!###############################################################################
!> Module for user inputs (direct and indirect).
!###############################################################################
module ovar
  !
  use module_kind_types, only : wp,lk,ldk,fals,true
  use module_kind_types, only : NavierStokes_Eqns
  use module_kind_types, only : Internally_Generated_Grid,Unknown_Format
  use module_kind_types, only : Legendre_Gauss,Equi_Distant_Points
  use module_kind_types, only : TVD3_RK
  use module_kind_types, only : invs_flux_ROE
  use module_kind_types, only : visc_flux_BR2
  use module_kind_types, only : visc_interface_fluxes
  use module_kind_types, only : laminar_flow
  use module_kind_types, only : one_work_unit
  use module_kind_types, only : iout,zero,one,half,two,five,pi,eps6,eps12
  use module_kind_types, only : mpi_wtime
  use module_kind_types, only : Correction_gDG
  use module_kind_types, only : TFM_OFF,IGN_OFF
  use module_kind_types, only : bc_unknown
  use module_kind_types, only : project_derivatives
  use module_kind_types, only : derivative_of_proj_fluxes
  use module_kind_types, only : WilliamsShunn_TriPoints
  use module_kind_types, only : ShunnHam_TetPoints
  use module_kind_types, only : LagrangePolynomial
  use module_kind_types, only : Global_Cell_Timestep
  use module_kind_types, only : debug_timer
  use module_kind_types, only : entering_procedure
  use module_kind_types, only : leaving_procedure
  !
  implicit none
  !
  private
  !
  public :: gfr_cpu_time
  public :: ovar_memory_usage
  public :: adjust_cfl_cycles
  public :: if_interpto_newgrid
  public :: move_bcin_to_device
  !
  !
  !
  ! ##############################
  ! ###  Module Derived Types  ###
  ! ##############################
  !
  !  relaxation_t :
  !> Derived type that provides relaxation parameters
  !> that can be used for certain boundary conditions.
  type, public :: relaxation_t
    real(wp) :: value = one
    real(wp) :: time  = zero
    integer  :: iter  = 0
  contains
    procedure :: is_not_default => is_not_default_relax
  end type relaxation_t
  !
  !  sem_bc_t :
  !> Derived type to provide information to use the
  !> synthetic eddy method (SEM) to generate an
  !> unsteady turbulent inflow boundary condition.
  type, public :: sem_bc_t
    logical(lk)        :: is_cylindrical      = fals
    real(wp)           :: length_scale        = zero
    real(wp)           :: convective_velocity = zero
    real(wp)           :: profile_height      = zero
    character(len=150) :: profile_file        = ""
    integer            :: sub_in_bc_grp       = 0
    integer            :: sub_in_bc_idx       = 0
  contains
    procedure :: is_not_default => is_not_default_sem
  end type sem_bc_t
  !
  !  periodic_bc_t :
  !> Derived type to provide information about a periodic
  !> boundary condition to allow it to be matched to another
  !> periodic boundary condition.
  !>
  !> Note:
  !> * At least one of matching_bc_name or matching_bc_idx needs to have
  !>   a non-default value. If both have been given non-default values,
  !>   the bc_input_t entity that they point to needs to be the same.
  !> * The rotation arrays define the geometrical transformation required
  !>   to find the matching faces from the connecting periodic boundary.
  !>   These are based off the implementation of periodic interfaces
  !>   within the CGNS standard. For more insight, see the
  !>   Periodic_t data structure in the CGNS SIDS.
  type, public :: periodic_bc_t
    !
    !  matching_bc_name :
    !> Name of the matching periodic boundary, i.e.:
    !> * matching_bc_name == bc_input(i)%name
    !> * matching_bc_idx = i
    character(len=32) :: matching_bc_name = ""
    !
    !  matching_bc_idx :
    !> Index of the matching periodic boundary interface, i.e.:
    !> * i = matching_bc_idx
    !> * matching_bc_name = bc_input(i)%name
    integer :: matching_bc_idx = 0
    !
    !  RotationCeter :
    !> The origin for defining the rotation angle between the periodic
    !> interfaces,
    real(wp) :: RotationCenter(1:3) = [zero,zero,zero]
    !
    !  RotationAngle :
    !> The axis-rotation angle from the current periodic
    !> interface. If rotating about more than one axis, the
    !> rotation is performed first about the x-axis, then the
    !> y-axis, then the z-axis.
    !>
    !> NOTE: EITHER RotationAngle OR RotationAxis can be used
    !>       to define rotated periodic boundary conditions,
    !>       but not both. If both are specified, RotationAxis
    !>       will be ignored and RotationAngle will be used.
    real(wp) :: RotationAngle(1:3)  = [zero,zero,zero]
    !
    !  RotationAxis :
    !> The axis of rotation as a vector and rotation angle about this axis
    !> * RotationAxis(1:3) = vector providing axis of rotation
    !> * RotationAxis(4)   = rotation angle about this axis
    !>
    !> NOTE: EITHER RotationAngle OR RotationAxis can be used
    !>       to define rotated periodic boundary conditions,
    !>       but not both. If both are specified, RotationAxis
    !>       will be ignored and RotationAngle will be used.
    real(wp) :: RotationAxis(1:4)  = [zero,zero,zero,zero]
    !
    !  Translation :
    !> The translation from the current periodic boundary interface to the
    !> connecting periodic boundary interface.
    real(wp) :: Translation(1:3)    = [zero,zero,zero]
    !
  contains
    procedure :: is_not_default => is_not_default_periodic
  end type periodic_bc_t
  !
  !  bc_input_t :
  !> Derived type used in the input namelist to specify
  !> the flow parameters for each boundary condition.
  type, public :: bc_input_t
    !
    character(len=32) :: name = ""
    !
    integer :: set_bc_value = bc_unknown
    !
    real(wp) :: t_static   = -one
    real(wp) :: p_static   = -one
    real(wp) :: rho_static = -one
    real(wp) :: mach       = -one
    real(wp) :: vx         = zero
    real(wp) :: vy         = zero
    real(wp) :: vz         = zero
    real(wp) :: t_total    = -one
    real(wp) :: p_total    = -one
    real(wp) :: rho_total  = -one
    real(wp) :: alpha_aoa  = zero
    real(wp) :: beta_aoa   = zero
    real(wp) :: wall_temp  = -one
    !
    ! NOTE: WE NEED TO MAKE SURE WE ADD LOGIC THAT SETS ROTATION_RATE(1)
    !       TO THE ROTATION FOR 2D. TECHNICALLY, ROTATION_RATE(3) SHOULD
    !       SPECIFY THE ROTATION FOR 2D AND YOU USE THE 3D CROSS PRODUCT
    !       (OMEGA x XYZ) WITH THE Z COORDINATE SET TO ZERO TO GET THE
    !       X AND Y ANGULAR VELOCITIES, BUT WE DONT WANT TO BE USING A
    !       'Z' VALUE FOR 2D.
    real(wp) :: rotation_rate(1:3) = [zero,zero,zero]
    !
    logical :: is_specified_in_mrf = fals
    !
#ifdef ALLOCATABLE_MASS_FRAC
    real(wp), allocatable :: mass_frac(:)
    character(len=10), allocatable :: species_name(:)
#else
    real(wp) :: mass_frac(7) = zero
    character(len=10) :: species_name(7) = ""
#endif
    !
    type(relaxation_t) :: relax
    !
    type(sem_bc_t) :: sem
    !
    type(periodic_bc_t) :: periodic
    !
  contains
    !
    procedure :: is_not_default => is_not_default_bc_input
    procedure :: has_default_flow_conditions
    procedure :: copy_nonflow_conditions
  !!procedure, private :: write_bc_input_t
  !!procedure, private :: read_bc_input_t
    !
#ifndef DISABLE_DTIO
  !!generic :: write(formatted) => write_bc_input_t
  !!generic :: read(formatted) => read_bc_input_t
#endif
    !
  end type bc_input_t
  !
  type, public :: pbc_rotation_t
    real(wp), allocatable :: mat(:,:)
  end type pbc_rotation_t
  !
  !  bc_in_t :
  !> Derived type that contains the boundary condition flow parameters
  !> used by the solver. This is a consolidated version of what is
  !> specified in the input namelist using bc_input_t.
  type, public :: bc_in_t
    real(wp), allocatable :: cv(:)
    real(wp), allocatable :: pv(:)
    real(wp), allocatable :: gp(:)
    real(wp) :: wall_temperature = -one
    real(wp) :: rotation_rate(1:3) = [zero,zero,zero]
    type(relaxation_t) :: relax
    type(sem_bc_t) :: sem
    type(pbc_rotation_t) :: pbc_rotation
    integer :: bc_input_idx
    integer :: bc_type = bc_unknown
  end type bc_in_t
  !
  !  mrf_spinup_t :
  !> Derived type for moving reference frame (MRF) spinup information.
  type, public :: mrf_spinup_t
    integer  :: iter = 0
    real(wp) :: tstar = zero
  contains
    procedure :: get_spinup_ratio => get_mrf_spinup_ratio
  end type mrf_spinup_t
  !
  !  mrf_t :
  !> Derived type for moving reference frame (MRF) information.
  type, public :: mrf_t
    real(wp) :: rotation(1:3) = [zero,zero,zero]
    real(wp) :: translation(1:3) = [zero,zero,zero]
    type(mrf_spinup_t) :: spinup
  contains
    procedure :: is_active => mrf_is_active
    procedure :: is_not_active => mrf_is_not_active
    procedure :: make_nondim => nondimensionalize_mrf
    procedure :: get_rotation => get_mrf_rotation
  end type mrf_t
  !
  !  thermofix_t
  !> Derived type to provide information to use the check/fix
  !> of bad thermo properties (pressure,temperature,density).
  type, public :: thermofix_t
    integer     :: method          = TFM_OFF
    integer     :: iterations      = 10
    logical(lk) :: check_cells     = true
    logical(lk) :: check_faces     = true
    logical(lk) :: check_nodes     = true
    logical(lk) :: verbose         = true
    real(wp)    :: scale_factor    = 0.5_wp     ! between [0,1]
    real(wp)    :: min_pressure    = 10000._wp  ! [Pa]
    real(wp)    :: min_temperature =   100._wp  ! [K]
    real(wp)    :: min_density     =  0.001_wp  ! [kg/m3]
  end type thermofix_t
  !
  ! point_probe_t : derived type containing point probe location
  !
  type, public :: point_probe_t
    !> point probe location in user grid units,
    !> must multiply by grid_scaling_units to get meters.
    real(wp) :: xyz(3) = [-999.0_wp,-999.0_wp,-999.0_wp]
  contains
    procedure :: is_not_default => is_not_default_point_probe
  end type point_probe_t
  !
  ! rake_probe_t : derived type containing rake probe location
  !
  type, public :: rake_probe_t
    !> number of probes in this rake
    integer  :: num_probes = 0
    !> start location of rake rake probes in user grid units,
    !> must multiply by grid_scaling_units to get meters.
    real(wp) :: xyz_beg(3) = [-999.0_wp,-999.0_wp,-999.0_wp]
    !> end   location of rake rake probes in user grid units,
    !> must multiply by grid_scaling_units to get meters.
    real(wp) :: xyz_end(3) = [-999.0_wp,-999.0_wp,-999.0_wp]
  contains
    procedure :: is_not_default => is_not_default_rake_probe
  end type rake_probe_t
  !
  !  geom_wts_t :
  !> Derived type to store geometry weights used for parallel load balancing.
  type, public :: geom_wts_t
    integer :: tria = 1
    integer :: quad = 1
    integer :: tetr = 1
    integer :: pyra = 1
    integer :: pris = 1
    integer :: hexa = 1
    integer :: face_tria = 1
    integer :: face_quad = 1
  end type geom_wts_t
  !
  !
  ! ################################
  ! ###  Module Local Variables  ###
  ! ################################
  !
  ! NOTE: Input defaults are set to run P3 48**3 DoF Taylor-Green vortex
  !       using gDG, Ga-LP, SSP-RK3, and Roe flux with entropy fix.
  !
  ! #############################
  ! ####   Input Variables   ####
  ! #############################
  !
  ! GEOMETRY WEIGHTS FOR METIS PARTITIONING
  type(geom_wts_t), public, save :: geom_wts(1:5)
  !
  ! GRID INFORMATION
  integer,          public, save :: grid_format = Internally_Generated_Grid
  character(len=:), public, save, allocatable :: gridfile
  character(len=:), public, save, allocatable :: plot3d_cutfile
  character(len=:), public, save, allocatable :: plot3d_bcsfile
  real(wp),         public, save :: grid_scaling_factor = 1.0_wp
  integer,          public, save :: internal_grid_order = -1
  !
  ! PLOT3D COARSENING OPTIONS
  integer, public, save :: plot3d_i_stride          = -1
  integer, public, save :: plot3d_j_stride          = -1
  integer, public, save :: plot3d_k_stride          = -1
  integer, public, save :: plot3d_all_strides       = -1
  integer, public, save :: plot3d_agglomerate_order = -1
  integer, public, save :: parametric_agglomeration =  0
  !
  ! These interp_* variables store inputs for the grid that
  ! will be interpolated onto.
  !
  ! INTERPOLATION UTILITY GRID INFORMATION
  integer,          public, save :: interpto_grid_format = Unknown_Format
  character(len=:), public, save, allocatable :: interpto_gridfile
  character(len=:), public, save, allocatable :: interpto_plot3d_cutfile
  character(len=:), public, save, allocatable :: interpto_plot3d_bcsfile
  real(wp),         public, save :: interpto_grid_scaling_factor = 1.0_wp
  !
  ! INTERPOLATION UTILITY PLOT3D COARSENING OPTIONS
  integer, public, save :: interpto_plot3d_i_stride          = -1
  integer, public, save :: interpto_plot3d_j_stride          = -1
  integer, public, save :: interpto_plot3d_k_stride          = -1
  integer, public, save :: interpto_plot3d_all_strides       = -1
  integer, public, save :: interpto_plot3d_agglomerate_order = -1
  !
  ! BASE FLUX RECONSTRUCTION (FR) OPTIONS
  integer, public, save :: flux_divergence     = LagrangePolynomial
  integer, public, save :: correction_function = Correction_gDG
  !
  ! LOCATION OF SOLUTION, FLUX, TRIANGLE, AND TETRAHEDRA POINTS
  integer, public, save :: loc_solution_pts   = Legendre_Gauss
  integer, public, save :: loc_flux_pts       = Legendre_Gauss
  integer, public, save :: loc_triangle_pts   = WilliamsShunn_TriPoints
  integer, public, save :: loc_tetrahedra_pts = ShunnHam_TetPoints
  !
  ! LOCATION OF INITIALIZATION POINTS
  integer, public, save :: loc_init_pts = 0 ! uses same as solution points
  !
  ! RUNGE-KUTTA OPTIONS
  integer, public, save :: Runge_Kutta_Scheme = TVD3_RK
  integer, public, save :: num_rk_stages      = 3
  !
  ! COMMON-FLUX OPTIONS
  integer, public, save :: invs_flux_method      = invs_flux_ROE
  integer, public, save :: visc_flux_method      = visc_flux_BR2
  integer, public, save :: visc_variables_to_ave = visc_interface_fluxes
  !
  ! TIME STEPPING INFORMATION
  integer,  public, save :: num_timesteps     = 25000
  real(wp), public, save :: Final_Time        = -10.0_wp
  real(wp), public, save :: constant_dt       = 1.0e-6_wp
  integer,  public, save :: Timestep_Type     = Global_Cell_Timestep
  integer,  public, save :: minimum_timesteps = 100
  !
  ! CFL CONTROL VARIABLES
  real(wp), public, save :: cfl            = 0.1_wp
  real(wp), public, save :: cfl_beg        = 0.1_wp
  real(wp), public, save :: cfl_end        = 0.1_wp
  integer,  public, save :: cfl_cycles     = 1
  integer,  public, save :: cfl_start_iter = 0
  !
  ! OVER-INTEGRATION OPTIONS
  integer, public, save :: loc_quadrature_pts              = Legendre_Gauss
  integer, public, save :: combine_common_flux_projections = 0
  !
  ! ERROR-INTEGRATION OPTIONS
  integer, public, save :: loc_error_pts = Legendre_Gauss
  !
  ! FORCE GLOBAL ORDER AND SOLUTION POINTS
  integer, public, save :: all_points = -1
  !
  ! PARALLEL OPTIONS
  integer, public, save :: metis_option = 3 ! k-way part, min comm volume
  integer, public, save :: use_parmetis = 0
  !
  ! INITIALIZATION / TEST CASE
  integer, public, save :: itestcase     = 9
  integer, public, save :: jet_init_func = 2
  !
  ! RESTART OPTIONS
  logical(lk),      public, save :: load_restart_file              = fals
  integer,          public, save :: restart_interval               = 0
  character(len=:), public, save, allocatable :: restart_file
  logical(lk),      public, save :: restart_output_uses_host_roots = true
  real(wp),         public, save :: CFL_reduction_for_NaN_restart  = 2.0_wp
  integer,          public, save :: num_saved_time_ave_restarts    = 2
  real(wp),         public, save :: restart_frequency              = -one
  logical(lk),      public, save :: completely_disable_restart_output = fals
  !
  ! OPTIONS TO USE SOLUTION FILE FOR INITIAL CONDITIONS
  integer,          public, save :: solution_format = Unknown_Format
  character(len=:), public, save, allocatable :: solution_file
  !
  ! OPTIONS FOR SEM
  logical(lk), public, save :: use_improved_SEM_normalization     = fals
  logical(lk), public, save :: use_SEM_profile_for_initialization = fals
  logical(lk), public, save :: SEM_eddy_volume_is_spherical       = fals
  real(wp),    public, save :: SEM_box_num_eddies_multiplier      = one
  logical(lk), public, save :: SEM_use_random_seed_from_restart   = fals
  logical(lk), public, save :: SEM_include_thermo_fluctuations    = fals
  logical(lk), public, save :: SEM_use_restart_eddies             = true
  !
  ! FLAG TO ACTIVATE A TEST OF THE SEM CODE
  logical(lk), public, save :: average_sem_inflow_plane = fals
  !
  ! STATISTICAL OUTPUT OPTIONS
  integer,     public, save :: iter_out_interval       = 1
  integer,     public, save :: results_interval        = 0
  logical(lk), public, save :: output_eperk_cfl        = true
  !
  !  disable_header_rewrites :
  !> Disables rewrites of the column headers for iteration statistics
  !> that are normally written to stdout every 100 iter_out_intervals.
  !> Does not affect the header written at the beginning of the run.
  logical(lk), public, save :: disable_header_rewrites = fals
  !
  logical(lk), public, save :: output_summary_file = fals
  !
  ! SOLUTION OUTPUT OPTIONS
  character(len=:), public, save, allocatable :: output_dir
  integer,          public, save :: output_interval   = -5
  real(wp),         public, save :: output_frequency  = -1.0_wp
  logical(lk),      public, save :: continuous_output = fals
  integer,          public, save :: loc_output_pts    = Equi_Distant_Points
  !
  ! SELECTION OF OUTPUT VARIABLES IN THE CGNS SOLUTION FILES
  logical(lk), public, save :: output_vorticity                      = fals
  logical(lk), public, save :: output_cylindrical_coordinates        = fals
  logical(lk), public, save :: output_cylindrical_velocity           = fals
  logical(lk), public, save :: output_relative_velocity              = fals
  logical(lk), public, save :: output_relative_conditions            = fals
  logical(lk), public, save :: output_qcriterion                     = fals
  logical(lk), public, save :: output_conserved_variables            = fals
  logical(lk), public, save :: output_momentum                       = fals
  logical(lk), public, save :: output_temperature                    = fals
  logical(lk), public, save :: output_thermo_variables               = fals
  logical(lk), public, save :: output_gradients                      = fals
  logical(lk), public, save :: output_density_gradient               = fals
  logical(lk), public, save :: output_density_gradient_magnitude     = fals
  logical(lk), public, save :: output_pressure_gradient              = fals
  logical(lk), public, save :: output_pressure_gradient_magnitude    = fals
  logical(lk), public, save :: output_pressure_dmat_gradient         = fals
  logical(lk), public, save :: output_temperature_gradient           = fals
  logical(lk), public, save :: output_temperature_gradient_magnitude = fals
  logical(lk), public, save :: output_temperature_dmat_gradient      = fals
  logical(lk), public, save :: output_schlieren                      = fals
  logical(lk), public, save :: output_modified_schlieren             = fals
  logical(lk), public, save :: output_transport_variables            = fals
  logical(lk), public, save :: output_stagnation_temperature         = fals
  logical(lk), public, save :: output_stagnation_pressure            = fals
  logical(lk), public, save :: output_stagnation_variables           = fals
  logical(lk), public, save :: output_marker_variables               = fals
  logical(lk), public, save :: output_grid_details                   = fals
  logical(lk), public, save :: output_local_timestep                 = fals
  logical(lk), public, save :: output_species_density                = fals
  logical(lk), public, save :: output_species_mass_fractions         = true
  logical(lk), public, save :: output_species_mole_fractions         = fals
  logical(lk), public, save :: output_extra_species_variables        = fals
  logical(lk), public, save :: output_debug_species_variables        = fals
  logical(lk), public, save :: output_grdnt_species_variables        = fals
  !
  ! GOVERNING EQUATIONS
  integer, public, save :: governing_equations = NavierStokes_Eqns
  integer, public, save :: turbulence_model    = laminar_flow
  !
  ! MOVING REFERENCE FRAME
  type(mrf_t), public, save :: moving_ref_frame
  !
  ! CHECK/FIX BAD THERMO PROPERTIES (PRES,TEMP,DENS)
  type(thermofix_t), public :: thermofix
  !
  ! CHEMICAL REACTION SOURCE TERMS
  logical(lk), public, save :: disable_chemical_reactions         = fals
  integer,     public, save :: semi_implicit_chemistry            = 0
  logical(lk), public, save :: traceback_unconverged_newton_temps = fals
 !integer,     public, save :: reactions_start_iter               = 0
  !
  ! IGNITION CONTROL SOURCE TERMS
  integer,  public, save :: &
      ignition_method = IGN_OFF, & ! method to use for the ignition source term
      ignition_start_timestep = 0, & ! timestep to begin ignition
      ignition_duration = 1   ! number of timesteps to maintain max temperature
  real(wp), public, save :: &
      ignition_prod2reac_L2_ratio_goal = -1.0_wp, &
      ignition_max_temperature         = 1600.0_wp, &   ! [K]
      ignition_spark_cap               = 1.0e+12_wp, &  ! [N/m^2]
      ! if > zero, then
      !   ignition_source_magnitude = delta temperature [K]
      ! if < zero, then
      !   abs(ignition_source_magnitude) = % of ref. total energy [no units]
      ignition_source_magnitude        = 20.0_wp, &
      ignition_source_radius           = zero, &
      ignition_source_location(1:3)    = [zero,zero,zero]
  !
  ! NONDIMENSIONAL REFERENCE CONDITIONS
  real(wp), public, save :: machref =  0.4_wp
  real(wp), public, save :: reyref  = -1.0_wp
  ! reyref default evaluates to 8.7319053606e+6
  real(wp), public, save :: gam     =  1.4_wp
  !
  ! REALIZABLE TOTAL TEMPERATURE LIMITER
  real(wp), public, save :: ideal_total_temperature = -1.0_wp
  real(wp), public, save :: limit_total_temperature = 0.2_wp
  !
  ! REFERENCE CONDITIONS
  real(wp), public, save :: rhoref  =     -1.0_wp
  ! rhoref default evaluates to 1.160833478437518_wp
  real(wp), public, save :: tref    =    300.0_wp
  real(wp), public, save :: pref    = 100000.0_wp
  real(wp), public, save :: rgasref =     -1.0_wp
  !
  ! REFERENCE TURBULENCE CONDITIONS
  real(wp), public, save :: tinref = 0.0_wp
  real(wp), public, save :: tvrref = 1.0e-6_wp
  !
  ! REFERENCE ANGLES OF ATTACK
  real(wp), public, save :: aoa_alpharef = 0.0_wp
  real(wp), public, save :: aoa_betaref  = 0.0_wp
  !
  ! LAMINAR AND TURBULENT PRANDTL NUMBERS
  real(wp), public, save :: Pr  = 0.72_wp
  real(wp), public, save :: Prt = 0.72_wp
  !
  ! LAMINAR AND TURBULENT SCHMIDT NUMBERS
  real(wp), public, save :: schmidt_lam  = 0.72_wp
  real(wp), public, save :: schmidt_turb = 0.72_wp
  !
  ! SUTHERLANDS LAW COEFFICIENTS
  real(wp), public, save :: suth_muref = 1.716e-5_wp
  real(wp), public, save :: suth_Tref  = 273.0_wp
  real(wp), public, save :: suth_Sref  = 110.4_wp
  !
  ! FILTER OPTIONS
  integer,  public, save :: Filter_Option = 0
  real(wp), public, save :: Filter_etac   = 0.0_wp
  integer,  public, save :: Filter_Order  = 16
  real(wp), public, save :: Filter_Alpha  = 36.0_wp
  real(wp), public, save :: Filter_eps    = 0.0_wp
  !
  integer,  public, save :: Pyramid_Filter_Option = 1
  !
  ! LIMITER FILTER OPTIONS
  integer,  public, save :: Limiter_Filter_Option = 0
  real(wp), public, save :: Limiter_Filter_etac   = 0.0_wp
  integer,  public, save :: Limiter_Filter_Order  = 16
  real(wp), public, save :: Limiter_Filter_Alpha  = 36.0_wp
  real(wp), public, save :: Limiter_Filter_eps    = 0.0_wp
  logical(lk), public, save :: Limiter_Use_Dynamic_Filter = fals
  logical(lk), public, save :: Track_Limiter_Filter_Dissipation = fals
  !
  ! SPONGE FILTER OPTIONS
  logical(lk), public, save :: use_sponge_filter    = fals
  real(wp),    public, save :: Sponge_Filter_etac   = 0.0_wp
  integer,     public, save :: Sponge_Filter_Order  = 8
  real(wp),    public, save :: Sponge_Filter_Alpha  = 36.0_wp
  real(wp),    public, save :: Sponge_Beyond_X      = huge(zero)
  !
  ! LIMITER OPTIONS
  integer, public, save :: Limiter_Option             = 0
  integer, public, save :: Resolution_Marker_Variable = 0
  integer, public, save ::       Jump_Marker_Variable = 0
  integer, public, save ::        MLP_Marker_Variable = 0
  !
  ! LIMITER MARKER SENSITIVITY OPTIONS
  real(wp), public, save :: Limiter_Marker_Min  =  1.0_wp
  real(wp), public, save :: Limiter_Marker_Max  =  1.0_wp
  real(wp), public, save :: Limiter_Marker_c1   = -1.0_wp
  real(wp), public, save :: Limiter_Marker_c2   = -1.0_wp
  real(wp), public, save :: Limiter_Marker_dpsi = -1.0_wp
  !
  ! CONTROL OVER FREQUENCY OF FILTERING/LIMITING
  integer, public, save :: Filtering_Interval        = 0
  integer, public, save :: Limiting_Interval         = 0
  integer, public, save :: Limiting_Interval_Output  = 1
  !< frequency to store limiter marker and scaled marker data:
  !<   =0  store each Runge-Kutta stage,
  !<   >0  store last Runge-Kutta stage,
  !<   <0  mimic the frequency of Limiting_Interval.
  !< Storing markers each stage increases memory and file size.
  !
  ! OUTPUT THE STATISTICS OF THE LIMITER
  logical(lk), public, save :: output_limiter_statistics = fals
  !
  ! OPTION TO OUTPUT FILTER DIAGONALS TO FILE
  logical(lk), public, save :: output_filter_diagonals = fals
  !
  ! ISENTROPIC VORTEX CONSTANTS
  real(wp), public, save :: vortex_bigR     = 1.0_wp
  real(wp), public, save :: vortex_bigGam   = 5.0_wp
  real(wp), public, save :: VortexGrid_Lref = 10.0_wp
  !
  ! TAYLOR-GREEN VORTEX CONSTANTS
  integer,     public, save :: VortexGrid_DOF             = 48
  logical(lk), public, save :: VortexGrid_random_perturb  = fals
  integer,     public, save :: VortexGrid_perturb_percent = 30
  logical(lk), public, save :: VortexGrid_Hex_to_Pyramids = fals
  logical(lk), public, save :: VortexGrid_is_3D           = fals
  !
  ! TAYLOR-COUETTE FLOW CONSTANTS
  real(wp), public, save :: TCF_Inner_Radius       = half
  real(wp), public, save :: TCF_Inner_Wall_Spacing = 0.002_wp
  real(wp), public, save :: TCF_Outer_Radius       = one
  real(wp), public, save :: TCF_Outer_Wall_Spacing = 0.0025_wp
  real(wp), public, save :: TCF_Aspect_Ratio       = two*pi
  integer,  public, save :: TCF_Axial_Cells        = 24
  integer,  public, save :: TCF_Radial_Cells       = 20
  integer,  public, save :: TCF_Azimuthal_Cells    = 48
  real(wp), public, save :: TCF_Theta_Shift        = zero
  integer,  public, save :: TCF_Periodic_Sectors   = 1
  ! TCF_Periodic_Sectors : periodic rotation angle = 360/TCF_Periodic_Sectors
  !
  ! FLAME SPEED TEST PARAMETERS
  integer, public, save :: flame_speed_test_dof = 0
  !
  ! SOD SHOCK TUBE TEST PARAMETERS
  integer,  public, save :: sod_test_dof = 0
    !< Sod shock tube resolution for 1d internal grid:
    !<   (<0) user specified number of cells,
    !<   (=0) 512 degrees of freedom,
    !<   (>0) user specified degrees of freedom.
  real(wp), public, save :: sod_test_right_vector(3) = [1.,0.,0.]
    !< Sod shock tube diaphram right side normal vector.
  !
  ! DENSITY PULSE CONVECTION TEST PARAMETERS
  integer,  public, save :: pulse_test_form = 0
    !< Form of the density pulse test shape:
    !<   0 = Combo (Gauss,Square,Triangle,Ellipse),
    !<   1 = Gauss,
    !<   2 = Square,
    !<   3 = Triangle,
    !<   4 = Ellipse.
  integer,  public, save :: pulse_test_form_dim = 1
    !< Dimension of the density test shape:
    !<   1 = 1d,
    !<   2 = 2d.
  integer,  public, save :: pulse_test_dof = 0
    !< Resolution for density pulse test internal grid:
    !<   (<0) user specified number of cells,
    !<   (=0) 512 degrees of freedom,
    !<   (>0) user specified degrees of freedom.
  real(wp), public, save :: pulse_test_length = 1.0
    !< 1d domain length for density pulse test internal grid.
  real(wp), public, save :: pulse_test_center(3) = [0.5,0.0,0.0]
    !< xyz at center of first density pulse.
  real(wp), public, save :: pulse_test_width  = 0.2
    !< Width of density pulse.
  real(wp), public, save :: pulse_test_gap    = 0.2
    !< Gap between density pulses.
  integer,  public, save :: pulse_test_repeat = 0
    !< Number of repeat density pulses.
  !
  ! SHU-OSHER SHOCK TUBE TEST PARAMETERS
  integer,  public, save :: shu_osher_test_dof = 0
    !< Shu-Osher shock tube resolution for 1d internal grid:
    !<   (<0) user specified number of cells,
    !<   (=0) 512 degrees of freedom,
    !<   (>0) user specified degrees of freedom.
  !
  ! CHANNEL FLOW PARAMETERS
  real(wp),         public, save :: channel_body_force          = 0.0_wp
  character(len=:), public, save, allocatable :: channel_init_file
  logical(lk),      public, save :: channel_spatial_ave_restart = fals
  logical(lk),      public, save :: channel_nondim_output_coord = fals
  !
  real(wp), public, save :: channel_yplus_of_one       = 8.42134e-6_wp
  real(wp), public, save :: channel_grid_ave_xplus     = -30.0_wp
  real(wp), public, save :: channel_grid_ave_yplus     = -15.0_wp
  real(wp), public, save :: channel_grid_ave_zplus     = -15.0_wp
  real(wp), public, save :: channel_grid_yplus_wall    = 2.0_wp
  integer,  public, save :: channel_rake_save_interval = 100
  !
  ! PIPE CHANNEL FLOW PARAMETERS
  real(wp), public, save :: pipe_channel_angle = five
  !
  ! OPTIONS FOR MODIFYING BC ALGORITHMS
  logical(lk), public, save :: walls_are_exact   = true
  integer,     public, save :: sub_inflow_method = 4
  logical(lk), public, save :: use_bc_conflicts  = fals
  !
  ! CGNS OUTPUT OPTIONS
  logical(lk), public, save :: convert_restart_to_cgns = fals
  logical(lk), public, save :: delete_converted_restart_files = fals
  logical(lk), public, save :: completely_disable_cgns = fals
  logical(lk), public, save :: interpolate_before_output_variable = true
  logical(lk), public, save :: cgns_high_order_output = fals
  logical(lk), public, save :: cgns_mixed_output = true
  logical(lk), public, save :: cgns_output_restart_time_ave_data = fals
  !
  ! CGNS INPUT OPTIONS
  logical(lk), public, save :: correct_cgns_grid_points = fals
  logical(lk), public, save :: force_cgns_linear_grid_points = fals
  !
  ! CONVERGENCE CRITERIA
  real(wp), public, save :: convergence_order_abs_res = real( -16 , kind=wp )
  real(wp), public, save :: convergence_order_max_res = real( -16 , kind=wp )
  !
  ! PROFILE PARALLEL IO
  logical(lk), public, save :: profile_io_all       = fals
  logical(lk), public, save :: profile_io_cgns      = fals
  logical(lk), public, save :: profile_io_cgns_subs = fals
  logical(lk), public, save :: profile_io_restart   = fals
  !
  ! OPTIONS TO OUTPUT BOUNDARY PROFILE
  logical(lk), public, save :: output_bnd_profiles = fals
  !
  ! CONTROL DOT_PRODUCT PROCEDURES BOUND TO MATRIX DERIVED-TYPE
  logical(lk), public, save :: use_unrolled_dot_products = true
  !
  ! LUSTRE FILE SYSTEM STRIPING PARAMETERS
  integer,          public, save :: lustre_stripe_count = 64
  character(len=:), public, save, allocatable :: lustre_stripe_size
  !
  ! OPENACC INPUT OPTIONS
  logical(lk), public, save :: use_async_for_error = fals
  !
  ! POINT PROBES
  !> probes sampling interval in iterations
  integer, public :: probes_collection_interval = 0
  integer, public :: probes_write_interval = 100
  !> list of probes
  type(point_probe_t), public, allocatable :: point_probes(:)
  type(rake_probe_t),  public, allocatable :: rake_probes(:)
  type(point_probe_t), public, allocatable :: cell_probes(:)
  !
  ! TIME AVERAGING
  logical(lk),      public, save :: output_time_averaging    = fals
  character(len=:), public, save, allocatable :: time_ave_file
  real(wp),         public, save :: time_scaling_factor      = 1.0_wp
  logical(lk),      public, save :: time_ave_vel_is_axisymm  = fals
  integer,          public, save :: start_averaging_timestep = -1
  real(wp),         public, save :: start_averaging_tstar    = -one
  !
  ! SPATIAL AVERAGE RESTART FILES
  logical(lk), public, save :: spatial_ave_restart_x_and_z_planes = fals
  logical(lk), public, save :: spatial_ave_restart_z_planes       = fals
  !
  ! UTILITY TO TIME AVERAGE EXISTING RESTART FILES
  logical(lk),      public, save :: time_average_restart_files = fals
  character(len=:), public, save, allocatable :: restart_list_file
  logical(lk),      public, save :: modify_time_ave_restart_history = fals
  !
  ! MEMORY USAGE
  logical(lk), public, save :: dump_memory_usage = fals
  !
  ! DUMP FLUX INFORMATION FOR DEBUGGING
  logical(lk), public, save :: dump_fluxes                  = fals
  logical(lk), public, save :: check_flux_point_coordinates = fals
  !
  ! ACTIVATE SPECIAL GRID FOR DEBUGGING THE POSTPROCESSOR
  logical(lk), public, save :: postproc_debug_grid = fals
  !
  ! OPTION TO OUTPUT PLOT3D GRID TO TECPLOT IN UNSTRUCTURED FORMAT
  logical(lk), public, save :: output_plot3d_to_tecplot = fals
  !
  ! OPTION TO OUTPUT BOUNDARY CONDITIONS TO TECPLOT FILE
  logical(lk), public, save :: output_bface_array = fals
  !
  logical(lk), public, save :: force_nd_subroutines = fals
  real(wp),    public, save :: newton_rlx = 1.0_wp
  real(wp),    public, save :: conv_temp = eps12
  integer,     public, save :: chemistry_sub_iterations = -1
  logical(lk), public, save :: ignore_thermo_newton_warnings = fals
  !
  ! INPUT BOUNDARY CONDITIONS
  type(bc_input_t), public, allocatable :: bc_input(:)
  !
  ! ###################################
  ! ####   Input Alias Variables   ####
  ! ###################################
  !
  ! Note:
  ! - Alias variables are only to make user input easier.
  ! - They should only be used by the namelist routines and
  !   during the post-input check to set other internal variables.
  ! - Never access these variables directly in the main code.
  !   Instead use the internal variable denoted in the trailing
  !   comments below.
  !
  ! POLYNOMIAL ORDER
  integer, public, save :: all_orders = -3                  ! sets [npq]_order
  integer, public, save :: solution_order = 3               ! sets n_order
  integer, public, save :: projection_order = -3            ! sets p_order
  integer, public, save :: quadrature_order = -3            ! sets q_order
  integer, public, save :: output_order = -3                ! sets o_order
  integer, public, save :: error_order = -3                 ! sets e_order
  !
  integer, public, save :: interpto_all_orders = -3         ! sets [npq]_order
  integer, public, save :: interpto_solution_order = 3      ! sets n_order
  integer, public, save :: interpto_projection_order = -3   ! sets p_order
  integer, public, save :: interpto_quadrature_order = -3   ! sets q_order
  integer, public, save :: interpto_output_order = -3       ! sets o_order
  integer, public, save :: interpto_error_order = -3        ! sets e_order
  !
  integer, public, save :: dmat_grad_temp_order = -1        ! sets t_order
  !
  !
  ! #################################
  ! ####   Non-Input Variables   ####
  ! #################################
  !
  real(wp), public, parameter :: Default_Marker_c1_RI = 4.0_wp
  real(wp), public, parameter :: Default_Marker_c2_RI = 4.25_wp
  !
  real(wp), public, parameter :: Default_Marker_c1_JI = 2.25_wp
  real(wp), public, parameter :: Default_Marker_c2_JI = 3.0_wp
  !
  real(wp), public, parameter :: Default_Marker_dpsi = 0.5_wp
  real(wp), public, save      :: Default_Marker_c1   = Default_Marker_c1_RI
  real(wp), public, save      :: Default_Marker_c2   = Default_Marker_c2_RI
  !
 !integer, public, parameter :: dflux_opt = project_derivatives
  integer, public, parameter :: dflux_opt = derivative_of_proj_fluxes
  !
  character(len=15), public, save :: lustre_stripe_factor = "64"
  character(len=15), public, save :: lustre_stripe_unit   = "4194304"
  !
  real(wp), public, save :: Lref = -1.0_wp
  real(wp), public, save :: aref
  real(wp), public, save :: uref
  real(wp), public, save :: cpref
  real(wp), public, save :: muref
  ! muref default evaluates to 1.84671531762365e-5_wp
  real(wp), public, save, dimension(3) :: velref
  !
  real(wp), public, save :: pref_nd
 !real(wp), public, save :: rgasref
  real(wp), public, save :: rgasref_nd
  real(wp), public, save :: muref_nd
  real(wp), public, save :: cpref_nd
  real(wp), public, save, dimension(3) :: velref_nd
  !
  real(wp), public, save :: tmin = 10.0_wp
  real(wp), public, save :: tmax = 15000.0_wp
  !
  real(wp), public, save :: tkeref
  real(wp), public, save :: tlsref
  !
  real(wp), public, save :: tkeref_nd
  real(wp), public, save :: tlsref_nd
  !
  logical(lk),  public, save :: turb_is_on = .false.
  logical(lk),  public, save :: include_tke_in_energy = .false.
  !
  real(wp), public, save :: dt_global
  !
  integer,  public, save :: itcur
  real(wp), public, save :: current_cfl
  real(wp), public, save :: d_t
  real(wp), public, save :: time
  real(wp), public, save :: time_ref = 1.0_wp
  integer,  public, save :: current_rk_stage = 1
  !
  real(wp),     public, save :: saved_ave_time = 0.0_wp
  real(wp),     public, save :: prev_ave_time = 0.0_wp
  real(wp),     public, save :: ave_start_time = 0.0_wp
  logical(ldk), public, save :: read_time_ave_restart = .false.
  !
  real(wp), public, save :: init_tstar = 0.0_wp
  !
  real(wp), public, save :: cfl_cycle_start
  real(wp), public, save :: cfl_cycle_end
  !
  integer,  public, save :: itrst = 0
  real(wp), public, save :: tmrst = zero
  real(wp), public, save :: dtrst = zero
  !
  real(wp), public, save :: marker_dpsi = Default_Marker_dpsi
  real(wp), public, save :: marker_c1   = Default_Marker_c1_RI
  real(wp), public, save :: marker_c2   = Default_Marker_c2_RI
  !
  ! jump indicator: weight and var-index for gas prop and usp
  real(wp), public, save :: ji_gwt
  real(wp), public, save :: ji_uwt
  integer,  public, save :: ji_gvar
  integer,  public, save :: ji_uvar
  !
  real(wp), public, save, allocatable :: cv_ref(:)
  real(wp), public, save, allocatable :: pv_ref(:)
  real(wp), public, save, allocatable :: gp_ref(:)
  !
  real(wp), public, save :: iso_wall_temp = 300.0_wp
  !
  real(wp), public, save :: p_exit
  real(wp), public, save :: local_volume
  real(wp), public, save :: total_volume
  !
  real(wp), public, save :: betark(4,4)
  !
  real(wp), public, save :: rl2mx = -huge(zero)
  !
  real(wp), public, save, allocatable :: uverr(:,:,:)
  !$ACC DECLARE CREATE(uverr)
  !
  integer, public, save, allocatable :: nlim(:)
  !
  real(wp), public, save, dimension(1:3,1:3,1:5) :: advec_error
  real(wp), public, save, dimension(1:5,1:2,1:5) :: vortex_error
  !
  integer, public, save, allocatable :: all_host_num(:)
  integer, public, save, allocatable :: all_cpu_num(:)
  !
  integer,  public, save :: Total_Periods
  !
  logical(lk),  public, save :: is_at_output_frequency = fals
  logical(lk),  public, save :: is_at_restart_frequency = fals
  !
  logical(lk),  public, save :: write_TaylorGreen_8s_solution = fals
  logical(lk),  public, save :: write_TaylorGreen_full_solution = fals
  !
  logical(lk),  public, save :: output_KEDR = fals
  !
  logical(lk),  public, save :: convergence_reached = fals
  logical(lk),  public, save :: this_is_final_timestep = fals
  !
  character(len=170), public, save :: cgns_dir = ""
  logical(lk),        public, save :: disable_cgns_output = fals
  logical(lk),        public, save :: cgns_grid_needs_correction = fals
  !
  real(wp),         save :: max_wall_time
  real(wp),         save :: beg_wall_time
  real(wp),         save :: end_wall_time
  real(wp), public, save :: elapsed_wall_time
  real(wp), public, save :: work_units
  real(wp), public, save :: total_sbu_cost = zero
  !
  real(wp), public, save :: max_restart_wall_time       = zero
  real(wp), public, save :: max_time_ave_rst_wall_time  = zero
  real(wp), public, save :: max_cgns_wall_time          = zero
  real(wp), public, save :: max_time_ave_cgns_wall_time = zero
  !
  real(wp), public, save :: Period_Time
  integer,  public, save :: Period_Timesteps
  real(wp), public, save, dimension(1:3) :: Grid_Max
  real(wp), public, save, dimension(1:3) :: Grid_Min
  real(wp), public, save, dimension(1:3) :: Grid_Length
  !
  real(wp), public, save :: suth_c1
  !
  real(wp), public, save :: LDG_beta(1:3) = [0.0_wp, 0.0_wp, 0.0_wp]
  real(wp), public, save :: LDG_tau = 0.0_wp
  !
  type(bc_in_t), public, save, allocatable :: bc_in(:)
  !
  character(len=32), public, save, allocatable :: bc_names(:)
  !
  integer, public, save :: interpto_input_mode = 1
  !
  !> Turn quadrature/matrix output on by setting matrix_output_unit to anything
  !> besides iout (which is a parameter variable set to the iso_fortran_env
  !> constant output_unit, which historically has the value 6).
  integer, public, save :: matrix_output_unit = iout
 !integer, public, save :: matrix_output_unit = 9203
  !
  logical(lk), save :: mrf_is_in_nondim = fals
  !
  type, public :: resource_list_t
    integer  :: n = 1
    real(wp) :: sbu_rate = zero
    character(len=:), allocatable :: model
  end type resource_list_t
  !
  type(resource_list_t), public, save, allocatable :: resource_list(:)
  !
  !> Interface to check if two boundary condition input derived types are equal.
  interface operator(==)
    module procedure bc_input_t_equals, &
                     relaxation_t_equals, &
                     sem_bc_t_equals, &
                     periodic_bc_t_equals, &
                     point_probe_t_equals, &
                     rake_probe_t_equals
  end interface
  !
  !> Interface to check if two boundary condition input derived types are not
  !> equal.
  interface operator(/=)
    module procedure bc_input_t_not_equal, &
                     relaxation_t_not_equal, &
                     sem_bc_t_not_equal, &
                     periodic_bc_t_not_equal, &
                     point_probe_t_not_equal, &
                     rake_probe_t_not_equal
  end interface
  !
  !$ACC DECLARE CREATE(walls_are_exact,sub_inflow_method, &
  !$ACC&               turb_is_on,include_tke_in_energy, &
  !$ACC&               time,itcur,gam,pr,aref,rhoref,tref,tmin,tmax, &
  !$ACC&               cpref_nd,rgasref_nd, &
  !$ACC&               suth_Tref,suth_Sref,suth_muref,suth_c1, &
  !$ACC&               LDG_beta,LDG_tau,visc_variables_to_ave, &
  !$ACC&               ideal_total_temperature,limit_total_temperature, &
  !$ACC&               invs_flux_method,Resolution_Marker_Variable, &
  !$ACC&               Jump_Marker_Variable,MLP_Marker_Variable)
contains
!
!###############################################################################
!> Move the bc_in derived type to GPU.
!###############################################################################
subroutine move_bcin_to_device
  !
  !.. Local Scalars ..
  integer :: i,i1,i2,m1,m2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "move_bcin_to_device"
  !
continue
  !
  if (allocated(bc_in)) then
    !$ACC ENTER DATA COPYIN(bc_in)
    i1 = lbound(bc_in,dim=1)
    i2 = ubound(bc_in,dim=1)
    !$ACC ENTER DATA COPYIN(bc_in(i1:i2))
    do i = i1,i2
      if (allocated(bc_in(i)%cv)) then
        !$ACC ENTER DATA COPYIN(bc_in(i)%cv)
        m1 = lbound(bc_in(i)%cv,dim=1)
        m2 = ubound(bc_in(i)%cv,dim=1)
        !$ACC ENTER DATA COPYIN(bc_in(i)%cv(m1:m2))
      end if
      if (allocated(bc_in(i)%pv)) then
        !$ACC ENTER DATA COPYIN(bc_in(i)%pv)
        m1 = lbound(bc_in(i)%pv,dim=1)
        m2 = ubound(bc_in(i)%pv,dim=1)
        !$ACC ENTER DATA COPYIN(bc_in(i)%pv(m1:m2))
      end if
      if (allocated(bc_in(i)%gp)) then
        !$ACC ENTER DATA COPYIN(bc_in(i)%gp)
        m1 = lbound(bc_in(i)%gp,dim=1)
        m2 = ubound(bc_in(i)%gp,dim=1)
        !$ACC ENTER DATA COPYIN(bc_in(i)%gp(m1:m2))
      end if
      !$ACC ENTER DATA COPYIN(bc_in(i)%relax)
    end do
  end if
  !
end subroutine move_bcin_to_device
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!> Check if a bc_input derived type is not set to default values.
!###############################################################################
pure function is_not_default_bc_input(this) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(bc_input_t), intent(in) :: this
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Allocatable Scalars ..
  type(bc_input_t), allocatable :: bc_input_default
  !
continue
  !
  ! Create a variable of type bc_input_t with default values
  !
  allocate ( bc_input_default )
  !
  ! Compare the input argument to the default values
  !
  return_value = ( this /= bc_input_default )
  !
end function is_not_default_bc_input
!
!###############################################################################
!> Check if the flow condition components of a bc_input derived type
!> have default values (within some tolerance).
!###############################################################################
pure function has_default_flow_conditions(this) result(return_value)
  !
  ! This function is slightly different than is_not_default_bc_input. This
  ! function is used to check only the flow condition components of the
  ! boundary condition against the default values. The function
  ! is_not_default_bc_input also checks bc_input%name and other non flow
  ! condition components against the default values.
  !
  !.. Passed Argument from Invoking Object ..
  class(bc_input_t), intent(in) :: this
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Allocatable Scalars ..
  type(bc_input_t), allocatable :: dflt
  !
continue
  !
  ! Create a variable of type bc_input_t with default values
  !
  allocate ( dflt )
  !
  ! Initialize the function result to false
  !
  return_value = fals
  !
  ! NOTE: Dont check this%name, this%set_bc_value, or this%periodic because
  !       those three components dont have anything to do with the flow
  !       conditions for this boundary condition
  !
  ! If any of the following conditionals are true, then the subcomponents
  ! of the input argument are not equal to the default values, so immediately
  ! return with a function result of false indicating they are not equal.
  !
  if (abs(this%t_static   - dflt%t_static  ) > eps6) return
  if (abs(this%p_static   - dflt%p_static  ) > eps6) return
  if (abs(this%rho_static - dflt%rho_static) > eps6) return
  if (abs(this%mach       - dflt%mach      ) > eps6) return
  if (abs(this%vx         - dflt%vx        ) > eps6) return
  if (abs(this%vy         - dflt%vy        ) > eps6) return
  if (abs(this%vz         - dflt%vz        ) > eps6) return
  if (abs(this%t_total    - dflt%t_total   ) > eps6) return
  if (abs(this%p_total    - dflt%p_total   ) > eps6) return
  if (abs(this%rho_total  - dflt%rho_total ) > eps6) return
  if (abs(this%alpha_aoa  - dflt%alpha_aoa ) > eps6) return
  if (abs(this%beta_aoa   - dflt%beta_aoa  ) > eps6) return
  if (abs(this%wall_temp  - dflt%wall_temp ) > eps6) return
  !
  do n = 1,size(this%rotation_rate)
    if (abs(this%rotation_rate(n) - dflt%rotation_rate(n)) > eps6) return
  end do
  !
  if (this%is_specified_in_mrf .neqv. dflt%is_specified_in_mrf) return
  !
  if ( this%relax%is_not_default() ) return
  !
  if ( this%sem%is_not_default() ) return
  !
  ! CHECK ON this%periodic%is_not_default() IS NOT DONE BECAUSE WHETHER OR
  ! NOT THE PERIODIC COMPONENT HAS DEFAULT CONDITIONS HAS NOTHING TO DO WITH
  ! WHETHER OR NOT THIS BC_INPUT ELEMENT HAS DEFAULT FLOW CONDITIONS WHICH
  ! IS WHAT THIS FUNCTION IS CHECKING.
  !
#ifdef ALLOCATABLE_MASS_FRAC
  !
  if (allocated(this%mass_frac   ) .neqv. allocated(dflt%mass_frac   )) return
  if (allocated(this%species_name) .neqv. allocated(dflt%species_name)) return
  !
  if (allocated(this%mass_frac)) then
    if (size(this%mass_frac) /= size(dflt%mass_frac)) return
    do n = 1,size(this%mass_frac)
      if (abs(this%mass_frac(n) - dflt%mass_frac(n)) > eps6) return
    end do
  end if
  !
  if (allocated(this%species_name)) then
    if (size(this%species_name) /= size(dflt%species_name)) return
    do n = 1,size(this%species_name)
      if (trim(adjustl(this%species_name(n))) /= &
          trim(adjustl(dflt%species_name(n)))) return
    end do
  end if
  !
#else
  !
  do n = 1,size(this%mass_frac)
    if (abs(this%mass_frac(n) - dflt%mass_frac(n)) > eps6) return
  end do
  !
  do n = 1,size(this%species_name)
    if (trim(adjustl(this%species_name(n))) /= &
        trim(adjustl(dflt%species_name(n)))) return
  end do
  !
#endif
  !
  return_value = true
  !
end function has_default_flow_conditions
!
!###############################################################################
!
pure subroutine copy_nonflow_conditions(this,copy_src)
  !
  ! This function copies the non-flow-condition parts of the input copy_src
  ! to the current element of type bc_input_t.
  !
  !.. Passed Argument from Invoking Object ..
  class(bc_input_t), intent(inout) :: this
  class(bc_input_t),    intent(in) :: copy_src
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Allocatable Scalars ..
  type(bc_input_t), allocatable :: dflt
  !
continue
  !
  ! Copy the name component
  !
  this%name = trim(adjustl(copy_src%name))
  !
  ! Copy the BC type
  !
  this%set_bc_value = copy_src%set_bc_value
  !
  ! Copy the relax component and its subcomponents
  !
  this%relax = copy_src%relax
  !
  ! Copy the sem component and its subcomponents
  !
  this%sem = copy_src%sem
  !
  ! Copy the periodic component and its subcomponents
  !
  this%periodic = copy_src%periodic
  !
end subroutine copy_nonflow_conditions
!
!###############################################################################
!> Check if two bc_input derived type variables are equal.
!###############################################################################
pure function bc_input_t_equals(a,b) result(return_value)
  !
  !.. Formal Arguments ..
  type(bc_input_t), intent(in) :: a
  type(bc_input_t), intent(in) :: b
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! Initialize the function result to false
  !
  return_value = fals
  !
  ! If any of the following conditionals are true, then the subcomponents
  ! of the input arguments are not equal, so immediately return with a
  ! function result of false indicating they are not equal.
  !
  if (trim(adjustl(a%name)) /= trim(adjustl(b%name))) return
  !
  if (a%set_bc_value /= b%set_bc_value) return
  !
  if (abs(a%t_static   - b%t_static  ) > eps6) return
  if (abs(a%p_static   - b%p_static  ) > eps6) return
  if (abs(a%rho_static - b%rho_static) > eps6) return
  if (abs(a%mach       - b%mach      ) > eps6) return
  if (abs(a%vx         - b%vx        ) > eps6) return
  if (abs(a%vy         - b%vy        ) > eps6) return
  if (abs(a%vz         - b%vz        ) > eps6) return
  if (abs(a%t_total    - b%t_total   ) > eps6) return
  if (abs(a%p_total    - b%p_total   ) > eps6) return
  if (abs(a%rho_total  - b%rho_total ) > eps6) return
  if (abs(a%alpha_aoa  - b%alpha_aoa ) > eps6) return
  if (abs(a%beta_aoa   - b%beta_aoa  ) > eps6) return
  if (abs(a%wall_temp  - b%wall_temp ) > eps6) return
  !
  do n = 1,size(a%rotation_rate)
    if (abs(a%rotation_rate(n) - b%rotation_rate(n)) > eps6) return
  end do
  !
  if (a%is_specified_in_mrf .neqv. b%is_specified_in_mrf) return
  !
#ifdef ALLOCATABLE_MASS_FRAC
  !
  if (allocated(a%mass_frac   ) .neqv. allocated(b%mass_frac   )) return
  if (allocated(a%species_name) .neqv. allocated(b%species_name)) return
  !
  if (allocated(a%mass_frac)) then
    if (size(a%mass_frac) /= size(b%mass_frac)) return
    do n = 1,size(a%mass_frac)
      if (abs(a%mass_frac(n) - b%mass_frac(n)) > eps6) return
    end do
  end if
  !
  if (allocated(a%species_name)) then
    if (size(a%species_name) /= size(b%species_name)) return
    do n = 1,size(a%species_name)
      if (trim(adjustl(a%species_name(n))) /= &
          trim(adjustl(b%species_name(n)))) return
    end do
  end if
  !
#else
  !
  do n = 1,size(a%mass_frac)
    if (abs(a%mass_frac(n) - b%mass_frac(n)) > eps6) return
  end do
  !
  do n = 1,size(a%species_name)
    if (trim(adjustl(a%species_name(n))) /= &
        trim(adjustl(b%species_name(n)))) return
  end do
  !
#endif
  !
  if (a%relax    /= b%relax   ) return
  if (a%sem      /= b%sem     ) return
  if (a%periodic /= b%periodic) return
  !
  ! If we have reached this point, all subcomponents of the input arguments
  ! are equal, so return true to indicate this is the case.
  !
  return_value = true
  !
end function bc_input_t_equals
!
!###############################################################################
!> Check if two bc_input derived type variables are not equal.
!###############################################################################
pure function bc_input_t_not_equal(a,b) result(return_value)
  !
  !.. Formal Arguments ..
  type(bc_input_t), intent(in) :: a
  type(bc_input_t), intent(in) :: b
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
continue
  !
  return_value = .not. (a == b)
  !
end function bc_input_t_not_equal
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!> Check if a relaxation derived type variable has non-default values.
!###############################################################################
pure function is_not_default_relax(this) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(relaxation_t), intent(in) :: this
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Allocatable Scalars ..
  type(relaxation_t), allocatable :: relax_default
  !
continue
  !
  ! Create a variable of type relaxation_t with default values
  !
  allocate ( relax_default )
  !
  ! Compare the input argument to the default values
  !
  return_value = ( this /= relax_default )
  !
end function is_not_default_relax
!
!###############################################################################
!> Check if two relaxation derived type variables are equal.
!###############################################################################
elemental function relaxation_t_equals(a,b) result(return_value)
  !
  !.. Formal Arguments ..
  type(relaxation_t), intent(in) :: a
  type(relaxation_t), intent(in) :: b
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
continue
  !
  ! Initialize the function result to false
  !
  return_value = fals
  !
  ! If any of the following conditionals are true, then the subcomponents
  ! of the input arguments are not equal, so immediately return with a
  ! function result of false indicating they are not equal.
  !
  if (abs(a%value - b%value) > eps6) return
  if (abs(a%time  - b%time ) > eps6) return
  if (a%iter /= b%iter) return
  !
  ! If we have reached this point, all subcomponents of the input arguments
  ! are equal, so return true to indicate this is the case.
  !
  return_value = true
  !
end function relaxation_t_equals
!
!###############################################################################
!> Check if two relaxation derived type variables are not equal.
!###############################################################################
elemental function relaxation_t_not_equal(a,b) result(return_value)
  !
  !.. Formal Arguments ..
  type(relaxation_t), intent(in) :: a
  type(relaxation_t), intent(in) :: b
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
continue
  !
  return_value = .not. (a == b)
  !
end function relaxation_t_not_equal
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!> Check if a sem_bc derived type variable has non-default values.
!###############################################################################
pure function is_not_default_sem(this) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(sem_bc_t), intent(in) :: this
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Allocatable Scalars ..
  type(sem_bc_t), allocatable :: sem_default
  !
continue
  !
  ! Create a variable of type sem_bc_t with default values
  !
  allocate ( sem_default )
  !
  ! Compare the input argument to the default values
  !
  return_value = ( this /= sem_default )
  !
end function is_not_default_sem
!
!###############################################################################
!> Check if two sem_bc derived type variables are equal.
!###############################################################################
elemental function sem_bc_t_equals(a,b) result(return_value)
  !
  !.. Formal Arguments ..
  type(sem_bc_t), intent(in) :: a
  type(sem_bc_t), intent(in) :: b
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
continue
  !
  ! Initialize the function result to false
  !
  return_value = fals
  !
  ! If any of the following conditionals are true, then the subcomponents
  ! of the input arguments are not equal, so immediately return with a
  ! function result of false indicating they are not equal.
  !
  if (a%is_cylindrical .neqv. b%is_cylindrical) return
  !
  if (abs(a%length_scale        - b%length_scale       ) > eps6) return
  if (abs(a%convective_velocity - b%convective_velocity) > eps6) return
  if (abs(a%profile_height      - b%profile_height     ) > eps6) return
  !
  if (trim(adjustl(a%profile_file)) /= trim(adjustl(b%profile_file))) return
  !
  if (a%sub_in_bc_grp /= b%sub_in_bc_grp) return
  if (a%sub_in_bc_idx /= b%sub_in_bc_idx) return
  !
  ! If we have reached this point, all subcomponents of the input arguments
  ! are equal, so return true to indicate this is the case.
  !
  return_value = true
  !
end function sem_bc_t_equals
!
!###############################################################################
!> Check if two sem_bc derived type variables are not equal.
!###############################################################################
elemental function sem_bc_t_not_equal(a,b) result(return_value)
  !
  !.. Formal Arguments ..
  type(sem_bc_t), intent(in) :: a
  type(sem_bc_t), intent(in) :: b
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
continue
  !
  return_value = .not. (a == b)
  !
end function sem_bc_t_not_equal
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!> Check if a periodic_bc derived type variable has non-default values.
!###############################################################################
pure function is_not_default_periodic(this) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(periodic_bc_t), intent(in) :: this
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Allocatable Scalars ..
  type(periodic_bc_t), allocatable :: periodic_default
  !
continue
  !
  ! Create a variable of type periodic_t with default values
  !
  allocate ( periodic_default )
  !
  ! Compare the input argument to the default values
  !
  return_value = ( this /= periodic_default )
  !
end function is_not_default_periodic
!
!###############################################################################
!> Check if two periodic_bc derived type variables are equal.
!###############################################################################
elemental function periodic_bc_t_equals(a,b) result(return_value)
  !
  !.. Formal Arguments ..
  type(periodic_bc_t), intent(in) :: a
  type(periodic_bc_t), intent(in) :: b
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! Initialize the function result to false
  !
  return_value = fals
  !
  ! If any of the following conditionals are true, then the subcomponents
  ! of the input arguments are not equal, so immediately return with a
  ! function result of false indicating they are not equal.
  !
  if (trim(adjustl(a%matching_bc_name)) /= &
      trim(adjustl(b%matching_bc_name))) return
  !
  if (a%matching_bc_idx /= b%matching_bc_idx) return
  !
  do n = 1,size(a%RotationCenter)
    if (abs(a%RotationCenter(n) - b%RotationCenter(n)) > eps6) return 
  end do
  !
  do n = 1,size(a%RotationAngle)
    if (abs(a%RotationAngle(n) - b%RotationAngle(n)) > eps6) return 
  end do
  !
  do n = 1,size(a%RotationAxis)
    if (abs(a%RotationAxis(n) - b%RotationAxis(n)) > eps6) return 
  end do
  !
  do n = 1,size(a%Translation)
    if (abs(a%Translation(n) - b%Translation(n)) > eps6) return 
  end do
  !
  ! If we have reached this point, all subcomponents of the input arguments
  ! are equal, so return true to indicate this is the case.
  !
  return_value = true
  !
end function periodic_bc_t_equals
!
!###############################################################################
!> Check if two periodic_bc derived type variables are not equal.
!###############################################################################
elemental function periodic_bc_t_not_equal(a,b) result(return_value)
  !
  !.. Formal Arguments ..
  type(periodic_bc_t), intent(in) :: a
  type(periodic_bc_t), intent(in) :: b
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
continue
  !
  return_value = .not. (a == b)
  !
end function periodic_bc_t_not_equal
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
!
!###############################################################################
!> Check if a point probe derived type is not set to default values.
!###############################################################################
pure function is_not_default_point_probe(this) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(point_probe_t), intent(in) :: this
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Allocatable Scalars ..
  type(point_probe_t), allocatable :: point_probe_default
  !
continue
  !
  ! Create a variable of type point_probe_t with default values
  !
  allocate ( point_probe_default )
  !
  ! Compare the input argument to the default values
  !
  return_value = ( this /= point_probe_default )
  !
end function is_not_default_point_probe
!
!###############################################################################
!> Check if two point probe derived type variables are equal.
!###############################################################################
elemental function point_probe_t_equals(a,b) result(return_value)
  !
  !.. Formal Arguments ..
  type(point_probe_t), intent(in) :: a
  type(point_probe_t), intent(in) :: b
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! Initialize the function result to false
  !
  return_value = fals
  !
  ! If any of the following conditionals are true, then the subcomponents
  ! of the are not equal, so immediately return with a function result of
  ! false indicating they are not equal.
  !
  do n = 1, size(a%xyz)
    if (abs(a%xyz(n) - b%xyz(n)) > eps6) return
  end do
  !
  ! Passed all checks, a and b are equal.
  !
  return_value = true
  !
end function point_probe_t_equals
!###############################################################################
!> Check if two point probe derived type variables are not equal.
!###############################################################################
elemental function point_probe_t_not_equal(a,b) result(return_value)
  !
  !.. Formal Arguments ..
  type(point_probe_t), intent(in) :: a
  type(point_probe_t), intent(in) :: b
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
continue
  !
  return_value = .not. (a == b)
  !
end function point_probe_t_not_equal
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
!
!###############################################################################
!> Check if a rake probe derived type is not set to default values.
!###############################################################################
pure function is_not_default_rake_probe(this) result(return_value)
  !
  !.. Use Statements ..
  use module_kind_types, only : rake_eps => eps4
  !
  !.. Passed Argument from Invoking Object ..
  class(rake_probe_t), intent(in) :: this
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Allocatable Scalars ..
  type(rake_probe_t), allocatable :: rake_probe_default
  !
continue
  !
  ! Create a variable of type rake_probe_t with default values
  !
  allocate ( rake_probe_default )
  !
  ! Compare the input argument to the default values
  !
  return_value = ( this /= rake_probe_default )
  !
  ! If this does not contain default values, make sure the values for this
  ! rake are valid, i.e., the distance between the rake beginning and end
  ! is reasonable, and 1 or more probes were requested for the rake.
  ! NOTE: If only one probe is requested for the rake, ignore the distance
  !       between the beginning and end of the rake since those will just
  !       get averaged to get the probe location.
  !
  ! If the values for this rake are invalid, make the return_value from this
  ! function indicate this rake has default values so this invalid rake
  ! gets deleted from the list of valid rakes.
  !
  if (return_value) then
    !
    if (this%num_probes <= 0) then
      !
      return_value = fals
      !
    else if (this%num_probes > 1 .and. &
             norm2( this%xyz_beg - this%xyz_end ) < rake_eps) then
      !
      return_value = fals
      !
    end if
    !
  end if
  !
end function is_not_default_rake_probe
!
!###############################################################################
!> Check if two rake probe derived type variables are equal.
!###############################################################################
elemental function rake_probe_t_equals(a,b) result(return_value)
  !
  !.. Formal Arguments ..
  type(rake_probe_t), intent(in) :: a
  type(rake_probe_t), intent(in) :: b
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! Initialize the function result to false
  !
  return_value = fals
  !
  ! If any of the following conditionals are true, then the subcomponents
  ! of the input arguments are not equal, so immediately return with a
  ! function result of false indicating they are not equal.
  !
  if (a%num_probes /= b%num_probes) return
  !
  do n = 1, size(a%xyz_beg)
    if (abs(a%xyz_beg(n) - b%xyz_beg(n)) > eps6) return
  end do
  do n = 1, size(a%xyz_end)
    if (abs(a%xyz_end(n) - b%xyz_end(n)) > eps6) return
  end do
  !
  ! Passed all checks, a and b are equal.
  !
  return_value = true
  !
end function rake_probe_t_equals
!###############################################################################
!> Check if two rake probe derived type variables are not equal.
!###############################################################################
elemental function rake_probe_t_not_equal(a,b) result(return_value)
  !
  !.. Formal Arguments ..
  type(rake_probe_t), intent(in) :: a
  type(rake_probe_t), intent(in) :: b
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
continue
  !
  return_value = .not. (a == b)
  !
end function rake_probe_t_not_equal
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
!
!###############################################################################
!> Check if a moving reference frame (MRF) derived type variable specifies
!> an active moving reference frame.
!###############################################################################
pure function mrf_is_active(this) result(return_value)
  !
  !.. Use Statements ..
  use module_kind_types, only : rpm_to_mps
  !
  !.. Passed Argument from Invoking Object ..
  class(mrf_t), intent(in) :: this
  !
  !.. Function Result ..
  logical :: return_value
  !
  !.. Local Scalars ..
  logical  :: in_trf,in_rrf
  real(wp) :: trf_dim,rrf_dim
  !
continue
  !
  if (mrf_is_in_nondim) then
    trf_dim = one/aref
    rrf_dim = rpm_to_mps/aref
  else
    trf_dim = one
    rrf_dim = one
  end if
  !
  in_rrf = any(abs(this%rotation   ) >= one*rrf_dim)  ! is in rotation
  in_trf = any(abs(this%translation) >= one*trf_dim)  ! is in translation
  !
  return_value = (in_rrf .or. in_trf)  ! active if rotating or translating
  !
end function mrf_is_active
!
!###############################################################################
!> Check if a moving reference frame (MRF) derived type variable specifies
!> an non-active moving reference frame.
!###############################################################################
pure function mrf_is_not_active(this) result(return_value)
  !.. Passed Argument from Invoking Object ..
  class(mrf_t), intent(in) :: this
  !.. Function Result ..
  logical :: return_value
continue
  return_value = ( .not. this%is_active() )
end function mrf_is_not_active
!
!###############################################################################
!> Nondimensionalize the information in the moving reference frame (MRF)
!> derived type, if not already nondimensional.
!###############################################################################
subroutine nondimensionalize_mrf(this)
  !
  !.. Use Statements ..
  use module_kind_types, only : rpm_to_mps
  !
  !.. Passed Argument from Invoking Object ..
  class(mrf_t), intent(inout) :: this
  !
continue
  !
  if (mrf_is_in_nondim) return
  !
  this%rotation(:) = this%rotation(:) * (rpm_to_mps/aref)
  !
  this%translation(:) = this%translation(:) / aref
  !
  mrf_is_in_nondim = true
  !
end subroutine nondimensionalize_mrf
!
!###############################################################################
!
pure &
function get_mrf_spinup_ratio(this) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(mrf_spinup_t), intent(in) :: this
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: current_tstar,spinup_iter,spinup_tstar,spinup_ratio
  !
continue
  !
  return_value = one
  !
  spinup_iter = one
  spinup_tstar = one
  !
  ! First check the current iteration versus the spinup%iter value
  !
  if (itcur < this%iter) then
    spinup_iter = real(itcur,kind=wp) / real(this%iter,kind=wp)
  end if
  !
  ! Next check the current tstar versus the spinup%tstar value
  !
  if (this%tstar > eps6) then
    current_tstar = time*time_ref*time_scaling_factor
    if (current_tstar < this%tstar) then
      spinup_tstar = current_tstar / this%tstar
    end if
  end if
  !
  ! Set the spinup ratio as the smaller of the spinup iteration ratio and
  ! the spinup tstar ratio, and make sure the value is between zero and one.
  !
  return_value = max( zero , min( spinup_tstar, spinup_iter, one ) )
  !
end function get_mrf_spinup_ratio
!
!###############################################################################
!> Get the moving reference frame (MRF) rotation speeds.
!###############################################################################
pure &
subroutine get_mrf_rotation(this,x_r,y_r,z_r,full_speed)
  !
  !.. Passed Argument from Invoking Object ..
  class(mrf_t), intent(in) :: this
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: x_r,y_r,z_r
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: full_speed
  !
  !.. Local Scalars ..
  real(wp)    :: spinup_ratio
  logical(lk) :: ignore_spinup
  !
continue
  !
  ignore_spinup = fals
  if (present(full_speed)) then
    ignore_spinup = full_speed
  end if
  !
  ! Get the spinup ratio
  !
  if (ignore_spinup) then
    spinup_ratio = one
  else
    spinup_ratio = this%spinup%get_spinup_ratio()
  end if
  !
  ! Multiply the final rotation speeds by the spinup ratio
  !
  x_r = spinup_ratio * this%rotation(1)
  y_r = spinup_ratio * this%rotation(2)
  z_r = spinup_ratio * this%rotation(3)
  !
end subroutine get_mrf_rotation
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
!
!###############################################################################
!> Write the data from a bc_input derived type variable to the specified
!> output unit number.
!###############################################################################
!!subroutine write_bc_input_t(dtv,iunit,iotype,vlist,my_iostat,my_iomsg)
!!  !
!!  !.. Passed Argument from Invoking Object ..
!!  class(bc_input_t), intent(in) :: dtv
!!  integer, intent(in) :: iunit
!!  character(len=*), intent(in) :: iotype
!!  integer, intent(in) :: vlist(:)
!!  integer, intent(out) :: my_iostat
!!  character(len=*), intent(inout) :: my_iomsg
!!  !
!!  integer :: m
!!  !
!!continue
!!  !
!!  write (iunit,1,iostat=my_iostat,iomsg=my_iomsg) &
!!    dtv%t_static, &
!!    dtv%p_static, &
!!    dtv%rho_static, &
!!    dtv%mach, &
!!    dtv%vx, &
!!    dtv%vy, &
!!    dtv%vz, &
!!    dtv%t_total, &
!!    dtv%p_total, &
!!    dtv%rho_total, &
!!    dtv%alpha_aoa, &
!!    dtv%beta_aoa, &
!!    dtv%wall_temp, &
!!    (dtv%rotation_rate(m),m=1,3), &
!!    dtv%relax%value, &
!!    dtv%relax%time, &
!!    dtv%relax%iter
!!#ifdef ALLOCATABLE_MASS_FRAC
!!  if (allocated(dtv%mass_frac)) then
!!    write (iunit,2,iostat=my_iostat,iomsg=my_iomsg) &
!!      (m,dtv%mass_frac(m), m=1,size(dtv%mass_frac))
!!  end if
!!#else
!!  write (iunit,2,iostat=my_iostat,iomsg=my_iomsg) &
!!    (m,dtv%mass_frac(m), m=1,size(dtv%mass_frac))
!!#endif
!!  !
!!  1 format ("t_static     = ",es14.6,/, &
!!            "p_static     = ",es14.6,/, &
!!            "rho_static     = ",es14.6,/, &
!!            "mach     = ",es14.6,/, &
!!            "vx     = ",es14.6,/, &
!!            "vy     = ",es14.6,/, &
!!            "vz     = ",es14.6,/, &
!!            "t_total     = ",es14.6,/, &
!!            "p_total     = ",es14.6,/, &
!!            "rho_total     = ",es14.6,/, &
!!            "alpha_aoa     = ",es14.6,/, &
!!            "beta_aoa     = ",es14.6,/, &
!!            "wall_temp     = ",es14.6,/, &
!!            "rotation_rate = [",es14.6,2(",",es14.6),"]",/, &
!!            "relax%value     = ",es14.6,/, &
!!            "relax%time     = ",es14.6,/, &
!!            "relax%iter     = ",es14.6,/)
!!  2 format (100("mass_frac(",i0,") = ",es14.6,:,/))
!!  !
!!end subroutine write_bc_input_t
!
!###############################################################################
!> Read the data into a bc_input derived type variable from the specified
!> input unit number.
!###############################################################################
!!subroutine read_bc_input_t(dtv,iunit,iotype,vlist,my_iostat,my_iomsg)
!!  !
!!  !.. Passed Argument from Invoking Object ..
!!  class(bc_input_t), intent(inout) :: dtv
!!  integer, intent(in) :: iunit
!!  character(len=*), intent(in) :: iotype
!!  integer, intent(in) :: vlist(:)
!!  integer, intent(out) :: my_iostat
!!  character(len=*), intent(inout) :: my_iomsg
!!  !
!!  integer :: m,ierr
!!  !
!!  !.. Local Parameters ..
!!  integer, parameter :: max_ncs = 7
!!  !
!!continue
!!  !
!!  read (iunit,*,iostat=my_iostat,iomsg=my_iomsg) &
!!    dtv%t_static, &
!!    dtv%p_static, &
!!    dtv%rho_static, &
!!    dtv%mach, &
!!    dtv%vx, &
!!    dtv%vy, &
!!    dtv%vz, &
!!    dtv%t_total, &
!!    dtv%p_total, &
!!    dtv%rho_total, &
!!    dtv%alpha_aoa, &
!!    dtv%beta_aoa, &
!!    dtv%wall_temp, &
!!    (dtv%rotation_rate(m),m=1,3), &
!!    dtv%relax%value, &
!!    dtv%relax%time, &
!!    dtv%relax%iter
!!#ifdef ALLOCATABLE_MASS_FRAC
!!  if (.not.allocated(dtv%mass_frac)) then
!!    allocate ( dtv%mass_frac(1:max_ncs) , source=zero )
!!  end if
!!#endif
!!  read (iunit,*,iostat=my_iostat,iomsg=my_iomsg) &
!!                                  (dtv%mass_frac(m), m=1,size(dtv%mass_frac))
!!  !
!!end subroutine read_bc_input_t
!
!###############################################################################
!> Function to determine whether code will interpolate onto new grid.
!###############################################################################
pure &
function if_interpto_newgrid() result(return_value)
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
continue
  !
  return_value = fals ! default
  !
  ! If user specifies interpto_gridfile filename,
  ! then code should interpolate onto it.
  !
  if (len_trim(adjustl(interpto_gridfile)) > 0) then
    return_value = true
  end if
  !
end function if_interpto_newgrid
!
!###############################################################################
!
subroutine adjust_cfl_cycles(new_cfl_inputs)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: new_cfl_inputs(:)
  !
  !.. Local Saved Scalars ..
  integer :: inp_cfl_cycles = 1
  integer :: inp_cfl_start_iter = 0
  !
  real(wp) :: inp_cfl     = -one
  real(wp) :: inp_cfl_beg = -one
  real(wp) :: inp_cfl_end = -one
  !
  logical(lk) :: cfl_adjustment_is_set = fals
  !
continue
  !
  if (present(new_cfl_inputs)) then
    !
    cfl_adjustment_is_set = true
    !
    inp_cfl = new_cfl_inputs(1)
    inp_cfl_beg = new_cfl_inputs(2)
    inp_cfl_end = new_cfl_inputs(3)
    inp_cfl_cycles = nint(new_cfl_inputs(4))
    inp_cfl_start_iter = nint(new_cfl_inputs(5))
    !
    if (inp_cfl_start_iter <= itcur) then
      cfl_start_iter = itcur
    else
      cfl_start_iter = inp_cfl_start_iter
    end if
    !
  else
    !
    if (cfl_adjustment_is_set) then
      !
      cfl = inp_cfl
      if (inp_cfl_cycles < 2) then
        cfl_beg = inp_cfl
        cfl_end = inp_cfl
        cfl_cycles = 1
      else
        cfl_beg = min(inp_cfl_beg,inp_cfl_end)
        cfl_end = max(inp_cfl_beg,inp_cfl_end)
        cfl_cycles = inp_cfl_cycles
      end if
      !
      cfl_adjustment_is_set = fals
      !
    end if
    !
    cfl_cycle_start = itcur
    cfl_cycle_end = itcur + cfl_cycles
    !
  end if
  !
end subroutine adjust_cfl_cycles
!
!###############################################################################
!
subroutine gfr_cpu_time(time_flag,io_unit)
  !
  !.. Use Statements ..
  use module_kind_types, only : gfr_timer
  !
  !.. Formal Arguments ..
  integer, intent(in) :: time_flag
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: io_unit
  !
  !.. Local Scalars ..
  integer :: hours,minutes,seconds,total_sec,io
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "gfr_cpu_time"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  io = iout
  if (present(io_unit)) then
    io = io_unit
  end if
  !
  if (time_flag == 0) then
    !
    call gfr_timer( beg_wall_time )
    work_units = zero
    !
  else if (time_flag == 1) then
    !
    call gfr_timer( end_wall_time )
    elapsed_wall_time = end_wall_time-beg_wall_time
    work_units = elapsed_wall_time/one_work_unit
    total_sbu_cost = get_sbu_cost(elapsed_wall_time)
    !
  else if (time_flag < 0) then
    !
    total_sec = nint(elapsed_wall_time)
    !
    hours   = total_sec / 3600
    minutes = (total_sec - hours*3600) / 60
    seconds = total_sec - hours*3600 - minutes*60
    !
    write (io,1) elapsed_wall_time, hours,minutes,seconds
    !
    if (time_flag == -99) then
      write (io,2) work_units
    else
      write (io,3) total_sbu_cost
    end if
    !
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  1 format (/,"TOTAL SIMULATION WALL TIME WAS ",f18.6," s.", &
              " (",i0,":",i2.2,":",i2.2,")")
  2 format (  "TOTAL COMPUTATIONAL COST WAS   ",f18.6," work units.")
  3 format (  "TOTAL SBU COST WAS             ",f18.6," SBUs.")
  !
end subroutine gfr_cpu_time
!
!###############################################################################
!
pure &
function get_sbu_cost(wall_time) result(return_value)
  !
  ! NOTE: WALL_TIME IS EXPECTED IN UNITS OF SECONDS
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: wall_time
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  integer  :: i
  real(wp) :: rnum_nodes,sbu_rate,sbu_time
  !
  !.. Local Parameters ..
  real(wp), parameter :: sec_to_hours = one/3600.0_wp
  !
continue
  !
  return_value = zero
  !
  ! Convert the input wall time from seconds to hours
  ! to get the time used to compute the total SBU cost
  !
  sbu_time = wall_time * sec_to_hours
  !
  if (allocated(resource_list)) then
    do i = 1,size(resource_list)
      rnum_nodes = real(resource_list(i)%n,kind=wp)
      sbu_rate = resource_list(i)%sbu_rate
      return_value = return_value + rnum_nodes*sbu_rate*sbu_time
    end do
  end if
  !
end function get_sbu_cost
!
!###############################################################################
!> Output memory requirement for ovar.
!###############################################################################
subroutine ovar_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
  !
  !.. Formal Arguments ..
  integer, optional, intent(in) :: iunit
  !
  !.. Local Scalars ..
  integer :: iou
  character(len=36) :: array_name
  type(memory) :: array
  type(memory) :: total
  !
#ifndef DISABLE_DTIO
continue
  !
  iou = iout
  if (present(iunit)) iou = iunit
  !
  call total%reset
  !
  write (iou,1)
  !
  if (allocated(cv_ref)) then
    call array%set( storage_size(cv_ref,kind=inttype) , &
                            size(cv_ref,kind=inttype) )
    call total%add(array)
    write (array_name,5) "cv_ref"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "cv_ref"
  end if
  if (allocated(pv_ref)) then
    call array%set( storage_size(pv_ref,kind=inttype) , &
                            size(pv_ref,kind=inttype) )
    call total%add(array)
    write (array_name,5) "pv_ref"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "pv_ref"
  end if
  if (allocated(gp_ref)) then
    call array%set( storage_size(gp_ref,kind=inttype) , &
                            size(gp_ref,kind=inttype) )
    call total%add(array)
    write (array_name,5) "gp_ref"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "gp_ref"
  end if
  if (allocated(uverr)) then
    call array%set( storage_size(uverr,kind=inttype) , &
                            size(uverr,kind=inttype) )
    call total%add(array)
    write (array_name,5) "uverr"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "uverr"
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module ovar")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine ovar_memory_usage
!
!###############################################################################
!
end module ovar
!
!###############################################################################
!
module input_namelist_mod
  !
  !.. Use Statements ..
  use ovar,      only : cfl, cfl_beg, cfl_end, cfl_cycles, cfl_start_iter
 !use ovar,      only : reactions_start_iter
  use ovar,      only : restart_file
  !
  implicit none
  !
  public
  !
  ! CFL ADJUSTMENT NAMELIST
  !
  namelist / cfl_adjustment / cfl, cfl_beg, cfl_end, cfl_cycles, cfl_start_iter
  !
  ! RESTART INTERUPTION NAMELIST
  !
  namelist / restart_interupt / restart_file
 !!
 !! START REACTIONS NAMELIST
 !!
 !namelist / start_reactions / reactions_start_iter
  !
contains
!
!###############################################################################
!
end module input_namelist_mod
!
!###############################################################################
!
