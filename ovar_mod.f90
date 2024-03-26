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
  use module_kind_types, only : iout,zero,one,half,two,pi,eps6,eps12
  use module_kind_types, only : mpi_wtime
  use module_kind_types, only : Correction_gDG
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
  ! relaxation_t : derived type that provides relaxation parameters
  !                that can be used for certain boundary conditions
  !
  type, public :: relaxation_t
    real(wp) :: value = one
    real(wp) :: time  = zero
    integer  :: iter  = 0
  contains
    procedure :: is_not_default => is_not_default_relax
  end type relaxation_t
  !
  ! sem_bc_t : derived type to provide information to use the
  !            synthetic eddy method (SEM) to generate an
  !            unsteady turbulent inflow boundary condition
  !
  type, public :: sem_bc_t
    logical(lk) :: is_cylindrical = fals
    real(wp) :: length_scale = zero
    real(wp) :: convective_velocity = zero
    real(wp) :: profile_height = zero
    character(len=150) :: profile_file = ""
    integer  :: sub_in_bc_grp = 0
    integer  :: sub_in_bc_idx = 0
  contains
    procedure :: is_not_default => is_not_default_sem
  end type sem_bc_t
  !
  ! periodic_bc_t : derived type to provide information about a periodic
  !                 boundary condition to allow it to be matched to another
  !                 periodic boundary condition
  !
  type, public :: periodic_bc_t
    !
    ! NOTE: At least one of matching_bc_name or matching_bc_idx needs to have
    !       a non-default value. If both have been given non-default values,
    !       the bc_input_t entity that they point to needs to be the same.
    !
    ! matching_bc_name : gives a character string that should match the value
    !                    of one of the other boundary condition names (that is
    !                    not itself), i.e.,
    !                       matching_bc_name == bc_input(i)%name
    !                       matching_bc_idx = i
    character(len=32) :: matching_bc_name = ""
    ! matching_bc_idx  : gives the index of bc_input(:) to find the matching
    !                    periodic boundary interface, i.e.,
    !                       i = matching_bc_idx
    !                       matching_bc_name = bc_input(i)%name
    integer :: matching_bc_idx = 0
    !
    ! The following three arrays define the geometrical transformation required
    ! to find the matching faces from the connecting periodic boundary.
    ! NOTE: These are based off the implementation of periodic interfaces
    !       within the CGNS standard. For more insight, see the
    !       Periodic_t data structure in the CGNS SIDS.
    !
    ! RotationCenter : the origin for defining the rotation angle
    !                  between the periodic interfaces
    real(wp) :: RotationCenter(1:3) = [zero,zero,zero]
    ! RotationAngle  : defines the axis-rotation angle from the current periodic
    !                  boundary interface to the connecting periodic boundary
    !                  interface. If rotating about more than one axis, the
    !                  rotation is performed first about the x-axis, then the
    !                  y-axis, then the z-axis.
    real(wp) :: RotationAngle(1:3)  = [zero,zero,zero]
    ! Translation    : defines the translation from the current periodic
    !                  boundary interface to the connecting periodic boundary
    !                  interface.
    real(wp) :: Translation(1:3)    = [zero,zero,zero]
    !
  contains
    procedure :: is_not_default => is_not_default_periodic
  end type periodic_bc_t
  !
  ! bc_input_t : derived type used in the input namelist to specify
  !              the flow parameters for each boundary condition
  !
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
  ! bc_in_t : derived type that contains the boundary condition flow parameters
  !           used by the solver. This is a consolidated version of what is
  !           specified in the input namelist using bc_input_t.
  !
  type, public :: bc_in_t
    real(wp), allocatable :: cv(:)
    real(wp), allocatable :: pv(:)
    real(wp), allocatable :: gp(:)
    real(wp) :: wall_temperature = -one
    real(wp) :: rotation_rate(1:3) = [zero,zero,zero]
    type(relaxation_t) :: relax
    type(sem_bc_t) :: sem
    integer :: bc_input_idx
    integer :: bc_type = bc_unknown
  end type bc_in_t
  !
  !
  type, public :: mrf_t
    real(wp) :: rotation(1:3) = [zero,zero,zero]
    real(wp) :: translation(1:3) = [zero,zero,zero]
  contains
    procedure :: is_active => mrf_is_active
    procedure :: is_not_active => mrf_is_not_active
    procedure :: make_nondim => nondimensionalize_mrf
  end type mrf_t
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
  ! BASE FR OPTIONS
  integer, public, save :: flux_divergence     = LagrangePolynomial
  integer, public, save :: correction_function = Correction_gDG
  !
  ! LOCATION OF SOLUTION, FLUX, AND TRIANGULAR POINTS
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
  ! COMMON FLUX OPTIONS
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
  !
  ! METHOD OF MANUFACTURED SOLUTIONS (MMS) OPTIONS
  integer,     public, save :: mms_opt             = 0
  real(wp),    public, save :: mms_init            = 0.1_wp
  logical(lk), public, save :: mms_output_solution = fals
  !
  ! INITIALIZATION / TEST CASE
  integer, public, save :: itestcase     = 9
  integer, public, save :: jet_init_func = 2
  !
  ! RESTART OPTIONS
  logical(lk),      public, save :: load_restart_file              = fals
  integer,          public, save :: restart_interval               = 0
  logical(lk),      public, save :: use_old_restart_format         = fals
  character(len=:), public, save, allocatable :: restart_file
  logical(lk),      public, save :: restart_output_uses_host_roots = true
  real(wp),         public, save :: CFL_reduction_for_NaN_restart  = 2.0_wp
  integer,          public, save :: num_saved_time_ave_restarts    = 2
  real(wp),         public, save :: restart_frequency              = -one
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
  !
  ! FLAG TO ACTIVATE A TEST OF THE SEM CODE
  logical(lk), public, save :: average_sem_inflow_plane = fals
  !
  ! STATISTICAL OUTPUT OPTIONS
  integer,     public, save :: iter_out_interval       = 1
  integer,     public, save :: results_interval        = 0
  logical(lk), public, save :: output_eperk_cfl        = true
  !     disable_header_rewrites : 
  !        The column headers for the iteration statistics are written to stdout
  !        before the first iteration, and by default are rewritten about once
  !        every 100 iter_out_intervals. This allows the user to disable the
  !        header from being rewritten after the first time it is written.
  logical(lk), public, save :: disable_header_rewrites = fals
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
  logical(lk), public, save :: output_vorticity               = fals
  logical(lk), public, save :: output_cylindrical_velocity    = fals
  logical(lk), public, save :: output_qcriterion              = fals
  logical(lk), public, save :: output_conserved_variables     = fals
  logical(lk), public, save :: output_temperature             = fals
  logical(lk), public, save :: output_thermo_variables        = fals
  logical(lk), public, save :: output_gradients               = fals
  logical(lk), public, save :: output_transport_variables     = fals
  logical(lk), public, save :: output_stagnation_temperature  = fals
  logical(lk), public, save :: output_stagnation_pressure     = fals
  logical(lk), public, save :: output_stagnation_variables    = fals
  logical(lk), public, save :: output_marker_variables        = fals
  logical(lk), public, save :: output_species_mass_fractions  = true
  logical(lk), public, save :: output_species_mole_fractions  = fals
  logical(lk), public, save :: output_extra_species_variables = fals
  logical(lk), public, save :: output_debug_species_variables = fals
  logical(lk), public, save :: output_grdnt_species_variables = fals
  !
  ! GOVERNING EQUATIONS
  integer, public, save :: governing_equations = NavierStokes_Eqns
  integer, public, save :: turbulence_model    = laminar_flow
  !
  ! MOVING REFERENCE FRAME
  type(mrf_t), public, save :: moving_ref_frame
  !
  ! CHEMICAL REACTION SOURCE TERMS
  logical(lk), public, save :: disable_chemical_reactions         = fals
  integer,     public, save :: semi_implicit_chemistry            = 0
  logical(lk), public, save :: traceback_unconverged_newton_temps = fals
 !integer,     public, save :: reactions_start_iter               = 0
  !
  ! IGNITION CONTROL SOURCE TERMS
  integer,  public, save :: &
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
  real(wp), public, save :: alpha_aoaref = 0.0_wp
  real(wp), public, save :: beta_aoaref  = 0.0_wp
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
  ! LIMITER FILTER OPTIONS
  integer,  public, save :: Limiter_Filter_Option = 0
  real(wp), public, save :: Limiter_Filter_etac   = 0.0_wp
  integer,  public, save :: Limiter_Filter_Order  = 16
  real(wp), public, save :: Limiter_Filter_Alpha  = 36.0_wp
  real(wp), public, save :: Limiter_Filter_eps    = 0.0_wp
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
  !
  ! LIMITER MARKER SENSITIVITY OPTIONS
  real(wp), public, save :: Limiter_Marker_Min  =  1.0_wp
  real(wp), public, save :: Limiter_Marker_Max  =  1.0_wp
  real(wp), public, save :: Limiter_Marker_c1   = -1.0_wp
  real(wp), public, save :: Limiter_Marker_c2   = -1.0_wp
  real(wp), public, save :: Limiter_Marker_dpsi = -1.0_wp
  !
  ! CONTROL OVER FREQUENCY OF FILTERING/LIMITING
  integer, public, save :: Filtering_Interval = 0
  integer, public, save :: Limiting_Interval  = 0
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
  !
  ! FLAME SPEED TEST PARAMETERS
  integer, public, save :: flame_speed_test_dof = 0
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
  !
  ! CGNS INPUT OPTIONS
  logical(lk), public, save :: correct_cgns_grid_points = fals
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
  !
  ! #################################
  ! ####   Non-Input Variables   ####
  ! #################################
  !
  real(wp), public, parameter :: Default_Marker_dpsi = 0.5_wp
  real(wp), public, parameter :: Default_Marker_c1   = 4.0_wp
  real(wp), public, parameter :: Default_Marker_c2   = 4.25_wp
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
  real(wp), public, save :: marker_c1   = Default_Marker_c1
  real(wp), public, save :: marker_c2   = Default_Marker_c2
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
  ! Turn quadrature/matrix output on by setting matrix_output_unit to anything
  ! besides iout (which is a parameter variable set to the iso_fortran_env
  ! constant output_unit, which historically has the value 6)
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
  ! Interface to check if two boundary condition input derived types are equal
  !
  interface operator(==)
    module procedure bc_input_t_equals, &
                     relaxation_t_equals, &
                     sem_bc_t_equals, &
                     periodic_bc_t_equals
  end interface
  !
  ! Interface to check if two boundary condition input derived types are equal
  !
  interface operator(/=)
    module procedure bc_input_t_not_equal, &
                     relaxation_t_not_equal, &
                     sem_bc_t_not_equal, &
                     periodic_bc_t_not_equal
  end interface
  !
  !$ACC DECLARE CREATE(walls_are_exact,sub_inflow_method, &
  !$ACC&               turb_is_on,include_tke_in_energy, &
  !$ACC&               time,itcur,gam,pr,aref,rhoref,tref,tmin,tmax, &
  !$ACC&               cpref_nd,rgasref_nd, &
  !$ACC&               suth_Tref,suth_Sref,suth_muref,suth_c1, &
  !$ACC&               LDG_beta,LDG_tau,visc_variables_to_ave, &
  !$ACC&               ideal_total_temperature,limit_total_temperature, &
  !$ACC&               invs_flux_method,Resolution_Marker_Variable)
contains
!
!###############################################################################
!
subroutine move_bcin_to_device
  !
  !.. Local Scalars ..
  integer :: i,i1,i2,m1,m2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "move_bcin_to_deviec"
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
!
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
!
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
!
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
!
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
!
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
!
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
!
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
!
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
!
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
!
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
!
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
!
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
!
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
  in_rrf = any(abs(this%rotation) >= one*rrf_dim)
  in_trf = any(abs(this%translation) >= one*trf_dim)
  !
  return_value = (in_rrf .or. in_trf)
  !
end function mrf_is_active
!
!###############################################################################
!
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
!
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
!
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
!
pure &
function if_interpto_newgrid() result(return_value)
  !
  ! Function to determine whether code will interpolate onto new grid.
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
subroutine adjust_cfl_cycles
  !
continue
  !
  cfl_cycle_start = itcur
  cfl_cycle_end = itcur + cfl_cycles
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
!
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
  use ovar,      only : grid_format, gridfile, plot3d_all_strides
  use ovar,      only : plot3d_i_stride, plot3d_j_stride, plot3d_k_stride
  use ovar,      only : plot3d_agglomerate_order,parametric_agglomeration
  use ovar,      only : plot3d_cutfile, plot3d_bcsfile, grid_scaling_factor
  use ovar,      only : internal_grid_order
  use ovar,      only : interpto_grid_format, interpto_gridfile
  use ovar,      only : interpto_plot3d_all_strides, interpto_plot3d_i_stride
  use ovar,      only : interpto_plot3d_j_stride, interpto_plot3d_k_stride
  use ovar,      only : interpto_plot3d_agglomerate_order
  use ovar,      only : interpto_plot3d_cutfile, interpto_plot3d_bcsfile
  use ovar,      only : interpto_grid_scaling_factor
  use ovar,      only : solution_format, solution_file
  use ovar,      only : loc_solution_pts, loc_flux_pts
  use ovar,      only : loc_triangle_pts, loc_tetrahedra_pts
  use ovar,      only : loc_quadrature_pts
  use ovar,      only : combine_common_flux_projections
  use ovar,      only : flux_divergence, correction_function
  use ovar,      only : Runge_Kutta_Scheme, num_rk_stages
  use ovar,      only : invs_flux_method, visc_flux_method
  use ovar,      only : visc_variables_to_ave, walls_are_exact
  use ovar,      only : num_timesteps, Final_Time, constant_dt
  use ovar,      only : cfl, cfl_beg, cfl_end, cfl_cycles, cfl_start_iter
  use ovar,      only : Timestep_Type, minimum_timesteps
  use ovar,      only : all_points, metis_option, use_old_restart_format
  use ovar,      only : mms_opt, mms_init, mms_output_solution
  use ovar,      only : itestcase, jet_init_func, governing_equations
  use ovar,      only : moving_ref_frame
  use ovar,      only : disable_chemical_reactions, semi_implicit_chemistry
  use ovar,      only : traceback_unconverged_newton_temps
 !use ovar,      only : reactions_start_iter
  use ovar,      only : ignition_start_timestep
  use ovar,      only : ignition_duration
  use ovar,      only : ignition_prod2reac_L2_ratio_goal
  use ovar,      only : ignition_source_magnitude
  use ovar,      only : ignition_max_temperature
  use ovar,      only : ignition_spark_cap
  use ovar,      only : ignition_source_radius
  use ovar,      only : ignition_source_location
  use ovar,      only : load_restart_file, restart_interval, restart_file
  use ovar,      only : restart_output_uses_host_roots
  use ovar,      only : CFL_reduction_for_NaN_restart
  use ovar,      only : num_saved_time_ave_restarts
  use ovar,      only : restart_frequency
  use ovar,      only : iter_out_interval, results_interval
  use ovar,      only : output_eperk_cfl, disable_header_rewrites
  use ovar,      only : output_summary_file
  use ovar,      only : output_dir, output_interval, output_frequency
  use ovar,      only : continuous_output
  use ovar,      only : postproc_debug_grid
  use ovar,      only : output_bnd_profiles,sub_inflow_method
  use ovar,      only : machref, reyref, gam, Pr, Prt
  use ovar,      only : schmidt_lam, schmidt_turb
  use ovar,      only : rhoref, tref, pref, rgasref
  use ovar,      only : turbulence_model, tinref, tvrref
  use ovar,      only : alpha_aoaref, beta_aoaref
  use ovar,      only : profile_io_all,profile_io_restart
  use ovar,      only : profile_io_cgns,profile_io_cgns_subs
  use ovar,      only : suth_muref, suth_Tref, suth_Sref
  use ovar,      only : Filter_Option, Limiter_Filter_Option
  use ovar,      only : Filter_etac,   Limiter_Filter_etac
  use ovar,      only : Filter_Order,  Limiter_Filter_Order
  use ovar,      only : Filter_Alpha,  Limiter_Filter_Alpha
  use ovar,      only : Filter_eps,    Limiter_Filter_eps
  use ovar,      only : use_sponge_filter
  use ovar,      only : Sponge_Filter_etac
  use ovar,      only : Sponge_Filter_Order
  use ovar,      only : Sponge_Filter_Alpha
  use ovar,      only : Sponge_Beyond_X
  use ovar,      only : Limiter_Option, Resolution_Marker_Variable
  use ovar,      only : Limiter_Marker_Min, Limiter_Marker_Max
  use ovar,      only : Limiter_Marker_dpsi
  use ovar,      only : Limiter_Marker_c1, Limiter_Marker_c2
  use ovar,      only : Filtering_Interval, Limiting_Interval
  use ovar,      only : loc_output_pts, loc_init_pts, loc_error_pts
  use ovar,      only : vortex_bigR, vortex_bigGam
  use ovar,      only : output_plot3d_to_tecplot, output_bface_array
  use ovar,      only : interpolate_before_output_variable
  use ovar,      only : cgns_high_order_output,cgns_mixed_output
  use ovar,      only : correct_cgns_grid_points
  use ovar,      only : convergence_order_abs_res, convergence_order_max_res
  use ovar,      only : dump_memory_usage,use_bc_conflicts
  use ovar,      only : dump_fluxes,check_flux_point_coordinates
  use ovar,      only : use_unrolled_dot_products
  use ovar,      only : VortexGrid_DOF,VortexGrid_random_perturb
  use ovar,      only : VortexGrid_perturb_percent,VortexGrid_Lref
  use ovar,      only : VortexGrid_Hex_to_Pyramids,VortexGrid_is_3D
  use ovar,      only : TCF_Inner_Radius,TCF_Inner_Wall_Spacing
  use ovar,      only : TCF_Outer_Radius,TCF_Outer_Wall_Spacing
  use ovar,      only : TCF_Aspect_Ratio,TCF_Axial_Cells
  use ovar,      only : TCF_Radial_Cells,TCF_Azimuthal_Cells
  use ovar,      only : flame_speed_test_dof
  use ovar,      only : channel_body_force,channel_init_file
  use ovar,      only : channel_spatial_ave_restart
  use ovar,      only : channel_nondim_output_coord
  use ovar,      only : channel_yplus_of_one, channel_rake_save_interval
  use ovar,      only : channel_grid_ave_xplus, channel_grid_ave_zplus
  use ovar,      only : channel_grid_ave_yplus, channel_grid_yplus_wall
  use ovar,      only : convert_restart_to_cgns, completely_disable_cgns
  use ovar,      only : delete_converted_restart_files
  use ovar,      only : lustre_stripe_count, lustre_stripe_size
  use ovar,      only : bc_input!,bc_input_t
  use ovar,      only : output_time_averaging, time_ave_file
  use ovar,      only : start_averaging_timestep, start_averaging_tstar
  use ovar,      only : time_ave_vel_is_axisymm
  use ovar,      only : time_scaling_factor
  use ovar,      only : spatial_ave_restart_x_and_z_planes
  use ovar,      only : spatial_ave_restart_z_planes
  use ovar,      only : time_average_restart_files, restart_list_file
  use ovar,      only : use_async_for_error
  use ovar,      only : ideal_total_temperature,limit_total_temperature
  use ovar,      only : output_limiter_statistics,output_filter_diagonals
  use ovar,      only : output_vorticity
  use ovar,      only : output_cylindrical_velocity
  use ovar,      only : output_qcriterion
  use ovar,      only : output_conserved_variables
  use ovar,      only : output_temperature
  use ovar,      only : output_thermo_variables
  use ovar,      only : output_gradients
  use ovar,      only : output_transport_variables
  use ovar,      only : output_stagnation_temperature
  use ovar,      only : output_stagnation_pressure
  use ovar,      only : output_stagnation_variables
  use ovar,      only : output_marker_variables
  use ovar,      only : output_species_mass_fractions
  use ovar,      only : output_species_mole_fractions
  use ovar,      only : output_extra_species_variables
  use ovar,      only : output_debug_species_variables
  use ovar,      only : output_grdnt_species_variables
  use ovar,      only : use_improved_SEM_normalization
  use ovar,      only : use_SEM_profile_for_initialization
  use ovar,      only : SEM_eddy_volume_is_spherical
  use ovar,      only : SEM_box_num_eddies_multiplier
  use ovar,      only : SEM_use_random_seed_from_restart
  use ovar,      only : average_sem_inflow_plane
  use ovar,      only : force_nd_subroutines, newton_rlx, conv_temp
  use ovar,      only : chemistry_sub_iterations
  use ovar,      only : ignore_thermo_newton_warnings
  use module_cgns_types, only : cgns_use_queue
  !
  implicit none
  !
  ! Set everything in the module to public
  !
  public
  !
  integer, save :: all_orders = -3
  integer, save :: solution_order = 3
  integer, save :: projection_order = -3
  integer, save :: quadrature_order = -3
  integer, save :: output_order = -3
  integer, save :: error_order = -3
  !
  integer, save :: interpto_all_orders = -3
  integer, save :: interpto_solution_order = 3
  integer, save :: interpto_projection_order = -3
  integer, save :: interpto_quadrature_order = -3
  integer, save :: interpto_output_order = -3
  integer, save :: interpto_error_order = -3
  !
  integer, save :: dmat_grad_temp_order = -1
  !
  ! ########################
  ! ###  Input Namelist  ###
  ! ########################
  !
  ! GRID INFORMATION
  namelist / input / grid_format, gridfile, plot3d_cutfile, plot3d_bcsfile, &
                     grid_scaling_factor, internal_grid_order
  ! PLOT3D COARSENING OPTIONS
  namelist / input / plot3d_i_stride, plot3d_j_stride, plot3d_k_stride, &
                     plot3d_all_strides, plot3d_agglomerate_order, &
                     parametric_agglomeration
  ! INTERPOLATION UTILITY GRID INFORMATION
  namelist / input / interpto_grid_format, interpto_gridfile, &
                     interpto_plot3d_cutfile, interpto_plot3d_bcsfile, &
                     interpto_grid_scaling_factor
  ! INTERPOLATION UTILITY PLOT3D COARSENING OPTIONS
  namelist / input / interpto_plot3d_i_stride, interpto_plot3d_j_stride, &
                     interpto_plot3d_k_stride, interpto_plot3d_all_strides, &
                     interpto_plot3d_agglomerate_order
  ! INTERPOLATION UTILITY SOLUTION ORDER
  namelist / input / interpto_solution_order, interpto_all_orders, &
                     interpto_projection_order, interpto_quadrature_order, &
                     interpto_output_order, interpto_error_order
  ! BASE FR OPTIONS
  namelist / input / solution_order, flux_divergence, correction_function
  ! LOCATION OF SOLUTION, FLUX, AND TRIANGULAR POINTS
  namelist / input / loc_solution_pts, loc_flux_pts, loc_triangle_pts, &
                     loc_tetrahedra_pts
  ! LOCATION OF INITIALIZATION POINTS
  namelist / input / loc_init_pts
  ! RUNGE-KUTTA OPTIONS
  namelist / input / Runge_Kutta_Scheme, num_rk_stages
  ! COMMON FLUX OPTIONS
  namelist / input / invs_flux_method, visc_flux_method, visc_variables_to_ave
  ! TIME STEPPING INFORMATION
  namelist / input / num_timesteps, Final_Time, constant_dt, &
                     Timestep_Type, minimum_timesteps
  ! CFL CONTROL VARIABLES
  namelist / input / cfl, cfl_beg, cfl_end, cfl_cycles, cfl_start_iter
  ! OVER-INTEGRATION OPTIONS
  namelist / input / projection_order, quadrature_order, loc_quadrature_pts, &
                     combine_common_flux_projections
  ! ERROR-INTEGRATION OPTIONS
  namelist / input / error_order, loc_error_pts
  ! FORCE GLOBAL ORDER AND SOLUTION POINTS
  namelist / input / all_orders, all_points
  ! PARALLEL OPTIONS
  namelist / input / metis_option
  ! METHOD OF MANUFACTURED SOLUTIONS (MMS) OPTIONS
  namelist / input / mms_opt, mms_init, mms_output_solution
  ! INITIALIZATION / TEST CASE
  namelist / input / itestcase, jet_init_func
  ! RESTART OPTIONS
  namelist / input / load_restart_file, restart_interval, &
                     use_old_restart_format, restart_file, &
                     restart_output_uses_host_roots, &
                     CFL_reduction_for_NaN_restart, &
                     num_saved_time_ave_restarts, &
                     restart_frequency
  ! OPTIONS TO USE SOLUTION FILE FOR INITIAL CONDITIONS
  namelist / input / solution_format, solution_file
  ! OPTIONS FOR SEM
  namelist / input / use_improved_SEM_normalization, &
                     use_SEM_profile_for_initialization, &
                     SEM_eddy_volume_is_spherical, &
                     SEM_box_num_eddies_multiplier, &
                     SEM_use_random_seed_from_restart
  ! FLAG TO ACTIVATE A TEST OF THE SEM CODE
  namelist / input / average_sem_inflow_plane
  ! STATISTICAL OUTPUT OPTIONS
  namelist / input / iter_out_interval, results_interval, output_eperk_cfl, &
                     disable_header_rewrites, output_summary_file
  ! SOLUTION OUTPUT OPTIONS
  namelist / input / output_dir, output_interval, continuous_output, &
                     output_order, loc_output_pts, output_frequency
  ! SELECTION OF OUTPUT VARIABLES IN THE CGNS SOLUTION FILES
  namelist / input / output_vorticity, output_conserved_variables, &
                     output_cylindrical_velocity, &
                     output_temperature, output_thermo_variables, &
                     output_gradients, output_transport_variables, &
                     output_stagnation_temperature, &
                     output_stagnation_pressure, &
                     output_stagnation_variables, output_marker_variables, &
                     output_qcriterion, &
                     output_species_mass_fractions, &
                     output_species_mole_fractions, &
                     output_extra_species_variables, &
                     output_debug_species_variables, &
                     output_grdnt_species_variables
  ! GOVERNING EQUATIONS
  namelist / input / governing_equations, turbulence_model
  ! MOVING REFERENCE FRAME
  namelist / input / moving_ref_frame
  ! CHEMICAL REACTION SOURCE TERMS
  namelist / input / disable_chemical_reactions, semi_implicit_chemistry, &
                     traceback_unconverged_newton_temps
                    !traceback_unconverged_newton_temps, &
                    !reactions_start_iter
  !
  ! IGNITION CONTROL SOURCE TERMS
  namelist / input / ignition_start_timestep, &
                     ignition_duration, &
                     ignition_prod2reac_L2_ratio_goal, &
                     ignition_max_temperature, &
                     ignition_spark_cap, &
                     ignition_source_magnitude, &
                     ignition_source_radius, &
                     ignition_source_location
  ! NONDIMENSIONAL REFERENCE CONDITIONS
  namelist / input / machref, reyref, gam
  ! REALIZABLE TOTAL TEMPERATURE LIMITER
  namelist / input / ideal_total_temperature, limit_total_temperature
  ! REFERENCE CONDITIONS
  namelist / input / rhoref, tref, pref, rgasref
  ! REFERENCE TURBULENCE CONDITIONS
  namelist / input / tinref, tvrref
  ! REFERENCE ANGLES OF ATTACK
  namelist / input / alpha_aoaref, beta_aoaref
  ! LAMINAR AND TURBULENT PRANDTL NUMBERS
  namelist / input / Pr, Prt
  ! LAMINAR AND TURBULENT SCHMIDT NUMBERS
  namelist / input / schmidt_lam, schmidt_turb
  ! SUTHERLANDS LAW COEFFICIENTS
  namelist / input / suth_muref, suth_Tref, suth_Sref
  ! FILTER OPTIONS
  namelist / input / Filter_Option, Filter_etac, Filter_Order, &
                     Filter_Alpha, Filter_eps
  ! LIMITER FILTER OPTIONS
  namelist / input / Limiter_Filter_Option, Limiter_Filter_etac, &
                     Limiter_Filter_Order, Limiter_Filter_Alpha, &
                     Limiter_Filter_eps
  ! SPONGE FILTER OPTIONS
  namelist / input / use_sponge_filter, &
                     Sponge_Filter_etac, &
                     Sponge_Filter_Order, &
                     Sponge_Filter_Alpha, &
                     Sponge_Beyond_X
  ! LIMITER OPTIONS
  namelist / input / Limiter_Option, Resolution_Marker_Variable
  ! LIMITER MARKER SENSITIVITY OPTIONS
  namelist / input / Limiter_Marker_Min, Limiter_Marker_Max, &
                     Limiter_Marker_c1, Limiter_Marker_c2, &
                     Limiter_Marker_dpsi
  ! CONTROL OVER FREQUENCY OF FILTERING/LIMITING
  namelist / input / Filtering_Interval, Limiting_Interval
  ! OUTPUT THE STATISTICS OF THE LIMITER
  namelist / input / output_limiter_statistics
  ! OPTION TO OUTPUT FILTER DIAGONALS TO FILE
  namelist / input / output_filter_diagonals
  ! ISENTROPIC VORTEX CONSTANTS
  namelist / input / vortex_bigR, vortex_bigGam, VortexGrid_Lref
  ! TAYLOR-GREEN VORTEX CONSTANTS
  namelist / input / VortexGrid_DOF, VortexGrid_random_perturb, &
                     VortexGrid_perturb_percent, VortexGrid_Hex_to_Pyramids, &
                     VortexGrid_is_3D
  ! TAYLOR-COUETTE FLOW CONSTANTS
  namelist / input / TCF_Inner_Radius, TCF_Inner_Wall_Spacing, &
                     TCF_Outer_Radius, TCF_Outer_Wall_Spacing, &
                     TCF_Aspect_Ratio, TCF_Axial_Cells, &
                     TCF_Radial_Cells, TCF_Azimuthal_Cells
  ! FLAME SPEED TEST PARAMETERS
  namelist / input / flame_speed_test_dof
  ! CHANNEL FLOW PARAMETERS
  namelist / input / channel_body_force, channel_init_file, &
                     channel_spatial_ave_restart, channel_nondim_output_coord, &
                     channel_grid_ave_xplus, channel_grid_ave_zplus, &
                     channel_grid_ave_yplus, channel_grid_yplus_wall, &
                     channel_yplus_of_one, channel_rake_save_interval
  ! OPTIONS FOR MODIFYING BC ALGORITHMS
  namelist / input / walls_are_exact, sub_inflow_method, use_bc_conflicts
  ! CGNS OUTPUT OPTIONS
  namelist / input / cgns_use_queue, interpolate_before_output_variable, &
                     cgns_high_order_output, cgns_mixed_output, &
                     convert_restart_to_cgns, completely_disable_cgns, &
                     delete_converted_restart_files
  ! CGNS INPUT OPTIONS
  namelist / input / correct_cgns_grid_points
  ! CONVERGENCE CRITERIA
  namelist / input / convergence_order_abs_res, convergence_order_max_res
  ! PROFILE PARALLEL IO
  namelist / input / profile_io_all, profile_io_restart, &
                     profile_io_cgns, profile_io_cgns_subs
  ! OPTIONS TO OUTPUT BOUNDARY PROFILE
  namelist / input / output_bnd_profiles
  ! CONTROL DOT_PRODUCT PROCEDURES BOUND TO MATRIX DERIVED-TYPE
  namelist / input / use_unrolled_dot_products
  ! LUSTRE FILE SYSTEM STRIPING PARAMETERS
  namelist / input / lustre_stripe_count, lustre_stripe_size
  ! OPENACC INPUT OPTIONS
  namelist / input / use_async_for_error
  ! TIME AVERAGING
  namelist / input / output_time_averaging, time_ave_file, &
                     time_scaling_factor, time_ave_vel_is_axisymm, &
                     start_averaging_timestep, start_averaging_tstar
  ! SPATIAL AVERAGE RESTART FILES
  namelist / input / spatial_ave_restart_x_and_z_planes, &
                     spatial_ave_restart_z_planes
  ! UTILITY TO TIME AVERAGE EXISTING RESTART FILES
  namelist / input / time_average_restart_files, restart_list_file
  ! MEMORY USAGE
  namelist / input / dump_memory_usage
  ! DUMP FLUX INFORMATION FOR DEBUGGING
  namelist / input / dump_fluxes, check_flux_point_coordinates
  ! ACTIVATE SPECIAL GRID FOR DEBUGGING THE POSTPROCESSOR
  namelist / input / postproc_debug_grid
  ! OPTION TO OUTPUT PLOT3D GRID TO TECPLOT IN UNSTRUCTURED FORMAT
  namelist / input / output_plot3d_to_tecplot
  ! OPTION TO OUTPUT CONTENTS OF BFACE ARRAY
  namelist / input / output_bface_array
  !
  namelist / input / force_nd_subroutines, dmat_grad_temp_order, &
                     newton_rlx, conv_temp, &
                     chemistry_sub_iterations, ignore_thermo_newton_warnings
  ! INPUT BOUNDARY CONDITIONS
  namelist / input / bc_input
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
subroutine allocate_namelist_variables(opt_n_bc_input,opt_char_len, &
                                       opt_n_species)
  !
  !.. Use Statements ..
  use module_kind_types, only : zero,alloc_error,error_message
  use module_kind_types, only : debug_timer
  use module_kind_types, only : entering_procedure
  use module_kind_types, only : leaving_procedure
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: opt_n_bc_input
  integer, optional, intent(in) :: opt_char_len
  integer, optional, intent(in) :: opt_n_species
  !
  !.. Local Scalars ..
  integer :: n,ierr
  integer :: default_len,n_species,n_bc_input
  character(len=20) :: blanks
  character(len=50) :: array_name
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "allocate_namelist_variables"
  integer, parameter :: max_default_len = 150
  integer, parameter :: max_n_species = 7
  integer, parameter :: max_n_bc_input = 20
  !
continue
  !
  n_bc_input = max_n_bc_input
  n_species = max_n_species
  default_len = max_default_len
  !
  if (present(opt_n_bc_input)) n_bc_input = opt_n_bc_input
  if (present(opt_n_species)) n_species = opt_n_species
  if (present(opt_char_len)) default_len = opt_char_len
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Make sure everything is deallocated first
  !
  if (allocated(bc_input)) then
    deallocate ( bc_input , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bc_input",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(gridfile)) then
    deallocate ( gridfile , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"gridfile",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(plot3d_cutfile)) then
    deallocate ( plot3d_cutfile , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"plot3d_cutfile",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(plot3d_bcsfile)) then
    deallocate ( plot3d_bcsfile , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"plot3d_bcsfile",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(interpto_gridfile)) then
    deallocate ( interpto_gridfile , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interpto_gridfile",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(interpto_plot3d_cutfile)) then
    deallocate ( interpto_plot3d_cutfile , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interpto_plot3d_cutfile",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(interpto_plot3d_bcsfile)) then
    deallocate ( interpto_plot3d_bcsfile , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interpto_plot3d_bcsfile",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(restart_file)) then
    deallocate ( restart_file , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"restart_file",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(solution_file)) then
    deallocate ( solution_file , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"solution_file",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(output_dir)) then
    deallocate ( output_dir , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"output_dir",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(channel_init_file)) then
    deallocate ( channel_init_file , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"channel_init_file",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(lustre_stripe_size)) then
    deallocate ( lustre_stripe_size , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"lustre_stripe_size",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(time_ave_file)) then
    deallocate ( time_ave_file , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"time_ave_file",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(restart_list_file)) then
    deallocate ( restart_list_file , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"restart_list_file",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  ! Allocate bc_input before reading the namelist
  !
  allocate ( bc_input(1:n_bc_input) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bc_input",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate the mass_frac and species_name components
  ! of bc_input before reading the namelist.
  !
#ifdef ALLOCATABLE_MASS_FRAC
  do n = 1,size(bc_input)
    allocate ( bc_input(n)%mass_frac(1:n_species) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,2) n,"mass_frac"
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    allocate ( bc_input(n)%species_name(1:n_species) , source="" , &
               stat=ierr , errmsg=error_message )
    write (array_name,2) n,"species_name"
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  end do
#endif
  !
  ! Allocate all the character input variables
  !
  allocate ( character(len=default_len) :: gridfile , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"gridfile",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( character(len=default_len) :: plot3d_cutfile , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"plot3d_cutfile",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( character(len=default_len) :: plot3d_bcsfile , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"plot3d_bcsfile",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( character(len=default_len) :: interpto_gridfile , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"interpto_gridfile",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( character(len=default_len) :: interpto_plot3d_cutfile , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"interpto_plot3d_cutfile",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( character(len=default_len) :: interpto_plot3d_bcsfile , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"interpto_plot3d_bcsfile",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( character(len=default_len) :: restart_file , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"restart_file",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( character(len=default_len) :: solution_file , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"solution_file",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( character(len=default_len) :: output_dir , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"output_dir",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( character(len=default_len) :: channel_init_file , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"channel_init_file",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( character(len=default_len) :: lustre_stripe_size , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"lustre_stripe_size",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( character(len=default_len) :: time_ave_file , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"time_ave_file",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( character(len=default_len) :: restart_list_file , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"restart_list_file",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  ! Set all the character input variables to their default values
  !
  if (default_len > 0) then
    !
    write (blanks,1) default_len
    !
    write (gridfile,blanks)
    write (plot3d_cutfile,blanks)
    write (plot3d_bcsfile,blanks)
    write (interpto_gridfile,blanks)
    write (interpto_plot3d_cutfile,blanks)
    write (interpto_plot3d_bcsfile,blanks)
    write (restart_file,blanks)
    write (solution_file,blanks)
    write (output_dir,blanks)
    write (channel_init_file,blanks)
    write (lustre_stripe_size,blanks)
    write (time_ave_file,blanks)
    write (restart_list_file,blanks)
    !
  end if
  !
  if (len(restart_file) >= 7) then
    restart_file(1:7) = "out.rst"
  else
    restart_file = "out.rst" ! F2003 AUTO-REALLOC
  end if
  !
  if (len(output_dir) >= 1) then
    output_dir(1:1) = "."
  else
    output_dir = "." ! F2003 AUTO-REALLOC
  end if
  !
  if (len(lustre_stripe_size) >= 2) then
    lustre_stripe_size(1:2) = "4m"
  else
    lustre_stripe_size = "4m" ! F2003 AUTO-REALLOC
  end if
  !
  if (len(time_ave_file) >= 7) then
    time_ave_file(1:7) = "out.ave"
  else
    time_ave_file = "out.ave" ! F2003 AUTO-REALLOC
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  1 format ('(',i0,'x)')
  2 format ("bc_input(",i0,")%",a)
  !
end subroutine allocate_namelist_variables
!
!###############################################################################
!
subroutine resize_namelist_allocatables
  !
  !.. Use Statements ..
  use module_kind_types, only : alloc_error,error_message
  use module_kind_types, only : debug_timer
  use module_kind_types, only : entering_procedure
  use module_kind_types, only : leaving_procedure
  use ovar,              only : bc_input_t
  !
  !.. Local Scalars ..
  integer :: i,n,ierr
  !
  !.. Local Allocatable Arrays ..
  type(bc_input_t), allocatable :: tmp_bc_input(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "resize_namelist_allocatables"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Reallocate bc_input to the actual size needed for the inputs
  !
  n = 0
  do i = 1,size(bc_input)
    if (bc_input(i)%is_not_default()) then
      n = n + 1
    end if
  end do
  !
  allocate ( tmp_bc_input(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_bc_input",1,__LINE__,__FILE__,ierr,error_message)
  !
  n = 0
  do i = 1,size(bc_input)
    if (bc_input(i)%is_not_default()) then
      n = n + 1
      tmp_bc_input(n) = bc_input(i)
    end if
  end do
  !
  call move_alloc( from=tmp_bc_input , to=bc_input )
  !
  do n = 1,size(bc_input)
    bc_input(n)%sem%profile_file = trim(adjustl(bc_input(n)%sem%profile_file))
  end do
  !
  ! Following are F2003 AUTO REALLOCATIONS
  !
  gridfile = trim(adjustl(gridfile))
  plot3d_cutfile = trim(adjustl(plot3d_cutfile))
  plot3d_bcsfile = trim(adjustl(plot3d_bcsfile))
  interpto_gridfile = trim(adjustl(interpto_gridfile))
  interpto_plot3d_cutfile = trim(adjustl(interpto_plot3d_cutfile))
  interpto_plot3d_bcsfile = trim(adjustl(interpto_plot3d_bcsfile))
  restart_file = trim(adjustl(restart_file))
  solution_file = trim(adjustl(solution_file))
  output_dir = trim(adjustl(output_dir))
  channel_init_file = trim(adjustl(channel_init_file))
  lustre_stripe_size = trim(adjustl(lustre_stripe_size))
  time_ave_file = trim(adjustl(time_ave_file))
  restart_list_file = trim(adjustl(restart_list_file))
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine resize_namelist_allocatables
!
!###############################################################################
!
end module input_namelist_mod
!
!###############################################################################
!
