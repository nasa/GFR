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
  use module_kind_types, only : iout,zero,one,eps6
  use module_kind_types, only : mpi_wtime
  use module_kind_types, only : Correction_gDG
  use module_kind_types, only : bc_unknown
  use module_kind_types, only : project_derivatives
  use module_kind_types, only : derivative_of_proj_fluxes
  use module_kind_types, only : BarycentricLobatto_TriPoints
  use module_kind_types, only : LagrangePolynomial
  use module_kind_types, only : Global_Cell_Timestep
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
  type, public :: relaxation_t
    real(wp) :: value = one
    real(wp) :: time  = zero
    integer  :: iter  = 0
  end type relaxation_t
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
  contains
    !
    procedure :: is_not_default
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
  type, public :: bc_in_t
    real(wp), allocatable :: cv(:)
    real(wp), allocatable :: pv(:)
    real(wp), allocatable :: gp(:)
    real(wp) :: wall_temperature = -one
    type(relaxation_t) :: relax
    integer :: bc_input_idx
    integer :: bc_type = bc_unknown
  end type bc_in_t
  !
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
  ! GRID INFORMATION
  integer,            public, save :: grid_format = Internally_Generated_Grid
  character(len=:), public, save, allocatable :: gridfile
  character(len=:), public, save, allocatable :: plot3d_cutfile
  character(len=:), public, save, allocatable :: plot3d_bcsfile
  real(wp),           public, save :: grid_scaling_factor = 1.0_wp
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
  integer,            public, save :: interpto_grid_format      = Unknown_Format
  character(len=:), public, save, allocatable :: interpto_gridfile
  character(len=:), public, save, allocatable :: interpto_plot3d_cutfile
  character(len=:), public, save, allocatable :: interpto_plot3d_bcsfile
  real(wp),           public, save :: interpto_grid_scaling_factor = 1.0_wp
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
  integer, public, save :: loc_solution_pts = Legendre_Gauss
  integer, public, save :: loc_flux_pts     = Legendre_Gauss
  integer, public, save :: loc_triangle_pts = BarycentricLobatto_TriPoints
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
  logical(lk),        public, save :: load_restart_file              = fals
  integer,            public, save :: restart_interval               = 0
  logical(lk),        public, save :: use_old_restart_format         = fals
  character(len=:), public, save, allocatable :: restart_file
  logical(lk),        public, save :: restart_output_uses_host_roots = true
  !
  ! OPTIONS TO USE SOLUTION FILE FOR INITIAL CONDITIONS
  integer,            public, save :: solution_format = Unknown_Format
  character(len=:), public, save, allocatable :: solution_file
  !
  ! PROFILE FILE FOR SEM TURBULENCE INFLOW GENERATOR
  character(len=:), public, save, allocatable :: sem_profile_file
  !
  ! STATISTICAL OUTPUT OPTIONS
  integer, public, save :: iter_out_interval = 1
  integer, public, save :: results_interval  = 0
  !
  ! SOLUTION OUTPUT OPTIONS
  character(len=:), public, save, allocatable :: output_dir
  integer,            public, save :: output_interval   = -5
  real(wp),           public, save :: output_frequency  = -1.0_wp
  logical(lk),        public, save :: continuous_output = fals
  integer,            public, save :: loc_output_pts    = Equi_Distant_Points
  !
  ! SELECTION OF OUTPUT VARIABLES IN THE CGNS SOLUTION FILES
  logical(lk), public, save :: output_vorticity               = fals
  logical(lk), public, save :: output_conserved_variables     = fals
  logical(lk), public, save :: output_temperature             = fals
  logical(lk), public, save :: output_thermo_variables        = fals
  logical(lk), public, save :: output_gradients               = fals
  logical(lk), public, save :: output_transport_variables     = fals
  logical(lk), public, save :: output_stagnation_variables    = fals
  logical(lk), public, save :: output_marker_variables        = fals
  logical(lk), public, save :: output_species_mass_fractions  = true
  logical(lk), public, save :: output_extra_species_variables = fals
  logical(lk), public, save :: output_debug_species_variables = fals
  logical(lk), public, save :: output_grdnt_species_variables = fals
  !
  ! GOVERNING EQUATIONS
  integer, public, save :: governing_equations = NavierStokes_Eqns
  integer, public, save :: turbulence_model    = laminar_flow
  !
  ! CHEMICAL REACTION SOURCE TERMS
  logical(lk), public, save :: disable_chemical_reactions         = fals
  integer,     public, save :: semi_implicit_chemistry            = 0
  logical(lk), public, save :: traceback_unconverged_newton_temps = fals
 !integer,     public, save :: reactions_start_iter               = 0
  !
  ! IGNITION CONTROL SOURCE TERMS
  integer,  public, save :: ignition_start_timestep = 0
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
  !
  ! CHANNEL FLOW PARAMETERS
  real(wp),           public, save :: channel_body_force          = 0.0_wp
  character(len=:), public, save, allocatable :: channel_init_file
  logical(lk),        public, save :: channel_spatial_ave_restart = fals
  !
  ! OPTIONS FOR MODIFYING BC ALGORITHMS
  logical(lk), public, save :: walls_are_exact   = true
  integer,     public, save :: sub_inflow_method = 1
  logical(lk), public, save :: use_bc_conflicts  = fals
  !
  ! CGNS OUTPUT OPTIONS
  logical(lk), public, save :: convert_restart_to_cgns = fals
  logical(lk), public, save :: completely_disable_cgns = fals
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
  integer,           public, save :: lustre_stripe_count = 64
  character(len=:), public, save, allocatable :: lustre_stripe_size
  !
  ! OPENACC INPUT OPTIONS
  logical(lk), public, save :: use_async_for_error = fals
  !
  ! TIME AVERAGING
  logical(lk),        public, save :: output_time_averaging    = fals
  character(len=:), public, save, allocatable :: time_ave_file
  real(wp),           public, save :: time_scaling_factor      = 1.0_wp
  logical(lk),        public, save :: time_ave_vel_is_axisymm  = fals
  integer,            public, save :: start_averaging_timestep = -1
  real(wp),           public, save :: start_averaging_tstar    = -one
  !
  ! UTILITY TO TIME AVERAGE EXISTING RESTART FILES
  logical(lk),        public, save :: time_average_restart_files = fals
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
  !
  logical(lk),  public, save :: write_TaylorGreen_8s_solution = fals
  logical(lk),  public, save :: write_TaylorGreen_full_solution = fals
  !
  logical(lk),  public, save :: convergence_reached = fals
  logical(lk),  public, save :: this_is_final_timestep = fals
  !
  real(wp),         save :: max_wall_time
  real(wp),         save :: beg_wall_time
  real(wp),         save :: end_wall_time
  real(wp),         save :: wall_time
  real(wp), public, save :: work_units
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
pure function is_not_default(this) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(bc_input_t), intent(in) :: this
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
continue
  !
  return_value = true
  !
  if (len_trim(adjustl(this%name)) > 0) return
  if (this%set_bc_value /= bc_unknown)  return
  !
  return_value = (.not. this%has_default_flow_conditions())
  !
end function is_not_default
!
!###############################################################################
!
pure function has_default_flow_conditions(this) result(return_value)
  !
  !.. Passed Argument from Invoking Object ..
  class(bc_input_t), intent(in) :: this
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
continue
  !
  return_value = fals
  !
  if (abs(this%t_static  +one) > eps6) return
  if (abs(this%p_static  +one) > eps6) return
  if (abs(this%rho_static+one) > eps6) return
  if (abs(this%mach      +one) > eps6) return
  if (abs(this%vx            ) > eps6) return
  if (abs(this%vy            ) > eps6) return
  if (abs(this%vz            ) > eps6) return
  if (abs(this%t_total   +one) > eps6) return
  if (abs(this%p_total   +one) > eps6) return
  if (abs(this%rho_total +one) > eps6) return
  if (abs(this%alpha_aoa     ) > eps6) return
  if (abs(this%beta_aoa      ) > eps6) return
  if (abs(this%wall_temp +one) > eps6) return
  !
  if (abs(this%relax%value- one) > eps6) return
  if (abs(this%relax%time -zero) > eps6) return
  if (    this%relax%iter        /=   0) return
  !
#ifdef ALLOCATABLE_MASS_FRAC
  if (allocated(this%mass_frac)) then
    if (any(this%mass_frac(:) > zero)) return
  end if
#else
  if (any(this%mass_frac(:) > zero)) return
#endif
  !
  return_value = true
  !
end function has_default_flow_conditions
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
subroutine gfr_cpu_time(time_flag)
  !
  !.. Use Statements ..
  use module_kind_types, only : gfr_timer
  !
  !.. Formal Arguments ..
  integer, intent(in) :: time_flag
  !
  !.. Local Scalars ..
  integer :: hours,minutes,seconds,total_sec
  !
continue
  !
  if (time_flag == 0) then
    !
    call gfr_timer( beg_wall_time )
    work_units = zero
    !
  else if (time_flag == 1) then
    !
    call gfr_timer( end_wall_time )
    wall_time = end_wall_time-beg_wall_time
    work_units = wall_time/one_work_unit
    !
  else if (time_flag == -1) then
    !
    total_sec = nint(wall_time)
    !
    hours   = total_sec / 3600
    minutes = (total_sec - hours*3600) / 60
    seconds = total_sec - hours*3600 - minutes*60
    !
    write (iout,1) wall_time,hours,minutes,seconds,work_units
    !
  end if
  !
  1 format (/,"TOTAL SIMULATION WALL TIME WAS ",f18.6," s.", &
              " (",i0,":",i2.2,":",i2.2,")",/, &
              "TOTAL COMPUTATIONAL COST WAS   ",f18.6," work units.")
  !
end subroutine gfr_cpu_time
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
  use ovar,      only : interpto_grid_format, interpto_gridfile
  use ovar,      only : interpto_plot3d_all_strides, interpto_plot3d_i_stride
  use ovar,      only : interpto_plot3d_j_stride, interpto_plot3d_k_stride
  use ovar,      only : interpto_plot3d_agglomerate_order
  use ovar,      only : interpto_plot3d_cutfile, interpto_plot3d_bcsfile
  use ovar,      only : interpto_grid_scaling_factor
  use ovar,      only : solution_format, solution_file
  use ovar,      only : loc_solution_pts, loc_flux_pts
  use ovar,      only : loc_triangle_pts, loc_quadrature_pts
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
  use ovar,      only : disable_chemical_reactions, semi_implicit_chemistry
  use ovar,      only : traceback_unconverged_newton_temps
 !use ovar,      only : reactions_start_iter
  use ovar,      only : ignition_start_timestep
  use ovar,      only : ignition_prod2reac_L2_ratio_goal
  use ovar,      only : ignition_source_magnitude
  use ovar,      only : ignition_max_temperature
  use ovar,      only : ignition_spark_cap
  use ovar,      only : ignition_source_radius
  use ovar,      only : ignition_source_location
  use ovar,      only : load_restart_file, restart_interval, restart_file
  use ovar,      only : restart_output_uses_host_roots
  use ovar,      only : iter_out_interval, results_interval
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
  use ovar,      only : Limiter_Option, Resolution_Marker_Variable
  use ovar,      only : Limiter_Marker_Min, Limiter_Marker_Max
  use ovar,      only : Limiter_Marker_dpsi
  use ovar,      only : Limiter_Marker_c1, Limiter_Marker_c2
  use ovar,      only : Filtering_Interval, Limiting_Interval
  use ovar,      only : loc_output_pts, loc_init_pts, loc_error_pts
  use ovar,      only : vortex_bigR, vortex_bigGam
  use ovar,      only : output_plot3d_to_tecplot, output_bface_array
  use ovar,      only : convergence_order_abs_res, convergence_order_max_res
  use ovar,      only : dump_memory_usage,use_bc_conflicts
  use ovar,      only : dump_fluxes,check_flux_point_coordinates
  use ovar,      only : use_unrolled_dot_products
  use ovar,      only : VortexGrid_DOF,VortexGrid_random_perturb
  use ovar,      only : VortexGrid_perturb_percent,VortexGrid_Lref
  use ovar,      only : channel_body_force,channel_init_file
  use ovar,      only : channel_spatial_ave_restart
  use ovar,      only : convert_restart_to_cgns, completely_disable_cgns
  use ovar,      only : lustre_stripe_count, lustre_stripe_size
  use ovar,      only : bc_input!,bc_input_t
  use ovar,      only : output_time_averaging, time_ave_file
  use ovar,      only : start_averaging_timestep, start_averaging_tstar
  use ovar,      only : time_ave_vel_is_axisymm
  use ovar,      only : time_scaling_factor
  use ovar,      only : time_average_restart_files, restart_list_file
  use ovar,      only : use_async_for_error
  use ovar,      only : ideal_total_temperature,limit_total_temperature
  use ovar,      only : output_limiter_statistics,output_filter_diagonals
  use ovar,      only : output_vorticity
  use ovar,      only : output_conserved_variables
  use ovar,      only : output_temperature
  use ovar,      only : output_thermo_variables
  use ovar,      only : output_gradients
  use ovar,      only : output_transport_variables
  use ovar,      only : output_stagnation_variables
  use ovar,      only : output_marker_variables
  use ovar,      only : output_species_mass_fractions
  use ovar,      only : output_extra_species_variables
  use ovar,      only : output_debug_species_variables
  use ovar,      only : output_grdnt_species_variables
  use ovar,      only : sem_profile_file
  use ovar,      only : force_nd_subroutines, newton_rlx
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
  integer, save :: projection_order = 3
  integer, save :: quadrature_order = 3
  integer, save :: output_order = -3
  integer, save :: error_order = -3
  integer, save :: interpto_solution_order = 3
  integer, save :: dmat_grad_temp_order = -1
  !
  ! ########################
  ! ###  Input Namelist  ###
  ! ########################
  !
  ! GRID INFORMATION
  namelist / input / grid_format, gridfile, plot3d_cutfile, plot3d_bcsfile, &
                     grid_scaling_factor
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
  namelist / input / interpto_solution_order
  ! BASE FR OPTIONS
  namelist / input / solution_order, flux_divergence, correction_function
  ! LOCATION OF SOLUTION, FLUX, AND TRIANGULAR POINTS
  namelist / input / loc_solution_pts, loc_flux_pts, loc_triangle_pts
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
                     restart_output_uses_host_roots
  ! OPTIONS TO USE SOLUTION FILE FOR INITIAL CONDITIONS
  namelist / input / solution_format, solution_file
  ! PROFILE FILE FOR SEM TURBULENCE INFLOW GENERATOR
  namelist / input / sem_profile_file
  ! STATISTICAL OUTPUT OPTIONS
  namelist / input / iter_out_interval, results_interval
  ! SOLUTION OUTPUT OPTIONS
  namelist / input / output_dir, output_interval, continuous_output, &
                     output_order, loc_output_pts, output_frequency
  ! SELECTION OF OUTPUT VARIABLES IN THE CGNS SOLUTION FILES
  namelist / input / output_vorticity, output_conserved_variables, &
                     output_temperature, output_thermo_variables, &
                     output_gradients, output_transport_variables, &
                     output_stagnation_variables, output_marker_variables, &
                     output_species_mass_fractions, &
                     output_extra_species_variables, &
                     output_debug_species_variables, &
                     output_grdnt_species_variables
  ! GOVERNING EQUATIONS
  namelist / input / governing_equations, turbulence_model
  ! CHEMICAL REACTION SOURCE TERMS
  namelist / input / disable_chemical_reactions, semi_implicit_chemistry, &
                     traceback_unconverged_newton_temps
                    !traceback_unconverged_newton_temps, &
                    !reactions_start_iter
  !
  ! IGNITION CONTROL SOURCE TERMS
  namelist / input / ignition_start_timestep, &
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
                     VortexGrid_perturb_percent
  ! CHANNEL FLOW PARAMETERS
  namelist / input / channel_body_force, channel_init_file, &
                     channel_spatial_ave_restart
  ! OPTIONS FOR MODIFYING BC ALGORITHMS
  namelist / input / walls_are_exact, sub_inflow_method, use_bc_conflicts
  ! CGNS OUTPUT OPTIONS
  namelist / input / cgns_use_queue, &
                     convert_restart_to_cgns, completely_disable_cgns
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
  namelist / input / force_nd_subroutines, dmat_grad_temp_order, newton_rlx, &
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
subroutine allocate_namelist_variables
  !
  !.. Use Statements ..
  use module_kind_types, only : zero,alloc_error,error_message
  !
  !.. Local Scalars ..
  integer :: n,ierr
  character(len=20) :: blanks
  character(len=50) :: array_name
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "allocate_namelist_variables"
  integer, parameter :: default_len = 150
  integer, parameter :: max_ncs = 7
  !
continue
  !
  ! Allocate bc_input before reading the namelist
  !
  allocate ( bc_input(1:20) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bc_input",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate the mass_frac and species_name components
  ! of bc_input before reading the namelist.
  !
#ifdef ALLOCATABLE_MASS_FRAC
  do n = 1,size(bc_input)
    allocate ( bc_input(n)%mass_frac(1:max_ncs) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,2) n,"mass_frac"
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    allocate ( bc_input(n)%species_name(1:max_ncs) , source="" , &
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
  allocate ( character(len=default_len) :: sem_profile_file , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sem_profile_file",1,__LINE__,__FILE__, &
                   ierr,error_message)
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
  write (sem_profile_file,blanks)
  write (output_dir,blanks)
  write (channel_init_file,blanks)
  write (lustre_stripe_size,blanks)
  write (time_ave_file,blanks)
  write (restart_list_file,blanks)
  !
  restart_file(1:7) = "out.rst"
  output_dir(1:1) = "."
  lustre_stripe_size(1:2)  = "4m"
  time_ave_file(1:7) = "out.ave"
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
  sem_profile_file = trim(adjustl(sem_profile_file))
  output_dir = trim(adjustl(output_dir))
  channel_init_file = trim(adjustl(channel_init_file))
  lustre_stripe_size = trim(adjustl(lustre_stripe_size))
  time_ave_file = trim(adjustl(time_ave_file))
  restart_list_file = trim(adjustl(restart_list_file))
  !
end subroutine resize_namelist_allocatables
!
!###############################################################################
!
end module input_namelist_mod
!
!###############################################################################
!
