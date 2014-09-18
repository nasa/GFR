module ovar
  !
  use module_kind_types, only : wp,fals,true
  use geovar, only : p_order,n_order,o_order,all_orders,q_order
  !
  implicit none
  !
  private
  !
  ! Public Input Namelist
  !
  public :: input
  !
  !
  ! ###########################
  ! ###  Module Parameters  ###
  ! ###########################
  !
  integer,  public, parameter :: Classic_RK = 1
  integer,  public, parameter :: TVD2_RK    = 2
  integer,  public, parameter :: TVD3_RK    = 3
  integer,  public, parameter :: TVD4_RK    = 4
  integer,  public, parameter :: CK4_RK     = 5
  !
  ! Using Intel 2013 compilers
  real(wp), public, parameter :: one_work_unit = 6.9205_wp
  ! Westmere nodes
 !real(wp), public, parameter :: one_work_unit = 25.551_wp/3.0_wp
  ! Sandy Bridge nodes
 !real(wp), public, parameter :: one_work_unit = 22.808_wp/3.0_wp
  ! Ivy Bridge nodes
 !real(wp), public, parameter :: one_work_unit = 19.096_wp/3.0_wp
  ! Using Intel version 10 compilers
 !real(wp), public, parameter :: one_work_unit = 7.026_wp
  !
  !
  !
  ! ################################
  ! ###  Module Local Variables  ###
  ! ################################
  !
  ! NOTE: Input defaults are set to run P3 80x80 Shu Vortex on [-10,10] domain
  !       for 1 period using Ga-LP, SSP-RK3, and Roe flux with entropy fix.
  !
  !.. Input Variables ..
  integer,  public, save :: all_points       = -1
  integer,  public, save :: iop_solpts       = 1
  integer,  public, save :: iop_flxpts       = 1
  integer,  public, save :: iop_tripts       = 2
  integer,  public, save :: iop_outpts       = 2
  integer,  public, save :: flux_divergence  = 0
  integer,  public, save :: iop_correction   = 3
  !
  integer,  public, save :: iRKtype     = 3
  integer,  public, save :: num_stages  = 3
  integer,  public, save :: invs_flux_method = 1
  integer,  public, save :: visc_flux_method = 1
  logical,  public, save :: walls_are_exact = true
  integer,  public, save :: visc_variables_to_ave = 1
  !
  integer,  public, save :: num_ts     = 25000
  real(wp), public, save :: Final_Time = -1.0_wp
  real(wp), public, save :: cfl        = 1.0_wp
  real(wp), public, save :: dt_in      = 0.0008_wp
  integer,  public, save :: Timestep_Type = 0
  !
  integer,  public, save :: iop_qdtpts = 1
  integer,  public, save :: dflux_opt  = 2
  !
  integer,  public, save :: metis_option = 2
  !
  integer,  public, save :: itestcase = 5
  !
  integer,  public, save :: irestart         = 0
  integer,  public, save :: restart_interval = 0
  !
  integer,  public, save :: iter_out_interval = 1
  integer,  public, save :: output_solution   = 1
  integer,  public, save :: output_interval   = 0
  integer,  public, save :: results_interval  = 0
  !
  real(wp), public, save :: rhoref  =     -1.00_wp ! should evaluate to 1.160833478437518_wp
  real(wp), public, save :: tref    =    300.00_wp
  real(wp), public, save :: pref    = 100000.00_wp
  real(wp), public, save :: machref =      0.40_wp
  real(wp), public, save :: aoaref  =     45.00_wp
  real(wp), public, save :: gam     =      1.40_wp
  real(wp), public, save :: rgasref =    287.15_wp
  real(wp), public, save :: muref   ! default evaluates to 1.84671531762365e-5_wp
  real(wp), public, save :: reyref  =     -1.00_wp ! default evaluates to 8.7319053606e+6
  !
  real(wp), public, save :: pref_nd
  real(wp), public, save :: vxref_nd
  real(wp), public, save :: vyref_nd
  real(wp), public, save :: rgasref_nd
  real(wp), public, save :: muref_nd
  real(wp), public, save :: cpref_nd
  !
  real(wp), public, save :: Pr = 0.72_wp
  real(wp), public, save :: Prt = 0.72_wp
  !
  real(wp), public, save :: suth_Sref = 110.4_wp
  real(wp), public, save :: suth_Tref = 273.0_wp
  real(wp), public, save :: suth_muref = 1.716e-5_wp
  !
  integer,  public, save :: Filter_Option = 0
  integer,  public, save :: Filter_Order  = 16
  real(wp), public, save :: Filter_Alpha  = 36.0_wp
  real(wp), public, save :: Filter_etac   = 0.0_wp
  !
  integer,  public, save :: init_loc  = 0
  integer,  public, save :: error_loc = 0
  !
  logical,  public, save :: bound_error       = .false.
  real(wp), public, save :: error_bound_limit = 2.0_wp
  !
  real(wp), public, save :: vortex_bigR = 1.0_wp
  real(wp), public, save :: vortex_bigGam = 5.0_wp
  !
  !.. Non-Input Variables ..
  real(wp), public, save :: aref
  real(wp), public, save :: uref
  real(wp), public, save :: vxref
  real(wp), public, save :: vyref
  real(wp), public, save :: cpref
  !
  real(wp), public, save :: dt_global
  !
  integer,  public, save :: itcur
  real(wp), public, save :: d_t
  real(wp), public, save :: time
  !
  integer,  public, save :: itrst
  real(wp), public, save :: tmrst
  real(wp), public, save :: dtrst
  !
  real(wp), public, save, allocatable, dimension(:) :: cv_in
  real(wp), public, save, allocatable, dimension(:) :: pv_in
  !
  real(wp), public, save :: iso_wall_temp = 300.0_wp
  !
  real(wp), public, save, dimension(4) :: bc_in
  real(wp), public, save, dimension(4) :: bc_top
  real(wp), public, save, dimension(4) :: bc_left
  real(wp), public, save, dimension(4) :: bc_right
  real(wp), public, save, dimension(3) :: bc_invar
  !
  real(wp), public, save :: p_exit
  real(wp), public, save :: local_volume
  real(wp), public, save :: total_volume
  !
  real(wp), public, save :: betark(4,4)
  !
  integer,  public, save :: nloc
  integer,  public, save :: kloc
  real(wp), public, save :: rlinf
  real(wp), public, save :: rlone
  real(wp), public, save :: rltwo
  !
  integer,  public, save, allocatable, dimension(:) :: unloc
  integer,  public, save, allocatable, dimension(:) :: ukloc
  real(wp), public, save, allocatable, dimension(:) :: ulinf
  real(wp), public, save, allocatable, dimension(:) :: ulone
  real(wp), public, save, allocatable, dimension(:) :: ultwo
  !
  integer,  public, save, allocatable, dimension(:) :: vnloc
  integer,  public, save, allocatable, dimension(:) :: vkloc
  real(wp), public, save, allocatable, dimension(:) :: vlinf
  real(wp), public, save, allocatable, dimension(:) :: vlone
  real(wp), public, save, allocatable, dimension(:) :: vltwo
  !
  real(wp), public, save, dimension(1:3,1:3,1:5) :: advec_error
  real(wp), public, save, dimension(1:5,1:2,1:5) :: vortex_error
  !
  integer,  public, save :: nentmax
  integer,  public, save :: kentmax
  real(wp), public, save :: xentmax
  real(wp), public, save :: yentmax
  real(wp), public, save :: entmax
  real(wp), public, save :: entL2
  !
  integer,  public, save :: nst
  !
  integer,  public, save :: Filter_Cutoff
  !
  integer,  public, save :: Total_Periods
  !
  logical,  public, save :: walltime_expiring
  real(wp), public, save :: max_wall_time
  real(wp), public, save :: wall_time_left
  real(wp), public, save :: beg_wall_time
  real(wp), public, save :: end_wall_time
  real(wp), public, save :: wall_time
  real(wp), public, save :: work_units
  !
  real(wp), public, save :: Period_Time
  integer,  public, save :: Period_Timesteps
  real(wp), public, save :: Grid_Max_X
  real(wp), public, save :: Grid_Max_Y
  real(wp), public, save :: Grid_Min_X
  real(wp), public, save :: Grid_Min_Y
  real(wp), public, save :: Grid_Length_X
  real(wp), public, save :: Grid_Length_Y
  !
  real(wp), public, save :: Density_Minimum_Residual
  real(wp), public, save :: Density_Maximum_Residual
  logical,  public, save :: vortex_has_collapsed
  integer,  public, save :: final_countdown
  !
  real(wp), public, save :: suth_c1
  !
  integer, public, save :: output_plot3d_to_tecplot = 0
  logical, public, save :: output_bface_array = fals
  !
  integer,            public, save :: grid_format = 2
  character(len=150), public, save :: gridfile = "grid_files/meilin/quad.40x40.2X.neu"
  character(len=150), public, save :: output_dir = "."
  character(len=150), public, save :: plot3d_cutfile = ""
  character(len=150), public, save :: plot3d_bcsfile = ""
  !
  !
  ! ########################
  ! ###  Input Namelist  ###
  ! ########################
  !
  ! GRID INFORMATION
  namelist / input / grid_format, gridfile, plot3d_cutfile, plot3d_bcsfile
  ! BASE FR OPTIONS
  namelist / input / n_order, iop_solpts, iop_flxpts, iop_tripts, iop_outpts, &
                     flux_divergence, iop_correction
  ! SOLVER OPTIONS
  namelist / input / iRKtype, num_stages, invs_flux_method, visc_flux_method
  ! VISCOUS OPTIONS
  namelist / input / visc_variables_to_ave
  ! TIME STEPPING INFORMATION
  namelist / input / num_ts, Final_Time, cfl, dt_in, Timestep_Type
  ! OVER-INTEGRATION OPTIONS
  namelist / input / p_order, q_order, iop_qdtpts, dflux_opt
  ! FORCE GLOBAL ORDER AND SOLUTION POINTS
  namelist / input / all_orders, all_points
  ! PARALLEL OPTIONS
  namelist / input / metis_option
  ! INITIALIZATION / TEST CASE
  namelist / input / itestcase
  ! RESTART OPTIONS
  namelist / input / irestart, restart_interval
  ! OUTPUT OPTIONS
  namelist / input / output_dir, iter_out_interval, output_solution, &
                     output_interval, results_interval, o_order
  ! NONDIMENSIONAL REFERENCE CONDITIONS
  namelist / input / machref, reyref, gam
  ! REFERENCE CONDITIONS
  namelist / input / rhoref, tref, pref, aoaref, rgasref
  ! LAMINAR AND TURBULENT PRANDTL NUMBERS
  namelist / input / Pr, Prt
  ! SUTHERLANDS LAW COEFFICIENTS
  namelist / input / suth_muref, suth_Tref, suth_Sref
  ! FILTER OPTIONS
  namelist / input / Filter_Option, Filter_etac, Filter_Order, Filter_Alpha
  ! BOUNDING ERROR OPTIONS
  namelist / input / bound_error, error_bound_limit
  ! INITIALIZATION/ERROR LOCATIONS
  namelist / input / init_loc, error_loc
  ! VORTEX CONSTANTS
  namelist / input / vortex_bigR, vortex_bigGam
  ! OPTION TO OUTPUT PLOT3D GRID TO TECPLOT IN UNSTRUCTURED FORMAT
  namelist / input / output_plot3d_to_tecplot
  ! LOGICAL IDENTIFYING HOW WALL SOLUTIONS ARE TREATED
  namelist / input / walls_are_exact
  ! VARIOUS DEBUG TYPE variables
  namelist / input / output_bface_array
  !
end module ovar
!
!###############################################################################
!
