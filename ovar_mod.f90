module ovar
  !
  use module_kind_types, only : wp,fals,true
  use geovar, only : p_order,n_order,o_order,all_orders,q_order
  !
  implicit none
  !
  private
  !
  public :: ovar_alloc_status
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
  real(wp), public, save :: gam     =      1.40_wp
  real(wp), public, save :: rgasref =    287.15_wp
  real(wp), public, save :: muref   ! default evaluates to 1.84671531762365e-5_wp
  real(wp), public, save :: reyref  =     -1.00_wp ! default evaluates to 8.7319053606e+6
  !
  real(wp), public, save :: alpha_aoaref = 45.00_wp
  real(wp), public, save :: beta_aoaref  =  0.00_wp
  !
  real(wp), public, save :: pref_nd
  real(wp), public, save :: rgasref_nd
  real(wp), public, save :: muref_nd
  real(wp), public, save :: cpref_nd
  real(wp), public, save, dimension(3) :: velref_nd
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
  real(wp), public, save :: cpref
  real(wp), public, save, dimension(3) :: velref
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
  real(wp), public, save :: rl2mx
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
  logical, public, save :: force_parallel_cgns = fals
  logical, public, save :: force_serial_cgns = fals
  logical, public, save :: cgns_use_queue = true
  logical, public, save :: cgns_is_parallel
  logical, public, save :: interpolate_before_output_variable = true
  !
  real(wp), public, save :: convergence_order_abs_res = real( -16 , kind=wp )
  real(wp), public, save :: convergence_order_max_res = real( -16 , kind=wp )
  logical,  public, save :: convergence_reached = fals
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
  namelist / input / rhoref, tref, pref, rgasref
  ! REFERENCE ANGLES OF ATTACK
  namelist / input / alpha_aoaref, beta_aoaref
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
  ! VARIOUS DEBUG TYPE VARIABLES
  namelist / input / output_bface_array
  ! CGNS OUTPUT OPTIONS
  namelist / input / force_parallel_cgns,force_serial_cgns,cgns_use_queue, &
                     interpolate_before_output_variable
  ! CONVERGENCE CRITERIA
  namelist / input / convergence_order_abs_res, convergence_order_max_res
  !
contains
!
!###############################################################################
!
subroutine ovar_alloc_status(iunit)
  !
  !.. Use Statements ..
  use module_memory
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iunit
  !
  !.. Local Scalars ..
  type(memory) :: array
  type(memory) :: total
  !
continue
  !
  total%siz = 0
  total%mem = 0
  !
  write (iunit,1)
  !
  if (allocated(cv_in)) then
    array%siz = size(cv_in)
    array%mem = array%siz*storage_size(cv_in)
    total%mem = total%mem + array%mem
    write (iunit,2) "'cv_in'",array%get_units()
  end if
  if (allocated(pv_in)) then
    array%siz = size(pv_in)
    array%mem = array%siz*storage_size(pv_in)
    total%mem = total%mem + array%mem
    write (iunit,2) "'pv_in'",array%get_units()
  end if
  if (allocated(unloc)) then
    array%siz = size(unloc)
    array%mem = array%siz*storage_size(unloc)
    total%mem = total%mem + array%mem
    write (iunit,2) "'unloc'",array%get_units()
  end if
  if (allocated(ukloc)) then
    array%siz = size(ukloc)
    array%mem = array%siz*storage_size(ukloc)
    total%mem = total%mem + array%mem
    write (iunit,2) "'ukloc'",array%get_units()
  end if
  if (allocated(ulinf)) then
    array%siz = size(ulinf)
    array%mem = array%siz*storage_size(ulinf)
    total%mem = total%mem + array%mem
    write (iunit,2) "'ulinf'",array%get_units()
  end if
  if (allocated(ulone)) then
    array%siz = size(ulone)
    array%mem = array%siz*storage_size(ulone)
    total%mem = total%mem + array%mem
    write (iunit,2) "'ulone'",array%get_units()
  end if
  if (allocated(ultwo)) then
    array%siz = size(ultwo)
    array%mem = array%siz*storage_size(ultwo)
    total%mem = total%mem + array%mem
    write (iunit,2) "'ultwo'",array%get_units()
  end if
  if (allocated(vnloc)) then
    array%siz = size(vnloc)
    array%mem = array%siz*storage_size(vnloc)
    total%mem = total%mem + array%mem
    write (iunit,2) "'vnloc'",array%get_units()
  end if
  if (allocated(vkloc)) then
    array%siz = size(vkloc)
    array%mem = array%siz*storage_size(vkloc)
    total%mem = total%mem + array%mem
    write (iunit,2) "'vkloc'",array%get_units()
  end if
  if (allocated(vlinf)) then
    array%siz = size(vlinf)
    array%mem = array%siz*storage_size(vlinf)
    total%mem = total%mem + array%mem
    write (iunit,2) "'vlinf'",array%get_units()
  end if
  if (allocated(vlone)) then
    array%siz = size(vlone)
    array%mem = array%siz*storage_size(vlone)
    total%mem = total%mem + array%mem
    write (iunit,2) "'vlone'",array%get_units()
  end if
  if (allocated(vltwo)) then
    array%siz = size(vltwo)
    array%mem = array%siz*storage_size(vltwo)
    total%mem = total%mem + array%mem
    write (iunit,2) "'vltwo'",array%get_units()
  end if
  !
  total = total%get_units()
  write (iunit,3) total%mem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module ovar")
  2 format (" Array: ",a7,dt)
  3 format (" Total memory = ",i0,1x,a)
  !
end subroutine ovar_alloc_status
!
!###############################################################################
!
end module ovar
!
!###############################################################################
!
module eqn_idx
  !
  use module_kind_types
  !
  implicit none
  !
  private
  !
  ! Public Procedures
  !
  public :: create_equation_indices
  !
  ! NEC : index for continuity equation
  !
  integer, public, protected, save :: nec
  !
  ! NMX : index for x-momentum equation
  ! NMY : index for y-momentum equation
  ! NMZ : index for z-momentum equation
  !
  integer, public, protected, save :: nmx
  integer, public, protected, save :: nmy
  integer, public, protected, save :: nmz
  !
  ! NMB : index for first momentum equation
  ! NME : index for last  momentum equation
  !
  integer, public, protected, save :: nmb
  integer, public, protected, save :: nme
  !
  ! NEC : index for energy equation
  !
  integer, public, protected, save :: nee
  !
  ! NQ : total number of equations
  !
  integer, public, protected, save :: nq = 4
  !
  ! NGOV : Equation control parameters 
  !
  !        NOTE: THIS IS COPIED FROM VULCAN MOST IS NOT USED
  !
  !        NGOV(1,1) = starting index of continuity equations
  !        NGOV(2,1) = ending   index of continuity equations
  !                    and the sum of all densities preceeding
  !        NGOV(3,1) = number of continuity equations
  !                  = 1 when thermally and calorically perfect
  !                  = no. of chemical species+1 when thermally perfect
  !                  = no. of chemical species+1 when thermally imperfect
  !        NGOV(4,1) = starting index of continuity equation source terms
  !        NGOV(5,1) = ending   index of continuity equation source terms
  !        NGOV(6,1) = no. of continuity equation source terms
  !                  = 0 for thermally and calorically perfect gas
  !        NGOV(7,1) = equation RHS solution procedure
  !                    0 = explicit
  !                    1 = point implicit
  !        NGOV(1,2) = starting index of momentum equations
  !        NGOV(2,2) = ending   index of momentum equations
  !        NGOV(3,2) = number of momentum equations
  !        NGOV(6,2) = no. of momentum equation source terms
  !                  = 0
  !        NGOV(7,2) = equation RHS solution procedure
  !                    0 = explicit
  !                    1 = point implicit
  !        NGOV(1,3) = starting index of energy equations
  !        NGOV(2,3) = ending   index of energy equations
  !        NGOV(3,3) = number of energy equations
  !                  = 1 when thermally and calorically perfect
  !                  = 1 when thermally perfect calorically imperfect
  !                  = 2 when thermally imperfect
  !        NGOV(6,3) = no. of energy equation source terms
  !                  = 0 when thermally and calorically perfect
  !                  = 0 when thermally and calorically imperfect
  !                  = N when thermally imperfect
  !        NGOV(1,4) = starting index of turbulence equations
  !        NGOV(2,4) = ending   index of turbulence equations
  !        NGOV(3,4) = number of turbulence equations
  !        NGOV(4,4) = starting index of turbulence equation source terms
  !        NGOV(5,4) = ending   index of turbulence equation source terms
  !        NGOV(6,4) = no. of turbulence equation source terms
  !        NGOV(7,4) = equation RHS solution procedure
  !                    0 = explicit
  !        NGOV(1,5) = starting index of "G" equations
  !        NGOV(2,5) = ending   index of "G" equations
  !        NGOV(3,5) = number of "G" equations
  !        NGOV(4,5) = starting index of "G" equation source terms
  !        NGOV(5,5) = ending   index of "G" equation source terms
  !        NGOV(6,5) = no. of "G" equation source terms
  !        NGOV(7,5) = equation RHS solution procedure
  !                    0 = explicit
  !                    1 = point implicit
  !
  integer, public, protected, save, dimension(3,4) :: ngov
  !
contains
!
!###############################################################################
!
subroutine create_equation_indices(num_dim)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: num_dim
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_equation_indices"
  !
continue
  !
  if (num_dim == 1) then
    !
    call stop_recon(stop_mpi,pname,"1D is not supported")
    !
  end if
  !
  ! Continuity equations
  !
  ngov(1:3,1) = 1
  !
  ! Momentum equations
  !
  ngov(3,2) = num_dim
  ngov(1,2) = ngov(2,1) + 1
  ngov(2,2) = ngov(1,2) + ngov(3,2) - 1
  !
  ! Energy equations
  !
  ngov(3,3) = 1
  ngov(1,3) = ngov(2,2) + 1
  ngov(2,3) = ngov(1,3) + ngov(3,3) - 1
  !
  ! Turbulence equations (NOT WORKING)
  !
  ngov(3,4) = 0
  ngov(1,4) = ngov(2,3) + 1
  ngov(2,4) = ngov(1,4) + ngov(3,4) - 1
  !
  ! Set the total number of equations
  !
  nq = ngov(3,1) + ngov(3,2) + ngov(3,3) + ngov(3,4)
  !
  ! Set the indices that will be used through the code
  !
  nec = ngov(2,1)
  nmx = ngov(1,2)
  nmy = ngov(1,2) + 1
  nmz = merge( ngov(1,2) + 2 , 0 , num_dim == 3 )
  nee = ngov(2,3)
  !
  nmb = nmx
  nme = nec + ngov(3,2)
  !
end subroutine create_equation_indices
!
!###############################################################################
!
end module eqn_idx
!
!###############################################################################
!
