module ovar
  !
  use module_kind_types, only : wp,lk,ldk,fals,true
  use module_kind_types, only : Euler_Eqns
  use module_kind_types, only : Gambit_Neutral
  use module_kind_types, only : Legendre_Gauss,Legendre_Gauss_Lobatto
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
  !
  implicit none
  !
  private
  !
  public :: recon_cpu_time
  public :: ovar_memory_usage
  public :: adjust_cfl_cycles
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
  ! #############################
  ! ####   Input Variables   ####
  ! #############################
  !
  integer,  public, save :: all_points       = -1
  integer,  public, save :: iop_solpts       = Legendre_Gauss
  integer,  public, save :: iop_flxpts       = Legendre_Gauss
  integer,  public, save :: iop_tripts       = 2
  integer,  public, save :: flux_divergence  = 1
  integer,  public, save :: correction_function = Correction_gDG
  !
  integer,  public, save :: governing_equations = Euler_Eqns
  !
  integer,  public, save :: iRKtype     = TVD3_RK
  integer,  public, save :: num_stages  = 3
  integer,  public, save :: invs_flux_method = invs_flux_ROE
  integer,  public, save :: visc_flux_method = visc_flux_BR2
  logical(lk),  public, save :: walls_are_exact = true
  integer,  public, save :: visc_variables_to_ave = visc_interface_fluxes
  !
  integer,  public, save :: num_ts     = 25000
  real(wp), public, save :: Final_Time = -1.0_wp
  real(wp), public, save :: cfl_beg    = 0.1_wp
  real(wp), public, save :: cfl_end    = 1.0_wp
  integer,  public, save :: cfl_cycles = 0
  real(wp), public, save :: cfl        = 1.0_wp
  real(wp), public, save :: dt_in      = 0.0008_wp
  integer,  public, save :: Timestep_Type = 0
  integer,  public, save :: minimum_timesteps = 100
  !
  integer,  public, save :: iop_qdtpts = 1
  integer,  public, save :: dflux_opt  = 2
  !
  integer,  public, save :: metis_option = 3
  !
  integer,  public, save :: mms_opt = 0
  real(wp), public, save :: mms_init = 0.1_wp
  logical(lk),  public, save :: mms_output_solution = fals
  !
  integer,  public, save :: itestcase = 5
  !
  integer,  public, save :: irestart         = 0
  integer,  public, save :: restart_interval = 0
  character(len=150), public, save :: restart_file = "out.rst"
  !
  integer,  public, save :: iter_out_interval =  1
  integer,  public, save :: results_interval  =  0
  !
  character(len=150), public, save :: output_dir = "."
  integer,            public, save :: output_interval = -5
  logical(lk),            public, save :: correct_output  = fals
  logical(lk),            public, save :: continuous_output  = fals
  logical(lk),            public, save :: postproc_debug_grid = fals
  !
  real(wp), public, save :: grid_scaling_factor = 1.00_wp
  !
  real(wp), public, save :: rhoref  =     -1.00_wp
  ! rhoref default evaluates to 1.160833478437518_wp
  real(wp), public, save :: tref    =    300.00_wp
  real(wp), public, save :: pref    = 100000.00_wp
  real(wp), public, save :: machref =      0.40_wp
  real(wp), public, save :: gam     =      1.40_wp
  real(wp), public, save :: rgasref =    287.15_wp
  real(wp), public, save :: reyref  =     -1.00_wp
  ! reyref default evaluates to 8.7319053606e+6
  !
  real(wp), public, save :: ptot2p_ratio = -1.00_wp
  real(wp), public, save :: ptotref = -1.00_wp
  real(wp), public, save :: ttotref = -1.00_wp
  !
  real(wp), public, save :: alpha_aoaref = 0.00_wp
  real(wp), public, save :: beta_aoaref  = 0.00_wp
  !
  real(wp), public, save :: Pr = 0.72_wp
  real(wp), public, save :: Prt = 0.72_wp
  !
  real(wp), public, save :: suth_Sref = 110.4_wp
  real(wp), public, save :: suth_Tref = 273.0_wp
  real(wp), public, save :: suth_muref = 1.716e-5_wp
  !
  integer,  public, save :: turbulence_model = laminar_flow
  real(wp), public, save :: tinref = 0.0_wp
  real(wp), public, save :: tvrref = 1.0e-6_wp
  !
  integer,  public, save :: Filter_Option = 0
  integer,  public, save :: Filter_Order  = 16
  real(wp), public, save :: Filter_Alpha  = 36.0_wp
  real(wp), public, save :: Filter_etac   = 0.0_wp
  real(wp), public, save :: Filter_eps    = 0.0_wp
  !
  ! Default to output solution at Legendre Gauss Lobatto points
  integer,  public, save :: loc_output_pts = Legendre_Gauss_Lobatto
  ! Default to initialize solution at solution points
  integer,  public, save :: loc_init_pts   = 0
  ! Default to compute solution error at Legendre Gauss points
  integer,  public, save :: loc_error_pts  = Legendre_Gauss
  !
  real(wp), public, save :: vortex_bigR = 1.0_wp
  real(wp), public, save :: vortex_bigGam = 5.0_wp
  !
  integer,  public, save :: VortexGrid_DOF = 64
  logical(lk),  public, save :: VortexGrid_random_perturb = fals
  integer,  public, save :: VortexGrid_perturb_percent = 30
  real(wp), public, save :: VortexGrid_Lref = 10.0_wp
  !
  real(wp),           public, save :: channel_body_force = 0.0_wp
  character(len=150), public, save :: channel_init_file = ""
  !
  logical(lk), public, save :: output_plot3d_to_tecplot = fals
  logical(lk), public, save :: output_bface_array = fals
  !
  integer,            public, save :: grid_format = Gambit_Neutral
  character(len=150), public, save :: gridfile = &
                                      "grid_files/meilin/quad.40x40.2X.neu"
  character(len=150), public, save :: plot3d_cutfile = ""
  character(len=150), public, save :: plot3d_bcsfile = ""
  !
  integer, public, save :: plot3d_all_strides = -1
  integer, public, save :: plot3d_i_stride    = -1
  integer, public, save :: plot3d_j_stride    = -1
  integer, public, save :: plot3d_k_stride    = -1
  !
  logical(lk), public, save :: interpolate_before_output_variable = true
  logical(lk), public, save :: completely_disable_cgns = fals
  !
  real(wp), public, save :: convergence_order_abs_res = real( -16 , kind=wp )
  real(wp), public, save :: convergence_order_max_res = real( -16 , kind=wp )
  !
  logical(lk), public, save :: use_bc_conflicts = true
  !
  logical(lk), public, save :: dump_memory_usage = fals
  logical(lk), public, save :: dump_fluxes       = fals
  !
  logical(lk), public, save :: check_flux_point_coordinates = fals
  !
  logical(lk), public, save :: profile_io_all     = fals
  logical(lk), public, save :: profile_io_cgns    = fals
  logical(lk), public, save :: profile_io_restart = fals
  !
  logical(lk), public, save :: output_bnd_profiles = fals
  !
  logical(lk), public, save :: use_unrolled_dot_products = true
  !
  logical(lk), public, save :: convert_restart_to_cgns = fals
  !
  logical(lk), public, save :: use_old_restart_format = fals
  !
  integer,           public, save :: lustre_stripe_count = 64
  character(len=10), public, save :: lustre_stripe_size  = "4m"
  !
  integer, public, save :: sub_inflow_method = 1
  !
  ! #################################
  ! ####   Non-Input Variables   ####
  ! #################################
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
  real(wp), public, save :: ptotref_nd
  real(wp), public, save :: ttotref_nd
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
  !
  real(wp), public, save :: cfl_cycle_start
  real(wp), public, save :: cfl_cycle_end
  !
  integer,  public, save :: itrst = 0
  real(wp), public, save :: tmrst = zero
  real(wp), public, save :: dtrst = zero
  !
  real(wp), public, save, allocatable, dimension(:) :: cv_in
  real(wp), public, save, allocatable, dimension(:) :: pv_in
  !
  real(wp), public, save :: iso_wall_temp = 300.0_wp
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
  integer,  public, save :: nst
  !
  integer,  public, save :: Total_Periods
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
  type, public :: relaxation_t
    real(wp) :: value = 1.0_wp
    real(wp) :: time  = 0.0_wp
    integer  :: iter  = 0
  end type relaxation_t
  !
  type, public :: bc_input_t
    !
    character(len=32) :: name = ""
    !
    integer :: set_bc_value = bc_unknown
    !
    real(wp) :: t_static   = -huge(0.0_wp)
    real(wp) :: p_static   = -huge(0.0_wp)
    real(wp) :: rho_static = -huge(0.0_wp)
    real(wp) :: mach       = -huge(0.0_wp)
    real(wp) :: vx         = 0.0_wp
    real(wp) :: vy         = 0.0_wp
    real(wp) :: vz         = 0.0_wp
    real(wp) :: t_total    = -huge(0.0_wp)
    real(wp) :: p_total    = -huge(0.0_wp)
    real(wp) :: rho_total  = -huge(0.0_wp)
    real(wp) :: alpha_aoa  = 0.0_wp
    real(wp) :: beta_aoa   = 0.0_wp
    real(wp) :: wall_temp  = -huge(0.0_wp)
    !
    type(relaxation_t) :: relax
    !
  contains
    !
    procedure :: has_default_flow_conditions
    !
  end type bc_input_t
  !
  type(bc_input_t), public, save :: bc_input(1:20)
  !
  type, public :: bc_in_t
    real(wp), allocatable :: cv(:)
    real(wp), allocatable :: pv(:)
    real(wp) :: wall_temperature
    type(relaxation_t) :: relax
  end type bc_in_t
  !
  type(bc_in_t), public, save, allocatable :: bc_in(:)
  !
  character(len=32), public, save, allocatable :: bc_names(:)
  !
contains
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
  if (abs(this%t_static  +huge(zero)) > eps6) return
  if (abs(this%p_static  +huge(zero)) > eps6) return
  if (abs(this%rho_static+huge(zero)) > eps6) return
  if (abs(this%mach      +huge(zero)) > eps6) return
  if (abs(this%vx                   ) > eps6) return
  if (abs(this%vy                   ) > eps6) return
  if (abs(this%vz                   ) > eps6) return
  if (abs(this%t_total   +huge(zero)) > eps6) return
  if (abs(this%p_total   +huge(zero)) > eps6) return
  if (abs(this%rho_total +huge(zero)) > eps6) return
  if (abs(this%alpha_aoa            ) > eps6) return
  if (abs(this%beta_aoa             ) > eps6) return
  if (abs(this%wall_temp +huge(zero)) > eps6) return
  !
  if (abs(this%relax%value- one) > eps6) return
  if (abs(this%relax%time -zero) > eps6) return
  if (    this%relax%iter        /=    0) return
  !
  return_value = true
  !
end function has_default_flow_conditions
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
subroutine recon_cpu_time(time_flag)
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
    beg_wall_time = mpi_wtime()
    work_units = zero
    !
  else if (time_flag == 1) then
    !
    end_wall_time = mpi_wtime()
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
end subroutine recon_cpu_time
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
continue
  !
  iou = iout
  if (present(iunit)) iou = iunit
  !
  call total%reset
  !
  write (iou,1)
  !
  if (allocated(cv_in)) then
    call array%set( storage_size(cv_in,kind=inttype) , &
                            size(cv_in,kind=inttype) )
    call total%add(array)
    write (array_name,5) "cv_in"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "cv_in"
  end if
  if (allocated(pv_in)) then
    call array%set( storage_size(pv_in,kind=inttype) , &
                            size(pv_in,kind=inttype) )
    call total%add(array)
    write (array_name,5) "pv_in"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "pv_in"
  end if
  if (allocated(unloc)) then
    call array%set( storage_size(unloc,kind=inttype) , &
                            size(unloc,kind=inttype) )
    call total%add(array)
    write (array_name,5) "unloc"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "unloc"
  end if
  if (allocated(ukloc)) then
    call array%set( storage_size(ukloc,kind=inttype) , &
                            size(ukloc,kind=inttype) )
    call total%add(array)
    write (array_name,5) "ukloc"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "ukloc"
  end if
  if (allocated(ulinf)) then
    call array%set( storage_size(ulinf,kind=inttype) , &
                            size(ulinf,kind=inttype) )
    call total%add(array)
    write (array_name,5) "ulinf"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "ulinf"
  end if
  if (allocated(ulone)) then
    call array%set( storage_size(ulone,kind=inttype) , &
                            size(ulone,kind=inttype) )
    call total%add(array)
    write (array_name,5) "ulone"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "ulone"
  end if
  if (allocated(ultwo)) then
    call array%set( storage_size(ultwo,kind=inttype) , &
                            size(ultwo,kind=inttype) )
    call total%add(array)
    write (array_name,5) "ultwo"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "ultwo"
  end if
  if (allocated(vnloc)) then
    call array%set( storage_size(vnloc,kind=inttype) , &
                            size(vnloc,kind=inttype) )
    call total%add(array)
    write (array_name,5) "vnloc"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "vnloc"
  end if
  if (allocated(vkloc)) then
    call array%set( storage_size(vkloc,kind=inttype) , &
                            size(vkloc,kind=inttype) )
    call total%add(array)
    write (array_name,5) "vkloc"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "vkloc"
  end if
  if (allocated(vlinf)) then
    call array%set( storage_size(vlinf,kind=inttype) , &
                            size(vlinf,kind=inttype) )
    call total%add(array)
    write (array_name,5) "vlinf"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "vlinf"
  end if
  if (allocated(vlone)) then
    call array%set( storage_size(vlone,kind=inttype) , &
                            size(vlone,kind=inttype) )
    call total%add(array)
    write (array_name,5) "vlone"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "vlone"
  end if
  if (allocated(vltwo)) then
    call array%set( storage_size(vltwo,kind=inttype) , &
                            size(vltwo,kind=inttype) )
    call total%add(array)
    write (array_name,5) "vltwo"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "vltwo"
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
  use ovar,      only : plot3d_cutfile, plot3d_bcsfile, grid_scaling_factor
  use ovar,      only : iop_solpts, iop_flxpts, iop_tripts, iop_qdtpts
  use ovar,      only : flux_divergence, correction_function, dflux_opt
  use ovar,      only : iRKtype, num_stages
  use ovar,      only : invs_flux_method, visc_flux_method
  use ovar,      only : visc_variables_to_ave, walls_are_exact
  use ovar,      only : num_ts, Final_Time, dt_in
  use ovar,      only : cfl, cfl_beg, cfl_end, cfl_cycles
  use ovar,      only : Timestep_Type, minimum_timesteps
  use ovar,      only : all_points, metis_option, use_old_restart_format
  use ovar,      only : mms_opt, mms_init, mms_output_solution
  use ovar,      only : itestcase, governing_equations
  use ovar,      only : irestart, restart_interval, restart_file
  use ovar,      only : iter_out_interval, results_interval
  use ovar,      only : output_dir, output_interval
  use ovar,      only : correct_output, continuous_output
  use ovar,      only : postproc_debug_grid
  use ovar,      only : output_bnd_profiles,sub_inflow_method
  use ovar,      only : machref, reyref, gam, Pr, Prt
  use ovar,      only : ptotref, ttotref, ptot2p_ratio
  use ovar,      only : rhoref, tref, pref, rgasref
  use ovar,      only : turbulence_model, tinref, tvrref
  use ovar,      only : alpha_aoaref, beta_aoaref
  use ovar,      only : profile_io_all,profile_io_cgns,profile_io_restart
  use ovar,      only : suth_muref, suth_Tref, suth_Sref
  use ovar,      only : Filter_Option, Filter_etac, Filter_Order
  use ovar,      only : Filter_Alpha, Filter_eps
  use ovar,      only : loc_output_pts, loc_init_pts, loc_error_pts
  use ovar,      only : vortex_bigR, vortex_bigGam
  use ovar,      only : output_plot3d_to_tecplot, output_bface_array
  use ovar,      only : interpolate_before_output_variable
  use ovar,      only : convergence_order_abs_res, convergence_order_max_res
  use ovar,      only : dump_memory_usage,use_bc_conflicts
  use ovar,      only : dump_fluxes,check_flux_point_coordinates
  use ovar,      only : use_unrolled_dot_products
  use ovar,      only : VortexGrid_DOF,VortexGrid_random_perturb
  use ovar,      only : VortexGrid_perturb_percent,VortexGrid_Lref
  use ovar,      only : channel_body_force,channel_init_file
  use ovar,      only : convert_restart_to_cgns, completely_disable_cgns
  use ovar,      only : lustre_stripe_count, lustre_stripe_size
  use ovar,      only : bc_input!,bc_input_t
  use order_mod, only : p_order,n_order,q_order
  use order_mod, only : o_order,e_order,all_orders
  use module_cgns_types, only : cgns_use_queue
  !
  implicit none
  !
  ! Set everything in the module to public
  !
  public
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
                     plot3d_all_strides
  ! BASE FR OPTIONS
  namelist / input / n_order, flux_divergence, correction_function
  ! LOCATION OF SOLUTION, FLUX, AND TRIANGULAR POINTS
  namelist / input / iop_solpts, iop_flxpts, iop_tripts
  ! LOCATION OF OUTPUT, INITIALIZATION, AND ERROR POINTS
  namelist / input / loc_output_pts, loc_init_pts
  ! SOLVER OPTIONS
  namelist / input / iRKtype, num_stages, invs_flux_method, visc_flux_method
  ! VISCOUS OPTIONS
  namelist / input / visc_variables_to_ave
  ! TIME STEPPING INFORMATION
  namelist / input / num_ts, Final_Time, dt_in, Timestep_Type, minimum_timesteps
  ! CFL CONTROL VARIABLES
  namelist / input / cfl, cfl_beg, cfl_end, cfl_cycles
  ! OVER-INTEGRATION OPTIONS
  namelist / input / p_order, q_order, iop_qdtpts, dflux_opt
  ! ERROR-INTEGRATION OPTIONS
  namelist / input / e_order, loc_error_pts
  ! FORCE GLOBAL ORDER AND SOLUTION POINTS
  namelist / input / all_orders, all_points
  ! PARALLEL OPTIONS
  namelist / input / metis_option
  ! METHOD OF MANUFACTURED SOLUTIONS (MMS) OPTIONS
  namelist / input / mms_opt, mms_init, mms_output_solution
  ! INITIALIZATION / TEST CASE
  namelist / input / itestcase
  ! RESTART OPTIONS
  namelist / input / irestart, restart_interval, use_old_restart_format, &
                     restart_file
  ! STATISTICAL OUTPUT OPTIONS
  namelist / input / iter_out_interval, results_interval
  ! SOLUTION OUTPUT OPTIONS
  namelist / input / output_dir, output_interval, o_order, &
                     correct_output, continuous_output
  ! ACTIVATE SPECIAL GRID FOR DEBUGGING THE POSTPROCESSOR
  namelist / input / postproc_debug_grid
  ! GOVERNING EQUATIONS
  namelist / input / governing_equations, turbulence_model
  ! NONDIMENSIONAL REFERENCE CONDITIONS
  namelist / input / machref, reyref, gam
  ! TOTAL REFERENCE CONDITIONS
  namelist / input / ptotref, ttotref, ptot2p_ratio
  ! REFERENCE CONDITIONS
  namelist / input / rhoref, tref, pref, rgasref
  ! REFERENCE TURBULENCE CONDITIONS
  namelist / input / tinref, tvrref
  ! REFERENCE ANGLES OF ATTACK
  namelist / input / alpha_aoaref, beta_aoaref
  ! LAMINAR AND TURBULENT PRANDTL NUMBERS
  namelist / input / Pr, Prt
  ! SUTHERLANDS LAW COEFFICIENTS
  namelist / input / suth_muref, suth_Tref, suth_Sref
  ! FILTER OPTIONS
  namelist / input / Filter_Option, Filter_etac, Filter_Order, &
                     Filter_Alpha, Filter_eps
  ! VORTEX CONSTANTS
  namelist / input / vortex_bigR, vortex_bigGam
  ! TAYLOR-GREEN VORTEX CONSTANTS
  namelist / input / VortexGrid_DOF, VortexGrid_random_perturb, &
                     VortexGrid_perturb_percent, VortexGrid_Lref
  ! CHANNEL FLOW PARAMETERS
  namelist / input / channel_body_force, channel_init_file
  ! OPTION TO OUTPUT PLOT3D GRID TO TECPLOT IN UNSTRUCTURED FORMAT
  namelist / input / output_plot3d_to_tecplot
  ! LOGICAL IDENTIFYING HOW WALL SOLUTIONS ARE TREATED
  namelist / input / walls_are_exact
  ! VARIOUS DEBUG TYPE VARIABLES
  namelist / input / output_bface_array
  ! CGNS OUTPUT OPTIONS
  namelist / input / cgns_use_queue, interpolate_before_output_variable, &
                     convert_restart_to_cgns, completely_disable_cgns
  ! CONVERGENCE CRITERIA
  namelist / input / convergence_order_abs_res, convergence_order_max_res
  ! MEMORY USAGE
  namelist / input / dump_memory_usage
  ! PROFILE PARALLEL IO
  namelist / input / profile_io_all, profile_io_cgns, profile_io_restart
  ! TEMPORARY DEBUGGING LOGICALS
  namelist / input / use_bc_conflicts,sub_inflow_method
  ! PROFILE OUTPUT OPTIONS
  namelist / input / output_bnd_profiles
  ! CONTROL DOT_PRODUCT PROCEDURES BOUND TO MATRIX DERIVED-TYPE
  namelist / input / use_unrolled_dot_products
  ! DUMP FLUX INFORMATION FOR DEBUGGING
  namelist / input / dump_fluxes, check_flux_point_coordinates
  ! LUSTRE FILE SYSTEM STRIPING PARAMETERS
  namelist / input / lustre_stripe_count, lustre_stripe_size
  ! INPUT BOUNDARY CONDITIONS
  namelist / input / bc_input
  !
  ! CFL ADJUSTMENT NAMELIST
  !
  namelist / cfl_adjustment / cfl, cfl_beg, cfl_end, cfl_cycles
  !
  ! RESTART INTERUPTION NAMELIST
  !
  namelist / restart_interupt / restart_file
  !
end module input_namelist_mod
!
!###############################################################################
!
