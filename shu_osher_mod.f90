!###############################################################################
!> Routines for Shu-Osher shock/entropy wave interaction problem.
!> Ref:
!> * Shu, C-W and Osher, S.,
!>   "Efficient implementation of essentially non-oscillatory
!>    shock-capturing schemes, II," J Comput Phys 1989; 83(1):32–78.
!>    http://dx.doi.org/10.1016/0021-9991(89)90222-2
!###############################################################################
module shu_osher_mod
  !
  use module_kind_types
  !
  implicit none
  !
  private
  !
  public :: initialize_shu_osher_problem
  public :: shu_osher_init
  !
  ! Module constants needed for the initialization functions.
  ! Flow conditions for the supersonic inflow and left domain initialization.
  !
  real(wp), parameter :: default_rho  = 3.857143_wp    ! [kg/m3]
  real(wp), parameter :: default_vx   = 2.629369_wp    ! [m/s]
  real(wp), parameter :: default_vy   = zero           ! [m/s]
  real(wp), parameter :: default_vz   = zero           ! [m/s]
  real(wp), parameter :: default_pres = 10.333333_wp   ! [Pa]
contains
!
!###############################################################################
!> Initialize parameters for the Shu-Osher problem.
!>
!> Sets the correct freestream reference conditions.
!> Sets the correct inflow boundary condition inputs for bc_input.
!> Does not initialize the solution field.
!###############################################################################
subroutine initialize_shu_osher_problem()
  !
  !.. Use Statements ..
  use ovar, only : bc_input_t, bc_input
  use ovar, only : rgasref,gam
  use ovar, only : alpha_aoaref,beta_aoaref
  use ovar, only : rhoref,pref,tref,uref,aref,machref,reyref
  use thermo_mod, only : rgas_cpg

  !
  !.. Local Variables ..
  type(bc_input_t ) :: default_bc_input
  integer :: n, bc_idx
  real(wp) :: muref_scale,tref_scale,reyn
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_shu_osher_problem"
  !
continue
  !
  ! Find the index in bc_input for the supersonic inflow.
  ! If missing, find the first empty name.
  !
  bc_idx = 0
  find_sup_inflow_bc: &
  do n = 1,size(bc_input)
    if (trim(adjustl(uppercase(bc_input(n)%name))) == "SUP_INFLOW") then
      bc_idx = n
      exit find_sup_inflow_bc
    else if (trim(adjustl(uppercase(bc_input(n)%name))) == "" .and. &
             bc_idx == 0 ) then
      bc_idx = n
    end if
  end do find_sup_inflow_bc
  !
  ! Override user input with correct supersonic inflow boundary condition.
  !
  bc_input(bc_idx)%rho_static = default_rho
  bc_input(bc_idx)%vx         = default_vx
  bc_input(bc_idx)%vy         = default_vy
  bc_input(bc_idx)%vz         = default_vz
  bc_input(bc_idx)%p_static   = default_pres
  !
  bc_input(bc_idx)%t_static   = default_bc_input%t_static
  bc_input(bc_idx)%t_total    = default_bc_input%t_total
  bc_input(bc_idx)%p_total    = default_bc_input%p_total
  bc_input(bc_idx)%rho_total  = default_bc_input%rho_total
  bc_input(bc_idx)%alpha_aoa  = default_bc_input%alpha_aoa
  bc_input(bc_idx)%beta_aoa   = default_bc_input%beta_aoa
  bc_input(bc_idx)%wall_temp  = default_bc_input%wall_temp
  !
  ! Override reference conditions with supersonic freestream condition.
  !
  if (gam <= zero) gam = 1.4_wp
  if (rgasref <= zero) rgasref = rgas_cpg
  !
  alpha_aoaref = zero
  beta_aoaref  = zero
  rhoref       = default_rho
  pref         = default_pres
  tref         = pref / (rhoref*rgasref)
  uref         = default_vx
  aref         = sqrt( max( gam*rgasref*tref , rndoff ) )
  machref      = uref / aref
  reyref       = -one
  !
end subroutine initialize_shu_osher_problem
!
!###############################################################################
!> Solution initialization for Shu-Osher shock/entropy wave interaction problem.
!> Ref:
!> * Shu, C-W and Osher, S.,
!>   "Efficient implementation of essentially non-oscillatory
!>    shock-capturing schemes, II," J Comput Phys 1989; 83(1):32–78.
!>    http://dx.doi.org/10.1016/0021-9991(89)90222-2
!###############################################################################
pure function shu_osher_init(xyz) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : aref,uref,Lref
  use ovar, only : gam,rhoref,tref,pref
  use ovar, only : rgasref
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xyz
  !
  !.. Function Return Value ..
  real(wp), dimension(1:2+size(xyz)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: nc,mb,me,ne,nx
  real(wp) :: gm1,rhoa,rhoa2
  real(wp) :: rho,pres
  !
  !.. Local Arrays ..
  real(wp), dimension(1:3) :: vel
  !
continue
  !
  nc = 1
  mb = nc + 1
  me = nc + size(xyz)
  ne = me + 1
  !
  ! Some needed constants
  !
  gm1 = gam - one
  rhoa = rhoref * aref
  rhoa2 = rhoa * aref
  !
  !
  ! ------------------------------------------
  ! Dimensional flow properties on left state.
  ! ------------------------------------------
  !
  rho    = default_rho
  vel(1) = default_vx
  vel(2) = default_vy
  vel(3) = default_vz
  pres   = default_pres
  !
  !
  ! -------------------------------------------
  ! Dimensional flow properties on right state.
  ! -------------------------------------------
  !
  if ( xyz(1) >= -4.0_wp ) then
    rho    = one + 0.2_wp * sin( five * xyz(1) )  ! [kg/m3]
    vel(1) = zero  ! [m/s]
    vel(2) = zero  ! [m/s]
    vel(3) = zero  ! [m/s]
    pres   = one   ! [Pa]
  end if
  !
  ! Evaluate the dimensional initial conditions using conserved variables
  !
  return_value(nc)    = rho
  return_value(mb:me) = rho*vel(1:size(xyz))
  return_value(ne)    = pres/gm1 + half*rho*selfdot(vel)
  !
  ! Non-dimensionalize the initial conditions
  !
  return_value(nc)    = return_value(nc)    / rhoref
  return_value(mb:me) = return_value(mb:me) / rhoa
  return_value(ne)    = return_value(ne)    / rhoa2
  !
end function shu_osher_init
!
!###############################################################################
!
end module shu_osher_mod
