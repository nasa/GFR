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
  use ovar,       only : bc_input_t,bc_input
  use ovar,       only : rgasref,gam
  use ovar,       only : alpha_aoaref,beta_aoaref
  use ovar,       only : rhoref,pref,tref,uref,aref,machref,reyref
  use thermo_mod, only : rgas_cpg
  !
  !.. Local Variables ..
  integer :: n,bc_idx,ierr
  !
  !.. Allocatable Scalars ...
  type(bc_input_t), allocatable :: tmp_bc_input
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_shu_osher_problem"
  !
continue
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
  ! Return if bc_input is not allocated
  !
  if (.not. allocated(bc_input)) then
    return
  end if
  !
  ! Find the index in bc_input for the supersonic inflow.
  ! If missing, find the first empty name.
  !
  bc_idx = 0
  find_sup_inflow_bc: &
  do n = 1,size(bc_input)
    if (trim(adjustl(uppercase(bc_input(n)%name))) == "INFLOW") then
      bc_idx = n
      exit find_sup_inflow_bc
    else if (trim(adjustl(uppercase(bc_input(n)%name))) == "" .and. &
             bc_idx == 0 ) then
      bc_idx = n
    end if
  end do find_sup_inflow_bc
  !
  ! Return if bc_idx is out of the bounds of bc_input,
  ! e.g., bc_idx has a value 0, or bc_input has a size of 0
  ! which should mean bc_input is allocated as bc_input(1:0).
  ! NOTE: If the flow conditions of bc_input are still the
  !       default values, it will set the boundary flow
  !       conditions to the freestream/reference conditions.
  !
  if (bc_idx < lbound(bc_input,dim=1) .or. &
      bc_idx > ubound(bc_input,dim=1)) then
    return
  end if
  !
  ! Allocate the tmp_bc_input scalar with default bc_input_t conditions
  !
  allocate ( tmp_bc_input , stat=ierr )
  !
  ! Save the non-flow-conditions part of bc_input to tmp_bc_input
  !
  call tmp_bc_input%copy_nonflow_conditions( bc_input(bc_idx) )
  !
  ! Set the inflow conditions for the Shu-Osher problem
  !
  tmp_bc_input%rho_static = default_rho
  tmp_bc_input%vx         = default_vx
  tmp_bc_input%vy         = default_vy
  tmp_bc_input%vz         = default_vz
  tmp_bc_input%p_static   = default_pres
  !
  ! Copy tmp_bc_input into bc_input(bc_idx)
  !
  bc_input(bc_idx) = tmp_bc_input
  !
  ! Deallocate tmp_bc_input
  !
  if (allocated(tmp_bc_input)) deallocate ( tmp_bc_input , stat=ierr )
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
 !if ( xyz(1) >= -3.0_wp ) then
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
