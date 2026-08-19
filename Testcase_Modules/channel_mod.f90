module channel_mod
  !
  ! Module containing procedures for the channel flow problem
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nee,nq
  !
  implicit none
  !
  private
  !
  !.. Public procedures
  public :: read_channel_init_file
  public :: channel_initialization
  public :: pipe_initialization
  public :: channel_source
  !
  !.. Local Module Variables ..
  real(wp), allocatable, save :: init_sol(:,:)
  !
contains
!
!###############################################################################
!
subroutine read_channel_init_file()
  !
  !.. Use Statements ..
  use ovar, only : channel_init_file
  !
  !.. Local Parameters..
  character(len=*), parameter :: pname = "read_channel_init_file"
  !
continue
  !
end subroutine read_channel_init_file
!
!###############################################################################
!
pure function channel_initialization( xyz ) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : bc_in
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: xyz
  !
  !.. Function Result ..
  real(wp), dimension(1:size(xyz,dim=2),1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer :: k,np
  !
continue
  !
  np = size(xyz,dim=2)
  !
  if (allocated(init_sol)) then
    !
    !
  else
    !
    do k = 1,np
      return_value(k,:) = bc_in(0)%cv(:)
    end do
    !
  end if
  !
end function channel_initialization
!
!###############################################################################
!
pure function pipe_initialization( xyz ) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : bc_in
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: xyz
  !
  !.. Function Result ..
  real(wp), dimension(1:size(xyz,dim=2),1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer  :: k,np
  real(wp) :: rho,rhou,rhov,rhow,rhoe
  real(wp) :: x,y,z,r,t
  !
continue
  !
  np = size(xyz,dim=2)
  !
  if (allocated(init_sol)) then
    !
    !
  else
    !
    do k = 1,np
      !
      x = xyz(1,k)
      y = xyz(2,k)
      z = xyz(3,k)
      r = y*y + z*z
      t = atan2(z,y)
      t = t + merge(two*pi,zero,t<zero)
      !
      rho  = bc_in(0)%cv(nec)
      rhou = bc_in(0)%cv(nmx)
      rhov = bc_in(0)%cv(nmy)
      rhow = bc_in(0)%cv(nmz)
      rhoe = bc_in(0)%cv(nee)
      !
      return_value(k,:) = bc_in(0)%cv(:)
      !
    end do
    !
  end if
  !
end function pipe_initialization
!
!###############################################################################
!
subroutine channel_source()
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_solpts
  use ovar,        only : channel_body_force
  use ovar,        only : rhoref,aref
  use flowvar_mod, only : residual,usp
  use metrics_mod, only : jac
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: bfx,hc,vx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "channel_source"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  bfx = channel_body_force / (rhoref*aref*aref)
  !
  do n = 1,n_solpts
    !
    vx = usp(n,nmx) / usp(n,nec)
    !
    hc = -bfx * jac(n)
    !
    residual(n,nmx) = residual(n,nmx) + hc
    residual(n,nee) = residual(n,nee) + hc * vx
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine channel_source
!
!###############################################################################
!
end module channel_mod
!
!###############################################################################
!
