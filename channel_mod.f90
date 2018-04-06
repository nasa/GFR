module channel_mod
  !
  ! Module containing procedures for the channel flow problem
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nee,nq
  !
  implicit none
  !
  private
  !
  !.. Public procedures
  public :: read_channel_init_file
  public :: channel_initialization
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
subroutine channel_source()
  !
  !.. Use Statements ..
  use geovar,      only : n_solpts
  use ovar,        only : channel_body_force
  use ovar,        only : rhoref,aref
  use flowvar,     only : residual,usp
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
