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
  use geovar,      only : ncell,cell
  use ovar,        only : channel_body_force
  use ovar,        only : rhoref,aref
  use flowvar,     only : residual,usp
  use metrics_mod, only : metrics_dt
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2
  !
  real(wp) :: bfx,hc,vx
  !
  !.. Local Parameters ..
  real(wp), parameter :: add_or_subtract = +one
 !real(wp), parameter :: add_or_subtract = -one
  !
 !integer, parameter :: jac_exp =  0
  integer, parameter :: jac_exp =  1
 !integer, parameter :: jac_exp = -1
  !
  character(len=*), parameter :: pname = "channel_source"
  !
continue
  !
  bfx = channel_body_force / (rhoref*aref*aref)
  !
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    !
    do n = n1,n2
      !
      vx = usp(n,nmx) / usp(n,nec)
      !
      hc = -bfx / (metrics_dt(nc)%jac(n-n1+1)**jac_exp)
      !
      residual(n,nmx) = residual(n,nmx) + add_or_subtract * hc
      residual(n,nee) = residual(n,nee) + add_or_subtract * hc * vx
      !
    end do
    !
  end do
  !
end subroutine channel_source
!
!###############################################################################
!
end module channel_mod
!
!###############################################################################
!
