module species_mod
  !
  use module_kind_types
  use eqn_idx, only : ncb,nce,nec
  !
  implicit none
  !
  private
  !
  public :: compute_mixture_density
  !
contains
!
!###############################################################################
!
subroutine compute_mixture_density(n1,n2,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nsp
  real(wp) :: rho
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_mixture_density"
  !
continue
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(usp)
  do n = n1,n2
    !
    rho = zero
    !$ACC LOOP VECTOR REDUCTION(+:rho)
    do m = ncb,nce
      rho = rho + usp(n,m)
    end do
   !rho = sum(usp(n,ncb:nce))
    !
    usp(nec,n) = rho
    !
  end do
  !
end subroutine compute_mixture_density
!
!###############################################################################
!
end module species_mod
