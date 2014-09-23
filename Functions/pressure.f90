pure function pressure_cv(uu) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: uu
  !
  !.. Function Result ..
  real(wp), dimension(1:size(uu,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer :: n,nq
  real(wp) :: gm1,q2
  !
continue
  !
  nq = size(uu,dim=1)
  gm1 = gam - one
  !
  do n = 1,size(uu,dim=2)
    !
    q2 = dot_product( uu(2:nq-1,n) , uu(2:nq-1,n) )
    !
    return_value(n) = gm1 * (uu(nq,n) - half*q2/uu(1,n))
    !
  end do
  !
end function pressure_cv
!
!###############################################################################
!
pure function pressure_cv_sp(uu) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: uu
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  integer  :: nq
  real(wp) :: gm1,q2
  !
continue
  !
  nq = size(uu)
  gm1 = gam - one
  !
  q2 = dot_product( uu(2:nq-1) , uu(2:nq-1) )
  !
  return_value = gm1 * (uu(nq) - half*q2/uu(1))
  !
end function pressure_cv_sp
!
!###############################################################################
!
