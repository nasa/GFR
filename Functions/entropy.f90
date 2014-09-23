pure function entropy_cv(uu) result(return_value)
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
  real(wp) :: gm1,q2,pres
  !
continue
  !
  nq = size(uu,dim=1)
  gm1 = gam - one
  !
  do n = 1,size(uu,dim=2)
    !
    q2 = dot_product( uu(2:nq-1,n) , uu(2:nq-1,n) )
    pres = gm1 * (uu(nq,n) - half*q2/uu(1,n))
    !
    return_value(n) = log(pres) - gam*log(uu(1,n))
    !
  end do
  !
end function entropy_cv
!
!###############################################################################
!
pure function entropy_cv_sp(uu) result(return_value)
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
  integer :: nq
  real(wp) :: gm1,q2,pres
  !
continue
  !
  nq = size(uu)
  gm1 = gam - one
  !
  q2 = dot_product( uu(2:nq-1) , uu(2:nq-1) )
  pres = gm1 * (uu(nq) - half*q2/uu(1))
  !
  return_value = log(pres) - gam*log(uu(1))
  !
end function entropy_cv_sp
!
!###############################################################################
!
pure function entropy_pv(vv) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vv
  !
  !.. Function Result ..
  real(wp), dimension(1:size(vv,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer :: n,nq
  !
continue
  !
  nq = size(vv,dim=1)
  !
  do n = 1,size(vv,dim=2)
    return_value(n) = log(vv(nq,n)) - gam*log(vv(1,n)) 
  end do
  !
end function entropy_pv
!
!###############################################################################
!
pure function entropy_pv_sp(vv) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: vv
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
continue
  !
  return_value = log(vv(size(vv))) - gam*log(vv(1))
  !
end function entropy_pv_sp
!
!###############################################################################
!
