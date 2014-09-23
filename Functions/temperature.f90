pure function temperature_cv(uu) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: uu
  !
  !.. Function Result ..
  real(wp), dimension(1:size(uu,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: n,nq
  real(wp) :: q2,cv_inv
  !
continue
  !
  nq = size(uu,dim=1)
  cv_inv = (gam-one)/rgasref_nd
  !
  do n = 1,size(uu,dim=2)
    !
    q2 = dot_product( uu(2:nq-1,n) , uu(2:nq-1,n) )
    !
    return_value(n) = cv_inv * (uu(nq,n) - half*q2/uu(1,n)) / uu(1,n)
    !
  end do
  !
end function temperature_cv
!
!###############################################################################
!
pure function temperature_cv_sp(uu) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: uu
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  integer :: nq
  real(wp) :: q2,cv_inv
  !
continue
  !
  nq = size(uu)
  cv_inv = (gam-one)/rgasref_nd
  !
  q2 = dot_product( uu(2:nq-1) , uu(2:nq-1) )
  !
  return_value = cv_inv * (uu(nq) - half*q2/uu(1)) / uu(1)
  !
end function temperature_cv_sp
!
!###############################################################################
!
pure function temperature_pv(vv) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vv
  !
  !.. Function Result ..
  real(wp), dimension(1:size(vv,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: n,nq
  real(wp) :: rgas_inv
  !
continue
  !
  nq = size(vv,dim=1)
  rgas_inv = one/rgasref_nd
  !
  do n = 1,size(vv,dim=2)
    return_value(n) = vv(nq,n) * rgas_inv / vv(1,n)
  end do
  !
end function temperature_pv
!
!###############################################################################
!
pure function temperature_pv_sp(vv) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: vv
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  integer :: nq
  !
continue
  !
  nq = size(vv)
  return_value = vv(nq)/vv(1)/rgasref_nd
  !
end function temperature_pv_sp
!
!###############################################################################
!
