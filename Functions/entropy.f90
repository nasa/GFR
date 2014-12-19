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
  integer :: n
  real(wp) :: gm1,q2,pres
  !
continue
  !
  gm1 = gam - one
  !
  do n = 1,size(uu,dim=2)
    !
    q2 = dot( uu(nmb:nme,n) , uu(nmb:nme,n) )
    pres = gm1 * (uu(nee,n) - half*q2/uu(nec,n))
    !
    return_value(n) = log(pres) - gam*log(uu(nec,n))
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
  real(wp) :: gm1,q2,pres
  !
continue
  !
  gm1 = gam - one
  !
  q2 = dot( uu(nmb:nme) , uu(nmb:nme) )
  pres = gm1 * (uu(nee) - half*q2/uu(nec))
  !
  return_value = log(pres) - gam*log(uu(nec))
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
  integer :: n
  !
continue
  !
  do n = 1,size(vv,dim=2)
    return_value(n) = log(vv(nee,n)) - gam*log(vv(nec,n)) 
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
  return_value = log(vv(nee)) - gam*log(vv(nec))
  !
end function entropy_pv_sp
!
!###############################################################################
!
