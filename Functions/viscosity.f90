pure function viscosity_cv(uu) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : tref,rhoref,aref
  use ovar, only : suth_c1,suth_sref
  use ovar, only : gam,rgasref_nd
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: uu
  !
  !.. Function Result ..
  real(wp), dimension(1:size(uu,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: cv_inv,q2,temp,tdim,mudim,muref_inv
  !
continue
  !
  muref_inv = one/rhoref/aref
  cv_inv = (gam-one)/rgasref_nd
  !
  do n = 1,size(uu,dim=2)
    !
    q2 = dot( uu(nmb:nme,n) , uu(nmb:nme,n) )
    !
    temp = cv_inv * (uu(nee,n) - half*q2/uu(nec,n)) / uu(nec,n)
    !
    tdim = temp * tref
    mudim = suth_c1 * tdim * sqrt(tdim) / (tdim + suth_sref)
    return_value(n) = mudim * muref_inv
    !
  end do
  !
end function viscosity_cv
!
!###############################################################################
!
pure function viscosity_cv_sp(uu) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : tref,rhoref,aref
  use ovar, only : suth_c1,suth_sref
  use ovar, only : gam,rgasref_nd
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: uu
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: cv_inv,q2,temp,tdim,mudim
  !
continue
  !
  cv_inv = (gam-one)/rgasref_nd
  !
  q2 = dot( uu(nmb:nme) , uu(nmb:nme) )
  !
  temp = cv_inv * (uu(nee) - half*q2/uu(nec)) / uu(nec)
  !
  tdim = temp * tref
  mudim = suth_c1 * tdim * sqrt(tdim) / (tdim + suth_sref)
  return_value = mudim / rhoref / aref
  !
end function viscosity_cv_sp
!
!###############################################################################
!
pure function viscosity_pv(vv) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : tref,rhoref,aref,rgasref_nd
  use ovar, only : suth_c1,suth_sref
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vv
  !
  !.. Function Result ..
  real(wp), dimension(1:size(vv,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: tdim,mudim,muref_inv,rgas_inv
  !
continue
  !
  muref_inv = one/rhoref/aref
  rgas_inv = one/rgasref_nd
  !
  do n = 1,size(vv,dim=2)
    tdim = (vv(nee,n) * rgas_inv / vv(nec,n)) * tref
    mudim = suth_c1 * tdim * sqrt(tdim) / (tdim + suth_sref)
    return_value(n) = mudim * muref_inv
  end do
  !
end function viscosity_pv
!
!###############################################################################
!
pure function viscosity_pv_sp(vv) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : tref,rhoref,aref,rgasref_nd
  use ovar, only : suth_c1,suth_sref
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: vv
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: tdim,mudim
  !
continue
  !
  tdim = (vv(nee) / rgasref_nd / vv(nec)) * tref
  mudim = suth_c1 * tdim * sqrt(tdim) / (tdim + suth_sref)
  return_value = mudim / rhoref / aref
  !
end function viscosity_pv_sp
!
!###############################################################################
!
