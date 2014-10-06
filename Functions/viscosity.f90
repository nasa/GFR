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
  integer  :: n,nq
  real(wp) :: cv_inv,q2,temp,tdim,mudim,muref_inv
  !
continue
  !
  nq = size(uu,dim=1)
  !
  muref_inv = one/rhoref/aref
  cv_inv = (gam-one)/rgasref_nd
  !
  do n = 1,size(uu,dim=2)
    !
    q2 = dot( uu(2:nq-1,n) , uu(2:nq-1,n) )
    !
    temp = cv_inv * (uu(nq,n) - half*q2/uu(1,n)) / uu(1,n)
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
  integer :: nq
  real(wp) :: cv_inv,q2,temp,tdim,mudim
  !
continue
  !
  nq = size(uu)
  cv_inv = (gam-one)/rgasref_nd
  !
  q2 = dot( uu(2:nq-1) , uu(2:nq-1) )
  !
  temp = cv_inv * (uu(nq) - half*q2/uu(1)) / uu(1)
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
  integer  :: n,nq
  real(wp) :: tdim,mudim,muref_inv,rgas_inv
  !
continue
  !
  nq = size(vv,dim=1)
  muref_inv = one/rhoref/aref
  rgas_inv = one/rgasref_nd
  !
  do n = 1,size(vv,dim=2)
    tdim = (vv(nq,n) * rgas_inv / vv(1,n)) * tref
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
  integer  :: nq
  real(wp) :: tdim,mudim
  !
continue
  !
  nq = size(vv)
  !
  tdim = (vv(nq) / rgasref_nd / vv(1)) * tref
  mudim = suth_c1 * tdim * sqrt(tdim) / (tdim + suth_sref)
  return_value = mudim / rhoref / aref
  !
end function viscosity_pv_sp
!
!###############################################################################
!
