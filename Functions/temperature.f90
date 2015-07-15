pure function temperature_cv(uu) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam, rgasref_nd, include_tke_in_energy
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: uu
  !
  !.. Function Result ..
  real(wp), dimension(1:size(uu,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: q2,cv_inv
  !
continue
  !
  cv_inv = (gam-one)/rgasref_nd
  !
  do n = 1,size(uu,dim=2)
    !
    q2 = dot( uu(nmb:nme,n) , uu(nmb:nme,n) )
    !
    return_value(n) = cv_inv * (uu(nee,n) - half*q2/uu(nec,n)) / uu(nec,n)
    !
  end do
  !
  if (include_tke_in_energy) then
    do n = 1,size(uu,dim=2)
      return_value(n) = return_value(n) - cv_inv * uu(ntk,n) / uu(nec,n)
    end do
  end if
  !
end function temperature_cv
!
!###############################################################################
!
pure function temperature_cv_sp(uu) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam, rgasref_nd, include_tke_in_energy
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: uu
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: q2,cv_inv
  !
continue
  !
  cv_inv = (gam-one)/rgasref_nd
  !
  q2 = dot( uu(nmb:nme) , uu(nmb:nme) )
  !
  return_value = cv_inv * (uu(nee) - half*q2/uu(nec)) / uu(nec)
  !
  if (include_tke_in_energy) then
    return_value = return_value - cv_inv * uu(ntk) / uu(nec)
  end if
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
  integer  :: n
  real(wp) :: rgas_inv
  !
continue
  !
  rgas_inv = one/rgasref_nd
  !
  do n = 1,size(vv,dim=2)
    return_value(n) = vv(nee,n) * rgas_inv / vv(nec,n)
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
continue
  !
  return_value = vv(nee)/vv(nec)/rgasref_nd
  !
end function temperature_pv_sp
!
!###############################################################################
!
