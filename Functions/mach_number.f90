pure function mach_number_cv(uu) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam, include_tke_in_energy
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: uu
  !
  !.. Function Result ..
  real(wp), dimension(1:size(uu,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer :: n
  real(wp) :: gm1,q2,pres,mach
  !
continue
  !
  gm1 = gam - one
  !
  do n = 1,size(uu,dim=2)
    !
    q2 = dot( uu(nmb:nme,n) , uu(nmb:nme,n) )
    !
    pres = gm1 * (uu(nee,n) - half*q2/uu(nec,n))
    !
    if (include_tke_in_energy) pres = pres - gm1*uu(ntk,n)
    !
    mach = q2 / (gam*pres*uu(nec,n))
    !
    return_value(n) = sqrt( max(mach,zero) )
    !
  end do
  !
end function mach_number_cv
!
!###############################################################################
!
pure function mach_number_cv_sp(uu) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam, include_tke_in_energy
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: uu
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: gm1,q2,pres,mach
  !
continue
  !
  gm1 = gam - one
  !
  q2 = dot( uu(nmb:nme) , uu(nmb:nme) )
  !
  pres = gm1 * (uu(nee) - half*q2/uu(nec))
  !
  if (include_tke_in_energy) pres = pres - gm1*uu(ntk)
  !
  mach = q2 / (gam*pres*uu(nec))
  !
  return_value = sqrt( max(mach,zero) )
  !
end function mach_number_cv_sp
!
!###############################################################################
!
pure function mach_number_pv(vv) result(return_value)
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
  real(wp) :: q2,mach
  !
continue
  !
  do n = 1,size(vv,dim=2)
    !
    q2 = dot( vv(nmb:nme,n) , vv(nmb:nme,n) )
    !
    mach = q2 / (gam*vv(nee,n)*vv(nec,n))
    !
    return_value(n) = sqrt( max(mach,zero) )
    !
  end do
  !
end function mach_number_pv
!
!###############################################################################
!
pure function mach_number_pv_sp(vv) result(return_value)
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
  !.. Local Scalars ..
  real(wp) :: q2,mach
  !
continue
  !
  q2 = dot( vv(nmb:nme) , vv(nmb:nme) )
  !
  mach = q2 / (gam*vv(nee)*vv(nec))
  !
  return_value = sqrt( max(mach,zero) )
  !
end function mach_number_pv_sp
!
!###############################################################################
!
