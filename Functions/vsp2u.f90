pure function vsp2u(vv) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: vv
  !
  !.. Function Result ..
  real(wp), dimension(1:size(vv,dim=1),1:size(vv,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  do n = 1,size(vv,dim=2)
    return_value(:,n) = vsp2u_sp( vv(:,n) )
  end do
  !
end function vsp2u
!
!###############################################################################
!
pure function vsp2u_sp(vv) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam, include_tke_in_energy
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: vv
  !
  !.. Function Result ..
  real(wp), dimension(1:size(vv)) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: gm1,q2
  !
continue
  !
  gm1 = gam - one
  !
  q2 = dot( vv(nmb:nme) , vv(nmb:nme) )
  !
  return_value(nec)     = vv(nec)
  return_value(nmb:nme) = vv(nec)*vv(nmb:nme)
  return_value(nee)     = vv(nee)/gm1 + half*vv(nec)*q2
  !
  if (include_tke_in_energy) then
    return_value(ntk) = vv(nec) * vv(ntk)
    return_value(ntl) = vv(nec) * vv(ntl)
    return_value(nee) = return_value(nee) + return_value(ntk)
  end if
  !
end function vsp2u_sp
!
!###############################################################################
!
