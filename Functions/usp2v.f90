pure function usp2v(uu,swap_t_for_rho) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: uu
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: swap_t_for_rho
  !
  !.. Function Result ..
  real(wp), dimension(1:size(uu,dim=1),1:size(uu,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  do n = 1,size(uu,dim=2)
    return_value(:,n) = usp2v_sp( uu(:,n) )
  end do
  !
  if (present(swap_t_for_rho)) then
    if (swap_t_for_rho) then
      return_value(nec,:) = return_value(nee,:)/(return_value(nec,:)*rgasref_nd)
    end if
  end if
  !
end function usp2v
!
!###############################################################################
!
pure function usp2v_sp(uu,swap_t_for_rho) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam, rgasref_nd, include_tke_in_energy
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: uu
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: swap_t_for_rho
  !
  !.. Function Result ..
  real(wp), dimension(1:size(uu)) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: gm1,q2
  !
continue
  !
  gm1 = gam - one
  !
  q2 = dot( uu(nmb:nme) , uu(nmb:nme) )
  !
  return_value(nec)     = uu(nec)
  return_value(nmb:nme) = uu(nmb:nme)/uu(nec)
  return_value(nee)     = gm1 * (uu(nee) - half*q2/uu(nec))
  !
  if (include_tke_in_energy) then
    return_value(nee) = return_value(nee) - gm1*uu(ntk)
    return_value(ntk) = uu(ntk)/uu(nec)
    return_value(ntl) = uu(ntl)/uu(nec)
  end if
  !
  if (present(swap_t_for_rho)) then
    if (swap_t_for_rho) then
      return_value(nec) = return_value(nee) / (return_value(nec) * rgasref_nd)
    end if
  end if
  !
end function usp2v_sp
!
!###############################################################################
!
