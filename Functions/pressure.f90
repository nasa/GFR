pure function pressure_cv(uu) result(return_value)
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
  real(wp) :: gm1,q2
  !
continue
  !
  gm1 = gam - one
  !
  do n = 1,size(uu,dim=2)
    !
    q2 = dot( uu(nmb:nme,n) , uu(nmb:nme,n) )
    !
    return_value(n) = gm1 * (uu(nee,n) - half*q2/uu(nec,n))
    !
  end do
  !
  if (include_tke_in_energy) then
    do n = 1,size(uu,dim=2)
      return_value(n) = return_value(n) - gm1*uu(ntk,n)
    end do
  end if
  !
end function pressure_cv
!
!###############################################################################
!
pure function pressure_cv_sp(uu) result(return_value)
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
  real(wp) :: gm1,q2
  !
continue
  !
  gm1 = gam - one
  !
  q2 = dot( uu(nmb:nme) , uu(nmb:nme) )
  !
  return_value = gm1 * (uu(nee) - half*q2/uu(nec))
  !
  if (include_tke_in_energy) return_value = return_value - gm1*uu(ntk)
  !
end function pressure_cv_sp
!
!###############################################################################
!
