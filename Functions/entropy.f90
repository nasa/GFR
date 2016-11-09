pure function entropy_cv(uu,log_opt,use_pref_nd) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam, pref_nd, include_tke_in_energy
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: uu(:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: log_opt
  logical(lk), optional, intent(in) :: use_pref_nd
  !
  !.. Function Result ..
  real(wp), dimension(1:size(uu,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: n
  logical(lk)  :: use_log
  real(wp) :: gm1,q2,pres
  !
continue
  !
  use_log = .true.
  if (present(log_opt)) use_log = log_opt
  !
  gm1 = gam - one
  !
  do n = 1,size(uu,dim=2)
    !
    q2 = dot( uu(nmb:nme,n) , uu(nmb:nme,n) )
    pres = gm1 * (uu(nee,n) - half*q2/uu(nec,n))
    !
    if (include_tke_in_energy) pres = pres - gm1*uu(ntk,n)
    !
    if (present(use_pref_nd)) then
      pres = pres / merge(pref_nd,one,use_pref_nd)
    end if
    !
    if (use_log) then
      return_value(n) = log(pres) - gam*log(uu(nec,n))
    else
      return_value(n) = pres / (uu(nec,n)**gam)
    end if
    !
  end do
  !
end function entropy_cv
!
!###############################################################################
!
pure function entropy_cv_sp(uu,log_opt,use_pref_nd) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam, pref_nd, include_tke_in_energy
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: uu(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: log_opt
  logical(lk), optional, intent(in) :: use_pref_nd
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: gm1,q2,pres
  logical(lk)  :: use_log
  !
continue
  !
  use_log = .true.
  if (present(log_opt)) use_log = log_opt
  !
  gm1 = gam - one
  !
  q2 = dot( uu(nmb:nme) , uu(nmb:nme) )
  pres = gm1 * (uu(nee) - half*q2/uu(nec))
  !
  if (include_tke_in_energy) pres = pres - gm1*uu(ntk)
  !
  if (present(use_pref_nd)) then
    pres = pres / merge(pref_nd,one,use_pref_nd)
  end if
  !
  if (use_log) then
    return_value = log(pres) - gam*log(uu(nec))
  else
    return_value = pres / (uu(nec)**gam)
  end if
  !
end function entropy_cv_sp
!
!###############################################################################
!
pure function entropy_pv(vv,log_opt,use_pref_nd) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam, pref_nd
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: vv(:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: log_opt
  logical(lk), optional, intent(in) :: use_pref_nd
  !
  !.. Function Result ..
  real(wp), dimension(1:size(vv,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: n
  logical(lk)  :: use_log
  real(wp) :: pres
  !
continue
  !
  use_log = .true.
  if (present(log_opt)) use_log = log_opt
  !
  do n = 1,size(vv,dim=2)
    !
    pres = vv(nee,n)
    if (present(use_pref_nd)) then
      pres = pres / merge(pref_nd,one,use_pref_nd)
    end if
    !
    if (use_log) then
      return_value(n) = log(pres) - gam*log(vv(nec,n))
    else
      return_value(n) = pres / (vv(nec,n)**gam)
    end if
    !
  end do
  !
end function entropy_pv
!
!###############################################################################
!
pure function entropy_pv_sp(vv,log_opt,use_pref_nd) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam, pref_nd
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: vv(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: log_opt
  logical(lk), optional, intent(in) :: use_pref_nd
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  logical(lk)  :: use_log
  real(wp) :: pres
  !
continue
  !
  use_log = .true.
  if (present(log_opt)) use_log = log_opt
  !
  pres = vv(nee)
  if (present(use_pref_nd)) then
    pres = pres / merge(pref_nd,one,use_pref_nd)
  end if
  !
  if (use_log) then
    return_value = log(pres) - gam*log(vv(nec))
  else
    return_value = pres / (vv(nec)**gam)
  end if
  !
end function entropy_pv_sp
!
!###############################################################################
!
