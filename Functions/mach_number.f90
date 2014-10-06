pure function mach_number_cv(uu) result(return_value)
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
  integer :: n,nq
  real(wp) :: gm1,q2,pres,mach
  !
continue
  !
  nq = size(uu,dim=1)
  gm1 = gam - one
  !
  do n = 1,size(uu,dim=2)
    !
    q2 = dot( uu(2:nq-1,n) , uu(2:nq-1,n) )
    !
    pres = gm1 * (uu(nq,n) - half*q2/uu(1,n))
    !
    mach = q2 / (gam*pres*uu(1,n))
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
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: uu
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  integer :: nq
  real(wp) :: gm1,q2,pres,mach
  !
continue
  !
  nq = size(uu)
  gm1 = gam - one
  !
  q2 = dot( uu(2:nq-1) , uu(2:nq-1) )
  !
  pres = gm1 * (uu(nq) - half*q2/uu(1))
  !
  mach = q2 / (gam*pres*uu(1))
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
  integer :: n,nq
  real(wp) :: q2,mach
  !
continue
  !
  nq = size(vv,dim=1)
  !
  do n = 1,size(vv,dim=2)
    !
    q2 = dot( vv(2:nq-1,n) , vv(2:nq-1,n) )
    !
    mach = q2 / (gam*vv(nq,n)*vv(1,n))
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
  integer :: nq
  real(wp) :: q2,mach
  !
continue
  !
  nq = size(vv)
  !
  q2 = dot( vv(2:nq-1) , vv(2:nq-1) )
  !
  mach = q2 / (gam*vv(nq)*vv(1))
  !
  return_value = sqrt( max(mach,zero) )
  !
end function mach_number_pv_sp
!
!###############################################################################
!
