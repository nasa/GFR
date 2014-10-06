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
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: vv
  !
  !.. Function Result ..
  real(wp), dimension(1:size(vv)) :: return_value
  !
  !.. Local Scalars ..
  integer :: nq
  real(wp) :: gm1,q2
  !
continue
  !
  nq = size(vv)
  gm1 = gam - one
  !
  q2 = dot( vv(2:nq-1) , vv(2:nq-1) )
  !
  return_value(1) = vv(1)
  return_value(2:nq-1) = vv(1)*vv(2:nq-1)
  return_value(nq) = vv(nq)/gm1 + half*vv(1)*q2
  !
end function vsp2u_sp
!
!###############################################################################
!
