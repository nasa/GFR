pure function usp2v(uu) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: uu
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
end function usp2v
!
!###############################################################################
!
pure function usp2v_sp(uu) result(return_value)
  !
  !.. Use Statements ..
  use ovar,   only : gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: uu
  !
  !.. Function Result ..
  real(wp), dimension(1:size(uu)) :: return_value
  !
  !.. Local Scalars ..
  integer :: nq
  real(wp) :: gm1,q2
  !
continue
  !
  nq = size(uu)
  gm1 = gam - one
  !
  q2 = dot( uu(2:nq-1) , uu(2:nq-1) )
  !
  return_value(1)      = uu(1)
  return_value(2:nq-1) = uu(2:nq-1)/uu(1)
  return_value(nq)     = gm1 * (uu(nq) - half*q2/uu(1))
  !
end function usp2v_sp
!
!###############################################################################
!
