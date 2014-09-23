pure function usp2utd(uu,jacobian) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:),   intent(in) :: jacobian
  real(wp), dimension(:,:), intent(in) :: uu
  !
  !.. Function Results ..
  real(wp), dimension(1:size(uu,dim=1),1:size(uu,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  do n = 1,size(uu,dim=2)
    return_value(:,n) = uu(:,n)*jacobian(n)
  end do
  !
end function usp2utd
!
!###############################################################################
!
