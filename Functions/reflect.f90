pure function reflect_DP(node,face_center,normal) result(return_value)
  !
  real(dp), dimension(:), intent(in) :: node
  real(dp), dimension(:), intent(in) :: face_center
  real(dp), dimension(:), intent(in) :: normal
  !
  real(dp), dimension(1:size(node)) :: return_value
  !
  real(dp), dimension(1:size(node)) :: ds
  !
continue
  !
  ds = face_center - node
  !
  return_value = node + 2.0_dp*dot(ds,normal)*normal
  !
end function reflect_DP
!
pure function reflect_QP(node,face_center,normal) result(return_value)
  !
  real(qp), dimension(:), intent(in) :: node
  real(qp), dimension(:), intent(in) :: face_center
  real(qp), dimension(:), intent(in) :: normal
  !
  real(qp), dimension(1:size(node)) :: return_value
  !
  real(qp), dimension(1:size(node)) :: ds
  !
continue
  !
  ds = face_center - node
  !
  return_value = node + 2.0_qp*dot_product(ds,normal)*normal
  !
end function reflect_QP
