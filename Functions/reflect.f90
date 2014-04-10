pure function reflect(node,face_center,normal) result(return_value)
  !
  real(vbp), dimension(:), intent(in) :: node
  real(vbp), dimension(:), intent(in) :: face_center
  real(vbp), dimension(:), intent(in) :: normal
  !
  real(vbp), dimension(1:size(node)) :: return_value
  !
  real(vbp), dimension(1:size(node)) :: ds
  !
continue
  !
  ds = face_center - node
  !
  return_value = node + two*dot_product(ds,normal)*normal
  !
end function reflect
