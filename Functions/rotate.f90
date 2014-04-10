pure function rotate(face_node1,face_node2,node_off_face) result(return_value)
  !
  real(vbp), dimension(:), intent(in) :: face_node1
  real(vbp), dimension(:), intent(in) :: face_node2
  real(vbp), dimension(:), intent(in) :: node_off_face
  !
  real(vbp), dimension(1:size(face_node1)) :: return_value
  !
continue
  !
  return_value = face_node1 + face_node2 - node_off_face
  !
end function rotate
