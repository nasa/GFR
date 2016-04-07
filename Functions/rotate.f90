pure function rotate_DP(face_node1,face_node2,node_off_face) &
                 result(return_value)
  !
  real(dp), dimension(:), intent(in) :: face_node1
  real(dp), dimension(:), intent(in) :: face_node2
  real(dp), dimension(:), intent(in) :: node_off_face
  !
  real(dp), dimension(1:size(face_node1)) :: return_value
  !
continue
  !
  return_value = face_node1 + face_node2 - node_off_face
  !
end function rotate_DP
!
pure function rotate_QP(face_node1,face_node2,node_off_face) &
                 result(return_value)
  !
  real(qp), dimension(:), intent(in) :: face_node1
  real(qp), dimension(:), intent(in) :: face_node2
  real(qp), dimension(:), intent(in) :: node_off_face
  !
  real(qp), dimension(1:size(face_node1)) :: return_value
  !
continue
  !
  return_value = face_node1 + face_node2 - node_off_face
  !
end function rotate_QP
