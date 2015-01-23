pure function node_center_SP(xyz_nodes) result(return_value)
  !.. Formal Arguments ..
  real(sp), dimension(:,:), intent(in) :: xyz_nodes
  !.. Function Result ..
  real(sp), dimension(1:size(xyz_nodes,dim=1)) :: return_value
  !.. Local Scalars ..
  real(sp) :: denom
  !.. Local Constants ..
  real(sp), parameter :: one = 1.0_sp
continue
  ! Count the number of nodes present
  denom = real( size(xyz_nodes,dim=2) , kind=sp )
  ! Prevent divide by zero if for some reason zero nodes are present
  denom = max(denom,one)
  ! Average the coordinates of the nodes
  return_value = sum(xyz_nodes,dim=2) / denom
end function node_center_SP
!
!###############################################################################
!
pure function node_center_DP(xyz_nodes) result(return_value)
  !.. Formal Arguments ..
  real(dp), dimension(:,:), intent(in) :: xyz_nodes
  !.. Function Result ..
  real(dp), dimension(1:size(xyz_nodes,dim=1)) :: return_value
  !.. Local Scalars ..
  real(dp) :: denom
  !.. Local Constants ..
  real(dp), parameter :: one = 1.0_dp
continue
  ! Count the number of nodes present
  denom = real( size(xyz_nodes,dim=2) , kind=dp )
  ! Prevent divide by zero if for some reason zero nodes are present
  denom = max(denom,one)
  ! Average the coordinates of the nodes
  return_value = sum(xyz_nodes,dim=2) / denom
end function node_center_DP
!
!###############################################################################
!
pure function node_center_QP(xyz_nodes) result(return_value)
  !.. Formal Arguments ..
  real(qp), dimension(:,:), intent(in) :: xyz_nodes
  !.. Function Result ..
  real(qp), dimension(1:size(xyz_nodes,dim=1)) :: return_value
  !.. Local Scalars ..
  real(qp) :: denom
  !.. Local Constants ..
  real(qp), parameter :: one = 1.0_qp
continue
  ! Count the number of nodes present
  denom = real( size(xyz_nodes,dim=2) , kind=qp )
  ! Prevent divide by zero if for some reason zero nodes are present
  denom = max(denom,one)
  ! Average the coordinates of the nodes
  return_value = sum(xyz_nodes,dim=2) / denom
end function node_center_QP
