pure function selfdot_DP(xyz) result(return_value)
  real(dp), dimension(:), intent(in) :: xyz
  real(dp) :: return_value
  return_value = dot(xyz,xyz)
end function selfdot_DP
!
pure function selfdot_QP(xyz) result(return_value)
  real(qp), dimension(:), intent(in) :: xyz
  real(qp) :: return_value
  return_value = dot_product(xyz,xyz)
end function selfdot_QP
!
pure function rnorm2_DP(xyz) result(return_value)
  real(dp), dimension(:), intent(in) :: xyz
  real(dp) :: return_value
  return_value = sqrt(selfdot(xyz))
end function rnorm2_DP
!
pure function rnorm2_QP(xyz) result(return_value)
  real(qp), dimension(:), intent(in) :: xyz
  real(qp) :: return_value
  return_value = sqrt(selfdot(xyz))
end function rnorm2_QP
!
pure function unit_vector_DP(xyz) result(return_value)
  real(dp), dimension(:), intent(in) :: xyz
  real(dp), dimension(1:size(xyz)) :: return_value
  real(dp) :: grads,grdsnz,grdsrt
  grads = rnorm2(xyz)
  grdsnz = max(grads,real(epsilon(0.0_dp),kind=dp))
  grdsrt = grads / grdsnz
  return_value(:) = xyz(:) / grdsnz*grdsrt
end function unit_vector_DP
!
pure function unit_vector_QP(xyz) result(return_value)
  real(qp), dimension(:), intent(in) :: xyz
  real(qp), dimension(1:size(xyz)) :: return_value
  real(qp) :: grads,grdsnz,grdsrt
  grads = rnorm2(xyz)
  grdsnz = max(grads,real(epsilon(0.0_dp),kind=qp))
  grdsrt = grads / grdsnz
  return_value(:) = xyz(:) / grdsnz*grdsrt
end function unit_vector_QP
