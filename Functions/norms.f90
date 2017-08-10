pure function selfdot_SP(xyz) result(return_value)
  real(SP), dimension(:), intent(in) :: xyz
  real(SP) :: return_value
  return_value = dot(xyz,xyz)
end function selfdot_SP
!
pure function selfdot_DP(xyz) result(return_value)
  real(DP), dimension(:), intent(in) :: xyz
  real(DP) :: return_value
  return_value = dot(xyz,xyz)
end function selfdot_DP
!
pure function selfdot_QP(xyz) result(return_value)
  real(QP), dimension(:), intent(in) :: xyz
  real(QP) :: return_value
  return_value = dot_product(xyz,xyz)
end function selfdot_QP
!
pure function rnorm2_SP(xyz) result(return_value)
  real(SP), dimension(:), intent(in) :: xyz
  real(SP) :: return_value
  return_value = sqrt(selfdot(xyz))
end function rnorm2_SP
!
pure function rnorm2_DP(xyz) result(return_value)
  real(DP), dimension(:), intent(in) :: xyz
  real(DP) :: return_value
  return_value = sqrt(selfdot(xyz))
end function rnorm2_DP
!
pure function rnorm2_QP(xyz) result(return_value)
  real(QP), dimension(:), intent(in) :: xyz
  real(QP) :: return_value
  return_value = sqrt(selfdot(xyz))
end function rnorm2_QP
!
pure function unit_vector_SP(xyz) result(return_value)
  real(SP), dimension(:), intent(in) :: xyz
  real(SP), dimension(1:size(xyz)) :: return_value
  real(SP) :: grads,grdsnz,grdsrt
  grads = norm2(xyz)
  grdsnz = max(grads,real(epsilon(0.0_SP),kind=SP))
  grdsrt = grads / grdsnz
  return_value(:) = xyz(:) / grdsnz*grdsrt
end function unit_vector_SP
!
pure function unit_vector_DP(xyz) result(return_value)
  real(DP), dimension(:), intent(in) :: xyz
  real(DP), dimension(1:size(xyz)) :: return_value
  real(DP) :: grads,grdsnz,grdsrt
  grads = norm2(xyz)
  grdsnz = max(grads,real(epsilon(0.0_DP),kind=DP))
  grdsrt = grads / grdsnz
  return_value(:) = xyz(:) / grdsnz*grdsrt
end function unit_vector_DP
!
pure function unit_vector_QP(xyz) result(return_value)
  real(QP), dimension(:), intent(in) :: xyz
  real(QP), dimension(1:size(xyz)) :: return_value
  real(QP) :: grads,grdsnz,grdsrt
  grads = norm2(xyz)
  grdsnz = max(grads,real(epsilon(0.0_DP),kind=QP))
  grdsrt = grads / grdsnz
  return_value(:) = xyz(:) / grdsnz*grdsrt
end function unit_vector_QP
