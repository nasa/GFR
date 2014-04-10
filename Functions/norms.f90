pure function selfdot(xyz) result(return_value)
  real(vbp), dimension(:), intent(in) :: xyz
  real(vbp) :: return_value
  return_value = dot_product(xyz,xyz)
end function selfdot
!
pure function rnorm2(xyz) result(return_value)
  real(vbp), dimension(:), intent(in) :: xyz
  real(vbp) :: return_value
  return_value = sqrt(selfdot(xyz))
end function rnorm2
!
pure function unit_vector(xyz) result(return_value)
  real(vbp), dimension(:), intent(in) :: xyz
  real(vbp), dimension(1:size(xyz)) :: return_value
  real(vbp) :: grads,grdsnz,grdsrt
  grads = rnorm2(xyz)
  grdsnz = max(grads,rndoff)
  grdsrt = grads / grdsnz
  return_value(:) = xyz(:) / grdsnz*grdsrt
end function unit_vector
