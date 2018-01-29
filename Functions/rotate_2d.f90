pure function rotate_2d_DP(xy,theta) result(return_value)
  real(dp),               intent(in) :: theta
  real(dp), dimension(:), intent(in) :: xy
  real(dp), dimension(1:size(xy)) :: return_value
  return_value(1) = xy(1)*cos(theta) - xy(2)*sin(theta)
  return_value(2) = xy(1)*sin(theta) + xy(2)*cos(theta)
end function rotate_2d_DP
!
pure function rotate_2d_QP(xy,theta) result(return_value)
  real(qp),               intent(in) :: theta
  real(qp), dimension(:), intent(in) :: xy
  real(qp), dimension(1:size(xy)) :: return_value
  return_value(1) = xy(1)*cos(theta) - xy(2)*sin(theta)
  return_value(2) = xy(1)*sin(theta) + xy(2)*cos(theta)
end function rotate_2d_QP
