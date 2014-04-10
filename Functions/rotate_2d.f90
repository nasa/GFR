pure function rotate_2d(xy,theta) result(return_value)
  real(vbp),               intent(in) :: theta
  real(vbp), dimension(:), intent(in) :: xy
  real(vbp), dimension(1:size(xy)) :: return_value
  return_value(1) = xy(1)*cos(theta) - xy(2)*sin(theta)
  return_value(2) = xy(1)*sin(theta) + xy(2)*cos(theta)
end function rotate_2d
