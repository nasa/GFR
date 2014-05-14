elemental function chop(a) result(return_value)
  real(vbp), intent(in) :: a
  real(vbp) :: return_value
  return_value = a * merge( one , zero , abs(a)>two*epsilon(a) )
end function chop
