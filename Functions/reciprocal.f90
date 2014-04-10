pure function reciprocal(input_value) result(return_value)
  real(vbp), intent(in) :: input_value
  real(vbp) :: return_value
  real(vbp) :: input_nonzero,input_ratio
  input_nonzero = max(abs(input_value),rndoff)
  input_ratio = input_value / input_nonzero
  return_value = one / input_nonzero*input_ratio
end function reciprocal
