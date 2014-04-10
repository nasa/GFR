pure function np2n(np) result(return_value)
  integer, intent(in) :: np
  integer :: return_value
  return_value = (int(sqrt(real(9-8*(1-np),vbp)))-3)/2
end function np2n
