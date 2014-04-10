elemental function deltafun(logic) result(return_value)
  logical, intent(in) :: logic
  real(vbp) :: return_value
  return_value = merge(one,zero,logic)
end function deltafun
