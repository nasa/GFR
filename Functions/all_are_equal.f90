pure function all_are_equal(array) result(return_value)
  integer, dimension(:), intent(in) :: array
  logical :: return_value
continue
  return_value = (minval(array) == maxval(array))
end function all_are_equal
