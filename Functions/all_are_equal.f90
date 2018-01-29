pure function all_are_equal(array) result(return_value)
  integer, dimension(:), intent(in) :: array
  logical(lk) :: return_value
continue
  if (size(array) == 0) then
    return_value = true
  else
    return_value = (minval(array) == maxval(array))
  end if
end function all_are_equal
!
pure function all_not_equal(array) result(return_value)
  integer, dimension(:), intent(in) :: array
  logical(lk) :: return_value
continue
  if (size(array) == 0) then
    return_value = fals
  else
    return_value = (minval(array) /= maxval(array))
  end if
end function all_not_equal
