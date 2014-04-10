pure function last(array) result(return_value)
  integer, dimension(:), intent(in) :: array
  integer :: return_value
  return_value = array(size(array))
end function last
