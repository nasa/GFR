function free_unit() result(return_value)
  integer :: return_value
  integer :: iunit
  logical(ldk) :: is_open
  integer, parameter :: min_unit = 100
  integer, parameter :: max_unit = 500
  do iunit = min_unit,max_unit
    return_value = iunit
    inquire(unit=iunit,opened=is_open)
    if (.not.is_open) exit
  end do
end function free_unit
