pure function bc_priority(boundary_conditions) result(return_value)
  integer, dimension(:), intent(in) :: boundary_conditions
  integer :: return_value
  integer :: n
continue
  do n = 1,size(bc_heirarchy)
#ifdef USE_INTRINSIC_FINDLOC
    return_value = findloc( boundary_conditions , bc_heirarchy(n) , dim=1 )
#else
    return_value = findloc( boundary_conditions , bc_heirarchy(n) )
#endif
    if (return_value /= 0) exit
  end do
end function bc_priority
