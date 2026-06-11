!
!###############################################################################
!###############################################################################
! 32-bit integer versions
!###############################################################################
!###############################################################################
!
pure function all_are_equal_i4(array) result(return_value)
  use, intrinsic :: iso_fortran_env, only : lp => int32
  integer(lp), dimension(:), intent(in) :: array
  logical(lk) :: return_value
continue
  if (size(array) == 0) then
    return_value = true
  else
    return_value = (minval(array) == maxval(array))
  end if
end function all_are_equal_i4
!
pure function all_not_equal_i4(array) result(return_value)
  use, intrinsic :: iso_fortran_env, only : lp => int32
  integer(lp), dimension(:), intent(in) :: array
  logical(lk) :: return_value
continue
  if (size(array) == 0) then
    return_value = fals
  else
    return_value = (minval(array) /= maxval(array))
  end if
end function all_not_equal_i4
!
!###############################################################################
!###############################################################################
! 64-bit integer versions
!###############################################################################
!###############################################################################
!
pure function all_are_equal_i8(array) result(return_value)
  use, intrinsic :: iso_fortran_env, only : lp => int64
  integer(lp), dimension(:), intent(in) :: array
  logical(lk) :: return_value
continue
  if (size(array) == 0) then
    return_value = true
  else
    return_value = (minval(array) == maxval(array))
  end if
end function all_are_equal_i8
!
pure function all_not_equal_i8(array) result(return_value)
  use, intrinsic :: iso_fortran_env, only : lp => int64
  integer(lp), dimension(:), intent(in) :: array
  logical(lk) :: return_value
continue
  if (size(array) == 0) then
    return_value = fals
  else
    return_value = (minval(array) /= maxval(array))
  end if
end function all_not_equal_i8
!
!###############################################################################
!###############################################################################
! 32-bit real versions
!###############################################################################
!###############################################################################
!
pure function all_are_equal_r4(array,acceptable_error) result(return_value)
  use, intrinsic :: iso_fortran_env, only : lp => real32
  real(lp), dimension(:), intent(in) :: array
  real(lp), optional,     intent(in) :: acceptable_error
  logical(lk) :: return_value
  real(lp) :: error_test
continue
  error_test = epsilon(1.0_lp) ! 1.192e-7 for real32
  if (present(acceptable_error)) then
    error_test = acceptable_error
  end if
  if (size(array) == 0) then
    return_value = true
  else
    return_value = ( abs((maxval(array)-minval(array))) <= error_test )
  end if
end function all_are_equal_r4
!
pure function all_not_equal_r4(array,acceptable_error) result(return_value)
  use, intrinsic :: iso_fortran_env, only : lp => real32
  real(lp), dimension(:), intent(in) :: array
  real(lp), optional,     intent(in) :: acceptable_error
  logical(lk) :: return_value
  real(lp) :: error_test
continue
  error_test = epsilon(1.0_lp) ! 1.192e-7 for real32
  if (present(acceptable_error)) then
    error_test = acceptable_error
  end if
  if (size(array) == 0) then
    return_value = fals
  else
    return_value = ( abs((maxval(array)-minval(array))) > error_test )
  end if
end function all_not_equal_r4
!
!###############################################################################
!###############################################################################
! 64-bit real versions
!###############################################################################
!###############################################################################
!
pure function all_are_equal_r8(array,acceptable_error) result(return_value)
  use, intrinsic :: iso_fortran_env, only : lp => real64
  real(lp), dimension(:), intent(in) :: array
  real(lp), optional,     intent(in) :: acceptable_error
  logical(lk) :: return_value
  real(lp) :: error_test
continue
  error_test = eps12
  if (present(acceptable_error)) then
    error_test = acceptable_error
  end if
  if (size(array) == 0) then
    return_value = true
  else
    return_value = ( abs((maxval(array)-minval(array))) <= error_test )
  end if
end function all_are_equal_r8
!
pure function all_not_equal_r8(array,acceptable_error) result(return_value)
  use, intrinsic :: iso_fortran_env, only : lp => real64
  real(lp), dimension(:), intent(in) :: array
  real(lp), optional,     intent(in) :: acceptable_error
  logical(lk) :: return_value
  real(lp) :: error_test
continue
  error_test = eps12
  if (present(acceptable_error)) then
    error_test = acceptable_error
  end if
  if (size(array) == 0) then
    return_value = fals
  else
    return_value = ( abs((maxval(array)-minval(array))) > error_test )
  end if
end function all_not_equal_r8
!
!###############################################################################
!###############################################################################
! 64-bit real versions
!###############################################################################
!###############################################################################
!
pure function all_are_equal_r16(array,acceptable_error) result(return_value)
  integer, parameter :: lp = QP
  real(lp), dimension(:), intent(in) :: array
  real(lp), optional,     intent(in) :: acceptable_error
  logical(lk) :: return_value
  real(lp) :: error_test
continue
  error_test = eps12
  if (present(acceptable_error)) then
    error_test = acceptable_error
  end if
  if (size(array) == 0) then
    return_value = true
  else
    return_value = ( abs((maxval(array)-minval(array))) <= error_test )
  end if
end function all_are_equal_r16
!
pure function all_not_equal_r16(array,acceptable_error) result(return_value)
  integer, parameter :: lp = QP
  real(lp), dimension(:), intent(in) :: array
  real(lp), optional,     intent(in) :: acceptable_error
  logical(lk) :: return_value
  real(lp) :: error_test
continue
  error_test = eps12
  if (present(acceptable_error)) then
    error_test = acceptable_error
  end if
  if (size(array) == 0) then
    return_value = fals
  else
    return_value = ( abs((maxval(array)-minval(array))) > error_test )
  end if
end function all_not_equal_r16
