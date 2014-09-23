pure function findloc_i4(array,value,mask,back) result(return_value)
  !.. Formal Arguments ..
  integer(i4), dimension(:), intent(in) :: array
  integer(i4),               intent(in) :: value
  !.. Optional Arguments ..
  logical, optional, dimension(:), intent(in) :: mask
  logical, optional,               intent(in) :: back
  !.. Function Result ..
  integer :: return_value
continue
  return_value = 0
  if (present(mask)) then
    if (count(array==value .and. mask) == 0) return
  else
    if (count(array==value) == 0) return
  end if
  if (present(mask) .and. present(back)) then
    return_value = maxloc(merge(1,0,array==value.and.mask),dim=1)
  else if (present(mask)) then
    return_value = maxloc(merge(1,0,array(size(array):1:-1)==value .and. &
                                     mask(size(array):1:-1)),dim=1)
  else if (present(back)) then
    return_value = maxloc(merge(1,0,array(size(array):1:-1)==value),dim=1)
  else
    return_value = maxloc(merge(1,0,array==value),dim=1)
  end if
end function findloc_i4
!
!###############################################################################
!
pure function findloc_i8(array,value,mask,back) result(return_value)
  !.. Formal Arguments ..
  integer(i8), dimension(:), intent(in) :: array
  integer(i8),               intent(in) :: value
  !.. Optional Arguments ..
  logical, optional, dimension(:), intent(in) :: mask
  logical, optional,               intent(in) :: back
  !.. Function Result ..
  integer :: return_value
continue
  return_value = 0
  if (present(mask)) then
    if (count(array==value .and. mask) == 0) return
  else
    if (count(array==value) == 0) return
  end if
  if (present(mask) .and. present(back)) then
    return_value = maxloc(merge(1,0,array==value.and.mask),dim=1)
  else if (present(mask)) then
    return_value = maxloc(merge(1,0,array(size(array):1:-1)==value .and. &
                                     mask(size(array):1:-1)),dim=1)
  else if (present(back)) then
    return_value = maxloc(merge(1,0,array(size(array):1:-1)==value),dim=1)
  else
    return_value = maxloc(merge(1,0,array==value),dim=1)
  end if
end function findloc_i8
!
!###############################################################################
!
pure function findloc_r4(array,value,mask,back) result(return_value)
  !.. Formal Arguments ..
  real(r4), dimension(:), intent(in) :: array
  real(r4),               intent(in) :: value
  !.. Optional Arguments ..
  logical, optional, dimension(:), intent(in) :: mask
  logical, optional,               intent(in) :: back
  !.. Function Result ..
  integer :: return_value
continue
  return_value = 0
  if (present(mask)) then
    if (count(array==value .and. mask) == 0) return
  else
    if (count(array==value) == 0) return
  end if
  if (present(mask) .and. present(back)) then
    return_value = maxloc(merge(1,0,array==value.and.mask),dim=1)
  else if (present(mask)) then
    return_value = maxloc(merge(1,0,array(size(array):1:-1)==value .and. &
                                     mask(size(array):1:-1)),dim=1)
  else if (present(back)) then
    return_value = maxloc(merge(1,0,array(size(array):1:-1)==value),dim=1)
  else
    return_value = maxloc(merge(1,0,array==value),dim=1)
  end if
end function findloc_r4
!
!###############################################################################
!
pure function findloc_r8(array,value,mask,back) result(return_value)
  !.. Formal Arguments ..
  real(r8), dimension(:), intent(in) :: array
  real(r8),               intent(in) :: value
  !.. Optional Arguments ..
  logical, optional, dimension(:), intent(in) :: mask
  logical, optional,               intent(in) :: back
  !.. Function Result ..
  integer :: return_value
continue
  return_value = 0
  if (present(mask)) then
    if (count(array==value .and. mask) == 0) return
  else
    if (count(array==value) == 0) return
  end if
  if (present(mask) .and. present(back)) then
    return_value = maxloc(merge(1,0,array==value.and.mask),dim=1)
  else if (present(mask)) then
    return_value = maxloc(merge(1,0,array(size(array):1:-1)==value .and. &
                                     mask(size(array):1:-1)),dim=1)
  else if (present(back)) then
    return_value = maxloc(merge(1,0,array(size(array):1:-1)==value),dim=1)
  else
    return_value = maxloc(merge(1,0,array==value),dim=1)
  end if
end function findloc_r8
!
!###############################################################################
!
