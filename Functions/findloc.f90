pure function findloc_i4(array,value,mask,back) result(return_value)
  !.. Formal Arguments ..
  integer(i4), dimension(:), intent(in) :: array
  integer(i4),               intent(in) :: value
  logical,     dimension(:), intent(in) :: mask
  !.. Optional Arguments ..
  logical, optional,         intent(in) :: back
  !.. Function Result ..
  integer :: return_value
  !.. Local Scalars ..
  logical :: reverse
continue
  return_value = 0
  if (count(array==value .and. mask) == 0) return
  reverse = .false.
  if (present(back)) reverse = back
  if (reverse) then
    return_value = maxloc(merge(1,0,array(size(array):1:-1)==value .and. &
                                     mask(size(array):1:-1)),dim=1)
  else
    return_value = maxloc(merge(1,0,array==value.and.mask),dim=1)
  end if
end function findloc_i4
!
!###############################################################################
!
pure function findloc_i8(array,value,mask,back) result(return_value)
  !.. Formal Arguments ..
  integer(i8), dimension(:), intent(in) :: array
  integer(i8),               intent(in) :: value
  logical,     dimension(:), intent(in) :: mask
  !.. Optional Arguments ..
  logical, optional,         intent(in) :: back
  !.. Function Result ..
  integer :: return_value
  !.. Local Scalars ..
  logical :: reverse
continue
  return_value = 0
  if (count(array==value .and. mask) == 0) return
  reverse = .false.
  if (present(back)) reverse = back
  if (reverse) then
    return_value = maxloc(merge(1,0,array(size(array):1:-1)==value .and. &
                                     mask(size(array):1:-1)),dim=1)
  else
    return_value = maxloc(merge(1,0,array==value.and.mask),dim=1)
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
  !.. Local Scalars ..
  logical :: reverse
continue
  !
  return_value = 0
  !
  if (present(mask)) then
    if (count(array==value .and. mask) == 0) return
  else
    if (count(array==value) == 0) return
  end if
  !
  reverse = .false.
  if (present(back)) reverse = back
  !
  if (present(mask)) then
    if (reverse) then
      return_value = maxloc(merge(1,0,array(size(array):1:-1)==value .and. &
                                       mask(size(array):1:-1)),dim=1)
    else
      return_value = maxloc(merge(1,0,array==value.and.mask),dim=1)
    end if
  else
    if (reverse) then
      return_value = maxloc(merge(1,0,array(size(array):1:-1)==value),dim=1)
    else
      return_value = maxloc(merge(1,0,array==value),dim=1)
    end if
  end if
  !
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
  !.. Local Scalars ..
  logical :: reverse
continue
  !
  return_value = 0
  !
  if (present(mask)) then
    if (count(array==value .and. mask) == 0) return
  else
    if (count(array==value) == 0) return
  end if
  !
  reverse = .false.
  if (present(back)) reverse = back
  !
  if (present(mask)) then
    if (reverse) then
      return_value = maxloc(merge(1,0,array(size(array):1:-1)==value .and. &
                                       mask(size(array):1:-1)),dim=1)
    else
      return_value = maxloc(merge(1,0,array==value.and.mask),dim=1)
    end if
  else
    if (reverse) then
      return_value = maxloc(merge(1,0,array(size(array):1:-1)==value),dim=1)
    else
      return_value = maxloc(merge(1,0,array==value),dim=1)
    end if
  end if
  !
end function findloc_r8
!
!###############################################################################
!
pure recursive function recursive_findloc_i4(array,value) result(return_value)
  !.. Formal Arguments ..
  integer(i4), dimension(:), intent(in) :: array
  integer(i4),               intent(in) :: value
  !.. Function Result ..
  integer :: return_value
  !.. Local Scalars ..
  integer :: i,imid
  !
continue
  !
  return_value = 0
  !
  imid = max( 1 , size(array)/2 )
  !
  if (value == array(imid)) then
    return_value = imid
  else if (size(array) == 1 .and. value /= array(imid)) then
    return_value = -huge(0)
  else if (value < array(imid)) then
    return_value = recursive_findloc_i4( array(:imid-1) , value )
  else
    return_value = imid + recursive_findloc_i4( array(imid+1:) , value )
  end if
  !
end function recursive_findloc_i4
!
!##############################################################################
!
pure recursive function recursive_findloc_i8(array,value) result(return_value)
  !.. Formal Arguments ..
  integer(i8), dimension(:), intent(in) :: array
  integer(i8),               intent(in) :: value
  !.. Function Result ..
  integer :: return_value
  !.. Local Scalars ..
  integer :: i,imid
  !
continue
  !
  return_value = 0
  !
  imid = max( 1 , size(array)/2 )
  !
  if (value == array(imid)) then
    return_value = imid
  else if (size(array) == 1 .and. value /= array(imid)) then
    return_value = -huge(0)
  else if (value < array(imid)) then
    return_value = recursive_findloc_i8( array(:imid-1) , value )
  else
    return_value = imid + recursive_findloc_i8( array(imid+1:) , value )
  end if
  !
end function recursive_findloc_i8
!
!##############################################################################
!
pure recursive function recursive_findloc_jeff_patch(array,value) &
                                               result(return_value)
  !
  ! Seths recursive version of find_loc
  ! with mods to prevent race condition
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: array
  integer,               intent(in) :: value
  !.. Function Result ..
  integer :: return_value
  !.. Local Scalars ..
  integer :: i,imid,asiz,amid
  !
continue
  !
  asiz = size(array)
  !
  if (asiz == 0) then
    return_value = 0 ; return
  endif
  !
  imid = max( 1 , asiz/2 )
  !
  amid = array(imid)
  !
  return_value = -1
  !
  !write(iout,*) " size(array), imid, array(imid), value =", &
  !                asiz, imid, array(imid), value
  if (value == amid) then
    return_value = imid
  else if (asiz == 1 .and. value /= amid) then
    return_value = 0
  else if (value < amid) then
    return_value = recursive_findloc_jeff_patch( array(:imid-1) , value )
  else
    return_value = recursive_findloc_jeff_patch( array(imid+1:) , value )
    if (return_value > 0) then
      return_value = imid  + return_value
    end if
  end if
  !write(iout,*) " imid, return_value = ", imid, return_value
  !
end function recursive_findloc_jeff_patch
!
!##############################################################################
!
