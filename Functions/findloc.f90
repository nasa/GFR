pure function findloc_i4(array,value,mask,back) result(return_value)
  !.. Kind Parameters ..
  integer, parameter :: lp = i4
  !.. Formal Arguments ..
  integer(lp), intent(in) :: array(:)
  integer(lp), intent(in) :: value
  logical(lk), intent(in) :: mask(:)
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: back
  !.. Function Result ..
  integer :: return_value
  !.. Local Scalars ..
  logical(lk) :: reverse
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
  !.. Kind Parameters ..
  integer, parameter :: lp = i8
  !.. Formal Arguments ..
  integer(lp), intent(in) :: array(:)
  integer(lp), intent(in) :: value
  logical(lk), intent(in) :: mask(:)
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: back
  !.. Function Result ..
  integer :: return_value
  !.. Local Scalars ..
  logical(lk) :: reverse
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
  !.. Kind Parameters ..
  integer, parameter :: lp = r4
  !.. Formal Arguments ..
  real(lp),    intent(in) :: array(:)
  real(lp),    intent(in) :: value
  logical(lk), intent(in) :: mask(:)
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: back
  !.. Function Result ..
  integer :: return_value
  !.. Local Scalars ..
  logical(lk) :: reverse
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
end function findloc_r4
!
!###############################################################################
!
pure function findloc_r8(array,value,mask,back) result(return_value)
  !.. Kind Parameters ..
  integer, parameter :: lp = r8
  !.. Formal Arguments ..
  real(lp),    intent(in) :: array(:)
  real(lp),    intent(in) :: value
  logical(lk), intent(in) :: mask(:)
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: back
  !.. Function Result ..
  integer :: return_value
  !.. Local Scalars ..
  logical(lk) :: reverse
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
end function findloc_r8
!
!###############################################################################
!
pure recursive function recursive_findloc_i4(array,value) result(return_value)
  !
  ! Recursive version of findloc with mods to prevent race condition
  !
  !.. Kind Parameters ..
  integer, parameter :: lp = i4
  !.. Formal Arguments ..
  integer(lp), intent(in) :: array(:)
  integer(lp), intent(in) :: value
  !.. Function Result ..
  integer :: return_value
  !.. Local Scalars ..
  integer     :: i,imid,asiz
  integer(lp) :: amid
  !
continue
  !
  asiz = size(array)
  !
  if (asiz == 0) then
    return_value = 0 ; return
  end if
  !
  imid = max( 1 , asiz/2 )
  !
  amid = array(imid)
  !
  return_value = -1
  !
  if (value == amid) then
    return_value = imid
  else if (asiz == 1 .and. value /= amid) then
    return_value = 0
  else if (value < amid) then
    return_value = findloc( array(:imid-1) , value )
  else
    return_value = findloc( array(imid+1:) , value )
    if (return_value > 0) then
      return_value = imid + return_value
    end if
  end if
  !
end function recursive_findloc_i4
!
!###############################################################################
!
pure recursive function recursive_findloc_i8(array,value) result(return_value)
  !
  ! Recursive version of findloc with mods to prevent race condition
  !
  !.. Kind Parameters ..
  integer, parameter :: lp = i8
  !.. Formal Arguments ..
  integer(lp), intent(in) :: array(:)
  integer(lp), intent(in) :: value
  !.. Function Result ..
  integer :: return_value
  !.. Local Scalars ..
  integer     :: i,imid,asiz
  integer(lp) :: amid
  !
continue
  !
  asiz = size(array)
  !
  if (asiz == 0) then
    return_value = 0 ; return
  end if
  !
  imid = max( 1 , asiz/2 )
  !
  amid = array(imid)
  !
  return_value = -1
  !
  if (value == amid) then
    return_value = imid
  else if (asiz == 1 .and. value /= amid) then
    return_value = 0
  else if (value < amid) then
    return_value = findloc( array(:imid-1) , value )
  else
    return_value = findloc( array(imid+1:) , value )
    if (return_value > 0) then
      return_value = imid + return_value
    end if
  end if
  !
end function recursive_findloc_i8
!
!###############################################################################
!
pure recursive function recursive_findloc_r4(array,value) result(return_value)
  !
  ! Recursive version of findloc with mods to prevent race condition
  !
  !.. Kind Parameters ..
  integer, parameter :: lp = r4
  !.. Formal Arguments ..
  real(lp), intent(in) :: array(:)
  real(lp), intent(in) :: value
  !.. Function Result ..
  integer :: return_value
  !.. Local Scalars ..
  integer  :: i,imid,asiz
  real(lp) :: amid
  !
continue
  !
  asiz = size(array)
  !
  if (asiz == 0) then
    return_value = 0 ; return
  end if
  !
  imid = max( 1 , asiz/2 )
  !
  amid = array(imid)
  !
  return_value = -1
  !
  if (value == amid) then
    return_value = imid
  else if (asiz == 1 .and. value /= amid) then
    return_value = 0
  else if (value < amid) then
    return_value = findloc( array(:imid-1) , value )
  else
    return_value = findloc( array(imid+1:) , value )
    if (return_value > 0) then
      return_value = imid + return_value
    end if
  end if
  !
end function recursive_findloc_r4
!
!###############################################################################
!
pure recursive function recursive_findloc_r8(array,value) result(return_value)
  !
  ! Recursive version of findloc with mods to prevent race condition
  !
  !.. Kind Parameters ..
  integer, parameter :: lp = r8
  !.. Formal Arguments ..
  real(lp), intent(in) :: array(:)
  real(lp), intent(in) :: value
  !.. Function Result ..
  integer :: return_value
  !.. Local Scalars ..
  integer  :: i,imid,asiz
  real(lp) :: amid
  !
continue
  !
  asiz = size(array)
  !
  if (asiz == 0) then
    return_value = 0 ; return
  end if
  !
  imid = max( 1 , asiz/2 )
  !
  amid = array(imid)
  !
  return_value = -1
  !
  if (value == amid) then
    return_value = imid
  else if (asiz == 1 .and. value /= amid) then
    return_value = 0
  else if (value < amid) then
    return_value = findloc( array(:imid-1) , value )
  else
    return_value = findloc( array(imid+1:) , value )
    if (return_value > 0) then
      return_value = imid + return_value
    end if
  end if
  !
end function recursive_findloc_r8
!
!###############################################################################
!
