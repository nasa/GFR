pure function uppercase(string) result(return_value)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: string
  !
  !.. Function Return Value ..
  character(len=len(string)) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,n
  !
continue
  !
  ! This uses the ASCII character sequence where
  ! A-Z are numbers 65-90 and a-z are 97-122
  !
  return_value = string
  !
  do i = 97,122
    forall ( n=1:len_trim(string) , string(n:n) == achar(i) ) &
                              return_value(n:n) = achar(i-32)
  end do
  !
end function uppercase
!
!###############################################################################
!
