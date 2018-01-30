pure function reformat_string(input_strng,line_length,leading_blanks) &
                       result(return_value)
!function reformat_string(input_strng,line_length,leading_blanks) &
!                  result(return_value)
  !
  character(len=*), intent(in) :: input_strng
  integer, optional, intent(in) :: line_length
  integer, optional, intent(in) :: leading_blanks
  !
  character(len=len(input_strng)) :: return_value
  !
 !integer :: j
  integer :: n,last_space,max_strng_len
  integer :: llen,nblanks
  character(len=len(input_strng)) :: strng
  character(len=100) :: frmt
  !
continue
  !
  return_value = ""
  !
  llen = 62
  nblanks = 17
  strng = trim(adjustl(input_strng))
  !
  if (present(line_length)) llen = min(line_length,100)
  if (present(leading_blanks)) nblanks = leading_blanks
  !
  write (frmt,1) nblanks
  !
  ! Change all tab characters to a single space
  !
  forall (n=1:len_trim(strng),strng(n:n)==achar(9)) strng(n:n) = " "
  !
 !j = 0
  interval_loop: do
    !
   !j = j + 1
   !write (*,2) j
    !
    ! Reset max_strng_len for what is left of strng
    !
    max_strng_len = len_trim(strng)
   !write (*,3) "max_strng_len",max_strng_len
    !
    ! Maximum value between first space character and line length, limited
    ! to the length of the what is left of strng (this is in case the
    ! interval 1:llen does not contain a space character)
    !
    n = min( max_strng_len , max(llen,scan(strng," ")) )
   !write (*,3) "n",n
    !
    ! Location of last space character on interval 1:n
    !
    last_space = scan(strng(1:n)," ",back=true)
   !write (*,3) "last_space",last_space
    !
    ! Make sure that we keep the same value for 'n' if the interval 1:n
    ! contains no spaces or 'n' is the same as the length of what is left of
    ! strng.
    !
    n = merge( n , last_space , last_space==0 .or. n==max_strng_len )
   !write (*,3) "n using last_space",n,"len(substrng)",len_trim(strng(1:n))
   !write (*,4) trim(strng(1:n))
    !
    if (len_trim(return_value) == 0)then
      write (return_value,fmt=trim(frmt)) trim(strng(1:n))
    else
      write (return_value,fmt=trim(frmt)) trim(return_value)//new_line('a'), &
                                          trim(strng(1:n))
    end if
    !
    if (len_trim(strng) == n) exit interval_loop
    !
    strng = adjustl(strng(n+1:max_strng_len))
    !
  end do interval_loop
  !
  ! Format Statements
  !
  1 format ('(a,:,',i0,'x,a)')
 !2 format (/,' ** Loop ',i0,' **',/)
 !3 format (2x,100(a," = ",i0,:,"; "))
 !4 format ('"',a,'"')
  !
end function reformat_string
