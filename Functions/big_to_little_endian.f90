elemental function big_to_little_endian_i2(n) result(return_value)
  !
  ! Function to convert a big_endian i_2 (32 bit) integer
  !               to a little_endian i_2 (32 bit) integer
  !
  use iso_fortran_env, only : i2 => int16
  !
  implicit none
  !
  integer(i2), intent(in) :: n
  !
  integer(i2) :: return_value
  !
  integer(i4) :: i
  !
continue
  !
  return_value = 0_i2
  !
  do i = 0_i4,8_i4,8_i4
    call mvbits(n,8_i4-i,8_i4,return_value,i)
  end do
  !
end function big_to_little_endian_i2
!
!###############################################################################
!
elemental function big_to_little_endian_i4(n) result(return_value)
  !
  ! Function to convert a big_endian i_4 (32 bit) integer
  !               to a little_endian i_4 (32 bit) integer
  !
  implicit none
  !
  integer(i4), intent(in) :: n
  !
  integer(i4) :: return_value
  !
  integer(i4) :: i
  !
continue
  !
  return_value = 0_i4
  !
  do i = 0_i4,24_i4,8_i4
    call mvbits(n,24_i4-i,8_i4,return_value,i)
  end do
  !
end function big_to_little_endian_i4
!
!###############################################################################
!
elemental function big_to_little_endian_i8(m) result(return_value)
  !
  ! Function to convert a big_endian i_8 (64 bit) integer
  !               to a little_endian i_8 (64 bit) integer
  !
  implicit none
  !
  integer(i8), intent(in) :: m
  !
  integer(i8) :: return_value
  !
  integer(i4) :: i
  !
continue
  !
  return_value = 0_i8
  !
  do i = 0_i4,56_i4,8_i4
    call mvbits(m,56_i4-i,8_i4,return_value,i)
  end do
  !
end function big_to_little_endian_i8
!
!###############################################################################
!
elemental function big_to_little_endian_r4(s) result(return_value)
  !
  ! Function to convert a big_endian r_4 (32 bit) real
  !               to a little_endian r_4 (32 bit) real
  !
  implicit none
  !
  real(r4), intent(in) :: s
  !
  real(r4) :: return_value
  !
  integer(i4) :: i
  integer(i4) :: int_s
  integer(i4) :: int_s_br
  !
continue
  !
  return_value = 0.0_r4
  !
  int_s = 0_i4
  int_s_br = 0_i4
  !
  int_s = transfer(s,int_s)
  do i = 0_i4,24_i4,8_i4
    call mvbits(int_s,24_i4-i,8_i4,int_s_br,i)
  end do
  return_value = transfer(int_s_br,s)
  !
end function big_to_little_endian_r4
!
!###############################################################################
!
elemental function big_to_little_endian_r8(d) result(return_value)
  !
  ! Function to convert a big_endian r_8 (64 bit) real
  !               to a little_endian r_8 (64 bit) real
  !
  implicit none
  !
  real(r8), intent(in) :: d
  !
  real(r8) :: return_value
  !
  integer(i4) :: i
  integer(i8) :: int_d
  integer(i8) :: int_d_br
  !
continue
  !
  return_value = 0.0_r8
  !
  int_d = 0_i8
  int_d_br = 0_i8
  !
  int_d = transfer(d,int_d)
  do i = 0_i4,56_i4,8_i4
    call mvbits(int_d,56_i4-i,8_i4,int_d_br,i)
  end do
  return_value = transfer(int_d_br,d)
  !
end function big_to_little_endian_r8
!
!###############################################################################
!
