!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
! INIT_ARRAY FUNCTIONS WITH SIZE1,SIZE2,SIZE3 INPUTS AS 32-BIT INTEGERS
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
!
!###############################################################################
!###############################################################################
!###############################################################################
! 32-BIT SIZE Functions with val input as 8-bit integer
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_i1_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int8
  use iso_fortran_env, only : ip => int32
  !
  integer(ip), intent(in) :: size1
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_i1_rank1
!
!###############################################################################
!
pure function init_array_i1_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int8
  use iso_fortran_env, only : ip => int32
  !
  integer(ip), intent(in) :: size1,size2
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_i1_rank2
!
!###############################################################################
!
pure function init_array_i1_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int8
  use iso_fortran_env, only : ip => int32
  !
  integer(ip), intent(in) :: size1,size2,size3
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_i1_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! 32-BIT SIZE Functions with val input as 16-bit integer
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_i2_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int16
  use iso_fortran_env, only : ip => int32
  !
  integer(ip), intent(in) :: size1
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_i2_rank1
!
!###############################################################################
!
pure function init_array_i2_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int16
  use iso_fortran_env, only : ip => int32
  !
  integer(ip), intent(in) :: size1,size2
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_i2_rank2
!
!###############################################################################
!
pure function init_array_i2_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int16
  use iso_fortran_env, only : ip => int32
  !
  integer(ip), intent(in) :: size1,size2,size3
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_i2_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! 32-BIT SIZE Functions with val input as 32-bit integer
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_i4_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int32
  use iso_fortran_env, only : ip => int32
  !
  integer(ip), intent(in) :: size1
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_i4_rank1
!
!###############################################################################
!
pure function init_array_i4_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int32
  use iso_fortran_env, only : ip => int32
  !
  integer(ip), intent(in) :: size1,size2
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_i4_rank2
!
!###############################################################################
!
pure function init_array_i4_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int32
  use iso_fortran_env, only : ip => int32
  !
  integer(ip), intent(in) :: size1,size2,size3
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_i4_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! 32-BIT SIZE Functions with val input as 64-bit integer
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_i8_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int64
  use iso_fortran_env, only : ip => int32
  !
  integer(ip), intent(in) :: size1
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_i8_rank1
!
!###############################################################################
!
pure function init_array_i8_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int64
  use iso_fortran_env, only : ip => int32
  !
  integer(ip), intent(in) :: size1,size2
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_i8_rank2
!
!###############################################################################
!
pure function init_array_i8_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int64
  use iso_fortran_env, only : ip => int32
  !
  integer(ip), intent(in) :: size1,size2,size3
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_i8_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! 32-BIT SIZE Functions with val input as 32-bit real
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_r4_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real32
  use iso_fortran_env, only : ip => int32
  !
  integer(ip),  intent(in) :: size1
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_r4_rank1
!
!###############################################################################
!
pure function init_array_r4_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real32
  use iso_fortran_env, only : ip => int32
  !
  integer(ip),  intent(in) :: size1,size2
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_r4_rank2
!
!###############################################################################
!
pure function init_array_r4_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real32
  use iso_fortran_env, only : ip => int32
  !
  integer(ip),  intent(in) :: size1,size2,size3
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_r4_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! 32-BIT SIZE Functions with val input as 64-bit real
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_r8_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real64
  use iso_fortran_env, only : ip => int32
  !
  integer(ip),  intent(in) :: size1
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_r8_rank1
!
!###############################################################################
!
pure function init_array_r8_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real64
  use iso_fortran_env, only : ip => int32
  !
  integer(ip),  intent(in) :: size1,size2
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_r8_rank2
!
!###############################################################################
!
pure function init_array_r8_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real64
  use iso_fortran_env, only : ip => int32
  !
  integer(ip),  intent(in) :: size1,size2,size3
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_r8_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! 32-BIT SIZE Functions with val input as 128-bit real
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_r16_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real128
  use iso_fortran_env, only : ip => int32
  !
  integer(ip),  intent(in) :: size1
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_r16_rank1
!
!###############################################################################
!
pure function init_array_r16_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real128
  use iso_fortran_env, only : ip => int32
  !
  integer(ip),  intent(in) :: size1,size2
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_r16_rank2
!
!###############################################################################
!
pure function init_array_r16_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real128
  use iso_fortran_env, only : ip => int32
  !
  integer(ip),  intent(in) :: size1,size2,size3
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_r16_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! 32-BIT SIZE Functions with val input as logical kind of least precision
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_lk_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : ip => int32
  !
  integer, parameter :: lp = lk
  !
  integer(ip), intent(in) :: size1
  logical(lp), intent(in) :: val
  !
  logical(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_lk_rank1
!
!###############################################################################
!
pure function init_array_lk_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : ip => int32
  !
  integer, parameter :: lp = lk
  !
  integer(ip), intent(in) :: size1,size2
  logical(lp), intent(in) :: val
  !
  logical(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_lk_rank2
!
!###############################################################################
!
pure function init_array_lk_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : ip => int32
  !
  integer, parameter :: lp = lk
  !
  integer(ip), intent(in) :: size1,size2,size3
  logical(lp), intent(in) :: val
  !
  logical(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_lk_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! 32-BIT SIZE Functions with val input as compiler-default logical kind
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_l4_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : ip => int32
  !
  integer, parameter :: lp = ldk
  !
  integer(ip), intent(in) :: size1
  logical(lp), intent(in) :: val
  !
  logical(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_l4_rank1
!
!###############################################################################
!
pure function init_array_l4_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : ip => int32
  !
  integer, parameter :: lp = ldk
  !
  integer(ip), intent(in) :: size1,size2
  logical(lp), intent(in) :: val
  !
  logical(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_l4_rank2
!
!###############################################################################
!
pure function init_array_l4_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : ip => int32
  !
  integer, parameter :: lp = ldk
  !
  integer(ip), intent(in) :: size1,size2,size3
  logical(lp), intent(in) :: val
  !
  logical(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_l4_rank3
!
!
!
!
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
! INIT_ARRAY FUNCTIONS WITH SIZE1,SIZE2,SIZE3 INPUTS AS 64-BIT INTEGERS
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
!
!
!
!
!###############################################################################
!###############################################################################
!###############################################################################
! 64-BIT SIZE Functions with val input as 8-bit integer
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_s8_i1_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int8
  use iso_fortran_env, only : ip => int64
  !
  integer(ip), intent(in) :: size1
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_s8_i1_rank1
!
!###############################################################################
!
pure function init_array_s8_i1_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int8
  use iso_fortran_env, only : ip => int64
  !
  integer(ip), intent(in) :: size1,size2
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_s8_i1_rank2
!
!###############################################################################
!
pure function init_array_s8_i1_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int8
  use iso_fortran_env, only : ip => int64
  !
  integer(ip), intent(in) :: size1,size2,size3
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_s8_i1_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! 64-BIT SIZE Functions with val input as 16-bit integer
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_s8_i2_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int16
  use iso_fortran_env, only : ip => int64
  !
  integer(ip), intent(in) :: size1
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_s8_i2_rank1
!
!###############################################################################
!
pure function init_array_s8_i2_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int16
  use iso_fortran_env, only : ip => int64
  !
  integer(ip), intent(in) :: size1,size2
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_s8_i2_rank2
!
!###############################################################################
!
pure function init_array_s8_i2_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int16
  use iso_fortran_env, only : ip => int64
  !
  integer(ip), intent(in) :: size1,size2,size3
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_s8_i2_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! 64-BIT SIZE Functions with val input as 32-bit integer
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_s8_i4_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int32
  use iso_fortran_env, only : ip => int64
  !
  integer(ip), intent(in) :: size1
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_s8_i4_rank1
!
!###############################################################################
!
pure function init_array_s8_i4_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int32
  use iso_fortran_env, only : ip => int64
  !
  integer(ip), intent(in) :: size1,size2
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_s8_i4_rank2
!
!###############################################################################
!
pure function init_array_s8_i4_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int32
  use iso_fortran_env, only : ip => int64
  !
  integer(ip), intent(in) :: size1,size2,size3
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_s8_i4_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! 64-BIT SIZE Functions with val input as 64-bit integer
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_s8_i8_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int64
  use iso_fortran_env, only : ip => int64
  !
  integer(ip), intent(in) :: size1
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_s8_i8_rank1
!
!###############################################################################
!
pure function init_array_s8_i8_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int64
  use iso_fortran_env, only : ip => int64
  !
  integer(ip), intent(in) :: size1,size2
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_s8_i8_rank2
!
!###############################################################################
!
pure function init_array_s8_i8_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : lp => int64
  use iso_fortran_env, only : ip => int64
  !
  integer(ip), intent(in) :: size1,size2,size3
  integer(lp), intent(in) :: val
  !
  integer(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_s8_i8_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! 64-BIT SIZE Functions with val input as 32-bit real
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_s8_r4_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real32
  use iso_fortran_env, only : ip => int64
  !
  integer(ip),  intent(in) :: size1
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_s8_r4_rank1
!
!###############################################################################
!
pure function init_array_s8_r4_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real32
  use iso_fortran_env, only : ip => int64
  !
  integer(ip),  intent(in) :: size1,size2
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_s8_r4_rank2
!
!###############################################################################
!
pure function init_array_s8_r4_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real32
  use iso_fortran_env, only : ip => int64
  !
  integer(ip),  intent(in) :: size1,size2,size3
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_s8_r4_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! 64-BIT SIZE Functions with val input as 64-bit real
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_s8_r8_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real64
  use iso_fortran_env, only : ip => int64
  !
  integer(ip),  intent(in) :: size1
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_s8_r8_rank1
!
!###############################################################################
!
pure function init_array_s8_r8_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real64
  use iso_fortran_env, only : ip => int64
  !
  integer(ip),  intent(in) :: size1,size2
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_s8_r8_rank2
!
!###############################################################################
!
pure function init_array_s8_r8_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real64
  use iso_fortran_env, only : ip => int64
  !
  integer(ip),  intent(in) :: size1,size2,size3
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_s8_r8_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! 64-BIT SIZE Functions with val input as 128-bit real
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_s8_r16_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real128
  use iso_fortran_env, only : ip => int64
  !
  integer(ip),  intent(in) :: size1
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_s8_r16_rank1
!
!###############################################################################
!
pure function init_array_s8_r16_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real128
  use iso_fortran_env, only : ip => int64
  !
  integer(ip),  intent(in) :: size1,size2
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_s8_r16_rank2
!
!###############################################################################
!
pure function init_array_s8_r16_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : lp => real128
  use iso_fortran_env, only : ip => int64
  !
  integer(ip),  intent(in) :: size1,size2,size3
  real(lp), intent(in) :: val
  !
  real(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_s8_r16_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! 64-BIT SIZE Functions with val input as logical kind of least precision
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_s8_lk_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : ip => int64
  !
  integer, parameter :: lp = lk
  !
  integer(ip), intent(in) :: size1
  logical(lp), intent(in) :: val
  !
  logical(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_s8_lk_rank1
!
!###############################################################################
!
pure function init_array_s8_lk_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : ip => int64
  !
  integer, parameter :: lp = lk
  !
  integer(ip), intent(in) :: size1,size2
  logical(lp), intent(in) :: val
  !
  logical(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_s8_lk_rank2
!
!###############################################################################
!
pure function init_array_s8_lk_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : ip => int64
  !
  integer, parameter :: lp = lk
  !
  integer(ip), intent(in) :: size1,size2,size3
  logical(lp), intent(in) :: val
  !
  logical(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_s8_lk_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
! 64-BIT SIZE Functions with val input as compiler-default logical kind
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function init_array_s8_l4_rank1(size1,val) result(return_value)
  !
  use iso_fortran_env, only : ip => int64
  !
  integer, parameter :: lp = ldk
  !
  integer(ip), intent(in) :: size1
  logical(lp), intent(in) :: val
  !
  logical(lp), allocatable :: return_value(:)
  !
continue
  return_value = spread(val,dim=1,ncopies=size1) ! F2003 AUTO-REALLOC
end function init_array_s8_l4_rank1
!
!###############################################################################
!
pure function init_array_s8_l4_rank2(size1,size2,val) result(return_value)
  !
  use iso_fortran_env, only : ip => int64
  !
  integer, parameter :: lp = ldk
  !
  integer(ip), intent(in) :: size1,size2
  logical(lp), intent(in) :: val
  !
  logical(lp), allocatable :: return_value(:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2) , &
 !                        shape=[size1,size2] )   ! F2003 AUTO-REALLOC
end function init_array_s8_l4_rank2
!
!###############################################################################
!
pure function init_array_s8_l4_rank3(size1,size2,size3,val) result(return_value)
  !
  use iso_fortran_env, only : ip => int64
  !
  integer, parameter :: lp = ldk
  !
  integer(ip), intent(in) :: size1,size2,size3
  logical(lp), intent(in) :: val
  !
  logical(lp), allocatable :: return_value(:,:,:)
  !
continue
  ! F2003 AUTO-REALLOC
  return_value = reshape( [val] , shape=[size1,size2,size3] , pad=[val] )
 !return_value = reshape( spread(val,dim=1,ncopies=size1*size2*size3) , &
 !                        shape=[size1,size2,size3] )   ! F2003 AUTO-REALLOC
end function init_array_s8_l4_rank3
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
