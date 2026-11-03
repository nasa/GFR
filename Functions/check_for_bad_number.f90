subroutine check_for_bad_number(line_num,file_name,ival)
  !
#define USE_GATHER
  !.. Use Statements ..
  use module_kind_types, only : ldk,int_mpi
  use module_kind_types, only : mypnum,ncpu,mpierr
  use module_kind_types, only : mpi_inttyp
#ifdef USE_MPI_F08
  use module_kind_types, only : mpi_allgather
  use module_kind_types, only : mpi_allreduce
#endif
  use module_kind_types, only : stop_gfr,stop_mpi
  use module_kind_types, only : src_file_with_bad_number
  use module_kind_types, only : line_num_with_bad_number
  use module_kind_types, only : MPI_MAX
  use module_kind_types, only : MPI_COMM_WORLD
! use module_kind_types, only : MPI_IN_PLACE
#ifdef USE_MPI_F08
  use mpi_f08, only : MPI_IN_PLACE
#else
  use mpi,     only : MPI_IN_PLACE
#endif
#ifndef SPECIAL_FOR_GNU
  use, intrinsic :: ieee_exceptions, only : ieee_usual
  use, intrinsic :: ieee_exceptions, only : ieee_get_flag
#endif
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: line_num
  character(len=*), intent(in) :: file_name
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: ival
  !
#ifndef SPECIAL_FOR_GNU
  !
  !.. Local Scalars ..
  integer :: n,m
  !
  character(len=200) :: overflow_fmt
  character(len=200) :: divide_by_zero_fmt
  character(len=200) :: invalid_fmt
  !
  character(len=:), allocatable :: overflow_error
  character(len=:), allocatable :: divide_by_zero_error
  character(len=:), allocatable :: invalid_error
  character(len=:), allocatable :: output_message
  !
  !.. Local Arrays ..
  integer      :: bad_number_found(1:size(ieee_usual))
  integer      :: bad_numbers_found(1:size(ieee_usual),1:ncpu)
  logical(ldk) :: signaling_exceptions(1:size(ieee_usual))
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: signaling_cpus(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_for_bad_number"
  !
continue
#ifndef DISABLE_CHECK_FOR_BAD_NUMBER
  !
  overflow_error = repeat(" ",1000000)       ! F2003 AUTO-REALLOC
  divide_by_zero_error = repeat(" ",1000000) ! F2003 AUTO-REALLOC
  invalid_error = repeat(" ",1000000)        ! F2003 AUTO-REALLOC
  output_message = repeat(" ",4000000)       ! F2003 AUTO-REALLOC
  !
  ! First check to make sure that no floating-point exceptions are signaling
  !
  call ieee_get_flag(ieee_usual,signaling_exceptions)
  !
#ifdef USE_GATHER
  bad_numbers_found(:,:) = 0
  bad_numbers_found(:,mypnum+1) = merge(1,0,signaling_exceptions)
#else
  bad_number_found = merge(1,0,signaling_exceptions)
#endif
  !
  if (ncpu > 1) then
#ifdef USE_GATHER
  call mpi_allgather(MPI_IN_PLACE,0_int_mpi,mpi_inttyp,bad_numbers_found, &
                     size(bad_numbers_found,dim=1,kind=int_mpi), &
                     mpi_inttyp,MPI_COMM_WORLD,mpierr)
#else
  call mpi_allreduce(MPI_IN_PLACE,bad_number_found, &
                     size(bad_number_found,kind=int_mpi), &
                     mpi_inttyp,MPI_MAX,MPI_COMM_WORLD,mpierr)
#endif
  end if
  !
#ifdef USE_GATHER
  if (sum(bad_numbers_found) > 0) then
#else
  if (sum(bad_number_found) > 0) then
#endif
    !
    if (len_trim(adjustl(src_file_with_bad_number)) == 0) then
      src_file_with_bad_number = file_name
      line_num_with_bad_number = line_num
    end if
    !
    ! A floating-point exception was detected.
    !
#ifdef USE_GATHER
    !
    if (any(bad_numbers_found(1,:) > 0)) then
      ! F2003 AUTO-REALLOCATION
      signaling_cpus = pack( intseq(0,ncpu-1) , (bad_numbers_found(1,:) > 0) )
      if (ncpu > 1) then
        write (overflow_fmt,201) size(signaling_cpus)-1
        write (overflow_error,overflow_fmt) size(signaling_cpus),ncpu, &
                  (signaling_cpus(n), n=1,size(signaling_cpus))
      else
        write (overflow_error,101)
      end if
    end if
    !
    if (any(bad_numbers_found(2,:) > 0)) then
      ! F2003 AUTO-REALLOCATION
      signaling_cpus = pack( intseq(0,ncpu-1) , (bad_numbers_found(2,:) > 0) )
      if (ncpu > 1) then
        write (divide_by_zero_fmt,202) size(signaling_cpus)-1
        write (divide_by_zero_error,divide_by_zero_fmt) &
                  size(signaling_cpus),ncpu, &
                  (signaling_cpus(n), n=1,size(signaling_cpus))
      else
        write (divide_by_zero_error,102)
      end if
    end if
    !
    if (any(bad_numbers_found(3,:) > 0)) then
      ! F2003 AUTO-REALLOCATION
      signaling_cpus = pack( intseq(0,ncpu-1) , (bad_numbers_found(3,:) > 0) )
      if (ncpu > 1) then
        write (invalid_fmt,203) size(signaling_cpus)-1
        write (invalid_error,invalid_fmt) size(signaling_cpus),ncpu, &
                  (signaling_cpus(n), n=1,size(signaling_cpus))
      else
        write (invalid_error,103)
      end if
    end if
    !
    if (present(ival)) then
      write (output_message,100) trim(adjustl(overflow_error)), &
                                trim(adjustl(divide_by_zero_error)), &
                                trim(adjustl(invalid_error)), &
                                ival
    else
      write (output_message,100) trim(adjustl(overflow_error)), &
                                trim(adjustl(divide_by_zero_error)), &
                                trim(adjustl(invalid_error))
    end if
    !
#else
    if (bad_number_found(1) > 0) write (overflow_error,101)
    if (bad_number_found(2) > 0) write (divide_by_zero_error,102)
    if (bad_number_found(3) > 0) write (invalid_error,103)
    if (present(ival)) then
      write (output_message,100) trim(adjustl(overflow_error)), &
                                trim(adjustl(divide_by_zero_error)), &
                                trim(adjustl(invalid_error)), &
                                ival
    else
      write (output_message,100) trim(adjustl(overflow_error)), &
                                trim(adjustl(divide_by_zero_error)), &
                                trim(adjustl(invalid_error))
    end if
#endif
    !
   !call stop_gfr(stop_mpi,pname,line_num,file_name,output_message)
    !
  end if
  !
  100 format ("A SIGNALING FLOATING POINT EXCEPTION WAS DETECTED!!!", &
             !3(2x,a,:))
              3(2x,a,:),2x,"ival = ",i0)
  !
  101 format ("OVERFLOW EXCEPTION DETECTED!!!")
  102 format ("DIVIDE-BY-ZERO EXCEPTION DETECTED!!!")
  103 format ("INVALID NUMBER (NaN) EXCEPTION DETECTED!!!")
  !
  201 format ("('OVERFLOW EXCEPTION DETECTED ON ',", &
              "i0,' OUT OF ',i0,' PROCESSORS ',", &
              i0,"(i0,','),i0,' !!!')")
  202 format ("('DIVIDE-BY-ZERO EXCEPTION DETECTED ON ',", &
              "i0,' OUT OF ',i0,' PROCESSORS ',", &
              i0,"(i0,','),i0,' !!!')")
  203 format ("('INVALID NUMBER (NaN) EXCEPTION DETECTED ON ',", &
              "i0,' OUT OF ',i0,' PROCESSORS ',", &
              i0,"(i0,','),i0,' !!!')")
  !
  300 format ("A SIGNALING FLOATING POINT EXCEPTION WAS DETECTED!!!", &
              3(2x,a,:),2x,"ival = ",i0)
#endif
#endif
  !
end subroutine check_for_bad_number
!
!###############################################################################
!
