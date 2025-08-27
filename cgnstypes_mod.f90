module cgns_types_mod
  !
#ifdef USE_C_NULL_CHAR
#define _CGNS_GOTO_END_ C_NULL_CHAR
#else
#define _CGNS_GOTO_END_ "end"
#endif
  !
  !.. Use Statements
  use, intrinsic :: iso_c_binding, only : C_CHAR,C_NULL_CHAR
  use module_kind_types, only : i4,i8,r4,r8,wp,lk,fals,true
  use module_kind_types, only : zero,eps6
  use module_kind_types, only : MPI_COMM_WORLD,mpierr
  use module_kind_types, only : mpi_barrier,mpi_wtime
  use module_kind_types, only : mypnum,glb_root,iout
  use module_kind_types, only : dbg_mpi_frmt
  !
  use cgns
  !
  implicit none
  !
  public
  !
  private :: i4,i8,r4,r8,lk,fals,true
  !
  ! Make cgsize_t and cgenum_t private to force the use of CBT, CST, and CET
  !
  private :: cgsize_t
  private :: cgenum_t
  !
  ! Define the base CGNS integer types
  !
  integer, parameter :: CBT = kind(0)  ! CGNS base integer type
  integer, parameter :: CST = cgsize_t ! CGNS data integer type
  integer, parameter :: CET = cgenum_t ! CGNS enumerator integer type
  !
  integer, parameter :: CGLEN = 32     ! maximum length for a CGNS string
  !
#ifdef SPECIAL_FOR_PGI
  integer, parameter :: CRT = r4       ! CGNS data real type
  integer(kind(REALSINGLE)), parameter :: cgns_kind_real = REALSINGLE
  !
 !integer, parameter :: CRT = r8       ! CGNS data real type
 !integer(kind(REALSINGLE)), parameter :: cgns_kind_real = REALDOUBLE
#else
  integer, parameter :: CRT = r4       ! CGNS data real type
 !integer, parameter :: CRT = r8       ! CGNS data real type
  !
  integer(kind(REALSINGLE)), parameter :: &
      cgns_kind_real = merge(REALDOUBLE,REALSINGLE,CRT == r8)
#endif
  !
  ! Last argument required for calls to cg_goto_f
  !
  character(len=*,kind=C_CHAR), parameter :: CGNS_GOTO_END = _CGNS_GOTO_END_
  !
  ! CGNS error variable
  !
  integer(CBT), save :: cgierr
  !
  ! Integer parameters for CGNS timer
  !
  integer, parameter :: start_cgns_file_timer      = 1
  integer, parameter :: stop_cgns_file_timer       = 2
  integer, parameter :: start_cgns_procedure_timer = 3
  integer, parameter :: stop_cgns_procedure_timer  = 4
  !
contains
!
!###############################################################################
!
subroutine cgns_error(procedure_name,file_name,cgns_routine,cgierr_val, &
                      line_number,src_file_name,optmessag,optinteger, &
                      not_parallel)
  !
  !.. Use Statements ..
  use module_kind_types, only : iout
  use module_kind_types, only : int_mpi
  use module_kind_types, only : mpi_inttyp
  use module_kind_types, only : mpierr
  use module_kind_types, only : mypnum
  use module_kind_types, only : error_message
  use module_kind_types, only : io_error
#ifdef USE_MPI_F08
  use mpi_f08,           only : MPI_IN_PLACE
  use mpi_f08,           only : MPI_MAX
  use mpi_f08,           only : MPI_COMM_WORLD
  use mpi_f08,           only : mpi_allreduce
#else
  use mpi,               only : MPI_IN_PLACE
  use mpi,               only : MPI_MAX
  use mpi,               only : MPI_COMM_WORLD
 !use mpi,               only : mpi_allreduce
#endif
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: procedure_name
  character(len=*), intent(in) :: file_name
  character(len=*), intent(in) :: cgns_routine
  integer(CBT),     intent(in) :: cgierr_val
  integer,          intent(in) :: line_number
  character(len=*), intent(in) :: src_file_name
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: optmessag
  integer,          optional, intent(in) :: optinteger
  logical(lk),      optional, intent(in) :: not_parallel
  !
  !.. Local Scalars ..
  integer :: ierror
  logical(lk) :: file_operation_is_serial
  !
continue
  !
  file_operation_is_serial = fals
  if (present(not_parallel)) then
    file_operation_is_serial = not_parallel
  end if
  !
  if (cgierr_val == CG_OK) then
    ierror = 0
  else
    ierror = 1
  end if
  if (.not. file_operation_is_serial) then
    call mpi_allreduce(MPI_IN_PLACE,ierror,1_int_mpi,mpi_inttyp, &
                       MPI_MAX,MPI_COMM_WORLD,mpierr)
 !else
 !  write (iout,20) trim(adjustl(cgns_routine))
  end if
  if (ierror == 0) return
  !
  if (mypnum == 0) then
    !
    write (iout,10) trim(adjustl(file_name))
    !
    if (present(optmessag) .and. present(optinteger)) then
      write (iout,1) trim(adjustl(optmessag)),optinteger
    else if (present(optmessag)) then
      write (iout,2) trim(adjustl(optmessag))
    else if (present(optinteger)) then
      write (iout,3) optinteger
    end if
    !
    write (iout,11) trim(adjustl(procedure_name)),cgierr_val, &
                    trim(adjustl(cgns_routine))
    !
    if (.not. file_operation_is_serial) then
      write (iout,12) line_number, trim(adjustl(src_file_name))
    end if
    !
    write (iout,13)
    !
  end if
  !
  if (file_operation_is_serial) then
    call cg_get_error_f(error_message)
    call io_error(procedure_name,file_name,5,line_number,src_file_name, &
                  cgierr_val,error_message)
  else
   !call cg_get_error_f(error_message)
    if (mypnum == glb_root) then
      write (iout,30)
      flush (iout)
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    call cg_error_print_f
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      write (iout,31)
      flush (iout)
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    call cgp_error_exit_f
  end if
  !
  ! Format Statements
  !
  1  format (" ",a," = ",i0)
  2  format (" ",a)
  3  format (" ",i0)
  10 format (/," ***************************************************** ",/, &
               " ERROR RETURNED ACCESSING CGNS FILE '",a,"'")
  11 format (" INSIDE PROCEDURE: ",a,/, &
             " ERROR INDICATOR = ",i0,/, &
             " ERROR RETURNED FROM CGNS PROCEDURE: ",a)
  12 format ("                       ON LINE #  ",i0,/, &
             "                     OF THE FILE '",a,"'")
  13 format (" ***************************************************** ",//)
  20 format (2x,"Checking error for: ",a)
  30 format (//,2x,72("#"),//, &
                6x,"BEFORE CALL TO CG_ERROR_PRINT_F",//, &
                2x,72("-"),/)
  31 format (/,2x,72("-"),//, &
               6x,"AFTER  CALL TO CG_ERROR_PRINT_F",//, &
               2x,72("#"),//, &
               6x,"BEFORE CALL TO CG_ERROR_EXIT_F",//, &
               2x,72("-"),/)
  !
end subroutine cgns_error
!
!###############################################################################
!
subroutine cgns_timer(time_flag,profile_io_cgns,message,override_barrier)
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: time_flag
  logical(lk),      intent(in) :: profile_io_cgns
  character(len=*), intent(in) :: message
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: override_barrier
  !
  !.. Local Scalars ..
  real(wp), save :: beg_file_time
  real(wp)       :: end_file_time
  !
  real(wp), save :: beg_procedure_time
  real(wp)       :: end_procedure_time
  !
  real(wp)       :: wall_time
  !
  ! Initialized to zero in case start_cgns_file_timer isnt first use
  ! of this subroutine.
  real(wp), save :: total_sub_time = zero
  !
  logical(lk) :: skip_barrier
  !
continue
  !
  if (.not. profile_io_cgns) return
  !
  skip_barrier = fals
  !
  if (present(override_barrier)) then
    skip_barrier = override_barrier
  end if
  !
  if (.not. skip_barrier) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  if (mypnum == glb_root) then
    !
    if (time_flag == start_cgns_file_timer) then
      !
      total_sub_time = zero
      !
      beg_file_time = mpi_wtime()
     !beg_file_time = dclock()
     !call cpu_time( beg_file_time )
      !
      write (iout,1) trim(adjustl(message))
      !
    else if (time_flag == stop_cgns_file_timer) then
      !
      end_file_time = mpi_wtime()
     !end_file_time = dclock()
     !call cpu_time( end_file_time )
      wall_time = end_file_time - beg_file_time
      !
      if (total_sub_time > eps6) then
        write (iout,5) trim(adjustl(message)), &
                       wall_time-total_sub_time, &
                       total_sub_time, wall_time
      else
        write (iout,2) trim(adjustl(message)),wall_time
      end if
      !
    else if (time_flag == start_cgns_procedure_timer) then
      !
      beg_procedure_time = mpi_wtime()
     !beg_procedure_time = dclock()
     !call cpu_time( beg_procedure_time )
     !if (skip_barrier) then
     !  write (iout,1) trim(adjustl(message))
     !end if
      !
      write (iout,4) trim(adjustl(message))
      !
    else if (time_flag == stop_cgns_procedure_timer) then
      !
      end_procedure_time = mpi_wtime()
      !end_procedure_time = dclock()
      !call cpu_time( end_procedure_time )
      wall_time = end_procedure_time - beg_procedure_time
      !
      total_sub_time = total_sub_time + wall_time
      !
      write (iout,3) trim(adjustl(message)),wall_time
      !
    end if
    !
    flush (iout)
    !
  end if
  !
  if (.not. skip_barrier) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Format Statements
  !
  1 format (/,2x,"Starting I/O timer for: ",a)
  2 format (2x,"Total wall time to write ",a," : ",f12.6," seconds",/)
  3 format (4x,"Time to complete call to ",a," : ",f12.6," seconds")
  4 format (10x,"Starting timer for ",a)
  5 format (2x,"Total wall time to write ",a,/, &
            20x,"Time spent on computations : ",f12.6," seconds",/, &
            20x,"Time spent    in CGNS subs : ",f12.6," seconds",/, &
            20x,"Total wall time   combined : ",f12.6," seconds",/)
  !
end subroutine cgns_timer
!
!###############################################################################
!
end module cgns_types_mod
!
!###############################################################################
!
