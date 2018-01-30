module module_cgns_types
  !
#ifdef USE_C_NULL_CHAR
#define _CGNS_GOTO_END_ C_NULL_CHAR
#else
#define _CGNS_GOTO_END_ "end"
#endif
  !
  !.. Use Statements
  use, intrinsic :: iso_c_binding, only : C_CHAR,C_NULL_CHAR
  use module_kind_types, only : i4,i8,r4,r8,lk,fals,true
#ifdef CGNS_3_3
  use cgns
#endif
  !
  implicit none
  !
  public
  !
  private :: i4,i8,r4,r8,lk,fals,true
  !
#ifndef CGNS_3_3
  include "cgnslib_f95.h"
  !
  enum, bind(c)
    enumerator :: ElementTypeNull = CG_NULL
    enumerator :: BAR_5 = HEXA_64+1
    enumerator :: TRI_12
    enumerator :: TRI_15
    enumerator :: QUAD_P4_16
    enumerator :: QUAD_25
    enumerator :: TETRA_22
    enumerator :: TETRA_34
    enumerator :: TETRA_35
    enumerator :: PYRA_P4_29
    enumerator :: PYRA_50
    enumerator :: PYRA_55
    enumerator :: PENTA_33
    enumerator :: PENTA_66
    enumerator :: PENTA_75
    enumerator :: HEXA_44
    enumerator :: HEXA_98
    enumerator :: HEXA_125
  end enum
  !
  enum, bind(c)
    enumerator :: BCTypeNull = CG_NULL
  end enum
  !
  enum, bind(c)
    enumerator :: test_enumerator = 0
  end enum
  private :: test_enumerator
  !
  integer, parameter :: cgsize_t = merge(i8,kind(0),CG_BUILD_64BIT==1)
 !integer, parameter :: cgsize_t = i8
  integer, parameter :: cgenum_t = kind(test_enumerator)
  !
#endif
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
  ! Constants to control the CGNS I/O queue
  !
  integer(CBT), parameter :: begin_cgns_queue = 1_cbt
  integer(CBT), parameter :: stop_cgns_queue  = 0_cbt
  !
  ! Last argument required for calls to cg_goto_f
  !
  character(len=*,kind=C_CHAR), parameter :: CGNS_GOTO_END = _CGNS_GOTO_END_
  !
  ! CGNS error variable
  !
  integer(CBT), save :: cgierr
  !
  !
  ! Logical to activate/deactivate using CGNS queue functions
  !
  logical(lk), save :: cgns_use_queue = .false.
  !
contains
!
!###############################################################################
!
subroutine disable_cgns_queue(pname,cgns_grid_file)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: pname
  character(len=*), intent(in) :: cgns_grid_file
  !
continue
  !
#ifndef CGNS_3_3
  if (.not. cgns_use_queue) then
    call cgp_queue_set_f(STOP_CGNS_QUEUE,cgierr)
    call cgns_error(pname,cgns_grid_file,"cgp_queue_set_f", &
                    cgierr,__LINE__,__FILE__)
  end if
#endif
  !
end subroutine disable_cgns_queue
!
!###############################################################################
!
subroutine set_cgns_queue(queue_operation,pname,cgns_grid_file)
  !
  !.. Formal Arguments ..
  integer(CBT),     intent(in) :: queue_operation
  character(len=*), intent(in) :: pname
  character(len=*), intent(in) :: cgns_grid_file
  !
continue
  !
#ifndef CGNS_3_3
  if (cgns_use_queue) then
    call cgp_queue_set_f(BEGIN_CGNS_QUEUE,cgierr)
    call cgns_error(pname,cgns_grid_file,"cgp_queue_set_f", &
                    cgierr,__LINE__,__FILE__)
  end if
#endif
  !
end subroutine set_cgns_queue
!
!###############################################################################
!
subroutine flush_cgns_queue(pname,cgns_grid_file)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: pname
  character(len=*), intent(in) :: cgns_grid_file
  !
continue
  !
#ifndef CGNS_3_3
  if (cgns_use_queue) then
    call cgp_queue_flush_f(cgierr)
    call cgns_error(pname,cgns_grid_file,"cgp_queue_flush_f", &
                    cgierr,__LINE__,__FILE__)
  end if
#endif
  !
end subroutine flush_cgns_queue
!
!###############################################################################
!
subroutine cgns_error(procedure_name,file_name,cgns_routine,iemessag, &
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
  integer(CBT),     intent(in) :: iemessag
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
  if (iemessag == CG_OK) then
    ierror = 0
  else
    ierror = 1
  end if
  if (.not. file_operation_is_serial) then
    call mpi_allreduce(MPI_IN_PLACE,ierror,1_int_mpi,mpi_inttyp, &
                       MPI_MAX,MPI_COMM_WORLD,mpierr)
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
    write (iout,11) trim(adjustl(procedure_name)),iemessag, &
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
                  iemessag,error_message)
  else
    call cg_error_print_f
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
  !
end subroutine cgns_error
!
!###############################################################################
!
end module module_cgns_types
!
!###############################################################################
!
