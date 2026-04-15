module restart_mod
  !
#include <mpi_defs.h>
!!#define PROFILE_ON
#define _DEBUG_WRITE_TIMER_ global_debug_timer
!#define _DEBUG_WRITE_TIMER_ debug_timer
!!!#define VERBOSE_RESTART
#ifdef VERBOSE_RESTART
#define ONLY_READ_RESTARTS_ARE_VERBOSE
!!!#define ONLY_WRITE_RESTARTS_ARE_VERBOSE
#endif
  !
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  use eqn_idx_mod, only : nq,nec,ncs
  !
  use iso_fortran_env, only : int8,int16,int32,int64
  use iso_fortran_env, only : real32,real64,real128
  !
  implicit none
  !
  private
  !
  !###########################
  !###  Public Procedures  ###
  !###########################
  !
  public :: read_restart_file
  public :: write_restart_file
  public :: check_restart_dir
  public :: create_restart_fname_formats
  public :: deallocate_solver_arrays
  public :: modify_time_average_restart_history
  !
  !###########################
  !###  Module Parameters  ###
  !###########################
  !
  ! Integer type to use for the ordered_geom array
  integer, parameter :: ogp = int32
  !
  integer, parameter :: restart_gauss_nodal_solution   = -Legendre_Gauss
  integer, parameter :: restart_lobatto_nodal_solution = -Legendre_Gauss_Lobatto
  integer, parameter :: restart_legendre_modal_solution = Legendre_Gauss
  integer, parameter :: restart_monomial_modal_solution = Legendre_Gauss_Lobatto
  !
  character(len=*), parameter :: file_datarep = "native"
  !
  !#################################
  !###  Module Scalar Variables  ###
  !#################################
  !
  logical(lk), save :: restart_needs_initialization = true
  logical(lk), save :: time_ave_needs_initialization = true
  !
  logical(lk), save :: restart_operation_is_read = fals
  !
  integer, save :: restart_solution_form
  !
#ifdef PBS_ENV
  _MPI_DATA_TYPE_, save :: ordered_geom_type
  _MPI_DATA_TYPE_, save :: solution_type
  _MPI_DATA_TYPE_, save :: sol_block_type
  _MPI_DATA_TYPE_, save :: averaging_type
  _MPI_DATA_TYPE_, save :: ave_block_type
#else
  _MPI_DATA_TYPE_, save :: ordered_geom_type = MPI_DATATYPE_NULL
  _MPI_DATA_TYPE_, save :: solution_type     = MPI_DATATYPE_NULL
  _MPI_DATA_TYPE_, save :: sol_block_type    = MPI_DATATYPE_NULL
  _MPI_DATA_TYPE_, save :: averaging_type    = MPI_DATATYPE_NULL
  _MPI_DATA_TYPE_, save :: ave_block_type    = MPI_DATATYPE_NULL
#endif
  !
  character(len=170), save :: restart_dir = ""
  character(len=170), save :: rstfile_fmt
  character(len=170), save :: lnkcmnd_fmt
  character(len=170), save :: statfile_fmt
  character(len=170), save :: tafile_fmt
  character(len=170), save :: talink_fmt
  character(len=170), save :: tastat_fmt
  character(len=170), save :: tadone_fmt
  character(len=170), save :: talkdn_fmt
  character(len=170), save :: rstfname_fmt
  character(len=170), save :: tafname_fmt
  character(len=170), save :: tadonefname_fmt
  character(len=170), save :: taiterlk_fmt
  character(len=170), save :: rmtalink_fmt
  !
  integer, save :: averaging_file_switch = 1
  !
  integer, allocatable, save :: output_iterations(:)
  !
  _MPI_DATA_TYPE_, save :: read_type
  _MPI_DATA_TYPE_, save :: read_block_type
  _MPI_DATA_TYPE_, save :: real_type
  !
  integer, save :: num_pts = 0
  integer, save :: num_global_pts = 0
  !
  !##################################
  !###  Module Derived Datatypes  ###
  !##################################
  !
  type :: rstrt_cell_t
    integer :: input_cell
    integer :: input_beg_sp
    integer :: input_end_sp
    integer :: input_geom
    integer :: input_order
  end type rstrt_cell_t
  !
  type :: sp_idx_t
    integer :: beg_sp
    integer :: end_sp
  end type sp_idx_t
  !
  type :: backup_t
    character(len=200) :: filename
    integer            :: iter
    real(wp)           :: cfl
  end type backup_t
  !
  type :: backup_restart_t
    character(len=200) :: filename
    integer            :: iter
    real(wp)           :: cfl
    type(backup_t)     :: backup
  end type backup_restart_t
  !
  ! TA_DATA_T : Derived type to store all the data from a restart file.
  !             Currently, this is only used by the utility to modify the
  !             history of a time-average restart file.
  type :: ta_data_t
    real(wp),         allocatable :: sol(:,:)
    integer(ogp),     allocatable :: ordered_geom(:)
    integer(wip),     allocatable :: int_header(:)
    real(wp),         allocatable :: real_header(:)
    integer(wip),     allocatable :: var_header(:)
    character(len=:), allocatable :: ta_filename
  end type ta_data_t
  !
  !###################################
  !###  Module Allocatable Arrays  ###
  !###################################
  !
  integer, save, allocatable :: cell_output_order(:)
  !
  integer(INT_MPI), save, allocatable :: og_counts(:)
  integer(INT_MPI), save, allocatable :: og_displs(:)
  !
  integer(INT_MPI), save, allocatable :: pts_counts(:)
  integer(INT_MPI), save, allocatable :: pts_displs(:)
  !
  integer(int_mpi), save, allocatable :: pts_offset(:)
  integer(int_mpi), save, allocatable :: pts_blklen(:)
  !
  type(sp_idx_t), save, allocatable :: host_output_order(:)
  !
  type(rstrt_cell_t), save, allocatable :: cell_input_data(:)
  !
  type(backup_restart_t), save :: backup_restart
  !
  ! ORIGINAL_TA_DATA : This is used by the utility for modifying the history
  !                    of a time-average restart file to store the original
  !                    time-average data that is ultimately being modified,
  !                    in case it will be modified by multiple other
  !                    time-average restart files.
  type(ta_data_t), allocatable :: original_ta_data
  !
  !###########################################################################
  !#  Generic interfaces for subroutines to read INTEGER data using MPI I/O  #
  !###########################################################################
  !
  interface read_rstrt_file_i1
    module procedure read_rstrt_file_int8_data_i4, &
                     read_rstrt_file_int8_data_i8
  end interface read_rstrt_file_i1
  !
  interface read_rstrt_file_i2
    module procedure read_rstrt_file_int16_data_i4, &
                     read_rstrt_file_int16_data_i8
  end interface read_rstrt_file_i2
  !
  interface read_rstrt_file_i4
    module procedure read_rstrt_file_int32_data_i4, &
                     read_rstrt_file_int32_data_i8
  end interface read_rstrt_file_i4
  !
  interface read_rstrt_file_i8
    module procedure read_rstrt_file_int64_data_i4, &
                     read_rstrt_file_int64_data_i8
  end interface read_rstrt_file_i8
  !
  !###########################################################################
  !# Generic interfaces for procedures to read INTEGER data using Serial I/O #
  !###########################################################################
  !
  interface read_rstrt_file_serial_i1
    module procedure read_rstrt_file_serial_int8_data_i4, &
                     read_rstrt_file_serial_int8_data_i8
  end interface read_rstrt_file_serial_i1
  !
  interface read_rstrt_file_serial_i2
    module procedure read_rstrt_file_serial_int16_data_i4, &
                     read_rstrt_file_serial_int16_data_i8
  end interface read_rstrt_file_serial_i2
  !
  interface read_rstrt_file_serial_i4
    module procedure read_rstrt_file_serial_int32_data_i4, &
                     read_rstrt_file_serial_int32_data_i8
  end interface read_rstrt_file_serial_i4
  !
  interface read_rstrt_file_serial_i8
    module procedure read_rstrt_file_serial_int64_data_i4, &
                     read_rstrt_file_serial_int64_data_i8
  end interface read_rstrt_file_serial_i8
  !
  !############################################################################
  !###  Generic interfaces for subroutines to read REAL data using MPI I/O  ###
  !############################################################################
  !
  interface read_rstrt_file_r4
    module procedure read_rstrt_file_r4_rank1, &
                     read_rstrt_file_r4_rank2
  end interface read_rstrt_file_r4
  !
  interface read_rstrt_file_r8
    module procedure read_rstrt_file_r8_rank1, &
                     read_rstrt_file_r8_rank2
  end interface read_rstrt_file_r8
  !
  interface read_rstrt_file_r16
    module procedure read_rstrt_file_r16_rank1, &
                     read_rstrt_file_r16_rank2
  end interface read_rstrt_file_r16
  !
  interface read_rstrt_file_wp
    module procedure read_rstrt_file_wp_rank1, &
                     read_rstrt_file_wp_rank2
  end interface read_rstrt_file_wp
  !
  !############################################################################
  !### Generic interfaces for procedures to read REAL data using Serial I/O ###
  !############################################################################
  !
  interface read_rstrt_file_serial_r4
    module procedure read_rstrt_file_serial_r4_rank1, &
                     read_rstrt_file_serial_r4_rank2
  end interface read_rstrt_file_serial_r4
  !
  interface read_rstrt_file_serial_r8
    module procedure read_rstrt_file_serial_r8_rank1, &
                     read_rstrt_file_serial_r8_rank2
  end interface read_rstrt_file_serial_r8
  !
  interface read_rstrt_file_serial_r16
    module procedure read_rstrt_file_serial_r16_rank1, &
                     read_rstrt_file_serial_r16_rank2
  end interface read_rstrt_file_serial_r16
  !
  interface read_rstrt_file_serial_wp
    module procedure read_rstrt_file_serial_wp_rank1, &
                     read_rstrt_file_serial_wp_rank2
  end interface read_rstrt_file_serial_wp
  !
contains
!
!###############################################################################
!
subroutine restart_from_backup()
  !
  !.. Use Statements ..
  use ovar, only : itcur,time,d_t
  use ovar, only : itrst,tmrst,dtrst
  use ovar, only : cfl,cfl_beg,cfl_end,cfl_cycles
  use ovar, only : restart_interval
  use ovar, only : CFL_reduction_for_NaN_restart
  use ovar, only : adjust_cfl_cycles
  !
  !.. Local Scalars ..
  integer  :: iter,min_iter
  real(wp) :: restart_cfl
  !
  character(len=200) :: restart_fname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "restart_from_backup"
  !
continue
  !
  call global_debug_timer(entering_procedure,pname)
  !
  iter = itcur - backup_restart%iter
  min_iter = min( 1000 , nint(real(restart_interval,kind=wp)/ten) )
  !
  if (iter < min_iter) then
    !
    backup_restart%filename = backup_restart%backup%filename
    backup_restart%iter     = backup_restart%backup%iter
    backup_restart%cfl      = backup_restart%backup%cfl
    !
    backup_restart%backup%cfl = backup_restart%cfl*CFL_reduction_for_NaN_restart
    !
  end if
  !
  restart_fname = backup_restart%filename
  !
  call read_restart_file(restart_fname)
  !
  d_t   = dtrst
  itcur = itrst
  time  = tmrst
  !
  cfl = backup_restart%cfl*CFL_reduction_for_NaN_restart
  !
  backup_restart%cfl = cfl
  !
  cfl_beg = cfl
  cfl_end = cfl
  cfl_cycles = 0
  call adjust_cfl_cycles
  !
  call global_debug_timer(leaving_procedure,pname)
  !
end subroutine restart_from_backup
!
!###############################################################################
!
subroutine read_restart_file(restart_fname,force_serial_read)
  !
  !.. Use Statements ..
  use ovar,        only : itrst,tmrst,dtrst,rl2mx
  use ovar,        only : prev_ave_time,ave_start_time
  use ovar,        only : restart_file,time_ave_file
  use ovar,        only : output_time_averaging
  use ovar,        only : read_time_ave_restart
  use ovar,        only : modify_time_ave_restart_history
  use flowvar_mod, only : usp,uavesp,time_ave_variables
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: restart_fname
  logical(lk),      optional, intent(in) :: force_serial_read
  !
  !.. Local Scalars ..
  integer :: n,ierr
  character(len=150) :: fname
  !
  !.. Local Allocatable Arrays ..
  integer(wip), allocatable :: int_header(:)
  integer(wip), allocatable :: var_header(:)
  real(wp),     allocatable :: real_header(:)
  integer(ogp), allocatable :: ordered_geom(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_restart_file"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  restart_operation_is_read = true
  !
  fname = trim(adjustl(restart_file))
  if (present(restart_fname)) fname = trim(adjustl(restart_fname))
  !
  !
  ! Read the restart solution file
  ! NOTE: usp SHOULD BE ALLOCATED AT THIS POINT
  !
  call read_solution_file(fname,int_header,real_header,usp,ordered_geom, &
                          force_serial_read=force_serial_read)
  !
  ! Extract the required information from the header arrays
  !
  itrst = int( int_header(4) , kind=kind(itrst) )
  !
  tmrst = real( real_header(1) , kind=kind(tmrst) )
  dtrst = real( real_header(2) , kind=kind(dtrst) )
  !
  if (size(real_header) >= 3) then
    rl2mx = real( real_header(3) , kind=kind(rl2mx) )
  end if
  !
  if (output_time_averaging .or. modify_time_ave_restart_history) then
    !
    if (read_time_ave_restart .or. modify_time_ave_restart_history) then
      !
      fname = trim(adjustl(time_ave_file))
      !
      ! NOTE: uavesp SHOULD NOT BE ALLOCATED AT THIS POINT
      !
      call read_solution_file(fname,int_header,real_header,uavesp, &
                              ordered_geom, &
                              var_header=var_header, &
                              force_serial_read=force_serial_read)
      !
      prev_ave_time = real( real_header(2) , kind=kind(prev_ave_time) )
      ave_start_time = real( real_header(3) , kind=kind(ave_start_time) )
      !
      if (size(real_header) >= 4) then
        rl2mx = real( real_header(4) , kind=kind(rl2mx) )
      end if
      !
      if (int_header(4) /= itrst) then
        write (error_message,1) int_header(4),itrst
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      !
      if (abs(real_header(1)-tmrst) > eps9) then
        write (error_message,2) real_header(1),tmrst
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      !
      if (allocated(time_ave_variables)) then
        deallocate ( time_ave_variables , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"time_ave_variables",2,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
      !
      allocate ( time_ave_variables(1:size(var_header)) , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"time_ave_variables",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      do n = 1,size(var_header)
        time_ave_variables(n) = transfer(var_header(n),time_ave_variables(n))
        time_ave_variables(n) = trim(adjustl(time_ave_variables(n)))
      end do
      !
      if (modify_time_ave_restart_history) then
        !
        ! If running the utility to modify the history of a time-average
        ! restart file, save all of the data from this time-average restart
        ! file to the module derived type original_ta_data so it can be
        ! used again later in the utility to write new modified
        ! time-average restart files.
        !
        call save_ta_data(original_ta_data, &
                          ta_filename=fname, &
                          int_header=int_header, &
                          real_header=real_header, &
                          var_header=var_header, &
                          ordered_geom=ordered_geom, &
                          solution=uavesp)
        !
      end if
      !
    else
      !
      ! Reset rl2mx since we are starting a new time averaging
      !
      rl2mx = -huge(zero)
      !
    end if
    !
  end if
  !
  ! Deallocate ordered_geom just to be sure it is done. The other local
  ! allocatable arrays are small enough to not care.
  !
  if (allocated(ordered_geom)) then
    deallocate ( ordered_geom , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ordered_geom",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  restart_operation_is_read = fals
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format("The current iteration number in the time averaging file does ", &
           "not match the iteration number in the restart file: ", &
           "time_ave_iter = ",i0,", restart_iter = ",i0)
  2 format("The current solution time in the time averaging file does ", &
           "not match the solution time in the restart file: ", &
           "time_ave_time = ",es14.6,", restart_time = ",es14.6)
  !
end subroutine read_restart_file
!
!###############################################################################
!
subroutine read_solution_file(file_name,int_header,real_header,solution, &
                              ordered_geom,var_header,ignore_sem, &
                              force_serial_read)
  !
  !.. Formal Arguments ..
  character(len=*),          intent(in   ) :: file_name
  integer(wip), allocatable, intent(  out) :: int_header(:)
  real(wp),     allocatable, intent(  out) :: real_header(:)
  real(wp),     allocatable, intent(inout) :: solution(:,:)
  integer(ogp), allocatable, intent(inout) :: ordered_geom(:)
  !
  !.. Optional Arguments ..
  integer(wip), optional, allocatable, intent(  out) :: var_header(:)
  logical(lk),  optional,              intent(in   ) :: ignore_sem
  logical(lk),  optional,              intent(in   ) :: force_serial_read
  !
  !.. Local Scalars ..
  logical(lk) :: force_serial
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_solution_file"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  force_serial = fals
  if (present(force_serial_read)) force_serial = force_serial_read
  !
  if (ncpu == 1 .or. force_serial) then
    call read_solution_file_serial(file_name,int_header,real_header,solution, &
                                   ordered_geom,var_header,ignore_sem)
  else
    call read_solution_file_mpi(file_name,int_header,real_header,solution, &
                                ordered_geom,var_header,ignore_sem)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_solution_file
!
!###############################################################################
!
subroutine read_solution_file_mpi(file_name,int_header,real_header,solution, &
                                  ordered_geom,var_header,ignore_sem)
  !
  !.. Use Statements ..
  use geomvar_mod, only : ncell,cell,n_solpts
  use geomvar_mod, only : n_global_cell,n_global_solpts
  use ovar,        only : time_average_restart_files
  use thermo_mod,  only : inert_species_map
  !
  !.. Formal Arguments ..
  character(len=*),          intent(in   ) :: file_name
  integer(wip), allocatable, intent(  out) :: int_header(:)
  real(wp),     allocatable, intent(  out) :: real_header(:)
  real(wp),     allocatable, intent(inout) :: solution(:,:)
  integer(ogp), allocatable, intent(inout) :: ordered_geom(:)
  !
  !.. Optional Arguments ..
  integer(wip), optional, allocatable, intent(  out) :: var_header(:)
  logical(lk),  optional,              intent(in   ) :: ignore_sem
  !
  !.. Local Scalars ..
  integer :: l0,m,mbeg,mm,n,np,n0,ierr
  integer :: sol_nq,ncs_diff
  integer :: inert_nq,inert_ncs
  integer :: restart_ip,restart_rp
  logical(lk) :: rst_form_matches_sol_form
  logical(lk) :: skip_sem_restart_data
  !
  _MPI_FILE_TYPE_ :: file_handle
  _MPI_INFO_TYPE_ :: file_info
  !.. Need to be MPI_OFFSET_KIND integers
  integer(MPI_OFFSET_KIND) :: file_disp
  integer(MPI_OFFSET_KIND) :: file_offset
  integer(MPI_OFFSET_KIND) :: real_bytes
  !
  integer :: solution_errors
  integer :: restart_nq,restart_ncs,restart_ncells
  integer :: n_header_int,n_header_real
  integer :: num_sem_boxes
  !
  integer :: inp_cell,inp_geom,inp_order
  integer :: this_cell,this_geom,this_order
  !
  ! Variables for restart status
  integer :: ios
  integer, save :: read_counts = 0
  !
  character(len=25) :: statfile
  !
  _MPI_COMM_TYPE_ :: file_comm
  _MPI_STATUS_TYPE_ :: istat(1:_MPI_STATUS_SIZE_)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: rst_sol(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_solution_file_mpi"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  skip_sem_restart_data = fals
  if (present(ignore_sem)) then
    skip_sem_restart_data = ignore_sem
  end if
  !
  num_sem_boxes = 0
  !
  solution_errors = 0
  !
  read_counts = read_counts + 1
  write (statfile,'("read_restart.",i0,".stat")') read_counts
  file_comm = MPI_COMM_WORLD
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            open_file=true, &
                            stat_strng="init stat, before mpi_info")
  !
  ! Initialize file_info to MPI_INFO_NULL
  !
  file_info = MPI_INFO_NULL
  !
#ifdef PBS_ENV
#ifndef IGNORE_MPI_INFO
  !
  ! Create the MPI info object for this restart file
  ! NOTE: For now, we will leave the MPI info object as is without any extra
  !       hints since this will just be for reading the restart file this one
  !       time before starting the simulation.
  !
  call debug_timer(start_timer,"calling mpi_info_create")
  call mpi_info_create(file_info,mpierr)
  call debug_timer(stop_timer,"calling mpi_info_create")
  !
  ! Added read_once to the MPI info object just in case this helps
  !
  call mpi_info_set(file_info,"access_style","read_once",mpierr)
  !
  ! Added disabling of romio_cb_read (reading with collective buffering) because
  ! this was causing hangs on Rome nodes while reading the restart header info.
  !
  call mpi_info_set(file_info,"romio_cb_read","disable",mpierr)
  !
#endif
#endif
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                     stat_strng="after mpi_info, before mpi_file_open")
  !
  ! Open the restart file for MPI I/O
  !
  call debug_timer(start_timer,"calling mpi_file_open")
  call mpi_file_open(MPI_COMM_WORLD,trim(adjustl(file_name)), &
                     MPI_MODE_RDONLY,file_info,file_handle,mpierr)
  call debug_timer(stop_timer,"calling mpi_file_open")
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            stat_strng="after mpi_file_open")
  !
  ! Have all processors read the header section of the restart file
  !
! ######################################
! ### Integer part of restart header ###
! ######################################
  !
  ! Read the integer type for the integer part of the restart header
  !
  call debug_timer(start_timer,"reading int type for int header")
  call mpi_file_read_all(file_handle,restart_ip,1_int_mpi, &
                         MPI_INTEGER,istat(1),mpierr)
  call debug_timer(stop_timer,"reading int type for int header")
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            mpistat=istat, &
                            stat_strng="after reading int type for int header")
  !
  ! Read the number of integer header items
  !
  call debug_timer(start_timer,"reading # of int header items")
  call mpi_file_read_all(file_handle,n_header_int,1_int_mpi, &
                         MPI_INTEGER,istat(1),mpierr)
  call debug_timer(stop_timer,"reading # of int header items")
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            mpistat=istat, &
                            stat_strng="after reading # of int header items")
  !
  ! Allocate the int_header array
  !
  if (allocated(int_header)) then
    deallocate ( int_header , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"int_header",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( int_header(1:n_header_int) , source=0_wip , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"int_header",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read the integer part of the restart header
  !
  call debug_timer(start_timer,"reading integer header")
  select case (restart_ip)
    case (int32) ; call read_rstrt_file_i4(file_handle,int_header, &
                                           mpistat=istat)
    case (int64) ; call read_rstrt_file_i8(file_handle,int_header, &
                                           mpistat=istat)
    case (int8)  ; call read_rstrt_file_i1(file_handle,int_header, &
                                           mpistat=istat)
    case (int16) ; call read_rstrt_file_i2(file_handle,int_header, &
                                           mpistat=istat)
    case default
      write (error_message,20)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end select
  call debug_timer(stop_timer,"reading integer header")
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            mpistat=istat, &
                            stat_strng="after reading integer header")
  !
  ! Extract the information from the integer section of the header
  !
  ! Extract the number of flow variables saved at each solution point
  ! within the restart file
  !
  restart_nq = int( int_header(1) , kind=kind(restart_nq) )
  !
  if (allocated(solution)) then
    !
    sol_nq = size(solution,dim=2)
    !
    if (ncs > 0) then
      !
      inert_ncs = size(inert_species_map)
      inert_nq = (nq-nec+1) + inert_ncs
      restart_ncs = restart_nq - (nq-nec+1)
      !
      if (all(restart_ncs /= [ncs,inert_ncs])) then
        write (error_message,50) restart_ncs,ncs,inert_ncs
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      !
      if (all(restart_nq /= [sol_nq,inert_nq])) then
        write (error_message,51) restart_nq,sol_nq,inert_nq
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      !
    else if (restart_nq /= sol_nq) then
      !
      write (error_message,10) restart_nq,sol_nq
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      !
    end if
    !
  end if
  !
  ! Extracted the total number of grid cells contained within the restart file
  !
  restart_ncells = int( int_header(2) , kind=kind(restart_ncells) )
  !
  if (restart_ncells /= n_global_cell) then
    write (error_message,11) restart_ncells,n_global_cell
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Extract the format that the solution data is in within the restart file
  !
  restart_solution_form = int( int_header(3) , &
                               kind=kind(restart_solution_form) )
!%!
!%call check_restart_solution_form(restart_solution_form)
  !
  ! Check if there is restart information for
  ! the synthetic inflow turbulence generators
  !
  if (size(int_header) >= 5) then
    num_sem_boxes = max(0,int_header(5))
  end if
  !
! ###################################
! ### Real part of restart header ###
! ###################################
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                          stat_strng="before reading real type for real header")
  !
  ! Read the real type for the real part of the restart header
  !
  call debug_timer(start_timer,"reading real type for real header")
  call mpi_file_read_all(file_handle,restart_rp,1_int_mpi, &
                         MPI_INTEGER,istat(1),mpierr)
  call debug_timer(stop_timer,"reading real type for real header")
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            mpistat=istat, &
                           stat_strng="after reading real type for real header")
  !
  ! Read the number of real header items
  !
  call debug_timer(start_timer,"reading # of real header items")
  call mpi_file_read_all(file_handle,n_header_real,1_int_mpi, &
                         MPI_INTEGER,istat(1),mpierr)
  call debug_timer(stop_timer,"reading # of real header items")
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            mpistat=istat, &
                            stat_strng="after reading # of real header items")
  !
  ! Allocate the real_header array
  !
  if (allocated(real_header)) then
    deallocate ( real_header , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"real_header",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( real_header(1:n_header_real) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"real_header",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read the real part of the restart header
  !
  call debug_timer(start_timer,"reading real header")
  if (restart_rp == wp) then
    call read_rstrt_file_wp(file_handle,real_header,mpistat=istat)
  else
    select case (restart_rp)
      case (real64)  ; call read_rstrt_file_r8(file_handle,real_header, &
                                               mpistat=istat)
      case (real32)  ; call read_rstrt_file_r4(file_handle,real_header, &
                                               mpistat=istat)
      case (real128) ; call read_rstrt_file_r16(file_handle,real_header, &
                                                mpistat=istat)
      case default
        write (error_message,21)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end select
  end if
  call debug_timer(stop_timer,"reading real header")
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            mpistat=istat, &
                            stat_strng="after reading real header")
  !
! #############################################
! ### Ordered_Geom part of the restart file ###
! #############################################
  !
  ! Read the integer type for the ordered_geom array
  !
  call debug_timer(start_timer,"reading int type for ordered_geom")
  call mpi_file_read_all(file_handle,restart_ip,1_int_mpi, &
                         MPI_INTEGER,istat(1),mpierr)
  call debug_timer(stop_timer,"reading int type for ordered_geom")
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            mpistat=istat, &
                           stat_strng="after reading int type for ordered_geom")
  !
  if (.not. allocated(ordered_geom)) then
    allocate ( ordered_geom(1:restart_ncells) , source=0_ogp , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ordered_geom",1,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  ! All processors have to read the ordered_geom data so that each
  ! processor can compute the data offsets for each of the cells
  ! within its partition.
  !
  call debug_timer(start_timer,"reading ordered_geom")
  select case (restart_ip)
    case (int32) ; call read_rstrt_file_i4(file_handle,ordered_geom, &
                                           mpistat=istat)
    case (int64) ; call read_rstrt_file_i8(file_handle,ordered_geom, &
                                           mpistat=istat)
    case (int8)  ; call read_rstrt_file_i1(file_handle,ordered_geom, &
                                           mpistat=istat)
    case (int16) ; call read_rstrt_file_i2(file_handle,ordered_geom, &
                                           mpistat=istat)
    case default
      write (error_message,22)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end select
  call debug_timer(stop_timer,"reading ordered_geom")
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            mpistat=istat, &
                            stat_strng="after reading ordered_geom")
  !
  ! Check to see if the form of the restart data matches the
  ! form of the solution data for the current simulation.
  !
  call check_restart_solution_form(restart_solution_form,ordered_geom, &
                                   rst_form_matches_sol_form)
  !
! ################################################
! ### Variable header part of the restart file ###
! ################################################
  !
  if (present(var_header)) then
    !
    ! Read the integer type for the variable header
    !
    call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            stat_strng="before reading int type for var_header")
    !
    call debug_timer(start_timer,"reading int type for variable header")
    call mpi_file_read_all(file_handle,restart_ip,1_int_mpi, &
                           MPI_INTEGER,istat(1),mpierr)
    call debug_timer(stop_timer,"reading int type for variable header")
    !
    call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                              mpistat=istat, &
                             stat_strng="after reading int type for var_header")
    !
    if (allocated(var_header)) then
      deallocate ( var_header , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"var_header",2,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
    !
    allocate ( var_header(1:restart_nq) , source=0_wip , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"var_header",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Read the variable header for the solution data
    !
    call debug_timer(start_timer,"reading variable header")
    select case (restart_ip)
      case (int32) ; call read_rstrt_file_i4(file_handle,var_header, &
                                             mpistat=istat)
      case (int64) ; call read_rstrt_file_i8(file_handle,var_header, &
                                             mpistat=istat)
      case (int8)  ; call read_rstrt_file_i1(file_handle,var_header, &
                                             mpistat=istat)
      case (int16) ; call read_rstrt_file_i2(file_handle,var_header, &
                                             mpistat=istat)
      case default
        write (error_message,24)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end select
    call debug_timer(stop_timer,"reading variable header")
    !
    call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                              mpistat=istat, &
                              stat_strng="after reading variable header")
    !
    if (.not. allocated(solution)) then
      allocate ( solution(1:n_solpts,1:restart_nq) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"solution",1,__LINE__,__FILE__,ierr,error_message)
    end if
    !
  end if
  !
! #########################################
! ### Solution part of the restart file ###
! #########################################
  !
  ! Read the real type for the solution data
  !
  call debug_timer(start_timer,"reading real type for solution data")
  call mpi_file_read_all(file_handle,restart_rp,1_int_mpi, &
                         MPI_INTEGER,istat(1),mpierr)
  call debug_timer(stop_timer,"reading real type for solution data")
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            mpistat=istat, &
                         stat_strng="after reading real type for solution data")
  !
  ! Get the updated file displacement
  !
  call debug_timer(start_timer,"getting file position to read solution data")
  call mpi_file_get_position(file_handle,file_offset,mpierr)
  call debug_timer(stop_timer,"getting file position to read solution data")
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
          stat_strng="after getting file position before reading solution data")
  !
  call debug_timer(start_timer,"getting byte offset to read solution data")
  call mpi_file_get_byte_offset(file_handle,file_offset,file_disp,mpierr)
  call debug_timer(stop_timer,"getting byte offset to read solution data")
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
            stat_strng="after getting byte offset before reading solution data")
  !
  ! Get the MPI derived types required to read the solution data section
  ! of the restart file if cell_input_data is not allocated
  !
  if (.not. allocated(cell_input_data)) then
    !
    call debug_timer(start_timer,"calling get_input_mpi_datatypes")
    call get_input_mpi_datatypes(restart_nq,restart_rp,real_type, &
                                 read_block_type,read_type, &
                                 num_pts,num_global_pts, &
                                 ordered_geom,cell_input_data)
    call debug_timer(stop_timer,"calling get_input_mpi_datatypes")
    !
  end if
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                  stat_strng="after getting MPI datatype to read solution data")
  !
!%!
!%! IF WE ARE TO MOVE THE POINT IN THE PRE-PROCESSING PHASE AT WHICH THE
!%! RESTART FILE IS READ, WE NEED TO MAKE SURE THAT N_SOLPTS HAS BEEN
!%! RECOMPUTED SO THAT IT IS LOCAL TO EACH PROCESSOR BECAUSE THIS IS NOT
!%! DONE WITHIN PARALLEL_MOD AFTER PARTITIONING THE GRID.
!%!
!%if (num_pts /= n_solpts) then
!%  write (iout,12) mypnum,num_pts,n_solpts
!%  solution_errors = 100
!%end if
!%!
!%! Sum up the value of solution_errors across all processors
!%!
!%call mpi_allreduce(MPI_IN_PLACE,solution_errors,1_int_mpi, &
!%                   mpi_inttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
!%!
!%if (solution_errors /= 0) then
!%  write (error_message,13)
!%  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
!%end if
  !
  ! Set the file view to to read in the solution restart data
  !
  call debug_timer(start_timer,"setting file view to read solution data")
  call mpi_file_set_view(file_handle,file_disp,read_block_type, &
                         read_type, trim(adjustl(file_datarep)), &
                         file_info, mpierr)
  call debug_timer(stop_timer,"setting file view to read solution data")
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
              stat_strng="after setting file view before reading solution data")
  !
  ! Allocate the rst_sol array to read in the solution restart data
  !
!%allocate ( rst_sol(1:restart_nq,1:n_solpts) , source=zero , &
  allocate ( rst_sol(1:restart_nq,1:num_pts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"rst_sol",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read in the solution restart data
  !
  call debug_timer(start_timer,"reading solution data")
  if (restart_rp == wp) then
    call read_rstrt_file_wp(file_handle,rst_sol,read_block_type,mpistat=istat)
  else
    select case (restart_rp)
      case (real64)  ; call read_rstrt_file_r8(file_handle,rst_sol, &
                                               read_block_type, &
                                               mpistat=istat)
      case (real32)  ; call read_rstrt_file_r4(file_handle,rst_sol, &
                                               read_block_type, &
                                               mpistat=istat)
      case (real128) ; call read_rstrt_file_r16(file_handle,rst_sol, &
                                                read_block_type, &
                                                mpistat=istat)
      case default
        write (error_message,23)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end select
  end if
  call debug_timer(stop_timer,"reading solution data")
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            mpistat=istat, &
                            stat_strng="after reading solution data")
  !
  ! Read in the restart information for the synthetic
  ! inflow turbulence generators
  !
  if (.not. skip_sem_restart_data) then
    !
    if (num_sem_boxes > 0) then
      !
      ! Update the file displacement since we know the total number of
      ! solution points, the number of variables per solution point, and
      ! the byte-size of for each element of the solution array.
      !
      select case (restart_rp)
        case (real32)  ; real_bytes = int(4,kind=MPI_OFFSET_KIND)
        case (real64)  ; real_bytes = int(8,kind=MPI_OFFSET_KIND)
        case (real128) ; real_bytes = int(16,kind=MPI_OFFSET_KIND)
      end select
      !
      file_disp = file_disp + real_bytes * &
                              int(int_header(1),kind=MPI_OFFSET_KIND) * &
                              int(num_global_pts,kind=MPI_OFFSET_KIND)
      !
      ! Update the file view to after all the solution data
      !
      call debug_timer(start_timer,"setting file view to read sem data")
      call mpi_file_set_view(file_handle,file_disp,MPI_INTEGER, &
                             MPI_INTEGER, trim(adjustl(file_datarep)), &
                             file_info, mpierr)
      call debug_timer(stop_timer,"setting file view to read sem data")
      !
      call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                  stat_strng="after setting file view before reading sem data")
      !
      ! Read in the SEM data from the restart file
      !
      call debug_timer(start_timer,"calling read_sem_restart_mpi")
      call read_sem_restart_mpi(num_sem_boxes,file_handle)
      call debug_timer(stop_timer,"calling read_sem_restart_mpi")
      !
      call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                  stat_strng="after reading all sem data")
      !
    end if
    !
  end if
  !
  ! Finally close the restart file
  !
  call debug_timer(start_timer,"calling mpi_file_close")
  call mpi_file_close(file_handle,mpierr)
  call debug_timer(stop_timer,"calling mpi_file_close")
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            close_file=true, &
                            stat_strng="after mpi_file_close")
  !
  !
  !
  if (rst_form_matches_sol_form) then
    !
    ! Copy the restart solution data over into the solution array
    !
    solution_errors = 0
    !
    do this_cell = 1,ncell
      !
      inp_cell = cell_input_data(this_cell)%input_cell
      !
      inp_geom = cell_input_data(this_cell)%input_geom
      inp_order = cell_input_data(this_cell)%input_order
      !
      this_geom = cell(this_cell)%geom
      this_order = cell(this_cell)%order
      !
      ! NEED TO ADD SOMETHING HERE TO CONVERT THE SOLUTION IF THE SOLUTION
      ! FORM DOESNT MATCH THE CURRENT SIMULATION.  FOR NOW, JUST CREATE AN
      ! ERROR IF IT DOESNT MATCH THE RESTART FORM.
      !
      if ((inp_geom /= this_geom) .or. (inp_order /= this_order)) then
        write (iout,30) mypnum, inp_cell, inp_geom, inp_order, &
                               this_cell,this_geom,this_order
        solution_errors = solution_errors + 1
      end if
      !
      n0 = cell(this_cell)%beg_sp - 1
      np = cell(this_cell)%end_sp - n0
      !
      l0 = cell_input_data(this_cell)%input_beg_sp - 1
      !
      ! Copy over the restart solutions for the continuity, momentum, and
      ! energy equations. Also copy over the restart solutions for the
      ! species equations if the number of chemical species in the restart
      ! file matches the number of chemical species in the current simulation.
      !
      ncs_diff = 0
      mbeg = 1
      if (ncs > 0 .and. restart_ncs /= ncs) then
        ncs_diff = ncs - restart_ncs
        mbeg = restart_ncs + 1
      end if
      !
      ! We need to flip the indices between the two arrays to match how the
      ! restart solution was read using the MPI derived datatype.
      !
      do n = 1,np
        do m = mbeg,restart_nq
          solution(n0+n,m+ncs_diff) = rst_sol(m,l0+n)
        end do
      end do
      !
      ! Copy over the restart solutions for the species equations if 
      ! the number of chemical species in the restart file does not
      ! match the number of chemical species in the current simulation.
      !
      if (ncs > 0 .and. restart_ncs /= ncs) then
        !
        if (restart_ncs == inert_ncs) then
          !
          ! We need to flip the indices between the two arrays to match how the
          ! restart solution was read using the MPI derived datatype.
          !
          do n = 1,np
            do mm = 1,restart_ncs
              m = inert_species_map(mm)
              solution(n0+n,m) = rst_sol(mm,l0+n)
            end do
          end do
          !
        else
          !
          write (error_message,50) restart_ncs,ncs,inert_ncs
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
          !
        end if
        !
      end if
      !
    end do
    !
  else
    !
    ! Interpolate/project the restart solution
    ! to the current solution points/order
    !
    call interp_restart_solution(restart_solution_form,rst_sol,solution)
    !
  end if
  !
  ! Deallocate rst_sol since it is no longer needed
  !
  deallocate ( rst_sol , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"rst_sol",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Sum up the value of solution_errors across all processors
  !
  call mpi_allreduce(MPI_IN_PLACE,solution_errors,1_int_mpi, &
                     mpi_inttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
  !
  ! Stop the simulation if there were any errors copying the restart
  ! solution data into solution.
  !
  if (solution_errors > 0) then
    write (error_message,31) solution_errors
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Deallocate cell_input_data if we are not time-averaging a bunch of
  ! restart files. Also reset the variables real_type, read_type, and
  ! num_pts to default values in case we ever need to come back here
  ! and read in a new restart file.
  !
  if (.not. time_average_restart_files) then
    !
    if (allocated(cell_input_data)) then
      deallocate ( cell_input_data , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"cell_input_data",2,__LINE__,__FILE__,ierr, &
                       error_message)
    end if
    !
    real_type = MPI_DATATYPE_NULL
    read_type = MPI_DATATYPE_NULL
    num_pts = 0
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  10 format("The number of variables in the restart file does not match the ", &
            "number of variables in the solution array: restart_nq=",i0, &
            ", current_nq=",i0,".")
  11 format("The number of global grid cells in the restart file does not ", &
            "match the current simulation: restart_n_global_cell=",i0, &
            ", current_n_global_cell=",i0,".")
  12 format(/,"  The number of solution points for the cells on CPU-",i0," ", &
              "does not match the",/,"   number of solution points for the ", &
              "corresponding cells in the restart file: ",/, &
              "        restart_n_solpts=",i0,/, &
              "        current_n_solpts=",i0)
  13 format("The number of solution points for cells in the restart file ", &
            "do not seem to match the corresponding cells in the current ", &
            "simulation.")
  !
  20 format("Unknown integer type when trying to read the integer section ", &
            "of the restart file header.")
  21 format("Unknown real type when trying to read the real section ", &
            "of the restart file header.")
  22 format("Unknown integer type when trying to read the ordered_geom data ", &
            "from the restart file.")
  23 format("Unknown real type when trying to read the solution data ", &
            "from the restart file.")
  24 format("Unknown integer type when trying to read the variable header ", &
            "for the solution data from the restart file.")
  !
  30 format(" CPU-",i0, &
            ": Error copying restart solution into solution array.",/, &
            "      restart cell: ",i0,":   geom = ",i0,", order = ",i0,/, &
            "      current cell: ",i0,":   geom = ",i0,", order = ",i0)
  31 format("The geometry and structure of the solution data in the restart ", &
            "file does not match the current simulation. There were ",i0," ", &
            "total errors encountered across all the processors when ", &
            "trying to copy the restart solution data into the current ", &
            "solution array.")
  !
  40 format ("The number of integer header items in the restart file",/, &
             "is greater than the size of the int_header array.",/, &
             "         n_header_int = ",i0,/, &
             "     size(int_header) = ",i0)
  41 format ("The number of real header items in the restart file",/, &
             "is greater than the size of the real_header array.",/, &
             "        n_header_real = ",i0,/, &
             "    size(real_header) = ",i0)
  !
  50 format("The number of chemical species in the restart file does not ", &
            "match the number of chemical species in the solution array: ", &
            "restart_ncs=",i0, &
            ", current_ncs=",i0, &
            ", inert_ncs=",i0,".")
  51 format("The number of variables in the restart file does not match the ", &
            "number of variables in the solution array: restart_nq=",i0, &
            ", current_nq=",i0, &
            ", inert_nq=",i0,".")
  !
end subroutine read_solution_file_mpi
!
!###############################################################################
!
subroutine read_sem_restart_mpi(num_sem_boxes,file_handle)
  !
  !.. Use Statements
  use random_mod,         only : iseed
  use synthetic_eddy_mod, only : sem_box
  !
  !.. Formal Arguments ..
  integer,         intent(in) :: num_sem_boxes
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  !
  !.. Local Scalars ..
  integer  :: restart_ip,restart_rp
  integer  :: restart_size,num_eddies
  integer  :: i,ierr,n,sem_idx
  real(wp) :: x,y,z
  !
  character(len=50) :: array_name
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: num_eddies_in_box(:)
  real(wp), allocatable :: eddy_data(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_sem_restart_mpi"
  !
continue
  !
  ! What needs to be read in:
  !   1. integer type of random seed array
  !   2. size of random seed array
  !   3. random seed array
  !   4. integer type of number of eddies in each box
  !   5. number of eddies in each box
  !   6. real type of eddy location data
  !   7. x,y,z eddy location data
  !         (one big array for all eddies across all SEM boxes)
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Read the integer type of the random seed array
  !
  call mpi_file_read_all(file_handle,restart_ip,1_int_mpi, &
                         MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Read the number of items in the random seed array
  !
  call mpi_file_read_all(file_handle,restart_size,1_int_mpi, &
                         MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  iseed = spread(0,dim=1,ncopies=restart_size) ! F2003 AUTO-REALLOC
  !
  ! Read the random seed array
  !
  select case (restart_ip)
    case (int8)  ; call read_rstrt_file_i1(file_handle,iseed)
    case (int16) ; call read_rstrt_file_i2(file_handle,iseed)
    case (int32) ; call read_rstrt_file_i4(file_handle,iseed)
    case (int64) ; call read_rstrt_file_i8(file_handle,iseed)
    case default
      write (error_message,1)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end select
  !
  ! #######
  ! ## 4 ##
  ! #######
  !
  ! Read the integer type of the array giving
  ! the number of eddies in each SEM box
  !
  call mpi_file_read_all(file_handle,restart_ip,1_int_mpi, &
                         MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
  !
  ! #######
  ! ## 5 ##
  ! #######
  !
  ! Read the number of eddies in each SEM box
  !
  num_eddies_in_box = spread(0,dim=1,ncopies=num_sem_boxes) ! F2003 AUTO-REALLOC
  !
  select case (restart_ip)
    case (int8)  ; call read_rstrt_file_i1(file_handle,num_eddies_in_box)
    case (int16) ; call read_rstrt_file_i2(file_handle,num_eddies_in_box)
    case (int32) ; call read_rstrt_file_i4(file_handle,num_eddies_in_box)
    case (int64) ; call read_rstrt_file_i8(file_handle,num_eddies_in_box)
    case default
      write (error_message,2)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end select
  !
  ! #######
  ! ## 6 ##
  ! #######
  !
  ! Read the real type of the eddy data
  !
  call mpi_file_read_all(file_handle,restart_rp,1_int_mpi, &
                         MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
  !
  ! #######
  ! ## 7 ##
  ! #######
  !
  ! Read the eddy data for all eddies across all boxes
  !
  num_eddies = sum(num_eddies_in_box)
  !
  ! F2003 AUTO-REALLOC
  eddy_data = reshape( spread(zero,dim=1,ncopies=9*num_eddies) , &
                       [9,num_eddies] )
  !
  if (restart_rp == wp) then
    call read_rstrt_file_wp(file_handle,eddy_data)
  else
    select case (restart_rp)
      case (real32)  ; call read_rstrt_file_r4(file_handle,eddy_data)
      case (real64)  ; call read_rstrt_file_r8(file_handle,eddy_data)
      case (real128) ; call read_rstrt_file_r16(file_handle,eddy_data)
      case default
        write (error_message,3)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end select
  end if
  !
  ! Then, loop over the SEM boxes and distribute the eddy data to the
  ! respective SEM box
  !
  if (allocated(sem_box)) then
    deallocate ( sem_box , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"sem_box",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( sem_box(1:num_sem_boxes) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sem_box",1,__LINE__,__FILE__,ierr,error_message)
  !
  i = 0
  !
  do sem_idx = 1,num_sem_boxes
    !
    num_eddies = num_eddies_in_box(sem_idx)
    !
    sem_box(sem_idx)%num_eddies = num_eddies
    !
    allocate ( sem_box(sem_idx)%eddy(1:num_eddies) , &
               stat=ierr , errmsg=error_message )
    write (array_name,10) sem_idx
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    do n = 1,num_eddies
      !
      i = i + 1
      !
      x = eddy_data(1,i)
      y = eddy_data(2,i)
      z = eddy_data(3,i)
      !
      sem_box(sem_idx)%eddy(n)%x = x
      sem_box(sem_idx)%eddy(n)%y = y
      sem_box(sem_idx)%eddy(n)%z = z
     !
     ! The polar components r and t are now computed in subroutine
     ! get_polar_coords_for_sem_fp which is called from initialize_sem.
     !
     !sem_box(sem_idx)%eddy(n)%r = hypot(y,z)
     !if (all([y,z] == zero)) then
     !  sem_box(sem_idx)%eddy(n)%t = zero
     !else
     !  sem_box(sem_idx)%eddy(n)%t = atan2(z,y)
     !end if
      !
      sem_box(sem_idx)%eddy(n)%x_intensity = eddy_data(4,i)
      sem_box(sem_idx)%eddy(n)%y_intensity = eddy_data(5,i)
      sem_box(sem_idx)%eddy(n)%z_intensity = eddy_data(6,i)
      !
      sem_box(sem_idx)%eddy(n)%x_frequency = eddy_data(7,i)
      sem_box(sem_idx)%eddy(n)%y_frequency = eddy_data(8,i)
      sem_box(sem_idx)%eddy(n)%z_frequency = eddy_data(9,i)
      !
    end do
    !
  end do
  !
  ! Deallocate local arrays before leaving
  !
  if (allocated(num_eddies_in_box)) then
    deallocate ( num_eddies_in_box , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"num_eddies_in_box",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(eddy_data)) then
    deallocate ( eddy_data , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"eddy_data",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Format Statements
  !
  1 format ("Unknown integer type when trying to read the random seed ", &
            "information from the SEM section of the restart file.")
  2 format ("Unknown integer type when trying to read the number of eddies ", &
            "in each SEM box from the SEM section of the restart file.")
  3 format ("Unknown real type when trying to read the eddy data ", &
            "from the SEM section of the restart file.")
  !
  10 format ("sem_box(",i0,")%eddy")
  !
end subroutine read_sem_restart_mpi
!
!###############################################################################
!
subroutine read_solution_file_serial(file_name,int_header,real_header, &
                                     solution,ordered_geom,var_header, &
                                     ignore_sem)
  !
  !.. Use Statements ..
  use geomvar_mod, only : ncell,cell,n_solpts
  use geomvar_mod, only : n_global_cell
  use ovar,        only : time_average_restart_files
  use thermo_mod,  only : inert_species_map
  !
  !.. Formal Arguments ..
  character(len=*),              intent(in   ) :: file_name
  integer(wip),     allocatable, intent(  out) :: int_header(:)
  real(wp),         allocatable, intent(  out) :: real_header(:)
  real(wp),         allocatable, intent(inout) :: solution(:,:)
  integer(ogp),     allocatable, intent(inout) :: ordered_geom(:)
  !
  !.. Optional Arguments ..
  integer(wip), optional, allocatable, intent(  out) :: var_header(:)
  logical(lk),  optional,              intent(in   ) :: ignore_sem
  !
  !.. Local Scalars ..
  integer :: io,l0,l2,m,mbeg,mm,n,np,n0,ierr
  integer :: sol_nq,ncs_diff,og
  integer :: inert_nq,inert_ncs
  integer :: restart_ip,restart_rp
  logical(lk) :: rst_form_matches_sol_form
  logical(lk) :: skip_sem_restart_data
  !
  integer :: solution_errors
  integer :: restart_nq,restart_ncs,restart_ncells
  integer :: n_header_int,n_header_real
  integer :: num_sem_boxes
  !
  integer :: inp_geom,inp_order
  integer :: this_cell,this_geom,this_order
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: rst_sol(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_solution_file_serial"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  skip_sem_restart_data = fals
  if (present(ignore_sem)) then
    skip_sem_restart_data = ignore_sem
  end if
  !
  num_sem_boxes = 0
  !
  solution_errors = 0
  !
  ! Open the restart file
  !
  open (newunit=io,file=trim(adjustl(file_name)),status="old",action="read", &
                   form="unformatted",access="stream",position="rewind", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,file_name,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read the header section of the restart file
  !
! ######################################
! ### Integer part of restart header ###
! ######################################
  !
  ! Read the integer type for the integer part of the restart header
  !
  read (io,iostat=ierr) restart_ip
  !
  ! Read the number of integer header items
  !
  read (io,iostat=ierr) n_header_int
  !
  ! Allocate the int_header array
  !
  if (allocated(int_header)) then
    deallocate ( int_header , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"int_header",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( int_header(1:n_header_int) , source=0_wip , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"int_header",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read the integer part of the restart header
  !
  select case (restart_ip)
    case (int32) ; call read_rstrt_file_serial_i4(io,int_header)
    case (int64) ; call read_rstrt_file_serial_i8(io,int_header)
    case (int8)  ; call read_rstrt_file_serial_i1(io,int_header)
    case (int16) ; call read_rstrt_file_serial_i2(io,int_header)
    case default
      write (error_message,20)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end select
  !
  ! Extract the information from the integer section of the header
  !
  ! Extract the number of flow variables saved at each solution point
  ! within the restart file
  !
  restart_nq = int( int_header(1) , kind=kind(restart_nq) )
  !
  if (allocated(solution)) then
    !
    sol_nq = size(solution,dim=2)
    !
    if (ncs > 0) then
      !
      inert_ncs = size(inert_species_map)
      inert_nq = (nq-nec+1) + inert_ncs
      restart_ncs = restart_nq - (nq-nec+1)
      !
      if (all(restart_ncs /= [ncs,inert_ncs])) then
        write (error_message,50) restart_ncs,ncs,inert_ncs
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      !
      if (all(restart_nq /= [sol_nq,inert_nq])) then
        write (error_message,51) restart_nq,sol_nq,inert_nq
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      !
    else if (restart_nq /= sol_nq) then
      !
      write (error_message,10) restart_nq,sol_nq
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      !
    end if
    !
  end if
  !
  ! Extracted the total number of grid cells contained within the restart file
  !
  restart_ncells = int( int_header(2) , kind=kind(restart_ncells) )
  !
  if (restart_ncells /= n_global_cell) then
    write (error_message,11) restart_ncells,n_global_cell
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Extract the format that the solution data is in within the restart file
  !
  restart_solution_form = int( int_header(3) , &
                               kind=kind(restart_solution_form) )
!%!
!%call check_restart_solution_form(restart_solution_form)
  !
  ! Check if there is restart information for
  ! the synthetic inflow turbulence generators
  !
  if (size(int_header) >= 5) then
    num_sem_boxes = max(0,int_header(5))
  end if
  !
! ###################################
! ### Real part of restart header ###
! ###################################
  !
  ! Read the real type for the real part of the restart header
  !
  read (io,iostat=ierr) restart_rp
  !
  ! Read the number of real header items
  !
  read (io,iostat=ierr) n_header_real
  !
  ! Allocate the real_header array
  !
  if (allocated(real_header)) then
    deallocate ( real_header , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"real_header",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( real_header(1:n_header_real) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"real_header",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read the real part of the restart header
  !
  if (restart_rp == wp) then
    call read_rstrt_file_serial_wp(io,real_header)
  else
    select case (restart_rp)
      case (real64)  ; call read_rstrt_file_serial_r8(io,real_header)
      case (real32)  ; call read_rstrt_file_serial_r4(io,real_header)
      case (real128) ; call read_rstrt_file_serial_r16(io,real_header)
      case default
        write (error_message,21)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end select
  end if
  !
! #############################################
! ### Ordered_Geom part of the restart file ###
! #############################################
  !
  ! Read the integer type for the ordered_geom array
  !
  read (io,iostat=ierr) restart_ip
  !
  if (.not. allocated(ordered_geom)) then
    allocate ( ordered_geom(1:restart_ncells) , source=0_ogp , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ordered_geom",1,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  ! Read the ordered_geom data so that each
  !
  select case (restart_ip)
    case (int32) ; call read_rstrt_file_serial_i4(io,ordered_geom)
    case (int64) ; call read_rstrt_file_serial_i8(io,ordered_geom)
    case (int8)  ; call read_rstrt_file_serial_i1(io,ordered_geom)
    case (int16) ; call read_rstrt_file_serial_i2(io,ordered_geom)
    case default
      write (error_message,22)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end select
  !
  ! Check to see if the form of the restart data matches the
  ! form of the solution data for the current simulation.
  !
  call check_restart_solution_form(restart_solution_form,ordered_geom, &
                                   rst_form_matches_sol_form)
  !
! ################################################
! ### Variable header part of the restart file ###
! ################################################
  !
  if (present(var_header)) then
    !
    ! Read the real type for the solution data
    !
    read (io,iostat=ierr) restart_ip
    !
    if (allocated(var_header)) then
      deallocate ( var_header , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"var_header",2,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
    !
    allocate ( var_header(1:restart_nq) , source=0_wip , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"var_header",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Read the variable header for the solution data
    !
    select case (restart_ip)
      case (int32) ; call read_rstrt_file_serial_i4(io,var_header)
      case (int64) ; call read_rstrt_file_serial_i8(io,var_header)
      case (int8)  ; call read_rstrt_file_serial_i1(io,var_header)
      case (int16) ; call read_rstrt_file_serial_i2(io,var_header)
      case default
        write (error_message,24)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end select
    !
    if (.not. allocated(solution)) then
      allocate ( solution(1:n_solpts,1:restart_nq) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"solution",1,__LINE__,__FILE__,ierr,error_message)
    end if
    !
  end if
  !
! #########################################
! ### Solution part of the restart file ###
! #########################################
  !
  ! Read the real type for the solution data
  !
  read (io,iostat=ierr) restart_rp
  !
  ! Compute the number of solution points in the restart file
  !
  num_pts = 0
  do this_cell = 1,ncell
    !
    og = int(ordered_geom(this_cell),kind=kind(og))
    !
    inp_order = og / 100
    inp_geom = og - 100*inp_order
    !
    num_pts = num_pts + Cell_Solpts(inp_geom,inp_order)
    !
  end do
!%!
!%! IF WE ARE TO MOVE THE POINT IN THE PRE-PROCESSING PHASE AT WHICH THE
!%! RESTART FILE IS READ, WE NEED TO MAKE SURE THAT N_SOLPTS HAS BEEN
!%! RECOMPUTED SO THAT IT IS LOCAL TO EACH PROCESSOR BECAUSE THIS IS NOT
!%! DONE WITHIN PARALLEL_MOD AFTER PARTITIONING THE GRID.
!%!
!%if (num_pts /= n_solpts) then
!%  write (iout,12) mypnum,num_pts,n_solpts
!%  write (error_message,13)
!%  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
!%end if
  !
  ! Allocate the rst_sol array to read in the solution restart data
  !
!%allocate ( rst_sol(1:restart_nq,1:n_solpts) , source=zero , &
  allocate ( rst_sol(1:restart_nq,1:num_pts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"rst_sol",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read in the solution restart data
  !
  if (restart_rp == wp) then
    call read_rstrt_file_serial_wp(io,rst_sol)
  else
    select case (restart_rp)
      case (real64)  ; call read_rstrt_file_serial_r8(io,rst_sol)
      case (real32)  ; call read_rstrt_file_serial_r4(io,rst_sol)
      case (real128) ; call read_rstrt_file_serial_r16(io,rst_sol)
      case default
        write (error_message,23)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end select
  end if
  !
  ! Read in the restart information for the synthetic
  ! inflow turbulence generators
  !
  if (.not. skip_sem_restart_data) then
    if (num_sem_boxes > 0) then
      call read_sem_restart_serial(num_sem_boxes,io)
    end if
  end if
  !
  ! Finally close the restart file
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,file_name,2,__LINE__,__FILE__,ierr,error_message)
  !
  !
  !
  if (rst_form_matches_sol_form) then
    !
    ! Copy the restart solution data over into the solution array
    !
    solution_errors = 0
    !
    l2 = 0
    do this_cell = 1,ncell
      !
      og = int(ordered_geom(this_cell),kind=kind(og))
      !
      inp_order = og / 100
      inp_geom = og - 100*inp_order
      !
      this_geom = cell(this_cell)%geom
      this_order = cell(this_cell)%order
      !
      ! NEED TO ADD SOMETHING HERE TO CONVERT THE SOLUTION IF THE SOLUTION
      ! FORM DOESNT MATCH THE CURRENT SIMULATION.  FOR NOW, JUST CREATE AN
      ! ERROR IF IT DOESNT MATCH THE RESTART FORM.
      !
      if ((inp_geom /= this_geom) .or. (inp_order /= this_order)) then
        write (iout,30) mypnum,this_cell, inp_geom, inp_order, &
                               this_cell,this_geom,this_order
        solution_errors = solution_errors + 1
      end if
      !
      n0 = cell(this_cell)%beg_sp - 1
      np = cell(this_cell)%end_sp - n0
      !
      ! Just in case inp_geom/=this_geom or inp_order/=this_order, we will use
      ! Cell_Solpts to compute the number of points for the input cell in order
      ! to increment l2. On the off chance either of those do not match, this
      ! will prevent any potential memory errors or segmentation faults that
      ! might occur if we assume  Cell_Solpts(inp_geom,inp_order) == np.
      !
      l0 = l2
      l2 = l2 + Cell_Solpts(inp_geom,inp_order)
      !
      ! Copy over the restart solutions for the continuity, momentum, and
      ! energy equations. Also copy over the restart solutions for the
      ! species equations if the number of chemical species in the restart
      ! file matches the number of chemical species in the current simulation.
      !
      ncs_diff = 0
      mbeg = 1
      if (ncs > 0 .and. restart_ncs /= ncs) then
        ncs_diff = ncs - restart_ncs
        mbeg = restart_ncs + 1
      end if
      !
      ! We need to flip the indices between the two arrays to match how the
      ! restart solution was read using the MPI derived datatype.
      !
      do n = 1,np
        do m = mbeg,restart_nq
          solution(n0+n,m+ncs_diff) = rst_sol(m,l0+n)
        end do
      end do
      !
      ! Copy over the restart solutions for the species equations if 
      ! the number of chemical species in the restart file does not
      ! match the number of chemical species in the current simulation.
      !
      if (ncs > 0 .and. restart_ncs /= ncs) then
        !
        if (restart_ncs == inert_ncs) then
          !
          ! We need to flip the indices between the two arrays to match how the
          ! restart solution was read using the MPI derived datatype.
          !
          do n = 1,np
            do mm = 1,restart_ncs
              m = inert_species_map(mm)
              solution(n0+n,m) = rst_sol(mm,l0+n)
            end do
          end do
          !
        else
          !
          write (error_message,50) restart_ncs,ncs,inert_ncs
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
          !
        end if
        !
      end if
      !
    end do
    !
  else
    !
    ! Interpolate/project the restart solution
    ! to the current solution points/order
    !
    call interp_restart_solution(restart_solution_form,rst_sol,solution, &
                                 ordered_geom=ordered_geom)
    !
  end if
  !
  ! Deallocate rst_sol since it is no longer needed
  !
  deallocate ( rst_sol , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"rst_sol",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Stop the simulation if there were any errors copying the restart
  ! solution data into solution.
  !
  if (solution_errors > 0) then
    write (error_message,31) solution_errors
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Reset num_pts to default values in case we ever need to come back here
  ! and read in a new restart file.
  !
  if (.not. time_average_restart_files) then
    num_pts = 0
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  10 format("The number of variables in the restart file does not match the ", &
            "number of variables in the solution array: restart_nq=",i0, &
            ", current_nq=",i0,".")
  11 format("The number of global grid cells in the restart file does not ", &
            "match the current simulation: restart_n_global_cell=",i0, &
            ", current_n_global_cell=",i0,".")
  12 format(/,"  The number of solution points for the cells on CPU-",i0," ", &
              "does not match the",/,"   number of solution points for the ", &
              "corresponding cells in the restart file: ",/, &
              "        restart_n_solpts=",i0,/, &
              "        current_n_solpts=",i0)
  13 format("The number of solution points for cells in the restart file ", &
            "do not seem to match the corresponding cells in the current ", &
            "simulation.")
  !
  20 format("Unknown integer type when trying to read the integer section ", &
            "of the restart file header.")
  21 format("Unknown real type when trying to read the real section ", &
            "of the restart file header.")
  22 format("Unknown integer type when trying to read the ordered_geom data ", &
            "from the restart file.")
  23 format("Unknown real type when trying to read the solution data ", &
            "from the restart file.")
  24 format("Unknown integer type when trying to read the variable header ", &
            "for the solution data from the restart file.")
  !
  30 format(" CPU-",i0, &
            ": Error copying restart solution into solution array.",/, &
            "      restart cell: ",i0,":   geom = ",i0,", order = ",i0,/, &
            "      current cell: ",i0,":   geom = ",i0,", order = ",i0)
  31 format("The geometry and structure of the solution data in the restart ", &
            "file does not match the current simulation. There were ",i0," ", &
            "total errors encountered across all the processors when ", &
            "trying to copy the restart solution data into the current ", &
            "solution array.")
  !
  50 format("The number of chemical species in the restart file does not ", &
            "match the number of chemical species in the solution array: ", &
            "restart_ncs=",i0, &
            ", current_ncs=",i0, &
            ", inert_ncs=",i0,".")
  51 format("The number of variables in the restart file does not match the ", &
            "number of variables in the solution array: restart_nq=",i0, &
            ", current_nq=",i0, &
            ", inert_nq=",i0,".")
  !
end subroutine read_solution_file_serial
!
!###############################################################################
!
subroutine read_sem_restart_serial(num_sem_boxes,io)
  !
  !.. Use Statements
  use random_mod,         only : iseed
  use synthetic_eddy_mod, only : sem_box
  !
  !.. Formal Arguments ..
  integer, intent(in) :: num_sem_boxes
  integer, intent(in) :: io
  !
  !.. Local Scalars ..
  integer  :: restart_ip,restart_rp
  integer  :: restart_size,num_eddies
  integer  :: i,ierr,n,sem_idx
  real(wp) :: x,y,z
  !
  character(len=50) :: array_name
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: num_eddies_in_box(:)
  real(wp), allocatable :: eddy_data(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_sem_restart_serial"
  !
continue
  !
  ! What needs to be read in:
  !   1. integer type of random seed array
  !   2. size of random seed array
  !   3. random seed array
  !   4. integer type of number of eddies in each box
  !   5. number of eddies in each box
  !   6. real type of eddy location data
  !   7. x,y,z eddy location data
  !         (one big array for all eddies across all SEM boxes)
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Read the integer type of the random seed array
  !
  read (io,iostat=ierr) restart_ip
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Read the number of items in the random seed array
  !
  read (io,iostat=ierr) restart_size
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  iseed = spread(0,dim=1,ncopies=restart_size) ! F2003 AUTO-REALLOC
  !
  ! Read the random seed array
  !
  select case (restart_ip)
    case (int8)  ; call read_rstrt_file_serial_i1(io,iseed)
    case (int16) ; call read_rstrt_file_serial_i2(io,iseed)
    case (int32) ; call read_rstrt_file_serial_i4(io,iseed)
    case (int64) ; call read_rstrt_file_serial_i8(io,iseed)
    case default
      write (error_message,1)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end select
  !
  ! #######
  ! ## 4 ##
  ! #######
  !
  ! Read the integer type of the array giving
  ! the number of eddies in each SEM box
  !
  read (io,iostat=ierr) restart_ip
  !
  ! #######
  ! ## 5 ##
  ! #######
  !
  ! Read the number of eddies in each SEM box
  !
  num_eddies_in_box = spread(0,dim=1,ncopies=num_sem_boxes) ! F2003 AUTO-REALLOC
  !
  select case (restart_ip)
    case (int8)  ; call read_rstrt_file_serial_i1(io,num_eddies_in_box)
    case (int16) ; call read_rstrt_file_serial_i2(io,num_eddies_in_box)
    case (int32) ; call read_rstrt_file_serial_i4(io,num_eddies_in_box)
    case (int64) ; call read_rstrt_file_serial_i8(io,num_eddies_in_box)
    case default
      write (error_message,2)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end select
  !
  ! #######
  ! ## 6 ##
  ! #######
  !
  ! Read the real type of the eddy data
  !
  read (io,iostat=ierr) restart_rp
  !
  ! #######
  ! ## 7 ##
  ! #######
  !
  ! Read the eddy data for all eddies across all boxes
  !
  num_eddies = sum(num_eddies_in_box)
  !
  ! F2003 AUTO-REALLOC
  eddy_data = reshape( spread(zero,dim=1,ncopies=9*num_eddies) , &
                       [9,num_eddies] )
  !
  if (restart_rp == wp) then
    call read_rstrt_file_serial_wp(io,eddy_data)
  else
    select case (restart_rp)
      case (real32)  ; call read_rstrt_file_serial_r4(io,eddy_data)
      case (real64)  ; call read_rstrt_file_serial_r8(io,eddy_data)
      case (real128) ; call read_rstrt_file_serial_r16(io,eddy_data)
      case default
        write (error_message,3)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end select
  end if
  !
  ! Then, loop over the SEM boxes and distribute the eddy data to the
  ! respective SEM box
  !
  if (allocated(sem_box)) then
    deallocate ( sem_box , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"sem_box",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( sem_box(1:num_sem_boxes) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sem_box",1,__LINE__,__FILE__,ierr,error_message)
  !
  i = 0
  !
  do sem_idx = 1,num_sem_boxes
    !
    num_eddies = num_eddies_in_box(sem_idx)
    !
    sem_box(sem_idx)%num_eddies = num_eddies
    !
    allocate ( sem_box(sem_idx)%eddy(1:num_eddies) , &
               stat=ierr , errmsg=error_message )
    write (array_name,10) sem_idx
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    do n = 1,num_eddies
      !
      i = i + 1
      !
      x = eddy_data(1,i)
      y = eddy_data(2,i)
      z = eddy_data(3,i)
      !
      sem_box(sem_idx)%eddy(n)%x = x
      sem_box(sem_idx)%eddy(n)%y = y
      sem_box(sem_idx)%eddy(n)%z = z
     !
     ! The polar components r and t are now computed in subroutine
     ! get_polar_coords_for_sem_fp which is called from initialize_sem.
     !
     !sem_box(sem_idx)%eddy(n)%r = hypot(y,z)
     !if (all([y,z] == zero)) then
     !  sem_box(sem_idx)%eddy(n)%t = zero
     !else
     !  sem_box(sem_idx)%eddy(n)%t = atan2(z,y)
     !end if
      !
      sem_box(sem_idx)%eddy(n)%x_intensity = eddy_data(4,i)
      sem_box(sem_idx)%eddy(n)%y_intensity = eddy_data(5,i)
      sem_box(sem_idx)%eddy(n)%z_intensity = eddy_data(6,i)
      !
      sem_box(sem_idx)%eddy(n)%x_frequency = eddy_data(7,i)
      sem_box(sem_idx)%eddy(n)%y_frequency = eddy_data(8,i)
      sem_box(sem_idx)%eddy(n)%z_frequency = eddy_data(9,i)
      !
    end do
    !
  end do
  !
  ! Deallocate local arrays before leaving
  !
  if (allocated(num_eddies_in_box)) then
    deallocate ( num_eddies_in_box , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"num_eddies_in_box",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(eddy_data)) then
    deallocate ( eddy_data , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"eddy_data",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Format Statements
  !
  1 format ("Unknown integer type when trying to read the random seed ", &
            "information from the SEM section of the restart file.")
  2 format ("Unknown integer type when trying to read the number of eddies ", &
            "in each SEM box from the SEM section of the restart file.")
  3 format ("Unknown real type when trying to read the eddy data ", &
            "from the SEM section of the restart file.")
  !
  10 format ("sem_box(",i0,")%eddy")
  !
end subroutine read_sem_restart_serial
!
!###############################################################################
!
subroutine write_restart_file(itnum)
  !
  !.. Use Statements ..
  use geomvar_mod,        only : n_global_cell
  use ovar,               only : itcur,time,d_t,rl2mx
  use ovar,               only : loc_solution_pts
  use ovar,               only : output_time_averaging
  use ovar,               only : restart_output_uses_host_roots
  use ovar,               only : max_restart_wall_time
  use ovar,               only : completely_disable_restart_output
  use flowvar_mod,        only : usp
  use synthetic_eddy_mod, only : sem_box
  use synthetic_eddy_mod, only : dump_sem_data
  !.. Use Statements for old restart files ..
  use geomvar_mod,        only : n_global_totpts
  !
  !.. Formal Arguments ..
  integer, intent(in) :: itnum
  !
  !.. Local Scalars ..
  integer :: ierr
  integer :: exitstat,cmndstat
  !
  logical(lk) :: involved_with_output
  logical(lk) :: restart_is_good_data
  !
  character(len=100) :: fname
  character(len=200) :: fullname
  character(len=100) :: statfile
  character(len=300) :: command
  character(len=300) :: cmndmsg
  !
  _MPI_COMM_TYPE_, save :: restart_comm
  !
  !.. Local Arrays ..
  integer(wip) :: int_header(1:5)
  real(wp)     :: real_header(1:3)
  !
  !.. Local Allocatable Arrays ..
  integer(ogp), allocatable :: ordered_geom(:)
  real(wp),     allocatable :: solution(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_restart_file"
  !
  logical(ldk), parameter :: wait_for_completion = .true.
  !
continue
  !
  if (completely_disable_restart_output) return
  !
  call _DEBUG_WRITE_TIMER_(entering_procedure,pname)
  !
  !
  !
  restart_is_good_data = true
  !
  ! Perform some initialization tasks if the restart module has not
  ! been initialized yet.
  !
  if (restart_needs_initialization) then
    !
    ! For now, write the nodal solution to the restart file. At some point we
    ! should switch so it writes the modal solution since that would give more
    ! flexibility regarding changes in the solution order between runs.
    !
    if (loc_solution_pts == Legendre_Gauss) then
      restart_solution_form = restart_gauss_nodal_solution
    else if (loc_solution_pts == Legendre_Gauss_Lobatto) then
      restart_solution_form = restart_lobatto_nodal_solution
    end if
    !
    ! Get the MPI derived types required to read/write the restart file
    !
    if (restart_output_uses_host_roots) then
      call get_host_output_mpi_datatypes
    else
      call get_output_mpi_datatypes
    end if
    !
#ifdef DEBUG_ON
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      write (iout,*) "After get_output_mpi_datatypes"
      flush (iout)
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
    !
    ! Unset the restart initialization flag since this has now been done
    !
    restart_needs_initialization = fals
    !
  end if
  !
  ! Get the name for the current restart file
  !
  if (itnum > 0) then
    write (fname,rstfname_fmt) itcur,""
    write (fullname,rstfile_fmt) itcur,""
    write (command,lnkcmnd_fmt) itcur,""
    write (statfile,statfile_fmt) itcur,""
  else if (itnum == 0) then
    write (fname,rstfname_fmt) itcur,".init"
    write (fullname,rstfile_fmt) itcur,".init"
    write (command,lnkcmnd_fmt) itcur,".init"
    write (statfile,statfile_fmt) itcur,".init"
  else if (itnum == -2) then
    write (fname,rstfname_fmt) itcur-1,".NaN"
    write (fullname,rstfile_fmt) itcur-1,".NaN"
    write (command,lnkcmnd_fmt) itcur-1,".NaN"
    write (statfile,statfile_fmt) itcur-1,".NaN"
    restart_is_good_data = fals
  else
    write (fname,rstfname_fmt) itcur,".done"
    write (fullname,rstfile_fmt) itcur,".done"
    write (command,lnkcmnd_fmt) itcur,".done"
    write (statfile,statfile_fmt) itcur,".done"
  end if
  !
  ! Create the cell information to write to the restart file
  !
  if (restart_output_uses_host_roots) then
    involved_with_output = i_am_host_root
    restart_comm = host_roots_comm
    call collect_restart_arrays(current_sol=usp, &
                                output_sol=solution, &
                                ordered_geom=ordered_geom)
  else
    involved_with_output = true
    restart_comm = MPI_COMM_WORLD
    call get_restart_arrays(current_sol=usp, &
                            output_sol=solution, &
                            ordered_geom=ordered_geom)
  end if
  !
  ! Store the size information for the simulation into the int_header array
  !
  int_header = [ int(nq,kind=wip) , &
                 int(n_global_cell,kind=wip) , &
                 int(restart_solution_form,kind=wip) , &
                 int(itcur,kind=wip) , &
                 size(sem_box,kind=wip) ]
  !
  ! Store the time data into the real_header array
  !
  real_header = [ time , d_t , rl2mx ]
  !
  ! Now write the actual restart file using all the
  ! MPI processes actually involved with the I/O
  !
  if (involved_with_output) then
    !
    if (ncpu > 1) then
      call write_solution_file_mpi(fullname,int_header,real_header, &
                                   restart_comm, &
                                   ordered_geom,ordered_geom_type, &
                                   solution,solution_type,sol_block_type, &
                                   max_restart_wall_time,statfile=statfile)
    else
      call write_solution_file_serial(fullname,int_header,real_header, &
                                      ordered_geom,solution, &
                                      max_restart_wall_time)
    end if
    !
  end if
  if (size(sem_box) > 0) then
    call dump_sem_data("write")
  end if
!!!#ifdef DEBUG_ON
#define USE_CREATE_SYMLINK
!!!#endif
#ifdef USE_CREATE_SYMLINK
  !
  ! Create a link to the most recent restart file for piggy-backing jobs
  !
  if (restart_is_good_data) then
    if (mypnum == glb_root) then
      call create_restart_link(fname)
    end if
  end if
#else
  !
  ! Create a link to the most recent restart file for piggy-backing jobs
  !
  if (restart_is_good_data) then
    call execute_command(command,wait_for_completion,exitstat, &
                         cmndstat,cmndmsg,pname)
  end if
#endif
  !
  ! Deallocate the solution array since its no longer needed
  !
  if (allocated(solution)) then
    deallocate ( solution , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"solution",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  ! Write the time averaged solution if we are time averaging
  !
  if (output_time_averaging .and. restart_is_good_data) then
    call write_time_averaged_solution(itnum,ordered_geom,int_header)
  end if
  !
  ! Deallocate ordered_geom since it is no longer needed
  !
  if (allocated(ordered_geom)) then
    deallocate ( ordered_geom , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ordered_geom",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  call _DEBUG_WRITE_TIMER_(leaving_procedure,pname)
  !
end subroutine write_restart_file
!
!###############################################################################
!
subroutine write_solution_file_mpi(file_name,int_header,real_header,file_comm, &
                                   ordered_geom,og_type,solution,sol_type, &
                                   sol_blktyp,max_wall_time,var_header,statfile)
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_global_solpts
  use ovar,        only : lustre_stripe_factor
  use ovar,        only : lustre_stripe_unit
  use ovar,        only : restart_output_uses_host_roots
  use ovar,        only : profile_io_restart
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: file_name
  integer,          intent(in) :: int_header(:)
  real(wp),         intent(in) :: real_header(:)
  _MPI_COMM_TYPE_,  intent(in) :: file_comm
  integer(ogp),     intent(in) :: ordered_geom(:)
  _MPI_DATA_TYPE_,  intent(in) :: og_type
  real(wp),         intent(in) :: solution(:,:)
  _MPI_DATA_TYPE_,  intent(in) :: sol_type
  _MPI_DATA_TYPE_,  intent(in) :: sol_blktyp
  real(wp),      intent(inout) :: max_wall_time
  !
  !.. Optional Arguments ..
  integer,          optional, intent(in) :: var_header(:)
  character(len=*), optional, intent(in) :: statfile
  !
  !.. Local Scalars ..
  integer      :: ios
  integer      :: num_sem_boxes
  integer(wip) :: disp,isize
  !
  _MPI_FILE_TYPE_ :: file_handle
  _MPI_INFO_TYPE_ :: file_info
  integer(MPI_OFFSET_KIND) :: file_disp
  integer(MPI_OFFSET_KIND) :: file_offset
  integer(MPI_OFFSET_KIND) :: int_bytes
  integer(MPI_OFFSET_KIND) :: real_bytes
  integer(kind(MPI_MODE_WRONLY)) :: file_mode
  !
  real(wp) :: beg_wall_time
  real(wp) :: end_wall_time
  real(wp) :: restart_wall_time
  !
  character(len=200) :: strng
  integer(int_mpi)   :: saved_mpierr
  !
  !.. Local Arrays ..
  _MPI_STATUS_TYPE_ :: istat(1:_MPI_STATUS_SIZE_)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_solution_file_mpi"
  !
continue
  !
  call _DEBUG_WRITE_TIMER_(entering_procedure,pname)
  !
  ! Check if there is restart information for
  ! the synthetic inflow turbulence generators
  !
  num_sem_boxes = 0
  if (size(int_header) >= 5) then
    num_sem_boxes = int_header(5)
  end if
  !
  call mpi_barrier(file_comm,mpierr)
  if (mypnum == glb_root) then
    write (iout,1) pname
    flush (iout)
    beg_wall_time = mpi_wtime()
  end if
  call mpi_barrier(file_comm,mpierr)
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            open_file=true, &
                            stat_strng="init stat, before mpi_info")
  !
  ! Initialize file_info to MPI_INFO_NULL
  !
  file_info = MPI_INFO_NULL
  !
#ifdef PBS_ENV
#ifndef IGNORE_MPI_INFO
  !
  ! Create the MPI info object for this restart file
  !
  call mpi_info_create(file_info,mpierr)
  strng = "after mpi_info_create, " // &
          "before mpi_info_set(access_style)"
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            stat_strng=strng,mpi_error=mpierr)
  !
  ! Now give MPI some I/O hints through the info object
  !
  !-----------------------  EXAMPLE  ------------------------
  !call mpi_info_set(file_info,"info key","key value",mpierr)
  !-----------------------  EXAMPLE  ------------------------
  !
  ! Specify the nubmer of I/O devices in the system
  !
 !call mpi_info_set(file_info,"num_io_nodes","1",mpierr)
  !
  call mpi_info_set(file_info,"access_style","write_once",mpierr)
  strng = "after mpi_info_set(access_style), " // &
          "before mpi_info_set(romio_cb_write)"
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            stat_strng=strng,mpi_error=mpierr)
  !
  ! Enable collective buffering for writes
  !
  call mpi_info_set(file_info,"romio_cb_write","enable",mpierr)
  strng = "after mpi_info_set(romio_cb_write), " // &
          "before mpi_info_set(striping_factor)"
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            stat_strng=strng,mpi_error=mpierr)
  !
  ! Specify the number of I/O devices that the
  ! restart file should be striped across
  !
  call mpi_info_set(file_info,"striping_factor", &
                    trim(adjustl(lustre_stripe_factor)),mpierr)
  strng = "after mpi_info_set(striping_factor), " // &
          "before mpi_info_set(striping_unit)"
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            stat_strng=strng,mpi_error=mpierr)
  !
  ! Specifiy the suggested striping unit to be used for this file. The
  ! striping unit is the amount of consecutive data assigned to one I/O
  ! device before progressing to the next device, when striping across
  ! a number of devices. It is expressed in bytes.
  !
  call mpi_info_set(file_info,"striping_unit", &
                    trim(adjustl(lustre_stripe_unit)),mpierr)
  strng = "after mpi_info_set(striping_unit)"
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            stat_strng=strng,mpi_error=mpierr)
  !
#endif
#endif
#ifdef DEBUG_ON
  call mpi_barrier(file_comm,mpierr)
  if (mypnum == glb_root) then
    write (iout,*) "Before mpi_file_open"
    flush (iout)
  end if
  call mpi_barrier(file_comm,mpierr)
#endif
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            stat_strng="before mpi_file_open")
  !
  ! Open the restart file for MPI I/O
  !
  file_mode = MPI_MODE_WRONLY + MPI_MODE_CREATE
  call mpi_file_open(file_comm,trim(adjustl(file_name)), &
                     file_mode,file_info,file_handle,mpierr)
  saved_mpierr = mpierr
#ifdef DEBUG_ON
  call mpi_barrier(file_comm,mpierr)
  if (mypnum == glb_root) then
    write (iout,*) "After mpi_file_open"
    flush (iout)
  end if
  call mpi_barrier(file_comm,mpierr)
#endif
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            stat_strng="after mpi_file_open", &
                            mpi_error=saved_mpierr)
  !
  ! Have the root processor write the header for the restart file
  !
  if (mypnum == glb_root) then
    !
    ! Write the integer type for the size information for the simulation
    !
    call mpi_file_write(file_handle,wip, &
                        1_int_mpi,MPI_INTEGER,istat(1),mpierr)
    !
    ! Write the number of size data items
    !
    isize = size(int_header,kind=wip)
    call mpi_file_write(file_handle,isize, &
                        1_int_mpi,MPI_INTEGER,istat(1),mpierr)
    !
    ! Write the size information for the simulation to the restart file
    !
    call mpi_file_write(file_handle,int_header, &
                        size(int_header,kind=int_mpi), &
                        mpi_inttyp,istat(1),mpierr)
    !
    ! Write the real type for the time data
    !
    call mpi_file_write(file_handle,wp, &
                        1_int_mpi,MPI_INTEGER,istat(1),mpierr)
    !
    ! Write the number of time data items
    !
    isize = size(real_header,kind=wip)
    call mpi_file_write(file_handle,isize, &
                        1_int_mpi,MPI_INTEGER,istat(1),mpierr)
    !
    ! Write the time data to the restart file
    !
    call mpi_file_write(file_handle,real_header, &
                        size(real_header,kind=int_mpi), &
                        mpi_flttyp,istat(1),mpierr)
    !
    ! Write the integer type for the ordered_geom array
    !
    call mpi_file_write(file_handle,i4, &
                        1_int_mpi,MPI_INTEGER,istat(1),mpierr)
    !
    ! Get the updated file displacement for the root processor
    !
    call mpi_file_get_position(file_handle,file_offset,mpierr)
    call mpi_file_get_byte_offset(file_handle,file_offset,file_disp,mpierr)
    !
    ! Convert the MPI_OFFSET_KIND file displacement to the working integer
    ! precision so that it can be broadcast to all processors
    !
    disp = int( file_disp , kind=wip )
    !
  end if
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            stat_strng="after root writes restart header")
  !
  ! Broadcast the current file displacement to all processors
  !
  call mpi_bcast(disp,1_int_mpi,mpi_inttyp,glb_root,file_comm,mpierr)
  !
  ! Convert the file displacement to MPI_OFFSET_KIND
  !
  file_disp = int( disp , kind=MPI_OFFSET_KIND )
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                     stat_strng="after bcast disp from restart header", &
                     mpi_error=mpierr)
  !
  ! Set the file view to write the ordered_geom array to the restart file
  !
  call mpi_file_set_view(file_handle,file_disp,MPI_INTEGER, &
                         og_type, trim(adjustl(file_datarep)), &
                         file_info, mpierr)
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                     stat_strng="after setting view for ordered_geom", &
                     mpi_error=mpierr)
  !
  ! Write the ordered_geom array to the restart file
  !
  call mpi_file_write_all(file_handle,ordered_geom, &
                          size(ordered_geom,kind=int_mpi), &
                          MPI_INTEGER,istat(1),mpierr)
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            mpistat=istat, &
                            stat_strng="after writing ordered_geom", &
                            mpi_error=mpierr)
  !
  ! Update the file displacement since we know the total number of cells
  ! and that the elements of the array ordered_geom are 4-byte integers.
  !
  int_bytes = storage_size(ordered_geom,kind=MPI_OFFSET_KIND) / &
              int(8,kind=MPI_OFFSET_KIND)
  file_disp = file_disp + int_bytes*int(int_header(2),kind=MPI_OFFSET_KIND)
  !
  ! Update the file view to after all the ordered_geom data
  !
  call mpi_file_set_view(file_handle,file_disp,MPI_INTEGER, &
                         MPI_INTEGER, trim(adjustl(file_datarep)), &
                         file_info, mpierr)
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            stat_strng="after setting view for wp", &
                            mpi_error=mpierr)
  !
  ! Have the root processor write the variable header and
  ! the byte size for the solution data to the restart file
  !
  if (mypnum == glb_root) then
    !
    if (present(var_header)) then
      !
      ! Write the integer type for the size information for the simulation
      !
      call mpi_file_write(file_handle,wip, &
                          1_int_mpi,MPI_INTEGER,istat(1),mpierr)
      !
      ! Write the variable header information for
      ! the solution data to the restart file
      !
      call mpi_file_write(file_handle,var_header, &
                          size(var_header,kind=int_mpi), &
                          MPI_INTEGER,istat(1),mpierr)
      !
    end if
    !
    ! Write the kind value for real working-precision
    !
    call mpi_file_write(file_handle,wp,1_int_mpi, &
                        MPI_INTEGER,istat(1),mpierr)
    !
    ! Get the updated file displacement for the root processor
    !
    call mpi_file_get_position(file_handle,file_offset,mpierr)
    call mpi_file_get_byte_offset(file_handle,file_offset,file_disp,mpierr)
    !
    ! Convert the MPI_OFFSET_KIND file displacement to the working integer
    ! precision so that it can be broadcast to all processors
    !
    disp = int( file_disp , kind=wip )
    !
  end if
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            stat_strng="after root writes wp")
  !
  ! Broadcast the current file displacement to all processors
  !
  call mpi_bcast(disp,1_int_mpi,mpi_inttyp,glb_root,file_comm,mpierr)
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            stat_strng="after bcast disp from wp", &
                            mpi_error=mpierr)
  !
  ! Convert the file displacement to MPI_OFFSET_KIND
  !
  file_disp = int( disp , kind=MPI_OFFSET_KIND )
  !
  ! Set the file view to write the solution data to the restart file
  !
  call mpi_file_set_view(file_handle,file_disp,sol_blktyp, &
                         sol_type, trim(adjustl(file_datarep)), &
                         file_info, mpierr)
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            stat_strng="after setting view for solution", &
                            mpi_error=mpierr)
  !
  ! Write the solution data to the restart file
  !
  call mpi_file_write_all(file_handle,solution, &
                          size(solution,dim=2,kind=int_mpi), &
                          sol_blktyp,istat(1),mpierr)
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            mpistat=istat, &
                            stat_strng="after writing solution", &
                            mpi_error=mpierr)
  !
  ! Update the file displacement since we know the total number of
  ! solution points, the number of variables per solution point, and
  ! the byte-size of for each element of the solution array.
  !
  real_bytes = storage_size(solution,kind=MPI_OFFSET_KIND) / &
               int(8,kind=MPI_OFFSET_KIND)
  file_disp = file_disp + real_bytes * &
                          int(int_header(1),kind=MPI_OFFSET_KIND) * &
                          int(n_global_solpts,kind=MPI_OFFSET_KIND)
  !
  ! Update the file view to after all the solution data
  !
  call mpi_file_set_view(file_handle,file_disp,MPI_INTEGER, &
                         MPI_INTEGER, trim(adjustl(file_datarep)), &
                         file_info, mpierr)
  !
  ! Write the synthetic inflow turbulence generator data to the restart file
  !
  if (num_sem_boxes > 0) then
    call write_sem_restart_mpi(file_handle)
  end if
  !
  ! Finally close the restart file
  !
  call mpi_file_close(file_handle,mpierr)
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            close_file=true, &
                            stat_strng="after mpi_file_close", &
                            mpi_error=mpierr)
  !
  ! Output the wall time taken to write the restart solution file if
  ! the input parameter profile_io_restart is true
  !
  call mpi_barrier(file_comm,mpierr)
  if (mypnum == glb_root) then
    end_wall_time = mpi_wtime()
    restart_wall_time = end_wall_time - beg_wall_time
    max_wall_time = max(max_wall_time,restart_wall_time)
    if (profile_io_restart) then
      write (iout,2) pname,trim(adjustl(file_name)), &
                     restart_wall_time, &
                     trim(adjustl(merge("ONLY HOST ROOTS  ", &
                                        "ALL MPI PROCESSES", &
                                        restart_output_uses_host_roots)))
    else
      write (iout,2) pname
    end if
    flush (iout)
  end if
  call mpi_barrier(file_comm,mpierr)
  !
  call _DEBUG_WRITE_TIMER_(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," Entering ",a)
  2 format (" Leaving ",a,/,:, &
          5x,"Wall time taken to write the file '",a, &
             "' = ",f12.6," seconds",/, &
         10x,"(USING ",a," FOR MPI-I/O)",/)
  !
end subroutine write_solution_file_mpi
!
!###############################################################################
!
subroutine write_sem_restart_mpi(file_handle)
  !
  !.. Use Statements
  use random_mod,         only : get_random_seed,iseed
  use synthetic_eddy_mod, only : sem_box
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  !
  !.. Local Scalars ..
  integer(int_mpi) :: restart_ip,restart_rp
  integer(int_mpi) :: restart_size
  !
  integer :: i,ierr,n,num_eddies,num_sem_boxes,sem_idx
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: num_eddies_in_box(:)
  real(wp), allocatable :: eddy_data(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_sem_restart_mpi"
  !
continue
  !
  ! Get the current random seed to save to the restart file
  !
  call get_random_seed()
  !
  ! Collect all the eddy data into a single array
  ! NOTE: ONLY THE GLOBAL ROOT PROCESS NEEDS TO DO THE REST OF THIS SUBROUTINE
  !
  if (mypnum == glb_root) then
    !
    num_sem_boxes = size(sem_box)
    !
    num_eddies_in_box = sem_box(:)%num_eddies ! F2003 AUTO-REALLOC
    !
    num_eddies = sum(num_eddies_in_box)
    !
    ! Collect all the eddy data into a single array
    !
    ! F2003 AUTO-REALLOC
    eddy_data = reshape( spread(zero,dim=1,ncopies=9*num_eddies) , &
                         [9,num_eddies] )
    !
    i = 0
    !
    do sem_idx = 1,num_sem_boxes
      !
      do n = 1,sem_box(sem_idx)%num_eddies
        !
        i = i + 1
        !
        eddy_data(1,i) = sem_box(sem_idx)%eddy(n)%x
        eddy_data(2,i) = sem_box(sem_idx)%eddy(n)%y
        eddy_data(3,i) = sem_box(sem_idx)%eddy(n)%z
        !
        eddy_data(4,i) = sem_box(sem_idx)%eddy(n)%x_intensity
        eddy_data(5,i) = sem_box(sem_idx)%eddy(n)%y_intensity
        eddy_data(6,i) = sem_box(sem_idx)%eddy(n)%z_intensity
        !
        eddy_data(7,i) = sem_box(sem_idx)%eddy(n)%x_frequency
        eddy_data(8,i) = sem_box(sem_idx)%eddy(n)%y_frequency
        eddy_data(9,i) = sem_box(sem_idx)%eddy(n)%z_frequency
        !
      end do
      !
    end do
    !
    ! What needs to be written to the restart file:
    !   1. integer type of random seed array
    !   2. size of random seed array
    !   3. random seed array
    !   4. integer type of number of eddies in each box
    !   5. number of eddies in each box
    !   6. real type of eddy location data
    !   7. x,y,z eddy location data
    !         (one big array for all eddies across all SEM boxes)
    !
    ! #######
    ! ## 1 ##
    ! #######
    !
    ! Write the integer type of the random seed array
    !
    restart_ip = kind(iseed)
    call mpi_file_write(file_handle,restart_ip,1_int_mpi, &
                        MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
    !
    ! #######
    ! ## 2 ##
    ! #######
    !
    ! Write the number of items in the random seed array
    !
    restart_size = size(iseed,kind=int_mpi)
    call mpi_file_write(file_handle,restart_size,1_int_mpi, &
                        MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
    !
    ! #######
    ! ## 3 ##
    ! #######
    !
    ! Write the random seed array
    !
    call mpi_file_write(file_handle,iseed,size(iseed,kind=int_mpi), &
                        mpi_inttyp,MPI_STATUS_IGNORE,mpierr)
    !
    ! #######
    ! ## 4 ##
    ! #######
    !
    ! Write the integer type of the array giving
    ! the number of eddies in each SEM box
    !
    restart_ip = kind(num_eddies_in_box)
    call mpi_file_write(file_handle,restart_ip,1_int_mpi, &
                        MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
    !
    ! #######
    ! ## 5 ##
    ! #######
    !
    ! Write the number of eddies in each SEM box
    !
    call mpi_file_write(file_handle,num_eddies_in_box, &
                        size(num_eddies_in_box,kind=int_mpi), &
                        mpi_inttyp,MPI_STATUS_IGNORE,mpierr)
    !
    ! Deallocate num_eddies_in_box since it is no longer needed
    !
    if (allocated(num_eddies_in_box)) then
      deallocate ( num_eddies_in_box , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"num_eddies_in_box",2,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
    !
    ! #######
    ! ## 6 ##
    ! #######
    !
    !
    ! Write the real type of the eddy data
    !
    restart_rp = kind(eddy_data)
    call mpi_file_write(file_handle,restart_rp,1_int_mpi, &
                        MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
    !
    ! #######
    ! ## 7 ##
    ! #######
    !
    ! Write the eddy data for all eddies across all boxes
    !
    call mpi_file_write(file_handle,eddy_data,size(eddy_data,kind=int_mpi), &
                        mpi_flttyp,MPI_STATUS_IGNORE,mpierr)
    !
    ! Deallocate eddy_data before leaving
    !
    if (allocated(eddy_data)) then
      deallocate ( eddy_data , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"eddy_data",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
  end if
  !
end subroutine write_sem_restart_mpi
!
!###############################################################################
!
subroutine write_solution_file_serial(file_name,int_header,real_header, &
                                      ordered_geom,solution,max_wall_time, &
                                      var_header)
  !
  !.. Use Statements ..
  use ovar, only : profile_io_restart
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: file_name
  integer,          intent(in) :: int_header(:)
  real(wp),         intent(in) :: real_header(:)
  integer(ogp),     intent(in) :: ordered_geom(:)
  real(wp),         intent(in) :: solution(:,:)
  real(wp),      intent(inout) :: max_wall_time
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: var_header(:)
  !
  !.. Local Scalars ..
  integer :: i,n,io,ierr
  integer :: num_sem_boxes
  !
  real(wp) :: beg_wall_time
  real(wp) :: end_wall_time
  real(wp) :: restart_wall_time
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_solution_file_serial"
  !
continue
  !
  call global_debug_timer(entering_procedure,pname)
  !
  ! Check if there is restart information for
  ! the synthetic inflow turbulence generators
  !
  num_sem_boxes = 0
  if (size(int_header) >= 5) then
    num_sem_boxes = int_header(5)
  end if
  !
  write (iout,1) pname
  flush (iout)
  beg_wall_time = mpi_wtime()
  !
  ! Open the restart file
  !
  open (newunit=io,file=trim(adjustl(file_name)), &
                   status="replace",action="write", &
                   form="unformatted",access="stream",position="rewind", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,file_name,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Write the header for the restart file
  !
  ! Write the integer type for the size information for the simulation
  !
  write (io,iostat=ierr) wip
  !
  ! Write the number of size data items
  !
  write (io,iostat=ierr) size(int_header,kind=wip)
  !
  ! Write the size information for the simulation to the restart file
  !
  write (io,iostat=ierr) (int_header(i), i=1,size(int_header))
  !
  ! Write the real type for the time data
  !
  write (io,iostat=ierr) wp
  !
  ! Write the number of time data items
  !
  write (io,iostat=ierr) size(real_header,kind=wip)
  !
  ! Write the time data to the restart file
  !
  write (io,iostat=ierr) (real_header(i), i=1,size(real_header))
  !
  ! Write the integer type for the ordered_geom array
  !
  write (io,iostat=ierr) i4
  !
  ! Write the ordered_geom array to the restart file
  !
  write (io,iostat=ierr) (ordered_geom(i), i=1,size(ordered_geom))
  !
  ! Write the variable header for the solution data to the restart file
  !
  if (present(var_header)) then
    !
    ! Write the integer type for the size information for the simulation
    !
    write (io,iostat=ierr) wip
    !
    ! Write the variable header information for
    ! the solution data to the restart file
    !
    write (io,iostat=ierr) (var_header(i), i=1,size(var_header))
    !
  end if
  !
  ! Write the kind value for real working-precision
  !
  write (io,iostat=ierr) wp
  !
  ! Write the solution data to the restart file
  !
  write (io,iostat=ierr) ((solution(i,n), i=1,size(solution,dim=1)), &
                                          n=1,size(solution,dim=2))
  !
  ! Write the synthetic inflow turbulence generator data to the restart file
  !
  if (num_sem_boxes > 0) then
    call write_sem_restart_serial(io)
  end if
  !
  ! Finally close the restart file
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,file_name,2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Output the wall time taken to write the restart solution file if
  ! the input parameter profile_io_restart is true
  !
  end_wall_time = mpi_wtime()
  restart_wall_time = end_wall_time - beg_wall_time
  max_wall_time = max(max_wall_time,restart_wall_time)
  if (profile_io_restart) then
    write (iout,2) pname,trim(adjustl(file_name)), restart_wall_time
  else
    write (iout,2) pname
  end if
  flush (iout)
  !
  call global_debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," Entering ",a)
  2 format (" Leaving ",a,/,:, &
          5x,"Wall time taken to write the file '",a, &
             "' = ",f12.6," seconds",/, &
         10x,"(USING SERIAL-I/O)",/)
  !
end subroutine write_solution_file_serial
!
!###############################################################################
!
subroutine write_sem_restart_serial(io)
  !
  !.. Use Statements
  use random_mod,         only : get_random_seed,iseed
  use synthetic_eddy_mod, only : sem_box
  !
  !.. Formal Arguments ..
  integer, intent(in) :: io
  !
  !.. Local Scalars ..
  integer :: i,ierr,n,num_eddies,num_sem_boxes,sem_idx
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: num_eddies_in_box(:)
  real(wp), allocatable :: eddy_data(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_sem_restart_serial"
  !
continue
  !
  ! Get the current random seed to save to the restart file
  !
  call get_random_seed()
  !
  ! Collect all the eddy data into a single array
  ! NOTE: ONLY THE GLOBAL ROOT PROCESS NEEDS TO DO THE REST OF THIS SUBROUTINE
  !
  if (mypnum == glb_root) then
    !
    num_sem_boxes = size(sem_box)
    !
    num_eddies_in_box = sem_box(:)%num_eddies ! F2003 AUTO-REALLOC
    !
    num_eddies = sum(num_eddies_in_box)
    !
    ! Collect all the eddy data into a single array
    !
    ! F2003 AUTO-REALLOC
    eddy_data = reshape( spread(zero,dim=1,ncopies=9*num_eddies) , &
                         [9,num_eddies] )
    !
    i = 0
    !
    do sem_idx = 1,num_sem_boxes
      !
      do n = 1,sem_box(sem_idx)%num_eddies
        !
        i = i + 1
        !
        eddy_data(1,i) = sem_box(sem_idx)%eddy(n)%x
        eddy_data(2,i) = sem_box(sem_idx)%eddy(n)%y
        eddy_data(3,i) = sem_box(sem_idx)%eddy(n)%z
        !
        eddy_data(4,i) = sem_box(sem_idx)%eddy(n)%x_intensity
        eddy_data(5,i) = sem_box(sem_idx)%eddy(n)%y_intensity
        eddy_data(6,i) = sem_box(sem_idx)%eddy(n)%z_intensity
        !
        eddy_data(7,i) = sem_box(sem_idx)%eddy(n)%x_frequency
        eddy_data(8,i) = sem_box(sem_idx)%eddy(n)%y_frequency
        eddy_data(9,i) = sem_box(sem_idx)%eddy(n)%z_frequency
        !
      end do
      !
    end do
    !
    ! What needs to be written to the restart file:
    !   1. integer type of random seed array
    !   2. size of random seed array
    !   3. random seed array
    !   4. integer type of number of eddies in each box
    !   5. number of eddies in each box
    !   6. real type of eddy location data
    !   7. x,y,z eddy location data
    !         (one big array for all eddies across all SEM boxes)
    !
    ! #######
    ! ## 1 ##
    ! #######
    !
    ! Write the integer type of the random seed array
    !
    write (io,iostat=ierr) kind(iseed)
    !
    ! #######
    ! ## 2 ##
    ! #######
    !
    ! Write the number of items in the random seed array
    !
    write (io,iostat=ierr) size(iseed,kind=int_mpi)
    !
    ! #######
    ! ## 3 ##
    ! #######
    !
    ! Write the random seed array
    !
    write (io,iostat=ierr) (iseed(i), i=1,size(iseed))
    !
    ! #######
    ! ## 4 ##
    ! #######
    !
    ! Write the integer type of the array giving
    ! the number of eddies in each SEM box
    !
    write (io,iostat=ierr) kind(num_eddies_in_box)
    !
    ! #######
    ! ## 5 ##
    ! #######
    !
    ! Write the number of eddies in each SEM box
    !
    write (io,iostat=ierr) (num_eddies_in_box(i), i=1,size(num_eddies_in_box))
    !
    ! Deallocate num_eddies_in_box since it is no longer needed
    !
    if (allocated(num_eddies_in_box)) then
      deallocate ( num_eddies_in_box , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"num_eddies_in_box",2,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
    !
    ! #######
    ! ## 6 ##
    ! #######
    !
    !
    ! Write the real type of the eddy data
    !
    write (io,iostat=ierr) kind(eddy_data)
    !
    ! #######
    ! ## 7 ##
    ! #######
    !
    ! Write the eddy data for all eddies across all boxes
    !
    write (io,iostat=ierr) ((eddy_data(i,n), i=1,size(eddy_data,dim=1)), &
                                             n=1,size(eddy_data,dim=2))
    !
    ! Deallocate eddy_data before leaving
    !
    if (allocated(eddy_data)) then
      deallocate ( eddy_data , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"eddy_data",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
  end if
  !
end subroutine write_sem_restart_serial
!
!###############################################################################
!
subroutine write_time_averaged_solution(itnum,ordered_geom,int_header)
  !
  !.. Use Statements ..
  use ovar,        only : time,prev_ave_time,ave_start_time,itcur,rl2mx
  use ovar,        only : restart_output_uses_host_roots
  use ovar,        only : max_time_ave_rst_wall_time
  use ovar,        only : num_saved_time_ave_restarts
  use ovar,        only : is_at_restart_frequency
  use flowvar_mod, only : uavesp,time_ave_variables
  !
  !.. Formal Arguments ..
  integer,         intent(in) :: itnum
  integer(ogp),    intent(in) :: ordered_geom(:)
  integer(wip), intent(inout) :: int_header(:)
  !
  !.. Local Scalars ..
  integer :: n,ierr
  integer :: exitstat,cmndstat
  !
  logical(lk) :: involved_with_output
  !
  character(len=100) :: fname
  character(len=200) :: fullname
  character(len=100) :: statfile
  character(len=300) :: link_current_command
  character(len=300) :: link_iter_command
  character(len=300) :: rm_iter_link_command
  character(len=300) :: cmndmsg
  !
  _MPI_COMM_TYPE_ :: timeave_comm
  !
  !.. Local Arrays ..
  real(wp) :: real_header(1:4)
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: var_header(:)
  real(wp), allocatable :: solution(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_time_averaged_solution"
  !
  logical(ldk), parameter :: wait_for_completion = .true.
  !
continue
  !
  call _DEBUG_WRITE_TIMER_(entering_procedure,pname)
  !
  if (time_ave_needs_initialization) then
    !
    ! Get the MPI derived types required to
    ! read/write the time-averaging restart file
    !
    if (restart_output_uses_host_roots) then
      call get_host_time_ave_mpi_datatypes
    else
      call get_time_ave_mpi_datatypes
    end if
    !
    ! F2003 AUTO-REALLOC
    output_iterations = init_array(num_saved_time_ave_restarts,0)
    !
    ! Unset the restart initialization flag since this has now been done
    !
    time_ave_needs_initialization = fals
    !
  end if
  !
  ! Get the name for the current time averaged solution file
  !
  ! RDIR = trim(adjustl(restart_dir))
  ! FID = final_iter_digits = # of digits for itrst + num_timesteps
  ! TAD = time_ave_digits   = # of digits for averaging_file_switch
  !
  ! tafname_fmt     =      "('out.',i0.TAD,a,'.ave')"
  ! tafile_fmt      = "('RDIR/out.',i0.TAD,a,'.ave')"
  ! tadonefname_fmt =      "('out.',i0.FID,a,'.ave')"
  ! tadone_fmt      = "('RDIR/out.',i0.FID,a,'.ave')"
  !
  if (itnum > 0) then
    if (is_at_restart_frequency) then
      write (fname,tadonefname_fmt) itcur,""   ! file name only
      write (fullname,tadone_fmt) itcur,""     ! file with path
      write (link_current_command,talkdn_fmt) itcur,""
      write (statfile,tadone_fmt) itcur,""
      statfile = trim(adjustl(statfile)) // ".stat"
    else
      write (fname,tafname_fmt) averaging_file_switch,""   ! file name only
      write (fullname,tafile_fmt) averaging_file_switch,"" ! file with path
      write (link_current_command,talink_fmt) averaging_file_switch,""
      write (statfile,tastat_fmt) averaging_file_switch,""
    end if
  else
    write (fname,tadonefname_fmt) itcur,".done"
    write (fullname,tadone_fmt) itcur,".done"
    write (link_current_command,talkdn_fmt) itcur,".done"
   !write (statfile,tastat_fmt) itcur,".done"
    write (statfile,tadone_fmt) itcur,".done"
    statfile = trim(adjustl(statfile)) // ".stat"
  end if
  !
  ! Create the time averaged information to write to the averaging file
  !
  if (restart_output_uses_host_roots) then
    involved_with_output = i_am_host_root
    timeave_comm = host_roots_comm
    call collect_restart_arrays(current_sol=uavesp, &
                                output_sol=solution)
  else
    involved_with_output = true
    timeave_comm = MPI_COMM_WORLD
    call get_restart_arrays(current_sol=uavesp, &
                            output_sol=solution)
  end if
  !
  ! Store the size information specific to the time averaged solution
  ! into the respective location within the int_header array
  !
  int_header(1) = size(uavesp,dim=2)
  !
  ! Store the time data specific to the time averaged solution
  ! into the respective location within the real_header array
  !
  real_header(:) = [ time , prev_ave_time , ave_start_time , rl2mx ]
  !
  ! Create the variable header for the time averaged solution
  !
  allocate ( var_header(1:size(time_ave_variables)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"var_header",1,__LINE__,__FILE__,ierr,error_message)
  !
  do n = 1,size(time_ave_variables)
    var_header(n) = transfer(time_ave_variables(n),var_header(n))
  end do
  !
  ! Now write the actual time averaged solution file using all the
  ! MPI processes actually involved with the I/O
  !
  if (involved_with_output) then
    !
    if (ncpu > 1) then
      call write_solution_file_mpi(fullname,int_header,real_header, &
                                   timeave_comm, &
                                   ordered_geom,ordered_geom_type, &
                                   solution,averaging_type,ave_block_type, &
                                   max_time_ave_rst_wall_time, &
                                   var_header=var_header,statfile=statfile)
    else
      call write_solution_file_serial(fullname,int_header,real_header, &
                                      ordered_geom,solution, &
                                      max_time_ave_rst_wall_time, &
                                      var_header=var_header)
    end if
    !
  end if
#ifdef USE_CREATE_SYMLINK
  !
  ! Create a link to the most recent restart file for piggy-backing jobs
  !
  if (mypnum == glb_root) then
    call create_restart_link(fname,iteration=itcur)
  end if
#else
  !
  ! Create a link to the most recent restart file for piggy-backing jobs
  !
  call execute_command(link_current_command,wait_for_completion,exitstat, &
                       cmndstat,cmndmsg,pname)
  !
  if (.not. is_at_restart_frequency) then
    !
    if (output_iterations(averaging_file_switch) /= 0) then
      n = output_iterations(averaging_file_switch)
      write (rm_iter_link_command,rmtalink_fmt) n
    else
      rm_iter_link_command = ''
    end if
    write (link_iter_command,taiterlk_fmt) averaging_file_switch, itcur
    output_iterations(averaging_file_switch) = itcur
    !
    ! Create a link to this time averaged restart file that contains the
    ! iteration number in the link files name.
    !
    call execute_command(link_iter_command,wait_for_completion,exitstat, &
                         cmndstat,cmndmsg,pname)
    !
    ! Remove the old iteration link if we have begun to start
    ! overwriting previous time average restart files.
    !
    if (len_trim(adjustl(rm_iter_link_command)) > 0) then
      call execute_command(rm_iter_link_command,wait_for_completion,exitstat, &
                           cmndstat,cmndmsg,pname)
    end if
    !
  end if
  !
#endif
  !
  ! Deallocate the solution array since its no longer needed
  !
  if (allocated(solution)) then
    deallocate ( solution , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"solution",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  ! Toggle the averaging for the next time we write the averaged solution
  ! NOTE: num_saved_time_ave_restarts should be > 0, but we will use abs
  !       just to make sure that is the case.
  !
  if (.not. is_at_restart_frequency) then
    if (averaging_file_switch >= abs(num_saved_time_ave_restarts)) then
      averaging_file_switch = 1
    else
      averaging_file_switch = averaging_file_switch + 1
    end if
  end if
  !
  call _DEBUG_WRITE_TIMER_(leaving_procedure,pname)
  !
end subroutine write_time_averaged_solution
!
!###############################################################################
!
subroutine write_modified_time_averaged_solution(file_name,int_header, &
                                                 real_header,var_header, &
                                                 ordered_geom,solution)
  !
  !.. Use Statements ..
  use ovar, only : restart_output_uses_host_roots
  use ovar, only : max_time_ave_rst_wall_time
  use ovar, only : loc_solution_pts
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: file_name
  integer(wip),     intent(in) :: int_header(:)
  real(wp),         intent(in) :: real_header(:)
  integer,          intent(in) :: var_header(:)
  integer(ogp),     intent(in) :: ordered_geom(:)
  real(wp),         intent(in) :: solution(:,:)
  !
  !.. Local Scalars ..
  integer     :: ierr
  logical(lk) :: involved_with_output
  !
  _MPI_COMM_TYPE_ :: timeave_comm
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: fullname,statfile
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: output_sol(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_modified_time_averaged_solution"
  !
continue
  !
  call _DEBUG_WRITE_TIMER_(entering_procedure,pname)
  !
  ! Perform some initialization tasks if the restart module has not
  ! been initialized yet.
  ! NOTE: This is usually done to prepare for writing the instantaneous
  !       restart file, but we need to do it here because a few module
  !       arrays are created that will be needed to get the MPI datatypes
  !       for writing the time-average restart files.
  !
  if (restart_needs_initialization) then
    !
    ! For now, write the nodal solution to the restart file. At some point we
    ! should switch so it writes the modal solution since that would give more
    ! flexibility regarding changes in the solution order between runs.
    !
    if (loc_solution_pts == Legendre_Gauss) then
      restart_solution_form = restart_gauss_nodal_solution
    else if (loc_solution_pts == Legendre_Gauss_Lobatto) then
      restart_solution_form = restart_lobatto_nodal_solution
    end if
    !
    ! Get the MPI derived types required to read/write the restart file
    !
    if (restart_output_uses_host_roots) then
      call get_host_output_mpi_datatypes
    else
      call get_output_mpi_datatypes
    end if
    !
#ifdef DEBUG_ON
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      write (iout,*) "After get_output_mpi_datatypes"
      flush (iout)
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
    !
    ! Unset the restart initialization flag since this has now been done
    !
    restart_needs_initialization = fals
    !
  end if
  !
  ! Perform some initialization tasks to get the MPI datatypes needed to
  ! write the time-average restart files
  !
  if (time_ave_needs_initialization) then
    !
    ! Get the MPI derived types required to
    ! read/write the time-averaging restart file
    !
    if (restart_output_uses_host_roots) then
      call get_host_time_ave_mpi_datatypes
    else
      call get_time_ave_mpi_datatypes
    end if
    !
    ! Unset the restart initialization flag since this has now been done
    !
    time_ave_needs_initialization = fals
    !
  end if
  !
  ! Get the name for the time
  !
  ! F2003 AUTO-REALLOC
  fullname = trim(adjustl(restart_dir)) // trim(adjustl(file_name)) 
  ! F2003 AUTO-REALLOC
  statfile = trim(adjustl(fullname)) // ".stat"
  !
  ! Create the time averaged information to write to the averaging file
  !
  if (restart_output_uses_host_roots) then
    involved_with_output = i_am_host_root
    timeave_comm = host_roots_comm
    call collect_restart_arrays(current_sol=solution, &
                                output_sol=output_sol)
  else
    involved_with_output = true
    timeave_comm = MPI_COMM_WORLD
    call get_restart_arrays(current_sol=solution, &
                            output_sol=output_sol)
  end if
  !
  ! Now write the actual time averaged solution file using all the
  ! MPI processes actually involved with the I/O
  !
  if (involved_with_output) then
    !
    if (ncpu > 1) then
      call write_solution_file_mpi(fullname,int_header,real_header, &
                                   timeave_comm, &
                                   ordered_geom,ordered_geom_type, &
                                   output_sol,averaging_type,ave_block_type, &
                                   max_time_ave_rst_wall_time, &
                                   var_header=var_header,statfile=statfile)
    else
      call write_solution_file_serial(fullname,int_header,real_header, &
                                      ordered_geom,output_sol, &
                                      max_time_ave_rst_wall_time, &
                                      var_header=var_header)
    end if
    !
  end if
  !
  ! Deallocate the output_sol array since its no longer needed
  !
  if (allocated(output_sol)) then
    deallocate ( output_sol , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"output_sol",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  call _DEBUG_WRITE_TIMER_(leaving_procedure,pname)
  !
end subroutine write_modified_time_averaged_solution
!
!###############################################################################
!
subroutine get_restart_arrays(current_sol,output_sol,ordered_geom)
  !
  !.. Use Statements ..
  use geomvar_mod, only : ncell,cell,n_solpts
  !
  !.. Formal Arguments ..
  real(wp),                 intent(in) :: current_sol(:,:)
  real(wp), allocatable, intent(inout) :: output_sol(:,:)
  !
  !.. Optional Arguments ..
  integer(ogp), optional, allocatable, intent(inout) :: ordered_geom(:)
  !
  !.. Local Scalars ..
  integer :: k,l0,l2,m,n,nc,np,nvar,n0,ierr
  integer :: this_geom,this_order
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_restart_arrays"
  !
continue
  !
  call _DEBUG_WRITE_TIMER_(entering_procedure,pname)
  !
  if (present(ordered_geom)) then
    !
    allocate ( ordered_geom(1:ncell) , source=0_ogp , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ordered_geom",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    do n = 1,ncell
      !
      nc = cell_output_order(n)
      !
      this_geom = cell(nc)%geom
      this_order = cell(nc)%order
      !
      ordered_geom(n) = int( 100*this_order + this_geom , kind=ogp )
      !
    end do
    !
  end if
  !
  ! Allocate a temporary array to copy over only the interior solution
  ! data for each processor. We also have to copy over the solution data
  ! in the same order that is expected by the MPI derived datatype.
  !
  nvar = size(current_sol,dim=2)
  !
  allocate ( output_sol(1:nvar,1:n_solpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"output_sol",1,__LINE__,__FILE__,ierr,error_message)
  !
  l2 = 0
  do n = 1,ncell
    !
    nc = cell_output_order(n)
    !
    n0 = cell(nc)%beg_sp-1
    np = cell(nc)%end_sp-n0
    !
    l0 = l2
    l2 = l2 + np
    !
    ! We need to flip the indices between the two arrays to match what
    ! is expected by the MPI derived datatype.
    !
    do k = 1,np
      do m = 1,nvar
        output_sol(m,l0+k) = current_sol(n0+k,m)
      end do
    end do
    !
  end do
  !
  call _DEBUG_WRITE_TIMER_(leaving_procedure,pname)
  !
end subroutine get_restart_arrays
!
!###############################################################################
!
subroutine collect_restart_arrays(current_sol,output_sol,ordered_geom)
  !
  !.. Use Statements ..
  use geomvar_mod, only : ncell,cell,n_solpts
  !
  !.. Formal Arguments ..
  real(wp),                 intent(in) :: current_sol(:,:)
  real(wp), allocatable, intent(inout) :: output_sol(:,:)
  !
  !.. Optional Arguments ..
  integer(ogp), optional, allocatable, intent(inout) :: ordered_geom(:)
  !
  !.. Local Scalars ..
  integer :: nvar,ierr,lcell,lpts
  integer :: ii
  !
  _MPI_DATA_TYPE_ :: og_mpi_type
  _MPI_COMM_TYPE_ :: debug_barrier_comm
  !
  !.. LOcal Allocatable Arrays ..
  integer(ogp), allocatable :: local_og(:)
  real(wp),     allocatable :: local_sol(:,:)
  !
  integer(INT_MPI), allocatable :: sol_counts(:)
  integer(INT_MPI), allocatable :: sol_displs(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "collect_restart_arrays"
  !
continue
  !
  call _DEBUG_WRITE_TIMER_(entering_procedure,pname)
  !
 !debug_barrier_comm = MPI_COMM_WORLD
  debug_barrier_comm = my_host_comm
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  if (present(ordered_geom)) then
    !
    lcell = 1
    if (i_am_host_root) then
      lcell = sum(og_counts)
    end if
    !
    allocate ( ordered_geom(1:lcell) , source=0_ogp , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ordered_geom",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    ! F2003 AUTO-REALLOCATION
#ifdef SPECIAL_FOR_GCC
    local_og = [int( 100*cell(1:ncell)%order + cell(1:ncell)%geom , kind=ogp )]
#else
    local_og = int( 100*cell(1:ncell)%order + cell(1:ncell)%geom , kind=ogp )
#endif
    !
    if (ogp == int64) then
      og_mpi_type = MPI_INTEGER8
    else
      og_mpi_type = MPI_INTEGER
    end if
    !
    call mpi_gatherv(local_og,int(ncell,kind=int_mpi),og_mpi_type, &
                     ordered_geom,og_counts,og_displs,og_mpi_type, &
                     host_root,my_host_comm,mpierr)
    !
    deallocate ( local_og , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"local_og",2,__LINE__,__FILE__,ierr,error_message)
    !
    ! Reorder the ordered_geom array to match its MPI derived datatype
    !
    if (i_am_host_root) then
      local_og = ordered_geom ! F2003 AUTO-REALLOCATION
      ordered_geom(:) = local_og(cell_output_order)
    end if
    !
  end if
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  nvar = size(current_sol,dim=2)
  !
  sol_counts = pts_counts * nvar ! F2003 AUTO-REALLOCATION
  sol_displs = pts_displs * nvar ! F2003 AUTO-REALLOCATION
  !
  lpts = 1
  if (i_am_host_root) then
    lpts = sum(pts_counts)
  end if
  !
  allocate ( output_sol(1:nvar,1:lpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"output_sol",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate a temporary array to copy over only the interior solution
  ! data for each processor. We also have to copy over the solution data
  ! in the same order that is expected by the MPI derived datatype.
  !
  ! F2003 AUTO-REALLOCATION
  local_sol = transpose(current_sol(1:n_solpts,1:nvar))
  !
  call mpi_gatherv(local_sol,size(local_sol,kind=int_mpi),mpi_flttyp, &
                   output_sol,sol_counts,sol_displs,mpi_flttyp, &
                   host_root,my_host_comm,mpierr)
  !
  deallocate ( local_sol , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"local_sol",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Reorder the solution array to match its MPI derived datatypes
  !
  if (i_am_host_root) then
   !local_sol = output_sol ! F2003 AUTO-REALLOCATION
   !output_sol(:,:) = reorder_host_solution( local_sol(:,:) )
    output_sol(:,:) = reorder_host_solution( output_sol(:,:) )
  end if
  !
  call _DEBUG_WRITE_TIMER_(leaving_procedure,pname)
  !
contains
!===============================================================================
subroutine debug_write(string,ival)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: string
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: ival
  !
  !.. Local Scalars ..
  integer :: icpu
  !
continue
  !
  ii = ii + 1
  !
  if (mypnum == glb_root) then
    write (iout,2)
    flush (iout)
  end if
  call mpi_barrier(debug_barrier_comm,mpierr)
  !
  flush (iout)
  call mpi_barrier(debug_barrier_comm,mpierr)
  !
  do icpu = 1,ncpu
    if (mypnum == icpu-1) then
      if (present(ival)) then
        write (iout,1) ii,mypnum,string,ival
      else
        write (iout,1) ii,mypnum,string
      end if
      flush (iout)
    end if
    call mpi_barrier(debug_barrier_comm,mpierr)
    flush (iout)
  end do
  !
  flush (iout)
  call mpi_barrier(debug_barrier_comm,mpierr)
  flush (iout)
  !
  1 format (i2.2,1x,"CPU # ",i3.3,": ",a,:," = ",i0)
  2 format (100("-"))
  !
end subroutine debug_write
!===============================================================================
end subroutine collect_restart_arrays
!
!###############################################################################
!
subroutine debug_write_module(ii,string,ival,comm)
  !
  !.. Formal Arguments ..
  integer, intent(inout) :: ii
  character(len=*), intent(in) :: string
  integer, intent(in) :: ival
  _MPI_COMM_TYPE_, intent(in) :: comm
  !
  !.. Local Scalars ..
  integer :: icpu
  !
continue
  !
  ii = ii + 1
  !
  flush (iout)
  call mpi_barrier(comm,mpierr)
  !
  do icpu = 1,ncpu
    if (mypnum == icpu-1) then
      write (iout,1) ii,mypnum,string,ival
      flush (iout)
    end if
    call mpi_barrier(comm,mpierr)
    flush (iout)
  end do
  !
  flush (iout)
  call mpi_barrier(comm,mpierr)
  flush (iout)
  !
  1 format (i2.2,1x,"CPU # ",i3.3,": ",a,:," = ",i0)
  !
end subroutine debug_write_module
!
!###############################################################################
!
pure function reorder_host_solution(solution) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: solution(:,:)
  !
  !.. Function Result ..
  ! NOTE: Changed the function result to allocatable because having it as an
  !       automatic array based on the size of the input solution caused the
  !       compiler to put return_value on the stack, which sometimes caused
  !       the code to crash due to it running out of stack memory. Changing
  !       return_value to allocatable causes the compiler to put the array
  !       in heap memory.
  real(wp), allocatable :: return_value(:,:)
 !real(wp) :: return_value(1:size(solution,dim=1), &
 !                         1:size(solution,dim=2))
  !
  !.. Local Scalars ..
  integer :: nc,n1,n2,l1,l2
  !
continue
  !
  n1 = size(solution,dim=1)
  n2 = size(solution,dim=2)
  allocate ( return_value(1:n1,1:n2) , source=zero )
  !
  l2 = 0
  do nc = 1,size(host_output_order)
    !
    n1 = host_output_order(nc)%beg_sp
    n2 = host_output_order(nc)%end_sp
    !
    l1 = l2 + 1
    l2 = l2 + n2-n1+1
    !
    return_value(:,l1:l2) = solution(:,n1:n2)
    !
  end do
  !
end function reorder_host_solution
!
!###############################################################################
!
subroutine write_restart_status(ios,file_comm,file_name,stat_strng, &
                                close_file,open_file,mpistat,mpi_error)
  !
  !.. Formal Arguments ..
  integer,       intent(inout) :: ios
  _MPI_COMM_TYPE_,  intent(in) :: file_comm
  !
  !.. Optional Arguments ..
  character(len=*),  optional, intent(in) :: file_name
  character(len=*),  optional, intent(in) :: stat_strng
  logical(lk),       optional, intent(in) :: open_file
  logical(lk),       optional, intent(in) :: close_file
  _MPI_STATUS_TYPE_, optional, intent(in) :: mpistat(:)
  integer(int_mpi),  optional, intent(in) :: mpi_error
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_restart_status"
  !
continue
  !
#ifdef VERBOSE_RESTART
  !
#ifdef ONLY_READ_RESTARTS_ARE_VERBOSE
  if (.not. restart_operation_is_read) then
    return
  end if
#endif
  !
#ifdef ONLY_WRITE_RESTARTS_ARE_VERBOSE
  if (restart_operation_is_read) then
    return
  end if
#endif
  !
  if (present(file_name)) then
    !
    if (present(open_file)) then
      if (open_file) then
        if (mypnum == glb_root) then
          open (newunit=ios,file=trim(adjustl(file_name)),status="replace", &
                            action="write",form="formatted",iostat=ierr, &
                            iomsg=error_message)
          call io_error(pname,file_name,1,__LINE__,__FILE__,ierr,error_message)
        end if
        call mpi_barrier(file_comm,mpierr)
        if (mypnum == glb_root) flush (ios)
        call mpi_barrier(file_comm,mpierr)
      end if
    end if
    !
    if (present(stat_strng).or.present(mpistat).or.present(mpi_error)) then
      !
      if (mypnum == glb_root) then
        write (ios,2)
        flush (ios)
      end if
      call mpi_barrier(file_comm,mpierr)
      if (mypnum == glb_root) then
        if (present(stat_strng)) write (ios,1) mypnum,trim(adjustl(stat_strng))
        if (present(mpistat)) then
#ifdef USE_MPI_F08
          write (ios,3) mypnum,mpistat(1)%mpi_source, &
                               mpistat(1)%mpi_tag, &
                               mpistat(1)%mpi_error
#else
          write (ios,3) mypnum,mpistat(1),mpistat(2),mpistat(3)
#endif
        end if
        if (present(mpi_error)) write (ios,4) mypnum,mpi_error
      end if
      call mpi_barrier(file_comm,mpierr)
      if (mypnum == glb_root) flush (ios)
      call mpi_barrier(file_comm,mpierr)
      if (mypnum == glb_root) then
        write (ios,2)
        flush (ios)
      end if
      !
    end if
    !
    if (present(close_file)) then
      if (close_file) then
        call mpi_barrier(file_comm,mpierr)
        if (mypnum == glb_root) flush (ios)
        call mpi_barrier(file_comm,mpierr)
        if (mypnum == glb_root) then
          close (ios,iostat=ierr,iomsg=error_message)
          call io_error(pname,file_name,2,__LINE__,__FILE__,ierr,error_message)
        end if
      end if
    end if
    !
  end if
  !
  ! Format Statements
  !
  1 format ("CPU # ",i5.5," : ",a)
  2 format (80("-"))
  3 format ("CPU # ",i5.5," : MPI_SOURCE = ",i0,/, &
                         14x,"   MPI_TAG = ",i0,/, &
                         14x," MPI_ERROR = ",i0)
  4 format ("CPU # ",i5.5," Error code returned from MPI function (mpierr): ",i0)
  !
#endif
  !
end subroutine write_restart_status
!
!###############################################################################
!
subroutine check_restart_dir
  !
  !.. Use Statements ..
  use ovar, only : output_dir,itrst,num_timesteps
  use ovar, only : lustre_stripe_count,lustre_stripe_size
  use ovar, only : load_restart_file,restart_file
  use ovar, only : output_time_averaging
  use ovar, only : time_ave_file,read_time_ave_restart
  use ovar, only : restart_output_uses_host_roots
  use ovar, only : matrix_output_unit
  !
  !.. Local Scalars ..
  integer :: clen,ierr,iodir
  integer :: exitstat,cmndstat
  integer :: final_iter_digits
  logical(lk)  :: using_restart_dir
  logical(ldk) :: file_exists
  character(len=300) :: command
  character(len=300) :: cmndmsg
  character(len=200) :: test_file
  character(len=200) :: test_dir
  character(len=100) :: char_num
  character(len=150) :: fname
  !
  !.. Local Arrays ..
  logical(ldk) :: restart_logicals(1:2)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_restart_dir"
  !
  logical(ldk), parameter :: wait_for_completion = .true.
  !
continue
  !
  if (matrix_output_unit /= iout) return
  !
  call global_debug_timer(entering_procedure,pname)
  !
  ! Initialize exitstat, cmndmsg, and restart_logicals
  !
  exitstat = 0
  cmndmsg = ""
  restart_logicals(:) = .false.
  !
  ! First, if a restart was requested, make
  ! sure that the input restart file exists
  !
  if (mypnum == glb_root) then
    !
    if (load_restart_file) then
      !
      fname = trim(adjustl(restart_file))
      !
      inquire (file=fname,exist=file_exists)
      !
      if (.not. file_exists) then
        write (iout,100) trim(adjustl(fname))
        flush (iout)
        load_restart_file = fals
       !call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
      end if
      !
    end if
    !
    read_time_ave_restart = .false.
    if (load_restart_file) then
      !
      if (output_time_averaging) then
        !
        fname = trim(adjustl(time_ave_file))
        !
        inquire (file=fname,exist=read_time_ave_restart)
        !
      end if
      !
    end if
    !
    restart_logicals(:) = [ logical(load_restart_file,kind=ldk) , &
                            read_time_ave_restart ]
    !
  end if
  !
  call mpi_bcast(restart_logicals,size(restart_logicals,kind=int_mpi), &
                 MPI_LOGICAL,glb_root,MPI_COMM_WORLD,mpierr)
  !
  load_restart_file = logical( restart_logicals(1) , kind=lk )
  read_time_ave_restart = restart_logicals(2)
  !
  ! Try to open a file in the directory named Restart_files
  !
  test_dir = trim(output_dir) // "Restart_files/"
  test_file = trim(test_dir) // "Restart_dir.test_file.dat"
  !
  ! Only perform the following tests using the root processor
  ! to minimize contention and prevent each MPI process from
  ! executing the same command.
  !
  if (mypnum == glb_root) then
    !
    open (newunit=iodir,file=trim(test_file),status="replace", &
                        action="write",iostat=ierr)
    !
    if (ierr == 0) then
      !
      ! The open was successful, so close and delete the test file
      !
      close (iodir,status="delete",iostat=ierr)
      !
      ! Set the logical saying that we are using a separate restart directory
      !
      using_restart_dir = true
      restart_dir = trim(test_dir)
      !
    else
      !
      ! The open failed, so more than likely the directory Restart_files
      ! doesnt exist. Try to use the intrinsic execute_command_line
      ! to create the directory.
      !
      command = "mkdir -p "//trim(test_dir)
      !
      call execute_command(command,wait_for_completion,exitstat, &
                           cmndstat,cmndmsg,pname)
     !call execute_command_line( trim(command) , wait_for_completion , &
     !                           exitstat , cmndstat , cmndmsg )
     !call write_command_results(pname,command,exitstat,cmndstat,cmndmsg)
      !
      if (exitstat /= 0 .or. cmndstat /= 0) then
        using_restart_dir = fals
        restart_dir = trim(output_dir)
      else
        using_restart_dir = true
        restart_dir = trim(test_dir)
      end if
      !
    end if
    !
#ifdef PBS_ENV
    !
    ! If we are using a separate restart directory and we are running
    ! this through PBS on NAS, try to modify the lfs striping for
    ! the restart directory for better I/O performance.
    !
    if (using_restart_dir) then
      !
     !command = "lfs setstripe -c 64 -S 4m "//trim(restart_dir)
      write (command,3) max(0,lustre_stripe_count), &
                        trim(adjustl(lustre_stripe_size)), &
                        trim(adjustl(restart_dir))
      write (iout,4) trim(adjustl(command))
      flush (iout)
      !
      call execute_command(command,wait_for_completion,exitstat, &
                           cmndstat,cmndmsg,pname)
     !call execute_command_line( trim(command) , wait_for_completion , &
     !                           exitstat , cmndstat , cmndmsg )
     !call write_command_results(pname,command,exitstat,cmndstat,cmndmsg)
      !
      ! We arent going to test the status of executing this command because
      ! it worked or it didnt. Either way, we will continue on and output
      ! the restart solution to the restart directory.
      !
    end if
    !
#endif
    !
  end if
  !
  ! Broadcast the value of restart_dir to all the non-root MPI processes
  !
  clen = len_trim(restart_dir)
  call mpi_bcast(clen,1_int_mpi,mpi_inttyp,glb_root,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(restart_dir(1:clen),int(clen,kind=int_mpi),MPI_CHARACTER, &
                 glb_root,MPI_COMM_WORLD,mpierr)
  !
  ! Create the various format strings needed for reading/writing restart files
  !
  call create_restart_fname_formats(itrst)
  !
  ! Make sure restart_output_uses_host_roots is
  ! false if only one processor is being used.
  !
  if (ncpu <= 1) then
    restart_output_uses_host_roots = fals
  end if
  !
  call global_debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  3  format ("lfs setstripe -c ",i0," -S ",a,1x,a)
  4  format ("Lustre striping command = '",a,"'")
  !
  100 format (/,1x,"A simulation restart was requested but", &
                1x,"the restart file given in the input file",/, &
                10x,"'",a,"'",/, &
                1x,"does not exist!!!",/, &
                1x,"NOTE: The file 'out.rst' is the default", &
                1x,"if no file name is explicitly given in",/, &
                1x,"the input file. The simulation will proceed", &
                1x,"as if a restart was not requested.")
  !
end subroutine check_restart_dir
!
!###############################################################################
!
subroutine create_restart_fname_formats(starting_iteration)
  !
  !.. Use Statements ..
  use ovar, only : num_timesteps
  use ovar, only : num_saved_time_ave_restarts
  !
  !.. Formal Arguments ..
  integer, intent(in) :: starting_iteration
  !
  !.. Local Scalars ..
  integer :: final_iter_digits,time_ave_digits
  character(len=100) :: char_num
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: rstdir
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_restart_fname_formats"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Get the final/maximum iteration number for this simulation and extract
  ! the number of digits needed for an integer of this size.
  !
  write (char_num,1) starting_iteration + num_timesteps
  final_iter_digits = max(4,len_trim(adjustl(char_num)))
  !
  ! Get the number of digits for time averaged restart file
  ! NOTE: num_saved_time_ave_restarts should be > 0, but we will use abs
  !       just to make sure that is the case.
  !
  write (char_num,1) abs(num_saved_time_ave_restarts)
  time_ave_digits = max(1,len_trim(adjustl(char_num)))
  !
  ! Create the format strings for creating the restart file names and commands
  !
  rstdir = trim(adjustl(restart_dir)) ! F2003 AUTO-REALLOC
  !
  write (rstfile_fmt,12) rstdir, final_iter_digits, "rst"
  write (lnkcmnd_fmt,13) final_iter_digits, "rst", rstdir, "rst"
  write (statfile_fmt,12) rstdir, final_iter_digits, "stat"
  write (tafile_fmt,12) rstdir, time_ave_digits, "ave"
  write (talink_fmt,13) time_ave_digits, "ave", rstdir, "ave"
  write (tastat_fmt,12) rstdir, time_ave_digits, "ave.stat"
  write (tadone_fmt,12) rstdir, final_iter_digits, "ave"
  write (talkdn_fmt,13) final_iter_digits, "ave", rstdir, "ave"
  write (taiterlk_fmt,14) time_ave_digits, rstdir, final_iter_digits
  write (rmtalink_fmt,15) rstdir, final_iter_digits
  !
  write (rstfname_fmt,12) "", final_iter_digits, "rst"
  write (tafname_fmt,12) "", time_ave_digits, "ave"
  write (tadonefname_fmt,12) "", final_iter_digits, "ave"
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (i0)
  12 format ("('",a,"out.',i0.",i0,",a,'.",a,"')")
  13 format ("('/bin/ln -sf out.',i0.",i0,",a,'.",a," ",a,"out.current.",a,"')")
  14 format ("('/bin/ln -sf out.',i0.",i0,",'.ave ",a,"out.',i0.",i0,".ave')")
  15 format ("('/bin/rm -f ",a,"out.',i0.",i0,",'.ave')")
  !
end subroutine create_restart_fname_formats
!
!###############################################################################
!
subroutine create_restart_link(fname,iteration)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  use ovar, only : is_at_restart_frequency
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: fname
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: iteration
  !
  !.. Local Scalars ..
  integer :: l,n
  integer(c_int) :: ierr
  logical(ldk) :: file_exists
  logical(lk) :: linking_ta_file
  character(len=300) :: ta_fname
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: existing_file
  character(len=:), allocatable :: link_file
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_restart_link"
  !
continue
  !
  linking_ta_file = fals
  !
  existing_file = trim(adjustl(fname)) ! F2003 AUTO-REALLOC
  !
  n = len_trim(existing_file)
  !
  if (uppercase(existing_file(n-2:n)) == "AVE") then
    ! F2003 AUTO-REALLOC
    link_file = trim(adjustl(restart_dir)) // "out.current.ave"
    if (uppercase(existing_file(n-7:n)) /= "DONE.AVE") then
      linking_ta_file = true
    end if
  else
    ! F2003 AUTO-REALLOC
    link_file = trim(adjustl(restart_dir)) // "out.current.rst"
  end if
  !
  call create_symlink(existing_file,link_file)
  !
  if (linking_ta_file .and. .not. is_at_restart_frequency) then
    if (present(iteration)) then
      !
      write (ta_fname,tadone_fmt) iteration,""
      link_file = trim(adjustl(ta_fname)) ! F2003 AUTO-REALLOC
      !
      call create_symlink(existing_file,link_file)
      !
      if (output_iterations(averaging_file_switch) /= 0) then
        !
        write (ta_fname,tadone_fmt) output_iterations(averaging_file_switch),""
        link_file = trim(adjustl(ta_fname)) ! F2003 AUTO-REALLOC
        !
        ierr = unlink_f(link_file)
        if (ierr /= 0) then
          write (iout,1) pname,link_file,ierr
          flush (iout)
        end if
        !
      end if
      !
      output_iterations(averaging_file_switch) = iteration
      !
    end if
  end if
  !
  1 format (/,2x,"ERROR RETURNED FROM C FUNCTION unlink!!!",/, &
              6x,"    procedure = '",a,"'",/, &
              6x,"    link_file = '",a,"'",/, &
              6x,"         ierr = ",i0,/)
  !
end subroutine create_restart_link
!
!###############################################################################
!
subroutine create_symlink(existing_file,link_file)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: existing_file
  character(len=*), intent(in) :: link_file
  !
  !.. Local Scalars ..
  integer :: l
  integer(c_int) :: ierr
  logical(ldk) :: file_exists
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: efile
  character(len=:), allocatable :: lfile
  character(len=:), allocatable :: efile_with_path
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_symlink"
  !
continue
  !
  efile = trim(adjustl(existing_file)) ! F2003 AUTO-REALLOC
  lfile = trim(adjustl(link_file))     ! F2003 AUTO-REALLOC
  !
  efile_with_path = trim(adjustl(restart_dir)) // efile ! F2003 AUTO-REALLOC
  !
  inquire (file=efile_with_path,exist=file_exists)
  if (.not. file_exists) then
    write (iout,1) pname,efile_with_path,lfile
    flush (iout)
    return
  end if
  !
  retry_loop: &
  do l = 1,5
    !
    ierr = unlink_f(lfile)
    if (ierr /= 0) then
      write (iout,2) pname,lfile,ierr
      flush (iout)
    end if
    !
    ierr = symlink_f(efile,lfile)
    !
    if (ierr /= 0) then
      write (iout,3) pname,efile,lfile,ierr
      flush (iout)
    else
      exit retry_loop
    end if
    !
  end do retry_loop
  !
  1 format (/,2x,"ERROR CREATING SYMLINK FOR RESTART FILE BECAUSE ", &
                 "RESTART FILE DOES NOT EXIST!!!",/, &
              6x,"    procedure = '",a,"'",/, &
              6x,"existing_file = '",a,"'",/, &
              6x,"    link_file = '",a,"'",/, &
              2x,"ABORTING THIS ATTEMPT TO CREATE A SYMLINK FILE AND ", & 
                 "CONTINUING ON!!!",/)
  2 format (/,2x,"ERROR RETURNED FROM C FUNCTION unlink!!!",/, &
              6x,"    procedure = '",a,"'",/, &
              6x,"    link_file = '",a,"'",/, &
              6x,"         ierr = ",i0,/)
  3 format (/,2x,"ERROR RETURNED FROM C FUNCTION symlink!!!",/, &
              6x,"    procedure = '",a,"'",/, &
              6x,"existing_file = '",a,"'",/, &
              6x,"    link_file = '",a,"'",/, &
              6x,"         ierr = ",i0,/)
  !
end subroutine create_symlink
!
!###############################################################################
!
function symlink_f(existing_file,link_file) result(return_value)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int,c_null_char
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: existing_file
  character(len=*), intent(in) :: link_file
  !
  !.. Function Result ..
  integer(c_int) :: return_value
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: efile_null
  character(len=:), allocatable :: lfile_null
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "symlink_f"
  !
  !.. C Function Interfaces ..
  interface
    integer(c_int) function symlink(existing_file,link_file) &
                               bind(c,name="symlink")
      use iso_c_binding, only : c_int,c_char
      implicit none
      character(kind=c_char), dimension(*), intent(in) :: existing_file
      character(kind=c_char), dimension(*), intent(in) :: link_file
    end function symlink
  end interface
  !
continue
  !
  efile_null = trim(adjustl(existing_file)) // c_null_char ! F2003 AUTO-REALLOC
  lfile_null = trim(adjustl(link_file)) // c_null_char     ! F2003 AUTO-REALLOC
  !
  return_value = symlink(efile_null,lfile_null)
  !
end function symlink_f
!
!###############################################################################
!
function unlink_f(link_file) result(return_value)
  !
  !.. Use Statements ..
  use, intrinsic :: iso_c_binding, only : c_int,c_null_char
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: link_file
  !
  !.. Function Result ..
  integer(c_int) :: return_value
  !
  !.. Local Scalars ..
  logical(ldk) :: file_exists
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: lfile
  character(len=:), allocatable :: lfile_null
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "unlink_f"
  !
  !.. C Function Interfaces ..
  interface
    integer(c_int) function unlink(link_file) &
                               bind(c,name="unlink")
      use iso_c_binding, only : c_int,c_char
      implicit none
      character(kind=c_char), dimension(*), intent(in) :: link_file
    end function unlink
  end interface
  !
continue
  !
  return_value = int(0,kind=c_int)
  !
  lfile = trim(adjustl(link_file))  ! F2003 AUTO-REALLOC
  lfile_null = lfile // c_null_char ! F2003 AUTO-REALLOC
  !
  inquire (file=lfile,exist=file_exists)
  if (file_exists) then
    return_value = unlink(lfile_null)
  end if
  !
end function unlink_f
!
!###############################################################################
!
subroutine get_host_output_mpi_datatypes()
  !
  !.. Use Statement ..
  use geomvar_mod, only : ncell,cell,n_solpts,cell_map
  use geomvar_mod, only : global_pinc_ptr
  use flowvar_mod, only : uavesp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,naq,l1,l2,ierr
  integer :: local_cell,lcell
  integer :: global_cell,this_part
  integer :: np,np_sum,np_offset
  !
  !.. Local Allocatable Arrays ..
  logical(lk),      allocatable :: msk(:)
  integer(int_mpi), allocatable :: og_offset(:)
  integer(int_mpi), allocatable :: tmp_copy(:)
  integer(int_mpi), allocatable :: recvcounts(:)
  integer(int_mpi), allocatable :: recvdispls(:)
  integer,          allocatable :: sp_idx(:)
  integer,          allocatable :: host_sp_idx(:)
  integer,          allocatable :: host_cells(:)
  integer,          allocatable :: host_pts(:)
 !type(sp_idx_t),   allocatable :: host_sol_idx(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_host_output_mpi_datatypes"
  !
continue
  !
  call _DEBUG_WRITE_TIMER_(entering_procedure,pname)
  !
  ! Get the number of cells for each of the
  ! processes that is in the current host group
  !
  allocate ( host_cells(1:size(my_host_grp)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"host_cells",1,__LINE__,__FILE__,ierr,error_message)
  !
  do n = 1,size(my_host_grp)
    this_part = my_host_grp(n)+1
    host_cells(n) = size(cell_map(this_part)%loc_to_glb)
  end do
  !
  lcell = sum(host_cells)
  !
  ! Get the number of solution points for each of
  ! the processes that is in the current host group
  !
  allocate ( host_pts(1:size(my_host_grp)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"host_pts",1,__LINE__,__FILE__,ierr,error_message)
  !
  call mpi_allgather(n_solpts,1_int_mpi,mpi_inttyp, &
                     host_pts,1_int_mpi,mpi_inttyp, &
                     my_host_comm,mpierr)
  !
  ! Allocate the arrays for the host root
  ! NOTE : Dont allocate these on every process to limit
  !        the amount of memory consumed on each host.
  !
  if (i_am_host_root) then
    !
    allocate ( og_offset(1:lcell) , source=0_int_mpi , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"og_offset",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( pts_offset(1:lcell) , source=0_int_mpi , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"pts_offset",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( pts_blklen(1:lcell) , source=0_int_mpi , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"pts_blklen",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( host_sp_idx(1:2*lcell) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"host_sp_idx",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( host_output_order(1:lcell) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"host_output_order",1,__LINE__,__FILE__, &
                     ierr,error_message)
    !
  else
    !
    ! Have all other processes allocate host_sp_idx so MPI doesnt complain
    ! that host_sp_idx isnt allocated when using mpi_gatherv
    !
    allocate ( host_sp_idx(1:1) , source=0 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"host_sp_idx",1,__LINE__,__FILE__,ierr,error_message)
    !
  end if
  !
  ! Create an array of the beginning and ending solution
  ! point indices for each cell on the local process
  !
  allocate ( sp_idx(1:2*ncell) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sp_idx",1,__LINE__,__FILE__,ierr,error_message)
  !
  n2 = 0
  do local_cell = 1,ncell
    n1 = n2 + 1
    n2 = n1 + 1
    sp_idx(n1) = cell(local_cell)%beg_sp
    sp_idx(n2) = cell(local_cell)%end_sp
  end do
  !
  ! Get the counts and displacements for collecting the ordered_geom
  ! array from each process in the current host group
  ! NOTE: We can create these arrays for all processes on this host since
  !       these are very small arrays and consume little to no memory.
  !
  ! F2003 AUTO-REALLOCATION
#ifdef SPECIAL_FOR_GCC
  og_counts = [ int(host_cells,kind=int_mpi) ]
#else
  og_counts = int(host_cells,kind=int_mpi)
#endif
  !
  ! F2003 AUTO-REALLOCATION
#ifdef SPECIAL_FOR_PGI
  allocate ( og_displs(1:size(og_counts)) , &
             source=eoshift( og_counts , shift=TO_END , boundary=0_int_mpi ) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"og_displs",1,__LINE__,__FILE__,ierr,error_message)
#else
  og_displs = eoshift( og_counts , shift=TO_END , boundary=0_int_mpi )
#endif
  !
  do n = 2,size(og_displs)
    og_displs(n) = og_displs(n) + og_displs(n-1)
  end do
  !
  recvcounts = 2_int_mpi*og_counts ! F2003 AUTO-REALLOCATION
  recvdispls = 2_int_mpi*og_displs ! F2003 AUTO-REALLOCATION
  !
  ! F2003 AUTO-REALLOCATION
#ifdef SPECIAL_FOR_GCC
  pts_counts = [ int(host_pts,kind=int_mpi) ]
#else
  pts_counts = int(host_pts,kind=int_mpi)
#endif
  !
  ! F2003 AUTO-REALLOCATION
#ifdef SPECIAL_FOR_PGI
  allocate ( pts_displs(1:size(pts_counts)) , &
             source=eoshift( pts_counts , shift=TO_END , boundary=0_int_mpi ), &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"pts_displs",1,__LINE__,__FILE__,ierr,error_message)
#else
  pts_displs = eoshift( pts_counts , shift=TO_END , boundary=0_int_mpi )
#endif
  !
  do n = 2,size(pts_displs)
    pts_displs(n) = pts_displs(n) + pts_displs(n-1)
  end do
  !
  ! Gather the solution point indices on the host root
  !
  call mpi_gatherv(sp_idx,size(sp_idx,kind=int_mpi),mpi_inttyp, &
                   host_sp_idx,recvcounts,recvdispls,mpi_inttyp, &
                   host_root,my_host_comm,mpierr)
  !
  deallocate ( sp_idx , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sp_idx",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( recvcounts , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recvcounts",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( recvdispls , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recvdispls",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( host_cells , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"host_cells",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( host_pts , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"host_pts",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the cell offsets as well as the block lengths
  ! and offsets for the solution data
  !
  if (i_am_host_root) then
    !
    nc = 0
    l2 = 0
    np_offset = 0
    !
    partition_loop: do n = 1,size(my_host_grp)
      !
      this_part = my_host_grp(n)+1
      !
      np_sum = 0 ! reset the running sum for the number of sol-pts on this_part
      !
      cell_loop: do local_cell = 1,size(cell_map(this_part)%loc_to_glb)
        !
        nc = nc + 1 ! cell counter for all partitions on this host
        !
        l1 = l2 + 1 ! location of beginning solpt idx within host_sp_idx
        l2 = l1 + 1 ! location of ending    solpt idx within host_sp_idx
        !
        ! Indices of the data for this cell in the solution array
        !
        n1 = host_sp_idx(l1) ! beginning solpt idx of local_cell on this_part
        n2 = host_sp_idx(l2) ! ending    solpt idx of local_cell on this_part
        np = n2-n1+1    ! total number of sol-pts for local_cell on this_part
        !
        ! Update the running sum of the total number of sol-pts on this_part
        !
        np_sum = np_sum + np
       !!
       !! Save these indices for later, relative to
       !! their location within the solution array
       !!
       !host_sol_idx(nc)%beg_sp = np_offset + n1
       !host_sol_idx(nc)%end_sp = np_offset + n2
        !
        ! Update these indices to their relative
        ! location within the solution array
        !
        host_sp_idx(l1) = np_offset + n1
        host_sp_idx(l2) = np_offset + n2
        !
        ! Store the block length for this local cell
        !
        pts_blklen(nc) = int( np , kind=int_mpi )
        !
        ! Get the global cell index for this local cell
        !
        global_cell = cell_map(this_part)%loc_to_glb(local_cell)
        !
        ! Store the global cell index for the ordered_geom offset
        !
        og_offset(nc) = int( global_cell-1 , kind=int_mpi )
        !
        ! Starting point of the solution data for this global cell
        ! NOTE: The value global_pinc_ptr(global_cell)+1 is actually the
        !       starting solution point for this cell. The offsets for MPI
        !       start at 0 so dropping the +1 is the actual offset we want.
        !
        pts_offset(nc) = int( global_pinc_ptr(global_cell) , kind=int_mpi )
        !
      end do cell_loop
      !
      ! Update the solution point offset for the next partition
      !
      np_offset = np_offset + np_sum
      !
    end do partition_loop
    !
    ! Allocate the cell_output_order array and
    ! the mask array for sorting the cell offsets
    !
    allocate ( cell_output_order(1:lcell) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_output_order",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    ! Sort the og_offset array to find the output order for the cells so that
    ! the offsets for the MPI datatypes are monotonically increasing.
    !
    call heap_sort_integer(lcell,og_offset,cell_output_order)
    !
    ! Get the solution point indices for each cell relative to their location
    ! within the solution array, and save them to their new reordered location
    !
    do nc = 1,lcell
      !
      n = cell_output_order(nc)
      !
      host_output_order(nc)%beg_sp = host_sp_idx(2*n-1)
      host_output_order(nc)%end_sp = host_sp_idx(2*n  )
      !
    end do
    !
   !allocate ( msk(1:lcell) , source=true , stat=ierr , errmsg=error_message )
   !call alloc_error(pname,"msk",1,__LINE__,__FILE__,ierr,error_message)
   !!
   !! Sort the og_offset array to find the output order for the cells so that
   !! the offsets for the MPI datatypes are monotonically increasing.
   !! NOTE: This array is assumed to be small enough so that it doesnt require
   !!       some fancy sorting algorithm and this simple minloc method will
   !!       more than likely be relatively inexpensive.
   !!
   !do nc = 1,lcell
   !  !
   !  n = minloc( og_offset , dim=1 , mask=msk )
   !  !
   !  msk(n) = fals ! mark this cell as found
   !  !
   !  cell_output_order(nc) = n
   !  !
   !  ! Get the solution point indices for this cell relative to their location
   !  ! within the solution array, and save them to their new reordered location
   !  !
   !  host_output_order(nc)%beg_sp = host_sp_idx(2*n-1)
   !  host_output_order(nc)%end_sp = host_sp_idx(2*n  )
   !  !
   !end do
    !
    ! Allocate the temporary array for applying the sorting to the offset
    ! and block length arrays.
    !
   !write (mypnum+1000,1) (cell_output_order(n),n=1,size(cell_output_order))
   !flush (mypnum+1000)
   !write (mypnum+1001,1) (og_offset(n),n=1,size(og_offset))
   !flush (mypnum+1001)
   !write (mypnum+1002,1) (pts_offset(n),n=1,size(pts_offset))
   !flush (mypnum+1002)
   !write (mypnum+1003,2) (host_sp_idx(2*n-1), &
   !                       host_sp_idx(2*n  ), &
   !                       n=1,lcell)
   !flush (mypnum+1003)
   !write (mypnum+1004,2) (cell_output_order(n), &
   !                       og_offset(cell_output_order(n)), &
   !                       n=1,size(og_offset))
   !flush (mypnum+1004)
   !write (mypnum+1005,2) (cell_output_order(n), &
   !                       pts_offset(cell_output_order(n)), &
   !                       n=1,size(pts_offset))
   !flush (mypnum+1005)
    !
    if (allocated(msk)) then
      deallocate ( msk , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    if (allocated(host_sp_idx)) then
      deallocate ( host_sp_idx , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"host_sp_idx",2,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
    !
    allocate ( tmp_copy(1:lcell) , source=0_int_mpi , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"tmp_copy",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Apply the sort to the og_offset array
    !
    tmp_copy(:) = og_offset
    og_offset(:) = tmp_copy(cell_output_order)
    !
    ! Apply the sort to the pts_offset array
    !
    tmp_copy(:) = pts_offset
    pts_offset(:) = tmp_copy(cell_output_order)
    !
    ! Apply the sort to the pts_blklen array
    !
    tmp_copy(:) = pts_blklen
    pts_blklen(:) = tmp_copy(cell_output_order)
   !write (mypnum+1004,1) (og_offset(n),n=1,size(og_offset))
   !flush (mypnum+1004)
   !write (mypnum+1005,1) (pts_offset(n),n=1,size(pts_offset))
   !flush (mypnum+1005)
   !write (mypnum+1006,2) (host_output_order(n)%beg_sp, &
   !                       host_output_order(n)%end_sp, &
   !                       n=1,size(host_output_order))
   !flush (mypnum+1006)
   !1 format (i0)
   !2 format (i9,1x,i9)
    !
    ! Deallocate the temporary array
    !
    deallocate ( tmp_copy , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"tmp_copy",2,__LINE__,__FILE__,ierr,error_message)
    !
    ! Create the MPI derived datatype for writing the
    ! ordered_geom array to the restart file
    !
    call mpi_type_create_indexed_block(size(og_offset,kind=int_mpi), &
                                       1_int_mpi,og_offset, &
                                       MPI_INTEGER,ordered_geom_type,mpierr)
    call mpi_type_commit(ordered_geom_type,mpierr)
    !
    ! Create the MPI derived datatype for writing the
    ! solution data to the restart file
    !
    call mpi_type_contiguous(int(nq,kind=int_mpi),mpi_flttyp, &
                             sol_block_type,mpierr)
    call mpi_type_commit(sol_block_type,mpierr)
    !
    !
    !
    call mpi_type_indexed(size(pts_offset,kind=int_mpi), &
                          pts_blklen,pts_offset, &
                          sol_block_type,solution_type,mpierr)
    call mpi_type_commit(solution_type,mpierr)
    !
  end if
  !
  ! Deallocate the local arrays before leaving
  !
  if (allocated(og_offset)) then
    deallocate ( og_offset , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"og_offset",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call _DEBUG_WRITE_TIMER_(leaving_procedure,pname)
  !
end subroutine get_host_output_mpi_datatypes
!
!###############################################################################
!
subroutine get_output_mpi_datatypes()
  !
  !.. Use Statement ..
  use geomvar_mod, only : ncell,cell,cell_map
  use geomvar_mod, only : global_pinc_ptr
  use flowvar_mod, only : uavesp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,naq,ierr
  integer :: local_cell
  integer :: global_cell
  !
  !.. Local Allocatable Arrays ..
  logical(lk),      allocatable :: msk(:)
  integer(int_mpi), allocatable :: og_offset(:)
  integer(int_mpi), allocatable :: sol_offset(:)
  integer(int_mpi), allocatable :: sol_blklen(:)
  integer(int_mpi), allocatable :: tmp_copy(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_output_mpi_datatypes"
  !
continue
  !
  call _DEBUG_WRITE_TIMER_(entering_procedure,pname)
  !
  ! Allocate the local arrays
  !
  allocate ( og_offset(1:ncell) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"og_offset",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( pts_offset(1:ncell) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"pts_offset",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( pts_blklen(1:ncell) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"pts_blklen",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the cell offsets as well as the block lengths
  ! and offsets for the solution data
  !
  do local_cell = 1,ncell
    !
    ! Indices of the data for this cell in the solution array
    !
    n1 = cell(local_cell)%beg_sp
    n2 = cell(local_cell)%end_sp
    !
    pts_blklen(local_cell) = int( n2-n1+1 , kind=int_mpi )
    !
    ! Get the global cell index for this local cell
    !
    global_cell = cell_map(mypnum+1)%loc_to_glb(local_cell)
    !
    ! Store the global cell index for the ordered_geom offset
    !
    og_offset(local_cell) = int( global_cell-1 , kind=int_mpi )
    !
    ! Starting point of the solution data for this global cell
    ! NOTE: The value global_pinc_ptr(global_cell)+1 is actually the starting
    !       solution point for this cell. The offsets for MPI start at 0 so
    !       dropping the +1 is the actual offset that we want.
    !
    pts_offset(local_cell) = int( global_pinc_ptr(global_cell) , kind=int_mpi )
    !
  end do
  !
  ! Allocate the cell_output_order array and
  ! the mask array for sorting the cell offsets
  !
  allocate ( cell_output_order(1:ncell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_output_order",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! Sort the og_offset array to find the output order for the cells so that
  ! the offsets for the MPI datatypes are monotonically increasing.
  !
  call heap_sort_integer(ncell,og_offset,cell_output_order)
 !!
 !allocate ( msk(1:ncell) , source=true , stat=ierr , errmsg=error_message )
 !call alloc_error(pname,"msk",1,__LINE__,__FILE__,ierr,error_message)
 !!
 !! Sort the og_offset array to find the output order for the cells so that
 !! the offsets for the MPI datatypes are monotonically increasing.
 !! NOTE: This array is assumed to be small enough so that it doesnt require
 !!       some fancy sorting algorithm and this simple minloc method will
 !!       more than likely be relatively inexpensive.
 !!
 !do nc = 1,ncell
 !  n = minloc( og_offset , dim=1 , mask=msk )
 !  msk(n) = fals
 !  cell_output_order(nc) = n
 !end do
 !do nc = 1,ncell
 !  if (cell_output_order(nc) /= nc) then
 !    write (iout,1) mypnum
 !    write (mypnum+1000,2) (n,cell_output_order(n),n=1,ncell)
 !    exit
 !  end if
 !end do
  1 format (" CPU #-",i4.4," : cell_output_order /= [1,ncell]")
  2 format ("cell_output_order(",i8,") = ",i0)
  !
  ! Deallocate the msk array
  !
  if (allocated(msk)) then
    deallocate ( msk , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Allocate the temporary array for applying the sorting to the offset
  ! and block length arrays.
  !
  allocate ( tmp_copy(1:ncell) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_copy",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Apply the sort to the og_offset array
  !
  tmp_copy = og_offset
  og_offset = tmp_copy(cell_output_order)
  !
  ! Apply the sort to the pts_offset array
  !
  tmp_copy = pts_offset
  pts_offset = tmp_copy(cell_output_order)
  !
  ! Apply the sort to the pts_blklen array
  !
  tmp_copy = pts_blklen
  pts_blklen = tmp_copy(cell_output_order)
  !
  ! Deallocate the temporary array
  !
  deallocate ( tmp_copy , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_copy",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the MPI derived datatype for writing the
  ! ordered_geom array to the restart file
  !
  call mpi_type_create_indexed_block(size(og_offset,kind=int_mpi), &
                                     1_int_mpi,og_offset, &
                                     MPI_INTEGER,ordered_geom_type,mpierr)
  call mpi_type_commit(ordered_geom_type,mpierr)
  !
  ! Create the MPI derived datatype for writing the
  ! solution data to the restart file
  !
  ! Multiply the offsets and block lengths by nq to account
  ! for the number of equations in the solution data
  !
  sol_offset = pts_offset(:) * int(nq,kind=int_mpi) ! F2003 AUTO-REALLOCATION
  sol_blklen = pts_blklen(:) * int(nq,kind=int_mpi) ! F2003 AUTO-REALLOCATION
  !
  call mpi_type_indexed(size(sol_offset,kind=int_mpi), &
                        sol_blklen,sol_offset, &
                        mpi_flttyp,solution_type,mpierr)
  call mpi_type_commit(solution_type,mpierr)
  !
  ! Deallocate the local arrays before leaving
  !
  if (allocated(og_offset)) then
    deallocate ( og_offset , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"og_offset",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(sol_offset)) then
    deallocate ( sol_offset , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"sol_offset",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(sol_blklen)) then
    deallocate ( sol_blklen , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"sol_blklen",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call _DEBUG_WRITE_TIMER_(leaving_procedure,pname)
  !
end subroutine get_output_mpi_datatypes
!
!###############################################################################
!
subroutine get_host_time_ave_mpi_datatypes()
  !
  !.. Use Statement ..
  use flowvar_mod, only : uavesp
  !
  !.. Local Scalars ..
  integer :: naq,ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_host_time_ave_mpi_datatypes"
  !
continue
  !
  call _DEBUG_WRITE_TIMER_(entering_procedure,pname)
  !
  ! Create the cell offsets as well as the block lengths
  ! and offsets for the solution data
  !
  if (i_am_host_root) then
    !
    ! Create the MPI derived datatype for writing the
    ! time averaged solution data to the averaging file
    !
    if (allocated(uavesp)) then
      !
      naq = minval( shape(uavesp) )
      !
      call mpi_type_contiguous(int(naq,kind=int_mpi),mpi_flttyp, &
                               ave_block_type,mpierr)
      call mpi_type_commit(ave_block_type,mpierr)
      !
      !
      !
      call mpi_type_indexed(size(pts_offset,kind=int_mpi), &
                            pts_blklen,pts_offset, &
                            ave_block_type,averaging_type,mpierr)
      call mpi_type_commit(averaging_type,mpierr)
      !
    end if
    !
  end if
  !
  ! Deallocate the local arrays before leaving
  !
  if (allocated(pts_offset)) then
    deallocate ( pts_offset , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"pts_offset",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(pts_blklen)) then
    deallocate ( pts_blklen , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"pts_blklen",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call _DEBUG_WRITE_TIMER_(leaving_procedure,pname)
  !
end subroutine get_host_time_ave_mpi_datatypes
!
!###############################################################################
!
subroutine get_time_ave_mpi_datatypes()
  !
  !.. Use Statement ..
  use flowvar_mod, only : uavesp
  !
  !.. Local Scalars ..
  integer :: naq,ierr
  !
  !.. Local Allocatable Arrays ..
  logical(lk),      allocatable :: msk(:)
  integer(int_mpi), allocatable :: ave_offset(:)
  integer(int_mpi), allocatable :: ave_blklen(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_time_ave_mpi_datatypes"
  !
continue
  !
  call _DEBUG_WRITE_TIMER_(entering_procedure,pname)
  !
  ! Create the MPI derived datatype for writing the
  ! time averaged solution data to the averaging file
  !
  if (allocated(uavesp)) then
    !
    naq = minval( shape(uavesp) )
    !
    ! Multiply the offsets and block lengths by naq to account for
    ! the number of time averaged quantities in the solution data
    !
    ave_offset = pts_offset(:)*int(naq,kind=int_mpi) ! F2003 AUTO-REALLOCATION
    ave_blklen = pts_blklen(:)*int(naq,kind=int_mpi) ! F2003 AUTO-REALLOCATION
    !
    call mpi_type_indexed(size(ave_offset,kind=int_mpi), &
                          ave_blklen,ave_offset, &
                          mpi_flttyp,averaging_type,mpierr)
    call mpi_type_commit(averaging_type,mpierr)
    !
  end if
  !
  ! Deallocate the local arrays before leaving
  !
  if (allocated(pts_offset)) then
    deallocate ( pts_offset , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"pts_offset",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(pts_blklen)) then
    deallocate ( pts_blklen , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"pts_blklen",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(ave_offset)) then
    deallocate ( ave_offset , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ave_offset",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(ave_blklen)) then
    deallocate ( ave_blklen , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ave_blklen",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call _DEBUG_WRITE_TIMER_(leaving_procedure,pname)
  !
end subroutine get_time_ave_mpi_datatypes
!
!###############################################################################
!
subroutine get_input_mpi_datatypes(restart_nq,restart_rp,real_type, &
                                   read_block_type,read_type, &
                                   num_pts,num_global_pts, &
                                   ordered_geom,cell_input_data)
  !
  !.. Use Statements ..
  use geomvar_mod, only : ncell ! ncell is recomputed after cell_map is created
  use geomvar_mod, only : cell_map
  !
  !.. Formal Arguments ..
  integer(wip),                     intent(in) :: restart_nq
  integer(wip),                     intent(in) :: restart_rp
  _MPI_DATA_TYPE_,               intent(inout) :: real_type
  _MPI_DATA_TYPE_,               intent(inout) :: read_block_type
  _MPI_DATA_TYPE_,               intent(inout) :: read_type
  integer(wip),                  intent(inout) :: num_pts
  integer(wip),                  intent(inout) :: num_global_pts
  integer(ogp),                     intent(in) :: ordered_geom(:)
  type(rstrt_cell_t), allocatable, intent(out) :: cell_input_data(:)
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,ierr,og
  integer :: local_cell,global_cell
  integer :: this_cell,this_geom,this_order
  integer :: inp_cell,inp_geom,inp_order
  !
  !.. Local Allocatable Arrays ..
  logical(lk),      allocatable :: msk(:)
  integer(int_mpi), allocatable :: sol_offset(:)
  integer(int_mpi), allocatable :: sol_blklen(:)
  integer(int_mpi), allocatable :: tmp_copy(:)
  integer,          allocatable :: pinc_ptr(:)
  integer,          allocatable :: cell_input_order(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_input_mpi_datatypes"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  read_block_type = MPI_DATATYPE_NULL
  read_type = MPI_DATATYPE_NULL
  !
  ! First determine the MPI real type for the restart solution data
  !
  if (restart_rp == wp) then
    real_type = mpi_flttyp
  else
    select case (restart_rp)
      case (real64)  ; real_type = MPI_REAL8
      case (real32)  ; real_type = MPI_REAL4
      case (real128) ; real_type = MPI_REAL16
      case default   ; real_type = MPI_REAL
    end select
  end if
  !
  ! Allocate the local arrays
  !
  allocate ( pinc_ptr(1:size(ordered_geom)+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"pinc_ptr",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( sol_offset(1:ncell) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sol_offset",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( sol_blklen(1:ncell) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sol_blklen",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create a pointer array for the solution points of each cell
  ! corresponding to how the cells are stored in the restart file.
  !
  do global_cell = 1,size(ordered_geom)
    !
    og = int(ordered_geom(global_cell),kind=kind(og))
    this_order = og / 100
    this_geom = og - 100*this_order
    !
    pinc_ptr(global_cell+1) = pinc_ptr(global_cell) + &
                              Cell_Solpts(this_geom,this_order)
    !
  end do
  !
  ! Get the total number of points in the restart file
  !
  num_global_pts = pinc_ptr(size(ordered_geom)+1)
  !
  ! Create the offset and block length arrays for the solution data
  !
  do local_cell = 1,ncell
    !
    ! Get the global cell index for this local cell
    !
    global_cell = cell_map(mypnum+1)%loc_to_glb(local_cell)
    !
    ! For now, store the value in ordered_geom within sol_blklen
    !
    sol_blklen(local_cell) = int( ordered_geom(global_cell) , kind=int_mpi )
    !
    ! Starting point of the solution data for this global cell
    ! NOTE: The value global_pinc_ptr(global_cell)+1 is actually the starting
    !       solution point for this cell. The offsets for MPI start at 0 so
    !       dropping the +1 is the actual offset that we want.
    !
    sol_offset(local_cell) = int( pinc_ptr(global_cell) , kind=int_mpi )
    !
  end do
  !
  ! Deallocate the pinc_ptr since it is no longer needed
  !
  deallocate ( pinc_ptr , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"pinc_ptr",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate the cell_input_data array and
  ! the mask array for sorting the cell offsets
  !
  if (allocated(cell_input_data)) then
    deallocate ( cell_input_data , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_input_data",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  allocate ( cell_input_data(1:ncell) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_input_data",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  allocate ( cell_input_order(1:ncell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_input_order",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  allocate ( msk(1:ncell) , source=true , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Sort the sol_offset array to find the output order for the cells so that
  ! the offsets for the MPI datatypes are monotonically increasing.
  ! NOTE: This array is assumed to be small enough so that it doesnt require
  !       some fancy sorting algorithm and this simple minloc method will
  !       more than likely be relatively inexpensive.
  !
  call heap_sort_integer(ncell,sol_offset,cell_input_order)
  !
  do nc = 1,ncell
    cell_input_data(nc)%input_cell = cell_input_order(nc)
  end do
  !
  deallocate ( cell_input_order , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_input_order",2,__LINE__,__FILE__, &
                   ierr,error_message)
  !
 !do nc = 1,ncell
 !  n = minloc( sol_offset , dim=1 , mask=msk )
 !  msk(n) = fals
 !  cell_input_data(nc)%input_cell = n
 !end do
  !
  ! Deallocate the msk array
  !
  deallocate ( msk , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate the temporary array for applying the sorting to the offset
  ! and block length arrays.
  !
  allocate ( tmp_copy(1:ncell) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_copy",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Apply the sort to the sol_offset array
  !
  tmp_copy(:) = sol_offset
  sol_offset(:) = tmp_copy(cell_input_data(:)%input_cell)
  !
  ! Apply the sort to the sol_blklen array
  !
  tmp_copy(:) = sol_blklen
  sol_blklen(:) = tmp_copy(cell_input_data(:)%input_cell)
  !
  ! Deallocate the temporary array
  !
  deallocate ( tmp_copy , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_copy",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the remaining information within cell_input_data and convert
  ! the ordered_geom values that are currently stored in the block length
  ! array to the actual block lengths for the restart data
  !
  do this_cell = 1,ncell
    !
    inp_cell = cell_input_data(this_cell)%input_cell
    !
    ! Extract the order and geometry for this cell from the ordered_geom
    ! value that was stored in sol_blklen.
    !
    inp_order = sol_blklen(inp_cell) / 100
    inp_geom = sol_blklen(inp_cell) - inp_order*100
    !
    ! Compute the actual block length for this cell
    !
    sol_blklen(inp_cell) = int( Cell_Solpts(inp_geom,inp_order) , kind=int_mpi )
    !
    ! Store all this information in cell_input_data for this_cell
    !
    cell_input_data(this_cell)%input_geom   = inp_geom
    cell_input_data(this_cell)%input_order  = inp_order
    !
  end do
  !
  ! Get the beginning and ending solution points for each cell relative
  ! to the ordering the cells are in within the restart file.
  ! NOTE: We have to do this as a separate loop from above so that the
  !       information in sol_blklen has been corrected for all cells.
  !
  do this_cell = 1,ncell
    !
    inp_cell = cell_input_data(this_cell)%input_cell
    !
    n1 = sum( sol_blklen(1:inp_cell-1) ) + 1
    n2 = n1-1 + sol_blklen(inp_cell)
    !
    cell_input_data(this_cell)%input_beg_sp = n1
    cell_input_data(this_cell)%input_end_sp = n2
    !
  end do
  !
  ! Get the total number of solution points that will be read from the
  ! restart file for this processor
  !
  num_pts = int( sum(sol_blklen) , kind=kind(num_pts) )
  !
  ! Dont forget to account for the number of equations in the offsets
  ! and block lengths for the solution data
  !
  call mpi_type_contiguous(int(restart_nq,kind=int_mpi),real_type, &
                           read_block_type,mpierr)
  call mpi_type_commit(read_block_type,mpierr)
  !
  ! Create the MPI derived datatype for reading the
  ! solution data to the restart file
  !
  call mpi_type_indexed(size(sol_offset,kind=int_mpi), &
                        sol_blklen,sol_offset, &
                        read_block_type,read_type,mpierr)
  call mpi_type_commit(read_type,mpierr)
  !
  ! Deallocate the local arrays before leaving
  !
  deallocate ( sol_offset , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sol_offset",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( sol_blklen , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sol_blklen",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  !
end subroutine get_input_mpi_datatypes
!
!###############################################################################
!
subroutine check_restart_solution_form(restart_solution_form,ordered_geom, &
                                       rst_form_matches_sol_form)
  !
  !.. Use Statements ..
  use ovar,        only : loc_solution_pts
  use geomvar_mod, only : global_cell_geom
  use geomvar_mod, only : global_cell_order
  !
  !.. Formal Arguments ..
  integer,        intent(in) :: restart_solution_form
  integer(ogp),   intent(in) :: ordered_geom(:)
  logical(lk), intent(inout) :: rst_form_matches_sol_form
  !
  !.. Local Scalars ..
  integer :: nc
  integer :: rst_geom,rst_order
  integer :: sol_geom,sol_order
  !
  logical(lk) :: found_geom_mismatch
  logical(lk) :: found_order_mismatch
  logical(lk) :: requires_interpolation
  logical(lk) :: requires_projection
  !
  character(len=10) :: inp_solpts
  character(len=50) :: restart_form
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_restart_solution_form"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  requires_interpolation = fals
  requires_projection = fals
  !
  rst_form_matches_sol_form = true
  !
  ! Check each of the restart cells to see if they
  ! match the cells for the current simulation.
  !
  found_geom_mismatch = fals
  found_order_mismatch = fals
  !
  do nc = 1,size(ordered_geom)
    !
    rst_order = int(ordered_geom(nc),kind=kind(rst_order)) / 100
    rst_geom = int(ordered_geom(nc),kind=kind(rst_geom)) - 100*rst_order
    !
    sol_order = global_cell_order(nc)
    sol_geom = global_cell_geom(nc)
    !
    if (sol_order /= rst_order) then
      found_order_mismatch = true
      if (rst_order > sol_order) then
        requires_projection = true
      else
        requires_interpolation = true
      end if
    end if
    !
    if (sol_geom /= rst_geom) found_geom_mismatch = true
    !
  end do
  !
  ! Stop execution if there is a geometry mismatch which means
  ! the grid file and restart file are not compatible.
  ! 
  if (found_geom_mismatch) then
    write (error_message,100)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Change rst_form_matches_sol_form if an order mismatch is found
  !
  if (found_order_mismatch) rst_form_matches_sol_form = fals
  !
  ! Check to see if the restart solution is stored
  ! in a different form than the current solution
  !
  if (restart_solution_form /= -loc_solution_pts) then
    !
    rst_form_matches_sol_form = fals
    !
    select case (loc_solution_pts)
      case (Legendre_Gauss)
        inp_solpts = "Gauss"
      case (Legendre_Gauss_Lobatto)
        inp_solpts = "Lobatto"
      case default
        write (error_message,101)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end select
    !
    select case (restart_solution_form)
      case (restart_gauss_nodal_solution)
        restart_form = "nodal form for the Gauss solution points"
      case (restart_lobatto_nodal_solution)
        restart_form = "nodal form for the Lobatto solution points"
      case (restart_legendre_modal_solution)
        write (error_message,102) "Legendre"
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      case (restart_monomial_modal_solution)
        write (error_message,102) "monomial"
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      case default
        write (error_message,103)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end select
    !
    if (mypnum == glb_root) then
      write (iout,1) trim(adjustl(restart_form)),trim(adjustl(inp_solpts))
      flush (iout)
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format(/,1x,93("#"),/, &
             40x,"!!! WARNING !!!",/, &
             2x,"The format of the solution data in the restart file does ", &
                "not match the current simulation.",/, &
             6x,"The restart data is in ",a,", and",/, &
             6x,"the current simulation requires the data to be defined ", &
                "at the ",a," points.",/, &
             2x,"Attempting to convert the restart solution into the form ", &
                "needed for this simulation.",/, &
             1x,93("#"),/)
  !
  100 format ("The geometry of one or more cells from the restart file ", &
              "does not match that of the current simulation! It appears ", &
              "that the current grid file is not compatible with the ", &
              "specified restart file!")
  101 format ("The solution for the current simulation seems to be in ", &
              "an unknown form when trying to read in the solution data ", &
              "from the restart file!")
  102 format ("The restart file appears to contain solution data stored ", &
              "in modal form using ",a," basis functions. Restarting ", &
              "from modal solution data is not supported at this time.")
  103 format ("The solution data in the restart file seems to be in an ", &
              "unknown form!")
  !
end subroutine check_restart_solution_form
!
!###############################################################################
!###############################################################################
! MPI read subroutine for reading 8-bit INTEGER data into a 32-bit array.
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_int8_data_i4(file_handle,int_data,mpistat)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  integer(i4),  intent(inout) :: int_data(:)
  !
  !.. Optional Arguments ..
  _MPI_STATUS_TYPE_, optional, intent(inout) :: mpistat(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
  !.. Local Allocatable Arrays ..
  integer(int8), allocatable :: read_int(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_int8_data_i4"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(int_data)
  !
  allocate ( read_int(1:np) , stat=ierr )
  !
  if (present(mpistat)) then
    call mpi_file_read_all(file_handle,read_int,int(np,kind=int_mpi), &
                           MPI_INTEGER1,mpistat(1),mpierr)
  else
    call mpi_file_read_all(file_handle,read_int,int(np,kind=int_mpi), &
                           MPI_INTEGER1,MPI_STATUS_IGNORE,mpierr)
  end if
  !
  int_data(:) = int( read_int(:) , kind=kind(int_data) )
  !
  deallocate ( read_int , stat=ierr )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_int8_data_i4
!
!###############################################################################
!###############################################################################
! MPI read subroutine for reading 16-bit INTEGER data into a 32-bit array.
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_int16_data_i4(file_handle,int_data,mpistat)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  integer(i4),  intent(inout) :: int_data(:)
  !
  !.. Optional Arguments ..
  _MPI_STATUS_TYPE_, optional, intent(inout) :: mpistat(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
  !.. Local Allocatable Arrays ..
  integer(int16), allocatable :: read_int(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_int16_data_i4"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(int_data)
  !
  allocate ( read_int(1:np) , stat=ierr )
  !
  if (present(mpistat)) then
    call mpi_file_read_all(file_handle,read_int,int(np,kind=int_mpi), &
                           MPI_INTEGER2,mpistat(1),mpierr)
  else
    call mpi_file_read_all(file_handle,read_int,int(np,kind=int_mpi), &
                           MPI_INTEGER2,MPI_STATUS_IGNORE,mpierr)
  end if
  !
  int_data(:) = int( read_int(:) , kind=kind(int_data) )
  !
  deallocate ( read_int , stat=ierr )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_int16_data_i4
!
!###############################################################################
!###############################################################################
! MPI read subroutine for reading 32-bit INTEGER data into a 32-bit array.
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_int32_data_i4(file_handle,int_data,mpistat)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  integer(i4),  intent(inout) :: int_data(:)
  !
  !.. Optional Arguments ..
  _MPI_STATUS_TYPE_, optional, intent(inout) :: mpistat(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
  !.. Local Allocatable Arrays ..
  integer(int32), allocatable :: read_int(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_int32_data_i4"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(int_data)
  !
  allocate ( read_int(1:np) , stat=ierr )
  !
  if (present(mpistat)) then
    call mpi_file_read_all(file_handle,read_int,int(np,kind=int_mpi), &
                           MPI_INTEGER4,mpistat(1),mpierr)
  else
    call mpi_file_read_all(file_handle,read_int,int(np,kind=int_mpi), &
                           MPI_INTEGER4,MPI_STATUS_IGNORE,mpierr)
  end if
  !
  int_data(:) = int( read_int(:) , kind=kind(int_data) )
  !
  deallocate ( read_int , stat=ierr )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_int32_data_i4
!
!###############################################################################
!###############################################################################
! MPI read subroutine for reading 64-bit INTEGER data into a 32-bit array.
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_int64_data_i4(file_handle,int_data,mpistat)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  integer(i4),  intent(inout) :: int_data(:)
  !
  !.. Optional Arguments ..
  _MPI_STATUS_TYPE_, optional, intent(inout) :: mpistat(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
  !.. Local Allocatable Arrays ..
  integer(int64), allocatable :: read_int(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_int64_data_i4"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(int_data)
  !
  allocate ( read_int(1:np) , stat=ierr )
  !
  if (present(mpistat)) then
    call mpi_file_read_all(file_handle,read_int,int(np,kind=int_mpi), &
                           MPI_INTEGER8,mpistat(1),mpierr)
  else
    call mpi_file_read_all(file_handle,read_int,int(np,kind=int_mpi), &
                           MPI_INTEGER8,MPI_STATUS_IGNORE,mpierr)
  end if
  !
  int_data(:) = int( read_int(:) , kind=kind(int_data) )
  !
  deallocate ( read_int , stat=ierr )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_int64_data_i4
!
!###############################################################################
!###############################################################################
! MPI read subroutine for reading 8-bit INTEGER data into a 64-bit array.
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_int8_data_i8(file_handle,int_data,mpistat)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  integer(i8),  intent(inout) :: int_data(:)
  !
  !.. Optional Arguments ..
  _MPI_STATUS_TYPE_, optional, intent(inout) :: mpistat(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
  !.. Local Allocatable Arrays ..
  integer(int8), allocatable :: read_int(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_int8_data_i8"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(int_data)
  !
  allocate ( read_int(1:np) , stat=ierr )
  !
  if (present(mpistat)) then
    call mpi_file_read_all(file_handle,read_int,int(np,kind=int_mpi), &
                           MPI_INTEGER1,mpistat(1),mpierr)
  else
    call mpi_file_read_all(file_handle,read_int,int(np,kind=int_mpi), &
                           MPI_INTEGER1,MPI_STATUS_IGNORE,mpierr)
  end if
  !
  int_data(:) = int( read_int(:) , kind=kind(int_data) )
  !
  deallocate ( read_int , stat=ierr )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_int8_data_i8
!
!###############################################################################
!###############################################################################
! MPI read subroutine for reading 16-bit INTEGER data into a 64-bit array.
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_int16_data_i8(file_handle,int_data,mpistat)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  integer(i8),  intent(inout) :: int_data(:)
  !
  !.. Optional Arguments ..
  _MPI_STATUS_TYPE_, optional, intent(inout) :: mpistat(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
  !.. Local Allocatable Arrays ..
  integer(int16), allocatable :: read_int(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_int16_data_i8"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(int_data)
  !
  allocate ( read_int(1:np) , stat=ierr )
  !
  if (present(mpistat)) then
    call mpi_file_read_all(file_handle,read_int,int(np,kind=int_mpi), &
                           MPI_INTEGER2,mpistat(1),mpierr)
  else
    call mpi_file_read_all(file_handle,read_int,int(np,kind=int_mpi), &
                           MPI_INTEGER2,MPI_STATUS_IGNORE,mpierr)
  end if
  !
  int_data(:) = int( read_int(:) , kind=kind(int_data) )
  !
  deallocate ( read_int , stat=ierr )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_int16_data_i8
!
!###############################################################################
!###############################################################################
! MPI read subroutine for reading 32-bit INTEGER data into a 64-bit array.
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_int32_data_i8(file_handle,int_data,mpistat)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  integer(i8),  intent(inout) :: int_data(:)
  !
  !.. Optional Arguments ..
  _MPI_STATUS_TYPE_, optional, intent(inout) :: mpistat(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
  !.. Local Allocatable Arrays ..
  integer(int32), allocatable :: read_int(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_int32_data_i8"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(int_data)
  !
  allocate ( read_int(1:np) , stat=ierr )
  !
  if (present(mpistat)) then
    call mpi_file_read_all(file_handle,read_int,int(np,kind=int_mpi), &
                           MPI_INTEGER4,mpistat(1),mpierr)
  else
    call mpi_file_read_all(file_handle,read_int,int(np,kind=int_mpi), &
                           MPI_INTEGER4,MPI_STATUS_IGNORE,mpierr)
  end if
  !
  int_data(:) = int( read_int(:) , kind=kind(int_data) )
  !
  deallocate ( read_int , stat=ierr )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_int32_data_i8
!
!###############################################################################
!###############################################################################
! MPI read subroutine for reading 64-bit INTEGER data into a 64-bit array.
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_int64_data_i8(file_handle,int_data,mpistat)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  integer(i8),  intent(inout) :: int_data(:)
  !
  !.. Optional Arguments ..
  _MPI_STATUS_TYPE_, optional, intent(inout) :: mpistat(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
  !.. Local Allocatable Arrays ..
  integer(int64), allocatable :: read_int(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_int64_data_i8"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(int_data)
  !
  allocate ( read_int(1:np) , stat=ierr )
  !
  if (present(mpistat)) then
    call mpi_file_read_all(file_handle,read_int,int(np,kind=int_mpi), &
                           MPI_INTEGER8,mpistat(1),mpierr)
  else
    call mpi_file_read_all(file_handle,read_int,int(np,kind=int_mpi), &
                           MPI_INTEGER8,MPI_STATUS_IGNORE,mpierr)
  end if
  !
  int_data(:) = int( read_int(:) , kind=kind(int_data) )
  !
  deallocate ( read_int , stat=ierr )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_int64_data_i8
!
!###############################################################################
!###############################################################################
! MPI read subroutines for 32-bit REAL data
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_r4_rank1(file_handle,real_data,mpistat)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  real(wp),     intent(inout) :: real_data(:)
  !
  !.. Optional Arguments ..
  _MPI_STATUS_TYPE_, optional, intent(inout) :: mpistat(:)
  !
  !.. Local Arrays ..
  real(real32) :: read_real(1:size(real_data))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_r4_rank1"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (present(mpistat)) then
    call mpi_file_read_all(file_handle,read_real, &
                           size(read_real,kind=int_mpi), &
                           MPI_REAL4,mpistat(1),mpierr)
  else
    call mpi_file_read_all(file_handle,read_real, &
                           size(read_real,kind=int_mpi), &
                           MPI_REAL4,MPI_STATUS_IGNORE,mpierr)
  end if
  !
  real_data = real( read_real , kind=wp )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_r4_rank1
!-------------------------------------------------------------------------------
subroutine read_rstrt_file_r4_rank2(file_handle,real_data, &
                                    read_block_type,mpistat)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  real(wp),     intent(inout) :: real_data(:,:)
  !
  !.. Optional Arguments ..
  _MPI_DATA_TYPE_,   optional,    intent(in) :: read_block_type
  _MPI_STATUS_TYPE_, optional, intent(inout) :: mpistat(:)
  !
  !.. Local Arrays ..
  real(real32) :: read_real(1:size(real_data,dim=1), &
                            1:size(real_data,dim=2))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_r4_rank2"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (present(mpistat)) then
    if (present(read_block_type)) then
      call mpi_file_read_all(file_handle,read_real, &
                             size(read_real,dim=2,kind=int_mpi), &
                             read_block_type,mpistat(1),mpierr)
    else
      call mpi_file_read_all(file_handle,read_real, &
                             size(read_real,kind=int_mpi), &
                             MPI_REAL4,mpistat(1),mpierr)
    end if
  else
    if (present(read_block_type)) then
      call mpi_file_read_all(file_handle,read_real, &
                             size(read_real,dim=2,kind=int_mpi), &
                             read_block_type,MPI_STATUS_IGNORE,mpierr)
    else
      call mpi_file_read_all(file_handle,read_real, &
                             size(read_real,kind=int_mpi), &
                             MPI_REAL4,MPI_STATUS_IGNORE,mpierr)
    end if
  end if
  !
  real_data = real( read_real , kind=wp )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_r4_rank2
!
!###############################################################################
!###############################################################################
! MPI read subroutines for 64-bit REAL data
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_r8_rank1(file_handle,real_data,mpistat)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  real(wp),     intent(inout) :: real_data(:)
  !
  !.. Optional Arguments ..
  _MPI_STATUS_TYPE_, optional, intent(inout) :: mpistat(:)
  !
  !.. Local Arrays ..
  real(real64) :: read_real(1:size(real_data))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_r8_rank1"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (present(mpistat)) then
    call mpi_file_read_all(file_handle,read_real, &
                           size(read_real,kind=int_mpi), &
                           MPI_REAL8,mpistat(1),mpierr)
  else
    call mpi_file_read_all(file_handle,read_real, &
                           size(read_real,kind=int_mpi), &
                           MPI_REAL8,MPI_STATUS_IGNORE,mpierr)
  end if
  !
  real_data = real( read_real , kind=wp )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_r8_rank1
!-------------------------------------------------------------------------------
subroutine read_rstrt_file_r8_rank2(file_handle,real_data, &
                                    read_block_type,mpistat)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  real(wp),     intent(inout) :: real_data(:,:)
  !
  !.. Optional Arguments ..
  _MPI_DATA_TYPE_,   optional,    intent(in) :: read_block_type
  _MPI_STATUS_TYPE_, optional, intent(inout) :: mpistat(:)
  !
  !.. Local Arrays ..
  real(real64) :: read_real(1:size(real_data,dim=1), &
                            1:size(real_data,dim=2))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_r8_rank2"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (present(mpistat)) then
    if (present(read_block_type)) then
      call mpi_file_read_all(file_handle,read_real, &
                             size(read_real,dim=2,kind=int_mpi), &
                             read_block_type,mpistat(1),mpierr)
    else
      call mpi_file_read_all(file_handle,read_real, &
                             size(read_real,kind=int_mpi), &
                             MPI_REAL8,mpistat(1),mpierr)
    end if
    if (present(read_block_type)) then
      call mpi_file_read_all(file_handle,read_real, &
                             size(read_real,dim=2,kind=int_mpi), &
                             read_block_type,MPI_STATUS_IGNORE,mpierr)
    else
      call mpi_file_read_all(file_handle,read_real, &
                             size(read_real,kind=int_mpi), &
                             MPI_REAL8,MPI_STATUS_IGNORE,mpierr)
    end if
  end if
  !
  real_data = real( read_real , kind=wp )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_r8_rank2
!
!###############################################################################
!###############################################################################
! MPI read subroutines for 128-bit REAL data
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_r16_rank1(file_handle,real_data,mpistat)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  real(wp),     intent(inout) :: real_data(:)
  !
  !.. Optional Arguments ..
  _MPI_STATUS_TYPE_, optional, intent(inout) :: mpistat(:)
  !
  !.. Local Arrays ..
  real(max(real64,real128)) :: read_real(1:size(real_data))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_r16_rank1"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (present(mpistat)) then
    call mpi_file_read_all(file_handle,read_real, &
                           size(read_real,kind=int_mpi), &
                           MPI_REAL16,mpistat(1),mpierr)
  else
    call mpi_file_read_all(file_handle,read_real, &
                           size(read_real,kind=int_mpi), &
                           MPI_REAL16,MPI_STATUS_IGNORE,mpierr)
  end if
  !
  real_data = real( read_real , kind=wp )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_r16_rank1
!-------------------------------------------------------------------------------
subroutine read_rstrt_file_r16_rank2(file_handle,real_data, &
                                     read_block_type,mpistat)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  real(wp),     intent(inout) :: real_data(:,:)
  !
  !.. Optional Arguments ..
  _MPI_DATA_TYPE_,   optional,    intent(in) :: read_block_type
  _MPI_STATUS_TYPE_, optional, intent(inout) :: mpistat(:)
  !
  !.. Local Arrays ..
  real(max(real64,real128)) :: read_real(1:size(real_data,dim=1), &
                                         1:size(real_data,dim=2))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_r16_rank2"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (present(mpistat)) then
    if (present(read_block_type)) then
      call mpi_file_read_all(file_handle,read_real, &
                             size(read_real,dim=2,kind=int_mpi), &
                             read_block_type,mpistat(1),mpierr)
    else
      call mpi_file_read_all(file_handle,read_real, &
                             size(read_real,kind=int_mpi), &
                             MPI_REAL16,mpistat(1),mpierr)
    end if
    if (present(read_block_type)) then
      call mpi_file_read_all(file_handle,read_real, &
                             size(read_real,dim=2,kind=int_mpi), &
                             read_block_type,MPI_STATUS_IGNORE,mpierr)
    else
      call mpi_file_read_all(file_handle,read_real, &
                             size(read_real,kind=int_mpi), &
                             MPI_REAL16,MPI_STATUS_IGNORE,mpierr)
    end if
  end if
  !
  real_data = real( read_real , kind=wp )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_r16_rank2
!
!###############################################################################
!###############################################################################
! MPI read subroutines for current working precision REAL data (wp)
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_wp_rank1(file_handle,real_data,mpistat)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  real(wp),     intent(inout) :: real_data(:)
  !
  !.. Optional Arguments ..
  _MPI_STATUS_TYPE_, optional, intent(inout) :: mpistat(:)
  !
  !.. Local Arrays ..
  real(wp) :: read_real(1:size(real_data))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_wp_rank1"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (present(mpistat)) then
    call mpi_file_read_all(file_handle,read_real, &
                           size(real_data,kind=int_mpi), &
                           mpi_flttyp,mpistat(1),mpierr)
  else
    call mpi_file_read_all(file_handle,read_real, &
                           size(real_data,kind=int_mpi), &
                           mpi_flttyp,MPI_STATUS_IGNORE,mpierr)
  end if
  !
  real_data = read_real
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_wp_rank1
!-------------------------------------------------------------------------------
subroutine read_rstrt_file_wp_rank2(file_handle,real_data, &
                                    read_block_type,mpistat)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  real(wp),     intent(inout) :: real_data(:,:)
  !
  !.. Optional Arguments ..
  _MPI_DATA_TYPE_,   optional,    intent(in) :: read_block_type
  _MPI_STATUS_TYPE_, optional, intent(inout) :: mpistat(:)
  !
  !.. Local Arrays ..
  real(wp) :: read_real(1:size(real_data,dim=1), &
                        1:size(real_data,dim=2))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_wp_rank2"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (present(mpistat)) then
    if (present(read_block_type)) then
      call mpi_file_read_all(file_handle,read_real, &
                             size(real_data,dim=2,kind=int_mpi), &
                             read_block_type,mpistat(1),mpierr)
    else
      call mpi_file_read_all(file_handle,read_real, &
                             size(real_data,kind=int_mpi), &
                             mpi_flttyp,mpistat(1),mpierr)
    end if
  else
    if (present(read_block_type)) then
      call mpi_file_read_all(file_handle,read_real, &
                             size(real_data,dim=2,kind=int_mpi), &
                             read_block_type,MPI_STATUS_IGNORE,mpierr)
    else
      call mpi_file_read_all(file_handle,read_real, &
                             size(real_data,kind=int_mpi), &
                             mpi_flttyp,MPI_STATUS_IGNORE,mpierr)
    end if
  end if
  !
  real_data = read_real
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_wp_rank2
!
!###############################################################################
!###############################################################################
! Serial read subroutine for reading 8-bit INTEGER data into a 32-bit array.
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_serial_int8_data_i4(io,int_data)
  !
  !.. Formal Arguments ..
  integer,        intent(in) :: io
  integer(i4), intent(inout) :: int_data(:)
  !
  !.. Local Scalars ..
  integer :: i,ierr,np
  !
  !.. Local Allocatable Arrays ..
  integer(int8), allocatable :: read_int(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_serial_int8_data_i4"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(int_data)
  !
  allocate ( read_int(1:np) , stat=ierr )
  !
  read (io,iostat=ierr) (read_int(i), i=1,np)
  !
  do i = 1,np
    int_data(i) = int( read_int(i) , kind=kind(int_data) )
  end do
  !
  deallocate ( read_int , stat=ierr )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_serial_int8_data_i4
!
!###############################################################################
!###############################################################################
! Serial read subroutine for reading 16-bit INTEGER data into a 32-bit array.
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_serial_int16_data_i4(io,int_data)
  !
  !.. Formal Arguments ..
  integer,        intent(in) :: io
  integer(i4), intent(inout) :: int_data(:)
  !
  !.. Local Scalars ..
  integer :: i,ierr,np
  !
  !.. Local Allocatable Arrays ..
  integer(int16), allocatable :: read_int(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_serial_int16_data_i4"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(int_data)
  !
  allocate ( read_int(1:np) , stat=ierr )
  !
  read (io,iostat=ierr) (read_int(i), i=1,np)
  !
  do i = 1,np
    int_data(i) = int( read_int(i) , kind=kind(int_data) )
  end do
  !
  deallocate ( read_int , stat=ierr )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_serial_int16_data_i4
!
!###############################################################################
!###############################################################################
! Serial read subroutine for reading 32-bit INTEGER data into a 32-bit array.
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_serial_int32_data_i4(io,int_data)
  !
  !.. Formal Arguments ..
  integer,        intent(in) :: io
  integer(i4), intent(inout) :: int_data(:)
  !
  !.. Local Scalars ..
  integer :: i,ierr,np
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_serial_int32_data_i4"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(int_data)
  !
  read (io,iostat=ierr) (int_data(i), i=1,np)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_serial_int32_data_i4
!
!###############################################################################
!###############################################################################
! Serial read subroutine for reading 64-bit INTEGER data into a 32-bit array.
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_serial_int64_data_i4(io,int_data)
  !
  !.. Formal Arguments ..
  integer,        intent(in) :: io
  integer(i4), intent(inout) :: int_data(:)
  !
  !.. Local Scalars ..
  integer :: i,ierr,np
  !
  !.. Local Allocatable Arrays ..
  integer(int64), allocatable :: read_int(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_serial_int64_data_i4"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(int_data)
  !
  allocate ( read_int(1:np) , stat=ierr )
  !
  read (io,iostat=ierr) (read_int(i), i=1,np)
  !
  do i = 1,np
    int_data(i) = int( read_int(i) , kind=kind(int_data) )
  end do
  !
  deallocate ( read_int , stat=ierr )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_serial_int64_data_i4
!
!###############################################################################
!###############################################################################
! Serial read subroutine for reading 8-bit INTEGER data into a 64-bit array.
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_serial_int8_data_i8(io,int_data)
  !
  !.. Formal Arguments ..
  integer,        intent(in) :: io
  integer(i8), intent(inout) :: int_data(:)
  !
  !.. Local Scalars ..
  integer :: i,ierr,np
  !
  !.. Local Allocatable Arrays ..
  integer(int8), allocatable :: read_int(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_serial_int8_data_i8"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(int_data)
  !
  allocate ( read_int(1:np) , stat=ierr )
  !
  read (io,iostat=ierr) (read_int(i), i=1,np)
  !
  do i = 1,np
    int_data(i) = int( read_int(i) , kind=kind(int_data) )
  end do
  !
  deallocate ( read_int , stat=ierr )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_serial_int8_data_i8
!
!###############################################################################
!###############################################################################
! Serial read subroutine for reading 16-bit INTEGER data into a 64-bit array.
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_serial_int16_data_i8(io,int_data)
  !
  !.. Formal Arguments ..
  integer,        intent(in) :: io
  integer(i8), intent(inout) :: int_data(:)
  !
  !.. Local Scalars ..
  integer :: i,ierr,np
  !
  !.. Local Allocatable Arrays ..
  integer(int16), allocatable :: read_int(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_serial_int16_data_i8"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(int_data)
  !
  allocate ( read_int(1:np) , stat=ierr )
  !
  read (io,iostat=ierr) (read_int(i), i=1,np)
  !
  do i = 1,np
    int_data(i) = int( read_int(i) , kind=kind(int_data) )
  end do
  !
  deallocate ( read_int , stat=ierr )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_serial_int16_data_i8
!
!###############################################################################
!###############################################################################
! Serial read subroutine for reading 32-bit INTEGER data into a 64-bit array.
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_serial_int32_data_i8(io,int_data)
  !
  !.. Formal Arguments ..
  integer,        intent(in) :: io
  integer(i8), intent(inout) :: int_data(:)
  !
  !.. Local Scalars ..
  integer :: i,ierr,np
  !
  !.. Local Allocatable Arrays ..
  integer(int32), allocatable :: read_int(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_serial_int32_data_i8"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(int_data)
  !
  allocate ( read_int(1:np) , stat=ierr )
  !
  read (io,iostat=ierr) (read_int(i), i=1,np)
  !
  do i = 1,np
    int_data(i) = int( read_int(i) , kind=kind(int_data) )
  end do
  !
  deallocate ( read_int , stat=ierr )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_serial_int32_data_i8
!
!###############################################################################
!###############################################################################
! Serial read subroutine for reading 64-bit INTEGER data into a 64-bit array.
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_serial_int64_data_i8(io,int_data)
  !
  !.. Formal Arguments ..
  integer,        intent(in) :: io
  integer(i8), intent(inout) :: int_data(:)
  !
  !.. Local Scalars ..
  integer :: i,ierr,np
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_serial_int64_data_i8"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = size(int_data)
  !
  read (io,iostat=ierr) (int_data(i), i=1,np)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_serial_int64_data_i8
!
!###############################################################################
!###############################################################################
! Serial read subroutines for 32-bit REAL data
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_serial_r4_rank1(io,real_data)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: io
  real(wp), intent(inout) :: real_data(:)
  !
  !.. Local Scalars ..
  integer :: i,ierr
  !
  !.. Local Arrays ..
  real(real32) :: read_real(1:size(real_data))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_serial_r4_rank1"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  read (io,iostat=ierr) (read_real(i), i=1,size(read_real))
  !
  do i = 1,size(real_data,dim=1)
    real_data(i) = real( read_real(i) , kind=wp )
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_serial_r4_rank1
!-------------------------------------------------------------------------------
subroutine read_rstrt_file_serial_r4_rank2(io,real_data)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: io
  real(wp), intent(inout) :: real_data(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,ierr
  !
  !.. Local Arrays ..
  real(real32) :: read_real(1:size(real_data,dim=1), &
                            1:size(real_data,dim=2))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_serial_r4_rank2"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  read (io,iostat=ierr) ((read_real(i,n), i=1,size(read_real,dim=1)), &
                                          n=1,size(read_real,dim=2))
  !
  do n = 1,size(real_data,dim=2)
    do i = 1,size(real_data,dim=1)
      real_data(i,n) = real( read_real(i,n) , kind=wp )
    end do
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_serial_r4_rank2
!
!###############################################################################
!###############################################################################
! Serial read subroutines for 64-bit REAL data
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_serial_r8_rank1(io,real_data)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: io
  real(wp), intent(inout) :: real_data(:)
  !
  !.. Local Scalars ..
  integer :: i,ierr
  !
  !.. Local Arrays ..
  real(real64) :: read_real(1:size(real_data))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_serial_r8_rank1"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  read (io,iostat=ierr) (read_real(i), i=1,size(read_real))
  !
  do i = 1,size(real_data,dim=1)
    real_data(i) = real( read_real(i) , kind=wp )
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_serial_r8_rank1
!-------------------------------------------------------------------------------
subroutine read_rstrt_file_serial_r8_rank2(io,real_data)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: io
  real(wp), intent(inout) :: real_data(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,ierr
  !
  !.. Local Arrays ..
  real(real64) :: read_real(1:size(real_data,dim=1), &
                            1:size(real_data,dim=2))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_serial_r8_rank2"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  read (io,iostat=ierr) ((read_real(i,n), i=1,size(read_real,dim=1)), &
                                          n=1,size(read_real,dim=2))
  !
  do n = 1,size(real_data,dim=2)
    do i = 1,size(real_data,dim=1)
      real_data(i,n) = real( read_real(i,n) , kind=wp )
    end do
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_serial_r8_rank2
!
!###############################################################################
!###############################################################################
! Serial read subroutines for 128-bit REAL data
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_serial_r16_rank1(io,real_data)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: io
  real(wp), intent(inout) :: real_data(:)
  !
  !.. Local Scalars ..
  integer :: i,ierr
  !
  !.. Local Arrays ..
  real(max(real64,real128)) :: read_real(1:size(real_data))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_serial_r16_rank1"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  read (io,iostat=ierr) (read_real(i), i=1,size(read_real))
  !
  do i = 1,size(real_data,dim=1)
    real_data(i) = real( read_real(i) , kind=wp )
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_serial_r16_rank1
!-------------------------------------------------------------------------------
subroutine read_rstrt_file_serial_r16_rank2(io,real_data)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: io
  real(wp), intent(inout) :: real_data(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,ierr
  !
  !.. Local Arrays ..
  real(max(real64,real128)) :: read_real(1:size(real_data,dim=1), &
                                         1:size(real_data,dim=2))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_serial_r16_rank2"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  read (io,iostat=ierr) ((read_real(i,n), i=1,size(read_real,dim=1)), &
                                          n=1,size(read_real,dim=2))
  !
  do n = 1,size(real_data,dim=2)
    do i = 1,size(real_data,dim=1)
      real_data(i,n) = real( read_real(i,n) , kind=wp )
    end do
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_serial_r16_rank2
!
!###############################################################################
!###############################################################################
! Serial read subroutines for current working precision REAL data (wp)
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_serial_wp_rank1(io,real_data)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: io
  real(wp), intent(inout) :: real_data(:)
  !
  !.. Local Scalars ..
  integer :: i,ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_serial_wp_rank1"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  read (io,iostat=ierr) (real_data(i), i=1,size(real_data))
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_serial_wp_rank1
!-------------------------------------------------------------------------------
subroutine read_rstrt_file_serial_wp_rank2(io,real_data)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: io
  real(wp), intent(inout) :: real_data(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_serial_wp_rank2"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  read (io,iostat=ierr) ((real_data(i,n), i=1,size(real_data,dim=1)), &
                                          n=1,size(real_data,dim=2))
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_serial_wp_rank2
!
!###############################################################################
!
subroutine interp_restart_solution(restart_solution_form,rst_sol,usp, &
                                   ordered_geom)
  !
  !.. Use Statements ..
  use order_mod,         only : n_min_geom,n_min_order
  use order_mod,         only : n_max_geom,n_max_order
  use geomvar_mod,       only : ncell,cell
  use ovar,              only : loc_solution_pts
  use interpolation_mod, only : compute_generic_interp_matrix
  use projection_mod,    only : compute_generic_projct_matrix
  !
  use matrix_types_mod,  only : matrix
  use matrix_types_mod,  only : create_matrix_sparsity
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: restart_solution_form
  real(wp),    intent(in) :: rst_sol(:,:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Optional Arguments ..
  integer(ogp), optional, intent(in) :: ordered_geom(:)
  !
  !.. Local Derived Types ..
  type :: rstrt_interp_t
    integer :: omin = huge(0)
    integer :: omax = -huge(0)
    type(matrix), allocatable :: fromCell(:)
  end type rstrt_interp_t
  !
  !.. Local Scalars ..
  integer :: i,j,n,nc,np,n0,lp,l0,l1,l2
  integer :: gmin,gmax,omin,omax
  integer :: rst_omin,rst_omax,ierr
  integer :: this_geom,this_order
  integer :: rst_geom,rst_order
  integer :: loc_rst_pts
  !
  real(wp) :: a
  !
  character(len=200) :: array_name
  !
  !.. Local Arrays ..
  real(wp) :: val(1:nq)
  !
  !.. Local Allocatable Arrays ..
  type(rstrt_cell_t),   allocatable :: rst_cell(:)
  type(rstrt_interp_t), allocatable :: rst_interp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interp_restart_solution"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  loc_rst_pts = abs(restart_solution_form)
  !
  ! If ordered_geom is present, this is for reading a restart file
  ! in serial and cell_input_data doesnt exist, so create rst_cell
  ! using ordered_geom.
  !
  if (present(ordered_geom)) then
    !
    allocate ( rst_cell(1:ncell) )
    !
    l2 = 0
    do nc = 1,ncell
      !
      rst_order = int(ordered_geom(nc),kind=kind(rst_order)) / 100
      rst_geom = int(ordered_geom(nc),kind=kind(rst_geom)) - 100*rst_order
      !
      l1 = l2 + 1
      l2 = l2 + Cell_Solpts(rst_geom,rst_order)
      !
      rst_cell(nc)%input_geom = rst_geom
      rst_cell(nc)%input_order = rst_order
      rst_cell(nc)%input_beg_sp = l1
      rst_cell(nc)%input_end_sp = l2
      rst_cell(nc)%input_cell = nc
      !
    end do
    !
  else
    !
    ! Copy cell_input_data into rst_cell
    !
    rst_cell = cell_input_data(:) ! F2003 AUTO-REALLOC
    !
  end if
  !
  ! Get the min/max order of the restart solution
  !
  rst_omin = minval(rst_cell(:)%input_order)
  rst_omax = maxval(rst_cell(:)%input_order)
  !
  gmin = n_min_geom
  gmax = n_max_geom
  omin = n_min_order
  omax = n_max_order
  !
  ! Allocate the rst_interp array
  !
  allocate ( rst_interp(gmin:gmax,omin:omax) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"rst_interp",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Find the max/min restart order for each geom/order cell combination
  !
  do nc = 1,ncell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    rst_order = rst_cell(nc)%input_order
    !
    if (rst_order < rst_interp(this_geom,this_order)%omin) then
      rst_interp(this_geom,this_order)%omin = rst_order
    end if
    !
    if (rst_order > rst_interp(this_geom,this_order)%omax) then
      rst_interp(this_geom,this_order)%omax = rst_order
    end if
    !
  end do
  !
  ! Allocate the fromCell components of the rst_interp array
  ! and create all the interpolation/projection matrices
  !
  do this_geom = gmin,gmax
    do this_order = omin,omax
      !
      ! Extract the min/max restart cell order for the current geom/order combo
      !
      rst_omin = rst_interp(this_geom,this_order)%omin
      rst_omax = rst_interp(this_geom,this_order)%omax
      !
      ! If the minimum restart order is greater than the maximum restart order
      ! then there are no restart cells of this geom/order combo so cycle to
      ! the next combo.
      !
      if (rst_omin > rst_omax) cycle
      !
      ! Allocate the fromCell component of rst_interp for this geom/order combo
      !
      allocate ( rst_interp(this_geom,this_order)%fromCell(rst_omin:rst_omax), &
                 stat=ierr, errmsg=error_message )
      write (array_name,1) Geom_Name(this_geom),this_order,"fromCell"
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Create interpolation matrices for all restart orders
      ! less than or equal to this_order
      !
      do rst_order = rst_omin,min(this_order,rst_omax)
        !
        ! Compute the interpolation matrix for this geometry to 
        ! interpolate from a restart cell of order rst_order to a
        ! solution cell of order this_order.
        !
        ! F2003 AUTO-REALLOC
        rst_interp(this_geom,this_order)%fromCell(rst_order)%mat = &
                          compute_generic_interp_matrix(this_geom, &
                                                        from_order=rst_order, &
                                                        from_loc=loc_rst_pts, &
                                                        to_order=this_order, &
                                                        to_loc=loc_solution_pts)
        !
      end do
      !
      ! Create projection matrices for all restart orders
      ! greater than this_order
      !
      do rst_order = max(this_order+1,rst_omin),rst_omax
        !
        ! Compute the projection matrix for this geometry to 
        ! project from a restart cell of order rst_order to a
        ! solution cell of order this_order.
        !
        ! F2003 AUTO-REALLOC
        rst_interp(this_geom,this_order)%fromCell(rst_order)%mat = &
                          compute_generic_projct_matrix(this_geom, &
                                                        from_order=rst_order, &
                                                        from_loc=loc_rst_pts, &
                                                        to_order=this_order, &
                                                        to_loc=loc_solution_pts)
        !
      end do
      !
    end do
  end do
  !
  ! Create the sparse matrix arrays for rst_interp
  !
  do this_geom = gmin,gmax
    do this_order = omin,omax
      if (allocated(rst_interp(this_geom,this_order)%fromCell)) then
        call create_matrix_sparsity( rst_interp(this_geom,this_order)%fromCell )
      end if
    end do
  end do
  !
  ! Now interpolate/project the restart solution
  ! to the solution points of each cell
  !
  do nc = 1,ncell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    rst_order = rst_cell(nc)%input_order
    !
    l0 = rst_cell(nc)%input_beg_sp - 1
    lp = rst_cell(nc)%input_end_sp - l0
    !
    do n = 1,np
      !
      val(1:nq) = zero
      !
      if (rst_order > this_order) then
        !
        ! Use the row sparse arrays for projecting
        !
        do i = rst_interp(this_geom,this_order)% &
                            fromCell(rst_order)%row_ptr(n)+1, &
               rst_interp(this_geom,this_order)% &
                            fromCell(rst_order)%row_ptr(n+1)
          !
          j = rst_interp(this_geom,this_order)%fromCell(rst_order)%row_idx(i)
          a = rst_interp(this_geom,this_order)%fromCell(rst_order)%row_val(i)
          !
          val(1:nq) = val(1:nq) + a * rst_sol(1:nq,l0+j)
          !
        end do
        !
      else
        !
        ! Use the row sparse arrays for interpolating
        !
        do i = rst_interp(this_geom,this_order)% &
                            fromCell(rst_order)%row_ptr(n)+1, &
               rst_interp(this_geom,this_order)% &
                            fromCell(rst_order)%row_ptr(n+1)
          !
          j = rst_interp(this_geom,this_order)%fromCell(rst_order)%row_idx(i)
          a = rst_interp(this_geom,this_order)%fromCell(rst_order)%row_val(i)
          !
          val(1:nq) = val(1:nq) + a * rst_sol(1:nq,l0+j)
          !
        end do
        !
      end if
      !
      ! Finally, save the interpolated/projected restart solution into usp
      !
      usp(n0+n,1:nq) = val(1:nq)
      !
    end do
    !
  end do
  !
  ! Deallocate local allocatable arrays before leaving
  !
  if (allocated(rst_interp)) deallocate ( rst_interp )
  if (allocated(rst_cell)) deallocate ( rst_cell )
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("rst_interp(",a,",",i0,")%",a,:,"(",i0,")%mat")
  !
end subroutine interp_restart_solution
!
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine save_ta_data(ta_data,ta_filename,int_header,real_header,var_header, &
                        ordered_geom,solution)
  !
  !.. Formal Arguments ..
  type(ta_data_t), allocatable, intent(out) :: ta_data
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: ta_filename
  integer(wip),     optional, intent(in) :: int_header(:)
  real(wp),         optional, intent(in) :: real_header(:)
  integer(wip),     optional, intent(in) :: var_header(:)
  integer(ogp),     optional, intent(in) :: ordered_geom(:)
  real(wp),         optional, intent(in) :: solution(:,:)
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "save_ta_data"
  !
continue
  !
  allocate ( ta_data , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ta_data",1,__LINE__,__FILE__,ierr,error_message)
  !
  if (present(solution)) then
    ! F2003 AUTO-REALLOC
    ta_data%sol = real( solution(:,:) , kind=kind(ta_data%sol) )
  end if
  !
  if (present(ordered_geom)) then
    ta_data%ordered_geom = int( ordered_geom(:) , & ! F2003 AUTO-REALLOC
                                kind=kind(ta_data%ordered_geom))
  end if
  !
  if (present(int_header)) then
    ! F2003 AUTO-REALLOC
    ta_data%int_header = int(int_header(:) , kind=kind(ta_data%int_header))
  end if
  !
  if (present(real_header)) then
    ! F2003 AUTO-REALLOC
    ta_data%real_header = real(real_header(:) , kind=kind(ta_data%real_header))
  end if
  !
  if (present(var_header)) then
    ! F2003 AUTO-REALLOC
    ta_data%var_header = int(var_header(:) , kind=kind(ta_data%var_header))
  end if
  !
  if (present(ta_filename)) then
    ta_data%ta_filename = trim(adjustl(ta_filename)) ! F2003 AUTO-REALLOC
  end if
  !
end subroutine save_ta_data
!
!###############################################################################
!
subroutine modify_time_average_restart_history()
  !
  !.. Use Statements ..
  use ovar, only : restart_list_file
  use ovar, only : time_ref,time_scaling_factor
  use ovar, only : delete_converted_restart_files
  !
  !.. Local Scalars ..
  integer  :: ierr,io,i,n
  real(wp) :: out_time,time_scale
  !
  logical(ldk) :: file_exists
  !
  character(len=150) :: text
  character(len=300) :: fname
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: list_file_name
  character(len=:), allocatable :: read_file_name
  !
  !.. Local Allocatable Arrays ..
  character(len=300), allocatable :: list_of_restart_files(:)
  character(len=300), allocatable :: temp_list(:)
  !
  !.. Local Allocatable Derived Types ..
  type(ta_data_t), allocatable :: new_start
  type(ta_data_t), allocatable :: modified_ta_data
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "modify_time_average_restart_history"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  time_scale = time_ref*time_scaling_factor
  !
  ! First deallocate unnecessary arrays that were allocated for the solver but
  ! wont be needed since execution will abort at the end of this subroutine
  !
  call deallocate_solver_arrays
  !
  ! ###########
  ! ###########
  ! ###  1  ###
  ! ###########
  ! ###########
  !
  ! Check to make sure the restart_list_file exists
  !
  list_file_name = trim(adjustl(restart_list_file)) ! F2003 AUTO-REALLOC
  inquire (file=list_file_name,exist=file_exists)
  !
  if (.not. file_exists .or. len_trim(adjustl(restart_list_file)) == 0) then
    write (error_message,100)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! ###########
  ! ###########
  ! ###  2  ###
  ! ###########
  ! ###########
  !
  ! Count the number of lines in the file restart_list_file
  !
  open (newunit=io,file=list_file_name,status="old",action="read", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,list_file_name,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Count the number of restart files given in the list of files
  !
  n = 0
  count_file_lines_loop: &
  do
    read (io,1,iostat=ierr,iomsg=error_message) text
    if (ierr < 0) then
      ! There was an end-of-file (or possibly an end-of record) error,
      ! so assume we have hit the end of the file and exit this loop.
      exit count_file_lines_loop
    else if (ierr > 0) then
      ! There was some other kind of error reading from the file 
      call io_error(pname,list_file_name,-1,__LINE__,__FILE__, &
                    ierr,error_message)
    else
      n = n + 1
    end if
  end do count_file_lines_loop
  !
  rewind (io)
  !
  ! ###########
  ! ###########
  ! ###  3  ###
  ! ###########
  ! ###########
  !
  ! Read in the list of time-average restart files that will be used
  ! to create new time-average restart files with modified histories.
  !
  ! Allocate the arrays for storing the file names from the list of files
  !
  allocate ( list_of_restart_files(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"list_of_restart_files",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( temp_list(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_list",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  ! Read the file names from the list of files
  !
  do n = 1,size(list_of_restart_files)
    read (io,1,iostat=ierr,iomsg=error_message) list_of_restart_files(n)
    call io_error(pname,list_file_name,-1,__LINE__,__FILE__,ierr,error_message)
  end do
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,list_file_name,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! ###########
  ! ###########
  ! ###  4  ###
  ! ###########
  ! ###########
  !
  ! Check the validity of the filenames in the restart file list.
  !   1) Ignore any filenames that start with !
  !   2) For those filenames not starting with !, use inquire to check whether
  !      or not the file exists. If the file does exist, add the filename
  !      to temp_list.
  !
  i = 0
  do n = 1,size(list_of_restart_files)
    ! F2003 AUTO-REALLOC
    read_file_name = trim(adjustl(list_of_restart_files(n)))
    if (mypnum == glb_root) write (iout,10) n,read_file_name
    if (read_file_name(1:1) /= "!") then
      if (mypnum == glb_root) write (iout,11) read_file_name
      inquire (file=read_file_name,exist=file_exists)
      if (file_exists) then
        i = i + 1
        if (mypnum == glb_root) write (iout,12) i
        temp_list(i) = read_file_name
      else
        if (mypnum == glb_root) write (iout,13)
      end if
    else
      if (mypnum == glb_root) write (iout,14)
    end if
  end do
  flush (iout)
  !
  ! Use F2003 AUTO-REALLOCATION to cutoff the unused parts of temp_list
  !
  temp_list = temp_list(1:i) ! F2003 AUTO-REALLOC
  !
  ! Use move_alloc to over-write list_of_restart_files with temp_list
  !
  call move_alloc( from=temp_list , to=list_of_restart_files )
  !
  ! ###########
  ! ###########
  ! ###  5  ###
  ! ###########
  ! ###########
  !
  ! Loop over the list of files to read each restart
  ! file and write out a CGNS solution file
  !
  loop_over_restart_files: &
  do n = 1,size(list_of_restart_files)
    !
    ! F2003 AUTO-REALLOC
    read_file_name = trim(adjustl(list_of_restart_files(n)))
    !
    if (mypnum == glb_root) then
      write (iout,2) trim(adjustl(read_file_name)), &
                     n, size(list_of_restart_files)
    end if
    !
    if (allocated(new_start)) then
      deallocate ( new_start , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"new_start",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    allocate ( new_start , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"new_start",1,__LINE__,__FILE__,ierr,error_message)
    !
    new_start%ta_filename = read_file_name ! F2003 AUTO-REALLOC
    !
    call read_solution_file(new_start%ta_filename, &
                            new_start%int_header, &
                            new_start%real_header, &
                            new_start%sol, &
                            new_start%ordered_geom, &
                            var_header=new_start%var_header, &
                            ignore_sem=true)
    !
    ! Modify the time-average data from the original restart file using
    ! this new time-average as the new starting point of the averaging.
    !
    call remove_time_average_history(original_ta_data,new_start, &
                                     modified_ta_data,ierr)
    !
    if (ierr /= 0) then
      cycle loop_over_restart_files
    end if
    !
    if (mypnum == glb_root) then
      out_time = (modified_ta_data%real_header(1) - &
                  modified_ta_data%real_header(3) )*time_scale
      write (iout,3) (new_start%real_header(1) - &
                      new_start%real_header(3))*time_scale, &
                     (original_ta_data%real_header(1) - &
                      original_ta_data%real_header(3))*time_scale, &
                     out_time
      flush (iout)
    end if
   !!
   !! Update the parts of the real header that need to altered with the
   !! change in time-average history
   !! NOTE: In actuality, real_header(4)=rl2mx could possibly need to be
   !!       updated as well, but it is impossible to know for sure or even
   !!       what the correct rl2mx would be since it is just the maximum
   !!       l2 norm of the density residual over all the iterations since
   !!       the time averaging began. We dont know what iteration rl2mx
   !!       occurred, and we dont have a history of rl2mx to know what
   !!       it should be from the new time-averaging interval. We will just
   !!       leave it as is since it is really just used for the residual.dat
   !!       file which is not all that useful for LES.
   !!
   !real_header(3) = new_ave_start_time
    !
    ! Write the new time-average restart file
    ! NOTE: We will just reuse the contents of the ordered_geom, int_header,
    !       and var_header arrays since those remain the same as they were
    !       before modifying the time-average history.
    !
    call write_modified_time_averaged_solution(modified_ta_data%ta_filename, &
                                               modified_ta_data%int_header, &
                                               modified_ta_data%real_header, &
                                               modified_ta_data%var_header, &
                                               modified_ta_data%ordered_geom, &
                                               modified_ta_data%sol)
    !
    ! If requested, delete the time-average restart file that was just used to
    ! modify the original time-average restart file.
    !
    if (delete_converted_restart_files) then
      call delete_file(read_file_name,pname)
    end if
    !
  end do loop_over_restart_files
  !
  call debug_timer(leaving_procedure,pname)
  !
  write (error_message,4)
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  !
  ! Format Statements
  !
  1 format (a)
  2 format (/,3(1x,100("#"),/),/, &
              3x,"Working on time-average restart file '",a,"' ", &
                 "(",i3," of ",i3,")",//, &
              3(1x,100("#"),/),/, &
              4x,"Reading the time-average restart file")
  3 format (/,1x,71("-"),/, &
              " -----  Creating new time-average history by ", &
                                             "removing the initial  -----",/, &
              " -----  ",f10.4," tstar of time-averaging read ", &
                                                "in from this       -----",/, &
              " -----  time-average restart file. The new ", &
                                           "time-average will be    -----",/, &
              " -----  reduced from ",f10.4," to ",f10.4," ", &
                                              "tstar of averaging.  -----",/, &
              1x,71("-"),/)
  4 format ("STOPPING EXECUTION AFTER MODIFYING THE TIME-AVERAGE RESTART ", &
            "HISTORY USING THE TIME-AVERAGE DATA FROM THE LIST OF ", &
            "TIME-AVERAGE RESTART FILES!")
  !
  10 format (/,2x,"Filename # ",i0," from restart file list = '",a,"'")
  11 format (4x,"Checking existance of time-average restart file '",a,"'")
  12 format (4x,"Time-average restart file found! Saving to temp_list(",i0,")")
  13 format (8x,"ERROR: TIME-AVERAGE RESTART FILE DOES NOT EXIST!",/, &
             8x,"       IGNORING THIS LINE FROM THE RESTART FILE LIST!")
  14 format (8x,"WARNING: IGNORING FILENAME BEGINNING WITH '!'...")
  !
  20 format ("out.",i0,"-",i0,".modified.ave")
  !
  100 format ("THE FILE SPECIFIED BY THE INPUT - RESTART_LIST_FILE - COULD ", &
              "NOT BE FOUND! THIS FILE MUST PROVIDE A LIST OF ONE OR MORE ", &
              "TIME-AVERAGE RESTART FILES FROM THE SAME LINEAGE AS THE ", &
              "CURRENT TIME-AVERAGE RESTART FILE THAT WAS SPECIFIED BY THE ", &
              "INPUT - TIME_AVE_FILE - IN ORDER TO MODIFY THE HISTORY OF ", &
              "THE CURRENT TIME-AVERAGED DATA.")
  !
end subroutine modify_time_average_restart_history
!
!###############################################################################
!
subroutine remove_time_average_history(original,remove,modified,return_error)
  !
  !.. Use Statements ..
  use flowvar_mod, only : time_ave_variables
  !
  !.. Formal Arguments ..
  type(ta_data_t),              intent(in   ) :: original
  type(ta_data_t),              intent(in   ) :: remove
  type(ta_data_t), allocatable, intent(  out) :: modified
  integer,                      intent(inout) :: return_error
  !
  !.. Local Scalars ..
  integer  :: i,j,iii,ierr
  real(wp) :: orig_final_time,orig_start_time,orig_tstar,orig_coef
  real(wp) :: drop_final_time,drop_start_time,drop_tstar,drop_coef
  !
  character(len=1000) :: fname
  !
  character(len=len(time_ave_variables)) :: o_ta_var,s_ta_var
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "remove_time_average_history"
  !
 !logical(lk), parameter :: simple_F2003_allocation = fals
  logical(lk), parameter :: simple_F2003_allocation = true
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  return_error = 1
  !
  single_pass_loop: &
  do iii = 1,1
    !
    !===========================================================================
    !
    ! #####
    ! # 1 #
    ! #####
    !
    ! First check that the original restart data is compatible with the remove
    ! restart data that will be used as the new starting point of the modified
    ! time-average data. If they arent compatible, there is no point in doing
    ! any of the rest of the stuff in this subroutine.
    !
    return_error = 1
    !
    ! Check that the number of time-averaged variables is the same
    !
    if (original%int_header(1) /= remove%int_header(1)) then
      if (mypnum == glb_root) then
        write (iout,100) "# of time-average variables"
        write (iout,101) original%int_header(1), remove%int_header(1)
      end if
      exit single_pass_loop
    end if
    !
    ! Check that the number of grid cells is the same
    !
    if (original%int_header(2) /= remove%int_header(2)) then
      if (mypnum == glb_root) then
        write (iout,100) "# of grid cells"
        write (iout,101) original%int_header(2), remove%int_header(2)
      end if
      exit single_pass_loop
    end if
    !
    ! Check that time averages were started at the same solution time. The
    ! only way the math behind this modification works is if both sets of
    ! time-average data were started at the same solution time.
    !
    if (abs(original%real_header(3)-remove%real_header(3)) > eps6) then
      if (mypnum == glb_root) then
        write (iout,100) "remove time of accumulating time-average data"
        write (iout,101) original%real_header(3), remove%real_header(3)
      end if
      exit single_pass_loop
    end if
    !
    ! Check that all the time-average variables are the same
    !
    ierr = 0
    do i = 1,size(original%var_header)
      if (original%var_header(i) /= remove%var_header(i)) then
        ierr = ierr + 1
        if (mypnum == glb_root) then
          if (ierr == 1) then
            write (iout,100) "actual variables being time averaged"
          end if
          o_ta_var = transfer(original%var_header(i),o_ta_var)
          s_ta_var = transfer(remove%var_header(i),s_ta_var)
          write (iout,102) ierr,trim(adjustl(o_ta_var)),trim(adjustl(s_ta_var))
        end if
      end if
    end do
    if (ierr > 0) then
      if (mypnum == glb_root) write (iout,103) ierr
      exit single_pass_loop
    end if
    !
    ! Check that ordered_geom matches for all cells
    !
    ierr = 0
    do i = 1,size(original%ordered_geom)
      if (original%ordered_geom(i) /= remove%ordered_geom(i)) then
        ierr = ierr + 1
      end if
    end do
    if (ierr > 0) then
      if (mypnum == glb_root) then
        write (iout,100) "contents of ORDERED_GEOM array"
        write (iout,103) ierr
      end if
      exit single_pass_loop
    end if
    !
    !===========================================================================
    !
    ! #####
    ! # 2 #
    ! #####
    !
    return_error = 2
    !
    ! Allocate the input argument modified as a copy of the argument original
    !
    if (simple_F2003_allocation) then
      !
      ! Copy the contents of the original data into modified.
      ! NOTE: This should automatically allocate modified, and for each
      !       component of original, automatically allocate the respective
      !       component of modified to the same size, and then copy the data
      !       from each compenent of original to the respective component of
      !       modified.
      !
      modified = original ! F2003 AUTO-REALLOC
      !
    else
      !
      ! Allocate modified
      ! NOTE: We dont need to check if modified is allocated already since
      !       it has the intent(out) specification, which means Fortran will
      !       automatically deallocate it upon entry into this subroutine.
      !
      allocate ( modified , stat=ierr , errmsg=cerr )
      if (ierr /= 0) then
        write (iout,200) "modified",ierr,trim(adjustl(cerr))
        exit single_pass_loop
      end if
      !
      ! Allocate the sol component of modified using the same component
      ! of original as the source to copy the data over to modified.
      !
      allocate ( modified%sol , source=original%sol , stat=ierr , errmsg=cerr )
      if (ierr /= 0) then
        write (iout,200) "modified%sol",ierr,trim(adjustl(cerr))
        exit single_pass_loop
      end if
      !
      ! Allocate the ordered_geom component of modified using the same component
      ! of original as the source to copy the data over to modified.
      !
      allocate ( modified%ordered_geom , source=original%ordered_geom , &
                 stat=ierr , errmsg=cerr )
      if (ierr /= 0) then
        write (iout,200) "modified%ordered_geom",ierr,trim(adjustl(cerr))
        exit single_pass_loop
      end if
      !
      ! Allocate the int_header component of modified using the same component
      ! of original as the source to copy the data over to modified.
      !
      allocate ( modified%int_header , source=original%int_header , &
                 stat=ierr , errmsg=cerr )
      if (ierr /= 0) then
        write (iout,200) "modified%int_header",ierr,trim(adjustl(cerr))
        exit single_pass_loop
      end if
      !
      ! Allocate the real_header component of modified using the same component
      ! of original as the source to copy the data over to modified.
      !
      allocate ( modified%real_header , source=original%real_header , &
                 stat=ierr , errmsg=cerr )
      if (ierr /= 0) then
        write (iout,200) "modified%real_header",ierr,trim(adjustl(cerr))
        exit single_pass_loop
      end if
      !
      ! Allocate the var_header component of modified using the same component
      ! of original as the source to copy the data over to modified.
      !
      allocate ( modified%var_header , source=original%var_header , &
                 stat=ierr , errmsg=cerr )
      if (ierr /= 0) then
        write (iout,200) "modified%var_header",ierr,trim(adjustl(cerr))
        exit single_pass_loop
      end if
      !
    end if
    !
    !===========================================================================
    !
    ! #####
    ! # 3 #
    ! #####
    !
    ! Remove the averaging coefficient (i.e, denominator term) from the
    ! original time-averaged data that was copied into modified and the data
    ! in remove that will be subtracted from the original data to shift the
    ! starting solution time of the modified time-average data.
    !
    orig_final_time = modified%real_header(1)
    orig_start_time = modified%real_header(3)
    !
    drop_final_time = remove%real_header(1)
    drop_start_time = remove%real_header(3) ! should be same as orig_start_time
    !
    orig_tstar = orig_final_time - orig_start_time
    drop_tstar = drop_final_time - drop_start_time
    !
    orig_coef = orig_tstar / (orig_tstar - drop_tstar)
    drop_coef = one - orig_coef ! = -drop_tstar / (orig_tstar - drop_tstar)
    !
    do j = 1,size(modified%sol,dim=2)
      do i = 1,size(modified%sol,dim=1)
        modified%sol(i,j) = orig_coef*modified%sol(i,j) + &
                            drop_coef*remove%sol(i,j)
      end do
    end do
    !
    !===========================================================================
    !
    ! #####
    ! # 4 #
    ! #####
    !
    ! Update the starting solution time for the modified time-average data
    ! in the real_header array
    !
    modified%real_header(3) = drop_final_time
    !
    !===========================================================================
    !
    ! #####
    ! # 5 #
    ! #####
    !
    ! Create the file name for the new modified time-average restart file
    !
    write (fname,300) original%int_header(4), remove%int_header(4)
    !
    modified%ta_filename = trim(adjustl(fname)) ! F2003 AUTO-REALLOC
    !
    !===========================================================================
    !
    ! If we have reached this point, there doesnt seem to be any errors with
    ! modifying the time average history, so set the return error to 0.
    !
    return_error = 0
    !
  end do single_pass_loop
  !
  flush (iout)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  100 format (/,1x,80("*"),/, &
                " ERROR: Incompatibility between original and remove data!",/, &
                10x,"Mismatch between  '",a,"'")
  101 format (14x,"Original = ",g0,/, &
              14x,"  Remove = ",g0,/, &
              1x,80("*"),/)
  102 format (14x,"Variable mismatch # ",i0,/, &
              18x,"Original = '",a,"'",/, &
              18x,"  Remove = '",a,"'")
  103 format (14x,"# of mismatches = ",i0, &
              1x,80("*"),/)
  !
  200 format (/,1x,80("*"),/, &
                " ERROR ALLOCATING : ",a,/, &
                "      ERROR INDICATOR = ",i0,/, &
                "      ERROR MESSAGE   = '",a,"'",/, &
                1x,80("*"))
  !
  300 format ("out.",i0,"-",i0,".modified.ave")
  !
end subroutine remove_time_average_history
!
!###############################################################################
!
subroutine delete_file(file_name,procedure_name)
  !
  !###################################################################
  !###################################################################
  ! NOTE: USING THIS SUBROUTINE IS PERMANENT, AND THE FILE CANNOT BE
  !       RECOVERED IF THIS IS USED, SO MAKE SURE THIS IS ONLY CALLED
  !       WHEN YOU ARE ABSOLUTELY SURE YOU WANT TO DELETE THE FILE.
  !###################################################################
  !###################################################################
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: file_name
  character(len=*), intent(in) :: procedure_name
  !
  !.. Local Scalars ..
  integer :: io,ierr
  logical :: file_exists
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: fname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "delete_file"
  !
continue
  !
  fname = trim(adjustl(file_name)) ! F2003 AUTO-REALLOC
  !
  if (mypnum == glb_root) then
    !
    ! First check to see if the file exists
    !
    inquire (file=fname,exist=file_exists)
    !
    ! If it exists, procede to delete the file
    !
    if (file_exists) then
      !
      ! Open the file as read-write
      !
      open (newunit=io,file=fname,status="old",action="readwrite", &
                       access="stream",form="unformatted", &
                       iostat=ierr,iomsg=error_message)
      call io_error(procedure_name,fname,1,__LINE__,__FILE__, &
                    ierr,error_message)
      !
      ! Close the file, using the delete status to delete it
      !
      !###################################################################
      !###################################################################
      ! NOTE: THIS IS PERMANENT, AND THE FILE CANNOT BE RECOVERED
      !       IF THIS IS USED, SO MAKE SURE THIS IS ONLY CALLED WHEN
      !       YOU ARE ABSOLUTELY SURE YOU WANT TO DELETE THE FILE.
      !###################################################################
      !###################################################################
      !
      close (io,status="delete",iostat=ierr,iomsg=error_message)
      call io_error(procedure_name,fname,2,__LINE__,__FILE__, &
                    ierr,error_message)
      !
    end if
    !
  end if
  !
end subroutine delete_file
!
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine deallocate_solver_arrays()
  !
  !.. Use Statements ..
  use geomvar_mod, only : face,grid
  use flowvar_mod, only : faceusp_dt,facedusp_dt,faceflx,facexyz_dt
  use flowvar_mod, only : residual,dtsp,uoldsp,rk_sum,rssprk
  use flowvar_mod, only : nodeusp_dt,nodedusp_dt
  use flowvar_mod, only : edgeusp_dt,edgedusp_dt,edgexyz_dt
  use metrics_mod, only : metrics_dt,geofa,facenrm
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "deallocate_solver_arrays"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (allocated(face)) then
    deallocate ( face , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"face",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(grid)) then
    deallocate ( grid , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(faceusp_dt)) then
    deallocate ( faceusp_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"faceusp_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(facedusp_dt)) then
    deallocate ( facedusp_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facedusp_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(faceflx)) then
    deallocate ( faceflx , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"faceflx",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(facexyz_dt)) then
    deallocate ( facexyz_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facexyz_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(residual)) then
    deallocate ( residual , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"residual",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(dtsp)) then
    deallocate ( dtsp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"dtsp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(uoldsp)) then
    deallocate ( uoldsp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"uoldsp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(rk_sum)) then
    deallocate ( rk_sum , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"rk_sum",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(rssprk)) then
    deallocate ( rssprk , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"rssprk",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(nodeusp_dt)) then
    deallocate ( nodeusp_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodeusp_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(nodedusp_dt)) then
    deallocate ( nodedusp_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodedusp_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edgeusp_dt)) then
    deallocate ( edgeusp_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edgeusp_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edgedusp_dt)) then
    deallocate ( edgedusp_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edgedusp_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edgexyz_dt)) then
    deallocate ( edgexyz_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edgexyz_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(metrics_dt)) then
    deallocate ( metrics_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"metrics_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(geofa)) then
    deallocate ( geofa , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"geofa",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(facenrm)) then
    deallocate ( facenrm , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facenrm",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine deallocate_solver_arrays
!
!###############################################################################
!
subroutine collect_mem(string,mem,total_mem)
  !
  !.. Formal Arguments ..
  character(len=*),  intent(in) :: string
  integer(int64), intent(inout) :: mem
  integer(int64), intent(inout) :: total_mem
  !
  !.. Local Scalars ..
  real(wp) :: rmem
  character(len=1) :: cmem
  !
continue
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,mem,1_int_mpi,MPI_INTEGER8, &
                       MPI_SUM,MPI_COMM_WORLD,mpierr)
  end if
  !
  call get_mem_scale(mem,rmem,cmem)
  !
  total_mem = total_mem + mem
  !
  if (mypnum == glb_root) then
    write (iout,1) string,rmem,cmem
  end if
  flush (iout)
  !
  ! Reset mem for next set of arrays to deallocate
  !
  mem = 0_int64
  !
  ! Format Statements
  !
  1 format (2x,"Deallocating: ",a29," => Freed up ",f10.4, &
               " ",a1,"Bytes of memory!")
  !
end subroutine collect_mem
!
!###############################################################################
!
subroutine get_mem_scale(mem,rmem,cmem)
  !
  !.. Formal Arguments ..
  integer(int64),      intent(in) :: mem
  real(wp),         intent(inout) :: rmem
  character(len=1), intent(inout) :: cmem
  !
  !.. Local Scalars ..
  integer :: n,mem_exp
  integer(int64) :: lmem
  !
continue
  !
  lmem = mem/8_int64
  !
  mem_exp = 0
  !
  mem_scale_loop: &
  do n = 0,4
    lmem = lmem / 1024_int64
    if (lmem == 0_int64) then
      mem_exp = n
      exit mem_scale_loop
    end if
  end do mem_scale_loop
  !
  if (mem_exp == 0) then
    cmem = ""
  else if (mem_exp == 1) then
    cmem = "K"
  else if (mem_exp == 2) then
    cmem = "M"
  else if (mem_exp == 3) then
    cmem = "G"
  else if (mem_exp == 4) then
    cmem = "T"
  end if
  !
  rmem = real(mem/8_int64,kind=wp) / (1024.0_wp**mem_exp)
  !
end subroutine get_mem_scale
!
!###############################################################################
!###############################################################################
!###############################################################################
!
include "Functions/cell_solpts.f90"
#include "Functions/sort.f90"
!
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module restart_mod
