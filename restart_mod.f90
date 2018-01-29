module restart_mod
  !
#include <mpi_defs.h>
!!!#define VERBOSE_RESTART
  !
  use module_kind_types
  use eqn_idx, only : nq
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
  !
  !###########################
  !###  Module Parameters  ###
  !###########################
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
  !
  integer, save :: averaging_file_switch = 1
  !
  _MPI_DATA_TYPE_, save :: read_type
  _MPI_DATA_TYPE_, save :: read_block_type
  _MPI_DATA_TYPE_, save :: real_type
  !
  integer, save :: num_pts = 0
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
  !############################################################################
  !###  Generic interfaces for subroutines to read REAL data using MPI I/O  ###
  !############################################################################
  !
  interface read_rstrt_file_r4
    module procedure read_rstrt_file_r4_rank1, read_rstrt_file_r4_rank2
  end interface read_rstrt_file_r4
  !
  interface read_rstrt_file_r8
    module procedure read_rstrt_file_r8_rank1, read_rstrt_file_r8_rank2
  end interface read_rstrt_file_r8
  !
  interface read_rstrt_file_r16
    module procedure read_rstrt_file_r16_rank1, read_rstrt_file_r16_rank2
  end interface read_rstrt_file_r16
  !
  interface read_rstrt_file_wp
    module procedure read_rstrt_file_wp_rank1, read_rstrt_file_wp_rank2
  end interface read_rstrt_file_wp
  !
contains
!
!###############################################################################
!
subroutine read_restart_file
  !
  !.. Use Statements ..
  use ovar,    only : itrst,tmrst,dtrst,rl2mx
  use ovar,    only : prev_ave_time,ave_start_time
  use ovar,    only : restart_file,time_ave_file
  use ovar,    only : output_time_averaging
  use ovar,    only : read_time_ave_restart
  use ovar,    only : use_old_restart_format
  use flowvar, only : usp,uavesp,time_ave_variables
  !
  !.. Local Scalars ..
  integer :: n,ierr
  character(len=150) :: fname
  !
  !.. Local Allocatable Arrays ..
  integer(wip), allocatable :: int_header(:)
  integer(wip), allocatable :: var_header(:)
  real(wp),     allocatable :: real_header(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_restart_file"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  fname = trim(adjustl(restart_file))
  !
  if (use_old_restart_format) then
    !
    call read_restart_old_format(fname,usp,sol_is_transposed=fals)
    !
  else
    !
    ! Read the restart solution file
    ! NOTE: usp SHOULD BE ALLOCATED AT THIS POINT
    !
    call read_solution_file(fname,int_header,real_header,usp, &
                            sol_is_transposed=fals)
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
    if (output_time_averaging) then
      !
      if (read_time_ave_restart) then
        !
        fname = trim(adjustl(time_ave_file))
        !
        ! NOTE: uavesp SHOULD NOT BE ALLOCATED AT THIS POINT
        !
        call read_solution_file(fname,int_header,real_header,uavesp, &
                                var_header,sol_is_transposed=fals)
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
  end if
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
                              var_header,sol_is_transposed)
  !
  !.. Use Statements ..
  use geovar, only : ncell,cell,n_solpts
  use geovar, only : n_global_cell
  use ovar,   only : time_average_restart_files
  !
  !.. Formal Arguments ..
  character(len=*),              intent(in   ) :: file_name
  integer(wip),     allocatable, intent(  out) :: int_header(:)
  real(wp),         allocatable, intent(  out) :: real_header(:)
  real(wp),         allocatable, intent(inout) :: solution(:,:)
  !
  !.. Optional Arguments ..
  integer(wip), optional, allocatable, intent(out) :: var_header(:)
  logical(lk),  optional,               intent(in) :: sol_is_transposed
  !
  !.. Local Scalars ..
  integer :: l0,l1,l2,m,n,np,n0,n1,n2,ierr
  integer :: sol_nq
  integer :: restart_ip,restart_rp
  logical(lk) :: flip_sol_indices
  !
  _MPI_FILE_TYPE_ :: file_handle
  _MPI_INFO_TYPE_ :: file_info
  !.. Need to be MPI_OFFSET_KIND integers
  integer(MPI_OFFSET_KIND) :: file_disp
  integer(MPI_OFFSET_KIND) :: file_offset
  !
  integer :: solution_errors
  integer :: restart_nq,restart_ncells
  integer :: n_header_int,n_header_real
  !
  integer :: inp_cell,inp_geom,inp_order
  integer :: this_cell,this_geom,this_order
  !
  !.. Local Allocatable Arrays ..
  integer(wip),       allocatable :: ordered_geom(:)
  real(wp),           allocatable :: rst_sol(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_solution_file"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  flip_sol_indices = fals
  if (present(sol_is_transposed)) then
    flip_sol_indices = sol_is_transposed
  end if
  !
  !
  solution_errors = 0
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
  call mpi_info_create(file_info,mpierr)
  !
#endif
#endif
  !
  ! Open the restart file for MPI I/O
  !
  call mpi_file_open(MPI_COMM_WORLD,trim(adjustl(file_name)), &
                     MPI_MODE_RDONLY,file_info,file_handle,mpierr)
  !
  ! Have all processors read the header section of the restart file
  !
! ######################################
! ### Integer part of restart header ###
! ######################################
  !
  ! Read the integer type for the integer part of the restart header
  !
  call mpi_file_read_all(file_handle,restart_ip,1_int_mpi, &
                         MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
  !
  ! Read the number of integer header items
  !
  call mpi_file_read_all(file_handle,n_header_int,1_int_mpi, &
                         MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
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
  if (restart_ip == wip) then
    call read_rstrt_file_wip(file_handle,int_header)
  else
    select case (restart_ip)
      case (int32) ; call read_rstrt_file_i4(file_handle,int_header)
      case (int64) ; call read_rstrt_file_i8(file_handle,int_header)
      case (int8)  ; call read_rstrt_file_i1(file_handle,int_header)
      case (int16) ; call read_rstrt_file_i2(file_handle,int_header)
      case default
        write (error_message,20)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end select
  end if
  !
  ! Extract the information from the integer section of the header
  !
  ! Extract the number of flow variables saved at each solution point
  ! within the restart file
  !
  restart_nq = int( int_header(1) , kind=kind(restart_nq) )
  !
  if (allocated(solution)) then
    if (flip_sol_indices) then
      sol_nq = size(solution,dim=2)
    else
      sol_nq = size(solution,dim=1)
    end if
    if (restart_nq /= sol_nq) then
      write (error_message,10) restart_nq,sol_nq
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
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
  !
  call check_restart_solution_form(restart_solution_form)
  !
! ###################################
! ### Real part of restart header ###
! ###################################
  !
  ! Read the real type for the real part of the restart header
  !
  call mpi_file_read_all(file_handle,restart_rp,1_int_mpi, &
                         MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
  !
  ! Read the number of real header items
  !
  call mpi_file_read_all(file_handle,n_header_real,1_int_mpi, &
                         MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
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
    call read_rstrt_file_wp(file_handle,real_header)
  else
    select case (restart_rp)
      case (real64)  ; call read_rstrt_file_r8(file_handle,real_header)
      case (real32)  ; call read_rstrt_file_r4(file_handle,real_header)
      case (real128) ; call read_rstrt_file_r16(file_handle,real_header)
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
  call mpi_file_read_all(file_handle,restart_ip,1_int_mpi, &
                         MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
  !
  allocate ( ordered_geom(1:restart_ncells) , source=0_wip , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ordered_geom",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! All processors have to read the ordered_geom data so that each
  ! processor can compute the data offsets for each of the cells
  ! within its partition.
  !
  if (restart_ip == wip) then
    call read_rstrt_file_wip(file_handle,ordered_geom)
  else
    select case (restart_ip)
      case (int32) ; call read_rstrt_file_i4(file_handle,ordered_geom)
      case (int64) ; call read_rstrt_file_i8(file_handle,ordered_geom)
      case (int8)  ; call read_rstrt_file_i1(file_handle,ordered_geom)
      case (int16) ; call read_rstrt_file_i2(file_handle,ordered_geom)
      case default
        write (error_message,22)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end select
  end if
  !
! ################################################
! ### Variable header part of the restart file ###
! ################################################
  !
  if (present(var_header)) then
    !
    ! Read the real type for the solution data
    !
    call mpi_file_read_all(file_handle,restart_ip,1_int_mpi, &
                           MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
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
    if (restart_ip == wip) then
      call read_rstrt_file_wip(file_handle,var_header)
    else
      select case (restart_ip)
        case (int32) ; call read_rstrt_file_i4(file_handle,var_header)
        case (int64) ; call read_rstrt_file_i8(file_handle,var_header)
        case (int8)  ; call read_rstrt_file_i1(file_handle,var_header)
        case (int16) ; call read_rstrt_file_i2(file_handle,var_header)
        case default
          write (error_message,24)
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end select
    end if
    !
    if (.not. allocated(solution)) then
      if (flip_sol_indices) then
        allocate ( solution(1:n_solpts,1:restart_nq) , source=zero , &
                   stat=ierr , errmsg=error_message )
      else
        allocate ( solution(1:restart_nq,1:n_solpts) , source=zero , &
                   stat=ierr , errmsg=error_message )
      end if
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
  call mpi_file_read_all(file_handle,restart_rp,1_int_mpi, &
                         MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
  !
  ! Get the updated file displacement
  !
  call mpi_file_get_position(file_handle,file_offset,mpierr)
  call mpi_file_get_byte_offset(file_handle,file_offset,file_disp,mpierr)
  !
  ! Get the MPI derived types required to read the solution data section
  ! of the restart file if cell_input_data is not allocated
  !
  if (.not. allocated(cell_input_data)) then
    !
    call get_input_mpi_datatypes(restart_nq,restart_rp,real_type, &
                                 read_block_type,read_type, &
                                 num_pts,ordered_geom,cell_input_data)
    !
  end if
  !
  ! Deallocate ordered_geom since it is no longer needed
  !
  if (allocated(ordered_geom)) then
    deallocate ( ordered_geom , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ordered_geom",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  ! IF WE ARE TO MOVE THE POINT IN THE PRE-PROCESSING PHASE AT WHICH THE
  ! RESTART FILE IS READ, WE NEED TO MAKE SURE THAT N_SOLPTS HAS BEEN
  ! RECOMPUTED SO THAT IT IS LOCAL TO EACH PROCESSOR BECAUSE THIS IS NOT
  ! DONE WITHIN PARALLEL_MOD AFTER PARTITIONING THE GRID.
  !
  if (num_pts /= n_solpts) then
    write (iout,12) mypnum,num_pts,n_solpts
    solution_errors = 100
  end if
  !
  ! Sum up the value of solution_errors across all processors
  !
  call mpi_allreduce(MPI_IN_PLACE,solution_errors,1_int_mpi, &
                     mpi_inttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
  !
  if (solution_errors /= 0) then
    write (error_message,13)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Set the file view to write the ordered_geom array to the restart file
  !
  call mpi_file_set_view(file_handle,file_disp,read_block_type, &
                         read_type, trim(adjustl(file_datarep)), &
                         file_info, mpierr)
  !
  ! Allocate the rst_sol array to read in the solution restart data
  !
  allocate ( rst_sol(1:restart_nq,1:n_solpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"rst_sol",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read in the solution restart data
  !
  if (restart_rp == wp) then
    call read_rstrt_file_wp(file_handle,rst_sol, &
                            read_block_type)
  else
    select case (restart_rp)
      case (real64)  ; call read_rstrt_file_r8(file_handle,rst_sol, &
                                               read_block_type)
      case (real32)  ; call read_rstrt_file_r4(file_handle,rst_sol, &
                                               read_block_type)
      case (real128) ; call read_rstrt_file_r16(file_handle,rst_sol, &
                                                read_block_type)
      case default
        write (error_message,23)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end select
  end if
  !
  ! Finally close the restart file
  !
  call mpi_file_close(file_handle,mpierr)
  !
  ! Copy the restart solution data over into the solution array
  !
  solution_errors = 0
  !
  do this_cell = 1,ncell
    !
    inp_cell = cell_input_data(this_cell)%input_cell
    !
    l1 = cell_input_data(this_cell)%input_beg_sp
    l2 = cell_input_data(this_cell)%input_end_sp
    l0 = l1-1
    !
    inp_geom = cell_input_data(this_cell)%input_geom
    inp_order = cell_input_data(this_cell)%input_order
    !
    n1 = cell(this_cell)%beg_sp
    n2 = cell(this_cell)%end_sp
    n0 = n1-1
    np = n2-n0
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
    if (flip_sol_indices) then
      do n = 1,np
        do m = 1,restart_nq
          solution(n0+n,m) = rst_sol(m,l0+n)
        end do
      end do
    else
      do n = 1,np
        do m = 1,restart_nq
          solution(m,n0+n) = rst_sol(m,l0+n)
        end do
      end do
    end if
   !solution(:,n1:n2) = rst_sol(:,l1:l2)
    !
  end do
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
end subroutine read_solution_file
!
!###############################################################################
!
subroutine read_restart_old_format(file_name,solution,sol_is_transposed)
  !
  !.. Use Statements ..
  use geovar,       only : ncell,cell
  use geovar,       only : n_global_totpts
  use geovar,       only : global_pinc_ptr
  use ovar,         only : itrst,tmrst,dtrst
  use flowvar,      only : global_usp
  use parallel_mod, only : cell_map
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: file_name
  real(wp),      intent(inout) :: solution(:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: sol_is_transposed
  !
  !.. Local Scalars ..
  integer :: iorst,g0,n,m,nc,np,n0,n1,n2,ierr
  logical(lk) :: flip_sol_indices
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_restart_old_format"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  flip_sol_indices = fals
  if (present(sol_is_transposed)) then
    flip_sol_indices = sol_is_transposed
  end if
  !
  allocate ( global_usp(1:nq,1:n_global_totpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"global_usp",1,__LINE__,__FILE__,ierr,error_message)
  !
  open(newunit=iorst,file=file_name,status="old",action="read", &
                     form="unformatted",iostat=ierr,iomsg=error_message)
  call io_error(pname,file_name,1,__LINE__,__FILE__,ierr,error_message)
  !
  read (iorst) itrst, tmrst, dtrst
  read (iorst) global_usp
  !
  close (iorst,iostat=ierr,iomsg=error_message)
  call io_error(pname,file_name,2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Copy the solution data from only those cells
  ! that are local to this processor into solution
  ! NOTE: For a serial job, cell_map(mypnum+1)%loc_to_glb = [1,ncell]
  !
  do nc = 1,ncell
    !
    ! Indices of the data for this cell in the solution array
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
   !n1 = cell(nc)%beg_sp
   !n2 = cell(nc)%end_sp
   !np = n2-n1+1
    !
    ! Starting point of the data for this cell in the global_usp array
    !
    g0 = global_pinc_ptr( cell_map(mypnum+1)%loc_to_glb(nc) )
    !
    ! Copy the solution data from the global_usp array for this cell
    !
    if (flip_sol_indices) then
      do n = 1,np
        do m = 1,nq
          solution(n0+n,m) = global_usp(m,g0+n)
        end do
      end do
    else
      do n = 1,np
        do m = 1,nq
          solution(m,n0+n) = global_usp(m,g0+n)
        end do
      end do
    end if
   !solution(:,n1:n2) = global_usp(:,g0+1:g0+np)
    !
  end do
  !
  deallocate ( global_usp , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"global_usp",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_restart_old_format
!
!###############################################################################
!
subroutine write_restart_file(itnum)
  !
  !.. Use Statements ..
  use geovar,       only : n_global_cell
  use ovar,         only : itcur,time,d_t,rl2mx
  use ovar,         only : loc_solution_pts
  use ovar,         only : output_time_averaging
  use ovar,         only : restart_output_uses_host_roots
  use flowvar,      only : usp
  !.. Use Statements for old restart files ..
  use geovar,       only : n_global_totpts
  use ovar,         only : use_old_restart_format
  use flowvar,      only : global_usp
  use parallel_mod, only : collect_global_solution
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
  logical(lk) :: sol_is_transposed
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
  integer(wip) :: int_header(1:4)
  real(wp)     :: real_header(1:3)
  !
  !.. Local Allocatable Arrays ..
  integer(i4), allocatable :: ordered_geom(:)
  real(wp),    allocatable :: solution(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_restart_file"
  !
  logical(ldk), parameter :: wait_for_completion = .true.
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  sol_is_transposed = fals
  !
  if (use_old_restart_format) then
    !
    if (itnum > 0) then
      write (fullname,'("out.",i0.8,".rst")') itcur
    else if (itnum == -2) then
      write (fullname,'("out.",i0.8,".NaN.rst")') itcur-1
    else
      write (fullname,'("out.",i0.8,".done.rst")') itcur
    end if
    !
    if (mypnum == glb_root) then
      if (sol_is_transposed) then
        allocate ( global_usp(1:n_global_totpts,1:nq) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"global_usp",1,__LINE__,__FILE__,ierr, &
                         error_message)
      else
        allocate ( global_usp(1:nq,1:n_global_totpts) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"global_usp",1,__LINE__,__FILE__,ierr, &
                         error_message)
      end if
    end if
    !
    ! Collect the global solution
    !
    call collect_global_solution
    !
    ! Only use the root partition to output the restart file
    !
    if (mypnum == glb_root) then
      if (sol_is_transposed) then
        global_usp = transpose(global_usp) ! F2003 AUTO-REALLOCATION
      end if
      call write_restart_old_format(fullname,global_usp)
    end if
    !
    if (allocated(global_usp)) then
      deallocate ( global_usp , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"global_usp",2,__LINE__,__FILE__,ierr, &
                       error_message)
    end if
    !
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    !
  else
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
                                  ordered_geom=ordered_geom, &
                                  sol_is_transposed=sol_is_transposed)
    else
      involved_with_output = true
      restart_comm = MPI_COMM_WORLD
      call get_restart_arrays(current_sol=usp, &
                              output_sol=solution, &
                              ordered_geom=ordered_geom, &
                              sol_is_transposed=sol_is_transposed)
    end if
    !
    ! Store the size information for the simulation into the int_header array
    !
    int_header = [ int(nq,kind=wip) , &
                   int(n_global_cell,kind=wip) , &
                   int(restart_solution_form,kind=wip) , &
                   int(itcur,kind=wip) ]
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
      call write_solution_file(fullname,int_header,real_header,restart_comm, &
                               ordered_geom,ordered_geom_type, &
                               solution,solution_type,sol_block_type, &
                               statfile=statfile)
      !
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
      call write_time_averaged_solution(itnum,ordered_geom,int_header, &
                                        sol_is_transposed=sol_is_transposed)
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
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine write_restart_file
!
!###############################################################################
!
subroutine write_solution_file(file_name,int_header,real_header,file_comm, &
                               ordered_geom,og_type,solution,sol_type, &
                               sol_blktyp,var_header,statfile)
  !
  !.. Use Statements ..
  use ovar,    only : lustre_stripe_factor
  use ovar,    only : lustre_stripe_unit
  use ovar,    only : restart_output_uses_host_roots
  use ovar,    only : profile_io_restart
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: file_name
  integer,          intent(in) :: int_header(:)
  real(wp),         intent(in) :: real_header(:)
  _MPI_COMM_TYPE_,  intent(in) :: file_comm
  integer(i4),      intent(in) :: ordered_geom(:)
  _MPI_DATA_TYPE_,  intent(in) :: og_type
  real(wp),         intent(in) :: solution(:,:)
  _MPI_DATA_TYPE_,  intent(in) :: sol_type
  _MPI_DATA_TYPE_,  intent(in) :: sol_blktyp
  !
  !.. Optional Arguments ..
  integer,          optional, intent(in) :: var_header(:)
  character(len=*), optional, intent(in) :: statfile
  !
  !.. Local Scalars ..
  integer      :: ios
  integer(wip) :: disp,isize
  !
  _MPI_FILE_TYPE_ :: file_handle
  _MPI_INFO_TYPE_ :: file_info
  integer(MPI_OFFSET_KIND) :: file_disp
  integer(MPI_OFFSET_KIND) :: file_offset
  integer(kind(MPI_MODE_WRONLY)) :: file_mode
  !
  real(wp) :: beg_wall_time
  real(wp) :: end_wall_time
  !
  !.. Local Arrays ..
  _MPI_STATUS_TYPE_ :: istat(1:_MPI_STATUS_SIZE_)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_solution_file"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  call mpi_barrier(file_comm,mpierr)
  if (mypnum == glb_root) then
    write (iout,1)
    flush (iout)
    if (profile_io_restart) beg_wall_time = mpi_wtime()
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
  !
  ! Enable collective buffering for writes
  !
  call mpi_info_set(file_info,"romio_cb_write","enable",mpierr)
  !
  ! Specify the number of I/O devices that the
  ! restart file should be striped across
  !
  call mpi_info_set(file_info,"striping_factor", &
                    trim(adjustl(lustre_stripe_factor)),mpierr)
  !
  ! Specifiy the suggested striping unit to be used for this file. The
  ! striping unit is the amount of consecutive data assigned to one I/O
  ! device before progressing to the next device, when striping across
  ! a number of devices. It is expressed in bytes.
  !
  call mpi_info_set(file_info,"striping_unit", &
                    trim(adjustl(lustre_stripe_unit)),mpierr)
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
                     stat_strng="after mpi_info, before mpi_file_open")
  !
  ! Open the restart file for MPI I/O
  !
  file_mode = MPI_MODE_WRONLY + MPI_MODE_CREATE
  call mpi_file_open(file_comm,trim(adjustl(file_name)), &
                     file_mode,file_info,file_handle,mpierr)
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
                            stat_strng="after mpi_file_open")
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
                     stat_strng="after bcast disp from restart header")
  !
  ! Set the file view to write the ordered_geom array to the restart file
  !
  call mpi_file_set_view(file_handle,file_disp,MPI_INTEGER, &
                         og_type, trim(adjustl(file_datarep)), &
                         file_info, mpierr)
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                     stat_strng="after setting view for ordered_geom")
  !
  ! Write the ordered_geom array to the restart file
  !
  call mpi_file_write_all(file_handle,ordered_geom, &
                          size(ordered_geom,kind=int_mpi), &
                          MPI_INTEGER,istat(1),mpierr)
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            mpistat=istat, &
                            stat_strng="after writing ordered_geom")
  !
  ! Update the file displacement since we know the total number of cells
  ! and that the elements of the array ordered_geom are 4-byte integers.
  !
  file_disp = file_disp + int( 4*int_header(2) , kind=MPI_OFFSET_KIND )
  !
  ! Update the file view to after all the ordered_geom data
  !
  call mpi_file_set_view(file_handle,file_disp,MPI_INTEGER, &
                         MPI_INTEGER, trim(adjustl(file_datarep)), &
                         file_info, mpierr)
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            stat_strng="after setting view for wp")
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
                            stat_strng="after bcast disp from wp")
  !
  ! Convert the file displacement to MPI_OFFSET_KIND
  !
  file_disp = int( disp , kind=MPI_OFFSET_KIND )
  !
  ! Set the file fiew to write the solution data to the restart file
  !
  call mpi_file_set_view(file_handle,file_disp,sol_blktyp, &
                         sol_type, trim(adjustl(file_datarep)), &
                         file_info, mpierr)
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            stat_strng="after setting view for solution")
  !
  ! Write the solution data to the restart file
  !
  call mpi_file_write_all(file_handle,solution, &
                          size(solution,dim=2,kind=int_mpi), &
                          sol_blktyp,istat(1),mpierr)
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            mpistat=istat, &
                            stat_strng="after writing solution")
  !
  ! Finally close the restart file
  !
  call mpi_file_close(file_handle,mpierr)
  !
  call WRITE_RESTART_STATUS(ios,file_comm,file_name=statfile, &
                            close_file=true, &
                            stat_strng="after mpi_file_close")
  !
  ! Output the wall time taken to write the restart solution file if
  ! the input parameter profile_io_restart is true
  !
  call mpi_barrier(file_comm,mpierr)
  if (mypnum == glb_root) then
    if (profile_io_restart) then
      end_wall_time = mpi_wtime()
      write (iout,2) trim(adjustl(file_name)), &
                     end_wall_time - beg_wall_time, &
                     trim(adjustl(merge("ONLY HOST ROOTS  ", &
                                        "ALL MPI PROCESSES", &
                                        restart_output_uses_host_roots)))
    else
      write (iout,2)
    end if
    flush (iout)
  end if
  call mpi_barrier(file_comm,mpierr)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," Entering write_solution_file")
  2 format (" Leaving write_solution_file",/,:, &
          5x,"Wall time taken to write the file '",a, &
             "' = ",f12.6," seconds",/, &
         10x,"(USING ",a," FOR MPI-I/O)",/)
  !
end subroutine write_solution_file
!
!###############################################################################
!
subroutine write_restart_old_format(file_name,solution)
  !
  !.. Use Statements ..
  use ovar, only : itcur,time,d_t
  use ovar, only : profile_io_restart
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: file_name
  real(wp),         intent(in) :: solution(:,:)
  !
  !.. Local Scalars ..
  integer  :: iorst,ierr
  real(wp) :: beg_wall_time
  real(wp) :: end_wall_time
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_restart_old_format"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  write (iout,1)
  flush (iout)
  if (profile_io_restart) beg_wall_time = mpi_wtime()
  !
  open (newunit=iorst,file=trim(adjustl(file_name)),status="replace", &
                      action="write",form='unformatted',iostat=ierr, &
                      iomsg=error_message)
  call io_error(pname,trim(adjustl(file_name)),1,__LINE__,__FILE__, &
                ierr,error_message)
  !
  write (iorst) itcur, time, d_t
  write (iorst) solution
  !
  close (iorst,iostat=ierr,iomsg=error_message)
  call io_error(pname,trim(adjustl(file_name)),2,__LINE__,__FILE__, &
                ierr,error_message)
  !
  ! Output the wall time taken to write the restart solution file if
  ! the input parameter profile_io_restart is true
  !
  if (profile_io_restart) then
    end_wall_time = mpi_wtime()
    write (iout,2) trim(adjustl(file_name)), &
                   end_wall_time - beg_wall_time
  else
    write (iout,2)
  end if
  flush (iout)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," Entering write_restart_old_format")
  2 format (" Leaving write_restart_old_format",/,:, &
          5x,"Wall time taken to write the file '",a, &
             "' = ",f12.6," seconds",/, &
         10x,"(USING GLOBAL ROOT FOR SERIAL I/O)",/)
  !
end subroutine write_restart_old_format
!
!###############################################################################
!
subroutine write_time_averaged_solution(itnum,ordered_geom,int_header, &
                                        sol_is_transposed)
  !
  !.. Use Statements ..
  use ovar,    only : time,prev_ave_time,ave_start_time,itcur,rl2mx
  use ovar,    only : restart_output_uses_host_roots
  use flowvar, only : uavesp,time_ave_variables
  !
  !.. Formal Arguments ..
  integer,         intent(in) :: itnum
  integer(i4),     intent(in) :: ordered_geom(:)
  integer(wip), intent(inout) :: int_header(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: sol_is_transposed
  !
  !.. Local Scalars ..
  integer :: n,ierr
  integer :: exitstat,cmndstat
  !
  logical(lk) :: involved_with_output
  logical(lk) :: flip_sol_indices
  !
  character(len=100) :: fname
  character(len=200) :: fullname
  character(len=100) :: statfile
  character(len=300) :: command
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
  call debug_timer(entering_procedure,pname)
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
  flip_sol_indices = fals
  if (present(sol_is_transposed)) then
    flip_sol_indices = sol_is_transposed
  end if
  !
  ! Get the name for the current time averaged solution file
  !
  if (itnum > 0) then
    write (fname,tafname_fmt) averaging_file_switch,""
    write (fullname,tafile_fmt) averaging_file_switch,""
    write (command,talink_fmt) averaging_file_switch,""
    write (statfile,tastat_fmt) averaging_file_switch,""
  else
    write (fname,tadonefname_fmt) itcur,".done"
    write (fullname,tadone_fmt) itcur,".done"
    write (command,talkdn_fmt) itcur,".done"
    write (statfile,tastat_fmt) itcur,".done"
  end if
  !
  ! Create the time averaged information to write to the averaging file
  !
  if (restart_output_uses_host_roots) then
    involved_with_output = i_am_host_root
    timeave_comm = host_roots_comm
    call collect_restart_arrays(current_sol=uavesp, &
                                output_sol=solution, &
                                sol_is_transposed=flip_sol_indices)
  else
    involved_with_output = true
    timeave_comm = MPI_COMM_WORLD
    call get_restart_arrays(current_sol=uavesp, &
                            output_sol=solution, &
                            sol_is_transposed=flip_sol_indices)
  end if
  !
  ! Store the size information specific to the time averaged solution
  ! into the respective location within the int_header array
  !
  if (flip_sol_indices) then
    int_header(1) = size(uavesp,dim=2)
  else
    int_header(1) = size(uavesp,dim=1)
  end if
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
    call write_solution_file(fullname,int_header,real_header,timeave_comm, &
                             ordered_geom,ordered_geom_type, &
                             solution,averaging_type,ave_block_type, &
                             var_header=var_header,statfile=statfile)
    !
  end if
#ifdef USE_CREATE_SYMLINK
  !
  ! Create a link to the most recent restart file for piggy-backing jobs
  !
  if (mypnum == glb_root) then
    call create_restart_link(fname)
  end if
#else
  !
  ! Create a link to the most recent restart file for piggy-backing jobs
  !
  call execute_command(command,wait_for_completion,exitstat, &
                       cmndstat,cmndmsg,pname)
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
  !
  averaging_file_switch = merge(2,1,averaging_file_switch==1)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (i0)
  2 format (" Leaving write_time_averaged_solution",/,:, &
            5x,"Wall time taken to write time averaged solution = ", &
            f12.6," seconds")
  3 format (10x,"(USING ",a," FOR MPI-I/O)",/)
  4 format (/," Entering write_time_averaged_solution")
  !
end subroutine write_time_averaged_solution
!
!###############################################################################
!
subroutine get_restart_arrays(current_sol,output_sol,ordered_geom, &
                              sol_is_transposed)
  !
  !.. Use Statements ..
  use geovar,  only : ncell,cell,n_solpts
  !
  !.. Formal Arguments ..
  real(wp),                 intent(in) :: current_sol(:,:)
  real(wp), allocatable, intent(inout) :: output_sol(:,:)
  !
  !.. Optional Arguments ..
  integer(i4), optional, allocatable, intent(inout) :: ordered_geom(:)
  logical(lk), optional,                 intent(in) :: sol_is_transposed
  !
  !.. Local Scalars ..
  integer :: k,l0,l1,l2,m,n,nc,np,nvar,n0,n1,n2,ierr
  integer :: this_geom,this_order
  logical(lk) :: flip_sol_indices
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_restart_arrays"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  flip_sol_indices = fals
  if (present(sol_is_transposed)) then
    flip_sol_indices = sol_is_transposed
  end if
  !
  if (present(ordered_geom)) then
    !
    allocate ( ordered_geom(1:ncell) , source=0_i4 , &
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
      ordered_geom(n) = int( 100*this_order + this_geom , kind=i4 )
      !
    end do
    !
  end if
  !
  ! Allocate a temporary array to copy over only the interior solution
  ! data for each processor. We also have to copy over the solution data
  ! in the same order that is expected by the MPI derived datatype.
  !
  if (flip_sol_indices) then
    nvar = size(current_sol,dim=2)
  else
    nvar = size(current_sol,dim=1)
  end if
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
   !n1 = cell(nc)%beg_sp
   !n2 = cell(nc)%end_sp
   !n0 = n1-1
   !np = n2-n0
    !
    l0 = l2
    l2 = l2 + np
   !l1 = l2 + 1
   !l2 = l2 + n2-n1+1
    !
    if (flip_sol_indices) then
      do k = 1,np
        do m = 1,nvar
          output_sol(m,l0+k) = current_sol(n0+k,m)
        end do
      end do
    else
      do k = 1,np
        do m = 1,nvar
          output_sol(m,l0+k) = current_sol(m,n0+k)
        end do
      end do
    end if
   !output_sol(1:np,l1:l2) = current_sol(1:np,n1:n2)
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_restart_arrays
!
!###############################################################################
!
subroutine collect_restart_arrays(current_sol,output_sol,ordered_geom, &
                                  sol_is_transposed)
  !
  !.. Use Statements ..
  use geovar,  only : ncell,cell,n_solpts
  !
  !.. Formal Arguments ..
  real(wp),                 intent(in) :: current_sol(:,:)
  real(wp), allocatable, intent(inout) :: output_sol(:,:)
  !
  !.. Optional Arguments ..
  integer(i4), optional, allocatable, intent(inout) :: ordered_geom(:)
  logical(lk), optional,                 intent(in) :: sol_is_transposed
  !
  !.. Local Scalars ..
  integer :: nvar,ierr,lcell,lpts
  logical(lk) :: flip_sol_indices
  !
  !.. LOcal Allocatable Arrays ..
  integer(i4), allocatable :: local_og(:)
  real(wp),    allocatable :: local_sol(:,:)
  !
  integer(INT_MPI), allocatable :: sol_counts(:)
  integer(INT_MPI), allocatable :: sol_displs(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "collect_restart_arrays"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  flip_sol_indices = fals
  if (present(sol_is_transposed)) then
    flip_sol_indices = sol_is_transposed
  end if
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
    allocate ( ordered_geom(1:lcell) , source=0_i4 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ordered_geom",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    ! F2003 AUTO-REALLOCATION
#ifdef SPECIAL_FOR_GCC
    local_og = [ int( 100*cell(1:ncell)%order + cell(1:ncell)%geom , kind=i4 ) ]
#else
    local_og = int( 100*cell(1:ncell)%order + cell(1:ncell)%geom , kind=i4 )
#endif
    !
    call mpi_gatherv(local_og,int(ncell,kind=int_mpi),MPI_INTEGER, &
                     ordered_geom,og_counts,og_displs,MPI_INTEGER, &
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
  if (flip_sol_indices) then
    nvar = size(current_sol,dim=2)
  else
    nvar = size(current_sol,dim=1)
  end if
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
  if (flip_sol_indices) then
    ! F2003 AUTO-REALLOCATION
    local_sol = transpose(current_sol(1:n_solpts,1:nvar))
  else
    ! F2003 AUTO-REALLOCATION
    local_sol = current_sol(1:nvar,1:n_solpts)
  end if
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
  call debug_timer(leaving_procedure,pname)
  !
end subroutine collect_restart_arrays
!
!###############################################################################
!
pure function reorder_host_solution(solution) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: solution(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(solution,dim=1), &
                           1:size(solution,dim=2))
  !
  !.. Local Scalars ..
  integer :: nc,n1,n2,l1,l2
  !
continue
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
                                close_file,open_file,mpistat)
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
    if (present(stat_strng).or.present(mpistat)) then
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
  call debug_timer(entering_procedure,pname)
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
  ! Create the format string for creating the restart file names
  !
  write (char_num,1) itrst + num_timesteps
  final_iter_digits = max(4,len_trim(adjustl(char_num)))
  !
  write (rstfile_fmt,2) trim(adjustl(restart_dir)),final_iter_digits
  write (lnkcmnd_fmt,5) final_iter_digits,trim(adjustl(restart_dir))
  write (statfile_fmt,6) trim(adjustl(restart_dir)),final_iter_digits
  write (tafile_fmt,7) trim(adjustl(restart_dir))
  write (talink_fmt,8) trim(adjustl(restart_dir))
  write (tastat_fmt,9) trim(adjustl(restart_dir))
  write (tadone_fmt,10) trim(adjustl(restart_dir)),final_iter_digits
  write (talkdn_fmt,11) final_iter_digits,trim(adjustl(restart_dir))
  !
  write (rstfname_fmt,20) final_iter_digits
  write (tafname_fmt,21)
  write (tadonefname_fmt,22) final_iter_digits
  !
  ! Make sure restart_output_uses_host_roots is
  ! false if only one processor is being used.
  !
  if (ncpu <= 1) then
    restart_output_uses_host_roots = fals
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1  format (i0)
  2  format ("('",a,"out.',i0.",i0,",a,'.rst')")
  3  format ("lfs setstripe -c ",i0," -S ",a,1x,a)
  4  format ("Lustre striping command = '",a,"'")
  5  format ("('/bin/ln -sf out.',i0.",i0,",a,'.rst ",a,"out.current.rst')")
  6  format ("('",a,"out.',i0.",i0,",a,'.stat')")
  7  format ("('",a,"out.',i0,a,'.ave')")
  8  format ("('/bin/ln -sf out.',i0,a,'.ave ",a,"out.current.ave')")
  9  format ("('",a,"out.',i0,a,'.ave.stat')")
  10 format ("('",a,"out.',i0.",i0,",a,'.ave')")
  11 format ("('/bin/ln -sf out.',i0.",i0,",a,'.ave ",a,"out.current.ave')")
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
  20 format ("('out.',i0.",i0,",a,'.rst')")  ! rstfile_fmt
  21 format ("('out.',i0,a,'.ave')")         ! tafile_fmt
  22 format ("('out.',i0.",i0,",a,'.ave')")  ! tadone_fmt
  !
end subroutine check_restart_dir
!
!###############################################################################
!
subroutine create_restart_fname_formats
  !
  !.. Use Statements ..
  use ovar, only : itcur,num_timesteps
  !
  !.. Local Scalars ..
  integer :: final_iter_digits
  character(len=100) :: char_num
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
  ! NOTE: This procedure should be called after initializing the solution from
  !       a restart file or by other means and itcur is given the value of
  !       itrst from the restart file or 0 if not using a restart file.
  !
  write (char_num,1) itcur + num_timesteps
  final_iter_digits = max(4,len_trim(adjustl(char_num)))
  !
  ! Create the format strings for creating the restart file names and commands
  !
  write (rstfile_fmt,2) trim(adjustl(restart_dir)),final_iter_digits
  write (lnkcmnd_fmt,3) final_iter_digits,trim(adjustl(restart_dir))
  write (statfile_fmt,4) trim(adjustl(restart_dir)),final_iter_digits
  write (tafile_fmt,5) trim(adjustl(restart_dir))
  write (talink_fmt,6) trim(adjustl(restart_dir))
  write (tastat_fmt,7) trim(adjustl(restart_dir))
  write (tadone_fmt,8) trim(adjustl(restart_dir)),final_iter_digits
  write (talkdn_fmt,9) final_iter_digits,trim(adjustl(restart_dir))
  !
  write (rstfname_fmt,20) final_iter_digits
  write (tafname_fmt,21)
  write (tadonefname_fmt,22) final_iter_digits
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (i0)
  2 format ("('",a,"out.',i0.",i0,",a,'.rst')")
  3 format ("('/bin/ln -sf out.',i0.",i0,",a,'.rst ",a,"out.current.rst')")
  4 format ("('",a,"out.',i0.",i0,",a,'.stat')")
  5 format ("('",a,"out.',i0,a,'.ave')")
  6 format ("('/bin/ln -sf out.',i0,a,'.ave ",a,"out.current.ave')")
  7 format ("('",a,"out.',i0,a,'.ave.stat')")
  8 format ("('",a,"out.',i0.",i0,",a,'.ave')")
  9 format ("('/bin/ln -sf out.',i0.",i0,",a,'.ave ",a,"out.current.ave')")
  !
  20 format ("('out.',i0.",i0,",a,'.rst')")  ! rstfile_fmt
  21 format ("('out.',i0,a,'.ave')")         ! tafile_fmt
  22 format ("('out.',i0.",i0,",a,'.ave')")  ! tadone_fmt
  !
end subroutine create_restart_fname_formats
!
!###############################################################################
!
subroutine create_restart_link(fname)
  !
  !.. Use Statements ..
  use iso_c_binding, only : c_int,c_null_char
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: fname
  !
  !.. Local Scalars ..
  integer :: l,n
  integer(c_int) :: ierr
  logical(ldk) :: file_exists
  character(len=300) :: existing_file
  character(len=300) :: existing_file_null
  character(len=300) :: link_file
  character(len=300) :: link_file_null
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_restart_link"
  !
  !.. C Function Interfaces ..
  interface
    integer(c_int) function unlink(link_file) &
                               bind(c,name="unlink")
      use iso_c_binding, only : c_int,c_char
      implicit none
      character(kind=c_char), dimension(*), intent(in) :: link_file
    end function unlink
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
  existing_file = trim(adjustl(fname))
  existing_file_null = trim(adjustl(existing_file)) // c_null_char
  existing_file_null = trim(adjustl(existing_file_null))
  !
  n = len_trim(existing_file)
  !
  if (uppercase(existing_file(n-2:n)) == "AVE") then
    link_file = trim(adjustl(restart_dir)) // "out.current.ave"
  else
    link_file = trim(adjustl(restart_dir)) // "out.current.rst"
  end if
  link_file = trim(adjustl(link_file))
  link_file_null = trim(adjustl(link_file)) // c_null_char
  link_file_null = trim(adjustl(link_file_null))
  !
  do l = 1,5
    !
    inquire (file=link_file,exist=file_exists)
    if (file_exists) then
      ierr = unlink(link_file_null)
      if (ierr /= 0) then
        write (iout,1) trim(adjustl(link_file)),ierr
        flush (iout)
      end if
    end if
    !
    ierr = symlink(existing_file_null,link_file_null)
    !
    if (ierr /= 0) then
      write (iout,2) trim(adjustl(existing_file)), &
                     trim(adjustl(link_file)),ierr
      flush (iout)
    else
      exit
    end if
    !
  end do
  !
  1 format (/,2x,"ERROR RETURNED FROM C FUNCTION unlink!!!",/, &
              6x,"    link_file = '",a,"'",/, &
              6x,"         ierr = ",i0,/)
  2 format (/,2x,"ERROR RETURNED FROM C FUNCTION symlink!!!",/, &
              6x,"existing_file = '",a,"'",/, &
              6x,"    link_file = '",a,"'",/, &
              6x,"         ierr = ",i0,/)
  !
end subroutine create_restart_link
!
!###############################################################################
!
subroutine get_host_output_mpi_datatypes()
  !
  !.. Use Statement ..
  use geovar,  only : ncell,cell,n_solpts
  use geovar,  only : global_pinc_ptr
  use flowvar, only : uavesp
  use parallel_mod, only : cell_map
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
  call debug_timer(entering_procedure,pname)
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
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_host_output_mpi_datatypes
!
!###############################################################################
!
subroutine get_output_mpi_datatypes()
  !
  !.. Use Statement ..
  use geovar,  only : ncell,cell
  use geovar,  only : global_pinc_ptr
  use flowvar, only : uavesp
  use parallel_mod, only : cell_map
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
  call debug_timer(entering_procedure,pname)
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
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_output_mpi_datatypes
!
!###############################################################################
!
subroutine get_host_time_ave_mpi_datatypes()
  !
  !.. Use Statement ..
  use flowvar, only : uavesp
  !
  !.. Local Scalars ..
  integer :: naq,ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_host_time_ave_mpi_datatypes"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
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
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_host_time_ave_mpi_datatypes
!
!###############################################################################
!
subroutine get_time_ave_mpi_datatypes()
  !
  !.. Use Statement ..
  use flowvar, only : uavesp
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
  call debug_timer(entering_procedure,pname)
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
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_time_ave_mpi_datatypes
!
!###############################################################################
!
subroutine get_input_mpi_datatypes(restart_nq,restart_rp,real_type, &
                                   read_block_type,read_type, &
                                   num_pts,ordered_geom,cell_input_data)
  !
  !.. Use Statements ..
  use geovar,       only : ncell ! ncell is recomputed after cell_map is created
  use parallel_mod, only : cell_map
  !
  !.. Formal Arguments ..
  integer(wip),                     intent(in) :: restart_nq
  integer(wip),                     intent(in) :: restart_rp
  _MPI_DATA_TYPE_,               intent(inout) :: real_type
  _MPI_DATA_TYPE_,               intent(inout) :: read_block_type
  _MPI_DATA_TYPE_,               intent(inout) :: read_type
  integer(wip),                  intent(inout) :: num_pts
  integer(wip),                     intent(in) :: ordered_geom(:)
  type(rstrt_cell_t), allocatable, intent(out) :: cell_input_data(:)
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,ierr
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
    this_order = ordered_geom(global_cell) / 100
    this_geom = ordered_geom(global_cell) - 100*this_order
    !
    pinc_ptr(global_cell+1) = pinc_ptr(global_cell) + &
                              Cell_Solpts(this_geom,this_order)
    !
  end do
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
  allocate ( msk(1:ncell) , source=true , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Sort the sol_offset array to find the output order for the cells so that
  ! the offsets for the MPI datatypes are monotonically increasing.
  ! NOTE: This array is assumed to be small enough so that it doesnt require
  !       some fancy sorting algorithm and this simple minloc method will
  !       more than likely be relatively inexpensive.
  !
  do nc = 1,ncell
    n = minloc( sol_offset , dim=1 , mask=msk )
    msk(n) = fals
    cell_input_data(nc)%input_cell = n
  end do
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
subroutine check_restart_solution_form(restart_solution_form)
  !
  !.. Use Statements ..
  use ovar, only : loc_solution_pts
  !
  !.. Formal Arguments ..
  integer, intent(in) :: restart_solution_form
  !
  !.. Local Scalars ..
  character(len=10) :: inp_solpts
  character(len=50) :: restart_form
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_restart_solution_form"
  !
continue
  !
  if (restart_solution_form /= -loc_solution_pts) then
    !
    select case (loc_solution_pts)
      case (Legendre_Gauss)         ; inp_solpts = "Gauss"
      case (Legendre_Gauss_Lobatto) ; inp_solpts = "Lobatto"
      case default                  ; inp_solpts = "Unknown"
    end select
    !
    select case (restart_solution_form)
      case (restart_gauss_nodal_solution)
        restart_form = "nodal form for the Gauss solution points"
      case (restart_lobatto_nodal_solution)
        restart_form = "nodal form for the Lobatto solution points"
      case (restart_legendre_modal_solution)
        restart_form = "modal form for Legendre basis functions"
      case (restart_monomial_modal_solution)
        restart_form = "modal form for monomial basis functions"
      case default
        restart_form = "Unknown"
    end select
    !
    write (error_message,1) trim(adjustl(restart_form)), &
                            trim(adjustl(inp_solpts))
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    !
  end if
  !
  ! Format Statements
  !
  1 format("The format of the solution data in the restart file does not ", &
           "match the current simulation. The restart data is in ",a," ", &
           "and the current simulation requires the data to be defined at ", &
           " the ",a," points.")
  !
end subroutine check_restart_solution_form
!
!###############################################################################
!###############################################################################
! MPI read subroutine for 8-bit INTEGER data
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_i1(file_handle,int_data)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  integer(wip), intent(inout) :: int_data(:)
  !
  !.. Local Arrays ..
  integer(int8) :: read_int(1:size(int_data))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_i1"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  call mpi_file_read_all(file_handle,read_int, &
                         size(read_int,kind=int_mpi), &
                         MPI_INTEGER1,MPI_STATUS_IGNORE,mpierr)
  !
  int_data = int( read_int , kind=wip )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_i1
!
!###############################################################################
!###############################################################################
! MPI read subroutine for 16-bit INTEGER data
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_i2(file_handle,int_data)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  integer(wip), intent(inout) :: int_data(:)
  !
  !.. Local Arrays ..
  integer(int16) :: read_int(1:size(int_data))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_i2"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  call mpi_file_read_all(file_handle,read_int, &
                         size(read_int,kind=int_mpi), &
                         MPI_INTEGER2,MPI_STATUS_IGNORE,mpierr)
  !
  int_data = int( read_int , kind=wip )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_i2
!
!###############################################################################
!###############################################################################
! MPI read subroutine for 32-bit INTEGER data
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_i4(file_handle,int_data)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  integer(wip), intent(inout) :: int_data(:)
  !
  !.. Local Arrays ..
  integer(int32) :: read_int(1:size(int_data))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_i4"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  call mpi_file_read_all(file_handle,read_int, &
                         size(read_int,kind=int_mpi), &
                         MPI_INTEGER4,MPI_STATUS_IGNORE,mpierr)
  !
  int_data = int( read_int , kind=wip )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_i4
!
!###############################################################################
!###############################################################################
! MPI read subroutine for 64-bit INTEGER data
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_i8(file_handle,int_data)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  integer(wip), intent(inout) :: int_data(:)
  !
  !.. Local Arrays ..
  integer(int64) :: read_int(1:size(int_data))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_i8"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  call mpi_file_read_all(file_handle,read_int, &
                         size(read_int,kind=int_mpi), &
                         MPI_INTEGER8,MPI_STATUS_IGNORE,mpierr)
  !
  int_data = int( read_int , kind=wip )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_i8
!
!###############################################################################
!###############################################################################
! MPI read subroutine for INTEGER data matching the working integer precision
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_wip(file_handle,int_data)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  integer(wip), intent(inout) :: int_data(:)
  !
  !.. Local Arrays ..
  integer(wip) :: read_int(1:size(int_data))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_rstrt_file_wip"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  call mpi_file_read_all(file_handle,read_int, &
                         size(int_data,kind=int_mpi), &
                         mpi_inttyp,MPI_STATUS_IGNORE,mpierr)
  !
  int_data = read_int
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_wip
!
!###############################################################################
!###############################################################################
! MPI read subroutines for 32-bit REAL data
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_r4_rank1(file_handle,real_data)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  real(wp),     intent(inout) :: real_data(:)
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
  call mpi_file_read_all(file_handle,read_real, &
                         size(read_real,kind=int_mpi), &
                         MPI_REAL4,MPI_STATUS_IGNORE,mpierr)
  !
  real_data = real( read_real , kind=wp )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_r4_rank1
!-------------------------------------------------------------------------------
subroutine read_rstrt_file_r4_rank2(file_handle,real_data,read_block_type)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  real(wp),     intent(inout) :: real_data(:,:)
  _MPI_DATA_TYPE_, intent(in) :: read_block_type
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
  call mpi_file_read_all(file_handle,read_real, &
                         size(read_real,dim=2,kind=int_mpi), &
                         read_block_type,MPI_STATUS_IGNORE,mpierr)
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
subroutine read_rstrt_file_r8_rank1(file_handle,real_data)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  real(wp),     intent(inout) :: real_data(:)
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
  call mpi_file_read_all(file_handle,read_real, &
                         size(read_real,kind=int_mpi), &
                         MPI_REAL8,MPI_STATUS_IGNORE,mpierr)
  !
  real_data = real( read_real , kind=wp )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_r8_rank1
!-------------------------------------------------------------------------------
subroutine read_rstrt_file_r8_rank2(file_handle,real_data,read_block_type)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  real(wp),     intent(inout) :: real_data(:,:)
  _MPI_DATA_TYPE_, intent(in) :: read_block_type
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
  call mpi_file_read_all(file_handle,read_real, &
                         size(read_real,dim=2,kind=int_mpi), &
                         read_block_type,MPI_STATUS_IGNORE,mpierr)
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
subroutine read_rstrt_file_r16_rank1(file_handle,real_data)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  real(wp),     intent(inout) :: real_data(:)
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
  call mpi_file_read_all(file_handle,read_real, &
                         size(read_real,kind=int_mpi), &
                         MPI_REAL16,MPI_STATUS_IGNORE,mpierr)
  !
  real_data = real( read_real , kind=wp )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_r16_rank1
!-------------------------------------------------------------------------------
subroutine read_rstrt_file_r16_rank2(file_handle,real_data,read_block_type)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  real(wp),     intent(inout) :: real_data(:,:)
  _MPI_DATA_TYPE_, intent(in) :: read_block_type
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
  call mpi_file_read_all(file_handle,read_real, &
                         size(read_real,dim=2,kind=int_mpi), &
                         read_block_type,MPI_STATUS_IGNORE,mpierr)
  !
  real_data = real( read_real , kind=wp )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_r16_rank2
!
!###############################################################################
!###############################################################################
! MPI read subroutines for REAL data matching the current working precision (wp)
!###############################################################################
!###############################################################################
!
subroutine read_rstrt_file_wp_rank1(file_handle,real_data)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  real(wp),     intent(inout) :: real_data(:)
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
  call mpi_file_read_all(file_handle,read_real, &
                         size(real_data,kind=int_mpi), &
                         mpi_flttyp,MPI_STATUS_IGNORE,mpierr)
  !
  real_data = read_real
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_wp_rank1
!-------------------------------------------------------------------------------
subroutine read_rstrt_file_wp_rank2(file_handle,real_data,read_block_type)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_, intent(in) :: file_handle
  real(wp),     intent(inout) :: real_data(:,:)
  _MPI_DATA_TYPE_, intent(in) :: read_block_type
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
  call mpi_file_read_all(file_handle,read_real, &
                         size(real_data,dim=2,kind=int_mpi), &
                         read_block_type,MPI_STATUS_IGNORE,mpierr)
  !
  real_data = read_real
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_rstrt_file_wp_rank2
!
!###############################################################################
!###############################################################################
!###############################################################################
!
include "Functions/cell_solpts.f90"
#include "Functions/sort.f90"
!
!###############################################################################
!
end module restart_mod
