module restart_mod
  !
#include <mpi_defs.h>
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
  logical, save :: restart_needs_initialization = true
  !
  integer, save :: restart_solution_form
  !
  _MPI_DATA_TYPE_, save :: ordered_geom_type = MPI_DATATYPE_NULL
  _MPI_DATA_TYPE_, save :: solution_type = MPI_DATATYPE_NULL
  !
  character(len=170), save :: restart_dir = ""
  character(len=170), save :: rstfile_fmt
  !
  !###################################
  !###  Module Allocatable Arrays  ###
  !###################################
  !
  integer, save, allocatable :: cell_output_order(:)
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
  use geovar,       only : ncell,cell,n_solpts
  use geovar,       only : n_global_cell
  use ovar,         only : itrst,tmrst,dtrst
  use ovar,         only : iop_solpts
  use flowvar,      only : usp
  !.. Use Statements for old restart files ..
  use geovar,       only : n_global_totpts
  use geovar,       only : global_pinc_ptr
  use ovar,         only : use_old_restart_format
  use flowvar,      only : global_usp
  use parallel_mod, only : cell_map
  !
  !.. Local Scalars ..
  integer :: n,nc,np,n1,n2,l1,l2,gp,ierr
  integer :: iorst
  !
  !.. Need to be default integers
  integer :: restart_ip,restart_rp
  !.. Need to be MPI integers
  _MPI_FILE_TYPE_ :: file_handle
  _MPI_INFO_TYPE_ :: file_info
  _MPI_DATA_TYPE_ :: read_type
  _MPI_DATA_TYPE_ :: real_type
  !.. Need to be MPI_OFFSET_KIND integers
  integer(MPI_OFFSET_KIND) :: file_disp
  integer(MPI_OFFSET_KIND) :: file_offset
  !
  integer :: solution_errors,num_pts
  integer :: restart_nq,restart_ncells
  integer :: n_header_int,n_header_real
  !
  integer :: inp_cell,inp_geom,inp_order
  integer :: this_cell,this_geom,this_order
  !
  character(len=10) :: inp_solpts
  character(len=50) :: restart_form
  !
  !.. Local Allocatable Arrays ..
  integer(wip),       allocatable :: ordered_geom(:)
  integer(wip),       allocatable :: int_header(:)
  real(wp),           allocatable :: real_header(:)
  real(wp),           allocatable :: solution(:,:)
  type(rstrt_cell_t), allocatable :: cell_input_data(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_restart_file"
  character(len=*), parameter :: fname = "out.rst"
  !
continue
  !
  solution_errors = 0
  !
  if (use_old_restart_format) then
    !
    allocate ( global_usp(1:nq,1:n_global_totpts) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_usp",1,ierr,error_message)
    !
    open(newunit=iorst,file=fname,status="old",action="read", &
                       form="unformatted",iostat=ierr)
    call io_error(pname,fname,1,ierr)
    !
    read (iorst) itrst, tmrst, dtrst
    read (iorst) global_usp
    !
    close (iorst,iostat=ierr)
    call io_error(pname,fname,2,ierr)
    !
    ! If using multiple processors, copy the solution data from 
    ! only those cells that are local to this processor into usp
    !
    if (ncpu > 1) then
      !
      do nc = 1,ncell
        !
        ! Indices of the data for this cell in the usp array
        !
        n1 = cell(nc)%beg_sp
        n2 = cell(nc)%end_sp
        np = n2-n1+1
        !
        ! Starting point of the data for this cell in the global_usp array
        !
        gp = global_pinc_ptr( cell_map(mypnum+1)%loc_to_glb(nc) )
        !
        ! Copy the solution data from the global_usp array for this cell
        !
        usp(:,n1:n2) = global_usp(:,gp+1:gp+np)
        !
      end do
      !
    else
      !
      ! For a serial job, just copy the global_usp array into usp
      !
      usp(:,:) = global_usp(:,:)
      !
    end if
    !
    deallocate ( global_usp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_usp",2,ierr,error_message)
    !
  else
    !
   !if (ncpu > 1) then
      !
#ifdef PBS_ENV
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
      ! Specify the number of I/O devices that the
      ! restart file should be striped across
      !
     !call mpi_info_set(file_info,"striping_factor","1",mpierr)
      !
      ! Specifiy the suggested striping unit to be used for this file. The
      ! striping unit is the amount of consecutive data assigned to one I/O
      ! device before progressing to the next device, when striping across
      ! a number of devices. It is expressed in bytes.
      !
     !call mpi_info_set(file_info,"striping_unit","1",mpierr)
      !
#else
      !
      ! Set file_info to MPI_INFO_NULL
      !
      file_info = MPI_INFO_NULL
      !
#endif
      !
      ! Open the restart file for MPI I/O
      !
      call mpi_file_open(MPI_COMM_WORLD, trim(adjustl(fname)), &
                         MPI_MODE_RDONLY,file_info,file_handle,mpierr)
      !
      ! Have all processors read the header section of the restart file
      !
      ! Read the integer type for the size information for the simulation
      !
      call mpi_file_read_all(file_handle,restart_ip,1_int_mpi, &
                             MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
      !
      ! Read the number of size data items
      !
      call mpi_file_read_all(file_handle,n_header_int,1_int_mpi, &
                             MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
      !
      allocate ( int_header(1:n_header_int) , source=0_wip , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"int_header",1,ierr,error_message)
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
            call stop_recon(stop_mpi,pname,error_message)
        end select
      end if
      !
      ! Extract the information from the integer section of the header
      !
      restart_nq = int( int_header(1) , kind=kind(restart_nq) )
      !
      if (restart_nq /= nq) then
        write (error_message,10) restart_nq,nq
        call stop_recon(stop_mpi,pname,error_message)
      end if
      !
      restart_ncells = int( int_header(2) , kind=kind(restart_ncells) )
      !
      if (restart_ncells /= n_global_cell) then
        write (error_message,11) restart_ncells,n_global_cell
        call stop_recon(stop_mpi,pname,error_message)
      end if
      !
      restart_solution_form = int( int_header(3) , &
                                   kind=kind(restart_solution_form) )
      !
      if (restart_solution_form /= -iop_solpts) then
        !
        select case (iop_solpts)
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
        write (error_message,12) trim(adjustl(restart_form)), &
                                 trim(adjustl(inp_solpts))
        call stop_recon(stop_mpi,pname,error_message)
        !
      end if
      !
      itrst = int( int_header(4) , kind=kind(itrst) )
      !
      deallocate ( int_header , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"int_header",2,ierr,error_message)
      !
      ! Read the real type for the time data
      !
      call mpi_file_read_all(file_handle,restart_rp,1_int_mpi, &
                             MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
      !
      ! Read the number of time data items
      !
      call mpi_file_read_all(file_handle,n_header_real,1_int_mpi, &
                             MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
      !
      allocate ( real_header(1:n_header_real) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"real_header",1,ierr,error_message)
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
            call stop_recon(stop_mpi,pname,error_message)
        end select
      end if
      !
      ! Extract the information from the real section of the header
      !
      tmrst  = real( real_header(1) , kind=kind(tmrst) )
      dtrst  = real( real_header(2) , kind=kind(dtrst) )
      !
      deallocate ( real_header , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"real_header",2,ierr,error_message)
      !
      ! Read the integer type for the ordered_geom array
      !
      call mpi_file_read_all(file_handle,restart_ip,1_int_mpi, &
                          MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
      !
      allocate ( ordered_geom(1:restart_ncells) , source=0_wip , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"ordered_geom",1,ierr,error_message)
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
            call stop_recon(stop_mpi,pname,error_message)
        end select
      end if
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
      ! of the restart file
      !
      call get_input_mpi_datatypes(restart_nq,restart_rp,real_type,read_type, &
                                   num_pts,ordered_geom,cell_input_data)
      !
      ! IF WE ARE TO MOVE THE POINT IN THE PRE-PROCESSING PHASE AT WHICH THE
      ! RESTART FILE IS READ, WE NEED TO MAKE SURE THAT N_SOLPTS HAS BEEN
      ! RECOMPUTED SO THAT IT IS LOCAL TO EACH PROCESSOR BECAUSE THIS IS NOT
      ! DONE WITHIN PARALLEL_MOD AFTER PARTITIONING THE GRID.
      !
      if (num_pts /= n_solpts) then
        write (iout,13) mypnum,num_pts,n_solpts
        solution_errors = 100
      end if
      !
      ! Sum up the value of solution_errors across all processors
      !
      call mpi_allreduce(MPI_IN_PLACE,solution_errors,1_int_mpi, &
                         mpi_inttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
      !
      if (solution_errors /= 0) then
        write (error_message,14)
        call stop_recon(stop_mpi,pname,error_message)
      end if
      !
      ! Set the file view to write the ordered_geom array to the restart file
      !
      call mpi_file_set_view(file_handle,file_disp,real_type, &
                             read_type, trim(adjustl(file_datarep)), &
                             file_info, mpierr)
      !
      ! Allocate the solution array to read in the solution restart data
      !
      allocate ( solution(1:nq,1:n_solpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"solution",1,ierr,error_message)
      !
      ! Read in the solution restart data
      !
      if (restart_rp == wp) then
        call read_rstrt_file_wp(file_handle,solution)
      else
        select case (restart_rp)
          case (real64)  ; call read_rstrt_file_r8(file_handle,solution)
          case (real32)  ; call read_rstrt_file_r4(file_handle,solution)
          case (real128) ; call read_rstrt_file_r16(file_handle,solution)
          case default
            write (error_message,23)
            call stop_recon(stop_mpi,pname,error_message)
        end select
      end if
      !
      ! Finally close the restart file
      !
      call mpi_file_close(file_handle,mpierr)
      !
      ! Copy the restart solution data over into usp
      !
      solution_errors = 0
      !
      do this_cell = 1,ncell
        !
        inp_cell = cell_input_data(this_cell)%input_cell
        !
        l1 = cell_input_data(this_cell)%input_beg_sp
        l2 = cell_input_data(this_cell)%input_end_sp
        !
        inp_geom = cell_input_data(this_cell)%input_geom
        inp_order = cell_input_data(this_cell)%input_order
        !
        n1 = cell(this_cell)%beg_sp
        n2 = cell(this_cell)%end_sp
        !
        this_geom = cell(this_cell)%geom
        this_order = cell(this_cell)%order
        !
        if ((inp_geom /= this_geom) .or. (inp_order /= this_order)) then
          write (iout,30) mypnum, inp_cell, inp_geom, inp_order, &
                                 this_cell,this_geom,this_order
          solution_errors = solution_errors + 1
        end if
        !
        usp(1:nq,n1:n2) = solution(1:nq,l1:l2)
        !
      end do
      !
      ! Sum up the value of solution_errors across all processors
      !
      call mpi_allreduce(MPI_IN_PLACE,solution_errors,1_int_mpi, &
                         mpi_inttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
      !
      ! Stop the simulation if there were any errors copying the restart
      ! solution data into usp.
      !
      if (solution_errors > 0) then
        write (error_message,31) solution_errors
        call stop_recon(stop_mpi,pname,error_message)
      end if
      !
      ! Deallocate any local arrays that are still allocated
      !
      if (allocated(cell_input_data)) then
        deallocate ( cell_input_data , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"cell_input_data",2,ierr,error_message)
      end if
      !
      if (allocated(ordered_geom)) then
        deallocate ( ordered_geom , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"ordered_geom",2,ierr,error_message)
      end if
      !
      if (allocated(solution)) then
        deallocate ( solution , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"solution",2,ierr,error_message)
      end if
      !
      !
   !else
      !
      ! Stuff to read the data in serial
      !
   !end if
    !
  end if
  !
  ! Format Statements
  !
  10 format("The number of equations in the restart file does not match the ", &
            "current simulation: restart_nq=",i0,", current_nq=",i0,".")
  11 format("The number of global grid cells in the restart file does not ", &
            "match the current simulation: restart_n_global_cell=",i0, &
            ", current_n_global_cell=",i0,".")
  12 format("The format of the solution data in the restart file does not ", &
            "match the current simulation. The restart data is in ",a," ", &
            "and the current simulation requires the data to be defined at ", &
            " the ",a," points.")
  13 format(/,"  The number of solution points for the cells on CPU-",i0," ", &
              "does not match the",/,"   number of solution points for the ", &
              "corresponding cells in the restart file: ",/, &
              "        restart_n_solpts=",i0,/, &
              "        current_n_solpts=",i0)
  14 format("The number of solution points for cells in the restart file ", &
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
end subroutine read_restart_file
!
!###############################################################################
!
subroutine write_restart_file(itnum)
  !
  !.. Use Statements ..
  use geovar,  only : cell,ncell,n_solpts
  use geovar,  only : n_global_cell
  use ovar,    only : itcur,time,d_t
  use ovar,    only : iop_solpts
  use flowvar, only : usp
  !.. Use Statements for old restart files ..
  use geovar,  only : n_global_totpts
  use ovar,    only : use_old_restart_format
  use flowvar, only : global_usp
  use parallel_mod, only : collect_global_solution
  !
  !.. Formal Arguments ..
  integer, intent(in) :: itnum
  !
  !.. Local Scalars ..
  integer :: ierr,iorst
  integer :: n,nc,n1,n2,l1,l2
  integer :: this_geom,this_order
  integer(wip) :: disp
  !
  _MPI_FILE_TYPE_ :: file_handle
  _MPI_INFO_TYPE_ :: file_info
  integer(MPI_OFFSET_KIND) :: file_disp
  integer(MPI_OFFSET_KIND) :: file_offset
  !
  character(len=40) :: fname
  !
  !.. Local Arrays ..
  integer(wip) :: int_header(1:4)
  real(wp)     :: real_header(1:2)
  !
  !.. Local Allocatable Arrays ..
  integer(i4), allocatable :: ordered_geom(:)
  real(wp),    allocatable :: solution(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_restart_file"
  !
continue
  !
  if (use_old_restart_format) then
    !
    if (mypnum == 0) then
      allocate ( global_usp(1:nq,1:n_global_totpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"global_usp",1,ierr,error_message)
    end if
    !
    ! Collect the global solution
    !
    call collect_global_solution
    !
    ! Only use the root partition to output the restart file
    !
    if (mypnum == 0) then
      !
      if (itnum > 0) then
        write (fname,'("out.",i0.8,".rst")') itcur
      else
        write (fname,'("out.",i0.8,".done.rst")') itcur
      end if
      !
      open (newunit=iorst,file=trim(adjustl(fname)),status="replace", &
                          action="write",form='unformatted',iostat=ierr)
      call io_error(pname,trim(adjustl(fname)),1,ierr)
      !
      write (iorst) itcur, time, d_t
      write (iorst) global_usp
      !
      close (iorst,iostat=ierr)
      call io_error(pname,trim(adjustl(fname)),2,ierr)
      !
      deallocate ( global_usp , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"global_usp",2,ierr,error_message)
      !
    end if
    !
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    !
  else
    !
    ! Perform some initialization tasks if the restart module has not
    ! been initialized yet.
    !
    if (restart_needs_initialization) then
      !
      ! Check for the restart directory and create it if it does not exist
      !
      call check_restart_dir
      !
      ! For now, write the nodal solution to the restart file. At some point we
      ! should switch so it writes the modal solution since that would give more
      ! flexibility regarding changes in the solution order between runs.
      !
      if (iop_solpts == Legendre_Gauss) then
        restart_solution_form = restart_gauss_nodal_solution
      else if (iop_solpts == Legendre_Gauss_Lobatto) then
        restart_solution_form = restart_lobatto_nodal_solution
      end if
      !
      ! Get the MPI derived types required to read/write the restart file
      !
      call get_output_mpi_datatypes
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
      write (fname,rstfile_fmt) itcur,""
    else
      write (fname,rstfile_fmt) itcur,".done"
    end if
    !
    ! Create the cell information to write to the restart file
    !
    allocate ( ordered_geom(1:ncell) , source=0_i4 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ordered_geom",1,ierr,error_message)
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
#ifdef PBS_ENV
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
    ! Specify the number of I/O devices that the
    ! restart file should be striped across
    !
   !call mpi_info_set(file_info,"striping_factor","1",mpierr)
    !
    ! Specifiy the suggested striping unit to be used for this file. The
    ! striping unit is the amount of consecutive data assigned to one I/O
    ! device before progressing to the next device, when striping across
    ! a number of devices. It is expressed in bytes.
    !
   !call mpi_info_set(file_info,"striping_unit","1",mpierr)
    !
#else
    !
    ! Set file_info to MPI_INFO_NULL
    !
    file_info = MPI_INFO_NULL
    !
#endif
    !
    ! Open the restart file for MPI I/O
    !
    call mpi_file_open(MPI_COMM_WORLD,trim(adjustl(fname)), &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       file_info,file_handle,mpierr)
    !
    ! Have the root processor write the header for the restart file
    !
    if (mypnum == 0) then
      !
      ! Write the integer type for the size information for the simulation
      !
      call mpi_file_write(file_handle,wip,1_int_mpi, &
                          MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
      !
      ! Write the number of size data items
      !
      call mpi_file_write(file_handle,size(int_header,kind=int_mpi), &
                          1_int_mpi,MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
      !
      ! Store the size information for the simulation into the temporary array
      ! int_header so that we dont have to call mpi_file_write multiple times.
      !
      int_header = [ int(nq,kind=wip) , &
                     int(n_global_cell,kind=wip) , &
                     int(restart_solution_form,kind=wip) , &
                     int(itcur,kind=wip) ]
      !
      ! Write the size information for the simulation to the restart file
      !
      call mpi_file_write(file_handle,int_header, &
                          size(int_header,kind=int_mpi), &
                          mpi_inttyp,MPI_STATUS_IGNORE,mpierr)
      !
      ! Write the real type for the time data
      !
      call mpi_file_write(file_handle,wp,1_int_mpi, &
                          MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
      !
      ! Write the number of time data items
      !
      call mpi_file_write(file_handle,size(real_header,kind=int_mpi), &
                          1_int_mpi,MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
      !
      ! Store the time data into the temporary array real_header so
      ! that we dont have to call mpi_file_write multiple times.
      !
      real_header = [ time , d_t ]
      !
      ! Write the time data to the restart file
      !
      call mpi_file_write(file_handle,real_header, &
                          size(real_header,kind=int_mpi), &
                          mpi_flttyp,MPI_STATUS_IGNORE,mpierr)
      !
      ! Write the integer type for the ordered_geom array
      !
      call mpi_file_write(file_handle,i4,1_int_mpi, &
                          MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
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
    ! Broadcast the current file displacement to all processors
    !
    call mpi_bcast(disp,1_int_mpi,mpi_inttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
    !
    ! Convert the file displacement to MPI_OFFSET_KIND
    !
    file_disp = int( disp , kind=MPI_OFFSET_KIND )
    !
    ! Set the file view to write the ordered_geom array to the restart file
    !
    call mpi_file_set_view(file_handle,file_disp,MPI_INTEGER, &
                           ordered_geom_type, trim(adjustl(file_datarep)), &
                           file_info, mpierr)
    !
    ! Write the ordered_geom array to the restart file
    !
    call mpi_file_write_all(file_handle,ordered_geom, &
                            size(ordered_geom,kind=int_mpi), &
                            MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
    !
    ! Update the file displacement since we know the total number of cells
    ! and that the elements of the array ordered_geom are 4-byte integers.
    !
    file_disp = file_disp + int( 4*n_global_cell , kind=MPI_OFFSET_KIND )
    !
    ! Update the file view to after all the ordered_geom data 
    !
    call mpi_file_set_view(file_handle,file_disp,MPI_INTEGER, &
                           MPI_INTEGER, trim(adjustl(file_datarep)), &
                           file_info, mpierr)
    !
    ! Have the root processor write the byte size for the solution data
    !
    if (mypnum == 0) then
      !
      ! Write the kind value for real working-precision
      !
      call mpi_file_write(file_handle,wp,1_int_mpi, &
                          MPI_INTEGER,MPI_STATUS_IGNORE,mpierr)
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
    ! Broadcast the current file displacement to all processors
    !
    call mpi_bcast(disp,1_int_mpi,mpi_inttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
    !
    ! Convert the file displacement to MPI_OFFSET_KIND
    !
    file_disp = int( disp , kind=MPI_OFFSET_KIND )
    !
    ! Set the file fiew to write the solution data to the restart file
    !
    call mpi_file_set_view(file_handle,file_disp,mpi_flttyp, &
                           solution_type, trim(adjustl(file_datarep)), &
                           file_info, mpierr)
    !
    ! Allocate a temporary array to copy over only the interior solution
    ! data for each processor. We also have to copy over the solution data
    ! in the same order that is expected by the MPI derived datatype.
    !
    allocate ( solution(1:nq,1:n_solpts) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"solution",1,ierr,error_message)
    ! 
    l2 = 0
    do n = 1,ncell
      !
      nc = cell_output_order(n)
      !
      n1 = cell(nc)%beg_sp
      n2 = cell(nc)%end_sp
      !
      l1 = l2 + 1
      l2 = l2 + n2-n1+1
      !
      solution(1:nq,l1:l2) = usp(1:nq,n1:n2)
      !
    end do
    !
    ! Write the solution data to the restart file
    !
    call mpi_file_write_all(file_handle,solution, &
                            size(solution,kind=int_mpi), &
                            mpi_flttyp,MPI_STATUS_IGNORE,mpierr)
    !
    ! Finally close the restart file
    !
    call mpi_file_close(file_handle,mpierr)
    !
    ! Deallocate ordered_geom and solution since they are no longer needed
    !
    deallocate ( ordered_geom , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ordered_geom",2,ierr,error_message)
    !
    deallocate ( solution , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"solution",2,ierr,error_message)
    !
  end if
  !
end subroutine write_restart_file
!
!###############################################################################
!
subroutine check_restart_dir
  !
  !.. Use Statements ..
  use ovar, only : output_dir,itrst,num_ts
  use ovar, only : lustre_stripe_count,lustre_stripe_size
  !
  !.. Local Scalars ..
  integer :: clen,ierr,iodir
  integer :: exitstat,cmndstat
  integer :: final_iter_digits
  logical :: using_restart_dir
  character(len=300) :: command
  character(len=300) :: cmndmsg
  character(len=200) :: test_file
  character(len=200) :: test_dir
  character(len=100) :: char_num
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_restart_dir"
  !
  logical, parameter :: wait_for_completion = true
  !
continue
  !
  ! Initialize exitstat and cmndmsg
  !
  exitstat = 0
  cmndmsg = ""
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
  if (mypnum == 0) then
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
      call execute_command_line( trim(command) , wait_for_completion , &
                                 exitstat , cmndstat , cmndmsg )
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
     !command = "lfs setstripe -c 64 -s 4m "//trim(restart_dir)
      write (command,3) max(0,lustre_stripe_count), &
                        trim(adjustl(lustre_stripe_size)), &
                        trim(adjustl(restart_dir))
      !
      call execute_command_line( trim(command) , wait_for_completion , &
                                 exitstat , cmndstat , cmndmsg )
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
  call mpi_bcast(clen,1_int_mpi,mpi_inttyp, &
                 0_int_mpi,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(restart_dir(1:clen),int(clen,kind=int_mpi),MPI_CHARACTER, &
                 0_int_mpi,MPI_COMM_WORLD,mpierr)
  !
  ! Create the format string for creating the restart file names
  !
  write (char_num,1) itrst + num_ts
  final_iter_digits = max(4,len_trim(adjustl(char_num)))
  !
  write (rstfile_fmt,2) trim(adjustl(restart_dir)),final_iter_digits
  !
  ! Format Statements
  !
  1 format (i0)
  2 format ("('",a,"out.',i0.",i0,",a,'.rst')")
  3 format ("lfs setstripe -c ",i0," -s ",a,1x,a)
  !
end subroutine check_restart_dir
!
!###############################################################################
!
subroutine get_output_mpi_datatypes()
  !
  !.. Use Statement ..
  use geovar,  only : ncell,cell
  use geovar,  only : global_pinc_ptr
  use parallel_mod, only : cell_map
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,ierr
  integer :: local_cell
  integer :: global_cell
  !
  !.. Local Allocatable Arrays ..
  logical,          allocatable :: msk(:)
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
  ! Allocate the local arrays
  !
  allocate ( og_offset(1:ncell) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"og_offset",1,ierr,error_message)
  !
  allocate ( sol_offset(1:ncell) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sol_offset",1,ierr,error_message)
  !
  allocate ( sol_blklen(1:ncell) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sol_blklen",1,ierr,error_message)
  !
  ! Create the cell offsets as well as the block lengths
  ! and offsets for the solution data
  !
  do local_cell = 1,ncell
    !
    ! Indices of the data for this cell in the usp array
    !
    n1 = cell(local_cell)%beg_sp
    n2 = cell(local_cell)%end_sp
    !
    sol_blklen(local_cell) = int( n2-n1+1 , kind=int_mpi )
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
    sol_offset(local_cell) = int( global_pinc_ptr(global_cell) , kind=int_mpi )
    !
  end do
  !
  ! Dont forget to account for the number of equations in the offsets
  ! and block lengths for the solution data
  !
  sol_offset = sol_offset * int( nq , kind=int_mpi )
  sol_blklen = sol_blklen * int( nq , kind=int_mpi )
  !
  ! Allocate the cell_output_order array and
  ! the mask array for sorting the cell offsets
  !
  allocate ( cell_output_order(1:ncell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_output_order",1,ierr,error_message)
  !
  allocate ( msk(1:ncell) , source=true , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",1,ierr,error_message)
  !
  ! Sort the og_offset array to find the output order for the cells so that
  ! the offsets for the MPI datatypes are monotonically increasing.
  ! NOTE: This array is assumed to be small enough so that it doesnt require
  !       some fancy sorting algorithm and this simple minloc method will
  !       more than likely be relatively inexpensive.
  !
  do nc = 1,ncell
    n = minloc( og_offset , dim=1 , mask=msk )
    msk(n) = fals
    cell_output_order(nc) = n
  end do
  !
  ! Deallocate the msk array
  !
  deallocate ( msk , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",2,ierr,error_message)
  !
  ! Allocate the temporary array for applying the sorting to the offset
  ! and block length arrays.
  !
  allocate ( tmp_copy(1:ncell) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_copy",1,ierr,error_message)
  !
  ! Apply the sort to the og_offset array
  !
  tmp_copy = og_offset
  og_offset = tmp_copy(cell_output_order)
  !
  ! Apply the sort to the sol_offset array
  !
  tmp_copy = sol_offset
  sol_offset = tmp_copy(cell_output_order)
  !
  ! Apply the sort to the sol_blklen array
  !
  tmp_copy = sol_blklen
  sol_blklen = tmp_copy(cell_output_order)
  !
  ! Deallocate the temporary array
  !
  deallocate ( tmp_copy , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_copy",2,ierr,error_message)
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
  call mpi_type_indexed(size(sol_offset,kind=int_mpi), &
                        sol_blklen,sol_offset, &
                        mpi_flttyp,solution_type,mpierr)
  call mpi_type_commit(solution_type,mpierr)
  !
  ! Deallocate the local arrays before leaving
  !
  deallocate ( og_offset , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"og_offset",2,ierr,error_message)
  !
  deallocate ( sol_offset , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sol_offset",2,ierr,error_message)
  !
  deallocate ( sol_blklen , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sol_blklen",2,ierr,error_message)
  !
end subroutine get_output_mpi_datatypes
!
!###############################################################################
!
subroutine get_input_mpi_datatypes(restart_nq,restart_rp,real_type,read_type, &
                                   num_pts,ordered_geom,cell_input_data)
  !
  !.. Use Statements ..
  use geovar,       only : ncell ! ncell is recomputed after cell_map is created
  use parallel_mod, only : cell_map
  !
  !.. Formal Arguments ..
  integer(wip),                                   intent(in) :: restart_nq
  integer(wip),                                   intent(in) :: restart_rp
  _MPI_DATA_TYPE_,                            intent(inout) :: real_type
  _MPI_DATA_TYPE_,                            intent(inout) :: read_type
  integer(wip),                                intent(inout) :: num_pts
  integer(wip),                    dimension(:),  intent(in) :: ordered_geom
  type(rstrt_cell_t), allocatable, dimension(:), intent(out) :: cell_input_data
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,ierr
  integer :: local_cell,global_cell
  integer :: this_cell,this_geom,this_order
  integer :: inp_cell,inp_geom,inp_order
  !
  !.. Local Allocatable Arrays ..
  logical,          allocatable :: msk(:)
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
  call alloc_error(pname,"pinc_ptr",1,ierr,error_message)
  !
  allocate ( sol_offset(1:ncell) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sol_offset",1,ierr,error_message)
  !
  allocate ( sol_blklen(1:ncell) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sol_blklen",1,ierr,error_message)
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
  call alloc_error(pname,"pinc_ptr",2,ierr,error_message)
  !
  ! Allocate the cell_input_data array and
  ! the mask array for sorting the cell offsets
  !
  if (allocated(cell_input_data)) then
    deallocate ( cell_input_data , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_input_data",2,ierr,error_message)
  end if
  allocate ( cell_input_data(1:ncell) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_input_data",1,ierr,error_message)
  !
  allocate ( msk(1:ncell) , source=true , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",1,ierr,error_message)
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
  call alloc_error(pname,"msk",2,ierr,error_message)
  !
  ! Allocate the temporary array for applying the sorting to the offset
  ! and block length arrays.
  !
  allocate ( tmp_copy(1:ncell) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_copy",1,ierr,error_message)
  !
  ! Apply the sort to the sol_offset array
  !
  tmp_copy = sol_offset
  sol_offset = tmp_copy(cell_input_data(:)%input_cell)
  !
  ! Apply the sort to the sol_blklen array
  !
  tmp_copy = sol_blklen
  sol_blklen = tmp_copy(cell_input_data(:)%input_cell)
  !
  ! Deallocate the temporary array
  !
  deallocate ( tmp_copy , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_copy",2,ierr,error_message)
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
  sol_offset = sol_offset * int( restart_nq , kind=int_mpi )
  sol_blklen = sol_blklen * int( restart_nq , kind=int_mpi )
  !
  ! Create the MPI derived datatype for reading the
  ! solution data to the restart file
  !
  call mpi_type_indexed(size(sol_offset,kind=int_mpi), &
                        sol_blklen,sol_offset, &
                        real_type,read_type,mpierr)
  call mpi_type_commit(read_type,mpierr)
  !
  ! Deallocate the local arrays before leaving
  !
  deallocate ( sol_offset , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sol_offset",2,ierr,error_message)
  !
  deallocate ( sol_blklen , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sol_blklen",2,ierr,error_message)
  !
  !
end subroutine get_input_mpi_datatypes
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
  _MPI_FILE_TYPE_,              intent(in) :: file_handle
  integer(wip), dimension(:), intent(inout) :: int_data
  !
  !.. Local Arrays ..
  integer(int8), dimension(1:size(int_data)) :: read_int
  !
continue
  !
  call mpi_file_read_all(file_handle,read_int, &
                         size(read_int,kind=int_mpi), &
                         MPI_INTEGER1,MPI_STATUS_IGNORE,mpierr)
  !
  int_data = int( read_int , kind=wip )
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
  _MPI_FILE_TYPE_,              intent(in) :: file_handle
  integer(wip), dimension(:), intent(inout) :: int_data
  !
  !.. Local Arrays ..
  integer(int16), dimension(1:size(int_data)) :: read_int
  !
continue
  !
  call mpi_file_read_all(file_handle,read_int, &
                         size(read_int,kind=int_mpi), &
                         MPI_INTEGER2,MPI_STATUS_IGNORE,mpierr)
  !
  int_data = int( read_int , kind=wip )
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
  _MPI_FILE_TYPE_,              intent(in) :: file_handle
  integer(wip), dimension(:), intent(inout) :: int_data
  !
  !.. Local Arrays ..
  integer(int32), dimension(1:size(int_data)) :: read_int
  !
continue
  !
  call mpi_file_read_all(file_handle,read_int, &
                         size(read_int,kind=int_mpi), &
                         MPI_INTEGER4,MPI_STATUS_IGNORE,mpierr)
  !
  int_data = int( read_int , kind=wip )
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
  _MPI_FILE_TYPE_,              intent(in) :: file_handle
  integer(wip), dimension(:), intent(inout) :: int_data
  !
  !.. Local Arrays ..
  integer(int64), dimension(1:size(int_data)) :: read_int
  !
continue
  !
  call mpi_file_read_all(file_handle,read_int, &
                         size(read_int,kind=int_mpi), &
                         MPI_INTEGER8,MPI_STATUS_IGNORE,mpierr)
  !
  int_data = int( read_int , kind=wip )
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
  _MPI_FILE_TYPE_,              intent(in) :: file_handle
  integer(wip), dimension(:), intent(inout) :: int_data
  !
continue
  !
  call mpi_file_read_all(file_handle,int_data, &
                         size(int_data,kind=int_mpi), &
                         mpi_inttyp,MPI_STATUS_IGNORE,mpierr)
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
  _MPI_FILE_TYPE_,          intent(in) :: file_handle
  real(wp), dimension(:), intent(inout) :: real_data
  !
  !.. Local Arrays ..
  real(real32), dimension(1:size(real_data)) :: read_real
  !
continue
  !
  call mpi_file_read_all(file_handle,read_real, &
                         size(read_real,kind=int_mpi), &
                         MPI_REAL4,MPI_STATUS_IGNORE,mpierr)
  !
  real_data = real( read_real , kind=wp )
  !
end subroutine read_rstrt_file_r4_rank1
!-------------------------------------------------------------------------------
subroutine read_rstrt_file_r4_rank2(file_handle,real_data)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_,            intent(in) :: file_handle
  real(wp), dimension(:,:), intent(inout) :: real_data
  !
  !.. Local Arrays ..
  real(real32) :: read_real(1:size(real_data,dim=1),1:size(real_data,dim=2))
  !
continue
  !
  call mpi_file_read_all(file_handle,read_real, &
                         size(read_real,kind=int_mpi), &
                         MPI_REAL4,MPI_STATUS_IGNORE,mpierr)
  !
  real_data = real( read_real , kind=wp )
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
  _MPI_FILE_TYPE_,          intent(in) :: file_handle
  real(wp), dimension(:), intent(inout) :: real_data
  !
  !.. Local Arrays ..
  real(real64), dimension(1:size(real_data)) :: read_real
  !
continue
  !
  call mpi_file_read_all(file_handle,read_real, &
                         size(read_real,kind=int_mpi), &
                         MPI_REAL8,MPI_STATUS_IGNORE,mpierr)
  !
  real_data = real( read_real , kind=wp )
  !
end subroutine read_rstrt_file_r8_rank1
!-------------------------------------------------------------------------------
subroutine read_rstrt_file_r8_rank2(file_handle,real_data)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_,            intent(in) :: file_handle
  real(wp), dimension(:,:), intent(inout) :: real_data
  !
  !.. Local Arrays ..
  real(real64) :: read_real(1:size(real_data,dim=1),1:size(real_data,dim=2))
  !
continue
  !
  call mpi_file_read_all(file_handle,read_real, &
                         size(read_real,kind=int_mpi), &
                         MPI_REAL8,MPI_STATUS_IGNORE,mpierr)
  !
  real_data = real( read_real , kind=wp )
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
  _MPI_FILE_TYPE_,          intent(in) :: file_handle
  real(wp), dimension(:), intent(inout) :: real_data
  !
  !.. Local Arrays ..
  real(real128), dimension(1:size(real_data)) :: read_real
  !
continue
  !
  call mpi_file_read_all(file_handle,read_real, &
                         size(read_real,kind=int_mpi), &
                         MPI_REAL16,MPI_STATUS_IGNORE,mpierr)
  !
  real_data = real( read_real , kind=wp )
  !
end subroutine read_rstrt_file_r16_rank1
!-------------------------------------------------------------------------------
subroutine read_rstrt_file_r16_rank2(file_handle,real_data)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_,            intent(in) :: file_handle
  real(wp), dimension(:,:), intent(inout) :: real_data
  !
  !.. Local Arrays ..
  real(real128) :: read_real(1:size(real_data,dim=1),1:size(real_data,dim=2))
  !
continue
  !
  call mpi_file_read_all(file_handle,read_real, &
                         size(read_real,kind=int_mpi), &
                         MPI_REAL16,MPI_STATUS_IGNORE,mpierr)
  !
  real_data = real( read_real , kind=wp )
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
  _MPI_FILE_TYPE_,          intent(in) :: file_handle
  real(wp), dimension(:), intent(inout) :: real_data
  !
continue
  !
  call mpi_file_read_all(file_handle,real_data, &
                         size(real_data,kind=int_mpi), &
                         mpi_flttyp,MPI_STATUS_IGNORE,mpierr)
  !
end subroutine read_rstrt_file_wp_rank1
!-------------------------------------------------------------------------------
subroutine read_rstrt_file_wp_rank2(file_handle,real_data)
  !
  !.. Formal Arguments ..
  _MPI_FILE_TYPE_,            intent(in) :: file_handle
  real(wp), dimension(:,:), intent(inout) :: real_data
  !
continue
  !
  call mpi_file_read_all(file_handle,real_data, &
                         size(real_data,kind=int_mpi), &
                         mpi_flttyp,MPI_STATUS_IGNORE,mpierr)
  !
end subroutine read_rstrt_file_wp_rank2
!
!###############################################################################
!###############################################################################
!###############################################################################
!
include "Functions/cell_solpts.f90"
!
!###############################################################################
!
end module restart_mod
