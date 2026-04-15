module write_cgns_mod
  !
#include <mpi_defs.h>
  !
  !.. Use Statements ..
  use module_kind_types
  use cgns_types_mod
  use eqn_idx_mod, only : nec,nmb,nme,nmx,nmy,nmz,nee,nq
  use eqn_idx_mod, only : ncs,npres,ntemp,ngamm,nenth,ncp,nvisc,ncond,ng
  use eqn_idx_mod, only : NCD_Axial,NCD_Vert,NCD_Span
  use eqn_idx_mod, only : NM_Axial,NM_Vert,NM_Span
  !
  implicit none
  !
  private
  !
  ! ##################################
  ! ###  Module Public Procedures  ###
  ! ##################################
  !
  public :: write_parallel_cgns
  public :: TaylorGreen_8s_vorticity_solution
  public :: convert_restart_files_to_cgns
  !
  ! ################################
  ! ###  Module Local Variables  ###
  ! ################################
  !
  !.. Output File Names ..
  character(len=*), parameter :: cgns_base_fname = "solution.grid.cgns"
  character(len=*), parameter :: cgns_grid_fname = "solution.grid.cgns"
  character(len=*), parameter :: cgns_time_fname = "solution.time.cgns"
  character(len=170),    save :: cgns_grid_file
  character(len=170),    save :: cgns_time_file
  character(len=170),    save :: cgnsfile_fmt
  !
  !.. CGNS Time Dependent Variables ..
  integer, save :: cgns_zone_counter
  !
  real(sp),             save, allocatable :: soltime(:)
  integer,              save, allocatable :: soliter(:)
  character(len=CGLEN), save, allocatable :: solname(:)
  character(len=CGLEN), save, allocatable :: solfile(:)
  !
  ! CGNS Base and Zone names (hard-coded for now)
  character(len=*), parameter :: base_name = "gfr_base"
  character(len=*), parameter :: zone_name = "Solution"
  !
  ! Derived type for CGNS sections
  !
  type :: cgns_sect_t
    ! NAME : This is the string name for this CGNS section
    character(len=CGLEN) :: name
    ! TYPE : Element type for this section
    integer(CET) :: type
    ! ISECT_G : CGNS section index
    integer(CBT) :: isect_g = 0_CBT
    ! ALL_CELLS_BEG : Index of first element in this section (global)
    integer(CST) :: all_cells_beg
    ! ALL_CELLS_END : Index of last element in this section (global)
    integer(CST) :: all_cells_end
    ! MY_CELLS_BEG : Index of first element in this section (local)
    integer(CST) :: my_cells_beg
    ! MY_CELLS_END : Index of last element in this section (local)
    integer(CST) :: my_cells_end
    ! MAXOFFSET : This is just the size of ElemConn or the last value in the
    !             array Offsets. This is just a helpful integer to have
    !             when writing a MIXED CGNS section instead of having to
    !             extract it from ElemConn or Offsets.
    integer(CST) :: MaxOffset
    ! ELEMCONN : This is the actual element connectivities.
    ! NOTE: For a MIXED section, each elements connectivity is prepended
    !       by the element type, e.g., the connectivity for a linear hex
    !       is 9 integers: [HEXA_8, node 1, ..., node 8].
    integer(CST), allocatable :: ElemConn(:)
    ! OFFSETS : This is only needed for a MIXED section. This gives the
    !           start/offset of each element within the ElemConn array.
    integer(CST), allocatable :: Offsets(:)
  end type cgns_sect_t
  !
  ! CGNS section information
  ! NOTE: This needs to be an array for when using a mixed element grid
  !       and writing out each element type as its own section.
  !
  type(cgns_sect_t), save, allocatable :: cgns_sect(:)
  !
  ! Derived type containing various information for
  ! each combination of cell geom/order
  !
  type :: geord_t
    ! USED,NOT_USED : logicals indicating if this geom/order combo exists
    logical(lk) :: used = fals
    logical(lk) :: not_used = true
    ! NCELL : number of global cells of this geom/order combo
    integer :: ncell = 0
    ! TYPE : CGNS element type for this geom/order combo
    integer(CET) :: type
    ! NPTS : number of output nodes/points of this element type
    integer :: npts
    ! PTS : element connectivity / point indices
    integer(CST), allocatable :: pts(:)
  end type geord_t
  !
  ! Array of derived types to store various information for
  ! the different combinations of cell geom/order
  !
  type(geord_t), save, allocatable :: geord(:,:)
  !
  !
  integer, save, allocatable :: cell_k0(:)
  !
  logical(lk), save :: interpolated_dusp_is_needed
  !
  character(len=*), parameter :: char_xyz(1:3) = ["X","Y","Z"]
  character(len=*), parameter :: char_lc_xyz(1:3) = ["x","y","z"]
  character(len=*), parameter :: char_cyl(1:3) = ["Axial","R    ","Theta"]
  !
  !.. CGNS Parallel Point Indices ..
  integer(CST), save :: my_pts_beg
  integer(CST), save :: my_pts_end
  !
  integer(CST), dimension(1:3) :: zone_size
  !
  ! Generic interface for interpolate_for_output
  !
  interface interpolate_for_output
    module procedure interpolate_for_output_1d, &
                     interpolate_for_output_2d, &
                     interpolate_for_output_3d
  end interface interpolate_for_output
  !
contains
!
!###############################################################################
!
subroutine convert_restart_files_to_cgns()
  !
  !.. Use Statements ..
  use ovar,          only : restart_list_file,restart_file
  use ovar,          only : output_time_averaging,read_time_ave_restart
  use ovar,          only : ave_start_time
  use ovar,          only : time,itcur,tmrst,itrst
  use ovar,          only : time_ref,time_scaling_factor
  use ovar,          only : delete_converted_restart_files
  use geomvar_mod,   only : n_solpts
  use flowvar_mod,   only : usp,gpsp
  use restart_mod,   only : read_restart_file
  use restart_mod,   only : deallocate_solver_arrays
  use flux_mod,      only : discontinuous_gradient_at_solpts
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Local Scalars ..
  integer  :: ierr,io,i,n
  real(wp) :: out_time
  !
  logical(ldk) :: file_exists
  !
  character(len=150) :: text
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: fname
  !
  !.. Local Allocatable Arrays ..
  character(len=300), allocatable :: restart_file_list(:)
  character(len=300), allocatable :: temp_list(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "convert_restart_files_to_cgns"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! First deallocate unnecessary arrays that were allocated for the solver but
  ! wont be needed since execution will abort at the end of this subroutine
  !
  call deallocate_solver_arrays(keep_cgns_arrays=true)
  !
  ! Write out the CGNS grid file and first solution file using
  ! using the restart file provided in the input file
  !
  if (mypnum == glb_root) then
    write (iout,2) trim(adjustl(restart_file)),1,1
  end if
  !
  if (mypnum == glb_root) then
    out_time = time*time_ref*time_scaling_factor
    if (output_time_averaging) then
      out_time = out_time - ave_start_time*time_ref*time_scaling_factor
    end if
    write (iout,3) out_time
  end if
  !
  call write_parallel_cgns(open_CGNS)
  !
  if (len_trim(adjustl(restart_list_file)) == 0) then
    write (error_message,5)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  !###############################################
  !################    FOR NOW    ################
  !###############################################
  !
  ! Disable reading in of a time averaged restart files and
  ! writing the CGNS time averaged solution files
  !
  read_time_ave_restart = fals
  output_time_averaging = fals
  !
  !###############################################
  !###############################################
  !###############################################
  !
  ! Read in the list of restart files that will be averaged
  !
  fname = trim(adjustl(restart_list_file)) ! F2003 AUTO-REALLOC
  !
  open (newunit=io,file=fname,status="old",action="read", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Count the number of restart files given in the list of files
  !
  n = 0
  do
    read (io,1,iostat=ierr,iomsg=error_message) text
    if (ierr < 0) then
      exit
    else if (ierr > 0) then
      call io_error(pname,fname,-1,__LINE__,__FILE__,ierr,error_message)
    else
      n = n + 1
    end if
  end do
  !
  rewind (io)
  !
  ! Allocate the arrays for storing the file names from the list of files
  !
  allocate ( restart_file_list(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"restart_file_list",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( temp_list(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_list",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  ! Read the file names from the list of files
  !
  do n = 1,size(restart_file_list)
    read (io,1,iostat=ierr,iomsg=error_message) restart_file_list(n)
    call io_error(pname,fname,-1,__LINE__,__FILE__,ierr,error_message)
  end do
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Check the validity of the filenames in the restart file list.
  !   1) Ignore any filenames that start with !
  !   2) For those filenames not starting with !, use inquire to check whether
  !      or not the file exists. If the file does exist, add the filename
  !      to temp_list.
  !
  i = 0
  do n = 1,size(restart_file_list)
    fname = trim(adjustl(restart_file_list(n))) ! F2003 AUTO-REALLOC
    if (mypnum == glb_root) write (iout,10) n,trim(adjustl(fname))
    if (fname(1:1) /= "!") then
      if (mypnum == glb_root) write (iout,11) trim(adjustl(fname))
      inquire (file=fname,exist=file_exists)
      if (file_exists) then
        i = i + 1
        if (mypnum == glb_root) write (iout,12) i
        temp_list(i) = fname
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
  ! Use move_alloc to over-write restart_file_list with temp_list
  !
  call move_alloc( from=temp_list , to=restart_file_list )
  !
  ! Loop over the list of files to read each restart
  ! file and write out a CGNS solution file
  !
  do n = 1,size(restart_file_list)
    !
    restart_file = trim(adjustl(restart_file_list(n)))
    !
    if (mypnum == glb_root) then
      write (iout,2) trim(adjustl(restart_file)),n,size(restart_file_list)
    end if
    !
    call read_restart_file
    !
    time = tmrst
    itcur = itrst
    !
    if (mypnum == glb_root) then
      out_time = time*time_ref*time_scaling_factor
      if (output_time_averaging) then
        out_time = out_time - ave_start_time*time_ref*time_scaling_factor
      end if
      write (iout,3) out_time
    end if
    !
    call discontinuous_gradient_at_solpts
    !
    call debug_timer(entering_procedure,"update_thermo_properties")
    call update_thermo_properties(1,n_solpts,usp,gpsp,caller_name=pname)
    call debug_timer(leaving_procedure,"update_thermo_properties")
    !
    call debug_timer(entering_procedure,"update_transport_properties")
    call update_transport_properties(1,n_solpts,usp,gpsp)
    call debug_timer(leaving_procedure,"update_transport_properties")
    !
    if (n == size(restart_file_list)) then
      call write_parallel_cgns(close_CGNS)
    else
      call write_parallel_cgns(add_to_CGNS)
    end if
    !
    if (delete_converted_restart_files) then
      if (mypnum == glb_root) then
        !
        open (newunit=io,file=restart_file,status="old",action="readwrite", &
                         access="stream",form="unformatted", &
                         iostat=ierr,iomsg=error_message)
        call io_error(pname,restart_file,1,__LINE__,__FILE__,ierr,error_message)
        !
        close (io,status="delete",iostat=ierr,iomsg=error_message)
        call io_error(pname,restart_file,2,__LINE__,__FILE__,ierr,error_message)
        !
      end if
    end if
    !
  end do
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
              3x,"Working on restart file '",a,"' (",i3," of ",i3,")",//, &
              3(1x,100("#"),/),/, &
              4x,"Reading the restart file")
  3 format (/,1x,81("-"),/, &
              1x,"-----  Creating CGNS solution file from restart file at ", &
                 "tstar = ",f10.4,"  -----",/, &
              1x,81("-"),/)
  4 format ("STOPPING EXECUTION AFTER CONVERTING THE LIST OF ", &
            "RESTART FILES TO CGNS SOLUTION FILES!")
  5 format ("STOPPING EXECUTION AFTER CONVERTING THE RESTART ", &
            "FILE TO A CGNS SOLUTION FILE!")
  !
  10 format (/,2x,"Filename # ",i0," from restart file list = '",a,"'")
  11 format (4x,"Checking existance of restart file '",a,"'")
  12 format (4x,"Restart file found! Saving to temp_list(",i0,")")
  13 format (8x,"ERROR: RESTART FILE DOES NOT EXIST!")
  14 format (8x,"WARNING: IGNORING FILENAME BEGINNING WITH '!'...")
  !
end subroutine convert_restart_files_to_cgns
!
!###############################################################################
!
subroutine write_parallel_cgns(iflag)
  !
  !.. Use Statements ..
  use ovar, only : itcur,itrst,num_timesteps,time,time_ref
  use ovar, only : write_TaylorGreen_full_solution
  use ovar, only : completely_disable_cgns
  use ovar, only : output_time_averaging
  use ovar, only : time_scaling_factor
  use ovar, only : disable_cgns_output
  use ovar, only : thermofix
  use ovar, only : cgns_output_restart_time_ave_data
  use thermofix_mod, only : check_thermo_properties
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iflag
  !
  !.. Local Scalars ..
  integer :: n,ierr,final_iter_digits
  !
  character(len=CGLEN) :: tname
  character(len=100)   :: char_num
  character(len=170)   :: sol_file
  !
  !.. Local Saved Variables ..
  logical(lk), save :: this_is_initial_output = true
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_parallel_cgns"
  !
continue
  !
  if (completely_disable_cgns) return
  !
  if (disable_cgns_output) then
    if (iflag == add_to_CGNS) then
      if (.not. write_TaylorGreen_full_solution) then
        return
      end if
    end if
  end if
  !
  call debug_timer(entering_procedure,pname)
  !
  if (this_is_initial_output) then
    !
    ! Create the cell_k0 array if running in debug mode and using OpenMP
    ! which will allow using OpenMP to speed up the interpolations
    ! computed in interpolate_for_output.
    !
    cell_k0 = create_cell_output_offsets() ! F2003 AUTO-REALLOC
    !
    ! Create the format string for creating the CGNS file names
    !
    write (char_num,4) itrst + num_timesteps
    final_iter_digits = max(4,len_trim(adjustl(char_num)))
    !
    write (cgnsfile_fmt,3) final_iter_digits
    !
    ! Create the CGNS grid file
    !
    call write_parallel_cgns_grid_file
    !
    ! Initialize the zone counter
    !
    cgns_zone_counter = 0
    !
    this_is_initial_output = fals
    !
  end if
  !
  ! Get the name for the current solution node
  !
  write (tname,'("FlowSolution",i0)') cgns_zone_counter
 !write (sol_file,'("solution.",i0,".cgns")') cgns_zone_counter
  if (iflag == -3) then
    write (sol_file,cgnsfile_fmt) max(0,itcur-1),".NaN"
  else
    write (sol_file,cgnsfile_fmt) itcur,""
  end if
  !
  ! Reallocate the soltime and solname arrays so we
  ! can concatenate this solution onto those arrays
  !
  call reallocate(soltime,cgns_zone_counter+1)
  call reallocate(soliter,cgns_zone_counter+1)
  call reallocate(solname,cgns_zone_counter+1)
  call reallocate(solfile,cgns_zone_counter+1)
  soltime(size(soltime)) = real(time*time_ref*time_scaling_factor,kind(soltime))
  soliter(size(soliter)) = itcur
  solname(size(solname)) = tname
  solfile(size(solfile)) = sol_file
  !
  ! Check/fix bad thermo properties
  !
 !call check_thermo_properties( thermofix%check_cells, thermofix%check_faces, &
 !                              thermofix%check_nodes, update_faces=true, &
 !                              verbose=thermofix%verbose )
  !
  ! Write the solution file
  !
  call write_parallel_cgns_solution_file(sol_file,tname)
  !
  ! Write the solution file
  !
  if (output_time_averaging .and. iflag /= -3) then
    call write_parallel_cgns_time_ave_file
    if (cgns_output_restart_time_ave_data) then
      call write_parallel_cgns_timeave_data_file
    end if
  end if
  !
  ! ################################################
  ! #####   FINISHED WRITING FIELD VARIABLES   #####
  ! ################################################
  !
  ! Write the time iterative solution file
  !
  if (any(iflag == [close_CGNS,-3])) then
    call write_parallel_cgns_time_file
  end if
  !
  ! Update the zone counter
  !
  cgns_zone_counter = cgns_zone_counter + 1
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("#############################################",/, &
            "    FOR PROCESSOR #",i0)
  2 format (a," = ",i0)
  3 format ("('solution.',i0.",i0,",a,'.cgns')")
  4 format (i0)
  !
end subroutine write_parallel_cgns
!
!###############################################################################
!
subroutine write_parallel_cgns_time_file()
  !
  !.. Use Statements ..
  use ovar, only : profile_io_cgns,profile_io_cgns_subs
  use ovar, only : cgns_dir
  !
  !.. Local Scalars ..
  integer :: n
  !
  integer(CBT) :: ifile_t,ibase_t,izone_t,iiter_t
  !
  real(wp) :: beg_wall_time
  real(wp) :: end_wall_time
  real(wp) :: cgns_wall_time
  !
  character(len=300) :: link_path
  !
  !.. Local Arrays ..
  integer(CST) :: idata(1:2)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_parallel_cgns_time_file"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
#ifndef SPECIAL_FOR_PGI
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      write (iout,1) pname
      flush (iout)
      beg_wall_time = mpi_wtime()
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  iiter_t = 1_CBT
  !
  call cgp_mpi_comm_f(mpi_comm_self_val,cgierr)
  !
  if (mypnum == glb_root) then
    !
    cgns_time_file = trim(adjustl(cgns_dir)) // cgns_time_fname
    !
    ! Start the timer for writing this time file that links all
    ! the previous CGNS solution files into a single CGNS file.
    !
    call cgns_timer(start_cgns_file_timer,profile_io_cgns_subs, &
                    cgns_time_file,override_barrier=true)
    !
    ! Open the time file linking all the previous CGNS solution files into
    ! a single CGNS file.
    ! NOTE: This opens the CGNS file and creates the base and zone sections.
    !
    call open_cgns_file(cgns_time_file,pname,ifile_t,ibase_t,izone_t, &
                        not_parallel=true)
    !
    ! Create the grid coordinates and element connectivity sections by
    ! using links to these respective sections within the CGNS grid file.
    !
    call link_cgns_grid_file(cgns_time_file,pname,ifile_t,ibase_t,izone_t, &
                             not_parallel=true)
    !
    ! Create a links to the solution node in each of the solution files
    !
   !call cg_goto_f(ifile_t,ibase_t,cgierr, &
   !               "Zone_t",izone_t, &
   !               CGNS_GOTO_END)
   !call cgns_error(pname,cgns_time_file,"cg_goto_f", &
   !                cgierr,__LINE__,__FILE__,not_parallel=true)
    !
    do n = 1,size(solname)
      !
      link_path = base_name//"/"//zone_name//"/"//trim(solname(n))
      !
      call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                      "cg_link_write_f",override_barrier=true)
      call cg_link_write_f(trim(solname(n)),trim(solfile(n)), &
                           trim(link_path),cgierr)
      call cgns_error(pname,cgns_time_file,"cg_link_write_f", &
                      cgierr,__LINE__,__FILE__,not_parallel=true)
      call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                      "cg_link_write_f",override_barrier=true)
      !
    end do
    !
    call write_iterative_data_to_cgns_file(cgns_time_file,pname, &
                                           ifile_t,ibase_t,izone_t, &
                                           soltime,soliter,solname, &
                                           not_parallel=true)
    !
    ! Close the CGNS FILE
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_close_f",override_barrier=true)
    call cgp_close_f(ifile_t,cgierr)
    call cgns_error(pname,cgns_time_file,"cgp_close_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_close_f",override_barrier=true)
    !
    call cgns_timer(stop_cgns_file_timer,profile_io_cgns_subs, &
                    cgns_time_file,override_barrier=true)
    !
  end if
  !
  call cgp_mpi_comm_f(mpi_comm_world_val,cgierr)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      end_wall_time = mpi_wtime()
      cgns_wall_time = end_wall_time - beg_wall_time
      write (iout,2) pname,trim(adjustl(cgns_time_file)),cgns_wall_time
      flush (iout)
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," Entering ",a)
  2 format (1x,"Leaving ",a,/,:, &
            5x,"Wall time taken to write the file '",a, &
            "' = ",f12.6," seconds",/)
  !
end subroutine write_parallel_cgns_time_file
!
!###############################################################################
!
subroutine write_iterative_data_to_cgns_file(cgns_file,procedure_name, &
                                             ifile,ibase,izone, &
                                             sol_times,sol_iters,sol_names, &
                                             not_parallel)
  !
  !.. Use Statements ..
  use ovar, only : profile_io_cgns_subs
  !
  !.. Formal Arguments ..
  character(len=*),     intent(in) :: cgns_file
  character(len=*),     intent(in) :: procedure_name
  integer(CBT),         intent(in) :: ifile
  integer(CBT),         intent(in) :: ibase
  integer(CBT),         intent(in) :: izone
  real(sp),             intent(in) :: sol_times(:)
  integer,              intent(in) :: sol_iters(:)
  character(len=CGLEN), intent(in) :: sol_names(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: not_parallel
  !
  !.. Local Scalars ..
  integer(CBT) :: iiter_t
  !
  character(len=300) :: pname
  !
  !.. Local Arrays ..
  integer(CST) :: idata(1:2)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: sub_pname = "write_iterative_data_to_cgns_file"
  !
continue
  !
#ifndef SPECIAL_FOR_PGI
  !
  write (pname,1) sub_pname,trim(adjustl(procedure_name))
  !
  iiter_t = 1_CBT
  !
  !=============================================================================
  !
  ! Create the BaseIterativeData_t node
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_biter_write_f",override_barrier=not_parallel)
  call cg_biter_write_f(ifile,ibase,'TimeIterValues', &
                        size(sol_times,kind=CBT),cgierr)
  call cgns_error(pname,cgns_file,"cg_biter_write_f", &
                  cgierr,__LINE__,__FILE__,not_parallel=not_parallel)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_biter_write_f",override_barrier=not_parallel)
   !
   ! Go to the end of the BaseIterativeData_t node so we can write the
   ! solution time values for each zone.
   !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_goto_f",override_barrier=not_parallel)
  call cg_goto_f(ifile,ibase,cgierr, &
                 "BaseIterativeData_t",iiter_t, &
                 CGNS_GOTO_END)
  call cgns_error(pname,cgns_file,"cg_goto_f", &
                  cgierr,__LINE__,__FILE__,not_parallel=not_parallel)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_goto_f",override_barrier=not_parallel)
  !
  ! Write the solution time values for each zone 
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_array_write_f",override_barrier=not_parallel)
  call cg_array_write_f("TimeValues",CGNS_KIND_REAL,1_CBT, &
                        size(sol_times,kind=CST),sol_times,cgierr)
  call cgns_error(pname,cgns_file,"cg_array_write_f", &
                  cgierr,__LINE__,__FILE__,not_parallel=not_parallel)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_array_write_f",override_barrier=not_parallel)
  !
  ! Write the solution iteration numbers for each zone 
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_array_write_f",override_barrier=not_parallel)
  call cg_array_write_f("IterationValues",CGNS_KIND_INT,1_CBT, &
                        size(sol_iters,kind=CST),sol_iters,cgierr)
  call cgns_error(pname,cgns_file,"cg_array_write_f", &
                  cgierr,__LINE__,__FILE__,not_parallel=not_parallel)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_array_write_f",override_barrier=not_parallel)
   !
  !=============================================================================
  !
  ! Create the ZoneIterativeData_t node
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_ziter_write_f",override_barrier=not_parallel)
  call cg_ziter_write_f(ifile,ibase,izone,"ZoneIterativeData",cgierr)
  call cgns_error(pname,cgns_file,"cg_ziter_write_f", &
                  cgierr,__LINE__,__FILE__,not_parallel=not_parallel)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_ziter_write_f",override_barrier=not_parallel)
  !
  ! Go to the end of the ZoneIterativeData_t node so we can write the
  ! FlowSolutionPointers data to link each zone to their respective
  ! data in the BaseIterativeData_t node.
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_goto_f",override_barrier=not_parallel)
  call cg_goto_f(ifile,ibase,cgierr, &
                 "Zone_t",izone, &
                 "ZoneIterativeData_t",iiter_t, &
                 CGNS_GOTO_END)
  call cgns_error(pname,cgns_file,"cg_goto_f", &
                  cgierr,__LINE__,__FILE__,not_parallel=not_parallel)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_goto_f",override_barrier=not_parallel)
  !
  ! Write the FlowSolutionPointers data to the ZoneIterativeData_t node.
  !
  idata(:) = int( [32,size(sol_names)] , kind=kind(idata) )
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_array_write_f",override_barrier=not_parallel)
  call cg_array_write_f("FlowSolutionPointers",CHARACTER,2_CBT, &
                        idata,sol_names,cgierr)
  call cgns_error(pname,cgns_file,"cg_array_write_f", &
                  cgierr,__LINE__,__FILE__,not_parallel=not_parallel)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_array_write_f",override_barrier=not_parallel)
  !
  !=============================================================================
  !
  ! Identify the flow solution within the CGNS file as time dependent
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_simulation_type_write_f",override_barrier=not_parallel)
  call cg_simulation_type_write_f(ifile,ibase,TIMEACCURATE,cgierr)
  call cgns_error(pname,cgns_file,"cg_simulation_type_write_f", &
                  cgierr,__LINE__,__FILE__,not_parallel=not_parallel)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_simulation_type_write_f",override_barrier=not_parallel)
  !
#endif
  !
  ! Format Statements
  !
  1 format (a," (called from ",a,")")
  !
end subroutine write_iterative_data_to_cgns_file
!
!###############################################################################
!
subroutine write_parallel_cgns_grid_file()
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr,n_solpts,xyz
  use ovar,        only : profile_io_cgns,profile_io_cgns_subs
  use ovar,        only : grid_scaling_factor
  use ovar,        only : itestcase,channel_nondim_output_coord
  use ovar,        only : cgns_dir
  use ovar,        only : output_cylindrical_coordinates
  !
  !.. Local Scalars ..
  integer :: l,n,ierr
  !
  integer(CBT) :: ifile_g,ibase_g,izone_g,isol_g,igrid_g
  !
  integer(CST) :: npts
  !
  integer(CET) :: cgpio_mode
  !
  real(wp) :: beg_wall_time
  real(wp) :: end_wall_time
  real(wp) :: cgns_wall_time
  real(wp) :: coord_output_scaling,theta
  !
  character(len=CGLEN) :: varname
  !
  !.. Local Arrays ..
  integer(CBT) :: ixyz_g(1:nr)
  integer(CBT) :: iarray_g(1:nr)
  !
  !.. Local Allocatable Arrays ..
  real(CRT), allocatable :: var(:)
  real(wp),  allocatable :: vartmp(:)
  real(wp),  allocatable :: xyztmp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_parallel_cgns_grid_file"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
#ifndef SPECIAL_FOR_PGI
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      write (iout,1) pname
      flush (iout)
      beg_wall_time = mpi_wtime()
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  cgns_grid_file = trim(adjustl(cgns_dir)) // cgns_grid_fname
  !
  ! Get the element connectivities and point/node and cell beginning and
  ! ending indices that are needed for CGNS output. This will basically
  ! create the cgns_sect(:) derived type, but it also gets the values for
  ! the module variables/arrays my_pts_beg, my_pts_end, and zone_size.
  !
  call get_connectivity_for_cgns()
  !
  ! Start the timer for writing the grid file
  !
  call cgns_timer(start_cgns_file_timer,profile_io_cgns_subs,cgns_grid_file)
  !
  ! Open the CGNS grid file
  ! NOTE: This opens the CGNS files and creates the base and zone sections
  ! NOTE: If using a mixed element grid and cgns_mixed_output is not being
  !       used (i.e., size(cgns_sect) > 1), we need to open the CGNS file
  !       for process-independent parallel IO instead of the default
  !       collective parallel IO. This is needed in case any of the MPI
  !       processes has zero cells matching the element type of one of the
  !       sections.
  !
  cgpio_mode = CGP_COLLECTIVE
  if (size(cgns_sect) > 1) then
    cgpio_mode = CGP_INDEPENDENT
  end if
  call open_cgns_file(cgns_grid_file,pname,ifile_g,ibase_g,izone_g, &
                      cgpio_mode=cgpio_mode)
  !
  ! Create data nodes for the coordinates
  !
  !
  ! Create a new GridCoordinates_t node for the Cartesian coordinates
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_grid_write_f")
  call cg_grid_write_f(ifile_g,ibase_g,izone_g, &
                       "GridCoordinates",igrid_g,cgierr)
  call cgns_error(pname,cgns_grid_file,"cg_grid_write_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_grid_write_f")
  !
  ! Move to the end of the CGNS file to write the Cartesian coordinate
  ! data using the cgp_array_write functions to this 
  !
  call cg_goto_f(ifile_g,ibase_g,cgierr, &
                 "Zone_t",izone_g, &
                 "GridCoordinates_t",igrid_g, &
                 CGNS_GOTO_END)
  call cgns_error(pname,cgns_grid_file,"cg_goto_f",cgierr,__LINE__,__FILE__)
  !
  do l = 1,nr
    varname = "Coordinate"//char_xyz(l)
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_coord_write_f")
    call cgp_coord_write_f(ifile_g,ibase_g,izone_g,CGNS_KIND_REAL, &
                           trim(varname),ixyz_g(l),cgierr)
    call cgns_error(pname,cgns_grid_file,"cgp_coord_write_f", &
                    cgierr,__LINE__,__FILE__,"Write node: Coordinate",l)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_coord_write_f")
  end do
  !
  ! Now create and write the coordinate data in parallel
  !
  ! Allocate the temporary variable arrays
  !
  npts = my_pts_end - my_pts_beg + 1_CST
  !
  xyztmp = init_array(npts,int(nr,kind=CST),zero)  ! F2003 AUTO-REALLOC
  vartmp = init_array(npts,zero)     ! F2003 AUTO-REALLOC
  var = init_array(npts,real(0,CRT)) ! F2003 AUTO-REALLOC
  !
  ! Interpolate the coordinates to the output locations and
  ! store in the temporary array xyztmp
  !
  do l = 1,nr
    !
    ! Interpolate this coordinate to the output locations
    !
    call interpolate_for_output( xyz(l,1:n_solpts) , xyztmp(1:npts,l) )
    !
  end do
  !
  ! Get the scaling constant used to remove the grid scaling factor or
  ! apply any other nondimensionalization to the coordinates
  !
  coord_output_scaling = one/grid_scaling_factor
  !
  if (itestcase == abs(Channel_Flow) .and. channel_nondim_output_coord) then
    ! Nondimensionalize the grid coordinates by the channel-half-height
    coord_output_scaling = coord_output_scaling / (half*channel_grid_Ly)
  end if
  !
  ! Apply the output coordinate scaling to the output coordinates
  !
  do l = 1,nr
    do n = 1,npts
      xyztmp(n,l) = xyztmp(n,l) * coord_output_scaling
    end do
  end do
  !
  ! Write the output coordinate data to the CGNS file
  !
  do l = 1,nr
    !
    ! Convert this output coordinate to the real type required for the CGNS file
    !
    var(:) = real( xyztmp(1:npts,l) , kind(var) )
    !
    ! Write the current coordinate variable to the CGNS file
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_coord_write_data_f")
    call cgp_coord_write_data_f(ifile_g,ibase_g,izone_g,ixyz_g(l), &
                                my_pts_beg,my_pts_end,var,cgierr)
    call cgns_error(pname,cgns_grid_file,"cgp_coord_write_data_f", &
                    cgierr,__LINE__,__FILE__,"Write data: Coordinate",l)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_coord_write_data_f")
    !
  end do
  !
  ! Write the cylindrical coordinates as flow variables to the grid file
  !
  if (output_cylindrical_coordinates .and. nr == 3) then
    !
    ! Create a new GridCoordinates_t node for the cylindrical coordinates
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_grid_write_f")
    call cg_grid_write_f(ifile_g,ibase_g,izone_g, &
                         "CylindricalCoordinates",igrid_g,cgierr)
    call cgns_error(pname,cgns_grid_file,"cg_grid_write_f", &
                    cgierr,__LINE__,__FILE__)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_grid_write_f")
    !
    ! Move to the end of the CGNS file to write the cylindrical coordinate
    ! data using the cgp_array_write functions to this GridCoordinates_t node.
    !
    call cg_goto_f(ifile_g,ibase_g,cgierr, &
                   "Zone_t",izone_g, &
                   "GridCoordinates_t",igrid_g, &
                   CGNS_GOTO_END)
    call cgns_error(pname,cgns_grid_file,"cg_goto_f",cgierr,__LINE__,__FILE__)
    !
    do l = 1,nr
      varname = "Coordinate"//char_cyl(l)
      call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                      "cgp_array_write_f")
      call cgp_array_write_f(trim(varname),CGNS_KIND_REAL,1_CBT,zone_size, &
                             iarray_g(l),cgierr)
      call cgns_error(pname,cgns_grid_file,"cgp_array_write_f", &
                      cgierr,__LINE__,__FILE__,"Write node: Coordinate Array",l)
      call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                      "cgp_array_write_f")
    end do
    !
    ! Compute the axial coordinate and write it to the CGNS grid file
    !
    do n = 1,npts
      vartmp(n) = xyztmp(n,NCD_Axial)
    end do
    var(:) = real( vartmp , kind=kind(var) )
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_array_write_data_f")
    call cgp_array_write_data_f(iarray_g(1),my_pts_beg,my_pts_end,var,cgierr)
    call cgns_error(pname,cgns_grid_file,"cgp_array_write_data_f", &
                    cgierr,__LINE__,__FILE__,"Write data: Coordinate Array",1)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_array_write_data_f")
    !
    ! Compute the radial coordinate and write it to the CGNS grid file
    !
    do n = 1,npts
      vartmp(n) = sqrt(xyztmp(n,NCD_Vert)**2 + xyztmp(n,NCD_Span)**2)
    end do
    var(:) = real( vartmp , kind=kind(var) )
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_array_write_data_f")
    call cgp_array_write_data_f(iarray_g(2),my_pts_beg,my_pts_end,var,cgierr)
    call cgns_error(pname,cgns_grid_file,"cgp_array_write_data_f", &
                    cgierr,__LINE__,__FILE__,"Write data: Coordinate Array",2)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_array_write_data_f")
    !
    ! Compute the theta coordinate and write it to the CGNS grid file
    !
    do n = 1,npts
      theta = atan2(xyztmp(n,NCD_Span),xyztmp(n,NCD_Vert))
      if (theta < zero) theta = theta + two*pi
      vartmp(n) = theta*radians2degrees
    end do
    var(:) = real( vartmp , kind=kind(var) )
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_array_write_data_f")
    call cgp_array_write_data_f(iarray_g(3),my_pts_beg,my_pts_end,var,cgierr)
    call cgns_error(pname,cgns_grid_file,"cgp_array_write_data_f", &
                    cgierr,__LINE__,__FILE__,"Write data: Coordinate Array",3)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_array_write_data_f")
    !
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  ! Deallocate the var and vartmp arrays
  !
  if (allocated(var)) then
    deallocate ( var , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"var",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(vartmp)) then
    deallocate ( vartmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vartmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(xyztmp)) then
    deallocate ( xyztmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyztmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Write all of the CGNS sections which contain the element connectivities
  !
  do n = 1,size(cgns_sect)
    call write_cgns_section(ifile_g,ibase_g,izone_g,cgns_sect(n),sect_idx=n)
  end do
  !
  ! Close the CGNS FILE
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cgp_close_f")
  call cgp_close_f(ifile_g,cgierr)
  call cgns_error(pname,cgns_grid_file,"cgp_close_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cgp_close_f")
  !
  call cgns_timer(stop_cgns_file_timer,profile_io_cgns_subs,cgns_grid_file)
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      end_wall_time = mpi_wtime()
      cgns_wall_time = end_wall_time - beg_wall_time
      write (iout,2) pname,trim(adjustl(cgns_grid_file)),cgns_wall_time
      flush (iout)
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," Entering ",a)
  2 format (1x,"Leaving ",a,/,:, &
            5x,"Wall time taken to write the file '",a, &
            "' = ",f12.6," seconds",/)
  !
end subroutine write_parallel_cgns_grid_file
!
!###############################################################################
!
subroutine open_cgns_file(cgns_file,procedure_name,ifile,ibase,izone, &
                          cgpio_mode,not_parallel)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  use ovar,        only : profile_io_cgns_subs
  use ovar,        only : lustre_stripe_factor
  use ovar,        only : lustre_stripe_unit
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: cgns_file
  character(len=*), intent(in) :: procedure_name
  integer(CBT),  intent(inout) :: ifile
  integer(CBT),  intent(inout) :: ibase
  integer(CBT),  intent(inout) :: izone
  !
  !.. Optional Arguments ..
  integer(CET), optional, intent(in) :: cgpio_mode
  logical(lk),  optional, intent(in) :: not_parallel
  !
  !.. Local Scalars ..
  integer(CET) :: cgns_pio_mode
  character(len=300) :: pname
  !
  _MPI_INFO_TYPE_ :: file_info
  !
  !.. Local Parameters ..
  character(len=*), parameter :: sub_pname = "open_cgns_file"
  !
continue
  !
  write (pname,1) sub_pname,trim(adjustl(procedure_name))
  !
  cgns_pio_mode = CGP_COLLECTIVE
  if (present(cgpio_mode)) then
    if (cgpio_mode == CGP_INDEPENDENT) then
      ! Make sure to only change this if it takes the other valid value
      ! of CGP_INDEPENDENT, otherwise keep it set to the default value.
      cgns_pio_mode = CGP_INDEPENDENT
    end if
  end if
  !
#ifdef PBS_ENV
  !
  ! Set the CGNS parallel IO mode to independent for PBS
  !
  cgns_pio_mode = CGP_INDEPENDENT
  !
  ! Create an MPI Info object for writing the CGNS file in the PBS environment
  call mpi_info_create(file_info,mpierr)
  ! Tell MPI this is a write once operation
  call mpi_info_set(file_info,"access_style","write_once",mpierr)
  ! Disable writing with collective buffering
  call mpi_info_set(file_info,"romio_cb_write","disable",mpierr)
  ! Add the striping factor to the MPI Info object
  call mpi_info_set(file_info,"striping_factor", &
                    trim(adjustl(lustre_stripe_factor)),mpierr)
  ! Add the striping unit to the MPI Info object
  call mpi_info_set(file_info,"striping_unit", &
                    trim(adjustl(lustre_stripe_unit)),mpierr)
  !
  ! Send the MPI Info object to the CGNS library before opening the file
  !
#ifdef USE_MPI_F08
  call cgp_mpi_info_f(file_info%mpi_val,cgierr)
#else
  call cgp_mpi_info_f(file_info,cgierr)
#endif
#endif
  !
  call cgp_pio_mode_f(cgns_pio_mode,cgierr)
  call cgns_error(pname,cgns_file,"cgp_pio_mode_f", &
                  cgierr,__LINE__,__FILE__,not_parallel=not_parallel)
  !
  ! Open the CGNS solution file
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cgp_open_f",override_barrier=not_parallel)
  call cgp_open_f(cgns_file,CG_MODE_WRITE,ifile,cgierr)
  call cgns_error(pname,cgns_file,"cgp_open_f", &
                  cgierr,__LINE__,__FILE__, &
                  "Mode = CG_MODE_WRITE",not_parallel=not_parallel)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cgp_open_f",override_barrier=not_parallel)
  !
  ! Create the CGNS base for this solution
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_base_write_f",override_barrier=not_parallel)
  call cg_base_write_f(ifile,base_name,int(nr,kind=CBT), &
                       int(nr,kind=CBT),ibase,cgierr)
  call cgns_error(pname,cgns_file,"cg_base_write_f", &
                  cgierr,__LINE__,__FILE__,not_parallel=not_parallel)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_base_write_f",override_barrier=not_parallel)
  !
  ! Write the zone information for this solution
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_zone_write_f",override_barrier=not_parallel)
  call cg_zone_write_f(ifile,ibase,zone_name,zone_size, &
                       UNSTRUCTURED,izone,cgierr)
  call cgns_error(pname,cgns_file,"cg_zone_write_f", &
                  cgierr,__LINE__,__FILE__,not_parallel=not_parallel)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_zone_write_f",override_barrier=not_parallel)
  !
  ! Format Statements
  !
  1 format (a," (called from ",a,")")
  !
end subroutine open_cgns_file
!
!###############################################################################
!
subroutine link_cgns_grid_file(cgns_file,procedure_name, &
                               ifile,ibase,izone,not_parallel)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  use ovar,        only : profile_io_cgns_subs
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: cgns_file
  character(len=*), intent(in) :: procedure_name
  integer(CBT),     intent(in) :: ifile
  integer(CBT),     intent(in) :: ibase
  integer(CBT),     intent(in) :: izone
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: not_parallel
  !
  !.. Local Scalars ..
  integer :: n
  !
  character(len=CGLEN) :: sect_name
  !
  character(len=300) :: link_path
  character(len=300) :: pname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: sub_pname = "link_cgns_grid_file"
  !
continue
  !
  write (pname,1) sub_pname,trim(adjustl(procedure_name))
  !
  ! Create a link to the grid coordinates in the grid file
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_goto_f",override_barrier=not_parallel)
  call cg_goto_f(ifile,ibase,cgierr, &
                 "Zone_t",izone, &
                 CGNS_GOTO_END)
  call cgns_error(pname,cgns_file,"cg_goto_f", &
                  cgierr,__LINE__,__FILE__,not_parallel=not_parallel)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_goto_f",override_barrier=not_parallel)
  !
  link_path = base_name//"/"//zone_name//"/GridCoordinates"
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f",override_barrier=not_parallel)
  call cg_link_write_f("GridCoordinates",trim(cgns_grid_fname), &
                       trim(link_path),cgierr)
  call cgns_error(pname,cgns_file,"cg_link_write_f(GridCoordinates)", &
                  cgierr,__LINE__,__FILE__,not_parallel=not_parallel)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f",override_barrier=not_parallel)
 !!
 !! Create a link to the element connectivity in the grid file
 !!
 !link_path = base_name//"/"//zone_name//"/"//trim(cgns_sect_name)
 !!
 !call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
 !                "cg_link_write_f",override_barrier=not_parallel)
 !call cg_link_write_f(trim(cgns_sect_name),trim(cgns_grid_fname), &
 !                     trim(link_path),cgierr)
 !call cgns_error(pname,cgns_file, &
 !                "cg_link_write_f("//trim(cgns_sect_name)//")", &
 !                cgierr,__LINE__,__FILE__,not_parallel=not_parallel)
 !call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
 !                "cg_link_write_f",override_barrier=not_parallel)
  !
  ! Create a link to the element connectivity in the grid file
  !
  do n = 1,size(cgns_sect)
    !
    sect_name = trim(adjustl(cgns_sect(n)%name))
    !
    link_path = base_name//"/"//zone_name//"/"//trim(sect_name)
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_link_write_f",override_barrier=not_parallel)
    call cg_link_write_f(trim(sect_name),trim(cgns_grid_fname), &
                         trim(link_path),cgierr)
    call cgns_error(pname,cgns_file, &
                    "cg_link_write_f("//trim(sect_name)//")", &
                    cgierr,__LINE__,__FILE__,not_parallel=not_parallel)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_link_write_f",override_barrier=not_parallel)
    !
  end do
  !
  ! Format Statements
  !
  1 format (a," (called from ",a,")")
  !
end subroutine link_cgns_grid_file
!
!###############################################################################
!
subroutine link_cylindrical_coordinates(cgns_file,procedure_name, &
                                        ifile,ibase,izone,isol)
  !
  !.. Use Statements ..
  use ovar, only : profile_io_cgns_subs
  use ovar, only : output_cylindrical_coordinates
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: cgns_file
  character(len=*), intent(in) :: procedure_name
  integer(CBT),     intent(in) :: ifile
  integer(CBT),     intent(in) :: ibase
  integer(CBT),     intent(in) :: izone
  integer(CBT),     intent(in) :: isol
  !
  !.. Local Scalars ..
  character(len=300) :: link_path
  character(len=300) :: pname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: sub_pname = "link_cylindrical_coordinates"
  !
continue
  !
  ! Return without doing anything if the cylindrical coordinates
  ! were not written to the CGNS grid file.
  !
  if (.not. output_cylindrical_coordinates) then
    return
  end if
  !
  !
  !
  write (pname,1) sub_pname,trim(adjustl(procedure_name))
  !
  ! Create a link to the axial grid coordinates in the grid file
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cg_goto_f")
  call cg_goto_f(ifile,ibase,cgierr, &
                 "Zone_t",izone, &
                 "FlowSolution_t",isol, &
                 CGNS_GOTO_END)
  call cgns_error(pname,cgns_file,"cg_goto_f",cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cg_goto_f")
  !
  link_path = base_name//"/"//zone_name// &
              "/CylindricalCoordinates/CoordinateAxial"
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  call cg_link_write_f("CoordinateAxial",trim(cgns_grid_fname), &
                       trim(link_path),cgierr)
  call cgns_error(pname,cgns_file,"cg_link_write_f(CoordinateAxial)", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  !
  ! Create a link to the radial grid coordinates in the grid file
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cg_goto_f")
  call cg_goto_f(ifile,ibase,cgierr, &
                 "Zone_t",izone, &
                 "FlowSolution_t",isol, &
                 CGNS_GOTO_END)
  call cgns_error(pname,cgns_file,"cg_goto_f",cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cg_goto_f")
  !
  link_path = base_name//"/"//zone_name// &
              "/CylindricalCoordinates/CoordinateR"
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  call cg_link_write_f("CoordinateR",trim(cgns_grid_fname), &
                       trim(link_path),cgierr)
  call cgns_error(pname,cgns_file,"cg_link_write_f(CoordinateR)", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  !
  ! Create a link to the theta grid coordinates in the grid file
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cg_goto_f")
  call cg_goto_f(ifile,ibase,cgierr, &
                 "Zone_t",izone, &
                 "FlowSolution_t",isol, &
                 CGNS_GOTO_END)
  call cgns_error(pname,cgns_file,"cg_goto_f",cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cg_goto_f")
  !
  link_path = base_name//"/"//zone_name// &
              "/CylindricalCoordinates/CoordinateTheta"
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  call cg_link_write_f("CoordinateTheta",trim(cgns_grid_fname), &
                       trim(link_path),cgierr)
  call cgns_error(pname,cgns_file,"cg_link_write_f(CoordinateTheta)", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  !
  ! Format Statements
  !
  1 format (a," (called from ",a,")")
  !
end subroutine link_cylindrical_coordinates
!
!###############################################################################
!
subroutine write_cgns_section(ifil,ibas,izon,s,sect_idx)
  !
  !.. Use Statements ..
  use iso_c_binding, only : c_null_ptr
  use ovar,          only : profile_io_cgns_subs
  !
  !.. Formal Arguments ..
  integer(CBT),         intent(in) :: ifil
  integer(CBT),         intent(in) :: ibas
  integer(CBT),         intent(in) :: izon
  type(cgns_sect_t), intent(inout) :: s
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: sect_idx
  !
  !.. Local Scalars ..
  integer(CBT) :: isec
  integer(CET) :: ityp
  integer(CST) :: glb1,glb2
  integer(CST) :: loc1,loc2
  integer(CST) :: mxoff
  !
  character(len=CGLEN) :: snam
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_cgns_section"
  !
continue
  !
  snam = trim(adjustl(s%name))
  !
  ityp = s%type
  !
  glb1 = s%all_cells_beg
  glb2 = s%all_cells_end
  !
  loc1 = s%my_cells_beg
  loc2 = s%my_cells_end
  !
  mxoff = s%MaxOffset
  !
  if (ityp == MIXED) then
    !
    ! For mixed sections, we need to use the cgp_poly_... functions
    !
    ! First, create the empty section to prepare for writing the element
    ! connectivities for this section in parallel
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_poly_section_write_f")
    call cgp_poly_section_write_f(ifil,ibas,izon,snam,ityp,glb1,glb2,mxoff, &
                                  0_CBT,isec,cgierr)
    call cgns_error(pname,cgns_grid_file,"cgp_poly_section_write_f", &
                    cgierr,__LINE__,__FILE__)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_poly_section_write_f")
    !
    ! Now write the actual element connectivities in parallel
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_poly_elements_write_data_f")
    call cgp_poly_elements_write_data_f(ifil,ibas,izon,isec,loc1,loc2, &
                                        s%ElemConn,s%Offsets,cgierr)
    call cgns_error(pname,cgns_grid_file,"cgp_poly_elements_write_data_f", &
                    cgierr,__LINE__,__FILE__)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_poly_elements_write_data_f")
    !
  else
    !
    ! For sections containing a single element type, just use the normal
    ! cgp_section_write_f and cgp_element_write_data_f functions
    !
    ! First, create the empty section to prepare for writing the element
    ! connectivities for this section in parallel
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_section_write_f")
    call cgp_section_write_f(ifil,ibas,izon,snam,ityp,glb1,glb2, &
                             0_CBT,isec,cgierr)
    call cgns_error(pname,cgns_grid_file,"cgp_section_write_f", &
                    cgierr,__LINE__,__FILE__)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_section_write_f")
    !
    ! Now write the actual element connectivities in parallel
    ! NOTE: Only have those paritions with cells matching the element type
    !       for this section (i.e., s%ElemConn has been allocated) write
    !       their element connectivity for this section to the CGNS grid
    !       file.
    !     - I tried passing c_null_ptr instead of s%ElemConn if it wasnt
    !       allocated, but this just overwrote the first element of the
    !       section with an element connectivity of all 0's. 
    !     - I also tried making loc1 and loc2 zero, but that just
    !       returned errors from CGNS. 
    !     - Sending in the unallocated s%ElemConn array itself results in
    !       runtime errors for trying to access an unallocated array.
    ! NOTE: Since not all MPI processes are guaranteed to call
    !       cgp_elements_write_data_f, we need to make sure that the
    !       parallel IO mode is set to independent instead of the default
    !       collective. This is done in the call to open_cgns_file.
    !
#ifdef DEBUG_ON
    if (ncpu == 1) then
      write (iout,3) "snam",trim(adjustl(snam))
      write (iout,1) "ifil",ifil
      write (iout,1) "ibas",ibas
      write (iout,1) "izon",izon
      write (iout,1) "n_total_pts",zone_size(1)
      write (iout,1) "n_total_cells",zone_size(2)
      write (iout,1) "isec",isec
      write (iout,1) "glb1",glb1
      write (iout,1) "glb2",glb2
      write (iout,1) "loc1",loc1
      write (iout,1) "loc2",loc2
      write (iout,1) "ityp",ityp
      write (iout,1) "mxoff",mxoff
      write (iout,2) "allocated(s%ElemConn)",allocated(s%ElemConn)
      if (allocated(s%ElemConn)) then
        write (iout,1) "size(s%ElemConn)",size(s%ElemConn)
        write (iout,1) "lbound(s%ElemConn)",lbound(s%ElemConn,dim=1)
        write (iout,1) "ubound(s%ElemConn)",ubound(s%ElemConn,dim=1)
      end if
      flush (iout)
      1 format (2x,a," = ",i0)
      2 format (2x,a," = ",l1)
      3 format (2x,a," = ",a)
    end if
    flush (iout)
#endif
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_elements_write_data_f")
    if (allocated(s%ElemConn)) then
      call cgp_elements_write_data_f(ifil,ibas,izon,isec,loc1,loc2, &
                                     s%ElemConn,cgierr)
    end if
    call cgns_error(pname,cgns_grid_file,"cgp_elements_write_data_f", &
                    cgierr,__LINE__,__FILE__)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_elements_write_data_f")
    !
  end if
  !
  ! Save isec back into cgns_sect in case it will be needed later
  !
  s%isect_g = isec
  !
end subroutine write_cgns_section
!
!###############################################################################
!
subroutine write_parallel_cgns_solution_file(sol_file,tname)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  use ovar,        only : itcur
  use ovar,        only : tref,rgasref_nd,tmin
  use ovar,        only : profile_io_cgns,profile_io_cgns_subs
  use ovar,        only : max_wall_time => max_cgns_wall_time
  use ovar,        only : cgns_dir
  use ovar,        only : output_cylindrical_coordinates
  use flowvar_mod, only : usp,dusp,gpsp
  !
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  use flux_mod, only : discontinuous_gradient_at_solpts
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: sol_file
  character(len=*), intent(in) :: tname
  !
  !.. Local Scalars ..
  integer :: l,n,npts,ierr
  !
  real(wp) :: beg_wall_time
  real(wp) :: end_wall_time
  real(wp) :: cgns_wall_time
  !
  character(len=170) :: cgns_sol_file
  !
  !.. Local Allocatable Arrays ..
#ifdef SPECIAL_FOR_INTEL
  character(len=:), allocatable :: char_sname(:)
#else
  character(len=CGLEN), allocatable :: char_sname(:)
#endif
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_parallel_cgns_solution_file"
  !
  !.. Local Scalars ..
  integer(CBT) :: ifile_s,ibase_s,izone_s,isol_s
  !
  real(sp)             :: sol_time(1)
  integer              :: sol_iter(1)
  character(len=CGLEN) :: sol_name(1)
  !
  !.. Local Allocatable Arrays ..
  integer(CBT), allocatable, dimension(:) :: isvar_s
  !
  real(CRT), allocatable :: var(:)
  real(wp),  allocatable :: vartmp(:)
  real(wp),  allocatable :: interp_usp(:,:)
  real(wp),  allocatable :: interp_gpsp(:,:)
  real(wp),  allocatable :: interp_dusp(:,:,:)
  real(wp),  allocatable :: dusp_copy(:,:,:)
  !
 !logical(lk), parameter :: output_discontinuous_gradients = true
  logical(lk), parameter :: output_discontinuous_gradients = fals
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
#ifndef SPECIAL_FOR_PGI
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      write (iout,1) pname
      flush (iout)
      beg_wall_time = mpi_wtime()
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Get the list of variables that will be written to the CGNS solution file
  !
  call get_cgns_sol_variables(char_sname)
  !
  ! Get the name of the current solution file
  !
  cgns_sol_file = trim(adjustl(cgns_dir)) // trim(adjustl(sol_file))
  !
  ! Start the timer for writing this solution file
  !
  call cgns_timer(start_cgns_file_timer,profile_io_cgns_subs,cgns_sol_file)
  !
  ! Open the solution file.
  ! NOTE: This opens the CGNS file and creates the base and zone sections.
  !
  call open_cgns_file(cgns_sol_file,pname,ifile_s,ibase_s,izone_s)
  !
  ! Create the grid coordinates and element connectivity sections by
  ! using links to these respective sections within the CGNS grid file.
  !
  call link_cgns_grid_file(cgns_sol_file,pname,ifile_s,ibase_s,izone_s)
  !
  ! Create the flow solution section of the CGNS file
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_sol_write_f")
  call cg_sol_write_f(ifile_s,ibase_s,izone_s,tname,VERTEX,isol_s,cgierr)
  call cgns_error(pname,cgns_sol_file,"cg_sol_write_f",cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_sol_write_f")
  !
  ! #######################################
  ! #####   WRITING FIELD VARIABLES   #####
  ! #######################################
  !
  npts = int( my_pts_end - my_pts_beg + 1_CST , kind=kind(npts) )
  !
  ! Allocate the local arrays used for output
  !
  var = init_array(npts,real(0,CRT))        ! F2003 AUTO-REALLOC
  vartmp = init_array(npts,zero)            ! F2003 AUTO-REALLOC
  interp_usp = init_array(npts,nq,zero)     ! F2003 AUTO-REALLOC
  interp_dusp = init_array(npts,nr,nq,zero) ! F2003 AUTO-REALLOC
  !
  ! Interpolate usp to the output points
  !
  call interpolate_for_output( usp , interp_usp )
  !
  if (allocated(dusp)) then
    !
    if (output_discontinuous_gradients) then
      !
      ! Calling discontinuous_gradient_at_solpts will replace the
      ! contents of dusp with the discontinuous gradient, so save
      ! the contents of dusp in dusp_copy.
      !
      dusp_copy = dusp(:,:,:) ! F2003 AUTOMATIC-REALLOCATION
      !
      call discontinuous_gradient_at_solpts()
      !
    end if
    !
    if (interpolated_dusp_is_needed) then
      call interpolate_for_output( dusp , interp_dusp )
    end if
    !
    if (output_discontinuous_gradients) then
      !
      ! Now that we have the discontinuous gradients in interp_dusp,
      ! copy the original contents of dusp that were saved in dusp_copy
      ! back into dusp, then deallocate dusp_copy.
      !
      dusp(:,:,:) = dusp_copy(:,:,:)
      !
      deallocate ( dusp_copy , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"dusp_copy",2,__LINE__,__FILE__,ierr,error_message)
      !
    end if
    !
  end if
  !
  ! ################################################
  ! #####   WRITING SOLUTION FIELD VARIABLES   #####
  ! ################################################
  !
  ! Count the total number of solution variables that will be written
  ! and allocate isvar_s accordingly
  !
  n = count_output_variables(char_sname)
  !
  ! Allocate the CGNS field indexing array
  !
  isvar_s = init_array(n,0_CBT) ! F2003 AUTO-REALLOC
  !
  ! Create data nodes for the solution variables
  !
  call write_field_headers(cgns_sol_file,ifile_s,ibase_s,izone_s,isol_s, &
                           char_sname,isvar_s)
  !
  ! Allocate the array for storing the gas properties
  ! evaluated at the output points
  !
  interp_gpsp = init_array(npts,ng,zero) ! F2003 AUTO-REALLOC
  !
  ! The most accurate gas properties (GP) we can output are those computed
  ! from the conserved variable (CV) values which have already been
  ! interpolated to the output points (OP). This is because the CV are an
  ! actual polynomial good for interpolating, whereas the GP at the solution
  ! points (SP) are not a true polynomial because they are computed from the
  ! CV values at the SP. However, for the special case of thermally perfect
  ! gas mixture (TPG) simulations, computing the temperature from the CV
  ! requires a Newton solve and thus an initial guess for the temperature.
  ! To get this initial guess, we need to interpolate the GP in gpsp to the
  ! OP and store in interp_gpsp. These interpolated GP in interp_gpsp will be
  ! overwritten once update_thermo_properties is called. NOTE, this only
  ! needs to be done for simulations using a TPG, since the GP for a
  ! calorically perfect gas are simply computed from the CV.
  !
  if (ncs > 0) then
    call interpolate_for_output( gpsp , interp_gpsp )
  end if
  !
  ! Now actually compute the updated thermodynamic and transport properties
  ! NOTE: If this is a TPG simulation, interp_gpsp has been initialized to
  !       the interpolated values of gpsp so the Newton solve has an initial
  !       guess for computing temperature.
  !
  call update_thermo_properties(1,npts,interp_usp,interp_gpsp, &
                                caller_name=pname)
  if (itcur == 0) then
    if (any(interp_gpsp(:,ntemp) <= zero) .or. &
        any(interp_usp(:,nec) <= zero)) then
      do n = 1,size(interp_gpsp,dim=1)
        if (interp_gpsp(n,ntemp) <= tmin/tref) then
          interp_gpsp(n,ntemp) = tmin/tref
          if (interp_usp(n,nec) < zero) then
            interp_usp(n,nec) = interp_gpsp(n,npres) / &
                                (rgasref_nd*interp_gpsp(n,ntemp))
          else
            interp_gpsp(n,npres) = interp_usp(n,nec) * &
                                   rgasref_nd*interp_gpsp(n,ntemp)
          end if
        end if
      end do
    end if
  end if
  call update_transport_properties(1,npts,interp_usp,interp_gpsp)
  !
  ! Now create the output data and write it to the CGNS file in parallel
  !
  n = 0
  do l = 1,size(char_sname)
    !
    call output_variable(ifile_s,ibase_s,izone_s,isol_s,isvar_s, &
                         cgns_sol_file,char_sname(l),n, &
                         interp_usp,interp_gpsp,interp_dusp,vartmp,var)
    !
  end do
  !
  ! Create a link to the grid coordinates in the grid file
  !
  if (output_cylindrical_coordinates .and. nr == 3) then
    call link_cylindrical_coordinates(cgns_sol_file,pname, &
                                      ifile_s,ibase_s,izone_s,isol_s)
  end if
  !
  ! ################################################
  ! #####   FINISHED WRITING FIELD VARIABLES   #####
  ! ################################################
  !
  ! Add the time information to the solution file
  !
  sol_time(1) = soltime(size(soltime))
  sol_name(1) = tname
 !sol_name(1) = solname(size(solname))
  sol_iter(1) = itcur
  call write_iterative_data_to_cgns_file(cgns_sol_file,pname, &
                                         ifile_s,ibase_s,izone_s, &
                                         sol_time,sol_iter,sol_name)
  !
  ! Close the CGNS FILE
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cgp_close_f")
  call cgp_close_f(ifile_s,cgierr)
  call io_error(pname,cgns_sol_file,2,__LINE__,__FILE__,cgierr)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cgp_close_f")
  !
  call cgns_timer(stop_cgns_file_timer,profile_io_cgns_subs,cgns_sol_file)
  !
  ! Deallocate the arrays var and vartmp
  !
  if (allocated(var)) then
    deallocate ( var , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"var",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(vartmp)) then
    deallocate ( vartmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vartmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(interp_usp)) then
    deallocate ( interp_usp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interp_usp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(interp_dusp)) then
    deallocate ( interp_dusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interp_dusp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(isvar_s)) then
    deallocate ( isvar_s , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"isvar_s",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      end_wall_time = mpi_wtime()
      cgns_wall_time = end_wall_time - beg_wall_time
      max_wall_time = max(max_wall_time,cgns_wall_time)
      write (iout,2) pname,trim(adjustl(cgns_sol_file)),cgns_wall_time
      flush (iout)
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," Entering ",a)
  2 format (1x,"Leaving ",a,/,:, &
            5x,"Wall time taken to write the file '",a, &
            "' = ",f12.6," seconds",/)
  !
end subroutine write_parallel_cgns_solution_file
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function count_output_variables(variable_names) result(return_value)
  !
  !.. Use Statements ..
  use geomvar_mod,     only : nr
  use ovar,            only : governing_equations
  use metrics_mod,     only : time_metrics
  use flowvar_mod,     only : nmarker
  !
  !.. Formal Arguments ..
#ifdef SPECIAL_FOR_INTEL
  character(len=:), allocatable, intent(in) :: variable_names(:)
#else
  character(len=CGLEN), intent(in) :: variable_names(:)
#endif
  !
  !.. Function Result ..
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: n,l
  !
continue
  !
  n = 0
  do l = 1,size(variable_names)
    if (trim(variable_names(l)) == "Velocity") then
      n = n + nr
    else if (trim(variable_names(l)) == "RelativeVelocity") then
      n = n + nr
    else if (trim(variable_names(l)) == "CylindricalVelocity") then
      n = n + 2
    else if (trim(variable_names(l)) == "RCylindricalVelocity") then
      n = n + 2
    else if (trim(variable_names(l)) == "Vorticity") then
      n = n + merge(nr,1,nr==3)
    else if (trim(variable_names(l)) == "Momentum") then
      n = n + nr
    else if (trim(variable_names(l)) == "RelativeMomentum") then
      n = n + nr
    else if (trim(variable_names(l)) == "GradDensity") then
      n = n + nr
    else if (trim(variable_names(l)) == "GradMomentum") then
      n = n + nr*nr
    else if (trim(variable_names(l)) == "GradVelocity") then
      n = n + nr*nr
    else if (trim(variable_names(l)) == "GradEnergy") then
      n = n + nr
    else if (trim(variable_names(l)) == "GradTemperature") then
      n = n + nr
    else if (trim(variable_names(l)) == "DmatTemperature") then
      n = n + nr
    else if (trim(variable_names(l)) == "GradMassFractions") then
      n = n + ncs*nr
    else if (trim(variable_names(l)) == "GradSpeciesDensities") then
      n = n + ncs*nr
    else if (trim(variable_names(l)) == "Viscosity") then
      if (governing_equations == NavierStokes_Eqns) then
        n = n + 1
      end if
    else if (trim(variable_names(l)) == "ResolutionMarker") then
      n = n + 2*nmarker
    else if (trim(variable_names(l)) == "JumpMarker") then
      n = n + 2*nmarker
    else if (trim(variable_names(l)) == "DucrosMarker") then
      n = n + 2*nmarker + 3
    else if (trim(variable_names(l)) == "MLPMarker") then
      n = n + 2*nmarker
    else if (trim(variable_names(l)) == "TimeMetrics") then
      if (allocated(time_metrics)) then
        n = n + nr + 1
      end if
    else if (trim(variable_names(l)) == "CoMetrics") then
   !  n = n + nr*nr
    else if (trim(variable_names(l)) == "ContraMetrics") then
      n = n + nr*nr
    else if (trim(variable_names(l)) == "MassFractions") then
      if (ncs > 0) then
        n = n + ncs
      end if
    else if (trim(variable_names(l)) == "MoleFractions") then
      if (ncs > 0) then
        n = n + ncs
      end if
    else if (trim(variable_names(l)) == "SpeciesDensities") then
      if (ncs > 0) then
        n = n + ncs
      end if
    else if (trim(variable_names(l)) == "SpeciesEnthalpy") then
      if (ncs > 0) then
        n = n + ncs
      end if
    else if (trim(variable_names(l)) == "SpeciesCp") then
      if (ncs > 0) then
        n = n + ncs
      end if
    else if (trim(variable_names(l)) == "cpi/R") then
      if (ncs > 0) then
        n = n + ncs
      end if
    else if (trim(variable_names(l)) == "hi/RT") then
      if (ncs > 0) then
        n = n + ncs
      end if
    else if (trim(variable_names(l)) == "hi/R") then
      if (ncs > 0) then
        n = n + ncs
      end if
    else if (trim(variable_names(l)) == "dTdx_Breakdown") then
      if (ncs > 0) then
        n = n + (2*nr + 2*ncs + 7 + ncs*nr)
      end if
    else if (trim(variable_names(l)) == "QCriterion") then
      n = n + 2
    else if (trim(variable_names(l)) == "LocalTimestep") then
      n = n + 2
    else if (trim(variable_names(l)) == "FilterDissipation") then
      n = n + nq + 1
    else if (trim(variable_names(l)) == "CellSize") then
      n = n + 2
    else
      n = n + 1
    end if
  end do
  !
  return_value = n
  !
end function count_output_variables
!
!###############################################################################
!
subroutine write_field_headers(cgns_sol_file,ifile,ibase,izone,isol, &
                               variable_names,ifld)
  !
  !.. Use Statements ..
  use geomvar_mod,     only : nr
  use ovar,            only : governing_equations
  use ovar,            only : profile_io_cgns_subs
  use thermo_mod,      only : species_name
  use metrics_mod,     only : time_metrics
  use flowvar_mod,     only : nmarker
  !
  !.. Formal Arguments ..
  character(len=*),     intent(in) :: cgns_sol_file
  integer(CBT),         intent(in) :: ifile
  integer(CBT),         intent(in) :: ibase
  integer(CBT),         intent(in) :: izone
  integer(CBT),         intent(in) :: isol
#ifdef SPECIAL_FOR_INTEL
  character(len=:), allocatable, intent(in) :: variable_names(:)
#else
  character(len=CGLEN), intent(in) :: variable_names(:)
#endif
  integer(CBT),      intent(inout) :: ifld(:)
  !
  !.. Local Scalars ..
  integer :: i,k,l,m,n
  integer :: imarker
  !
  character(len=CGLEN) :: varname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_field_headers"
  character(len=*), parameter :: char_ijk(1:3) = ["Xi  ","Eta ","Zeta"]
  character(len=*), parameter :: char_cyl(1:3) = ["R    ","Theta","Phi  "]
  !
continue
  !
  n = 0
  do k = 1,size(variable_names)
    if (trim(variable_names(k)) == "Velocity") then
      do m = 1,nr
        n = n + 1
        varname = "Velocity"//char_xyz(m)
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
    else if (trim(variable_names(k)) == "RelativeVelocity") then
      do m = 1,nr
        n = n + 1
        varname = "RelativeVelocity"//char_xyz(m)
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
    else if (trim(variable_names(k)) == "CylindricalVelocity") then
      do m = 1,2
        n = n + 1
        varname = "Velocity"//trim(adjustl(char_cyl(m)))
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
    else if (trim(variable_names(k)) == "RCylindricalVelocity") then
      do m = 1,2
        n = n + 1
        varname = "RelativeVelocity"//trim(adjustl(char_cyl(m)))
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
    else if (trim(variable_names(k)) == "Vorticity") then
      do m = 1,nr
        n = n + 1
        varname = "Vorticity"//char_xyz(merge(m,3,nr==3))
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        if (nr < 3) exit
      end do
    else if (trim(variable_names(k)) == "Momentum") then
      do m = 1,nr
        n = n + 1
        varname = "Momentum"//char_xyz(m)
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
    else if (trim(variable_names(k)) == "RelativeMomentum") then
      do m = 1,nr
        n = n + 1
        varname = "RelativeMomentum"//char_xyz(m)
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
    else if (trim(variable_names(k)) == "GradDensity") then
      do l = 1,nr
        n = n + 1
        varname = "DensityGradient"//char_xyz(l)
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
    else if (trim(variable_names(k)) == "GradMomentum") then
      do m = 1,nr
        do l = 1,nr
          n = n + 1
          varname = "Momentum"//char_xyz(m)//"Gradient"//char_xyz(l)
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
      end do
    else if (trim(variable_names(k)) == "GradVelocity") then
      do m = 1,nr
        do l = 1,nr
          n = n + 1
          varname = "Velocity"//char_xyz(m)//"Gradient"//char_xyz(l)
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
      end do
    else if (trim(variable_names(k)) == "GradEnergy") then
      do l = 1,nr
        n = n + 1
        varname = "EnergyGradient"//char_xyz(l)
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
    else if (trim(variable_names(k)) == "GradTemperature") then
      do l = 1,nr
        n = n + 1
        varname = "TemperatureGradient"//char_xyz(l)
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
    else if (trim(variable_names(k)) == "DmatTemperature") then
      do l = 1,nr
        n = n + 1
        varname = "TemperatureDMatrix"//char_xyz(l)
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
    else if (trim(variable_names(k)) == "GradMassFractions") then
      do m = 1,ncs
        do l = 1,nr
          n = n + 1
          varname = trim(adjustl(species_name(m))) // &
                    "MassFractionGradient"//char_xyz(l)
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
      end do
    else if (trim(variable_names(k)) == "GradSpeciesDensities") then
      do m = 1,ncs
        do l = 1,nr
          n = n + 1
          varname = trim(adjustl(species_name(m))) // &
                    "DensityGradient"//char_xyz(l)
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
      end do
    else if (trim(variable_names(k)) == "Viscosity") then
      if (governing_equations == NavierStokes_Eqns) then
        n = n + 1
        varname = trim(variable_names(k))
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end if
    else if (trim(variable_names(k)) == "ResolutionMarker") then
      !
      do imarker = 1, nmarker
        n = n + 1
        varname = trim(variable_names(k))
        if ( imarker > 1 ) write(varname,'(a,i2.2)') trim(varname),imarker
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
      !
      do imarker = 1, nmarker
        n = n + 1
        varname = "Scaled"//trim(variable_names(k))
        if ( imarker > 1 ) write(varname,'(a,i2.2)') trim(varname),imarker
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
      !
    else if (trim(variable_names(k)) == "JumpMarker") then
      !
      do imarker = 1, nmarker
        n = n + 1
        varname = trim(variable_names(k))
        if ( imarker > 1 ) write(varname,'(a,i2.2)') trim(varname),imarker
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
      !
      do imarker = 1, nmarker
        n = n + 1
        varname = "Scaled"//trim(variable_names(k))
        if ( imarker > 1 ) write(varname,'(a,i2.2)') trim(varname),imarker
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
      !
    else if (trim(variable_names(k)) == "DucrosMarker") then
      !
      n = n + 1
      varname = "DucrosFunction"
      call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      !
      n = n + 1
      varname = "DucrosDivergence"
      call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      !
      n = n + 1
      varname = "DucrosSensor"
      call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      !
      do imarker = 1, nmarker
        n = n + 1
        varname = "JumpMarker"
       !varname = trim(variable_names(k))
        if ( imarker > 1 ) write(varname,'(a,i2.2)') trim(varname),imarker
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
      !
      do imarker = 1, nmarker
        n = n + 1
        varname = "Scaled"//trim(variable_names(k))
        if ( imarker > 1 ) write(varname,'(a,i2.2)') trim(varname),imarker
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
      !
    else if (trim(variable_names(k)) == "MLPMarker") then
      !
      do imarker = 1, nmarker
        n = n + 1
        varname = trim(variable_names(k))
        if ( imarker > 1 ) write(varname,'(a,i2.2)') trim(varname),imarker
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
      !
      do imarker = 1, nmarker
        n = n + 1
        varname = "Scaled"//trim(variable_names(k))
        if ( imarker > 1 ) write(varname,'(a,i2.2)') trim(varname),imarker
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
      !
    else if (trim(variable_names(k)) == "FilterDissipation") then
      !
      n = n + 1
      varname = "DensityFilterDissipation"
      call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      !
      do m = 1,nr
        n = n + 1
        varname = "Momentum"//char_xyz(m)//"FilterDissipation"
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end do
      !
      n = n + 1
      varname = "EnergyFilterDissipation"
      call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      !
      n = n + 1
      varname = "MarkerCounter"
      call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      !
    else if (trim(variable_names(k)) == "CellSize") then
      !
      n = n + 1
      varname = "CellVolume"
      call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      !
      n = n + 1
      varname = "CellLengthScale"
      call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      !
    else if (trim(variable_names(k)) == "CoMetrics") then
   !  do m = 1,nr
   !    do i = 1,nr
   !      n = n + 1
   !      varname = "Covariant Metric d" // &
   !                trim(adjustl(char_xyz(m))) // "d" // char_ijk(i)
   !      call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
   !    end do
   !  end do
    else if (trim(variable_names(k)) == "ContraMetrics") then
      do m = 1,nr
        do i = 1,nr
          n = n + 1
          varname = "Contravariant Metric d" // &
                    trim(adjustl(char_ijk(m))) // "d" // char_xyz(i)
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
      end do
    else if (trim(variable_names(k)) == "TimeMetrics") then
      if (allocated(time_metrics)) then
        do m = 1,nr
          n = n + 1
          varname = "Time Metric d" // char_xyz(m)
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
        n = n + 1
        varname = "Time Metric dA"
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      end if
    else if (trim(variable_names(k)) == "MassFractions") then
      if (ncs > 0) then
        do m = 1,ncs
          n = n + 1
          varname = "MassFraction"//trim(adjustl(species_name(m)))
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
      end if
    else if (trim(variable_names(k)) == "MoleFractions") then
      if (ncs > 0) then
        do m = 1,ncs
          n = n + 1
          varname = "MoleFraction"//trim(adjustl(species_name(m)))
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
      end if
    else if (trim(variable_names(k)) == "SpeciesDensities") then
      if (ncs > 0) then
        do m = 1,ncs
          n = n + 1
          varname = "SpeciesDensity"//trim(adjustl(species_name(m)))
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
      end if
    else if (trim(variable_names(k)) == "SpeciesEnthalpy") then
      if (ncs > 0) then
        do m = 1,ncs
          n = n + 1
          varname = "SpeciesEnthalpy"//trim(adjustl(species_name(m)))
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
      end if
    else if (trim(variable_names(k)) == "SpeciesCp") then
      if (ncs > 0) then
        do m = 1,ncs
          n = n + 1
          varname = "SpeciesCp"//trim(adjustl(species_name(m)))
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
      end if
    else if (trim(variable_names(k)) == "cpi/R") then
      if (ncs > 0) then
        do m = 1,ncs
          n = n + 1
          varname = "cpi_over_R_"//trim(adjustl(species_name(m)))
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
      end if
    else if (trim(variable_names(k)) == "hi/R") then
      if (ncs > 0) then
        do m = 1,ncs
          n = n + 1
          varname = "hi_over_R_"//trim(adjustl(species_name(m)))
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
      end if
    else if (trim(variable_names(k)) == "hi/RT") then
      if (ncs > 0) then
        do m = 1,ncs
          n = n + 1
          varname = "hi_over_RT_"//trim(adjustl(species_name(m)))
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
      end if
    else if (trim(variable_names(k)) == "dTdx_Breakdown") then
      if (ncs > 0) then
        do l = 1,nr
          n = n + 1
          varname = "CPG_dTd"//char_xyz(l)
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
        do m = 1,ncs
          n = n + 1
          varname = "Cv("//trim(adjustl(species_name(m)))//")"
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
        do m = 1,ncs
          n = n + 1
          varname = "dCp("//trim(adjustl(species_name(m)))//")dT"
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
        !
        n = n + 1
        varname = "dCp(mix)dT"
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        !
        n = n + 1
        varname = "T_from_gradcv2gradpv"
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        !
        n = n + 1
        varname = "Cp(mix)"
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        !
        n = n + 1
        varname = "Gamma(Mix)"
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        !
        n = n + 1
        varname = "dTdx_lead_coefficient"
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        !
        n = n + 1
        varname = "CPG_constant_Cv"
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        !
        n = n + 1
        varname = "CPG_1_over_constant_Cv"
        call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        !
        do l = 1,nr
          n = n + 1
          varname = "TPG_dTd"//char_xyz(l)
          call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
        end do
        !
        do m = 1,ncs
          do l = 1,nr
            n = n + 1
            varname = "dY("//trim(adjustl(species_name(m)))// &
                      ")d"//char_xyz(l)//"_from_dTd"//char_xyz(l)
            call cgns_write_field_header(n,varname,ifld(n), &
                                         __LINE__,__FILE__)
          end do
        end do
      end if
    else if (trim(variable_names(k)) == "QCriterion") then
      !
      n = n + 1
      varname = trim(variable_names(k))
      call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      !
      n = n + 1
      varname = trim(variable_names(k))//"_Positive"
      call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      !
    else if (trim(variable_names(k)) == "LocalTimestep") then
      !
      n = n + 1
      varname = trim(variable_names(k))
      call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      !
      n = n + 1
      varname = trim(variable_names(k))//"_Scaled"
      call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
      !
    else
      n = n + 1
      if (trim(variable_names(k)) == "Energy") then
        varname = "EnergyStagnationDensity"
      else if (trim(variable_names(k)) == "T0") then
        varname = "TemperatureStagnation"
      else if (trim(variable_names(k)) == "P0") then
        varname = "PressureStagnation"
      else if (trim(variable_names(k)) == "H0") then
        varname = "EnthalpyStagnation"
      else if (trim(variable_names(k)) == "Viscosity") then
        varname = "ViscosityMolecular"
      else if (trim(variable_names(k)) == "Conductivity") then
        varname = "ThermalConductivity"
      else if (trim(variable_names(k)) == "Cp") then
        varname = "SpecificHeatPressure"
      else if (trim(variable_names(k)) == "Gamma") then
        varname = "SpecificHeatRatio"
      else if (trim(variable_names(k)) == "Rgas") then
        varname = "IdealGasConstant"
      else
        varname = trim(variable_names(k))
      end if
      call cgns_write_field_header(n,varname,ifld(n),__LINE__,__FILE__)
    end if
  end do
  !
  !-----------------------------------------------------------------------------
  !
  contains
  !
#include "Functions/cgns_write_field_header.f90"
  !
  !-----------------------------------------------------------------------------
  !
end subroutine write_field_headers
!
!###############################################################################
!
subroutine output_variable(ifile_s,ibase_s,izone_s,isol_s,isvar_s, &
                           cgns_sol_file,vname,nv, &
                           cv,gp,dcvdx,vartmp,var)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr,n_solpts,ncell,cell
  use geomvar_mod, only : xyz
  use ovar,        only : aref,rhoref,time_ref,tref
  use ovar,        only : gam,pref,pref_nd,machref,d_t
  use ovar,        only : rgasref,cpref,Pr
  use ovar,        only : itestcase,governing_equations
  use ovar,        only : profile_io_cgns_subs
  use ovar,        only : Limiter_Option
  use ovar,        only : moving_ref_frame
  use flowvar_mod, only : dtsp,dt_save,cell_diss,cell_num_marks
  !
  use metrics_mod, only : time_metrics,metrics_array,jac,wtdjac
  !
  use flowvar_mod,  only : marker,scaled_marker,nmarker
  use limiters_mod, only : resolution_indicator
  use limiters_mod, only : jump_indicator
  use limiters_mod, only : ducros_indicator
  use limiters_mod, only : cell_lengths
  use flowvar_mod,  only : deallocate_marker
  use flowvar_mod,  only : deallocate_scaled_marker
  use mlp_mod,      only : limit_by_mlp
  !
  use thermo_mod, only : rspec
  use thermo_mod, only : compute_mixture_rgas
  use thermo_mod, only : compute_total_temperature
  use thermo_mod, only : compute_total_pressure
  use thermo_mod, only : compute_species_enthalpy
  use thermo_mod, only : compute_species_cp
  !
  use transport_mod, only : compute_mole_fractions
  !
  use gpu_mod, only : grad_cv_to_grad_pv_gpu
  !
  use taylor_couette_mod, only : remove_freestream_tcf_init
  !
  !.. Formal Arguments ..
  integer(CBT),     intent(in) :: ifile_s
  integer(CBT),     intent(in) :: ibase_s
  integer(CBT),     intent(in) :: izone_s
  integer(CBT),     intent(in) :: isol_s
  integer(CBT),     intent(in) :: isvar_s(:)
  character(len=*), intent(in) :: cgns_sol_file
  character(len=*), intent(in) :: vname
  integer,       intent(inout) :: nv
  real(wp),      intent(inout) :: cv(:,:)
  real(wp),      intent(inout) :: gp(:,:)
  real(wp),         intent(in) :: dcvdx(:,:,:)
  real(wp),      intent(inout) :: vartmp(:)
  real(CRT),     intent(inout) :: var(:)
  !
  !.. Local Scalars ..
  integer  :: i,j,l,m,m1,m2,n,nc,np,n0,n1,n2,npts,ierr
  integer  :: nmi,nmj,l_comp,l_phys,nvar
  integer  :: this_order,limit_order
  integer  :: imarker
  real(wp) :: freestream_ref,rhoa,rhoa2,time_ref_inv
  real(wp) :: cp_dim,cond_dim,dim_val,val,eps
  real(wp) :: rho_inv,vx,vy,vz,a2,q2,s,fval
  real(wp) :: x1,x2,v1,v2,theta,vr,vt
  real(wp) :: gamm,gm1,gam_over_gm1,mach,pres,temp
  real(wp) :: mi,mj,drhodi,drhodj,dmjdi,dmidj
  real(wp) :: w12,w13,w23,wijwij
  real(wp) :: s12,s13,s23,sijsij,s11,s22,s33,skk
  real(wp) :: dudx,dudy,dudz
  real(wp) :: dvdx,dvdy,dvdz
  real(wp) :: dwdx,dwdy,dwdz
  real(wp) :: qmin,qmax,qave,qpts,qsum,qstddev
  real(wp) :: x_r,y_r,z_r,vel_r,r1,r2,vx_r,vy_r,vz_r
  !
  integer(CBT) :: ivar_s
  !
  logical(lk) :: marker_was_not_allocated
  logical(lk) :: scaled_marker_was_not_allocated
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: vartmp_at_sp(:)
  real(wp), allocatable :: xyztmp(:,:)
  real(wp), allocatable :: spcvar(:,:)
  real(wp), allocatable :: dtdx(:,:)
  real(wp), allocatable :: dpvdx(:,:,:)
  real(wp), allocatable :: other_variables(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_variable"
  !
 !real(wp), parameter :: qcriterion_ref_vel = one
  real(wp), parameter :: qcriterion_ref_vel = zero
  !
continue
  !
  npts = size(vartmp)
  !
  freestream_ref = zero
  rhoa = rhoref*aref
  rhoa2 = rhoa*aref
  cp_dim = aref*aref/tref
  cond_dim = rhoa*cp_dim
  !
  ! Pre-load vartmp with Mach number if the current variable
  ! is Mach number, total pressure, or total temperature.
  !
  select case (trim(adjustl(vname)))
    !
    case ("Mach","P0","T0")
      !
      if (ncs == 0) then
        do n = 1,npts
          gp(n,ngamm) = gam
        end do
      end if
      !
      if (nr == 3) then
        do n = 1,npts
          rho_inv = one/cv(n,nec)
          vx = cv(n,nmx)*rho_inv
          vy = cv(n,nmy)*rho_inv
          vz = cv(n,nmz)*rho_inv
          a2 = max( gp(n,ngamm)*gp(n,npres)*rho_inv , rndoff )
          q2 = vx*vx + vy*vy + vz*vz
          vartmp(n) = sqrt( q2 / a2 ) - freestream_ref
        end do
      else if (nr == 2) then
        do n = 1,npts
          rho_inv = one/cv(n,nec)
          vx = cv(n,nmx)*rho_inv
          vy = cv(n,nmy)*rho_inv
          a2 = max( gp(n,ngamm)*gp(n,npres)*rho_inv , rndoff )
          q2 = vx*vx + vy*vy
          vartmp(n) = sqrt( q2 / a2 ) - freestream_ref
        end do
      else
        do n = 1,npts
          rho_inv = one/cv(n,nec)
          vx = cv(n,nmx)*rho_inv
          a2 = max( gp(n,ngamm)*gp(n,npres)*rho_inv , rndoff )
          vartmp(n) = sqrt( vx*vx / a2 ) - freestream_ref
        end do
      end if
      !
    case ("GradVelocity","VorticityMagnitude","StrainRateMagnitude", &
          "QCriterion","DivergenceOfVelocity","DucrosMarker")
      !
      allocate ( dpvdx(1:npts,1:nr,1:nq) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"dpvdx",1,__LINE__,__FILE__,ierr,error_message)
      !
      call grad_cv_to_grad_pv_gpu(1,npts,cv,gp,dcvdx,dpvdx,velocity_only=true)
      !
    case ("GradTemperature","GradMassFractions","GradTempMagnitude")
      !
      allocate ( dpvdx(1:npts,1:nr,1:nq) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"dpvdx",1,__LINE__,__FILE__,ierr,error_message)
      !
      call grad_cv_to_grad_pv_gpu(1,npts,cv,gp,dcvdx,dpvdx)
      !
    case ("dTdx_Breakdown")
      !
      if (ncs > 0) then
        !
        allocate ( dpvdx(1:npts,1:nr,1:nq) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"dpvdx",1,__LINE__,__FILE__,ierr,error_message)
        !
        nvar = nr + 2*ncs + 7
        allocate ( other_variables(1:npts,1:nvar) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"other_variables",1,__LINE__,__FILE__, &
                         ierr,error_message)
        !
        call grad_cv_to_grad_pv_gpu(1,npts,cv,gp,dcvdx,dpvdx, &
                                    output_variables=other_variables)
        !
      end if
      !
    case ("SpeciesEnthalpy","SpeciesCp","hi/RT","hi/R","cpi/R")
      !
      if (ncs > 0) then
        !
        allocate ( spcvar(1:npts,1:ncs) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"spcvar",1,__LINE__,__FILE__,ierr,error_message)
        !
        select case (trim(adjustl(vname)))
          case ("SpeciesCp","cpi/R")
            call compute_species_cp(1,npts,gp(1:npts,ntemp),spcvar, &
                                    dimensional=true)
          case default
            call compute_species_enthalpy(1,npts,gp(1:npts,ntemp),spcvar, &
                                          dimensional=true)
        end select
        !
      end if
      !
    case ("MoleFractions")
      !
      if (ncs > 0) then
        !
        allocate ( spcvar(1:npts,1:ncs) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"spcvar",1,__LINE__,__FILE__,ierr,error_message)
        !
        call compute_mole_fractions(1,npts,cv,spcvar)
        !
      end if
      !
  end select
  !
  if (abs(itestcase) == Freestream_TCF) then
    xyztmp = init_array(npts,nr,zero) ! F2003 AUTO-REALLOC
    do m = 1,nr
      call interpolate_for_output( xyz(m,1:n_solpts) , xyztmp(1:npts,m) )
    end do
  end if
  !
  ! Compute the current output variable and write it to the CGNS file
  !
  variable_select: &
  select case (trim(adjustl(vname)))
    !
    ! #####################
    ! ###  NEED TO ADD  ###
    ! #####################
    ! ---------------------
    !   FOR GAS MIXTURES
    ! ---------------------
    ! sensible enthalpy
    ! ---------------------
    !   FOR CPG AND TPG
    ! ---------------------
    ! pitot pressure
    ! strain rate
    ! pressure gradient
    ! magnitude of pressure gradient
    !
    case ("Density") ! Density
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (itestcase == Freestream_Preservation) freestream_ref = rhoref
      !
      do n = 1,npts
        vartmp(n) = cv(n,nec) * rhoref - freestream_ref
      end do
      !
      if (abs(itestcase) == Freestream_TCF) then
        call remove_freestream_tcf_init(nec,xyztmp,vartmp)
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("Momentum") ! Momentum
      !!{{{
      do m = nmb,nme
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,npts
          vartmp(n) = cv(n,m) * rhoa
        end do
        !
        if (abs(itestcase) == Freestream_TCF) then
          call remove_freestream_tcf_init(m,xyztmp,vartmp)
        end if
        !
        call cgns_write_field_data
        !
      end do
      !!}}}
    case ("RelativeMomentum") ! Relative Momentum
      !!{{{
      !
      ! First we need to get the coordinates for the output points
      !
      if (.not. allocated(xyztmp)) then
        xyztmp = init_array(npts,nr,zero) ! F2003 AUTO-REALLOC
        do m = 1,nr
          call interpolate_for_output( xyz(m,1:n_solpts) , xyztmp(1:npts,m) )
        end do
      end if
      !
      call moving_ref_frame%get_rotation(x_r,y_r,z_r)
      !
      do m = nmb,nme
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        if (m == nmx) then
          r1 = y_r ; m1 = 3 ; r2 = z_r ; m2 = 2
        else if (m == nmy) then
          r1 = z_r ; m1 = 1 ; r2 = x_r ; m2 = 3
        else if (m == nmz) then
          r1 = x_r ; m1 = 2 ; r2 = y_r ; m2 = 1
        end if
        !
        do n = 1,npts
          vel_r = r1*xyztmp(n,m1) - r2*xyztmp(n,m2)
          vartmp(n) = (cv(n,m) - cv(n,nec)*vel_r) * rhoa
        end do
        !
        if (abs(itestcase) == Freestream_TCF) then
          call remove_freestream_tcf_init(m,xyztmp,vartmp)
        end if
        !
        call cgns_write_field_data
        !
      end do
      !!}}}
    case ("Energy") ! Energy
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,npts
        vartmp(n) = cv(n,nee) * rhoa2
      end do
      !
      if (abs(itestcase) == Freestream_TCF) then
        call remove_freestream_tcf_init(nee,xyztmp,vartmp)
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("SpeciesDensities") ! Species densities
      !!{{{
      do m = 1,ncs
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,npts
          vartmp(n) = cv(n,m) * rhoref
        end do
        !
        call cgns_write_field_data
        !
      end do
      !!}}}
    case ("MassFractions") ! Species mass fractions
      !!{{{
      do m = 1,ncs
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,npts
          vartmp(n) = cv(n,m)/cv(n,nec)
        end do
        !
        call cgns_write_field_data(flush_to_zero=fals)
        !
      end do
      !!}}}
    case ("MoleFractions") ! Species mole fractions
      !!{{{
      do m = 1,ncs
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        vartmp(1:npts) = spcvar(1:npts,m)
        !
        call cgns_write_field_data
        !
      end do
      !!}}}
    case ("GradSpeciesDensities") ! Species densities gradients
      !!{{{
      do m = 1,ncs
        do l = 1,nr
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          do n = 1,npts
            vartmp(n) = dcvdx(n,l,m) * rhoref
          end do
          !
          call cgns_write_field_data
          !
        end do
      end do
      !!}}}
    case ("GradDensity") ! Density gradient
      !!{{{
      do l = 1,nr
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,npts
          vartmp(n) = dcvdx(n,l,nec) * rhoref
        end do
        !
        call cgns_write_field_data
        !
      end do
      !!}}}
    case ("GradDensityMagnitude") ! Magnitude of density gradient
      !!{{{
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,npts
        vartmp(n) = norm2( dcvdx(n,1:nr,nec) ) * rhoref
      end do
      !
      call cgns_write_field_data
      !
      !!}}}
    case ("GradMomentum") ! Momentum gradients
      !!{{{
      do m = nmb,nme
        do l = 1,nr
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          do n = 1,npts
            vartmp(n) = dcvdx(n,l,m) * rhoa
          end do
          !
          call cgns_write_field_data
          !
        end do
      end do
      !!}}}
    case ("GradEnergy") ! Total energy gradient
      !!{{{
      do l = 1,nr
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,npts
          vartmp(n) = dcvdx(n,l,nee) * rhoa2
        end do
        !
        call cgns_write_field_data
        !
      end do
      !!}}}
    case ("GradMassFractions") ! Species mass fraction gradients
      !!{{{
      do m = 1,ncs
        do l = 1,nr
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          do n = 1,npts
            vartmp(n) = dpvdx(n,l,m)
          end do
          !
          call cgns_write_field_data
          !
        end do
      end do
      !!}}}
    case ("GradVelocity") ! Velocity gradients
      !!{{{
      do m = nmb,nme
        do l = 1,nr
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          do n = 1,npts
            vartmp(n) = dpvdx(n,l,m) * aref
          end do
          !
          call cgns_write_field_data
          !
        end do
      end do
      !!}}}
    case ("GradTemperature") ! Temperature gradient
      !!{{{
      do l = 1,nr
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,npts
          vartmp(n) = dpvdx(n,l,nee) * tref
        end do
        !
        call cgns_write_field_data
        !
      end do
      !!}}}
    case ("GradTempMagnitude") ! Magnitude of temperature gradient
      !!{{{
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,npts
        vartmp(n) = norm2( dpvdx(n,1:nr,nee) ) * tref
      end do
      !
      call cgns_write_field_data
      !
      !!}}}
    case ("dTdx_Breakdown") ! Breakdown of terms to compute dT/dx
      !!{{{
      if (ncs > 0) then
        do m = 1,nvar
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          do n = 1,npts
            vartmp(n) = other_variables(n,m)
          end do
          !
          call cgns_write_field_data
          !
        end do
        !
        do l = 1,nr
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          do n = 1,npts
            vartmp(n) = dpvdx(n,l,nee)
          end do
          !
          call cgns_write_field_data
          !
        end do
        !
        do m = 1,ncs
          do l = 1,nr
            !
            nv = nv + 1
            ivar_s = isvar_s(nv)
            !
            do n = 1,npts
              vartmp(n) = dpvdx(n,l,m)
            end do
            !
            call cgns_write_field_data
            !
          end do
        end do
      end if
      !!}}}
    case ("DmatTemperature") ! Temperature gradient using derivative matrices
      !!{{{
      !
      allocate ( dtdx(1:n_solpts,1:nr) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"dtdx",1,__LINE__,__FILE__,ierr,error_message)
      !
      call compute_dmatrix_grad_temperature(dtdx)
      !
      if (.not. allocated(vartmp_at_sp)) then
        allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
      !
      do l = 1,nr
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,n_solpts
          vartmp_at_sp(n) = dtdx(n,l) * tref
        end do
        !
        call interpolate_for_output( vartmp_at_sp , vartmp )
        !
        call cgns_write_field_data
        !
      end do
      !
      if (allocated(dtdx)) then
        deallocate ( dtdx , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"dtdx",2,__LINE__,__FILE__,ierr,error_message)
      end if
      !
      !!}}}
    case ("Velocity") ! Velocity
      !!{{{
      do m = nmb,nme
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        if (itestcase == Freestream_Preservation) then
          freestream_ref = merge( aref*machref , zero , m==nmb )
        end if
        !
        do n = 1,npts
          vartmp(n) = cv(n,m)/cv(n,nec) * aref - freestream_ref
        end do
        !
        if (abs(itestcase) == Freestream_TCF) then
          call remove_freestream_tcf_init(-m,xyztmp,vartmp)
        end if
        !
        call cgns_write_field_data
        !
      end do
      !!}}}
    case ("RelativeVelocity") ! Relative Velocity
      !!{{{
      !
      ! First we need to get the coordinates for the output points
      !
      if (.not. allocated(xyztmp)) then
        xyztmp = init_array(npts,nr,zero) ! F2003 AUTO-REALLOC
        do m = 1,nr
          call interpolate_for_output( xyz(m,1:n_solpts) , xyztmp(1:npts,m) )
        end do
      end if
      !
      call moving_ref_frame%get_rotation(x_r,y_r,z_r)
      !
      do m = nmb,nme
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        if (itestcase == Freestream_Preservation) then
          freestream_ref = merge( aref*machref , zero , m==nmb )
        end if
        !
        if (m == nmx) then
          r1 = y_r ; m1 = 3 ; r2 = z_r ; m2 = 2
        else if (m == nmy) then
          r1 = z_r ; m1 = 1 ; r2 = x_r ; m2 = 3
        else if (m == nmz) then
          r1 = x_r ; m1 = 2 ; r2 = y_r ; m2 = 1
        end if
        !
        do n = 1,npts
          vel_r = r1*xyztmp(n,m1) - r2*xyztmp(n,m2)
          vartmp(n) = (cv(n,m)/cv(n,nec) - vel_r) * aref - freestream_ref
        end do
        !
        if (abs(itestcase) == Freestream_TCF) then
          call remove_freestream_tcf_init(-m,xyztmp,vartmp)
        end if
        !
        call cgns_write_field_data
        !
      end do
      !!}}}
    case ("CylindricalVelocity") ! Cylindrical Velocity
      !!{{{
      !
      ! First we need to get the coordinates for the output points
      !
      if (.not. allocated(xyztmp)) then
        xyztmp = init_array(npts,nr,zero) ! F2003 AUTO-REALLOC
        do m = 1,nr
          call interpolate_for_output( xyz(m,1:n_solpts) , xyztmp(1:npts,m) )
        end do
      end if
      !
      ! Get the indices for the primary and secondary directions
      !
      m2 = nr   ! z for 3D, y for 2D
      m1 = m2-1 ! y for 3D, x for 2D
      !
      ! Use vartmp_at_sp to temporarily store Vtheta while Vr is saved
      ! to vartmp. This is so we dont have to do the following loop twice.
      ! After the loop, we will write Vr (in vartmp) to the CGNS solution
      ! file, then we can just copy Vtheta (in vartmp_at_sp) into vartmp
      ! to output Vtheta to the CGNS solution file.
      !
      vartmp_at_sp = init_array(npts,zero) ! F2003 AUTO-REALLOC
      !
      do n = 1,npts
        !
        v1 = cv(n,NM_Vert)/cv(n,nec) * aref
        v2 = cv(n,NM_Span)/cv(n,nec) * aref
        !
        x1 = xyztmp(n,NCD_Vert)
        x2 = xyztmp(n,NCD_Span)
        !
        theta = atan2(x2,x1)
        !
        vr = v2*sin(theta) + v1*cos(theta)
        vt = v2*cos(theta) - v1*sin(theta)
        !
        vartmp(n) = vr
        vartmp_at_sp(n) = vt
        !
      end do
      !
      ! Deallocate xyztmp since it is no longer needed
      !
      if (allocated(xyztmp)) then
        deallocate ( xyztmp , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"xyztmp",2,__LINE__,__FILE__,ierr,error_message)
      end if
      !
      ! Write Vr which is in vartmp to the CGNS solution file
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      call cgns_write_field_data
      !
      ! Copy Vtheta which is in vartmp_to_sp to vartmp
      !
      vartmp(1:npts) = vartmp_at_sp(1:npts)
      !
      ! Write Vtheta which is now in vartmp to the CGNS solution file
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      call cgns_write_field_data
      !!}}}
    case ("RCylindricalVelocity") ! Relative Cylindrical Velocity
      !!{{{
      !
      ! First we need to get the coordinates for the output points
      !
      if (.not. allocated(xyztmp)) then
        xyztmp = init_array(npts,nr,zero) ! F2003 AUTO-REALLOC
        do m = 1,nr
          call interpolate_for_output( xyz(m,1:n_solpts) , xyztmp(1:npts,m) )
        end do
      end if
      !
      ! Get the indices for the primary and secondary directions
      !
      m2 = nr   ! z for 3D, y for 2D
      m1 = m2-1 ! y for 3D, x for 2D
      !
      call moving_ref_frame%get_rotation(x_r,y_r,z_r)
      !
      ! Use vartmp_at_sp to temporarily store Vtheta while Vr is saved
      ! to vartmp. This is so we dont have to do the following loop twice.
      ! After the loop, we will write Vr (in vartmp) to the CGNS solution
      ! file, then we can just copy Vtheta (in vartmp_at_sp) into vartmp
      ! to output Vtheta to the CGNS solution file.
      !
      vartmp_at_sp = init_array(npts,zero) ! F2003 AUTO-REALLOC
      !
      do n = 1,npts
        !
        vy_r = z_r*xyztmp(n,1) - x_r*xyztmp(n,3)
        vz_r = x_r*xyztmp(n,2) - y_r*xyztmp(n,1)
        !
        v1 = (cv(n,nec+m1)/cv(n,nec) - vy_r) * aref
        v2 = (cv(n,nec+m2)/cv(n,nec) - vz_r) * aref
        !
        x1 = xyztmp(n,m1)
        x2 = xyztmp(n,m2)
        !
        theta = atan2(x2,x1)
        !
        vr = v2*sin(theta) + v1*cos(theta)
        vt = v2*cos(theta) - v1*sin(theta)
        !
        vartmp(n) = vr
        vartmp_at_sp(n) = vt
        !
      end do
      !
      ! Deallocate xyztmp since it is no longer needed
      !
      if (allocated(xyztmp)) then
        deallocate ( xyztmp , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"xyztmp",2,__LINE__,__FILE__,ierr,error_message)
      end if
      !
      ! Write Vr which is in vartmp to the CGNS solution file
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      call cgns_write_field_data
      !
      ! Copy Vtheta which is in vartmp_to_sp to vartmp
      !
      vartmp(1:npts) = vartmp_at_sp(1:npts)
      !
      ! Write Vtheta which is now in vartmp to the CGNS solution file
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      call cgns_write_field_data
      !!}}}
    case ("Pressure") ! Pressure
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (itestcase == Freestream_Preservation) freestream_ref = pref
      !
      do n = 1,npts
        vartmp(n) = gp(n,npres) * rhoa2 - freestream_ref
      end do
      !
      if (abs(itestcase) == Freestream_TCF) then
        call remove_freestream_tcf_init(nq+npres,xyztmp,vartmp)
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("P0") ! Total Pressure
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (ncs > 0) then
        !
        call compute_total_pressure(1,npts,cv,gp,vartmp,dimensional=true)
        !
      else
        !
        ! NOTE: vartmp has been pre-loaded with the Mach number, and the
        !       auxiliary variable in gp has already been switched to gamma
        !
        do n = 1,npts
          mach = vartmp(n)
          pres = gp(n,npres) * rhoa2
          gamm = gp(n,ngamm)
          gm1 = gamm - one
          gam_over_gm1 = gamm / gm1
          vartmp(n) = pres*( (one + half*gm1*mach*mach)**gam_over_gm1 )
        end do
        !
        if (abs(itestcase) == Freestream_TCF) then
          call remove_freestream_tcf_init(-nq-npres,xyztmp,vartmp)
        end if
        !
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("Mach") ! Mach number
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      ! NOTE: vartmp has been pre-loaded with the Mach number
      !
      if (itestcase == Freestream_Preservation) then
        do n = 1,npts
          vartmp(n) = vartmp(n) - machref
        end do
      end if
      !
      if (abs(itestcase) == Freestream_TCF) then
        call remove_freestream_tcf_init(-nee,xyztmp,vartmp)
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("RelativeMach") ! Relative Mach number
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      ! NOTE: vartmp has been pre-loaded with the Mach number
      !
      if (nr == 3) then
        !
        ! First we need to get the coordinates for the output points
        !
        if (.not. allocated(xyztmp)) then
          xyztmp = init_array(npts,nr,zero) ! F2003 AUTO-REALLOC
          do m = 1,nr
            call interpolate_for_output( xyz(m,1:n_solpts) , xyztmp(1:npts,m) )
          end do
        end if
        !
        call moving_ref_frame%get_rotation(x_r,y_r,z_r)
        !
        do n = 1,npts
          !
          rho_inv = one/cv(n,nec)
          !
          vx_r = y_r*xyztmp(n,3) - z_r*xyztmp(n,2)
          vy_r = z_r*xyztmp(n,1) - x_r*xyztmp(n,3)
          vz_r = x_r*xyztmp(n,2) - y_r*xyztmp(n,1)
          !
          vx = cv(n,nmx)*rho_inv - vx_r
          vy = cv(n,nmy)*rho_inv - vy_r
          vz = cv(n,nmz)*rho_inv - vz_r
          
          a2 = max( gp(n,ngamm)*gp(n,npres)*rho_inv , rndoff )
          q2 = vx*vx + vy*vy + vz*vz
          vartmp(n) = sqrt( q2 / a2 ) - freestream_ref
          !
        end do
        !
      end if
      !
      if (itestcase == Freestream_Preservation) then
        do n = 1,npts
          vartmp(n) = vartmp(n) - machref
        end do
      end if
      !
      if (abs(itestcase) == Freestream_TCF) then
        call remove_freestream_tcf_init(-nee,xyztmp,vartmp)
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("Entropy") ! Entropy
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (ncs == 0) then
        do n = 1,npts
          gp(n,ngamm) = gam
        end do
      end if
      !
      if (itestcase == Inviscid_Gaussian_Bump) then
        do n = 1,npts
          s = (gp(n,npres)/pref_nd) / (cv(n,nec)**gp(n,ngamm))
          vartmp(n) = abs( one - s )
        end do
      else
        do n = 1,npts
          vartmp(n) = log(gp(n,npres)) - gp(n,ngamm)*log(cv(n,nec))
        end do
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("Temperature") ! Temperature
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (itestcase == Freestream_Preservation) freestream_ref = tref
      !
      do n = 1,npts
        vartmp(n) = gp(n,ntemp) * tref - freestream_ref
      end do
      !
      if (abs(itestcase) == Freestream_TCF) then
        call remove_freestream_tcf_init(nq+ntemp,xyztmp,vartmp)
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("T0") ! Total Temperature
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (ncs > 0) then
        !
        call compute_total_temperature(1,npts,cv,gp,vartmp,dimensional=true)
        !
      else
        !
        ! NOTE: vartmp has been pre-loaded with the Mach number, and the
        !       auxiliary variable in gp has already been switched to gamma
        !
        do n = 1,npts
          mach = vartmp(n)
          temp = gp(n,ntemp) * tref
          gm1 = gp(n,ngamm) - one
          vartmp(n) = temp*(one + half*gm1*mach*mach)
        end do
        !
        if (abs(itestcase) == Freestream_TCF) then
          call remove_freestream_tcf_init(-nq-ntemp,xyztmp,vartmp)
        end if
        !
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("Enthalpy") ! Enthalpy
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (ncs > 0) then
        do n = 1,npts
          vartmp(n) = gp(n,nenth) * aref*aref
        end do
      else
        gam_over_gm1 = gam / (gam-one)
        do n = 1,npts
          vartmp(n) = gam_over_gm1*gp(n,npres)/cv(n,nec) * aref*aref
        end do
      end if
      !
      if (abs(itestcase) == Freestream_TCF) then
        call remove_freestream_tcf_init(nq+ngamm,xyztmp,vartmp)
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("H0") ! Total Enthalpy
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (ncs > 0) then
        do n = 1,npts
          vartmp(n) = gp(n,nenth) * aref*aref
        end do
      else
        gam_over_gm1 = gam / (gam-one)
        do n = 1,npts
          vartmp(n) = gam_over_gm1*gp(n,npres)/cv(n,nec) * aref*aref
        end do
      end if
      !
      if (nr == 3) then
        do n = 1,npts
          rho_inv = one/cv(n,nec)
          vx = cv(n,nmx)*rho_inv * aref
          vy = cv(n,nmy)*rho_inv * aref
          vz = cv(n,nmz)*rho_inv * aref
          q2 = vx*vx + vy*vy + vz*vz
          vartmp(n) = vartmp(n) + half*q2
        end do
      else if (nr == 2) then
        do n = 1,npts
          rho_inv = one/cv(n,nec)
          vx = cv(n,nmx)*rho_inv * aref
          vy = cv(n,nmy)*rho_inv * aref
          q2 = vx*vx + vy*vy
          vartmp(n) = vartmp(n) + half*q2
        end do
      else
        do n = 1,npts
          vx = cv(n,nmx)/cv(n,nec) * aref
          vartmp(n) = vartmp(n) + half*vx*vx
        end do
      end if
      !
      if (abs(itestcase) == Freestream_TCF) then
        call remove_freestream_tcf_init(-nq-ngamm,xyztmp,vartmp)
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("Rgas") ! Mixture gas constant
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (ncs > 0) then
        do m = 1,ncs
          do n = 1,npts
            cv(n,m) = cv(n,m)/cv(n,nec)
          end do
        end do
        call compute_mixture_rgas(1,npts,cv(1:npts,1:ncs),vartmp, &
                                  dimensional=true)
        do m = 1,ncs
          do n = 1,npts
            cv(n,m) = cv(n,m)*cv(n,nec)
          end do
        end do
      else
        do n = 1,npts
          vartmp(n) = rgasref
        end do
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("Gamma") ! Ratio of specific heats
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (ncs > 0) then
        do n = 1,npts
          vartmp(n) = gp(n,ngamm)
        end do
      else
        do n = 1,npts
          vartmp(n) = gam
        end do
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("Cp") ! Specific heat at constant pressure
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (ncs > 0) then
        do n = 1,npts
          vartmp(n) = gp(n,ncp) * cp_dim
        end do
      else
        do n = 1,npts
          vartmp(n) = cpref
        end do
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("Viscosity") ! Viscosity
      !!{{{
      if (governing_equations == NavierStokes_Eqns) then
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,npts
          vartmp(n) = gp(n,nvisc) * rhoa
        end do
        !
        call cgns_write_field_data
        !
      end if
      !!}}}
    case ("Conductivity") ! Conductivity
      !!{{{
      if (governing_equations == NavierStokes_Eqns) then
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,npts
          vartmp(n) = gp(n,ncond) * cond_dim
        end do
        !
        call cgns_write_field_data
        !
      end if
      !!}}}
    case ("Prandtl") ! Mixture Prandtl number
      !!{{{
      if (governing_equations == NavierStokes_Eqns) then
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        if (ncs > 0) then
          do n = 1,npts
            vartmp(n) = gp(n,nvisc) * gp(n,ncp) / gp(n,ncond)
          end do
        else
          do n = 1,npts
            vartmp(n) = Pr
          end do
        end if
        !
        call cgns_write_field_data
        !
      end if
      !!}}}
    case ("SpeciesEnthalpy") ! Species enthalpies
      !!{{{
      if (ncs > 0) then
        do m = 1,ncs
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          vartmp(1:npts) = spcvar(1:npts,m)
          !
          call cgns_write_field_data
          !
        end do
      end if
      !!}}}
    case ("SpeciesCp") ! Species specific heats at constant pressure
      !!{{{
      if (ncs > 0) then
        do m = 1,ncs
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          vartmp(1:npts) = spcvar(1:npts,m)
          !
          call cgns_write_field_data
          !
        end do
      end if
      !!}}}
    case ("cpi/R") ! Species Cp / Species gas constant
      !!{{{
      if (ncs > 0) then
        do m = 1,ncs
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          dim_val = one/(rspec(m)*cp_dim)
          vartmp(1:npts) = spcvar(1:npts,m)*dim_val
          !
          call cgns_write_field_data
          !
        end do
      end if
      !!}}}
    case ("hi/RT") ! Species enthalpies / (Species gas constant * temperature)
      !!{{{
      if (ncs > 0) then
        do m = 1,ncs
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          dim_val = rspec(m)*cp_dim*tref
          do n = 1,npts
            vartmp(n) = spcvar(n,m)/(gp(n,ntemp)*dim_val)
          end do
          !
          call cgns_write_field_data
          !
        end do
      end if
      !!}}}
    case ("hi/R") ! Species enthalpies / Species gas constant
      !!{{{
      if (ncs > 0) then
        do m = 1,ncs
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          dim_val = one/(rspec(m)*cp_dim)
          vartmp(1:npts) = spcvar(1:npts,m)*dim_val
          !
          call cgns_write_field_data
          !
        end do
      end if
      !!}}}
    case ("Vorticity") ! Vorticity
      !!{{{
      time_ref_inv = one/time_ref
      !
      vorticity_loop: &
      do l = 1,nr
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        ! [2-(2/4)*3,3-(3/4)*3,4-(4/4)*3] = [2-0*3,3-0*3,4-1*3] = [2,3,1]
        i = (l+1) - ( (l+1)/4 )*3
        ! [3-(3/4)*3,4-(4/4)*3,5-(5/4)*3] = [3-0*3,4-1*3,5-1*3] = [3,1,2]
        j = (l+2) - ( (l+2)/4 )*3
        !
        nmi = nmb-1+i ! [y,z,x]
        nmj = nmb-1+j ! [z,x,y]
        !
        do n = 1,npts
          !
          rho_inv = one/cv(n,nec)
          !
          mi = cv(n,nmi)
          mj = cv(n,nmj)
          !
          drhodi = dcvdx(n,i,nec) * time_ref_inv
          drhodj = dcvdx(n,j,nec) * time_ref_inv
          dmjdi = dcvdx(n,i,nmj) * time_ref_inv
          dmidj = dcvdx(n,j,nmi) * time_ref_inv
          !
          vartmp(n) = rho_inv*( (dmjdi-dmidj) - rho_inv*(mj*drhodi-mi*drhodj) )
          !
        end do
        !
        call cgns_write_field_data
        !
        if (nr < 3) exit vorticity_loop
        !
      end do vorticity_loop
      !!}}}
    case ("VorticityMagnitude") ! Vorticity magnitude
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,npts
        !
        w12 = dpvdx(n,2,nmx) - dpvdx(n,1,nmy)
        w13 = dpvdx(n,3,nmx) - dpvdx(n,1,nmz)
        w23 = dpvdx(n,3,nmy) - dpvdx(n,2,nmz)
        !
        vartmp(n) = sqrt( w12*w12 + w13*w13 + w23*w23 ) * aref
        !
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("StrainRateMagnitude") ! Strain rate magnitude
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,npts
        !
        dudx = dpvdx(n,1,nmx)
        dudy = dpvdx(n,2,nmx)
        dudz = dpvdx(n,3,nmx)
        !
        dvdx = dpvdx(n,1,nmy)
        dvdy = dpvdx(n,2,nmy)
        dvdz = dpvdx(n,3,nmy)
        !
        dwdx = dpvdx(n,1,nmz)
        dwdy = dpvdx(n,2,nmz)
        dwdz = dpvdx(n,3,nmz)
        !
        skk = one3 * (dudx + dvdy + dwdz)
        !
        s11 = dudx - skk
        s22 = dvdy - skk
        s33 = dwdz - skk
        !
        s12 = half * (dudy + dvdx)
        s13 = half * (dudz + dwdx)
        s23 = half * (dvdz + dwdy)
        !
        sijsij = s11**2 + s22**2 + s33**2 + two*(s12**2 + s13**2 + s23**2)
        !
        vartmp(n) = sqrt( two * sijsij ) * aref
        !
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("QCriterion") ! Q-criterion
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      ! Q = 0.5*( div(u)**2 + ||rotation||^2 - ||strain||^2 )
      !
      ! div(u) = dudx + dvdy + dwdz
      ! rotation[x,y,z] = 0.5*[dwdy-dvdz,dudz-dwdx,dvdx-dudy]
      !
      call moving_ref_frame%get_rotation(x_r,y_r,z_r)
      !
      do n = 1,npts
        !
        dudx = dpvdx(n,1,nmx)
        dudy = dpvdx(n,2,nmx) + z_r*qcriterion_ref_vel
        dudz = dpvdx(n,3,nmx) - y_r*qcriterion_ref_vel
        !
        dvdx = dpvdx(n,1,nmy) - z_r*qcriterion_ref_vel
        dvdy = dpvdx(n,2,nmy)
        dvdz = dpvdx(n,3,nmy) + x_r*qcriterion_ref_vel
        !
        dwdx = dpvdx(n,1,nmz) + y_r*qcriterion_ref_vel
        dwdy = dpvdx(n,2,nmz) - x_r*qcriterion_ref_vel
        dwdz = dpvdx(n,3,nmz)
        !
        w12 = half * (dudy - dvdx)
        w13 = half * (dudz - dwdx)
        w23 = half * (dvdz - dwdy)
        !
        wijwij = two * (w12**2 + w13**2 + w23**2)
        !
        skk = one3 * (dudx + dvdy + dwdz)
        !
        s11 = dudx - skk
        s22 = dvdy - skk
        s33 = dwdz - skk
        !
        s12 = half * (dudy + dvdx)
        s13 = half * (dudz + dwdx)
        s23 = half * (dvdz + dwdy)
        !
        sijsij = s11**2 + s22**2 + s33**2 + two*(s12**2 + s13**2 + s23**2)
        !
        vartmp(n) = half*(wijwij - sijsij)
       !vartmp(n) = vartmp(n) * aref*aref
        !
      end do
      !
      ! Save the original Q-criterion
      !
      vartmp_at_sp = vartmp(:) ! F2003 AUTO-REALLOC
      !
      ! Get the min/max values of the Q-criterion on this partition
      !
      qmin = minval(vartmp(:))
      qmax = maxval(vartmp(:))
      !
      ! Get the collective min/max Q-criterion values across all partitions
      !
      !##################################################################
      !##################################################################
      !##################################################################
      ! NOTE: IF WE EVER CHANGE CGNS OUTPUT TO SOMETHING LIKE RESTART
      !       FILES WHERE ONLY THE HOST MPI PROCESS ON EACH NODE
      !       PARTICIPATES IN WRITING THE CGNS FILE (I.E., SEE INPUT
      !       RESTART_OUTPUT_USES_HOST_ROOTS), WE WILL HAVE TO CHANGE
      !       THE COMMUNICATOR USED IN THE FOLLOWING MPI CALLS.
      !##################################################################
      !##################################################################
      !##################################################################
      !
      call mpi_allreduce(MPI_IN_PLACE,qmin,1_int_mpi,mpi_flttyp, &
                         MPI_MIN,MPI_COMM_WORLD,mpierr)
      !
      call mpi_allreduce(MPI_IN_PLACE,qmax,1_int_mpi,mpi_flttyp, &
                         MPI_MAX,MPI_COMM_WORLD,mpierr)
      !
      ! Normalize the Q-criterion to the interval [0,1]
      !
      if (abs(qmax-qmin) > eps9) then
        vartmp(:) = (vartmp(:) - qmin)/(qmax-qmin)
      else
        vartmp(:) = vartmp(:) - qmin
      end if
      !
#define STDDEV_NORMALIZED_QCRITERION
#ifdef STDDEV_NORMALIZED_QCRITERION
      !
      ! Compute the average of the normalized Q-criterion
      !
      qave = sum(vartmp(:))
      qpts = real(size(vartmp),kind=kind(qpts))
      !
      call mpi_allreduce(MPI_IN_PLACE,qave,1_int_mpi,mpi_flttyp, &
                         MPI_SUM,MPI_COMM_WORLD,mpierr)
      !
      call mpi_allreduce(MPI_IN_PLACE,qpts,1_int_mpi,mpi_flttyp, &
                         MPI_SUM,MPI_COMM_WORLD,mpierr)
      !
      qave = qave / qpts
      !
      ! Compute the standard deviation
      !
      qsum = zero
      do n = 1,npts
        qsum = qsum + (vartmp(n)-qave)**2
      end do
      !
      call mpi_allreduce(MPI_IN_PLACE,qsum,1_int_mpi,mpi_flttyp, &
                         MPI_SUM,MPI_COMM_WORLD,mpierr)
      !
      qstddev = sqrt( qsum / (qpts-one) )
      !
      ! Scale the Q-criterion so that the average value is 0, and plus/minus
      ! one is plus/minus one/100 standard deviation from the average.
      !
      if (qstddev > eps14) then
        vartmp(:) = (vartmp(:) - qave) / (eps2*qstddev)
      else
        vartmp(:) = vartmp(:) - qave
      end if
#endif
      !
      call cgns_write_field_data
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      ! Copy the original positive values of Q-Criterion back to vartmp
      ! and make the negative values zero
      !
      vartmp(:) = merge( vartmp_at_sp(:) , zero , vartmp_at_sp(:) >= zero )
      !
      call cgns_write_field_data
      !
      !!}}}
    case ("DivergenceOfVelocity") ! Divergence of velocity
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,npts
        vartmp(n) = (dpvdx(n,1,nmx) + dpvdx(n,2,nmy) + dpvdx(n,3,nmz)) * aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("ResolutionMarker")
      !!{{{
      if (.not. allocated(vartmp_at_sp)) then
        allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
      !
      ! Compute the resolution marker if it has yet to be computed
      ! or was computed using a different marker type.
      !
      marker_was_not_allocated        = ( .not. allocated(marker) )
      scaled_marker_was_not_allocated = ( .not. allocated(scaled_marker) )
      !
      if (  all(Limiter_Option /= [-1,1]) .or. &
                 marker_was_not_allocated .or. &
          scaled_marker_was_not_allocated ) then
        ! NOTE: Passing orders_to_cutoff=0 ensures that only the highest mode
        !       is cutoff to compute the resolution marker. We want it to use
        !       the P order cutoff filter, and without this, it will instead
        !       use the P-1 order cutoff filter. This was causing a crash when
        !       P=1, since the cutoff filter isnt allocated with a P=0 cutoff
        !       filter.
        call resolution_indicator(nst=1,apply_projection=fals, &
                                  orders_to_cutoff=0)
      end if
      !
      ! Actual resolution marker value
      !
      do imarker = 1, nmarker
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do nc = 1,ncell
          n0 = cell(nc)%beg_sp - 1
          np = cell(nc)%end_sp - n0
          do n = 1,np
            vartmp_at_sp(n0+n) = marker(nc,imarker)
          end do
        end do
        !
        call interpolate_for_output( vartmp_at_sp , vartmp )
        !
        call cgns_write_field_data
        !
      end do
      !
      ! Scaled resolution marker value
      !
      do imarker = 1, nmarker
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        if (Limiter_Option <= -1 .and. Limiter_Option /= -2) then
          do nc = 1,ncell
            n0 = cell(nc)%beg_sp - 1
            np = cell(nc)%end_sp - n0
            this_order = cell(nc)%order
            limit_order = int( real(this_order,kind=wp) * &
                               (one-scaled_marker(nc,imarker)) )
            limit_order = max(1,min(this_order,limit_order))
            do n = 1,np
              vartmp_at_sp(n0+n) = limit_order
            end do
          end do
        else
          do nc = 1,ncell
            n0 = cell(nc)%beg_sp - 1
            np = cell(nc)%end_sp - n0
            do n = 1,np
              vartmp_at_sp(n0+n) = scaled_marker(nc,imarker)
            end do
          end do
        end if
        !
        call interpolate_for_output( vartmp_at_sp , vartmp )
        !
        call cgns_write_field_data
        !
      end do
      !
      ! Deallocate marker and/or scaled_marker if they were not originally
      ! allocated when this case construct was entered.
      !
      if (marker_was_not_allocated) then
        call deallocate_marker(pname)
      end if
      !
      if (scaled_marker_was_not_allocated) then
        call deallocate_scaled_marker(pname)
      end if
      !
      !!}}}
    case ("JumpMarker")
      !!{{{
      allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      ! Compute the pressure jump marker if it has yet to be computed
      !
      marker_was_not_allocated        = ( .not. allocated(marker) )
      scaled_marker_was_not_allocated = ( .not. allocated(scaled_marker) )
      !
      if (       abs(Limiter_Option) /= 2 .or. &
                 marker_was_not_allocated .or. &
          scaled_marker_was_not_allocated ) then
        call jump_indicator(nst=1)
      end if
      !
      ! Actual pressure jump marker value
      !
      do imarker = 1, nmarker
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do nc = 1,ncell
          n0 = cell(nc)%beg_sp - 1
          np = cell(nc)%end_sp - n0
          do n = 1,np
            vartmp_at_sp(n0+n) = marker(nc,imarker)
          end do
        end do
        !
        call interpolate_for_output( vartmp_at_sp , vartmp )
        !
        call cgns_write_field_data
      end do
      !
      ! Scaled pressure jump marker value
      !
      do imarker = 1, nmarker
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do nc = 1,ncell
          n0 = cell(nc)%beg_sp - 1
          np = cell(nc)%end_sp - n0
          do n = 1,np
            vartmp_at_sp(n0+n) = scaled_marker(nc,imarker)
          end do
        end do
        !
        call interpolate_for_output( vartmp_at_sp , vartmp )
        !
        call cgns_write_field_data
      end do
      !
      ! Deallocate marker and/or scaled_marker if they were not originally
      ! allocated when this case construct was entered.
      !
      if (marker_was_not_allocated) then
        call deallocate_marker(pname)
      end if
      !
      if (scaled_marker_was_not_allocated) then
        call deallocate_scaled_marker(pname)
      end if
      !
      !!}}}
    case ("DucrosMarker")
      !!{{{
      allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      ! Compute the Ducros marker if it has yet to be computed
      !
      marker_was_not_allocated        = ( .not. allocated(marker) )
      scaled_marker_was_not_allocated = ( .not. allocated(scaled_marker) )
      !
      if (       abs(Limiter_Option) /= 4 .or. &
                 marker_was_not_allocated .or. &
          scaled_marker_was_not_allocated ) then
        call ducros_indicator(nst=1,update_gradients=fals)
      end if
      !
      ! Compute the value of the Ducros function at each solution point
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do nc = 1,ncell
        n1 = cell(nc)%beg_sp
        n2 = cell(nc)%end_sp
        eps = eps8 * machref / cell_lengths(nc)
        do n = n1,n2
          !
          skk = dpvdx(n,1,nmx) + dpvdx(n,2,nmy) + dpvdx(n,3,nmz)
          !
          w12 = dpvdx(n,2,nmx) - dpvdx(n,1,nmy)
          w13 = dpvdx(n,3,nmx) - dpvdx(n,1,nmz)
          w23 = dpvdx(n,3,nmy) - dpvdx(n,2,nmz)
          !
          wijwij = w12*w12 + w13*w13 + w23*w23
          !
          vartmp_at_sp(n) = (skk*skk + eps) / (skk*skk + wijwij + eps)
          !
        end do
      end do
      !
      call interpolate_for_output( vartmp_at_sp , vartmp )
      !
      call cgns_write_field_data
      !
      ! Compute the value of the Ducros 'S' value at each solution point
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do nc = 1,ncell
        n1 = cell(nc)%beg_sp
        n2 = cell(nc)%end_sp
        val = half*cell_lengths(nc)
        do n = n1,n2
          !
          skk = dpvdx(n,1,nmx) + dpvdx(n,2,nmy) + dpvdx(n,3,nmz)
          !
          a2 = gp(n,ngamm)*gp(n,npres)/max(cv(n,nec),rndoff)
          !
          vartmp_at_sp(n) = abs(skk) * val / sqrt(a2)
          !
        end do
      end do
      !
      call interpolate_for_output( vartmp_at_sp , vartmp )
      !
      call cgns_write_field_data
      !
      ! Compute the value of the actual Ducros sensor at each solution point
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do nc = 1,ncell
        n1 = cell(nc)%beg_sp
        n2 = cell(nc)%end_sp
        val = half*cell_lengths(nc)
        eps = eps8 * machref / cell_lengths(nc)
        do n = n1,n2
          !
          skk = dpvdx(n,1,nmx) + dpvdx(n,2,nmy) + dpvdx(n,3,nmz)
          !
          w12 = dpvdx(n,2,nmx) - dpvdx(n,1,nmy)
          w13 = dpvdx(n,3,nmx) - dpvdx(n,1,nmz)
          w23 = dpvdx(n,3,nmy) - dpvdx(n,2,nmz)
          !
          wijwij = w12*w12 + w13*w13 + w23*w23
          !
          a2 = gp(n,ngamm)*gp(n,npres)/max(cv(n,nec),rndoff)
          !
          fval = (skk*skk + eps) / (skk*skk + wijwij + eps)
          fval = min( one , fval )

          s = min( one , abs(skk)*val/sqrt(a2) )
          !
          vartmp_at_sp(n) = min( one , max( s*s , min(0.01_wp,fval) ) )
          !
        end do
      end do
      !
      call interpolate_for_output( vartmp_at_sp , vartmp )
      !
      call cgns_write_field_data
      !
      ! Value of the jump marker
      !
      do imarker = 1, nmarker
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do nc = 1,ncell
          n0 = cell(nc)%beg_sp - 1
          np = cell(nc)%end_sp - n0
          do n = 1,np
            vartmp_at_sp(n0+n) = marker(nc,imarker)
          end do
        end do
        !
        call interpolate_for_output( vartmp_at_sp , vartmp )
        !
        call cgns_write_field_data
      end do
      !
      ! Scaled Ducros marker which is the Ducros sensor multiplied
      ! by the scaled jump_marker
      !
      do imarker = 1, nmarker
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do nc = 1,ncell
          n0 = cell(nc)%beg_sp - 1
          np = cell(nc)%end_sp - n0
          do n = 1,np
            vartmp_at_sp(n0+n) = scaled_marker(nc,imarker)
          end do
        end do
        !
        call interpolate_for_output( vartmp_at_sp , vartmp )
        !
        call cgns_write_field_data
      end do
      !
      ! Deallocate marker and/or scaled_marker if they were not originally
      ! allocated when this case construct was entered.
      !
      if (marker_was_not_allocated) then
        call deallocate_marker(pname)
      end if
      !
      if (scaled_marker_was_not_allocated) then
        call deallocate_scaled_marker(pname)
      end if
      !
      !!}}}
    case ("MLPMarker")
      !!{{{
      if (.not. allocated(vartmp_at_sp)) then
        allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
      !
      ! Compute the mlp marker if it has yet to be computed
      !
      marker_was_not_allocated        = ( .not. allocated(marker) )
      scaled_marker_was_not_allocated = ( .not. allocated(scaled_marker) )
      !
      if (            Limiter_Option /= 3 .or. &
                 marker_was_not_allocated .or. &
          scaled_marker_was_not_allocated ) then
        ! Use relax=0 so solution is not updated with limited values.
        call limit_by_mlp(nst=1,relax=zero)
      end if
      !
      ! Actual marker value
      !
      do imarker = 1, nmarker
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do nc = 1,ncell
          n0 = cell(nc)%beg_sp - 1
          np = cell(nc)%end_sp - n0
          do n = 1,np
            vartmp_at_sp(n0+n) = marker(nc,imarker)
          end do
        end do
        !
        call interpolate_for_output( vartmp_at_sp , vartmp )
        !
        call cgns_write_field_data
        !
      end do
      !
      ! Scaled marker value
      !
      do imarker = 1, nmarker
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do nc = 1,ncell
          n0 = cell(nc)%beg_sp - 1
          np = cell(nc)%end_sp - n0
          do n = 1,np
            vartmp_at_sp(n0+n) = scaled_marker(nc,imarker)
          end do
        end do
        !
        call interpolate_for_output( vartmp_at_sp , vartmp )
        !
        call cgns_write_field_data
        !
      end do
      !
      ! Deallocate marker and/or scaled_marker if they were not originally
      ! allocated when this case construct was entered.
      !
      if (marker_was_not_allocated) then
        call deallocate_marker(pname)
      end if
      !
      if (scaled_marker_was_not_allocated) then
        call deallocate_scaled_marker(pname)
      end if
      !!}}}
    case ("FilterDissipation")
      !!{{{
      !
      if (.not. allocated(vartmp_at_sp)) then
        allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
      !
      do m = 1,nq
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        ! Save the limiter filter dissipation for this equation to vartmp_at_sp
        !
        do nc = 1,ncell
          n0 = cell(nc)%beg_sp - 1
          np = cell(nc)%end_sp - n0
          do n = 1,np
            vartmp_at_sp(n0+n) = cell_diss(nc,m)
          end do
        end do
        !
        call interpolate_for_output( vartmp_at_sp , vartmp )
        !
        call cgns_write_field_data
        !
      end do
      !
      ! Save the number of times each cell has been marked to vartmp_at_sp
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do nc = 1,ncell
        n0 = cell(nc)%beg_sp - 1
        np = cell(nc)%end_sp - n0
        do n = 1,np
          vartmp_at_sp(n0+n) = cell_num_marks(nc)
        end do
      end do
      !
      call interpolate_for_output( vartmp_at_sp , vartmp )
      !
      call cgns_write_field_data
      !
      !!}}}
    case ("Partition")
      !!{{{
      allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,n_solpts
        vartmp_at_sp(n) = mypnum+1
      end do
      !
      call interpolate_for_output( vartmp_at_sp , vartmp )
      !
      call cgns_write_field_data
      !!}}}
    case ("CellSize")
      !!{{{
      allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      ! Compute and output the cell volumes
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do nc = 1,ncell
        n1 = cell(nc)%beg_sp
        n2 = cell(nc)%end_sp
        val = sum( wtdjac(n1:n2) )
        do n = n1,n2
          vartmp_at_sp(n) = val
        end do
      end do
      !
      call interpolate_for_output( vartmp_at_sp , vartmp )
      !
      ! Compute and output the cell edge lengths (vol**(1/3))
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,n_solpts
        vartmp_at_sp(n) = vartmp_at_sp(n)**one3
      end do
      !
      call interpolate_for_output( vartmp_at_sp , vartmp )
      !
      call cgns_write_field_data
      !!}}}
    case ("Timestep")
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      call interpolate_for_output( dtsp(1:n_solpts) , vartmp , nondim=time_ref )
      !
      call cgns_write_field_data
      !!}}}
    case ("LocalTimestep")
      !!{{{
      !
      ! Write out the actual dimensional time step size
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      call interpolate_for_output( dt_save(1:n_solpts), vartmp, nondim=time_ref)
      !
      call cgns_write_field_data(flush_to_zero=fals)
      !
      ! Write out the time step sizes scaled by the minimum time step
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      dt_save(:) = dt_save(:) / d_t
      call interpolate_for_output( dt_save(1:n_solpts), vartmp, nondim=time_ref)
      !
      call cgns_write_field_data
      !!}}}
    case ("CoMetrics")
     !!!{{{
     !allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
     !           stat=ierr , errmsg=error_message )
     !call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
     !                 ierr,error_message)
     !!
     !do l_phys = 1,nr
     !  do l_comp = 1,nr
     !    !
     !    nv = nv + 1
     !    ivar_s = isvar_s(nv)
     !    !
     !    do n = 1,n_solpts
     !      vartmp_at_sp(n) = metrics_array(n,l_phys,l_comp)
     !    end do
     !    !
     !    call interpolate_for_output( vartmp_at_sp , vartmp )
     !    !
     !    call cgns_write_field_data
     !    !
     !  end do
     !end do
     !!!}}}
    case ("ContraMetrics")
      !!{{{
      allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      do l_comp = 1,nr
        do l_phys = 1,nr
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          do n = 1,n_solpts
            vartmp_at_sp(n) = jac(n)*metrics_array(n,l_comp,l_phys)
          end do
          !
          call interpolate_for_output( vartmp_at_sp , vartmp )
          !
          call cgns_write_field_data
          !
        end do
      end do
      !!}}}
    case ("TimeMetrics")
      !!{{{
      allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      if (allocated(time_metrics)) then
        !
        do l = 1,nr
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          do n = 1,n_solpts
            vartmp_at_sp(n) = one/time_metrics(l,n)
          end do
          !
          call interpolate_for_output( vartmp_at_sp , vartmp )
          !
          call cgns_write_field_data
          !
        end do
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,n_solpts
          vartmp_at_sp(n) = time_metrics(nr+1,n) * time_ref
        end do
        !
        call interpolate_for_output( vartmp_at_sp , vartmp )
        !
        call cgns_write_field_data
        !
      end if
      !!}}}
    case default
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,npts
        vartmp(n) = cv(n,nec) * rhoref
      end do
      !
      call cgns_write_field_data
      !!}}}
  end select variable_select
  !
  if (allocated(vartmp_at_sp)) then
    deallocate ( vartmp_at_sp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vartmp_at_sp",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(dpvdx)) then
    deallocate ( dpvdx , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"dpvdx",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(spcvar)) then
    deallocate ( spcvar , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"spcvar",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(other_variables)) then
    deallocate ( other_variables , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"other_variables",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  !-----------------------------------------------------------------------------
  !
  contains
  !
#include "Functions/cgns_write_field_data.f90"
  !
  !-----------------------------------------------------------------------------
  !
end subroutine output_variable
!
!###############################################################################
!
subroutine get_cgns_sol_variables(variable_list)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  use ovar,        only : itestcase,governing_equations
  use ovar,        only : limiter_option
  use ovar,        only : output_vorticity
  use ovar,        only : output_cylindrical_coordinates
  use ovar,        only : output_cylindrical_velocity
  use ovar,        only : output_relative_velocity
  use ovar,        only : output_relative_conditions
  use ovar,        only : output_qcriterion
  use ovar,        only : output_conserved_variables
  use ovar,        only : output_momentum
  use ovar,        only : output_temperature
  use ovar,        only : output_thermo_variables
  use ovar,        only : output_gradients
  use ovar,        only : output_transport_variables
  use ovar,        only : output_stagnation_temperature
  use ovar,        only : output_stagnation_pressure
  use ovar,        only : output_stagnation_variables
  use ovar,        only : output_marker_variables
  use ovar,        only : output_grid_details
  use ovar,        only : output_local_timestep
  use ovar,        only : Resolution_Marker_Variable
  use ovar,        only : Jump_Marker_Variable
  use ovar,        only : MLP_Marker_Variable
  use ovar,        only : output_species_density
  use ovar,        only : output_species_mass_fractions
  use ovar,        only : output_species_mole_fractions
  use ovar,        only : output_extra_species_variables
  use ovar,        only : output_debug_species_variables
  use ovar,        only : output_grdnt_species_variables
  use ovar,        only : Track_Limiter_Filter_Dissipation
  !
  !.. Formal Arguments ..
#ifdef SPECIAL_FOR_INTEL
  character(len=:), allocatable, intent(out) :: variable_list(:)
#else
  character(len=CGLEN), allocatable, intent(out) :: variable_list(:)
#endif
  !
  !.. Local Scalars ..
  integer :: i,j,n,m
  !
  character(len=CGLEN) :: varname
  character(len=CGLEN) :: varcheck
  !
  !.. Local Allocatable Arrays ..
#ifdef SPECIAL_FOR_INTEL
  character(len=:), allocatable :: temp_variables(:)
#else
  character(len=CGLEN), allocatable :: temp_variables(:)
#endif
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_cgns_sol_variables"
  !
  character(len=*), parameter :: default_list(1:4) = ["Density ", &
                                                      "Velocity", &
                                                      "Pressure", &
                                                      "Mach    "]
  !
  ! AVAILABLE QUANTITIES TO OUTPUT AND ORDER IN WHICH THEY WILL BE OUTPUT
  !
  character(len=*), parameter :: var_order(1:59) = ["Density             ", &
                                                    "Momentum            ", &
                                                    "Energy              ", &
                                                    "Velocity            ", &
                                                    "Pressure            ", &
                                                    "Temperature         ", &
                                                    "Mach                ", &
                                                    "MassFractions       ", &
                                                    "MoleFractions       ", &
                                                    "SpeciesDensities    ", &
                                                    "Enthalpy            ", &
                                                    "Cp                  ", &
                                                    "Entropy             ", &
                                                    "P0                  ", &
                                                    "T0                  ", &
                                                    "H0                  ", &
                                                    "Viscosity           ", &
                                                    "Conductivity        ", &
                                                    "Prandtl             ", &
                                                    "Gamma               ", &
                                                    "Rgas                ", &
                                                    "SpeciesEnthalpy     ", &
                                                    "SpeciesCp           ", &
                                                    "cpi/R               ", &
                                                    "hi/RT               ", &
                                                    "hi/R                ", &
                                                    "GradSpeciesDensities", &
                                                    "GradDensity         ", &
                                                    "GradMomentum        ", &
                                                    "GradEnergy          ", &
                                                    "GradMassFractions   ", &
                                                    "GradVelocity        ", &
                                                    "Vorticity           ", &
                                                    "CylindricalVelocity ", &
                                                    "RelativeMomentum    ", &
                                                    "RelativeVelocity    ", &
                                                    "RCylindricalVelocity", &
                                                    "RelativeMach        ", &
                                                    "GradTemperature     ", &
                                                    "ResolutionMarker    ", &
                                                    "JumpMarker          ", &
                                                    "DucrosMarker        ", &
                                                    "MLPMarker           ", &
                                                    "Timestep            ", &
                                                    "LocalTimestep       ", &
                                                    "Partition           ", &
                                                    "CellSize            ", &
                                                    "ContraMetrics       ", &
                                                    "TimeMetrics         ", &
                                                    "dTdx_Breakdown      ", &
                                                    "DmatTemperature     ", &
                                                    "CoMetrics           ", &
                                                    "GradTempMagnitude   ", &
                                                    "GradDensityMagnitude", &
                                                    "VorticityMagnitude  ", &
                                                    "StrainRateMagnitude ", &
                                                    "QCriterion          ", &
                                                    "DivergenceOfVelocity", &
                                                    "FilterDissipation   "]
  !
  !.. Variables requiring dusp for output
  character(len=*), parameter :: grad_vars(1:17) = ["GradSpeciesDensities", &
                                                    "GradDensity         ", &
                                                    "GradDensityMagnitude", &
                                                    "GradMomentum        ", &
                                                    "GradEnergy          ", &
                                                    "GradMassFractions   ", &
                                                    "GradVelocity        ", &
                                                    "GradTemperature     ", &
                                                    "GradTempMagnitude   ", &
                                                    "dTdx_Breakdown      ", &
                                                    "DmatTemperature     ", &
                                                    "Vorticity           ", &
                                                    "VorticityMagnitude  ", &
                                                    "StrainRateMagnitude ", &
                                                    "QCriterion          ", &
                                                    "DivergenceOfVelocity", &
                                                    "DucrosMarker        "]
  !
  ! Conserved Variables
  character(len=*), parameter :: consv_var(1:3) = ["Density ", &
                                                   "Momentum", &
                                                   "Energy  "]
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate variable_list to the default_list
  !
  call add_sol_variables(variable_list,default_list)
  !
  ! Make sure vorticity is one of the variables that will
  ! be output if running the Taylor-Green vortex problem.
  !
  if (any(itestcase == Taylor_Problems)) then
    call add_sol_variable(variable_list,"Vorticity")
    call add_sol_variable(variable_list,"VorticityMagnitude")
  end if
  !
  ! Output vorticity and other variables useful for identifying
  ! vortical structures
  !
  if (output_vorticity) then
    call add_sol_variable(variable_list,"Vorticity")
    call add_sol_variable(variable_list,"DivergenceOfVelocity")
    if (nr > 2) then
      call add_sol_variable(variable_list,"VorticityMagnitude")
      call add_sol_variable(variable_list,"StrainRateMagnitude")
      call add_sol_variable(variable_list,"QCriterion")
    end if
  end if
  !
  ! Output the cylindrical velocity components (vr,vtheta)
  if (output_cylindrical_coordinates .or. output_cylindrical_velocity) then
    if (nr >= 2) then
      call add_sol_variable(variable_list,"CylindricalVelocity")
      if (output_relative_conditions .or. output_relative_velocity) then
        call add_sol_variable(variable_list,"RCylindricalVelocity")
      end if
    end if
  end if
  !
  ! Output the relative velocities
  !
  if (output_relative_velocity) then
    call add_sol_variable(variable_list,"RelativeVelocity")
  else if (output_relative_conditions) then
    call add_sol_variable(variable_list,"RelativeVelocity")
    call add_sol_variable(variable_list,"RelativeMomentum")
    call add_sol_variable(variable_list,"RelativeMach")
  end if
  !
  ! Output Q-criterion without all the other vorticity variables
  !
  if (output_qcriterion) then
    if (nr > 2) then
      call add_sol_variable(variable_list,"QCriterion")
    end if
  end if
  !
  ! Make sure entropy is one of the variables that will be
  ! output if running any of the isentropic vortex problems.
  !
  if (any(itestcase == Transport_Problems)) then
    call add_sol_variable(variable_list,"Entropy")
  end if
  !
  if (output_conserved_variables) then
    call add_sol_variables(variable_list,consv_var)
    if (ncs > 1) then
      call add_sol_variable(variable_list,"SpeciesDensities")
    end if
  end if
  !
  if (output_momentum) then
    call add_sol_variable(variable_list,"Momentum")
  end if
  !
  if (output_temperature) then
    call add_sol_variable(variable_list,"Temperature")
  end if
  !
  if (output_thermo_variables) then
    call add_sol_variable(variable_list,"Pressure")
    call add_sol_variable(variable_list,"Temperature")
    call add_sol_variable(variable_list,"Enthalpy")
    call add_sol_variable(variable_list,"Cp")
    call add_sol_variable(variable_list,"Entropy")
  end if
  !
  if (output_stagnation_temperature) then
    call add_sol_variable(variable_list,"T0")
  end if
  !
  if (output_stagnation_pressure) then
    call add_sol_variable(variable_list,"P0")
  end if
  !
  if (output_stagnation_variables) then
    call add_sol_variable(variable_list,"P0")
    call add_sol_variable(variable_list,"T0")
    call add_sol_variable(variable_list,"H0")
  end if
  !
  if (output_local_timestep) then
    call add_sol_variable(variable_list,"LocalTimestep")
  end if
  !
  if (output_grid_details) then
    call add_sol_variable(variable_list,"Partition")
    call add_sol_variable(variable_list,"CellSize")
  end if
  !
  ! Output the marker quantities if requested
  !
  if (output_marker_variables) then
    if ( limiter_option == 3 ) then
      call add_sol_variable(variable_list,"MLPMarker")
      if (MLP_Marker_Variable > 0) then
        if (MLP_Marker_Variable == nec .or. ncs > 1) then
          m = 1
        else if (MLP_Marker_Variable == nee) then
          m = 3
        else
          m = 2
        end if
        call add_sol_variable(variable_list,trim(consv_var(m)))
      else
        call add_sol_variables(variable_list,consv_var)
      end if
    else if ( abs(limiter_option) == 2 ) then
      call add_sol_variable(variable_list,"JumpMarker")
      if (Jump_Marker_Variable > 0) then
        if (Jump_Marker_Variable == nec .or. ncs > 1) then
          m = 1
        else if (Jump_Marker_Variable == nee) then
          m = 3
        else
          m = 2
        end if
        call add_sol_variable(variable_list,trim(consv_var(m)))
      else
        call add_sol_variable(variable_list,'Pressure')
      end if
    else if ( abs(limiter_option) == 4 ) then
      call add_sol_variable(variable_list,"DucrosMarker")
      if (Jump_Marker_Variable > 0) then
        if (Jump_Marker_Variable == nec .or. ncs > 1) then
          m = 1
        else if (Jump_Marker_Variable == nee) then
          m = 3
        else
          m = 2
        end if
        call add_sol_variable(variable_list,trim(consv_var(m)))
      else
        call add_sol_variable(variable_list,'Pressure')
      end if
    else
      call add_sol_variable(variable_list,"ResolutionMarker")
      if (Resolution_Marker_Variable > 0) then
        if (Resolution_Marker_Variable == nec .or. ncs > 1) then
          m = 1
        else if (Resolution_Marker_Variable == nee) then
          m = 3
        else
          m = 2
        end if
        call add_sol_variable(variable_list,trim(consv_var(m)))
      else
        call add_sol_variables(variable_list,consv_var)
      end if
    end if
  end if
  !
  if (Track_Limiter_Filter_Dissipation) then
    call add_sol_variable(variable_list,"FilterDissipation")
  end if
  !
  ! Add species mass fractions to the variable list if using a gas mixture.
  ! Also output any other requested variables for gas mixtures.
  !
  if (ncs > 0) then
    !
    if (output_species_density) then
      call add_sol_variable(variable_list,"SpeciesDensities")
    end if
    !
    if (output_species_mass_fractions) then
      call add_sol_variable(variable_list,"MassFractions")
    end if
    !
    if (output_species_mole_fractions) then
      call add_sol_variable(variable_list,"MoleFractions")
    end if
    !
    if (output_extra_species_variables) then
      call add_sol_variable(variable_list,"SpeciesDensities")
      call add_sol_variable(variable_list,"Enthalpy")
      call add_sol_variable(variable_list,"H0")
      call add_sol_variable(variable_list,"Gamma")
      call add_sol_variable(variable_list,"Cp")
      call add_sol_variable(variable_list,"Rgas")
      call add_sol_variable(variable_list,"Momentum")
      call add_sol_variable(variable_list,"Energy")
      call add_sol_variable(variable_list,"P0")
      call add_sol_variable(variable_list,"T0")
    end if
    !
    if (output_grdnt_species_variables) then
      call add_sol_variable(variable_list,"GradMassFractions")
      call add_sol_variable(variable_list,"GradSpeciesDensities")
      call add_sol_variable(variable_list,"GradDensity")
      call add_sol_variable(variable_list,"GradMomentum")
      call add_sol_variable(variable_list,"GradEnergy")
      call add_sol_variable(variable_list,"GradVelocity")
      call add_sol_variable(variable_list,"GradTemperature")
    end if
    !
    if (output_debug_species_variables) then
      call add_sol_variable(variable_list,"SpeciesEnthalpy")
      call add_sol_variable(variable_list,"SpeciesCp")
      call add_sol_variable(variable_list,"hi/RT")
      call add_sol_variable(variable_list,"hi/R")
      call add_sol_variable(variable_list,"cpi/R")
      call add_sol_variable(variable_list,"dTdx_Breakdown")
    end if
    !
  end if
  !
  ! Output transport properties if solving the Navier-Stokes equations
  !
  if (output_transport_variables) then
    if (governing_equations == NavierStokes_Eqns) then
      !
      call add_sol_variable(variable_list,"Viscosity")
      call add_sol_variable(variable_list,"Conductivity")
      !
      if (ncs > 0) then
        if (output_extra_species_variables) then
          call add_sol_variable(variable_list,"Prandtl")
        end if
      end if
      !
    end if
  end if
  !
  ! If outputing gradients was requested, add the gradients of whatever
  ! variables are in the variable list.
  !
  if (output_gradients) then
    !
    allocate ( temp_variables , source=variable_list )
    !
    do i = 1,size(temp_variables)
      select case (uppercase(trim(adjustl(temp_variables(i)))))
        case ("DENSITY")
          call add_sol_variable(variable_list,"GradDensity")
        case ("VELOCITY")
          call add_sol_variable(variable_list,"GradVelocity")
        case ("ENERGY")
          call add_sol_variable(variable_list,"GradEnergy")
        case ("TEMPERATURE")
          call add_sol_variable(variable_list,"GradTemperature")
        case ("MOMENTUM")
          call add_sol_variable(variable_list,"GradMomentum")
        case ("MASSFRACTIONS")
          call add_sol_variable(variable_list,"GradMassFractions")
        case ("SPECIESDENSITIES")
          call add_sol_variable(variable_list,"GradSpeciesDensities")
      end select
    end do
    !
    deallocate ( temp_variables )
    !
  end if
  !
  ! Make sure there are no duplicate variables being output
  !
  allocate ( temp_variables , source=variable_list )
  !
  n = 0
  copy_loop: &
  do i = 1,size(variable_list)
    varcheck = uppercase( trim(adjustl(variable_list(i))) )
    do j = 1,n
      varname = uppercase( trim(adjustl(temp_variables(j))) )
      if (varname == varcheck) then
        ! This variable already exists in temp_variables, so dont add this
        ! to temp_variables and go on to the next variable in variable_list.
        cycle copy_loop
      end if
    end do
    n = n + 1
    temp_variables(n) = trim(adjustl(variable_list(i)))
  end do copy_loop
  !
  ! Use F2003 AUTO-REALLOCATION to cutoff the duplicate variables
  ! still contained in temp_variables(n+1:size(temp_variables)).
  !
  temp_variables = temp_variables(1:n)
  !
  ! Use move_alloc to over-write variable_list with temp_variables
  !
  call move_alloc( from=temp_variables , to=variable_list )
  !
  ! Put the output variables in order
  !
  allocate ( temp_variables , source=variable_list )
  !
  n = 0
  do i = 1,size(var_order)
    varname = uppercase( trim(adjustl(var_order(i))) )
    reorder_loop: &
    do j = 1,size(variable_list)
      varcheck = uppercase( trim(adjustl(variable_list(j))) )
      if (varname == varcheck) then
        n = n + 1
        temp_variables(n) = trim(adjustl(variable_list(j)))
        exit reorder_loop
      end if
    end do reorder_loop
  end do
  !
  ! Use move_alloc to over-write variable_list with temp_variables
  !
  call move_alloc( from=temp_variables , to=variable_list )
  !
  ! Check to see if any of the gradients are needed for an output variable.
  ! If not, we wont need to call interpolate_for_output on dusp.
  !
  interpolated_dusp_is_needed = fals
  dusp_check_loop: &
  do i = 1,size(variable_list)
    varname = uppercase( trim(adjustl(variable_list(i))) )
    do j = 1,size(grad_vars)
      varcheck = uppercase( trim(adjustl(grad_vars(j))) )
      if (varname == varcheck) then
        interpolated_dusp_is_needed = true
        exit dusp_check_loop
      end if
    end do
  end do dusp_check_loop
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_cgns_sol_variables
!
!###############################################################################
!
subroutine add_sol_variable(variable_list,new_variable)
  !
  !.. Formal Arguments ..
  character(len=*),                     intent(in) :: new_variable
#ifdef SPECIAL_FOR_INTEL
  character(len=:),     allocatable, intent(inout) :: variable_list(:)
#else
  character(len=CGLEN), allocatable, intent(inout) :: variable_list(:)
#endif
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "add_sol_variable"
  !
continue
  !
  if (.not. allocated(variable_list)) then
    !
#ifdef SPECIAL_FOR_INTEL
    n = len_trim(adjustl(new_variable))
    ! USING F2003 AUTO-REALLOCATION
    variable_list = [ character(len=n) :: trim(adjustl(new_variable)) ]
#else
    ! USING F2003 AUTO-REALLOCATION
    variable_list = [ adjustl(new_variable) ]
#endif
    !
  else
    !
#ifdef SPECIAL_FOR_INTEL
    n = max( len_trim(adjustl(new_variable)) , len(variable_list) )
    ! USING F2003 AUTO-REALLOCATION
    variable_list = [ character(len=n) :: variable_list , &
                                          trim(adjustl(new_variable)) ]
#else
    ! USING F2003 AUTO-REALLOCATION
    variable_list = [ variable_list , adjustl(new_variable) ]
#endif
    !
  end if
  !
end subroutine add_sol_variable
!
!###############################################################################
!
subroutine add_sol_variables(variable_list,new_variables)
  !
  !.. Formal Arguments ..
  character(len=*),                     intent(in) :: new_variables(:)
#ifdef SPECIAL_FOR_INTEL
  character(len=:),     allocatable, intent(inout) :: variable_list(:)
#else
  character(len=CGLEN), allocatable, intent(inout) :: variable_list(:)
#endif
  !
  !.. Local Scalars ..
  integer :: i,ibeg,n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "add_sol_variables"
  !
continue
  !
  if (.not. allocated(variable_list)) then
    !
#ifdef SPECIAL_FOR_INTEL
    n = len_trim(adjustl(new_variables(1)))
    ! USING F2003 AUTO-REALLOCATION
    variable_list = [ character(len=n) :: trim(adjustl(new_variables(1))) ]
#else
    ! USING F2003 AUTO-REALLOCATION
    variable_list = [ adjustl(new_variables(1)) ]
#endif
    !
    ibeg = 2
    !
  else
    !
    ibeg = 1
    !
  end if
  !
  do i = ibeg,size(new_variables)
#ifdef SPECIAL_FOR_INTEL
    n = max( len_trim(adjustl(new_variables(i))) , len(variable_list) )
    ! USING F2003 AUTO-REALLOCATION
    variable_list = [ character(len=n) :: variable_list , &
                                          trim(adjustl(new_variables(i))) ]
#else
    ! USING F2003 AUTO-REALLOCATION
    variable_list = [ variable_list , adjustl(new_variables(i)) ]
#endif
  end do
  !
end subroutine add_sol_variables
!
!###############################################################################
!
subroutine write_parallel_cgns_time_ave_file()
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  use ovar,        only : profile_io_cgns,profile_io_cgns_subs
  use ovar,        only : max_wall_time => max_time_ave_cgns_wall_time
  use ovar,        only : cgns_dir,itcur,itrst,num_timesteps
  use ovar,        only : time_ave_vel_is_axisymm
  use ovar,        only : time,time_ref,time_scaling_factor,ave_start_time
  use flowvar_mod, only : uavesp
  !
  !.. Local Scalars ..
  integer :: l,n,npts,ierr
  integer :: final_iter_digits
  !
  real(wp) :: beg_wall_time
  real(wp) :: end_wall_time
  real(wp) :: cgns_wall_time
  real(wp) :: cur_ave_time
  !
  character(len=100) :: char_num
  character(len=170) :: cgns_sol_file
  !
  integer(CBT) :: ifile_s,ibase_s,izone_s,isol_s
  !
  real(sp)             :: sol_time(1)
  integer              :: sol_iter(1)
  character(len=CGLEN) :: sol_name(1)
  !
  !.. Local Allocatable Arrays ..
  integer(CBT), allocatable :: isvar_s(:)
  !
  real(CRT), allocatable :: var(:)
  real(wp),  allocatable :: vartmp(:)
  !
#ifdef SPECIAL_FOR_INTEL
  character(len=:), allocatable :: char_sname(:)
#else
  character(len=CGLEN), allocatable :: char_sname(:)
#endif
  !
  !.. Local Saved Variables ..
  logical(lk),        save :: this_is_initial_time_ave_output = true
  character(len=170), save :: tafile_fmt
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_parallel_cgns_time_ave_file"
  character(len=*), parameter :: tname = "TimeAveragedSolution"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
#ifndef SPECIAL_FOR_PGI
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      write (iout,1) pname
      flush (iout)
      beg_wall_time = mpi_wtime()
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Create the format for the file name
  !
  if (this_is_initial_time_ave_output) then
    write (char_num,4) itrst + num_timesteps
    final_iter_digits = max(4,len_trim(adjustl(char_num)))
    write (tafile_fmt,3) trim(adjustl(cgns_dir)), final_iter_digits
    this_is_initial_time_ave_output = fals
  end if
  !
  ! Get the list of variables that will be written
  ! to the CGNS time averaged solution file
  !
  call get_cgns_time_ave_variables(char_sname)
  !
  ! Get the name of the current time averaged solution file
  !
  write (cgns_sol_file,tafile_fmt) itcur
  !
  ! Start the timer for writing this time averaged solution file
  !
  call cgns_timer(start_cgns_file_timer,profile_io_cgns_subs,cgns_sol_file)
  !
  ! Open the time averaged solution file.
  ! NOTE: This opens the CGNS file and creates the base and zone sections.
  !
  call open_cgns_file(cgns_sol_file,pname,ifile_s,ibase_s,izone_s)
  !
  ! Create the grid coordinates and element connectivity sections by
  ! using links to these respective sections within the CGNS grid file.
  !
  call link_cgns_grid_file(cgns_sol_file,pname,ifile_s,ibase_s,izone_s)
  !
  ! Create the flow solution section of the CGNS file
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_sol_write_f")
  call cg_sol_write_f(ifile_s,ibase_s,izone_s,tname,VERTEX,isol_s,cgierr)
  call cgns_error(pname,cgns_sol_file,"cg_sol_write_f",cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_sol_write_f")
  !
  !
  ! ##############################################################
  ! #####   WRITING TIME AVERAGED SOLUTION FIELD VARIABLES   #####
  ! ##############################################################
  !
  ! Count the total number of time averaged variables that will be written
  ! and allocate isvar_s accordingly
  !
  n = count_output_variables(char_sname)
  !
  isvar_s = init_array(n,0_CBT) ! F2003 AUTO-REALLOC
  !
  ! Create data nodes for the time averaged solution variables
  !
  call write_field_headers(cgns_sol_file,ifile_s,ibase_s,izone_s,isol_s, &
                           char_sname,isvar_s)
  !
  ! Preallocate var for writing variables in the correct precision
  ! to the CGNS file
  !
  npts = int( my_pts_end - my_pts_beg + 1_CST , kind=kind(npts) )
  !
  vartmp = init_array(npts,zero)          ! F2003 AUTO-REALLOC
  var = init_array(npts,real(0,kind=CRT)) ! F2003 AUTO-REALLOC
  !
  ! Create the output data and write it to the CGNS file in parallel
  !
  n = 0
  do l = 1,size(char_sname)
    call output_time_ave_variable(ifile_s,ibase_s,izone_s,isol_s,isvar_s, &
                                  cgns_sol_file,char_sname(l),n, &
                                  uavesp,vartmp,var)
  end do
  !
  ! Create a link to the grid coordinates in the grid file
  !
  if (time_ave_vel_is_axisymm .and. nr == 3) then
    call link_cylindrical_coordinates(cgns_sol_file,pname, &
                                      ifile_s,ibase_s,izone_s,isol_s)
  end if
  !
  ! ################################################
  ! #####   FINISHED WRITING FIELD VARIABLES   #####
  ! ################################################
  !
  ! Add the time information to the solution file
  !
  cur_ave_time = (time - ave_start_time)*time_ref*time_scaling_factor
  sol_time(1) = real(cur_ave_time,kind=kind(sol_time))
  sol_name(1) = tname
 !sol_name(1) = solname(size(solname))
  sol_iter(1) = itcur
  call write_iterative_data_to_cgns_file(cgns_sol_file,pname, &
                                         ifile_s,ibase_s,izone_s, &
                                         sol_time,sol_iter,sol_name)
  !
  ! Close the CGNS FILE
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cgp_close_f")
  call cgp_close_f(ifile_s,cgierr)
  call io_error(pname,cgns_sol_file,2,__LINE__,__FILE__,cgierr)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cgp_close_f")
  !
  call cgns_timer(stop_cgns_file_timer,profile_io_cgns_subs,cgns_sol_file)
  !
  ! Deallocate the arrays var and vartmp
  !
  if (allocated(var)) then
    deallocate ( var , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"var",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(vartmp)) then
    deallocate ( vartmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vartmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(isvar_s)) then
    deallocate ( isvar_s , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"isvar_s",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      end_wall_time = mpi_wtime()
      cgns_wall_time = end_wall_time - beg_wall_time
      max_wall_time = max(max_wall_time,cgns_wall_time)
      write (iout,2) pname,trim(adjustl(cgns_sol_file)),cgns_wall_time
      flush (iout)
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," Entering ",a)
  2 format (1x,"Leaving ",a,/,:, &
            5x,"Wall time taken to write the file '",a, &
            "' = ",f12.6," seconds",/)
  3 format ("('",a,"solution.time_ave.',i0.",i0,",'.cgns')")
  4 format (i0)
  !
end subroutine write_parallel_cgns_time_ave_file
!
!###############################################################################
!
subroutine get_cgns_time_ave_variables(variable_list)
  !
  !.. Use Statements ..
  use ovar,        only : time_ave_vel_is_axisymm
  use mrf_mod,     only : is_in_mrf
  use flowvar_mod, only : time_ave_variables
  !
  !.. Formal Arguments ..
#ifdef SPECIAL_FOR_INTEL
  character(len=:), allocatable, intent(out) :: variable_list(:)
#else
  character(len=CGLEN), allocatable, intent(out) :: variable_list(:)
#endif
  !
  !.. Local Scalars ..
  integer :: i,j,n,m
  !
  character(len=CGLEN) :: varname
  character(len=CGLEN) :: varcheck
  !
  !.. Local Allocatable Arrays ..
#ifdef SPECIAL_FOR_INTEL
  character(len=:), allocatable :: temp_variables(:)
  character(len=:), allocatable :: output_variables(:)
#else
  character(len=CGLEN), allocatable :: temp_variables(:)
  character(len=CGLEN), allocatable :: output_variables(:)
#endif
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_cgns_time_ave_variables"
  !
  character(len=*), parameter :: default_list(1:15) = ["Density    ", &
                                                       "Momentum   ", &
                                                       "Velocity   ", &
                                                       "Pressure   ", &
                                                       "Temperature", &
                                                       "uu         ", &
                                                       "vv         ", &
                                                       "ww         ", &
                                                       "uv         ", &
                                                       "uw         ", &
                                                       "vw         ", &
                                                       "TT         ", &
                                                       "uT         ", &
                                                       "vT         ", &
                                                       "wT         "]
  !
  character(len=*), parameter :: axi_list(1:15) = ["Density            ", &
                                                   "Velocity           ", &
                                                   "CylindricalVelocity", &
                                                   "Pressure           ", &
                                                   "Temperature        ", &
                                                   "uu                 ", &
                                                   "rr                 ", &
                                                   "ff                 ", &
                                                   "ur                 ", &
                                                   "uf                 ", &
                                                   "rf                 ", &
                                                   "TT                 ", &
                                                   "uT                 ", &
                                                   "rT                 ", &
                                                   "fT                 "]
  !
  ! AVAILABLE QUANTITIES TO OUTPUT AND ORDER IN WHICH THEY WILL BE OUTPUT
  !
  character(len=*), parameter :: var_order(1:33) = ["Density            ", &
                                                    "Momentum           ", &
                                                    "Velocity           ", &
                                                    "CylindricalVelocity", &
                                                    "Pressure           ", &
                                                    "Temperature        ", &
                                                    "P0                 ", &
                                                    "T0                 ", &
                                                    "uu                 ", &
                                                    "vv                 ", &
                                                    "ww                 ", &
                                                    "rr                 ", &
                                                    "ff                 ", &
                                                    "uv                 ", &
                                                    "uw                 ", &
                                                    "vw                 ", &
                                                    "ur                 ", &
                                                    "uf                 ", &
                                                    "rf                 ", &
                                                    "TT                 ", &
                                                    "uT                 ", &
                                                    "vT                 ", &
                                                    "wT                 ", &
                                                    "rT                 ", &
                                                    "fT                 ", &
                                                    "dd                 ", &
                                                    "du                 ", &
                                                    "dv                 ", &
                                                    "dw                 ", &
                                                    "dr                 ", &
                                                    "df                 ", &
                                                    "dT                 ", &
                                                    "PP                 "]
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate variable_list to the default_list
  !
  if (time_ave_vel_is_axisymm) then
    call add_sol_variables(variable_list,axi_list)
  else
    call add_sol_variables(variable_list,default_list)
  end if
 !!
 !! Add total temperature and pressure if in a moving reference frame
 !!
 !if ( is_in_mrf() ) then
 !  call add_sol_variable(variable_list,"T0")
 !  call add_sol_variable(variable_list,"P0")
 !end if
  !
  ! Add total temperature if it exists in the time averaged variables
  !
  t0_check: &
  do m = 1,size(time_ave_variables)
    if (trim(adjustl( time_ave_variables(m) )) == "T") then
      call add_sol_variable(variable_list,"T0")
      exit t0_check
    end if
  end do t0_check
  !
  ! Add total pressure if it exists in the time averaged variables
  !
  p0_check: &
  do m = 1,size(time_ave_variables)
    if (trim(adjustl( time_ave_variables(m) )) == "P") then
      call add_sol_variable(variable_list,"P0")
      exit p0_check
    end if
  end do p0_check
  !
  ! Make sure there are no duplicate variables being output
  !
  allocate ( temp_variables , source=variable_list )
  !
  n = 0
  copy_loop: &
  do i = 1,size(variable_list)
    varcheck = uppercase( trim(adjustl(variable_list(i))) )
    do j = 1,n
      varname = uppercase( trim(adjustl(temp_variables(j))) )
      if (varname == varcheck) then
        ! This variable already exists in temp_variables, so dont add this
        ! to temp_variables and go on to the next variable in variable_list.
        cycle copy_loop
      end if
    end do
    n = n + 1
    temp_variables(n) = trim(adjustl(variable_list(i)))
  end do copy_loop
  !
  ! Use F2003 AUTO-REALLOCATION to cutoff the duplicate variables
  ! still contained in temp_variables(n+1:size(temp_variables)).
  !
  temp_variables = temp_variables(1:n)
  !
  ! Use move_alloc to over-write variable_list with temp_variables
  !
  call move_alloc( from=temp_variables , to=variable_list )
  !
  ! Put the output variables in order
  !
  allocate ( temp_variables , source=variable_list )
  !
  n = 0
  do i = 1,size(var_order)
    varname = uppercase( trim(adjustl(var_order(i))) )
    reorder_loop: &
    do j = 1,size(variable_list)
      varcheck = uppercase( trim(adjustl(variable_list(j))) )
      if (varname == varcheck) then
        n = n + 1
        temp_variables(n) = trim(adjustl(variable_list(j)))
        exit reorder_loop
      end if
    end do reorder_loop
  end do
  !
  ! Use move_alloc to over-write variable_list with temp_variables
  !
  call move_alloc( from=temp_variables , to=variable_list )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_cgns_time_ave_variables
!
!###############################################################################
!
subroutine output_time_ave_variable(ifile_s,ibase_s,izone_s,isol_s,isvar_s, &
                                    cgns_sol_file,vname,nv, &
                                    uavesp,vartmp,var)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  use ovar,        only : aref,rhoref,tref
  use ovar,        only : profile_io_cgns_subs
  use ovar,        only : interpolate_before_output_variable
  use flowvar_mod, only : time_ave_variables
  !
  !.. Formal Arguments ..
  integer(CBT),     intent(in) :: ifile_s
  integer(CBT),     intent(in) :: ibase_s
  integer(CBT),     intent(in) :: izone_s
  integer(CBT),     intent(in) :: isol_s
  integer(CBT),     intent(in) :: isvar_s(:)
  character(len=*), intent(in) :: cgns_sol_file
  character(len=*), intent(in) :: vname
  integer,       intent(inout) :: nv
  real(wp),         intent(in) :: uavesp(:,:)
  real(wp),      intent(inout) :: vartmp(:)
  real(CRT),     intent(inout) :: var(:)
  !
  !.. Local Scalars ..
  integer :: l,m,n,npts,ierr
  integer :: idx1,idx2,idx3
  !
  integer(CBT) :: ivar_s
  !
  character(len=1) :: char_axi
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: vartmp_at_sp(:)
  real(wp), allocatable :: v1(:)
  real(wp), allocatable :: v2(:)
  real(wp), allocatable :: v1v2(:)
  !
  character(len=:), allocatable :: ta_var
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_time_ave_variable"
  !
continue
  !
  npts = size(vartmp)
  !
  idx1 = 0
  idx2 = 0
  idx3 = 0
  !
  ! Initialize vartmp to zero in case variable indices were not found
  !
  vartmp(:) = zero
  !
  ! Compute the current time averaged output variable
  ! and write it to the CGNS file
  !
  select case (trim(adjustl(vname)))
    !
    case ("Density") ! Density
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "d") then
          idx1 = m
          exit
        end if
      end do
      !
      if (idx1 > 0) then
        call interpolate_for_output( uavesp(:,idx1) , vartmp )
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * rhoref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("Velocity") ! Mean Velocity
      !!{{{
      do l = 1,nr
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do m = 1,size(time_ave_variables)
          ! F2003 AUTO-REALLOC
          ta_var = trim(adjustl( time_ave_variables(m) ))
          if (ta_var == char_lc_xyz(l)) then
            idx1 = m
            exit
          end if
        end do
        !
        if (idx1 > 0) then
          call interpolate_for_output( uavesp(:,idx1) , vartmp )
        end if
        !
        do n = 1,npts
          vartmp(n) = vartmp(n) * aref
        end do
        !
        call cgns_write_field_data
        !
        idx1 = 0
        !
      end do
      !!}}}
    case ("Momentum") ! Mean momentum
      !!{{{
      do l = 1,nr
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do m = 1,size(time_ave_variables)
          ! F2003 AUTO-REALLOC
          ta_var = trim(adjustl( time_ave_variables(m) ))
          if (ta_var == "d"//char_lc_xyz(l)) then
            idx1 = m
            exit
          end if
        end do
        !
        if (idx1 > 0) then
          call interpolate_for_output( uavesp(:,idx1) , vartmp )
        end if
        !
        do n = 1,npts
          vartmp(n) = vartmp(n) * rhoref * aref
        end do
        !
        call cgns_write_field_data
        !
        idx1 = 0
        !
      end do
      !!}}}
    case ("CylindricalVelocity") ! Cylindrical Velocities (Vr and Vtheta)
      !!{{{
      do l = 1,2
        !
        char_axi = merge("r","f",l==1)
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do m = 1,size(time_ave_variables)
          ! F2003 AUTO-REALLOC
          ta_var = trim(adjustl( time_ave_variables(m) ))
          if (ta_var == char_axi) then
            idx1 = m
            exit
          end if
        end do
        !
        if (idx1 > 0) then
          call interpolate_for_output( uavesp(:,idx1) , vartmp )
        end if
        !
        do n = 1,npts
          vartmp(n) = vartmp(n) * aref
        end do
        !
        call cgns_write_field_data
        !
        idx1 = 0
        !
      end do
      !!}}}
    case ("Pressure") ! Pressure
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "p") then
          idx1 = m
          exit
        end if
      end do
      !
      if (idx1 > 0) then
        call interpolate_for_output( uavesp(:,idx1) , vartmp )
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * rhoref*aref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("Temperature") ! Temperature
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "t") then
          idx1 = m
          exit
        end if
      end do
      !
      if (idx1 > 0) then
        call interpolate_for_output( uavesp(:,idx1) , vartmp )
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * tref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("P0") ! Pressure
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "P") then
          idx1 = m
          exit
        end if
      end do
      !
      if (idx1 > 0) then
        call interpolate_for_output( uavesp(:,idx1) , vartmp )
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * rhoref*aref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("T0") ! Temperature
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "T") then
          idx1 = m
          exit
        end if
      end do
      !
      if (idx1 > 0) then
        call interpolate_for_output( uavesp(:,idx1) , vartmp )
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * tref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("uu") ! axial velocity fluctuation
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "x") then
          idx1 = m
        else if (ta_var == "xx") then
          idx2 = m
        end if
        if (all([idx1,idx2] > 0)) exit
      end do
      !
      if (all([idx1,idx2] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx2) - uavesp(:,idx1)*uavesp(:,idx1)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          do n = 1,npts
            vartmp(n) = v2(n) - v1(n)*v1(n)
          end do
          deallocate ( v1 , v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * aref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("vv") ! vertical velocity fluctuation
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "y") then
          idx1 = m
        else if (ta_var == "yy") then
          idx2 = m
        end if
        if (all([idx1,idx2] > 0)) exit
      end do
      !
      if (all([idx1,idx2] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx2) - uavesp(:,idx1)*uavesp(:,idx1)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          do n = 1,npts
            vartmp(n) = v2(n) - v1(n)*v1(n)
          end do
          deallocate ( v1 , v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * aref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("ww") ! spanwise velocity fluctuation
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "z") then
          idx1 = m
        else if (ta_var == "zz") then
          idx2 = m
        end if
        if (all([idx1,idx2] > 0)) exit
      end do
      !
      if (all([idx1,idx2] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx2) - uavesp(:,idx1)*uavesp(:,idx1)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          do n = 1,npts
            vartmp(n) = v2(n) - v1(n)*v1(n)
          end do
          deallocate ( v1 , v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * aref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("rr") ! radial velocity fluctuation
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "r") then
          idx1 = m
        else if (ta_var == "rr") then
          idx2 = m
        end if
        if (all([idx1,idx2] > 0)) exit
      end do
      !
      if (all([idx1,idx2] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx2) - uavesp(:,idx1)*uavesp(:,idx1)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          do n = 1,npts
            vartmp(n) = v2(n) - v1(n)*v1(n)
          end do
          deallocate ( v1 , v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * aref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("ff") ! azimuthal velocity fluctuation
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "f") then
          idx1 = m
        else if (ta_var == "ff") then
          idx2 = m
        end if
        if (all([idx1,idx2] > 0)) exit
      end do
      !
      if (all([idx1,idx2] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx2) - uavesp(:,idx1)*uavesp(:,idx1)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          do n = 1,npts
            vartmp(n) = v2(n) - v1(n)*v1(n)
          end do
          deallocate ( v1 , v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * aref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("uv") ! X-Y shear stress
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "x") then
          idx1 = m
        else if (ta_var == "y") then
          idx2 = m
        else if (ta_var == "xy") then
          idx3 = m
        end if
        if (all([idx1,idx2,idx3] > 0)) exit
      end do
      !
      if (all([idx1,idx2,idx3] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx3) - uavesp(:,idx1)*uavesp(:,idx2)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v1v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          call interpolate_for_output( uavesp(:,idx3) , v1v2 )
          do n = 1,npts
            vartmp(n) = v1v2(n) - v1(n)*v2(n)
          end do
          deallocate ( v1 , v2 , v1v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * aref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("uw") ! X-Z shear stress
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "x") then
          idx1 = m
        else if (ta_var == "z") then
          idx2 = m
        else if (ta_var == "xz") then
          idx3 = m
        end if
        if (all([idx1,idx2,idx3] > 0)) exit
      end do
      !
      if (all([idx1,idx2,idx3] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx3) - uavesp(:,idx1)*uavesp(:,idx2)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v1v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          call interpolate_for_output( uavesp(:,idx3) , v1v2 )
          do n = 1,npts
            vartmp(n) = v1v2(n) - v1(n)*v2(n)
          end do
          deallocate ( v1 , v2 , v1v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * aref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("vw") ! Y-Z shear stress
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "y") then
          idx1 = m
        else if (ta_var == "z") then
          idx2 = m
        else if (ta_var == "yz") then
          idx3 = m
        end if
        if (all([idx1,idx2,idx3] > 0)) exit
      end do
      !
      if (all([idx1,idx2,idx3] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx3) - uavesp(:,idx1)*uavesp(:,idx2)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v1v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          call interpolate_for_output( uavesp(:,idx3) , v1v2 )
          do n = 1,npts
            vartmp(n) = v1v2(n) - v1(n)*v2(n)
          end do
          deallocate ( v1 , v2 , v1v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * aref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("ur") ! X-R shear stress
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "x") then
          idx1 = m
        else if (ta_var == "r") then
          idx2 = m
        else if (ta_var == "xr") then
          idx3 = m
        end if
        if (all([idx1,idx2,idx3] > 0)) exit
      end do
      !
      if (all([idx1,idx2,idx3] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx3) - uavesp(:,idx1)*uavesp(:,idx2)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v1v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          call interpolate_for_output( uavesp(:,idx3) , v1v2 )
          do n = 1,npts
            vartmp(n) = v1v2(n) - v1(n)*v2(n)
          end do
          deallocate ( v1 , v2 , v1v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * aref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("uf") ! X-Theta shear stress
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "x") then
          idx1 = m
        else if (ta_var == "f") then
          idx2 = m
        else if (ta_var == "xf") then
          idx3 = m
        end if
        if (all([idx1,idx2,idx3] > 0)) exit
      end do
      !
      if (all([idx1,idx2,idx3] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx3) - uavesp(:,idx1)*uavesp(:,idx2)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v1v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          call interpolate_for_output( uavesp(:,idx3) , v1v2 )
          do n = 1,npts
            vartmp(n) = v1v2(n) - v1(n)*v2(n)
          end do
          deallocate ( v1 , v2 , v1v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * aref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("rf") ! R-Theta shear stress
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "r") then
          idx1 = m
        else if (ta_var == "f") then
          idx2 = m
        else if (ta_var == "rf") then
          idx3 = m
        end if
        if (all([idx1,idx2,idx3] > 0)) exit
      end do
      !
      if (all([idx1,idx2,idx3] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx3) - uavesp(:,idx1)*uavesp(:,idx2)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v1v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          call interpolate_for_output( uavesp(:,idx3) , v1v2 )
          do n = 1,npts
            vartmp(n) = v1v2(n) - v1(n)*v2(n)
          end do
          deallocate ( v1 , v2 , v1v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * aref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("TT") ! Temperature Fluctuation
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "t") then
          idx1 = m
        else if (ta_var == "tt") then
          idx2 = m
        end if
        if (all([idx1,idx2] > 0)) exit
      end do
      !
      if (all([idx1,idx2] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx2) - uavesp(:,idx1)*uavesp(:,idx1)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          do n = 1,npts
            vartmp(n) = v2(n) - v1(n)*v1(n)
          end do
          deallocate ( v1 , v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * tref*tref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("uT") ! Streamwise heat flux
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "x") then
          idx1 = m
        else if (ta_var == "t") then
          idx2 = m
        else if (ta_var == "xt") then
          idx3 = m
        end if
        if (all([idx1,idx2,idx3] > 0)) exit
      end do
      !
      if (all([idx1,idx2,idx3] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx3) - uavesp(:,idx1)*uavesp(:,idx2)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v1v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          call interpolate_for_output( uavesp(:,idx3) , v1v2 )
          do n = 1,npts
            vartmp(n) = v1v2(n) - v1(n)*v2(n)
          end do
          deallocate ( v1 , v2 , v1v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * aref*tref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("vT") ! Vertical heat flux
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "y") then
          idx1 = m
        else if (ta_var == "t") then
          idx2 = m
        else if (ta_var == "yt") then
          idx3 = m
        end if
        if (all([idx1,idx2,idx3] > 0)) exit
      end do
      !
      if (all([idx1,idx2,idx3] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx3) - uavesp(:,idx1)*uavesp(:,idx2)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v1v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          call interpolate_for_output( uavesp(:,idx3) , v1v2 )
          do n = 1,npts
            vartmp(n) = v1v2(n) - v1(n)*v2(n)
          end do
          deallocate ( v1 , v2 , v1v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * aref*tref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("wT") ! Spanwise heat flux
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "z") then
          idx1 = m
        else if (ta_var == "t") then
          idx2 = m
        else if (ta_var == "zt") then
          idx3 = m
        end if
        if (all([idx1,idx2,idx3] > 0)) exit
      end do
      !
      if (all([idx1,idx2,idx3] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx3) - uavesp(:,idx1)*uavesp(:,idx2)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v1v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          call interpolate_for_output( uavesp(:,idx3) , v1v2 )
          do n = 1,npts
            vartmp(n) = v1v2(n) - v1(n)*v2(n)
          end do
          deallocate ( v1 , v2 , v1v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * aref*tref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("rT") ! Radial heat flux
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "r") then
          idx1 = m
        else if (ta_var == "t") then
          idx2 = m
        else if (ta_var == "rt") then
          idx3 = m
        end if
        if (all([idx1,idx2,idx3] > 0)) exit
      end do
      !
      if (all([idx1,idx2,idx3] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx3) - uavesp(:,idx1)*uavesp(:,idx2)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v1v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          call interpolate_for_output( uavesp(:,idx3) , v1v2 )
          do n = 1,npts
            vartmp(n) = v1v2(n) - v1(n)*v2(n)
          end do
          deallocate ( v1 , v2 , v1v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * aref*tref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("fT") ! Azimuthal heat flux
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "f") then
          idx1 = m
        else if (ta_var == "t") then
          idx2 = m
        else if (ta_var == "ft") then
          idx3 = m
        end if
        if (all([idx1,idx2,idx3] > 0)) exit
      end do
      !
      if (all([idx1,idx2,idx3] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx3) - uavesp(:,idx1)*uavesp(:,idx2)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v1v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          call interpolate_for_output( uavesp(:,idx3) , v1v2 )
          do n = 1,npts
            vartmp(n) = v1v2(n) - v1(n)*v2(n)
          end do
          deallocate ( v1 , v2 , v1v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * aref*tref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("dd") ! Density Fluctuation
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "d") then
          idx1 = m
        else if (ta_var == "dd") then
          idx2 = m
        end if
        if (all([idx1,idx2] > 0)) exit
      end do
      !
      if (all([idx1,idx2] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx2) - uavesp(:,idx1)*uavesp(:,idx1)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          do n = 1,npts
            vartmp(n) = v2(n) - v1(n)*v1(n)
          end do
          deallocate ( v1 , v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * rhoref*rhoref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("du") ! Axial momentum fluctuation
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "d") then
          idx1 = m
        else if (ta_var == "x") then
          idx2 = m
        else if (ta_var == "dx") then
          idx3 = m
        end if
        if (all([idx1,idx2,idx3] > 0)) exit
      end do
      !
      if (all([idx1,idx2,idx3] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx3) - uavesp(:,idx1)*uavesp(:,idx2)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v1v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          call interpolate_for_output( uavesp(:,idx3) , v1v2 )
          do n = 1,npts
            vartmp(n) = v1v2(n) - v1(n)*v2(n)
          end do
          deallocate ( v1 , v2 , v1v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * rhoref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("dv") ! Vertical momentum fluctuation
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "d") then
          idx1 = m
        else if (ta_var == "y") then
          idx2 = m
        else if (ta_var == "dy") then
          idx3 = m
        end if
        if (all([idx1,idx2,idx3] > 0)) exit
      end do
      !
      if (all([idx1,idx2,idx3] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx3) - uavesp(:,idx1)*uavesp(:,idx2)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v1v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          call interpolate_for_output( uavesp(:,idx3) , v1v2 )
          do n = 1,npts
            vartmp(n) = v1v2(n) - v1(n)*v2(n)
          end do
          deallocate ( v1 , v2 , v1v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * rhoref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("dw") ! Spanwise momentum fluctuation
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "d") then
          idx1 = m
        else if (ta_var == "z") then
          idx2 = m
        else if (ta_var == "dz") then
          idx3 = m
        end if
        if (all([idx1,idx2,idx3] > 0)) exit
      end do
      !
      if (all([idx1,idx2,idx3] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx3) - uavesp(:,idx1)*uavesp(:,idx2)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v1v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          call interpolate_for_output( uavesp(:,idx3) , v1v2 )
          do n = 1,npts
            vartmp(n) = v1v2(n) - v1(n)*v2(n)
          end do
          deallocate ( v1 , v2 , v1v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * rhoref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("dr") ! Radial momentum fluctuation
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "d") then
          idx1 = m
        else if (ta_var == "r") then
          idx2 = m
        else if (ta_var == "dr") then
          idx3 = m
        end if
        if (all([idx1,idx2,idx3] > 0)) exit
      end do
      !
      if (all([idx1,idx2,idx3] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx3) - uavesp(:,idx1)*uavesp(:,idx2)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v1v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          call interpolate_for_output( uavesp(:,idx3) , v1v2 )
          do n = 1,npts
            vartmp(n) = v1v2(n) - v1(n)*v2(n)
          end do
          deallocate ( v1 , v2 , v1v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * rhoref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("df") ! Azimuthal momentum fluctuation
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "d") then
          idx1 = m
        else if (ta_var == "f") then
          idx2 = m
        else if (ta_var == "df") then
          idx3 = m
        end if
        if (all([idx1,idx2,idx3] > 0)) exit
      end do
      !
      if (all([idx1,idx2,idx3] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx3) - uavesp(:,idx1)*uavesp(:,idx2)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v1v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          call interpolate_for_output( uavesp(:,idx3) , v1v2 )
          do n = 1,npts
            vartmp(n) = v1v2(n) - v1(n)*v2(n)
          end do
          deallocate ( v1 , v2 , v1v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * rhoref*aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("dT") ! Density-temperature fluctuation
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "d") then
          idx1 = m
        else if (ta_var == "t") then
          idx2 = m
        else if (ta_var == "dt") then
          idx3 = m
        end if
        if (all([idx1,idx2,idx3] > 0)) exit
      end do
      !
      if (all([idx1,idx2,idx3] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx3) - uavesp(:,idx1)*uavesp(:,idx2)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v1v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          call interpolate_for_output( uavesp(:,idx3) , v1v2 )
          do n = 1,npts
            vartmp(n) = v1v2(n) - v1(n)*v2(n)
          end do
          deallocate ( v1 , v2 , v1v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * rhoref*tref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("PP") ! Pressure Fluctuation
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do m = 1,size(time_ave_variables)
        ! F2003 AUTO-REALLOC
        ta_var = trim(adjustl( time_ave_variables(m) ))
        if (ta_var == "p") then
          idx1 = m
        else if (ta_var == "pp") then
          idx2 = m
        end if
        if (all([idx1,idx2] > 0)) exit
      end do
      !
      if (all([idx1,idx2] > 0)) then
        if (interpolate_before_output_variable) then
          ! F2003 AUTO-REALLOC
          vartmp_at_sp = uavesp(:,idx2) - uavesp(:,idx1)*uavesp(:,idx1)
          call interpolate_for_output( vartmp_at_sp , vartmp )
          deallocate ( vartmp_at_sp , stat=ierr )
        else
          v1 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          v2 = init_array(npts,zero) ! F2003 AUTO-REALLOC
          call interpolate_for_output( uavesp(:,idx1) , v1 )
          call interpolate_for_output( uavesp(:,idx2) , v2 )
          do n = 1,npts
            vartmp(n) = v2(n) - v1(n)*v1(n)
          end do
          deallocate ( v1 , v2 , stat=ierr )
        end if
      end if
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * (rhoref*aref*aref)**2
      end do
      !
      call cgns_write_field_data
      !!}}}
    case default
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      idx1 = 1
      call interpolate_for_output( uavesp(:,idx1) , vartmp )
      !
      do n = 1,npts
        vartmp(n) = vartmp(n) * rhoref
      end do
      !
      call cgns_write_field_data
      !!}}}
  end select
  !
  if (allocated(vartmp_at_sp)) then
    deallocate ( vartmp_at_sp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vartmp_at_sp",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(v1)) then
    deallocate ( v1 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"v1",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(v2)) then
    deallocate ( v2 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"v2",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(v1v2)) then
    deallocate ( v1v2 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"v1v2",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  !-----------------------------------------------------------------------------
  !
  contains
  !
#include "Functions/cgns_write_field_data.f90"
  !
  !-----------------------------------------------------------------------------
  !
end subroutine output_time_ave_variable
!
!###############################################################################
!
subroutine write_parallel_cgns_timeave_data_file()
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  use ovar,        only : profile_io_cgns,profile_io_cgns_subs
  use ovar,        only : max_wall_time => max_time_ave_cgns_wall_time
  use ovar,        only : cgns_dir,itcur,itrst,num_timesteps
  use ovar,        only : time_ave_vel_is_axisymm
  use ovar,        only : time,time_ref,time_scaling_factor,ave_start_time
  use flowvar_mod, only : uavesp
  !
  !.. Local Scalars ..
  integer :: l,n,nv,nvar,npts,ierr
  integer :: final_iter_digits
  !
  real(wp) :: beg_wall_time
  real(wp) :: end_wall_time
  real(wp) :: cgns_wall_time
  real(wp) :: cur_ave_time
  !
  character(len=100) :: char_num
  character(len=170) :: cgns_sol_file
  !
  character(len=CGLEN) :: varname
  !
  integer(CBT) :: ifile_s,ibase_s,izone_s,isol_s,ivar_s
  integer(CBT) :: ifile,ibase,izone,isol
  !
  real(sp)             :: sol_time(1)
  integer              :: sol_iter(1)
  character(len=CGLEN) :: sol_name(1)
  !
  !.. Local Allocatable Arrays ..
  integer(CBT), allocatable :: isvar_s(:)
  !
  real(CRT), allocatable :: var(:)
  real(wp),  allocatable :: vartmp(:)
  real(wp),  allocatable :: variable_ref(:)
  !
#ifdef SPECIAL_FOR_INTEL
  character(len=:), allocatable :: char_sname(:)
#else
  character(len=CGLEN), allocatable :: char_sname(:)
#endif
  !
  !.. Local Saved Variables ..
  logical(lk),        save :: this_is_initial_time_ave_output = true
  character(len=170), save :: tafile_fmt
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_parallel_cgns_timeave_data_file"
  character(len=*), parameter :: tname = "TimeAveragedData"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
#ifndef SPECIAL_FOR_PGI
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      write (iout,1) pname
      flush (iout)
      beg_wall_time = mpi_wtime()
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Create the format for the file name
  !
  if (this_is_initial_time_ave_output) then
    write (char_num,4) itrst + num_timesteps
    final_iter_digits = max(4,len_trim(adjustl(char_num)))
    write (tafile_fmt,3) trim(adjustl(cgns_dir)), final_iter_digits
    this_is_initial_time_ave_output = fals
  end if
  !
  ! Get the list of variables that will be written
  ! to the CGNS time averaged solution file
  !
  call create_time_ave_data_names(char_sname,variable_ref)
  !
  ! Get the name of the current time averaged solution file
  !
  write (cgns_sol_file,tafile_fmt) itcur
  !
  ! Start the timer for writing this time averaged solution file
  !
  call cgns_timer(start_cgns_file_timer,profile_io_cgns_subs,cgns_sol_file)
  !
  ! Open the time averaged solution file.
  ! NOTE: This opens the CGNS file and creates the base and zone sections.
  !
  call open_cgns_file(cgns_sol_file,pname,ifile_s,ibase_s,izone_s)
  !
  ! Create the grid coordinates and element connectivity sections by
  ! using links to these respective sections within the CGNS grid file.
  !
  call link_cgns_grid_file(cgns_sol_file,pname,ifile_s,ibase_s,izone_s)
  !
  ! Create the flow solution section of the CGNS file
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_sol_write_f")
  call cg_sol_write_f(ifile_s,ibase_s,izone_s,tname,VERTEX,isol_s,cgierr)
  call cgns_error(pname,cgns_sol_file,"cg_sol_write_f",cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_sol_write_f")
  !
  !
  ! ##############################################################
  ! #####   WRITING TIME AVERAGED SOLUTION FIELD VARIABLES   #####
  ! ##############################################################
  !
  ! Count the total number of time averaged variables that will be written
  ! and allocate isvar_s accordingly
  !
  nvar = count_output_variables(char_sname)
  !
  isvar_s = init_array(nvar,0_CBT) ! F2003 AUTO-REALLOC
  !
  ! Create data nodes for the time averaged solution variables
  !
  do nv = 1,nvar
    varname = trim(adjustl(char_sname(nv)))
    ifile = ifile_s
    ibase = ibase_s
    izone = izone_s
    isol = isol_s
    call cgns_write_field_header(nv,varname,isvar_s(nv),__LINE__,__FILE__)
  end do
 !call write_field_headers(cgns_sol_file,ifile_s,ibase_s,izone_s,isol_s, &
 !                         char_sname,isvar_s)
  !
  ! Preallocate var for writing variables in the correct precision
  ! to the CGNS file
  !
  npts = int( my_pts_end - my_pts_beg + 1_CST , kind=kind(npts) )
  !
  vartmp = init_array(npts,zero)          ! F2003 AUTO-REALLOC
  var = init_array(npts,real(0,kind=CRT)) ! F2003 AUTO-REALLOC
  !
  ! Create the output data and write it to the CGNS file in parallel
  !
  do nv = 1,nvar
    !
    ivar_s = isvar_s(nv)
    !
    call interpolate_for_output( uavesp(:,nv) , vartmp )
    !
    do n = 1,npts
      vartmp(n) = vartmp(n) * variable_ref(nv)
    end do
    !
    call cgns_write_field_data
    !
  end do
 !n = 0
 !do l = 1,size(char_sname)
 !  call output_time_ave_variable(ifile_s,ibase_s,izone_s,isol_s,isvar_s, &
 !                                cgns_sol_file,char_sname(l),n, &
 !                                uavesp,vartmp,var)
 !end do
  !
  ! Create a link to the grid coordinates in the grid file
  !
  if (time_ave_vel_is_axisymm .and. nr == 3) then
    call link_cylindrical_coordinates(cgns_sol_file,pname, &
                                      ifile_s,ibase_s,izone_s,isol_s)
  end if
  !
  ! ################################################
  ! #####   FINISHED WRITING FIELD VARIABLES   #####
  ! ################################################
  !
  ! Add the time information to the solution file
  !
  cur_ave_time = (time - ave_start_time)*time_ref*time_scaling_factor
  sol_time(1) = real(cur_ave_time,kind=kind(sol_time))
  sol_name(1) = tname
 !sol_name(1) = solname(size(solname))
  sol_iter(1) = itcur
  call write_iterative_data_to_cgns_file(cgns_sol_file,pname, &
                                         ifile_s,ibase_s,izone_s, &
                                         sol_time,sol_iter,sol_name)
  !
  ! Close the CGNS FILE
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cgp_close_f")
  call cgp_close_f(ifile_s,cgierr)
  call io_error(pname,cgns_sol_file,2,__LINE__,__FILE__,cgierr)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cgp_close_f")
  !
  call cgns_timer(stop_cgns_file_timer,profile_io_cgns_subs,cgns_sol_file)
  !
  ! Deallocate the arrays var and vartmp
  !
  if (allocated(var)) then
    deallocate ( var , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"var",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(vartmp)) then
    deallocate ( vartmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vartmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(isvar_s)) then
    deallocate ( isvar_s , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"isvar_s",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      end_wall_time = mpi_wtime()
      cgns_wall_time = end_wall_time - beg_wall_time
      max_wall_time = max(max_wall_time,cgns_wall_time)
      write (iout,2) pname,trim(adjustl(cgns_sol_file)),cgns_wall_time
      flush (iout)
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," Entering ",a)
  2 format (1x,"Leaving ",a,/,:, &
            5x,"Wall time taken to write the file '",a, &
            "' = ",f12.6," seconds",/)
  3 format ("('",a,"solution.ta_data.',i0.",i0,",'.cgns')")
  4 format (i0)
  !
  !-----------------------------------------------------------------------------
  !
  contains
  !
#include "Functions/cgns_write_field_header.f90"
#include "Functions/cgns_write_field_data.f90"
  !
  !-----------------------------------------------------------------------------
  !
end subroutine write_parallel_cgns_timeave_data_file
!
!###############################################################################
!
subroutine create_time_ave_data_names(variable_list,variable_ref)
  !
  !.. Use Statements ..
  use ovar,    only : rhoref,aref,tref
  use flowvar_mod, only : time_ave_variables
  !
  !.. Formal Arguments ..
#ifdef SPECIAL_FOR_INTEL
  character(len=:), allocatable, intent(out) :: variable_list(:)
#else
  character(len=CGLEN), allocatable, intent(out) :: variable_list(:)
#endif
  real(wp), allocatable, intent(out) :: variable_ref(:)
  !
  !.. Local Scalars ..
  integer  :: l,n,nl,nvar
  real(wp) :: ta_ref
  !
  character(len=CGLEN) :: ta_var
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: tmp_var
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_time_ave_data_names"
  !
continue
  !
  if (.not. allocated(time_ave_variables)) return
  !
  call debug_timer(entering_procedure,pname)
  !
  nvar = size(time_ave_variables)
  !
  variable_ref = init_array(nvar,one) ! F2003 AUTO-REALLOC
  !
  do n = 1,nvar
    !
    ta_ref = one
    ta_var = "ave("
    !
    tmp_var = trim(adjustl(time_ave_variables(n))) ! F2003 AUTO-REALLOC
    nl = len_trim(tmp_var)
    !
    do l = 1,nl
      !
      select case ( tmp_var(l:l) )
        !
        case ("d")
          ta_var = trim(adjustl(ta_var)) // "rho"
          ta_ref = ta_ref * rhoref
        case ("e")
          ta_var = trim(adjustl(ta_var)) // "rho*e0"
          ta_ref = ta_ref * rhoref*aref*aref
        case ("p")
          ta_var = trim(adjustl(ta_var)) // "P"
          ta_ref = ta_ref * rhoref*aref*aref
        case ("P")
          ta_var = trim(adjustl(ta_var)) // "P0"
          ta_ref = ta_ref * rhoref*aref*aref
        case ("t")
          ta_var = trim(adjustl(ta_var)) // "T"
          ta_ref = ta_ref * tref
        case ("T")
          ta_var = trim(adjustl(ta_var)) // "T0"
          ta_ref = ta_ref * tref
        case ("x")
          ta_var = trim(adjustl(ta_var)) // "vx"
          ta_ref = ta_ref * aref
        case ("y")
          ta_var = trim(adjustl(ta_var)) // "vy"
          ta_ref = ta_ref * aref
        case ("z")
          ta_var = trim(adjustl(ta_var)) // "vz"
          ta_ref = ta_ref * aref
        case ("r")
          ta_var = trim(adjustl(ta_var)) // "vr"
          ta_ref = ta_ref * aref
        case ("f")
          ta_var = trim(adjustl(ta_var)) // "vt"
          ta_ref = ta_ref * aref
          !
      end select
      !
      ta_var = trim(adjustl(ta_var)) // merge(")","*",l==nl)
      !
    end do
    !
    ! Save the final nondimensionalization value
    ! for this time-averaged variable to variable_ref
    !
    variable_ref(n) = ta_ref
    !
    ! Save the final variable name for this
    ! time-averaged variable to variable_list
    !
    call add_sol_variable(variable_list,trim(adjustl(ta_var)))
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine create_time_ave_data_names
!
!###############################################################################
!
subroutine get_connectivity_for_cgns()
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  use geomvar_mod, only : global_cell_geom
  use geomvar_mod, only : global_cell_order
  use ovar,        only : cgns_high_order_output
  use ovar,        only : cgns_mixed_output
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  logical(lk) :: mixed_element_grid
  logical(lk) :: mixed_order_grid
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_connectivity_for_cgns"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (cgns_high_order_output) then
    !
   !!
   !! Use a MIXED CGNS section only if this is a mixed element/order grid
   !! and the cgns_mixed_output flag is enabled.
   !!
   !mixed_element_grid = all_not_equal(global_cell_geom)
   !mixed_order_grid = all_not_equal(global_cell_order)
   !!
   !if ((mixed_element_grid .or. mixed_order_grid) .and. cgns_mixed_output) then
    if (cgns_mixed_output) then
      !
      call get_mixed_connectivity()
      !
    else
      !
      call get_unmixed_connectivities()
      !
    end if
    !
  else
    !
    ! Use linear subcell connectivity to output the high order solution
    !
    allocate ( cgns_sect(1:1) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cgns_sect",1,__LINE__,__FILE__,ierr,error_message)
    !
    if (nr == 1) then
      cgns_sect(1)%name = "Edges"
      cgns_sect(1)%type = BAR_2
    else if (nr == 2) then
      cgns_sect(1)%name = "Quadrilaterals"
      cgns_sect(1)%type = QUAD_4
    else if (nr == 3) then
      cgns_sect(1)%name = "Hexahedrals"
      cgns_sect(1)%type = HEXA_8
    end if
    !
    call get_linear_subcell_connectivity(cgns_sect(1)%all_cells_beg, &
                                         cgns_sect(1)%all_cells_end, &
                                         cgns_sect(1)%my_cells_beg, &
                                         cgns_sect(1)%my_cells_end, &
                                         cgns_sect(1)%ElemConn)
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_connectivity_for_cgns
!
!###############################################################################
!
subroutine get_linear_subcell_connectivity(all_cells_beg,all_cells_end, &
                                           my_cells_beg,my_cells_end,ElemConn)
  !
  !.. Use Statements ..
  use geomvar_mod,    only : ncell,n_global_cell,cell
  use geomvar_mod,    only : pst_geom
  use geomvar_mod,    only : global_cell_geom
  use geomvar_mod,    only : global_cell_order
  use geomvar_mod,    only : cell_map
  use order_mod,      only : o_order,geom_solpts
  !
  !.. Formal Arguments ..
  integer(CST),              intent(inout) :: all_cells_beg
  integer(CST),              intent(inout) :: all_cells_end
  integer(CST),              intent(inout) :: my_cells_beg
  integer(CST),              intent(inout) :: my_cells_end
  integer(CST), allocatable, intent(inout) :: ElemConn(:)
  !
  !.. Local Scalars ..
  integer :: n,nc,np,ierr,icpu
  integer :: gmin,gmax,this_subcells
  integer :: this_cell,this_geom
  integer :: this_order,output_order
  !
  integer(CST) :: p1,this_npts
  integer(CST) :: ibeg,iend,npts
  integer(CST) :: n_total_cells
  integer(CST) :: n_total_pts
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_linear_subcell_connectivity"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize n_total_cells and n_total_pts which give:
  !
  ! n_total_cells: total number of output cells for the global grid
  ! n_total_pts: total number of output points for the global grid
  !
  n_total_pts = 0_CST
  n_total_cells = 0_CST
  !
  do icpu = 1,ncpu
    !
    if (icpu-1 == mypnum) then
      my_pts_beg = n_total_pts + 1_CST
      my_cells_beg = n_total_cells + 1_CST
    end if
    !
    do nc = 1,size(cell_map(icpu)%loc_to_glb)
      !
      ! Index for grid cell on the current processor
      !
      this_cell = cell_map(icpu)%loc_to_glb(nc)
      !
      ! Geometry of this grid cell
      !
      this_geom = global_cell_geom(this_cell)
      !
      ! Polynomial order for this grid cell
      !
      this_order = global_cell_order(this_cell)
      !
      ! Number of solution points for this cell geometry at this order
      !
      output_order = max(this_order,o_order)
      !
      this_npts = int(geom_solpts(this_geom,output_order),kind=kind(this_npts))
      !
      ! Get the number of subcells for this cell geometry at this order.
      ! The number of subcells is basically:
      !     this_order**(spatial dimension of this_geom) or
      ! NOTE: Make sure that the number of subcells is at least 1,
      !       particularly for the case of this_order=0
      !
      this_subcells = &
            max( 1 , size(pst_geom(this_geom,output_order)%connectivity,dim=2) )
     !this_subcells = max( 1 , output_order**geom_dimen(this_geom) )
      !
      ! Add the number of solution points for this cell geometry at
      ! this order to the running count of the global total number of
      ! output solution points.
      !
      n_total_pts = n_total_pts + this_npts
      !
      ! Add the number of subcells for this cell geometry at this order
      ! to the running count of the global total number of output
      ! subcells.
      !
      n_total_cells = n_total_cells + &
                      int(this_subcells,kind=kind(n_total_cells))
      !
    end do
    !
    if (icpu-1 == mypnum) then
      my_pts_end = n_total_pts
      my_cells_end = n_total_cells
    end if
    !
  end do
  !
  ! Set the beginning and ending cell indices for the total number of linear
  ! sub-cells that will be output
  !
  all_cells_beg = 1_CST
  all_cells_end = n_total_cells
  !
  ! Save n_total_pts and n_total_cells to zone_size
  !
  zone_size(1) = n_total_pts
  zone_size(2) = n_total_cells
  zone_size(3) = 0_CST
  !
  ! Some quick stuff needed before creating the connectivity array
  !
  gmin = minval( global_cell_geom(1:n_global_cell) )
  gmax = maxval( global_cell_geom(1:n_global_cell) )
  !
  if (any(Geom_Unknown == [gmin,gmax])) then
    write (error_message,1)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  else if (geom_dimen(gmin) /= geom_dimen(gmax)) then
    write (error_message,2)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Get the number of nodes needed to define grid cells of this dimension
  ! NOTE: We will be creating degenerate grid cells for unstructured type
  !       cells such as Geom_Tria,Geom_Tetr,Geom_Pyra,Geom_Pris.
  ! Therefore, the number of nodes needed are:
  !         Geom_0D  => 2**0 = 1
  !         Geom_1D  => 2**1 = 2
  !         Geom_2D  => 2**2 = 4
  !         Geom_3D  => 2**3 = 8
  !
  np = 2**geom_dimen(gmin)
  !
  ! Create the connectivity array
  !
  npts = int(np,kind=kind(npts)) * (my_cells_end - my_cells_beg + 1_CST)
  !
  ElemConn = init_array(npts,0_CST) ! F2003 AUTO-REALLOC
  !
  nc = 0
  p1 = my_pts_beg
  !
  iend = 0_CST
  !
  do this_cell = 1,ncell
    !
    ! Geometry of this grid cell
    !
    this_geom = cell(this_cell)%geom
    !
    ! Polynomial order for this grid cell
    !
    this_order = cell(this_cell)%order
    !
    ! Number of solution points for this cell geometry at this order
    !
    output_order = max(this_order,o_order)
    !
    this_npts = int(geom_solpts(this_geom,output_order),kind=kind(this_npts))
    !
    ! Get the number of subcells for this cell geometry at this order.
    ! The number of subcells is basically:
    !     this_order**(spatial dimension of this_geom) or
    ! NOTE: Make sure that the number of subcells is at least 1,
    !       particularly for the case of this_order=0
    !
    this_subcells = &
          max( 1 , size(pst_geom(this_geom,output_order)%connectivity,dim=2) )
   !this_subcells = max( 1 , output_order**geom_dimen(this_geom) )
    !
    ! Loop through the number of subcells and add the subconnectivity
    ! to the first solution point of the current grid cell to get the
    ! connectivity for the current subcell.
    !
    do n = 1,this_subcells
      !
      nc = nc + 1
      !
      ibeg = iend + 1_CST
      iend = iend + int(np,kind=kind(iend))
      !
      ElemConn(ibeg:iend) =  p1 + int( pst_geom(this_geom,output_order)% &
                                                      connectivity(:,n) , &
                                       kind=kind(ElemConn) )
      !
    end do
    !
    p1 = p1 + this_npts
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (" A grid cell of an unknown geometry type was found!")
  2 format (" The minimum and maximum dimensions of the grid cells do not", &
            " match. For example, this would happen if there is a 2D cell", &
            " within a 3D grid. This configuration is currently not supported.")
  !
end subroutine get_linear_subcell_connectivity
!
!###############################################################################
!
subroutine get_mixed_connectivity()
  !
  ! USE FROM MODULE HEADER
  ! my_pts_beg
  ! my_pts_end
  ! zone_size
  ! cgns_sect
  ! geord
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_global_cell,ncell
  use geomvar_mod, only : global_cell_geom,cell_geom
  use geomvar_mod, only : global_cell_order,cell_order
  use geomvar_mod, only : cell_map
  !
  !.. Formal Arguments ..
  !
  !.. Local Scalars ..
  integer :: gmin,gmax,omin,omax,ierr,icpu,nc
  integer :: this_cell,this_geom,this_order,this_npts
  !
  integer(CST) :: n_total_pts
  integer(CST) :: n_total_cells
  integer(CST) :: my_cells_beg
  integer(CST) :: my_cells_end
  !
  logical(lk) :: mixed_geom
  logical(lk) :: mixed_order
  logical(lk) :: mixed_grid
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_mixed_connectivity"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Create the geord derived type that contains various information
  ! for the different combinations of cell geom/order. This also returns
  ! the min/max cell geometry and min/max cell order across all the cells
  ! in the current simulation.
  ! NOTE: We have to explicitly specify the indices of the geom and order
  !       argument arrays because these arrays are allocated to the number
  !       of (cells + boundary faces) and not just to the number of cells.
  !
  call create_geord(global_cell_geom(1:n_global_cell), &
                    global_cell_order(1:n_global_cell), &
                    gmin,gmax,omin,omax)
  !
  mixed_geom  = ( gmin /= gmax )
  ! Be careful here, because the order can only be mixed if (omin < 4)
  ! since everything P5 and higher will be output as P4
  mixed_order = ( (omin < 4) .and. (omin /= omax) )
  mixed_grid  = ( mixed_geom .or. mixed_order )
  !
  ! Allocate the cgns_sect derived type
  !
  allocate ( cgns_sect(1:1) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cgns_sect",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Get the name for the CGNS section.
  !
  cgns_sect(1)%name = get_mixed_section_name()
  !
  ! Get the CGNS element type for the CGNS section
  ! NOTE: We only want to use a MIXED CGNS element type if the grid actually
  !       contains more than one geom/order combo, because the element
  !       connectivity + element offset arrays for a MIXED section take up
  !       2*N_GLOBAL_CELL*CST additional bytes of storage compared to the
  !       same data stored in a regular single-element-type section.
  !
  if (mixed_grid) then
    cgns_sect(1)%type = MIXED
  else
    cgns_sect(1)%type = get_cell_type(gmin,omin)
  end if
  !
  ! Set the beginning and ending cell indices for this section
  !
  cgns_sect(1)%all_cells_beg = 1_CST
  cgns_sect(1)%all_cells_end = int(n_global_cell,kind=CST)
  !
  ! Get the total number of output points/nodes and output cells, and the
  ! beginning and ending indices for the points/nodes and cells output
  ! by this process.
  !
  n_total_pts = 0_CST
  n_total_cells = 0_CST
  !
  do icpu = 1,ncpu
    !
    ! If this process/CPU number matches my process/CPU number, save the
    ! current sum of the output points/nodes+1 and output cells+1 as the
    ! beginning index of my output points/nodes and output cells, respectively.
    !
    if (icpu-1 == mypnum) then
      my_pts_beg = n_total_pts + 1_CST
      my_cells_beg = n_total_cells + 1_CST
    end if
    !
    ! Count the number of output points/nodes and output cells for this
    ! partition.
    ! NOTE: For a MIXED section, the number of output cells is really
    !       just the number of cells on the partition.
    !
    do nc = 1,size(cell_map(icpu)%loc_to_glb)
      !
      ! Get the global cell index for the current local cell index (nc)
      !
      this_cell = cell_map(icpu)%loc_to_glb(nc)
      !
      this_geom = global_cell_geom(this_cell)   ! cell geometry
      this_order = global_cell_order(this_cell) ! cell order
      !
      ! Get the number of output points for this cell geom/order combo
      !
      this_npts = geord(this_geom,this_order)%npts
      !
      ! Add the number of points/nodes for this cell to the running sum
      ! of the total number of output points/nodes.
      !
      n_total_pts = n_total_pts + int(this_npts,kind=kind(n_total_pts))
      !
      ! Add this cell to the section cell count
      ! 
      n_total_cells = n_total_cells + 1_CST
      !
    end do
    !
    ! If this process/CPU number matches my process/CPU number, save the
    ! current sum of the output points/nodes and output cells as the ending
    ! index of my output points/nodes and output cells, respectively.
    !
    if (icpu-1 == mypnum) then
      my_pts_end = n_total_pts
      my_cells_end = n_total_cells
    end if
    !
  end do
  !
  ! Save the beginning and ending indices of my cells to cgns_sect
  !
  cgns_sect(1)%my_cells_beg = my_cells_beg
  cgns_sect(1)%my_cells_end = my_cells_end
  !
  ! Save the total number of output points and cells to zone_size
  !
  zone_size(1) = n_total_pts
  zone_size(2) = n_total_cells
  zone_size(3) = 0_CST
  !
  ! Create the CGNS element connectivity and offsets for this partition
  !
  call create_cgns_connectivity(my_pts_beg,my_pts_end, &
                                my_cells_beg,my_cells_end, &
                                cgns_sect(1)%ElemConn, &
                                cgns_sect(1)%Offsets, &
                                mixed_grid=mixed_grid)
  !
  ! Save the total size of ElemConn to MaxOffset so we dont have to extract
  ! this information when creating the CGNS section.
  !
 !cgns_sect(1)%MaxOffset = cgns_sect(1)%Offsets(size(cell_geom)+1)
 !cgns_sect(1)%MaxOffset = size(cgns_sect(1)%ElemConn,kind=CST)
  cgns_sect(1)%MaxOffset = n_total_pts + n_total_cells
  !
  if (mixed_grid) then
    !
    ! We only need to do a sanity check on MaxOffset if this is a MIXED
    ! section, because MaxOffset is only used in the call to
    ! cgp_poly_section_write_f which is only used for MIXED sections.
    !
    nc = size(cgns_sect(1)%ElemConn)
    if (ncpu > 1) then
      call mpi_allreduce(MPI_IN_PLACE,nc,1_int_mpi, &
                         mpi_inttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
    end if
    !
    if (nc /= cgns_sect(1)%MaxOffset) then
      write (error_message,1) nc,cgns_sect(1)%MaxOffset
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
    !
  end if
  !
  ! If this is not a mixed element section, deallocate cgns_sect(1)%Offsets
  ! since we dont need it.
  !
  if (.not. mixed_grid) then
    if (allocated(cgns_sect(1)%Offsets)) then
      deallocate ( cgns_sect(1)%Offsets , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"cgns_sect(1)%Offsets",2,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
  end if
  !
  if (allocated(geord)) then
    deallocate ( geord , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"geord",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Format statements
  !
  1 format ("The MaxOffset or elementDataSize input for the CGNS function/&", &
            "cgp_poly_section_write_f does not seem to match the sum of/&", &
            "the total number of output points and cells which should be/&", &
            "the case for a MIXED CGNS section./&", &
            "  the sum of the size of all local ElemConn arrays   = ",i0,'/&', &
            "  the sum of total number of output points and cells = ",i0)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_mixed_connectivity
!
!###############################################################################
!
subroutine get_unmixed_connectivities()
  !
  !.. Use Statements ..
  use geomvar_mod, only : ncell,n_global_cell
  use geomvar_mod, only : global_cell_geom,cell_geom
  use geomvar_mod, only : global_cell_order,cell_order
  use geomvar_mod, only : cell_map
  !
  !.. Local Scalars ..
  integer :: nc,icpu,gmin,gmax,omin,omax,ierr
  integer :: nsect,this_sect,this_geom,this_order
  integer :: num_cells,this_cell,this_npts,npts
  integer :: ibeg,iend
  !
  integer(CET) :: this_type
  !
  integer(CST) :: my_cells_beg,n_cells_beg
  integer(CST) :: my_cells_end,n_cells_end
  integer(CST) :: n_total_cells
  integer(CST) :: n_total_pts
  integer(CST) :: sect_total_cells
  integer(CST) :: sect_cells
  integer(CST) :: jbeg,jend,joff
  !
  !.. Allocatable Arrays ..
  integer(CET), allocatable :: global_cell_type(:)
  integer(CET), allocatable :: cell_type(:)
  integer(CST), allocatable :: ElemConn(:)
  integer(CST), allocatable :: Offsets(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_unmixed_connectivities"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Create the geord derived type that contains various information
  ! for the different combinations of cell geom/order. This also returns
  ! the min/max cell geometry and min/max cell order across all the cells
  ! in the current simulation.
  ! NOTE: We have to explicitly specify the indices of the geom and order
  !       argument arrays because these arrays are allocated to the number
  !       of (cells + boundary faces) and not just to the number of cells.
  !
  call create_geord(global_cell_geom(1:n_global_cell), &
                    global_cell_order(1:n_global_cell), &
                    gmin,gmax,omin,omax)
  !
  ! Create local and global arrays that contain the CGNS element type
  ! for all the local and global cells, respectively. This is just to
  ! make things easier instead of having to keep getting the CGNS element
  ! type from the combination of the cell geometry and cell order.
  ! NOTE: We have to explicitly specify the indices of the geom and order
  !       argument arrays because these arrays are allocated to the number
  !       of (cells + boundary faces) and not just to the number of cells.
  !
  ! F2003 AUTO-REALLOC
  global_cell_type = get_cell_type(global_cell_geom(1:n_global_cell), &
                                   global_cell_order(1:n_global_cell))
  ! F2003 AUTO-REALLOC
  cell_type = get_cell_type(cell_geom(1:ncell),cell_order(1:ncell))
  !
  ! Count the number of CGNS sections that will be needed
  !
  nsect = count( geord(:,:)%used )
  !
  ! Allocate the cgns_sect derived type
  !
  allocate ( cgns_sect(1:nsect) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cgns_sect",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the basic information for each section: section name, element type,
  ! beginning/ending indices for all elements in this section
  !
  n_cells_end = 0_CST
  !
  this_sect = 0
  do this_geom = gmin,gmax
    do this_order = omin,min(omax,4)
      !
      if (geord(this_geom,this_order)%not_used) cycle
      !
      this_sect = this_sect + 1
      !
      ! Get the section name
      !
      write (cgns_sect(this_sect)%name,'("P",i0,1x,a)') &
                            this_order, trim(adjustl(Geom_Full_Name(this_geom)))
      !
      ! Get the section element type
      !
      cgns_sect(this_sect)%type = geord(this_geom,this_order)%type
      !
      ! Set the beginning and ending element index for this section
      ! NOTE: This gives the beginning and ending index for the cells
      !       in this section within the range of [1:n_global_cell]
      !
      num_cells = geord(this_geom,this_order)%ncell
      !
      n_cells_beg = n_cells_end + 1_CST
      n_cells_end = n_cells_end + int(num_cells,kind=CST)
      !
      ! Save the beginning and ending element indices to cgns_sect
      !
      cgns_sect(this_sect)%all_cells_beg = n_cells_beg
      cgns_sect(this_sect)%all_cells_end = n_cells_end
      !
    end do
  end do
  !
  ! For each section, count the number of cells for each partition using
  ! cell_map. This also creates a running sum of the total number of cells
  ! across all partitions to make sure it matches up with the counts we
  ! just found.
  !
  n_total_cells = 0_CST
  !
  do this_sect = 1,nsect
    !
    sect_total_cells = 0_CST
    !
    do icpu = 1,ncpu
      !
      ! If this process/CPU number matches my process/CPU number, save the
      ! current sum of the section cells +1 as the beginning index of my
      ! cells for this section.
      !
      if (icpu-1 == mypnum) then
        my_cells_beg = sect_total_cells + 1_CST
      end if
      !
      ! Count the number of cells on this partition that belong to this section
      !
      do nc = 1,size(cell_map(icpu)%loc_to_glb)
        !
        ! Get the global cell index for the current local cell index (nc)
        !
        this_cell = cell_map(icpu)%loc_to_glb(nc)
        !
        ! Add this cell to the section cell count if the CGNS element type
        ! for this cell geom/order combo matches the CGNS element type
        ! for this section.
        !
        if (global_cell_type(this_cell) == cgns_sect(this_sect)%type) then
          sect_total_cells = sect_total_cells + 1_CST
        end if
        !
      end do
      !
      ! If this process/CPU number matches my process/CPU number, save the
      ! current sum of the section cells as the ending index of my cells
      ! for this section.
      !
      if (icpu-1 == mypnum) then
        my_cells_end = sect_total_cells
      end if
      !
    end do
    !
    ! Save the beginning and ending indices of my cells 
    ! for this section to cgns_sect.
    ! NOTE: We need to add an offset for the cells from the previous sections
    !
    joff = cgns_sect(this_sect)%all_cells_beg - 1_CST
    !
    cgns_sect(this_sect)%my_cells_beg = joff + my_cells_beg
    cgns_sect(this_sect)%my_cells_end = joff + my_cells_end
    !
    ! Compare the global number of cells found for this section to the count
    ! found previously, and abort if they dont match.
    !
    sect_cells = cgns_sect(this_sect)%all_cells_end - &
                 cgns_sect(this_sect)%all_cells_beg + 1_CST
    !
    if (sect_total_cells /= sect_cells) then
      write (error_message,1) this_sect,sect_total_cells,sect_cells
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
    !
    ! Update the count for the total number of output cells to include
    ! the cells found for this section.
    !
    n_total_cells = n_total_cells + sect_total_cells
    !
  end do
  !
  ! Compare the total number of output cells across all sections to the
  ! number of global cells in the simulation, and abort if they dont match.
  !
  if (n_total_cells /= int(n_global_cell,kind=kind(n_total_cells))) then
    write (error_message,2) n_total_cells,n_global_cell
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! We can deallocate global_cell_type since it is no longer needed
  !
  if (allocated(global_cell_type)) then
    deallocate ( global_cell_type , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_cell_type",1,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  ! Get the total number of output points/nodes and the beginning and ending
  ! index for the points/nodes output by this process
  ! NOTE: The output points/nodes are independent of sections
  !
  n_total_pts = 0_CST
  !
  do icpu = 1,ncpu
    !
    ! If this process/CPU number matches my process/CPU number, save the
    ! current sum of the output points/nodes +1 as the beginning index of
    ! my output points/nodes.
    !
    if (icpu-1 == mypnum) then
      my_pts_beg = n_total_pts + 1_CST
    end if
    !
    ! Count the number of output points/nodes for this partition
    !
    do nc = 1,size(cell_map(icpu)%loc_to_glb)
      !
      ! Get the global cell index for the current local cell index (nc)
      !
      this_cell = cell_map(icpu)%loc_to_glb(nc)
      !
      this_geom = global_cell_geom(this_cell)   ! cell geometry
      this_order = global_cell_order(this_cell) ! cell order
      !
      ! Get the number of output points for this cell geom/order combo
      !
      this_npts = geord(this_geom,this_order)%npts
      !
      ! Add the number of points/nodes for this cell to the running sum
      ! of the total number of output points/nodes.
      !
      n_total_pts = n_total_pts + int(this_npts,kind=kind(n_total_pts))
      !
    end do
    !
    ! If this process/CPU number matches my process/CPU number, save the
    ! current sum of the output points/nodes as the ending index of
    ! my output points/nodes.
    !
    if (icpu-1 == mypnum) then
      my_pts_end = n_total_pts
    end if
    !
  end do
  !
  ! Save the total number of output points and cells to zone_size
  !
  zone_size(1) = n_total_pts
  zone_size(2) = n_total_cells
  zone_size(3) = 0_CST
  !
  ! Create the CGNS element connectivity and offsets for this partition
  !
  call create_cgns_connectivity(my_pts_beg,my_pts_end, &
                                my_cells_beg,my_cells_end, &
                                ElemConn,Offsets)
 !!
 !! Get the total number of output points/nodes for all the cells on this
 !! partition to get the size of the element connectivity array ElemConn.
 !!
 !my_npts = int( my_pts_end - my_pts_beg + 1_CST , kind=kind(my_npts) )
 !!
 !! Create cell connectivities and offsets for all the cells on this partition
 !!
 !Offsets = init_array(ncell+1,0_CST) ! F2003 AUTO-REALLOC
 !ElemConn = init_array(my_npts,0_CST) ! F2003 AUTO-REALLOC
 !!
 !! Set the initial ending point/node index to the ending point/node index
 !! of the partition before this one, which is my_pts_beg-1.
 !!
 !iend = my_pts_beg - 1_CST
 !!
 !do nc = 1,ncell
 !  !
 !  this_geom = cell_geom(nc)
 !  this_order = cell_order(nc)
 !  this_npts = geord(this_geom,this_order)%npts
 !  !
 !  ! Save the index offset for this cell which is the final index of the
 !  ! previous cell
 !  !
 !  ioff = iend
 !  !
 !  Offsets(nc) = ioff
 !  !
 !  ! Get the beginning and ending indices for the connectivity for this cell
 !  !
 !  ibeg = iend + 1_CST
 !  iend = iend + int(this_npts,kind=kind(iend))
 !  !
 !  ! Save the element connectivity for this cell and add in the offset to
 !  ! make the connectivity indices match the point/node data for this cell
 !  !
 !  ElemConn(ibeg:iend) = ioff + &
 !         int(geord(this_geom,this_order)%pts(1:this_npts),kind=kind(ElemConn))
 !  !
 !end do
  !
  ! Now go through the sections and save the element connectivities for all
  ! the cells of a given element type to the element connectivity for the
  ! respective section in cgns_sect.
  !
  do this_sect = 1,nsect
    !
    ! Get the CGNS element type for this section
    !
    this_type = cgns_sect(this_sect)%type
    !
    ! Get the number of output points/nodes for this element type
    !
    call cg_npe_f(this_type,this_npts,cgierr)
    !
    ! Get the size of the element connectivity array for this section,
    ! which is just the number of output points/nodes for this element type
    ! multiplied by the number of cells on this partition of that type.
    !
    npts = this_npts * count(cell_type(1:ncell) == this_type)
    !
    ! If the current partition has no cells of this element type, first,
    ! set the beginning and ending cell indices for this partition to the
    ! global index for the first cell of this section, the continue to the
    ! next section so that ElemConn doesnt get allocated.
    !
    if (npts == 0) then
      joff = cgns_sect(this_sect)%all_cells_beg
      cgns_sect(this_sect)%my_cells_beg = joff
      cgns_sect(this_sect)%my_cells_end = joff
      cycle
    end if
    !
    ! Allocate the element connectivity array for the cells on this
    ! partition matching the element type of the current section.
    !
    cgns_sect(this_sect)%ElemConn = init_array(npts,0_CST) ! F2003 AUTO-REALLOC
    !
    ! Loop through the cells on this partition and add the connectivities
    ! for the cells that match the element type of the current section to
    ! the element connectivity for this section.
    !
    iend = 0
    do nc = 1,ncell
      !
      ! Skip to the next cell if the current cell is not
      ! the same element type as the current section.
      !
      if (cell_type(nc) /= this_type) cycle
      !
      ! Get the beginning and ending indices within the section
      ! element connectivity for the current cell.
      !
      ibeg = iend + 1
      iend = iend + this_npts
      !
      ! Get the beginning and ending indices within the partition
      ! element connectivity for the current cell.
      !
     !joff = Offsets(nc) - (my_pts_beg-1_CST)
      joff = Offsets(nc)
      jbeg = joff + 1_CST
      jend = joff + int(this_npts,kind=kind(jend))
      !
      ! Save the element connectivity for this cell
      ! to the element connectivity for this section.
      !
      cgns_sect(this_sect)%ElemConn(ibeg:iend) = ElemConn(jbeg:jend)
      !
    end do
    !
  end do
  !
  ! Deallocate geord and the local arrays before leaving
  !
  if (allocated(geord)) then
    deallocate ( geord , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"geord",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(cell_type)) then
    deallocate ( cell_type , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_type",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(ElemConn)) then
    deallocate ( ElemConn , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ElemConn",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(Offsets)) then
    deallocate ( Offsets , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"Offsets",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("The number of cells found for CGNS section #",i0," across/&", &
            "all partitions does not match the expected number of cells/&", &
            "determined from the global grid information./&", &
            "  # of cells found from across all partitions = ",i0,'/&', &
            "  # of cells determined from the global grid  = ",i0)
  2 format ("The total number of cells across all CGNS sections does/&", &
            "not match the total number of cells in the global grid./&", &
            "  # of cells across all CGNS sections = ",i0,'/&', &
            "  # of cells in the global grid       = ",i0)
  !
end subroutine get_unmixed_connectivities
!
!###############################################################################
!
subroutine create_geord(cell_geom,cell_order,gmin,gmax,omin,omax)
  !
  !.. Use Statements ..
  use cgns_mod, only : reorder_to_cgns_pts
  !
  !.. Formal Arguments ..
  integer,    intent(in) :: cell_geom(:)
  integer,    intent(in) :: cell_order(:)
  integer, intent(inout) :: gmin
  integer, intent(inout) :: gmax
  integer, intent(inout) :: omin
  integer, intent(inout) :: omax
  !
  !.. Local Scalars ..
  integer      :: nc,igeo,iord,npts,ierr
  integer(CET) :: ityp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_geord"
  !
continue
  !
  gmin = minval(cell_geom)
  gmax = maxval(cell_geom)
  omin = minval(cell_order)
  omax = maxval(cell_order)
  !
  ! Make sure the minimum order is not greater than 4 because CGNS
  ! does not support anything higher than P4
  !
  omin = min(4,omin)
  !
  ! Allocate the geord derived type
  !
  allocate ( geord(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"geord",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Count the number of cells for each geom/order combination
  !
  do nc = 1,size(cell_geom)
    igeo = cell_geom(nc)
    iord = cell_order(nc)
    geord(igeo,iord)%ncell = geord(igeo,iord)%ncell + 1
  end do
  !
  ! If there are any cells P5 or higher, add these cells to the P4 count for
  ! the respective geometry and set the higher-than-P4 cell count to zero.
  ! This is because CGNS only supports up to P4 cells, so any cells P5 or
  ! higher will have to be interpolated to the P4 output points.
  !
  do igeo = gmin,gmax
    if (omax >= 5) then
      geord(igeo,4)%ncell = sum(geord(igeo,4:omax)%ncell)
      geord(igeo,5:omax)%ncell = 0
    end if
  end do
  !
  ! Go through each geom/order combo and indicate if the combo exists
  ! in the current simulation
  !
  do igeo = gmin,gmax
    do iord = omin,omax
      if (geord(igeo,iord)%ncell > 0) then
        geord(igeo,iord)%used = true
        geord(igeo,iord)%not_used = fals
      end if
    end do
  end do
  !
  ! Get the element type and number of output points/nodes for each
  ! geom/order combo
  !
  do igeo = gmin,gmax
    do iord = omin,omax
      geord(igeo,iord)%type = get_cell_type(igeo,iord)
      call cg_npe_f(geord(igeo,iord)%type, &
                    geord(igeo,iord)%npts, &
                    cgierr)
    end do
  end do
  !
  ! Create the the point connectivity for each present geom/order combo. The
  ! indices in the point connectivity make it so that the current point data
  ! can be written directly to the CGNS file without reordering, and the indices
  ! themselves reorder the point data to the order defined for the element type
  ! by CGNS. The indices are the same for all cells of a particular geom/order
  ! combo, except for an offset based on the number of points from all previous
  ! cells. Creating this derived type allows us to create these indices once,
  ! and we can just add the offset as the full cell connectivity is created.
  !
  do igeo = gmin,gmax
    do iord = omin,omax
      !
      ! Skip to the next geom/order combo if there are no cells for this combo
      !
      if (geord(igeo,iord)%not_used) cycle
      !
      ityp = geord(igeo,iord)%type ! CGNS element type for this combo
      npts = geord(igeo,iord)%npts ! num of output points/nodes for this combo
      !
      ! Get the point connectivity for this geom/order combo
      !
      ! F2003 AUTO-REALLOC
      geord(igeo,iord)%pts = int( reorder_to_cgns_pts(ityp,intseq(1,npts)) , &
                                  kind=kind(geord(igeo,iord)%pts) )
      !
    end do
  end do
  !
end subroutine create_geord
!
!###############################################################################
!
elemental function get_cell_type(geom,order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: geom
  integer, intent(in) :: order
  !
  !.. Function Result ..
  integer(CET) :: return_value
  !
  !.. Local Scalars ..
  integer :: cap_order
  !
  !.. Local Parameters ..
  integer(CET), parameter :: edge(1:4) = [BAR_2,  BAR_3,   BAR_4,   BAR_5]
  integer(CET), parameter :: tria(1:4) = [TRI_3,  TRI_6,   TRI_10,  TRI_15]
  integer(CET), parameter :: quad(1:4) = [QUAD_4, QUAD_9,  QUAD_16, QUAD_25]
  integer(CET), parameter :: tetr(1:4) = [TETRA_4,TETRA_10,TETRA_20,TETRA_35]
  integer(CET), parameter :: pyra(1:4) = [PYRA_5, PYRA_14, PYRA_30, PYRA_55]
  integer(CET), parameter :: pris(1:4) = [PENTA_6,PENTA_18,PENTA_40,PENTA_75]
  integer(CET), parameter :: hexa(1:4) = [HEXA_8, HEXA_27, HEXA_64, HEXA_125]
  !
continue
  !
  ! We need to cap order at 4 since CGNS only supports up to P4 elements
  ! Also make sure order is at least P1
  !
  cap_order = max( 1 , min(order,4) )
  !
  if (geom == Geom_Node) then
    return_value = NODE ! I doubt this will be needed
  else if (geom == Geom_Edge) then
    return_value = edge( cap_order )
  else if (geom == Geom_Tria) then
    return_value = tria( cap_order )
  else if (geom == Geom_Quad) then
    return_value = quad( cap_order )
  else if (geom == Geom_Tetr) then
    return_value = tetr( cap_order )
  else if (geom == Geom_Pyra) then
    return_value = pyra( cap_order )
  else if (geom == Geom_Pris) then
    return_value = pris( cap_order )
  else if (geom == Geom_Hexa) then
    return_value = hexa( cap_order )
  else
    return_value = -1_CET
  end if
  !
end function get_cell_type
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine get_geom_order_from_cgns_type(type,geom,order)
  !
  !.. Formal Arguments ..
  integer(CET), intent(in) :: type
  integer,     intent(out) :: geom
  integer,     intent(out) :: order
  !
  !.. Local Parameters ..
  integer(CET), parameter :: edge(1:4) = [BAR_2,  BAR_3,   BAR_4,   BAR_5]
  integer(CET), parameter :: tria(1:4) = [TRI_3,  TRI_6,   TRI_10,  TRI_15]
  integer(CET), parameter :: quad(1:4) = [QUAD_4, QUAD_9,  QUAD_16, QUAD_25]
  integer(CET), parameter :: tetr(1:4) = [TETRA_4,TETRA_10,TETRA_20,TETRA_35]
  integer(CET), parameter :: pyra(1:4) = [PYRA_5, PYRA_14, PYRA_30, PYRA_55]
  integer(CET), parameter :: pris(1:4) = [PENTA_6,PENTA_18,PENTA_40,PENTA_75]
  integer(CET), parameter :: hexa(1:4) = [HEXA_8, HEXA_27, HEXA_64, HEXA_125]
  !
continue
  !
  if (any(type == edge(:))) then
    geom = Geom_Edge
    order = findloc(edge,type,dim=1)
  else if (any(type == tria(:))) then
    geom = Geom_Tria
    order = findloc(tria,type,dim=1)
  else if (any(type == quad(:))) then
    geom = Geom_Quad
    order = findloc(quad,type,dim=1)
  else if (any(type == tetr(:))) then
    geom = Geom_Tetr
    order = findloc(tetr,type,dim=1)
  else if (any(type == pyra(:))) then
    geom = Geom_Pyra
    order = findloc(pyra,type,dim=1)
  else if (any(type == pris(:))) then
    geom = Geom_Pris
    order = findloc(pris,type,dim=1)
  else if (any(type == hexa(:))) then
    geom = Geom_Hexa
    order = findloc(hexa,type,dim=1)
  else
    geom = Geom_Unknown
    order = -1
  end if
  !
end subroutine get_geom_order_from_cgns_type
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_mixed_section_name() result(return_value)
  !
  !.. Function Result ..
  character(len=CGLEN) :: return_value
  !
  !.. Local Scalars ..
  integer :: igeo,iord
  integer :: found_geom,found_order
  !
  logical(lk) :: mixed_geom,geom_not_set
  logical(lk) :: mixed_order,order_not_set
  !
  character(len=13) :: cord
  character(len=15) :: cgeo
  !
continue
  !
  cgeo = "Elements"
  cord = "Unknown-Order"
  !
  mixed_geom = fals
  mixed_order = fals
  !
  geom_not_set = true
  order_not_set = true
  !
  do igeo = lbound(geord,dim=1),ubound(geord,dim=1)
    do iord = lbound(geord,dim=2),ubound(geord,dim=2)
      !
      ! Skip to the next geom/order combo if there are no cells for this combo
      !
      if (geord(igeo,iord)%not_used) cycle
      !
      if (.not. mixed_geom) then
        if (geom_not_set) then
          cgeo = trim(adjustl(Geom_Full_Name(igeo)))
          geom_not_set = fals
          found_geom = igeo
        else
          if (igeo /= found_geom) then
            cgeo = "Mixed-Element"
            mixed_geom = true
          end if
        end if
      end if
      !
      if (.not. mixed_order) then
        if (order_not_set) then
          write (cord,'("P",i0)') iord
          order_not_set = fals
          found_order = iord
        else
          if (iord /= found_order) then
            cord = "Mixed-Order"
            mixed_order = true
          end if
        end if
      end if
      !
     !if (mixed_geom) then
     !  ! There are multiple geometries present so set the geometry part of
     !  ! the section name to indicate its a mixed element section
     !  cgeo = "Mixed-Element"
     !else
     !  found_geom = igeo
     !  ! This is the first geometry found, so set the geometry part of the
     !  ! section name to the name of this geometry
     !  cgeo = trim(adjustl(Geom_Full_Name(igeo)))
     !  ! Set mixed_geom to 1 so that we know a geometry has already been found
     !  mixed_geom = true
     !end if
     !!
     !if (mixed_order) then
     !  ! There are multiple orders present so set the order part of
     !  ! the section name to indicate its a mixed order section
     !  cord = "Mixed-Order"
     !else
     !  ! This is the first order found, so set the order part of the
     !  ! section name to this order
     !  write (cord,'("P",i0)') iord
     !  ! Set mixed_order to 1 so that we know an order has already been found
     !  mixed_order = true
     !end if
      !
    end do
  end do
  !
  ! Write the section name to the function result
  !
  write (return_value,'(a,1x,a)') trim(adjustl(cord)), trim(adjustl(cgeo))
  !
  return_value = trim(adjustl(return_value))
 !return_value = "Hexahedrals"
  !
end function get_mixed_section_name
!
!###############################################################################
!
subroutine create_cgns_connectivity(ptsbeg,ptsend,celbeg,celend, &
                                    ElemConn,Offsets,mixed_grid)
  !
  !.. Use Statements ..
  use geomvar_mod, only : ncell,cell_geom,cell_order
  !
  !.. Formal Arguments ..
  integer(CST),               intent(in) :: ptsbeg
  integer(CST),               intent(in) :: ptsend
  integer(CST),               intent(in) :: celbeg
  integer(CST),               intent(in) :: celend
  integer(CST), allocatable, intent(out) :: ElemConn(:)
  integer(CST), allocatable, intent(out) :: Offsets(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: mixed_grid
  !
  !.. Local Scalars ..
  integer :: nc,igeo,iord
  integer(CST) :: ibeg,iend,ioff,npts,my_npts
  integer(CST) :: pts_offset,offset_for_elem_type
  logical(lk) :: mixed_elem_grid
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_cgns_connectivity"
  !
continue
  !
  mixed_elem_grid = fals
  if (present(mixed_grid)) then
    mixed_elem_grid = mixed_grid
  end if
  !
  ! Get the total number of output points/nodes and cells for this partition
  !
  my_npts = ptsend - ptsbeg + 1_CST
  !
  ! If this is actually a MIXED section, add ncell to my_npts because the
  ! element connectivity for a MIXED section must include the element type
  ! before the actual connectivity for each cell.
  !
  if (mixed_elem_grid) then
    my_npts = my_npts + int(ncell,kind=kind(my_npts))
  end if
  !
  ! Create cell connectivities and offsets for all the cells on this partition
  !
  Offsets = init_array(ncell+1,0_CST) ! F2003 AUTO-REALLOC
  ElemConn = init_array(my_npts,0_CST) ! F2003 AUTO-REALLOC
  !
  ! Create an integer binary to use for a mixed element section that is used
  ! to add the +1 offset needed for storing the element type in the element
  ! connectivity array.
  !
  offset_for_elem_type = merge(1_CST,0_CST,mixed_elem_grid)
  !
  ! Initialize the points/node offset counting index for this partition
  !
  pts_offset = ptsbeg - 1_CST
  !
  ! Initialize the ioff counting index which is basically the offset value
  ! for the Offsets array. The Offsets array is needed for 2 situations:
  !  1) Writing a MIXED section, in which ioff needs to be
  !     pts_offset + celbeg - 1
  !  2) Writing multiple sections each having a fixed element type. For this
  !     situation, Offsets is only used to copy only the connectivity for cells
  !     of a certain type from the ElemConn array into the respective
  !     cgns_sect(i)%ElemConn array. For this use, Offsets should start at
  !     zero since it is only for cells local to the partition.
  !  3) The last situation is when using the cgns_mixed_output for a grid that
  !     only contains a single element type. For this situation, the Offsets
  !     array wont be used, so it doesnt matter what ioff is initialized to.
  !
 !ioff = pts_offset + offset_for_elem_type*(celbeg - 1_CST)
  ioff = offset_for_elem_type*(pts_offset + celbeg - 1_CST)
  !
  ! Initialize the counting index for the ElemConn array
  !
  iend = 0_CST
  !
  ! Loop over the cells and create the element connectivity for each cell
  ! and the offset array pointing to starting location of the connectivity
  ! for a given cell within the connectivity array.
  !
  do nc = 1,ncell
    !
    igeo = cell_geom(nc)
    iord = cell_order(nc)
    npts = int(geord(igeo,iord)%npts,kind=kind(npts))
    !
    ! Save the index offset for this cell which is the final index of the
    ! previous cell
    !
    Offsets(nc) = ioff
    !
    ! Get the beginning and ending indices for the connectivity for this cell
    !
    ibeg = iend + 1_CST
    iend = iend + npts + offset_for_elem_type
    !
    ! Update the offset value for the next cell
    !
    ioff = ioff + npts + offset_for_elem_type
    !
    ! Add the element type for this cell to the ibeg location of ElemConn
    ! if this is for a mixed element section.
    !
    if (mixed_elem_grid) then
      ElemConn(ibeg) = int(geord(igeo,iord)%type,kind=kind(ElemConn))
    end if
    !
    ! Add offset_for_elem_type to ibeg to increase it by one if this is
    ! for a mixed element section so that the array indices (ibeg:iend)
    ! match the size of geord(igeo,iord)%pts(1:npts).
    !
    ibeg = ibeg + offset_for_elem_type
    !
    ! Save the element connectivity for this cell and add in the offset to
    ! make the connectivity indices match the point/node data for this cell
    !
    ElemConn(ibeg:iend) = pts_offset + &
           int(geord(igeo,iord)%pts(1:npts),kind=kind(ElemConn))
    !
    ! Update the points/node offset index for the next cell.
    ! We need to do this separately from iend because these
    ! two offsets will be different for a mixed element section.
    !
    pts_offset = pts_offset + npts
    !
  end do
  !
  ! The last element of Offsets should give the total size of the ElemConn array
  !
  Offsets(ncell+1) = ioff
  !
  ! Check to make sure the final ioff counter matches the size of ElemConn
  !
  if (Offsets(ncell+1)-Offsets(1) /= my_npts) then
    write (error_message,1) Offsets(ncell+1)-Offsets(1),my_npts
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Format Statements
  !
  1 format ("The final offset counter does not match the size of the/&", &
            "connectivity array!/&", &
            "  final offset counter = ",i0,'/&', &
            "        size(ElemConn) = ",i0)
  !
end subroutine create_cgns_connectivity
!
!###############################################################################
!
subroutine TaylorGreen_8s_vorticity_solution()
  !
  !.. Use Statements ..
  use ovar, only : write_TaylorGreen_8s_solution
  use ovar, only : profile_io_cgns
  use ovar, only : cgns_dir
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  integer(CBT) :: ifile_v,ibase_v,izone_v
  integer(CBT) :: isect_v,isol_v,isvar_v
  !
  integer(CST) :: my_beg,my_end,my_pts
  !
  real(wp) :: beg_wall_time
  real(wp) :: end_wall_time
  real(wp) :: cgns_wall_time
  !
  character(len=170) :: cgns_vort_file
  !
  integer(INT_MPI) :: color,key
  _MPI_COMM_TYPE_ :: vort_comm
  logical(lk) :: i_am_vort_cpu
  !
  !.. Local Arrays ..
  integer(CST), dimension(1:3) :: my_vort_pts
  integer(CST), dimension(1:3) :: my_vort_cells
  !
  integer(CBT), dimension(1:2) :: ixyz_v
  integer(CST), dimension(1:3) :: vort_zone_size
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:) :: my_vort_faces
  !
  integer(CST), allocatable, dimension(:,:) :: ipelem
  !
  real(CRT), allocatable, dimension(:) :: vary
  real(CRT), allocatable, dimension(:) :: varz
  !
  !.. Local Parameters ..
  character(len=*), parameter :: bname = "TGV_vort_8s_base"
  character(len=*), parameter :: zname = "TGV_vort_8s_zone"
  character(len=*), parameter :: fname = "TGV_vorticity_8s.cgns"
  character(len=*), parameter :: pname = "TaylorGreen_8s_vorticity_solution"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
#ifndef SPECIAL_FOR_PGI
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      write (iout,1) pname
      flush (iout)
      beg_wall_time = mpi_wtime()
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  cgns_vort_file = trim(adjustl(cgns_dir)) // fname
  !
  ! Get/collect all the data necessary to output the 8s vorticity solution
  !
  call get_data_for_8s_vorticity_solution(my_vort_pts,my_vort_cells, &
                                          vort_zone_size,my_vort_faces,ipelem)
  !
  my_pts = my_vort_pts(1)
  my_beg = my_vort_pts(2)
  my_end = my_vort_pts(3)
  !
  i_am_vort_cpu = (my_pts > 0)
  color = merge( 1_INT_MPI , MPI_UNDEFINED , i_am_vort_cpu )
  key = 0_INT_MPI
  !
  call mpi_comm_split(MPI_COMM_WORLD,color,key,vort_comm,mpierr)
  !
  if (i_am_vort_cpu) then
    !
#ifdef USE_MPI_F08
    call cgp_mpi_comm_f(vort_comm%mpi_val,cgierr)
#else
    call cgp_mpi_comm_f(vort_comm,cgierr)
#endif
    !
    ! Open the CGNS file
    !
    call cgp_open_f(cgns_vort_file,CG_MODE_WRITE,ifile_v,cgierr)
    call io_error(pname,cgns_vort_file,1,__LINE__,__FILE__,cgierr)
    !
    ! Set the CGNS base (who knows what this is for)
    ! NOTE: We are creating this base 2-dimensional with surface (2D) cells
    !       since this file will contain a 2D solution along a boundary face
    !
    call cg_base_write_f(ifile_v,bname,2_CBT,2_CBT,ibase_v,cgierr)
    call cgns_error(pname,cgns_vort_file,"cg_base_write_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    !
    ! Write the title for this solution
    !
    call cg_zone_write_f(ifile_v,ibase_v,zname,vort_zone_size, &
                         UNSTRUCTURED,izone_v,cgierr)
    call cgns_error(pname,cgns_vort_file,"cg_zone_write_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    !
    ! Create data nodes for the coordinates
    !
    call cgp_coord_write_f(ifile_v,ibase_v,izone_v,CGNS_KIND_REAL, &
                           "CoordinateY",ixyz_v(1),cgierr)
    call cgns_error(pname,cgns_vort_file,"cgp_coord_write_f", &
                    cgierr,__LINE__,__FILE__,"Write node: Coordinate = 1 (Y)", &
                    not_parallel=true)
    !
    call cgp_coord_write_f(ifile_v,ibase_v,izone_v,CGNS_KIND_REAL, &
                           "CoordinateZ",ixyz_v(2),cgierr)
    call cgns_error(pname,cgns_vort_file,"cgp_coord_write_f", &
                    cgierr,__LINE__,__FILE__,"Write node: Coordinate = 2 (Z)", &
                    not_parallel=true)
    !
    ! Now create and write the coordinate data in parallel
    !
    ! For each boundary face within the grid region of interest, interpolate the
    ! coordinates of the solution points for the host cell to the boundary face
    ! and this interpolate those face points to the output face points.
    !
    allocate ( vary(1:my_pts) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vary",1,__LINE__,__FILE__,ierr,error_message)
    vary = real( zero , kind=kind(vary) )
    !
    allocate ( varz(1:my_pts) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"varz",1,__LINE__,__FILE__,ierr,error_message)
    varz = real( zero , kind=kind(varz) )
    !
    ! Interpolate the coordinates to the output points
    !
    call get_vort_face_coordinates(my_vort_faces,vary,varz)
    !
    ! Add the data for the y-coordinate to the output queue
    ! NOTE: It is an error for a processor to try an write data if it doesnt
    !       have any data to write, so make sure that only those processors
    !       that have data to write call the CGNS write subroutine.
    !
    cgierr = 0
    call cgp_coord_write_data_f(ifile_v,ibase_v,izone_v,ixyz_v(1), &
                                my_beg,my_end,vary,cgierr)
    call cgns_error(pname,cgns_vort_file,"cgp_coord_write_data_f", &
                    cgierr,__LINE__,__FILE__,"Write data: Coordinate = 1 (Y)", &
                    not_parallel=true)
    !
    ! Add the data for the z-coordinate to the output queue
    ! NOTE: It is an error for a processor to try an write data if it doesnt
    !       have any data to write, so make sure that only those processors
    !       that have data to write call the CGNS write subroutine.
    !
    cgierr = 0
    call cgp_coord_write_data_f(ifile_v,ibase_v,izone_v,ixyz_v(2), &
                                my_beg,my_end,varz,cgierr)
    call cgns_error(pname,cgns_vort_file,"cgp_coord_write_data_f", &
                    cgierr,__LINE__,__FILE__,"Write data: Coordinate = 2 (Z)", &
                    not_parallel=true)
    !
    ! Create the data node for the grid connectivity
    !
    call cgp_section_write_f(ifile_v,ibase_v,izone_v,"Quadrilaterals", &
                             QUAD_4,1_CST,vort_zone_size(2), &
                             0_CBT,isect_v,cgierr)
    call cgns_error(pname,cgns_vort_file,"cgp_section_write_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    !
    ! Write the element connectivity in parallel
    !
    my_beg = my_vort_cells(2)
    my_end = my_vort_cells(3)
    !
    ! Add the element connectivity data to the output queue.
    ! NOTE: It is an error for a processor to try an write data if it doesnt
    !       have any data to write, so make sure that only those processors
    !       that have data to write call the CGNS write subroutine.
    !
    cgierr = 0
    call cgp_elements_write_data_f(ifile_v,ibase_v,izone_v,isect_v, &
                                   my_beg,my_end,ipelem,cgierr)
    call cgns_error(pname,cgns_vort_file,"cgp_elements_write_data_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    !
    ! Deallocate the arrays varz, and ipelem
    ! NOTE: Keep vary allocated and use it to store the vorticity variable
    !
    if (allocated(varz)) then
      deallocate ( varz , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"varz",2,__LINE__,__FILE__,ierr,error_message)
    end if
    if (allocated(ipelem)) then
      deallocate ( ipelem , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"ipelem",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    ! Evaluate and write the vorticity solution
    !
    my_beg = my_vort_pts(2)
    my_end = my_vort_pts(3)
    !
    ! Create the solution node
    !
    call cg_sol_write_f(ifile_v,ibase_v,izone_v,"VorticitySolution", &
                        VERTEX,isol_v,cgierr)
    call cgns_error(pname,cgns_vort_file,"cg_sol_write_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    !
    ! Create a data node for the vorticity solution
    !
    call cgp_field_write_f(ifile_v,ibase_v,izone_v,isol_v,CGNS_KIND_REAL, &
                           "VorticityMagnitude",isvar_v,cgierr)
    call cgns_error(pname,cgns_vort_file,"cgp_field_write_f", &
                    cgierr,__LINE__,__FILE__, &
                    "Write node: Solution variable = 1 (VorticityMagnitude)", &
                    not_parallel=true)
    !
    ! Compute the vorticity magnitude at each face point
    !
    call get_vort_face_solution(my_vort_faces,vary)
    !
    ! Add the data for current variable to the output queue.
    ! NOTE: It is an error for a processor to try an write data if it doesnt
    !       have any data to write, so make sure that only those processors
    !       that have data to write call the CGNS write subroutine.
    !
    cgierr = 0
    call cgp_field_write_data_f(ifile_v,ibase_v,izone_v,isol_v,isvar_v, &
                                my_beg,my_end,vary,cgierr)
    call cgns_error(pname,cgns_vort_file,"cgp_field_write_data_f", &
                    cgierr,__LINE__,__FILE__, &
                    "Write data: Solution variable = 1 (VorticityMagnitude)", &
                    not_parallel=true)
    !
    ! Deallocate vary now that we are finished with it
    !
    if (allocated(vary)) then
      deallocate ( vary , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"vary",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    ! Close the CGNS file
    !
    call cgp_close_f(ifile_v,cgierr)
    call io_error(pname,cgns_vort_file,2,__LINE__,__FILE__,cgierr)
    !
  end if
  !
  call cgp_mpi_comm_f(mpi_comm_world_val,cgierr)
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      end_wall_time = mpi_wtime()
      cgns_wall_time = end_wall_time - beg_wall_time
      write (iout,2) pname,trim(adjustl(cgns_vort_file)),cgns_wall_time
      flush (iout)
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
#endif
  !
  write_TaylorGreen_8s_solution = fals
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," Entering ",a)
  2 format (1x,"Leaving ",a,/,:, &
            5x,"Wall time taken to write the file '",a, &
            "' = ",f12.6," seconds",/)
  !
end subroutine TaylorGreen_8s_vorticity_solution
!
!###############################################################################
!
subroutine get_data_for_8s_vorticity_solution(my_vort_pts, &
                                              my_vort_cells, &
                                              vort_zone_size, &
                                              my_vort_faces, &
                                              ipelem)
  !
  !.. Use Statements ..
  use ovar,        only : Lref,VortexGrid_DOF
  use geomvar_mod, only : nfbnd,face
  use geomvar_mod, only : xyz_nodes,pst_geom
  use order_mod,   only : geom_solpts,n_order
  !
  !.. Formal Arguments ..
  integer(CST),              intent(inout) :: my_vort_pts(1:3)
  integer(CST),              intent(inout) :: my_vort_cells(1:3)
  integer(CST),              intent(inout) :: vort_zone_size(1:3)
  integer,      allocatable, intent(inout) :: my_vort_faces(:)
  integer(CST), allocatable, intent(inout) :: ipelem(:,:)
  !
  !.. Local Scalars ..
  integer :: nc,nf,n,ierr,ip
  integer :: this_face,face_geom,face_order
  integer :: face_subcells
  !
  integer(INT_MPI) :: int_size,cvd_size
  _MPI_DATA_TYPE_ :: mpi_vorttyp
  !
  integer(CST) :: p1,face_npts
  integer(CST) :: n_my_vort_pts,n_my_vort_cells
  integer(CST) :: my_vort_pts_beg,my_vort_cells_beg
  integer(CST) :: my_vort_pts_end,my_vort_cells_end
  integer(CST) :: n_total_vort_cells,n_total_vort_pts
  !
  real(wp) :: LPI,half_LPI,dx_eps,cell_dx,x,y,z
  !
  !.. Local Allocatable Arrays ..
  logical(lk),      allocatable, dimension(:)   :: maskf
  integer(CST), allocatable, dimension(:,:) :: cpu_vort_data
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_data_for_8s_vorticity"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Compute a few initial constants relating to the grid length
  !
  LPI = Lref * PI
  half_LPI = half * LPI
  !
  ! Compute the average cell height and then use that to get a small
  ! percentage of the cell height to use a limit to identify boundary
  ! faces within the grid region of interest
  !
  cell_dx = two * LPI / anint( real(VortexGrid_DOF,kind=wp) / &
                               real(n_order+1,kind=wp) , kind=wp )
  dx_eps = cell_dx * eps3 / ten
  !
  ! Initialize the variables giving the total number of vorticity points
  ! and subcells as well as the variables giving the same information for
  ! the respective information local to the current processor
  !
  n_total_vort_pts = 0_CST
  n_total_vort_cells = 0_CST
  !
  n_my_vort_pts = 0_CST
  n_my_vort_cells = 0_CST
  !
  allocate ( maskf(1:nfbnd) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"maskf",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Find the boundary faces that lie on the plane x = -PI*Lref within the region
  ! y = [0,0.5*PI*Lref] and z = [0.5*PI*Lref,PI*Lref]
  !
  maskf_loop: do nf = 1,nfbnd
    !
    bface_node_loop: do n = 1,size(face(nf)%nodes)
      !
      face_geom = face(nf)%geom
      face_order = face(nf)%order ! order of boundary face
      !
      ! NOTE: We dont care here if any of the nodes stored in face(nf)%nodes
      !       for this face are duplicates. The if statements will just produce
      !       a false value for the second time.  We only care if the value
      !       ip that is returned is zero because this will cause the xyz_nodes
      !       array to go out-of-bounds so cycle the inner loop if this
      !       situation arrises.
      !
      ip = face(nf)%nodes(n)
      !
      if (ip == 0) cycle bface_node_loop
      !
      x = xyz_nodes(1,ip)
      y = xyz_nodes(2,ip)
      z = xyz_nodes(3,ip)
      !
      if (abs(x + LPI) <= eps6) then
        if (z >= half_LPI + dx_eps) then
          if ( (y >= dx_eps) .and. (y <= half_LPI - dx_eps) ) then
            !
            ! Mark this face for quick identification later
            !
            maskf(nf) = true
            !
            ! Number of solution points for this face geometry at this order
            !
            face_npts = int(geom_solpts(face_geom,face_order), &
                            kind=kind(face_npts))
            !
            ! Get the number of subcells for this face geometry at this order.
            ! The number of subcells is basically:
            !     face_order**(spatial dimension of face_geom) or
            ! NOTE: Make sure that the number of subcells is at least 1,
            !       particularly for the case of face_order=0
            !
            face_subcells = max( 1 , face_order**geom_dimen(face_geom) )
            !
            ! Add the number of solution points and subcells to the total
            ! number local to this processor
            !
            n_my_vort_pts   = n_my_vort_pts   + face_npts
            n_my_vort_cells = n_my_vort_cells + &
                              int(face_subcells,kind=kind(n_my_vort_cells))
            !
            cycle maskf_loop
            !
          end if
        end if
      end if
      !
    end do bface_node_loop
    !
  end do maskf_loop
  !
  ! Create the array that stores the index for all boundary faces on this
  ! processor that lie within the grid region of interest
  !
  if (allocated(my_vort_faces)) then
    deallocate ( my_vort_faces , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"my_vort_faces",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  allocate ( my_vort_faces(1:count(maskf)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"my_vort_faces",1,__LINE__,__FILE__,ierr,error_message)
  !
  my_vort_faces = pack( intseq(1,nfbnd) , mask=maskf )
  !
  ! Deallocate the maskf array since it is no longer needed
  !
  deallocate ( maskf , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"maskf",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! We need to collect this information for all processors
  !
  allocate ( cpu_vort_data(1:2,1:ncpu) , source=0_CST , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cpu_vort_data",1,__LINE__,__FILE__,ierr,error_message)
  !
  cpu_vort_data(1,mypnum+1) = n_my_vort_pts
  cpu_vort_data(2,mypnum+1) = n_my_vort_cells
  !
  if (ncpu > 1) then
    !
#ifdef USE_MPI_F08
    call mpi_sizeof(reshape(cpu_vort_data,[size(cpu_vort_data)]), &
                    int_size,mpierr)
    call mpi_type_match_size(MPI_TYPECLASS_INTEGER,int_size,mpi_vorttyp,mpierr)
#else
    if (CST == I8) then
      mpi_vorttyp = MPI_INTEGER8
    else
      mpi_vorttyp = MPI_INTEGER
    end if
#endif
    !
    cvd_size = size(cpu_vort_data,dim=1,kind=int_mpi)
    call mpi_allgather(MPI_IN_PLACE, cvd_size,mpi_vorttyp, &
                       cpu_vort_data,cvd_size,mpi_vorttyp, &
                       MPI_COMM_WORLD,mpierr)
    !
  end if
  !
  ! Sum up all the vorticity points and subcells from all processors
  !
  n_total_vort_pts   = sum( cpu_vort_data(1,:) )
  n_total_vort_cells = sum( cpu_vort_data(2,:) )
  !
  ! Save the total number of vorticity points and subcells to vort_zone_size
  !
  vort_zone_size = [n_total_vort_pts,n_total_vort_cells,0_CST]
  !
  ! Get the beginning index for the vorticity points and subcells local to this
  ! processor by summing up all the vorticity points and subcells from all
  ! processors with a smaller rank than the this processor.
  ! NOTE: The second index of cpu_vort_data goes from 1:ncpu whereas the value
  !       given by mypnum is from 0:ncpu-1. Therefore, summing from 1:mypnum
  !       gives the range of all processors with a smaller rank.
  !
  my_vort_pts_beg   = sum( cpu_vort_data(1,1:mypnum) ) + 1_CST
  my_vort_cells_beg = sum( cpu_vort_data(2,1:mypnum) ) + 1_CST
  !
  ! Get the ending index for the vorticity points and subcells local to this
  ! processor by adding the numbers found in fmask_loop to the beginning indices
  !
  my_vort_pts_end   = (my_vort_pts_beg-1)   + n_my_vort_pts
  my_vort_cells_end = (my_vort_cells_beg-1) + n_my_vort_cells
  !
  ! Save these values in the my_vort_pts and my_vort_cells arrays
  !
  my_vort_pts   = [n_my_vort_pts,  my_vort_pts_beg,  my_vort_pts_end]
  my_vort_cells = [n_my_vort_cells,my_vort_cells_beg,my_vort_cells_end]
  !
  ! Deallocate the cpu_vort_data array since it is no longer needed
  !
  deallocate ( cpu_vort_data , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cpu_vort_data",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the connectivity array
  ! NOTE: The Taylor-Green vortex case is 3D so we know the boundary faces
  !       will always be 2D and they will require 4 nodes to define each
  !       subcell within the CGNS file
  !
  if (allocated(ipelem)) then
    deallocate ( ipelem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ipelem",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( ipelem(1:4,1:n_my_vort_cells) , source=0_CST , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ipelem",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Initialize the counters nc and p1 for this processor
  !
  nc = 0
  p1 = my_vort_pts_beg
  !
  ! Loop through the boundary faces on this processor that are within the
  ! grid region of interest and get the subcell connectivity for each face
  !
  do nf = 1,size(my_vort_faces)
    !
    this_face = my_vort_faces(nf)
    !
    face_geom = face(this_face)%geom   ! geometry of boundary face
    face_order = face(this_face)%order ! order of boundary face
    !
    ! Number of solution points for this face geometry at this order
    !
    face_npts = int(geom_solpts(face_geom,face_order),kind=kind(face_npts))
    !
    ! Get the number of subcells for this face geometry at this order.
    ! The number of subcells is basically:
    !     face_order**(spatial dimension of face_geom) or
    ! NOTE: Make sure that the number of subcells is at least 1,
    !       particularly for the case of face_order=0
    !
    face_subcells = max( 1 , face_order**geom_dimen(face_geom) )
    !
    ! Loop through the number of subcells and add the subconnectivity
    ! to the first solution point of the current grid cell to get the
    ! connectivity for the current subcell.
    !
    do n = 1,face_subcells
      !
      nc = nc + 1
      !
      ipelem(:,nc) = p1 + int( pst_geom(face_geom,face_order)% &
                                            connectivity(:,n) , &
                               kind=kind(ipelem) )
      !
    end do
    !
    p1 = p1 + face_npts
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_data_for_8s_vorticity_solution
!
!###############################################################################
!
subroutine get_vort_face_coordinates(my_vort_faces,y,z)
  !
  !.. Use Statements ..
  use order_mod,         only : maxSP,maxFP
  use order_mod,         only : geom_solpts
  use geomvar_mod,       only : face,cell,xyz,fp
  use interpolation_mod, only : interp
  use interpolation_mod, only : outerp,outerp_mv_mult
  !
  !.. Formal Arguments ..
  integer,   dimension(:),    intent(in) :: my_vort_faces
  real(CRT), dimension(:), intent(inout) :: y
  real(CRT), dimension(:), intent(inout) :: z
  !
  !.. Local Scalars ..
  integer :: i,j,nf,np,n1,n2,nfp,f1,f2,k,kl
  integer :: this_face,face_geom,face_order
  integer :: host_cell,host_geom,host_order
  real(wp) :: a
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxFP) :: face_y
  real(wp), dimension(1:maxFP) :: face_z
  real(wp), dimension(1:maxSP,1:2) :: txyz
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_vort_face_coordinates"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  f2 = 0
  !
  do nf = 1,size(my_vort_faces)
    !
    this_face = my_vort_faces(nf) ! boundary face of current profile
    !
    face_geom = face(this_face)%geom   ! geometry of boundary face
    face_order = face(this_face)%order ! geometry of boundary face
    !
    host_cell = face(this_face)%left%cell ! host cell on current boundary face
    !
    host_geom = cell(host_cell)%geom   ! host cell geometry type
    host_order = cell(host_cell)%order ! host cell solution order
    !
    n1 = cell(host_cell)%beg_sp ! beginning index for sol-pts in host cell
    n2 = cell(host_cell)%end_sp ! ending    index for sol-pts in host cell
    np = n2-n1+1                ! number of sol-pts in host cell
    !
    nfp = geom_solpts(face_geom,face_order) ! number of sol-pts on boundary face
    f1 = f2 + 1                  ! beginning index for sol-pts on boundary face
    f2 = f2 + nfp                 ! ending    index for sol-pts on boundary face
    !
    ! Create transpose of the xyz(2:3,n1:n2) to reduce cache misses
    !
    txyz(1:np,1:2) = transpose( xyz(2:3,n1:n2) )
    !
    ! Compute the profile quantities at each data point for this face
    !
    do k = 1,nfp
      !
      ! Index of flx-pt local to the current host cell
      !
      kl = face(this_face)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                        rot(face(this_face)%left%rotation)% &
                                            pts(k)
      !
      ! First interpolate to the face flux points
      !
      face_y(k) = zero
      face_z(k) = zero
      do i = interp(host_geom,host_order)%toFace(face_order)%row_ptr(kl)+1, &
             interp(host_geom,host_order)%toFace(face_order)%row_ptr(kl+1)
        !
        j = interp(host_geom,host_order)%toFace(face_order)%row_idx(i)
        a = interp(host_geom,host_order)%toFace(face_order)%row_val(i)
        !
        face_y(k) = face_y(k) + a * txyz(j,1)
        face_z(k) = face_z(k) + a * txyz(j,2)
        !
      end do
      !
    end do
    !
    ! Finally, interpolate the face flux points
    ! to the output face solution points
    !
    if (allocated(outerp(face_geom,face_order)%output)) then
      !
      call outerp_mv_mult( outerp(face_geom,face_order)%output , face_y(1:nfp) )
      call outerp_mv_mult( outerp(face_geom,face_order)%output , face_z(1:nfp) )
      !
    end if
    !
    y(f1:f2) = real( face_y(1:nfp) , kind=kind(y) )
    z(f1:f2) = real( face_z(1:nfp) , kind=kind(z) )
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_vort_face_coordinates
!
!###############################################################################
!
subroutine get_vort_face_solution(my_vort_faces,vorticity)
  !
  !.. Use Statements ..
  use order_mod, only : maxSP,maxFP
  use order_mod, only : geom_solpts
  use geomvar_mod, only : nr,face,cell,fp
  use ovar, only : time_ref
  use interpolation_mod, only : interp
  use interpolation_mod, only : outerp,outerp_mv_mult
  use flowvar_mod, only : usp,dusp
  !
  !.. Formal Arguments ..
  integer,   dimension(:),    intent(in) :: my_vort_faces
  real(CRT), dimension(:), intent(inout) :: vorticity
  !
  !.. Local Scalars ..
  integer :: i,j,k,kl,l,m,nf,n0,nfp,f1,f2
  integer :: this_face,face_geom,face_order
  integer :: host_cell,host_geom,host_order
  real(wp) :: a
  real(wp) :: rho_inv,vx,vy,vz
  real(wp) :: decdx,dmydx,dmzdx
  real(wp) :: decdy,dmxdy,dmzdy
  real(wp) :: decdz,dmxdz,dmydz
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: vort_vec
  real(wp), dimension(1:nq,1:maxFP) :: cvfp
  real(wp), dimension(1:maxFP,1:nq) :: cv_tmp
  real(wp), dimension(1:nr,1:nq,1:maxFP) :: dcvdx
  real(wp), dimension(1:maxFP,1:nr,1:nq) :: dcvdx_tmp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_vort_face_solution"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  f2 = 0
  !
  do nf = 1,size(my_vort_faces)
    !
    this_face = my_vort_faces(nf) ! boundary face of current profile
    !
    face_geom = face(this_face)%geom   ! geometry of boundary face
    face_order = face(this_face)%order ! geometry of boundary face
    !
    host_cell = face(this_face)%left%cell ! host cell on current boundary face
    !
    host_geom = cell(host_cell)%geom   ! host cell geometry type
    host_order = cell(host_cell)%order ! host cell solution order
    !
    n0 = cell(host_cell)%beg_sp-1 ! beginning index for sol-pts in host cell
    !
    nfp = geom_solpts(face_geom,face_order) ! number of sol-pts on boundary face
    f1 = f2 + 1                  ! beginning index for sol-pts on boundary face
    f2 = f2 + nfp                 ! ending    index for sol-pts on boundary face
    !
    ! Compute the profile quantities at each data point for this face
    !
    do k = 1,nfp
      !
      ! Index of flx-pt local to the current host cell
      !
      kl = face(this_face)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                        rot(face(this_face)%left%rotation)% &
                                            pts(k)
      !
      ! First interpolate to the face flux points
      !
      do m = 1,nq
        cv_tmp(k,m) = zero
        dcvdx_tmp(k,1:nr,m) = zero
        do i = interp(host_geom,host_order)%toFace(face_order)%row_ptr(kl)+1, &
               interp(host_geom,host_order)%toFace(face_order)%row_ptr(kl+1)
          !
          j = interp(host_geom,host_order)%toFace(face_order)%row_idx(i)
          a = interp(host_geom,host_order)%toFace(face_order)%row_val(i)
          !
          cv_tmp(k,m) = cv_tmp(k,m) + a * usp(n0+j,m)
          do l = 1,nr
            dcvdx_tmp(k,l,m) = dcvdx_tmp(k,l,m) + a * dusp(n0+j,l,m)
          end do
          !
        end do
      end do
      !
    end do
    !
    ! Now, interpolate the vorticity at the face flux points
    ! to the output face solution points
    !
    if (allocated(outerp(face_geom,face_order)%output)) then
      !
      do m = 1,nq
        call outerp_mv_mult( outerp(face_geom,face_order)%output , &
                             cv_tmp(1:nfp,m) , cvfp(m,1:nfp) )
        do l = 1,nr
          call outerp_mv_mult( outerp(face_geom,face_order)%output , &
                               dcvdx_tmp(1:nfp,l,m) , dcvdx(l,m,1:nfp) )
        end do
      end do
      !
    else
      !
      cvfp(1:nq,1:nfp) = transpose( cv_tmp(1:nfp,1:nq) )
      !
      do k = 1,nfp
        do m = 1,nq
          do l = 1,nr
            dcvdx(l,m,k) = dcvdx_tmp(k,l,m)
          end do
        end do
      end do
      !
    end if
    !
    ! Finally, compute the magnitude of the vorticity vector at each face point
    !
    do k = 1,nfp
      !
      ! Compute the vorticity vector at this face point using the
      ! gradients of the momentums and the product rule
      !
      rho_inv = one/cvfp(nec,k)
      vx = cvfp(nmx,k)*rho_inv
      vy = cvfp(nmy,k)*rho_inv
      vz = cvfp(nmz,k)*rho_inv
      !
      decdx = dcvdx(1,nec,k)
      decdy = dcvdx(2,nec,k)
      decdz = dcvdx(3,nec,k)
      !
      dmxdy = dcvdx(2,nmx,k)
      dmxdz = dcvdx(3,nmx,k)
      !
      dmydx = dcvdx(1,nmy,k)
      dmydz = dcvdx(3,nmy,k)
      !
      dmzdx = dcvdx(1,nmz,k)
      dmzdy = dcvdx(2,nmz,k)
      !
      vort_vec(1) = rho_inv*( (dmzdy-dmydz) - (vz*decdy-vy*decdz) )
      vort_vec(2) = rho_inv*( (dmxdz-dmzdx) - (vx*decdz-vz*decdx) )
      vort_vec(3) = rho_inv*( (dmydx-dmxdy) - (vy*decdx-vx*decdy) )
      !
      ! Compute the magnitude of the vorticity vector
      ! NOTE: Dont forget to adjust the nondimensionalization of the
      !       vorticity by dividing through by time_ref to get vorticity*
      !
      vorticity(f1-1+k) = real( norm2(vort_vec) / time_ref , &
                                kind=kind(vorticity) )
      !
    end do
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_vort_face_solution
!
!###############################################################################
!
pure function create_cell_output_offsets() result(return_value)
  !
  !.. Use Statements ..
  use geomvar_mod,       only : ncell,cell
  use interpolation_mod, only : outerp
  !
  !.. Function Result ..
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: nc,this_geom,this_order,k0,kp
  !
continue
  !
  return_value = init_array(ncell,0) ! F2003 AUTO-REALLOC
  !
  k0 = 0
  kp = 0
  !
  do nc = 1,ncell
    !
    k0 = k0 + kp
    !
    return_value(nc) = k0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    if (allocated(outerp(this_geom,this_order)%output)) then
      kp = size(outerp(this_geom,this_order)%output%mat,dim=1)
    else
      kp = cell(nc)%end_sp - cell(nc)%beg_sp + 1
    end if
    !
  end do
  !
end function create_cell_output_offsets
!
!###############################################################################
!
!pure &
subroutine interpolate_for_output_1d(array,return_value,nondim)
  !
  !.. Use Statements ..
  use geomvar_mod,       only : ncell,cell
  use interpolation_mod, only : outerp
  !
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  use omp_lib, only : omp_get_thread_num
  use omp_lib, only : omp_get_num_threads
#endif
#endif
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: array(:)
  real(wp), intent(inout) :: return_value(:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: nondim
  !
  !.. Local Scalars ..
  integer  :: i,k,n,k0,kp,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interpolate_for_output_1d"
  !
continue
  !
  if (.not. allocated(cell_k0)) then
    cell_k0 = create_cell_output_offsets() ! F2003 AUTO-REALLOC
  end if
! !
  ! Interpolate the data
  !
  k0 = 0
  kp = 0
  !
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP PARALLEL &
  !$OMP&         DEFAULT(SHARED) &
  !$OMP&         PRIVATE(nc,this_geom,this_order,n0,np,k0,kp,k,val,i,n,a,tid)
  !
  tid = omp_get_thread_num()
  nthreads = omp_get_num_threads()
  !
  if (mypnum == glb_root) then
    if (tid == glb_root) then
      write (iout,thread_report) pname,nthreads
    end if
  end if
  !
  chunk_size = ceiling( real(ncell,kind=wp) / real(nthreads,kind=wp) )
  chunk_size = chunk_size / 10
  !
  !$OMP DO
 !!$OMP DO SCHEDULE(RUNTIME)
 !!$OMP DO SCHEDULE(GUIDED,chunk_size)
 !!$OMP DO SCHEDULE(STATIC,chunk_size)
#endif
#endif
  do nc = 1,ncell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
   !k0 = k0 + kp
    k0 = cell_k0(nc)
    !
    if (allocated(outerp(this_geom,this_order)%output)) then
      !
      kp = size(outerp(this_geom,this_order)%output%mat,dim=1)
      !
      do k = 1,kp
        !
        val = zero
        !
        do i = outerp(this_geom,this_order)%output%col_ptr(k)+1, &
               outerp(this_geom,this_order)%output%col_ptr(k+1)
          !
          n = outerp(this_geom,this_order)%output%col_idx(i)
          a = outerp(this_geom,this_order)%output%col_val(i)
          !
          val = val + a * array(n0+n)
          !
        end do
        !
        return_value(k0+k) = val
        !
      end do
      !
    else
      !
      kp = np
      !
      do n = 1,np
        return_value(n0+n) = array(n0+n)
      end do
      !
    end if
    !
  end do
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP END DO
#endif
#endif
  !
  ! Dimensionalize the variable if the nondimensionalization factor is present
  !
  if (present(nondim)) then
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
    !$OMP DO
   !!$OMP DO SCHEDULE(RUNTIME)
   !!$OMP DO SCHEDULE(GUIDED,chunk_size)
   !!$OMP DO SCHEDULE(STATIC,chunk_size)
#endif
#endif
    do n = 1,size(return_value)
      return_value(n) = return_value(n) * nondim
    end do
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
    !$OMP END DO
#endif
#endif
  end if
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP END PARALLEL
#endif
#endif
  !
end subroutine interpolate_for_output_1d
!
!###############################################################################
!
!pure &
subroutine interpolate_for_output_2d(array,return_value,nondim)
  !
  !.. Use Statements ..
  use geomvar_mod,       only : ncell,cell
  use interpolation_mod, only : outerp
  !
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  use omp_lib, only : omp_get_thread_num
  use omp_lib, only : omp_get_num_threads
#endif
#endif
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: array(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: nondim
  !
  !.. Local Scalars ..
  integer  :: i,k,m,n,k0,kp,nc,np,n0,nvar
  integer  :: this_geom,this_order
  real(wp) :: a,val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interpolate_for_output_2d"
  !
continue
  !
  if (.not. allocated(cell_k0)) then
    cell_k0 = create_cell_output_offsets() ! F2003 AUTO-REALLOC
  end if
  !
  nvar = size(array,dim=2)
  !
  ! Interpolate the data
  !
  k0 = 0
  kp = 0
  !
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP PARALLEL &
  !$OMP&         DEFAULT(SHARED) &
  !$OMP&         PRIVATE(nc,this_geom,this_order,n0,np,k0,kp,k,val,i,n,m,a,tid)
  !
  tid = omp_get_thread_num()
  nthreads = omp_get_num_threads()
  !
  if (mypnum == glb_root) then
    if (tid == glb_root) then
      write (iout,thread_report) pname,nthreads
    end if
  end if
  !
  chunk_size = ceiling( real(ncell,kind=wp) / real(nthreads,kind=wp) )
  chunk_size = chunk_size / 10
  !
  !$OMP DO
 !!$OMP DO SCHEDULE(RUNTIME)
 !!$OMP DO SCHEDULE(GUIDED,chunk_size)
 !!$OMP DO SCHEDULE(STATIC,chunk_size)
#endif
#endif
  do nc = 1,ncell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
   !k0 = k0 + kp
    k0 = cell_k0(nc)
    !
    if (allocated(outerp(this_geom,this_order)%output)) then
      !
      kp = size(outerp(this_geom,this_order)%output%mat,dim=1)
      !
      do m = 1,nvar
        do k = 1,kp
          !
          val = zero
          !
          do i = outerp(this_geom,this_order)%output%col_ptr(k)+1, &
                 outerp(this_geom,this_order)%output%col_ptr(k+1)
            !
            n = outerp(this_geom,this_order)%output%col_idx(i)
            a = outerp(this_geom,this_order)%output%col_val(i)
            !
            val = val + a * array(n0+n,m)
            !
          end do
          !
          return_value(k0+k,m) = val
          !
        end do
      end do
      !
    else
      !
      kp = np
      !
      do m = 1,nvar
        do n = 1,np
          return_value(k0+n,m) = array(n0+n,m)
        end do
      end do
      !
    end if
    !
  end do
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP END DO
#endif
#endif
  !
  ! Dimensionalize the variables if the nondimensionalization factor is present
  !
  if (present(nondim)) then
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
    !$OMP DO
   !!$OMP DO SCHEDULE(RUNTIME)
   !!$OMP DO SCHEDULE(GUIDED,chunk_size)
   !!$OMP DO SCHEDULE(STATIC,chunk_size)
#endif
#endif
    do m = 1,size(return_value,dim=2)
      do n = 1,size(return_value,dim=1)
        return_value(n,m) = return_value(n,m) * nondim
      end do
    end do
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
    !$OMP END DO
#endif
#endif
  end if
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP END PARALLEL
#endif
#endif
  !
end subroutine interpolate_for_output_2d
!
!###############################################################################
!
!pure &
subroutine interpolate_for_output_3d(array,return_value,nondim)
  !
  !.. Use Statements ..
  use geomvar_mod,       only : ncell,cell
  use interpolation_mod, only : outerp
  !
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  use omp_lib, only : omp_get_thread_num
  use omp_lib, only : omp_get_num_threads
#endif
#endif
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: array(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: nondim
  !
  !.. Local Scalars ..
  integer  :: i,l,k,m,n,k0,kp,nc,np,n0,nvar,ndim
  integer  :: this_geom,this_order
  real(wp) :: a,val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interpolate_for_output_3d"
  !
continue
  !
  if (.not. allocated(cell_k0)) then
    cell_k0 = create_cell_output_offsets() ! F2003 AUTO-REALLOC
  end if
  !
  ndim = size(array,dim=2)
  nvar = size(array,dim=3)
  !
  ! Interpolate the data
  !
  k0 = 0
  kp = 0
  !
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP PARALLEL &
  !$OMP&         DEFAULT(SHARED) &
  !$OMP&         PRIVATE(nc,this_geom,this_order,n0,np,k0,kp,k, &
  !$OMP&                 val,i,n,m,l,a,tid)
  !
  tid = omp_get_thread_num()
  nthreads = omp_get_num_threads()
  !
  if (mypnum == glb_root) then
    if (tid == glb_root) then
      write (iout,thread_report) pname,nthreads
    end if
  end if
  !
  chunk_size = ceiling( real(ncell,kind=wp) / real(nthreads,kind=wp) )
  chunk_size = chunk_size / 10
  !
  !$OMP DO
 !!$OMP DO SCHEDULE(RUNTIME)
 !!$OMP DO SCHEDULE(GUIDED,chunk_size)
 !!$OMP DO SCHEDULE(STATIC,chunk_size)
#endif
#endif
  do nc = 1,ncell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
   !k0 = k0 + kp
    k0 = cell_k0(nc)
    !
    if (allocated(outerp(this_geom,this_order)%output)) then
      !
      kp = size(outerp(this_geom,this_order)%output%mat,dim=1)
      !
      do m = 1,nvar
        do l = 1,ndim
          do k = 1,kp
            !
            val = zero
            !
            do i = outerp(this_geom,this_order)%output%col_ptr(k)+1, &
                   outerp(this_geom,this_order)%output%col_ptr(k+1)
              !
              n = outerp(this_geom,this_order)%output%col_idx(i)
              a = outerp(this_geom,this_order)%output%col_val(i)
              !
              val = val + a * array(n0+n,l,m)
              !
            end do
            !
            return_value(k0+k,l,m) = val
            !
          end do
        end do
      end do
      !
    else
      !
      kp = np
      !
      do m = 1,nvar
        do l = 1,ndim
          do n = 1,np
            return_value(k0+n,l,m) = array(n0+n,l,m)
          end do
        end do
      end do
      !
    end if
    !
  end do
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP END DO
#endif
#endif
  !
  ! Dimensionalize the variables if the nondimensionalization factor is present
  !
  if (present(nondim)) then
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
    !$OMP DO
   !!$OMP DO SCHEDULE(RUNTIME)
   !!$OMP DO SCHEDULE(GUIDED,chunk_size)
   !!$OMP DO SCHEDULE(STATIC,chunk_size)
#endif
#endif
    do m = 1,size(return_value,dim=3)
      do l = 1,size(return_value,dim=2)
        do n = 1,size(return_value,dim=1)
          return_value(n,l,m) = return_value(n,l,m) * nondim
        end do
      end do
    end do
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
    !$OMP END DO
#endif
#endif
  end if
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP END PARALLEL
#endif
#endif
  !
end subroutine interpolate_for_output_3d
!
!###############################################################################
!
subroutine compute_dmatrix_grad_temperature(dtdx)
  !
  !.. Use Statements
  use geomvar_mod,       only : nr,n_solpts,ncell,cell
  use order_mod,         only : n_order,t_order,geom_solpts
  use metrics_mod,       only : metrics_array
  use interpolation_mod, only : interp
  use derivatives_mod,   only : deriv
  use projection_mod,    only : projct
  use flowvar_mod,       only : usp,gpsp
  !
  use thermo_mod, only : update_thermo_properties
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: dtdx(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,l,m,n,nc,np,n0
  integer  :: this_geom,this_order
  integer  :: nqp,num_qpts,q0,ierr
  real(wp) :: a,t,val
  !
  !.. Local Arrays ..
  real(wp) :: dtdi(1:3)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: cv_qp(:,:)
  real(wp), allocatable :: gp_qp(:,:)
  real(wp), allocatable :: temp(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_dmatrix_grad_temperature"
  !
continue
  !
  allocate ( temp(1:n_solpts) , source=zero , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp",1,__LINE__,__FILE__,ierr,error_message)
  !
  if (t_order <= n_order) then
    !
    do n = 1,n_solpts
      temp(n) = gpsp(n,ntemp)
    end do
    !
  else
    !
    ! 1. Interpolate conserved variables to quadrature points and compute
    !    temperature at quadrature points
    !
    nqp = geom_solpts(cell(1)%geom,t_order)
    num_qpts = ncell*nqp
    !
    allocate ( cv_qp(1:num_qpts,1:nq) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cv_qp",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( gp_qp(1:num_qpts,1:ng) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"gp_qp",1,__LINE__,__FILE__,ierr,error_message)
    !
    !
    !
    q0 = 0
    !
    do nc = 1,ncell
      !
      n0 = cell(nc)%beg_sp - 1
      !
      this_geom = cell(nc)%geom   ! geometry type of this cell
      this_order = cell(nc)%order ! solution order of this cell
      !
      ! Now interpolate the conserved variables from the solution points
      ! to the quadrature points
      !
      do n = 1,nqp
        do m = 1,nq
          !
          val = zero
          !
          do i = interp(this_geom,this_order)%toCell(t_order)%row_ptr(n)+1, &
                 interp(this_geom,this_order)%toCell(t_order)%row_ptr(n+1)
            !
            j = interp(this_geom,this_order)%toCell(t_order)%row_idx(i)
            a = interp(this_geom,this_order)%toCell(t_order)%row_val(i)
            !
            val = val + a * usp(n0+j,m)
            !
          end do
          !
          cv_qp(q0+n,m) = val
          !
        end do
        !
      end do
      !
      q0 = q0 + nqp
      !
    end do
    !
    call update_thermo_properties(1,num_qpts,cv_qp,gp_qp,caller_name=pname)
    !
    if (allocated(cv_qp)) then
      deallocate ( cv_qp , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"cv_qp",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    ! 2. Project temperature at quadrature points to solution order polynomial
    !    at solution points
    !
    q0 = 0
    !
    do nc = 1,ncell
      !
      n0 = cell(nc)%beg_sp - 1
      np = cell(nc)%end_sp - n0
      !
      this_geom = cell(nc)%geom   ! geometry type of this cell
      this_order = cell(nc)%order ! solution order of this cell
      !
      ! Now interpolate the derivatives of the conserved variables
      ! from the solution points to the quadrature points
      !
      do n = 1,np
        !
        t = zero
        !
        do i = projct(this_geom,this_order)%CellOI(t_order)%row_ptr(n)+1, &
               projct(this_geom,this_order)%CellOI(t_order)%row_ptr(n+1)
          !
          j = projct(this_geom,this_order)%CellOI(t_order)%row_idx(i)
          a = projct(this_geom,this_order)%CellOI(t_order)%row_val(i)
          !
          t = t + a * gp_qp(q0+j,ntemp)
          !
        end do
        !
        temp(n0+n) = t
        !
      end do
      !
      q0 = q0 + nqp
      !
    end do
    !
    if (allocated(gp_qp)) then
      deallocate ( gp_qp , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"gp_qp",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
  end if
  !
  ! 4. Use derivative matrices to compute gradient of temperature 
  !    in computation space
  !
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1  ! beginning index for solpts in cell nc
    np = cell(nc)%end_sp-n0 ! number of solution points in this cell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Compute the gradient of the conservative variables at each solution point
    !
    do n = 1,np
      !
      do l = 1,nr
        !
        dtdi(l) = zero
        !
        do i = deriv(this_geom,this_order)%d(l)%row_ptr(n)+1, &
               deriv(this_geom,this_order)%d(l)%row_ptr(n+1)
          !
          j = deriv(this_geom,this_order)%d(l)%row_idx(i)
          a = deriv(this_geom,this_order)%d(l)%row_val(i)
          !
          dtdi(l) = dtdi(l) + a * temp(n0+j)
          !
        end do
        !
      end do
      !
      do l = 1,nr
        dtdx(n0+n,l) = dot( dtdi(1:nr) , metrics_array(n0+n,1:nr,l) )
      end do
      !
    end do
    !
  end do
  !
  if (allocated(temp)) then
    deallocate ( temp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"temp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
end subroutine compute_dmatrix_grad_temperature
!
!###############################################################################
!
end module write_cgns_mod
!
!###############################################################################
!
