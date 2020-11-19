module postprocessor_mod
  !
#ifdef DEBUG_ON
!!#define DEBUG_POST_PROCESSOR
#endif
!#define VERBOSE_MPI
  !
  !.. Global Use Statements ..
  use module_kind_types, global_debug_timer => debug_timer
  use module_cgns_types
  use eqn_idx, only : nec,nmb,nme,nmx,nmy,nmz,nee,nq
  use eqn_idx, only : ncs,npres,ntemp,ngamm,nenth,ncp,nvisc,ncond,ng
  !
  implicit none
  !
  private
  !
  ! ##################################
  ! ###  Module Public Procedures  ###
  ! ##################################
  !
  public :: new_write_parallel_cgns
  public :: convert_restart_files_to_continuous_cgns
  !
  !
  ! ################################
  ! ###  Module Local Variables  ###
  ! ################################
  !
  type :: geom_pnt_map_t
    integer, allocatable :: sp(:)
    integer, allocatable :: fp(:)
    integer, allocatable :: ep(:)
    integer, allocatable :: np(:)
  end type geom_pnt_map_t
  !
  type :: output_face_pts_t
    integer              :: face_idx
    integer, allocatable :: out_idx(:)
  end type output_face_pts_t
  !
  type :: output_edge_pts_t
    integer              :: edge_idx
    integer, allocatable :: out_idx(:)
  end type output_edge_pts_t
  !
  type :: output_node_pts_t
    integer :: node_idx
    integer :: out_idx
  end type output_node_pts_t
  !
  integer,                 allocatable :: output_sol_pts(:)
  !
  type(output_face_pts_t), allocatable :: output_face_pts(:)
  !
  type(output_edge_pts_t), allocatable :: output_edge_pts(:)
  !
  type(output_node_pts_t), allocatable :: output_node_pts(:)
  !
  type :: plot_func_t
    integer :: func
    character(len=CGLEN) :: varname = ""
  contains
    procedure :: get_function_name
  end type plot_func_t
  !
  type(plot_func_t), public, parameter :: &
              plot_static_density     =  plot_func_t( 1),    &
              plot_velocity_x         =  plot_func_t( 2),    &
              plot_velocity_y         =  plot_func_t( 3),    &
              plot_velocity_z         =  plot_func_t( 4),    &
              plot_static_pressure    =  plot_func_t( 5),    &
              plot_static_temperature =  plot_func_t( 6),    &
              plot_momentum_x         =  plot_func_t( 7),    &
              plot_momentum_y         =  plot_func_t( 8),    &
              plot_momentum_z         =  plot_func_t( 9),    &
              plot_total_energy       =  plot_func_t(10),    &
              plot_entropy            =  plot_func_t(11),    &
              plot_mach_number        =  plot_func_t(12),    &
              plot_vorticity_x        =  plot_func_t(13),    &
              plot_vorticity_y        =  plot_func_t(14),    &
              plot_vorticity_z        =  plot_func_t(15),    &
              plot_grad_density(3)    = [plot_func_t(16), &
                                         plot_func_t(17), &
                                         plot_func_t(18)],   &
              plot_grad_momentum_x(3) = [plot_func_t(19), &
                                         plot_func_t(20), &
                                         plot_func_t(21)],   &
              plot_grad_momentum_y(3) = [plot_func_t(22), &
                                         plot_func_t(23), &
                                         plot_func_t(24)],   &
              plot_grad_momentum_z(3) = [plot_func_t(25), &
                                         plot_func_t(26), &
                                         plot_func_t(27)],   &
              plot_grad_energy(3)     = [plot_func_t(28), &
                                         plot_func_t(29), &
                                         plot_func_t(30)]
  !
  type(plot_func_t), allocatable :: plot_var(:)
  !
  integer, save, allocatable :: ipelem_of_cells(:,:,:)
  !
  ! INITIALIZATION_NEEDED: Private module variable that indicates whether the
  !                        initializations required for CGNS have been completed
  !                        or not.
  !
  logical(lk), save :: initialization_needed = true
  !
  character(len=*), parameter :: char_xyz(1:3) = ["X","Y","Z"]
  !
  ! ###########################
  ! ###  Module Parameters  ###
  ! ###########################
  !
  logical(lk), parameter :: check_gradients = fals
 !logical(lk), parameter :: check_gradients = true
  !
  ! ################################
  ! ###  Module Local Variables  ###
  ! ################################
  !
  !.. Output File Names ..
  character(len=*),      parameter :: cgns_grid_fname = "solution.grid.cgns"
  character(len=*),      parameter :: cgns_time_fname = "solution.time.cgns"
  character(len=170),         save :: cgns_dir = ""
  character(len=170),         save :: cgns_grid_file
  !
  character(len=170),         save :: cgnsfile_fmt
  !
  !.. CGNS Path Variables ..
  integer(CBT),      save :: ifile_n
  integer(CBT),      save :: ibase_n
  integer(CBT),      save :: izone_n
  integer(CBT),      save :: isect_n
  integer(CBT), parameter :: iiter_n = 1_CBT
  !
  !.. CGNS Time Dependent Variables ..
  integer, save :: cgns_zone_counter
  !
  real(sp),             save, allocatable :: soltime(:)
  character(len=CGLEN), save, allocatable :: solname(:)
  character(len=CGLEN), save, allocatable :: solfile(:)
  !
  character(len=*), parameter :: base_name = "gfr_base"
  character(len=*), parameter :: zone_name = "Solution"
  !
  !.. CGNS Parallel Point Indices ..
  integer(CST), save :: n_my_output_pts
  integer(CST), save :: my_pts_beg
  integer(CST), save :: my_pts_end
  !
  integer(CST), dimension(1:3) :: zone_size
  !
  ! DISABLE_CGNS_OUTPUT: Private module variable to disable CGNS output. It
  !                      is currently only relevant to the Taylor-Green Vortex
  !                      problem and is more of a temporary addition until the
  !                      output logic/code becomes more robust.
  !
  logical(lk), save :: disable_cgns_output = fals
  !
!!" FINISHED SUBROUTINE MAP_TO_OUTPUT_POINTS AND THE RELATED FUNCTIONS FOR"
!!" HEXAHEDRA CELLS. NOW NEED TO GO BACK AND FIGURE OUT HOW TO APPLY THIS"
!!" FOR THE CONNECTIVITY AND FOR CREATING THE SOLUTION ARRAY SENT TO CGNS."
!!!
!!" WHEN CREATING THE SOLUTION ARRAY, NEED TO FIGURE OUT WHAT TO DO ABOUT"
!!" POINTS FOR WHICH AN ADJACENT PARTITION IS RESPONSIBLE. IN THESE CASES"
!!" THERE IS A POINT IN THE SOLUTION ARRAY THAT NEEDS TO BE SKIPPED FOR A"
!!" PARTITION THAT SHARES THIS POINT BUT IS NOT RESPONSIBLE FOR ITS SOLUTION."
!!!
!!!
!!" NEED TO FIGURE OUT CONNECTIVITY OF THE GLOBAL GRID AND WHICH CELLS ARE"
!!" RESPONSIBLE FOR WHICH OUTPUT POINTS."
!!"  - NEED TO FIGURE THIS OUT FOR BOTH SERIAL AND PARALLEL CASES."
!!"  - USE CELL(:)%NODE(:)%CELL_IDX, CELL(:)%EDGE(:)%CELL_IDX(:), AND"
!!"    CELL(:)%FACE(:)%CELL_IDX(:) TO MAP THE LOCATIONS OF THE NODE, EDGE,"
!!"    AND FACE VALUES WITHIN NODEUSP, NODEDUSP, EDGEUSP, EDGEDUSP, FACEUSP,"
!!"    AND FACEDUSP TO THE CORRECT LOCATIONS WITHIN THE ARRAY THATS CREATED"
!!"    TO COLLECT THE SOLUTIONS FOR SENDING TO CGNS."
!!"  - TO REMOVE THE DUPLICATE POINT INDICES ON PARTITION BOUNDARIES IN THE"
!!"    CGNS FILE, A POSSIBLE PATH COULD BE TO USE THE NODEUSP, EDGEUSP, AND"
!!"    FACEUSP ARRAYS TO STORE A REAL VALUE OF A GIVEN POINT INDEX ON THE"
!!"    PROCESSOR THAT IS RESPONSIBLE FOR WRITING THE SOLUTION AT THIS POINT,"
!!"    AND ALL OTHER PROCESSORS STORE A NEGATIVE VALUE (OR SOMETHING THAT CAN"
!!"    BE IDENTIFIED AS A POINT INDEX FROM THE WRONG PROCESSOR) AND USE THE"
!!"    EXISTING COMMUNICATION/EXCHANGE CODE TO SYNCRONIZE THE CORRECT POINT"
!!"    INDEX ON ALL PROCESSORS. THERE COULD BE SOME LOGIC THAT ONLY ADDS"
!!"    POSITIVE VALUES FROM THE RECEIVE BUFFER BACK INTO NODEUSP, EDGEUSP,"
!!"    OR FACEUSP. AFTER THAT, ITS JUST A MATTER OF CONVERTING THE REAL"
!!"    VALUES TO INTEGERS AND USE THEM TO OVERWRITE THE PROPER POINT INDICES"
!!"    IN THE CONNECTIVITY ARRAY THAT SENT INTO CGNS. JUST MAKE SURE NOT"
!!"    TO INCLUDE ANY AVERAGING COEFFICIENTS."
!!!
!!" NEED TO FIGURE OUT HOW TO OVERRIDE THE INTERIOR SOLUTIONS FOR THE"
!!" SOLUTIONS ON BOUNDARY FACES, EDGES, AND NODES."
!!"   - NEED TO ADJUST THE LOCAL CELLS_WITH_NODE AND CELLS_WITH_EDGE ARRAYS"
!!"     ON EACH PARTITION SO THAT IF A NODE OR AN EDGE ARE ON A BOUNDARY"
!!"     FACES THE ONLY CELLS REMAINING WITHIN THESE ARRAYS ARE THE GHOST"
!!"     CELLS/FACES THAT CONTRIBUTE TO THE SOLUTION AT THE NODE/EDGE."
!!"       - NEED TO TAKE INTO ACCOUNT LOCAL BC PRIORITIES AS WELL."
!!"   - NEED TO CREATE SOMETHING TO INTERPOLATE THE BOUNDARY FACE SOLUTION"
!!"     FOR THE GHOST CELLS/FACES TO THE NODES/EDGES FOR THE BOUNDARY FACES."
!!"     INTERPOLATING THE INTERIOR SOLUTION TO THOSE NODE/EDGE POINTS WOULD"
!!"     NOT BE CORRECT, FOR EXAMPLE INTERPOLATING THE HOST CELL SOLUTION"
!!"     TO A BOUNDARY FACE THAT HAS A NO-SLIP WALL BOUNDARY CONDITION."
!!"   - THIS WILL PROBABLY NEED TO BE DONE BEFORE EXCHANGING SINCE THE"
!!"     SOLUTION FROM AN ADJACENT PARTITION MIGHT OVERRIDE THE LOCAL"
!!"     SOLUTION OR VICE-VERSA."
!!"   + IF THIS IS DONE BEFORE THE EXCHANGE, WE NEED TO MAKE SURE THE"
!!"     AVERAGING CORRECTION IS BASED ON THE NUMBER OF BOUNDARY FACES"
!!"     THAT CONTAIN THE NODE/EDGE AND NOT THE NUMBER OF INTERNAL GRID"
!!"     CELLS."
!!!
!!" NEED TO FIGURE OUT HOW TO CONSTRUCT THE GLOBAL EDGE AND NODE"
!!" CONNECTIVITY SO THAT WE CAN CREATE THE MPI DATA TYPES FOR EXCHANGING"
!!" EDGE AND NODE SOLUTIONS."
!!"   - NEED TO CREATE THE EDGE DERIVED TYPE ARRAY"
!!"       + MAKE SURE THAT THE INDICES IN EDGE(:)%CPU(:)%EDGPTS STARTS AT 0"
!!"       - MIGHT NEED TO FIX THE ORDERING OF THE EDGE POINT INDICES THAT"
!!"         ARE STORED IN EDGE(:)%CPU(:)%EDGPTS. CURRENTLY IT JUST ASSUMES"
!!"         THE ORDERING IS [0:NUMBER OF EDGE POINTS-1]."
!!"   + NEED TO CREATE THE NODE DERIVED TYPE ARRAY"
!!"   + NEED TO CREATE THE EDGE_MAP ARRAY"
!!"   + NEED TO CREATE THE EXCH_NEUSP ARRAY"
!!"   - "
!!"   - "
!!!
!!" NEED TO FIGURE OUT HOW TO COLLECT THE INTERIOR, FACE, EDGE, AND NODE"
!!" SOLUTIONS TOGETHER AS ONE SOLUTION FOR A CELL BEFORE WRITING IT TO"
!!" THE CGNS SOLUTION FILE."
!!"   - USE THE SAME SUBCONNECTIVITY FROM BEFORE INCLUDING THE ADD POINTS"
!!"     VALUE DEPENDING ON THE CONTINUOUS SOLUTION INPUT PARAMETER. WE JUST"
!!"     NEED A SIMPLE INTEGER ARRAY THAT MAPS THE VARIOUS FACE, EDGE, AND"
!!"     NODE POINTS TO THE CORRECT LOCATIONS WITHIN THE ARRAY THAT WILL"
!!"     STORE THE SOLUTION FOR THE CELL."
!!!
!!" NEED TO DEALLOCATE THE EDGEXYZ AND FACEXYZ (EXCEPT WHEN USING MMS)"
!!" ARRAYS AFTER WRITING THE FIRST GRID CGNS FILE SINCE THEY WILL NO"
!!" LONGER BE NEEDED."
!!!
!!" THE LOCAL AVERAGING COEFFICIENTS SHOULD ACCOUNT ONLY FOR LOCAL CELLS."
!!" THE BND_NODES AND BND_EDGES ARRAYS WILL CORRECT THE LOCAL AVERAGING"
!!" COEFFICIENTS TO THE CORRECT GLOBAL AVERAGING COEFFICIENT."
  !
  !
  interface operator(==)
    module procedure plot_func_type_equals
  end interface
  !
contains
!
!###############################################################################
!
subroutine convert_restart_files_to_continuous_cgns()
  !
  !.. Use Statements ..
  use ovar,           only : restart_list_file,restart_file
  use ovar,           only : output_time_averaging,read_time_ave_restart
  use ovar,           only : ave_start_time,prev_ave_time
  use ovar,           only : time,itcur,tmrst,itrst
  use ovar,           only : time_ref,time_scaling_factor
  use restart_mod,    only : read_restart_file
  !
  !.. Local Scalars ..
  integer  :: ierr,io,n
  real(wp) :: out_time
  !
  character(len=150) :: fname
  character(len=150) :: text
  !
  !.. Local Allocatable Arrays ..
  character(len=100), allocatable :: restart_file_list(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: &
                              pname = "convert_restart_files_to_continuous_cgns"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Write out the CGNS grid file and first solution file using
  ! using the restart file provided in the input file
  !
  if (mypnum == glb_root) then
    out_time = time*time_ref*time_scaling_factor
    if (output_time_averaging) then
      out_time = out_time - ave_start_time*time_ref*time_scaling_factor
    end if
    write (iout,3) out_time
  end if
  !
  call new_write_parallel_cgns(open_CGNS)
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
  fname = trim(adjustl(restart_list_file))
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
  ! Loop over the list of files to read each restart
  ! file and write out a CGNS solution file
  !
  do n = 1,size(restart_file_list)
    !
    restart_file = trim(adjustl(restart_file_list(n)))
    !
    if (mypnum == glb_root) then
      write (iout,2) trim(adjustl(restart_file))
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
    if (n == size(restart_file_list)) then
      call new_write_parallel_cgns(close_CGNS)
    else
      call new_write_parallel_cgns(add_to_CGNS)
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
  2 format (/,4x,"Working on restart file '",a,"'",/, &
              8x,"Reading the restart file")
  3 format (/,2x,"Creating CGNS solution file from restart file at ", &
                 "tstar = ",f10.4,/)
  4 format ("STOPPING EXECUTION AFTER CONVERTING THE LIST OF ", &
            "RESTART FILES TO CGNS SOLUTION FILES!")
  5 format ("STOPPING EXECUTION AFTER CONVERTING THE RESTART ", &
            "FILE TO A CGNS SOLUTION FILE!")
  !
  !
end subroutine convert_restart_files_to_continuous_cgns
!
!###############################################################################
!
subroutine new_write_parallel_cgns(iflag)
  !
  !.. Use Statements ..
  use ovar, only : itcur,time,time_ref
  use ovar, only : write_TaylorGreen_full_solution
  use ovar, only : completely_disable_cgns
  use ovar, only : output_time_averaging
  use ovar, only : time_scaling_factor
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iflag
  !
  !.. Local Scalars ..
  integer :: n,ierr
  !
  character(len=CGLEN) :: tname

  character(len=170)   :: sol_file
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "new_write_parallel_cgns"
  !
continue
  !
  if (completely_disable_cgns) return
  !
  if (disable_cgns_output) then
    if (.not. initialization_needed) then
      if (.not. write_TaylorGreen_full_solution) then
        return
      end if
    end if
  end if
  !
  call debug_timer(entering_procedure,pname)
  !
  if (initialization_needed) then
    !
    ! Check for the CGNS output directory and create it if not found
    !
    call check_cgns_dir




    !
    ! Write the CGNS grid file
    !
    call write_grid_parallel_cgns
    !
    ! Initialize the zone counter
    !
    cgns_zone_counter = 0
    !
    ! Now that the CGNS grid file has been successfully created, set the
    ! initialization flag to false so this isnt done again.
    !
    initialization_needed = fals
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
  call reallocate(solname,cgns_zone_counter+1)
  call reallocate(solfile,cgns_zone_counter+1)
  soltime(size(soltime)) = real(time*time_ref*time_scaling_factor,kind(soltime))
  solname(size(solname)) = tname
  solfile(size(solfile)) = sol_file
  !
  ! Write the solution file
  !
  call write_solution_parallel_cgns(sol_file,tname)
  !
  ! Write the solution file
  !
  if (output_time_averaging) then
 !  call write_time_ave_parallel_cgns
  end if
  !
  ! ################################################
  ! #####   FINISHED WRITING FIELD VARIABLES   #####
  ! ################################################
  !
  call set_cgns_queue(STOP_CGNS_QUEUE,pname,cgns_grid_file)
  !
  ! Write the time iterative solution file
  !
  if (any(iflag == [close_CGNS,-3])) then
 !  call write_time_file_parallel_cgns
  end if
  !
  ! Update the zone counter
  !
  cgns_zone_counter = cgns_zone_counter + 1
  !
  call debug_timer(leaving_procedure,pname)
  !
  write (error_message,100)
  100 format ("Temporary stop after writing initial CGNS files using the ", &
              "new continuous output from the post-processor module.")
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  !
  ! Format Statements
  !
  1 format ("#############################################",/, &
            "    FOR PROCESSOR #",i0)
  2 format (a," = ",i0)


  !
end subroutine new_write_parallel_cgns
!
!###############################################################################
!
subroutine write_grid_parallel_cgns()
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : profile_io_cgns,profile_io_cgns_subs
  use ovar,   only : grid_scaling_factor
  use ovar,   only : channel_nondim_output_coord
  !
  !.. Local Scalars ..
  integer :: i,l,k,n,ierr
  !
  integer(CST) :: my_cells_beg,my_cells_end
  integer(CST) :: n_total_output_cells
  !
  real(wp) :: beg_wall_time
  real(wp) :: end_wall_time
  real(wp) :: cgns_wall_time
  !
  character(len=CGLEN) :: varname
  !
  !.. Local Arrays ..
  integer(CBT) :: ixyz_n(1:nr)
  !
  !.. Local Allocatable Arrays ..
  integer(CST), allocatable :: ipelem(:,:)
  !
  real(CRT), allocatable :: var(:)
  real(wp),  allocatable :: vartmp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_grid_parallel_cgns"
  character(len=*), parameter :: char_xyz(1:3) = ["X","Y","Z"]
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
#ifndef SPECIAL_FOR_PGI
  !
  if (profile_io_cgns) then
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum,"before mpi_barrier",__LINE__,__FILE__
#endif
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum," after mpi_barrier",__LINE__,__FILE__
#endif
    if (mypnum == glb_root) then
      write (iout,1) pname
      flush (iout)
      beg_wall_time = mpi_wtime()
    end if
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum,"before mpi_barrier",__LINE__,__FILE__
#endif
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum," after mpi_barrier",__LINE__,__FILE__
#endif
  end if
  !
  cgns_grid_file = trim(adjustl(cgns_dir)) // cgns_grid_fname
  !
  ! Initialize the variables and data arrays needed for CGNS output
  ! NOTE: This call to routine create_cell_subconnectivity will initialize
  !       the module variables my_pts_beg, my_pts_end, and n_my_output_pts.
  !
  call create_cell_subconnectivity(my_cells_beg,my_cells_end,ipelem)
  !
  ! Extract the total number of cells from the zone_size array
  !
  n_total_output_cells = zone_size(2)
  !
  ! Open the cgns grid file
  !
  call cgns_timer(start_cgns_file_timer,profile_io_cgns_subs,cgns_grid_file)
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cgp_open_f")
  call cgp_open_f(cgns_grid_file,CG_MODE_WRITE,ifile_n,cgierr)
  call cgns_error(pname,cgns_grid_file,"cgp_open_f", &
                  cgierr,__LINE__,__FILE__, &
                  "Mode = CG_MODE_WRITE")
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cgp_open_f")
  !
  ! Set the CGNS base (who knows what this is for)
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_base_write_f")
  call cg_base_write_f(ifile_n,base_name,int(nr,kind=CBT), &
                       int(nr,kind=CBT),ibase_n,cgierr)
  call cgns_error(pname,cgns_grid_file,"cg_base_write_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_base_write_f")
  !
  ! Write the title for this solution
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_zone_write_f")
  call cg_zone_write_f(ifile_n,ibase_n,zone_name,zone_size, &
                       UNSTRUCTURED,izone_n,cgierr)
  call cgns_error(pname,cgns_grid_file,"cg_zone_write_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_zone_write_f")
  !
  ! Create data nodes for the coordinates
  !
  do n = 1,nr
    varname = "Coordinate"//char_xyz(n)
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_coord_write_f")
    call cgp_coord_write_f(ifile_n,ibase_n,izone_n,CGNS_KIND_REAL, &
                           trim(varname),ixyz_n(n),cgierr)
    call cgns_error(pname,cgns_grid_file,"cgp_coord_write_f", &
                    cgierr,__LINE__,__FILE__,"Write node: Coordinate",n)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_coord_write_f")
  end do
  !
  ! Now create and write the coordinate data in parallel
  !
  ! Allocate the temporary variable arrays
  !
  allocate ( vartmp(1:nr,1:n_my_output_pts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"vartmp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( var(1:n_my_output_pts) , source=real(0,CRT) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"var",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Stop/Disable the CGNS queue if it is not to be used
  !
  call disable_cgns_queue(pname,cgns_grid_file)
  !
  ! Collect the coordinates for all the output points that this processor
  ! is responsible for writing to the CGNS files
  !
  call collect_coordinates_for_output(vartmp)
  !
#ifdef DEBUG_POST_PROCESSOR
  if (ncpu == 1) then
    call write_tecplot_file(vartmp,ipelem)
  end if
#endif
  !
  ! Write the coordinates for the output points to the CGNS grid file
  !
  write (700+mypnum,20) my_pts_beg,my_pts_end,my_pts_end-my_pts_beg+1
  do n = 1,size(vartmp,dim=2)
    write (700+mypnum,21) n+my_pts_beg-1,(vartmp(i,n), i=1,size(vartmp,dim=1))
  end do
  20 format (/,2x,"my_pts_beg = ",i0,/, &
               2x,"my_pts_end = ",i0,/, &
               2x,"my_pts_tot = ",i0,/)
  21 format (i10," : ",3es15.6)
  !
  do n = 1,nr
    !
    ! Copy the current coordinate variable into the array var
    ! and convert to the real type required for the CGNS file.
    !
    if (channel_nondim_output_coord) then
      vartmp(n,:) = vartmp(n,:) / 0.003048_wp
    end if
    !
    var(:) = real( vartmp(n,:)/grid_scaling_factor , kind=kind(var) )
    !
    ! Reset the queue.
    !
    call set_cgns_queue(BEGIN_CGNS_QUEUE,pname,cgns_grid_file)
    !
    ! Add the data for current coordinate variable to the output queue.
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_coord_write_data_f")
    call cgp_coord_write_data_f(ifile_n,ibase_n,izone_n,ixyz_n(n), &
                                my_pts_beg,my_pts_end,var,cgierr)
    call cgns_error(pname,cgns_grid_file,"cgp_coord_write_data_f", &
                    cgierr,__LINE__,__FILE__,"Write data: Coordinate",n)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_coord_write_data_f")
    !
    ! Finally, flush the queue containing the
    ! current coordinate variable to the CGNS file.
    !
    call flush_cgns_queue(pname,cgns_grid_file)
    !
  end do
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
  !
  ! Create the data node for the grid connectivity
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cgp_section_write_f")
  if (nr == 1) then
    call cgp_section_write_f(ifile_n,ibase_n,izone_n,"Edges", &
                             BAR_2,1_CST,n_total_output_cells, &
                             0_CBT,isect_n,cgierr)
  else if (nr == 2) then
    call cgp_section_write_f(ifile_n,ibase_n,izone_n,"Quadrilaterals", &
                             QUAD_4,1_CST,n_total_output_cells, &
                             0_CBT,isect_n,cgierr)
  else if (nr == 3) then
    call cgp_section_write_f(ifile_n,ibase_n,izone_n,"Hexahedrals", &
                             HEXA_8,1_CST,n_total_output_cells, &
                             0_CBT,isect_n,cgierr)
  end if
  call cgns_error(pname,cgns_grid_file,"cgp_section_write_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cgp_section_write_f")
  !
  ! Write the element connectivity in parallel
  !
  ! Reset the queue.
  !
  call set_cgns_queue(BEGIN_CGNS_QUEUE,pname,cgns_grid_file)
  !
  ! Add the element connectivity data to the output queue.
  !
  write (600+mypnum,10) my_cells_beg,my_cells_end,my_cells_end-my_cells_beg+1
  do n = 1,size(ipelem,dim=2)
    write (600+mypnum,11) n+my_cells_beg-1,(ipelem(i,n), i=1,size(ipelem,dim=1))
  end do
  10 format (/,2x,"my_cells_beg = ",i0,/, &
               2x,"my_cells_end = ",i0,/, &
               2x,"my_cells_tot = ",i0,/)
  11 format (i10," : ",8i9)
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cgp_elements_write_data_f")
  call cgp_elements_write_data_f(ifile_n,ibase_n,izone_n,isect_n, &
                                 my_cells_beg,my_cells_end,ipelem,cgierr)
  call cgns_error(pname,cgns_grid_file,"cgp_elements_write_data_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cgp_elements_write_data_f")
  !
  ! Finally, flush the queue containing the
  ! element connectivity to the CGNS file.
  !
  call flush_cgns_queue(pname,cgns_grid_file)
  !
  ! Deallocate ipelem
  !
  if (allocated(ipelem)) then
    deallocate ( ipelem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ipelem",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Close the CGNS FILE
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cgp_close_f")
  call cgp_close_f(ifile_n,cgierr)
  call cgns_error(pname,cgns_grid_file,"cgp_close_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cgp_close_f")
  !
  call cgns_timer(stop_cgns_file_timer,profile_io_cgns_subs,cgns_grid_file)
  !
  if (profile_io_cgns) then
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum,"before mpi_barrier",__LINE__,__FILE__
#endif
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum," after mpi_barrier",__LINE__,__FILE__
#endif
    if (mypnum == glb_root) then
      end_wall_time = mpi_wtime()
      cgns_wall_time = end_wall_time - beg_wall_time
      write (iout,2) pname,trim(adjustl(cgns_grid_file)),cgns_wall_time
      flush (iout)
    end if
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum,"before mpi_barrier",__LINE__,__FILE__
#endif
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum," after mpi_barrier",__LINE__,__FILE__
#endif
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
end subroutine write_grid_parallel_cgns
!
!###############################################################################
!
subroutine write_solution_parallel_cgns(sol_file,tname)
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : profile_io_cgns,profile_io_cgns_subs
  use ovar,   only : max_wall_time => max_cgns_wall_time
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: sol_file
  character(len=*), intent(in) :: tname
  !
  !.. Local Scalars ..
  integer :: l,n,npts,npvar,ierr
  !
  real(wp) :: beg_wall_time
  real(wp) :: end_wall_time
  real(wp) :: cgns_wall_time
  !
  character(len=CGLEN) :: sect_name
  character(len=170)   :: cgns_sol_file
  character(len=300)   :: link_path
  !
  !.. Local Allocatable Arrays ..
#ifdef SPECIAL_FOR_INTEL
  character(len=:), allocatable :: char_sname(:)
#else
  character(len=CGLEN), allocatable :: char_sname(:)
#endif
  !
  !.. Local Arrays ..
  integer(CBT) :: ifile_s,ibase_s,izone_s,isol_s
  !
  !.. Local Allocatable Arrays ..
  integer(CBT), allocatable :: isvar_s(:)
  !
  !.. Local Allocatable Arrays ..
  real(CRT), allocatable :: var(:)
  real(wp),  allocatable :: vartmp(:)
  real(wp),  allocatable :: output_usp(:,:)
  real(wp),  allocatable :: output_gpsp(:,:)
  real(wp),  allocatable :: output_dusp(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_solution_parallel_cgns"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
#ifndef SPECIAL_FOR_PGI
  !
  if (profile_io_cgns) then
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum,"before mpi_barrier",__LINE__,__FILE__
#endif
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum," after mpi_barrier",__LINE__,__FILE__
#endif
    if (mypnum == glb_root) then
      write (iout,1) pname
      flush (iout)
      beg_wall_time = mpi_wtime()
    end if
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum,"before mpi_barrier",__LINE__,__FILE__
#endif
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum," after mpi_barrier",__LINE__,__FILE__
#endif
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
  ! Open the cgns output file if iflag is equal to open_CGNS
  !
  call cgns_timer(start_cgns_file_timer,profile_io_cgns_subs,cgns_sol_file)
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cgp_open_f")
  call cgp_open_f(cgns_sol_file,CG_MODE_WRITE,ifile_s,cgierr)
  call io_error(pname,cgns_sol_file,1,__LINE__,__FILE__,cgierr)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cgp_open_f")
  !
  ! Set the CGNS base (who knows what this is for)
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_base_write_f")
  call cg_base_write_f(ifile_s,base_name,int(nr,kind=CBT), &
                       int(nr,kind=CBT),ibase_s,cgierr)
  call cgns_error(pname,cgns_sol_file,"cg_base_write_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_base_write_f")
  !
  ! Write the title for this solution
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_zone_write_f")
  call cg_zone_write_f(ifile_s,ibase_s,zone_name,zone_size, &
                       UNSTRUCTURED,izone_s,cgierr)
  call cgns_error(pname,cgns_sol_file,"cg_zone_write_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_zone_write_f")
  !
  ! Create a link to the grid coordinates in the grid file
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cg_goto_f")
  call cg_goto_f(ifile_s,ibase_s,cgierr, &
                 "Zone_t",izone_s, &
                 CGNS_GOTO_END)
  call cgns_error(pname,cgns_sol_file,"cg_goto_f",cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cg_goto_f")
  !
  link_path = base_name//"/"//zone_name//"/GridCoordinates"
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  call cg_link_write_f("GridCoordinates",trim(cgns_grid_fname), &
                       trim(link_path),cgierr)
  call cgns_error(pname,cgns_sol_file,"cg_link_write_f(GridCoordinates)", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  !
  ! Create a link to the element connectivity in the grid file
  !
  if (nr == 1) then
    sect_name = "Edges"
  else if (nr == 2) then
    sect_name = "Quadrilaterals"
  else if (nr == 3) then
    sect_name = "Hexahedrals"
  end if
  link_path = base_name//"/"//zone_name//"/"//trim(sect_name)
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  call cg_link_write_f(trim(sect_name),trim(cgns_grid_fname), &
                       trim(link_path),cgierr)
  call cgns_error(pname,cgns_sol_file, &
                  "cg_link_write_f("//trim(sect_name)//")", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  !
  ! Write the flow variables
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_sol_write_f")
  call cg_sol_write_f(ifile_s,ibase_s,izone_s,tname,VERTEX,isol_s,cgierr)
  call cgns_error(pname,cgns_sol_file,"cg_sol_write_f",cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_sol_write_f")
  !
  ! Count the total number of variables that will be written to the CGNS file
  !
  call get_num_cgns_variables(char_sname,npvar)
  !
  ! Allocate isvar_s
  !
  allocate ( isvar_s(1:npvar) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"isvar_s",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate var and vartmp if they arent already allocated
  !
  npts = n_my_output_pts
  !
  allocate ( vartmp(1:npts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"vartmp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( var(1:npts) , source=real(0,CRT) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"var",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( output_usp(1:npts,1:nq) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"output_usp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( output_gpsp(1:npts,1:ng) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"output_gpsp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( output_dusp(1:npts,1:nr,1:nq) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"output_dusp",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Collect the solution into output_usp, output_gpsp, and output_dusp
  !
  call collect_solution_for_output(output_usp,output_gpsp,output_dusp)
  !
  ! ################################################
  ! #####   WRITING SOLUTION FIELD VARIABLES   #####
  ! ################################################
  !
  ! Create data nodes for each of the solution variables
  !
  call write_field_headers(ifile_s,ibase_s,izone_s,isol_s,isvar_s, &
                           cgns_sol_file,char_sname)
  !
  ! Now create the output data and write it to the CGNS file in parallel
  !
  n = 0
  do l = 1,size(char_sname)
    !
    call get_output_variable(ifile_s,ibase_s,izone_s,isol_s,isvar_s, &
                             cgns_sol_file,char_sname(l),n, &
                             output_usp,output_gpsp,output_dusp,vartmp,var)
    !
  end do
  !
  ! ################################################
  ! #####   FINISHED WRITING FIELD VARIABLES   #####
  ! ################################################
  !
  call set_cgns_queue(STOP_CGNS_QUEUE,pname,cgns_sol_file)
  !
  ! Close the CGNS FILE
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cgp_close_f")
  call cgp_close_f(ifile_s,cgierr)
  call io_error(pname,cgns_grid_file,2,__LINE__,__FILE__,cgierr)
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
  if (allocated(output_usp)) then
    deallocate ( output_usp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"output_usp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(output_gpsp)) then
    deallocate ( output_gpsp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"output_gpsp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(output_dusp)) then
    deallocate ( output_dusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"output_dusp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(isvar_s)) then
    deallocate ( isvar_s , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"isvar_s",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (profile_io_cgns) then
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum,"before mpi_barrier",__LINE__,__FILE__
#endif
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum," after mpi_barrier",__LINE__,__FILE__
#endif
    if (mypnum == glb_root) then
      end_wall_time = mpi_wtime()
      cgns_wall_time = end_wall_time - beg_wall_time
      max_wall_time = max(max_wall_time,cgns_wall_time)
      write (iout,2) pname,trim(adjustl(cgns_sol_file)),cgns_wall_time
      flush (iout)
    end if
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum,"before mpi_barrier",__LINE__,__FILE__
#endif
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum," after mpi_barrier",__LINE__,__FILE__
#endif
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
end subroutine write_solution_parallel_cgns
!
!###############################################################################
!
!pure &
subroutine get_num_cgns_variables(char_sname,npvar)
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : governing_equations
  !
  use metrics_mod, only : time_metrics
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: char_sname(:)
  integer,       intent(inout) :: npvar
  !
  !.. Local Scalars ..
  integer :: l
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_num_cgns_variables"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Count the total number of solution variables that will be written
  ! and allocate isvar_s accordingly
  !
  npvar = 0
  do l = 1,size(char_sname)
    if (trim(char_sname(l)) == "Velocity") then
      npvar = npvar + nr
    else if (trim(char_sname(l)) == "Vorticity") then
      npvar = npvar + merge(nr,1,nr==3)
    else if (trim(char_sname(l)) == "Momentum") then
      npvar = npvar + nr
    else if (trim(char_sname(l)) == "GradDensity") then
      npvar = npvar + nr
    else if (trim(char_sname(l)) == "GradMomentum") then
      npvar = npvar + nr*nr
    else if (trim(char_sname(l)) == "GradVelocity") then
      npvar = npvar + nr*nr
    else if (trim(char_sname(l)) == "GradEnergy") then
      npvar = npvar + nr
    else if (trim(char_sname(l)) == "GradTemperature") then
      npvar = npvar + nr
  ! else if (trim(char_sname(l)) == "DmatTemperature") then
  !   npvar = npvar + nr
    else if (trim(char_sname(l)) == "GradMassFractions") then
      npvar = npvar + ncs*nr
    else if (trim(char_sname(l)) == "GradSpeciesDensities") then
      npvar = npvar + ncs*nr
    else if (trim(char_sname(l)) == "Viscosity") then
      if (governing_equations == NavierStokes_Eqns) then
        npvar = npvar + 1
      end if
  ! else if (trim(char_sname(l)) == "ResolutionMarker") then
  !   npvar = npvar + 2
  ! else if (trim(char_sname(l)) == "JumpMarker") then
  !   npvar = npvar + 2
  ! else if (trim(char_sname(l)) == "TimeMetrics") then
  !   if (allocated(time_metrics)) then
  !     npvar = npvar + nr + 1
  !   end if
  ! else if (trim(char_sname(l)) == "CoMetrics") then
  !!  npvar = npvar + nr*nr
  ! else if (trim(char_sname(l)) == "ContraMetrics") then
  !   npvar = npvar + nr*nr
    else if (trim(char_sname(l)) == "MassFractions") then
      if (ncs > 0) then
        npvar = npvar + ncs
      end if
    else if (trim(char_sname(l)) == "MoleFractions") then
      if (ncs > 0) then
        npvar = npvar + ncs
      end if
    else if (trim(char_sname(l)) == "SpeciesDensities") then
      if (ncs > 0) then
        npvar = npvar + ncs
      end if
    else if (trim(char_sname(l)) == "SpeciesEnthalpy") then
      if (ncs > 0) then
        npvar = npvar + ncs
      end if
    else if (trim(char_sname(l)) == "SpeciesCp") then
      if (ncs > 0) then
        npvar = npvar + ncs
      end if
    else if (trim(char_sname(l)) == "cpi/R") then
      if (ncs > 0) then
        npvar = npvar + ncs
      end if
    else if (trim(char_sname(l)) == "hi/RT") then
      if (ncs > 0) then
        npvar = npvar + ncs
      end if
    else if (trim(char_sname(l)) == "hi/R") then
      if (ncs > 0) then
        npvar = npvar + ncs
      end if
    else if (trim(char_sname(l)) == "dTdx_Breakdown") then
      if (ncs > 0) then
        npvar = npvar + (2*nr + 2*ncs + 7 + ncs*nr)
      end if
    else if (trim(char_sname(l)) == "QCriterion") then
      npvar = npvar + 2
    else
      npvar = npvar + 1
    end if
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_num_cgns_variables
!
!###############################################################################
!
subroutine write_field_headers(ifile_s,ibase_s,izone_s,isol_s,isvar_s, &
                               cgns_sol_file,char_sname)
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : governing_equations
  use ovar,   only : profile_io_cgns_subs
  !
  use metrics_mod, only : time_metrics
  !
  use thermo_mod, only : species_name
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: cgns_sol_file
  character(len=*), intent(in) :: char_sname(:)
  integer,          intent(in) :: ifile_s
  integer,          intent(in) :: ibase_s
  integer,          intent(in) :: izone_s
  integer,          intent(in) :: isol_s
  integer,       intent(inout) :: isvar_s(:)
  !
  !.. Local Scalars ..
  integer :: i,k,l,n,m
  !
  character(len=CGLEN) :: varname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_field_headers"
  character(len=*), parameter :: char_xyz(1:3) = ["X","Y","Z"]
  character(len=*), parameter :: char_ijk(1:3) = ["Xi  ","Eta ","Zeta"]
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Create data nodes for the solution variables
  !
  n = 0
  do k = 1,size(char_sname)
    if (trim(char_sname(k)) == "Velocity") then
      do m = 1,nr
        n = n + 1
        varname = "Velocity"//char_xyz(m)
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
      end do
    else if (trim(char_sname(k)) == "Vorticity") then
      do m = 1,nr
        n = n + 1
        varname = "Vorticity"//char_xyz(merge(m,3,nr==3))
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        if (nr < 3) exit
      end do
    else if (trim(char_sname(k)) == "Momentum") then
      do m = 1,nr
        n = n + 1
        varname = "Momentum"//char_xyz(m)
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
      end do
    else if (trim(char_sname(k)) == "GradDensity") then
      do l = 1,nr
        n = n + 1
        varname = "DensityGradient"//char_xyz(l)
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
      end do
    else if (trim(char_sname(k)) == "GradMomentum") then
      do m = 1,nr
        do l = 1,nr
          n = n + 1
          varname = "Momentum"//char_xyz(m)//"Gradient"//char_xyz(l)
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      end do
    else if (trim(char_sname(k)) == "GradVelocity") then
      do m = 1,nr
        do l = 1,nr
          n = n + 1
          varname = "Velocity"//char_xyz(m)//"Gradient"//char_xyz(l)
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      end do
    else if (trim(char_sname(k)) == "GradEnergy") then
      do l = 1,nr
        n = n + 1
        varname = "EnergyGradient"//char_xyz(l)
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
      end do
    else if (trim(char_sname(k)) == "GradTemperature") then
      do l = 1,nr
        n = n + 1
        varname = "TemperatureGradient"//char_xyz(l)
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
      end do
  ! else if (trim(char_sname(k)) == "DmatTemperature") then
  !   do l = 1,nr
  !     n = n + 1
  !     varname = "TemperatureDMatrix"//char_xyz(l)
  !     call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
  !   end do
    else if (trim(char_sname(k)) == "GradMassFractions") then
      do m = 1,ncs
        do l = 1,nr
          n = n + 1
          varname = trim(adjustl(species_name(m))) // &
                    "MassFractionGradient"//char_xyz(l)
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      end do
    else if (trim(char_sname(k)) == "GradSpeciesDensities") then
      do m = 1,ncs
        do l = 1,nr
          n = n + 1
          varname = trim(adjustl(species_name(m))) // &
                    "DensityGradient"//char_xyz(l)
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      end do
    else if (trim(char_sname(k)) == "Viscosity") then
      if (governing_equations == NavierStokes_Eqns) then
        n = n + 1
        varname = trim(char_sname(k))
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
      end if
  ! else if (trim(char_sname(k)) == "ResolutionMarker") then
  !   !
  !   n = n + 1
  !   varname = trim(char_sname(k))
  !   call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
  !   !
  !   n = n + 1
  !   varname = "Scaled"//trim(char_sname(k))
  !   call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
  !   !
  ! else if (trim(char_sname(k)) == "JumpMarker") then
  !   !
  !   n = n + 1
  !   varname = trim(char_sname(k))
  !   call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
  !   !
  !   n = n + 1
  !   varname = "Scaled"//trim(char_sname(k))
  !   call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
  !   !
  ! else if (trim(char_sname(k)) == "CoMetrics") then
  !!  do m = 1,nr
  !!    do i = 1,nr
  !!      n = n + 1
  !!      varname = "Covariant Metric d" // &
  !!                trim(adjustl(char_xyz(m))) // "d" // char_ijk(i)
  !!      call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
  !!    end do
  !!  end do
  ! else if (trim(char_sname(k)) == "ContraMetrics") then
  !   do m = 1,nr
  !     do i = 1,nr
  !       n = n + 1
  !       varname = "Contravariant Metric d" // &
  !                 trim(adjustl(char_ijk(m))) // "d" // char_xyz(i)
  !       call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
  !     end do
  !   end do
  ! else if (trim(char_sname(k)) == "TimeMetrics") then
  !   if (allocated(time_metrics)) then
  !     do m = 1,nr
  !       n = n + 1
  !       varname = "Time Metric d" // char_xyz(m)
  !       call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
  !     end do
  !     n = n + 1
  !     varname = "Time Metric dA"
  !     call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
  !   end if
    else if (trim(char_sname(k)) == "MassFractions") then
      if (ncs > 0) then
        do m = 1,ncs
          n = n + 1
          varname = "MassFraction"//trim(adjustl(species_name(m)))
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      end if
    else if (trim(char_sname(k)) == "MoleFractions") then
      if (ncs > 0) then
        do m = 1,ncs
          n = n + 1
          varname = "MoleFraction"//trim(adjustl(species_name(m)))
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      end if
    else if (trim(char_sname(k)) == "SpeciesDensities") then
      if (ncs > 0) then
        do m = 1,ncs
          n = n + 1
          varname = "SpeciesDensity"//trim(adjustl(species_name(m)))
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      end if
    else if (trim(char_sname(k)) == "SpeciesEnthalpy") then
      if (ncs > 0) then
        do m = 1,ncs
          n = n + 1
          varname = "SpeciesEnthalpy"//trim(adjustl(species_name(m)))
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      end if
    else if (trim(char_sname(k)) == "SpeciesCp") then
      if (ncs > 0) then
        do m = 1,ncs
          n = n + 1
          varname = "SpeciesCp"//trim(adjustl(species_name(m)))
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      end if
    else if (trim(char_sname(k)) == "cpi/R") then
      if (ncs > 0) then
        do m = 1,ncs
          n = n + 1
          varname = "cpi_over_R_"//trim(adjustl(species_name(m)))
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      end if
    else if (trim(char_sname(k)) == "hi/R") then
      if (ncs > 0) then
        do m = 1,ncs
          n = n + 1
          varname = "hi_over_R_"//trim(adjustl(species_name(m)))
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      end if
    else if (trim(char_sname(k)) == "hi/RT") then
      if (ncs > 0) then
        do m = 1,ncs
          n = n + 1
          varname = "hi_over_RT_"//trim(adjustl(species_name(m)))
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      end if
    else if (trim(char_sname(k)) == "dTdx_Breakdown") then
      if (ncs > 0) then
        do l = 1,nr
          n = n + 1
          varname = "CPG_dTd"//char_xyz(l)
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
        do m = 1,ncs
          n = n + 1
          varname = "Cv("//trim(adjustl(species_name(m)))//")"
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
        do m = 1,ncs
          n = n + 1
          varname = "dCp("//trim(adjustl(species_name(m)))//")dT"
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
        !
        n = n + 1
        varname = "dCp(mix)dT"
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        !
        n = n + 1
        varname = "T_from_gradcv2gradpv"
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        !
        n = n + 1
        varname = "Cp(mix)"
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        !
        n = n + 1
        varname = "Gamma(Mix)"
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        !
        n = n + 1
        varname = "dTdx_lead_coefficient"
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        !
        n = n + 1
        varname = "CPG_constant_Cv"
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        !
        n = n + 1
        varname = "CPG_1_over_constant_Cv"
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        !
        do l = 1,nr
          n = n + 1
          varname = "TPG_dTd"//char_xyz(l)
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
        !
        do m = 1,ncs
          do l = 1,nr
            n = n + 1
            varname = "dY("//trim(adjustl(species_name(m)))// &
                      ")d"//char_xyz(l)//"_from_dTd"//char_xyz(l)
            call cgns_write_field_header(n,varname,isvar_s(n), &
                                         __LINE__,__FILE__)
          end do
        end do
      end if
    else if (trim(char_sname(k)) == "QCriterion") then
      !
      n = n + 1
      varname = trim(char_sname(k))
      call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
      !
      n = n + 1
      varname = trim(char_sname(k))//"_Positive"
      call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
      !
    else
      n = n + 1
      if (trim(char_sname(k)) == "Energy") then
        varname = "EnergyStagnationDensity"
      else if (trim(char_sname(k)) == "T0") then
        varname = "TemperatureStagnation"
      else if (trim(char_sname(k)) == "P0") then
        varname = "PressureStagnation"
      else if (trim(char_sname(k)) == "H0") then
        varname = "EnthalpyStagnation"
      else if (trim(char_sname(k)) == "Viscosity") then
        varname = "ViscosityMolecular"
      else if (trim(char_sname(k)) == "Conductivity") then
        varname = "ThermalConductivity"
      else if (trim(char_sname(k)) == "Cp") then
        varname = "SpecificHeatPressure"
      else if (trim(char_sname(k)) == "Gamma") then
        varname = "SpecificHeatRatio"
      else if (trim(char_sname(k)) == "Rgas") then
        varname = "IdealGasConstant"
      else
        varname = trim(char_sname(k))
      end if
      call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
    end if
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  !
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
subroutine collect_coordinates_for_output(var)
  !
  !.. Use Statements ..
  use geovar,  only : nfbnd,ncell,cell
  use geovar,  only : xyz,xyz_nodes
  use flowvar, only : edgexyz_dt,facexyz_dt
  use flowvar, only : face_xyz
  use ovar,    only : mms_opt
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(out) :: var
  !
  !.. Local Scalars ..
  integer :: k,l,n,nc,ne,nf,nn,np1,np2,ierr
  integer :: face_idx,edge_idx,node_idx
  !
  !.. Local Allocatable Arrays ..
  type(face_xyz), allocatable :: temp_facexyz(:)
  integer, allocatable :: edge_pts(:)
  integer, allocatable :: face_pts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "collect_coordinates_for_output"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
 !write (iout,1000) "size(xyz,dim=1) = ",size(xyz,dim=1)
 !write (iout,1000) "size(xyz,dim=2) = ",size(xyz,dim=2)
 !write (iout,1000) "size(var,dim=1) = ",size(var,dim=1)
 !write (iout,1000) "size(var,dim=2) = ",size(var,dim=2)
 !write (iout,1000) "size(output_sol_pts) = ",size(output_sol_pts)
 !write (iout,1000) "minval(output_sol_pts) = ",minval(output_sol_pts)
 !write (iout,1000) "maxval(output_sol_pts) = ",maxval(output_sol_pts)
 !write (iout,1000) "n_my_output_pts = ",n_my_output_pts
 !write (iout,1000) "my_pts_beg = ",my_pts_beg
 !write (iout,1000) "my_pts_end = ",my_pts_end
 !1000 format (2x,a,i0)
  cell_loop: do nc = 1,ncell
    !
    np1 = cell(nc)%beg_sp
    np2 = cell(nc)%end_sp
    !
#ifdef DEBUG_POST_PROCESSOR
    write (1000+mypnum,1) nc
#endif
    do n = np1,np2
      !
      l = output_sol_pts(n) - my_pts_beg + 1
#ifdef DEBUG_POST_PROCESSOR
      write (1000+mypnum,2) n,l,(xyz(nn,n), nn=1,size(var,dim=1))
#endif
      !
      var(:,l) = xyz(:,n)
      !
    end do
    !
  end do cell_loop
  1 format (" Cell ",i0)
  2 format ("    SP ",i0," => ",i0,"  : ",3es14.6)
  3 format (" Face ",i0)
  4 format ("    FP ",i0," => ",i0,"  : ",3es14.6)
  5 format (" Edge ",i0)
  6 format ("    EP ",i0," => ",i0,"  : ",3es14.6)
  7 format (" Node")
  8 format ("    NP ",i0," => ",i0,"  : ",3es14.6)
  !
  !
 !l = 0
 !do nf = 1,size(output_face_pts)
 !  l = l + size(output_face_pts(nf)%out_idx)
 !end do
 !face_pts = init_array(l,0) ! F2003 AUTO-REALLOC
 !l = 0
 !do nf = 1,size(output_face_pts)
 !  do k = 1,size(output_face_pts(nf)%out_idx)
 !    l = l + 1
 !    face_pts(l) = output_face_pts(nf)%out_idx(k)
 !  end do
 !end do
 !write (iout,1000) "size(output_face_pts) = ",size(face_pts)
 !write (iout,1000) "minval(output_face_pts) = ",minval(face_pts)
 !write (iout,1000) "maxval(output_face_pts) = ",maxval(face_pts)
 !deallocate ( face_pts )
  !
  face_loop: do nf = 1,size(output_face_pts)
    !
    face_idx = output_face_pts(nf)%face_idx
    !
#ifdef DEBUG_POST_PROCESSOR
    write (1000+mypnum,3) face_idx
#endif
    do k = 1,size(output_face_pts(nf)%out_idx)
      !
      l = output_face_pts(nf)%out_idx(k) - my_pts_beg + 1
      !
#ifdef DEBUG_POST_PROCESSOR
      write (1000+mypnum,4) k,l,(facexyz_dt(face_idx)%v(n,k), n=1,size(var,dim=1))
#endif
      var(:,l) = facexyz_dt(face_idx)%v(:,k)
      !
    end do
    !
  end do face_loop
  !
  !
 !l = 0
 !do ne = 1,size(output_edge_pts)
 !  l = l + size(output_edge_pts(ne)%out_idx)
 !end do
 !edge_pts = init_array(l,0) ! F2003 AUTO-REALLOC
 !l = 0
 !do ne = 1,size(output_edge_pts)
 !  do k = 1,size(output_edge_pts(ne)%out_idx)
 !    l = l + 1
 !    edge_pts(l) = output_edge_pts(ne)%out_idx(k)
 !  end do
 !end do
 !write (iout,1000) "size(output_edge_pts) = ",size(edge_pts)
 !write (iout,1000) "minval(output_edge_pts) = ",minval(edge_pts)
 !write (iout,1000) "maxval(output_edge_pts) = ",maxval(edge_pts)
 !deallocate ( edge_pts )
  !
  edge_loop: do ne = 1,size(output_edge_pts)
    !
    edge_idx = output_edge_pts(ne)%edge_idx
    !
#ifdef DEBUG_POST_PROCESSOR
    write (1000+mypnum,5) edge_idx
#endif
    do k = 1,size(output_edge_pts(ne)%out_idx)
      !
     !l = output_edge_pts(ne)%out_idx(k)
      l = output_edge_pts(ne)%out_idx(k) - my_pts_beg + 1
      !
#ifdef DEBUG_POST_PROCESSOR
      write (1000+mypnum,6) k,l,(edgexyz_dt(edge_idx)%v(n,k), n=1,size(var,dim=1))
#endif
      var(:,l) = edgexyz_dt(edge_idx)%v(:,k)
      !
    end do
    !
  end do edge_loop
  !
  !
 !write (iout,1000) "size(output_node_pts) = ",size(output_node_pts(:))
 !write (iout,1000) "minval(output_node_pts) = ",minval(output_node_pts(:)%out_idx)
 !write (iout,1000) "maxval(output_node_pts) = ",maxval(output_node_pts(:)%out_idx)
  !
#ifdef DEBUG_POST_PROCESSOR
  write (1000+mypnum,7)
#endif
  node_loop: do nn = 1,size(output_node_pts)
    !
    node_idx = output_node_pts(nn)%node_idx
    !
   !l = output_node_pts(nn)%out_idx
    l = output_node_pts(nn)%out_idx - my_pts_beg + 1
    !
#ifdef DEBUG_POST_PROCESSOR
    write (1000+mypnum,8) node_idx,l,(xyz_nodes(n,node_idx), n=1,size(var,dim=1))
#endif
    var(:,l) = xyz_nodes(:,node_idx)
    !
  end do node_loop
  !
  ! Deallocate edgexyz since it is no longer needed
  !
  deallocate ( edgexyz_dt , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"edgexyz_dt",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Deallocate or shrink facexyz depending on if MMS is enabled
  !
  if (mms_opt /= 0) then
    !
    ! MMS is enabled, so shrink facexyz to only the boundary faces because the
    ! coordinates for the face points on the interior faces are no longer needed
    !
    allocate ( temp_facexyz(1:nfbnd) , source=facexyz_dt(1:nfbnd) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"temp_facexyz",1,__LINE__,__FILE__,ierr, &
                     error_message)
    call move_alloc( from=temp_facexyz , to=facexyz_dt )
    !
  else
    !
    ! MMS is not enabled, so deallocate facexyz_dt
    ! entirely since it is no longer needed
    !
    deallocate ( facexyz_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facexyz_dt",2,__LINE__,__FILE__,ierr,error_message)
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine collect_coordinates_for_output
!
!###############################################################################
!
subroutine collect_solution_for_output(cv,gp,dcvdx)
  !
  !.. Use Statements ..
  use geovar,  only : nr,ncell,cell
  use flowvar, only : usp,gpsp,dusp
  use flowvar, only : faceusp_dt,facegpsp_dt,facedusp_dt
  use flowvar, only : nodeusp_dt,nodegpsp_dt,nodedusp_dt
  use flowvar, only : edgeusp_dt,edgegpsp_dt,edgedusp_dt
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: cv(:,:)
  real(wp), intent(inout) :: gp(:,:)
  real(wp), intent(inout) :: dcvdx(:,:,:)
  !
  !.. Local Scalars ..
  integer :: k,l,m,n,nc,ne,nf,nn,np1,np2
  integer :: face_idx,edge_idx,node_idx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "collect_solution_for_output"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  call get_face_solutions
  !
  call get_edge_and_node_solutions
  !
  call pp_update_gas_properties
  !
  !
  cell_loop: do nc = 1,ncell
    !
    np1 = cell(nc)%beg_sp
    np2 = cell(nc)%end_sp
    !
    do n = np1,np2
      !
      k = output_sol_pts(n)
      !
      do m = 1,nq
        cv(k,m) = usp(n,m)
      end do
      !
      do m = 1,ng
        gp(k,m) = gpsp(n,m)
      end do
      !
      do m = 1,nq
        do l = 1,nr
          dcvdx(k,l,m) = dusp(n,l,m)
        end do
      end do
     !var(:,k) = compute_output_variables(usp(n,:),gpsp(n,:),dusp(n,:,:))
      !
    end do
    !
  end do cell_loop
  !
  !
  !
  face_loop: do nf = 1,size(output_face_pts)
    !
    face_idx = output_face_pts(nf)%face_idx
    !
    do n = 1,size(output_face_pts(nf)%out_idx)
      !
      k = output_face_pts(nf)%out_idx(n)
      !
      do m = 1,nq
        cv(k,m) = faceusp_dt(face_idx)%v(m,n,1)
      end do
      !
      do m = 1,ng
        gp(k,m) = facegpsp_dt(face_idx)%v(m,n,1)
      end do
      !
      do m = 1,nq
        do l = 1,nr
          dcvdx(k,l,m) = facedusp_dt(face_idx)%v(l,m,n,1)
        end do
      end do
     !var(:,k) = compute_output_variables(faceusp_dt(face_idx)%v(:,n,1), &
     !                                   facegpsp_dt(face_idx)%v(:,n,1), &
     !                                   facedusp_dt(face_idx)%v(:,:,n,1))
      !
    end do
    !
  end do face_loop
  !
  !
  !
  edge_loop: do ne = 1,size(output_edge_pts)
    !
    edge_idx = output_edge_pts(ne)%edge_idx
    !
    do n = 1,size(output_edge_pts(ne)%out_idx)
      !
      k = output_edge_pts(ne)%out_idx(n)
      !
      do m = 1,nq
        cv(k,m) = edgeusp_dt(edge_idx)%v(m,n)
      end do
      !
      do m = 1,ng
        gp(k,m) = edgegpsp_dt(edge_idx)%v(m,n)
      end do
      !
      do m = 1,nq
        do l = 1,nr
          dcvdx(k,l,m) = edgedusp_dt(edge_idx)%v(l,m,n)
        end do
      end do
     !var(:,k) = compute_output_variables(edgeusp_dt(edge_idx)%v(:,n), &
     !                                   edgegpsp_dt(edge_idx)%v(:,n), &
     !                                   edgedusp_dt(edge_idx)%v(:,:,n))
      !
    end do
    !
  end do edge_loop
  !
  !
  !
  node_loop: do nn = 1,size(output_node_pts)
    !
    node_idx = output_node_pts(nn)%node_idx
    !
    k = output_node_pts(nn)%out_idx
    !
    do m = 1,nq
      cv(k,m) = nodeusp_dt(node_idx)%v(m)
    end do
    !
    do m = 1,ng
      gp(k,m) = nodegpsp_dt(node_idx)%v(m)
    end do
    !
    do m = 1,nq
      do l = 1,nr
        dcvdx(k,l,m) = nodedusp_dt(node_idx)%v(l,m)
      end do
    end do
   !var(:,k) = compute_output_variables(nodeusp_dt(node_idx)%v(:), &
   !                                   nodegpsp_dt(node_idx)%v(:), &
   !                                   nodedusp_dt(node_idx)%v(:,:))
    !
  end do node_loop
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine collect_solution_for_output
!
!###############################################################################
!
subroutine get_output_variable(ifile_s,ibase_s,izone_s,isol_s,isvar_s, &
                               cgns_sol_file,vname,nv, &
                               cv,gp,dcvdx,vartmp,var)
  !
  !.. Use Statements ..
  use geovar,  only : nr,n_solpts,ncell,cell
  use ovar,    only : aref,rhoref,time_ref,tref
  use ovar,    only : gam,pref,pref_nd,machref
  use ovar,    only : rgasref,cpref,Pr
  use ovar,    only : itestcase,governing_equations
  use ovar,    only : profile_io_cgns_subs
  use ovar,    only : Limiter_Option
  use flowvar, only : dtsp
  !
  use metrics_mod, only : time_metrics,metrics_array,jac
  !
  use module_limiters, only : marker,scaled_marker
  use module_limiters, only : resolution_indicator
  use module_limiters, only : jump_indicator
  use module_limiters, only : deallocate_marker
  use module_limiters, only : deallocate_scaled_marker
  !
  use thermo_mod, only : rspec
  use thermo_mod, only : compute_mixture_rgas
  use thermo_mod, only : compute_species_enthalpy
  use thermo_mod, only : compute_species_cp
  !
  use transport_mod, only : compute_mole_fractions
  !
  use gpu_mod, only : grad_cv_to_grad_pv_gpu
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
  integer  :: i,j,l,m,n,nc,np,n0,npts,ierr
  integer  :: nmi,nmj,l_comp,l_phys,nvar
  integer  :: this_order,limit_order
  real(wp) :: freestream_ref,rhoa,rhoa2,time_ref_inv
  real(wp) :: cp_dim,cond_dim,dim_val
  real(wp) :: rho_inv,vx,vy,vz,a2,q2,s
  real(wp) :: gamm,gm1,gam_over_gm1,mach,pres,temp
  real(wp) :: mi,mj,drhodi,drhodj,dmjdi,dmidj
  real(wp) :: w12,w13,w23,wijwij
  real(wp) :: s12,s13,s23,sijsij,s11,s22,s33,skk
  real(wp) :: dudx,dudy,dudz
  real(wp) :: dvdx,dvdy,dvdz
  real(wp) :: dwdx,dwdy,dwdz
  real(wp) :: qmin,qmax
  !
  integer(CBT) :: ivar_s
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: vartmp_at_sp(:)
  real(wp), allocatable :: spcvar(:,:)
  real(wp), allocatable :: dtdx(:,:)
  real(wp), allocatable :: dpvdx(:,:,:)
  real(wp), allocatable :: other_variables(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_output_variable"
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
    case ("GradVelocity","GradTemperature","GradMassFractions", &
          "GradTempMagnitude","VorticityMagnitude","StrainRateMagnitude", &
          "QCriterion","DivergenceOfVelocity")
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
  ! Compute the current output variable and write it to the CGNS file
  !
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
          vartmp(n) = cv(n,m)
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
        call cgns_write_field_data
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
  ! case ("DmatTemperature") ! Temperature gradient using derivative matrices
  !   !!{{{
  !   !
  !   allocate ( dtdx(1:n_solpts,1:nr) , source=zero , &
  !              stat=ierr , errmsg=error_message )
  !   call alloc_error(pname,"dtdx",1,__LINE__,__FILE__,ierr,error_message)
  !   !
  !   call compute_dmatrix_grad_temperature(dtdx)
  !   !
  !   if (.not. allocated(vartmp_at_sp)) then
  !     allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
  !                stat=ierr , errmsg=error_message )
  !     call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
  !                      ierr,error_message)
  !   end if
  !   !
  !   do l = 1,nr
  !     !
  !     nv = nv + 1
  !     ivar_s = isvar_s(nv)
  !     !
  !     do n = 1,n_solpts
  !       vartmp_at_sp(n) = dtdx(n,l) * tref
  !     end do
  !     !
  !     call interpolate_for_output( vartmp_at_sp , vartmp )
  !     !
  !     call cgns_write_field_data
  !     !
  !   end do
  !   !
  !   if (allocated(dtdx)) then
  !     deallocate ( dtdx , stat=ierr , errmsg=error_message )
  !     call alloc_error(pname,"dtdx",2,__LINE__,__FILE__,ierr,error_message)
  !   end if
  !   !
  !   !!}}}
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
        call cgns_write_field_data
        !
      end do
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
      call cgns_write_field_data
      !!}}}
    case ("P0") ! Total Pressure
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      ! NOTE: vartmp has been pre-loaded with the Mach number, and the auxiliary
      !       variable in gp has already been switched to gamma
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
      call cgns_write_field_data
      !!}}}
    case ("T0") ! Total Temperature
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      ! NOTE: vartmp has been pre-loaded with the Mach number, and the auxiliary
      !       variable in gp has already been switched to gamma
      !
      do n = 1,npts
        mach = vartmp(n)
        temp = gp(n,ntemp) * tref
        gm1 = gp(n,ngamm) - one
        vartmp(n) = temp*(one + half*gm1*mach*mach)
      end do
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
        i = (l+1) - ( (l+1)/4 )*3
        j = (l+2) - ( (l+2)/4 )*3
        !
        nmi = nmb-1+i
        nmj = nmb-1+j
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
      ! Normalize the whole Q-criterion to the interval 0->1
      !
      qmin = minval(vartmp(:))
      qmax = maxval(vartmp(:))
      !
      vartmp(:) = (vartmp(:) - qmin)/(qmax-qmin)
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
  ! case ("ResolutionMarker")
  !   !!{{{
  !   if (.not. allocated(vartmp_at_sp)) then
  !     allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
  !                stat=ierr , errmsg=error_message )
  !     call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
  !                      ierr,error_message)
  !   end if
  !   !
  !   ! Compute the resolution marker if it has yet to be computed
  !   !
  !   if (Limiter_Option >= 0 .or. &
  !       .not. allocated(marker) .or. &
  !       .not. allocated(scaled_marker)) then
  !     call resolution_indicator(apply_projection=fals)
  !   end if
  !   !
  !   ! Actual resolution marker value
  !   !
  !   nv = nv + 1
  !   ivar_s = isvar_s(nv)
  !   !
  !   do nc = 1,ncell
  !     n0 = cell(nc)%beg_sp - 1
  !     np = cell(nc)%end_sp - n0
  !     do n = 1,np
  !       vartmp_at_sp(n0+n) = marker(nc)
  !     end do
  !   end do
  !   !
  !   call interpolate_for_output( vartmp_at_sp , vartmp )
  !   !
  !   call cgns_write_field_data
  !   !
  !   ! Scaled resolution marker value
  !   !
  !   nv = nv + 1
  !   ivar_s = isvar_s(nv)
  !   !
  !   if (Limiter_Option <= -1) then
  !     do nc = 1,ncell
  !       n0 = cell(nc)%beg_sp - 1
  !       np = cell(nc)%end_sp - n0
  !       this_order = cell(nc)%order
  !       limit_order = int( real(this_order,kind=wp) * &
  !                          (one-scaled_marker(nc)) )
  !       limit_order = max(1,min(this_order,limit_order))
  !       do n = 1,np
  !         vartmp_at_sp(n0+n) = limit_order
  !       end do
  !     end do
  !   else
  !     do nc = 1,ncell
  !       n0 = cell(nc)%beg_sp - 1
  !       np = cell(nc)%end_sp - n0
  !       do n = 1,np
  !         vartmp_at_sp(n0+n) = scaled_marker(nc)
  !       end do
  !     end do
  !   end if
  !   !
  !   call interpolate_for_output( vartmp_at_sp , vartmp )
  !   !
  !   call cgns_write_field_data
  !   !!}}}
  ! case ("JumpMarker")
  !   !!{{{
  !   allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
  !              stat=ierr , errmsg=error_message )
  !   call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
  !                    ierr,error_message)
  !   !
  !   ! Compute the pressure jump marker if it has yet to be computed
  !   !
  !   if (Limiter_Option <= 1 .or. &
  !       .not. allocated(marker) .or. &
  !       .not. allocated(scaled_marker)) then
  !     call jump_indicator
  !   end if
  !   !
  !   ! Actual pressure jump marker value
  !   !
  !   nv = nv + 1
  !   ivar_s = isvar_s(nv)
  !   !
  !   do nc = 1,ncell
  !     n0 = cell(nc)%beg_sp - 1
  !     np = cell(nc)%end_sp - n0
  !     do n = 1,np
  !       vartmp_at_sp(n0+n) = marker(nc)
  !     end do
  !   end do
  !   !
  !   call interpolate_for_output( vartmp_at_sp , vartmp )
  !   !
  !   call cgns_write_field_data
  !   !
  !   ! Scaled pressure jump marker value
  !   !
  !   nv = nv + 1
  !   ivar_s = isvar_s(nv)
  !   !
  !   do nc = 1,ncell
  !     n0 = cell(nc)%beg_sp - 1
  !     np = cell(nc)%end_sp - n0
  !     do n = 1,np
  !       vartmp_at_sp(n0+n) = scaled_marker(nc)
  !     end do
  !   end do
  !   !
  !   call interpolate_for_output( vartmp_at_sp , vartmp )
  !   !
  !   call cgns_write_field_data
  !   !!}}}
  ! case ("Partition")
  !   !!{{{
  !   allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
  !              stat=ierr , errmsg=error_message )
  !   call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
  !                    ierr,error_message)
  !   !
  !   nv = nv + 1
  !   ivar_s = isvar_s(nv)
  !   !
  !   do n = 1,n_solpts
  !     vartmp_at_sp(n) = mypnum+1
  !   end do
  !   !
  !   call interpolate_for_output( vartmp_at_sp , vartmp )
  !   !
  !   call cgns_write_field_data
  !   !!}}}
  ! case ("Timestep")
  !   !!{{{
  !   nv = nv + 1
  !   ivar_s = isvar_s(nv)
  !   !
  !   call interpolate_for_output( dtsp(1:n_solpts) , vartmp , nondim=time_ref )
  !   !
  !   call cgns_write_field_data
  !   !!}}}
  ! case ("CoMetrics")
  !  !!!{{{
  !  !allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
  !  !           stat=ierr , errmsg=error_message )
  !  !call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
  !  !                 ierr,error_message)
  !  !!
  !  !do l_phys = 1,nr
  !  !  do l_comp = 1,nr
  !  !    !
  !  !    nv = nv + 1
  !  !    ivar_s = isvar_s(nv)
  !  !    !
  !  !    do n = 1,n_solpts
  !  !      vartmp_at_sp(n) = metrics_array(n,l_phys,l_comp)
  !  !    end do
  !  !    !
  !  !    call interpolate_for_output( vartmp_at_sp , vartmp )
  !  !    !
  !  !    call cgns_write_field_data
  !  !    !
  !  !  end do
  !  !end do
  !  !!!}}}
  ! case ("ContraMetrics")
  !   !!{{{
  !   allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
  !              stat=ierr , errmsg=error_message )
  !   call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
  !                    ierr,error_message)
  !   !
  !   do l_comp = 1,nr
  !     do l_phys = 1,nr
  !       !
  !       nv = nv + 1
  !       ivar_s = isvar_s(nv)
  !       !
  !       do n = 1,n_solpts
  !         vartmp_at_sp(n) = jac(n)*metrics_array(n,l_comp,l_phys)
  !       end do
  !       !
  !       call interpolate_for_output( vartmp_at_sp , vartmp )
  !       !
  !       call cgns_write_field_data
  !       !
  !     end do
  !   end do
  !   !!}}}
  ! case ("TimeMetrics")
  !   !!{{{
  !   allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
  !              stat=ierr , errmsg=error_message )
  !   call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
  !                    ierr,error_message)
  !   !
  !   if (allocated(time_metrics)) then
  !     !
  !     do l = 1,nr
  !       !
  !       nv = nv + 1
  !       ivar_s = isvar_s(nv)
  !       !
  !       do n = 1,n_solpts
  !         vartmp_at_sp(n) = one/time_metrics(l,n)
  !       end do
  !       !
  !       call interpolate_for_output( vartmp_at_sp , vartmp )
  !       !
  !       call cgns_write_field_data
  !       !
  !     end do
  !     !
  !     nv = nv + 1
  !     ivar_s = isvar_s(nv)
  !     !
  !     do n = 1,n_solpts
  !       vartmp_at_sp(n) = time_metrics(nr+1,n) * time_ref
  !     end do
  !     !
  !     call interpolate_for_output( vartmp_at_sp , vartmp )
  !     !
  !     call cgns_write_field_data
  !     !
  !   end if
  !   !!}}}
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
  end select
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
end subroutine get_output_variable
!
!###############################################################################
!
subroutine get_cgns_sol_variables(variable_list)
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : itestcase,governing_equations
  use ovar,   only : output_vorticity
  use ovar,   only : output_qcriterion
  use ovar,   only : output_conserved_variables
  use ovar,   only : output_temperature
  use ovar,   only : output_thermo_variables
  use ovar,   only : output_gradients
  use ovar,   only : output_transport_variables
  use ovar,   only : output_stagnation_temperature
  use ovar,   only : output_stagnation_pressure
  use ovar,   only : output_stagnation_variables
  use ovar,   only : output_marker_variables
  use ovar,   only : Resolution_Marker_Variable
  use ovar,   only : output_species_mass_fractions
  use ovar,   only : output_species_mole_fractions
  use ovar,   only : output_extra_species_variables
  use ovar,   only : output_debug_species_variables
  use ovar,   only : output_grdnt_species_variables
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
  character(len=*), parameter :: var_order(1:41) = ["Density             ", &
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
                                                    "GradTemperature     ", &
                                                  ! "ResolutionMarker    ", &
                                                  ! "JumpMarker          ", &
                                                  ! "Timestep            ", &
                                                  ! "Partition           ", &
                                                  ! "ContraMetrics       ", &
                                                  ! "TimeMetrics         ", &
                                                    "dTdx_Breakdown      ", &
                                                  ! "DmatTemperature     ", &
                                                  ! "CoMetrics           ", &
                                                    "GradTempMagnitude   ", &
                                                    "GradDensityMagnitude", &
                                                    "VorticityMagnitude  ", &
                                                    "StrainRateMagnitude ", &
                                                    "QCriterion          ", &
                                                    "DivergenceOfVelocity"]
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
  if (itestcase == Taylor_Green_Vortex) then
    call add_sol_variable(variable_list,"Vorticity")
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
  ! Output the marker quantities if requested
  !
  if (output_marker_variables) then
  ! call add_sol_variable(variable_list,"ResolutionMarker")
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
  !
  ! Add species mass fractions to the variable list if using a gas mixture.
  ! Also output any other requested variables for gas mixtures.
  !
  if (ncs > 0) then
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
!pure function compute_output_variables(upt,gppt,dupt) result(return_value)
!  !
!  !.. Formal Arguments ..
!  real(wp), intent(in) :: upt(:)
!  real(wp), intent(in) :: gppt(:)
!  real(wp), intent(in) :: dupt(:,:)
!  !
!  !.. Function Result ..
!  real(wp) :: return_value(1:size(plot_var))
!  !
!  !.. Local Scalars ..
!  integer :: i
!  !
!continue
!  !
!  do i = 1,size(plot_var)
!    return_value(i) = compute_output_variable(upt,gppt,dupt,plot_var(i))
!  end do
!  !
!end function compute_output_variables
!
!###############################################################################
!
!pure function compute_output_variable(upt,gppt,dupt,plot_variable) &
!                               result(return_value)
!  !
!  !.. Use Statements ..
!  use ovar,    only : rhoref,aref,tref,time_ref,gam,pref_nd,itestcase
!  !
!  !.. Formal Arguments ..
!  real(wp),          intent(in) :: upt(:)
!  real(wp),          intent(in) :: gppt(:)
!  real(wp),          intent(in) :: dupt(:,:)
!  type(plot_func_t), intent(in) :: plot_variable
!  !
!  !.. Function Result ..
!  real(wp) :: return_value
!  !
!  !.. Local Scalars ..
!  real(wp) :: gamm
!  !
!continue
!  !
!  if (plot_variable == plot_static_density) then
!    !
!    return_value = upt(nec)                           * rhoref
!    !
!  else if (plot_variable == plot_velocity_x) then
!    !
!    return_value = upt(nmx) / upt(nec)                * aref
!    !
!  else if (plot_variable == plot_velocity_y) then
!    !
!    return_value = upt(nmy) / upt(nec)                * aref
!    !
!  else if (plot_variable == plot_velocity_z) then
!    !
!    return_value = upt(nmz) / upt(nec)                * aref
!    !
!  else if (plot_variable == plot_momentum_x) then
!    !
!    return_value = upt(nmx)                           * rhoref * aref
!    !
!  else if (plot_variable == plot_momentum_y) then
!    !
!    return_value = upt(nmy)                           * rhoref * aref
!    !
!  else if (plot_variable == plot_momentum_z) then
!    !
!    return_value = upt(nmz)                           * rhoref * aref
!    !
!  else if (plot_variable == plot_static_pressure) then
!    !
!    return_value = gppt(npres)                        * rhoref * aref * aref
!    !
!  else if (plot_variable == plot_static_temperature) then
!    !
!    return_value = gppt(ntemp)                        * tref
!    !
!  else if (plot_variable == plot_total_energy) then
!    !
!    return_value = upt(nee)                           * rhoref * aref * aref
!    !
!  else if (plot_variable == plot_mach_number) then
!    !
!    if (ncs == 0) then
!      gamm = gam
!    else
!      gamm = gppt(ngamm)
!    end if
!    !
!    return_value = norm2(upt(nmb:nme)) / sqrt( upt(nec)*gamm*gppt(npres) )
!    !
!  else if (plot_variable == plot_entropy) then
!    !
!    if (ncs == 0) then
!      gamm = gam
!    else
!      gamm = gppt(ngamm)
!    end if
!    !
!    if (itestcase == Inviscid_Gaussian_Bump) then
!      return_value = (gppt(npres)/pref_nd) / (upt(nec)**gamm)
!    else
!      return_value = log(gppt(npres)) - gamm*log(upt(nec))
!    end if
!    !
!  else if (plot_variable == plot_vorticity_x) then
!    !
!    return_value = vorticity_dusp( upt , gppt , dupt , 1 )   / time_ref
!    !
!  else if (plot_variable == plot_vorticity_y) then
!    !
!    return_value = vorticity_dusp( upt , gppt , dupt , 2 )   / time_ref
!    !
!  else if (plot_variable == plot_vorticity_z) then
!    !
!    return_value = vorticity_dusp( upt , gppt , dupt , 3 )   / time_ref
!    !
!  else if (plot_variable == plot_grad_density(1)) then
!    !
!    return_value = dupt(1,nec)                        * rhoref
!    !
!  else if (plot_variable == plot_grad_density(2)) then
!    !
!    return_value = dupt(2,nec)                        * rhoref
!    !
!  else if (plot_variable == plot_grad_density(3)) then
!    !
!    return_value = dupt(3,nec)                        * rhoref
!    !
!  else if (plot_variable == plot_grad_momentum_x(1)) then
!    !
!    return_value = dupt(1,nmx)                        * rhoref * aref
!    !
!  else if (plot_variable == plot_grad_momentum_x(2)) then
!    !
!    return_value = dupt(2,nmx)                        * rhoref * aref
!    !
!  else if (plot_variable == plot_grad_momentum_x(3)) then
!    !
!    return_value = dupt(3,nmx)                        * rhoref * aref
!    !
!  else if (plot_variable == plot_grad_momentum_y(1)) then
!    !
!    return_value = dupt(1,nmy)                        * rhoref * aref
!    !
!  else if (plot_variable == plot_grad_momentum_y(2)) then
!    !
!    return_value = dupt(2,nmy)                        * rhoref * aref
!    !
!  else if (plot_variable == plot_grad_momentum_y(3)) then
!    !
!    return_value = dupt(3,nmy)                        * rhoref * aref
!    !
!  else if (plot_variable == plot_grad_momentum_z(1)) then
!    !
!    return_value = dupt(1,nmz)                        * rhoref * aref
!    !
!  else if (plot_variable == plot_grad_momentum_z(2)) then
!    !
!    return_value = dupt(2,nmz)                        * rhoref * aref
!    !
!  else if (plot_variable == plot_grad_momentum_z(3)) then
!    !
!    return_value = dupt(3,nmz)                        * rhoref * aref
!    !
!  else if (plot_variable == plot_grad_energy(1)) then
!    !
!    return_value = dupt(1,nee)                        * rhoref * aref * aref
!    !
!  else if (plot_variable == plot_grad_energy(2)) then
!    !
!    return_value = dupt(2,nee)                        * rhoref * aref * aref
!    !
!  else if (plot_variable == plot_grad_energy(3)) then
!    !
!    return_value = dupt(3,nee)                        * rhoref * aref * aref
!    !
!  end if
!  !
!end function compute_output_variable
!
!###############################################################################
!
subroutine pp_update_gas_properties()
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  use geovar,    only : face,edge
  use geovar,    only : n_solpts,n_flxpts
  use flowvar,   only : usp,gpsp
  use flowvar,   only : faceusp_dt,facegpsp_dt
  use flowvar,   only : edgeusp_dt,edgegpsp_dt
  use flowvar,   only : nodeusp_dt,nodegpsp_dt
  !
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Local Scalars ..
  integer  :: nf,ne,nn,n0,np,k,m,ierr
  integer  :: n_edgpts,n_nodpts
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: utmp(:,:)
  real(wp), allocatable :: gptmp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "pp_update_gas_properties"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Count the total number of edge points
  !
  n_edgpts = 0
  do ne = 1,size(edgeusp_dt)
    n_edgpts = n_edgpts + size(edgeusp_dt(ne)%v,dim=2)
  end do
  !
  n_nodpts = size(nodeusp_dt)
  !
  ! 1. Update the gas properties for the interior solution points
  !
  call update_thermo_properties(1,n_solpts,usp,gpsp, &
                                loc="cell",caller_name=pname)
  call update_transport_properties(1,n_solpts,usp,gpsp)
  !
  ! 2. Update the gas properties for the face points
  !
  ! Allocate utmp and gptmp for all the face points
  !
  if (allocated(utmp)) then
    deallocate ( utmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"utmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(gptmp)) then
    deallocate ( gptmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"gptmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( utmp(1:n_flxpts,1:nq) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"utmp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( gptmp(1:n_flxpts,1:ng) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"gptmp",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Copy faceusp_dt and facegpsp_dt into utmp and gptmp
  !
  do nf = 1,size(faceusp_dt)
    n0 = face(nf)%beg_fp - 1
    np = face(nf)%end_fp - n0
    do k = 1,np
      do m = 1,nq
        utmp(n0+k,m) = faceusp_dt(nf)%v(m,k,1)
      end do
      do m = 1,ng
        gptmp(n0+k,m) = facegpsp_dt(nf)%v(m,k,1)
      end do
    end do
  end do
  !
  ! Use utmp to update the gas properties in gptmp
  !
  call update_thermo_properties(1,n_flxpts,utmp,gptmp, &
                                loc="face",caller_name=pname)
  call update_transport_properties(1,n_flxpts,utmp,gptmp)
  !
  ! Copy the gas properties in gptmp back into facegpsp_dt
  !
  do nf = 1,size(faceusp_dt)
    n0 = face(nf)%beg_fp - 1
    np = face(nf)%end_fp - n0
    do m = 1,ng
      do k = 1,np
        facegpsp_dt(nf)%v(m,k,1) = gptmp(n0+k,m)
      end do
    end do
  end do
  !
  ! 3. Update the gas properties for the edge points
  !
  ! Allocate utmp and gptmp for all the edge points
  !
  if (allocated(utmp)) then
    deallocate ( utmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"utmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(gptmp)) then
    deallocate ( gptmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"gptmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( utmp(1:n_edgpts,1:nq) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"utmp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( gptmp(1:n_edgpts,1:ng) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"gptmp",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Copy edgeusp_dt and edgegpsp_dt into utmp and gptmp
  !
  n0 = 0
  np = 0
  do ne = 1,size(edgeusp_dt)
    n0 = n0 + np
    np = geom_solpts(Geom_Edge,edge(ne)%order)
    do k = 1,np
      do m = 1,nq
        utmp(n0+k,m) = edgeusp_dt(ne)%v(m,k)
      end do
      do m = 1,ng
        gptmp(n0+k,m) = edgegpsp_dt(ne)%v(m,k)
      end do
    end do
  end do
  !
  ! Use utmp to update the gas properties in gptmp
  !
  call update_thermo_properties(1,n_edgpts,utmp,gptmp,caller_name=pname)
  call update_transport_properties(1,n_edgpts,utmp,gptmp)
  !
  ! Copy the gas properties in gptmp back into edgegpsp_dt
  !
  n0 = 0
  np = 0
  do ne = 1,size(edgeusp_dt)
    n0 = n0 + np
    np = geom_solpts(Geom_Edge,edge(ne)%order)
    do m = 1,ng
      do k = 1,np
        edgegpsp_dt(ne)%v(m,k) = gptmp(n0+k,m)
      end do
    end do
  end do
  !
  ! 4. Update the gas properties for the node points
  !
  ! Allocate utmp and gptmp for all the node points
  !
  if (allocated(utmp)) then
    deallocate ( utmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"utmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(gptmp)) then
    deallocate ( gptmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"gptmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( utmp(1:n_nodpts,1:nq) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"utmp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( gptmp(1:n_nodpts,1:ng) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"gptmp",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Copy nodeusp_dt and nodegpsp_dt into utmp and gptmp
  !
  do nn = 1,size(nodeusp_dt)
    do m = 1,nq
      utmp(nn,m) = nodeusp_dt(nn)%v(m)
    end do
    do m = 1,ng
      gptmp(nn,m) = nodegpsp_dt(nn)%v(m)
    end do
  end do
  !
  ! Use utmp to update the gas properties in gptmp
  !
  call update_thermo_properties(1,n_nodpts,utmp,gptmp,caller_name=pname)
  call update_transport_properties(1,n_nodpts,utmp,gptmp)
  !
  ! Copy the gas properties in gptmp back into nodegpsp_dt
  !
  do nn = 1,size(nodeusp_dt)
    do m = 1,ng
      nodegpsp_dt(nn)%v(m) = gptmp(nn,m)
    end do
  end do
  !
  ! Deallocate utmp and gptmp before leaving
  !
  if (allocated(utmp)) then
    deallocate ( utmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"utmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(gptmp)) then
    deallocate ( gptmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"gptmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine pp_update_gas_properties
!
!###############################################################################
!
subroutine get_face_solutions
  !
  !.. Use Statements ..
  use geovar,    only : n_bnd_flxpts,n_flxpts
  use flowvar,   only : faceusp_array,facegpsp
  !
  use module_bc, only : get_faceusp_bc
  use module_bc, only : get_facedusp_bc
  !
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  use parallel_mod, only : exchange_faceusp
  use parallel_mod, only : wait_faceusp
  use parallel_mod, only : exchange_facedusp
  use parallel_mod, only : wait_facedusp
  !
  !.. Local Scalars ..
  integer :: n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_face_solutions"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  call pp_interp_cv_to_flxpts
  !
  ! Update the thermodynamic and transport gas properties
  ! at all the flux points of the interior faces
  !
  n1 = n_bnd_flxpts+1
  n2 = n_flxpts
  !
  call update_thermo_properties(n1,n2,faceusp_array,facegpsp,caller_name=pname)
  call update_transport_properties(n1,n2,faceusp_array,facegpsp)
  !
  ! Exchange the left and right states for faces on communication boundaries
  !
  if (ncpu > 1) then
    call exchange_faceusp
    call wait_faceusp
  end if
  !
  call get_faceusp_bc
  !
  ! Update the thermodynamic and transport gas properties for
  ! all the exterior flux points of all the boundary faces
  ! NOTE: This was done for the interior flux points
  !       of all boundary faces within get_faceusp_bc
  !
  n1 = 1
  n2 = n_bnd_flxpts
  !
  call update_thermo_properties(n1,n2,faceusp_array(:,:,2),facegpsp(:,:,2), &
                                      loc="face",caller_name=pname)
  call update_transport_properties(n1,n2,faceusp_array(:,:,2),facegpsp(:,:,2))
  !
  call pp_interp_dusp_to_flxpts
  !
  ! Exchange the left and right gradients for faces on communication boundaries
  !
  if (ncpu > 1) then
    call exchange_facedusp
    call wait_facedusp
  end if
  !
  call get_facedusp_bc
  !
  call copy_and_average_face_arrays_to_dt
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_face_solutions
!
!###############################################################################
!
subroutine copy_and_average_face_arrays_to_dt()
  !
  !.. Use Statements ..
  use geovar,    only : nfbnd,bface,nface,face,n_flxpts
  use flowvar,   only : faceusp_dt,faceusp_array
  use flowvar,   only : facegpsp_dt,facegpsp
  use flowvar,   only : facedusp_dt,facedusp_array
  !
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Local Scalars ..
  integer  :: n1,n2,nf,f0,nfp,k,l,m
  real(wp) :: val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "copy_and_average_face_arrays_to_dt"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Zero out the component arrays in the face derived types
  !
#ifdef SPECIAL_FOR_PGI
  call faceusp_dt%zero_face_var
  call facegpsp_dt%zero_face_var
  call facedusp_dt%zero_face_vec
#else
  call faceusp_dt%zero
  call facegpsp_dt%zero
  call facedusp_dt%zero
#endif
  !
  ! Compute a common state for all boundary flux points based
  ! on the type of boundary condition for the boundary face
  !
  do nf = 1,nfbnd
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    if (any(bface(1,nf) == bc_comm)) then
      !
      ! Average left and right states and gradients for
      ! exterior communication and periodic boundary faces
      !
      do m = 1,size(faceusp_array,dim=2)
        do k = 1,nfp
          !
          val = half*(faceusp_array(f0+k,m,1) + faceusp_array(f0+k,m,2))
          faceusp_array(f0+k,m,1) = val
          faceusp_array(f0+k,m,2) = val
          faceusp_dt(nf)%v(m,k,1) = val
          faceusp_dt(nf)%v(m,k,2) = val
          !
          do l = 1,size(facedusp_array,dim=2)
            val = half*(facedusp_array(f0+k,l,m,1) + facedusp_array(f0+k,l,m,2))
            facedusp_array(f0+k,l,m,1) = val
            facedusp_array(f0+k,l,m,2) = val
            facedusp_dt(nf)%v(l,m,k,1) = val
            facedusp_dt(nf)%v(l,m,k,2) = val
          end do
          !
        end do
      end do
      !
    else
      !
      ! Use the right (exterior) state for exterior
      ! boundary faces with physical boundary conditions
      !
      do m = 1,size(faceusp_array,dim=2)
        do k = 1,nfp
          !
          val = faceusp_array(f0+k,m,2)
          faceusp_array(f0+k,m,1) = val
          faceusp_array(f0+k,m,2) = val
          faceusp_dt(nf)%v(m,k,1) = val
          faceusp_dt(nf)%v(m,k,2) = val
          !
          do l = 1,size(facedusp_array,dim=2)
            val = facedusp_array(f0+k,l,m,2)
            facedusp_array(f0+k,l,m,1) = val
            facedusp_array(f0+k,l,m,2) = val
            facedusp_dt(nf)%v(l,m,k,1) = val
            facedusp_dt(nf)%v(l,m,k,2) = val
          end do
          !
        end do
      end do
      !
    end if
    !
  end do
  !
  ! Average left and right states and gradients for interior faces
  !
  do nf = nfbnd+1,nface
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    do m = 1,size(faceusp_array,dim=2)
      do k = 1,nfp
        !
        val = half*(faceusp_array(f0+k,m,1) + faceusp_array(f0+k,m,2))
        faceusp_array(f0+k,m,1) = val
        faceusp_array(f0+k,m,2) = val
        faceusp_dt(nf)%v(m,k,1) = val
        faceusp_dt(nf)%v(m,k,2) = val
        !
        do l = 1,size(facedusp_array,dim=2)
          val = half*(facedusp_array(f0+k,l,m,1) + facedusp_array(f0+k,l,m,2))
          facedusp_array(f0+k,l,m,1) = val
          facedusp_array(f0+k,l,m,2) = val
          facedusp_dt(nf)%v(l,m,k,1) = val
          facedusp_dt(nf)%v(l,m,k,2) = val
        end do
        !
      end do
    end do
    !
  end do
  !
  ! Update the gas properties using the newly updated faceusp_array
  !
  n1 = 1
  n2 = n_flxpts
  !
  call update_thermo_properties(n1,n2,faceusp_array(:,:,1),facegpsp(:,:,1), &
                                      loc="face",caller_name=pname)
  call update_transport_properties(n1,n2,faceusp_array(:,:,1),facegpsp(:,:,1))
  !
  ! Copy the updated gas properties to facegpsp_dt
  !
  do nf = 1,nface
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    do m = 1,size(facegpsp,dim=2)
      do k = 1,nfp
        val = facegpsp(f0+k,m,1)
        facegpsp(f0+k,m,1) = val
        facegpsp(f0+k,m,2) = val
        facegpsp_dt(nf)%v(m,k,1) = val
        facegpsp_dt(nf)%v(m,k,2) = val
      end do
    end do
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine copy_and_average_face_arrays_to_dt
!
!###############################################################################
!
subroutine pp_interp_cv_to_flxpts()
  !
  use order_mod,         only : maxSP,maxFP
  use geovar,            only : ncell,cell,face,fp
  use interpolation_mod, only : interp
  use flowvar,           only : usp
  use flowvar,           only : faceusp_array
  !
  !.. Local Scalars ..
  integer  :: n,nc,n0,n1,n2,k,m,f0
  integer  :: nf,kf,lf,sf,i,j,nfp
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  integer  :: fp_offset,face_rot
  integer  :: num_cell_faces
  real(wp) :: a,val
  !
  !.. Local Arrays ..
  integer, dimension(1:maxFP) :: fpidx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "pp_interp_cv_to_flxpts_nd"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1 ! beginning index of sol-pts for this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Now interpolate the interior solution gradients to the flux points
    !
    ! Loop over the faces of this cell
    !
    num_cell_faces = size(cell(nc)%face)
    !
    do n = 1,num_cell_faces
      !
      nf = cell(nc)%face(n)%idx  ! index of the current face
      sf = cell(nc)%face(n)%side ! side of face that cell nc is on
      !
      face_geom = cell(nc)%face(n)%geom   ! geom  of current face
      face_order = cell(nc)%face(n)%order ! order of current face
      !
      fp_offset = cell(nc)%face(n)%fp_offset
      face_rot = cell(nc)%face(n)%rotation
      !
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      fpidx(1:nfp) = fp_offset + fp%geord(face_geom,face_order)% &
                                    rot(face_rot)%pts(1:nfp)
      !
      do kf = 1,nfp
        do m = 1,nq
          !
          ! kf: index of the flx-pt local to this face
          ! lf: index of flx-pt within interp%mat
          !
          lf = fpidx(kf)
          val = zero
          !
          do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1,&
                 interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
            !
            j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
            a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
            !
            val = val + a * usp(n0+j,m)
            !
          end do
          !
          faceusp_array(f0+kf,m,sf) = val
          !
        end do
        !
      end do
      !
    end do
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine pp_interp_cv_to_flxpts
!
!###############################################################################
!
subroutine pp_interp_dusp_to_flxpts()
  !
  !.. Use Statements ..
  use order_mod,         only : maxFP
  use geovar,            only : nr,ncell,cell,face,fp
  use interpolation_mod, only : interp
  use flowvar,           only : dusp,facedusp_array
  !
  !.. Local Scalars ..
  integer  :: n,nc,n0,k,l,m,f0
  integer  :: nf,kf,lf,sf,i,j,nfp
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  integer  :: fp_offset,face_rot
  real(wp) :: a,val
  !
  !.. Local Arrays ..
  integer, dimension(1:maxFP) :: fpidx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "pp_interp_dusp_to_flxpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1 ! beginning index of sol-pts for this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Now interpolate the interior solution gradients to the flux points
    !
    ! Loop over the faces of this cell
    !
    do n = 1,size(cell(nc)%face)
      !
      nf = cell(nc)%face(n)%idx  ! index of the current face
      sf = cell(nc)%face(n)%side ! side of face that cell nc is on
      !
      face_geom = cell(nc)%face(n)%geom   ! geom  of current face
      face_order = cell(nc)%face(n)%order ! order of current face
      !
      fp_offset = cell(nc)%face(n)%fp_offset
      face_rot = cell(nc)%face(n)%rotation
      !
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      fpidx(1:nfp) = fp_offset + fp%geord(face_geom,face_order)% &
                                    rot(face_rot)%pts(1:nfp)
      !
      do kf = 1,nfp
        !
        ! kf: index of the flx-pt local to this face
        ! lf: index of flx-pt within interp%mat
        !
        do m = 1,nq
          do l = 1,nr
            !
            lf = fpidx(kf)
            val = zero
            !
            do i = interp(this_geom,this_order)% &
                             toFace(face_order)%row_ptr(lf)+1, &
                   interp(this_geom,this_order)% &
                             toFace(face_order)%row_ptr(lf+1)
              !
              j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
              a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
              !
              val = val + a * dusp(n0+j,l,m)
              !
            end do
            !
            facedusp_array(f0+kf,l,m,sf) = val
            !
          end do
        end do
        !
      end do
      !
    end do
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine pp_interp_dusp_to_flxpts
!
!###############################################################################
!
subroutine get_edge_and_node_solutions
  !
  !.. Use Statements ..
  use order_mod,  only : maxSP
  use geovar,     only : nr,ncell,cell
  use flowvar,    only : usp,dusp
  use flowvar,    only : edgeusp_dt,edgedusp_dt
  use flowvar,    only : nodeusp_dt,nodedusp_dt
  use interpolation_mod, only : interp
  !
  use parallel_mod, only : exch_edge_and_node_solutions
  !
  !.. Local Scalars ..
  integer :: i,j,nc,n0,k,l,m,n,ne,nn,kf,lf
  integer :: this_geom,this_order,edge_order
  real(wp) :: a,ave
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq)      :: upt
  real(wp), dimension(1:nr,1:nq) :: dupt
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_edge_and_node_solutions"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Zero out the nodeusp, nodedusp, edgeusp, and edgedusp arrays
  !
#ifdef SPECIAL_FOR_PGI
  call nodeusp_dt%zero_ne1d_t
  call nodedusp_dt%zero_ne2d_t
  call edgeusp_dt%zero_ne2d_t
  call edgedusp_dt%zero_ne3d_t
#else
  call nodeusp_dt%zero
  call nodedusp_dt%zero
  call edgeusp_dt%zero
  call edgedusp_dt%zero
#endif
  !
  ! Loop over all the cells on the current processor and compute the
  ! contribution to the edge and node solutions from each cell
  !
  cell_loop: do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1 ! beginning index of sol-pts for this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Interpolate the solution and gradients to the edge points
    ! of the cell and add to the edgeusp and edgedusp arrays
    !
    cell_edges: do n = 1,size(cell(nc)%edge)
      !
      ne = cell(nc)%edge(n)%idx
      ave = cell(nc)%edge(n)%ave_coef
      edge_order = cell(nc)%edge(n)%order
      !
      edge_pts: do kf = 1,edge_order+1
        !
        lf = cell(nc)%edge(n)%ep_offset + merge( edge_order+2-kf , kf , &
                                                cell(nc)%edge(n)%rev_dir)
        !
        upt(1:nq) = zero
        dupt(1:nr,1:nq) = zero
        !
        do i = interp(this_geom,this_order)%toEdge(edge_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)%toEdge(edge_order)%row_ptr(lf+1)
          !
          j = interp(this_geom,this_order)%toEdge(edge_order)%row_idx(i)
          a = interp(this_geom,this_order)%toEdge(edge_order)%row_val(i)
          !
          do m = 1,nq
            upt(m) = upt(m) + a * usp(n0+j,m)
            do l = 1,nr
              dupt(l,m) = dupt(l,m) + a * dusp(n0+j,l,m)
            end do
          end do
          !
        end do
        !
        edgeusp_dt(ne)%v(1:nq,kf) = edgeusp_dt(ne)%v(1:nq,kf) + ave*upt(1:nq)
        !
        edgedusp_dt(ne)%v(1:nr,1:nq,kf) = edgedusp_dt(ne)%v(1:nr,1:nq,kf) + &
                                       ave*dupt(1:nr,1:nq)
        !
      end do edge_pts
      !
    end do cell_edges
    !
    ! Interpolate the solution and gradients to the nodes
    ! of the cell and add to the nodeusp and nodedusp arrays
    !
    cell_nodes: do n = 1,size(cell(nc)%node)
      !
      nn  = cell(nc)%node(n)%idx
      ave = cell(nc)%node(n)%ave_coef
      lf  = n
      !
      upt(1:nq) = zero
      dupt(1:nr,1:nq) = zero
      !
      do i = interp(this_geom,this_order)%toNode%row_ptr(lf)+1, &
             interp(this_geom,this_order)%toNode%row_ptr(lf+1)
        !
        j = interp(this_geom,this_order)%toNode%row_idx(i)
        a = interp(this_geom,this_order)%toNode%row_val(i)
        !
        do m = 1,nq
          upt(m) = upt(m) + a * usp(n0+j,m)
          do l = 1,nr
            dupt(l,m) = dupt(l,m) + a * dusp(n0+j,l,m)
          end do
        end do
        !
      end do
      !
      nodeusp_dt(nn)%v(1:nq) = nodeusp_dt(nn)%v(1:nq) + ave*upt(1:nq)
      !
      nodedusp_dt(nn)%v(1:nr,1:nq) = nodedusp_dt(nn)%v(1:nr,1:nq) + &
                                  ave*dupt(1:nr,1:nq)
      !
    end do cell_nodes
    !
  end do cell_loop
  !
  !
  !
  if (ncpu > 1) then
    call exch_edge_and_node_solutions(nodeusp_dt,nodedusp_dt, &
                                      edgeusp_dt,edgedusp_dt)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_edge_and_node_solutions
!
!###############################################################################
!
subroutine create_cell_subconnectivity(my_cells_beg,my_cells_end,ipelem)
  !
  !.. Use Statements ..
  use order_mod, only : gmin => n_min_geom
  use order_mod, only : gmax => n_max_geom
  use order_mod, only : omin => n_min_order
  use order_mod, only : omax => n_max_order
  use order_mod, only : geom_solpts
  use geovar,    only : n_solpts,ncell,cell
  use geovar,    only : pst_geom,nr,fp
  use flowvar,   only : faceusp_dt
  use flowvar,   only : edgeusp_dt,edgedusp_dt
  use flowvar,   only : nodeusp_dt,nodedusp_dt
  use ovar,      only : postproc_debug_grid
  !
  use parallel_mod, only : exch_connectivity
  !
  !.. Formal Arguments ..
  integer(CST),                              intent(inout) :: my_cells_beg
  integer(CST),                              intent(inout) :: my_cells_end
  integer(CST), allocatable, dimension(:,:), intent(inout) :: ipelem
  !
  !.. Local Scalars ..
  integer :: nop,nsp,nfp,nep,nnp
  integer :: mop,msp,mfp,mep,mnp
  integer :: n,m,nc,ne,nf,nn,n1,n2,np
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: mxpt,my_pts_offset,edge_order
  integer :: this_subcells,n_my_output_cells
  integer :: bad_indices,sf,kf,lf
  integer :: ifc,iec,inc,ierr,nop1,nop2
  integer :: n_total_output_cells
  integer :: n_total_output_pts
  integer :: running_offset
  integer :: ibeg,iend,i,j
  integer :: loc_n_my_output_pts
  integer :: loc_my_pts_beg
  integer :: loc_my_pts_end
  integer :: debug_unit
  !
  integer(int_mpi) :: this_cpu
  !
  character(len=100) :: array_name
  !
  !.. Local Arrays ..
  integer, dimension(1:ncpu) :: cpu_mxpts
  integer, dimension(1:ncpu) :: cpu_subcells
  !
  !.. Local Allocatable Arrays ..
  type(geom_pnt_map_t), allocatable :: map(:,:)
  !
  integer, allocatable :: cell_points_ptr(:)
  integer, allocatable :: cell_points(:)
  integer, allocatable :: tmp_cell_points(:)
  integer, allocatable :: cell_pts(:)
  integer, allocatable :: rev_pts(:)
  integer, allocatable :: face_pt(:)
  integer, allocatable :: edge_pt(:)
  integer, allocatable :: node_pt(:)
  integer, allocatable :: remaining_pts(:)
  !
  logical(lk), allocatable :: my_faces(:)
  logical(lk), allocatable :: my_edges(:)
  logical(lk), allocatable :: my_nodes(:)
  logical(lk), allocatable :: node_already_done(:)
  logical(lk), allocatable :: edge_already_done(:)
  logical(lk), allocatable :: msk(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_cell_subconnectivity"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  debug_unit = 160
  !
  ! Initialize the maximum number of output points to zero for all processors
  !
  cpu_mxpts = 0
  !
  allocate ( map(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"map",1,__LINE__,__FILE__,ierr,error_message)
  !
  !##########################################################################
  !##########################################################################
  !###                                                                    ###
  !###  1. First create the integer pointer array for the indices of all  ###
  !###     the points within each cell on this partition.                 ###
  !###                                                                    ###
  !##########################################################################
  !##########################################################################
  !
  allocate ( cell_points_ptr(1:ncell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_points_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! Loop through each cell and count the number of output points
  ! for each cell on this partition
  !
  msp = 0
  mfp = 0
  mep = 0
  mnp = 0
  mop = 0
  !
  do nc = 1,ncell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Create the mapping arrays if they havent already been created for
    ! this combination of cell geometry and solution order
    !
    if (.not. allocated(map(this_geom,this_order)%sp)) then
      call map_to_output_points(this_geom,this_order,map(this_geom,this_order))
    end if
    !
    ! Total number of solution points for this cell
    !
    nsp = size( map(this_geom,this_order)%sp )
    msp = max(msp,nsp)
    !
    ! Total number of face points for this cell
    !
    nfp = size( map(this_geom,this_order)%fp )
    mfp = max(mfp,nfp)
    !
    ! Total number of edge points for this cell
    !
    nep = size( map(this_geom,this_order)%ep )
    mep = max(mep,nep)
    !
    ! Total number of node points for this cell
    !
    nnp = size( map(this_geom,this_order)%np )
    mnp = max(mnp,nnp)
    !
    ! Total number of output points for this cell
    !
    nop = nsp + nfp + nep + nnp
    mop = max(mop,nop)
    !
    ! Add the total number of output points for this cell to
    ! the cell_points_ptr array
    !
    cell_points_ptr(nc+1) = cell_points_ptr(nc) + nop
    !
  end do
  !
  ! Save the total number of output points for this processor
  !
  cpu_mxpts(mypnum+1) = last(cell_points_ptr)
  !
  ! Collect the total number of output points from all processors
  !
  if (ncpu > 1) then
   !call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_DATATYPE_NULL, &
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum,"before mpi_allgather",__LINE__,__FILE__
#endif
    call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_INTEGER, &
                       cpu_mxpts,1_int_mpi,mpi_inttyp, &
                       MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum," after mpi_allgather",__LINE__,__FILE__
#endif
  end if
  !
  ! Create the point offset for this processor
  !
  my_pts_offset = sum( cpu_mxpts(1:mypnum) )
  !
  !#########################################################################
  !#########################################################################
  !###                                                                   ###
  !###  2. Perform the necessary allocations and initializations before  ###
  !###     creating the initial point indices for each cell on this      ###
  !###     partition.                                                    ###
  !###                                                                   ###
  !#########################################################################
  !#########################################################################
  !
  ! Allocate the cell_points array which will contain the indices for the output
  ! points within each cell. This will initially contain the indices for each
  ! cell with all colocated node, edge, and face points on partition boundaries
  ! having a unique index from each cell containing the node, edge, or face.
  ! After exchanging the unique indices between adjacent partitions, the indices
  ! for these colocated points will be reduced to a single index that is shared
  ! by all cells that contain the point.
  !
  allocate ( cell_points(1:last(cell_points_ptr)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_points",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate the arrays to mark nodes and edges that have already been completed
  !
  nn = size(nodeusp_dt)
  allocate ( node_already_done(1:nn) , source=fals , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"node_already_done",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ne = size(edgeusp_dt)
  allocate ( edge_already_done(1:ne) , source=fals , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"edge_already_done",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! Allocate the cell_pts array to temporarily store the point indices of
  ! the current cell thats being operated on
  !
  allocate ( cell_pts(1:mop) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_pts",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( face_pt(1:mfp) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_pt",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( edge_pt(1:mep) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"edge_pt",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( node_pt(1:mnp) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"node_pt",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Zero the nodeusp, nodedusp, edgeusp, edgedusp, and faceusp array
  ! before we start storing point indices in these arrays
  !
#ifdef SPECIAL_FOR_PGI
  call nodeusp_dt%zero_ne1d_t
  call nodedusp_dt%zero_ne2d_t
  call edgeusp_dt%zero_ne2d_t
  call edgedusp_dt%zero_ne3d_t
  call faceusp_dt%zero_face_var
#else
  call nodeusp_dt%zero
  call nodedusp_dt%zero
  call edgeusp_dt%zero
  call edgedusp_dt%zero
  call faceusp_dt%zero
#endif
  !
  !###########################################################################
  !###########################################################################
  !###                                                                     ###
  !###  3. Loop through the cells on this partition and find the indices   ###
  !###     for the output points of each cell and use the node and edge    ###
  !###     mask arrays that were just allocated to make sure that we only  ###
  !###     assign a point index for the first node or edge point found     ###
  !###     for the nodes and edges with multiple colocated points and      ###
  !###     that are only found on this partition.                          ###
  !###                                                                     ###
  !###########################################################################
  !###########################################################################
  !
  pack_cell_pts: do nc = 1,ncell
    !
    ! First initialize cell_pts to 0
    !
    cell_pts = 0
    face_pt  = 0
    edge_pt  = 0
    node_pt  = 0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    nop1 = cell_points_ptr(nc)+1
    nop2 = cell_points_ptr(nc+1)
    nop  = nop2-nop1+1
    !
    ! Total number of solution points for this cell
    !
    nsp = size( map(this_geom,this_order)%sp )
    !
    ! Find the total number of face points for this cell and create the
    ! map from the face point lf to the correct location in cell_pts
    !
    nfp = size( map(this_geom,this_order)%fp )
    !
    face_pt(1:nfp) = nsp + intseq(1,nfp)
    !
    ! Find the total number of edge points for this cell and create the
    ! map from the edge point lf to the correct location in cell_pts
    !
    nep = size( map(this_geom,this_order)%ep )
    !
    edge_pt(1:nep) = nsp + nfp + intseq(1,nep)
    !
    ! Find the total number of node points for this cell and create the
    ! map from the node point lf to the correct location in cell_pts
    !
    nnp = size( map(this_geom,this_order)%np )
    !
    node_pt(1:nnp) = nsp + nfp + nep + intseq(1,nnp)
    !
    ! Create the initial indices for the solution points in this cell
    !
    cell_pts(1:nsp) = nop1 + map(this_geom,this_order)%sp
    !
    ! Create the initial indices for the face points of this cell
    !
    cell_pts(face_pt) = nop1 + map(this_geom,this_order)%fp
    !
    ! Create the initial indices for the edge points of this cell
    !
    cell_pts(edge_pt) = nop1 + map(this_geom,this_order)%ep
    !
    ! Create the initial indices for the node points of this cell
    !
    cell_pts(node_pt) = nop1 + map(this_geom,this_order)%np
    !
    ! Adjust cell_pts for the offset of this processor
    !
    cell_pts(1:nop) = cell_pts(1:nop) + my_pts_offset
    !
    ! Pack up the face point indices of this cell into faceusp
    !
    pack_cell_face_pts: do n = 1,size(cell(nc)%face)
      !
      nf = cell(nc)%face(n)%idx
      sf = cell(nc)%face(n)%side
      !
      face_geom = cell(nc)%face(n)%geom
      face_order = cell(nc)%face(n)%order
      !
      do kf = 1,geom_solpts(face_geom,face_order)
        !
        lf = cell(nc)%face(n)%fp_offset + fp%geord(face_geom,face_order)% &
                                          rot(cell(nc)%face(n)%rotation)% &
                                          pts(kf)
        !
        faceusp_dt(nf)%v(1,kf,sf) = real( cell_pts( face_pt(lf) ) , kind=wp )
        !
      end do
      !
    end do pack_cell_face_pts
    !
    ! Pack up the edge point indices of this cell into edgeusp
    !
    pack_cell_edge_pts: do n = 1,size(cell(nc)%edge)
      !
      ne = cell(nc)%edge(n)%idx
      !
      if (edge_already_done(ne)) cycle pack_cell_edge_pts
      !
      edge_order = cell(nc)%edge(n)%order
      !
      do kf = 1,edge_order+1
        !
        lf = cell(nc)%edge(n)%ep_offset + merge( edge_order+2-kf , kf , &
                                                cell(nc)%edge(n)%rev_dir )
        !
        edgeusp_dt(ne)%v(1,kf) = real( cell_pts( edge_pt(lf) ) , kind=wp )
        !
      end do
      !
      edge_already_done(ne) = true
      !
    end do pack_cell_edge_pts
    !
    ! Pack up the node point indices of this cell into nodeusp
    !
    pack_cell_node_pts: do n = 1,size(cell(nc)%node)
      !
      nn = cell(nc)%node(n)%idx
      !
      if (node_already_done(nn)) cycle pack_cell_node_pts
      !
      lf = n
      !
      nodeusp_dt(nn)%v(1) = real( cell_pts( node_pt(lf) ) , kind=wp )
      !
      node_already_done(nn) = true
      !
    end do pack_cell_node_pts
    !
    ! Save the output points for this cell, in correct order, and
    ! before any corrections for duplicates
    !
    cell_points(nop1:nop2) = cell_pts(1:nop)
    !
  end do pack_cell_pts
  !
  if (allocated(node_already_done)) then
    deallocate ( node_already_done , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"node_already_done",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (allocated(edge_already_done)) then
    deallocate ( edge_already_done , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edge_already_done",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (ncpu == 1) then
    ! debug_unit = 160
    call write_cell_points(debug_unit,omin+1,cell_points_ptr,cell_points)
  end if
  !
  !###########################################################################
  !###########################################################################
  !###                                                                     ###
  !###  4. If using multiple partitions/processors, exchange the point     ###
  !###     indices for any nodes, edges, and faces, that are on partition  ###
  !###     boundaries.                                                     ###
  !###                                                                     ###
  !###########################################################################
  !###########################################################################
  !
  if (ncpu > 1) then
    call exch_connectivity(nodeusp_dt,nodedusp_dt, &
                           edgeusp_dt,edgedusp_dt,faceusp_dt)
  end if
  !
  !##########################################################################
  !##########################################################################
  !###                                                                    ###
  !###  5. Loop through the cells on this partition and unpack the point  ###
  !###     indices from the nodeusp, edgeusp, and faceusp arrays.         ###
  !###                                                                    ###
  !##########################################################################
  !##########################################################################
  !
  do nc = 1,ncell
    !
    cell_pts = 0
    face_pt  = 0
    edge_pt  = 0
    node_pt  = 0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    nop1 = cell_points_ptr(nc)+1
    nop2 = cell_points_ptr(nc+1)
    nop  = nop2-nop1+1
    !
    ! Total number of solution points for this cell
    !
    nsp = size( map(this_geom,this_order)%sp )
    !
    ! Find the total number of face points for this cell and create the
    ! map from the face point lf to the correct location in cell_pts
    !
    nfp = size( map(this_geom,this_order)%fp )
    !
    face_pt(1:nfp) = nsp + intseq(1,nfp)
    !
    ! Find the total number of edge points for this cell and create the
    ! map from the edge point lf to the correct location in cell_pts
    !
    nep = size( map(this_geom,this_order)%ep )
    !
    edge_pt(1:nep) = nsp + nfp + intseq(1,nep)
    !
    ! Find the total number of node points for this cell and create the
    ! map from the node point lf to the correct location in cell_pts
    !
    nnp = size( map(this_geom,this_order)%np )
    !
    node_pt(1:nnp) = nsp + nfp + nep + intseq(1,nnp)
    !
    ! Copy the output points for this cell, in correct order, and
    ! before any corrections for duplicates
    !
    cell_pts(1:nop) = cell_points(nop1:nop2)
    !
    ! Unpack the face point indices of this cell from faceusp
    !
    unpack_cell_face_pts: do n = 1,size(cell(nc)%face)
      !
      nf = cell(nc)%face(n)%idx
      !
      face_geom = cell(nc)%face(n)%geom
      face_order = cell(nc)%face(n)%order
      !
      do kf = 1,geom_solpts(face_geom,face_order)
        !
        lf = cell(nc)%face(n)%fp_offset + fp%geord(face_geom,face_order)% &
                                          rot(cell(nc)%face(n)%rotation)% &
                                          pts(kf)
        !
        cell_pts( face_pt(lf) ) = nint( faceusp_dt(nf)%v(1,kf,1) , &
                                        kind=kind(cell_pts) )
        !
      end do
      !
    end do unpack_cell_face_pts
    !
    ! Unpack the edge point indices of this cell from edgeusp
    !
    unpack_cell_edge_pts: do n = 1,size(cell(nc)%edge)
      !
      ne = cell(nc)%edge(n)%idx
      edge_order = cell(nc)%edge(n)%order
      !
      do kf = 1,edge_order+1
        !
        lf = cell(nc)%edge(n)%ep_offset + merge( edge_order+2-kf , kf , &
                                                cell(nc)%edge(n)%rev_dir)
        !
        cell_pts( edge_pt(lf) ) = nint( edgeusp_dt(ne)%v(1,kf) , &
                                        kind=kind(cell_pts) )
        !
      end do
      !
    end do unpack_cell_edge_pts
    !
    ! Unpack the node point indices of this cell from nodeusp
    !
    unpack_cell_node_pts: do n = 1,size(cell(nc)%node)
      !
      nn = cell(nc)%node(n)%idx
      !
      lf = n
      !
      cell_pts( node_pt(lf) ) = nint( nodeusp_dt(nn)%v(1), kind=kind(cell_pts) )
      !
    end do unpack_cell_node_pts
    !
    ! Copy the unpacked point indices into cell_points
    !
    cell_points(nop1:nop2) = cell_pts(1:nop)
    !
  end do
  !
  if (ncpu == 1) then
    ! debug_unit = 161
    call write_cell_points(debug_unit,omin+1,cell_points_ptr,cell_points)
  end if
  !
  !#############################################################################
  !#############################################################################
  !###                                                                       ###
  !###  6. Find the nodes, edges, and faces of the grid that this processor  ###
  !###     is responsible for writing to the CGNS solution file.             ###
  !###                                                                       ###
  !#############################################################################
  !#############################################################################
  !
  nf = size(faceusp_dt)
  allocate ( my_faces(1:nf) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"my_faces",1,__LINE__,__FILE__,ierr,error_message)
  !
  ne = size(edgeusp_dt)
  allocate ( my_edges(1:ne) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"my_edges",1,__LINE__,__FILE__,ierr,error_message)
  !
  nn = size(nodeusp_dt)
  allocate ( my_nodes(1:nn) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"my_nodes",1,__LINE__,__FILE__,ierr,error_message)
  !
  bad_indices = 0
  !
  ! Find the face points for which this processor is responsible
  !
  find_my_faces: do nf = 1,size(faceusp_dt)
    !
    n = size(faceusp_dt(nf)%v,dim=2)
    m = count(faceusp_dt(nf)%v(1,:,1) > zero)
    !
    if (n == m) then
      !
      my_faces(nf) = true
      !
    else if (m > 0) then
      !
      if (bad_indices == 0) write (iout,11)
      !
      bad_indices = bad_indices + 1
      write (iout,12) "Face",nf,m,n
      !
    end if
    !
  end do find_my_faces
  !
  ! Find the edge points for which this processor is responsible
  !
  find_my_edges: do ne = 1,size(edgeusp_dt)
    !
    n = size(edgeusp_dt(ne)%v,dim=2)
    m = count(edgeusp_dt(ne)%v(1,:) > zero)
    !
    if (n == m) then
      !
      my_edges(ne) = true
      !
    else if (m > 0) then
      !
      if (bad_indices == 0) write (iout,11)
      !
      bad_indices = bad_indices + 1
      write (iout,12) "Edge",ne,m,n
      !
    end if
    !
  end do find_my_edges
  !
  ! Find the node points for which this processor is responsible
  !
  find_my_nodes: do nn = 1,size(nodeusp_dt)
    !
    if (nodeusp_dt(nn)%v(1) > zero) then
      my_nodes(nn) = true
    end if
    !
  end do find_my_nodes
  !
  ! Reduce the maximum value of bad_indices across all processors
  !
  if (ncpu > 1) then
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum,"before mpi_allreduce",__LINE__,__FILE__
#endif
    call mpi_allreduce(MPI_IN_PLACE,bad_indices,1_int_mpi,mpi_inttyp, &
                       MPI_MAX,MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum," after mpi_allreduce",__LINE__,__FILE__
#endif
  end if
  !
  ! Abort execution if bad_indices was non-zero on any processor
  !
  if (bad_indices /= 0) then
    write (error_message,13)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  !############################################################################
  !############################################################################
  !###                                                                      ###
  !###  7. Since all the colocated points with multiple unique indices for  ###
  !###     the same point were reduced to a single index, there is now a    ###
  !###     large number of point indices which are no longer used,          ###
  !###     creating 'holes' in the indexing. Now we will go through and     ###
  !###     remove these holes in the indexing so the point indices are a    ###
  !###     continuous interval.                                             ###
  !###                                                                      ###
  !############################################################################
  !############################################################################
  !
  ! Find the maximum point index
  !
  mxpt = maxval(abs(cell_points))
  !
  ! Collectively find the maximum point index from all processors
  !
  if (ncpu > 1) then
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum,"before mpi_allreduce",__LINE__,__FILE__
#endif
    call mpi_allreduce(MPI_IN_PLACE,mxpt,1_int_mpi,mpi_inttyp, &
                       MPI_MAX,MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum," after mpi_allreduce",__LINE__,__FILE__
#endif
  end if
  !
  ! Allocate the array to mark the points on each processor that still remain
  !
  allocate ( remaining_pts(1:mxpt) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"remaining_pts",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Mark all the remaining points on this processor
  ! NOTE: We need to use the absolute value of cell_points because some
  !       of the indices are still negative to mark points for which an
  !       adjacent processor is responsible for writing the solution at
  !       the point.
  !
  remaining_pts( abs(cell_points(:)) ) = 1
!#define DEBUG_CREATE_CELL_SUBCONNECTIVITY
#ifdef DEBUG_CREATE_CELL_SUBCONNECTIVITY
  write (iout,100) "                     mxpt = ",mxpt
  write (iout,100) "count(  cell_points <= 0) = ",count(cell_points<=0)
  write (iout,100) "count(remaining_pts == 1) = ",count(remaining_pts==1)
#endif
  !
  ! Collectively reduce the remaining points on all processors
  !
  if (ncpu > 1) then
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum,"before mpi_allreduce",__LINE__,__FILE__
#endif
    call mpi_allreduce(MPI_IN_PLACE,remaining_pts,int(mxpt,kind=int_mpi), &
                       mpi_inttyp,MPI_MAX,MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum," after mpi_allreduce",__LINE__,__FILE__
#endif
  end if
  !
  ! Convert remaining_pts into an array that maps a point index in cell_points
  ! into the equivalent point index if all non-used point indices were removed
  !
  mxpt = 0
  do n = 1,size(remaining_pts)
    if (remaining_pts(n) /= 0) then
      mxpt = mxpt + 1
      remaining_pts(n) = mxpt
    end if
  end do
  !
  ! Save the total number of output points across all processors
  !
  n_total_output_pts = mxpt
  !
  ! Now redo the point indices in cell_points to account
  ! for only the remaining point indices
  ! NOTE: We need to make sure that we keep the negative sign for indices
  !       within cell_points that mark points for which an adjacent partition
  !       is responsible for writing the solution at the point.
  !
  do n = 1,size(cell_points)
    cell_points(n) = sign( remaining_pts( abs(cell_points(n)) ) , &
                           cell_points(n) )
  end do
  !
  if (ncpu == 1) then
    ! debug_unit = 162
    call write_cell_points(debug_unit,omin+1,cell_points_ptr,cell_points)
  end if
  !
  allocate ( msk(1:n_total_output_pts) , source=fals , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",1,__LINE__,__FILE__,ierr,error_message)
  !
 !where (cell_points > 0) msk(cell_points) = true
  do n = 1,size(cell_points)
    if (cell_points(n) > 0) then
      msk(cell_points(n)) = true
    end if
  end do
  loc_n_my_output_pts = count(msk)
  !
  if (ncpu > 1) then
    !
    ! Reallocate remaining_pts to the new maximum point index
    ! and reset it to 0
    !
    call reallocate( remaining_pts , n_total_output_pts )
    remaining_pts(:) = 0
    !
    ! USING F2003 AUTO-REALLOCATION
    tmp_cell_points = cell_points
    !
    running_offset = 0
    !
    do n = 1,ncpu
      !
      this_cpu = int(n-1,kind=int_mpi)
      !
      ! Reset remaining_pts to 0 on all processors
      !
      remaining_pts(:) = 0
      !
      !
      if (mypnum == this_cpu) then
        ibeg = running_offset + 1
        iend = running_offset + loc_n_my_output_pts
        loc_my_pts_beg = ibeg
        loc_my_pts_end = iend
        remaining_pts(:) = unpack( intseq(ibeg,iend) , mask=msk , field=0 )
      end if
      !
#ifdef VERBOSE_MPI
      write (iout,dbg_mpi_frmt) mypnum,"before mpi_bcast",__LINE__,__FILE__
#endif
      call mpi_bcast(remaining_pts,size(remaining_pts,kind=int_mpi), &
                     mpi_inttyp,this_cpu,MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
      write (iout,dbg_mpi_frmt) mypnum," after mpi_bcast",__LINE__,__FILE__
#endif
      !
      if (mypnum == this_cpu) then
        do i = 1,size(cell_points)
          if (cell_points(i) > 0) then
            tmp_cell_points(i) = remaining_pts( cell_points(i) )
          end if
        end do
      else
        do i = 1,size(cell_points)
          if (cell_points(i) < 0) then
            j = abs(cell_points(i))
            if (remaining_pts(j) > 0) then
              tmp_cell_points(i) = -remaining_pts(j)
            end if
          end if
        end do
      end if
      !
      ! Update the running offset needed to get the
      ! starting point index for the next processor
      !
      running_offset = running_offset + count(remaining_pts /= 0)
      !
    end do
    !
    ! Transfer tmp_cell_points to cell_points
    ! NOTE: This should leave tmp_cell_pts deallocated
    !
    call move_alloc( from=tmp_cell_points , to=cell_points )
    !
    remaining_pts(:) = 0
   !where (cell_points > 0) remaining_pts(cell_points) = cell_points
    do n = 1,size(cell_points)
      if (cell_points(n) > 0) then
        remaining_pts(cell_points(n)) = cell_points(n)
      end if
    end do
    !
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum,"before mpi_barrier",__LINE__,__FILE__
#endif
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum," after mpi_barrier",__LINE__,__FILE__
#endif
    do n = 1,ncpu
      !
      flush (iout)
#ifdef VERBOSE_MPI
      write (iout,dbg_mpi_frmt) mypnum,"before mpi_barrier",__LINE__,__FILE__
#endif
      call mpi_barrier(MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
      write (iout,dbg_mpi_frmt) mypnum," after mpi_barrier",__LINE__,__FILE__
#endif
      if (mypnum == int(n-1,kind=int_mpi)) then
        write (iout,31) mypnum,minval(remaining_pts)
        write (iout,32) mypnum,maxval(remaining_pts)
        write (iout,33) mypnum,maxval(remaining_pts)-minval(remaining_pts)+1
        write (iout,34) mypnum,count(remaining_pts /= 0)
        write (iout,35) mypnum,loc_my_pts_beg
        write (iout,36) mypnum,loc_my_pts_end
        write (iout,37) mypnum,loc_n_my_output_pts
        if (n == ncpu) write (iout,38)
      end if
      flush (iout)
#ifdef VERBOSE_MPI
      write (iout,dbg_mpi_frmt) mypnum,"before mpi_barrier",__LINE__,__FILE__
#endif
      call mpi_barrier(MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
      write (iout,dbg_mpi_frmt) mypnum," after mpi_barrier",__LINE__,__FILE__
#endif
      !
    end do
    !
  end if
  !
  if (allocated(msk)) then
    deallocate ( msk , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(remaining_pts)) then
    deallocate ( remaining_pts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"remaining_pts",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  ! TEMPORARY DEBUG STUFF
  ! TEMPORARY DEBUG STUFF
  ! TEMPORARY DEBUG STUFF
  if (ncpu == 1) then
    nsp = (omin+1)**3
    nfp = (omin+1)**2
    nep = (omin+1)
    loc_n_my_output_pts = ncell*nsp + size(faceusp_dt)*nfp + &
                                      size(edgeusp_dt)*nep + &
                                      size(nodeusp_dt)
    write (iout,40) loc_n_my_output_pts
    ! debug_unit = 163
    call write_cell_points(debug_unit,omin+1,cell_points_ptr,cell_points)
  end if
  n1 = minval( cell_points , mask=(cell_points>0) )
  n2 = maxval( cell_points , mask=(cell_points>0) )
  allocate ( msk(n1:n2) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",1,__LINE__,__FILE__,ierr,error_message)
  do n = 1,size(cell_points)
    if (cell_points(n) > 0) then
      msk(cell_points(n)) = true
    end if
  end do
  loc_n_my_output_pts = count(msk)
  deallocate ( msk , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
  write (iout,41) loc_n_my_output_pts
  ! TEMPORARY DEBUG STUFF
  ! TEMPORARY DEBUG STUFF
  ! TEMPORARY DEBUG STUFF
  !
  !###########################################################################
  !###########################################################################
  !###                                                                     ###
  !###  8. Create the arrays that map all the points that the current      ###
  !###     processor is responsible for writing to the CGNS solution file  ###
  !###     to the correct location within the array that collects the      ###
  !###     plot function values at each of these points and is then sent   ###
  !###     into the CGNS function for writing the solution.                ###
  !###                                                                     ###
  !###########################################################################
  !###########################################################################
  !
  allocate ( output_sol_pts(1:n_solpts) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"output_sol_pts",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  nf = count( my_faces )
  allocate ( output_face_pts(1:nf) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"output_face_pts",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ne = count( my_edges )
  allocate ( output_edge_pts(1:ne) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"output_edge_pts",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  nn = count( my_nodes )
  allocate ( output_node_pts(1:nn) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"output_node_pts",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  cpu_subcells(:) = 0
  !
  ifc = 0
  iec = 0
  inc = 0
  !
  do nc = 1,ncell
    !
    cell_pts = 0
    face_pt  = 0
    edge_pt  = 0
    node_pt  = 0
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Get the number of subcells for this cell geometry at this order.
    ! The number of subcells is basically:
    !     (this_order+2)**(spatial dimension of this_geom)
    !
    this_subcells = (this_order+2)**geom_dimen(this_geom)
    !
    ! Add the number of subcells for this cell geometry to the running
    ! count of total number of output subcells for this processor
    !
    cpu_subcells(mypnum+1) = cpu_subcells(mypnum+1) + this_subcells
    !
    ! Get the beginning and ending indices for the output point indices
    !
    nop1 = cell_points_ptr(nc)+1
    nop2 = cell_points_ptr(nc+1)
    nop  = nop2-nop1+1
    !
    ! Total number of solution points for this cell
    !
    nsp = size( map(this_geom,this_order)%sp )
    !
    ! Find the total number of face points for this cell and create the
    ! map from the face point lf to the correct location in cell_pts
    !
    nfp = size( map(this_geom,this_order)%fp )
    !
    face_pt(1:nfp) = nsp + intseq(1,nfp)
    !
    ! Find the total number of edge points for this cell and create the
    ! map from the edge point lf to the correct location in cell_pts
    !
    nep = size( map(this_geom,this_order)%ep )
    !
    edge_pt(1:nep) = nsp + nfp + intseq(1,nep)
    !
    ! Find the total number of node points for this cell and create the
    ! map from the node point lf to the correct location in cell_pts
    !
    nnp = size( map(this_geom,this_order)%np )
    !
    node_pt(1:nnp) = nsp + nfp + nep + intseq(1,nnp)
    !
    ! Make a copy of the output points for this cell
    !
    cell_pts(1:nop) = cell_points(nop1:nop2)
    !
    ! Get the output point indices for the solution points of this cell
    !
    output_sol_pts(n1:n2) = cell_pts(1:nsp)
    !
    ! Get the output point indices for the face points of this cell
    !
    create_output_face_pts: do n = 1,size(cell(nc)%face)
      !
      nf = cell(nc)%face(n)%idx
      !
      ! my_faces(nf) will be true if the current processor is responsible
      ! for the solution at the points on this face AND this face hasnt
      ! already been added to output_face_pts
      !
      if (.not. my_faces(nf)) cycle create_output_face_pts
      !
      ! Increment the face counter by one and then set
      ! my_faces(nf) to false so that this face is not
      ! found a second time by the other cell on this face
      !
      ifc = ifc + 1
      my_faces(nf) = fals
      !
      ! Save this face index to output_face_pts
      !
      output_face_pts(ifc)%face_idx = nf
      !
      ! Get the number of face points for this face
      !
      face_geom = cell(nc)%face(n)%geom
      face_order = cell(nc)%face(n)%order
      !
      nfp = geom_solpts(face_geom,face_order)
      !
      ! Allocate the out_idx component of output_face_pts for this face
      !
      allocate ( output_face_pts(ifc)%out_idx(1:nfp) , source=0 , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) "face",ifc
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      do kf = 1,nfp
        !
        lf = cell(nc)%face(n)%fp_offset + fp%geord(face_geom,face_order)% &
                                          rot(cell(nc)%face(n)%rotation)% &
                                          pts(kf)
        !
        output_face_pts(ifc)%out_idx(kf) = cell_pts( face_pt(lf) )
        !
      end do
      !
     !if (any(n == [Xmin3D_Face,Ymax3D_Face,Zmin3D_Face])) then
     !  output_face_pts(ifc)%out_idx = &
     !                transpose_face(output_face_pts(ifc)%out_idx)
     !end if
      !
    end do create_output_face_pts
    !
    ! Get the output point indices for the edge points of this cell
    !
    create_output_edge_pts: do n = 1,size(cell(nc)%edge)
      !
      ne = cell(nc)%edge(n)%idx
      !
      ! my_edges(ne) will be true if the current processor is responsible
      ! for the solution at the points on this edge AND this edge hasnt
      ! already been added to output_edge_pts
      !
      if (.not. my_edges(ne)) cycle create_output_edge_pts
      !
      ! Increment the edge counter by one and then set
      ! my_edges(ne) to false so that this edge is not
      ! found a second time by other cells with this edge
      !
      iec = iec + 1
      my_edges(ne) = fals
      !
      ! Save this edge index to output_edge_pts
      !
      output_edge_pts(iec)%edge_idx = ne
      !
      ! Get the number of edge points for this edge
      !
      edge_order = cell(nc)%edge(n)%order
      nep = edge_order+1
      !
      ! Allocate the out_idx component of output_edge_pts for this edge
      !
      allocate ( output_edge_pts(iec)%out_idx(1:nep) , source=0 , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) "edge",iec
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      do kf = 1,nep
        !
        lf = cell(nc)%edge(n)%ep_offset + merge( edge_order+2-kf , kf , &
                                                cell(nc)%edge(n)%rev_dir)
        !
        output_edge_pts(iec)%out_idx(kf) = cell_pts( edge_pt(lf) )
        !
      end do
      !
    end do create_output_edge_pts
    !
    ! Get the output point indices for the nodes of this cell
    !
    create_output_node_pts: do n = 1,size(cell(nc)%node)
      !
      nn = cell(nc)%node(n)%idx
      !
      ! my_nodes(nn) will be true if the current processor is responsible
      ! for the solution at this node AND this node hasnt already been
      ! added to output_node_pts
      !
      if (.not. my_nodes(nn)) cycle create_output_node_pts
      !
      ! Increment the node counter by one and then set
      ! my_nodes(nn) to false so that this node is not
      ! found a second time by other cells with this node
      !
      inc = inc + 1
      my_nodes(nn) = fals
      !
      ! Save this node index to output_node_pts
      !
      output_node_pts(inc)%node_idx = nn
      !
      ! Save the output index for this node to output_node_pts
      !
      lf = n
      !
      output_node_pts(inc)%out_idx = cell_pts( node_pt(lf) )
      !
      !
    end do create_output_node_pts
    !
  end do
  !
  ! Check to make sure that the indices for all faces, edges, and nodes
  ! that this processor is responsible for were found and saved.
  !
  bad_indices = merge(1,0,any([my_faces,my_edges,my_nodes]))
  !
  ! Reduce the maximum value of bad_indices across all processors
  !
  if (ncpu > 1) then
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum,"before mpi_allreduce",__LINE__,__FILE__
#endif
    call mpi_allreduce(MPI_IN_PLACE,bad_indices,1_int_mpi,mpi_inttyp, &
                       MPI_MAX,MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum," after mpi_allreduce",__LINE__,__FILE__
#endif
  end if
  !
  ! Dump some error messages and abort execution
  ! if bad_indices was non-zero on any processor
  !
  if (bad_indices /= 0) then
    !
    if (mypnum == 0) write (iout,21)
    flush (iout)
    !
    ! Write an error if any faces for this processor were not found and saved
    !
    if (any(my_faces)) then
      write (iout,22) mypnum,count(my_faces),"Faces"
      flush (iout)
    end if
    !
    ! Write an error if any edges for this processor were not found and saved
    !
    if (any(my_edges)) then
      write (iout,22) mypnum,count(my_edges),"Edges"
      flush (iout)
    end if
    !
    ! Write an error if any nodes for this processor were not found and saved
    !
    if (any(my_nodes)) then
      write (iout,22) mypnum,count(my_nodes),"Nodes"
      flush (iout)
    end if
    !
    ! Abort execution if bad_indices was non-zero on any processor
    !
    write (error_message,23)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    !
  end if
  !
  ! Count the updated total number of output points for this processor
  !
  cpu_mxpts(:) = 0
  cpu_mxpts(mypnum+1) = size(output_sol_pts) + size(output_node_pts)
  do n = 1,size(output_face_pts)
    cpu_mxpts(mypnum+1) = cpu_mxpts(mypnum+1) + size(output_face_pts(n)%out_idx)
  end do
  do n = 1,size(output_edge_pts)
    cpu_mxpts(mypnum+1) = cpu_mxpts(mypnum+1) + size(output_edge_pts(n)%out_idx)
  end do
  !
  ! TEMPORARY DEBUG STUFF
  ! TEMPORARY DEBUG STUFF
  ! TEMPORARY DEBUG STUFF
  n_my_output_pts = cpu_mxpts(mypnum+1)
  write (iout,42) n_my_output_pts
  n_my_output_cells = cpu_subcells(mypnum+1)
  write (iout,43) n_my_output_cells
  ! TEMPORARY DEBUG STUFF
  ! TEMPORARY DEBUG STUFF
  ! TEMPORARY DEBUG STUFF
  !
  ! Collect the updated total number of output points
  ! and subcells from all processors
  !
  if (ncpu > 1) then
   !call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_DATATYPE_NULL, &
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum,"before mpi_allgather",__LINE__,__FILE__
#endif
    call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_INTEGER, &
                       cpu_mxpts,1_int_mpi,mpi_inttyp, &
                       MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum," after mpi_allgather",__LINE__,__FILE__
#endif
   !call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_DATATYPE_NULL, &
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum,"before mpi_allgather",__LINE__,__FILE__
#endif
    call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_INTEGER, &
                       cpu_subcells,1_int_mpi,mpi_inttyp, &
                       MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
    write (iout,dbg_mpi_frmt) mypnum," after mpi_allgather",__LINE__,__FILE__
#endif
  end if
  !
  my_pts_beg = sum( cpu_mxpts(1:mypnum) ) + 1
  my_pts_end = sum( cpu_mxpts(1:mypnum+1) )
  !
  my_cells_beg = sum( cpu_subcells(1:mypnum) ) + 1
  my_cells_end = sum( cpu_subcells(1:mypnum+1) )
  !
  n_total_output_pts = sum( cpu_mxpts )
  n_total_output_cells = sum( cpu_subcells )
  !
  ! Save n_total_output_pts and n_total_output_cells to zone_size
  !
  zone_size = [n_total_output_pts,n_total_output_cells,0]
  !
  ! n_my_output_pts: number of output points for this processor
  !
  n_my_output_pts = my_pts_end - my_pts_beg + 1
  !
  ! n_my_output_cells: number of output cell for this processor
  !
  n_my_output_cells = my_cells_end - my_cells_beg + 1
  !
  ! TEMPORARY DEBUG STUFF
  ! TEMPORARY DEBUG STUFF
  ! TEMPORARY DEBUG STUFF
  write (iout,44) n_my_output_pts
  write (iout,45) n_my_output_cells
  ! TEMPORARY DEBUG STUFF
  ! TEMPORARY DEBUG STUFF
  ! TEMPORARY DEBUG STUFF
  !
  !#######################################################
  !#######################################################
  !###                                                 ###
  !###  9.                                             ###
  !###                                                 ###
  !#######################################################
  !#######################################################
  !
  np = 2**nr
  !
  ! Create the connectivity array
  !
  if (allocated(ipelem)) then
    deallocate ( ipelem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ipelem",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( ipelem(1:np,1:n_my_output_cells) , source=0_CST , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ipelem",1,__LINE__,__FILE__,ierr,error_message)
  !
  if (postproc_debug_grid .and. ncpu == 1) then
    if (mod(n_my_output_cells,ncell) == 0) then
      n = n_my_output_cells / ncell
      allocate ( ipelem_of_cells(1:np,1:n,1:ncell) , source=0 , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"ipelem_of_cells",1,__LINE__,__FILE__,ierr, &
                       error_message)
    else
      write (error_message,'(a)') &
            "n_my_output_cells should be a multiple of ncell!"
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
  end if
  !
  allocate ( rev_pts(1:size(cell_pts)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"rev_pts",1,__LINE__,__FILE__,ierr,error_message)
  !
  inc = 0
  !
  do nc = 1,ncell
    !
    rev_pts(:) = 0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    nop1 = cell_points_ptr(nc)+1
    nop2 = cell_points_ptr(nc+1)
    nop  = nop2-nop1+1
    !
    ! Any points in this cell that another processor is responsible for
    ! will be identified with a negative value in cell_points, so take
    ! the absolute value of cell_points so that the subcell connectivities
    ! are correct when stored in ipelem.
    !
    cell_pts(1:nop) = abs( cell_points(nop1:nop2) )
    !
    ! Total number of solution points for this cell
    !
    nsp = size( map(this_geom,this_order)%sp )
    !
    ! Find the total number of face points for this cell and create the
    ! map from the face point lf to the correct location in cell_pts
    !
    nfp = size( map(this_geom,this_order)%fp )
    !
    face_pt(1:nfp) = nsp + intseq(1,nfp)
    !
    ! Find the total number of edge points for this cell and create the
    ! map from the edge point lf to the correct location in cell_pts
    !
    nep = size( map(this_geom,this_order)%ep )
    !
    edge_pt(1:nep) = nsp + nfp + intseq(1,nep)
    !
    ! Find the total number of node points for this cell and create the
    ! map from the node point lf to the correct location in cell_pts
    !
    nnp = size( map(this_geom,this_order)%np )
    !
    node_pt(1:nnp) = nsp + nfp + nep + intseq(1,nnp)
    !
    rev_pts(map(this_geom,this_order)%sp+1) = cell_pts(1:nsp)
    rev_pts(map(this_geom,this_order)%fp+1) = cell_pts(face_pt(1:nfp))
    rev_pts(map(this_geom,this_order)%ep+1) = cell_pts(edge_pt(1:nep))
    rev_pts(map(this_geom,this_order)%np+1) = cell_pts(node_pt(1:nnp))
    !
    ! Loop through the number of subcells for this cell and get the
    ! connectivity for each one using the output indices in cell_pts
    !
    do n = 1,size(pst_geom(this_geom,this_order)%connectivity,dim=2)
      !
      inc = inc + 1
      !
      ipelem(:,inc) = rev_pts(pst_geom(this_geom,this_order)% &
                                           connectivity(:,n) + 1)
      if (allocated(ipelem_of_cells)) then
        ipelem_of_cells(:,n,nc) = rev_pts(pst_geom(this_geom,this_order)% &
                                             connectivity(:,n) + 1)
      end if
      !
    end do
    !
  end do
  !
  !#######################################################
  !#######################################################
  !###                                                 ###
  !###  9. Deallocate the local arrays before leaving  ###
  !###                                                 ###
  !#######################################################
  !#######################################################
  !
  if (allocated(cell_points)) then
    deallocate ( cell_points , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_points",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(cell_points_ptr)) then
    deallocate ( cell_points_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_points_ptr",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (allocated(my_faces)) then
    deallocate ( my_faces , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"my_faces",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(my_edges)) then
    deallocate ( my_edges , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"my_edges",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(my_nodes)) then
    deallocate ( my_nodes , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"my_nodes",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(map)) then
    deallocate ( map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"map",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(cell_pts)) then
    deallocate ( cell_pts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_pts",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(face_pt)) then
    deallocate ( face_pt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"face_pt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edge_pt)) then
    deallocate ( edge_pt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edge_pt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(node_pt)) then
    deallocate ( node_pt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"node_pt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(rev_pts)) then
    deallocate ( rev_pts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"rev_pts",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1  format ("output_",a,"_pts(",i0,")%out_idx")
  11 format ("Error with inconsistent output point indices!")
  12 format (1x,a,1x,i0,1x,i0,1x,i0)
  13 format ("Error with inconsistent output point indices!")
  21 format ("Error: There appears to be faces, edges, or nodes, that",/, &
             "       were not found and saved on a processor!")
  22 format (4x,"on CPU ",i0," : ",i0," missed ",a)
  23 format ("There appears to be faces, edges, or nodes, that were not ", &
             "found and saved on a processor!")
  31 format (60("-"),/, &
             " (31) CPU-",i0," : minval(remaining_pts)     = ",i0)
  32 format (" (32) CPU-",i0," : maxval(remaining_pts)     = ",i0)
  33 format (" (33) CPU-",i0," : maxval - minval + 1       = ",i0)
  34 format (" (34) CPU-",i0," : count(remaining_pts /= 0) = ",i0)
  35 format (" (35) CPU-",i0," : loc_n_my_pts_beg          = ",i0)
  36 format (" (36) CPU-",i0," : loc_n_my_pts_end          = ",i0)
  37 format (" (37) CPU-",i0," : loc_n_my_output_pts       = ",i0)
  38 format (60("-"))
  40 format (" (40) Correct value : loc_n_my_output_pts   = ",i0)
  41 format (" (41) Using msk: loc_n_my_output_pts        = ",i0)
  42 format (" (42) Using cpu_mxpts: n_my_output_pts      = ",i0)
  43 format (" (43) Using cpu_subcells: n_my_output_cells = ",i0)
  44 format (" (44) Final n_my_output_pts                 = ",i0)
  45 format (" (45) Final n_my_output_cells               = ",i0)
  !
  100 format (1x,a,i0)
  !
end subroutine create_cell_subconnectivity
!
!###############################################################################
!
pure function transpose_face(face_pts) result(return_value)
  integer, dimension(:), intent(in) :: face_pts
  integer, dimension(1:size(face_pts)) :: return_value
  integer :: np1d
continue
  np1d = nint( sqrt(real(size(face_pts),kind=wp)) )
  return_value = reshape( &
                   transpose( reshape(face_pts,[np1d,np1d]) ) , &
                   shape(face_pts) )
end function transpose_face
!
!###############################################################################
!
elemental function plot_func_type_equals(a,b) result(return_value)
  !
  type(plot_func_t), intent(in) :: a
  type(plot_func_t), intent(in) :: b
  !
  logical(ldk) :: return_value
  !
continue
  !
  return_value = (a%func == b%func)
  !
end function plot_func_type_equals
!
!###############################################################################
!
pure function compute_output_variables(upt,gppt,dupt) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: upt(:)
  real(wp), intent(in) :: gppt(:)
  real(wp), intent(in) :: dupt(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(plot_var))
  !
  !.. Local Scalars ..
  integer :: i
  !
continue
  !
  do i = 1,size(plot_var)
    return_value(i) = compute_output_variable(upt,gppt,dupt,plot_var(i))
  end do
  !
end function compute_output_variables
!
!###############################################################################
!
pure function compute_output_variable(upt,gppt,dupt,plot_variable) &
                               result(return_value)
  !
  !.. Use Statements ..
  use ovar,    only : rhoref,aref,tref,time_ref,gam,pref_nd,itestcase
  !
  !.. Formal Arguments ..
  real(wp),          intent(in) :: upt(:)
  real(wp),          intent(in) :: gppt(:)
  real(wp),          intent(in) :: dupt(:,:)
  type(plot_func_t), intent(in) :: plot_variable
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: gamm
  !
continue
  !
  if (plot_variable == plot_static_density) then
    !
    return_value = upt(nec)                           * rhoref
    !
  else if (plot_variable == plot_velocity_x) then
    !
    return_value = upt(nmx) / upt(nec)                * aref
    !
  else if (plot_variable == plot_velocity_y) then
    !
    return_value = upt(nmy) / upt(nec)                * aref
    !
  else if (plot_variable == plot_velocity_z) then
    !
    return_value = upt(nmz) / upt(nec)                * aref
    !
  else if (plot_variable == plot_momentum_x) then
    !
    return_value = upt(nmx)                           * rhoref * aref
    !
  else if (plot_variable == plot_momentum_y) then
    !
    return_value = upt(nmy)                           * rhoref * aref
    !
  else if (plot_variable == plot_momentum_z) then
    !
    return_value = upt(nmz)                           * rhoref * aref
    !
  else if (plot_variable == plot_static_pressure) then
    !
    return_value = gppt(npres)                        * rhoref * aref * aref
    !
  else if (plot_variable == plot_static_temperature) then
    !
    return_value = gppt(ntemp)                        * tref
    !
  else if (plot_variable == plot_total_energy) then
    !
    return_value = upt(nee)                           * rhoref * aref * aref
    !
  else if (plot_variable == plot_mach_number) then
    !
    if (ncs == 0) then
      gamm = gam
    else
      gamm = gppt(ngamm)
    end if
    !
    return_value = norm2(upt(nmb:nme)) / sqrt( upt(nec)*gamm*gppt(npres) )
    !
  else if (plot_variable == plot_entropy) then
    !
    if (ncs == 0) then
      gamm = gam
    else
      gamm = gppt(ngamm)
    end if
    !
    if (itestcase == Inviscid_Gaussian_Bump) then
      return_value = (gppt(npres)/pref_nd) / (upt(nec)**gamm)
    else
      return_value = log(gppt(npres)) - gamm*log(upt(nec))
    end if
    !
  else if (plot_variable == plot_vorticity_x) then
    !
    return_value = vorticity_dusp( upt , gppt , dupt , 1 )   / time_ref
    !
  else if (plot_variable == plot_vorticity_y) then
    !
    return_value = vorticity_dusp( upt , gppt , dupt , 2 )   / time_ref
    !
  else if (plot_variable == plot_vorticity_z) then
    !
    return_value = vorticity_dusp( upt , gppt , dupt , 3 )   / time_ref
    !
  else if (plot_variable == plot_grad_density(1)) then
    !
    return_value = dupt(1,nec)                        * rhoref
    !
  else if (plot_variable == plot_grad_density(2)) then
    !
    return_value = dupt(2,nec)                        * rhoref
    !
  else if (plot_variable == plot_grad_density(3)) then
    !
    return_value = dupt(3,nec)                        * rhoref
    !
  else if (plot_variable == plot_grad_momentum_x(1)) then
    !
    return_value = dupt(1,nmx)                        * rhoref * aref
    !
  else if (plot_variable == plot_grad_momentum_x(2)) then
    !
    return_value = dupt(2,nmx)                        * rhoref * aref
    !
  else if (plot_variable == plot_grad_momentum_x(3)) then
    !
    return_value = dupt(3,nmx)                        * rhoref * aref
    !
  else if (plot_variable == plot_grad_momentum_y(1)) then
    !
    return_value = dupt(1,nmy)                        * rhoref * aref
    !
  else if (plot_variable == plot_grad_momentum_y(2)) then
    !
    return_value = dupt(2,nmy)                        * rhoref * aref
    !
  else if (plot_variable == plot_grad_momentum_y(3)) then
    !
    return_value = dupt(3,nmy)                        * rhoref * aref
    !
  else if (plot_variable == plot_grad_momentum_z(1)) then
    !
    return_value = dupt(1,nmz)                        * rhoref * aref
    !
  else if (plot_variable == plot_grad_momentum_z(2)) then
    !
    return_value = dupt(2,nmz)                        * rhoref * aref
    !
  else if (plot_variable == plot_grad_momentum_z(3)) then
    !
    return_value = dupt(3,nmz)                        * rhoref * aref
    !
  else if (plot_variable == plot_grad_energy(1)) then
    !
    return_value = dupt(1,nee)                        * rhoref * aref * aref
    !
  else if (plot_variable == plot_grad_energy(2)) then
    !
    return_value = dupt(2,nee)                        * rhoref * aref * aref
    !
  else if (plot_variable == plot_grad_energy(3)) then
    !
    return_value = dupt(3,nee)                        * rhoref * aref * aref
    !
  end if
  !
end function compute_output_variable
!
!###############################################################################
!
elemental subroutine get_function_name(this)
  !
  !.. Formal Arguments ..
  class(plot_func_t), intent(inout) :: this
  !
continue
  !
  if (this == plot_static_density) then
    this%varname = "Density"
  else if (this == plot_velocity_x) then
    this%varname = "VelocityX"
  else if (this == plot_velocity_y) then
    this%varname = "VelocityY"
  else if (this == plot_velocity_z) then
    this%varname = "VelocityZ"
  else if (this == plot_momentum_x) then
    this%varname = "MomentumX"
  else if (this == plot_momentum_y) then
    this%varname = "MomentumY"
  else if (this == plot_momentum_z) then
    this%varname = "MomentumZ"
  else if (this == plot_static_pressure) then
    this%varname = "Pressure"
  else if (this == plot_static_temperature) then
    this%varname = "Temperature"
  else if (this == plot_total_energy) then
    this%varname = "EnergyStagnationDensity"
  else if (this == plot_mach_number) then
    this%varname = "Mach"
  else if (this == plot_entropy) then
    this%varname = "Entropy"
  else if (this == plot_vorticity_x) then
    this%varname = "VorticityX"
  else if (this == plot_vorticity_y) then
    this%varname = "VorticityY"
  else if (this == plot_vorticity_z) then
    this%varname = "VorticityZ"
  else if (this == plot_grad_density(1)) then
    this%varname = "DensityGradientX"
  else if (this == plot_grad_density(2)) then
    this%varname = "DensityGradientY"
  else if (this == plot_grad_density(3)) then
    this%varname = "DensityGradientZ"
  else if (this == plot_grad_momentum_x(1)) then
    this%varname = "MomentumXGradientX"
  else if (this == plot_grad_momentum_x(2)) then
    this%varname = "MomentumXGradientY"
  else if (this == plot_grad_momentum_x(3)) then
    this%varname = "MomentumXGradientZ"
  else if (this == plot_grad_momentum_y(1)) then
    this%varname = "MomentumYGradientX"
  else if (this == plot_grad_momentum_y(2)) then
    this%varname = "MomentumYGradientY"
  else if (this == plot_grad_momentum_y(3)) then
    this%varname = "MomentumYGradientZ"
  else if (this == plot_grad_momentum_z(1)) then
    this%varname = "MomentumZGradientX"
  else if (this == plot_grad_momentum_z(2)) then
    this%varname = "MomentumZGradientY"
  else if (this == plot_grad_momentum_z(3)) then
    this%varname = "MomentumZGradientZ"
  else if (this == plot_grad_energy(1)) then
    this%varname = "EnergyGradientX"
  else if (this == plot_grad_energy(2)) then
    this%varname = "EnergyGradientY"
  else if (this == plot_grad_energy(3)) then
    this%varname = "EnergyGradientZ"
  end if
  !
  this%varname = trim(adjustl(this%varname))
  !
end subroutine get_function_name
!
!###############################################################################
!
subroutine get_plot_functions
  !
  !.. Local Scalars ..
  integer :: ierr,npvar
  integer, allocatable :: isvar_s(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_plot_functions"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  npvar = 20
  !
  allocate ( plot_var(1:npvar) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"plot_var",1,__LINE__,__FILE__,ierr,error_message)
  !
  plot_var = [ plot_static_density, &
               plot_momentum_x, &
               plot_momentum_y, &
               plot_momentum_z, &
               plot_total_energy, &
               plot_grad_density, &
               plot_grad_momentum_x, &
               plot_grad_momentum_y, &
               plot_grad_momentum_z, &
               plot_grad_energy ]
  !
  ! Get the name for each of the functions that will be plotted
  !
  call plot_var%get_function_name
  !
  ! Allocate the index array for the CGNS field solutions
  !
  allocate ( isvar_s(1:npvar) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"isvar_s",1,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_plot_functions
!
!###############################################################################
!
subroutine map_to_output_points(this_geom,this_order,map)
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: this_geom
  integer,                 intent(in) :: this_order
  type(geom_pnt_map_t), intent(inout) :: map
  !
  !.. Local Scalars ..
  integer :: ne,nf,nep,nfp,nnp,nsp,ierr
  integer :: i,j,k,l,m,n,np1d,np2d,nt1d,nt2d
  !
  character(len=100) :: map_fmt
  !
  !.. Local Arrays ..
  integer, dimension(1:2,1:this_order+1,1:12) :: edge_values
  integer, dimension(1:2,1:(this_order+1)**2,1:6) :: face_values
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "map_to_output_points"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  !
  !###################################
  !###################################
  !#####          FACES          #####
  !###################################
  !###################################
  !
  ! Number of solution points for this combination
  ! of cell geometry and solution order
  !
  nsp = Cell_Solpts(this_geom,this_order)
  !
  ! Number of faces for this cell
  !
  nfp = 0
  do nf = 1,geom_faces(this_geom)
    nfp = nfp + Cell_Solpts( cell_face_geom(nf,this_geom) , this_order )
  end do
  !
  ! Number of edges for this cell
  !
  nep = 0
  if (any(this_geom == Geom_3D)) then
    do ne = 1,geom_edges(this_geom)
      nep = nep + Cell_Solpts(Geom_Edge,this_order)
    end do
  end if
  !
  ! Number of nodes for this cell
  !
  nnp = geom_nodes(this_geom)
  !
  ! Allocate the different components of the map derived type
  !
  allocate ( map%sp(1:nsp) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"map%sp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( map%fp(1:nfp) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"map%sp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( map%ep(1:nep) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"map%sp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( map%np(1:nnp) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"map%sp",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Compute the maps from the local indices to the output indices
  !
  select case (this_geom)
      !
    case (Geom_Tria)
      !
      map%sp = get_tria_solpt_indices(this_order)
      map%fp = get_tria_face_indices(this_order)
      map%np = get_tria_node_indices(this_order)
      !
    case (Geom_Quad)
      !
      map%sp = get_quad_solpt_indices(this_order)
      map%fp = get_quad_face_indices(this_order)
      map%np = get_quad_node_indices(this_order)
      !
    case (Geom_Tetr)
      !
      map%sp = get_tetr_solpt_indices(this_order)
      map%fp = get_tetr_face_indices(this_order)
      map%ep = get_tetr_edge_indices(this_order)
      map%np = get_tetr_node_indices(this_order)
      !
    case (Geom_Pyra)
      !
      map%sp = get_pyra_solpt_indices(this_order)
      map%fp = get_pyra_face_indices(this_order)
      map%ep = get_pyra_edge_indices(this_order)
      map%np = get_pyra_node_indices(this_order)
      !
    case (Geom_Pris)
      !
      map%sp = get_pris_solpt_indices(this_order)
      map%fp = get_pris_face_indices(this_order)
      map%ep = get_pris_edge_indices(this_order)
      map%np = get_pris_node_indices(this_order)
      !
    case (Geom_Hexa)
      !
      map%sp = get_hexa_solpt_indices(this_order)
      map%fp = get_hexa_face_indices(this_order)
      map%ep = get_hexa_edge_indices(this_order)
      map%np = get_hexa_node_indices(this_order)
      !
    case default
      !
      write (error_message,1)
      call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
      !
  end select
  !
  ! Finally, adjust the mapping indices so that they start at 0 instead of 1
  !
  map%sp = map%sp - 1
  map%fp = map%fp - 1
  map%ep = map%ep - 1
  map%np = map%np - 1
  !
  np1d = this_order+1
  np2d = np1d*np1d
  nt1d = np1d + 2
  nt2d = nt1d*nt1d
  !
  edge_values = 0
  face_values = 0
  !
  l = 0
  do n = 1,size(face_values,dim=3)
    do m = 1,size(face_values,dim=2)
      l = l + 1
      face_values(:,m,n) = [l,map%fp(l)+1]
    end do
  end do
  !
!#define DEBUG_MAP_TO_OUTPUT_POINTS
#ifdef DEBUG_MAP_TO_OUTPUT_POINTS
  do n = 1,2
    write (iout,14) (j,j=3*n-2,3*n)
    do m = 1,np2d
      write (iout,15) ((face_values(i,m,j), i=1,2), j=3*n-2,3*n)
    end do
  end do
  !
  write (iout,21)
  do l = 1,size(map%np)
    write (iout,22) l,map%np(l)+1
  end do
#endif
  !
  l = 0
  do n = 1,size(edge_values,dim=3)
    do m = 1,size(edge_values,dim=2)
      l = l + 1
      edge_values(:,m,n) = [l,map%ep(l)+1]
    end do
  end do
  !
#ifdef DEBUG_MAP_TO_OUTPUT_POINTS
  do n = 1,4
    write (iout,35) (j,j=3*n-2,3*n)
    do m = 1,np1d
      write (iout,36) ((edge_values(i,m,j), i=1,2), j=3*n-2,3*n)
    end do
  end do
  !
  write (map_fmt,43) np1d,np1d,np1d
  write (iout,41)
  l = 0
  do k = 2,np1d+1
    write (iout,42) k-1
    do j = 2,np1d+1
      write (iout,map_fmt) j-1, (l+i, i=1,np1d), &
                           ((k-1)*nt2d+(j-1)*nt1d+i, i=2,np1d+1), &
                           (map%sp(l+i)+1, i=1,np1d)
      l = l + np1d
    end do
  end do
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("ERROR! Encountered invalid geometry when trying to compute ", &
            "the arrays to map local point indices to the corresponding ", &
            "output indices.")
  14 format (/,3(2x,"Face",i2,:,17x))
  15 format (3(4x,"fp",i4," =>",i4,:,4x))
  21 format (//,"  Nodes")
  22 format ("    np ",i2," =>",2(1x,i3),";  [i,j,k] = [",i0,",",i0,",",i0,"]")
  35 format (/,3(2x,"Edge",i2,:,17x))
  36 format (3(4x,"ep",i4," =>",i4,:,4x))
  41 format (//,"  Solution Points")
  42 format (" SP  k = ",i0)
  43 format ("(7x,'j = ',i0,3x,",i0,"i3,3x,",i0,"i4,3x,",i0,"i4)")
  !
end subroutine map_to_output_points
!
!###############################################################################
!
pure function get_hexa_solpt_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, dimension(1:(order+1)**3) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,n,np1d,nt1d
  !
continue
  !
  return_value = 0
  !
  np1d = order + 1
  nt1d = np1d + 2
  !
  l = 0
  n = nt1d*nt1d
  !
  do k = 1,np1d
    n = n + nt1d ! skip all points along the j-min face
    do j = 1,np1d
      n = n + 1 ! skip the point for the i-min face
      do i = 1,np1d
        l = l + 1
        n = n + 1
        return_value(l) = n
      end do
      n = n + 1 ! skip the point for the i-max face
    end do
    n = n + nt1d ! skip all points along the j-max face
  end do
  !
end function get_hexa_solpt_indices
!
!###############################################################################
!
pure function get_hexa_node_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, dimension(1:8) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,ip,jp,kp,np1d,nt1d,nt2d
  !
continue
  !
  return_value = 0
  !
  np1d = order + 1
  nt1d = np1d + 2
  nt2d = nt1d**2
  !
  l = 0
  do k = 1,2
    kp = (k-1)*(nt1d-1)
    do j = 1,2
      jp = (j-1)*(nt1d-1)
      do i = 1,2
        ip = (i-1)*(nt1d-1) - (2*i-3)*jp
        l = l + 1
        return_value(l) = kp*nt2d + jp*nt1d + ip + 1
      end do
    end do
  end do
 !return_value(1) = 1
 !return_value(2) = nt1d
 !return_value(3) = nt1d*nt1d
 !return_value(4) = nt1d*nt1d - nt1d + 1
 !return_value(5) = nt1d*nt1d*(nt1d-1) + 1
 !return_value(6) = nt1d*nt1d*(nt1d-1) + nt1d
 !return_value(7) = nt1d*nt1d*(nt1d-1) + nt1d*nt1d
 !return_value(8) = nt1d*nt1d*(nt1d-1) + nt1d*nt1d - nt1d + 1

  !
end function get_hexa_node_indices
!
!###############################################################################
!
pure function get_hexa_face_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, dimension(1:6*(order+1)*(order+1)) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,ip,iexp,istride
  integer :: j,jp,jexp,jstride
  integer :: k,kp,kexp,kstride
  integer :: l,l1,l2,nf,np1d,nt1d
  !
  !.. Local Parameters ..
  integer, parameter :: dx = 0 ! The i/x/xi   direction
  integer, parameter :: dy = 1 ! The j/y/eta  direction
  integer, parameter :: dz = 2 ! The k/z/zeta direction
  !
  ! fexp(1,:) :  gives the stride exponent for the direction in which the local
  !              face indices are consecutive
  !
  ! fexp(2,:) :  gives the stride exponent for the direction in which the local
  !              face indices are not consecutive
  !
  integer, parameter :: fexp(2,6) = reshape( [dy,dx, dx,dy, &
                                              dy,dz, dz,dy, &
                                              dx,dz, dz,dx], shape(fexp) )
  !
continue
  !
  return_value = 0
  !
  np1d = order + 1
  nt1d = np1d + 2
  !
  l = 0
  do nf = 1,6
    !
    l1 = l + 1
    !
    k = merge( 1 , nt1d , any(nf == [Xmin3D_Face,Ymin3D_Face,ZMin3D_Face]) )
    !
    iexp = fexp(1,nf)
    jexp = fexp(2,nf)
    kexp = 3 - iexp - jexp
    !
    istride = nt1d**iexp
    jstride = nt1d**jexp
    kstride = nt1d**kexp
    !
    kp = (k-1)*kstride
    do j = 2,np1d+1
      jp = (j-1)*jstride
      do i = 2,np1d+1
        ip = (i-1)*istride
        l = l + 1
        return_value(l) = kp + jp + ip + 1
      end do
    end do
    !
    l2 = l
    !
    if (any(nf == [Xmin3D_Face,Ymax3D_Face,Zmin3D_Face])) then
      return_value(l1:l2) = transpose_face(return_value(l1:l2))
    end if
    !
  end do
  !
end function get_hexa_face_indices
!
!###############################################################################
!
pure function get_hexa_edge_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result ..
  integer, dimension(1:12*(order+1)) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,l,ne,idir,np1d,nt1d
  integer :: ic,ip,iexp,istride
  integer :: jc,jp,jexp,jstride
  integer :: kc,kp,kexp,kstride
  !
  !.. Local Parameters ..
  integer, parameter :: dx = 0 ! The i/x/xi   direction
  integer, parameter :: dy = 1 ! The j/y/eta  direction
  integer, parameter :: dz = 2 ! The k/z/zeta direction
  !
  integer, parameter :: fwd = 0 ! The edge-pt indices are in increasing order
  integer, parameter :: rev = 1 ! The edge-pt indices are in decreasing order
  !
  integer, parameter :: dmin = 0 ! The direction is a min face
  integer, parameter :: dmax = 1 ! The direction is a max face
  !
  ! eexp(1,:) : gives the stride exponent for the direction of
  !             lowest  dimension that is constant for the edge
  !
  ! eexp(2,:) : gives the stride exponent for the direction of
  !             highest dimension that is constant for the edge
  !
  integer, parameter :: eexp(2,12) = reshape( [dy,dz, dx,dz, dy,dz, dx,dz,  &
                                               dx,dy, dx,dy, dx,dy, dx,dy,  &
                                               dy,dz, dx,dz, dy,dz, dx,dz], &
                                               shape(eexp) )
  !
  ! ecoef(1,:) : gives the coefficient that puts the edge-pt indices in
  !              increasing or decreasing order for the non-constant
  !              direction of the edge
  !
  ! ecoef(2,:) : gives the coefficient that adjusts the edge-pt
  !              indices based on whether the edge is on a
  !              min or max face for the direction eexp(1,:)
  !
  ! ecoef(3,:) : gives the coefficient that adjusts the edge-pt
  !              indices based on whether the edge is on a
  !              min or max face for the direction eexp(2,:)
  !
  integer, parameter :: ecoef(3,12) = reshape( [fwd,dmin,dmin, fwd,dmax,dmin, &
                                                rev,dmax,dmin, rev,dmin,dmin, &
                                                fwd,dmin,dmin, fwd,dmax,dmin, &
                                                fwd,dmax,dmax, fwd,dmin,dmax, &
                                                fwd,dmin,dmax, fwd,dmax,dmax, &
                                                rev,dmax,dmax, rev,dmin,dmax], &
                                                shape(ecoef) )
  !
continue
  !
  return_value = 0
  !
  np1d = order + 1
  nt1d = np1d + 2
  !
  l = 0
  do ne = 1,12
    !
    ic = ecoef(1,ne)
    jc = ecoef(2,ne)
    kc = ecoef(3,ne)
    !
    jexp = eexp(1,ne)
    kexp = eexp(2,ne)
    iexp = 3 - kexp - jexp
    !
    istride = nt1d**iexp
    jstride = nt1d**jexp
    kstride = nt1d**kexp
    !
    ip = ic * (nt1d-1) * istride
    jp = jc * (nt1d-1) * jstride
    kp = kc * (nt1d-1) * kstride
    !
    do i = 2,np1d+1
      l = l + 1
      idir = (2*ic-1) * (1-i) * istride
      return_value(l) = kp + jp + ip + idir + 1
    end do
    !
  end do
  !
end function get_hexa_edge_indices
!
!###############################################################################
!
pure function get_tria_solpt_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Tria,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_tria_solpt_indices
!
!###############################################################################
!
pure function get_tria_face_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Tria,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_tria_face_indices
!
!###############################################################################
!
pure function get_tria_node_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Tria,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_tria_node_indices
!
!###############################################################################
!
pure function get_quad_solpt_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Quad,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_quad_solpt_indices
!
!###############################################################################
!
pure function get_quad_face_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Quad,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_quad_face_indices
!
!###############################################################################
!
pure function get_quad_node_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Quad,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_quad_node_indices
!
!###############################################################################
!
pure function get_tetr_solpt_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Tetr,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_tetr_solpt_indices
!
!###############################################################################
!
pure function get_tetr_face_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Tetr,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_tetr_face_indices
!
!###############################################################################
!
pure function get_tetr_edge_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Tetr,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_tetr_edge_indices
!
!###############################################################################
!
pure function get_tetr_node_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Tetr,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_tetr_node_indices
!
!###############################################################################
!
pure function get_pyra_solpt_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Pyra,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_pyra_solpt_indices
!
!###############################################################################
!
pure function get_pyra_face_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Pyra,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_pyra_face_indices
!
!###############################################################################
!
pure function get_pyra_edge_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Pyra,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_pyra_edge_indices
!
!###############################################################################
!
pure function get_pyra_node_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Pyra,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_pyra_node_indices
!
!###############################################################################
!
pure function get_pris_solpt_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Pris,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_pris_solpt_indices
!
!###############################################################################
!
pure function get_pris_face_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Pris,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_pris_face_indices
!
!###############################################################################
!
pure function get_pris_edge_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Pris,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_pris_edge_indices
!
!###############################################################################
!
pure function get_pris_node_indices(order) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: order
  !
  !.. Function Result .
  integer, allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: np,ierr
  !
continue
  !
  np = Cell_Solpts(Geom_Pris,order)
  !
  allocate ( return_value(1:np) , source=0 , stat=ierr )
  !
end function get_pris_node_indices
!
!###############################################################################
!
subroutine get_output_connectivity(bface,cell_geom, &
                                   nodes_of_cell,nodes_of_cell_ptr)
  !
  !.. Use Statements ..
  use connectivity_mod, only : find_cells_surrounding_each_node
  use connectivity_mod, only : find_cells_surrounding_each_cell
  use connectivity_mod, only : find_nodes_on_each_edge
  use connectivity_mod, only : find_cells_surrounding_each_edge
  use connectivity_mod, only : get_unit_cell_subconnectivity
  !
  !.. Formal Arguments ..
  integer, dimension(:,:), intent(in) :: bface
  integer, dimension(:),   intent(in) :: cell_geom
  integer, dimension(:),   intent(in) :: nodes_of_cell
  integer, dimension(:),   intent(in) :: nodes_of_cell_ptr
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:)   :: cells_surr_cell
  integer, allocatable, dimension(:)   :: cells_surr_cell_ptr
  integer, allocatable, dimension(:,:) :: nodes_on_edge
  !
  integer, allocatable, dimension(:)   :: cells_with_node
  integer, allocatable, dimension(:)   :: cells_with_node_ptr
  integer, allocatable, dimension(:)   :: cells_with_edge
  integer, allocatable, dimension(:)   :: cells_with_edge_ptr
  integer, allocatable, dimension(:,:) :: cells_on_face
  !
  !.. Local Scalars ..
  integer :: lnode,lfbnd,lcell
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_output_connectivity"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Set the module size parameters
  !
  lnode = maxval(nodes_of_cell)
  lfbnd = size(bface,dim=2)
  lcell = size(cell_geom) - lfbnd
  !
  ! Get the cells surrounding each node
  !
  call find_cells_surrounding_each_node(bface, &
                                        nodes_of_cell,nodes_of_cell_ptr, &
                                        cells_with_node,cells_with_node_ptr)
  !
  ! Get the cells surrounding each cell
  !
  call find_cells_surrounding_each_cell(bface,cell_geom, &
                                        nodes_of_cell,nodes_of_cell_ptr, &
                                        cells_with_node,cells_with_node_ptr, &
                                        cells_surr_cell,cells_surr_cell_ptr)
  !
  ! Get the nodes defining each edge as well as the edges surrounding each node
  !
  call find_nodes_on_each_edge(cell_geom, &
                               nodes_of_cell,nodes_of_cell_ptr, &
                               cells_with_node,cells_with_node_ptr, &
                               nodes_on_edge)
  !
  ! Get the cells that surround each edge
  !
  call find_cells_surrounding_each_edge(nodes_on_edge, &
                                        nodes_of_cell,nodes_of_cell_ptr, &
                                        cells_with_node,cells_with_node_ptr, &
                                        cells_with_edge,cells_with_edge_ptr)
  !
  ! Get the cells on each face
  !
  call find_cells_on_each_face(bface,cells_surr_cell,cells_surr_cell_ptr, &
                               cells_on_face)
  !
  ! Get the unit cell subconnectivities for use when writing solution file
  !
  call get_unit_cell_subconnectivity
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_output_connectivity
!
!###############################################################################
!
subroutine find_cells_on_each_face(bface,cells_surr_cell,cells_surr_cell_ptr, &
                                   cells_on_face)
  !
  ! Procedure that creates the face arrays. Look in
  ! the module geovar for definitions of each array.
  !
  !.. Formal Arguments ..
  integer,              dimension(:,:),    intent(in) :: bface
  integer, allocatable, dimension(:),   intent(inout) :: cells_surr_cell
  integer, allocatable, dimension(:),   intent(inout) :: cells_surr_cell_ptr
  integer, allocatable, dimension(:,:), intent(inout) :: cells_on_face
  !
  !.. Local Scalars ..
  integer :: c1,c2,n1,n2,nf,nc,ierr
  integer :: face_counter
  integer :: lfbnd,lface,lcell
  integer :: this_cell,adj_cell
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_cells_on_each_face"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  lfbnd = size(bface,dim=2)
  lcell = size(cells_surr_cell_ptr)-1
  !
  ! ###########
  ! ###########
  ! ###  1  ###
  ! ###########
  ! ###########
  !
  ! Find the total number of faces
  ! and then allocate the face array
  !
  lface = lfbnd
  do nc = 1,lcell
    !
    n1 = cells_surr_cell_ptr(nc)+1
    n2 = cells_surr_cell_ptr(nc+1)
    !
    lface = lface + count( cells_surr_cell(n1:n2)>nc .and. &
                           cells_surr_cell(n1:n2)<=lcell )
    !
  end do
  !
  allocate ( cells_on_face(1:2,1:lface) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_on_face",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! ###########
  ! ###########
  ! ###  2  ###
  ! ###########
  ! ###########
  !
  ! Get the cells on each boundary face
  !
  do nf = 1,lfbnd
    !
    cells_on_face(1,nf) = bface( 2,nf) ! Host  cell on boundary face
    cells_on_face(2,nf) = bface(10,nf) ! Ghost cell on boundary face
    !
  end do
  !
  ! ###########
  ! ###########
  ! ###  3  ###
  ! ###########
  ! ###########
  !
  ! Get the cells for the interior faces
  !
  face_counter = lfbnd
  !
  do this_cell = 1,lcell
    !
    c1 = cells_surr_cell_ptr(this_cell)+1
    c2 = cells_surr_cell_ptr(this_cell+1)
    !
    do nc = c1,c2
      !
      adj_cell = cells_surr_cell(nc)
      !
      if (adj_cell>this_cell .and. adj_cell<=lcell) then
        !
        face_counter = face_counter + 1
        !
        ! Store the cells on this interior face
        !
        cells_on_face(1,face_counter) = this_cell
        cells_on_face(2,face_counter) = adj_cell
        !
      end if
      !
    end do
    !
  end do
  !
  ! We are finished with cells_surr_cell
  ! and cells_surr_cell_ptr so deallocate
  !
  deallocate ( cells_surr_cell , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_surr_cell",2,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  deallocate ( cells_surr_cell_ptr , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_surr_cell_ptr",2,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine find_cells_on_each_face
!
!###############################################################################
!
pure function vorticity_dusp(cv,gp,dcvdx,dir) result(return_value)
  !
  use gpu_mod, only : grad_cv_to_grad_pv_gpu
  !
  integer,  intent(in) :: dir
  real(wp), intent(in) :: cv(:)
  real(wp), intent(in) :: gp(:)
  real(wp), intent(in) :: dcvdx(:,:)
  !
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,n1,n2
  !
  !.. Local Arrays ..
  real(wp) :: dpvdx(1:size(dcvdx,dim=1),1:size(dcvdx,dim=2))
  !
  !.. Local Parameters ..
  integer, parameter, dimension(3) :: idx = [2,3,1]
  integer, parameter, dimension(3) :: jdx = [3,1,2]
  !
continue
  !
  n1 = size(dcvdx,dim=1)
  n2 = size(dcvdx,dim=2)
  !
  i = idx(dir)
  j = jdx(dir)
  !
  call grad_cv_to_grad_pv_gpu( cv , gp , dcvdx , dpvdx )
  !
  return_value = dpvdx(i,nmb-1+j) - dpvdx(j,nmb-1+i)
  !
end function vorticity_dusp
!
!###############################################################################
!
subroutine check_cgns_dir
  !
  !.. Use Statements ..
  use ovar, only : output_dir,itrst,num_timesteps
  use ovar, only : lustre_stripe_count
  use ovar, only : lustre_stripe_size
  !
  !.. Local Scalars ..
  integer :: clen,ierr,iodir
  integer :: exitstat,cmndstat
  integer :: final_iter_digits
  logical(lk) :: using_cgns_dir
  character(len=300) :: command
  character(len=300) :: cmndmsg
  character(len=200) :: test_file
  character(len=200) :: test_dir
  character(len=100) :: char_num
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_cgns_dir"
  !
  logical(ldk), parameter :: wait_for_completion = .true.
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize exitstat and cmndmsg
  !
  exitstat = 0
  cmndmsg = ""
  !
  ! Try to open a file in the directory named CGNS_files
  !
  test_dir = trim(output_dir) // "CGNS_files/"
  test_file = trim(test_dir) // "cgns_dir.test_file.dat"
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
      ! Set the logical saying that we are using a separate CGNS directory
      !
      using_cgns_dir = true
      cgns_dir = trim(test_dir)
      !
    else
      !
      ! The open failed, so more than likely the directory CGNS_files
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
        using_cgns_dir = fals
        cgns_dir = trim(output_dir)
      else
        using_cgns_dir = true
        cgns_dir = trim(test_dir)
      end if
      !
    end if
    !
#ifdef PBS_ENV
    !
    ! If we are using a separate CGNS directory and we are running
    ! this through PBS on NAS, try to modify the lfs striping for
    ! the CGNS directory for better I/O performance.
    !
    if (using_cgns_dir) then
      !
     !command = "lfs setstripe -c 64 -S 4m "//trim(cgns_dir)
      write (command,1) max(0,lustre_stripe_count), &
                        trim(adjustl(lustre_stripe_size)), &
                        trim(adjustl(cgns_dir))
      write (iout,2) trim(adjustl(command))
      !
      call execute_command(command,wait_for_completion,exitstat, &
                           cmndstat,cmndmsg,pname)
     !call execute_command_line( trim(command) , wait_for_completion , &
     !                           exitstat , cmndstat , cmndmsg )
     !call write_command_results(pname,command,exitstat,cmndstat,cmndmsg)
      !
      ! We arent going to test the status of executing this command because
      ! it worked or it didnt. Either way, we will continue on and output
      ! the CGNS solution to the CGNS directory.
      !
    end if
    !
#endif
    !
  end if
  !
  ! Broadcast the value of cgns_dir to all the non-root MPI processes
  !
  clen = len_trim(cgns_dir)
#ifdef VERBOSE_MPI
  write (iout,dbg_mpi_frmt) mypnum,"before mpi_bcast",__LINE__,__FILE__
#endif
  call mpi_bcast(clen,1_int_mpi,mpi_inttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
  write (iout,dbg_mpi_frmt) mypnum," after mpi_bcast",__LINE__,__FILE__
#endif
#ifdef VERBOSE_MPI
  write (iout,dbg_mpi_frmt) mypnum,"before mpi_bcast",__LINE__,__FILE__
#endif
  call mpi_bcast(cgns_dir(1:clen),int(clen,kind=int_mpi),MPI_CHARACTER, &
                 0_int_mpi,MPI_COMM_WORLD,mpierr)
#ifdef VERBOSE_MPI
  write (iout,dbg_mpi_frmt) mypnum," after mpi_bcast",__LINE__,__FILE__
#endif
  !
  ! Create the format string for creating the CGNS file names
  !
  write (char_num,3) itrst + num_timesteps
  final_iter_digits = max(4,len_trim(adjustl(char_num)))
  !
  write (cgnsfile_fmt,4) final_iter_digits
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("lfs setstripe -c ",i0," -S ",a,1x,a)
  2 format ("Lustre striping command = '",a,"'")
  3 format (i0)
  4 format ("('solution.',i0.",i0,",a,'.cgns')")
  !
end subroutine check_cgns_dir
!
!###############################################################################
!
subroutine write_tecplot_file(xyz_pts,ipelem)
  !
  !.. Formal Arguments ..
  real(wp),     dimension(:,:), intent(in) :: xyz_pts
  integer(CST), dimension(:,:), intent(in) :: ipelem
  !
  !.. Local Scalars ..
  integer :: l,n,nc,iotec,ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_tecplot_file"
  character(len=*), parameter :: fname = "postproc.tec.dat"
  character(len=*), parameter :: cname = "postproc.cells.tec.dat"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  open (newunit=iotec,file=fname,form="formatted",status="replace", &
                      action="write",position="rewind",iostat=ierr, &
                      iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (iotec,1) size(xyz_pts,dim=2),size(ipelem,dim=2)
  !
  do l = 1,size(xyz_pts,dim=1)
    write (iotec,2) (xyz_pts(l,n), n=1,size(xyz_pts,dim=2))
  end do
  !
  do n = 1,size(ipelem,dim=2)
    write (iotec,3) (ipelem(l,n), l=1,size(ipelem,dim=1))
  end do
  !
  if (allocated(ipelem_of_cells)) then
    !
    do nc = 1,size(ipelem_of_cells,dim=3)
      !
      write (iotec,4) nc,size(xyz_pts,dim=2),size(ipelem_of_cells,dim=2)
      !
      do n = 1,size(ipelem_of_cells,dim=2)
        write (iotec,3) (ipelem_of_cells(l,n,nc), &
                                l=1,size(ipelem_of_cells,dim=1))
      end do
      !
    end do
    !
    deallocate ( ipelem_of_cells , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ipelem_of_cells",2,__LINE__,__FILE__,ierr, &
                     error_message)
    !
  end if
  !
  close (iotec,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ('VARIABLES = "X", "Y", "Z"',/, &
            'ZONE T="Continuous Output", N=',i0,', E=',i0, &
            ', DATAPACKING=BLOCK, ZONETYPE=FEBRICK')
  2 format (5e18.10)
  3 format (8(1x,i0))
  4 format ('ZONE T="Grid Cell ',i0,'", N=',i0,', E=',i0, &
            ', DATAPACKING=BLOCK, ZONETYPE=FEBRICK',/, &
            'VARSHARELIST=([1-3]=1)')
  !
end subroutine write_tecplot_file
!
!###############################################################################
!
subroutine write_cell_points(iunit,nep,cell_points_ptr,cell_points)
  !
  !.. Formal Arguments ..
  integer,            intent(inout) :: iunit
  integer,               intent(in) :: nep
  integer, dimension(:), intent(in) :: cell_points_ptr
  integer, dimension(:), intent(in) :: cell_points
  !
  !.. Local Scalars ..
  integer :: i,n,n1,n2
  !
  character(len=300) :: map_fmt
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_cell_points"
  !
continue
  !
#ifdef DEBUG_POST_PROCESSOR
  call debug_timer(entering_procedure,pname)
  !
  write (map_fmt,1) (nep,i=1,6)
  !
  do n = 1,size(cell_points_ptr)-1
    n1 = cell_points_ptr(n)+1
    n2 = cell_points_ptr(n+1)
   !write (iunit,1) n,(cell_points(i),i=n1,n2)
    write (iunit,map_fmt) n,(cell_points(i),i=n1,n2)
  end do
  !
  iunit = iunit + 1
  !
  call debug_timer(leaving_procedure,pname)
#endif
  !
  ! Format Statements
  !
  1 format ("('  CELL ',i0,/,60('-'),/,", &
            "'Solution Points',/,", &
            i0,"(",i0,"(",i0,"i5,/),/),/,", &
            "'Face Points',/,", &
              "6(",i0,"(",i0,"i5,/),/),/,", &
            "'Edge Points',/,", &
                    "12(",i0,"i5,/),/,", &
            "'Node Points',/,", &
                           "8i5,/,", &
            "60('-'),/)")
 !2 format ("  CELL ",i0,/,60("-"),/, &
 !          "Solution Points",/, &
 !          <nep>(<nep>(<nep>i5,/),/),/, &
 !          "Face Points",/, &
 !              6(<nep>(<nep>i5,/),/),/, &
 !          "Edge Points",/, &
 !                   12(<nep>i5,/),/, &
 !          "Node Points",/, &
 !                          8i5,/, &
 !          60("-"),/)
  !
end subroutine write_cell_points
!
!###############################################################################
!
include "Functions/last.f90"
include "Functions/cell_solpts.f90"
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module postprocessor_mod
