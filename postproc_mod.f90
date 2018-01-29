module postprocessor_mod
  !
  !.. Global Use Statements ..
  use module_kind_types
  use module_cgns_types
  use eqn_idx, only : nq,nec,nmx,nmy,nmz,nee,nmb,nme,ntk
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
  integer, save :: npvar = 20
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
  character(len=*),      parameter :: cgns_base_fname = "solution.grid.cgns"
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
  integer(CBT), save, allocatable, dimension(:) :: isvar_s
  !
  !.. CGNS Time Dependent Variables ..
  integer, save :: cgns_zone_counter
  !
  real(sp),             save, allocatable, dimension(:) :: soltime
  character(len=CGLEN), save, allocatable, dimension(:) :: solname
  !
  character(len=*), parameter :: base_name = "gfr_base"
  character(len=*), parameter :: zone_name = "Solution"
  !
  !.. CGNS Parallel Point Indices ..
  integer, save :: n_my_output_pts
  integer, save :: my_pts_beg
  integer, save :: my_pts_end
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
subroutine new_write_parallel_cgns()
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : itcur,time,time_ref
  use ovar,   only : completely_disable_cgns
  use ovar,   only : write_TaylorGreen_full_solution
  !
  use, intrinsic :: iso_c_binding, only : C_NULL_CHAR
  !
  !.. Local Scalars ..
  integer :: n,ierr
  !
  character(len=CGLEN) :: tname
  character(len=170)   :: sol_file
  character(len=300)   :: link_path
  !
  !.. Local Arrays ..
  integer, dimension(1:2) :: idata
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
    ! Initialize the functions that will be written to the CGNS solution files
    !
    call get_plot_functions
    !
    ! Write the CGNS grid file
    !
    call write_grid_parallel_cgns
    !
    ! Now that the CGNS grid file has been successfully created, set the
    ! initialization flag to false so this isnt done again.
    !
    initialization_needed = fals
    !
  else
    !
    ! The CGNS grid file has already been created so we need to reopen
    ! it in modification mode.
    !
    call cgp_open_f(cgns_grid_file,CG_MODE_MODIFY,ifile_n,cgierr)
    call cgns_error(pname,cgns_grid_file,"cgp_open_f", &
                    cgierr,__LINE__,__FILE__,"Mode = CG_MODE_MODIFY")
    !
  end if
  !
  ! Get the name for the current solution node
  !
  write (tname,'("FlowSolution",i0)') cgns_zone_counter
 !write (sol_file,'("solution.",i0,".cgns")') cgns_zone_counter
  write (sol_file,cgnsfile_fmt) itcur
  !
  ! Reallocate the soltime and solname arrays so we
  ! can concatenate this solution onto those arrays
  !
  call reallocate(soltime,cgns_zone_counter+1)
  call reallocate(solname,cgns_zone_counter+1)
  soltime(size(soltime)) = real(time*time_ref,kind(soltime))
  solname(size(solname)) = tname
  !
  ! Create a link to the solution node in the current solution file
  !
  call cg_goto_f(ifile_n,ibase_n,cgierr, &
                 "Zone_t",izone_n, &
                 CGNS_GOTO_END)
  call cgns_error(pname,cgns_grid_file,"cg_goto_f",cgierr,__LINE__,__FILE__)
  !
  link_path = base_name//"/"//zone_name//"/"//trim(tname)
  !
  call cg_link_write_f(trim(tname),trim(sol_file),trim(link_path),cgierr)
  call cgns_error(pname,cgns_grid_file,"cg_link_write_f", &
                  cgierr,__LINE__,__FILE__)
  !
  ! Write the solution file
  !
  call write_solution_parallel_cgns(sol_file,tname)
  !
  ! ################################################
  ! #####   FINISHED WRITING FIELD VARIABLES   #####
  ! ################################################
  !
  call set_cgns_queue(STOP_CGNS_QUEUE,pname,cgns_grid_file)
  !
  ! Update the solution times within the base iterative node
  !
  call cg_biter_write_f(ifile_n,ibase_n,'TimeIterValues', &
                        size(soltime,kind=CBT),cgierr)
  call cgns_error(pname,cgns_grid_file,"cg_biter_write_f", &
                  cgierr,__LINE__,__FILE__)
  !
  call cg_goto_f(ifile_n,ibase_n,cgierr, &
                 "BaseIterativeData_t",iiter_n, &
                 CGNS_GOTO_END)
  call cgns_error(pname,cgns_grid_file,"cg_goto_f", &
                  cgierr,__LINE__,__FILE__)
  !
  call cg_array_write_f("TimeValues",CGNS_KIND_REAL,1_CBT, &
                        size(soltime,kind=CST),soltime,ierr)
  call cgns_error(pname,cgns_grid_file,"cg_array_write_f", &
                  cgierr,__LINE__,__FILE__)
  !
  ! Update the names of the zone iterative nodes
  !
  call cg_ziter_write_f(ifile_n,ibase_n,izone_n,"ZoneIterativeData",cgierr)
  call cgns_error(pname,cgns_grid_file,"cg_ziter_write_f", &
                  cgierr,__LINE__,__FILE__)
  !
  call cg_goto_f(ifile_n,ibase_n,cgierr, &
                 "Zone_t",izone_n, &
                 "ZoneIterativeData_t",iiter_n, &
                 CGNS_GOTO_END)
  call cgns_error(pname,cgns_grid_file,"cg_goto_f", &
                  cgierr,__LINE__,__FILE__)
  !
  idata(:) = int( [32,size(solname)] , kind=CST )
  call cg_array_write_f("FlowSolutionPointers",CHARACTER,2_CBT, &
                        idata,solname,cgierr)
  call cgns_error(pname,cgns_grid_file,"cg_array_write_f", &
                  cgierr,__LINE__,__FILE__)
  !
  ! Identify the flow solution within the CGNS file as time dependent
  !
  call cg_simulation_type_write_f(ifile_n,ibase_n,TIMEACCURATE,cgierr)
  call cgns_error(pname,cgns_grid_file,"cg_simulation_type_write_f", &
                  cgierr,__LINE__,__FILE__)
  !
  ! Close the CGNS FILE
  !
  call cgp_close_f(ifile_n,cgierr)
  call cgns_error(pname,cgns_grid_file,"cgp_close_f",cgierr,__LINE__,__FILE__)
  !
  ! Update the zone counter
  !
  cgns_zone_counter = cgns_zone_counter + 1
  !
  write (error_message,100)
  100 format ("Temporary stop after writing initial CGNS files using the ", &
              "new continuous output from the post-processor module.")
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  !
  call debug_timer(leaving_procedure,pname)
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
  !
  !.. Local Scalars ..
  integer :: n,ierr
  !
  integer(CST) :: my_cells_beg,my_cells_end
  integer(CST) :: n_total_output_pts
  integer(CST) :: n_total_output_cells
  !
  character(len=CGLEN) :: varname
  !
  !.. Local Arrays ..
  integer(CBT), dimension(1:nr) :: ixyz_n
  !
  !.. Local Allocatable Arrays ..
  integer(CST), allocatable :: ipelem(:,:)
  real(wp),     allocatable :: vartmp(:,:)
  real(CRT),    allocatable :: var(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_grid_parallel_cgns"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  cgns_grid_file = trim(adjustl(cgns_dir)) // cgns_base_fname
  !
  ! Initialize the variables and data arrays needed for CGNS output
  ! NOTE: This call to routine create_cell_subconnectivity will initialize
  !       the module variables my_pts_beg, my_pts_end, and n_my_output_pts.
  !
  call create_cell_subconnectivity(my_cells_beg,my_cells_end,ipelem)
  !
  ! Extract the total number of pts and cells from the zone_size array
  !
  n_total_output_pts = zone_size(1)
  n_total_output_cells = zone_size(2)
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
  ! Initialize the zone counter
  !
  cgns_zone_counter = 0
  !
  call cgp_open_f(cgns_grid_file,CG_MODE_WRITE,ifile_n,cgierr)
  call cgns_error(pname,cgns_grid_file,"cgp_open_f", &
                  cgierr,__LINE__,__FILE__,"Mode = CG_MODE_WRITE")
  !
  ! Set the CGNS base (who knows what this is for)
  !
  call cg_base_write_f(ifile_n,base_name,int(nr,kind=CBT), &
                       int(nr,kind=CBT),ibase_n,cgierr)
  call cgns_error(pname,cgns_grid_file,"cg_base_write_f", &
                  cgierr,__LINE__,__FILE__)
  !
  ! Write the title for this solution
  !
  call cg_zone_write_f(ifile_n,ibase_n,zone_name,zone_size, &
                       UNSTRUCTURED,izone_n,ierr)
  call cgns_error(pname,cgns_grid_file,"cg_zone_write_f", &
                  cgierr,__LINE__,__FILE__)
  !
  ! Create data nodes for the coordinates
  !
  do n = 1,nr
    varname = "Coordinate"//char_xyz(n)
    call cgp_coord_write_f(ifile_n,ibase_n,izone_n,CGNS_KIND_REAL, &
                           trim(varname),ixyz_n(n),ierr)
    call cgns_error(pname,cgns_grid_file,"cgp_coord_write_f", &
                    cgierr,__LINE__,__FILE__,"Write node: Coordinate",n)
  end do
  !
  ! Now create and write the coordinate data in parallel
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
  if (ncpu == 1) then
    call write_tecplot_file(vartmp,ipelem)
  end if
  !
  ! Write the coordinates for the output points to the CGNS grid file
  !
  do n = 1,nr
    !
    ! Copy the current coordinate variable into the array var
    ! and convert to the real type required for the CGNS file.
    !
    var(:) = real( vartmp(n,:) , kind=kind(var) )
    !
    ! Reset the queue.
    !
    call set_cgns_queue(BEGIN_CGNS_QUEUE,pname,cgns_grid_file)
    !
    ! Add the data for current coordinate variable to the output queue.
    !
    call cgp_coord_write_data_f(ifile_n,ibase_n,izone_n,ixyz_n(n), &
                                my_pts_beg,my_pts_end,var,ierr)
    call cgns_error(pname,cgns_grid_file,"cgp_coord_write_data_f", &
                    cgierr,__LINE__,__FILE__,"Write data: Coordinate",n)
    !
    ! Finally, flush the queue containing the
    ! current coordinate variable to the CGNS file.
    !
    call flush_cgns_queue(pname,cgns_grid_file)
    !
  end do
  !
  ! Create the data node for the grid connectivity
  !
  if (nr == 2) then
    call cgp_section_write_f(ifile_n,ibase_n,izone_n,"Quadrilaterals", &
                             QUAD_4,1_CST,n_total_output_cells, &
                             0_CBT,isect_n,ierr)
    call cgns_error(pname,cgns_grid_file,"cgp_section_write_f", &
                    cgierr,__LINE__,__FILE__,"Element Type = QUAD_4")
  else if (nr == 3) then
    call cgp_section_write_f(ifile_n,ibase_n,izone_n,"Hexahedrals", &
                             HEXA_8,1_CST,n_total_output_cells, &
                             0_CBT,isect_n,ierr)
    call cgns_error(pname,cgns_grid_file,"cgp_section_write_f", &
                    cgierr,__LINE__,__FILE__,"Element Type = HEXA_8")
  end if
  !
  ! Write the element connectivity in parallel
  !
  ! Reset the queue.
  !
  call set_cgns_queue(BEGIN_CGNS_QUEUE,pname,cgns_grid_file)
  !
  ! Add the element connectivity data to the output queue.
  !
  call cgp_elements_write_data_f(ifile_n,ibase_n,izone_n,isect_n, &
                                 my_cells_beg,my_cells_end,ipelem,ierr)
  call cgns_error(pname,cgns_grid_file,"cgp_elements_write_data_f", &
                  cgierr,__LINE__,__FILE__)
  !
  ! Finally, flush the queue containing the
  ! element connectivity to the CGNS file.
  !
  call flush_cgns_queue(pname,cgns_grid_file)
  !
  ! Deallocate the arrays var, vartmp, and ipelem
  !
  if (allocated(var)) then
    deallocate ( var , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"var",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(vartmp)) then
    deallocate ( vartmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vartmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(ipelem)) then
    deallocate ( ipelem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ipelem",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine write_grid_parallel_cgns
!
!###############################################################################
!
subroutine write_solution_parallel_cgns(sol_file,tname)
  !
  !.. Use Statements ..
  use geovar,  only : nr
  !
  use, intrinsic :: iso_c_binding, only : C_NULL_CHAR
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: sol_file
  character(len=*), intent(in) :: tname
  !
  !.. Local Scalars ..
  integer :: n,ierr
  !
  character(len=CGLEN) :: sect_name
  character(len=170)   :: cgns_sol_file
  character(len=300)   :: link_path
  !
  !.. Local Arrays ..
  integer(CBT) :: ifile_s,ibase_s,izone_s,isol_s
  !
  !.. Local Allocatable Arrays ..
  real(wp),       allocatable :: vartmp(:,:)
  real(CRT), allocatable :: var(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_solution_parallel_cgns"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  cgns_sol_file = trim(adjustl(cgns_dir)) // trim(adjustl(sol_file))
  !
  ! Open the cgns output file if iflag is equal to open_CGNS
  !
  ! Initialize the zone counter
  !
  call cgp_open_f(cgns_sol_file,CG_MODE_WRITE,ifile_s,cgierr)
  call io_error(pname,cgns_sol_file,1,__LINE__,__FILE__,cgierr)
  !
  ! Set the CGNS base (who knows what this is for)
  !
  call cg_base_write_f(ifile_s,base_name,int(nr,kind=CBT), &
                       int(nr,kind=CBT),ibase_s,cgierr)
  call cgns_error(pname,cgns_sol_file,"cg_base_write_f", &
                  cgierr,__LINE__,__FILE__)
  !
  ! Write the title for this solution
  !
  call cg_zone_write_f(ifile_s,ibase_s,zone_name,zone_size, &
                       UNSTRUCTURED,izone_s,ierr)
  call cgns_error(pname,cgns_sol_file,"cg_zone_write_f", &
                  cgierr,__LINE__,__FILE__)
  !
  ! Create a link to the grid coordinates in the grid file
  !
  call cg_goto_f(ifile_s,ibase_s,cgierr, &
                 "Zone_t",izone_s, &
                 CGNS_GOTO_END)
  call cgns_error(pname,cgns_sol_file,"cg_goto_f",cgierr,__LINE__,__FILE__)
  !
  link_path = base_name//"/"//zone_name//"/GridCoordinates"
  !
  call cg_link_write_f("GridCoordinates",trim(cgns_base_fname), &
                       trim(link_path),ierr)
  call cgns_error(pname,cgns_sol_file,"cg_link_write_f(GridCoordinates)", &
                  cgierr,__LINE__,__FILE__)
  !
  ! Create a link to the element connectivity in the grid file
  !
  if (nr == 2) then
    sect_name = "Quadrilaterals"
  else if (nr == 3) then
    sect_name = "Hexahedrals"
  end if
  link_path = base_name//"/"//zone_name//"/"//trim(sect_name)
  !
  call cg_link_write_f(trim(sect_name),trim(cgns_base_fname), &
                       trim(link_path),ierr)
  call cgns_error(pname,cgns_sol_file, &
                  "cg_link_write_f("//trim(sect_name)//")", &
                  cgierr,__LINE__,__FILE__)
  !
  ! Allocate var and vartmp if they arent already allocated
  !
  allocate ( vartmp(1:npvar,1:n_my_output_pts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"vartmp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( var(1:n_my_output_pts) , source=real(0,CRT) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"var",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Collect the solution into vartmp
  !
  call collect_solution_for_output(vartmp)
  !
  ! Write the flow variables
  !
  call cg_sol_write_f(ifile_s,ibase_s,izone_s,tname,VERTEX,isol_s,cgierr)
  call cgns_error(pname,cgns_sol_file,"cg_sol_write_f",cgierr,__LINE__,__FILE__)
  !
  ! ################################################
  ! #####   WRITING SOLUTION FIELD VARIABLES   #####
  ! ################################################
  !
  ! Create data nodes for the solution variables
  !
  do n = 1,size(plot_var)
    call cgp_field_write_f(ifile_s,ibase_s,izone_s,isol_s,CGNS_KIND_REAL, &
                           trim(plot_var(n)%varname),isvar_s(n),ierr)
    call cgns_error(pname,cgns_sol_file,"cgp_field_write_f", &
                    cgierr,__LINE__,__FILE__,"Write node: Solution variable",n)
  end do
  !
  ! Now create the output data and write it to the CGNS file in parallel
  !
  do n = 1,size(plot_var)
    !
    ! Save this plot function into a separate contiguous array
    !
    var(:) = real( vartmp(n,:) , kind=kind(var) )
    !
    ! Reset the queue.
    !
    call set_cgns_queue(BEGIN_CGNS_QUEUE,pname,cgns_sol_file)
    !
    ! Add the data for current variable to the output queue.
    !
    call cgp_field_write_data_f(ifile_s,ibase_s,izone_s,isol_s,isvar_s(n), &
                                my_pts_beg,my_pts_end,var,ierr)
    call cgns_error(pname,cgns_sol_file,"cgp_field_write_data_f", &
                    cgierr,__LINE__,__FILE__,"Write data: Solution variable",n)
    !
    ! Finally, flush the queue containing
    ! the current variable to the CGNS file.
    !
    call flush_cgns_queue(pname,cgns_sol_file)
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
  call cgp_close_f(ifile_s,cgierr)
  call io_error(pname,cgns_grid_file,2,__LINE__,__FILE__,cgierr)
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
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("#############################################",/, &
            "    FOR PROCESSOR #",i0)
  2 format (a," = ",i0)
  !
end subroutine write_solution_parallel_cgns
!
!###############################################################################
!
subroutine collect_coordinates_for_output(var)
  !
  !.. Use Statements ..
  use geovar,  only : nfbnd,ncell,cell
  use geovar,  only : xyz,xyz_nodes
  use flowvar, only : edgexyz,facexyz
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
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "collect_coordinates_for_output"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  cell_loop: do nc = 1,ncell
    !
    np1 = cell(nc)%beg_sp
    np2 = cell(nc)%end_sp
    !
    write (170,1) nc
    do n = np1,np2
      !
      l = output_sol_pts(n)
      write (170,2) n,l,(xyz(nn,n), nn=1,size(var,dim=1))
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
  !
  face_loop: do nf = 1,size(output_face_pts)
    !
    face_idx = output_face_pts(nf)%face_idx
    !
    write (170,3) face_idx
    do k = 1,size(output_face_pts(nf)%out_idx)
      !
      l = output_face_pts(nf)%out_idx(k)
      !
      write (170,4) k,l,(facexyz(face_idx)%v(n,k), n=1,size(var,dim=1))
      var(:,l) = facexyz(face_idx)%v(:,k)
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
    write (170,5) edge_idx
    do k = 1,size(output_edge_pts(ne)%out_idx)
      !
      l = output_edge_pts(ne)%out_idx(k)
      !
      write (170,6) k,l,(edgexyz(edge_idx)%v(n,k), n=1,size(var,dim=1))
      var(:,l) = edgexyz(edge_idx)%v(:,k)
      !
    end do
    !
  end do edge_loop
  !
  !
  !
  write (170,7)
  node_loop: do nn = 1,size(output_node_pts)
    !
    node_idx = output_node_pts(nn)%node_idx
    !
    l = output_node_pts(nn)%out_idx
    !
    write (170,8) node_idx,l,(xyz_nodes(n,node_idx), n=1,size(var,dim=1))
    var(:,l) = xyz_nodes(:,node_idx)
    !
  end do node_loop
  !
  ! Deallocate edgexyz since it is no longer needed
  !
  deallocate ( edgexyz , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"edgexyz",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Deallocate or shrink facexyz depending on if MMS is enabled
  !
  if (mms_opt /= 0) then
    !
    ! MMS is enabled, so shrink facexyz to only the boundary faces because the
    ! coordinates for the face points on the interior faces are no longer needed
    !
    allocate ( temp_facexyz(1:nfbnd) , source=facexyz(1:nfbnd) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"temp_facexyz",1,__LINE__,__FILE__,ierr, &
                     error_message)
    call move_alloc( from=temp_facexyz , to=facexyz )
    !
  else
    !
    ! MMS is not enabled, so deallocate facexyz
    ! entirely since it is no longer needed
    !
    deallocate ( facexyz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facexyz",2,__LINE__,__FILE__,ierr,error_message)
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine collect_coordinates_for_output
!
!###############################################################################
!
subroutine collect_solution_for_output(var)
  !
  !.. Use Statements ..
  use geovar,  only : ncell,cell
  use flowvar, only : usp,dusp
  use flowvar, only : faceusp,facedusp
  use flowvar, only : nodeusp,nodedusp
  use flowvar, only : edgeusp,edgedusp
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(out) :: var
  !
  !.. Local Scalars ..
  integer :: k,l,n,nc,ne,nf,nn,np1,np2
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
  !
  cell_loop: do nc = 1,ncell
    !
    np1 = cell(nc)%beg_sp
    np2 = cell(nc)%end_sp
    !
    do n = np1,np2
      !
      l = output_sol_pts(n)
      !
      var(:,l) = compute_output_variables(usp(:,n),dusp(:,:,n))
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
    do k = 1,size(output_face_pts(nf)%out_idx)
      !
      l = output_face_pts(nf)%out_idx(k)
      !
      var(:,l) = compute_output_variables(faceusp(face_idx)%v(:,k,1), &
                                         facedusp(face_idx)%v(:,:,k,1))
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
    do k = 1,size(output_edge_pts(ne)%out_idx)
      !
      l = output_edge_pts(ne)%out_idx(k)
      !
      var(:,l) = compute_output_variables(edgeusp(edge_idx)%v(:,k), &
                                         edgedusp(edge_idx)%v(:,:,k))
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
    l = output_node_pts(nn)%out_idx
    !
    var(:,l) = compute_output_variables(nodeusp(node_idx)%v(:), &
                                       nodedusp(node_idx)%v(:,:))
    !
  end do node_loop
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine collect_solution_for_output
!
!###############################################################################
!
subroutine get_face_solutions
  !
  !.. Use Statements ..
  use geovar,    only : nfbnd,nface
  use flowvar,   only : faceusp,facedusp
  !
  use flux_mod,  only : interp_cv_to_flxpts
  use flux_mod,  only : interp_dusp_to_flxpts
  use module_bc, only : get_faceusp_bc
  use module_bc, only : get_facedusp_bc
  !
  !.. Local Scalars ..
  integer :: nf
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_face_solutions"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
#ifdef SPECIAL_FOR_PGI
  call faceusp%zero_face_var
  call facedusp%zero_face_vec
#else
  call faceusp%zero
  call facedusp%zero
#endif
  !
  call interp_cv_to_flxpts
  !
  call get_faceusp_bc
  !
  do nf = 1,nfbnd
    faceusp(nf)%v(:,:,1) = faceusp(nf)%v(:,:,2)
  end do
  !
  do nf = nfbnd+1,nface
    faceusp(nf)%v(:,:,1) = half * (faceusp(nf)%v(:,:,1) + faceusp(nf)%v(:,:,2))
  end do
  !
  !
  !
  call interp_dusp_to_flxpts
  !
  call get_facedusp_bc
  !
  do nf = 1,nfbnd
    facedusp(nf)%v(:,:,:,1) = facedusp(nf)%v(:,:,:,2)
  end do
  !
  do nf = nfbnd+1,nface
    facedusp(nf)%v(:,:,:,1) = half * (facedusp(nf)%v(:,:,:,1) + &
                                      facedusp(nf)%v(:,:,:,2))
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_face_solutions
!
!###############################################################################
!
subroutine get_edge_and_node_solutions
  !
  !.. Use Statements ..
  use order_mod,  only : maxSP
  use geovar,     only : nr,ncell,cell
  use flowvar,    only : usp,dusp
  use flowvar,    only : edgeusp,edgedusp
  use flowvar,    only : nodeusp,nodedusp
  use interpolation_mod, only : interp
  !
  use parallel_mod, only : exch_edge_and_node_solutions
  !
  !.. Local Scalars ..
  integer :: nc,n1,n2,np,k,l,m,n,ne,nn,kf,lf
  integer :: this_geom,this_order,edge_order
  real(wp) :: ave
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq)              :: upt
  real(wp), dimension(1:nr,1:nq)         :: dupt
  real(wp), dimension(1:maxSP,1:nq)      :: tusp
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: tdusp
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
  call nodeusp%zero_ne1d_t
  call nodedusp%zero_ne2d_t
  call edgeusp%zero_ne2d_t
  call edgedusp%zero_ne3d_t
#else
  call nodeusp%zero
  call nodedusp%zero
  call edgeusp%zero
  call edgedusp%zero
#endif
  !
  ! Loop over all the cells on the current processor and compute the
  ! contribution to the edge and node solutions from each cell
  !
  cell_loop: do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    np = n2-n1+1         ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    ! Create transpose of the usp array to reduce cache misses
    !
    do k = 1,np
      do m = 1,nq
        tusp(k,m) = usp(m,k+n1-1)
      end do
    end do
    !
    ! Create transpose of the dusp array to reduce cache misses
    !
    do k = 1,np
      do m = 1,nq
        do l = 1,nr
          tdusp(k,l,m) = dusp(l,m,k+n1-1)
        end do
      end do
    end do
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
        lf = cell(nc)%edge(n)%get_ep_idx(kf)
        !
        do m = 1,nq
          do l = 1,nr
            dupt(l,m) = interp(this_geom,this_order)% &
                                  toEdge(edge_order)% &
                        dot(lf,tdusp(1:np,l,m))
          end do
        end do
        !
        do m = 1,nq
          upt(m) = interp(this_geom,this_order)% &
                             toEdge(edge_order)% &
                   dot(lf,tusp(1:np,m))
        end do
        !
        edgeusp(ne)%v(1:nq,kf) = edgeusp(ne)%v(1:nq,kf) + ave*upt(1:nq)
        !
        edgedusp(ne)%v(1:nr,1:nq,kf) = edgedusp(ne)%v(1:nr,1:nq,kf) + &
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
      do m = 1,nq
        do l = 1,nr
          dupt(l,m) = interp(this_geom,this_order)% &
                                            toNode% &
                      dot(lf,tdusp(1:np,l,m))
        end do
      end do
      !
      do m = 1,nq
        upt(m) = interp(this_geom,this_order)% &
                                       toNode% &
                 dot(lf,tusp(1:np,m))
      end do
      !
      nodeusp(nn)%v(1:nq) = nodeusp(nn)%v(1:nq) + ave*upt(1:nq)
      !
      nodedusp(nn)%v(1:nr,1:nq) = nodedusp(nn)%v(1:nr,1:nq) + &
                                  ave*dupt(1:nr,1:nq)
      !
    end do cell_nodes
    !
  end do cell_loop
  !
  !
  !
  if (ncpu > 1) then
    call exch_edge_and_node_solutions(nodeusp,nodedusp,edgeusp,edgedusp)
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
  use geovar,    only : pst_geom,nr
  use flowvar,   only : faceusp
  use flowvar,   only : edgeusp,edgedusp
  use flowvar,   only : nodeusp,nodedusp
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
  integer :: mxpt,my_pts_offset
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
    call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_INTEGER, &
                       cpu_mxpts,1_int_mpi,mpi_inttyp, &
                       MPI_COMM_WORLD,mpierr)
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
  nn = size(nodeusp)
  allocate ( node_already_done(1:nn) , source=fals , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"node_already_done",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ne = size(edgeusp)
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
  call nodeusp%zero_ne1d_t
  call nodedusp%zero_ne2d_t
  call edgeusp%zero_ne2d_t
  call edgedusp%zero_ne3d_t
  call faceusp%zero_face_var
#else
  call nodeusp%zero
  call nodedusp%zero
  call edgeusp%zero
  call edgedusp%zero
  call faceusp%zero
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
        lf = cell(nc)%face(n)%get_fp_idx(kf)
        !
        faceusp(nf)%v(1,kf,sf) = real( cell_pts( face_pt(lf) ) , kind=wp )
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
      do kf = 1,cell(nc)%edge(n)%order+1
        !
        lf = cell(nc)%edge(n)%get_ep_idx(kf)
        !
        edgeusp(ne)%v(1,kf) = real( cell_pts( edge_pt(lf) ) , kind=wp )
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
      nodeusp(nn)%v(1) = real( cell_pts( node_pt(lf) ) , kind=wp )
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
    call exch_connectivity(nodeusp,nodedusp,edgeusp,edgedusp,faceusp)
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
        lf = cell(nc)%face(n)%get_fp_idx(kf)
        !
        cell_pts( face_pt(lf) ) = nint( faceusp(nf)%v(1,kf,1) , &
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
      !
      do kf = 1,cell(nc)%edge(n)%order+1
        !
        lf = cell(nc)%edge(n)%get_ep_idx(kf)
        !
        cell_pts( edge_pt(lf) ) = nint( edgeusp(ne)%v(1,kf) , &
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
      cell_pts( node_pt(lf) ) = nint( nodeusp(nn)%v(1) , kind=kind(cell_pts) )
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
  nf = size(faceusp)
  allocate ( my_faces(1:nf) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"my_faces",1,__LINE__,__FILE__,ierr,error_message)
  !
  ne = size(edgeusp)
  allocate ( my_edges(1:ne) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"my_edges",1,__LINE__,__FILE__,ierr,error_message)
  !
  nn = size(nodeusp)
  allocate ( my_nodes(1:nn) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"my_nodes",1,__LINE__,__FILE__,ierr,error_message)
  !
  bad_indices = 0
  !
  ! Find the face points for which this processor is responsible
  !
  find_my_faces: do nf = 1,size(faceusp)
    !
    n = size(faceusp(nf)%v,dim=2)
    m = count(faceusp(nf)%v(1,:,1) > zero)
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
  find_my_edges: do ne = 1,size(edgeusp)
    !
    n = size(edgeusp(ne)%v,dim=2)
    m = count(edgeusp(ne)%v(1,:) > zero)
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
  find_my_nodes: do nn = 1,size(nodeusp)
    !
    if (nodeusp(nn)%v(1) > zero) then
      my_nodes(nn) = true
    end if
    !
  end do find_my_nodes
  !
  ! Reduce the maximum value of bad_indices across all processors
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,bad_indices,1_int_mpi,mpi_inttyp, &
                       MPI_MAX,MPI_COMM_WORLD,mpierr)
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
    call mpi_allreduce(MPI_IN_PLACE,mxpt,1_int_mpi,mpi_inttyp, &
                       MPI_MAX,MPI_COMM_WORLD,mpierr)
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
  write (iout,100) "                     mxpt = ",mxpt
  write (iout,100) "count(  cell_points <= 0) = ",count(cell_points<=0)
  write (iout,100) "count(remaining_pts == 1) = ",count(remaining_pts==1)
  !
  ! Collectively reduce the remaining points on all processors
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,remaining_pts,int(mxpt,kind=int_mpi), &
                       mpi_inttyp,MPI_MAX,MPI_COMM_WORLD,mpierr)
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
      call mpi_bcast(remaining_pts,size(remaining_pts,kind=int_mpi), &
                     mpi_inttyp,this_cpu,MPI_COMM_WORLD,mpierr)
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
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    do n = 1,ncpu
      !
      flush (iout)
      call mpi_barrier(MPI_COMM_WORLD,mpierr)
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
      call mpi_barrier(MPI_COMM_WORLD,mpierr)
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
    loc_n_my_output_pts = ncell*nsp + size(faceusp)*nfp + &
                                      size(edgeusp)*nep + &
                                      size(nodeusp)
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
        lf = cell(nc)%face(n)%get_fp_idx(kf)
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
      nep = cell(nc)%edge(n)%order+1
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
        lf = cell(nc)%edge(n)%get_ep_idx(kf)
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
    call mpi_allreduce(MPI_IN_PLACE,bad_indices,1_int_mpi,mpi_inttyp, &
                       MPI_MAX,MPI_COMM_WORLD,mpierr)
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
    call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_INTEGER, &
                       cpu_mxpts,1_int_mpi,mpi_inttyp, &
                       MPI_COMM_WORLD,mpierr)
   !call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_DATATYPE_NULL, &
    call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_INTEGER, &
                       cpu_subcells,1_int_mpi,mpi_inttyp, &
                       MPI_COMM_WORLD,mpierr)
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
pure function compute_output_variables(upt,dupt) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:),   intent(in) :: upt
  real(wp), dimension(:,:), intent(in) :: dupt
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
    return_value(i) = compute_output_variable(upt,dupt,plot_var(i))
  end do
  !
end function compute_output_variables
!
!###############################################################################
!
pure function compute_output_variable(upt,dupt,plot_variable) &
                               result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : rhoref,aref,tref,time_ref
  !
  !.. Formal Arguments ..
  real(wp), dimension(:),   intent(in) :: upt
  real(wp), dimension(:,:), intent(in) :: dupt
  type(plot_func_t),        intent(in) :: plot_variable
  !
  !.. Function Result ..
  real(wp) :: return_value
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
    return_value = pressure_cv_sp( upt )              * rhoref * aref * aref
    !
  else if (plot_variable == plot_static_temperature) then
    !
    return_value = temperature_cv_sp( upt )           * tref
    !
  else if (plot_variable == plot_total_energy) then
    !
    return_value = upt(nee)                           * rhoref * aref * aref
    !
  else if (plot_variable == plot_mach_number) then
    !
    return_value = mach_number_cv_sp( upt )
    !
  else if (plot_variable == plot_entropy) then
    !
    return_value = entropy_cv_sp( upt )
    !
  else if (plot_variable == plot_vorticity_x) then
    !
    return_value = vorticity_dusp( upt , dupt , 1 )   / time_ref
    !
  else if (plot_variable == plot_vorticity_y) then
    !
    return_value = vorticity_dusp( upt , dupt , 2 )   / time_ref
    !
  else if (plot_variable == plot_vorticity_z) then
    !
    return_value = vorticity_dusp( upt , dupt , 3 )   / time_ref
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
  integer :: ierr
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
  !
  l = 0
  do n = 1,size(edge_values,dim=3)
    do m = 1,size(edge_values,dim=2)
      l = l + 1
      edge_values(:,m,n) = [l,map%ep(l)+1]
    end do
  end do
  !
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
pure function vorticity_dusp(usp,dusp,dir) result(return_value)
  !
  integer,                  intent(in) :: dir
  real(wp), dimension(:),   intent(in) :: usp
  real(wp), dimension(:,:), intent(in) :: dusp
  !
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,n1,n2
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(dusp,dim=1),1:size(dusp,dim=2)) :: dpvdx
  !
  !.. Local Parameters ..
  integer, parameter, dimension(3) :: idx = [2,3,1]
  integer, parameter, dimension(3) :: jdx = [3,1,2]
  !
continue
  !
  n1 = size(dusp,dim=1)
  n2 = size(dusp,dim=2)
  !
  i = idx(dir)
  j = jdx(dir)
  !
  dpvdx(1:n1,1:n2) =  grad_cv_to_grad_pv_sp( usp(:) , dusp(:,:) )
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
  use ovar, only : lustre_stripe_count,lustre_stripe_size
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
  !
  call mpi_bcast(clen,1_int_mpi, &
                 mpi_inttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
  !
  call mpi_bcast(cgns_dir(1:clen),int(clen,kind=int_mpi), &
                 MPI_CHARACTER,0_int_mpi,MPI_COMM_WORLD,mpierr)
  !
  ! Create the format string for creating the CGNS file names
  !
  write (char_num,1) itrst + num_timesteps
  final_iter_digits = max(4,len_trim(adjustl(char_num)))
  !
  write (cgnsfile_fmt,2) final_iter_digits
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (i0)
  2 format ("('solution.',i0.",i0,",'.cgns)'")
  3 format ("lfs setstripe -c ",i0," -S ",a,1x,a)
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
  3 format (8i6)
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
include "Functions/pressure.f90"
include "Functions/temperature.f90"
include "Functions/mach_number.f90"
include "Functions/entropy.f90"
include "Functions/grad_cv_to_grad_pv.f90"
include "Functions/cell_solpts.f90"
!
!###############################################################################
!
end module postprocessor_mod
