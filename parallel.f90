module parallel_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nq
  !
  implicit none
  !
  private
#ifdef METIS_5
  include "metis.fh"
#endif
  !
  !.. Public Module Procedures ..
  !
  public :: create_serial_cell_map
  public :: partition_grid
  public :: collect_global_solpts
  public :: collect_global_solution
  public :: wait_faceusp
  public :: exchange_faceusp
  public :: wait_facedusp
  public :: exchange_facedusp
  public :: exch_edge_and_node_solutions
  public :: parallel_memory_usage
 !public :: exchange_modes
  !
  ! IDEA: 1. all processors send in their global versions of bface, xyz_nodes, 
  !          nodes_of_cell, nodes_of_cell_ptr, etc. to a single subroutine
  !       2. subroutine makes a copy of the global information necessary for
  !          output within the module
  !       3. Subroutine uses METIS to partition the grid
  !       4. Each processor then creates local versions of the arrays input
  !          to the subroutine which is then output from the subroutine
  !       5. Datatypes are needed to perform a parallel reduction of data 
  !          from the local partitions to a single solution array for output
  !       6. Datatypes are needed to exchange information along partition
  !          boundaries between processors
  !       7. We will probably have to create the solution points within each
  !          after the partitions are created and localized to the processors
  !
  integer, save :: nr
  integer, save :: nbfai
  !
  ! Dimensions of the global grid
  !
  integer, save :: ngcell
  integer, save :: ngnode
  integer, save :: ngfbnd
  integer, save :: ngface
  integer, save :: ngtotcel
  integer, save :: ngtotnod
  !
  ! Number of communication partners
  !
  integer, save :: ncomm
  !
  ! Exchange partners
  !
  integer, save, allocatable, dimension(:) :: exch_cpu
  !
  ! Test MPI datatypes for exchanging information
  !
  integer(int_mpi), save, allocatable, dimension(:) :: testf_send
  integer(int_mpi), save, allocatable, dimension(:) :: testf_recv
  !
  ! Request and status arrays for non-blocking exchanges
  !
  integer(int_mpi), save, allocatable, dimension(:)   :: fusp_rqst
  integer(int_mpi), save, allocatable, dimension(:,:) :: fusp_istat
  !
  integer(int_mpi), save, allocatable, dimension(:)   :: fdusp_rqst
  integer(int_mpi), save, allocatable, dimension(:,:) :: fdusp_istat
  !
  ! Exchange buffer
  !
  real(r8), save, allocatable, dimension(:) :: fusp_send_buffer
  real(r8), save, allocatable, dimension(:) :: fusp_recv_buffer
  real(r8), save, allocatable, dimension(:) :: fdusp_send_buffer
  real(r8), save, allocatable, dimension(:) :: fdusp_recv_buffer
  !
  integer, save :: eusp_size
  integer, save :: edusp_size
  integer, save :: nusp_size
  integer, save :: ndusp_size
  !
  real(r8), save, allocatable, dimension(:) :: neusp_buffer
  !
  ! MPI datatype for collecting the global
  ! solution point coordinates on all processors
  !
  integer(int_mpi), save, allocatable, dimension(:) :: recv_xyz
  !
  ! MPI datatype for collecting the global solution on all processors
  !
  integer(int_mpi), save, allocatable, dimension(:) :: recv_sol
  !
  ! MPI datatype for collecting a global variable on all processors
  !
  integer(int_mpi), save, allocatable, dimension(:) :: recv_var
  !
  !
  !
  ! ##########################################################
  ! #####   MAPPING DERIVED TYPE AND ASSOCIATED ARRAYS   #####
  ! ##########################################################
  !
  ! MAP_T : derived type to store mappings between local and global
  !         indices for cells, faces, and nodes
  !
  type :: map_t
    ! loc_to_glb : local index to global index mapping array
    integer, allocatable :: loc_to_glb(:)
  end type map_t
  !
  ! CELL_MAP : mappings between local and global cell indices
  !            NOTE: I believe this mapping array is the only one needed outside
  !                  this module. For example, it is needed for initializing the
  !                  local solution from the global solution that is read in 
  !                  from the restart file.
  !
  type(map_t), public, save, allocatable :: cell_map(:)
  !
  !
  !
  !
  ! #####################################################################
  ! #####   FACE CONNECTIVITY DERIVED TYPES AND ASSOCIATED ARRAYS   #####
  ! #####################################################################
  !
  ! CELL_ON_FACE_T : derived type to store the information about
  !                  a cell on one side of a given face
  !                  NOTE: This is similar to the datatype of the same name
  !                        in the module geovar but with some slight changes.
  !
  type :: cell_on_face_t
    ! cell : index of the grid cell on the current side of the face
    integer :: cell
    ! cell_face : index for the side/surface of the cell that comprises
    !             the current face
    !               Examples: for a triangle, in the range 1:3, quad 1:4,
    !                         tetr 1:4, pyra 1:5, pris 1:5, hexa 1:6
    integer :: cell_face
    ! geom : geometry of this cell
    integer :: geom
    ! order : order of this cell
    integer :: order
    ! partition : grid partition/processor on which this cell resides
    integer :: partition
    ! flxpts : indices for the flux points of this cell on the current face
    !          NOTE: In the version of the datatype cell_on_face_t in the module
    !                geovar, flxpts gives an index that is relative to all of
    !                the flux points for a given cell.
    !                Where flxpts is used in this module, we only care about the
    !                ordering on the current face so flxpts will contain indices
    !                from 1 to size(flxpts).
    integer, allocatable :: flxpts(:)
  end type cell_on_face_t
  !
  ! FACE_T : derived type to store information for each grid face
  !          NOTE: This is similar to the datatype of the same name
  !                in the module geovar but with some slight changes.
  !
  type :: face_t
    ! geom : geometry type of the current face
    integer :: geom
    ! order : solution order of the current face
    !         NOTE: This should be the maximum solution order between the two
    !               cells attached to the current face.
    integer :: order
    ! nodes : array containing the indices for the grid nodes that define this
    !         face. At most this array will contain four node indices for a
    !         quadrilateral face. For triangular or edge faces, the remaining
    !         unused nodes will remain zero. For example:
    !            edge_face => face(i)%nodes(1:4) = [1,2,0,0]
    !            tria_face => face(i)%nodes(1:4) = [1,2,3,0]
    !            quad_face => face(i)%nodes(1:4) = [1,2,3,4]
    integer, dimension(1:4) :: nodes = 0
    ! left : Information for the cell on the left side of the current face.
    !        NOTE: The normal to the current face always points AWAY FROM the
    !              left cell on the face.
    !        NOTE: If the current face is a boundary face, the left side is the
    !              host/interior cell.
    type(cell_on_face_t) :: left
    ! right : Information for the cell on the right side of the current face.
    !         NOTE: The normal to the current face always points TOWARDS / INTO
    !               the right cell on the face.
    !         NOTE: If the current face is a boundary face, the right side is
    !               the ghost cell.
    type(cell_on_face_t) :: right
    ! side : Possible replacement for left and right components above. This
    !        would allow for the cell_idx component of face_of_cell_t to be
    !        removed since it is a duplicate of the flxpts component of
    !        cell_on_face_t. The side component of face_of_cell_t would give
    !        either the value 1 or 2 which would correspond to the location in
    !        this side component of face_t.
    !        side(1) is the same as the left  component
    !        side(2) is the same as the right component
   !type(cell_on_face_t) :: side(1:2)
  end type face_t
  !
  ! FACE : Array of type face_t that contains information for all grid faces
  !
  type(face_t), save, allocatable :: face(:)
  !
  !
  !
  ! #####################################################################
  ! #####   MPI COMMUNICATION DERIVED TYPES AND ASSOCIATED ARRAYS   #####
  ! #####################################################################
  !
  ! EXCH_DATA_T : Derived type containing the information needed to communicate
  !               data between two processes
  !
  type :: exch_data_t
    ! adj_cpu : process/CPU ID for the communication partner
    integer(int_mpi) :: adj_cpu
    ! send_tag : message tag for matching the MPI communication
    integer(int_mpi) :: send_tag
    ! recv_tag : message tag for matching the MPI communication
    integer(int_mpi) :: recv_tag
    ! send_datatype : handle for the MPI user-defined datatype to send data
    integer(int_mpi) :: send_datatype = MPI_DATATYPE_NULL
    ! recv_datatype : handle for the MPI user-defined datatype to receive data
    integer(int_mpi) :: recv_datatype = MPI_DATATYPE_NULL
  end type exch_data_t
  !
  ! EXCH_FUSP : Array containing all the information needed to communicate the
  !             face flux point solution variables between adjacent processes.
  !
  type(exch_data_t), save, allocatable :: exch_fusp(:)
  !
  ! EXCH_FDUSP : Array containing all the information needed to communicate
  !              the face flux point gradients between adjacent processes.
  !
  type(exch_data_t), save, allocatable :: exch_fdusp(:)
  !
  ! CPU_T : Derived type containing the information for one CPU that
  !         constitutes half of a communication pair
  !
  type :: cpu_t
    ! cpu_pnum : process/CPU ID for one of a communication pair
    integer(int_mpi) :: cpu_pnum
    ! send_datatype : handle for the MPI user-defined datatype to send data
    integer(int_mpi) :: send_datatype = MPI_DATATYPE_NULL
    ! recv_datatype : handle for the MPI user-defined datatype to receive data
    integer(int_mpi) :: recv_datatype = MPI_DATATYPE_NULL
  end type cpu_t
  !
  ! EXCH_PAIR_T : Derived type containing information for each of the
  !               communication pairs
  !
  type :: exch_pair_t
    ! tag1 : tag for the cpu1->cpu2 communication
    integer(int_mpi) :: tag1
    ! tag2 : tag for the cpu2->cpu1 communication
    integer(int_mpi) :: tag2
    ! cpu : information for the two cpus making up this pair
    type(cpu_t) :: cpu(1:2)
  end type exch_pair_t
  !
  ! EXCH_NEUSP : Array containing the information for all the communication
  !              pairs that need to exchange edge and node solutions
  !              
  type(exch_pair_t), save, allocatable :: exch_neusp(:)
  !
  !
  ! #########################################
  ! #####   MODULE GENERIC INTERFACES   #####
  ! #########################################
  !
 !interface exchange_modes
 !  module procedure exchange_modes_r4, exchange_modes_r8
 !end interface exchange_modes
  !
contains
!
!###############################################################################
!
subroutine create_serial_cell_map(ncell)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: ncell
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_serial_cell_map"
  !
continue
  !
  allocate ( cell_map(1:1) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_map",1,ierr,error_message,skip_alloc_pause)
  !
  allocate ( cell_map(1)%loc_to_glb(1:ncell) , source=intseq(1,ncell) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_map(1)%loc_to_glb",1,ierr,error_message)
  !
end subroutine create_serial_cell_map
!
!###############################################################################
!
subroutine partition_grid(npart,metis_option_requested,bface,nodes_of_cell, &
                          nodes_of_cell_ptr,xyz_nodes,cell_geom,cell_order)
  !
  !.. Formal Arguments ..
  integer,                                  intent(in) :: npart
  integer,                                  intent(in) :: metis_option_requested
  integer,  allocatable, dimension(:),   intent(inout) :: cell_geom
  integer,  allocatable, dimension(:),   intent(inout) :: cell_order
  integer,  allocatable, dimension(:),   intent(inout) :: nodes_of_cell
  integer,  allocatable, dimension(:),   intent(inout) :: nodes_of_cell_ptr
  integer,  allocatable, dimension(:,:), intent(inout) :: bface
  real(wp), allocatable, dimension(:,:), intent(inout) :: xyz_nodes
  !
  !.. Local Scalars ..
  integer :: n,ierr
  character(len=200) :: array_name
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:) :: xadj
  integer, allocatable, dimension(:) :: adjncy
  integer, allocatable, dimension(:) :: epart
  !
  !.. Local Allocatable Derived-Type Arrays ..
  type(map_t),  allocatable, dimension(:) :: node_map
  type(map_t),  allocatable, dimension(:) :: face_map
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "partition_grid"
  !
continue
  !
  ! Get the size of the global grid
  !
  nr = size(xyz_nodes,dim=1)
  nbfai = size(bface,dim=1)
  !
  ! The total number of nodes and cells, including ghost cells
  !
  ngtotcel = size(nodes_of_cell_ptr) - 1
  ngtotnod = size(xyz_nodes,dim=2)
  !
  ! The number of boundary faces and interior cells
  !
  ngfbnd = size(bface,dim=2)
  ngcell = ngtotcel - ngfbnd
  !
  ! Get the number of nodes that define the interior cells
  !
  n = nodes_of_cell_ptr(ngcell+1)
  ngnode = maxval(nodes_of_cell(1:n))
  !
  ! Create the dual graph for the global grid not including ghost cells
  !
  call memory_pause("Before: calling create_dual_graph")
  call create_dual_graph(cell_geom,cell_order,nodes_of_cell, &
                         nodes_of_cell_ptr,bface,xadj,adjncy)
  !
  ! Allocate the epart array to store the cell partition numbers
  !
  allocate ( epart(1:ngcell) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"epart",1,ierr,error_message)
 !if (mypnum == 0) then
 !  call output_partitions(99,nodes_of_cell,nodes_of_cell_ptr,xyz_nodes)
 !end if
  !
  ! Use METIS to partition the grid
  !
#ifdef METIS_4
  !
  call memory_pause("Before: calling partition_using_metis_4")
  call partition_using_metis_4(metis_option_requested,npart, &
                               ngcell,cell_geom,cell_order, &
                               xadj,adjncy,epart)
  !
#elif METIS_5
  !
  call memory_pause("Before: calling partition_using_metis_5")
  call partition_using_metis_5(metis_option_requested,npart, &
                               ngcell,xadj,adjncy,epart)
  !
#else
  !
  write (error_message,2)
  call stop_recon(stop_mpi,pname,error_message)
  !
#endif
  !
  ! Deallocate no longer needed METIS arrays
  !
  deallocate ( xadj , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xadj",2,ierr,error_message)
  deallocate ( adjncy , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"adjncy",2,ierr,error_message)
  !
  ! Distribute partitions across all processors
  !
  call mpi_bcast(epart,int(ngcell,kind=kind(int_mpi)), &
                 mpi_inttyp,0,mpi_comm_world,mpierr)
  !
  ! Output the partition grid to a tecplot file
  !
 !if (mypnum == 0) then
 !  call output_partitions(mypnum,nodes_of_cell,nodes_of_cell_ptr, &
 !                         xyz_nodes,epart)
 !end if
 !call stop_recon(stop_mpi,pname,"outputing partitions")
  !
  ! Find the faces that separate two adjacent partitions
  ! and add these faces to boundary faces array
  !
  call memory_pause("Before: calling find_partition_boundaries")
  call find_partition_boundaries(epart,bface)
  !
  ! Now that we have created boundary conditions for faces
  ! located on partition boundaries, use epart to create
  ! mappings from local to global block orderings.
  ! This also includes creating the reverse mapping arrays
  ! or the maps from global to local block orderings
  !
  call memory_pause("Before: calling create_grid_mappings")
  call create_grid_mappings(npart,epart, &
                            nodes_of_cell,nodes_of_cell_ptr, &
                            node_map,face_map,cell_map)
  !
 !call parallel_memory_usage(iout)
  !
  ! Now that the mappings have been created, we
  ! need to localize the grid on each processor
  !
  call memory_pause("Before: calling localize_grid")
  call localize_grid(nodes_of_cell,nodes_of_cell_ptr,bface, &
                     xyz_nodes,epart,cell_geom,cell_order, &
                     node_map,face_map,cell_map)
  !
  ! Deallocate epart since it is no longer needed
  !
  if (allocated(epart)) then
    deallocate ( epart , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"epart",2,ierr,error_message)
  end if
  !
  ! Deallocate the node and edge arrays and the node and edge
  ! mapping arrays since they are no longer needed
  !
  if (allocated(node_map)) then
    deallocate ( node_map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"node_map",2,ierr,error_message)
  end if
  !
  ! Finally, create MPI datatypes for exchanging face data between processes
  !
  call memory_pause("Before: calling create_face_nonblocking_datatypes")
  call create_face_nonblocking_datatypes(bface,face_map)
  !
  ! Deallocate the face array because I dont think it is needed anymore
  !
  if (allocated(face)) then
    deallocate ( face , stat=ierr ,errmsg=error_message )
    call alloc_error(pname,"face",2,ierr,error_message)
  end if
  !
  if (allocated(face_map)) then
    deallocate ( face_map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"face_map",2,ierr,error_message)
  end if
  !
  ! Format Statements
  !
  2 format(" It seems that the code was compiled without either of the METIS", &
           " pre-processor options being defined. One of these options needs", &
           " to be defined in order to partition the grid and allow the use", &
           " of multiple processors.")
  !
end subroutine partition_grid
!
!###############################################################################
!
subroutine partition_using_metis_4(metis_option_requested,npart, &
                                   ngcell,cell_geom,cell_order, &
                                   xadj,adjncy,epart)
  !
  !.. Use Statements ..
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  !
  !.. Formal Arguments ..
  integer,                  intent(in) :: metis_option_requested
  integer,                  intent(in) :: npart
  integer,                  intent(in) :: ngcell
  integer, dimension(:),    intent(in) :: cell_geom
  integer, dimension(:),    intent(in) :: cell_order
  integer, dimension(:),    intent(in) :: xadj
  integer, dimension(:),    intent(in) :: adjncy
  integer, dimension(:), intent(inout) :: epart
  !
  !.. Local Scalars ..
  integer :: n,ip,ierr
  integer :: wgtflag,numflag,vol
  integer :: metis_option
  logical :: passes_inspection
  character(len=19) :: method
  character(len=12) :: weight
  character(len=20) :: minmzd
  !
  !.. Local Arrays ..
  integer,       dimension(1:5)     :: metis_option_hierarchy
  integer,       dimension(1:5)     :: options
  real(c_float), dimension(1:npart) :: tpwgts
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:) :: vwgt
  integer, allocatable, dimension(:) :: adjwgt
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "partition_using_metis_4"
  !
continue
#ifdef METIS_4
  !
  ! Create the priority list for METIS options to try
  !
  metis_option_hierarchy = [metis_option_requested,3,2,1,4]
  !
  ! Allocate the weighting arrays
  !
  allocate ( vwgt(1:ngcell) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"vwgt",1,ierr,error_message)
  !
  allocate ( adjwgt(1:last(xadj)-1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"adjwgt",1,ierr,error_message)
  !
  ! One-pass loop to allow for easy error handling if none of the
  ! partitioning methods produce a quality grid partitioning
  !
  part_loop: do ip = 1,1
    !
    do n = 1,size(metis_option_hierarchy)
      !
      ! Get the METIS option to try based on the priority list
      !
      metis_option = metis_option_hierarchy(n)
      !
      ! Create the weights for the vertices of the dual graph
      !
      call create_grid_weights(metis_option,vwgt,adjwgt, &
                               cell_geom,cell_order)
      !
      ! Define other METIS variables
      !
      numflag = 1
      options(:) = 0
      vol = 0
      !
      ! Determine whether weights are being used and set wgtflag accordingly
      !
      wgtflag = 0
      if (minval(adjwgt) /= maxval(adjwgt)) wgtflag = 1
      if (minval(  vwgt) /= maxval(  vwgt)) wgtflag = wgtflag + 2
      !
      ! Partition the grid using one of the METIS routines
      !
      if (metis_option == 1) then
        !
        ! This METIS routine computes non-weighted partitions using a
        ! multilevel recursive bisection while minimizing the edge cut.
        !
        method = "recursive bisection"
        weight = "non-weighted"
        minmzd = "      edge cut      "
        !
        call metis_partgraphrecursive(ngcell,xadj,adjncy,vwgt,adjwgt,wgtflag, &
                                      numflag,npart,options,vol,epart)
        !
      else if (metis_option == 2) then
        !
        ! This METIS routine computes non-weighted partitions using
        ! a multilevel k-way partitioning algorithm while minimizing
        ! the edge cut.
        !
        method = "k-way partitioning "
        weight = "non-weighted"
        minmzd = "      edge cut      "
        !
        call metis_partgraphkway(ngcell,xadj,adjncy,vwgt,adjwgt,wgtflag, &
                                 numflag,npart,options,vol,epart)
        !
      else if (metis_option == 3) then
        !
        ! This METIS routine computes non-weighted partitions using a
        ! a multilevel k-way partitioning algorithm while minimizing
        ! the total communication volume.
        !
        minmzd = "communication volume"
        weight = "non-weighted"
        method = "k-way partitioning "
        !
        call metis_partgraphvkway(ngcell,xadj,adjncy,vwgt,adjwgt,wgtflag, &
                                  numflag,npart,options,vol,epart)
        !
      else if (metis_option == 4) then
        !
        ! This METIS routine computes weighted partitions using a
        ! multilevel k-way partitioning algorithm while minimizing
        ! the total communication volume.
        !
        method = "k-way partitioning "
        weight = "  weighted  "
        minmzd = "communication volume"
        !
        ! Compute the prescribed partition weights
        ! NOTE: The array tpwgts needs to be the Fortran real
        !       kind that corresponds to the C float type
        !
        !##############################################################
        !##############################################################
        !
        ! For now, set the partition weights equal and make sure they
        ! sum to one. After looking at the METIS source code, this is
        ! actually equivalent to just calling metis_partgraphvkway.
        !
        tpwgts(:) = real(one,c_float)/real(npart,c_float)
        tpwgts(:) = tpwgts(:)/sum(tpwgts)
        !
        !##############################################################
        !##############################################################
        !
        call metis_wpartgraphvkway(ngcell,xadj,adjncy,vwgt,adjwgt,wgtflag, &
                                   numflag,npart,tpwgts,options,vol,epart)
        !
      else
        !
        ! The METIS option requested is invalid!
        !
        write (error_message,103)
        call stop_recon(stop_mpi,pname,error_message)
        !
      end if
      !
      ! I dont trust that the METIS routines will compute identical
      ! partitioning on each processor so have the root process check
      ! its local partitioning for errors and then broadcast that
      ! partitioning to the other processes if none are found.
      !
      if (mypnum == 0) then
        call check_quality_of_partitioning(npart,epart, &
                                           passes_inspection)
      end if
      !
      ! Broadcast the logical passes_inspection to all processors so
      ! that they all will know whether to exit part_loop.
      !
      call mpi_bcast(passes_inspection,1,mpi_logical,0,mpi_comm_world,mpierr)
      !
      ! If the partitioning passes inspection, have the root processor
      ! output the results of the partitioning and exit from part_loop.
      !
      if (passes_inspection) then
        if (mypnum == 0) write (iout,101) method,npart,weight,minmzd,vol
        exit part_loop
      end if
      !
    end do
    !
    write (error_message,102)
    call stop_recon(stop_mpi,pname,error_message)
    !
  end do part_loop
  !
  ! Deallocate the weight arrays before we leave
  !
  deallocate ( vwgt , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"vwgt",2,ierr,error_message)
  deallocate ( adjwgt , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"adjwgt",2,ierr,error_message)
  !
  ! Format statements
  !
  101 format (/,'-------------------------------------------------------',//, &
                '         Used METIS 4.0 to partition the graph         ',/, &
                '    using a multilevel ',a,' algorithm    ',/, &
                '           into ',i0,' ',a,' partitions.          ',/, &
                '     The minimized ',a,' was ',i0,'.    ',//, &
                '-------------------------------------------------------',/)
  102 format (" None of the partitioning methods available in METIS 4.0 were", &
              " able to find a quality partitioning of the grid!")
  103 format (" The METIS partitioning method requested in the input file is", &
              " invalid!")
  !
#endif
end subroutine partition_using_metis_4
!
!###############################################################################
!
subroutine partition_using_metis_5(metis_option_requested,npart, &
                                   ngcell,xadj,adjncy,epart)
  !
  !.. Use Statements ..
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_null_ptr
  !
  !.. Formal Arguments ..
  integer,                  intent(in) :: metis_option_requested
  integer,                  intent(in) :: npart
  integer,                  intent(in) :: ngcell
  integer, dimension(:),    intent(in) :: xadj
  integer, dimension(:),    intent(in) :: adjncy
  integer, dimension(:), intent(inout) :: epart
  !
  !.. Local Scalars ..
  integer :: n,ip,ierr
  integer :: ncon,objval
  integer :: metis_option
  logical :: passes_inspection
  character(len=19) :: method
  character(len=12) :: weight
  character(len=20) :: minmzd
  !
  !.. Local Arrays ..
  integer, dimension(1:5) :: metis_option_hierarchy
#ifdef METIS_5
  integer, dimension(0:METIS_NOPTIONS-1) :: options
#endif
  !
  !.. Local Allocatable Arrays ..
  integer,       allocatable, dimension(:) :: vwgt
  integer,       allocatable, dimension(:) :: vsize
  integer,       allocatable, dimension(:) :: adjwgt
  real(c_float), allocatable, dimension(:) :: tpwgts
  real(c_float), allocatable, dimension(:) :: ubvec
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "partition_using_metis_5"
  !
continue
#ifdef METIS_5
  !
  ! Set the default options for METIS
  !
  call METIS_SetDefaultOptions(options)
  !
  ! Customize some of these options
  !
  ! Set numbering to Fortran style (begins at 1)
  options(METIS_OPTION_NUMBERING) = 1
  ! Force contiguous partitions
  options(METIS_OPTION_CONTIG) = 1
  ! Explicitly minimize the maximum connectivity
  options(METIS_OPTION_MINCONN) = 1
  !
  ! Define other METIS variables
  !
  ncon = 1 ! number of weights associated with each vertex (>= 1)
  !
  ! Allocate the local arrays
  !
#ifdef DEBUG_ON
  allocate ( vwgt(1:ncon*ngcell) , source=1 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"vwgt",1,ierr,error_message)
  allocate ( vsize(1:ngcell) , source=1 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"vsize",1,ierr,error_message)
  allocate ( adjwgt(1:size(adjncy)) , source=1 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"adjwgt",1,ierr,error_message)
  allocate ( tpwgts(0:ncon*npart-1) , &
             source=1.0_c_float/real(npart,c_float) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tpwgts",1,ierr,error_message)
  tpwgts = real( one/real(npart,kind=wp) , kind=kind(tpwgts) )
  tpwgts = tpwgts / sum(tpwgts)
  allocate ( ubvec(1:ncon) , source=1.01_c_float , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ubvec",1,ierr,error_message)
  write (iout,*) mypnum,sum(tpwgts)
#endif
  !
  ! Create the hierarchy of METIS algorithms to try
  !
  metis_option_hierarchy = [metis_option_requested,3,2,1,4]
  !
  ! One-pass loop to allow for easy error handling if none of the
  ! partitioning methods produce a quality grid partitioning
  !
  part_loop: do ip = 1,1
    !
    do n = 1,size(metis_option_hierarchy)
      !
      ! Get the METIS option to try based on the priority list
      !
      metis_option = metis_option_hierarchy(n)
      !
      ! Partition the grid using one of the METIS routines
      !
      if (metis_option == 1) then
        !
        ! This METIS routine computes non-weighted partitions using a
        ! multilevel recursive bisection while minimizing the edge cut.
        !
        method = "recursive bisection"
        weight = "non-weighted"
        minmzd = "      edge cut      "
        !
        call metis_partgraphrecursive(ngcell,ncon,xadj,adjncy,vwgt,vsize, &
                                      adjwgt,npart,tpwgts,ubvec,options, &
                                      objval,epart)
        !
      else if (metis_option == 2) then
        !
        ! This METIS routine computes non-weighted partitions using
        ! a multilevel k-way partitioning algorithm while minimizing
        ! the edge cut.
        !
        method = "k-way partitioning "
        weight = "non-weighted"
        minmzd = "      edge cut      "
        !
        options(METIS_OPTION_OBJTYPE) = METIS_OBJTYPE_CUT
        !
        call metis_partgraphkway(ngcell,ncon,xadj,adjncy,vwgt,vsize, &
                                 adjwgt,npart,tpwgts,ubvec,options, &
                                 objval,epart)
        !
      else if (metis_option == 3) then
        !
        ! This METIS routine computes non-weighted partitions using a
        ! a multilevel k-way partitioning algorithm while minimizing
        ! the total communication volume.
        !
        minmzd = "communication volume"
        weight = "non-weighted"
        method = "k-way partitioning "
        !
        options(METIS_OPTION_OBJTYPE) = METIS_OBJTYPE_VOL
        !
        call metis_partgraphkway(ngcell,ncon,xadj,adjncy,vwgt,vsize, &
                                 adjwgt,npart,tpwgts,ubvec,options, &
                                 objval,epart)
        !
      else if (metis_option == 4) then
        !
        ! This METIS routine computes weighted partitions using a
        ! multilevel k-way partitioning algorithm while minimizing
        ! the total communication volume.
        !
        method = "k-way partitioning "
        weight = "  weighted  "
        minmzd = "communication volume"
        !
        !##############################################################
        !##############################################################
        !
        ! NEED TO SET SOMETHING UP TO COMPUTE PARTITION WEIGHTS.
        ! FOR NOW, THE PARTITION WEIGHTS ARE EQUAL AND THIS IS
        ! JUST EQUIVALENT TO (metis_option == 3).
        !
        !##############################################################
        !##############################################################
        !
        options(METIS_OPTION_OBJTYPE) = METIS_OBJTYPE_VOL
        !
        call metis_partgraphkway(ngcell,ncon,xadj,adjncy,vwgt,vsize, &
                                 adjwgt,npart,tpwgts,ubvec,options, &
                                 objval,epart)
        !
      else
        !
        ! The METIS option requested is invalid!
        !
        write (error_message,3)
        call stop_recon(stop_mpi,pname,error_message)
        !
      end if
      !
      ! I dont trust that the METIS routines will compute identical
      ! partitioning on each processor so have the root process check
      ! its local partitioning for errors and then broadcast that
      ! partitioning to the other processes if none are found.
      !
      if (mypnum == 0) then
        call check_quality_of_partitioning(npart,epart, &
                                           passes_inspection)
      end if
      !
      ! Broadcast the logical passes_inspection to all processors so
      ! that they all will know whether to exit part_loop.
      !
      call mpi_bcast(passes_inspection,1,mpi_logical,0,mpi_comm_world,mpierr)
      !
      ! If the partitioning passes inspection, have the root processor
      ! output the results of the partitioning and exit from part_loop.
      !
      if (passes_inspection) then
        if (mypnum == 0) write (iout,1) method,npart,weight,minmzd,objval
        exit part_loop
      end if
      !
      if (mypnum == 0) write (iout,4) method,npart,weight,minmzd
      !
    end do
    !
    write (error_message,2)
    call stop_recon(stop_mpi,pname,error_message)
    !
  end do part_loop
  !
  ! Deallocate all the local arrays before leaving
  !
  if (allocated(vwgt)) then
    deallocate ( vwgt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vwgt",2,ierr,error_message)
  end if
  if (allocated(vsize)) then
    deallocate ( vsize , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vsize",2,ierr,error_message)
  end if
  if (allocated(adjwgt)) then
    deallocate ( adjwgt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"adjwgt",2,ierr,error_message)
  end if
  if (allocated(tpwgts)) then
    deallocate ( tpwgts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"tpwgts",2,ierr,error_message)
  end if
  if (allocated(ubvec)) then
    deallocate ( ubvec , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ubvec",2,ierr,error_message)
  end if
  !
  ! Format statements
  !
  1 format (/,'-------------------------------------------------------',//, &
              '         Used METIS 5.1 to partition the graph         ',/, &
              '    using a multilevel ',a,' algorithm    ',/, &
              '           into ',i0,' ',a,' partitions.          ',/, &
              '     The minimized ',a,' was ',i0,'.    ',//, &
              '-------------------------------------------------------',/)
  2 format (" None of the partitioning methods available in METIS 5.1 were", &
            " able to find a quality partitioning of the grid!")
  3 format (" The METIS partitioning method requested in the input file is", &
            " invalid!")
  4 format (///,'########################################################',/, &
                '########## !! PARTITIONING ATTEMPT FAILED !! ###########',/, &
                '########################################################',//, &
                '  METIS 5.1 failed to successfully partition the graph  ',/, &
                '    using a multilevel ',a,' algorithm    ',/, &
                '       into ',i0,' ',a,' partitions',/, &
                '   while minimizing the ',a,'!',/, &
                '         Switching to the next METIS algorithm         ',/, &
                '               to partition the graph...               ',//, &
                '########################################################',/, &
                '########## !! PARTITIONING ATTEMPT FAILED !! ###########',/, &
                '########################################################',///)
  !
#endif
end subroutine partition_using_metis_5
!
!###############################################################################
!
subroutine create_dual_graph(cell_geom,cell_order,nodes_of_cell, &
                             nodes_of_cell_ptr,bface,xadj,adjncy)
  !
  !.. Use Statements ..
  use order_mod,        only : geom_solpts
  use connectivity_mod, only : GetFaceNodes_for_QuadFace
  use connectivity_mod, only : GetFaceNodes_for_TriaFace
  use connectivity_mod, only : RotateFluxPoints_Quad
  use connectivity_mod, only : RotateFluxPoints_Tria
  !
  !.. Formal Arguments ..
  integer,              dimension(:),      intent(in) :: cell_geom
  integer,              dimension(:),      intent(in) :: cell_order
  integer,              dimension(:),      intent(in) :: nodes_of_cell
  integer,              dimension(:),      intent(in) :: nodes_of_cell_ptr
  integer,              dimension(:,:),    intent(in) :: bface
  integer, allocatable, dimension(:),   intent(inout) :: xadj
  integer, allocatable, dimension(:),   intent(inout) :: adjncy
  !
  !.. Local Scalars ..
  integer :: ic,jc,n,nc,nf,ncf,nfp,nfnods,ifa,ieadj,ierr
  integer :: i1,i2,j1,j2,l1,l2,n1,n2,r1,r2,c1,c2,csc_ptr
  integer :: this_cell,this_side,this_geom,this_order
  integer :: ghost_cell,adj_cell,adj_geom,adj_order
  integer :: host_cell,host_side,host_geom,host_order
  integer :: left_cell,left_side,left_geom,left_order
  integer :: right_cell,right_side,right_geom,right_order
  integer :: face_geom,face_order,face_counter
  logical :: error_finding_rotation
  character(len=200) :: array_name
  !
  !.. Local Arrays ..
  integer, dimension(1:4) :: ip
  integer, dimension(1:4) :: face_nodes
  integer, dimension(1:4) :: left_nodes
  integer, dimension(1:4) :: right_nodes
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:)   :: lpoin
  integer, allocatable, dimension(:)   :: cells_surr_cell
  integer, allocatable, dimension(:)   :: cells_surr_cell_ptr
  integer, allocatable, dimension(:)   :: cells_with_node
  integer, allocatable, dimension(:)   :: cells_with_node_ptr
  integer, allocatable, dimension(:,:) :: faces_of_cell
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_dual_graph"
  !
continue
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Find the cells that surround each node
  !
  allocate ( cells_with_node_ptr(1:ngnode+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_with_node_ptr",1,ierr,error_message)
  !
  ! Loop over the cells, storing 'ahead'
  !
  do n = 1,ngcell
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    do ic = n1,n2
      i1 = nodes_of_cell(ic) + 1
      cells_with_node_ptr(i1) = cells_with_node_ptr(i1) + 1
    end do
  end do
  !
  ! Reshuffle cells_with_node_ptr
  !
  do n = 2,size(cells_with_node_ptr)
    cells_with_node_ptr(n) = cells_with_node_ptr(n) + cells_with_node_ptr(n-1)
  end do
  !
  allocate ( cells_with_node(1:last(cells_with_node_ptr)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_with_node",1,ierr,error_message)
  !
  ! Now store the cells surrounding each node in cells_with_node
  !
  do n = 1,ngcell
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    do ic = n1,n2
      i1 = nodes_of_cell(ic)
      i2 = cells_with_node_ptr(i1) + 1
      cells_with_node_ptr(i1) = i2
      cells_with_node(i2) = n
    end do
  end do
  !
  ! Finally, reorder cells_with_node_ptr
  !
  do n = size(cells_with_node_ptr),2,-1
    cells_with_node_ptr(n) = cells_with_node_ptr(n-1)
  end do
  cells_with_node_ptr(1) = 0
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Find the cells surrounding each cell which
  ! defines the dual graph of the nodal grid 
  !
  allocate ( cells_surr_cell_ptr(1:ngcell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_surr_cell_ptr",1,ierr,error_message)
  !
  ! Loop to set up cells_surr_cell_ptr
  !
  do n = 1,ngcell
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    cells_surr_cell_ptr(n+1) = cells_surr_cell_ptr(n) + n2-n1+1
  end do
 !!
 !! Reshuffle cells_surr_cell_ptr
 !!
 !do n = 2,size(cells_surr_cell_ptr)
 !  cells_surr_cell_ptr(n) = cells_surr_cell_ptr(n) + cells_surr_cell_ptr(n-1)
 !end do
  !
  allocate ( cells_surr_cell(1:last(cells_surr_cell_ptr)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_surr_cell",1,ierr,error_message)
  !
  ! Loop to find cells_surr_cells
  !
  allocate ( lpoin(1:ngnode) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"lpoin",1,ierr,error_message)
  !
  do n = 1,ngcell
    !
    ! Get the surrounding cells based on the geometry of the current cell
    !
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    !
    ! Get the pointer index of cells_surr_cell for this cell
    !
    csc_ptr = cells_surr_cell_ptr(n)
    !
    ! Geometry of the current cell
    !
    this_geom = cell_geom(n)
    !
    ! Loop over the number of cell faces for the current geometry
    !
    do this_side = 1,geom_faces(this_geom)
      !
      ! Find the number of nodes that need to be matched for this face
      !
      nfnods = num_face_nodes(this_side,this_geom)
      !
      ! Node offsets for this face of the current cell
      !
      ip(1:nfnods) = side_nodes(1:nfnods,this_side,this_geom)
      !
      ! Store nodes of the current face for the current cell
      !
      face_nodes(1:nfnods) = nodes_of_cell( n1 + ip(1:nfnods) )
      !
      ! Mark each of these nodes on the current face
      !
      lpoin ( face_nodes(1:nfnods) ) = 1
      !
      ! Get the beginning and ending indices for the
      ! cells surrounding the first node on the face
      !
      i1 = cells_with_node_ptr( face_nodes(1) )+1
      i2 = cells_with_node_ptr( face_nodes(1) +1)
      !
      ! Initialize the adjacent face marker to zero
      !
      ieadj = 0
      !
      ! Loop over the cells that contain the first node on this face
      !
      find_cell: do ic = i1,i2
        !
        ! Current cell containing this node
        !
        jc = cells_with_node(ic)
        !
        ! Skip to the next cell if this is the current cell we are working on
        !
        if (jc == n) cycle
        !
        ! Get the beginning and ending indices of the nodes defining the
        ! possible adjacent cell
        !
        j1 = nodes_of_cell_ptr(jc)+1
        j2 = nodes_of_cell_ptr(jc+1)
        !
        ! If summing the lpoin array using the current nodes of the possible
        ! adjacent cell is equal to the number of nodes we need to match, then
        ! this is the adjacent cell to the current face.
        !
        if (sum(lpoin(nodes_of_cell(j1:j2))) == nfnods) then
          ieadj = jc
          exit find_cell
        end if
        !
      end do find_cell
      !
      ! Store the adjacent cell index
      !
      cells_surr_cell(csc_ptr+this_side) = ieadj
      !
      ! Re-mark the nodes in lpoin to zero
      !
      lpoin( face_nodes(1:nfnods) ) = 0
      !
    end do
    !
  end do
  !
  ! We are finished with cells_with_node and cells_with_node_ptr so deallocate
  !
  deallocate ( cells_with_node , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_with_node",2,ierr,error_message)
  deallocate ( cells_with_node_ptr , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_with_node_ptr",2,ierr,error_message)
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! We cant just use cells_surr_cell_ptr and cells_surr_cell as xadj
  ! and adjncy because cells_surr_cell is zero for faces of cells
  ! that correspond to boundary faces. In the loop "find_cell" above,
  ! if the sum(lpoin(nodes_of_cell(j1:j2))) term never equals 2 then
  ! this loop exits with ieadj still equal to 0. We cannot send adjncy
  ! containing zeros to the METIS routines so these need to be removed
  ! and the locations in xadj also have to be adjusted accordingly.
  ! Additionally, have xadj to start the indexing for adjncy at 1,
  ! i.e., xadj(1)=1 (cells_surr_cell_ptr indexing for cells_surr_cell
  ! begins at 0, i.e. cells_surr_cell_ptr(1)=0) since this is what
  ! METIS expects when using Fortran array ordering.
  !
 !do n = 1,ngcell
 !  n1 = cells_surr_cell_ptr(n)+1
 !  n2 = cells_surr_cell_ptr(n+1)
 !  write (101,'(4(3x,i8))') cells_surr_cell(n1:n2)
 !end do
  !
  allocate ( xadj(1:ngcell+1) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xadj",1,ierr,error_message)
  xadj(1) = 1
  !
  n = count( cells_surr_cell(:) /= 0 )
  allocate ( adjncy(1:n) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"adjncy",1,ierr,error_message)
  !
  ! Loop through the cells to copy only the non-zero
  ! locations of cells_surr_cell into the adjncy array and
  ! keep count of the number copied for each cell in xadj
  !
  ic = 0
  do n = 1,ngcell
    n1 = cells_surr_cell_ptr(n)+1
    n2 = cells_surr_cell_ptr(n+1)
    do nf = n1,n2
      if (cells_surr_cell(nf) /= 0) then
        xadj(n+1) = xadj(n+1) + 1
        ic = ic + 1
        adjncy(ic) = cells_surr_cell(nf)
      end if
    end do
  end do
  !
  ! Finally, reshuffle xadj
  !
  do n = 2,size(xadj)
    xadj(n) = xadj(n) + xadj(n-1)
  end do
  !
 !do n = 1,ngcell
 !  n1 = xadj(n)
 !  n2 = xadj(n+1)-1
 !  write (102,'(4(3x,i8))') adjncy(n1:n2)
 !end do
  !
  ! #######
  ! ## 4 ##
  ! #######
  !
  ! Find the cells on each face. Once the grid is partitioned,
  ! this will be needed to adjust the boundary faces to include
  ! former interior faces that now separate adjacent partitions
  !
  ! Find the total number of faces and then
  ! allocate the face and faces_of_cell arrays
  !
  ngface = ngfbnd
  do n = 1,ngcell
    n1 = cells_surr_cell_ptr(n)+1
    n2 = cells_surr_cell_ptr(n+1)
    ngface = ngface + count( cells_surr_cell(n1:n2) > n )
  end do
  !
  ! Create the face and faces_of_cell arrays
  !
  allocate ( face(1:ngface) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face",1,ierr,error_message)
  !
  allocate ( faces_of_cell(1:6,1:ngcell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faces_of_cell",1,ierr,error_message)
  !
  ! ##########
  ! ##      ##
  ! ## 4(A) ##
  ! ##      ##
  ! ##########
  !
  ! Get the cells on each boundary face
  !
  do nf = 1,ngfbnd
    !
    host_cell = bface(2,nf)   ! Host cell on boundary face
    ghost_cell = bface(10,nf) ! Ghost cell on boundary face
    !
    host_geom = cell_geom(host_cell)             ! geometry of host cell
    host_order = cell_order(host_cell)           ! order of host cell
    host_side = bface(4,nf)                      ! face of host cell on boundary
    nfnods = num_face_nodes(host_side,host_geom) ! # of nodes defining bnd face
    !
    ! beginning index of nodes defining host cell
    n1 = nodes_of_cell_ptr(host_cell)+1
    ! ending    index of nodes defining host cell
    n2 = nodes_of_cell_ptr(host_cell+1)
    !
    ! Node offsets for this face
    ip(1:nfnods) = side_nodes(1:nfnods,host_side,host_geom)
    !
    ! Store the cells on this boundary face
    !
    face(nf)%left%cell = host_cell
    face(nf)%left%cell_face = host_side
    face(nf)%left%geom = host_geom
    face(nf)%left%order = host_order
    !
    face(nf)%right%cell = ghost_cell
    face(nf)%right%cell_face = 1
    face(nf)%right%geom = geom_of_face(nfnods)
    face(nf)%right%order = host_order
    !
    face(nf)%geom = geom_of_face(nfnods)
    face(nf)%order = host_order
    !
    ! Store the nodes on this boundary face
    !
    face(nf)%nodes(1:nfnods) = nodes_of_cell( n1 + ip(1:nfnods) )
    !
    ! Store this face for the host cell
    !
    faces_of_cell(host_side,host_cell) = nf
    !
  end do
  !
  ! ##########
  ! ##      ##
  ! ## 4(B) ##
  ! ##      ##
  ! ##########
  !
  ! Get the cells for the interior faces
  !
  face_counter = ngfbnd
  do this_cell = 1,ngcell
    !
    c1 = cells_surr_cell_ptr(this_cell)+1
    c2 = cells_surr_cell_ptr(this_cell+1)
    !
    n1 = nodes_of_cell_ptr(this_cell)+1
    n2 = nodes_of_cell_ptr(this_cell+1)
    !
    this_geom = cell_geom(this_cell)
    this_order = cell_order(this_cell)
    !
    do nc = c1,c2
      !
      adj_cell = cells_surr_cell(nc)
      !
      if (adj_cell>this_cell .and. adj_cell<=ngcell) then
        !
        adj_geom = cell_geom(adj_cell)   ! geometry of adjacent cell
        adj_order = cell_order(adj_cell) ! order    of adjacent cell
        !
        face_counter = face_counter + 1
        !
        this_side = nc-c1+1
        !
        ! # of nodes defining this face
        nfnods = num_face_nodes(this_side,this_geom)
        ! node offsets for this face
        ip(1:nfnods) = side_nodes(1:nfnods,this_side,this_geom)
        !
        ! Store the cells on this interior face
        !
        face(face_counter)%left%cell = this_cell
        face(face_counter)%right%cell = adj_cell
        !
        face(face_counter)%left%geom = this_geom
        face(face_counter)%right%geom = adj_geom
        !
        face(face_counter)%left%order = this_order
        face(face_counter)%right%order = adj_order
        !
        ! Store the nodes on this interior face
        !
        face(face_counter)%nodes(1:nfnods) = nodes_of_cell( n1 + ip(1:nfnods) )
        !
        ! Store the geometry of this face
        !
        face(face_counter)%geom = geom_of_face(nfnods)
        face(face_counter)%order = max( this_order, adj_order )
        !
        ! Store the side of the left cell that is on this face since we have
        ! this information available right now. We will figure out the side of
        ! the right cell in section 4(C).
        !
        face(face_counter)%left%cell_face = this_side
        !
        ! Store this face for the left cell
        !
        faces_of_cell(this_side,this_cell) = face_counter
        !
      end if
      !
    end do
    !
  end do
  !
  ! We no longer need cells_surr_cell and cells_surr_cell_ptr so deallocate
  !
  deallocate ( cells_surr_cell , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_surr_cell",2,ierr,error_message)
  deallocate ( cells_surr_cell_ptr , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_surr_cell_ptr",2,ierr,error_message)
  !
  ! ##########
  ! ##      ##
  ! ## 4(C) ##
  ! ##      ##
  ! ##########
  !
  ! Find the side of the right cell that is on each interior face
  !
  lpoin = 0
  !
  loop_right_cell_face: do nf = ngfbnd+1,ngface
    !
    nfnods = geom_nodes(face(nf)%geom) ! # of nodes defining this face
    !
    ! Mark the nodes on the current face
    !
    lpoin( face(nf)%nodes(1:nfnods) ) = 1
    !
    ! Information for the right cell on the face
    !
    right_cell = face(nf)%right%cell
    right_geom = face(nf)%right%geom
    !
    n1 = nodes_of_cell_ptr(right_cell)+1
    n2 = nodes_of_cell_ptr(right_cell+1)
    !
    loop_cell_faces: do ncf = 1,geom_faces(right_geom)
      !
      ! Skip this cell face if the number of face nodes doesnt match
      !
      if (num_face_nodes(ncf,right_geom) /= nfnods) cycle loop_cell_faces
      !
      ! Store the nodes of the cell face in question
      !
      ip(1:nfnods) = nodes_of_cell(n1+side_nodes(1:nfnods,ncf,right_geom))
      !
      ! If this face subset of lpoin sums up to nfnods, this is the side
      ! of the right cell on the current face. Store the value, reset all
      ! of lpoin to zero, and skip to the next face.
      !
      if (sum(lpoin(ip(1:nfnods))) == nfnods) then
        face(nf)%right%cell_face = ncf
        faces_of_cell(ncf,right_cell) = nf
        lpoin( ip(1:nfnods) ) = 0
        cycle loop_right_cell_face
      end if
      !
    end do loop_cell_faces
    !
  end do loop_right_cell_face
  !
  ! We dont need lpoin anymore so deallocate
  !
  deallocate ( lpoin , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"lpoin",2,ierr,error_message)
  !
  ! ##########
  ! ##      ##
  ! ## 4(D) ##
  ! ##      ##
  ! ##########
  !
  !
  ! For each cell on the face, get the indices for the face/flux points of the
  ! cell that are on this face. These indices are used in the correct and
  ! interp(:)%toFace matrices.
  !
  error_finding_rotation = fals
  !
  do nf = 1,ngface
    !
    left_cell = face(nf)%left%cell
    right_cell = face(nf)%right%cell
    !
    left_side = face(nf)%left%cell_face
    right_side = face(nf)%right%cell_face
    !
    left_geom = face(nf)%left%geom
    right_geom = face(nf)%right%geom
    !
    left_order = face(nf)%left%order
    right_order = face(nf)%right%order
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    nfp = geom_solpts(face_geom,face_order)
    !
    nfnods = geom_nodes(face_geom)
    !
    l1 = nodes_of_cell_ptr(left_cell)+1
    l2 = nodes_of_cell_ptr(left_cell+1)
    !
    r1 = nodes_of_cell_ptr(right_cell)+1
    r2 = nodes_of_cell_ptr(right_cell+1)
    !
    ip(1:nfnods) = side_nodes(1:nfnods,left_side,left_geom)
    left_nodes(1:nfnods) = nodes_of_cell(l1 + ip(1:nfnods))
    !
    ! #############################
    ! ### LEFT CELL ON THE FACE ###
    ! #############################
    !
    ! Allocate the flxpts component of the face array for the left cell
    !
    allocate ( face(nf)%left%flxpts(1:nfp) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,2) nf,"left"
    call alloc_error(pname,array_name,1,ierr,error_message,skip_alloc_pause)
    !
    ! The indices of the flux points for the left cell are just the 
    ! sequential values of [1:nfp]
    !
    face(nf)%left%flxpts = intseq(0,nfp-1)!+ 1
    !
    ! ##############################
    ! ### RIGHT CELL ON THE FACE ###
    ! ##############################
    !
    ! For the right cell, it depends on whether
    ! this is a boundary or interior face
    !
    ! Allocate the flxpts component of the face array for the left cell
    !
    allocate ( face(nf)%right%flxpts(1:nfp) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,2) nf,"right"
    call alloc_error(pname,array_name,1,ierr,error_message,skip_alloc_pause)
    !
    if (nf <= ngfbnd) then
      !
      ! This is a boundary face so get the boundary nodes
      !
      ip(1:nfnods) = intseq(0,nfnods-1)
      right_nodes(1:nfnods) = nodes_of_cell(r1 + ip(1:nfnods))
      !
    else
      !
      ! This is an interior face so get the nodes
      ! for this side of the right cell
      !
      ip(1:nfnods) = side_nodes(1:nfnods,right_side,right_geom)
      right_nodes(1:nfnods) = nodes_of_cell(r1 + ip(1:nfnods))
      !
    end if
    !
    ! Get the flux point rotation for the right cell
    !
    if (face_geom == Geom_Edge) then
      !
      ! If this is a 1D edge between two 2D cells, the ordering of the flux
      ! points for the right cell is just the reverse ordering of the flux
      ! points for the left cell.
      !
      face(nf)%right%flxpts = intseq(nfp-1,0,-1)!+ 1
      !
    else if (face_geom == Geom_Quad) then
      !
      ! Get the nodes from the left and right cells that are on the current face
      ! and in the correct ABCD ordering needed for the RotateFluxPoints_Quad
      ! function.
      !
      call GetFaceNodes_for_QuadFace(nodes_of_cell(l1:l2),left_side, &
                                     nodes_of_cell(r1:r2),right_side, &
                                     left_nodes(1:nfnods),right_nodes(1:nfnods))
      !
      ! Find the face rotation between the two cells on the interface and use
      ! this rotation to get the correct ordering of the flux points for the
      ! right cell relative to the flux points of the left cell.
      !
      face(nf)%right%flxpts = RotateFluxPoints_Quad(left_nodes(1:nfnods), &
                                                    right_nodes(1:nfnods), &
                                                    face_order)!+ 1
      !
      ! If the previous function made no changes, mark that there was an error
      ! finding the correct face rotation between cells.
      !
      if (all(face(nf)%right%flxpts == 1)) error_finding_rotation = true
      !
    else if (face_geom == Geom_Tria) then
      !
      call GetFaceNodes_for_TriaFace(nodes_of_cell(l1:l2),left_side, &
                                     nodes_of_cell(r1:r2),right_side, &
                                     left_nodes(1:nfnods),right_nodes(1:nfnods))
      !
      face(nf)%right%flxpts = RotateFluxPoints_Tria(left_nodes(1:nfnods), &
                                                    right_nodes(1:nfnods), &
                                                    nfp)!+ 1
      !
      if (all(face(nf)%right%flxpts == 1)) error_finding_rotation = true
      !
    end if
    !
  end do
  !
  ! Deallocate faces_of_cell before we leave
  !
  deallocate ( faces_of_cell , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faces_of_cell",2,ierr,error_message)
  !
  ! Report and stop execution if there was an error finding
  ! the rotation of the flux points on any of the faces.
  !
  if (error_finding_rotation) then
    write (iout,*) "There was an error rotating flux points!"
    call stop_recon(stop_mpi,pname)
  end if
  !
  ! Format Statements
  !
  1 format (" Error on face - ",i0,",  nodes ",i0,3(:,2x,i0))
  2 format ("face(",i0,")%",a,"%flxpts")
  !
end subroutine create_dual_graph
!
!###############################################################################
!
subroutine create_grid_weights(metis_option,vwgt,adjwgt,cell_geom,cell_order)
  !
  !.. Formal Arguments ..
  integer,                  intent(in) :: metis_option
  integer, dimension(:),    intent(in) :: cell_geom
  integer, dimension(:),    intent(in) :: cell_order
  integer, dimension(:), intent(inout) :: vwgt
  integer, dimension(:), intent(inout) :: adjwgt
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_grid_weights"
  !
continue
  !
  ! Initialize vwgt
  ! NOTE: vwgt gives the computational weight of each cell
  !
  vwgt = 0
  !
  !############################################################
  !############################################################
  ! We need to figure out the computational difference between
  ! triangle and quad cells and populate vwgt accordingly 
  !############################################################
  !############################################################
  !
  ! Initialize adjwgt
  !
  if ( any(metis_option == (/1,2/)) ) then
    !
    ! If minimizing the edgecut, adjwgt gives the weights of the faces
    !
    !###############################################################
    !###############################################################
    ! Not sure if we need any weighting on the faces so set to zero
    adjwgt = 0
    !###############################################################
    !###############################################################
    !
  else if ( any(metis_option == (/3,4/)) ) then
    !
    ! If minimizing the communication volume, adjwgt gives the communication
    ! weight of the cells (i.e., how much data has to be sent to/from this cell)
    !
    ! Set the communication weight of the cells to
    ! the number of solution points for each cell
    !
    adjwgt = Cell_Solpts( cell_geom , cell_order )
    !
  end if
  !
end subroutine create_grid_weights
!
!###############################################################################
!
subroutine check_quality_of_partitioning(npart,epart, &
                                         passes_inspection)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: npart
  integer, dimension(:), intent(in) :: epart
  logical, intent(out) :: passes_inspection
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_quality_of_partitioning"
  !
continue
  !
  passes_inspection = true
  !
  ! Check that all processors got part of the global grid
  !
  do n = 1,npart
    if (any(epart == n)) cycle
    write (iout,100) n
    passes_inspection = fals
   !call stop_recon(abort,pname)
  end do
  !
  !############################################################################
  !############################################################################
  !############################################################################
  !#####   NEED TO ADD CODE TO CHECK THAT ALL PARTITIONS ARE CONTIGUOUS   #####
  !############################################################################
  !############################################################################
  !############################################################################
  !
  ! Format statements
  !
  100 format(/,' Error: Processor ',i0,' was excluded from the          ',/, &
               ' partitioning of the elements. More than likely too many',/, &
               ' processors were used in the paritioning of this mesh.  ',/, &
               ' Please decrease number of processors and try again.    ',/)
  !
end subroutine check_quality_of_partitioning
!
!###############################################################################
!
subroutine find_partition_boundaries(epart,bface)
  !
  !.. Formal Arguments ..
  integer,              dimension(:),      intent(in) :: epart
  integer, allocatable, dimension(:,:), intent(inout) :: bface
  !
  !.. Local Scalars ..
  integer :: new_ngfbnd,new_ngface,ierr
  integer :: n,nf,n1,n2,np
  integer :: i,ic,l1,l2,lp,ifa,pf
  integer :: j,jc,r1,r2,rp,jfa,pc
  integer :: host_cell,host_geom,host_side
  integer :: ghost_cell,nfnods,switch_index
  integer :: left_cell,left_geom,left_side
  integer :: right_cell,right_geom,right_side
  !
  !.. Local Arrays ..
  integer, dimension(1:4) :: ip
  !
  !.. Local Allocatable Arrays ..
  logical,      allocatable, dimension(:) :: msk
  type(face_t), allocatable, dimension(:) :: temp_face
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_partition_boundaries"
  !
continue
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  allocate ( msk(1:ngface) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",1,ierr,error_message)
  !
  ! Find the interior faces that divide two partitions
  !
  do nf = ngfbnd+1,ngface
    msk(nf) = epart( face(nf)%left%cell ) /= epart( face(nf)%right%cell )
  end do
  !
  ! Copy the partition information into the face array. For now, dont
  ! do anything about the partition information for the ghost cells
  ! on the boundary faces.
  !
  do nf = 1,ngface
    face(nf)%left%partition = epart(face(nf)%left%cell)
    if (nf > ngfbnd) then
      face(nf)%right%partition = epart(face(nf)%right%cell)
    end if
  end do
  !
  ! Recompute the number of boundary faces and the total number of faces.
  !   Each internal face dividing unstructured blocks will be turned
  !   into two boundary faces, one for each side of the face meaning:
  !   1. the number of internal faces will shrink by count(msk)
  !   2. the number of boundary faces will increase by 2*count(msk)
  !   3. the total number of faces will increase by count(msk)
  !
  new_ngfbnd = ngfbnd + 2*count(msk)
  new_ngface = ngface + 1*count(msk)
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Reallocate bface to make room for new partition boundary faces
  !
  call reallocate(bface,new_size2=new_ngfbnd)
  !
  ! Replace bface(10,nf) with the location within the face array to extract data
  !
  do nf = 1,ngfbnd
    bface(10,nf) = nf
  end do
  !
  ! Collect and add the new boundary faces to bface
  !
  nf = ngfbnd + 1
  do n = ngfbnd+1,ngface
    !
    if (.not. msk(n)) cycle
    !
    left_cell = face(n)%left%cell
    left_side = face(n)%left%cell_face
    left_geom = face(n)%left%geom
    !
    right_cell = face(n)%right%cell
    right_side = face(n)%right%cell_face
    right_geom = face(n)%right%geom
    !
    ! Now add this information the boundary face array
    !
    ! Left side of the face
    !
    bface( 1,nf) = bc_cpu_bnd  ! boundary condition
    bface( 2,nf) = left_cell   ! host cell
    bface( 3,nf) = left_geom   ! host cell geometry type
    bface( 4,nf) = left_side   ! face of host cell on boundary
    bface( 5,nf) = nf + 1      ! adjacent partition boundary face
    bface( 6,nf) = right_cell  ! adjacent partition host cell
    bface( 7,nf) = epart(right_cell) - 1  ! adjacent partition processor ID
    bface( 8,nf) = right_geom  ! adjacent partition host cell geometry type
    bface( 9,nf) = right_side  ! face of adjacent partition host cell on boundary
    bface(10,nf) = n           ! location within face array to extract data
    !
    nf = nf + 1
    !
    ! Right side of the face
    !
    bface( 1,nf) = bc_cpu_bnd  ! boundary condition
    bface( 2,nf) = right_cell  ! host cell
    bface( 3,nf) = right_geom  ! host cell geometry type
    bface( 4,nf) = right_side  ! face of host cell on boundary
    bface( 5,nf) = nf - 1      ! adjacent partition boundary face
    bface( 6,nf) = left_cell   ! adjacent partition host cell
    bface( 7,nf) = epart(left_cell) - 1  ! adjacent partition processor ID
    bface( 8,nf) = left_geom   ! adjacent partition host cell geometry type
    bface( 9,nf) = left_side   ! face of adjacent partition host cell on boundary
    bface(10,nf) = -n          ! location within face array to extract data
    !                          ! NOTE: This is negative to mark that this host
    nf = nf + 1                !       cell used to be a right cell on the face
    !
  end do
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Allocate a temporary array to resize and reorder face
  !
  allocate ( temp_face(1:new_ngface) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_face",1,ierr,error_message)
  !
  ! Redo the boundary face section of face
  ! because of the new additions to bface
  !
  do nf = 1,new_ngfbnd
    !
    ! Copy everything over from the old location in the face array
    !
    temp_face(nf) = face(abs(bface(10,nf)))
    !
    if (bface(1,nf) == bc_cpu_bnd) then
      !
      ! If bface(10,nf) is negative, that means the host cell for this boundary
      ! face used to be the right cell when this partition boundary previously
      ! was an interior face. Therefore, we need to flip the left and right
      ! cells on the boundary face in order to make the right cell the host
      ! cell for this boundary face. We also need to reverse the 
      ! rotation/ordering of the face nodes so that the nodes of the
      ! ghost cell are defined such that the normal of the face points
      ! outside the grid domain.
      !
      if (bface(10,nf) < 0) then
        !
        ! Switch the left and right cell information for this face
        !
        temp_face(nf)%left = face(abs(bface(10,nf)))%right
        temp_face(nf)%right = face(abs(bface(10,nf)))%left
        !
        ! Switch the order of the face nodes so that the face normal points
        ! outside the grid domain. See subroutine adjust_for_ghost_cells for
        ! further details on whats being done here.
        !
        nfnods = geom_nodes(temp_face(nf)%geom) ! # of nodes defining this face
        !
        ! Get the index of the node to switch with node 1
        !
        switch_index = min( nfnods , 3 )
        !
        ! Switch the two face nodes required to reverse the face node ordering
        !
        n = temp_face(nf)%nodes(1)
        temp_face(nf)%nodes(1) = temp_face(nf)%nodes(switch_index)
        temp_face(nf)%nodes(switch_index) = n
        !
      end if
      !
    end if
    !
    ! Reset bface(10,nf) back to the index for the ghost cell
    !
    bface(10,nf) = ngcell + nf
    !
  end do
  !
  ! Copy the interior face information to the temporary array
  ! but skip those interior faces that are now boundary faces
  !
  nf = new_ngfbnd + 1
  do n = ngfbnd+1,ngface
    if ( msk(n) ) cycle
    temp_face(nf) = face(n)
    nf = nf + 1
  end do
  !
  ! Use move_alloc to reallocate face. This does the following
  !     1. Deallocates face
  !     2. Allocates face with the same type, type
  !        parameters, array bounds, and value as temp_face
  !     3. Deallocates temp_face
  !
  call move_alloc( from=temp_face , to=face )
  !
  ! We need to update the global grid parameters
  !
  ngfbnd = new_ngfbnd
  ngface = new_ngface
  !
  ! Deallocate msk since we are finished with it
  !
  deallocate ( msk , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",2,ierr,error_message)
  !
  ! #######
  ! ## 4 ##
  ! #######
  ! 
  ! Since we have all the information here, we need to make an addition to
  ! the bface array for periodic boundary faces. For the periodic boundary
  ! face 'nf' (with host cell 'nc') which is linked to another periodic
  ! boundary face 'pf' (with host cell 'pc'), we need to store in
  ! bface(7,nf) the processor that contains the adjacent partition. We
  ! will also set the processor value in bface(7,nf) to the current
  ! processor ID for all other non-communication boundaries.
  !
  do nf = 1,ngfbnd
    if (bface(1,nf) == bc_periodic) then
      pf = bface(5,nf)
      pc = bface(2,pf)
      bface(7,nf) = epart(pc) - 1  ! adjacent partition processor ID
    else if (bface(1,nf) /= bc_cpu_bnd) then
      bface(7,nf) = mypnum
     !bface(7,nf) = epart(bface(2,nf))-1
    end if
  end do
  !
  ! Now copy the partition information for boundary faces into the face array
  !
  do nf = 1,ngfbnd
    if (bface(1,nf) /= bc_cpu_bnd) then
      face(nf)%right%partition = bface(7,nf)
    end if
  end do
  !
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  ! #####                                                                #####
  ! #####  AT SOME POINT WE WILL NEED TO MAKE SOME ALTERATIONS TO THE    #####
  ! #####  FLXPTS, CELL_FACE, AND OTHER COMPONENTS OF THE FACE(:)        #####
  ! #####  ARRAY FOR COMMUNICATION BOUNDARY FACES.                       #####
  ! #####                                                                #####
  ! #####                  THIS IS TO MAKE SURE THAT                     #####
  ! #####                                                                #####
  ! #####  THE HOST CELL FLUX POINTS AND BOUNDARY FACE FLUX POINTS ARE   #####
  ! #####  CONSISTENT LOCALLY ON THEIR PARTITION                         #####
  ! #####                                                                #####
  ! #####                       AS WELL AS                               #####
  ! #####                                                                #####
  ! #####  THE HOST CELL FLUX POINT DATA AND THE ADJACENT COMMUNICATION  #####
  ! #####  CELL FLUX POINT DATA IS CONSISTENT DURING THE COMMUNICATION   #####
  ! #####  SENDS AND RECEIVES BETWEEN THE TWO PARTITIONS.                #####
  ! #####                                                                #####
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  !
  ! Perform one last quick test to make sure that the attached cells
  ! for each interior face are located within the same partition
  !
  if ( any( face(ngfbnd+1:ngface)%left%partition /= &
            face(ngfbnd+1:ngface)%right%partition ) ) then
    write (iout,1)
    call stop_recon(stop_mpi,pname)
  end if
  !
  ! Format Statements
  !
  1 format (//,"Error: Found two cells attached to an interior face",/, &
               "       that are not located within the same partition!",//)
  !
end subroutine find_partition_boundaries
!
!###############################################################################
!
subroutine create_grid_mappings(npart,epart, &
                                nodes_of_cell,nodes_of_cell_ptr, &
                                node_map,face_map,cell_map)
  !
  !.. Formal Arguments ..
  integer,                                   intent(in) :: npart
  integer,                  dimension(:),    intent(in) :: epart
  integer,                  dimension(:),    intent(in) :: nodes_of_cell
  integer,                  dimension(:),    intent(in) :: nodes_of_cell_ptr
  type(map_t), allocatable, dimension(:), intent(inout) :: node_map
  type(map_t), allocatable, dimension(:), intent(inout) :: face_map
  type(map_t), allocatable, dimension(:), intent(inout) :: cell_map
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,ncount,ierr
  character(len=200) :: array_name 
  !
  !.. Local Allocatable Arrays ..
  logical, allocatable, dimension(:) :: msk
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_grid_mappings"
  !
continue
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Create the cell mappings
  !
  allocate ( cell_map(1:npart) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_map",1,ierr,error_message)
  !
  ! Loop through the partitions to find the
  ! local number of cells for each partition
  !
  do n = 1,npart
    !
    ! Count the cells on this partition
    !
    ncount = count( epart == n )
    !
    ! Allocate the loc_to_glb component of cell_map for this partition
    !
    allocate ( cell_map(n)%loc_to_glb(1:ncount) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "cell_map",n
    call alloc_error(pname,array_name,1,ierr,error_message)
    !
    ! Add the cells of the current partition to cell_map
    !
    cell_map(n)%loc_to_glb = pack( intseq(1,ngcell) , epart(1:ngcell) == n )
   !cell_map(n)%loc_to_glb = find_subset(epart,n)
    !
  end do
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Create the node mappings
  !
  allocate ( node_map(1:npart) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"node_map",1,ierr,error_message)
  !
  ! Loop through the partitions to find the
  ! local number of nodes for each partition
  !
  allocate ( msk(1:ngnode) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",1,ierr,error_message)
  !
  do n = 1,npart
    !
    ! Initialize/Reset msk
    !
    msk(:) = fals
    !
    ! Mark the nodes that are within this partition
    !
    do nc = 1,ngcell
      if (epart(nc) /= n) cycle
      n1 = nodes_of_cell_ptr(nc)+1
      n2 = nodes_of_cell_ptr(nc+1)
      msk(nodes_of_cell(n1:n2)) = true
    end do
    !
    ! Count the nodes on this partition
    !
    ncount = count( msk(1:ngnode) )
    !
    ! Allocate the loc_to_glb component of node_map for this partition
    !
    allocate ( node_map(n)%loc_to_glb(1:ncount) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "node_map",n
    call alloc_error(pname,array_name,1,ierr,error_message)
    !
    ! Add the nodes of the current partition to node_map
    !
    node_map(n)%loc_to_glb = pack( intseq(1,ngnode) , msk(1:ngnode) )
    !
  end do
  !
  ! Deallocate msk because it is no longer needed
  !
  deallocate ( msk , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",2,ierr,error_message)
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Create the face mappings
  !
  allocate ( face_map(1:npart) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_map",1,ierr,error_message)
  !
  ! Loop through the partitions to find the
  ! local number of faces for each partition
  !
  do n = 1,npart
    !
    ! Count the faces that are within this partition
    !
    ncount = count( face(1:ngface)%left%partition == n )
    !
    ! Allocate the loc_to_glb component of face_map for this partition
    !
    allocate ( face_map(n)%loc_to_glb(1:ncount) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "face_map",n
    call alloc_error(pname,array_name,1,ierr,error_message)
    !
    ! Add the faces of the current partition to face_map
    !
    face_map(n)%loc_to_glb = pack( intseq(1,ngface) , &
                                   face(1:ngface)%left%partition == n )
   !face_map(n)%loc_to_glb = find_subset(face(1:ngface)%left%partition,n)
    !
  end do
  !
  ! Format Statements
  !
  1 format (a,"(",i0,")%loc_to_glb")
  !
end subroutine create_grid_mappings
!
!###############################################################################
!
pure function find_subset(array,find_value) result(return_value)
  !
  !.. Formal Arguments ..
  integer,               intent(in) :: find_value
  integer, dimension(:), intent(in) :: array
  !
  !.. Function Return Value ..
  integer, allocatable, dimension(:) :: return_value
  !
  !.. Local Scalars ..
  integer :: ncount,ierr
  !
continue
  !
  ncount = count( array == find_value )
  !
  allocate ( return_value(1:ncount) , source=0 , stat=ierr )
  !
  return_value = pack( intseq(1,size(array)) , array == find_value )
  !
end function find_subset
!
!###############################################################################
!
subroutine localize_grid(nodes_of_cell,nodes_of_cell_ptr,bface, &
                         xyz_nodes,epart,cell_geom,cell_order, &
                         node_map,face_map,cell_map)
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  use geovar,    only : bc_conflict
  !
  !.. Formal Arguments ..
  integer,                dimension(:),      intent(in) :: epart
  type(map_t),            dimension(:),      intent(in) :: node_map
  type(map_t),            dimension(:),      intent(in) :: face_map
  type(map_t),            dimension(:),      intent(in) :: cell_map
  integer,   allocatable, dimension(:),   intent(inout) :: cell_geom
  integer,   allocatable, dimension(:),   intent(inout) :: cell_order
  integer,   allocatable, dimension(:),   intent(inout) :: nodes_of_cell
  integer,   allocatable, dimension(:),   intent(inout) :: nodes_of_cell_ptr
  integer,   allocatable, dimension(:,:), intent(inout) :: bface
  real(wp),  allocatable, dimension(:,:), intent(inout) :: xyz_nodes
  !
  !.. Local Scalars ..
  integer :: n,nf,np,nmin,nmax,ierr
  integer :: nfbnd,ncell,nnode,nface
  integer :: n_totcel,n_totnod,nfnods
  integer :: rev_glb_node,glb_geom,glb_order
  integer :: l1,l2,g1,g2,loc_ghost_cell
  integer :: loc_face,loc_cell,loc_node
  integer :: glb_face,glb_cell,glb_node
  integer :: adj_face,adj_cell,adj_part
  integer :: my_part,this_part,master_cpu,slave_cpu
  !
  !.. Local Allocatable Arrays ..
  logical,   allocatable, dimension(:)   :: msk
  integer,   allocatable, dimension(:)   :: recv_blklen
  integer,   allocatable, dimension(:)   :: recv_offset
  integer,   allocatable, dimension(:)   :: cell_offset
  integer,   allocatable, dimension(:)   :: ghostcell_map
  integer,   allocatable, dimension(:)   :: ghostnode_map
  integer,   allocatable, dimension(:)   :: ghostnode_rmap
  integer,   allocatable, dimension(:)   :: global_geom
  integer,   allocatable, dimension(:)   :: global_order
  integer,   allocatable, dimension(:)   :: global_nofc
  integer,   allocatable, dimension(:)   :: global_nofc_ptr
  integer,   allocatable, dimension(:,:) :: global_bface
  real(wp),  allocatable, dimension(:,:) :: global_xyz
  !
  !.. Objects of MPI Integer Type ..
  integer(int_mpi) :: np_int_mpi
  integer(int_mpi), allocatable, dimension(:) :: recv_blklen_int_mpi
  integer(int_mpi), allocatable, dimension(:) :: recv_offset_int_mpi
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "localize_grid"
  !
continue
  !
  ! Before we start localizing everything, localize bc_conflict.
  ! Go through each bc conflict and get the partition to which the
  ! boundary face belongs and then localize the boundary face indices.
  ! If we find a face that is on our partition, mark the conflict so
  ! that we know we are involved in the communication for this conflict.
  !
  if (allocated(bc_conflict)) then
    !
    do n = 1,size(bc_conflict)
      !
      ! Reinitialize bc_conflict(n)%not_my_conflict to true and only
      ! change if we find a face on this conflict for this partition.
      !
      bc_conflict(n)%not_my_conflict = true
      !
      ! First get the partition for the master face
      !
      glb_face = bc_conflict(n)%master%face
      !
      master_cpu = epart( bface(2,glb_face) )
      !
      loc_face = findloc( face_map(master_cpu)%loc_to_glb , glb_face)
      !
      bc_conflict(n)%master%face = loc_face
      bc_conflict(n)%master%cpu = master_cpu
      !
      ! Mark this conflict if this processor is involved with it
      !
      if (master_cpu == mypnum) bc_conflict(n)%not_my_conflict = fals
      !
      ! Now do the same for the slave faces
      !
      do nf = 1,size(bc_conflict(n)%slave)
        !
        glb_face = bc_conflict(n)%slave(nf)%face
        !
        slave_cpu = epart( bface(2,glb_face) )
        !
        loc_face = findloc( face_map(slave_cpu)%loc_to_glb , glb_face)
        !
        bc_conflict(n)%slave(nf)%face = loc_face
        bc_conflict(n)%slave(nf)%cpu = slave_cpu
        !
        ! Mark this conflict if this processor is involved with it
        !
        if (slave_cpu == mypnum) bc_conflict(n)%not_my_conflict = fals
        !
        ! Check if this conflict is local to only one partition.
        ! NOTE: We first check bc_conflict(n)%is_local to see if it is true
        !       before checking if the master and slave processors are the 
        !       same. This way, if we already determined that this conflict
        !       is not local, we dont have to worry about the following
        !       logical assignment overwriting the value in is_local with
        !       the value 'true'. For example:
        !       master_cpu = 1
        !       slave_cpu #1 = 2
        !       bc_conflict(n)%is_local = (master_cpu == slave_cpu #1) = .false.
        !       slave_cpu #2 = 1
        !       bc_conflict(n)%is_local = (master_cpu == slave_cpu #2) = .true.
        !
        if (bc_conflict(n)%is_local) then
          bc_conflict(n)%is_local = (master_cpu == slave_cpu)
        end if
        !
      end do
      !
    end do
    !
   !do n = 1,size(bc_conflict)
   !  write (mypnum+200,1) n,bc_conflict(n)%is_local, &
   !                         bc_conflict(n)%not_my_conflict, &
   !                         bc_conflict(n)%master%face, &
   !                         bc_conflict(n)%master%flux_point, &
   !                         bc_conflict(n)%master%cpu, &
   !          trim(bc_integer_to_string(bface(1,bc_conflict(n)%master%face)))
   !  do nf = 1,size(bc_conflict(n)%slave)
   !    write (mypnum+200,2) nf,bc_conflict(n)%slave(nf)%face, &
   !                            bc_conflict(n)%slave(nf)%flux_point, &
   !                            bc_conflict(n)%slave(nf)%cpu, &
   !          trim(bc_integer_to_string(bface(1,bc_conflict(n)%slave(nf)%face)))
   !  end do
   !end do
   !1 format ("Conflict #",i0," - is_local=",l1,", not_my_conflict=",l1,/, &
   !          "MF :   Face=",i5,",  FluxPoint=",i2,",  CPU=",i0, &
   !          "    BC=",a)
   !2 format ("SF",i0,":   Face=",i5,",  FluxPoint=",i2,",  CPU=",i0, &
   !          "    BC=",a)
    !
  end if
  !
  !
  !
  ! We need to make copies of the global grid arrays nodes_of_cell,
  ! nodes_of_cell_ptr, bface, cell_geom, cell_order, and xyz_nodes
  ! before we localize them. Since we need to reallocate the original
  ! arrays to the new local grid sizes, use the intrinsic subroutine
  ! move_alloc which does:
  !     1. Deallocates 'to=' array if it is allocated
  !     2. Allocates 'to=' array with the same type, type parameters,
  !        array bounds, and with the same value as 'from=' array
  !     3. Deallocates 'from=' array
  !
  call move_alloc( from=nodes_of_cell     , to=global_nofc )
  call move_alloc( from=nodes_of_cell_ptr , to=global_nofc_ptr )
  call move_alloc( from=bface             , to=global_bface )
  call move_alloc( from=cell_geom         , to=global_geom )
  call move_alloc( from=cell_order        , to=global_order )
  call move_alloc( from=xyz_nodes         , to=global_xyz )
  !
  ! The process ID, mypnum, starts numbering at 0 so
  ! identify which partition is local to each processor
  !
  my_part = mypnum + 1
  !
  ! Count the local number of cells, nodes, faces, and boundary faces
  !
  ncell = size( cell_map(my_part)%loc_to_glb )
  nnode = size( node_map(my_part)%loc_to_glb )
  nface = size( face_map(my_part)%loc_to_glb )
  nfbnd = count( face_map(my_part)%loc_to_glb <= ngfbnd )
  !
  ! Now that all ghost cells are actually faces and not cells,
  ! n_totnod is just the same as nnode
  !
  n_totnod = nnode
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Localize the boundary face array, but not the information it contains.
  !
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,ierr,error_message)
  !
  ! First extract the information that is local to the current
  ! processor from the global_bface array to the bface array.
  ! Also create a mapping array to identify the global ghost cells
  ! needed for the local grid.
  !
  do loc_face = 1,nfbnd
    !
    ! Global face ID for this local boundary face
    !
    glb_face = face_map(my_part)%loc_to_glb(loc_face)
    !
    ! Copy everything to the local bface array
    !
    bface(:,loc_face) = global_bface(:,glb_face)
    !
    ! Now localize the information for this face
    !
    ! local host cell
    bface(2,loc_face) = findloc( cell_map(my_part)%loc_to_glb , &
                                 bface(2,loc_face) )
    !
    if ( any(bface(1,loc_face) == bc_comm) ) then
      !
      adj_face = bface(5,loc_face)      ! connecting boundary face
      adj_cell = bface(6,loc_face)      ! connecting host cell
      adj_part = bface(7,loc_face) + 1  ! partition of connecting cell
      !
      bface(5,loc_face) = findloc( face_map(adj_part)%loc_to_glb , adj_face)
      bface(6,loc_face) = findloc( cell_map(adj_part)%loc_to_glb , adj_cell)
      !
    end if
    !
    ! Localize the ghost cell for this face
    ! NOTE: Since ghost cells are actually faces of an interior cell they are
    !       one dimensions lower that the interior cells, 
    !          e.g. 2D tria/quad face for a 3D grid or 
    !               1D edge face for a 2D grid
    !       each ghost cell is uniquely defined so the ghost cell index is just
    !        (the number of interior cells) + (index of current boundary face)
    !
    bface(10,loc_face) = ncell + loc_face
    !
  end do
  !
  n_totcel = ncell + nfbnd
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Localize the cell_geom and cell_order arrays
  !
  allocate ( cell_geom(1:n_totcel) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_geom",1,ierr,error_message)
  !
  allocate ( cell_order(1:n_totcel) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_order",1,ierr,error_message)
  !
  ! Get the geometry and order for the interior cells local to this partition
  !
  do loc_cell = 1,ncell
    glb_cell = cell_map(my_part)%loc_to_glb(loc_cell)
    cell_geom(loc_cell) = global_geom(glb_cell)
    cell_order(loc_cell) = global_order(glb_cell)
  end do
  !
  ! Get the geometry and order for the ghost cells of this partition
  !
  do loc_face = 1,nfbnd
    glb_face = face_map(my_part)%loc_to_glb(loc_face)
    loc_cell = ncell + loc_face
    cell_geom(loc_cell) = face(glb_face)%geom
    cell_order(loc_cell) = face(glb_face)%order
  end do
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Localize the node connectivity pointer for each cell
  !
  ! Create the local nodes_of_cell_ptr
  !
  allocate ( nodes_of_cell_ptr(1:n_totcel+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell_ptr",1,ierr,error_message)
  !
  ! First set the node connectivity pointer to the number
  ! of nodes that define each of the interior cells
  !
  do loc_cell = 1,ncell
    glb_cell = cell_map(my_part)%loc_to_glb(loc_cell)
    nodes_of_cell_ptr(loc_cell+1) = global_nofc_ptr(glb_cell+1) - &
                                    global_nofc_ptr(glb_cell)
  end do
  !
  ! Next set the node connectivity pointer to the number of nodes that
  ! define each of the ghost cells by using the face(:)%nodes array.
  !
  do loc_face = 1,nfbnd
    loc_ghost_cell = bface(10,loc_face)
    glb_face = face_map(my_part)%loc_to_glb(loc_face)
    nfnods = geom_nodes(face(glb_face)%geom) ! # of nodes defining this face
    nodes_of_cell_ptr(loc_ghost_cell+1) = nfnods
  end do
  !
  ! Reshuffle nodes_of_cell_ptr to make the pointer cumulative 
  !
  do loc_cell = 2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(loc_cell) = nodes_of_cell_ptr(loc_cell) + &
                                  nodes_of_cell_ptr(loc_cell-1)
  end do
  !
  ! #######
  ! ## 4 ##
  ! #######
  !
  ! Create the local nodes_of_cell array for each cell
  ! NOTE: The node indices copied into this array are actually the global
  !       indices and they will not be localized until section 6.
  !
  allocate ( nodes_of_cell(1:last(nodes_of_cell_ptr)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell",1,ierr,error_message)
  !
  ! First extract the information that is local to the current
  ! processor from the global_nofc array to the nodes_of_cell array.
  !
  do loc_cell = 1,ncell
    !
    l1 = nodes_of_cell_ptr(loc_cell)+1
    l2 = nodes_of_cell_ptr(loc_cell+1)
    !
    glb_cell = cell_map(my_part)%loc_to_glb(loc_cell)
    !
    g1 = global_nofc_ptr(glb_cell)+1
    g2 = global_nofc_ptr(glb_cell+1)
    !
    nodes_of_cell(l1:l2) = global_nofc(g1:g2)
    !
  end do
  !
  ! Now go through the boundary faces and extract the nodal information for
  ! the ghost cells from the faces(:)%nodes array.
  !
  do loc_face = 1,nfbnd
    !
    glb_face = face_map(my_part)%loc_to_glb(loc_face)
    !
    loc_ghost_cell = bface(10,loc_face)
    !
    l1 = nodes_of_cell_ptr(loc_ghost_cell)+1
    l2 = nodes_of_cell_ptr(loc_ghost_cell+1)
    !
    nfnods = geom_nodes(face(glb_face)%geom) ! # of nodes defining this face
    !
    nodes_of_cell(l1:l2) = face(glb_face)%nodes(1:nfnods)
    !
  end do
  !
 !if (mypnum == 0) then
 !  write (mypnum+130,'(a,i0)') "ncell = ",ngcell
 !  do loc_cell = 1,size(global_nofc_ptr)-1
 !    write (mypnum+130,*)
 !    l1 = global_nofc_ptr(loc_cell)+1
 !    l2 = global_nofc_ptr(loc_cell+1)
 !    write (mypnum+130,40) loc_cell,l1,l2
 !    write (mypnum+130,41) (global_nofc(n),n=l1,l2)
 !  end do
 !end if
 !!
 !flush (mypnum+130)
 !call mpi_barrier(mpi_comm_world,mpierr)
 !write (mypnum+120,'(a,i0)') "ncell = ",ncell
 !do loc_cell = 1,size(nodes_of_cell_ptr)-1
 !  write (mypnum+120,*)
 !  l1 = nodes_of_cell_ptr(loc_cell)+1
 !  l2 = nodes_of_cell_ptr(loc_cell+1)
 !  write (mypnum+120,40) loc_cell,l1,l2
 !  write (mypnum+120,41) (nodes_of_cell(n),n=l1,l2)
 !end do
 !!
 !40 format ("loc_cell=",i0,";   l1:l2=",i0,":",i0)
 !41 format ("nodes_of_cell(l1:l2) = ",*(i0,:,", "))
 !!
 !call stop_recon(stop_mpi,pname)
  !
  ! We no longer need global_bface or global_nofc so deallocate
  !
  deallocate ( global_bface , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"global_bface",2,ierr,error_message)
  deallocate ( global_nofc , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"global_nofc",2,ierr,error_message)
  !
  ! #######
  ! ## 5 ##
  ! #######
  !
  ! We need to create collective receiving MPI datatypes to collect
  ! localized data from each processor into a shared global form.
  !
  ! First, we need to create an array containing the
  ! solution point offsets for every global cell
  !
  allocate ( cell_offset(1:ngcell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_offset",1,ierr,error_message)
  !
  do glb_cell = 1,ngcell-1
    !
    glb_geom = global_geom(glb_cell)
    glb_order = global_order(glb_cell)
    !
    cell_offset(glb_cell+1) = cell_offset(glb_cell) + &
                              geom_solpts(glb_geom,glb_order)
    !
  end do
  !
  ! Loop over all the processors to create the collective receive MPI datatypes
  !
  np = 0
  do this_part = 1,size(cell_map)
    np = max( np , size(cell_map(this_part)%loc_to_glb) )
  end do
  !
  allocate ( recv_offset(1:np) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_offset",1,ierr,error_message)
  !
  allocate ( recv_offset_int_mpi(1:np) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_offset_int_mpi",1,ierr,error_message)
  !
  allocate ( recv_blklen(1:np) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_blklen",1,ierr,error_message)
  !
  allocate ( recv_blklen_int_mpi(1:np) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_blklen_int_mpi",1,ierr,error_message)
  !
  allocate ( recv_xyz(1:ncpu) , source=mpi_datatype_null , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_xyz",1,ierr,error_message)
  !
  allocate ( recv_sol(1:ncpu) , source=mpi_datatype_null , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_sol",1,ierr,error_message)
  !
  allocate ( recv_var(1:ncpu) , source=mpi_datatype_null , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_var",1,ierr,error_message)
  !
  do this_part = 1,ncpu
    !
    np = size(cell_map(this_part)%loc_to_glb)
    np_int_mpi = int( np , kind=int_mpi )
    !
    do loc_cell = 1,np
      !
      glb_cell = cell_map(this_part)%loc_to_glb(loc_cell)
      glb_geom = global_geom(glb_cell)
      glb_order = global_order(glb_cell)
      !
      recv_offset(loc_cell) = cell_offset(glb_cell)
      recv_blklen(loc_cell) = geom_solpts(glb_geom,glb_order)
      !
    end do
    !
    ! Now create the receive MPI datatype for the local
    ! solution point coordinates on the current processor
    !
    ! adjust the offset for nr dimensions
    recv_offset(1:np) = recv_offset(1:np)*nr
    recv_blklen(1:np) = recv_blklen(1:np)*nr
    !
    recv_offset_int_mpi(1:np) = int( recv_offset(1:np) , kind=int_mpi )
    recv_blklen_int_mpi(1:np) = int( recv_blklen(1:np) , kind=int_mpi )
    !
    call mpi_type_indexed(np_int_mpi, &
                          recv_blklen_int_mpi(1:np), &
                          recv_offset_int_mpi(1:np), &
                          mpi_flttyp,recv_xyz(this_part),mpierr)
    call mpi_type_commit(recv_xyz(this_part),mpierr)
    !
    ! cancel out the nr adjustment above
    recv_offset(1:np) = recv_offset(1:np)/nr
    recv_blklen(1:np) = recv_blklen(1:np)/nr
    !
    ! Now create the receive MPI datatype for the local
    ! solution at the solution points on the current processor
    !
    ! adjust the offset for nq flow equations
    recv_offset(1:np) = recv_offset(1:np)*nq
    recv_blklen(1:np) = recv_blklen(1:np)*nq
    !
    recv_offset_int_mpi(1:np) = int( recv_offset(1:np) , kind=int_mpi )
    recv_blklen_int_mpi(1:np) = int( recv_blklen(1:np) , kind=int_mpi )
    !
    call mpi_type_indexed(np_int_mpi, &
                          recv_blklen_int_mpi(1:np), &
                          recv_offset_int_mpi(1:np), &
                          mpi_flttyp,recv_sol(this_part),mpierr)
    call mpi_type_commit(recv_sol(this_part),mpierr)
    !
    ! cancel out the nq adjustment above
    recv_offset(1:np) = recv_offset(1:np)/nq
    recv_blklen(1:np) = recv_blklen(1:np)/nq
    !
    ! Now create the receive MPI datatype for an arbitrary
    ! variable at the solution points on the current processor
    !
    recv_offset_int_mpi(1:np) = int( recv_offset(1:np) , kind=int_mpi )
    recv_blklen_int_mpi(1:np) = int( recv_blklen(1:np) , kind=int_mpi )
    !
    call mpi_type_indexed(np_int_mpi, &
                          recv_blklen_int_mpi(1:np), &
                          recv_offset_int_mpi(1:np), &
                          mpi_flttyp,recv_var(this_part),mpierr)
    call mpi_type_commit(recv_var(this_part),mpierr)
    !
  end do
  !
  ! We no longer need global_nofc_ptr, recv_offset,
  ! recv_blklen, or cell_offset so deallocate
  !
  deallocate ( global_nofc_ptr , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"global_nofc_ptr",2,ierr,error_message)
  deallocate ( recv_offset , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_offset",2,ierr,error_message)
  deallocate ( recv_blklen , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_blklen",2,ierr,error_message)
  deallocate ( recv_offset_int_mpi , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_offset_int_mpi",2,ierr,error_message)
  deallocate ( recv_blklen_int_mpi , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_blklen_int_mpi",2,ierr,error_message)
  deallocate ( cell_offset , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_offset",2,ierr,error_message)
  !
  ! #######
  ! ## 6 ##
  ! #######
  !
  ! Reverse map the interior nodes to the local ordering
  ! NOTE: The node indices stored in nodes of cells are currently the
  !       global node indices
  !
  do loc_node = 1,size(nodes_of_cell)
    glb_node = nodes_of_cell(loc_node)
    nodes_of_cell(loc_node) = findloc( node_map(my_part)%loc_to_glb , glb_node )
  end do
  !
  ! #######
  ! ## 7 ##
  ! #######
  !
  ! Localize the node coordinates
  !
  allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_nodes",1,ierr,error_message)
  !
  ! Localize the interior nodes
  !
  do loc_node = 1,nnode
    glb_node = node_map(my_part)%loc_to_glb(loc_node)
    xyz_nodes(:,loc_node) = global_xyz(:,glb_node)
  end do
  !
  ! We no longer need global_xyz so deallocate
  !
  deallocate ( global_xyz , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"global_xyz",2,ierr,error_message)
  !
  ! #######
  ! ## 8 ##
  ! #######
  !
  ! Now that we have localized the grid, change any periodic boundary conditions
  ! that require communication across processors to communication boundaries
  !
  do nf = 1,nfbnd
    if (bface(1,nf) == bc_periodic .and. bface(7,nf) /= mypnum) then
      bface(1,nf) = bc_cpu_bnd
    end if
  end do
 !!
 !! #######
 !! ## 9 ##
 !! #######
 !!
 !! Output the individual partitions to separate Tecplot files
 !!
 !call output_partitions(my_part,nodes_of_cell,nodes_of_cell_ptr,xyz_nodes)
 !!
 !call output_partitions(my_part+ncpu+1,nodes_of_cell, &
 !                       nodes_of_cell_ptr(1:ncell+1), &
 !                       xyz_nodes(1:nr,1:nnode))
  !
end subroutine localize_grid
!
!###############################################################################
!
subroutine create_face_nonblocking_datatypes(bface,face_map)
  !
  !.. Use Statements ..
  use geovar, only : n_global_node
  use geovar, only : n_global_edge
  !
  !.. Formal Arguments ..
  integer,     dimension(:,:), intent(in) :: bface
  type(map_t), dimension(:),   intent(in) :: face_map
  !
  !.. Local Scalars ..
  logical :: flag
  integer :: n,nfp,nfbnd,i,ierr
  integer :: rmin,rmax,scnt,rcnt,p1,p2
  integer :: loc_face,glb_face,my_part
  integer :: adj_cpu,adj_part,adj_face
  integer(int_mpi) :: block_size,np
  integer(int_mpi) :: this_flxpt,tag_ub
  character(len=200) :: array_name
  !
  !.. Local Derived Data Types
  type :: comm_pts_t
    integer(int_mpi), allocatable :: send(:)
    integer(int_mpi), allocatable :: recv(:)
  end type
  !
  !.. Local Arrays ..
  integer, dimension(1:ncpu) :: cpu_id
  !
  !.. Local Allocatable Arrays ..
  logical,          allocatable :: send_msk(:)
  integer,          allocatable :: recv_map(:)
  integer(int_mpi), allocatable :: flxpt_offset(:)
  type(comm_pts_t), allocatable :: comm_pts(:)
  !
 !real(r8), dimension(1:nr,1:n_global_totpts,1:ncpu) :: gxyz
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_face_nonblocking_datatypes"
  !
continue
  !
  my_part = mypnum + 1
  !
  nfbnd = size(bface,dim=2)
  !
 !call stop_recon(stop_mpi,"NEED TO SET UP FACEDUSP COMMUNICATION FIRST!!!")
 !!
 !call get_solution_points_2d
 !do n = 1,ncpu
 !  if (mypnum == n-1) np = n_totpts
 !  call mpi_bcast(np,1,mpi_inttyp,n-1,mpi_comm_world,mpierr)
 !  if (mypnum == n-1) gxyz(1:nr,1:np,n) = real(xyz(1:nr,1:n_totpts),kind=r8)
 !  call mpi_bcast(gxyz(1:nr,1:np,n),np*nr,mpi_flttyp,n-1,mpi_comm_world,mpierr)
 !end do
  !
  ! Create an offset array for the number of flux
  ! points on each communication boundary faces
  !
  allocate ( flxpt_offset(1:nfbnd) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"flxpt_offset",1,ierr,error_message)
  !
  do loc_face = 1,nfbnd
    !
    if (bface(1,loc_face) /= bc_cpu_bnd) cycle
    !
    glb_face = face_map(my_part)%loc_to_glb(loc_face)
    !
    flxpt_offset(loc_face) = size( face(glb_face)%left%flxpts , kind=int_mpi )
    !
  end do
  !
  ! Reshuffle flxpt_offset by summing the total number of flux points for
  ! all the communication boundary faces before the current one.
  !
  do loc_face = nfbnd,1,-1
    !
    if (bface(1,loc_face) /= bc_cpu_bnd) cycle
    !
    flxpt_offset(loc_face) = sum( flxpt_offset(1:loc_face-1) )
    !
  end do
  !
  ! Now onto creating the communication datatypes
  !
  ! First we need to find which processes the
  ! current process needs to exchange with
  !
  cpu_id(:) = -99
  do loc_face = 1,nfbnd
    if (bface(1,loc_face) == bc_cpu_bnd) then
      cpu_id( bface(7,loc_face)+1 ) = bface(7,loc_face)
    end if
  end do
  !
  ! Make sure that the current process doesnt try to communicate with itself
  !
  cpu_id(my_part) = -99
  !
  ! Count the number of processes the current process needs to communicate with
  !
  ncomm = count(cpu_id >= 0)
  !
  ! Allocate the arrays that will contain the communication information
  !
  allocate ( exch_fusp(1:ncomm) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"exch_fusp",1,ierr,error_message)
  !
  allocate ( exch_fdusp(1:ncomm) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"exch_fusp",1,ierr,error_message)
  !
  ! Copy the process IDs that the current process needs to communicate
  ! with into these arrays.
  !
  exch_fusp(1:ncomm)%adj_cpu = int(pack( cpu_id , cpu_id >= 0 ),kind=int_mpi)
  exch_fdusp(1:ncomm)%adj_cpu = exch_fusp(1:ncomm)%adj_cpu
  !
  ! Find the range of face indices that this process will be
  ! receiving from all of its adjacent partitions
  !
  rmin = minval(bface(5,:)) ! min face index this process will be receiving from
  rmax = maxval(bface(5,:)) ! max face index this process will be receiving from
  !
  ! Allocate all the arrays needed for creating the send and receive face maps
  !
  allocate ( send_msk(1:nfbnd) , source=fals , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"send_msk",1,ierr,error_message)
  !
  allocate ( recv_map(rmin:rmax) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_map",1,ierr,error_message)
  !
  allocate ( comm_pts(1:ncomm) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"comm_pts",1,ierr,error_message)
  !
  ! Loop through each of the processes this one will be communicating with
  ! and create the send and receive face maps
  !
  do i = 1,ncomm
    !
    adj_cpu = int( exch_fusp(i)%adj_cpu , kind=kind(adj_cpu) )
    adj_part = adj_cpu + 1
    !
    ! Loop to find the communication boundary faces between the
    ! current processor and the current exchange partner, adj_cpu.
    ! We will create a mask array, send_msk, to identify all the
    ! local host cells that need to be sent to the ghost cells
    ! of adj_cpu. We will also create an array, recv_map, that maps
    ! the host cells of adj_cpu to the local ghost cells whose data
    ! needs to be received from adj_cpu.
    !
    do loc_face = 1,nfbnd
      !
      if (bface(1,loc_face) /= bc_cpu_bnd) cycle
      !
      if (bface(7,loc_face) == adj_cpu) then
        adj_face = bface(5,loc_face)
        send_msk(loc_face) = true
        recv_map(adj_face) = loc_face
      end if
      !
    end do
    !
    ! Now create the send and receive face maps. These maps will contain
    ! in order the faces that need to be sent to and received from the
    ! current adjacent partition adj_cpu
    ! NOTE: Each processor will send its face data in order of the face IDs.
    !       The receiving datatypes need to account for this so that the
    !       received data is placed in the correct location within an array.
    !       This is why we have created the array recv_map so that we can
    !       loop in order through the faces of the adjacent partition to
    !       create the receive offsets in the same order of the adjacent
    !       partitions face IDs.
    !
    ! Create the send face maps
    !
    !
    ! Count how many flux points in total will need to be communicated.
    ! NOTE: This should be the same for both sending and receiving
    !
    scnt = 0
    do loc_face = 1,nfbnd
      !
      if (.not. send_msk(loc_face)) cycle
      !
      glb_face = face_map(my_part)%loc_to_glb(loc_face)
      !
      scnt = scnt + size( face(glb_face)%left%flxpts )
      !
    end do
    !
    ! Allocate the send and recv components of comm_pts now that we know
    ! the total number of flux points that will need to be communicated.
    !
    allocate ( comm_pts(i)%send(1:scnt) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,100) i,"send"
    call alloc_error(pname,array_name,1,ierr,error_message)
    !
    allocate ( comm_pts(i)%recv(1:scnt) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,100) i,"recv"
    call alloc_error(pname,array_name,1,ierr,error_message)
    !
    ! Create the send point offsets
    !
    scnt = 0
    do loc_face = 1,nfbnd
      !
      if (.not. send_msk(loc_face)) cycle
      !
      glb_face = face_map(my_part)%loc_to_glb(loc_face)
      !
      do n = 1,size( face(glb_face)%left%flxpts )
        !
        scnt = scnt + 1
        !
        this_flxpt = int( face(glb_face)%left%flxpts(n) , kind=int_mpi )
        !
        comm_pts(i)%send(scnt) = flxpt_offset(loc_face) + this_flxpt
        !
      end do
      !
    end do
    !
    ! Create the receive point offsets
    !
    rcnt = 0
    do adj_face = rmin,rmax
      !
      if (recv_map(adj_face) == 0) cycle
      !
      loc_face = recv_map(adj_face)
      !
      glb_face = face_map(my_part)%loc_to_glb(loc_face)
      !
      do n = 1,size( face(glb_face)%left%flxpts )
        !
        rcnt = rcnt + 1
        !
        this_flxpt = int( face(glb_face)%left%flxpts(n) , kind=int_mpi )
        !
        comm_pts(i)%recv(rcnt) = flxpt_offset(loc_face) + this_flxpt
        !
      end do
      !
    end do
    !
    send_msk(:) = fals
    recv_map(:) = 0
    !
  end do
  !
  ! Deallocate send_msk, recv_map, and flxpt_offset
  ! because they are no longer needed
  !
  deallocate ( send_msk , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"send_msk",2,ierr,error_message)
  deallocate ( recv_map , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_map",2,ierr,error_message)
  deallocate ( flxpt_offset , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"flxpt_offset",2,ierr,error_message)
  !
  !##########################################################
  !##########################################################
  ! Create the MPI datatypes for exchanging
  ! solution data between adjacent partitions
  !##########################################################
  !##########################################################
  !
  ! Set the block size for this MPI datatype
  !
  block_size = int( nq , kind=int_mpi )
  !
  do i = 1,ncomm
    !
    ! Adjust comm_pts(i)%send and comm_pts(i)%recv for the solution arrays
    !
    comm_pts(i)%send = comm_pts(i)%send * block_size
    comm_pts(i)%recv = comm_pts(i)%recv * block_size
    !
    np = size( comm_pts(i)%send , kind=int_mpi )
    !
    call mpi_type_create_indexed_block(np,block_size,comm_pts(i)%send, &
                                       mpi_flttyp,exch_fusp(i)%send_datatype, &
                                       mpierr)
    call mpi_type_commit(exch_fusp(i)%send_datatype,mpierr)
    !
    np = size( comm_pts(i)%recv , kind=int_mpi )
    !
    call mpi_type_create_indexed_block(np,block_size,comm_pts(i)%recv, &
                                       mpi_flttyp,exch_fusp(i)%recv_datatype, &
                                       mpierr)
    call mpi_type_commit(exch_fusp(i)%recv_datatype,mpierr)
    !
    ! Remove the solution array adjustment from
    ! comm_pts(i)%send and comm_pts(i)%recv
    !
    comm_pts(i)%send = comm_pts(i)%send / block_size
    comm_pts(i)%recv = comm_pts(i)%recv / block_size
    !
  end do
  !
  !##########################################################
  !##########################################################
  ! Create the MPI datatypes for exchanging solution
  ! gradients data between adjacent partitions
  !##########################################################
  !##########################################################
  !
  ! Set the block size for this MPI datatype
  !
  block_size = int( nr*nq , kind=int_mpi )
  !
  do i = 1,ncomm
    !
    ! Adjust comm_pts(i)%send and comm_pts(i)%recv for the solution arrays
    !
    comm_pts(i)%send = comm_pts(i)%send * block_size
    comm_pts(i)%recv = comm_pts(i)%recv * block_size
    !
    np = size( comm_pts(i)%send , kind=int_mpi )
    !
    call mpi_type_create_indexed_block(np,block_size,comm_pts(i)%send, &
                                       mpi_flttyp,exch_fdusp(i)%send_datatype, &
                                       mpierr)
    call mpi_type_commit(exch_fdusp(i)%send_datatype,mpierr)
    !
    np = size( comm_pts(i)%recv , kind=int_mpi )
    !
    call mpi_type_create_indexed_block(np,block_size,comm_pts(i)%recv, &
                                       mpi_flttyp,exch_fdusp(i)%recv_datatype, &
                                       mpierr)
    call mpi_type_commit(exch_fdusp(i)%recv_datatype,mpierr)
    !
    ! Remove the solution array adjustment from
    ! comm_pts(i)%send and comm_pts(i)%recv
    !
    comm_pts(i)%send = comm_pts(i)%send / block_size
    comm_pts(i)%recv = comm_pts(i)%recv / block_size
    !
  end do
  !
  !
  !##########################################################
  !##########################################################
  ! Output all the solution point offsets for debugging
  !##########################################################
  !##########################################################
  !
 !do i = 1,ncomm
 !  !
 !  adj_cpu = exch_fusp(i)%adj_cpu
 !  !
 !  write (mypnum+180,1) mypnum,adj_cpu
 !  !
 ! !do n = size(comm_pts(i)%send)
 ! !  p1 = comm_pts(i)%send(n)+1
 ! !  p2 = comm_pts(i)%recv(n)+1
 ! !  write (mypnum+180,7) n,p1,gxyz(1,p1,mypnum+1),gxyz(2,p1,mypnum+1), &
 ! !                         p2,gxyz(1,p2,adj_cpu+1),gxyz(2,p2,adj_cpu+1)
 ! !end do
 !  !
 !end do
  !
  ! Create the test MPI datatypes
  !
  allocate ( testf_send(1:ncomm) , source=mpi_datatype_null , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"testf_send",1,ierr,error_message)
  !
  allocate ( testf_recv(1:ncomm) , source=mpi_datatype_null , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"testf_recv",1,ierr,error_message)
  !
  ! Set the block size for this MPI datatype
  !
  block_size = int( nr , kind=int_mpi )
  !
  do i = 1,ncomm
    !
    ! Adjust comm_pts(i)%send and comm_pts(i)%recv for the solution arrays
    !
    comm_pts(i)%send = comm_pts(i)%send * block_size
    comm_pts(i)%recv = comm_pts(i)%recv * block_size
    !
    np = size( comm_pts(i)%send , kind=int_mpi )
    !
    call mpi_type_create_indexed_block(np,block_size,comm_pts(i)%send, &
                                       mpi_flttyp,testf_send(i),mpierr)
    call mpi_type_commit(testf_send(i),mpierr)
    !
    np = size( comm_pts(i)%recv , kind=int_mpi )
    !
    call mpi_type_create_indexed_block(np,block_size,comm_pts(i)%recv, &
                                       mpi_flttyp,testf_recv(i),mpierr)
    call mpi_type_commit(testf_recv(i),mpierr)
    !
    ! Remove the solution array adjustment from
    ! comm_pts(i)%send and comm_pts(i)%recv
    !
    comm_pts(i)%send = comm_pts(i)%send / block_size
    comm_pts(i)%recv = comm_pts(i)%recv / block_size
    !
  end do
  !
  ! We no longer need comm_pts so deallocate
  !
  deallocate ( comm_pts , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"comm_pts",1,ierr,error_message)
  !
  ! Set the tag values for MPI communication
  !
  ! Get the upper bound for tag from MPI
  !
 !call mpi_comm_get_attr(MPI_COMM_WORLD,MPI_TAG_UB,tag_ub,flag,mpierr)
 !write (iout,*) mypnum,huge(tag_ub)
 !write (iout,*) mypnum,mpierr
 !write (iout,*) mypnum,flag
 !write (iout,*) mypnum,tag_ub
 !write (iout,*)
 !flush (iout)
 !call mpi_barrier(mpi_comm_world,mpierr)
 !!
 !! Divide the tag upper bound by two because it will be the starting
 !! value of tag for the gradient communications. This way the solution
 !! variable communication has the first half of valid tag values to work
 !! with and the gradient communication has the second half.
 !!
 !tag_ub = tag_ub / 2_int_mpi
 !!
 !write (iout,*) mypnum,tag_ub
 !flush (iout)
 !call mpi_barrier(mpi_comm_world,mpierr)
  !
  do i = 1,ncomm
    !
    exch_fusp(i)%send_tag = mypnum
    exch_fusp(i)%recv_tag = exch_fusp(i)%adj_cpu
    !
    exch_fdusp(i)%send_tag = mypnum!+ tag_ub
    exch_fdusp(i)%recv_tag = exch_fdusp(i)%adj_cpu!+ tag_ub
    !
  end do
  !
  ! Count the total number of flux points on all boundary faces
  !
  nfp = 0
  do loc_face = 1,nfbnd
    glb_face = face_map(my_part)%loc_to_glb(loc_face)
    nfp = nfp + size( face(glb_face)%left%flxpts )
  end do
  !
  ! Lets allocate the module arrays needed to
  ! exchange face data between processors
  !
  ! Arrays for exchanging solution variables at face flux points
  !
  ! Receive buffer for communicating solution variables
  !
  allocate ( fusp_recv_buffer(1:nq*nfp) , source=0.0_r8 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fusp_recv_buffer",1,ierr,error_message)
  !
  ! Send buffer for communicating solution variables
  !
  allocate ( fusp_send_buffer(1:nq*nfp) , source=0.0_r8 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fusp_send_buffer",1,ierr,error_message)
  !
  ! MPI request objects for communicating solution variables
  !
  allocate ( fusp_rqst(1:2*ncomm) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fusp_rqst",1,ierr,error_message)
  !
  ! MPI status objects returned after communicating solution variables
  !
  allocate ( fusp_istat(1:mpi_status_size,1:2*ncomm) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fusp_istat",1,ierr,error_message)
  !
  ! Arrays for exchanging gradients of solution variables at face flux points
  !
  ! Receive buffer for communicating gradients
  !
  allocate ( fdusp_recv_buffer(1:nr*nq*nfp) , source=0.0_r8 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fdusp_recv_buffer",1,ierr,error_message)
  !
  ! Send buffer for communicating gradients
  !
  allocate ( fdusp_send_buffer(1:nr*nq*nfp) , source=0.0_r8 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fdusp_send_buffer",1,ierr,error_message)
  !
  ! MPI request objects for communicating gradients
  !
  allocate ( fdusp_rqst(1:2*ncomm) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fdusp_rqst",1,ierr,error_message)
  !
  ! MPI status objects returned after communicating gradients
  !
  allocate ( fdusp_istat(1:mpi_status_size,1:2*ncomm) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fdusp_istat",1,ierr,error_message)
  !
  ! Buffer for communicating node/edge solutions and gradients
  !
  eusp_size  =    nq*n_global_edge
  edusp_size = nr*nq*n_global_edge
  nusp_size  =    nq*n_global_node
  ndusp_size = nr*nq*n_global_node
  !
  n = max(edusp_size,ndusp_size)
  !
  allocate ( neusp_buffer(1:n) , source=0.0_r8 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"neusp_buffer",1,ierr,error_message)
  !
  ! Format Statements
  !
  1 format (/," CPU ",i0," -> CPU ",i0," :: Sending")
  2 format (/," CPU ",i0," <- CPU ",i0," :: Receiving")
  3 format (*(12i6,/))
  4 format (" Cell = ",i0,*(/,12i6))
  5 format (" CPU ",i0," -> CPU ",i0," :: Sending   - ",3i6)
  6 format (" CPU ",i0," <- CPU ",i0," :: Receiving - ",3i6)
  7 format (2x,i6,": Sending - ",i0,4x,2f12.4,/, &
            "        Receiving - ",i0,4x,2f12.4)  
  100 format ("comm_pts(",i0,")%",a)
  !
end subroutine create_face_nonblocking_datatypes
!
!###############################################################################
!
subroutine exchange_faceusp(faceusp)
  !
  !.. Use Statements ..
#ifndef USE_TYPEBOUND_FACE_PACK
  use geovar,  only : nfbnd,bface
#endif
  use flowvar, only : face_var
  !
  !.. Formal Arguments ..
  type(face_var), dimension(:), intent(in) :: faceusp
  !
  !.. Local Scalars..
  integer :: n,nf,nsz,ibeg,iend
  !
  integer(int_mpi) :: send_cpu,send_dtype
  integer(int_mpi) :: recv_cpu,recv_dtype
  integer(int_mpi) :: mesg_tag
  !
continue
  !
#ifdef USE_TYPEBOUND_FACE_PACK
  call faceusp%pack(fusp_send_buffer)
#else
  iend = 0
  do nf = 1,nfbnd
    if (bface(1,nf) /= bc_cpu_bnd) cycle
    nsz = size(faceusp(nf)%v(:,:,1))
    ibeg = iend + 1
    iend = iend + nsz
    fusp_send_buffer(ibeg:iend) = &
                 reshape( real(faceusp(nf)%v(:,:,1),kind=r8) , [nsz] )
  end do
#endif
  !
  ! Initialize all the receives so each process
  ! is ready once the send is posted
  !
  do n = 1,ncomm
    !
    ! Process sending this message
    send_cpu   = exch_fusp(n)%adj_cpu
    ! Message tag
    mesg_tag   = exch_fusp(n)%recv_tag
    ! MPI user-defined datatype for receive
    recv_dtype = exch_fusp(n)%recv_datatype
    !
    call mpi_irecv(fusp_recv_buffer,1,recv_dtype,send_cpu,mesg_tag, &
                   mpi_comm_world,fusp_rqst(n),mpierr)
    !
  end do
  !
  ! Send out the face flux point solution variables to the adjacent processes
  !
  do n = 1,ncomm
    !
    ! Process receiving this message
    recv_cpu   = exch_fusp(n)%adj_cpu
    ! Message tag
    mesg_tag   = exch_fusp(n)%send_tag
    ! MPI user-defined datatype for send
    send_dtype = exch_fusp(n)%send_datatype
    !
    call mpi_isend(fusp_send_buffer,1,send_dtype,recv_cpu,mesg_tag, &
                   mpi_comm_world,fusp_rqst(n+ncomm),mpierr)
    !
  end do
  !
end subroutine exchange_faceusp
!
!###############################################################################
!
subroutine wait_faceusp(faceusp)
  !
  !.. Use Statements ..
#ifndef USE_TYPEBOUND_FACE_PACK
  use geovar,  only : nfbnd,bface
#endif
  use flowvar, only : face_var
  !
  !.. Formal Arguments ..
  type(face_var), dimension(:), intent(inout) :: faceusp
  !
  !.. Local Scalars ..
  integer :: ibeg,iend,nf
  !
  !.. Local Arrays ..
  integer, dimension(1:2) :: nsz
  !
continue
  !
  ! Wait for all faceusp sends and receives to complete
  !
  call mpi_waitall(int(2*ncomm,kind=int_mpi),fusp_rqst,fusp_istat,mpierr)
  !
  ! Update the right states for only the communication boundary faces
  !
#ifdef USE_TYPEBOUND_FACE_PACK
  call faceusp%unpack(fusp_recv_buffer)
#else
  iend = 0
  do nf = 1,nfbnd
    if (bface(1,nf) /= bc_cpu_bnd) cycle
    nsz(1) = size(faceusp(nf)%v,dim=1)
    nsz(2) = size(faceusp(nf)%v,dim=2)
    ibeg = iend + 1
    iend = iend + product(nsz)
    faceusp(nf)%v(:,:,2) = &
                 reshape( real(fusp_recv_buffer(ibeg:iend),kind=wp) , nsz )
  end do
#endif
  !
end subroutine wait_faceusp
!
!###############################################################################
!
subroutine exchange_facedusp(facedusp)
  !
  !.. Use Statements ..
#ifndef USE_TYPEBOUND_FACE_PACK
  use geovar,  only : nfbnd,bface
#endif
  use flowvar, only : face_vec
  !
  !.. Formal Arguments ..
  type(face_vec), dimension(:), intent(in) :: facedusp
  !
  !.. Local Scalars..
  integer :: n,nf,nsz,ibeg,iend
  !
  integer(int_mpi) :: send_cpu,send_dtype
  integer(int_mpi) :: recv_cpu,recv_dtype
  integer(int_mpi) :: mesg_tag
  !
continue
  !
#ifdef USE_TYPEBOUND_FACE_PACK
  call facedusp%pack(fdusp_send_buffer)
#else
  iend = 0
  do nf = 1,nfbnd
    if (bface(1,nf) /= bc_cpu_bnd) cycle
    nsz = size(facedusp(nf)%v(:,:,:,1))
    ibeg = iend + 1
    iend = iend + nsz
    fdusp_send_buffer(ibeg:iend) = &
                 reshape( real(facedusp(nf)%v(:,:,:,1),kind=r8) , [nsz] )
  end do
#endif
  !
  ! Initialize all the receives so each process
  ! is ready once the send is posted
  !
  do n = 1,ncomm
    !
    ! Process sending this message
    send_cpu   = exch_fdusp(n)%adj_cpu
    ! Message tag
    mesg_tag   = exch_fdusp(n)%recv_tag
    ! MPI user-defined datatype for receive
    recv_dtype = exch_fdusp(n)%recv_datatype
    !
    call mpi_irecv(fdusp_recv_buffer,1,recv_dtype,send_cpu,mesg_tag, &
                   mpi_comm_world,fdusp_rqst(n),mpierr)
    !
  end do
  !
  ! Send out the adjacent cells
  !
  do n = 1,ncomm
    !
    ! Process receiving this message
    recv_cpu   = exch_fdusp(n)%adj_cpu
    ! Message tag
    mesg_tag   = exch_fdusp(n)%send_tag
    ! MPI user-defined datatype for send
    send_dtype = exch_fdusp(n)%send_datatype
    !
    call mpi_isend(fdusp_send_buffer,1,send_dtype,recv_cpu,mesg_tag, &
                   mpi_comm_world,fdusp_rqst(n+ncomm),mpierr)
    !
  end do
  !
end subroutine exchange_facedusp
!
!###############################################################################
!
subroutine wait_facedusp(facedusp)
  !
  !.. Use Statements ..
#ifndef USE_TYPEBOUND_FACE_PACK
  use geovar,  only : nfbnd,bface
#endif
  use flowvar, only : face_vec
  !
  !.. Formal Arguments ..
  type(face_vec), dimension(:), intent(inout) :: facedusp
  !
  !.. Local Scalars ..
  integer :: ibeg,iend,nf
  !.. Local Arrays ..
  !
  integer, dimension(1:3) :: nsz
  !
continue
  !
  ! Wait for all facedusp sends and receives to complete
  !
  call mpi_waitall(int(2*ncomm,kind=int_mpi),fdusp_rqst,fdusp_istat,mpierr)
  !
  ! Update the right states for only the communication boundary faces
  !
#ifdef USE_TYPEBOUND_FACE_PACK
  call facedusp%unpack(fdusp_recv_buffer)
#else
  iend = 0
  do nf = 1,nfbnd
    if (bface(1,nf) /= bc_cpu_bnd) cycle
    nsz(1) = size(facedusp(nf)%v,dim=1)
    nsz(2) = size(facedusp(nf)%v,dim=2)
    nsz(3) = size(facedusp(nf)%v,dim=3)
    ibeg = iend + 1
    iend = iend + product(nsz)
    facedusp(nf)%v(:,:,:,2) = &
                 reshape( real(fdusp_recv_buffer(ibeg:iend),kind=wp) , nsz )
  end do
#endif
  !
end subroutine wait_facedusp
!
!###############################################################################
!
subroutine exch_edge_and_node_solutions(nodeusp,nodedusp,edgeusp,edgedusp)
  !
  !.. Use Statements ..
  use flowvar, only : ne1d_t,ne2d_t,ne3d_t
  !
  !.. Formal Arguments ..
  type(ne1d_t), dimension(:), intent(inout) :: nodeusp
  type(ne2d_t), dimension(:), intent(inout) :: nodedusp
  type(ne2d_t), dimension(:), intent(inout) :: edgeusp
  type(ne3d_t), dimension(:), intent(inout) :: edgedusp
  !
  !.. Local Scalars..
  integer :: n,n1,n2,n3,nsz
  integer :: ibeg,iend,ierr
  integer :: this_pair,npts
  !
  !.. MPI integers ..
  integer(int_mpi) :: cpu1,tag1,send_dtype
  integer(int_mpi) :: cpu2,tag2,recv_dtype
  integer(int_mpi) :: istat(1:mpi_status_size)
  !
  !.. Local Allocatable Arrays ..
  real(r8), allocatable, dimension(:) :: send_buffer
  real(r8), allocatable, dimension(:) :: recv_buffer
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "exch_edge_and_node_solutions"
  !
continue
  !
  ! Get the total size of all the edge and node data
  !
  npts = 0
  !
  do n = 1,size(nodeusp)
    npts = npts + size(nodeusp(n)%v)
  end do
  !
  do n = 1,size(nodedusp)
    npts = npts + size(nodedusp(n)%v)
  end do
  !
  do n = 1,size(edgeusp)
    npts = npts + size(edgeusp(n)%v)
  end do
  !
  do n = 1,size(edgedusp)
    npts = npts + size(edgedusp(n)%v)
  end do
  !
  ! Allocate the send and receive buffer arrays
  !
  allocate ( send_buffer(1:npts) , source=0.0_r8 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"send_buffer",1,ierr,error_message)
  !
  allocate ( recv_buffer(1:npts) , source=0.0_r8 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_buffer",1,ierr,error_message)
  !
  ! Create the send buffer
  !
  iend = 0
  !
  do n = 1,size(nodeusp)
    nsz = size(nodeusp(n)%v(:))
   !ibeg = nodeusp(n)%beg_pt
   !iend = ibeg-1 + nsz
    ibeg = iend + 1
    iend = iend + nsz
    send_buffer(ibeg:iend) = reshape( real(nodeusp(n)%v(:),kind=r8) , [nsz] )
  end do
  !
  do n = 1,size(nodedusp)
    nsz = size(nodedusp(n)%v(:,:))
   !ibeg = nodedusp(n)%beg_pt
   !iend = ibeg-1 + nsz
    ibeg = iend + 1
    iend = iend + nsz
    send_buffer(ibeg:iend) = &
                 reshape( real(nodedusp(n)%v(:,:),kind=r8) , [nsz] )
  end do
  !
  do n = 1,size(edgeusp)
    nsz = size(edgeusp(n)%v(:,:))
   !ibeg = edgeusp(n)%beg_pt
   !iend = ibeg-1 + nsz
    ibeg = iend + 1
    iend = iend + nsz
    send_buffer(ibeg:iend) = reshape( real(edgeusp(n)%v(:,:),kind=r8) , [nsz] )
  end do
  !
  do n = 1,size(edgedusp)
    nsz = size(edgedusp(n)%v(:,:,:))
   !ibeg = edgedusp(n)%beg_pt
   !iend = ibeg-1 + nsz
    ibeg = iend + 1
    iend = iend + nsz
    send_buffer(ibeg:iend) = &
                 reshape( real(edgedusp(n)%v(:,:,:),kind=r8) , [nsz] )
  end do
  !
  ! Loop through all the communication pairs and exchange solutions
  ! for the nodes and edge points, adding the received values to the
  ! current solutions before going to the next pair
  !
  exchange_loop: do this_pair = 1,size(exch_neusp)
    !
    tag1 = exch_neusp(this_pair)%tag1
    tag2 = exch_neusp(this_pair)%tag2
    !
    cpu1 = exch_neusp(this_pair)%cpu(1)%cpu_pnum
    cpu2 = exch_neusp(this_pair)%cpu(2)%cpu_pnum
    !
    if (mypnum == cpu1) then
      !
      send_dtype = exch_neusp(this_pair)%cpu(1)%send_datatype
      recv_dtype = exch_neusp(this_pair)%cpu(1)%recv_datatype
      !
      call mpi_send(send_buffer,1,send_dtype,cpu2,tag1,mpi_comm_world,mpierr)
      !
      call mpi_recv(recv_buffer,1,recv_dtype,cpu2,tag2,mpi_comm_world,istat,mpierr)
      !
    else if (mypnum == cpu2) then
      !
      send_dtype = exch_neusp(this_pair)%cpu(2)%send_datatype
      recv_dtype = exch_neusp(this_pair)%cpu(2)%recv_datatype
      !
      call mpi_recv(recv_buffer,1,recv_dtype,cpu2,tag1,mpi_comm_world,istat,mpierr)
      !
      call mpi_send(send_buffer,1,send_dtype,cpu2,tag2,mpi_comm_world,mpierr)
      !
    end if
    !
    if (any(mypnum == [cpu1,cpu2])) then
      !
      ! Add the received data to the solutions at the nodes and edges
      ! NOTE: This will unfortunately add a bunch of zeros because of
      !       how its implemented
      !
      iend = 0
      !
      do n = 1,size(nodeusp)
        n1 = size(nodeusp(n)%v,dim=1)
       !ibeg = nodeusp(n)%beg_pt
       !iend = ibeg-1 + n1
        ibeg = iend + 1
        iend = iend + n1
        nodeusp(n)%v(:) = nodeusp(n)%v(:) + &
                   reshape( real(recv_buffer(ibeg:iend),kind=wp) , [n1] )
      end do
      !
      do n = 1,size(nodedusp)
        n1 = size(nodedusp(n)%v,dim=1)
        n2 = size(nodedusp(n)%v,dim=2)
       !ibeg = nodedusp(n)%beg_pt
       !iend = ibeg-1 + n1*n2
        ibeg = iend + 1
        iend = iend + n1*n2
        nodedusp(n)%v(:,:) = nodedusp(n)%v(:,:) + &
                   reshape( real(recv_buffer(ibeg:iend),kind=wp) , [n1,n2] )
      end do
      !
      do n = 1,size(edgeusp)
        n1 = size(edgeusp(n)%v,dim=1)
        n2 = size(edgeusp(n)%v,dim=2)
       !ibeg = edgeusp(n)%beg_pt
       !iend = ibeg-1 + n1*n2
        ibeg = iend + 1
        iend = iend + n1*n2
        edgeusp(n)%v(:,:) = edgeusp(n)%v(:,:) + &
                   reshape( real(recv_buffer(ibeg:iend),kind=wp) , [n1,n2] )
      end do
      !
      do n = 1,size(edgedusp)
        n1 = size(edgedusp(n)%v,dim=1)
        n2 = size(edgedusp(n)%v,dim=2)
        n3 = size(edgedusp(n)%v,dim=3)
       !ibeg = edgedusp(n)%beg_pt
       !iend = ibeg-1 + n1*n2*n3
        ibeg = iend + 1
        iend = iend + n1*n2*n3
        edgedusp(n)%v(:,:,:) = edgedusp(n)%v(:,:,:) + &
                   reshape( real(recv_buffer(ibeg:iend),kind=wp) , [n1,n2,n3] )
      end do
      !
      ! Reset the receiver buffer to zero for the next communication
      !
      recv_buffer = 0.0_r8
      !
    end if
    !
    ! Have all processors wait for the current communication to complete
    ! before going on to the next communication pair.  This is to ensure
    ! there are no hangups and shouldnt be too much of a slowdown since
    ! this is only done for writing the solution files.
    !
    call mpi_barrier(mpi_comm_world,mpierr)
    !
  end do exchange_loop
  !
  if (allocated(send_buffer)) then
    deallocate ( send_buffer , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"send_buffer",2,ierr,error_message)
  end if
  !
  if (allocated(recv_buffer)) then
    deallocate ( recv_buffer , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_buffer",2,ierr,error_message)
  end if
  !
end subroutine exch_edge_and_node_solutions
!
!###############################################################################
!
subroutine exchange_edgeusp(edgeusp)
  !
  !.. Use Statements ..
  use flowvar, only : ne2d_t
  !
  !.. Formal Arguments ..
  type(ne2d_t), dimension(:), intent(inout) :: edgeusp
  !
  !.. Local Scalars..
  integer :: ne,vsz,ibeg,iend
  !
  !.. Local Arrays ..
  integer, dimension(1:2) :: nsz
  !
continue
  !
  ! Initialize neusp_buffer to zero
  !
  neusp_buffer = 0.0_r8
  !
#ifdef USE_TYPEBOUND_FACE_PACK
  call edgeusp%pack(neusp_buffer(1:eusp_size))
#else
  do ne = 1,size(edgeusp)
    vsz = size(edgeusp(ne)%v(:,:))
    ibeg = edgeusp(ne)%beg_pt
    iend = ibeg-1 + vsz
    neusp_buffer(ibeg:iend) = reshape( real(edgeusp(ne)%v(:,:),kind=r8) , [vsz] )
  end do
#endif
  !
  ! Sum all the edge point values from all the processors
  !
  call mpi_allreduce(MPI_IN_PLACE,neusp_buffer(1:eusp_size),eusp_size, &
                     mpi_flttyp,MPI_SUM,mpi_comm_world,mpierr)
  !
  ! Update the right states for only the communication boundary faces
  !
#ifdef USE_TYPEBOUND_FACE_PACK
  call edgeusp%unpack(neusp_buffer(1:eusp_size))
#else
  do ne = 1,size(edgeusp)
    nsz(1) = size(edgeusp(ne)%v,dim=1)
    nsz(2) = size(edgeusp(ne)%v,dim=2)
    ibeg = edgeusp(ne)%beg_pt
    iend = ibeg-1 + product(nsz)
    edgeusp(ne)%v(:,:) = reshape( real(neusp_buffer(ibeg:iend),kind=wp) , nsz )
  end do
#endif
  !
end subroutine exchange_edgeusp
!
!###############################################################################
!
subroutine exchange_edgedusp(edgedusp)
  !
  !.. Use Statements ..
  use flowvar, only : ne3d_t
  !
  !.. Formal Arguments ..
  type(ne3d_t), dimension(:), intent(inout) :: edgedusp
  !
  !.. Local Scalars..
  integer :: ne,vsz,ibeg,iend
  !
  !.. Local Arrays ..
  integer, dimension(1:3) :: nsz
  !
continue
  !
  ! Initialize neusp_buffer to zero
  !
  neusp_buffer = 0.0_r8
  !
#ifdef USE_TYPEBOUND_FACE_PACK
  call edgedusp%pack(neusp_buffer(1:edusp_size))
#else
  do ne = 1,size(edgedusp)
    vsz = size(edgedusp(ne)%v(:,:,:))
    ibeg = edgedusp(ne)%beg_pt
    iend = ibeg-1 + vsz
    neusp_buffer(ibeg:iend) = &
                 reshape( real(edgedusp(ne)%v(:,:,:),kind=r8) , [vsz] )
  end do
#endif
  !
  ! Sum all the edge point values from all the processors
  !
  call mpi_allreduce(MPI_IN_PLACE,neusp_buffer(1:edusp_size),edusp_size, &
                     mpi_flttyp,MPI_SUM,mpi_comm_world,mpierr)
  !
  ! Update the edge point values using the reduced values
  !
#ifdef USE_TYPEBOUND_FACE_PACK
  call edgedusp%unpack(neusp_buffer(1:edusp_size))
#else
  do ne = 1,size(edgedusp)
    nsz(1) = size(edgedusp(ne)%v,dim=1)
    nsz(2) = size(edgedusp(ne)%v,dim=2)
    nsz(3) = size(edgedusp(ne)%v,dim=3)
    ibeg = edgedusp(ne)%beg_pt
    iend = ibeg-1 + product(nsz)
    edgedusp(ne)%v(:,:,:) = &
                 reshape( real(neusp_buffer(ibeg:iend),kind=wp) , nsz )
  end do
#endif
  !
end subroutine exchange_edgedusp
!
!###############################################################################
!
subroutine exchange_nodeusp(nodeusp)
  !
  !.. Use Statements ..
  use flowvar, only : ne1d_t
  !
  !.. Formal Arguments ..
  type(ne1d_t), dimension(:), intent(inout) :: nodeusp
  !
  !.. Local Scalars..
  integer :: nn,vsz,ibeg,iend
  !
  !.. Local Arrays ..
  integer, dimension(1) :: nsz
  !
continue
  !
  ! Initialize neusp_buffer to zero
  !
  neusp_buffer = 0.0_r8
  !
#ifdef USE_TYPEBOUND_FACE_PACK
  call nodeusp%pack(neusp_buffer(1:nusp_size))
#else
  do nn = 1,size(nodeusp)
    vsz = size(nodeusp(nn)%v(:))
    ibeg = nodeusp(nn)%beg_pt
    iend = ibeg-1 + vsz
    neusp_buffer(ibeg:iend) = reshape( real(nodeusp(nn)%v(:),kind=r8) , [vsz] )
  end do
#endif
  !
  ! Sum all the node point values from all the processors
  !
  call mpi_allreduce(MPI_IN_PLACE,neusp_buffer(1:nusp_size),nusp_size, &
                     mpi_flttyp,MPI_SUM,mpi_comm_world,mpierr)
  !
  ! Update the right states for only the communication boundary faces
  !
#ifdef USE_TYPEBOUND_FACE_PACK
  call nodeusp%unpack(neusp_buffer(1:nusp_size))
#else
  do nn = 1,size(nodeusp)
    nsz(1) = size(nodeusp(nn)%v,dim=1)
    ibeg = nodeusp(nn)%beg_pt
    iend = ibeg-1 + product(nsz)
    nodeusp(nn)%v(:) = reshape( real(neusp_buffer(ibeg:iend),kind=wp) , nsz )
  end do
#endif
  !
end subroutine exchange_nodeusp
!
!###############################################################################
!
subroutine exchange_nodedusp(nodedusp)
  !
  !.. Use Statements ..
  use flowvar, only : ne2d_t
  !
  !.. Formal Arguments ..
  type(ne2d_t), dimension(:), intent(inout) :: nodedusp
  !
  !.. Local Scalars..
  integer :: nn,vsz,ibeg,iend
  !
  !.. Local Arrays ..
  integer, dimension(1:2) :: nsz
  !
continue
  !
  ! Initialize neusp_buffer to zero
  !
  neusp_buffer = 0.0_r8
  !
#ifdef USE_TYPEBOUND_FACE_PACK
  call nodedusp%pack(neusp_buffer(1:ndusp_size))
#else
  do nn = 1,size(nodedusp)
    vsz = size(nodedusp(nn)%v(:,:))
    ibeg = nodedusp(nn)%beg_pt
    iend = ibeg-1 + vsz
    neusp_buffer(ibeg:iend) = &
                 reshape( real(nodedusp(nn)%v(:,:),kind=r8) , [vsz] )
  end do
#endif
  !
  ! Sum all the node point values from all the processors
  !
  call mpi_allreduce(MPI_IN_PLACE,neusp_buffer(1:ndusp_size),ndusp_size, &
                     mpi_flttyp,MPI_SUM,mpi_comm_world,mpierr)
  !
  ! Update the node point values using the reduced values
  !
#ifdef USE_TYPEBOUND_FACE_PACK
  call nodedusp%unpack(neusp_buffer(1:ndusp_size))
#else
  do nn = 1,size(nodedusp)
    nsz(1) = size(nodedusp(nn)%v,dim=1)
    nsz(2) = size(nodedusp(nn)%v,dim=2)
    ibeg = nodedusp(nn)%beg_pt
    iend = ibeg-1 + product(nsz)
    nodedusp(nn)%v(:,:) = &
                 reshape( real(neusp_buffer(ibeg:iend),kind=wp) , nsz )
  end do
#endif
  !
end subroutine exchange_nodedusp
!
!###############################################################################
!
!subroutine exchange_modes_r4(umodes)
!  !
!  !.. Formal Arguments ..
!  real(r4), dimension(:,:), intent(inout) :: umodes
!  !
!  !.. Local Scalars..
!  integer :: n
!  integer(int_mpi) :: s,r,dtype
!  !
!  !.. Local Arrays ..
!  real(r4), dimension(1:size(umodes,1),1:size(umodes,2)) :: exch_send
!  integer(int_mpi), dimension(1:2*ncomm) :: rqst
!  integer(int_mpi), dimension(1:mpi_status_size,1:2*ncomm) :: istat
!  !
!continue
!  !
!  exch_send = umodes
!  !
!  ! Initialize all the receives so each process
!  ! is ready once the send is posted
!  !
!  do n = 1,ncomm
!    s = int( exch_cpu(n) , kind=int_mpi )
!    dtype = recv_fdusp(n)
!    call mpi_irecv(umodes,1,dtype,s,s,mpi_comm_world,rqst(n),mpierr)
!  end do
!  !
!  ! Send out the adjacent cells
!  !
!  do n = 1,ncomm
!    r = int( exch_cpu(n) , kind=int_mpi )
!    dtype = send_fdusp(n)
!    call mpi_isend(exch_send,1,dtype,r,mypnum, &
!                   mpi_comm_world,rqst(n+ncomm),mpierr)
!  end do
!  !
!  ! Wait for all sends and receives to complete
!  !
!  call mpi_waitall(int(2*ncomm,kind=int_mpi),rqst,istat,mpierr)
!  !
!end subroutine exchange_modes_r4
!
!###############################################################################
!
!subroutine exchange_modes_r8(umodes)
!  !
!  !.. Formal Arguments ..
!  real(r8), dimension(:,:), intent(inout) :: umodes
!  !
!  !.. Local Scalars..
!  integer :: n
!  integer(int_mpi) :: s,r,dtype
!  !
!  !.. Local Arrays ..
!  real(r8), dimension(1:size(umodes,1),1:size(umodes,2)) :: exch_send
!  integer(int_mpi), dimension(1:2*ncomm) :: rqst
!  integer(int_mpi), dimension(1:mpi_status_size,1:2*ncomm) :: istat
!  !
!continue
!  !
!  exch_send = umodes
!  !
!  ! Initialize all the receives so each process
!  ! is ready once the send is posted
!  !
!  do n = 1,ncomm
!    s = int( exch_cpu(n) , kind=int_mpi )
!    dtype = recv_fdusp(n)
!    call mpi_irecv(umodes,1,dtype,s,s,mpi_comm_world,rqst(n),mpierr)
!  end do
!  !
!  ! Send out the adjacent cells
!  !
!  do n = 1,ncomm
!    r = int( exch_cpu(n) , kind=int_mpi )
!    dtype = send_fdusp(n)
!    call mpi_isend(exch_send,1,dtype,r,mypnum, &
!                   mpi_comm_world,rqst(n+ncomm),mpierr)
!  end do
!  !
!  ! Wait for all sends and receives to complete
!  !
!  call mpi_waitall(int(2*ncomm,kind=int_mpi),rqst,istat,mpierr)
!  !
!end subroutine exchange_modes_r8
!
!###############################################################################
!
subroutine collect_global_solpts
  !
  !.. Use Statements ..
  use geovar, only : n_solpts,xyz,global_xyz
  !
  !.. Local Scalars ..
  integer :: i,icpu,n,n1,n2,ierr
  !
  !.. Local Arrays ..
  integer, dimension(1:mpi_status_size) :: istat
  !
  !.. Local Allocatable Arrays ..
  real(r8), allocatable, dimension(:,:) :: xyz_r8
  real(r8), allocatable, dimension(:,:) :: global_xyz_r8
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "collect_global_solpts"
  !
continue
  !
  if (ncpu > 1) then
    !
    n1 = size(xyz,dim=1)
    n2 = size(xyz,dim=2)
    !
    allocate ( xyz_r8(1:n1,1:n2) , source=real(xyz,r8) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_r8",1,ierr,error_message)
    !
    n1 = size(global_xyz,dim=1)
    n2 = size(global_xyz,dim=2)
    !
    allocate ( global_xyz_r8(1:n1,1:n2) , source=0.0_r8 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_xyz_r8",1,ierr,error_message)
    !
    n = size(xyz(:,1:n_solpts))
    !
    if (mypnum /= 0) then
      !
      call mpi_send(xyz_r8(:,1:n_solpts),n,mpi_flttyp,0, &
                    mypnum,mpi_comm_world,mpierr)
      !
    else
      !
      call mpi_sendrecv(xyz_r8(:,1:n_solpts),n,mpi_flttyp,0,n, &
                        global_xyz_r8,1,recv_xyz(1),0,n, &
                        mpi_comm_world,istat,mpierr)
      !
      do i = 2,ncpu
        !
        icpu = i-1
        !
        call mpi_recv(global_xyz_r8,1,recv_xyz(i),icpu,icpu, &
                      mpi_comm_world,istat,mpierr)
        !
      end do
      !
    end if
    !
    call mpi_bcast(global_xyz_r8,size(global_xyz_r8),mpi_flttyp, &
                   0,mpi_comm_world,mpierr)
    !
    global_xyz(:,:) = real( global_xyz_r8(:,:) , kind=wp )
    !
    deallocate ( xyz_r8 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_r8",2,ierr,error_message)
    deallocate ( global_xyz_r8 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_xyz_r8",2,ierr,error_message)
    !
  else
    !
    global_xyz(:,:) = xyz(:,:)
    !
  end if
  !
end subroutine collect_global_solpts
!
!###############################################################################
!
subroutine collect_global_solution
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : usp,global_usp
  !
  !.. Local Scalars ..
  integer :: i,icpu,n,n1,n2,ierr
  !
  !.. Local Arrays ..
  integer, dimension(1:mpi_status_size) :: istat
  !
  !.. Local Allocatable Arrays ..
  real(r8), allocatable, dimension(:,:) :: sp_r8
  real(r8), allocatable, dimension(:,:) :: global_r8
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "collect_global_solution"
  !
continue
  !
  if (ncpu > 1) then
    !
    n1 = size(usp,dim=1)
    n2 = size(usp,dim=2)
    !
    allocate ( sp_r8(1:n1,1:n2) , source=real(usp,kind=r8) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"sp_r8",1,ierr,error_message)
    !
    ! Collect the global usp array on the root processor
    !
    n = size(sp_r8(:,1:n_solpts))
    !
    if (mypnum /= 0) then
      !
      ! Have all the non-root processors send their
      ! local solutions to the root processor.
      !
      call mpi_send(sp_r8(:,1:n_solpts),n,mpi_flttyp,0, &
                    mypnum,mpi_comm_world,mpierr)
      !
    else
      !
      ! Only collect the global solution on the root processor. It is the
      ! only processor that will be writing the restart file so there is
      ! no need for the other processors to also have a copy of the global
      ! solution.
      !
      n1 = size(global_usp,dim=1)
      n2 = size(global_usp,dim=2)
      !
      allocate ( global_r8(1:n1,1:n2) , source=0.0_r8 , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"global_r8",1,ierr,error_message)
      !
      ! Have the root processor send its own solution to itself so
      ! it gets put in the right location within the global array.
      !
      call mpi_sendrecv(sp_r8(:,1:n_solpts),n,mpi_flttyp,0,n, &
                        global_r8,1,recv_sol(1),0,n, &
                        mpi_comm_world,istat,mpierr)
      !
      ! Now have the root processor receive the local
      ! solutions from all the other processors.
      !
      do i = 2,ncpu
        !
        icpu = i-1
        !
        call mpi_recv(global_r8,1,recv_sol(i),icpu,icpu, &
                      mpi_comm_world,istat,mpierr)
        !
      end do
      !
      ! Convert the global_r8 array to working precision
      !
      global_usp(:,:) = real( global_r8(:,:) , kind=wp )
      !
      ! Deallocate the global_r8 array
      !
      deallocate ( global_r8 , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"global_r8",2,ierr,error_message)
      !
    end if
    !
    deallocate ( sp_r8 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"sp_r8",2,ierr,error_message)
    !
  else
    !
    global_usp(:,:) = usp(:,:)
    !
  end if
  !
end subroutine collect_global_solution
!
!###############################################################################
!
subroutine collect_global_variable(var)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts,n_global_solpts
  use flowvar, only : global_var
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: var
  !
  !.. Local Scalars ..
  integer :: i,icpu,npts,ierr
  !
  !.. Local Arrays ..
  integer, dimension(1:mpi_status_size) :: istat
  !
  !.. Local Allocatable Arrays ..
  real(r8), allocatable, dimension(:) :: sp_r8
  real(r8), allocatable, dimension(:) :: global_r8
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "collect_global_variable"
  !
continue
  !
  if (.not. allocated(global_var)) then
    allocate ( global_var(1:n_global_solpts) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_var",1,ierr,error_message)
  end if
  !
  if (ncpu > 1) then
    !
    allocate ( sp_r8(1:size(var)) , source=real(var,r8) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"sp_r8",1,ierr,error_message)
    allocate ( global_r8(1:size(global_var)) , source=0.0_r8 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_r8",1,ierr,error_message)
    !
    ! Collect the global variable array on all processors
    !
    npts = n_solpts
    !
    if (mypnum /= 0) then
      !
      call mpi_send(sp_r8(1:npts),npts,mpi_flttyp,0, &
                    mypnum,mpi_comm_world,mpierr)
      !
    else
      !
      call mpi_sendrecv(sp_r8(1:n_solpts),npts,mpi_flttyp,0,npts, &
                        global_r8,1,recv_var(1),0,npts, &
                        mpi_comm_world,istat,mpierr)
      !
      do i = 2,ncpu
        !
        icpu = i-1
        !
        call mpi_recv(global_r8,1,recv_var(i),icpu,icpu, &
                      mpi_comm_world,istat,mpierr)
        !
      end do
      !
    end if
    !
    call mpi_bcast(global_r8,size(global_r8),mpi_flttyp,0,mpi_comm_world,mpierr)
    !
    global_var(:) = real( global_r8(:) , kind=wp )
    !
    deallocate ( sp_r8 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"sp_r8",2,ierr,error_message)
    deallocate ( global_r8 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_r8",2,ierr,error_message)
    !
  else
    !
    global_var(:) = var(:)
    !
  end if
  !
end subroutine collect_global_variable
!
!###############################################################################
!
subroutine output_partitions(id,nodes_of_cell,nodes_of_cell_ptr,xyz_nodes,epart)
  !
  !.. Formal Arguments ..
  integer,                          intent(in) :: id
  integer,  dimension(:),           intent(in) :: nodes_of_cell
  integer,  dimension(:),           intent(in) :: nodes_of_cell_ptr
  real(wp), dimension(:,:),         intent(in) :: xyz_nodes
  integer,  dimension(:), optional, intent(in) :: epart
  !
  !.. Local Scalars ..
  integer :: i,n,n1,n2,ierr,np,nc,ndim
  character(len=100) :: fname
  character(len=20)  :: zonetype
  character(len=50)  :: variables
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:,:) :: ipelem
  !
  !.. Local Constants ..
  integer, parameter :: iot = 103
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_partitions"
  !
continue
  !
  ndim = size(xyz_nodes,dim=1)
  !
  zonetype = "FEQUADRILATERAL"
  variables = 'VARIABLES = "X", "Y"'
  !
  if (ndim == 3) then
    variables = trim(adjustl(variables)) // ', "Z"'
    zonetype = "FEBRICK"
  end if
  !
  if (present(epart)) then
    fname = "partitions.tec"
    np = ngnode
    nc = ngcell
    variables = trim(adjustl(variables)) // ', "Partition"'
  else
    write (fname,'("partitions.",i0,".tec")') id
    np = size(xyz_nodes,dim=2)
    nc = size(nodes_of_cell_ptr)-1
  end if
  fname = adjustl(fname)
  !
  if (ndim == 3) then
    !
    allocate ( ipelem(1:8,1:nc) , source=0 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ipelem",1,ierr,error_message)
    !
    do n = 1,nc
      !
      n1 = nodes_of_cell_ptr(n)+1
      n2 = nodes_of_cell_ptr(n+1)
      !
      if (n2-n1+1 == 4) then
        !
        ipelem(1,n) = nodes_of_cell(n1+0)
        ipelem(2,n) = nodes_of_cell(n1+1)
        ipelem(3,n) = nodes_of_cell(n1+2)
        ipelem(4,n) = nodes_of_cell(n1+2)
        ipelem(5,n) = nodes_of_cell(n1+3)
        ipelem(6,n) = nodes_of_cell(n1+3)
        ipelem(7,n) = nodes_of_cell(n1+3)
        ipelem(8,n) = nodes_of_cell(n1+3)
        !
      else if (n2-n1+1 == 5) then
        !
        ipelem(1,n) = nodes_of_cell(n1+0)
        ipelem(2,n) = nodes_of_cell(n1+1)
        ipelem(3,n) = nodes_of_cell(n1+2)
        ipelem(4,n) = nodes_of_cell(n1+3)
        ipelem(5,n) = nodes_of_cell(n1+4)
        ipelem(6,n) = nodes_of_cell(n1+4)
        ipelem(7,n) = nodes_of_cell(n1+4)
        ipelem(8,n) = nodes_of_cell(n1+4)
        !
      else if (n2-n1+1 == 6) then
        !
        ipelem(1,n) = nodes_of_cell(n1+0)
        ipelem(2,n) = nodes_of_cell(n1+1)
        ipelem(3,n) = nodes_of_cell(n1+2)
        ipelem(4,n) = nodes_of_cell(n1+2)
        ipelem(5,n) = nodes_of_cell(n1+3)
        ipelem(6,n) = nodes_of_cell(n1+4)
        ipelem(7,n) = nodes_of_cell(n1+5)
        ipelem(8,n) = nodes_of_cell(n1+5)
        !
      else if (n2-n1+1 == 8) then
        !
        ipelem(1,n) = nodes_of_cell(n1+0)
        ipelem(2,n) = nodes_of_cell(n1+1)
        ipelem(3,n) = nodes_of_cell(n1+2)
        ipelem(4,n) = nodes_of_cell(n1+3)
        ipelem(5,n) = nodes_of_cell(n1+4)
        ipelem(6,n) = nodes_of_cell(n1+5)
        ipelem(7,n) = nodes_of_cell(n1+6)
        ipelem(8,n) = nodes_of_cell(n1+7)
        !
      end if
      !
    end do
    !
  else
    !
    allocate ( ipelem(1:4,1:nc) , source=0 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ipelem",1,ierr,error_message)
    !
    do n = 1,nc
      !
      n1 = nodes_of_cell_ptr(n)+1
      n2 = nodes_of_cell_ptr(n+1)
      !
      ipelem(1,n) = nodes_of_cell(n1  )
      ipelem(2,n) = nodes_of_cell(n1+1)
      ipelem(3,n) = nodes_of_cell(n1+2)
      ipelem(4,n) = nodes_of_cell(n2  )
      !
    end do
    !
  end if 
  !
  open (iot,file=trim(fname),status="replace",action="write", &
        iostat=ierr,iomsg=error_message)
  call io_error(pname,trim(fname),1,ierr,error_message)
  !
  ! Tecplot header with grid size parameters
  !
  write (iot,6) trim(adjustl(variables)),np,nc,trim(adjustl(zonetype))
  if (present(epart)) write (iot,7) ndim + 1
  !
  ! Node coordinates
  !
  do n = 1,size(xyz_nodes,dim=1)
    write (iot,12) (xyz_nodes(n,i),i=1,np)
  end do
  !
  ! Cell partition 
  !
  if (present(epart)) write (iot,13) (epart(i),i=1,nc)
  !
  ! Cell Connectivity
  !
  do n = 1,nc
    write (iot,13) (ipelem(i,n), i=1,size(ipelem,dim=1))
  end do
  !
  close (iot,iostat=ierr,iomsg=error_message)
  call io_error(pname,trim(fname),2,ierr,error_message)
  !
  if (allocated(ipelem)) then
    deallocate ( ipelem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ipelem",2,ierr,error_message)
  end if
  !
 !1 format (a,/,)
 !2 format ('VARIABLES = "X", "Y", "Partition"')
 !3 format ('VARIABLES = "X", "Y", "Z"')
 !4 format ('VARIABLES = "X", "Y", "Z", "Partition"')
  5 format ('VARIABLES = "X" "Y"',/, &
            'ZONE T="Partitioned Grid", N=',i0,', E=',i0, &
            ', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL')
  1 format ('VARIABLES = "X" "Y" "Partition"',/, &
            'ZONE T="Partitioned Grid", N=',i0,', E=',i0, &
            ', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL',/, &
            'VARLOCATION = ([3]=CellCentered)')
  6 format (a,/,'ZONE T="Partitioned Grid", N=',i0,', E=',i0, &
            ', DATAPACKING=BLOCK, ZONETYPE=',a:,/,a)
  7 format ("VARLOCATION = ([",i0,"]=CellCentered)")
  12 format (6(1x,es15.8))
  13 format (8(4x,i0))
  14 format (4(3x,i0))
  !
end subroutine output_partitions
!
!###############################################################################
!
subroutine parallel_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
  !
  !.. Formal Arguments ..
  integer, optional, intent(in) :: iunit
  !
  !.. Local Scalars ..
  integer :: iou,i
  integer :: loc_siz,glb_siz
  character(len=36) :: array_name
  type(memory) :: array,dt
  type(memory) :: total
  !
continue
  !
  iou = iout
  if (present(iunit)) iou = iunit
  !
  call total%reset
  !
  write (iou,1)
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(cell_map)) then
    !
    call dt%set( storage_size(cell_map,kind=inttype) , &
                         size(cell_map,kind=inttype) )
    !
    do i = lbound(cell_map,dim=1),ubound(cell_map,dim=1)
      if (allocated(cell_map(i)%loc_to_glb)) then
        call dt%add( storage_size(cell_map(i)%loc_to_glb,kind=inttype) , &
                             size(cell_map(i)%loc_to_glb,kind=inttype) )
      end if
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "cell_map"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "cell_map"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
!!  if (allocated(node_map)) then
!!    !
!!    call dt%set( storage_size(node_map,kind=inttype) , &
!!                         size(node_map,kind=inttype) )
!!   !write (iout,*)
!!   !write (iout,*) "node_map(:)"
!!   !write (iout,*) "size = ",size(node_map)
!!   !write (iout,*) "storage_size = ",storage_size(node_map)
!!   !!
!!   !loc_siz = 0
!!   !glb_siz = 0
!!    !
!!    do i = lbound(node_map,dim=1),ubound(node_map,dim=1)
!!      !
!!      if (allocated(node_map(i)%loc_to_glb)) then
!!        call dt%add( storage_size(node_map(i)%loc_to_glb,kind=inttype) , &
!!                             size(node_map(i)%loc_to_glb,kind=inttype) )
!!   !    loc_siz = loc_siz + size(node_map(i)%loc_to_glb)
!!      end if
!!      !
!!    end do
!!    !
!!   !i = lbound(node_map,dim=1)
!!   !if (allocated(node_map(i)%loc_to_glb)) then
!!   !  write (iout,*)
!!   !  write (iout,*) "node_map(:)%loc_to_glb"
!!   !  write (iout,*) "size = ",loc_siz
!!   !  write (iout,*) "storage_size = ",storage_size(node_map(i)%loc_to_glb)
!!   !end if
!!    !
!!    call total%add(dt)
!!    write (array_name,5) "node_map"
!!    write (iou,2) array_name,dt%get_units()
!!    !
!!  else
!!    !
!!    write (iou,7) "node_map"
!!    !
!!  end if
!!  !
!!  !-----------------------------------------------------------------------------
!!  !
!!  if (allocated(face_map)) then
!!    !
!!    call dt%set( storage_size(face_map,kind=inttype) , &
!!                         size(face_map,kind=inttype) )
!!   !write (iout,*)
!!   !write (iout,*) "face_map(:)"
!!   !write (iout,*) "size = ",size(face_map)
!!   !write (iout,*) "storage_size = ",storage_size(face_map)
!!   !!
!!   !loc_siz = 0
!!   !glb_siz = 0
!!    !
!!    do i = lbound(face_map,dim=1),ubound(face_map,dim=1)
!!      !
!!      if (allocated(face_map(i)%loc_to_glb)) then
!!        call dt%add( storage_size(face_map(i)%loc_to_glb,kind=inttype) , &
!!                             size(face_map(i)%loc_to_glb,kind=inttype) )
!!   !    loc_siz = loc_siz + size(face_map(i)%loc_to_glb)
!!      end if
!!      !
!!    end do
!!    !
!!   !i = lbound(face_map,dim=1)
!!   !if (allocated(face_map(i)%loc_to_glb)) then
!!   !  write (iout,*)
!!   !  write (iout,*) "face_map(:)%loc_to_glb"
!!   !  write (iout,*) "size = ",loc_siz
!!   !  write (iout,*) "storage_size = ",storage_size(face_map(i)%loc_to_glb)
!!   !end if
!!    !
!!    call total%add(dt)
!!    write (array_name,5) "face_map"
!!    write (iou,2) array_name,dt%get_units()
!!    !
!!  else
!!    !
!!    write (iou,7) "face_map"
!!    !
!!  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(face)) then
    !
    call dt%set( storage_size(face,kind=inttype) , &
                         size(face,kind=inttype) )
    !
    do i = lbound(face,dim=1),ubound(face,dim=1)
      !
      call dt%add( storage_size(face(i)%geom,kind=inttype) )
      call dt%add( storage_size(face(i)%order,kind=inttype) )
      !
      call dt%add( storage_size(face(i)%nodes,kind=inttype) , &
                           size(face(i)%nodes,kind=inttype) )
      !
      call dt%add( storage_size(face(i)%left,kind=inttype) )
      call dt%add( storage_size(face(i)%left%cell,kind=inttype) )
      call dt%add( storage_size(face(i)%left%cell_face,kind=inttype) )
      call dt%add( storage_size(face(i)%left%geom,kind=inttype) )
      call dt%add( storage_size(face(i)%left%order,kind=inttype) )
      call dt%add( storage_size(face(i)%left%partition,kind=inttype) )
      !
      if (allocated(face(i)%left%flxpts)) then
        call dt%add( storage_size(face(i)%left%flxpts,kind=inttype) , &
                             size(face(i)%left%flxpts,kind=inttype) )
      end if
      !
      call dt%add( storage_size(face(i)%right,kind=inttype) )
      call dt%add( storage_size(face(i)%right%cell,kind=inttype) )
      call dt%add( storage_size(face(i)%right%cell_face,kind=inttype) )
      call dt%add( storage_size(face(i)%right%geom,kind=inttype) )
      call dt%add( storage_size(face(i)%right%order,kind=inttype) )
      call dt%add( storage_size(face(i)%right%partition,kind=inttype) )
      !
      if (allocated(face(i)%right%flxpts)) then
        call dt%add( storage_size(face(i)%right%flxpts,kind=inttype) , &
                             size(face(i)%right%flxpts,kind=inttype) )
      end if
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "face"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "face"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(exch_fusp)) then
    !
    call dt%set( storage_size(exch_fusp,kind=inttype) , &
                         size(exch_fusp,kind=inttype) )
    !
    do i = lbound(exch_fusp,dim=1),ubound(exch_fusp,dim=1)
      !
      call dt%add( storage_size(exch_fusp(i)%adj_cpu,kind=inttype) )
      call dt%add( storage_size(exch_fusp(i)%send_tag,kind=inttype) )
      call dt%add( storage_size(exch_fusp(i)%recv_tag,kind=inttype) )
      call dt%add( storage_size(exch_fusp(i)%send_datatype,kind=inttype) )
      call dt%add( storage_size(exch_fusp(i)%recv_datatype,kind=inttype) )
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "exch_fusp"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "exch_fusp"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(exch_fdusp)) then
    !
    call dt%set( storage_size(exch_fdusp,kind=inttype) , &
                         size(exch_fdusp,kind=inttype) )
    !
    do i = lbound(exch_fdusp,dim=1),ubound(exch_fdusp,dim=1)
      !
      call dt%add( storage_size(exch_fdusp(i)%adj_cpu,kind=inttype) )
      call dt%add( storage_size(exch_fdusp(i)%send_tag,kind=inttype) )
      call dt%add( storage_size(exch_fdusp(i)%recv_tag,kind=inttype) )
      call dt%add( storage_size(exch_fdusp(i)%send_datatype,kind=inttype) )
      call dt%add( storage_size(exch_fdusp(i)%recv_datatype,kind=inttype) )
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "exch_fdusp"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "exch_fdusp"
    !
  end if
  !
  !
  ! Check the memory of the arrays within this module
  !
  if (allocated(exch_cpu)) then
    call array%set( storage_size(exch_cpu,kind=inttype) , &
                            size(exch_cpu,kind=inttype) )
    call total%add(array)
    write (array_name,5) "exch_cpu"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "exch_cpu"
  end if
  if (allocated(testf_send)) then
    call array%set( storage_size(testf_send,kind=inttype) , &
                            size(testf_send,kind=inttype) )
    call total%add(array)
    write (array_name,5) "testf_send"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "testf_send"
  end if
  if (allocated(testf_recv)) then
    call array%set( storage_size(testf_recv,kind=inttype) , &
                            size(testf_recv,kind=inttype) )
    call total%add(array)
    write (array_name,5) "testf_recv"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "testf_recv"
  end if
  if (allocated(fusp_rqst)) then
    call array%set( storage_size(fusp_rqst,kind=inttype) , &
                            size(fusp_rqst,kind=inttype) )
    call total%add(array)
    write (array_name,5) "fusp_rqst"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "fusp_rqst"
  end if
  if (allocated(fusp_istat)) then
    call array%set( storage_size(fusp_istat,kind=inttype) , &
                            size(fusp_istat,kind=inttype) )
    call total%add(array)
    write (array_name,5) "fusp_istat"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "fusp_istat"
  end if
  if (allocated(fdusp_rqst)) then
    call array%set( storage_size(fdusp_rqst,kind=inttype) , &
                            size(fdusp_rqst,kind=inttype) )
    call total%add(array)
    write (array_name,5) "fdusp_rqst"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "fdusp_rqst"
  end if
  if (allocated(fdusp_istat)) then
    call array%set( storage_size(fdusp_istat,kind=inttype) , &
                            size(fdusp_istat,kind=inttype) )
    call total%add(array)
    write (array_name,5) "fdusp_istat"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "fdusp_istat"
  end if
  if (allocated(fusp_send_buffer)) then
    call array%set( storage_size(fusp_send_buffer,kind=inttype) , &
                            size(fusp_send_buffer,kind=inttype) )
    call total%add(array)
    write (array_name,5) "fusp_send_buffer"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "fusp_send_buffer"
  end if
  if (allocated(fusp_recv_buffer)) then
    call array%set( storage_size(fusp_recv_buffer,kind=inttype) , &
                            size(fusp_recv_buffer,kind=inttype) )
    call total%add(array)
    write (array_name,5) "fusp_recv_buffer"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "fusp_recv_buffer"
  end if
  if (allocated(fdusp_send_buffer)) then
    call array%set( storage_size(fdusp_send_buffer,kind=inttype) , &
                            size(fdusp_send_buffer,kind=inttype) )
    call total%add(array)
    write (array_name,5) "fdusp_send_buffer"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "fdusp_send_buffer"
  end if
  if (allocated(fdusp_recv_buffer)) then
    call array%set( storage_size(fdusp_recv_buffer,kind=inttype) , &
                            size(fdusp_recv_buffer,kind=inttype) )
    call total%add(array)
    write (array_name,5) "fdusp_recv_buffer"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "fdusp_recv_buffer"
  end if
  if (allocated(recv_xyz)) then
    call array%set( storage_size(recv_xyz,kind=inttype) , &
                            size(recv_xyz,kind=inttype) )
    call total%add(array)
    write (array_name,5) "recv_xyz"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "recv_xyz"
  end if
  if (allocated(recv_sol)) then
    call array%set( storage_size(recv_sol,kind=inttype) , &
                            size(recv_sol,kind=inttype) )
    call total%add(array)
    write (array_name,5) "recv_sol"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "recv_sol"
  end if
  if (allocated(recv_var)) then
    call array%set( storage_size(recv_var,kind=inttype) , &
                            size(recv_var,kind=inttype) )
    call total%add(array)
    write (array_name,5) "recv_var"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "recv_var"
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module parallel_mod")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  4 format (a,"(",i0,")%",a)
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
end subroutine parallel_memory_usage
!
!###############################################################################
!
include "Functions/last.f90"
include "Functions/cell_solpts.f90"
!
!###############################################################################
!
end module parallel_mod
