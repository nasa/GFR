module parallel_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nq
  use geovar,  only : face_t
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
  public :: exch_connectivity
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
  integer, save :: n_node_edge_pts
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
  ! ######################################################################
  ! #####   BOUNDARY NODE/EDGE DERIVED TYPES AND ASSOCIATED ARRAYS   #####
  ! ######################################################################
  !
  ! BND_DATA_T : derived type to compactly store the local indices of either
  !              the nodes or edges of a partition that are located on a
  !              partition boundary
  !
  type :: bnd_data_t
    ! idx : local index of the node or edge on a partition boundary
    integer :: idx
    ! ave_correction : the averaging coefficient used to compute the local node
    !                  or edge point values is based on the number of cells on
    !                  the local partition that contribute to this node/edge.
    !                  If the number of local cells containing this node/edge
    !                  is lc, these values are locally averaged using the
    !                  coefficient (1/lc). If the number of global cells
    !                  containing this node/edge is gc, this averaging
    !                  correction is simply multiplying these node/edge
    !                  values by (lc/gc), resulting in the final averaging
    !                  coefficient of (1/gc) for this node/edge.
    real(wp) :: ave_correction
    ! my_responsibility : logical flag that indicates if the node or edge point
    !                     indices for the subcell connectivities involving this
    !                     this node or edge are the responsibility of this
    !                     processor
    logical :: my_responsibility = fals
  end type bnd_data_t
  !
  ! BND_NODES : Array of type bnd_data_t to store the local index and
  !             averaging correction for all nodes of the partition that
  !             are located on a partition boundary
  !
  type(bnd_data_t), save, allocatable, dimension(:) :: bnd_nodes
  !
  ! BND_EDGES : Array of type bnd_data_t to store the local index and
  !             averaging correction for all edges of the partition that
  !             are located on a partition boundary
  !
  type(bnd_data_t), save, allocatable, dimension(:) :: bnd_edges
  !
  ! BND_FACES : Array of type bnd_data_t to store the local index and
  !             averaging correction for all faces of the partition that
  !             are located on a partition boundary
  !             NOTE: The averaging correction is not needed for faces as
  !                   this is only really used to identify the partition
  !                   that is responsible for the face point indices
  !
  type(bnd_data_t), save, allocatable, dimension(:) :: bnd_faces
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
  ! #####################################################################
  ! #####   EDGE CONNECTIVITY DERIVED TYPES AND ASSOCIATED ARRAYS   #####
  ! #####################################################################
  !
  ! CPU_WITH_EDGE_T : derived type to store the information about
  !                   a partition/processor that contains a given edge
  !
  type :: cpu_with_edge_t
    ! partition : grid partition/processor containing this edge
    integer :: partition
    ! loc_edge : index of this edge local this partition
    integer :: loc_edge
    ! num_cells : number of cells on this partition containing this edge
    integer :: num_cells
    ! edgpts : indices for the edge points on the current edge
    integer, allocatable :: edgpts(:)
  end type cpu_with_edge_t
  !
  logical, parameter :: use_edgpts = true
 !logical, parameter :: use_edgpts = fals
  !
  ! EDGE_T : Derived type to store information for each grid edge
  !          NOTE: This is different from the derived type
  !                of the same name in the module geovar
  !
  type :: edge_t
    ! is_on_boundary : logical that identifies if this edge is on a
    !                  non-communication boundary face
    logical :: is_on_boundary = fals
    ! controlling_partition : partition that is responsible for the index
    !                         of this node for the CGNS connectivity
    integer :: controlling_partition
    ! order : solution order of the current edge
    integer :: order
    ! cpu : information for each cpu that contains the current edge
    type(cpu_with_edge_t), allocatable :: cpu(:)
  end type edge_t
  !
  !
  !
  ! #####################################################################
  ! #####   NODE CONNECTIVITY DERIVED TYPES AND ASSOCIATED ARRAYS   #####
  ! #####################################################################
  !
  ! CPU_WITH_NODE_T : derived type to store the information about
  !                   a partition/processor that contains a given node
  !
  type :: cpu_with_node_t
    ! partition : grid partition/processor containing this node
    integer :: partition
    ! loc_node : index of this node local this partition
    integer :: loc_node
    ! num_cells : number of cells on this partition containing this edge
    integer :: num_cells
  end type cpu_with_node_t
  !
  ! NODE_T : Derived type to store information for each grid node
  !          NOTE: This is different from the derived type
  !                of the same name in the module geovar
  !
  type :: node_t
    ! is_on_boundary : logical that identifies if this edge is on a
    !                  non-communication boundary face
    logical :: is_on_boundary = fals
    ! controlling_partition : partition that is responsible for the index
    !                         of this node for the CGNS connectivity
    integer :: controlling_partition
    ! cpu : information for each cpu that contains the current node
    type(cpu_with_node_t), allocatable :: cpu(:)
  end type node_t
  !
  !
  !
  ! #####################################################################
  ! #####   MPI COMMUNICATION DERIVED TYPES AND ASSOCIATED ARRAYS   #####
  ! #####################################################################
  !
  ! EXCH_DATA_T : Derived type containing the information needed to communicate
  !               data between two processes
  !               NOTE: adj_cpu, send_tag, and recv_tag all default to -1
  !                     because this represents an invalid value if these
  !                     components are used without explicitly redefining
  !                     them to a valid value
  !
  type :: exch_data_t
    ! adj_cpu : process/CPU ID for the communication partner
    integer(int_mpi) :: adj_cpu = -1_int_mpi
    ! send_tag : message tag for matching the MPI communication
    integer(int_mpi) :: send_tag = -1_int_mpi
    ! recv_tag : message tag for matching the MPI communication
    integer(int_mpi) :: recv_tag = -1_int_mpi
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
  ! EXCH_NEUSP : Array containing the information for all the communication
  !              pairs that need to exchange edge and node solutions
  !              
  type(exch_data_t), save, allocatable :: exch_neusp(:)
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
subroutine create_serial_cell_map(lcell)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: lcell
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
  allocate ( cell_map(1)%loc_to_glb(1:lcell) , source=intseq(1,lcell) , &
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
  integer :: lcell,ierr
  character(len=200) :: array_name
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:)   :: xadj
  integer, allocatable, dimension(:)   :: adjncy
  integer, allocatable, dimension(:)   :: epart
  integer, allocatable, dimension(:)   :: fpart
  integer, allocatable, dimension(:)   :: cells_with_node
  integer, allocatable, dimension(:)   :: cells_with_node_ptr
  integer, allocatable, dimension(:)   :: cells_surr_cell
  integer, allocatable, dimension(:)   :: cells_surr_cell_ptr
  integer, allocatable, dimension(:)   :: cells_with_edge
  integer, allocatable, dimension(:)   :: cells_with_edge_ptr
  integer, allocatable, dimension(:,:) :: nodes_on_edge
  !
  !.. Local Allocatable Derived-Type Arrays ..
  type(node_t), allocatable, dimension(:) :: node
  type(map_t),  allocatable, dimension(:) :: node_map
  type(edge_t), allocatable, dimension(:) :: edge
  type(map_t),  allocatable, dimension(:) :: edge_map
  type(face_t), allocatable, dimension(:) :: face
  type(map_t),  allocatable, dimension(:) :: face_map
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "partition_grid"
  logical, parameter :: use_old_dual_graph = true
 !logical, parameter :: use_old_dual_graph = fals
  !
continue
  !
  ! Get the size of the global grid
  !
  nr = size(xyz_nodes,dim=1)
  !
  ! Get the number of interior grid cells
  !
  lcell = size(nodes_of_cell_ptr)-1 - size(bface,dim=2)
  !
  ! Create the global grid connectivity
  !
  call get_global_connectivity(cell_geom,cell_order,bface,xyz_nodes, &
                               nodes_of_cell,nodes_of_cell_ptr, &
                               cells_with_node,cells_with_node_ptr, &
                               cells_surr_cell,cells_surr_cell_ptr, &
                               cells_with_edge,cells_with_edge_ptr, &
                               nodes_on_edge,face)
  !
  if (use_old_dual_graph) then
    !
    ! Create the dual graph for the global grid not including ghost cells
    !
    call memory_pause("Before: calling create_dual_graph")
    call create_dual_graph(cell_geom,cell_order,nodes_of_cell, &
                           nodes_of_cell_ptr,bface,xadj,adjncy,face)
  else
    !
    ! Create the dual graph for the global grid not including ghost cells
    !
    call memory_pause("Before: calling create_metis_dual_graph")
    call create_metis_dual_graph(cells_surr_cell,cells_surr_cell_ptr, &
                                 xadj,adjncy)
    !
  end if
  !
  ! Deallocate cells_surr_cell and cells_surr_cell_ptr since they
  ! are no longer needed
  !
  if (allocated(cells_surr_cell)) then
    deallocate ( cells_surr_cell , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cells_surr_cell",2,ierr,error_message)
  end if
  !
  if (allocated(cells_surr_cell_ptr)) then
    deallocate ( cells_surr_cell_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cells_surr_cell_ptr",2,ierr,error_message)
  end if
  !
  ! Allocate the epart array to store the cell partition numbers
  !
  allocate ( epart(1:lcell) , source=0 , stat=ierr , errmsg=error_message )
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
                               cell_geom,cell_order, &
                               xadj,adjncy,epart)
  !
#elif METIS_5
  !
  call memory_pause("Before: calling partition_using_metis_5")
  call partition_using_metis_5(metis_option_requested,npart, &
                               xadj,adjncy,epart)
  !
#else
  !
  write (error_message,1)
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
  call mpi_bcast(epart,int(lcell,kind=kind(int_mpi)), &
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
  ! Create the node and edge arrays
  !
  call memory_pause("Before: calling create_node_and_edge_arrays")
  call create_node_and_edge_arrays(nodes_of_cell,nodes_of_cell_ptr, &
                                   cells_with_node,cells_with_node_ptr, &
                                   cells_with_edge,cells_with_edge_ptr, &
                                  !nodes_on_edge, &
                                   cell_order,epart,bface,node,edge)
  !
  ! Deallocate nodes_on_edge, cells_with_node, and cells_with_node_ptr
  ! since they are no longer needed
  !
  if (allocated(nodes_on_edge)) then
    deallocate ( nodes_on_edge , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodes_on_edge",2,ierr,error_message)
  end if
  !
  if (allocated(cells_with_node)) then
    deallocate ( cells_with_node , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cells_with_node",2,ierr,error_message)
  end if
  !
  if (allocated(cells_with_node_ptr)) then
    deallocate ( cells_with_node_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cells_with_node_ptr",2,ierr,error_message)
  end if
  !
  ! Find the faces that separate two adjacent partitions
  ! and add these faces to boundary faces array
  !
  call memory_pause("Before: calling find_partition_boundaries")
  call find_partition_boundaries(cell_geom,epart,bface,fpart,face)
  !
  ! Now that we have created boundary conditions for faces
  ! located on partition boundaries, use epart to create
  ! mappings from local to global block orderings.
  ! This also includes creating the reverse mapping arrays
  ! or the maps from global to local block orderings
  !
  call memory_pause("Before: calling create_grid_mappings")
  call create_grid_mappings(epart,face,node,edge, &
                            nodes_of_cell,nodes_of_cell_ptr, &
                            cells_with_edge,cells_with_edge_ptr, &
                            node_map,edge_map,face_map,cell_map)
  !
  ! Deallocate cells_with_edge and cells_with_edge_ptr since they
  ! are no longer needed
  !
  if (allocated(cells_with_edge)) then
    deallocate ( cells_with_edge , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cells_with_edge",2,ierr,error_message)
  end if
  !
  if (allocated(cells_with_edge_ptr)) then
    deallocate ( cells_with_edge_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cells_with_edge_ptr",2,ierr,error_message)
  end if
  !
 !call parallel_memory_usage(iout)
  !
  ! Now that the mappings have been created, we
  ! need to localize the grid on each processor
  !
  call memory_pause("Before: calling localize_grid")
  call localize_grid(nodes_of_cell,nodes_of_cell_ptr,bface, &
                     xyz_nodes,epart,cell_geom,cell_order, &
                     face,node_map,face_map,cell_map)
  !
  ! Deallocate epart since it is no longer needed
  !
  if (allocated(epart)) then
    deallocate ( epart , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"epart",2,ierr,error_message)
  end if
  !
  ! Create MPI datatypes for exchanging node and edge data between processes
  !
  call memory_pause("Before: calling create_node_and_edge_datatypes")
  call create_node_and_edge_datatypes(node,node_map,edge,edge_map)
  !
  ! Deallocate the node and edge arrays and the node and edge
  ! mapping arrays since they are no longer needed
  !
  if (allocated(node)) then
    deallocate ( node , stat=ierr ,errmsg=error_message )
    call alloc_error(pname,"node",2,ierr,error_message)
  end if
  !
  if (allocated(edge)) then
    deallocate ( edge , stat=ierr ,errmsg=error_message )
    call alloc_error(pname,"edge",2,ierr,error_message)
  end if
  !
  if (allocated(node_map)) then
    deallocate ( node_map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"node_map",2,ierr,error_message)
  end if
  !
  if (allocated(edge_map)) then
    deallocate ( edge_map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edge_map",2,ierr,error_message)
  end if
  !
  ! Finally, create MPI datatypes for exchanging face data between processes
  !
  call memory_pause("Before: calling create_face_nonblocking_datatypes")
  call create_face_nonblocking_datatypes(bface,fpart,face,face_map)
  !
  ! Deallocate face, face_map, and fpart since they are no longer needed
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
  if (allocated(fpart)) then
    deallocate ( fpart , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"fpart",2,ierr,error_message)
  end if
  !
  ! Format Statements
  !
  1 format(" It seems that the code was compiled without either of the METIS", &
           " pre-processor options being defined. One of these options needs", &
           " to be defined in order to partition the grid and allow the use", &
           " of multiple processors.")
  !
end subroutine partition_grid
!
!###############################################################################
!
subroutine partition_using_metis_4(metis_option_requested,npart, &
                                   cell_geom,cell_order, &
                                   xadj,adjncy,epart)
  !
  !.. Use Statements ..
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  !
  !.. Formal Arguments ..
  integer,                  intent(in) :: metis_option_requested
  integer,                  intent(in) :: npart
  integer, dimension(:),    intent(in) :: cell_geom
  integer, dimension(:),    intent(in) :: cell_order
  integer, dimension(:),    intent(in) :: xadj
  integer, dimension(:),    intent(in) :: adjncy
  integer, dimension(:), intent(inout) :: epart
  !
  !.. Local Scalars ..
  integer :: n,ip,ierr,lcell
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
  integer, allocatable, dimension(:) :: vsize
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "partition_using_metis_4"
  !
continue
#ifdef METIS_4
  !
  lcell = size(epart)
  !
  ! Create the priority list for METIS options to try
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
      ! Create the weights for the vertices of the dual graph
      !
      call create_grid_weights(metis_option,1,lcell,size(adjncy), &
                               vwgt,adjwgt,vsize,wgtflag,cell_geom,cell_order)
      !
      ! Define other METIS variables
      !
      numflag = 1
      options(:) = 0
      vol = 0
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
        call metis_partgraphrecursive(lcell,xadj,adjncy,vwgt,adjwgt,wgtflag, &
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
        call metis_partgraphkway(lcell,xadj,adjncy,vwgt,adjwgt,wgtflag, &
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
        call metis_partgraphvkway(lcell,xadj,adjncy,vwgt,vsize,wgtflag, &
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
        call metis_wpartgraphvkway(lcell,xadj,adjncy,vwgt,vsize,wgtflag, &
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
  if (allocated(vwgt)) then
    deallocate ( vwgt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vwgt",2,ierr,error_message)
  end if
  !
  if (allocated(adjwgt)) then
    deallocate ( adjwgt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"adjwgt",2,ierr,error_message)
  end if
  !
  if (allocated(vsize)) then
    deallocate ( vsize , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vsize",2,ierr,error_message)
  end if
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
                                   xadj,adjncy,epart)
  !
  !.. Use Statements ..
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  use iso_c_binding, only : c_null_ptr
  !
  !.. Formal Arguments ..
  integer,                  intent(in) :: metis_option_requested
  integer,                  intent(in) :: npart
  integer, dimension(:),    intent(in) :: xadj
  integer, dimension(:),    intent(in) :: adjncy
  integer, dimension(:), intent(inout) :: epart
  !
  !.. Local Scalars ..
  integer :: n,ip,ierr,lcell
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
  lcell = size(epart)
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
  allocate ( vwgt(1:ncon*lcell) , source=1 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"vwgt",1,ierr,error_message)
  allocate ( vsize(1:lcell) , source=1 , stat=ierr , errmsg=error_message )
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
        call metis_partgraphrecursive(lcell,ncon,xadj,adjncy,vwgt,vsize, &
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
        call metis_partgraphkway(lcell,ncon,xadj,adjncy,vwgt,vsize, &
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
        call metis_partgraphkway(lcell,ncon,xadj,adjncy,vwgt,vsize, &
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
        call metis_partgraphkway(lcell,ncon,xadj,adjncy,vwgt,vsize, &
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
subroutine get_global_connectivity(cell_geom,cell_order,bface,xyz_nodes, &
                                   nodes_of_cell,nodes_of_cell_ptr, &
                                   cells_with_node,cells_with_node_ptr, &
                                   cells_surr_cell,cells_surr_cell_ptr, &
                                   cells_with_edge,cells_with_edge_ptr, &
                                   nodes_on_edge,face)
  !
  !.. Use Statements ..
  use connectivity_mod, only : find_cells_surrounding_each_node
  use connectivity_mod, only : find_cells_surrounding_each_cell
  use connectivity_mod, only : find_nodes_on_each_edge
  use connectivity_mod, only : find_cells_surrounding_each_edge
  use connectivity_mod, only : create_face_derived_type
  !
  !.. Formal Arguments ..
  integer,              dimension(:),      intent(in) :: cell_geom
  integer,              dimension(:),      intent(in) :: cell_order
  integer,              dimension(:,:),    intent(in) :: bface
  integer,              dimension(:),      intent(in) :: nodes_of_cell
  integer,              dimension(:),      intent(in) :: nodes_of_cell_ptr
  real(wp),             dimension(:,:),    intent(in) :: xyz_nodes
  integer, allocatable, dimension(:),   intent(inout) :: cells_with_node
  integer, allocatable, dimension(:),   intent(inout) :: cells_with_node_ptr
  integer, allocatable, dimension(:),   intent(inout) :: cells_surr_cell
  integer, allocatable, dimension(:),   intent(inout) :: cells_surr_cell_ptr
  integer, allocatable, dimension(:),   intent(inout) :: cells_with_edge
  integer, allocatable, dimension(:),   intent(inout) :: cells_with_edge_ptr
  integer, allocatable, dimension(:,:), intent(inout) :: nodes_on_edge
  !
  type(face_t), allocatable, dimension(:), intent(inout) :: face
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_global_connectivity"
  !
continue
  !
  ! Find the cells that surround each node
  !
  call find_cells_surrounding_each_node(bface, &
                                        nodes_of_cell,nodes_of_cell_ptr, &
                                        cells_with_node,cells_with_node_ptr, &
                                        include_bnd_faces=true, &
                                        make_bnd_faces_negative=true)
  !
  ! Find the cells surrounding each cell which
  ! defines the dual graph of the nodal grid 
  !
  call find_cells_surrounding_each_cell(bface,cell_geom, &
                                        nodes_of_cell,nodes_of_cell_ptr, &
                                        cells_with_node,cells_with_node_ptr, &
                                        cells_surr_cell,cells_surr_cell_ptr)
  !
  ! Find the nodes defining each edge as well as the edges surrounding each node
  !
  call find_nodes_on_each_edge(cell_geom, &
                               nodes_of_cell,nodes_of_cell_ptr, &
                               cells_with_node,cells_with_node_ptr, &
                               nodes_on_edge)
  !
  ! Find the cells that surround each edge
  !
  call find_cells_surrounding_each_edge(nodes_on_edge, &
                                        nodes_of_cell,nodes_of_cell_ptr, &
                                        cells_with_node,cells_with_node_ptr, &
                                        cells_with_edge,cells_with_edge_ptr)
  !
  ! Find the cells on each face. Once the grid is partitioned,
  ! this will be needed to adjust the boundary faces to include
  ! former interior faces that now separate adjacent partitions
  !
  call create_face_derived_type(bface,cell_geom,cell_order,xyz_nodes, &
                                nodes_of_cell,nodes_of_cell_ptr, &
                                cells_surr_cell,cells_surr_cell_ptr, &
                                face,force_flxpt_offset_to_zero=true)
  !
end subroutine get_global_connectivity
!
!###############################################################################
!
subroutine create_metis_dual_graph(cells_surr_cell,cells_surr_cell_ptr, &
                                   xadj,adjncy)
  !
  ! This subroutine takes the grid cells and nodes and creates the
  ! dual graph arrays required by METIS to partition the global grid
  !
  !.. Formal Arguments ..
  integer,              dimension(:),    intent(in) :: cells_surr_cell
  integer,              dimension(:),    intent(in) :: cells_surr_cell_ptr
  integer, allocatable, dimension(:), intent(inout) :: xadj
  integer, allocatable, dimension(:), intent(inout) :: adjncy
  !
  !.. Local Scalars ..
  integer :: n,nf,n1,n2,ic,ierr,lcell
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_metis_dual_graph"
  !
continue
  !
  lcell = size(cells_surr_cell_ptr)-1
  !
  ! We cant just use cells_surr_cell_ptr and cells_surr_cell as xadj and adjncy
  ! because cells_surr_cell is:
  !
  !   1. negative for faces of cells that correspond to boundary faces if the
  !      optional arguments include_bnd_faces and make_bnd_faces_negative are
  !      both sent into subroutine find_cells_surrounding_each_node with the
  !      value .TRUE.
  !
  !   OR
  !
  !   2. zero for faces of cells that correspond to boundary faces if the
  !      optional argument include_bnd_faces was either absent or present with
  !      the value .FALSE. in the call to the subroutine
  !      find_cells_surrounding_each_node above. In the loop "find_cell" within
  !      subroutine find_cells_surrounding_each_cell, if the
  !      sum(lpoin(nodes_of_cell(j1:j2))) term never equals nfnods then this
  !      loop exits with ieadj still equal to 0.
  !
  ! We cannot send adjncy containing values less than 1 to the METIS routines so
  ! these need to be removed and the locations in xadj also have to be adjusted
  ! accordingly. Additionally, have xadj to start the indexing for adjncy at 1,
  ! i.e., xadj(1)=1 (cells_surr_cell_ptr indexing for cells_surr_cell begins
  ! at 0, i.e. cells_surr_cell_ptr(1)=0) since this is what METIS expects when
  ! using Fortran array ordering.
  !
  ! Allocate the xadj and adjncy arrays
  !
  allocate ( xadj(1:lcell+1) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xadj",1,ierr,error_message)
  xadj(1) = 1
  !
  n = count( cells_surr_cell(:) > 0 )
  allocate ( adjncy(1:n) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"adjncy",1,ierr,error_message)
  !
  ! Loop through the cells to copy only the greater than zero
  ! locations of cells_surr_cell into the adjncy array and
  ! keep count of the number copied for each cell in xadj
  !
  ic = 0
  do n = 1,lcell
    n1 = cells_surr_cell_ptr(n)+1
    n2 = cells_surr_cell_ptr(n+1)
    do nf = n1,n2
      if (cells_surr_cell(nf) > 0) then
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
end subroutine create_metis_dual_graph
!
!###############################################################################
!
subroutine create_grid_weights(metis_option,ncon,n_dual_vertices,n_dual_edges, &
                               vwgt,adjwgt,vsize,wgtflag,cell_geom,cell_order)
  !
  !.. Formal Arguments ..
  integer,                               intent(in) :: metis_option
  integer,                               intent(in) :: ncon
  integer,                               intent(in) :: n_dual_vertices
  integer,                               intent(in) :: n_dual_edges
  integer,              dimension(:),    intent(in) :: cell_geom
  integer,              dimension(:),    intent(in) :: cell_order
  integer,                            intent(inout) :: wgtflag
  integer, allocatable, dimension(:), intent(inout) :: vwgt
  integer, allocatable, dimension(:), intent(inout) :: adjwgt
  integer, allocatable, dimension(:), intent(inout) :: vsize
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_grid_weights"
  !
continue
  !
  ! Initialize the weighting flag
  !
  wgtflag = 0
  !
  ! Initialize vwgt
  ! NOTE: vwgt gives the computational weight of each cell
  !
  if (allocated(vwgt)) then
    deallocate ( vwgt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vwgt",2,ierr,error_message)
  end if
  !
  allocate ( vwgt(1:ncon*n_dual_vertices) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"vwgt",1,ierr,error_message)
  !
 !vwgt = 0
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
    ! Allocate the adjwgt array
    !
    if (allocated(adjwgt)) then
      deallocate ( adjwgt , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"adjwgt",2,ierr,error_message)
    end if
    !
    allocate ( adjwgt(1:n_dual_edges) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"adjwgt",1,ierr,error_message)
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
    ! Determine whether edge weights are being used
    ! and adjust wgtflag accordingly
    !
    if (minval(adjwgt) /= maxval(adjwgt)) wgtflag = wgtflag + 1
    !
    ! Determine whether vertex weights are being used
    ! and adjust wgtflag accordingly
    !
    if (minval(vwgt) /= maxval(vwgt)) wgtflag = wgtflag + 2
    !
  else if ( any(metis_option == (/3,4/)) ) then
    !
    ! Allocate the vsize array
    !
    if (allocated(vsize)) then
      deallocate ( vsize , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"vsize",2,ierr,error_message)
    end if
    !
    allocate ( vsize(1:n_dual_vertices) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vsize",1,ierr,error_message)
    !
    ! If minimizing the communication volume, vsize gives the communication
    ! weight of the cells (i.e., how much data has to be sent to/from this cell)
    !
    ! Set the communication weight of the cells to
    ! the number of solution points for each cell
    !
    vsize = Cell_Solpts( cell_geom , cell_order )
    !
    ! Determine whether communication weights are being used
    ! and adjust wgtflag accordingly
    !
    if (minval(vsize) /= maxval(vsize)) wgtflag = wgtflag + 1
    !
    ! Determine whether computation weights are being used
    ! and adjust wgtflag accordingly
    !
    if (minval(vwgt) /= maxval(vwgt)) wgtflag = wgtflag + 2
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
subroutine output_value(char_input,val_input,io_unit,extra_line,surpress_cpu)
  !
  character(len=*),  intent(in) :: char_input
  integer,           intent(in) :: val_input
  integer, optional, intent(in) :: io_unit
  logical, optional, intent(in) :: extra_line
  logical, optional, intent(in) :: surpress_cpu
  !
  integer :: n,io
  logical :: add_extra_line
  logical :: include_cpu
  character(len=10) :: cval
  character(len=50) :: frmt
  !
continue
  !
  io = iout
  if (present(io_unit)) then
    io = io_unit
  end if
  !
  add_extra_line = true
  if (present(extra_line)) then
    add_extra_line = extra_line
  end if
  !
  include_cpu = true
  if (present(surpress_cpu)) then
    include_cpu = surpress_cpu
  end if
  !
  if (include_cpu) then
    !
    write (cval,'(i0)') ncpu
    write (frmt,1) len_trim(adjustl(cval))
    !
    flush (io)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    !
    do n = 0,ncpu-1
      if (mypnum == n) then
        write (io,frmt) mypnum,trim(adjustl(char_input)),val_input
        flush (io)
      end if
      call mpi_barrier(MPI_COMM_WORLD,mpierr)
    end do
    !
    flush (io)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    !
    if (mypnum == 0) write (io,*)
    !
    flush (io)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    !
  else
    !
    write (io,2) trim(adjustl(char_input)),val_input
    if (add_extra_line) then
      write (io,*)
    end if
    !
  end if
  !
  1 format ('(" CPU ",i',i0,',": ",a," = ",i0)')
  2 format (a," = ",i0)
  !
end subroutine output_value
!
!###############################################################################
!
subroutine create_node_and_edge_arrays(nodes_of_cell,nodes_of_cell_ptr, &
                                       cells_with_node,cells_with_node_ptr, &
                                       cells_with_edge,cells_with_edge_ptr, &
                                      !nodes_on_edge, &
                                       cell_order,epart,bface,node,edge)
  !
  !.. Formal Arguments ..
  integer,                   dimension(:),    intent(in) :: nodes_of_cell
  integer,                   dimension(:),    intent(in) :: nodes_of_cell_ptr
  integer,                   dimension(:),    intent(in) :: cells_with_node
  integer,                   dimension(:),    intent(in) :: cells_with_node_ptr
  integer,                   dimension(:),    intent(in) :: cells_with_edge
  integer,                   dimension(:),    intent(in) :: cells_with_edge_ptr
 !integer,                   dimension(:,:),  intent(in) :: nodes_on_edge
  integer,                   dimension(:),    intent(in) :: epart
  integer,                   dimension(:),    intent(in) :: cell_order
  integer,                   dimension(:,:),  intent(in) :: bface
  type(node_t), allocatable, dimension(:), intent(inout) :: node
  type(edge_t), allocatable, dimension(:), intent(inout) :: edge
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,nf,np,i,ierr
  integer :: lnode,ledge,lcell
  integer :: lcpu,ucpu,num_cpus
  integer :: host_cell,ghost_cell
  character(len=100) :: array_name
  !
  !.. Local Arrays ..
  logical, dimension(1:ncpu) :: cpu_mask
  integer, dimension(1:ncpu) :: cpu_cells
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:) :: bc_list
  integer, allocatable, dimension(:) :: cell_list
  integer, allocatable, dimension(:) :: part_list
  !
  !.. Local Paramters ..
  character(len=*), parameter :: pname = "create_node_and_edge_arrays"
  !
continue
  !
  lnode = size(cells_with_node_ptr)-1
  ledge = size(cells_with_edge_ptr)-1
  lcell = size(epart)
  !
  ! Find the maximum number of cells with a single node
  !
  n1 = maxval( cells_with_node_ptr(2:lnode+1) - &
               cells_with_node_ptr(1:lnode  ) )
  !
  ! Find the maximum number of cells with a single edge
  !
  n2 = maxval( cells_with_edge_ptr(2:ledge+1) - &
               cells_with_edge_ptr(1:ledge  ) )
  !
  ! Allocate the array to store the partitions numbers for all cells with a
  ! given node/edge, and the array to store the BCs those boundary faces
  ! with a given node/edge.
  !
  n = max(n1,n2)
  !
  allocate ( cell_list(1:n) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_list",1,ierr,error_message)
  !
  allocate ( part_list(1:n) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"part_list",1,ierr,error_message)
  !
  allocate ( bc_list(1:n) , source=not_a_bc , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bc_list",1,ierr,error_message)
  !
  ! Create the node array
  !
  allocate ( node(1:lnode) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"node",1,ierr,error_message)
  !
  ! Loop through the grid nodes and find any nodes that exist on multiple
  ! partitions. Save these partition numbers to the node array
  !
  do n = 1,lnode
    !
    ! Reset the local arrays for this node
    !
    cpu_mask(:) = fals
    cpu_cells(:) = 0
    part_list(:) = 0
    cell_list(:) = 0
    bc_list(:) = not_a_bc
    !
    n1 = cells_with_node_ptr(n)+1 ! index of first cell containing this node
    n2 = cells_with_node_ptr(n+1) ! index of last  cell containing this node
    np = n2-n1+1                  ! number of cells containing this node
    !
    ! We first need to check if this node lies
    ! on a non-communication boundary face
    !
    if (any( cells_with_node(n1:n2) < 0 )) then
      !
      ! One or more of the cells containing this node are actually a
      ! non-communication boundary face so we need to treat this node
      ! a little differently
      !
      ! Set the logical flag that identifies
      ! this node is part of a boundary face
      !
      node(n)%is_on_boundary = true
      !
      np = 0
      !
      do i = n1,n2
        !
        if (cells_with_node(i) < 0) then
          !
          np = np + 1
          !
          nf = abs(cells_with_node(i)) - lcell
          !
          host_cell  = bface(2,nf)
          ghost_cell = bface(10,nf)
          !
          if (abs(cells_with_node(i)) /= ghost_cell) then
            write (error_message,2) n,n,nf,cells_with_node(i), &
                                    ghost_cell,host_cell
            call stop_recon(stop_mpi,pname,error_message)
          end if
          !
          ! Store the host cell for this boundary face
          !
          cell_list(np) = host_cell
          !
          ! Store the B.C. for this boundary face
          !
          bc_list(np) = bface(1,nf)
          !
        end if
        !
      end do
      !
      ! Save the partitions numbers of the cells with this edge to a local array
      !
      part_list(1:np) = epart( cell_list(1:np) )
      !
      ! Find the boundary faces with this node that have the highest
      ! BC priority, and set the BC to bc_unknown for all the other
      ! boundary faces with a BC of lesser priority.
      !
      call find_priority_bc( bc_list(1:np) , cell_list(1:np) )
      !
    else
      !
      ! This node is not part of a boundary face
      !
      ! Set bc_list to not_a_bc to indicate that we want
      ! to add each of these cells to cpu_cells
      !
      bc_list(1:np) = not_a_bc
      !
      ! Save the cells with this node to a local array
      !
      cell_list(1:np) = cells_with_node(n1:n2)
      !
      ! Save the partitions numbers of the cells with this edge to a local array
      !
      part_list(1:np) = epart( cell_list(1:np) )
      !
    end if
    !
    ! Get the lower and upper bounds of the partitions with this node
    !
    lcpu = minval( part_list(1:np) )
    ucpu = maxval( part_list(1:np) )
    !
    ! Skip this node if the partition number is
    ! the same for all cells with this node
    !
    if (lcpu == ucpu) cycle
    !
    ! Sum the number of cells containing this node on each partition
    !
    do i = 1,np
      !
      ! Mark this partition as having a cell containing this node
      !
      cpu_mask( part_list(i) ) = true
      !
      ! Only add another cell to this partition if this is not
      ! a boundary face (i.e., bc_list(i) == not_a_bc) or this
      ! is a boundary face with the BC of highest priority for
      ! this node (i.e., bc_list(i) /= bc_unknown).
      !
      if (bc_list(i) /= bc_unknown) then
        cpu_cells( part_list(i) ) = cpu_cells( part_list(i) ) + 1
      end if
      !
    end do
    !
    ! Count the number of cpus that contain this node
    !
    num_cpus = count( cpu_mask(lcpu:ucpu) )
    !
    ! Allocate the cpu component of the node array for this node
    !
    allocate ( node(n)%cpu(1:num_cpus) , stat=ierr , errmsg=error_message )
    write (array_name,1) "node",n
    call alloc_error(pname,array_name,1,ierr,error_message,skip_alloc_pause)
    !
    node(n)%cpu(1:num_cpus)%partition = pack( intseq(lcpu,ucpu) , &
                                              cpu_mask(lcpu:ucpu) )
    !
    ! Save the number of cells that contain this node on each partition
    ! NOTE: This may be zero for a partition if:
    !         1) this node is on a boundary face, and
    !         2) on this partition, all the boundary faces with this node
    !            have BCs that are of lesser priority than a higher
    !            priority BC for a boundary face on another partition.
    !       In such a case, we want this to be zero because that means the
    !       averaging correction for this node (computed later in subroutine
    !       create_node_and_edge_datatypes) will also be zero for this
    !       partition, and therefore the contribution to the solution for
    !       this node from this partition will be canceled out.
    !
    do i = 1,size(node(n)%cpu)
      node(n)%cpu(i)%num_cells = cpu_cells( node(n)%cpu(i)%partition )
    end do
    !
    ! Assign the partition that contains the cell with the lowest index
    ! value with this node as the controlling partition for this node
    !
    i = minval( cell_list , mask=(cell_list>0 .and. cell_list<lcell) )
    !
    node(n)%controlling_partition = epart(i)
    !
  end do
  !
  ! Create the edge array
  !
  allocate ( edge(1:ledge) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"edge",1,ierr,error_message)
  !
  ! Loop through the grid edges and find any edges that exist on multiple
  ! partitions. Save these partition numbers to the edge array
  !
  do n = 1,ledge
    !
    ! Reset the local arrays for this edge
    !
    cpu_mask(:) = fals
    cpu_cells(:) = 0
    part_list(:) = 0
    cell_list(:) = 0
    bc_list(:) = not_a_bc
    !
    n1 = cells_with_edge_ptr(n)+1 ! index of first cell containing this edge
    n2 = cells_with_edge_ptr(n+1) ! index of last  cell containing this edge
    np = n2-n1+1                  ! number of cells containing this edge
    !
    ! We first need to check if this edge lies
    ! on a non-communication boundary face
    !
    if (any( cells_with_edge(n1:n2) < 0 )) then
      !
      ! One or more of the cells containing this edge are actually a
      ! non-communication boundary face so we need to treat this edge
      ! a little differently
      !
      ! Set the logical flag that identifies
      ! this edge is part of a boundary face
      !
      edge(n)%is_on_boundary = true
      !
      edge(n)%order = 0
      np = 0
      !
      do i = n1,n2
        !
        if (cells_with_edge(i) < 0) then
          !
          np = np + 1
          !
          nf = abs(cells_with_edge(i)) - lcell
          !
          host_cell  = bface(2,nf)
          ghost_cell = bface(10,nf)
          !
          if (abs(cells_with_edge(i)) /= ghost_cell) then
            write (error_message,3) n,n,nf,cells_with_edge(i), &
                                    ghost_cell,host_cell
            call stop_recon(stop_mpi,pname,error_message)
          end if
          !
          ! Store the host cell for this boundary face
          !
          cell_list(np) = host_cell
          !
          ! Store the B.C. for this boundary face
          !
          bc_list(np) = bface(1,nf)
          !
        end if
        !
      end do
      !
      ! Temporarily store the host cells in part_list
      !
      part_list(1:np) = cell_list(1:np)
      !
      ! Find the boundary faces with this edge that have the highest
      ! BC priority, and set the BC to bc_unknown for all the other
      ! boundary faces with a BC of lesser priority.
      !
      call find_priority_bc( bc_list(1:np) , cell_list(1:np) )
      !
      ! Now that we have removed all the boundary faces with BC of lower
      ! priority, find the maximum order of all the host cells for the
      ! remaining boundary faces. Also replace the index for the host cells
      ! that are temporarily stored in part_list with the partition numbers
      ! for each of these host cells.
      !
      do i = 1,np
        host_cell = part_list(i)
        if (bc_list(i) /= bc_unknown) then
          edge(n)%order = max( edge(n)%order , cell_order(host_cell) )
        end if
        part_list(i) = epart( host_cell )
      end do
      !
    else
      !
      ! This edge is not part of a boundary face
      !
      ! Set bc_list to not_a_bc to indicate that we want
      ! to add each of these cells to cpu_cells
      !
      bc_list(1:np) = not_a_bc
      !
      ! Save the cells with this edge to a local array
      !
      cell_list(1:np) = cells_with_edge(n1:n2)
      !
      ! Save the partitions numbers of the cells with this edge to a local array
      !
      part_list(1:np) = epart( cell_list(1:np) )
      !
      ! Get the maximum order of the cells with this edge
      ! and save this to the edge array
      !
      edge(n)%order = maxval( cell_order(cell_list(1:np)) )
      !
    end if
    !
    ! Get the lower and upper bounds of the partitions with this edge
    !
    lcpu = minval( part_list(1:np) )
    ucpu = maxval( part_list(1:np) )
    !
    ! Skip this edge if the partition number is
    ! the same for all cells with this edge
    !
    if (lcpu == ucpu) cycle
    !
    ! Sum the number of cells containing this edge on each partition
    !
    do i = 1,np
      !
      ! Mark this partition as having a cell containing this edge
      !
      cpu_mask( part_list(i) ) = true
      !
      ! Only add another cell to this partition if this is not
      ! a boundary face (i.e., bc_list(i) == not_a_bc) or this
      ! is a boundary face with the BC of highest priority for
      ! this edge (i.e., bc_list(i) /= bc_unknown).
      !
      if (bc_list(i) /= bc_unknown) then
        cpu_cells( part_list(i) ) = cpu_cells( part_list(i) ) + 1
      end if
      !
    end do
    !
    ! Count the number of cpus that contain this edge
    !
    num_cpus = count( cpu_mask(lcpu:ucpu) )
    !
    ! Allocate the cpu component of the edge array for this edge
    !
    allocate ( edge(n)%cpu(1:num_cpus) , stat=ierr , errmsg=error_message )
    write (array_name,1) "edge",n
    call alloc_error(pname,array_name,1,ierr,error_message,skip_alloc_pause)
    !
    edge(n)%cpu(1:num_cpus)%partition = pack( intseq(lcpu,ucpu) , &
                                              cpu_mask(lcpu:ucpu) )
    !
    ! Save the number of cells that contain this edge on each partition
    ! NOTE: This may be zero for a partition if:
    !         1) this edge is on a boundary face, and
    !         2) on this partition, all the boundary faces with this edge
    !            have BCs that are of lesser priority than a higher
    !            priority BC for a boundary face on another partition.
    !       In such a case, we want this to be zero because that means the
    !       averaging correction for this edge (computed later in subroutine
    !       create_node_and_edge_datatypes) will also be zero for this
    !       partition, and therefore the contribution to the solution for
    !       this edge from this partition will be canceled out.
    !
    do i = 1,size(edge(n)%cpu)
      edge(n)%cpu(i)%num_cells = cpu_cells( edge(n)%cpu(i)%partition )
    end do
    !
    ! Assign the partition that contains the cell with the lowest index
    ! value with this edge as the controlling partition for this edge
    !
    i = minval( cell_list , mask=(cell_list>0 .and. cell_list<lcell) )
    !
    edge(n)%controlling_partition = epart(i)
    !
    ! Finally, we need to create the ordering of the edge points
    !
    if (use_edgpts) then
      !
      np = Cell_Solpts( Geom_Edge , edge(n)%order )
      !
      do i = 1,size(edge(n)%cpu)
        !
        allocate ( edge(n)%cpu(i)%edgpts(1:np) , source=0 , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) "edge",n,i
        call alloc_error(pname,array_name,1,ierr,error_message,skip_alloc_pause)
        !
        ! NOT SURE IF THIS IS CORRECT. WE MIGHT NEED TO CHECK THE NODES ON
        ! THE EDGE AND THE ORDERING AND USE intseq(np,1,-1) IF THE NODES ARE
        ! REVERSED
        !
        edge(n)%cpu(i)%edgpts(1:np) = intseq(1,np)
        !
      end do
      !
    end if
    !
  end do
  !
  ! Deallocate cell_list, part_list, and bc_list before we leave
  !
  deallocate ( cell_list , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_list",2,ierr,error_message)
  !
  deallocate ( part_list , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"part_list",2,ierr,error_message)
  !
  deallocate ( bc_list , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bc_list",2,ierr,error_message)
  !
  ! Format Statements
  !
  1 format(a,"(",i0,")%cpu",:,"(",i0,")%edgpts")
  2 format("ERROR! A cell containing the node ",i0," was found to ", &
           "be a boundary",/,8x,"face, but the ghost cell value in ", &
           "cells_with_node was not the",/,8x,"negative of the ghost ", &
           "cell value in the face array.",/, &
           11x,"node index = ",i0,/, &
           11x,"face index = ",i0,/, &
           11x,"cells_with_node(i)  = ",i0,/, &
           11x,"face(nf)%right%cell = ",i0,/, &
           11x,"face(nf)%left%cell  = ",i0,/)
  3 format("ERROR! A cell containing the edge ",i0," was found to ", &
           "be a boundary",/,8x,"face, but the ghost cell value in ", &
           "cells_with_edge was not the",/,8x,"negative of the ghost ", &
           "cell value in the face array.",/, &
           11x,"edge index = ",i0,/, &
           11x,"face index = ",i0,/, &
           11x,"cells_with_edge(i)  = ",i0,/, &
           11x,"face(nf)%right%cell = ",i0,/, &
           11x,"face(nf)%left%cell  = ",i0,/)
  !
end subroutine create_node_and_edge_arrays
!
!###############################################################################
!
subroutine find_priority_bc(bc_list,cell_list)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(inout) :: bc_list
  !
  !.. Optional Arguments ..
  integer, optional, dimension(:), intent(inout) :: cell_list
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Arrays ..
 !integer, dimension(1:size(bc_list)) :: array_of_zeros
  !
continue
  !
 !array_of_zeros(:) = 0
  !
  ! Loop through the priority heirarchy of the different BC
  ! until we find the BC of highest priority within bc_list
  !
  do n = 1,size(bc_heirarchy)
    !
    if (any(bc_list == bc_heirarchy(n))) then
      !
      ! We found the BC of highest priority within bc_list.
      ! Now set all other BCs within bc_list to the parameter
      ! bc_unknown so that only the current BC remains.
      !
      where (bc_list /= bc_heirarchy(n)) bc_list = bc_unknown
      !
      ! Exit now that we have collected the boundary faces
      ! that have the BC of highest priority
      !
      exit
      !
    end if
    !
  end do
  !
  ! Now pack the host cells in cell_list that correspond to the
  ! non-bc_unknown boundary faces in bc_list up to the front of cell_list
  ! and pad the remaining elements of cell_list with zeros.
  !
 !if (present(cell_list)) then
 !  cell_list = pack( cell_list , bc_list /= bc_unknown , array_of_zeros )
 !end if
  !
end subroutine find_priority_bc
!
!###############################################################################
!
subroutine find_partition_boundaries(cell_geom,epart,bface,fpart,face)
  !
  !.. Formal Arguments ..
  integer,                   dimension(:),      intent(in) :: cell_geom
  integer,                   dimension(:),      intent(in) :: epart
  integer,      allocatable, dimension(:,:), intent(inout) :: bface
  integer,      allocatable, dimension(:),   intent(inout) :: fpart
  type(face_t), allocatable, dimension(:),   intent(inout) :: face
  !
  !.. Local Scalars ..
  integer :: new_lfbnd,new_lface
  integer :: lfbnd,lface,lcell
  integer :: n,nf,n1,n2,np,ierr
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
  lfbnd = size(bface,dim=2)
  lface = size(face)
  lcell = size(epart)
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  allocate ( msk(1:lface) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",1,ierr,error_message)
  !
  ! Find the interior faces that divide two partitions
  !
  do nf = lfbnd+1,lface
    msk(nf) = epart( face(nf)%left%cell ) /= epart( face(nf)%right%cell )
  end do
  !
  ! Recompute the number of boundary faces and the total number of faces.
  !   Each internal face dividing unstructured blocks will be turned
  !   into two boundary faces, one for each side of the face meaning:
  !   1. the number of internal faces will shrink by count(msk)
  !   2. the number of boundary faces will increase by 2*count(msk)
  !   3. the total number of faces will increase by count(msk)
  !
  new_lfbnd = lfbnd + 2*count(msk)
  new_lface = lface + 1*count(msk)
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Reallocate bface to make room for new partition boundary faces
  !
  call reallocate(bface,new_size2=new_lfbnd)
  !
  ! Allocate the fpart array that identifies which partition/processor is
  ! is responsible for the face/flux point indices for each boundary face
  ! on a partition boundary thats needed for creating the CGNS connectivity
  !
  allocate ( fpart(1:new_lfbnd) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fpart",1,ierr,error_message)
  !
  ! Replace bface(10,nf) with the location within the face array to extract data
  ! and set the responsible partition to the partition for the host cell
  !
  do nf = 1,lfbnd
    bface(10,nf) = nf
    fpart(nf) = epart( bface(2,nf) )
  end do
  !
  ! Collect and add the new boundary faces to bface
  !
  nf = lfbnd + 1
  do n = lfbnd+1,lface
    !
    if (.not. msk(n)) cycle
    !
    left_cell = face(n)%left%cell
    left_side = face(n)%left%cell_face
    left_geom = cell_geom(left_cell)
    !
    right_cell = face(n)%right%cell
    right_side = face(n)%right%cell_face
    right_geom = cell_geom(right_cell)
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
    bface( 9,nf) = right_side  ! face of adj partition host cell on boundary
    bface(10,nf) = n           ! location within face array to extract data
    !
    ! Set the responsible partition for this new boundary face to the
    ! partition for what was originally the left cell on this interior face
    !
    fpart(nf) = epart(left_cell)
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
    bface( 9,nf) = left_side   ! face of adj partition host cell on boundary
    bface(10,nf) = -n          ! location within face array to extract data
    !                          ! NOTE: This is negative to mark that this host
    !                          !       cell used to be a right cell on the face
    !
    ! Set the responsible partition for this new boundary face to the
    ! partition for what was originally the left cell on this interior face
    !
    fpart(nf) = epart(left_cell)
    !
    nf = nf + 1
    !
  end do
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Allocate a temporary array to resize and reorder face
  !
  allocate ( temp_face(1:new_lface) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_face",1,ierr,error_message)
  !
  ! Redo the boundary face section of face
  ! because of the new additions to bface
  !
  do nf = 1,new_lfbnd
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
    bface(10,nf) = lcell + nf
    !
  end do
  !
  ! Copy the interior face information to the temporary array
  ! but skip those interior faces that are now boundary faces
  !
  nf = new_lfbnd + 1
  do n = lfbnd+1,lface
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
  lfbnd = new_lfbnd
  lface = new_lface
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
  do nf = 1,lfbnd
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
  if ( any( epart(face(lfbnd+1:lface)%left%cell) /= &
            epart(face(lfbnd+1:lface)%right%cell) ) ) then
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
subroutine output_node_arrays(node,node_map,bnd_nodes)
  !
  type(node_t),     optional, dimension(:), intent(in) :: node
  type(map_t),      optional, dimension(:), intent(in) :: node_map
  type(bnd_data_t), optional, dimension(:), intent(in) :: bnd_nodes
  !
  integer :: i,n,io
  !
continue
  !
  if (present(node)) then
    io = mypnum+140
    do n = 1,size(node)
      write (io,10) n
      write (io,11) node(n)%is_on_boundary
      write (io,12) node(n)%controlling_partition
      if (allocated(node(n)%cpu)) then
        do i = 1,size(node(n)%cpu)
          write (io,13) i
          write (io,14) node(n)%cpu(i)%partition
          write (io,15) node(n)%cpu(i)%loc_node
          write (io,16) node(n)%cpu(i)%num_cells
        end do
      end if
    end do
  end if
  !
  10 format (" Contents of node(",i0,")")
  11 format (5x,"%is_on_boundary        = ",l1)
  12 format (5x,"%controlling_partition = ",i0)
  13 format (5x,"%cpu(",i0,")")
  14 format (10x,"%partition = ",i0)
  15 format (10x,"%loc_node  = ",i0)
  16 format (10x,"%num_cells = ",i0)
  !
  if (present(node_map)) then
    io = mypnum+150
    do n = 1,size(node_map)
      write (io,20) n
      do i = 1,size(node_map(n)%loc_to_glb)
        write (io,21) i,node_map(n)%loc_to_glb(i)
      end do
      write (io,22)
    end do
  end if
  !
  20 format (" Map for CPU ",i0)
  21 format ("   loc = ",i5," => glb = ",i5)
  22 format (30("#"))
  !
  if (present(bnd_nodes)) then
    io = mypnum+160
    do n = 1,size(bnd_nodes)
      write (io,30) n
      write (io,31) bnd_nodes(n)%idx
      write (io,32) bnd_nodes(n)%ave_correction
      write (io,33) bnd_nodes(n)%my_responsibility
    end do
  end if
  !
  30 format (" Contents of bnd_nodes(",i0,")")
  31 format (5x,"%idx               = ",i0)
  32 format (5x,"%ave_correction    = ",f9.6)
  33 format (5x,"%my_responsibility = ",l1)
  !
end subroutine output_node_arrays
!
!###############################################################################
!
subroutine create_grid_mappings(epart,face,node,edge, &
                                nodes_of_cell,nodes_of_cell_ptr, &
                                cells_with_edge,cells_with_edge_ptr, &
                                node_map,edge_map,face_map,cell_map)
  !
  !.. Formal Arguments ..
  integer,                  dimension(:),    intent(in) :: epart
  integer,                  dimension(:),    intent(in) :: nodes_of_cell
  integer,                  dimension(:),    intent(in) :: nodes_of_cell_ptr
  integer,                  dimension(:),    intent(in) :: cells_with_edge
  integer,                  dimension(:),    intent(in) :: cells_with_edge_ptr
  type(face_t),             dimension(:),    intent(in) :: face
  type(node_t),             dimension(:), intent(inout) :: node
  type(edge_t),             dimension(:), intent(inout) :: edge
  type(map_t), allocatable, dimension(:), intent(inout) :: node_map
  type(map_t), allocatable, dimension(:), intent(inout) :: edge_map
  type(map_t), allocatable, dimension(:), intent(inout) :: face_map
  type(map_t), allocatable, dimension(:), intent(inout) :: cell_map
  !
  !.. Local Scalars ..
  integer :: n,nc,ne,n1,n2,ierr
  integer :: npart,this_part,ncount
  integer :: loc_edge,glb_edge
  integer :: loc_node,glb_node
  integer :: lcell,lnode,lface,ledge
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
  npart = maxval(epart) - minval(epart) + 1
  lcell = size(epart)
  lnode = size(node)
  ledge = size(edge)
  lface = size(face)
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Allocate the msk array
  !
  ncount = max( lcell , lnode , lface , ledge )
  !
  allocate ( msk(1:ncount) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",1,ierr,error_message)
  !
  ! Create the cell mappings
  !
  allocate ( cell_map(1:npart) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_map",1,ierr,error_message)
  !
  ! Loop through the partitions to find the
  ! local number of cells for each partition
  !
  do this_part = 1,npart
    !
    ! Mark the cells on this partition
    !
    msk(1:lcell) = ( epart(1:lcell) == this_part )
    !
    ! Count the cells on this partition
    !
    ncount = count( msk(1:lcell) )
    !
    ! Allocate the loc_to_glb component of cell_map for this partition
    !
    allocate ( cell_map(this_part)%loc_to_glb(1:ncount) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "cell_map",this_part
    call alloc_error(pname,array_name,1,ierr,error_message)
    !
    ! Add the cells of the current partition to cell_map
    !
    cell_map(this_part)%loc_to_glb = pack( intseq(1,lcell) , msk(1:lcell) )
   !cell_map(this_part)%loc_to_glb = find_subset(epart,this_part)
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
  do this_part = 1,npart
    !
    ! Initialize/Reset msk
    !
    msk(:) = fals
    !
    ! Mark the nodes that are within this partition
    !
    do nc = 1,lcell
      if (epart(nc) /= this_part) cycle
      n1 = nodes_of_cell_ptr(nc)+1
      n2 = nodes_of_cell_ptr(nc+1)
      msk(nodes_of_cell(n1:n2)) = true
    end do
    !
    ! Count the nodes on this partition
    !
    ncount = count( msk(1:lnode) )
    !
    ! Allocate the loc_to_glb component of node_map for this partition
    !
    allocate ( node_map(this_part)%loc_to_glb(1:ncount) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "node_map",this_part
    call alloc_error(pname,array_name,1,ierr,error_message)
    !
    ! Add the nodes of the current partition to node_map
    !
    node_map(this_part)%loc_to_glb = pack( intseq(1,lnode) , msk(1:lnode) )
    !
    ! Save the local node index for each node of this partition that
    ! also exists on a different partition
    ! NOTE: The cpu component of the node array should only be allocated
    !       if the global node exists on multiple partitions
    !
    local_nodes: do loc_node = 1,size(node_map(this_part)%loc_to_glb)
      !
      glb_node = node_map(this_part)%loc_to_glb(loc_node)
      !
      if (.not. allocated(node(glb_node)%cpu)) cycle local_nodes
      !
      do n = 1,size(node(glb_node)%cpu)
        if (this_part == node(glb_node)%cpu(n)%partition) then
          node(glb_node)%cpu(n)%loc_node = loc_node
          cycle local_nodes
        end if
      end do
      !
    end do local_nodes
    !
  end do
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
  do this_part = 1,npart
    !
    ! Mark the faces that are within this partition
    !
    msk(1:lface) = ( epart(face(1:lface)%left%cell) == this_part )
    !
    ! Count the faces that are within this partition
    !
    ncount = count( msk(1:lface) )
    !
    ! Allocate the loc_to_glb component of face_map for this partition
    !
    allocate ( face_map(this_part)%loc_to_glb(1:ncount) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "face_map",this_part
    call alloc_error(pname,array_name,1,ierr,error_message)
    !
    ! Add the faces of the current partition to face_map
    !
    face_map(this_part)%loc_to_glb = pack( intseq(1,lface) , msk(1:lface) )
   !face_map(this_part)%loc_to_glb = &
   !                  find_subset(epart(face(1:lface)%left%cell),this_part)
    !
  end do
  !
  ! #######
  ! ## 4 ##
  ! #######
  !
  ! Create the edge mappings
  !
  allocate ( edge_map(1:npart) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"edge_map",1,ierr,error_message)
  !
  ! Loop through the partitions to find the
  ! local number of edges for each partition
  !
  do this_part = 1,npart
    !
    ! Initialize/Reset msk
    !
    msk(:) = fals
    !
    ! Mark the edges that are within this partition
    !
    do ne = 1,ledge
      n1 = cells_with_edge_ptr(ne)+1
      n2 = cells_with_edge_ptr(ne+1)
      msk(ne) = any( epart(cells_with_edge(n1:n2)) == this_part )
    end do
    !
    ! Count the edges on this partition
    !
    ncount = count( msk(1:ledge) )
    !
    ! Allocate the loc_to_glb component of edge_map for this partition
    !
    allocate ( edge_map(this_part)%loc_to_glb(1:ncount) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "edge_map",this_part
    call alloc_error(pname,array_name,1,ierr,error_message)
    !
    ! Add the edges of the current partition to edge_map
    !
    edge_map(this_part)%loc_to_glb = pack( intseq(1,ledge) , msk(1:ledge) )
    !
    ! Save the local edge index for each edge of this partition that
    ! also exists on a different partition
    ! NOTE: The cpu component of the edge array should only be allocated
    !       if the global edge exists on multiple partitions
    ! 
    local_edges: do loc_edge = 1,size(edge_map(this_part)%loc_to_glb)
      !
      glb_edge = edge_map(this_part)%loc_to_glb(loc_edge)
      !
      if (.not. allocated(edge(glb_edge)%cpu)) cycle local_edges
      ! 
      do n = 1,size(edge(glb_edge)%cpu)
        if (this_part == edge(glb_edge)%cpu(n)%partition) then
          edge(glb_edge)%cpu(n)%loc_edge = loc_edge
          cycle local_edges
        end if
      end do
      !
    end do local_edges
    !
  end do
  !
  ! Deallocate msk before we leave
  !
  deallocate ( msk , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",2,ierr,error_message)
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
                         face,node_map,face_map,cell_map)
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  use geovar,    only : bc_conflict
  !
  !.. Formal Arguments ..
  integer,                dimension(:),      intent(in) :: epart
  type(face_t),           dimension(:),      intent(in) :: face
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
  integer :: gnode,gfbnd,gface,gcell
  integer :: lnode,lfbnd,lface,lcell
  integer :: n_totcel,nfnods,nbfai
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
  ! Before we deallocate any of the input arrays, get the size of the global
  ! grid from these arrays
  !
  gnode = size(xyz_nodes,dim=2)
  gfbnd = size(bface,dim=2)
  gface = size(face)
  gcell = size(epart)
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
  lcell = size( cell_map(my_part)%loc_to_glb )
  lnode = size( node_map(my_part)%loc_to_glb )
  lface = size( face_map(my_part)%loc_to_glb )
  lfbnd = count( face_map(my_part)%loc_to_glb <= gfbnd )
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Localize the boundary face array, but not the information it contains.
  !
  nbfai = size(global_bface,dim=1)
  !
  allocate ( bface(1:nbfai,1:lfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,ierr,error_message)
  !
  ! First extract the information that is local to the current
  ! processor from the global_bface array to the bface array.
  ! Also create a mapping array to identify the global ghost cells
  ! needed for the local grid.
  !
  do loc_face = 1,lfbnd
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
    bface(10,loc_face) = lcell + loc_face
    !
  end do
  !
  n_totcel = lcell + lfbnd
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
  do loc_cell = 1,lcell
    glb_cell = cell_map(my_part)%loc_to_glb(loc_cell)
    cell_geom(loc_cell) = global_geom(glb_cell)
    cell_order(loc_cell) = global_order(glb_cell)
  end do
  !
  ! Get the geometry and order for the ghost cells of this partition
  !
  do loc_face = 1,lfbnd
    glb_face = face_map(my_part)%loc_to_glb(loc_face)
    loc_cell = lcell + loc_face
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
  do loc_cell = 1,lcell
    glb_cell = cell_map(my_part)%loc_to_glb(loc_cell)
    nodes_of_cell_ptr(loc_cell+1) = global_nofc_ptr(glb_cell+1) - &
                                    global_nofc_ptr(glb_cell)
  end do
  !
  ! Next set the node connectivity pointer to the number of nodes that
  ! define each of the ghost cells by using the face(:)%nodes array.
  !
  do loc_face = 1,lfbnd
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
  do loc_cell = 1,lcell
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
  do loc_face = 1,lfbnd
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
 !  write (mypnum+130,'(a,i0)') "gcell = ",gcell
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
 !write (mypnum+120,'(a,i0)') "lcell = ",lcell
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
  allocate ( cell_offset(1:gcell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_offset",1,ierr,error_message)
  !
  do glb_cell = 1,gcell-1
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
  allocate ( xyz_nodes(1:nr,1:lnode) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_nodes",1,ierr,error_message)
  !
  ! Localize the interior nodes
  !
  do loc_node = 1,lnode
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
  do nf = 1,lfbnd
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
 !                       nodes_of_cell_ptr(1:lcell+1), &
 !                       xyz_nodes(1:nr,1:lnode))
  !
end subroutine localize_grid
!
!###############################################################################
!
subroutine create_face_nonblocking_datatypes(bface,fpart,face,face_map)
  !
  !.. Formal Arguments ..
  integer,      dimension(:,:), intent(in) :: bface
  integer,      dimension(:),   intent(in) :: fpart
  type(face_t), dimension(:),   intent(in) :: face
  type(map_t),  dimension(:),   intent(in) :: face_map
  !
  !.. Local Scalars ..
  logical :: flag
  integer :: n,nfp,lfbnd,i,ierr
  integer :: rmin,rmax,scnt,rcnt,p1,p2
  integer :: loc_face,glb_face,my_part
  integer :: adj_face,adj_part,adj_cpu
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
  lfbnd = size(bface,dim=2)
  !
  ! Count the number of faces located on partition boundaries
  !
  n = count(bface(1,:) == bc_cpu_bnd)
  !
  ! Allocate the bnd_faces array
  !
  allocate ( bnd_faces(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bnd_faces",1,ierr,error_message)
  !
  ! Now fill in the components within the bnd_faces array
  !
  n = 0
  fill_bnd_faces: do loc_face = 1,lfbnd
    !
    if (bface(1,loc_face) /= bc_cpu_bnd) cycle fill_bnd_faces
    !
    n = n + 1
    !
    bnd_faces(n)%idx = loc_face
    !
    glb_face = face_map(my_part)%loc_to_glb(loc_face)
    !
    if (my_part == fpart(glb_face)) then
      bnd_faces(n)%my_responsibility = true
    end if
    !
  end do fill_bnd_faces
  !
  ! Create an offset array for the number of flux
  ! points on each communication boundary faces
  !
  allocate ( flxpt_offset(1:lfbnd) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"flxpt_offset",1,ierr,error_message)
  !
  do loc_face = 1,lfbnd
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
  do loc_face = lfbnd,1,-1
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
  do loc_face = 1,lfbnd
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
  allocate ( send_msk(1:lfbnd) , source=fals , &
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
    do loc_face = 1,lfbnd
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
    do loc_face = 1,lfbnd
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
    do loc_face = 1,lfbnd
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
  do loc_face = 1,lfbnd
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
subroutine create_node_and_edge_datatypes(node,node_map,edge,edge_map)
  !
  !.. Formal Arguments ..
  type(node_t), dimension(:), intent(in) :: node
  type(map_t),  dimension(:), intent(in) :: node_map
  type(edge_t), dimension(:), intent(in) :: edge
  type(map_t),  dimension(:), intent(in) :: edge_map
  !
  !.. Local Scalars ..
  integer :: i,ierr,n,n1,n2,nep,nvar,scnt,rcnt,rmin,rmax
  integer :: my_part,adj_part,adj_cpu,npair,this_pair
  integer :: lnode,loc_node,glb_node,adj_node,node_part
  integer :: ledge,loc_edge,glb_edge,adj_edge,edge_part
  integer(int_mpi) :: this_edgpt,tag12,np
  integer(int_mpi) :: block_size,tag21
  real(wp) :: rd,rn
  !
  !.. Local Arrays ..
  integer, dimension(1:ncpu) :: rnmin
  integer, dimension(1:ncpu) :: rnmax
  integer, dimension(1:ncpu) :: remin
  integer, dimension(1:ncpu) :: remax
  integer, dimension(1:ncpu) :: cpu_nums
  integer, dimension(1:ncpu) :: cpu_disp
  logical, dimension(1:ncpu) :: cpu_msk
  !
  !.. Local Allocatable Arrays ..
  logical,          allocatable :: send_node_msk(:)
  integer,          allocatable :: recv_node_map(:)
  logical,          allocatable :: send_edge_msk(:)
  integer,          allocatable :: recv_edge_map(:)
  integer,          allocatable :: comm_pair(:,:)
  integer(int_mpi), allocatable :: send_disp(:)
  integer(int_mpi), allocatable :: recv_disp(:)
  integer(int_mpi), allocatable :: edgpt_offset(:)
  integer(int_mpi), allocatable :: nodpt_offset(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_node_and_edge_datatypes"
  !
continue
  !
  my_part = mypnum + 1
  !
  lnode = size( node_map(my_part)%loc_to_glb )
  ledge = size( edge_map(my_part)%loc_to_glb )
  !
  ! Count the number of nodes located on multiple processors
  !
  n = 0
  do loc_node = 1,lnode
    glb_node = node_map(my_part)%loc_to_glb(loc_node)
    if (allocated(node(glb_node)%cpu)) then
      if (size(node(glb_node)%cpu) > 1) then
        n = n + 1
      end if
    end if
  end do
  !
  ! Allocate the bnd_nodes array
  !
  allocate ( bnd_nodes(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bnd_nodes",1,ierr,error_message)
  !
  ! Now fill in the components within the bnd_nodes array
  !
  n = 0
  fill_bnd_nodes: do loc_node = 1,lnode
    !
    glb_node = node_map(my_part)%loc_to_glb(loc_node)
    !
    if (allocated(node(glb_node)%cpu)) then
      if (size(node(glb_node)%cpu) > 1) then
        !
        n = n + 1
        !
        bnd_nodes(n)%idx = loc_node
        !
        if (my_part == node(glb_node)%controlling_partition) then
          bnd_nodes(n)%my_responsibility = true
        end if
        !
        rd = real( sum( node(glb_node)%cpu(:)%num_cells ) , kind=wp )
        !
        do i = 1,size(node(glb_node)%cpu)
          if (my_part == node(glb_node)%cpu(i)%partition) then
            rn = real( node(glb_node)%cpu(i)%num_cells , kind=wp ) 
            bnd_nodes(n)%ave_correction = rn / rd
            cycle fill_bnd_nodes
          end if
        end do
        !
      end if
    end if
    !
  end do fill_bnd_nodes
  !
  ! Count the number of edges located on multiple processors
  !
  n = 0
  do loc_edge = 1,ledge
    glb_edge = edge_map(my_part)%loc_to_glb(loc_edge)
    if (allocated(edge(glb_edge)%cpu)) then
      if (size(edge(glb_edge)%cpu) > 1) then
        n = n + 1
      end if
    end if
  end do
  !
  ! Allocate the bnd_edges array
  !
  allocate ( bnd_edges(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bnd_edges",1,ierr,error_message)
  !
  ! Now fill in the components within the bnd_edges array
  !
  n = 0
  fill_bnd_edges: do loc_edge = 1,ledge
    !
    glb_edge = edge_map(my_part)%loc_to_glb(loc_edge)
    !
    if (allocated(edge(glb_edge)%cpu)) then
      if (size(edge(glb_edge)%cpu) > 1) then
        !
        n = n + 1
        !
        bnd_edges(n)%idx = loc_edge
        !
        if (my_part == edge(glb_edge)%controlling_partition) then
          bnd_edges(n)%my_responsibility = true
        end if
        !
        rd = real( sum( edge(glb_edge)%cpu(:)%num_cells ) , kind=wp )
        !
        do i = 1,size(edge(glb_edge)%cpu)
          if (my_part == edge(glb_edge)%cpu(i)%partition) then
            rn = real( edge(glb_edge)%cpu(i)%num_cells , kind=wp ) 
            bnd_edges(n)%ave_correction = rn / rd
            cycle fill_bnd_edges
          end if
        end do
        !
      end if
    end if
    !
  end do fill_bnd_edges
  !
  ! Number of variables being stored for each node or edge point
  !
  nvar = nq * (nr + 1)
  !
  ! Allocate the node and edge point offset array
  !
  allocate ( nodpt_offset(1:lnode) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodpt_offset",1,ierr,error_message)
  !
  allocate ( edgpt_offset(1:ledge) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"edgpt_offset",1,ierr,error_message)
  !
  ! Compute the total size of the buffer for
  ! communicating node/edge solutions and gradients
  !
  ! Node solution points are 1-to-1 with the grid nodes so initialize
  ! n_node_edge_pts to just nvar times the number of grid nodes located
  ! on partition boundaries
  !
  n_node_edge_pts = nvar * size(bnd_nodes)
  !
  ! Loop through the edges on partition boundaries and add nvar times
  ! the number of edge points for each edge to n_node_edge_pts
  !
  do n = 1,size(bnd_edges)
    !
    loc_edge = bnd_edges(n)%idx
    glb_edge = edge_map(my_part)%loc_to_glb(loc_edge)
    !
    !
    if (use_edgpts) then
      !
      do i = 1,size(edge(glb_edge)%cpu)
        if (edge(glb_edge)%cpu(i)%partition == my_part) then
          nep = size(edge(glb_edge)%cpu(i)%edgpts)
          n_node_edge_pts = n_node_edge_pts +  nvar*nep
        end if
      end do
      !
    else
      !
      nep = Cell_Solpts( Geom_Edge , edge(glb_edge)%order )
      n_node_edge_pts = n_node_edge_pts +  nvar*nep
      !
    end if
    !
  end do
  !
  ! Initialize the min/max edge/node indices
  !
  rnmin =  huge(0)
  remin =  huge(0)
  rnmax = -huge(0)
  remax = -huge(0)
  !
  ! Create an array that identifies all processors that the current processor
  ! needs to communicate with and are of a higher CPU rank
  !
  cpu_msk = fals
  !
  ! First, find the communication required for nodes while also creating an
  ! offset array for the number of node points that lie on a partition boundary
  !
  do n = 1,size(bnd_nodes)
    !
    loc_node = bnd_nodes(n)%idx
    glb_node = node_map(my_part)%loc_to_glb(loc_node)
    !
    do i = 1,size(node(glb_node)%cpu)
      !
      ! Mark this partition for communication
      !
      cpu_msk( node(glb_node)%cpu(i)%partition ) = true
      !
      node_part = node(glb_node)%cpu(i)%partition
      !
      if (node_part == my_part) then
        !
        ! This is the node for this processor so
        ! create the node point offset for this node
        !
        nodpt_offset(loc_node) = 1_int_mpi
        !
      else
        !
        ! This is a node is on an adjacent processor so update the
        ! min and max node index used to allocate recv_node_map
        !
        rnmin(node_part) = min( rnmin(node_part) , &
                                node(glb_node)%cpu(i)%loc_node )
        rnmax(node_part) = max( rnmax(node_part) , &
                                node(glb_node)%cpu(i)%loc_node )
        !
      end if
      !
    end do
    !
  end do
  !
  ! Reshuffle nodpt_offset by summing the total number of nodes for
  ! all the communication nodes before the current one
  !
  do n = size(bnd_nodes),1,-1
    !
    loc_node = bnd_nodes(n)%idx
    !
    nodpt_offset(loc_node) = sum( nodpt_offset(1:loc_node-1) )
    !
  end do
  !
  if (maxval(nodpt_offset,dim=1) /= size(bnd_nodes)) then
    write (error_message,1) "Something is wrong here with the node offset"
    call stop_recon(stop_mpi,pname,error_message)
  end if
  !
  ! Second, find the communication required for edges while also creating an
  ! offset array for the number of edge points for each edge that lies on a
  ! partition boundary
  !
  do n = 1,size(bnd_edges)
    !
    loc_edge = bnd_edges(n)%idx
    glb_edge = edge_map(my_part)%loc_to_glb(loc_edge)
    !
    do i = 1,size(edge(glb_edge)%cpu)
      !
      ! Mark this partition for communication
      !
      cpu_msk( edge(glb_edge)%cpu(i)%partition ) = true
      !
      edge_part = edge(glb_edge)%cpu(i)%partition
      !
      if (edge_part == my_part) then
        !
        ! This is the edge for this processor so
        ! create the edge point offset for this edge
        !
        nep = Cell_Solpts( Geom_Edge , edge(glb_edge)%order )
        edgpt_offset(loc_edge) = int( nep , kind=int_mpi )
        !
      else
        !
        ! This is an edge on an adjacent processor so update the
        ! min and max edge index used to allocate recv_edge_map
        !
        remin(edge_part) = min( remin(edge_part) , &
                                edge(glb_edge)%cpu(i)%loc_edge )
        remax(edge_part) = max( remax(edge_part) , &
                                edge(glb_edge)%cpu(i)%loc_edge )
        !
      end if
      !
    end do
    !
  end do
  !
  ! Reshuffle edgpt_offset by summing the total number of edge points for
  ! all the communication edges before the current one
  !
  do n = size(bnd_edges),1,-1
    !
    loc_edge = bnd_edges(n)%idx
    !
    edgpt_offset(loc_edge) = sum( edgpt_offset(1:loc_edge-1) )
    !
  end do
  !
  ! Finally, add the total node point offset into
  ! the offset for each of the edge points
  !
  edgpt_offset = edgpt_offset + size(bnd_nodes,kind=int_mpi)
  !
  ! Count the number of communication pairs of the current processor
  ! with all other processors that have a higher rank
  !
  cpu_nums = 0
  cpu_nums(my_part) = count( cpu_msk(my_part+1:ncpu) )
  !
  ! Collect the number of communication pairs for all processors
  !
  call mpi_allgather(MPI_IN_PLACE,0,MPI_DATATYPE_NULL, &
                     cpu_nums,1,mpi_inttyp,mpi_comm_world,mpierr)
  !
  ! Count the total number of communication pairs for nodes and edges
  !
  npair = sum(cpu_nums)
  !
  ! Allocate the comm_pair array
  !
  allocate ( comm_pair(1:2,1:npair) , source=-1 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"comm_pair",1,ierr,error_message)
  !
  ! Create the displacement array for gathering the comm_pair array
  ! NOTE: cpu_nums is multiplied by 2 because each pair contains two
  !       pieces of information when collecting comm_pair
  !
  cpu_disp = 0
  do n = 2,ncpu
    cpu_disp(n) = cpu_disp(n-1) + cpu_nums(n-1)*size(comm_pair,dim=1)
  end do
  !
  ! Fill in the section of comm_pair for the
  ! current processor before communicating
  !
  ! beginning index for the second dim of comm_pair for the current processor
  n1 = cpu_disp(my_part)+1
  ! ending    index for the second dim of comm_pair for the current processor
  n2 = cpu_disp(my_part)+cpu_nums(my_part)
  !
  ! processor of lower  rank that belongs to this communication pair
  comm_pair(1,n1:n2) = mypnum
  ! processor of higher rank that belongs to this communication pair
  comm_pair(2,n1:n2) = pack( intseq(mypnum+1,ncpu-1) , &
                             mask=cpu_msk(my_part+1:ncpu) )
  !
  ! Collect the individual comm_pair arrays for all processors
  !
  call mpi_allgatherv(MPI_IN_PLACE,0,MPI_DATATYPE_NULL, &
                      comm_pair,cpu_nums,cpu_disp,mpi_inttyp, &
                      mpi_comm_world,mpierr)
  !
  ! Allocate the array containing the communication information
  !
  allocate ( exch_neusp(1:npair) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"exch_neusp",1,ierr,error_message)
  !
  ! Fill in the adj_cpu and tag parts of exch_neusp
  !
  do this_pair = 1,size(exch_neusp)
    !
    tag12 = int( this_pair                    , kind=int_mpi )
    tag21 = int( this_pair + size(exch_neusp) , kind=int_mpi )
    !
    if (comm_pair(1,this_pair) == mypnum) then
      !
      exch_neusp(this_pair)%adj_cpu  = int(comm_pair(2,this_pair),kind=int_mpi)
      exch_neusp(this_pair)%send_tag = tag12
      exch_neusp(this_pair)%recv_tag = tag21
     !exch_neusp(this_pair)%send_tag = mypnum
     !exch_neusp(this_pair)%recv_tag = int(comm_pair(2,this_pair),kind=int_mpi)
      !
    else if (comm_pair(2,this_pair) == mypnum) then
      !
      exch_neusp(this_pair)%adj_cpu  = int(comm_pair(1,this_pair),kind=int_mpi)
      exch_neusp(this_pair)%send_tag = tag21
      exch_neusp(this_pair)%recv_tag = tag12
     !exch_neusp(this_pair)%send_tag = mypnum
     !exch_neusp(this_pair)%recv_tag = int(comm_pair(1,this_pair),kind=int_mpi)
      !
    end if
    !
  end do
  !
  ! Deallocate comm_pair since it is no longer needed
  !
  deallocate ( comm_pair , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"comm_pair",2,ierr,error_message)
  !
  ! Allocate all the arrays needed for creating
  ! the send and receive node and edge maps
  !
  allocate ( send_node_msk(1:lnode) , source=fals , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"send_node_msk",1,ierr,error_message)
  !
  allocate ( send_edge_msk(1:ledge) , source=fals , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"send_edge_msk",1,ierr,error_message)
  !
  ! Loop through each of the node/edge communication pairs and create the
  ! send and receive MPI datatypes for the two processors of this pair
  !
  pair_loop: do this_pair = 1,size(exch_neusp)
    !
    if (exch_neusp(this_pair)%adj_cpu < 0_int_mpi) cycle pair_loop
    !
    adj_cpu = exch_neusp(this_pair)%adj_cpu
    !
    adj_part = adj_cpu + 1
    !
    ! Reset the mask arrays for this communication pair
    !
    send_node_msk(:) = fals
    send_edge_msk(:) = fals
    !
    ! Allocate the mapping arrays for this communication pair
    !
    rmin = rnmin(adj_part)
    rmax = rnmax(adj_part)
    allocate ( recv_node_map(rmin:rmax) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_node_map",1,ierr,error_message)
    !
    rmin = remin(adj_part)
    rmax = remax(adj_part)
    allocate ( recv_edge_map(rmin:rmax) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_edge_map",1,ierr,error_message)
    !
    ! Loop to find the communication boundary nodes between the
    ! current processor and the current exchange partner, adj_cpu.
    ! 1. We will create a mask array, send_node_msk, to identify all
    !    the local boundary nodes that need to be sent to adj_cpu.
    ! 2. We will also create an array, recv_node_map, that maps
    !    the boundary nodes of adj_cpu to the local boundary nodes
    !    whose data needs to be received from adj_cpu.
    !
    find_nodes: do n = 1,size(bnd_nodes)
      !
      loc_node = bnd_nodes(n)%idx
      glb_node = node_map(my_part)%loc_to_glb(loc_node)
      !
      find_nodes_find_cpu: do i = 1,size(node(glb_node)%cpu)
        !
        if (node(glb_node)%cpu(i)%partition == adj_part) then
          !
          adj_node = node(glb_node)%cpu(i)%loc_node
          send_node_msk(loc_node) = true
          recv_node_map(adj_node) = loc_node
          !
          exit find_nodes_find_cpu
          !
        end if
        !
      end do find_nodes_find_cpu
      !
    end do find_nodes
    !
    ! Loop to find the communication boundary edges between the
    ! current processor and the current exchange partner, adj_cpu.
    ! 1. We will also create a mask array, send_edge_msk, to identify
    !    all the local boundary edges that need to be sent to adj_cpu.
    ! 2. We will also create an array, recv_edge_map, that maps
    !    the boundary edges of adj_cpu to the local boundary edges
    !    whose data needs to be received from adj_cpu.
    !
    find_edges: do n = 1,size(bnd_edges)
      !
      loc_edge = bnd_edges(n)%idx
      glb_edge = edge_map(my_part)%loc_to_glb(loc_edge)
      !
      find_edges_find_cpu: do i = 1,size(edge(glb_edge)%cpu)
        !
        if (edge(glb_edge)%cpu(i)%partition == adj_part) then
          !
          adj_edge = edge(glb_edge)%cpu(i)%loc_edge
          send_edge_msk(loc_edge) = true
          recv_edge_map(adj_edge) = loc_edge
          !
          exit find_edges_find_cpu
          !
        end if
        !
      end do find_edges_find_cpu
      !
    end do find_edges
    !
    ! Now create the send and receive maps. These maps will contain,
    ! in order, the faces that need to be sent to and received from the
    ! current adjacent partition adj_cpu
    ! NOTE: Each processor will send its node and edge data in the order of
    !       the local node and edge IDs. The receiving datatypes need to
    !       account for this so that the received data is placed in the
    !       correct location within an array. This is why we have created
    !       the arrays recv_node_map and recv_edge_map so that we can loop,
    !       in order, through the nodes and edges of the adjacent partition
    !       to create the receive offsets in the same order of the adjacent
    !       partitions edge and node IDs.
    !
    ! Create the send node maps
    !
    !
    ! Count how many node and edge points in total will need to be communicated.
    ! NOTE: This should be the same for both sending and receiving
    !
    ! Initialize the counter to the number of node points
    ! that need to be communicated
    !
    scnt = count( send_node_msk )
    !
    ! Loop through the edges and add to the counter the number of edge points
    ! for each edge that needs to be communicated
    !
    count_edges: do loc_edge = 1,ledge
      !
      if (.not. send_edge_msk(loc_edge)) cycle count_edges
      !
      glb_edge = edge_map(my_part)%loc_to_glb(loc_edge)
      !
      if (use_edgpts) then
        !
        count_edges_find_cpu: do i = 1,size(edge(glb_edge)%cpu)
          if (edge(glb_edge)%cpu(i)%partition == my_part) then
            scnt = scnt + size( edge(glb_edge)%cpu(i)%edgpts )
            exit count_edges_find_cpu
          end if
        end do count_edges_find_cpu
        !
      else
        !
        scnt = scnt + Cell_Solpts( Geom_Edge , edge(glb_edge)%order )
        !
      end if
      !
    end do count_edges
    !
    ! Allocate the send and recv arrays needed to create the MPI datatypes now
    ! that we know the total number of points that will need to be communicated.
    !
    allocate ( send_disp(1:scnt) , source=0_int_mpi , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"send_disp",1,ierr,error_message)
    !
    allocate ( recv_disp(1:scnt) , source=0_int_mpi , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_disp",1,ierr,error_message)
    !
    !###########################################################################
    !###########################################################################
    ! Create the send displacement offsets for the edge and node points that
    ! need to be communicated for this communication pair so that we can create
    ! the MPI datatypes to easily exchange this data with adj_cpu
    !###########################################################################
    !###########################################################################
    !
    ! Create the send displacement offsets for the node points
    !
    scnt = 0
    node_send_disp: do loc_node = 1,lnode
      !
      if (.not. send_node_msk(loc_node)) cycle node_send_disp
      !
      scnt = scnt + 1
      !
      send_disp(scnt) = nodpt_offset(loc_node)
      !
    end do node_send_disp
    !
    ! Create the send displacement offsets for the edge points
    !
    edge_send_disp: do loc_edge = 1,ledge
      !
      if (.not. send_edge_msk(loc_edge)) cycle edge_send_disp
      !
      glb_edge = edge_map(my_part)%loc_to_glb(loc_edge)
      !
      if (use_edgpts) then
        !
        edge_send_disp_find_cpu: do i = 1,size(edge(glb_edge)%cpu)
          !
          if (edge(glb_edge)%cpu(i)%partition == my_part) then
            !
            do n = 1,size( edge(glb_edge)%cpu(i)%edgpts )
              !
              scnt = scnt + 1
              !
              this_edgpt = int( edge(glb_edge)%cpu(i)%edgpts(n) , kind=int_mpi )
              !
              send_disp(scnt) = edgpt_offset(loc_edge) + this_edgpt
              !
            end do
            !
            exit edge_send_disp_find_cpu
            !
          end if
          !
        end do edge_send_disp_find_cpu
        !
      else
        !
        nep = Cell_Solpts( Geom_Edge , edge(glb_edge)%order )
        !
        do n = 0,nep-1
          !
          scnt = scnt + 1
          !
          this_edgpt = int( n , kind=int_mpi )
          !
          send_disp(scnt) = edgpt_offset(loc_edge) + this_edgpt
          !
        end do
        !
      end if
      !
    end do edge_send_disp
    !
    !###########################################################################
    !###########################################################################
    ! Create the receive displacement offsets for the edge and node points that
    ! need to be communicated for this communication pair so that we can create
    ! the MPI datatypes to easily exchange this data with adj_cpu
    !###########################################################################
    !###########################################################################
    !
    ! Initialize the receive counter to zero
    !
    rcnt = 0
    !
    ! Create the receive displacement offsets for the node points
    !
    node_recv_disp: do adj_node = lbound(recv_node_map,dim=1), &
                                  ubound(recv_node_map,dim=1)
      !
      if (recv_node_map(adj_node) == 0) cycle node_recv_disp
      !
      loc_node = recv_node_map(adj_node)
      !
      rcnt = rcnt + 1
      !
      recv_disp(rcnt) = nodpt_offset(loc_node)
      !
    end do node_recv_disp
    !
    ! Create the receive displacement offsets for the edge points
    !
    edge_recv_disp: do adj_edge = lbound(recv_edge_map,dim=1), &
                                  ubound(recv_edge_map,dim=1)
      !
      if (recv_edge_map(adj_edge) == 0) cycle edge_recv_disp
      !
      loc_edge = recv_edge_map(adj_edge)
      !
      glb_edge = edge_map(my_part)%loc_to_glb(loc_edge)
      !
      if (use_edgpts) then
        !
        edge_recv_disp_find_cpu: do i = 1,size(edge(glb_edge)%cpu)
          !
          if (edge(glb_edge)%cpu(i)%partition == my_part) then
            !
            do n = 1,size( edge(glb_edge)%cpu(i)%edgpts )
              !
              rcnt = rcnt + 1
              !
              this_edgpt = int( edge(glb_edge)%cpu(i)%edgpts(n) , kind=int_mpi )
              !
              recv_disp(rcnt) = edgpt_offset(loc_edge) + this_edgpt
              !
            end do
            !
            exit edge_recv_disp_find_cpu
            !
          end if
          !
        end do edge_recv_disp_find_cpu
        !
      else
        !
        nep = Cell_Solpts( Geom_Edge , edge(glb_edge)%order )
        !
        do n = 0,nep-1
          !
          rcnt = rcnt + 1
          !
          this_edgpt = int( n , kind=int_mpi )
          !
          recv_disp(rcnt) = edgpt_offset(loc_edge) + this_edgpt
          !
        end do
        !
      end if
      !
    end do edge_recv_disp
    !
    ! Deallocate recv_node_map and recv_edge_map since they are
    ! no longer needed for the current communication pair
    !
    deallocate ( recv_node_map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_node_map",2,ierr,error_message)
    !
    deallocate ( recv_edge_map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_edge_map",2,ierr,error_message)
    !
    !##########################################################
    !##########################################################
    ! Create the MPI datatypes for exchanging node and edge
    ! solution and gradient data between adjacent partitions
    !##########################################################
    !##########################################################
    !
    ! The block size is the number of data values for each node or edge point
    ! i.e., nq + nr*nq
    !
    block_size = int( nvar , kind=int_mpi )
    !
    ! Adjust send_disp and recv_disp to account for the block size
    !
    send_disp = send_disp * block_size
    recv_disp = recv_disp * block_size
    !
    ! Create and commit the send MPI datatype for this communication pair
    !
    np = size( send_disp , kind=int_mpi )
    call mpi_type_create_indexed_block(np,block_size,send_disp,mpi_flttyp, &
                                       exch_neusp(this_pair)%send_datatype, &
                                       mpierr)
    call mpi_type_commit(exch_neusp(this_pair)%send_datatype,mpierr)
    !
    ! Create and commit the receive MPI datatype for this communication pair
    !
    np = size( recv_disp , kind=int_mpi )
    call mpi_type_create_indexed_block(np,block_size,recv_disp,mpi_flttyp, &
                                       exch_neusp(this_pair)%recv_datatype, &
                                       mpierr)
    call mpi_type_commit(exch_neusp(this_pair)%recv_datatype,mpierr)
    !
    ! Deallocate the arrays for this communication
    ! pair before continuing to the next pair
    !
    deallocate ( send_disp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"send_disp",2,ierr,error_message)
    !
    deallocate ( recv_disp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_disp",2,ierr,error_message)
    !
  end do pair_loop
  !
  ! Deallocate send_node_msk, send_edge_msk, nodpt_offset,
  ! and edgpt_offset before we leave
  !
  deallocate ( send_node_msk , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"send_node_msk",2,ierr,error_message)
  !
  deallocate ( send_edge_msk , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"send_edge_msk",2,ierr,error_message)
  !
  deallocate ( nodpt_offset , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodpt_offset",2,ierr,error_message)
  !
  deallocate ( edgpt_offset , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"edgpt_offset",2,ierr,error_message)
  !
  ! Format Statements
  !
  1 format (a)
  !
end subroutine create_node_and_edge_datatypes
!
!###############################################################################
!
subroutine exchange_faceusp(faceusp)
  !
  !.. Use Statements ..
#ifndef USE_TYPEBOUND_FACE_PACK
  use geovar,  only : bface
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
  do nf = 1,size(bface,dim=2)
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
  use geovar,  only : bface
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
  do nf = 1,size(bface,dim=2)
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
  use geovar,  only : bface
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
  do nf = 1,size(bface,dim=2)
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
  use geovar,  only : bface
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
  do nf = 1,size(bface,dim=2)
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
  integer  :: n,n1,n2,n3,nsz
  integer  :: i,ibeg,iend,ierr
  integer  :: this_pair
  real(wp) :: ave_correction
  !
  !.. MPI integers ..
  integer(int_mpi) :: adj_cpu
  integer(int_mpi) :: send_tag,send_dtype
  integer(int_mpi) :: recv_tag,recv_dtype
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
  ! Allocate the send and receive buffer arrays
  !
  allocate ( send_buffer(1:n_node_edge_pts) , source=0.0_r8 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"send_buffer",1,ierr,error_message)
  !
  allocate ( recv_buffer(1:n_node_edge_pts) , source=0.0_r8 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_buffer",1,ierr,error_message)
  !
  ! Create the send buffer
  !
  iend = 0
  !
  ! Gather and pack up the node solutions and gradients
  ! NOTE: This data will be stored in the 1D send_buffer array such that
  !       it would be equivalent to send_buffer(1:nr+1,1:nq,:) where
  !       send_buffer(1,1:nq,:) are the conserved solution variables and
  !       send_buffer(2:nr+1,1:nq,:) are the solution gradients
  !       
  pack_nodes: do i = 1,size(bnd_nodes)
    !
    n = bnd_nodes(i)%idx
    !
    ! First correct the values in nodeusp and nodedusp since the averaging
    ! coefficient used to compute the node values was based on the number
    ! of cells on the local partition that contribute to this node. If the
    ! number of local cells containing this node is lc, these values are
    ! currently averaged using the coefficient (1/lc). If the number of global
    ! cells containing this node is gc, this averaging correction is simply
    ! multiplying these nodal values by (lc/gc) resulting in the final
    ! averaging coefficient of (1/gc).
    !
    ave_correction = bnd_nodes(i)%ave_correction
    !
    nodeusp(n)%v(:)    = ave_correction * nodeusp(n)%v(:)
    nodedusp(n)%v(:,:) = ave_correction * nodedusp(n)%v(:,:)
    !
    ! Pack the conserved solution variables for this node point
    !
    nsz = size(nodeusp(n)%v(:))
    ibeg = iend + 1
    iend = iend + nsz
    send_buffer(ibeg:iend) = reshape( real(nodeusp(n)%v(:),kind=r8) , [nsz] )
    !
    ! Pack the solution gradients for this node point
    !
    nsz = size(nodedusp(n)%v(:,:))
    ibeg = iend + 1
    iend = iend + nsz
    send_buffer(ibeg:iend) = &
                 reshape( real(nodedusp(n)%v(:,:),kind=r8) , [nsz] )
    !
  end do pack_nodes
  !
  ! Gather and pack up the edge solutions and gradients
  ! NOTE: This data will be stored in the 1D send_buffer array such that
  !       it would be equivalent to send_buffer(1:nr+1,1:nq,:) where
  !       send_buffer(1,1:nq,:) are the conserved solution variables and
  !       send_buffer(2:nr+1,1:nq,:) are the solution gradients
  !
  pack_edges: do i = 1,size(bnd_edges)
    !
    n = bnd_edges(i)%idx
    !
    ! First correct the values in edgeusp and edgedusp since the averaging
    ! coefficient used to compute the edge values was based on the number
    ! of cells on the local partition that contribute to this edge. If the
    ! number of local cells containing this edge is lc, these values are
    ! currently averaged using the coefficient (1/lc). If the number of global
    ! cells containing this edge is gc, this averaging correction is simply
    ! multiplying these nodal values by (lc/gc) resulting in the final
    ! averaging coefficient of (1/gc).
    !
    ave_correction = bnd_edges(i)%ave_correction
    !
    edgeusp(n)%v(:,:)    = ave_correction * edgeusp(n)%v(:,:)
    edgedusp(n)%v(:,:,:) = ave_correction * edgedusp(n)%v(:,:,:)
    !
    ! Pack the conserved solution variables for these edge points
    !
    nsz = size(edgeusp(n)%v(:,:))
    ibeg = iend + 1
    iend = iend + nsz
    send_buffer(ibeg:iend) = reshape( real(edgeusp(n)%v(:,:),kind=r8) , [nsz] )
    !
    ! Pack the solution gradients for these edge points
    !
    nsz = size(edgedusp(n)%v(:,:,:))
    ibeg = iend + 1
    iend = iend + nsz
    send_buffer(ibeg:iend) = &
                 reshape( real(edgedusp(n)%v(:,:,:),kind=r8) , [nsz] )
    !
  end do pack_edges
  !
  ! Loop through all the communication pairs and exchange solutions
  ! for the nodes and edge points, adding the received values to the
  ! current solutions before going to the next pair
  !
  exchange_loop: do this_pair = 1,size(exch_neusp)
    !
    ! If exch_neusp(this_pair)%adj_cpu is not the default value, this processor
    ! is one of the two processors making up this communication pair.
    !
    if (exch_neusp(this_pair)%adj_cpu >= 0) then
      !
      adj_cpu = exch_neusp(this_pair)%adj_cpu
      !
      send_tag = exch_neusp(this_pair)%send_tag
      recv_tag = exch_neusp(this_pair)%recv_tag
      !
      send_dtype = exch_neusp(this_pair)%send_datatype
      recv_dtype = exch_neusp(this_pair)%recv_datatype
      !
      ! Use MPI_SENDRECV to concurrently send and receive the necessary
      ! solution data from the adjacent processor. This is essentially
      ! the same as processor 1 calling MPI_SEND then MPI_RECV while
      ! processor 2 calls MPI_RECV then MPI_SEND.
      !
      call mpi_sendrecv(send_buffer,1,send_dtype,adj_cpu,send_tag, &
                        recv_buffer,1,recv_dtype,adj_cpu,recv_tag, &
                        mpi_comm_world,istat,mpierr)
      !
      ! Add the received data to the solutions at the nodes and edges
      ! NOTE: This will unfortunately add a bunch of zeros because of
      !       how its implemented
      !
      iend = 0
      !
      ! Unpack the received node solutions and gradients
      !
      unpack_nodes: do i = 1,size(bnd_nodes)
        !
        n = bnd_nodes(i)%idx
        !
        ! Unpack the conserved solution variables for this node point
        ! and add them to the local node values
        !
        n1 = size(nodeusp(n)%v,dim=1)
        ibeg = iend + 1
        iend = iend + n1
        nodeusp(n)%v(:) = nodeusp(n)%v(:) + &
                   reshape( real(recv_buffer(ibeg:iend),kind=wp) , [n1] )
        !
        ! Unpack the solution gradients for this node point
        ! and add them to the local node values
        !
        n1 = size(nodedusp(n)%v,dim=1)
        n2 = size(nodedusp(n)%v,dim=2)
        ibeg = iend + 1
        iend = iend + n1*n2
        nodedusp(n)%v(:,:) = nodedusp(n)%v(:,:) + &
                   reshape( real(recv_buffer(ibeg:iend),kind=wp) , [n1,n2] )
        !
      end do unpack_nodes
      !
      ! Unpack the received edge solutions and gradients
      !
      unpack_edges: do i = 1,size(bnd_edges)
        !
        n = bnd_edges(i)%idx
        !
        ! Unpack the conserved solution variables for these edge points
        ! and add them to the local edge values
        !
        n1 = size(edgeusp(n)%v,dim=1)
        n2 = size(edgeusp(n)%v,dim=2)
        ibeg = iend + 1
        iend = iend + n1*n2
        edgeusp(n)%v(:,:) = edgeusp(n)%v(:,:) + &
                   reshape( real(recv_buffer(ibeg:iend),kind=wp) , [n1,n2] )
        !
        ! Unpack the solution gradients for these edge points
        ! and add them to the local edge values
        !
        n1 = size(edgedusp(n)%v,dim=1)
        n2 = size(edgedusp(n)%v,dim=2)
        n3 = size(edgedusp(n)%v,dim=3)
        ibeg = iend + 1
        iend = iend + n1*n2*n3
        edgedusp(n)%v(:,:,:) = edgedusp(n)%v(:,:,:) + &
                   reshape( real(recv_buffer(ibeg:iend),kind=wp) , [n1,n2,n3] )
        !
      end do unpack_edges
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
subroutine exch_connectivity(nodeusp,nodedusp,edgeusp,edgedusp,faceusp)
  !
  !.. Use Statements ..
  use flowvar, only : ne1d_t,ne2d_t,ne3d_t,face_var
  !
  !.. Formal Arguments ..
  type(ne1d_t),   dimension(:), intent(inout) :: nodeusp
  type(ne2d_t),   dimension(:), intent(inout) :: nodedusp
  type(ne2d_t),   dimension(:), intent(inout) :: edgeusp
  type(ne3d_t),   dimension(:), intent(inout) :: edgedusp
  type(face_var), dimension(:), intent(inout) :: faceusp
  !
  !.. Local Scalars..
  integer  :: n,n1,n2,k
  integer  :: i,ibeg,iend,ierr
  integer  :: this_pair
  !
  !.. MPI integers ..
  integer(int_mpi) :: adj_cpu,mesg_tag
  integer(int_mpi) :: send_cpu,send_tag,send_dtype
  integer(int_mpi) :: recv_cpu,recv_tag,recv_dtype
  integer(int_mpi) :: istat(1:mpi_status_size)
  !
  !.. Local Allocatable Arrays ..
  real(r8), allocatable, dimension(:) :: send_buffer
  real(r8), allocatable, dimension(:) :: recv_buffer
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "exch_connectivity"
  !
continue
  !
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  ! #####                                                                #####
  ! #####  1) Exchange the point indices for the duplicated node and     #####
  ! #####     edge points on partition boundaries in order to establish  #####
  ! #####     the global cell subconnectivity thats needed to create     #####
  ! #####     the CGNS solution files.                                   #####
  ! #####                                                                #####
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  !
  ! Allocate the send and receive buffer arrays. Initialize these arrays
  ! to +1 so the correct point index from the responsible processor can
  ! be identified by making their point index negative
  !
  allocate ( send_buffer(1:n_node_edge_pts) , &
             source=real(one,kind=kind(send_buffer)) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"send_buffer",1,ierr,error_message)
  !
  allocate ( recv_buffer(1:n_node_edge_pts) , &
             source=real(one,kind=kind(recv_buffer)) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_buffer",1,ierr,error_message)
  !
  !##################################################################
  !##################################################################
  !###                                                            ###
  !###  1.1) Create the send buffer for the node and edge points  ###
  !###                                                            ###
  !##################################################################
  !##################################################################
  !
  iend = 0
  !
  ! 1.1.1) Gather and pack up the node points for
  !        the nodes on boundary partitions
  !       
  pack_nodes: do i = 1,size(bnd_nodes)
    !
    n = bnd_nodes(i)%idx
    !
    ibeg = iend + 1
    iend = iend + size(nodeusp(n)%v(:))
    !
    ! Make the index for the node point negative if this partition
    ! is responsible for the point index of this node
    !
    if (bnd_nodes(i)%my_responsibility) then
      send_buffer(ibeg) = real( -nodeusp(n)%v(1) , kind=kind(send_buffer) )
    end if
    !
    ! There is no information to pack up from the nodedusp array, but we do
    ! need to update the ibeg and iend indices so we can continue to pack up
    ! the next node point from nodeusp and store it in the correct location
    ! in send_buffer
    !
    ibeg = iend + 1
    iend = iend + size(nodedusp(n)%v)
    !
  end do pack_nodes
  !
  ! 1.1.2) Gather and pack up the edge points for
  !        the edges on partition boundaries
  !
  pack_edges: do i = 1,size(bnd_edges)
    !
    n = bnd_edges(i)%idx
    !
    n1 = size(edgeusp(n)%v,dim=1)
    n2 = size(edgeusp(n)%v,dim=2)
    !
    do k = 1,n2
      !
      ibeg = iend + 1
      iend = iend + n1
      !
      ! Make the index for the edge point negative if this partition
      ! is responsible for the point indices of this edge
      !
      if (bnd_edges(i)%my_responsibility) then
        send_buffer(ibeg) = real( -edgeusp(n)%v(1,k) , kind=kind(send_buffer) )
      end if
      !
    end do
    !
    ! There is no information to pack up from the edgedusp array, but we do
    ! need to update the ibeg and iend indices so we can continue to pack up
    ! the edge points of the next edge from edgeusp and store them in the
    ! correct location in send_buffer
    !
    ibeg = iend + 1
    iend = iend + size(edgedusp(n)%v)
    !
  end do pack_edges
  !
  !#######################################################################
  !#######################################################################
  !###                                                                 ###
  !### 1.2) Loop through all the communication pairs and exchange the  ###
  !###      indices for the node and edge points, and save the index   ###
  !###      for each point if the received index value is negative     ###
  !###      which indicates the adjacent partition is responsible for  ###
  !###      the index of that point.                                   ###
  !###                                                                 ###
  !#######################################################################
  !#######################################################################
  !
  exchange_loop: do this_pair = 1,size(exch_neusp)
    !
    ! If exch_neusp(this_pair)%adj_cpu is not the default value, this processor
    ! is one of the two processors making up this communication pair.
    !
    if (exch_neusp(this_pair)%adj_cpu >= 0) then
      !
      adj_cpu = exch_neusp(this_pair)%adj_cpu
      !
      send_tag = exch_neusp(this_pair)%send_tag
      recv_tag = exch_neusp(this_pair)%recv_tag
      !
      send_dtype = exch_neusp(this_pair)%send_datatype
      recv_dtype = exch_neusp(this_pair)%recv_datatype
      !
      ! Use MPI_SENDRECV to concurrently send and receive the necessary
      ! solution data from the adjacent processor. This is essentially
      ! the same as processor 1 calling MPI_SEND then MPI_RECV while
      ! processor 2 calls MPI_RECV then MPI_SEND.
      !
      call mpi_sendrecv(send_buffer,1,send_dtype,adj_cpu,send_tag, &
                        recv_buffer,1,recv_dtype,adj_cpu,recv_tag, &
                        mpi_comm_world,istat,mpierr)
      !
      ! Add the received data to the solutions at the nodes and edges
      ! NOTE: This will unfortunately add a bunch of zeros because of
      !       how its implemented
      !
      iend = 0
      !
      ! Unpack the received node solutions and gradients
      !
      unpack_nodes: do i = 1,size(bnd_nodes)
        !
        n = bnd_nodes(i)%idx
        !
        ! Unpack the node point index for this node
        !
        ibeg = iend + 1
        iend = iend + size(nodeusp(n)%v)
        !
        ! If the value of recv_buffer(ibeg) is negative, the adjacent
        ! partition is responsible for the point index of this node
        ! so copy the negative value of this point index into nodeusp
        ! for the current node.
        ! NOTE: We need to keep this negative so that later we can identify
        !       that this processor is not responsible for writing the
        !       solution variables for this node point.
        !
        if (recv_buffer(ibeg) < real(zero,kind=kind(recv_buffer))) then
          nodeusp(n)%v(1) = real( recv_buffer(ibeg) , &
                                  kind=kind(nodeusp(n)%v) )
        end if
        !
        ! We dont need to unpack anything for nodedusp, but we do need
        ! to update the ibeg and iend indices so we can continue to unpack
        ! the next node point for nodeusp from the correct location in
        ! recv_buffer
        !
        ibeg = iend + 1
        iend = iend + size(nodedusp(n)%v)
        !
      end do unpack_nodes
      !
      ! Unpack the received edge points in edgeusp
      !
      unpack_edges: do i = 1,size(bnd_edges)
        !
        n = bnd_edges(i)%idx
        !
        ! Unpack the edge point indices for these edge points
        !
        n1 = size(edgeusp(n)%v,dim=1)
        n2 = size(edgeusp(n)%v,dim=2)
        !
        do k = 1,n2
          !
          ibeg = iend + 1
          iend = iend + n1
          !
          ! If the value of recv_buffer(ibeg) is negative, the adjacent
          ! partition is responsible for the point indices of this edge
          ! so copy the negative value of this point index into edgeusp
          ! for this edge point on the current edge.
          ! NOTE: We need to keep this negative so that later we can identify
          !       that this processor is not responsible for writing the
          !       solution variables for this edge point.
          !
          if (recv_buffer(ibeg) < real(zero,kind=kind(recv_buffer))) then
            edgeusp(n)%v(1,k) = real( recv_buffer(ibeg) , &
                                      kind=kind(edgeusp(n)%v) )
          end if
          !
        end do
        !
        ! We dont need to unpack anything for edgedusp, but we do need
        ! to update the ibeg and iend indices so we can continue to unpack
        ! the edge points of the next edge for edgeusp from the correct
        ! location in recv_buffer
        !
        ibeg = iend + 1
        iend = iend + size(edgedusp(n)%v)
        !
      end do unpack_edges
      !
      ! Reset the receiver buffer to one for the next communication
      !
      recv_buffer = real( one , kind=kind(recv_buffer) )
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
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  ! #####                                                                #####
  ! #####  2) Exchange the point indices for the duplicated face points  #####
  ! #####     on partition boundaries in order to establish the global   #####
  ! #####     cell subconnectivity thats needed to create the CGNS       #####
  ! #####     solution files.                                            #####
  ! #####                                                                #####
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  !
  ! Reset the fusp_send_buffer and fusp_recv_buffer arrays to +1
  !
  fusp_send_buffer(:) = real( one , kind=kind(fusp_send_buffer) )
  fusp_recv_buffer(:) = real( one , kind=kind(fusp_recv_buffer) )
  !
  ! 2.1) Pack up the face points
  !
  iend = 0
  pack_faces: do i = 1,size(bnd_faces)
    !
    n = bnd_faces(i)%idx
    !
    n1 = size(faceusp(n)%v,dim=1)
    n2 = size(faceusp(n)%v,dim=2)
    !
    do k = 1,n2
      !
      ibeg = iend + 1
      iend = iend + n1
      !
      ! Make the index for the face point negative if this partition
      ! is responsible for the point indices of this face
      !
      if (bnd_faces(i)%my_responsibility) then
        fusp_send_buffer(ibeg) = real( -faceusp(n)%v(1,k,1) , &
                                       kind=kind(fusp_send_buffer) )
      end if
      !
    end do
    !
  end do pack_faces
  !
  ! 2.2) Initialize all the receives so each process
  !      is ready once the send is posted
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
  ! 2.3) Send out the face flux point solution variables
  !      to the adjacent processes
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
  ! 2.4) Wait for all faceusp sends and receives to complete
  !
  call mpi_waitall(int(2*ncomm,kind=int_mpi),fusp_rqst,fusp_istat,mpierr)
  !
  ! 2.5) Update the right states for only the communication boundary faces
  !
  iend = 0
  unpack_faces: do i = 1,size(bnd_faces)
    !
    n = bnd_faces(i)%idx
    !
    n1 = size(faceusp(n)%v,dim=1)
    n2 = size(faceusp(n)%v,dim=2)
    !
    do k = 1,n2
      !
      ibeg = iend + 1
      iend = iend + n1
      !
      ! If the value of fusp_recv_buffer(ibeg) is negative, the adjacent
      ! partition is responsible for the point indices of this face
      ! so copy the negative value of this point index into faceusp
      ! for this face point on the current face.
      ! NOTE: We need to keep this negative so that later we can identify
      !       that this processor is not responsible for writing the
      !       solution variables for this face point.
      !
      if (fusp_recv_buffer(ibeg) < real(zero,kind=kind(fusp_recv_buffer))) then
        faceusp(n)%v(1,k,1) = real( fusp_recv_buffer(ibeg) , &
                                    kind=kind(faceusp(n)%v) )
      end if
      !
    end do
    !
  end do unpack_faces
  !
  ! The bnd_faces array should no longer be needed so deallocate before leaving
  !
  deallocate ( bnd_faces , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bnd_faces",2,ierr,error_message)
  !
end subroutine exch_connectivity
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
    np = size(xyz_nodes,dim=2)
    nc = size(epart)
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
!!if (allocated(node_map)) then
!!  !
!!  call dt%set( storage_size(node_map,kind=inttype) , &
!!                       size(node_map,kind=inttype) )
!! !write (iout,*)
!! !write (iout,*) "node_map(:)"
!! !write (iout,*) "size = ",size(node_map)
!! !write (iout,*) "storage_size = ",storage_size(node_map)
!! !!
!! !loc_siz = 0
!! !glb_siz = 0
!!  !
!!  do i = lbound(node_map,dim=1),ubound(node_map,dim=1)
!!    !
!!    if (allocated(node_map(i)%loc_to_glb)) then
!!      call dt%add( storage_size(node_map(i)%loc_to_glb,kind=inttype) , &
!!                           size(node_map(i)%loc_to_glb,kind=inttype) )
!! !    loc_siz = loc_siz + size(node_map(i)%loc_to_glb)
!!    end if
!!    !
!!  end do
!!  !
!! !i = lbound(node_map,dim=1)
!! !if (allocated(node_map(i)%loc_to_glb)) then
!! !  write (iout,*)
!! !  write (iout,*) "node_map(:)%loc_to_glb"
!! !  write (iout,*) "size = ",loc_siz
!! !  write (iout,*) "storage_size = ",storage_size(node_map(i)%loc_to_glb)
!! !end if
!!  !
!!  call total%add(dt)
!!  write (array_name,5) "node_map"
!!  write (iou,2) array_name,dt%get_units()
!!  !
!!else
!!  !
!!  write (iou,7) "node_map"
!!  !
!!end if
!!!
!!!-----------------------------------------------------------------------------
!!!
!!if (allocated(face_map)) then
!!  !
!!  call dt%set( storage_size(face_map,kind=inttype) , &
!!                       size(face_map,kind=inttype) )
!! !write (iout,*)
!! !write (iout,*) "face_map(:)"
!! !write (iout,*) "size = ",size(face_map)
!! !write (iout,*) "storage_size = ",storage_size(face_map)
!! !!
!! !loc_siz = 0
!! !glb_siz = 0
!!  !
!!  do i = lbound(face_map,dim=1),ubound(face_map,dim=1)
!!    !
!!    if (allocated(face_map(i)%loc_to_glb)) then
!!      call dt%add( storage_size(face_map(i)%loc_to_glb,kind=inttype) , &
!!                           size(face_map(i)%loc_to_glb,kind=inttype) )
!! !    loc_siz = loc_siz + size(face_map(i)%loc_to_glb)
!!    end if
!!    !
!!  end do
!!  !
!! !i = lbound(face_map,dim=1)
!! !if (allocated(face_map(i)%loc_to_glb)) then
!! !  write (iout,*)
!! !  write (iout,*) "face_map(:)%loc_to_glb"
!! !  write (iout,*) "size = ",loc_siz
!! !  write (iout,*) "storage_size = ",storage_size(face_map(i)%loc_to_glb)
!! !end if
!!  !
!!  call total%add(dt)
!!  write (array_name,5) "face_map"
!!  write (iou,2) array_name,dt%get_units()
!!  !
!!else
!!  !
!!  write (iou,7) "face_map"
!!  !
!!end if
!!!
!!!-----------------------------------------------------------------------------
!!!
!!if (allocated(face)) then
!!  !
!!  call dt%set( storage_size(face,kind=inttype) , &
!!                       size(face,kind=inttype) )
!!  !
!!  do i = lbound(face,dim=1),ubound(face,dim=1)
!!    !
!!    call dt%add( storage_size(face(i)%geom,kind=inttype) )
!!    call dt%add( storage_size(face(i)%order,kind=inttype) )
!!    !
!!    call dt%add( storage_size(face(i)%nodes,kind=inttype) , &
!!                         size(face(i)%nodes,kind=inttype) )
!!    !
!!    call dt%add( storage_size(face(i)%left,kind=inttype) )
!!    call dt%add( storage_size(face(i)%left%cell,kind=inttype) )
!!    call dt%add( storage_size(face(i)%left%cell_face,kind=inttype) )
!!    !
!!    if (allocated(face(i)%left%flxpts)) then
!!      call dt%add( storage_size(face(i)%left%flxpts,kind=inttype) , &
!!                           size(face(i)%left%flxpts,kind=inttype) )
!!    end if
!!    !
!!    call dt%add( storage_size(face(i)%right,kind=inttype) )
!!    call dt%add( storage_size(face(i)%right%cell,kind=inttype) )
!!    call dt%add( storage_size(face(i)%right%cell_face,kind=inttype) )
!!    !
!!    if (allocated(face(i)%right%flxpts)) then
!!      call dt%add( storage_size(face(i)%right%flxpts,kind=inttype) , &
!!                           size(face(i)%right%flxpts,kind=inttype) )
!!    end if
!!    !
!!  end do
!!  !
!!  call total%add(dt)
!!  write (array_name,5) "face"
!!  write (iou,2) array_name,dt%get_units()
!!  !
!!else
!!  !
!!  write (iou,7) "face"
!!  !
!!end if
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
subroutine create_dual_graph(cell_geom,cell_order,nodes_of_cell, &
                             nodes_of_cell_ptr,bface,xadj,adjncy,face)
  !
  !.. Use Statements ..
  use order_mod,        only : geom_solpts
  use connectivity_mod, only : GetFaceNodes_for_QuadFace
  use connectivity_mod, only : GetFaceNodes_for_TriaFace
  use connectivity_mod, only : RotateFluxPoints_Quad
  use connectivity_mod, only : RotateFluxPoints_Tria
  !
  !.. Formal Arguments ..
  integer,                   dimension(:),    intent(in) :: cell_geom
  integer,                   dimension(:),    intent(in) :: cell_order
  integer,                   dimension(:),    intent(in) :: nodes_of_cell
  integer,                   dimension(:),    intent(in) :: nodes_of_cell_ptr
  integer,                   dimension(:,:),  intent(in) :: bface
  integer,      allocatable, dimension(:), intent(inout) :: xadj
  integer,      allocatable, dimension(:), intent(inout) :: adjncy
  type(face_t), allocatable, dimension(:), intent(inout) :: face
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
  integer :: lnode,lfbnd,lface,lcell
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
  lnode = maxval(nodes_of_cell)
  lfbnd = size(bface,dim=2)
  lcell = size(nodes_of_cell_ptr)-1 - lfbnd
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Find the cells that surround each node
  !
  allocate ( cells_with_node_ptr(1:lnode+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_with_node_ptr",1,ierr,error_message)
  !
  ! Loop over the cells, storing 'ahead'
  !
  do n = 1,lcell
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
  do n = 1,lcell
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
  allocate ( cells_surr_cell_ptr(1:lcell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_surr_cell_ptr",1,ierr,error_message)
  !
  ! Loop to set up cells_surr_cell_ptr
  !
  do n = 1,lcell
    cells_surr_cell_ptr(n+1) = cells_surr_cell_ptr(n) + geom_faces(cell_geom(n))
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
  allocate ( lpoin(1:lnode) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"lpoin",1,ierr,error_message)
  !
  do n = 1,lcell
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
  ! if the sum(lpoin(nodes_of_cell(j1:j2))) term never equals nfnods then
  ! this loop exits with ieadj still equal to 0. We cannot send adjncy
  ! containing zeros to the METIS routines so these need to be removed
  ! and the locations in xadj also have to be adjusted accordingly.
  ! Additionally, have xadj to start the indexing for adjncy at 1,
  ! i.e., xadj(1)=1 (cells_surr_cell_ptr indexing for cells_surr_cell
  ! begins at 0, i.e. cells_surr_cell_ptr(1)=0) since this is what
  ! METIS expects when using Fortran array ordering.
  !
 !do n = 1,lcell
 !  n1 = cells_surr_cell_ptr(n)+1
 !  n2 = cells_surr_cell_ptr(n+1)
 !  write (101,'(4(3x,i8))') cells_surr_cell(n1:n2)
 !end do
  !
  allocate ( xadj(1:lcell+1) , source=0 , stat=ierr , errmsg=error_message )
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
  do n = 1,lcell
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
 !do n = 1,lcell
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
  lface = lfbnd
  do n = 1,lcell
    n1 = cells_surr_cell_ptr(n)+1
    n2 = cells_surr_cell_ptr(n+1)
    lface = lface + count( cells_surr_cell(n1:n2) > n )
  end do
  !
  ! Create the face and faces_of_cell arrays
  !
  allocate ( face(1:lface) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face",1,ierr,error_message)
  !
  allocate ( faces_of_cell(1:6,1:lcell) , source=0 , &
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
  do nf = 1,lfbnd
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
    !
    face(nf)%right%cell = ghost_cell
    face(nf)%right%cell_face = 1
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
  face_counter = lfbnd
  do this_cell = 1,lcell
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
      if (adj_cell>this_cell .and. adj_cell<=lcell) then
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
  loop_right_cell_face: do nf = lfbnd+1,lface
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
    right_geom = cell_geom(right_cell)
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
  do nf = 1,lface
    !
    left_cell = face(nf)%left%cell
    right_cell = face(nf)%right%cell
    !
    left_side = face(nf)%left%cell_face
    right_side = face(nf)%right%cell_face
    !
    left_geom = cell_geom(left_cell)
    right_geom = cell_geom(right_cell)
    !
    left_order = cell_order(left_cell)
    right_order = cell_order(right_cell)
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
    if (nf <= lfbnd) then
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
include "Functions/last.f90"
include "Functions/cell_solpts.f90"
!
!###############################################################################
!
end module parallel_mod
