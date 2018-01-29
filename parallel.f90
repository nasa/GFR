module parallel_mod
  !
#include <mpi_defs.h>
#define PROFILE_ON
  !
#ifdef DEBUG_ON
!!!#define DEBUG_PARALLEL
#endif
  !
  !.. Use Statements ..
  use module_kind_types, global_debug_timer => debug_timer
  use eqn_idx, only : nq
  use geovar,  only : face_t
  use geovar,  only : fp_t
  !
#ifdef METIS_5
  use metis5_mod
#else
  use intrinsic :: iso_fortran_env, only : idx_t => int32
  use intrinsic :: iso_fortran_env, only : real_t => real32
#endif
  !
  implicit none
  !
  private
  !
  !.. Public Module Procedures ..
  !
  public :: create_serial_cell_map
  public :: partition_grid
 !public :: collect_global_solpts
  public :: collect_global_solution
 !public :: collect_global_variable
  public :: wait_faceusp
  public :: exchange_faceusp
  public :: wait_facedusp
  public :: exchange_facedusp
  public :: exch_edge_and_node_solutions
  public :: parallel_memory_usage
  public :: exch_connectivity
 !public :: exchange_modes
  !
  integer, save :: nr
  !
  ! Number of communication partners
  !
  integer, save :: ncomm
  !
  ! Request and status arrays for non-blocking exchanges
  !
  _MPI_REQUEST_TYPE_, save, allocatable :: fusp_rqst(:)
  _MPI_REQUEST_TYPE_, save, allocatable :: fdusp_rqst(:)
  !
  _MPI_STATUS_TYPE_, save, allocatable :: fusp_istat(:,:)
  _MPI_STATUS_TYPE_, save, allocatable :: fdusp_istat(:,:)
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
  ! MPI datatype for collecting the global solution on all processors
  !
  _MPI_DATA_TYPE_, save, allocatable, dimension(:) :: recv_sol
  !
  !
  !
  type :: flux_pts_t
    integer, allocatable :: pts(:)
  end type flux_pts_t
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
    logical(lk) :: my_responsibility = fals
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
    integer :: partition = 0
    ! loc_edge : index of this edge local this partition
    integer :: loc_edge = 0
    ! num_cells : number of cells on this partition containing this edge
    integer :: num_cells = 0
    ! edgpts : indices for the edge points on the current edge
    integer, allocatable :: edgpts(:)
  end type cpu_with_edge_t
  !
  logical(lk), parameter :: use_edgpts = true
 !logical(lk), parameter :: use_edgpts = fals
  !
  ! EDGE_T : Derived type to store information for each grid edge
  !          NOTE: This is different from the derived type
  !                of the same name in the module geovar
  !
  type :: edge_t
    ! is_on_boundary : logical that identifies if this edge is on a
    !                  non-communication boundary face
    logical(lk) :: is_on_boundary = fals
    ! controlling_partition : partition that is responsible for the index
    !                         of this node for the CGNS connectivity
    integer :: controlling_partition = 0
    ! order : solution order of the current edge
    integer :: order = 0
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
    integer :: partition = 0
    ! loc_node : index of this node local this partition
    integer :: loc_node = 0
    ! num_cells : number of cells on this partition containing this edge
    integer :: num_cells = 0
  end type cpu_with_node_t
  !
  ! NODE_T : Derived type to store information for each grid node
  !          NOTE: This is different from the derived type
  !                of the same name in the module geovar
  !
  type :: node_t
    ! is_on_boundary : logical that identifies if this edge is on a
    !                  non-communication boundary face
    logical(lk) :: is_on_boundary = fals
    ! controlling_partition : partition that is responsible for the index
    !                         of this node for the CGNS connectivity
    integer :: controlling_partition = 0
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
#ifdef PBS_ENV
    _MPI_DATA_TYPE_ :: send_datatype
#else
    _MPI_DATA_TYPE_ :: send_datatype = MPI_DATATYPE_NULL
#endif
    ! recv_datatype : handle for the MPI user-defined datatype to receive data
#ifdef PBS_ENV
    _MPI_DATA_TYPE_ :: recv_datatype
#else
    _MPI_DATA_TYPE_ :: recv_datatype = MPI_DATATYPE_NULL
#endif
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
 !logical(lk), parameter :: use_mpi_pack = fals
 !logical(lk), parameter :: use_mpi_pack = true
 !!
 !logical(lk), parameter :: send_using_large_buffer = fals
 !logical(lk), parameter :: send_using_large_buffer = true
  !
 !logical(lk), parameter :: send_grid_elem_tags = fals
  logical(lk), parameter :: send_grid_elem_tags = true
  !
  logical(lk), parameter :: send_grid_elem_host_cell = fals
 !logical(lk), parameter :: send_grid_elem_host_cell = true
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
  call debug_timer(entering_procedure,pname)
  !
  allocate ( cell_map(1:1) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_map",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  !
  allocate ( cell_map(1)%loc_to_glb(1:lcell) , source=intseq(1,lcell) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_map(1)%loc_to_glb",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine create_serial_cell_map
!
!###############################################################################
!
subroutine partition_grid(npart,metis_option_requested,bface,nodes_of_cell, &
                          nodes_of_cell_ptr,xyz_nodes,cell_geom,cell_order)
  !
  !.. Use Statements ..
  use ovar, only : continuous_output
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
  integer, allocatable :: cells_with_node(:)
  integer, allocatable :: cells_with_node_ptr(:)
  integer, allocatable :: cells_surr_cell(:)
  integer, allocatable :: cells_surr_cell_ptr(:)
  integer, allocatable :: cells_with_edge(:)
  integer, allocatable :: cells_with_edge_ptr(:)
  integer, allocatable :: nodes_on_edge(:,:)
  integer, allocatable :: face_rotation(:,:)
  !
  integer(idx_t), allocatable :: xadj(:)
  integer(idx_t), allocatable :: adjncy(:)
  integer(idx_t), allocatable :: epart(:)
  !
  !.. Local Allocatable Derived-Type Arrays ..
  type(node_t),     allocatable :: node(:)
  type(map_t),      allocatable :: node_map(:)
  type(edge_t),     allocatable :: edge(:)
  type(map_t),      allocatable :: edge_map(:)
  type(face_t),     allocatable :: face(:)
  type(map_t),      allocatable :: face_map(:)
  type(flux_pts_t), allocatable :: flx(:)
  type(fp_t),       allocatable :: fp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "partition_grid"
 !logical(lk), parameter :: use_old_dual_graph = true
  logical(lk), parameter :: use_old_dual_graph = fals
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Get the size of the global grid
  !
  nr = size(xyz_nodes,dim=1)
  !
  ! Get the number of interior grid cells
  !
  lcell = size(nodes_of_cell_ptr)-1 - size(bface,dim=2)
  !
  if (i_am_host_root) then
    !
    ! Create the global grid connectivity
    !
    call memory_pause("Before: calling get_global_connectivity")
    call get_global_connectivity(cell_geom,cell_order,bface,xyz_nodes, &
                                 nodes_of_cell,nodes_of_cell_ptr, &
                                 cells_with_node,cells_with_node_ptr, &
                                 cells_surr_cell,cells_surr_cell_ptr, &
                                 cells_with_edge,cells_with_edge_ptr, &
                                 nodes_on_edge,face,fp)
    !
    ! Create the dual graph for the global grid not including ghost cells
    !
    call memory_pause("Before: calling create_metis_dual_graph")
    call create_metis_dual_graph(cells_surr_cell,cells_surr_cell_ptr, &
                                 xadj,adjncy)
    !
    ! Deallocate cells_surr_cell and cells_surr_cell_ptr since they
    ! are no longer needed
    !
    if (allocated(cells_surr_cell)) then
      deallocate ( cells_surr_cell , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"cells_surr_cell",2,__LINE__,__FILE__,ierr, &
                       error_message)
    end if
    !
    if (allocated(cells_surr_cell_ptr)) then
      deallocate ( cells_surr_cell_ptr , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"cells_surr_cell_ptr",2,__LINE__,__FILE__,ierr, &
                       error_message)
    end if
    !
    ! Allocate the epart array to store the cell partition numbers
    !
    allocate ( epart(1:lcell) , source=0_idx_t , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"epart",1,__LINE__,__FILE__,ierr,error_message)
   !if (mypnum == 0) then
   !  call output_partitions(99,nodes_of_cell,nodes_of_cell_ptr,xyz_nodes)
   !end if
    !
    ! Use METIS to partition the grid
    ! NOTE: ONLY USE THE GLOBAL ROOT PROCESS TO PARTITION THE GRID
    !
    if (mypnum == glb_root) then
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
      call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
      !
#endif
      !
    end if
    !
    ! Have the global root process broadcast
    ! the partitioning to all host roots
    !
    call mpi_bcast(epart,int(lcell,kind=int_mpi),mpi_inttyp, &
                   glb_root,host_roots_comm,mpierr)
    !
  end if
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  ! Deallocate no longer needed METIS arrays
  !
  if (allocated(xadj)) then
    deallocate ( xadj , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xadj",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(adjncy)) then
    deallocate ( adjncy , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"adjncy",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Distribute partitions across all processors
  !
 !call mpi_bcast(epart,int(lcell,kind=int_mpi), &
 !               mpi_inttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
  !
  ! Output the partition grid to a tecplot file
  !
 !if (mypnum == 0) then
 !  call output_partitions(mypnum,nodes_of_cell,nodes_of_cell_ptr, &
 !                         xyz_nodes,epart)
 !end if
 !call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,"outputing partitions")
  !
  if (i_am_host_root) then
    !
    ! Create the node and edge arrays
    !
    if (continuous_output) then
      call memory_pause("Before: calling create_node_and_edge_arrays")
      call create_node_and_edge_arrays(nodes_of_cell,nodes_of_cell_ptr, &
                                       cells_with_node,cells_with_node_ptr, &
                                       cells_with_edge,cells_with_edge_ptr, &
                                       cell_order,epart,bface,node,edge)
    end if
    !
    ! Deallocate nodes_on_edge since it is no longer needed
    !
    if (allocated(nodes_on_edge)) then
      deallocate ( nodes_on_edge , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"nodes_on_edge",2,__LINE__,__FILE__,ierr, &
                       error_message)
    end if
    !
    ! Find the faces that separate two adjacent partitions
    ! and add these faces to boundary faces array
    !
    call memory_pause("Before: calling find_partition_boundaries")
    call find_partition_boundaries(epart,bface,face,cell_geom,cell_order)
    !
    ! Create the flx and face_rotation arrays from fp and the newly
    ! updated face array that now contains communication boundary faces
    !
    call create_face_rotations(face,fp,flx,face_rotation)
    !
    ! Deallocate fp since it is no longer needed
    !
    if (allocated(fp)) then
      deallocate ( fp , stat=ierr ,errmsg=error_message )
      call alloc_error(pname,"fp",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    ! Now that we have created boundary conditions for faces
    ! located on partition boundaries, use epart to create
    ! mappings from local to global block orderings.
    ! This also includes creating the reverse mapping arrays
    ! or the maps from global to local block orderings
    !
    call memory_pause("Before: calling create_cell_map")
    call create_cell_map(npart,epart,cell_map)
    !
    call memory_pause("Before: calling create_node_map")
    call create_node_map(npart,epart,node,cells_with_node, &
                         cells_with_node_ptr,node_map)
    !
    call memory_pause("Before: calling create_face_map")
    call create_face_map(npart,epart,face,face_map)
    !
    if (allocated(edge)) then
      call memory_pause("Before: calling create_edge_map")
      call create_edge_map(npart,epart,edge,cells_with_edge, &
                           cells_with_edge_ptr,edge_map)
    end if
    !
    ! Deallocate cells_with_node and cells_with_node_ptr
    ! since they are no longer needed
    !
    if (allocated(cells_with_node)) then
      deallocate ( cells_with_node , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"cells_with_node",2,__LINE__,__FILE__,ierr, &
                       error_message)
    end if
    !
    if (allocated(cells_with_node_ptr)) then
      deallocate ( cells_with_node_ptr , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"cells_with_node_ptr",2,__LINE__,__FILE__,ierr, &
                       error_message)
    end if
    !
#ifdef DEBUG_ON
   !call output_mapping_arrays(node_map,edge_map,face_map,cell_map)
#endif
    !
  end if
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  ! Deallocate epart, cells_with_edge, and cells_with_edge_ptr
  ! since they are no longer needed
  !
  if (allocated(epart)) then
    deallocate ( epart , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"epart",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(cells_with_edge)) then
    deallocate ( cells_with_edge , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cells_with_edge",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (allocated(cells_with_edge_ptr)) then
    deallocate ( cells_with_edge_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cells_with_edge_ptr",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  ! Have the root processor broadcast the flx, face_rotation
  ! and the mapping arrays
  !
 !if (use_mpi_pack) then
 !  call pack_fp_and_mapping_arrays(npart,node_map,edge_map,face_map, &
 !                                  cell_map,flx,face_rotation)
 !else
 !  if (send_using_large_buffer) then
 !    call bcast_fp_and_mapping_arrays(npart,node_map,edge_map,face_map, &
 !                                     cell_map,flx,face_rotation)
 !  else
      call bcast_each_fp_and_map_array(npart,node_map,edge_map,face_map, &
                                       cell_map,flx,face_rotation)
 !  end if
 !end if
  !
 !call parallel_memory_usage(iout)
  !
  ! Before we localize everything, we need to create MPI datatypes
  ! to collectively receive localized data from each processor into
  ! a shared global form.
  ! NOTE: The root processor will be the only processor receiving this global
  !       data so it is the only processor that needs to create these datatypes
  !
  if (mypnum == glb_root) then
    call create_global_collective_datatypes(cell_map,cell_geom,cell_order)
  end if
  !
  ! Before we start localizing everything, localize bc_conflict.
  ! Go through each bc conflict and get the partition to which the
  ! boundary face belongs and then localize the boundary face indices.
  ! If we find a face that is on our partition, mark the conflict so
  ! that we know we are involved in the communication for this conflict.
  !
 !call localize_bc_conflict(bface,face_map)
  !
  ! Now that the mappings have been created, we
  ! need to localize the grid on each processor
  !
 !call memory_pause("Before: calling localize_grid")
 !call localize_grid(nodes_of_cell,nodes_of_cell_ptr,bface, &
 !                   xyz_nodes,cell_geom,cell_order, &
 !                   face,node_map,face_map,cell_map)
  call memory_pause("Before: calling root_localize_grid")
  call root_localizes_grid(nodes_of_cell,nodes_of_cell_ptr,bface, &
                           xyz_nodes,cell_geom,cell_order, &
                           face,node_map,face_map,cell_map)
  !
  ! Deallocate face since it is no longer needed
  !
  if (allocated(face)) then
    deallocate ( face , stat=ierr ,errmsg=error_message )
    call alloc_error(pname,"face",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Create MPI datatypes for exchanging node and edge data between processes
  !
  call memory_pause("Before: calling create_node_and_edge_datatypes")
  call bcast_node_and_edge_arrays(node,edge)
 !call report_parallel_memory_usage(xadj,adjncy,epart,node,edge,face, &
 !                                  cells_with_node,cells_with_node_ptr, &
 !                                  cells_surr_cell,cells_surr_cell_ptr, &
 !                                  cells_with_edge,cells_with_edge_ptr, &
 !                                  nodes_on_edge,face_rotation,flx,fp, &
 !                                  node_map,edge_map,face_map)
  !
  if (continuous_output) then
    call create_node_and_edge_datatypes(node,node_map,edge,edge_map)
  end if
  !
  ! Deallocate the node and edge arrays and the node and edge
  ! mapping arrays since they are no longer needed
  !
  if (allocated(node)) then
    deallocate ( node , stat=ierr ,errmsg=error_message )
    call alloc_error(pname,"node",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edge)) then
    deallocate ( edge , stat=ierr ,errmsg=error_message )
    call alloc_error(pname,"edge",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(node_map)) then
    deallocate ( node_map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"node_map",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edge_map)) then
    deallocate ( edge_map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edge_map",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Finally, create MPI datatypes for exchanging face data between processes
  !
  call memory_pause("Before: calling create_face_nonblocking_datatypes")
  call create_face_nonblocking_datatypes(bface,face_map,flx,face_rotation)
 !call create_face_nonblocking_datatypes(bface,face,face_map,flx,face_rotation)
  !
  ! Deallocate face, face_map, flx, and face_rotation
  ! since they are no longer needed
  !
 !if (mypnum == 0) then
 !  write (iout,*)
 !  write (iout,*) "PAUSING BEFORE DEALLOCATING ANY OF THE FACE ARRAYS!"
 !  write (iout,*) "PRESS ANY KEY TO EXECUTE ABORT!"
 !  read (*,*)
 !end if
 !call mpi_barrier(MPI_COMM_WORLD,mpierr)
  if (allocated(face)) then
    deallocate ( face , stat=ierr ,errmsg=error_message )
    call alloc_error(pname,"face",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
 !if (mypnum == 0) then
 !  write (iout,*)
 !  write (iout,*) "PAUSING AFTER DEALLOCATING THE FACE ARRAY!"
 !  write (iout,*) "PRESS ANY KEY TO EXECUTE ABORT!"
 !  read (*,*)
 !end if
 !call mpi_barrier(MPI_COMM_WORLD,mpierr)
  if (allocated(face_map)) then
    deallocate ( face_map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"face_map",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
 !if (mypnum == 0) then
 !  write (iout,*)
 !  write (iout,*) "PAUSING AFTER DEALLOCATING THE FACE_MAP ARRAY!"
 !  write (iout,*) "PRESS ANY KEY TO EXECUTE ABORT!"
 !  read (*,*)
 !end if
 !call mpi_barrier(MPI_COMM_WORLD,mpierr)
  if (allocated(flx)) then
    deallocate ( flx , stat=ierr ,errmsg=error_message )
    call alloc_error(pname,"flx",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
 !if (mypnum == 0) then
 !  write (iout,*)
 !  write (iout,*) "PAUSING AFTER DEALLOCATING THE FLX ARRAY!"
 !  write (iout,*) "PRESS ANY KEY TO EXECUTE ABORT!"
 !  read (*,*)
 !end if
 !call mpi_barrier(MPI_COMM_WORLD,mpierr)
  if (allocated(face_rotation)) then
    deallocate ( face_rotation , stat=ierr ,errmsg=error_message )
    call alloc_error(pname,"face_rotation",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
 !if (mypnum == 0) then
 !  write (iout,*)
 !  write (iout,*) "PAUSING AFTER DEALLOCATING THE FACE_ROTATION ARRAY!"
 !  write (iout,*) "PRESS ANY KEY TO EXECUTE ABORT!"
 !  read (*,*)
 !end if
  !
  call debug_timer(leaving_procedure,pname)
  !
 !if (mypnum == 0) then
 !  write (iout,*)
 !  write (iout,*) "PAUSING BEFORE ABORTING! "
 !  write (iout,*) "PRESS ANY KEY TO EXECUTE ABORT!"
 !  read (*,*)
 !end if
 !call mpi_barrier(MPI_COMM_WORLD,mpierr)
 !call mpi_barrier(MPI_COMM_WORLD,mpierr)
 !call stop_gfr(stop_mpi,pname,__LINE__,__FILE__, &
 !                "TEMPORARY STOP AFTER PARTITIONING THE GRID!!")
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
  integer,           intent(in) :: metis_option_requested
  integer,           intent(in) :: npart
  integer,           intent(in) :: cell_geom(:)
  integer,           intent(in) :: cell_order(:)
  integer(idx_t),    intent(in) :: xadj(:)
  integer(idx_t),    intent(in) :: adjncy(:)
  integer(idx_t), intent(inout) :: epart(:)
  !
  !.. Local Scalars ..
  integer :: n,ip,ierr,lcell
  integer :: wgtflag,numflag,vol
  integer :: metis_option
  logical(ldk) :: passes_inspection
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
  call debug_timer(entering_procedure,pname)
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
        call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
        !
      end if
      !
      ! I dont trust that the METIS routines will compute identical
      ! partitioning on each processor so have the root process check
      ! its local partitioning for errors and then broadcast that
      ! partitioning to the other processes if none are found.
      !
     !if (mypnum == 0) then
        call check_quality_of_partitioning(npart,epart, &
                                           passes_inspection)
     !end if
      !
      ! Broadcast the logical passes_inspection to all processors so
      ! that they all will know whether to exit part_loop.
      !
     !call mpi_bcast(passes_inspection,1_int_mpi,MPI_LOGICAL, &
     !               0_int_mpi,MPI_COMM_WORLD,mpierr)
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
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
    !
  end do part_loop
  !
  ! Deallocate the weight arrays before we leave
  !
  if (allocated(vwgt)) then
    deallocate ( vwgt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vwgt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(adjwgt)) then
    deallocate ( adjwgt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"adjwgt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(vsize)) then
    deallocate ( vsize , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vsize",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
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
  !.. Formal Arguments ..
  integer,           intent(in) :: metis_option_requested
  integer,           intent(in) :: npart
  integer(idx_t),    intent(in) :: xadj(:)
  integer(idx_t),    intent(in) :: adjncy(:)
  integer(idx_t), intent(inout) :: epart(:)
  !
  !.. Local Scalars ..
  integer        :: n,ip,ierr,lcell,metis_option
  integer(idx_t) :: nvtxs,ncon,nparts,objval
  logical(ldk)   :: passes_inspection
  character(len=19) :: method
  character(len=12) :: weight
  character(len=20) :: minmzd
  !
  !.. Local Arrays ..
  integer :: metis_option_hierarchy(1:5)
#ifdef METIS_5
  integer(idx_t) :: options(0:METIS_NOPTIONS-1)
#endif
  !
  !.. Local Allocatable Arrays ..
  integer(idx_t), allocatable :: vwgt(:)
  integer(idx_t), allocatable :: vsize(:)
  integer(idx_t), allocatable :: adjwgt(:)
  real(real_t),   allocatable :: tpwgts(:)
  real(real_t),   allocatable :: ubvec(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "partition_using_metis_5"
  !
continue
#ifdef METIS_5
  !
  call debug_timer(entering_procedure,pname)
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
  nvtxs = int(lcell,kind=idx_t)
  ncon = 1_idx_t ! number of weights associated with each vertex (>= 1)
  nparts = int(npart,kind=idx_t)
  !
  ! Allocate the local arrays
  !
#ifdef DEBUG_ON
  allocate ( vwgt(1:ncon*lcell) , source=1_idx_t , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"vwgt",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( vsize(1:lcell) , source=1_idx_t , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"vsize",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( adjwgt(1:size(adjncy)) , source=1_idx_t , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"adjwgt",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( tpwgts(0:ncon*npart-1) , &
             source=1.0_real_t/real(npart,kind=real_t) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tpwgts",1,__LINE__,__FILE__,ierr,error_message)
  !
  tpwgts(:) = real( one/real(npart,kind=wp) , kind=kind(tpwgts) )
  tpwgts(:) = tpwgts(:) / sum(tpwgts(:))
  !
  allocate ( ubvec(1:ncon) , source=1.01_real_t , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ubvec",1,__LINE__,__FILE__,ierr,error_message)
  !
  if (mypnum == glb_root) write (iout,5) mypnum,sum(tpwgts)
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
        call metis_partgraphrecursive(nvtxs,ncon,xadj,adjncy,vwgt,vsize, &
                                      adjwgt,nparts,tpwgts,ubvec,options, &
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
        call metis_partgraphkway(nvtxs,ncon,xadj,adjncy,vwgt,vsize, &
                                 adjwgt,nparts,tpwgts,ubvec,options, &
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
        call metis_partgraphkway(nvtxs,ncon,xadj,adjncy,vwgt,vsize, &
                                 adjwgt,nparts,tpwgts,ubvec,options, &
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
        call metis_partgraphkway(nvtxs,ncon,xadj,adjncy,vwgt,vsize, &
                                 adjwgt,nparts,tpwgts,ubvec,options, &
                                 objval,epart)
        !
      else
        !
        ! The METIS option requested is invalid!
        !
        write (error_message,3)
        call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
        !
      end if
      !
      ! I dont trust that the METIS routines will compute identical
      ! partitioning on each processor so have the root process check
      ! its local partitioning for errors and then broadcast that
      ! partitioning to the other processes if none are found.
      !
     !if (mypnum == 0) then
        call check_quality_of_partitioning(npart,epart, &
                                           passes_inspection)
     !end if
      !
      ! Broadcast the logical passes_inspection to all processors so
      ! that they all will know whether to exit part_loop.
      !
     !call mpi_bcast(passes_inspection,1_int_mpi,MPI_LOGICAL, &
     !               0_int_mpi,MPI_COMM_WORLD,mpierr)
      !
      ! If the partitioning passes inspection, have the root processor
      ! output the results of the partitioning and exit from part_loop.
      !
      if (passes_inspection) then
        write (iout,1) method,npart,weight,minmzd,objval
        exit part_loop
      end if
      !
      write (iout,4) method,npart,weight,minmzd
      !
    end do
    !
    write (error_message,2)
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
    !
  end do part_loop
  !
  ! Deallocate all the local arrays before leaving
  !
  if (allocated(vwgt)) then
    deallocate ( vwgt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vwgt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(vsize)) then
    deallocate ( vsize , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vsize",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(adjwgt)) then
    deallocate ( adjwgt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"adjwgt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(tpwgts)) then
    deallocate ( tpwgts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"tpwgts",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(ubvec)) then
    deallocate ( ubvec , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ubvec",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
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
  5 format (" CPU # ",i0," : sum(tpwgts) = ",es14.6)
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
                                   nodes_on_edge,face,fp)
  !
  !.. Use Statements ..
  use connectivity_mod, only : find_cells_surrounding_each_node
  use connectivity_mod, only : find_cells_surrounding_each_cell
  use connectivity_mod, only : find_nodes_on_each_edge
  use connectivity_mod, only : find_cells_surrounding_each_edge
  use connectivity_mod, only : create_face_derived_type
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: cell_geom(:)
  integer,                 intent(in) :: cell_order(:)
  integer,                 intent(in) :: bface(:,:)
  integer,                 intent(in) :: nodes_of_cell(:)
  integer,                 intent(in) :: nodes_of_cell_ptr(:)
  real(wp),                intent(in) :: xyz_nodes(:,:)
  integer, allocatable, intent(inout) :: cells_with_node(:)
  integer, allocatable, intent(inout) :: cells_with_node_ptr(:)
  integer, allocatable, intent(inout) :: cells_surr_cell(:)
  integer, allocatable, intent(inout) :: cells_surr_cell_ptr(:)
  integer, allocatable, intent(inout) :: cells_with_edge(:)
  integer, allocatable, intent(inout) :: cells_with_edge_ptr(:)
  integer, allocatable, intent(inout) :: nodes_on_edge(:,:)
  !
  type(face_t), allocatable, intent(inout) :: face(:)
  type(fp_t),   allocatable, intent(inout) :: fp
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Allocatable Arrays ..
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_global_connectivity"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
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
                                face,fp,force_flxpt_offset_to_zero=true)
  !
  call debug_timer(leaving_procedure,pname)
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
  integer,                        intent(in) :: cells_surr_cell(:)
  integer,                        intent(in) :: cells_surr_cell_ptr(:)
  integer(idx_t), allocatable, intent(inout) :: xadj(:)
  integer(idx_t), allocatable, intent(inout) :: adjncy(:)
  !
  !.. Local Scalars ..
  integer :: n,nf,n1,n2,ic,ierr,lcell
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_metis_dual_graph"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
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
  allocate ( xadj(1:lcell+1) , source=0_idx_t , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xadj",1,__LINE__,__FILE__,ierr,error_message)
  xadj(1) = 1_idx_t
  !
  n = count( cells_surr_cell(:) > 0 )
  allocate ( adjncy(1:n) , source=0_idx_t , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"adjncy",1,__LINE__,__FILE__,ierr,error_message)
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
        xadj(n+1) = xadj(n+1) + 1_idx_t
        ic = ic + 1
        adjncy(ic) = int( cells_surr_cell(nf) , kind=idx_t )
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
  call debug_timer(leaving_procedure,pname)
  !
end subroutine create_metis_dual_graph
!
!###############################################################################
!
subroutine create_grid_weights(metis_option,ncon,n_dual_vertices,n_dual_edges, &
                               vwgt,adjwgt,vsize,wgtflag,cell_geom,cell_order)
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: metis_option
  integer,                 intent(in) :: ncon
  integer,                 intent(in) :: n_dual_vertices
  integer,                 intent(in) :: n_dual_edges
  integer,                 intent(in) :: cell_geom(:)
  integer,                 intent(in) :: cell_order(:)
  integer,              intent(inout) :: wgtflag
  integer, allocatable, intent(inout) :: vwgt(:)
  integer, allocatable, intent(inout) :: adjwgt(:)
  integer, allocatable, intent(inout) :: vsize(:)
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_grid_weights"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
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
    call alloc_error(pname,"vwgt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( vwgt(1:ncon*n_dual_vertices) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"vwgt",1,__LINE__,__FILE__,ierr,error_message)
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
      call alloc_error(pname,"adjwgt",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    allocate ( adjwgt(1:n_dual_edges) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"adjwgt",1,__LINE__,__FILE__,ierr,error_message)
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
      call alloc_error(pname,"vsize",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    allocate ( vsize(1:n_dual_vertices) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vsize",1,__LINE__,__FILE__,ierr,error_message)
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
  call debug_timer(leaving_procedure,pname)
  !
end subroutine create_grid_weights
!
!###############################################################################
!
subroutine check_quality_of_partitioning(npart,epart,passes_inspection)
  !
  !.. Formal Arguments ..
  integer,        intent(in) :: npart
  integer(idx_t), intent(in) :: epart(:)
  logical(ldk),  intent(out) :: passes_inspection
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_quality_of_partitioning"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  passes_inspection = true
  !
  ! Check that all processors got part of the global grid
  !
  do n = 1,npart
    if (any(epart == n)) cycle
    write (iout,100) n
    passes_inspection = fals
   !call stop_gfr(abort,pname,__LINE__,__FILE__)
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
  call debug_timer(leaving_procedure,pname)
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
  logical(lk), optional, intent(in) :: extra_line
  logical(lk), optional, intent(in) :: surpress_cpu
  !
  integer :: n,io
  logical(lk) :: add_extra_line
  logical(lk) :: include_cpu
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
                                       cell_order,epart,bface,node,edge)
  !
  !.. Formal Arguments ..
  integer,                      intent(in) :: nodes_of_cell(:)
  integer,                      intent(in) :: nodes_of_cell_ptr(:)
  integer,                      intent(in) :: cells_with_node(:)
  integer,                      intent(in) :: cells_with_node_ptr(:)
  integer,                      intent(in) :: cells_with_edge(:)
  integer,                      intent(in) :: cells_with_edge_ptr(:)
  integer(idx_t),               intent(in) :: epart(:)
  integer,                      intent(in) :: cell_order(:)
  integer,                      intent(in) :: bface(:,:)
  type(node_t), allocatable, intent(inout) :: node(:)
  type(edge_t), allocatable, intent(inout) :: edge(:)
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,nf,np,i,ierr
  integer :: lnode,ledge,lcell
  integer :: lcpu,ucpu,num_cpus
  integer :: host_cell,ghost_cell
  character(len=100) :: array_name
  !
  !.. Local Arrays ..
  logical(lk), dimension(1:ncpu) :: cpu_mask
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
  call debug_timer(entering_procedure,pname)
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
  n2 = 0
  if (ledge > 0) then
    n2 = maxval( cells_with_edge_ptr(2:ledge+1) - &
                 cells_with_edge_ptr(1:ledge  ) )
  end if
  !
  !
  ! Allocate the array to store the partitions numbers for all cells with a
  ! given node/edge, and the array to store the BCs those boundary faces
  ! with a given node/edge.
  !
  n = max(n1,n2)
  !
  allocate ( cell_list(1:n) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_list",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( part_list(1:n) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"part_list",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( bc_list(1:n) , source=not_a_bc , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bc_list",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the node array
  !
  allocate ( node(1:lnode) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"node",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Loop through the grid nodes and find any nodes that exist on multiple
  ! partitions. Save these partition numbers to the node array
  !
  create_node_loop: do n = 1,lnode
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
            call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
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
    if (lcpu == ucpu) then
      node(n)%controlling_partition = lcpu
      cycle create_node_loop
    end if
    !
    ! Sum the number of cells containing this node on each partition
    !
    do i = 1,np
      !
      ! Mark this partition as having a cell containing this node
      ! NOTE: For a boundary node, we want to keep all partitions
      !       that contain this node, not just the ones with
      !       bc_list /= bc_unknown, because all partitions will
      !       still have to provide the connectivity for this node.
      !       To make sure that those partitions with
      !       bc_list == bc_unknown dont contribute to the solution
      !       at this node, the num_cells component will be 0 for
      !       these partitions.
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
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
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
    i = minval( cell_list(1:np) , &
                mask=(cell_list(1:np)>0 .and. cell_list(1:np)<=lcell) )
    !
    node(n)%controlling_partition = epart(i)
    !
  end do create_node_loop
  !
  ! Create the edge array
  !
  allocate ( edge(1:ledge) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"edge",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Loop through the grid edges and find any edges that exist on multiple
  ! partitions. Save these partition numbers to the edge array
  !
  create_edge_loop: do n = 1,ledge
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
            call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
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
    if (lcpu == ucpu) then
      edge(n)%controlling_partition = lcpu
      cycle create_edge_loop
    end if
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
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
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
    i = minval( cell_list(1:np) , &
                mask=(cell_list(1:np)>0 .and. cell_list(1:np)<=lcell) )
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
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
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
  end do create_edge_loop
  !
  ! Deallocate cell_list, part_list, and bc_list before we leave
  !
  if (allocated(cell_list)) then
    deallocate ( cell_list , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_list",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(part_list)) then
    deallocate ( part_list , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"part_list",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(bc_list)) then
    deallocate ( bc_list , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bc_list",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
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
subroutine find_partition_boundaries(epart,bface,face,cell_geom,cell_order)
  !
  !.. Formal Arguments ..
  integer(idx_t),               intent(in) :: epart(:)
  integer,      allocatable, intent(inout) :: bface(:,:)
  type(face_t), allocatable, intent(inout) :: face(:)
  integer,      allocatable, intent(inout) :: cell_geom(:)
  integer,      allocatable, intent(inout) :: cell_order(:)
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
  logical(lk),  allocatable :: msk(:)
  type(face_t), allocatable :: temp_face(:)
  integer,      allocatable :: temp_geom(:)
  integer,      allocatable :: temp_order(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_partition_boundaries"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
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
  call alloc_error(pname,"msk",1,__LINE__,__FILE__,ierr,error_message)
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
  ! Replace bface(10,nf) with the location within the face array to extract data
  ! and set the responsible partition to the partition for the host cell
  !
  do nf = 1,lfbnd
    bface(10,nf) = nf
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
    !                          ! NOTE: This is positive to indicate that this
    nf = nf + 1                !       partition is responsible for this new
    !                          !       boundary face
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
    !                          ! NOTE: This is negative to indicate that this
    nf = nf + 1                !       host cell used to be the right cell on
    !                          !       the original interior face and that this
  end do                       !       partition is not responsible for this
  !                            !       new boundary face
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Allocate a temporary array to resize and reorder face
  !
  allocate ( temp_face(1:new_lface) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_face",1,__LINE__,__FILE__,ierr,error_message)
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
  call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! #######
  ! ## 4 ##
  ! #######
  !
  ! Update the boundary face sections of cell_geom and cell_order
  ! to include the newly created boundary faces
  !
  call reallocate(cell_geom,lcell+new_lfbnd)
  call reallocate(cell_order,lcell+new_lfbnd)
  !
  do nf = 1,new_lfbnd
    cell_geom(lcell+nf) = face(nf)%geom
    cell_order(lcell+nf) = face(nf)%order
  end do
  !
  ! #######
  ! ## 5 ##
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
     !bface(7,nf) = mypnum
      bface(7,nf) = epart(bface(2,nf))-1
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
    call stop_gfr(abort,pname,__LINE__,__FILE__)
  end if
  !
  call debug_timer(leaving_procedure,pname)
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
subroutine create_face_rotations(face,fp,flx,face_rotation)
  !
  !.. Formal Arguments ..
  type(face_t),                     intent(in) :: face(:)
  type(fp_t),                       intent(in) :: fp
  type(flux_pts_t), allocatable, intent(inout) :: flx(:)
  integer,          allocatable, intent(inout) :: face_rotation(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,n,nf,np,ierr
  integer :: ibeg,jbeg,kbeg
  integer :: iend,jend,kend
  integer :: face_geom,right_rotation
  integer :: face_order,left_rotation
  integer :: face_pts
  !
  character(len=100) :: array_name
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: fpidx(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_face_rotation"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ibeg = lbound(fp%geord,dim=1)
  iend = ubound(fp%geord,dim=1)
  !
  jbeg = lbound(fp%geord,dim=2)
  jend = ubound(fp%geord,dim=2)
  !
  kbeg = huge(0)
  kend = -huge(0)
  !
  do j = jbeg,jend
    do i = ibeg,iend
      kbeg = min(kbeg,lbound(fp%geord(i,j)%rot,dim=1))
      kend = max(kend,ubound(fp%geord(i,j)%rot,dim=1))
    end do
  end do
  !
  allocate ( fpidx(ibeg:iend,jbeg:jend,kbeg:kend) , source=-huge(0) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fpidx",1,__LINE__,__FILE__,ierr,error_message)
  !
  n = 0
  do j = jbeg,jend
    do i = ibeg,iend
      do k = lbound(fp%geord(i,j)%rot,dim=1),ubound(fp%geord(i,j)%rot,dim=1)
        if (any(fp%geord(i,j)%rot(k)%pts(:) /= 0)) then
          n = n + 1
          fpidx(i,j,k) = n
        end if
      end do
    end do
  end do
  !
  allocate ( flx(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"flx",1,__LINE__,__FILE__,ierr,error_message)
  !
  n = 0
  do j = jbeg,jend
    do i = ibeg,iend
      do k = lbound(fp%geord(i,j)%rot,dim=1),ubound(fp%geord(i,j)%rot,dim=1)
        if (any(fp%geord(i,j)%rot(k)%pts(:) /= 0)) then
          !
          n = n + 1
          !
          np = size(fp%geord(i,j)%rot(k)%pts)
          !
          allocate ( flx(n)%pts(1:np) , &
                     source=fp%geord(i,j)%rot(k)%pts(1:np) , &
                     stat=ierr , errmsg=error_message )
          write (array_name,1) n
          call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                           error_message)
          !
        end if
      end do
    end do
  end do
  !
  allocate ( face_rotation(1:2,1:size(face)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_rotation",1,__LINE__,__FILE__,ierr,error_message)
  !
  do nf = 1,size(face)
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    left_rotation = face(nf)%left%rotation
    right_rotation = face(nf)%right%rotation
    !
    face_rotation(1,nf) = fpidx(face_geom,face_order,left_rotation)
    face_rotation(2,nf) = fpidx(face_geom,face_order,right_rotation)
    !
  end do
  !
  deallocate ( fpidx , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fpidx",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("flx(",i0,")%pts")
  2 format (a,i0)
  !
end subroutine create_face_rotations
!
!###############################################################################
!
subroutine output_node_arrays(node,node_map,bnd_nodes)
  !
  type(node_t),     allocatable, optional, intent(in) :: node(:)
  type(map_t),      allocatable, optional, intent(in) :: node_map(:)
  type(bnd_data_t), allocatable, optional, intent(in) :: bnd_nodes(:)
  !
  integer :: i,n,io
  !
continue
  !
  if (present(node)) then
    if (allocated(node)) then
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
    else
      write (io,17)
    end if
  end if
  !
  10 format (" Contents of node(",i0,")")
  11 format (5x,"%is_on_boundary        = ",l1)
  12 format (5x,"%controlling_partition = ",i0)
  13 format (5x,"%cpu(",i0,")")
  14 format (10x,"%partition = ",i0)
  15 format (10x,"%loc_node  = ",i0)
  16 format (10x,"%num_cells = ",i0)
  17 format (" THE node ARRAY IS NOT ALLOCATED!!!")
  !
  if (present(node_map)) then
    if (allocated(node_map)) then
      io = mypnum+150
      do n = 1,size(node_map)
        if (allocated(node_map(n)%loc_to_glb)) then
          write (io,20) n
          do i = 1,size(node_map(n)%loc_to_glb)
            write (io,21) i,node_map(n)%loc_to_glb(i)
          end do
        else
          write (io,20) n, "IS NOT ALLOCATED!!!"
        end if
        write (io,22)
      end do
    else
      write (io,23)
    end if
  end if
  !
  20 format (" Map for CPU ",i0,:,1x,a)
  21 format ("   loc = ",i5," => glb = ",i5)
  22 format (30("#"))
  23 format (" THE node_map ARRAY IS NOT ALLOCATED!!!")
  !
  if (present(bnd_nodes)) then
    if (allocated(bnd_nodes)) then
      io = mypnum+160
      do n = 1,size(bnd_nodes)
        write (io,30) n
        write (io,31) bnd_nodes(n)%idx
        write (io,32) bnd_nodes(n)%ave_correction
        write (io,33) bnd_nodes(n)%my_responsibility
      end do
    else
      write (io,34)
    end if
  end if
  !
  30 format (" Contents of bnd_nodes(",i0,")")
  31 format (5x,"%idx               = ",i0)
  32 format (5x,"%ave_correction    = ",f9.6)
  33 format (5x,"%my_responsibility = ",l1)
  34 format (" THE bnd_nodes ARRAY IS NOT ALLOCATED!!!")
  !
end subroutine output_node_arrays
!
!###############################################################################
!
subroutine create_cell_map(npart,epart,cell_map)
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: npart
  integer(idx_t),              intent(in) :: epart(:)
  type(map_t), allocatable, intent(inout) :: cell_map(:)
  !
  !.. Local Scalars ..
  integer :: i,j,nc,ierr,this_part
  integer :: gcell,cell_div,init_cell
  character(len=100) :: array_name
  !
  !.. Local Arrays ..
  integer :: cell_count(1:npart)
  !
  !.. Local Allocatable Arrays ..
  type(map_t), allocatable :: tmp_map(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_cell_map"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  gcell = size(epart)
  !
  cell_div = gcell / npart
  !
  init_cell = min( gcell , cell_div*2 )
  !
  ! Create the cell mappings
  !
  allocate ( tmp_map(1:npart) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_map",1,__LINE__,__FILE__,ierr,error_message)
  !
  do this_part = 1,size(tmp_map)
    allocate ( tmp_map(this_part)%loc_to_glb(1:init_cell) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) this_part
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  end do
  !
  cell_count(:) = 0
  do nc = 1,gcell
    this_part = epart(nc)
    i = cell_count(this_part) + 1
    if (i > size(tmp_map(this_part)%loc_to_glb)) then
      tmp_map(this_part)%loc_to_glb = [tmp_map(this_part)%loc_to_glb(:), &
                                       (0,j=1,cell_div)]
    end if
    tmp_map(this_part)%loc_to_glb(i) = nc
    cell_count(this_part) = i
  end do
  !
  ! Copy tmp_map over to cell_map
  !
  allocate ( cell_map(1:npart) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_map",1,__LINE__,__FILE__,ierr,error_message)
  !
  do this_part = 1,size(cell_map)
    nc = cell_count(this_part)
    ! F2003 AUTO-REALLOCATION
    cell_map(this_part)%loc_to_glb = tmp_map(this_part)%loc_to_glb(1:nc)
  end do
  !
  deallocate ( tmp_map , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_map",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("tmp_map(",i0,")%loc_to_glb")
  !
end subroutine create_cell_map
!
!###############################################################################
!
subroutine create_node_map(npart,epart,node,cells_with_node, &
                           cells_with_node_ptr,node_map)
  !
  !.. Formal Arguments ..
  integer,                      intent(in) :: npart
  integer(idx_t),               intent(in) :: epart(:)
  integer,                      intent(in) :: cells_with_node(:)
  integer,                      intent(in) :: cells_with_node_ptr(:)
  type(node_t), allocatable, intent(inout) :: node(:)
  type(map_t),  allocatable, intent(inout) :: node_map(:)
  !
  !.. Local Scalars ..
  integer :: i,j,n,nn,n1,n2,ierr
  integer :: this_part,loc_node,glb_node
  integer :: gnode,node_div,init_node
  character(len=200) :: array_name
  !
  !.. Local Arrays ..
  integer     :: node_count(1:npart)
  logical(lk) :: part_mask(1:npart)
  !
  !.. Local Allocatable Arrays ..
  type(map_t), allocatable :: tmp_map(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_node_map"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  gnode = size(cells_with_node_ptr)-1
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  node_div = gnode / npart
  !
  init_node = min( gnode , node_div*4 )
  !
  allocate ( tmp_map(1:npart) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_map",1,__LINE__,__FILE__,ierr,error_message)
  !
  do this_part = 1,size(tmp_map)
    allocate ( tmp_map(this_part)%loc_to_glb(1:init_node) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) this_part
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  end do
  !
  node_count(:) = 0
  do nn = 1,gnode
    part_mask(:) = fals
    n1 = cells_with_node_ptr(nn)+1
    n2 = cells_with_node_ptr(nn+1)
    do n = n1,n2
      if (cells_with_node(n) > 0) then
        this_part = epart(cells_with_node(n))
        if (part_mask(this_part)) cycle
        i = node_count(this_part) + 1
        if (i > size(tmp_map(this_part)%loc_to_glb)) then
          tmp_map(this_part)%loc_to_glb = [tmp_map(this_part)%loc_to_glb(:), &
                                           (0,j=1,node_div)]
        end if
        tmp_map(this_part)%loc_to_glb(i) = nn
        node_count(this_part) = i
        part_mask(this_part) = true
      end if
    end do
  end do
  !
  ! Create the node mappings
  !
  allocate ( node_map(1:npart) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"node_map",1,__LINE__,__FILE__,ierr,error_message)
  !
  do this_part = 1,size(node_map)
    !
    nn = node_count(this_part)
    !
    ! Copy tmp_map over to node_map
    !
    ! F2003 AUTO-REALLOCATION
    node_map(this_part)%loc_to_glb = tmp_map(this_part)%loc_to_glb(1:nn)
    !
  end do
  !
  deallocate ( tmp_map , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_map",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Save the local node index for each node of this partition that
  ! also exists on a different partition
  ! NOTE: The cpu component of the node array should only be allocated
  !       if the global node exists on multiple partitions
  !
  if (allocated(node)) then
    !
    do this_part = 1,size(node_map)
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
  end if
  !
 !call output_node_arrays(node,node_map)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("tmp_map(",i0,")%loc_to_glb")
  !
end subroutine create_node_map
!
!###############################################################################
!
subroutine create_face_map(npart,epart,face,face_map)
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: npart
  integer(idx_t),              intent(in) :: epart(:)
  type(face_t),                intent(in) :: face(:)
  type(map_t), allocatable, intent(inout) :: face_map(:)
  !
  !.. Local Scalars ..
  integer :: i,j,nf,ierr,this_part
  integer :: gface,face_div,init_face
  character(len=100) :: array_name
  !
  !.. Local Arrays ..
  integer :: face_count(1:npart)
  !
  !.. Local Allocatable Arrays ..
  type(map_t), allocatable :: tmp_map(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_face_map"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Create the face mappings
  !
  gface = size(face)
  !
  face_div = gface / npart
  !
  init_face = min( gface , face_div*4 )
  !
  allocate ( tmp_map(1:npart) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_map",1,__LINE__,__FILE__,ierr,error_message)
  !
  do this_part = 1,size(tmp_map)
    allocate ( tmp_map(this_part)%loc_to_glb(1:init_face) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) this_part
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  end do
  !
  face_count(:) = 0
  do nf = 1,gface
    this_part = epart( face(nf)%left%cell )
    i = face_count(this_part) + 1
    if (i > size(tmp_map(this_part)%loc_to_glb)) then
      tmp_map(this_part)%loc_to_glb = [tmp_map(this_part)%loc_to_glb(:), &
                                       (0,j=1,face_div)]
    end if
    tmp_map(this_part)%loc_to_glb(i) = nf
    face_count(this_part) = i
  end do
  !
  ! Copy tmp_map over to face_map
  !
  allocate ( face_map(1:npart) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_map",1,__LINE__,__FILE__,ierr,error_message)
  !
  do this_part = 1,size(face_map)
    nf = face_count(this_part)
    ! F2003 AUTO-REALLOCATION
    face_map(this_part)%loc_to_glb = tmp_map(this_part)%loc_to_glb(1:nf)
  end do
  !
  deallocate ( tmp_map , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_map",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("tmp_map(",i0,")%loc_to_glb")
  !
end subroutine create_face_map
!
!###############################################################################
!
subroutine create_edge_map(npart,epart,edge,cells_with_edge, &
                           cells_with_edge_ptr,edge_map)
  !
  !.. Formal Arguments ..
  integer,                      intent(in) :: npart
  integer(idx_t),               intent(in) :: epart(:)
  integer,                      intent(in) :: cells_with_edge(:)
  integer,                      intent(in) :: cells_with_edge_ptr(:)
  type(edge_t), allocatable, intent(inout) :: edge(:)
  type(map_t),  allocatable, intent(inout) :: edge_map(:)
  !
  !.. Local Scalars ..
  integer :: i,j,n,ne,n1,n2,ierr
  integer :: this_part,loc_edge,glb_edge
  integer :: gedge,edge_div,init_edge
  character(len=200) :: array_name
  !
  !.. Local Arrays ..
  integer     :: edge_count(1:npart)
  logical(lk) :: part_mask(1:npart)
  !
  !.. Local Allocatable Arrays ..
  type(map_t), allocatable :: tmp_map(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_edge_map"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  gedge = size(cells_with_edge_ptr)-1
  !
 !if (mypnum == glb_root) write (iout,15) gedge
 !15 format (" size(edge) = ",i0)
  !
  edge_div = gedge / npart
  !
  init_edge = min( gedge , edge_div*4 )
  !
  allocate ( tmp_map(1:npart) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_map",1,__LINE__,__FILE__,ierr,error_message)
  !
  do this_part = 1,size(tmp_map)
    allocate ( tmp_map(this_part)%loc_to_glb(1:init_edge) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) this_part
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  end do
  !
  edge_count(:) = 0
  do ne = 1,gedge
    part_mask(:) = fals
    n1 = cells_with_edge_ptr(ne)+1
    n2 = cells_with_edge_ptr(ne+1)
    do n = n1,n2
      if (cells_with_edge(n) > 0) then
        this_part = epart(cells_with_edge(n))
        if (part_mask(this_part)) cycle
        i = edge_count(this_part) + 1
        if (i > size(tmp_map(this_part)%loc_to_glb)) then
          tmp_map(this_part)%loc_to_glb = [tmp_map(this_part)%loc_to_glb(:), &
                                           (0,j=1,edge_div)]
        end if
        tmp_map(this_part)%loc_to_glb(i) = ne
        edge_count(this_part) = i
        part_mask(this_part) = true
      end if
    end do
  end do
  !
  ! Create the edge mappings
  !
  allocate ( edge_map(1:npart) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"edge_map",1,__LINE__,__FILE__,ierr,error_message)
  !
  do this_part = 1,size(edge_map)
    !
    ne = edge_count(this_part)
    !
    ! Copy tmp_map over to edge_map
    !
    ! F2003 AUTO-REALLOCATION
    edge_map(this_part)%loc_to_glb = tmp_map(this_part)%loc_to_glb(1:ne)
    !
  end do
  !
  deallocate ( tmp_map , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_map",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Save the local edge index for each edge of this partition that
  ! also exists on a different partition
  ! NOTE: The cpu component of the edge array should only be allocated
  !       if the global edge exists on multiple partitions
  !
  if (allocated(edge)) then
    !
    do this_part = 1,size(edge_map)
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
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("tmp_map(",i0,")%loc_to_glb")
  !
end subroutine create_edge_map
!
!###############################################################################
!
subroutine output_mapping_arrays(node_map,edge_map,face_map,cell_map)
  !
  !.. Formal Arguments ..
  type(map_t), allocatable, intent(in) :: node_map(:)
  type(map_t), allocatable, intent(in) :: edge_map(:)
  type(map_t), allocatable, intent(in) :: face_map(:)
  type(map_t), allocatable, intent(in) :: cell_map(:)
  !
  !.. Local Scalars ..
  integer :: ierr,io,l,this_part
  character(len=50) :: fname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_mapping_arrays"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Output node_map to file
  !
  write (fname,1) "node",mypnum
  !
  open (newunit=io,file=fname,status="replace",action="write", &
        iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  if (allocated(node_map)) then
    write (io,11) "node",size(node_map)
    do this_part = 1,size(node_map)
      if (allocated(node_map(this_part)%loc_to_glb)) then
        write (io,12) "node",this_part,size(node_map(this_part)%loc_to_glb)
        do l = 1,size(node_map(this_part)%loc_to_glb)
          write (io,13) "node",this_part,l,node_map(this_part)%loc_to_glb(l)
        end do
      else
        write (io,22) "node",this_part
      end if
    end do
  else
    write (io,21) "node"
  end if
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Output edge_map to file
  !
  write (fname,1) "edge",mypnum
  !
  open (newunit=io,file=fname,status="replace",action="write", &
        iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  if (allocated(edge_map)) then
    write (io,11) "edge",size(edge_map)
    do this_part = 1,size(edge_map)
      if (allocated(edge_map(this_part)%loc_to_glb)) then
        write (io,12) "edge",this_part,size(edge_map(this_part)%loc_to_glb)
        do l = 1,size(edge_map(this_part)%loc_to_glb)
          write (io,13) "edge",this_part,l,edge_map(this_part)%loc_to_glb(l)
        end do
      else
        write (io,22) "edge",this_part
      end if
    end do
  else
    write (io,21) "edge"
  end if
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Output face_map to file
  !
  write (fname,1) "face",mypnum
  !
  open (newunit=io,file=fname,status="replace",action="write", &
        iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  if (allocated(face_map)) then
    write (io,11) "face",size(face_map)
    do this_part = 1,size(face_map)
      if (allocated(face_map(this_part)%loc_to_glb)) then
        write (io,12) "face",this_part,size(face_map(this_part)%loc_to_glb)
        do l = 1,size(face_map(this_part)%loc_to_glb)
          write (io,13) "face",this_part,l,face_map(this_part)%loc_to_glb(l)
        end do
      else
        write (io,22) "face",this_part
      end if
    end do
  else
    write (io,21) "face"
  end if
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
  !
  ! #######
  ! ## 4 ##
  ! #######
  !
  ! Output cell_map to file
  !
  write (fname,1) "cell",mypnum
  !
  open (newunit=io,file=fname,status="replace",action="write", &
        iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  if (allocated(cell_map)) then
    write (io,11) "cell",size(cell_map)
    do this_part = 1,size(cell_map)
      if (allocated(cell_map(this_part)%loc_to_glb)) then
        write (io,12) "cell",this_part,size(cell_map(this_part)%loc_to_glb)
        do l = 1,size(cell_map(this_part)%loc_to_glb)
          write (io,13) "cell",this_part,l,cell_map(this_part)%loc_to_glb(l)
        end do
      else
        write (io,22) "cell",this_part
      end if
    end do
  else
    write (io,21) "cell"
  end if
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
  !
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__, &
                  "TEMPORARY STOP AFTER WRITING OUT THE MAPPING ARRAYS!!")
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1  format (a,"_map.",i3.3,".dat")
  11 format ("size(",a,"_map) = ",i0)
  12 format (4x,"size(",a,"_map(",i3.3,")%loc_to_glb) = ",i0)
  13 format (8x,a,"_map(",i3.3,")%loc_to_glb(",i8,") = ",i0)
  21 format ("WARNING!!! ",a,"_map IS UNALLOCATED!!!")
  22 format (4x,"WARNING!!! ",a,"_map(",i0,")%loc_to_glb IS UNALLOCATED")
  !
end subroutine output_mapping_arrays
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
  allocate ( return_value(1:ncount) , &
             source=pack( intseq(1,size(array)) , array == find_value ) , &
             stat=ierr )
  !
end function find_subset
!
!###############################################################################
!
subroutine bcast_each_fp_and_map_array(npart,node_map,edge_map,face_map, &
                                       cell_map,flx,face_rotation)
  !
  !.. Formal Arguments ..
  integer,                          intent(in) :: npart
  type(map_t),      allocatable, intent(inout) :: node_map(:)
  type(map_t),      allocatable, intent(inout) :: edge_map(:)
  type(map_t),      allocatable, intent(inout) :: face_map(:)
  type(map_t),      allocatable, intent(inout) :: cell_map(:)
  type(flux_pts_t), allocatable, intent(inout) :: flx(:)
  integer,          allocatable, intent(inout) :: face_rotation(:,:)
  !
  !.. Local Scalars ..
  integer :: n,np,i0,i1,i2,ierr
  !
  character(len=100) :: array_name
  !
  !.. Local Arrays ..
  integer :: array_sizes(1:4*npart+3)
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: buffer(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "bcast_each_fp_and_map_array"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (i_am_host_root) then
    !
    np = 0
    do n = 1,size(flx)
      np = np + size(flx(n)%pts) + 1
    end do
    !
    do n = 1,npart
      array_sizes(n+0*npart) = size(node_map(n)%loc_to_glb)
      array_sizes(n+1*npart) = -1
      if (allocated(edge_map)) then
        if (allocated(edge_map(n)%loc_to_glb)) then
          array_sizes(n+1*npart) = size(edge_map(n)%loc_to_glb)
        end if
      end if
      array_sizes(n+2*npart) = size(face_map(n)%loc_to_glb)
      array_sizes(n+3*npart) = size(cell_map(n)%loc_to_glb)
    end do
    array_sizes(1+4*npart) = size(face_rotation,dim=2)
    array_sizes(2+4*npart) = size(flx)
    array_sizes(3+4*npart) = np
    !
  else
    !
    allocate ( node_map(1:npart) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"node_map",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( edge_map(1:npart) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edge_map",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( face_map(1:npart) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"face_map",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( cell_map(1:npart) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_map",1,__LINE__,__FILE__,ierr,error_message)
    !
  end if
  !
  call mpi_bcast(array_sizes,size(array_sizes,kind=int_mpi),mpi_inttyp, &
                 host_root,my_host_comm,mpierr)
  !
  ! Have all non-root processors pre-allocate the loc_to_glb components
  ! of the mapping arrays, and the face_rotation array
  !
  if (.not. i_am_host_root) then
    !
    do n = 1,size(node_map)
      np = array_sizes(n+0*npart)
      allocate ( node_map(n)%loc_to_glb(1:np) , source=0 , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) "node",n
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    end do
    !
    if (all(array_sizes(npart+1:2*npart) == -1)) then
      if (allocated(edge_map)) then
        deallocate ( edge_map , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"edge_map",2,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
    end if
    !
    if (allocated(edge_map)) then
      do n = 1,size(edge_map)
        np = array_sizes(n+1*npart)
        if (np > 0) then
          allocate ( edge_map(n)%loc_to_glb(1:np) , source=0 , &
                     stat=ierr , errmsg=error_message )
          write (array_name,1) "edge",n
          call alloc_error(pname,array_name,1,__LINE__,__FILE__, &
                           ierr,error_message)
        end if
      end do
    end if
    !
    do n = 1,size(face_map)
      np = array_sizes(n+2*npart)
      allocate ( face_map(n)%loc_to_glb(1:np) , source=0 , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) "face",n
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    end do
    !
    do n = 1,size(cell_map)
      np = array_sizes(n+3*npart)
      allocate ( cell_map(n)%loc_to_glb(1:np) , source=0 , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) "cell",n
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    end do
    !
    np = array_sizes(1+4*npart)
    allocate ( face_rotation(1:2,1:np) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"face_rotation",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
  end if
  !
  ! Now broadcast all these arrays individually (to prevent having to
  ! allocate a bunch of extra memory for buffers) to all processors
  !
  ! Broadcast each node_map(:)%loc_to_glb array to all processors
  !
  do n = 1,size(node_map)
    call mpi_bcast(node_map(n)%loc_to_glb, &
                   size(node_map(n)%loc_to_glb,kind=int_mpi), &
                   mpi_inttyp,host_root,my_host_comm,mpierr)
  end do
  !
  ! Broadcast each edge_map(:)%loc_to_glb array to all processors
  !
  if (allocated(edge_map)) then
    do n = 1,size(edge_map)
      if (allocated(edge_map(n)%loc_to_glb)) then
        call mpi_bcast(edge_map(n)%loc_to_glb, &
                       size(edge_map(n)%loc_to_glb,kind=int_mpi), &
                       mpi_inttyp,host_root,my_host_comm,mpierr)
      end if
    end do
  end if
  !
  ! Broadcast each face_map(:)%loc_to_glb array to all processors
  !
  do n = 1,size(face_map)
    call mpi_bcast(face_map(n)%loc_to_glb, &
                   size(face_map(n)%loc_to_glb,kind=int_mpi), &
                   mpi_inttyp,host_root,my_host_comm,mpierr)
  end do
  !
  ! Broadcast each cell_map(:)%loc_to_glb array to all processors
  !
  do n = 1,size(cell_map)
    call mpi_bcast(cell_map(n)%loc_to_glb, &
                   size(cell_map(n)%loc_to_glb,kind=int_mpi), &
                   mpi_inttyp,host_root,my_host_comm,mpierr)
  end do
  !
  ! Broadcast the face_rotation array to all processors
  !
  call mpi_bcast(face_rotation,size(face_rotation,kind=int_mpi), &
                 mpi_inttyp,host_root,my_host_comm,mpierr)
  !
  ! Allocate the buffer for broadcasting the components of the flx array
  !
  np = array_sizes(3+4*npart)
  allocate ( buffer(1:np) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"buffer",1,__LINE__,__FILE__,ierr,error_message)
  !
  if (i_am_host_root) then
    !
    ! Have the root processor create the buffer for the flux array
    !
    i2 = 0
    !
    do n = 1,size(flx)
      np = size(flx(n)%pts)
      i1 = i2 + 1
      i2 = i2 + 1 + np
      buffer(i1:i2) = [ np , flx(n)%pts(:) ]
    end do
    !
  end if
  !
  ! Broadcast the buffer for the flx array to all processors
  !
  call mpi_bcast(buffer,size(buffer,kind=int_mpi),mpi_inttyp, &
                 host_root,my_host_comm,mpierr)
  !
  ! Have all non-root processors unpack the buffer for the flx array
  !
  if (.not. i_am_host_root) then
    !
    ! First, allocate the flux array
    !
    np = array_sizes(2+4*npart)
    allocate ( flx(1:np) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"flx",1,__LINE__,__FILE__,ierr,error_message)
    !
    i2 = 0
    !
    do n = 1,size(flx)
      !
      i0 = i2 + 1
      np = buffer(i0)
      !
      i1 = i0 + 1
      i2 = i0 + np
#ifdef SPECIAL_FOR_PGI
     !flx(n)%pts = buffer(i1:i2) ! using F2003 auto-reallocation
      flx(n)%pts = buffer(i1-1 + intseq(1,np)) ! using F2003 auto-reallocation
#else
      flx(n)%pts = buffer(i1:i2) ! using F2003 auto-reallocation
#endif
      !
    end do
    !
  end if
  !
  ! Deallocate the buffer array before we leave
  !
  deallocate ( buffer , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"buffer",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (a,"_map(",i0,")%loc_to_glb")
  !
end subroutine bcast_each_fp_and_map_array
!
!###############################################################################
!
subroutine root_localizes_grid(nodes_of_cell,nodes_of_cell_ptr,bface, &
                               xyz_nodes,cell_geom,cell_order, &
                               face,node_map,face_map,cell_map)
  !
  !.. Use Statements ..
  use geovar, only : grid_t,grid
  !
  !.. Formal Arguments ..
  type(face_t), allocatable,    intent(in) :: face(:)
  type(map_t),                  intent(in) :: node_map(:)
  type(map_t),                  intent(in) :: face_map(:)
  type(map_t),                  intent(in) :: cell_map(:)
  integer,      allocatable, intent(inout) :: cell_geom(:)
  integer,      allocatable, intent(inout) :: cell_order(:)
  integer,      allocatable, intent(inout) :: nodes_of_cell(:)
  integer,      allocatable, intent(inout) :: nodes_of_cell_ptr(:)
  integer,      allocatable, intent(inout) :: bface(:,:)
  real(wp),     allocatable, intent(inout) :: xyz_nodes(:,:)
  !
  !.. Local Scalars ..
  integer :: n,my_part,this_part,lcell,lnode,gnode,ierr
  !
  character(len=100) :: part_message
  !
  !.. Local Allocatable Arrays ..
  integer,      allocatable :: global_geom(:)
  integer,      allocatable :: global_order(:)
  integer,      allocatable :: global_nofc(:)
  integer,      allocatable :: global_nofc_ptr(:)
  integer,      allocatable :: global_bface(:,:)
  real(wp),     allocatable :: global_xyz(:,:)
  integer,      allocatable :: grid_elem(:)
  real(wp),     allocatable :: grid_xyz(:,:)
  type(grid_t), allocatable :: global_grid
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "root_localizes_grid"
  !
continue
  !
  part_message = ""
  !
  call debug_timer(entering_procedure,pname)
  !
  if (i_am_host_root) then
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
    call move_alloc( from=grid              , to=global_grid )
    !
    gnode = size(global_xyz,dim=2)
    !
  end if
  !
  ! Have all processors deallocate these global
  ! arrays if they are currently allocated.
  !
  ! NOTE: The above calls to move_alloc should have already
  !       accomplished this task for the root processor
  !
  if (allocated(cell_geom)) then
    deallocate ( cell_geom , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_geom",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(cell_order)) then
    deallocate ( cell_order , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_order",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (allocated(nodes_of_cell)) then
    deallocate ( nodes_of_cell , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodes_of_cell",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (allocated(nodes_of_cell_ptr)) then
    deallocate ( nodes_of_cell_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodes_of_cell_ptr",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (allocated(bface)) then
    deallocate ( bface , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bface",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(xyz_nodes)) then
    deallocate ( xyz_nodes , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_nodes",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(grid)) then
    deallocate ( grid , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! The process ID, mypnum, starts numbering at 0 so
  ! identify which partition is local to each processor
  !
  my_part = mypnum + 1
  !
  ! Have the root processor localize each partition and
  ! then send the localized arrays to the matching processor
  !
  do n = 2,size(my_host_grp)
    !
    this_part = my_host_grp(n)+1
    !
    write (part_message,11) this_part
    call debug_timer(start_timer,part_message)
    !
    if (i_am_host_root) then
      !
#ifdef DEBUG_PARALLEL
      write (iout,10) this_part
#endif
      !
      ! Localize the bface array
      !
      call root_localizes_bface(this_part,cell_map,face_map,global_bface,bface)
      !
      ! Localize the cell_geom and cell_order arrays
      !
      call root_localizes_cell_geom_and_order(this_part,face, &
                                              face_map,cell_map, &
                                              global_geom,global_order, &
                                              cell_geom,cell_order)
      !
      ! Localize the nodes_of_cell_ptr and nodes_of_cell arrays
      !
      call root_localizes_nodes_of_cell(this_part,gnode,bface,face, &
                                        node_map,face_map,cell_map, &
                                        global_nofc,global_nofc_ptr, &
                                        nodes_of_cell,nodes_of_cell_ptr)
      !
      ! Localize the xyz_nodes array
      !
      call root_localizes_xyz_nodes(this_part,node_map,global_xyz,xyz_nodes)
      !
      ! Localize the grid derived type
      !
      call root_localizes_grid_dt(this_part,cell_map,global_grid, &
                                  grid_elem,grid_xyz)
      !
    end if
    !
    if (i_am_host_root .or. my_part == this_part) then
      !
     !if (use_mpi_pack) then
     !  call pack_localized_partition(this_part,bface,cell_geom,cell_order, &
     !                                nodes_of_cell,nodes_of_cell_ptr, &
     !                                xyz_nodes,grid_elem,grid_xyz)
     !else
     !  if (send_using_large_buffer) then
     !    call bsend_localized_partition(this_part,bface,cell_geom,cell_order, &
     !                                   nodes_of_cell,nodes_of_cell_ptr, &
     !                                   xyz_nodes,grid_elem,grid_xyz)
     !  else
          call send_localized_partition(this_part,bface,cell_geom,cell_order, &
                                        nodes_of_cell,nodes_of_cell_ptr, &
                                        xyz_nodes,grid_elem,grid_xyz)
     !  end if
     !end if
      !
    end if
    !
    ! Unpack the grid_elem and grid_xyz arrays into grid
    !
    if (my_part == this_part) then
      call unpack_localized_grid_dt(size(cell_map(this_part)%loc_to_glb), &
                                    grid_elem,grid_xyz,grid)
    end if
    !
    call debug_timer(stop_timer,part_message)
    !
  end do
  !
  ! Now have the root processor localize its own partition
  !
  if (i_am_host_root) then
    !
    this_part = mypnum + 1
    !
    write (part_message,11) this_part
    call debug_timer(start_timer,part_message)
    !
#ifdef DEBUG_PARALLEL
    write (iout,10) this_part
#endif
    !
    ! Localize the bface array
    !
    call root_localizes_bface(this_part,cell_map,face_map,global_bface,bface)
    !
    ! Localize the cell_geom and cell_order arrays
    !
    call root_localizes_cell_geom_and_order(this_part,face, &
                                            face_map,cell_map, &
                                            global_geom,global_order, &
                                            cell_geom,cell_order)
    !
    ! Localize the nodes_of_cell_ptr and nodes_of_cell arrays
    !
    call root_localizes_nodes_of_cell(this_part,gnode,bface,face, &
                                      node_map,face_map,cell_map, &
                                      global_nofc,global_nofc_ptr, &
                                      nodes_of_cell,nodes_of_cell_ptr)
    !
    ! Localize the xyz_nodes array
    !
    call root_localizes_xyz_nodes(this_part,node_map,global_xyz,xyz_nodes)
    !
    ! Localize the grid derived type
    !
    call root_localizes_grid_dt(this_part,cell_map,global_grid, &
                                grid_elem,grid_xyz)
    !
    ! Unpack the grid_elem and grid_xyz arrays into grid
    !
    call unpack_localized_grid_dt(size(cell_map(this_part)%loc_to_glb), &
                                  grid_elem,grid_xyz,grid)
    !
    call debug_timer(stop_timer,part_message)
    !
  end if
  !
  ! We no longer need the global arrays so deallocate
  !
  if (allocated(global_grid)) then
    deallocate ( global_grid , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_grid",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(global_bface)) then
    deallocate ( global_bface , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_bface",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (allocated(global_geom)) then
    deallocate ( global_geom , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_geom",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(global_order)) then
    deallocate ( global_order , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_order",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (allocated(global_nofc)) then
    deallocate ( global_nofc , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_nofc",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(global_nofc_ptr)) then
    deallocate ( global_nofc_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_nofc_ptr",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (allocated(global_xyz)) then
    deallocate ( global_xyz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_xyz",2,__LINE__,__FILE__,ierr,error_message)
  end if
 !!
 !! Output the individual partitions to separate Tecplot files
 !!
 !lcell = size(cell_map(my_part)%loc_to_glb)
 !lnode = size(node_map(my_part)%loc_to_glb)
 !!
 !call output_partitions(my_part,nodes_of_cell,nodes_of_cell_ptr,xyz_nodes)
 !!
 !call output_partitions(my_part+ncpu+1,nodes_of_cell, &
 !                       nodes_of_cell_ptr(1:lcell+1), &
 !                       xyz_nodes(1:nr,1:lnode))
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  10 format (/,"For partition # ",i0)
  11 format ("Root localizing partition # ",i0)
  !
end subroutine root_localizes_grid
!
!###############################################################################
!
subroutine unpack_localized_grid_dt(lcell,grid_elem,grid_xyz,grid)
  !
  !.. Use Statements ..
  use geovar,        only : grid_t
  use ovar,          only : grid_format,plot3d_agglomerate_order
  use module_cgns,   only : cgns_element_properties
  use module_gmsh,   only : gmsh_element_properties
  use module_plot3d, only : plot3d_element_properties
  !
  !.. Formal Arguments ..
  integer,                      intent(in) :: lcell
  integer,      allocatable, intent(inout) :: grid_elem(:)
  real(wp),     allocatable, intent(inout) :: grid_xyz(:,:)
  type(grid_t), allocatable, intent(inout) :: grid
  !
  !.. Local Scalars ..
  integer :: n,nc,nn,np,n1,n2,ierr
  integer :: i,ntags,elem_type,plot3d_order
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "unpack_localized_grid_dt"
  !
continue
  !
  ! First allocate the grid derived type
  !
  if (allocated(grid)) then
    deallocate ( grid , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( grid , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the grid%xyz array
  !
  ! USE MOVE_ALLOC (AVOIDS COPYING DATA TO NEW MEMORY)
  call move_alloc( from=grid_xyz , to=grid%xyz )
  ! USE STANDARD ALLOCATION (COPIES DATA TO NEW MEMORY)
 !allocate ( grid%xyz , source=grid_xyz , stat=ierr , errmsg=error_message )
 !call alloc_error(pname,"grid%xyz",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the grid%elem array
  !
  allocate ( grid%elem(1:lcell) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid%elem",1,__LINE__,__FILE__,ierr,error_message)
  !
  n = 0
  do nc = 1,lcell
    !
    ! Unpack the elem_type/ntags combo
    !
    i = 1
    elem_type = mod( grid_elem(n+i) , 1000 )
    ntags = grid_elem(n+i) / 1000
    !
    ! Get the element properties for this type of element
    !
    if (grid_format == Gmsh_Format) then
      grid%elem(nc)%prop = gmsh_element_properties(elem_type)
    else if (grid_format == CGNS_Format) then
      grid%elem(nc)%prop = cgns_element_properties(elem_type)
    else if (grid_format == Plot3D_Format) then
      plot3d_order = max(1,plot3d_agglomerate_order)
      grid%elem(nc)%prop = plot3d_element_properties(elem_type,plot3d_order)
    else
      grid%elem(nc)%prop = linear_element_properties(elem_type)
    end if
    !
    ! Unpack the host cell if this is being used
    !
    if (send_grid_elem_host_cell) then
      i = i + 1
      grid%elem(nc)%host_cell = grid_elem(n+i)
    end if
    !
    ! Unpack the point indices from grid_elem if the next value in grid_elem
    ! is not 0; otherwise skip this step and increase the counter for this
    ! element by 1
    !
    if (grid_elem(n+i+1) /= 0) then
      np = grid%elem(nc)%prop%npts
      n1 = n + i + 1
      n2 = n + i + np
#ifdef SPECIAL_FOR_PGI
      grid%elem(nc)%pts = grid_elem(n1-1+intseq(1,np)) ! F2003 AUTO-REALLOC
#else
      grid%elem(nc)%pts = grid_elem(n1:n2) ! using F2003 auto-reallocation
#endif
      i = i + np
    else
      i = i + 1
    end if
    !
    ! Unpack the tags from grid_elem
    !
    n1 = n + i + 1
    n2 = n + i + ntags
#ifdef SPECIAL_FOR_PGI
    grid%elem(nc)%tags = grid_elem(n1-1+intseq(1,ntags)) ! F2003 AUTO-REALLOC
#else
    grid%elem(nc)%tags = grid_elem(n1:n2) ! using F2003 auto-reallocation
#endif
    i = i + ntags
    !
    ! Unpack the node indices from grid_elem if the next value in grid_elem
    ! is not 0; otherwise skip this step and increase the counter for this
    !
    if (grid_elem(n+i+1) /= 0) then
      nn = geom_nodes( grid%elem(nc)%prop%geom )
      n1 = n + i + 1
      n2 = n + i + nn
#ifdef SPECIAL_FOR_PGI
      grid%elem(nc)%nodes = grid_elem(n1-1+intseq(1,nn)) ! F2003 AUTO-REALLOC
#else
      grid%elem(nc)%nodes = grid_elem(n1:n2) ! using F2003 auto-reallocation
#endif
      i = i + nn
    else
      i = i + 1
    end if
    !
    ! Add the final counter for this element to the grid_elem
    ! counter before going to the next element
    !
    n = n + i
    !
  end do
  !
  ! Check to make sure the final counter matches the size of the grid_elem array
  !
  if (n /= size(grid_elem)) then
    write (error_message,1) mypnum+1,n,size(grid_elem)
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Deallocate grid_elem and grid_xyz because they are no longer needed
  !
  if (allocated(grid_elem)) then
    deallocate ( grid_elem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid_elem",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(grid_xyz)) then
    deallocate ( grid_xyz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid_xyz",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Format Statements
  !
  1 format ("ERROR CREATING grid_elem FOR PARTITION #",i0,"!!! ", &
            "THE FINAL COUNTER n=",i0,", WHEN IT SHOULD BE ", &
            "size(grid_elem)=",i0,"!!!")
  !
end subroutine unpack_localized_grid_dt
!
!###############################################################################
!
pure function linear_element_properties(elem_type) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : prop_t
  !
  !.. Formal Arguments ..
  integer, intent(in) :: elem_type
  !
  !.. Function Result ..
  type(prop_t) :: return_value
  !
continue
  !
 !return_value = prop_t(npts=geom_nodes(elem_type), &
 !                      order=1, &
 !                      geom=elem_type, &
 !                      elem_type=elem_type)
  return_value%npts = geom_nodes(elem_type)
  return_value%order = 1
  return_value%geom = elem_type
  return_value%elem_type = elem_type
  !
end function linear_element_properties
!
!###############################################################################
!
subroutine send_localized_partition(this_part,bface,cell_geom,cell_order, &
                                    nodes_of_cell,nodes_of_cell_ptr, &
                                    xyz_nodes,grid_elem,grid_xyz)
  !
  !.. Formal Arguments ..
  integer,                  intent(in) :: this_part
  integer,  allocatable, intent(inout) :: cell_geom(:)
  integer,  allocatable, intent(inout) :: cell_order(:)
  integer,  allocatable, intent(inout) :: nodes_of_cell(:)
  integer,  allocatable, intent(inout) :: nodes_of_cell_ptr(:)
  integer,  allocatable, intent(inout) :: bface(:,:)
  real(wp), allocatable, intent(inout) :: xyz_nodes(:,:)
  integer,  allocatable, intent(inout) :: grid_elem(:)
  real(wp), allocatable, intent(inout) :: grid_xyz(:,:)
  !
  !.. Local Scalars ..
  integer :: n1,n2,ierr
  integer(INT_MPI) :: send_cpu,recv_cpu
  integer(INT_MPI) :: itag1,itag2,itag3,itag4,itag5
  integer(INT_MPI) :: itag6,itag7,itag8,itag9,nsz
  !
  !.. Local Arrays ..
  integer :: array_sizes(1:13)
  _MPI_STATUS_TYPE_ :: istat(1:_MPI_STATUS_SIZE_)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "send_localized_partition"
  !
continue
  !
  call debug_timer(entering_procedure,pname,processor=0,part_num=this_part)
  !
  send_cpu = my_host_root
  recv_cpu = int(this_part-1,kind=int_mpi)
  itag1 = int(this_part,kind=int_mpi)
  itag2 = itag1 + ncpu
  itag3 = itag2 + ncpu
  itag4 = itag3 + ncpu
  itag5 = itag4 + ncpu
  itag6 = itag5 + ncpu
  itag7 = itag6 + ncpu
  itag8 = itag7 + ncpu
  itag9 = itag8 + ncpu
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Create and send array_sizes to the processor handling this partition
  !
  if (mypnum == send_cpu) then
    !
    array_sizes( 1) = size(cell_geom)
    array_sizes( 2) = merge( 1 , size(cell_geom), &
                            all(cell_geom(1)==cell_geom(2:)) )
    array_sizes( 3) = size(cell_order)
    array_sizes( 4) = merge( 1 , size(cell_order) , &
                            all(cell_order(1)==cell_order(2:)) )
    array_sizes( 5) = size(nodes_of_cell)
    array_sizes( 6) = size(nodes_of_cell_ptr)
    array_sizes( 7) = size(grid_elem)
    array_sizes( 8) = size(bface,dim=1)
    array_sizes( 9) = size(bface,dim=2)
    array_sizes(10) = size(xyz_nodes,dim=1)
    array_sizes(11) = size(xyz_nodes,dim=2)
    array_sizes(12) = size(grid_xyz,dim=1)
    array_sizes(13) = size(grid_xyz,dim=2)
    !
#ifdef DEBUG_PARALLEL
    write (iout,11) "cell_geom",         size(cell_geom)
    write (iout,11) "cell_order",        size(cell_order)
    write (iout,11) "nodes_of_cell",     size(nodes_of_cell)
    write (iout,11) "nodes_of_cell_ptr", size(nodes_of_cell_ptr)
    write (iout,11) "grid_elem",         size(grid_elem)
    write (iout,12) "bface",1,           size(bface,dim=1)
    write (iout,12) "bface",2,           size(bface,dim=2)
    write (iout,11) "bface",             size(bface)
    write (iout,12) "xyz_nodes",1,       size(xyz_nodes,dim=1)
    write (iout,12) "xyz_nodes",2,       size(xyz_nodes,dim=2)
    write (iout,11) "xyz_nodes",         size(xyz_nodes)
    write (iout,12) "grid_xyz",1,        size(grid_xyz,dim=1)
    write (iout,12) "grid_xyz",2,        size(grid_xyz,dim=2)
    write (iout,11) "grid_xyz",          size(grid_xyz)
#endif
    !
    call mpi_send(array_sizes,size(array_sizes,kind=int_mpi), &
                  mpi_inttyp,recv_cpu,itag1,MPI_COMM_WORLD,mpierr)
    !
  else if (mypnum == recv_cpu) then
    !
    call mpi_recv(array_sizes,size(array_sizes,kind=int_mpi), &
                  mpi_inttyp,send_cpu,itag1,MPI_COMM_WORLD,istat(1),mpierr)
    !
  end if
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Send the cell_geom array to the receiving processor
  !
  nsz = int(array_sizes(2),kind=int_mpi)
  !
  if (mypnum == send_cpu) then
    !
    call mpi_send(cell_geom(1:nsz),nsz,mpi_inttyp, &
                  recv_cpu,itag2,MPI_COMM_WORLD,mpierr)
#ifdef DEBUG_PARALLEL
    write (iout,11) "cell_geom",size(cell_geom), &
                    minval(cell_geom,dim=1),maxval(cell_geom),"send"
#endif
    !
  else if (mypnum == recv_cpu) then
    !
    n1 = array_sizes(1)
    allocate ( cell_geom(1:n1) , source=0 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_geom",1,__LINE__,__FILE__,ierr,error_message)
    !
    call mpi_recv(cell_geom(1:nsz),nsz,mpi_inttyp, &
                  send_cpu,itag2,MPI_COMM_WORLD,istat(1),mpierr)
    !
    if (nsz == 1) cell_geom(2:) = cell_geom(1)
#ifdef DEBUG_PARALLEL
    write (iout,11) "cell_geom",size(cell_geom), &
                    minval(cell_geom,dim=1),maxval(cell_geom),"recv"
#endif

    !
  end if
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Send the cell_order array to the receiving processor
  !
  nsz = int(array_sizes(4),kind=int_mpi)
  !
  if (mypnum == send_cpu) then
    !
    call mpi_send(cell_order(1:nsz),nsz,mpi_inttyp, &
                  recv_cpu,itag3,MPI_COMM_WORLD,mpierr)
#ifdef DEBUG_PARALLEL
    write (iout,11) "cell_order",size(cell_order), &
                    minval(cell_order,dim=1),maxval(cell_order),"send"
#endif
    !
  else if (mypnum == recv_cpu) then
    !
    n1 = array_sizes(3)
    allocate ( cell_order(1:n1) , source=0 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_order",1,__LINE__,__FILE__,ierr,error_message)
    !
    call mpi_recv(cell_order(1:nsz),nsz,mpi_inttyp, &
                  send_cpu,itag3,MPI_COMM_WORLD,istat(1),mpierr)
    !
    if (nsz == 1) cell_order(2:) = cell_order(1)
#ifdef DEBUG_PARALLEL
    write (iout,11) "cell_order",size(cell_order), &
                    minval(cell_order,dim=1),maxval(cell_order),"recv"
#endif
    !
  end if
  !
  ! #######
  ! ## 4 ##
  ! #######
  !
  ! Send the nodes_of_cell array to the receiving processor
  !
  if (mypnum == send_cpu) then
    !
    call mpi_send(nodes_of_cell,size(nodes_of_cell,kind=int_mpi), &
                  mpi_inttyp,recv_cpu,itag4,MPI_COMM_WORLD,mpierr)
    !
  else if (mypnum == recv_cpu) then
    !
    n1 = array_sizes(5)
    allocate ( nodes_of_cell(1:n1) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    call mpi_recv(nodes_of_cell,size(nodes_of_cell,kind=int_mpi), &
                  mpi_inttyp,send_cpu,itag4,MPI_COMM_WORLD,istat(1),mpierr)
    !
  end if
  !
  ! #######
  ! ## 5 ##
  ! #######
  !
  ! Send the nodes_of_cell_ptr array to the receiving processor
  !
  if (mypnum == send_cpu) then
    !
    call mpi_send(nodes_of_cell_ptr,size(nodes_of_cell_ptr,kind=int_mpi), &
                  mpi_inttyp,recv_cpu,itag5,MPI_COMM_WORLD,mpierr)
    !
  else if (mypnum == recv_cpu) then
    !
    n1 = array_sizes(6)
    allocate ( nodes_of_cell_ptr(1:n1) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    call mpi_recv(nodes_of_cell_ptr,size(nodes_of_cell_ptr,kind=int_mpi), &
                  mpi_inttyp,send_cpu,itag5,MPI_COMM_WORLD,istat(1),mpierr)
    !
  end if
  !
  ! #######
  ! ## 6 ##
  ! #######
  !
  ! Send the grid_elem array to the receiving processor
  !
  if (mypnum == send_cpu) then
    !
    call mpi_send(grid_elem,size(grid_elem,kind=int_mpi), &
                  mpi_inttyp,recv_cpu,itag6,MPI_COMM_WORLD,mpierr)
    !
  else if (mypnum == recv_cpu) then
    !
    n1 = array_sizes(7)
    allocate ( grid_elem(1:n1) , source=0 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid_elem",1,__LINE__,__FILE__,ierr,error_message)
    !
    call mpi_recv(grid_elem,size(grid_elem,kind=int_mpi),mpi_inttyp, &
                  send_cpu,itag6,MPI_COMM_WORLD,istat(1),mpierr)
    !
  end if
  !
  ! #######
  ! ## 7 ##
  ! #######
  !
  ! Send the bface array to the receiving processor
  !
  if (mypnum == send_cpu) then
    !
    call mpi_send(bface,size(bface,kind=int_mpi), &
                  mpi_inttyp,recv_cpu,itag7,MPI_COMM_WORLD,mpierr)
    !
  else if (mypnum == recv_cpu) then
    !
    n1 = array_sizes(8)
    n2 = array_sizes(9)
    allocate ( bface(1:n1,1:n2) , source=0 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
    !
    call mpi_recv(bface,size(bface,kind=int_mpi),mpi_inttyp, &
                  send_cpu,itag7,MPI_COMM_WORLD,istat(1),mpierr)
    !
  end if
  !
  ! #######
  ! ## 8 ##
  ! #######
  !
  ! Send the xyz_nodes array to the receiving processor
  !
  if (mypnum == send_cpu) then
    !
    call mpi_send(xyz_nodes,size(xyz_nodes,kind=int_mpi), &
                  mpi_flttyp,recv_cpu,itag8,MPI_COMM_WORLD,mpierr)
    !
  else if (mypnum == recv_cpu) then
    !
    n1 = array_sizes(10)
    n2 = array_sizes(11)
    allocate ( xyz_nodes(1:n1,1:n2) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_nodes",1,__LINE__,__FILE__,ierr,error_message)
    !
    call mpi_recv(xyz_nodes,size(xyz_nodes,kind=int_mpi),mpi_flttyp, &
                  send_cpu,itag8,MPI_COMM_WORLD,istat(1),mpierr)
    !
  end if
  !
  ! #######
  ! ## 9 ##
  ! #######
  !
  ! Send the grid_xyz array to the receiving processor
  !
  if (mypnum == send_cpu) then
    !
    call mpi_send(grid_xyz,size(grid_xyz,kind=int_mpi), &
                  mpi_flttyp,recv_cpu,itag9,MPI_COMM_WORLD,mpierr)
    !
  else if (mypnum == recv_cpu) then
    !
    n1 = array_sizes(12)
    n2 = array_sizes(13)
    allocate ( grid_xyz(1:n1,1:n2) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid_xyz",1,__LINE__,__FILE__,ierr,error_message)
    !
    call mpi_recv(grid_xyz,size(grid_xyz,kind=int_mpi),mpi_flttyp, &
                  send_cpu,itag9,MPI_COMM_WORLD,istat(1),mpierr)
    !
  end if
  !
  call debug_timer(leaving_procedure,pname,processor=0,part_num=this_part)
  !
  ! Format Statements
  !
  11 format (5x,"size(",a,")",t29," = ",i11,:," = [",i0,",",i0,"]  - ",a)
  12 format (5x,"size(",a,",dim=",i0,")",t29," = ",i11)
  13 format (5x,"Adding ",a," to the buffer")
  14 format (5x,"Sending ",a,:," = ",i0)
  !
end subroutine send_localized_partition
!
!###############################################################################
!
subroutine root_localizes_bface(this_part,cell_map,face_map,global_bface,bface)
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: this_part
  type(map_t),             intent(in) :: face_map(:)
  type(map_t),             intent(in) :: cell_map(:)
  integer,                 intent(in) :: global_bface(:,:)
  integer, allocatable, intent(inout) :: bface(:,:)
  !
  !.. Local Scalars ..
  integer :: ierr
  integer :: gfbnd,lfbnd,lcell,nbfai
  integer :: loc_face,glb_face
  integer :: adj_face,adj_cell,adj_part
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "root_localizes_bface"
  !
continue
  !
  call debug_timer(entering_procedure,pname,processor=0,part_num=this_part)
  !
  ! Before we deallocate any of the input arrays, get the size of the global
  ! grid from these arrays
  !
  gfbnd = size(global_bface,dim=2)
  !
  ! Count the local number of cells and boundary faces
  !
  lcell = size( cell_map(this_part)%loc_to_glb )
  lfbnd = count( face_map(this_part)%loc_to_glb <= gfbnd )
  !
#ifdef DEBUG_PARALLEL
  write (iout,10) pname,gfbnd,lcell,lfbnd
#endif
  10 format (3x,"In ",a," :",/,10x,"gfbnd = ",i0,/, &
                               10x,"lcell = ",i0,/, &
                               10x,"lfbnd = ",i0)
  !
  ! Localize the boundary face array, but not the information it contains.
  !
  nbfai = size(global_bface,dim=1)
  !
  if (allocated(bface)) then
    deallocate ( bface , stat=ierr , errmsg=error_message)
    call alloc_error(pname,"bface",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( bface(1:nbfai,1:lfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
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
    glb_face = face_map(this_part)%loc_to_glb(loc_face)
    !
    ! Copy everything to the local bface array
    !
    bface(:,loc_face) = global_bface(:,glb_face)
    !
    ! Now localize the information for this face
    !
    ! local host cell
    bface(2,loc_face) = findloc( cell_map(this_part)%loc_to_glb , &
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
    ! NOTE: Make sure this is negative if global_bface(10,glb_face) is
    !       also negative
    !
    bface(10,loc_face) = sign( lcell+loc_face , global_bface(10,glb_face) )
    !
    ! Now that we have localized the information for this boundary face,
    ! change any periodic boundary conditions that require communication
    ! across processors to communication boundaries.
    ! NOTE: Set the value to -bc_cpu_bnd to identify this face as a former
    !       periodic boundary condition face.
    !
    if (bface(1,loc_face) == bc_periodic) then
      if (bface(7,loc_face) /= this_part-1) then
        bface(1,loc_face) = -bc_cpu_bnd
      end if
    end if
    !
  end do
  !
  call debug_timer(leaving_procedure,pname,processor=0,part_num=this_part)
  !
end subroutine root_localizes_bface
!
!###############################################################################
!
subroutine root_localizes_cell_geom_and_order(this_part,face, &
                                              face_map,cell_map, &
                                              global_geom,global_order, &
                                              cell_geom,cell_order)
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: this_part
  type(face_t),            intent(in) :: face(:)
  type(map_t),             intent(in) :: face_map(:)
  type(map_t),             intent(in) :: cell_map(:)
  integer,                 intent(in) :: global_geom(:)
  integer,                 intent(in) :: global_order(:)
  integer, allocatable, intent(inout) :: cell_geom(:)
  integer, allocatable, intent(inout) :: cell_order(:)
  !
  !.. Local Scalars ..
  integer :: ierr
  integer :: gcell,lcell
  integer :: gfbnd,lfbnd
  integer :: loc_face,loc_cell
  integer :: glb_face,glb_cell
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "root_localizes_cell_geom_and_order"
  !
continue
  !
  call debug_timer(entering_procedure,pname,processor=0,part_num=this_part)
  !
  ! Before we deallocate any of the input arrays, get the size of the global
  ! grid from these arrays
  !
  gcell = sum_map_sizes(cell_map)
  gfbnd = size(global_geom) - gcell
  !
  ! Count the local number of cells and boundary faces
  !
  lcell = size( cell_map(this_part)%loc_to_glb )
  lfbnd = count( face_map(this_part)%loc_to_glb <= gfbnd )
  !
#ifdef DEBUG_PARALLEL
  write (iout,10) pname, &
                  size(global_geom), &
                  minval(global_geom,dim=1),maxval(global_geom,dim=1), &
                  size(global_order), &
                  minval(global_order,dim=1),maxval(global_order,dim=1), &
                  gcell,gfbnd,lcell,lfbnd
#endif
  10 format (3x,"In ",a," :",/, &
             10x,"size(global_geom) = ",i0," = [",i0,",",i0,"]",/, &
             10x,"size(global_order) = ",i0," = [",i0,",",i0,"]",/, &
             10x,"gcell = ",i0,/,10x,"gfbnd = ",i0,/, &
             10x,"lcell = ",i0,/,10x,"lfbnd = ",i0)
  11 format (3x,"In ",a," :",/, &
             10x,"size(cell_geom) = ",i0," = [",i0,",",i0,"]",/, &
             10x,"size(cell_order) = ",i0," = [",i0,",",i0,"]")
  !
  ! Localize the cell_geom and cell_order arrays
  !
  if (allocated(cell_geom)) then
    deallocate ( cell_geom , stat=ierr , errmsg=error_message)
    call alloc_error(pname,"cell_geom",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( cell_geom(1:lcell+lfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_geom",1,__LINE__,__FILE__,ierr,error_message)
  !
  if (allocated(cell_order)) then
    deallocate ( cell_order , stat=ierr , errmsg=error_message)
    call alloc_error(pname,"cell_order",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( cell_order(1:lcell+lfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_order",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Get the geometry and order for the interior cells local to this partition
  !
  do loc_cell = 1,lcell
    glb_cell = cell_map(this_part)%loc_to_glb(loc_cell)
    cell_geom(loc_cell) = global_geom(glb_cell)
    cell_order(loc_cell) = global_order(glb_cell)
  end do
  !
  ! Get the geometry and order for the ghost cells of this partition
  !
  do loc_face = 1,lfbnd
    glb_face = face_map(this_part)%loc_to_glb(loc_face)
    loc_cell = lcell + loc_face
    cell_geom(loc_cell) = face(glb_face)%geom
    cell_order(loc_cell) = face(glb_face)%order
  end do
#ifdef DEBUG_PARALLEL
  write (iout,11) pname, &
                  size(cell_geom), &
                  minval(cell_geom,dim=1),maxval(cell_geom,dim=1), &
                  size(cell_order), &
                  minval(cell_order,dim=1),maxval(cell_order,dim=1)
#endif
  !
  call debug_timer(leaving_procedure,pname,processor=0,part_num=this_part)
  !
end subroutine root_localizes_cell_geom_and_order
!
!###############################################################################
!
subroutine root_localizes_nodes_of_cell(this_part,gnode,bface,face, &
                                        node_map,face_map,cell_map, &
                                        global_nofc,global_nofc_ptr, &
                                        nodes_of_cell,nodes_of_cell_ptr)
  !
  !.. Formal Arguments ..
  integer,                 intent(in) :: this_part
  integer,                 intent(in) :: gnode
  integer,                 intent(in) :: bface(:,:)
  type(face_t),            intent(in) :: face(:)
  type(map_t),             intent(in) :: node_map(:)
  type(map_t),             intent(in) :: face_map(:)
  type(map_t),             intent(in) :: cell_map(:)
  integer,                 intent(in) :: global_nofc(:)
  integer,                 intent(in) :: global_nofc_ptr(:)
  integer, allocatable, intent(inout) :: nodes_of_cell(:)
  integer, allocatable, intent(inout) :: nodes_of_cell_ptr(:)
  !
  !.. Local Scalars ..
  integer :: ierr,l,l1,l2,g1,g2
  integer :: lfbnd,lcell,nfnods
  integer :: loc_face,loc_cell,loc_node
  integer :: glb_face,glb_cell,glb_node
  integer :: loc_ghost_cell
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: rev_node_map(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "root_localizes_nodes_of_cell"
  !
continue
  !
  call debug_timer(entering_procedure,pname,processor=0,part_num=this_part)
  !
  ! Count the local number of cells and boundary faces
  !
  lcell = size( cell_map(this_part)%loc_to_glb )
  lfbnd = size( bface , dim=2 )
  !
#ifdef DEBUG_PARALLEL
  write (iout,10) pname,lcell,lfbnd
#endif
  10 format (3x,"In ",a," :",/,10x,"lcell = ",i0,/,10x,"lfbnd = ",i0)
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Localize the node connectivity pointer for each cell
  !
  ! Create the local nodes_of_cell_ptr
  !
  if (allocated(nodes_of_cell_ptr)) then
    deallocate ( nodes_of_cell_ptr , stat=ierr , errmsg=error_message)
    call alloc_error(pname,"nodes_of_cell_ptr",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  allocate ( nodes_of_cell_ptr(1:lcell+lfbnd+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! First set the node connectivity pointer to the number
  ! of nodes that define each of the interior cells
  !
  do loc_cell = 1,lcell
    glb_cell = cell_map(this_part)%loc_to_glb(loc_cell)
    nodes_of_cell_ptr(loc_cell+1) = global_nofc_ptr(glb_cell+1) - &
                                    global_nofc_ptr(glb_cell)
  end do
  !
  ! Next set the node connectivity pointer to the number of nodes that
  ! define each of the ghost cells by using the face(:)%nodes array.
  !
  ! NOTE: bface(10,loc_face) could be negative to indicate this partition
  !       is not responsible for this boundary face
  !
  do loc_face = 1,lfbnd
    loc_ghost_cell = abs( bface(10,loc_face) )
    glb_face = face_map(this_part)%loc_to_glb(loc_face)
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
  ! ## 2 ##
  ! #######
  !
  ! Create the local nodes_of_cell array for each cell
  ! NOTE: The node indices copied into this array are actually the global
  !       indices and they will not be localized until the very end.
  !
  if (allocated(nodes_of_cell)) then
    deallocate ( nodes_of_cell , stat=ierr , errmsg=error_message)
    call alloc_error(pname,"nodes_of_cell",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  allocate ( nodes_of_cell(1:last(nodes_of_cell_ptr)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! First extract the information that is local to the current
  ! processor from the global_nofc array to the nodes_of_cell array.
  !
  do loc_cell = 1,lcell
    !
    l1 = nodes_of_cell_ptr(loc_cell)+1
    l2 = nodes_of_cell_ptr(loc_cell+1)
    !
    glb_cell = cell_map(this_part)%loc_to_glb(loc_cell)
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
    glb_face = face_map(this_part)%loc_to_glb(loc_face)
    !
    ! NOTE: bface(10,loc_face) could be negative to indicate this partition
    !       is not responsible for this boundary face
    !
    loc_ghost_cell = abs( bface(10,loc_face) )
    !
    l1 = nodes_of_cell_ptr(loc_ghost_cell)+1
    l2 = nodes_of_cell_ptr(loc_ghost_cell+1)
    !
    nfnods = geom_nodes(face(glb_face)%geom) ! # of nodes defining this face
    !
    nodes_of_cell(l1:l2) = face(glb_face)%nodes(1:nfnods)
    !
  end do
 !!
 !! Debugging localization of nodes_of_cell_ptr and nodes_of_cell by having
 !! the root processor output the global versions, and then having all
 !! processors output their local versions before the node indices that are
 !! stored in nodes_of_cell are reverse mapped to their local node indices.
 !!
 !if (mypnum == 0) then
 !  do loc_cell = 1,size(global_nofc_ptr)-1
 !    write (mypnum+130,*)
 !    l1 = global_nofc_ptr(loc_cell)+1
 !    l2 = global_nofc_ptr(loc_cell+1)
 !    write (mypnum+130,40) loc_cell,l1,l2
 !    write (mypnum+130,41) (global_nofc(l),l=l1,l2)
 !  end do
 !end if
 !!
 !flush (mypnum+130)
 !call mpi_barrier(MPI_COMM_WORLD,mpierr)
 !write (mypnum+120,'(a,i0)') "lcell = ",lcell
 !do loc_cell = 1,size(nodes_of_cell_ptr)-1
 !  write (mypnum+120,*)
 !  l1 = nodes_of_cell_ptr(loc_cell)+1
 !  l2 = nodes_of_cell_ptr(loc_cell+1)
 !  write (mypnum+120,40) loc_cell,l1,l2
 !  write (mypnum+120,41) (nodes_of_cell(l),l=l1,l2)
 !end do
 !!
 !call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Reverse map the interior nodes to the local ordering
  ! NOTE: The node indices stored in nodes of cells are currently the
  !       global node indices
  !
  allocate ( rev_node_map(1:gnode) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"rev_node_map",1,__LINE__,__FILE__,ierr,error_message)
  !
  rev_node_map(node_map(this_part)%loc_to_glb(:)) = &
                              intseq(1,size(node_map(this_part)%loc_to_glb))
  !
  do loc_node = 1,size(nodes_of_cell)
    glb_node = nodes_of_cell(loc_node)
    nodes_of_cell(loc_node) = rev_node_map(glb_node)
   !nodes_of_cell(loc_node) = findloc( node_map(this_part)%loc_to_glb , &
   !                                   glb_node )
  end do
  !
  deallocate ( rev_node_map , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"rev_node_map",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname,processor=0,part_num=this_part)
  !
  ! Format Statements
  !
 !40 format ("loc_cell=",i0,";   l1:l2=",i0,":",i0)
 !41 format ("nodes_of_cell(l1:l2) = ",*(i0,:,", "))
  !
end subroutine root_localizes_nodes_of_cell
!
!###############################################################################
!
subroutine root_localizes_xyz_nodes(this_part,node_map,global_xyz,xyz_nodes)
  !
  !.. Formal Arguments ..
  integer,                  intent(in) :: this_part
  type(map_t),              intent(in) :: node_map(:)
  real(wp),                 intent(in) :: global_xyz(:,:)
  real(wp), allocatable, intent(inout) :: xyz_nodes(:,:)
  !
  !.. Local Scalars ..
  integer :: ierr,lnode,loc_node,glb_node
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "root_localizes_xyz_nodes"
  !
continue
  !
  call debug_timer(entering_procedure,pname,processor=0,part_num=this_part)
  !
  ! Count the local number of nodes
  !
  lnode = size( node_map(this_part)%loc_to_glb )
  !
  ! Localize the node coordinates
  !
  if (allocated(xyz_nodes)) then
    deallocate ( xyz_nodes , stat=ierr , errmsg=error_message)
    call alloc_error(pname,"xyz_nodes",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( xyz_nodes(1:nr,1:lnode) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_nodes",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Localize the interior nodes
  !
  do loc_node = 1,lnode
    glb_node = node_map(this_part)%loc_to_glb(loc_node)
    xyz_nodes(:,loc_node) = global_xyz(:,glb_node)
  end do
  !
  call debug_timer(leaving_procedure,pname,processor=0,part_num=this_part)
  !
end subroutine root_localizes_xyz_nodes
!
!###############################################################################
!
subroutine root_localizes_grid_dt(this_part,cell_map,global_grid, &
                                  grid_elem,grid_xyz)
  !
  !.. Use Statements ..
  use geovar, only : grid_t
  !
  !.. Formal Arguments ..
  integer,                  intent(in) :: this_part
  type(map_t),              intent(in) :: cell_map(:)
  type(grid_t),             intent(in) :: global_grid
  integer,  allocatable, intent(inout) :: grid_elem(:)
  real(wp), allocatable, intent(inout) :: grid_xyz(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,n,ngp,nge,nlp,ierr
  integer :: loc_cell,glb_cell
  integer :: n1,n2,loc_pnt,glb_pnt
  !
  !.. Local Allocatable Arrays ..
  logical(lk), allocatable :: msk(:)
  integer,     allocatable :: l2g_map(:)
  integer,     allocatable :: g2l_map(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "root_localizes_grid_dt"
  !
continue
  !
  call debug_timer(entering_procedure,pname,processor=0,part_num=this_part)
  !
  ! Deallocate the grid_elem and grid_xyz arrays
  !
  if (allocated(grid_elem)) then
    deallocate ( grid_elem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid_elem",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(grid_xyz)) then
    deallocate ( grid_xyz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid_xyz",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Create the mask array to identify all the grid points for this partition
  !
  ngp = size(global_grid%xyz,dim=2)
  !
  allocate ( msk(1:ngp) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Find all the local cells for this partition and perform two tasks:
  !   1. Mark all the grid points using the msk array
  !   2. Find the size needed for grid_elem array
  !
  ! NOTE: RIGHT NOW, WE ARE ASSUMING THAT WE DONT NEED TO INCLUDE THE
  !       ELEMENTS IN global_grid%elem THAT ARE BOUNDARY FACES, SINCE THE FACE
  !       POINTS SHOULD JUST BE THE SAME AS THE POINTS ON THE FACE OF
  !       THE HOST CELL THAT CORRESPONDS TO THE BOUNDARY FACE.
  !
  call debug_timer(start_timer,"marking the grid points for this partition")
  nge = 0
  do loc_cell = 1,size(cell_map(this_part)%loc_to_glb)
    !
    glb_cell = cell_map(this_part)%loc_to_glb(loc_cell)
    !
    i = 1       ! 1 for elem_type, size(tags) combo
    if (send_grid_elem_host_cell) then
      i = i + 1 ! 1 for host cell
    end if
    !
    ! Add the size of global_grid%elem(:)%pts to the counter for this element
    !
    if (allocated(global_grid%elem(glb_cell)%pts)) then
      ! NOTE: Make sure that the counter is increased by at least 1 if the
      !       pts array is allocated but its size is 0.
      i = i + max( 1 , size(global_grid%elem(glb_cell)%pts) )
      msk( global_grid%elem(glb_cell)%pts ) = true
    else
      i = i + 1
    end if
    !
    ! Add the size of global_grid%elem(:)%tags to the counter for this element
    ! if we are sending this data to the other processors
    !
    if (send_grid_elem_tags) then
      if (allocated(global_grid%elem(glb_cell)%tags)) then
        i = i + size(global_grid%elem(glb_cell)%tags)
      else
        ! We will know if tags data exists for this element from the
        ! elem_type, size(tags) combo, so we dont need to add 1 here
      end if
    end if
    !
    ! Add the size of global_grid%elem(:)%nodes to the counter for this element
    !
    if (allocated(global_grid%elem(glb_cell)%nodes)) then
      ! NOTE: Make sure that the counter is increased by at least 1 if the
      !       nodes array is allocated but its size is 0.
      i = i + max( 1 , size(global_grid%elem(glb_cell)%nodes) )
      msk( global_grid%elem(glb_cell)%nodes ) = true
    else
      i = i + 1
    end if
    !
    nge = nge + i
    !
  end do
  call debug_timer(stop_timer,"marking the grid points for this partition")
  call debug_timer(start_timer, &
                   "creating grid_xyz and reverse point mapping array")
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Create grid_xyz
  !
  ! SLOWER TO INCLUDE THE PACK WITHIN THE AUTOMATIC REALLOCATION
 !grid_xyz = global_grid%xyz(:, pack(intseq(1,ngp),msk) ) ! F2003 auto-realloc
  !
  ! Create the points map to map the grid points local to this partition
  ! to their global point indices
  !
 !l2g_map = pack(intseq(1,ngp),msk) ! using F2003 automatic reallocation
  !
  nlp = count(msk)
  allocate ( l2g_map(1:nlp) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"l2g_map",1,__LINE__,__FILE__,ierr,error_message)
  n = 0
  do glb_pnt = 1,size(msk)
    if (msk(glb_pnt)) then
      n = n + 1
      l2g_map(n) = glb_pnt
    end if
  end do
  !
  ! Deallocate the msk array since it is no longer needed
  !
  deallocate ( msk , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate the grid_xyz array and copy data into it from global_grid%xyz
  !
 !grid_xyz = global_grid%xyz(:,l2g_map) ! F2003 auto-realloc
  !
  n1 = size(global_grid%xyz,dim=1)
  !
  allocate ( grid_xyz(1:n1,1:nlp), source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid_xyz",1,__LINE__,__FILE__,ierr,error_message)
  !
  do loc_pnt = 1,nlp
    grid_xyz(:,loc_pnt) = global_grid%xyz(:,l2g_map(loc_pnt))
  end do
  !
  ! Create the global point index to local point index reverse mapping array
  ! NOTE: I TRIED TO COMBINE THE PREVIOUS LOOP WITH THE FOLLOWING LOOP INTO A
  !       SINGLE LOOP, AS IN:
  !          do loc_pnt = 1,nlp
  !            glb_pnt = l2g_map(loc_pnt)
  !            grid_xyz(:,loc_pnt) = global_grid%xyz(:,glb_pnt)
  !            g2l_map(glb_pnt) = loc_pnt
  !          end do
  !       BUT THIS WAS ACTUALLY SLOWER THAN HAVING TWO SEPARATE LOOPS. IT MUST
  !       HAVE SOMETHING TO DO WITH ACCESSING ELEMENTS OF global_grid%xyz,
  !       grid_xyz, l2g_map, AND g2l_map FOR EACH LOOP ITERATION PREVENTS
  !       THE OPTIMIZER FROM TAKING ADVANTAGE OF DATA LOCALITY TO PROPERLY
  !       VECTORIZE THESE TASKS.
  !
  ! SLOWER TO USE UNPACK AND AUTOMATIC REALLOCATION
 !g2l_map = unpack(intseq(1,count(msk)),mask=msk,field=0) ! F2003 auto-realloc
  !
  allocate ( g2l_map(1:ngp) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"g2l_map",1,__LINE__,__FILE__,ierr,error_message)
  !
  do loc_pnt = 1,nlp
    g2l_map(l2g_map(loc_pnt)) = loc_pnt
  end do
  !
  !
  ! Deallocate l2g_map since it is no longer needed
  !
  if (allocated(l2g_map)) then
    deallocate ( l2g_map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"l2g_map",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Create grid_elem
  !
  ! Allocate the grid_elem array
  !
  call debug_timer(stop_timer, &
                   "creating grid_xyz and reverse point mapping array")
  call debug_timer(start_timer,"creating grid_elem")
  allocate ( grid_elem(1:nge) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid_elem",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Loop over the elements for this partition and store the data for
  ! each element in the grid_elem array
  !
  n = 0
  do loc_cell = 1,size(cell_map(this_part)%loc_to_glb)
    !
    glb_cell = cell_map(this_part)%loc_to_glb(loc_cell)
    !
    ! Since the values of global_grid%elem(glb_cell)%prop%elem_type and
    ! size(global_grid%elem(glb_cell)%tags) should both be at most 16-bit
    ! integers, combine them into a single default sized integer (probably
    ! 32-bit) to reduce the size needed for the grid_elem array.
    !
    ! FOR ELEM_TYPE, SIZE(TAGS) COMBO:
    ! combo = 1000*size(tags) + elem_type
    ! elem_type = mod(combo,1000)
    ! ntags = combo/1000
    !
    i = 1 ! 1 for elem_type, size(tags) combo
    j = global_grid%elem(glb_cell)%prop%elem_type
    if (allocated(global_grid%elem(glb_cell)%tags)) then
      j = j + merge(1000,0,send_grid_elem_tags) * &
              size(global_grid%elem(glb_cell)%tags)
    end if
    grid_elem(n+i) = j
    !
    ! Add global_grid%elem(:)%host_cell if needed
    ! NOTE: WE NEED TO LOCALIZE HOST CELL BEFORE ADDING IT TO grid_elem
    !
    if (send_grid_elem_host_cell) then
      i = i + 1 ! 1 for host cell
      grid_elem(n+i) = findloc( cell_map(this_part)%loc_to_glb , &
                                global_grid%elem(glb_cell)%host_cell )
    end if
    !
    ! Add the localized indices for global_grid%elem(:)%pts
    !
    ! Start by setting the next value of grid_elem to 0 in case the pts
    ! component isnt allocated or its allocated with a size of 0. This will
    ! indicate to the receiving processor that the global_grid%elem(:)%pts
    ! array didnt contain any data for this element.
    !
    i = i + 1
    grid_elem(n+i) = 0
    !
    ! Check to make sure global_grid%elem(glb_cell)%pts is
    ! allocated before trying to add it to grid_elem
    !
    if (allocated(global_grid%elem(glb_cell)%pts)) then
      !
      ! The global_grid%elem(:)%pts array is allocated, now make sure
      ! its size is greater than zero and contains data
      ! NOTE: We need to perform this check so that we can decrease the
      !       counter for this element before we start adding stuff to
      !       grid_elem.
      !
      if (size(global_grid%elem(glb_cell)%pts) > 0) then
        !
        ! Its size is greater than 0, so first decrease the counter for this
        ! element by 1, so the previous assignment of 0 to what was the next
        ! value of grid_elem can be overwritten
        !
        i = i - 1
        !
        ! Now add the localized indices from
        ! global_grid%elem(:)%pts to grid_elem
        !
        do j = lbound(global_grid%elem(glb_cell)%pts,dim=1), &
               ubound(global_grid%elem(glb_cell)%pts,dim=1)
          i = i + 1
          grid_elem(n+i) = g2l_map( global_grid%elem(glb_cell)%pts(j) )
        end do
        !
      end if
      !
    end if
    !
    ! Add global_grid%elem(:)%tags if its allocated and we
    ! are sending this data to the other processors
    !
    if (send_grid_elem_tags) then
      if (allocated(global_grid%elem(glb_cell)%tags)) then
        do j = lbound(global_grid%elem(glb_cell)%tags,dim=1), &
               ubound(global_grid%elem(glb_cell)%tags,dim=1)
          i = i + 1
          grid_elem(n+i) = global_grid%elem(glb_cell)%tags(j)
        end do
      end if
    end if
    !
    ! Finally, add the localized indices for global_grid%elem(:)%nodes
    !
    ! Start by setting the next value of grid_elem to 0 in case the nodes
    ! component isnt allocated or its allocated with a size of 0. This will
    ! indicate to the receiving processor that the global_grid%elem(:)%nodes
    ! array didnt contain any data for this element.
    !
    i = i + 1
    grid_elem(n+i) = 0
    !
    ! Check to make sure global_grid%elem(glb_cell)%nodes is
    ! allocated before trying to add it to grid_elem
    !
    if (allocated(global_grid%elem(glb_cell)%nodes)) then
      !
      ! The global_grid%elem(:)%nodes array is allocated, now make sure
      ! its size is greater than zero and contains data
      ! NOTE: We need to perform this check so that we can decrease the
      !       counter for this element before we start adding stuff to
      !       grid_elem.
      !
      if (size(global_grid%elem(glb_cell)%nodes) > 0) then
        !
        ! Its size is greater than 0, so first decrease the counter for this
        ! element by 1, so the previous assignment of 0 to what was the next
        ! value of grid_elem can be overwritten
        !
        i = i - 1
        !
        ! Now add the localized indices from
        ! global_grid%elem(:)%nodes to grid_elem
        !
        do j = lbound(global_grid%elem(glb_cell)%nodes,dim=1), &
               ubound(global_grid%elem(glb_cell)%nodes,dim=1)
          i = i + 1
          grid_elem(n+i) = g2l_map( global_grid%elem(glb_cell)%nodes(j) )
        end do
        !
      end if
      !
    end if
    !
    ! Add the final counter for this element to the grid_elem
    ! counter before going to the next element
    !
    n = n + i
    !
  end do
  call debug_timer(stop_timer,"creating grid_elem")
  !
  ! Check to make sure the final counter matches the size of the grid_elem array
  !
  if (n /= size(grid_elem)) then
    write (error_message,1) this_part,n,size(grid_elem)
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Deallocate g2l_map before we leave
  !
  if (allocated(g2l_map)) then
    deallocate ( g2l_map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"g2l_map",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname,processor=0,part_num=this_part)
  !
  ! Format Statements
  !
  1 format ("ERROR CREATING grid_elem FOR PARTITION #",i0,"!!! ", &
            "THE FINAL COUNTER n=",i0,", WHEN IT SHOULD BE ", &
            "size(grid_elem)=",i0,"!!!")
  !
end subroutine root_localizes_grid_dt
!
!###############################################################################
!
subroutine localize_bc_conflict(global_bface,face_map,part_id)
  !
  !.. Use Statements ..
  use geovar, only : bc_conflict
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: global_bface(:,:)
  type(map_t), intent(in) :: face_map(:)
  !
  !.. Optional Arguments
  integer, optional, intent(in) :: part_id
  !
  !.. Local Scalars ..
  integer :: n,master_cpu,glb_face
  integer :: nf,slave_cpu,loc_face
  integer :: my_part,adj_face
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "localize_bc_conflict"
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
    call debug_timer(entering_procedure,pname)
    !
    my_part = mypnum + 1
    if (present(part_id)) my_part = part_id
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
      ! Get the index for the matching face on the adjacent partition
      !
      adj_face = global_bface(5,glb_face)
      !
      ! For the adjacent face, get the adjacent partition processor ID
      ! which is just the processor ID for the current glb_face
      ! NOTE: The value in global_bface(7,:) is the processor ID, but we
      !       want the partition number which is just the processor ID + 1
      !
      master_cpu = global_bface(7,adj_face) + 1
      ! Equivalent to above, but using epart
     !master_cpu = epart( global_bface(2,glb_face) )
      !
      loc_face = findloc( face_map(master_cpu)%loc_to_glb , glb_face)
      !
      bc_conflict(n)%master%face = loc_face
      bc_conflict(n)%master%cpu = master_cpu
      !
      ! Mark this conflict if this processor is involved with it
      !
      if (master_cpu == my_part) bc_conflict(n)%not_my_conflict = fals
      !
      ! Now do the same for the slave faces
      !
      do nf = 1,size(bc_conflict(n)%slave)
        !
        glb_face = bc_conflict(n)%slave(nf)%face
        !
        ! Get the index for the matching face on the adjacent partition
        !
        adj_face = global_bface(5,glb_face)
        !
        ! For the adjacent face, get the adjacent partition processor ID
        ! which is just the processor ID for the current glb_face
        ! NOTE: The value in global_bface(7,:) is the processor ID, but we
        !       want the partition number which is just the processor ID + 1
        !
        slave_cpu = global_bface(7,adj_face) + 1
        ! Equivalent to above, but using epart
       !slave_cpu = epart( global_bface(2,glb_face) )
        !
        loc_face = findloc( face_map(slave_cpu)%loc_to_glb , glb_face)
        !
        bc_conflict(n)%slave(nf)%face = loc_face
        bc_conflict(n)%slave(nf)%cpu = slave_cpu
        !
        ! Mark this conflict if this processor is involved with it
        !
        if (slave_cpu == my_part) bc_conflict(n)%not_my_conflict = fals
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
   !          trim(bc_integer_to_string( &
   !                   global_bface(1,bc_conflict(n)%master%face) &
   !                                   ))
   !  do nf = 1,size(bc_conflict(n)%slave)
   !    write (mypnum+200,2) nf,bc_conflict(n)%slave(nf)%face, &
   !                            bc_conflict(n)%slave(nf)%flux_point, &
   !                            bc_conflict(n)%slave(nf)%cpu, &
   !          trim(bc_integer_to_string( &
   !                   global_bface(1,bc_conflict(n)%slave(nf)%face) &
   !                                   ))
   !  end do
   !end do
    !
    call debug_timer(leaving_procedure,pname)
    !
  end if
  !
  ! Format Statements
  !
 !1 format ("Conflict #",i0," - is_local=",l1,", not_my_conflict=",l1,/, &
 !          "MF :   Face=",i5,",  FluxPoint=",i2,",  CPU=",i0, &
 !          "    BC=",a)
 !2 format ("SF",i0,":   Face=",i5,",  FluxPoint=",i2,",  CPU=",i0, &
 !          "    BC=",a)
  !
end subroutine localize_bc_conflict
!
!###############################################################################
!
subroutine create_global_collective_datatypes(cell_map,global_geom,global_order)
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  use ovar,      only : use_old_restart_format
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: global_geom(:)
  integer,     intent(in) :: global_order(:)
  type(map_t), intent(in) :: cell_map(:)
  !
  !.. Local Scalars ..
  integer :: loc_cell,gcell,ierr,np,this_part
  integer :: glb_cell,glb_geom,glb_order
  integer(int_mpi) :: np_int_mpi
  !
  !.. Local Allocatable Arrays ..
  integer,          allocatable :: cell_offset(:)
  integer,          allocatable :: recv_offset(:)
  integer(int_mpi), allocatable :: recv_offset_int_mpi(:)
  integer,          allocatable :: recv_blklen(:)
  integer(int_mpi), allocatable :: recv_blklen_int_mpi(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_global_collective_datatypes"
  !
continue
  !
  ! Only create the collective datatypes if the old restart format is being used
  !
  if (use_old_restart_format) then
    !
    ! We need to create collective receiving MPI datatypes to collect
    ! localized data from each processor into a shared global form.
    !
    ! Before we start, find the total number of grid cells
    !
    gcell = sum_map_sizes(cell_map)
    !
    ! First, we need to create an array containing the
    ! solution point offsets for every global cell
    !
    allocate ( cell_offset(1:gcell) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_offset",1,__LINE__,__FILE__,ierr,error_message)
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
    ! Find the maximum number of cells for a single partition
    np = 0
    do this_part = 1,size(cell_map)
      np = max( np , size(cell_map(this_part)%loc_to_glb) )
    end do
    !
    ! Allocate the local arrays needed for storing the offsets and block lengths
    ! required as input into the MPI procedures for creating MPI datatypes
    !
    allocate ( recv_offset(1:np) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_offset",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( recv_offset_int_mpi(1:np) , source=0_int_mpi , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_offset_int_mpi",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    allocate ( recv_blklen(1:np) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_blklen",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( recv_blklen_int_mpi(1:np) , source=0_int_mpi , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_blklen_int_mpi",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    ! Allocate the module arrays that will contained
    ! the collective receive MPI datatypes
    !
   !allocate ( recv_xyz(1:ncpu) , source=MPI_DATATYPE_NULL , &
   !           stat=ierr , errmsg=error_message )
   !call alloc_error(pname,"recv_xyz",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( recv_sol(1:ncpu) , source=MPI_DATATYPE_NULL , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_sol",1,__LINE__,__FILE__,ierr,error_message)
    !
   !allocate ( recv_var(1:ncpu) , source=MPI_DATATYPE_NULL , &
   !           stat=ierr , errmsg=error_message )
   !call alloc_error(pname,"recv_var",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Loop over all the processors to create
    ! the collective receive MPI datatypes
    !
    do this_part = 1,ncpu
      !
      np = size(cell_map(this_part)%loc_to_glb)
      np_int_mpi = int( np , kind=int_mpi )
      !
      ! Get the offset and block length for each cell of this partition
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
     !!
     !! Now create the receive MPI datatype for the local
     !! solution point coordinates on the current processor
     !!
     !! adjust the offset for nr dimensions
     !recv_offset(1:np) = recv_offset(1:np)*nr
     !recv_blklen(1:np) = recv_blklen(1:np)*nr
     !!
     !recv_offset_int_mpi(1:np) = int( recv_offset(1:np) , kind=int_mpi )
     !recv_blklen_int_mpi(1:np) = int( recv_blklen(1:np) , kind=int_mpi )
     !!
     !call mpi_type_indexed(np_int_mpi, &
     !                      recv_blklen_int_mpi(1:np), &
     !                      recv_offset_int_mpi(1:np), &
     !                      mpi_flttyp,recv_xyz(this_part),mpierr)
     !call mpi_type_commit(recv_xyz(this_part),mpierr)
     !!
     !! cancel out the nr adjustment above
     !recv_offset(1:np) = recv_offset(1:np)/nr
     !recv_blklen(1:np) = recv_blklen(1:np)/nr
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
     !! Now create the receive MPI datatype for an arbitrary
     !! variable at the solution points on the current processor
     !!
     !recv_offset_int_mpi(1:np) = int( recv_offset(1:np) , kind=int_mpi )
     !recv_blklen_int_mpi(1:np) = int( recv_blklen(1:np) , kind=int_mpi )
     !!
     !call mpi_type_indexed(np_int_mpi, &
     !                      recv_blklen_int_mpi(1:np), &
     !                      recv_offset_int_mpi(1:np), &
     !                      mpi_flttyp,recv_var(this_part),mpierr)
     !call mpi_type_commit(recv_var(this_part),mpierr)
     !!
    end do
    !
    ! Deallocate the local arrays before leaving
    !
    deallocate ( cell_offset , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_offset",2,__LINE__,__FILE__,ierr,error_message)
    !
    deallocate ( recv_offset , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_offset",2,__LINE__,__FILE__,ierr,error_message)
    !
    deallocate ( recv_offset_int_mpi , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_offset_int_mpi",2,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    deallocate ( recv_blklen , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_blklen",2,__LINE__,__FILE__,ierr,error_message)
    !
    deallocate ( recv_blklen_int_mpi , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_blklen_int_mpi",2,__LINE__,__FILE__,ierr, &
                     error_message)
    !
  end if
  !
end subroutine create_global_collective_datatypes
!
!###############################################################################
!
subroutine create_face_nonblocking_datatypes(bface,face_map,flx,face_rotation)
!subroutine create_face_nonblocking_datatypes(bface,face,face_map, &
!                                             flx,face_rotation)
  !
  !.. Formal Arguments ..
  integer,       intent(inout) :: bface(:,:)
 !type(face_t),     intent(in) :: face(:)
  type(map_t),      intent(in) :: face_map(:)
  type(flux_pts_t), intent(in) :: flx(:)
  integer,          intent(in) :: face_rotation(:,:)
  !
  !.. Local Scalars ..
  logical(lk) :: flag
  integer :: n,nfp,lfbnd,i,ibeg,iend,interval,ierr
  integer :: rmin,rmax,scnt,rcnt,p1,p2
  integer :: loc_face,glb_face,my_part
  integer :: adj_face,adj_cpu
  integer :: left_rotation
  integer(int_mpi) :: block_size,np
  integer(int_mpi) :: this_flxpt
 !integer(int_mpi) :: tag_ub
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
  logical(lk),      allocatable :: send_msk(:)
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
  call debug_timer(entering_procedure,pname)
  !
  my_part = mypnum + 1
  !
  lfbnd = size(bface,dim=2)
  !
  ! Count the number of faces located on partition boundaries
  !
  n = count(abs(bface(1,:)) == bc_cpu_bnd)
  !
  ! Allocate the bnd_faces array
  !
  allocate ( bnd_faces(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bnd_faces",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Now fill in the components within the bnd_faces array
  !
  n = 0
  fill_bnd_faces: do loc_face = 1,lfbnd
    !
    if (abs(bface(1,loc_face)) /= bc_cpu_bnd) cycle fill_bnd_faces
    !
    n = n + 1
    !
    bnd_faces(n)%idx = loc_face
    !
    glb_face = face_map(my_part)%loc_to_glb(loc_face)
    !
    if (bface(10,loc_face) > 0) then
      bnd_faces(n)%my_responsibility = true
    end if
    !
  end do fill_bnd_faces
  !
  ! We no longer need bface(10,:) to be positive/negative to differentiate
  ! between partitions that are responsible for the boundary face, so make
  ! bface(10,:) positive for all boundary faces
  !
  bface(10,:) = abs(bface(10,:))
  !
  ! Create an offset array for the number of flux
  ! points on each communication boundary faces
  !
  allocate ( flxpt_offset(1:lfbnd) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"flxpt_offset",1,__LINE__,__FILE__,ierr,error_message)
  !
  do loc_face = 1,lfbnd
    !
    if (abs(bface(1,loc_face)) /= bc_cpu_bnd) cycle
    !
    glb_face = face_map(my_part)%loc_to_glb(loc_face)
    !
    left_rotation = face_rotation(1,glb_face)
    !
    flxpt_offset(loc_face) = size( flx(left_rotation)%pts , kind=int_mpi )
   !flxpt_offset(loc_face) = size( face(glb_face)%left%flxpts , kind=int_mpi )
    !
  end do
  !
  ! Reshuffle flxpt_offset by summing the total number of flux points for
  ! all the communication boundary faces before the current one.
  !
  do loc_face = lfbnd,1,-1
    !
    if (abs(bface(1,loc_face)) /= bc_cpu_bnd) cycle
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
    if (abs(bface(1,loc_face)) == bc_cpu_bnd) then
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
  call alloc_error(pname,"exch_fusp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( exch_fdusp(1:ncomm) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"exch_fusp",1,__LINE__,__FILE__,ierr,error_message)
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
  call alloc_error(pname,"send_msk",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( recv_map(rmin:rmax) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_map",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( comm_pts(1:ncomm) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"comm_pts",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Loop through each of the processes this one will be communicating with
  ! and create the send and receive face maps
  !
  do i = 1,ncomm
    !
    adj_cpu = int( exch_fusp(i)%adj_cpu , kind=kind(adj_cpu) )
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
      if (abs(bface(1,loc_face)) /= bc_cpu_bnd) cycle
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
      left_rotation = face_rotation(1,glb_face)
      !
      scnt = scnt + size( flx(left_rotation)%pts )
     !scnt = scnt + size( face(glb_face)%left%flxpts )
      !
    end do
    !
    ! Allocate the send and recv components of comm_pts now that we know
    ! the total number of flux points that will need to be communicated.
    !
    allocate ( comm_pts(i)%send(1:scnt) , source=0_int_mpi , &
               stat=ierr , errmsg=error_message )
    write (array_name,100) i,"send"
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( comm_pts(i)%recv(1:scnt) , source=0_int_mpi , &
               stat=ierr , errmsg=error_message )
    write (array_name,100) i,"recv"
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
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
      left_rotation = face_rotation(1,glb_face)
      !
      ibeg = lbound( flx(left_rotation)%pts , dim=1 )
      iend = ubound( flx(left_rotation)%pts , dim=1 )
      interval = 1
      !
      ! Reverse the flux point ordering if this
      ! used to be a periodic boundary face.
      !
      if (bface(1,loc_face) == -bc_cpu_bnd) then
        if (nr < 3) then
          ibeg = ubound( flx(left_rotation)%pts , dim=1 )
          iend = lbound( flx(left_rotation)%pts , dim=1 )
          interval = -1
        end if
        bface(1,loc_face) = bc_cpu_bnd
      end if
      !
     !do n = 1,size( face(glb_face)%left%flxpts )
     !do n = 1,size( flx(left_rotation)%pts )
      do n = ibeg,iend,interval
        !
        scnt = scnt + 1
        !
        this_flxpt = int( flx(left_rotation)%pts(n) , kind=int_mpi )
       !this_flxpt = int( face(glb_face)%left%flxpts(n) , kind=int_mpi )
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
      left_rotation = face_rotation(1,glb_face)
      !
      do n = lbound( flx(left_rotation)%pts , dim=1 ), &
             ubound( flx(left_rotation)%pts , dim=1 )
     !do n = 1,size( face(glb_face)%left%flxpts )
        !
        rcnt = rcnt + 1
        !
        this_flxpt = int( flx(left_rotation)%pts(n) , kind=int_mpi )
       !this_flxpt = int( face(glb_face)%left%flxpts(n) , kind=int_mpi )
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
  call alloc_error(pname,"send_msk",2,__LINE__,__FILE__,ierr,error_message)
  deallocate ( recv_map , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_map",2,__LINE__,__FILE__,ierr,error_message)
  deallocate ( flxpt_offset , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"flxpt_offset",2,__LINE__,__FILE__,ierr,error_message)
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
    comm_pts(i)%send(:) = comm_pts(i)%send(:) * block_size
    comm_pts(i)%recv(:) = comm_pts(i)%recv(:) * block_size
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
    comm_pts(i)%send(:) = comm_pts(i)%send(:) / block_size
    comm_pts(i)%recv(:) = comm_pts(i)%recv(:) / block_size
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
    comm_pts(i)%send(:) = comm_pts(i)%send(:) * block_size
    comm_pts(i)%recv(:) = comm_pts(i)%recv(:) * block_size
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
    comm_pts(i)%send(:) = comm_pts(i)%send(:) / block_size
    comm_pts(i)%recv(:) = comm_pts(i)%recv(:) / block_size
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
  ! We no longer need comm_pts so deallocate
  !
  deallocate ( comm_pts , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"comm_pts",1,__LINE__,__FILE__,ierr,error_message)
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
 !call mpi_barrier(MPI_COMM_WORLD,mpierr)
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
 !call mpi_barrier(MPI_COMM_WORLD,mpierr)
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
    left_rotation = face_rotation(1,glb_face)
    nfp = nfp + size( flx(left_rotation)%pts )
   !nfp = nfp + size( face(glb_face)%left%flxpts )
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
  call alloc_error(pname,"fusp_recv_buffer",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! Send buffer for communicating solution variables
  !
  allocate ( fusp_send_buffer(1:nq*nfp) , source=0.0_r8 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fusp_send_buffer",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! MPI request objects for communicating solution variables
  !
  allocate ( fusp_rqst(1:2*ncomm) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fusp_rqst",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! MPI status objects returned after communicating solution variables
  !
  allocate ( fusp_istat(1:_MPI_STATUS_SIZE_,1:2*ncomm) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fusp_istat",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Arrays for exchanging gradients of solution variables at face flux points
  !
  ! Receive buffer for communicating gradients
  !
  allocate ( fdusp_recv_buffer(1:nr*nq*nfp) , source=0.0_r8 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fdusp_recv_buffer",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! Send buffer for communicating gradients
  !
  allocate ( fdusp_send_buffer(1:nr*nq*nfp) , source=0.0_r8 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fdusp_send_buffer",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! MPI request objects for communicating gradients
  !
  allocate ( fdusp_rqst(1:2*ncomm) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fdusp_rqst",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! MPI status objects returned after communicating gradients
  !
  allocate ( fdusp_istat(1:_MPI_STATUS_SIZE_,1:2*ncomm) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fdusp_istat",1,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," CPU ",i0," -> CPU ",i0," :: Sending")
  2 format (/," CPU ",i0," <- CPU ",i0," :: Receiving")
  3 format (100(12i6,/))
 !3 format (*(12i6,/))
  4 format (" Cell = ",i0,100(/,12i6))
 !4 format (" Cell = ",i0,*(/,12i6))
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
subroutine bcast_node_and_edge_arrays(node,edge)
  !
  !.. Formal Arguments ..
  type(node_t), allocatable, intent(inout) :: node(:)
  type(edge_t), allocatable, intent(inout) :: edge(:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,n,nn,ne,nep,npd,num_cpus,ierr
  !
  character(len=100) :: array_name
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: packed_data(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "bcast_node_and_edge_arrays"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (mypnum == 0) then
    !
    nn = 0
    if (allocated(node)) then
      do i = 1,size(node)
        nn = nn + 3
        if (allocated(node(i)%cpu)) then
          nn = nn + 3*size(node(i)%cpu)
        end if
      end do
    end if
    !
    ne = 0
    if (allocated(edge)) then
      do i = 1,size(edge)
        ne = ne + 4
        if (allocated(edge(i)%cpu)) then
          do j = 1,size(edge(i)%cpu)
            ne = ne + 4
            if (allocated(edge(i)%cpu(j)%edgpts)) then
              ne = ne + size(edge(i)%cpu(j)%edgpts)
            end if
          end do
        end if
      end do
    end if
    !
    npd = nn + ne + 2
    !
  else
    !
    if (allocated(node)) then
      deallocate ( node , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"node",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    if (allocated(edge)) then
      deallocate ( edge , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"edge",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
  end if
  !
  call mpi_bcast(npd,1_int_mpi,mpi_inttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
  !
  ! Leave this subroutine if npd == 2, which means that both the
  ! node and edge arrays are not allocated on the root process
  !
  if (npd == 2) then
    call debug_timer(leaving_procedure,pname)
    return
  end if
  !
  allocate ( packed_data(1:npd) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"packed_data",1,__LINE__,__FILE__,ierr,error_message)
  !
  if (mypnum == 0) then
    !
    n = 0
    !
    if (allocated(node)) then
      n = n+1 ; packed_data(n) = size(node)
      do i = 1,size(node)
        n = n+1 ; packed_data(n) = merge(1,0,node(i)%is_on_boundary)
        n = n+1 ; packed_data(n) = node(i)%controlling_partition
        if (allocated(node(i)%cpu)) then
          n = n+1 ; packed_data(n) = size(node(i)%cpu)
          do j = 1,size(node(i)%cpu)
            n = n+1 ; packed_data(n) = node(i)%cpu(j)%partition
            n = n+1 ; packed_data(n) = node(i)%cpu(j)%loc_node
            n = n+1 ; packed_data(n) = node(i)%cpu(j)%num_cells
          end do
        else
          n = n+1 ; packed_data(n) = -1
        end if
      end do
    else
      n = n+1 ; packed_data(n) = -1
    end if
    !
    if (allocated(edge)) then
      n = n+1 ; packed_data(n) = size(edge)
      do i = 1,size(edge)
        n = n+1 ; packed_data(n) = merge(1,0,edge(i)%is_on_boundary)
        n = n+1 ; packed_data(n) = edge(i)%controlling_partition
        n = n+1 ; packed_data(n) = edge(i)%order
        if (allocated(edge(i)%cpu)) then
          n = n+1 ; packed_data(n) = size(edge(i)%cpu)
          do j = 1,size(edge(i)%cpu)
            n = n+1 ; packed_data(n) = edge(i)%cpu(j)%partition
            n = n+1 ; packed_data(n) = edge(i)%cpu(j)%loc_edge
            n = n+1 ; packed_data(n) = edge(i)%cpu(j)%num_cells
            if (allocated(edge(i)%cpu(j)%edgpts)) then
              n = n+1 ; packed_data(n) = size(edge(i)%cpu(j)%edgpts)
              do k = 1,size(edge(i)%cpu(j)%edgpts)
                n = n+1 ; packed_data(n) = edge(i)%cpu(j)%edgpts(k)
              end do
            else
              n = n+1 ; packed_data(n) = -1
            end if
          end do
        else
          n = n+1 ; packed_data(n) = -1
        end if
      end do
    else
      n = n+1 ; packed_data(n) = -1
    end if
    !
  end if
  !
  call mpi_bcast(packed_data,size(packed_data,kind=int_mpi),mpi_inttyp, &
                 0_int_mpi,MPI_COMM_WORLD,mpierr)
  !
  if (mypnum /= 0) then
    !
    n = 0
    !
    n = n+1 ; nn = packed_data(n)
    !
    if (nn > 0) then
      !
      allocate ( node(1:nn) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"node",1,__LINE__,__FILE__,ierr,error_message)
      !
      do i = 1,size(node)
        !
        n = n+1 ; node(i)%is_on_boundary = (packed_data(n) == 1)
        n = n+1 ; node(i)%controlling_partition = packed_data(n)
        n = n+1 ; num_cpus = packed_data(n)
        !
        if (num_cpus >= 0) then
          !
          allocate ( node(i)%cpu(1:num_cpus) ,  &
                     stat=ierr , errmsg=error_message )
          write (array_name,1) "node",i
          call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                           error_message)
          !
          do j = 1,size(node(i)%cpu)
            n = n+1 ; node(i)%cpu(j)%partition = packed_data(n)
            n = n+1 ; node(i)%cpu(j)%loc_node = packed_data(n)
            n = n+1 ; node(i)%cpu(j)%num_cells = packed_data(n)
          end do
          !
        end if
      end do
      !
    end if
    !
    n = n+1 ; ne = packed_data(n)
    !
    if (ne > 0) then
      !
      allocate ( edge(1:ne) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"edge",1,__LINE__,__FILE__,ierr,error_message)
      !
      do i = 1,size(edge)
        !
        n = n+1 ; edge(i)%is_on_boundary = (packed_data(n) == 1)
        n = n+1 ; edge(i)%controlling_partition = packed_data(n)
        n = n+1 ; edge(i)%order = packed_data(n)
        n = n+1 ; num_cpus = packed_data(n)
        !
        if (num_cpus >= 0) then
          !
          allocate ( edge(i)%cpu(1:num_cpus) ,  &
                     stat=ierr , errmsg=error_message )
          write (array_name,1) "edge",i
          call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                           error_message)
          !
          do j = 1,size(edge(i)%cpu)
            !
            n = n+1 ; edge(i)%cpu(j)%partition = packed_data(n)
            n = n+1 ; edge(i)%cpu(j)%loc_edge = packed_data(n)
            n = n+1 ; edge(i)%cpu(j)%num_cells = packed_data(n)
            n = n+1 ; nep = packed_data(n)
            !
            if (nep >= 0) then
              !
              allocate ( edge(i)%cpu(j)%edgpts(1:nep) , source=0 , &
                         stat=ierr , errmsg=error_message )
              write (array_name,1) "edge",i,j
              call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                               error_message)
              !
              do k = 1,size(edge(i)%cpu(j)%edgpts)
                n = n+1 ; edge(i)%cpu(j)%edgpts(k) = packed_data(n)
              end do
              !
            end if
            !
          end do
          !
        end if
      end do
      !
    end if
    !
  end if
  !
  if (allocated(packed_data)) then
    deallocate ( packed_data , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"packed_data",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (a,"(",i0,")%cpu",:,"(",i0,")%edgpts")
  !
end subroutine bcast_node_and_edge_arrays
!
!###############################################################################
!
subroutine create_node_and_edge_datatypes(node,node_map,edge,edge_map)
  !
  !.. Formal Arguments ..
  type(node_t), allocatable, intent(in) :: node(:)
  type(map_t),               intent(in) :: node_map(:)
  type(edge_t), allocatable, intent(in) :: edge(:)
  type(map_t),               intent(in) :: edge_map(:)
  !
  !.. Local Scalars ..
  integer :: i,ierr,n,n1,n2,nep,nvar,scnt,rcnt,rmin,rmax
  integer :: my_part,adj_part,npair,this_pair,last_offset
  integer :: lnode,loc_node,glb_node,adj_node,node_part
  integer :: ledge,loc_edge,glb_edge,adj_edge,edge_part
  integer(int_mpi) :: adj_cpu
  integer(int_mpi) :: this_edgpt,tag12,np
  integer(int_mpi) :: block_size,tag21
  integer(int_mpi) :: max_nodpt_offset
  real(wp) :: rd,rn
  !
  !.. Local Arrays ..
  integer     :: rnmin(1:ncpu)
  integer     :: rnmax(1:ncpu)
  integer     :: remin(1:ncpu)
  integer     :: remax(1:ncpu)
  logical(lk) :: cpu_msk(1:ncpu)
  !
  integer(int_mpi) :: cpu_nums(1:ncpu)
  integer(int_mpi) :: cpu_disp(1:ncpu)
  !
  !.. Local Allocatable Arrays ..
  logical(lk),      allocatable :: send_node_msk(:)
  integer,          allocatable :: recv_node_map(:)
  logical(lk),      allocatable :: send_edge_msk(:)
  integer,          allocatable :: recv_edge_map(:)
  integer,          allocatable :: comm_pair(:,:)
  integer(int_mpi), allocatable :: send_disp(:)
  integer(int_mpi), allocatable :: recv_disp(:)
  integer(int_mpi), allocatable :: edgpt_offset(:)
  integer(int_mpi), allocatable :: nodpt_offset(:)
  integer(int_mpi), allocatable :: nodpt_offset_save(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_node_and_edge_datatypes"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  call debug_timer(start_timer,"section 1")
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
    !
    ! Check to make sure node%cpu is allocated for this node
    !
    if (allocated(node(glb_node)%cpu)) then
      !
      ! Check to see if this node is on multiple partitions which
      ! identifies that this node is on a communication boundary
      !
      if (size(node(glb_node)%cpu) > 1) then
        !
        ! Check to make sure that the current partition is one of the
        ! partitions that needs to communicate data for this node.
        !
        ! NOTE: This check prevents the situation where
        !       1. This partition contains this node but it is only found in
        !          the definition of an interior cell
        !       and
        !       2. This node is found in the definition of boundary faces on
        !          two or more additional partitions
        !
        !                          \|
        ! For example:     BOUNDARY\| Partition B
        !               ///////////\|
        !               ------------*---------------
        !                           |
        !               Partition A | This partition
        !                           |
        !
        if (any(node(glb_node)%cpu%partition == my_part)) then
          n = n + 1
        end if
      end if
    end if
  end do
  !
  ! Allocate the bnd_nodes array
  !
  allocate ( bnd_nodes(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bnd_nodes",1,__LINE__,__FILE__,ierr,error_message)
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
        do i = 1,size(node(glb_node)%cpu)
          if (my_part == node(glb_node)%cpu(i)%partition) then
            !
            n = n + 1
            !
            bnd_nodes(n)%idx = loc_node
            !
            rd = real( sum( node(glb_node)%cpu(:)%num_cells ) , kind=wp )
            rn = real( node(glb_node)%cpu(i)%num_cells , kind=wp )
            !
            bnd_nodes(n)%ave_correction = rn / rd
            bnd_nodes(n)%my_responsibility = &
                   ( node(glb_node)%controlling_partition == my_part )
            !
            cycle fill_bnd_nodes
            !
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
        if (any(edge(glb_edge)%cpu%partition == my_part)) then
          n = n + 1
        end if
      end if
    end if
  end do
  !
  ! Allocate the bnd_edges array
  !
  allocate ( bnd_edges(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bnd_edges",1,__LINE__,__FILE__,ierr,error_message)
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
        do i = 1,size(edge(glb_edge)%cpu)
          if (my_part == edge(glb_edge)%cpu(i)%partition) then
            !
            n = n + 1
            !
            bnd_edges(n)%idx = loc_edge
            !
            rd = real( sum( edge(glb_edge)%cpu(:)%num_cells ) , kind=wp )
            rn = real( edge(glb_edge)%cpu(i)%num_cells , kind=wp )
            !
            bnd_edges(n)%ave_correction = rn / rd
            bnd_edges(n)%my_responsibility = &
                   ( edge(glb_edge)%controlling_partition == my_part )
            !
            cycle fill_bnd_edges
            !
          end if
        end do
        !
      end if
    end if
    !
  end do fill_bnd_edges
  call debug_timer(stop_timer,"section 1")
  call debug_timer(start_timer,"section 2")
  !
  ! Number of variables being stored for each node or edge point
  !
  nvar = nq * (nr + 1)
  !
  ! Allocate the node and edge point offset array
  !
  allocate ( nodpt_offset(1:lnode) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodpt_offset",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( edgpt_offset(1:ledge) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"edgpt_offset",1,__LINE__,__FILE__,ierr,error_message)
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
  call debug_timer(stop_timer,"section 2")
  call debug_timer(start_timer,"section 3")
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
      node_part = node(glb_node)%cpu(i)%partition
      !
      cpu_msk(node_part) = true
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
 !call output_nodpt_offsets(nodpt_offset,0,size(bnd_nodes))
 !call output_bnd_nodes(bnd_nodes,node,node_map)
 !nodpt_offset_save = nodpt_offset
  !
  ! Reshuffle nodpt_offset by summing the total number of nodes for
  ! all the communication nodes before the current one
  !
  ! SLOWEST METHOD
  !
 !do n = size(bnd_nodes),1,-1
 !  loc_node = bnd_nodes(n)%idx
 !  nodpt_offset(loc_node) = sum( nodpt_offset(1:loc_node-1) )
 !end do
 !call output_nodpt_offsets(nodpt_offset,1,size(bnd_nodes))
 !nodpt_offset = nodpt_offset_save
  !
  ! MUCH FASTER
  !
 !do n = size(bnd_nodes),1,-1
 !  loc_node = bnd_nodes(n)%idx
 !  nodpt_offset(loc_node) = sum( nodpt_offset(bnd_nodes(1:n-1)%idx) )
 !end do
 !call output_nodpt_offsets(nodpt_offset,2,size(bnd_nodes))
 !nodpt_offset = nodpt_offset_save
  !
  ! FASTEST METHOD
  !
  n1 = 0
  last_offset = 0
  do n = 1,size(bnd_nodes)
    n1 = n1 + last_offset
    loc_node = bnd_nodes(n)%idx
    last_offset = nodpt_offset(loc_node)
    nodpt_offset(loc_node) = n1
  end do
 !call output_nodpt_offsets(nodpt_offset,3,size(bnd_nodes))
  !
  if (maxval(nodpt_offset,dim=1)+1 /= size(bnd_nodes)) then
    write (error_message,1) "Something is wrong here with the node offset"
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
  end if
  call debug_timer(stop_timer,"section 3")
  call debug_timer(start_timer,"section 4")
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
      edge_part = edge(glb_edge)%cpu(i)%partition
      !
      cpu_msk(edge_part) = true
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
  ! Extract the last non-zero nodpt_offset
  !
  max_nodpt_offset = maxval(nodpt_offset,dim=1)
  !
  ! Reshuffle edgpt_offset by summing the total number of edge points for
  ! all the communication edges before the current one, making sure to
  ! also account for the last node offset by adding in max_nodpt_offset
 !!
 !! SLOWEST METHOD
 !!
 !do n = size(bnd_edges),1,-1
 !  loc_edge = bnd_edges(n)%idx
 !  edgpt_offset(loc_edge) = sum(edgpt_offset(1:loc_edge-1)) + max_nodpt_offset
 !end do
 !!
 !! MUCH FASTER
 !!
 !do n = size(bnd_edges),1,-1
 !  loc_edge = bnd_edges(n)%idx
 !  edgpt_offset(loc_edge) = sum( edgpt_offset(bnd_edges(1:n-1)%idx) ) + &
 !                           max_nodpt_offset
 !end do
  !
  ! FASTEST METHOD
  !
  n1 = max_nodpt_offset
  last_offset = 0
  do n = 1,size(bnd_edges)
    n1 = n1 + last_offset
    loc_edge = bnd_edges(n)%idx
    last_offset = edgpt_offset(loc_edge)
    edgpt_offset(loc_edge) = n1
  end do
  call debug_timer(stop_timer,"section 4")
  call debug_timer(start_timer,"section 5a")
 !write (iout,823) mypnum,"a"
  823 format (4x,"CPU-",i4.4," is here (",a,") !!",:,i10)
  !
  ! Finally, add the total node point offset into
  ! the offset for each of the edge points
  !
 !edgpt_offset = edgpt_offset + (size(bnd_nodes,kind=int_mpi)-1)
  !
  ! Count the number of communication pairs of the current processor
  ! with all other processors that have a higher rank
  !
  cpu_nums(:) = 0_int_mpi
#ifdef SPECIAL_FOR_PGI
  cpu_nums(my_part) = int( count( cpu_msk(my_part+1:ncpu) ) , kind=int_mpi )
#else
  cpu_nums(my_part) = count( cpu_msk(my_part+1:ncpu) , kind=int_mpi )
#endif
  !
  ! Collect the number of communication pairs for all processors
  !
 !call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_DATATYPE_NULL, &
  call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_INTEGER, &
                     cpu_nums,1_int_mpi,MPI_INTEGER,MPI_COMM_WORLD,mpierr)
  call debug_timer(stop_timer,"section 5a")
  call debug_timer(start_timer,"section 5b")
  !
  ! Count the total number of communication pairs for nodes and edges
  !
  npair = sum(cpu_nums)
  !
  ! Allocate the comm_pair array
  !
  allocate ( comm_pair(1:2,1:npair) , source=-1 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"comm_pair",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Fill in the section of comm_pair for the
  ! current processor before communicating
  !
  ! beginning index for the second dim of comm_pair for the current processor
  n1 = sum( cpu_nums(1:my_part-1) ) + 1
  ! ending    index for the second dim of comm_pair for the current processor
  n2 = sum( cpu_nums(1:my_part) )
  !
  ! processor of lower  rank that belongs to this communication pair
  comm_pair(1,n1:n2) = mypnum
  ! processor of higher rank that belongs to this communication pair
  comm_pair(2,n1:n2) = pack( intseq(my_part+1,ncpu) - 1 , &
                       mask=cpu_msk(my_part+1:ncpu) )
  !
  ! Create the displacement and count arrays for gathering the comm_pair array
  ! NOTE: We first need to multiply cpu_nums by the size of the first
  !       dimension of comm_pair because each pair contains that many
  !       pieces of information when collecting comm_pair
  !
  cpu_nums = cpu_nums * size(comm_pair,dim=1,kind=int_mpi)
  !
  cpu_disp = 0_int_mpi
  do n = 2,ncpu
    cpu_disp(n) = cpu_disp(n-1) + cpu_nums(n-1)
  end do
  call debug_timer(stop_timer,"section 5b")
  call debug_timer(start_timer,"section 5c")
  !
  ! Collect the individual comm_pair arrays for all processors
  !
 !call mpi_allgatherv(MPI_IN_PLACE,0_int_mpi,MPI_DATATYPE_NULL, &
  call mpi_allgatherv(MPI_IN_PLACE,0_int_mpi,MPI_INTEGER, &
                      comm_pair,cpu_nums,cpu_disp,mpi_inttyp, &
                      MPI_COMM_WORLD,mpierr)
  call debug_timer(stop_timer,"section 5c")
  call debug_timer(start_timer,"section 5d")
  !
  ! Allocate the array containing the communication information
  !
  allocate ( exch_neusp(1:npair) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"exch_neusp",1,__LINE__,__FILE__,ierr,error_message)
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
  call alloc_error(pname,"comm_pair",2,__LINE__,__FILE__,ierr,error_message)
  call debug_timer(stop_timer,"section 5d")
  call debug_timer(start_timer,"section 6: pair_loop")
  !
  ! Allocate all the arrays needed for creating
  ! the send and receive node and edge maps
  !
  allocate ( send_node_msk(1:lnode) , source=fals , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"send_node_msk",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( send_edge_msk(1:ledge) , source=fals , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"send_edge_msk",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Loop through each of the node/edge communication pairs and create the
  ! send and receive MPI datatypes for the two processors of this pair
  !
  pair_loop: do this_pair = 1,size(exch_neusp)
    !
    adj_cpu = exch_neusp(this_pair)%adj_cpu
    !
    if (adj_cpu < 0_int_mpi) cycle pair_loop
    !
    adj_part = int(adj_cpu,kind=kind(adj_part)) + 1
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
    call alloc_error(pname,"recv_node_map",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    rmin = remin(adj_part)
    rmax = remax(adj_part)
    allocate ( recv_edge_map(rmin:rmax) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_edge_map",1,__LINE__,__FILE__,ierr, &
                     error_message)
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
    call alloc_error(pname,"send_disp",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( recv_disp(1:scnt) , source=0_int_mpi , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_disp",1,__LINE__,__FILE__,ierr,error_message)
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
    call alloc_error(pname,"recv_node_map",2,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    deallocate ( recv_edge_map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_edge_map",2,__LINE__,__FILE__,ierr, &
                     error_message)
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
    send_disp(:) = send_disp(:) * block_size
    recv_disp(:) = recv_disp(:) * block_size
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
    call alloc_error(pname,"send_disp",2,__LINE__,__FILE__,ierr,error_message)
    !
    deallocate ( recv_disp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_disp",2,__LINE__,__FILE__,ierr,error_message)
    !
  end do pair_loop
  call debug_timer(stop_timer,"section 6: pair_loop")
  !
  ! Deallocate send_node_msk, send_edge_msk, nodpt_offset,
  ! and edgpt_offset before we leave
  !
  deallocate ( send_node_msk , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"send_node_msk",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( send_edge_msk , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"send_edge_msk",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( nodpt_offset , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodpt_offset",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( edgpt_offset , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"edgpt_offset",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (a)
  !
contains
subroutine output_nodpt_offsets(nodpt_offset,method,size_bnd_nodes)
  !
  integer, intent(in) :: nodpt_offset(:)
  integer, intent(in) :: method
  integer, intent(in) :: size_bnd_nodes
  !
  integer :: n,io,ierr
  character(len=100) :: fname
  !
continue
  !
  write (fname,1) mypnum,method
  !
  open (newunit=io,file=fname,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  !
  do n = 1,size(nodpt_offset)
    write (io,2) n,nodpt_offset(n)
  end do
  !
  write (io,3) maxval(nodpt_offset,dim=1)+1,size_bnd_nodes
  !
  close (io,iostat=ierr,iomsg=error_message)
  !
  ! Format Statements
  !
  1 format ("nodpt_offsets/",i4.4,".",i1,".dat")
  2 format ("nodpt_offset(",i4,") = ",i0)
  3 format ("maxval(nodpt_offset)+1 = ",i0,/, &
            "       size(bnd_nodes) = ",i0)
  !
end subroutine output_nodpt_offsets
subroutine output_bnd_nodes(bnd_nodes,node,node_map)
  !
  type(bnd_data_t), intent(in) :: bnd_nodes(:)
  type(node_t),     intent(in) :: node(:)
  type(map_t),      intent(in) :: node_map(:)
  !
  integer :: i,n,io,ierr
  integer :: loc_node,glb_node
  character(len=100) :: fname
  !
continue
  !
  write (fname,1) mypnum
  !
  open (newunit=io,file=fname,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  !
  do n = 1,size(bnd_nodes)
    loc_node = bnd_nodes(n)%idx
    glb_node = node_map(mypnum+1)%loc_to_glb(loc_node)
    if (allocated(node(glb_node)%cpu)) then
      write (io,2) n,loc_node,glb_node, &
            (node(glb_node)%cpu(i)%partition,i=1,size(node(glb_node)%cpu))
    else
      write (io,3) n,loc_node,glb_node
    end if
  end do
  !
  close (io,iostat=ierr,iomsg=error_message)
  !
  ! Format Statements
  !
  1 format ("nodpt_offsets/",i4.4,".4.dat")
  2 format ("bnd_nodes(",i4,")%idx = ",i0,"  ; node(gn=",i0, &
            ")%cpu(:)%partition = [",100(i0,:,","))
           !")%cpu(:)%partition = [",*(i0,:,","))
  3 format ("bnd_nodes(",i4,")%idx = ",i0,"  ; node(gn=",i0, &
            ")%cpu IS NOT ALLOCATED")
  !
end subroutine output_bnd_nodes
end subroutine create_node_and_edge_datatypes
!
!###############################################################################
!
subroutine exchange_faceusp(faceusp)
  !
  !.. Use Statements ..
  use geovar,  only : bface
  use flowvar, only : face_var
  !
  !.. Formal Arguments ..
!!!#define DEBUG_SEND_FACEUSP
#ifdef DEBUG_SEND_FACEUSP
  type(face_var), allocatable, intent(inout) :: faceusp(:)
#else
  type(face_var), dimension(:), intent(in) :: faceusp
#endif
  !
  !.. Local Scalars..
  integer :: n,nf,nsz,ibeg,iend
  !
  integer(int_mpi) :: send_cpu
  integer(int_mpi) :: recv_cpu
  integer(int_mpi) :: mesg_tag
  _MPI_DATA_TYPE_ :: send_dtype,recv_dtype
#ifdef DEBUG_SEND_FACEUSP
  integer :: i,j,k,m,ll
  type(face_var), allocatable :: tmp_faceusp(:)
#endif
  !
continue
#ifdef DEBUG_SEND_FACEUSP
  do j = 1,2
    fusp_send_buffer(:) = real(zero,kind=r8)
    if (j == 1) then
      tmp_faceusp = faceusp ! F2003 AUTO-REALLOCATION
      i = 0
      do ll = 1,2
        do m = 1,nq
          do nf = 1,size(bface,dim=2)
            do k = 1,size(faceusp(nf)%v,dim=2)
              i = i + 1
              faceusp(nf)%v(m,k,ll) = real(i,kind=wp)
            end do
          end do
        end do
      end do
    else
      call move_alloc( from=tmp_faceusp , to=faceusp )
    end if
#endif
  !
  iend = 0
  do nf = 1,size(bface,dim=2)
    if (bface(1,nf) /= bc_cpu_bnd) cycle
    nsz = size(faceusp(nf)%v(:,:,1))
    ibeg = iend + 1
    iend = iend + nsz
    fusp_send_buffer(ibeg:iend) = &
                 reshape( real(faceusp(nf)%v(:,:,1),kind=r8) , [nsz] )
  end do
#ifdef DEBUG_SEND_FACEUSP
    if (j == 1) then
      do i = 1,size(fusp_send_buffer)
        write (mypnum+3000,*) i,nint(fusp_send_buffer(i))
      end do
      flush (mypnum+3000)
    end if
  end do
#endif
!!!#define DEBUG_RECV_FACEUSP
#ifdef DEBUG_RECV_FACEUSP
  fusp_send_buffer(:) = real( intseq(1,size(fusp_send_buffer)) , kind=r8 )
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
    call mpi_irecv(fusp_recv_buffer,1_int_mpi,recv_dtype,send_cpu,mesg_tag, &
                   MPI_COMM_WORLD,fusp_rqst(n),mpierr)
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
    call mpi_isend(fusp_send_buffer,1_int_mpi,send_dtype,recv_cpu,mesg_tag, &
                   MPI_COMM_WORLD,fusp_rqst(n+ncomm),mpierr)
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
  use geovar,  only : bface
  use flowvar, only : face_var
  !
  !.. Formal Arguments ..
  type(face_var), dimension(:), intent(inout) :: faceusp
  !
  !.. Local Scalars ..
  integer :: ibeg,iend,nf
  !
  !.. Local Arrays ..
  integer :: nsz(1:2)
  !
  !.. Local Allocatable Arrays ..
  _MPI_STATUS_TYPE_, allocatable :: fstat(:)
  !
continue
  !
  ! Wait for all faceusp sends and receives to complete
  !
  fstat = reshape(fusp_istat,[size(fusp_istat)]) ! F2003 AUTO-REALLOCATION
  call mpi_waitall(size(fusp_rqst,kind=int_mpi),fusp_rqst,fstat,mpierr)
 !call mpi_waitall(size(fusp_rqst,kind=int_mpi),fusp_rqst, &
 !                 reshape(fusp_istat,[size(fusp_istat)]),mpierr)
#ifdef DEBUG_RECV_FACEUSP
  do nf = 1,size(fusp_recv_buffer)
    write (mypnum+7000,*) nf,nint(fusp_recv_buffer(nf))
  end do
  flush (mypnum+7000)
#endif
  !
  ! Update the right states for only the communication boundary faces
  !
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
  !
end subroutine wait_faceusp
!
!###############################################################################
!
subroutine exchange_facedusp(facedusp)
  !
  !.. Use Statements ..
  use geovar,  only : bface
  use flowvar, only : face_vec
  !
  !.. Formal Arguments ..
#ifdef DEBUG_SEND_FACEUSP
  type(face_vec), allocatable, intent(inout) :: facedusp(:)
#else
  type(face_vec), dimension(:), intent(in) :: facedusp
#endif
  !
  !.. Local Scalars..
  integer :: n,nf,nsz,ibeg,iend
  !
  integer(int_mpi) :: send_cpu
  integer(int_mpi) :: recv_cpu
  integer(int_mpi) :: mesg_tag
  _MPI_DATA_TYPE_ :: send_dtype,recv_dtype
#ifdef DEBUG_SEND_FACEUSP
  integer :: i,j,k,l,m,ll
  type(face_vec), allocatable :: tmp_facedusp(:)
#endif
  !
continue
#ifdef DEBUG_SEND_FACEUSP
  do j = 1,2
    fdusp_send_buffer(:) = real(zero,kind=r8)
    if (j == 1) then
      tmp_facedusp = facedusp ! F2003 AUTO-REALLOCATION
      i = 0
      do ll = 1,2
        do m = 1,nq
          do l = 1,nr
            do nf = 1,size(bface,dim=2)
              do k = 1,size(facedusp(nf)%v,dim=3)
                i = i + 1
                facedusp(nf)%v(l,m,k,ll) = real(i,kind=wp)
              end do
            end do
          end do
        end do
      end do
    else
      call move_alloc( from=tmp_facedusp , to=facedusp )
    end if
#endif
  !
  iend = 0
  do nf = 1,size(bface,dim=2)
    if (bface(1,nf) /= bc_cpu_bnd) cycle
    nsz = size(facedusp(nf)%v(:,:,:,1))
    ibeg = iend + 1
    iend = iend + nsz
    fdusp_send_buffer(ibeg:iend) = &
                 reshape( real(facedusp(nf)%v(:,:,:,1),kind=r8) , [nsz] )
  end do
#ifdef DEBUG_SEND_FACEUSP
    if (j == 1) then
      do i = 1,size(fdusp_send_buffer)
        write (mypnum+4000,*) i,nint(fdusp_send_buffer(i))
      end do
      flush (mypnum+4000)
    end if
  end do
#endif
#ifdef DEBUG_RECV_FACEUSP
  fdusp_send_buffer(:) = real( intseq(1,size(fdusp_send_buffer)) , kind=r8 )
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
    call mpi_irecv(fdusp_recv_buffer,1_int_mpi,recv_dtype,send_cpu,mesg_tag, &
                   MPI_COMM_WORLD,fdusp_rqst(n),mpierr)
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
    call mpi_isend(fdusp_send_buffer,1_int_mpi,send_dtype,recv_cpu,mesg_tag, &
                   MPI_COMM_WORLD,fdusp_rqst(n+ncomm),mpierr)
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
  use geovar,  only : bface
  use flowvar, only : face_vec
  !
  !.. Formal Arguments ..
  type(face_vec), dimension(:), intent(inout) :: facedusp
  !
  !.. Local Scalars ..
  integer :: ibeg,iend,nf
  !
  !.. Local Arrays
  integer :: nsz(1:3)
  !
  !.. Local Allocatable Arrays ..
  _MPI_STATUS_TYPE_, allocatable :: fstat(:)
  !
continue
  !
  ! Wait for all facedusp sends and receives to complete
  !
  fstat = reshape(fdusp_istat,[size(fdusp_istat)]) ! F2003 AUTO-REALLOCATION
  call mpi_waitall(size(fdusp_rqst,kind=int_mpi),fdusp_rqst,fstat,mpierr)
 !call mpi_waitall(size(fdusp_rqst,kind=int_mpi),fdusp_rqst, &
 !                 reshape(fdusp_istat,[size(fdusp_istat)]),mpierr)
#ifdef DEBUG_RECV_FACEUSP
  do nf = 1,size(fdusp_recv_buffer)
    write (mypnum+8000,*) nf,nint(fdusp_recv_buffer(nf))
  end do
  flush (mypnum+8000)
  call stop_gfr(stop_mpi,"wait_facedusp",__LINE__,__FILE__)
#endif
  !
  ! Update the right states for only the communication boundary faces
  !
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
  integer(int_mpi) :: send_tag
  integer(int_mpi) :: recv_tag
  _MPI_DATA_TYPE_ :: send_dtype,recv_dtype
  _MPI_STATUS_TYPE_ :: istat(1:_MPI_STATUS_SIZE_)
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
  call alloc_error(pname,"send_buffer",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( recv_buffer(1:n_node_edge_pts) , source=0.0_r8 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_buffer",1,__LINE__,__FILE__,ierr,error_message)
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
    if (exch_neusp(this_pair)%adj_cpu >= 0_int_mpi) then
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
      call mpi_sendrecv(send_buffer,1_int_mpi,send_dtype,adj_cpu,send_tag, &
                        recv_buffer,1_int_mpi,recv_dtype,adj_cpu,recv_tag, &
                        MPI_COMM_WORLD,istat(1),mpierr)
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
      recv_buffer(:) = 0.0_r8
      !
    end if
    !
    ! Have all processors wait for the current communication to complete
    ! before going on to the next communication pair.  This is to ensure
    ! there are no hangups and shouldnt be too much of a slowdown since
    ! this is only done for writing the solution files.
    !
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    !
  end do exchange_loop
  !
  if (allocated(send_buffer)) then
    deallocate ( send_buffer , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"send_buffer",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(recv_buffer)) then
    deallocate ( recv_buffer , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_buffer",2,__LINE__,__FILE__,ierr,error_message)
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
  integer(int_mpi) :: send_cpu,send_tag
  integer(int_mpi) :: recv_cpu,recv_tag
  _MPI_DATA_TYPE_ :: send_dtype,recv_dtype
  _MPI_STATUS_TYPE_ :: istat(1:_MPI_STATUS_SIZE_)
  !
  !.. Local Allocatable Arrays ..
  real(r8), allocatable, dimension(:) :: send_buffer
  real(r8), allocatable, dimension(:) :: recv_buffer
  !
  _MPI_STATUS_TYPE_, allocatable :: fstat(:)
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
  call alloc_error(pname,"send_buffer",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( recv_buffer(1:n_node_edge_pts) , &
             source=real(one,kind=kind(recv_buffer)) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"recv_buffer",1,__LINE__,__FILE__,ierr,error_message)
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
      call mpi_sendrecv(send_buffer,1_int_mpi,send_dtype,adj_cpu,send_tag, &
                        recv_buffer,1_int_mpi,recv_dtype,adj_cpu,recv_tag, &
                        MPI_COMM_WORLD,istat(1),mpierr)
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
      recv_buffer(:) = real( one , kind=kind(recv_buffer) )
      !
    end if
    !
    ! Have all processors wait for the current communication to complete
    ! before going on to the next communication pair.  This is to ensure
    ! there are no hangups and shouldnt be too much of a slowdown since
    ! this is only done for writing the solution files.
    !
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    !
  end do exchange_loop
  !
  if (allocated(send_buffer)) then
    deallocate ( send_buffer , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"send_buffer",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(recv_buffer)) then
    deallocate ( recv_buffer , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"recv_buffer",2,__LINE__,__FILE__,ierr,error_message)
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
    call mpi_irecv(fusp_recv_buffer,1_int_mpi,recv_dtype,send_cpu,mesg_tag, &
                   MPI_COMM_WORLD,fusp_rqst(n),mpierr)
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
    call mpi_isend(fusp_send_buffer,1_int_mpi,send_dtype,recv_cpu,mesg_tag, &
                   MPI_COMM_WORLD,fusp_rqst(n+ncomm),mpierr)
    !
  end do
  !
  ! 2.4) Wait for all faceusp sends and receives to complete
  !
  fstat = reshape(fusp_istat,[size(fusp_istat)]) ! F2003 AUTO-REALLOCATION
  call mpi_waitall(size(fusp_rqst,kind=int_mpi),fusp_rqst,fstat,mpierr)
 !call mpi_waitall(size(fusp_rqst,kind=int_mpi),fusp_rqst, &
 !                 reshape(fusp_istat,[size(fusp_istat)]),mpierr)
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
  call alloc_error(pname,"bnd_faces",2,__LINE__,__FILE__,ierr,error_message)
  !
end subroutine exch_connectivity
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
  _MPI_STATUS_TYPE_ :: istat(1:_MPI_STATUS_SIZE_)
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
    call alloc_error(pname,"sp_r8",1,__LINE__,__FILE__,ierr,error_message)
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
      call mpi_send(sp_r8(:,1:n_solpts),n,mpi_flttyp,0_int_mpi, &
                    mypnum,MPI_COMM_WORLD,mpierr)
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
      call alloc_error(pname,"global_r8",1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Have the root processor send its own solution to itself so
      ! it gets put in the right location within the global array.
      !
      call mpi_sendrecv(sp_r8(:,1:n_solpts),n,mpi_flttyp,0_int_mpi,n, &
                        global_r8,1_int_mpi,recv_sol(1),0_int_mpi,n, &
                        MPI_COMM_WORLD,istat(1),mpierr)
      !
      ! Now have the root processor receive the local
      ! solutions from all the other processors.
      !
      do i = 2,ncpu
        !
        icpu = i-1
        !
        call mpi_recv(global_r8,1_int_mpi,recv_sol(i),icpu,icpu, &
                      MPI_COMM_WORLD,istat(1),mpierr)
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
      call alloc_error(pname,"global_r8",2,__LINE__,__FILE__,ierr,error_message)
      !
    end if
    !
    deallocate ( sp_r8 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"sp_r8",2,__LINE__,__FILE__,ierr,error_message)
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
subroutine output_partitions(id,nodes_of_cell,nodes_of_cell_ptr,xyz_nodes,epart)
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: id
  integer,  intent(in) :: nodes_of_cell(:)
  integer,  intent(in) :: nodes_of_cell_ptr(:)
  real(wp), intent(in) :: xyz_nodes(:,:)
  !
  !.. Optional Arguments ..
  integer(idx_t), optional, intent(in) :: epart(:)
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
    call alloc_error(pname,"ipelem",1,__LINE__,__FILE__,ierr,error_message)
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
    call alloc_error(pname,"ipelem",1,__LINE__,__FILE__,ierr,error_message)
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
  call io_error(pname,trim(fname),1,__LINE__,__FILE__,ierr,error_message)
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
  call io_error(pname,trim(fname),2,__LINE__,__FILE__,ierr,error_message)
  !
  if (allocated(ipelem)) then
    deallocate ( ipelem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ipelem",2,__LINE__,__FILE__,ierr,error_message)
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
#ifndef DISABLE_DTIO
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
  if (allocated(bnd_nodes)) then
    !
    call dt%set( storage_size(bnd_nodes,kind=inttype) , &
                         size(bnd_nodes,kind=inttype) )
    !
    do i = lbound(bnd_nodes,dim=1),ubound(bnd_nodes,dim=1)
      !
      call dt%add( storage_size(bnd_nodes(i)%idx,kind=inttype) )
      call dt%add( storage_size(bnd_nodes(i)%ave_correction,kind=inttype) )
      call dt%add( storage_size(bnd_nodes(i)%my_responsibility,kind=inttype) )
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "bnd_nodes"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "bnd_nodes"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(bnd_edges)) then
    !
    call dt%set( storage_size(bnd_edges,kind=inttype) , &
                         size(bnd_edges,kind=inttype) )
    !
    do i = lbound(bnd_edges,dim=1),ubound(bnd_edges,dim=1)
      !
      call dt%add( storage_size(bnd_edges(i)%idx,kind=inttype) )
      call dt%add( storage_size(bnd_edges(i)%ave_correction,kind=inttype) )
      call dt%add( storage_size(bnd_edges(i)%my_responsibility,kind=inttype) )
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "bnd_edges"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "bnd_edges"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(bnd_faces)) then
    !
    call dt%set( storage_size(bnd_faces,kind=inttype) , &
                         size(bnd_faces,kind=inttype) )
    !
    do i = lbound(bnd_faces,dim=1),ubound(bnd_faces,dim=1)
      !
      call dt%add( storage_size(bnd_faces(i)%idx,kind=inttype) )
      call dt%add( storage_size(bnd_faces(i)%ave_correction,kind=inttype) )
      call dt%add( storage_size(bnd_faces(i)%my_responsibility,kind=inttype) )
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "bnd_faces"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "bnd_faces"
    !
  end if
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
  !-----------------------------------------------------------------------------
  !
  if (allocated(exch_neusp)) then
    !
    call dt%set( storage_size(exch_neusp,kind=inttype) , &
                         size(exch_neusp,kind=inttype) )
    !
    do i = lbound(exch_neusp,dim=1),ubound(exch_neusp,dim=1)
      !
      call dt%add( storage_size(exch_neusp(i)%adj_cpu,kind=inttype) )
      call dt%add( storage_size(exch_neusp(i)%send_tag,kind=inttype) )
      call dt%add( storage_size(exch_neusp(i)%recv_tag,kind=inttype) )
      call dt%add( storage_size(exch_neusp(i)%send_datatype,kind=inttype) )
      call dt%add( storage_size(exch_neusp(i)%recv_datatype,kind=inttype) )
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "exch_neusp"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "exch_neusp"
    !
  end if
  !
  !
  ! Check the memory of the arrays within this module
  !
  !--------------------------------------------------------------------
  if (allocated(fusp_rqst)) then
    call array%set( storage_size(fusp_rqst,kind=inttype) , &
                            size(fusp_rqst,kind=inttype) )
    call total%add(array)
    write (array_name,5) "fusp_rqst"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "fusp_rqst"
  end if
  !--------------------------------------------------------------------
  if (allocated(fdusp_rqst)) then
    call array%set( storage_size(fdusp_rqst,kind=inttype) , &
                            size(fdusp_rqst,kind=inttype) )
    call total%add(array)
    write (array_name,5) "fdusp_rqst"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "fdusp_rqst"
  end if
  !--------------------------------------------------------------------
  if (allocated(fusp_istat)) then
    call array%set( storage_size(fusp_istat,kind=inttype) , &
                            size(fusp_istat,kind=inttype) )
    call total%add(array)
    write (array_name,5) "fusp_istat"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "fusp_istat"
  end if
  !--------------------------------------------------------------------
  if (allocated(fdusp_istat)) then
    call array%set( storage_size(fdusp_istat,kind=inttype) , &
                            size(fdusp_istat,kind=inttype) )
    call total%add(array)
    write (array_name,5) "fdusp_istat"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "fdusp_istat"
  end if
  !--------------------------------------------------------------------
  if (allocated(fusp_send_buffer)) then
    call array%set( storage_size(fusp_send_buffer,kind=inttype) , &
                            size(fusp_send_buffer,kind=inttype) )
    call total%add(array)
    write (array_name,5) "fusp_send_buffer"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "fusp_send_buffer"
  end if
  !--------------------------------------------------------------------
  if (allocated(fusp_recv_buffer)) then
    call array%set( storage_size(fusp_recv_buffer,kind=inttype) , &
                            size(fusp_recv_buffer,kind=inttype) )
    call total%add(array)
    write (array_name,5) "fusp_recv_buffer"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "fusp_recv_buffer"
  end if
  !--------------------------------------------------------------------
  if (allocated(fdusp_send_buffer)) then
    call array%set( storage_size(fdusp_send_buffer,kind=inttype) , &
                            size(fdusp_send_buffer,kind=inttype) )
    call total%add(array)
    write (array_name,5) "fdusp_send_buffer"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "fdusp_send_buffer"
  end if
  !--------------------------------------------------------------------
  if (allocated(fdusp_recv_buffer)) then
    call array%set( storage_size(fdusp_recv_buffer,kind=inttype) , &
                            size(fdusp_recv_buffer,kind=inttype) )
    call total%add(array)
    write (array_name,5) "fdusp_recv_buffer"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "fdusp_recv_buffer"
  end if
  !--------------------------------------------------------------------
  if (allocated(recv_sol)) then
    call array%set( storage_size(recv_sol,kind=inttype) , &
                            size(recv_sol,kind=inttype) )
    call total%add(array)
    write (array_name,5) "recv_sol"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "recv_sol"
  end if
  !--------------------------------------------------------------------
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
#endif
end subroutine parallel_memory_usage
!
!###############################################################################
!
subroutine partition_memory_usage(iou,xadj,adjncy,epart,node,edge,face, &
                                  cells_with_node,cells_with_node_ptr, &
                                  cells_surr_cell,cells_surr_cell_ptr, &
                                  cells_with_edge,cells_with_edge_ptr, &
                                  nodes_on_edge,face_rotation,flx,fp, &
                                  node_map,edge_map,face_map)
  !
  !.. Use Statements ..
  use module_memory
  !
  !.. Formal Arguments ..
  integer,                       intent(in) :: iou
  integer(idx_t),   allocatable, intent(in) :: xadj(:)
  integer(idx_t),   allocatable, intent(in) :: adjncy(:)
  integer(idx_t),   allocatable, intent(in) :: epart(:)
  integer,          allocatable, intent(in) :: cells_with_node(:)
  integer,          allocatable, intent(in) :: cells_with_node_ptr(:)
  integer,          allocatable, intent(in) :: cells_surr_cell(:)
  integer,          allocatable, intent(in) :: cells_surr_cell_ptr(:)
  integer,          allocatable, intent(in) :: cells_with_edge(:)
  integer,          allocatable, intent(in) :: cells_with_edge_ptr(:)
  integer,          allocatable, intent(in) :: nodes_on_edge(:,:)
  integer,          allocatable, intent(in) :: face_rotation(:,:)
  type(map_t),      allocatable, intent(in) :: node_map(:)
  type(map_t),      allocatable, intent(in) :: edge_map(:)
  type(map_t),      allocatable, intent(in) :: face_map(:)
  type(node_t),     allocatable, intent(in) :: node(:)
  type(edge_t),     allocatable, intent(in) :: edge(:)
  type(face_t),     allocatable, intent(in) :: face(:)
  type(flux_pts_t), allocatable, intent(in) :: flx(:)
  type(fp_t),       allocatable, intent(in) :: fp
  !
  !.. Local Scalars ..
  integer :: i,j,k,loc_siz,glb_siz
  character(len=36) :: array_name
  type(memory) :: array,dt
  type(memory) :: total
  type(memory) :: edge_cpu
  type(memory) :: flxpts
  type(memory) :: left_fp
  type(memory) :: right_fp
  !
#ifndef DISABLE_DTIO
continue
  !
  call total%reset
  call edge_cpu%reset
  call flxpts%reset
  call left_fp%reset
  call right_fp%reset
  !
  write (iou,1)
  !
  !-----------------------------------------------------------------------------
  !
  ! Check the memory of the local allocatable derived type
  ! arrays within the subroutine partition_grid
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(node_map)) then
    !
    call dt%set( storage_size(node_map,kind=inttype) , &
                         size(node_map,kind=inttype) )
    !
    do i = lbound(node_map,dim=1),ubound(node_map,dim=1)
      if (allocated(node_map(i)%loc_to_glb)) then
        call dt%add( storage_size(node_map(i)%loc_to_glb,kind=inttype) , &
                             size(node_map(i)%loc_to_glb,kind=inttype) )
      end if
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "node_map"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "node_map"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(face_map)) then
    !
    call dt%set( storage_size(face_map,kind=inttype) , &
                         size(face_map,kind=inttype) )
    !
    do i = lbound(face_map,dim=1),ubound(face_map,dim=1)
      if (allocated(face_map(i)%loc_to_glb)) then
        call dt%add( storage_size(face_map(i)%loc_to_glb,kind=inttype) , &
                             size(face_map(i)%loc_to_glb,kind=inttype) )
      end if
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "face_map"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "face_map"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(edge_map)) then
    !
    call dt%set( storage_size(edge_map,kind=inttype) , &
                         size(edge_map,kind=inttype) )
    !
    do i = lbound(edge_map,dim=1),ubound(edge_map,dim=1)
      if (allocated(edge_map(i)%loc_to_glb)) then
        call dt%add( storage_size(edge_map(i)%loc_to_glb,kind=inttype) , &
                             size(edge_map(i)%loc_to_glb,kind=inttype) )
      end if
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "edge_map"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "edge_map"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(fp)) then
    !
    call dt%set( storage_size(fp,kind=inttype) , 1_inttype )
    !
    if (allocated(fp%geord)) then
      !
      call array%set( storage_size(fp%geord,kind=inttype), &
                              size(fp%geord,kind=inttype) )
      !
      call dt%add(array)
      !
      do j = lbound(fp%geord,dim=2),ubound(fp%geord,dim=2)
        do i = lbound(fp%geord,dim=2),ubound(fp%geord,dim=2)
          if (allocated(fp%geord(i,j)%rot)) then
            !
            call array%set( storage_size(fp%geord(i,j)%rot,kind=inttype) , &
                                    size(fp%geord(i,j)%rot,kind=inttype) )
            call dt%add(array)
            !
            do k = lbound(fp%geord(i,j)%rot,dim=1), &
                   ubound(fp%geord(i,j)%rot,dim=1)
              if (allocated(fp%geord(i,j)%rot(k)%pts)) then
                !
                call array%set( &
                       storage_size(fp%geord(i,j)%rot(k)%pts,kind=inttype) , &
                               size(fp%geord(i,j)%rot(k)%pts) )
                !
                call dt%add(array)
                !
              end if
            end do
            !
          end if
        end do
      end do
      !
    end if
    !
    call total%add(dt)
    write (array_name,5) "fp"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "fp"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(face)) then
    !
    call array%set( storage_size(face,kind=inttype) , &
                            size(face,kind=inttype) )
    !
    write (array_name,5) "face(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    do i = lbound(face,dim=1),ubound(face,dim=1)
      !
      call dt%add( storage_size(face(i)%geom,kind=inttype) )
      call dt%add( storage_size(face(i)%order,kind=inttype) )
      call dt%add( storage_size(face(i)%nodes,kind=inttype) , &
                           size(face(i)%nodes,kind=inttype) )
      call dt%add( storage_size(face(i)%left,kind=inttype) )
      call dt%add( storage_size(face(i)%left%cell,kind=inttype) )
      call dt%add( storage_size(face(i)%left%cell_face,kind=inttype) )
      call dt%add( storage_size(face(i)%right,kind=inttype) )
      call dt%add( storage_size(face(i)%right%cell,kind=inttype) )
      call dt%add( storage_size(face(i)%right%cell_face,kind=inttype) )
      !
      call array%set( storage_size(face(i)%left%fp_offset,kind=inttype) )
      call array%add( storage_size(face(i)%left%rotation,kind=inttype) )
      call left_fp%add(array)
      call flxpts%add(array)
      call dt%add(array)
      !
      call array%set( storage_size(face(i)%right%fp_offset,kind=inttype) )
      call array%add( storage_size(face(i)%right%rotation,kind=inttype) )
      call right_fp%add(array)
      call flxpts%add(array)
      call dt%add(array)
      !
    end do
    !
    write (array_name,5) "face%left%flxpts"
    write (iou,6) array_name,left_fp%get_units()
    !
    write (array_name,5) "face%right%flxpts"
    write (iou,6) array_name,right_fp%get_units()
    !
    write (array_name,5) "face%[left,right]%flxpts"
    write (iou,6) array_name,flxpts%get_units()
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
  if (allocated(node)) then
    !
    call array%set( storage_size(node,kind=inttype) , &
                            size(node,kind=inttype) )
    !
    write (array_name,5) "node(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    do i = lbound(node,dim=1),ubound(node,dim=1)
      !
      call dt%add( storage_size(node(i)%is_on_boundary,kind=inttype) )
      call dt%add( storage_size(node(i)%controlling_partition,kind=inttype) )
      !
      if (allocated(node(i)%cpu)) then
        !
        call array%set( storage_size(node(i)%cpu,kind=inttype) , &
                                size(node(i)%cpu,kind=inttype) )
        call dt%add(array)
        !
        do j = lbound(node(i)%cpu,dim=1),ubound(node(i)%cpu,dim=1)
          !
          call dt%add( storage_size(node(i)%cpu(j)%partition,kind=inttype) )
          call dt%add( storage_size(node(i)%cpu(j)%loc_node,kind=inttype) )
          call dt%add( storage_size(node(i)%cpu(j)%num_cells,kind=inttype) )
          !
        end do
        !
      end if
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "node"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "node"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(edge)) then
    !
    call array%set( storage_size(edge,kind=inttype) , &
                            size(edge,kind=inttype) )
    !
    write (array_name,5) "edge(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    do i = lbound(edge,dim=1),ubound(edge,dim=1)
      !
      call dt%add( storage_size(edge(i)%is_on_boundary,kind=inttype) )
      call dt%add( storage_size(edge(i)%controlling_partition,kind=inttype) )
      call dt%add( storage_size(edge(i)%order,kind=inttype) )
      !
      if (allocated(edge(i)%cpu)) then
        !
        call array%set( storage_size(edge(i)%cpu,kind=inttype) , &
                                size(edge(i)%cpu,kind=inttype) )
        call dt%add(array)
        !
        do j = lbound(edge(i)%cpu,dim=1),ubound(edge(i)%cpu,dim=1)
          !
          call dt%add( storage_size(edge(i)%cpu(j)%partition,kind=inttype) )
          call dt%add( storage_size(edge(i)%cpu(j)%loc_edge,kind=inttype) )
          call dt%add( storage_size(edge(i)%cpu(j)%num_cells,kind=inttype) )
          !
          if (allocated(edge(i)%cpu(j)%edgpts)) then
            call array%set( storage_size(edge(i)%cpu(j)%edgpts,kind=inttype), &
                                    size(edge(i)%cpu(j)%edgpts,kind=inttype) )
            call edge_cpu%add(array)
            call dt%add(array)
          end if
          !
        end do
        !
      end if
      !
    end do
    !
    write (array_name,5) "edge%cpu%edgpts"
    write (iou,6) array_name,edge_cpu%get_units()
    !
    call total%add(dt)
    write (array_name,5) "edge"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "edge"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(flx)) then
    !
    call dt%set( storage_size(flx,kind=inttype) , &
                         size(flx,kind=inttype) )
    !
    do i = lbound(flx,dim=1),ubound(flx,dim=1)
      if (allocated(flx(i)%pts)) then
        call dt%add( storage_size(flx(i)%pts,kind=inttype) , &
                             size(flx(i)%pts,kind=inttype) )
      end if
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "flx"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "flx"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  ! Check the memory of the local allocatable arrays
  ! within the subroutine partition_grid
  !
  !--------------------------------------------------------------------
  if (allocated(xadj)) then
    call array%set( storage_size(xadj,kind=inttype) , &
                            size(xadj,kind=inttype) )
    call total%add(array)
    write (array_name,5) "xadj"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "xadj"
  end if
  !--------------------------------------------------------------------
  if (allocated(adjncy)) then
    call array%set( storage_size(adjncy,kind=inttype) , &
                            size(adjncy,kind=inttype) )
    call total%add(array)
    write (array_name,5) "adjncy"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "adjncy"
  end if
  !--------------------------------------------------------------------
  if (allocated(epart)) then
    call array%set( storage_size(epart,kind=inttype) , &
                            size(epart,kind=inttype) )
    call total%add(array)
    write (array_name,5) "epart"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "epart"
  end if
  !--------------------------------------------------------------------
  if (allocated(cells_with_node)) then
    call array%set( storage_size(cells_with_node,kind=inttype) , &
                            size(cells_with_node,kind=inttype) )
    call total%add(array)
    write (array_name,5) "cells_with_node"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "cells_with_node"
  end if
  !--------------------------------------------------------------------
  if (allocated(cells_with_node_ptr)) then
    call array%set( storage_size(cells_with_node_ptr,kind=inttype) , &
                            size(cells_with_node_ptr,kind=inttype) )
    call total%add(array)
    write (array_name,5) "cells_with_node_ptr"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "cells_with_node_ptr"
  end if
  !--------------------------------------------------------------------
  if (allocated(cells_surr_cell)) then
    call array%set( storage_size(cells_surr_cell,kind=inttype) , &
                            size(cells_surr_cell,kind=inttype) )
    call total%add(array)
    write (array_name,5) "cells_surr_cell"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "cells_surr_cell"
  end if
  !--------------------------------------------------------------------
  if (allocated(cells_surr_cell_ptr)) then
    call array%set( storage_size(cells_surr_cell_ptr,kind=inttype) , &
                            size(cells_surr_cell_ptr,kind=inttype) )
    call total%add(array)
    write (array_name,5) "cells_surr_cell_ptr"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "cells_surr_cell_ptr"
  end if
  !--------------------------------------------------------------------
  if (allocated(cells_with_edge)) then
    call array%set( storage_size(cells_with_edge,kind=inttype) , &
                            size(cells_with_edge,kind=inttype) )
    call total%add(array)
    write (array_name,5) "cells_with_edge"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "cells_with_edge"
  end if
  !--------------------------------------------------------------------
  if (allocated(cells_with_edge_ptr)) then
    call array%set( storage_size(cells_with_edge_ptr,kind=inttype) , &
                            size(cells_with_edge_ptr,kind=inttype) )
    call total%add(array)
    write (array_name,5) "cells_with_edge_ptr"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "cells_with_edge_ptr"
  end if
  !--------------------------------------------------------------------
  if (allocated(nodes_on_edge)) then
    call array%set( storage_size(nodes_on_edge,kind=inttype) , &
                            size(nodes_on_edge,kind=inttype) )
    call total%add(array)
    write (array_name,5) "nodes_on_edge"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "nodes_on_edge"
  end if
  !--------------------------------------------------------------------
  if (allocated(face_rotation)) then
    call array%set( storage_size(face_rotation,kind=inttype) , &
                            size(face_rotation,kind=inttype) )
    call total%add(array)
    write (array_name,5) "face_rotation"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "face_rotation"
  end if
  !--------------------------------------------------------------------
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of local allocatable derived-types/arrays ", &
              "in partition_grid subroutine")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a,:,a)
  4 format (a,"(",a,",",i0,")%",a,:,"(",i0,")%mat")
  5 format ("'",a,"'")
  6 format (" Array Component: ",a,dt)
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine partition_memory_usage
!
!###############################################################################
!
subroutine report_parallel_memory_usage(xadj,adjncy,epart,node,edge,face, &
                                        cells_with_node,cells_with_node_ptr, &
                                        cells_surr_cell,cells_surr_cell_ptr, &
                                        cells_with_edge,cells_with_edge_ptr, &
                                        nodes_on_edge,face_rotation,flx,fp, &
                                        node_map,edge_map,face_map, &
                                        n,request_stop)
  !
  !.. Use Statements ..
  use geovar,        only : geovar_memory_usage
  use ovar,          only : ovar_memory_usage
  use flowvar,       only : flowvar_memory_usage
  use module_plot3d, only : plot3d_memory_usage
  use module_gmsh,   only : gmsh_memory_usage
  use module_cgns,   only : cgns_memory_usage
  !
  !.. Formal Arguments ..
  integer(idx_t),   allocatable, intent(in) :: xadj(:)
  integer(idx_t),   allocatable, intent(in) :: adjncy(:)
  integer(idx_t),   allocatable, intent(in) :: epart(:)
  integer,          allocatable, intent(in) :: cells_with_node(:)
  integer,          allocatable, intent(in) :: cells_with_node_ptr(:)
  integer,          allocatable, intent(in) :: cells_surr_cell(:)
  integer,          allocatable, intent(in) :: cells_surr_cell_ptr(:)
  integer,          allocatable, intent(in) :: cells_with_edge(:)
  integer,          allocatable, intent(in) :: cells_with_edge_ptr(:)
  integer,          allocatable, intent(in) :: nodes_on_edge(:,:)
  integer,          allocatable, intent(in) :: face_rotation(:,:)
  type(node_t),     allocatable, intent(in) :: node(:)
  type(map_t),      allocatable, intent(in) :: node_map(:)
  type(edge_t),     allocatable, intent(in) :: edge(:)
  type(map_t),      allocatable, intent(in) :: edge_map(:)
  type(face_t),     allocatable, intent(in) :: face(:)
  type(map_t),      allocatable, intent(in) :: face_map(:)
  type(flux_pts_t), allocatable, intent(in) :: flx(:)
  type(fp_t),       allocatable, intent(in) :: fp
  !
  !.. Optional Arguments ..
  integer,     optional, intent(in) :: n
  logical(lk), optional, intent(in) :: request_stop
  !
  !.. Local Scalars ..
  integer :: i,ierr,ncpu_digits
  !
  character(len=100) :: memfile
  character(len=100) :: frmt
  character(len=100) :: char_num
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "report_parallel_memory_usage"
  !
continue
  !
  frmt = "('memory_usage'"
  !
  if (ncpu > 1) then
    !
    write (char_num,1) ncpu
    ncpu_digits = max(1,len_trim(adjustl(char_num)))
    !
    write (frmt,10) trim(adjustl(frmt)),ncpu_digits,ncpu_digits
    !
  end if
  !
  if (present(n)) then
    write (frmt,11) trim(adjustl(frmt)),n
  end if
  !
  frmt = trim(adjustl(frmt)) // ",'.dat')"
  !
  if (ncpu > 1) then
    write (memfile,frmt) mypnum
  else
    write (memfile,frmt)
  end if
  !
  open (newunit=i,file=trim(adjustl(memfile)),status="replace", &
        action="write",form="formatted",position="rewind", &
        iostat=ierr,iomsg=error_message)
  call io_error(pname,memfile,1,__LINE__,__FILE__,ierr,error_message)
  !
  call geovar_memory_usage(i)
 !call ovar_memory_usage(i)
 !call plot3d_memory_usage(i)
 !call gmsh_memory_usage(i)
  call cgns_memory_usage(i)
 !call flowvar_memory_usage(i)
  call parallel_memory_usage(i)
  call partition_memory_usage(i,xadj,adjncy,epart,node,edge,face, &
                              cells_with_node,cells_with_node_ptr, &
                              cells_surr_cell,cells_surr_cell_ptr, &
                              cells_with_edge,cells_with_edge_ptr, &
                              nodes_on_edge,face_rotation,flx,fp, &
                              node_map,edge_map,face_map)
  !
  flush (i)
  !
  close (i,iostat=ierr,iomsg=error_message)
  call io_error(pname,memfile,2,__LINE__,__FILE__,ierr,error_message)
  !
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  if (present(request_stop)) then
    if (request_stop) then
      write (error_message,20)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
  end if
  !
  ! Format Statements
  !
  1  format (i0)
  10 format (a,",'.CPU-',i",i0,".",i0)
  11 format (a,",'.CP-",i0,"'")
  20 format ("After dumping memory usage")
  !
end subroutine report_parallel_memory_usage
!
!###############################################################################
!
pure function sum_map_sizes(map) result(return_value)
  !.. Formal Arguments ..
  type(map_t), intent(in) :: map(:)
  !.. Function Result ..
  integer :: return_value
  !.. Local Scalars ..
  integer :: n
continue
  return_value = 0
  do n = lbound(map,dim=1),ubound(map,dim=1)
    return_value = return_value + size(map(n)%loc_to_glb)
  end do
end function sum_map_sizes
!
!###############################################################################
!
subroutine debug_timer(time_flag,message,processor,part_num,acc)
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: time_flag
  !
  character(len=*), intent(in) :: message
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: processor
  integer, optional, intent(in) :: part_num
  integer, optional, intent(in) :: acc
  !
#ifdef PROFILE_ON
  !.. Local Scalars ..
  real(wp), save :: beg_wall_time
  real(wp)       :: end_wall_time
  real(wp)       :: wall_time
  !
  integer, save :: procedure_counter = 0
  integer       :: i,hours,minutes,seconds,total_sec
  integer       :: output_cpu
  !
  character(len=300) :: out_message
  !
  !.. Local Arrays ..
  real(wp), save :: beg_procedure_time(10)
  real(wp), save :: accumulated_time(10)
  real(wp), save :: accumulation_beg_time(10)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: leading_space = "                    "
  !
continue
  !
  output_cpu = 0
  if (present(processor)) output_cpu = processor
  !
  if (time_flag == start_timer) then
    !
    if (present(part_num)) then
      write (out_message,6) message,part_num
    else
      write (out_message,5) message
    end if
    !
    i = 3 + procedure_counter*2
    !
    beg_wall_time = mpi_wtime()
   !beg_wall_time = dclock()
   !call cpu_time( beg_wall_time )
    !
    if (mypnum == output_cpu) then
      write (iout,1) leading_space(1:i),trim(adjustl(out_message))
    end if
    !
  else if (time_flag == stop_timer) then
    !
    end_wall_time = mpi_wtime()
   !end_wall_time = dclock()
   !call cpu_time( end_wall_time )
    wall_time = end_wall_time - beg_wall_time
    !
    total_sec = nint(wall_time)
    !
    hours   = total_sec / 3600
    minutes = (total_sec - hours*3600) / 60
    seconds = total_sec - hours*3600 - minutes*60
    !
    if (present(part_num)) then
      write (out_message,6) message,part_num
    else
      write (out_message,5) message
    end if
    !
    i = 3 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,2) leading_space(1:i),wall_time,trim(adjustl(out_message))
    end if
    !
    ! Reset beg_wall_time in case this procedure is called twice
    ! in a row with time_flag set to stop_timer
    !
    beg_wall_time = mpi_wtime()
   !beg_wall_time = dclock()
   !call cpu_time( beg_wall_time )
    !
  else if (time_flag == reset_accumulated_time) then
    !
    if (present(acc)) then
      accumulated_time(acc) = zero
    end if
    !
  else if (time_flag == report_accumulated_time) then
    !
    if (present(part_num)) then
      write (out_message,6) message,part_num
    else
      write (out_message,5) message
    end if
    !
    i = 3 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,2) leading_space(1:i),accumulated_time(acc), &
                     trim(adjustl(out_message))
    end if
    !
  else if (time_flag == start_accumulation) then
    !
    if (present(acc)) then
      accumulation_beg_time(acc) = mpi_wtime()
    end if
    !
  else if (time_flag == stop_accumulation) then
    !
    if (present(acc)) then
      wall_time = mpi_wtime() - accumulation_beg_time(acc)
      accumulated_time(acc) = accumulated_time(acc) + wall_time
    end if
    !
  else if (time_flag == entering_procedure) then
    !
    i = 1 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,3) leading_space(1:i),trim(adjustl(message))
    end if
    !
    procedure_counter = procedure_counter + 1
    !
    beg_procedure_time(procedure_counter) = mpi_wtime()
   !beg_procedure_time(procedure_counter) = dclock()
   !call cpu_time( beg_procedure_time(procedure_counter) )
    !
  else if (time_flag == leaving_procedure) then
    !
    end_wall_time = mpi_wtime()
   !end_wall_time = dclock()
   !call cpu_time( end_wall_time )
    wall_time = end_wall_time - beg_procedure_time(procedure_counter)
    !
    procedure_counter = max( 0 , procedure_counter-1 )
    !
    i = 1 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,4) leading_space(1:i),trim(adjustl(message)),wall_time
    end if
    !
  end if
  !
  flush (iout)
  !
  ! Format Statements
  !
  1 format (a,"## Starting timer for: ",a)
  2 format (a,"     -- TIME WAS",f12.6," seconds for ",a)
  3 format (a,"=> Entering procedure : ",a)
  4 format (a," <= Leaving procedure : ",a,"   after",f12.6," seconds")
  5 format (a)
  6 format (a," for partition # ",i0)
  !
#endif
  !
end subroutine debug_timer
!
!##############################################################################
!
include "Functions/last.f90"
include "Functions/cell_solpts.f90"
!
!###############################################################################
!
end module parallel_mod
