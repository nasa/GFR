module parallel_mod
  !
  !.. Use Statements ..
  use module_kind_types
 !use mpivar
  !
  implicit none
  !
  private
  !
 !include "kind_defs.h"
 !include "constants.h"
  !
  public :: create_serial_cell_map
  public :: partition_grid
  public :: blocking_exchange
  public :: nonblocking_exchange
  public :: test_nonblocking_exchange
  public :: test_nonblocking_fexchange
  public :: collect_global_solpts
  public :: collect_global_solution
  public :: wait_faceusp
  public :: exchange_faceusp
  public :: wait_facedusp
  public :: exchange_facedusp
  !
  ! IDEA: 1. all processors send in their global versions of bface, xyz_nodes, 
  !          nodes_of_cell, nodes_of_cell_ptr, etc. to a single subroutine
  !       2. subroutine makes a copy of the global information necessary for
  !          output within the module
  !       3. Subroutine uses MeTiS to partition the grid
  !       4. Each processor then creates local versions of the arrays input
  !          to the subroutine which is then output from the subroutine
  !       5. Datatypes are needed to perform a parallel reduction of data 
  !          from the local partitions to a single solution array for output
  !       6. Datatypes are needed to exchange information along partition
  !          boundaries between processors
  !       7. We will probably have to create the solution points within each
  !          after the partitions are created and localized to the processors
  !
  integer, parameter :: nifac = 4
  !
  integer, save :: nr
  integer, save :: nq
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
  ! Number of solution points based on solution polynomial order
  !
  integer, save :: n_sp1de
  integer, save :: n_sp2dt
  integer, save :: n_sp2dq
  !
  ! MPI user-defined datatypes
  !
  integer, public, allocatable, dimension(:) :: snd_usp
  integer, public, allocatable, dimension(:) :: rcv_usp
  !
  ! Definitions of communication data
  !
  !         partition boundary faces shared
  !            by left_cpu and right_cpu
  !                        |
  !                        V
  !    __________________________________________
  !   |                    |                     |
  !   |   Left Partition   |   Right Partition   |
  !   |    on left_cpu     |    on right_cpu     |
  !   |____________________|_____________________|
  !                           
  !
  !     Left Partition ( comm_links_ptr(1) ) view
  !
  !           partition boundary face (lfa)
  !                        |
  !                        v
  !            x-----------o------------o
  !            |           |/           |
  !            |    lhc    |/    lgc    |
  !            |           |/           |
  !            x-----------o------------o
  !            <-interior->|<---ghost--->
  !
  !
  !    Right Partition ( comm_links_ptr(2) ) view
  !
  !           partition boundary face (rfa)
  !                        |
  !                        v
  !            x-----------o------------o
  !            |           |/           |
  !            |    rgc    |/    rhc    |
  !            |           |/           |
  !            x-----------o------------o
  !            <---ghost-->|<--interior->
  !
  ! n_com_bnd : number of unique communication boundaries
  ! n_com_lnk : number of individual communication "links" on a processor
  !             NOTE: A "link" is defined as a boundary face that
  !                   requires communication to another partition.
  !
  integer, save :: n_com_bnd
  integer, save :: n_com_lnk
  !
  ! comm_links_ptr : defines the communication parameters for a
  !                  partition boundary i_bnd. The data in this
  !                  array is identical across all processors.
  !     comm_links_ptr(1,i_bnd) = left_cpu  = left  partition number - 1
  !     comm_links_ptr(2,i_bnd) = right_cpu = right partition number - 1
  !     comm_links_ptr(3,i_bnd) = partition boundary faces shared
  !                               by left_cpu and right_cpu
  !
  integer, save, allocatable, dimension(:,:) :: comm_links_ptr
  !
  !  Example: 4 partitions
  !
  !      o---o---o---o   o---o---o---o---o
  !      |   |   |   |   |   |   |   |   |
  !      o---o---o---o   o---o---o---o---o
  !      |   | 1 |   |   |   | 2 |   |   |
  !      o---o---o---o   o---o---o---o---o
  !      |   |   |   |   |   |   |   |   |
  !      o---o---o---o   o---o---o---o---o
  !
  !      o---o---o---o   o---o---o---o---o
  !      |   |   |   |   |   |   |   |   |
  !      o---o---o---o   o---o---o---o---o
  !      |   | 3 |   |   |   | 4 |   |   |
  !      o---o---o---o   o---o---o---o---o
  !      |   |   |   |   |   |   |   |   |
  !      o---o---o---o   o---o---o---o---o
  !      |   |   |   |   |   |   |   |   |
  !      o---o---o---o   o---o---o---o---o
  !
  ! n_com_bnd = 4  [(1-2),(3-4),(1-3),(2-4)]
  !
  ! comm_links_ptr(1,1) = 0
  ! comm_links_ptr(2,1) = 1
  ! comm_links_ptr(3,1) = 3
  !
  ! comm_links_ptr(1,2) = 2
  ! comm_links_ptr(2,2) = 3
  ! comm_links_ptr(3,2) = 4
  !
  ! comm_links_ptr(1,3) = 0
  ! comm_links_ptr(2,3) = 2
  ! comm_links_ptr(3,3) = 3
  !
  ! comm_links_ptr(1,4) = 1
  ! comm_links_ptr(2,4) = 3
  ! comm_links_ptr(3,4) = 4
  !
  ! For processor 0 -> n_com_lnk = 6
  ! For processor 1 -> n_com_lnk = 7
  ! For processor 2 -> n_com_lnk = 7
  ! For processor 3 -> n_com_lnk = 8
  !
  ! comm_links : defines the connectivity of each partition boundary
  !              face, i_link. Only data relative to each processor
  !              is stored in this array. Therefore, for partition
  !              boundary face i_link which is shared by left_cpu and
  !              right_cpu, the data for i_link is identical between
  !              these two processors but the location of this data
  !              within the comm_links array is processor independent.
  !
  ! comm_links(1,1,i_link) = lfa = left face
  ! comm_links(2,1,i_link) = lhc = left host cell
  ! comm_links(3,1,i_link) = lgc = left ghost cell
  ! comm_links(1,2,i_link) = rfa = right face
  ! comm_links(2,2,i_link) = rhc = right host cell
  ! comm_links(3,2,i_link) = rgc = right ghost cell
  !
  integer, save, allocatable, dimension(:,:,:) :: comm_links
  !
  ! Number of communication partners
  !
  integer, save :: ncomm
  !
  ! Exchange partners
  !
  integer, save, allocatable, dimension(:) :: exch_cpu
  integer, save, allocatable, dimension(:) :: exchf_cpu
  !
  ! Send and receive cell maps and pointers
  !
  integer, save, allocatable, dimension(:) :: send_cell_ptr
  integer, save, allocatable, dimension(:) :: send_cell
  integer, save, allocatable, dimension(:) :: recv_cell_ptr
  integer, save, allocatable, dimension(:) :: recv_cell
  !
  ! Send and receive solution point maps and pointers
  !
  integer, save, allocatable, dimension(:) :: send_pts_ptr
  integer, save, allocatable, dimension(:) :: send_pts
  integer, save, allocatable, dimension(:) :: recv_pts_ptr
  integer, save, allocatable, dimension(:) :: recv_pts
  !
  ! Array that contains the offsets for the solution points
  ! within a boundary cell depending on cell geometry type
  ! and the face of the cell on the boundary face.
  !
  integer, save, allocatable, dimension(:,:,:) :: geom_offset
  !
  ! Array to store the orientation of each local ghost
  ! cell with respect to its corresponding global cell
  !
  integer, save, allocatable, dimension(:) :: ghost_face
  !
  ! Test MPI datatypes for exchanging information
  !
  integer, save, allocatable, dimension(:) :: test_send
  integer, save, allocatable, dimension(:) :: test_recv
  integer, save, allocatable, dimension(:) :: testf_send
  integer, save, allocatable, dimension(:) :: testf_recv
  !
  ! MPI datatypes for exchanging solution data between adjacent partitions
  !
  integer, save, allocatable, dimension(:) :: send_usp
  integer, save, allocatable, dimension(:) :: recv_usp
  !
  integer, save, allocatable, dimension(:) :: send_fusp
  integer, save, allocatable, dimension(:) :: recv_fusp
  !
  integer, save, allocatable, dimension(:) :: send_fdusp
  integer, save, allocatable, dimension(:) :: recv_fdusp
  !
  ! Request and status arrays for non-blocking exchanges
  !
  integer, save, allocatable, dimension(:)   :: fusp_rqst
  integer, save, allocatable, dimension(:,:) :: fusp_istat
  !
  integer, save, allocatable, dimension(:)   :: fdusp_rqst
  integer, save, allocatable, dimension(:,:) :: fdusp_istat
  !
  ! Exchange buffer
  !
  real(r8), save, allocatable, dimension(:,:)       :: exch_buffer
  real(r8), save, allocatable, dimension(:,:,:,:)   :: fusp_buffer
  real(r8), save, allocatable, dimension(:,:,:,:,:) :: fdusp_buffer
  !
  ! MPI datatype for collecting the global
  ! solution point coordinates on all processors
  !
  integer, save, allocatable, dimension(:) :: recv_xyz
  !
  ! MPI datatype for collecting the global solution on all processors
  !
  integer, save, allocatable, dimension(:) :: recv_sol
  !
  ! Local to global cell mappings - These need to be saved
  ! for use at various points (e.g., reading a restart file)
  !
  integer, public, save, allocatable, dimension(:) :: cell_map
  integer, public, save, allocatable, dimension(:) :: cell_map_ptr
 !!
 !! Generic interface for routine to reallocate an array and increase its size
 !!
 !interface reallocate
 !  module procedure reallocate_1, &
 !                   reallocate_2
 !end interface reallocate
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
continue
  !
  allocate ( cell_map_ptr(1:2) , source=[0,ncell] , stat=ierr )
  call alloc_error("create_serial_cell_map","cell_map_ptr",1,ierr)
  allocate ( cell_map(1:ncell) , source=[1:ncell] , stat=ierr )
  call alloc_error("create_serial_cell_map","cell_map",1,ierr)
  !
end subroutine create_serial_cell_map
!
!###############################################################################
!
subroutine partition_grid(npart,n_order,metis_option,bface,nodes_of_cell, &
                          nodes_of_cell_ptr,xyz_nodes)
  !
  !.. Use Statements ..
  use iso_c_binding, only : c_float
  use iso_c_binding, only : c_double
  !
  !.. Formal Arguments ..
  integer,                                   intent(in) :: npart
  integer,                                   intent(in) :: n_order
  integer,                                   intent(in) :: metis_option
  integer,   allocatable, dimension(:),   intent(inout) :: nodes_of_cell
  integer,   allocatable, dimension(:),   intent(inout) :: nodes_of_cell_ptr
  integer,   allocatable, dimension(:,:), intent(inout) :: bface
  real(vbp), allocatable, dimension(:,:), intent(inout) :: xyz_nodes
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,ierr
  integer :: wgtflag,numflag,vol
  character(len=19) :: method
  character(len=12) :: weight
  character(len=20) :: minmzd
  !
  !.. Local Arrays ..
  integer,       dimension(1:5)     :: options
  real(c_float), dimension(1:npart) :: tpwgts
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:)   :: xadj
  integer, allocatable, dimension(:)   :: vwgt
  integer, allocatable, dimension(:)   :: adjncy
  integer, allocatable, dimension(:)   :: adjwgt
  integer, allocatable, dimension(:)   :: epart
  integer, allocatable, dimension(:,:) :: cell_rmap
  integer, allocatable, dimension(:)   :: node_map
  integer, allocatable, dimension(:)   :: node_map_ptr
  integer, allocatable, dimension(:,:) :: node_rmap
  integer, allocatable, dimension(:)   :: face_map
  integer, allocatable, dimension(:)   :: face_map_ptr
  integer, allocatable, dimension(:,:) :: face_rmap
  !
  ! Unlike the cells_on_face array defined in the module
  ! geovar, cells_on_face is defined in this routine as:
  !
  ! For face i_face, with face vector going from node 1 to node 2, ...
  !   cells_on_face(1,i_face) = attached cell to the left  of i_face
  !   cells_on_face(2,i_face) = attached cell to the right of i_face
  !   cells_on_face(3,i_face) = node 1 defining i_face
  !   cells_on_face(4,i_face) = node 2 defining i_face
  !
  integer, allocatable, dimension(:,:) :: cells_on_face
  !
continue
  !
  ! Get the number of solution points based on n_order
  !
  n_sp1de = n_order+1
  n_sp2dt = n_sp1de*(n_sp1de+1)/2
  n_sp2dq = n_sp1de*n_sp1de
  !
  ! Get the size of the global grid
  !
  nr = size(xyz_nodes,dim=1)
  nq = nr + 2
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
  call create_dual_graph(nodes_of_cell,nodes_of_cell_ptr, &
                         bface,xadj,adjncy,cells_on_face)
  !
  ! Create the weights for the vertices of the dual graph
  !
  call create_grid_weights(metis_option,last(xadj)-1,vwgt,adjwgt, &
                           nodes_of_cell_ptr)
  !
  ! Allocate the epart array to store the cell partition numbers
  !
  allocate ( epart(1:ngcell) , source=0 , stat=ierr )
  call alloc_error("partition_grid","epart",1,ierr)
 !if (mypnum == 0) then
 !  call output_partitions(99,nodes_of_cell,nodes_of_cell_ptr,xyz_nodes)
 !end if
  !
  ! Define other MeTiS variables
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
  ! Partition the grid using one of the MeTiS routines
  !
  if (metis_option == 1) then
    !
    ! This MeTiS routine computes non-weighted partitions using a
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
    ! This MeTiS routine computes non-weighted partitions using
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
    ! This MeTiS routine computes non-weighted partitions using a
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
    ! This MeTiS routine computes weighted partitions using a
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
    ! sum to one. After looking at the MeTiS source code, this is
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
  end if
  !
  ! Deallocate no longer needed MeTiS arrays
  !
  deallocate ( xadj , stat=ierr )
  call alloc_error("partition_grid","xadj",2,ierr)
  deallocate ( adjncy , stat=ierr )
  call alloc_error("partition_grid","adjncy",2,ierr)
  deallocate ( vwgt , stat=ierr )
  call alloc_error("partition_grid","vwgt",2,ierr)
  deallocate ( adjwgt , stat=ierr )
  call alloc_error("partition_grid","adjwgt",2,ierr)
  !
  ! I dont trust that the MeTiS routines will compute identical
  ! partitioning on each processor so have the root process check
  ! its local partitioning for errors and then broadcast that
  ! partitioning to the other processes if none are found.
  !
  if (mypnum == 0) then
    !
    ! Check that all processors got part of the global grid
    !
    do n = 1,npart
      if (any(epart == n)) cycle
      write (iout,100) n
      call stop_recon(abort,"partition_grid")
    end do
    !
    ! Output the results of partitioning
    !
    write (iout,101) method,npart,weight,minmzd,vol
    !
  end if
  !
  ! Distribute partitions across all processors
  !
  call mpi_bcast(epart,ngcell,mpi_inttyp,0,mpi_comm_world,mpierr)
  !
  ! Output the partition grid to a tecplot file
  !
 !if (mypnum == 0) then
 !  call output_partitions(mypnum,nodes_of_cell,nodes_of_cell_ptr, &
 !                         xyz_nodes,epart)
 !end if
  !
  ! Find the faces that separate two adjacent partitions
  ! and add these faces to boundary faces array
  !
  call find_partition_boundaries(nodes_of_cell,nodes_of_cell_ptr,epart, &
                                 cells_on_face,bface)
  !
  ! Now that we have created boundary conditions for faces
  ! located on partition boundaries, use epart to create
  ! mappings from local to global block orderings.
  ! This also includes creating the reverse mapping arrays
  ! or the maps from global to local block orderings
  ! NOTE: cell_map and cell_map_ptr are passed here as if they are
  !       local arrays but they are actually module arrays because
  !       they are needed later after this routine is completed.
  !
  call create_grid_mappings(npart,epart,cells_on_face, &
                            nodes_of_cell,nodes_of_cell_ptr, &
                            cell_map,cell_map_ptr,cell_rmap, &
                            node_map,node_map_ptr,node_rmap, &
                            face_map,face_map_ptr,face_rmap)
  !
  ! Now that the mappings have been created, we
  ! need to localize the grid on each processor
  !
  call localize_grid(cell_map,cell_map_ptr,cell_rmap, &
                     node_map,node_map_ptr,node_rmap, &
                     face_map,face_map_ptr,face_rmap, &
                     nodes_of_cell,nodes_of_cell_ptr, &
                     bface,xyz_nodes)
  !
  ! Before creating the MPI datatypes, we need to create
  ! the offset array for the solution points within a cell
  !
  call create_geometry_offsets
  !
  ! Finally, create MPI datatypes for exchanging data between processes
  !
  call create_face_nonblocking_datatypes(bface)
  call create_nonblocking_exchange_datatypes(nodes_of_cell, &
                                             nodes_of_cell_ptr,bface)
  !
 !call create_exchange_datatypes()
  !
  ! Format Statements
  !
  100 format(/,' Error: Processor ',i0,' was excluded from the          ',/, &
               ' partitioning of the elements. More than likely too many',/, &
               ' processors were used in the paritioning of this mesh.  ',/, &
               ' Please decrease number of processors and try again.    ',/)
  101 format (/,'-------------------------------------------------------',//, &
                '         Used MeTiS 4.0 to partition the graph         ',/, &
                '    using a multilevel ',a,' algorithm    ',/, &
                '           into ',i0,' ',a,' partitions.          ',/, &
                '     The minimized ',a,' was ',i0,'.    ',//, &
                '-------------------------------------------------------',/)
  !
end subroutine partition_grid
!
!###############################################################################
!
subroutine create_dual_graph(nodes_of_cell,nodes_of_cell_ptr, &
                             bface,xadj,adjncy,cells_on_face)
  !
  !.. Formal Arguments ..
  integer,              dimension(:),      intent(in) :: nodes_of_cell
  integer,              dimension(:),      intent(in) :: nodes_of_cell_ptr
  integer,              dimension(:,:),    intent(in) :: bface
  integer, allocatable, dimension(:),   intent(inout) :: xadj
  integer, allocatable, dimension(:),   intent(inout) :: adjncy
  integer, allocatable, dimension(:,:), intent(inout) :: cells_on_face
  !
  !.. Local Scalars ..
  integer :: ic,i1,i2,ip,ip1,ip2,is,ieadj
  integer :: jc,j1,j2,n,n1,n2,nf,nadj,ierr
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:) :: lpoin
  integer, allocatable, dimension(:) :: cells_surr_cell
  integer, allocatable, dimension(:) :: cells_surr_cell_ptr
  integer, allocatable, dimension(:) :: cells_with_node
  integer, allocatable, dimension(:) :: cells_with_node_ptr
  !
continue
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Find the cells that surround each node
  !
  allocate ( cells_with_node_ptr(1:ngnode+1) , source=0 , stat=ierr )
  call alloc_error("create_dual_graph","cells_with_node_ptr",1,ierr)
  !
  ! Loop over the cells, storing 'ahead'
  !
  do n = 1,ngcell
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    do ic = n1,n2
      ip = nodes_of_cell(ic) + 1
      cells_with_node_ptr(ip) = cells_with_node_ptr(ip) + 1
    end do
  end do
  !
  ! Reshuffle cells_with_node_ptr
  !
  do n = 2,size(cells_with_node_ptr)
    cells_with_node_ptr(n) = cells_with_node_ptr(n) + cells_with_node_ptr(n-1)
  end do
  !
  allocate ( cells_with_node(1:last(cells_with_node_ptr)) , source=0 , stat=ierr )
  call alloc_error("create_dual_graph","cells_with_node",1,ierr)
  !
  ! Now store the cells surrounding each node in cells_with_node
  !
  do n = 1,ngcell
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    do ic = n1,n2
      ip = nodes_of_cell(ic)
      is = cells_with_node_ptr(ip) + 1
      cells_with_node_ptr(ip) = is
      cells_with_node(is) = n
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
  allocate ( cells_surr_cell_ptr(1:ngcell+1) , source=0 , stat=ierr )
  call alloc_error("create_dual_graph","cells_surr_cell_ptr",1,ierr)
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
  allocate ( cells_surr_cell(1:last(cells_surr_cell_ptr)) , source=0 , stat=ierr )
  call alloc_error("create_dual_graph","cells_surr_cell",1,ierr)
  !
  ! Loop to find cells_surr_cells
  !
  allocate ( lpoin(1:ngnode) , source=0 , stat=ierr )
  call alloc_error("create_dual_graph","lpoin",1,ierr)
  !
  do n = 1,ngcell
    !
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    !
    ! Loop over the nodes (representative of faces) for the current cell
    !
    do nf = n1,n2
      !
      ! Nodes of the current face
      !
      ip1 = nodes_of_cell(nf)
      ip2 = nodes_of_cell(merge(nf+1,n1,nf/=n2))
      !
      ! Mark these nodes
      !
      lpoin(ip1) = 1
      lpoin(ip2) = 1
      !
      ! Loop over the elements that surround node ip1
      !
      i1 = cells_with_node_ptr(ip1)+1
      i2 = cells_with_node_ptr(ip1+1)
      !
      ieadj = 0
      !
      find_cell: do ic = i1,i2
        jc = cells_with_node(ic)
        if (jc == n) cycle
        j1 = nodes_of_cell_ptr(jc)+1
        j2 = nodes_of_cell_ptr(jc+1)
        if (sum(lpoin(nodes_of_cell(j1:j2))) == 2) then
          ieadj = jc
          exit find_cell
        end if
      end do find_cell
      !
      ! Found the adjacent cell
      !
      cells_surr_cell(cells_surr_cell_ptr(n)+nf-n1+1) = ieadj
      !
      ! Reset lpoin for the next cell
      !
      lpoin(ip1) = 0
      lpoin(ip2) = 0
      !
    end do
    !
  end do
  !
  ! We are finished with lpoin, cells_with_node,
  ! and cells_with_node_ptr so deallocate
  !
  deallocate ( lpoin , stat=ierr )
  call alloc_error("create_dual_graph","lpoin",2,ierr)
  deallocate ( cells_with_node , stat=ierr )
  call alloc_error("create_dual_graph","cells_with_node",2,ierr)
  deallocate ( cells_with_node_ptr , stat=ierr )
  call alloc_error("create_dual_graph","cells_with_node_ptr",2,ierr)
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! We cant just use cells_surr_cell_ptr and cells_surr_cell as xadj
  ! and adjncy because cells_surr_cell is zero for faces of cells that
  ! that correspond to boundary faces. In the loop "find_cell" above,
  ! if the sum(lpoin(nodes_of_cell(j1:j2))) term never equals 2 then
  ! this loop exits with ieadj still equal to 0. We cannot send adjncy
  ! containing zeros to the MeTiS routines so these need to be removed
  ! and the locations in xadj also have to be adjusted accordingly.
  ! Additionally, have xadj to start the indexing for adjncy at 1
  ! instead of 0 like cells_surr_cell_ptr (i.e., xadj(1)=1) since
  ! this is what MeTiS expects when using Fortran array ordering.
  !
 !do n = 1,ngcell
 !  n1 = cells_surr_cell_ptr(n)+1
 !  n2 = cells_surr_cell_ptr(n+1)
 !  write (101,'(4(3x,i8))') cells_surr_cell(n1:n2)
 !end do
  !
  allocate ( xadj(1:ngcell+1) , source=0 , stat=ierr )
  call alloc_error("create_dual_graph","xadj",1,ierr)
  xadj(1) = 1
  !
  nadj = count( cells_surr_cell(:) /= 0 )
  allocate ( adjncy(1:nadj) , source=0 , stat=ierr )
  call alloc_error("create_dual_graph","adjncy",1,ierr)
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
  ! Find the total number of faces and then allocate cells_on_face
  !
  ngface = ngfbnd
  do n = 1,ngcell
    n1 = cells_surr_cell_ptr(n)+1
    n2 = cells_surr_cell_ptr(n+1)
    ngface = ngface + count( cells_surr_cell(n1:n2) > n )
  end do
  !
  allocate ( cells_on_face(1:nifac,1:ngface) , source=0 , stat=ierr )
  call alloc_error("create_dual_graph","cells_on_face",1,ierr)
  !
  ! Create the cells_on_face array
  !
  ! Get the cells on each boundary face
  !
  do nf = 1,ngfbnd
    !
    ic = bface(2,nf)
    !
    n1 = cells_surr_cell_ptr(ic)+1
    n2 = cells_surr_cell_ptr(ic+1)
    !
    ip = bface(4,nf) + n1 - 1
    !
    cells_on_face(1,nf) = ic                                     ! host  cell
    cells_on_face(2,nf) = ngcell + nf                            ! ghost cell
    cells_on_face(3,nf) = nodes_of_cell(ip)                      ! node 1
    cells_on_face(4,nf) = nodes_of_cell( merge(ip+1,n1,ip/=n2) ) ! node 2
    !
  end do
  !
  ! Get the cells for the interior faces
  !
  nf = ngfbnd
  do ic = 1,ngcell
    !
    n1 = cells_surr_cell_ptr(ic)+1
    n2 = cells_surr_cell_ptr(ic+1)
    !
    do ip = n1,n2
      !
      jc = cells_surr_cell(ip)
      !
      if (jc > ic .and. jc <= ngcell) then
        !
        nf = nf + 1
        !
        cells_on_face(1,nf) = ic                                     ! left  cell
        cells_on_face(2,nf) = jc                                     ! right cell
        cells_on_face(3,nf) = nodes_of_cell(ip)                      ! node 1
        cells_on_face(4,nf) = nodes_of_cell( merge(ip+1,n1,ip/=n2) ) ! node 2
        !
      end if
      !
    end do
    !
  end do
  !
  ! Deallocate cells_surr_cell and cells_surr_cell_ptr before we leave
  !
  deallocate ( cells_surr_cell , stat=ierr )
  call alloc_error("create_dual_graph","cells_surr_cell",2,ierr)
  deallocate ( cells_surr_cell_ptr , stat=ierr )
  call alloc_error("create_dual_graph","cells_surr_cell_ptr",2,ierr)
  !
end subroutine create_dual_graph
!
!###############################################################################
!
subroutine create_grid_weights(metis_option,nadjfac, &
                               vwgt,adjwgt,nodes_of_cell_ptr)
  !
  !.. Formal Arguments ..
  integer,                               intent(in) :: nadjfac
  integer,                               intent(in) :: metis_option
  integer,              dimension(:),    intent(in) :: nodes_of_cell_ptr
  integer, allocatable, dimension(:), intent(inout) :: vwgt
  integer, allocatable, dimension(:), intent(inout) :: adjwgt
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,np,ierr
  !
continue
  !
  ! Allocate and initialize vwgt
  ! NOTE: vwgt gives the computational weight of each cell
  !
  allocate ( vwgt(1:ngcell) , source=0 , stat=ierr )
  call alloc_error("create_grid_weights","vwgt",1,ierr)
  !
  !############################################################
  !############################################################
  ! We need to figure out the computational difference between
  ! triangle and quad cells and populate vwgt accordingly 
  !############################################################
  !############################################################
  !
  ! Allocate and initialize adjwgt
  !
  if ( any(metis_option == (/1,2/)) ) then
    !
    ! If minimizing the edgecut, adjwgt gives the weights of the faces
    !
    !###############################################################
    !###############################################################
    ! Not sure if we need any weighting on the faces so set to zero
    allocate ( adjwgt(1:nadjfac) , source=0 , stat=ierr )
    call alloc_error("create_grid_weights","adjwgt",1,ierr)
    !###############################################################
    !###############################################################
    !
  else if ( any(metis_option == (/3,4/)) ) then
    !
    ! If minimizing the communication volume, adjwgt gives the communication
    ! weight of the cells (i.e., how much data has to be sent to/from this cell)
    !
    allocate ( adjwgt(1:ngcell) , source=0 , stat=ierr )
    call alloc_error("create_grid_weights","adjwgt",1,ierr)
    !
    ! Set the communication weight of the cells to
    ! the number of solution points for each cell
    !
    do n = 1,ngcell
      !
      n1 = nodes_of_cell_ptr(n)+1
      n2 = nodes_of_cell_ptr(n+1)
      np = n2-n1+1
      !
      adjwgt(n) = merge(n_sp2dt,n_sp2dq,np==3)
      !
    end do
    !
  end if
  !
end subroutine create_grid_weights
!
!###############################################################################
!
subroutine find_partition_boundaries(nodes_of_cell,nodes_of_cell_ptr,epart, &
                                     cells_on_face,bface)
  !
  !.. Formal Arguments ..
  integer,              dimension(:),      intent(in) :: nodes_of_cell
  integer,              dimension(:),      intent(in) :: nodes_of_cell_ptr
  integer,              dimension(:),      intent(in) :: epart
  integer, allocatable, dimension(:,:), intent(inout) :: cells_on_face
  integer, allocatable, dimension(:,:), intent(inout) :: bface
  !
  !.. Local Scalars ..
  integer :: new_ngfbnd,new_ngface,ierr
  integer :: n,nf,n1,n2,np,np1,np2
  integer :: i,ic,i1,i2,ip,ifa,pf
  integer :: j,jc,j1,j2,jp,jfa,pc
  !
  !.. Local Arrays ..
  integer, dimension(1:4) :: lnode
  !
  !.. Local Allocatable Arrays ..
  logical, allocatable, dimension(:)   :: msk
  integer, allocatable, dimension(:,:) :: temp_array
  !
continue
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  allocate ( msk(1:ngface) , source=fals , stat=ierr )
  call alloc_error("find_partition_boundaries","msk",1,ierr)
  !
  ! Find the interior faces that divide two partitions
  !
  do nf = ngfbnd+1,ngface
    msk(nf) = epart( cells_on_face(1,nf) ) .ne. epart( cells_on_face(2,nf) )
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
  ! Collect and add the new boundary faces to bface
  !
  nf = ngfbnd + 1
  do n = ngfbnd+1,ngface
    !
    if (.not. msk(n)) cycle
    !
    ic = cells_on_face(1,n)
    jc = cells_on_face(2,n)
    np1 = cells_on_face(3,n)
    np2 = cells_on_face(4,n)
    !
    ! 1. The left side of the face
    !
    i1 = nodes_of_cell_ptr(ic)+1
    i2 = nodes_of_cell_ptr(ic+1)
    ip = i2-i1+1
    !
    ! Find which face of the left cell is on the boundary
    !
    lnode(:) = 0
    do i = i1,i2
      if (nodes_of_cell(i) == np1) lnode(i-i1+1) = 1
      if (nodes_of_cell(i) == np2) lnode(i-i1+1) = 1
    end do
    !
    if (sum(lnode(1:2)) == 2) then
      ifa = 1
    else if (sum(lnode(2:3)) == 2) then
      ifa = 2
    else if (sum(lnode(3:4)) == 2) then
      ifa = 3
    else
      ifa = merge(3,4,ip==3)
    end if
    !
    ! 2. The right side of the face
    !
    j1 = nodes_of_cell_ptr(jc)+1
    j2 = nodes_of_cell_ptr(jc+1)
    jp = j2-j1+1
    !
    ! Find which face of the right cell is on the boundary
    !
    lnode(:) = 0
    do j = j1,j2
      if (nodes_of_cell(j) == np1) lnode(j-j1+1) = 1
      if (nodes_of_cell(j) == np2) lnode(j-j1+1) = 1
    end do
    !
    if (sum(lnode(1:2)) == 2) then
      jfa = 1
    else if (sum(lnode(2:3)) == 2) then
      jfa = 2
    else if (sum(lnode(3:4)) == 2) then
      jfa = 3
    else
      jfa = merge(3,4,jp==3)
    end if
    !
    ! 3. Now add this information the boundary face array
    !
    ! Left side of the face
    !
    bface(1,nf) = bc_cpu_bnd       ! boundary condition
    bface(2,nf) = ic               ! host cell
    bface(3,nf) = merge(3,2,ip==3) ! host cell geometry type
    bface(4,nf) = ifa              ! face of host cell on boundary
    bface(5,nf) = nf + 1           ! adjacent partition boundary face
    bface(6,nf) = jc               ! adjacent partition host cell
    bface(7,nf) = epart(jc) - 1    ! adjacent partition processor ID
    bface(8,nf) = merge(3,2,jp==3) ! adjacent partition host cell geometry type
    bface(9,nf) = jfa              ! face of adjacent partition host cell on boundary
    !
    nf = nf + 1
    !
    ! Right side of the face
    !
    bface(1,nf) = bc_cpu_bnd       ! boundary condition
    bface(2,nf) = jc               ! host cell
    bface(3,nf) = merge(3,2,jp==3) ! host cell geometry type
    bface(4,nf) = jfa              ! face of host cell on boundary
    bface(5,nf) = nf - 1           ! adjacent partition boundary face
    bface(6,nf) = ic               ! adjacent partition host cell
    bface(7,nf) = epart(ic) - 1    ! adjacent partition processor ID
    bface(8,nf) = merge(3,2,ip==3) ! adjacent partition host cell geometry type
    bface(9,nf) = ifa              ! face of adjacent partition host cell on boundary
    !
    nf = nf + 1
    !
  end do
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Allocate a temporary array to resize and reorder cells_on_face
  !
  allocate ( temp_array(1:nifac,1:new_ngface) , source=0 , stat=ierr )
  call alloc_error("find_partition_boundaries","temp_array",1,ierr)
  !
  ! Redo the boundary face section of cells_on_face
  ! because of the new additions to bface
  !
  do nf = 1,new_ngfbnd
    !
    n = bface(2,nf)
    !
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    !
    np = bface(4,nf) + n1 - 1
    !
    temp_array(1,nf) = bface(2,nf)
    temp_array(2,nf) = ngcell + nf
    temp_array(3,nf) = nodes_of_cell(np)
    temp_array(4,nf) = nodes_of_cell( merge(np+1,n1,np/=n2) )
    !
  end do
  !
  ! Copy the interior face information to the temporary array
  ! but skip those interior faces that are now boundary faces
  !
  nf = new_ngfbnd + 1
  do n = ngfbnd+1,ngface
    if ( msk(n) ) cycle
    temp_array(:,nf) = cells_on_face(:,n)
    nf = nf + 1
  end do
  !
  ! Use move_alloc to reallocate cells_on_face. This does the following
  !     1. Deallocates cells_on_face
  !     2. Allocates cells_on_face with the same type, type
  !        parameters, array bounds, and value as temp_array
  !     3. Deallocates temp_array
  !
  call move_alloc( from=temp_array , to=cells_on_face )
  !
  ! We need to update the global grid parameters
  !
  ngfbnd = new_ngfbnd
  ngface = new_ngface
  !
  ! Deallocate msk since we are finished with it
  !
  deallocate ( msk , stat=ierr )
  call alloc_error("find_partition_boundaries","msk",2,ierr)
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
  ! Perform one last quick test to make sure that the attached cells
  ! for each interior face are located within the same partition
  !
  if ( any( epart(cells_on_face(1,ngfbnd+1:ngface)) /= &
            epart(cells_on_face(2,ngfbnd+1:ngface)) ) ) then
    write (iout,1)
    call stop_recon(stop_mpi,"find_partition_boundaries")
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
subroutine create_grid_mappings(npart,epart,cells_on_face, &
                                nodes_of_cell,nodes_of_cell_ptr, &
                                cell_map,cell_map_ptr,cell_rmap, &
                                node_map,node_map_ptr,node_rmap, &
                                face_map,face_map_ptr,face_rmap)
  !
  !.. Formal Arguments ..
  integer,                                 intent(in) :: npart
  integer,              dimension(:),      intent(in) :: epart
  integer,              dimension(:),      intent(in) :: nodes_of_cell
  integer,              dimension(:),      intent(in) :: nodes_of_cell_ptr
  integer,              dimension(:,:),    intent(in) :: cells_on_face
  integer, allocatable, dimension(:),   intent(inout) :: cell_map
  integer, allocatable, dimension(:),   intent(inout) :: cell_map_ptr
  integer, allocatable, dimension(:,:), intent(inout) :: cell_rmap
  integer, allocatable, dimension(:),   intent(inout) :: node_map
  integer, allocatable, dimension(:),   intent(inout) :: node_map_ptr
  integer, allocatable, dimension(:,:), intent(inout) :: node_rmap
  integer, allocatable, dimension(:),   intent(inout) :: face_map
  integer, allocatable, dimension(:),   intent(inout) :: face_map_ptr
  integer, allocatable, dimension(:,:), intent(inout) :: face_rmap
  !
  !.. Local Scalars ..
  integer :: i,n,nc,n1,n2,nf,ierr
  !
  !.. Local Allocatable Arrays ..
  logical, allocatable, dimension(:) :: msk
  !
continue
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Create the local to global cell mappings
  ! NOTE: No need to preallocate cell_map because it will
  !       be allocated upon the first call to reallocate.
  !
  allocate ( cell_map_ptr(1:npart+1) , source=0 , stat=ierr )
  call alloc_error("create_grid_mappings","cell_map_ptr",1,ierr)
  allocate ( cell_rmap(1:ngcell,1:npart) , source=0 , stat=ierr )
  call alloc_error("create_grid_mappings","cell_rmap",1,ierr)
  !
  ! Loop through the partitions to find the local number
  ! of cells for each partition and dynamically create the
  ! cell mappings instead of having to do this loop twice.
  !
  do n = 2,size(cell_map_ptr)
    !
    ! Count the cells on this partition
    ! NOTE: partition being worked on is actually partition 'n-1'
    !
    cell_map_ptr(n) = cell_map_ptr(n-1) + count( epart == (n-1) )
    !
    ! The extents of cell_map for the current partition
    !
    n1 = cell_map_ptr(n-1)+1
    n2 = cell_map_ptr(n)
    !
    ! Expand cell_map to make room for the cells of the current partition
    !
    call reallocate(cell_map,n2)
    !
    ! Add the cells of the current partition to cell_map
    !
    cell_map(n1:n2) = pack( [1:ngcell] , epart(1:ngcell) == (n-1) )
   !cell_map = [cell_map,find_subset(epart,n-1)]
    !
  end do
  !
  ! Create the reverse maps for cells
  !
  cell_rmap(:,:) = create_reverse_map(ngcell,cell_map,cell_map_ptr)
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Create the local to global node mappings
  ! NOTE: No need to preallocate node_map because it will
  !       be allocated upon the first call to reallocate.
  !
  allocate ( node_map_ptr(1:npart+1) , source=0 , stat=ierr )
  call alloc_error("create_grid_mappings","node_map_ptr",1,ierr)
  allocate ( node_rmap(1:ngnode,1:npart) , source=0 , stat=ierr )
  call alloc_error("create_grid_mappings","node_rmap",1,ierr)
  !
  ! Loop through the partitions to find the local number of
  ! nodes for each partition and dynamically create the node
  ! mappings instead of having to do this double loop twice.
  !
  allocate ( msk(1:ngnode) , source=fals , stat=ierr )
  call alloc_error("create_grid_mappings","msk",1,ierr)
  !
  do n = 2,size(node_map_ptr)
    !
    ! Initialize/Reset msk
    !
    msk(:) = fals
    !
    ! Mark the nodes that are within this partition
    ! NOTE: partition being worked on is actually partition 'n-1'
    !
    do nc = 1,ngcell
      if (epart(nc) /= n-1) cycle
      n1 = nodes_of_cell_ptr(nc)+1
      n2 = nodes_of_cell_ptr(nc+1)
      msk(nodes_of_cell(n1:n2)) = true
    end do
    !
    ! Count the nodes on this partition
    !
    node_map_ptr(n) = node_map_ptr(n-1) + count( msk(1:ngnode) )
    !
    ! The extents of node_map for the current partition
    !
    n1 = node_map_ptr(n-1)+1
    n2 = node_map_ptr(n)
    !
    ! Expand node_map to make room for the nodes of the current partition
    !
    call reallocate(node_map,n2)
    !
    ! Add the nodes of the current partition to node_map
    !
    node_map(n1:n2) = pack( [1:ngnode] , msk(1:ngnode) )
    !
  end do
  !
  deallocate ( msk , stat=ierr )
  call alloc_error("create_grid_mappings","msk",2,ierr)
  !
  ! Create the reverse maps for nodes
  !
  node_rmap(:,:) = create_reverse_map(ngnode,node_map,node_map_ptr)
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Create the local to global face mappings
  ! NOTE: No need to preallocate face_map because it will
  !       be allocated upon the first call to reallocate.
  !
  allocate ( face_map_ptr(1:npart+1) , source=0 , stat=ierr )
  call alloc_error("create_grid_mappings","face_map_ptr",1,ierr)
  allocate ( face_rmap(1:ngface,1:npart) , source=0 , stat=ierr )
  call alloc_error("create_grid_mappings","face_rmap",1,ierr)
  !
  ! Loop through the partitions to find the local number of
  ! faces for each partition and dynamically create the face
  ! mappings instead of having to do this double loop twice.
  !
  do n = 2,size(face_map_ptr)
    !
    ! Count the faces that are within this partition
    ! NOTE: partition being worked on is actually partition 'n-1'
    !
    face_map_ptr(n) = face_map_ptr(n-1) + &
                      count( epart(cells_on_face(1,1:ngface)) == (n-1) )
    !
    ! The extents of face_map for the current partition
    !
    n1 = face_map_ptr(n-1)+1
    n2 = face_map_ptr(n)
    !
    ! Expand face_map to make room for the faces of the current partition
    !
    call reallocate(face_map,n2)
    !
    ! Add the faces of the current partition to face_map
    !
    face_map(n1:n2) = pack( [1:ngface] , &
                            epart(cells_on_face(1,1:ngface)) == (n-1) )
   !face_map = [face_map,find_subset( epart(cells_on_face(1,1:ngface)) , n-1 )]
    !
  end do
  !
  ! Create the reverse maps for faces
  !
  face_rmap(:,:) = create_reverse_map(ngface,face_map,face_map_ptr)
  !
end subroutine create_grid_mappings
!
!###############################################################################
!
pure function create_reverse_map(nsize,map,map_ptr) result(return_value)
  !
  !.. Formal Arguments ..
  integer,               intent(in) :: nsize
  integer, dimension(:), intent(in) :: map
  integer, dimension(:), intent(in) :: map_ptr
  !
  !.. Function Return Value ..
  integer, allocatable, dimension(:,:) :: return_value
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,np,ierr
  !
continue
  !
  np = size(map_ptr)-1
  !
  ! Allocate the function result
  !
  allocate ( return_value(1:nsize,1:np) , source=0 , stat=ierr )
 !call alloc_error("create_reverse_map","return_value",1,ierr)
  !
  ! Create the reverse mapping array for each partition
  !
  do n = 1,np
    n1 = map_ptr(n)+1
    n2 = map_ptr(n+1)
    return_value(map(n1:n2),n) = [1:n2-n1+1]
  end do
  !
end function create_reverse_map
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
  integer :: i
  !
continue
  !
  return_value = pack( [(i,i=1,size(array))] , array == find_value )
  !
end function find_subset
!
!###############################################################################
!
subroutine localize_grid(cell_map,cell_map_ptr,cell_rmap, &
                         node_map,node_map_ptr,node_rmap, &
                         face_map,face_map_ptr,face_rmap, &
                         nodes_of_cell,nodes_of_cell_ptr, &
                         bface,xyz_nodes)
  !
  !.. Formal Arguments ..
  integer,                dimension(:),      intent(in) :: cell_map
  integer,                dimension(:),      intent(in) :: cell_map_ptr
  integer,                dimension(:,:),    intent(in) :: cell_rmap
  integer,                dimension(:),      intent(in) :: node_map
  integer,                dimension(:),      intent(in) :: node_map_ptr
  integer,                dimension(:,:),    intent(in) :: node_rmap
  integer,                dimension(:),      intent(in) :: face_map
  integer,                dimension(:),      intent(in) :: face_map_ptr
  integer,                dimension(:,:),    intent(in) :: face_rmap
  integer,   allocatable, dimension(:),   intent(inout) :: nodes_of_cell
  integer,   allocatable, dimension(:),   intent(inout) :: nodes_of_cell_ptr
  integer,   allocatable, dimension(:,:), intent(inout) :: bface
  real(vbp), allocatable, dimension(:,:), intent(inout) :: xyz_nodes
  !
  !.. Local Scalars ..
  integer :: n,nf,nc,np,n1,n2,g1,g2,gc,gp
  integer :: np1,np2,nc1,nc2,nf1,nf2
  integer :: nfbnd,ncell,nnode,nface
  integer :: n_totcel,n_totnod,nmin,nmax,gn,rev_gn,gbf
  integer :: inbeg,inend,igbeg,igend
  integer :: lhc,lhg,lgg,lfa,lgc,lh1,lg1,gfa
  integer :: my_part,ierr
  !
  !.. Local Arrays ..
  integer, dimension(1:1) :: li
  !
  !.. Local Allocatable Arrays ..
  logical,   allocatable, dimension(:)   :: msk
  integer,   allocatable, dimension(:)   :: recv_blklen
  integer,   allocatable, dimension(:)   :: recv_offset
  integer,   allocatable, dimension(:)   :: cell_offset
  integer,   allocatable, dimension(:)   :: ghostcell_map
  integer,   allocatable, dimension(:)   :: ghostnode_map
  integer,   allocatable, dimension(:)   :: ghostnode_rmap
  integer,   allocatable, dimension(:)   :: global_nofc
  integer,   allocatable, dimension(:)   :: global_nofc_ptr
  integer,   allocatable, dimension(:,:) :: global_bface
  real(vbp), allocatable, dimension(:,:) :: global_xyz
  !
continue
  !
  ! We need to make copies of the global grid arrays nodes_of_cell,
  ! nodes_of_cell_ptr, bface, and xyz_nodes before we localize them.
  ! Since we need to reallocate the original arrays to the new local
  ! grid sizes, use the intrinsic subroutine move_alloc which does:
  !     1. Deallocates 'to=' array if it is allocated
  !     2. Allocates 'to=' array with the same type, type parameters,
  !        array bounds, and with the same value as 'from=' array
  !     3. Deallocates 'from=' array
  !
  call move_alloc( from=nodes_of_cell     , to=global_nofc )
  call move_alloc( from=nodes_of_cell_ptr , to=global_nofc_ptr )
  call move_alloc( from=bface             , to=global_bface )
  call move_alloc( from=xyz_nodes         , to=global_xyz )
  !
  ! The process ID, mypnum, starts numbering at 0 so
  ! identify which partition is local to each procesor
  !
  my_part = mypnum + 1
  !
  ! Get the offset indices for this processor
  !
  nc1 = cell_map_ptr(my_part)
  nc2 = cell_map_ptr(my_part+1)
  !
  np1 = node_map_ptr(my_part)
  np2 = node_map_ptr(my_part+1)
  !
  nf1 = face_map_ptr(my_part)
  nf2 = face_map_ptr(my_part+1)
  !
  ! Count the local number of cells, nodes, faces, and boundary faces
  !
  ncell = nc2 - nc1
  nnode = np2 - np1
  nface = nf2 - nf1
  nfbnd = count( face_rmap(1:ngfbnd,my_part) /= 0 )
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Localize the boundary face array, but not the information it contains.
  !
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , stat=ierr )
  call alloc_error("localize_grid","bface",1,ierr)
  allocate ( ghostcell_map(1:ngtotcel) , source=0 , stat=ierr )
  call alloc_error("localize_grid","ghostcell_map",1,ierr)
  !
  ! First extract the information that is local to the current
  ! processor from the global_bface array to the bface array.
  ! Also create a mapping array to identify the global ghost cells
  ! needed for the local grid.
  !
  nc = ncell
  do nf = 1,nfbnd
    !
    ! Global face ID for this local boundary face
    !
    gbf = face_map(nf1+nf)
    !
    ! Copy everything to the local bface array
    !
    bface(:,nf) = global_bface(:,gbf)
    !
    ! Now localize the information for this face
    !
    bface(2,nf) = cell_rmap( bface(2,nf) , my_part ) ! local host cell
    !
    if ( any(bface(1,nf) == bc_comm) ) then
      bface(5,nf) = face_rmap( bface(5,nf) , bface(7,nf)+1 ) ! connecting boundary face
      bface(6,nf) = cell_rmap( bface(6,nf) , bface(7,nf)+1 ) ! connecting host cell
    end if
    !
    ! The global cell ID for the local ghost cell on this local boundary face.
    ! If this boundary face is a non-periodic communication boundary, the ghost
    ! cell corresponds to an interior cell of the global grid. Otherwise, the
    ! ghost cell corresponds to a ghost cell of the global grid.
    !
    gc = merge( ngcell+gbf , global_bface(6,gbf) , bface(1,nf) /= bc_cpu_bnd )
    !
    ! Create a mapping array, ghostcell_map, just for the ghost cells of the
    ! local grid. If ghostcell_map is given the global cell ID for the local
    ! ghost cell, it will return the local ghost cell ID. Also store the
    ! local ghost cell ID in the bface(10,:) location.
    ! EXPLANATION: If no communication boundaries existed, none of this would
    !              be necessary since each boundary face should correspond
    !              to its own unique ghost cell. However, when communication
    !              boundaries exist, more than one communication boundary on
    !              the same partition can correspond to the same cell on an
    !              adjacent partition. Therefore, the local ghost cell for
    !              these multiple communication boundaries is shared and thus
    !              not unique. The following if statement ensures that only a
    !              single ghost cell exists if shared by multiple boundary
    !              faces.
    !
    if (ghostcell_map(gc) /= 0) then
      bface(10,nf) = ghostcell_map(gc)
    else
      nc = nc + 1
      ghostcell_map(gc) = nc
      bface(10,nf) = nc
    end if
    !
  end do
  !
  n_totcel = nc
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Localize the node connectivity pointer for each cell
  !
  ! Create the local nodes_of_cell_ptr
  !
  allocate ( nodes_of_cell_ptr(1:n_totcel+1) , source=0 , stat=ierr )
  call alloc_error("localize_grid","nodes_of_cell_ptr",1,ierr)
  !
  ! First set the node connectivity pointer to the number
  ! of nodes that define each of the interior cells
  ! NOTE: Instead of initially creating this pointer array in a commulative
  !       manner as it will be once it is finished, we will just set the
  !       array to the number of nodes for each cell. Once this is finished
  !       we will then reshuffle the array to be commulative. This is done
  !       because of the ghost cell loop below which relies on the value nc
  !       returned by ghostcell_map to update the pointer corresponding to
  !       the ghost cell nc. These values returned by ghostcell_map are
  !       generally not ordered consecutively and therfore the commulative
  !       pointer nc-1 is not guaranteed to be defined yet.
  !
  do nc = 1,ncell
    gc = cell_map(nc1+nc)
    nodes_of_cell_ptr(nc+1) = global_nofc_ptr(gc+1) - global_nofc_ptr(gc)
  end do
  !
  ! Next set the node connectivity pointer to the number
  ! of nodes that define each of the ghost cells
  !
  do gc = 1,ngtotcel
    if (ghostcell_map(gc) == 0) cycle
    nc = ghostcell_map(gc)
    nodes_of_cell_ptr(nc+1) = global_nofc_ptr(gc+1) - global_nofc_ptr(gc)
  end do
  !
  ! Reshuffle the ghost cell section of nodes_of_cell_ptr
  ! to make the pointer commulative 
  !
  do nc = 2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(nc) = nodes_of_cell_ptr(nc) + nodes_of_cell_ptr(nc-1)
  end do
  !
  ! We are now finished with ghostcell_map so deallocate
  !
  deallocate ( ghostcell_map , stat=ierr )
  call alloc_error("localize_grid","ghostcell_map",2,ierr)
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Localize the node connectivity for each cell
  !
  allocate ( nodes_of_cell(1:last(nodes_of_cell_ptr)) , source=0 , stat=ierr )
  call alloc_error("localize_grid","nodes_of_cell",1,ierr)
  !
  ! First extract the information that is local to the current
  ! processor from the global_nofc array to the nodes_of_cell array.
  !
  do nc = 1,ncell
    !
    n1 = nodes_of_cell_ptr(nc)+1
    n2 = nodes_of_cell_ptr(nc+1)
    !
    gc = cell_map(nc1+nc)
    !
    g1 = global_nofc_ptr(gc)+1
    g2 = global_nofc_ptr(gc+1)
    !
    nodes_of_cell(n1:n2) = global_nofc(g1:g2)
    !
  end do
  !
  ! Now go through the boundary faces and extract the nodal information for
  ! the ghost cells from the global_nofc array to the nodes_of_cell array.
  ! NOTE: This loop is made more complicated because we require that the
  !       first two nodes of the ghost cell are the two nodes that lie on
  !       the boundary face. However, if a ghost cell is connected to
  !       multiple boundary faces, the first two nodes of the ghost cell
  !       are the two nodes on the face for only one boundary face.
  !
  allocate ( ghost_face(ncell+1:n_totcel) , source=0 , stat=ierr )
  call alloc_error("localize_grid","ghost_face",1,ierr)
  !
  do nf = 1,nfbnd
    !
    gbf = face_map(nf1+nf) ! global boundary face
    !
    ! On a non-periodic communication boundary, the global ghost
    ! cell is the host cell on the adjacent partition for the
    ! current boundary face.  For all other boundary conditions,
    ! the global ghost cell is just a ghost cell from the global grid.
    !
    gc = merge( ngcell+gbf , global_bface(6,gbf) , bface(1,nf) /= bc_cpu_bnd )
    !
    g1 = global_nofc_ptr(gc)+1 ! beginning node index for global ghost cell
    g2 = global_nofc_ptr(gc+1) ! ending    node index for global ghost cell
    !
    lhc = bface(2,nf)               ! local host cell
    lhg = merge(3,4,bface(3,nf)==3) ! geometry of local host cell
    lgg = merge(3,4,bface(8,nf)==3) ! geometry of local ghost cell
    lfa = bface(4,nf)               ! face of local host cell on boundary
    lgc = bface(10,nf)              ! local ghost cell
    !
    lh1 = nodes_of_cell_ptr(lhc)+1 ! beginning node index for local host  cell
    lg1 = nodes_of_cell_ptr(lgc)+1 ! beginning node index for local ghost cell
    !
    ! Make the first two nodes for the ghost cell
    ! the two nodes that lie on the boundary face
    !
    nodes_of_cell(lg1+0) = nodes_of_cell(lh1+mod(lfa,lhg))
    nodes_of_cell(lg1+1) = nodes_of_cell(lh1+lfa-1)
    !
    ! Find the offset from g1 of the global node that corresponds to
    ! the local ghost cell node lg1. For non-periodic communication
    ! boundaries, bface(9,nf) corresponds to the local ghost cell node
    ! lg1+1. For periodic communication boundaries, the value bface(9,nf)
    ! gives the face of the opposite periodic host cell which does not
    ! necessarily align with the face of the global ghost cell. Therfore,
    ! we need to match the local ghost cell nodes with those of the global
    ! ghost cell nodes. For all other boundaries, there are no alignment
    ! issues and we just need to match the global ghost cell nodes to the
    ! local ghost cell nodes. Doing this for non-periodic communication
    ! boundaries will result in the same results as if using the rotation
    ! value found in bface(9,nf) so just do the same for all boundaries.
    ! NOTE: Because g1 is the first node index and li(1) is a value from 1
    !       to 3 (for triangles) or 4 (for quads), li(1) is actually the
    !       offset+1 and therfore the value g1+li(1) actually corresponds
    !       to the local ghost cell node lg1+1.
    !
    li = pack( [1:lgg] , nodes_of_cell(lg1) == global_nofc(g1:g2) )
   !write (*,'(7i5)') mypnum,nf,lfa,bface(9,nf),li(1),lhc,lgc
    !
    if (nodes_of_cell(lg1+1) /= global_nofc(g1+mod(li(1),lgg))) then
      write (*,*) "Error with ghost cell node ordering ", &
                  "in subroutine localize grid!"
      call stop_recon(abort,"localize_grid")
    end if
    !
    ! Now copy the remaining nodes from the global ghost cell
    ! to the local ghost cell while maintaining node ordering
    !
    do n = 1,lgg-2
      nodes_of_cell(lg1+1+n) = global_nofc(g1+mod(li(1)+n,lgg))
    end do
    !
    ! While we have the information available, we need to store the
    ! orientation of the ghost cell in relation to its corresponding
    ! global cell. If the first two nodes of the ghost cell and global
    ! cell match then ghost_face(ghost_cell)=1. If the first two nodes
    ! of the ghost cell match the second and third nodes of the global
    ! cell then ghost_face(ghost_cell)=2. And so on... This is the same
    ! as the face of the adjacent cell that is on the boundary.
    ! NOTE: When finding the value li above, the value in bface(9,nf) was
    !       not valid for periodic communication boundaries because the
    !       alignment of the ghost cell nodes and the opposite periodic
    !       cell are not necessarily the same. However, when we are
    !       looking to communicate the information from the opposite
    !       periodic cell to the local ghost cell, the value in bface(9,nf)
    !       does give the correct orientation between those two cells.
    !
    ghost_face(lgc) = bface(9,nf)
    !
  end do
  !
  ! We no longer need global_bface or global_nofc so deallocate
  !
  deallocate ( global_bface , stat=ierr )
  call alloc_error("localize_grid","global_bface",2,ierr)
  deallocate ( global_nofc , stat=ierr )
  call alloc_error("localize_grid","global_nofc",2,ierr)
  !
  ! #######
  ! ## 4 ##
  ! #######
  !
  ! We need to create collective receiving MPI datatypes to collect
  ! localized data from each processor into a shared global form.
  !
  ! First, we need to create an array containing the
  ! solution point offsets for every global cell
  !
  allocate ( cell_offset(1:ngcell) , source=0 , stat=ierr )
  call alloc_error("localize_grid","cell_offset",1,ierr)
  !
  do n = 1,ngcell-1
    np = global_nofc_ptr(n+1) - global_nofc_ptr(n)
    cell_offset(n+1) = cell_offset(n) + merge(n_sp2dt,n_sp2dq,np==3)
  end do
  !
  ! Loop over all the processors to create the collective receive MPI datatypes
  !
  n = size(cell_map_ptr)
  np = maxval( cell_map_ptr(2:size(cell_map_ptr)  ) - &
               cell_map_ptr(1:size(cell_map_ptr)-1) )
  allocate ( recv_offset(1:np) , source=0 , stat=ierr )
  call alloc_error("localize_grid","recv_offset",1,ierr)
  allocate ( recv_blklen(1:np) , source=0 , stat=ierr )
  call alloc_error("localize_grid","recv_blklen",1,ierr)
  allocate ( recv_xyz(1:ncpu) , source=mpi_datatype_null , stat=ierr )
  call alloc_error("localize_grid","recv_xyz",1,ierr)
  allocate ( recv_sol(1:ncpu) , source=mpi_datatype_null , stat=ierr )
  call alloc_error("localize_grid","recv_sol",1,ierr)
  !
  do n = 1,ncpu
    !
    n1 = cell_map_ptr(n)
    n2 = cell_map_ptr(n+1)
    np = n2-n1
    !
    do nc = 1,np
      !
      gc = cell_map(n1+nc)
      gp = global_nofc_ptr(gc+1) - global_nofc_ptr(gc)
      !
      recv_offset(nc) = cell_offset(gc)
      recv_blklen(nc) = merge(n_sp2dt,n_sp2dq,gp==3)
      !
    end do
    !
    ! Now create the receive MPI datatype for the local
    ! solution point coordinates on the current processor
    !
    recv_offset(1:np) = recv_offset(1:np)*nr ! adjust the offset for nr dimensions
    recv_blklen(1:np) = recv_blklen(1:np)*nr ! adjust the offset for nr dimensions
    !
    call mpi_type_indexed(np,recv_blklen(1:np),recv_offset(1:np), &
                          mpi_flttyp,recv_xyz(n),mpierr)
    call mpi_type_commit(recv_xyz(n),mpierr)
    !
    recv_offset(1:np) = recv_offset(1:np)/nr ! cancel out the nr adjustment above
    recv_blklen(1:np) = recv_blklen(1:np)/nr ! cancel out the nr adjustment above
    !
    ! Now create the receive MPI datatype for the local
    ! solution at the solution points on the current processor
    !
    recv_offset(1:np) = recv_offset(1:np)*nq ! adjust the offset for nq flow equations
    recv_blklen(1:np) = recv_blklen(1:np)*nq ! adjust the offset for nq flow equations
    !
    call mpi_type_indexed(np,recv_blklen(1:np),recv_offset(1:np), &
                          mpi_flttyp,recv_sol(n),mpierr)
    call mpi_type_commit(recv_sol(n),mpierr)
    !
  end do
  !
  ! We no longer need global_nofc_ptr, recv_offset,
  ! recv_blklen, or cell_offset so deallocate
  !
  deallocate ( global_nofc_ptr , stat=ierr )
  call alloc_error("localize_grid","global_nofc_ptr",2,ierr)
  deallocate ( recv_offset , stat=ierr )
  call alloc_error("localize_grid","recv_offset",2,ierr)
  deallocate ( recv_blklen , stat=ierr )
  call alloc_error("localize_grid","recv_blklen",2,ierr)
  deallocate ( cell_offset , stat=ierr )
  call alloc_error("localize_grid","cell_offset",2,ierr)
  !
  ! #######
  ! ## 5 ##
  ! #######
  !
  inbeg = 1                          ! beginning index of interior cell nodes
  inend = nodes_of_cell_ptr(ncell+1) ! ending    index of interior cell nodes
  igbeg = inend + 1                  ! beginning index of ghost    cell nodes
  igend = size(nodes_of_cell)        ! ending    index of ghost    cell nodes
  !
  ! Reverse map the interior nodes to the local ordering
  !
  do n = inbeg,inend
    nodes_of_cell(n) = node_rmap( nodes_of_cell(n) , my_part )
  end do
  !
  ! Now go through the ghost cells and reverse map those nodes that
  ! are interior nodes of the global grid and mark the remaining
  ! nodes by making them the negative of their global node ID
  !
  do n = igbeg,igend
    !
    gn = nodes_of_cell(n)
    if (gn <= ngnode) then
      rev_gn = node_rmap( gn , my_part )
      nodes_of_cell(n) = merge( rev_gn , -gn , rev_gn /= 0 )
    else
      nodes_of_cell(n) = -gn
    end if
    !
  end do
  !
  ! Now find the negative node values and add them to the total number of nodes
  !
  nmin = abs(maxval(nodes_of_cell, mask=nodes_of_cell<0))
  nmax = abs(minval(nodes_of_cell, mask=nodes_of_cell<0))
  !
  allocate ( msk(nmin:nmax) , source=fals , stat=ierr )
  call alloc_error("localize_grid","msk",1,ierr)
  !
  ! Create a mask array for all the negative nodes in nodes_of_cell
  !
  where (nodes_of_cell < 0) msk(abs(nodes_of_cell)) = true
  !
  n_totnod = nnode + count(msk) ! total number of nodes for the local partition
  !
  ! Create the forward and reverse mapping arrays for the new local ghost nodes
  !
  allocate ( ghostnode_map(nnode+1:n_totnod) , source=0 , stat=ierr )
  call alloc_error("localize_grid","ghostnode_map",1,ierr)
  allocate ( ghostnode_rmap(nmin:nmax) , source=0 , stat=ierr )
  call alloc_error("localize_grid","ghostnode_rmap",1,ierr)
  !
  ghostnode_map(:) = pack( [nmin:nmax] , msk )          ! forward map
  ghostnode_rmap(ghostnode_map(:)) = [nnode+1:n_totnod] ! reverse map
  !
  ! Now loop through and replace the negative node
  ! values with the new localized ghost node values
  !
  where (nodes_of_cell < 0) nodes_of_cell = ghostnode_rmap( abs(nodes_of_cell) )
  !
  ! We no longer need ghostnode_rmap or msk so deallocate
  !
  deallocate ( msk , stat=ierr )
  call alloc_error("localize_grid","msk",2,ierr)
  deallocate ( ghostnode_rmap , stat=ierr )
  call alloc_error("localize_grid","ghostnode_rmap",2,ierr)
  !
  ! #######
  ! ## 6 ##
  ! #######
  !
  ! Localize the node coordinates
  !
  allocate ( xyz_nodes(1:nr,1:n_totnod) , source=zero , stat=ierr )
  call alloc_error("localize_grid","xyz_nodes",1,ierr)
  !
  ! Localize the interior nodes
  !
  do np = 1,nnode
    xyz_nodes(:,np) = global_xyz(:,node_map(np1+np))
  end do
  !
  ! Localize the remaining nodes of the ghost
  ! cells that are not on the boundary faces
  !
  do np = nnode+1,n_totnod
    xyz_nodes(:,np) = global_xyz(:,ghostnode_map(np))
  end do
  !
  ! We no longer need ghostnode_map or global_xyz so deallocate
  !
  deallocate ( ghostnode_map , stat=ierr )
  call alloc_error("localize_grid","ghostnode_map",2,ierr)
  deallocate ( global_xyz , stat=ierr )
  call alloc_error("localize_grid","global_xyz",2,ierr)
  !
  ! #######
  ! ## 7 ##
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
 !! ## 8 ##
 !! #######
 !!
 !! Output the individual partitions to separate tecplot files
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
! WHEN CREATING THE EXCHANGE DATATYPES, A POSSIBLE METHOD WOULD BE TO USE MPI
! COLLECTIVE COMMUNICATION TO GATHER IDENTICAL ARRAYS ACROSS ALL PROCESSORS
! CONTAINING CONCATENATED COPIES OF THE ARRAYS LOCAL TO EACH PROCESSOR.
! E.G., GATHER A GLOBAL BFACE ARRAY, gbface, ON ALL PROCESSORS CONTAINING
! gbface = [bface from cpu 0,bface from cpu 1,...,bface from cpu ncpu]
!
!###############################################################################
!
subroutine create_nonblocking_exchange_datatypes(nodes_of_cell, &
                                                 nodes_of_cell_ptr,bface)
  !
  !.. Formal Arguments ..
  integer, dimension(:),   intent(in) :: nodes_of_cell
  integer, dimension(:),   intent(in) :: nodes_of_cell_ptr
  integer, dimension(:,:), intent(in) :: bface
  !
  !.. Local Scalars ..
  integer :: n,nf,np,scnt,rcnt,i,icpu,ierr
  integer :: smin,smax,rmin,rmax,s1,s2,r1,r2
  integer :: p1,p2,nc,rfa,rge
  integer :: nfbnd,ncell,n_totcel,n_totpts
  !
  !.. Local Arrays ..
  integer, dimension(1:ncpu) :: cpu_id
  !
  !.. Local Allocatable Arrays ..
  logical, allocatable, dimension(:) :: send_msk
  integer, allocatable, dimension(:) :: recv_map
  integer, allocatable, dimension(:) :: send_temp
  integer, allocatable, dimension(:) :: recv_temp
  integer, allocatable, dimension(:) :: cell_offset
  !
continue
  !
  nfbnd = size(bface,dim=2)
  ncell = minval(bface(10,:)) - 1
  n_totcel = size(nodes_of_cell_ptr)-1
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
 !write (150+mypnum,*) "   bc    hc  htype  hloc  adjf  adjh  adjp  adjt  adjl   gc"
 !do nf = 1,nfbnd
 !  write (150+mypnum,'(*(i6))') (bface(n,nf),n=1,size(bface,1))
 !end do
  ! Create the offset for each cell of the local
  ! grid due to the solution points within the cell
  !
  allocate ( cell_offset(1:n_totcel) , source=0 , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","cell_offset",1,ierr)
  !
  do n = 1,n_totcel-1
    np = nodes_of_cell_ptr(n+1) - nodes_of_cell_ptr(n)
    cell_offset(n+1) = cell_offset(n) + merge(n_sp2dt,n_sp2dq,np==3)
  end do
  np = nodes_of_cell_ptr(n_totcel+1) - nodes_of_cell_ptr(n_totcel)
  n_totpts = last(cell_offset) + merge(n_sp2dt,n_sp2dq,np==3)
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Now onto creating the communication datatypes
  !
  ! First we need to find which processes the
  ! current process needs to exchange with
  !
  cpu_id(:) = -99
  do nf = 1,nfbnd
    if ( any(bface(1,nf) == bc_comm) ) cpu_id(bface(7,nf)+1) = bface(7,nf)
  end do
  !
  ! Make sure that the current process doesnt try to communicate with itself
  !
  cpu_id(mypnum+1) = -99
  !
  ! Count the number of processes the current process needs to communicate with
  !
  ncomm = count(cpu_id >= 0)
  !
  ! Create a condensed array containing the process
  ! IDs the current process needs to communicate with
  !
  allocate ( exch_cpu(1:ncomm) , source=pack(cpu_id,cpu_id>=0) , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","exch_cpu",1,ierr)
  !
  ! Now loop through the processes to exchange with
  ! and create the send and receive cell maps
  !
  rmin = minval(bface(6,:))
  rmax = maxval(bface(6,:))
  smin = minval(bface(2,:))
  smax = maxval(bface(2,:))
  allocate ( send_msk(smin:smax) , source=fals , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","send_msk",1,ierr)
  allocate ( recv_map(rmin:rmax) , source=0 , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","recv_map",1,ierr)
  allocate ( send_temp(1:smax-smin+1) , source=0 , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","send_temp",1,ierr)
  allocate ( recv_temp(1:rmax-rmin+1) , source=0 , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","recv_temp",1,ierr)
  allocate ( send_cell_ptr(1:ncomm+1) , source=0 , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","send_cell_ptr",1,ierr)
  allocate ( recv_cell_ptr(1:ncomm+1) , source=0 , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","recv_cell_ptr",1,ierr)
  !
  do i = 1,ncomm
    !
    icpu = exch_cpu(i)
    !
    ! Loop to find the communication boundary faces between the
    ! current processor and the current exchange partner, icpu.
    ! We will create a mask array, send_msk, to identify all the
    ! local host cells that need to be sent to the ghost cells
    ! of icpu. We will also create an array, recv_map, that maps
    ! the host cells of icpu to the local ghost cells whose data
    ! needs to be received from icpu.
    !
    ! #################################################################
    ! #################################################################
    ! #################################################################
    !            WARNING! WARNING! WARNING! WARNING! WARNING!
    ! #################################################################
    ! #################################################################
    ! #################################################################
    !
    ! THERE IS A GIANT FLAW IN THE NEXT LOOP IF A HOST CELL ON ONE 
    ! PROCESSOR CORRESPONDS TO TWO OR MORE SEPARATE GHOST CELLS ON
    ! ANOTHER PROCESSOR
    !
    ! #################################################################
    ! #################################################################
    ! #################################################################
    !            WARNING! WARNING! WARNING! WARNING! WARNING!
    ! #################################################################
    ! #################################################################
    ! #################################################################
    !
    do nf = 1,nfbnd
      !
      if (any(bface(1,nf) == bc_comm) .and. bface(7,nf) == icpu) then
        send_msk(bface(2,nf)) = true
        recv_map(bface(6,nf)) = bface(10,nf)
      end if
      !
    end do
   !write (160+mypnum,*) icpu
   !write (170+mypnum,*) icpu
   !do n = smin,smax
   !  write (160+mypnum,*) n,send_msk(n)
   !end do
   !do n = rmin,rmax
   !  write (170+mypnum,*) n,recv_map(n)
   !end do
   !write (160+mypnum,*)
   !write (170+mypnum,*)
    !
    ! Now create the send and receive cell maps. These maps will contain
    ! in order the cells that need to sent and received from the current
    ! adjacent partition icpu
    ! NOTE: Each processor will send its cell data in order of the cell IDs.
    !       The receiving datatypes need to account for this so that the
    !       received data is placed in the correct location within an array.
    !       This is why we have created the array recv_map so that we can
    !       loop in order through the cells of the adjacent partition to
    !       create the receive offsets in the same order of the adjacent
    !       partitions cell IDs.
    !
    ! Create the send cell maps
    !
    scnt = 0
    do n = smin,smax
      if (.not. send_msk(n)) cycle
      scnt = scnt + 1
      send_temp(scnt) = n
    end do
    !
    send_cell_ptr(i+1) = send_cell_ptr(i) + scnt
    !
    s1 = send_cell_ptr(i)+1
    s2 = send_cell_ptr(i+1)
    !
    call reallocate(send_cell,s2)
    send_cell(s1:s2) = send_temp(1:scnt)
    !
    ! Create the receive cell maps
    !
    rcnt = 0
    do n = rmin,rmax
      if (recv_map(n) == 0) cycle
      rcnt = rcnt + 1
      recv_temp(rcnt) = recv_map(n)
    end do
    !
    recv_cell_ptr(i+1) = recv_cell_ptr(i) + rcnt
    !
    r1 = recv_cell_ptr(i)+1
    r2 = recv_cell_ptr(i+1)
    !
    call reallocate(recv_cell,r2)
    recv_cell(r1:r2) = recv_temp(1:rcnt)
    !
    ! Finally, reset send_msk, recv_map, send_temp,
    ! and recv_temp before next communication loop
    !
    send_msk(:) = fals
    recv_map(:) = 0
    send_temp(:) = 0
    recv_temp(:) = 0
    !
  end do
  !
  deallocate ( send_msk , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","send_msk",2,ierr)
  deallocate ( recv_map , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","recv_map",2,ierr)
  deallocate ( send_temp , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","send_temp",2,ierr)
  deallocate ( recv_temp , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","recv_temp",2,ierr)
  !
  ! Now create pointers for the send and receive solution point maps
  !
  allocate ( send_pts_ptr(1:ncomm+1) , source=0 , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","send_pts_ptr",1,ierr)
  allocate ( recv_pts_ptr(1:ncomm+1) , source=0 , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","recv_pts_ptr",1,ierr)
  !
  do i = 1,ncomm
    !
    ! 1. Cells that need to be sent
    !
    s1 = send_cell_ptr(i)+1
    s2 = send_cell_ptr(i+1)
    !
    do n = s1,s2
      nc = send_cell(n)
      np = nodes_of_cell_ptr(nc+1) - nodes_of_cell_ptr(nc)
      send_pts_ptr(i+1) = send_pts_ptr(i+1) + merge(n_sp2dt,n_sp2dq,np==3)
    end do
    send_pts_ptr(i+1) = send_pts_ptr(i+1) + send_pts_ptr(i)
    !
    ! 2. Cells that need to be received
    !
    r1 = recv_cell_ptr(i)+1
    r2 = recv_cell_ptr(i+1)
    !
    do n = r1,r2
      nc = recv_cell(n)
      np = nodes_of_cell_ptr(nc+1) - nodes_of_cell_ptr(nc)
      recv_pts_ptr(i+1) = recv_pts_ptr(i+1) + merge(n_sp2dt,n_sp2dq,np==3)
    end do
    recv_pts_ptr(i+1) = recv_pts_ptr(i+1) + recv_pts_ptr(i)
    !
  end do
  !
  ! Create the send and receive solution point maps
  !
  allocate ( send_pts(1:last(send_pts_ptr)) , source=0 , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","send_pts",1,ierr)
  allocate ( recv_pts(1:last(recv_pts_ptr)) , source=0 , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","recv_pts",1,ierr)
  !
  do i = 1,ncomm
    !
    ! 1. Cells that need to be sent
    !
    p1 = send_pts_ptr(i)+1
    p2 = send_pts_ptr(i+1)
    !
    s1 = send_cell_ptr(i)+1
    s2 = send_cell_ptr(i+1)
    !
    scnt = send_pts_ptr(i)
    do n = s1,s2
      !
      nc = send_cell(n)
      np = nodes_of_cell_ptr(nc+1) - nodes_of_cell_ptr(nc)
      np = merge(n_sp2dt,n_sp2dq,np==3)
      !
      send_pts(scnt+1:scnt+np) = cell_offset(nc) + [0:np-1]
      !
      scnt = scnt + np
      !
    end do
    !
    ! 2. Cells that need to be received
    !
    p1 = recv_pts_ptr(i)+1
    p2 = recv_pts_ptr(i+1)
    !
    r1 = recv_cell_ptr(i)+1
    r2 = recv_cell_ptr(i+1)
    !
    rcnt = recv_pts_ptr(i)
    do n = r1,r2
      !
      nc = recv_cell(n)
      np = nodes_of_cell_ptr(nc+1) - nodes_of_cell_ptr(nc)
      np = merge(n_sp2dt,n_sp2dq,np==3)
      !
      rfa = ghost_face(nc)
      rge = merge(1,2,np==n_sp2dt)
      !
      recv_pts(rcnt+1:rcnt+np) = cell_offset(nc) + geom_offset(1:np,rfa,rge)
      !
      rcnt = rcnt + np
      !
    end do
    !
  end do
  !
  ! We no longer need ghost_face so deallocate
  !
  deallocate ( ghost_face , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","ghost_face",2,ierr)
  !
  ! Finally create the MPI datatypes for exchanging
  ! solution data between adjacent partitions
  !
  allocate ( send_usp(1:ncomm) , source=mpi_datatype_null , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","send_usp",1,ierr)
  allocate ( recv_usp(1:ncomm) , source=mpi_datatype_null , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","recv_usp",1,ierr)
  !
  ! Adjust send_pts and recv_pts for the solution arrays
  !
  send_pts(:) = send_pts(:)*nq
  recv_pts(:) = recv_pts(:)*nq
  !
  do i = 1,ncomm
    !
    s1 = send_pts_ptr(i)+1
    s2 = send_pts_ptr(i+1)
    np = s2-s1+1
    !
    call mpi_type_create_indexed_block(np,nq,send_pts(s1:s2),mpi_flttyp, &
                                       send_usp(i),mpierr)
    call mpi_type_commit(send_usp(i),mpierr)
    !
    r1 = recv_pts_ptr(i)+1
    r2 = recv_pts_ptr(i+1)
    np = r2-r1+1
    !
    call mpi_type_create_indexed_block(np,nq,recv_pts(r1:r2),mpi_flttyp, &
                                       recv_usp(i),mpierr)
    call mpi_type_commit(recv_usp(i),mpierr)
    !
  end do
  !
  ! Remove the solution array adjustment from send_pts and recv_pts
  !
  send_pts(:) = send_pts(:)/nq
  recv_pts(:) = recv_pts(:)/nq
 !!
 !! Output all the solution point offsets for debugging
 !!
 !do i = 1,ncomm
 !  !
 !  icpu = exch_cpu(i)
 !  !
 !  p1 = send_cell_ptr(i)+1
 !  p2 = send_cell_ptr(i+1)
 !  !
 !  s1 = send_pts_ptr(i)+1
 !  s2 = send_pts_ptr(i+1)
 !  !
 !  write (mypnum+180,1) mypnum,icpu
 !  !
 !  do n = p1,p2
 !    nc = send_cell(n)
 !    np = nodes_of_cell_ptr(nc+1) - nodes_of_cell_ptr(nc)
 !    np = merge(n_sp2dt,n_sp2dq,np==3)
 !    s2 = s1+np-1
 !    write (mypnum+180,4) nc,send_pts(s1:s2)+1
 !    s1 = s2+1
 !  end do
 !  write (*,5) mypnum,icpu,s2,send_pts_ptr(i+1),send_pts_ptr(i+1)-send_pts_ptr(i)
 !  !
 !  p1 = recv_cell_ptr(i)+1
 !  p2 = recv_cell_ptr(i+1)
 !  !
 !  r1 = recv_pts_ptr(i)+1
 !  r2 = recv_pts_ptr(i+1)
 !  !
 !  write (mypnum+190,2) mypnum,icpu
 !  !
 !  do n = p1,p2
 !    nc = recv_cell(n)
 !    np = nodes_of_cell_ptr(nc+1) - nodes_of_cell_ptr(nc)
 !    np = merge(n_sp2dt,n_sp2dq,np==3)
 !    r2 = r1+np-1
 !    write (mypnum+190,4) nc,recv_pts(r1:r2)+1
 !    r1 = r2+1
 !  end do
 !  write (*,6) mypnum,icpu,r2,recv_pts_ptr(i+1),recv_pts_ptr(i+1)-recv_pts_ptr(i)
 !  !
 !end do
  !
  ! I don't think we need geom_offset anymore so deallocate
  !
  deallocate ( geom_offset , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","geom_offset",2,ierr)
  !
  ! Create the test MPI datatypes
  !
  allocate ( test_send(1:ncomm) , source=mpi_datatype_null , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","test_send",1,ierr)
  allocate ( test_recv(1:ncomm) , source=mpi_datatype_null , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","test_recv",1,ierr)
  !
  do i = 1,ncomm
    !
    s1 = send_pts_ptr(i)+1
    s2 = send_pts_ptr(i+1)
    np = s2-s1+1
    !
    call mpi_type_create_indexed_block(np,1,send_pts(s1:s2),mpi_flttyp, &
                                       test_send(i),mpierr)
    call mpi_type_commit(test_send(i),mpierr)
    !
    r1 = recv_pts_ptr(i)+1
    r2 = recv_pts_ptr(i+1)
    np = r2-r1+1
    !
    call mpi_type_create_indexed_block(np,1,recv_pts(r1:r2),mpi_flttyp, &
                                       test_recv(i),mpierr)
    call mpi_type_commit(test_recv(i),mpierr)
    !
  end do
  !
  ! We dont need send_pts, recv_pts, send_pts_ptr, recv_pts_ptr, send_cell,
  ! send_cell_ptr, recv_cell, or recv_cell_ptr anymore so deallocate
  !
  deallocate ( send_pts , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","send_pts",2,ierr)
  deallocate ( send_pts_ptr , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","send_pts_ptr",2,ierr)
  deallocate ( recv_pts , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","recv_pts",2,ierr)
  deallocate ( recv_pts_ptr , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","recv_pts_ptr",2,ierr)
  deallocate ( send_cell , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","send_cell",2,ierr)
  deallocate ( send_cell_ptr , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","send_cell_ptr",2,ierr)
  deallocate ( recv_cell , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","recv_cell",2,ierr)
  deallocate ( recv_cell_ptr , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","recv_cell_ptr",2,ierr)
  !
  ! These offsets are only to the grid cells of each partition. We need to
  ! account for the fact that each cell contains multiple solution points.
  ! All of the solution point data for a cell is lumped together so all
  ! that is necessary to exchange this block of information is to know the
  ! geometry type of the cell. However, the ordering of the solution points
  ! within a cell is not guaranteed to be the same across partitions. The
  ! ordering of the solution points is only dependent on the ordering of the
  ! nodes defining the element. Therefore, if we can find the ordering of the
  ! cell nodes on each partition, we can find the offsets for each solution
  ! point within the cell. The data in bface(4,:) gives the face of the local
  ! host cell that is on the boundary face and bface(9,:) gives the face of
  ! the adjacent host cell that is on the boundary face. Using these two
  ! pieces of data we should be able to find the solution point offsets.
  !
  !
  ! Deallocate before leaving
  !
  deallocate ( cell_offset , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","cell_offset",2,ierr)
  !
  ! Lets allocate the exchange buffer for the solver
  !
  allocate ( exch_buffer(1:nq,1:n_totpts) , source=0.0_r8 , stat=ierr )
  call alloc_error("create_nonblocking_exchange_datatypes","exch_buffer",1,ierr)
  !
  1 format (/," CPU ",i0," -> CPU ",i0," :: Sending")
  2 format (/," CPU ",i0," <- CPU ",i0," :: Receiving")
  3 format (*(12i6,/))
  4 format (" Cell = ",i0,*(/,12i6))
  5 format (" CPU ",i0," -> CPU ",i0," :: Sending   - ",3i6)
  6 format (" CPU ",i0," <- CPU ",i0," :: Receiving - ",3i6)
  !
end subroutine create_nonblocking_exchange_datatypes
!
!###############################################################################
!
subroutine create_face_nonblocking_datatypes(bface)
  !
  use geovar, only : n_fp1de
 !use metrics_mod, only : get_solution_points_2d
 !use geovar, only : xyz,n_global_totpts,nr,n_totpts
  !
  !.. Formal Arguments ..
  integer, dimension(:,:), intent(inout) :: bface
  !
  !.. Local Scalars ..
  integer :: n,nf,np,scnt,rcnt,i,icpu,ierr
  integer :: s1,s2,r1,r2,block_size
  integer :: rmin,rmax,adj_face,loc_face
  integer :: p1,p2,nc,nfbnd
  !
  !.. Local Arrays ..
  integer, dimension(1:ncpu) :: cpu_id
  !
  !.. Local Allocatable Arrays ..
  logical, allocatable, dimension(:) :: send_msk
  integer, allocatable, dimension(:) :: recv_map
  integer, allocatable, dimension(:) :: send_temp
  integer, allocatable, dimension(:) :: recv_temp
  !
 !real(r8), dimension(1:nr,1:n_global_totpts,1:ncpu) :: gxyz
  !
continue
  !
  !
 !call stop_recon(stop_mpi,"NEED TO SET UP FACEDUSP COMMUNICATION FIRST!!!")
  !
 !call get_solution_points_2d
 !do n = 1,ncpu
 !  if (mypnum == n-1) np = n_totpts
 !  call mpi_bcast(np,1,mpi_inttyp,n-1,mpi_comm_world,mpierr)
 !  if (mypnum == n-1) gxyz(1:nr,1:np,n) = real(xyz(1:nr,1:n_totpts),kind=r8)
 !  call mpi_bcast(gxyz(1:nr,1:np,n),np*nr,mpi_flttyp,n-1,mpi_comm_world,mpierr)
 !end do
  !
  nfbnd = size(bface,dim=2)
  !
  ! Now onto creating the communication datatypes
  !
  ! First we need to find which processes the
  ! current process needs to exchange with
  !
  cpu_id(:) = -99
  do nf = 1,nfbnd
    if ( any(bface(1,nf) == bc_comm) ) cpu_id(bface(7,nf)+1) = bface(7,nf)
  end do
  !
  ! Make sure that the current process doesnt try to communicate with itself
  !
  cpu_id(mypnum+1) = -99
  !
  ! Count the number of processes the current process needs to communicate with
  !
  ncomm = count(cpu_id >= 0)
  !
  ! Create a condensed array containing the process
  ! IDs the current process needs to communicate with
  !
  allocate ( exchf_cpu(1:ncomm) , source=pack(cpu_id,cpu_id>=0) , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","exchf_cpu",1,ierr)
  !
  ! Now loop through the processes to exchange with
  ! and create the send and receive cell maps
  !
  rmin = minval(bface(5,:))
  rmax = maxval(bface(5,:))
  allocate ( send_msk(1:nfbnd) , source=fals , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","send_msk",1,ierr)
  allocate ( recv_map(rmin:rmax) , source=0 , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","recv_map",1,ierr)
  allocate ( send_temp(1:n_fp1de*nfbnd) , source=0 , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","send_temp",1,ierr)
  allocate ( recv_temp(1:n_fp1de*nfbnd) , source=0 , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","recv_temp",1,ierr)
  allocate ( send_pts_ptr(1:ncomm+1) , source=0 , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","send_pts_ptr",1,ierr)
  allocate ( recv_pts_ptr(1:ncomm+1) , source=0 , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","recv_pts_ptr",1,ierr)
  !
  do i = 1,ncomm
    !
    icpu = exchf_cpu(i)
    !
    ! Loop to find the communication boundary faces between the
    ! current processor and the current exchange partner, icpu.
    ! We will create a mask array, send_msk, to identify all the
    ! local host cells that need to be sent to the ghost cells
    ! of icpu. We will also create an array, recv_map, that maps
    ! the host cells of icpu to the local ghost cells whose data
    ! needs to be received from icpu.
    !
    do nf = 1,nfbnd
      if (any(bface(1,nf) == bc_comm) .and. bface(7,nf) == icpu) then
        adj_face = bface(5,nf)
        send_msk(nf) = true
        recv_map(adj_face) = nf
      end if
    end do
    !
    ! Now create the send and receive cell maps. These maps will contain
    ! in order the cells that need to sent and received from the current
    ! adjacent partition icpu
    ! NOTE: Each processor will send its cell data in order of the cell IDs.
    !       The receiving datatypes need to account for this so that the
    !       received data is placed in the correct location within an array.
    !       This is why we have created the array recv_map so that we can
    !       loop in order through the cells of the adjacent partition to
    !       create the receive offsets in the same order of the adjacent
    !       partitions cell IDs.
    !
    ! Create the send cell maps
    !
    scnt = 0
    do nf = 1,nfbnd
      !
      if (.not. send_msk(nf)) cycle
      !
      ! Create the send point offsets
      !
      do n = 1,n_fp1de
        scnt = scnt + 1
        send_temp(scnt) = (nf-1)*2*n_fp1de + n-1
      end do
      !
    end do
    !
    ! Create the receive point offsets
    ! NOTE: We want to receive these in reverse order
    !       to match the host cell face points
    !
    rcnt = 0
    do nf = rmin,rmax
      !
      if (recv_map(nf) == 0) cycle
      !
      adj_face = nf
      loc_face = recv_map(nf)
      !
      do n = n_fp1de,1,-1
        rcnt = rcnt + 1
        recv_temp(rcnt) = (loc_face-1)*2*n_fp1de + n_fp1de + n-1
      end do
      !
    end do
    !
    ! Create the send_face and send_face_ptr
    ! arrays for this communication couple
    !
    send_pts_ptr(i+1) = send_pts_ptr(i) + scnt
    !
    s1 = send_pts_ptr(i)+1
    s2 = send_pts_ptr(i+1)
    !
    call reallocate(send_pts,s2)
    send_pts(s1:s2) = send_temp(1:scnt)
    !
    ! Create the recv_face and recv_face_ptr
    ! arrays for this communication couple
    !
    recv_pts_ptr(i+1) = recv_pts_ptr(i) + rcnt
    !
    r1 = recv_pts_ptr(i)+1
    r2 = recv_pts_ptr(i+1)
    !
    call reallocate(recv_pts,r2)
    recv_pts(r1:r2) = recv_temp(1:rcnt)
    !
    ! Finally, reset send_msk, recv_map, send_temp,
    ! and recv_temp before next communication loop
    !
    send_msk(:) = fals
    recv_map(:) = 0
    send_temp(:) = 0
    recv_temp(:) = 0
    !
  end do
  !
  deallocate ( send_msk , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","send_msk",2,ierr)
  deallocate ( recv_map , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","recv_map",2,ierr)
  deallocate ( send_temp , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","send_temp",2,ierr)
  deallocate ( recv_temp , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","recv_temp",2,ierr)
  !
  !##########################################################
  !##########################################################
  ! Create the MPI datatypes for exchanging
  ! solution data between adjacent partitions
  !##########################################################
  !##########################################################
  !
  allocate ( send_fusp(1:ncomm) , source=mpi_datatype_null , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","send_fusp",1,ierr)
  allocate ( recv_fusp(1:ncomm) , source=mpi_datatype_null , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","recv_fusp",1,ierr)
  !
  ! Set the block size for this MPI datatype
  !
  block_size = nq
  !
  ! Adjust send_pts and recv_pts for the solution arrays
  !
  send_pts(:) = send_pts(:)*block_size
  recv_pts(:) = recv_pts(:)*block_size
  !
  do i = 1,ncomm
    !
    s1 = send_pts_ptr(i)+1
    s2 = send_pts_ptr(i+1)
    np = s2-s1+1
    !
    call mpi_type_create_indexed_block(np,block_size,send_pts(s1:s2), &
                                       mpi_flttyp,send_fusp(i),mpierr)
    call mpi_type_commit(send_fusp(i),mpierr)
    !
    r1 = recv_pts_ptr(i)+1
    r2 = recv_pts_ptr(i+1)
    np = r2-r1+1
    !
    call mpi_type_create_indexed_block(np,block_size,recv_pts(r1:r2), &
                                       mpi_flttyp,recv_fusp(i),mpierr)
    call mpi_type_commit(recv_fusp(i),mpierr)
    !
  end do
  !
  ! Remove the solution array adjustment from send_pts and recv_pts
  !
  send_pts(:) = send_pts(:)/block_size
  recv_pts(:) = recv_pts(:)/block_size
  !
  !##########################################################
  !##########################################################
  ! Create the MPI datatypes for exchanging solution
  ! gradients data between adjacent partitions
  !##########################################################
  !##########################################################
  !
  allocate ( send_fdusp(1:ncomm) , source=mpi_datatype_null , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","send_fdusp",1,ierr)
  allocate ( recv_fdusp(1:ncomm) , source=mpi_datatype_null , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","recv_fdusp",1,ierr)
  !
  ! Set the block size for this MPI datatype
  !
  block_size = nr*nq
  !
  ! Adjust send_pts and recv_pts for the solution arrays
  !
  send_pts(:) = send_pts(:)*block_size
  recv_pts(:) = recv_pts(:)*block_size
  !
  do i = 1,ncomm
    !
    s1 = send_pts_ptr(i)+1
    s2 = send_pts_ptr(i+1)
    np = s2-s1+1
    !
    call mpi_type_create_indexed_block(np,block_size,send_pts(s1:s2), &
                                       mpi_flttyp,send_fdusp(i),mpierr)
    call mpi_type_commit(send_fdusp(i),mpierr)
    !
    r1 = recv_pts_ptr(i)+1
    r2 = recv_pts_ptr(i+1)
    np = r2-r1+1
    !
    call mpi_type_create_indexed_block(np,block_size,recv_pts(r1:r2), &
                                       mpi_flttyp,recv_fdusp(i),mpierr)
    call mpi_type_commit(recv_fdusp(i),mpierr)
    !
  end do
  !
  ! Remove the solution array adjustment from send_pts and recv_pts
  !
  send_pts(:) = send_pts(:)/block_size
  recv_pts(:) = recv_pts(:)/block_size
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
 !  icpu = exchf_cpu(i)
 !  !
 !  s1 = send_pts_ptr(i)+1
 !  s2 = send_pts_ptr(i+1)
 !  !
 !  r1 = recv_pts_ptr(i)+1
 !  r2 = recv_pts_ptr(i+1)
 !  !
 !  write (mypnum+180,1) mypnum,icpu
 !  !
 !  do n = s1,s2
 !    p1 = send_pts(n)+1
 !    p2 = recv_pts(n)+1
 !    write (mypnum+180,7) n,p1,gxyz(1,p1,mypnum+1),gxyz(2,p1,mypnum+1), &
 !                           p2,gxyz(1,p2,icpu+1),gxyz(2,p2,icpu+1)
 !  end do
 !  !
 !  write (*,5) mypnum,icpu,s2,send_pts_ptr(i+1),send_pts_ptr(i+1)-send_pts_ptr(i)
 !  write (*,6) mypnum,icpu,r2,recv_pts_ptr(i+1),recv_pts_ptr(i+1)-recv_pts_ptr(i)
 !  !
 !end do
  !
  block_size = nr
  !
  send_pts(:) = send_pts(:) * block_size
  recv_pts(:) = recv_pts(:) * block_size
  !
  ! Create the test MPI datatypes
  !
  allocate ( testf_send(1:ncomm) , source=mpi_datatype_null , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","testf_send",1,ierr)
  allocate ( testf_recv(1:ncomm) , source=mpi_datatype_null , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","testf_recv",1,ierr)
  !
  do i = 1,ncomm
    !
    s1 = send_pts_ptr(i)+1
    s2 = send_pts_ptr(i+1)
    np = s2-s1+1
    !
    call mpi_type_create_indexed_block(np,block_size,send_pts(s1:s2),mpi_flttyp, &
                                       testf_send(i),mpierr)
    call mpi_type_commit(testf_send(i),mpierr)
    !
    r1 = recv_pts_ptr(i)+1
    r2 = recv_pts_ptr(i+1)
    np = r2-r1+1
    !
    call mpi_type_create_indexed_block(np,block_size,recv_pts(r1:r2),mpi_flttyp, &
                                       testf_recv(i),mpierr)
    call mpi_type_commit(testf_recv(i),mpierr)
    !
  end do
  !
  ! We dont need send_pts, recv_pts, send_pts_ptr,
  ! or recv_pts_ptr anymore so deallocate
  !
  deallocate ( send_pts , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","send_pts",2,ierr)
  deallocate ( send_pts_ptr , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","send_pts_ptr",2,ierr)
  deallocate ( recv_pts , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","recv_pts",2,ierr)
  deallocate ( recv_pts_ptr , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","recv_pts_ptr",2,ierr)
  !
  ! These offsets are only to the grid cells of each partition. We need to
  ! account for the fact that each cell contains multiple solution points.
  ! All of the solution point data for a cell is lumped together so all
  ! that is necessary to exchange this block of information is to know the
  ! geometry type of the cell. However, the ordering of the solution points
  ! within a cell is not guaranteed to be the same across partitions. The
  ! ordering of the solution points is only dependent on the ordering of the
  ! nodes defining the element. Therefore, if we can find the ordering of the
  ! cell nodes on each partition, we can find the offsets for each solution
  ! point within the cell. The data in bface(4,:) gives the face of the local
  ! host cell that is on the boundary face and bface(9,:) gives the face of
  ! the adjacent host cell that is on the boundary face. Using these two
  ! pieces of data we should be able to find the solution point offsets.
  !
  ! Lets allocate the exchange buffer for the solver
  !
  allocate ( fusp_buffer(1:nq,1:n_fp1de,1:2,1:nfbnd) , source=0.0_r8 , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","fusp_buffer",1,ierr)
  allocate ( fusp_rqst(1:2*ncomm) , source=0 , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","fusp_rqst",1,ierr)
  allocate ( fusp_istat(1:mpi_status_size,1:2*ncomm) , source=0 , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","fusp_istat",1,ierr)
  !
  allocate ( fdusp_buffer(1:nr,1:nq,1:n_fp1de,1:2,1:nfbnd) , source=0.0_r8 , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","fdusp_buffer",1,ierr)
  allocate ( fdusp_rqst(1:2*ncomm) , source=0 , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","fdusp_rqst",1,ierr)
  allocate ( fdusp_istat(1:mpi_status_size,1:2*ncomm) , source=0 , stat=ierr )
  call alloc_error("create_face_nonblocking_datatypes","fdusp_istat",1,ierr)
  !
  1 format (/," CPU ",i0," -> CPU ",i0," :: Sending")
  2 format (/," CPU ",i0," <- CPU ",i0," :: Receiving")
  3 format (*(12i6,/))
  4 format (" Cell = ",i0,*(/,12i6))
  5 format (" CPU ",i0," -> CPU ",i0," :: Sending   - ",3i6)
  6 format (" CPU ",i0," <- CPU ",i0," :: Receiving - ",3i6)
  7 format (2x,i6,": Sending - ",i0,4x,2f12.4,/, &
            "        Receiving - ",i0,4x,2f12.4)  
  !
end subroutine create_face_nonblocking_datatypes
!
!###############################################################################
!
subroutine create_geometry_offsets()
  !
  !.. Local Scalars ..
  integer :: i,j,n,lp,gp,ierr
  !
continue
  !
  allocate ( geom_offset(1:n_sp2dq,1:4,1:2) , source=0 , stat=ierr )
  call alloc_error("create_geometry_offsets","geom_offset",1,ierr)
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Create the offsets for triangle cells
  !
  ! 1a. Create the triangle offsets for bface(9,:)=1
  !
  lp = 0
  do j = 1,n_sp1de
    do i = 1,n_sp1de-j+1
      lp = lp + 1
      gp = lp
      geom_offset(lp,1,1) = gp-1
    end do
  end do
  !
  ! 1b. Create the triangle offsets for bface(9,:)=2
  !
  lp = 0
  do j = 1,n_sp1de
    lp = lp + 1
    gp = n_sp2dt - sum([0:j-1])
    geom_offset(lp,2,1) = gp-1
    do i = 1,n_sp1de-j
      lp = lp + 1
      gp = gp - i - j
      geom_offset(lp,2,1) = gp-1
    end do
  end do
  !
  ! 1c. Create the triangle offsets for bface(9,:)=3
  !
  lp = 0
  do j = 1,n_sp1de
    gp = 1 - j
    do i = 1,n_sp1de-j+1
      lp = lp + 1
      gp = gp + (n_sp1de-i+1)
      geom_offset(lp,3,1) = gp-1
    end do
  end do
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Create the offsets for quad cells
  !
  ! 2a. Create the quad offsets for bface(9,:)=1
  !
  lp = 0
  do j = 1,n_sp1de
    do i = 1,n_sp1de
      lp = lp + 1
      gp = (j-1)*n_sp1de + i
      geom_offset(lp,1,2) = gp-1
    end do
  end do
  !
  ! 2b. Create the quad offsets for bface(9,:)=2
  !
  lp = 0
  do j = 1,n_sp1de
    do i = 1,n_sp1de
      lp = lp + 1
      gp = j + n_sp1de*(n_sp1de-i)
      geom_offset(lp,2,2) = gp-1
    end do
  end do
  !
  ! 2c. Create the quad offsets for bface(9,:)=3
  !
  lp = 0
  do j = 1,n_sp1de
    do i = 1,n_sp1de
      lp = lp + 1
      gp = (n_sp1de-j)*n_sp1de + (n_sp1de-i+1)
      geom_offset(lp,3,2) = gp-1
    end do
  end do
  !
  ! 2d. Create the quad offsets for bface(9,:)=4
  !
  lp = 0
  do j = 1,n_sp1de
    do i = 1,n_sp1de
      lp = lp + 1
      gp = (n_sp1de-j+1) + n_sp1de*(i-1)
      geom_offset(lp,4,2) = gp-1
    end do
  end do
 !!
 !! Output geometry offsets for debugging
 !!
 !if (mypnum == 0) then
 !  !
 !  write (99,1) "Triangle"
 !  !
 !  do i = 1,3
 !    !
 !    write (99,2) i
 !    lp = n_sp2dt
 !    do j = n_sp1de,1,-1
 !      write (99,3) (1+geom_offset(lp+gp,i,1),gp=0,n_sp1de-j)
 !      lp = lp - (n_sp1de-j+2)
 !    end do
 !    !
 !  end do
 !  !
 !  write (99,1) "Quad"
 !  !
 !  do i = 1,4
 !    !
 !    write (99,2) i
 !    do j = n_sp1de,1,-1
 !      write (99,3) (1+geom_offset(n_sp1de*(j-1)+gp,i,2),gp=1,n_sp1de)
 !    end do
 !    !
 !  end do
 !  !
 !end if
 !!
 !1 format (/," ",a," offsets",/)
 !2 format ("Face ",i0," of cell on boundary")
 !3 format (<n_sp1de>i5)
  !
end subroutine create_geometry_offsets
!
!###############################################################################
!
!subroutine create_nonblocking_exchange_datatypes()
!  !
!  !.. Formal Arguments ..
!  !
!  !.. Local Scalars ..
!  !
!  !.. Local Allocatable Arrays ..
!  logical, allocatable, dimension(:) :: send_msk
!  logical, allocatable, dimension(:) :: recv_msk
!  integer, allocatable, dimension(:) :: send_offsets
!  integer, allocatable, dimension(:) :: recv_offsets
!  !
!continue
!  !
!  ! Allocate the send and receive mask arrays
!  !
!  mxcell = maxval( comm_cpu(2,:) )
!  !
!  rmin = minval( comm_links(2,1,:) )
!  rmax = maxval( comm_links(2,1,:) )
!  smin = minval( comm_links(2,2,:) )
!  smax = maxval( comm_links(2,2,:) )
!  !
!  allocate ( send_msk(smin:smax) , source=fals , stat=ierr )
!  allocate ( recv_msk(rmin:rmax) , source=fals , stat=ierr )
!  !
!  ! Loop through the communication partners
!  !
!  iend = 0
!  send_count = 0
!  recv_count = 0
!  !
!  do npb = 1,ncomm
!    !
!    rcpu = comm_cpu(1,npb)
!    !
!    ibeg = iend + 1
!    iend = iend + comm_cpu(2,npb)
!    !
!    ! Mark the send and receive data locations
!    !
!    do i = ibeg,iend
!      !
!      if (comm_links(3,2,i) == rcpu) then
!        send_msk(comm_links(2,2,i)) = true
!        recv_msk(comm_links(2,1,i)) = true
!      end if
!      !
!    end do
!    !
!    ! Create the sending offsets
!    !
!    do i = smin,smax
!      if (.not. send_msk(i)) cycle
!      send_count = send_count + 1
!      send_offsets(send_count) = i - 1
!    end do
!    !
!    ! Create the receiving offsets
!    !
!    do i = rmin,rmax
!      if (.not. recv_msk(i)) cycle
!      recv_count = recv_count + 1
!      recv_offsets(recv_count) = i - 1
!    end do
!    !
!  end do
!  !
!end subroutine create_nonblocking_exchange_datatypes
!!
!!###############################################################################
!!
!subroutine create_blocking_exchange_datatypes()
!  !
!  !.. Local Scalars ..
!  integer :: nbp,lcpu,rcpu,block,nsnd,nrcv
!  !
!  !.. Local Allocatable Arrays ..
!  integer, allocatable, dimension(:) :: send_offsets
!  integer, allocatable, dimension(:) :: recv_offsets
!  !
!continue
!  !
!  ! Block size assuming <nr> momentum equations
!  ! plus one density plus one energy equation
!  !
!  block = nq
!  ! 
!  ! Allocate the sending and receiving datatypes
!  !
!  allocate ( snd_usp(1:n_com_bnd) , source=mpi_datatype_null , stat=ierr )
!  allocate ( rcv_usp(1:n_com_bnd) , source=mpi_datatype_null , stat=ierr )
!  !
!  ! Allocate and initialize the offset arrays
!  !
!  offmax = n_sp2dq*maxval(comm_links_ptr(3,:))
!  !
!  allocate ( send_offsets(1:offmax) , source=0 , stat=ierr )
!  allocate ( recv_offsets(1:offmax) , source=0 , stat=ierr )
!  !
!  ! Create the MPI datatypes to communicate data across partition boundaries
!  !
!  do npb = 1,n_com_bnd
!    !
!    lcpu = comm_links_ptr(1,npb)
!    rcpu = comm_links_ptr(2,npb)
!    !
!    ! For this partition boundary, find the number of solution points
!    ! that require the exchange of information on each processor.
!    !
!    if (mypnum == lcpu) then
!      !
!      ! Processor for the left partition on the current boundary
!      !
!    else if (mypnum == rcpu) then
!      !
!      ! Processor for the right partition on the current boundary
!      !
!    end if
!    !
!    ! Create the send datatype for the solution array usp
!    !
!    call mpi_type_create_indexed_block(nsnd,block,send_offsets(1:nsnd), &
!                                       mpi_fltval,snd_usp(npb),mpierr)
!    call mpi_type_commit(snd_usp(npb),mpierr)
!    !
!    ! Create the receive datatype for the solution array usp
!    !
!    call mpi_type_create_indexed_block(nrcv,block,recv_offsets(1:nrcv), &
!                                       mpi_fltval,rcv_usp(npb),mpierr)
!    call mpi_type_commit(rcv_usp(npb),mpierr)
!    !
!  end do
!  !
!end subroutine create_blocking_exchange_datatypes
!
!###############################################################################
!
!subroutine create_communication_links2()
!  !
!  ! Create only those communication links relevant to each processor
!  ! separately and combine into a global array once completed.
!  !
!  !.. Formal Arguments ..
!  integer,                 intent(in) :: npart
!  integer, dimension(:),   intent(in) :: epart
!  integer, dimension(:,:), intent(in) :: bface
!  !
!  !.. Local Scalars ..
!  integer :: nf
!  !
!  !.. Local Allocatable Arrays ..
!  integer, allocatable, dimension(:) :: link_count
!  !
!continue
!  !
!  ! Loop through the boundary faces and find the number of
!  ! processors the current processor has to communicate with.
!  !
!  allocate ( link_count(1:npart) , source=0 , stat=ierr )
!  !
!  do nf = 1,ngfbnd
!    !
!    ! Cycle if not a communication boundary face
!    !
!    if (.not. any(bface(1,nf) == bc_comm)) cycle
!    !
!    ! If applicable to the current processor, add one to the
!    ! number of links for the adjacent processor.
!    ! NOTE: Make sure that the communication is to an external
!    !       processor and not for two cells on the current
!    !       processor linked by a periodic boundary condition.
!    !
!    lcpu = epart(bface(2,nf))-1
!    rcpu = bface(7,nf)
!    !
!    if (lcpu == mypnum .and. rcpu /= mypnum) &
!      link_count(rcpu+1) = link_count(rcpu+1) + 1
!    !
!  end do
!  !
!  ! Count the number of processors to communicate with and the
!  ! total number of communication links for all processors
!  !
!  ncomm = count( link_count(:) /= 0 )
!  nlink = sum( link_count(:) )
!  !
!  ! Create the list of processors to communicate with
!  !
!  allocate ( comm_cpu(1:2,1:ncomm) , source=0 , stat=ierr )
!  !
!  i = 0
!  do n = 1,npart
!    if (link_count(n) /= 0) then
!      i = i + 1
!      comm_cpu(1,i) = n - 1
!      comm_cpu(2,i) = link_count(n)
!    end if
!  end do
!  !
!  ! Now create the list of communication links for each processor
!  !
!  allocate ( comm_links(1:4,1:2,1:nlink) , source=0 , stat=ierr )
!  !
!  i = 0
!  !
!  do npb = 1,ncomm
!    !
!    rcpu = comm_cpu(1,npb)
!    !
!    do nf = 1,ngfbnd
!      !
!      ! Cycle if not a communication boundary face
!      !
!      if (.not. any(bface(1,nf) == bc_comm)) cycle
!      !
!      lcpu = epart(bface(2,nf))-1
!      !
!      if (lcpu == mypnum .and. rcpu == bface(7,nf)) then
!        !
!        ! Create the link information for this communication boundary face
!        !
!        ! for communication boudnary face i_lnk:
!        !
!        ! comm_links(1,1,i_lnk) = for lcpu, global boundary face
!        ! comm_links(2,1,i_lnk) = for lcpu, global host  cell
!        ! comm_links(3,1,i_lnk) = for lcpu, processor ID
!        ! comm_links(4,1,i_lnk) = for lcpu, number of nodes defining host cell
!        !
!        ! comm_links(1,2,i_lnk) = for rcpu, global boundary face
!        ! comm_links(2,2,i_lnk) = for rcpu, global host  cell
!        ! comm_links(3,2,i_lnk) = for rcpu, processor ID
!        ! comm_links(4,2,i_lnk) = for rcpu, number of nodes defining host cell
!        !
!        i = i + 1
!        !
!        ! Link information for left cell
!        !
!        comm_links(1,1,i) = nf
!        comm_links(2,1,i) = bface(2,nf)
!        comm_links(3,1,i) = lcpu
!        comm_links(4,1,i) = nodes_of_cell_ptr(bface(2,nf)+1) - &
!                            nodes_of_cell_ptr(bface(2,nf))
!        !
!        ! Link information for right cell
!        !
!        comm_links(1,2,i) = bface(5,nf)
!        comm_links(2,2,i) = bface(6,nf)
!        comm_links(3,2,i) = rcpu
!        comm_links(4,2,i) = nodes_of_cell_ptr(bface(6,nf)+1) - &
!                            nodes_of_cell_ptr(bface(6,nf))
!        !
!      end if
!      !
!    end do
!    !
!  end do
!  !
!  ! 
!
!
!
!end subroutine create_communication_links2
!!
!!###############################################################################
!!
!subroutine create_communication_links(npart,bface)
!  !
!  !.. Formal Arguments ..
!  integer, intent(in) :: npart
!  integer, dimension(:,:), intent(in) :: bface
!  !
!  !.. Local Scalars ..
!  integer :: llnk,nlnk,mlnk,i,ierr,np1,np2,n,nf,lcpu,rcpu
!  !
!  !.. Local Arrays ..
!  integer, dimension(1:npart,1:npart) :: link_count
!  !
!continue
!  !
!  ! Mark the number of links for each processor
!  !
!  link_count(:,:) = 0
!  do n = 1,npart
!    do nf = 1,ngfbnd
!      if ( any(bface(1,nf) == bc_comm) ) then
!        link_count(n,bface(7,nf)+1) = link_count(n,bface(7,nf)+1) + 1
!      end if
!    end do
!  end do
!  !
!  ! Count the total number of communication boundaries excluding duplicates
!  !
!  n_com_bnd = 0
!  do np1 = npart,1,-1
!    do np2 = 1,np1-1
!      if (link_count(np1,np2) /= 0) n_com_bnd = n_com_bnd + 1
!    end do
!  end do
!  !
!  ! Create the pointer for communication links
!  !
!  allocate ( comm_links_ptr(1:3,1:n_com_bnd) , source=0 , stat=ierr )
!  !
!  llnk = 0
!  do np1 = npart,1,-1
!    do np2 = 1,np1-1
!      if (link_count(np1,np2) /= 0) then
!        llnk = llnk + 1
!        comm_links_ptr(1,llnk) = np1 - 1
!        comm_links_ptr(2,llnk) = np2 - 1
!        comm_links_ptr(3,llnk) = link_count(np1,np2)
!      end if
!    end do
!  end do
!  !
!  ! Count the total number of communication links on each local processor
!  !
! !n_com_lnk = 0
! !do n = 1,n_com_bnd
! !  if ( any(mypnum == comm_links_ptr(1:2,n)) ) then
! !    n_com_lnk = n_com_lnk + comm_links_ptr(3,n)
! !  end if
! !end do
!  n_com_lnk = sum( comm_links_ptr(3,:) , &
!                   mask = any(mypnum==comm_links_ptr(1:2,:),dim=1) )
!  !
!  ! Create the communication links on each local processor
!  !
!  allocate ( comm_links(1:3,1:2,1:n_com_lnk) , source=0 , stat=ierr )
!  !
!  nlnk = 0  ! count of total communication links on local processor
!  mlnk = 0  ! count of total communication links on all   processors
!  !
!  !########################################################################
!  !########################################################################
!  !########################################################################
!  !########################################################################
!  !########################################################################
!  !########################################################################
!  !
!  ! THIS NEEDS TO BE REWORKED BECAUSE I THINK THERE WILL BE DUPLICATE
!  ! COMMUNICATION LINKS SINCE PARTITION BOUNDARY FACES WERE SPLIT FROM
!  ! ONE INTERIOR FACE TO TWO BOUNDARY FACES WITH ONE DEFINED FOR EACH
!  ! PARTITION. THIS NEXT LOOP NEEDS TO BE THOUGHT OUT A LITTLE MORE.
!  !
!  !########################################################################
!  !########################################################################
!  !########################################################################
!  !########################################################################
!  !########################################################################
!  !########################################################################
!  do n = 1,n_com_bnd
!    !
!    lcpu = comm_links_ptr(1,n)
!    rcpu = comm_links_ptr(2,n)
!    !
!    llnk = 0  ! count of total communication links for current boundary
!    !
!    if ( any(mypnum == [lcpu,rcpu]) ) then
!      !
!      do nf = 1,ngfbnd
!        !
!        if ( any(bface(1,nf) == bc_comm) .and. &
!             any(bface(7,nf) == [lcpu,rcpu]) ) then
!          llnk = llnk + 1
!          nlnk = nlnk + 1
!          if (bface(7,nf) == lcpu) then
!            comm_links(1,1,nlnk) = bface(5,nf)
!            comm_links(2,1,nlnk) = bface(6,nf)
!            comm_links(3,1,nlnk) = 
!            comm_links(1,2,nlnk) = nf
!            comm_links(2,2,nlnk) = bface(2,nf)
!            comm_links(3,2,nlnk) = 
!          else
!            comm_links(1,2,nlnk) = nf
!            comm_links(2,2,nlnk) = bface(2,nf)
!            comm_links(3,1,nlnk) = 
!            comm_links(1,2,nlnk) = bface(5,nf)
!            comm_links(2,2,nlnk) = bface(6,nf)
!            comm_links(3,2,nlnk) = 
!          end if
!        end if
!        !
!      end do
!      !
!      if (llnk /= comm_links_ptr(3,n)) &
!        write (iout,1) mypnum,n,llnk,comm_links_ptr(3,n)
!      !
!    end if
!    !
!    call mpi_bcast(llnk,1,mpi_inttyp,lcpu,mpi_comm_world,mpierr)
!    !
!    mlnk = mlnk + llnk
!    !
!  end do
!  !
!  ! Check that the number of communication
!  ! links on each processor is as expected
!  !
!  i = 0
!  do n = 0,npart-1
!    if (mypnum == n) then
!      if (nlnk /= n_com_lnk) then
!        write (iout,5) mypnum,nlnk,n_com_lnk
!        i = 1
!      end if
!    end if
!    call mpi_barrier(mpi_comm_world,mpierr)
!  end do
!  !
!  ! Check that the global number of links is as expected
!  !
!  if (mlnk /= sum(comm_links_ptr(3,:))) then
!    if (mypnum == 0) then
!      write (iout,6) mlnk,sum(comm_links_ptr(3,:))
!    end if
!    i = 1
!  end if
!  !
!  ! Sync the error flag across all processors and kill if error detected
!  !
!  ierr = i
!  call mpi_allreduce(i,ierr,1,mpi_inttyp,mpi_max,mpi_comm_world,mpierr)
!  if (ierr /= 0) call stop_recon(stop_mpi,"create_communication_links")
!  !
!  ! Format Statements
!  !
!  4  format (/,"ERROR: On processor ",i0,",",/, &
!               "The number of links found for communication",/, &
!               "boundary ",i0," does not match what was expected!",/, &
!               "Number of links found   : ",i0,/, &
!               "Number of links expected: ",i0,/)
!  5  format (/,"ERROR: On processor ",i0,",",/, &
!               "The total number of communication links found on the",/, &
!               "current processor does not match what was expected!",/, &
!               "Total number of links found    on processor: ",i0,/, &
!               "Total number of links expected on processor: ",i0,/)
!  6  format (/,"ERROR COUNTING COMMUNICATION LINKS!",/, &
!               "GLOBAL NUMBER OF LINKS FOUND   : ",i0,/, &
!               "GLOBAL NUMBER OF LINKS EXPECTED: ",i0,/)
!  !
!end subroutine create_communication_links
!
!###############################################################################
!
subroutine blocking_exchange(usp)
  !
  !.. Formal Arguments ..
  real(vbp), dimension(:,:), intent(inout) :: usp
  !
  !.. Local Scalars ..
  integer :: i,n,lcpu,rcpu,ltag,rtag
  !
  !.. Local Arrays ..
  integer, dimension(1:mpi_status_size) :: istat
  !
continue
  !
  ! Loop over the communication boundaries between partitions
  !
  i = 0
  !
  do n = 1,n_com_bnd
    !
    lcpu = comm_links_ptr(1,n)
    rcpu = comm_links_ptr(2,n)
    !
    ! Put a barrier here to ensure against deadlock
    ! NOTE: Not sure if this is needed
    !
    call mpi_barrier(mpi_comm_world,mpierr)
    !
    ! Have the processors associated with these two partitions
    ! on the current partition boundary exchange variables
    !
    ltag = n + lcpu
    rtag = n + rcpu
    !
    if (mypnum == lcpu) then
      !
      i = i + 1
      !
      call mpi_send(usp,1,snd_usp(i),rcpu,ltag,mpi_comm_world,mpierr)
      !
      call mpi_recv(usp,1,rcv_usp(i),rcpu,rtag,mpi_comm_world,istat,mpierr)
      !
    else if (mypnum == rcpu) then
      !
      i = i + 1
      !
      call mpi_recv(usp,1,rcv_usp(i),lcpu,ltag,mpi_comm_world,istat,mpierr)
      !
      call mpi_send(usp,1,snd_usp(i),lcpu,rtag,mpi_comm_world,mpierr)
      !
    end if
    !
  end do
  !
end subroutine blocking_exchange
!
!###############################################################################
!
subroutine exchange_faceusp(faceusp)
  !
  !.. Use Statements ..
  use geovar, only : nfbnd
  !
  !.. Formal Arguments ..
  real(vbp), dimension(:,:,:,:), intent(in) :: faceusp
  !
  !.. Local Scalars..
  integer :: n,s,r,dtype
  !
  !.. Local Arrays ..
  real(r8), dimension(1:size(faceusp,1),1:size(faceusp,2), &
                      1:size(faceusp,3),1:nfbnd) :: exch_send
  !
continue
  !
  fusp_buffer(:,:,:,:) = real(faceusp(:,:,:,1:nfbnd),kind=r8)
  exch_send(:,:,:,:)   = real(faceusp(:,:,:,1:nfbnd),kind=r8)
  !
  ! Initialize all the receives so each process
  ! is ready once the send is posted
  !
  do n = 1,ncomm
    s = exchf_cpu(n)
    dtype = recv_fusp(n)
    call mpi_irecv(fusp_buffer,1,dtype,s,s,mpi_comm_world,fusp_rqst(n),mpierr)
  end do
  !
  ! Send out the adjacent cells
  !
  do n = 1,ncomm
    r = exchf_cpu(n)
    dtype = send_fusp(n)
    call mpi_isend(exch_send,1,dtype,r,mypnum, &
                   mpi_comm_world,fusp_rqst(n+ncomm),mpierr)
  end do
  !
end subroutine exchange_faceusp
!
!###############################################################################
!
subroutine wait_faceusp(faceusp)
  !
  !.. Use Statements ..
  use geovar, only : nfbnd
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:,:,:), intent(inout) :: faceusp
  !
continue
  !
  ! Wait for all faceusp sends and receives to complete
  !
  call mpi_waitall(2*ncomm,fusp_rqst,fusp_istat,mpierr)
  !
  ! Update the right states for all the boundary faces
  !
  faceusp(:,:,2,1:nfbnd) = real(fusp_buffer(:,:,2,1:nfbnd),kind=wp)
  !
end subroutine wait_faceusp
!
!###############################################################################
!
subroutine exchange_facedusp(facedusp)
  !
  !.. Use Statements ..
  use geovar, only : nfbnd
  !
  !.. Formal Arguments ..
  real(vbp), dimension(:,:,:,:,:), intent(in) :: facedusp
  !
  !.. Local Scalars..
  integer :: n,s,r,dtype
  !
  !.. Local Arrays ..
  real(r8), dimension(1:size(facedusp,1),1:size(facedusp,2), &
                      1:size(facedusp,3),1:size(facedusp,4), &
                      1:nfbnd) :: exch_send
  !
continue
  !
  fdusp_buffer(:,:,:,:,:) = real(facedusp(:,:,:,:,1:nfbnd),kind=r8)
  exch_send(:,:,:,:,:)    = real(facedusp(:,:,:,:,1:nfbnd),kind=r8)
  !
  ! Initialize all the receives so each process
  ! is ready once the send is posted
  !
  do n = 1,ncomm
    s = exchf_cpu(n)
    dtype = recv_fdusp(n)
    call mpi_irecv(fdusp_buffer,1,dtype,s,s,mpi_comm_world,fdusp_rqst(n),mpierr)
  end do
  !
  ! Send out the adjacent cells
  !
  do n = 1,ncomm
    r = exchf_cpu(n)
    dtype = send_fdusp(n)
    call mpi_isend(exch_send,1,dtype,r,mypnum, &
                   mpi_comm_world,fdusp_rqst(n+ncomm),mpierr)
  end do
  !
end subroutine exchange_facedusp
!
!###############################################################################
!
subroutine wait_facedusp(facedusp)
  !
  !.. Use Statements ..
  use geovar, only : nfbnd
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:,:,:,:), intent(inout) :: facedusp
  !
continue
  !
  ! Wait for all facedusp sends and receives to complete
  !
  call mpi_waitall(2*ncomm,fdusp_rqst,fdusp_istat,mpierr)
  !
  ! Update the right states for all the boundary faces
  !
  facedusp(:,:,:,2,1:nfbnd) = real(fdusp_buffer(:,:,:,2,1:nfbnd),kind=wp)
  !
end subroutine wait_facedusp
!
!###############################################################################
!
subroutine nonblocking_exchange(usp)
  !
  use geovar, only : n_totpts,n_solpts
  !
  !.. Formal Arguments ..
  real(vbp), dimension(:,:), intent(inout) :: usp
  !
  !.. Local Scalars..
  integer :: n,s,r,dtype
  !
  !.. Local Arrays ..
  integer, dimension(1:2*ncomm) :: rqst
  integer, dimension(1:mpi_status_size,1:2*ncomm) :: istat
  real(r8), dimension(1:size(usp,1),1:size(usp,2)) :: exch_send
  !
continue
  !
  exch_buffer(:,:) = real(usp(:,:),r8)
  exch_send(:,:) = real(usp(:,:),r8)
  !
  ! Initialize all the receives so each process
  ! is ready once the send is posted
  !
  do n = 1,ncomm
    s = exch_cpu(n)
    dtype = recv_usp(n)
    call mpi_irecv(exch_buffer,1,dtype,s,s,mpi_comm_world,rqst(n),mpierr)
  end do
  !
  ! Send out the adjacent cells
  !
  do n = 1,ncomm
    r = exch_cpu(n)
    dtype = send_usp(n)
    call mpi_isend(exch_send,1,dtype,r,mypnum,mpi_comm_world,rqst(n+ncomm),mpierr)
  end do
  !
  ! Wait for all sends and receives to complete
  !
  call mpi_waitall(2*ncomm,rqst,istat,mpierr)
  !
  ! Update the adjacent elements
  !
  usp(:,n_solpts+1:n_totpts) = real(exch_buffer(:,n_solpts+1:n_totpts),vbp)
  !
end subroutine nonblocking_exchange
!
!###############################################################################
!
subroutine test_nonblocking_fexchange()
  !
  use geovar, only : nr,xyz,nfbnd,cells_on_face,pts_in_cell_ptr,interp_s2f,bface
  !
  !.. Local Scalars..
  integer :: i,n,s,r,dtype,ierr,nf,jl,kbl,l1,l2,nl,wl,k,kl
  real(wp) :: tmp
  !
  !.. Local Arrays ..
  integer, dimension(1:2*ncomm) :: rqst
  integer, dimension(1:mpi_status_size,1:2*ncomm) :: istat
  !
  !.. Local Allocatable Arrays ..
  real(r8), dimension(1:nr,1:n_sp1de,1:2,1:nfbnd) :: recv_buffer
  real(r8), dimension(1:nr,1:n_sp1de,1:2,1:nfbnd) :: send_buffer
  !
continue
  !
  send_buffer(:,:,:,:) = real(0.0,kind=r8)
  recv_buffer(:,:,:,:) = real(0.0,kind=r8)
  do nf = 1,nfbnd
    !
    jl = cells_on_face(1,nf)
    kbl = cells_on_face(5,nf)
    !
    l1 = pts_in_cell_ptr(jl)+1
    l2 = pts_in_cell_ptr(jl+1)
    !
    nl = l2-l1+1
    wl = merge(1,2,nl==n_sp2dt)
    !
    do k = 1,n_sp1de
      kl = kbl-1 + k
      do n = 1,nr
        tmp = dot_product(interp_s2f(1:nl,kl,wl),xyz(n,l1:l2))
        send_buffer(n,k,1,nf) = real(tmp,kind=r8)
      end do
    end do
    !
  end do
  !
  recv_buffer(:,:,:,:) = send_buffer(:,:,:,:)
  !
  ! Initialize all the receives so each process
  ! is ready once the send is posted
  !
  do n = 1,ncomm
    s = exchf_cpu(n)
    dtype = testf_recv(n)
    call mpi_irecv(recv_buffer,1,dtype,s,s,mpi_comm_world,rqst(n),mpierr)
  end do
  !
  ! Send out the adjacent cells
  !
  do n = 1,ncomm
    r = exchf_cpu(n)
    dtype = testf_send(n)
    call mpi_isend(send_buffer,1,dtype,r,mypnum,mpi_comm_world,rqst(n+ncomm),mpierr)
  end do
  !
  ! Wait for all sends and receives to complete
  !
  call mpi_waitall(2*ncomm,rqst,istat,mpierr)
  !
  ! Update the adjacent elements
  !
  send_buffer(:,:,:,:) = recv_buffer(:,:,:,:)
  !
  do nf = 1,nfbnd
    if (bface(1,nf) /= bc_cpu_bnd) cycle
    write (150+mypnum,1) nf,cells_on_face(1:2,nf)
    do k = 1,n_sp1de
      write (150+mypnum,2) k,send_buffer(1,k,1,nf),send_buffer(2,k,1,nf), &
                             send_buffer(1,k,2,nf),send_buffer(2,k,2,nf)
    end do
  end do
  !
 !call stop_recon(stop_mpi,"test_nonblocking_fexchange")
  !
  1 format (" Face = ",i6,4x,"=> Left cell = ",i5,"  |  Right Cell = ",i5)
  2 format ("    k = ",i2," -   xl,yl = ",g12.4,",",g12.4, &
                          "  |  xr,yr = ",g12.4,",",g12.4)
 !2 format ("    k = ",i2," -   xl = ",g12.4,"  |  xr = ",g12.4,/, &
 !          "               yl = ",g12.4,"  |  yr = ",g12.4)
  !
end subroutine test_nonblocking_fexchange
!
!###############################################################################
!
subroutine test_nonblocking_exchange(array)
  !
  !.. Formal Arguments ..
  real(vbp), dimension(:), intent(inout) :: array
  !
  !.. Local Scalars..
  integer :: i,n,s,r,dtype,ierr
  !
  !.. Local Arrays ..
  integer, dimension(1:2*ncomm) :: rqst
  integer, dimension(1:mpi_status_size,1:2*ncomm) :: istat
  !
  !.. Local Allocatable Arrays ..
  real(vbp), allocatable, dimension(:) :: buffer
  !
continue
  !
  allocate ( buffer(1:size(array)) , source=array , stat=ierr )
  call alloc_error("test_nonblocking_exchange","buffer",1,ierr)
  !
  ! Initialize all the receives so each process
  ! is ready once the send is posted
  !
  do n = 1,ncomm
    s = exch_cpu(n)
    dtype = test_recv(n)
    call mpi_irecv(buffer,1,dtype,s,s,mpi_comm_world,rqst(n),mpierr)
  end do
  !
  ! Send out the adjacent cells
  !
  do n = 1,ncomm
    r = exch_cpu(n)
    dtype = test_send(n)
    call mpi_isend(array,1,dtype,r,mypnum,mpi_comm_world,rqst(n+ncomm),mpierr)
  end do
  !
  ! Wait for all sends and receives to complete
  !
  call mpi_waitall(2*ncomm,rqst,istat,mpierr)
  !
  ! Update the adjacent elements
  !
  array(:) = buffer(:)
  !
  ! Deallocate buffer before we leave
  !
  deallocate ( buffer , stat=ierr )
  call alloc_error("test_nonblocking_exchange","buffer",2,ierr)
  !
end subroutine test_nonblocking_exchange
!
!###############################################################################
!
subroutine collect_global_solpts
  !
  !.. Use Statements ..
  use geovar, only : n_solpts,xyz,global_xyz
  !
  !.. Local Scalars ..
  integer :: i,icpu,n,ierr
  !
  !.. Local Arrays ..
  integer, dimension(1:mpi_status_size) :: istat
  !
  !.. Local Allocatable Arrays ..
  real(r8), allocatable, dimension(:,:) :: xyz_r8
  real(r8), allocatable, dimension(:,:) :: global_xyz_r8
  !
continue
  !
  if (ncpu > 1) then
    !
    allocate ( xyz_r8(1:size(xyz,1),1:size(xyz,2)) , source=real(xyz,r8) , stat=ierr )
    call alloc_error("collect_global_solpts","xyz_r8",1,ierr)
    !
    allocate ( global_xyz_r8(1:size(global_xyz,1),1:size(global_xyz,2)) , &
               source=0.0_r8 , stat=ierr )
    call alloc_error("collect_global_solpts","global_xyz_r8",1,ierr)
    !
    n = size(xyz(:,1:n_solpts))
    !
    if (mypnum /= 0) then
      !
      call mpi_send(xyz_r8(:,1:n_solpts),n,mpi_flttyp,0,mypnum,mpi_comm_world,mpierr)
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
    global_xyz(:,:) = real( global_xyz_r8(:,:) , vbp )
    !
    deallocate ( xyz_r8 , stat=ierr )
    call alloc_error("collect_global_solpts","xyz_r8",2,ierr)
    deallocate ( global_xyz_r8 , stat=ierr )
    call alloc_error("collect_global_solpts","global_xyz_r8",2,ierr)
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
  use flowvar, only : vsp,global_vsp
  !
  !.. Local Scalars ..
  integer :: i,icpu,n,ierr
  !
  !.. Local Arrays ..
  integer, dimension(1:mpi_status_size) :: istat
  !
  !.. Local Allocatable Arrays ..
  real(r8), allocatable, dimension(:,:) :: sp_r8
  real(r8), allocatable, dimension(:,:) :: global_r8
  !
continue
  !
  if (ncpu > 1) then
    !
    allocate ( sp_r8(1:size(usp,1),1:size(usp,2)) , source=real(usp,r8) , stat=ierr )
    call alloc_error("collect_global_solution","sp_r8",1,ierr)
    allocate ( global_r8(1:size(global_usp,1),1:size(global_usp,2)) , &
               source=0.0_r8 , stat=ierr )
    call alloc_error("collect_global_solution","global_r8",1,ierr)
    !
    ! Collect the global usp array on all processors
    !
    n = size(sp_r8(:,1:n_solpts))
    !
    if (mypnum /= 0) then
      !
      call mpi_send(sp_r8(:,1:n_solpts),n,mpi_flttyp,0,mypnum,mpi_comm_world,mpierr)
      !
    else
      !
      call mpi_sendrecv(sp_r8(:,1:n_solpts),n,mpi_flttyp,0,n, &
                        global_r8,1,recv_sol(1),0,n, &
                        mpi_comm_world,istat,mpierr)
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
    end if
    !
    call mpi_bcast(global_r8,size(global_r8),mpi_flttyp,0,mpi_comm_world,mpierr)
    !
    global_usp(:,:) = real( global_r8(:,:) , vbp )
    !
    ! Collect the global vsp array on all processors
    !
    sp_r8(:,:) = real( vsp(:,:) , r8 )
    global_r8(:,:) = 0.0_r8
    !
    n = size(sp_r8(:,1:n_solpts))
    !
    if (mypnum /= 0) then
      !
      call mpi_send(sp_r8(:,1:n_solpts),n,mpi_flttyp,0,mypnum,mpi_comm_world,mpierr)
      !
    else
      !
      call mpi_sendrecv(sp_r8(:,1:n_solpts),n,mpi_flttyp,0,n, &
                        global_r8,1,recv_sol(1),0,n, &
                        mpi_comm_world,istat,mpierr)
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
    end if
    !
    call mpi_bcast(global_r8,size(global_r8),mpi_flttyp,0,mpi_comm_world,mpierr)
    !
    global_vsp(:,:) = real( global_r8(:,:) , vbp )
    !
    deallocate ( sp_r8 , stat=ierr )
    call alloc_error("collect_global_solution","sp_r8",2,ierr)
    deallocate ( global_r8 , stat=ierr )
    call alloc_error("collect_global_solution","global_r8",2,ierr)
    !
  else
    !
    global_usp(:,:) = usp(:,:)
    global_vsp(:,:) = vsp(:,:)
    !
  end if
  !
end subroutine collect_global_solution
!
!###############################################################################
!###############################################################################
!###############################################################################
!
! Routines to expand an array while keeping its data intact. This module
! contains a generic interface for these two routines to allow use through
! a single routine name for either a rank one or rank two array.
!
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine reallocate_1(array,new_size)
  !
  !.. Formal Arguments ..
  integer,                               intent(in) :: new_size
  integer, allocatable, dimension(:), intent(inout) :: array
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:) :: temp_array
  !
  !.. Local Scalars ..
  integer :: ierr
  !
continue
  !
  ! Allocate the temporary array to store the old contents of array
  !
  allocate ( temp_array(1:new_size) , source=0 , stat=ierr )
  call alloc_error("reallocate_1","temp_array",1,ierr)
  !
  ! If array is allocated, copy the old contents
  ! of array into the temporary array
  !
  if (allocated(array)) temp_array(1:size(array)) = array
  !
  ! move_alloc does:
  !     1. deallocates array if it is allocated
  !     2. allocates array with the same type, type parameters,
  !        array bounds, and value as temp_array
  !     3. deallocates temp_array
  !
  call move_alloc( from=temp_array , to=array )
  !
end subroutine reallocate_1
!
!###############################################################################
!
subroutine reallocate_2(array,new_size1,new_size2)
  !
  !.. Formal Arguments ..
  integer, optional,                       intent(in) :: new_size1
  integer, optional,                       intent(in) :: new_size2
  integer, allocatable, dimension(:,:), intent(inout) :: array
  !
  !.. Local Scalars ..
  integer :: old_size1,old_size2,ierr
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:,:) :: temp_array
  !
continue
  !
  ! Get the old dimensions of array
  !
  old_size1 = 1
  old_size2 = 1
  if (allocated(array)) then
    old_size1 = size(array,1)
    old_size2 = size(array,2)
  end if
  !
  ! Allocate the temporary array to store the old contents of array
  !
  if (present(new_size1) .and. present(new_size2)) then
    allocate ( temp_array(1:new_size1,1:new_size2) , source=0 , stat=ierr )
  else if (present(new_size1)) then
    allocate ( temp_array(1:new_size1,1:old_size2) , source=0 , stat=ierr )
  else if (present(new_size2)) then
    allocate ( temp_array(1:old_size1,1:new_size2) , source=0 , stat=ierr )
  else
    allocate ( temp_array(1:old_size1,1:old_size2) , source=0 , stat=ierr )
  end if
  call alloc_error("reallocate_2","temp_array",1,ierr)
  !
  ! Copy the old contents of array into the temporary array
  !
  if (allocated(array)) temp_array(1:old_size1,1:old_size2) = array
  !
  ! move_alloc does:
  !     1. deallocates array
  !     2. allocates array with the same type, type parameters,
  !        array bounds, and value as temp_array
  !     3. deallocates temp_array
  !
  call move_alloc( from=temp_array , to=array )
  !
end subroutine reallocate_2
!
!###############################################################################
!
subroutine output_partitions(id,nodes_of_cell,nodes_of_cell_ptr,xyz_nodes,epart)
  !
  !.. Formal Arguments ..
  integer,                           intent(in) :: id
  integer,   dimension(:),           intent(in) :: nodes_of_cell
  integer,   dimension(:),           intent(in) :: nodes_of_cell_ptr
  real(vbp), dimension(:,:),         intent(in) :: xyz_nodes
  integer,   dimension(:), optional, intent(in) :: epart
  !
  !.. Local Scalars ..
  integer :: i,n,n1,n2,ierr,np,nc
  character(len=100) :: fname
  !
  !.. Local Constants ..
  integer, parameter :: iot = 103
  !
continue
  !
  if (present(epart)) then
    fname = "partitions.tec"
    np = ngnode
    nc = ngcell
  else
    write (fname,'("partitions.",i0,".tec")') id
    np = size(xyz_nodes,dim=2)
    nc = size(nodes_of_cell_ptr)-1
  end if
  fname = adjustl(fname)
  !
  open (iot,file=trim(fname),status="replace",action="write",iostat=ierr)
  call io_error("output_partition",trim(fname),1,ierr)
  !
  ! Tecplot header with grid size parameters
  !
  if (present(epart)) then
    write (iot,1) np,nc
  else
    write (iot,5) np,nc
  end if
  !
  ! Node coordinates
  !
  do n = 1,size(xyz_nodes,dim=1)
    write (iot,2) (xyz_nodes(n,i),i=1,np)
  end do
  !
  ! Cell partition 
  !
  if (present(epart)) write (iot,3) (epart(i),i=1,nc)
  !
  ! Cell Connectivity
  !
  do n = 1,nc
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    if (n2-n1+1 == 3) then
      write (iot,4) (nodes_of_cell(i),i=n1,n2),nodes_of_cell(n2)
    else
      write (iot,4) (nodes_of_cell(i),i=n1,n2)
    end if
  end do
  !
  close (iot,iostat=ierr)
  call io_error("output_partitions",trim(fname),2,ierr)
  !
  5 format ('VARIABLES = "X" "Y"',/, &
            'ZONE T="Partitioned Grid", N=',i0,', E=',i0, &
            ', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL')
  1 format ('VARIABLES = "X" "Y" "Partition"',/, &
            'ZONE T="Partitioned Grid", N=',i0,', E=',i0, &
            ', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL',/, &
            'VARLOCATION = ([3]=CellCentered)')
  2 format (6(1x,es15.8))
  3 format (8(4x,i0))
  4 format (4(3x,i0))
  !
end subroutine output_partitions
!
!###############################################################################
!
include "Functions/last.f90"
!
!###############################################################################
!
end module parallel_mod
