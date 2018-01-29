module geovar
  !
  use module_kind_types, only : wp,lk,iout,int_mpi
  use module_kind_types, only : Geom_Min,Geom_Max,Geom_Name
  use module_kind_types, only : geom_family_t,Complete
  !
  implicit none
  !
  private
  !
  public :: geovar_memory_usage
  !
  !  c - cells
  !  v - vertices
  !  nbr - neighbors
  !  ic - interior cells
  !  e - edges
  !  be - boundary edges
  !  idbe - identity of boundary edge
  !  pc - periodic cell
  !  eunl - edge unit normal and (edge) length (rotate edge 90 degrees).
  !
  !  allocate in subroutine ncc2cese
  !
  ! n_solpts : total number of solution points across all cells
  ! n_flxpts : total number of flux points across all cells
  !
  integer, public, save :: n_solpts
  integer, public, save :: n_flxpts
  !
  ! over-integration parameters
  !
  logical(lk), public, save :: using_quadrature
  !
  ! nbfai : number of storage locations needed for
  !         each boundary face in array bface
  !
  integer, parameter, public :: nbfai = 10
  !
  ! nr    : number of dimensions for computation
  ! nnode : total number of grid nodes
  ! ncell : total number of grid cells
  ! nface : total number of edges/faces
  ! nfbnd : total number of boundary edges/faces
  !
  integer, public, save :: nr
  integer, public, save :: nnode
  integer, public, save :: ncell
  integer, public, save :: nface
  integer, public, save :: nedge
  integer, public, save :: nfbnd
  integer, public, save :: n_totcel
  integer, public, save :: n_totnod
  integer, public, save :: n_totpts
  !
  ! cell_geom : this gives the geometry of each cell
  !
  integer, public, save, allocatable, dimension(:) :: cell_geom
  !
  ! cell_order : gives the order of the solution for each cell
  !              This if for future use if we want the solution order
  !              to not be constant across all cells
  !
  integer, public, save, allocatable, dimension(:) :: cell_order
  !
  ! xyz_nodes : x,y coordinates of the vertices
  ! xyz       : x,y coordinates of all solution points
  !
  real(wp), public, save, allocatable, dimension(:,:) :: xyz_nodes
  real(wp), public, save, allocatable, dimension(:,:) :: xyz
  !
  ! Grid size parameters for global grid
  !
  integer, public, save :: n_global_node
  integer, public, save :: n_global_totnod
  integer, public, save :: n_global_cell
  integer, public, save :: n_global_totcel
  integer, public, save :: n_global_edge
  integer, public, save :: n_global_fbnd
  integer, public, save :: n_global_solpts
  integer, public, save :: n_global_totpts
  !
  ! Arrays for the global grid
  !
  real(wp), public, save, allocatable, dimension(:,:) :: global_xyz_nodes
  real(wp), public, save, allocatable, dimension(:,:) :: global_xyz
  integer,  public, save, allocatable, dimension(:)   :: global_nofc_ptr
  integer,  public, save, allocatable, dimension(:)   :: global_nofc
  integer,  public, save, allocatable, dimension(:)   :: global_pinc_ptr
  integer,  public, save, allocatable, dimension(:)   :: global_cell_geom
  integer,  public, save, allocatable, dimension(:)   :: global_cell_order
  !
  ! bface : Boundary face information
  !        bface( 1,:) = Type of boundary condition
  !        bface( 2,:) = Host cell of boundary face
  !        bface( 3,:) = Type of cell on boundary
  !                       = 1 -> edge
  !                       = 2 -> quad
  !                       = 3 -> triangle
  !                       = 4 -> tetrahedra
  !                       = 5 -> pyramid
  !                       = 6 -> prism
  !                       = 8 -> hexahedra
  !        bface( 4,:) = Face of host cell on boundary
  !                for triangles: 1 -> nodes 1 and 2 of cell are on the boundary
  !                for triangles: 2 -> nodes 2 and 3 of cell are on the boundary
  !                for triangles: 3 -> nodes 3 and 1 of cell are on the boundary
  !                for quads: 1 -> nodes 1 and 2 of cell are on the boundary
  !                for quads: 2 -> nodes 2 and 3 of cell are on the boundary
  !                for quads: 3 -> nodes 3 and 4 of cell are on the boundary
  !                for quads: 4 -> nodes 4 and 1 of cell are on the boundary
  !        bface( 5,:) = A) For a non-communication (i.e., physical boundary)
  !                         boundary face i_face, this gives the boundary
  !                         condition group to which i_face belongs. This is
  !                         used to group together boundary faces for analysis
  !                         during the output/post-processing parts of the
  !                         simulation.
  !                      B) For a periodic or partition boundary face i_face
  !                         with host cell i_cell, this gives the boundary
  !                         face j_face that links to i_face. To find host cell
  !                         j_cell for periodic face j_face, use
  !                         bface(2,bface(5,i_face)).
  !        bface( 6,:) = A) For a non-communication (i.e., physical boundary)
  !                         boundary face i_face, this gives the index of
  !                         the array bc_in(:) that provides the flow
  !                         conditions for this boundary condition.
  !                      B) For a periodic or partition boundary face i_face
  !                         with host cell i_cell, this gives the host cell
  !                         j_cell (i.e., bface(2,bface(5,i_face)) if they are
  !                         located on the same partition). This is needed if
  !                         they arent on the same partition because the local
  !                         process wont have access to host cell information
  !                         in bface(2,j_face) for the adjacent partition
  !                         which is actually stored on another processor.
  !        bface( 7,:) = For a periodic or partition boundary face i_face
  !                      with host cell i_cell, this gives the partition that
  !                      j_cell belongs to (i.e., bface(2,bface(5,i_face))
  !                      if they are on the same partition).
  !        bface( 8,:) = Geometry type of host cell on adjacent partition
  !                      Equivalent to bface(3,:) but for adjacent partition
  !                      host cell
  !        bface( 9,:) = Face of adjacent partition host cell on boundary.
  !                      Equivalent to bface(4,:) but for adjacent partition
  !                      host cell.
  !        bface(10,:) = Index of ghost cell for current boundary face
  !                      NOTE: For a grid with only one partition, this is
  !                            simply bface(10,i_face) = ncell + i_face.
  !                            However, for multiple partitions, two boundary
  !                            faces can connect to the same cell on an adjacent
  !                            partition and therefore connect to the same ghost
  !                            cell.
  !
  integer, public, save, allocatable, dimension(:,:) :: bface
  !
  ! CONFLICT : A data type that provides the node and face indices for any
  !            boundary condition conflicts. A boundary condition conflict
  !            is where two adjacent boundary faces have non-matching /
  !            conflicting boundary conditions. These conflicts are located
  !            at a single node for a 2D grid and along an edge for a 3D grid.
  !            NOTE: This data type might need to be tweaked for 3D. The whole
  !                  point of it is to exchange the computed solution at a
  !                  boundary flux point so that they are consistent. This is
  !                  straight forward for 2D because we only need to exchange
  !                  this information for a single point. However, for 3D,
  !                  we will need to exchange information along the entire
  !                  edge meaning we will have multiple points instead of
  !                  just the one. Example, add the allocatable component
  !                  flxpts as a possible solution to this potential flaw.
  !                  For 2D, flxpts is allocated to flxpts(1:2,1:1) where
  !                  the first rank is for which face the flux point is
  !                  from while the second rank gives the actual flux point.
  !                  ??? For 3D, this could be allocated to something like
  !                  flxpts(1:size(faces),1:n_sp1de) ???
  !
  type, public :: point
    integer          :: flux_point
    integer          :: face
    integer(int_mpi) :: cpu
  end type point
  !
  type, public :: conflict
    ! is_local : true if conflict not on a partition boundary
    logical(lk) :: is_local
    ! not_my_conflict :  true if processor has no involvement with conflict
    logical(lk) :: not_my_conflict
    ! master :
    type(point) :: master
    ! slave :
    type(point), allocatable, dimension(:) :: slave
  end type conflict
  !
  ! BC_CONFLICT : Array for storing all the boundary condition conflicts.
  !
  type(conflict), public, save, allocatable, dimension(:) :: bc_conflict
  !
  ! For cell i_cell,
  ! nodes_of_cell(ibeg:iend) gives the grid nodes that define i_cell
  ! where : ibeg = nodes_of_cell_ptr(i_cell)+1
  !         iend = nodes_of_cell_ptr(i_cell+1)
  !
  integer, public, save, allocatable, dimension(:) :: nodes_of_cell_ptr
  integer, public, save, allocatable, dimension(:) :: nodes_of_cell
  !
  !
  !
  ! GEOM_CONN :
  !
  type, public :: geom_conn
    integer, allocatable, dimension(:,:) :: connectivity
  end type geom_conn
  !
  type(geom_conn), public, save, allocatable, dimension(:,:) :: pst_geom
  !
  type, public :: fp_rot_pts_t
    integer, allocatable :: pts(:)
  end type fp_rot_pts_t
  !
  type, public :: fp_rot_t
    type(fp_rot_pts_t), allocatable :: rot(:)
  end type fp_rot_t
  !
  type, public :: fp_t
    type(fp_rot_t), allocatable :: geord(:,:)
  end type fp_t
  !
  type(fp_t), public, save, allocatable :: fp
  !
  ! CELL_ON_FACE_T : derived type to store the information about
  !                  a cell on one side of a given face
  !
  type, public :: cell_on_face_t
    ! cell : index of the grid cell on the current side of the face
    integer :: cell
    ! cell_face : index for the side/surface of the cell that comprises
    !             the current face
    !               Examples: for a triangle, in the range 1:3, quad 1:4,
    !                         tetr 1:4, pyra 1:5, pris 1:5, hexa 1:6
    integer :: cell_face
    ! fp_offset : offset to add to the flux point indices
    !             given by the rotation index
    integer :: fp_offset = 0
    ! rotation : index giving the rotation of the flux points
    !            for the current element on this face
    integer :: rotation
  contains
    private
    procedure :: get_fp_index => get_face_fp_index
    procedure :: get_fp_indices => get_face_fp_indices
    generic, public :: get_fp_idx => get_fp_index, get_fp_indices
  end type cell_on_face_t
  !
  ! FACE_T : derived type to store information for each grid face
  !
  type, public :: face_t
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
  type(face_t), public, save, allocatable :: face(:)
  !
  !
  !
  !
  ! FACE_OF_CELL_T : derived type to store the information about the flux points
  !                  on a given cell face
  !
  type, public :: face_of_cell_t
    ! idx : index of the current face relative to all grid faces (1:nface)
    integer :: idx
    ! side : side of face that the current cell is on (1=left,2=right)
    integer :: side
    ! geom : geometry type of the current face
    integer :: geom
    ! order : polynomial order of the current face
    integer :: order
    ! fp_offset : offset to add to the flux point indices
    !             given by the rotation index
    integer :: fp_offset
    ! rotation : index giving the rotation of the flux points
    !            for the current face of this element
    integer :: rotation
  contains
    private
    procedure :: get_fp_index => get_cell_fp_index
    procedure :: get_fp_indices => get_cell_fp_indices
    generic, public :: get_fp_idx => get_fp_index, get_fp_indices
  end type face_of_cell_t
  !
  !
  ! EDGE_OF_CELL_T : derived type to store the information about the solution
  !                  points on a given edge (this is only used for output)
  !
  type, public :: edge_of_cell_t
    ! idx : index of this edge within the edgeusp and edgedusp arrays
    integer :: idx
    ! order : polynomial order of the edge
    integer :: order
    ! ave_coef : gives the coefficient needed to compute an average
    !            i.e., ave_coef = 1/(number of cells that contain this edge)
    real(wp) :: ave_coef
    ! ep_offset : offset to add to the edge point indices
    integer :: ep_offset
    ! rev_dir : logical to indicate if the direction/rotation of the
    !           edge point indices needs to be reversed
    !  = fals : edge point indices are in increasing order and are given by
    !              intseq(order+1,rev_dir) = 1,...,order+1
    !  = true : edge point indices are in decreasing order and are given by
    !              intseq(order+1,rev_dir) = order+1,...,1
    logical(lk) :: rev_dir
  contains
    private
    procedure :: get_ep_index
    procedure :: get_ep_indices
    generic, public :: get_ep_idx => get_ep_index, get_ep_indices
  end type edge_of_cell_t
  !
  !
  ! NODE_OF_CELL_T : derived type to store the information about a node
  !                  (this is only used for output)
  !
  type, public :: node_of_cell_t
    ! idx : index of this node within the nodeusp and nodedusp arrays
    integer :: idx
    ! ave_coef : gives the coefficient needed to compute an average
    !            i.e., ave_coef = 1/(number of cells that contain this node)
    real(wp) :: ave_coef
    ! cell_idx : index of the node in the cell
    !            i.e., 1 <= cell_index <= geom_nodes(cell_geometry)
    integer :: cell_idx
  end type node_of_cell_t
  !
  !
  ! CELL_T : derived type to store information for each grid cell
  !
  type, public :: cell_t
    ! geom : geometry type of the current cell
    integer :: geom
    ! order : solution order of the current cell
    integer :: order
    ! beg_sp : index of first solution point in the current cell
    integer :: beg_sp
    ! end_sp : index of last  solution point in the current cell
    integer :: end_sp
    ! face : array giving the flx-pt information
    !        for all faces of the current cell
    type(face_of_cell_t), allocatable :: face(:)
    ! edge : array giving the edge-pt information
    !        for all edges of the current cell
    type(edge_of_cell_t), allocatable :: edge(:)
    ! node : array giving the information for all nodes of the current cell
    type(node_of_cell_t), allocatable :: node(:)
  end type cell_t
  !
  ! CELL : Array of type cell_t that contains information for all grid cells
  !
  type(cell_t), public, save, allocatable :: cell(:)
  !
  !
  ! EDGE_T : derived type to store information for each grid edge
  !
  type, public :: edge_t
    ! nodes : the two nodes that define this edge
    integer :: nodes(1:2)
    ! order : the maximum order of all cells containing this edge
    integer :: order
  end type edge_t
  !
  ! EDGE : Array of type edge_t that contains information for all grid edges
  !
  type(edge_t), public, save, allocatable :: edge(:)
  !
  ! CGNS_PATH_t
  !
 !type, public :: cgns_idx_t
 !  integer :: base
 !  integer :: zone
 !  integer :: section
 !  integer :: boco = 0
 !end type cgns_idx_t
  !
  ! PROP_T
  !
  type, public :: prop_t
    integer :: npts
    integer :: order
    integer :: geom
    integer :: elem_type
    type(geom_family_t) :: family = Complete
    logical(lk) :: warning = .false.
  end type prop_t
  !
  ! ELEM_T
  !
  type, public :: elem_t
    integer :: host_cell = 0
    type(prop_t) :: prop
    integer, allocatable :: pts(:)
    integer, allocatable :: tags(:)
    integer, allocatable :: nodes(:)
   !type(cgns_idx_t), allocatable :: idx
  end type elem_t
  !
  ! GRID_T
  !
  type, public :: grid_t
    type(elem_t), allocatable :: elem(:)
    real(wp),     allocatable :: xyz(:,:)
  end type grid_t
  !
  type(grid_t), public, save, allocatable :: grid
  !
contains
!
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function get_ep_indices(this) result(return_value)
  !.. Use Statements ..
  use module_kind_types, only : intseq
  !.. Passed Argument from Invoking Object ..
  class(edge_of_cell_t), intent(in) :: this
  !.. Function Result ..
  integer :: return_value(1:this%order+1)
continue
  return_value = this%ep_offset + intseq(this%order+1,this%rev_dir)
end function get_ep_indices
!
pure function get_ep_index(this,k) result(return_value)
  !.. Passed Argument from Invoking Object ..
  class(edge_of_cell_t), intent(in) :: this
  !.. Formal Arguments ..
  integer, intent(in) :: k
  !.. Function Result ..
  integer :: return_value
continue
  return_value = this%ep_offset + merge(this%order+2-k,k,this%rev_dir)
end function get_ep_index
!
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function get_cell_fp_indices(this) result(return_value)
  !.. Passed Argument from Invoking Object ..
  class(face_of_cell_t), intent(in) :: this
  !.. Function Result ..
  integer, allocatable :: return_value(:)
continue
  ! Using F2003 auto-reallocation
  return_value = this%fp_offset + fp%geord(this%geom,this%order)% &
                                              rot(this%rotation)% &
                                  pts(:)
end function get_cell_fp_indices
!
pure function get_cell_fp_index(this,k) result(return_value)
  !.. Passed Argument from Invoking Object ..
  class(face_of_cell_t), intent(in) :: this
  !.. Formal Arguments ..
  integer, intent(in) :: k
  !.. Function Result ..
  integer :: return_value
continue
  return_value = this%fp_offset + fp%geord(this%geom,this%order)% &
                                              rot(this%rotation)% &
                                  pts(k)
end function get_cell_fp_index
!
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function get_face_fp_indices(this,face_geom,face_order) &
                           result(return_value)
  !.. Passed Argument from Invoking Object ..
  class(cell_on_face_t), intent(in) :: this
  !.. Formal Arguments ..
  integer, intent(in) :: face_geom
  integer, intent(in) :: face_order
  !.. Function Result ..
  integer, allocatable :: return_value(:)
continue
  ! Using F2003 auto-reallocation
  return_value = this%fp_offset + fp%geord(face_geom,face_order)% &
                                              rot(this%rotation)% &
                                  pts(:)
end function get_face_fp_indices
!
pure function get_face_fp_index(this,face_geom,face_order,k) &
                         result(return_value)
  !.. Passed Argument from Invoking Object ..
  class(cell_on_face_t), intent(in) :: this
  !.. Formal Arguments ..
  integer, intent(in) :: face_geom
  integer, intent(in) :: face_order
  integer, intent(in) :: k
  !.. Function Result ..
  integer :: return_value
continue
  return_value = this%fp_offset + fp%geord(face_geom,face_order)% &
                                              rot(this%rotation)% &
                                  pts(k)
end function get_face_fp_index
!
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine geovar_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
  !
  !.. Formal Arguments ..
  integer, optional, intent(in) :: iunit
  !
  !.. Local Scalars ..
  integer :: iou,i,j,k
  character(len=36) :: array_name
  type(memory) :: array,dt
  type(memory) :: face_idx
  type(memory) :: edge_idx
  type(memory) :: flxpts
  type(memory) :: left_fp
  type(memory) :: right_fp
  type(memory) :: total
  !
#ifndef DISABLE_DTIO
continue
  !
  iou = iout
  if (present(iunit)) iou = iunit
  !
  call total%reset
  call face_idx%reset
  call edge_idx%reset
  call flxpts%reset
  call left_fp%reset
  call right_fp%reset
  !
  write (iou,1)
  !
  !
  ! Check the memory of the arrays within this module
  !
  if (allocated(face)) then
    !
    !##########################################################
   !type, public :: cell_on_face_t
   !  integer :: cell
   !  integer :: cell_face
   !  integer, allocatable :: flxpts(:)
   !end type cell_on_face_t
   !type, public :: face_t
   !  integer :: geom
   !  integer :: order
   !  integer :: nodes(1:4)
   !  type(cell_on_face_t) :: left
   !  type(cell_on_face_t) :: right
   !end type face_t
   !type(face_t), public, save, allocatable :: face(:)
    !##########################################################
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
  if (allocated(cell)) then
    !
    !##########################################################
   !type, public :: face_of_cell_t
   !  integer :: idx
   !  integer :: side
   !  integer :: geom
   !  integer :: order
   !  integer, allocatable :: cell_idx(:)
   !end type face_of_cell_t
   !type, public :: edge_of_cell_t
   !  integer :: idx
   !  integer :: order
   !  real(wp) :: ave_coef
   !  integer, allocatable :: cell_idx(:)
   !end type edge_of_cell_t
   !type, public :: node_of_cell_t
   !  integer :: idx
   !  real(wp) :: ave_coef
   !  integer :: cell_idx
   !end type node_of_cell_t
   !type, public :: cell_t
   !  integer :: geom
   !  integer :: order
   !  integer :: beg_sp
   !  integer :: end_sp
   !  type(face_of_cell_t), allocatable :: face(:)
   !  type(edge_of_cell_t), allocatable :: edge(:)
   !  type(node_of_cell_t), allocatable :: node(:)
   !end type cell_t
   !type(cell_t), public, save, allocatable :: cell(:)
    !##########################################################
    !
    call array%set( storage_size(cell,kind=inttype) , &
                            size(cell,kind=inttype) )
    !
    write (array_name,5) "cell(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    do i = lbound(cell,dim=1),ubound(cell,dim=1)
      !
      call dt%add( storage_size(cell(i)%geom,kind=inttype) )
      call dt%add( storage_size(cell(i)%order,kind=inttype) )
      call dt%add( storage_size(cell(i)%beg_sp,kind=inttype) )
      call dt%add( storage_size(cell(i)%end_sp,kind=inttype) )
      !
      if (allocated(cell(i)%face)) then
        !
        call array%set( storage_size(cell(i)%face,kind=inttype) , &
                                size(cell(i)%face,kind=inttype) )
        call dt%add(array)
        !
        do j = lbound(cell(i)%face,dim=1),ubound(cell(i)%face,dim=1)
          !
          call dt%add( storage_size(cell(i)%face(j)%idx,kind=inttype) )
          call dt%add( storage_size(cell(i)%face(j)%side,kind=inttype) )
          call dt%add( storage_size(cell(i)%face(j)%geom,kind=inttype) )
          call dt%add( storage_size(cell(i)%face(j)%order,kind=inttype) )
          !
          call array%set( storage_size(cell(i)%face(j)%fp_offset,kind=inttype) )
          call array%add( storage_size(cell(i)%face(j)%rotation,kind=inttype) )
          call face_idx%add(array)
          call dt%add(array)
          !
        end do
        !
      end if
      !
      if (allocated(cell(i)%edge)) then
        !
        call array%set( storage_size(cell(i)%edge,kind=inttype) , &
                                size(cell(i)%edge,kind=inttype) )
        call dt%add(array)
        !
        do j = lbound(cell(i)%edge,dim=1),ubound(cell(i)%edge,dim=1)
          !
          call dt%add( storage_size(cell(i)%edge(j)%idx,kind=inttype) )
          call dt%add( storage_size(cell(i)%edge(j)%order,kind=inttype) )
          call dt%add( storage_size(cell(i)%edge(j)%ave_coef,kind=inttype) )
          !
          call array%set( storage_size(cell(i)%edge(j)%ep_offset,kind=inttype) )
          call array%add( storage_size(cell(i)%edge(j)%rev_dir,kind=inttype) )
          call edge_idx%add(array)
          call dt%add(array)
          !
        end do
        !
      end if
      !
      if (allocated(cell(i)%node)) then
        !
        call array%set( storage_size(cell(i)%node,kind=inttype) , &
                                size(cell(i)%node,kind=inttype) )
        call dt%add(array)
        !
        do j = lbound(cell(i)%node,dim=1),ubound(cell(i)%node,dim=1)
          !
          call dt%add( storage_size(cell(i)%node(j)%idx,kind=inttype) )
          call dt%add( storage_size(cell(i)%node(j)%ave_coef,kind=inttype) )
          call dt%add( storage_size(cell(i)%node(j)%cell_idx,kind=inttype) )
          !
        end do
        !
      end if
      !
    end do
    !
    write (array_name,5) "cell%face%cell_idx"
    write (iou,6) array_name,face_idx%get_units()
    !
    write (array_name,5) "cell%edge%cell_idx"
    write (iou,6) array_name,edge_idx%get_units()
    !
    call total%add(dt)
    write (array_name,5) "cell"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "cell"
    !
  end if
  !
  ! Memory usage for the fp array of type fp_t
  !
  if (allocated(fp)) then
    !
    !##########################################################
   !type, public :: fp_rot_pts_t
   !  integer, allocatable :: pts
   !end type fp_rot_pts_t
   !type, public :: fp_rot_t
   !  type(fp_rot_pts_t), allocatable :: rot(:)
   !end type fp_rot_t
   !type, public :: fp_t
   !  type(fp_rot_t), allocatable :: geord(:,:)
   !end type fp_t
   !type(fp_t), public, save, allocatable :: fp
    !##########################################################
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
  ! Memory usage for the edge array of type edge_t
  !
  if (allocated(edge)) then
    !
    call array%set( storage_size(edge,kind=inttype) , &
                            size(edge,kind=inttype) )
    write (array_name,5) "edge(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    do i = lbound(edge,dim=1),ubound(edge,dim=1)
      !
      call dt%add( storage_size(edge(i)%order,kind=inttype) )
      !
      call array%set( storage_size(edge(i)%nodes,kind=inttype) , &
                              size(edge(i)%nodes,kind=inttype) )
      !
      call dt%add(array)
      !
    end do
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
  ! Memory usage for the grid variable of type grid_t
  !
  if (allocated(grid)) then
    !
    call dt%set( storage_size(grid,kind=inttype) , 1_inttype )
    !
    if (allocated(grid%elem)) then
      !
      call array%set( storage_size(grid%elem,kind=inttype) , &
                              size(grid%elem,kind=inttype) )
      !
      call dt%add(array)
      !
      do i = lbound(grid%elem,dim=1),ubound(grid%elem,dim=1)
        !
        call dt%add( storage_size(grid%elem(i)%host_cell,kind=inttype) )
        call dt%add( storage_size(grid%elem(i)%prop,kind=inttype) )
        call dt%add( storage_size(grid%elem(i)%prop%npts,kind=inttype) )
        call dt%add( storage_size(grid%elem(i)%prop%order,kind=inttype) )
        call dt%add( storage_size(grid%elem(i)%prop%geom,kind=inttype) )
        call dt%add( storage_size(grid%elem(i)%prop%elem_type,kind=inttype) )
        call dt%add( storage_size(grid%elem(i)%prop%family,kind=inttype) )
        call dt%add( storage_size(grid%elem(i)%prop%family%v,kind=inttype) )
        call dt%add( storage_size(grid%elem(i)%prop%warning,kind=inttype) )
        !
        if (allocated(grid%elem(i)%pts)) then
          call dt%add( storage_size(grid%elem(i)%pts,kind=inttype) , &
                               size(grid%elem(i)%pts,kind=inttype) )
        end if
        if (allocated(grid%elem(i)%tags)) then
          call dt%add( storage_size(grid%elem(i)%tags,kind=inttype) , &
                               size(grid%elem(i)%tags,kind=inttype) )
        end if
        if (allocated(grid%elem(i)%nodes)) then
          call dt%add( storage_size(grid%elem(i)%nodes,kind=inttype) , &
                               size(grid%elem(i)%nodes,kind=inttype) )
        end if
        !
      end do
      !
    end if
    !
    if (allocated(grid%xyz)) then
      call dt%add( storage_size(grid%xyz,kind=inttype) , &
                           size(grid%xyz,kind=inttype) )
    end if
    !
    call total%add(dt)
    write (array_name,5) "grid"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "grid"
    !
  end if

  !
  ! Memory usage for the pst_geom array of type geom_conn
  !
  if (allocated(pst_geom)) then
    !
    call array%set( storage_size(pst_geom,kind=inttype) , &
                            size(pst_geom,kind=inttype) )
    write (array_name,5) "pst_geom(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    do j = lbound(pst_geom,dim=2),ubound(pst_geom,dim=2)
      do i = lbound(pst_geom,dim=1),ubound(pst_geom,dim=1)
        !
        call array%set( &
                   storage_size(pst_geom(i,j)%connectivity,kind=inttype) , &
                           size(pst_geom(i,j)%connectivity,kind=inttype) )
        !
        write (array_name,8) Geom_Name(i),j
        write (iou,2) array_name,array%get_units()
        !
        call dt%add(array)
        !
      end do
    end do
    !
    call total%add(dt)
    write (array_name,5) "pst_geom"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "pst_geom"
    !
  end if
  !
  ! Memory usage for the pst_geom array of type geom_conn
  !
  if (allocated(bc_conflict)) then
    !
    call dt%set( storage_size(bc_conflict,kind=inttype) , &
                         size(bc_conflict,kind=inttype) )
    !
    do i = lbound(bc_conflict,dim=1),ubound(bc_conflict,dim=1)
      !
      call dt%add( storage_size(bc_conflict(i)%is_local,kind=inttype) )
      call dt%add( storage_size(bc_conflict(i)%not_my_conflict,kind=inttype) )
      call dt%add( storage_size(bc_conflict(i)%master,kind=inttype) )
      call dt%add( storage_size(bc_conflict(i)%master%flux_point,kind=inttype) )
      call dt%add( storage_size(bc_conflict(i)%master%face,kind=inttype) )
      call dt%add( storage_size(bc_conflict(i)%master%cpu,kind=inttype) )
      !
      if (allocated(bc_conflict(i)%slave)) then
        !
        call dt%add( storage_size(bc_conflict(i)%slave,kind=inttype) , &
                             size(bc_conflict(i)%slave,kind=inttype) )
        !
        do j = lbound(bc_conflict(i)%slave,dim=1), &
               ubound(bc_conflict(i)%slave,dim=1)
          call dt%add( storage_size(bc_conflict(i)%slave(j)%flux_point, &
                                    kind=inttype) )
          call dt%add( storage_size(bc_conflict(i)%slave(j)%face,kind=inttype) )
          call dt%add( storage_size(bc_conflict(i)%slave(j)%cpu,kind=inttype) )
        end do
        !
      end if
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "bc_conflict"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "bc_conflict"
    !
  end if
  !
  ! Memory usage for all the remaining arrays of intrinsic datatypes
  !
  !
  !
  if (allocated(cell_geom)) then
    call array%set( storage_size(cell_geom,kind=inttype) , &
                            size(cell_geom,kind=inttype) )
    call total%add(array)
    write (array_name,5) "cell_geom"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "cell_geom"
  end if
  !
  !
  !
  if (allocated(cell_order)) then
    call array%set( storage_size(cell_order,kind=inttype) , &
                            size(cell_order,kind=inttype) )
    call total%add(array)
    write (array_name,5) "cell_order"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "cell_order"
  end if
  !
  !
  !
  if (allocated(xyz_nodes)) then
    call array%set( storage_size(xyz_nodes,kind=inttype) , &
                            size(xyz_nodes,kind=inttype) )
    call total%add(array)
    write (array_name,5) "xyz_nodes"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "xyz_nodes"
  end if
  !
  !
  !
  if (allocated(xyz)) then
    call array%set( storage_size(xyz,kind=inttype) , &
                            size(xyz,kind=inttype) )
    call total%add(array)
    write (array_name,5) "xyz"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "xyz"
  end if
  !
  !
  !
  if (allocated(global_xyz_nodes)) then
    call array%set( storage_size(global_xyz_nodes,kind=inttype) , &
                            size(global_xyz_nodes,kind=inttype) )
    call total%add(array)
    write (array_name,5) "global_xyz_nodes"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "global_xyz_nodes"
  end if
  !
  !
  !
  if (allocated(global_xyz)) then
    call array%set( storage_size(global_xyz,kind=inttype) , &
                            size(global_xyz,kind=inttype) )
    call total%add(array)
    write (array_name,5) "global_xyz"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "global_xyz"
  end if
  !
  !
  !
  if (allocated(global_nofc_ptr)) then
    call array%set( storage_size(global_nofc_ptr,kind=inttype) , &
                            size(global_nofc_ptr,kind=inttype) )
    call total%add(array)
    write (array_name,5) "global_nofc_ptr"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "global_nofc_ptr"
  end if
  !
  !
  !
  if (allocated(global_nofc)) then
    call array%set( storage_size(global_nofc,kind=inttype) , &
                            size(global_nofc,kind=inttype) )
    call total%add(array)
    write (array_name,5) "global_nofc"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "global_nofc"
  end if
  !
  !
  !
  if (allocated(global_pinc_ptr)) then
    call array%set( storage_size(global_pinc_ptr,kind=inttype) , &
                            size(global_pinc_ptr,kind=inttype) )
    call total%add(array)
    write (array_name,5) "global_pinc_ptr"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "global_pinc_ptr"
  end if
  !
  !
  !
  if (allocated(global_cell_geom)) then
    call array%set( storage_size(global_cell_geom,kind=inttype) , &
                            size(global_cell_geom,kind=inttype) )
    call total%add(array)
    write (array_name,5) "global_cell_geom"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "global_cell_geom"
  end if
  !
  !
  !
  if (allocated(global_cell_order)) then
    call array%set( storage_size(global_cell_order,kind=inttype) , &
                            size(global_cell_order,kind=inttype) )
    call total%add(array)
    write (array_name,5) "global_cell_order"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "global_cell_order"
  end if
  !
  !
  !
  if (allocated(bface)) then
    call array%set( storage_size(bface,kind=inttype) , &
                            size(bface,kind=inttype) )
    call total%add(array)
    write (array_name,5) "bface"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "bface"
  end if
  !
  !
  !
  if (allocated(nodes_of_cell_ptr)) then
    call array%set( storage_size(nodes_of_cell_ptr,kind=inttype) , &
                            size(nodes_of_cell_ptr,kind=inttype) )
    call total%add(array)
    write (array_name,5) "nodes_of_cell_ptr"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "nodes_of_cell_ptr"
  end if
  !
  !
  !
  if (allocated(nodes_of_cell)) then
    call array%set( storage_size(nodes_of_cell,kind=inttype) , &
                            size(nodes_of_cell,kind=inttype) )
    call total%add(array)
    write (array_name,5) "nodes_of_cell"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "nodes_of_cell"
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module geovar")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a,:,a)
  4 format (a,"(",a,",",i0,")%",a,:,"(",i0,")%mat")
  5 format ("'",a,"'")
  6 format (" Array Component: ",a,dt)
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  8 format ("pst_geom(",a,",",i0,")%connectivity")
  !
#endif
end subroutine geovar_memory_usage
!
!###############################################################################
!
end module geovar
!
!###############################################################################
!
