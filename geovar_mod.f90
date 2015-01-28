module geovar
  !
  use module_kind_types, only : wp,Geom_Min,Geom_Max
  !
  implicit none
  !
  private
  !
  public :: geovar_alloc_status
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
  logical, public, save :: using_quadrature
  !
  ! nbfai : number of storage locations needed for each boundary face in array bface
  !
  integer, public :: nbfai = 10
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
  integer, public, save :: nfbnd
  integer, public, save :: n_totcel
  integer, public, save :: n_totnod
  integer, public, save :: n_totpts
  !
  ! For cell i_cell
  ! pts_in_cell_ptr(i_cell)+1 gives the first solution point within i_cell
  ! pts_in_cell_ptr(i_cell+1) gives the last  solution point within i_cell
  !
  integer, public, save, allocatable, dimension(:) :: pts_in_cell_ptr
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
  !                      for triangles: 1 -> nodes 1 and 2 of cell are on the boundary
  !                      for triangles: 2 -> nodes 2 and 3 of cell are on the boundary
  !                      for triangles: 3 -> nodes 3 and 1 of cell are on the boundary
  !                      for quads: 1 -> nodes 1 and 2 of cell are on the boundary
  !                      for quads: 2 -> nodes 2 and 3 of cell are on the boundary
  !                      for quads: 3 -> nodes 3 and 4 of cell are on the boundary
  !                      for quads: 4 -> nodes 4 and 1 of cell are on the boundary
  !        bface( 5,:) = A) For a non-communication (i.e., physical boundary)
  !                         boundary face i_face, this gives the boundary condition
  !                         group to which i_face belongs. This is used to group
  !                         together boundary faces for analysis during the
  !                         output/post-processing parts of the simulation.
  !                      B) For a periodic or partition boundary face i_face with
  !                         host cell i_cell, this gives the boundary face j_face
  !                         that links to i_face. To find host cell j_cell for
  !                         periodic face j_face, use bface(2,bface(5,i_face)).
  !        bface( 6,:) = For a periodic or partition boundary face i_face
  !                      with host cell i_cell, this gives the host cell
  !                      j_cell (i.e., bface(2,bface(5,i_face)) if they are
  !                      located on the same partition). This is needed if
  !                      they arent on the same partition because the local
  !                      process wont have access to the host cell information
  !                      in bface(2,j_face) for the adjacent partition
  !                      which is actually stored on another processor.
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
    integer :: flux_point
    integer :: face
    integer :: cpu
  end type point
  !
  type, public :: conflict
    logical :: is_local ! true if conflict not on a partition boundary
    logical :: not_my_conflict ! true if processor has no involvement with conflict
    type(point) :: master
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
  ! For face i_face, with the face vector going from
  ! nodes_on_face(1,i_face) to nodes_on_face(2,i_face), ...
  !
  integer, public, save, allocatable, dimension(:,:) :: nodes_on_face
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
  !
  !
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
    integer, allocatable :: flxpts(:)
  end type cell_on_face_t
  !
  ! FACE_T : derived type to store information for each grid face
  !
  type, public :: face_t
    ! geom : geometry type of the current face
    integer :: geom
    ! order : solution order of the current face
    !         NOTE: This is should be the maximum solution order between the two
    !               cells attached to the current face.
    integer :: order
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
    ! face_idx : index of flx-pt relative to all flx-pts on the current face
    integer, allocatable :: face_idx(:)
    ! cell_idx : index of flx-pt relative to all flx-pts of the current cell
    integer, allocatable :: cell_idx(:)
  end type face_of_cell_t
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
    ! face : array giving the flx-pt information for all faces of the current cell
    type(face_of_cell_t), allocatable :: face(:)
  end type cell_t
  !
  ! CELL : Array of type cell_t that contains information for all grid cells
  !
  type(cell_t), public, save, allocatable :: cell(:)
  !
contains
!
!###############################################################################
!
subroutine geovar_alloc_status(iunit)
  !
  !.. Use Statements ..
  use module_memory
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iunit
  !
  !.. Local Scalars ..
  type(memory) :: array
  type(memory) :: total
  !
continue
  !
  total%siz = 0
  total%mem = 0
  !
  write (iunit,1)
  !
  if (allocated(pts_in_cell_ptr)) then
    array%siz = size(pts_in_cell_ptr)
    array%mem = array%siz*storage_size(pts_in_cell_ptr)
    total%mem = total%mem + array%mem
    write (iunit,2) "'pts_in_cell_ptr'",array%get_units()
  end if
  if (allocated(xyz_nodes)) then
    array%siz = size(xyz_nodes)
    array%mem = array%siz*storage_size(xyz_nodes)
    total%mem = total%mem + array%mem
    write (iunit,2) "'xyz_nodes'",array%get_units()
  end if
  if (allocated(xyz)) then
    array%siz = size(xyz)
    array%mem = array%siz*storage_size(xyz)
    total%mem = total%mem + array%mem
    write (iunit,2) "'xyz'",array%get_units()
  end if
  if (allocated(global_xyz_nodes)) then
    array%siz = size(global_xyz_nodes)
    array%mem = array%siz*storage_size(global_xyz_nodes)
    total%mem = total%mem + array%mem
    write (iunit,2) "'global_xyz_nodes'",array%get_units()
  end if
  if (allocated(global_xyz)) then
    array%siz = size(global_xyz)
    array%mem = array%siz*storage_size(global_xyz)
    total%mem = total%mem + array%mem
    write (iunit,2) "'global_xyz'",array%get_units()
  end if
  if (allocated(global_nofc_ptr)) then
    array%siz = size(global_nofc_ptr)
    array%mem = array%siz*storage_size(global_nofc_ptr)
    total%mem = total%mem + array%mem
    write (iunit,2) "'global_nofc_ptr'",array%get_units()
  end if
  if (allocated(global_nofc)) then
    array%siz = size(global_nofc)
    array%mem = array%siz*storage_size(global_nofc)
    total%mem = total%mem + array%mem
    write (iunit,2) "'global_nofc'",array%get_units()
  end if
  if (allocated(global_pinc_ptr)) then
    array%siz = size(global_pinc_ptr)
    array%mem = array%siz*storage_size(global_pinc_ptr)
    total%mem = total%mem + array%mem
    write (iunit,2) "'global_pinc_ptr'",array%get_units()
  end if
  if (allocated(global_cell_order)) then
    array%siz = size(global_cell_order)
    array%mem = array%siz*storage_size(global_cell_order)
    total%mem = total%mem + array%mem
    write (iunit,2) "'global_cell_order'",array%get_units()
  end if
  if (allocated(global_cell_geom)) then
    array%siz = size(global_cell_geom)
    array%mem = array%siz*storage_size(global_cell_geom)
    total%mem = total%mem + array%mem
    write (iunit,2) "'global_cell_geom'",array%get_units()
  end if
  if (allocated(bface)) then
    array%siz = size(bface)
    array%mem = array%siz*storage_size(bface)
    total%mem = total%mem + array%mem
    write (iunit,2) "'bface'",array%get_units()
  end if
  if (allocated(bc_conflict)) then
    array%siz = size(bc_conflict)
    array%mem = array%siz*storage_size(bc_conflict)
    total%mem = total%mem + array%mem
    write (iunit,2) "'bc_conflict'",array%get_units()
  end if
  if (allocated(cell_order)) then
    array%siz = size(cell_order)
    array%mem = array%siz*storage_size(cell_order)
    total%mem = total%mem + array%mem
    write (iunit,2) "'cell_order'",array%get_units()
  end if
  if (allocated(cell_geom)) then
    array%siz = size(cell_geom)
    array%mem = array%siz*storage_size(cell_geom)
    total%mem = total%mem + array%mem
    write (iunit,2) "'cell_geom'",array%get_units()
  end if
  if (allocated(nodes_of_cell_ptr)) then
    array%siz = size(nodes_of_cell_ptr)
    array%mem = array%siz*storage_size(nodes_of_cell_ptr)
    total%mem = total%mem + array%mem
    write (iunit,2) "'nodes_of_cell_ptr'",array%get_units()
  end if
  if (allocated(nodes_of_cell)) then
    array%siz = size(nodes_of_cell)
    array%mem = array%siz*storage_size(nodes_of_cell)
    total%mem = total%mem + array%mem
    write (iunit,2) "'nodes_of_cell'",array%get_units()
  end if
  if (allocated(nodes_on_face)) then
    array%siz = size(nodes_on_face)
    array%mem = array%siz*storage_size(nodes_on_face)
    total%mem = total%mem + array%mem
    write (iunit,2) "'nodes_on_face'",array%get_units()
  end if
  !
  total = total%get_units()
  write (iunit,3) total%mem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module geovar")
  2 format (" Array: ",a21,dt)
  3 format (" Total memory = ",i0,1x,a)
  !
end subroutine geovar_alloc_status
!
!###############################################################################
!
end module geovar
!
!###############################################################################
!
