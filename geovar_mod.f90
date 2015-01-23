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
  type, public :: cell_type
    integer :: geom
    integer :: order
  end type cell_type
  !
  type(cell_type), public, save, allocatable :: cell(:)
  !
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
  ! cells_on_face(1,i_face) : attached cell to the left  of i_face
  ! cells_on_face(2,i_face) : attached cell to the right of i_face
  ! cells_on_face(3,i_face) : node of left  triangular cell not on i_face (quads = 0)
  ! cells_on_face(4,i_face) : node of right triangular cell not on i_face (quads = 0)
  ! cells_on_face(5,i_face) : index of correction function array needed for left  cell
  ! cells_on_face(6,i_face) : index of correction function array needed for right cell
  !
  ! For Triangles : 
  !   for cells_on_face(3,:) = (    1    ,     2     ,3)
  !       cells_on_face(5,:) = (n_sp1de+1,2*n_sp1de+1,1)
  !   for cells_on_face(4,:) = (    1    ,    2    ,   3   )
  !       cells_on_face(6,:) = (2*n_sp1de,3*n_sp1de,n_sp1de)
  ! For Quadrilaterals : 
  !   for cells_on_face(3,:) = (    1    ,     2     ,     3     ,4)
  !       cells_on_face(5,:) = (n_sp1de+1,2*n_sp1de+1,3*n_sp1de+1,1)
  !   for cells_on_face(4,:) = (    1    ,    2    ,    3    ,   4   )
  !       cells_on_face(6,:) = (2*n_sp1de,3*n_sp1de,4*n_sp1de,n_sp1de)
  !
  integer, public, save, allocatable, dimension(:,:) :: nodes_on_face
  integer, public, save, allocatable, dimension(:,:) :: cells_on_face
  !
  ! For cell i_cell, 
  !    flxpts_of_cell_ptr(i_cell)+1 = first flux point for i_cell
  !    flxpts_of_cell_ptr(i_cell+1) = last  flux point for i_cell
  ! For flux point i_flxpt of i_cell
  !    flxpts_of_cell(1,i_flxpt) = grid face that i_flxpt is on
  !    flxpts_of_cell(2,i_flxpt) = index of i_flxpt on this face from 1:n_sp1de
  !    flxpts_of_cell(3,i_flxpt) = index of lagrange polynomial for interpolating
  !                                i_cell solution point data to i_flxpt
  !    flxpts_of_cell(4,i_flxpt) = gives the side of the face i_cell lies on
  !
  integer, public, parameter :: nflxprm = 4
  !
  integer, public, save, allocatable, dimension(:)   :: flxpts_of_cell_ptr
  integer, public, save, allocatable, dimension(:,:) :: flxpts_of_cell
  !
  !
  type, public :: geom_conn
    integer, allocatable, dimension(:,:) :: connectivity
  end type geom_conn
  !
  type(geom_conn), public, save, allocatable, dimension(:,:) :: pst_geom
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
  if (allocated(cells_on_face)) then
    array%siz = size(cells_on_face)
    array%mem = array%siz*storage_size(cells_on_face)
    total%mem = total%mem + array%mem
    write (iunit,2) "'cells_on_face'",array%get_units()
  end if
  if (allocated(flxpts_of_cell_ptr)) then
    array%siz = size(flxpts_of_cell_ptr)
    array%mem = array%siz*storage_size(flxpts_of_cell_ptr)
    total%mem = total%mem + array%mem
    write (iunit,2) "'flxpts_of_cell_ptr'",array%get_units()
  end if
  if (allocated(flxpts_of_cell)) then
    array%siz = size(flxpts_of_cell)
    array%mem = array%siz*storage_size(flxpts_of_cell)
    total%mem = total%mem + array%mem
    write (iunit,2) "'flxpts_of_cell'",array%get_units()
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
