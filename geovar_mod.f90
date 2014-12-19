module geovar
  !
  use module_kind_types, only : wp,Geom_Min,Geom_Max
  !
  implicit none
  !
  private
  !
  public :: write_variable
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
  ! gauss_1d       : coordinates of 1D Gauss solution points of order n_order
  ! gauss_wts_1d   : weights     of 1D Gauss solution points of order n_order
  ! gauss_wts_2d   : weights     of 2D Gauss solution points of order n_order
  !
  real(wp), public, save, allocatable, dimension(:) :: gauss_1d
  real(wp), public, save, allocatable, dimension(:) :: gauss_wts_1d
  real(wp), public, save, allocatable, dimension(:) :: gauss_wts_2d
  !
  ! lobatto_1d     : coordinates of 1D Lobatto solution points of order n_order
  ! lobatto_wts_1d : weights     of 1D Lobatto solution points of order n_order
  ! lobatto_wts_2d : weights     of 2D Lobatto solution points of order n_order
  !
  real(wp), public, save, allocatable, dimension(:) :: lobatto_1d
  real(wp), public, save, allocatable, dimension(:) :: lobatto_wts_1d
  real(wp), public, save, allocatable, dimension(:) :: lobatto_wts_2d
  !
  ! solpts_1d : coordinates of 1D solution points of order n_order
  ! flxpts_1d : coordinates of 1D flux points of order n_order
  ! solwts_1d : weights of the 1D solution points of order n_order
  ! flxwts_1d : weights of the 1D flux points of order n_order
  ! rs        : 2D Lobatto nodes of triangle in r-s space
  ! 
  real(wp), public, save, allocatable, dimension(:)   :: solpts_1d
  real(wp), public, save, allocatable, dimension(:)   :: flxpts_1d
  real(wp), public, save, allocatable, dimension(:)   :: solwts_1d
  real(wp), public, save, allocatable, dimension(:)   :: flxwts_1d
  real(wp), public, save, allocatable, dimension(:,:) :: rs
  !
  real(wp), public, save, allocatable, dimension(:)   :: qdtpts_1d
  real(wp), public, save, allocatable, dimension(:)   :: qdtwts_1d
  !
  ! gauss2numer   : matrix to move from Gauss   SP to the numerical SP
  ! lobatto2numer : matrix to move from Lobatto SP to the numerical SP
  !
  real(wp), public, save, allocatable, dimension(:,:) :: gauss2numer
  real(wp), public, save, allocatable, dimension(:,:) :: lobatto2numer
  !
  ! numer2gauss   : matrix to move from the numerical SP to Gauss   SP
  ! numer2lobatto : matrix to move from the numerical SP to Lobatto SP
  !
  real(wp), public, save, allocatable, dimension(:,:) :: numer2gauss
  real(wp), public, save, allocatable, dimension(:,:) :: numer2lobatto
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
  ! xyz_flxpts : x,y coordinates of all flux points
  !
  real(wp), public, save, allocatable, dimension(:,:,:) :: xyz_flxpts
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
  ! faces_of_cell : gives the global face index for each local face of each cell
  !
  integer, public, save, allocatable, dimension(:,:) :: faces_of_cell
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
  ! Various Vandermonde matrices
  !
  real(wp), public, save, allocatable, dimension(:,:,:) :: vand_inv
  !
  real(wp), public, save, allocatable, dimension(:,:) :: interp_f2g
  real(wp), public, save, allocatable, dimension(:,:) :: projct_g2f
  !
  real(wp), public, save, allocatable, dimension(:,:,:) :: interp_s2f
  real(wp), public, save, allocatable, dimension(:,:,:) :: interp_s2g
  real(wp), public, save, allocatable, dimension(:,:,:) :: interp_s2q
  real(wp), public, save, allocatable, dimension(:,:,:) :: projct_q2s
  real(wp), public, save, allocatable, dimension(:,:,:,:) :: projct_s2p
  !
  real(wp), public, save, allocatable, dimension(:,:,:,:) :: projct_limiter
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
subroutine write_variable(iw,var)
  !
  integer, intent(in) :: iw
  real(wp), dimension(:,:), intent(in) :: var
  !
  integer :: i,i1,i2,j,m
  !
continue
  !
  do i = 1,ncell
    !
    i1 = pts_in_cell_ptr(i)+1
    i2 = pts_in_cell_ptr(i+1)
    !
    write (iw,"(29x,i3)") i
    !
    do j = i1,i2
      write (iw,1) j,(var(m,j),m=1,size(var,1))
    end do
    !
  end do
  !
 !1 format (i6,3x,100(es14.6))
  1 format (i6,3x,*(es14.6))
  !
end subroutine write_variable
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
  if (allocated(gauss_1d)) then
    array%siz = size(gauss_1d)
    array%mem = array%siz*storage_size(gauss_1d)
    total%mem = total%mem + array%mem
    write (iunit,2) "'gauss_1d'",array%get_units()
  end if
  if (allocated(gauss_wts_1d)) then
    array%siz = size(gauss_wts_1d)
    array%mem = array%siz*storage_size(gauss_wts_1d)
    total%mem = total%mem + array%mem
    write (iunit,2) "'gauss_wts_1d'",array%get_units()
  end if
  if (allocated(gauss_wts_2d)) then
    array%siz = size(gauss_wts_2d)
    array%mem = array%siz*storage_size(gauss_wts_2d)
    total%mem = total%mem + array%mem
    write (iunit,2) "'gauss_wts_2d'",array%get_units()
  end if
  if (allocated(lobatto_1d)) then
    array%siz = size(lobatto_1d)
    array%mem = array%siz*storage_size(lobatto_1d)
    total%mem = total%mem + array%mem
    write (iunit,2) "'lobatto_1d'",array%get_units()
  end if
  if (allocated(lobatto_wts_1d)) then
    array%siz = size(lobatto_wts_1d)
    array%mem = array%siz*storage_size(lobatto_wts_1d)
    total%mem = total%mem + array%mem
    write (iunit,2) "'lobatto_wts_1d'",array%get_units()
  end if
  if (allocated(lobatto_wts_2d)) then
    array%siz = size(lobatto_wts_2d)
    array%mem = array%siz*storage_size(lobatto_wts_2d)
    total%mem = total%mem + array%mem
    write (iunit,2) "'lobatto_wts_2d'",array%get_units()
  end if
  if (allocated(solpts_1d)) then
    array%siz = size(solpts_1d)
    array%mem = array%siz*storage_size(solpts_1d)
    total%mem = total%mem + array%mem
    write (iunit,2) "'solpts_1d'",array%get_units()
  end if
  if (allocated(flxpts_1d)) then
    array%siz = size(flxpts_1d)
    array%mem = array%siz*storage_size(flxpts_1d)
    total%mem = total%mem + array%mem
    write (iunit,2) "'flxpts_1d'",array%get_units()
  end if
  if (allocated(solwts_1d)) then
    array%siz = size(solwts_1d)
    array%mem = array%siz*storage_size(solwts_1d)
    total%mem = total%mem + array%mem
    write (iunit,2) "'solwts_1d'",array%get_units()
  end if
  if (allocated(flxwts_1d)) then
    array%siz = size(flxwts_1d)
    array%mem = array%siz*storage_size(flxwts_1d)
    total%mem = total%mem + array%mem
    write (iunit,2) "'flxwts_1d'",array%get_units()
  end if
  if (allocated(rs)) then
    array%siz = size(rs)
    array%mem = array%siz*storage_size(rs)
    total%mem = total%mem + array%mem
    write (iunit,2) "'rs'",array%get_units()
  end if
  if (allocated(qdtpts_1d)) then
    array%siz = size(qdtpts_1d)
    array%mem = array%siz*storage_size(qdtpts_1d)
    total%mem = total%mem + array%mem
    write (iunit,2) "'qdtpts_1d'",array%get_units()
  end if
  if (allocated(qdtwts_1d)) then
    array%siz = size(qdtwts_1d)
    array%mem = array%siz*storage_size(qdtwts_1d)
    total%mem = total%mem + array%mem
    write (iunit,2) "'qdtwts_1d'",array%get_units()
  end if
  if (allocated(gauss2numer)) then
    array%siz = size(gauss2numer)
    array%mem = array%siz*storage_size(gauss2numer)
    total%mem = total%mem + array%mem
    write (iunit,2) "'gauss2numer'",array%get_units()
  end if
  if (allocated(lobatto2numer)) then
    array%siz = size(lobatto2numer)
    array%mem = array%siz*storage_size(lobatto2numer)
    total%mem = total%mem + array%mem
    write (iunit,2) "'lobatto2numer'",array%get_units()
  end if
  if (allocated(numer2gauss)) then
    array%siz = size(numer2gauss)
    array%mem = array%siz*storage_size(numer2gauss)
    total%mem = total%mem + array%mem
    write (iunit,2) "'numer2gauss'",array%get_units()
  end if
  if (allocated(numer2lobatto)) then
    array%siz = size(numer2lobatto)
    array%mem = array%siz*storage_size(numer2lobatto)
    total%mem = total%mem + array%mem
    write (iunit,2) "'numer2lobatto'",array%get_units()
  end if
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
  if (allocated(xyz_flxpts)) then
    array%siz = size(xyz_flxpts)
    array%mem = array%siz*storage_size(xyz_flxpts)
    total%mem = total%mem + array%mem
    write (iunit,2) "'xyz_flxpts'",array%get_units()
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
  if (allocated(faces_of_cell)) then
    array%siz = size(faces_of_cell)
    array%mem = array%siz*storage_size(faces_of_cell)
    total%mem = total%mem + array%mem
    write (iunit,2) "'faces_of_cell'",array%get_units()
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
  if (allocated(interp_f2g)) then
    array%siz = size(interp_f2g)
    array%mem = array%siz*storage_size(interp_f2g)
    total%mem = total%mem + array%mem
    write (iunit,2) "'interp_f2g'",array%get_units()
  end if
  if (allocated(projct_g2f)) then
    array%siz = size(projct_g2f)
    array%mem = array%siz*storage_size(projct_g2f)
    total%mem = total%mem + array%mem
    write (iunit,2) "'projct_g2f'",array%get_units()
  end if
  if (allocated(interp_s2f)) then
    array%siz = size(interp_s2f)
    array%mem = array%siz*storage_size(interp_s2f)
    total%mem = total%mem + array%mem
    write (iunit,2) "'interp_s2f'",array%get_units()
  end if
  if (allocated(interp_s2g)) then
    array%siz = size(interp_s2g)
    array%mem = array%siz*storage_size(interp_s2g)
    total%mem = total%mem + array%mem
    write (iunit,2) "'interp_s2g'",array%get_units()
  end if
  if (allocated(interp_s2q)) then
    array%siz = size(interp_s2q)
    array%mem = array%siz*storage_size(interp_s2q)
    total%mem = total%mem + array%mem
    write (iunit,2) "'interp_s2q'",array%get_units()
  end if
  if (allocated(projct_q2s)) then
    array%siz = size(projct_q2s)
    array%mem = array%siz*storage_size(projct_q2s)
    total%mem = total%mem + array%mem
    write (iunit,2) "'projct_q2s'",array%get_units()
  end if
  if (allocated(projct_s2p)) then
    array%siz = size(projct_s2p)
    array%mem = array%siz*storage_size(projct_s2p)
    total%mem = total%mem + array%mem
    write (iunit,2) "'projct_s2p'",array%get_units()
  end if
  if (allocated(projct_limiter)) then
    array%siz = size(projct_limiter)
    array%mem = array%siz*storage_size(projct_limiter)
    total%mem = total%mem + array%mem
    write (iunit,2) "'projct_limiter'",array%get_units()
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
module qpvar
  !
  ! Module containing all the temporary quad precision arrays
  ! used during the preprocessing/startup phase for computing
  ! the various arrays and matrices use throughout the solver
  !
  use module_kind_types, only : qp,alloc_error
  !
  implicit none
  !
  private
  !
  public :: deallocate_all_qp_arrays
  public :: qpvar_alloc_status
  !
  ! gauss_1d       : coordinates of 1D Gauss solution points of order n_order
  ! gauss_wts_1d   : weights     of 1D Gauss solution points of order n_order
  ! gauss_wts_2d   : weights     of 2D Gauss solution points of order n_order
  !
  real(qp), public, save, allocatable, dimension(:) :: gauss_1d_qp
  real(qp), public, save, allocatable, dimension(:) :: gauss_wts_1d_qp
  real(qp), public, save, allocatable, dimension(:) :: gauss_wts_2d_qp
  !
  ! lobatto_1d     : coordinates of 1D Lobatto solution points of order n_order
  ! lobatto_wts_1d : weights     of 1D Lobatto solution points of order n_order
  ! lobatto_wts_2d : weights     of 2D Lobatto solution points of order n_order
  !
  real(qp), public, save, allocatable, dimension(:) :: lobatto_1d_qp
  real(qp), public, save, allocatable, dimension(:) :: lobatto_wts_1d_qp
  real(qp), public, save, allocatable, dimension(:) :: lobatto_wts_2d_qp
  !
  ! solpts_1d : coordinates of 1D solution points of order n_order
  ! flxpts_1d : coordinates of 1D flux points of order n_order
  ! solwts_1d : weights of the 1D solution points of order n_order
  ! flxwts_1d : weights of the 1D flux points of order n_order
  ! rs        : 2D Lobatto nodes of triangle in r-s space
  !
  real(qp), public, save, allocatable, dimension(:)   :: vp_qp
  real(qp), public, save, allocatable, dimension(:)   :: solpts_1d_qp
  real(qp), public, save, allocatable, dimension(:)   :: flxpts_1d_qp
  real(qp), public, save, allocatable, dimension(:)   :: solwts_1d_qp
  real(qp), public, save, allocatable, dimension(:)   :: flxwts_1d_qp
  real(qp), public, save, allocatable, dimension(:,:) :: rs_qp
  !
  real(qp), public, save, allocatable, dimension(:)   :: qdtpts_1d_qp
  real(qp), public, save, allocatable, dimension(:)   :: qdtwts_1d_qp
  !
  ! gauss2numer   : matrix to move from Gauss   SP to the numerical SP
  ! lobatto2numer : matrix to move from Lobatto SP to the numerical SP
  !
  real(qp), public, save, allocatable, dimension(:,:) :: gauss2numer_qp
  real(qp), public, save, allocatable, dimension(:,:) :: lobatto2numer_qp
  !
  ! numer2gauss   : matrix to move from the numerical SP to Gauss   SP
  ! numer2lobatto : matrix to move from the numerical SP to Lobatto SP
  !
  real(qp), public, save, allocatable, dimension(:,:) :: numer2gauss_qp
  real(qp), public, save, allocatable, dimension(:,:) :: numer2lobatto_qp
  !
  ! Various projection/interpolation matrices
  !
  real(qp), public, save, allocatable, dimension(:,:) :: interp_f2g_qp
  real(qp), public, save, allocatable, dimension(:,:) :: projct_g2f_qp
  !
  real(qp), public, save, allocatable, dimension(:,:,:) :: interp_s2f_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: interp_s2g_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: interp_s2q_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: projct_q2s_qp
  real(qp), public, save, allocatable, dimension(:,:,:,:) :: projct_s2p_qp
  !
  real(qp), public, save, allocatable, dimension(:,:,:,:) :: projct_limiter_qp
  !
  ! Collocation derivative matrices
  !
  real(qp), public, save, allocatable, dimension(:,:,:) :: dr_mat_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: ds_mat_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: cormat_qp
  !
  ! Quadrature derivative matrices
  !
  real(qp), public, save, allocatable, dimension(:,:,:) :: dqr_LPmat_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: dqs_LPmat_qp
  !
contains
!
!###############################################################################
!
subroutine deallocate_all_qp_arrays()
  !
  !..Local Scalars ..
  integer :: ierr
  !
continue
  !
  if (allocated(gauss_1d_qp)) then
    deallocate ( gauss_1d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","gauss_1d_qp",2,ierr)
  end if
  !
  if (allocated(gauss_wts_1d_qp)) then
    deallocate ( gauss_wts_1d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","gauss_wts_1d_qp",2,ierr)
  end if
  !
  if (allocated(gauss_wts_2d_qp)) then
    deallocate ( gauss_wts_2d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","gauss_wts_2d_qp",2,ierr)
  end if
  !
  if (allocated(lobatto_1d_qp)) then
    deallocate ( lobatto_1d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","lobatto_1d_qp",2,ierr)
  end if
  !
  if (allocated(lobatto_wts_1d_qp)) then
    deallocate ( lobatto_wts_1d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","lobatto_wts_1d_qp",2,ierr)
  end if
  !
  if (allocated(lobatto_wts_2d_qp)) then
    deallocate ( lobatto_wts_2d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","lobatto_wts_2d_qp",2,ierr)
  end if
  !
  if (allocated(vp_qp)) then
    deallocate ( vp_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","vp_qp",2,ierr)
  end if
  !
  if (allocated(solpts_1d_qp)) then
    deallocate ( solpts_1d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","solpts_1d_qp",2,ierr)
  end if
  !
  if (allocated(flxpts_1d_qp)) then
    deallocate ( flxpts_1d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","flxpts_1d_qp",2,ierr)
  end if
  !
  if (allocated(solwts_1d_qp)) then
    deallocate ( solwts_1d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","solwts_1d_qp",2,ierr)
  end if
  !
  if (allocated(flxwts_1d_qp)) then
    deallocate ( flxwts_1d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","flxwts_1d_qp",2,ierr)
  end if
  !
  if (allocated(rs_qp)) then
    deallocate ( rs_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","rs_qp",2,ierr)
  end if
  !
  if (allocated(qdtpts_1d_qp)) then
    deallocate ( qdtpts_1d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","qdtpts_1d_qp",2,ierr)
  end if
  !
  if (allocated(qdtwts_1d_qp)) then
    deallocate ( qdtwts_1d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","qdtwts_1d_qp",2,ierr)
  end if
  !
  if (allocated(gauss2numer_qp)) then
    deallocate ( gauss2numer_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","gauss2numer_qp",2,ierr)
  end if
  !
  if (allocated(lobatto2numer_qp)) then
    deallocate ( lobatto2numer_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","lobatto2numer_qp",2,ierr)
  end if
  !
  if (allocated(numer2gauss_qp)) then
    deallocate ( numer2gauss_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","numer2gauss_qp",2,ierr)
  end if
  !
  if (allocated(numer2lobatto_qp)) then
    deallocate ( numer2lobatto_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","numer2lobatto_qp",2,ierr)
  end if
  !
  if (allocated(interp_f2g_qp)) then
    deallocate ( interp_f2g_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","interp_f2g_qp",2,ierr)
  end if
  !
  if (allocated(projct_g2f_qp)) then
    deallocate ( projct_g2f_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","projct_g2f_qp",2,ierr)
  end if
  !
  if (allocated(interp_s2f_qp)) then
    deallocate ( interp_s2f_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","interp_s2f_qp",2,ierr)
  end if
  !
  if (allocated(interp_s2g_qp)) then
    deallocate ( interp_s2g_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","interp_s2g_qp",2,ierr)
  end if
  !
  if (allocated(interp_s2q_qp)) then
    deallocate ( interp_s2q_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","interp_s2q_qp",2,ierr)
  end if
  !
  if (allocated(projct_q2s_qp)) then
    deallocate ( projct_q2s_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","projct_q2s_qp",2,ierr)
  end if
  !
  if (allocated(projct_s2p_qp)) then
    deallocate ( projct_s2p_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","projct_s2p_qp",2,ierr)
  end if
  !
  if (allocated(projct_limiter_qp)) then
    deallocate ( projct_limiter_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","projct_limiter_qp",2,ierr)
  end if
  !
  if (allocated(dr_mat_qp)) then
    deallocate ( dr_mat_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","dr_mat_qp",2,ierr)
  end if
  !
  if (allocated(ds_mat_qp)) then
    deallocate ( ds_mat_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","ds_mat_qp",2,ierr)
  end if
  !
  if (allocated(cormat_qp)) then
    deallocate ( cormat_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","cormat_qp",2,ierr)
  end if
  !
  if (allocated(dqr_LPmat_qp)) then
    deallocate ( dqr_LPmat_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","dqr_LPmat_qp",2,ierr)
  end if
  !
  if (allocated(dqs_LPmat_qp)) then
    deallocate ( dqs_LPmat_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","dqs_LPmat_qp",2,ierr)
  end if
  !
end subroutine deallocate_all_qp_arrays
!
!###############################################################################
!
subroutine qpvar_alloc_status(iunit)
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
  if (allocated(gauss_1d_qp)) then
    array%siz = size(gauss_1d_qp)
    array%mem = array%siz*storage_size(gauss_1d_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'gauss_1d_qp'",array%get_units()
  end if
  if (allocated(gauss_wts_1d_qp)) then
    array%siz = size(gauss_wts_1d_qp)
    array%mem = array%siz*storage_size(gauss_wts_1d_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'gauss_wts_1d_qp'",array%get_units()
  end if
  if (allocated(gauss_wts_2d_qp)) then
    array%siz = size(gauss_wts_2d_qp)
    array%mem = array%siz*storage_size(gauss_wts_2d_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'gauss_wts_2d_qp'",array%get_units()
  end if
  if (allocated(lobatto_1d_qp)) then
    array%siz = size(lobatto_1d_qp)
    array%mem = array%siz*storage_size(lobatto_1d_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'lobatto_1d_qp'",array%get_units()
  end if
  if (allocated(lobatto_wts_1d_qp)) then
    array%siz = size(lobatto_wts_1d_qp)
    array%mem = array%siz*storage_size(lobatto_wts_1d_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'lobatto_wts_1d_qp'",array%get_units()
  end if
  if (allocated(lobatto_wts_2d_qp)) then
    array%siz = size(lobatto_wts_2d_qp)
    array%mem = array%siz*storage_size(lobatto_wts_2d_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'lobatto_wts_2d_qp'",array%get_units()
  end if
  if (allocated(vp_qp)) then
    array%siz = size(vp_qp)
    array%mem = array%siz*storage_size(vp_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'vp_qp'",array%get_units()
  end if
  if (allocated(solpts_1d_qp)) then
    array%siz = size(solpts_1d_qp)
    array%mem = array%siz*storage_size(solpts_1d_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'solpts_1d_qp'",array%get_units()
  end if
  if (allocated(flxpts_1d_qp)) then
    array%siz = size(flxpts_1d_qp)
    array%mem = array%siz*storage_size(flxpts_1d_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'flxpts_1d_qp'",array%get_units()
  end if
  if (allocated(solwts_1d_qp)) then
    array%siz = size(solwts_1d_qp)
    array%mem = array%siz*storage_size(solwts_1d_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'solwts_1d_qp'",array%get_units()
  end if
  if (allocated(flxwts_1d_qp)) then
    array%siz = size(flxwts_1d_qp)
    array%mem = array%siz*storage_size(flxwts_1d_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'flxwts_1d_qp'",array%get_units()
  end if
  if (allocated(rs_qp)) then
    array%siz = size(rs_qp)
    array%mem = array%siz*storage_size(rs_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'rs_qp'",array%get_units()
  end if
  if (allocated(qdtpts_1d_qp)) then
    array%siz = size(qdtpts_1d_qp)
    array%mem = array%siz*storage_size(qdtpts_1d_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'qdtpts_1d_qp'",array%get_units()
  end if
  if (allocated(qdtwts_1d_qp)) then
    array%siz = size(qdtwts_1d_qp)
    array%mem = array%siz*storage_size(qdtwts_1d_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'qdtwts_1d_qp'",array%get_units()
  end if
  if (allocated(gauss2numer_qp)) then
    array%siz = size(gauss2numer_qp)
    array%mem = array%siz*storage_size(gauss2numer_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'gauss2numer_qp'",array%get_units()
  end if
  if (allocated(lobatto2numer_qp)) then
    array%siz = size(lobatto2numer_qp)
    array%mem = array%siz*storage_size(lobatto2numer_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'lobatto2numer_qp'",array%get_units()
  end if
  if (allocated(numer2gauss_qp)) then
    array%siz = size(numer2gauss_qp)
    array%mem = array%siz*storage_size(numer2gauss_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'numer2gauss_qp'",array%get_units()
  end if
  if (allocated(numer2lobatto_qp)) then
    array%siz = size(numer2lobatto_qp)
    array%mem = array%siz*storage_size(numer2lobatto_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'numer2lobatto_qp'",array%get_units()
  end if
  if (allocated(interp_f2g_qp)) then
    array%siz = size(interp_f2g_qp)
    array%mem = array%siz*storage_size(interp_f2g_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'interp_f2g_qp'",array%get_units()
  end if
  if (allocated(projct_g2f_qp)) then
    array%siz = size(projct_g2f_qp)
    array%mem = array%siz*storage_size(projct_g2f_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'projct_g2f_qp'",array%get_units()
  end if
  if (allocated(interp_s2f_qp)) then
    array%siz = size(interp_s2f_qp)
    array%mem = array%siz*storage_size(interp_s2f_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'interp_s2f_qp'",array%get_units()
  end if
  if (allocated(interp_s2g_qp)) then
    array%siz = size(interp_s2g_qp)
    array%mem = array%siz*storage_size(interp_s2g_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'interp_s2g_qp'",array%get_units()
  end if
  if (allocated(interp_s2q_qp)) then
    array%siz = size(interp_s2q_qp)
    array%mem = array%siz*storage_size(interp_s2q_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'interp_s2q_qp'",array%get_units()
  end if
  if (allocated(projct_q2s_qp)) then
    array%siz = size(projct_q2s_qp)
    array%mem = array%siz*storage_size(projct_q2s_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'projct_q2s_qp'",array%get_units()
  end if
  if (allocated(projct_s2p_qp)) then
    array%siz = size(projct_s2p_qp)
    array%mem = array%siz*storage_size(projct_s2p_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'projct_s2p_qp'",array%get_units()
  end if
  if (allocated(projct_limiter_qp)) then
    array%siz = size(projct_limiter_qp)
    array%mem = array%siz*storage_size(projct_limiter_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'projct_limiter_qp'",array%get_units()
  end if
  if (allocated(dr_mat_qp)) then
    array%siz = size(dr_mat_qp)
    array%mem = array%siz*storage_size(dr_mat_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'dr_mat_qp'",array%get_units()
  end if
  if (allocated(ds_mat_qp)) then
    array%siz = size(ds_mat_qp)
    array%mem = array%siz*storage_size(ds_mat_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'ds_mat_qp'",array%get_units()
  end if
  if (allocated(cormat_qp)) then
    array%siz = size(cormat_qp)
    array%mem = array%siz*storage_size(cormat_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'cormat_qp'",array%get_units()
  end if
  if (allocated(dqr_LPmat_qp)) then
    array%siz = size(dqr_LPmat_qp)
    array%mem = array%siz*storage_size(dqr_LPmat_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'dqr_LPmat_qp'",array%get_units()
  end if
  if (allocated(dqs_LPmat_qp)) then
    array%siz = size(dqs_LPmat_qp)
    array%mem = array%siz*storage_size(dqs_LPmat_qp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'dqs_LPmat_qp'",array%get_units()
  end if
  !
  total = total%get_units()
  write (iunit,3) total%mem
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module qpvar")
  2 format (" Array: ",a19,dt)
  3 format (" Total memory = ",i0,1x,a)
  !
end subroutine qpvar_alloc_status
!
!###############################################################################
!
end module qpvar
!
!###############################################################################
!
