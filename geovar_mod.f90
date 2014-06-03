module geovar
  !
  use module_kind_types, only : wp
  !
  implicit none
  !
  private
  !
  public :: write_variable
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
  ! n_order  : order of the solution polynomial requested in input file
  ! n_sp1dl  : number of solution points on a 1D edge/face
  ! n_sp2dt  : number of solution points within a 2D triangular cell
  ! n_sp2dq  : number of solution points within a 2d quadrilateral cell
  ! n_solpts : total number of solution points across all cells
  ! n_flxpts : total number of flux points across all cells
  !
  integer, public, save :: o_order
  integer, public, save :: p_order
  integer, public, save :: n_order
  integer, public, save :: n_sp1de
  integer, public, save :: n_fp1de
  integer, public, save :: n_sp2dt
  integer, public, save :: n_sp2dq
  integer, public, save :: n_solpts
  integer, public, save :: n_flxpts
  !
  ! over-integration parameters
  !
  logical, public, save :: using_quadrature
  integer, public, save :: q_order
  integer, public, save :: n_qp1de
  integer, public, save :: n_qp2dt
  integer, public, save :: n_qp2dq
  !
  ! neqns : number of governing equations
  !
  integer, public, parameter :: neqns = 4
  !
  ! nbfai : number of storage locations needed for each boundary face in array bface
  !
  integer, public, parameter :: nbfai = 10
  !
  ! nr    : number of dimensions for computation
  ! nnode : total number of grid nodes
  ! ncell : total number of grid cells
  ! nface : total number of edges/faces
  ! nfbnd : total number of boundary edges/faces
  ! nptpf : number of solution points per face
  !
  integer, public, save :: nr
  integer, public, save :: nnode
  integer, public, save :: ncell
  integer, public, save :: nface
  integer, public, save :: nfbnd
  integer, public, save :: nptpf
  integer, public, save :: nghost
  integer, public, save :: n_totcel
  integer, public, save :: n_totnod
  integer, public, save :: n_totpts
  !
  ! radau_1d       : coordinates of 1D Radau solution points of order n_order
  ! radau_wts_1d   : weights     of 1D Radau solution points of order n_order
  ! radau_wts_2d   : weights     of 2D Radau solution points of order n_order
  !
  real(wp), public, save, allocatable, dimension(:) :: radau_1d
  real(wp), public, save, allocatable, dimension(:) :: radau_wts_1d
  real(wp), public, save, allocatable, dimension(:) :: radau_wts_2d
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
  ! weights2d : weights of the solution points of order n_order
  ! rs        : 2D Lobatto nodes of triangle in r-s space
  ! 
  real(wp), public, save, allocatable, dimension(:)   :: solpts_1d
  real(wp), public, save, allocatable, dimension(:)   :: flxpts_1d
  real(wp), public, save, allocatable, dimension(:)   :: solwts_1d
  real(wp), public, save, allocatable, dimension(:)   :: flxwts_1d
  real(wp), public, save, allocatable, dimension(:,:) :: weights2d
  real(wp), public, save, allocatable, dimension(:,:) :: rs
  !
  real(wp), public, save, allocatable, dimension(:)   :: qdtpts_1d
  real(wp), public, save, allocatable, dimension(:)   :: qdtwts_1d
  real(wp), public, save, allocatable, dimension(:,:) :: qdtwts_2d
  !
  ! radau2numer   : matrix to move from Radau   SP to the numerical SP
  ! gauss2numer   : matrix to move from Gauss   SP to the numerical SP
  ! lobatto2numer : matrix to move from Lobatto SP to the numerical SP
  !
  real(wp), public, save, allocatable, dimension(:,:) :: radau2numer
  real(wp), public, save, allocatable, dimension(:,:) :: gauss2numer
  real(wp), public, save, allocatable, dimension(:,:) :: lobatto2numer
  !
  ! numer2radau   : matrix to move from the numerical SP to Radau   SP
  ! numer2gauss   : matrix to move from the numerical SP to Gauss   SP
  ! numer2lobatto : matrix to move from the numerical SP to Lobatto SP
  !
  real(wp), public, save, allocatable, dimension(:,:) :: numer2radau
  real(wp), public, save, allocatable, dimension(:,:) :: numer2gauss
  real(wp), public, save, allocatable, dimension(:,:) :: numer2lobatto
  !
  ! exact2numer : matrix to move from the exact SP to the numerical SP
  ! numer2exact : matrix to move from the numerical SP to the exact SP
  !
  real(wp), public, save, allocatable, dimension(:,:) :: exact2numer
  real(wp), public, save, allocatable, dimension(:,:) :: numer2exact
  ! 
  ! For cell i_cell
  ! pts_in_cell_ptr(i_cell)+1 gives the first solution point within i_cell
  ! pts_in_cell_ptr(i_cell+1) gives the last  solution point within i_cell
  !
  integer, public, save, allocatable, dimension(:) :: pts_in_cell_ptr
  !
  ! xyz_nodes : x,y coordinates of the vertices
  ! xyz       : x,y coordinates of all solution points
  ! xyz_error : x,y coordinates of all error quadrature points
  !
  real(wp), public, save, allocatable, dimension(:,:) :: xyz_nodes
  real(wp), public, save, allocatable, dimension(:,:) :: xyz
  real(wp), public, save, allocatable, dimension(:,:) :: xyz_error
  !
  ! xyz_radau   : x,y coordinates of all the Radau   SP
  ! xyz_gauss   : x,y coordinates of all the Gauss   SP
  ! xyz_lobatto : x,y coordinates of all the Lobatto SP
  !
  real(wp), public, save, allocatable, dimension(:,:) :: xyz_radau
  real(wp), public, save, allocatable, dimension(:,:) :: xyz_gauss
  real(wp), public, save, allocatable, dimension(:,:) :: xyz_lobatto
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
  integer,   public, save, allocatable, dimension(:)   :: global_nofc_ptr
  integer,   public, save, allocatable, dimension(:)   :: global_nofc
  integer,   public, save, allocatable, dimension(:)   :: global_pinc_ptr
  integer,   public, save, allocatable, dimension(:)   :: global_cell_offset
  !
  ! bface : Boundary face information
  !         bface(1,:) = type of boundary condition
  !         bface(2,:) = host cell of boundary face
  !         bface(3,:) = type of cell on boundary (1-edge,2-quad,3-triangle,etc.)
  !         bface(4,:) = face of host cell on boundary
  !                      for triangles: 1 -> nodes 1 and 2 of cell are on the boundary
  !                      for triangles: 2 -> nodes 2 and 3 of cell are on the boundary
  !                      for triangles: 3 -> nodes 3 and 1 of cell are on the boundary
  !                      for quads: 1 -> nodes 1 and 2 of cell are on the boundary
  !                      for quads: 2 -> nodes 2 and 3 of cell are on the boundary
  !                      for quads: 3 -> nodes 3 and 4 of cell are on the boundary
  !                      for quads: 4 -> nodes 4 and 1 of cell are on the boundary
  !         bface(5,:) = for a periodic or partition boundary face i_face with
  !                      host cell i_cell, this gives the boundary face j_face
  !                      that links to i_face. To find host cell j_cell for
  !                      periodic face j_face, use bface(2,bface(5,i_face)).
  !         bface(6,:) = for a periodic or partition boundary face i_face
  !                      with host cell i_cell, this gives the host cell
  !                      j_cell (i.e., bface(2,bface(5,i_face)) if they are
  !                      located on the same partition). This is needed if
  !                      they arent on the same partition because the local
  !                      process wont have access to the host cell information
  !                      in bface(2,j_face) for the adjacent partition
  !                      which is actually stored on another processor.
  !         bface(7,:) = for a periodic or partition boundary face i_face
  !                      with host cell i_cell, this gives the partition that
  !                      j_cell (i.e., bface(2,bface(5,i_face)) if they are
  !                      on the same partition) belongs to.
  !
  integer, public, save, allocatable, dimension(:,:) :: bface
  !
  ! For face i_face,
  ! pts_on_face(1,1:n_sp1de,i_face) : solution points of the left  cell on i_face
  ! pts_on_face(2,1:n_sp1de,i_face) : solution points of the right cell on i_face
  !
  integer, public, save, allocatable, dimension(:,:,:) :: pts_on_face
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
  ! For boundary face i_face, ...
  ! pts_of_periodic(ibeg:iend) gives the matching periodic solution points for
  !                            the solution points within the ghost cell on
  !                            the periodic boundary face i_face
  ! where : ibeg = pts_of_periodic_ptr(i_face)+1
  !         iend = pts_of_periodic_ptr(i_face+1)
  !
  integer, public, save, allocatable, dimension(:) :: pts_of_periodic_ptr
  integer, public, save, allocatable, dimension(:) :: pts_of_periodic
  !
  integer, public, save, allocatable, dimension(:,:) :: tria_connectivity
  integer, public, save, allocatable, dimension(:,:) :: quad_connectivity
  !
  ! Exponential filter matrices
  !
  real(wp), public, save, allocatable, dimension(:,:,:) :: filter
  !
  ! Various Vandermonde matrices
  !
  real(wp), public, save, allocatable, dimension(:,:) :: vand_1d
  real(wp), public, save, allocatable, dimension(:,:) :: vand_1d_inv
  real(wp), public, save, allocatable, dimension(:,:) :: vand_tria
  real(wp), public, save, allocatable, dimension(:,:) :: vand_tria_inv
  real(wp), public, save, allocatable, dimension(:,:) :: vand_quad
  real(wp), public, save, allocatable, dimension(:,:) :: vand_quad_inv
  !
  real(wp), public, save, allocatable, dimension(:,:,:) :: invM_mat
  !
  real(wp), public, save, allocatable, dimension(:,:,:) :: vand_qdtr
  real(wp), public, save, allocatable, dimension(:,:,:) :: vand_qdtr_inv
  !
  real(wp), public, save, allocatable, dimension(:,:) :: interp_f2g
  real(wp), public, save, allocatable, dimension(:,:) :: projct_g2f
  !
  real(wp), public, save, allocatable, dimension(:,:,:) :: interp_s2f
  real(wp), public, save, allocatable, dimension(:,:,:) :: interp_s2g
  real(wp), public, save, allocatable, dimension(:,:,:) :: interp_s2q
  real(wp), public, save, allocatable, dimension(:,:,:,:) :: interp_s2o
  real(wp), public, save, allocatable, dimension(:,:,:) :: projct_q2s
  real(wp), public, save, allocatable, dimension(:,:,:,:) :: projct_s2p
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
  public :: deallocate_all_qp_arrays
  !
  ! radau_1d       : coordinates of 1D Radau solution points of order n_order
  ! radau_wts_1d   : weights     of 1D Radau solution points of order n_order
  ! radau_wts_2d   : weights     of 2D Radau solution points of order n_order
  !
  real(qp), public, save, allocatable, dimension(:) :: radau_1d_qp
  real(qp), public, save, allocatable, dimension(:) :: radau_wts_1d_qp
  real(qp), public, save, allocatable, dimension(:) :: radau_wts_2d_qp
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
  ! weights2d : weights of the solution points of order n_order
  ! rs        : 2D Lobatto nodes of triangle in r-s space
  !
  real(qp), public, save, allocatable, dimension(:)   :: vp_qp
  real(qp), public, save, allocatable, dimension(:)   :: solpts_1d_qp
  real(qp), public, save, allocatable, dimension(:)   :: flxpts_1d_qp
  real(qp), public, save, allocatable, dimension(:)   :: solwts_1d_qp
  real(qp), public, save, allocatable, dimension(:)   :: flxwts_1d_qp
  real(qp), public, save, allocatable, dimension(:,:) :: weights2d_qp
  real(qp), public, save, allocatable, dimension(:,:) :: rs_qp
  !
  real(qp), public, save, allocatable, dimension(:)   :: qdtpts_1d_qp
  real(qp), public, save, allocatable, dimension(:)   :: qdtwts_1d_qp
  real(qp), public, save, allocatable, dimension(:,:) :: qdtwts_2d_qp
  !
  ! radau2numer   : matrix to move from Radau   SP to the numerical SP
  ! gauss2numer   : matrix to move from Gauss   SP to the numerical SP
  ! lobatto2numer : matrix to move from Lobatto SP to the numerical SP
  !
  real(qp), public, save, allocatable, dimension(:,:) :: radau2numer_qp
  real(qp), public, save, allocatable, dimension(:,:) :: gauss2numer_qp
  real(qp), public, save, allocatable, dimension(:,:) :: lobatto2numer_qp
  !
  ! numer2radau   : matrix to move from the numerical SP to Radau   SP
  ! numer2gauss   : matrix to move from the numerical SP to Gauss   SP
  ! numer2lobatto : matrix to move from the numerical SP to Lobatto SP
  !
  real(qp), public, save, allocatable, dimension(:,:) :: numer2radau_qp
  real(qp), public, save, allocatable, dimension(:,:) :: numer2gauss_qp
  real(qp), public, save, allocatable, dimension(:,:) :: numer2lobatto_qp
  !
  ! exact2numer : matrix to move from the exact SP to the numerical SP
  ! numer2exact : matrix to move from the numerical SP to the exact SP
  !
  real(qp), public, save, allocatable, dimension(:,:) :: exact2numer_qp
  real(qp), public, save, allocatable, dimension(:,:) :: numer2exact_qp
  ! 
  ! Exponential filter matrices
  !
  real(qp), public, save, allocatable, dimension(:,:,:) :: filter_qp
  !
  ! Various Vandermonde matrices
  !
  real(qp), public, save, allocatable, dimension(:,:) :: vand_1d_qp
  real(qp), public, save, allocatable, dimension(:,:) :: vand_1d_inv_qp
  real(qp), public, save, allocatable, dimension(:,:) :: vand_tria_qp
  real(qp), public, save, allocatable, dimension(:,:) :: vand_tria_inv_qp
  real(qp), public, save, allocatable, dimension(:,:) :: vand_quad_qp
  real(qp), public, save, allocatable, dimension(:,:) :: vand_quad_inv_qp
  !
  ! Inverse mass matrix
  !
  real(qp), public, save, allocatable, dimension(:,:,:) :: invM_mat_qp
  !
  ! Quadrature Vandermonde matrices
  !
  real(qp), public, save, allocatable, dimension(:,:,:) :: vand_qdtr_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: vand_qdtr_inv_qp
  !
  ! Various projection/interpolation matrices
  !
  real(qp), public, save, allocatable, dimension(:,:) :: interp_f2g_qp
  real(qp), public, save, allocatable, dimension(:,:) :: projct_g2f_qp
  !
  real(qp), public, save, allocatable, dimension(:,:,:) :: interp_s2f_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: interp_s2g_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: interp_s2q_qp
  real(qp), public, save, allocatable, dimension(:,:,:,:) :: interp_s2o_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: projct_q2s_qp
  real(qp), public, save, allocatable, dimension(:,:,:,:) :: projct_s2p_qp
  !
  ! Collocation derivative matrices
  !
  real(qp), public, save, allocatable, dimension(:,:)   :: d1_mat_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: dr_mat_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: ds_mat_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: dr_hat_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: ds_hat_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: cormat_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: cormat_dg_qp
  !
  ! Quadrature derivative matrices
  !
  real(qp), public, save, allocatable, dimension(:,:)   :: dq1D_mat_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: dqr_mat_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: dqs_mat_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: dqr_LPmat_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: dqs_LPmat_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: dqr_CRmat_qp
  real(qp), public, save, allocatable, dimension(:,:,:) :: dqs_CRmat_qp
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
  if (allocated(radau_1d_qp)) then
    deallocate ( radau_1d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","radau_1d_qp",2,ierr)
  end if
  !
  if (allocated(radau_wts_1d_qp)) then
    deallocate ( radau_wts_1d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","radau_wts_1d_qp",2,ierr)
  end if
  !
  if (allocated(radau_wts_2d_qp)) then
    deallocate ( radau_wts_2d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","radau_wts_2d_qp",2,ierr)
  end if
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
  if (allocated(weights2d_qp)) then
    deallocate ( weights2d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","weights2d_qp",2,ierr)
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
  if (allocated(qdtwts_2d_qp)) then
    deallocate ( qdtwts_2d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","qdtwts_2d_qp",2,ierr)
  end if
  !
  if (allocated(radau2numer_qp)) then
    deallocate ( radau2numer_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","radau2numer_qp",2,ierr)
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
  if (allocated(numer2radau_qp)) then
    deallocate ( numer2radau_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","numer2radau_qp",2,ierr)
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
  if (allocated(exact2numer_qp)) then
    deallocate ( exact2numer_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","exact2numer_qp",2,ierr)
  end if
  !
  if (allocated(numer2exact_qp)) then
    deallocate ( numer2exact_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","numer2exact_qp",2,ierr)
  end if
  !
  if (allocated(filter_qp)) then
    deallocate ( filter_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","filter_qp",2,ierr)
  end if
  !
  if (allocated(vand_1d_qp)) then
    deallocate ( vand_1d_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","vand_1d_qp",2,ierr)
  end if
  !
  if (allocated(vand_1d_inv_qp)) then
    deallocate ( vand_1d_inv_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","vand_1d_inv_qp",2,ierr)
  end if
  !
  if (allocated(vand_tria_qp)) then
    deallocate ( vand_tria_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","vand_tria_qp",2,ierr)
  end if
  !
  if (allocated(vand_tria_inv_qp)) then
    deallocate ( vand_tria_inv_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","vand_tria_inv_qp",2,ierr)
  end if
  !
  if (allocated(vand_quad_qp)) then
    deallocate ( vand_quad_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","vand_quad_qp",2,ierr)
  end if
  !
  if (allocated(vand_quad_inv_qp)) then
    deallocate ( vand_quad_inv_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","vand_quad_inv_qp",2,ierr)
  end if
  !
  if (allocated(invM_mat_qp)) then
    deallocate ( invM_mat_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","invM_mat_qp",2,ierr)
  end if
  !
  if (allocated(vand_qdtr_qp)) then
    deallocate ( vand_qdtr_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","vand_qdtr_qp",2,ierr)
  end if
  !
  if (allocated(vand_qdtr_inv_qp)) then
    deallocate ( vand_qdtr_inv_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","vand_qdtr_inv_qp",2,ierr)
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
  if (allocated(interp_s2o_qp)) then
    deallocate ( interp_s2o_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","interp_s2o_qp",2,ierr)
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
  if (allocated(d1_mat_qp)) then
    deallocate ( d1_mat_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","d1_mat_qp",2,ierr)
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
  if (allocated(dr_hat_qp)) then
    deallocate ( dr_hat_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","dr_hat_qp",2,ierr)
  end if
  !
  if (allocated(ds_hat_qp)) then
    deallocate ( ds_hat_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","ds_hat_qp",2,ierr)
  end if
  !
  if (allocated(cormat_qp)) then
    deallocate ( cormat_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","cormat_qp",2,ierr)
  end if
  !
  if (allocated(cormat_dg_qp)) then
    deallocate ( cormat_dg_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","cormat_dg_qp",2,ierr)
  end if
  !
  if (allocated(dq1D_mat_qp)) then
    deallocate ( dq1D_mat_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","dq1D_mat_qp",2,ierr)
  end if
  !
  if (allocated(dqr_mat_qp)) then
    deallocate ( dqr_mat_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","dqr_mat_qp",2,ierr)
  end if
  !
  if (allocated(dqs_mat_qp)) then
    deallocate ( dqs_mat_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","dqs_mat_qp",2,ierr)
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
  if (allocated(dqr_CRmat_qp)) then
    deallocate ( dqr_CRmat_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","dqr_CRmat_qp",2,ierr)
  end if
  !
  if (allocated(dqs_CRmat_qp)) then
    deallocate ( dqs_CRmat_qp , stat=ierr )
    call alloc_error("deallocate_all_qp_arrays","dqs_CRmat_qp",2,ierr)
  end if
  !
  !
end subroutine deallocate_all_qp_arrays
!
!###############################################################################
!
end module qpvar
!
!###############################################################################
!
module geovar_test
  !
  use module_kind_types, only : wp, dflux_opt_ktmod, project_derivatives
  !
  use module_kind_types, only : p_order => p_order_ktmod
  use module_kind_types, only : n_order => n_order_ktmod
  use module_kind_types, only : q_order => q_order_ktmod
  use module_kind_types, only : using_quadrature => using_quadrature_ktmod
  !
  implicit none
  !
  private
  !
  public :: write_variable
  !
  public :: p_order
  public :: n_order
  public :: q_order
  public :: using_quadrature
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
  ! n_order  : order of the solution polynomial requested in input file
  ! n_sp1dl  : number of solution points on a 1D edge/face
  ! n_sp2dt  : number of solution points within a 2D triangular cell
  ! n_sp2dq  : number of solution points within a 2d quadrilateral cell
  ! n_solpts : total number of solution points across all cells
  ! n_flxpts : total number of flux points across all cells
  !
  ! Parameters controlled by pre-processor directives
  !
 !integer, public, parameter :: p_order = p_order_fpp
 !integer, public, parameter :: n_order = n_order_fpp
 !integer, public, parameter :: q_order = q_order_fpp
 !logical, public, parameter :: using_quadrature = using_quadrature_fpp
  !
  integer, public, parameter :: n_sp1de = n_order + 1
  integer, public, parameter :: n_sp2dt = n_sp1de*(n_sp1de+1)/2
  integer, public, parameter :: n_sp2dq = n_sp1de*n_sp1de
  !
  ! over-integration parameters
  !
  integer, public, parameter :: n_qp1de = q_order + 1
  integer, public, parameter :: n_qp2dt = n_qp1de*(n_qp1de+1)/2
  integer, public, parameter :: n_qp2dq = n_qp1de*n_qp1de
  !
  integer, public, save :: o_order
  integer, public, save :: n_fp1de
 !integer, public, parameter :: n_fp1de = merge( n_qp1de , n_sp1de , &
 !                using_quadrature .and. dflux_opt_ktmod == project_derivatives )
  !
  ! neqns : number of governing equations
  !
  integer, public, parameter :: neqns = 4
  !
  ! nbfai : number of storage locations needed for each boundary face in array bface
  !
  integer, public, parameter :: nbfai = 10
  !
  ! nr    : number of dimensions for computation
  ! nnode : total number of grid nodes
  ! ncell : total number of grid cells
  ! nface : total number of edges/faces
  ! nfbnd : total number of boundary edges/faces
  ! nptpf : number of solution points per face
  !
  integer, public, save :: n_solpts
  integer, public, save :: n_flxpts
  !
  integer, public, save :: nr
  integer, public, save :: nnode
  integer, public, save :: ncell
  integer, public, save :: nface
  integer, public, save :: nfbnd
  integer, public, save :: nptpf
  integer, public, save :: nghost
  integer, public, save :: n_totcel
  integer, public, save :: n_totnod
  integer, public, save :: n_totpts
  !
  ! radau_1d       : coordinates of 1D Radau solution points of order n_order
  ! radau_wts_1d   : weights     of 1D Radau solution points of order n_order
  ! radau_wts_2d   : weights     of 2D Radau solution points of order n_order
  !
  real(wp), public, save, allocatable, dimension(:) :: radau_1d
  real(wp), public, save, allocatable, dimension(:) :: radau_wts_1d
  real(wp), public, save, allocatable, dimension(:) :: radau_wts_2d
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
  ! weights2d : weights of the solution points of order n_order
  ! rs        : 2D Lobatto nodes of triangle in r-s space
  ! 
  real(wp), public, save, allocatable, dimension(:)   :: solpts_1d
  real(wp), public, save, allocatable, dimension(:)   :: flxpts_1d
  real(wp), public, save, allocatable, dimension(:)   :: solwts_1d
  real(wp), public, save, allocatable, dimension(:)   :: flxwts_1d
  real(wp), public, save, allocatable, dimension(:,:) :: weights2d
  real(wp), public, save, allocatable, dimension(:,:) :: rs
  !
  real(wp), public, save, allocatable, dimension(:)   :: qdtpts_1d
  real(wp), public, save, allocatable, dimension(:)   :: qdtwts_1d
  real(wp), public, save, allocatable, dimension(:,:) :: qdtwts_2d
  !
  ! radau2numer   : matrix to move from Radau   SP to the numerical SP
  ! gauss2numer   : matrix to move from Gauss   SP to the numerical SP
  ! lobatto2numer : matrix to move from Lobatto SP to the numerical SP
  !
  real(wp), public, save, allocatable, dimension(:,:) :: radau2numer
  real(wp), public, save, allocatable, dimension(:,:) :: gauss2numer
  real(wp), public, save, allocatable, dimension(:,:) :: lobatto2numer
  !
  ! numer2radau   : matrix to move from the numerical SP to Radau   SP
  ! numer2gauss   : matrix to move from the numerical SP to Gauss   SP
  ! numer2lobatto : matrix to move from the numerical SP to Lobatto SP
  !
  real(wp), public, save, allocatable, dimension(:,:) :: numer2radau
  real(wp), public, save, allocatable, dimension(:,:) :: numer2gauss
  real(wp), public, save, allocatable, dimension(:,:) :: numer2lobatto
  !
  ! exact2numer : matrix to move from the exact SP to the numerical SP
  ! numer2exact : matrix to move from the numerical SP to the exact SP
  !
  real(wp), public, save, allocatable, dimension(:,:) :: exact2numer
  real(wp), public, save, allocatable, dimension(:,:) :: numer2exact
  ! 
  ! For cell i_cell
  ! pts_in_cell_ptr(i_cell)+1 gives the first solution point within i_cell
  ! pts_in_cell_ptr(i_cell+1) gives the last  solution point within i_cell
  !
  integer, public, save, allocatable, dimension(:) :: pts_in_cell_ptr
  !
  ! xyz_nodes : x,y coordinates of the vertices
  ! xyz       : x,y coordinates of all solution points
  ! xyz_error : x,y coordinates of all error quadrature points
  !
  real(wp), public, save, allocatable, dimension(:,:) :: xyz_nodes
  real(wp), public, save, allocatable, dimension(:,:) :: xyz
  real(wp), public, save, allocatable, dimension(:,:) :: xyz_error
  !
  ! xyz_radau   : x,y coordinates of all the Radau   SP
  ! xyz_gauss   : x,y coordinates of all the Gauss   SP
  ! xyz_lobatto : x,y coordinates of all the Lobatto SP
  !
  real(wp), public, save, allocatable, dimension(:,:) :: xyz_radau
  real(wp), public, save, allocatable, dimension(:,:) :: xyz_gauss
  real(wp), public, save, allocatable, dimension(:,:) :: xyz_lobatto
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
  integer,   public, save, allocatable, dimension(:)   :: global_nofc_ptr
  integer,   public, save, allocatable, dimension(:)   :: global_nofc
  integer,   public, save, allocatable, dimension(:)   :: global_pinc_ptr
  integer,   public, save, allocatable, dimension(:)   :: global_cell_offset
  !
  ! bface : Boundary face information
  !         bface(1,:) = type of boundary condition
  !         bface(2,:) = host cell of boundary face
  !         bface(3,:) = type of cell on boundary (1-edge,2-quad,3-triangle,etc.)
  !         bface(4,:) = face of host cell on boundary
  !                      for triangles: 1 -> nodes 1 and 2 of cell are on the boundary
  !                      for triangles: 2 -> nodes 2 and 3 of cell are on the boundary
  !                      for triangles: 3 -> nodes 3 and 1 of cell are on the boundary
  !                      for quads: 1 -> nodes 1 and 2 of cell are on the boundary
  !                      for quads: 2 -> nodes 2 and 3 of cell are on the boundary
  !                      for quads: 3 -> nodes 3 and 4 of cell are on the boundary
  !                      for quads: 4 -> nodes 4 and 1 of cell are on the boundary
  !         bface(5,:) = for a periodic or partition boundary face i_face with
  !                      host cell i_cell, this gives the boundary face j_face
  !                      that links to i_face. To find host cell j_cell for
  !                      periodic face j_face, use bface(2,bface(5,i_face)).
  !         bface(6,:) = for a periodic or partition boundary face i_face
  !                      with host cell i_cell, this gives the host cell
  !                      j_cell (i.e., bface(2,bface(5,i_face)) if they are
  !                      located on the same partition). This is needed if
  !                      they arent on the same partition because the local
  !                      process wont have access to the host cell information
  !                      in bface(2,j_face) for the adjacent partition
  !                      which is actually stored on another processor.
  !         bface(7,:) = for a periodic or partition boundary face i_face
  !                      with host cell i_cell, this gives the partition that
  !                      j_cell (i.e., bface(2,bface(5,i_face)) if they are
  !                      on the same partition) belongs to.
  !
  integer, public, save, allocatable, dimension(:,:) :: bface
  !
  ! For face i_face,
  ! pts_on_face(1,1:n_sp1de,i_face) : solution points of the left  cell on i_face
  ! pts_on_face(2,1:n_sp1de,i_face) : solution points of the right cell on i_face
  !
  integer, public, save, allocatable, dimension(:,:,:) :: pts_on_face
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
  ! For boundary face i_face, ...
  ! pts_of_periodic(ibeg:iend) gives the matching periodic solution points for
  !                            the solution points within the ghost cell on
  !                            the periodic boundary face i_face
  ! where : ibeg = pts_of_periodic_ptr(i_face)+1
  !         iend = pts_of_periodic_ptr(i_face+1)
  !
  integer, public, save, allocatable, dimension(:) :: pts_of_periodic_ptr
  integer, public, save, allocatable, dimension(:) :: pts_of_periodic
  !
  integer, public, save, allocatable, dimension(:,:) :: tria_connectivity
  integer, public, save, allocatable, dimension(:,:) :: quad_connectivity
  !
  ! Exponential filter matrices
  !
  real(wp), public, save, allocatable, dimension(:,:,:) :: filter
  !
  ! Various Vandermonde matrices
  !
  real(wp), public, save, allocatable, dimension(:,:) :: vand_1d
  real(wp), public, save, allocatable, dimension(:,:) :: vand_1d_inv
  real(wp), public, save, allocatable, dimension(:,:) :: vand_tria
  real(wp), public, save, allocatable, dimension(:,:) :: vand_tria_inv
  real(wp), public, save, allocatable, dimension(:,:) :: vand_quad
  real(wp), public, save, allocatable, dimension(:,:) :: vand_quad_inv
  !
  real(wp), public, save, allocatable, dimension(:,:,:) :: invM_mat
  !
  real(wp), public, save, allocatable, dimension(:,:,:) :: vand_qdtr
  real(wp), public, save, allocatable, dimension(:,:,:) :: vand_qdtr_inv
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
end module geovar_test
!
!###############################################################################
!
