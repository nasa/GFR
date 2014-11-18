module metrics_mod
  !
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
  public :: get_solution_points_1d
  public :: get_solution_points_2d
  public :: get_grid_connectivity
  public :: get_metrics
  !
  integer, parameter :: mindex = 15
  integer, parameter :: nstack = 50
  integer, parameter :: jn     = 2000
  !
  integer, parameter :: DG_TriPoints = 1
  integer, parameter :: HT_TriPoints = 2
  !
  character(len=*), parameter :: CG_error = &
                   '(/,a," points at the Chebyshev-Gauss'// &
                         ' nodes is not currently supported!",/)'
  character(len=*), parameter :: CGL_error = &
                   '(/,a," points at the Chebyshev-Gauss-Lobatto'// &
                         ' nodes is not currently supported!",/)'
  character(len=*), parameter :: PNT_error = &
                   '(/,"Option for location of ",a,"'// &
                       ' points is not valid!",/)'
  !
  interface get_quad_weights_at_solpts
    module procedure get_quad_weights_at_solpts_DP, &
                     get_quad_weights_at_solpts_QP
  end interface get_quad_weights_at_solpts
  !
  interface get_triangle_weights_at_solpts
    module procedure get_triangle_weights_at_solpts_DP, &
                     get_triangle_weights_at_solpts_QP
  end interface get_triangle_weights_at_solpts
  !
contains
!
!###############################################################################
!
pure function get_solution_points_1d(xi,xnodes) result(return_value)
  !
  real(wp), dimension(:), intent(in) :: xi
  real(wp), dimension(:), intent(in) :: xnodes
  !
  real(wp), dimension(size(xi),size(xnodes)-1) :: return_value
  !
  integer :: np,nc,i,j
  !
continue
  !
  np = size(xi)
  nc = size(xnodes)-1
  !
  return_value(:,:) = zero
  !
  do j = 1,nc
    do i = 1,np
      return_value(i,j) = xnodes(j) + half*(xi(i)+one)*(xnodes(j+1)-xnodes(j))
    end do
  end do
  !
end function get_solution_points_1d
!
!###############################################################################
!
subroutine get_solution_points_2d()
  !
  use geovar
  use qpvar,          only : radau_1d_qp,radau_wts_1d_qp,radau_wts_2d_qp
  use qpvar,          only : gauss_1d_qp,gauss_wts_1d_qp,gauss_wts_2d_qp
  use qpvar,          only : lobatto_1d_qp,lobatto_wts_1d_qp,lobatto_wts_2d_qp
  use qpvar,          only : solpts_1d_qp,solwts_1d_qp,vp_qp,rs_qp
  use qpvar,          only : flxpts_1d_qp,flxwts_1d_qp,weights2d_qp
  use ovar,           only : iop_solpts,iop_flxpts,iop_tripts,error_loc
  use metvar,         only : d1_mat
  use polynomial_mod, only : nodes_legendre_gauss, &
                             nodes_legendre_gauss_lobatto, &
                             nodes_chebyshev_gauss, &
                             nodes_chebyshev_gauss_lobatto, &
                             nodes_right_radau, &
                             weights_legendre_gauss, &
                             weights_legendre_gauss_lobatto, &
                             weights_chebyshev_gauss, &
                             weights_chebyshev_gauss_lobatto, &
                             weights_right_radau, &
                             derivative_matrix_legendre_gauss, &
                             derivative_matrix_legendre_gauss_lobatto, &
                             derivative_matrix_chebyshev_gauss_lobatto
  use triangle_mod,   only : TriNodes2D_HT, &
                             TriNodes2D_DG
  use parallel_mod,   only : collect_global_solpts
  !
  !.. Local Scalars ..
  integer   :: i,j,n,n1,n2,p1,p2,np,ierr
  real(wp) :: r,s
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr,1:4) :: xy
  real(wp), dimension(1:n_sp1de) :: vp
  !
continue
  !
  ! We need to reset various grid size parameters in case partitioning
  ! was done for use with multiple processors. If only using with a
  ! single processor, this should just give the values we already have.
  !
  nfbnd = size(bface,dim=2)
  n_totcel = size(nodes_of_cell_ptr) - 1
  ncell = minval(bface(10,:)) - 1
  nghost = n_totcel - ncell
  n_totnod = size(xyz_nodes,dim=2)
  nnode = maxval(nodes_of_cell(1:nodes_of_cell_ptr(ncell+1)))
  !
  ! Get the 1D coordinates and weights of various node types
  !
  ! Allocate the working precision arrays
  !
  allocate ( radau_1d(1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","radau_1d",1,ierr)
  allocate ( radau_wts_1d(1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","radau_wts_1d",1,ierr)
  allocate ( radau_wts_2d(1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","radau_wts_2d",1,ierr)
  !
  allocate ( gauss_1d(1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","gauss_1d",1,ierr)
  allocate ( gauss_wts_1d(1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","gauss_wts_1d",1,ierr)
  allocate ( gauss_wts_2d(1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","gauss_wts_2d",1,ierr)
  !
  allocate ( lobatto_1d(1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","lobatto_1d",1,ierr)
  allocate ( lobatto_wts_1d(1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","lobatto_wts_1d",1,ierr)
  allocate ( lobatto_wts_2d(1:n_sp2dq) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","lobatto_wts_2d",1,ierr)
  !
  allocate ( solpts_1d(1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","solpts_1d",1,ierr)
  allocate ( solwts_1d(1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","solwts_1d",1,ierr)
  !
  allocate ( flxpts_1d(1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","flxpts_1d",1,ierr)
  allocate ( flxwts_1d(1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","flxwts_1d",1,ierr)
  !
  allocate ( rs(1:nr,1:n_sp2dt) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","rs",1,ierr)
  allocate ( weights2d(1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","weights2d",1,ierr)
  !
  ! #####################################################################
  ! #####################################################################
  ! ##     USE QUAD PRECISION ARRAYS TO CREATE ALL OF THE ARRAYS       ##
  ! ##   PREVIOUSLY ALLOCATED. WE WILL THEN PERFORM A CHOP TO ZERO     ##
  ! ##  EVERYTHING WITH AN ABSOLUTE VALUE SMALLER THAN EPSILON(0_WP).  ##
  ! #####################################################################
  ! #####################################################################
  !
  ! Allocate the quad precision arrays
  !
  allocate ( radau_1d_qp(1:n_sp1de) , source=qzero , stat=ierr )
  call alloc_error("get_solution_points_2d","radau_1d_qp",1,ierr)
  allocate ( radau_wts_1d_qp(1:n_sp1de) , source=qzero , stat=ierr )
  call alloc_error("get_solution_points_2d","radau_wts_1d_qp",1,ierr)
  allocate ( radau_wts_2d_qp(1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("get_solution_points_2d","radau_wts_2d_qp",1,ierr)
  !
  allocate ( gauss_1d_qp(1:n_sp1de) , source=qzero , stat=ierr )
  call alloc_error("get_solution_points_2d","gauss_1d_qp",1,ierr)
  allocate ( gauss_wts_1d_qp(1:n_sp1de) , source=qzero , stat=ierr )
  call alloc_error("get_solution_points_2d","gauss_wts_1d_qp",1,ierr)
  allocate ( gauss_wts_2d_qp(1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("get_solution_points_2d","gauss_wts_2d_qp",1,ierr)
  !
  allocate ( lobatto_1d_qp(1:n_sp1de) , source=qzero , stat=ierr )
  call alloc_error("get_solution_points_2d","lobatto_1d_qp",1,ierr)
  allocate ( lobatto_wts_1d_qp(1:n_sp1de) , source=qzero , stat=ierr )
  call alloc_error("get_solution_points_2d","lobatto_wts_1d_qp",1,ierr)
  allocate ( lobatto_wts_2d_qp(1:n_sp2dq) , source=qzero , stat=ierr )
  call alloc_error("get_solution_points_2d","lobatto_wts_2d_qp",1,ierr)
  !
  allocate ( solpts_1d_qp(1:n_sp1de) , source=qzero , stat=ierr )
  call alloc_error("get_solution_points_2d","solpts_1d_qp",1,ierr)
  allocate ( solwts_1d_qp(1:n_sp1de) , source=qzero , stat=ierr )
  call alloc_error("get_solution_points_2d","solwts_1d_qp",1,ierr)
  allocate ( vp_qp(1:n_sp1de) , source=qzero , stat=ierr )
  call alloc_error("get_solution_points_2d","vp_qp",1,ierr)
  !
  allocate ( flxpts_1d_qp(1:n_sp1de) , source=qzero , stat=ierr )
  call alloc_error("get_solution_points_2d","flxpts_1d_qp",1,ierr)
  allocate ( flxwts_1d_qp(1:n_sp1de) , source=qzero , stat=ierr )
  call alloc_error("get_solution_points_2d","flxwts_1d_qp",1,ierr)
  !
  allocate ( rs_qp(1:nr,1:n_sp2dt) , source=qzero , stat=ierr )
  call alloc_error("get_solution_points_2d","rs_qp",1,ierr)
  allocate ( weights2d_qp(1:n_sp2dq,1:2) , source=qzero , stat=ierr )
  call alloc_error("get_solution_points_2d","weights2d_qp",1,ierr)
  !
  ! Create the 1D arrays for the locations of the Gauss, Lobatto,
  ! and Radau solution points and their corresponding weights
  !
  call nodes_legendre_gauss(n_sp1de,gauss_1d_qp,vp_qp)
  call weights_legendre_gauss(n_sp1de,gauss_1d_qp,vp_qp,gauss_wts_1d_qp)
  !
  call nodes_legendre_gauss_lobatto(n_order,lobatto_1d_qp,vp_qp)
  call weights_legendre_gauss_lobatto(n_order,vp_qp,lobatto_wts_1d_qp)
  !
  call nodes_right_radau(n_order,radau_1d_qp)
  call weights_right_radau(n_order,radau_1d_qp,radau_wts_1d_qp)
  !
  ! Create the 2D weights for the Gauss and Lobatto solution points on quads
  !
  gauss_wts_2d_qp = get_quad_weights_at_solpts(gauss_wts_1d_qp)
  lobatto_wts_2d_qp = get_quad_weights_at_solpts(lobatto_wts_1d_qp)
  radau_wts_2d_qp = get_quad_weights_at_solpts(radau_wts_1d_qp)
  !
  ! Now get the 1D solution points for the computation
  ! NOTE: We will save into vp whichever of gauss_1d and lobatto_1d isnt solpts_1d
  !
  if (iop_solpts == Legendre_Gauss) then
    !
    solpts_1d_qp(:) = gauss_1d_qp(:)
    solwts_1d_qp(:) = gauss_wts_1d_qp(:)
    !
  else if (iop_solpts == Legendre_Gauss_Lobatto) then
    !
    solpts_1d_qp(:) = lobatto_1d_qp(:)
    solwts_1d_qp(:) = lobatto_wts_1d_qp(:)
    !
  else if (iop_solpts == Chebyshev_Gauss) then
    !
    call nodes_chebyshev_gauss(n_sp1de,solpts_1d_qp)
    call weights_chebyshev_gauss(n_sp1de,solwts_1d_qp)
    write (iout,CG_error) "Solution"
    call stop_recon(stop_mpi,"get_solution_points_2d")
    !
  else if (iop_solpts == Chebyshev_Gauss_Lobatto) then
    !
    call nodes_chebyshev_gauss_lobatto(n_order,solpts_1d_qp)
    call weights_chebyshev_gauss_lobatto(n_order,solwts_1d_qp)
    write (iout,CGL_error) "Solution"
    call stop_recon(stop_mpi,"get_solution_points_2d")
    !
  else
    !
    write (iout,PNT_error) "solution"
    call stop_recon(stop_mpi,"get_solution_points_2d")
    !
  end if
  !
  ! Now get the 1D flux points
  !
  if (iop_flxpts == Legendre_Gauss) then
    !
    flxpts_1d_qp(:) = gauss_1d_qp(:)
    flxwts_1d_qp(:) = gauss_wts_1d_qp(:)
    !
  else if (iop_flxpts == Legendre_Gauss_Lobatto) then
    !
    flxpts_1d_qp(:) = lobatto_1d_qp(:)
    flxwts_1d_qp(:) = lobatto_wts_1d_qp(:)
    !
  else if (iop_flxpts == Chebyshev_Gauss) then
    !
    call nodes_chebyshev_gauss(n_sp1de,flxpts_1d_qp)
    call weights_chebyshev_gauss(n_sp1de,flxwts_1d_qp)
    write (iout,CG_error) "Flux"
    call stop_recon(stop_mpi,"get_solution_points_2d")
    !
  else if (iop_flxpts == Chebyshev_Gauss_Lobatto) then
    !
    call nodes_chebyshev_gauss_lobatto(n_order,flxpts_1d_qp)
    call weights_chebyshev_gauss_lobatto(n_order,flxwts_1d_qp)
    write (iout,CGL_error) "Flux"
    call stop_recon(stop_mpi,"get_solution_points_2d")
    !
  else
    !
    write (iout,PNT_error) "flux"
    call stop_recon(stop_mpi,"get_solution_points_2d")
    !
  end if
  !
  ! Get the coordinates of the base triangle in r-s space
  !
  if (iop_tripts == DG_TriPoints) then
    rs_qp(:,:) = TriNodes2D_DG(n_order,qzero)
  else if (iop_tripts == HT_TriPoints) then
    rs_qp(:,:) = TriNodes2D_HT(n_order,qzero)
  end if
  !
  ! r and s come out defined on the interval [0,1]
  !
  ! Get the triangle and quadrilateral weights at the solution points
  !
  weights2d_qp(1:n_sp2dt,1) = get_triangle_weights_at_solpts(rs_qp)
  !
  if (abs(error_loc) == 1) then
    !
    vp_qp(:) = merge( gauss_1d_qp(:) , radau_1d_qp(:) , error_loc > 0 )
    !
    weights2d_qp(1:n_sp2dq,2) = merge( gauss_wts_2d_qp(1:n_sp2dq) , &
                                       radau_wts_2d_qp(1:n_sp2dq) , &
                                       error_loc > 0 )
    !
  else
    !
    vp_qp(:) = solpts_1d_qp(:)
    !
    weights2d_qp(1:n_sp2dq,2) = get_quad_weights_at_solpts(solwts_1d_qp)
    !
  end if
  !
  ! Chop the quad precision arrays
  !
  radau_1d_qp       = chop_qp( radau_1d_qp )
  radau_wts_1d_qp   = chop_qp( radau_wts_1d_qp )
  radau_wts_2d_qp   = chop_qp( radau_wts_2d_qp )
  gauss_1d_qp       = chop_qp( gauss_1d_qp )
  gauss_wts_1d_qp   = chop_qp( gauss_wts_1d_qp )
  gauss_wts_2d_qp   = chop_qp( gauss_wts_2d_qp )
  lobatto_1d_qp     = chop_qp( lobatto_1d_qp )
  lobatto_wts_1d_qp = chop_qp( lobatto_wts_1d_qp )
  lobatto_wts_2d_qp = chop_qp( lobatto_wts_2d_qp )
  solpts_1d_qp      = chop_qp( solpts_1d_qp )
  solwts_1d_qp      = chop_qp( solwts_1d_qp )
  flxpts_1d_qp      = chop_qp( flxpts_1d_qp )
  flxwts_1d_qp      = chop_qp( flxwts_1d_qp )
  vp_qp             = chop_qp( vp_qp )
  rs_qp             = chop_qp( rs_qp )
  weights2d_qp      = chop_qp( weights2d_qp )
  !
  ! ##############################################
  ! ##############################################
  ! ##   CHOP THE QUAD PRECISION ARRAYS AND     ##
  ! ##  STORE IN THE WORKING PRECISION ARRAYS.  ##
  ! ##############################################
  ! ##############################################
  !
  radau_1d(:)       = chop( radau_1d_qp(:) )
  radau_wts_1d(:)   = chop( radau_wts_1d_qp(:) )
  radau_wts_2d(:)   = chop( radau_wts_2d_qp(:) )
  gauss_1d(:)       = chop( gauss_1d_qp(:) )
  gauss_wts_1d(:)   = chop( gauss_wts_1d_qp(:) )
  gauss_wts_2d(:)   = chop( gauss_wts_2d_qp(:) )
  lobatto_1d(:)     = chop( lobatto_1d_qp(:) )
  lobatto_wts_1d(:) = chop( lobatto_wts_1d_qp(:) )
  lobatto_wts_2d(:) = chop( lobatto_wts_2d_qp(:) )
  solpts_1d(:)      = chop( solpts_1d_qp(:) )
  solwts_1d(:)      = chop( solwts_1d_qp(:) )
  flxpts_1d(:)      = chop( flxpts_1d_qp(:) )
  flxwts_1d(:)      = chop( flxwts_1d_qp(:) )
  vp(:)             = chop( vp_qp(:) )
  rs(:,:)           = chop( rs_qp(:,:) )
  weights2d(:,:)    = chop( weights2d_qp(:,:) )
  !
  !
  ! ######################################################################
  ! ######################################################################
  ! ##  CONTINUE ON AND COMPUTE THE COORDINATES OF THE SOLUTION POINTS  ##
  ! ######################################################################
  ! ######################################################################
  !
  ! Save the number of solution points within a cell to pts_in_cell_ptr
  !
  allocate ( pts_in_cell_ptr(1:n_totcel+1) , source=0 , stat=ierr )
  call alloc_error("get_solution_points_2d","pts_in_cell_ptr",1,ierr)
  !
  do n = 1,n_totcel
    np = nodes_of_cell_ptr(n+1) - nodes_of_cell_ptr(n)
    pts_in_cell_ptr(n+1) = pts_in_cell_ptr(n) + merge(n_sp2dt,n_sp2dq,np==3)
  end do
  !
  ! Total number of interior solution points
  !
  n_solpts = pts_in_cell_ptr(ncell+1)
  !
  ! Total number of interior plus ghost cell solution points
  !
  n_totpts = last(pts_in_cell_ptr)
  !
  ! Now calculate the solution points and the
  ! Gauss solution points within each cell
  !
  allocate ( xyz(1:nr,1:n_totpts) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","xyz",1,ierr)
  allocate ( xyz_error(1:nr,1:n_totpts) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","xyz_error",1,ierr)
  allocate ( xyz_radau(1:nr,1:n_totpts) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","xyz_radau",1,ierr)
  allocate ( xyz_gauss(1:nr,1:n_totpts) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","xyz_gauss",1,ierr)
  allocate ( xyz_lobatto(1:nr,1:n_totpts) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","xyz_lobatto",1,ierr)
  !
  do i = 1,n_totcel
    !
    n1 = nodes_of_cell_ptr(i)+1
    n2 = nodes_of_cell_ptr(i+1)
    !
    p1 = pts_in_cell_ptr(i)+1
    p2 = pts_in_cell_ptr(i+1)
    !
    do j = n1,n2
      xy(:,j-n1+1) = xyz_nodes(:,nodes_of_cell(j))
    end do
    !
    xyz(:,p1:p2) = merge(get_solution_points_tria(rs       ,xy(1:nr,1:3)), &
                         get_solution_points_quad(solpts_1d,xy(1:nr,1:4)), &
                         n2-n1+1 == 3)
    !
    xyz_error(:,p1:p2) = get_solution_points_quad(vp,xy(1:nr,1:4))
    !
    xyz_gauss(:,p1:p2) = get_solution_points_quad(gauss_1d,xy(1:nr,1:4))
    xyz_radau(:,p1:p2) = get_solution_points_quad(radau_1d,xy(1:nr,1:4))
    xyz_lobatto(:,p1:p2) = get_solution_points_quad(lobatto_1d,xy(1:nr,1:4))
    !
  end do
  !
  ! Collect the solution point coordinates for the global grid
  !
  call collect_global_solpts
  !
end subroutine get_solution_points_2d
!
!###############################################################################
!
pure function get_solution_points_tria(rs,xy) result(return_value)
  !
  real(wp), dimension(:,:), intent(in) :: rs
  real(wp), dimension(:,:), intent(in) :: xy
  !
  real(wp), dimension(1:size(rs,1),1:size(rs,2)) :: return_value
  !
  integer   :: i
  real(wp) :: r,s
  !
continue
  !
  ! On entry, r and s both need to be defined on the interval [0,1]
  !
  do i = 1,size(rs,2)
    !
    r = rs(1,i)
    s = rs(2,i)
    !
    return_value(:,i) = (one-r-s)*xy(:,1) + r*xy(:,2) + s*xy(:,3)
    !
  end do
  !
end function get_solution_points_tria
!
!###############################################################################
!
pure function get_solution_points_quad(xi,xy) result(return_value)
  !
  real(wp), dimension(:),   intent(in) :: xi
  real(wp), dimension(:,:), intent(in) :: xy
  !
  real(wp), dimension(1:size(xy,1),1:size(xi)**2) :: return_value
  !
  integer   :: i,j,l
  real(wp) :: etam1,etap1,xim1,xip1
  !
continue
  !
  l = 1
  do i = 1,size(xi)
    !
    etam1 = one - xi(i)
    etap1 = one + xi(i)
    !
    do j = 1,size(xi)
      !
      xim1 = one - xi(j)
      xip1 = one + xi(j)
      !
      return_value(:,l) = c0_25 * ( etam1*xim1*xy(:,1) + &
                                    etam1*xip1*xy(:,2) + &
                                    etap1*xip1*xy(:,3) + &
                                    etap1*xim1*xy(:,4) )
      !
      l = l + 1
      !
    end do
    !
  end do
  !
end function get_solution_points_quad
!
!###############################################################################
!
pure function get_triangle_weights_at_solpts_DP(rs) result(return_value)
  !
  ! NOTE : This seems to fail for large n_order (>15)
  !
  use polynomial_mod, only : nodes_legendre_gauss, &
                             weights_legendre_gauss
  !
  real(dp), dimension(:,:), intent(in) :: rs
  !
  real(dp), dimension(1:size(rs,2)) :: return_value
  !
  integer  :: i,j,k,n,m,np
  integer  :: nspts,nfpts
  real(dp) :: g,h,x,y
  !
  integer, parameter :: nquad = 20
  !
  real(dp), dimension(1:nquad)                   :: xi
  real(dp), dimension(1:nquad)                   :: dxi
  real(dp), dimension(1:nquad)                   :: weights
  real(dp), dimension(1:size(rs,2),1:size(rs,2)) :: V
  !
  integer, parameter :: lp = dp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  ! On entry, r and s need to be defined on the interval [0,1]
  !
  nspts = size(rs,2)
  nfpts = np2nfp(nspts)
  !
  ! Compute the coefficients of the Vandermonde matrix for
  ! the Lagrange polynomials using monomial basis functions
  !
  do n = 1,nspts
    k = 0
    do i = 1,nfpts
      do j = 1,i
        k = k + 1
        V(k,n) = rs(1,n)**(i-j) * rs(2,n)**(j-1)
      end do
    end do
  end do
  !
  ! Now get the inverse of this matrix for use later
  !
  V = invert_matrix(V)
  !
  ! Get the Gauss quadrature points on the interval [-1,1] and
  ! the corresponding weights necessary to perform the integration
  !
  call nodes_legendre_gauss(nquad,xi,dxi)
  call weights_legendre_gauss(nquad,xi,dxi,weights)
  !
  ! Reset the quadrature points and weights to the interval [0,1]
  !
  xi(:) = half*(one+xi(:))
  weights(:) = half*weights(:)
  !
  ! Initialize the triangle weights to zero
  !
  return_value(:) = zero
  !
  ! Loop over each solution point performing a double integral on the
  ! corresponding Lagrange polynomial over the domain of the triangle
  ! i.e., ipt    : current solution point
  !       L(ipt) : Lagrange polynomial for ipt
  !
  !   /\ 1    /\ 1-xi
  !   |       |
  !   | d(xi) | d(eta) [L(ipt)]
  !   |       |
  !  \/ 0    \/ 0
  !
  do np = 1,nspts
    !
    ! Outer integral loop
    !
    outer_integral: do n = 1,nquad
      !
      ! Inner integral loop
      !
      h = zero
      inner_integral: do m = 1,nquad
        !
        x = xi(n)
        y = (one-xi(n))*xi(m)
        !
        ! Evaluate the polynomial at the point x,y
        !
        k = 0
        g = zero
        do i = 1,nfpts
          do j = 1,i
            k = k + 1
            g = g + V(np,k) * x**(i-j) * y**(j-1)
          end do
        end do
        !
        ! Add the weighted polynomial value to the inner integral result
        !
        h = h + (one-xi(n)) * weights(m) * g
        !
      end do inner_integral
      !
      ! Add the weighted inner integral result to the outer integral result
      !
      return_value(np) = return_value(np) + weights(n) * h
      !
    end do outer_integral
    !
  end do
  !
end function get_triangle_weights_at_solpts_DP
!
!###############################################################################
!
pure function get_triangle_weights_at_solpts_QP(rs) result(return_value)
  !
  ! NOTE : This seems to fail for large n_order (>15)
  !
  use polynomial_mod, only : nodes_legendre_gauss, &
                             weights_legendre_gauss
  !
  real(qp), dimension(:,:), intent(in) :: rs
  !
  real(qp), dimension(1:size(rs,2)) :: return_value
  !
  integer  :: i,j,k,n,m,np
  integer  :: nspts,nfpts
  real(qp) :: g,h,x,y
  !
  integer, parameter :: nquad = 20
  !
  real(qp), dimension(1:nquad)                   :: xi
  real(qp), dimension(1:nquad)                   :: dxi
  real(qp), dimension(1:nquad)                   :: weights
  real(qp), dimension(1:size(rs,2),1:size(rs,2)) :: V
  !
  integer, parameter :: lp = qp
  !
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: half = 0.5_lp
  real(lp), parameter :: one  = 1.0_lp
  !
continue
  !
  ! On entry, r and s need to be defined on the interval [0,1]
  !
  nspts = size(rs,2)
  nfpts = np2nfp(nspts)
  !
  ! Compute the coefficients of the Vandermonde matrix for
  ! the Lagrange polynomials using monomial basis functions
  !
  do n = 1,nspts
    k = 0
    do i = 1,nfpts
      do j = 1,i
        k = k + 1
        V(k,n) = rs(1,n)**(i-j) * rs(2,n)**(j-1)
      end do
    end do
  end do
  !
  ! Now get the inverse of this matrix for use later
  !
  V = invert_matrix(V)
  !
  ! Get the Gauss quadrature points on the interval [-1,1] and
  ! the corresponding weights necessary to perform the integration
  !
  call nodes_legendre_gauss(nquad,xi,dxi)
  call weights_legendre_gauss(nquad,xi,dxi,weights)
  !
  ! Reset the quadrature points and weights to the interval [0,1]
  !
  xi(:) = half*(one+xi(:))
  weights(:) = half*weights(:)
  !
  ! Initialize the triangle weights to zero
  !
  return_value(:) = zero
  !
  ! Loop over each solution point performing a double integral on the
  ! corresponding Lagrange polynomial over the domain of the triangle
  ! i.e., ipt    : current solution point
  !       L(ipt) : Lagrange polynomial for ipt
  !
  !   /\ 1    /\ 1-xi
  !   |       |
  !   | d(xi) | d(eta) [L(ipt)]
  !   |       |
  !  \/ 0    \/ 0
  !
  do np = 1,nspts
    !
    ! Outer integral loop
    !
    outer_integral: do n = 1,nquad
      !
      ! Inner integral loop
      !
      h = zero
      inner_integral: do m = 1,nquad
        !
        x = xi(n)
        y = (one-xi(n))*xi(m)
        !
        ! Evaluate the polynomial at the point x,y
        !
        k = 0
        g = zero
        do i = 1,nfpts
          do j = 1,i
            k = k + 1
            g = g + V(np,k) * x**(i-j) * y**(j-1)
          end do
        end do
        !
        ! Add the weighted polynomial value to the inner integral result
        !
        h = h + (one-xi(n)) * weights(m) * g
        !
      end do inner_integral
      !
      ! Add the weighted inner integral result to the outer integral result
      !
      return_value(np) = return_value(np) + weights(n) * h
      !
    end do outer_integral
    !
  end do
  !
end function get_triangle_weights_at_solpts_QP
!
!###############################################################################
!
pure function get_quad_weights_at_solpts_DP(weights) result(return_value)
  !
  real(dp), dimension(:), intent(in) :: weights
  !
  real(dp), dimension(1:size(weights)**2) :: return_value
  !
  integer :: i,j,k
  !
continue
  !
  k = 0
  do i = 1,size(weights)
    do j = 1,size(weights)
      k = k + 1
      return_value(k) = weights(i)*weights(j)
    end do
  end do
  !
end function get_quad_weights_at_solpts_DP
!
!###############################################################################
!
pure function get_quad_weights_at_solpts_QP(weights) result(return_value)
  !
  real(qp), dimension(:), intent(in) :: weights
  !
  real(qp), dimension(1:size(weights)**2) :: return_value
  !
  integer :: i,j,k
  !
continue
  !
  k = 0
  do i = 1,size(weights)
    do j = 1,size(weights)
      k = k + 1
      return_value(k) = weights(i)*weights(j)
    end do
  end do
  !
end function get_quad_weights_at_solpts_QP
!
!###############################################################################
!
subroutine get_grid_connectivity()
  !
  !.. Use Statements ..
  use geovar, only : bface
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:)   :: cells_with_node
  integer, allocatable, dimension(:)   :: cells_with_node_ptr
  integer, allocatable, dimension(:)   :: cells_surr_cell
  integer, allocatable, dimension(:)   :: cells_surr_cell_ptr
  !
continue
  !
  ! Get the cells surrounding each node
  !
  call find_cells_surrounding_each_node(cells_with_node,cells_with_node_ptr)
  !
  ! Get the cells surrounding each cell
  !
  call find_cells_surrounding_each_cell(cells_with_node,cells_with_node_ptr, &
                                        cells_surr_cell,cells_surr_cell_ptr)
  !
  ! Get the cells on each face
  !
  call find_cells_on_each_face(cells_surr_cell,cells_surr_cell_ptr)
  !
  ! Get the solution points on each face
  !
  call find_points_on_each_face
  !
  ! Get the indices for the flux points for each cell
  !
  call find_flxpts_of_each_cell
  !
  ! Get the periodic connectivity if there is any
  !
  if (any(bface(1,:) == bc_periodic)) &
    call get_periodic_connectivity
  !
  ! Get the unit cell subconnectivities for use when writing solution file
  !
  call get_unit_cell_subconnectivity
  !
  ! Get the grid metrics
  !
 !call get_metrics
  !
end subroutine get_grid_connectivity
!
!###############################################################################
!
subroutine find_cells_surrounding_each_node(cells_with_node,cells_with_node_ptr)
  !
  ! Find the cells that surround each node
  !
  ! For node i_node,
  ! cells_with_node(ibeg:iend) gives all of the grid cells that contain i_node
  ! where : ibeg = cells_with_node_ptr(i_node)+1
  !         iend = cells_with_node_ptr(i_node+1)
  !
  !.. Use Statements ..
  use geovar
  !
  !.. Formal Arguments ..
  integer, allocatable, dimension(:), intent(inout) :: cells_with_node
  integer, allocatable, dimension(:), intent(inout) :: cells_with_node_ptr
  !
  !.. Local Scalars ..
  integer :: i,ip,is,ibeg,iend,n,nc,ierr
  !
  !.. Local Allocatable Arrays ..
  logical, allocatable, dimension(:,:) :: ghost_msk
  !
continue
  !
  allocate ( cells_with_node_ptr(1:nnode+1) , source=0 , stat=ierr )
  call alloc_error("find_cells_surrounding_each_cell","cells_with_node_ptr",1,ierr)
  !
  ! Loop over the cells, storing 'ahead'
  !
  do n = 1,ncell
    !
    ibeg = nodes_of_cell_ptr(n)+1
    iend = nodes_of_cell_ptr(n+1)
    !
    do i = ibeg,iend
      ip = nodes_of_cell(i) + 1
      cells_with_node_ptr(ip) = cells_with_node_ptr(ip) + 1
    end do
    !
  end do
  !
  ! Add in the ghost cells. Make sure to only use the nodes
  ! of the ghost cell which lie on the boundary face.
  ! NOTE: We need a mask array here for the ghost cells to
  !       prevent this loop from finding a particular ghost
  !       cell more than once on an individual node.
  !
  allocate ( ghost_msk(ncell+1:n_totcel,1:nnode) , source=fals , stat=ierr )
  call alloc_error("find_cells_surrounding_each_cell","ghost_msk",1,ierr)
  !
  do n = 1,nfbnd
    !
    nc = bface(10,n)
    !
    ibeg = nodes_of_cell_ptr(nc)+1
    iend = nodes_of_cell_ptr(nc+1)
    !
    do i = ibeg,iend
      ip = nodes_of_cell(i)
      if (ip <= nnode) then
        if (.not. ghost_msk(nc,ip)) then
          cells_with_node_ptr(ip+1) = cells_with_node_ptr(ip+1) + 1
          ghost_msk(nc,ip) = true
        end if
      end if
    end do
    !
  end do
  !
  ! Reshuffle cells_with_node_ptr
  !
  do n = 2,size(cells_with_node_ptr)
    cells_with_node_ptr(n) = cells_with_node_ptr(n) + cells_with_node_ptr(n-1)
  end do
  !
  allocate ( cells_with_node(1:last(cells_with_node_ptr)) , source=0 , stat=ierr )
  call alloc_error("find_cells_surrounding_each_cell","cells_with_node",1,ierr)
  !
  ! Now store the cells surrounding each node in cells_with_node
  !
  do n = 1,ncell
    !
    ibeg = nodes_of_cell_ptr(n)+1
    iend = nodes_of_cell_ptr(n+1)
    !
    do i = ibeg,iend
      ip = nodes_of_cell(i)
      is = cells_with_node_ptr(ip) + 1
      cells_with_node_ptr(ip) = is
      cells_with_node(is) = n
    end do
    !
  end do
  !
  ! Add in the ghost cells. Make sure to only use the nodes
  ! of the ghost cell which lie on the boundary face. 
  ! NOTE: We need a mask array here for the ghost cells to
  !       prevent this loop from finding a particular ghost
  !       cell more than once on an individual node.
  !
  ghost_msk(:,:) = fals
  !
  do n = 1,nfbnd
    !
    nc = bface(10,n)
    !
    ibeg = nodes_of_cell_ptr(nc)+1
    iend = nodes_of_cell_ptr(nc+1)
    !
    do i = ibeg,iend
      ip = nodes_of_cell(i)
      if (ip <= nnode) then
        if (.not. ghost_msk(nc,ip)) then
          is = cells_with_node_ptr(ip) + 1
          cells_with_node_ptr(ip) = is
          cells_with_node(is) = nc
          ghost_msk(nc,ip) = true
        end if
      end if
    end do
    !
  end do
  !
  ! We are finished with ghost_msk so deallocate
  !
  deallocate ( ghost_msk , stat=ierr )
  call alloc_error("find_cells_surrounding_each_node","ghost_msk",2,ierr)
  !
  ! Finally, reorder cells_with_node_ptr
  !
  do n = size(cells_with_node_ptr),2,-1
    cells_with_node_ptr(n) = cells_with_node_ptr(n-1)
  end do
  cells_with_node_ptr(1) = 0
  !
end subroutine find_cells_surrounding_each_node
!
!###############################################################################
!
subroutine find_cells_surrounding_each_cell(cells_with_node,cells_with_node_ptr, &
                                            cells_surr_cell,cells_surr_cell_ptr)
  !
  ! Find the cells that surround each cell
  !
  ! For cell i_cell, ...
  ! cells_surr_cell(ibeg:iend) gives all of the grid cells
  !                            that are adjacent to i_cell
  ! where : ibeg = cells_surr_cell_ptr(i_cell)+1
  !         iend = cells_surr_cell_ptr(i_cell+1)
  !
  !.. Use Statements ..
  use geovar
  !
  !.. Formal Arguments ..
  integer, allocatable, dimension(:), intent(inout) :: cells_with_node
  integer, allocatable, dimension(:), intent(inout) :: cells_with_node_ptr
  integer, allocatable, dimension(:), intent(inout) :: cells_surr_cell
  integer, allocatable, dimension(:), intent(inout) :: cells_surr_cell_ptr
  !
  !.. Local Scalars ..
  integer :: i,i1,i2,ip1,ip2,ieadj
  integer :: j,j1,j2,n,n1,n2,nf,ierr
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:) :: lpoin
  !
continue
  !
  allocate ( lpoin(1:n_totnod) , source=0 , stat=ierr )
  call alloc_error("find_cells_surrounding_each_cell","lpoin",1,ierr)
  allocate ( cells_surr_cell_ptr(1:ncell+1) , source=0 , stat=ierr )
  call alloc_error("find_cells_surrounding_each_cell","cells_surr_cell_ptr",1,ierr)
  !
  ! Loop to set up cells_surr_cell_ptr
  !
  do n = 1,ncell
    cells_surr_cell_ptr(n+1) = cells_surr_cell_ptr(n) + nodes_of_cell_ptr(n+1) &
                                                      - nodes_of_cell_ptr(n)
  end do
  !
  allocate ( cells_surr_cell(1:last(cells_surr_cell_ptr)) , source=0 , stat=ierr )
  call alloc_error("find_cells_surrounding_each_cell","cells_surr_cell",1,ierr)
  !
  ! Loop to find cells_surr_cell
  !
  do n = 1,ncell
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
      do i = i1,i2
        j = cells_with_node(i)
        if (j == n) cycle
        j1 = nodes_of_cell_ptr(j)+1
        j2 = nodes_of_cell_ptr(j+1)
        if (sum(lpoin(nodes_of_cell(j1:j2))) == 2) then
          ieadj = j
          exit
        end if
      end do
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
  deallocate ( lpoin , stat=ierr )
  call alloc_error("find_cells_surrounding_each_cell","lpoin",2,ierr)
  !
  ! We are also finished with cells_with_node
  ! and cells_with_node_ptr so deallocate
  !
  deallocate ( cells_with_node , stat=ierr )
  call alloc_error("find_cells_surrounding_each_cell","cells_with_node",2,ierr)
  deallocate ( cells_with_node_ptr , stat=ierr )
  call alloc_error("find_cells_surrounding_each_cell","cells_with_node_ptr",2,ierr)
  !
end subroutine find_cells_surrounding_each_cell
!
!###############################################################################
!
subroutine find_cells_on_each_face(cells_surr_cell,cells_surr_cell_ptr)
  !
  ! For face i_face, ...
  ! nodes_on_face(1:2,i_face) : grid nodes defining i_face
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
  !       cells_on_face(5,:) = (n_qp1de+1,2*n_qp1de+1,1)
  !   for cells_on_face(4,:) = (    1    ,    2    ,   3   )
  !       cells_on_face(6,:) = (2*n_qp1de,3*n_qp1de,n_qp1de)
  ! For Quadrilaterals : 
  !   for cells_on_face(3,:) = (    1    ,     2     ,     3     ,4)
  !       cells_on_face(5,:) = (n_qp1de+1,2*n_qp1de+1,3*n_qp1de+1,1)
  !   for cells_on_face(4,:) = (    1    ,    2    ,    3    ,   4   )
  !       cells_on_face(6,:) = (2*n_qp1de,3*n_qp1de,4*n_qp1de,n_qp1de)
  !
  !.. Use Statements ..
  use geovar
  !
  !.. Formal Arguments ..
  integer, allocatable, dimension(:),   intent(inout) :: cells_surr_cell
  integer, allocatable, dimension(:),   intent(inout) :: cells_surr_cell_ptr
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,nf,ierr
  integer :: i,ic,ifa,ip1,ip2,ip3,jc
  !
  real(wp) :: x1,y1,x2,y2,dl,unx,uny,x13,y13,dot
  !
  !.. Local Arrays ..
  integer, dimension(1:4) :: lnode
  !
continue
  !
  ! Find the total number of faces
  ! and then allocate cells_on_face
  !
  nface = nfbnd
  do n = 1,ncell
    !
    n1 = cells_surr_cell_ptr(n)+1
    n2 = cells_surr_cell_ptr(n+1)
    !
    nface = nface + count( cells_surr_cell(n1:n2)>n .and. &
                           cells_surr_cell(n1:n2)<=ncell )
   !do nf = n1,n2
   !  i = cells_surr_cell(nf)
   !  if (i>n .and. i<=ncell) nface = nface + 1
   !end do
    !
  end do
  !
  allocate ( cells_on_face(1:6,1:nface) , source=0 , stat=ierr )
  call alloc_error("find_cells_on_each_face","cells_on_face",1,ierr)
  allocate ( nodes_on_face(1:2,1:nface) , source=0 , stat=ierr )
  call alloc_error("find_cells_on_each_face","nodes_on_face",1,ierr)
  !
  ! Get the cells on each boundary face
  !
  do nf = 1,nfbnd
    !
    ic = bface(2,nf)
    !
    n1 = cells_surr_cell_ptr(ic)+1
    n2 = cells_surr_cell_ptr(ic+1)
    !
    ifa = bface(4,nf)+n1-1
    !
    ip1 = nodes_of_cell(ifa)
    ip2 = nodes_of_cell(merge(ifa+1,n1,ifa/=n2))
    !
    cells_on_face(1,nf) = bface(2,nf)
    cells_on_face(2,nf) = bface(10,nf)
    !
    nodes_on_face(1,nf) = ip1
    nodes_on_face(2,nf) = ip2
    !
  end do
  !
  ! Get the cells for the interior faces
  !
  ifa = nfbnd
  do ic = 1,ncell
    !
    n1 = cells_surr_cell_ptr(ic)+1
    n2 = cells_surr_cell_ptr(ic+1)
    !
    do nf = n1,n2
      !
      jc = cells_surr_cell(nf)
      !
      if (jc>ic .and. jc<=ncell) then
        !
        ifa = ifa + 1
        !
        ip1 = nodes_of_cell(nf)
        ip2 = nodes_of_cell(merge(nf+1,n1,nf/=n2))
        !
        cells_on_face(1,ifa) = ic
        cells_on_face(2,ifa) = jc
        !
        nodes_on_face(1,ifa) = ip1
        nodes_on_face(2,ifa) = ip2
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
  deallocate ( cells_surr_cell , stat=ierr )
  call alloc_error("find_cells_on_each_face","cells_surr_cell",2,ierr)
  deallocate ( cells_surr_cell_ptr , stat=ierr )
  call alloc_error("find_cells_on_each_face","cells_surr_cell_ptr",2,ierr)
  !
  ! For the face vector going from nodes_on_face(1,nf) to nodes_on_face(2,nf),
  ! this is a test to make sure that cells_on_face(1,nf) is to the left of
  ! this face vector and cells_on_face(2,nf) is on the right
  !
  ic = 0
  do nf = nfbnd+1,nface
    !
    ip1 = nodes_on_face(1,nf)
    ip2 = nodes_on_face(2,nf)
    !
    x1 = xyz_nodes(1,ip1)
    y1 = xyz_nodes(2,ip1)
    x2 = xyz_nodes(1,ip2)
    y2 = xyz_nodes(2,ip2)
    !
    dl = sqrt( (x2-x1)**2 + (y2-y1)**2 )
    !
    ! Get the negative unit normal for this face
    !
    unx = (y2-y1)/dl
    uny = (x1-x2)/dl
    !
    jc = cells_on_face(2,nf)
    n1 = nodes_of_cell_ptr(jc)+1
    n2 = nodes_of_cell_ptr(jc+1)
    !
    ! Find a cell node not on the face
    !
    if (n2-n1+1 == 3) then
      !
      ip3 = sum( nodes_of_cell(n1:n2) ) - (ip1+ip2)
      !
    else if (n2-n1+1 == 4) then
      !
      if (nodes_of_cell(n1)/=ip1 .and. nodes_of_cell(n1)/=ip2) then
        ip3 = nodes_of_cell(n1)
      else if (nodes_of_cell(n2)/=ip1 .and. nodes_of_cell(n2)/=ip2) then
        ip3 = nodes_of_cell(n2)
      else
        ip3 = nodes_of_cell(n1+1)
      end if
      !
    else
      !
      write (*,*) "Error: Should not be here in sub find_cells_on_each_face!"
      !
    end if
    !
    ! Vector from face node to non-face node
    !
    x13 = xyz_nodes(1,ip3) - xyz_nodes(1,ip1)
    y13 = xyz_nodes(2,ip3) - xyz_nodes(2,ip1)
    !
    ! Dot product of previous vector with negative face normal
    !
    dot = unx*x13 + uny*y13
    !
    ! If dot is negative, then the cell stored in cells_on_face(2,nf)
    ! is actually the cell to the left of the face vector and needs
    ! to be stored in cells_on_face(1,nf)
    !
    if (dot < zero) then
      write (*,'(a,i0,a,i0,2x,i0)') "Error on face - ",nf,",  nodes ",ip1,ip2
      ic = 1
    end if
    !
  end do
  if (ic == 1) call stop_recon(stop_mpi,"find_cells_on_each_face")
  !
  ! On each face, we need to find the vertex of triangular
  ! cells that is not located on the face. For quadrilateral
  ! cells, this will be set to zero for now.
  !
  do nf = 1,nface
    !
    ! Nodes on the current face
    !
    ip1 = nodes_on_face(1,nf)
    ip2 = nodes_on_face(2,nf)
    !
    ! For the left cell on the face
    !
    ic = cells_on_face(1,nf)
    n1 = nodes_of_cell_ptr(ic)+1
    n2 = nodes_of_cell_ptr(ic+1)
    !
    ! Find which face of the left cell is on the current face
    ! and the starting index used in the correction procedure
    !
    lnode = 0
    do i = n1,n2
      if (nodes_of_cell(i) == ip1) lnode(i-n1+1) = 1
      if (nodes_of_cell(i) == ip2) lnode(i-n1+1) = 1
    end do
    !
    if (sum(lnode(1:2)) == 2) then
      cells_on_face(3,nf) = 1
      cells_on_face(5,nf) = 1
    else if (sum(lnode(2:3)) == 2) then
      cells_on_face(3,nf) = 2
      cells_on_face(5,nf) = n_sp1de + 1
    else if (sum(lnode(3:4)) == 2) then
      cells_on_face(3,nf) = 3
      cells_on_face(5,nf) = 2*n_sp1de + 1
    else
      cells_on_face(3,nf) = 4
      cells_on_face(5,nf) = merge(2,3,n2-n1+1==3)*n_sp1de + 1
    end if
    !
    ! Now the right cell on the current face
    !
    if (nf <= nfbnd) then
      !
     !! This is a boundary face so the right cell
     !! needs to be handled a little differently
     !!
     !cells_on_face(4,nf) = 1
     !cells_on_face(6,nf) = n_sp1de
     !cells_on_face(6,nf) = n_qp1de
      !
      ! For boundary faces, I think everything should be set up
      ! properly so that we shouldnt have to differentiate here
      ! between boundary and interior faces for the right cell.
      ! However, we will keep them separate for the time being...
      !
      jc = cells_on_face(2,nf)
      n1 = nodes_of_cell_ptr(jc)+1
      n2 = nodes_of_cell_ptr(jc+1)
      !
      ! Find which face of the ghost cell is on the current face
      ! and the starting index used in the correction procedure.
      ! NOTE: For boundary faces, all but non-periodic communication
      !       boundaries should result in:
      !         cells_on_face(4,nf) = 1
      !         cells_on_face(6,nf) = 1
      !       For non-periodic communication boundaries, the face of the
      !       ghost cell on the boundary can vary due to the possibility
      !       of multiple boundary faces sharing a single ghost cell.
      !       We will explicitly match the boundary face nodes to determine
      !       the face orientation for all boundary face types in an effort
      !       to be thorough.
      !
      lnode = 0
      do i = n1,n2
        if (nodes_of_cell(i) == ip1) lnode(i-n1+1) = 1
        if (nodes_of_cell(i) == ip2) lnode(i-n1+1) = 1
      end do
      !
      if (sum(lnode(1:2)) == 2) then
        cells_on_face(4,nf) = 1
        cells_on_face(6,nf) = n_sp1de
      else if (sum(lnode(2:3)) == 2) then
        cells_on_face(4,nf) = 2
        cells_on_face(6,nf) = 2*n_sp1de
      else if (sum(lnode(3:4)) == 2) then
        cells_on_face(4,nf) = 3
        cells_on_face(6,nf) = 3*n_sp1de
      else
        cells_on_face(4,nf) = 4
        cells_on_face(6,nf) = merge(3,4,n2-n1+1==3)*n_sp1de
      end if
      !
    else
      !
      ! For the right cell on an interior face
      !
      jc = cells_on_face(2,nf)
      n1 = nodes_of_cell_ptr(jc)+1
      n2 = nodes_of_cell_ptr(jc+1)
      !
      ! Find which face of the right cell is on the current face
      ! and the starting index used in the correction procedure
      !
      lnode = 0
      do i = n1,n2
        if (nodes_of_cell(i) == ip1) lnode(i-n1+1) = 1
        if (nodes_of_cell(i) == ip2) lnode(i-n1+1) = 1
      end do
      !
      if (sum(lnode(1:2)) == 2) then
        cells_on_face(4,nf) = 1
        cells_on_face(6,nf) = n_sp1de
      else if (sum(lnode(2:3)) == 2) then
        cells_on_face(4,nf) = 2
        cells_on_face(6,nf) = 2*n_sp1de
      else if (sum(lnode(3:4)) == 2) then
        cells_on_face(4,nf) = 3
        cells_on_face(6,nf) = 3*n_sp1de
      else
        cells_on_face(4,nf) = 4
        cells_on_face(6,nf) = merge(3,4,n2-n1+1==3)*n_sp1de
      end if
      !
    end if
    !
  end do
  !
 !write (110,*)
 !write (110,*) "Boundary Faces"
 !write (110,*)
 !do nf = 1,nface
 !  write (110,1000) nf, cells_on_face(1:5:2,nf),cells_on_face(2:6:2,nf)
 !  if (nf == nfbnd) then
 !    write (110,*)
 !    write (110,*) "Interior Faces"
 !    write (110,*)
 !  end if
 !end do
 !1000 format ("Face ",i4,": ",2(6x,i4,1x,i1,1x,i2))
  !
end subroutine find_cells_on_each_face
!
!###############################################################################
!
subroutine find_points_on_each_face()
  !
  ! For face i_face,
  ! pts_on_face(1,1:n_sp1de,i_face) : solution points of the left  cell on i_face
  ! pts_on_face(2,1:n_sp1de,i_face) : solution points of the right cell on i_face
  !
  !.. Use Statements ..
  use geovar
  !
  !.. Local Scalars ..
  integer :: i,j,k,ic,jc,ierr
  integer :: kp,n1,n2,nf,p1,p2
  !
  !.. Local Arrays ..
  real(wp), dimension(1:2) :: dxy,xyfp,xykp
  !
continue
  !
  ! Allocate pts_on_face array and initialize to zero
  !
  allocate ( pts_on_face(1:2,1:n_sp1de,1:nface) , source=0 , stat=ierr )
  call alloc_error("find_points_on_each_face","pts_on_face",1,ierr)
  !
  ! Identify which solution points of a cell are on each edge
  !
  do nf = 1,nface
    !
    ic = cells_on_face(1,nf) ! cell i on face
    jc = cells_on_face(2,nf) ! cell j on face
    !
    n1 = nodes_on_face(1,nf) ! index of node 1 on face
    n2 = nodes_on_face(2,nf) ! index of node 2 on face
    !
    ! Vector from node 1 to node 2 of face
    !
    dxy(:) = xyz_nodes(:,n2) - xyz_nodes(:,n1)
    !
    face_solpts_loop: do k = 1,n_sp1de
      !
      ! Compute the coordinates of the k-th solution point on this face.
      ! Start at node 1 on the face, and move in the direction of the
      ! face vector, dxy, to the k-th Lobatto point.
      !
      xyfp(:) = xyz_nodes(:,n1) + (solpts_1d(k)+one)/two * dxy
      !
      ! Now loop over the solution points within cell ic to find the
      ! one that corresponds to the k-th solution point on this face.
      !
      p1 = pts_in_cell_ptr(ic)+1
      p2 = pts_in_cell_ptr(ic+1)
      !
      ic_solpts_loop: do kp = p1,p2
        !
        xykp(:) = xyz(:,kp)  ! x,y coordinates of solution point kp in cell ic
        !
        ! Test if dot product of xykp-xyfp with itself is zero to
        ! determine if this is the solution point on the face
        !
        if (selfdot(xykp-xyfp) < eps9) then
          pts_on_face(1,k,nf) = kp
          exit ic_solpts_loop
        end if
        !
      end do ic_solpts_loop
      !
      ! Now loop over the solution points within cell jc to find the
      ! one that corresponds to the k-th solution point on this face.
      !
      p1 = pts_in_cell_ptr(jc)+1
      p2 = pts_in_cell_ptr(jc+1)
      !
      jc_solpts_loop: do kp = p1,p2
        !
        xykp = xyz(:,kp)  ! x,y coordinates of solution point kp in cell jc
        !
        ! Test if dot product of xykp-xyfp with itself is zero to
        ! determine if this is the solution point on the face
        !
        if (selfdot(xykp-xyfp) < eps9) then
          pts_on_face(2,k,nf) = kp
          exit jc_solpts_loop
        end if
        !
      end do jc_solpts_loop
      !
    end do face_solpts_loop
    !
  end do
  !
end subroutine find_points_on_each_face
!
!###############################################################################
!
subroutine find_flxpts_of_each_cell()
  !
  !.. Use Statements ..
  use geovar, only : nfbnd,ncell,n_totcel,nface,n_sp1de,n_fp1de,nflxprm
  use geovar, only : flxpts_of_cell_ptr,flxpts_of_cell,bface
  use geovar, only : cells_on_face,nodes_of_cell_ptr
  !
  !.. Local Scalars ..
  integer :: n,nf,ierr,k
  integer :: jl,lf1,lf2,lk1,lk,lfp,lkbg
  integer :: jr,rf1,rf2,rk1,rk,rfp,rkbg
  !
  !.. Local Parameters ..
  integer, parameter :: left_side = 1
  integer, parameter :: right_side = 2
  !
continue
  !
  ! Create the pointer for the flxpts_of_cell array
  !
  allocate ( flxpts_of_cell_ptr(1:ncell+1) , source=0 , stat=ierr )
  call alloc_error("find_flxpts_of_each_cell","flxpts_of_cell_ptr",1,ierr)
  !
  ! Set up the flux points pointer array for each interior cell
  ! NOTE: This array is only used for interior cells so it does
  !       not need to account for ghost cells
  !
  do n = 1,ncell
    k = nodes_of_cell_ptr(n+1) - nodes_of_cell_ptr(n)
    flxpts_of_cell_ptr(n+1) = flxpts_of_cell_ptr(n) + k*n_fp1de
  end do
  !
  ! Allocate the flxpts_of_cell array
  !
  allocate ( flxpts_of_cell(1:nflxprm,1:last(flxpts_of_cell_ptr)) , &
             source=0 , stat=ierr )
  call alloc_error("find_flxpts_of_each_cell","flxpts_of_cell",1,ierr)
  !
  ! Loop through the boundary faces and fill
  ! the left part of the flxpts_of_cell array
  !
  do nf = 1,nfbnd
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    !
    lf1 = flxpts_of_cell_ptr(jl)+1 ! beginning flxpt index for left  cell
    lf2 = flxpts_of_cell_ptr(jl+1) ! ending    flxpt index for left  cell
    !
    lkbg = (cells_on_face(5,nf)-1)/n_sp1de
    !
    lk1 = lkbg*n_fp1de + 1 ! beginning flxpt on current face for left  cell
    !
    ! Loop through each flux point on the current face and get
    ! the flux point information for both cells on this face
    !
    do k = 1,n_fp1de
      !
      lk = lk1-1 + k ! idx of current flxpt for left  cell (from 1 to lf2-lf1+1)
      !
      lfp = lf1-1 + lk ! idx of current flxpt for left  cell (from lf1 to lf2)
      !
      flxpts_of_cell(1,lfp) = nf
      flxpts_of_cell(2,lfp) = k
      flxpts_of_cell(3,lfp) = lk
      flxpts_of_cell(4,lfp) = left_side
      !
    end do
    !
  end do
  !
  ! Loop through each interior face to fill the
  ! remaining parts of the flxpts_of_cell array
  !
  do nf = nfbnd+1,nface
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    lf1 = flxpts_of_cell_ptr(jl)+1 ! beginning flxpt index for left  cell
    lf2 = flxpts_of_cell_ptr(jl+1) ! ending    flxpt index for left  cell
    !
    rf1 = flxpts_of_cell_ptr(jr)+1 ! beginning flxpt index for right cell
    rf2 = flxpts_of_cell_ptr(jr+1) ! ending    flxpt index for right cell
    !
    lkbg = (cells_on_face(5,nf)-1)/n_sp1de
    rkbg = (cells_on_face(6,nf)-0)/n_sp1de
    !
    lk1 = lkbg*n_fp1de + 1 ! beginning flxpt on current face for left  cell
    rk1 = rkbg*n_fp1de     ! beginning flxpt on current face for right cell
    !
    ! Loop through each flux point on the current face and get
    ! the flux point information for both cells on this face
    !
    do k = 1,n_fp1de
      !
      lk = lk1-1 + k ! idx of current flxpt for left  cell (from 1 to lf2-lf1+1)
      rk = rk1+1 - k ! idx of current flxpt for right cell (from 1 to rf2-rf1+1)
      !
      lfp = lf1-1 + lk ! idx of current flxpt for left  cell (from lf1 to lf2)
      rfp = rf1-1 + rk ! idx of current flxpt for right cell (from rf1 to rf2)
      !
      flxpts_of_cell(1,lfp) = nf
      flxpts_of_cell(2,lfp) = k
      flxpts_of_cell(3,lfp) = lk
      flxpts_of_cell(4,lfp) = left_side
      !
      flxpts_of_cell(1,rfp) = nf
      flxpts_of_cell(2,rfp) = k
      flxpts_of_cell(3,rfp) = rk
      flxpts_of_cell(4,rfp) = right_side
      !
    end do
    !
  end do
  !
 !write (121+mypnum,6)
 !do nf = 1,nfbnd
 !  write (121+mypnum,7) nf,(bface(k,nf),k=1,size(bface,1))
 !end do
 !write (101+mypnum,1)
 !do nf = 1,ncell
 !  !
 !  lf1 = flxpts_of_cell_ptr(nf)+1
 !  lf2 = flxpts_of_cell_ptr(nf+1)
 !  !
 !  do k = lf1,lf2
 !    if (flxpts_of_cell(3,k) /= k-lf1+1) write (iout,*) "error",nf,k
 !    if (k == lf1) then
 !      if (flxpts_of_cell(1,k) <= nfbnd) then
 !        write (101+mypnum,8) nf,flxpts_of_cell(1:4,k)
 !      else
 !       !write (101+mypnum,2) nf,flxpts_of_cell(1:4,k)
 !      end if
 !    else if (mod(k-lf1+n_fp1de,n_fp1de) == 0) then
 !      if (flxpts_of_cell(1,k) <= nfbnd) then
 !        write (101+mypnum,8) nf,flxpts_of_cell(1:4,k)
 !      else
 !       !write (101+mypnum,3) flxpts_of_cell(1:4,k)
 !      end if
 !    else
 !      if (flxpts_of_cell(1,k) <= nfbnd) then
 !        write (101+mypnum,4) flxpts_of_cell(2:3,k)
 !      else
 !       !write (101+mypnum,4) flxpts_of_cell(2:3,k)
 !      end if
 !    end if
 !  end do
 !  write (101+mypnum,5)
 !end do
  !
  1 format (" Cell  Face  FaceFlxPt  CellFlxPt  FaceSide",/, &
              "----------------------------------------------")
  2 format (i5,i6,i8,i11,i10)
  3 format (5x,i6,i8,i11,i10)
  4 format (5x,6x,i8,i11,10x)
  5 format ("----------------------------------------------")
  6 format (" nf   bc        hc htype  hloc      adjf  adjh  adjp  adjt  adjl        gc")
  7 format (2i4,4x,3i6,4x,5i6,4x,i6)
  8 format (i5,"*",i4,"*",i8,i11,i10)
  9 format (5x,"*",i4,"*",i8,i11,i10)
  !
end subroutine find_flxpts_of_each_cell
!
!###############################################################################
!
subroutine get_periodic_connectivity()
  !
  ! Routine to establish the connectivity between boundary
  ! cells connected by a periodic boundary condition
  !
  !.. Use Statements ..
  use geovar
  use metvar
  !
  !.. Local Scalars ..
  integer   :: g1,g2,gc,nf,kg,gbeg,ierr
  integer   :: p1,p2,pc,pf,kp,kpts
  real(wp) :: theta
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: rgfc,rnrmg
  real(wp), dimension(1:nr) :: rpfc,rnrmp
  real(wp), dimension(1:nr,1:n_sp2dq) :: xyg,xyp
  !
continue
  !
  !
  ! Allocate and initialize the pointer array for pts_of_periodic
  !
  allocate ( pts_of_periodic_ptr(1:nfbnd+1) , source=0 , stat=ierr )
  call alloc_error("get_periodic_connectivity","pts_of_periodic_ptr",1,ierr)
  !
  ! Loop over the boundary faces and for each periodic boundary condition,
  ! store ahead the number of solution points for the corresponding cell
  !
  do nf = 1,nfbnd
    !
    ! Only interested in periodic boundary faces
    !
    if (bface(1,nf) /= bc_periodic) cycle
    !
    pf = bface(5,nf) ! index of connecting periodic boundary face
    pts_of_periodic_ptr(nf+1) = merge(n_sp2dt,n_sp2dq,bface(3,pf)==3)
    !
  end do
  !
  ! Reshuffle pts_of_periodic_ptr and allocate pts_of_periodic
  !
  do nf = 2,size(pts_of_periodic_ptr)
    pts_of_periodic_ptr(nf) = pts_of_periodic_ptr(nf) + pts_of_periodic_ptr(nf-1)
  end do
  !
  allocate ( pts_of_periodic(1:last(pts_of_periodic_ptr)) , source=0 , stat=ierr )
  call alloc_error("get_periodic_connectivity","pts_of_periodic",1,ierr)
  !
  ! Loop over the boundary faces and for each periodic boundary condition,
  ! find the matching solution points within the corresponding
  ! periodic cell for this ghost cell and store in pts_of_periodic
  !
  do nf = 1,nfbnd
    !
    ! Only do this for periodic boundary faces
    !
    if (bface(1,nf) /= bc_periodic) cycle
    !
    ! Beginning index of pts_of_periodic for face nf
    !
    gbeg = pts_of_periodic_ptr(nf)+1
    !
    ! Cells associated with this periodic boundary condition
    !
    pf = bface(5,nf)  ! index of corresponding boundary face to nf
    pc = bface(2,pf)  ! index of corresponding periodic cell
    gc = bface(10,nf) ! index of ghost cell
    !
    ! Beginning and ending indices for solution points in ghost cell
    !
    g1 = pts_in_cell_ptr(gc)+1
    g2 = pts_in_cell_ptr(gc+1)
    !
    ! Beginning and ending indices for solution points in periodic cell
    !
    p1 = pts_in_cell_ptr(pc)+1
    p2 = pts_in_cell_ptr(pc+1)
    !
    ! Centers of the boundary faces
    !
    rgfc(:) = half * sum( xyz_nodes(:,nodes_on_face(:,nf)) , dim=2 )
    rpfc(:) = half * sum( xyz_nodes(:,nodes_on_face(:,pf)) , dim=2 )
    !
    ! Number of solution points in cells gc and pc
    !
    kpts = p2 - p1 + 1
    !
    ! Coordinates of solution points in cells gc and pc with
    ! the periodic boundary face center located at the origin
    !
    do kp = 1,nr
      xyg(kp,1:kpts) = xyz(kp,g1:g2) - rgfc(kp)
      xyp(kp,1:kpts) = xyz(kp,p1:p2) - rpfc(kp)
    end do
    !
    ! Outward unit normal for boundary face nf
    !
    rnrmg(1) = xyz_nodes(2,nodes_on_face(2,nf)) - xyz_nodes(2,nodes_on_face(1,nf))
    rnrmg(2) = xyz_nodes(1,nodes_on_face(1,nf)) - xyz_nodes(1,nodes_on_face(2,nf))
    rnrmg(:) = unit_vector(rnrmg)
    !
    ! Outward unit normal for boundary face pf
    rnrmp(1) = xyz_nodes(2,nodes_on_face(2,pf)) - xyz_nodes(2,nodes_on_face(1,pf))
    rnrmp(2) = xyz_nodes(1,nodes_on_face(1,pf)) - xyz_nodes(1,nodes_on_face(2,pf))
    rnrmp(:) = unit_vector(rnrmp)
    !
    ! The angle between the two unit normals
    !
    theta = acos(dot(rnrmg,rnrmp))
    !
    ! Rotate the ghost cell solution points (theta+180) degrees
    ! around the origin (also the boundary face center). After
    ! this, the periodic and ghost cells should be oriented so
    ! that all the solution points match up between cells.
    !
    do kp = 1,kpts
      xyg(:,kp) = rotate_2d( xyg(:,kp) , theta+pi )
    end do
    !
    ! Find the matching solution points in the periodic cell
    ! for each of the solution points in the ghost cell
    !
    ghost_cell: do kg = 1,kpts
      periodic_cell: do kp = 1,kpts
        !
        if (selfdot( xyp(:,kp)-xyg(:,kg) ) <= eps9) then
          pts_of_periodic(gbeg+kg-1) = p1 + kp - 1
          cycle ghost_cell
        end if
        !
      end do periodic_cell
    end do ghost_cell
    !
  end do 
  !
end subroutine get_periodic_connectivity
!
!###############################################################################
!
subroutine get_unit_cell_subconnectivity()
  !
  !.. Use Statements ..
  use geovar, only : n_order,o_order,tria_connectivity,quad_connectivity
  !
  !.. Local Scalars ..
  integer :: j1,j2,j3,j4,n,m,sk,ierr,no1d,noc
  !
  !.. Local Arrays ..
  integer, allocatable, dimension(:,:) :: counter
  !
continue
  !
  if (o_order <= 0) then
    no1d = n_order
  else
    no1d = o_order
  end if
  !
  no1d = max(1,no1d)
  !
  noc = no1d*no1d
  !
  ! Allocate tria_connectivity and quad_connectivity
  !
  allocate ( counter(1:no1d+1,1:no1d+1) , source=0 , stat=ierr )
  call alloc_error("get_unit_cell_subconnectivity","counter",1,ierr)
  allocate ( tria_connectivity(1:4,1:noc) , source=0 , stat=ierr )
  call alloc_error("get_unit_cell_subconnectivity","tria_connectivity",1,ierr)
  allocate ( quad_connectivity(1:4,1:noc) , source=0 , stat=ierr )
  call alloc_error("get_unit_cell_subconnectivity","quad_connectivity",1,ierr)
  !
  ! Get the subconnectivity for the unit triangle
  !
  sk = 0
  do n = 1,no1d+1
    do m = 1,no1d+2-n
      sk = sk + 1
      counter(n,m) = sk
    end do
  end do
  !
  sk = 0
  do n = 1,no1d+1
    do m = 1,no1d+1-n
      !
      j1 = counter(n,m)
      j2 = counter(n,m+1)
      j3 = counter(n+1,m)
      j4 = counter(n+1,m+1)
      !
      sk = sk + 1
      tria_connectivity(1,sk) = j1
      tria_connectivity(2,sk) = j2
      tria_connectivity(3,sk) = j3
      tria_connectivity(4,sk) = j3
      !
      if (j4 /= 0) then
        sk = sk + 1
        tria_connectivity(1,sk) = j2
        tria_connectivity(2,sk) = j4
        tria_connectivity(3,sk) = j3
        tria_connectivity(4,sk) = j3
      end if
      !
    end do
  end do
  !
  ! Get the subconnectivity for the unit square
  !
  sk = 0
  do n = 1,no1d
    do m = 1,no1d
      !
      sk = sk + 1
      quad_connectivity(1,sk) = (no1d+1)*(n-1) + m
      quad_connectivity(2,sk) = (no1d+1)*(n-1) + m + 1
      quad_connectivity(3,sk) = (no1d+1)*(n  ) + m + 1
      quad_connectivity(4,sk) = (no1d+1)*(n  ) + m
      !
    end do
  end do
  !
  ! Now subtract one from all values of the subconnectivities. This
  ! is done so that the values go from (0:n_sp2dt-1) or (0:n_sp2dq-1)
  ! depending on the cell type and then this can simply be added
  ! to the beginning solution point index of pts_in_cell_ptr(icell)+1.
  !
  tria_connectivity(:,:) = tria_connectivity(:,:) - 1
  quad_connectivity(:,:) = quad_connectivity(:,:) - 1
  !
end subroutine get_unit_cell_subconnectivity
!
!###############################################################################
!
subroutine get_metrics()
  !
  !.. Use Statements ..
  use geovar
  use ovar, only : local_volume,total_volume
  use metvar
  !
  !.. Local Scalars ..
  integer   :: i,j,k,l,ierr
  integer   :: n1,n2,nf,np,p1,p2
  real(wp) :: tvol
  real(r8)  :: loc_vol,tot_vol
  real(wp) :: val
 !real(wp)  :: dxdr,dxds,dydr,dyds
  !
  !.. Local Arrays ..
  real(wp), dimension(1:2,1:4) :: xy
 !real(wp), dimension(1:n_qp2dq) :: qdtr_jac
 !real(wp), dimension(1:n_qp2dq) :: qdtr_wts
 !real(wp), dimension(1:n_qp2dq) :: dxdr
 !real(wp), dimension(1:n_qp2dq) :: dxds
 !real(wp), dimension(1:n_qp2dq) :: dydr
 !real(wp), dimension(1:n_qp2dq) :: dyds
  !
continue
  !
 !qdtr_wts(1:n_qp2dq) = get_quad_weights_at_solpts(qdtwts_1d)
 !!
 !do i = 110,150,10
 !  write (i+3,*) "using matmul"
 !  write (i+4,*) "using metrics loop"
 !end do
  !
  ! Allocate metrics arrays
  !
  allocate ( dx(1:n_totpts) , stat=ierr )
  call alloc_error("get_metrics","dx",1,ierr)
  allocate ( dy(1:n_totpts) , stat=ierr )
  call alloc_error("get_metrics","dy",1,ierr)
  allocate ( geofa(1:nr,1:n_sp1de,1:nface) , source=zero , stat=ierr )
  call alloc_error("get_metrics","geofa",1,ierr)
  allocate ( jacobian(1:n_totpts) , source=zero , stat=ierr )
  call alloc_error("get_metrics","jacobian",1,ierr)
  allocate ( dJ(1:n_totpts) , stat=ierr )
  call alloc_error("get_metrics","dJ",1,ierr)
  allocate ( volj(1:n_totcel) , source=zero , stat=ierr )
  call alloc_error("get_metrics","volj",1,ierr)
  !
  !
  local_volume = zero
  total_volume = zero
  !
  ! Loop over the cells to get the metrics at each solution point
  !
  metrics_loop: do i = 1,n_totcel
    !
    n1 = nodes_of_cell_ptr(i)+1
    n2 = nodes_of_cell_ptr(i+1)
    !
    p1 = pts_in_cell_ptr(i)+1
    p2 = pts_in_cell_ptr(i+1)
    !
    ! Load the nodes defining the cell into a temporary array
    !
    do j = n1,n2
      xy(:,j-n1+1) = xyz_nodes(:,nodes_of_cell(j))
    end do
    !
    ! Compute at each solution point (depending on the cell type)
    ! x metrics - dx/dr,dx/ds,d2x/dr2,d2x/ds2,d2x/drds,
    ! y metrics - dy/dr,dy/ds,d2y/dr2,d2y/ds2,d2y/drds,
    ! Jacobians - J,dJ/dr,dJ/ds
    !
    if (n2-n1+1 == 3) then
      !
      ! Metrics for triangular cell
      !
      tvol = zero
      !
      do j = p1,p2
        !
        ! These are correct if the reference triangle is defined
        ! by the three nodes (-1,-1) - (1,-1) - (-1,1)
        !
       !dx(j)%dr = half * (xy(1,2) - xy(1,1))
       !dy(j)%dr = half * (xy(2,2) - xy(2,1))
        !
       !dx(j)%ds = half * (xy(1,3) - xy(1,1))
       !dy(j)%ds = half * (xy(2,3) - xy(2,1))
        !
        ! These are correct if the reference triangle is defined
        ! by the three nodes (0,0) - (1,0) - (0,1)
        !
        dx(j)%dr = xy(1,2) - xy(1,1)
        dy(j)%dr = xy(2,2) - xy(2,1)
        !
        dx(j)%ds = xy(1,3) - xy(1,1)
        dy(j)%ds = xy(2,3) - xy(2,1)
        !
        dx(j)%dr2 = zero
        dy(j)%dr2 = zero
        !
        dx(j)%ds2 = zero
        dy(j)%ds2 = zero
        !
        dx(j)%drds = zero
        dy(j)%drds = zero
        !
        ! J = ( dx/dr * dy/ds ) - ( dx/ds * dy/dr )
        !
        jacobian(j) = dx(j)%dr*dy(j)%ds - dx(j)%ds*dy(j)%dr
        !
        ! dJ/dr = [ ( d2x/dr2 * dy/ds ) + ( dx/dr * d2y/drds ) ] - 
        !         [ ( d2x/drds * dy/dr ) + ( dx/ds * d2y/dr2 ) ]
        !
        dJ(j)%dr = ( dx(j)%dr2 *dy(j)%ds + dx(j)%dr*dy(j)%drds ) - &
                   ( dx(j)%drds*dy(j)%dr + dx(j)%ds*dy(j)%dr2  )
        !
        ! dJ/ds = [ ( d2x/dsdr * dy/ds ) + ( dx/dr * d2y/ds2 ) ] - 
        !         [ ( d2x/ds2 * dy/dr ) + ( dx/ds * d2y/dsdr ) ]
        !
        dJ(j)%ds = ( dx(j)%drds*dy(j)%ds + dx(j)%dr*dy(j)%ds2  ) - &
                   ( dx(j)%ds2 *dy(j)%dr + dx(j)%ds*dy(j)%drds )
        !
        tvol = tvol + jacobian(j)*weights2d(j-p1+1,1)
        !
      end do
      !
      volj(i) = tvol
      !
    else if (n2-n1+1 == 4) then
      !
      ! Metrics for quadrilateral cell
      !
      tvol = zero
      !
      do j = p1,p2
        !
        np = j-p1+1
        !
        l = (np-1)/n_sp1de + 1
        k = np - (l-1)*n_sp1de
        !
        ! These are correct if the reference quad is defined
        ! by the four nodes (0,0) - (1,0) - (1,1) - (0,1)
        !
       !dx(j)%dr = (one-solpts_1d(l))*(xy(1,2)-xy(1,1)) + &
       !                solpts_1d(l) *(xy(1,3)-xy(1,4))
       !dy(j)%dr = (one-solpts_1d(l))*(xy(2,2)-xy(2,1)) + &
       !                solpts_1d(l) *(xy(2,3)-xy(2,4))
       !!
       !dx(j)%ds = (one-solpts_1d(k))*(xy(1,4)-xy(1,1)) + &
       !                solpts_1d(k) *(xy(1,3)-xy(1,2))
       !dy(j)%ds = (one-solpts_1d(k))*(xy(2,4)-xy(2,1)) + &
       !                solpts_1d(k) *(xy(2,3)-xy(2,2))
       !!
       !dx(j)%drds = ( xy(1,1) - xy(1,2) + xy(1,3) - xy(1,4) )
       !dy(j)%drds = ( xy(2,1) - xy(2,2) + xy(2,3) - xy(2,4) )
        !
        ! These are correct if the reference quad is defined
        ! by the four nodes (-1,-1) - (1,-1) - (1,1) - (-1,1)
        !
        dx(j)%dr = c0_25 * ( (solpts_1d(l)-one)*(xy(1,1)-xy(1,2)) + &
                             (solpts_1d(l)+one)*(xy(1,3)-xy(1,4)) )
        dy(j)%dr = c0_25 * ( (solpts_1d(l)-one)*(xy(2,1)-xy(2,2)) + &
                             (solpts_1d(l)+one)*(xy(2,3)-xy(2,4)) )
        !
        dx(j)%ds = c0_25 * ( (solpts_1d(k)-one)*(xy(1,1)-xy(1,4)) + &
                             (solpts_1d(k)+one)*(xy(1,3)-xy(1,2)) )
        dy(j)%ds = c0_25 * ( (solpts_1d(k)-one)*(xy(2,1)-xy(2,4)) + &
                             (solpts_1d(k)+one)*(xy(2,3)-xy(2,2)) )
        !
        dx(j)%drds = c0_25 * ( xy(1,1) - xy(1,2) + xy(1,3) - xy(1,4) )
        dy(j)%drds = c0_25 * ( xy(2,1) - xy(2,2) + xy(2,3) - xy(2,4) )
        !
        dx(j)%dr2 = zero
        dy(j)%dr2 = zero
        !
        dx(j)%ds2 = zero
        dy(j)%ds2 = zero
        !
        ! J = ( dx/dr * dy/ds ) - ( dx/ds * dy/dr )
        !
        jacobian(j) = dx(j)%dr*dy(j)%ds - dx(j)%ds*dy(j)%dr
        !
        ! dJ/dr = [ ( d2x/dr2 * dy/ds ) + ( dx/dr * d2y/drds ) ] - 
        !         [ ( d2x/drds * dy/dr ) + ( dx/ds * d2y/dr2 ) ]
        !
        dJ(j)%dr = ( dx(j)%dr2 *dy(j)%ds + dx(j)%dr*dy(j)%drds ) - &
                   ( dx(j)%drds*dy(j)%dr + dx(j)%ds*dy(j)%dr2  )
        !
        ! dJ/ds = [ ( d2x/dsdr * dy/ds ) + ( dx/dr * d2y/ds2 ) ] - 
        !         [ ( d2x/ds2 * dy/dr ) + ( dx/ds * d2y/dsdr ) ]
        !
        dJ(j)%ds = ( dx(j)%drds*dy(j)%ds + dx(j)%dr*dy(j)%ds2  ) - &
                   ( dx(j)%ds2 *dy(j)%dr + dx(j)%ds*dy(j)%drds )
        !
        tvol = tvol + jacobian(j)*weights2d(j-p1+1,2)
        !
      end do
      !
     !qdtr_jac(:) = zero
     !dxdr(:) = zero
     !dxds(:) = zero
     !dydr(:) = zero
     !dyds(:) = zero
     !do j = 1,n_qp2dq
     !  l = (j-1)/n_qp1de + 1
     !  k = j - (l-1)*n_qp1de
     !  dxdr(j) = c0_25 * ( (qdtpts_1d(l)-one)*(xy(1,1)-xy(1,2)) + &
     !                      (qdtpts_1d(l)+one)*(xy(1,3)-xy(1,4)) )
     !  dydr(j) = c0_25 * ( (qdtpts_1d(l)-one)*(xy(2,1)-xy(2,2)) + &
     !                      (qdtpts_1d(l)+one)*(xy(2,3)-xy(2,4)) )
     !  dxds(j) = c0_25 * ( (qdtpts_1d(k)-one)*(xy(1,1)-xy(1,4)) + &
     !                      (qdtpts_1d(k)+one)*(xy(1,3)-xy(1,2)) )
     !  dyds(j) = c0_25 * ( (qdtpts_1d(k)-one)*(xy(2,1)-xy(2,4)) + &
     !                      (qdtpts_1d(k)+one)*(xy(2,3)-xy(2,2)) )
     !  qdtr_jac(j) = dxdr(j)*dyds(j) - dxds(j)*dydr(j)
     !end do
     !!
     !!------------------------------------------------------------------
     !write (112,*) "Cell = ",i
     !write (112,*) " Jacobian at solution points"
     !write (112,*) jacobian(p1:p2)
     !write (112,*) "sum = ",sum(jacobian(p1:p2)*weights2d(:,2))
     !write (112,*)
     !!
     !write (113,*) "Cell = ",i
     !write (113,*) " Jacobian at quadrature points"
     !write (113,*) matmul( interp_s2q , jacobian(p1:p2) )
     !write (113,*) "sum = ",sum(matmul( interp_s2q , jacobian(p1:p2) )*qdtr_wts(:))
     !write (113,*)
     !!
     !write (114,*) "Cell = ",i
     !write (114,*) " Jacobian at quadrature points"
     !write (114,*) qdtr_jac(:)
     !write (114,*) "sum = ",sum(qdtr_jac(:)*qdtr_wts(:))
     !write (114,*)
     !!------------------------------------------------------------------
     !write (122,*) "Cell = ",i
     !write (122,*) " dxdr at solution points"
     !write (122,*) dx(p1:p2)%dr
     !write (122,*)
     !!
     !write (123,*) "Cell = ",i
     !write (123,*) " dxdr at quadrature points"
     !write (123,*) matmul( interp_s2q , dx(p1:p2)%dr )
     !write (123,*)
     !!
     !write (124,*) "Cell = ",i
     !write (124,*) " dxdr at quadrature points"
     !write (124,*) dxdr(:)
     !write (124,*)
     !!------------------------------------------------------------------
     !write (132,*) "Cell = ",i
     !write (132,*) " dxds at solution points"
     !write (132,*) dx(p1:p2)%ds
     !write (132,*)
     !!
     !write (133,*) "Cell = ",i
     !write (133,*) " dxds at quadrature points"
     !write (133,*) matmul( interp_s2q , dx(p1:p2)%ds )
     !write (133,*)
     !!
     !write (134,*) "Cell = ",i
     !write (134,*) " dxds at quadrature points"
     !write (134,*) dxds(:)
     !write (134,*)
     !!------------------------------------------------------------------
     !write (142,*) "Cell = ",i
     !write (142,*) " dydr at solution points"
     !write (142,*) dy(p1:p2)%dr
     !write (142,*)
     !!
     !write (143,*) "Cell = ",i
     !write (143,*) " dydr at quadrature points"
     !write (143,*) matmul( interp_s2q , dy(p1:p2)%dr )
     !write (143,*)
     !!
     !write (144,*) "Cell = ",i
     !write (144,*) " dydr at quadrature points"
     !write (144,*) dydr(:)
     !write (144,*)
     !!------------------------------------------------------------------
     !write (152,*) "Cell = ",i
     !write (152,*) " dyds at solution points"
     !write (152,*) dy(p1:p2)%ds
     !write (152,*)
     !!
     !write (153,*) "Cell = ",i
     !write (153,*) " dyds at quadrature points"
     !write (153,*) matmul( interp_s2q , dy(p1:p2)%ds )
     !write (153,*)
     !!
     !write (154,*) "Cell = ",i
     !write (154,*) " dyds at quadrature points"
     !write (154,*) dyds(:)
     !write (154,*)
      !------------------------------------------------------------------
      !
      volj(i) = tvol
      !
    end if
    !
    if (i <= ncell) local_volume = local_volume + volj(i)
    !
  end do metrics_loop
 !write (*,*) mypnum,local_volume,sum(volj(1:ncell))
  !
  ! Get the face normal at the solution points on each face
  !
  do nf = 1,nface
    !
    ! Nodes of the current face
    !
    n1 = nodes_on_face(1,nf)
    n2 = nodes_on_face(2,nf)
    !
    ! Outward normal for all solution points on this face
    ! (assuming the face is linear between the two nodes)
    !
    ! NOTE: 1. This is correct for the reference triangle defined by
    !          the three nodes [(0,0)-(1,0)-(0,1)] and for the reference
    !          quad defined by the four nodes [(0,0)-(1,0)-(1,1)-(0,1)].
    !       2. This is incorrect (it needs to be multiplied by 0.5) for
    !          the reference triangle defined by the three nodes
    !          [(-1,-1)-(1,-1)-(-1,1)] and for the reference quad defined
    !          by the four nodes [(-1,-1)-(1,-1)-(1,1)-(-1,1)].
    !       3. Since the reference triangle is defined by (1) and the
    !          reference quad is defined by (2), the variable "cl" defined
    !          in the header of the module flux_mod is used to adjust the
    !          face normals between triangles and quads.
    !
    geofa(1,:,nf) = xyz_nodes(2,n2) - xyz_nodes(2,n1)
    geofa(2,:,nf) = xyz_nodes(1,n1) - xyz_nodes(1,n2)
    !
  end do
  !
  ! Reduce the volumes of each partition to the root
  ! processor to determine the volume of the global grid.
  !
  tot_vol = real(local_volume,r8)
  loc_vol = real(local_volume,r8)
 !total_volume = local_volume
  !
  if (ncpu > 1) then
    call mpi_allreduce(loc_vol,tot_vol,1,mpi_flttyp, &
                       mpi_sum,mpi_comm_world,mpierr)
  end if
  !
  total_volume = real(tot_vol,kind=wp)
  local_volume = real(loc_vol,kind=wp)
  !
 !tot_vol = 1.0e-7_wp
 !do
 !  val = zero
 !  do i = 1,ncell
 !    n1 = pts_in_cell_ptr(i)+1
 !    n2 = pts_in_cell_ptr(i+1)
 !    np = n2-n1+1
 !    p1= merge(1,2,np==n_sp2dt)
 !    do k = n1,n2
 !      val = val + tot_vol*tot_vol*jacobian(k)*weights2d(k-n1+1,p1)
 !    end do
 !  end do
 !  val = val/total_volume
 !  if (log10(val) > -3.0_wp) then
 !    exit
 !  else
 !    tot_vol = tot_vol*1.1_wp
 !  end if
 !end do
 !write (iout,88) tot_vol,val,log10(val)
 !!
 !88 format (/,3es14.6,/)
 !!
 !call stop_recon(stop_mpi,"end of get_metrics")
 !if (mypnum == 0) write (*,*) "Total volume = ",total_volume
  call mpi_barrier(mpi_comm_world,mpierr)
  !
end subroutine get_metrics
!
!###############################################################################
!
include "Functions/np2nfp.f90"
include "Functions/last.f90"
include "Functions/chop.f90"
!
!###############################################################################
!
end module metrics_mod
