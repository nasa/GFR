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
 !public :: get_flux_points_2d
  public :: get_grid_connectivity
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
contains
!
!###############################################################################
!
pure function get_solution_points_1d(xi,xnodes) result(return_value)
  !
  real(vbp), dimension(:), intent(in) :: xi
  real(vbp), dimension(:), intent(in) :: xnodes
  !
  real(vbp), dimension(size(xi),size(xnodes)-1) :: return_value
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
subroutine get_flux_points_2d()
  !
  !.. Use Statements ..
  use geovar,   only : n_order,n_sp1de,nr,nface
  use geovar,   only : flxpts_1d,flxwts_1d
  use geovar,   only : xyz_flxpts,xyz_nodes,nodes_on_face
  use ovar,     only : iop_flxpts
  use polynomial_mod, only : nodes_legendre_gauss, &
                             nodes_legendre_gauss_lobatto, &
                             nodes_chebyshev_gauss, &
                             nodes_chebyshev_gauss_lobatto, &
                             weights_legendre_gauss, &
                             weights_legendre_gauss_lobatto, &
                             weights_chebyshev_gauss, &
                             weights_chebyshev_gauss_lobatto
  !
  !.. Local Scalars ..
  integer :: i,n,nf,n1,n2,ierr
  !
  !.. Local Arrays ..
  real(vbp), dimension(1:n_sp1de) :: vp
  real(vbp), dimension(1:nr)   :: pt1
  real(vbp), dimension(1:nr)   :: pt2
  real(vbp), dimension(1:nr)   :: delpt
  !
continue
  !
  ! Allocate the 1D flux arrays
  !
  allocate ( flxpts_1d(1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("get_flux_points_2d","flxpts_1d",1,ierr)
  allocate ( flxwts_1d(1:n_sp1de) , source=zero , stat=ierr )
  call alloc_error("get_flux_points_2d","flxwts_1d",1,ierr)
  !
  ! Create the 1D array for the location of the flux points
  !
  if (iop_flxpts == Legendre_Gauss) then
    !
    call nodes_legendre_gauss(n_sp1de,flxpts_1d,vp)
    call weights_legendre_gauss(n_sp1de,flxpts_1d,vp,flxwts_1d)
    !
  else if (iop_flxpts == Legendre_Gauss_Lobatto) then
    !
    call nodes_legendre_gauss_lobatto(n_order,flxpts_1d,vp)
    call weights_legendre_gauss_lobatto(n_order,vp,flxwts_1d)
    !
  else if (iop_flxpts == Chebyshev_Gauss) then
    !
    call nodes_chebyshev_gauss(n_sp1de,flxpts_1d)
    call weights_chebyshev_gauss(n_sp1de,flxwts_1d)
    write (iout,CG_error) "Flux"
    call stop_recon(stop_mpi,"get_flux_points_2d")
    !
  else if (iop_flxpts == Chebyshev_Gauss_Lobatto) then
    !
    call nodes_chebyshev_gauss_lobatto(n_order,flxpts_1d)
    call weights_chebyshev_gauss_lobatto(n_order,flxwts_1d)
    write (iout,CGL_error) "Flux"
    call stop_recon(stop_mpi,"get_flux_points_2d")
    !
  else
    !
    write (iout,PNT_error) "flux"
    call stop_recon(stop_mpi,"get_flux_points_2d")
    !
  end if
 !!
 !! Now create the coordinates of the flux points on each face
 !!
 !allocate ( xyz_flxpts(1:nr,1:n_sp1de,1:nface) , source=zero , stat=ierr )
 !call alloc_error("get_flux_points_2d","xyz_flxpts",1,ierr)
 !!
 !do nf = 1,nface
 !  !
 !  n1 = nodes_on_face(1,nf)
 !  n2 = nodes_on_face(2,nf)
 !  !
 !  pt1(:) = xyz_nodes(:,n1)
 !  pt2(:) = xyz_nodes(:,n2)
 !  !
 !  delpt(:) = pt2(:) - pt1(:)
 !  !
 !  do n = 1,n_sp1de
 !    xyz_flxpts(:,n,nf) = pt1(:) + delpt(:)*(flxpts_1d(n)+one)/two
 !  end do
 !  !
 !end do
  !
end subroutine get_flux_points_2d
!
!###############################################################################
!
subroutine get_solution_points_2d()
  !
  use geovar
  use ovar,           only : iop_solpts,iop_tripts,error_loc
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
  real(vbp) :: r,s
  !
  !.. Local Arrays ..
  real(vbp), dimension(1:nr,1:4) :: xy
  real(vbp), dimension(1:n_order+1) :: vp
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
  ! Set some of the base constants needed
  !
  n_sp1de = n_order + 1                ! number of 1D solution points
  n_sp2dt = n_sp1de * (n_sp1de+1) / 2  ! number of 2D triangle solution points
  n_sp2dq = n_sp1de * n_sp1de          ! number of 2D quad     solution points
  !
  ! Get the 1D coordinates and weights of the Lobatto nodes
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
  ! Create the 1D arrays for the location of the Gauss, Lobatto,
  ! and Radau solution points and their corresponding weights
  !
  call nodes_legendre_gauss(n_sp1de,gauss_1d,vp)
  call weights_legendre_gauss(n_sp1de,gauss_1d,vp,gauss_wts_1d)
  !
  call nodes_legendre_gauss_lobatto(n_order,lobatto_1d,vp)
  call weights_legendre_gauss_lobatto(n_order,vp,lobatto_wts_1d)
  !
  call nodes_right_radau(n_order,radau_1d)
  call weights_right_radau(n_order,radau_1d,radau_wts_1d)
  !
 !write (*,*)
 !write (*,'(100f10.6)') radau_1d
 !write (*,'(100f10.6)') radau_wts_1d
 !write (*,*)
  !
  ! Create the 2D weights for the Gauss and Lobatto solution points on quads
  !
  gauss_wts_2d = get_quad_weights_at_solpts(gauss_wts_1d)
  lobatto_wts_2d = get_quad_weights_at_solpts(lobatto_wts_1d)
  radau_wts_2d = get_quad_weights_at_solpts(radau_wts_1d)
  !
  ! Now get the 1D solution points for the computation
  ! NOTE: We will save into vp whichever of gauss_1d and lobatto_1d isnt solpts_1d
  !
  if (iop_solpts == Legendre_Gauss) then
    !
    solpts_1d(:) = gauss_1d(:)
    solwts_1d(:) = gauss_wts_1d(:)
    !
  else if (iop_solpts == Legendre_Gauss_Lobatto) then
    !
    solpts_1d(:) = lobatto_1d(:)
    solwts_1d(:) = lobatto_wts_1d(:)
    !
  else if (iop_solpts == Chebyshev_Gauss) then
    !
    call nodes_chebyshev_gauss(n_order+1,solpts_1d)
    call weights_chebyshev_gauss(n_order+1,solwts_1d)
    write (iout,CG_error) "Solution"
    call stop_recon(stop_mpi,"get_solution_points_2d")
    !
  else if (iop_solpts == Chebyshev_Gauss_Lobatto) then
    !
    call nodes_chebyshev_gauss_lobatto(n_order,solpts_1d)
    call weights_chebyshev_gauss_lobatto(n_order,solwts_1d)
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
  ! Get the coordinates of the base triangle in r-s space
  !
  allocate ( rs(1:nr,1:n_sp2dt) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","rs",1,ierr)
  allocate ( weights2d(1:n_sp2dq,1:2) , source=zero , stat=ierr )
  call alloc_error("get_solution_points_2d","weights2d",1,ierr)
  !
  if (iop_tripts == DG_TriPoints) then
    rs(:,:) = TriNodes2D_DG(n_order)
  else if (iop_tripts == HT_TriPoints) then
    rs(:,:) = TriNodes2D_HT(n_order)
  end if
  !
  ! r and s come out defined on the interval [0,1]
  !
  ! Get the triangle and quadrilateral weights at the solution points
  !
  weights2d(1:n_sp2dt,1) = get_triangle_weights_at_solpts(rs)
  !
  if (abs(error_loc) == 1) then
    !
    vp(:) = merge( gauss_1d(:) , radau_1d(:) , error_loc > 0 )
    !
    weights2d(1:n_sp2dq,2) = merge( gauss_wts_2d(1:n_sp2dq) , &
                                    radau_wts_2d(1:n_sp2dq) , &
                                    error_loc > 0 )
    !
  else
    !
    vp(:) = solpts_1d(:)
    !
    weights2d(1:n_sp2dq,2) = get_quad_weights_at_solpts(solwts_1d)
    !
  end if
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
  real(vbp), dimension(:,:), intent(in) :: rs
  real(vbp), dimension(:,:), intent(in) :: xy
  !
  real(vbp), dimension(1:size(rs,1),1:size(rs,2)) :: return_value
  !
  integer   :: i
  real(vbp) :: r,s
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
  real(vbp), dimension(:),   intent(in) :: xi
  real(vbp), dimension(:,:), intent(in) :: xy
  !
  real(vbp), dimension(1:size(xy,1),1:size(xi)**2) :: return_value
  !
  integer   :: i,j,l
  real(vbp) :: etam1,etap1,xim1,xip1
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
pure function get_triangle_weights_at_solpts(rs) result(return_value)
  !
  ! NOTE : This seems to fail for large n_order (>15)
  !
  use polynomial_mod, only : nodes_legendre_gauss, &
                             weights_legendre_gauss
  use startup_mod, only : invert_matrix
  !
  real(vbp), dimension(:,:), intent(in) :: rs
  !
  real(vbp), dimension(1:size(rs,2)) :: return_value
  !
  integer   :: i,j,k,n,m,np
  integer   :: nspts,nfpts
  real(vbp) :: g,h,x,y
  !
  integer, parameter :: nquad = 20
  !
  real(vbp), dimension(1:nquad)                   :: xi
  real(vbp), dimension(1:nquad)                   :: dxi
  real(vbp), dimension(1:nquad)                   :: weights
  real(vbp), dimension(1:size(rs,2),1:size(rs,2)) :: V
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
end function get_triangle_weights_at_solpts
!
!###############################################################################
!
pure function get_quad_weights_at_solpts(weights) result(return_value)
  !
  real(vbp), dimension(:), intent(in) :: weights
  !
  real(vbp), dimension(1:size(weights)**2) :: return_value
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
end function get_quad_weights_at_solpts
!
!###############################################################################
!
!subroutine adjust_for_ghost_cells()
!  !
!  ! Routine to make some adjustments to the data from the grid file
!  ! to account for ghost cells and new nodes associated with them 
!  !
!  !.. Use Statements ..
!  use geovar, only : nr,ncell,nfbnd,nnode,n_totcel,n_totnod,n_totpts, &
!                     bface,nodes_of_cell_ptr,nodes_of_cell,xyz_nodes
!  !
!  !.. Local Scalars ..
!  integer   :: g1,g2,i,i1,i2,n1,n2,nf,ng,ni,np
!  integer   :: gcell,icell,igeom,inode,iside
!  integer   :: pf,pcell,p1,p2,pgeom,pside
!  integer   :: np1,np2,ni1,ni2,ierr
!  real(vbp) :: theta
!  !
!  !.. Local Arrays ..
!  real(vbp), dimension(1:nr)     :: ptf,rnrm,ds
!  real(vbp), dimension(1:nr)     :: ds1,ds2,pnrm,inrm
!  real(vbp), dimension(1:nr,1:4) :: xyp
!  !
!  !.. Local Allocatable Arrays ..
!  integer,   allocatable, dimension(:)   :: tmp_int
!  real(vbp), allocatable, dimension(:,:) :: tmp_real
!  !
!  integer, parameter :: rotation = 0
!  !
!continue
!  !
!  ! #######
!  ! ## 1 ##
!  ! #######
!  !
!  ! Total number of interior plus ghost cells
!  !
!  n_totcel = ncell + nfbnd
!  !
!  ! Reallocate nodes_of_cell_ptr to account for ghost cells
!  !
!  allocate ( tmp_int(1:n_totcel+1) , source=0 , stat=ierr )
!  call alloc_error("get_solution_points_2d","tmp_int",1,ierr)
!  ! 
!  tmp_int(1:size(nodes_of_cell_ptr)) = nodes_of_cell_ptr
!  call move_alloc( from=tmp_int , to=nodes_of_cell_ptr )
!  !
!  ! Set the number of nodes for each ghost cell
!  !
!  do nf = 1,nfbnd
!    !
!    if (bface(1,nf) /= bc_periodic) then
!      !
!      ! Not a periodic boundary, number of nodes is the same as interior
!      !
!      inode = merge(3,4,bface(3,nf)==3)
!      !
!    else
!      !
!      pf = bface(5,nf)
!      pcell = bface(2,pf)
!      !
!      n1 = nodes_of_cell_ptr(pcell)+1
!      n2 = nodes_of_cell_ptr(pcell+1)
!      !
!      inode = merge(3,4,n2-n1+1==3)
!      !
!    end if
!    !
!    nodes_of_cell_ptr(ncell+nf+1) = inode
!    !
!  end do
!  !
!  ! Reshuffle the ghost cell portion of nodes_of_cell_ptr
!  !
!  do i = ncell+2,size(nodes_of_cell_ptr)
!    nodes_of_cell_ptr(i) = nodes_of_cell_ptr(i) + nodes_of_cell_ptr(i-1)
!  end do
!  !
!  ! #######
!  ! ## 2 ##
!  ! #######
!  !
!  ! Reallocate nodes_of_cell to account for ghost cells
!  ! 
!  allocate ( tmp_int(1:last(nodes_of_cell_ptr)) , source=0 , stat=ierr )
!  call alloc_error("get_solution_points_2d","tmp_int",1,ierr)
!  !
!  tmp_int(1:size(nodes_of_cell)) = nodes_of_cell
!  call move_alloc( from=tmp_int , to=nodes_of_cell )
!  !
!  ! Loop over the boundary faces and find the nodes of each ghost cell.
!  ! The first two nodes are on the boundary face while the other nodes
!  ! are newly created and to be defined in step 3. The nodes are in
!  ! counter-clockwise ordering.
!  !
!  inode = nnode
!  !
!  do nf = 1,nfbnd
!    !
!    gcell = nf + ncell
!    !
!    ! Node indices for the ghost cell
!    !
!    g1 = nodes_of_cell_ptr(gcell)+1
!    g2 = nodes_of_cell_ptr(gcell+1)
!    !
!    ! Information about boundary face
!    !
!    icell = bface(2,nf)
!    igeom = merge(3,4,bface(3,nf)==3)
!    iside = bface(4,nf)
!    !
!    ! Starting node index for the interior cell
!    !
!    i1 = nodes_of_cell_ptr(icell)+1
!    !
!    ! Assign nodes on boundary from interior cell
!    !
!    nodes_of_cell(g1+0) = nodes_of_cell(i1+mod(iside,igeom))
!    nodes_of_cell(g1+1) = nodes_of_cell(i1+iside-1)
!    !
!    ! Create new nodes for the remaining ghost cell nodes
!    !
!    do i = g1+2,g2
!      inode = inode + 1
!      nodes_of_cell(i) = inode
!    end do
!    !
!  end do
!  !
!  ! Total number of nodes after the addition of
!  ! the newly created nodes for the ghost cells
!  !
!  n_totnod = inode
!  !
!  ! #######
!  ! ## 3 ##
!  ! #######
!  !
!  ! Reallocate xyz_nodes to account for ghost nodes
!  !
!  allocate ( tmp_real(1:nr,1:n_totnod) , source=zero , stat=ierr )
!  call alloc_error("get_solution_points_2d","tmp_real",1,ierr)
!  !
!  tmp_real(1:nr,1:nnode) = xyz_nodes
!  call move_alloc( from=tmp_real , to=xyz_nodes )
!  !
!  ! Loop over boundary faces to compute the
!  ! coordinates of the newly created ghost nodes
!  !
!  do nf = 1,nfbnd
!    !
!    ! Ghost cell and the corresponding node indices
!    !
!    gcell = nf + ncell
!    !
!    g1 = nodes_of_cell_ptr(gcell)+1
!    g2 = nodes_of_cell_ptr(gcell+1)
!    !
!    ! Interior cell and the corresponding node indices
!    !
!    icell = bface(2,nf)
!    !
!    i1 = nodes_of_cell_ptr(icell)+1
!    i2 = nodes_of_cell_ptr(icell+1)
!    !
!    ! Information about boundary face
!    !
!    igeom = merge(3,4,bface(3,nf)==3)
!    iside = bface(4,nf)
!    !
!    ! Conditional to treat periodic boundary conditions separately
!    !
!    if (bface(1,nf) /= bc_periodic) then
!      !
!      ! This is not a periodic boundary face
!      !
!      ! Nodes on face
!      !
!      n1 = nodes_of_cell(g1)
!      n2 = nodes_of_cell(g1+1)
!      !
!      ! Center of boundary face by averaging the coordinates of the two boundary nodes
!      !
!      ptf(:) = half*(xyz_nodes(:,n1)+xyz_nodes(:,n2))
!      !
!      ! Outward unit normal to boundary face. NOTE, we are using the
!      ! ordering of the boundary nodes for the ghost cell which is
!      ! opposite that of the interior cell so the normal is the negative
!      ! of the regular normal equation. (We are essentially finding the
!      ! inward unit normal of the boundary face within the ghost cell)
!      !
!      rnrm(1) = xyz_nodes(2,n1) - xyz_nodes(2,n2)
!      rnrm(2) = xyz_nodes(1,n2) - xyz_nodes(1,n1)
!      rnrm = unit_vector(rnrm)
!      !
!      ! Coordinates for nodes g1 and g1+1 of the ghost cell have already been defined
!      ! so we need to define nodes g1+2 and g2 (g1+2==g2 for triangular cells).
!      !
!      if (igeom == 3) then
!        !
!        ! Ghost node 3 for triangle
!        !
!        ni = nodes_of_cell( i1 + mod(iside+1,igeom) )
!        ng = nodes_of_cell( g1 + 2 )
!        !
!        ! Get the coordinates of ghost node 3
!        !
!        xyz_nodes(:,ng) = merge(rotate(xyz_nodes(:,n1),xyz_nodes(:,n2),xyz_nodes(:,ni)), &
!                                reflect(xyz_nodes(:,ni),ptf,rnrm), &
!                                rotation == 1)
!        !
!      else if (igeom == 4) then
!        !
!        ! Ghost node 3 for quadrilateral
!        !
!        ng = nodes_of_cell( g1 + 2 )
!        ni = nodes_of_cell( i1 + mod(iside+2-max(0,min(1,rotation)),igeom) )
!        !
!        ! Get the coordinates of ghost node 3
!        !
!        xyz_nodes(:,ng) = merge(rotate(xyz_nodes(:,n1),xyz_nodes(:,n2),xyz_nodes(:,ni)), &
!                                reflect(xyz_nodes(:,ni),ptf,rnrm), &
!                                rotation == 1)
!        !
!        ! Ghost node 4 for quadrilateral
!        !
!        ng = nodes_of_cell( g1 + 3 )
!        ni = nodes_of_cell( i1 + mod(iside+1+max(0,min(1,rotation)),igeom) )
!        !
!        ! Get the coordinates of ghost node 4
!        !
!        xyz_nodes(:,ng) = merge(rotate(xyz_nodes(:,n1),xyz_nodes(:,n2),xyz_nodes(:,ni)), &
!                                reflect(xyz_nodes(:,ni),ptf,rnrm), &
!                                rotation == 1)
!        !
!      end if
!      !
!    else
!      !
!      ! This is a periodic boundary face
!      !
!      ! Boundary face linked to this periodic boundary face
!      !
!      pf = bface(5,nf)
!      !
!      ! Linked cell and the corresponding node indices
!      !
!      pcell = bface(2,pf)
!      !
!      p1 = nodes_of_cell_ptr(pcell)+1
!      p2 = nodes_of_cell_ptr(pcell+1)
!      !
!      ! Information about the linked face
!      !
!      pgeom = merge(3,4,bface(3,pf)==3)
!      pside = bface(4,pf)
!      !
!      ! Node coordinates for periodic cell
!      !
!      xyp(:,1:pgeom) = xyz_nodes(:,nodes_of_cell(p1:p2))
!      !
!      ! Nodes on the boundary face for pcell
!      !
!      np1 = nodes_of_cell( p1 + pside - 1 )
!      np2 = nodes_of_cell( p1 + mod(pside,pgeom) )
!      !
!      ! Nodes on the boundary face for icell
!      !
!      ni1 = nodes_of_cell( i1 + iside - 1)
!      ni2 = nodes_of_cell( i1 + mod(iside,igeom) )
!      !
!      ! Center of boundary face for pcell
!      !
!      ptf(:) = half * (xyz_nodes(:,np1) + xyz_nodes(:,np2))
!      !
!      ! Outward unit normal for boundary face of pcell
!      !
!      pnrm(1) = xyz_nodes(2,np2) - xyz_nodes(2,np1)
!      pnrm(2) = xyz_nodes(1,np1) - xyz_nodes(1,np2)
!      pnrm = unit_vector(pnrm)
!      !
!      ! Outward unit normal for boundary face of icell
!      !
!      inrm(1) = xyz_nodes(2,ni2) - xyz_nodes(2,ni1)
!      inrm(2) = xyz_nodes(1,ni1) - xyz_nodes(1,ni2)
!      inrm = unit_vector(inrm)
!      !
!      ! Get the angle between the two outward unit normals
!      !
!      theta = acos(dot_product(pnrm,inrm))
!      !
!      ! Translation of the pcell nodes so that ptf is at the origin
!      !
!      do i = 1,pgeom
!        xyp(:,i) = xyp(:,i) - ptf(:)
!      end do
!      !
!      ! Rotate xyp about ptf theta+180 degrees
!      !
!      do i = 1,pgeom
!        xyp(:,i) = rotate_2d( xyp(:,i) , theta+pi )
!      end do
!      !
!      ! Translation of the pcell nodes back to original ptf
!      !
!      do i = 1,pgeom
!        xyp(:,i) = xyp(:,i) + ptf(:)
!      end do
!      !
!      ! Find the distances between the two corresponding nodes of
!      ! the two boundary faces. The two boundary faces should now be
!      ! parallel to each other but with opposite outward unit normals
!      ! and therefore the distances should be equal. 
!      !
!      ds1(:) = xyz_nodes(:,ni1) - xyp(:,1+mod(pside,pgeom))
!      ds2(:) = xyz_nodes(:,ni2) - xyp(:,pside)
!      !
!      if (abs(ds1(1)-ds2(1))>eps6 .or. abs(ds1(2)-ds2(2))>eps6) then
!        write (*,1) icell,pcell,abs(ds1-ds2)
!      end if
!      !
!      ! Move the interior nodes of pcell to the ghost nodes of icell
!      !
!      do i = 3,pgeom
!        !
!        ng = nodes_of_cell( g1 + i-1 )
!        np = nodes_of_cell( p1 + mod(pside+i-2,pgeom) )
!        !
!        xyz_nodes(:,ng) = xyz_nodes(:,np) + ds1(:)
!        !
!      end do
!      !
!    end if
!    !
!  end do
!  !
!  1 format ("ICELL = ",i0,"  PCELL = ",i0,"   (delta x,delta y) = ",2es24.11)
!  !
!end subroutine adjust_for_ghost_cells
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
  ! Now that we have the face connectivity, get the flux points on each face
  !
  call get_flux_points_2d
  !
  ! Get the solution points on each face
  !
  call find_points_on_each_face
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
  call get_metrics
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
  !       cells_on_face(5,:) = (n_sp1de+1,2*n_sp1de+1,1)
  !   for cells_on_face(4,:) = (    1    ,    2    ,   3   )
  !       cells_on_face(6,:) = (2*n_sp1de,3*n_sp1de,n_sp1de)
  ! For Quadrilaterals : 
  !   for cells_on_face(3,:) = (    1    ,     2     ,     3     ,4)
  !       cells_on_face(5,:) = (n_sp1de+1,2*n_sp1de+1,3*n_sp1de+1,1)
  !   for cells_on_face(4,:) = (    1    ,    2    ,    3    ,   4   )
  !       cells_on_face(6,:) = (2*n_sp1de,3*n_sp1de,4*n_sp1de,n_sp1de)
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
  real(vbp) :: x1,y1,x2,y2,dl,unx,uny,x13,y13,dot
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
  real(vbp), dimension(1:2) :: dxy,xyfp,xykp
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
  real(vbp) :: theta
  !
  !.. Local Arrays ..
  real(vbp), dimension(1:nr) :: rgfc,rnrmg
  real(vbp), dimension(1:nr) :: rpfc,rnrmp
  real(vbp), dimension(1:nr,1:n_sp2dq) :: xyg,xyp
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
    theta = acos(dot_product(rnrmg,rnrmp))
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
  use geovar
  !
  !.. Local Scalars ..
  integer :: j1,j2,j3,j4,n,m,sk,ierr
  !
  !.. Local Arrays ..
  integer, dimension(1:n_sp1de,1:n_sp1de) :: counter
  !
continue
  !
  ! Allocate tria_connectivity and quad_connectivity
  !
  allocate ( tria_connectivity(1:4,1:n_sp2dq) , source=0 , stat=ierr )
  call alloc_error("get_unit_cell_subconnectivity","tria_connectivity",1,ierr)
  allocate ( quad_connectivity(1:4,1:n_sp2dq) , source=0 , stat=ierr )
  call alloc_error("get_unit_cell_subconnectivity","quad_connectivity",1,ierr)
  counter(:,:) = 0
  !
  ! Get the subconnectivity for the unit triangle
  !
  sk = 0
  do n = 1,n_order+1
    do m = 1,n_order+2-n
      sk = sk + 1
      counter(n,m) = sk
    end do
  end do
  !
  sk = 0
  do n = 1,n_order+1
    do m = 1,n_order+1-n
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
  do n = 1,n_order
    do m = 1,n_order
      !
      sk = sk + 1
      quad_connectivity(1,sk) = (n_order+1)*(n-1) + m
      quad_connectivity(2,sk) = (n_order+1)*(n-1) + m + 1
      quad_connectivity(3,sk) = (n_order+1)*(n  ) + m + 1
      quad_connectivity(4,sk) = (n_order+1)*(n  ) + m
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
  real(vbp) :: tvol
  real(r8)  :: loc_vol,tot_vol
  !
  !.. Local Arrays ..
  real(vbp), dimension(1:2,1:4) :: xy
  !
continue
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
  total_volume = real(tot_vol,vbp)
  local_volume = real(loc_vol,vbp)
  !
 !if (mypnum == 0) write (*,*) "Total volume = ",total_volume
  call mpi_barrier(mpi_comm_world,mpierr)
  !
end subroutine get_metrics
!
!###############################################################################
!
pure function GeometricFactors2D(xy,Drs) result(return_value)
  !
  real(vbp), dimension(:,:,:), intent(in) :: xy
  real(vbp), dimension(:,:,:), intent(in) :: Drs
  !
  real(vbp), dimension(1:2,1:2,1:size(xy,2),1:size(xy,3)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,j
  real(vbp) :: xr,xs,yr,ys
  real(vbp) :: jacobian
  !
continue
  !
  do j = 1,size(xy,3)
    do i = 1,size(xy,2)
      !
      xr = dot_product(Drs(1,i,:),xy(1,:,j)) ! dx/dr
      xs = dot_product(Drs(2,i,:),xy(1,:,j)) ! dx/ds
      yr = dot_product(Drs(1,i,:),xy(2,:,j)) ! dy/dr
      ys = dot_product(Drs(2,i,:),xy(2,:,j)) ! dy/ds
      !
      jacobian = xr*ys - xs*yr
      !
      return_value(1,1,i,j) =  ys/jacobian ! dr/dx
      return_value(2,1,i,j) = -yr/jacobian ! ds/dx
      return_value(1,2,i,j) = -xs/jacobian ! dr/dy
      return_value(2,2,i,j) =  xr/jacobian ! ds/dy
      !
    end do
  end do
  !
end function GeometricFactors2D
!
!###############################################################################
!
pure function FaceNormals2D(xy,Drs) result(return_value)
  !
  real(vbp), dimension(:,:,:), intent(in) :: xy
  real(vbp), dimension(:,:,:), intent(in) :: Drs
  !
  real(vbp), dimension(1:2,1:np2nfp(size(xy,2)),1:size(xy,3)) :: return_value
  !
  !.. Local Scalars ..
  integer   :: i,j
  real(vbp) :: xr,xs,yr,ys
  real(vbp) :: jacobian
  !
continue
  !
  do j = 1,size(xy,3)
    do i = 1,size(xy,2)
      !
      xr = dot_product(Drs(i,:,1),xy(1,:,j)) ! dx/dr
      xs = dot_product(Drs(i,:,2),xy(1,:,j)) ! dx/ds
      yr = dot_product(Drs(i,:,1),xy(2,:,j)) ! dy/dr
      ys = dot_product(Drs(i,:,2),xy(2,:,j)) ! dy/ds
      !
      jacobian = xr*ys - xs*yr
      !
      return_value(1,i,j) = -xs
      return_value(2,i,j) =  xr
      !
    end do
  end do
  !
end function FaceNormals2D
!
!###############################################################################
!
!
!###############################################################################
!
!
!###############################################################################
!
!
!###############################################################################
!
include "Functions/np2nfp.f90"
include "Functions/norms.f90"
include "Functions/last.f90"
include "Functions/rotate.f90"
include "Functions/reflect.f90"
include "Functions/rotate_2d.f90"
!
!###############################################################################
!
end module metrics_mod
