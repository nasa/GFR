module mappings_mod
  !
  !.. Global Use Statements ..
  use module_kind_types
  !
  implicit none
  !
  private
  !
  !
  ! ##################################
  ! ###  Module Public Procedures  ###
  ! ##################################
  !
  public :: get_solution_points
  !
  !
  ! ###################################
  ! ###  Module Generic Interfaces  ###
  ! ###################################
  !
  interface get_quad_weights_at_solpts
    module procedure get_quad_weights_at_solpts_DP, &
                     get_quad_weights_at_solpts_QP
  end interface get_quad_weights_at_solpts
  !
contains
!
!###############################################################################
!
subroutine get_solution_points()
  !
  use geovar,         only : nr,nfbnd,ncell,nnode,n_totcel,n_totnod
  use geovar,         only : n_solpts,n_totpts,pts_in_cell_ptr
  use geovar,         only : nodes_of_cell_ptr,nodes_of_cell
  use geovar,         only : bface,xyz_nodes,xyz
  use geovar,         only : gauss_1d,gauss_wts_1d,gauss_wts_2d
  use geovar,         only : lobatto_1d,lobatto_wts_1d,lobatto_wts_2d
  use geovar,         only : solpts_1d,solwts_1d,flxpts_1d,flxwts_1d
  use geovar,         only : rs
  use geovar,         only : cell_geom,cell_order
  !
  use qpvar,          only : gauss_1d_qp,gauss_wts_1d_qp,gauss_wts_2d_qp
  use qpvar,          only : lobatto_1d_qp,lobatto_wts_1d_qp,lobatto_wts_2d_qp
  use qpvar,          only : vp_qp
  !
  use ovar,           only : iop_solpts,iop_flxpts,iop_tripts,error_loc
  !
  use polynomial_mod, only : nodes_legendre_gauss
  use polynomial_mod, only : nodes_legendre_gauss_lobatto
  use polynomial_mod, only : weights_legendre_gauss
  use polynomial_mod, only : weights_legendre_gauss_lobatto
  !
  use parallel_mod,   only : collect_global_solpts
  !
  use order_mod,      only : n_order,geom_solpts
  !
  use quadrature_mod, only : std_elem,face_elem
  !
  !.. Local Scalars ..
  integer  :: i,j,n,n1,n2,p1,p2,np,ierr
  integer  :: npts,this_geom,this_order
  integer  :: nsp1de,nsp2dt,nsp2dq,nsphex
  real(wp) :: r,s
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr,1:8) :: cell_xyz
  !
  !.. Local Pointers ..
  real(wp), pointer, contiguous :: edge_pts(:)
  real(wp), pointer, contiguous :: tria_pts(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_solution_points"
  !
continue
  !
  ! Initialize the local pointers
  !
  edge_pts => null()
  tria_pts => null()
  !
  ! We need to reset various grid size parameters in case partitioning
  ! was done for use with multiple processors. If only using with a
  ! single processor, this should just give the values we already have.
  !
  nfbnd = size(bface,dim=2)
  n_totcel = size(nodes_of_cell_ptr) - 1
  ncell = minval(bface(10,:)) - 1
  n_totnod = size(xyz_nodes,dim=2)
  nnode = maxval(nodes_of_cell(1:nodes_of_cell_ptr(ncell+1)))
  !
  nsp1de = geom_solpts(Geom_Edge,n_order)
  nsp2dt = geom_solpts(Geom_Tria,n_order)
  nsp2dq = geom_solpts(Geom_Quad,n_order)
  nsphex = geom_solpts(Geom_Hexa,n_order)
  !
  ! Get the 1D coordinates and weights of various node types
  !
  ! Allocate the working precision arrays
  !
  allocate ( gauss_1d(1:nsp1de) , source=zero , stat=ierr )
  call alloc_error(pname,"gauss_1d",1,ierr)
  allocate ( gauss_wts_1d(1:nsp1de) , source=zero , stat=ierr )
  call alloc_error(pname,"gauss_wts_1d",1,ierr)
  allocate ( gauss_wts_2d(1:nsp2dq) , source=zero , stat=ierr )
  call alloc_error(pname,"gauss_wts_2d",1,ierr)
  !
  allocate ( lobatto_1d(1:nsp1de) , source=zero , stat=ierr )
  call alloc_error(pname,"lobatto_1d",1,ierr)
  allocate ( lobatto_wts_1d(1:nsp1de) , source=zero , stat=ierr )
  call alloc_error(pname,"lobatto_wts_1d",1,ierr)
  allocate ( lobatto_wts_2d(1:nsp2dq) , source=zero , stat=ierr )
  call alloc_error(pname,"lobatto_wts_2d",1,ierr)
  !
  allocate ( solpts_1d(1:nsp1de) , source=zero , stat=ierr )
  call alloc_error(pname,"solpts_1d",1,ierr)
  allocate ( solwts_1d(1:nsp1de) , source=zero , stat=ierr )
  call alloc_error(pname,"solwts_1d",1,ierr)
  !
  allocate ( flxpts_1d(1:nsp1de) , source=zero , stat=ierr )
  call alloc_error(pname,"flxpts_1d",1,ierr)
  allocate ( flxwts_1d(1:nsp1de) , source=zero , stat=ierr )
  call alloc_error(pname,"flxwts_1d",1,ierr)
  !
  allocate ( rs(1:nr,1:nsp2dt) , source=zero , stat=ierr )
  call alloc_error(pname,"rs",1,ierr)
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
  allocate ( gauss_1d_qp(1:nsp1de) , source=qzero , stat=ierr )
  call alloc_error(pname,"gauss_1d_qp",1,ierr)
  allocate ( gauss_wts_1d_qp(1:nsp1de) , source=qzero , stat=ierr )
  call alloc_error(pname,"gauss_wts_1d_qp",1,ierr)
  allocate ( gauss_wts_2d_qp(1:nsp2dq) , source=qzero , stat=ierr )
  call alloc_error(pname,"gauss_wts_2d_qp",1,ierr)
  !
  allocate ( lobatto_1d_qp(1:nsp1de) , source=qzero , stat=ierr )
  call alloc_error(pname,"lobatto_1d_qp",1,ierr)
  allocate ( lobatto_wts_1d_qp(1:nsp1de) , source=qzero , stat=ierr )
  call alloc_error(pname,"lobatto_wts_1d_qp",1,ierr)
  allocate ( lobatto_wts_2d_qp(1:nsp2dq) , source=qzero , stat=ierr )
  call alloc_error(pname,"lobatto_wts_2d_qp",1,ierr)
  !
  allocate ( vp_qp(1:nsp1de) , source=qzero , stat=ierr )
  call alloc_error(pname,"vp_qp",1,ierr)
  !
  ! Create the 1D arrays for the locations of the Gauss,
  ! and Lobatto solution points and their corresponding weights
  !
  call nodes_legendre_gauss(nsp1de,gauss_1d_qp,vp_qp)
  call weights_legendre_gauss(nsp1de,gauss_1d_qp,vp_qp,gauss_wts_1d_qp)
  !
  call nodes_legendre_gauss_lobatto(n_order,lobatto_1d_qp,vp_qp)
  call weights_legendre_gauss_lobatto(n_order,vp_qp,lobatto_wts_1d_qp)
  !
  ! Create the 2D weights for the Gauss and Lobatto solution points on quads
  !
  gauss_wts_2d_qp = get_quad_weights_at_solpts(gauss_wts_1d_qp)
  lobatto_wts_2d_qp = get_quad_weights_at_solpts(lobatto_wts_1d_qp)
  !
  ! Chop the quad precision arrays
  !
  gauss_1d_qp       = chop( gauss_1d_qp , gauss_1d_qp )
  gauss_wts_1d_qp   = chop( gauss_wts_1d_qp , gauss_wts_1d_qp )
  gauss_wts_2d_qp   = chop( gauss_wts_2d_qp , gauss_wts_2d_qp )
  lobatto_1d_qp     = chop( lobatto_1d_qp , lobatto_1d_qp )
  lobatto_wts_1d_qp = chop( lobatto_wts_1d_qp , lobatto_wts_1d_qp )
  lobatto_wts_2d_qp = chop( lobatto_wts_2d_qp , lobatto_wts_2d_qp )
  !
  ! ##############################################
  ! ##############################################
  ! ##   CHOP THE QUAD PRECISION ARRAYS AND     ##
  ! ##  STORE IN THE WORKING PRECISION ARRAYS.  ##
  ! ##############################################
  ! ##############################################
  !
  gauss_1d       = chop( gauss_1d_qp , gauss_1d )
  gauss_wts_1d   = chop( gauss_wts_1d_qp , gauss_wts_1d )
  gauss_wts_2d   = chop( gauss_wts_2d_qp , gauss_wts_2d )
  lobatto_1d     = chop( lobatto_1d_qp , lobatto_1d )
  lobatto_wts_1d = chop( lobatto_wts_1d_qp , lobatto_wts_1d )
  lobatto_wts_2d = chop( lobatto_wts_2d_qp , lobatto_wts_2d )
  !
  ! Copy the 1D edge and face quadrature information into solpts_1d,
  ! solwts_1d, flxpts_1d, and flxwts_1d. Eventually these will be
  ! replaced all together.
  !
  if (allocated(std_elem(Geom_Edge,n_order)%pts)) then
    solpts_1d = std_elem(Geom_Edge,n_order)%pts(1,:)
  end if
  if (allocated(std_elem(Geom_Edge,n_order)%wts)) then
    solwts_1d = std_elem(Geom_Edge,n_order)%wts
  end if
  if (allocated(face_elem(Geom_Edge,n_order)%pts)) then
    flxpts_1d = face_elem(Geom_Edge,n_order)%pts(1,:)
  end if
  if (allocated(face_elem(Geom_Edge,n_order)%wts)) then
    flxwts_1d = face_elem(Geom_Edge,n_order)%wts
  end if
  !
  ! Copy the 2D triangle solution points into rs. Eventually this
  ! will be replaced all together.
  !
  if (allocated(std_elem(Geom_Tria,n_order)%pts)) then
    rs = std_elem(Geom_Tria,n_order)%pts
  end if
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
  call alloc_error(pname,"pts_in_cell_ptr",1,ierr)
  !
  do n = 1,n_totcel
    pts_in_cell_ptr(n+1) = pts_in_cell_ptr(n) + &
                           geom_solpts( cell_geom(n) , cell_order(n) )
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
  ! Now calculate the coordinates of the solution points within each cell
  !
  allocate ( xyz(1:nr,1:n_totpts) , source=zero , stat=ierr )
  call alloc_error(pname,"xyz",1,ierr)
  !
  !
 !do this_order = lbound(std_elem,dim=2),ubound(std_elem,dim=2)
 !  do this_geom = lbound(std_elem,dim=1),ubound(std_elem,dim=1)
 !    if (allocated(std_elem(this_geom,this_order)%pts)) then
 !      write (iout,3) Geom_Name(this_geom),this_order,"pts", &
 !      trim(merge("IS    ","is NOT",is_contiguous(std_elem(this_geom,this_order)%pts)))
 !    end if
 !    if (allocated(std_elem(this_geom,this_order)%wts)) then
 !      write (iout,3) Geom_Name(this_geom),this_order,"wts", &
 !      trim(merge("IS    ","is NOT",is_contiguous(std_elem(this_geom,this_order)%wts)))
 !    end if
 !  end do
 !end do
  !
  do i = 1,n_totcel
    !
    n1 = nodes_of_cell_ptr(i)+1
    n2 = nodes_of_cell_ptr(i+1)
    !
    p1 = pts_in_cell_ptr(i)+1
    p2 = pts_in_cell_ptr(i+1)
    !
    this_geom = cell_geom(i)
    this_order = cell_order(i)
    !
    do j = n1,n2
      cell_xyz(1:nr,j-n1+1) = xyz_nodes(1:nr,nodes_of_cell(j))
    end do
    !
    select case (this_geom)
      !
      case (Geom_Edge)
        !
        edge_pts => std_elem(Geom_Edge,this_order)%pts(1,:)
        !
        xyz(1:nr,p1:p2) = get_solpts_edge( edge_pts , cell_xyz(1:nr,1:2) )
        !
      case (Geom_Tria)
        !
        tria_pts => std_elem(Geom_Tria,this_order)%pts
        !
        xyz(1:nr,p1:p2) = get_solpts_tria( tria_pts , cell_xyz(1:nr,1:3) )
        !
      case (Geom_Quad)
        !
        edge_pts => std_elem(Geom_Edge,this_order)%pts(1,:)
        !
        xyz(1:nr,p1:p2) = get_solpts_quad( edge_pts , cell_xyz(1:nr,1:4) )
        !
      case (Geom_Tetr)
        !
        write (iout,2) "Tetrahedral"
        call stop_recon(stop_mpi,pname)
        !
        tria_pts => std_elem(Geom_Tria,this_order)%pts
        !
        xyz(1:nr,p1:p2) = get_solpts_tetr( tria_pts , cell_xyz(1:nr,1:4) )
        !
      case (Geom_Pyra)
        !
        write (iout,2) "Pyramid"
        call stop_recon(stop_mpi,pname)
        !
        edge_pts => std_elem(Geom_Edge,this_order)%pts(1,:)
        tria_pts => std_elem(Geom_Tria,this_order)%pts
        !
        xyz(1:nr,p1:p2) = get_solpts_pyra( edge_pts , tria_pts , cell_xyz(1:nr,1:5) )
        !
      case (Geom_Pris)
        !
        write (iout,2) "Prism"
        call stop_recon(stop_mpi,pname)
        !
        edge_pts => std_elem(Geom_Edge,this_order)%pts(1,:)
        tria_pts => std_elem(Geom_Tria,this_order)%pts
        !
        xyz(1:nr,p1:p2) = get_solpts_pris( edge_pts , tria_pts , cell_xyz(1:nr,1:6) )
        !
      case (Geom_Hexa)
        !
        edge_pts => std_elem(Geom_Edge,this_order)%pts(1,:)
        !
        xyz(1:nr,p1:p2) = get_solpts_hexa( edge_pts , cell_xyz(1:nr,1:8) )
        !
      case default
        !
        write (iout,1) i,cell_geom(i)
        call stop_recon(stop_mpi,pname)
        !
    end select
    !
    if (associated(edge_pts)) edge_pts => null()
    if (associated(tria_pts)) tria_pts => null()
    !
  end do
  !
  ! Collect the solution point coordinates for the global grid
  !
  call collect_global_solpts
  !
  ! Format Statements
  !
  1 format (" ERROR: Invalid cell geometry!",/, &
            "        Cell = ",i0,/, &
            "        Geom = ",i0,/)
  2 format ("ERROR: ",a," cells are currently not supported!")
  3 format (" std_elem(",a,",",i0,")%",a,1x,a," contiguous")
  !
end subroutine get_solution_points
!
!###############################################################################
!
pure function get_solpts_edge(xi,xy) result(return_value)
  !
  real(wp), dimension(:),   intent(in) :: xi
  real(wp), dimension(:,:), intent(in) :: xy
  !
  real(wp), dimension(1:size(xy,1),1:size(xi)) :: return_value
  !
  integer  :: i
  real(wp) :: xim1,xip1
  !
continue
  !
  do i = 1,size(xi)
    !
    xim1 = one - xi(i)
    xip1 = one + xi(i)
    !
    return_value(:,i) = half * ( xim1*xy(:,1) + xip1*xy(:,2) )
    !
  end do
  !
end function get_solpts_edge
!
!###############################################################################
!
pure function get_solpts_tria(rs,xy) result(return_value)
  !
  real(wp), dimension(:,:), intent(in) :: rs
  real(wp), dimension(:,:), intent(in) :: xy
  !
  real(wp), dimension(1:size(rs,1),1:size(rs,2)) :: return_value
  !
  integer  :: i
  real(wp) :: r,s,a,b,c
  !
continue
  !
  ! If tria_is_biunit is true:  a=0 ; b=1 ; c=1/2
  ! If tria_is_biunit is false: a=1 ; b=0 ; c=1
  !
  b = deltafun(tria_is_biunit)
  a = one - b
  c = one - half*b
  !
  ! On entry, r and s both need to be defined on the interval [0,1]
  !
  do i = 1,size(rs,2)
    !
    r = rs(1,i)
    s = rs(2,i)
    !
    return_value(:,i) = c * ((a-r-s)*xy(:,1) + (r+b)*xy(:,2) + (s+b)*xy(:,3))
    !
  end do
  !
end function get_solpts_tria
!
!###############################################################################
!
pure function get_solpts_quad(xi,xy) result(return_value)
  !
  real(wp), dimension(:),   intent(in) :: xi
  real(wp), dimension(:,:), intent(in) :: xy
  !
  real(wp), dimension(1:size(xy,1),1:size(xi)**2) :: return_value
  !
  integer  :: i,j,l
  real(wp) :: etam1,etap1,xim1,xip1
  !
continue
  !
  l = 1
  do j = 1,size(xi)
    !
    etam1 = one - xi(j)
    etap1 = one + xi(j)
    !
    do i = 1,size(xi)
      !
      xim1 = one - xi(i)
      xip1 = one + xi(i)
      !
      return_value(:,l) = one4 * ( etam1*xim1*xy(:,1) + &
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
end function get_solpts_quad
!
!###############################################################################
!
pure function get_solpts_tetr(rs,xyz) result(return_value)
  !
  real(wp), dimension(:,:), intent(in) :: rs
  real(wp), dimension(:,:), intent(in) :: xyz
  !
  real(wp), allocatable :: return_value(:,:)
  !
  integer :: ndim,npts,ierr
  integer :: n
  !
continue
  !
  n = size(rs,dim=2)
  !
  ndim = size(xyz,dim=1)
  npts = n*(3+nint(sqrt(real(1+8*n,kind=wp))))/6
  !
  allocate ( return_value(1:ndim,1:npts) , source=zero , stat=ierr )
  !
end function get_solpts_tetr
!
!###############################################################################
!
pure function get_solpts_pyra(xi,rs,xyz) result(return_value)
  !
  real(wp), dimension(:),   intent(in) :: xi
  real(wp), dimension(:,:), intent(in) :: rs
  real(wp), dimension(:,:), intent(in) :: xyz
  !
  real(wp), allocatable :: return_value(:,:)
  !
  integer :: ndim,npts,ierr
  integer :: n
  !
continue
  !
  n = size(xi)-1
  !
  ndim = size(xyz,dim=1)
  npts = (n+1)*(n+2)*(2*n+3)/6
  !
  allocate ( return_value(1:ndim,1:npts) , source=zero , stat=ierr )
  !
end function get_solpts_pyra
!
!###############################################################################
!
pure function get_solpts_pris(xi,rs,xyz) result(return_value)
  !
  real(wp), dimension(:),   intent(in) :: xi
  real(wp), dimension(:,:), intent(in) :: rs
  real(wp), dimension(:,:), intent(in) :: xyz
  !
  real(wp), allocatable :: return_value(:,:)
  !
  integer  :: ndim,npts,ierr
  integer  :: i,j,k
  real(wp) :: r,s,a,b,c,zetam1,zetap1
  !
  real(wp), dimension(1:size(xyz,1)) :: p1
  real(wp), dimension(1:size(xyz,1)) :: p2
  real(wp), dimension(1:size(xyz,1)) :: p3
  !
continue
  !
  ndim = size(xyz,dim=1)
  npts = size(xi)*size(rs,dim=2)
  !
  allocate ( return_value(1:ndim,1:npts) , source=zero , stat=ierr )
  !
  ! If tria_is_biunit is true:  a=0 ; b=1 ; c=1/2
  ! If tria_is_biunit is false: a=1 ; b=0 ; c=1
  !
  b = deltafun(tria_is_biunit)
  a = one - b
  c = one - half*b
  !
  k = 0
  !
  do j = 1,size(xi)
    !
    zetam1 = one - xi(j)
    zetap1 = one + xi(j)
    !
    ! Find the nodes of the triangle for the current xi-eta plane.
    !
    p1 = half * (zetam1*xyz(:,1) + zetap1*xyz(:,4))
    p2 = half * (zetam1*xyz(:,2) + zetap1*xyz(:,5))
    p3 = half * (zetam1*xyz(:,3) + zetap1*xyz(:,6))
    !
    do i = 1,size(rs,dim=2)
      !
      k = k + 1
      !
      r = rs(1,i)
      s = rs(2,i)
      !
      return_value(:,k) = c * ((a-r-s)*p1 + (r+b)*p2 + (s+b)*p3)
      !
    end do
    !
  end do
  !
end function get_solpts_pris
!
!###############################################################################
!
pure function get_solpts_hexa(xi,xyz) result(return_value)
  !
  real(wp), dimension(:),   intent(in) :: xi
  real(wp), dimension(:,:), intent(in) :: xyz
  !
  real(wp), dimension(1:size(xyz,1),1:size(xi)**3) :: return_value
  !
  integer  :: i,j,k,l
  real(wp) :: xim1,xip1,etam1,etap1,zetam1,zetap1
  !
continue
  !
  l = 1
  do k = 1,size(xi)
    !
    zetam1 = one - xi(k)
    zetap1 = one + xi(k)
    !
    do j = 1,size(xi)
      !
      etam1 = one - xi(j)
      etap1 = one + xi(j)
      !
      do i = 1,size(xi)
        !
        xim1 = one - xi(i)
        xip1 = one + xi(i)
        !
        return_value(:,l) = one8 * ( xim1*etam1*zetam1*xyz(:,1) + &
                                     xip1*etam1*zetam1*xyz(:,2) + &
                                     xip1*etap1*zetam1*xyz(:,3) + &
                                     xim1*etap1*zetam1*xyz(:,4) + &
                                     xim1*etam1*zetap1*xyz(:,5) + &
                                     xip1*etam1*zetap1*xyz(:,6) + &
                                     xip1*etap1*zetap1*xyz(:,7) + &
                                     xim1*etap1*zetap1*xyz(:,8) )
        !
        l = l + 1
        !
      end do
      !
    end do
    !
  end do
  !
end function get_solpts_hexa
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
include "Functions/last.f90"
!
!###############################################################################
!
end module mappings_mod
