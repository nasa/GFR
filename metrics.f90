module metrics_mod
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
  public :: get_metrics
  public :: metrics_alloc_status
  !
  type, public :: face_nrm
    real(wp), allocatable :: v(:,:)
  end type face_nrm
  !
  ! geofa : normal of the face at each of the solution points on the face
  !
  type(face_nrm), public, save, allocatable :: geofa(:)
  !
  ! volj : Volume of the grid cells evaluated by integrating the jacobian at the
  !        solution points using the quadrature rule
  !
  real(wp), public, save, allocatable, dimension(:) :: volj
  !
  ! jacobian : jacobian of the solution points
  !
  real(wp), public, save, allocatable, dimension(:) :: jacobian
  !
  ! dJ : derivative of the jacobian at the solution points
  !
  real(wp), public, save, allocatable, dimension(:,:) :: dJ
  !
  ! metrics : first derivative metrics terms at the solution points
  !     For solution point 'n', metrics(l,i,n,1) gives the partial derivative
  !                             d (x_l)
  !                             -------
  !                             d(xi_i)
  !                    whereas, metrics(l,i,n,2) gives the partial derivative
  !                             d(xi_l)
  !                             ------- * jacobian(n)
  !                             d(x_i)
  !     where ( x_1, x_2, x_3) = (x,y,z) -> the physical coordinates
  !     where (xi_1,xi_2,xi_3) = (xi,eta,zeta) or (r,s,t) -> the computational
  !                           coordinates within the standard elemental region
  !           NOTE: the first index l gives the index of the numerator and
  !                 the second index i gives the index of the denominator
  !           NOTE: For 2D z and zeta (or t) do not exist
  !           NOTE: The jacobian(n) term is implicitly included in the
  !                 metrics(:,:,:,2) part of the array because multiplication
  !                 by the jacobian is desired where this metric term is
  !                 predominantly used.
  !
  real(wp), public, save, allocatable, dimension(:,:,:,:) :: metrics
  !
  ! metrics2: second derivative metrics terms at the solution points
  !     For solution point 'n', metrics2(l,i,j,n) gives the partial derivative 
  !                                  2
  !                                 d (x_l)
  !                             ---------------
  !                             d(xi_i) d(xi_j)
  !     where ( x_1, x_2, x_3) = (x,y,z) -> the physical coordinates
  !     where (xi_1,xi_2,xi_3) = (xi,eta,zeta) or (r,s,t) -> the computational
  !                           coordinates within the standard elemental region
  !           NOTE: For 2D z and zeta (or t) do not exist
  !
  real(wp), public, save, allocatable, dimension(:,:,:,:) :: metrics2
  !
  ! metrics2: second derivative metrics terms at the solution points
  !     For solution point 'n', metrics3(l,i,j,k,n) gives the partial derivative 
  !                                      3
  !                                     d (x_l)
  !                             -----------------------
  !                             d(xi_i) d(xi_j) d(xi_k)
  !     where ( x_1, x_2, x_3) = (x,y,z) -> the physical coordinates
  !     where (xi_1,xi_2,xi_3) = (xi,eta,zeta) or (r,s,t) -> the computational
  !           NOTE: For 2D z and zeta (or t) do not exist
  !
  real(wp), public, save, allocatable, dimension(:,:,:,:,:) :: metrics3
  !
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! #####  ||||||||||||||||   ATTEMPT AT METRICS USING    |||||||||||||||  #####
  ! #####  VVVVVVVVVVVVVVVV  PARAMETERIZED DERIVED TYPES  VVVVVVVVVVVVVVV  #####
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  !
  ! NOTE: THERE IS A BUG IN INTEL FORTRAN 15.0.0 and 15.0.1 THAT PREVENTS THE
  !       USE OF LENGTH TYPE PARAMETERS IN DERIVED TYPES. CHECK BACK IN FEBRUARY
  !       2015 TO SEE IF THIS IS FIXED IN 15.0.2.
  !
  type, public :: sp_metrics
  end type sp_metrics
  !
  ! SP_METRICS_1: Parameterized derived type containing all the first
  !               derivative metric terms at the solution points
  !      COMPONENTS
  !           n : length parameter
  !           d : first derivative metrics
  !                d(l,i) gives the partial derivative 
  !                             d (x_l)
  !                             -------
  !                             d(xi_i)
  !     where ( x_1, x_2, x_3) = (x,y,z) -> the physical coordinates
  !     where (xi_1,xi_2,xi_3) = (xi,eta,zeta) or (r,s,t) -> the computational
  !                           coordinates within the standard elemental region
  !
#ifdef DEBUG_ON
  type, public, extends(sp_metrics) :: sp_metrics_1_2D
    real(wp)     :: d(2,2)
  end type sp_metrics_1_2D
  type, public, extends(sp_metrics) :: sp_metrics_1_3D
    real(wp)     :: d(3,3)
  end type sp_metrics_1_3D
#else
  type, public, extends(sp_metrics) :: sp_metrics_1
 !type, public, extends(sp_metrics) :: sp_metrics_1(n)
   !integer, len :: n = 2
   !real(wp)     :: d(n,n)
    real(wp), allocatable :: d(:,:)
  end type sp_metrics_1
#endif
  !
  !
  ! SP_METRICS_2: Parameterized derived type that extends the type SP_METRICS_1
  !               to also contain all the second derivative metrics terms at
  !               the solution points. It also contains all the derived type
  !               components contained within SP_METRICS_1.
  !      COMPONENTS (all the components of SP_METRICS_1 + ...)
  !           d2 : second derivative metrics
  !                d2(l,i,j) gives the partial derivative 
  !                                  2
  !                                 d (x_l)
  !                             ---------------
  !                             d(xi_i) d(xi_j)
  !     where ( x_1, x_2, x_3) = (x,y,z) -> the physical coordinates
  !     where (xi_1,xi_2,xi_3) = (xi,eta,zeta) or (r,s,t) -> the computational
  !                           coordinates within the standard elemental region
  !
#ifdef DEBUG_ON
  type, public, extends(sp_metrics_1_2D) :: sp_metrics_2_2D
    real(wp) :: d2(2,2,2)
  end type sp_metrics_2_2D
  type, public, extends(sp_metrics_1_3D) :: sp_metrics_2_3D
    real(wp) :: d2(3,3,3)
  end type sp_metrics_2_3D
#else
  type, public, extends(sp_metrics_1) :: sp_metrics_2
   !real(wp) :: d2(n,n,n)
    real(wp), allocatable :: d2(:,:,:)
  end type sp_metrics_2
#endif
  !
  ! SP_METRICS_3: Parameterized derived type that extends the type SP_METRICS_2
  !               to also contain all the third derivative metrics terms at
  !               the solution points. It also contains all the derived type
  !               components contained within SP_METRICS_1 and SP_METRICS_2.
  !      COMPONENTS (all the components of SP_METRICS_1 and SP_METRICS_2 + ...)
  !           d3 : third derivative metrics
  !                d3(l,i,j,k) gives the partial derivative 
  !                                      3
  !                                     d (x_l)
  !                             -----------------------
  !                             d(xi_i) d(xi_j) d(xi_k)
  !     where ( x_1, x_2, x_3) = (x,y,z) -> the physical coordinates
  !     where (xi_1,xi_2,xi_3) = (xi,eta,zeta) or (r,s,t) -> the computational
  !                           coordinates within the standard elemental region
  !
#ifdef DEBUG_ON
  type, public, extends(sp_metrics_2_3D) :: sp_metrics_3_3D
    real(wp) :: d3(3,3,3,3)
  end type sp_metrics_3_3D
#else
  type, public, extends(sp_metrics_2) :: sp_metrics_3
   !real(wp) :: d3(n,n,n,n)
    real(wp), allocatable :: d3(:,:,:,:)
  end type sp_metrics_3
#endif
  !
  ! METRICS : This array contains the grid metrics at each solution point
  !           The type declaration of class(sp_metrics_1(:)) is used to specify
  !           that metrics is an array whose type is in the class of types
  !           consisting of sp_metrics_1 and all of its extensions (i.e.,
  !           sp_metrics_2 and sp_metrics_3). This way we can allocate it to the
  !           correct type and size to prevent wasted memory when the second or
  !           third derivative metrics arent needed. The (:) after sp_metrics_1
  !           means that the length parameter, n, in the type definition is a
  !           deferred parameter that is to be specified at allocation of the
  !           array. If no value is passed to the length parameter, n, at
  !           allocation, it defaults to the value 2 given in the type definition.
  !
  class(sp_metrics), public, allocatable, dimension(:) :: pdt_metrics
  !
  class(sp_metrics), public, save, allocatable, dimension(:) :: pdt_met_dJ
  !
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! #####  AAAAAAAAAAAAAAAA   ATTEMPT AT METRICS USING    AAAAAAAAAAAAAAA  #####
  ! #####  ||||||||||||||||  PARAMETERIZED DERIVED TYPES  |||||||||||||||  #####
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  !
  !
  ! OLDER METRICS TYPE DEFINITIONS
  !
  !
  type, public :: derivative
    real(wp) :: dr
    real(wp) :: ds
    real(wp) :: dt
  end type derivative
  !
  type, public, extends(derivative) :: second_derivative
    real(wp) :: dr2
    real(wp) :: ds2
    real(wp) :: dt2
    real(wp) :: drds
    real(wp) :: drdt
    real(wp) :: dsdr
    real(wp) :: dsdt
    real(wp) :: dtdr
    real(wp) :: dtds
  end type second_derivative
  !
  type, public, extends(second_derivative) :: third_derivative
    real(wp) :: dr3,dr2ds,dr2dt
    real(wp) :: ds3,ds2dr,ds2dt
    real(wp) :: dt3,dt2dr,dt2ds
    real(wp) :: drdsdr,drdsds,drdsdt
    real(wp) :: drdtdr,drdtds,drdtdt
    real(wp) :: dsdrdr,dsdrds,dsdrdt
    real(wp) :: dsdtdr,dsdtds,dsdtdt
    real(wp) :: dtdrdr,dtdrds,dtdrdt
    real(wp) :: dtdsdr,dtdsds,dtdsdt
  end type third_derivative
  !
  ! dx : metrics of the x-direction at the solution points
  ! dy : metrics of the y-direction at the solution points
  ! dz : metrics of the z-direction at the solution points
  !
  type(second_derivative), public, save, allocatable, dimension(:) :: dx
  type(second_derivative), public, save, allocatable, dimension(:) :: dy
  type(second_derivative), public, save, allocatable, dimension(:) :: dz
  !
  ! dJ : derivative of the jacobian at the solution points
  !
  type(derivative), public, save, allocatable, dimension(:) :: dt_dJ
  !
  ! spmat : solution points matrix
  ! dximat : dx matrix
  ! detmat : dy matrix
  ! spcemat : curved edge CE (at bottom), solution points
  ! spceximat : CE metrics d_xi
  ! spceetmat : CE metrics d_eta
  ! spcexi2mat : CE metrics d_xi2
  ! spceet2mat : CE metrics d_eta2
  ! spcexietmat : CE metrics d_xieta
  !
  ! lagrange : gives the Lagrange polynomial coefficients to project
  !            variables/fluxes from the solution points in the cell
  !            interior to the flux points on the cell faces
  !
 !real(wp), public, save, allocatable, dimension(:,:) :: spmat       !(3,ksp2d)
 !real(wp), public, save, allocatable, dimension(:,:) :: spcemat     !(ksp1d+1,ksp2d)
 !real(wp), public, save, allocatable, dimension(:,:) :: spceximat   !(ksp1d+1,ksp2d)
 !real(wp), public, save, allocatable, dimension(:,:) :: spceetmat   !(ksp1d+1,ksp2d)
 !real(wp), public, save, allocatable, dimension(:,:) :: spcexi2mat  !(ksp1d+1,ksp2d)
 !real(wp), public, save, allocatable, dimension(:,:) :: spceet2mat  !(ksp1d+1,ksp2d)
 !real(wp), public, save, allocatable, dimension(:,:) :: spcexietmat !(ksp1d+1,ksp2d)
  !
  ! wsp  : quadrature weights
  ! xiSP : xi coordinates of the edge solution points
  !        !! xisp and xisp1d seemed to be duplicates
  ! xiCE : 
  !
  ! Hesthaven and Warburton equivalencies
  ! wsp = weights
  !
 !real(wp), public, save, allocatable, dimension(:) :: wsp  !(ksp2d)
 !real(wp), public, save, allocatable, dimension(:) :: xiSP !(ksp1d)
 !real(wp), public, save, allocatable, dimension(:) :: xiCE !(ksp1d)
  !
  ! XYSP    : x,y coordinates of the cell solution points
  ! XYSPxi  : d_(x,y)/d_xi  of the cell solution points (i.e. d_x/d_xi and d_y/d_xi)
  ! XYSPet  : d_(x,y)/d_eta of the cell solution points
  ! XYeSP   : x,y coordinates of the edge solution points !! not used
  ! XYeSPxi : d_(x,y)/d_xi of the edge solution points 
  !
  ! Hesthaven and Warburton equivalencies
  ! XYSP(1,:,:)   = x
  ! XYSP(2,:,:)   = y
  ! XYSPxi(1,:,:) = xr = matmul(Dr,x) -> rx =  ys/J
  ! XYSPxi(2,:,:) = yr = matmul(Dr,y) -> ry = -xs/J
  ! XYSPet(1,:,:) = xs = matmul(Ds,x) -> sx = -yr/J
  ! XYSPet(2,:,:) = ys = matmul(Ds,y) -> sy =  xr/J
  !
 !real(wp), public, save, allocatable, dimension(:,:,:) :: XYSP    !(2,ksp2d,-nfbnd:ncell)
 !real(wp), public, save, allocatable, dimension(:,:,:) :: XYSPxi  !(2,ksp2d,-nfbnd:ncell)
 !real(wp), public, save, allocatable, dimension(:,:,:) :: XYSPet  !(2,ksp2d,-nfbnd:ncell)
 !real(wp), public, save, allocatable, dimension(:,:,:) :: XYeSP   !(2,ksp1d,1:nface)
 !real(wp), public, save, allocatable, dimension(:,:,:) :: XYeSPxi !(2,ksp1d,1:nface)
  !
  ! DxyDrs : d_(x,y)/d_(r,s) of the cell solution points
  !
 !real(wp), public, save, allocatable, dimension(:,:,:) :: DxyDrs
  !
  ! rjaSP   : jacobian of the solution points
  ! rjaSPxi : d_jacobian/d_xi of the solution points
  ! rjaSPet : d_jacobian/d_eta of the solution points
  !
  ! Hesthaven and Warburthon equivalencies
  ! rjaSP = J = -xs*yr + xr*ys
  !
 !real(wp), public, save, allocatable, dimension(:,:) :: rjaSP   !(ksp2d,-nfbnd:ncell)
 !real(wp), public, save, allocatable, dimension(:,:) :: rjaSPxi !(ksp2d,-nfbnd:ncell)
 !real(wp), public, save, allocatable, dimension(:,:) :: rjaSPet !(ksp2d,-nfbnd:ncell)
  !
  !... for boundary cells, if id_be(j) = 0, then the cell has a curved edge
  !
  ! be_cSPind(1,k,j) k = 1,ksp: the solution point index of cell be_ic(j)
  ! be_cSPind(2,k,j) k = 1,ksp: the solution point index of exterior cell
  ! be_cSPind(3,k,j) k = 1,ksp: the solution point index of cell be_pc(j)
  ! e_cSPind(1,k,j) k = 1,ksp: the solution point index for the left cell
  ! e_cSPind(2,k,j) k = 1,ksp: the solution point index for the right cell
 !integer, public, save, allocatable, dimension(:,:,:) :: e_cSPind  !(2,ksp1d,nface)
 !integer, public, save, allocatable, dimension(:,:,:) :: be_cSPind !(3,ksp1d,nfbnd)
  !
contains
!
!###############################################################################
!
subroutine get_metrics()
  !
  !.. Use Statements ..
  use geovar, only : nr,nface,n_totpts,n_totcel,ncell
  use geovar, only : nodes_of_cell_ptr,nodes_of_cell
  use geovar, only : nodes_on_face,pts_in_cell_ptr
  use geovar, only : xyz_nodes,cell_geom,cell_order
  use geovar, only : face
  !
  use quadrature_mod, only : std_elem
  !
  use ovar, only : local_volume,total_volume
  use ovar, only : using_curved_boundaries
  use ovar, only : use_3rd_derivative_metrics
  !
  use order_mod, only : geom_solpts
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,ierr
  integer  :: n,nc,n1,n2,nf,np,p1,p2
  integer  :: n_edge_pts,nfnods,nfpts
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  real(wp) :: tvol,tria_coef,tc2,r,s,t
  real(r8) :: loc_vol,tot_vol
  real(wp) :: val,rm1,rp1,sm1,sp1,tm1,tp1,zm1,zp1
  character(len=200) :: array_name
  character(len=1000) :: error_message
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr,1:8) :: nodes
  !
  !.. Local Pointers ..
  real(wp), pointer, contiguous :: weights(:)
  real(wp), pointer, contiguous :: edge_pts(:)
  real(wp), pointer, contiguous :: rs_pts(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_metrics"
  !
continue
  !
  ! Initialize the local pointers to disassociated
  !
  weights => null()
  edge_pts => null()
  rs_pts => null()
  !
  ! Coefficients for triangle metrics depending on
  ! the definition of the standard triangular element
  !
  tria_coef = one - half*deltafun(tria_is_biunit)
  tc2 = half * tria_coef
  !
  ! Allocate the required metrics related arrays
  !
  allocate ( metrics(1:nr,1:nr,1:n_totpts,1:2) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"metrics",1,ierr,error_message)
  !
  allocate ( jacobian(1:n_totpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"jacobian",1,ierr,error_message)
  !
  allocate ( volj(1:n_totcel) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"volj",1,ierr,error_message)
  !
  ! Allocate the outer part of the geofa array
  !
  allocate ( geofa(1:nface) , stat=ierr )
  call alloc_error(pname,"geofa",1,ierr,error_message)
  !
  ! Allocate the data components of the geofa array
  !
  do nf = 1,nface
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    nfpts = geom_solpts(face_geom,face_order)
    !
    allocate ( geofa(nf)%v(1:nr,1:nfpts) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) nf
    call alloc_error(pname,array_name,1,ierr,error_message)
    !
  end do
  !
  ! Allocate some arrays needed for curved boundaries
  !
  if (using_curved_boundaries) then
    !
    allocate ( metrics2(1:nr,1:nr,1:nr,1:n_totpts) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"metrics2",1,ierr,error_message)
    !
    allocate ( dJ(1:nr,1:n_totpts) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"dJ",1,ierr,error_message)
    !
    if (use_3rd_derivative_metrics) then
      !
      allocate ( metrics3(1:nr,1:nr,1:nr,1:nr,1:n_totpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"metrics3",1,ierr,error_message)
      !
    end if
    !
  end if
  !
  local_volume = zero
  total_volume = zero
  !
  ! Loop over the cells to get the metrics at each solution point
  !
  metrics_loop: do nc = 1,ncell
    !
    n1 = nodes_of_cell_ptr(nc)+1
    n2 = nodes_of_cell_ptr(nc+1)
    !
    p1 = pts_in_cell_ptr(nc)+1
    p2 = pts_in_cell_ptr(nc+1)
    !
    this_geom = cell_geom(nc)
    this_order = cell_order(nc)
    !
    ! Associate weights with the quadrature weights of the
    ! current cells geometry and order
    !
    weights => std_elem(this_geom,this_order)%wts
    !
    ! Load the nodes defining the cell into a temporary array
    !
    do n = n1,n2
      nodes(:,n-n1+1) = xyz_nodes(:,nodes_of_cell(n))
    end do
    !
    ! Initialize the volume for this cell to zero
    !
    tvol = zero
    !
    ! Compute at each solution point (depending on the cell type)
    ! x metrics - dx/dr,dx/ds,d2x/dr2,d2x/ds2,d2x/drds,
    ! y metrics - dy/dr,dy/ds,d2y/dr2,d2y/ds2,d2y/drds,
    ! Jacobians - J,dJ/dr,dJ/ds
    !
    if (cell_geom(nc) == Geom_Edge) then
      !
      ! ###################################
      ! #####  Metrics for edge cell  #####
      ! ###################################
      !
      do n = p1,p2
        !
        np = n-p1+1
        !
        ! First derivative metrics
        !
        metrics(1:nr,1,n,1) = half * (nodes(1:nr,2) - nodes(1:nr,1))
        !                  _____________________________________
        ! jacobian : J = \/ (dx/dr)**2 + (dy/dr)**2 + (dz/dr)**2
        !
        jacobian(n) = rnorm2( metrics(1:nr,1,n,1) )
        !
        ! Add the jacobian and quadrature weight to this cells volume
        !
        tvol = tvol + jacobian(n)*weights(np)
        !
      end do
      !
   else if (cell_geom(nc) == Geom_Tria) then
      !
      ! ########################################
      ! #####  Metrics for triangular cell #####
      ! ########################################
      !
      do n = p1,p2
        !
        np = n-p1+1
        !
        ! First derivative metrics
        !
        metrics(1:nr,1,n,1) = tria_coef * (nodes(1:nr,2) - nodes(1:nr,1))
        !
        metrics(1:nr,2,n,1) = tria_coef * (nodes(1:nr,3) - nodes(1:nr,1))
        !
        ! Second derivative metrics
        !
        if (allocated(metrics2)) then
          metrics2(:,:,:,n) = zero
        end if
        !
        ! jacobian : J = ( dx/dr * dy/ds ) - ( dx/ds * dy/dr )
        !
        jacobian(n) = compute_jacobian( metrics(1:nr,1:nr,n,1) )
        !
        ! Compute the inverse metrics
        !
        metrics(1:nr,1:nr,n,2) = reverse_metrics( metrics(1:nr,1:nr,n,1) )
        !
        ! First derivative of jacobian
        !
        if (allocated(dJ)) then
          !
          ! dJ/dr = [ ( d2x/dr2 * dy/ds ) + ( dx/dr * d2y/drds ) ] - 
          !         [ ( d2x/drds * dy/dr ) + ( dx/ds * d2y/dr2 ) ]
          !
          dJ(1,n) = metrics2(1,1,1,n)*metrics(2,2,n,1) + &
                    metrics2(2,1,2,n)*metrics(1,1,n,1) - &
                    metrics2(1,1,2,n)*metrics(2,1,n,1) - &
                    metrics2(2,1,1,n)*metrics(1,2,n,1)
          !
          ! dJ/ds = [ ( d2x/dsdr * dy/ds ) + ( dx/dr * d2y/ds2 ) ] - 
          !         [ ( d2x/ds2 * dy/dr ) + ( dx/ds * d2y/dsdr ) ]
          !
          dJ(2,n) = metrics2(1,1,2,n)*metrics(2,2,n,1) + &
                    metrics2(2,2,2,n)*metrics(1,1,n,1) - &
                    metrics2(1,2,2,n)*metrics(2,1,n,1) - &
                    metrics2(2,1,2,n)*metrics(1,2,n,1)
          !
        end if
        !
        ! Add the jacobian and quadrature weight to this cells volume
        !
        tvol = tvol + jacobian(n)*weights(np)
        !
      end do
      !
    else if (cell_geom(nc) == Geom_Quad) then
      !
      ! ############################################
      ! #####  Metrics for quadrilateral cell  #####
      ! ############################################
      !
      n_edge_pts = geom_solpts(Geom_Edge,cell_order(nc))
      !
      edge_pts => std_elem(Geom_Edge,cell_order(nc))%pts(1,:)
      !
      do n = p1,p2
        !
        np = n-p1+1
        !
        j = (np-1)/n_edge_pts + 1
        i = np - (j-1)*n_edge_pts
        !
        rm1 = edge_pts(i) - one
        rp1 = edge_pts(i) + one
        !
        sm1 = edge_pts(j) - one
        sp1 = edge_pts(j) + one
        !
        ! First derivative metrics
        !
        metrics(1:nr,1,n,1) = one4 * ( sm1*(nodes(1:nr,1)-nodes(1:nr,2)) + &
                                       sp1*(nodes(1:nr,3)-nodes(1:nr,4)) )
        !
        metrics(1:nr,2,n,1) = one4 * ( rm1*(nodes(1:nr,1)-nodes(1:nr,4)) + &
                                       rp1*(nodes(1:nr,3)-nodes(1:nr,2)) )
        !
        ! Second derivative metrics
        !
        if (allocated(metrics2)) then
          !
          metrics2(1:nr,1,2,n) = one4 * ( nodes(1:nr,1) - nodes(1:nr,2) + &
                                          nodes(1:nr,3) - nodes(1:nr,4) )
          !
          metrics2(1:nr,2,1,n) = metrics2(1:nr,1,2,n)
          !
        end if
        !
        ! jacobian : J = ( dx/dr * dy/ds ) - ( dx/ds * dy/dr )
        !
        jacobian(n) = compute_jacobian( metrics(1:nr,1:nr,n,1) )
        !
        ! Compute the inverse metrics
        !
        metrics(1:nr,1:nr,n,2) = reverse_metrics( metrics(1:nr,1:nr,n,1) )
        !
        ! First derivative of jacobian
        !
        if (allocated(dJ)) then
          !
          ! dJ/dr = [ ( d2x/dr2 * dy/ds ) + ( dx/dr * d2y/drds ) ] - 
          !         [ ( d2x/drds * dy/dr ) + ( dx/ds * d2y/dr2 ) ]
          !
          dJ(1,n) = metrics2(1,1,1,n)*metrics(2,2,n,1) + &
                    metrics2(2,1,2,n)*metrics(1,1,n,1) - &
                    metrics2(1,1,2,n)*metrics(2,1,n,1) - &
                    metrics2(2,1,1,n)*metrics(1,2,n,1)
          !
          ! dJ/ds = [ ( d2x/dsdr * dy/ds ) + ( dx/dr * d2y/ds2 ) ] - 
          !         [ ( d2x/ds2 * dy/dr ) + ( dx/ds * d2y/dsdr ) ]
          !
          dJ(2,n) = metrics2(1,1,2,n)*metrics(2,2,n,1) + &
                    metrics2(2,2,2,n)*metrics(1,1,n,1) - &
                    metrics2(1,2,2,n)*metrics(2,1,n,1) - &
                    metrics2(2,1,2,n)*metrics(1,2,n,1)
          !
        end if
        !
        ! Add the jacobian and quadrature weight to this cells volume
        !
        tvol = tvol + jacobian(n)*weights(np)
        !
      end do
      !
      ! ##########################################################################
      ! ##########################################################################
      ! ##########################################################################
      ! ###############  ||||||||||                   ||||||||||  ################
      ! ###############  VVVVVVVVVV  NEED TO WORK ON  VVVVVVVVVV  ################
      ! ##########################################################################
      ! ##########################################################################
      ! ##########################################################################
      !
    else if (cell_geom(nc) == Geom_Tetr) then
      !
      ! ##########################################
      ! #####  Metrics for tetrahedral cell  #####
      ! ##########################################
      !
      write (error_message,101)
      call stop_recon(stop_mpi,pname,error_message)
      !
    else if (cell_geom(nc) == Geom_Pyra) then
      !
      ! ######################################
      ! #####  Metrics for pyramid cell  #####
      ! ######################################
      !
      write (error_message,102)
      call stop_recon(stop_mpi,pname,error_message)
      !
      ! ##########################################################################
      ! ##########################################################################
      ! ##########################################################################
      ! ###############  AAAAAAAAAA  NEED TO WORK ON  AAAAAAAAAA  ################
      ! ###############  ||||||||||                   ||||||||||  ################
      ! ##########################################################################
      ! ##########################################################################
      ! ##########################################################################
      !
    else if (cell_geom(nc) == Geom_Pris) then
      !
      ! ####################################
      ! #####  Metrics for prism cell  #####
      ! ####################################
      !
      n_edge_pts = geom_solpts(Geom_Edge,cell_order(nc))
      !
      edge_pts => std_elem(Geom_Edge,cell_order(nc))%pts(1,:)
      rs_pts => std_elem(Geom_Tria,cell_order(nc))%pts
      !
      do n = p1,p2
        !
        np = n-p1+1
        !
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        ! ###############  ||||||||||                   ||||||||||  ################
        ! ###############  VVVVVVVVVV   NEEDS WORK!!!   VVVVVVVVVV  ################
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        !
        k = (np-1)/(n_edge_pts**2) + 1
        j = (np-1)/n_edge_pts - (k-1)*n_edge_pts + 1
        i = np - (j-1)*n_edge_pts - (k-1)*n_edge_pts*n_edge_pts
        !
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        ! ###############  AAAAAAAAAA   NEEDS WORK!!!   AAAAAAAAAA  ################
        ! ###############  ||||||||||                   ||||||||||  ################
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        !
        r = tria_coef * (rs_pts(1,i) + deltafun(tria_is_biunit))
        s = tria_coef * (rs_pts(2,i) + deltafun(tria_is_biunit))
        t = tria_coef * (deltafun(tria_is_biunit) - r - s)
        !
        zm1 = tria_coef * (edge_pts(j) - one)
        zp1 = tria_coef * (edge_pts(j) + one)
        !
        metrics(1:nr,1,n,1) = half * ( zm1*(nodes(1:nr,1) - nodes(1:nr,2)) + &
                                       zp1*(nodes(1:nr,5) - nodes(1:nr,4)) )
        !
        metrics(1:nr,2,n,1) = half * ( zm1*(nodes(1:nr,1) - nodes(1:nr,3)) + &
                                       zp1*(nodes(1:nr,6) - nodes(1:nr,4)) )
        !
        metrics(1:nr,3,n,1) = half * ( r*(nodes(1:nr,5) - nodes(1:nr,2)) + &
                                       s*(nodes(1:nr,6) - nodes(1:nr,3)) + &
                                       t*(nodes(1:nr,4) - nodes(1:nr,1)) )
        !
        ! Second derivative metrics
        !
        if (allocated(metrics2)) then
          !
          metrics2(1:nr,2,3,n) = tc2 * ( (nodes(1:nr,1) - nodes(1:nr,3)) + &
                                         (nodes(1:nr,6) - nodes(1:nr,4)) )
          !
          metrics2(1:nr,3,2,n) = tc2 * ( (nodes(1:nr,1) - nodes(1:nr,4)) + &
                                         (nodes(1:nr,6) - nodes(1:nr,3)) )
          !
          metrics2(1:nr,3,1,n) = tc2 * ( (nodes(1:nr,1) - nodes(1:nr,4)) + &
                                         (nodes(1:nr,5) - nodes(1:nr,2)) )
          !
          metrics2(1:nr,1,3,n) = tc2 * ( (nodes(1:nr,1) - nodes(1:nr,2)) + &
                                         (nodes(1:nr,5) - nodes(1:nr,4)) )
          !
        end if
        !
        ! Third derivative metrics
        !
        if (allocated(metrics3)) then
          metrics3(:,:,:,:,n) = zero
        end if
        !
        ! jacobian : J = ( dx/dr * dy/ds ) - ( dx/ds * dy/dr )
        !
        jacobian(n) = compute_jacobian( metrics(1:nr,1:nr,n,1) )
        !
        ! Compute the inverse metrics
        !
        metrics(1:nr,1:nr,n,2) = reverse_metrics( metrics(1:nr,1:nr,n,1) )
        !
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        ! ###############  ||||||||||                   ||||||||||  ################
        ! ###############  VVVVVVVVVV  NEED TO WORK ON  VVVVVVVVVV  ################
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        !
        ! First derivative of jacobian
        !
        if (allocated(dJ)) then
          !
          ! dJ/dr = [ ( d2x/dr2 * dy/ds ) + ( dx/dr * d2y/drds ) ] - 
          !         [ ( d2x/drds * dy/dr ) + ( dx/ds * d2y/dr2 ) ]
          !
          dJ(1,n) = metrics2(1,1,1,n)*metrics(2,2,n,1) + &
                    metrics2(2,1,2,n)*metrics(1,1,n,1) - &
                    metrics2(1,1,2,n)*metrics(2,1,n,1) - &
                    metrics2(2,1,1,n)*metrics(1,2,n,1)
          !
          ! dJ/ds = [ ( d2x/dsdr * dy/ds ) + ( dx/dr * d2y/ds2 ) ] - 
          !         [ ( d2x/ds2 * dy/dr ) + ( dx/ds * d2y/dsdr ) ]
          !
          dJ(2,n) = metrics2(1,1,2,n)*metrics(2,2,n,1) + &
                    metrics2(2,2,2,n)*metrics(1,1,n,1) - &
                    metrics2(1,2,2,n)*metrics(2,1,n,1) - &
                    metrics2(2,1,2,n)*metrics(1,2,n,1)
          !
        end if
        !
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        ! ###############  AAAAAAAAAA  NEED TO WORK ON  AAAAAAAAAA  ################
        ! ###############  ||||||||||                   ||||||||||  ################
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        !
        ! Add the jacobian and quadrature weight to this cells volume
        !
        tvol = tvol + jacobian(n)*weights(np)
        !
      end do
      !
      !
      write (error_message,103)
      call stop_recon(stop_mpi,pname,error_message)
      !
    else if (cell_geom(nc) == Geom_Hexa) then
      !
      ! #########################################
      ! #####  Metrics for hexahedral cell  #####
      ! #########################################
      !
      n_edge_pts = geom_solpts(Geom_Edge,cell_order(nc))
      !
      edge_pts => std_elem(Geom_Edge,cell_order(nc))%pts(1,:)
      !
      do n = p1,p2
        !
        np = n-p1+1
        !
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        ! ###############  ||||||||||                   ||||||||||  ################
        ! ###############  VVVVVVVVVV DOUBLE CHECK THIS VVVVVVVVVV  ################
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        !
        k = (np-1)/(n_edge_pts**2) + 1
        j = (np-1)/n_edge_pts - (k-1)*n_edge_pts + 1
        i = np - (j-1)*n_edge_pts - (k-1)*n_edge_pts*n_edge_pts
        !
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        ! ###############  AAAAAAAAAA DOUBLE CHECK THIS AAAAAAAAAA  ################
        ! ###############  ||||||||||                   ||||||||||  ################
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        !
        rm1 = edge_pts(i) - one
        rp1 = edge_pts(i) + one
        !
        sm1 = edge_pts(j) - one
        sp1 = edge_pts(j) + one
        !
        tm1 = edge_pts(k) - one
        tp1 = edge_pts(k) + one
        !
        metrics(1:nr,1,n,1) = one8 * ( sm1*tm1*(nodes(1:nr,2) - nodes(1:nr,1)) + &
                                       sm1*tp1*(nodes(1:nr,5) - nodes(1:nr,6)) + &
                                       sp1*tm1*(nodes(1:nr,4) - nodes(1:nr,3)) + &
                                       sp1*tp1*(nodes(1:nr,7) - nodes(1:nr,8)) )
        !
        metrics(1:nr,2,n,1) = one8 * ( rm1*tm1*(nodes(1:nr,4) - nodes(1:nr,1)) + &
                                       rm1*tp1*(nodes(1:nr,5) - nodes(1:nr,8)) + &
                                       rp1*tm1*(nodes(1:nr,2) - nodes(1:nr,3)) + &
                                       rp1*tp1*(nodes(1:nr,7) - nodes(1:nr,6)) )
        !
        metrics(1:nr,3,n,1) = one8 * ( rm1*sm1*(nodes(1:nr,5) - nodes(1:nr,1)) + &
                                       rm1*sp1*(nodes(1:nr,4) - nodes(1:nr,8)) + &
                                       rp1*sm1*(nodes(1:nr,2) - nodes(1:nr,6)) + &
                                       rp1*sp1*(nodes(1:nr,7) - nodes(1:nr,3)) )
        !
        ! Second derivative metrics
        !
        if (allocated(metrics2)) then
          !
          metrics2(1:nr,2,1,n) = one8 * ( tm1*(nodes(1:nr,2) - nodes(1:nr,1)) + &
                                          tm1*(nodes(1:nr,4) - nodes(1:nr,3)) + &
                                          tp1*(nodes(1:nr,5) - nodes(1:nr,6)) + &
                                          tp1*(nodes(1:nr,7) - nodes(1:nr,8)) )
          metrics2(1:nr,1,2,n) = metrics2(1:nr,2,1,n)
          !
          metrics2(1:nr,3,2,n) = one8 * ( rm1*(nodes(1:nr,4) - nodes(1:nr,1)) + &
                                          rm1*(nodes(1:nr,5) - nodes(1:nr,8)) + &
                                          rp1*(nodes(1:nr,2) - nodes(1:nr,3)) + &
                                          rp1*(nodes(1:nr,7) - nodes(1:nr,6)) )
          metrics2(1:nr,2,3,n) = metrics2(1:nr,3,2,n)
          !
          metrics2(1:nr,1,3,n) = one8 * ( sm1*(nodes(1:nr,5) - nodes(1:nr,1)) + &
                                          sm1*(nodes(1:nr,2) - nodes(1:nr,6)) + &
                                          sp1*(nodes(1:nr,4) - nodes(1:nr,8)) + &
                                          sp1*(nodes(1:nr,7) - nodes(1:nr,3)) )
          metrics2(1:nr,3,1,n) = metrics2(1:nr,1,3,n)
          !
        end if
        !
        ! Third derivative metrics
        !
        if (allocated(metrics3)) then
          !
          metrics3(1:nr,3,2,1,n) = one8 * ( nodes(1:nr,2) - nodes(1:nr,1) + &
                                            nodes(1:nr,4) - nodes(1:nr,3) + &
                                            nodes(1:nr,5) - nodes(1:nr,6) + &
                                            nodes(1:nr,7) - nodes(1:nr,8) )
          metrics3(1:nr,2,3,1,n) = metrics3(1:nr,3,2,1,n)
          metrics3(1:nr,3,1,2,n) = metrics3(1:nr,3,2,1,n)
          metrics3(1:nr,1,3,2,n) = metrics3(1:nr,3,2,1,n)
          metrics3(1:nr,2,1,3,n) = metrics3(1:nr,3,2,1,n)
          metrics3(1:nr,1,2,3,n) = metrics3(1:nr,3,2,1,n)
          !
        end if
        !
        ! jacobian : J = ( dx/dr * dy/ds ) - ( dx/ds * dy/dr )
        !
        jacobian(n) = compute_jacobian( metrics(1:nr,1:nr,n,1) )
        !
        ! Compute the inverse metrics
        !
        metrics(1:nr,1:nr,n,2) = reverse_metrics( metrics(1:nr,1:nr,n,1) )
        !
        ! First derivative of jacobian
        !
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        ! ###############  ||||||||||                   ||||||||||  ################
        ! ###############  VVVVVVVVVV  NEED TO WORK ON  VVVVVVVVVV  ################
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        !
        if (allocated(dJ)) then
          !
          ! dJ/dr = [ ( d2x/dr2 * dy/ds ) + ( dx/dr * d2y/drds ) ] - 
          !         [ ( d2x/drds * dy/dr ) + ( dx/ds * d2y/dr2 ) ]
          !
          dJ(1,n) = metrics2(1,1,1,n)*metrics(2,2,n,1) + &
                    metrics2(2,1,2,n)*metrics(1,1,n,1) - &
                    metrics2(1,1,2,n)*metrics(2,1,n,1) - &
                    metrics2(2,1,1,n)*metrics(1,2,n,1)
          !
          ! dJ/ds = [ ( d2x/dsdr * dy/ds ) + ( dx/dr * d2y/ds2 ) ] - 
          !         [ ( d2x/ds2 * dy/dr ) + ( dx/ds * d2y/dsdr ) ]
          !
          dJ(2,n) = metrics2(1,1,2,n)*metrics(2,2,n,1) + &
                    metrics2(2,2,2,n)*metrics(1,1,n,1) - &
                    metrics2(1,2,2,n)*metrics(2,1,n,1) - &
                    metrics2(2,1,2,n)*metrics(1,2,n,1)
          !
        end if
        !
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        ! ###############  AAAAAAAAAA  NEED TO WORK ON  AAAAAAAAAA  ################
        ! ###############  ||||||||||                   ||||||||||  ################
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        !
        ! Add the jacobian and quadrature weight to this cells volume
        !
        tvol = tvol + jacobian(n)*weights(np)
        !
      end do
      !
    else
      !
      ! Unknown or unsupported geometry
      !
      write (error_message,100)
      call stop_recon(stop_mpi,pname,error_message)
      !
    end if
    !
    ! Store the computed volume of the current cell
    !
    volj(nc) = tvol
    !
    ! Add the volume of the current cell to the total volume
    ! of the current local partition
    !
    if (nc <= ncell) local_volume = local_volume + volj(nc)
    !
    ! Disassociate local pointers before moving to the next cell
    !
    if (associated(weights)) weights => null()
    if (associated(edge_pts)) edge_pts => null()
    if (associated(rs_pts)) rs_pts => null()
    !
  end do metrics_loop
 !write (*,*) mypnum,local_volume,sum(volj(1:ncell))
  !
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  ! ###############  ||||||||||                   ||||||||||  ################
  ! ###############  VVVVVVVVVV  NEED TO WORK ON  VVVVVVVVVV  ################
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  !
  ! Get the face normal at the solution points on each face
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
  do nf = 1,nface
    !
    nfnods = geom_nodes(face(nf)%geom) ! # of nodes defining this face
    !
    ! Coordinates of the nodes on this face
    !
    nodes(1:nr,1:nfnods) = xyz_nodes(1:nr,nodes_on_face(1:nfnods,nf))
    !
    ! !!!!!!!!!!!!!! WARNING: POSSIBLY BAD ASSUMPTIONS MADE HERE !!!!!!!!!!!!!!
    ! Assumption 1: the normal is the same for all flux points on the face
    ! !!!!!!!!!!!!!! WARNING: POSSIBLY BAD ASSUMPTIONS MADE HERE !!!!!!!!!!!!!!
    !
    ! Compute the face normal at the first flux point
    !
    geofa(nf)%v(1:nr,1) = face_normal( nodes(1:nr,1:nfnods) )
    !
    ! Copy the face normal from the first flux point to all the other
    ! face flux points (assuming the normal is the same for all of them)
    !
    do n = 2,size(geofa(nf)%v,dim=2)
      geofa(nf)%v(1:nr,n) = geofa(nf)%v(1:nr,1)
    end do
    !
  end do
  !
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  ! ###############  AAAAAAAAAA  NEED TO WORK ON  AAAAAAAAAA  ################
  ! ###############  ||||||||||                   ||||||||||  ################
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  !
  ! Reduce the volumes of each partition to the root
  ! processor to determine the volume of the global grid.
  !
  tot_vol = real(local_volume,kind=r8)
  loc_vol = real(local_volume,kind=r8)
  !
  if (ncpu > 1) then
    call mpi_allreduce(loc_vol,tot_vol,1,mpi_flttyp, &
                       mpi_sum,mpi_comm_world,mpierr)
  end if
  !
  total_volume = real(tot_vol,kind=wp)
  local_volume = real(loc_vol,kind=wp)
  !
  call mpi_barrier(mpi_comm_world,mpierr)
  !
  ! Format Statements
  !
  1   format ("geofa(",i0,")%v")
  100 format (" Encountered an unknown or unsupported cell geometry while", &
              " trying to compute the grid metrics!")
  101 format (" Encountered a tetrahedral grid cell while trying to compute", &
              " the grid metrics! This grid geometry is not currently supported!")
  102 format (" Encountered a pyramid grid cell while trying to compute", &
              " the grid metrics! This grid geometry is not currently supported!")
  103 format (" Encountered a prism grid cell while trying to compute", &
              " the grid metrics! This grid geometry is not currently supported!")
  104 format (" Encountered an invalid metrics derived type!")
  !
end subroutine get_metrics
!
!###############################################################################
!
subroutine metrics_alloc_status(iunit)
  !
  !.. Use Statements ..
  use module_memory
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iunit
  !
  !.. Local Scalars ..
  type(memory) :: total
  type(memory) :: array
  !
continue
  !
  total%siz = 0
  total%mem = 0
  !
  write (iunit,1)
  !
  if (allocated(jacobian)) then
    array%siz = size(jacobian)
    array%mem = array%siz*storage_size(jacobian)
    total%mem = total%mem + array%mem
    write (iunit,2) "'jacobian'",array%get_units()
  end if
  if (allocated(dJ)) then
    array%siz = size(dJ)
    array%mem = array%siz*storage_size(dJ)
    total%mem = total%mem + array%mem
    write (iunit,2) "'dJ'",array%get_units()
  end if
  if (allocated(volj)) then
    array%siz = size(volj)
    array%mem = array%siz*storage_size(volj)
    total%mem = total%mem + array%mem
    write (iunit,2) "'volj'",array%get_units()
  end if
  !
  total = total%get_units()
  write (iunit,3) total%mem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module metrics_mod")
  2 format (" Array: ",a11,dt)
  3 format (" Total memory = ",i0,1x,a)
  !
end subroutine metrics_alloc_status
!
!###############################################################################
!
end module metrics_mod
