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
  public :: metrics_memory_usage
  !
  type, public :: face_nrm
    real(wp), allocatable :: v(:,:)
  end type face_nrm
  !
  ! geofa : normal of the face at each of the solution points on the face
  !
  type(face_nrm), public, save, allocatable :: geofa(:)
  type(face_nrm), public, save, allocatable :: geofa2(:)
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
  ! metrics3: third derivative metrics terms at the solution points
  !     For solution point 'n', metrics3(l,i,j,k,n) gives the partial derivative
  !                                      3
  !                                     d (x_l)
  !                             -----------------------
  !                             d(xi_i) d(xi_j) d(xi_k)
  !     where ( x_1, x_2, x_3) = (x,y,z) -> the physical coordinates
  !     where (xi_1,xi_2,xi_3) = (xi,eta,zeta) or (r,s,t) -> the computational
  !                           coordinates within the standard elemental region
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
  !           allocation, it defaults to the value 2 given in the type
  !           definition.
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
 !real(wp), public, save, allocatable :: spmat(:,:)       !(3,ksp2d)
 !real(wp), public, save, allocatable :: spcemat(:,:)     !(ksp1d+1,ksp2d)
 !real(wp), public, save, allocatable :: spceximat(:,:)   !(ksp1d+1,ksp2d)
 !real(wp), public, save, allocatable :: spceetmat(:,:)   !(ksp1d+1,ksp2d)
 !real(wp), public, save, allocatable :: spcexi2mat(:,:)  !(ksp1d+1,ksp2d)
 !real(wp), public, save, allocatable :: spceet2mat(:,:)  !(ksp1d+1,ksp2d)
 !real(wp), public, save, allocatable :: spcexietmat(:,:) !(ksp1d+1,ksp2d)
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
  ! XYSPxi  : d_(x,y)/d_xi  of the cell solution points, i.e. d_x/d_xi, d_y/d_xi
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
 !real(wp), public, save, allocatable :: XYSP(:,:,:)    !(2,ksp2d,-nfbnd:ncell)
 !real(wp), public, save, allocatable :: XYSPxi(:,:,:)  !(2,ksp2d,-nfbnd:ncell)
 !real(wp), public, save, allocatable :: XYSPet(:,:,:)  !(2,ksp2d,-nfbnd:ncell)
 !real(wp), public, save, allocatable :: XYeSP(:,:,:)   !(2,ksp1d,1:nface)
 !real(wp), public, save, allocatable :: XYeSPxi(:,:,:) !(2,ksp1d,1:nface)
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
 !real(wp), public, save, allocatable :: rjaSP(:,:)   !(ksp2d,-nfbnd:ncell)
 !real(wp), public, save, allocatable :: rjaSPxi(:,:) !(ksp2d,-nfbnd:ncell)
 !real(wp), public, save, allocatable :: rjaSPet(:,:) !(ksp2d,-nfbnd:ncell)
  !
  !... for boundary cells, if id_be(j) = 0, then the cell has a curved edge
  !
  ! be_cSPind(1,k,j) k = 1,ksp: the solution point index of cell be_ic(j)
  ! be_cSPind(2,k,j) k = 1,ksp: the solution point index of exterior cell
  ! be_cSPind(3,k,j) k = 1,ksp: the solution point index of cell be_pc(j)
  ! e_cSPind(1,k,j) k = 1,ksp: the solution point index for the left cell
  ! e_cSPind(2,k,j) k = 1,ksp: the solution point index for the right cell
 !integer, public, save, allocatable :: e_cSPind(:,:,:)  !(2,ksp1d,nface)
 !integer, public, save, allocatable :: be_cSPind(:,:,:) !(3,ksp1d,nfbnd)
  !
contains
!
!###############################################################################
!
subroutine get_metrics()
  !
 !use ifport
  !.. Use Statements ..
  use geovar, only : nr,nface,n_totpts,n_totcel,ncell
  use geovar, only : nodes_of_cell_ptr,nodes_of_cell
  use geovar, only : xyz_nodes,xyz
  use geovar, only : face,cell
  !
  use derivatives_mod, only : deriv
  !
  use quadrature_mod, only : std_elem
  !
  use ovar, only : local_volume,total_volume
  use ovar, only : using_curved_boundaries
  use ovar, only : use_3rd_derivative_metrics
  use ovar, only : use_advanced_face_metrics
  use ovar, only : VortexGrid_random_perturb
  !
  use order_mod, only : geom_solpts,maxSP
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,ierr,j1,j2,ibeg,iend
  integer  :: n,nc,n1,n2,nf,np,p1,p2,p3
  integer  :: n_edge_pts,nfnods,nfpts
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  integer  :: ref_dir,phys_dir
  real(wp) :: tvol,tria_coef,tc2,r,s,t
  real(r8) :: loc_vol,tot_vol
  real(wp) :: val,rm1,rp1,sm1,sp1,tm1,tp1,zm1,zp1
  character(len=200) :: array_name
  !
 !integer :: exitstat,cmndstat
 !character(len=1000) :: cmndmsg,command
  !
  !.. Local Arrays ..
  integer, dimension(1:3) :: ksize
  integer, dimension(1:4,1:136,1:3) :: knodes
  real(wp), dimension(1:nr,1:16) :: diff
  real(wp), dimension(1:nr,1:8) :: nodes
  real(wp), dimension(1:maxSP,1:nr) :: txyz
  real(wp), dimension(1:nr,1:nr) :: outmet
  !
  integer, dimension(1:4) :: iom
  !
  !.. Local Pointers ..
  real(wp), pointer, contiguous :: weights(:)
  real(wp), pointer, contiguous :: edge_pts(:)
  real(wp), pointer, contiguous :: rs_pts(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_metrics"
  character(len=*), parameter :: fname(4) = ["metrics.fwd.dat", &
                                             "metrics.rev.dat", &
                                             "uniform.fwd.dat", &
                                             "uniform.rev.dat"]
  !
 !logical(lk), parameter :: output_metrics = true
  logical(lk), parameter :: output_metrics = fals
  !
  logical(lk), parameter :: use_derivative_matrices = true
 !logical(lk), parameter :: use_derivative_matrices = fals
  !
  character(len=*), parameter :: cxyz(3,2) = reshape( ["X","Y","Z", &
                                                       "R","S","T"], &
                                                      shape(cxyz) )
  !
continue
  !
  if (VortexGrid_random_perturb) then
    ibeg = 1
    iend = 2
  else
    ibeg = 3
    iend = 4
  end if
  if (output_metrics) then
    do i = ibeg,iend
      open(newunit=iom(i),file=fname(i),status="replace",action="write", &
                          form="formatted",position="rewind", &
                          iostat=ierr,iomsg=error_message)
      call io_error(pname,fname(i),1,ierr,error_message)
    end do
  end if
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
  call alloc_error(pname,"metrics",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( jacobian(1:n_totpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"jacobian",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( volj(1:n_totcel) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"volj",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate the outer part of the geofa array
  !
  allocate ( geofa(1:nface) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"geofa",1,__LINE__,__FILE__,ierr,error_message)
 !allocate ( geofa2(1:nface) , stat=ierr , errmsg=error_message )
 !call alloc_error(pname,"geofa2",1,__LINE__,__FILE__,ierr,error_message)
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
    write (array_name,1) "geofa",nf
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !
   !allocate ( geofa2(nf)%v(1:nr,1:nfpts) , source=zero , &
   !           stat=ierr , errmsg=error_message )
   !write (array_name,1) "geofa2",nf
   !call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
  end do
  !
  ! Allocate some arrays needed for curved boundaries
  !
  if (using_curved_boundaries) then
    !
    allocate ( metrics2(1:nr,1:nr,1:nr,1:n_totpts) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"metrics2",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( dJ(1:nr,1:n_totpts) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"dJ",1,__LINE__,__FILE__,ierr,error_message)
    !
    if (use_3rd_derivative_metrics) then
      !
      allocate ( metrics3(1:nr,1:nr,1:nr,1:nr,1:n_totpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"metrics3",1,__LINE__,__FILE__,ierr,error_message)
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
    p1 = cell(nc)%beg_sp
    p2 = cell(nc)%end_sp
    p3 = p2-p1+1
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    if (output_metrics) then
      do i = ibeg,iend
        write (iom(i),21) nc
      end do
    end if
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
    if (this_geom == Geom_Edge) then
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
   else if (this_geom == Geom_Tria) then
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
    else if (this_geom == Geom_Quad) then
      !
      ! ############################################
      ! #####  Metrics for quadrilateral cell  #####
      ! ############################################
      !
      n_edge_pts = geom_solpts(Geom_Edge,this_order)
      !
      edge_pts => std_elem(Geom_Edge,this_order)%pts(1,:)
      !
      if (use_derivative_matrices) then
        do k = p1,p2
          do l = 1,nr
            txyz(k-p1+1,l) = xyz(l,k)
          end do
        end do
      end if
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
        if (use_derivative_matrices) then
          !
          do ref_dir = 1,nr
            do phys_dir = 1,nr
              metrics(phys_dir,ref_dir,n,1) = deriv(this_geom,this_order)% &
                                                               d(ref_dir)% &
                                              dot( np , txyz(1:p3,phys_dir) )
            end do
          end do
          !
        else
          !
          do l = 1,nr
            metrics(l,1,n,1) = one4 * ( sm1*(nodes(l,1)-nodes(l,2)) + &
                                        sp1*(nodes(l,3)-nodes(l,4)) )
          end do
          !
          do l = 1,nr
            metrics(l,2,n,1) = one4 * ( rm1*(nodes(l,1)-nodes(l,4)) + &
                                        rp1*(nodes(l,3)-nodes(l,2)) )
          end do
          !
        end if
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
      ! ########################################################################
      ! ########################################################################
      ! ########################################################################
      ! ##############  ||||||||||                   ||||||||||  ###############
      ! ##############  VVVVVVVVVV  NEED TO WORK ON  VVVVVVVVVV  ###############
      ! ########################################################################
      ! ########################################################################
      ! ########################################################################
      !
    else if (this_geom == Geom_Tetr) then
      !
      ! ##########################################
      ! #####  Metrics for tetrahedral cell  #####
      ! ##########################################
      !
      write (error_message,2)
      call stop_recon(stop_mpi,pname,__LINE__,__FILE__,error_message)
      !
    else if (this_geom == Geom_Pyra) then
      !
      ! ######################################
      ! #####  Metrics for pyramid cell  #####
      ! ######################################
      !
      write (error_message,3)
      call stop_recon(stop_mpi,pname,__LINE__,__FILE__,error_message)
      !
      ! ########################################################################
      ! ########################################################################
      ! ########################################################################
      ! ##############  AAAAAAAAAA  NEED TO WORK ON  AAAAAAAAAA  ###############
      ! ##############  ||||||||||                   ||||||||||  ###############
      ! ########################################################################
      ! ########################################################################
      ! ########################################################################
      !
    else if (this_geom == Geom_Pris) then
      !
      ! ####################################
      ! #####  Metrics for prism cell  #####
      ! ####################################
      !
      n_edge_pts = geom_solpts(Geom_Edge,this_order)
      !
      edge_pts => std_elem(Geom_Edge,this_order)%pts(1,:)
      rs_pts => std_elem(Geom_Tria,this_order)%pts
      !
      do n = p1,p2
        !
        np = n-p1+1
        !
        ! ######################################################################
        ! ######################################################################
        ! ######################################################################
        ! #############  ||||||||||                   ||||||||||  ##############
        ! #############  VVVVVVVVVV   NEEDS WORK!!!   VVVVVVVVVV  ##############
        ! ######################################################################
        ! ######################################################################
        ! ######################################################################
        !
        k = (np-1)/(n_edge_pts**2) + 1
        j = (np-1)/n_edge_pts - (k-1)*n_edge_pts + 1
        i = np - (j-1)*n_edge_pts - (k-1)*n_edge_pts*n_edge_pts
        !
        ! ######################################################################
        ! ######################################################################
        ! ######################################################################
        ! #############  AAAAAAAAAA   NEEDS WORK!!!   AAAAAAAAAA  ##############
        ! #############  ||||||||||                   ||||||||||  ##############
        ! ######################################################################
        ! ######################################################################
        ! ######################################################################
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
        ! ######################################################################
        ! ######################################################################
        ! ######################################################################
        ! #############  ||||||||||                   ||||||||||  ##############
        ! #############  VVVVVVVVVV  NEED TO WORK ON  VVVVVVVVVV  ##############
        ! ######################################################################
        ! ######################################################################
        ! ######################################################################
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
        ! ######################################################################
        ! ######################################################################
        ! ######################################################################
        ! #############  AAAAAAAAAA  NEED TO WORK ON  AAAAAAAAAA  ##############
        ! #############  ||||||||||                   ||||||||||  ##############
        ! ######################################################################
        ! ######################################################################
        ! ######################################################################
        !
        ! Add the jacobian and quadrature weight to this cells volume
        !
        tvol = tvol + jacobian(n)*weights(np)
        !
      end do
      !
      !
      write (error_message,4)
      call stop_recon(stop_mpi,pname,__LINE__,__FILE__,error_message)
      !
    else if (this_geom == Geom_Hexa) then
      !
      ! #########################################
      ! #####  Metrics for hexahedral cell  #####
      ! #########################################
      !
      n_edge_pts = geom_solpts(Geom_Edge,this_order)
      !
      edge_pts => std_elem(Geom_Edge,this_order)%pts(1,:)
      !
      if (use_derivative_matrices) then
        do k = p1,p2
          do l = 1,nr
            txyz(k-p1+1,l) = xyz(l,k)
          end do
        end do
      end if
      !
      do n = p1,p2
        !
        np = n-p1+1
        !
        ! ######################################################################
        ! ######################################################################
        ! ######################################################################
        ! #############  ||||||||||                   ||||||||||  ##############
        ! #############  VVVVVVVVVV DOUBLE CHECK THIS VVVVVVVVVV  ##############
        ! ######################################################################
        ! ######################################################################
        ! ######################################################################
        !
        k = (np-1)/(n_edge_pts**2) + 1
        j = (np-1)/n_edge_pts - (k-1)*n_edge_pts + 1
        i = np - (j-1)*n_edge_pts - (k-1)*n_edge_pts*n_edge_pts
        !
        ! ######################################################################
        ! ######################################################################
        ! ######################################################################
        ! #############  AAAAAAAAAA DOUBLE CHECK THIS AAAAAAAAAA  ##############
        ! #############  ||||||||||                   ||||||||||  ##############
        ! ######################################################################
        ! ######################################################################
        ! ######################################################################
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
        if (use_derivative_matrices) then
          !
          do ref_dir = 1,nr
            do phys_dir = 1,nr
              metrics(phys_dir,ref_dir,n,1) = deriv(this_geom,this_order)% &
                                                               d(ref_dir)% &
                                              dot( np , txyz(1:p3,phys_dir) )
            end do
          end do
          !
        else
          !
          do l = 1,nr
            metrics(l,1,n,1) = one8*( sm1*tm1*(nodes(l,2) - nodes(l,1)) + &
                                      sm1*tp1*(nodes(l,5) - nodes(l,6)) + &
                                      sp1*tm1*(nodes(l,4) - nodes(l,3)) + &
                                      sp1*tp1*(nodes(l,7) - nodes(l,8)) )
          end do
          !
          do l = 1,nr
            metrics(l,2,n,1) = one8*( rm1*tm1*(nodes(l,4) - nodes(l,1)) + &
                                      rm1*tp1*(nodes(l,5) - nodes(l,8)) + &
                                      rp1*tm1*(nodes(l,2) - nodes(l,3)) + &
                                      rp1*tp1*(nodes(l,7) - nodes(l,6)) )
          end do
          !
          do l = 1,nr
            metrics(l,3,n,1) = one8*( rm1*sm1*(nodes(l,5) - nodes(l,1)) + &
                                      rm1*sp1*(nodes(l,4) - nodes(l,8)) + &
                                      rp1*sm1*(nodes(l,2) - nodes(l,6)) + &
                                      rp1*sp1*(nodes(l,7) - nodes(l,3)) )
          end do
          !
        end if
       !!
       !metrics(1:nr,1,n,1) = one8*( sm1*tm1*(nodes(1:nr,2) - nodes(1:nr,1)) + &
       !                             sm1*tp1*(nodes(1:nr,5) - nodes(1:nr,6)) + &
       !                             sp1*tm1*(nodes(1:nr,4) - nodes(1:nr,3)) + &
       !                             sp1*tp1*(nodes(1:nr,7) - nodes(1:nr,8)) )
       !!
       !metrics(1:nr,2,n,1) = one8*( rm1*tm1*(nodes(1:nr,4) - nodes(1:nr,1)) + &
       !                             rm1*tp1*(nodes(1:nr,5) - nodes(1:nr,8)) + &
       !                             rp1*tm1*(nodes(1:nr,2) - nodes(1:nr,3)) + &
       !                             rp1*tp1*(nodes(1:nr,7) - nodes(1:nr,6)) )
       !!
       !metrics(1:nr,3,n,1) = one8*( rm1*sm1*(nodes(1:nr,5) - nodes(1:nr,1)) + &
       !                             rm1*sp1*(nodes(1:nr,4) - nodes(1:nr,8)) + &
       !                             rp1*sm1*(nodes(1:nr,2) - nodes(1:nr,6)) + &
       !                             rp1*sp1*(nodes(1:nr,7) - nodes(1:nr,3)) )
        !
        ! Second derivative metrics
        !
        if (allocated(metrics2)) then
          !
          metrics2(1:nr,2,1,n) = one8*( tm1*(nodes(1:nr,2) - nodes(1:nr,1)) + &
                                        tm1*(nodes(1:nr,4) - nodes(1:nr,3)) + &
                                        tp1*(nodes(1:nr,5) - nodes(1:nr,6)) + &
                                        tp1*(nodes(1:nr,7) - nodes(1:nr,8)) )
          metrics2(1:nr,1,2,n) = metrics2(1:nr,2,1,n)
          !
          metrics2(1:nr,3,2,n) = one8*( rm1*(nodes(1:nr,4) - nodes(1:nr,1)) + &
                                        rm1*(nodes(1:nr,5) - nodes(1:nr,8)) + &
                                        rp1*(nodes(1:nr,2) - nodes(1:nr,3)) + &
                                        rp1*(nodes(1:nr,7) - nodes(1:nr,6)) )
          metrics2(1:nr,2,3,n) = metrics2(1:nr,3,2,n)
          !
          metrics2(1:nr,1,3,n) = one8*( sm1*(nodes(1:nr,5) - nodes(1:nr,1)) + &
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
        if (output_metrics) then
          do i = ibeg,iend
            j1 = mod(i+1,2) + 1
            j2 = mod(i  ,2) + 1
            outmet = metrics(1:nr,1:nr,n,j1)
           !if (j1 == 2) then
           !  outmet = outmet / jacobian(n)
           !end if
            where (abs(outmet) < eps14) outmet = zero
            write (iom(i),22) n
            write (iom(i),23) ((cxyz(k,j1),k=1,nr), &
                                cxyz(l,j2), &
                               (outmet(k,l), k=1,nr), l=1,nr)
          end do
        end if
        21 format (/,60("-"),/," Cell ",i0,/,60("-"))
        22 format (5x,"Solution Point ",i4)
        23 format (3(2x,"d(",a,",",a,",",a,")/d",a,"  :  ",3es14.6,:,/))
        !
        ! First derivative of jacobian
        !
        ! ######################################################################
        ! ######################################################################
        ! ######################################################################
        ! #############  ||||||||||                   ||||||||||  ##############
        ! #############  VVVVVVVVVV  NEED TO WORK ON  VVVVVVVVVV  ##############
        ! ######################################################################
        ! ######################################################################
        ! ######################################################################
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
        ! ######################################################################
        ! ######################################################################
        ! ######################################################################
        ! #############  AAAAAAAAAA  NEED TO WORK ON  AAAAAAAAAA  ##############
        ! #############  ||||||||||                   ||||||||||  ##############
        ! ######################################################################
        ! ######################################################################
        ! ######################################################################
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
      write (error_message,6)
      call stop_recon(stop_mpi,pname,__LINE__,__FILE__,error_message)
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
 !if (output_metrics) then
 !  exitstat = 0
 !  cmndstat = 0
 !  do i = ibeg,iend
 !        command = "sed -i " // &
 !                      "-e 's/-0\.000000E+00/      -      /g' " // &
 !                      "-e 's/0\.000000E+00/     -      /g' " // &
 !                  trim(adjustl(fname(i)))
 !    write (iout,194) "Executing: ",trim(adjustl(command))
 !    !
 !    exitstat = system( trim(command) )
 !   !call execute_command_line(trim(command),true,exitstat,cmndstat,cmndmsg)
 !    !
 !    if (exitstat /= 0 .or. cmndstat /= 0) then
 !      write (iout,199) exitstat,cmndstat,trim(adjustl(cmndmsg))
 !    end if
 !  end do
 !  write (error_message,*) "After outputting the metrics because ", &
 !                          "the output_metrics parameter is active!"
 !  call stop_recon(stop_mpi,pname,__LINE__,__FILE__,error_message)
 !  194 format (a,a)
 !  199 format (/," Error attempting to use command line!",/, &
 !                " exit status    = ",i0,/, &
 !                " command status = ",i0,/, &
 !                " error message  = ",a,//)
 !end if
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
  if (use_advanced_face_metrics) then
    !
    call compute_face_normals
    !
  else
    !
    do nf = 1,nface
      !
      nfnods = geom_nodes(face(nf)%geom) ! # of nodes defining this face
      !
      ! Coordinates of the nodes on this face
      !
      nodes(1:nr,1:nfnods) = xyz_nodes(1:nr,face(nf)%nodes(1:nfnods))
      !
      ! !!!!!!!!!!!!!! WARNING: POSSIBLY BAD ASSUMPTIONS MADE HERE !!!!!!!!!!!!!
      ! Assumption 1: the normal is the same for all flux points on the face
      ! !!!!!!!!!!!!!! WARNING: POSSIBLY BAD ASSUMPTIONS MADE HERE !!!!!!!!!!!!!
      !
      ! Compute the face normal at the first flux point
      !
     !geofa2(nf)%v(1:nr,1) = face_normal( nodes(1:nr,1:nfnods) )
      geofa(nf)%v(1:nr,1) = face_normal( nodes(1:nr,1:nfnods) )
      !
      ! Copy the face normal from the first flux point to all the other
      ! face flux points (assuming the normal is the same for all of them)
      !
     !do n = 2,size(geofa2(nf)%v,dim=2)
     !  geofa2(nf)%v(1:nr,n) = geofa2(nf)%v(1:nr,1)
     !end do
      do n = 2,size(geofa(nf)%v,dim=2)
        geofa(nf)%v(1:nr,n) = geofa(nf)%v(1:nr,1)
      end do
      !
     !write (180,10) nf
     !do n = 1,size(geofa(nf)%v,dim=2)
     !  write (180,11) n,(one4*geofa(nf)%v(i,n),i=1,size(geofa(nf)%v,dim=1))
     !end do
    end do
    !
  end if
  !
 !do nf = 1,nface
 !  write (182,10) nf
 !  do n = 1,size(geofa(nf)%v,dim=2)
 !    write (182,11) n,(geofa(nf)%v(i,n),i=1,size(geofa(nf)%v,dim=1))
 !  end do
 !end do
 !!
 !do nf = 1,nface
 !  write (182,10) nf
 !  do n = 1,size(geofa2(nf)%v,dim=2)
 !    write (182,11) n,( abs( (one4*geofa(nf)%v(i,n)-geofa2(nf)%v(i,n)) * &
 !                  reciprocal(one4*geofa(nf)%v(i,n)) ), &
 !                          i=1,size(geofa2(nf)%v,dim=1))
 !  end do
 !end do
 !!
 !open (newunit=nfpts,file="faces.tec.dat",status="replace",action="write", &
 !                    form="formatted",position="rewind",iostat=ierr)
 !write (nfpts,20) size(xyz_nodes,dim=2),nface
 !do l = 1,nr
 !  write (nfpts,21) (xyz_nodes(l,n),n=1,size(xyz_nodes,dim=2))
 !end do
 !do nf = 1,nface
 !  nfnods = geom_nodes(face(nf)%geom)
 !  write (nfpts,22) (face(nf)%nodes(l),l=1,nfnods)
 !end do
 !ksize = 0
 !do nf = 1,nface
 !  l = size(geofa(nf)%v,dim=1)
 !  n = size(geofa(nf)%v,dim=2)
 !  diff(1:l,1:n) = abs( (one4*geofa2(nf)%v - geofa(nf)%v) * &
 !                             reciprocal(one4*geofa2(nf)%v) )
 !  if (any(diff > 1.0e-02_wp)) then
 !    nfnods = geom_nodes(face(nf)%geom)
 !    nodes(1:nr,1:nfnods) = xyz_nodes(1:nr,face(nf)%nodes(1:nfnods))
 !   !write (183,12) nf,((nodes(l,n),l=1,nr),n=1,nfnods)
 !    if (any(abs(nodes(3,:)-0.0_wp) < eps9)) then
 !      ksize(1) = ksize(1) + 1
 !      knodes(1:nfnods,ksize(1),1) = face(nf)%nodes(1:nfnods)
 !    else if (any(abs(nodes(3,:)-0.075_wp) < eps9)) then
 !      ksize(2) = ksize(2) + 1
 !      knodes(1:nfnods,ksize(2),2) = face(nf)%nodes(1:nfnods)
 !    else if (any(abs(nodes(3,:)-0.15_wp) < eps9)) then
 !      ksize(3) = ksize(3) + 1
 !      knodes(1:nfnods,ksize(3),3) = face(nf)%nodes(1:nfnods)
 !    end if
 !  end if
 !end do
 !do k = 1,3
 !  write (nfpts,23) k,size(xyz_nodes,dim=2),ksize(k)
 !  do n = 1,ksize(k)
 !    write (nfpts,22) (knodes(l,n,k),l=1,4)
 !  end do
 !end do
 !20 format ('VARIABLES = "X", "Y", "Z"',/, &
 !           'ZONE T="All Faces", N=',i0,', E=',i0,', DATAPACKING=BLOCK', &
 !           ', ZONETYPE=FEQUADRILATERAL')
 !21 format (4es22.14)
 !22 format (*(i5))
 !23 format ('ZONE T="K',i0,' Faces", N=',i0,', E=',i0,', DATAPACKING=BLOCK', &
 !           ', ZONETYPE=FEQUADRILATERAL',/,'VARSHARELIST=([1-3]=1)')
 !close (nfpts,iostat=ierr)
 !call stop_recon(stop_mpi,pname,__LINE__,__FILE__)
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
    call mpi_allreduce(loc_vol,tot_vol,1_int_mpi,mpi_flttyp, &
                       MPI_SUM,MPI_COMM_WORLD,mpierr)
  end if
  !
  total_volume = real(tot_vol,kind=wp)
  local_volume = real(loc_vol,kind=wp)
  !
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  if (allocated(geofa2)) then
    deallocate ( geofa2 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"geofa2",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Format Statements
  !
  1 format (a,"(",i0,")%v")
  2 format (" Encountered a tetrahedral grid cell while trying to compute", &
            " the grid metrics! This grid geometry is not currently supported!")
  3 format (" Encountered a pyramid grid cell while trying to compute", &
            " the grid metrics! This grid geometry is not currently supported!")
  4 format (" Encountered a prism grid cell while trying to compute", &
            " the grid metrics! This grid geometry is not currently supported!")
  5 format (" Encountered an invalid metrics derived type!")
  6 format (" Encountered an unknown or unsupported cell geometry while", &
            " trying to compute the grid metrics!")
  !
  10 format (/,"######################################################",/, &
               "#################  FACE ",i0,"  ##########################",/, &
               "######################################################",/)
  11 format (i4,*(es18.10))
  12 format (i6,/,*(3(es18.10),/,:))
  !
end subroutine get_metrics
!
!###############################################################################
!
subroutine compute_face_normals
  !
  !.. Use Statements ..
  use order_mod, only : maxSP,maxFP,geom_solpts,n_order
  use geovar, only : nr,ncell,cell,xyz,xyz_nodes,nface
  use interpolation_mod, only : interp
  use derivatives_mod, only : deriv
  !
  !.. Local Scalars ..
  integer :: n,nc,nf,ncf,n1,n2,np
  integer :: kf,lf,l,k,nfp,ndim
  integer :: left_geom,left_order
  integer :: face_geom,face_order
  integer :: inrm,ierr,npcount,nfpts
  real(wp) :: nrm_coef
  logical(lk) :: dmat_has_been_output
  !
  !.. Local Arrays ..
  real(wp) :: dvec(1:3,1:3)
  real(wp) :: txyz(1:maxSP,1:nr)
  real(wp) :: xyz_face(1:maxFP,1:nr)
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: ipface(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_face_normals"
  character(len=*), parameter :: fname = "face_normals.tec.dat"
 !logical(lk), parameter :: output_face_normals = true
  logical(lk), parameter :: output_face_normals = fals
  !
continue
  !
 !dmat_has_been_output = fals
  !
  if (output_face_normals) then
    !
    npcount = 0
    nfpts = 0
    do nf = 1,nface
      nfpts = nfpts + size(geofa(nf)%v,dim=2)
    end do
    !
    open (newunit=inrm,file=fname,status="replace",action="write", &
                       form="formatted",position="rewind",iostat=ierr)
    call io_error(pname,fname,1,ierr)
    !
    if (nr == 2) then
      !
      write (inrm,20) nfpts,nface
      !
      allocate ( ipface(1:2,1:nface) , source=0 , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"ipface",1,__LINE__,__FILE__,ierr,error_message)
      !
    else if (nr == 3) then
      !
      write (inrm,21) nfpts,nface
      !
      allocate ( ipface(1:4,1:nface) , source=0 , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"ipface",1,__LINE__,__FILE__,ierr,error_message)
      !
    end if
    !
  end if
  !
  cell_loop: do nc = 1,ncell
    !
    left_geom  = cell(nc)%geom
    left_order = cell(nc)%order
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    np = n2-n1+1
    !
    do k = 1,np
      do l = 1,nr
        txyz(k,l) = xyz(l,k+n1-1)
      end do
    end do
   !!
   !! Compute the coordinates for the center of
   !! the current cell and save for later
   !!
   !i1 = nodes_of_cell_ptr(nc)+1
   !i2 = nodes_of_cell_ptr(nc+1)
   !ip = i2-i1+1
   !!
   !cell_center(1:nr) = sum(xyz_nodes(1:nr,nodes_of_cell(i1:i2))) / &
   !                                                   real(ip,kind=wp)
    !
    ! Loop through the faces of the current cell
    ! to get the normals at each face point
    !
    cell_face_loop: do ncf = 1,size(cell(nc)%face)
      !
      if (cell(nc)%face(ncf)%side == 2) cycle cell_face_loop
      !
      nf = cell(nc)%face(ncf)%idx
      !
      face_geom  = cell(nc)%face(ncf)%geom
      face_order = cell(nc)%face(ncf)%order
     !face_side  = cell(nc)%face(ncf)%side
      !
      nfp = geom_solpts(face_geom,face_order)
      !
      ! Interpolate the coordinates of the solution points to get the 
      ! coordinates of the flux points on the current face.
      !
      do l = 1,nr
        !
        do kf = 1,size(cell(nc)%face(ncf)%cell_idx)
          !
          lf = cell(nc)%face(ncf)%cell_idx(kf)
          !
          xyz_face(kf,l) = interp(left_geom,left_order)% &
                                     toFace(face_order)% &
                           dot( lf , txyz(1:np,l) )
          !
        end do
        !
      end do
     !do k = 1,size(cell(nc)%face(ncf)%cell_idx)
     !  write (184,'(3es12.4)') (xyz_face(k,l),l=1,nr)
     !end do
     !write (184,*)
      !
      ! If the current face is an edge, initialize the k components of the
      ! two vectors the dvec array to one and everything else to zero
      !
      if (face_geom == Geom_Edge) then
        dvec = zero
        dvec(3,:) = one
      end if
      !
      ! Now compute the derivative of the flux point coordinates
      !
      ndim = size(deriv(face_geom,face_order)%d)
      !
      do k = 1,size(cell(nc)%face(ncf)%cell_idx)
        !
        do n = 1,ndim
          !
          ! If the current cell face is a x-min, y-max, or z-min face, we need
          ! to reverse the order of the derivatives so that the cross_product
          ! computed below results in a outward facing normal for the left cell
          ! on the face.
          !
          kf = n
          if (face_geom /= Geom_Edge) then
            if (any(ncf == [Xmin3D_Face,Ymax3D_Face,Zmin3D_Face])) then
              kf = ndim+1 - n
            end if
          end if
          !
          do l = 1,nr
            dvec(l,kf) = deriv(face_geom,face_order)% &
                                                d(n)% &
                         dot( k , xyz_face(1:nfp,l) )
          end do
          !
        end do
        !
        ! Finally, the length/area normal to the current face at the current
        ! flux point is just the cross product of the two vectors in dvec
        !
       !write (183,'(3es12.4,4x,3es12.4)') (dvec(l,1),l=1,size(dvec,dim=1)), &
       !                                   (dvec(l,2),l=1,size(dvec,dim=1))
        dvec(:,3) = cross_product( dvec(:,1) , dvec(:,2) )
        !
        ! Zero/chop all the normals less than 10*roundoff
        !
        geofa(nf)%v(1:nr,k) = chop( real(dvec(1:nr,3),kind=qp)/qten ) * ten
       !geofa2(nf)%v(1:nr,k) = chop( real(dvec(1:nr,3),kind=qp)/qten ) * ten
        !
        if (output_face_normals) then
          !
          npcount = npcount + 1
          !
          if (nr == 2) then
            !
            if (k == 1) then
              ipface(1,nf) = npcount
            else if (k == size(cell(nc)%face(ncf)%cell_idx)) then
              ipface(2,nf) = npcount
            end if
            !
            write (inrm,22) (xyz_face(k,l),l=1,nr), &
                            (geofa(nf)%v(l,k),l=1,nr)!,npcount
            !
          else if (nr == 3) then
            !
            if (k == 1) then
              ipface(1,nf) = npcount
            else if (k == face_order+1) then
              ipface(2,nf) = npcount
            else if (k == size(cell(nc)%face(ncf)%cell_idx)) then
              ipface(3,nf) = npcount
            else if (face_geom == Geom_Quad) then
              if (k == (face_order+1)*face_order) then
                ipface(4,nf) = npcount
              end if
            end if
            !
            if (face_geom == Geom_Tria) ipface(4,nf) = ipface(3,nf)
            !
            write (inrm,23) (xyz_face(k,l),l=1,nr), &
                            (geofa(nf)%v(l,k),l=1,nr)!,npcount
            !
          end if
          !
        end if
        !
      end do
     !write (183,*)
      !
    end do cell_face_loop
    !
  end do cell_loop
  !
  if (output_face_normals) then
    !
    do nf = 1,nface
      write (inrm,24) (ipface(l,nf),l=1,size(ipface,dim=1))
    end do
    !
    close (inrm,iostat=ierr)
    call io_error(pname,"face_normals.tec.dat",2,ierr)
    !
    if (allocated(ipface)) then
      deallocate ( ipface , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"ipface",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
  end if
  !
  ! Format Statements
  !
  10 format (/,"######################################################",/, &
               "#################  FACE ",i0,"  ##########################",/, &
               "######################################################",/)
  11 format (2i4,/,*(3(es18.10),/,:))
 !12 format (/,"############",/, &
 !             "###  ",i0,"  ###",/, &
 !             "############",/, &
 !             *(<n_order+1>(<n_order+1>(es14.6,:),/),/))
  13 format (/,"##########################################################",/)
  !
  20 format ('VARIABLES = "X", "Y", "n_x", "n_y"',/, &
             'ZONE T="Face Normals", N=',i0,', E=',i0, &
             ', DATAPACKING=POINT, ZONETYPE=FELINESEG')
  21 format ('VARIABLES = "X", "Y", "Z", "n_x", "n_y", "n_z"',/, &
             'ZONE T="Face Normals", N=',i0,', E=',i0, &
             ', DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL')
  22 format (2es23.15,1x,2es23.15,:,1x,i0)
  23 format (3es14.6,1x,3es14.6,:,1x,i0)
  24 format (4(1x,i0))
  !
end subroutine compute_face_normals
!
!###############################################################################
!
subroutine metrics_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
  !
  !.. Formal Arguments ..
  integer, optional, intent(in) :: iunit
  !
  !.. Local Scalars ..
  integer :: i,iou
  character(len=36) :: array_name
  type(memory) :: array,dt
  type(memory) :: total
  !
continue
  !
  iou = iout
  if (present(iunit)) iou = iunit
  !
  call total%reset
  !
  write (iou,1)
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(geofa)) then
    !
    call dt%set( storage_size(geofa,kind=inttype) , &
                         size(geofa,kind=inttype) )
    !
    call array%reset
    !
    do i = lbound(geofa,dim=1),ubound(geofa,dim=1)
      !
      call array%add( storage_size(geofa(i)%v,kind=inttype) , &
                              size(geofa(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "geofa(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "geofa"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "geofa"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(geofa2)) then
    !
    call dt%set( storage_size(geofa2,kind=inttype) , &
                         size(geofa2,kind=inttype) )
    !
    call array%reset
    !
    do i = lbound(geofa2,dim=1),ubound(geofa2,dim=1)
      !
      call array%add( storage_size(geofa2(i)%v,kind=inttype) , &
                              size(geofa2(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "geofa2(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "geofa2"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "geofa2"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(dx)) then
    !
    call dt%set( storage_size(dx,kind=inttype) , &
                         size(dx,kind=inttype) )
    !
    do i = lbound(dx,dim=1),ubound(dx,dim=1)
      !
      call dt%add( storage_size(dx(i)%dr,kind=inttype) )
      call dt%add( storage_size(dx(i)%ds,kind=inttype) )
      call dt%add( storage_size(dx(i)%dt,kind=inttype) )
      call dt%add( storage_size(dx(i)%dr2,kind=inttype) )
      call dt%add( storage_size(dx(i)%ds2,kind=inttype) )
      call dt%add( storage_size(dx(i)%dt2,kind=inttype) )
      call dt%add( storage_size(dx(i)%drds,kind=inttype) )
      call dt%add( storage_size(dx(i)%drdt,kind=inttype) )
      call dt%add( storage_size(dx(i)%dsdr,kind=inttype) )
      call dt%add( storage_size(dx(i)%dsdt,kind=inttype) )
      call dt%add( storage_size(dx(i)%dtdr,kind=inttype) )
      call dt%add( storage_size(dx(i)%dtds,kind=inttype) )
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "dx"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "dx"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(dy)) then
    !
    call dt%set( storage_size(dy,kind=inttype) , &
                         size(dy,kind=inttype) )
    !
    do i = lbound(dy,dim=1),ubound(dy,dim=1)
      !
      call dt%add( storage_size(dy(i)%dr,kind=inttype) )
      call dt%add( storage_size(dy(i)%ds,kind=inttype) )
      call dt%add( storage_size(dy(i)%dt,kind=inttype) )
      call dt%add( storage_size(dy(i)%dr2,kind=inttype) )
      call dt%add( storage_size(dy(i)%ds2,kind=inttype) )
      call dt%add( storage_size(dy(i)%dt2,kind=inttype) )
      call dt%add( storage_size(dy(i)%drds,kind=inttype) )
      call dt%add( storage_size(dy(i)%drdt,kind=inttype) )
      call dt%add( storage_size(dy(i)%dsdr,kind=inttype) )
      call dt%add( storage_size(dy(i)%dsdt,kind=inttype) )
      call dt%add( storage_size(dy(i)%dtdr,kind=inttype) )
      call dt%add( storage_size(dy(i)%dtds,kind=inttype) )
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "dy"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "dy"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(dz)) then
    !
    call dt%set( storage_size(dz,kind=inttype) , &
                         size(dz,kind=inttype) )
    !
    do i = lbound(dz,dim=1),ubound(dz,dim=1)
      !
      call dt%add( storage_size(dz(i)%dr,kind=inttype) )
      call dt%add( storage_size(dz(i)%ds,kind=inttype) )
      call dt%add( storage_size(dz(i)%dt,kind=inttype) )
      call dt%add( storage_size(dz(i)%dr2,kind=inttype) )
      call dt%add( storage_size(dz(i)%ds2,kind=inttype) )
      call dt%add( storage_size(dz(i)%dt2,kind=inttype) )
      call dt%add( storage_size(dz(i)%drds,kind=inttype) )
      call dt%add( storage_size(dz(i)%drdt,kind=inttype) )
      call dt%add( storage_size(dz(i)%dsdr,kind=inttype) )
      call dt%add( storage_size(dz(i)%dsdt,kind=inttype) )
      call dt%add( storage_size(dz(i)%dtdr,kind=inttype) )
      call dt%add( storage_size(dz(i)%dtds,kind=inttype) )
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "dz"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "dz"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(dt_dJ)) then
    !
    call dt%set( storage_size(dt_dJ,kind=inttype) , &
                         size(dt_dJ,kind=inttype) )
    !
    do i = lbound(dt_dJ,dim=1),ubound(dt_dJ,dim=1)
      !
      call dt%add( storage_size(dt_dJ(i)%dr,kind=inttype) )
      call dt%add( storage_size(dt_dJ(i)%ds,kind=inttype) )
      call dt%add( storage_size(dt_dJ(i)%dt,kind=inttype) )
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "dt_dJ"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "dt_dJ"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(pdt_metrics)) then
    !
    call dt%set( storage_size(pdt_metrics,kind=inttype) , &
                         size(pdt_metrics,kind=inttype) )
    !
    select type (pdt_metrics)
      !
#ifdef DEBUG_ON
      type is (sp_metrics_3_3D)
        !
        do i = lbound(pdt_metrics,dim=1),ubound(pdt_metrics,dim=1)
          !
          call dt%add( storage_size(pdt_metrics(i)%d ,kind=inttype) , &
                               size(pdt_metrics(i)%d ,kind=inttype) )
          call dt%add( storage_size(pdt_metrics(i)%d2,kind=inttype) , &
                               size(pdt_metrics(i)%d2,kind=inttype) )
          call dt%add( storage_size(pdt_metrics(i)%d3,kind=inttype) , &
                               size(pdt_metrics(i)%d3,kind=inttype) )
          !
        end do
        !
      type is (sp_metrics_2_3D)
        !
        do i = lbound(pdt_metrics,dim=1),ubound(pdt_metrics,dim=1)
          !
          call dt%add( storage_size(pdt_metrics(i)%d ,kind=inttype) , &
                               size(pdt_metrics(i)%d ,kind=inttype) )
          call dt%add( storage_size(pdt_metrics(i)%d2,kind=inttype) , &
                               size(pdt_metrics(i)%d2,kind=inttype) )
          !
        end do
        !
      type is (sp_metrics_1_3D)
        !
        do i = lbound(pdt_metrics,dim=1),ubound(pdt_metrics,dim=1)
          !
          call dt%add( storage_size(pdt_metrics(i)%d ,kind=inttype) , &
                               size(pdt_metrics(i)%d ,kind=inttype) )
          !
        end do
        !
      type is (sp_metrics_2_2D)
        !
        do i = lbound(pdt_metrics,dim=1),ubound(pdt_metrics,dim=1)
          !
          call dt%add( storage_size(pdt_metrics(i)%d ,kind=inttype) , &
                               size(pdt_metrics(i)%d ,kind=inttype) )
          call dt%add( storage_size(pdt_metrics(i)%d2,kind=inttype) , &
                               size(pdt_metrics(i)%d2,kind=inttype) )
          !
        end do
        !
      type is (sp_metrics_1_2D)
        !
        do i = lbound(pdt_metrics,dim=1),ubound(pdt_metrics,dim=1)
          !
          call dt%add( storage_size(pdt_metrics(i)%d ,kind=inttype) , &
                               size(pdt_metrics(i)%d ,kind=inttype) )
          !
        end do
        !
#else
      type is (sp_metrics_3)
        !
        do i = lbound(pdt_metrics,dim=1),ubound(pdt_metrics,dim=1)
          !
          if (allocated(pdt_metrics(i)%d)) then
            call dt%add( storage_size(pdt_metrics(i)%d ,kind=inttype) , &
                                 size(pdt_metrics(i)%d ,kind=inttype) )
          end if
          if (allocated(pdt_metrics(i)%d2)) then
            call dt%add( storage_size(pdt_metrics(i)%d2,kind=inttype) , &
                                 size(pdt_metrics(i)%d2,kind=inttype) )
          end if
          if (allocated(pdt_metrics(i)%d3)) then
            call dt%add( storage_size(pdt_metrics(i)%d3,kind=inttype) , &
                                 size(pdt_metrics(i)%d3,kind=inttype) )
          end if
          !
        end do
        !
      type is (sp_metrics_2)
        !
        do i = lbound(pdt_metrics,dim=1),ubound(pdt_metrics,dim=1)
          !
          if (allocated(pdt_metrics(i)%d)) then
            call dt%add( storage_size(pdt_metrics(i)%d ,kind=inttype) , &
                                 size(pdt_metrics(i)%d ,kind=inttype) )
          end if
          if (allocated(pdt_metrics(i)%d2)) then
            call dt%add( storage_size(pdt_metrics(i)%d2,kind=inttype) , &
                                 size(pdt_metrics(i)%d2,kind=inttype) )
          end if
          !
        end do
        !
      type is (sp_metrics_1)
        !
        do i = lbound(pdt_metrics,dim=1),ubound(pdt_metrics,dim=1)
          !
          if (allocated(pdt_metrics(i)%d)) then
            call dt%add( storage_size(pdt_metrics(i)%d ,kind=inttype) , &
                                 size(pdt_metrics(i)%d ,kind=inttype) )
          end if
          !
        end do
        !
#endif
    end select
    !
    call total%add(dt)
    write (array_name,5) "pdt_metrics"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "pdt_metrics"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(pdt_met_dJ)) then
    !
    call dt%set( storage_size(pdt_met_dJ,kind=inttype) , &
                         size(pdt_met_dJ,kind=inttype) )
    !
    select type (pdt_met_dJ)
      !
#ifdef DEBUG_ON
      type is (sp_metrics_3_3D)
        !
        do i = lbound(pdt_met_dJ,dim=1),ubound(pdt_met_dJ,dim=1)
          !
          call dt%add( storage_size(pdt_met_dJ(i)%d ,kind=inttype) , &
                               size(pdt_met_dJ(i)%d ,kind=inttype) )
          call dt%add( storage_size(pdt_met_dJ(i)%d2,kind=inttype) , &
                               size(pdt_met_dJ(i)%d2,kind=inttype) )
          call dt%add( storage_size(pdt_met_dJ(i)%d3,kind=inttype) , &
                               size(pdt_met_dJ(i)%d3,kind=inttype) )
          !
        end do
        !
      type is (sp_metrics_2_3D)
        !
        do i = lbound(pdt_met_dJ,dim=1),ubound(pdt_met_dJ,dim=1)
          !
          call dt%add( storage_size(pdt_met_dJ(i)%d ,kind=inttype) , &
                               size(pdt_met_dJ(i)%d ,kind=inttype) )
          call dt%add( storage_size(pdt_met_dJ(i)%d2,kind=inttype) , &
                               size(pdt_met_dJ(i)%d2,kind=inttype) )
          !
        end do
        !
      type is (sp_metrics_1_3D)
        !
        do i = lbound(pdt_met_dJ,dim=1),ubound(pdt_met_dJ,dim=1)
          !
          call dt%add( storage_size(pdt_met_dJ(i)%d ,kind=inttype) , &
                               size(pdt_met_dJ(i)%d ,kind=inttype) )
          !
        end do
        !
      type is (sp_metrics_2_2D)
        !
        do i = lbound(pdt_met_dJ,dim=1),ubound(pdt_met_dJ,dim=1)
          !
          call dt%add( storage_size(pdt_met_dJ(i)%d ,kind=inttype) , &
                               size(pdt_met_dJ(i)%d ,kind=inttype) )
          call dt%add( storage_size(pdt_met_dJ(i)%d2,kind=inttype) , &
                               size(pdt_met_dJ(i)%d2,kind=inttype) )
          !
        end do
        !
      type is (sp_metrics_1_2D)
        !
        do i = lbound(pdt_met_dJ,dim=1),ubound(pdt_met_dJ,dim=1)
          !
          call dt%add( storage_size(pdt_met_dJ(i)%d ,kind=inttype) , &
                               size(pdt_met_dJ(i)%d ,kind=inttype) )
          !
        end do
        !
#else
      type is (sp_metrics_3)
        !
        do i = lbound(pdt_met_dJ,dim=1),ubound(pdt_met_dJ,dim=1)
          !
          if (allocated(pdt_met_dJ(i)%d)) then
            call dt%add( storage_size(pdt_met_dJ(i)%d ,kind=inttype) , &
                                 size(pdt_met_dJ(i)%d ,kind=inttype) )
          end if
          if (allocated(pdt_met_dJ(i)%d2)) then
            call dt%add( storage_size(pdt_met_dJ(i)%d2,kind=inttype) , &
                                 size(pdt_met_dJ(i)%d2,kind=inttype) )
          end if
          if (allocated(pdt_met_dJ(i)%d3)) then
            call dt%add( storage_size(pdt_met_dJ(i)%d3,kind=inttype) , &
                                 size(pdt_met_dJ(i)%d3,kind=inttype) )
          end if
          !
        end do
        !
      type is (sp_metrics_2)
        !
        do i = lbound(pdt_met_dJ,dim=1),ubound(pdt_met_dJ,dim=1)
          !
          if (allocated(pdt_met_dJ(i)%d)) then
            call dt%add( storage_size(pdt_met_dJ(i)%d ,kind=inttype) , &
                                 size(pdt_met_dJ(i)%d ,kind=inttype) )
          end if
          if (allocated(pdt_met_dJ(i)%d2)) then
            call dt%add( storage_size(pdt_met_dJ(i)%d2,kind=inttype) , &
                                 size(pdt_met_dJ(i)%d2,kind=inttype) )
          end if
          !
        end do
        !
      type is (sp_metrics_1)
        !
        do i = lbound(pdt_met_dJ,dim=1),ubound(pdt_met_dJ,dim=1)
          !
          if (allocated(pdt_met_dJ(i)%d)) then
            call dt%add( storage_size(pdt_met_dJ(i)%d ,kind=inttype) , &
                                 size(pdt_met_dJ(i)%d ,kind=inttype) )
          end if
          !
        end do
        !
#endif
    end select
    !
    call total%add(dt)
    write (array_name,5) "pdt_met_dJ"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "pdt_met_dJ"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(metrics)) then
    call array%set( storage_size(metrics,kind=inttype) , &
                            size(metrics,kind=inttype) )
    call total%add(array)
    write (array_name,5) "metrics"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "metrics"
  end if
  if (allocated(metrics2)) then
    call array%set( storage_size(metrics2,kind=inttype) , &
                            size(metrics2,kind=inttype) )
    call total%add(array)
    write (array_name,5) "metrics2"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "metrics2"
  end if
  if (allocated(metrics3)) then
    call array%set( storage_size(metrics3,kind=inttype) , &
                            size(metrics3,kind=inttype) )
    call total%add(array)
    write (array_name,5) "metrics3"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "metrics3"
  end if
  if (allocated(jacobian)) then
    call array%set( storage_size(jacobian,kind=inttype) , &
                            size(jacobian,kind=inttype) )
    call total%add(array)
    write (array_name,5) "jacobian"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "jacobian"
  end if
  if (allocated(dJ)) then
    call array%set( storage_size(dJ,kind=inttype) , &
                            size(dJ,kind=inttype) )
    call total%add(array)
    write (array_name,5) "dJ"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "dJ"
  end if
  if (allocated(volj)) then
    call array%set( storage_size(volj,kind=inttype) , &
                            size(volj,kind=inttype) )
    call total%add(array)
    write (array_name,5) "volj"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "volj"
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module metrics_mod")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  5 format ("'",a,"'")
  6 format (" Array Component: ",a,dt)
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
end subroutine metrics_memory_usage
!
!###############################################################################
!
end module metrics_mod
