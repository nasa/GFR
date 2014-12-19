module metvar
  !
  use module_kind_types, only : wp
  !
  implicit none
  !
  private
  !
  public :: metvar_alloc_status
  !
  !
  ! geofa : normal of the face at each of the solution points on the face
  !
  real(wp), public, save, allocatable, dimension(:,:,:) :: geofa
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
  !     For solution point 'n', metrics(l,i,n) gives the partial derivative
  !                             d (x_l)
  !                             -------
  !                             d(xi_i)
  !     where ( x_1, x_2, x_3) = (x,y,z) -> the physical coordinates
  !     where (xi_1,xi_2,xi_3) = (xi,eta,zeta) or (r,s,t) -> the computational
  !                           coordinates within the standard elemental region
  !           NOTE: For 2D z and zeta (or t) do not exist
  !
  real(wp), public, save, allocatable, dimension(:,:,:) :: metrics
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
  ! cormat : correction matrix
  ! spcemat : curved edge CE (at bottom), solution points
  ! spceximat : CE metrics d_xi
  ! spceetmat : CE metrics d_eta
  ! spcexi2mat : CE metrics d_xi2
  ! spceet2mat : CE metrics d_eta2
  ! spcexietmat : CE metrics d_xieta
  !
  ! Hesthaven and Warburton equivalencies
  ! spmat  = rs
  ! dximat = Dr
  ! detmat = Ds
  ! cormat = LIFT
  !
  real(wp), public, save, allocatable, dimension(:,:,:) :: dr_mat
  real(wp), public, save, allocatable, dimension(:,:,:) :: ds_mat
  real(wp), public, save, allocatable, dimension(:,:,:) :: cormat
  !
  ! lagrange : gives the Lagrange polynomial coefficients to project
  !            variables/fluxes from the solution points in the cell
  !            interior to the flux points on the cell faces
  !
  real(wp), public, save, allocatable, dimension(:,:,:) :: dqr_LPmat
  real(wp), public, save, allocatable, dimension(:,:,:) :: dqs_LPmat
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
subroutine metvar_alloc_status(iunit)
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
  if (allocated(dr_mat)) then
    array%siz = size(dr_mat)
    array%mem = array%siz*storage_size(dr_mat)
    total%mem = total%mem + array%mem
    write (iunit,2) "'dr_mat'",array%get_units()
  end if
  if (allocated(ds_mat)) then
    array%siz = size(ds_mat)
    array%mem = array%siz*storage_size(ds_mat)
    total%mem = total%mem + array%mem
    write (iunit,2) "'ds_mat'",array%get_units()
  end if
  if (allocated(cormat)) then
    array%siz = size(cormat)
    array%mem = array%siz*storage_size(cormat)
    total%mem = total%mem + array%mem
    write (iunit,2) "'cormat'",array%get_units()
  end if
  if (allocated(dqr_LPmat)) then
    array%siz = size(dqr_LPmat)
    array%mem = array%siz*storage_size(dqr_LPmat)
    total%mem = total%mem + array%mem
    write (iunit,2) "'dqr_LPmat'",array%get_units()
  end if
  if (allocated(dqs_LPmat)) then
    array%siz = size(dqs_LPmat)
    array%mem = array%siz*storage_size(dqs_LPmat)
    total%mem = total%mem + array%mem
    write (iunit,2) "'dqs_LPmat'",array%get_units()
  end if
  if (allocated(dx)) then
    array%siz = size(dx)
    array%mem = array%siz*storage_size(dx)
    total%mem = total%mem + array%mem
    write (iunit,2) "'dx'",array%get_units()
  end if
  if (allocated(dy)) then
    array%siz = size(dy)
    array%mem = array%siz*storage_size(dy)
    total%mem = total%mem + array%mem
    write (iunit,2) "'dy'",array%get_units()
  end if
  if (allocated(geofa)) then
    array%siz = size(geofa)
    array%mem = array%siz*storage_size(geofa)
    total%mem = total%mem + array%mem
    write (iunit,2) "'geofa'",array%get_units()
  end if
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
  1 format (/,"Status of arrays in module metvar")
  2 format (" Array: ",a11,dt)
  3 format (" Total memory = ",i0,1x,a)
  !
end subroutine metvar_alloc_status
!
!###############################################################################
!
end module metvar
!
!###############################################################################
!
