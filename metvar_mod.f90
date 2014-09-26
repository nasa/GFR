module metvar
  !
  use module_kind_types, only : wp
  !
  implicit none
  !
  private
  !
  type, public :: derivative
    real(wp) :: dr
    real(wp) :: ds
  end type derivative
  !
  ! NEED NEWER COMPILER FOR THIS
  !
 !type, extends(derivative) :: second_derivative
 !  real(wp) :: dr2
 !  real(wp) :: ds2
 !  real(wp) :: drds
 !end type second_derivative
  !
  ! FOR OLDER COMPILER
  !
  type, public :: second_derivative
    real(wp) :: dr
    real(wp) :: ds
    real(wp) :: dr2
    real(wp) :: ds2
    real(wp) :: drds
  end type second_derivative
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
  real(wp), public, save, allocatable, dimension(:,:)   :: d1_mat
  real(wp), public, save, allocatable, dimension(:,:,:) :: dr_mat
  real(wp), public, save, allocatable, dimension(:,:,:) :: ds_mat
  real(wp), public, save, allocatable, dimension(:,:,:) :: dr_hat
  real(wp), public, save, allocatable, dimension(:,:,:) :: ds_hat
  real(wp), public, save, allocatable, dimension(:,:,:) :: cormat
  real(wp), public, save, allocatable, dimension(:,:,:) :: cormat_dg
  !
  ! lagrange : gives the Lagrange polynomial coefficients to project
  !            variables/fluxes from the solution points in the cell
  !            interior to the flux points on the cell faces
  !
  real(wp), public, save, allocatable, dimension(:,:)   :: dq1D_mat
  real(wp), public, save, allocatable, dimension(:,:,:) :: dqr_mat
  real(wp), public, save, allocatable, dimension(:,:,:) :: dqs_mat
  real(wp), public, save, allocatable, dimension(:,:,:) :: dqr_LPmat
  real(wp), public, save, allocatable, dimension(:,:,:) :: dqs_LPmat
  real(wp), public, save, allocatable, dimension(:,:,:) :: dqr_CRmat
  real(wp), public, save, allocatable, dimension(:,:,:) :: dqs_CRmat
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
  type(second_derivative), public, save, allocatable, dimension(:) :: dx
  type(second_derivative), public, save, allocatable, dimension(:) :: dy
  !
  ! geofa : normal of the face at each of the solution points on the face
  !
  real(wp), public, save, allocatable, dimension(:,:,:) :: geofa
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
  ! jacobian : jacobian of the solution points
  ! dJdr     : d_jacobian/d_r of the solution points
  ! dJds     : d_jacobian/d_s of the solution points
  !
  real(wp), public, save, allocatable, dimension(:) :: jacobian
 !real(wp), public, save, allocatable, dimension(:) :: dJdr
 !real(wp), public, save, allocatable, dimension(:) :: dJds
  type(derivative), public, save, allocatable, dimension(:) :: dJ
  !
  ! volj : Volume of the grid cells evaluated using the jacobian of the solution points
  !
  real(wp), public, save, allocatable, dimension(:) :: volj  !(-nfbnd:ncell)
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
end module metvar
!
!###############################################################################
!
