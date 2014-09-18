module flowvar
  !
  use module_kind_types, only : wp
  !
  implicit none
  !
  private
  !
  ! USP  : conservative variables at the solution points (@SP)
  ! VSP  : primitive variables @SP
  ! DUSP : gradient of primitive variables @SP
  !
  real(wp), public, save, allocatable, dimension(:,:)   :: usp
  real(wp), public, save, allocatable, dimension(:,:)   :: vsp
  real(wp), public, save, allocatable, dimension(:,:,:) :: dusp
  !
  ! UOLDSP   : conservative variables of previous time step @SP 
  ! RESIDUAL : RHS residual @SP for current time/RK step
  ! 
  real(wp), public, save, allocatable, dimension(:,:) :: uoldsp
  real(wp), public, save, allocatable, dimension(:,:) :: residual
  !
  ! GLOBAL_USP : conservative variables @SP across all processors
  ! GLOBAL_VSP : primitive    variables @SP across all processors
  !
  real(wp), public, save, allocatable, dimension(:,:) :: global_usp
  real(wp), public, save, allocatable, dimension(:,:) :: global_vsp
  !
  ! FACEFLX  : flux normal to cell faces at each interface flux point
  ! FACEUSP  : conservative variables at the interface flux points
  ! FACEDUSP : conservative variable gradients at the interface flux points
  !
  ! faceflx(1:neqns,1:2,1:n_sp1de,1:nface) ! this is the normal face flux
  ! faceusp(1:neqns,1:2,1:n_sp1de,1:nface)
  ! facedusp(1:2,1:neqns,1:2,1:n_sp1de,1:nface)
  !
  real(wp), public, save, allocatable, dimension(:,:,:,:)   :: faceusp
  real(wp), public, save, allocatable, dimension(:,:,:,:,:) :: facedusp
  real(wp), public, save, allocatable, dimension(:,:,:,:)   :: faceflx
  !
  ! USP_SUM : auxilary array needed for SSP-RK4 method
  ! RSSPRK  : auxilary array containing old residuals, needed for some RK methods
  !
  real(wp), public, save, allocatable, dimension(:,:)   :: usp_sum
  real(wp), public, save, allocatable, dimension(:,:,:) :: rssprk
  !
  ! DTSP : Time step for each solution point
  !
  real(wp), public, save, allocatable, dimension(:) :: dtsp
  !
  ! PROFILE_QUANTITY: Derived type containing the name of a profile
  !                   quantity and the corresponding data
  !
  type, public :: profile_quantity
    character(len=100) :: qty_name
    real(wp), allocatable, dimension(:) :: qty_data
  end type profile_quantity
  !
  ! PROFILE: Derived type of a single profile containing multiple quantities
  !
  type, public :: profile
    character(len=100) :: prof_name
    integer :: group_num
    integer, allocatable, dimension(:) :: faces
    type(profile_quantity), allocatable, dimension(:) :: qty
  end type profile
  !
  ! BND_PROF : Solution profiles along each boundary
  ! SOL_PROF : Solution profiles within the interior
  !
  type(profile), public, save, allocatable, dimension(:) :: bnd_prof
  type(profile), public, save, allocatable, dimension(:) :: sol_prof
  !
end module flowvar
!
!###############################################################################
!
