module flowvar
  !
  use module_kind_types, only : wp
  !
  implicit none
  !
  private
  !
  public :: flowvar_alloc_status
  !
  ! USP  : conservative variables at the solution points (@SP)
  ! DUSP : gradient of primitive variables @SP
  !
  real(wp), public, save, allocatable, dimension(:,:)   :: usp
  real(wp), public, save, allocatable, dimension(:,:,:) :: dusp
  !
  ! UOLDSP   : conservative variables of previous time step @SP 
  ! RESIDUAL : RHS residual @SP for current time/RK step
  ! 
  real(wp), public, save, allocatable, dimension(:,:) :: uoldsp
  real(wp), public, save, allocatable, dimension(:,:) :: residual
  !
  ! GLOBAL_USP : conservative variables @SP across all processors
  ! GLOBAL_VAR : arbitrary    variable  @SP across all processors
  !
  real(wp), public, save, allocatable, dimension(:,:) :: global_usp
  real(wp), public, save, allocatable, dimension(:)   :: global_var
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
contains
!
!###############################################################################
!
subroutine flowvar_alloc_status(iunit)
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
  total%mem = 0
  !
  write (iunit,1)
  !
  if (allocated(usp)) then
    array%siz = size(usp)
    array%mem = array%siz*storage_size(usp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'usp'",array%get_units()
  end if
  if (allocated(dusp)) then
    array%siz = size(dusp)
    array%mem = array%siz*storage_size(dusp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'dusp'",array%get_units()
  end if
  if (allocated(uoldsp)) then
    array%siz = size(uoldsp)
    array%mem = array%siz*storage_size(uoldsp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'uoldsp'",array%get_units()
  end if
  if (allocated(residual)) then
    array%siz = size(residual)
    array%mem = array%siz*storage_size(residual)
    total%mem = total%mem + array%mem
    write (iunit,2) "'residual'",array%get_units()
  end if
  if (allocated(global_usp)) then
    array%siz = size(global_usp)
    array%mem = array%siz*storage_size(global_usp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'global_usp'",array%get_units()
  end if
  if (allocated(global_var)) then
    array%siz = size(global_var)
    array%mem = array%siz*storage_size(global_var)
    total%mem = total%mem + array%mem
    write (iunit,2) "'global_var'",array%get_units()
  end if
  if (allocated(faceusp)) then
    array%siz = size(faceusp)
    array%mem = array%siz*storage_size(faceusp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'faceusp'",array%get_units()
  end if
  if (allocated(facedusp)) then
    array%siz = size(facedusp)
    array%mem = array%siz*storage_size(facedusp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'facedusp'",array%get_units()
  end if
  if (allocated(faceflx)) then
    array%siz = size(faceflx)
    array%mem = array%siz*storage_size(faceflx)
    total%mem = total%mem + array%mem
    write (iunit,2) "'faceflx'",array%get_units()
  end if
  if (allocated(usp_sum)) then
    array%siz = size(usp_sum)
    array%mem = array%siz*storage_size(usp_sum)
    total%mem = total%mem + array%mem
    write (iunit,2) "'usp_sum'",array%get_units()
  end if
  if (allocated(rssprk)) then
    array%siz = size(rssprk)
    array%mem = array%siz*storage_size(rssprk)
    total%mem = total%mem + array%mem
    write (iunit,2) "'rssprk'",array%get_units()
  end if
  if (allocated(dtsp)) then
    array%siz = size(dtsp)
    array%mem = array%siz*storage_size(dtsp)
    total%mem = total%mem + array%mem
    write (iunit,2) "'dtsp'",array%get_units()
  end if
  if (allocated(bnd_prof)) then
    array%siz = size(bnd_prof)
    array%mem = array%siz*storage_size(bnd_prof)
    total%mem = total%mem + array%mem
    write (iunit,2) "'bnd_prof'",array%get_units()
  end if
  if (allocated(sol_prof)) then
    array%siz = size(sol_prof)
    array%mem = array%siz*storage_size(sol_prof)
    total%mem = total%mem + array%mem
    write (iunit,2) "'sol_prof'",array%get_units()
  end if
  !
  total = total%get_units()
  write (iunit,3) total%mem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module flowvar")
  2 format (" Array: ",a12,dt)
  3 format (" Total memory = ",i0,1x,a)
  !
end subroutine flowvar_alloc_status
!
!###############################################################################
!
end module flowvar
!
!###############################################################################
!
