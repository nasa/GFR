module flowvar
  !
  use module_kind_types, only : r8,wp
  use module_kind_types, only : zero
  use module_kind_types, only : Classic_RK
  use module_kind_types, only : CK4_RK
  use module_kind_types, only : TVD4_RK
  use module_kind_types, only : alloc_error
  use module_kind_types, only : bc_cpu_bnd
  !
  implicit none
  !
  private
  !
  public :: flowvar_alloc_status
  public :: allocate_flowvar_arrays
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
  ! TROUBLED_CELL : Mask/Marker array identifying troubled cells
  !
  logical, public, save, allocatable, dimension(:) :: troubled_cell
  !
  !
  !
  ! FACE_XYZ : Derived type containing a set of variables at the flux/face
  !            points on each side of a face between cells
  !
  type, public :: face_xyz
    real(wp), allocatable :: v(:,:)
  end type face_xyz
  !
  !
  ! FACE_VAR : Derived type containing a set of variables at the flux/face
  !            points on each side of a face between cells
  !
  type, public :: face_var
    real(wp), allocatable :: v(:,:,:)
  contains
    procedure :: zero => zero_face_var
    procedure :: pack => pack_face_var_data
    procedure :: unpack => unpack_face_var_data
  end type face_var
  !
  !
  ! FACE_VEC : Derived type containing a vector of a set of variables at the
  !            flux/face points on each side of a face between cells. This is
  !            basically the same thing as face_var but with an additional
  !            rank for v to allow for directional quantities for each
  !            variable in the set.
  !
  type, public :: face_vec
    real(wp), allocatable :: v(:,:,:,:)
  contains
    procedure :: zero => zero_face_vec
    procedure :: pack => pack_face_vec_data
    procedure :: unpack => unpack_face_vec_data
  end type face_vec
  !
  ! FACEFLX  : flux normal to cell faces at each interface flux point
  ! FACEUSP  : conservative variables at the interface flux points
  ! FACEDUSP : conservative variable gradients at the interface flux points
  ! FACEXYZ  : coordinates of the flux points on the boundary faces (used only for MMS)
  !
  ! faceflx(1:nface)%v(1:nq,1:n_sp1de,1:2) ! this is the normal face flux
  ! faceusp(1:nface)%v(1:nq,1:n_sp1de,1:2)
  ! facedusp(1:nface)%v(1:2,1:nq,1:n_sp1de,1:2)
  !
  type(face_var), public, save, allocatable, dimension(:) :: faceusp
  type(face_vec), public, save, allocatable, dimension(:) :: facedusp
  type(face_var), public, save, allocatable, dimension(:) :: faceflx
  type(face_xyz), public, save, allocatable, dimension(:) :: facexyz
  !
contains
!
!###############################################################################
!
subroutine allocate_flowvar_arrays()
  !
  !.. Use Statements ..
  use eqn_idx,   only : nq
  use geovar,    only : nr,n_solpts,n_totpts,nface,nfbnd,face
  use geovar,    only : n_global_totpts
  use order_mod, only : geom_solpts
  use order_mod, only : n_fp1de
  use ovar,      only : iRKtype,num_stages,mms_opt
  !
  !.. Local Scalars ..
  integer :: n,ierr,face_geom,face_order,nfpts,nf
  character(len=200) :: array_name
  character(len=200) :: error_message
  !
  !.. Local Constants ..
  character(len=*), parameter :: pname = "allocate_flowvar_arrays"
  !
continue
  !
  allocate ( usp(1:nq,1:n_totpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"usp",1,ierr,error_message)
  !
  !
  allocate ( dusp(1:nr,1:nq,1:n_totpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"dusp",1,ierr,error_message)
  !
  !
  allocate ( residual(1:nq,1:n_solpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"residual",1,ierr,error_message)
  !
  !
  allocate ( dtsp(1:n_solpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"dtsp",1,ierr,error_message)
  !
  !
  allocate ( uoldsp(1:nq,1:n_solpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"uoldsp",1,ierr,error_message)
  !
  !
  ! Allocate the outer part of the face solution/flux arrays
  !
  allocate ( faceusp(1:nface) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faceusp",1,ierr,error_message)
  !
  allocate ( facedusp(1:nface) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"facedusp",1,ierr,error_message)
  !
  allocate ( faceflx(1:nface) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faceflx",1,ierr,error_message)
  !
  ! If using MMS, allocate the facexyz array
  !
  if (mms_opt /= 0) then
    allocate ( facexyz(1:nfbnd) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facexyz",1,ierr,error_message)
  end if
  !
  ! Allocate the data component of these arrays based on the face geometry
  !
  do nf = 1,nface
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    nfpts = geom_solpts(face_geom,face_order)
    !
    allocate ( faceusp(nf)%v(1:nq,1:nfpts,1:2) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "faceusp",nf
    call alloc_error(pname,array_name,1,ierr,error_message)
    !
    allocate ( facedusp(nf)%v(1:nr,1:nq,1:nfpts,1:2) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "facedusp",nf
    call alloc_error(pname,array_name,1,ierr,error_message)
    !
    allocate ( faceflx(nf)%v(1:nq,1:nfpts,1:2) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "faceflx",nf
    call alloc_error(pname,array_name,1,ierr,error_message)
    !
    ! If using MMS and this is a boundary face,
    ! allocate the v component of facexyz
    !
    if (mms_opt /= 0 .and. nf <= nfbnd) then
      allocate ( facexyz(nf)%v(1:nr,1:nfpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) "facexyz",nf
      call alloc_error(pname,array_name,1,ierr,error_message)
    end if
    !
  end do
  !
  ! Allocate the storage for the Runge-Kutta methods
  !
  n = 0
  if (iRKtype == Classic_RK) then
    n = num_stages
  else if (iRKtype == CK4_RK .or. iRKtype == TVD4_RK) then
    n = 1
    if (iRKtype == TVD4_RK) then
      !
      allocate ( usp_sum(1:nq,1:n_solpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"usp_sum",1,ierr,error_message)
      !
    end if
  end if
  !
  if (n /= 0) then
    !
    allocate ( rssprk(1:nq,1:n_solpts,1:n) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"rssprk",1,ierr,error_message)
    !
  end if
  !
  ! Allocate the global_usp array
  !
  allocate ( global_usp(1:nq,1:n_global_totpts) , source=zero , stat=ierr )
  call alloc_error(pname,"global_usp",1,ierr)
  !
  ! Format Statements
  !
  1 format (a,"(",i0,")%v")
  !
end subroutine allocate_flowvar_arrays
!
!###############################################################################
!
elemental subroutine zero_face_var(this)
  !
  !.. Formal Arguments ..
  class(face_var), intent(inout) :: this
  !
continue
  !
  if (allocated(this%v)) this%v = zero
  !
end subroutine zero_face_var
!
!###############################################################################
!
elemental subroutine zero_face_vec(this)
  !
  !.. Formal Arguments ..
  class(face_vec), intent(inout) :: this
  !
continue
  !
  if (allocated(this%v)) this%v = zero
  !
end subroutine zero_face_vec
!
!###############################################################################
!
subroutine pack_face_var_data(this,package)
  !
  !.. Use Statements ..
  use geovar, only : nfbnd,bface
  !
  !.. Formal Arguments ..
  class(face_var), dimension(:),    intent(in) :: this
  real(r8),        dimension(:), intent(inout) :: package
  !
  !.. Local Scalars ..
  integer :: nf,ibeg,iend,nsz
  !
continue
  !
  iend = 0
  !
  do nf = 1,nfbnd
    !
    if (bface(1,nf) /= bc_cpu_bnd) cycle
    !
    nsz = size(this(nf)%v(:,:,1))
    !
    ibeg = iend + 1
    iend = iend + nsz
    !
    package(ibeg:iend) = reshape( real(this(nf)%v(:,:,1),kind=r8) , [nsz] )
    !
  end do
  !
end subroutine pack_face_var_data
!
!###############################################################################
!
subroutine pack_face_vec_data(this,package)
  !
  !.. Use Statements ..
  use geovar, only : nfbnd,bface
  !
  !.. Formal Arguments ..
  class(face_vec), dimension(:),    intent(in) :: this
  real(r8),        dimension(:), intent(inout) :: package
  !
  !.. Local Scalars ..
  integer :: nf,ibeg,iend,nsz
  !
continue
  !
  iend = 0
  !
  do nf = 1,nfbnd
    !
    if (bface(1,nf) /= bc_cpu_bnd) cycle
    !
    nsz = size(this(nf)%v(:,:,:,1))
    !
    ibeg = iend + 1
    iend = iend + nsz
    !
    package(ibeg:iend) = reshape( real(this(nf)%v(:,:,:,1),kind=r8) , [nsz] )
    !
  end do
  !
end subroutine pack_face_vec_data
!
!###############################################################################
!
subroutine unpack_face_var_data(this,package)
  !
  !.. Use Statements ..
  use geovar, only : nfbnd,bface
  !
  !.. Formal Arguments ..
  class(face_var), dimension(:), intent(inout) :: this
  real(r8),        dimension(:),    intent(in) :: package
  !
  !.. Local Scalars ..
  integer :: nf,ibeg,iend
  !
  !.. Local Arrays ..
  integer, dimension(1:2) :: nsz
  !
continue
  !
  iend = 0
  !
  do nf = 1,nfbnd
    !
    if (bface(1,nf) /= bc_cpu_bnd) cycle
    !
    nsz(1) = size(this(nf)%v,dim=1)
    nsz(2) = size(this(nf)%v,dim=2)
    !
    ibeg = iend + 1
    iend = iend + product(nsz)
    !
    this(nf)%v(:,:,2) = reshape( real(package(ibeg:iend),kind=wp) , nsz )
    !
  end do
  !
end subroutine unpack_face_var_data
!
!###############################################################################
!
subroutine unpack_face_vec_data(this,package)
  !
  !.. Use Statements ..
  use geovar, only : nfbnd,bface
  !
  !.. Formal Arguments ..
  class(face_vec), dimension(:), intent(inout) :: this
  real(r8),        dimension(:),    intent(in) :: package
  !
  !.. Local Scalars ..
  integer :: nf,ibeg,iend
  !
  !.. Local Arrays ..
  integer, dimension(1:3) :: nsz
  !
continue
  !
  iend = 0
  !
  do nf = 1,nfbnd
    !
    if (bface(1,nf) /= bc_cpu_bnd) cycle
    !
    nsz(1) = size(this(nf)%v,dim=1)
    nsz(2) = size(this(nf)%v,dim=2)
    nsz(3) = size(this(nf)%v,dim=3)
    !
    ibeg = iend + 1
    iend = iend + product(nsz)
    !
    this(nf)%v(:,:,:,2) = reshape( real(package(ibeg:iend),kind=wp) , nsz )
    !
  end do
  !
end subroutine unpack_face_vec_data
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
