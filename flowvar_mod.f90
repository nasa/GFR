module flowvar
  !
  use module_kind_types, only : r8,wp
  use module_kind_types, only : zero,iout
  use module_kind_types, only : Classic_RK
  use module_kind_types, only : CK4_RK
  use module_kind_types, only : TVD4_RK
  use module_kind_types, only : alloc_error
  use module_kind_types, only : bc_cpu_bnd
  use module_kind_types, only : skip_alloc_pause
  use module_kind_types, only : error_message
  !
  implicit none
  !
  private
  !
  public :: flowvar_memory_usage
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
#ifdef USE_TYPEBOUND_FACE_PACK
    ! THESE PACK AND UNPACK PROCEDURES WORK WITH INTEL 2015.3.187 BUT
    ! ARE NOT ALLOWED BY THE FORTRAN STANDARD. THEY WILL PRODUCE A
    ! COMPILE-TIME ERROR IF USED WITH 2016.0.109
    procedure :: pack => pack_face_var_data
    procedure :: unpack => unpack_face_var_data
#endif
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
#ifdef USE_TYPEBOUND_FACE_PACK
    ! THESE PACK AND UNPACK PROCEDURES WORK WITH INTEL 2015.3.187 BUT
    ! ARE NOT ALLOWED BY THE FORTRAN STANDARD. THEY WILL PRODUCE A
    ! COMPILE-TIME ERROR IF USED WITH 2016.0.109
    procedure :: pack => pack_face_vec_data
    procedure :: unpack => unpack_face_vec_data
#endif
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
  use order_mod, only : geom_solpts
  use order_mod, only : n_fp1de
  use ovar,      only : iRKtype,num_stages,mms_opt
  use ovar,      only : itestcase,num_ts
  use ovar,      only : convert_restart_to_cgns
  !
  !.. Local Scalars ..
  integer :: n,ierr,face_geom,face_order,nfpts,nf
  character(len=200) :: array_name
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
  if (.not. convert_restart_to_cgns) then
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
      call alloc_error(pname,array_name,1,ierr,error_message,skip_alloc_pause)
      !
      allocate ( facedusp(nf)%v(1:nr,1:nq,1:nfpts,1:2) , source=zero , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) "facedusp",nf
      call alloc_error(pname,array_name,1,ierr,error_message,skip_alloc_pause)
      !
      allocate ( faceflx(nf)%v(1:nq,1:nfpts,1:2) , source=zero , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) "faceflx",nf
      call alloc_error(pname,array_name,1,ierr,error_message,skip_alloc_pause)
      !
      ! If using MMS and this is a boundary face,
      ! allocate the v component of facexyz
      !
      if (mms_opt /= 0 .and. nf <= nfbnd) then
        allocate ( facexyz(nf)%v(1:nr,1:nfpts) , source=zero , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) "facexyz",nf
        call alloc_error(pname,array_name,1,ierr,error_message,skip_alloc_pause)
      end if
      !
    end do
    !
  end if
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
subroutine flowvar_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
  !
  !.. Formal Arguments ..
  integer, optional, intent(in) :: iunit
  !
  !.. Local Scalars ..
  integer :: iou,i,j,k
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
  if (allocated(faceusp)) then
    !
    write (iou,8) "faceusp",storage_size(faceusp),storage_size(faceusp)/8
    call dt%set( storage_size(faceusp,kind=inttype) , &
                         size(faceusp,kind=inttype) )
    !
    call array%reset
    !
    do i = lbound(faceusp,dim=1),ubound(faceusp,dim=1)
      !
      call array%add( storage_size(faceusp(i)%v,kind=inttype) , &
                              size(faceusp(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "faceusp(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "faceusp"
    write (iou,2) array_name,dt%get_units()
    !
  else
    write (iou,7) "faceusp"
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(facedusp)) then
    !
    write (iou,8) "facedusp",storage_size(facedusp),storage_size(facedusp)/8
    call dt%set( storage_size(facedusp,kind=inttype) , &
                         size(facedusp,kind=inttype) )
    !
    call array%reset
    !
    do i = lbound(facedusp,dim=1),ubound(facedusp,dim=1)
      !
      call array%add( storage_size(facedusp(i)%v,kind=inttype) , &
                              size(facedusp(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "facedusp(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "facedusp"
    write (iou,2) array_name,dt%get_units()
    !
  else
    write (iou,7) "facedusp"
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(faceflx)) then
    !
    write (iou,8) "faceflx",storage_size(faceflx),storage_size(faceflx)/8
    call dt%set( storage_size(faceflx,kind=inttype) , &
                         size(faceflx,kind=inttype) )
    !
    call array%reset
    !
    do i = lbound(faceflx,dim=1),ubound(faceflx,dim=1)
      !
      call array%add( storage_size(faceflx(i)%v,kind=inttype) , &
                              size(faceflx(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "faceflx(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "faceflx"
    write (iou,2) array_name,dt%get_units()
    !
  else
    write (iou,7) "faceflx"
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(facexyz)) then
    !
    write (iou,8) "facexyz",storage_size(facexyz),storage_size(facexyz)/8
    call dt%set( storage_size(facexyz,kind=inttype) , &
                         size(facexyz,kind=inttype) )
    !
    call array%reset
    !
    do i = lbound(facexyz,dim=1),ubound(facexyz,dim=1)
      !
      call array%add( storage_size(facexyz(i)%v,kind=inttype) , &
                              size(facexyz(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "facexyz(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "facexyz"
    write (iou,2) array_name,dt%get_units()
    !
  else
    write (iou,7) "facexyz"
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(usp)) then
    call array%set( storage_size(usp,kind=inttype) , &
                            size(usp,kind=inttype) )
    call total%add(array)
    write (array_name,5) "usp"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "usp"
  end if
  if (allocated(dusp)) then
    call array%set( storage_size(dusp,kind=inttype) , &
                            size(dusp,kind=inttype) )
    call total%add(array)
    write (array_name,5) "dusp"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "dusp"
  end if
  if (allocated(uoldsp)) then
    call array%set( storage_size(uoldsp,kind=inttype) , &
                            size(uoldsp,kind=inttype) )
    call total%add(array)
    write (array_name,5) "uoldsp"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "uoldsp"
  end if
  if (allocated(residual)) then
    call array%set( storage_size(residual,kind=inttype) , &
                            size(residual,kind=inttype) )
    call total%add(array)
    write (array_name,5) "residual"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "residual"
  end if
  if (allocated(global_usp)) then
    call array%set( storage_size(global_usp,kind=inttype) , &
                            size(global_usp,kind=inttype) )
    call total%add(array)
    write (array_name,5) "global_usp"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "global_usp"
  end if
  if (allocated(global_var)) then
    call array%set( storage_size(global_var,kind=inttype) , &
                            size(global_var,kind=inttype) )
    call total%add(array)
    write (array_name,5) "global_var"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "global_var"
  end if
  if (allocated(usp_sum)) then
    call array%set( storage_size(usp_sum,kind=inttype) , &
                            size(usp_sum,kind=inttype) )
    call total%add(array)
    write (array_name,5) "usp_sum"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "usp_sum"
  end if
  if (allocated(rssprk)) then
    call array%set( storage_size(rssprk,kind=inttype) , &
                            size(rssprk,kind=inttype) )
    call total%add(array)
    write (array_name,5) "rssprk"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "rssprk"
  end if
  if (allocated(dtsp)) then
    call array%set( storage_size(dtsp,kind=inttype) , &
                            size(dtsp,kind=inttype) )
    call total%add(array)
    write (array_name,5) "dtsp"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "dtsp"
  end if
  if (allocated(troubled_cell)) then
    call array%set( storage_size(troubled_cell,kind=inttype) , &
                            size(troubled_cell,kind=inttype) )
    call total%add(array)
    write (array_name,5) "troubled_cell"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "troubled_cell"
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module flowvar")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  4 format (a,"(",i0,")%",a)
  5 format ("'",a,"'")
  6 format (" Array Component: ",a,dt)
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  8 format (" Storage size for ",a," = ",i0," Bits, ",i0," Bytes")
  !
end subroutine flowvar_memory_usage
!
!###############################################################################
!
end module flowvar
!
!###############################################################################
!
