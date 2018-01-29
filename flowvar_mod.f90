module flowvar
  !
  use module_kind_types, only : wp,lk
  use module_kind_types, only : zero,iout
  use module_kind_types, only : true,fals
  use module_kind_types, only : debug_timer
  use module_kind_types, only : entering_procedure
  use module_kind_types, only : leaving_procedure
  use module_kind_types, only : Classic_RK
  use module_kind_types, only : CK4_RK
  use module_kind_types, only : TVD4_RK
  use module_kind_types, only : alloc_error
  use module_kind_types, only : skip_alloc_pause
  use module_kind_types, only : error_message
  use module_kind_types, only : Geom_Edge
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
  ! RSSPRK  : auxilary array containing old residuals,
  !           needed for some RK methods
  !
  real(wp), public, save, allocatable, dimension(:,:)   :: usp_sum
  real(wp), public, save, allocatable, dimension(:,:,:) :: rssprk
  !
  ! DTSP : Time step for each solution point
  !
  real(wp), public, save, allocatable, dimension(:) :: dtsp
  !
  ! UAVESP : time averaged variables at the solution points
  !
  real(wp), public, save, allocatable :: uavesp(:,:)
  !
  ! TIME_AVE_VARIABLES : character identifier for variables included
  !                      in the time averaging
  !
  character(len=4), public, save, allocatable :: time_ave_variables(:)
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
#ifdef SPECIAL_FOR_PGI
    procedure :: zero_face_var
#else
    procedure :: zero => zero_face_var
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
#ifdef SPECIAL_FOR_PGI
    procedure :: zero_face_vec
#else
    procedure :: zero => zero_face_vec
#endif
  end type face_vec
  !
  ! FACEFLX  : flux normal to cell faces at each interface flux point
  ! FACEUSP  : conservative variables at the interface flux points
  ! FACEDUSP : conservative variable gradients at the interface flux points
  ! FACEXYZ  : coordinates of the flux points on the boundary faces
  !
  ! faceflx(1:nface)%v(1:nq,1:n_sp1de,1:2) ! this is the normal face flux
  ! faceusp(1:nface)%v(1:nq,1:n_sp1de,1:2)
  ! facedusp(1:nface)%v(1:nr,1:nq,1:n_sp1de,1:2)
  !
  type(face_var), public, save, allocatable, dimension(:) :: faceusp
  type(face_vec), public, save, allocatable, dimension(:) :: facedusp
  type(face_var), public, save, allocatable, dimension(:) :: faceflx
  type(face_xyz), public, save, allocatable, dimension(:) :: facexyz
  !
  !
  !
  ! NE1D_T : Derived type containing a 1D array for storing edge/node
  !          solutions for a processor
  !
  type, public :: ne1d_t
    real(wp), allocatable :: v(:)
  contains
#ifdef SPECIAL_FOR_PGI
    procedure :: zero_ne1d_t
#else
    procedure :: zero => zero_ne1d_t
#endif
  end type ne1d_t
  !
  ! NE2D_T : Derived type containing a 2D array for storing edge/node
  !          solutions for a processor
  !
  type, public :: ne2d_t
    real(wp), allocatable :: v(:,:)
  contains
#ifdef SPECIAL_FOR_PGI
    procedure :: zero_ne2d_t
#else
    procedure :: zero => zero_ne2d_t
#endif
  end type ne2d_t
  !
  ! NE3D_T : Derived type containing a 3D array for storing edge/node
  !          solutions for a processor
  !
  type, public :: ne3d_t
    real(wp), allocatable :: v(:,:,:)
  contains
#ifdef SPECIAL_FOR_PGI
    procedure :: zero_ne3d_t
#else
    procedure :: zero => zero_ne3d_t
#endif
  end type ne3d_t
  !
  ! NODEUSP  : conservative variables at the nodes
  !                nodeusp(1:nnode)%v(1:nq)
  !
  type(ne1d_t), public, save, allocatable, dimension(:) :: nodeusp
  !
  ! NODEDUSP : conservative variable gradients at the nodes
  !                nodedusp(1:nnode)%v(1:nr,1:nq)
  !
  type(ne2d_t), public, save, allocatable, dimension(:) :: nodedusp
  !
  ! EDGEUSP  : conservative variables at the edge points
  !                edgeusp(1:nedge)%v(1:nq,1:n_sp1de)
  !
  type(ne2d_t), public, save, allocatable, dimension(:) :: edgeusp
  !
  ! EDGEDUSP : conservative variable gradients at the edge points
  !                edgedusp(1:nedge)%v(1:nr,1:nq,1:n_sp1de)
  !
  type(ne3d_t), public, save, allocatable, dimension(:) :: edgedusp
  !
  ! EDGEXYZ : conservative variable gradients at the edge points
  !               edgexyz(1:nedge)%v(1:nr,1:n_sp1de)
  !
  type(ne2d_t), public, save, allocatable, dimension(:) :: edgexyz
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
  use geovar,    only : nr,n_solpts,n_totpts
  use geovar,    only : nfbnd,nface,face
  use geovar,    only : nnode,nedge,edge
  use order_mod, only : geom_solpts
  use ovar,      only : Runge_Kutta_Scheme,num_rk_stages,mms_opt
  use ovar,      only : convert_restart_to_cgns
  use ovar,      only : continuous_output
  use ovar,      only : time_average_restart_files
  !
  !.. Local Scalars ..
  integer :: n,ne,nf,nn,nepts,nfpts
  integer :: face_geom,face_order,ierr
  character(len=200) :: array_name
  !
  !.. Local Constants ..
  character(len=*), parameter :: pname = "allocate_flowvar_arrays"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  allocate ( usp(1:nq,1:n_totpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"usp",1,__LINE__,__FILE__,ierr,error_message)
  !
  !
  allocate ( dusp(1:nr,1:nq,1:n_totpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"dusp",1,__LINE__,__FILE__,ierr,error_message)
  !
  !
  ! If we are only time averaging restart files, we dont need any of
  ! these other arrays so
  !
  if (time_average_restart_files) return
  !
  allocate ( residual(1:nq,1:n_solpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"residual",1,__LINE__,__FILE__,ierr,error_message)
  !
  !
  allocate ( dtsp(1:n_solpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"dtsp",1,__LINE__,__FILE__,ierr,error_message)
  !
  !
  allocate ( uoldsp(1:nq,1:n_solpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"uoldsp",1,__LINE__,__FILE__,ierr,error_message)
  !
  !
  ! Allocate the outer part of the face solution/flux arrays
  !
  if (.not. convert_restart_to_cgns .or. continuous_output) then
    !
    allocate ( faceusp(1:nface) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"faceusp",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( facedusp(1:nface) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facedusp",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( faceflx(1:nface) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"faceflx",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! If using the continuous output option or MMS, allocate the facexyz array
    !
    if (continuous_output) then
      !
      ! If using continuous output, facexyz needs to be defined
      ! for both interior and boundary faces
      !
      allocate ( facexyz(1:nface) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"facexyz",1,__LINE__,__FILE__,ierr,error_message)
      !
    else if (mms_opt /= 0) then
      !
      ! If using MMS but not using continuous output, facexyz only needs
      ! to be defined for the boundary faces
      !
      allocate ( facexyz(1:nfbnd) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"facexyz",1,__LINE__,__FILE__,ierr,error_message)
      !
    end if
    !
    ! Allocate the data component of these arrays based on the face geometry
    !
    !##########################################################################
    !##########################################################################
    !##  NOTE : THE LOOPS ALLOCATING THE INNER COMPONENTS ARE KEPT SEPARATE  ##
    !##         FROM EACH OTHER TO TRY AND KEEP THE ALLOCATABLE COMPONENTS   ##
    !##         OF EACH ARRAY CONTIGUOUS IN PHYSICAL MEMORY                  ##
    !##########################################################################
    !##########################################################################
    !
    if (allocated(faceusp)) then
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
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
      end do
    end if
    !
    if (allocated(facedusp)) then
      do nf = 1,nface
        !
        face_geom = face(nf)%geom
        face_order = face(nf)%order
        !
        nfpts = geom_solpts(face_geom,face_order)
        !
        allocate ( facedusp(nf)%v(1:nr,1:nq,1:nfpts,1:2) , source=zero , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) "facedusp",nf
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
      end do
    end if
    !
    if (allocated(faceflx)) then
      do nf = 1,nface
        !
        face_geom = face(nf)%geom
        face_order = face(nf)%order
        !
        nfpts = geom_solpts(face_geom,face_order)
        !
        allocate ( faceflx(nf)%v(1:nq,1:nfpts,1:2) , source=zero , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) "faceflx",nf
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
      end do
    end if
    !
    if (allocated(facexyz)) then
      do nf = 1,size(facexyz) ! could be nfbnd or nface
        !
        face_geom = face(nf)%geom
        face_order = face(nf)%order
        !
        nfpts = geom_solpts(face_geom,face_order)
        !
        allocate ( facexyz(nf)%v(1:nr,1:nfpts) , source=zero , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) "facexyz",nf
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
      end do
    end if
    !
    ! If using continuous output, allocate the nodeusp, nodedusp, edgeusp,
    ! and edgedusp arrays
    !
    !##########################################################################
    !##########################################################################
    !##  NOTE : THE LOOPS ALLOCATING THE INNER COMPONENTS ARE KEPT SEPARATE  ##
    !##         FROM EACH OTHER TO TRY AND KEEP THE ALLOCATABLE COMPONENTS   ##
    !##         OF EACH ARRAY CONTIGUOUS IN PHYSICAL MEMORY                  ##
    !##########################################################################
    !##########################################################################
    !
    if (continuous_output) then
      !
      ! Allocate the nodeusp array and its data components
      !
      allocate ( nodeusp(1:nnode) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"nodeusp",1,__LINE__,__FILE__,ierr,error_message)
      !
      do nn = 1,nnode
        allocate ( nodeusp(nn)%v(1:nq) , source=zero , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) "nodeusp",nn
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
      end do
      !
      ! Allocate the nodedusp array and its data components
      !
      allocate ( nodedusp(1:nnode) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"nodedusp",1,__LINE__,__FILE__,ierr,error_message)
      !
      do nn = 1,nnode
        allocate ( nodedusp(nn)%v(1:nr,1:nq) , source=zero , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) "nodedusp",nn
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
      end do
      !
      ! Allocate the edgeusp array and its data components
      !
      allocate ( edgeusp(1:nedge) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"edgeusp",1,__LINE__,__FILE__,ierr,error_message)
      !
      do ne = 1,nedge
        !
        nepts = geom_solpts(Geom_Edge,edge(ne)%order)
        !
        allocate ( edgeusp(ne)%v(1:nq,1:nepts) , source=zero , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) "edgeusp",ne
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
      end do
      !
      ! Allocate the edgedusp array and its data components
      !
      allocate ( edgedusp(1:nedge) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"edgedusp",1,__LINE__,__FILE__,ierr,error_message)
      !
      do ne = 1,nedge
        !
        nepts = geom_solpts(Geom_Edge,edge(ne)%order)
        !
        allocate ( edgedusp(ne)%v(1:nr,1:nq,1:nepts) , source=zero , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) "edgedusp",ne
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
      end do
      !
      ! Allocate the edgexyz array and its data components
      !
      allocate ( edgexyz(1:nedge) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"edgexyz",1,__LINE__,__FILE__,ierr,error_message)
      !
      do ne = 1,nedge
        !
        nepts = geom_solpts(Geom_Edge,edge(ne)%order)
        !
        allocate ( edgexyz(ne)%v(1:nr,1:nepts) , source=zero , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) "edgexyz",ne
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message,skip_alloc_pause)
        !
      end do
      !
    end if
    !
  end if
  !
  ! Allocate the storage for the Runge-Kutta methods
  !
  n = 0
  !
  if (Runge_Kutta_Scheme == Classic_RK) then
    n = num_rk_stages
  else if ( any(Runge_Kutta_Scheme == [CK4_RK,TVD4_RK]) ) then
    n = 1
    if (Runge_Kutta_Scheme == TVD4_RK) then
      allocate ( usp_sum(1:nq,1:n_solpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"usp_sum",1,__LINE__,__FILE__,ierr,error_message)
    end if
  end if
  !
  if (n /= 0) then
    allocate ( rssprk(1:nq,1:n_solpts,1:n) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"rssprk",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (a,"(",i0,")%v")
  !
end subroutine allocate_flowvar_arrays
!
!###############################################################################
!###############################################################################
!
!   FACE RELATED PROCEDURES FOR THE FACE_VAR DERIVED TYPE
!
!###############################################################################
!###############################################################################
!
elemental subroutine zero_face_var(this)
  class(face_var), intent(inout) :: this
  if (allocated(this%v)) this%v = zero
end subroutine zero_face_var
!
!###############################################################################
!###############################################################################
!
!   FACE RELATED PROCEDURES FOR THE FACE_VEC DERIVED TYPE
!
!###############################################################################
!###############################################################################
!
elemental subroutine zero_face_vec(this)
  class(face_vec), intent(inout) :: this
  if (allocated(this%v)) this%v = zero
end subroutine zero_face_vec
!
!###############################################################################
!###############################################################################
!
!   EDGE/NODE RELATED PROCEDURES FOR THE NE1D_T DERIVED TYPE
!
!###############################################################################
!###############################################################################
!
elemental subroutine zero_ne1d_t(this)
  class(ne1d_t), intent(inout) :: this
  if (allocated(this%v)) this%v = zero
end subroutine zero_ne1d_t
!
!###############################################################################
!###############################################################################
!
!   EDGE/NODE RELATED PROCEDURES FOR THE NE2D_T DERIVED TYPE
!
!###############################################################################
!###############################################################################
!
elemental subroutine zero_ne2d_t(this)
  class(ne2d_t), intent(inout) :: this
  if (allocated(this%v)) this%v = zero
end subroutine zero_ne2d_t
!
!###############################################################################
!###############################################################################
!
!   EDGE/NODE RELATED PROCEDURES FOR THE NE3D_T DERIVED TYPE
!
!###############################################################################
!###############################################################################
!
elemental subroutine zero_ne3d_t(this)
  class(ne3d_t), intent(inout) :: this
  if (allocated(this%v)) this%v = zero
end subroutine zero_ne3d_t
!
!###############################################################################
!###############################################################################
!
!   MEMORY USAGE
!
!###############################################################################
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
#ifndef DISABLE_DTIO
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
    call array%set( storage_size(faceusp,kind=inttype) )
    !
    write (array_name,5) "type(face_var)"
    write (iou,6) array_name,array%get_units()
    !
    call array%set( storage_size(faceusp,kind=inttype) , &
                            size(faceusp,kind=inttype) )
    !
    write (array_name,5) "faceusp(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
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
    call array%set( storage_size(facedusp,kind=inttype) )
    !
    write (array_name,5) "type(face_vec)"
    write (iou,6) array_name,array%get_units()
    !
    call array%set( storage_size(facedusp,kind=inttype) , &
                            size(facedusp,kind=inttype) )
    !
    write (array_name,5) "facedusp(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
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
    call array%set( storage_size(faceflx,kind=inttype) )
    !
    write (array_name,5) "type(face_var)"
    write (iou,6) array_name,array%get_units()
    !
    call array%set( storage_size(faceflx,kind=inttype) , &
                            size(faceflx,kind=inttype) )
    !
    write (array_name,5) "faceflx(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
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
    call array%set( storage_size(facexyz,kind=inttype) )
    !
    write (array_name,5) "type(face_xyz)"
    write (iou,6) array_name,array%get_units()
    !
    call array%set( storage_size(facexyz,kind=inttype) , &
                            size(facexyz,kind=inttype) )
    !
    write (array_name,5) "facexyz(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
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
  if (allocated(edgeusp)) then
    !
    call array%set( storage_size(edgeusp,kind=inttype) )
    !
    write (array_name,5) "type(ne2d_t)"
    write (iou,6) array_name,array%get_units()
    !
    call array%set( storage_size(edgeusp,kind=inttype) , &
                            size(edgeusp,kind=inttype) )
    !
    write (array_name,5) "edgeusp(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    call array%reset
    !
    do i = lbound(edgeusp,dim=1),ubound(edgeusp,dim=1)
      !
      call array%add( storage_size(edgeusp(i)%v,kind=inttype) , &
                              size(edgeusp(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "edgeusp(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "edgeusp"
    write (iou,2) array_name,dt%get_units()
    !
  else
    write (iou,7) "edgeusp"
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(edgedusp)) then
    !
    call array%set( storage_size(edgedusp,kind=inttype) )
    !
    write (array_name,5) "type(ne3d_t)"
    write (iou,6) array_name,array%get_units()
    !
    call array%set( storage_size(edgedusp,kind=inttype) , &
                            size(edgedusp,kind=inttype) )
    !
    write (array_name,5) "edgedusp(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    call array%reset
    !
    do i = lbound(edgedusp,dim=1),ubound(edgedusp,dim=1)
      !
      call array%add( storage_size(edgedusp(i)%v,kind=inttype) , &
                              size(edgedusp(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "edgedusp(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "edgedusp"
    write (iou,2) array_name,dt%get_units()
    !
  else
    write (iou,7) "edgedusp"
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(edgexyz)) then
    !
    call array%set( storage_size(edgexyz,kind=inttype) )
    !
    write (array_name,5) "type(ne2d_t)"
    write (iou,6) array_name,array%get_units()
    !
    call array%set( storage_size(edgexyz,kind=inttype) , &
                            size(edgexyz,kind=inttype) )
    !
    write (array_name,5) "edgexyz(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    call array%reset
    !
    do i = lbound(edgexyz,dim=1),ubound(edgexyz,dim=1)
      !
      call array%add( storage_size(edgexyz(i)%v,kind=inttype) , &
                              size(edgexyz(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "edgexyz(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "edgexyz"
    write (iou,2) array_name,dt%get_units()
    !
  else
    write (iou,7) "edgexyz"
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(nodeusp)) then
    !
    call array%set( storage_size(nodeusp,kind=inttype) )
    !
    write (array_name,5) "type(ne1d_t)"
    write (iou,6) array_name,array%get_units()
    !
    call array%set( storage_size(nodeusp,kind=inttype) , &
                            size(nodeusp,kind=inttype) )
    !
    write (array_name,5) "nodeusp(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    call array%reset
    !
    do i = lbound(nodeusp,dim=1),ubound(nodeusp,dim=1)
      !
      call array%add( storage_size(nodeusp(i)%v,kind=inttype) , &
                              size(nodeusp(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "nodeusp(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "nodeusp"
    write (iou,2) array_name,dt%get_units()
    !
  else
    write (iou,7) "nodeusp"
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(nodedusp)) then
    !
    call array%set( storage_size(nodedusp,kind=inttype) )
    !
    write (array_name,5) "type(ne2d_t)"
    write (iou,6) array_name,array%get_units()
    !
    call array%set( storage_size(nodedusp,kind=inttype) , &
                            size(nodedusp,kind=inttype) )
    !
    write (array_name,5) "nodedusp(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    call array%reset
    !
    do i = lbound(nodedusp,dim=1),ubound(nodedusp,dim=1)
      !
      call array%add( storage_size(nodedusp(i)%v,kind=inttype) , &
                              size(nodedusp(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "nodedusp(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "nodedusp"
    write (iou,2) array_name,dt%get_units()
    !
  else
    write (iou,7) "nodedusp"
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
#endif
end subroutine flowvar_memory_usage
!
!###############################################################################
!
end module flowvar
!
!###############################################################################
!
