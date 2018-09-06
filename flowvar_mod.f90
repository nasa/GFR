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
  ! SOL_NODES : conservation variables at plot3d nodes
  !
  real(wp), public, save, allocatable :: sol_nodes(:,:)
  !
  ! USP  : conservative variables at the solution points (@SP)
  ! DUSP : gradient of conservative variables @SP
  ! GPSP : gas properties at the solution points
  !
  real(wp), public, save, allocatable, dimension(:,:)   :: usp
  !$ACC DECLARE CREATE(usp)
  real(wp), public, save, allocatable, dimension(:,:,:) :: dusp
  !$ACC DECLARE CREATE(dusp)
  real(wp), public, save, allocatable, dimension(:,:) :: gpsp
  !$ACC DECLARE CREATE(gpsp)
  !
  ! UOLDSP   : conservative variables of previous time step @SP
  ! RESIDUAL : RHS residual @SP for current time/RK step
  !
  real(wp), public, save, allocatable, dimension(:,:) :: uoldsp
  !$ACC DECLARE CREATE(uoldsp)
  real(wp), public, save, allocatable, dimension(:,:) :: residual
  !$ACC DECLARE CREATE(residual)
  !
  ! GLOBAL_USP : conservative variables @SP across all processors
  ! GLOBAL_VAR : arbitrary    variable  @SP across all processors
  !
  real(wp), public, save, allocatable, dimension(:,:) :: global_usp
  real(wp), public, save, allocatable, dimension(:)   :: global_var
  !
  ! RK_SUM : auxilary array needed for SSP-RK4 method
  ! RSSPRK : auxilary array containing old residuals,
  !          needed for some RK methods
  !
  real(wp), public, save, allocatable, dimension(:,:)   :: rk_sum
  !$ACC DECLARE CREATE(rk_sum)
  real(wp), public, save, allocatable, dimension(:,:,:) :: rssprk
  !$ACC DECLARE CREATE(rssprk)
  !
  ! DTSP : Time step for each solution point
  !
  real(wp), public, save, allocatable, dimension(:) :: dtsp
  !$ACC DECLARE CREATE(dtsp)
  !
  ! UAVESP : time averaged variables at the solution points
  !
  real(wp), public, save, allocatable :: uavesp(:,:)
  !$ACC DECLARE CREATE(uavesp)
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
  ! FACEUSP_DT  : conservative variables at the interface flux points
  ! FACEGPSP_DT : gas properties at the interface flux points
  ! FACEDUSP_DT : conservative variable gradients at the interface flux points
  ! FACEXYZ_DT  : coordinates of the flux points on the boundary faces
  !
  ! faceusp_dt(1:nface)%v(1:nq,1:n_sp1de,1:2)
  ! facedusp_dt(1:nface)%v(1:nr,1:nq,1:n_sp1de,1:2)
  ! facexyz_dt(1:nface)%v(1:nr,1:n_sp1de)
  !
  type(face_var), public, save, allocatable, dimension(:) :: faceusp_dt
  type(face_var), public, save, allocatable, dimension(:) :: facegpsp_dt
  type(face_vec), public, save, allocatable, dimension(:) :: facedusp_dt
  type(face_xyz), public, save, allocatable, dimension(:) :: facexyz_dt
  !
  ! FACEUSP_ARRAY  : conservative variables at the interface flx-pts
  ! FACEDUSP_ARRAY : conservative variable gradients at the interface flx-pts
  !
  real(wp), public, save, allocatable :: faceusp_array(:,:,:)
  !$ACC DECLARE CREATE(faceusp_array)
  real(wp), public, save, allocatable :: facedusp_array(:,:,:,:)
  !$ACC DECLARE CREATE(facedusp_array)
  !
  ! FACEFLX  : flux normal to cell faces at each interface flux point
  !
  real(wp), public, save, allocatable :: faceflx(:,:)
  !$ACC DECLARE CREATE(faceflx)
  !
  ! FACEXYZ_ARRAY  : coordinates of the interface flux points
  !
  real(wp), public, save, allocatable :: facexyz_array(:,:)
  !$ACC DECLARE CREATE(facexyz_array)
  !
  ! FACEGPSP : gas properties at the interface flux points
  !
  real(wp), public, save, allocatable :: facegpsp(:,:,:)
  !$ACC DECLARE CREATE(facegpsp)
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
  ! NODEUSP_DT  : conservative variables at the nodes
  !                   nodeusp_dt(1:nnode)%v(1:nq)
  !
  type(ne1d_t), public, save, allocatable, dimension(:) :: nodeusp_dt
  !
  ! NODEGPSP_DT : gas properties at the nodes
  !                   nodegpsp_dt(1:nnode)%v(1:ng)
  !
  type(ne1d_t), public, save, allocatable, dimension(:) :: nodegpsp_dt
  !
  ! NODEDUSP_DT : conservative variable gradients at the nodes
  !                   nodedusp_dt(1:nnode)%v(1:nr,1:nq)
  !
  type(ne2d_t), public, save, allocatable, dimension(:) :: nodedusp_dt
  !
  ! EDGEUSP_DT  : conservative variables at the edge points
  !                   edgeusp_dt(1:nedge)%v(1:nq,1:n_sp1de)
  !
  type(ne2d_t), public, save, allocatable, dimension(:) :: edgeusp_dt
  !
  ! EDGEGPSP_DT : gas properties at the edge points
  !                   edgegpsp_dt(1:nedge)%v(1:ng,1:n_sp1de)
  !
  type(ne2d_t), public, save, allocatable, dimension(:) :: edgegpsp_dt
  !
  ! EDGEDUSP_DT : conservative variable gradients at the edge points
  !                   edgedusp_dt(1:nedge)%v(1:nr,1:nq,1:n_sp1de)
  !
  type(ne3d_t), public, save, allocatable, dimension(:) :: edgedusp_dt
  !
  ! EDGEXYZ_DT  : conservative variable gradients at the edge points
  !                   edgexyz_dt(1:nedge)%v(1:nr,1:n_sp1de)
  !
  type(ne2d_t), public, save, allocatable, dimension(:) :: edgexyz_dt
  !
  !
  !
  ! ###################################
  ! ###  MODULE GENERIC INTERFACES  ###
  ! ###################################
  !
  interface alloc_face_dt_components
    module procedure alloc_face_xyz_dt_components, &
                     alloc_face_var_dt_components, &
                     alloc_face_vec_dt_components
  end interface alloc_face_dt_components
  !
  interface alloc_ne_dt_components
    module procedure alloc_node_var_dt_components, &
                     alloc_node_vec_dt_components, &
                     alloc_edge_var_dt_components, &
                     alloc_edge_vec_dt_components
  end interface alloc_ne_dt_components
  !
contains
!
!###############################################################################
!
subroutine allocate_flowvar_arrays()
  !
  !.. Use Statements ..
  use eqn_idx,   only : nq,ng
  use geovar,    only : nr,n_solpts,n_totpts
  use geovar,    only : nfbnd,nface,face
  use geovar,    only : nnode,nedge,edge
  use order_mod, only : geom_solpts
  use ovar,      only : Runge_Kutta_Scheme,num_rk_stages,mms_opt
  use ovar,      only : convert_restart_to_cgns
  use ovar,      only : continuous_output
  use ovar,      only : time_average_restart_files
  use ovar,      only : if_interpto_newgrid
  !
  !.. Local Scalars ..
  integer :: n,ne,nf,nn,nepts,nfp,nfpts
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
  if (.not. allocated(usp)) then
    allocate ( usp(1:n_solpts,1:nq) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"usp",1,__LINE__,__FILE__,ierr,error_message)
    !$ACC UPDATE DEVICE(usp)
  end if
  !
  !
  allocate ( dusp(1:n_solpts,1:nr,1:nq) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"dusp",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(dusp)
  !
  allocate ( gpsp(1:n_solpts,1:ng) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"gpsp",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(gpsp)
  !
  ! If we are only time averaging restart files, we dont need any of
  ! these other arrays so return out of here
  !
  if (time_average_restart_files) return
  !
  !
  ! If we are only interpolating the solution to a new grid,
  ! then only need facexyz_dt before returning out of here.
  !
  if ( if_interpto_newgrid() ) then
    !
    ! Deallocate dusp because we do not need it when interpolating the
    ! solution to a new grid and it takes up too much memory.
    !
    deallocate ( dusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"dusp",2,__LINE__,__FILE__,ierr,error_message)
    !
    ! The interpolation utility needs facexyz_array to
    ! be defined for both interior and boundary faces.
    !
    allocate ( facexyz_dt(1:nface) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facexyz_dt",1,__LINE__,__FILE__,ierr,error_message)
    !
    call alloc_face_dt_components(facexyz_dt,nfpts,"facexyz_dt")
    !
    allocate ( facexyz_array(1:nr,1:nfpts) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facexyz_array",1,__LINE__,__FILE__, &
                     ierr,error_message)
    !
    ! Procedure interpsol needs faceusp_dt to print debug info
    !
    allocate ( faceusp_dt(1:nface) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"faceusp_dt",1,__LINE__,__FILE__,ierr,error_message)
    !
    call alloc_face_dt_components(faceusp_dt,nfpts,"faceusp_dt")
    !
    allocate ( faceusp_array(1:nfpts,1:nq,1:2) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"faceusp_array",1,__LINE__,__FILE__, &
                     ierr,error_message)
    !
    return
    !
  end if
  !
  !
  allocate ( residual(1:n_solpts,1:nq) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"residual",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(residual)
  !
  !
  allocate ( dtsp(1:n_solpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"dtsp",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(dtsp)
  !
  !
  allocate ( uoldsp(1:n_solpts,1:nq) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"uoldsp",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(uoldsp)
  !
  !
  if (convert_restart_to_cgns) then
    !
    nfpts = face(nface)%end_fp ! end fp index for last face
    !
    allocate ( faceusp_array(1:nfpts,1:nq,1:2) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"faceusp_array",1,__LINE__,__FILE__, &
                     ierr,error_message)
    !$ACC UPDATE DEVICE(faceusp_array)
   !!
   !allocate ( facedusp_array(1:nfpts,1:nr,1:nq,1:2) , source=zero , &
   !           stat=ierr , errmsg=error_message )
   !call alloc_error(pname,"facedusp_array",1,__LINE__,__FILE__, &
   !                 ierr,error_message)
   !!$ACC UPDATE DEVICE(facedusp_array)
    !
  end if
  !
  ! Allocate the outer part of the face solution/flux arrays
  !
  if (.not. convert_restart_to_cgns .or. continuous_output) then
    !
    ! If using the continuous output option or MMS,
    ! allocate the facexyz_dt array
    !
    if (continuous_output) then
      !
      ! If using continuous output, facexyz_dt needs to be defined
      ! for both interior and boundary faces
      !
      allocate ( facexyz_dt(1:nface) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"facexyz_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      if (.not. allocated(faceusp_dt)) then
        allocate ( faceusp_dt(1:nface) , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"faceusp_dt",1,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
      !
      if (.not. allocated(facegpsp_dt)) then
        allocate ( facegpsp_dt(1:nface) , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"facegpsp_dt",1,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
      !
      if (.not. allocated(facedusp_dt)) then
        allocate ( facedusp_dt(1:nface) , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"facedusp_dt",1,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
      !
    else if (mms_opt /= 0) then
      !
      ! If using MMS but not using continuous output, facexyz_dt only needs
      ! to be defined for the boundary faces
      !
      allocate ( facexyz_dt(1:nfbnd) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"facexyz_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
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
    if (allocated(faceusp_dt)) then
      call alloc_face_dt_components(faceusp_dt,nfpts,"faceusp_dt")
    end if
    !
    if (allocated(facegpsp_dt)) then
      call alloc_face_dt_components(facegpsp_dt,nfpts,"facegpsp_dt",num_var=ng)
    end if
    !
    if (allocated(facedusp_dt)) then
      call alloc_face_dt_components(facedusp_dt,nfpts,"facedusp_dt")
    end if
    !
    !
    nfpts = face(nface)%end_fp ! end fp index for last face
    !
    allocate ( faceflx(1:nfpts,1:nq) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"faceflx",1,__LINE__,__FILE__, &
                     ierr,error_message)
    !$ACC UPDATE DEVICE(faceflx)
    !
    if (.not. allocated(faceusp_array)) then
      allocate ( faceusp_array(1:nfpts,1:nq,1:2) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"faceusp_array",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !$ACC UPDATE DEVICE(faceusp_array)
    end if
    !
    if (.not. allocated(facedusp_array)) then
      allocate ( facedusp_array(1:nfpts,1:nr,1:nq,1:2) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"facedusp_array",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !$ACC UPDATE DEVICE(facedusp_array)
    end if
    !
    if (.not. allocated(facegpsp)) then
      allocate ( facegpsp(1:nfpts,1:ng,1:2) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"facegpsp",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !$ACC UPDATE DEVICE(facegpsp)
    end if
    !
    !
    if (allocated(facexyz_dt)) then
      !
      call alloc_face_dt_components(facexyz_dt,nfpts,"facexyz_dt")
      !
      allocate ( facexyz_array(1:nr,1:nfpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"facexyz_array",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      !$ACC UPDATE DEVICE(facexyz_array)
      !
    end if
    !
    ! If using continuous output, allocate the facegpsp_dt, nodeusp_dt,
    ! nodegpsp_dt, nodedusp_dt, edgeusp_dt, edgegpsp_dt, and edgedusp_dt arrays
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
      ! Allocate the nodeusp_dt array and its data components
      !
      allocate ( nodeusp_dt(1:nnode) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"nodeusp_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      call alloc_ne_dt_components(nodeusp_dt,"nodeusp_dt")
      !
      ! Allocate the nodegpsp_dt array and its data components
      !
      allocate ( nodegpsp_dt(1:nnode) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"nodegpsp_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      call alloc_ne_dt_components(nodegpsp_dt,"nodegpsp_dt",num_var=ng)
      !
      ! Allocate the nodedusp_dt array and its data components
      !
      allocate ( nodedusp_dt(1:nnode) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"nodedusp_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      call alloc_ne_dt_components(nodedusp_dt,"nodedusp_dt")
      !
      ! Allocate the edgeusp_dt array and its data components
      !
      allocate ( edgeusp_dt(1:nedge) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"edgeusp_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      call alloc_ne_dt_components(edgeusp_dt,nepts,"edgeusp_dt")
      !
      ! Allocate the edgegpsp_dt array and its data components
      !
      allocate ( edgegpsp_dt(1:nedge) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"edgegpsp_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      call alloc_ne_dt_components(edgegpsp_dt,nepts,"edgegpsp_dt",num_var=ng)
      !
      ! Allocate the edgedusp_dt array and its data components
      !
      allocate ( edgedusp_dt(1:nedge) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"edgedusp_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      call alloc_ne_dt_components(edgedusp_dt,nepts,"edgedusp_dt")
      !
      ! Allocate the edgexyz_dt array and its data components
      !
      allocate ( edgexyz_dt(1:nedge) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"edgexyz_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      call alloc_ne_dt_components(edgexyz_dt,nepts,"edgexyz_dt",num_var=nr)
      !
    end if
    !
  end if
  !
  ! Allocate the storage for the Runge-Kutta methods
  !
  if (Runge_Kutta_Scheme == Classic_RK) then
    allocate ( rssprk(1:n_solpts,1:nq,1:num_rk_stages) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"rssprk",1,__LINE__,__FILE__,ierr,error_message)
    !$ACC UPDATE DEVICE(rssprk)
  end if
  !
  if ( any(Runge_Kutta_Scheme == [CK4_RK,TVD4_RK]) ) then
    allocate ( rk_sum(1:n_solpts,1:nq) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"rk_sum",1,__LINE__,__FILE__,ierr,error_message)
    !$ACC UPDATE DEVICE(rk_sum)
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
!
subroutine alloc_face_xyz_dt_components(array,npts,aname,num_var)
  !
  !.. Use Statements ..
  use geovar,    only : nr,face
  use order_mod, only : geom_solpts
  !
  !.. Formal Arguments ..
  type(face_xyz), intent(inout) :: array(:)
  integer,        intent(inout) :: npts
  character(len=*),  intent(in) :: aname
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: num_var
  !
  !.. Local Scalars ..
  integer :: nf,nfp,nvar,ierr
  integer :: face_geom,face_order
  !
  character(len=200) :: array_name
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "alloc_face_xyz_dt_components"
  !
continue
  !
  nvar = nr
  if (present(num_var)) then
    nvar = num_var
  end if
  !
  npts = 0
  !
  do nf = lbound(array,dim=1),ubound(array,dim=1)
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    nfp = geom_solpts(face_geom,face_order)
    !
    npts = npts + nfp
    !
    allocate ( array(nf)%v(1:nvar,1:nfp) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) trim(adjustl(aname)),nf
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !
  end do
  !
  ! Format Statements
  !
  1 format (a,"(",i0,")%v")
  !
end subroutine alloc_face_xyz_dt_components
!
!###############################################################################
!
subroutine alloc_face_var_dt_components(array,npts,aname,num_var)
  !
  !.. Use Statements ..
  use eqn_idx,   only : nq
  use geovar,    only : face
  use order_mod, only : geom_solpts
  !
  !.. Formal Arguments ..
  type(face_var), intent(inout) :: array(:)
  integer,        intent(inout) :: npts
  character(len=*),  intent(in) :: aname
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: num_var
  !
  !.. Local Scalars ..
  integer :: nf,nfp,nvar,ierr
  integer :: face_geom,face_order
  !
  character(len=200) :: array_name
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "alloc_face_var_dt_components"
  !
continue
  !
  nvar = nq
  if (present(num_var)) then
    nvar = num_var
  end if
  !
  npts = 0
  !
  do nf = lbound(array,dim=1),ubound(array,dim=1)
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    nfp = geom_solpts(face_geom,face_order)
    !
    npts = npts + nfp
    !
    allocate ( array(nf)%v(1:nvar,1:nfp,1:2) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) trim(adjustl(aname)),nf
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !
  end do
  !
  ! Format Statements
  !
  1 format (a,"(",i0,")%v")
  !
end subroutine alloc_face_var_dt_components
!
!###############################################################################
!
subroutine alloc_face_vec_dt_components(array,npts,aname,num_var)
  !
  !.. Use Statements ..
  use eqn_idx,   only : nq
  use geovar,    only : nr,face
  use order_mod, only : geom_solpts
  !
  !.. Formal Arguments ..
  type(face_vec), intent(inout) :: array(:)
  integer,        intent(inout) :: npts
  character(len=*),  intent(in) :: aname
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: num_var
  !
  !.. Local Scalars ..
  integer :: nf,nfp,nvar,ierr
  integer :: face_geom,face_order
  !
  character(len=200) :: array_name
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "alloc_face_vec_dt_components"
  !
continue
  !
  nvar = nq
  if (present(num_var)) then
    nvar = num_var
  end if
  !
  npts = 0
  !
  do nf = lbound(array,dim=1),ubound(array,dim=1)
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    nfp = geom_solpts(face_geom,face_order)
    !
    npts = npts + nfp
    !
    allocate ( array(nf)%v(1:nr,1:nvar,1:nfp,1:2) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) trim(adjustl(aname)),nf
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !
  end do
  !
  ! Format Statements
  !
  1 format (a,"(",i0,")%v")
  !
end subroutine alloc_face_vec_dt_components
!
!###############################################################################
!
subroutine alloc_node_var_dt_components(array,aname,num_var)
  !
  !.. Use Statements ..
  use eqn_idx, only : nq
  !
  !.. Formal Arguments ..
  type(ne1d_t),  intent(inout) :: array(:)
  character(len=*), intent(in) :: aname
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: num_var
  !
  !.. Local Scalars ..
  integer :: nn,nvar,ierr
  !
  character(len=200) :: array_name
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "alloc_node_var_dt_components"
  !
continue
  !
  nvar = nq
  if (present(num_var)) then
    nvar = num_var
  end if
  !
  do nn = lbound(array,dim=1),ubound(array,dim=1)
    !
    allocate ( array(nn)%v(1:nvar) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) trim(adjustl(aname)),nn
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !
  end do
  !
  ! Format Statements
  !
  1 format (a,"(",i0,")%v")
  !
end subroutine alloc_node_var_dt_components
!
!###############################################################################
!
subroutine alloc_node_vec_dt_components(array,aname,num_var)
  !
  !.. Use Statements ..
  use eqn_idx, only : nq
  use geovar,  only : nr
  !
  !.. Formal Arguments ..
  type(ne2d_t),  intent(inout) :: array(:)
  character(len=*), intent(in) :: aname
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: num_var
  !
  !.. Local Scalars ..
  integer :: nn,nvar,ierr
  !
  character(len=200) :: array_name
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "alloc_node_vec_dt_components"
  !
continue
  !
  nvar = nq
  if (present(num_var)) then
    nvar = num_var
  end if
  !
  do nn = lbound(array,dim=1),ubound(array,dim=1)
    !
    allocate ( array(nn)%v(1:nr,1:nvar) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) trim(adjustl(aname)),nn
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !
  end do
  !
  ! Format Statements
  !
  1 format (a,"(",i0,")%v")
  !
end subroutine alloc_node_vec_dt_components
!
!###############################################################################
!
subroutine alloc_edge_var_dt_components(array,npts,aname,num_var)
  !
  !.. Use Statements ..
  use eqn_idx,   only : nq
  use geovar,    only : edge
  use order_mod, only : geom_solpts
  !
  !.. Formal Arguments ..
  type(ne2d_t),  intent(inout) :: array(:)
  integer,       intent(inout) :: npts
  character(len=*), intent(in) :: aname
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: num_var
  !
  !.. Local Scalars ..
  integer :: ne,nepts,nvar,ierr
  !
  character(len=200) :: array_name
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "alloc_edge_var_dt_components"
  !
continue
  !
  nvar = nq
  if (present(num_var)) then
    nvar = num_var
  end if
  !
  npts = 0
  !
  do ne = lbound(array,dim=1),ubound(array,dim=1)
    !
    nepts = geom_solpts(Geom_Edge,edge(ne)%order)
    !
    npts = npts + nepts
    !
    allocate ( array(ne)%v(1:nvar,1:nepts) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) trim(adjustl(aname)),ne
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !
  end do
  !
  ! Format Statements
  !
  1 format (a,"(",i0,")%v")
  !
end subroutine alloc_edge_var_dt_components
!
!###############################################################################
!
subroutine alloc_edge_vec_dt_components(array,npts,aname,num_var)
  !
  !.. Use Statements ..
  use eqn_idx,   only : nq
  use geovar,    only : nr
  use geovar,    only : edge
  use order_mod, only : geom_solpts
  !
  !.. Formal Arguments ..
  type(ne3d_t),  intent(inout) :: array(:)
  integer,       intent(inout) :: npts
  character(len=*), intent(in) :: aname
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: num_var
  !
  !.. Local Scalars ..
  integer :: ne,nepts,nvar,ierr
  !
  character(len=200) :: array_name
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "alloc_edge_vec_dt_components"
  !
continue
  !
  nvar = nq
  if (present(num_var)) then
    nvar = num_var
  end if
  !
  npts = 0
  !
  do ne = lbound(array,dim=1),ubound(array,dim=1)
    !
    nepts = geom_solpts(Geom_Edge,edge(ne)%order)
    !
    npts = npts + nepts
    !
    allocate ( array(ne)%v(1:nr,1:nvar,1:nepts) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) trim(adjustl(aname)),ne
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !
  end do
  !
  ! Format Statements
  !
  1 format (a,"(",i0,")%v")
  !
end subroutine alloc_edge_vec_dt_components
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
  if (allocated(faceusp_dt)) then
    !
    call array%set( storage_size(faceusp_dt,kind=inttype) )
    !
    write (array_name,5) "type(face_var)"
    write (iou,6) array_name,array%get_units()
    !
    call array%set( storage_size(faceusp_dt,kind=inttype) , &
                            size(faceusp_dt,kind=inttype) )
    !
    write (array_name,5) "faceusp_dt(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    call array%reset
    !
    do i = lbound(faceusp_dt,dim=1),ubound(faceusp_dt,dim=1)
      !
      call array%add( storage_size(faceusp_dt(i)%v,kind=inttype) , &
                              size(faceusp_dt(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "faceusp_dt(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "faceusp_dt"
    write (iou,2) array_name,dt%get_units()
    !
  else
    write (iou,7) "faceusp_dt"
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(facedusp_dt)) then
    !
    call array%set( storage_size(facedusp_dt,kind=inttype) )
    !
    write (array_name,5) "type(face_vec)"
    write (iou,6) array_name,array%get_units()
    !
    call array%set( storage_size(facedusp_dt,kind=inttype) , &
                            size(facedusp_dt,kind=inttype) )
    !
    write (array_name,5) "facedusp_dt(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    call array%reset
    !
    do i = lbound(facedusp_dt,dim=1),ubound(facedusp_dt,dim=1)
      !
      call array%add( storage_size(facedusp_dt(i)%v,kind=inttype) , &
                              size(facedusp_dt(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "facedusp_dt(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "facedusp_dt"
    write (iou,2) array_name,dt%get_units()
    !
  else
    write (iou,7) "facedusp_dt"
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(facexyz_dt)) then
    !
    call array%set( storage_size(facexyz_dt,kind=inttype) )
    !
    write (array_name,5) "type(face_xyz)"
    write (iou,6) array_name,array%get_units()
    !
    call array%set( storage_size(facexyz_dt,kind=inttype) , &
                            size(facexyz_dt,kind=inttype) )
    !
    write (array_name,5) "facexyz_dt(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    call array%reset
    !
    do i = lbound(facexyz_dt,dim=1),ubound(facexyz_dt,dim=1)
      !
      call array%add( storage_size(facexyz_dt(i)%v,kind=inttype) , &
                              size(facexyz_dt(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "facexyz_dt(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "facexyz_dt"
    write (iou,2) array_name,dt%get_units()
    !
  else
    write (iou,7) "facexyz_dt"
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(edgeusp_dt)) then
    !
    call array%set( storage_size(edgeusp_dt,kind=inttype) )
    !
    write (array_name,5) "type(ne2d_t)"
    write (iou,6) array_name,array%get_units()
    !
    call array%set( storage_size(edgeusp_dt,kind=inttype) , &
                            size(edgeusp_dt,kind=inttype) )
    !
    write (array_name,5) "edgeusp_dt(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    call array%reset
    !
    do i = lbound(edgeusp_dt,dim=1),ubound(edgeusp_dt,dim=1)
      !
      call array%add( storage_size(edgeusp_dt(i)%v,kind=inttype) , &
                              size(edgeusp_dt(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "edgeusp_dt(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "edgeusp_dt"
    write (iou,2) array_name,dt%get_units()
    !
  else
    write (iou,7) "edgeusp_dt"
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(edgedusp_dt)) then
    !
    call array%set( storage_size(edgedusp_dt,kind=inttype) )
    !
    write (array_name,5) "type(ne3d_t)"
    write (iou,6) array_name,array%get_units()
    !
    call array%set( storage_size(edgedusp_dt,kind=inttype) , &
                            size(edgedusp_dt,kind=inttype) )
    !
    write (array_name,5) "edgedusp_dt(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    call array%reset
    !
    do i = lbound(edgedusp_dt,dim=1),ubound(edgedusp_dt,dim=1)
      !
      call array%add( storage_size(edgedusp_dt(i)%v,kind=inttype) , &
                              size(edgedusp_dt(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "edgedusp_dt(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "edgedusp_dt"
    write (iou,2) array_name,dt%get_units()
    !
  else
    write (iou,7) "edgedusp_dt"
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(edgexyz_dt)) then
    !
    call array%set( storage_size(edgexyz_dt,kind=inttype) )
    !
    write (array_name,5) "type(ne2d_t)"
    write (iou,6) array_name,array%get_units()
    !
    call array%set( storage_size(edgexyz_dt,kind=inttype) , &
                            size(edgexyz_dt,kind=inttype) )
    !
    write (array_name,5) "edgexyz_dt(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    call array%reset
    !
    do i = lbound(edgexyz_dt,dim=1),ubound(edgexyz_dt,dim=1)
      !
      call array%add( storage_size(edgexyz_dt(i)%v,kind=inttype) , &
                              size(edgexyz_dt(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "edgexyz_dt(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "edgexyz_dt"
    write (iou,2) array_name,dt%get_units()
    !
  else
    write (iou,7) "edgexyz_dt"
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(nodeusp_dt)) then
    !
    call array%set( storage_size(nodeusp_dt,kind=inttype) )
    !
    write (array_name,5) "type(ne1d_t)"
    write (iou,6) array_name,array%get_units()
    !
    call array%set( storage_size(nodeusp_dt,kind=inttype) , &
                            size(nodeusp_dt,kind=inttype) )
    !
    write (array_name,5) "nodeusp_dt(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    call array%reset
    !
    do i = lbound(nodeusp_dt,dim=1),ubound(nodeusp_dt,dim=1)
      !
      call array%add( storage_size(nodeusp_dt(i)%v,kind=inttype) , &
                              size(nodeusp_dt(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "nodeusp_dt(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "nodeusp_dt"
    write (iou,2) array_name,dt%get_units()
    !
  else
    write (iou,7) "nodeusp_dt"
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(nodedusp_dt)) then
    !
    call array%set( storage_size(nodedusp_dt,kind=inttype) )
    !
    write (array_name,5) "type(ne2d_t)"
    write (iou,6) array_name,array%get_units()
    !
    call array%set( storage_size(nodedusp_dt,kind=inttype) , &
                            size(nodedusp_dt,kind=inttype) )
    !
    write (array_name,5) "nodedusp_dt(:)"
    write (iou,6) array_name,array%get_units()
    !
    call dt%set(array)
    !
    call array%reset
    !
    do i = lbound(nodedusp_dt,dim=1),ubound(nodedusp_dt,dim=1)
      !
      call array%add( storage_size(nodedusp_dt(i)%v,kind=inttype) , &
                              size(nodedusp_dt(i)%v,kind=inttype) )
      !
    end do
    !
    write (array_name,5) "nodedusp_dt(:)%v"
    write (iou,6) array_name,array%get_units()
    !
    call dt%add(array)
    !
    call total%add(dt)
    write (array_name,5) "nodedusp_dt"
    write (iou,2) array_name,dt%get_units()
    !
  else
    write (iou,7) "nodedusp_dt"
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
  if (allocated(gpsp)) then
    call array%set( storage_size(gpsp,kind=inttype) , &
                            size(gpsp,kind=inttype) )
    call total%add(array)
    write (array_name,5) "gpsp"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "gpsp"
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
  if (allocated(rk_sum)) then
    call array%set( storage_size(rk_sum,kind=inttype) , &
                            size(rk_sum,kind=inttype) )
    call total%add(array)
    write (array_name,5) "rk_sum"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "rk_sum"
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
  if (allocated(faceusp_array)) then
    call array%set( storage_size(faceusp_array,kind=inttype) , &
                            size(faceusp_array,kind=inttype) )
    call total%add(array)
    write (array_name,5) "faceusp_array"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "faceusp_array"
  end if
  if (allocated(facedusp_array)) then
    call array%set( storage_size(facedusp_array,kind=inttype) , &
                            size(facedusp_array,kind=inttype) )
    call total%add(array)
    write (array_name,5) "facedusp_array"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "facedusp_array"
  end if
  if (allocated(faceflx)) then
    call array%set( storage_size(faceflx,kind=inttype) , &
                            size(faceflx,kind=inttype) )
    call total%add(array)
    write (array_name,5) "faceflx"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "faceflx"
  end if
  if (allocated(facegpsp)) then
    call array%set( storage_size(facegpsp,kind=inttype) , &
                            size(facegpsp,kind=inttype) )
    call total%add(array)
    write (array_name,5) "facegpsp"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "facegpsp"
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
