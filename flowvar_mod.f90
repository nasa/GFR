module flowvar
  !
  use module_kind_types, only : wp,lk
  use module_kind_types, only : zero,eps1,ten2,iout
  use module_kind_types, only : true,fals
  use module_kind_types, only : debug_timer
  use module_kind_types, only : entering_procedure
  use module_kind_types, only : leaving_procedure
  use module_kind_types, only : Classic_RK
  use module_kind_types, only : CK4_RK
  use module_kind_types, only : TVD4_RK
  use module_kind_types, only : OERK
  use module_kind_types, only : EPERK
  use module_kind_types, only : alloc_error
  use module_kind_types, only : skip_alloc_pause
  use module_kind_types, only : error_message
  use module_kind_types, only : Geom_Edge
  use module_kind_types, only : bc_sem
  use module_kind_types, only : start_timer
  use module_kind_types, only : stop_timer
  use module_kind_types, only : mypnum
  use module_kind_types, only : glb_root
  use module_kind_types, only : stop_gfr
  use module_kind_types, only : stop_mpi
  use module_kind_types, only : THX_Phase_3
  use module_kind_types, only : bc_walls
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
  real(wp), public, save, allocatable, dimension(:,:) :: gpoldsp
  !$ACC DECLARE CREATE(gpoldsp)
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
  use eqn_idx,   only : nq,ng,npres,ntemp,ngamm,ncs
  use geovar,    only : nr,n_solpts,n_totpts
  use geovar,    only : nfbnd,nface,face
  use geovar,    only : nnode,nedge,edge
  use order_mod, only : geom_solpts
  use ovar,      only : Runge_Kutta_Scheme,num_rk_stages,mms_opt
  use ovar,      only : convert_restart_to_cgns
  use ovar,      only : continuous_output
  use ovar,      only : time_average_restart_files
  use ovar,      only : if_interpto_newgrid
  use ovar,      only : bc_in
  use ovar,      only : itestcase
  !
  !.. Local Scalars ..
  integer :: i,iwt,l,m,n,ne,nf,nn,nepts,nfp,nfpts
  integer :: face_geom,face_order,ierr
  character(len=200) :: array_name
  logical(lk) :: need_facexyz
  !
  !.. Local Arrays ..
  real(wp)           :: wall_time(1:50)
  character(len=100) :: message(1:50)
  !
  !.. Local Constants ..
  character(len=*), parameter :: pname = "allocate_flowvar_arrays"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  iwt = 0
  wall_time(:) = zero
  message(:) = ""
  !
  iwt = iwt + 1 ; message(iwt) = "Allocating usp, dusp, and gpsp"
call debug_timer(start_timer,message(iwt))
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
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
  !
  ! If we are only time averaging restart files, we dont need any of
  ! these other arrays so return out of here
  !
  if (time_average_restart_files) then
    call debug_timer(leaving_procedure,pname)
    return
  end if
  !
  if (any(abs(bc_in(:)%rotation_speed) > eps1)) then
    need_facexyz = fals
    do i = lbound(bc_in,dim=1),ubound(bc_in,dim=1)
      if (abs(bc_in(i)%rotation_speed) <= eps1) cycle
      if (any(bc_in(i)%bc_type == bc_walls)) then
        need_facexyz = true
        exit
      end if
    end do
    if (need_facexyz) then
  iwt = iwt + 1 ; message(iwt) = "Allocating facexyz_dt for rotating wall"
call debug_timer(start_timer,message(iwt))
      if (.not. allocated(facexyz_dt)) then
        allocate ( facexyz_dt(1:nfbnd) , stat=ierr , errmsg=error_message )
       !allocate ( facexyz_dt(1:nfbnd) , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"facexyz_dt",1,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
    end if
  end if
  !
  if (any(bc_in(:)%bc_type == bc_sem) .or. &
      abs(itestcase) == THX_Phase_3) then
  iwt = iwt + 1 ; message(iwt) = "Allocating facexyz_dt for SEM"
call debug_timer(start_timer,message(iwt))
    if (.not. allocated(facexyz_dt)) then
      allocate ( facexyz_dt(1:nface) , stat=ierr , errmsg=error_message )
     !allocate ( facexyz_dt(1:nfbnd) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"facexyz_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
  end if
  !
  if (ncs > 0) then
  iwt = iwt + 1 ; message(iwt) = "Allocating facexyz_dt for ncs > 0"
call debug_timer(start_timer,message(iwt))
    if (allocated(facexyz_dt)) then
      if (size(facexyz_dt) < nface) then
        deallocate ( facexyz_dt , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"facexyz_dt",2,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
    end if
    if (.not. allocated(facexyz_dt)) then
      allocate ( facexyz_dt(1:nface) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"facexyz_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
  end if
  !
  !
  ! If we are only interpolating the solution to a new grid,
  ! then only need facexyz_dt before returning out of here.
  !
  if ( if_interpto_newgrid() ) then
  iwt = iwt + 1 ; message(iwt) = "Allocating arrays for interpolation utility"
call debug_timer(start_timer,message(iwt))
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
    if (.not. allocated(facexyz_dt)) then
      allocate ( facexyz_dt(1:nface) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"facexyz_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
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
    allocate ( facegpsp_dt(1:nface) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facegpsp_dt",1,__LINE__,__FILE__, &
                     ierr,error_message)
    !
    call alloc_face_dt_components(facegpsp_dt,nfpts,"facegpsp_dt",num_var=ng)
    !
    allocate ( facegpsp(1:nfpts,1:ng,1:2) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facegpsp",1,__LINE__,__FILE__, &
                     ierr,error_message)
    !
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
    !
    call debug_timer(leaving_procedure,pname)
    !
    return
    !
  end if
  !
  !
  iwt = iwt + 1
  message(iwt) = "Allocating residual, dtsp, uoldsp, and gpoldsp"
call debug_timer(start_timer,message(iwt))
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
  allocate ( gpoldsp(1:n_solpts,1:max(ntemp,npres,ngamm)) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"gpoldsp",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(gpoldsp)
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
  !
  !
  if (convert_restart_to_cgns) then
  iwt = iwt + 1
  message(iwt) = "Allocating faceusp_array for convert_restart_to_cgns"
call debug_timer(start_timer,message(iwt))
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
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
  end if
  !
  ! Allocate the outer part of the face solution/flux arrays
  !
 !if (.not. convert_restart_to_cgns .or. continuous_output) then
    !
    ! If using the continuous output option or MMS,
    ! allocate the facexyz_dt array
    !
    if (continuous_output) then
      !
      ! If using continuous output, facexyz_dt needs to be defined
      ! for both interior and boundary faces
      !
  iwt = iwt + 1 
  message(iwt) = "Allocating face{xyz,usp,gpsp,dusp}_dt for continuous output"
call debug_timer(start_timer,message(iwt))
      if (allocated(facexyz_dt)) then
        if (size(facexyz_dt) < nface) then
          deallocate ( facexyz_dt , stat=ierr , errmsg=error_message )
          call alloc_error(pname,"facexyz_dt",2,__LINE__,__FILE__, &
                           ierr,error_message)
        end if
      end if
      if (.not. allocated(facexyz_dt)) then
        allocate ( facexyz_dt(1:nface) , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"facexyz_dt",1,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
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
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
      !
    else if (mms_opt /= 0) then
      !
      ! If using MMS but not using continuous output, facexyz_dt only needs
      ! to be defined for the boundary faces
      !
  iwt = iwt + 1 ; message(iwt) = "Allocating facexyz_dt for MMS"
call debug_timer(start_timer,message(iwt))
      if (.not. allocated(facexyz_dt)) then
        allocate ( facexyz_dt(1:nfbnd) , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"facexyz_dt",1,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
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
  iwt = iwt + 1 ; message(iwt) = "Allocating faceusp_dt components"
call debug_timer(start_timer,message(iwt))
      call alloc_face_dt_components(faceusp_dt,nfpts,"faceusp_dt")
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
    end if
    !
    if (allocated(facegpsp_dt)) then
  iwt = iwt + 1 ; message(iwt) = "Allocating facegpsp_dt components"
call debug_timer(start_timer,message(iwt))
      call alloc_face_dt_components(facegpsp_dt,nfpts,"facegpsp_dt",num_var=ng)
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
    end if
    !
    if (allocated(facedusp_dt)) then
  iwt = iwt + 1 ; message(iwt) = "Allocating facedusp_dt components"
call debug_timer(start_timer,message(iwt))
      call alloc_face_dt_components(facedusp_dt,nfpts,"facedusp_dt")
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
    end if
    !
    !
  iwt = iwt + 1
  message(iwt) = "Allocating faceflx and face{usp,gpsp,dusp}_array"
call debug_timer(start_timer,message(iwt))
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
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
    !
    !
    if (allocated(facexyz_dt)) then
      !
  iwt = iwt + 1 ; message(iwt) = "Allocating facexyz_dt components"
call debug_timer(start_timer,message(iwt))
      call alloc_face_dt_components(facexyz_dt,nfpts,"facexyz_dt")
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
      !
  iwt = iwt + 1 ; message(iwt) = "Allocating facexyz_array"
call debug_timer(start_timer,message(iwt))
      allocate ( facexyz_array(1:nr,1:nfpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"facexyz_array",1,__LINE__,__FILE__, &
                       ierr,error_message)
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
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
  iwt = iwt + 1 ; message(iwt) = "Allocating nodeusp_dt and its components"
call debug_timer(start_timer,message(iwt))
      allocate ( nodeusp_dt(1:nnode) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"nodeusp_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      call alloc_ne_dt_components(nodeusp_dt,"nodeusp_dt")
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
      !
      ! Allocate the nodegpsp_dt array and its data components
      !
  iwt = iwt + 1 ; message(iwt) = "Allocating nodegpsp_dt and its components"
call debug_timer(start_timer,message(iwt))
      allocate ( nodegpsp_dt(1:nnode) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"nodegpsp_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      call alloc_ne_dt_components(nodegpsp_dt,"nodegpsp_dt",num_var=ng)
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
      !
      ! Allocate the nodedusp_dt array and its data components
      !
  iwt = iwt + 1 ; message(iwt) = "Allocating nodedusp_dt and its components"
call debug_timer(start_timer,message(iwt))
      allocate ( nodedusp_dt(1:nnode) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"nodedusp_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      call alloc_ne_dt_components(nodedusp_dt,"nodedusp_dt")
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
      !
      ! Allocate the edgeusp_dt array and its data components
      !
  iwt = iwt + 1 ; message(iwt) = "Allocating edgeusp_dt and its components"
call debug_timer(start_timer,message(iwt))
      allocate ( edgeusp_dt(1:nedge) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"edgeusp_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      call alloc_ne_dt_components(edgeusp_dt,nepts,"edgeusp_dt")
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
      !
      ! Allocate the edgegpsp_dt array and its data components
      !
  iwt = iwt + 1 ; message(iwt) = "Allocating edgegpsp_dt and its components"
call debug_timer(start_timer,message(iwt))
      allocate ( edgegpsp_dt(1:nedge) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"edgegpsp_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      call alloc_ne_dt_components(edgegpsp_dt,nepts,"edgegpsp_dt",num_var=ng)
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
      !
      ! Allocate the edgedusp_dt array and its data components
      !
  iwt = iwt + 1 ; message(iwt) = "Allocating edgedusp_dt and its components"
call debug_timer(start_timer,message(iwt))
      allocate ( edgedusp_dt(1:nedge) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"edgedusp_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      call alloc_ne_dt_components(edgedusp_dt,nepts,"edgedusp_dt")
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
      !
      ! Allocate the edgexyz_dt array and its data components
      !
  iwt = iwt + 1 ; message(iwt) = "Allocating edgexyz_dt and its components"
call debug_timer(start_timer,message(iwt))
      allocate ( edgexyz_dt(1:nedge) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"edgexyz_dt",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      call alloc_ne_dt_components(edgexyz_dt,nepts,"edgexyz_dt",num_var=nr)
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
      !
    end if
    !
 !end if
  !
  ! Allocate the storage for the Runge-Kutta methods
  !
  iwt = iwt + 1 ; message(iwt) = "Allocating Runge Kutta arrays"
call debug_timer(start_timer,message(iwt))
  if ( any(Runge_Kutta_Scheme == [Classic_RK,-EPERK]) ) then
    allocate ( rssprk(1:n_solpts,1:nq,1:num_rk_stages) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"rssprk",1,__LINE__,__FILE__,ierr,error_message)
    !$ACC UPDATE DEVICE(rssprk)
  end if
  !
  if ( any(Runge_Kutta_Scheme == [CK4_RK,TVD4_RK,OERK,EPERK]) ) then
    allocate ( rk_sum(1:n_solpts,1:nq) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"rk_sum",1,__LINE__,__FILE__,ierr,error_message)
    !$ACC UPDATE DEVICE(rk_sum)
  end if
call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
  !
  if (mypnum == glb_root) then
    call sort_and_output_wall_times(wall_time(1:iwt),message(1:iwt))
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
subroutine sort_and_output_wall_times(wall_time,message)
  !
  !.. Formal Arguments ..
  real(wp),         intent(in) :: wall_time(:)
  character(len=*), intent(in) :: message(:)
  !
  !.. Local Scalars ..
  integer  :: i,n,nwt
  real(wp) :: total_time
  !
  !.. Local Arrays ..
  integer  :: sort_order(1:size(wall_time))
  real(wp) :: temp_time(1:size(wall_time))
  real(wp) :: percent_time(1:size(wall_time))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "sort_and_output_wall_times"
  !
continue
  !
#ifdef PROFILE_ON
  nwt = size(wall_time)
  !
  do n = 1,nwt
    temp_time(n) = max(zero,wall_time(n))
  end do
  !
  total_time = sum(temp_time(:))
  !
  do n = 1,nwt
    !
    i = minloc(temp_time,dim=1)
    !
    sort_order(n) = i
    !
    percent_time(i) = ten2 * temp_time(i) / total_time
    !
    temp_time(i) = huge(zero)
    !
  end do
  !
  write (iout,1)
  do n = 1,nwt
    i = sort_order(n)
    write (iout,2) wall_time(i),percent_time(i),trim(adjustl(message(i)))
  end do
  write (iout,*)
  flush (iout)
#endif
  !
  ! Format Statements
  !
  1 format (/,10x,"Sorted allocation times from allocate_flowvar_arrays")
  2 format (12x,f12.6," ( ",f5.2,"% ) seconds for ",a)
  !
end subroutine sort_and_output_wall_times
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
#ifndef DEBUG_ON
    write (array_name,1) trim(adjustl(aname)),nf
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
#endif
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
#ifndef DEBUG_ON
    write (array_name,1) trim(adjustl(aname)),nf
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
#endif
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
#ifndef DEBUG_ON
    write (array_name,1) trim(adjustl(aname)),nf
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
#endif
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
#ifndef DEBUG_ON
    write (array_name,1) trim(adjustl(aname)),nn
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
#endif
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
#ifndef DEBUG_ON
    write (array_name,1) trim(adjustl(aname)),nn
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
#endif
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
#ifndef DEBUG_ON
    write (array_name,1) trim(adjustl(aname)),ne
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
#endif
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
#ifndef DEBUG_ON
    write (array_name,1) trim(adjustl(aname)),ne
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
#endif
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
  if (allocated(gpoldsp)) then
    call array%set( storage_size(gpoldsp,kind=inttype) , &
                            size(gpoldsp,kind=inttype) )
    call total%add(array)
    write (array_name,5) "gpoldsp"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "gpoldsp"
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
