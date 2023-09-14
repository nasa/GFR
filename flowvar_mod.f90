module flowvar
  !
  use module_kind_types, only : wp
  use module_kind_types, only : zero
  use module_kind_types, only : debug_timer
  use module_kind_types, only : entering_procedure
  use module_kind_types, only : leaving_procedure
  use module_kind_types, only : alloc_error
  use module_kind_types, only : error_message
  use module_kind_types, only : start_timer
  use module_kind_types, only : stop_timer
  use module_kind_types, only : mypnum
  use module_kind_types, only : glb_root
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
  ! FACEUSP  : conservative variables at the interface flx-pts
  ! FACEDUSP : conservative variable gradients at the interface flx-pts
  !
  real(wp), public, save, allocatable :: faceusp(:,:,:)
  !$ACC DECLARE CREATE(faceusp)
  real(wp), public, save, allocatable :: facedusp(:,:,:,:)
  !$ACC DECLARE CREATE(facedusp)
  !
  ! FACEFLX  : flux normal to cell faces at each interface flux point
  !
  real(wp), public, save, allocatable :: faceflx(:,:)
  !$ACC DECLARE CREATE(faceflx)
  !
  ! FACEXYZ  : coordinates of the interface flux points
  !
  real(wp), public, save, allocatable :: facexyz(:,:)
  !$ACC DECLARE CREATE(facexyz)
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
  interface alloc_array
    module procedure alloc_rank1_array, &
                     alloc_rank2_array, &
                     alloc_rank3_array, &
                     alloc_rank4_array
  end interface alloc_array
  !
  interface alloc_fen_dt
    module procedure alloc_face_xyz_dt, &
                     alloc_face_var_dt, &
                     alloc_face_vec_dt, &
                     alloc_ne1d_dt, &
                     alloc_ne2d_dt, &
                     alloc_ne3d_dt
  end interface alloc_fen_dt
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
  integer,            save :: iwt = 0
  real(wp),           save :: wall_time(1:50)
  character(len=100), save :: message(1:50)
  !
contains
!
!###############################################################################
!
subroutine allocate_flowvar_arrays()
  !
  !.. Use Statements ..
  use module_kind_types, only : Classic_RK,CK4_RK,TVD4_RK,OERK,EPERK
  use eqn_idx,           only : nq,ng,npres,ntemp,ngamm
  use geovar,            only : nr,n_solpts,n_flxpts
  use geovar,            only : nface,nedge,nnode
  use ovar,              only : Runge_Kutta_Scheme,num_rk_stages
  use ovar,              only : continuous_output
  use ovar,              only : time_average_restart_files
  use ovar,              only : if_interpto_newgrid
  !
  !.. Local Scalars ..
  integer :: nvar,nf,ne,nn,nfpts,nepts,nspts
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
  ! For the special case of interpolating the solution to a new grid
  !
  if ( if_interpto_newgrid() ) then
    call alloc_interpto_arrays()
    call debug_timer(leaving_procedure,pname)
    return
  end if
  !
  ! Allocate usp, dusp, and gpsp since these should always be needed
  !
  call alloc_array(usp,n_solpts,nq,"usp",pname,__LINE__)
  !
  call alloc_array(dusp,n_solpts,nr,nq,"dusp",pname,__LINE__)
  !
  call alloc_array(gpsp,n_solpts,ng,"gpsp",pname,__LINE__)
  !
  ! If we are only time averaging restart files, we dont need any of
  ! these other arrays so return out of here
  !
  if (time_average_restart_files) then
    call debug_timer(leaving_procedure,pname)
    return
  end if
  !
  ! Allocate the other solution point arrays needed by the solver
  !
  call alloc_array(residual,n_solpts,nq,"residual",pname,__LINE__)
  !
  call alloc_array(dtsp,n_solpts,"dtsp",pname,__LINE__)
  !
  call alloc_array(uoldsp,n_solpts,nq,"uoldsp",pname,__LINE__)
  !
  ! NOTE: We dont want to use ng for the number of variables in gpoldsp
  !       because we only care about the values of temp/pres/gamma from the
  !       previous time step. The number of variables in gpsp, ng, can also
  !       include viscosity and thermal conductivity (for viscous flows) and
  !       possibly enthalpy and cp (for multi-species gas / reacting flows).
  !       There is no need to waste memory here, so just allocate gpoldsp
  !       to the max variable index out of temp/pres/gamma. This should be
  !       gamma, since (unless something has changed) these indices default
  !       to [ntemp,npres,ngamm]=[1,2,3].
  nvar = max(ntemp,npres,ngamm)
  call alloc_array(gpoldsp,n_solpts,nvar,"gpoldsp",pname,__LINE__)
  !
  ! Allocate the face point arrays needed by the solver 
  !
  call alloc_array(faceflx,n_flxpts,nq,"faceflx",pname,__LINE__)
  !
  call alloc_array(faceusp,n_flxpts,nq,2,"faceusp",pname,__LINE__)
  !
  call alloc_array(facedusp,n_flxpts,nr,nq,2,"facedusp",pname,__LINE__)
  !
  call alloc_array(facegpsp,n_flxpts,ng,2,"facegpsp",pname,__LINE__)
  !
  ! Allocate some extra storage needed for certain Runge-Kutta schemes
  !
  nspts = merge(n_solpts,0,any(Runge_Kutta_Scheme == [Classic_RK,-EPERK]))
  call alloc_array(rssprk,nspts,nq,num_rk_stages,"rssprk",pname,__LINE__)
  !
  nspts = merge(n_solpts,0,any(Runge_Kutta_Scheme==[CK4_RK,TVD4_RK,OERK,EPERK]))
  call alloc_array(rk_sum,nspts,nq,"rk_sum",pname,__LINE__)
  !
  ! Find the size needed for the facexyz array
  !
  call find_facexyz_size(nf,nfpts)
  !
  ! Allocate the facexyz and facexyz_dt arrays
  !
  call alloc_array(facexyz,nr,nfpts,"facexyz",pname,__LINE__)
  !
  call alloc_fen_dt(facexyz_dt,nf,"facexyz_dt",pname,__LINE__)
  !
  ! Allocate the face/edge/node derived type arrays
  !
  nf = 0
  if (continuous_output) nf = nface
  !
  call alloc_fen_dt(faceusp_dt,nf,"faceusp_dt",pname,__LINE__)
  !
  call alloc_fen_dt(facedusp_dt,nf,"facedusp_dt",pname,__LINE__)
  !
  call alloc_fen_dt(facegpsp_dt,nf,"facegpsp_dt",pname,__LINE__,num_var=ng)
  !
  ! Allocate the edge derived type arrays
  !
  ne = 0
  if (continuous_output .and. nr>=3) ne = nedge
  !
  call alloc_fen_dt(edgexyz_dt,ne,"edgexyz_dt",pname,__LINE__, &
                    num_var=nr,nepts=nepts)
  !
  call alloc_fen_dt(edgeusp_dt,ne,"edgeusp_dt",pname,__LINE__,nepts=nepts)
  !
  call alloc_fen_dt(edgedusp_dt,ne,"edgedusp_dt",pname,__LINE__)
  !
  call alloc_fen_dt(edgegpsp_dt,ne,"edgegpsp_dt",pname,__LINE__, &
                    num_var=ng,nepts=nepts)
  !
  ! Allocate the node derived type arrays
  !
  nn = 0
  if (continuous_output) nn = nnode
  !
  call alloc_fen_dt(nodeusp_dt,nf,"nodeusp_dt",pname,__LINE__)
  !
  call alloc_fen_dt(nodedusp_dt,nf,"nodedusp_dt",pname,__LINE__)
  !
  call alloc_fen_dt(nodegpsp_dt,nf,"nodegpsp_dt",pname,__LINE__,num_var=ng)
  !
  ! Finally, sort and output the time taken to allocate all these arrays
  !
  if (mypnum == glb_root) then
    call sort_and_output_wall_times(wall_time(1:iwt),message(1:iwt))
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine allocate_flowvar_arrays
!
!###############################################################################
!
subroutine alloc_interpto_arrays()
  !
  !.. Use Statements ..
  use eqn_idx, only : nq,ng
  use geovar,  only : nr,n_solpts,n_flxpts,nface
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "alloc_interpto_arrays"
  !
continue
  !
  ! Allocate all the arrays necessary for interpolating a solution to a new grid
  !
  call debug_timer(entering_procedure,pname)
  !
  !##################################
  !###  1. Allocate usp and gpsp  ###
  !##################################
  !
  call alloc_array(usp,n_solpts,nq,"usp",pname,__LINE__)
  !
  call alloc_array(gpsp,n_solpts,ng,"gpsp",pname,__LINE__)
  !
  !#############################################################
  !###  2. The interpolation utility needs facexyz to  ###
  !###     be defined for both interior and boundary faces.  ###
  !#############################################################
  !
  ! Allocate the derived type array facexyz_dt and
  ! the derived type components for each face
  !
  call alloc_fen_dt(facexyz_dt,nface,"facexyz_dt",pname,__LINE__)
  !
  ! Allocate the traditional array version of facexyz_dt
  !
  call alloc_array(facexyz,nr,n_flxpts,"facexyz",pname,__LINE__)
  !
  !#####################################################################
  !###  3. Procedure interpsol needs faceusp_dt to print debug info  ###
  !#####################################################################
  !
  ! Allocate the derived type array faceusp_dt and
  ! the derived type components for each face
  !
  call alloc_fen_dt(faceusp_dt,nface,"faceusp_dt",pname,__LINE__)
  !
  ! Allocate the traditional array version of faceusp_dt
  !
  call alloc_array(faceusp,n_flxpts,nq,2,"faceusp",pname,__LINE__)
  !
  !#################################################
  !###  4. The facegpsp_dt array is also needed  ###
  !#################################################
  !
  ! Allocate the derived type array facegpsp_dt and
  ! the derived type components for each face
  !
  call alloc_fen_dt(facegpsp_dt,nface,"facegpsp_dt",pname,__LINE__,num_var=ng)
  !
  ! Allocate the traditional array version of facegpsp_dt
  ! NOTE: nfp should be equal to n_flxpts
  !
  call alloc_array(facegpsp,n_flxpts,ng,2,"facegpsp",pname,__LINE__)
  !
  !########################################################################
  !###  5. Sort and output the time taken to allocate all these arrays  ###
  !########################################################################
  !
  if (mypnum == glb_root) then
    call sort_and_output_wall_times(wall_time(1:iwt),message(1:iwt))
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine alloc_interpto_arrays
!
!###############################################################################
!
subroutine find_facexyz_size(nf,nfpts)
  !
  !.. Use Statements ..
  use module_kind_types, only : one,rpm_to_mps,bc_sem,bc_walls,THX_Phase_3
  use eqn_idx,           only : ncs
  use geovar,            only : nfbnd,nface,n_bnd_flxpts,n_flxpts
  use ovar,              only : bc_in,moving_ref_frame,aref
  use ovar,              only : itestcase,continuous_output,mms_opt
  !
  !.. Formal Arguments ..
  integer, intent(inout) :: nf,nfpts
  !
  !.. Local Scalars ..
  integer :: i
  !
continue
  !
  ! Initialize nf and nfpts to 0, so that facexyz and facexyz_dt are only
  ! allocated if nf and nfpts are modified which would mean that the arrays
  ! are actually needed.
  !
  nf = 0
  nfpts = 0
  !
  ! If rotation_rate for a wall BC is nonzero, we need
  ! facexyz for at least the boundary face points.
  !
  rotating_wall_bc_check: &
  do i = lbound(bc_in,dim=1),ubound(bc_in,dim=1)
    if (any(bc_in(i)%bc_type == bc_walls)) then
      if (any(abs(bc_in(i)%rotation_rate(:)) >= one*rpm_to_mps/aref)) then
        nf = max(nf,nfbnd)
        nfpts = max(nfpts,n_bnd_flxpts)
        exit rotating_wall_bc_check
      end if
    end if
  end do rotating_wall_bc_check
  !
  ! If using an SEM BC, we need facexyz for at least the boundary face points.
  ! NOTE: facexyz can be deallocated after initializing the SEM module for
  !       assuming it isnt needed for another reason.
  !
  if (any(bc_in(:)%bc_type == bc_sem)) then
    nf = max(nf,nfbnd)
    nfpts = max(nfpts,n_bnd_flxpts)
  end if
  !
  ! If using MMS, we need facexyz for at least the boundary face points.
  !
  if (mms_opt /= 0) then
    nf = max(nf,nfbnd)
    nfpts = max(nfpts,n_bnd_flxpts)
  end if
  !
  ! If using a rotating reference frame, we need facexyz for all face points.
  !
  if ( moving_ref_frame%is_active() ) then
    nf = nface
    nfpts = n_flxpts
    ! Just return since these are the max values possible
    return
  end if
  !
  ! If running the THX3 case, we need facexyz for all face points.
  ! NOTE: facexyz can be deallocated after initializing the code
  !       for computing the blowing ratio for the THX3 case,
  !       assuming it isnt needed for another reason.
  !
  if (abs(itestcase) == THX_Phase_3) then
    nf = nface
    nfpts = n_flxpts
    ! Just return since these are the max values possible
    return
  end if
  !
  ! If running a multi-species gas flow (with or without reactions), we need
  ! facexyz for all face points in case the Newton solver for temperature 
  ! fails to converge, so that the coordinates where it failed can be output
  ! to help with troubleshooting.
  !
  if (ncs > 0) then
    nf = nface
    nfpts = n_flxpts
    ! Just return since these are the max values possible
    return
  end if
  !
  ! If using continuous output, we need facexyz for all face points.
  !
  if (continuous_output) then
    nf = nface
    nfpts = n_flxpts
    ! Just return since these are the max values possible
    return
  end if
  !
end subroutine find_facexyz_size
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine alloc_rank1_array(array,n1,aname,pname,line_num,imsg)
  !
  !.. Formal Arguments ..
  real(wp), allocatable, intent(inout) :: array(:)
  integer,                  intent(in) :: n1,line_num
  character(len=*),         intent(in) :: aname,pname
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: imsg
  !
  !.. Local Scalars ..
  integer :: ierr
  !
continue
  !
  ! If n1 is <= 0, deallocate array and return
  !
  if (n1 <= 0) then
    if (allocated(array)) then
      deallocate ( array , stat=ierr , errmsg=error_message )
      call alloc_error(pname,aname,2,line_num,__FILE__,ierr,error_message)
    end if
    return
  end if
  !
  iwt = iwt + 1
  message(iwt) = "Allocating " // trim(adjustl(aname))
  if (present(imsg)) then
    message(iwt) = trim(adjustl(message(iwt))) // " " // trim(adjustl(imsg))
  end if
  call debug_timer(start_timer,message(iwt))
  !
  if (.not. allocated(array)) then
    allocate ( array(1:n1) , source=zero , stat=ierr , errmsg=error_message )
    call alloc_error(pname,aname,1,line_num,__FILE__,ierr,error_message)
    !$ACC UPDATE DEVICE(array)
  end if
  !
  call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
  !
end subroutine alloc_rank1_array
!
!###############################################################################
!
subroutine alloc_rank2_array(array,n1,n2,aname,pname,line_num,imsg)
  !
  !.. Formal Arguments ..
  real(wp), allocatable, intent(inout) :: array(:,:)
  integer,                  intent(in) :: n1,n2,line_num
  character(len=*),         intent(in) :: aname,pname
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: imsg
  !
  !.. Local Scalars ..
  integer :: ierr
  !
continue
  !
  ! If any of [n1,n2] is <= 0, deallocate array and return
  !
  if (any([n1,n2] <= 0)) then
    if (allocated(array)) then
      deallocate ( array , stat=ierr , errmsg=error_message )
      call alloc_error(pname,aname,2,line_num,__FILE__,ierr,error_message)
    end if
    return
  end if
  !
  iwt = iwt + 1
  message(iwt) = "Allocating " // trim(adjustl(aname))
  if (present(imsg)) then
    message(iwt) = trim(adjustl(message(iwt))) // " " // trim(adjustl(imsg))
  end if
  call debug_timer(start_timer,message(iwt))
  !
  if (.not. allocated(array)) then
    allocate ( array(1:n1,1:n2) , source=zero, stat=ierr, errmsg=error_message )
    call alloc_error(pname,aname,1,line_num,__FILE__,ierr,error_message)
    !$ACC UPDATE DEVICE(array)
  end if
  !
  call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
  !
end subroutine alloc_rank2_array
!
!###############################################################################
!
subroutine alloc_rank3_array(array,n1,n2,n3,aname,pname,line_num,imsg)
  !
  !.. Formal Arguments ..
  real(wp), allocatable, intent(inout) :: array(:,:,:)
  integer,                  intent(in) :: n1,n2,n3,line_num
  character(len=*),         intent(in) :: aname,pname
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: imsg
  !
  !.. Local Scalars ..
  integer :: ierr
  !
continue
  !
  ! If any of [n1,n2,n3] is <= 0, deallocate array and return
  !
  if (any([n1,n2,n3] <= 0)) then
    if (allocated(array)) then
      deallocate ( array , stat=ierr , errmsg=error_message )
      call alloc_error(pname,aname,2,line_num,__FILE__,ierr,error_message)
    end if
    return
  end if
  !
  iwt = iwt + 1
  message(iwt) = "Allocating " // trim(adjustl(aname))
  if (present(imsg)) then
    message(iwt) = trim(adjustl(message(iwt))) // " " // trim(adjustl(imsg))
  end if
  call debug_timer(start_timer,message(iwt))
  !
  if (.not. allocated(array)) then
    allocate ( array(1:n1,1:n2,1:n3) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,aname,1,line_num,__FILE__,ierr,error_message)
    !$ACC UPDATE DEVICE(array)
  end if
  !
  call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
  !
end subroutine alloc_rank3_array
!
!###############################################################################
!
subroutine alloc_rank4_array(array,n1,n2,n3,n4,aname,pname,line_num,imsg)
  !
  !.. Formal Arguments ..
  real(wp), allocatable, intent(inout) :: array(:,:,:,:)
  integer,                  intent(in) :: n1,n2,n3,n4,line_num
  character(len=*),         intent(in) :: aname,pname
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: imsg
  !
  !.. Local Scalars ..
  integer :: ierr
  !
continue
  !
  ! If any of [n1,n2,n3,n4] is <= 0, deallocate array and return
  !
  if (any([n1,n2,n3,n4] <= 0)) then
    if (allocated(array)) then
      deallocate ( array , stat=ierr , errmsg=error_message )
      call alloc_error(pname,aname,2,line_num,__FILE__,ierr,error_message)
    end if
    return
  end if
  !
  iwt = iwt + 1
  message(iwt) = "Allocating " // trim(adjustl(aname))
  if (present(imsg)) then
    message(iwt) = trim(adjustl(message(iwt))) // " " // trim(adjustl(imsg))
  end if
  call debug_timer(start_timer,message(iwt))
  !
  if (.not. allocated(array)) then
    allocate ( array(1:n1,1:n2,1:n3,1:n4) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,aname,1,line_num,__FILE__,ierr,error_message)
    !$ACC UPDATE DEVICE(array)
  end if
  !
  call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
  !
end subroutine alloc_rank4_array
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine alloc_face_xyz_dt(dt_array,n1,aname,pname,line_num,imsg)
  !
  !.. Formal Arguments ..
  type(face_xyz), allocatable, intent(inout) :: dt_array(:)
  integer,                        intent(in) :: n1,line_num
  character(len=*),               intent(in) :: aname,pname
  !
  !.. Optional Arguments ..
  character(len=*), optional, intent(in) :: imsg
  !
  !.. Local Scalars ..
  integer :: ierr
  !
continue
  !
  ! If n1 is <= 0, deallocate dt_array and return
  !
  if (n1 <= 0) then
    if (allocated(dt_array)) then
      deallocate ( dt_array , stat=ierr , errmsg=error_message )
      call alloc_error(pname,aname,2,line_num,__FILE__,ierr,error_message)
    end if
    return
  end if
  !
  iwt = iwt + 1
  message(iwt) = "Allocating " // trim(adjustl(aname)) // " and its components"
  if (present(imsg)) then
    message(iwt) = trim(adjustl(message(iwt))) // " " // trim(adjustl(imsg))
  end if
  call debug_timer(start_timer,message(iwt))
  !
  ! First allocate the derived type array
  !
  if (.not. allocated(dt_array)) then
    allocate ( dt_array(1:n1) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,aname,1,line_num,__FILE__,ierr,error_message)
  end if
  !
  ! Now allocate the components of each element of the array dt_array
  !
  call alloc_face_dt_components(dt_array,aname,pname,line_num)
  !
  call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
  !
end subroutine alloc_face_xyz_dt
!
!###############################################################################
!
subroutine alloc_face_var_dt(dt_array,n1,aname,pname,line_num,num_var,imsg)
  !
  !.. Formal Arguments ..
  type(face_var), allocatable, intent(inout) :: dt_array(:)
  integer,                        intent(in) :: n1,line_num
  character(len=*),               intent(in) :: aname,pname
  !
  !.. Optional Arguments ..
  integer,          optional, intent(in) :: num_var
  character(len=*), optional, intent(in) :: imsg
  !
  !.. Local Scalars ..
  integer :: ierr
  !
continue
  !
  ! If n1 is <= 0, deallocate dt_array and return
  !
  if (n1 <= 0) then
    if (allocated(dt_array)) then
      deallocate ( dt_array , stat=ierr , errmsg=error_message )
      call alloc_error(pname,aname,2,line_num,__FILE__,ierr,error_message)
    end if
    return
  end if
  !
  iwt = iwt + 1
  message(iwt) = "Allocating " // trim(adjustl(aname)) // " and its components"
  if (present(imsg)) then
    message(iwt) = trim(adjustl(message(iwt))) // " " // trim(adjustl(imsg))
  end if
  call debug_timer(start_timer,message(iwt))
  !
  ! First allocate the derived type array
  !
  if (.not. allocated(dt_array)) then
    allocate ( dt_array(1:n1) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,aname,1,line_num,__FILE__,ierr,error_message)
  end if
  !
  ! Now allocate the components of each element of the array dt_array
  !
  call alloc_face_dt_components(dt_array,aname,pname,line_num,num_var)
  !
  call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
  !
end subroutine alloc_face_var_dt
!
!###############################################################################
!
subroutine alloc_face_vec_dt(dt_array,n1,aname,pname,line_num,num_var,imsg)
  !
  !.. Formal Arguments ..
  type(face_vec), allocatable, intent(inout) :: dt_array(:)
  integer,                        intent(in) :: n1,line_num
  character(len=*),               intent(in) :: aname,pname
  !
  !.. Optional Arguments ..
  integer,          optional, intent(in) :: num_var
  character(len=*), optional, intent(in) :: imsg
  !
  !.. Local Scalars ..
  integer :: ierr
  !
continue
  !
  ! If n1 is <= 0, deallocate dt_array and return
  !
  if (n1 <= 0) then
    if (allocated(dt_array)) then
      deallocate ( dt_array , stat=ierr , errmsg=error_message )
      call alloc_error(pname,aname,2,line_num,__FILE__,ierr,error_message)
    end if
    return
  end if
  !
  iwt = iwt + 1
  message(iwt) = "Allocating " // trim(adjustl(aname)) // " and its components"
  if (present(imsg)) then
    message(iwt) = trim(adjustl(message(iwt))) // " " // trim(adjustl(imsg))
  end if
  call debug_timer(start_timer,message(iwt))
  !
  ! First allocate the derived type array
  !
  if (.not. allocated(dt_array)) then
    allocate ( dt_array(1:n1) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,aname,1,line_num,__FILE__,ierr,error_message)
  end if
  !
  ! Now allocate the components of each element of the array dt_array
  !
  call alloc_face_dt_components(dt_array,aname,pname,line_num,num_var)
  !
  call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
  !
end subroutine alloc_face_vec_dt
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine alloc_ne1d_dt(dt_array,n1,aname,pname,line_num,num_var,imsg)
  !
  !.. Formal Arguments ..
  type(ne1d_t), allocatable, intent(inout) :: dt_array(:)
  integer,                      intent(in) :: n1,line_num
  character(len=*),             intent(in) :: aname,pname
  !
  !.. Optional Arguments ..
  integer,          optional, intent(in) :: num_var
  character(len=*), optional, intent(in) :: imsg
  !
  !.. Local Scalars ..
  integer :: ierr
  !
continue
  !
  ! If n1 is <= 0, deallocate dt_array and return
  !
  if (n1 <= 0) then
    if (allocated(dt_array)) then
      deallocate ( dt_array , stat=ierr , errmsg=error_message )
      call alloc_error(pname,aname,2,line_num,__FILE__,ierr,error_message)
    end if
    return
  end if
  !
  iwt = iwt + 1
  message(iwt) = "Allocating " // trim(adjustl(aname)) // " and its components"
  if (present(imsg)) then
    message(iwt) = trim(adjustl(message(iwt))) // " " // trim(adjustl(imsg))
  end if
  call debug_timer(start_timer,message(iwt))
  !
  ! First allocate the derived type array
  !
  if (.not. allocated(dt_array)) then
    allocate ( dt_array(1:n1) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,aname,1,line_num,__FILE__,ierr,error_message)
  end if
  !
  ! Now allocate the components of each element of the array dt_array
  !
  call alloc_ne_dt_components(dt_array,aname,pname,line_num,num_var)
  !
  call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
  !
end subroutine alloc_ne1d_dt
!
!###############################################################################
!
subroutine alloc_ne2d_dt(dt_array,n1,aname,pname,line_num,num_var,nepts,imsg)
  !
  !.. Formal Arguments ..
  type(ne2d_t), allocatable, intent(inout) :: dt_array(:)
  integer,                      intent(in) :: n1,line_num
  character(len=*),             intent(in) :: aname,pname
  !
  !.. Optional Arguments ..
  integer,          optional,    intent(in) :: num_var
  integer,          optional, intent(inout) :: nepts
  character(len=*), optional,    intent(in) :: imsg
  !
  !.. Local Scalars ..
  integer :: ierr,nep
  !
continue
  !
  ! If n1 is <= 0, deallocate dt_array and return
  !
  if (n1 <= 0) then
    if (allocated(dt_array)) then
      deallocate ( dt_array , stat=ierr , errmsg=error_message )
      call alloc_error(pname,aname,2,line_num,__FILE__,ierr,error_message)
    end if
    return
  end if
  !
  iwt = iwt + 1
  message(iwt) = "Allocating " // trim(adjustl(aname)) // " and its components"
  if (present(imsg)) then
    message(iwt) = trim(adjustl(message(iwt))) // " " // trim(adjustl(imsg))
  end if
  call debug_timer(start_timer,message(iwt))
  !
  ! First allocate the derived type array
  !
  if (.not. allocated(dt_array)) then
    allocate ( dt_array(1:n1) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,aname,1,line_num,__FILE__,ierr,error_message)
  end if
  !
  ! Now allocate the components of each element of the array dt_array
  !
  if (present(nepts)) then
    call alloc_ne_dt_components(dt_array,nep,aname,pname,line_num,num_var)
    nepts = nep
  else
    call alloc_ne_dt_components(dt_array,aname,pname,line_num,num_var)
  end if
  !
  call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
  !
end subroutine alloc_ne2d_dt
!
!###############################################################################
!
subroutine alloc_ne3d_dt(dt_array,n1,aname,pname,line_num,num_var,imsg)
  !
  !.. Formal Arguments ..
  type(ne3d_t), allocatable, intent(inout) :: dt_array(:)
  integer,                      intent(in) :: n1,line_num
  character(len=*),             intent(in) :: aname,pname
  !
  !.. Optional Arguments ..
  integer,          optional, intent(in) :: num_var
  character(len=*), optional, intent(in) :: imsg
  !
  !.. Local Scalars ..
  integer :: ierr,nep
  !
continue
  !
  ! If n1 is <= 0, deallocate dt_array and return
  !
  if (n1 <= 0) then
    if (allocated(dt_array)) then
      deallocate ( dt_array , stat=ierr , errmsg=error_message )
      call alloc_error(pname,aname,2,line_num,__FILE__,ierr,error_message)
    end if
    return
  end if
  !
  iwt = iwt + 1
  message(iwt) = "Allocating " // trim(adjustl(aname)) // " and its components"
  if (present(imsg)) then
    message(iwt) = trim(adjustl(message(iwt))) // " " // trim(adjustl(imsg))
  end if
  call debug_timer(start_timer,message(iwt))
  !
  ! First allocate the derived type array
  !
  if (.not. allocated(dt_array)) then
    allocate ( dt_array(1:n1) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,aname,1,line_num,__FILE__,ierr,error_message)
  end if
  !
  ! Now allocate the components of each element of the array dt_array
  !
  call alloc_ne_dt_components(dt_array,nep,aname,pname,line_num,num_var)
  !
  call debug_timer(stop_timer,message(iwt),return_time=wall_time(iwt))
  !
end subroutine alloc_ne3d_dt
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine sort_and_output_wall_times(wt,msg)
  !
  !.. Use Statements ..
  use module_kind_types, only : iout,ten2
  !
  !.. Formal Arguments ..
  real(wp),         intent(in) :: wt(:)
  character(len=*), intent(in) :: msg(:)
  !
  !.. Local Scalars ..
  integer  :: i,n,nwt
  real(wp) :: total_time
  !
  !.. Local Arrays ..
  integer  :: sort_order(1:size(wt))
  real(wp) :: temp_time(1:size(wt))
  real(wp) :: percent_time(1:size(wt))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "sort_and_output_wall_times"
  !
continue
  !
#ifdef PROFILE_ON
  nwt = size(wt)
  !
  do n = 1,nwt
    temp_time(n) = max(zero,wt(n))
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
    write (iout,2) wt(i),percent_time(i),trim(adjustl(msg(i)))
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
subroutine alloc_face_xyz_dt_components(array,aname,pname,line_num,num_var)
  !
  !.. Use Statements ..
  use module_kind_types, only : skip_alloc_pause
  use geovar,            only : nr,face
  use order_mod,         only : geom_solpts
  !
  !.. Formal Arguments ..
  type(face_xyz), intent(inout) :: array(:)
  character(len=*),  intent(in) :: aname,pname
  integer,           intent(in) :: line_num
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
continue
  !
  nvar = nr
  if (present(num_var)) then
    nvar = num_var
  end if
  !
  do nf = lbound(array,dim=1),ubound(array,dim=1)
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    nfp = geom_solpts(face_geom,face_order)
    !
    allocate ( array(nf)%v(1:nvar,1:nfp) , source=zero , &
               stat=ierr , errmsg=error_message )
#ifndef DEBUG_ON
    write (array_name,1) trim(adjustl(aname)),nf
    call alloc_error(pname,array_name,1,line_num,__FILE__,ierr, &
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
subroutine alloc_face_var_dt_components(array,aname,pname,line_num,num_var)
  !
  !.. Use Statements ..
  use module_kind_types, only : skip_alloc_pause
  use eqn_idx,           only : nq
  use geovar,            only : face
  use order_mod,         only : geom_solpts
  !
  !.. Formal Arguments ..
  type(face_var), intent(inout) :: array(:)
  character(len=*),  intent(in) :: aname,pname
  integer,           intent(in) :: line_num
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
continue
  !
  nvar = nq
  if (present(num_var)) then
    nvar = num_var
  end if
  !
  do nf = lbound(array,dim=1),ubound(array,dim=1)
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    nfp = geom_solpts(face_geom,face_order)
    !
    allocate ( array(nf)%v(1:nvar,1:nfp,1:2) , source=zero , &
               stat=ierr , errmsg=error_message )
#ifndef DEBUG_ON
    write (array_name,1) trim(adjustl(aname)),nf
    call alloc_error(pname,array_name,1,line_num,__FILE__,ierr, &
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
subroutine alloc_face_vec_dt_components(array,aname,pname,line_num,num_var)
  !
  !.. Use Statements ..
  use module_kind_types, only : skip_alloc_pause
  use eqn_idx,           only : nq
  use geovar,            only : nr,face
  use order_mod,         only : geom_solpts
  !
  !.. Formal Arguments ..
  type(face_vec), intent(inout) :: array(:)
  character(len=*),  intent(in) :: aname,pname
  integer,           intent(in) :: line_num
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
continue
  !
  nvar = nq
  if (present(num_var)) then
    nvar = num_var
  end if
  !
  do nf = lbound(array,dim=1),ubound(array,dim=1)
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    nfp = geom_solpts(face_geom,face_order)
    !
    allocate ( array(nf)%v(1:nr,1:nvar,1:nfp,1:2) , source=zero , &
               stat=ierr , errmsg=error_message )
#ifndef DEBUG_ON
    write (array_name,1) trim(adjustl(aname)),nf
    call alloc_error(pname,array_name,1,line_num,__FILE__,ierr, &
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
subroutine alloc_node_var_dt_components(array,aname,pname,line_num,num_var)
  !
  !.. Use Statements ..
  use module_kind_types, only : skip_alloc_pause
  use eqn_idx,           only : nq
  !
  !.. Formal Arguments ..
  type(ne1d_t),  intent(inout) :: array(:)
  character(len=*), intent(in) :: aname,pname
  integer,          intent(in) :: line_num
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: num_var
  !
  !.. Local Scalars ..
  integer :: nn,nvar,ierr
  !
  character(len=200) :: array_name
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
    call alloc_error(pname,array_name,1,line_num,__FILE__,ierr, &
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
subroutine alloc_node_vec_dt_components(array,aname,pname,line_num,num_var)
  !
  !.. Use Statements ..
  use module_kind_types, only : skip_alloc_pause
  use eqn_idx,           only : nq
  use geovar,            only : nr
  !
  !.. Formal Arguments ..
  type(ne2d_t),  intent(inout) :: array(:)
  character(len=*), intent(in) :: aname,pname
  integer,          intent(in) :: line_num
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: num_var
  !
  !.. Local Scalars ..
  integer :: nn,nvar,ierr
  !
  character(len=200) :: array_name
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
    call alloc_error(pname,array_name,1,line_num,__FILE__,ierr, &
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
subroutine alloc_edge_var_dt_components(array,npts,aname,pname,line_num,num_var)
  !
  !.. Use Statements ..
  use module_kind_types, only : Geom_Edge,skip_alloc_pause
  use eqn_idx,           only : nq
  use geovar,            only : edge
  use order_mod,         only : geom_solpts
  !
  !.. Formal Arguments ..
  type(ne2d_t),  intent(inout) :: array(:)
  integer,       intent(inout) :: npts
  character(len=*), intent(in) :: aname,pname
  integer,          intent(in) :: line_num
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: num_var
  !
  !.. Local Scalars ..
  integer :: ne,nepts,nvar,ierr
  !
  character(len=200) :: array_name
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
    call alloc_error(pname,array_name,1,line_num,__FILE__,ierr, &
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
subroutine alloc_edge_vec_dt_components(array,npts,aname,pname,line_num,num_var)
  !
  !.. Use Statements ..
  use module_kind_types, only : Geom_Edge,skip_alloc_pause
  use eqn_idx,           only : nq
  use geovar,            only : nr,edge
  use order_mod,         only : geom_solpts
  !
  !.. Formal Arguments ..
  type(ne3d_t),  intent(inout) :: array(:)
  integer,       intent(inout) :: npts
  character(len=*), intent(in) :: aname,pname
  integer,          intent(in) :: line_num
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: num_var
  !
  !.. Local Scalars ..
  integer :: ne,nepts,nvar,ierr
  !
  character(len=200) :: array_name
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
    call alloc_error(pname,array_name,1,line_num,__FILE__,ierr, &
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
  use module_kind_types, only : iout
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
  if (allocated(faceusp)) then
    call array%set( storage_size(faceusp,kind=inttype) , &
                            size(faceusp,kind=inttype) )
    call total%add(array)
    write (array_name,5) "faceusp"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "faceusp"
  end if
  if (allocated(facedusp)) then
    call array%set( storage_size(facedusp,kind=inttype) , &
                            size(facedusp,kind=inttype) )
    call total%add(array)
    write (array_name,5) "facedusp"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "facedusp"
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
  if (allocated(facexyz)) then
    call array%set( storage_size(facexyz,kind=inttype) , &
                            size(facexyz,kind=inttype) )
    call total%add(array)
    write (array_name,5) "facexyz"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "facexyz"
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
