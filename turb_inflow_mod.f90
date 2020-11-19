module synthetic_eddy_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nee,nq,ntemp,npres
  !
  implicit none
  !
  private
  !
  public :: sem
  public :: separate_sem_boundary_faces
  public :: output_sem_averages
  public :: dump_sem_data
  public :: initialize_from_sem_profile
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !###                                                                       ###
  !###   MODULE DERIVED-TYPES, PARAMETERS, VARIABLES, AND ALLOCATABLE        ###
  !###         ARRAYS FOR USING THE SYNTHETIC EDDY METHOD (SEM)              ###
  !###                                                                       ###
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! SEM DERIVED-TYPES
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  type :: sem_conn_t
    integer, allocatable :: ipelem(:,:)
  end type sem_conn_t
  !
  type(sem_conn_t), allocatable :: sem_conn(:)
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  type :: eddy_list_t
    integer :: num_eddies = 0
    integer, allocatable :: eddy_list(:)
  end type eddy_list_t
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  ! COORD_T : Derived-type to store the x,y,z location for an object
  !
  type :: coord_t
    real(wp) :: x = zero
    real(wp) :: y = zero
    real(wp) :: z = zero
  end type coord_t
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  ! WALL_T : Derived-type to store the location of a
  !          wall point and distance to that point
  !
  type :: wall_t
    real(wp) :: x,y,z
    real(wp) :: xnrm = zero
    real(wp) :: ynrm = zero
    real(wp) :: znrm = zero
    real(wp) :: distance = huge(zero)
  end type wall_t
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  ! SEM_FP_T : Derived-type to store the local SEM data for each SEM face point.
  !            This data includes
  !              - the index of the face point
  !              - the x,y,z coordinates of the face point
  !              - the list of eddies affecting the face point
  !              - the values of the SEM profile data interpolated to
  !                the face point.
  !
  type :: sem_fp_t
    ! IDX : gives the face point index for this SEM point
    !          i.e.,   rho = faceusp_array(idx,nec,2)
    integer :: idx
    ! X : gives the x coordinate for this SEM point
    ! Y : gives the y coordinate for this SEM point
    ! Z : gives the z coordinate for this SEM point
    ! R : gives the radial coordinate for this SEM point
    ! T : gives the angular coordinate for this SEM point
    !     These components are to eliminate the need for storing the array
    !     facexyz_dt or facexyz_array, because those are only really needed
    !     for MMS [facexyz_array(1:nr,1:n_bnd_flxpts)] or for the continuous
    !     output feature [facexyz_array(1:nr,1:n_flxpts)]. For SEM, we only
    !     need the subset of (1:n_bnd_flxpts) that are the flux/face points
    !     for the SEM boundary conditions, so we dont want to waste memory
    !     for the flux/face points for which we will never require their xyz
    !     coordinates.
    real(wp) :: x,y,z,r,t
    ! SIN_THETA : value of sin(atan2(z,y)) for this SEM point
    real(wp) :: sin_theta = zero
    ! COS_THETA : value of cos(atan2(z,y)) for this SEM point
    real(wp) :: cos_theta = zero
    ! RUNNING_AVE : Running average of the eddy concentration
    !               for this face point
    real(wp) :: running_ave = zero
    ! NUM_EDDIES : current number of eddies in eddy_list
    integer :: num_eddies = 0
    ! EDDY_LIST : list of eddies that contribute to this SEM point
    integer, allocatable :: eddy_list(:)
    ! WALL : derived type storing the location of the nearest
    !        wall point and the distance to that point
    type(wall_t) :: wall
    ! RAVG : averaged inflow density profile
    real(wp) :: ravg
    ! UAVG : averaged inflow velocity profile, x-direction
    ! VAVG : averaged inflow velocity profile, y-direction
    ! WAVG : averaged inflow velocity profile, z-direction
    real(wp) :: uavg
    real(wp) :: vavg
    real(wp) :: wavg
    ! TAVG : averaged inflow temperature profile
    real(wp) :: tavg
    ! R11,R21,R22,R31,R32,R33 : components of the prescribed Reynolds stress
    !                           tensor profile
    real(wp) :: r11 = zero
    real(wp) :: r21 = zero
    real(wp) :: r22 = zero
    real(wp) :: r31 = zero
    real(wp) :: r32 = zero
    real(wp) :: r33 = zero
    ! A11,A21,A22,A31,A32,A33 : components of the Cholesky decomposition of
    !                           the precribed Reynolds stress tensor
    real(wp) :: a11 = zero
    real(wp) :: a21 = zero
    real(wp) :: a22 = zero
    real(wp) :: a31 = zero
    real(wp) :: a32 = zero
    real(wp) :: a33 = zero
    !
!#define ENABLE_OUTPUT_TIME_HISTORY
#ifdef ENABLE_OUTPUT_TIME_HISTORY
    real(wp), allocatable :: time_history(:,:)
#endif
  end type sem_fp_t
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  ! SEM_PTS_T : Derived-type to store information about all the boundary
  !             flux points on a partition with an SEM boundary condition.
  !             NOTE: These will be referred to as SEM face points,
  !                   SEM flux points, or just SEM points.
  !
  type :: sem_pts_t
    ! BC_IDX : index within bc_in that these SEM points correspond to
    integer :: bc_idx
    ! GRP_IDX : index within bc_grps that these SEM points correspond to
    integer :: grp_idx
    ! NUM_FP : Number of SEM points on this partition for this boundary
    !          condition group.
    integer :: num_fp = 0
    ! SEM_FP : An array of all the SEM points on this partition for this
    !          boundary condition group. Each array element contains the
    !          face point index and coodinates, as well as the inflow
    !          profile data for the SEM point.
    type(sem_fp_t), allocatable :: sem_fp(:)
  end type sem_pts_t
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  ! SEM_PROFILE_T : Derived-type to store the data read
  !                 in from an SEM profile file.
  !
  type :: sem_profile_t
    !
    character(len=100) :: fname = ""
    !
    integer :: npts = 0
    !
    real(wp), allocatable :: y(:)
    !
    real(wp), allocatable :: ravg(:)
    real(wp), allocatable :: uavg(:)
    real(wp), allocatable :: vavg(:)
    real(wp), allocatable :: wavg(:)
    real(wp), allocatable :: tavg(:)
    !
    real(wp), allocatable :: tau11(:)
    real(wp), allocatable :: tau21(:)
    real(wp), allocatable :: tau22(:)
    real(wp), allocatable :: tau31(:)
    real(wp), allocatable :: tau32(:)
    real(wp), allocatable :: tau33(:)
    !
    real(wp), allocatable :: a11(:)
    real(wp), allocatable :: a21(:)
    real(wp), allocatable :: a22(:)
    real(wp), allocatable :: a31(:)
    real(wp), allocatable :: a32(:)
    real(wp), allocatable :: a33(:)
    !
  end type sem_profile_t
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  ! EDDY_T : Derived-type to store the properties each SEM eddy.
  !          These properties include:
  !            - x,y,z,r,theta coordinates of the eddy
  !            - the eddy intensity (+1 or -1) for each x,y,z direction
  !            - the eddy length scale (1/frequency) for each x,y,z direction
  !            - the effective eddy volume
  !            - the list of face points that are affected by the eddy
  !
  type, public :: eddy_t
    ! X : gives the x coordinate for this eddy
    ! Y : gives the y coordinate for this eddy
    ! Z : gives the z coordinate for this eddy
    ! R : gives the radial coordinate for this eddy
    ! T : gives the azimuthal (theta) coordinate for this eddy
    real(wp) :: x = zero
    real(wp) :: y = zero
    real(wp) :: z = zero
    real(wp) :: r = zero
    real(wp) :: t = zero
    ! X_INTENSITY : eddy intensity in the x direction
    ! Y_INTENSITY : eddy intensity in the y direction
    ! Z_INTENSITY : eddy intensity in the z direction
    real(wp) :: x_intensity = zero
    real(wp) :: y_intensity = zero
    real(wp) :: z_intensity = zero
    ! X_FREQUENCY : eddy turbulence frequency in the x direction
    ! Y_FREQUENCY : eddy turbulence frequency in the y direction
    ! Z_FREQUENCY : eddy turbulence frequency in the z direction
    !               NOTE: These are really just the eddy turbulence length
    !                     scales in each direction. We are just storing these
    !                     quantities as turbulence frequencies
    !                       (turb. freq.) = 1.0 / (turb. length scale)
    !                     because the length scales are only ever appear in the
    !                     denominators of equations in which they are used, and
    !                     pre-computing the reciprocals of the length scales
    !                     will help with computational efficiencies.
    real(wp) :: x_frequency = zero
    real(wp) :: y_frequency = zero
    real(wp) :: z_frequency = zero
    ! VOLUME : eddy effective volume
    !          NOTE: The eddy volume is stored as [1/(eddy volume)] because
    !                that is the only way it appears in the equations in which
    !                it is used.
    real(wp) :: volume = zero
    ! NUM_FP : number of face points affected by this eddy
    integer :: num_fp = 0
    ! FP_LIST : list of face points affected by this eddy
    integer, allocatable :: fp_list(:)
  end type eddy_t
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  ! SEM_BOX_T : Derived-type to store information about all the eddies for
  !             each SEM box.
  !
  type, public :: sem_box_t
    ! NUM_EDDIES : number of eddies in this SEM box
    integer :: num_eddies = 0
    ! BC_IDX : index within bc_in that corresponds to this SEM box
    integer :: bc_idx = 0
    ! GRP_IDX : index within bc_grps that corresponds to this SEM box
    integer :: grp_idx = 0
    ! LENGTH_SCALE : turbulence length scale for this SEM box
    real(wp) :: length_scale
    ! CONVECTIVE_VELOCITY : velocity used to convect eddies in this SEM box
    real(wp) :: convective_velocity
    ! BOXMIN : coordinate minimums for this SEM box
    real(wp) :: boxmin(1:5) ! [x,y,z,r,theta]
    ! BOXMAX : coordinate maximums for this SEM box
    real(wp) :: boxmax(1:5)
    ! BOX_VOLUME : volume of this SEM box
    real(wp) :: box_volume
    ! ALPHA_AOA : alpha angle of attack for the SEM inflow boundary
    real(wp) :: alpha_aoa = zero
    ! BETA_AOA : beta angle of attack for the SEM inflow boundary
    real(wp) :: beta_aoa = zero
    ! IS_CYLINDRICAL : logical flag indicating if the box extents are
    !                  best defined in cylindrical coordinates
    logical(lk) :: is_cylindrical = fals
    ! CENTER : coordinates for the center of this SEM box
    type(coord_t) :: center
    ! RNRM : collective inward facing normal for the inflow plane associated
    !        with this SEM box. This is needed for the coordinate transform
    !        required when the inflow plane is now perpendicular to x.
    type(coord_t) :: rnrm
    ! PROFILE : inflow profile related to this SEM box
    type(sem_profile_t), allocatable :: profile
    ! EDDY : information about each eddy in this SEM box
    type(eddy_t), allocatable :: eddy(:)
    !
  end type sem_box_t
  !
  !#############################################################################
  !
  type :: sem_fp_avg_t
    real(wp) :: r,ru,rv,rw,e
    real(wp) :: t,u,v,w,p
    real(wp) :: uu,vv,ww,uv,uw,vw
    real(wp) :: ut,vt,wt,rt
    real(wp) :: ruu,rvv,rww,ruv,ruw,rvw
    real(wp) :: rut,rvt,rwt
  end type sem_fp_avg_t
  !
  type :: sem_avg_t
    type(sem_fp_avg_t), allocatable :: fp_avg(:)
  end type sem_avg_t
  !
  type(sem_avg_t), allocatable :: sem_avg(:)
  !
  real(wp), save :: sem_prev_ave_time = -one
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! SEM PARAMETERS
  !
  integer, parameter :: min_eddies = 5000
  integer, parameter :: max_eddies = 500000
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! SEM VARIABLES
  !
  ! NUM_SEM_BOXES : number of SEM boxes
  !
  integer, save :: num_sem_boxes = 0
  !
  ! SEM_IS_INITIALIZED : logical flag that indicates whether or not the SEM
  !                      module has been initialized.
  !
  logical(lk), save :: sem_is_initialized = fals
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! SEM ALLOCATABLE ARRAYS
  !
  ! SEM_BOX : Derived-type array to store information about all eddies for all
  !           SEM boxes. Aside from the eddy face-point connectivity lists:
  !               sem_box(:)%eddy(:)%num_fp, and
  !               sem_box(:)%eddy(:)%fp_list,
  !           the sem_box(:) array should contain exactly the same information
  !           across all partitions/MPI-processes.
  !           NOTE: This needs to be public so that the eddy information can
  !                 be read from a restart file.
  !
  type(sem_box_t), public, save, allocatable :: sem_box(:)
  !
  ! MY_SEM_PTS : Derived-type array that contains the SEM information
  !              for all face points of the current partition that
  !              are on an SEM boundary. The size of this array will
  !              be the total number of SEM boundaries (1:num_sem_boxes)
  !              which is the same size as the array sem_box(:). For
  !              a given SEM box (given by index 'sem_idx'), if a
  !              partition does not have any boundary face points
  !              belonging to this SEM boundary/box, then:
  !                  my_sem_pts(sem_idx)%num_fp = 0, and
  !                  my_sem_pts(sem_idx)%sem_fp(:) is either unallocated or
  !                     is allocated with zero size [array bounds of (1:0)]
  !
  type(sem_pts_t), public, save, allocatable :: my_sem_pts(:)
  !
  ! XYZ_WALL : Array for storing the coordinates for all boundary
  !            flux points on physical wall boundaries and the
  !            interior facing unit normal at each point.
  !               xyz_wall(1,ifp) : x coordinate for flux point ifp
  !               xyz_wall(2,ifp) : y coordinate for flux point ifp
  !               xyz_wall(3,ifp) : z coordinate for flux point ifp
  !               xyz_wall(4,ifp) : x-component of the interior facing
  !                                 unit normal at flux point ifp
  !               xyz_wall(5,ifp) : y-component of the interior facing
  !                                 unit normal at flux point ifp
  !               xyz_wall(6,ifp) : z-component of the interior facing
  !                                 unit normal at flux point ifp
  !
 !real(wp), allocatable :: xyz_wall(:,:)
  !
  ! XYZ_SEM : Array for storing the coordinates for all boundary flux points
  !           on SEM inflow boundaries.
  !               xyz_sem(1,ifp) : x coordinate for flux point ifp
  !               xyz_sem(2,ifp) : y coordinate for flux point ifp
  !               xyz_sem(3,ifp) : z coordinate for flux point ifp
  !
 !real(wp), allocatable :: xyz_sem(:,:)
  !
contains
!
!###############################################################################
!
subroutine initialize_from_sem_profile(xyz,cv)
  !
  !.. Use Statements ..
  use ovar, only : bc_in,bc_input
  use ovar, only : rhoref,aref,tref
  use ovar, only : itestcase
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: xyz(:,:)
  real(wp), intent(inout) :: cv(:,:)
  !
  !.. Local Scalars ..
  integer :: k,n,npts,idx,iprof,nprof,bc_idx,ierr
  !
  !.. Local Allocatable Arrays ..
  type(sem_profile_t), allocatable :: prof(:)
  real(wp),            allocatable :: pvavg(:,:)
  integer,             allocatable :: bc_idx_prof(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_from_sem_profile"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate prof
  !
  nprof = count(bc_in(:)%bc_type == bc_sem)
  !
  allocate ( prof(1:nprof) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"prof",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate pvavg and initialize it to reference conditions
  !
  pvavg = init_array(nq,nprof,zero) ! F2003 AUTO-REALLOC
  !
  pvavg(nec,:) = bc_in(0)%cv(nec)
  pvavg(nmx,:) = bc_in(0)%pv(nmx)
  pvavg(nmy,:) = bc_in(0)%pv(nmy)
  pvavg(nmz,:) = bc_in(0)%pv(nmz)
  pvavg(nee,:) = bc_in(0)%gp(ntemp)
  !
  ! Allocate bc_idx_prof
  !
  bc_idx_prof = init_array(nprof,0) ! F2003 AUTO-REALLOC
  !
  ! Read in the profile file
  ! NOTE: This should be called before initialize_sem, so sem_box will not yet
  !       be allocated, and therefore the profile file will not yet have been
  !       read in and saved to sem_box(:)%profile.
  !
  iprof = 0
  one_pass_loop: &
  do k = 1,1
    !
    do bc_idx = 1,ubound(bc_in,dim=1)
      if (bc_in(bc_idx)%bc_type == bc_sem) then
        iprof = iprof + 1
        bc_idx_prof(iprof) = bc_idx
        prof(iprof)%fname = trim(adjustl( bc_in(bc_idx)%sem%profile_file ))
        if (len_trim(prof(iprof)%fname) > 0) then
          pvavg(nec,iprof) = bc_in(bc_idx)%cv(nec)
          pvavg(nmx,iprof) = bc_in(bc_idx)%pv(nmx)
          pvavg(nmy,iprof) = bc_in(bc_idx)%pv(nmy)
          pvavg(nmz,iprof) = bc_in(bc_idx)%pv(nmz)
          pvavg(nee,iprof) = bc_in(bc_idx)%gp(ntemp)
        end if
        if (iprof == nprof) exit one_pass_loop
      end if
    end do
    !
    write (error_message,1)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    !
  end do one_pass_loop
  !
  ! Read in the SEM profiles
  !
  do iprof = 1,nprof
    !
    call read_sem_profile(prof(iprof))
    !
    npts = prof(iprof)%npts
    !
    ! If the averaged density was not included in the profile
    ! file, set the density for all points in the profile to
    ! the density specified for this boundary condition.
    !
    if (.not. allocated(prof(iprof)%ravg)) then
      prof(iprof)%ravg = init_array(npts,pvavg(nec,iprof)) ! F2003 AUTO-REALLOC
    end if
    !
    ! If the averaged x-velocity was not included in the profile
    ! file, set the x-velocity for all points in the profile to
    ! the x-velocity specified for this boundary condition.
    !
    if (.not. allocated(prof(iprof)%uavg)) then
      prof(iprof)%uavg = init_array(npts,pvavg(nmx,iprof)) ! F2003 AUTO-REALLOC
    end if
    !
    ! If the averaged y-velocity was not included in the profile
    ! file, set the y-velocity for all points in the profile to
    ! the y-velocity specified for this boundary condition.
    !
    if (.not. allocated(prof(iprof)%vavg)) then
      prof(iprof)%vavg = init_array(npts,pvavg(nmy,iprof)) ! F2003 AUTO-REALLOC
    end if
    !
    ! If the averaged z-velocity was not included in the profile
    ! file, set the z-velocity for all points in the profile to
    ! the z-velocity specified for this boundary condition.
    !
    if (.not. allocated(prof(iprof)%wavg)) then
      prof(iprof)%wavg = init_array(npts,pvavg(nmz,iprof)) ! F2003 AUTO-REALLOC
    end if
    !
    ! If the averaged temperature was not included in the profile
    ! file, set the temperature for all points in the profile to
    ! the temperature specified for this boundary condition.
    !
    if (.not. allocated(prof(iprof)%tavg)) then
      prof(iprof)%tavg = init_array(npts,pvavg(nee,iprof)) ! F2003 AUTO-REALLOC
    end if
    !
  end do
  !
  if (mypnum == glb_root) then
    write (iout,*)
    do iprof = 1,nprof
      write (iout,2) bc_idx_prof(iprof),iprof, &
                     pvavg(nec,iprof)*rhoref, &
                     pvavg(nmx,iprof)*aref, &
                     pvavg(nmy,iprof)*aref, &
                     pvavg(nmz,iprof)*aref, &
                     pvavg(nee,iprof)*tref
      npts = prof(iprof)%npts
      write (iout,2) bc_idx_prof(iprof),iprof, &
                     prof(iprof)%ravg(npts)*rhoref, &
                     prof(iprof)%uavg(npts)*aref, &
                     prof(iprof)%vavg(npts)*aref, &
                     prof(iprof)%wavg(npts)*aref, &
                     prof(iprof)%tavg(npts)*tref
    end do
    write (iout,*)
  end if
  flush (iout)
  !
  ! Interpolate the flow profile to the solution points
  !
  if (itestcase == -THX_Phase_3) then
    !
    iprof = 1
    !
    find_pipe_sem: &
    do n = 1,nprof
      idx = bc_in(bc_idx_prof(iprof))%bc_input_idx
      if (trim(adjustl(uppercase(bc_input(idx)%name))) == "INFLOWPIPE") then
        iprof = n
        exit find_pipe_sem
      end if
    end do find_pipe_sem
    !
    call thx_pipe_sem_initialization(prof(iprof),pvavg(:,iprof),xyz,cv)
    !
  else if (itestcase == THX_Phase_3) then
    !
    if (size(bc_idx_prof) >= 2) then
      call thx_sem_initialization(prof,pvavg,bc_idx_prof,xyz,cv)
    else
      call thx_sem_plenum_initialization(prof,pvavg,bc_idx_prof,xyz,cv)
    end if
    !
  else
    !
    iprof = 1
    !
    call init_from_sem_profile(prof(iprof),pvavg(:,iprof),xyz,cv)
    !
  end if
  !
  !
  ! Deallocate all local arrays before leaving
  !
  if (allocated(prof)) then
    deallocate ( prof , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"prof",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(pvavg)) then
    deallocate ( pvavg , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"pvavg",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(bc_idx_prof)) then
    deallocate ( bc_idx_prof , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bc_idx_prof",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("ERROR!! TRYING TO INITIALIZE CHANNEL FLOW CASE USING ", &
            "THE SEM INFLOW PROFILE, BUT A FILE NAME FOR THE SEM ", &
            "PROFILE FILE COULD NOT BE FOUND FROM THE INPUT BOUNDARY ", &
            "CONDITIONS!")
  !
  2 format (/,2x,"Primitive variables from the SEM BC group # ",i0,/, &
              2x,"that will be used for initialization group # ",i0," : ",/, &
              6x,"ravg = ",f18.10," kg/(m^3)",/, &
              6x,"uavg = ",f18.10," m/s",/, &
              6x,"vavg = ",f18.10," m/s",/, &
              6x,"wavg = ",f18.10," m/s",/, &
              6x,"tavg = ",f18.10," K",/)
  !
end subroutine initialize_from_sem_profile
!
!###############################################################################
!
subroutine thx_sem_plenum_initialization(prof,pvavg,bc_idx_prof,xyz,cv)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
  use ovar, only : grid_scaling_factor
  use ovar, only : bc_in,bc_input
  !
  !.. Formal Arguments ..
  type(sem_profile_t), intent(in) :: prof(:)
  real(wp),            intent(in) :: pvavg(:,:)
  integer,             intent(in) :: bc_idx_prof(:)
  real(wp),            intent(in) :: xyz(:,:)
  real(wp),         intent(inout) :: cv(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,n1,n2,npts,idx,iprof
  integer  :: nozzle_prof,mywall
  real(wp) :: c1,c2,y,y1,y2,rlx_vel
  real(wp) :: r,u,v,w,t,q,q2,uu,vv,ww
  real(wp) :: rgas_over_gm1,ujet,uinf
  real(wp) :: xp,yp,zp
  real(wp) :: cos_theta,sin_theta
  !
  !.. Local Allocatable Arrays ..
  character(len=:), allocatable :: prof_name
  integer,          allocatable :: init_grp(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "thx_sem_plenum_initialization"
  !
  real(wp), parameter :: height = 2.68_wp
  real(wp), parameter :: half_height = half*height
  !
  logical(lk), parameter :: using_constant_area_nozzle = true
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nozzle_prof = 0
  !
  do iprof = 1,size(bc_idx_prof)
    !
    idx = bc_in(bc_idx_prof(iprof))%bc_input_idx
    !
    prof_name = trim(adjustl(uppercase(bc_input(idx)%name))) ! F03 AUTO-REALLOC
    !
    if (prof_name == "INFLOWNOZZLE") then
      nozzle_prof = iprof
    end if
    !
  end do
  !
  if (nozzle_prof == 0) then
    error_message = "Unable to find SEM profile for the nozzle!"
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Start by initializing everything to freestream
  !
  do k = 1,size(cv,dim=1)
    cv(k,nec) = bc_in(0)%cv(nec)
    cv(k,nmx) = bc_in(0)%cv(nmx)
    cv(k,nmy) = bc_in(0)%cv(nmy)
    cv(k,nmz) = bc_in(0)%cv(nmz)
    cv(k,nee) = bc_in(0)%cv(nee)
  end do
  !
  ! Allocate the array for identifying which section
  ! of the grid each solution point is in.
  !
  init_grp = init_array(size(xyz,dim=2),0) ! F2003 AUTO-REALLOC
  !
  ! Find all the points that will be initialized to the nozzle SEM profile
  !
  do k = 1,size(xyz,dim=2)
    !
    yp = xyz(2,k)/grid_scaling_factor
    zp = xyz(3,k)/grid_scaling_factor
    !
    if (yp >= -eps6 .and. yp <= (height+eps4)) then
      if (zp > (-half_height-eps4) .and. zp < (half_height+eps4)) then
        init_grp(k) = nozzle_prof
      end if
    end if
    !
  end do
  !
  ! Precompute Rgas/(gamma-1)
  !
  rgas_over_gm1 = rgasref_nd / (gam - one)
  !
  uinf = bc_in(0)%pv(nmx)
  ujet = prof(nozzle_prof)%uavg( prof(nozzle_prof)%npts )
  !
  init_loop: &
  do k = 1,size(xyz,dim=2)
    !
    iprof = init_grp(k)
    !
    ! If iprof is zero, this point should be initialized to free stream
    ! conditions which was already done at the beginning of this subroutine.
    !
    if (iprof == 0) cycle init_loop
    !
    ! Determine the wall distance based on which
    ! section of the grid this point is in
    !
    ! Get the location of this solution point in native grid units
    !
    xp = xyz(1,k)/grid_scaling_factor
    yp = xyz(2,k)/grid_scaling_factor
    zp = xyz(3,k)/grid_scaling_factor
    !
    ! Determine where this point is axially 
    !
    if (xp >= twelve) then
      !
      ! This point is beyond the plate
      !
      y = xyz(2,k)
      !
      cos_theta = one
      sin_theta = zero
      !
      ! This sets the relaxation factor between the values 1 and uinf/ujet,
      ! with xp=12 -> rlx_vel=1 and xp=48 -> rlx_vel=uinf/ujet.
      !
      rlx_vel = ( uinf*(xp-twelve) - ujet*(xp-four*twelve) ) / &
                (ujet*three*twelve)
      !
    else if (xp >= zero) then
      !
      ! This point is above the plate
      !
      y = xyz(2,k)
      !
      cos_theta = one
      sin_theta = zero
      !
      rlx_vel = one
      !
    else
      !
      ! This point is in the nozzle
      !
      ! Find which wall of the nozzle the current point is closest to
      !
      mywall = minloc(abs([yp,half_height-zp,height-yp,zp+half_height]),dim=1)
      !
      if (mywall == 1) then
        ! Closest to bottom (y-min) wall
        y = abs(xyz(2,k))
        cos_theta = one
        sin_theta = zero
      else if (mywall == 2) then
        ! Closest to right (z-max) wall
        y = abs(half_height - zp)*grid_scaling_factor
        cos_theta = zero
        sin_theta = one
      else if (mywall == 3) then
        ! Closest to top (y-max) wall
        y = abs(height - yp)*grid_scaling_factor
        cos_theta = -one
        sin_theta = zero
      else if (mywall == 4) then
        ! Closest to left (z-min) wall
        y = abs(zp + half_height)*grid_scaling_factor
        cos_theta = zero
        sin_theta = -one
      end if
      !
      rlx_vel = one
      !
    end if
    !
    ! Interpolate the primitive variables from the SEM profile
    ! based on the wall distance of the current point
    !
    npts = prof(iprof)%npts
    !
    n1 = 1
    n2 = 2
    search_profile_loop: &
    do n = 2,npts
      if (prof(iprof)%y(n) > y) exit search_profile_loop
      n1 = n
      n2 = max(n+1,npts)
    end do search_profile_loop
    !
    if (n1 == npts) then
      !
      r  = prof(iprof)%ravg(npts)
      uu = prof(iprof)%uavg(npts)
      vv = prof(iprof)%vavg(npts)
      ww = prof(iprof)%wavg(npts)
      t  = prof(iprof)%tavg(npts)
      !
    else
      !
      y1 = prof(iprof)%y(n1)
      y2 = prof(iprof)%y(n2)
      !
      c1 = (y-y2)/(y2-y1)
      c2 = (y-y1)/(y2-y1)
      !
      r  = c2*prof(iprof)%ravg(n2) - c1*prof(iprof)%ravg(n1)
      uu = c2*prof(iprof)%uavg(n2) - c1*prof(iprof)%uavg(n1)
      vv = c2*prof(iprof)%vavg(n2) - c1*prof(iprof)%vavg(n1)
      ww = c2*prof(iprof)%wavg(n2) - c1*prof(iprof)%wavg(n1)
      t  = c2*prof(iprof)%tavg(n2) - c1*prof(iprof)%tavg(n1)
      !
    end if
    !
    u = uu
    v =  vv*cos_theta + ww*sin_theta
    w = -vv*sin_theta + ww*cos_theta
    !
    u = u * rlx_vel
    v = v * rlx_vel
    w = w * rlx_vel
    !
    q2 = u*u + v*v + w*w
    !
    cv(k,nec) = r
    cv(k,nmx) = r*u
    cv(k,nmy) = r*v
    cv(k,nmz) = r*w
    cv(k,nee) = r*( rgas_over_gm1*t + half*q2 )
    !
  end do init_loop
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine thx_sem_plenum_initialization
!
!###############################################################################
!
subroutine thx_sem_initialization(prof,pvavg,bc_idx_prof,xyz,cv)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
  use ovar, only : grid_scaling_factor
  use ovar, only : bc_in,bc_input
  !
  !.. Formal Arguments ..
  type(sem_profile_t), intent(in) :: prof(:)
  real(wp),            intent(in) :: pvavg(:,:)
  integer,             intent(in) :: bc_idx_prof(:)
  real(wp),            intent(in) :: xyz(:,:)
  real(wp),         intent(inout) :: cv(:,:)
  !
  !.. Local Scalars ..
  integer  :: k,n,n1,n2,npts,idx,iprof
  integer  :: pipe_prof,nozzle_prof
  real(wp) :: c1,c2,y,y1,y2
  real(wp) :: r,u,v,w,t,q,q2,uu,vv,ww
  real(wp) :: rgas_over_gm1
  real(wp) :: dx,dy,dz,dl,unx,uny,unz
  real(wp) :: xp,yp,zp,xc,yc,angle
  real(wp) :: cn,cni,v11,v12,v13,v21,v22,v23,v31,v32,v33
  !
  !.. Local Allocatable Arrays ..
  character(len=:), allocatable :: prof_name
  integer,          allocatable :: init_grp(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "thx_sem_initialization"
  !
  real(wp), parameter :: height = 2.68_wp
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  pipe_prof = 0
  nozzle_prof = 0
  !
  do iprof = 1,size(bc_idx_prof)
    !
    idx = bc_in(bc_idx_prof(iprof))%bc_input_idx
    !
    prof_name = trim(adjustl(uppercase(bc_input(idx)%name))) ! F03 AUTO-REALLOC
    !
    if (prof_name == "INFLOWPIPE") then
      pipe_prof = iprof
    else if (prof_name == "INFLOWNOZZLE") then
      nozzle_prof = iprof
    end if
    !
  end do
  !
  if (pipe_prof == 0) then
    error_message = "Unable to find SEM profile for the pipe/injector!"
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  if (nozzle_prof == 0) then
    error_message = "Unable to find SEM profile for the nozzle!"
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Allocate the array for identifying which section
  ! of the grid each solution point is in.
  !
  init_grp = init_array(size(xyz,dim=2),0) ! F2003 AUTO-REALLOC
  !
  do k = 1,size(xyz,dim=2)
    !
    xp = xyz(1,k)/grid_scaling_factor
    yp = xyz(2,k)/grid_scaling_factor
    zp = xyz(3,k)/grid_scaling_factor
    !
    if ( (xp <= zero .and. yp < -half) .or. &
         (xp > zero .and. yp < zero) ) then
      !
      ! Find all solution points that are definitely
      ! in the pipe/injector section
      !
      init_grp(k) = pipe_prof
      !
    else if (xp > zero .and. xp <= twelve .and. yp >= zero) then
   !else if (xp > zero .and. xp <= twelve .and. &
   !         yp >= zero .and. yp < half) then
      !
      ! Find all solution points above the plate
      !
      init_grp(k) = nozzle_prof
      !
    else if (xp <= zero .and. yp > -half) then
      !
      ! Find all solution points in the nozzle area
      !
      if (xp >= -one .and. xp <= zero) then
        if ( (yp >= zero .and. yp <= height) .and. &
             (zp >= -half*height .and. zp <= half*height) ) then
          init_grp(k) = nozzle_prof
        end if
      else if (xp >= -five/four .and. xp < -one) then
        if ( (yp > -one4 .and. yp < (two+three4)) .and. &
             (zp > -seven/five .and. zp < seven/five) ) then
          init_grp(k) = nozzle_prof
        end if
      else if (xp < -five/four) then
        if ( (yp > -half .and. yp < three) .and. &
             (zp > -seven/four .and. zp < seven/four) ) then
          init_grp(k) = nozzle_prof
        end if
      end if
      !
      init_grp(k) = nozzle_prof
      !
    else
      !
      ! For now, just have everything else default to the plate section
      !
      init_grp(k) = nozzle_prof
      !
    end if
    !
  end do
  !
  ! Precompute Rgas/(gamma-1)
  !
  rgas_over_gm1 = rgasref_nd / (gam - one)
  !
  init_loop: &
  do k = 1,size(xyz,dim=2)
    !
    iprof = init_grp(k)
    !
    ! Determine the wall distance based on which
    ! section of the grid this point is in
    !
    if (iprof == pipe_prof) then
      !
      ! This point is in the pipe section
      !
      xp = xyz(1,k)/grid_scaling_factor
      yp = xyz(2,k)/grid_scaling_factor
      zp = xyz(3,k)/grid_scaling_factor
      !
      ! Get the turn angle of the pipe
      !
      angle = thx_pipe_turn_angle(xp,yp)
      !
      ! Get the x and y coordinates of the point on the pipe centerline
      ! that is perpendicular the current solution point
      !
      yc = thx_pipe_centerline_y(angle,xp,yp)
      xc = thx_pipe_centerline_x(angle,yc,xp)
      !
      ! Get the distance from the current solution point to the centerline point
      ! NOTE: If (angle <= zero), then (xp = xc), so explicitly set dx to zero
      !
      dx = merge( xc - xp , zero , angle > zero )
      dy = yc - yp
      dz = -zp     ! NOTE: zc = 0
      dl = sqrt( dx*dx + dy*dy + dz*dz )
      !
      ! Get the inward facing unit normal at the nearest wall point
      !
      unx = dx/dl
      uny = dy/dl
      unz = dz/dl
      !
      ! Get the wall distance
      ! NOTE: This is just [the radius of the pipe (0.25 in)] - [dl]
      ! NOTE: Need to use grid scaling factor to convert to meters
      !
      y = (one4 - dl)*grid_scaling_factor
      !
      ! Create the coefficients for the vector transformation
      !
      cn = sqrt(uny*uny + unz*unz)
      cni = one/cn
      !
      v11 = cn           ; v12 = unx ; v13 =  zero
      v21 = -unx*uny*cni ; v22 = uny ; v23 = -unz*cni
      v31 = -unx*unz*cni ; v32 = unz ; v33 =  uny*cni
      !
    else
      !
      ! This point is in the nozzle or above the plate
      !
      y = xyz(2,k)
      !
    end if
    !
    ! Interpolate the primitive variables from the SEM profile
    ! based on the wall distance of the current point
    !
    npts = prof(iprof)%npts
    !
    n1 = 1
    n2 = 2
    search_profile_loop: &
    do n = 2,npts
      if (prof(iprof)%y(n) > y) exit search_profile_loop
      n1 = n
      n2 = max(n+1,npts)
    end do search_profile_loop
    !
    if (n1 == npts) then
      !
     !r  = pvavg(nec,iprof)
     !uu = pvavg(nmx,iprof)
     !vv = pvavg(nmy,iprof)
     !ww = pvavg(nmz,iprof)
     !t  = pvavg(nee,iprof)
      r  = prof(iprof)%ravg(npts)
      uu = prof(iprof)%uavg(npts)
      vv = prof(iprof)%vavg(npts)
      ww = prof(iprof)%wavg(npts)
      t  = prof(iprof)%tavg(npts)
      !
    else
      !
      y1 = prof(iprof)%y(n1)
      y2 = prof(iprof)%y(n2)
      !
      c1 = (y-y2)/(y2-y1)
      c2 = (y-y1)/(y2-y1)
      !
      r  = c2*prof(iprof)%ravg(n2) - c1*prof(iprof)%ravg(n1)
      uu = c2*prof(iprof)%uavg(n2) - c1*prof(iprof)%uavg(n1)
      vv = c2*prof(iprof)%vavg(n2) - c1*prof(iprof)%vavg(n1)
      ww = c2*prof(iprof)%wavg(n2) - c1*prof(iprof)%wavg(n1)
      t  = c2*prof(iprof)%tavg(n2) - c1*prof(iprof)%tavg(n1)
      !
    end if
    !
    if (iprof == pipe_prof) then
      !
      if (angle > zero) then
        q = sqrt(uu*uu + vv*vv + ww*ww)
        if (angle >= 30.0_wp) then
          u = q*half*sqrt(three) ! q*cos(30)
          v = q*half             ! q*sin(30)
        else
          u = q*cos(angle*degrees2radians)
          v = q*sin(angle*degrees2radians)
        end if
        w = zero
      else
        u = uu
        v = vv
        w = ww
      end if
      !
    else
      !
      u = uu
      v = vv
      w = ww
      !
    end if
    !
    q2 = u*u + v*v + w*w
    !
    cv(k,nec) = r
    cv(k,nmx) = r*u
    cv(k,nmy) = r*v
    cv(k,nmz) = r*w
    cv(k,nee) = r*( rgas_over_gm1*t + half*q2 )
    !
  end do init_loop
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine thx_sem_initialization
!
!###############################################################################
!
subroutine thx_pipe_sem_initialization(prof,pvavg,xyz,cv)
  !
  !.. Use Statements ..
  use ovar, only : rhoref,aref,tref
  use ovar, only : gam,rgasref_nd
  use ovar, only : grid_scaling_factor
  !
  !.. Formal Arguments ..
  type(sem_profile_t), intent(in) :: prof
  real(wp),            intent(in) :: pvavg(:)
  real(wp),            intent(in) :: xyz(:,:)
  real(wp),         intent(inout) :: cv(:,:)
  !
  !.. Local Scalars ..
  integer  :: k,n,n1,n2,npts
  real(wp) :: c1,c2,cn,cni,y,y1,y2
  real(wp) :: r,u,v,w,t,q,q2,uu,vv,ww
  real(wp) :: rgas_over_gm1
  real(wp) :: dx,dy,dz,dl,unx,uny,unz
  real(wp) :: xp,yp,zp,xc,yc,angle
  real(wp) :: v11,v12,v13,v21,v22,v23,v31,v32,v33
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "thx_pipe_sem_initialization"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  npts = prof%npts
  !
  ! Precompute Rgas/(gamma-1)
  !
  rgas_over_gm1 = rgasref_nd / (gam - one)
  !
  do k = 1,size(xyz,dim=2)
    !
    xp = xyz(1,k)/grid_scaling_factor
    yp = xyz(2,k)/grid_scaling_factor
    zp = xyz(3,k)/grid_scaling_factor
    !
    ! Get the turn angle of the pipe
    !
    angle = thx_pipe_turn_angle(xp,yp)
    !
    ! Get the x and y coordinates of the point on the pipe centerline
    ! that is perpendicular the current solution point
    !
    yc = thx_pipe_centerline_y(angle,xp,yp)
    xc = thx_pipe_centerline_x(angle,yc,xp)
    !
    ! Get the distance from the current solution point to the centerline point
    ! NOTE: If (angle <= zero), then (xp = xc), so explicitly set dx to zero
    !
    dx = merge( xc - xp , zero , angle > zero )
    dy = yc - yp
    dz = -zp     ! NOTE: zc = 0
    dl = sqrt( dx*dx + dy*dy + dz*dz )
    !
    ! Get the inward facing unit normal at the nearest wall point
    !
    unx = dx/dl
    uny = dy/dl
    unz = dz/dl
    !
    ! Get the wall distance
    ! NOTE: This is just [the radius of the pipe (0.25 in)] - [dl]
    ! NOTE: Need to use grid scaling factor to convert to meters
    !
    y = (one4 - dl)*grid_scaling_factor
    !
    ! Create the coefficients for the vector transformation
    !
    cn = sqrt(uny*uny + unz*unz)
    cni = one/cn
    !
    v11 = cn           ; v12 = unx ; v13 =  zero
    v21 = -unx*uny*cni ; v22 = uny ; v23 = -unz*cni
    v31 = -unx*unz*cni ; v32 = unz ; v33 =  uny*cni
    !
    n1 = 1
    n2 = 2
    do n = 2,npts
      if (prof%y(n) > y) exit
      n1 = n
      n2 = max(n+1,npts)
    end do
    !
    if (n1 == npts) then
      !
     !r  = pvavg(nec)
     !uu = pvavg(nmx)
     !vv = pvavg(nmy)
     !ww = pvavg(nmz)
     !t  = pvavg(nee)
      r  = prof%ravg(npts)
      uu = prof%uavg(npts)
      vv = prof%vavg(npts)
      ww = prof%wavg(npts)
      t  = prof%tavg(npts)
      !
    else
      !
      y1 = prof%y(n1)
      y2 = prof%y(n2)
      !
      c1 = (y-y2)/(y2-y1)
      c2 = (y-y1)/(y2-y1)
      !
      r  = c2*prof%ravg(n2) - c1*prof%ravg(n1)
      uu = c2*prof%uavg(n2) - c1*prof%uavg(n1)
      vv = c2*prof%vavg(n2) - c1*prof%vavg(n1)
      ww = c2*prof%wavg(n2) - c1*prof%wavg(n1)
      t  = c2*prof%tavg(n2) - c1*prof%tavg(n1)
      !
    end if
    !
    ! Transform the SEM velocity vector to align it with the
    ! new coodinate system defined by the wall unit normal
    !
#define USE_SIMPLE_ROTATION_FOR_THX_PIPE
#ifdef USE_SIMPLE_ROTATION_FOR_THX_PIPE
    if (angle > zero) then
      q = sqrt(uu*uu + vv*vv + ww*ww)
      if (angle >= 30.0_wp) then
        u = q*half*sqrt(three) ! q*cos(30)
        v = q*half             ! q*sin(30)
      else
        u = q*cos(angle*degrees2radians)
        v = q*sin(angle*degrees2radians)
      end if
      w = zero
    else
      u = uu
      v = vv
      w = ww
    end if
#else
    if (angle >= 30.0_wp) then
      q = sqrt(uu*uu + vv*vv + ww*ww)
      u = q*sqrt(three)*half
      v = half*q
      w = zero
    else
      u = uu*v11 + vv*v12 + ww*v13
      v = uu*v21 + vv*v22 + ww*v23
      w = uu*v31 + vv*v32 + ww*v33
    end if
#endif
    !
    q2 = u*u + v*v + w*w
    !
    cv(k,nec) = r
    cv(k,nmx) = r*u
    cv(k,nmy) = r*v
    cv(k,nmz) = r*w
    cv(k,nee) = r*( rgas_over_gm1*t + half*q2 )
    !
!#define DEBUG_THX_PIPE_INIT
#ifdef DEBUG_THX_PIPE_INIT
    if (sqrt(q2)*aref > 215.0_wp) then
      write (10000+mypnum,1) k,xp,yp,zp, &
                             angle,xc,yc, &
                             unx,uny,unz, &
                             dl,y,y/grid_scaling_factor
      write (20000+mypnum,2) k,cn,cni, &
                             v11,v12,v13, &
                             v21,v22,v23, &
                             v31,v32,v33
      write (30000+mypnum,3) k, &
                             uu*aref, &
                             vv*aref, &
                             ww*aref, &
                             r*rhoref, &
                             u*aref, &
                             v*aref, &
                             w*aref, &
                             t*tref, &
                             sqrt(q2)*aref
    end if
    1 format (i8," : ",3(1x,f11.6), &
                    4x,f9.3,2(1x,f9.5), &
                    3x,3(1x,f9.5), &
                    4x,f9.5,1x,f9.5,2x,f11.7)
    2 format (i8," : ",f11.6,1x,f15.6, &
                    3(3x,3(1x,f15.6)))
    3 format (i8," : ",3(1x,f15.6), &
                      4x,6f15.6)
#endif
   !
  end do
  !
#ifdef DEBUG_THX_PIPE_INIT
  flush (10000+mypnum)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  error_message = "After dumping THX pipe SEM initialization"
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine thx_pipe_sem_initialization
!
!###############################################################################
!
pure &
function thx_pipe_turn_angle(x,y) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: x,y
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: numer,denom
  !
  !.. Local Constants ..
  real(wp), parameter :: c1 = sqrt(three) - 11.0_wp
  real(wp), parameter :: c2 = 635.0_wp - 74.0_wp*sqrt(three)
  real(wp), parameter :: c3 = 13.0_wp*sqrt(three) - two
  real(wp), parameter :: cx = half*(two - 13.0_wp*sqrt(three))
  !
continue
  !
  ! Compute the turn angle
  !
  ! NOTE: Using mathematica, I was able to prove that the argument to the
  !       sqrt in the denominator should always be positive and greater
  !       than zero if the coordinates [x,y] are a valid point within the
  !       pipe/injector geometry.
  !
  numer = c1 - two*y
  denom = sqrt( c2 + four*x*(c3+x) + four*y*(y-c1) )
  !
  return_value = acos( numer / denom ) * radians2degrees
  !
  ! Set the angle negative if it is in the SEM inflow region
  ! of the pipe/injector that is aligned with the x-axis.
  !
  if (x < cx) return_value = -return_value
  !
end function thx_pipe_turn_angle
!
!###############################################################################
!
pure &
function thx_pipe_centerline_y(angle,x,y) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: angle
  real(wp), intent(in) :: x
  real(wp), intent(in) :: y
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Constants ..
  real(wp), parameter :: c1  = -half*13.0_wp
  real(wp), parameter :: c2  = -one4*25.0_wp
  real(wp), parameter :: c3  = half*(two + sqrt(three))
  real(wp), parameter :: c4  = half*(sqrt(three) - 11.0_wp)
  real(wp), parameter :: c5  = -one8*three*sqrt(three)
  real(wp), parameter :: c6  = one4*sqrt(three)
  real(wp), parameter :: c30 = 30.0_wp
  !
continue
  !
  if (angle <= zero) then
    return_value = c1
  else if (angle >= c30) then
    return_value = c5 + c6*x + one4*y
   !return_value = c2 + half*c3*tan( (angle-c30)*degrees2radians )
  else
    return_value = c4 - c3*cos(angle*degrees2radians)
  end if
  !
end function thx_pipe_centerline_y
!
!###############################################################################
!
pure &
function thx_pipe_centerline_x(angle,yc,x) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: angle
  real(wp), intent(in) :: yc
  real(wp), intent(in) :: x
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Constants ..
  real(wp), parameter :: c1  = three/two
  real(wp), parameter :: c2  = sqrt(three)
  real(wp), parameter :: c3  = half*(two - 13.0_wp*sqrt(three))
  real(wp), parameter :: c4  = half*(two + sqrt(three))
  real(wp), parameter :: c5  = half*13.0_wp
  real(wp), parameter :: c30 = 30.0_wp
  !
continue
  !
  if (angle <= zero) then
    return_value = x
   !return_value = c3 - c5*tan(abs(angle)*degrees2radians)
  else if (angle >= c30) then
    return_value = c1 + c2*yc
  else
    return_value = c3 + c4*sin(angle*degrees2radians)
  end if
  !
end function thx_pipe_centerline_x
!
!###############################################################################
!
subroutine init_from_sem_profile(prof,pvavg,xyz,cv)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
  !
  !.. Formal Arguments ..
  type(sem_profile_t), intent(in) :: prof
  real(wp),            intent(in) :: pvavg(:)
  real(wp),            intent(in) :: xyz(:,:)
  real(wp),         intent(inout) :: cv(:,:)
  !
  !.. Local Scalars ..
  integer  :: k,n,n1,n2,npts
  real(wp) :: c1,c2,y,y1,y2
  real(wp) :: r,u,v,w,t,q2
  real(wp) :: rgas_over_gm1
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_from_sem_profile"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  npts = prof%npts
  !
  ! Precompute Rgas/(gamma-1)
  !
  rgas_over_gm1 = rgasref_nd / (gam - one)
  !
  do k = 1,size(xyz,dim=2)
    !
    y = xyz(2,k)
    !
    n1 = 1
    n2 = 2
    do n = 2,npts
      if (prof%y(n) > y) exit
      n1 = n
      n2 = max(n+1,npts)
    end do
    !
    if (n1 == npts) then
      !
      r = pvavg(nec)
      u = pvavg(nmx)
      v = pvavg(nmy)
      w = pvavg(nmz)
      t = pvavg(nee)
      !
    else
      !
      y1 = prof%y(n1)
      y2 = prof%y(n2)
      !
      c1 = (y-y2)/(y2-y1)
      c2 = (y-y1)/(y2-y1)
      !
      r = c2*prof%ravg(n2) - c1*prof%ravg(n1)
      u = c2*prof%uavg(n2) - c1*prof%uavg(n1)
      v = c2*prof%vavg(n2) - c1*prof%vavg(n1)
      w = c2*prof%wavg(n2) - c1*prof%wavg(n1)
      t = c2*prof%tavg(n2) - c1*prof%tavg(n1)
      !
    end if
    !
    q2 = u*u + v*v + w*w
    !
    cv(k,nec) = r
    cv(k,nmx) = r*u
    cv(k,nmy) = r*v
    cv(k,nmz) = r*w
    cv(k,nee) = r*( rgas_over_gm1*t + half*q2 )
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine init_from_sem_profile
!
!###############################################################################
!
subroutine sem(face_cv,face_gp)
  !
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  ! Synthetic Eddy Method
  !
  ! Based on
  ! Jarrin et al - Journal of Heat and Fluid Flow
  !
  ! Creates unsteady velocity fluctuations in a nozzle boundary layer
  ! based on a supplied turbulent kinetic energy profile.
  !
  ! Velocity and Reynolds stress profiles are read in.
  ! Convective velocity and turbulent length scale are input
  ! as variables in the bc specification.
  !
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  !.. Use Statements ..
  use ovar, only : d_t,gam,rgasref_nd,time_ref,itcur
  use ovar, only : use_improved_SEM_normalization
  use ovar, only : jet_init_func
  !
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,ierr,sem_idx,eddy_idx,num_eddies
  real(wp) :: ave,running_ave,time_scale,dt
  real(wp) :: c0,c1,eddy_vol,boxvol,lturb,dx,fsig,one_over_gm1
  real(wp) :: a11,a21,a22,a31,a32,a33,rho,temp,pres,ravg,tavg,q2
  real(wp) :: x,epx,elx,eix,tfx,ufluc,ufluc2,uavg,vx
  real(wp) :: y,epy,ely,eiy,tfy,vfluc,vfluc2,vavg,vy
  real(wp) :: z,epz,elz,eiz,tfz,wfluc,wfluc2,wavg,vz
  real(wp) :: st,ct,vr,vt,aoa,boa
  real(wp) :: boxmaxx,new_x
  !
  !.. Local Allocatable Arrays ..
  type(eddy_list_t), allocatable :: new_eddies(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "sem"
  !
  real(wp), parameter :: sqrt32 = sqrt(three*half)
  real(wp), parameter :: time_fraction = one5*one5
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (.not. sem_is_initialized) then
    call initialize_sem()
    sem_is_initialized = true
  end if
  !
  if (num_sem_boxes < 1) then
    call debug_timer(leaving_procedure,pname)
    return
  end if
  !
  ! Get the time step size for the synthetic eddies
  !
  dt = d_t
  if (allocated(sem_avg) .and. dt <= eps12) then
    !
    ! If d_t is zero, make the time step size for the synthetic eddies
    ! some random value, say 'a' such that eps6<=a<1, times a fraction
    ! (defined by parameter time_fraction) of the minimum convective
    ! time scale of all the SEM boxes.
    !
    call random_number(dt)
    dt = max(eps6,dt)
    !
    time_scale = minval( sem_box(:)%length_scale / &
                         sem_box(:)%convective_velocity )
    !
    dt = dt * time_scale * time_fraction
   !dt = dt * time_scale * time_fraction * time_ref
    !
  end if
  !
  one_over_gm1 = one / (gam - one)
  !
  ! ###########
  ! ###########
  ! ###  1  ###
  ! ###########
  ! ###########
  !
  !-----------------------------------------------------------------------
  ! Loop over each SEM boundary condition to apply fluctuations
  ! to the inflow velocity profile
  !-----------------------------------------------------------------------
  !
#ifdef DEBUG_SEM_EDDY_INFO
  do sem_idx = 1,num_sem_boxes
    !
    boxvol  = sem_box(sem_idx)%box_volume
   !lturb   = sem_box(sem_idx)%length_scale
    num_eddies = sem_box(sem_idx)%num_eddies
    !
    lturb = sem_box(sem_idx)%length_scale
    !
    c0 = boxvol/real(num_eddies,kind=wp)
    !
    do eddy_idx = 1,num_eddies
      !
      eix = sem_box(sem_idx)%eddy(eddy_idx)%x_intensity
      eiy = sem_box(sem_idx)%eddy(eddy_idx)%y_intensity
      eiz = sem_box(sem_idx)%eddy(eddy_idx)%z_intensity
      !
      elx = one/sem_box(sem_idx)%eddy(eddy_idx)%x_frequency
      ely = one/sem_box(sem_idx)%eddy(eddy_idx)%y_frequency
      elz = one/sem_box(sem_idx)%eddy(eddy_idx)%z_frequency
      !
      eddy_vol = sem_box(sem_idx)%eddy(eddy_idx)%volume
      !
      write (iout,61) sem_idx,eddy_idx, &
                      boxvol,num_eddies,c0, &
                      lturb,elx,ely,elz, &
                      elx*ely*elz,one/(eddy_vol*eddy_vol), &
                      eix,eiy,eiz
      !
    end do
    !
  end do
  !
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,"After writing eddy info.")
  61 format ("SEM Box # ",i0,", Eddy # ",i4,/, &
             2x,"Box volume = ",es14.6,/, &
             2x,"num_eddies = ",i0,/, &
             2x,"        c0 = ",es14.6,/, &
             2x,"Length Scale [Box,Eddy-x,y,z]: ",4es14.6,/, &
             2x,"Eddy Volume [l^3,eddy_vol]: ",2es14.6,/, &
             2x,"Eddy Intensities [x,y,z]: ",3es14.6,/)
#endif
  !
  loop_over_my_sem_boxes: &
  do sem_idx = 1,num_sem_boxes
    !
    boxvol  = sem_box(sem_idx)%box_volume
   !lturb   = sem_box(sem_idx)%length_scale
    num_eddies = sem_box(sem_idx)%num_eddies
    !
    time_scale = sem_box(sem_idx)%length_scale / &
                 sem_box(sem_idx)%convective_velocity
    !
    aoa = sem_box(sem_idx)%alpha_aoa * degrees2radians
    boa = sem_box(sem_idx)%beta_aoa  * degrees2radians
    !
    ! Pre-compute the parts of the normalization constant and
    ! shape/tent function that can be factored out of the summation
    ! over the eddies.
   !!
   !c0 = boxvol*sqrt(half*three)*half*three/(lturb**6)/real(num_eddies,kind=wp)
    !
    ! Pre-compute the part from the normalization constant
    !
    c0 = boxvol/real(num_eddies,kind=wp)
    !
    ! Loop through each face point on this SEM boundary and sum up the
    ! contributions of each eddy to the velocity fluctuation at that point
    !
    loop_over_sem_pts: &
    do n = 1,my_sem_pts(sem_idx)%num_fp
      !
      ! Get the index of this face point within the solver face arrays
      !
      k = my_sem_pts(sem_idx)%sem_fp(n)%idx
      !
      x = my_sem_pts(sem_idx)%sem_fp(n)%x
      y = my_sem_pts(sem_idx)%sem_fp(n)%y
      z = my_sem_pts(sem_idx)%sem_fp(n)%z
      !
      running_ave = my_sem_pts(sem_idx)%sem_fp(n)%running_ave
      !
      ufluc = zero
      vfluc = zero
      wfluc = zero
      !
      ave = zero
      !
      ! Loop over the eddies that contribute to the fluctuation at the
      ! current face point, and sum up these contributions to get the
      ! total velocity fluctuation at this face point.
      !
      loop_over_eddy_list: &
      do i = 1,my_sem_pts(sem_idx)%sem_fp(n)%num_eddies
        !
        eddy_idx = my_sem_pts(sem_idx)%sem_fp(n)%eddy_list(i)
        !
        epx = sem_box(sem_idx)%eddy(eddy_idx)%x
        epy = sem_box(sem_idx)%eddy(eddy_idx)%y
        epz = sem_box(sem_idx)%eddy(eddy_idx)%z
        !
        eix = sem_box(sem_idx)%eddy(eddy_idx)%x_intensity
        eiy = sem_box(sem_idx)%eddy(eddy_idx)%y_intensity
        eiz = sem_box(sem_idx)%eddy(eddy_idx)%z_intensity
        !
        elx = sem_box(sem_idx)%eddy(eddy_idx)%x_frequency
        ely = sem_box(sem_idx)%eddy(eddy_idx)%y_frequency
        elz = sem_box(sem_idx)%eddy(eddy_idx)%z_frequency
        !
        eddy_vol = sem_box(sem_idx)%eddy(eddy_idx)%volume
        !
       !tfx = sqrt32 * tentfun( abs(x-epx) , one/elx )
       !tfy = sqrt32 * tentfun( abs(y-epy) , one/ely )
       !tfz = sqrt32 * tentfun( abs(z-epz) , one/elz )
        tfx = sqrt32 * max( one - abs(x-epx)*elx , zero )
        tfy = sqrt32 * max( one - abs(y-epy)*ely , zero )
        tfz = sqrt32 * max( one - abs(z-epz)*elz , zero )
        !
        fsig = tfx*tfy*tfz
        !
        if (use_improved_SEM_normalization) then
          ave = ave + fsig*fsig
        else
          fsig = fsig * eddy_vol
        end if
        !
        ufluc = ufluc + fsig*eix
        vfluc = vfluc + fsig*eiy
        wfluc = wfluc + fsig*eiz
        !
      end do loop_over_eddy_list
      !
      ! Normalize the velocity fluctuations from the synthetic eddies
      !
      if (use_improved_SEM_normalization) then
        !
        c1 = dt/time_scale
        !
        running_ave = c1*ave + (one-c1)*running_ave
        !
        ufluc2 = zero
        vfluc2 = zero
        wfluc2 = zero
        !
        if (ufluc /= zero) ufluc2 = ufluc / sqrt(running_ave)
        if (vfluc /= zero) vfluc2 = vfluc / sqrt(running_ave)
        if (wfluc /= zero) wfluc2 = wfluc / sqrt(running_ave)
        !
        my_sem_pts(sem_idx)%sem_fp(n)%running_ave = running_ave
        !
      else
        !
        ufluc2 = c0*ufluc
        vfluc2 = c0*vfluc
        wfluc2 = c0*wfluc
        !
      end if
      !
      ! Compute the final velocity fluctuations for this face point
      !
      a11 = my_sem_pts(sem_idx)%sem_fp(n)%a11
      a21 = my_sem_pts(sem_idx)%sem_fp(n)%a21
      a22 = my_sem_pts(sem_idx)%sem_fp(n)%a22
      a31 = my_sem_pts(sem_idx)%sem_fp(n)%a31
      a32 = my_sem_pts(sem_idx)%sem_fp(n)%a32
      a33 = my_sem_pts(sem_idx)%sem_fp(n)%a33
      !
      ufluc = ufluc2*a11
      vfluc = ufluc2*a21 + vfluc2*a22
      wfluc = ufluc2*a31 + vfluc2*a32 + wfluc2*a33
      !
      ravg = my_sem_pts(sem_idx)%sem_fp(n)%ravg
      uavg = my_sem_pts(sem_idx)%sem_fp(n)%uavg
      vavg = my_sem_pts(sem_idx)%sem_fp(n)%vavg
      wavg = my_sem_pts(sem_idx)%sem_fp(n)%wavg
      tavg = my_sem_pts(sem_idx)%sem_fp(n)%tavg
      !
      vx = ufluc + uavg
      vy = vfluc + vavg
      vz = wfluc + wavg
      !
      q2 = vx*vx + vy*vy + vz*vz
      !
      ! One version of Jims code uses rho_inf and temp_inf, and the other
      ! uses rho_avg and temp_avg from the inflow profile
      !
      if (jet_init_func == -1) then
        !
        ! Use interior pressure, SEM density, and SEM velocities
        !
        pres = face_gp(k,npres,1)
        rho = ravg
        temp = pres / (rho*rgasref_nd)
        !
      else
        !
        ! Use SEM density, SEM velocities, and SEM temperature/presssure
        !
        rho  = ravg
        temp = tavg
       !rho  = one  ! = rhoref / rhoref
       !temp = one  ! = tref / tref
        pres = rho*rgasref_nd*temp
        !
      end if
      !
      face_cv(k,nec,2) = rho
      face_cv(k,nmx,2) = rho*vx
      face_cv(k,nmy,2) = rho*vy
      face_cv(k,nmz,2) = rho*vz
      face_cv(k,nee,2) = pres*one_over_gm1 + half*rho*q2
      !
      face_gp(k,npres,2) = pres
      face_gp(k,ntemp,2) = temp
      !
      if (jet_init_func == -2) then
        !
        ! Copy the SEM boundary states to the interior face points
        !
        face_cv(k,nec,1) = face_cv(k,nec,2)
        face_cv(k,nmx,1) = face_cv(k,nmx,2)
        face_cv(k,nmy,1) = face_cv(k,nmy,2)
        face_cv(k,nmz,1) = face_cv(k,nmz,2)
        face_cv(k,nee,1) = face_cv(k,nee,2)
        !
        face_gp(k,npres,1) = face_gp(k,npres,2)
        face_gp(k,ntemp,1) = face_gp(k,ntemp,2)
        !
        call update_thermo_properties(face_cv(k,:,1),face_gp(k,:,1), &
                                      caller_name=pname)
        !
        call update_transport_properties(face_cv(k,:,1),face_gp(k,:,1))
        !
      end if
      !
#ifdef ENABLE_OUTPUT_TIME_HISTORY
      my_sem_pts(sem_idx)%sem_fp(n)%time_history(nec ,itcur) = rho
      my_sem_pts(sem_idx)%sem_fp(n)%time_history(nmx ,itcur) = vx
      my_sem_pts(sem_idx)%sem_fp(n)%time_history(nmy ,itcur) = vy
      my_sem_pts(sem_idx)%sem_fp(n)%time_history(nmz ,itcur) = vz
      my_sem_pts(sem_idx)%sem_fp(n)%time_history(nee ,itcur) = face_cv(k,nee,2)
      my_sem_pts(sem_idx)%sem_fp(n)%time_history(nq+1,itcur) = temp
#endif
      !
    end do loop_over_sem_pts
    !
  end do loop_over_my_sem_boxes
  !
  if (allocated(sem_avg)) then
    call update_sem_ave(dt,face_cv(:,:,2))
  end if
  !
!#define TROUBLESHOOT_PROFILE_INTERPOLATION
#ifdef TROUBLESHOOT_PROFILE_INTERPOLATION
  call output_sem_averages()
  error_message = "#define TROUBLESHOOT_PROFILE_INTERPOLATION"
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
#endif
  !
  ! ###########
  ! ###########
  ! ###  2  ###
  ! ###########
  ! ###########
  !
  !-----------------------------------------------------------------------
  ! Convect the eddies through the box
  !-----------------------------------------------------------------------
  !
  allocate ( new_eddies(1:num_sem_boxes) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"new_eddies",2,__LINE__,__FILE__,ierr,error_message)
  !
  loop_over_all_sem_boxes: &
  do sem_idx = 1,num_sem_boxes
    !
    ! Get the axial location of the outflow plane of the current SEM box
    !
    boxmaxx = sem_box(sem_idx)%boxmax(1)
    !
    ! Get the axial distance traveled by the eddies over the past time step
    !
    dx = sem_box(sem_idx)%convective_velocity * dt
    !
    ! Loop over the eddies in this SEM box
    !
    loop_over_eddies_in_box: &
    do eddy_idx = 1,sem_box(sem_idx)%num_eddies
      !
      x = sem_box(sem_idx)%eddy(eddy_idx)%x
      !
      ! Compute the new axial location of the current eddy
      ! NOTE: THIS NEEDS TO BE ALTERED IF THE INFLOW PROFILE VELOCITY IS
      !       NOT X-DIRECTION ONLY.
      !
      new_x = x + dx
      !
      if (new_x > boxmaxx) then
        !
        ! This eddy has left the current SEM box, so we need to replace it
        ! with an eddy at a new random y-z location along the SEM box inflow
        ! plane.
        !
        ! First, we need to save the indices for this eddy
        ! to the list of new eddies.
        !
        call add_new_eddy_to_list(new_eddies(sem_idx),eddy_idx)
        !
      else
        !
        ! This eddy is still within the SEM box, so we just need to update
        ! its axial location.
        !
        sem_box(sem_idx)%eddy(eddy_idx)%x = new_x
        !
      end if
      !
    end do loop_over_eddies_in_box
    !
  end do loop_over_all_sem_boxes
  !
  ! ###########
  ! ###########
  ! ###  3  ###
  ! ###########
  ! ###########
  !
  !-----------------------------------------------------------------------
  ! Generate new random eddies for all eddies that have left the SEM box
  ! NOTE: All eddies are regenerated at the inflow plane within the solver
  !       so the argument 'generate_at_xmin' needs to be true.
  !-----------------------------------------------------------------------
  !
  call generate_new_eddies(new_eddies,generate_at_xmin=true,dt=dt)
  !
  ! Now we need to update the face point eddy-lists and the
  ! eddy face-point-lists due to the "new eddies" created to
  ! replace the "old eddies" that convected out of their SEM box.
  !
  ! Loop over the list of new eddies, and update the list of contributing eddies
  ! for each face point (of this partition on an SEM boundary) that contained
  ! an eddy that left its SEM box or will contain one of the new eddies.
  !
  ! ###########
  ! ###########
  ! ###  4  ###
  ! ###########
  ! ###########
  !
  !-----------------------------------------------------------------------
  ! Remove all the "old eddies" that exited their SEM box
  !-----------------------------------------------------------------------
  !
  loop_to_remove_old_eddies: &
  do sem_idx = 1,num_sem_boxes
    !
    do i = 1,new_eddies(sem_idx)%num_eddies
      !
      eddy_idx = new_eddies(sem_idx)%eddy_list(i)
      !
      do n = 1,sem_box(sem_idx)%eddy(eddy_idx)%num_fp
        !
        ! Get the index of the next SEM face point affected by this eddy
        !
        k = sem_box(sem_idx)%eddy(eddy_idx)%fp_list(n)
        !
        ! Remove the "old eddy" from the list of
        ! eddies affecting this face point
        !
        call remove_eddy_from_fp_eddy_list(my_sem_pts(sem_idx)%sem_fp(k), &
                                           eddy_idx)
        !
      end do
      !
      ! Finally, reset/delete the list of face points affected by this eddy.
      ! NOTE : The next loop will recreate this list for the new set of face
      !        points that are affected by the new random location of this eddy.
      !
      call delete_eddy_fp_list(sem_box(sem_idx)%eddy(eddy_idx))
      !
    end do
    !
  end do loop_to_remove_old_eddies
  !
  ! ###########
  ! ###########
  ! ###  5  ###
  ! ###########
  ! ###########
  !
  !-----------------------------------------------------------------------
  ! Find any SEM face points that will be affected
  ! by each of the "new eddies", then:
  !    1) add this eddy to the face point eddy-list, and
  !    2) add this face point to the eddy face-point-list.
  !-----------------------------------------------------------------------
  !
  call update_eddy_lists(new_eddies)
  !
  ! Deallocate new_eddies since we no longer need it
  !
  deallocate ( new_eddies , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"new_eddies",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine sem
!
!###############################################################################
!
subroutine update_sem_ave(dt,cv,alpha,beta)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd,itcur
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: dt
  real(wp), intent(in) :: cv(:,:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: alpha
  real(wp), optional, intent(in) :: beta
  !
  !.. Local Scalars ..
  integer  :: k,n,sem_idx
  real(wp) :: current_time,a,b
  real(wp) :: u,v,w,q2,t,p,gm1,rgas_inv
  real(wp) :: r,e,r_inv
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_sem_ave"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  gm1 = gam - one
  rgas_inv = one / rgasref_nd
  !
  if (present(alpha) .and. present(beta)) then
    a = alpha
    b = beta
  else
    if (sem_prev_ave_time < zero) then
      current_time = zero
      a = zero
      b = one
    else
      current_time = sem_prev_ave_time + dt
      a = sem_prev_ave_time / current_time
      b = one - a
    end if
  end if
  !
  do sem_idx = 1,num_sem_boxes
    !
    do n = 1,my_sem_pts(sem_idx)%num_fp
      !
      k = my_sem_pts(sem_idx)%sem_fp(n)%idx
      !
      r = cv(k,nec)
      r_inv = one / r
      !
      u = cv(k,nmx) * r_inv
      v = cv(k,nmy) * r_inv
      w = cv(k,nmz) * r_inv
      e = cv(k,nee)
      !
      q2 = u*u + v*v + w*w
      !
      p = gm1*(e - half*r*q2)
      t = p * rgas_inv * r_inv
      !
      call update_sem_ave_fp(a,b,r,u,v,w,e,t,p, &
                             sem_avg(sem_idx)%fp_avg(n))
      !
    end do
    !
  end do
  !
  ! Update the previous average time
  !
  sem_prev_ave_time = current_time
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_sem_ave
!
!###############################################################################
!
pure &
subroutine update_sem_ave_fp(a,b,r,u,v,w,e,t,p,this)
  !
  !.. Formal Arguments ..
  real(wp),              intent(in) :: a,b,r,u,v,w,e,t,p
  type(sem_fp_avg_t), intent(inout) :: this
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_sem_ave_fp"
  !
continue
  !
  this%r   = a * this%r   + b * r
  this%ru  = a * this%ru  + b * r*u
  this%rv  = a * this%rv  + b * r*v
  this%rw  = a * this%rw  + b * r*w
  this%e   = a * this%e   + b * e
  this%t   = a * this%t   + b * t
  this%u   = a * this%u   + b * u
  this%v   = a * this%v   + b * v
  this%w   = a * this%w   + b * w
  this%p   = a * this%p   + b * p
  this%uu  = a * this%uu  + b * u*u
  this%vv  = a * this%vv  + b * v*v
  this%ww  = a * this%ww  + b * w*w
  this%uv  = a * this%uv  + b * u*v
  this%uw  = a * this%uw  + b * u*w
  this%vw  = a * this%vw  + b * v*w
  this%ut  = a * this%ut  + b * u*t
  this%vt  = a * this%vt  + b * v*t
  this%wt  = a * this%wt  + b * w*t
  this%rt  = a * this%rt  + b * r*t
  this%ruu = a * this%ruu + b * r*u*u
  this%rvv = a * this%rvv + b * r*v*v
  this%rww = a * this%rww + b * r*w*w
  this%ruv = a * this%ruv + b * r*u*v
  this%ruw = a * this%ruw + b * r*u*w
  this%rvw = a * this%rvw + b * r*v*w
  this%rut = a * this%rut + b * r*u*t
  this%rvt = a * this%rvt + b * r*v*t
  this%rwt = a * this%rwt + b * r*w*t
  !
end subroutine update_sem_ave_fp
!
!###############################################################################
!
subroutine generate_new_eddies(new_eddies,generate_at_xmin,dt)
  !
  !.. Use Statements ..
  use ovar, only : SEM_eddy_volume_is_spherical
  !
  !.. Formal Arguments ..
  type(eddy_list_t), intent(in) :: new_eddies(:)
  logical(lk),       intent(in) :: generate_at_xmin
  real(wp),          intent(in) :: dt
  !
  !.. Local Scalars ..
  integer  :: n,num_eddies,sem_idx,eddy_idx
  real(wp) :: eddy_vol,x_coef,vol_coef,x_offset
  real(wp) :: x,xc,eix,elx,boxminx,boxmaxx,box_dx
  real(wp) :: y,yc,eiy,ely,boxminy,boxmaxy,box_dy
  real(wp) :: z,zc,eiz,elz,boxminz,boxmaxz,box_dz
  real(wp) :: r,t,rmin2,rmax2
  !
  !.. Local Arrays ..
  real(wp) :: new_xyz(1:3)
  real(wp) :: new_int(1:3)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "generate_new_eddies"
  !
  ! USE_X_OFFSET : constant to switch on/off the use of a tiny x-offset when
  !                regenerating eddies that have left an SEM box so as to not
  !                create "groups" of eddies that have the same x-location
  !                and regenerate together.
 !real(wp), parameter :: use_x_offset = zero
  real(wp), parameter :: use_x_offset = one
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  x_coef = merge( zero , one , generate_at_xmin )
  !
  vol_coef = merge( pi/six , one , SEM_eddy_volume_is_spherical )
  !
  loop_over_sem_boxes: &
  do sem_idx = 1,num_sem_boxes
    !
    xc = sem_box(sem_idx)%center%x
    yc = sem_box(sem_idx)%center%y
    zc = sem_box(sem_idx)%center%z
    !
    num_eddies = new_eddies(sem_idx)%num_eddies
    !
    x_offset = half*dt*sem_box(sem_idx)%convective_velocity * use_x_offset
    !
    if (.not. allocated(sem_box(sem_idx)%eddy)) then
      allocate ( sem_box(sem_idx)%eddy(1:num_eddies) )
    end if
    !
    ! Get the min/max extents of the current SEM box
    !
    boxminx = sem_box(sem_idx)%boxmin(1)
    boxmaxx = sem_box(sem_idx)%boxmax(1)
    !
    if (sem_box(sem_idx)%is_cylindrical) then
      !
      boxminy = sem_box(sem_idx)%boxmin(4)
      boxminz = sem_box(sem_idx)%boxmin(5)
      !
      boxmaxy = sem_box(sem_idx)%boxmax(4)
      boxmaxz = sem_box(sem_idx)%boxmax(5)
      !
    else
      !
      boxminy = sem_box(sem_idx)%boxmin(2)
      boxminz = sem_box(sem_idx)%boxmin(3)
      !
      boxmaxy = sem_box(sem_idx)%boxmax(2)
      boxmaxz = sem_box(sem_idx)%boxmax(3)
      !
    end if
    !
    box_dx = (boxmaxx - boxminx)*x_coef
    box_dy =  boxmaxy - boxminy
    box_dz =  boxmaxz - boxminz
    !
    loop_over_new_eddies_in_box: &
    do n = 1,num_eddies
      !
      eddy_idx = new_eddies(sem_idx)%eddy_list(n)
      !
      ! Create new random values to get the x,y,z coordinates for the new eddy
      !
      call random_number(new_xyz)
      !
      x = boxminx + new_xyz(1)*box_dx + x_offset*(new_xyz(1)-one)
      !
      if (sem_box(sem_idx)%is_cylindrical) then
        !
        ! Random theta coordinate is same as x,y,z
        !
        t = boxminz + new_xyz(3)*box_dz
        !
        ! Random radial coordinate need to be handled a little differently.
        ! Naively picking a random point between boxminy and boxmaxy, i.e.,
        !     r = boxminy + new_xyz(2)*box_dy
        ! will result in more eddies near boxminy and fewer eddies near boxmaxy.
        ! To generate evenly distributed eddies within the SEM annulus box, we
        ! need to get our "random" radius using the equation
        !     r = sqrt( new_xyz(2)*(R2^2 + R1^2) + R1^2 )
        ! where R2=boxmaxy and R1=boxminy.
        !
        rmin2 = boxminy*boxminy
        rmax2 = boxmaxy*boxmaxy
        !
        r = sqrt( new_xyz(2)*(rmax2 - rmin2) + rmin2 )
        !
        ! Convert r and theta to y and z
        !
        y = yc + r*cos(t)
        z = zc + r*sin(t)
        !
      else
        !
        y = boxminy + new_xyz(2)*box_dy
        z = boxminz + new_xyz(3)*box_dz
        !
        ! Convert y and z to r and theta
        !
        call yz_to_rt(y,z,r,t,yc,zc)
        !
      end if
      !
      ! Save the new random coordinates for the new eddy in sem_box
      !
      sem_box(sem_idx)%eddy(eddy_idx)%x = x
      sem_box(sem_idx)%eddy(eddy_idx)%y = y
      sem_box(sem_idx)%eddy(eddy_idx)%z = z
      sem_box(sem_idx)%eddy(eddy_idx)%r = r
      sem_box(sem_idx)%eddy(eddy_idx)%t = t
      !
      ! Create new random intensities for the new eddy
      !
      call random_number(new_int)
      !
      eix = merge(-one,one,new_int(1)<=half)
      eiy = merge(-one,one,new_int(2)<=half)
      eiz = merge(-one,one,new_int(3)<=half)
      !
      sem_box(sem_idx)%eddy(eddy_idx)%x_intensity = eix
      sem_box(sem_idx)%eddy(eddy_idx)%y_intensity = eiy
      sem_box(sem_idx)%eddy(eddy_idx)%z_intensity = eiz
      !
      !####################################################################
      !####################################################################
      !###                                                              ###
      !###                     TEMPORARY !!!!!!!                        ###
      !###                                                              ###
      !### |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ###
      !### VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV ###
      !###                                                              ###
      !####################################################################
      !####################################################################
      !
      ! For now, each SEM box has its own specified length scale, so set
      ! the eddy length scales to the length scale of the SEM box
      !
      elx = sem_box(sem_idx)%length_scale
      ely = sem_box(sem_idx)%length_scale
      elz = sem_box(sem_idx)%length_scale
      !
      !####################################################################
      !####################################################################
      !###                                                              ###
      !### AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA ###
      !### |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ###
      !###                                                              ###
      !####################################################################
      !####################################################################
      !
      ! Store the eddy length scales in each direction as the eddy frequencies
      !
      sem_box(sem_idx)%eddy(eddy_idx)%x_frequency = one / elx
      sem_box(sem_idx)%eddy(eddy_idx)%y_frequency = one / ely
      sem_box(sem_idx)%eddy(eddy_idx)%z_frequency = one / elz
      !
      ! Compute the volume of this eddy
      ! NOTE: The input parameter SEM_eddy_volume_is_spherical determines
      !       whether the variable vol_coef is equal to either
      !       (4/3)*PI*(1/2)**2 or 1 depending on whether the eddy is
      !       treated as sphere-like or box-like in volume, respectively.
      !       If a box-like eddy has a volume of l**3, the equivalent
      !       volume of a sphere-like eddy is (4/3)*pi*((l/2)**3).
      !
      eddy_vol = vol_coef*elx*ely*elz
      !
      ! Store the eddy volume as the reciprocal of the volume
      !
      sem_box(sem_idx)%eddy(eddy_idx)%volume = one / eddy_vol
      !
    end do loop_over_new_eddies_in_box
    !
  end do loop_over_sem_boxes
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine generate_new_eddies
!
!###############################################################################
!
subroutine initialize_eddy_lists()
  !
  !.. Local Scalars ..
  integer  :: i,ierr,n,sem_idx,eddy_idx
  real(wp) :: lturb,xe,ye,ze,xp,yp,zp,lx,ly,lz,ds
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_eddy_lists"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Make sure none of the sem_box(:)%eddy(:)%fp_list or
  ! my_sem_pts(:)%sem_fp(:)%eddy_list allocatable arrays
  ! are allocated.
  !
  do sem_idx = 1,num_sem_boxes
    !
    do eddy_idx = 1,sem_box(sem_idx)%num_eddies
      sem_box(sem_idx)%eddy(eddy_idx)%num_fp = 0
      if (allocated(sem_box(sem_idx)%eddy(eddy_idx)%fp_list)) then
        deallocate ( sem_box(sem_idx)%eddy(eddy_idx)%fp_list , stat=ierr )
      end if
    end do
    !
    do n = 1,my_sem_pts(sem_idx)%num_fp
      my_sem_pts(sem_idx)%sem_fp(n)%num_eddies = 0
      if (allocated(my_sem_pts(sem_idx)%sem_fp(n)%eddy_list)) then
        deallocate ( my_sem_pts(sem_idx)%sem_fp(n)%eddy_list , stat=ierr )
      end if
    end do
    !
  end do
  !
  ! Now create all the eddy face-point-lists and face point eddy-lists
  !
  loop_over_boxes: &
  do sem_idx = 1,num_sem_boxes
    !
    loop_over_eddies_in_box: &
    do eddy_idx = 1,sem_box(sem_idx)%num_eddies
      !
      ! Get the coordinates for this eddy
      !
      xe = sem_box(sem_idx)%eddy(eddy_idx)%x
      ye = sem_box(sem_idx)%eddy(eddy_idx)%y
      ze = sem_box(sem_idx)%eddy(eddy_idx)%z
      !
      ! Get the directional length scales for this eddy
      !
      lx = one / sem_box(sem_idx)%eddy(eddy_idx)%x_frequency
      ly = one / sem_box(sem_idx)%eddy(eddy_idx)%y_frequency
      lz = one / sem_box(sem_idx)%eddy(eddy_idx)%z_frequency
      !
      lturb = max(lx,ly,lz)
      !
      loop_over_fp_in_box: &
      do n = 1,my_sem_pts(sem_idx)%num_fp
        !
        ! Get the coordinates of this face point
        !
        xp = my_sem_pts(sem_idx)%sem_fp(n)%x
        yp = my_sem_pts(sem_idx)%sem_fp(n)%y
        zp = my_sem_pts(sem_idx)%sem_fp(n)%z
        !
        ! Get the minimum distance from this face point
        ! to the line that will be traversed by the eddy.
        ! NOTE: With the current assumptions that:
        !         1) the SEM inflow boundary is a constant-x plane and
        !         2) the eddy convective velocity is only in the x-direction,
        !       the minimum distance between the face point and the line
        !       traversed by the eddy is simply the distance between the
        !       face point and the eddy when the x-coordinate of the eddy is
        !       equal to the x-coordinate of the face point, i.e., xe == xp.
        !
       !ds = hypot( yp-ye , zp-ze )
        ds = norm2( [ yp-ye , zp-ze ] )
       !ds = sqrt( (yp-ye)**2 + (zp-ze)**2 )
       !ds = sqrt( (yp-ye)*(yp-ye) + (zp-ze)*(zp-ze) )
        !
        if (ds <= lturb) then
          !
          ! Add this face point to the list of face points affected by this eddy
          !
          call add_fp_to_eddy_fp_list(sem_box(sem_idx)%eddy(eddy_idx),n)
          !
          !
          ! Add this eddy to the list of eddies affecting this face point
          !
          call add_eddy_to_fp_eddy_list(my_sem_pts(sem_idx)%sem_fp(n), &
                                        eddy_idx)
          !
        end if
        !
      end do loop_over_fp_in_box
      !
    end do loop_over_eddies_in_box
    !
  end do loop_over_boxes
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine initialize_eddy_lists
!
!###############################################################################
!
subroutine update_eddy_lists(new_eddies)
  !
  !.. Formal Arguments ..
  type(eddy_list_t), intent(in) :: new_eddies(:)
  !
  !.. Local Scalars ..
  integer  :: i,n,sem_idx,eddy_idx
  real(wp) :: lturb,xe,ye,ze,xp,yp,zp,lx,ly,lz,ds
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_eddy_lists"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  loop_to_add_new_eddies: &
  do sem_idx = 1,num_sem_boxes
    !
    loop_over_new_eddies_in_box: &
    do i = 1,new_eddies(sem_idx)%num_eddies
      !
      eddy_idx = new_eddies(sem_idx)%eddy_list(i)
      !
      ! Get the coordinates for this eddy
      !
      xe = sem_box(sem_idx)%eddy(eddy_idx)%x
      ye = sem_box(sem_idx)%eddy(eddy_idx)%y
      ze = sem_box(sem_idx)%eddy(eddy_idx)%z
      !
      ! Get the directional length scales for this eddy
      !
      lx = one / sem_box(sem_idx)%eddy(eddy_idx)%x_frequency
      ly = one / sem_box(sem_idx)%eddy(eddy_idx)%y_frequency
      lz = one / sem_box(sem_idx)%eddy(eddy_idx)%z_frequency
      !
      lturb = max(lx,ly,lz)
      !
      do n = 1,my_sem_pts(sem_idx)%num_fp
        !
        ! Get the coordinates of this face point
        !
        xp = my_sem_pts(sem_idx)%sem_fp(n)%x
        yp = my_sem_pts(sem_idx)%sem_fp(n)%y
        zp = my_sem_pts(sem_idx)%sem_fp(n)%z
        !
        ! Get the minimum distance from this face point
        ! to the line that will be traversed by the eddy.
        ! NOTE: With the current assumptions that:
        !         1) the SEM inflow boundary is a constant-x plane and
        !         2) the eddy convective velocity is only in the x-direction,
        !       the minimum distance between the face point and the line
        !       traversed by the eddy is simply the distance between the
        !       face point and the eddy when the x-coordinate of the eddy is
        !       equal to the x-coordinate of the face point, i.e., xe == xp.
        !
       !ds = hypot( yp-ye , zp-ze )
        ds = norm2( [ yp-ye , zp-ze ] )
       !ds = sqrt( (yp-ye)**2 + (zp-ze)**2 )
       !ds = sqrt( (yp-ye)*(yp-ye) + (zp-ze)*(zp-ze) )
        !
        if (ds <= lturb) then
          !
          ! Add this face point to the list of face points affected by this eddy
          !
          call add_fp_to_eddy_fp_list(sem_box(sem_idx)%eddy(eddy_idx),n)
          !
          !
          ! Add this eddy to the list of eddies affecting this face point
          !
          call add_eddy_to_fp_eddy_list(my_sem_pts(sem_idx)%sem_fp(n), &
                                        eddy_idx)
          !
        end if
        !
      end do
      !
    end do loop_over_new_eddies_in_box
    !
  end do loop_to_add_new_eddies
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_eddy_lists
!
!###############################################################################
!
pure &
subroutine add_fp_to_eddy_fp_list(this,fp_idx)
  !
  !.. Formal Arguments ..
  type(eddy_t), intent(inout) :: this
  integer,         intent(in) :: fp_idx
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "add_fp_to_eddy_fp_list"
  !
continue
  !
  if (.not. allocated(this%fp_list)) then
    ! Allocate fp_list to an initial array of 20 zeros
    this%fp_list = init_array(20,0) ! F2003 AUTO-REALLOC
    ! Initialize num_fp to zero
    this%num_fp = 0
  end if
  !
  n = this%num_fp ! current number of face points in fp_list
  !
  ! Make sure that fp_list is large enough to accommodate one more face point
  !
  if (n+1 > size(this%fp_list)) then
    ! Use F2003 auto-reallocation to increase the size of fp_list
    ! by appending 20 more zeros to the end of the array
    ! F2003 AUTO-REALLOC
    this%fp_list = [ this%fp_list , init_array(20,0) ]
  end if
  !
  ! Add the current face point to the list
  ! of face points affected by this eddy.
  !
  this%fp_list(n+1) = fp_idx
  this%num_fp = n+1
  !
end subroutine add_fp_to_eddy_fp_list
!
!###############################################################################
!
pure &
subroutine delete_eddy_fp_list(this)
  !
  !.. Formal Arguments ..
  type(eddy_t), intent(inout) :: this
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "delete_eddy_fp_list"
  !
continue
  !
  ! Delete/reset the old fp_list array by re-allocating
  ! it to an initial array of 20 zeros
  !
  this%fp_list = init_array(20,0) ! F2003 AUTO-REALLOC
  !
  ! Initialize num_fp to zero
  !
  this%num_fp = 0
  !
end subroutine delete_eddy_fp_list
!
!###############################################################################
!
pure &
subroutine add_eddy_to_fp_eddy_list(this,eddy_idx)
  !
  !.. Formal Arguments ..
  type(sem_fp_t), intent(inout) :: this
  integer,           intent(in) :: eddy_idx
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "add_eddy_to_fp_eddy_list"
  !
continue
  !
  if (.not. allocated(this%eddy_list)) then
    ! Allocate eddy_list to an initial array of 20 zeros
    this%eddy_list = init_array(20,0) ! F2003 AUTO-REALLOC
    ! Initialize num_eddies to zero
    this%num_eddies = 0
  end if
  !
  n = this%num_eddies ! current number of eddies in eddy_list
  !
  ! Make sure that eddy_list is large enough to accommodate one more eddy
  !
  if (n+1 > size(this%eddy_list)) then
    ! Use F2003 auto-reallocation to increase the size of eddy_list
    ! by appending 20 more zeros to the end of the array
    ! F2003 AUTO-REALLOC
    this%eddy_list = [ this%eddy_list , init_array(20,0) ]
  end if
  !
  ! Add the current eddy to the list of eddies affecting this face point
  !
  this%eddy_list(n+1) = eddy_idx
  this%num_eddies = n+1
  !
end subroutine add_eddy_to_fp_eddy_list
!
!###############################################################################
!
pure &
subroutine remove_eddy_from_fp_eddy_list(this,eddy_idx)
  !
  !.. Formal Arguments ..
  type(sem_fp_t), intent(inout) :: this
  integer,           intent(in) :: eddy_idx
  !
  !.. Local Scalars ..
  integer :: i,n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "remove_eddy_from_fp_eddy_list"
  !
continue
  !
  ! Find the location of this eddy within the list of eddies for this face point
  !
  i = findloc( this%eddy_list , eddy_idx , dim=1 )
  !
  ! Get the index for the last eddy in the eddy list
  !
  n = this%num_eddies
  !
  ! Replace the old eddy that left the box with the last eddy in the list
  !
  this%eddy_list(i) = this%eddy_list(n)
  !
  ! Zero out what was the last eddy in the list
  !
  this%eddy_list(n) = 0
  !
  ! Subtract 1 from the size of the eddy list
  !
  this%num_eddies = n-1
  !
end subroutine remove_eddy_from_fp_eddy_list
!
!###############################################################################
!
pure &
subroutine add_new_eddy_to_list(this,eddy_idx)
  !
  !.. Formal Arguments ..
  type(eddy_list_t), intent(inout) :: this
  integer,              intent(in) :: eddy_idx
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "add_new_eddy_to_list"
  !
continue
  !
  ! First, make sure that eddy_list is allocated.
  ! If it isnt, allocate it and reset num_eddies to zero.
  !
  if (.not. allocated(this%eddy_list)) then
    ! Allocate eddy_list to an initial array of 20 zeros
    this%eddy_list = init_array(20,0) ! F2003 AUTO-REALLOC
    ! Initialize num_eddies to zero
    this%num_eddies = 0
  end if
  !
  n = this%num_eddies ! current number of new eddies for this SEM box
  !
  ! Make sure that the list used to store the indices for new eddies
  ! is large enough to accommodate one more new eddy before adding
  ! this eddy to that list. If it isnt large enough, increase the size
  ! of the array for an additional 20 eddies.
  !
  if (n+1 > size(this%eddy_list)) then
    ! Use F2003 auto-reallocation to increase the size of eddy_list
    ! by appending 20 more zeros to the end of the array
    ! F2003 AUTO-REALLOC
    this%eddy_list = [ this%eddy_list , init_array(20,0) ]
  end if
  !
  ! Add the current eddy to the list of new eddies for this SEM box
  !
  this%eddy_list(n+1) = eddy_idx
  this%num_eddies = n+1
  !
end subroutine add_new_eddy_to_list
!
!###############################################################################
!
subroutine initialize_sem()
  !
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  ! Synthetic Eddy Method initialization
  !
  ! Based on
  ! Jarrin et al - Journal of Heat and Fluid Flow
  !
  ! Creates unsteady velocity fluctuations in a nozzle boundary layer
  ! based on a supplied turbulent kinetic energy profile.
  !
  ! Velocity and Reynolds stress profiles are read in.
  ! Convective velocity and turbulent length scale are input
  ! as variables in the bc specification.
  !
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  !.. Use Statements ..
  use ovar,       only : bc_input,bc_in
  use ovar,       only : average_sem_inflow_plane,num_timesteps
  use geovar,     only : bc_grps
  use flowvar,    only : facexyz_array
  use random_mod, only : set_random_seed
  !
  !.. Local Scalars ..
  integer  :: ierr,k,n,np,n1,n2
  integer  :: sem_idx,grp_idx,bc_idx,num_eddies
  integer  :: num_restart_sem_boxes
  integer  :: bc_input_idx
  real(wp) :: lturb,uconv,x,y,z,r,t
  !
  character(len=50)  :: array_name
  character(len=100) :: fname
  !
  !.. Local Allocatable Arrays ..
  type(eddy_list_t), allocatable :: new_eddies(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_sem"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize the random number generator
  !
#ifdef TROUBLESHOOT_PROFILE_INTERPOLATION
  call set_random_seed(seed=1)
#else
  call set_random_seed
#endif
  !
  ! Count the number of BC groups on this processor that are of the SEM type
  !
  num_restart_sem_boxes = 0
  num_sem_boxes = count( bc_in(:)%bc_type == bc_sem )
  !
  if (num_sem_boxes < 1) then
    call debug_timer(leaving_procedure,pname)
    return
  end if
  !
  ! Allocate the derived-type sem_box if needed.
  !
  if (allocated(sem_box)) then
    !
    ! sem_box has already been allocated which means it was read in from the
    ! restart file. Check to make sure the number of SEM boxes in the restart
    ! file matches the number of SEM boxes specified in the input file
    ! boundary conditions.
    !
    num_restart_sem_boxes = size(sem_box)
    !
    if (num_sem_boxes /= num_restart_sem_boxes) then
      write (error_message,1) num_sem_boxes,num_restart_sem_boxes
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
    !
  else
    !
    allocate ( sem_box(1:num_sem_boxes) )
    !
  end if
  !
  ! Allocate a derived-type sem_pts_t for each of the SEM BC groups
  !
  allocate ( my_sem_pts(1:num_sem_boxes) )
  !
  ! Find the SEM boundary conditions and use the BC specifications to
  ! begin populating the my_sem_pts and sem_box derived-types.
  !
  ! Create all the SEM boxes.
  ! NOTE: All MPI processes will contain all SEM boxes and their respective
  !       eddies, and this information should be identical across all MPI
  !       processes.
  !
  sem_idx = 0
  do bc_idx = 1,ubound(bc_in,dim=1)
    !
    if (bc_in(bc_idx)%bc_type == bc_sem) then
      !
      sem_idx = sem_idx + 1
      !
      lturb = bc_in(bc_idx)%sem%length_scale
      uconv = bc_in(bc_idx)%sem%convective_velocity
      fname = trim(adjustl(bc_in(bc_idx)%sem%profile_file))
      !
      bc_input_idx = bc_in(bc_idx)%bc_input_idx
      !
      sem_box(sem_idx)%bc_idx = bc_idx
      sem_box(sem_idx)%length_scale = lturb
      sem_box(sem_idx)%convective_velocity = uconv
      !
      sem_box(sem_idx)%alpha_aoa = bc_input( bc_input_idx )%alpha_aoa
      sem_box(sem_idx)%beta_aoa  = bc_input( bc_input_idx )%beta_aoa
      !
      if (mypnum == glb_root) then
        n = merge(1,0,bc_in(bc_idx)%sem%is_cylindrical)
      end if
      !
      call mpi_bcast(n,1_int_mpi,mpi_inttyp,glb_root,MPI_COMM_WORLD,mpierr)
      !
      sem_box(sem_idx)%is_cylindrical = (n == 1)
      !
      if (.not. allocated(sem_box(sem_idx)%profile)) then
        allocate ( sem_box(sem_idx)%profile )
      end if
      sem_box(sem_idx)%profile%fname = trim(adjustl(fname))
      !
    end if
    !
  end do
  !
  ! Create my_sem_pts on each MPI process.
  ! NOTE: Each element of this derived-type array represents a single SEM box,
  !       and the information within each element is localized and unique to
  !       each MPI process. If the current MPI process does not have any
  !       boundary faces for a given SEM box, the array element respective to
  !       that SEM box will be mostly empty.
  !
  create_my_sem_pts: &
  do grp_idx = 1,size(bc_grps)
    !
    if (bc_grps(grp_idx)%bc_type == bc_sem) then
      !
      bc_idx = bc_grps(grp_idx)%bc_idx
      !
      find_sem_idx: &
      do n = 1,num_sem_boxes
        if (sem_box(n)%bc_idx == bc_idx) then
          sem_idx = n
          exit find_sem_idx
        end if
      end do find_sem_idx
      !
     !write (iout,90) grp_idx,bc_idx,sem_idx
     !write (iout,91) trim(adjustl(fname))
     !90 format (/,2x,"grp_idx = ",i0,/, &
     !             2x," bc_idx = ",i0,/, &
     !             2x,"sem_idx = ",i0)
     !91 format (2x,"fname = '",a,"'",/)
      !
      ! Save this group index to the sem_box array
      !
      sem_box(sem_idx)%grp_idx = grp_idx
      !
      ! Save the beginning and ending SEM face point indices for
      ! the BC group to my_sem_pts for the respective SEM box
      !
      n1 = bc_grps(grp_idx)%beg_fp
      n2 = bc_grps(grp_idx)%end_fp
      np = n2-n1+1
      !
      my_sem_pts(sem_idx)%bc_idx = bc_idx ! Index for this BC
      my_sem_pts(sem_idx)%grp_idx = grp_idx ! Index for this BC group
      my_sem_pts(sem_idx)%num_fp = np      ! # of points in this SEM set
      !
      ! Allocate the SEM face points array
      !
      allocate ( my_sem_pts(sem_idx)%sem_fp(1:np) )
      !
      ! Save the face point indices for all the face points in this SEM set
      !
      my_sem_pts(sem_idx)%sem_fp(1:np)%idx = intseq(n1,n2)
      !
      ! Save the coordinates for each face point in this SEM set
      !
      do n = 1,my_sem_pts(sem_idx)%num_fp
        !
        k = my_sem_pts(sem_idx)%sem_fp(n)%idx
        !
        x = facexyz_array(1,k)
        y = facexyz_array(2,k)
        z = facexyz_array(3,k)
        !
        my_sem_pts(sem_idx)%sem_fp(n)%x = x
        my_sem_pts(sem_idx)%sem_fp(n)%y = y
        my_sem_pts(sem_idx)%sem_fp(n)%z = z
        !
#ifdef ENABLE_OUTPUT_TIME_HISTORY
        allocate ( my_sem_pts(sem_idx)% &
                             sem_fp(n)%time_history(1:nq+1,0:num_timesteps) , &
                             source=zero )
#endif
        !
      end do
      !
    end if
    !
  end do create_my_sem_pts
  !
  ! Initialize the SEM averaging arrays if needed
  !
  if (average_sem_inflow_plane) then
    call initialize_sem_avg()
  end if
  !
  ! For each SEM box:
  !   1) get the coordinates for the center of the box,
  !   2) get the inward facing normal for the inflow
  !      plane associated with the box,
  !   3) get the min/max extents of the box,
  !   4) compute the volume of the box, and
  !   5) compute the number of eddies to use for the box.
  !
  ! Get the coordinates for the center of each SEM box
  !
  call get_center_of_sem_boxes()
  !
  ! Get the polar/cylindrical coordinates for each SEM face point
  ! NOTE: This will also compute the polar/cylindrical coodinates for all the
  !       eddies if sem_box(:)%eddy is allocated (restart containing SEM).
  !
  call get_polar_coords_for_sem_fp()
  !
  ! Get the inward facing normal of each SEM box
  !
  call get_normal_of_sem_boxes()
  !
  ! Get the min/max extents of each SEM box
  !
  call get_size_of_sem_boxes()
  !
  ! If the eddy information was read in from a restart file, check that the
  ! location of the eddies are within the extents of their respective SEM box.
  !
  call check_restart_eddies()
  !
  ! Read in velocity profiles and Reynolds stresses
  !
  call read_sem_box_profiles()
  !
  ! Recompute the volume of the SEM boxes and the number of eddies for
  ! each box if this is not a restart simulation
  !
  call recompute_sem_box_volume()
  !
  ! Write out parameters of the SEM boxes
  !
  call output_details_of_sem_boxes()
  !
  ! Create new_eddies so that we can generate the new eddies and create
  ! the face point eddy-lists and the eddy face-point-lists.
  ! NOTE: This will use the same subroutine that updates these lists during
  !       the solver due to "new" eddies that replace "old" eddies that have
  !       convected outside the SEM box. To utilize this same subroutine, we
  !       basically need to make the same new_eddies array that is used
  !       during the solver, but here the new_eddies will contain all
  !       eddies across all SEM boxes since they are all really "new" eddies.
  !
  allocate ( new_eddies(1:num_sem_boxes) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"new_eddies",1,__LINE__,__FILE__,ierr,error_message)
  !
  do sem_idx = 1,num_sem_boxes
    !
    num_eddies = sem_box(sem_idx)%num_eddies
    !
    if (.not. allocated(sem_box(sem_idx)%eddy)) then
      !
      ! If sem_box(sem_idx)%eddy is not allocated, it needs to
      ! be allocated and initialized with new random eddies.
      !
      new_eddies(sem_idx)%num_eddies = num_eddies
      new_eddies(sem_idx)%eddy_list = intseq(1,num_eddies) ! F2003 AUTO-REALLOC
      !
    else
      !
      ! If sem_box(sem_idx)%eddy is already allocated, the eddy
      ! information was probably read in from a restart file. Setting
      ! new_eddies(sem_idx)%num_eddies to zero will prevent the loop within
      ! generate_new_eddies from altering the sem_box(sem_idx)%eddy array.
      !
      new_eddies(sem_idx)%num_eddies = 0
      !
    end if
    !
  end do
  !
  ! Generate initial distribution of eddies. If not initializing, eddy
  ! information for each SEM box is read in from the restart file and
  ! new_eddies(sem_idx)%num_eddies should be zero from the above loop,
  ! so the loop within generate_new_eddies will not alter the
  ! sem_box(sem_idx)%eddy(:) array.
  !
  ! NOTE: Since we are initializing eddies for the SEM, we want the
  !       x-location of the eddies to be randomly located throughout
  !       the SEM box. Therefore, the argument 'generate_at_xmin'
  !       needs to be false.
  !
  call generate_new_eddies(new_eddies,generate_at_xmin=fals,dt=zero)
  !
  ! Deallocate new_eddies since we no longer need it
  !
  deallocate ( new_eddies , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"new_eddies",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Now initialize the face point eddy-lists and the
  ! eddy face-point-lists using new_eddies
  !
  call initialize_eddy_lists()
  !
  ! Compute the distance to the nearest wall for each of the SEM face points
  !
  call compute_sem_fp_wall_distance()
  !
  ! Finally, interpolate the profile data to the location of each face point
  !
  do sem_idx = 1,num_sem_boxes
    !
    do n = 1,my_sem_pts(sem_idx)%num_fp
      !
      call interpolate_profile_data_to_fp(my_sem_pts(sem_idx)%sem_fp(n), &
                                          sem_box(sem_idx)%profile, &
                                          sem_idx,n)
      !
    end do
    !
  end do
  !
  call check_for_bad_number(__LINE__,__FILE__)
  !
  if (average_sem_inflow_plane) then
    call output_sem_averages(output_inflow_profile=true)
  end if
  !
 !call output_interpolated_profiles()
  !
  ! Now that we have interpolated all the profile data to the SEM
  ! face points, we no longer need the actual profile data so deallocate
  ! all the sem_box(:)%profile components to free up that memory.
  !
 !do sem_idx = 1,num_sem_boxes
 !  if (allocated(sem_box(sem_idx)%profile)) then
 !    deallocate ( sem_box(sem_idx)%profile , stat=ierr , errmsg=error_message )
 !    write (array_name,'("sem_box(",i0,")%profile")') sem_idx
 !    call alloc_error(pname,array_name,2,__LINE__,__FILE__,ierr,error_message)
 !  end if
 !end do
  !
  if (num_restart_sem_boxes > 0) then
    call dump_sem_data("read")
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format("ERROR!!! THE NUMBER OF SEM BOXES IN THE RESTART FILE ", &
           "DO NOT MATCH THE NUMBER OF SEM BOXES IN THE INPUT FILE! ", &
           "# of SEM boxes in the input file = ",i4, &
           "# of SEM boxes in the restart file = ",i4)
  !
end subroutine initialize_sem
!
!###############################################################################
!
subroutine output_interpolated_profiles()
  !
  !.. Local Scalars ..
  integer :: n,sem_idx,io,ierr
  !
  character(len=100) :: fname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_interpolated_profiles"
  !
continue
  !
  if (sum(my_sem_pts(:)%num_fp) > 0) then
    !
    write (fname,1) mypnum
    1 format ("sem_profile.",i4.4,".dat")
    !
    open (newunit=io,file=fname,status="replace",action="write", &
                     form="formatted",position="rewind", &
                     iostat=ierr,iomsg=error_message)
    call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
    !
    do sem_idx = 1,num_sem_boxes
      do n = 1,my_sem_pts(sem_idx)%num_fp
        call output_fp_interpolated_profiles(my_sem_pts(sem_idx)%sem_fp(n), &
                                             sem_idx,n,io)
      end do
    end do
    !
    close (io,iostat=ierr,iomsg=error_message)
    call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
    !
  end if
  !
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
end subroutine output_interpolated_profiles
!
!###############################################################################
!
subroutine output_fp_interpolated_profiles(this,sem_idx,fp_idx,io)
  !
  !.. Use Statements ..
  use ovar, only : grid_scaling_factor,rhoref,aref,tref
  !
  !.. Formal Arguments ..
  type(sem_fp_t), intent(in) :: this
  integer,        intent(in) :: sem_idx
  integer,        intent(in) :: fp_idx
  integer,        intent(in) :: io
  !
  !.. Local Scalars ..
  real(wp) :: y,z,yw,unx,uny,unz
  real(wp) :: ravg,uavg,vavg,wavg,tavg
  real(wp) :: r11,r21,r22,r31,r32,r33
  real(wp) :: a11,a21,a22,a31,a32,a33
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_fp_interpolated_profiles"
  !
continue
  !
  y = this%y
  z = this%z
  !
  yw = this%wall%distance
  !
  unx = this%wall%xnrm
  uny = this%wall%ynrm
  unz = this%wall%znrm
  !
  ravg = this%ravg*rhoref
  uavg = this%uavg*aref
  vavg = this%vavg*aref
  wavg = this%wavg*aref
  tavg = this%tavg*tref
  !
  r11 = this%r11*aref*aref
  r21 = this%r21*aref*aref
  r22 = this%r22*aref*aref
  r31 = this%r31*aref*aref
  r32 = this%r32*aref*aref
  r33 = this%r33*aref*aref
  !
  a11 = this%a11*aref
  a21 = this%a21*aref
  a22 = this%a22*aref
  a31 = this%a31*aref
  a32 = this%a32*aref
  a33 = this%a33*aref
  !
  write (io,1) sem_idx,fp_idx,this%idx
  !
  write (io,7) " Y",y," m"
  write (io,7) " Y",y/grid_scaling_factor," in"
  write (io,7) " Z",z," m"
  write (io,7) " Z",z/grid_scaling_factor," in",""
  !
  write (io,7) "  WALL DISTANCE",yw," m"
  write (io,7) "  WALL DISTANCE",yw/grid_scaling_factor," in",""
  !
  write (io,2) " UNX",unx
  write (io,2) " UNY",uny
  write (io,2) " UNZ",unz,""
  !
  write (io,3)
  write (io,2) "  RHO",ravg
  write (io,2) "    U",uavg
  write (io,2) "    V",vavg
  write (io,2) "    W",wavg
  write (io,2) "    T",tavg,""
  !
  write (io,4)
  write (io,2) "  <UU>",r11
  write (io,2) "  <VV>",r22
  write (io,2) "  <WW>",r33
  write (io,2) "  <UV>",r21
  write (io,2) "  <UW>",r31
  write (io,2) "  <VW>",r32,""
  !
  write (io,5)
  write (io,2) "  A11",a11
  write (io,2) "  A22",a22
  write (io,2) "  A33",a33
  if (abs(a21) < 1000.0_wp) then
    write (io,2) "  A21",a21
  else
    write (io,6) "  A21",a21
  end if
  if (abs(a31) < 1000.0_wp) then
    write (io,2) "  A31",a31
  else
    write (io,6) "  A31",a31
  end if
  if (abs(a32) < 1000.0_wp) then
    write (io,2) "  A32",a32
  else
    write (io,6) "  A32",a32
  end if
  write (io,*)
  write (io,2) "  A11            ",a11
  if (abs(a21+a22) < 1000.0_wp) then
    write (io,2) "  A21 + A22      ",a21+a22
  else
    write (io,6) "  A21 + A22      ",a21+a22
  end if
  if (abs(a31+a32+a33) < 1000.0_wp) then
    write (io,2) "  A31 + A32 + A33",a31+a32+a33,""
  else
    write (io,6) "  A31 + A32 + A33",a31+a32+a33,""
  end if
  !
  1 format (/,80("#"),/, &
            2x,"SEM Box       = ",i0,/, &
            2x,"SEM FP Index  = ",i0,/, &
            2x,"Grid FP Index = ",i0,/, &
            80("-"),/)
  2 format (2x,a," = ",f18.13,:,/,a)
  3 format (2x,"Interpolated Transformed Primitive Variables Profile")
  4 format (2x,"Interpolated Transformed Reynolds Stress Tensor Profile")
  5 format (2x,"Interpolated Transformed Cholesky Decomposition ", &
               "of Reynolds Stress Tensor Profile")
  6 format (2x,a," = ",es18.10,:,/,a)
  7 format (2x,a," = ",f18.13,a,:,/,a)
  !
end subroutine output_fp_interpolated_profiles
!
!###############################################################################
!
subroutine initialize_sem_avg()
  !
  !.. Local Scalars ..
  integer :: n,nfp,sem_idx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_sem_avg"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  allocate ( sem_avg(1:num_sem_boxes) )
  !
  do sem_idx = 1,num_sem_boxes
    !
    nfp = my_sem_pts(sem_idx)%num_fp
    !
    allocate ( sem_avg(sem_idx)%fp_avg(1:nfp) )
    !
    do n = 1,nfp
      !
      sem_avg(sem_idx)%fp_avg(n)%r   = zero
      sem_avg(sem_idx)%fp_avg(n)%ru  = zero
      sem_avg(sem_idx)%fp_avg(n)%rv  = zero
      sem_avg(sem_idx)%fp_avg(n)%rw  = zero
      sem_avg(sem_idx)%fp_avg(n)%e   = zero
      sem_avg(sem_idx)%fp_avg(n)%t   = zero
      sem_avg(sem_idx)%fp_avg(n)%u   = zero
      sem_avg(sem_idx)%fp_avg(n)%v   = zero
      sem_avg(sem_idx)%fp_avg(n)%w   = zero
      sem_avg(sem_idx)%fp_avg(n)%p   = zero
      sem_avg(sem_idx)%fp_avg(n)%uu  = zero
      sem_avg(sem_idx)%fp_avg(n)%vv  = zero
      sem_avg(sem_idx)%fp_avg(n)%ww  = zero
      sem_avg(sem_idx)%fp_avg(n)%uv  = zero
      sem_avg(sem_idx)%fp_avg(n)%uw  = zero
      sem_avg(sem_idx)%fp_avg(n)%vw  = zero
      sem_avg(sem_idx)%fp_avg(n)%ut  = zero
      sem_avg(sem_idx)%fp_avg(n)%vt  = zero
      sem_avg(sem_idx)%fp_avg(n)%wt  = zero
      sem_avg(sem_idx)%fp_avg(n)%rt  = zero
      sem_avg(sem_idx)%fp_avg(n)%ruu = zero
      sem_avg(sem_idx)%fp_avg(n)%rvv = zero
      sem_avg(sem_idx)%fp_avg(n)%rww = zero
      sem_avg(sem_idx)%fp_avg(n)%ruv = zero
      sem_avg(sem_idx)%fp_avg(n)%ruw = zero
      sem_avg(sem_idx)%fp_avg(n)%rvw = zero
      sem_avg(sem_idx)%fp_avg(n)%rut = zero
      sem_avg(sem_idx)%fp_avg(n)%rvt = zero
      sem_avg(sem_idx)%fp_avg(n)%rwt = zero
      !
    end do
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine initialize_sem_avg
!
!###############################################################################
!
subroutine compute_sem_fp_wall_distance()
  !
  !.. Local Scalars ..
  integer  :: i,ierr,n,sem_idx
  real(wp) :: c1,c2,y,y1,y2
  !
  !.. Local Allocatable Arrays ..
  real(wp),     allocatable :: xyz_wall(:,:)
  real(wp),     allocatable :: xyz_sem(:,:)
  type(wall_t), allocatable :: wall(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_sem_fp_wall_distance"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! First, collect the coordinates of all the
  ! wall face points across all processors.
  !
  call collect_wall_fp_coordinates(xyz_wall)
  !
  ! Create an array containing the coordinates
  ! of all SEM face points on this processor.
  !
  call collect_sem_fp_coordinates(xyz_sem)
  !
  ! Now, do a search for the minimum wall distance of each SEM face point
  !
  call brute_force_wall_distance(xyz_wall,xyz_sem,wall)
  !
  ! Deallocate xyz_wall and xyz_sem since they are no longer needed
  !
  deallocate ( xyz_wall , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_wall",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( xyz_sem , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_sem",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Finally, copy the computed wall distances into my_sem_pts
  !
  i = 0
  do sem_idx = 1,num_sem_boxes
    !
    do n = 1,my_sem_pts(sem_idx)%num_fp
      i = i + 1
      my_sem_pts(sem_idx)%sem_fp(n)%wall = wall(i)
    end do
    !
  end do
  !
  ! Deallocate wall since it is no longer needed
  !
  deallocate ( wall , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"wall",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_sem_fp_wall_distance
!
!###############################################################################
!
subroutine collect_wall_fp_coordinates(xyz_wall)
  !
  !.. Use Statements ..
  use geovar,      only : bc_grps
  use flowvar,     only : facexyz_array
  use metrics_mod, only : facenrm
  !
  !.. Formal Arguments ..
  real(wp), allocatable, intent(out) :: xyz_wall(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,ierr,n,np,n1,n2,grp_idx
  real(wp) :: rnx,rny,rnz,ds
  !
  !.. Local Arrays ..
  integer(INT_MPI) :: cpu_nums(1:ncpu)
  integer(INT_MPI) :: cpu_disp(1:ncpu)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "collect_wall_fp_coordinates"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  cpu_nums = 0_int_mpi
  cpu_disp = 0_int_mpi
  !
  ! Count the number of wall face points on the current processor
  !
  do grp_idx = 1,size(bc_grps)
    !
    if (any(bc_grps(grp_idx)%bc_type == bc_physical_walls)) then
      !
      n1 = bc_grps(grp_idx)%beg_fp
      n2 = bc_grps(grp_idx)%end_fp
      np = n2-n1+1
      !
      cpu_nums(mypnum+1) = cpu_nums(mypnum+1) + int(np,kind=int_mpi)
      !
    end if
    !
  end do
  !
  ! Gather the number of wall face points for each processor
  !
  if (ncpu > 1) then
    call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_DATATYPE_NULL, &
                       cpu_nums,1_int_mpi,MPI_INTEGER, &
                       MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Allocate the xyz_wall array
  !
  np = sum(cpu_nums)
  !
  allocate ( xyz_wall(1:6,1:np) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_wall",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Get the offsets for gathering the face point coordinates
  !
  cpu_disp = 0_int_mpi
  do n = 2,ncpu
    cpu_disp(n) = cpu_disp(n-1) + cpu_nums(n-1)
  end do
  !
  ! Fill in the section of xyz_wall for the current processor
  ! with the coordinates of the wall face points.
  !
  i = cpu_disp(mypnum+1) ! beginning offset for this processor
  !
  do grp_idx = 1,size(bc_grps)
    !
    if (any(bc_grps(grp_idx)%bc_type == bc_physical_walls)) then
      !
      n1 = bc_grps(grp_idx)%beg_fp
      n2 = bc_grps(grp_idx)%end_fp
      !
      do n = n1,n2
        !
        i = i + 1
        !
        ! Copy the coordinates of this wall face point into xyz_wall.
        !
        xyz_wall(1:3,i) = facexyz_array(1:3,n)
        !
        ! Copy the interior-pointing unit normal to the wall into xyz_wall.
        ! However, facenrm contains the exterior-pointing area normal to the
        ! wall, so we need to multiply each normal component by -1 and convert
        ! to a unit normal.
        !
        rnx = -facenrm(1,n)
        rny = -facenrm(2,n)
        rnz = -facenrm(3,n)
        !
        ds = one/sqrt(rnx*rnx + rny*rny + rnz*rnz)
        !
        xyz_wall(4,i) = rnx*ds
        xyz_wall(5,i) = rny*ds
        xyz_wall(6,i) = rnz*ds
        !
      end do
      !
    end if
    !
  end do
  !
  ! Adjust the displacement and cpu_nums because
  ! each face point will have 3 real values
  !
  cpu_nums(:) = size(xyz_wall,dim=1,kind=int_mpi) * cpu_nums(:)
  cpu_disp(:) = size(xyz_wall,dim=1,kind=int_mpi) * cpu_disp(:)
  !
  ! Gather the face point coordinates from all processors
  !
  if (ncpu > 1) then
    call mpi_allgatherv(MPI_IN_PLACE,0_int_mpi,MPI_DATATYPE_NULL, &
                        xyz_wall,cpu_nums,cpu_disp,mpi_flttyp, &
                        MPI_COMM_WORLD,mpierr)
  end if
  !
#ifdef DEBUG_WALL_NORMALS
  if (mypnum == glb_root) then
    do n = 1,size(xyz_wall,dim=2)
      write (200,10) xyz_wall(4,n),xyz_wall(5,n),xyz_wall(6,n)
    end do
  end if
  10 format (3f15.10)
  !
  call debug_normals
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine collect_wall_fp_coordinates
!
!###############################################################################
!
subroutine debug_normals()
  !
  !.. Use Statements ..
  use geovar, only : nfbnd,ncell,bface,face,grid
  use flowvar, only : facexyz_array
  !
  !.. Local Scalars ..
  integer  :: i,n,nc,nf
  real(wp) :: x,y,z,r,t
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: pts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "debug_normals"
  !
continue
  !
  do nf = 1,nfbnd
    !
    if (all(bface(1,nf) /= bc_walls)) cycle
    !
    nc = ncell + nf
    !
    pts = grid%elem(nc)%pts(:) ! F2003 AUTO-REALLOC
    !
    write (1000,1) nf
    write (1001,1) nf
    do n = 1,size(pts)
      i = pts(n)
      x = grid%xyz(1,i)
      y = grid%xyz(2,i)
      z = grid%xyz(3,i)
      r = hypot(y,z)
      if (r /= zero) then
        t = atan2(z,y)
        if (t < zero) t = t + two*pi
        t = t * radians2degrees
      else
        t = zero
      end if
      write (1000,2) i,x,y,z,r,t
    end do
    do i = face(nf)%beg_fp,face(nf)%end_fp
      x = facexyz_array(1,i)
      y = facexyz_array(2,i)
      z = facexyz_array(3,i)
      r = hypot(y,z)
      if (r /= zero) then
        t = atan2(z,y)
        if (t < zero) t = t + two*pi
        t = t * radians2degrees
      else
        t = zero
      end if
      write (1001,2) i,x,y,z,r,t
    end do
    !
  end do
  !
  1 format (/,2x,"BOUNDARY FACE ",i0)
  2 format (4x,i6,4f15.10,f19.10)
  !
end subroutine debug_normals
!
!###############################################################################
!
subroutine collect_sem_fp_coordinates(xyz_sem)
  !
  !.. Formal Arguments ..
  real(wp), allocatable, intent(out) :: xyz_sem(:,:)
  !
  !.. Local Scalars ..
  integer :: i,ierr,n,np,sem_idx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "collect_sem_fp_coordinates"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  np = sum(my_sem_pts(:)%num_fp)
  !
  allocate ( xyz_sem(1:3,1:np) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_sem",1,__LINE__,__FILE__,ierr,error_message)
  !
  i = 0
  do sem_idx = 1,num_sem_boxes
    !
    do n = 1,my_sem_pts(sem_idx)%num_fp
      !
      i = i + 1
      !
      xyz_sem(1,i) = my_sem_pts(sem_idx)%sem_fp(n)%x
      xyz_sem(2,i) = my_sem_pts(sem_idx)%sem_fp(n)%y
      xyz_sem(3,i) = my_sem_pts(sem_idx)%sem_fp(n)%z
      !
    end do
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine collect_sem_fp_coordinates
!
!###############################################################################
!
subroutine brute_force_wall_distance(xyz_wall,xyz_pts,wall)
  !
  !.. Formal Arguments ..
  real(wp),                   intent(in) :: xyz_wall(:,:)
  real(wp),                   intent(in) :: xyz_pts(:,:)
  type(wall_t), allocatable, intent(out) :: wall(:)
  !
  !.. Local Scalars ..
  integer  :: i,ierr,n,nfp,nwp
  real(wp) :: xf,yf,zf,dx,dy,dz,dw2,dw2_min
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "brute_force_wall_distance"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nwp = size(xyz_wall,dim=2)
  nfp = size(xyz_pts,dim=2)
  !
  allocate ( wall(1:nfp) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"wall",1,__LINE__,__FILE__,ierr,error_message)
  !
 !!$OMP PARALLEL DO &
 !!$OMP&         DEFAULT(NONE) &
 !!$OMP&         SHARED(nfp,nwp,xyz_pts,xyz_wall,wall) &
 !!$OMP&         PRIVATE(n,i,xf,yf,zf,dx,dy,dz,dw2,dw2_min)
  !
  do n = 1,nfp
    !
    xf = xyz_pts(1,n)
    yf = xyz_pts(2,n)
    zf = xyz_pts(3,n)
    !
    dw2_min = huge(dw2_min)
    !
    do i = 1,nwp
      !
      dx = xf - xyz_wall(1,i)
      dy = yf - xyz_wall(2,i)
      dz = zf - xyz_wall(3,i)
      !
      dw2 = dx*dx + dy*dy + dz*dz
      !
      if (dw2 < dw2_min) then
        !
        wall(n)%x = xyz_wall(1,i)
        wall(n)%y = xyz_wall(2,i)
        wall(n)%z = xyz_wall(3,i)
        !
        if (size(xyz_wall,dim=1) >= 6) then
          wall(n)%xnrm = xyz_wall(4,i)
          wall(n)%ynrm = xyz_wall(5,i)
          wall(n)%znrm = xyz_wall(6,i)
        end if
        !
        dw2_min = dw2
        !
      end if
      !
    end do
    !
    wall(n)%distance = sqrt( dw2_min )
    !
    if (ncpu == 1) then
      if (mypnum == glb_root) then
        if (mod(n,100) == 0) then
          write (iout,1) n,nfp,100.0_wp*real(n,kind=wp)/real(nfp,kind=wp)
        end if
      end if
    end if
    !
  end do
  1 format (2x,"Finished ",i8," of ",i8," (",f9.4," %)")
  !
 !!$OMP END PARALLEL DO
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine brute_force_wall_distance
!
!###############################################################################
!
#ifdef TROUBLESHOOT_PROFILE_INTERPOLATION
#define EXTREME_TROUBLESHOOT_PROFILE_INTERPOLATION
#endif
#ifndef EXTREME_TROUBLESHOOT_PROFILE_INTERPOLATION
pure &
#endif
subroutine interpolate_profile_data_to_fp(this,prof,sem_idx,fp_idx)
  !
#ifdef EXTREME_TROUBLESHOOT_PROFILE_INTERPOLATION
  !.. Use Statements ..
  use ovar, only : aref
#endif
  !
  !.. Formal Arguments ..
  type(sem_fp_t),   intent(inout) :: this
  type(sem_profile_t), intent(in) :: prof
  integer,             intent(in) :: sem_idx
  integer,             intent(in) :: fp_idx
  !
  !.. Local Scalars ..
  integer  :: n,n1,n2,np
  real(wp) :: unx,uny,unz,unx2,uny2,unz2
  real(wp) :: c1,c2,c1i,c2i
  real(wp) :: uu,vv,ww,y,y1,y2
  real(wp) :: s11,s22,s33,s21,s31,s32
  real(wp) :: v11,v12,v13
  real(wp) :: v21,v22,v23
  real(wp) :: v31,v32,v33
  real(wp) :: t11,t12,t13,t14,t15,t16
  real(wp) :: t21,t22,t23,t24,t25,t26
  real(wp) :: t31,t32,t33,t34,t35,t36
  real(wp) :: t41,t42,t43,t44,t45,t46
  real(wp) :: t51,t52,t53,t54,t55,t56
  real(wp) :: t61,t62,t63,t64,t65,t66
  !
#ifdef EXTREME_TROUBLESHOOT_PROFILE_INTERPOLATION
  real(wp) :: radius,theta,aref2
  integer  :: io,ierr
  logical(lk) :: output_profile
  integer, save :: nti = 0 ! top    (   0 < theta <  11.25 ) inside
  integer, save :: nto = 0 ! top    (   0 < theta <  11.25 ) outside
  integer, save :: nbi = 0 ! bottom ( 180 < theta < 191.25 ) inside
  integer, save :: nbo = 0 ! bottom ( 180 < theta < 191.25 ) outside
  integer, save :: ni1 = 0 ! sem_idx == 1
  integer, save :: ni2 = 0 ! sem_idx == 2
#endif
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: r(:),u(:),v(:),w(:),t(:)
  real(wp), allocatable :: r11(:),r22(:),r33(:)
  real(wp), allocatable :: r21(:),r31(:),r32(:)
  real(wp), allocatable :: a11(:),a22(:),a33(:)
  real(wp), allocatable :: a21(:),a31(:),a32(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interpolate_profile_data_to_fp"
  !
continue
  !
  y = this%wall%distance
  !
#ifdef EXTREME_TROUBLESHOOT_PROFILE_INTERPOLATION
  call yz_to_rt(this%y,this%z,radius,theta, &
                yc=sem_box(sem_idx)%center%y, &
                zc=sem_box(sem_idx)%center%z, &
                return_degrees=true)
  !
  io = 0
  if (sem_idx == 1) then
    ni1 = ni1 + 1
    io = 10000 + ni1
  else if (sem_idx == 2) then
    ni2 = ni2 + 1
    io = 20000 + ni2
  end if
 !if (theta < 11.25_wp) then
 !  if (sem_idx == 1) then
 !    nti = nti + 1
 !    io = 600 + nti
 !  else if (sem_idx == 2) then
 !    nto = nto + 1
 !    io = 800 + nto
 !  end if
 !else if (theta > 180.0_wp .and. theta < 191.25_wp) then
 !  if (sem_idx == 1) then
 !    nbi = nbi + 1
 !    io = 700 + nbi
 !  else if (sem_idx == 2) then
 !    nbo = nbo + 1
 !    io = 900 + nbo
 !  end if
 !end if
  output_profile = (io /= 0)
  !
  aref2 = aref*aref
#endif
  !
  unx = this%wall%xnrm
  uny = this%wall%ynrm
  unz = this%wall%znrm
  !
  unx2 = unx*unx
  uny2 = uny*uny
  unz2 = unz*unz
  !
  c2 = uny*uny + unz*unz
  c1 = sqrt(c2)
  !
  c2i = one / c2 ! normal should never be only in unx so this should be fine
  c1i = one / c1
  !
  ! Create the coefficients for the vector transformation
  !
  v11 = c1           ; v12 = unx ; v13 =  zero
  v21 = -unx*uny*c1i ; v22 = uny ; v23 = -unz*c1i
  v31 = -unx*unz*c1i ; v32 = unz ; v33 =  uny*c1i
  !
  ! Create the coefficients for the tensor transformation
  !
  t11 = c2
  t21 = unx2*uny2*c2i
  t31 = unx2*unz2*c2i
  t41 = -unx*uny
  t51 = -unx*unz
  t61 = unx2*uny*unz*c2i
  !
  t12 = unx2
  t22 = uny2
  t32 = unz2
  t42 = unx*uny
  t52 = unx*unz
  t62 = uny*unz
  !
  t13 = zero
  t23 = unz2*c2i
  t33 = uny2*c2i
  t43 = zero
  t53 = zero
  t63 = -uny*unz*c2i
  !
  t14 = two*c1*unx
  t24 = -two*unx*uny2*c1i
  t34 = -two*unx*unz2*c1i
  t44 = c1*uny - unx2*uny*c1i
  t54 = c1*unz - unx2*unz*c1i
  t64 = -two*unx*uny*unz*c1i
  !
  t15 = zero
  t25 =  two*unx*uny*unz*c2i
  t35 = -two*unx*uny*unz*c2i
  t45 = -unz
  t55 =  uny
 !t65 = unx*(unz2 - uny2)*c2i
  t65 = -unx*uny2*c2i + unx*unz2*c2i
  !
  t16 = zero
  t26 = -two*uny*unz*c1i
  t36 =  two*uny*unz*c1i
  t46 = -unx*unz*c1i
  t56 =  unx*uny*c1i
 !t66 = (uny2 - unz2)*c1i
  t66 = uny2*c1i - unz2*c1i
  !
#ifdef EXTREME_TROUBLESHOOT_PROFILE_INTERPOLATION
  if (output_profile) then
    !
    write (io,1) "    y",this%y
    write (io,1) "    z",this%z
    write (io,4) "    r",radius,ten2*(radius-0.0125_wp)/(0.045_wp-0.0125_wp)
   !write (io,1) "    r",radius
    write (io,1) "theta",theta,""
    write (io,1) "wall_distance",y,""
    write (io,1) "unx",unx
    write (io,1) "uny",uny
    write (io,1) "unz",unz,""
    write (io,1) "unx2",unx2
    write (io,1) "uny2",uny2
    write (io,1) "unz2",unz2,""
    write (io,1) "c1",c1
    write (io,1) "c2",c2
    write (io,1) "c1i",c1i
    write (io,1) "c2i",c2i,""
    !
    write (io,2) v11,v12,v13,v21,v22,v23,v31,v32,v33
    !
    write (io,3) t11,t12,t13,t14,t15,t16, &
                 t21,t22,t23,t24,t25,t26, &
                 t31,t32,t33,t34,t35,t36, &
                 t41,t42,t43,t44,t45,t46, &
                 t51,t52,t53,t54,t55,t56, &
                 t61,t62,t63,t64,t65,t66
    !
  end if
  1 format (2x,a," = ",f18.13,:,/,a)
  2 format (2x,"Vector transformation matrix",/,3(4x,3f18.13,/))
  3 format (2x,"Tensor transformation matrix",/,6(4x,6f18.13,/))
  4 format (2x,a," = ",f18.13," ( ",f7.3," %)",:,/,a)
#endif
  !
  ! First, initialize all the local allocatable arrays
  !
  np = prof%npts
  !
  r = prof%ravg(:)    ! F2003 AUTO-REALLOC
  u = prof%uavg(:)    ! F2003 AUTO-REALLOC
  v = prof%vavg(:)    ! F2003 AUTO-REALLOC
  w = prof%wavg(:)    ! F2003 AUTO-REALLOC
  t = prof%tavg(:)    ! F2003 AUTO-REALLOC
  !
  r11 = prof%tau11(:) ! F2003 AUTO-REALLOC
  r22 = prof%tau22(:) ! F2003 AUTO-REALLOC
  r33 = prof%tau33(:) ! F2003 AUTO-REALLOC
  r21 = prof%tau21(:) ! F2003 AUTO-REALLOC
  r31 = prof%tau31(:) ! F2003 AUTO-REALLOC
  r32 = prof%tau32(:) ! F2003 AUTO-REALLOC
  !
  a11 = prof%a11(:)   ! F2003 AUTO-REALLOC
  a22 = prof%a22(:)   ! F2003 AUTO-REALLOC
  a33 = prof%a33(:)   ! F2003 AUTO-REALLOC
  a21 = prof%a21(:)   ! F2003 AUTO-REALLOC
  a31 = prof%a31(:)   ! F2003 AUTO-REALLOC
  a32 = prof%a32(:)   ! F2003 AUTO-REALLOC
  !
  ! Use the interior-facing wall unit-normal to convert the profile velocity
  ! vector and stress tensor from wall-aligned coordinates to Cartesian
  ! coordinates.
  !
  do n = 1,np
    !
    uu = prof%uavg(n)
    vv = prof%vavg(n)
    ww = prof%wavg(n)
    !
    u(n) = uu*v11 + vv*v12 + ww*v13
    v(n) = uu*v21 + vv*v22 + ww*v23
    w(n) = uu*v31 + vv*v32 + ww*v33
    !
    s11 = prof%tau11(n)
    s22 = prof%tau22(n)
    s33 = prof%tau33(n)
    s21 = prof%tau21(n)
    s31 = prof%tau31(n)
    s32 = prof%tau32(n)
    !
    r11(n) = t11*s11 + t12*s22 + t13*s33 + t14*s21 + t15*s31 + t16*s32
    r22(n) = t21*s11 + t22*s22 + t23*s33 + t24*s21 + t25*s31 + t26*s32
    r33(n) = t31*s11 + t32*s22 + t33*s33 + t34*s21 + t35*s31 + t36*s32
    r21(n) = t41*s11 + t42*s22 + t43*s33 + t44*s21 + t45*s31 + t46*s32
    r31(n) = t51*s11 + t52*s22 + t53*s33 + t54*s21 + t55*s31 + t56*s32
    r32(n) = t61*s11 + t62*s22 + t63*s33 + t64*s21 + t65*s31 + t66*s32
    !
    a11(n) = sqrt( max( zero , r11(n) ) )
    !
    if (abs(a11(n)) < rndoff) then
      a11(n) = zero
      a21(n) = zero
      a31(n) = zero
    else
      a21(n) = r21(n) / (a11(n)+rndoff)
      a31(n) = r31(n) / (a11(n)+rndoff)
    end if
    !
    a22(n) = sqrt( max( zero , r22(n) - (a21(n)**2) ) )
    !
    if (abs(a22(n)) < rndoff) then
      a22(n) = zero
      a32(n) = zero
    else
      a32(n) = (r32(n) - a21(n)*a31(n)) / (a22(n)+rndoff)
    end if
    !
    a33(n) = sqrt( max( zero , r33(n) - (a31(n)**2) -(a32(n)**2) ) )
    !
    if (abs(a33(n)) < rndoff) a33(n) = zero
    !
  end do
  !
  ! Find the indices of the profile data that
  ! bound the wall distance of this face point.
  !
  n1 = 1
  n2 = 2
  do n = 2,np
    if (prof%y(n) > y) exit
    n1 = n
    n2 = n+1
  end do
  !
  ! Compute the profile data at this face point depending on whether or not
  ! the index of the upper bound is greater than the sizes of the profile data.
  !
  if (n2 > np) then
    !
    ! The index of the upper bound IS GREATER THAN the size of
    ! the profile data, so just set the profile data for the
    ! face point equal to the profile at the last data point.
    !
    n1 = np
    n2 = np
    !
    y1 = prof%y(n1)
    y2 = prof%y(n2)
    !
    c1 = zero
    c2 = one
    !
  else
    !
    ! The index of the upper bound IS LESS THAN OR EQUAL to the
    ! size of the profile data, so perform a linear interpolation
    ! of the bounding profile data to the face point.
    !
    y1 = prof%y(n1)
    y2 = prof%y(n2)
    !
    c1 = (y-y2)/(y2-y1)
    c2 = (y-y1)/(y2-y1)
    !
  end if
  !
  this%ravg = c2*r(n2) - c1*r(n1)
  this%uavg = c2*u(n2) - c1*u(n1)
  this%vavg = c2*v(n2) - c1*v(n1)
  this%wavg = c2*w(n2) - c1*w(n1)
  this%tavg = c2*t(n2) - c1*t(n1)
  !
  this%r11 = c2*r11(n2) - c1*r11(n1)
  this%r21 = c2*r21(n2) - c1*r21(n1)
  this%r22 = c2*r22(n2) - c1*r22(n1)
  this%r31 = c2*r31(n2) - c1*r31(n1)
  this%r32 = c2*r32(n2) - c1*r32(n1)
  this%r33 = c2*r33(n2) - c1*r33(n1)
  !
  this%a11 = c2*a11(n2) - c1*a11(n1)
  this%a21 = c2*a21(n2) - c1*a21(n1)
  this%a22 = c2*a22(n2) - c1*a22(n1)
  this%a31 = c2*a31(n2) - c1*a31(n1)
  this%a32 = c2*a32(n2) - c1*a32(n1)
  this%a33 = c2*a33(n2) - c1*a33(n1)
  !
#ifdef EXTREME_TROUBLESHOOT_PROFILE_INTERPOLATION
  if (output_profile) then
    !
    write (io,*)
    write (io,30) "n1",n1
    write (io,30) "n2",n2,""
    !
    write (io,31) "y1",y1
    write (io,31) "y2",y2,""
    !
    write (io,31) "c1",c1
    write (io,31) "c2",c2,""
    !
    write (io,32)
    write (io,31) "  U",this%uavg*aref
    write (io,31) "  V",this%vavg*aref
    write (io,31) "  W",this%wavg*aref,""
    !
    write (io,33)
    write (io,31) "  <UU>",this%r11*aref2
    write (io,31) "  <VV>",this%r22*aref2
    write (io,31) "  <WW>",this%r33*aref2
    write (io,31) "  <UV>",this%r21*aref2
    write (io,31) "  <UW>",this%r31*aref2
    write (io,31) "  <VW>",this%r32*aref2,""
    !
    write (io,34)
    write (io,31) "  A11",this%a11*aref
    write (io,31) "  A22",this%a22*aref
    write (io,31) "  A33",this%a33*aref
    write (io,31) "  A21",this%a21*aref
    write (io,31) "  A31",this%a31*aref
    write (io,31) "  A32",this%a32*aref
    write (io,*)
    write (io,31) "  A11            ",(this%a11)*aref
    write (io,31) "  A21 + A22      ",(this%a21+this%a22)*aref
    write (io,31) "  A31 + A32 + A33",(this%a31+this%a32+this%a33)*aref,""
    !
!#define SUPER_EXTREME_TROUBLESHOOT_PROFILE_INTERPOLATION
#ifdef SUPER_EXTREME_TROUBLESHOOT_PROFILE_INTERPOLATION
    write (io,11) "Original"
    do n = 1,np
      write (io,21) n,prof%uavg(n)*aref,prof%vavg(n)*aref,prof%wavg(n)*aref
    end do
    !
    write (io,12) "Original"
    do n = 1,np
      write (io,21) n,prof%tau11(n)*aref2, &
                      prof%tau22(n)*aref2, &
                      prof%tau33(n)*aref2, &
                      prof%tau21(n)*aref2, &
                      prof%tau31(n)*aref2, &
                      prof%tau32(n)*aref2
    end do
    !
    write (io,13) "Original"
    do n = 1,np
      write (io,21) n,prof%a11(n)*aref,prof%a22(n)*aref,prof%a33(n)*aref, &
                      prof%a21(n)*aref,prof%a31(n)*aref,prof%a32(n)*aref
    end do
    !
    write (io,11) "Transformed"
    do n = 1,np
      write (io,21) n,u(n)*aref,v(n)*aref,w(n)*aref
    end do
    !
    write (io,12) "Transformed"
    do n = 1,np
      write (io,21) n,r11(n)*aref2,r22(n)*aref2,r33(n)*aref2, &
                      r21(n)*aref2,r31(n)*aref2,r32(n)*aref2
    end do
    !
    write (io,13) "Transformed"
    do n = 1,np
      write (io,21) n,a11(n)*aref,a22(n)*aref,a33(n)*aref, &
                      a21(n)*aref,a31(n)*aref,a32(n)*aref
    end do
#endif
    !
  end if
  !
  11 format (/, &
             2x,a," Velocity Profile",/, &
             11x,"X-Velocity",3x,"Y-Velocity",3x,"Z-Velocity")
  12 format (/, &
             2x,a," Reynolds Stress Tensor Profile",/, &
             14x,"<UU>",9x,"<VV>",9x,"<WW>",9x,"<UV>",9x,"<UW>",9x,"<VW>")
  13 format (/, &
             2x,a," Cholesky Decomposition of", &
                  " Reynolds Stress Tensor Profile",/, &
             15x,"A11",10x,"A22",10x,"A33",10x,"A21",10x,"A31",10x,"A32")
  !
  21 format (3x,i4,2x,6f13.8)
  !
  30 format (2x,a," = ",i0,:,/,a)
  31 format (2x,a," = ",f18.13,:,/,a)
  32 format (2x,"Interpolated Transformed Velocity Profile")
  33 format (2x,"Interpolated Transformed Reynolds Stress Tensor Profile")
  34 format (2x,"Interpolated Transformed Cholesky Decomposition ", &
                "of Reynolds Stress Tensor Profile")
  !
  close (unit=io,iostat=ierr)
#endif
  !
end subroutine interpolate_profile_data_to_fp
!
!###############################################################################
!
subroutine get_center_of_sem_boxes()
  !
  !.. Local Scalars ..
  integer  :: n,num_fp,sem_idx,ierr
  real(wp) :: val
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: box_num(:)
  real(wp), allocatable :: box_sum(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_center_of_sem_boxes"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  box_num = init_array(num_sem_boxes,zero)   ! F2003 AUTO-REALLOC
  box_sum = init_array(3,num_sem_boxes,zero) ! F2003 AUTO-REALLOC
  !
  ! Get the average contribution to each SEM box from this MPI process
  !
  do sem_idx = 1,num_sem_boxes
    !
    num_fp = my_sem_pts(sem_idx)%num_fp
    !
    box_num(sem_idx) = num_fp
    !
    do n = 1,num_fp
      !
      box_sum(1,sem_idx) = box_sum(1,sem_idx) + my_sem_pts(sem_idx)%sem_fp(n)%x
      box_sum(2,sem_idx) = box_sum(2,sem_idx) + my_sem_pts(sem_idx)%sem_fp(n)%y
      box_sum(3,sem_idx) = box_sum(3,sem_idx) + my_sem_pts(sem_idx)%sem_fp(n)%z
      !
    end do
   !box_sum(1,sem_idx) = sum( my_sem_pts(sem_idx)%sem_fp(1:num_fp)%x )
   !box_sum(2,sem_idx) = sum( my_sem_pts(sem_idx)%sem_fp(1:num_fp)%y )
   !box_sum(3,sem_idx) = sum( my_sem_pts(sem_idx)%sem_fp(1:num_fp)%z )
    !
  end do
  !
  ! Sum the average contributions across all MPI processes
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,box_num,size(box_num,kind=int_mpi), &
                       mpi_inttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
    call mpi_allreduce(MPI_IN_PLACE,box_sum,size(box_sum,kind=int_mpi), &
                       mpi_flttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Finish computing the average and save to the SEM box
  !
  do sem_idx = 1,num_sem_boxes
    !
    val = one / real(box_num(sem_idx),kind=wp)
    !
    sem_box(sem_idx)%center%x = box_sum(1,sem_idx) * val
    sem_box(sem_idx)%center%y = box_sum(2,sem_idx) * val
    sem_box(sem_idx)%center%z = box_sum(3,sem_idx) * val
    !
  end do
  !
  ! Deallocate box_num and box_sum before we leave
  !
  if (allocated(box_num)) then
    deallocate ( box_num , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"box_num",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(box_sum)) then
    deallocate ( box_sum , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"box_sum",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_center_of_sem_boxes
!
!###############################################################################
!
subroutine get_polar_coords_for_sem_fp()
  !
  !.. Local Scalars ..
  integer  :: n,sem_idx,eddy_idx
  real(wp) :: xc,yc,zc,t
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_polar_coords_for_sem_fp"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  do sem_idx = 1,num_sem_boxes
    !
    xc = sem_box(sem_idx)%center%x
    yc = sem_box(sem_idx)%center%y
    zc = sem_box(sem_idx)%center%z
    !
    ! Compute the polar coordinates for each SEM face point
    ! relative to the center of the current SEM box.
    !
    do n = 1,my_sem_pts(sem_idx)%num_fp
      !
      call yz_to_rt(my_sem_pts(sem_idx)%sem_fp(n)%y, &
                    my_sem_pts(sem_idx)%sem_fp(n)%z, &
                    my_sem_pts(sem_idx)%sem_fp(n)%r, &
                    my_sem_pts(sem_idx)%sem_fp(n)%t, &
                    yc,zc)
      !
      if (my_sem_pts(sem_idx)%sem_fp(n)%r /= zero) then
        t = my_sem_pts(sem_idx)%sem_fp(n)%t
        my_sem_pts(sem_idx)%sem_fp(n)%sin_theta = sin(t)
        my_sem_pts(sem_idx)%sem_fp(n)%cos_theta = cos(t)
      else
        my_sem_pts(sem_idx)%sem_fp(n)%sin_theta = zero
        my_sem_pts(sem_idx)%sem_fp(n)%cos_theta = zero
      end if
      !
    end do
    !
    ! If the sem_box(sem_idx)%eddy array is allocated (eddies were read in
    ! from the restart file), recompute the polar coordinates for all the
    ! eddies in this box.
    !
    if (allocated(sem_box(sem_idx)%eddy)) then
      !
      do eddy_idx = 1,sem_box(sem_idx)%num_eddies
        !
        call yz_to_rt(sem_box(sem_idx)%eddy(eddy_idx)%y, &
                      sem_box(sem_idx)%eddy(eddy_idx)%z, &
                      sem_box(sem_idx)%eddy(eddy_idx)%r, &
                      sem_box(sem_idx)%eddy(eddy_idx)%t, &
                      yc,zc)
        !
      end do
      !
    end if
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_polar_coords_for_sem_fp
!
!###############################################################################
!
subroutine get_normal_of_sem_boxes()
  !
  !.. Use Statements ..
  use metrics_mod, only : facenrm
  !
  !.. Local Scalars ..
  integer  :: n,k,num_fp,sem_idx,ierr
  real(wp) :: rnx,rny,rnz,ds
  real(wp) :: unx,uny,unz
  real(wp) :: val,aoa,boa
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: box_num(:)
  real(wp), allocatable :: box_sum(:,:)
  real(wp), allocatable :: box_nrm(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_normal_of_sem_boxes"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize all the SEM box normals to only the x-direction
  !
  sem_box(1:num_sem_boxes)%rnrm%x = one
  sem_box(1:num_sem_boxes)%rnrm%y = zero
  sem_box(1:num_sem_boxes)%rnrm%z = zero
  !
  ! Allocate the averaging arrays
  !
  box_num = init_array(num_sem_boxes,zero)   ! F2003 AUTO-REALLOC
  box_sum = init_array(3,num_sem_boxes,zero) ! F2003 AUTO-REALLOC
  box_nrm = init_array(3,num_sem_boxes,zero) ! F2003 AUTO-REALLOC
  !
  ! Get the average contribution to each SEM box from this MPI process
  !
  do sem_idx = 1,num_sem_boxes
    !
    num_fp = my_sem_pts(sem_idx)%num_fp
    !
    box_num(sem_idx) = num_fp
    !
    do n = 1,num_fp
      !
      k = my_sem_pts(sem_idx)%sem_fp(n)%idx
      !
      ! We need the inward facing unit normal for this face point
      ! NOTE: facenrm stores the outward facing area-weighted normal. Taking
      !       the negative value of facenrm gives us the inward facing
      !       area-weighted normal. We dont want to average the area-weighted
      !       normals, so we need to convert these to unit normals.
      !
      rnx = -facenrm(1,k)
      rny = -facenrm(2,k)
      rnz = -facenrm(3,k)
      ds = sqrt(rnx*rnx + rny*rny + rnz*rnz)
      unx = rnx/ds
      uny = rny/ds
      unz = rnz/ds
      !
      ! Add the unit normals to the average summation
      !
      box_sum(1,sem_idx) = box_sum(1,sem_idx) + unx
      box_sum(2,sem_idx) = box_sum(2,sem_idx) + uny
      box_sum(3,sem_idx) = box_sum(3,sem_idx) + unz
      !
    end do
    !
  end do
  !
  ! Sum the average contributions across all MPI processes
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,box_num,size(box_num,kind=int_mpi), &
                       mpi_inttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
    call mpi_allreduce(MPI_IN_PLACE,box_sum,size(box_sum,kind=int_mpi), &
                       mpi_flttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Finish computing the average and save to the SEM box
  !
  do sem_idx = 1,num_sem_boxes
    !
    val = one / real(box_num(sem_idx),kind=wp)
   !!
   !rnx = box_sum(1,sem_idx) * val
   !rny = box_sum(2,sem_idx) * val
   !rnz = box_sum(3,sem_idx) * val
   !ds = sqrt(rnx*rnx + rny*rny + rnz*rnz)
   !unx = rnx/ds
   !uny = rny/ds
   !unz = rnz/ds
   !!
   !box_nrm(1,sem_idx) = unx
   !box_nrm(2,sem_idx) = unxuny
   !box_nrm(3,sem_idx) = unxunz
    !
    box_nrm(1,sem_idx) = box_sum(1,sem_idx) * val
    box_nrm(2,sem_idx) = box_sum(2,sem_idx) * val
    box_nrm(3,sem_idx) = box_sum(3,sem_idx) * val
    !
  end do
  !
  ! Reset box_sum back to 0
  !
  box_sum(:,:) = zero
  !
  ! Check the averaged norms against all the flux point normals
  !
  do sem_idx = 1,num_sem_boxes
    !
    do n = 1,my_sem_pts(sem_idx)%num_fp
      !
      k = my_sem_pts(sem_idx)%sem_fp(n)%idx
      !
      ! We need the inward facing unit normal for this face point
      !
      rnx = -facenrm(1,k)
      rny = -facenrm(2,k)
      rnz = -facenrm(3,k)
      ds = sqrt(rnx*rnx + rny*rny + rnz*rnz)
      unx = rnx/ds
      uny = rny/ds
      unz = rnz/ds
      !
      ! Check the averaged unit normal agains this unit normal
      !
      if (abs(box_nrm(1,sem_idx)-unx) > eps3) then
        box_sum(1,sem_idx) = box_sum(1,sem_idx) + one
      end if
      !
      if (abs(box_nrm(2,sem_idx)-uny) > eps3) then
        box_sum(2,sem_idx) = box_sum(2,sem_idx) + one
      end if
      !
      if (abs(box_nrm(3,sem_idx)-unz) > eps3) then
        box_sum(3,sem_idx) = box_sum(3,sem_idx) + one
      end if
      !
    end do
    !
  end do
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,box_sum,size(box_sum,kind=int_mpi), &
                       mpi_flttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
  end if
  !
  do sem_idx = 1,num_sem_boxes
    !
    if (any( nint(box_sum(:,sem_idx)) > 0)) then
      !
     !if (mypnum == glb_root) then
     !  !
     !  write (iout,1) sem_idx
     !  !
     !  if (nint(box_sum(1,sem_idx)) > 0) then
     !    write (iout,2) "X",nint(box_sum(1,sem_idx))
     !  end if
     !  if (nint(box_sum(2,sem_idx)) > 0) then
     !    write (iout,2) "Y",nint(box_sum(2,sem_idx))
     !  end if
     !  if (nint(box_sum(3,sem_idx)) > 0) then
     !    write (iout,2) "Z",nint(box_sum(3,sem_idx))
     !  end if
     !  !
     !end if
      !
      ! Set the SEM box unit normal in only the x direction
      !
      sem_box(sem_idx)%rnrm%x = one
      sem_box(sem_idx)%rnrm%y = zero
      sem_box(sem_idx)%rnrm%z = zero
      !
    else
      !
      unx = box_nrm(1,sem_idx)
      uny = box_nrm(2,sem_idx)
      unz = box_nrm(3,sem_idx)
      !
      sem_box(sem_idx)%rnrm%x = unx
      sem_box(sem_idx)%rnrm%y = uny
      sem_box(sem_idx)%rnrm%z = unz
      !
      aoa = atan2(uny,unx)
      boa = atan2(unz,sqrt(unx*unx+uny*uny))
      !
      sem_box(sem_idx)%alpha_aoa = aoa
      sem_box(sem_idx)%beta_aoa  = boa
      !
    end if
    !
  end do
  !
  !
  !
  ! Deallocate box_num, box_sum, and box_nrm before we leave
  !
  if (allocated(box_num)) then
    deallocate ( box_num , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"box_num",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(box_sum)) then
    deallocate ( box_sum , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"box_sum",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(box_nrm)) then
    deallocate ( box_nrm , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"box_nrm",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_normal_of_sem_boxes
!
!###############################################################################
!
subroutine get_size_of_sem_boxes()
  !
  !.. Use Statements ..
  use ovar, only : SEM_eddy_volume_is_spherical
  use ovar, only : SEM_box_num_eddies_multiplier
  !
  !.. Local Scalars ..
  integer  :: ierr,n,num_eddies,sem_idx
  real(wp) :: lturb,t_lturb,box_volume
  real(wp) :: xmin,xmax,x
  real(wp) :: ymin,ymax,y
  real(wp) :: zmin,zmax,z
  real(wp) :: rmin,rmax,r
  real(wp) :: tmin,tmax,t
  real(wp) :: mean_eddy_volume
  real(wp) :: rnum_eddies
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: box_min(:,:)
  real(wp), allocatable :: box_max(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_size_of_sem_boxes"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate box_min and box_max to store the min/max extents
  ! (local to this partition/process) of the SEM boxes
  !
  box_min = init_array(5,num_sem_boxes, huge(zero)) ! F2003 AUTO-REALLOC
  box_max = init_array(5,num_sem_boxes,-huge(zero)) ! F2003 AUTO-REALLOC
  !
  ! Compute the min/max extents (local to this
  ! partition/process) of the SEM boxes
  !
  do sem_idx = 1,num_sem_boxes
    !
    lturb = sem_box(sem_idx)%length_scale
    !
    ! Compute the x,y,z min/max for the current SEM box
    !
    do n = 1,my_sem_pts(sem_idx)%num_fp
      !
      x = my_sem_pts(sem_idx)%sem_fp(n)%x
      y = my_sem_pts(sem_idx)%sem_fp(n)%y
      z = my_sem_pts(sem_idx)%sem_fp(n)%z
      r = my_sem_pts(sem_idx)%sem_fp(n)%r
      t = my_sem_pts(sem_idx)%sem_fp(n)%t
      !
      xmin = x - lturb
      xmax = x + lturb
      ymin = y - lturb
      ymax = y + lturb
      zmin = z - lturb
      zmax = z + lturb
      rmin = r - lturb
      rmax = r + lturb
      !
      ! The turbulence length scale is an arclength for theta, so determine
      ! the angle that corresponds to an arclength of lturb,
      ! i.e. theta_length_scale = lturb/r
      !
      if (r == zero) then
        t_lturb = zero
      else
        t_lturb = lturb/r
      end if
      !
      tmin = t - t_lturb
      tmax = t + t_lturb
      !
      ! Make sure rmin, tmin, and tmax are realizable
      !
      rmin = max(zero,rmin)
      tmin = max(zero,tmin)
      tmax = min(two*pi,tmax)
     !tmin = max(zero,min(two*pi,tmin))
     !tmax = min(two*pi,max(zero,tmax))
      !
      box_min(1,sem_idx) = min(xmin,box_min(1,sem_idx))
      box_min(2,sem_idx) = min(ymin,box_min(2,sem_idx))
      box_min(3,sem_idx) = min(zmin,box_min(3,sem_idx))
      box_min(4,sem_idx) = min(rmin,box_min(4,sem_idx))
      box_min(5,sem_idx) = min(tmin,box_min(5,sem_idx))
      !
      box_max(1,sem_idx) = max(xmax,box_max(1,sem_idx))
      box_max(2,sem_idx) = max(ymax,box_max(2,sem_idx))
      box_max(3,sem_idx) = max(zmax,box_max(3,sem_idx))
      box_max(4,sem_idx) = max(rmax,box_max(4,sem_idx))
      box_max(5,sem_idx) = max(tmax,box_max(5,sem_idx))
      !
    end do
   !!
   !! Compute the r,theta min/max for the current box from the y,z min/max
   !!
   !if (my_sem_pts(sem_idx)%num_fp > 0) then
   !  !
   !  ymin = box_min(2,sem_idx)
   !  zmin = box_min(3,sem_idx)
   !  !
   !  ymax = box_max(2,sem_idx)
   !  zmax = box_max(3,sem_idx)
   !  !
   !  r1 = hypot( ymin , zmin )
   !  r2 = hypot( ymax , zmin )
   !  r3 = hypot( ymin , zmax )
   !  r4 = hypot( ymax , zmax )
   !  !
   ! !r1 = sqrt( ymin*ymin + zmin*zmin )
   ! !r2 = sqrt( ymax*ymax + zmin*zmin )
   ! !r3 = sqrt( ymin*ymin + zmax*zmax )
   ! !r4 = sqrt( ymax*ymax + zmax*zmax )
   !  !
   !  rmin = max( zero , min(r1,r2,r3,r3) )
   !  rmax = max( zero , max(r1,r2,r3,r4) )
   !  !
   !  t1 = zero
   !  t2 = zero
   !  t3 = zero
   !  t4 = zero
   !  !
   !  if (ymin /= zero .or. zmin /= zero) t1 = atan2( zmin , ymin )
   !  if (ymax /= zero .or. zmin /= zero) t2 = atan2( zmin , ymax )
   !  if (ymin /= zero .or. zmax /= zero) t3 = atan2( zmax , ymin )
   !  if (ymax /= zero .or. zmax /= zero) t4 = atan2( zmax , ymax )
   !  !
   !  ! The atan2 function returns a value on the interval [-Pi,Pi].
   !  ! Converting an angle from this interval to the interval [0,2*Pi]
   !  ! simply requires adding 2*Pi to any angle value less than zero.
   !  ! For example:
   !  !    theta[-Pi,Pi] = -Pi    => theta[0,2*Pi] = Pi     = 180 degrees
   !  !    theta[-Pi,Pi] = -Pi/2  => theta[0,2*Pi] = 3*Pi/2 = 270 degrees
   !  !    theta[-Pi,Pi] = -0     => theta[0,2*Pi] = 2*Pi   = 360 degrees
   !  !
   !  if (t1 < zero) t1 = t1 + two*pi
   !  if (t2 < zero) t2 = t2 + two*pi
   !  if (t3 < zero) t3 = t3 + two*pi
   !  if (t4 < zero) t4 = t4 + two*pi
   !  !
   !  tmin = min(t1,t2,t3,t4)
   !  tmax = max(t1,t2,t3,t4)
   !  !
   !  box_min(4,sem_idx) = rmin
   !  box_min(5,sem_idx) = tmin
   !  !
   !  box_max(4,sem_idx) = rmax
   !  box_max(5,sem_idx) = tmax
   !  !
   !end if
    !
  end do
  !
!#define OUTPUT_SEM_BOX_SIZES
#ifdef OUTPUT_SEM_BOX_SIZES
  do n = 0,ncpu-1
    if (mypnum == n) then
      write (iout,*)
      do sem_idx = 1,num_sem_boxes
        write (iout,1) mypnum,sem_idx,"x",box_min(1,sem_idx),box_max(1,sem_idx)
        write (iout,1) mypnum,sem_idx,"y",box_min(2,sem_idx),box_max(2,sem_idx)
        write (iout,1) mypnum,sem_idx,"z",box_min(3,sem_idx),box_max(3,sem_idx)
        write (iout,1) mypnum,sem_idx,"r",box_min(4,sem_idx),box_max(4,sem_idx)
        write (iout,1) mypnum,sem_idx,"t",box_min(5,sem_idx),box_max(5,sem_idx)
      end do
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    flush (iout)
  end do
  if (mypnum == glb_root) write (iout,*)
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  flush (iout)
  1 format (2x,"CPU #",i2," : SEM Box # ",i2, &
            2x,a1,"[min,max] = [",es14.6,",",es14.6,"]")
#endif
  !
  ! Reduce the box_min and box_max arrays across all processes so that
  ! every process has the true global min/max extents for each SEM box
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,box_min,size(box_min,kind=int_mpi), &
                       mpi_flttyp,MPI_MIN,MPI_COMM_WORLD,mpierr)
    call mpi_allreduce(MPI_IN_PLACE,box_max,size(box_max,kind=int_mpi), &
                       mpi_flttyp,MPI_MAX,MPI_COMM_WORLD,mpierr)
  end if
  !
#ifdef OUTPUT_SEM_BOX_SIZES
  do n = 0,ncpu-1
    if (mypnum == n) then
      write (iout,*)
      do sem_idx = 1,num_sem_boxes
        write (iout,1) mypnum,sem_idx,"x",box_min(1,sem_idx),box_max(1,sem_idx)
        write (iout,1) mypnum,sem_idx,"y",box_min(2,sem_idx),box_max(2,sem_idx)
        write (iout,1) mypnum,sem_idx,"z",box_min(3,sem_idx),box_max(3,sem_idx)
        write (iout,1) mypnum,sem_idx,"r",box_min(4,sem_idx),box_max(4,sem_idx)
        write (iout,1) mypnum,sem_idx,"t",box_min(5,sem_idx),box_max(5,sem_idx)
      end do
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    flush (iout)
  end do
  if (mypnum == glb_root) write (iout,*)
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  flush (iout)
#endif
  !
  ! Copy the global min/max extents for each SEM box into sem_box.
  ! Also compute the volume and the number of eddies for each box.
  !
  do sem_idx = 1,num_sem_boxes
    !
    lturb = sem_box(sem_idx)%length_scale
    !
    xmin = box_min(1,sem_idx)
    ymin = box_min(2,sem_idx)
    zmin = box_min(3,sem_idx)
    rmin = box_min(4,sem_idx)
    tmin = box_min(5,sem_idx)
    !
    xmax = box_max(1,sem_idx)
    ymax = box_max(2,sem_idx)
    zmax = box_max(3,sem_idx)
    rmax = box_max(4,sem_idx)
    tmax = box_max(5,sem_idx)
    !
    ! Compute the volume of this SEM box
    !
    if (sem_box(sem_idx)%is_cylindrical) then
      box_volume = half*(xmax-xmin)*(tmax-tmin)*(rmax**2-rmin**2)
    else
      box_volume = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
    end if
    !
    ! If the SEM eddies were not read in from a restart file,
    ! i.e., sem_box(sem_idx)%eddy is not yet allocated,
    ! compute the number of eddies to use for this box based on the box
    ! volume and the turbulence length scale specified for this SEM box.
    !
    if (.not. allocated(sem_box(sem_idx)%eddy)) then
      !
      ! NOTE: The input parameter SEM_eddy_volume_is_spherical determines
      !       whether the variable mean_eddy_volume is equal to either
      !       (4/3)*PI*(1/2)**2 or 1 depending on whether the eddy is
      !       treated as sphere-like or box-like in volume, respectively.
      !       If a box-like eddy has a volume of l**3, the equivalent
      !       volume of a sphere-like eddy is (4/3)*pi*((l/2)**3).
      !
      if (SEM_eddy_volume_is_spherical) then
        mean_eddy_volume = (four/three)*pi*((half*lturb)**3)
      else
        mean_eddy_volume = lturb**3
      end if
      !
      rnum_eddies = SEM_box_num_eddies_multiplier*box_volume / mean_eddy_volume
      !
      num_eddies = min( max_eddies , max( min_eddies , int( rnum_eddies ) ) )
      !
      ! Save the number of eddies for this box into sem_box
      !
      sem_box(sem_idx)%num_eddies = num_eddies
      !
    end if
    !
    ! Save the box min/max extents and the box volume into sem_box.
    !
    sem_box(sem_idx)%boxmin(1) = xmin
    sem_box(sem_idx)%boxmin(2) = ymin
    sem_box(sem_idx)%boxmin(3) = zmin
    sem_box(sem_idx)%boxmin(4) = rmin
    sem_box(sem_idx)%boxmin(5) = tmin
    !
    sem_box(sem_idx)%boxmax(1) = xmax
    sem_box(sem_idx)%boxmax(2) = ymax
    sem_box(sem_idx)%boxmax(3) = zmax
    sem_box(sem_idx)%boxmax(4) = rmax
    sem_box(sem_idx)%boxmax(5) = tmax
    !
    sem_box(sem_idx)%box_volume = box_volume
    !
  end do
  call dump_sem_data("preread")
  !
  ! Deallocate box_min and box_max before we leave
  !
  if (allocated(box_min)) then
    deallocate ( box_min , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"box_min",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(box_max)) then
    deallocate ( box_max , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"box_max",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_size_of_sem_boxes
!
!###############################################################################
!
subroutine recompute_sem_box_volume()
  !
  !.. Use Statements ..
  use ovar, only : SEM_eddy_volume_is_spherical
  use ovar, only : SEM_box_num_eddies_multiplier
  !
  !.. Local Scalars ..
  integer  :: sem_idx,num_eddies
  real(wp) :: xmin,ymin,zmin,rmin,tmin
  real(wp) :: xmax,ymax,zmax,rmax,tmax
  real(wp) :: lturb,box_volume
  real(wp) :: mean_eddy_volume
  real(wp) :: rnum_eddies
  !
  ! Copy the global min/max extents for each SEM box into sem_box.
  ! Also compute the volume and the number of eddies for each box.
  !
  do sem_idx = 1,num_sem_boxes
    !
    lturb = sem_box(sem_idx)%length_scale
    !
    xmin = sem_box(sem_idx)%boxmin(1)
    ymin = sem_box(sem_idx)%boxmin(2)
    zmin = sem_box(sem_idx)%boxmin(3)
    rmin = sem_box(sem_idx)%boxmin(4)
    tmin = sem_box(sem_idx)%boxmin(5)
    !
    xmax = sem_box(sem_idx)%boxmax(1)
    ymax = sem_box(sem_idx)%boxmax(2)
    zmax = sem_box(sem_idx)%boxmax(3)
    rmax = sem_box(sem_idx)%boxmax(4)
    tmax = sem_box(sem_idx)%boxmax(5)
    !
    ! Compute the volume of this SEM box
    !
    if (sem_box(sem_idx)%is_cylindrical) then
      box_volume = half*(xmax-xmin)*(tmax-tmin)*(rmax**2-rmin**2)
    else
      box_volume = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
    end if
    !
    ! If the SEM eddies were not read in from a restart file,
    ! i.e., sem_box(sem_idx)%eddy is not yet allocated,
    ! compute the number of eddies to use for this box based on the box
    ! volume and the turbulence length scale specified for this SEM box.
    !
    if (.not. allocated(sem_box(sem_idx)%eddy)) then
      !
      ! NOTE: The input parameter SEM_eddy_volume_is_spherical determines
      !       whether the variable mean_eddy_volume is equal to either
      !       (4/3)*PI*(1/2)**2 or 1 depending on whether the eddy is
      !       treated as sphere-like or box-like in volume, respectively.
      !       If a box-like eddy has a volume of l**3, the equivalent
      !       volume of a sphere-like eddy is (4/3)*pi*((l/2)**3).
      !
      if (SEM_eddy_volume_is_spherical) then
        mean_eddy_volume = (four/three)*pi*((half*lturb)**3)
      else
        mean_eddy_volume = lturb**3
      end if
      !
      rnum_eddies = SEM_box_num_eddies_multiplier*box_volume / mean_eddy_volume
      !
      num_eddies = min( max_eddies , max( min_eddies , int( rnum_eddies ) ) )
      !
      ! Save the number of eddies for this box into sem_box
      !
      sem_box(sem_idx)%num_eddies = num_eddies
      !
    end if
    !
    sem_box(sem_idx)%box_volume = box_volume
    !
  end do
  !
end subroutine recompute_sem_box_volume
!
!###############################################################################
!
subroutine check_restart_eddies()
  !
  !.. Use Statements ..
  use ovar, only : SEM_eddy_volume_is_spherical
  use ovar, only : dtrst
  !
  !.. Local Scalars ..
  integer  :: ierr,n,sem_idx,eddy_idx
  real(wp) :: vol_coef,eddy_vol,elx,ely,elz
  real(wp) :: x_offset
  !
  logical(lk) :: error_found
  !
  !.. Local Arrays ..
  real(wp) :: eddy_loc(1:3)
  !
  character(len=1) :: cxyz(1:3)
  !
  !.. Local Allocatable Arrays ..
  integer,         allocatable :: num_bad_eddies(:)
  real(wp),        allocatable :: box_min(:,:)
  real(wp),        allocatable :: box_max(:,:)
  real(wp),        allocatable :: eddy_min(:,:)
  real(wp),        allocatable :: eddy_max(:,:)
  logical(lk),     allocatable :: is_inside(:,:)
  type(sem_box_t), allocatable :: tmp_sem_box(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_restart_eddies"
  character(len=*), parameter :: char_xyz(1:3) = ["x","y","z"]
  character(len=*), parameter :: char_xrt(1:3) = ["x","r","t"]
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  num_bad_eddies = init_array(num_sem_boxes,0) ! F2003 AUTO-REALLOC
  !
  box_min = init_array(3,num_sem_boxes, huge(zero)) ! F2003 AUTO-REALLOC
  box_max = init_array(3,num_sem_boxes,-huge(zero)) ! F2003 AUTO-REALLOC
  !
  ! Compute the min/max extents (local to this
  ! partition/process) of the SEM boxes
  !
  loop_over_sem_boxes: &
  do sem_idx = 1,num_sem_boxes
    !
    x_offset = half*dtrst*sem_box(sem_idx)%convective_velocity
    !
    ! Get the extents of the current SEM box
    !
   !box_min(1,sem_idx) = sem_box(sem_idx)%boxmin(1) - eps6
   !box_max(1,sem_idx) = sem_box(sem_idx)%boxmax(1) + eps6
    box_min(1,sem_idx) = sem_box(sem_idx)%boxmin(1) - x_offset
    box_max(1,sem_idx) = sem_box(sem_idx)%boxmax(1) + x_offset
    !
    if (sem_box(sem_idx)%is_cylindrical) then
      !
      box_min(2,sem_idx) = sem_box(sem_idx)%boxmin(4)
      box_min(3,sem_idx) = sem_box(sem_idx)%boxmin(5)
      !
      box_max(2,sem_idx) = sem_box(sem_idx)%boxmax(4)
      box_max(3,sem_idx) = sem_box(sem_idx)%boxmax(5)
      !
    else
      !
      box_min(2,sem_idx) = sem_box(sem_idx)%boxmin(2)
      box_min(3,sem_idx) = sem_box(sem_idx)%boxmin(3)
      !
      box_max(2,sem_idx) = sem_box(sem_idx)%boxmax(2)
      box_max(3,sem_idx) = sem_box(sem_idx)%boxmax(3)
      !
    end if
    !
    if (.not. allocated(sem_box(sem_idx)%eddy)) cycle loop_over_sem_boxes
    !
    n = 0
    !
    loop_over_eddies: &
    do eddy_idx = 1,sem_box(sem_idx)%num_eddies
      !
      eddy_loc(1) = sem_box(sem_idx)%eddy(eddy_idx)%x
      !
      if (sem_box(sem_idx)%is_cylindrical) then
        eddy_loc(2) = sem_box(sem_idx)%eddy(eddy_idx)%r
        eddy_loc(3) = sem_box(sem_idx)%eddy(eddy_idx)%t
      else
        eddy_loc(2) = sem_box(sem_idx)%eddy(eddy_idx)%y
        eddy_loc(3) = sem_box(sem_idx)%eddy(eddy_idx)%z
      end if
      !
      if (any(eddy_loc(1:3) < box_min(1:3,sem_idx)) .or. &
          any(eddy_loc(1:3) > box_max(1:3,sem_idx))) then
        n = n + 1
      end if
      !
    end do loop_over_eddies
    !
    num_bad_eddies(sem_idx) = n
    !
  end do loop_over_sem_boxes
  !
  ! Check if there were any eddies that found outside their respective SEM box.
  ! If any were found, write out a warning. Then, try to troubleshoot the boxes
  ! to see if the restart SEM boxes and input SEM boxes might simply be out of
  ! order. Finally, rearrange the restart SEM boxes if doing so matches them
  ! up with the input SEM boxes.
  !
  if (sum(num_bad_eddies) > 0) then
    !
    ! Output a warning that some of the eddies read from the restart file
    ! have been found to be outside of their respective SEM box.
    !
    if (mypnum == glb_root) then
      write (iout,1)
      do sem_idx = 1,num_sem_boxes
        if (num_bad_eddies(sem_idx) > 0) then
          write (iout,2) num_bad_eddies(sem_idx),sem_idx
        end if
      end do
      write (iout,3)
    end if
    flush (iout)
    !
    ! Deallocate num_bad_eddies since it is no longer needed
    !
    if (allocated(num_bad_eddies)) then
      deallocate ( num_bad_eddies , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"num_bad_eddies",2,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
    !
    ! Find the min/max extents of the eddy locations for each SEM box
    !
    eddy_min = init_array(3,num_sem_boxes, huge(zero)) ! F2003 AUTO-REALLOC
    eddy_max = init_array(3,num_sem_boxes,-huge(zero)) ! F2003 AUTO-REALLOC
    !
    do sem_idx = 1,num_sem_boxes
      !
      if (.not. allocated(sem_box(sem_idx)%eddy)) then
        !
        ! If the eddy array is not allocated, just make the min/max extents
        ! of the eddies equal to the box dimensions +/- 10**-6
        !
        eddy_min(1,sem_idx) = sem_box(sem_idx)%boxmin(1) + eps6
        eddy_max(1,sem_idx) = sem_box(sem_idx)%boxmax(1) - eps6
        !
        if (sem_box(sem_idx)%is_cylindrical) then
          !
          eddy_min(2,sem_idx) = sem_box(sem_idx)%boxmin(4) + eps6
          eddy_max(2,sem_idx) = sem_box(sem_idx)%boxmax(4) - eps6
          !
          eddy_min(3,sem_idx) = sem_box(sem_idx)%boxmin(5) + eps6
          eddy_max(3,sem_idx) = sem_box(sem_idx)%boxmax(5) - eps6
          !
        else
          !
          eddy_min(2,sem_idx) = sem_box(sem_idx)%boxmin(2) + eps6
          eddy_max(2,sem_idx) = sem_box(sem_idx)%boxmax(2) - eps6
          !
          eddy_min(3,sem_idx) = sem_box(sem_idx)%boxmin(3) + eps6
          eddy_max(3,sem_idx) = sem_box(sem_idx)%boxmax(3) - eps6
          !
        end if
        !
      else
        !
        ! Get the minimum value of each of the coordinates
        ! over all the eddies in this SEM box
        !
        eddy_min(1,sem_idx) = minval(sem_box(sem_idx)%eddy(:)%x)
        !
        if (sem_box(sem_idx)%is_cylindrical) then
          eddy_min(2,sem_idx) = minval(sem_box(sem_idx)%eddy(:)%r)
          eddy_min(3,sem_idx) = minval(sem_box(sem_idx)%eddy(:)%t)
        else
          eddy_min(2,sem_idx) = minval(sem_box(sem_idx)%eddy(:)%y)
          eddy_min(3,sem_idx) = minval(sem_box(sem_idx)%eddy(:)%z)
        end if
        !
        ! Get the maximum value of each of the coordinates
        ! over all the eddies in this SEM box
        !
        eddy_max(1,sem_idx) = maxval(sem_box(sem_idx)%eddy(:)%x)
        !
        if (sem_box(sem_idx)%is_cylindrical) then
          eddy_max(2,sem_idx) = maxval(sem_box(sem_idx)%eddy(:)%r)
          eddy_max(3,sem_idx) = maxval(sem_box(sem_idx)%eddy(:)%t)
        else
          eddy_max(2,sem_idx) = maxval(sem_box(sem_idx)%eddy(:)%y)
          eddy_max(3,sem_idx) = maxval(sem_box(sem_idx)%eddy(:)%z)
        end if
        !
      end if
      !
    end do
    !
    ! Check the eddy extents against the SEM box extents to see if we can
    ! rearrange the eddies to get a unique 1-1 matching of eddies to SEM boxes.
    !
    is_inside = init_array(num_sem_boxes,num_sem_boxes,fals) ! F03 AUTO-REALLOC
    !
    do sem_idx = 1,num_sem_boxes
      do eddy_idx = 1,num_sem_boxes
        !
        if (all(eddy_min(1:3,eddy_idx) >= box_min(1:3,sem_idx)) .and. &
            all(eddy_max(1:3,eddy_idx) <= box_max(1:3,sem_idx))) then
          is_inside(eddy_idx,sem_idx) = true
        end if
        !
      end do
    end do
    !
    ! Check to make sure that each group of eddies was found to be in only
    ! one SEM box, and that each SEM box only contains one group of eddies.
    !
    error_found = fals
    !
    ! Check to make sure at most one group
    ! of eddies is contained by each SEM box
    !
    do sem_idx = 1,num_sem_boxes
      n = count( is_inside(1:num_sem_boxes,sem_idx) )
      ! NOTE: n IS allowed to be 0 here, which means that none of
      !       the restart eddy groups is contained by this SEM box;
      !       therefore, the eddies for this SEM box will need to
      !       be initialized to new eddies.
      if (n > 1) then
        if (mypnum == glb_root) write (iout,10) sem_idx,n
        error_found = true
      end if
    end do
    !
    ! Check to make sure that each group of
    ! eddies is contained by only one SEM box.
    !
    do eddy_idx = 1,num_sem_boxes
      n = count( is_inside(eddy_idx,1:num_sem_boxes) )
      ! NOTE: n is NOT allowed to be 0 here, which would mean that the
      !       current group of eddies from the restart file is not
      !       contained by any of the current SEM boxes. Each group of
      !       eddies from the restart file MUST be contained by one,
      !       and only one, SEM box.
      if (n /= 1) then
        if (mypnum == glb_root) write (iout,11) eddy_idx,n
        error_found = true
      end if
    end do
    !
    if (error_found) then
      if (mypnum == glb_root) then
        write (iout,12) (n, n=1,num_sem_boxes)
        do sem_idx = 1,num_sem_boxes
          write (iout,13) merge("SEM Box # ","          ",sem_idx==1), &
                          sem_idx,(is_inside(n,sem_idx), n=1,num_sem_boxes)
        end do
        write (iout,'(//)')
        !
        do sem_idx = 1,num_sem_boxes
          if (sem_box(sem_idx)%is_cylindrical) then
            cxyz(1:3) = char_xrt(1:3)
          else
            cxyz(1:3) = char_xyz(1:3)
          end if
          write (iout,20) sem_idx,(cxyz(n), &
                                   box_min(n,sem_idx), &
                                   box_max(n,sem_idx), n=1,3)
          write (iout,21) sem_idx,(cxyz(n), &
                                   eddy_min(n,sem_idx), &
                                   eddy_max(n,sem_idx), n=1,3)
        end do
      end if
      flush (iout)
      write (error_message,14)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
    !
    ! If we passed the previous tests, there is a 1-to-1 matching of SEM boxes
    ! to groups of eddies. Shuffle the eddy components of the sem_box array
    ! so that the groups of eddies are matched with their respective SEM box.
    !
    ! First create a duplicate of sem_box
    ! NOTE: This will deep-copy everything within sem_box and its underlying
    !       components over to tmp_sem_box.
    !
    tmp_sem_box = sem_box ! F2003 AUTO-REALLOC
    !
    shuffle_eddies: &
    do sem_idx = 1,num_sem_boxes
      !
      ! Get the index for the group of eddies that matches this SEM box
      ! NOTE: The sum function is used here on an integer sequence from
      !       1 to num_sem_boxes, because the mask array should contain
      !       only one true value thus simply making the sum equal to
      !       the index of the first dimension of is_inside whose value
      !       is true.
      !
      eddy_idx = sum( intseq(1,num_sem_boxes) , &
                      mask=is_inside(1:num_sem_boxes,sem_idx) )
      !
      ! Use move_alloc to copy tmp_sem_box(eddy_idx)%eddy over to
      ! sem_box(sem_idx)%eddy. This will do the following:
      !   1. deallocate sem_box(sem_idx)%eddy if already allocated
      !   2. allocate sem_box(sem_idx)%eddy to the same size as
      !      tmp_sem_box(eddy_idx)%eddy
      !   3. copy the contents of tmp_sem_box(eddy_idx)%eddy over to
      !      sem_box(sem_idx)%eddy
      !   4. deallocate tmp_sem_box(eddy_idx)
      !
      if (eddy_idx > 0) then
        call move_alloc( from=tmp_sem_box(eddy_idx)%eddy , &
                         to=sem_box(sem_idx)%eddy )
      end if
      !
    end do shuffle_eddies
    !
    ! Deallocate box_min, box_max, eddy_min, and eddy_max
    ! since they are no longer needed
    !
    if (allocated(box_min)) then
      deallocate ( box_min , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"box_min",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    if (allocated(box_max)) then
      deallocate ( box_max , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"box_max",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    if (allocated(eddy_min)) then
      deallocate ( eddy_min , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"eddy_min",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    if (allocated(eddy_max)) then
      deallocate ( eddy_max , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"eddy_max",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
  end if
  !
  ! Deallocate tmp_sem_box and is_inside because they are no longer needed.
  !
  if (allocated(tmp_sem_box)) then
    deallocate ( tmp_sem_box , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"tmp_sem_box",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(is_inside)) then
    deallocate ( is_inside , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"is_inside",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Finally, compute the volume of each eddy using the eddy frequencies
  ! (i.e., one/length-scales) since this information is not read in
  ! from the restart file.
  !
  vol_coef = merge( pi/six , one , SEM_eddy_volume_is_spherical )
  !
  do sem_idx = 1,num_sem_boxes
    !
    if (.not. allocated(sem_box(sem_idx)%eddy)) cycle
    !
    do eddy_idx = 1,sem_box(sem_idx)%num_eddies
      !
      elx = one / sem_box(sem_idx)%eddy(eddy_idx)%x_frequency
      ely = one / sem_box(sem_idx)%eddy(eddy_idx)%y_frequency
      elz = one / sem_box(sem_idx)%eddy(eddy_idx)%z_frequency
      !
      eddy_vol = vol_coef*elx*ely*elz
      !
      sem_box(sem_idx)%eddy(eddy_idx)%volume = one / eddy_vol
      !
    end do
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1  format (/, &
       2x,78("#"),/, &
       2x,5("#"),26x, &
         "!!! WARNING !!!",27x,5("#"),/, &
       2x,5("#"),10x, &
         "A MISMATCH WAS DETECTED BETWEEN THE LOCATION OF",11x,5("#"),/, &
       2x,5("#"),10x, &
         "THE RESTART EDDIES AND THE SEM BOUNDING BOXES!!!",10x,5("#"))
  2  format (2x,5("#"),10x, &
             i8," restart eddies are outside SEM box # ",i2,10x,5("#"))
  3  format (2x,78("#"),/)
  !
  10 format (/,2x,"SEM box # ",i0," appears to contain ",i0," groups of ", &
                  "restart eddies!")
  11 format (/,2x,"Restart eddy group # ",i0," appears to be contained ", &
                  "by ",i0," SEM boxes!")
  12 format (//,20x,"Restart Eddy Group #",/, &
                2x,"Contained by",6x,i2,100(2x,i2,:))
  13 format (5x,a10,i2," [",100(2x,l1," ]",:,tl2))
  14 format ("ERROR!!! An error occurred while trying to create a 1-to-1 ", &
             "matching between the groups of eddies read from the restart ", &
             "file and the current set of SEM bounding boxes!")
  !
  20 format (/,2x,"Coordinate range for SEM box # ",i0," : ",/, &
               6x,a1," = [",es14.6,",",es14.6,"]",/, &
               6x,a1," = [",es14.6,",",es14.6,"]",/, &
               6x,a1," = [",es14.6,",",es14.6,"]")
  21 format (/,2x,"Coordinate range for Restart Eddy Group # ",i0," : ",/, &
               6x,a1," = [",es14.6,",",es14.6,"]",/, &
               6x,a1," = [",es14.6,",",es14.6,"]",/, &
               6x,a1," = [",es14.6,",",es14.6,"]")
  !
end subroutine check_restart_eddies
!
!###############################################################################
!
subroutine read_sem_box_profiles()
  !
  !.. Use Statements
  use ovar, only : bc_in
  !
  !.. Local Scalars ..
  integer  :: n,npts,sem_idx
  real(wp) :: tau11,tau21,tau22,tau31,tau32,tau33
  real(wp) :: a11,a21,a22,a31,a32,a33
  real(wp) :: ravg,uavg,vavg,wavg,tavg
  !
  logical(lk) :: profile_data_is_allocated
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_sem_box_profiles"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! First, check to see if the profile data has already been allocated.
  ! If it has been allocated, it is assumed that the profile data also
  ! exists in these allocated arrays.
  !
  call check_sem_profile_allocations(profile_data_is_allocated)
  !
  ! If the profile data has already been allocated, we do not need to read
  ! the data in so exit from this subroutine.
  !
  if (profile_data_is_allocated) then
    call debug_timer(leaving_procedure,pname)
    return
  end if
  !
  ! We need to create the profile data, so read the data from the profile files.
  !
  do sem_idx = 1,num_sem_boxes
    !
    call read_sem_profile(sem_box(sem_idx)%profile, &
                          is_cylindrical=sem_box(sem_idx)%is_cylindrical)
    !
    npts = sem_box(sem_idx)%profile%npts
    !
    ! Allocate any Reynolds stress tensor arrays
    ! that were not read from the profile file.
    !
    if (.not. allocated(sem_box(sem_idx)%profile%tau11)) &
    sem_box(sem_idx)%profile%tau11 = init_array(npts,zero) ! F2003 AUTO-REALLOC
    !
    if (.not. allocated(sem_box(sem_idx)%profile%tau21)) &
    sem_box(sem_idx)%profile%tau21 = init_array(npts,zero) ! F2003 AUTO-REALLOC
    !
    if (.not. allocated(sem_box(sem_idx)%profile%tau22)) &
    sem_box(sem_idx)%profile%tau22 = init_array(npts,zero) ! F2003 AUTO-REALLOC
    !
    if (.not. allocated(sem_box(sem_idx)%profile%tau31)) &
    sem_box(sem_idx)%profile%tau31 = init_array(npts,zero) ! F2003 AUTO-REALLOC
    !
    if (.not. allocated(sem_box(sem_idx)%profile%tau32)) &
    sem_box(sem_idx)%profile%tau32 = init_array(npts,zero) ! F2003 AUTO-REALLOC
    !
    if (.not. allocated(sem_box(sem_idx)%profile%tau33)) &
    sem_box(sem_idx)%profile%tau33 = init_array(npts,zero) ! F2003 AUTO-REALLOC
    !
    ! Allocate the components of the Cholesky decomposition of
    ! the prescribed Reynolds stress tensor for this profile.
    !
    sem_box(sem_idx)%profile%a11 = init_array(npts,zero) ! F2003 AUTO-REALLOC
    sem_box(sem_idx)%profile%a21 = init_array(npts,zero) ! F2003 AUTO-REALLOC
    sem_box(sem_idx)%profile%a22 = init_array(npts,zero) ! F2003 AUTO-REALLOC
    sem_box(sem_idx)%profile%a31 = init_array(npts,zero) ! F2003 AUTO-REALLOC
    sem_box(sem_idx)%profile%a32 = init_array(npts,zero) ! F2003 AUTO-REALLOC
    sem_box(sem_idx)%profile%a33 = init_array(npts,zero) ! F2003 AUTO-REALLOC
    !
    ! Compute the Reynolds stress tensor for the profile
    !
    do n = 1,npts
      !
      tau11 = sem_box(sem_idx)%profile%tau11(n)
      tau21 = sem_box(sem_idx)%profile%tau21(n)
      tau22 = sem_box(sem_idx)%profile%tau22(n)
      tau31 = sem_box(sem_idx)%profile%tau31(n)
      tau32 = sem_box(sem_idx)%profile%tau32(n)
      tau33 = sem_box(sem_idx)%profile%tau33(n)
      !
      a11 = sqrt(    tau11 )
     !a11 = sqrt(abs(tau11))
      !
      a21 = tau21/(a11+rndoff)
      a22 = sqrt( max( zero ,     tau22 -(a21**2) ) )
     !a22 = sqrt( max( zero , abs(tau22)-(a21**2) ) )
      !
      a31 = tau31/(a11+rndoff)
      a32 = (tau32-a21*a31)/(a22+rndoff)
      a33 = sqrt( max( zero ,     tau33 -(a31**2)-(a32**2) ) )
     !a33 = sqrt( max( zero , abs(tau33)-(a31**2)-(a32**2) ) )
      !
      sem_box(sem_idx)%profile%a11(n) = a11
      sem_box(sem_idx)%profile%a21(n) = a21
      sem_box(sem_idx)%profile%a22(n) = a22
      sem_box(sem_idx)%profile%a31(n) = a31
      sem_box(sem_idx)%profile%a32(n) = a32
      sem_box(sem_idx)%profile%a33(n) = a33
      !
    end do
    !
    ! If the averaged density was not included in the profile
    ! file, set the density for all points in the profile to
    ! the density specified for this boundary condition.
    !
    if (.not. allocated(sem_box(sem_idx)%profile%ravg)) then
      ravg = bc_in( sem_box(sem_idx)%bc_idx )%cv(nec)
      sem_box(sem_idx)%profile%ravg = init_array(npts,ravg) ! F2003 AUTO-REALLOC
    end if
    !
    ! If the averaged x-velocity was not included in the profile
    ! file, set the x-velocity for all points in the profile to
    ! the x-velocity specified for this boundary condition.
    !
    if (.not. allocated(sem_box(sem_idx)%profile%uavg)) then
      uavg = bc_in( sem_box(sem_idx)%bc_idx )%cv(nmx)
      sem_box(sem_idx)%profile%uavg = init_array(npts,uavg) ! F2003 AUTO-REALLOC
    end if
    !
    ! If the averaged y-velocity was not included in the profile
    ! file, set the y-velocity for all points in the profile to
    ! the y-velocity specified for this boundary condition.
    !
    if (.not. allocated(sem_box(sem_idx)%profile%vavg)) then
      vavg = bc_in( sem_box(sem_idx)%bc_idx )%cv(nmy)
      sem_box(sem_idx)%profile%vavg = init_array(npts,vavg) ! F2003 AUTO-REALLOC
    end if
    !
    ! If the averaged z-velocity was not included in the profile
    ! file, set the z-velocity for all points in the profile to
    ! the z-velocity specified for this boundary condition.
    !
    if (.not. allocated(sem_box(sem_idx)%profile%wavg)) then
      wavg = bc_in( sem_box(sem_idx)%bc_idx )%cv(nmz)
      sem_box(sem_idx)%profile%wavg = init_array(npts,wavg) ! F2003 AUTO-REALLOC
    end if
    !
    ! If the averaged temperature was not included in the profile
    ! file, set the temperature for all points in the profile to
    ! the temperature specified for this boundary condition.
    !
    if (.not. allocated(sem_box(sem_idx)%profile%tavg)) then
      tavg = bc_in( sem_box(sem_idx)%bc_idx )%gp(ntemp)
      sem_box(sem_idx)%profile%tavg = init_array(npts,tavg) ! F2003 AUTO-REALLOC
    end if
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_sem_box_profiles
!
!###############################################################################
!
subroutine read_sem_profile(prof,is_cylindrical)
  !
  !.. Formal Arguments ..
  type(sem_profile_t), intent(inout) :: prof
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(inout) :: is_cylindrical
  !
  !.. Local Scalars ..
  integer :: i,ierr,io,i1,i2,n,npts,nvar
  !
  logical(lk) :: r_instead_of_y
  !
  character(len=150) :: fname
  character(len=150) :: profile_variables
  !
  !.. Local Arrays ..
  integer :: temp_size(1:2)
  !
  !.. Local Allocatable Arrays ..
  integer,           allocatable :: int_char(:,:)
  real(wp),          allocatable :: temp_data(:,:)
  character(len=10), allocatable :: var_names(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_sem_profile"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! We need to create the profile data, so read the data from the profile files.
  !
  r_instead_of_y = fals
  !
  fname = trim(adjustl( prof%fname ))
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%  OPEN THE PROFILE FILE AND READ THE PROFILE VARIABLES  %%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  ! Have only the root MPI process read the profile file, and then
  ! broadcast the read data to all the other MPI processes. We want
  ! to do this to prevent file access contention for simulations
  ! involving a very large number of MPI processes.
  !
  if (mypnum == glb_root) then
    !
    open (newunit=io,file=fname,status="old",action="read", &
                     position="rewind",iostat=ierr,iomsg=error_message)
    call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Read the number of points in this profile file
    ! and allocate all the profile data arrays
    !
    read (io,'(a)') profile_variables
    !
    call identify_profile_variables(profile_variables,var_names)
    !
    nvar = size(var_names)
    !
    read (io,*) npts
    !
    prof%npts = npts
    !
    ! Read in all the profile data
    !
    temp_data = init_array(nvar,npts,zero) ! F2003 AUTO-REALLOC
    !
    do n = 1,npts
      read (io,*) (temp_data(i,n), i=1,nvar)
    end do
    !
    ! Close the profile file for this SEM box
    !
    close (io,iostat=ierr,iomsg=error_message)
    call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
    !
    temp_size(1) = nvar
    temp_size(2) = npts
    !
  end if
  !
  if (ncpu > 1) then
    !
    call mpi_bcast(temp_size,size(temp_size,kind=int_mpi),mpi_inttyp, &
                   glb_root,MPI_COMM_WORLD,mpierr)
    !
    nvar = temp_size(1)
    npts = temp_size(2)
    !
    prof%npts = npts
    !
    if (.not. allocated(temp_data)) then
      temp_data = init_array(nvar,npts,zero) ! F2003 AUTO-REALLOC
    end if
    !
    call mpi_bcast(temp_data,size(temp_data,kind=int_mpi),mpi_flttyp, &
                   glb_root,MPI_COMM_WORLD,mpierr)
    !
    ! Broadcast var_names to all the non-root MPI processes
    !
    if (.not. allocated(var_names)) then
      allocate ( var_names(1:nvar) , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"var_names",1,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    int_char = init_array( (len(var_names)/4)+1 , nvar , 0) ! F2003 AUTO-REALLOC
    !
    if (mypnum == glb_root) then
      do n = 1,nvar
        do i = 1,size(int_char,dim=1)
          i1 = (i-1)*4 + 1
          i2 = min(i*4,len(var_names))
          int_char(i,n) = transfer( var_names(n)(i1:i2) , int_char(i,n) )
        end do
      end do
    end if
    !
    call mpi_bcast(int_char,size(int_char,kind=int_mpi),mpi_inttyp, &
                   glb_root,MPI_COMM_WORLD,mpierr)
    !
    if (mypnum /= glb_root) then
      do n = 1,nvar
        do i = 1,size(int_char,dim=1)
          i1 = (i-1)*4 + 1
          i2 = min(i*4,len(var_names))
          var_names(n)(i1:i2) = transfer( int_char(i,n) , var_names(n)(i1:i2) )
        end do
      end do
    end if
    !
  end if
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%  NOW SAVE THE READ DATA INTO THE PROFILE DERIVED TYPE  %%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  do i = 1,nvar
    select case (uppercase(trim(adjustl(var_names(i)))))
      !
      case ("Y")
        prof%y = temp_data(i,1:npts)      ! F2003 AUTO-REALLOC
        !
      case ("R")
        r_instead_of_y = true
        prof%y = temp_data(i,1:npts)      ! F2003 AUTO-REALLOC
        !
      case ("RHO","DENSITY")
        prof%ravg = temp_data(i,1:npts)   ! F2003 AUTO-REALLOC
      case ("VX","U")
        prof%uavg = temp_data(i,1:npts)   ! F2003 AUTO-REALLOC
      case ("VY","V")
        prof%vavg = temp_data(i,1:npts)   ! F2003 AUTO-REALLOC
      case ("VZ","W")
        prof%wavg = temp_data(i,1:npts)   ! F2003 AUTO-REALLOC
      case ("T","TEMPERATURE")
        prof%tavg = temp_data(i,1:npts)   ! F2003 AUTO-REALLOC
        !
      case ("TAUXX","TAU11","RXX","R11")
        prof%tau11 = temp_data(i,1:npts)  ! F2003 AUTO-REALLOC
      case ("TAUYX","TAUXY","TAU21","TAU12","RYX","RXY","R21","R12")
        prof%tau21 = temp_data(i,1:npts)  ! F2003 AUTO-REALLOC
      case ("TAUYY","TAU22","RYY","R22")
        prof%tau22 = temp_data(i,1:npts)  ! F2003 AUTO-REALLOC
      case ("TAUZX","TAUXZ","TAU31","TAU13","RZX","RXZ","R31","R13")
        prof%tau31 = temp_data(i,1:npts)  ! F2003 AUTO-REALLOC
      case ("TAUZY","TAUYZ","TAU32","TAU23","RZY","RYZ","R32","R23")
        prof%tau32 = temp_data(i,1:npts)  ! F2003 AUTO-REALLOC
      case ("TAUZZ","TAU33","RZZ","R33")
        prof%tau33 = temp_data(i,1:npts)  ! F2003 AUTO-REALLOC
        !
    end select
    !
  end do
  !
  deallocate ( temp_data , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_data",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Nondimensionalize all the profile quantities
  !
  call nondimensionalize_sem_profile(prof)
  !
  ! If is_cylindrical is present, return true if the variable R
  ! was found in the profile file instead of variable Y.
  !
  ! NOTE: We do NOT want to return the value of r_instead_of_y to is_cylindrical
  !       because is_cylindrical could already be true through other means of
  !       detecting a cylindrical based SEM box. We do not want to change this
  !       to false because the variable Y was in the profile file instead of the
  !       variable R. In these situations, Y is just assumed to be R.
  !
  if (present(is_cylindrical)) then
    if (r_instead_of_y) then
      is_cylindrical = true
    end if
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_sem_profile
!
!###############################################################################
!
pure &
subroutine nondimensionalize_sem_profile(this)
  !
  !.. Use Statements ..
  use ovar, only : rhoref,aref,tref
  !
  !.. Formal Arguments ..
  type(sem_profile_t), intent(inout) :: this
  !
  !.. Local Scalars ..
  real(wp) :: rhoref_inv,aref_inv,tref_inv,aref2_inv
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "nondimensionalize_sem_profile"
  !
continue
  !
  rhoref_inv = one / rhoref
  aref_inv   = one / aref
  tref_inv   = one / tref
  !
  aref2_inv = aref_inv * aref_inv
  !
  if (allocated(this%ravg)) this%ravg(:) = this%ravg(:) * rhoref_inv
  if (allocated(this%uavg)) this%uavg(:) = this%uavg(:) * aref_inv
  if (allocated(this%vavg)) this%vavg(:) = this%vavg(:) * aref_inv
  if (allocated(this%wavg)) this%wavg(:) = this%wavg(:) * aref_inv
  if (allocated(this%tavg)) this%tavg(:) = this%tavg(:) * tref_inv
  !
  if (allocated(this%tau11)) this%tau11(:) = this%tau11(:) * aref2_inv
  if (allocated(this%tau21)) this%tau21(:) = this%tau21(:) * aref2_inv
  if (allocated(this%tau22)) this%tau22(:) = this%tau22(:) * aref2_inv
  if (allocated(this%tau31)) this%tau31(:) = this%tau31(:) * aref2_inv
  if (allocated(this%tau32)) this%tau32(:) = this%tau32(:) * aref2_inv
  if (allocated(this%tau33)) this%tau33(:) = this%tau33(:) * aref2_inv
  !
end subroutine nondimensionalize_sem_profile
!
!###############################################################################
!
subroutine identify_profile_variables(text,variable_names)
  !
  ! Routine to parse a character string to extract a list variable names
  !
  ! TEXT           : character string to be parsed
  ! VARIABLE_NAMES : list of variable names
  !
  !.. Formal Arguments ..
  character(len=*),               intent(in) :: text
  character(len=*), allocatable, intent(out) :: variable_names(:)
  !
  !.. Local Scalars ..
  integer :: ibeg,ierr,is,m,n,nn,nvar
  !
  character(len=len(text)) :: strng
  !
  !.. Local Allocatable Arrays ..
  character(len=len(variable_names)), allocatable :: clist(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: delimiters(1:4) = ["'",'"',",",achar(9)]
  character(len=*), parameter :: pname = "identify_profile_variables"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  strng = trim(adjustl(text))
  !
  ! Replace all single quotes, double quotes, commas, and tabs with spaces
  !
  do n = 1,len_trim(strng)
    if (any(strng(n:n) == delimiters)) then
      strng(n:n) = " "
    end if
  end do
  !
  m = 0
  !
  if (len_trim(strng) > 0) then
    is = 1
    find_next_word: do
      !
      if (strng(is:is) /= " ") then
        !
        ibeg = is
        !
        find_end_of_word: do
          if (is == len_trim(strng)) then
           !nn = len(clist) - len(strng(ibeg:is))
           !if (.not. allocated(clist)) then
           !  ! F2003 AUTO-REALLOC
           !  clist = [ strng(ibeg:is) // repeat(" ",ncopies=nn) ]
           !else
           !  ! F2003 AUTO-REALLOC
           !  clist = [ clist , strng(ibeg:is) // repeat(" ",ncopies=nn) ]
           !end if
            call reallocate( clist , m+1 )
            clist(m+1) = strng(ibeg:is)
          else if (strng(is:is) == " ") then
           !nn = len(clist) - len(strng(ibeg:is-1))
           !if (.not. allocated(clist)) then
           !  ! F2003 AUTO-REALLOC
           !  clist = [ strng(ibeg:is-1) // repeat(" ",ncopies=nn) ]
           !else
           !  ! F2003 AUTO-REALLOC
           !  clist = [ clist , strng(ibeg:is-1) // repeat(" ",ncopies=nn) ]
           !end if
            m = m + 1
            call reallocate( clist , m )
            clist(m) = strng(ibeg:is-1)
            ibeg = is
            exit find_end_of_word
          end if
          is = is + 1
        end do find_end_of_word
        !
      else
        !
        is = is + 1
        !
      end if
      !
      if (is > len_trim(strng)) exit find_next_word
      !
    end do find_next_word
  end if
  !
  ! Get the final number of sub-strings found in strng
  !
  nvar = size(clist)
  !
  allocate ( variable_names(1:nvar) , stat=ierr )
  !
  ! Make sure that everything in clist is adjusted left and trimmed.
  ! Probably not necessary but the extra overhead is basically zero
  ! and worth the extra piece of mind.
  !
  do n = 1,nvar
    variable_names(n) = trim(adjustl(clist(n)))
  end do
  !
  deallocate ( clist , stat=ierr )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine identify_profile_variables
!
!###############################################################################
!
subroutine check_sem_profile_allocations(profile_data_is_allocated)
  !
  !.. Formal Arguments ..
  logical(lk), intent(inout) :: profile_data_is_allocated
  !
  !.. Local Scalars ..
  integer :: num_allocs,num_checks,sem_idx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_sem_profile_allocations"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  profile_data_is_allocated = fals
  !
  num_allocs = 0
  num_checks = 0
  !
  ! Check the allocation status of all the profile arrays
  ! for each of the SEM boxes.
  !
  do sem_idx = 1,num_sem_boxes
    !
    ! Make sure that sem_box(sem_idx)%profile is allocated for this SEM box.
    !
    if (.not. allocated(sem_box(sem_idx)%profile)) then
      allocate ( sem_box(sem_idx)%profile )
    end if
    !
    ! Now check the allocation status of all the
    ! component arrays of sem_box(sem_idx)%profile.
    !
    call check_array_alloc(sem_box(sem_idx)%profile%y, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    !
    call check_array_alloc(sem_box(sem_idx)%profile%ravg, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    call check_array_alloc(sem_box(sem_idx)%profile%uavg, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    call check_array_alloc(sem_box(sem_idx)%profile%vavg, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    call check_array_alloc(sem_box(sem_idx)%profile%wavg, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    call check_array_alloc(sem_box(sem_idx)%profile%tavg, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    !
    call check_array_alloc(sem_box(sem_idx)%profile%tau11, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    call check_array_alloc(sem_box(sem_idx)%profile%tau21, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    call check_array_alloc(sem_box(sem_idx)%profile%tau22, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    call check_array_alloc(sem_box(sem_idx)%profile%tau31, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    call check_array_alloc(sem_box(sem_idx)%profile%tau32, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    call check_array_alloc(sem_box(sem_idx)%profile%tau33, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    !
    call check_array_alloc(sem_box(sem_idx)%profile%a11, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    call check_array_alloc(sem_box(sem_idx)%profile%a21, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    call check_array_alloc(sem_box(sem_idx)%profile%a22, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    call check_array_alloc(sem_box(sem_idx)%profile%a31, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    call check_array_alloc(sem_box(sem_idx)%profile%a32, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    call check_array_alloc(sem_box(sem_idx)%profile%a33, &
                           profile_data_is_allocated, &
                           num_checks,num_allocs)
    !
  end do
  !
  ! Check to make sure that either none or all of the profile data arrays
  ! have been allocated.
  !
  if (profile_data_is_allocated) then
    if (num_allocs /= num_checks) then
      write (error_message,1)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("ERROR!!! CHECKING THE SEM PROFILE DATA BEFORE ATTEMPTING ", &
            "READ THE PROFILE FILES HAS FOUND THAT SOME, BUT NOT ALL, ", &
            "OF THE PROFILE DATA ARRAYS HAVE BEEN ALLOCATED. IT SHOULD ", &
            "BE THAT EITHER NONE OR ALL OF THE PROFILE DATA ARRAYS ARE ", &
            "ALLOCATED. SOMETHING IS NOT CORRECT HERE.")
  !
end subroutine check_sem_profile_allocations
!
!###############################################################################
!
pure &
subroutine check_array_alloc(array,profile_data_is_allocated, &
                             num_checks,num_allocs)
  !
  !.. Formal Arguments ..
  real(wp), allocatable, intent(in) :: array(:)
  logical(lk),        intent(inout) :: profile_data_is_allocated
  integer,            intent(inout) :: num_checks
  integer,            intent(inout) :: num_allocs
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_array_alloc"
  !
continue
  !
  num_checks = num_checks + 1
  !
  if (allocated(array)) then
    profile_data_is_allocated = true
    num_allocs = num_allocs + 1
  end if
  !
end subroutine check_array_alloc
!
!###############################################################################
!
subroutine output_details_of_sem_boxes()
  !
  !.. Use Statements ..
  use ovar, only : aref
  use ovar, only : grid_scaling_factor
  !
  !.. Local Scalars ..
  integer  :: sem_idx
  real(wp) :: gsf
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_details_of_sem_boxes"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  gsf = one
  if (abs(grid_scaling_factor) > eps6) gsf = one/grid_scaling_factor 
  !
  if (mypnum == glb_root) then
    do sem_idx = 1,num_sem_boxes
      write (iout,*)
      write (iout,1) sem_idx
      write (iout,*) 'Uconv     = ',sem_box(sem_idx)%convective_velocity*aref
      write (iout,*) 'Lturb     = ',sem_box(sem_idx)%length_scale*gsf
      write (iout,*) 'Neddies   = ',sem_box(sem_idx)%num_eddies
      write (iout,*) 'Volbox    = ',sem_box(sem_idx)%box_volume*gsf*gsf*gsf
      write (iout,*) 'Min x     = ',sem_box(sem_idx)%boxmin(1)*gsf
      write (iout,*) 'Max x     = ',sem_box(sem_idx)%boxmax(1)*gsf
      write (iout,*) 'Min y     = ',sem_box(sem_idx)%boxmin(2)*gsf
      write (iout,*) 'Max y     = ',sem_box(sem_idx)%boxmax(2)*gsf
      write (iout,*) 'Min z     = ',sem_box(sem_idx)%boxmin(3)*gsf
      write (iout,*) 'Max z     = ',sem_box(sem_idx)%boxmax(3)*gsf
      if (sem_box(sem_idx)%is_cylindrical) then
        write (iout,*) 'Min r     = ',sem_box(sem_idx)%boxmin(4)*gsf
        write (iout,*) 'Max r     = ',sem_box(sem_idx)%boxmax(4)*gsf
        write (iout,*) 'Min theta = ',sem_box(sem_idx)%boxmin(5)*radians2degrees
        write (iout,*) 'Max theta = ',sem_box(sem_idx)%boxmax(5)*radians2degrees
      end if
    end do
    write (iout,*)
  end if
  flush (iout)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format(1x,'Synthetic Eddy Method Initialization, SEM Box = ',i4,/, &
           1x,'----------------------------------------------------')
  !
end subroutine output_details_of_sem_boxes
!
!###############################################################################
!
pure &
function tentfun(dl,sigma) result(return_value)
  !
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Simplified tent function for Jarrin's turb. inflow boundary
  ! NOTE: This is a simplified version of the original function:
  !         - the sqrt(3/2) constant has been pulled out for computational
  !           efficiency and is multiplied to the result after performing
  !           the summation of the tent functions over all the eddies
  !         - The original functions single argument of
  !              ( x(face-point-idx) - x(eddy-idx) ) / length-scale
  !           has been replaced with two arguments of
  !              abs( x(face-point-idx) - x(eddy-idx) ), and
  !              length-scale
  !           because this simplifies the conditional logic within the
  !           function as well as the function result.
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  !    DL : abs[ x(point on SEM plane) - x(location of synthetic eddy) ]
  ! SIGMA : length scale input
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: dl
  real(wp), intent(in) :: sigma
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
continue
  !
  if (dl < sigma) then
    return_value = sigma - dl
  else
    return_value = zero
  end if
  !
end function tentfun
!
!###############################################################################
!
pure &
function old_tentfun(arg) result(return_value)
  !
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Original tent function for Jarrin's turb. inflow boundary
  ! NOTE: The single input argument should be the value:
  !        ( x(face-point-idx) - x(eddy-idx) ) / length-scale
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  ! ARG : [x(point on SEM plane) - x(loc. of synthetic eddy)] / length-scale
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: arg
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
continue
  !
  if (abs(arg) < one) then
    return_value = sqrt(three/two)*(one-abs(arg))
  else
    return_value = zero
  end if
  !
end function old_tentfun
!
!###############################################################################
!
subroutine get_sem_profile_heights()
  !
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  ! Extract the maximum heights from the profiles
  ! associated with each synthetic eddy method boundary
  !
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  !.. Use Statements ..
  use ovar, only : bc_input
  !
  !.. Local Scalars ..
  integer  :: i,io,ierr,n,npts,bc_idx,iskip
  real(wp) :: ymax,vskip
  !
  character(len=100) :: fname
  character(len=150) :: profile_variables
  !
  !.. Local Allocatable Arrays ..
  real(wp),          allocatable :: y(:)
  character(len=10), allocatable :: var_names(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_sem_profile_heights"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Only have the global root processor extract the maximum heights
  !
  if (mypnum == glb_root) then
    !
    ! Find the SEM boundary conditions and extract
    ! the maximum heights of each SEM profile.
    !
    do bc_idx = 1,size(bc_input)
      !
      if (bc_input(bc_idx)%set_bc_value == bc_sem) then
        !
        fname = bc_input(bc_idx)%sem%profile_file
        !
        ! Open the profile file for this SEM boundary
        !
        open (newunit=io,file=trim(adjustl(fname)),status="old",action="read", &
                         position="rewind",iostat=ierr,iomsg=error_message)
        call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
        !
        ! Read the name of the profile variables
        !
        read (io,'(a)') profile_variables
        !
        ! Identify the variable names
        !
        call identify_profile_variables(profile_variables,var_names)
        !
        ! Count the number of variables to skip
        ! to get to the y/wall-distance variable
        !
        iskip = 0
        find_iskip_loop: &
        do i = 1,size(var_names)
          if (uppercase(trim(adjustl(var_names(i)))) == "Y") then
            exit find_iskip_loop
          else if (uppercase(trim(adjustl(var_names(i)))) == "R") then
            ! Mark this SEM box as cylindrical if the variable name is "R"
            bc_input(bc_idx)%sem%is_cylindrical = true
            exit find_iskip_loop
          end if
          iskip = iskip + 1
        end do find_iskip_loop
        !
        ! Read in the number of data points in this SEM profile
        !
        read (io,*) npts
        !
        ! Allocate the array for reading in the wall distances
        !
        y = init_array(npts,zero) ! F2003 AUTO-REALLOC
        !
        ! Read the wall distance of each SEM profile data point
        !
        do n = 1,npts
          read (io,*) (vskip, i=1,iskip), y(n)
        end do
        !
        ! Close the profile file for this SEM boundary
        !
        close (io,iostat=ierr,iomsg=error_message)
        call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
        !
        ! Compute the maximum wall distance for this profile
        !
        ymax = maxval(y(1:npts))
        !
        ! Save the maximum wall distance into bc_input
        !
        bc_input(bc_idx)%sem%profile_height = ymax
        !
      end if
      !
    end do
    !
    ! Save the maximum wall distances into the y array so they can be
    ! broadcast to all the other processors
    !
    ! F2003 AUTO-REALLOC
    y = [ bc_input(:)%sem%profile_height ]
    !
  else
    !
    ! Allocate y to the combined size of bc_input so that the
    ! non-global-root processes can receive the broadcast of the maximum
    ! wall distances of each SEM boundary profile.
    !
    npts = size(bc_input)
    !
    y = init_array(npts,zero) ! F2003 AUTO-REALLOC
    !
  end if
  !
  if (ncpu > 1) then
    call mpi_bcast(y,size(y,kind=int_mpi),mpi_flttyp, &
                   glb_root,MPI_COMM_WORLD,mpierr)
  end if
  !
  ! All processes copy the maximum wall distances back into bc_input
  !
  do bc_idx = 1,size(bc_input)
    bc_input(bc_idx)%sem%profile_height = y(bc_idx)
  end do
  !
  ! Deallocate y before we leave
  !
  deallocate ( y , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"y",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_sem_profile_heights
!
!###############################################################################
!
subroutine separate_sem_boundary_faces()
  !
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  ! Separate any boundary faces included with an SEM boundary that do not
  ! have a face node within the SEM profile height, and make those boundary
  ! faces a part of a new subsonic inflow boundary condition.
  !
  ! The following need to be defined/populated before calling this procedure
  !   - bface(1,:)
  !   - bface(5,:)
  !   - grid%xyz(:,:)
  !   - grid%elem(ncell+1:ncell+nfbnd)%nodes(:)
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  !.. Use Statements ..
  use geovar, only : nfbnd,ncell,nnode
  use geovar, only : bface,grid,xyz_nodes
  use ovar,   only : bc_input,bc_input_t,bc_names
  use ovar,   only : output_bface_array
  !
  !.. Local Scalars ..
  integer  :: i,ierr,n,nc,nf,nn,npts,bc_idx,grp_idx
  integer  :: ibc,nbc,ibcn,nbcn
  integer  :: n_sem_pts,n_sem_bc
  real(wp) :: profile_height
  !
  character(len=32) :: cbcn
  character(len=32) :: cbci
  !
  !.. Local Allocatable Arrays ..
  real(wp),     allocatable :: xyz_wall(:,:)
  real(wp),     allocatable :: xyz_sem(:,:)
  type(wall_t), allocatable :: wall(:)
  !
  logical(lk), allocatable :: msk(:)
  !
  integer, allocatable :: sem2node(:)
  integer, allocatable :: node2sem(:)
  !
  type(bc_input_t), allocatable :: tmp_bc_input(:)
  !
  character(len=len(bc_names)), allocatable :: tmp_bc_names(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "separate_sem_boundary_faces"
  !
continue
  !
  if (output_bface_array) call output_bface
  !
  ! Make sure that there are SEM boundary conditions. If not,
  ! dont do anything and just return to the calling procedure.
  !
  if (all(bface(1,:) /= bc_sem)) return
  !
  call debug_timer(entering_procedure,pname)
  !
  ! First, read all the SEM profile files and extract the maximum profile
  ! height from each SEM profile.
  !
  call get_sem_profile_heights
  !
  ! Find all the boundary nodes that are on wall boundary conditions
  !
  if (i_am_host_root) then
    !
    msk = init_array(nnode,fals) ! F2003 AUTO-REALLOC
    !
    do nf = 1,nfbnd
      if (any(bface(1,nf) == bc_physical_walls)) then
        msk( grid%elem(ncell+nf)%nodes(:) ) = true
      end if
    end do
    !
    ! F2003 AUTO-REALLOC
    xyz_wall = xyz_nodes( 1:3 , pack(intseq(1,nnode),mask=msk(1:nnode)) )
    !
    !
    ! Reinitialize the msk array to false
    !
    msk = init_array(nnode,fals) ! F2003 AUTO-REALLOC
    !
    do nf = 1,nfbnd
      if (bface(1,nf) == bc_sem) then
        msk( grid%elem(ncell+nf)%nodes(:) ) = true
      end if
    end do
    !
    n_sem_pts = count( msk(1:nnode) )
    !
    sem2node = pack( intseq(1,nnode) , mask=msk(1:nnode) ) ! F2003 AUTO-REALLOC
    node2sem = unpack( intseq(1,n_sem_pts) , &             ! F2003 AUTO-REALLOC
                       mask=msk(1:nnode) , field=0 )
    !
    ! F2003 AUTO-REALLOC
    xyz_sem = xyz_nodes( 1:3 , sem2node(1:n_sem_pts) )
    !
    ! Compute the minimum wall distance for all
    ! the boundary nodes on SEM boundaries
    !
    call brute_force_wall_distance(xyz_wall,xyz_sem,wall)
    !
  end if
  !
  ! Deallocate msk, xyz_wall, sem2node, xyz_sem since they are no longer needed
  !
  if (allocated(msk)) then
    deallocate ( msk , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(xyz_wall)) then
    deallocate ( xyz_wall , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_wall",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(sem2node)) then
    deallocate ( sem2node , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"sem2node",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(xyz_sem)) then
    deallocate ( xyz_sem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_sem",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Reallocate bc_input to create additional entries for the subsonic
  ! inflow boundary conditions that are separated from each SEM boundary.
  !
  n_sem_bc = 0
  do bc_idx = 1,size(bc_input)
    if (len_trim(adjustl(bc_input(bc_idx)%sem%profile_file)) > 0) then
      n_sem_bc = n_sem_bc + 1
    end if
  end do
  !
  nbc = size(bc_input)
  !
  allocate ( tmp_bc_input(1:nbc+n_sem_bc) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_bc_input",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Copy bc_input into tmp_bc_input
  !
  do bc_idx = 1,size(bc_input)
    tmp_bc_input(bc_idx) = bc_input(bc_idx)
  end do
  !
  nbcn = size(bc_names)
  !
  allocate ( tmp_bc_names(1:nbcn+n_sem_bc) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_bc_names",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Copy bc_names into tmp_bc_names
  !
  do bc_idx = 1,size(bc_names)
    tmp_bc_names(bc_idx) = bc_names(bc_idx)
  end do
  !
  ! Create the counters for the new elements of bc_input and bc_names
  !
  ibc = nbc
  ibcn = nbcn
  !
  ! Create the separate subsonic inflow boundary condition for each SEM boundary
  !
  do bc_idx = 1,size(bc_input)
    if (len_trim(adjustl(bc_input(bc_idx)%sem%profile_file)) > 0) then
      !
      ibc = ibc + 1
      ibcn = ibcn + 1
      !
      tmp_bc_input(ibc) = bc_input(bc_idx)
      tmp_bc_input(ibc)%set_bc_value = bc_sub_inflow
      tmp_bc_input(ibc)%name = trim(adjustl(bc_input(bc_idx)%name)) // &
                               " (Subsonic Inflow)"
      !
      tmp_bc_input(bc_idx)%sem%sub_in_bc_idx = ibc
      tmp_bc_input(bc_idx)%sem%sub_in_bc_grp = ibcn
      !
      tmp_bc_names(ibcn) = trim(adjustl(tmp_bc_input(ibc)%name))
      !
    end if
  end do
  !
  ! Transfer tmp_bc_input to bc_input
  !
  call move_alloc( from=tmp_bc_input , to=bc_input )
  !
  ! Transfer tmp_bc_names to bc_names
  !
  call move_alloc( from=tmp_bc_names , to=bc_names )
  !
  ! Loop through the boundary faces to find all the boundary faces that are
  ! on an SEM boundary. Check the wall distance of every node for each SEM
  ! boundary face, and change the boundary condition for the boundary face
  ! to a subsonic inflow boundary condition if all the nodes are further from
  ! the nearest wall than the profile height of the SEM boundary profile.
  !
 !do bc_idx = 1,size(bc_input)
 !  write (iout,91) bc_idx,bc_input(bc_idx)%sem%profile_height
 !end do
 !write (iout,*)
 !91 format ("bc_input(",i2,")%sem%profile_height = ",es14.6)
  !
  if (i_am_host_root) then
    !
    bnd_face_loop: &
    do nf = 1,nfbnd
      !
      if (bface(1,nf) /= bc_sem) cycle bnd_face_loop
      !
      nc = ncell + nf
      !
      bc_idx = bface(6,nf)
      !
      profile_height = bc_input(bc_idx)%sem%profile_height
     !write (iout,92) nf,bc_idx,profile_height
     !92 format ("nf = ",i4,"; bc_idx = ",i2,"; profile_height = ",es14.6)
      !
      one_pass_loop: &
      do i = 1,1
        !
        node_loop: &
        do n = 1,size(grid%elem(nc)%nodes)
          !
          ! Get the index for the current node of this boundary face
          ! relative to all other nodes on SEM boundary faces,
          ! i.e., the index of this node within the array wall
          !
          nn = node2sem( grid%elem(nc)%nodes(n) )
          !
          ! If the distance to the nearest wall of the current node is
          ! within the profile height of the SEM boundary profile, exit
          ! from the one-pass-loop so that the boundary conditions for
          ! this boundary face are unchanged and remain an SEM BC.
          !
          if (wall(nn)%distance <= profile_height) exit one_pass_loop
          !
        end do node_loop
        !
        ! If the loop over the nodes completed without triggering an exit
        ! from the one-pass-loop, this boundary face is outside the profile
        ! height of the SEM boundary, so adjust the boundary conditions for
        ! this boundary face to a subsonic inflow boundary.
        !
        bface(1,nf) = bc_sub_inflow
        bface(5,nf) = bc_input(bc_idx)%sem%sub_in_bc_grp
        bface(6,nf) = bc_input(bc_idx)%sem%sub_in_bc_idx
        !
        if (allocated(grid%elem(nc)%tags)) then
          grid%elem(nc)%tags(1) = bc_sub_inflow
        else
          grid%elem(nc)%tags = [ bc_sub_inflow ] ! F2003 AUTO-REALLOC
        end if
        !
      end do one_pass_loop
      !
    end do bnd_face_loop
    !
  end if
  !
  ! Broadcast the updated bface array to all the non-host-root MPI processes
  ! on each host.
  !
  if (ncpu > 1) then
    call mpi_bcast(bface,size(bface,kind=int_mpi),mpi_inttyp, &
                   host_root,my_host_comm,mpierr)
  end if
  !
  !#######################################################################
  !#######################################################################
  !
  ! MIGHT NEED TO ADD CODE TO MAKE SURE THAT EACH NEW SUBSONIC BC ELEMENT
  ! OF BC_INPUT AND BC_NAMES HAS AT LEAST ONE BOUNDARY FACE ASSOCIATED
  ! WITH IT. IF THERE ARE NO BOUNDARY FACES ASSOCIATED WITH IT, DELETE
  ! THAT ELEMENT FROM BC_INPUT AND BC_NAMES, SHRINK THE ARRAYS, AND MAKE
  ! THE NECESSARY ADJUSTMENTS TO BFACE, BC_INPUT(:)%SEM%SUB_IN_BC_IDX,
  ! BC_INPUT(:)%SEM%SUB_IN_BC_GRP, ETC.
  !
  ! HOWEVER, I THINK FOR NOW, JUST LEAVING THESE "EMPTY" ELEMENTS OF
  ! BC_INPUT AND BC_NAMES (I.E., NO BOUNDARY FACES ARE ASSOCIATED WITH
  ! THE ELEMENT) WILL NOT CREATE ANY PROBLEMS, AND THESE ELEMENTS WILL
  ! JUST TAKE UP A NEGLIGIBLE AMOUNT OF MEMORY.
  !
 !tmp_bc_input = bc_input ! F2003 AUTO-REALLOC
 !!
 !bc_input_check: &
 !do ibc = 1,nbc
 !  !
 !  if (len_trim(adjustl(bc_input(bc_idx)%sem%profile_file)) > 0) then
 !    !
 !    grp_idx = bc_input(ibc)%sem%sub_in_bc_grp
 !    bc_idx  = bc_input(ibc)%sem%sub_in_bc_grp
 !    !
 !    do nf = 1,nfbnd
 !      if (bface(1,nf) == bc_sub_inflow) then
 !        if (bface(5,nf) == grp_idx) then
 !          if (bface(6,nf) == bc_idx) then
 !            cycle bc_input_check
 !          end if
 !        end if
 !      end if
 !    end do
 !    !
 !    ! Now do something if the above loop did not find any boundary faces
 !    ! associated with the subsonic inflow group created for this SEM
 !    ! boundary condition.
 !    !
 !    ...
 !    !
 !  end if
 !  !
 !end do bc_input_check
 !!
 !call move_alloc( from=tmp_bc_input , to=bc_input )
 !!
 !tmp_bc_names = bc_names ! F2003 AUTO-REALLOC
 !!
 !bc_names_check: &
 !do ibcn = nbcn+1,size(bc_names)
 !  !
 !  cbcn = uppercase( trim(adjustl(bc_names(ibcn))) )
 !  !
 !  do bc_idx = 1,size(bc_input)
 !    !
 !    cbci = uppercase( trim(adjustl(bc_input(bc_idx)%name)) )
 !    !
 !    if (cbci == cbcn) then
 !      cycle bc_names_check
 !    end if
 !    !
 !  end do
 !  !
 !  ! Now do something if the above loop did not find a matching
 !  ! bc_input(:)%name to the current element of bc_names
 !  !
 !  ...
 !  !
 !end do
 !!
 !call move_alloc( from=tmp_bc_names , to=bc_names )
  !
  !#######################################################################
  !#######################################################################
  !
  ! Deallocate local arrays before leaving if still allocated
  !
  if (allocated(wall)) then
    deallocate ( wall , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"wall",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(node2sem)) then
    deallocate ( node2sem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"node2sem",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine separate_sem_boundary_faces
!
!###############################################################################
!
subroutine output_sem_averages(output_inflow_profile)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: output_inflow_profile
  !
  !.. Local Scalars ..
  integer :: i,ierr,n,nvar,num_fp,num_total_sem_pts,sem_idx
  !
  logical(lk) :: output_interpolated_inflow_profile
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: num_sem_pts(:,:)
  integer,  allocatable :: num_box_sem_pts(:)
  integer,  allocatable :: num_my_sem_pts(:)
  real(wp), allocatable :: my_ave_data(:,:)
  real(wp), allocatable :: all_ave_data(:,:)
  !
  integer(int_mpi), allocatable :: counts(:)
  integer(int_mpi), allocatable :: offset(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_sem_averages"
  !
continue
  !
  if (.not. allocated(sem_avg)) return
  !
  call debug_timer(entering_procedure,pname)
  !
  output_interpolated_inflow_profile = fals
  if (present(output_inflow_profile)) then
    output_interpolated_inflow_profile = output_inflow_profile
  end if

  !
  ! Have each MPI process count the number of SEM face points it has
  ! in each SEM box and save that to its location within num_sem_pts.
  !
  num_sem_pts = init_array(num_sem_boxes,ncpu,0) ! F2003 AUTO-REALLOC
  !
  do sem_idx = 1,num_sem_boxes
    num_sem_pts(sem_idx,mypnum+1) = my_sem_pts(sem_idx)%num_fp
  end do
  !
  ! Now gather the number of SEM face points for all MPI processes
  !
  if (ncpu > 1) then
    call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_DATATYPE_NULL, &
                       num_sem_pts,int(num_sem_boxes,kind=INT_MPI), &
                       mpi_inttyp,MPI_COMM_WORLD,mpierr)
  end if
  !
!#define DEBUG_SEM
#ifdef DEBUG_SEM
  if (mypnum == glb_root) then
    write (iout,*)
  end if
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  do i = 0,ncpu-1
    if (mypnum == i) then
      do sem_idx = 1,num_sem_boxes
        do n = 1,ncpu
          write (iout,1) mypnum,sem_idx,n-1,num_sem_pts(sem_idx,n)
        end do
      end do
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end do
  1 format (2x,"CPU #",i2," : num_sem_pts(SEM BOX #",i0,",CPU #",i0,") = ",i0)
#endif
  !
  ! Count the total number of SEM face points
  ! in each box across all MPI processes.
  !
  num_box_sem_pts = sum(num_sem_pts,dim=2) ! F2003 AUTO-REALLOC
  !
#ifdef DEBUG_SEM
  if (mypnum == glb_root) then
    write (iout,*)
  end if
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  do i = 0,ncpu-1
    if (mypnum == i) then
      do sem_idx = 1,num_sem_boxes
        write (iout,2) mypnum,sem_idx,num_box_sem_pts(sem_idx)
      end do
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end do
  2 format (2x,"CPU #",i2," : num_box_sem_pts(SEM BOX #",i0,") = ",i0)
#endif
  !
  ! Count the total number of SEM face points
  ! across all SEM boxes for each MPI process
  !
  num_my_sem_pts = sum(num_sem_pts,dim=1) ! F2003 AUTO-REALLOC
  !
#ifdef DEBUG_SEM
  if (mypnum == glb_root) then
    write (iout,*)
  end if
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  do i = 0,ncpu-1
    if (mypnum == i) then
      do n = 1,ncpu
        write (iout,3) mypnum,n-1,num_my_sem_pts(n)
      end do
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end do
  3 format (2x,"CPU #",i2," : num_my_sem_pts(CPU #",i0,") = ",i0)
#endif
  !
  ! Count the total number of SEM face points
  ! across all SEM boxes and all MPI processes
  !
  num_total_sem_pts = sum(num_sem_pts)
  !
#ifdef DEBUG_SEM
  if (mypnum == glb_root) then
    write (iout,*)
  end if
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  do i = 0,ncpu-1
    if (mypnum == i) then
      write (iout,4) mypnum,num_total_sem_pts
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end do
  4 format (2x,"CPU #",i2," : num_total_sem_pts = ",i0)
#endif
  !
  ! Have all MPI processes collect all their SEM averaged values into
  ! a local array.
  !
  if (output_interpolated_inflow_profile) then
    nvar = 24
  else
    nvar = 43
#ifdef ENABLE_OUTPUT_TIME_HISTORY
    nvar = nvar + 6
#endif
  end if
  !
  my_ave_data = init_array(nvar,num_my_sem_pts(mypnum+1),zero) ! F03 AUTO-REALOC
  !
  if (output_interpolated_inflow_profile) then
    call copy_my_sem_fp_profiles(my_prof_data=my_ave_data)
  else
    call copy_my_sem_ave_data(my_ave_data)
  end if
  !
  ! Gather all the SEM averaged values on global root process
  !
  if (ncpu > 1) then
    !
    ! Create the offset and counts for each MPI process
    !
    counts = int( nvar*num_my_sem_pts , kind=int_mpi ) ! F2003 AUTO-REALLOC
    offset = init_array(ncpu,0_int_mpi) ! F2003 AUTO-REALLOC
    !
    do n = 2,ncpu
      offset(n) = offset(n-1) + counts(n-1)
    end do
    !
    ! Allocate the array to collect all the SEM averaged values
    ! from all SEM face points across all MPI processes.
    ! NOTE: Only the global root process will allocate this array to the
    !       actually needed size; all other processes will allocate it such
    !       that there is only a single SEM face point.
    !
    num_fp = merge( num_total_sem_pts , 1 , mypnum == glb_root )
    !
    all_ave_data = init_array(nvar,num_fp,zero) ! F2003 AUTO-REALLOC
    !
    call mpi_gatherv(my_ave_data,size(my_ave_data,kind=int_mpi),mpi_flttyp, &
                     all_ave_data,counts,offset,mpi_flttyp, &
                     glb_root,MPI_COMM_WORLD,mpierr)
    !
  else
    !
    ! For a single MPI process, simply move the
    ! local data to the global data array
    !
    call move_alloc( from=my_ave_data , to=all_ave_data )
    !
  end if
  !
  call collect_output_sem_connectivities(num_sem_pts)
  !
  ! Now have the global root process group all the SEM face points by
  ! their respective SEM boxes so that the averaged data can be output.
  !
  if (mypnum == glb_root) then
    !
    if (output_interpolated_inflow_profile) then
      call write_sem_profile_data(num_sem_pts,all_prof_data=all_ave_data)
    else
      call write_sem_ave_data(num_sem_pts,all_ave_data)
    end if
    !
  end if
  !
  ! Deallocate local arrays before leaving
  !
  if (allocated(all_ave_data)) then
    deallocate ( all_ave_data , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"all_ave_data",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(my_ave_data)) then
    deallocate ( my_ave_data , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"my_ave_data",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(counts)) then
    deallocate ( counts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"counts",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(offset)) then
    deallocate ( offset , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"offset",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(num_my_sem_pts)) then
    deallocate ( num_my_sem_pts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"num_my_sem_pts",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(num_box_sem_pts)) then
    deallocate ( num_box_sem_pts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"num_box_sem_pts",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(num_sem_pts)) then
    deallocate ( num_sem_pts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"num_sem_pts",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (mypnum == glb_root) then
    call output_sem_time_history()
  end if
  !
  ! Since only the global root is writing the SEM averages output file,
  ! throw in an MPI barrier here just in case
  !
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine output_sem_averages
!
!###############################################################################
!
subroutine write_sem_ave_data(num_sem_pts,all_ave_data)
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: num_sem_pts(:,:)
  real(wp), intent(in) :: all_ave_data(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,io,ierr,m,n,n0,nvar,num_fp,sem_idx,cpu_idx
  real(wp) :: y,z,r,u,v,w,e,t,p,ru,rv,rw,rt
  real(wp) :: uu,vv,ww,uv,uw,vw,ut,vt,wt
  real(wp) :: ruu,rvv,rww,ruv,ruw,rvw,rut,rvt,rwt
  real(wp) :: er,eu,ev,ew,et,ep
  real(wp) :: r11,r22,r33,r21,r31,r32
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: offset(:,:)
  real(wp), allocatable :: this_ave(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_sem_ave_data"
  !
  character(len=*), parameter :: fname = "sem_ave_inflow.tec.dat"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nvar = 20
#ifdef ENABLE_OUTPUT_TIME_HISTORY
  nvar = nvar + 6
#endif
  !
  ! Create the offset array for each box on each MPI process
  !
  offset = init_array(ncpu,num_sem_boxes,0) ! F2003 AUTO-REALLOC
  !
  i = 0
  do cpu_idx = 1,ncpu
    do sem_idx = 1,num_sem_boxes
      offset(cpu_idx,sem_idx) = i
      i = i + num_sem_pts(sem_idx,cpu_idx)
    end do
  end do
  !
  ! Open the output file
  !
  open (newunit=io,file=fname,status="replace",action="write", &
                   form="formatted",position="rewind", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Write the Tecplot header
  !
  write (io,1)
  !
  ! Loop over the SEM boxes and write the SEM averaged data for each
  ! box into a separate Tecplot zone
  !
  loop_over_boxes: &
  do sem_idx = 1,num_sem_boxes
    !
    num_fp = sum( num_sem_pts(sem_idx,:) )
    !
    ! Allocate the array to collect the averaged data for the current SEM box
    !
    this_ave = init_array(nvar,num_fp,zero) ! F2003 AUTO-REALLOC
    !
    !############################################################
    !############################################################
    !#####  OUTPUT THE AVERAGED SEM DATA TO A TECPLOT ZONE  #####
    !############################################################
    !############################################################
    !
    i = 0
    !
    loop_over_cpus_ave: &
    do cpu_idx = 1,ncpu
      !
      n0 = offset(cpu_idx,sem_idx)
      !
      loop_over_cpu_pts_ave: &
      do n = 1,num_sem_pts(sem_idx,cpu_idx)
        !
        i = i + 1
        !
        y   = all_ave_data( 1,n0+n)
        z   = all_ave_data( 2,n0+n)
        !
        r   = all_ave_data( 3,n0+n)
        ru  = all_ave_data( 4,n0+n)
        rv  = all_ave_data( 5,n0+n)
        rw  = all_ave_data( 6,n0+n)
        e   = all_ave_data( 7,n0+n)
        !
        t   = all_ave_data( 8,n0+n)
        u   = all_ave_data( 9,n0+n)
        v   = all_ave_data(10,n0+n)
        w   = all_ave_data(11,n0+n)
        p   = all_ave_data(12,n0+n)
        !
        uu  = all_ave_data(13,n0+n)
        vv  = all_ave_data(14,n0+n)
        ww  = all_ave_data(15,n0+n)
        uv  = all_ave_data(16,n0+n)
        uw  = all_ave_data(17,n0+n)
        vw  = all_ave_data(18,n0+n)
        !
        ut  = all_ave_data(19,n0+n)
        vt  = all_ave_data(20,n0+n)
        wt  = all_ave_data(21,n0+n)
        rt  = all_ave_data(22,n0+n)
        !
        ruu = all_ave_data(23,n0+n)
        rvv = all_ave_data(24,n0+n)
        rww = all_ave_data(25,n0+n)
        ruv = all_ave_data(26,n0+n)
        ruw = all_ave_data(27,n0+n)
        rvw = all_ave_data(28,n0+n)
        !
        rut = all_ave_data(29,n0+n)
        rvt = all_ave_data(30,n0+n)
        rwt = all_ave_data(31,n0+n)
        !
        this_ave( 1,i) = z
        this_ave( 2,i) = y
        !
        this_ave( 3,i) = r
        this_ave( 4,i) = u
        this_ave( 5,i) = v
        this_ave( 6,i) = w
        this_ave( 7,i) = t
        this_ave( 8,i) = p
        !
        this_ave( 9,i) = uu - u*u
        this_ave(10,i) = vv - v*v
        this_ave(11,i) = ww - w*w
        this_ave(12,i) = uv - u*v
        this_ave(13,i) = uw - u*w
        this_ave(14,i) = vw - v*w
        !
        this_ave(15,i) = sqrt( max( zero , uu - u*u ) )
        this_ave(16,i) = sqrt( max( zero , vv - v*v ) )
        this_ave(17,i) = sqrt( max( zero , ww - w*w ) )
        !
        this_ave(18,i) = all_ave_data(43,n0+n)
        !
        call yz_to_rt(y,z,r,t, &
                      yc=sem_box(sem_idx)%center%y, &
                      zc=sem_box(sem_idx)%center%z, &
                      return_degrees=true)
        this_ave(19,i) = r
        this_ave(20,i) = t
        !
#ifdef ENABLE_OUTPUT_TIME_HISTORY
        this_ave(21,i) = all_ave_data(44,n0+n)
        this_ave(22,i) = all_ave_data(45,n0+n)
        this_ave(23,i) = all_ave_data(46,n0+n)
        this_ave(24,i) = all_ave_data(47,n0+n)
        this_ave(25,i) = all_ave_data(48,n0+n)
        this_ave(26,i) = all_ave_data(49,n0+n)
#endif
        !
      end do loop_over_cpu_pts_ave
      !
    end do loop_over_cpus_ave
    !
    ! Dimensionalize the averaged data
    !
    call dimensionalize_sem_ave(this_ave)
    !
    ! Write the zone header for this SEM box
    !
    write (io,2) sem_idx,"Averaged",num_fp,size(sem_conn(sem_idx)%ipelem,dim=2)
    !
    ! Write the SEM averaged data for this SEM box
    !
    do m = 1,nvar
      write (io,3) (this_ave(m,i), i=1,num_fp)
    end do
    !
    do n = 1,size(sem_conn(sem_idx)%ipelem,dim=2)
      write (io,4) (sem_conn(sem_idx)%ipelem(i,n), &
                    i=1,size(sem_conn(sem_idx)%ipelem,dim=1))
    end do
    !
    !############################################################
    !############################################################
    !#####  OUTPUT THE EXACT SEM PROFILE TO A TECPLOT ZONE  #####
    !############################################################
    !############################################################
    !
    i = 0
    !
    loop_over_cpus_prof: &
    do cpu_idx = 1,ncpu
      !
      n0 = offset(cpu_idx,sem_idx)
      !
      loop_over_cpu_pts_prof: &
      do n = 1,num_sem_pts(sem_idx,cpu_idx)
        !
        i = i + 1
        !
        y   = all_ave_data( 1,n0+n)
        z   = all_ave_data( 2,n0+n)
        !
        r   = all_ave_data(32,n0+n)
        u   = all_ave_data(33,n0+n)
        v   = all_ave_data(34,n0+n)
        w   = all_ave_data(35,n0+n)
        t   = all_ave_data(36,n0+n)
        p   = r*rgasref_nd*t
        !
        uu  = all_ave_data(37,n0+n)
        vv  = all_ave_data(38,n0+n)
        ww  = all_ave_data(39,n0+n)
        uv  = all_ave_data(40,n0+n)
        uw  = all_ave_data(41,n0+n)
        vw  = all_ave_data(42,n0+n)
        !
        this_ave( 1,i) = z
        this_ave( 2,i) = y
        !
        this_ave( 3,i) = r
        this_ave( 4,i) = u
        this_ave( 5,i) = v
        this_ave( 6,i) = w
        this_ave( 7,i) = t
        this_ave( 8,i) = p
        !
        this_ave( 9,i) = uu
        this_ave(10,i) = vv
        this_ave(11,i) = ww
        this_ave(12,i) = uv
        this_ave(13,i) = uw
        this_ave(14,i) = vw
        !
        this_ave(15,i) = sqrt( max( zero , uu ) )
        this_ave(16,i) = sqrt( max( zero , vv ) )
        this_ave(17,i) = sqrt( max( zero , ww ) )
        !
        this_ave(18,i) = all_ave_data(43,n0+n)
        !
        call yz_to_rt(y,z,r,t, &
                      yc=sem_box(sem_idx)%center%y, &
                      zc=sem_box(sem_idx)%center%z, &
                      return_degrees=true)
        this_ave(19,i) = r
        this_ave(20,i) = t
        !
#ifdef ENABLE_OUTPUT_TIME_HISTORY
        this_ave(21,i) = all_ave_data(44,n0+n)
        this_ave(22,i) = all_ave_data(45,n0+n)
        this_ave(23,i) = all_ave_data(46,n0+n)
        this_ave(24,i) = all_ave_data(47,n0+n)
        this_ave(25,i) = all_ave_data(48,n0+n)
        this_ave(26,i) = all_ave_data(49,n0+n)
#endif
        !
      end do loop_over_cpu_pts_prof
      !
    end do loop_over_cpus_prof
    !
    ! Dimensionalize the averaged data
    !
    call dimensionalize_sem_ave(this_ave)
    !
    ! Write the zone header for this SEM box
    !
    write (io,2) sem_idx,"Profile",num_fp,size(sem_conn(sem_idx)%ipelem,dim=2)
    !
    ! Write the SEM averaged data for this SEM box
    !
    do m = 1,nvar
      write (io,3) (this_ave(m,i), i=1,num_fp)
    end do
    !
    do n = 1,size(sem_conn(sem_idx)%ipelem,dim=2)
      write (io,4) (sem_conn(sem_idx)%ipelem(i,n), &
                    i=1,size(sem_conn(sem_idx)%ipelem,dim=1))
    end do
    !
    !##############################################################
    !##############################################################
    !#####   OUTPUT THE DIFFERENCE BETWEEN THE AVERAGED SEM   #####
    !#####  DATA AND THE EXACT SEM PROFILE TO A TECPLOT ZONE  #####
    !##############################################################
    !##############################################################
    !
    i = 0
    !
    loop_over_cpus_diff: &
    do cpu_idx = 1,ncpu
      !
      n0 = offset(cpu_idx,sem_idx)
      !
      loop_over_cpu_pts_diff: &
      do n = 1,num_sem_pts(sem_idx,cpu_idx)
        !
        i = i + 1
        !
        y   = all_ave_data( 1,n0+n)
        z   = all_ave_data( 2,n0+n)
        !
        r   = all_ave_data( 3,n0+n)
        ru  = all_ave_data( 4,n0+n)
        rv  = all_ave_data( 5,n0+n)
        rw  = all_ave_data( 6,n0+n)
        e   = all_ave_data( 7,n0+n)
        !
        t   = all_ave_data( 8,n0+n)
        u   = all_ave_data( 9,n0+n)
        v   = all_ave_data(10,n0+n)
        w   = all_ave_data(11,n0+n)
        p   = all_ave_data(12,n0+n)
        !
        uu  = all_ave_data(13,n0+n)
        vv  = all_ave_data(14,n0+n)
        ww  = all_ave_data(15,n0+n)
        uv  = all_ave_data(16,n0+n)
        uw  = all_ave_data(17,n0+n)
        vw  = all_ave_data(18,n0+n)
        !
        er  = all_ave_data(32,n0+n)
        eu  = all_ave_data(33,n0+n)
        ev  = all_ave_data(34,n0+n)
        ew  = all_ave_data(35,n0+n)
        et  = all_ave_data(36,n0+n)
        ep  = er*rgasref_nd*et
        !
        r11 = all_ave_data(37,n0+n)
        r22 = all_ave_data(38,n0+n)
        r33 = all_ave_data(39,n0+n)
        r21 = all_ave_data(40,n0+n)
        r31 = all_ave_data(41,n0+n)
        r32 = all_ave_data(42,n0+n)
        !
        this_ave( 1,i) = z
        this_ave( 2,i) = y
        !
        this_ave( 3,i) = r - er
        this_ave( 4,i) = u - eu
        this_ave( 5,i) = v - ev
        this_ave( 6,i) = w - ew
        this_ave( 7,i) = t - et
        this_ave( 8,i) = p - ep
        !
        this_ave( 9,i) = (uu - u*u) - r11
        this_ave(10,i) = (vv - v*v) - r22
        this_ave(11,i) = (ww - w*w) - r33
        this_ave(12,i) = (uv - u*v) - r21
        this_ave(13,i) = (uw - u*w) - r31
        this_ave(14,i) = (vw - v*w) - r32
        !
        this_ave(15,i) = sqrt(max( zero , uu - u*u )) - sqrt(max(zero,r11))
        this_ave(16,i) = sqrt(max( zero , vv - v*v )) - sqrt(max(zero,r22))
        this_ave(17,i) = sqrt(max( zero , ww - w*w )) - sqrt(max(zero,r33))
        !
        this_ave(18,i) = all_ave_data(43,n0+n)
        !
        call yz_to_rt(y,z,r,t, &
                      yc=sem_box(sem_idx)%center%y, &
                      zc=sem_box(sem_idx)%center%z, &
                      return_degrees=true)
        this_ave(19,i) = r
        this_ave(20,i) = t
        !
#ifdef ENABLE_OUTPUT_TIME_HISTORY
        this_ave(21,i) = all_ave_data(44,n0+n)
        this_ave(22,i) = all_ave_data(45,n0+n)
        this_ave(23,i) = all_ave_data(46,n0+n)
        this_ave(24,i) = all_ave_data(47,n0+n)
        this_ave(25,i) = all_ave_data(48,n0+n)
        this_ave(26,i) = all_ave_data(49,n0+n)
#endif
        !
      end do loop_over_cpu_pts_diff
      !
    end do loop_over_cpus_diff
    !
    ! Dimensionalize the averaged data
    !
    call dimensionalize_sem_ave(this_ave)
    !
    ! Write the zone header for this SEM box
    !
    write (io,2) sem_idx,"Averaged-Profile",num_fp, &
                 size(sem_conn(sem_idx)%ipelem,dim=2)
    !
    ! Write the SEM averaged data for this SEM box
    !
    do m = 1,nvar
      write (io,3) (this_ave(m,i), i=1,num_fp)
    end do
    !
    do n = 1,size(sem_conn(sem_idx)%ipelem,dim=2)
      write (io,4) (sem_conn(sem_idx)%ipelem(i,n), &
                    i=1,size(sem_conn(sem_idx)%ipelem,dim=1))
    end do
    !
  end do loop_over_boxes
  !
  ! Close the output file
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Deallocate local arrays before leaving
  !
  if (allocated(this_ave)) then
    deallocate ( this_ave , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"this_ave",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(offset)) then
    deallocate ( offset , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"offset",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ('VARIABLES = "Z"',/, &
            '"Y"',/, &
            '"<greek>r</greek>"',/, &
            '"U"',/, &
            '"V"',/, &
            '"W"',/, &
            '"T"',/, &
            '"P"',/, &
            '"u<math>\\034</math>u<math>\\034</math>"',/, &
            '"v<math>\\034</math>v<math>\\034</math>"',/, &
            '"w<math>\\034</math>w<math>\\034</math>"',/, &
            '"u<math>\\034</math>v<math>\\034</math>"',/, &
            '"u<math>\\034</math>w<math>\\034</math>"',/, &
            '"v<math>\\034</math>w<math>\\034</math>"',/, &
            '"u<math>\\034</math>"',/, &
            '"v<math>\\034</math>"',/, &
            '"w<math>\\034</math>"',/, &
            '"Wall Distance"',/, &
            '"Radius"',/, &
            '"Theta"' &
#ifdef ENABLE_OUTPUT_TIME_HISTORY
         ,/,'"min(u-U)"',/ &
            '"max(u-U)"',/ &
            '"min(v-V)"',/ &
            '"max(v-V)"',/ &
            '"min(w-W)"',/ &
            '"max(w-W)"' &
#endif
            )
  2 format ('ZONE T="SEM BC #',i0,' (',a,')", N=',i0,', E=',i0, &
#ifdef USE_DATAPACKING_POINT
            ', DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL')
#else
            ', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL')
#endif
  3 format (6es14.6)
  4 format (8i8)
  !
end subroutine write_sem_ave_data
!
!###############################################################################
!
subroutine write_sem_profile_data(num_sem_pts,all_prof_data)
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: num_sem_pts(:,:)
  real(wp), intent(in) :: all_prof_data(:,:)
  !
  !.. Local Scalars ..
  integer :: i,io,ierr,m,n,n0,nvar,num_fp,sem_idx,cpu_idx
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: offset(:,:)
  real(wp), allocatable :: this_prof(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_sem_profile_data"
  !
  character(len=*), parameter :: fname = "sem_profile_inflow.tec.dat"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nvar = 24
  !
  ! Create the offset array for each box on each MPI process
  !
  offset = init_array(ncpu,num_sem_boxes,0) ! F2003 AUTO-REALLOC
  !
  i = 0
  do cpu_idx = 1,ncpu
    do sem_idx = 1,num_sem_boxes
      offset(cpu_idx,sem_idx) = i
      i = i + num_sem_pts(sem_idx,cpu_idx)
    end do
  end do
  !
  ! Open the output file
  !
  open (newunit=io,file=fname,status="replace",action="write", &
                   form="formatted",position="rewind", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Write the Tecplot header
  !
  write (io,1)
  !
  ! Loop over the SEM boxes and write the SEM averaged data for each
  ! box into a separate Tecplot zone
  !
  loop_over_boxes: &
  do sem_idx = 1,num_sem_boxes
    !
    num_fp = sum( num_sem_pts(sem_idx,:) )
    !
    ! Allocate the array to collect the averaged data for the current SEM box
    !
    this_prof = init_array(nvar,num_fp,zero) ! F2003 AUTO-REALLOC
    !
    !############################################################
    !############################################################
    !#####  OUTPUT THE AVERAGED SEM DATA TO A TECPLOT ZONE  #####
    !############################################################
    !############################################################
    !
    i = 0
    !
    loop_over_cpus_ave: &
    do cpu_idx = 1,ncpu
      !
      n0 = offset(cpu_idx,sem_idx)
      !
      loop_over_cpu_pts_ave: &
      do n = 1,num_sem_pts(sem_idx,cpu_idx)
        !
        i = i + 1
        !
        do m = 1,nvar
          this_prof(m,i) = all_prof_data(m,n0+n)
        end do
        !
      end do loop_over_cpu_pts_ave
      !
    end do loop_over_cpus_ave
    !
    ! Dimensionalize the averaged data
    !
    call dimensionalize_sem_prof(this_prof)
    !
    ! Write the zone header for this SEM box
    !
    write (io,2) sem_idx,"Profile",num_fp,size(sem_conn(sem_idx)%ipelem,dim=2)
    !
    ! Write the SEM averaged data for this SEM box
    !
    do m = 1,nvar
      write (io,3) (this_prof(m,i), i=1,num_fp)
    end do
    !
    do n = 1,size(sem_conn(sem_idx)%ipelem,dim=2)
      write (io,4) (sem_conn(sem_idx)%ipelem(i,n), &
                    i=1,size(sem_conn(sem_idx)%ipelem,dim=1))
    end do
    !
  end do loop_over_boxes
  !
  ! Close the output file
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Deallocate local arrays before leaving
  !
  if (allocated(this_prof)) then
    deallocate ( this_prof , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"this_prof",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(offset)) then
    deallocate ( offset , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"offset",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ('VARIABLES = "Z"',/,'"Y"',/,'"X"',/, &
            '"<greek>r</greek>"',/,'"U"',/,'"V"',/,'"W"',/,'"T"',/, &
            '"uu"',/,'"vv"',/,'"ww"',/, &
            '"uv"',/,'"uw"',/,'"vw"',/, &
            '"a11"',/,'"a22"',/,'"a33"',/, &
            '"a21"',/,'"a31"',/,'"a32"',/, &
            '"Wall Distance"',/,'"nx"',/,'"ny"',/,'"nz"')
  2 format ('ZONE T="SEM BC #',i0,' (',a,')", N=',i0,', E=',i0, &
            ', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL')
  3 format (6es14.6)
  4 format (8i8)
  !
end subroutine write_sem_profile_data
!
!###############################################################################
!
subroutine collect_output_sem_connectivities(num_sem_pts)
  !
  !.. Use Statements ..
  use geovar, only : nfbnd,bface,face
  !
  !.. Formal Arguments ..
  integer, intent(in) :: num_sem_pts(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,ierr,n,nvar,num_faces,num_total_sem_faces,sem_idx
  integer :: grp_idx,bc_idx
  integer :: n5,n6,nfp5,nfp6
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: num_sem_faces(:,:)
  integer, allocatable :: num_box_sem_faces(:)
  integer, allocatable :: num_my_sem_faces(:)
  integer, allocatable :: my_ipelem(:,:)
  integer, allocatable :: all_ipelem(:,:)
  !
  type(sem_conn_t), allocatable :: my_sem_conn(:)
  !
  integer(int_mpi), allocatable :: counts(:)
  integer(int_mpi), allocatable :: offset(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "collect_output_sem_connectivities"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Have each MPI process count the number of SEM face points it has
  ! in each SEM box and save that to its location within num_sem_faces.
  !
  num_sem_faces = init_array(num_sem_boxes,ncpu,0) ! F2003 AUTO-REALLOC
  !
  do sem_idx = 1,num_sem_boxes
    !
   !grp_idx = sem_box(sem_idx)%grp_idx
    bc_idx = sem_box(sem_idx)%bc_idx
    !
   !n5 = count( bface(5,1:nfbnd) == grp_idx )
   !n6 = count( bface(6,1:nfbnd) == bc_idx )
   !!
   !nfp5 = sum( (face(1:nfbnd)%order+1)**2 , &
   !            mask=(bface(5,1:nfbnd) == grp_idx) )
   !nfp6 = sum( (face(1:nfbnd)%order+1)**2 , &
   !            mask=(bface(6,1:nfbnd) == bc_idx) )
   !!
   !if (n5 /= n6) then
   !  write (iout,900) "sem_idx",sem_idx,"n5",n5,"n6",n6
   !  call stop_gfr(abort,pname,__LINE__,__FILE__, &
   !                "ERROR WITH SEM CONNECTIVITIES")
   !else if (nfp5 /= nfp6) then
   !  write (iout,900) "sem_idx",sem_idx,"nfp5",nfp5,"nfp6",nfp6
   !  call stop_gfr(abort,pname,__LINE__,__FILE__, &
   !                "ERROR WITH SEM CONNECTIVITIES")
   !else if (nfp5 /= my_sem_pts(sem_idx)%num_fp) then
   !  write (iout,900) "sem_idx",sem_idx,"nfp5",nfp5, &
   !                   "my_sem_pts(sem_idx)%num_fp",my_sem_pts(sem_idx)%num_fp
   !  call stop_gfr(abort,pname,__LINE__,__FILE__, &
   !                "ERROR WITH SEM CONNECTIVITIES")
   !end if
   !900 format (/,3(2x,a," = ",i0,/))
    !
    num_faces = sum( face(1:nfbnd)%order**2 , &
                     mask=(bface(6,1:nfbnd) == bc_idx .and. &
                           bface(1,1:nfbnd) == bc_sem) )
    !
    num_sem_faces(sem_idx,mypnum+1) = num_faces
    !
  end do
  !
  ! Now gather the number of SEM face points for all MPI processes
  !
  if (ncpu > 1) then
    call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_DATATYPE_NULL, &
                       num_sem_faces,int(num_sem_boxes,kind=INT_MPI), &
                       mpi_inttyp,MPI_COMM_WORLD,mpierr)
  end if
  !
!#define DEBUG_SEM_CONNECTIVITY
#ifdef DEBUG_SEM_CONNECTIVITY
  if (mypnum == glb_root) then
    write (iout,*)
  end if
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  do i = 0,ncpu-1
    if (mypnum == i) then
      do sem_idx = 1,num_sem_boxes
        do n = 1,ncpu
          write (iout,1) mypnum,sem_idx,n-1,num_sem_faces(sem_idx,n)
        end do
      end do
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end do
  1 format (2x,"CPU #",i2," : num_sem_faces(SEM BOX #",i0,",CPU #",i0,") = ",i0)
#endif
  !
  ! Count the total number of SEM face points
  ! in each box across all MPI processes.
  !
  num_box_sem_faces = sum(num_sem_faces,dim=2) ! F2003 AUTO-REALLOC
  !
#ifdef DEBUG_SEM_CONNECTIVITY
  if (mypnum == glb_root) then
    write (iout,*)
  end if
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  do i = 0,ncpu-1
    if (mypnum == i) then
      do sem_idx = 1,num_sem_boxes
        write (iout,2) mypnum,sem_idx,num_box_sem_faces(sem_idx)
      end do
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end do
  2 format (2x,"CPU #",i2," : num_box_sem_faces(SEM BOX #",i0,") = ",i0)
#endif
  !
  ! Count the total number of SEM face points
  ! across all SEM boxes for each MPI process
  !
  num_my_sem_faces = sum(num_sem_faces,dim=1) ! F2003 AUTO-REALLOC
  !
#ifdef DEBUG_SEM_CONNECTIVITY
  if (mypnum == glb_root) then
    write (iout,*)
  end if
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  do i = 0,ncpu-1
    if (mypnum == i) then
      do n = 1,ncpu
        write (iout,3) mypnum,n-1,num_my_sem_faces(n)
      end do
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end do
  3 format (2x,"CPU #",i2," : num_my_sem_faces(CPU #",i0,") = ",i0)
#endif
  !
  ! Count the total number of SEM face points
  ! across all SEM boxes and all MPI processes
  !
  num_total_sem_faces = sum(num_sem_faces)
  !
#ifdef DEBUG_SEM_CONNECTIVITY
  if (mypnum == glb_root) then
    write (iout,*)
  end if
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  do i = 0,ncpu-1
    if (mypnum == i) then
      write (iout,4) mypnum,num_total_sem_faces
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end do
  4 format (2x,"CPU #",i2," : num_total_sem_faces = ",i0)
#endif
  !
  ! Have all MPI processes collect all their SEM averaged values into
  ! a local array.
  !
  call create_my_sem_connectivity(num_sem_faces,my_sem_conn)
  !
  my_ipelem = init_array(4,num_my_sem_faces(mypnum+1),0) ! F03 AUTO-REALOC
  !
  i = 0
  do sem_idx = 1,num_sem_boxes
    do n = 1,size(my_sem_conn(sem_idx)%ipelem,dim=2)
      i = i + 1
      do j = 1,size(my_sem_conn(sem_idx)%ipelem,dim=1)
        my_ipelem(j,i) = my_sem_conn(sem_idx)%ipelem(j,n)
      end do
    end do
  end do
  !
  ! Gather all the SEM averaged values on global root process
  !
  if (ncpu > 1) then
    !
    ! Create the offset and counts for each MPI process
    !
    counts = int( 4*num_my_sem_faces , kind=int_mpi ) ! F2003 AUTO-REALLOC
    offset = init_array(ncpu,0_int_mpi) ! F2003 AUTO-REALLOC
    !
    do n = 2,ncpu
      offset(n) = offset(n-1) + counts(n-1)
    end do
    !
    ! Allocate the array to collect all the SEM averaged values
    ! from all SEM face points across all MPI processes.
    ! NOTE: Only the global root process will allocate this array to the
    !       actually needed size; all other processes will allocate it such
    !       that there is only a single SEM face point.
    !
    num_faces = merge( num_total_sem_faces , 1 , mypnum == glb_root )
    !
    all_ipelem = init_array(4,num_faces,0) ! F2003 AUTO-REALLOC
    !
    call mpi_gatherv(my_ipelem,size(my_ipelem,kind=int_mpi),mpi_inttyp, &
                     all_ipelem,counts,offset,mpi_inttyp, &
                     glb_root,MPI_COMM_WORLD,mpierr)
    !
  else
    !
    ! For a single MPI process, simply move the
    ! local data to the global data array
    !
    call move_alloc( from=my_ipelem , to=all_ipelem )
    !
  end if
  !
  ! Now have the global root process group all the SEM face points by
  ! their respective SEM boxes so that the averaged data can be output.
  !
  if (mypnum == glb_root) then
    !
    call create_all_sem_connectivity(num_sem_pts,num_sem_faces,all_ipelem)
    !
  end if
  !
  ! Deallocate local arrays before leaving
  !
  if (allocated(all_ipelem)) then
    deallocate ( all_ipelem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"all_ipelem",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(my_ipelem)) then
    deallocate ( my_ipelem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"my_ipelem",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(counts)) then
    deallocate ( counts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"counts",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(offset)) then
    deallocate ( offset , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"offset",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(num_my_sem_faces)) then
    deallocate ( num_my_sem_faces , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"num_my_sem_faces",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(num_box_sem_faces)) then
    deallocate ( num_box_sem_faces , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"num_box_sem_faces",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(num_sem_faces)) then
    deallocate ( num_sem_faces , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"num_sem_faces",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  ! Since only the global root is writing the SEM averages output file,
  ! throw in an MPI barrier here just in case
  !
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine collect_output_sem_connectivities
!
!###############################################################################
!
subroutine create_my_sem_connectivity(num_sem_faces,my_sem_conn)
  !
  !.. Use Statements ..
  use geovar, only : nfbnd,bface,face,n_bnd_flxpts
  !
  !.. Formal Arguments ..
  integer,                        intent(in) :: num_sem_faces(:,:)
  type(sem_conn_t), allocatable, intent(out) :: my_sem_conn(:)
  !
  !.. Local Scalars ..
  integer :: i,k,n,nc,nf,f1,ierr
  integer :: num_fp,num_faces
  integer :: sem_idx,bc_idx,grp_idx,this_geom
  integer :: this_order,min_order,max_order
  !
  !.. Local Derived Types ..
  type :: conn_t
    integer, allocatable :: conn(:,:)
  end type conn_t
  !
  !.. Local Allocatable Arrays ..
  integer,      allocatable :: ip2fp(:)
  integer,      allocatable :: fp2ip(:)
  type(conn_t), allocatable :: tria(:)
  type(conn_t), allocatable :: quad(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_my_sem_connectivity"
  !
continue
  !
  !
  !
  min_order = minval( face(1:nfbnd)%order , dim=1 )
  max_order = maxval( face(1:nfbnd)%order , dim=1 )
  !
  allocate ( tria(min_order:max_order) )
  allocate ( quad(min_order:max_order) )
  !
  do this_order = min_order,max_order
    call create_tria_connectivity(this_order,tria(this_order)%conn)
    call create_quad_connectivity(this_order,quad(this_order)%conn)
  end do
  !
  !
  !
  if (allocated(my_sem_conn)) then
    deallocate ( my_sem_conn , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"my_sem_conn",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( my_sem_conn(1:num_sem_boxes) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"my_sem_conn",1,__LINE__,__FILE__,ierr,error_message)
  !
  !
  do sem_idx = 1,num_sem_boxes
    !
    num_fp = my_sem_pts(sem_idx)%num_fp
    num_faces = num_sem_faces(sem_idx,mypnum+1)
    !
    my_sem_conn(sem_idx)%ipelem = init_array(4,num_faces,0) ! F03 AUTO-REALLOC
    !
   !grp_idx = sem_box(sem_idx)%grp_idx
    bc_idx = sem_box(sem_idx)%bc_idx
    !
    i = 0
    !
    do nf = 1,nfbnd
      !
      if (bface(1,nf) /= bc_sem) cycle
      !
     !if (bface(5,nf) /= grp_idx) cycle
      if (bface(6,nf) /= bc_idx) cycle
      !
      f1 = face(nf)%beg_fp
      !
      this_geom = face(nf)%geom
      this_order = face(nf)%order
      nc = this_order*this_order
      !
      if (this_geom == Geom_Tria) then
        do n = 1,nc
          i = i + 1
          my_sem_conn(sem_idx)%ipelem(:,i) = f1 + tria(this_order)%conn(:,n)
        end do
      else
        do n = 1,nc
          i = i + 1
          my_sem_conn(sem_idx)%ipelem(:,i) = f1 + quad(this_order)%conn(:,n)
        end do
      end if
      !
    end do
    !
    ip2fp = init_array(num_fp,0)       ! F2003 AUTO-REALLOC
    fp2ip = init_array(n_bnd_flxpts,0) ! F2003 AUTO-REALLOC
    !
    do n = 1,num_fp
      k = my_sem_pts(sem_idx)%sem_fp(n)%idx
      ip2fp(n) = k  ! 1 <= n <= num_fp
      fp2ip(k) = n  ! 1 <= k <= n_bnd_flxpts
    end do
    !
#ifdef DEBUG_SEM_CONNECTIVITY
    k = 100*(sem_idx+2) + mypnum
    do n = 1,num_fp
      write (k,2) "ip2fp",n,ip2fp(n)
    end do
    !
    k = 100*(sem_idx+4) + mypnum
    do n = 1,n_bnd_flxpts
      write (k,2) "fp2ip",n,fp2ip(n)
    end do
    !
    do nf = 1,num_faces
      write (sem_idx*100+mypnum,1) nf,(my_sem_conn(sem_idx)%ipelem(n,nf), &
                                    n=1,size(my_sem_conn(sem_idx)%ipelem,dim=1))
    end do
    !
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    flush (100*(sem_idx+0)+mypnum)
    flush (100*(sem_idx+2)+mypnum)
    flush (100*(sem_idx+4)+mypnum)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    !
    1 format (2x,i6,4x,4i6)
    2 format (a,"(",i5,") = ",i0)
#endif
    !
    do nf = 1,num_faces
      do n = 1,size(my_sem_conn(sem_idx)%ipelem,dim=1)
        k = my_sem_conn(sem_idx)%ipelem(n,nf)
        my_sem_conn(sem_idx)%ipelem(n,nf) = fp2ip(k)
      end do
    end do
    !
  end do
  !
end subroutine create_my_sem_connectivity
!
!###############################################################################
!
subroutine create_all_sem_connectivity(num_sem_pts,num_sem_faces,all_ipelem)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: num_sem_pts(:,:)
  integer, intent(in) :: num_sem_faces(:,:)
  integer, intent(in) :: all_ipelem(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,n,nf,np,n0,f0,ierr
  integer :: num_faces,sem_idx,cpu_idx
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: np_offset(:,:)
  integer, allocatable :: nf_offset(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_my_sem_connectivity"
  !
continue
  !
  if (allocated(sem_conn)) then
    deallocate ( sem_conn , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"sem_conn",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( sem_conn(1:num_sem_boxes) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sem_conn",1,__LINE__,__FILE__,ierr,error_message)
  !
  !
  !
  np_offset = init_array(ncpu,num_sem_boxes,0) ! F2003 AUTO-REALLOC
  !
  do sem_idx = 1,num_sem_boxes
    np = 0
    do cpu_idx = 1,ncpu
      np_offset(cpu_idx,sem_idx) = np
      np = np + num_sem_pts(sem_idx,cpu_idx)
    end do
  end do
  !
  !
  !
  nf_offset = init_array(ncpu,num_sem_boxes,0) ! F2003 AUTO-REALLOC
  !
  nf = 0
  do cpu_idx = 1,ncpu
    do sem_idx = 1,num_sem_boxes
      nf_offset(cpu_idx,sem_idx) = nf
      nf = nf + num_sem_faces(sem_idx,cpu_idx)
    end do
  end do
  !
  !
  !
  do sem_idx = 1,num_sem_boxes
    !
    num_faces = sum( num_sem_faces(sem_idx,:) )
    sem_conn(sem_idx)%ipelem = init_array(4,num_faces,0) ! F2003 AUTO-REALLOC
    !
    i = 0
    !
    do cpu_idx = 1,ncpu
      !
      n0 = nf_offset(cpu_idx,sem_idx)
      f0 = np_offset(cpu_idx,sem_idx)
      !
      do nf = 1,num_sem_faces(sem_idx,cpu_idx)
        !
        i = i + 1
        !
        do j = 1,size(all_ipelem,dim=1)
          sem_conn(sem_idx)%ipelem(j,i) = f0 + all_ipelem(j,n0+nf)
        end do
        !
      end do
      !
    end do
    !
    if (i /= num_faces) then
      write (iout,1) sem_idx,i,num_faces
      flush (iout)
    end if
    !
  end do
  !
  ! Format Statements
  !
  1 format (/,2x,"Error copying all_ipelem into ", &
                 "sem_conn(",i0,")%ipelem while trying",/, &
              2x,"to create connectivity for the SEM ", &
                 "time-averaged boundary states!",/, &
              9x,"number of SEM faces found = ",i0,/, &
              6x,"number of SEM faces expected = ",i0,/)

  !
end subroutine create_all_sem_connectivity
!
!###############################################################################
!
subroutine create_tria_connectivity(no1d,conn)
  !
  !.. Formal Arguments ..
  integer,               intent(in) :: no1d
  integer, allocatable, intent(out) :: conn(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,n,no2d
  integer :: p1,p2,p3,p4
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: counter(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_tria_connectivity"
  !
continue
  !
  no2d = no1d*no1d
  !
  conn = init_array(4,no2d,0) ! F2003 AUTO-REALLOC
  !
  counter = init_array(no1d+1,no1d+1,0) ! F2003 AUTO-REALLOC
  !
  n = 0
  do j = 1,no1d+1
    do i = 1,no1d+2-j
      n = n + 1
      counter(j,i) = n
    end do
  end do
  !
  n = 0
  do j = 1,no1d+1
    do i = 1,no1d+1-j
      !
      p1 = counter(j,i)
      p2 = counter(j,i+1)
      p3 = counter(j+1,i)
      p4 = counter(j+1,i+1)
      !
      n = n + 1
      conn(1,n) = p1
      conn(2,n) = p2
      conn(3,n) = p3
      conn(4,n) = p3
      !
      if (p4 /= 0) then
        n = n + 1
        conn(1,n) = p2
        conn(2,n) = p4
        conn(3,n) = p3
        conn(4,n) = p3
      end if
      !
    end do
  end do
  !
  conn(:,:) = conn(:,:) - 1
  !
end subroutine create_tria_connectivity
!
!###############################################################################
!
subroutine create_quad_connectivity(no1d,conn)
  !
  !.. Formal Arguments ..
  integer,               intent(in) :: no1d
  integer, allocatable, intent(out) :: conn(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,n,no2d
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_quad_connectivity"
  !
continue
  !
  no2d = no1d*no1d
  !
  conn = init_array(4,no2d,0) ! F2003 AUTO-REALLOC
  !
  n = 0
  do j = 1,no1d
    do i = 1,no1d
      !
      n = n + 1
      conn(1,n) = (no1d+1)*(j-1) + i
      conn(2,n) = (no1d+1)*(j-1) + i + 1
      conn(3,n) = (no1d+1)*(j  ) + i + 1
      conn(4,n) = (no1d+1)*(j  ) + i
      !
    end do
  end do
  !
  conn(:,:) = conn(:,:) - 1
  !
end subroutine create_quad_connectivity
!
!###############################################################################
!
subroutine output_sem_time_history()
  !
  !.. Use Statements ..
  use ovar, only : num_timesteps,rhoref,aref,tref
  !
  !.. Local Scalars ..
  integer  :: i,io,ierr,n,sem_idx
  real(wp) :: r,u,v,w,e,t
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_sem_time_history"
  character(len=*), parameter :: fname = "sem_fp_time_histories.tec.dat"
  !
continue
#ifdef ENABLE_OUTPUT_TIME_HISTORY
  !
  open (newunit=io,file=fname,status="replace",action="write", &
                   form="formatted",position="rewind", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (io,1)
  !
  do sem_idx = 1,num_sem_boxes
    !
    do n = 1,my_sem_pts(sem_idx)%num_fp
      !
      write (io,2) sem_idx,n, &
                   my_sem_pts(sem_idx)%sem_fp(n)%y, &
                   my_sem_pts(sem_idx)%sem_fp(n)%z, &
                   num_timesteps
      !
      do i = lbound(my_sem_pts(sem_idx)%sem_fp(n)%time_history,dim=2), &
             ubound(my_sem_pts(sem_idx)%sem_fp(n)%time_history,dim=2)
        !
        r = my_sem_pts(sem_idx)%sem_fp(n)%time_history(nec,i)*rhoref
        u = my_sem_pts(sem_idx)%sem_fp(n)%time_history(nmx,i)*aref
        v = my_sem_pts(sem_idx)%sem_fp(n)%time_history(nmy,i)*aref
        w = my_sem_pts(sem_idx)%sem_fp(n)%time_history(nmz,i)*aref
        e = my_sem_pts(sem_idx)%sem_fp(n)%time_history(nee,i)*rhoref*aref*aref
        t = my_sem_pts(sem_idx)%sem_fp(n)%time_history(nq+1,i)*tref
        !
        write (io,3) i,r,u,v,w,e,t
        !
      end do
      !
    end do
    !
  end do
  !
  ! Format Statements
  !
  1 format ('VARIABLES = "Timestep" "rho" "u" "v" "w" "e" "T"')
  2 format ('ZONE T="Box-',i0,',FP-',i0,': y=',f9.6,',z=',f9.6, &
            '", I=',i0,", DATAPACKING=POINT")
  3 format (i8,6es14.6)
  !
#endif
end subroutine output_sem_time_history
!
!###############################################################################
!
pure &
subroutine copy_my_sem_ave_data(my_ave_data)
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: my_ave_data(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,n,sem_idx
  real(wp) :: umin,umax,vmin,vmax,wmin,wmax
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "copy_my_sem_ave_data"
  !
continue
  !
  i = 0
  do sem_idx = 1,num_sem_boxes
    do n = 1,my_sem_pts(sem_idx)%num_fp
      !
      i = i + 1
      !
      my_ave_data( 1,i) = my_sem_pts(sem_idx)%sem_fp(n)%y
      my_ave_data( 2,i) = my_sem_pts(sem_idx)%sem_fp(n)%z
      !
      my_ave_data( 3,i) = sem_avg(sem_idx)%fp_avg(n)%r
      my_ave_data( 4,i) = sem_avg(sem_idx)%fp_avg(n)%ru
      my_ave_data( 5,i) = sem_avg(sem_idx)%fp_avg(n)%rv
      my_ave_data( 6,i) = sem_avg(sem_idx)%fp_avg(n)%rw
      my_ave_data( 7,i) = sem_avg(sem_idx)%fp_avg(n)%e
      !
      my_ave_data( 8,i) = sem_avg(sem_idx)%fp_avg(n)%t
      my_ave_data( 9,i) = sem_avg(sem_idx)%fp_avg(n)%u
      my_ave_data(10,i) = sem_avg(sem_idx)%fp_avg(n)%v
      my_ave_data(11,i) = sem_avg(sem_idx)%fp_avg(n)%w
      my_ave_data(12,i) = sem_avg(sem_idx)%fp_avg(n)%p
      !
      my_ave_data(13,i) = sem_avg(sem_idx)%fp_avg(n)%uu
      my_ave_data(14,i) = sem_avg(sem_idx)%fp_avg(n)%vv
      my_ave_data(15,i) = sem_avg(sem_idx)%fp_avg(n)%ww
      my_ave_data(16,i) = sem_avg(sem_idx)%fp_avg(n)%uv
      my_ave_data(17,i) = sem_avg(sem_idx)%fp_avg(n)%uw
      my_ave_data(18,i) = sem_avg(sem_idx)%fp_avg(n)%vw
      !
      my_ave_data(19,i) = sem_avg(sem_idx)%fp_avg(n)%ut
      my_ave_data(20,i) = sem_avg(sem_idx)%fp_avg(n)%vt
      my_ave_data(21,i) = sem_avg(sem_idx)%fp_avg(n)%wt
      my_ave_data(22,i) = sem_avg(sem_idx)%fp_avg(n)%rt
      !
      my_ave_data(23,i) = sem_avg(sem_idx)%fp_avg(n)%ruu
      my_ave_data(24,i) = sem_avg(sem_idx)%fp_avg(n)%rvv
      my_ave_data(25,i) = sem_avg(sem_idx)%fp_avg(n)%rww
      my_ave_data(26,i) = sem_avg(sem_idx)%fp_avg(n)%ruv
      my_ave_data(27,i) = sem_avg(sem_idx)%fp_avg(n)%ruw
      my_ave_data(28,i) = sem_avg(sem_idx)%fp_avg(n)%rvw
      !
      my_ave_data(29,i) = sem_avg(sem_idx)%fp_avg(n)%rut
      my_ave_data(30,i) = sem_avg(sem_idx)%fp_avg(n)%rvt
      my_ave_data(31,i) = sem_avg(sem_idx)%fp_avg(n)%rwt
      !
      my_ave_data(32,i) = my_sem_pts(sem_idx)%sem_fp(n)%ravg
      my_ave_data(33,i) = my_sem_pts(sem_idx)%sem_fp(n)%uavg
      my_ave_data(34,i) = my_sem_pts(sem_idx)%sem_fp(n)%vavg
      my_ave_data(35,i) = my_sem_pts(sem_idx)%sem_fp(n)%wavg
      my_ave_data(36,i) = my_sem_pts(sem_idx)%sem_fp(n)%tavg
      !
      my_ave_data(37,i) = my_sem_pts(sem_idx)%sem_fp(n)%r11
      my_ave_data(38,i) = my_sem_pts(sem_idx)%sem_fp(n)%r22
      my_ave_data(39,i) = my_sem_pts(sem_idx)%sem_fp(n)%r33
      my_ave_data(40,i) = my_sem_pts(sem_idx)%sem_fp(n)%r21
      my_ave_data(41,i) = my_sem_pts(sem_idx)%sem_fp(n)%r31
      my_ave_data(42,i) = my_sem_pts(sem_idx)%sem_fp(n)%r32
      !
      my_ave_data(43,i) = my_sem_pts(sem_idx)%sem_fp(n)%wall%distance
      !
#ifdef ENABLE_OUTPUT_TIME_HISTORY
      !
      umin = minval( my_sem_pts(sem_idx)%sem_fp(n)%time_history(nmx,:) , dim=1 )
      vmin = minval( my_sem_pts(sem_idx)%sem_fp(n)%time_history(nmy,:) , dim=1 )
      wmin = minval( my_sem_pts(sem_idx)%sem_fp(n)%time_history(nmz,:) , dim=1 )
      !
      umax = maxval( my_sem_pts(sem_idx)%sem_fp(n)%time_history(nmx,:) , dim=1 )
      vmax = maxval( my_sem_pts(sem_idx)%sem_fp(n)%time_history(nmy,:) , dim=1 )
      wmax = maxval( my_sem_pts(sem_idx)%sem_fp(n)%time_history(nmz,:) , dim=1 )
      !
      my_ave_data(44,i) = umin - my_ave_data(33,i)
      my_ave_data(45,i) = umax - my_ave_data(33,i)
      my_ave_data(46,i) = vmin - my_ave_data(34,i)
      my_ave_data(47,i) = vmax - my_ave_data(34,i)
      my_ave_data(48,i) = wmin - my_ave_data(35,i)
      my_ave_data(49,i) = wmax - my_ave_data(35,i)
      !
#endif
      !
    end do
  end do
  !
end subroutine copy_my_sem_ave_data
!
!###############################################################################
!
pure &
subroutine copy_my_sem_fp_profiles(my_prof_data)
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: my_prof_data(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,sem_idx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "copy_my_sem_fp_profiles"
  !
continue
  !
  i = 0
  do sem_idx = 1,num_sem_boxes
    do n = 1,my_sem_pts(sem_idx)%num_fp
      !
      i = i + 1
      !
      my_prof_data( 1,i) = my_sem_pts(sem_idx)%sem_fp(n)%z
      my_prof_data( 2,i) = my_sem_pts(sem_idx)%sem_fp(n)%y
      my_prof_data( 3,i) = my_sem_pts(sem_idx)%sem_fp(n)%x
      !
      my_prof_data( 4,i) = my_sem_pts(sem_idx)%sem_fp(n)%ravg
      my_prof_data( 5,i) = my_sem_pts(sem_idx)%sem_fp(n)%uavg
      my_prof_data( 6,i) = my_sem_pts(sem_idx)%sem_fp(n)%vavg
      my_prof_data( 7,i) = my_sem_pts(sem_idx)%sem_fp(n)%wavg
      my_prof_data( 8,i) = my_sem_pts(sem_idx)%sem_fp(n)%tavg
      !
      my_prof_data( 9,i) = my_sem_pts(sem_idx)%sem_fp(n)%r11
      my_prof_data(10,i) = my_sem_pts(sem_idx)%sem_fp(n)%r22
      my_prof_data(11,i) = my_sem_pts(sem_idx)%sem_fp(n)%r33
      my_prof_data(12,i) = my_sem_pts(sem_idx)%sem_fp(n)%r21
      my_prof_data(13,i) = my_sem_pts(sem_idx)%sem_fp(n)%r31
      my_prof_data(14,i) = my_sem_pts(sem_idx)%sem_fp(n)%r32
      !
      my_prof_data(15,i) = my_sem_pts(sem_idx)%sem_fp(n)%a11
      my_prof_data(16,i) = my_sem_pts(sem_idx)%sem_fp(n)%a22
      my_prof_data(17,i) = my_sem_pts(sem_idx)%sem_fp(n)%a33
      my_prof_data(18,i) = my_sem_pts(sem_idx)%sem_fp(n)%a21
      my_prof_data(19,i) = my_sem_pts(sem_idx)%sem_fp(n)%a31
      my_prof_data(20,i) = my_sem_pts(sem_idx)%sem_fp(n)%a32
      !
      my_prof_data(21,i) = my_sem_pts(sem_idx)%sem_fp(n)%wall%distance
      my_prof_data(22,i) = my_sem_pts(sem_idx)%sem_fp(n)%wall%xnrm
      my_prof_data(23,i) = my_sem_pts(sem_idx)%sem_fp(n)%wall%ynrm
      my_prof_data(24,i) = my_sem_pts(sem_idx)%sem_fp(n)%wall%znrm
      !
    end do
  end do
  !
end subroutine copy_my_sem_fp_profiles
!
!###############################################################################
!
pure &
subroutine dimensionalize_sem_ave(this_ave)
  !
  !.. Use Statements ..
  use ovar, only : rhoref,aref,tref
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: this_ave(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,npts
  real(wp) :: aref2,rhoa2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "dimensionalize_sem_ave"
  !
continue
  !
  npts = size(this_ave,dim=2)
  !
  aref2 = aref * aref
  rhoa2 = rhoref * aref2
  !
  do n = 1,npts
    !
    this_ave( 3,n) = this_ave( 3,n) * rhoref
    this_ave( 4,n) = this_ave( 4,n) * aref
    this_ave( 5,n) = this_ave( 5,n) * aref
    this_ave( 6,n) = this_ave( 6,n) * aref
    this_ave( 7,n) = this_ave( 7,n) * tref
    this_ave( 8,n) = this_ave( 8,n) * rhoa2
    !
    this_ave( 9,n) = this_ave( 9,n) * aref2
    this_ave(10,n) = this_ave(10,n) * aref2
    this_ave(11,n) = this_ave(11,n) * aref2
    this_ave(12,n) = this_ave(12,n) * aref2
    this_ave(13,n) = this_ave(13,n) * aref2
    this_ave(14,n) = this_ave(14,n) * aref2
    !
    this_ave(15,n) = this_ave(15,n) * aref
    this_ave(16,n) = this_ave(16,n) * aref
    this_ave(17,n) = this_ave(17,n) * aref
    !
#ifdef ENABLE_OUTPUT_TIME_HISTORY
    this_ave(21,n) = this_ave(21,n) * aref
    this_ave(22,n) = this_ave(22,n) * aref
    this_ave(23,n) = this_ave(23,n) * aref
    this_ave(24,n) = this_ave(24,n) * aref
    this_ave(25,n) = this_ave(25,n) * aref
    this_ave(26,n) = this_ave(26,n) * aref
#endif
    !
  end do
  !
end subroutine dimensionalize_sem_ave
!
!###############################################################################
!
pure &
subroutine dimensionalize_sem_prof(this_prof)
  !
  !.. Use Statements ..
  use ovar, only : rhoref,aref,tref,grid_scaling_factor
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: this_prof(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,npts
  real(wp) :: aref2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "dimensionalize_sem_prof"
  !
continue
  !
  npts = size(this_prof,dim=2)
  !
  aref2 = aref * aref
  !
  do n = 1,npts
    !
    this_prof( 1,n) = this_prof( 1,n) / grid_scaling_factor
    this_prof( 2,n) = this_prof( 2,n) / grid_scaling_factor
    this_prof( 3,n) = this_prof( 3,n) / grid_scaling_factor
    !
    this_prof( 4,n) = this_prof( 4,n) * rhoref
    this_prof( 5,n) = this_prof( 5,n) * aref
    this_prof( 6,n) = this_prof( 6,n) * aref
    this_prof( 7,n) = this_prof( 7,n) * aref
    this_prof( 8,n) = this_prof( 8,n) * tref
    !
    this_prof( 9,n) = this_prof( 9,n) * aref2
    this_prof(10,n) = this_prof(10,n) * aref2
    this_prof(11,n) = this_prof(11,n) * aref2
    this_prof(12,n) = this_prof(12,n) * aref2
    this_prof(13,n) = this_prof(13,n) * aref2
    this_prof(14,n) = this_prof(14,n) * aref2
    !
    this_prof(15,n) = this_prof(15,n) * aref
    this_prof(16,n) = this_prof(16,n) * aref
    this_prof(17,n) = this_prof(17,n) * aref
    this_prof(18,n) = this_prof(18,n) * aref
    this_prof(19,n) = this_prof(19,n) * aref
    this_prof(20,n) = this_prof(20,n) * aref
    !
    this_prof(21,n) = this_prof(21,n) / grid_scaling_factor
    !
  end do
  !
end subroutine dimensionalize_sem_prof
!
!###############################################################################
!
subroutine dump_sem_data(cloc)
  !
  !.. Use Statements ..
  use random_mod, only : iseed
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: cloc
  !
  !.. Local Scalars ..
  integer  :: i,io,ierr,sem_idx,eddy_idx
  !
  character(len=50) :: fname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "dump_sem_data"
  !
continue
!!#define ENABLE_DUMP_SEM_DATA
#ifdef ENABLE_DUMP_SEM_DATA
  !
  write (fname,1) cloc,mypnum
  !
  !
  open (newunit=io,file=fname,status="replace",action="write", &
                   form="formatted",position="rewind", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (io,4) (iseed(i), i=1,size(iseed))
  write (io,5)
  !
  do sem_idx = 1,num_sem_boxes
    !
    write (io,2) sem_idx, &
                 sem_box(sem_idx)%num_eddies, &
                 sem_box(sem_idx)%bc_idx, &
                 sem_box(sem_idx)%grp_idx, &
                 sem_box(sem_idx)%length_scale, &
                 sem_box(sem_idx)%convective_velocity, &
                 sem_box(sem_idx)%boxmin(1),sem_box(sem_idx)%boxmax(1), &
                 sem_box(sem_idx)%boxmin(2),sem_box(sem_idx)%boxmax(2), &
                 sem_box(sem_idx)%boxmin(3),sem_box(sem_idx)%boxmax(3), &
                 sem_box(sem_idx)%boxmin(4),sem_box(sem_idx)%boxmax(4), &
                 sem_box(sem_idx)%boxmin(5),sem_box(sem_idx)%boxmax(5), &
                 sem_box(sem_idx)%box_volume, &
                 sem_box(sem_idx)%is_cylindrical
    !
    do eddy_idx = 1,sem_box(sem_idx)%num_eddies
      !
      write (io,3) eddy_idx, &
                   sem_box(sem_idx)%eddy(eddy_idx)%x, &
                   sem_box(sem_idx)%eddy(eddy_idx)%y, &
                   sem_box(sem_idx)%eddy(eddy_idx)%z, &
                   sem_box(sem_idx)%eddy(eddy_idx)%r, &
                   sem_box(sem_idx)%eddy(eddy_idx)%t, &
                   sem_box(sem_idx)%eddy(eddy_idx)%x_intensity, &
                   sem_box(sem_idx)%eddy(eddy_idx)%y_intensity, &
                   sem_box(sem_idx)%eddy(eddy_idx)%z_intensity, &
                   sem_box(sem_idx)%eddy(eddy_idx)%x_frequency, &
                   sem_box(sem_idx)%eddy(eddy_idx)%y_frequency, &
                   sem_box(sem_idx)%eddy(eddy_idx)%z_frequency, &
                   sem_box(sem_idx)%eddy(eddy_idx)%volume
      !
    end do
    !
  end do
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Format Statements
  !
  1 format ("sem_dump.",a,".cpu-",i0,".dat")
  2 format (/, 2x,"SEM Box # ",i0,/, &
              13x,"num_eddies = ",i0,/, &
              17x,"bc_idx = ",i0,/, &
              16x,"grp_idx = ",i0,/, &
              11x,"length_scale = ",es14.6,/, &
               4x,"convective_velocity = ",es14.6,/, &
              13x,"x[min,max] = [",es14.6,", ",es14.6,"]",/, &
              13x,"y[min,max] = [",es14.6,", ",es14.6,"]",/, &
              13x,"z[min,max] = [",es14.6,", ",es14.6,"]",/, &
              13x,"r[min,max] = [",es14.6,", ",es14.6,"]",/, &
              13x,"t[min,max] = [",es14.6,", ",es14.6,"]",/, &
              13x,"box volume = ",es14.6,/, &
              10x,"is cylindrical = ",l1)
  3 format (/, 4x,"Eddy # ",i0,/, &
              17x,"x = ",es14.6,/, &
              17x,"y = ",es14.6,/, &
              17x,"z = ",es14.6,/, &
              17x,"r = ",es14.6,/, &
              17x,"t = ",es14.6,/, &
               7x,"x_intensity = ",es14.6,/, &
               7x,"y_intensity = ",es14.6,/, &
               7x,"z_intensity = ",es14.6,/, &
               7x,"x_frequency = ",es14.6,/, &
               7x,"y_frequency = ",es14.6,/, &
               7x,"z_frequency = ",es14.6,/, &
              12x,"volume = ",es14.6)
  4 format (/,"Random seed = ",100(i0,:,", "))
  5 format (//)
  !
#endif
end subroutine dump_sem_data
!
!###############################################################################
!
subroutine output_bface
  !
  !.. Use Statements ..
  use geovar, only : bface,face
  !
  !.. Local Scalars ..
  integer :: i,ii,n,nf,maxint,myio,ierr
  character(len=100) :: write_format
  character(len=100) :: file_format
  character(len=100) :: my_file
  !
  !.. Local Arrays ..
  integer :: intlen(1:size(bface,dim=1)+1)
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: tmp_bface(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_bface"
  !
continue
  !
  tmp_bface = bface ! F2003 AUTO-REALLOC
  !
 !do nf = 1,size(bface,dim=2)
 !  !
 !  ii = huge(0)
 !  do i = 1,size(bface,dim=2)
 !    if (bface(1,i) < 0) cycle
 !    if (bface(2,i) < ii) then
 !      ii = bface(2,i)
 !      n = i
 !    end if
 !  end do
 ! !n = minloc(bface(2,:),dim=1,mask=bface(1,:)<0)
 !  !
 !  tmp_bface(:,nf) = bface(:,n)
 !  bface(1,n) = -abs(bface(1,n))
 !  write (iout,'("nf,n = ",i0,", ",i0,", ",i0)') nf,n,bface(2,n)
 !  !
 !end do
 !!
 !bface(1,:) = abs(bface(1,:))
  !
  ! Get the maximum length of the string containing the name of the BC types
  !
  intlen(1) = maxval( len_trim( bc_integer_to_string(tmp_bface(1,:)) ) )
  !
  do i = 1,size(intlen)-1
    n = 1
    maxint = maxval(tmp_bface(i,:))
    check_bface_size: do
      if (maxint/(10**n) == 0) then
        intlen(i+1) = n+1
        exit check_bface_size
      else
        n = n + 1
        cycle check_bface_size
      end if
    end do check_bface_size
  end do
  !
  ! Primary write of the format statement which has an unlimited
  ! '*' repetition multiplier on the last quantity to allow for
  ! an arbitrarily large number of values contained in intlen.
  ! The ":" within the part that is unlimitedly repeated serves
  ! to terminate the format statement at that point preventing
  ! it from going back and adding on an extra ',3x,i' at the end.
  !
  write (write_format,1) (intlen(i),i=1,size(intlen))
 !write (iout,'(a,a,a)') "write_format = '",trim(write_format),"'"
  !
  ! Finally, add on the parentheses to the format statement
  !
  write (write_format,2) trim(write_format)
 !write (iout,'(a,a,a)') "write_format = '",trim(write_format),"'"
  !
  ! Get the number of digits of the integer for the number of cpus
  !
  n = 1
  check_ncpu_size: do
    if (ncpu/(10**n) == 0) then
      write (file_format,3) n
      exit check_ncpu_size
    else
      n = n + 1
      cycle check_ncpu_size
    end if
  end do check_ncpu_size
  !
  ! Create the name of the output file
  !
  write (my_file,fmt=trim(file_format)) mypnum
  !
  ! Open the output file
  !
  open (newunit=myio,file=trim(my_file),status="replace",action="write", &
                     form="formatted",position="rewind", &
                     iostat=ierr,iomsg=error_message)
  call io_error(pname,trim(my_file),1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Output the bface array
  !
  do n = 1,size(tmp_bface,dim=2)
    write (myio,fmt=trim(write_format)) &
                               trim(bc_integer_to_string(tmp_bface(1,n))), &
                                   (tmp_bface(i,n),i=1,size(tmp_bface,dim=1))
  end do
  !
  ! Close the output file
  !
  close (myio,iostat=ierr,iomsg=error_message)
  call io_error(pname,trim(my_file),2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Format Statements
  !
  1 format ('a',i0,',i',i0,',": ",1x,i',i0,',1x,i',i0,',1x,i',i0, &
            ',5x,1x,i',i0,',1x,i',i0,',1x,i',i0,',1x,i',i0,',1x,i',i0, &
            100(',3x,i',i0,:))
           !*(',3x,i',i0,:))
  2 format ('(',a,')')
  3 format ('("sem_bface.",i0.',i0,',".dat")')
  !
end subroutine output_bface
!
!###############################################################################
!
pure &
subroutine yz_to_rt(y,z,r,t,yc,zc,return_degrees)
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: y
  real(wp),    intent(in) :: z
  real(wp), intent(inout) :: r
  real(wp), intent(inout) :: t
  !
  !.. Optional Arguments ..
  real(wp),    optional, intent(in) :: yc
  real(wp),    optional, intent(in) :: zc
  logical(lk), optional, intent(in) :: return_degrees
  !
  !.. Local Scalars ..
  real(wp)    :: y0,z0,dy,dz
  logical(lk) :: returns_degrees
  !
continue
  !
  y0 = zero
  z0 = zero
  returns_degrees = fals
  !
  if (present(yc)) y0 = yc
  if (present(zc)) z0 = zc
  if (present(return_degrees)) returns_degrees = return_degrees
  !
  dy = y - y0
  dz = z - z0
  !
  r = hypot(dy,dz)
  !
  if (r /= zero) then
    t = atan2(dz,dy)
    if (t < zero) t = t + two*pi
    if (returns_degrees) t = t*radians2degrees
  else
    t = zero
  end if
  !
end subroutine yz_to_rt
!
!###############################################################################
!
#include "Functions/check_for_bad_number.f90"
!
!###############################################################################
!
end module synthetic_eddy_mod
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
module digital_filtering_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nee
  !
  implicit none
  !
  private
  !
  public :: df
 !public :: dfcg
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !###                                                                       ###
  !###   MODULE DERIVED-TYPES, PARAMETERS, VARIABLES, AND ALLOCATABLE        ###
  !###        ARRAYS FOR USING THE DIGITAL FILTERING METHOD (DFM)            ###
  !###                                                                       ###
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! DIGITAL FILTERING VARIABLES AND ARRAYS
  !
  integer, save :: nfilx
  integer, save :: nfily
  integer, save :: nfilz
  integer, save :: npydf
  integer, save :: npzdf
  integer, save :: mxnfily
  integer, save :: mxnfilz
  !
  real(wp), save :: ymin
  real(wp), save :: ymax
  real(wp), save :: zmin
  real(wp), save :: zmax
  real(wp), save :: dx
  real(wp), save :: dy
  real(wp), save :: dz
  real(wp), save :: tchar
  real(wp), save :: cm1,cm2,cm3,cm4,cm5,cm6
  !
  real(wp), save, allocatable :: cmbc(:,:)
  real(wp), save, allocatable :: by(:)
  real(wp), save, allocatable :: bz(:)
  real(wp), save, allocatable :: bycg(:,:,:)
  real(wp), save, allocatable :: bzcg(:,:,:)
  real(wp), save, allocatable :: xdf(:,:,:)
  real(wp), save, allocatable :: qdf(:,:,:)
  real(wp), save, allocatable :: phi(:,:,:)
  real(wp), save, allocatable :: phiold(:,:,:)
  real(wp), save, allocatable :: r(:,:,:)
  !
  integer, save, allocatable :: nfilycg(:,:)
  integer, save, allocatable :: nfilzcg(:,:)
  !
!###############################################################################
  !
  integer, save :: nit
  integer, save :: nb
 !integer, save :: ntbl
  integer, save :: iinit
  integer, save :: tininit
  !
  real(wp), save, allocatable :: ravgin(:)
  real(wp), save, allocatable :: uavgin(:)
  real(wp), save, allocatable :: tavgin(:)
  !
  real(wp), save, allocatable :: a11(:)
  real(wp), save, allocatable :: a21(:)
  real(wp), save, allocatable :: a22(:)
  real(wp), save, allocatable :: a31(:)
  real(wp), save, allocatable :: a32(:)
  real(wp), save, allocatable :: a33(:)
  !
  !
  real(wp), save, allocatable :: dt(:,:,:)
 !real(wp), save, allocatable :: amutbc(:,:,:)
  !
  real(wp), save, allocatable :: x(:,:,:,:)
  real(wp), save, allocatable :: q(:,:,:,:)
 !real(wp), save, allocatable :: qbc(:,:,:,:)
  !
  ! DF_IS_INITIALIZED : logical flag that indicates whether or not the DF
  !                     module has been initialized.
  !
  logical(lk), save :: df_is_initialized = fals
  !
contains
!
!###############################################################################
!
subroutine df(iside,icnst,min1,max1,min2,max2,uconv,lsx,lsy,lsz)
  !
  !***********************************************************************
  !
  ! Digital filtering Turbulent inflow routine.
  ! Ref: Xie & Castro, Flow Turbulence and Combustion, 2008, 81:449-470
  !
  ! This version performs the digital filtering on a uniform grid
  ! and the results are interpolated onto the computational grid.
  !
  ! Note: this only works for an i-min/constant-x plane
  !
  !***********************************************************************
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
 !use data
 !use dataturbin
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: iside
  integer,  intent(in) :: icnst
  integer,  intent(in) :: min1
  integer,  intent(in) :: max1
  integer,  intent(in) :: min2
  integer,  intent(in) :: max2
  real(wp), intent(in) :: uconv
  real(wp), intent(in) :: lsx
  real(wp), intent(in) :: lsy
  real(wp), intent(in) :: lsz
  !
  !.. Local Scalars ..
  integer  :: i,j,k,m,n,l,isize
  integer  :: jdf,kdf
  real(wp) :: rgas_over_gm1
  real(wp) :: facold,facnew,uinst,vinst,winst
  real(wp) :: y0,z0,y1,z1,y2,z2,c11,c12,c21,c22,c00,dtdf
  real(wp) :: ufluc,ufluc2,vx
  real(wp) :: vfluc,vfluc2,vy
  real(wp) :: wfluc,wfluc2,vz
  real(wp) :: rho,temp,q2,ravg,uavg,tavg
  !
continue
  !
  ! initialize on first iteration
  !
  if (.not. df_is_initialized) then
    call initialize_df(iside,icnst,min1,max1,min2,max2,uconv,lsx,lsy,lsz)
    df_is_initialized = true
  end if
  !
  rgas_over_gm1 = rgasref_nd / (gam - one)
  !
  ! Populate random numbers
  !
  call random_number(r)
  isize = (npydf+2*nfily+1)*(npzdf+2*nfilz+1)*3
  call boxmuller(r(-nfily,-nfilz,1),isize)
  !
  ! Put last time step into old
  !
  if (nit /= 1) then
    do k = 1,npzdf
      do j = 1,npydf
        phiold(j,k,1) = phi(j,k,1)
        phiold(j,k,2) = phi(j,k,2)
        phiold(j,k,3) = phi(j,k,3)
      end do
    end do
  end if
  !
  ! Filter random numbers
  !
  phi(:,:,:) = zero
  do k = 1,npzdf
    do j = 1,npydf
      do m = -nfilz,nfilz
        do l = -nfily,nfily
          phi(j,k,1) = phi(j,k,1) + by(l)*bz(m)*r(j+l,k+m,1)
          phi(j,k,2) = phi(j,k,2) + by(l)*bz(m)*r(j+l,k+m,2)
          phi(j,k,3) = phi(j,k,3) + by(l)*bz(m)*r(j+l,k+m,3)
        end do
      end do
    end do
  end do
  !
  ! Compute new phi based on old and new phi
  !
  dtdf = dt(1,1,1)
  !
  if (nit == 1 .and. iinit == 1) then
    facold = zero
    facnew = one
  else if (nit == 1 .and. tininit == 1) then
    facold = zero
    facnew = one
  else
    facold = exp(-pi*dtdf/(two*tchar))
    facnew = sqrt(one-exp(-pi*dtdf/tchar))
  end if
  !
  do k = 1,npzdf
    do j = 1,npydf
      phi(j,k,1) = phiold(j,k,1)*facold+phi(j,k,1)*facnew
      phi(j,k,2) = phiold(j,k,2)*facold+phi(j,k,2)*facnew
      phi(j,k,3) = phiold(j,k,3)*facold+phi(j,k,3)*facnew
    end do
  end do
  !
  ! project amplitude tensor onto phi & average df values
  !
  do k = 1,npzdf
    do j = 1,npydf
      !
      ufluc2 = phi(j,k,1)
      vfluc2 = phi(j,k,2)
      wfluc2 = phi(j,k,3)
      !
      ufluc = ufluc2*a11(j)
      vfluc = ufluc2*a21(j) + vfluc2*a22(j)
      wfluc = ufluc2*a31(j) + vfluc2*a32(j) + wfluc2*a33(j)
      !
      ravg = ravgin(j)
      uavg = uavgin(j)
      tavg = tavgin(j)
      !
      vx = ufluc + uavg
      vy = vfluc
      vz = wfluc
      !
      q2 = vx*vx + vy*vy + vz*vz
      !
      rho  = ravg
      temp = tavg
      !
      qdf(j,k,nec) = rho
      qdf(j,k,nmx) = rho*vx
      qdf(j,k,nmy) = rho*vy
      qdf(j,k,nmz) = rho*vz
      qdf(j,k,nee) = rho*( rgas_over_gm1*temp + half*q2 )
      !
    end do
  end do
  !
  ! interpolate from df grid to cfd grid
  !
  do k = min2,max2
    do j = min1,max1
      !
      y0  = x(icnst,j,k,2)
      z0  = x(icnst,j,k,3)
      !
      jdf = int((y0-ymin)/dy) + 1
      kdf = int((z0-zmin)/dz) + 1
      !
      y1  = xdf(jdf  ,kdf  ,1)
      z1  = xdf(jdf  ,kdf  ,2)
      y2  = xdf(jdf+1,kdf+1,1)
      z2  = xdf(jdf+1,kdf+1,2)
      !
      c00 = one/((y2-y1)*(z2-z1))
      !
      c11 = (y2-y0)*(z2-z0) * c00
      c21 = (y2-y0)*(z0-z1) * c00
      c12 = (y0-y1)*(z2-z0) * c00
      c22 = (y0-y1)*(z0-z1) * c00
      !
      q(icnst,j,k,nec) = qdf(jdf  ,kdf  ,nec)*c11 + &
                         qdf(jdf+1,kdf  ,nec)*c21 + &
                         qdf(jdf  ,kdf+1,nec)*c12 + &
                         qdf(jdf+1,kdf+1,nec)*c22
      !
      q(icnst,j,k,nmx) = qdf(jdf  ,kdf  ,nmx)*c11 + &
                         qdf(jdf+1,kdf  ,nmx)*c21 + &
                         qdf(jdf  ,kdf+1,nmx)*c12 + &
                         qdf(jdf+1,kdf+1,nmx)*c22
      !
      q(icnst,j,k,nmy) = qdf(jdf  ,kdf  ,nmy)*c11 + &
                         qdf(jdf+1,kdf  ,nmy)*c21 + &
                         qdf(jdf  ,kdf+1,nmy)*c12 + &
                         qdf(jdf+1,kdf+1,nmy)*c22
      !
      q(icnst,j,k,nmz) = qdf(jdf  ,kdf  ,nmz)*c11 + &
                         qdf(jdf+1,kdf  ,nmz)*c21 + &
                         qdf(jdf  ,kdf+1,nmz)*c12 + &
                         qdf(jdf+1,kdf+1,nmz)*c22
      !
      q(icnst,j,k,nee) = qdf(jdf  ,kdf  ,nee)*c11 + &
                         qdf(jdf+1,kdf  ,nee)*c21 + &
                         qdf(jdf  ,kdf+1,nee)*c12 + &
                         qdf(jdf+1,kdf+1,nee)*c22
      !
    end do
  end do
  !
end subroutine df
!
!###############################################################################
!
subroutine initialize_df(iside,icnst,min1,max1,min2,max2,uconv,lsx,lsy,lsz)
  !
  !***********************************************************************
  !
  ! Initializes Digital Filter on 1st iteration
  !
  ! Need
  ! lsx,lsy,lsz - length scales in all 3 directions
  ! uconv       - convective vel.
  ! npx,nyp     - dimensions of d.f. grid
  !
  !***********************************************************************
  !
  !.. Use Statements ..
  use random_mod, only : set_random_seed
 !use data
 !use dataturbin
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: iside
  integer,  intent(in) :: icnst
  integer,  intent(in) :: min1
  integer,  intent(in) :: max1
  integer,  intent(in) :: min2
  integer,  intent(in) :: max2
  real(wp), intent(in) :: uconv
  real(wp), intent(in) :: lsx
  real(wp), intent(in) :: lsy
  real(wp), intent(in) :: lsz
  !
  !.. Local Scalars ..
  integer  :: i,j,k,m,n,l
  integer  :: nlsx,nlsy,nlsz,npin,dfunit
  real(wp) :: yout,denom,facto
  real(wp) :: ydim,zdim
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: yin(:)
  real(wp), allocatable :: qin(:,:)
  real(wp), allocatable :: qout(:,:)
  !
continue
  !
  ! Initialize the random number generator
  !
  call set_random_seed
  !
  tchar = lsx/uconv
  !
  ! setup digital filtering grid
  !
  allocate ( xdf(1:npydf,1:npzdf,1:2) )
  allocate ( qout(1:npydf,1:7) )
  allocate ( qdf(1:npydf,1:npzdf,1:5) )
  !
  ymin = minval(x(icnst,min1:max1,min2:max2,2))
  ymax = maxval(x(icnst,min1:max1,min2:max2,2))
  zmin = minval(x(icnst,min1:max1,min2:max2,3))
  zmax = maxval(x(icnst,min1:max1,min2:max2,3))
  ydim = ymax - ymin
  zdim = zmax - zmin
  !
  dy = ydim/real(npydf-1,kind=wp)
  dz = zdim/real(npzdf-1,kind=wp)
  !
  ! extend domain 1/4 cell beyond boundaries
  !
  ydim = ydim + half*dy
  zdim = zdim + half*dz
  ymin = ymin - one4*dy
  ymax = ymax + one4*dy
  zmin = zmin - one4*dz
  zmax = zmax + one4*dz
  !
  ! recompute df grid spacings based on larger boundaries
  !
  dy = ydim/real(npydf-1,kind=wp)
  dz = zdim/real(npzdf-1,kind=wp)
  dx = dy
  !
  ! Generate grid for digital filtering
  !
  do k = 1,npzdf
    do j = 1,npydf
      xdf(j,k,1) = real(j-1,kind=wp)*dy + ymin ! y-coordinate
      xdf(j,k,2) = real(k-1,kind=wp)*dz + zmin ! z-coordinate
    end do
  end do
  !
  nlsx   = nint(lsx/dx)
  nfilx  = two*nlsx
  nlsy   = nint(lsy/dy)
  nfily  = two*nlsy
  nlsz   = nint(lsz/dz)
  nfilz  = two*nlsz
  !
  ! Read specified velocity and Reynolds stress profiles
  ! qin, 1 = u_avg
  !      2 = <uu>
  !      3 = <vv>
  !      4 = <ww>
  !      5 = <uv>
  !      6 = rho_avg
  !      7 = T_avg
  ! note: dimensional quantities!
  !
  dfunit = 100 + nb
  read (dfunit,*) npin
  !
  allocate (yin(npin),qin(npin,7))
  do j = 1,npin
    read (dfunit,*) yin(j),(qin(j,i), i=1,7)
   !read (dfunit,*) yin(j),qin(j,1),qin(j,2),qin(j,3), &
   !              qin(j,4),qin(j,5),qin(j,6),qin(j,7)
  end do
  !
  ! interpolate inputs onto the digital filter grid
  !
  m = 1
  do n = 2,npydf
    yout = xdf(n,1,1)
    while_loop: &
    do while (m <= npin)
      if (yout >= yin(m) .and. yout < yin(m+1)) then
        denom = yin(m+1)-yin(m)
        facto = yout - yin(m)
        do l = 1,7
          qout(n,l) = (qin(m+1,l)-qin(m,l))/denom*facto + qin(m,l)
        end do
       !goto 20
        exit while_loop
      end if
      m = m + 1
    end do while_loop
   !20 continue
  end do
  !
  ! Filter Coefficients
  allocate ( by(-nfily:nfily) )
  allocate ( bz(-nfilz:nfilz) )
  !
  ! Random numbers
  allocate ( r(-nfily:npydf+nfily,-nfilz:npzdf+nfilz,1:3) )
  !
  ! Filtered function
  allocate ( phi(1:npydf,1:npzdf,1:3) )
  !
  ! Amplitude tensor
  allocate ( a11(1:npydf) )
  allocate ( a21(1:npydf) )
  allocate ( a22(1:npydf) )
  allocate ( a33(1:npydf) )
  !
  ! Mean inflow values
  allocate ( uavgin(1:npydf) )
  allocate ( ravgin(1:npydf) )
  allocate ( tavgin(1:npydf) )
  !
  ! Filtered function @ previous time
  ! allocated in solnm/w if from a restart
  if (iinit == 1 .or. tininit == 1) then
    allocate ( phiold(1:npydf,1:npzdf,1:3) )
  end if
  !
  phi(    :,:,:) = zero
  if (iinit == 1 .or. tininit == 1) then
    phiold( :,:,:) = zero
  end if
  !
  ! Compute filter coefficients
  !
  ! y-direction
  !
  denom = zero
  do k = -nfily,nfily
    by(k) = dffilter(k,nlsy)
    denom = denom + by(k)**2
  end do
  denom = sqrt(denom)
  do k = -nfily,nfily
    by(k) = by(k)/denom
  end do
  !
  ! z-direction
  !
  denom = zero
  do k = -nfilz,nfilz
    bz(k) = dffilter(k,nlsz)
    denom = denom + bz(k)**2
  end do
  denom = sqrt(denom)
  do k = -nfilz,nfilz
    bz(k) = bz(k)/denom
  end do
  !
  ! Compute amplitude tensor
  !
  do j = 1,npydf
    a11(j) = sqrt(qout(j,2))
    a21(j) = qout(j,5)/(a11(j)+rndoff)
    a22(j) = sqrt(qout(j,3)-a21(j)**2)
    a33(j) = sqrt(qout(j,4))
    uavgin(j) = qout(j,1)
    ravgin(j) = qout(j,6)
   !ravgin(j) = 2.377e-03_wp
    tavgin(j) = qout(j,7)
   !tavgin(j) = 519.0_wp
  end do
  !
  deallocate (qin,qout,yin)
  !
  write (iout,*)
  write (iout,55) nb
  write (iout,*) 'Nlines = ',npin
  write (iout,*) 'ymin   = ',ymin
  write (iout,*) 'ymax   = ',ymax
  write (iout,*) 'zmin   = ',zmin
  write (iout,*) 'zmax   = ',zmax
  write (iout,*) 'dx     = ',dx
  write (iout,*) 'dy     = ',dy
  write (iout,*) 'dz     = ',dz
  write (iout,*) 'lsx    = ',lsx
  write (iout,*) 'uconv  = ',uconv
  write (iout,*) 'tchar  = ',tchar
  write (iout,*) 'dx     = ',dx
  write (iout,*) 'nlsx   = ',nlsx
  write (iout,*) 'nlsy   = ',nlsy
  write (iout,*) 'nlsz   = ',nlsz
  write (iout,*) 'nfilx  = ',nfilx
  write (iout,*) 'nfily  = ',nfily
  write (iout,*) 'nfilz  = ',nfilz
  write (iout,*)
  !
  55 format(1x,'Digital Filtering Initialization, Blk = ',i4,/, &
            1x,'--------------------------------------------')
  56 format(1x,'Lines in vel and Re stress file = ',i4)
  !
end subroutine initialize_df
!
!###############################################################################
!
subroutine dfcg(iside,icnst,min1,max1,min2,max2,uconv,lsx,lsy,lsz)
  !
  !***********************************************************************
  !
  ! Digital Filtering Turbulent Inflow
  ! Ref: Xie & Castro, Flow Turbulence and Combustion, 2008, 81:449-470
  !
  ! This version (cg) applies DF on the computational grid.  The cell
  ! cellsizes at each grid point are computed and used to scale the
  ! filter width to provide the correct lenght scales.
  !
  ! Note: this only works for an i-min/constant-x plane
  !
  !***********************************************************************
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
 !use data
 !use dataturbin
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: iside
  integer,  intent(in) :: icnst
  integer,  intent(in) :: min1
  integer,  intent(in) :: max1
  integer,  intent(in) :: min2
  integer,  intent(in) :: max2
  real(wp), intent(in) :: uconv
  real(wp), intent(in) :: lsx
  real(wp), intent(in) :: lsy
  real(wp), intent(in) :: lsz
  !
  !.. Local Scalars ..
  integer  :: i,j,k,m,n,l,isize
  integer  :: jdf,kdf
  real(wp) :: rgas_over_gm1
  real(wp) :: facold,facnew,uinst,vinst,winst
  real(wp) :: y0,z0,y1,z1,y2,z2,c11,c12,c21,c22,c00,dtdf
  !
continue
  !
  ! initialize on first iteration
  !
  if (.not. df_is_initialized) then
    call initialize_dfcg(iside,icnst,min1,max1,min2,max2,uconv,lsx,lsy,lsz)
    df_is_initialized = true
  end if
  !
  rgas_over_gm1 = rgasref_nd / (gam - one)
  !
  ! Populate random numbers
  !
  call random_number(r)
  isize = 3 * ((max1+mxnfily)-(min1-mxnfily)+1) &
            * ((max2+mxnfilz)-(min2-mxnfilz)+1)
  call boxmuller(r(min1-mxnfily,min2-mxnfilz,1),isize)
  !
  ! Put last time step into old
  !
  if (nit /= 1) then
    do k = min2,max2
      do j = min1,max1
        phiold(j,k,1) = phi(j,k,1)
        phiold(j,k,2) = phi(j,k,2)
        phiold(j,k,3) = phi(j,k,3)
      end do
    end do
  end if
  !
  ! Filter random numbers
  !
  phi(:,:,:) = zero
  do k = min2,max2
    do j = min1,max1
      do m = -nfilzcg(j,k),nfilzcg(j,k)
        do l = -nfilycg(j,k),nfilycg(j,k)
          phi(j,k,1) = phi(j,k,1) + bycg(j,k,l)*bzcg(j,k,m)*r(j+l,k+m,1)
          phi(j,k,2) = phi(j,k,2) + bycg(j,k,l)*bzcg(j,k,m)*r(j+l,k+m,2)
          phi(j,k,3) = phi(j,k,3) + bycg(j,k,l)*bzcg(j,k,m)*r(j+l,k+m,3)
        end do
      end do
    end do
  end do
  !
  ! Compute new phi based on old and new phi
  !
  dtdf = dt(1,1,1)
  !
  if (nit == 1 .and. iinit == 1) then
    facold = zero
    facnew = one
  else if (nit == 1 .and. tininit == 1) then
    facold = zero
    facnew = one
  else
    facold = exp(-pi*dtdf/(two*tchar))
    facnew = sqrt(one-exp(-pi*dtdf/tchar))
  end if
  !
  do k = min2,max2
    do j = min1,max1
      phi(j,k,1) = phiold(j,k,1)*facold+phi(j,k,1)*facnew
      phi(j,k,2) = phiold(j,k,2)*facold+phi(j,k,2)*facnew
      phi(j,k,3) = phiold(j,k,3)*facold+phi(j,k,3)*facnew
    end do
  end do
  !
  ! project amplitude tensor onto phi & average df values
  !
  do k = min2,max2
    do j = min1,max1
      uinst = uavgin(j)+a11(j)*phi(j,k,1)
      vinst = a21(j)*phi(j,k,1) + a22(j)*phi(j,k,2)
      winst = a33(j)*phi(j,k,3)
      q(icnst,j,k,1) = ravgin(j)
      q(icnst,j,k,2) = ravgin(j)*uinst
      q(icnst,j,k,3) = ravgin(j)*vinst
      q(icnst,j,k,4) = ravgin(j)*winst
      q(icnst,j,k,5) = ravgin(j)*(rgas_over_gm1*tavgin(j) + &
                       half*(uinst**2+ vinst**2 + winst**2))
    end do
  end do
  !
end subroutine dfcg
!
!###############################################################################
!
subroutine initialize_dfcg(iside,icnst,min1,max1,min2,max2,uconv,lsx,lsy,lsz)
  !
  !***********************************************************************
  !
  ! Initializes Digital Filter on 1st iteration
  ! Version for filtering on computational grid
  !
  ! Need
  ! lsx,lsy,lsz - length scales in all 3 directions
  ! uconv       - convective vel.
  !
  !***********************************************************************
  !
  !.. Use Statements ..
  use random_mod, only : set_random_seed
 !use data
 !use dataturbin
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: iside
  integer,  intent(in) :: icnst
  integer,  intent(in) :: min1
  integer,  intent(in) :: max1
  integer,  intent(in) :: min2
  integer,  intent(in) :: max2
  real(wp), intent(in) :: uconv
  real(wp), intent(in) :: lsx
  real(wp), intent(in) :: lsy
  real(wp), intent(in) :: lsz
  !
  !.. Local Scalars ..
  integer  :: i,j,k,m,n,l
  integer  :: npin,dfunit
  real(wp) :: yout,denom,facto
  real(wp) :: ydim,zdim
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: nlsy(:,:)
  integer,  allocatable :: nlsz(:,:)
  real(wp), allocatable :: yin(:)
  real(wp), allocatable :: qin(:,:)
  real(wp), allocatable :: dely(:,:)
  real(wp), allocatable :: delz(:,:)
  !
continue
  !
  ! Initialize the random number generator
  !
  call set_random_seed
  !
  allocate (    dely(min1:max1,min2:max2) )
  allocate (    delz(min1:max1,min2:max2) )
  allocate (    nlsy(min1:max1,min2:max2) )
  allocate (    nlsz(min1:max1,min2:max2) )
  allocate ( nfilycg(min1:max1,min2:max2) )
  allocate ( nfilzcg(min1:max1,min2:max2) )
  !
  tchar = lsx/uconv
  !
  ! compute local grid spacing
  !
  call dfcellsize(icnst,min1,max1,min2,max2,dely,delz)
  !
  ! use grid spacing to scale filter widths
  !
  do k = min2,max2
    do j = min1,max1
      nlsy( j,k) = nint(lsy/dely(j,k))
      nfilycg(j,k) = two*nlsy(j,k)
      nlsz( j,k) = nint(lsz/delz(j,k))
      nfilzcg(j,k) = two*nlsz(j,k)
    end do
  end do
  !
  mxnfily = maxval(nfilycg)         ! max. filter width in y
  mxnfilz = maxval(nfilzcg)         ! max. filter width in z
  !
  npydf = max1-min1+1               ! define these variables for use
  npzdf = max2-min2+1               ! in DF i/o
  !
  ! Read specified velocity and Reynolds stress profiles
  ! qin, 1 = u_avg
  !      2 = <uu>
  !      3 = <vv>
  !      4 = <ww>
  !      5 = <uv>
  !      6 = rho_avg
  !      7 = T_avg
  ! note: dimensional quantities!
  !
  dfunit = 100 + nb
  allocate ( yin(min1:max1) )
  allocate ( qin(min1:max1,1:7) )
  !
  do j = min1,max1
    read (dfunit,*) yin(j),(qin(j,i), i=1,7)
   !read (dfunit,*) yin(j),qin(j,1),qin(j,2),qin(j,3), &
   !              qin(j,4),qin(j,5),qin(j,6),qin(j,7)
  end do
  !
  ! Filter Coefficients
  allocate ( bycg(min1:max1,min2:max2,-mxnfily:mxnfily) )
  allocate ( bzcg(min1:max1,min2:max2,-mxnfilz:mxnfilz) )
  !
  ! Random numbers
  allocate ( r(min1-mxnfily:max1+mxnfily,min2-mxnfilz:max2+mxnfilz,1:3))
  !
  ! Filtered function
  allocate ( phi(min1:max1,min2:max2,1:3))
  !
  ! Amplitude tensor
  allocate ( a11(min1:max1) )
  allocate ( a21(min1:max1) )
  allocate ( a22(min1:max1) )
  allocate ( a33(min1:max1) )
  !
  ! Mean inflow values
  allocate ( uavgin(min1:max1) )
  allocate ( ravgin(min1:max1) )
  allocate ( tavgin(min1:max1) )
  !
  if (iinit == 1 .or. tininit == 1) then
    ! Filtered function @ previous time
    ! allocated in solnm/w if from a restart
    allocate ( phiold(min1:max1,min2:max2,1:3) )
  end if
  !
  phi(:,:,:) = zero
  ! initialize if not read in
  if (iinit == 1 .or. tininit == 1) then
    phiold(:,:,:) = zero
  end if
  !
  ! Compute filter coefficients
  !
  ! y-direction
  !
  do k = min2,max2
    do j = min1,max1
      denom = zero
      do l = -nfilycg(j,k),nfilycg(j,k)
        bycg(j,k,l) = dffilter(l,nlsy(j,k))
        denom = denom + bycg(j,k,l)**2
      end do
      denom = sqrt(denom)
      do l = -nfilycg(j,k),nfilycg(j,k)
        bycg(j,k,l) = bycg(j,k,l)/denom
      end do
    end do
  end do
  !
  ! z-direction
  !
  do k = min2,max2
    do j = min1,max1
      denom = zero
      do l = -nfilzcg(j,k),nfilzcg(j,k)
        bzcg(j,k,l) = dffilter(l,nlsz(j,k))
        denom = denom + bzcg(j,k,l)**2
      end do
      denom = sqrt(denom)
      do l = -nfilzcg(j,k),nfilzcg(j,k)
        bzcg(j,k,l) = bzcg(j,k,l)/denom
      end do
    end do
  end do
  !
  ! Compute amplitude tensor
  !
  do j = min1,max1
    a11(j) = sqrt(qin(j,2))
    a21(j) = qin(j,5)/(a11(j)+rndoff)
    a22(j) = sqrt(qin(j,3)-a21(j)**2)
    a33(j) = sqrt(qin(j,4))
    uavgin(j) = qin(j,1)
    ravgin(j) = qin(j,6)
   !ravgin(j) = 2.377e-03_wp
    tavgin(j) = qin(j,7)
   !tavgin(j) = 519.0_wp
  end do
  !
  deallocate ( yin )
  !
  write (iout,*)
  write (iout,55) nb
  write (iout,*) 'lsx    = ',lsx
  write (iout,*) 'lsy    = ',lsy
  write (iout,*) 'lsz    = ',lsz
  write (iout,*) 'uconv  = ',uconv
  write (iout,*) 'tchar  = ',tchar
  write (iout,*) 'mxnfily  = ',mxnfily
  write (iout,*) 'mxnfilz  = ',mxnfilz
  !
  55 format(1x,'Digital Filtering Initialization, Blk = ',i4,/, &
            1x,'--------------------------------------------')
  !
end subroutine initialize_dfcg
!
!###############################################################################
!
subroutine dfcellsize(i,min1,max1,min2,max2,dely,delz)
  !
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! This routine computes the cell size for the dfcg routines
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
 !!.. Use Statements ..
 !use data
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: i
  integer,     intent(in) :: min1
  integer,     intent(in) :: max1
  integer,     intent(in) :: min2
  integer,     intent(in) :: max2
  real(wp), intent(inout) :: dely(min1:max1,min2:max2)
  real(wp), intent(inout) :: delz(min1:max1,min2:max2)
  !
  !.. Local Scalars ..
  real(wp) :: dxi12,det12,dze12
  integer  :: j,k,nbc,nhol,size2,ibcs,ic
  !
  !.. Local Arrays ..
  real(wp) :: yet(min1:max1,min2:max2)
  real(wp) :: yze(min1:max1,min2:max2)
  real(wp) :: zet(min1:max1,min2:max2)
  real(wp) :: zze(min1:max1,min2:max2)
  real(wp) :: dxi(1:3)
  real(wp) :: rdxi(1:3)
  !
continue
  !
  dxi(1) = one
  dxi(2) = one
  dxi(3) = one
  rdxi(1) = one/dxi(1)
  rdxi(2) = one/dxi(2)
  rdxi(3) = one/dxi(3)
  dxi12 = twelve*dxi(1)
  det12 = twelve*dxi(2)
  dze12 = twelve*dxi(3)
  !
  ! eta direction
  !
  do k = min2,max2
    do j = min1+6,max1-6
      yet(j,k) = cm1*(x(i,j+1,k,2)-x(i,j-1,k,2)) + &
                 cm2*(x(i,j+2,k,2)-x(i,j-2,k,2)) + &
                 cm3*(x(i,j+3,k,2)-x(i,j-3,k,2)) + &
                 cm4*(x(i,j+4,k,2)-x(i,j-4,k,2)) + &
                 cm5*(x(i,j+5,k,2)-x(i,j-5,k,2)) + &
                 cm6*(x(i,j+6,k,2)-x(i,j-6,k,2))
      zet(j,k) = cm1*(x(i,j+1,k,3)-x(i,j-1,k,3)) + &
                 cm2*(x(i,j+2,k,3)-x(i,j-2,k,3)) + &
                 cm3*(x(i,j+3,k,3)-x(i,j-3,k,3)) + &
                 cm4*(x(i,j+4,k,3)-x(i,j-4,k,3)) + &
                 cm5*(x(i,j+5,k,3)-x(i,j-5,k,3)) + &
                 cm6*(x(i,j+6,k,3)-x(i,j-6,k,3))
    end do
  end do
  !
  ! zeta direction
  !
  do k = min2+6,max2-6
    do j = min1,max1
      yze(j,k) = cm1*(x(i,j,k+1,2)-x(i,j,k-1,2)) + &
                 cm2*(x(i,j,k+2,2)-x(i,j,k-2,2)) + &
                 cm3*(x(i,j,k+3,2)-x(i,j,k-3,2)) + &
                 cm4*(x(i,j,k+4,2)-x(i,j,k-4,2)) + &
                 cm5*(x(i,j,k+5,2)-x(i,j,k-5,2)) + &
                 cm6*(x(i,j,k+6,2)-x(i,j,k-6,2))
      zze(j,k) = cm1*(x(i,j,k+1,3)-x(i,j,k-1,3)) + &
                 cm2*(x(i,j,k+2,3)-x(i,j,k-2,3)) + &
                 cm3*(x(i,j,k+3,3)-x(i,j,k-3,3)) + &
                 cm4*(x(i,j,k+4,3)-x(i,j,k-4,3)) + &
                 cm5*(x(i,j,k+5,3)-x(i,j,k-5,3)) + &
                 cm6*(x(i,j,k+6,3)-x(i,j,k-6,3))
    end do
  end do
  !
  ! jmin boundary
  !
  do k = min2,max2
    ibcs = 1
    do j = min1,min1+5
      ic = min1
      yet(j,k) = cmbc(13,ibcs)*x(i,ic+12,k,2) + &
                 cmbc(12,ibcs)*x(i,ic+11,k,2) + &
                 cmbc(11,ibcs)*x(i,ic+10,k,2) + &
                 cmbc(10,ibcs)*x(i,ic+ 9,k,2) + &
                 cmbc(9 ,ibcs)*x(i,ic+ 8,k,2) + &
                 cmbc(8 ,ibcs)*x(i,ic+ 7,k,2) + &
                 cmbc(7 ,ibcs)*x(i,ic+ 6,k,2) + &
                 cmbc(6 ,ibcs)*x(i,ic+ 5,k,2) + &
                 cmbc(5 ,ibcs)*x(i,ic+ 4,k,2) + &
                 cmbc(4 ,ibcs)*x(i,ic+ 3,k,2) + &
                 cmbc(3 ,ibcs)*x(i,ic+ 2,k,2) + &
                 cmbc(2 ,ibcs)*x(i,ic+ 1,k,2) + &
                 cmbc(1 ,ibcs)*x(i,ic   ,k,2)
      zet(j,k) = cmbc(13,ibcs)*x(i,ic+12,k,3) + &
                 cmbc(12,ibcs)*x(i,ic+11,k,3) + &
                 cmbc(11,ibcs)*x(i,ic+10,k,3) + &
                 cmbc(10,ibcs)*x(i,ic+ 9,k,3) + &
                 cmbc(9 ,ibcs)*x(i,ic+ 8,k,3) + &
                 cmbc(8 ,ibcs)*x(i,ic+ 7,k,3) + &
                 cmbc(7 ,ibcs)*x(i,ic+ 6,k,3) + &
                 cmbc(6 ,ibcs)*x(i,ic+ 5,k,3) + &
                 cmbc(5 ,ibcs)*x(i,ic+ 4,k,3) + &
                 cmbc(4 ,ibcs)*x(i,ic+ 3,k,3) + &
                 cmbc(3 ,ibcs)*x(i,ic+ 2,k,3) + &
                 cmbc(2 ,ibcs)*x(i,ic+ 1,k,3) + &
                 cmbc(1 ,ibcs)*x(i,ic   ,k,3)
      ibcs = ibcs + 1
    end do
  end do
  !
  ! jmax boundary
  !
  do k = min2,max2
    ibcs = 1
    do j = max1,max1-5,-1
      ic = max1
      yet(j,k) = -one*(cmbc(13,ibcs)*x(i,ic-12,k,2) + &
                       cmbc(12,ibcs)*x(i,ic-11,k,2) + &
                       cmbc(11,ibcs)*x(i,ic-10,k,2) + &
                       cmbc(10,ibcs)*x(i,ic- 9,k,2) + &
                       cmbc(9 ,ibcs)*x(i,ic- 8,k,2) + &
                       cmbc(8 ,ibcs)*x(i,ic- 7,k,2) + &
                       cmbc(7 ,ibcs)*x(i,ic- 6,k,2) + &
                       cmbc(6 ,ibcs)*x(i,ic- 5,k,2) + &
                       cmbc(5 ,ibcs)*x(i,ic- 4,k,2) + &
                       cmbc(4 ,ibcs)*x(i,ic- 3,k,2) + &
                       cmbc(3 ,ibcs)*x(i,ic- 2,k,2) + &
                       cmbc(2 ,ibcs)*x(i,ic- 1,k,2) + &
                       cmbc(1 ,ibcs)*x(i,ic   ,k,2))
      zet(j,k) = -one*(cmbc(13,ibcs)*x(i,ic-12,k,3) + &
                       cmbc(12,ibcs)*x(i,ic-11,k,3) + &
                       cmbc(11,ibcs)*x(i,ic-10,k,3) + &
                       cmbc(10,ibcs)*x(i,ic- 9,k,3) + &
                       cmbc(9 ,ibcs)*x(i,ic- 8,k,3) + &
                       cmbc(8 ,ibcs)*x(i,ic- 7,k,3) + &
                       cmbc(7 ,ibcs)*x(i,ic- 6,k,3) + &
                       cmbc(6 ,ibcs)*x(i,ic- 5,k,3) + &
                       cmbc(5 ,ibcs)*x(i,ic- 4,k,3) + &
                       cmbc(4 ,ibcs)*x(i,ic- 3,k,3) + &
                       cmbc(3 ,ibcs)*x(i,ic- 2,k,3) + &
                       cmbc(2 ,ibcs)*x(i,ic- 1,k,3) + &
                       cmbc(1 ,ibcs)*x(i,ic   ,k,3))
      ibcs = ibcs + 1
    end do
  end do
  !
  ! kmin boundary
  !
  do j = min1,max1
    ibcs = 1
    do k = min2,min2+5
      ic = min2
      yze(j,k) = cmbc(13,ibcs)*x(i,j,ic+12,2) + &
                 cmbc(12,ibcs)*x(i,j,ic+11,2) + &
                 cmbc(11,ibcs)*x(i,j,ic+10,2) + &
                 cmbc(10,ibcs)*x(i,j,ic+ 9,2) + &
                 cmbc(9 ,ibcs)*x(i,j,ic+ 8,2) + &
                 cmbc(8 ,ibcs)*x(i,j,ic+ 7,2) + &
                 cmbc(7 ,ibcs)*x(i,j,ic+ 6,2) + &
                 cmbc(6 ,ibcs)*x(i,j,ic+ 5,2) + &
                 cmbc(5 ,ibcs)*x(i,j,ic+ 4,2) + &
                 cmbc(4 ,ibcs)*x(i,j,ic+ 3,2) + &
                 cmbc(3 ,ibcs)*x(i,j,ic+ 2,2) + &
                 cmbc(2 ,ibcs)*x(i,j,ic+ 1,2) + &
                 cmbc(1 ,ibcs)*x(i,j,ic   ,2)
      zze(j,k) = cmbc(13,ibcs)*x(i,j,ic+12,3) + &
                 cmbc(12,ibcs)*x(i,j,ic+11,3) + &
                 cmbc(11,ibcs)*x(i,j,ic+10,3) + &
                 cmbc(10,ibcs)*x(i,j,ic+ 9,3) + &
                 cmbc(9 ,ibcs)*x(i,j,ic+ 8,3) + &
                 cmbc(8 ,ibcs)*x(i,j,ic+ 7,3) + &
                 cmbc(7 ,ibcs)*x(i,j,ic+ 6,3) + &
                 cmbc(6 ,ibcs)*x(i,j,ic+ 5,3) + &
                 cmbc(5 ,ibcs)*x(i,j,ic+ 4,3) + &
                 cmbc(4 ,ibcs)*x(i,j,ic+ 3,3) + &
                 cmbc(3 ,ibcs)*x(i,j,ic+ 2,3) + &
                 cmbc(2 ,ibcs)*x(i,j,ic+ 1,3) + &
                 cmbc(1 ,ibcs)*x(i,j,ic   ,3)
      ibcs = ibcs + 1
    end do
  end do
  !
  ! kmax boundary
  !
  do j = min1,max1
    ibcs = 1
    do k = max2,max2-5,-1
      ic = max2
      yze(j,k) = -one*(cmbc(13,ibcs)*x(i,j,ic-12,2) + &
                       cmbc(12,ibcs)*x(i,j,ic-11,2) + &
                       cmbc(11,ibcs)*x(i,j,ic-10,2) + &
                       cmbc(10,ibcs)*x(i,j,ic- 9,2) + &
                       cmbc(9 ,ibcs)*x(i,j,ic- 8,2) + &
                       cmbc(8 ,ibcs)*x(i,j,ic- 7,2) + &
                       cmbc(7 ,ibcs)*x(i,j,ic- 6,2) + &
                       cmbc(6 ,ibcs)*x(i,j,ic- 5,2) + &
                       cmbc(5 ,ibcs)*x(i,j,ic- 4,2) + &
                       cmbc(4 ,ibcs)*x(i,j,ic- 3,2) + &
                       cmbc(3 ,ibcs)*x(i,j,ic- 2,2) + &
                       cmbc(2 ,ibcs)*x(i,j,ic- 1,2) + &
                       cmbc(1 ,ibcs)*x(i,j,ic   ,2))
      zze(j,k) = -one*(cmbc(13,ibcs)*x(i,j,ic-12,3) + &
                       cmbc(12,ibcs)*x(i,j,ic-11,3) + &
                       cmbc(11,ibcs)*x(i,j,ic-10,3) + &
                       cmbc(10,ibcs)*x(i,j,ic- 9,3) + &
                       cmbc(9 ,ibcs)*x(i,j,ic- 8,3) + &
                       cmbc(8 ,ibcs)*x(i,j,ic- 7,3) + &
                       cmbc(7 ,ibcs)*x(i,j,ic- 6,3) + &
                       cmbc(6 ,ibcs)*x(i,j,ic- 5,3) + &
                       cmbc(5 ,ibcs)*x(i,j,ic- 4,3) + &
                       cmbc(4 ,ibcs)*x(i,j,ic- 3,3) + &
                       cmbc(3 ,ibcs)*x(i,j,ic- 2,3) + &
                       cmbc(2 ,ibcs)*x(i,j,ic- 1,3) + &
                       cmbc(1 ,ibcs)*x(i,j,ic   ,3))
      ibcs = ibcs + 1
    end do
  end do
  !
  ! compute grid spacings
  !
  do k = min2,max2
    do j = min1,max1
      dely(j,k) = yet(j,k)*dxi(2) + yze(j,k)*dxi(3)
      delz(j,k) = zet(j,k)*dxi(2) + zze(j,k)*dxi(3)
    end do
  end do
  !
 !deallocate (yet,yze,zet,zze)
  !
end subroutine dfcellsize
!
!###############################################################################
!
pure &
function dffilter(a,b) result(return_value)
  !
  !***********************************************************************
  ! Filter Function
  !***********************************************************************
  !
  !.. Formal Arguments ..
  integer, intent(in) :: a,b
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
continue
  !
  return_value = exp(-pi*abs(real(a,kind=wp))/real(b,kind=wp))
  !
end function dffilter
!
!###############################################################################
!
pure &
subroutine boxmuller(a,n)
  !
  !***********************************************************************
  ! Box-Muller: Computes a normally distributed set of random numbers
  !             from a uniformly distributed set.
  !***********************************************************************
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n
  real(wp), intent(inout) :: a(1:n+1)
  !
  !.. Local Scalars ..
  integer  :: nmax,i
  real(wp) :: u1,u2
  !
continue
  !
  if (mod(n,2) /= 0) then
    a(n+1) = a(1)
    nmax = n
  else
    nmax = n-1
  end if
  !
  do i = 1,nmax,2
    u1     = a(i)
    u2     = a(i+1)
    a(i)   = sqrt(-two*log(u1))*cos(two*pi*u2)
    a(i+1) = sqrt(-two*log(u1))*sin(two*pi*u2)
  end do
  !
end subroutine boxmuller
!
!###############################################################################
!
end module digital_filtering_mod
!
!###############################################################################
!
