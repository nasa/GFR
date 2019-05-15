module synthetic_eddy_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nee,nq,ntemp
  !
  implicit none
  !
  private
  !
  public :: sem
  public :: separate_sem_boundary_faces
  public :: output_sem_averages
  public :: dump_sem_data
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
  type :: eddy_list_t
    integer :: num_eddies = 0
    integer, allocatable :: eddy_list(:)
  end type eddy_list_t
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
    !     These components are to eliminate the need for storing the array
    !     facexyz_dt or facexyz_array, because those are only really needed
    !     for MMS [facexyz_array(1:nr,1:n_bnd_flxpts)] or for the continuous
    !     output feature [facexyz_array(1:nr,1:n_flxpts)]. For SEM, we only
    !     need the subset of (1:n_bnd_flxpts) that are the flux/face points
    !     for the SEM boundary conditions, so we dont want to waste memory
    !     for the flux/face points for which we will never require their xyz
    !     coordinates.
    real(wp) :: x,y,z
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
    ! WALL_DISTANCE : distance of this SEM point from the nearest wall
    real(wp) :: wall_distance
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
    ! IS_CYLINDRICAL : logical flag indicating if the box extents are
    !                  best defined in cylindrical coordinates
    logical(lk) :: is_cylindrical = fals
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
  integer, parameter :: max_eddies = 5000
  !
  ! TREAT EDDY VOLUMES AS SPHERE-LIKE
 !logical(lk), parameter :: eddy_volume_is_spherical = true
  ! TREAT EDDY VOLUMES AS BOX-LIKE
  logical(lk), parameter :: eddy_volume_is_spherical = fals
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
  type(sem_pts_t), save, allocatable :: my_sem_pts(:)
  !
  ! XYZ_WALL : Array to store the xyz coordinates for all grid nodes located
  !            on a wall boundary condition
  !
  real(wp), save, allocatable :: xyz_wall(:,:)
  !
contains
!
!###############################################################################
!
subroutine sem(cv)
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
  use ovar, only : d_t,gam,rgasref_nd
  use ovar, only : use_improved_SEM_normalization
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: cv(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,ierr,sem_idx,eddy_idx,num_eddies
  real(wp) :: ave,running_ave,time_scale,dt
  real(wp) :: c0,c1,eddy_vol,boxvol,lturb,dx,fsig,rgas_over_gm1
  real(wp) :: a11,a21,a22,a31,a32,a33,rho,temp,ravg,tavg,q2
  real(wp) :: x,epx,elx,eix,tfx,ufluc,ufluc2,uavg,vx
  real(wp) :: y,epy,ely,eiy,tfy,vfluc,vfluc2,vavg,vy
  real(wp) :: z,epz,elz,eiz,tfz,wfluc,wfluc2,wavg,vz
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
    !
  end if
  !
  rgas_over_gm1 = rgasref_nd / (gam - one)
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
      rho  = ravg
      temp = tavg
     !rho  = one  ! = rhoref / rhoref
     !temp = one  ! = tref / tref
      !
      cv(k,nec) = rho
      cv(k,nmx) = rho*vx
      cv(k,nmy) = rho*vy
      cv(k,nmz) = rho*vz
      cv(k,nee) = rho*( rgas_over_gm1*temp + half*q2 )
      !
#ifdef ENABLE_OUTPUT_TIME_HISTORY
      my_sem_pts(sem_idx)%sem_fp(n)%time_history(nec ,itcur) = rho
      my_sem_pts(sem_idx)%sem_fp(n)%time_history(nmx ,itcur) = vx
      my_sem_pts(sem_idx)%sem_fp(n)%time_history(nmy ,itcur) = vy
      my_sem_pts(sem_idx)%sem_fp(n)%time_history(nmz ,itcur) = vy
      my_sem_pts(sem_idx)%sem_fp(n)%time_history(nee ,itcur) = cv(k,nee)
      my_sem_pts(sem_idx)%sem_fp(n)%time_history(nq+1,itcur) = temp
#endif
      !
    end do loop_over_sem_pts
    !
  end do loop_over_my_sem_boxes
  !
  if (allocated(sem_avg)) then
    call update_sem_ave(dt,cv)
  end if
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
  call generate_new_eddies(new_eddies,generate_at_xmin=true)
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
subroutine update_sem_ave(dt,cv)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd,itcur
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: dt
  real(wp), intent(in) :: cv(:,:)
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
  if (sem_prev_ave_time < zero) then
    current_time = zero
    a = zero
    b = one
  else
    current_time = sem_prev_ave_time + dt
    a = sem_prev_ave_time / current_time
    b = one - a
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
subroutine generate_new_eddies(new_eddies,generate_at_xmin)
  !
  !.. Formal Arguments ..
  type(eddy_list_t), intent(in) :: new_eddies(:)
  logical(lk),       intent(in) :: generate_at_xmin
  !
  !.. Local Scalars ..
  integer  :: n,num_eddies,sem_idx,eddy_idx
  real(wp) :: eddy_vol,x_coef,vol_coef
  real(wp) :: x,eix,elx,boxminx,boxmaxx
  real(wp) :: y,eiy,ely,boxminy,boxmaxy
  real(wp) :: z,eiz,elz,boxminz,boxmaxz
  !
  !.. Local Arrays ..
  real(wp) :: new_xyz(1:3)
  real(wp) :: new_int(1:3)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "generate_new_eddies"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  x_coef = merge( zero , one , generate_at_xmin )
  !
  vol_coef = merge( (four/three)*pi , one , eddy_volume_is_spherical )
  !
  loop_over_sem_boxes: &
  do sem_idx = 1,num_sem_boxes
    !
    num_eddies = new_eddies(sem_idx)%num_eddies
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
    loop_over_new_eddies_in_box: &
    do n = 1,num_eddies
      !
      eddy_idx = new_eddies(sem_idx)%eddy_list(n)
      !
      ! Create new random values to get the x,y,z coordinates for the new eddy
      !
      call random_number(new_xyz)
      !
      x = boxminx + new_xyz(1)*(boxmaxx-boxminx)*x_coef
      y = boxminy + new_xyz(2)*(boxmaxy-boxminy)
      z = boxminz + new_xyz(3)*(boxmaxz-boxminz)
      !
      ! Save the new random coordinates for the new eddy in sem_box
      !
      if (sem_box(sem_idx)%is_cylindrical) then
        sem_box(sem_idx)%eddy(eddy_idx)%x = x
        sem_box(sem_idx)%eddy(eddy_idx)%y = y*cos(z) ! y
        sem_box(sem_idx)%eddy(eddy_idx)%z = y*sin(z) ! z
        sem_box(sem_idx)%eddy(eddy_idx)%r = y        ! r
        sem_box(sem_idx)%eddy(eddy_idx)%t = z        ! theta
      else
        sem_box(sem_idx)%eddy(eddy_idx)%x = x
        sem_box(sem_idx)%eddy(eddy_idx)%y = y          ! y
        sem_box(sem_idx)%eddy(eddy_idx)%z = z          ! z
        sem_box(sem_idx)%eddy(eddy_idx)%r = hypot(y,z) ! r
        sem_box(sem_idx)%eddy(eddy_idx)%t = atan2(z,y) ! theta
      end if
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
      ! NOTE: The constant eddy_volume_is_spherical is equal to either
      !       (4/3)*PI or 1 depending on whether the eddy is treated as
      !       sphere-like or box-like in volume, respectively.
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
  use ovar,       only : bc_in,average_sem_inflow_plane
  use geovar,     only : bc_grps
  use flowvar,    only : facexyz_array
  use random_mod, only : set_random_seed
  !
  !.. Local Scalars ..
  integer  :: ierr,k,n,np,n1,n2
  integer  :: sem_idx,grp_idx,bc_idx,num_eddies
  integer  :: num_restart_sem_boxes
  real(wp) :: lturb,uconv,x,y,z,theta
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
  call set_random_seed
  !
  ! Count the number of BC groups on this processor that are of the SEM type
  !
  num_restart_sem_boxes = 0
  num_sem_boxes = count( bc_in(:)%bc_type == bc_sem )
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
      sem_box(sem_idx)%bc_idx = bc_idx
      sem_box(sem_idx)%length_scale = lturb
      sem_box(sem_idx)%convective_velocity = uconv
      !
      sem_box(sem_idx)%is_cylindrical = bc_in(bc_idx)%sem%is_cylindrical
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
        if (y /= zero .or. z /= zero) then
          !
          theta = atan2( z , y )
          if (theta < zero) theta = theta + two*pi
          !
          my_sem_pts(sem_idx)%sem_fp(n)%sin_theta = sin(theta)
          my_sem_pts(sem_idx)%sem_fp(n)%cos_theta = cos(theta)
          !
        else
          !
          my_sem_pts(sem_idx)%sem_fp(n)%sin_theta = zero
          my_sem_pts(sem_idx)%sem_fp(n)%cos_theta = zero
          !
        end if
        !
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
  !   1) get the min/max extents of the box,
  !   2) compute the volume of the box, and
  !   3) compute the number of eddies to use for the box.
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
  call read_sem_profiles()
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
  call generate_new_eddies(new_eddies,generate_at_xmin=fals)
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
                                          sem_box(sem_idx)%profile)
      !
    end do
    !
  end do
  !
  ! Now that we have interpolated all the profile data to the SEM
  ! face points, we no longer need the actual profile data so deallocate
  ! all the sem_box(:)%profile components to free up that memory.
  !
  do sem_idx = 1,num_sem_boxes
    if (allocated(sem_box(sem_idx)%profile)) then
      deallocate ( sem_box(sem_idx)%profile , stat=ierr , errmsg=error_message )
      write (array_name,'("sem_box(",i0,")%profile")') sem_idx
      call alloc_error(pname,array_name,2,__LINE__,__FILE__,ierr,error_message)
    end if
  end do
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
 !real(wp), allocatable :: xyz_wall(:,:)
  real(wp), allocatable :: xyz_sem(:,:)
  real(wp), allocatable :: wall_distance(:)
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
  if (.not. allocated(xyz_wall)) then
    call collect_wall_fp_coordinates(xyz_wall)
  end if
  !
  ! Create an array containing the coordinates
  ! of all SEM face points on this processor.
  !
  call collect_sem_fp_coordinates(xyz_sem)
  !
  ! Now, do a search for the minimum wall distance of each SEM face point
  !
  call brute_force_wall_distance(xyz_wall,xyz_sem,wall_distance)
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
      my_sem_pts(sem_idx)%sem_fp(n)%wall_distance = wall_distance(i)
    end do
    !
  end do
  !
  ! Deallocate wall_distance since it is no longer needed
  !
  deallocate ( wall_distance , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"wall_distance",2,__LINE__,__FILE__,ierr,error_message)
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
  use geovar,  only : bc_grps
  use flowvar, only : facexyz_array
  !
  !.. Formal Arguments ..
  real(wp), allocatable, intent(out) :: xyz_wall(:,:)
  !
  !.. Local Scalars ..
  integer :: i,ierr,n,np,n1,n2,grp_idx
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
  allocate ( xyz_wall(1:3,1:np) , source=zero , &
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
        i = i + 1
        xyz_wall(1:3,i) = facexyz_array(1:3,n)
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
  call debug_timer(leaving_procedure,pname)
  !
end subroutine collect_wall_fp_coordinates
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
subroutine brute_force_wall_distance(xyz_wall,xyz_sem,wall_distance)
  !
  !.. Formal Arguments ..
  real(wp),               intent(in) :: xyz_wall(:,:)
  real(wp),               intent(in) :: xyz_sem(:,:)
  real(wp), allocatable, intent(out) :: wall_distance(:)
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
  nfp = size(xyz_sem,dim=2)
  !
  allocate ( wall_distance(1:nfp) , source=huge(zero) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"wall_distance",1,__LINE__,__FILE__,ierr,error_message)
  !
  !$OMP PARALLEL DO &
  !$OMP&         DEFAULT(NONE) &
  !$OMP&         SHARED(nfp,nwp,xyz_sem,xyz_wall,wall_distance) &
  !$OMP&         PRIVATE(n,i,xf,yf,zf,dx,dy,dz,dw2,dw2_min)
  !
  do n = 1,nfp
    !
    xf = xyz_sem(1,n)
    yf = xyz_sem(2,n)
    zf = xyz_sem(3,n)
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
      dw2_min = min( dw2_min , dw2 )
      !
    end do
    !
    wall_distance(n) = sqrt( dw2_min )
    !
  end do
  !
  !$OMP END PARALLEL DO
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine brute_force_wall_distance
!
!###############################################################################
!
pure &
subroutine interpolate_profile_data_to_fp(this,prof)
  !
  !.. Formal Arguments ..
  type(sem_fp_t),   intent(inout) :: this
  type(sem_profile_t), intent(in) :: prof
  !
  !.. Local Scalars ..
  integer  :: n,n1,n2
  real(wp) :: c1,c2,y,y1,y2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "interpolate_profile_data_to_fp"
  !
continue
  !
  y = this%wall_distance
  !
  ! Find the indices of the profile data that
  ! bound the wall distance of this face point.
  !
  n1 = 1
  n2 = 2
  do n = 2,prof%npts
    if (prof%y(n) > y) exit
    n1 = n
    n2 = n+1
  end do
  !
  ! Compute the profile data at this face point depending on whether or not
  ! the index of the upper bound is greater than the sizes of the profile data.
  !
  if (n2 > prof%npts) then
    !
    ! The index of the upper bound IS GREATER THAN the size of
    ! the profile data, so just set the profile data for the
    ! face point equal to the profile at the last data point.
    !
    n = prof%npts
    !
    this%ravg = prof%ravg(n)
    this%uavg = prof%uavg(n)
    this%vavg = prof%vavg(n)
    this%wavg = prof%wavg(n)
    this%tavg = prof%tavg(n)
    !
    this%r11 = prof%tau11(n)
    this%r21 = prof%tau21(n)
    this%r22 = prof%tau22(n)
    this%r31 = prof%tau31(n)
    this%r32 = prof%tau32(n)
    this%r33 = prof%tau33(n)
    !
    this%a11 = prof%a11(n)
    this%a21 = prof%a21(n)
    this%a22 = prof%a22(n)
    this%a31 = prof%a31(n)
    this%a32 = prof%a32(n)
    this%a33 = prof%a33(n)
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
    this%ravg = c2*prof%ravg(n2) - c1*prof%ravg(n1)
    this%uavg = c2*prof%uavg(n2) - c1*prof%uavg(n1)
    this%vavg = c2*prof%vavg(n2) - c1*prof%vavg(n1)
    this%wavg = c2*prof%wavg(n2) - c1*prof%wavg(n1)
    this%tavg = c2*prof%tavg(n2) - c1*prof%tavg(n1)
    !
    this%r11 = c2*prof%tau11(n2) - c1*prof%tau11(n1)
    this%r21 = c2*prof%tau21(n2) - c1*prof%tau21(n1)
    this%r22 = c2*prof%tau22(n2) - c1*prof%tau22(n1)
    this%r31 = c2*prof%tau31(n2) - c1*prof%tau31(n1)
    this%r32 = c2*prof%tau32(n2) - c1*prof%tau32(n1)
    this%r33 = c2*prof%tau33(n2) - c1*prof%tau33(n1)
    !
    this%a11 = c2*prof%a11(n2) - c1*prof%a11(n1)
    this%a21 = c2*prof%a21(n2) - c1*prof%a21(n1)
    this%a22 = c2*prof%a22(n2) - c1*prof%a22(n1)
    this%a31 = c2*prof%a31(n2) - c1*prof%a31(n1)
    this%a32 = c2*prof%a32(n2) - c1*prof%a32(n1)
    this%a33 = c2*prof%a33(n2) - c1*prof%a33(n1)
    !
  end if
  !
end subroutine interpolate_profile_data_to_fp
!
!###############################################################################
!
subroutine get_size_of_sem_boxes()
  !
  !.. Local Scalars ..
  integer  :: ierr,n,num_eddies,sem_idx
  real(wp) :: lturb,t_lturb,box_volume
  real(wp) :: xmin,xmax,x
  real(wp) :: ymin,ymax,y
  real(wp) :: zmin,zmax,z
  real(wp) :: rmin,rmax,r
  real(wp) :: tmin,tmax,t
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
      r = hypot(y,z)
      if (r == zero) then
        t = zero
      else
        t = atan2( z , y )
        if (t < zero) t = t + two*pi
      end if
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
      box_volume = half*(xmax-xmin)*(tmax-tmin)*((rmax-rmin)**2)
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
      num_eddies = min( max_eddies , int( box_volume/(lturb**3) ) )
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
subroutine check_restart_eddies()
  !
  !.. Local Scalars ..
  integer  :: ierr,n,sem_idx,eddy_idx
  real(wp) :: vol_coef,eddy_vol,elx,ely,elz
  !
  logical(lk) :: error_found
  !
  !.. Local Arrays ..
  real(wp) :: eddy_loc(1:3)
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
    ! Get the extents of the current SEM box
    !
    box_min(1,sem_idx) = sem_box(sem_idx)%boxmin(1) - eps6
    box_max(1,sem_idx) = sem_box(sem_idx)%boxmax(1) + eps6
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
    is_inside = init_array(num_sem_boxes,num_sem_boxes,fals)
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
  vol_coef = merge( (four/three)*pi , one , eddy_volume_is_spherical )
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
  13 format (5x,a10,i2," [  ",l1,100(3x,l1," ]",:,tl2))
  14 format ("ERROR!!! An error occurred while trying to create a 1-to-1 ", &
             "matching between the groups of eddies read from the restart ", &
             "file and the current set of SEM bounding boxes!")
  !
end subroutine check_restart_eddies
!
!###############################################################################
!
subroutine read_sem_profiles()
  !
  !.. Use Statements
  use ovar, only : bc_in
  !
  !.. Local Scalars ..
  integer  :: i,ierr,io,n,npts,nvar,sem_idx
  real(wp) :: tau11,tau21,tau22,tau31,tau32,tau33
  real(wp) :: a11,a21,a22,a31,a32,a33
  real(wp) :: ravg,vavg,wavg,tavg
  !
  logical(lk) :: profile_data_is_allocated
  logical(lk) :: density_was_read
  logical(lk) :: vy_was_read
  logical(lk) :: vz_was_read
  logical(lk) :: temperature_was_read
  !
  character(len=150) :: fname
  character(len=150) :: profile_variables
  !
  !.. Local Allocatable Arrays ..
  real(wp),          allocatable :: temp_data(:,:)
  character(len=10), allocatable :: var_names(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_sem_profiles"
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
  density_was_read = fals
  vy_was_read = fals
  vz_was_read = fals
  temperature_was_read = fals
  !
  do sem_idx = 1,num_sem_boxes
    !
    fname = sem_box(sem_idx)%profile%fname
    !
    ! Open the profile file for this SEM box
    !
    open (newunit=io,file=trim(adjustl(fname)),status="old",action="read", &
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
    allocate ( sem_box(sem_idx)%profile%y(1:npts) , source=zero )
    !
    allocate ( sem_box(sem_idx)%profile%ravg(1:npts) , source=zero )
    allocate ( sem_box(sem_idx)%profile%uavg(1:npts) , source=zero )
    allocate ( sem_box(sem_idx)%profile%vavg(1:npts) , source=zero )
    allocate ( sem_box(sem_idx)%profile%wavg(1:npts) , source=zero )
    allocate ( sem_box(sem_idx)%profile%tavg(1:npts) , source=zero )
    !
    allocate ( sem_box(sem_idx)%profile%tau11(1:npts) , source=zero )
    allocate ( sem_box(sem_idx)%profile%tau21(1:npts) , source=zero )
    allocate ( sem_box(sem_idx)%profile%tau22(1:npts) , source=zero )
    allocate ( sem_box(sem_idx)%profile%tau31(1:npts) , source=zero )
    allocate ( sem_box(sem_idx)%profile%tau32(1:npts) , source=zero )
    allocate ( sem_box(sem_idx)%profile%tau33(1:npts) , source=zero )
    !
    allocate ( sem_box(sem_idx)%profile%a11(1:npts) , source=zero )
    allocate ( sem_box(sem_idx)%profile%a21(1:npts) , source=zero )
    allocate ( sem_box(sem_idx)%profile%a22(1:npts) , source=zero )
    allocate ( sem_box(sem_idx)%profile%a31(1:npts) , source=zero )
    allocate ( sem_box(sem_idx)%profile%a32(1:npts) , source=zero )
    allocate ( sem_box(sem_idx)%profile%a33(1:npts) , source=zero )
    !
    sem_box(sem_idx)%profile%npts = npts
    !
    ! Read in all the profile data
    !
    allocate ( temp_data(1:nvar,1:npts) , source=zero )
    !
    do n = 1,npts
      read (io,*) (temp_data(i,n), i=1,nvar)
     !read (io,*) sem_box(sem_idx)%profile%y(n), &
     !            sem_box(sem_idx)%profile%uavg(n), &
     !            sem_box(sem_idx)%profile%tau11(n), &
     !            sem_box(sem_idx)%profile%tau22(n), &
     !            sem_box(sem_idx)%profile%tau33(n), &
     !            sem_box(sem_idx)%profile%tau21(n), &
     !            sem_box(sem_idx)%profile%ravg(n), &
     !            sem_box(sem_idx)%profile%tavg(n)
    end do
    !
    ! Close the profile file for this SEM box
    !
    close (io,iostat=ierr,iomsg=error_message)
    call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
    !
    do i = 1,nvar
      select case (uppercase(trim(adjustl(var_names(i)))))
        !
        case ("Y")
          sem_box(sem_idx)%profile%y(1:npts) = temp_data(i,1:npts)
        !
        case ("R")
          sem_box(sem_idx)%is_cylindrical = true
          sem_box(sem_idx)%profile%y(1:npts) = temp_data(i,1:npts)
          !
        case ("RHO","DENSITY")
          density_was_read = true
          sem_box(sem_idx)%profile%ravg(1:npts) = temp_data(i,1:npts)
        case ("VX","U")
          sem_box(sem_idx)%profile%uavg(1:npts) = temp_data(i,1:npts)
        case ("VY","V")
          vy_was_read = true
          sem_box(sem_idx)%profile%vavg(1:npts) = temp_data(i,1:npts)
        case ("VZ","W")
          vz_was_read = true
          sem_box(sem_idx)%profile%wavg(1:npts) = temp_data(i,1:npts)
        case ("T","TEMPERATURE")
          temperature_was_read = true
          sem_box(sem_idx)%profile%tavg(1:npts) = temp_data(i,1:npts)
          !
        case ("TAUXX","TAU11","RXX","R11")
          sem_box(sem_idx)%profile%tau11(1:npts) = temp_data(i,1:npts)
        case ("TAUYX","TAUXY","TAU21","TAU12","RYX","RXY","R21","R12")
          sem_box(sem_idx)%profile%tau21(1:npts) = temp_data(i,1:npts)
        case ("TAUYY","TAU22","RYY","R22")
          sem_box(sem_idx)%profile%tau22(1:npts) = temp_data(i,1:npts)
        case ("TAUZX","TAUXZ","TAU31","TAU13","RZX","RXZ","R31","R13")
          sem_box(sem_idx)%profile%tau31(1:npts) = temp_data(i,1:npts)
        case ("TAUZY","TAUYZ","TAU32","TAU23","RZY","RYZ","R32","R23")
          sem_box(sem_idx)%profile%tau32(1:npts) = temp_data(i,1:npts)
        case ("TAUZZ","TAU33","RZZ","R33")
          sem_box(sem_idx)%profile%tau33(1:npts) = temp_data(i,1:npts)
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
    call nondimensionalize_sem_profile(sem_box(sem_idx)%profile)
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
    if (.not. density_was_read) then
      ravg = bc_in( sem_box(sem_idx)%bc_idx )%cv(nec)
      do n = 1,npts
        sem_box(sem_idx)%profile%ravg(n) = ravg
      end do
    end if
    !
    ! If the averaged y-velocity was not included in the profile
    ! file, set the y-velocity for all points in the profile to
    ! the y-velocity specified for this boundary condition.
    !
    if (.not. vy_was_read) then
      vavg = bc_in( sem_box(sem_idx)%bc_idx )%pv(nmy)
      do n = 1,npts
        sem_box(sem_idx)%profile%vavg(n) = vavg
      end do
    end if
    !
    ! If the averaged z-velocity was not included in the profile
    ! file, set the z-velocity for all points in the profile to
    ! the z-velocity specified for this boundary condition.
    !
    if (.not. vz_was_read) then
      wavg = bc_in( sem_box(sem_idx)%bc_idx )%pv(nmz)
      do n = 1,npts
        sem_box(sem_idx)%profile%wavg(n) = wavg
      end do
    end if
    !
    ! If the averaged temperature was not included in the profile
    ! file, set the temperature for all points in the profile to
    ! the temperature specified for this boundary condition.
    !
    if (.not. temperature_was_read) then
      tavg = bc_in( sem_box(sem_idx)%bc_idx )%gp(ntemp)
      do n = 1,npts
        sem_box(sem_idx)%profile%tavg(n) = tavg
      end do
    end if
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_sem_profiles
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
  integer  :: n
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
  do n = 1,this%npts
    !
    this%ravg(n) = this%ravg(n) * rhoref_inv
    this%uavg(n) = this%uavg(n) * aref_inv
    this%vavg(n) = this%vavg(n) * aref_inv
    this%wavg(n) = this%wavg(n) * aref_inv
    this%tavg(n) = this%tavg(n) * tref_inv
    !
    this%tau11(n) = this%tau11(n) * aref2_inv
    this%tau21(n) = this%tau21(n) * aref2_inv
    this%tau22(n) = this%tau22(n) * aref2_inv
    this%tau31(n) = this%tau31(n) * aref2_inv
    this%tau32(n) = this%tau32(n) * aref2_inv
    this%tau33(n) = this%tau33(n) * aref2_inv
    !
  end do
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
  !
  !.. Local Scalars ..
  integer :: sem_idx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_details_of_sem_boxes"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (mypnum == glb_root) then
    do sem_idx = 1,num_sem_boxes
      write (iout,*)
      write (iout,1) sem_idx
      write (iout,*) 'Uconv     = ',sem_box(sem_idx)%convective_velocity*aref
      write (iout,*) 'Lturb     = ',sem_box(sem_idx)%length_scale
      write (iout,*) 'Neddies   = ',sem_box(sem_idx)%num_eddies
      write (iout,*) 'Volbox    = ',sem_box(sem_idx)%box_volume
      write (iout,*) 'Min x     = ',sem_box(sem_idx)%boxmin(1)
      write (iout,*) 'Max x     = ',sem_box(sem_idx)%boxmax(1)
      write (iout,*) 'Min y     = ',sem_box(sem_idx)%boxmin(2)
      write (iout,*) 'Max y     = ',sem_box(sem_idx)%boxmax(2)
      write (iout,*) 'Min z     = ',sem_box(sem_idx)%boxmin(3)
      write (iout,*) 'Max z     = ',sem_box(sem_idx)%boxmax(3)
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
  integer  :: io,ierr,n,npts,bc_idx
  real(wp) :: ymax
  !
  character(len=100) :: fname
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: y(:)
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
        ! Skip reading the variable names
        ! NOTE: This assumes that y/wall-distance is the first variable
        !       in the profile file.
        !
        read (io,*)
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
          read (io,*) y(n)
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
  integer  :: i,ierr,n,nc,nf,nn,nbc,nbcn,npts,bc_idx
  integer  :: n_sem_pts,n_sem_bc,n_wall_pts
  real(wp) :: profile_height
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: xyz_sem(:,:)
  real(wp), allocatable :: wall_distance(:)
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
    n_wall_pts = count( msk(1:nnode) )
    !
    ! F2003 AUTO-REALLOC
    xyz_wall = xyz_nodes( 1:3 , pack(intseq(1,nnode),mask=msk(1:nnode)) )
    !
  end if
  !
  ! Broadcast the number of wall boundary nodes
  ! to all the non-host-root MPI processes
  !
  if (ncpu > 1) then
    call mpi_bcast(n_wall_pts,1_int_mpi,mpi_inttyp, &
                   host_root,my_host_comm,mpierr)
  end if
  !
  if (.not. allocated(xyz_wall)) then
    allocate ( xyz_wall(1:3,1:n_wall_pts) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_wall",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Broadcast the coordinates of all the grid nodes on wall BCs.
  ! NOTE: All MPI processes need a copy of xyz_wall because it will be used
  !       later in compute_sem_fp_wall_distance to compute the distance to
  !       the nearest wall for all boundary face points on SEM boundaries.
  !
  if (ncpu > 1) then
    call mpi_bcast(xyz_wall,size(xyz_wall,kind=int_mpi),mpi_flttyp, &
                   host_root,my_host_comm,mpierr)
  end if
  !
  ! Find all the boundary nodes that are on SEM boundaries
  !
  if (i_am_host_root) then
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
    call brute_force_wall_distance(xyz_wall,xyz_sem,wall_distance)
    !
  end if
  !
  ! Deallocate msk, sem2node, and xyz_sem since they are no longer needed
  !
  if (allocated(msk)) then
    deallocate ( msk , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
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
  ! Create the separate subsonic inflow boundary condition for each SEM boundary
  !
  do bc_idx = 1,size(bc_input)
    if (len_trim(adjustl(bc_input(bc_idx)%sem%profile_file)) > 0) then
      !
      nbc = nbc + 1
      nbcn = nbcn + 1
      !
      tmp_bc_input(nbc) = bc_input(bc_idx)
      tmp_bc_input(nbc)%set_bc_value = bc_sub_inflow
      tmp_bc_input(nbc)%name = trim(adjustl(bc_input(bc_idx)%name)) // &
                               " (Subsonic Inflow)"
      !
      tmp_bc_input(bc_idx)%sem%sub_in_bc_idx = nbc
      tmp_bc_input(bc_idx)%sem%sub_in_bc_grp = nbcn
      !
      tmp_bc_names(nbcn) = trim(adjustl(tmp_bc_input(nbc)%name))
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
          ! i.e., the index of this node within the array wall_distance
          !
          nn = node2sem( grid%elem(nc)%nodes(n) )
          !
          ! If the distance to the nearest wall of the current node is
          ! within the profile height of the SEM boundary profile, exit
          ! from the one-pass-loop so that the boundary conditions for
          ! this boundary face are unchanged and remain an SEM BC.
          !
          if (wall_distance(nn) <= profile_height) exit one_pass_loop
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
  ! Deallocate local arrays before leaving if still allocated
  !
  if (allocated(wall_distance)) then
    deallocate ( wall_distance , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"wall_distance",2,__LINE__,__FILE__,ierr,error_message)
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
subroutine output_sem_averages()
  !
  !.. Local Scalars ..
  integer :: i,ierr,n,nvar,num_fp,num_total_sem_pts,sem_idx
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
  ! Have all MPI processes fill collect all their SEM averaged values into
  ! a local array.
  !
  nvar = 42
  !
  my_ave_data = init_array(nvar,num_my_sem_pts(mypnum+1),zero) ! F03 AUTO-REALOC
  !
  call copy_my_sem_ave_data(my_ave_data)
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
    ! NOTE: Only the global root process will allocate this array to the actually
    !       needed size; all other processes will allocate it such that there is
    !       only a single SEM face point.
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
  ! Now have the global root process group all the SEM face points by
  ! their respective SEM boxes so that the averaged data can be output.
  !
  if (mypnum == glb_root) then
    !
    call write_sem_ave_data(num_sem_pts,all_ave_data)
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
  nvar = 17
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
    write (io,2) sem_idx,"Averaged",num_fp
    !
    ! Write the SEM averaged data for this SEM box
    !
    do m = 1,nvar
      write (io,3) (this_ave(m,i), i=1,num_fp)
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
    write (io,2) sem_idx,"Profile",num_fp
    !
    ! Write the SEM averaged data for this SEM box
    !
    do m = 1,nvar
      write (io,3) (this_ave(m,i), i=1,num_fp)
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
    write (io,2) sem_idx,"Averaged-Profile",num_fp
    !
    ! Write the SEM averaged data for this SEM box
    !
    do m = 1,nvar
      write (io,3) (this_ave(m,i), i=1,num_fp)
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
            '"u<math>\\034</math>u<math>\\034</math>" ',/, &
            '"v<math>\\034</math>v<math>\\034</math>" ',/, &
            '"w<math>\\034</math>w<math>\\034</math>" ',/, &
            '"u<math>\\034</math>v<math>\\034</math>" ',/, &
            '"u<math>\\034</math>w<math>\\034</math>" ',/, &
            '"v<math>\\034</math>w<math>\\034</math>" ',/, &
            '"u<math>\\034</math>" ',/, &
            '"v<math>\\034</math>" ',/, &
            '"w<math>\\034</math>" ')
  2 format ('ZONE T="SEM BC #',i0,' (',a,')", I=',i0,', DATAPACKING=BLOCK')
  3 format (6es14.6)
  !
end subroutine write_sem_ave_data
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
  integer :: i,n,sem_idx
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
    end do
  end do
  !
end subroutine copy_my_sem_ave_data
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
  end do
  !
end subroutine dimensionalize_sem_ave
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
