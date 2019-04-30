module synthetic_eddy_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nee
  !
  implicit none
  !
  private
  !
  public :: sem
  public :: separate_sem_boundary_faces
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
   !! SKIP_THIS_FP : Logical flag to indicate that this SEM face point is
   !!                outside the SEM inflow profile/boundary-layer, so it
   !!                should be skipped when computing the SEM inflow
   !!                conditions. The inflow conditions for these skipped
   !!                points should be computed using a regular inflow BC
   !!                once the SEM inflow conditions for the not-skipped
   !!                points is completed.
   !logical(lk) :: skip_this_fp = fals
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
    ! BC_GRP : Boundary condition group these SEM points correspond to
    integer :: bc_grp
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
    !          NOTE: The eddy volume is stored as [1/sqrt(eddy volume)] because
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
    ! BC_GRP : BC group that corresponds to this SEM box
    integer :: bc_grp = 0
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
    ! PROFILE : inflow profile related to this SEM box
    type(sem_profile_t) :: profile
    ! EDDY : information about each eddy in this SEM box
    type(eddy_t), allocatable :: eddy(:)
    !
  end type sem_box_t
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! SEM PARAMETERS
  !
  integer, parameter :: max_eddies = 5000
  !
  logical(lk), parameter :: sem_box_is_cylindrical = fals
 !logical(lk), parameter :: sem_box_is_cylindrical = true
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
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: cv(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,ierr,sem_idx,eddy_idx,num_eddies
  real(wp) :: ave,running_ave,time_scale
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
  logical(lk), parameter :: use_improved_SEM_normalization = fals
 !logical(lk), parameter :: use_improved_SEM_normalization = true
  !
  real(wp), parameter :: sqrt32 = sqrt(three*half)
  !
continue
  !
  if (.not. sem_is_initialized) then
    call initialize_sem()
    sem_is_initialized = true
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
  loop_over_my_sem_boxes: &
  do sem_idx = 1,size(my_sem_pts)
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
   !c0 = sqrt(boxvol*half*three)*half*three/(lturb**6)/real(num_eddies,kind=wp)
    !
    ! Pre-compute the part from the normalization constant
    !
    c0 = sqrt(boxvol)/real(num_eddies,kind=wp)
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
        c1 = d_t/time_scale
        !
        running_ave = c1*ave + (one-c1)*running_ave
        !
        ufluc2 = ufluc / sqrt(running_ave)
        vfluc2 = vfluc / sqrt(running_ave)
        wfluc2 = wfluc / sqrt(running_ave)
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
    end do loop_over_sem_pts
    !
  end do loop_over_my_sem_boxes
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
  allocate ( new_eddies(1:size(sem_box)) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"new_eddies",2,__LINE__,__FILE__,ierr,error_message)
  !
  loop_over_all_sem_boxes: &
  do sem_idx = 1,size(sem_box)
    !
    ! Get the axial location of the outflow plane of the current SEM box
    !
    boxmaxx = sem_box(sem_idx)%boxmax(1)
    !
    ! Get the axial distance traveled by the eddies over the past time step
    !
    dx = sem_box(sem_idx)%convective_velocity * d_t
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
  do sem_idx = 1,size(new_eddies)
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
end subroutine sem
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
  ! TREAT EDDY VOLUMES AS SPHERE-LIKE
 !logical(lk), parameter :: eddy_volume_is_spherical = true
  ! TREAT EDDY VOLUMES AS BOX-LIKE
  logical(lk), parameter :: eddy_volume_is_spherical = fals
  !
continue
  !
  x_coef = merge( zero , one , generate_at_xmin )
  !
  vol_coef = merge( (four/three)*pi , one , eddy_volume_is_spherical )
  !
  loop_over_sem_boxes: &
  do sem_idx = 1,size(new_eddies)
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
    if (sem_box_is_cylindrical) then
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
      if (sem_box_is_cylindrical) then
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
      ! Store the eddy volume as the reciprocal
      ! of the squareroot of the volume
      !
      sem_box(sem_idx)%eddy(eddy_idx)%volume = one / sqrt(eddy_vol)
      !
    end do loop_over_new_eddies_in_box
    !
  end do loop_over_sem_boxes
    !
end subroutine generate_new_eddies
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
  loop_to_add_new_eddies: &
  do sem_idx = 1,size(new_eddies)
    !
    lturb = sem_box(sem_idx)%length_scale
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
    this%fp_list = spread(0,dim=1,ncopies=20) ! F2003 AUTO-REALLOC
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
    this%fp_list = [ this%fp_list , spread(0,dim=1,ncopies=20) ]
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
  this%fp_list = spread(0,dim=1,ncopies=20) ! F2003 AUTO-REALLOC
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
    this%eddy_list = spread(0,dim=1,ncopies=20) ! F2003 AUTO-REALLOC
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
    this%eddy_list = [ this%eddy_list , spread(0,dim=1,ncopies=20) ]
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
    this%eddy_list = spread(0,dim=1,ncopies=20) ! F2003 AUTO-REALLOC
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
    this%eddy_list = [ this%eddy_list , spread(0,dim=1,ncopies=20) ]
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
  use ovar,       only : bc_in
  use geovar,     only : bc_grps
  use flowvar,    only : facexyz_array
  use random_mod, only : set_random_seed
  !
  !.. Local Scalars ..
  integer  :: ierr,k,n,np,n1,n2
  integer  :: sem_idx,grp_idx,bc_idx,num_eddies
  integer  :: num_restart_sem_boxes
  real(wp) :: lturb,uconv
  !
  logical(lk) :: sem_from_restart
  !
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
  ! Initialize the random number generator
  !
  call set_random_seed
  !
  ! Count the number of BC groups on this processor that are of the SEM type
  !
  num_sem_boxes = count( bc_grps(:)%bc_type == bc_sem )
  !
  ! Allocate the derived-type sem_box if needed.
  !
  if (allocated(sem_box)) then
    !
    sem_from_restart = true
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
    sem_from_restart = fals
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
  sem_idx = 0
  do grp_idx = 1,size(bc_grps)
    !
    if (bc_grps(grp_idx)%bc_type == bc_sem) then
      !
      bc_idx = bc_grps(grp_idx)%bc_idx
      !
      sem_idx = sem_idx + 1
      !
      lturb = bc_in(bc_idx)%sem%length_scale
      uconv = bc_in(bc_idx)%sem%convective_velocity
      fname = trim(adjustl(bc_in(bc_idx)%sem%profile_file))
      !
      if (.not. sem_from_restart) then
        sem_box(sem_idx)%bc_grp = grp_idx
        sem_box(sem_idx)%length_scale = lturb
        sem_box(sem_idx)%convective_velocity = uconv
        sem_box(sem_idx)%profile%fname = trim(adjustl(fname))
      end if
      !
      n1 = bc_grps(grp_idx)%beg_fp
      n2 = bc_grps(grp_idx)%end_fp
      np = n2-n1+1
      !
      my_sem_pts(sem_idx)%bc_grp = grp_idx ! Index for this BC group
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
        k = my_sem_pts(sem_idx)%sem_fp(n)%idx
        my_sem_pts(sem_idx)%sem_fp(n)%x = facexyz_array(1,k)
        my_sem_pts(sem_idx)%sem_fp(n)%y = facexyz_array(2,k)
        my_sem_pts(sem_idx)%sem_fp(n)%z = facexyz_array(3,k)
      end do
      !
    end if
    !
  end do
  !
  ! For each SEM box:
  !   1) get the min/max extents of the box,
  !   2) compute the volume of the box, and
  !   3) compute the number of eddies to use for the box.
  !
  call get_size_of_sem_boxes(sem_from_restart)
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
  do sem_idx = 1,size(new_eddies)
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
  ! Now update the face point eddy-lists and the 
  ! eddy face-point-lists using new_eddies
  !
  call update_eddy_lists(new_eddies)
  !
  ! Deallocate new_eddies since we no longer need it
  !
  deallocate ( new_eddies , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"new_eddies",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Compute the distance to the nearest wall for each of the SEM face points
  !
  call compute_sem_fp_wall_distance()
  !
  ! Finally, interpolate the profile data to the location of each face point
  !
  do sem_idx = 1,size(sem_box)
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
  do sem_idx = 1,size(my_sem_pts)
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
  call mpi_allgather(MPI_IN_PLACE,0_int_mpi,MPI_DATATYPE_NULL, &
                     cpu_nums,1_int_mpi,MPI_INTEGER, &
                     MPI_COMM_WORLD,mpierr)
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
  call mpi_allgatherv(MPI_IN_PLACE,0_int_mpi,MPI_DATATYPE_NULL, &
                      xyz_wall,cpu_nums,cpu_disp,mpi_flttyp, &
                      MPI_COMM_WORLD,mpierr)
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
  np = sum(my_sem_pts(:)%num_fp)
  !
  allocate ( xyz_sem(1:3,1:np) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_sem",1,__LINE__,__FILE__,ierr,error_message)
  !
  i = 0
  do sem_idx = 1,size(my_sem_pts)
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
subroutine get_size_of_sem_boxes(sem_from_restart)
  !
  !.. Formal Arguments ..
  logical(lk), intent(in) :: sem_from_restart
  !
  !.. Local Scalars ..
  integer  :: ierr,n,num_eddies,sem_idx
  real(wp) :: lturb,box_volume
  real(wp) :: xmin,xmax,x
  real(wp) :: ymin,ymax,y
  real(wp) :: zmin,zmax,z
  real(wp) :: rmin,rmax,r1,r2,r3,r4
  real(wp) :: tmin,tmax,t1,t2,t3,t4
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: boxmin(:,:)
  real(wp), allocatable :: boxmax(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_size_of_sem_boxes"
  !
continue
  !
  ! Allocate boxmin and boxmax to store the min/max extents
  ! (local to this partition/process) of the SEM boxes
  !
  allocate ( boxmin(1:5,1:num_sem_boxes) , source=huge(zero) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"boxmin",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( boxmax(1:5,1:num_sem_boxes) , source=-huge(zero) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"boxmax",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Compute the min/max extents (local to this
  ! partition/process) of the SEM boxes
  !
  do sem_idx = 1,size(my_sem_pts)
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
      !
      xmin = x - lturb
      xmax = x + lturb
      ymin = y - lturb
      ymax = y + lturb
      zmin = z - lturb
      zmax = z + lturb
      !
      boxmin(1,sem_idx) = min(xmin,boxmin(1,sem_idx))
      boxmin(2,sem_idx) = min(ymin,boxmin(2,sem_idx))
      boxmin(3,sem_idx) = min(zmin,boxmin(3,sem_idx))
      !
      boxmax(1,sem_idx) = max(xmax,boxmax(1,sem_idx))
      boxmax(2,sem_idx) = max(ymax,boxmax(2,sem_idx))
      boxmax(3,sem_idx) = max(zmax,boxmax(3,sem_idx))
      !
    end do
    !
    ! Compute the r,theta min/max for the current box from the y,z min/max
    !
    ymin = boxmin(2,sem_idx)
    zmin = boxmin(3,sem_idx)
    !
    ymax = boxmax(2,sem_idx)
    zmax = boxmax(3,sem_idx)
    !
    r1 = sqrt( ymin*ymin + zmin*zmin )
    r2 = sqrt( ymax*ymax + zmin*zmin )
    r3 = sqrt( ymin*ymin + zmax*zmax )
    r4 = sqrt( ymax*ymax + zmax*zmax )
    !
    rmin = max( zero , min(r1,r2,r3,r3) )
    rmax = max( zero , max(r1,r2,r3,r4) )
    !
    t1 = zero
    t2 = zero
    t3 = zero
    t4 = zero
    !
    if (ymin /= zero .or. zmin /= zero) t1 = atan2( zmin , ymin )
    if (ymax /= zero .or. zmin /= zero) t2 = atan2( zmin , ymax )
    if (ymin /= zero .or. zmax /= zero) t3 = atan2( zmax , ymin )
    if (ymax /= zero .or. zmax /= zero) t4 = atan2( zmax , ymax )
    !
    ! The atan2 function returns a value on the interval [-Pi,Pi].
    ! Converting an angle from this interval to the interval [0,2*Pi]
    ! simply requires adding 2*Pi to any angle value less than zero.
    ! For example:
    !    theta[-Pi,Pi] = -Pi    => theta[0,2*Pi] = Pi     = 180 degrees
    !    theta[-Pi,Pi] = -Pi/2  => theta[0,2*Pi] = 3*Pi/2 = 270 degrees
    !    theta[-Pi,Pi] = -0     => theta[0,2*Pi] = 2*Pi   = 360 degrees
    !
    if (t1 < zero) t1 = t1 + two*pi
    if (t2 < zero) t2 = t2 + two*pi
    if (t3 < zero) t3 = t3 + two*pi
    if (t4 < zero) t4 = t4 + two*pi
    !
    tmin = min(t1,t2,t3,t4)
    tmax = max(t1,t2,t3,t4)
    !
    boxmin(4,sem_idx) = rmin
    boxmin(5,sem_idx) = tmin
    !
    boxmax(4,sem_idx) = rmax
    boxmax(5,sem_idx) = tmax
    !
  end do
  !
  ! Reduce the boxmin and boxmax arrays across all processes so that
  ! every process has the true global min/max extents for each SEM box
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,boxmin,size(boxmin,kind=int_mpi), &
                       mpi_flttyp,MPI_MIN,MPI_COMM_WORLD,mpierr)
    call mpi_allreduce(MPI_IN_PLACE,boxmax,size(boxmax,kind=int_mpi), &
                       mpi_flttyp,MPI_MAX,MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Copy the global min/max extents for each SEM box into sem_box.
  ! Also compute the volume and the number of eddies for each box.
  !
  do sem_idx = 1,size(sem_box)
    !
    lturb = sem_box(sem_idx)%length_scale
    !
    xmin = boxmin(1,sem_idx)
    ymin = boxmin(2,sem_idx)
    zmin = boxmin(3,sem_idx)
    rmin = boxmin(4,sem_idx)
    tmin = boxmin(5,sem_idx)
    !
    xmax = boxmax(1,sem_idx)
    ymax = boxmax(2,sem_idx)
    zmax = boxmax(3,sem_idx)
    rmax = boxmax(4,sem_idx)
    tmax = boxmax(5,sem_idx)
    !
    ! Compute the volume of this SEM box
    !
    if (sem_box_is_cylindrical) then
      box_volume = half*(xmax-xmin)*(tmax-tmin)*((rmax-rmin)**2)
    else
      box_volume = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
    end if
    !
    ! If the SEM eddies were not read in from a restart file,
    ! compute the number of eddies to use for this box based on the box
    ! volume and the turbulence length scale specified for this SEM box.
    !
    if (.not. sem_from_restart) then
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
  !
  ! Deallocate boxmin and boxmax before we leave
  !
  deallocate ( boxmin , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"boxmin",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( boxmax , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"boxmax",2,__LINE__,__FILE__,ierr,error_message)
  !
end subroutine get_size_of_sem_boxes
!
!###############################################################################
!
subroutine read_sem_profiles()
  !
  !.. Local Scalars ..
  integer  :: i,ierr,io,n,npts,nvar,sem_idx
  real(wp) :: tau11,tau21,tau22,tau31,tau32,tau33
  real(wp) :: a11,a21,a22,a31,a32,a33
  !
  logical(lk) :: profile_data_is_allocated
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
  ! First, check to see if the profile data has already been allocated.
  ! If it has been allocated, it is assumed that the profile data also
  ! exists in these allocated arrays.
  !
  call check_sem_profile_allocations(profile_data_is_allocated)
  !
  ! If the profile data has already been allocated, we do not need to read
  ! the data in so exit from this subroutine.
  !
  if (profile_data_is_allocated) return
  !
  ! We need to create the profile data, so read the data from the profile files.
  !
  do sem_idx = 1,size(sem_box)
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
        case ("RHO","DENSITY")
          sem_box(sem_idx)%profile%ravg(1:npts) = temp_data(i,1:npts)
        case ("VX","U")
          sem_box(sem_idx)%profile%uavg(1:npts) = temp_data(i,1:npts)
        case ("VY","V")
          sem_box(sem_idx)%profile%vavg(1:npts) = temp_data(i,1:npts)
        case ("VZ","W")
          sem_box(sem_idx)%profile%wavg(1:npts) = temp_data(i,1:npts)
        case ("T","TEMPERATURE")
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
      a22 = sqrt(    tau22 -(a21**2))
     !a22 = sqrt(abs(tau22)-(a21**2))
      !
      a31 = tau31/(a11+rndoff)
      a32 = (tau32-a21*a31)/(a22+rndoff)
      a33 = sqrt(    tau33 -(a31**2)-(a32**2))
     !a33 = sqrt(abs(tau33)-(a31**2)-(a32**2))
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
  end do
  !
end subroutine read_sem_profiles
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
 !m = 0
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
            nn = len(clist) - len(strng(ibeg:is))
            ! F2003 AUTO-REALLOC
            clist = [ clist , strng(ibeg:is) // repeat(" ",ncopies=nn) ]
           !call reallocate( clist , m+1 )
           !clist(m+1) = strng(ibeg:is)
          else if (strng(is:is) == " ") then
            nn = len(clist) - len(strng(ibeg:is-1))
            ! F2003 AUTO-REALLOC
            clist = [ clist , strng(ibeg:is-1) // repeat(" ",ncopies=nn) ]
           !m = m + 1
           !call reallocate( clist , m )
           !clist(m) = strng(ibeg:is-1)
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
  profile_data_is_allocated = fals
  !
  num_allocs = 0
  num_checks = 0
  !
  ! Check the allocation status of all the profile arrays
  ! for each of the SEM boxes.
  !
  do sem_idx = 1,size(sem_box)
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
  do sem_idx = 1,size(sem_box)
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
        y = spread(zero,dim=1,ncopies=npts) ! F2003 AUTO-REALLOC
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
    y = spread(zero,dim=1,ncopies=npts) ! F2003 AUTO-REALLOC
    !
  end if
  !
  call mpi_bcast(y,size(y,kind=int_mpi),mpi_flttyp, &
                 glb_root,MPI_COMM_WORLD,mpierr)
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
  use geovar, only : bface,grid
  use ovar,   only : bc_input,bc_input_t
  !
  !.. Local Scalars ..
  integer  :: i,ierr,n,nc,nf,nn,nbc,npts,bc_idx
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
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "separate_sem_boundary_faces"
  !
continue
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
    msk = spread(fals,dim=1,ncopies=nnode) ! F2003 AUTO-REALLOC
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
    xyz_wall = grid%xyz( 1:3 , pack(intseq(1,nnode),mask=msk(1:nnode)) )
    !
  end if
  !
  ! Broadcast the number of wall boundary nodes
  ! to all the non-host-root MPI processes
  !
  call mpi_bcast(n_wall_pts,1_int_mpi,mpi_inttyp, &
                 host_root,my_host_comm,mpierr)
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
  call mpi_bcast(xyz_wall,size(xyz_wall,kind=int_mpi),mpi_flttyp, &
                 host_root,my_host_comm,mpierr)
  !
  ! Find all the boundary nodes that are on SEM boundaries
  !
  if (i_am_host_root) then
    !
    msk = spread(fals,dim=1,ncopies=nnode) ! F2003 AUTO-REALLOC
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
    xyz_sem = grid%xyz( 1:3 , sem2node(1:n_sem_pts) )
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
  do n = 1,size(bc_input)
    if (len_trim(adjustl(bc_input(n)%sem%profile_file)) > 0) then
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
  do n = 1,size(bc_input)
    tmp_bc_input(n) = bc_input(n)
  end do
  !
  ! Create the separate subsonic inflow boundary condition for each SEM boundary
  !
  do bc_idx = 1,size(bc_input)
    if (len_trim(adjustl(bc_input(n)%sem%profile_file)) > 0) then
      !
      nbc = nbc + 1
      !
      tmp_bc_input(nbc) = bc_input(bc_idx)
      tmp_bc_input(nbc)%set_bc_value = bc_sub_inflow
      tmp_bc_input(nbc)%name = trim(adjustl(bc_input(bc_idx)%name)) // &
                               " (Subsonic Inflow)"
      !
      tmp_bc_input(bc_idx)%sem%sub_in_bc_grp = nbc
      !
    end if
  end do
  !
  ! Transfer tmp_bc_input to bc_input
  !
  call move_alloc( from=tmp_bc_input , to=bc_input )
  !
  ! Loop through the boundary faces to find all the boundary faces that are
  ! on an SEM boundary. Check the wall distance of every node for each SEM
  ! boundary face, and change the boundary condition for the boundary face
  ! to a subsonic inflow boundary condition if all the nodes are further from
  ! the nearest wall than the profile height of the SEM boundary profile.
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
      bc_idx = bface(5,nf)
      !
      profile_height = bc_input(bc_idx)%sem%profile_height
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
  call mpi_bcast(bface,size(bface,kind=int_mpi),mpi_inttyp, &
                 host_root,my_host_comm,mpierr)
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
end subroutine separate_sem_boundary_faces
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
