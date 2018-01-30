module initialization_mod
  !
  ! Module containing the functions and routines to
  ! initialize the solution
  !
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,ntk,ntl,nq
  !
  implicit none
  !
  private
  !
  public :: initialize_flow
  public :: initialize_reference
  !
contains
!
!###############################################################################
!
subroutine initialize_flow()
  !
  ! Routine to initialize the flow
  !
  !.. Use Statements ..
  use generic_types_mod, only : matrix
  !
  use geovar, only : nr,n_solpts,ncell,cell,xyz,grid
  !
  use ovar, only : rhoref,aref,itestcase
  use ovar, only : cv_in
  use ovar, only : load_restart_file,d_t,constant_dt,itcur,time
  use ovar, only : dtrst,itrst,tmrst
  use ovar, only : rhoref,pref,tref,gam,rgasref,cpref
  use ovar, only : aref,machref,velref,uref,muref,reyref
  use ovar, only : alpha_aoaref,beta_aoaref
  use ovar, only : pref_nd,velref_nd,muref_nd
  use ovar, only : rgasref_nd,cpref_nd
  use ovar, only : adjust_cfl_cycles
 !use ovar, only : output_time_averaging
  !
  use flowvar, only : usp,uoldsp
  use flowvar, only : allocate_flowvar_arrays
  !
  use order_mod, only : maxpts
  !
  use restart_mod, only : read_restart_file
  use restart_mod, only : create_restart_fname_formats
  !
  use projection_mod, only : project_to_porder
  !
  use interpolation_mod, only : outerp
  !
 !use time_mod, only : initialize_time_ave
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,np,l,l1,l2,m
  integer :: gmin,gmax,omin,omax,ierr
  integer :: this_geom,this_order
  integer :: init_cond
  real(wp) :: rhoa2
  character(len=200) :: array_name
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: xyz_init(:,:)
  type(matrix), allocatable, target :: outerp_inv(:,:)
  !
  !.. Local Pointers ..
 !real(wp), pointer, contiguous :: numer2init(:,:)
 !real(wp), pointer, contiguous :: init2numer(:,:)
  real(wp), pointer :: numer2init(:,:)
  real(wp), pointer :: init2numer(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_flow"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize the bc_in array
  !
  call create_bc_in
  !
  ! Initialize all local pointers to disassociated
  !
  numer2init => null()
  init2numer => null()
  !
  ! Some constants from the free-stream conditions
  !
  rhoa2 = rhoref * aref * aref
  !
  ! Allocate the arrays within the module flowvar
  !
  call allocate_flowvar_arrays
  !
  ! Check for a restart job, otherwise initialize usp
  !
  if (load_restart_file) then
    !
    ! Read in usp and the current iteration, time step,
    ! and elapsed time from the restart file
    !
    call read_restart_file
    !
    ! Initialize the iteration counter, time step,
    ! and elapsed time to that in the restart file
    !
    d_t   = dtrst
    itcur = itrst
    time  = tmrst
    !
  else
    !
    if (itestcase == Nozzle_Jet) then
      call find_elem_for_initial_jet_profile(point_is_in_smc000_nozzle)
    else if (itestcase == -Nozzle_Jet) then
      call find_elem_for_initial_jet_profile(point_is_in_arn2_nozzle)
    end if
    !
    allocate ( xyz_init(1:nr,1:maxpts) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_init",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Allocate the local array that stores the inverse of the matrix used for
    ! interpolating from the solution points to the initialization points.
    ! NOTE: This is array is used so that we only invert the matrix once for a
    !       given cell geometry/order combination.  We dont want to keep
    !       inverting the same matrix over and over for cells with the same
    !       geometry/order combination.
    !
    gmin = lbound(outerp,dim=1)
    gmax = ubound(outerp,dim=1)
    omin = lbound(outerp,dim=2)
    omax = ubound(outerp,dim=2)
    !
    allocate ( outerp_inv(gmin:gmax,omin:omax) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"outerp_inv",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Initialize the solution in each cell
    !
    do n = 1,ncell
      !
      this_geom = cell(n)%geom
      this_order = cell(n)%order
      !
      n1 = cell(n)%beg_sp
      n2 = cell(n)%end_sp
      np = n2-n1+1
      !
      init_cond = 0
      if (allocated(grid)) then
        if (allocated(grid%elem)) then
          if (allocated(grid%elem(n)%tags)) then
            if (size(grid%elem(n)%tags) >= 2) then
              init_cond = grid%elem(n)%tags(2)
            end if
          end if
        end if
      end if
      !
      ! Check whether the init component of outerp is associated with
      ! an interpolation matrix to determine where we need to compute
      ! the initial conditions.
      !
      if (.not. allocated(outerp(this_geom,this_order)%init)) then
        !
        ! The init component of outerp is not associated so just compute
        ! the initial conditions at the solution points
        !
        usp(1:nq,n1:n2) = initial_conditions( xyz(1:nr,n1:n2) , init_cond )
        !
      else
        !
        ! The init component of outerp is associated so we need to compute
        ! the initial conditions at the initialization points before
        ! interpolating it back to the solution points.
        !
        ! Assign local pointer to this matrix as an alias to simplify the code
        !
        numer2init => outerp(this_geom,this_order)%init%mat
        !
        ! Check if we have already created an inverse of this interpolation
        ! matrix. If we havent, allocate its location in outerp_inv and
        ! compute this inverse matrix.
        !
        if (.not. allocated(outerp_inv(this_geom,this_order)%mat)) then
          !
          l1 = size(numer2init,dim=1)
          l2 = size(numer2init,dim=2)
          !
          allocate ( outerp_inv(this_geom,this_order)%mat(1:l1,1:l2) , &
                     source=zero , stat=ierr , errmsg = error_message )
          write (array_name,1) Geom_Name(this_geom),this_order
          call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                           error_message)
          !
          outerp_inv(this_geom,this_order)%mat = invert_matrix(numer2init)
          !
        end if
        !
        ! Assign local pointer to this matrix as an alias to simplify the code
        !
        init2numer => outerp_inv(this_geom,this_order)%mat
        !
        ! Interpolate the cell solution points to the initialization points
        !
        do l = 1,nr
          xyz_init(l,1:np) = matmul( numer2init , xyz(l,n1:n2) )
        end do
        !
        ! Evaluate the initial conditions at the initialization points
        !
        usp(1:nq,n1:n2) = initial_conditions( xyz_init(1:nr,1:np) , init_cond )
        !
        ! Interpolate the initial conditions from the initialization points
        ! to the solution points
        !
        do m = 1,nq
          usp(m,n1:n2) = matmul( init2numer , usp(m,n1:n2) )
        end do
        !
        ! Disassociate all local pointers before continuing
        !
        if (associated(numer2init)) numer2init => null()
        if (associated(init2numer)) init2numer => null()
        !
      end if
      !
    end do
    !
    ! Project the initial solution from the
    ! initialization points to the solution points
    !
   !call project_to_porder(cell,usp)
    !
    ! Initialize the iteration counter, time step, and elapsed time
    !
    d_t   = constant_dt
    itcur = 0
    time  = zero
    !
  end if
  !
  ! Create the character format strings for the
  ! restart file names and the linking commands
  !
  call create_restart_fname_formats
  !
  ! Initialize the time-averaging information if outputting
  ! time-averaged variables to solution/restart files
  !
 !if (output_time_averaging) then
 !  call initialize_time_ave
 !end if
  !
  ! Initialize the CFL cycles
  !
  call adjust_cfl_cycles
  !
  ! Deallocate the grid derived type
  !
  if (allocated(grid)) then
    deallocate ( grid , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Initialize the solution at previous timestep array
  !
  if (allocated(uoldsp)) then
    uoldsp(1:nq,1:n_solpts) = usp(1:nq,1:n_solpts)
  end if
  !
  ! Initialize the solution gradient array
  !
  call initialize_dusp
  !
  ! Output the boundary conditions
  !
  call write_bc_conditions
  !
  ! Output the free-stream conditions
  !
  if (mypnum == glb_root) then
    if (all(itestcase /= [Density_Transport,Entropy_Transport])) then
      !
      write (iout,10) rhoref,pref,tref,machref,aref,velref(1)
      if (nr >= 2) write (iout,11) velref(2)
      if (nr == 3) write (iout,12) velref(3)
      write (iout,13) alpha_aoaref,beta_aoaref,gam,rgasref,cpref,muref,reyref
      !
      write (iout,20) rhoref/rhoref,pref_nd,tref/tref,aref/aref,velref_nd(1)
      if (nr >= 2) write (iout,21) velref_nd(2)
      if (nr == 3) write (iout,22) velref_nd(3)
      write (iout,23) rgasref_nd,cpref_nd,muref_nd
      !
    else
      !
      write (iout,30) rhoref,pref,tref,machref,uref/machref, &
                      velref(1),velref(2),alpha_aoaref,beta_aoaref,gam, &
                      rgasref,muref,reyref
      !
    end if
    !
  end if
  !
  ! Deallocate local arrays and disassociate local pointers before leaving
  !
  if (allocated(outerp_inv)) then
    deallocate ( outerp_inv , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"outerp_inv",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(xyz_init)) then
    deallocate ( xyz_init , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_init",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (associated(numer2init)) numer2init => null()
  if (associated(init2numer)) init2numer => null()
  !
  ! Create the facexyz and edgexyz arrays if
  ! using the continuous output option or MMS
  !
  call get_face_and_edge_point_coordinates
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("outerp_inv(",a,",",i0,")%mat")
  !
  10 format (/, &
             " Dimensional Free-stream conditions",/, &
             " Density          = ",es14.6," kg/m^3",/, &
             " Pressure         = ",es14.6," N/m^2",/, &
             " Temperature      = ",es14.6," K",/, &
             " Mach Number      = ",es14.6,/, &
             " Speed of Sound   = ",es14.6," m/s",/, &
             " X-Velocity       = ",es14.6," m/s")
  11 format (" Y-Velocity       = ",es14.6," m/s")
  12 format (" Z-Velocity       = ",es14.6," m/s")
  13 format (" Alpha A.o.Attack = ",es14.6," degrees"/, &
             " Beta  A.o.Attack = ",es14.6," degrees"/, &
             " Gamma            = ",es14.6,/, &
             " Gas Constant     = ",es14.6," J/(kg*K)",/, &
             " Cp               = ",es14.6," J/(kg*K)",/, &
             " Viscosity        = ",es14.6," kg/(m*s)",/, &
             " Reynolds #       = ",es14.6,/)
  !
  20 format (/, &
             " Non-dimensional Free-stream conditions",/, &
             " Density         = ",es14.6,/, &
             " Pressure        = ",es14.6,/, &
             " Temperature     = ",es14.6,/, &
             " Speed of Sound  = ",es14.6,/, &
             " X-Velocity      = ",es14.6)
  21 format (" Y-Velocity      = ",es14.6)
  22 format (" Z-Velocity      = ",es14.6)
  23 format (" Gas Constant    = ",es14.6,/, &
             " Cp              = ",es14.6,/, &
             " Viscosity       = ",es14.6,/)
  !
  30 format (/, &
             " The density vortex transport problem is being",/, &
             " used to model the linear 2D advection equation.",//, &
             " The Non-dimensional Free-stream conditions are",/, &
             " Density          = ",es14.6,/, &
             " Pressure         = ",es14.6,/, &
             " Temperature      = ",es14.6,/, &
             " Mach Number      = ",es14.6,/, &
             " Speed of Sound   = ",es14.6,/, &
             " X-Velocity       = ",es14.6,/, &
             " Y-Velocity       = ",es14.6,/, &
             " Alpha A.o.Attack = ",es14.6," degrees"/, &
             " Beta  A.o.Attack = ",es14.6," degrees"/, &
             " Gamma            = ",es14.6,/, &
             " Gas Constant     = ",es14.6,/)
  !
end subroutine initialize_flow
!
!###############################################################################
!
pure &
function initial_conditions(xyz,init_cond) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : ncell
  use ovar, only : bc_in
  use ovar, only : itestcase
  use ovar, only : mms_opt,mms_init
  use ovar, only : postproc_debug_grid
  use mms_mod, only : mms_solution_at_xyz
  use channel_mod, only : channel_initialization
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:,:)
  integer,  intent(in) :: init_cond
  !
  !.. Function Result ..
  real(wp), dimension(1:nq,1:size(xyz,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: k,np,ierr
  real(wp) :: div,dl,r,r0,nozzle_len
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: adj_xyz(:,:)
  !
continue
  !
  np = size(xyz,dim=2)
  !
  if (any(itestcase == Transport_Problems)) then
    !
    do k = 1,np
      return_value(:,k) = vortex_solution( xyz(:,k) )
    end do
    !
  else if (itestcase == Taylor_Green_Vortex) then
    !
    adj_xyz = xyz ! F2003 AUTO-REALLOCATION
    !
    if (postproc_debug_grid) then
      dl = real(ncell,kind=wp)**(-one3)
      adj_xyz = PI * (two*dl*xyz - one)
    end if
    !
    do k = 1,np
      return_value(:,k) = taylor_green_initialization( adj_xyz(:,k) )
    end do
    !
  else if (itestcase == Channel_Flow) then
    !
    return_value(:,:) = channel_initialization( xyz(:,:) )
    !
  else if (mms_opt /= 0) then
    !
    do k = 1,np
      return_value(:,k) = mms_init * mms_solution_at_xyz( xyz(:,k) , &
                                                          return_nondim=true, &
                                                          return_cv=true )
    end do
    !
  else if (abs(itestcase) == Nozzle_Jet) then
    !
    if (init_cond == 2) then
      !
      nozzle_len = merge(9.875_wp,7.74_wp,itestcase > 0)
      !
      do k = 1,np
        return_value(:,k) = jet_profile2(nozzle_len,xyz(1,k),bc_in(0)%pv(:))
      end do
      !
    else if (init_cond == 3) then
      !
      nozzle_len = merge(9.875_wp,7.74_wp,itestcase > 0)
      !
      if (itestcase == Nozzle_Jet) then
        div = eps1*smc000_interior_wall_radius(one)
      else if (itestcase == -Nozzle_Jet) then
        div = eps1*arn2_interior_wall_radius(one)
      end if
      !
      do k = 1,np
        !
        if (itestcase == Nozzle_Jet) then
          r0 = smc000_interior_wall_radius(xyz(1,k))
        else if (itestcase == -Nozzle_Jet) then
          r0 = arn2_interior_wall_radius(xyz(1,k))
        end if
        r = norm2(xyz(2:,k))
        return_value(:,k) = jet_profile3(nozzle_len,xyz(1,k),r,r0,div, &
                                         bc_in(1)%pv(:), &
                                         bc_in(0)%pv(:))
        !
      end do
      !
    else if (init_cond /= 0) then
      !
      ! div = thickness parameter for shear layer
      !     = (radius of nozzle at the lip) * 0.1
      ! NOTE: The nozzle lip should be located at x=0.0 so sending any value
      !       greater than 0.0 into the function smc000_interior_wall_radius
      !       should return the radius of the nozzle at the lip.
      !
      if (itestcase == Nozzle_Jet) then
        div = eps1*smc000_interior_wall_radius(one)
      else if (itestcase == -Nozzle_Jet) then
        div = eps1*arn2_interior_wall_radius(one)
      end if
      !
      do k = 1,np
        !
        if (itestcase == Nozzle_Jet) then
          r0 = smc000_interior_wall_radius(xyz(1,k))
        else if (itestcase == -Nozzle_Jet) then
          r0 = arn2_interior_wall_radius(xyz(1,k))
        end if
        r = norm2(xyz(2:,k))
        return_value(:,k) = jet_profile(r,r0,div,bc_in(0)%pv(:), &
                                        bc_in(init_cond)%pv(:))
        !
      end do
      !
    else
      !
      do k = 1,np
        return_value(:,k) = bc_in(init_cond)%cv(:)
      end do
      !
    end if
    !
  else if (itestcase == Infinite_Cylinder) then
    !
    do k = 1,np
      return_value(:,k) = initialize_cylinder(xyz(:,k),bc_in(0)%pv(:))
    end do
    !
  else
    !
    do k = 1,np
      return_value(:,k) = bc_in(init_cond)%cv(:)
    end do
    !
  end if
  !
end function initial_conditions
!
!###############################################################################
!
pure function jet_profile(r,r0,div,fs_pv,bc_pv) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: r
  real(wp), intent(in) :: r0
  real(wp), intent(in) :: div
  real(wp), intent(in) :: fs_pv(:)
  real(wp), intent(in) :: bc_pv(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(fs_pv))
  !
  !.. Local Scalars ..
  real(wp) :: func
  !
continue
 !!
 !! Initial profile for a jet
 !!
 !!  r0     = radius of jet profile
 !!  div    = thickness parameter for shear layer
 !!  ujet   = jet velocity
 !!  rhojet = jet density
 !!
 !r0   = 1._wp/12._wp
 !div  = 0.1*r0
 !ujet = 558.31_wp       ! (feet / second)
 !rhojet = 2.5019e-03_wp ! (slugs / cubic-feet)
 !!
 !do k = 1,mk(nb)
 !  do j = 1,mj(nb)
 !    do i = 1,mi(nb)
 !      if (nb < 2) then
 !        r0 = sqrt(x(i,mj(nb),k,2)**2+x(i,mj(nb),k,3)**2)
 !      end if
 !      r = sqrt(x(i,j,k,2)**2+x(i,j,k,3)**2)
 !      func  = half*(one-erf((r-r0)/(div)))
 !      uinit   = ujet*func   + (one-func)*fsvel
 !      rhoinit = rhojet*func + (one-func)*fsrho
 !      q(i,j,k,1) = rhoinit
 !      q(i,j,k,2) = rhoinit*uinit
 !      q(i,j,k,3) = zero
 !      q(i,j,k,4) = zero
 !      q(i,j,k,5) = rhoinit*(pinf/(rhoinit*(gam-one)) + half*uinit**2)
 !    end do
 !  end do
 !end do
  !
  func = half * (one - erf((r-r0)/div))
  !
  return_value(nec) = func*bc_pv(nec) + (one-func)*fs_pv(nec)
  return_value(nmx) = func*bc_pv(nmx) + (one-func)*fs_pv(nmx)
  return_value(nmx+1:nme) = zero
 !return_value(nmz) = zero
  return_value(nee) = fs_pv(nee)
  !
  return_value(:) = vsp2u_sp( return_value(:) )
  !
end function jet_profile
!
!###############################################################################
!
pure function jet_profile3(nozzle_len,x,r,r0,div,bc_pv,fs_pv) &
                    result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : gam,grid_scaling_factor
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: nozzle_len
  real(wp), intent(in) :: x
  real(wp), intent(in) :: r
  real(wp), intent(in) :: r0
  real(wp), intent(in) :: div
  real(wp), intent(in) :: bc_pv(:)
  real(wp), intent(in) :: fs_pv(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(fs_pv))
  !
  !.. Local Scalars ..
  real(wp) :: xbeg,func,xfunc,rfunc
  !
  !.. Local Arrays ..
  real(wp) :: vel(1:nr)
  !
continue
  !
  ! Initialize to the free-stream primitive variables
  !
  return_value(:) = fs_pv(:)
  !
  ! If within the nozzle, initialize the velocity mach number to
  ! a minimum of 0.05 at the inflow plane and use the error function
  ! to smoothly transition this to the free-stream mach number within
  ! the interior of the nozzle.
  !
  if (x < zero) then
    !
    xbeg = nozzle_len*grid_scaling_factor
    !
    rfunc = half * (one - erf((r-r0)/div))
    xfunc = half * (one + erf((ten/xbeg)*abs(x)-seven))
    func = rfunc * xfunc
    !
    return_value(nec) = func*bc_pv(nec) + (one-func)*fs_pv(nec)
    return_value(nmx) = func*bc_pv(nmx) + (one-func)*fs_pv(nmx)
    return_value(nmx+1:nme) = zero
    return_value(nee) = fs_pv(nee)
    !
  end if
  !
  return_value(:) = vsp2u_sp( return_value(:) )
  !
end function jet_profile3
!
!###############################################################################
!
pure function jet_profile2(nozzle_len,x,fs_pv) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : gam,grid_scaling_factor
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: nozzle_len
  real(wp), intent(in) :: x
  real(wp), intent(in) :: fs_pv(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(fs_pv))
  !
  !.. Local Scalars ..
  real(wp) :: aspd,mach,xbeg,func
  !
  !.. Local Arrays ..
  real(wp) :: vel(1:nr)
  !
continue
  !
  ! Initialize to the free-stream primitive variables
  !
  return_value(:) = fs_pv(:)
  !
  ! If within the nozzle, initialize the velocity mach number to
  ! a minimum of 0.05 at the inflow plane and use the error function
  ! to smoothly transition this to the free-stream mach number within
  ! the interior of the nozzle.
  !
  if (x < zero) then
    !
    aspd = sqrt(max(gam*fs_pv(nee)/fs_pv(nec),rndoff))
    mach = norm2(fs_pv(nmb:nme)) / aspd
    !
    if (mach < 0.05_wp) then
      !
      xbeg = nozzle_len*grid_scaling_factor
      !
      vel(:) = 0.05_wp * aspd * unit_vector( fs_pv(nmb:nme) )
      !
      if (x <= -xbeg) then
        !
        return_value(nmb:nme) = vel(:)
        !
      else
        !
        func = half * (one - erf((ten/xbeg)*abs(x)-seven))
        !
        return_value(nmb:nme) = func*fs_pv(nmb:nme) + (one-func)*vel(:)
        !
      end if
      !
    end if
    !
  end if
  !
  return_value(:) = vsp2u_sp( return_value(:) )
  !
end function jet_profile2
!
!###############################################################################
!
subroutine find_elem_for_initial_jet_profile(point_is_in_nozzle_interior)
  !
  !.. Use Statements ..
  use geovar, only : grid
  !
  !.. Dummy Function Interface ..
  interface
    function point_is_in_nozzle_interior(xyz) result(return_value)
      import :: wp,lk
      real(wp), intent(in) :: xyz(:)
      logical(lk) :: return_value
    end function point_is_in_nozzle_interior
  end interface
  !
  !.. Local Scalars ..
  integer :: l,n,init
  !
  !.. Local Arrays ..
  real(wp) :: xyz_ave(1:3)
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: pts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_elem_for_initial_jet_profile"
  !
 !integer, parameter :: init_val = 1 ! radial error function profile
  integer, parameter :: init_val = 2 ! inflow mach number limiter
 !integer, parameter :: init_val = 3 ! options 1 and 2 combined
  !
continue
  !
  xyz_ave(:) = zero
  !
  do n = 1,size(grid%elem)
    !
    if (grid%elem(n)%tags(1) == 0) then
      !
      pts = grid%elem(n)%pts(:) ! USING F2003 AUTO-REALLOCATION
      !
      do l = 1,size(grid%xyz,dim=1)
        xyz_ave(l) = sum(grid%xyz(l,pts(:))) / real(size(pts),kind=wp)
      end do
      !
      if (xyz_ave(1) > zero) then
        if (any(init_val == [2,3])) then
          init = 0
        else
          init = init_val
        end if
      else
        init = merge(init_val,0,point_is_in_nozzle_interior(xyz_ave))
      end if
      !
      if (size(grid%elem(n)%tags) < 2) then
        grid%elem(n)%tags = [grid%elem(n)%tags(1),init] ! F2003 AUTO-REALLOC
      else
        grid%elem(n)%tags(2) = init
      end if
      !
    end if
    !
  end do
  !
end subroutine find_elem_for_initial_jet_profile
!
!###############################################################################
!
pure &
function initialize_cylinder(xyz,fs_pv) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:)
  real(wp), intent(in) :: fs_pv(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(fs_pv))
  !
  !.. Local Scalars ..
  real(wp) :: a,b,c,u,v,w,rho,pres
  real(wp) :: r,theta,uinf
  !
  !.. Local Parameters ..
  real(wp), parameter :: sig = one/ten
  real(wp), parameter :: three2 = three/two
  real(wp), parameter :: five4 = five/four
  real(wp), parameter :: halfpi = half*pi
  !
continue
  !
  r = norm2(xyz(1:2))
  theta = acos(xyz(1)/r)
  !
  rho  = fs_pv(nec)
  uinf = fs_pv(nmx)
  u    = uinf
  v    = fs_pv(nmy)
  w    = fs_pv(nmz)
  pres = fs_pv(nee)
  !
  if (r < three2) then
    !
    u = uinf * cylinder_function(r,one,half)
    !
    a = sin( halfpi * xyz(3) )
    !
    if (r > half .and. r <= one) then
      b = cylinder_function(r,three4,one4)
    else
      b = one - cylinder_function(r,five4,one4)
    end if
    !
    if (theta > zero .and. theta < pi) then
      c = cylinder_function(theta,halfpi,halfpi)
    else
      c = one - cylinder_function(theta,three*halfpi,halfpi)
    end if
    !
    v = sig * uinf * b * c
    w = sig * uinf * b * a*a
    !
  end if
  !
  return_value(nec) = rho
  return_value(nmx) = u
  return_value(nmy) = v
  return_value(nmz) = w
  return_value(nee) = pres
  !
  return_value(:) = vsp2u_sp( return_value(:) )
  !
end function initialize_cylinder
!
!###############################################################################
!
pure function cylinder_function(r,rm,ra) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: r,rm,ra
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: r_min_rm,numer,denom
  !
continue
  !
  r_min_rm = r - rm
  !
  numer = two * ra * r_min_rm
  denom = ra*ra - r_min_rm*r_min_rm
  !
  return_value = half * (one + tanh(numer/denom))
  !
end function cylinder_function
!
!###############################################################################
!
subroutine initialize_dusp()
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use order_mod,  only : maxSP
  use geovar,  only : nr,ncell,cell
  use metrics_mod,  only : metrics_dt
  use derivatives_mod, only : deriv
  use flowvar, only : usp,dusp
  !
  !.. Local Scalars ..
  integer :: k,l,m,n,nc,n1,n2,np
  integer :: this_geom,this_order
  integer :: i,i1,i2,i3
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: dudr
  real(wp), dimension(1:maxSP,1:nq) :: tusp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_dusp"
  !
 !logical(lk), parameter :: use_conservative_form = true
  logical(lk), parameter :: use_conservative_form = fals
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Loop over the cells and get the derivative of
  ! the primitive variables at each solution point
  !
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index for solpts in cell nc
    n2 = cell(nc)%end_sp ! ending    index for solpts in cell nc
    np = n2-n1+1 ! number of solution points in this cell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Store the transpose of usp for this cell to reduce cache misses
    !
    do k = 1,np
      do m = 1,nq
        tusp(k,m) = usp(m,k+n1-1)
      end do
    end do
    !
    ! Compute the gradient of the conservative variables at each solution point
    !
    do k = 1,np
      !
      n = k+n1-1 ! total index of solution point on current partition
      !
      i = this_order+1
      i3 = (k-1)/(i*i) + 1
      i2 = (k-1)/i - (i3-1)*i + 1
      i1 = k - (i2-1)*i - (i3-1)*i*i
      !
      do m = 1,nq
        !
        if (use_conservative_form) then
          !
          do l = 1,nr
            !
            dudr(1) = deriv(this_geom,this_order)%d(1)% &
                      dot( k , tusp(1:np,m)*metrics_dt(nc)%met(:,1,l) )
            !
            dudr(2) = deriv(this_geom,this_order)%d(2)% &
                      dot( k , tusp(1:np,m)*metrics_dt(nc)%met(:,2,l) )
            !
            dudr(3) = deriv(this_geom,this_order)%d(3)% &
                      dot( k , tusp(1:np,m)*metrics_dt(nc)%met(:,3,l) )
            !
            dusp(l,m,n) = sum(dudr)!* jacobian(n)
            !
          end do
          !
        else
          !
          ! Compute the derivatives in r-s space
          !
          dudr(1:nr) = deriv(this_geom,this_order)%grad( k , tusp(1:np,m) )
          !
          ! Transform the derivatives back to x-y space
          !
          do l = 1,nr
            dusp(l,m,n) = dot( dudr(1:nr) , metrics_dt(nc)%met(k,1:nr,l) )
          end do
          !
        end if
        !
      end do
      !
    end do
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine initialize_dusp
!
!###############################################################################
!
subroutine get_face_and_edge_point_coordinates
  !
  use order_mod,         only : maxSP,geom_solpts
  use geovar,            only : nr,cell,face,xyz
  use flowvar,           only : edgexyz,facexyz
  use interpolation_mod, only : interp
  !
  !.. Local Scalars ..
  integer :: n,nc,ne,nf,np,n1,n2
  integer :: k,kf,l,lf,ierr
  integer :: this_geom,this_order
  integer :: face_geom,face_order
  integer :: edge_order
  !
  !.. Local Arrays ..
  real(wp) :: txyz(1:maxSP,1:nr)
  !
  !.. Local Allocatable Arrays ..
  logical(lk), allocatable :: msk(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_face_and_edge_point_coordinates"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Only do this if facexyz is allocated which means that either
  ! the continuous output option or MMS is being used
  !
  if (allocated(facexyz)) then
    !
    do nf = 1,size(facexyz)
      !
      face_geom = face(nf)%geom
      face_order = face(nf)%order
      !
      nc = face(nf)%left%cell
      !
      n1 = cell(nc)%beg_sp
      n2 = cell(nc)%end_sp
      np = n2-n1+1
      !
      this_geom = cell(nc)%geom
      this_order = cell(nc)%order
      !
      ! Create transpose of xyz array to reduce cache misses
      !
      do k = 1,np
        do l = 1,nr
          txyz(k,l) = xyz(l,k+n1-1)
        end do
      end do
      !
      ! Interpolate the solution point coordinates to each flux point
      ! on the current face
      !
      do kf = 1,geom_solpts(face_geom,face_order)
        !
        lf = face(nf)%left%get_fp_idx(face_geom,face_order,kf)
        !
        do l = 1,nr
          facexyz(nf)%v(l,kf) = interp(this_geom,this_order)% &
                                          toFace(face_order)% &
                                dot( lf , txyz(1:np,l) )
        end do
        !
      end do
      !
    end do
    !
  end if
  !
  ! Only do this if edgexyz is allocated which means
  ! that the continuous output option is being used
  !
  if (allocated(edgexyz)) then
    !
    ! Initialize edgexyz to zero
    !
#ifdef SPECIAL_FOR_PGI
    call edgexyz%zero_ne2d_t
#else
    call edgexyz%zero
#endif
    !
    ! Allocate the mask array to mark which edges have been finished
    !
    allocate ( msk(1:size(edgexyz)) , source=fals , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"msk",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Loop through the cells and compute the coordinates
    ! of the edge points on each grid edge by interpolating
    ! the coordinates of the cell solution points
    !
    do nc = 1,size(cell)
      !
      n1 = cell(nc)%beg_sp
      n2 = cell(nc)%end_sp
      np = n2-n1+1
      !
      this_geom = cell(nc)%geom
      this_order = cell(nc)%order
      !
      ! Create transpose of xyz array to reduce cache misses
      !
      do k = 1,np
        do l = 1,nr
          txyz(k,l) = xyz(l,k+n1-1)
        end do
      end do
      !
      ! Interpolate the solution point coordinates to each edge point
      ! on the current edge
      !
      edges_of_cell: do n = 1,size(cell(nc)%edge)
        !
        ne = cell(nc)%edge(n)%idx
        !
        ! Skip to the next edge of this cell if this edge has already been done
        !
        if (msk(ne)) cycle edges_of_cell
        !
        edge_order = cell(nc)%edge(n)%order
        !
        ! Loop over the edge points on this edge and interpolate the
        ! coordinates of the cell solution points to each edge point
        !
        do kf = 1,edge_order+1
          !
          lf = cell(nc)%edge(n)%get_ep_idx(kf)
          !
          do l = 1,nr
            edgexyz(ne)%v(l,kf) = interp(this_geom,this_order)% &
                                            toEdge(edge_order)% &
                                  dot(lf,txyz(1:np,l))
          end do
          !
        end do
        !
        ! Mark this edge as finished
        !
        msk(ne) = true
        !
      end do edges_of_cell
      !
    end do
    !
    ! Deallocate the msk array before leaving
    !
    deallocate ( msk , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_face_and_edge_point_coordinates
!
!###############################################################################
!
subroutine initialize_reference
  !
  !.. Use Statements ..
  use ovar,   only : rhoref,aref,tref,gam
  use ovar,   only : velref,velref_nd,uref
  use ovar,   only : pref,pref_nd,tmin,tmax
  use ovar,   only : rgasref,rgasref_nd
  use ovar,   only : cpref,cpref_nd,reyref
  use ovar,   only : cv_in,pv_in,machref
  use ovar,   only : ttotref,ttotref_nd
  use ovar,   only : ptotref,ptotref_nd
  use ovar,   only : ptot2p_ratio
  use ovar,   only : suth_Sref,suth_Tref
  use ovar,   only : alpha_aoaref,beta_aoaref
  use ovar,   only : itestcase,Lref
  !
  !.. Local Scalars ..
  integer  :: ierr
  real(wp) :: rhoa,rhoa2,tda2
  real(wp) :: gam_over_gm1
  real(wp) :: gm1_over_gam
  real(wp) :: total_cond
  real(wp) :: pres_check
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_reference"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Set up the base reference conditions
  !
  if (any(itestcase == Advec_Problems) .or. &
      abs(itestcase) == Shu_Vortex .or. &
      itestcase == -Vortex_Transport) then
    !
    ! Set up the base reference conditions for different
    ! versions of the isentropic Euler vortex problem
    !
    ! Set some reference conditions common to all of these vortex problems
    !
    beta_aoaref = zero
    velref(:) = zero
    !
    if (itestcase == -Shu_Vortex) then
      !
      ! Stationary vortex problem
      !
      alpha_aoaref = 90.0_wp
      itestcase = abs(itestcase)
      !
    else if (itestcase == -Vortex_Transport) then
      !
      ! Alteration to make the Shu case propagate in the
      ! y-direction only instead of in both x and y directions
      !
      alpha_aoaref = 90.0_wp
      velref(2) = one
      itestcase = Shu_Vortex
      !
    else
      !
      ! Original Shu vortex problem with diagonal propagation
      !
      alpha_aoaref = 45.0_wp
      velref(1) = one
      velref(2) = one
      !
    end if
    !
    machref = one/sqrt(gam)
    aref = one
    uref = one
    rhoref = one
    pref = one
    tref = one
    rgasref = one
    !
  else if (itestcase == Taylor_Green_Vortex) then
    !
    ! If Lref is negative, the grid was not created on the fly and
    ! we need to add code to compute Lref and make sure the input
    ! grid is valid for the Taylor-Green vortex problem
    !
    if (Lref <= zero) then
      write (error_message,101)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
    !
    ! Set up the base reference conditions for the Taylor-Green Vortex problem
    !
    alpha_aoaref = zero
    beta_aoaref = zero
    !
    gam = 1.4_wp
    tref = 300.0_wp ! (K)
    pref = 101325.0_wp ! (Pa)
    ! ( universal gas constant / molecular weight of air )
    rgasref = 8314.34_wp / 28.9651159_wp ! (J)/(kg*K)
    rhoref = pref / rgasref / tref ! (kg/m^3) ~ 1.17663794
    !
    aref = sqrt( gam*pref/rhoref ) ! (m/s) ~ 347.2169357
    machref = 0.1_wp
    uref = aref * machref ! (m/s) ~ 34.72169357
    !
    if (reyref <= zero) then
      reyref = 1600.0_wp
    end if
    !
    suth_Sref = 110.5556_wp
    suth_Tref = 273.111_wp
    !
  else
    !
    ! Set up the base reference conditions based
    ! on the input flow reference conditions
    !
    ! Check the reference conditions for consistency
    !
    if (rhoref < zero) then
      rhoref = pref/(rgasref*tref)
    else if (pref < zero) then
      pref = rhoref*rgasref*tref
    else if (tref < zero) then
      tref = pref/(rgasref*rhoref)
    else
      pres_check = rhoref*rgasref*tref
      if (abs(pres_check-pref) > eps6) then
       !write (iout,102) pref,pres_check
        write (error_message,102) pref,pres_check
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
    end if
    !
    ! Compute some of the remaining reference conditions
    !
    aref = sqrt( gam * pref / rhoref )
    uref = machref * aref
    velref(1) = uref * cos(beta_aoaref  * degrees2radians) * &
                       cos(alpha_aoaref * degrees2radians)
    velref(2) = uref * cos(beta_aoaref  * degrees2radians) * &
                       sin(alpha_aoaref * degrees2radians)
    velref(3) = uref * sin(beta_aoaref  * degrees2radians)
    !
  end if
  !
  ! Pre-compute some reference constants that will be used
  !
  rhoa = rhoref * aref
  rhoa2 = rhoa * aref
  tda2 = tref / aref / aref
  gam_over_gm1 = gam / (gam-one)
  gm1_over_gam = one / gam_over_gm1
  !
  ! Allocate the free stream inflow arrays
  !
  allocate ( cv_in(1:nq) , source=zero , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cv_in",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  allocate ( pv_in(1:nq) , source=zero , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"pv_in",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  !
  ! Initialize the inflow boundary conditions to the reference conditions
  ! NOTE: This is hard-wired for 2D and needs to be generalized for 3D
  !
  pv_in(nec) = rhoref
  pv_in(nmb:nme) = velref
  pv_in(nee) = pref
  !
  cv_in(nec) = pv_in(nec)
  cv_in(nmb:nme) = pv_in(nec) * pv_in(nmb:nme)
  cv_in(nee) = pv_in(nee)/(gam-one) + half*pv_in(nec)*selfdot(pv_in(nmb:nme))
  !
  ! Nondimensional reference velocities
  velref_nd = velref / aref
  !
  ! Nondimensional reference pressure
  pref_nd = pref / rhoa2
  !
  ! Nondimensional reference gas constant
  rgasref_nd = rgasref * tda2
  !
  ! Dimensional reference specific heat at constant pressure
  cpref    = rgasref * gam / (gam - one)
  ! Nondimensional reference specific heat at constant pressure
  cpref_nd = cpref * tda2
  !
  ! Finally, nondimensionalize the inflow boundary conditions
  !
  cv_in(nec)     = cv_in(nec) / rhoref
  cv_in(nmb:nme) = cv_in(nmb:nme) / rhoa
  cv_in(nee)     = cv_in(nee) / rhoa2
  !
  pv_in(nec)     = pv_in(nec) / rhoref
  pv_in(nmb:nme) = pv_in(nmb:nme) / aref
  pv_in(nee)     = pv_in(nee) / rhoa2
  !
  ! Total pressure and temperature
  !
  if (ptot2p_ratio < zero) then
    !
    total_cond = one + half*(gam-one)*machref*machref
    !
    ttotref = tref * total_cond
    ptotref = pref * (total_cond**gam_over_gm1)
    !
    ptot2p_ratio = ptotref/pref
    !
  else
    !
    ttotref = tref * (ptot2p_ratio**gm1_over_gam)
    ptotref = pref * ptot2p_ratio
    !
  end if
  !
  ttotref_nd = ttotref / tref
  ptotref_nd = ptotref / rhoa2
  !
  tmin = 10.0_wp
  tmax = 15000.0_wp
  !
  ! Initialize the viscosity reference conditions
  !
  call initialize_viscosity
  !
  ! Initialize the turbulence reference conditions
  !
  call initialize_turbulence
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  101 format ( &
     " The grid being used for the Taylor-Green vortex problem was not", &
     " generated on the fly and seems to be input from a grid file. Code", &
     " still needs to be added to compute Lref and make sure the input", &
     " grid is valid for the Taylor-Green vortex problem.")
  102 format ( &
     " Free-stream conditions in the input file are not consistent!",/, &
     " Reference pressure given    : ",es14.6,/, &
     " Reference pressure computed : ",es14.6,/, &
     " Please correct the input file or set one of density, pressure,", &
     " and temperature to a negative value so it is explicitly computed.")
  !
end subroutine initialize_reference
!
!###############################################################################
!
subroutine initialize_viscosity
  !
  !.. Use Statements ..
  use ovar, only : muref, muref_nd, reyref
  use ovar, only : pv_in, rhoref, aref, machref, tref
  use ovar, only : suth_muref, suth_c1
  use ovar, only : suth_Tref, suth_Sref
  !
  !.. Local Scalars ..
  real(wp) :: muref_scale,tref_scale
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_viscosity"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  muref_scale = (one + suth_Sref/suth_Tref) / sqrt(suth_Tref)
  !
  if (reyref <= zero) then
    !
    suth_c1 = suth_muref * muref_scale
    !
    muref_nd = viscosity_pv_sp(pv_in)
    muref    = muref_nd * rhoref * aref
    !
    reyref = rhoref * machref * aref * reciprocal( muref )
    !
  else
    !
    muref = rhoref * machref * aref / reyref
    !
    muref_nd = muref / (rhoref * aref)
    !
    tref_scale = tref * sqrt(tref) / (tref + suth_Sref)
    !
    suth_muref = muref / muref_scale / tref_scale
    !
    suth_c1 = suth_muref * muref_scale
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine initialize_viscosity
!
!###############################################################################
!
subroutine initialize_turbulence()
  !
  !.. Use Statements ..
  use ovar, only : tkeref, tkeref_nd
  use ovar, only : tlsref, tlsref_nd
  use ovar, only : tvrref, tinref
  use ovar, only : rhoref, muref
  use ovar, only : aref, machref
  use ovar, only : cv_in, pv_in
  use ovar, only : turbulence_model
  use ovar, only : turb_is_on
  use ovar, only : include_tke_in_energy
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_turbulence"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Turbulence equations
  !
  if (turbulence_model /= laminar_flow) then
    !
    ! Set the flag indicating that turbulence is active
    !
    turb_is_on = true
    !
    ! Set flag that says whether TKE is included in the total energy
    !
    if (abs(turbulence_model) == k_omega) then
      include_tke_in_energy = true
    end if
    !
    ! Dimensional reference conditions
    !
    tkeref = 1.5_wp * (tinref*machref*aref)**2
    tlsref = rhoref * tkeref / (tvrref*muref)
    !
    ! Nondimensional reference conditions
    !
    tkeref_nd = tkeref / (rhoref * aref * aref)
    tlsref_nd = tlsref / (rhoref * aref)
    !
    ! Set the turbulence conserved variables inflow conditions
    !
    cv_in(ntk) = cv_in(nec) * tkeref_nd
    cv_in(ntl) = cv_in(nec) * tlsref_nd
    !
    ! Add the tke contribution to the energy equation
    !
    cv_in(nee) = cv_in(nee) + cv_in(ntk)
    !
    ! Recompute the primitive variable inflow conditions
    !
    pv_in(1:nq) =  usp2v_sp( cv_in(1:nq) )
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine initialize_turbulence
!
!###############################################################################
!
subroutine create_bc_in
  !
  !.. Use Statements ..
  use geovar, only : nr,bface
  use ovar,   only : bc_input,bc_in,bc_names
  use ovar,   only : tmin,tmax,tref,pref,rgasref
  use ovar,   only : velref,rhoref,aref,gam
  !
  !.. Local Scalars ..
  integer     :: i,ii,n,nf,nbc,ierr
  real(wp)    :: r,p,t,r0,p0,t0,mach,aspd,visc,reyn
  real(wp)    :: alpha,beta,rhoa2,total_cond
  real(wp)    :: gm1,gam_over_gm1,one_over_gm1
  logical(lk) :: static_is_set
  logical(lk) :: total_is_set
  logical(lk) :: mach_is_set
  logical(lk) :: vel_is_set
  !
  character(len=32) :: bcnam,bcinp
  !
  !.. Local Arrays ..
  real(wp) :: vel(1:nr)
  real(wp) :: pv(1:nq)
  real(wp) :: cv(1:nq)
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: bc_list(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_bc_in"
  !
continue
  !
  rhoa2 = rhoref * aref * aref
  gm1 = gam - one
  one_over_gm1 = one / gm1
  gam_over_gm1 = gam*one_over_gm1
  !
  ! Find any boundary conditions that were set in the
  ! input namelist file using the bc_input array
  !
  do n = 1,size(bc_input)
    !
    bc_input(n)%name = adjustl(bc_input(n)%name)
    !
    if (len_trim(bc_input(n)%name) > 0) then
      if (.not. allocated(bc_list)) then
        allocate ( bc_list(1:1) , source=n , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"bc_list",1,__LINE__,__FILE__,ierr,error_message)
      else
        bc_list = [ bc_list , n ] ! Using F2003 auto-reallocation
      end if
    end if
    !
  end do
  !
  nbc = 0
  if (allocated(bc_list)) nbc = size(bc_list)
  !
  allocate ( bc_in(0:nbc) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bc_in",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Set the reference boundary conditions
  !
  pv(nec) = rhoref / rhoref
  pv(nmb:nme) = velref(:) / aref
  pv(nee) = pref / rhoa2
  !
  cv(:) = vsp2u_sp(pv)
  !
  bc_in(0)%pv = pv(:) ! Using F2003 auto-reallocation
  bc_in(0)%cv = cv(:) ! Using F2003 auto-reallocation
  !
  ! Set the other boundary conditions
  !
  do i = 1,nbc
    !
    n = bc_list(i)
    !
    ! Save the index within bc_input for this bc_in
    !
    bc_in(i)%bc_input_idx = n
    !
    if ( bc_input(n)%has_default_flow_conditions() ) then
      !
      if (mypnum == glb_root) then
        write (iout,2) i,trim(adjustl(bc_input(n)%name)), &
          trim(adjustl(bc_integer_to_string(bc_input(n)%set_bc_value)))
      end if
      !
      ! Set the boundary conditions for this BC to the reference conditions
      !
      bc_in(i)%pv = bc_in(0)%pv(:) ! Using F2003 auto-reallocation
      bc_in(i)%cv = bc_in(0)%cv(:) ! Using F2003 auto-reallocation
      !
    else
      !
      static_is_set = fals
      total_is_set = fals
      mach_is_set = fals
      vel_is_set = fals
      !
      if (mypnum == glb_root) then
        write (iout,3) i,trim(adjustl(bc_input(n)%name))
      end if
      !
      ! Copy over the relaxation information
      !
      bc_in(i)%relax%value = min(one,max(bc_input(n)%relax%value,zero))
      bc_in(i)%relax%time = max(bc_input(n)%relax%time,zero)
      bc_in(i)%relax%iter = max(bc_input(n)%relax%iter,0)
      !
      if ( (bc_input(n)%wall_temp >= tmin) .and. &
           (bc_input(n)%wall_temp <= tmax) ) then
        !
        bc_in(i)%wall_temperature = bc_input(n)%wall_temp/tref
        !
        if (size(vel) == 1) then
          vel(:) = [ bc_input(n)%vx ]
        else if (size(vel) == 2) then
          vel(:) = [ bc_input(n)%vx , bc_input(n)%vy ]
        else if (size(vel) == 3) then
          vel(:) = [ bc_input(n)%vx , bc_input(n)%vy , bc_input(n)%vz ]
        end if
        !
        pv(nec) = bc_in(0)%pv(nec)
        pv(nmb:nme) = vel(1:nr) / aref
        pv(nee) = bc_in(0)%pv(nee)
        !
        cv(:) = vsp2u_sp(pv)
        !
        bc_in(i)%pv = pv(:) ! USING F2003 AUTO-REALLOCATION
        bc_in(i)%cv = cv(:) ! USING F2003 AUTO-REALLOCATION
        !
      else
        !
        r = bc_input(n)%rho_static
        p = bc_input(n)%p_static
        t = bc_input(n)%t_static
        !
        r0 = bc_input(n)%rho_total
        p0 = bc_input(n)%p_total
        t0 = bc_input(n)%t_total
        !
        mach  = bc_input(n)%mach
        alpha = bc_input(n)%alpha_aoa
        beta  = bc_input(n)%beta_aoa
        !
        if (size(vel) == 1) then
          vel(:) = [ bc_input(n)%vx ]
        else if (size(vel) == 2) then
          vel(:) = [ bc_input(n)%vx , bc_input(n)%vy ]
        else if (size(vel) == 3) then
          vel(:) = [ bc_input(n)%vx , bc_input(n)%vy , bc_input(n)%vz ]
        end if
        !
        iterate_bc_conditions: do ii = 1,5
          !
          mach_is_set = ( abs(mach) < five*five .and. mach > zero )
          vel_is_set  = ( norm2(vel) > zero )
          !
          call igl(r,p,t,static_is_set,n)
         !call dump_bc_conditions(ii,r,p,t,r0,p0,t0,vel,mach,alpha,beta, &
         !                        static_is_set,total_is_set, &
         !                        vel_is_set,mach_is_set,"igl(r,p,t)")
          call igl(r0,p0,t0,total_is_set,n)
         !call dump_bc_conditions(ii,r,p,t,r0,p0,t0,vel,mach,alpha,beta, &
         !                        static_is_set,total_is_set, &
         !                        vel_is_set,mach_is_set,"igl(r0,p0,t0)")
          !
          call mach_total_pres(p,p0,mach,mach_is_set,n)
         !call dump_bc_conditions(ii,r,p,t,r0,p0,t0,vel,mach,alpha,beta, &
         !                        static_is_set,total_is_set, &
         !                        vel_is_set,mach_is_set,"mach_total_pres")
          call mach_total_temp(t,t0,mach,mach_is_set,n)
         !call dump_bc_conditions(ii,r,p,t,r0,p0,t0,vel,mach,alpha,beta, &
         !                        static_is_set,total_is_set, &
         !                        vel_is_set,mach_is_set,"mach_total_temp")
          call mach_total_dens(r,r0,mach,mach_is_set,n)
         !call dump_bc_conditions(ii,r,p,t,r0,p0,t0,vel,mach,alpha,beta, &
         !                        static_is_set,total_is_set, &
         !                        vel_is_set,mach_is_set,"mach_total_dens")
          !
          call vel_from_static(r,p,t,vel,mach,alpha,beta,vel_is_set,mach_is_set)
         !call dump_bc_conditions(ii,r,p,t,r0,p0,t0,vel,mach,alpha,beta, &
         !                        static_is_set,total_is_set, &
         !                        vel_is_set,mach_is_set,"vel_from_static")
          !
          if (static_is_set .and. vel_is_set) exit iterate_bc_conditions
          !
          if (ii == 5) then
            write (error_message,1) n
            call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
          end if
          !
        end do iterate_bc_conditions
        !
        pv(nec) = r
        pv(nmb:nme) = vel
        pv(nee) = p
        !
        call check_bc_input(bc_input(n),pv,n)
        !
        pv(nec) = pv(nec) / rhoref
        pv(nmb:nme) = pv(nmb:nme) / aref
        pv(nee) = pv(nee) / rhoa2
        !
        cv(:) = vsp2u_sp(pv)
        !
        bc_in(i)%pv = pv(:) ! Using F2003 auto-reallocation
        bc_in(i)%cv = cv(:) ! Using F2003 auto-reallocation
        !
        if (mypnum == glb_root) then
          !
          r = bc_in(i)%pv(nec)*rhoref
          p = bc_in(i)%pv(nee)*rhoa2
          t = temperature_pv_sp(bc_in(i)%pv)*tref
          !
          aspd = sqrt(gam*rgasref*t)
          mach = aref*norm2(bc_in(i)%pv(nmb:nme)) / aspd
          !
          total_cond = one + half*gm1*mach*mach
          !
          r0 = r * total_cond**one_over_gm1
          p0 = p * total_cond**gam_over_gm1
          t0 = t * total_cond
          !
          alpha = zero
          beta = zero
          if (nr >= 2) then
            alpha = atan2(bc_in(i)%pv(nmy),bc_in(i)%pv(nmx))
            if (nr == 3) then
              beta = asin(bc_in(i)%pv(nmz)/norm2(bc_in(i)%pv(nmb:nme)))
            end if
          end if
          alpha = alpha * radians2degrees
          beta = beta * radians2degrees
          !
          visc = viscosity_pv_sp( [r,vel,p] )
          !
          reyn = norm2(vel)*r/visc
          !
          write (iout,4) r,(bc_in(i)%pv(ii)*aref,ii=nmb,nme)
          write (iout,5) p,t,mach,aspd,r0,p0,t0,alpha,beta,visc,reyn
          !
        end if
        !
      end if
      !
    end if
    !
    ii = 0
    bface_loop: do nf = 1,size(bface,dim=2)
      !
      if (any(bface(1,nf) == bc_comm) .or. bface(5,nf) == 0) cycle bface_loop
      !
      bcnam = trim(adjustl( uppercase( bc_names(bface(5,nf)) ) ))
      bcinp = trim(adjustl( uppercase( bc_input(n)%name      ) ))
      !
      if (bcnam == bcinp) then
        !
        ii = ii + 1
        !
        ! Save the index for this boundary condition
        ! within the bc_in array to bface(6,:)
        !
        bface(6,nf) = i
        !
        ! Force the boundary condition to that specified in the input file
        ! if the value of set_bc_value is not equal to bc_unknown
        !
        if (bc_input(n)%set_bc_value /= bc_unknown) then
          bface(1,nf) = bc_input(n)%set_bc_value
        end if
        !
      end if
      !
    end do bface_loop
    !
    call mpi_allreduce(MPI_IN_PLACE,ii,1_int_mpi,mpi_inttyp, &
                       MPI_SUM,MPI_COMM_WORLD,mpierr)
    !
    if (mypnum == glb_root) then
      write (iout,6) ii
    end if
    !
  end do
  !
  if (ncpu > 1) then
    if (mypnum == glb_root) flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Format Statements
  !
  1 format ("Error iterating on BC # ",i0," !!!")
  2 format (/," Setting boundary condition for BC # ",i0," - '",a,"'",/, &
              "    to ",a)
  3 format (/," Flow conditions for BC # ",i0," - '",a,"'")
  4 format (4x,"Static Density        = ",es14.6," kg/m^3",/, &
            4x,"X-velocity            = ",es14.6," m/s",:,/, &
            4x,"Y-velocity            = ",es14.6," m/s",:,/, &
            4x,"Z-velocity            = ",es14.6," m/s")
  5 format (4x,"Static Pressure       = ",es14.6," N/m^2",/, &
            4x,"Static Temperature    = ",es14.6," K",/, &
            4x,"Mach Number           = ",es14.6,/, &
            4x,"Speed of Sound        = ",es14.6," m/s",/, &
            4x,"Total Density         = ",es14.6," kg/m^3",/, &
            4x,"Total Pressure        = ",es14.6," N/m^2",/, &
            4x,"Total Temperature     = ",es14.6," K",/, &
            4x,"Alpha Angle of Attack = ",es14.6," degrees",/, &
            4x,"Beta  Angle of Attack = ",es14.6," degrees",/, &
            4x,"Viscosity             = ",es14.6," kg/(m*s)",/, &
            4x,"Reynolds #            = ",es14.6,/)
  6 format (4x,"Number of boundary faces for this BC = ",i0)
  !
end subroutine create_bc_in
!
!###############################################################################
!
subroutine write_bc_conditions
  !
  !.. Use Statements ..
  use geovar, only : nr,bface
  use ovar,   only : bc_input,bc_in
  use ovar,   only : tmin,tmax,tref,rgasref
  use ovar,   only : rhoref,aref,gam
  !
  !.. Local Scalars ..
  integer     :: i,ii,nf,bc_type
  real(wp)    :: r,p,t,r0,p0,t0,mach,aspd,visc,reyn
  real(wp)    :: alpha,beta,rhoa2,total_cond
  real(wp)    :: gm1,gam_over_gm1,one_over_gm1
  !
  !.. Local Arrays ..
  real(wp) :: vel(1:nr)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_bc_in"
  !
continue
  !
  rhoa2 = rhoref * aref * aref
  gm1 = gam - one
  one_over_gm1 = one / gm1
  gam_over_gm1 = gam*one_over_gm1
  !
  ! Output all the boundary conditions
  !
  do i = lbound(bc_in,dim=1),ubound(bc_in,dim=1)
    !
    if (mypnum == glb_root) then
      !
      if (i == 0) then
        write (iout,1)
      else
        !
        bc_type_loop: do nf = 1,size(bface,dim=2)
          if (any(bface(1,nf) == bc_comm) .or. bface(5,nf) == 0) then
            cycle bc_type_loop
          end if
          if (bface(6,nf) == i) then
            bc_type = bface(1,nf)
            exit bc_type_loop
          end if
        end do bc_type_loop
        !
        write (iout,2) i,trim(adjustl(bc_input(bc_in(i)%bc_input_idx)%name)), &
                         trim(adjustl(bc_integer_to_string(bc_type)))
        !
      end if
      !
      if ( (bc_in(i)%wall_temperature*tref >= tmin) .and. &
           (bc_in(i)%wall_temperature*tref <= tmax) ) then
        !
        vel(1:nr) = bc_in(i)%pv(nmb:nme)*aref
        !
        write (iout,3) bc_in(i)%wall_temperature*tref, &
                       (vel(ii),ii=1,nr)
        !
      else
        !
        r = bc_in(i)%pv(nec)*rhoref
        p = bc_in(i)%pv(nee)*rhoa2
        t = temperature_pv_sp(bc_in(i)%pv)*tref
        !
        vel(1:nr) = bc_in(i)%pv(nmb:nme)*aref
        !
        aspd = sqrt(gam*rgasref*t)
        mach = aref*norm2(bc_in(i)%pv(nmb:nme)) / aspd
        !
        total_cond = one + half*gm1*mach*mach
        !
        r0 = r * total_cond**one_over_gm1
        p0 = p * total_cond**gam_over_gm1
        t0 = t * total_cond
        !
        alpha = zero
        beta = zero
        if (nr >= 2) then
          alpha = atan2(bc_in(i)%pv(nmy),bc_in(i)%pv(nmx))
          if (nr == 3) then
            beta = asin(bc_in(i)%pv(nmz)/norm2(bc_in(i)%pv(nmb:nme)))
          end if
        end if
        alpha = alpha * radians2degrees
        beta = beta * radians2degrees
        !
        visc = viscosity_pv_sp( [r,vel,p] )
        !
        reyn = norm2(vel)*r/visc
        !
        write (iout,4) r,(vel(ii),ii=1,nr)
        write (iout,5) p,t,mach,aspd,r0,p0,t0,alpha,beta,visc,reyn
        !
      end if
      !
      flush (iout)
      !
    end if
    !
    ii = 0
    bface_loop: do nf = 1,size(bface,dim=2)
      !
      if (any(bface(1,nf) == bc_comm) .or. bface(5,nf) == 0) cycle bface_loop
      !
      if (bface(6,nf) == i) then
        ii = ii + 1
      end if
      !
    end do bface_loop
    !
    call mpi_allreduce(MPI_IN_PLACE,ii,1_int_mpi,mpi_inttyp, &
                       MPI_SUM,MPI_COMM_WORLD,mpierr)
    !
    if (mypnum == glb_root) then
      write (iout,6) ii
      flush (iout)
    end if
    !
  end do
  !
  if (ncpu > 1) then
    if (mypnum == glb_root) flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Format Statements
  !
  1 format (/," Flow conditions for BC # 0 - Reference/Freestream BCs")
  2 format (/," Flow conditions for BC # ",i0," - '",a,"'",/, &
              4x,"BC Type = ",a)
  3 format (4x,"Wall Temperature      = ",es14.6," K",/, &
            4x,"X-velocity            = ",es14.6," m/s",/,:, &
            4x,"Y-velocity            = ",es14.6," m/s",/,:, &
            4x,"Z-velocity            = ",es14.6," m/s",/)
  4 format (4x,"Static Density        = ",es14.6," kg/m^3",/, &
            4x,"X-velocity            = ",es14.6," m/s",:,/, &
            4x,"Y-velocity            = ",es14.6," m/s",:,/, &
            4x,"Z-velocity            = ",es14.6," m/s")
  5 format (4x,"Static Pressure       = ",es14.6," N/m^2",/, &
            4x,"Static Temperature    = ",es14.6," K",/, &
            4x,"Mach Number           = ",es14.6,/, &
            4x,"Speed of Sound        = ",es14.6," m/s",/, &
            4x,"Total Density         = ",es14.6," kg/m^3",/, &
            4x,"Total Pressure        = ",es14.6," N/m^2",/, &
            4x,"Total Temperature     = ",es14.6," K",/, &
            4x,"Alpha Angle of Attack = ",es14.6," degrees",/, &
            4x,"Beta  Angle of Attack = ",es14.6," degrees",/, &
            4x,"Viscosity             = ",es14.6," kg/(m*s)",/, &
            4x,"Reynolds #            = ",es14.6,/)
  6 format (4x,"Number of boundary faces for this BC = ",i0)
  !
end subroutine write_bc_conditions
!
!###############################################################################
!
subroutine dump_bc_conditions(ii,r,p,t,r0,p0,t0,vel,mach,alpha,beta, &
                              static_is_set,total_is_set, &
                              vel_is_set,mach_is_set,string)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: ii
  real(wp),         intent(in) :: r,r0,p,p0,t,t0
  real(wp),         intent(in) :: mach,alpha,beta
  real(wp),         intent(in) :: vel(:)
  logical(lk),      intent(in) :: static_is_set,total_is_set
  logical(lk),      intent(in) :: vel_is_set,mach_is_set
  character(len=*), intent(in) :: string
  !
  !.. Local Scalars ..
  integer :: i
  !
continue
  !
  if (mypnum /= 0) return
  !
  write (iout,1) trim(adjustl(string)),ii, &
                 r,p,t,static_is_set, &
                 r0,p0,t0,total_is_set, &
                 sqrt(merge(gam*p/r,zero,gam*p/r>zero)), &
                 mach,mach_is_set
  write (iout,2) (vel(i),i=1,size(vel))
  write (iout,3) alpha,beta,vel_is_set
  flush (iout)
  !
  ! Format Statements
  !
  1 format (/, &
            3x,"Dumping variables at location '",a, &
               "' for BC iteration # ",i0,/, &
            5x,"Static Density        = ",es14.6,/, &
            5x,"Static Pressure       = ",es14.6,/, &
            5x,"Static Temperature    = ",es14.6,/, &
            5x,"static_is_set         = ",l1,/, &
            5x,"Total Density         = ",es14.6,/, &
            5x,"Total Pressure        = ",es14.6,/, &
            5x,"Total Temperature     = ",es14.6,/, &
            5x,"total_is_set          = ",l1,/, &
            5x,"Speed of Sound        = ",es14.6,/, &
            5x,"Mach Number           = ",es14.6,/, &
            5x,"mach_is_set           = ",l1)
  2 format (5x,"X-velocity            = ",es14.6,:,/, &
            5x,"Y-velocity            = ",es14.6,:,/, &
            5x,"Z-velocity            = ",es14.6)
  3 format (5x,"Alpha Angle of Attack = ",es14.6,/, &
            5x,"Beta  Angle of Attack = ",es14.6,/, &
            5x,"vel_is_set            = ",l1,/)
  !
end subroutine dump_bc_conditions
!
!###############################################################################
!
subroutine check_bc_input(bc,pv,n)
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : rgasref,gam
  use ovar,   only : bc_input_t
  !
  !.. Formal Arguments ..
  type(bc_input_t), intent(in) :: bc
  real(wp),         intent(in) :: pv(:)
  integer,          intent(in) :: n
  !
  !.. Local Scalars ..
  real(wp)    :: r,p,t,r0,p0,t0,vmag,bc_vel,mach,a,alpha,beta
  real(wp)    :: gm1,one_over_gm1,gam_over_gm1,err_tol
  real(wp)    :: total_cond
  logical(lk) :: errors_found
  !
  !.. Local Arrays ..
  real(wp) :: vel(1:nr)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_bc_input"
  !
continue
  !
  errors_found = fals
  err_tol = eps6
  !
  gm1 = gam - one
  one_over_gm1 = one / gm1
  gam_over_gm1 = gam / gm1
  !
  r = pv(nec)
  vel(:) = pv(nmb:nme)
  vmag = norm2(vel)
  p = pv(nee)
  t = p / (rgasref*r)
  a = sqrt(gam*rgasref*t)
  mach = vmag/a
  !
  total_cond = one + half*gm1*mach*mach
  r0 = r * (total_cond**one_over_gm1)
  p0 = p * (total_cond**gam_over_gm1)
  t0 = t *  total_cond
  !
  alpha = zero
  beta = zero
  if (size(vel) == 3) beta = asin(vel(3)/vmag) * radians2degrees
  if (size(vel) >= 2) alpha = atan2(vel(2),vel(1)) * radians2degrees
  !
  ! Check static temperature
  !
  if (bc%t_static > zero) then
    if (abs(bc%t_static - t) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,1) bc%t_static,t
      end if
      errors_found = true
    end if
  end if
  !
  ! Check static pressure
  !
  if (bc%p_static > zero) then
    if (abs(bc%p_static - p) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,2) bc%p_static,p
      end if
      errors_found = true
    end if
  end if
  !
  ! Check static density
  !
  if (bc%rho_static > zero) then
    if (abs(bc%rho_static - r) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,3) bc%rho_static,r
      end if
      errors_found = true
    end if
  end if
  !
  ! Check total temperature
  !
  if (bc%t_total > zero) then
    if (abs(bc%t_total - t0) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,4) bc%t_total,t0
      end if
      errors_found = true
    end if
  end if
  !
  ! Check total pressure
  !
  if (bc%p_total > zero) then
    if (abs(bc%p_total - p0) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,5) bc%p_total,p0
      end if
      errors_found = true
    end if
  end if
  !
  ! Check total density
  !
  if (bc%rho_total > zero) then
    if (abs(bc%rho_total - r0) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,6) bc%rho_total,r0
      end if
      errors_found = true
    end if
  end if
  !
  ! Check mach number
  !
  if (bc%mach > zero) then
    if (abs(bc%mach - mach) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,7) bc%mach,mach
      end if
      errors_found = true
    end if
  end if
  !
  ! Check alpha angle of attack
  !
  if (bc%alpha_aoa /= zero) then
    if (abs(bc%alpha_aoa - alpha) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,8) bc%alpha_aoa,alpha
      end if
      errors_found = true
    end if
  end if
  !
  ! Check beta angle of attack
  !
  if (bc%beta_aoa /= zero) then
    if (abs(bc%beta_aoa - beta) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,9) bc%beta_aoa,beta
      end if
      errors_found = true
    end if
  end if
  !
  ! Check velocity
  !
  bc_vel = norm2( [ bc%vx , bc%vy , bc%vz ] )
  !
  if (bc_vel > zero) then
    if (abs(bc_vel - vmag) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,10) bc_vel,vmag
      end if
      errors_found = true
    end if
    if (abs(bc%vx - vel(1)) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,11) bc%vx,vel(1)
      end if
      errors_found = true
    end if
    if (size(vel) >= 2) then
      if (abs(bc%vy - vel(2)) > err_tol) then
        if (mypnum == glb_root) then
          if (.not.errors_found) write (iout,100) n
          write (iout,12) bc%vy,vel(2)
        end if
        errors_found = true
      end if
    end if
    if (size(vel) == 3) then
      if (abs(bc%vz - vel(3)) > err_tol) then
        if (mypnum == glb_root) then
          if (.not.errors_found) write (iout,100) n
          write (iout,13) bc%vz,vel(3)
        end if
        errors_found = true
      end if
    end if
  end if
  !
  if (errors_found) then
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  !
  ! Format Statements
  !
  1   format (5x,"Input       static temperature = ",es14.6,/, &
              5x,"Computed    static temperature = ",es14.6)
  2   format (5x,"Input          static pressure = ",es14.6,/, &
              5x,"Computed       static pressure = ",es14.6)
  3   format (5x,"Input           static density = ",es14.6,/, &
              5x,"Computed        static density = ",es14.6)
  4   format (5x,"Input        total temperature = ",es14.6,/, &
              5x,"Computed     total temperature = ",es14.6)
  5   format (5x,"Input           total pressure = ",es14.6,/, &
              5x,"Computed        total pressure = ",es14.6)
  6   format (5x,"Input            total density = ",es14.6,/, &
              5x,"Computed         total density = ",es14.6)
  7   format (5x,"Input              Mach number = ",es14.6,/, &
              5x,"Computed           Mach number = ",es14.6)
  8   format (5x,"Input    alpha angle of attack = ",es14.6,/, &
              5x,"Computed alpha angle of attack = ",es14.6)
  9   format (5x,"Input    beta  angle of attack = ",es14.6,/, &
              5x,"Computed beta  angle of attack = ",es14.6)
  10  format (5x,"Input       velocity magnitude = ",es14.6,/, &
              5x,"Computed    velocity magnitude = ",es14.6)
  11  format (5x,"Input               x-velocity = ",es14.6,/, &
              5x,"Computed            x-velocity = ",es14.6)
  12  format (5x,"Input               y-velocity = ",es14.6,/, &
              5x,"Computed            y-velocity = ",es14.6)
  13  format (5x,"Input               z-velocity = ",es14.6,/, &
              5x,"Computed            z-velocity = ",es14.6)
  100 format (/,3x,"An error was found checking the final ", &
                   "flow conditions for BC # ",i0," !!!",/)
  !
end subroutine check_bc_input
!
!###############################################################################
!
subroutine vel_from_static(r,p,t,vel,mach,alpha,beta,vel_is_set,mach_is_set)
  !
  !.. Use Statements ..
  use ovar, only : rgasref,gam
  !
  !.. Formal Arguments ..
  real(wp),       intent(in) :: r,p,t
  real(wp),    intent(inout) :: vel(:)
  real(wp),    intent(inout) :: mach
  real(wp),    intent(inout) :: alpha
  real(wp),    intent(inout) :: beta
  logical(lk), intent(inout) :: vel_is_set
  logical(lk), intent(inout) :: mach_is_set
  !
  !.. Local Scalars ..
  real(wp) :: aspd
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "vel_from_static"
  !
continue
  !
  aspd = -huge(zero)
  if (t > zero) then
    aspd = sqrt(gam*rgasref*t)
  else if (r > zero .and. p > zero) then
    aspd = sqrt(gam*p/r)
  else if (vel_is_set .and. mach_is_set) then
    aspd = norm2(vel) / mach
  end if
  !
  if (vel_is_set) then
    if (size(vel) == 3) beta = asin(vel(3)/norm2(vel)) * radians2degrees
    if (size(vel) >= 2) alpha = atan2(vel(2),vel(1)) * radians2degrees
    if (aspd > zero .and. .not. mach_is_set) then
      mach = norm2(vel)/aspd
      mach_is_set = true
    end if
  else
    if (aspd > zero .and. mach_is_set) then
      vel(:) = mach*aspd
      if (size(vel) >= 2) then
        vel(1) = vel(1)*cos(alpha*degrees2radians)
        vel(2) = vel(2)*sin(alpha*degrees2radians)
      end if
      if (size(vel) == 3) then
        vel(1) = vel(1)*cos(beta *degrees2radians)
        vel(2) = vel(2)*cos(beta *degrees2radians)
        vel(3) = vel(3)*sin(beta *degrees2radians)
      end if
      vel_is_set = true
    end if
  end if
  !
end subroutine vel_from_static
!
!###############################################################################
!
subroutine mach_total_temp(t,t0,mach,mach_is_set,n)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp),    intent(inout) :: t,t0
  real(wp),    intent(inout) :: mach
  logical(lk), intent(inout) :: mach_is_set
  integer,        intent(in) :: n
  !
  !.. Local Scalars ..
  real(wp) :: tmach,total_cond
  real(wp) :: gm1,two_over_gm1
  real(wp) :: tratio
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "mach_total_temp"
  !
continue
  !
  gm1 = gam - one
  !
  if (t > zero .and. t0 > zero) then
    !
    two_over_gm1 = two / gm1
    !
    tratio = t0/t
    !
    tmach = sqrt( two_over_gm1 * (tratio - one) )
    !
    if (mach_is_set) then
      if (abs(mach-tmach) > eps6) then
        write (error_message,1) n,mach,tmach
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
    else
      mach = tmach
      mach_is_set = true
    end if
    !
  else
    !
    if (mach_is_set) then
      total_cond = one + half*gm1*mach*mach
      if (t > zero) then
        t0 = t * total_cond
      else if (t0 > zero) then
        t = t0 / total_cond
      end if
    end if
    !
  end if
  !
  ! Format Statements
  !
  1 format ( &
     " Error with BC # ",i0," !",/, &
     " Mach number specified for boundary condition = ",es14.6,/, &
     " Mach number computed from the stagnation temperature = ",es14.6,/, &
     " Please correct the input values for this BC!")
  !
end subroutine mach_total_temp
!
!###############################################################################
!
subroutine mach_total_pres(p,p0,mach,mach_is_set,n)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp),    intent(inout) :: p,p0
  real(wp),    intent(inout) :: mach
  logical(lk), intent(inout) :: mach_is_set
  integer,        intent(in) :: n
  !
  !.. Local Scalars ..
  real(wp) :: pmach,total_cond
  real(wp) :: two_over_gm1,gm1
  real(wp) :: gam_over_gm1,pratio
  real(wp) :: gm1_over_gam
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "mach_total_pres"
  !
continue
  !
  gm1 = gam - one
  !
  if (p > zero .and. p0 > zero) then
    !
    gm1_over_gam = gm1 / gam
    two_over_gm1 = two / gm1
    !
    pratio = (p0/p)**gm1_over_gam
    !
    pmach = sqrt( two_over_gm1 * (pratio - one) )
    !
    if (mach_is_set) then
      if (abs(mach-pmach) > eps6) then
        write (error_message,1) n,mach,pmach
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
    else
      mach = pmach
      mach_is_set = true
    end if
    !
  else
    !
    if (mach_is_set) then
      gam_over_gm1 = gam / gm1
      total_cond = (one + half*gm1*mach*mach)**gam_over_gm1
      if (p > zero) then
        p0 = p * total_cond
      else if (p0 > zero) then
        p = p0 / total_cond
      end if
    end if
    !
  end if
  !
  ! Format Statements
  !
  1 format ( &
     " Error with BC # ",i0," !",/, &
     " Mach number specified for boundary condition = ",es14.6,/, &
     " Mach number computed from the stagnation pressure = ",es14.6,/, &
     " Please correct the input values for this BC!")
  !
end subroutine mach_total_pres
!
!###############################################################################
!
subroutine mach_total_dens(r,r0,mach,mach_is_set,n)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp),    intent(inout) :: r,r0
  real(wp),    intent(inout) :: mach
  logical(lk), intent(inout) :: mach_is_set
  integer,        intent(in) :: n
  !
  !.. Local Scalars ..
  real(wp) :: rmach,total_cond
  real(wp) :: two_over_gm1,gm1
  real(wp) :: one_over_gm1,rratio
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "mach_total_dens"
  !
continue
  !
  gm1 = gam - one
  !
  if (r > zero .and. r0 > zero) then
    !
    two_over_gm1 = two / gm1
    rratio = (r0/r)**gm1
    !
    rmach = sqrt( two_over_gm1 * (rratio - one) )
    !
    if (mach_is_set) then
      if (abs(mach-rmach) > eps6) then
        write (error_message,1) n,mach,rmach
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
    else
      mach = rmach
      mach_is_set = true
    end if
    !
  else
    !
    if (mach_is_set) then
      one_over_gm1 = one / gm1
      total_cond = (one + half*gm1*mach*mach)**one_over_gm1
      if (r > zero) then
        r0 = r * total_cond
      else if (r0 > zero) then
        r = r0 / total_cond
      end if
    end if
    !
  end if
  !
  ! Format Statements
  !
  1 format ( &
     " Error with BC # ",i0," !",/, &
     " Mach number specified for boundary condition = ",es14.6,/, &
     " Mach number computed from the stagnation density = ",es14.6,/, &
     " Please correct the input values for this BC!")
  !
end subroutine mach_total_dens
!
!###############################################################################
!
subroutine igl(rho,pres,temp,success,n)
  !
  !.. Use Statements ..
  use ovar, only : rgasref
  !
  !.. Formal Arguments ..
  real(wp),    intent(inout) :: rho
  real(wp),    intent(inout) :: pres
  real(wp),    intent(inout) :: temp
  logical(lk), intent(inout) :: success
  integer,        intent(in) :: n
  !
  !.. Local Scalars ..
  real(wp) :: pres_check
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "igl"
  !
continue
  !
  success = fals
  !
  if (count([rho,pres,temp] > zero) >= 2) then
    !
    if ( rho > zero .and. pres > zero .and. temp > zero ) then
      !
      pres_check = rho * rgasref * temp
      !
      if (abs(pres_check-pres) > eps6) then
        write (error_message,1) n,pres,pres_check
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      !
    else if ( rho < zero .and. pres > zero .and. temp > zero ) then
      rho = pres / (rgasref * temp)
    else if ( pres < zero .and. rho > zero .and. temp > zero ) then
      pres = rho * rgasref * temp
    else if ( temp < zero .and. rho > zero .and. pres > zero ) then
      temp = pres / (rgasref * rho)
    end if
    !
    success = true
    !
  end if
  !
  ! Format Statements
  !
  1 format (" Error with BC # ",i0," !",/, &
            " Pressure specified for boundary condition = ",es14.6,/, &
            " Pressure computed from other specified conditions = ",es14.6,/, &
            " Please correct the flow conditions for this BC in the input",/, &
            " file or set one of density, pressure, and temperature to a",/, &
            " negative value so it is explicitly computed.")
  !
end subroutine igl
!
!###############################################################################
!
include "Functions/vortex.f90"
include "Functions/taylor_green.f90"
include "Functions/usp2v.f90"
include "Functions/vsp2u.f90"
include "Functions/viscosity.f90"
include "Functions/temperature.f90"
include "Functions/smc000_interior_wall_radius.f90"
include "Functions/arn2_interior_wall_radius.f90"
!
!###############################################################################
!
end module initialization_mod
