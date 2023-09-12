module taylor_couette_mod
  !
  ! Module containing variables and routines
  ! specific to the Taylor-Couette problem.
  !
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nee,nq
  !
  implicit none
  !
  private
  !
  public :: initialize_taylor_couette_problem
  public :: laminar_taylor_couette_init
  public :: turbulent_taylor_couette_init
  !
  ! Module constants needed for the initialization functions
  !
  real(wp), parameter :: sig = 0.1_wp
  real(wp), parameter :: sec_per_min = 60.0_wp
  !
  real(wp), parameter :: default_Re = 4000.0_wp
  real(wp), parameter :: default_Mach = 0.1_wp
  real(wp), parameter :: default_wall_temperature = 300.0_wp ! Kelvin
  real(wp), parameter :: default_pressure = 101325.0_wp ! Pa
  !
  real(wp), save :: ri ! radius of inner cylinder
  real(wp), save :: ro ! radius of outer cylinder
  real(wp), save :: lc ! length of cylinders
  !
  real(wp), save :: ro2  ! ro*ro
  !
  real(wp), save :: pconst ! pressure constant term due to wall temperature
  real(wp), save :: pcoeff ! leading coefficient for pressure term
  !
  real(wp), save :: sigui ! sig*(tangential wall velocity)
  !
  real(wp), save :: bigA ! constant for laminar solution
  real(wp), save :: bigB ! constant for laminar solution
  !
contains
!
!###############################################################################
!
subroutine initialize_taylor_couette_problem(ref_length)
  !
  !.. Use Statements ..
  use geovar,     only : xyz_nodes
  use ovar,       only : rhoref,tref,aref,pref,tmin
  use ovar,       only : gam,rgasref,reyref,machref
  use ovar,       only : bc_input
  use thermo_mod, only : rgas_cpg
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: ref_length
  !
  !.. Local Scalars ..
  integer  :: n,bc_idx
  real(wp) :: y,z,r2,ri2,eta2,xmin,xmax,rmin,rmax
  real(wp) :: rotation_rate,omega,wall_temperature
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_taylor_couette_problem"
  !
continue
  !
  ! Find the index in bc_input for the inner cylinder boundary conditions
  !
  bc_idx = 0
  find_inner_wall_bc: &
  do n = 1,size(bc_input)
    if (trim(adjustl(uppercase(bc_input(n)%name))) == "INNERWALL") then
      bc_idx = n
      exit find_inner_wall_bc
    end if
  end do find_inner_wall_bc
  !
  if (bc_idx == 0) then
    write (error_message,101)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Find the min/max values of x to determine the length of the cylinders
  !
  xmin = minval(xyz_nodes(1,:))
  xmax = maxval(xyz_nodes(1,:))
  !
  ! Save the length of the cylinders to its module variable
  !
  lc = xmax - xmin
  !
  ! Go through xyz_nodes to determine the radii of the two cylinders by
  ! finding minimum and maximum radii for the grid.
  ! NOTE: We will actually find the minimum and maximum radii^2 since that is
  !       cheaper to compute in the search loop. We will squareroot these
  !       after the loop to get the actual radii.
  ! 
  rmin = huge(zero)
  rmax = -huge(zero)
  !
  do n = 1,size(xyz_nodes,dim=2)
    !
    y = xyz_nodes(2,n)
    z = xyz_nodes(3,n)
    !
    r2 = y*y + z*z
    !
    if (r2 < rmin) rmin = r2
    if (r2 > rmax) rmax = r2
    !
  end do
  !
  ! Save the inner and outer radii to their module variables
  !
  ri = sqrt(rmin)
  ro = sqrt(rmax)
  !
  ! Set the reference length for the Taylor-Couette problem to the difference
  ! between the radii of the two cylinders, d=ro-ri. This reference length is
  ! only used to compute time_ref after returning from this subroutine.
  !
  ref_length = ro - ri
  !
  ! Get the wall temperature from the boundary condition. If the wall
  ! temperature hasnt been provided, first use the reference temperature
  ! if it is > tref, otherwise use the default wall temperature value
  ! from the module header.
  !
  wall_temperature = bc_input(bc_idx)%wall_temp
  if (wall_temperature <= tmin) then
    if (tref > tmin) then
      wall_temperature = tref
    else
      wall_temperature = default_wall_temperature
    end if
  end if
  !
  ! For this problem, set the reference temperature to
  ! the wall temperature for the inner rotating cylinder
  !
  tref = wall_temperature
  !
  ! Compute gamma, the reference gas constant, and the reference speed of sound
  !
  if (gam <= zero) gam = 1.4_wp
  if (rgasref <= zero) rgasref = rgas_cpg
  !
  aref = sqrt( max( gam*rgasref*tref , rndoff ) )
  !
  ! Set the reference pressure to 1 atm
  !
  pref = default_pressure
  !
  ! Compute the reference pressure using the ideal gas law
  !
  rhoref = pref/(rgasref*tref)
  !
  ! Check the input reyref for the Reynolds number to use for this problem.
  ! If reyref hasnt been specified (reyref < zero), use the default Reynolds
  ! number value from the module header (also set reyref to this value).
  ! Setting reyref will adjust the leading viscosity coefficient to get the
  ! correct Reynolds number for the simulation.
  !
  if (reyref <= zero) then
    reyref = default_Re
  end if
  !
  ! Get the rotation speed from the boundary condition.
  ! NOTE: This should be in RPMs
  !
  rotation_rate = bc_input(bc_idx)%rotation_rate(1)
  !
  if (rotation_rate <= eps1) then
    !
    machref = default_Mach
    !
    rotation_rate = machref*(aref*sec_per_min)/(ri*two*pi)
    !
  else
    !
    machref = rotation_rate*(ri*two*pi)/(aref*sec_per_min)
    !
  end if
  !
  ! Save the wall temperature and rotation speed back into the boundary
  ! condition in case either of these were updated
  !
  bc_input(bc_idx)%wall_temp = wall_temperature
  bc_input(bc_idx)%rotation_rate(1) = rotation_rate
  !
  ! Save the inner cylinder rotation speed.
  ! NOTE: This converts it from rpm to 1/s and when multiplied
  !       by a radius becomes m/s
  !
  omega = rotation_rate*two*pi/sec_per_min
  !
  ! Compute the remaining module variables that
  ! will be used in the initialization functions
  !
  ri2 = ri*ri
  ro2 = ro*ro
  !
  eta2 = ri2/ro2
  !
  bigA = -omega*eta2/(one-eta2)
  bigB =  omega* ri2/(one-eta2)
  !
  pcoeff = half*rhoref*bigB*bigB/ro2
  !
  pconst = pref - pcoeff*(eta2 - four*log(ri) - one/eta2)
  !
  sigui = sig*omega*ri
  !
  ! Format Statements
  !
  101 format ( &
    " ERROR: No inner wall boundary condition found! When running the", &
    " Taylor-Couette flow problem, one of the boundary conditions must", &
    " have the name 'InnerWall'. This should correspond to the inner", &
    " cylinder wall boundary, and should ideally specify the rotation", &
    " speed and wall temperature for the boundary.")
  !
end subroutine initialize_taylor_couette_problem
!
!###############################################################################
!
subroutine get_tcf_rpm(aspd,mach,tcf_rpm,iw_bc,ow_bc)
  !
  !.. Use Statements ..
  use ovar, only : moving_ref_frame,bc_input_t
  !
  !.. Formal Arguments ..
  real(wp),            intent(in) :: aspd
  real(wp),         intent(inout) :: mach
  real(wp),         intent(inout) :: tcf_rpm
  type(bc_input_t), intent(inout) :: iw_bc
  type(bc_input_t), intent(inout) :: ow_bc
  !
  !.. Local Scalars ..
  real(wp) :: iw_abs_rpm,ow_abs_rpm
  real(wp) :: mrf_rpm
  !
continue
  !
  ! The y- and z-rotations are zero for the TCF problem
  !
  moving_ref_frame%rotation(2:3) = zero
  iw_bc%rotation_rate(2:3) = zero
  ow_bc%rotation_rate(2:3) = zero
  !
  ! Set the MRF translation to zero
  moving_ref_frame%translation(:) = zero
  !
  if ( moving_ref_frame%is_active() ) then
    !
    ! Get the MRF RPMs
    !
    mrf_rpm = moving_ref_frame%rotation(1)
    !
    ! For now, if a rotating reference frame is being used, use the default
    ! Mach number for this problem to compute the rotation rates for the
    ! inner and outer wall boundary conditions.
    !
    mach = default_Mach
    !
    tcf_rpm = mach*aspd/(ri*rpm_to_mps)
    !
    ! Set the BC rotations depending on value of
    ! the flag is_specified_in_mrf for each BC.
    !
    if (iw_bc%is_specified_in_mrf) then
      iw_bc%rotation_rate(1) = tcf_rpm - mrf_rpm
    else
      iw_bc%rotation_rate(1) = tcf_rpm
    end if
    !
    if (ow_bc%is_specified_in_mrf) then
      ow_bc%rotation_rate(1) = -mrf_rpm
    else
      ow_bc%rotation_rate(1) = zero
    end if
    !
  else
    !
    iw_bc%is_specified_in_mrf = fals
    ow_bc%is_specified_in_mrf = fals
    !
    iw_abs_rpm = iw_bc%rotation_rate(1)
    ow_abs_rpm = ow_bc%rotation_rate(1)
    !
    if (abs(iw_abs_rpm-ow_abs_rpm) >= one) then
      !
      tcf_rpm = abs(iw_abs_rpm - ow_abs_rpm)
      !
      mach = tcf_rpm*ri*rpm_to_mps/aspd
      !
    else
      !
      mach = default_Mach
      !
      tcf_rpm = mach*aspd/(ri*rpm_to_mps)
      !
    end if
    !
    iw_bc%rotation_rate(1) = tcf_rpm ! Reset the IW rotation to TCF
    ow_bc%rotation_rate(1) = zero    ! Reset the OW rotation to zero
    !
  end if
  !
end subroutine get_tcf_rpm
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function laminar_taylor_couette_init(xyz) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam,rhoref,aref
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:nq)
  !
  !.. Local Scalars ..
  real(wp) :: r,r2,y,z,theta
  real(wp) :: vt,pres,rho,vx,vy,vz
  real(wp) :: gm1,rhoa,rhoa2
  !
continue
  !
  gm1 = gam - one
  rhoa = rhoref * aref
  rhoa2 = rhoa * aref
  !
  y = xyz(2)
  z = xyz(3)
  !
  r2 = y*y + z*z
  r  = sqrt(r2)
  !
  theta = atan2(z,y)
  !
  vt = bigA*r + bigB/r ! r >= ri > zero
  !
  pres = pcoeff*(r2/ro2 - four*log(r) - ro2/r2)
  !
  ! Compute the primitive variable initial conditions
  !
  rho  = rhoref
  vx   = zero
  vy   = -vt * sin(theta)
  vz   =  vt * cos(theta)
  pres = pres + pconst
  !
  ! Convert to conserved variables
  !
  return_value(nec) = rho
  return_value(nmx) = rho*vx
  return_value(nmy) = rho*vy
  return_value(nmz) = rho*vz
  return_value(nee) = pres/gm1 + half*rho*(vx*vx + vy*vy + vz*vz)
  !
  ! Non-dimensionalize the initial conditions
  !
  return_value(nec) = return_value(nec) / rhoref
  return_value(nmx) = return_value(nmx) / rhoa
  return_value(nmy) = return_value(nmy) / rhoa
  return_value(nmz) = return_value(nmz) / rhoa
  return_value(nee) = return_value(nee) / rhoa2
  !
end function laminar_taylor_couette_init
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function turbulent_taylor_couette_init(xyz) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam,rhoref,aref
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:nq)
  !
  !.. Local Scalars ..
  real(wp) :: r,r2,x,y,z,a,l,theta
  real(wp) :: vr,vt,vt_lam,p_lam
  real(wp) :: pres,rho,vx,vy,vz
  real(wp) :: gm1,rhoa,rhoa2
  !
continue
  !
  gm1 = gam - one
  rhoa = rhoref * aref
  rhoa2 = rhoa * aref
  !
  x = xyz(1)
  y = xyz(2)
  z = xyz(3)
  !
  r2 = y*y + z*z
  r  = sqrt(r2)
  !
  theta = atan2(z,y)
  !
  ! Compute the laminar v-theta and pressure
  !
  vt_lam = bigA*r + bigB/r ! r >= ri > zero
  !
  p_lam = pcoeff*(r2/ro2 - four*log(r) - ro2/r2) + pconst
  !
  ! Compute the radial and theta velocities with additional terms designed
  ! to induce turbulence within the flow. The additional terms are based off
  ! the initial conditions for the Taylor-Green vortex problem.
  !
  a = pi * (r/ri - one)
  l = x / (ro - ri)
  !
  vr = sigui*cos(theta)*sin(a)*sin(l)
  vt = sigui*sin(theta)*sin(a)*sin(l) + vt_lam
  !
  ! Compute the pressure with a similar additional term
  !
  pres = half*rhoref*sigui*sigui*cos(two*theta)*sin(two*a)*sin(two*l)
  !
  ! Set the initial conditions for this point
  !
  rho  = rhoref
  vx   = zero
  vy   = vr*cos(theta) - vt*sin(theta)
  vz   = vr*sin(theta) + vt*cos(theta)
  pres = p_lam + pres
  !
  ! Convert to conserved variables
  !
  return_value(nec) = rho
  return_value(nmx) = rho*vx
  return_value(nmy) = rho*vy
  return_value(nmz) = rho*vz
  return_value(nee) = pres/gm1 + half*rho*(vx*vx + vy*vy + vz*vz)
  !
  ! Non-dimensionalize the initial conditions
  !
  return_value(nec) = return_value(nec) / rhoref
  return_value(nmx) = return_value(nmx) / rhoa
  return_value(nmy) = return_value(nmy) / rhoa
  return_value(nmz) = return_value(nmz) / rhoa
  return_value(nee) = return_value(nee) / rhoa2
  !
end function turbulent_taylor_couette_init
!
!###############################################################################
!
function freestream_taylor_couette_init(xyz) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam,bc_in,itestcase
  use ovar, only : moving_ref_frame
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:nq)
  !
  !.. Local Scalars ..
  integer  :: n,ibc
  real(wp) :: x_r,y_r,z_r
  real(wp) :: bc_rot,bc_max_rot
  real(wp) :: gm1,pres,rho,vx,vy,vz
  !
continue
  !
  gm1 = gam - one
  !
  ! Always initialize to freestream conditions
  !
  return_value(:) = bc_in(0)%cv(:)
  !
  if (.not. allocated(tcf_fs_cv)) then
    !
    tcf_fs_cv = bc_in(0)%cv(:) ! F2003 AUTO-REALLOC
    tcf_fs_pv = bc_in(0)%pv(:) ! F2003 AUTO-REALLOC
    tcf_fs_gp = bc_in(0)%gp(:) ! F2003 AUTO-REALLOC
    !
    rho = bc_in(0)%cv(nec)
    pres = bc_in(0)%gp(npres)
    !
    tcf_fs_pv(nec) = sqrt(gam*pres/rho) ! speed of sound
    tcf_fs_pv(nee) = pres/gm1 ! internal energy
    tcf_fs_gp(ngamm) = gam*pres/(gm1*rho) ! enthalpy
    !
  end if
  !
  ! If itestcase is negative, initialize the y and z velocities to either
  !   1. If using a MRF, the local MRF y and z velocities, otherwise
  !   2. The local y and z velocity using the maximum rotation rate of
  !      all the boundary conditions.
  !
  if (itestcase < 0) then
    !
    if ( moving_ref_frame%is_active() ) then
      !
      x_r = moving_ref_frame%rotation(1)
      y_r = moving_ref_frame%rotation(2)
      z_r = moving_ref_frame%rotation(3)
      !
    else
      !
      bc_max_rot = -huge(zero)
      do n = 0,ubound(bc_in,dim=1)
        bc_rot = norm2(bc_in(n)%rotation_rate(:))
        if (bc_rot > bc_max_rot) then
          bc_max_rot = bc_rot
          ibc = n
        end if
      end do
      !
      x_r = bc_in(ibc)%rotation_rate(1)
      y_r = bc_in(ibc)%rotation_rate(2)
      z_r = bc_in(ibc)%rotation_rate(3)
      !
    end if
    !
    if (.not. allocated(tcf_rotation)) then
      tcf_rotation = [x_r,y_r,z_r] ! F2003 AUTO-REALLOC
    end if
    !
    ! Get freestream density and pressure
    !
    rho = bc_in(0)%cv(nec)
    pres = bc_in(0)%gp(npres)
    !
    ! Get freestream velocity plus rotation
    !
    vx = bc_in(0)%pv(nmx) + (y_r*xyz(3) - z_r*xyz(2))
    vy = bc_in(0)%pv(nmy) + (z_r*xyz(1) - x_r*xyz(3))
    vz = bc_in(0)%pv(nmz) + (x_r*xyz(2) - y_r*xyz(1))
    !
    ! Convert to conserved variables
    !
    return_value(nec) = rho
    return_value(nmx) = rho*vx
    return_value(nmy) = rho*vy
    return_value(nmz) = rho*vz
    return_value(nee) = pres/gm1 + half*rho*(vx*vx + vy*vy + vz*vz)
    !
  end if
  !
  if (.not. allocated(tcf_rotation)) then
    tcf_rotation = [zero,zero,zero] ! F2003 AUTO-REALLOC
  end if
  !
end function freestream_taylor_couette_init
!
!###############################################################################
!
subroutine remove_freestream_tcf_init(nvar,xyz,var)
  !
  !.. Use Statements ..
  use ovar, only : gam,aref,rhoref,tref
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nvar
  real(wp),    intent(in) :: xyz(:,:)
  real(wp), intent(inout) :: var(:)
  !
  !.. Local Scalars ..
  integer  :: n,npts
  real(wp) :: x_r,y_r,z_r,varref
  real(wp) :: gm1,pres,rho,vx,vy,vz,vr
  !
continue
  !
 !return
  !
  ! The input argument nvar selects which variable we are going to subtract
  ! out the freestream TCF initial conditions from. The following is the list
  ! of variables for nvar:
  !
  !   value of nvar :  variable
  !   ----------------------------------------------
  !   nec         : density
  !   nmx         : x-momentum
  !   nmy         : y-momentum
  !   nmz         : z-momentum
  !   nee         : total energy per unit volume
  !   nq+npres    : static pressure
  !   nq+ntemp    : static temperature
  !   nq+ngamm    : static enthalpy
  !   -nec        : total density
  !   -nmx        : x-velocity
  !   -nmy        : y-velocity
  !   -nmz        : z-velocity
  !   -nee        : Mach number
  !   -(nq+npres) : total pressure
  !   -(nq+ntemp) : total temperature
  !   -(nq+ngamm) : total enthalpy
  !
  npts = size(var)
  !
  gm1 = gam - one
  !
  x_r = tcf_rotation(1)
  y_r = tcf_rotation(2)
  z_r = tcf_rotation(3)
  !
  if (nvar == nec) then
    ! Static density
    varref = rhoref
    do n = 1,npts
      var(n) = var(n) - varref*tcf_fs_cv(nec)
    end do
    !
  else if (nvar == nmx) then
    ! x-momentum
    varref = rhoref*aref
    do n = 1,npts
      vr = y_r*xyz(n,3) - z_r*xyz(n,2)
      var(n) = var(n) - varref*(tcf_fs_cv(nmx) + tcf_fs_cv(nec)*vr)
    end do
    !
  else if (nvar == nmy) then
    ! y-momentum
    varref = rhoref*aref
    do n = 1,npts
      vr = z_r*xyz(n,1) - x_r*xyz(n,3)
      var(n) = var(n) - varref*(tcf_fs_cv(nmy) + tcf_fs_cv(nec)*vr)
    end do
    !
  else if (nvar == nmz) then
    ! z-momentum
    varref = rhoref*aref
    do n = 1,npts
      vr = x_r*xyz(n,2) - y_r*xyz(n,1)
      var(n) = var(n) - varref*(tcf_fs_cv(nmz) + tcf_fs_cv(nec)*vr)
    end do
    !
  else if (nvar == nee) then
    ! Total energy per unit volume
    varref = rhoref*aref*aref
    do n = 1,npts
      rho = tcf_fs_cv(nec)
      pres = tcf_fs_gp(npres)
      vx = tcf_fs_pv(nmx) + (y_r*xyz(n,3) - z_r*xyz(n,2))
      vy = tcf_fs_pv(nmy) + (z_r*xyz(n,1) - x_r*xyz(n,3))
      vz = tcf_fs_pv(nmz) + (x_r*xyz(n,2) - y_r*xyz(n,1))
      var(n) = var(n) - varref*(pres/gm1 + half*rho*(vx*vx + vy*vy + vz*vz))
    end do
    !
  else if (nvar == nq+npres) then
    ! Static pressure
    varref = rhoref*aref*aref
    do n = 1,npts
      var(n) = var(n) - varref*tcf_fs_gp(npres)
    end do
    !
  else if (nvar == nq+ntemp) then
    ! Static temperature
    varref = tref
    do n = 1,npts
      var(n) = var(n) - varref*tcf_fs_gp(ntemp)
    end do
    !
  else if (nvar == nq+ngamm) then
    ! Static enthalpy
    varref = aref*aref
    do n = 1,npts
      var(n) = var(n) - varref*tcf_fs_gp(ngamm)
    end do
    !
  else if (nvar == -nec) then
    ! Total density
    varref = rhoref
    do n = 1,npts
      vx = tcf_fs_pv(nmx) + (y_r*xyz(n,3) - z_r*xyz(n,2))
      vy = tcf_fs_pv(nmy) + (z_r*xyz(n,1) - x_r*xyz(n,3))
      vz = tcf_fs_pv(nmz) + (x_r*xyz(n,2) - y_r*xyz(n,1))
      var(n) = var(n) - varref*(tcf_fs_cv(nec) + half*(vx*vx+vy*vy+vz*vz))
    end do
    !
  else if (nvar == -nmx) then
    ! x-velocity
    varref = aref
    do n = 1,npts
      vr = y_r*xyz(n,3) - z_r*xyz(n,2)
      var(n) = var(n) - varref*(tcf_fs_pv(nmx) + vr)
    end do
    !
  else if (nvar == -nmy) then
    ! y-velocity
    varref = aref
    do n = 1,npts
      vr = z_r*xyz(n,1) - x_r*xyz(n,3)
      var(n) = var(n) - varref*(tcf_fs_pv(nmy) + vr)
    end do
    !
  else if (nvar == -nmz) then
    ! z-velocity
    varref = aref
    do n = 1,npts
      vr = x_r*xyz(n,2) - y_r*xyz(n,1)
      var(n) = var(n) - varref*(tcf_fs_pv(nmz) + vr)
    end do
    !
  else if (nvar == -nee) then
    ! Mach number
    ! NOTE: Freestream speed of sound was stored in tcf_fs_pv(nec)
    do n = 1,npts
      vx = tcf_fs_pv(nmx) + (y_r*xyz(n,3) - z_r*xyz(n,2))
      vy = tcf_fs_pv(nmy) + (z_r*xyz(n,1) - x_r*xyz(n,3))
      vz = tcf_fs_pv(nmz) + (x_r*xyz(n,2) - y_r*xyz(n,1))
      var(n) = var(n) - (sqrt(vx*vx + vy*vy + vz*vz)/tcf_fs_pv(nec))
    end do
    !
  else if (nvar == -nq-npres) then
    ! Total pressure
    varref = rhoref*aref*aref
    do n = 1,npts
      vx = tcf_fs_pv(nmx) + (y_r*xyz(n,3) - z_r*xyz(n,2))
      vy = tcf_fs_pv(nmy) + (z_r*xyz(n,1) - x_r*xyz(n,3))
      vz = tcf_fs_pv(nmz) + (x_r*xyz(n,2) - y_r*xyz(n,1))
      var(n) = var(n) - varref*(tcf_fs_gp(npres) + half*(vx*vx+vy*vy+vz*vz))
    end do
    !
  else if (nvar == -nq-ntemp) then
    ! Total temperature
    varref = tref
    do n = 1,npts
      vx = tcf_fs_pv(nmx) + (y_r*xyz(n,3) - z_r*xyz(n,2))
      vy = tcf_fs_pv(nmy) + (z_r*xyz(n,1) - x_r*xyz(n,3))
      vz = tcf_fs_pv(nmz) + (x_r*xyz(n,2) - y_r*xyz(n,1))
      var(n) = var(n) - varref*(tcf_fs_gp(ntemp) + half*(vx*vx+vy*vy+vz*vz))
    end do
    !
  else if (nvar == -nq-ngamm) then
    ! Total enthalpy
    varref = aref*aref
    do n = 1,npts
      vx = tcf_fs_pv(nmx) + (y_r*xyz(n,3) - z_r*xyz(n,2))
      vy = tcf_fs_pv(nmy) + (z_r*xyz(n,1) - x_r*xyz(n,3))
      vz = tcf_fs_pv(nmz) + (x_r*xyz(n,2) - y_r*xyz(n,1))
      var(n) = var(n) - varref*(tcf_fs_gp(ngamm) + half*(vx*vx+vy*vy+vz*vz))
    end do
    !
  else
    !
    ! Static density
    varref = rhoref
    do n = 1,npts
      var(n) = var(n) - varref*tcf_fs_cv(nec)
    end do
    !
  end if
  !
end subroutine remove_freestream_tcf_init
!
!###############################################################################
!
end module taylor_couette_mod
