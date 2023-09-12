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
end module taylor_couette_mod
