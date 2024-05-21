!> Initialization for Sod's shock tube problem.
!> Intended for domain on x=[-0.5,0.5] with initial jump at x=0.0.
!>   Left  state is reference freestream state.
!>   Left  state is:  P/Pref=1.00, rho/rhoref=1.000.
!>   Right state is:  P/Pref=0.10, rho/rhoref=0.125.
!>   Velocity is zero for both states.
pure function sod_init(xyz) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : aref,uref,Lref
  use ovar, only : gam,rhoref,tref,pref
  use ovar, only : rgasref
  use ovar, only : sod_test_right_vector
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xyz
  !
  !.. Function Return Value ..
  real(wp), dimension(1:2+size(xyz)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: nc,mb,me,ne,nx
  real(wp) :: gm1,rhoa,rhoa2
  real(wp) :: rho,pres
  real(wp) :: dist_right
  !
  !.. Local Arrays ..
  real(wp), dimension(1:3) :: vel
  !
continue
  !
  nc = 1
  mb = nc + 1
  me = nc + size(xyz)
  ne = me + 1
  !
  ! Some needed constants
  !
  gm1 = gam - one
  rhoa = rhoref * aref
  rhoa2 = rhoa * aref
  !
  !
  ! ------------------------------
  ! Flow properties on left state.
  ! ------------------------------
  !
  ! Compute initial velocities
  !
  vel(1) = zero
  vel(2) = zero
  vel(3) = zero
  !
  ! Compute initial pressure
  !
  pres = pref
  !
  ! Compute density assuming constant temperature
  !
  rho = pres / (rgasref*tref)
  !
  !
  ! -------------------------------
  ! Flow properties on right state.
  ! -------------------------------
  !
  ! Distance to the right of the diaphram by projecting xyz
  ! onto the diaphram right side normal vector.
  ! This method supports grids not aligned with x-direction.
  nx = size(xyz,dim=1)  ! grid dim: 1d,2d,or 3d
  dist_right = sum( xyz(1:nx) * sod_test_right_vector(1:nx) )
  !
  if ( dist_right > 0.0_wp ) then
    pres = pres * 0.10_wp
    rho  = rho  * 0.125_wp
  end if
  !
  ! Evaluate the dimensional initial conditions using conserved variables
  !
  return_value(nc)    = rho
  return_value(mb:me) = rho*vel(1:size(xyz))
  return_value(ne)    = pres/gm1 + half*rho*selfdot(vel)
  !
  ! Non-dimensionalize the initial conditions
  !
  return_value(nc)    = return_value(nc)    / rhoref
  return_value(mb:me) = return_value(mb:me) / rhoa
  return_value(ne)    = return_value(ne)    / rhoa2
  !
end function sod_init
