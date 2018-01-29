pure function taylor_green_initialization(xyz) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : aref,uref,Lref
  use ovar, only : gam,rhoref,tref,pref
  use ovar, only : rgasref
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xyz
  !
  !.. Function Return Value ..
  real(wp), dimension(1:2+size(xyz)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: nc,mb,me,ne
  real(wp) :: x,y,z,gm1,rhoa,rhoa2
  real(wp) :: rho,pres
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
  x = xyz(1) / Lref
  y = xyz(2) / Lref
  z = xyz(3) / Lref
  !
  ! Compute initial velocities
  !
  vel(1) =  uref * sin(x) * cos(y) * cos(z)
  vel(2) = -uref * cos(x) * sin(y) * cos(z)
  vel(3) = zero
  !
  ! Compute initial pressure
  !
  pres = pref + rhoref*uref*uref/sixteen * (cos(two*x) + cos(two*y)) &
                                         * (cos(two*z) + two)
  !
  ! Compute density assuming constant temperature
  !
  rho = pres / (rgasref*tref)
  !
  ! Evaluate the dimensional initial conditions using conserved variables
  !
  return_value(nc)    = rho
  return_value(mb:me) = rho*vel
  return_value(ne)    = pres/gm1 + half*rho*selfdot(vel)
  !
  ! Non-dimensionalize the initial conditions
  !
  return_value(nc)    = return_value(nc)    / rhoref
  return_value(mb:me) = return_value(mb:me) / rhoa
  return_value(ne)    = return_value(ne)    / rhoa2
  !
end function taylor_green_initialization
