pure function vortex2D(x,y) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : velref,aref,uref,machref
  use ovar, only : gam,rhoref,tref,pref,itestcase
  use ovar, only : vortex_bigR,vortex_bigGam
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: x
  real(wp), intent(in) :: y
  !
  !.. Function Return Value ..
  real(wp), dimension(1:4) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: gm1,rhoa,rhoa2
  real(wp) :: c1,c2,c3,func,du,dv,dT,vx,vy,temp,rho,pres
  real(wp) :: bigR,bigGam
  !
continue
  !
  ! Some needed constants
  !
  gm1 = gam - one
  rhoa = rhoref * aref
  rhoa2 = rhoa * aref
  !
  ! For density/entropy transport or Shu vortex -> bigR = 1;   bigGam = 5
  ! For vortex          transport               -> bigR = 1.5; bigGam = 13.5
  !
  bigR   = merge(vortex_bigR,  one, any(itestcase==Vortex_Problems))
  bigGam = merge(vortex_bigGam,five,any(itestcase==Vortex_Problems))
  !
  c2 = half*machref*machref*gm1
  !
  ! If doing density transport, make c3 = 0 to zero out pressure part of energy
  !
  c3 = one - deltafun(itestcase==Density_Transport)
  !
  func = (one - x*x - y*y)/(two*bigR*bigR)
  !
  c1 = bigGam*exp(func)/(two*pi)
  !
  du = -c1*y/bigR * uref
  dv =  c1*x/bigR * uref
  dT = -c2*c1*c1 * tref
  !
 !du = merge(du,zero,abs(du)>=eps6)
 !dv = merge(dv,zero,abs(dv)>=eps6)
  !
  ! Only add on the velocity pertubations if doing vortex transport
  ! i.e., the velocity is constant if doing density/entropy transport
  !
  vx   = velref(1) + du*deltafun(any(itestcase==Vortex_Problems))
  vy   = velref(2) + dv*deltafun(any(itestcase==Vortex_Problems))
  temp = tref  + dT
  !
  rho = rhoref * (temp/tref)**(one/gm1)
  pres = pref  * (temp/tref)**(gam/gm1)
  !
  ! If doing density/entropy transport, make the pressure
  ! constant everywhere and set it to the free-stream pressure 
  !
  pres = merge(pres,pref,any(itestcase==Vortex_Problems))
  !
  ! Evaluate the dimensional initial conditions using conserved variables
  !
  return_value(1) = rho
  return_value(2) = rho*vx
  return_value(3) = rho*vy
  return_value(4) = pres/gm1*c3 + half*rho*(vx*vx + vy*vy)
  !
  ! Non-dimensionalize the initial conditions
  !
  return_value(1) = return_value(1) / rhoref
  return_value(2) = return_value(2) / rhoa
  return_value(3) = return_value(3) / rhoa
  return_value(4) = return_value(4) / rhoa2
  !
end function vortex2D
!
!###############################################################################
!
pure subroutine vortex_position(delx,dely)
  !
  !.. Use Statements ..
  use ovar, only : itcur,velref,aref
  use ovar, only : Period_Timesteps,Period_Time
  use ovar, only : Grid_Length_X,Grid_Length_Y
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: delx
  real(wp), intent(out) :: dely
  !
  !.. Local Scalars ..
  real(wp) :: delt_since_lastperiod
  !
continue
  !
  ! Compute the time that has elapsed since the start of the current period
  !
  delt_since_lastperiod = real(mod(itcur,Period_Timesteps),wp) &
                        / real(Period_Timesteps,wp)
  delt_since_lastperiod = delt_since_lastperiod * Period_Time
  !
  ! Find the x and y distances the vortex has
  ! traveled since the start of the current period
  !
  delx = deltafun(velref(1)>eps3) * (velref(1)/aref) * delt_since_lastperiod
  dely = deltafun(velref(2)>eps3) * (velref(2)/aref) * delt_since_lastperiod
  !
  ! Now translate these into a periodic domain
  ! with the grid lengths for each direction
  !
  delx = mod(delx,Grid_Length_X)
  dely = mod(dely,Grid_Length_Y)
  !
end subroutine vortex_position
!
!###############################################################################
!
!pure function vortex3D(xyz) result(return_value)
!end function vortex3D
!
!###############################################################################
!
