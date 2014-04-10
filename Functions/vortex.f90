pure function vortex2D(x,y) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : vxref,vyref,aref,uref,machref
  use ovar, only : gam,rhoref,tref,pref,itestcase
  !
  !.. Formal Arguments ..
  real(vbp), intent(in) :: x
  real(vbp), intent(in) :: y
  !
  !.. Function Return Value ..
  real(vbp), dimension(1:4) :: return_value
  !
  !.. Local Scalars ..
  real(vbp) :: gm1,rhoa,rhoa2
  real(vbp) :: c1,c2,c3,func,du,dv,dT,vx,vy,temp,rho,pres
  real(vbp) :: bigR,bigGam
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
  bigR   = merge(vortex_bigR,  density_bigR,  any(itestcase==Vortex_Problems))
  bigGam = merge(vortex_bigGam,density_bigGam,any(itestcase==Vortex_Problems))
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
  vx   = vxref + du*deltafun(any(itestcase==Vortex_Problems))
  vy   = vyref + dv*deltafun(any(itestcase==Vortex_Problems))
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
!pure function vortex3D(xyz) result(return_value)
!end function vortex3D
