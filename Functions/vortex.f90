pure function vortex_solution(xyz) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : velref,aref,uref,machref
  use ovar, only : gam,rhoref,tref,pref,itestcase
  use ovar, only : vortex_bigR,vortex_bigGam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: xyz
  !
  !.. Function Return Value ..
  real(wp), dimension(1:2+size(xyz)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: nc,mb,me,ne
  real(wp) :: x,y,gm1,rhoa,rhoa2
  real(wp) :: c1,c2,c3,func,du,dv,dT,vx,vy,temp,rho,pres
  real(wp) :: bigR,bigGam
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
  x = xyz(1)
  y = xyz(2)
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
  vel(1) = velref(1) + du*deltafun(any(itestcase==Vortex_Problems))
  vel(2) = velref(2) + dv*deltafun(any(itestcase==Vortex_Problems))
  vel(3) = velref(3)
  temp   = tref  + dT
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
  return_value(nc)    = rho
  return_value(mb:me) = rho*vel
  return_value(ne)    = pres/gm1*c3 + half*rho*selfdot(vel)
  !
  ! Non-dimensionalize the initial conditions
  !
  return_value(nc)    = return_value(nc)    / rhoref
  return_value(mb:me) = return_value(mb:me) / rhoa
  return_value(ne)    = return_value(ne)    / rhoa2
  !
end function vortex_solution
!
!###############################################################################
!
pure subroutine vortex_position(delx,dely)
  !
  !.. Use Statements ..
  use ovar, only : itcur,velref,aref
  use ovar, only : Period_Timesteps,Period_Time
  use ovar, only : Grid_Length
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
  delx = mod(delx,Grid_Length(1))
  dely = mod(dely,Grid_Length(2))
  !
end subroutine vortex_position
!
!###############################################################################
!
pure function vortex_offset(xyz) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : itcur,velref,aref
  use ovar, only : Period_Timesteps,Period_Time
  use ovar, only : Grid_Length,Grid_Min
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: xyz
  !
  !.. Function Result ..
  real(wp), dimension(1:size(xyz,dim=1),1:size(xyz,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer :: nr,np,l,k
  real(wp) :: delt_since_lastperiod
  !
  !.. Local Arrays ..
  real(wp) :: del(1:size(xyz,dim=1))
  !
continue
  !
  nr = size(xyz,dim=1)
  np = size(xyz,dim=2)
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
  do l = 1,nr
    del(l) = deltafun(velref(l)>eps3) * (velref(l)/aref) * delt_since_lastperiod
  end do
  !
  ! Now translate these into a periodic domain
  ! with the grid lengths for each direction
  !
  do l = 1,nr
    del(l) = mod(del(l),Grid_Length(l))
  end do
  !
  do k = 1,np
    do l = 1,nr
      return_value(l,k) = xyz(l,k) - merge( del(l) , &
                                            del(l) - Grid_Length(l) , &
                                            xyz(l,k) - del(l) >= Grid_Min(l) )
    end do
  end do
  !
end function vortex_offset
!
!###############################################################################
!
!pure function vortex3D(xyz) result(return_value)
!end function vortex3D
!
!###############################################################################
!
