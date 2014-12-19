pure function normal_invs_flx_cv(uu,rnrm) result(return_value)
  !
  use geovar, only : nr
  !
  real(wp), dimension(1:nq), intent(in) :: uu
  real(wp), dimension(1:nr),    intent(in) :: rnrm
  !
  real(wp), dimension(1:nq) :: return_value
  !
  integer :: n
  real(wp), dimension(1:nr) :: unrm
  !
continue
  !
  unrm(:) = unit_vector( rnrm(:) )
  !
  return_value = zero
  !
  do n = 1,nr
    return_value = return_value + unrm(n) * invs_flx_cv(uu,n)
  end do
  !
end function normal_invs_flx_cv
!
!###############################################################################
!
pure function invs_flx_cv(uu,dir) result(return_value) 
  !
  use ovar, only : gam
  !
  integer,                   intent(in) :: dir
  real(wp), dimension(1:nq), intent(in) :: uu
  !
  real(wp), dimension(1:nq) :: return_value
  !
  integer :: meq
  real(wp) :: pres
  !
continue
  !
  meq = nec + dir
  !
  pres = (gam-one)*(uu(nee) - half*selfdot(uu(nmb:nme))/uu(nec))
  !
  return_value(nec) = uu(meq)
  return_value(nmb:nme) = uu(meq)*uu(nmb:nme)/uu(nec)
  return_value(meq) = return_value(meq) + pres
  return_value(nee) = uu(meq)/uu(nec) * (uu(nee) + pres)
  !
end function invs_flx_cv
!
!###############################################################################
!
pure function invs_flx_pv(vv,dir) result(return_value) 
  !
  use ovar, only : gam
  !
  integer,                   intent(in) :: dir
  real(wp), dimension(1:nq), intent(in) :: vv
  !
  real(wp), dimension(1:nq) :: return_value
  !
  integer :: meq
  real(wp) :: ener
  !
continue
  !
  meq = nec + dir
  !
  ener = vv(nee)/(gam-one) + half*vv(nec)*selfdot(vv(nmb:nme))
  !
  return_value(nec) = vv(nec)*vv(meq)
  return_value(nmb:nme) = vv(nec)*vv(nmb:nme) * vv(meq)
  return_value(meq) = return_value(meq) + vv(nee)
  return_value(nee) = (ener + vv(nee)) * vv(meq)
  !
end function invs_flx_pv
!
!###############################################################################
!
pure function invs_flx_jacobian_cv(uu,dir) result(return_value)
  !
  use geovar, only : nr
  use ovar, only : gam
  !
  integer,                   intent(in) :: dir
  real(wp), dimension(1:nq), intent(in) :: uu
  !
  real(wp), dimension(1:nq,1:nq) :: return_value
  !
  integer :: meq,n
  real(wp) :: a,a2,q2,etot,gm1
  !
  real(wp), dimension(1:nr) :: vel
  !
continue
  !
  return_value = zero
  !
  meq = nec + dir
  !
  gm1 = gam - one
  !
  etot  = uu(nee)/uu(nec)
  !
  vel = uu(nmb:nme)/uu(nec)
  q2  = selfdot(vel)
  a   = vel(dir) ! Velocity in the flux direction, NOT speed of sound
  a2  = a*a
  !
  ! return_value(column,row)
  !
  ! DENSITY ROW
  return_value(meq,nec) = one
  !
  ! ENERGY COLUMN
  return_value(nee,meq) = gm1
  return_value(nee,nee) = gam*a
  !
  ! INITIALIZE MOMENTUM PARTS THAT ARE NONZERO IGNORING DIAGONAL CONTRIBUTION
  return_value(meq,nmb:nme) = vel
  return_value(nmb:nme,meq) = -gm1*vel
  !
  ! DIAGONAL CONTRIBUTION TO MOMENTUM
  do n = nmb,nme
    return_value(n,n) = return_value(n,n) + a
  end do
  ! DIAGONAL CONTRIBUTION TO MOMENTUM FOR CURRENT DIRECTION IS ACTUALLY AN
  ! ADDICTIONAL 2*a SO ADD ANOTHER a TO THIS TERM
  return_value(meq,meq) = return_value(meq,meq) + a
  !
  ! DENSITY COLUMN, ENERGY ROW
  return_value(nec,nee) = (gm1*q2 - gam*etot) * a
  !
  ! DENSITY COLUMN, MOMENTUM ROWS
  return_value(nec,nmb:nme) = -a * vel
  ! ADD ADDITIONAL TERM TO MOMENTUM ROW FOR CURRENT DIRECTION
  return_value(nec,meq) = return_value(nec,meq) + half*gm1*q2
  !
  ! ENERGY ROW, MOMENTUM COLUMNS
  return_value(nmb:nme,nee) = -a * gm1 * vel
  ! ADD ADDITIONAL TERM TO MOMENTUM COLUMN FOR CURRENT DIRECTION
  return_value(meq,nee) = return_value(meq,nee) + gam*etot - half*gm1*q2
  !
end function invs_flx_jacobian_cv
!
!###############################################################################
!
pure function invs_wall_flx_cv(unrm,uu) result(return_value) 
  !
  use geovar, only : nr
  use ovar, only : gam
  !
  real(wp), dimension(1:nr),    intent(in) :: unrm
  real(wp), dimension(1:nq), intent(in) :: uu
  !
  real(wp), dimension(1:nq) :: return_value
  !
  real(wp) :: pres
  !
continue
  !
  ! Initialize all fluxes to zero since only the pressure
  ! part of the momentum fluxes is all that remains of the
  ! interface fluxes at a wall boundary.
  !
  return_value(:) = zero
  !
  ! Get the interior pressure at the wall
  ! NOTE: The interior velocity/momentum at the wall is NOT necessarily zero
  !
  pres = (gam-one)*(uu(nee) - half*selfdot(uu(nmb:nme))/uu(nec))
  !
  ! All momentum fluxes are just the pressure normal to the wall
  !
  return_value(nmb:nme) = unrm * pres
  !
end function invs_wall_flx_cv
!
!###############################################################################
!
!pure function f_invs_flx_cv(uu) result(return_value) 
!  !
!  use ovar, only : gam
!  !
!  real(wp), dimension(1:nq), intent(in) :: uu
!  !
!  real(wp), dimension(1:nq) :: return_value
!  !
!  real(wp) :: pres
!  !
!continue
!  !
!  pres = (gam-one)*(uu(nq) - half*selfdot(uu(nmb:nme))/uu(1))
!  !
!  return_value(1) = uu(2)
!  return_value(nmb:nme) = uu(2)*uu(nmb:nme)/uu(1)
!  return_value(2) = return_value(2) + pres
!  return_value(nq) = uu(2)/uu(1) * (uu(nq) + pres)
!  !
!end function f_invs_flx_cv
!!
!!###############################################################################
!!
!pure function g_invs_flx_cv(uu) result(return_value)
!  !
!  use ovar, only : gam
!  !
!  real(wp), dimension(1:nq), intent(in) :: uu
!  !
!  real(wp), dimension(1:nq) :: return_value
!  !
!  real(wp) :: pres
!  !
!continue
!  !
!  pres = (gam-one)*(uu(nq) - half*selfdot(uu(nmb:nme))/uu(1))
!  !
!  return_value(1) = uu(3)
!  return_value(nmb:nme) = uu(3)*uu(nmb:nme)/uu(1)
!  return_value(3) = return_value(3) + pres
!  return_value(nq) = uu(3)/uu(1) * (uu(nq) + pres)
!  !
!end function g_invs_flx_cv
!!
!!###############################################################################
!!
!pure function f_invs_flx_pv(vv) result(return_value) 
!  !
!  use ovar, only : gam
!  !
!  real(wp), dimension(1:nq), intent(in) :: vv
!  !
!  real(wp), dimension(1:nq) :: return_value
!  !
!  real(wp) :: ener
!  !
!continue
!  !
!  ener = vv(nq)/(gam-one) + half*vv(1)*selfdot(vv(nmb:nme))
!  !
!  return_value(1) = vv(1)*vv(2)
!  return_value(nmb:nme) = vv(1)*vv(nmb:nme) * vv(2)
!  return_value(2) = return_value(2) + vv(nq)
!  return_value(nq) = (ener + vv(nq)) * vv(2)
!  !
!end function f_invs_flx_pv
!!
!!###############################################################################
!!
!pure function g_invs_flx_pv(vv) result(return_value)
!  !
!  use ovar, only : gam
!  !
!  real(wp), dimension(1:nq), intent(in) :: vv
!  !
!  real(wp), dimension(1:nq) :: return_value
!  !
!  real(wp) :: ener
!  !
!continue
!  !
!  ener = vv(nq)/(gam-one) + half*vv(1)*selfdot(vv(nmb:nme))
!  !
!  return_value(1) = vv(1)*vv(3)
!  return_value(nmb:nme) = vv(1)*vv(nmb:nme) * vv(3)
!  return_value(3) = return_value(3) + vv(nq)
!  return_value(nq) = (ener + vv(nq)) * vv(3)
!  !
!end function g_invs_flx_pv
!!
!!###############################################################################
!!
!pure function f_invs_flx_jacobian_cv(uu) result(return_value)
!  !
!  use ovar, only : gam
!  !
!  real(wp), dimension(1:nq), intent(in) :: uu
!  !
!  real(wp), dimension(1:nq,1:nq) :: return_value
!  !
!  real(wp) :: r,u,u2,v,v2,e,gm1,gm3
!  !
!continue
!  !
!  r  = uu(1)
!  u  = uu(2)/r
!  v  = uu(3)/r
!  e  = uu(4)/r
!  !
!  u2 = u*u
!  v2 = v*v
!  !
!  gm1 = gam - one
!  gm3 = gam - three
!  !
!  return_value(1,1) = zero
!  return_value(2,1) = one
!  return_value(3,1) = zero
!  return_value(4,1) = zero
!  !
!  return_value(1,2) = half * (gm1*v2 + gm3*u2)
!  return_value(2,2) = -gm3*u
!  return_value(3,2) = -gm1*v
!  return_value(4,2) = gm1
!  !
!  return_value(1,3) = -u*v
!  return_value(2,3) = v
!  return_value(3,3) = u
!  return_value(4,3) = zero
!  !
!  return_value(1,4) = gm1*u*(u2+v2) - gam*u*e
!  return_value(2,4) = gam*e - half*gm1*(three*u2+v2)
!  return_value(3,4) = -gm1*u*v
!  return_value(4,4) = gam*u
!  !
!end function f_invs_flx_jacobian_cv
!!
!!###############################################################################
!!
!pure function g_invs_flx_jacobian_cv(uu) result(return_value)
!  !
!  use ovar, only : gam
!  !
!  real(wp), dimension(1:nq), intent(in) :: uu
!  !
!  real(wp), dimension(1:nq,1:nq) :: return_value
!  !
!  real(wp) :: r,u,u2,v,v2,e,gm1,gm3
!  !
!continue
!  !
!  gm1 = gam - one
!  gm3 = gam - three
!  !
!  r  = uu(1)
!  u  = uu(2)/r
!  v  = uu(3)/r
!  e  = uu(4)/r
!  !
!  u2 = u*u
!  v2 = v*v
!  !
!  return_value(1,1) = zero
!  return_value(2,1) = zero
!  return_value(3,1) = one
!  return_value(4,1) = zero
!  !
!  return_value(1,2) = -u*v
!  return_value(2,2) = v
!  return_value(3,2) = u
!  return_value(4,2) = zero
!  !
!  return_value(1,3) = half * (gm1*u2 + gm3*v2)
!  return_value(2,3) = -gm1*u
!  return_value(3,3) = -gm3*v
!  return_value(4,3) = gm1
!  !
!  return_value(1,4) = gm1*v*(u2+v2) - gam*v*e
!  return_value(2,4) = -gm1*u*v
!  return_value(3,4) = gam*e - half*gm1*(u2+three*v2)
!  return_value(4,4) = gam*v
!  !
!end function g_invs_flx_jacobian_cv
!
!###############################################################################
!
!pure function f_invs_flx_jacobian_pv(vv) result(return_value)
!  !
!  use ovar, only : gam
!  !
!  real(wp), dimension(1:nq), intent(in) :: vv
!  !
!  real(wp), dimension(1:nq,1:nq) :: return_value
!  !
!  real(wp) :: r,u,u2,v,v2,e,gm1,gm3
!  !
!continue
!  !
!  gm1 = gam - one
!  gm3 = gam - three
!  !
!  r  = vv(1)
!  u  = vv(2)
!  v  = vv(3)
!  e  = (vv(4)/gm1 + half*r*(vv(2)*vv(2)+vv(3)*vv(3))) / r
!  !
!  u2 = u*u
!  v2 = v*v
!  !
!  return_value(1,1) = zero
!  return_value(1,2) = one
!  return_value(1,3) = zero
!  return_value(1,4) = zero
!  !
!  return_value(2,1) = half * (gm1*v2 + gm3*u2)
!  return_value(2,2) = -gm3*u
!  return_value(2,3) = -gm1*v
!  return_value(2,4) = gm1
!  !
!  return_value(3,1) = -u*v
!  return_value(3,2) = v
!  return_value(3,3) = u
!  return_value(3,4) = zero
!  !
!  return_value(4,1) = gm1*u*(u2+v2) - gam*u*e
!  return_value(4,2) = gam*e - half*gm1*(three*u2+v2)
!  return_value(4,3) = -gm1*u*v
!  return_value(4,4) = gam*u
!  !
!end function f_invs_flx_jacobian_pv
!!
!!###############################################################################
!!
!pure function g_invs_flx_jacobian_pv(vv) result(return_value)
!  !
!  use ovar, only : gam
!  !
!  real(wp), dimension(1:nq), intent(in) :: vv
!  !
!  real(wp), dimension(1:nq,1:nq) :: return_value
!  !
!  real(wp) :: r,u,u2,v,v2,e,gm1,gm3
!  !
!continue
!  !
!  gm1 = gam - one
!  gm3 = gam - three
!  !
!  r  = vv(1)
!  u  = vv(2)
!  v  = vv(3)
!  e  = (vv(4)/gm1 + half*r*(vv(2)*vv(2)+vv(3)*vv(3))) / r
!  !
!  u2 = u*u
!  v2 = v*v
!  !
!  return_value(1,1) = zero
!  return_value(1,2) = zero
!  return_value(1,3) = one
!  return_value(1,4) = zero
!  !
!  return_value(2,1) = -u*v
!  return_value(2,2) = v
!  return_value(2,3) = u
!  return_value(2,4) = zero
!  !
!  return_value(3,1) = half * (gm1*u2 + gm3*v2)
!  return_value(3,2) = -gm1*u
!  return_value(3,3) = -gm3*v
!  return_value(3,4) = gm1
!  !
!  return_value(4,1) = gm1*v*(u2+v2) - gam*v*e
!  return_value(4,2) = -gm1*u*v
!  return_value(4,3) = gam*e - half*gm1*(u2+three*v2)
!  return_value(4,4) = gam*v
!  !
!end function g_invs_flx_jacobian_pv
!!
!!###############################################################################
!!
