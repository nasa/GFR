pure function normal_invs_flx_cv(uu,rnrm) result(return_value)
  !
  use geovar, only : nr
  !
  real(wp), dimension(1:nq), intent(in) :: uu
  real(wp), dimension(1:nr),    intent(in) :: rnrm
  !
  real(wp), dimension(1:nq) :: return_value
  !
  real(wp), dimension(1:nr)    :: unrm
  real(wp), dimension(1:nq) :: f_flx,g_flx
  !
continue
  !
  unrm(:) = unit_vector( rnrm(:) )
  !
  f_flx(:) = f_invs_flx_cv( uu(:) )
  g_flx(:) = g_invs_flx_cv( uu(:) )
  !
  return_value(:) = f_flx(:)*unrm(1) + g_flx(:)*unrm(2)
  !
end function normal_invs_flx_cv
!
!###############################################################################
!
pure function f_invs_flx_cv(uu) result(return_value) 
  !
  use ovar, only : gam
  !
  real(wp), dimension(1:nq), intent(in) :: uu
  !
  real(wp), dimension(1:nq) :: return_value
  !
  real(wp) :: pres
  !
continue
  !
  pres = (gam-one)*(uu(nq) - half*selfdot(uu(2:nq-1))/uu(1))
  !
  return_value(1) = uu(2)
  return_value(2:nq-1) = uu(2)*uu(2:nq-1)/uu(1)
  return_value(2) = return_value(2) + pres
  return_value(nq) = uu(2)/uu(1) * (uu(nq) + pres)
  !
end function f_invs_flx_cv
!
!###############################################################################
!
pure function g_invs_flx_cv(uu) result(return_value)
  !
  use ovar, only : gam
  !
  real(wp), dimension(1:nq), intent(in) :: uu
  !
  real(wp), dimension(1:nq) :: return_value
  !
  real(wp) :: pres
  !
continue
  !
  pres = (gam-one)*(uu(nq) - half*selfdot(uu(2:nq-1))/uu(1))
  !
  return_value(1) = uu(3)
  return_value(2:nq-1) = uu(3)*uu(2:nq-1)/uu(1)
  return_value(3) = return_value(3) + pres
  return_value(nq) = uu(3)/uu(1) * (uu(nq) + pres)
  !
end function g_invs_flx_cv
!
!###############################################################################
!
pure function f_invs_flx_pv(vv) result(return_value) 
  !
  use ovar, only : gam
  !
  real(wp), dimension(1:nq), intent(in) :: vv
  !
  real(wp), dimension(1:nq) :: return_value
  !
  real(wp) :: ener
  !
continue
  !
  ener = vv(nq)/(gam-one) + half*vv(1)*selfdot(vv(2:nq-1))
  !
  return_value(1) = vv(1)*vv(2)
  return_value(2:nq-1) = vv(1)*vv(2:nq-1) * vv(2)
  return_value(2) = return_value(2) + vv(nq)
  return_value(nq) = (ener + vv(nq)) * vv(2)
  !
end function f_invs_flx_pv
!
!###############################################################################
!
pure function g_invs_flx_pv(vv) result(return_value)
  !
  use ovar, only : gam
  !
  real(wp), dimension(1:nq), intent(in) :: vv
  !
  real(wp), dimension(1:nq) :: return_value
  !
  real(wp) :: ener
  !
continue
  !
  ener = vv(nq)/(gam-one) + half*vv(1)*selfdot(vv(2:nq-1))
  !
  return_value(1) = vv(1)*vv(3)
  return_value(2:nq-1) = vv(1)*vv(2:nq-1) * vv(3)
  return_value(3) = return_value(3) + vv(nq)
  return_value(nq) = (ener + vv(nq)) * vv(3)
  !
end function g_invs_flx_pv
!
!###############################################################################
!
pure function f_invs_flx_jacobian_cv(uu) result(return_value)
  !
  use ovar, only : gam
  !
  real(wp), dimension(1:nq), intent(in) :: uu
  !
  real(wp), dimension(1:nq,1:nq) :: return_value
  !
  real(wp) :: r,u,u2,v,v2,e,gm1,gm3
  !
continue
  !
  r  = uu(1)
  u  = uu(2)/r
  v  = uu(3)/r
  e  = uu(4)/r
  !
  u2 = u*u
  v2 = v*v
  !
  gm1 = gam - one
  gm3 = gam - three
  !
  return_value(1,1) = zero
  return_value(2,1) = one
  return_value(3,1) = zero
  return_value(4,1) = zero
  !
  return_value(1,2) = half * (gm1*v2 + gm3*u2)
  return_value(2,2) = -gm3*u
  return_value(3,2) = -gm1*v
  return_value(4,2) = gm1
  !
  return_value(1,3) = -u*v
  return_value(2,3) = v
  return_value(3,3) = u
  return_value(4,3) = zero
  !
  return_value(1,4) = gm1*u*(u2+v2) - gam*u*e
  return_value(2,4) = gam*e - half*gm1*(three*u2+v2)
  return_value(3,4) = -gm1*u*v
  return_value(4,4) = gam*u
  !
end function f_invs_flx_jacobian_cv
!
!###############################################################################
!
pure function g_invs_flx_jacobian_cv(uu) result(return_value)
  !
  use ovar, only : gam
  !
  real(wp), dimension(1:nq), intent(in) :: uu
  !
  real(wp), dimension(1:nq,1:nq) :: return_value
  !
  real(wp) :: r,u,u2,v,v2,e,gm1,gm3
  !
continue
  !
  gm1 = gam - one
  gm3 = gam - three
  !
  r  = uu(1)
  u  = uu(2)/r
  v  = uu(3)/r
  e  = uu(4)/r
  !
  u2 = u*u
  v2 = v*v
  !
  return_value(1,1) = zero
  return_value(2,1) = zero
  return_value(3,1) = one
  return_value(4,1) = zero
  !
  return_value(1,2) = -u*v
  return_value(2,2) = v
  return_value(3,2) = u
  return_value(4,2) = zero
  !
  return_value(1,3) = half * (gm1*u2 + gm3*v2)
  return_value(2,3) = -gm1*u
  return_value(3,3) = -gm3*v
  return_value(4,3) = gm1
  !
  return_value(1,4) = gm1*v*(u2+v2) - gam*v*e
  return_value(2,4) = -gm1*u*v
  return_value(3,4) = gam*e - half*gm1*(u2+three*v2)
  return_value(4,4) = gam*v
  !
end function g_invs_flx_jacobian_cv
!
!###############################################################################
!
pure function f_invs_flx_jacobian_pv(vv) result(return_value)
  !
  use ovar, only : gam
  !
  real(wp), dimension(1:nq), intent(in) :: vv
  !
  real(wp), dimension(1:nq,1:nq) :: return_value
  !
  real(wp) :: r,u,u2,v,v2,e,gm1,gm3
  !
continue
  !
  gm1 = gam - one
  gm3 = gam - three
  !
  r  = vv(1)
  u  = vv(2)
  v  = vv(3)
  e  = (vv(4)/gm1 + half*r*(vv(2)*vv(2)+vv(3)*vv(3))) / r
  !
  u2 = u*u
  v2 = v*v
  !
  return_value(1,1) = zero
  return_value(1,2) = one
  return_value(1,3) = zero
  return_value(1,4) = zero
  !
  return_value(2,1) = half * (gm1*v2 + gm3*u2)
  return_value(2,2) = -gm3*u
  return_value(2,3) = -gm1*v
  return_value(2,4) = gm1
  !
  return_value(3,1) = -u*v
  return_value(3,2) = v
  return_value(3,3) = u
  return_value(3,4) = zero
  !
  return_value(4,1) = gm1*u*(u2+v2) - gam*u*e
  return_value(4,2) = gam*e - half*gm1*(three*u2+v2)
  return_value(4,3) = -gm1*u*v
  return_value(4,4) = gam*u
  !
end function f_invs_flx_jacobian_pv
!
!###############################################################################
!
pure function g_invs_flx_jacobian_pv(vv) result(return_value)
  !
  use ovar, only : gam
  !
  real(wp), dimension(1:nq), intent(in) :: vv
  !
  real(wp), dimension(1:nq,1:nq) :: return_value
  !
  real(wp) :: r,u,u2,v,v2,e,gm1,gm3
  !
continue
  !
  gm1 = gam - one
  gm3 = gam - three
  !
  r  = vv(1)
  u  = vv(2)
  v  = vv(3)
  e  = (vv(4)/gm1 + half*r*(vv(2)*vv(2)+vv(3)*vv(3))) / r
  !
  u2 = u*u
  v2 = v*v
  !
  return_value(1,1) = zero
  return_value(1,2) = zero
  return_value(1,3) = one
  return_value(1,4) = zero
  !
  return_value(2,1) = -u*v
  return_value(2,2) = v
  return_value(2,3) = u
  return_value(2,4) = zero
  !
  return_value(3,1) = half * (gm1*u2 + gm3*v2)
  return_value(3,2) = -gm1*u
  return_value(3,3) = -gm3*v
  return_value(3,4) = gm1
  !
  return_value(4,1) = gm1*v*(u2+v2) - gam*v*e
  return_value(4,2) = -gm1*u*v
  return_value(4,3) = gam*e - half*gm1*(u2+three*v2)
  return_value(4,4) = gam*v
  !
end function g_invs_flx_jacobian_pv
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
  pres = (gam-one)*(uu(nq) - half*selfdot(uu(2:nq-1))/uu(1))
  !
  ! All momentum fluxes are just the pressure normal to the wall
  !
  return_value(2:nq-1) = unrm(1:nr)*pres
  !
end function invs_wall_flx_cv
!
!###############################################################################
!
