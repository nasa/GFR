pure function normal_invs_flx_cv(uu,rnrm) result(return_value)
  !
  use geovar, only : nr
  !
  real(wp), dimension(1:nq), intent(in) :: uu
  real(wp), dimension(1:nr), intent(in) :: rnrm
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
pure function invs_face_flx_cv(uu,rnrm,cell_geom,cell_order,face) &
                        result(return_value)
  !
  !.. Use Statements ..
  use geovar,            only : face_of_cell_t
  use interpolation_mod, only : interp
  !
  !.. Formal Arguments ..
  integer,              intent(in) :: cell_geom
  integer,              intent(in) :: cell_order
  real(wp),             intent(in) :: uu(:,:)
  real(wp),             intent(in) :: rnrm(:,:)
  type(face_of_cell_t), intent(in) :: face
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(uu,dim=1),1:size(rnrm,dim=2))
  !
  !.. Local Scalars ..
  integer :: k,kf,l,lf,m,face_order
  !
  !.. Local Arrays ..
  real(wp) :: flxsp(1:size(uu,dim=2), &
                    1:size(uu,dim=1), &
                    1:size(rnrm,dim=1))
  real(wp) :: vflx(1:size(rnrm,dim=1))
  real(wp) :: unrm(1:size(rnrm,dim=1))
  !
continue
  !
  face_order = face%order
  !
  do k = 1,size(uu,dim=2)
    flxsp(k,:,:) = invs_flx_vec_cv( uu(:,k) )
  end do
  !
  do kf = 1,size(rnrm,dim=2)
    !
    lf = face%get_fp_idx(kf)
    !
    unrm(:) = unit_vector( rnrm(:,kf) )
    !
    do m = 1,size(uu,dim=1)
      do l = 1,size(rnrm,dim=1)
        vflx(l) = interp(cell_geom,cell_order)%toFace(face_order)% &
                  dot( lf , flxsp(:,m,l) )
      end do
      !
      return_value(m,kf) = dot( vflx(:) , unrm(:) )
      !
    end do
    !
  end do
  !
end function invs_face_flx_cv
!
!###############################################################################
!
pure function invs_flx_vec_cv(uu) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  real(wp), intent(in) :: uu(:)
  !
  real(wp) :: return_value(1:size(uu),1:nr)
  !
  integer :: l
  !
continue
  !
  do l = 1,nr
    return_value(:,l) = invs_flx_cv( uu(:) , l )
  end do
  !
end function invs_flx_vec_cv
!
!###############################################################################
!
pure function invs_flx_cv1(uu) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : gam
  use ovar,   only : turb_is_on
  use ovar,   only : include_tke_in_energy
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: uu(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:nr,1:nq)
  !
  !.. Local Scalars ..
  integer  :: l,m
  real(wp) :: pres
  !
continue
  !
  pres = (gam-one)*(uu(nee) - half*selfdot(uu(nmb:nme))/uu(nec))
  !
  return_value(1:nr,nec) = uu(nmb:nme)
  !
  do m = nmb,nme
    do l = 1,nr
      return_value(l,m) = uu(m)*uu(nec+l)/uu(nec)
    end do
    return_value(m-nec,m) = return_value(m-nec,m) + pres
  end do
  !
  return_value(1:nr,nee) = (uu(nee) + pres) * uu(nmb:nme)/uu(nec)
  !
  if (turb_is_on) then
    return_value(1:nr,ntk) = uu(ntk) * uu(nmb:nme)/uu(nec)
    return_value(1:nr,ntl) = uu(ntl) * uu(nmb:nme)/uu(nec)
    if (include_tke_in_energy) then
      do l = 1,nr
        return_value(l,nec+l) = return_value(l,nec+l) - (gam-one)*uu(ntk)
        return_value(l,nee) = return_value(l,nee) - &
                              (gam-one)*uu(ntk)*uu(nec+l)/uu(nec)
      end do
    end if
  end if
  !
end function invs_flx_cv1
!
!###############################################################################
!
pure function invs_flx_cv(uu,dir) result(return_value)
  !
  use ovar, only : gam
  use ovar, only : turb_is_on
  use ovar, only : include_tke_in_energy
  !
  integer,                   intent(in) :: dir
  real(wp), dimension(1:nq), intent(in) :: uu
  !
  real(wp), dimension(1:nq) :: return_value
  !
  integer :: meq
  real(wp) :: pres,flxvel
  !
continue
  !
  meq = nec + dir
  !
  flxvel = uu(meq)/uu(nec)
  !
  pres = (gam-one)*(uu(nee) - half*selfdot(uu(nmb:nme))/uu(nec))
  !
  return_value(nec) = uu(meq)
  return_value(nmb:nme) = flxvel * uu(nmb:nme)
  return_value(meq) = return_value(meq) + pres
  return_value(nee) = flxvel * (uu(nee) + pres)
  !
  if (turb_is_on) then
    return_value(ntk) = flxvel * uu(ntk)
    return_value(ntl) = flxvel * uu(ntl)
    if (include_tke_in_energy) then
      return_value(meq) = return_value(meq) - (gam-one)*uu(ntk)
      return_value(nee) = return_value(nee) - (gam-one)*uu(ntk)*flxvel
    end if
  end if
  !
end function invs_flx_cv
!
!###############################################################################
!
pure function invs_flx_pv(vv,dir) result(return_value)
  !
  use ovar, only : gam
  use ovar, only : turb_is_on
  use ovar, only : include_tke_in_energy
  !
  integer,                   intent(in) :: dir
  real(wp), dimension(1:nq), intent(in) :: vv
  !
  real(wp), dimension(1:nq) :: return_value
  !
  integer :: meq
  real(wp) :: ener,flxvel
  !
continue
  !
  meq = nec + dir
  !
  flxvel = vv(meq)
  !
  ener = vv(nee)/(gam-one) + half*vv(nec)*selfdot(vv(nmb:nme))
  !
  return_value(nec) = flxvel * vv(nec)
  return_value(nmb:nme) = flxvel * vv(nec)*vv(nmb:nme)
  return_value(meq) = return_value(meq) + vv(nee)
  return_value(nee) = flxvel * (ener + vv(nee))
  !
  if (turb_is_on) then
    return_value(ntk) = flxvel * vv(nec)*vv(ntk)
    return_value(ntl) = flxvel * vv(nec)*vv(ntl)
    if (include_tke_in_energy) then
      return_value(meq) = return_value(meq) - (gam-one)*vv(nec)*vv(ntk)
      return_value(nee) = return_value(nee) - (gam-one)*vv(nec)*vv(ntk)*flxvel
    end if
  end if
  !
end function invs_flx_pv
!
!###############################################################################
!
pure function invs_flx_jacobian_cv(uu,dir) result(return_value)
  !
  use geovar, only : nr
  use ovar,   only : gam
  use ovar,   only : turb_is_on
  use ovar,   only : include_tke_in_energy
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
  ! ADDITIONAL 2*a SO ADD ANOTHER a TO THIS TERM
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
  if (turb_is_on) then
    ! TKE COLUMN, MOMENTUM ROW FOR CURRENT DIRECTION
    return_value(ntk,meq) = -gm1
    ! TKE COLUMN, TKE ROW
    return_value(ntk,ntk) = a
    ! TLS COLUMN, TLS ROW
    return_value(ntl,ntl) = a
    ! DENSITY COLUMN, TURBULENCE ROWS
    return_value(nec,ntk) = -a * uu(ntk)/uu(nec)
    return_value(nec,ntl) = -a * uu(ntl)/uu(nec)
    ! MOMENTUM COLUMN FOR CURRENT DIRECTION, TURBULENCE ROWS
    return_value(meq,ntk) = uu(ntk)/uu(nec)
    return_value(meq,ntl) = uu(ntl)/uu(nec)
    ! ADDITIONAL TERMS FOR ENERGY ROW IF TKE CONTRIBUTES TO ENERGY
    if (include_tke_in_energy) then
      ! DENSITY COLUMN, ENERGY ROW
      return_value(nec,nee) = return_value(nec,nee) + gm1*uu(ntk)
      ! MOMENTUM COLUMN FOR CURRENT DIRECTION, ENERGY ROW
      return_value(meq,nee) = return_value(meq,nee) - gm1*uu(ntk)/uu(nec)
      ! TKE COLUMN, ENERGY ROW
      return_value(ntk,nee) = -a * gm1
    end if
  end if
  !
end function invs_flx_jacobian_cv
!
!###############################################################################
!
pure function invs_wall_flx_cv(unrm,uu) result(return_value)
  !
  use geovar, only : nr
  use ovar,   only : gam
  use ovar,   only : turb_is_on
  use ovar,   only : include_tke_in_energy
  !
  real(wp), dimension(1:nr), intent(in) :: unrm
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
  if (include_tke_in_energy) pres = pres - (gam-one)*uu(ntk)
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
!
!###############################################################################
!
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
!
!###############################################################################
!
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
!
!###############################################################################
!
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
!
!###############################################################################
!
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
!
!###############################################################################
!
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
!
!###############################################################################
!
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
!
!###############################################################################
!
