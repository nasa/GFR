pure subroutine invs_flx_cv(uu,dir,return_value) 
! !$ACC ROUTINE SEQ
  !
  use ovar, only : gam
  use ovar, only : turb_is_on
  use ovar, only : include_tke_in_energy
  !
  integer,  value,              intent(in) :: dir
  real(wp), dimension(1:nq),    intent(in) :: uu
  !
  real(wp), dimension(1:nq), intent(inout) :: return_value
  !
  integer  :: n,meq
  real(wp) :: pres,flxvel
  !
continue
  !
 !meq = 0
  return_value(1:nq) = zero
 !return_value(:) = zero
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
end subroutine invs_flx_cv
!
!###############################################################################
!
pure subroutine invs_flx_cv_gpu_3d(usp,return_value) 
  !
  !.. Use Statements ..
  use ovar, only : gam
  use ovar, only : turb_is_on
  use ovar, only : include_tke_in_energy
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: usp(:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,np
  real(wp) :: velx,vely,velz,q2,gm1
  real(wp) :: rhou,rhov,rhow,rhoq2
  real(wp) :: rho,rho_inv,rhoe,pres
  real(wp) :: rhok,rhol,kflx
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: uu(:,:)
  !
continue
  !
  np = size(return_value,dim=1)
  gm1 = gam - one
  !
  allocate ( uu(1:np,1:nq) )
  !$ACC DATA CREATE(uu) &
  !$ACC&     PRESENT(return_value) &
  !$ACC&     COPYIN(usp)
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = 1,np
    uu(n,nec) = usp(nec,n)
    uu(n,nmx) = usp(nmx,n)
    uu(n,nmy) = usp(nmy,n)
    uu(n,nmz) = usp(nmz,n)
    uu(n,nee) = usp(nee,n)
  end do
  !
  ! Continuity equation
  !
  !$ACC PARALLEL LOOP
  do n = 1,np
    !
    rho = uu(n,nec)
    rho_inv = one/rho
    !
    velx = uu(n,nmx)*rho_inv
    vely = uu(n,nmy)*rho_inv
    velz = uu(n,nmz)*rho_inv
    rhoe = uu(n,nee)
    !
    q2 = velx*velx + vely*vely + velz*velz
    !
    pres = gm1*(rhoe - half*rho*q2)
    !
    return_value(n,1,nec) = rho*velx
    return_value(n,2,nec) = rho*vely
    return_value(n,3,nec) = rho*velz
    !
    return_value(n,1,nmx) = rho*velx*velx + pres
    return_value(n,2,nmx) = rho*velx*vely
    return_value(n,3,nmx) = rho*velx*velz
    !
    return_value(n,1,nmy) = rho*vely*velx
    return_value(n,2,nmy) = rho*vely*vely + pres
    return_value(n,3,nmy) = rho*vely*velz
    !
    return_value(n,1,nmz) = rho*velz*velx
    return_value(n,2,nmz) = rho*velz*vely
    return_value(n,3,nmz) = rho*velz*velz + pres
    !
    return_value(n,1,nee) = velx*(rhoe + pres)
    return_value(n,2,nee) = vely*(rhoe + pres)
    return_value(n,3,nee) = velz*(rhoe + pres)
    !
  end do
  !
  ! Turbulence equations
  !
  if (turb_is_on) then
    !$ACC PARALLEL LOOP
    do n = 1,np
      !
      rho_inv = one/uu(n,nec)
      velx = uu(n,nmx)*rho_inv
      vely = uu(n,nmy)*rho_inv
      velz = uu(n,nmz)*rho_inv
      rhok = uu(n,ntk)
      rhol = uu(n,ntl)
      !
      return_value(n,1,ntk) = velx * rhok
      return_value(n,2,ntk) = vely * rhok
      return_value(n,3,ntk) = velz * rhok
      !
      return_value(n,1,ntl) = velx * rhol
      return_value(n,2,ntl) = vely * rhol
      return_value(n,3,ntl) = velz * rhol
      !
    end do
    if (include_tke_in_energy) then
      !$ACC PARALLEL LOOP
      do n = 1,np
        !
        rho_inv = one/uu(n,nec)
        velx = uu(n,nmx)*rho_inv
        vely = uu(n,nmy)*rho_inv
        velz = uu(n,nmz)*rho_inv
        kflx = gm1*uu(n,ntk)
        !
        return_value(n,1,nmx) = return_value(n,1,nmx) - kflx
        return_value(n,2,nmy) = return_value(n,2,nmy) - kflx
        return_value(n,3,nmz) = return_value(n,3,nmz) - kflx
        !
        return_value(n,1,nee) = return_value(n,1,nee) - kflx*velx
        return_value(n,2,nee) = return_value(n,2,nee) - kflx*vely
        return_value(n,3,nee) = return_value(n,3,nee) - kflx*velz
        !
      end do
    end if
  end if
  !
  !$ACC END DATA
  deallocate ( uu )
  !
end subroutine invs_flx_cv_gpu_3d
!
!###############################################################################
!
pure subroutine invs_flx_cv_gpu_2d(usp,return_value) 
  !
  use ovar, only : gam
  use ovar, only : turb_is_on
  use ovar, only : include_tke_in_energy
  !
  real(wp),    intent(in) :: usp(:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  integer  :: l,m,n,np,meq
  real(wp) :: pres,gm1,val,flxvel,q2
  !
  real(wp), allocatable :: uu(:,:)
  !
continue
  !
  np = size(return_value,dim=1)
  gm1 = gam - one
  !
  allocate ( uu(1:np,1:nq) )
  !$ACC DATA CREATE(uu) &
  !$ACC&     PRESENT(return_value) &
  !$ACC&     COPYIN(usp)
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do n = 1,np
    do m = 1,nq
      uu(n,m) = usp(m,n)
    end do
  end do
  !
  ! Continuity equation
  !
  !$ACC PARALLEL LOOP COLLAPSE(2)
  do l = 1,2
    do n = 1,np
      return_value(n,l,nec) = uu(n,nec+l)
    end do
  end do
  !
  ! Momentum equations (without pressure flux)
  !
  !$ACC PARALLEL LOOP COLLAPSE(3)
  do m = nmb,nme
    do l = 1,2
      do n = 1,np
        return_value(n,l,m) = uu(n,m)*uu(n,nec+l)/uu(n,nec)
      end do
    end do
  end do
  !
  ! Energy equation and pressure flux
  !
  !$ACC PARALLEL LOOP COLLAPSE(2)
  do l = 1,2
    do n = 1,np
      meq = nec + l
      q2 = zero
      !$ACC LOOP SEQ
      do m = nmb,nme
        q2 = q2 + uu(n,m)*uu(n,m)
      end do
      pres = gm1*(uu(n,nee) - half*q2/uu(n,nec))
      return_value(n,l,meq) = return_value(n,l,meq) + pres
      return_value(n,l,nee) = (uu(n,nee) + pres)*uu(n,meq)/uu(n,nec)
    end do
  end do
  !
  ! Turbulence equations
  !
  if (turb_is_on) then
    !$ACC PARALLEL LOOP COLLAPSE(2)
    do l = 1,2
      do n = 1,np
        flxvel = uu(n,nec+l) / uu(n,nec)
        return_value(n,l,ntk) = flxvel * uu(n,ntk)
        return_value(n,l,ntl) = flxvel * uu(n,ntl)
      end do
    end do
    if (include_tke_in_energy) then
      !$ACC PARALLEL LOOP COLLAPSE(2)
      do l = 1,2
        do n = 1,np
          meq = nec + l
          val = gm1*uu(n,ntk)
          flxvel = uu(n,meq)/uu(n,nec)
          return_value(n,l,meq) = return_value(n,l,meq) - val
          return_value(n,l,nee) = return_value(n,l,nee) - val*flxvel
        end do
      end do
    end if
  end if
  !
  !$ACC END DATA
  deallocate ( uu )
  !
end subroutine invs_flx_cv_gpu_2d
!
!###############################################################################
!
pure subroutine invs_flx_cv_gpu_nd(usp,return_value) 
  !
  use ovar, only : gam
  use ovar, only : turb_is_on
  use ovar, only : include_tke_in_energy
  !
  real(wp),    intent(in) :: usp(:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  integer :: l,m,n,np,nr,meq
  real(wp) :: pres,gm1,val,flxvel,q2
  !
  real(wp), allocatable :: uu(:,:)
  !
continue
  !
  np = size(return_value,dim=1)
  nr = size(return_value,dim=2)
  gm1 = gam - one
  !
  allocate ( uu(1:np,1:nq) )
  !$ACC DATA CREATE(uu) &
  !$ACC&     PRESENT(return_value) &
  !$ACC&     COPYIN(usp)
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do n = 1,np
    do m = 1,nq
      uu(n,m) = usp(m,n)
    end do
  end do
  !
  ! Continuity equation
  !
  !$ACC PARALLEL LOOP COLLAPSE(2)
  do l = 1,nr
    do n = 1,np
      return_value(n,l,nec) = uu(n,nec+l)
    end do
  end do
  !
  ! Momentum equations (without pressure flux)
  !
  !$ACC PARALLEL LOOP COLLAPSE(3)
  do m = nmb,nme
    do l = 1,nr
      do n = 1,np
        return_value(n,l,m) = uu(n,m)*uu(n,nec+l)/uu(n,nec)
      end do
    end do
  end do
  !
  ! Energy equation and pressure flux
  !
  !$ACC PARALLEL LOOP COLLAPSE(2)
  do l = 1,nr
    do n = 1,np
      meq = nec + l
      q2 = zero
      !$ACC LOOP SEQ
      do m = nmb,nme
        q2 = q2 + uu(n,m)*uu(n,m)
      end do
      pres = gm1*(uu(n,nee) - half*q2/uu(n,nec))
      return_value(n,l,meq) = return_value(n,l,meq) + pres
      return_value(n,l,nee) = (uu(n,nee) + pres)*uu(n,meq)/uu(n,nec)
    end do
  end do
  !
  ! Turbulence equations
  !
  if (turb_is_on) then
    !$ACC PARALLEL LOOP COLLAPSE(2)
    do l = 1,nr
      do n = 1,np
        flxvel = uu(n,nec+l) / uu(n,nec)
        return_value(n,l,ntk) = flxvel * uu(n,ntk)
        return_value(n,l,ntl) = flxvel * uu(n,ntl)
      end do
    end do
    if (include_tke_in_energy) then
      !$ACC PARALLEL LOOP COLLAPSE(2)
      do l = 1,nr
        do n = 1,np
          meq = nec + l
          val = gm1*uu(n,ntk)
          flxvel = uu(n,meq)/uu(n,nec)
          return_value(n,l,meq) = return_value(n,l,meq) - val
          return_value(n,l,nee) = return_value(n,l,nee) - val*flxvel
        end do
      end do
    end if
  end if
  !
  !$ACC END DATA
  deallocate ( uu )
  !
end subroutine invs_flx_cv_gpu_nd
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
  return_value(:,:) = zero
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
pure subroutine invs_wall_flx_gpu_3d(n1,n2,facenrm,face_cv,return_value) 
  !
  !.. Use Statements ..
  use ovar,   only : gam
  use ovar,   only : include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,rho,rhou,rhov,rhow,rhoe,rhok,rhoq2,pres
  !
continue
  !
  gm1 = gam - one
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,return_value)
  do n = n1,n2
    !
    rho  = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    rhow = face_cv(n,nmz,1)
    rhoe = face_cv(n,nee,1)
    !
    rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
    !
    pres = gm1*(rhoe - half*rhoq2/rho)
    !
    return_value(n,nec) = zero
    return_value(n,nmx) = facenrm(1,n)*pres
    return_value(n,nmy) = facenrm(2,n)*pres
    return_value(n,nmz) = facenrm(3,n)*pres
    return_value(n,nee) = zero
    !
  end do
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,return_value)
    do n = n1,n2
      rhok = face_cv(n,ntk,1)
      return_value(n,nmx) = return_value(n,nmx) - gm1*rhok*facenrm(1,n)
      return_value(n,nmy) = return_value(n,nmy) - gm1*rhok*facenrm(2,n)
      return_value(n,nmz) = return_value(n,nmz) - gm1*rhok*facenrm(3,n)
    end do
  end if
  !
end subroutine invs_wall_flx_gpu_3d
!
!###############################################################################
!
pure subroutine invs_wall_flx_gpu_2d(n1,n2,facenrm,face_cv,return_value) 
  !
  !.. Use Statements ..
  use ovar,   only : gam
  use ovar,   only : include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,rho,rhou,rhov,rhoe,rhok,rhoq2,pres
  !
continue
  !
  gm1 = gam - one
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,return_value)
  do n = n1,n2
    !
    rho  = face_cv(n,nec,1)
    rhou = face_cv(n,nmx,1)
    rhov = face_cv(n,nmy,1)
    rhoe = face_cv(n,nee,1)
    !
    rhoq2 = rhou*rhou + rhov*rhov
    !
    pres = gm1*(rhoe - half*rhoq2/rho)
    !
    return_value(n,nec) = zero
    return_value(n,nmx) = facenrm(1,n)*pres
    return_value(n,nmy) = facenrm(2,n)*pres
    return_value(n,nee) = zero
    !
  end do
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,return_value)
    do n = n1,n2
      rhok = face_cv(n,ntk,1)
      return_value(n,nmx) = return_value(n,nmx) - gm1*rhok*facenrm(1,n)
      return_value(n,nmy) = return_value(n,nmy) - gm1*rhok*facenrm(2,n)
    end do
  end if
  !
end subroutine invs_wall_flx_gpu_2d
!
!###############################################################################
!
pure subroutine invs_wall_flx_gpu_nd(n1,n2,facenrm,face_cv,return_value) 
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : gam
  use ovar,   only : include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: l,m,n
  real(wp) :: gm1,pres,rhoq2
  !
continue
  !
  gm1 = gam - one
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,return_value)
  do n = n1,n2
    !
    rhoq2 = zero
    !$ACC LOOP SEQ
    do m = nmb,nme
      rhoq2 = rhoq2 + face_cv(n,m,1)*face_cv(n,m,1)
    end do
    !
    pres = gm1*(face_cv(n,nee,1) - half*rhoq2/face_cv(n,nec,1))
    !
    return_value(n,nec) = zero
    !$ACC LOOP SEQ
    do l = 1,nr
      return_value(n,nec+l) = facenrm(l,n)*pres
    end do
    return_value(n,nee) = zero
    !
  end do
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP COLLAPSE(2) &
    !$ACC&              PRESENT(facenrm,face_cv,return_value)
    do n = n1,n2
      do l = 1,nr
        return_value(n,nec+l) = return_value(n,nec+l) - &
                                gm1*face_cv(n,ntk,1)*facenrm(l,n)
      end do
    end do
  end if
  !
end subroutine invs_wall_flx_gpu_nd
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
