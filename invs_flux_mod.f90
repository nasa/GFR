module invs_flux_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,ntk,ntl,nq
  use eqn_idx, only : ncs,npres,ngamm,ng
  use timer_mod, only : solver_timer
  !
  implicit none
  !
  private
  !
  public :: invs_flx_cv
  public :: invs_flx_cv_gpu
  public :: invs_wall_flx_gpu
  public :: invs_wall_flx_cv
  !
contains
!
!###############################################################################
!
!pure &
subroutine invs_flx_cv_gpu(usp,gpsp,return_value,xyz)
  !
  !.. Use Statements ..
  use mrf_mod, only : mrf_flx_cv
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: usp(:,:)
  real(wp), intent(inout) :: gpsp(:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  real(wp),    intent(in) :: xyz(:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
  !. Local Saved Scalars ..
  integer, save :: idx = -1
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "invs_flx_cv_gpu"
  !
continue
  !
  call solver_timer(entering_procedure,pname,idx)
  !
  nr = size(return_value,dim=2)
  !
  if (nr == 3) then      ! this should also work for ncs > 0
    call invs_flx_cv_gpu_3d(usp,gpsp,return_value)
  else if (nr == 2) then ! this should also work for ncs > 0
    call invs_flx_cv_gpu_2d(usp,gpsp,return_value)
  else
    call invs_flx_cv_gpu_nd(usp,gpsp,return_value)
  end if
  !
  ! Add the flux terms due to the moving reference frame.
  ! NOTE: If not using a moving reference frame, this will just
  !       return without doing anything.
  !
  call mrf_flx_cv(usp,xyz,return_value)
  !
  call solver_timer(leaving_procedure,pname,idx)
  !
end subroutine invs_flx_cv_gpu
!
!###############################################################################
!
!pure &
subroutine invs_flx_cv_gpu_3d(usp,gpsp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  use ovar, only : turb_is_on
  use ovar, only : include_tke_in_energy
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: usp(:,:)
  real(wp), intent(inout) :: gpsp(:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,np
  real(wp) :: velx,vely,velz,gm1
  real(wp) :: rho,rho_inv,rhoe,pres
  real(wp) :: rhok,rhol,kflx,mass_frac
  !
continue
  !
  np = size(return_value,dim=1)
  gm1 = gam - one
  !
  ! Species continuity equations
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(return_value,usp)
    do m = 1,ncs
      do n = 1,np
        !
        mass_frac = usp(n,m)/usp(n,nec)
        !
        return_value(n,1,m) = mass_frac*usp(n,nmx)
        return_value(n,2,m) = mass_frac*usp(n,nmy)
        return_value(n,3,m) = mass_frac*usp(n,nmz)
        !
      end do
    end do
  end if
  !
  ! Continuity equation
  !
  !$ACC PARALLEL LOOP PRESENT(return_value,usp,gpsp)
  do n = 1,np
    !
    rho = usp(n,nec)
    rho_inv = one/rho
    !
    velx = usp(n,nmx)*rho_inv
    vely = usp(n,nmy)*rho_inv
    velz = usp(n,nmz)*rho_inv
    rhoe = usp(n,nee)
    !
    pres = gpsp(n,npres)
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
    !$ACC PARALLEL LOOP PRESENT(return_value,usp)
    do n = 1,np
      !
      rho_inv = one/usp(n,nec)
      velx = usp(n,nmx)*rho_inv
      vely = usp(n,nmy)*rho_inv
      velz = usp(n,nmz)*rho_inv
      rhok = usp(n,ntk)
      rhol = usp(n,ntl)
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
      if (ncs == 0) then
        !$ACC PARALLEL LOOP PRESENT(gpsp)
        do n = 1,np
          gpsp(n,ngamm) = gam
        end do
      end if
      !$ACC PARALLEL LOOP PRESENT(return_value,usp,gpsp)
      do n = 1,np
        !
        rho_inv = one/usp(n,nec)
        velx = usp(n,nmx)*rho_inv
        vely = usp(n,nmy)*rho_inv
        velz = usp(n,nmz)*rho_inv
        gm1 = gpsp(n,ngamm) - one
        kflx = gm1*usp(n,ntk)
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
end subroutine invs_flx_cv_gpu_3d
!
!###############################################################################
!
!pure &
subroutine invs_flx_cv_gpu_2d(usp,gpsp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  use ovar, only : turb_is_on
  use ovar, only : include_tke_in_energy
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: usp(:,:)
  real(wp), intent(inout) :: gpsp(:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,np
  real(wp) :: velx,vely,gm1
  real(wp) :: rho,rho_inv,rhoe,pres
  real(wp) :: rhok,rhol,kflx,mass_frac
  !
continue
  !
  np = size(return_value,dim=1)
  gm1 = gam - one
  !
  ! Species continuity equations
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(return_value,usp)
    do m = 1,ncs
      do n = 1,np
        !
        mass_frac = usp(n,m)/usp(n,nec)
        !
        return_value(n,1,m) = mass_frac*usp(n,nmx)
        return_value(n,2,m) = mass_frac*usp(n,nmy)
        !
      end do
    end do
  end if
  !
  ! Continuity equation
  !
  !$ACC PARALLEL LOOP PRESENT(return_value,usp,gpsp)
  do n = 1,np
    !
    rho = usp(n,nec)
    rho_inv = one/rho
    !
    velx = usp(n,nmx)*rho_inv
    vely = usp(n,nmy)*rho_inv
    rhoe = usp(n,nee)
    !
    pres = gpsp(n,npres)
    !
    return_value(n,1,nec) = rho*velx
    return_value(n,2,nec) = rho*vely
    !
    return_value(n,1,nmx) = rho*velx*velx + pres
    return_value(n,2,nmx) = rho*velx*vely
    !
    return_value(n,1,nmy) = rho*vely*velx
    return_value(n,2,nmy) = rho*vely*vely + pres
    !
    return_value(n,1,nee) = velx*(rhoe + pres)
    return_value(n,2,nee) = vely*(rhoe + pres)
    !
  end do
  !
  ! Turbulence equations
  !
  if (turb_is_on) then
    !$ACC PARALLEL LOOP PRESENT(return_value,usp)
    do n = 1,np
      !
      rho_inv = one/usp(n,nec)
      velx = usp(n,nmx)*rho_inv
      vely = usp(n,nmy)*rho_inv
      rhok = usp(n,ntk)
      rhol = usp(n,ntl)
      !
      return_value(n,1,ntk) = velx * rhok
      return_value(n,2,ntk) = vely * rhok
      !
      return_value(n,1,ntl) = velx * rhol
      return_value(n,2,ntl) = vely * rhol
      !
    end do
    if (include_tke_in_energy) then
      if (ncs == 0) then
        !$ACC PARALLEL LOOP PRESENT(gpsp)
        do n = 1,np
          gpsp(n,ngamm) = gam
        end do
      end if
      !$ACC PARALLEL LOOP PRESENT(return_value,usp,gpsp)
      do n = 1,np
        !
        rho_inv = one/usp(n,nec)
        velx = usp(n,nmx)*rho_inv
        vely = usp(n,nmy)*rho_inv
        gm1 = gpsp(n,ngamm) - one
        kflx = gm1*usp(n,ntk)
        !
        return_value(n,1,nmx) = return_value(n,1,nmx) - kflx
        return_value(n,2,nmy) = return_value(n,2,nmy) - kflx
        !
        return_value(n,1,nee) = return_value(n,1,nee) - kflx*velx
        return_value(n,2,nee) = return_value(n,2,nee) - kflx*vely
        !
      end do
    end if
  end if
  !
end subroutine invs_flx_cv_gpu_2d
!
!###############################################################################
!
!pure &
subroutine invs_flx_cv_gpu_nd(usp,gpsp,return_value)
  !
  use ovar, only : gam
  use ovar, only : turb_is_on
  use ovar, only : include_tke_in_energy
  !
  real(wp),    intent(in) :: usp(:,:)
  real(wp), intent(inout) :: gpsp(:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  integer  :: l,m,n,np,nr,meq
  real(wp) :: pres,gm1,val,flxvel
  !
continue
  !
  np = size(return_value,dim=1)
  nr = size(return_value,dim=2)
  gm1 = gam - one
  !
  ! Continuity equations and momentum equations (without pressure flux)
  !
  !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(return_value,usp)
  do m = 1,nme
    do l = 1,nr
      do n = 1,np
        return_value(n,l,m) = usp(n,m)*usp(n,nec+l)/usp(n,nec)
      end do
    end do
  end do
  !
  ! Energy equation and pressure flux
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(return_value,usp,gpsp)
  do l = 1,nr
    do n = 1,np
      meq = nec + l
      pres = gpsp(n,npres)
      return_value(n,l,meq) = return_value(n,l,meq) + pres
      return_value(n,l,nee) = (usp(n,nee) + pres)*usp(n,meq)/usp(n,nec)
    end do
  end do
  !
  ! Turbulence equations
  !
  if (turb_is_on) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(return_value,usp)
    do l = 1,nr
      do n = 1,np
        flxvel = usp(n,nec+l) / usp(n,nec)
        return_value(n,l,ntk) = flxvel * usp(n,ntk)
        return_value(n,l,ntl) = flxvel * usp(n,ntl)
      end do
    end do
    if (include_tke_in_energy) then
      if (ncs == 0) then
        !$ACC PARALLEL LOOP PRESENT(gpsp)
        do n = 1,np
          gpsp(n,ngamm) = gam
        end do
      end if
      !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(return_value,usp,gpsp)
      do l = 1,nr
        do n = 1,np
          meq = nec + l
          gm1 = gpsp(n,ngamm) - one
          val = gm1*usp(n,ntk)
          flxvel = usp(n,meq)/usp(n,nec)
          return_value(n,l,meq) = return_value(n,l,meq) - val
          return_value(n,l,nee) = return_value(n,l,nee) - val*flxvel
        end do
      end do
    end if
  end if
  !
end subroutine invs_flx_cv_gpu_nd
!
!###############################################################################
!
!pure &
subroutine invs_wall_flx_gpu(n1,n2,facenrm,facexyz,face_cv,face_gp,return_value)
  !
  !.. Use Statements ..
  use mrf_mod, only : mrf_wall_flx_cv
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: facexyz(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
  !. Local Saved Scalars ..
  integer, save :: idx = -1
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "invs_wall_flx_gpu"
  !
continue
  !
  call solver_timer(entering_procedure,pname,idx)
  !
  nr = size(facenrm,dim=1)
  !
  if (nr == 3) then      ! this should also work for ncs > 0
    call invs_wall_flx_gpu_3d(n1,n2,facenrm,face_cv,face_gp,return_value)
  else if (nr == 2) then ! this should also work for ncs > 0
    call invs_wall_flx_gpu_2d(n1,n2,facenrm,face_cv,face_gp,return_value)
  else
    call invs_wall_flx_gpu_nd(n1,n2,facenrm,face_cv,face_gp,return_value)
  end if
  !
  ! Add the flux terms due to the moving reference frame.
  ! NOTE: If not using a moving reference frame, this will just
  !       return without doing anything.
  !
  call mrf_wall_flx_cv(n1,n2,facenrm,facexyz,face_cv,return_value)
  !
  call solver_timer(leaving_procedure,pname,idx)
  !
end subroutine invs_wall_flx_gpu
!
!###############################################################################
!
!pure &
subroutine invs_wall_flx_gpu_3d(n1,n2,facenrm,face_cv,face_gp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  use ovar, only : include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: gm1,kflx,pres
  !
continue
  !
  gm1 = gam - one
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(return_value)
    do m = 1,ncs
      do n = n1,n2
        return_value(n,m) = zero
      end do
    end do
  end if
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,return_value)
  do n = n1,n2
    !
    pres = face_gp(n,npres,1)
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
    if (ncs == 0) then
      !$ACC PARALLEL LOOP PRESENT(face_gp)
      do n = n1,n2
        face_gp(n,ngamm,1) = gam
      end do
    end if
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,return_value)
    do n = n1,n2
      gm1 = face_gp(n,ngamm,1) - one
      kflx = gm1*face_cv(n,ntk,1)
      return_value(n,nmx) = return_value(n,nmx) - kflx*facenrm(1,n)
      return_value(n,nmy) = return_value(n,nmy) - kflx*facenrm(2,n)
      return_value(n,nmz) = return_value(n,nmz) - kflx*facenrm(3,n)
    end do
  end if
  !
end subroutine invs_wall_flx_gpu_3d
!
!###############################################################################
!
!pure &
subroutine invs_wall_flx_gpu_2d(n1,n2,facenrm,face_cv,face_gp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  use ovar, only : include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: gm1,pres,kflx
  !
continue
  !
  gm1 = gam - one
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(return_value)
    do m = 1,ncs
      do n = n1,n2
        return_value(n,m) = zero
      end do
    end do
  end if
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,return_value)
  do n = n1,n2
    !
    pres = face_gp(n,npres,1)
    !
    return_value(n,nec) = zero
    return_value(n,nmx) = facenrm(1,n)*pres
    return_value(n,nmy) = facenrm(2,n)*pres
    return_value(n,nee) = zero
    !
  end do
  !
  if (include_tke_in_energy) then
    if (ncs == 0) then
      !$ACC PARALLEL LOOP PRESENT(face_gp)
      do n = n1,n2
        face_gp(n,ngamm,1) = gam
      end do
    end if
    !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,return_value)
    do n = n1,n2
      gm1 = face_gp(n,ngamm,1) - one
      kflx = gm1*face_cv(n,ntk,1)
      return_value(n,nmx) = return_value(n,nmx) - kflx*facenrm(1,n)
      return_value(n,nmy) = return_value(n,nmy) - kflx*facenrm(2,n)
    end do
  end if
  !
end subroutine invs_wall_flx_gpu_2d
!
!###############################################################################
!
!pure &
subroutine invs_wall_flx_gpu_nd(n1,n2,facenrm,face_cv,face_gp,return_value)
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
  real(wp), intent(inout) :: face_gp(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: l,m,n
  real(wp) :: gm1,pres,kflx
  !
continue
  !
  gm1 = gam - one
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(return_value)
  do m = 1,nq
    do n = n1,n2
      return_value(n,m) = zero
    end do
  end do
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,return_value)
  do n = n1,n2
    !
    pres = face_gp(n,npres,1)
    !
    !$ACC LOOP SEQ
    do l = 1,nr
      return_value(n,nec+l) = facenrm(l,n)*pres
    end do
    !
  end do
  !
  if (include_tke_in_energy) then
    if (ncs == 0) then
      !$ACC PARALLEL LOOP PRESENT(face_gp)
      do n = n1,n2
        face_gp(n,ngamm,1) = gam
      end do
    end if
    !$ACC PARALLEL LOOP COLLAPSE(2) &
    !$ACC&              PRESENT(facenrm,face_cv,face_gp,return_value)
    do n = n1,n2
      do l = 1,nr
        gm1 = face_gp(n,ngamm,1) - one
        kflx = gm1*face_cv(n,ntk,1)
        return_value(n,nec+l) = return_value(n,nec+l) - kflx*facenrm(l,n)
      end do
    end do
  end if
  !
end subroutine invs_wall_flx_gpu_nd
!
!###############################################################################
!
pure subroutine invs_flx_cv(uu,gp,dir,return_value)
! !$ACC ROUTINE SEQ
  !
  use ovar, only : gam
  use ovar, only : turb_is_on
  use ovar, only : include_tke_in_energy
  !
  integer,  value,              intent(in) :: dir
  real(wp), dimension(1:nq),    intent(in) :: uu
  real(wp), dimension(1:ng),    intent(in) :: gp
  !
  real(wp), dimension(1:nq), intent(inout) :: return_value
  !
  integer  :: n,meq
  real(wp) :: flxvel
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
  return_value(1:nme) = flxvel * uu(1:nme)
  return_value(meq) = return_value(meq) + gp(npres)
  return_value(nee) = flxvel * (uu(nee) + gp(npres))
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
  return_value(1:nme) = flxvel * vv(nec)*vv(1:nme)
  return_value(nec) = flxvel * vv(nec)
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
end module invs_flux_mod
