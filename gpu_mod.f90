module gpu_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,ntk,ntl,nq
  !
  implicit none
  !
  private
  !
  public :: viscosity_cv_gpu
  public :: usp2v_gpu
  public :: tusp2v_gpu
  public :: vsp2u_gpu
  public :: tvsp2u_gpu
  public :: grad_cv_to_grad_pv_gpu
  !
contains
!
!###############################################################################
!
pure subroutine viscosity_cv_gpu(n1,n2,tusp,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: tusp(:,:)
  real(wp), intent(inout) :: return_value(:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call viscosity_cv_gpu_3d(n1,n2,tusp,return_value)
  else if (nr == 2) then
    call viscosity_cv_gpu_2d(n1,n2,tusp,return_value)
  else
    call viscosity_cv_gpu_nd(n1,n2,tusp,return_value)
  end if
  !
end subroutine viscosity_cv_gpu
!
!###############################################################################
!
pure subroutine viscosity_cv_gpu_3d(n1,n2,tusp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : tref,rhoref,aref
  use ovar, only : suth_c1,suth_sref
  use ovar, only : gam,rgasref_nd
  use ovar, only : include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: tusp(:,:)
  real(wp), intent(inout) :: return_value(:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: cv_inv,rhoaref_inv
  real(wp) :: rho_inv,rhou,rhov,rhow,rhoq2
  real(wp) :: rhoe,rhoe_int,tdim,mudim
  !
continue
  !
  cv_inv = (gam-one)/rgasref_nd
  rhoaref_inv = one / rhoref / aref
  !
  !$ACC DATA PRESENT(tusp,return_value)
  !
  !$ACC PARALLEL LOOP
  do n = n1,n2
    !
    rho_inv = one/tusp(n,nec)
    !
    rhou = tusp(n,nmx)
    rhov = tusp(n,nmy)
    rhow = tusp(n,nmz)
    rhoe = tusp(n,nee)
    !
    rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
    !
    rhoe_int = rhoe - half*rhoq2*rho_inv
    !
    ! Temporarily store the temperature in the return_value array
    !
    return_value(n) = cv_inv * rhoe_int * rho_inv
    !
  end do
  !
  ! Subtract out the turbulence kinetic energy part
  ! of the internal energy from the temperature
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP
    do n = n1,n2
      return_value(n) = return_value(n) - cv_inv*tusp(n,ntk)/tusp(n,nec)
    end do
  end if
  !
  ! Compute the viscosity using Sutherlands law
  !
  !$ACC PARALLEL LOOP
  do n = n1,n2
    tdim = return_value(n)*tref
    mudim = suth_c1 * tdim * sqrt(tdim) / (tdim + suth_sref)
    return_value(n) = mudim * rhoaref_inv
  end do
  !
  !$ACC END DATA
  !
end subroutine viscosity_cv_gpu_3d
!
!###############################################################################
!
pure subroutine viscosity_cv_gpu_2d(n1,n2,tusp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : tref,rhoref,aref
  use ovar, only : suth_c1,suth_sref
  use ovar, only : gam,rgasref_nd
  use ovar, only : include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: tusp(:,:)
  real(wp), intent(inout) :: return_value(:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: cv_inv,rhoaref_inv
  real(wp) :: rho_inv,rhou,rhov,rhoq2
  real(wp) :: rhoe,rhoe_int,tdim,mudim
  !
continue
  !
  cv_inv = (gam-one)/rgasref_nd
  rhoaref_inv = one / rhoref / aref
  !
  !$ACC DATA PRESENT(tusp,return_value)
  !
  !$ACC PARALLEL LOOP
  do n = n1,n2
    !
    rho_inv = one/tusp(n,nec)
    !
    rhou = tusp(n,nmx)
    rhov = tusp(n,nmy)
    rhoe = tusp(n,nee)
    !
    rhoq2 = rhou*rhou + rhov*rhov
    !
    rhoe_int = rhoe - half*rhoq2*rho_inv
    !
    ! Temporarily store the temperature in the return_value array
    !
    return_value(n) = cv_inv * rhoe_int * rho_inv
    !
  end do
  !
  ! Subtract out the turbulence kinetic energy part
  ! of the internal energy from the temperature
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP
    do n = n1,n2
      return_value(n) = return_value(n) - cv_inv*tusp(n,ntk)/tusp(n,nec)
    end do
  end if
  !
  ! Compute the viscosity using Sutherlands law
  !
  !$ACC PARALLEL LOOP
  do n = n1,n2
    tdim = return_value(n)*tref
    mudim = suth_c1 * tdim * sqrt(tdim) / (tdim + suth_sref)
    return_value(n) = mudim * rhoaref_inv
  end do
  !
  !$ACC END DATA
  !
end subroutine viscosity_cv_gpu_2d
!
!###############################################################################
!
pure subroutine viscosity_cv_gpu_nd(n1,n2,tusp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : tref,rhoref,aref
  use ovar, only : suth_c1,suth_sref
  use ovar, only : gam,rgasref_nd
  use ovar, only : include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: tusp(:,:)
  real(wp), intent(inout) :: return_value(:)
  !
  !.. Local Scalars ..
  integer  :: nr,n
  real(wp) :: cv_inv,rhoaref_inv
  real(wp) :: rho_inv,rhou,rhov,rhow,rhoq2
  real(wp) :: rhoe,rhoe_int,tdim,mudim
  !
continue
  !
  nr = nme-nmb+1
  !
  cv_inv = (gam-one)/rgasref_nd
  rhoaref_inv = one / rhoref / aref
  !
  !$ACC DATA PRESENT(tusp,return_value)
  !
  if (nr == 3) then
    !
    !$ACC PARALLEL LOOP
    do n = n1,n2
      !
      rho_inv = one/tusp(n,nec)
      !
      rhou = tusp(n,nmx)
      rhov = tusp(n,nmy)
      rhow = tusp(n,nmz)
      rhoe = tusp(n,nee)
      !
      rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
      !
      rhoe_int = rhoe - half*rhoq2*rho_inv
      !
      ! Temporarily store the temperature in the return_value array
      !
      return_value(n) = cv_inv * rhoe_int * rho_inv
      !
    end do
    !
  else if (nr == 2) then
    !
    !$ACC PARALLEL LOOP
    do n = n1,n2
      !
      rho_inv = one/tusp(n,nec)
      !
      rhou = tusp(n,nmx)
      rhov = tusp(n,nmy)
      rhoe = tusp(n,nee)
      !
      rhoq2 = rhou*rhou + rhov*rhov
      !
      rhoe_int = rhoe - half*rhoq2*rho_inv
      !
      ! Temporarily store the temperature in the return_value array
      !
      return_value(n) = cv_inv * rhoe_int * rho_inv
      !
    end do
    !
  else
    !
    !$ACC PARALLEL LOOP
    do n = n1,n2
      !
      rho_inv = one/tusp(n,nec)
      !
      rhou = tusp(n,nmx)
      rhoe = tusp(n,nee)
      !
      rhoq2 = rhou*rhou
      !
      rhoe_int = rhoe - half*rhoq2*rho_inv
      !
      ! Temporarily store the temperature in the return_value array
      !
      return_value(n) = cv_inv * rhoe_int * rho_inv
      !
    end do
    !
  end if
  !
  ! Subtract out the turbulence kinetic energy part
  ! of the internal energy from the temperature
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP
    do n = n1,n2
      return_value(n) = return_value(n) - cv_inv*tusp(n,ntk)/tusp(n,nec)
    end do
  end if
  !
  ! Compute the viscosity using Sutherlands law
  !
  !$ACC PARALLEL LOOP
  do n = n1,n2
    tdim = return_value(n)*tref
    mudim = suth_c1 * tdim * sqrt(tdim) / (tdim + suth_sref)
    return_value(n) = mudim * rhoaref_inv
  end do
  !
  !$ACC END DATA
  !
end subroutine viscosity_cv_gpu_nd
!
!###############################################################################
!
pure subroutine usp2v_gpu(n1,n2,usp,return_value,swap_t_for_rho)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: usp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: swap_t_for_rho
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call usp2v_gpu_3d(n1,n2,usp,return_value,swap_t_for_rho)
  else if (nr == 2) then
    call usp2v_gpu_2d(n1,n2,usp,return_value,swap_t_for_rho)
  else
    call usp2v_gpu_nd(n1,n2,usp,return_value,swap_t_for_rho)
  end if
  !
end subroutine usp2v_gpu
!
!###############################################################################
!
pure subroutine usp2v_gpu_3d(n1,n2,usp,return_value,swap_t_for_rho)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd,include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: usp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: swap_t_for_rho
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,rgasref_nd_inv
  real(wp) :: rho,rho_inv,velx,vely,velz,rhoe,q2,pres
  !
continue
  !
  gm1 = gam - one
  rgasref_nd_inv = one/rgasref_nd
  !
  !$ACC DATA PRESENT(usp,return_value)
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = n1,n2
    !
    rho = usp(nec,n)
    rho_inv = one/rho
    !
    velx = usp(nmx,n)*rho_inv
    vely = usp(nmy,n)*rho_inv
    velz = usp(nmz,n)*rho_inv
    rhoe = usp(nee,n)
    !
    q2 = velx*velx + vely*vely + velz*velz
    !
    pres = gm1*(rhoe - half*rho*q2)
    !
    return_value(nec,n) = rho
    return_value(nmx,n) = velx
    return_value(nmy,n) = vely
    return_value(nmz,n) = velz
    return_value(nee,n) = pres
    !
  end do
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
    do n = n1,n2
      rho_inv = one/usp(nec,n)
      return_value(nee,n) = return_value(nee,n) - gm1*usp(ntk,n)
      return_value(ntk,n) = usp(ntk,n)*rho_inv
      return_value(ntl,n) = usp(ntl,n)*rho_inv
    end do
  end if
  !
  if (present(swap_t_for_rho)) then
    if (swap_t_for_rho) then
      !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
      do n = n1,n2
        rho = return_value(nec,n)
        pres = return_value(nee,n)
        return_value(nec,n) = pres * rgasref_nd_inv / rho 
      end do
    end if
  end if
  !
  !$ACC END DATA
  !
end subroutine usp2v_gpu_3d
!
!###############################################################################
!
pure subroutine usp2v_gpu_2d(n1,n2,usp,return_value,swap_t_for_rho)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd,include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: usp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: swap_t_for_rho
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,rgasref_nd_inv
  real(wp) :: rho,rho_inv,velx,vely,rhoe,q2,pres
  !
continue
  !
  gm1 = gam - one
  rgasref_nd_inv = one/rgasref_nd
  !
  !$ACC DATA PRESENT(usp,return_value)
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = n1,n2
    !
    rho = usp(nec,n)
    rho_inv = one/rho
    !
    velx = usp(nmx,n)*rho_inv
    vely = usp(nmy,n)*rho_inv
    rhoe = usp(nee,n)
    !
    q2 = velx*velx + vely*vely
    !
    pres = gm1*(rhoe - half*rho*q2)
    !
    return_value(nec,n) = rho
    return_value(nmx,n) = velx
    return_value(nmy,n) = vely
    return_value(nee,n) = pres
    !
  end do
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
    do n = n1,n2
      rho_inv = one/usp(nec,n)
      return_value(nee,n) = return_value(nee,n) - gm1*usp(ntk,n)
      return_value(ntk,n) = usp(ntk,n)*rho_inv
      return_value(ntl,n) = usp(ntl,n)*rho_inv
    end do
  end if
  !
  if (present(swap_t_for_rho)) then
    if (swap_t_for_rho) then
      !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
      do n = n1,n2
        rho = return_value(nec,n)
        pres = return_value(nee,n)
        return_value(nec,n) = pres * rgasref_nd_inv / rho 
      end do
    end if
  end if
  !
  !$ACC END DATA
  !
end subroutine usp2v_gpu_2d
!
!###############################################################################
!
pure subroutine usp2v_gpu_nd(n1,n2,usp,return_value,swap_t_for_rho)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd,include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: usp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: swap_t_for_rho
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: gm1,rgasref_nd_inv
  real(wp) :: rho,rho_inv,velx,vely,velz,rhoe,rhoq2,pres
  !
continue
  !
  gm1 = gam - one
  rgasref_nd_inv = one/rgasref_nd
  !
  !$ACC DATA PRESENT(usp,return_value)
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = n1,n2
    !
    rhoq2 = zero
    !$ACC LOOP SEQ
    do m = nmb,nme
      rhoq2 = rhoq2 + usp(m,n)*usp(m,n)
    end do
    !
    return_value(nec,n)     = usp(nec,n)
    return_value(nmb:nme,n) = usp(nmb:nme,n)/usp(nec,n)
    return_value(nee,n)     = gm1 * (usp(nee,n) - half*rhoq2/usp(nec,n))
    !
  end do
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
    do n = n1,n2
      rho_inv = one/usp(nec,n)
      return_value(nee,n) = return_value(nee,n) - gm1*usp(ntk,n)
      return_value(ntk,n) = usp(ntk,n)*rho_inv
      return_value(ntl,n) = usp(ntl,n)*rho_inv
    end do
  end if
  !
  if (present(swap_t_for_rho)) then
    if (swap_t_for_rho) then
      !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
      do n = n1,n2
        rho = return_value(nec,n)
        pres = return_value(nee,n)
        return_value(nec,n) = pres * rgasref_nd_inv / rho 
      end do
    end if
  end if
  !
  !$ACC END DATA
  !
end subroutine usp2v_gpu_nd
!
!###############################################################################
!
pure subroutine tusp2v_gpu(n1,n2,tusp,return_value,swap_t_for_rho)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: tusp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: swap_t_for_rho
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call tusp2v_gpu_3d(n1,n2,tusp,return_value,swap_t_for_rho)
  else if (nr == 2) then
    call tusp2v_gpu_2d(n1,n2,tusp,return_value,swap_t_for_rho)
  else
    call tusp2v_gpu_nd(n1,n2,tusp,return_value,swap_t_for_rho)
  end if
  !
end subroutine tusp2v_gpu
!
!###############################################################################
!
pure subroutine tusp2v_gpu_3d(n1,n2,tusp,return_value,swap_t_for_rho)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd,include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: tusp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: swap_t_for_rho
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,rgasref_nd_inv
  real(wp) :: rho,rho_inv,velx,vely,velz,rhoe,q2,pres
  !
continue
  !
  gm1 = gam - one
  rgasref_nd_inv = one/rgasref_nd
  !
  !$ACC DATA PRESENT(tusp,return_value)
  !
  !$ACC PARALLEL LOOP
  do n = n1,n2
    !
    rho = tusp(n,nec)
    rho_inv = one/rho
    !
    velx = tusp(n,nmx)*rho_inv
    vely = tusp(n,nmy)*rho_inv
    velz = tusp(n,nmz)*rho_inv
    rhoe = tusp(n,nee)
    !
    q2 = velx*velx + vely*vely + velz*velz
    !
    pres = gm1*(rhoe - half*rho*q2)
    !
    return_value(n,nec) = rho
    return_value(n,nmx) = velx
    return_value(n,nmy) = vely
    return_value(n,nmz) = velz
    return_value(n,nee) = pres
    !
  end do
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP
    do n = n1,n2
      rho_inv = one/tusp(n,nec)
      return_value(n,nee) = return_value(n,nee) - gm1*tusp(n,ntk)
      return_value(n,ntk) = tusp(n,ntk)*rho_inv
      return_value(n,ntl) = tusp(n,ntl)*rho_inv
    end do
  end if
  !
  if (present(swap_t_for_rho)) then
    if (swap_t_for_rho) then
      !$ACC PARALLEL LOOP
      do n = n1,n2
        rho = return_value(n,nec)
        pres = return_value(n,nee)
        return_value(n,nec) = pres * rgasref_nd_inv / rho 
      end do
    end if
  end if
  !
  !$ACC END DATA
  !
end subroutine tusp2v_gpu_3d
!
!###############################################################################
!
pure subroutine tusp2v_gpu_2d(n1,n2,tusp,return_value,swap_t_for_rho)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd,include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: tusp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: swap_t_for_rho
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,rgasref_nd_inv
  real(wp) :: rho,rho_inv,velx,vely,rhoe,q2,pres
  !
continue
  !
  gm1 = gam - one
  rgasref_nd_inv = one/rgasref_nd
  !
  !$ACC DATA PRESENT(tusp,return_value)
  !
  !$ACC PARALLEL LOOP
  do n = n1,n2
    !
    rho = tusp(n,nec)
    rho_inv = one/rho
    !
    velx = tusp(n,nmx)*rho_inv
    vely = tusp(n,nmy)*rho_inv
    rhoe = tusp(n,nee)
    !
    q2 = velx*velx + vely*vely
    !
    pres = gm1*(rhoe - half*rho*q2)
    !
    return_value(n,nec) = rho
    return_value(n,nmx) = velx
    return_value(n,nmy) = vely
    return_value(n,nee) = pres
    !
  end do
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP
    do n = n1,n2
      rho_inv = one/tusp(n,nec)
      return_value(n,nee) = return_value(n,nee) - gm1*tusp(n,ntk)
      return_value(n,ntk) = tusp(n,ntk)*rho_inv
      return_value(n,ntl) = tusp(n,ntl)*rho_inv
    end do
  end if
  !
  if (present(swap_t_for_rho)) then
    if (swap_t_for_rho) then
      !$ACC PARALLEL LOOP
      do n = n1,n2
        rho = return_value(n,nec)
        pres = return_value(n,nee)
        return_value(n,nec) = pres * rgasref_nd_inv / rho 
      end do
    end if
  end if
  !
  !$ACC END DATA
  !
end subroutine tusp2v_gpu_2d
!
!###############################################################################
!
pure subroutine tusp2v_gpu_nd(n1,n2,tusp,return_value,swap_t_for_rho)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd,include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: tusp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: swap_t_for_rho
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: gm1,rgasref_nd_inv
  real(wp) :: rho,rho_inv,rhoq2,pres
  !
continue
  !
  gm1 = gam - one
  rgasref_nd_inv = one/rgasref_nd
  !
  !$ACC DATA PRESENT(tusp,return_value)
  !
  !$ACC PARALLEL LOOP
  do n = n1,n2
    !
    rhoq2 = zero
    !$ACC LOOP SEQ
    do m = nmb,nme
      rhoq2 = rhoq2 + tusp(n,m)*tusp(n,m)
    end do
    !
    return_value(n,nec)     = tusp(n,nec)
    return_value(n,nmb:nme) = tusp(n,nmb:nme)/tusp(n,nec)
    return_value(n,nee)     = gm1 * (tusp(n,nee) - half*rhoq2/tusp(n,nec))
    !
  end do
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP
    do n = n1,n2
      rho_inv = one/tusp(n,nec)
      return_value(n,nee) = return_value(n,nee) - gm1*tusp(n,ntk)
      return_value(n,ntk) = tusp(n,ntk)*rho_inv
      return_value(n,ntl) = tusp(n,ntl)*rho_inv
    end do
  end if
  !
  if (present(swap_t_for_rho)) then
    if (swap_t_for_rho) then
      !$ACC PARALLEL LOOP
      do n = n1,n2
        rho = return_value(n,nec)
        pres = return_value(n,nee)
        return_value(n,nec) = pres * rgasref_nd_inv / rho 
      end do
    end if
  end if
  !
  !$ACC END DATA
  !
end subroutine tusp2v_gpu_nd
!
!###############################################################################
!
pure subroutine vsp2u_gpu(n1,n2,vsp,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: vsp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call vsp2u_gpu_3d(n1,n2,vsp,return_value)
  else if (nr == 2) then
    call vsp2u_gpu_2d(n1,n2,vsp,return_value)
  else
    call vsp2u_gpu_nd(n1,n2,vsp,return_value)
  end if
  !
end subroutine vsp2u_gpu
!
!###############################################################################
!
pure subroutine vsp2u_gpu_3d(n1,n2,vsp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam, include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: vsp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1_inv,q2,rho,velx,vely,velz,pres,tke,tls
  !
continue
  !
  gm1_inv = one/(gam-one)
  !
  !$ACC DATA PRESENT(vsp,return_value)
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = n1,n2
    !
    rho = vsp(nec,n)
    velx = vsp(nmx,n)
    vely = vsp(nmy,n)
    velz = vsp(nmz,n)
    pres = vsp(nee,n)
    !
    q2 = velx*velx + vely*vely + velz*velz
    !
    return_value(nec,n) = rho
    return_value(nmx,n) = rho*velx
    return_value(nmy,n) = rho*vely
    return_value(nmz,n) = rho*velz
    return_value(nee,n) = pres*gm1_inv + half*rho*q2
    !
  end do
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
    do n = n1,n2
      rho = vsp(nec,n)
      tke = vsp(ntk,n)
      tls = vsp(ntl,n)
      return_value(ntk,n) = rho*tke
      return_value(ntl,n) = rho*tls
      return_value(nee,n) = return_value(nee,n) + rho*tke
    end do
  end if
  !
  !$ACC END DATA
  !
end subroutine vsp2u_gpu_3d
!
!###############################################################################
!
pure subroutine vsp2u_gpu_2d(n1,n2,vsp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam, include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: vsp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1_inv,q2,rho,velx,vely,pres,tke,tls
  !
continue
  !
  gm1_inv = one/(gam-one)
  !
  !$ACC DATA PRESENT(vsp,return_value)
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = n1,n2
    !
    rho = vsp(nec,n)
    velx = vsp(nmx,n)
    vely = vsp(nmy,n)
    pres = vsp(nee,n)
    !
    q2 = velx*velx + vely*vely
    !
    return_value(nec,n) = rho
    return_value(nmx,n) = rho*velx
    return_value(nmy,n) = rho*vely
    return_value(nee,n) = pres*gm1_inv + half*rho*q2
    !
  end do
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
    do n = n1,n2
      rho = vsp(nec,n)
      tke = vsp(ntk,n)
      tls = vsp(ntl,n)
      return_value(ntk,n) = rho*tke
      return_value(ntl,n) = rho*tls
      return_value(nee,n) = return_value(nee,n) + rho*tke
    end do
  end if
  !
  !$ACC END DATA
  !
end subroutine vsp2u_gpu_2d
!
!###############################################################################
!
pure subroutine vsp2u_gpu_nd(n1,n2,vsp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam, include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: vsp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: gm1_inv,q2,rho,tke,tls
  !
continue
  !
  gm1_inv = one/(gam-one)
  !
  !$ACC DATA PRESENT(vsp,return_value)
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = n1,n2
    !
    rho = vsp(nec,n)
    !
    q2 = zero
    !$ACC LOOP SEQ
    do m = nmb,nme
      q2 = q2 + vsp(m,n)*vsp(m,n)
    end do
    !
    return_value(nec,n)     = rho
    return_value(nmb:nme,n) = rho*vsp(nmb:nme,n)
    return_value(nee,n)     = vsp(nee,n)*gm1_inv + half*rho*q2
    !
  end do
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
    do n = n1,n2
      rho = vsp(nec,n)
      tke = vsp(ntk,n)
      tls = vsp(ntl,n)
      return_value(ntk,n) = rho*tke
      return_value(ntl,n) = rho*tls
      return_value(nee,n) = return_value(nee,n) + rho*tke
    end do
  end if
  !
  !$ACC END DATA
  !
end subroutine vsp2u_gpu_nd
!
!###############################################################################
!
pure subroutine tvsp2u_gpu(n1,n2,tvsp,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: tvsp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call tvsp2u_gpu_3d(n1,n2,tvsp,return_value)
  else if (nr == 2) then
    call tvsp2u_gpu_2d(n1,n2,tvsp,return_value)
  else
    call tvsp2u_gpu_nd(n1,n2,tvsp,return_value)
  end if
  !
end subroutine tvsp2u_gpu
!
!###############################################################################
!
pure subroutine tvsp2u_gpu_3d(n1,n2,tvsp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam, include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: tvsp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1_inv,q2,rho,velx,vely,velz,pres,tke,tls
  !
continue
  !
  gm1_inv = one/(gam-one)
  !
  !$ACC DATA PRESENT(tvsp,return_value)
  !
  !$ACC PARALLEL LOOP
  do n = n1,n2
    !
    rho = tvsp(n,nec)
    velx = tvsp(n,nmx)
    vely = tvsp(n,nmy)
    velz = tvsp(n,nmz)
    pres = tvsp(n,nee)
    !
    q2 = velx*velx + vely*vely + velz*velz
    !
    return_value(n,nec) = rho
    return_value(n,nmx) = rho*velx
    return_value(n,nmy) = rho*vely
    return_value(n,nmz) = rho*velz
    return_value(n,nee) = pres*gm1_inv + half*rho*q2
    !
  end do
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP
    do n = n1,n2
      rho = tvsp(n,nec)
      tke = tvsp(n,ntk)
      tls = tvsp(n,ntl)
      return_value(n,ntk) = rho*tke
      return_value(n,ntl) = rho*tls
      return_value(n,nee) = return_value(n,nee) + rho*tke
    end do
  end if
  !
  !$ACC END DATA
  !
end subroutine tvsp2u_gpu_3d
!
!###############################################################################
!
pure subroutine tvsp2u_gpu_2d(n1,n2,tvsp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam, include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: tvsp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1_inv,q2,rho,velx,vely,pres,tke,tls
  !
continue
  !
  gm1_inv = one/(gam-one)
  !
  !$ACC DATA PRESENT(tvsp,return_value)
  !
  !$ACC PARALLEL LOOP
  do n = n1,n2
    !
    rho = tvsp(n,nec)
    velx = tvsp(n,nmx)
    vely = tvsp(n,nmy)
    pres = tvsp(n,nee)
    !
    q2 = velx*velx + vely*vely
    !
    return_value(n,nec) = rho
    return_value(n,nmx) = rho*velx
    return_value(n,nmy) = rho*vely
    return_value(n,nee) = pres*gm1_inv + half*rho*q2
    !
  end do
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP
    do n = n1,n2
      rho = tvsp(n,nec)
      tke = tvsp(n,ntk)
      tls = tvsp(n,ntl)
      return_value(n,ntk) = rho*tke
      return_value(n,ntl) = rho*tls
      return_value(n,nee) = return_value(n,nee) + rho*tke
    end do
  end if
  !
  !$ACC END DATA
  !
end subroutine tvsp2u_gpu_2d
!
!###############################################################################
!
pure subroutine tvsp2u_gpu_nd(n1,n2,tvsp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam, include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: tvsp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: gm1_inv,q2,rho,tke,tls
  !
continue
  !
  gm1_inv = one/(gam-one)
  !
  !$ACC DATA PRESENT(tvsp,return_value)
  !
  !$ACC PARALLEL LOOP
  do n = n1,n2
    !
    rho = tvsp(n,nec)
    !
    q2 = zero
    !$ACC LOOP SEQ
    do m = nmb,nme
      q2 = q2 + tvsp(n,m)*tvsp(n,m)
    end do
    !
    return_value(n,nec)     = rho
    return_value(n,nmb:nme) = rho*tvsp(n,nmb:nme)
    return_value(n,nee)     = tvsp(n,nee)*gm1_inv + half*rho*q2
    !
  end do
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP
    do n = n1,n2
      rho = tvsp(n,nec)
      tke = tvsp(n,ntk)
      tls = tvsp(n,ntl)
      return_value(n,ntk) = rho*tke
      return_value(n,ntl) = rho*tls
      return_value(n,nee) = return_value(n,nee) + rho*tke
    end do
  end if
  !
  !$ACC END DATA
  !
end subroutine tvsp2u_gpu_nd
!
!###############################################################################
!
pure subroutine grad_cv_to_grad_pv_gpu(n1,n2,tusp,tdcvdx,return_value, &
                                       tusp_includes_rho_inv)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: tusp(:,:)
  real(wp),    intent(in) :: tdcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: tusp_includes_rho_inv
  !
  !.. Local Scalars ..
  integer  :: nr
  !
continue
  !
  nr = size(tdcvdx,dim=2)
  !
  if (nr == 3) then
    call grad_cv_to_grad_pv_gpu_3d(n1,n2,tusp,tdcvdx,return_value, &
                                   tusp_includes_rho_inv)
  else if (nr == 2) then
    call grad_cv_to_grad_pv_gpu_2d(n1,n2,tusp,tdcvdx,return_value, &
                                   tusp_includes_rho_inv)
  else
    call grad_cv_to_grad_pv_gpu_nd(n1,n2,tusp,tdcvdx,return_value, &
                                   tusp_includes_rho_inv)
  end if
  !
end subroutine grad_cv_to_grad_pv_gpu
!
!###############################################################################
!
pure subroutine grad_cv_to_grad_pv_gpu_3d(n1,n2,tusp,tdcvdx,return_value, &
                                          tusp_includes_rho_inv)
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: tusp(:,:)
  real(wp),    intent(in) :: tdcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: tusp_includes_rho_inv
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: cv_inv,rho,rho_inv
  real(wp) :: velx,vely,velz,etot
  real(wp) :: drhodx,drhody,drhodz
  !
  logical(lk) :: tusp_needs_rho_inv
  !
continue
  !
  tusp_needs_rho_inv = true
  if (present(tusp_includes_rho_inv)) then
    tusp_needs_rho_inv = .not. tusp_includes_rho_inv
  end if
  !
  !$ACC DATA PRESENT(tusp,tdcvdx,return_value)
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  ! Embed 1/rho in the all variables except density
  !
  if (tusp_needs_rho_inv) then
    !$ACC PARALLEL LOOP
    do n = n1,n2
      rho_inv = one/tusp(n,nec)
      tusp(n,nmx) = tusp(n,nmx)*rho_inv
      tusp(n,nmy) = tusp(n,nmy)*rho_inv
      tusp(n,nmz) = tusp(n,nmz)*rho_inv
      tusp(n,nee) = tusp(n,nee)*rho_inv
    end do
  end if
  !
  ! Copy the density gradient
  !
  !$ACC PARALLEL LOOP
  do n = n1,n2
    return_value(n,1,nec) = tdcvdx(n,1,nec)
    return_value(n,2,nec) = tdcvdx(n,2,nec)
    return_value(n,3,nec) = tdcvdx(n,3,nec)
  end do
  !
  ! Using the chain rule:
  !
  ! 1. Compute the velocity gradients (nmx,nmy,nmz)
  ! 2. Compute the total energy gradient and store in the
  !    location for the temperature gradient (nee)
  !
  !$ACC PARALLEL LOOP
  do n = n1,n2
    !
    rho_inv = one/tusp(n,nec)  ! tusp already includes 1/rho meaning:
    velx = tusp(n,nmx)*rho_inv ! this is actually vx/rho
    vely = tusp(n,nmy)*rho_inv ! this is actually vy/rho
    velz = tusp(n,nmz)*rho_inv ! this is actually vz/rho
    etot = tusp(n,nee)*rho_inv ! this is actually (total energy)/rho
    !
    drhodx = tdcvdx(n,1,nec)
    drhody = tdcvdx(n,2,nec)
    drhodz = tdcvdx(n,3,nec)
    !
    return_value(n,1,nmx) = rho_inv*tdcvdx(n,1,nmx) - drhodx*velx
    return_value(n,2,nmx) = rho_inv*tdcvdx(n,2,nmx) - drhody*velx
    return_value(n,3,nmx) = rho_inv*tdcvdx(n,3,nmx) - drhodz*velx
    !
    return_value(n,1,nmy) = rho_inv*tdcvdx(n,1,nmy) - drhodx*vely
    return_value(n,2,nmy) = rho_inv*tdcvdx(n,2,nmy) - drhody*vely
    return_value(n,3,nmy) = rho_inv*tdcvdx(n,3,nmy) - drhodz*vely
    !
    return_value(n,1,nmz) = rho_inv*tdcvdx(n,1,nmz) - drhodx*velz
    return_value(n,2,nmz) = rho_inv*tdcvdx(n,2,nmz) - drhody*velz
    return_value(n,3,nmz) = rho_inv*tdcvdx(n,3,nmz) - drhodz*velz
    !
    return_value(n,1,nee) = rho_inv*tdcvdx(n,1,nee) - drhodx*etot
    return_value(n,2,nee) = rho_inv*tdcvdx(n,2,nee) - drhody*etot
    return_value(n,3,nee) = rho_inv*tdcvdx(n,3,nee) - drhodz*etot
    !
   !return_value(n,l,m) = rho_inv*tdcvdx(n,l,m) - &
   !                      rho_inv*tdcvdx(n,l,nec)*tusp(n,m)
    !
  end do
  !
  ! The location for the temperature gradient currently contains
  ! the total energy gradient from the previous loop.
  ! Now use the velocity gradients to subtract the gradient of
  ! kinetic-energy-per-unit-mass from the total energy gradient.
  ! 
  !$ACC PARALLEL LOOP
  do n = n1,n2
    !
    velx = tusp(n,nmx)
    vely = tusp(n,nmy)
    velz = tusp(n,nmz)
    !
    return_value(n,1,nee) = return_value(n,1,nee) &
                          - return_value(n,1,nmx)*velx &
                          - return_value(n,1,nmy)*vely &
                          - return_value(n,1,nmz)*velz
    !
    return_value(n,2,nee) = return_value(n,2,nee) &
                          - return_value(n,2,nmx)*velx &
                          - return_value(n,2,nmy)*vely &
                          - return_value(n,2,nmz)*velz
    !
    return_value(n,3,nee) = return_value(n,3,nee) &
                          - return_value(n,3,nmx)*velx &
                          - return_value(n,3,nmy)*vely &
                          - return_value(n,3,nmz)*velz
    !
  end do
  !
  ! Finish computing the temperature gradient by dividing through by
  ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
  !
  !$ACC PARALLEL LOOP
  do n = n1,n2
    return_value(n,1,nee) = return_value(n,1,nee) * cv_inv
    return_value(n,2,nee) = return_value(n,2,nee) * cv_inv
    return_value(n,3,nee) = return_value(n,3,nee) * cv_inv
  end do
  !
  ! Remove the embedded 1/rho in the all variables except density
  !
  if (tusp_needs_rho_inv) then
    !$ACC PARALLEL LOOP
    do n = n1,n2
      rho = tusp(n,nec)
      tusp(n,nmx) = tusp(n,nmx)*rho
      tusp(n,nmy) = tusp(n,nmy)*rho
      tusp(n,nmz) = tusp(n,nmz)*rho
      tusp(n,nee) = tusp(n,nee)*rho
    end do
  end if
  !$ACC END DATA
  !
end subroutine grad_cv_to_grad_pv_gpu_3d
!
!###############################################################################
!
pure subroutine grad_cv_to_grad_pv_gpu_2d(n1,n2,tusp,tdcvdx,return_value, &
                                          tusp_includes_rho_inv)
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: tusp(:,:)
  real(wp),    intent(in) :: tdcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: tusp_includes_rho_inv
  !
  !.. Local Scalars ..
  integer  :: l,m,n
  real(wp) :: cv_inv,rho_inv
  !
  logical(lk) :: tusp_needs_rho_inv
  !
continue
  !
  tusp_needs_rho_inv = true
  if (present(tusp_includes_rho_inv)) then
    tusp_needs_rho_inv = .not. tusp_includes_rho_inv
  end if
  !
  !$ACC DATA PRESENT(tusp,tdcvdx,return_value)
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  ! Embed 1/rho in the all variables except density
  !
  if (tusp_needs_rho_inv) then
    !$ACC PARALLEL LOOP COLLAPSE(2)
    do m = nmb,nq
      do n = n1,n2
        tusp(n,m) = tusp(n,m)/tusp(n,nec)
      end do
    end do
  end if
  !
  ! Copy the density gradient
  !
  !$ACC PARALLEL LOOP COLLAPSE(2)
  do l = 1,2
    do n = n1,n2
      return_value(n,l,nec) = tdcvdx(n,l,nec)
    end do
  end do
  !
  ! Using the chain rule:
  !
  ! 1. Compute the velocity gradients (m = 2,nq-1)
  ! 2. Compute the total energy gradient and store in the
  !    location for the temperature gradient (m = nq)
  !
  !$ACC PARALLEL LOOP COLLAPSE(3)
  do m = nmb,nq
    do l = 1,2
      do n = n1,n2
        rho_inv = one/tusp(n,nec)
        return_value(n,l,m) = rho_inv*tdcvdx(n,l,m) - &
                              rho_inv*tdcvdx(n,l,nec)*tusp(n,m)
      end do
    end do
  end do
  !
  ! The location for the temperature gradient currently contains
  ! the total energy gradient from the previous loop.
  ! Now use the velocity gradients to subtract the gradient of
  ! kinetic-energy-per-unit-mass from the total energy gradient.
  ! 
  !
#ifdef NOT_SURE_ABOUT_LOOP
  do m = nmb,nme
    !$ACC PARALLEL LOOP COLLAPSE(2)
    do l = 1,2
      do n = n1,n2
        return_value(n,l,nee) = return_value(n,l,nee) &
                              - return_value(n,l,m)*tusp(n,m)
      end do
    end do
  end do
#else
  !$ACC PARALLEL LOOP COLLAPSE(2)
  do l = 1,2
    do n = n1,n2
      return_value(n,l,nee) = return_value(n,l,nee) &
                            - return_value(n,l,nmx)*tusp(n,nmx) &
                            - return_value(n,l,nmy)*tusp(n,nmy)
    end do
  end do
#endif
  !
  ! Finish computing the temperature gradient by dividing through by
  ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
  !
  !$ACC PARALLEL LOOP COLLAPSE(2)
  do l = 1,2
    do n = n1,n2
      return_value(n,l,nee) = return_value(n,l,nee) * cv_inv
    end do
  end do
  !
  ! Remove the embedded 1/rho in the all variables except density
  !
  if (tusp_needs_rho_inv) then
    !$ACC PARALLEL LOOP COLLAPSE(2)
    do m = nmb,nq
      do n = n1,n2
        tusp(n,m) = tusp(n,m)*tusp(n,nec)
      end do
    end do
  end if
  !$ACC END DATA
  !
end subroutine grad_cv_to_grad_pv_gpu_2d
!
!###############################################################################
!
pure subroutine grad_cv_to_grad_pv_gpu_nd(n1,n2,tusp,tdcvdx,return_value, &
                                          tusp_includes_rho_inv)
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: tusp(:,:)
  real(wp),    intent(in) :: tdcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: tusp_includes_rho_inv
  !
  !.. Local Scalars ..
  integer  :: l,m,n,nr
  real(wp) :: cv_inv,rho_inv
  !
  logical(lk) :: tusp_needs_rho_inv
  !
continue
  !
  tusp_needs_rho_inv = true
  if (present(tusp_includes_rho_inv)) then
    tusp_needs_rho_inv = .not. tusp_includes_rho_inv
  end if
  !
  !$ACC DATA PRESENT(tusp,tdcvdx,return_value)
  !
  nr = size(tdcvdx,dim=2)
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  ! Embed 1/rho in the all variables except density
  !
  if (tusp_needs_rho_inv) then
    !$ACC PARALLEL LOOP COLLAPSE(2)
    do m = nmb,nq
      do n = n1,n2
        tusp(n,m) = tusp(n,m)/tusp(n,nec)
      end do
    end do
  end if
  !
  ! Copy the density gradient
  !
  !$ACC PARALLEL LOOP COLLAPSE(2)
  do l = 1,nr
    do n = n1,n2
      return_value(n,l,nec) = tdcvdx(n,l,nec)
    end do
  end do
  !
  ! Using the chain rule:
  !
  ! 1. Compute the velocity gradients (m = 2,nq-1)
  ! 2. Compute the total energy gradient and store in the
  !    location for the temperature gradient (m = nq)
  !
  !$ACC PARALLEL LOOP COLLAPSE(3)
  do m = nmb,nq
    do l = 1,nr
      do n = n1,n2
        rho_inv = one/tusp(n,nec)
        return_value(n,l,m) = rho_inv*tdcvdx(n,l,m) - &
                              rho_inv*tdcvdx(n,l,nec)*tusp(n,m)
      end do
    end do
  end do
  !
  ! The location for the temperature gradient currently contains
  ! the total energy gradient from the previous loop.
  ! Now use the velocity gradients to subtract the gradient of
  ! kinetic-energy-per-unit-mass from the total energy gradient.
  ! 
  !
#ifdef NOT_SURE_ABOUT_LOOP
  do m = nmb,nme
    !$ACC PARALLEL LOOP COLLAPSE(2)
    do l = 1,nr
      do n = n1,n2
        return_value(n,l,nee) = return_value(n,l,nee) &
                              - return_value(n,l,m)*tusp(n,m)
      end do
    end do
  end do
#else
  if (nr == 3) then
    !$ACC PARALLEL LOOP COLLAPSE(2)
    do l = 1,3
      do n = n1,n2
        return_value(n,l,nee) = return_value(n,l,nee) &
                              - return_value(n,l,nmx)*tusp(n,nmx) &
                              - return_value(n,l,nmy)*tusp(n,nmy) &
                              - return_value(n,l,nmz)*tusp(n,nmz)
      end do
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP COLLAPSE(2)
    do l = 1,2
      do n = n1,n2
        return_value(n,l,nee) = return_value(n,l,nee) &
                              - return_value(n,l,nmx)*tusp(n,nmx) &
                              - return_value(n,l,nmy)*tusp(n,nmy)
      end do
    end do
  else
    !$ACC PARALLEL LOOP
    do n = n1,n2
      return_value(n,1,nee) = return_value(n,1,nee) &
                            - return_value(n,1,nmx)*tusp(n,nmx)
    end do
  end if
#endif
  !
  ! Finish computing the temperature gradient by dividing through by
  ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
  !
  !$ACC PARALLEL LOOP COLLAPSE(2)
  do l = 1,nr
    do n = n1,n2
      return_value(n,l,nee) = return_value(n,l,nee) * cv_inv
    end do
  end do
  !
  ! Remove the embedded 1/rho in the all variables except density
  !
  if (tusp_needs_rho_inv) then
    !$ACC PARALLEL LOOP COLLAPSE(2)
    do m = nmb,nq
      do n = n1,n2
        tusp(n,m) = tusp(n,m)*tusp(n,nec)
      end do
    end do
  end if
  !$ACC END DATA
  !
end subroutine grad_cv_to_grad_pv_gpu_nd
!
!###############################################################################
!
end module gpu_mod
