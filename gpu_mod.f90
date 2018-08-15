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
  public :: usp2v_gpu
  public :: tusp2v_gpu
  public :: vsp2u_gpu
  public :: tvsp2u_gpu
  public :: grad_cv_to_grad_pv_gpu
  public :: vortex_solution_gpu
  !
  interface tvsp2u_gpu
    module procedure tvsp2u_gpu_rv, tvsp2u_gpu_inout
  end interface tvsp2u_gpu
  !
contains
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(usp,return_value)
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
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
    !$ACC&              PRESENT(usp,return_value)
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
      !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
      !$ACC&              PRESENT(return_value)
      do n = n1,n2
        rho = return_value(nec,n)
        pres = return_value(nee,n)
        return_value(nec,n) = pres * rgasref_nd_inv / rho
      end do
    end if
  end if
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(usp,return_value)
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
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
    !$ACC&              PRESENT(usp,return_value)
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
      !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
      !$ACC&              PRESENT(return_value)
      do n = n1,n2
        rho = return_value(nec,n)
        pres = return_value(nee,n)
        return_value(nec,n) = pres * rgasref_nd_inv / rho
      end do
    end if
  end if
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(usp,return_value)
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
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
    !$ACC&              PRESENT(usp,return_value)
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
      !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
      !$ACC&              PRESENT(return_value)
      do n = n1,n2
        rho = return_value(nec,n)
        pres = return_value(nee,n)
        return_value(nec,n) = pres * rgasref_nd_inv / rho
      end do
    end if
  end if
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
  !$ACC PARALLEL LOOP PRESENT(tusp,return_value)
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
    !$ACC PARALLEL LOOP PRESENT(tusp,return_value)
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
      !$ACC PARALLEL LOOP PRESENT(return_value)
      do n = n1,n2
        rho = return_value(n,nec)
        pres = return_value(n,nee)
        return_value(n,nec) = pres * rgasref_nd_inv / rho
      end do
    end if
  end if
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
  !$ACC PARALLEL LOOP PRESENT(tusp,return_value)
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
    !$ACC PARALLEL LOOP PRESENT(tusp,return_value)
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
      !$ACC PARALLEL LOOP PRESENT(return_value)
      do n = n1,n2
        rho = return_value(n,nec)
        pres = return_value(n,nee)
        return_value(n,nec) = pres * rgasref_nd_inv / rho
      end do
    end if
  end if
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
  !$ACC PARALLEL LOOP PRESENT(tusp,return_value)
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
    !$ACC PARALLEL LOOP PRESENT(tusp,return_value)
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
      !$ACC PARALLEL LOOP PRESENT(return_value)
      do n = n1,n2
        rho = return_value(n,nec)
        pres = return_value(n,nee)
        return_value(n,nec) = pres * rgasref_nd_inv / rho
      end do
    end if
  end if
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(vsp,return_value)
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
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
    !$ACC&              PRESENT(vsp,return_value)
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(vsp,return_value)
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
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
    !$ACC&              PRESENT(vsp,return_value)
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(vsp,return_value)
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
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
    !$ACC&              PRESENT(vsp,return_value)
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
end subroutine vsp2u_gpu_nd
!
!###############################################################################
!
pure subroutine tvsp2u_gpu_rv(n1,n2,tvsp,return_value)
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
    call tvsp2u_gpu_rv_3d(n1,n2,tvsp,return_value)
  else if (nr == 2) then
    call tvsp2u_gpu_rv_2d(n1,n2,tvsp,return_value)
  else
    call tvsp2u_gpu_rv_nd(n1,n2,tvsp,return_value)
  end if
  !
end subroutine tvsp2u_gpu_rv
!
!###############################################################################
!
pure subroutine tvsp2u_gpu_rv_3d(n1,n2,tvsp,return_value)
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
  !$ACC PARALLEL LOOP PRESENT(tvsp,return_value)
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
    !$ACC PARALLEL LOOP PRESENT(tvsp,return_value)
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
end subroutine tvsp2u_gpu_rv_3d
!
!###############################################################################
!
pure subroutine tvsp2u_gpu_rv_2d(n1,n2,tvsp,return_value)
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
  !$ACC PARALLEL LOOP PRESENT(tvsp,return_value)
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
    !$ACC PARALLEL LOOP PRESENT(tvsp,return_value)
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
end subroutine tvsp2u_gpu_rv_2d
!
!###############################################################################
!
pure subroutine tvsp2u_gpu_rv_nd(n1,n2,tvsp,return_value)
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
  !$ACC PARALLEL LOOP PRESENT(tvsp,return_value)
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
    !$ACC PARALLEL LOOP PRESENT(tvsp,return_value)
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
end subroutine tvsp2u_gpu_rv_nd
!
!###############################################################################
!
pure subroutine tvsp2u_gpu_inout(n1,n2,tvsp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: tvsp(:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call tvsp2u_gpu_inout_3d(n1,n2,tvsp)
  else if (nr == 2) then
    call tvsp2u_gpu_inout_2d(n1,n2,tvsp)
  else
    call tvsp2u_gpu_inout_nd(n1,n2,tvsp)
  end if
  !
end subroutine tvsp2u_gpu_inout
!
!###############################################################################
!
pure subroutine tvsp2u_gpu_inout_3d(n1,n2,tvsp)
  !
  !.. Use Statements ..
  use ovar, only : gam, include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: tvsp(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1_inv,q2,rho,rhoe,velx,vely,velz,pres,tke,tls
  !
continue
  !
  gm1_inv = one/(gam-one)
  !
  !$ACC PARALLEL LOOP PRESENT(tvsp)
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
    tvsp(n,nec) = rho
    tvsp(n,nmx) = rho*velx
    tvsp(n,nmy) = rho*vely
    tvsp(n,nmz) = rho*velz
    tvsp(n,nee) = pres*gm1_inv + half*rho*q2
    !
  end do
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP PRESENT(tvsp)
    do n = n1,n2
      rho  = tvsp(n,nec)
      rhoe = tvsp(n,nee)
      tke  = tvsp(n,ntk)
      tls  = tvsp(n,ntl)
      tvsp(n,ntk) = rho*tke
      tvsp(n,ntl) = rho*tls
      tvsp(n,nee) = rhoe + rho*tke
    end do
  end if
  !
end subroutine tvsp2u_gpu_inout_3d
!
!###############################################################################
!
pure subroutine tvsp2u_gpu_inout_2d(n1,n2,tvsp)
  !
  !.. Use Statements ..
  use ovar, only : gam, include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: tvsp(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1_inv,q2,rho,rhoe,velx,vely,pres,tke,tls
  !
continue
  !
  gm1_inv = one/(gam-one)
  !
  !$ACC PARALLEL LOOP PRESENT(tvsp)
  do n = n1,n2
    !
    rho = tvsp(n,nec)
    velx = tvsp(n,nmx)
    vely = tvsp(n,nmy)
    pres = tvsp(n,nee)
    !
    q2 = velx*velx + vely*vely
    !
    tvsp(n,nec) = rho
    tvsp(n,nmx) = rho*velx
    tvsp(n,nmy) = rho*vely
    tvsp(n,nee) = pres*gm1_inv + half*rho*q2
    !
  end do
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP PRESENT(tvsp)
    do n = n1,n2
      rho  = tvsp(n,nec)
      rhoe = tvsp(n,nee)
      tke  = tvsp(n,ntk)
      tls  = tvsp(n,ntl)
      tvsp(n,ntk) = rho*tke
      tvsp(n,ntl) = rho*tls
      tvsp(n,nee) = rhoe + rho*tke
    end do
  end if
  !
end subroutine tvsp2u_gpu_inout_2d
!
!###############################################################################
!
pure subroutine tvsp2u_gpu_inout_nd(n1,n2,tvsp)
  !
  !.. Use Statements ..
  use ovar, only : gam, include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: tvsp(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: gm1_inv,q2,rho,rhoe,velx,vely,velz,pres,tke,tls
  !
continue
  !
  nr = nme-nmb+1
  !
  gm1_inv = one/(gam-one)
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(tvsp)
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
     !tvsp(n,nec) = rho
      tvsp(n,nmx) = rho*velx
      tvsp(n,nmy) = rho*vely
      tvsp(n,nmz) = rho*velz
      tvsp(n,nee) = pres*gm1_inv + half*rho*q2
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(tvsp)
    do n = n1,n2
      !
      rho = tvsp(n,nec)
      velx = tvsp(n,nmx)
      vely = tvsp(n,nmy)
      pres = tvsp(n,nee)
      !
      q2 = velx*velx + vely*vely
      !
     !tvsp(n,nec) = rho
      tvsp(n,nmx) = rho*velx
      tvsp(n,nmy) = rho*vely
      tvsp(n,nee) = pres*gm1_inv + half*rho*q2
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(tvsp)
    do n = n1,n2
      !
      rho = tvsp(n,nec)
      velx = tvsp(n,nmx)
      pres = tvsp(n,nee)
      !
     !tvsp(n,nec) = rho
      tvsp(n,nmx) = rho*velx
      tvsp(n,nee) = pres*gm1_inv + half*rho*velx*velx
      !
    end do
  end if
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP PRESENT(tvsp)
    do n = n1,n2
      rho  = tvsp(n,nec)
      rhoe = tvsp(n,nee)
      tke  = tvsp(n,ntk)
      tls  = tvsp(n,ntl)
      tvsp(n,ntk) = rho*tke
      tvsp(n,ntl) = rho*tls
      tvsp(n,nee) = rhoe + rho*tke
    end do
  end if
  !
end subroutine tvsp2u_gpu_inout_nd
!
!###############################################################################
!
pure subroutine grad_cv_to_grad_pv_gpu(n1,n2,tusp,tdcvdx,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: tusp(:,:)
  real(wp),    intent(in) :: tdcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: nr
  !
continue
  !
  nr = size(tdcvdx,dim=2)
  !
  if (nr == 3) then
    call grad_cv_to_grad_pv_gpu_3d(n1,n2,tusp,tdcvdx,return_value)
  else if (nr == 2) then
    call grad_cv_to_grad_pv_gpu_2d(n1,n2,tusp,tdcvdx,return_value)
  else
    call grad_cv_to_grad_pv_gpu_nd(n1,n2,tusp,tdcvdx,return_value)
  end if
  !
end subroutine grad_cv_to_grad_pv_gpu
!
!###############################################################################
!
pure subroutine grad_cv_to_grad_pv_gpu_3d(n1,n2,tusp,tdcvdx,return_value)
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: tusp(:,:)
  real(wp),    intent(in) :: tdcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: cv_inv,rho_inv
  real(wp) :: decdx,decdy,decdz
  real(wp) :: rho2_inv,a
  !
continue
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  ! Using the chain rule:
  !
  ! 1. Copy the density gradient
  ! 2. Compute the velocity gradients (nmx,nmy,nmz)
  ! 3. Compute the total energy gradient and store in the
  !    location for the temperature gradient (nee)
  !
  !$ACC PARALLEL LOOP PRESENT(tusp,tdcvdx,return_value)
  do n = n1,n2
    !
    rho_inv = one/tusp(n,nec)
    !
    decdx = tdcvdx(n,1,nec)
    decdy = tdcvdx(n,2,nec)
    decdz = tdcvdx(n,3,nec)
    !
    return_value(n,1,nec) = decdx
    return_value(n,2,nec) = decdy
    return_value(n,3,nec) = decdz
    !
    decdx = decdx*rho_inv*rho_inv
    decdy = decdy*rho_inv*rho_inv
    decdz = decdz*rho_inv*rho_inv
    !
    return_value(n,1,nmx) = rho_inv*tdcvdx(n,1,nmx) - decdx*tusp(n,nmx)
    return_value(n,2,nmx) = rho_inv*tdcvdx(n,2,nmx) - decdy*tusp(n,nmx)
    return_value(n,3,nmx) = rho_inv*tdcvdx(n,3,nmx) - decdz*tusp(n,nmx)
    !
    return_value(n,1,nmy) = rho_inv*tdcvdx(n,1,nmy) - decdx*tusp(n,nmy)
    return_value(n,2,nmy) = rho_inv*tdcvdx(n,2,nmy) - decdy*tusp(n,nmy)
    return_value(n,3,nmy) = rho_inv*tdcvdx(n,3,nmy) - decdz*tusp(n,nmy)
    !
    return_value(n,1,nmz) = rho_inv*tdcvdx(n,1,nmz) - decdx*tusp(n,nmz)
    return_value(n,2,nmz) = rho_inv*tdcvdx(n,2,nmz) - decdy*tusp(n,nmz)
    return_value(n,3,nmz) = rho_inv*tdcvdx(n,3,nmz) - decdz*tusp(n,nmz)
    !
    return_value(n,1,nee) = rho_inv*tdcvdx(n,1,nee) - decdx*tusp(n,nee)
    return_value(n,2,nee) = rho_inv*tdcvdx(n,2,nee) - decdy*tusp(n,nee)
    return_value(n,3,nee) = rho_inv*tdcvdx(n,3,nee) - decdz*tusp(n,nee)
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
  !$ACC PARALLEL LOOP PRESENT(tusp,return_value)
  do n = n1,n2
    !
    rho_inv = one/tusp(n,nec)
    !
    return_value(n,1,nee) = return_value(n,1,nee) &
                          - return_value(n,1,nmx)*tusp(n,nmx)*rho_inv &
                          - return_value(n,1,nmy)*tusp(n,nmy)*rho_inv &
                          - return_value(n,1,nmz)*tusp(n,nmz)*rho_inv
    !
    return_value(n,2,nee) = return_value(n,2,nee) &
                          - return_value(n,2,nmx)*tusp(n,nmx)*rho_inv &
                          - return_value(n,2,nmy)*tusp(n,nmy)*rho_inv &
                          - return_value(n,2,nmz)*tusp(n,nmz)*rho_inv
    !
    return_value(n,3,nee) = return_value(n,3,nee) &
                          - return_value(n,3,nmx)*tusp(n,nmx)*rho_inv &
                          - return_value(n,3,nmy)*tusp(n,nmy)*rho_inv &
                          - return_value(n,3,nmz)*tusp(n,nmz)*rho_inv
    !
  end do
  !
  ! Finish computing the temperature gradient by dividing through by
  ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
  !
  !$ACC PARALLEL LOOP PRESENT(return_value)
  do n = n1,n2
    return_value(n,1,nee) = return_value(n,1,nee) * cv_inv
    return_value(n,2,nee) = return_value(n,2,nee) * cv_inv
    return_value(n,3,nee) = return_value(n,3,nee) * cv_inv
  end do
  !
end subroutine grad_cv_to_grad_pv_gpu_3d
!
!###############################################################################
!
pure subroutine grad_cv_to_grad_pv_gpu_2d(n1,n2,tusp,tdcvdx,return_value)
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: tusp(:,:)
  real(wp),    intent(in) :: tdcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: cv_inv,rho_inv
  !
continue
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  ! Copy the density gradient
  !
  !$ACC PARALLEL LOOP PRESENT(tdcvdx,return_value)
  do n = n1,n2
    return_value(n,1,nec) = tdcvdx(n,1,nec)
    return_value(n,2,nec) = tdcvdx(n,2,nec)
  end do
  !
  ! Using the chain rule:
  !
  ! 1. Compute the velocity gradients (nmx,nmy)
  ! 2. Compute the total energy gradient and store in the
  !    location for the temperature gradient (nee)
  !
  !$ACC PARALLEL LOOP PRESENT(tusp,tdcvdx,return_value)
  do n = n1,n2
    !
    rho_inv = one/tusp(n,nec)
    !
    return_value(n,1,nmx) = rho_inv*tdcvdx(n,1,nmx) &
                          - rho_inv*tdcvdx(n,1,nec)*tusp(n,nmx)*rho_inv
    return_value(n,2,nmx) = rho_inv*tdcvdx(n,2,nmx) &
                          - rho_inv*tdcvdx(n,2,nec)*tusp(n,nmx)*rho_inv
    !
    return_value(n,1,nmy) = rho_inv*tdcvdx(n,1,nmy) &
                          - rho_inv*tdcvdx(n,1,nec)*tusp(n,nmy)*rho_inv
    return_value(n,2,nmy) = rho_inv*tdcvdx(n,2,nmy) &
                          - rho_inv*tdcvdx(n,2,nec)*tusp(n,nmy)*rho_inv
    !
    return_value(n,1,nee) = rho_inv*tdcvdx(n,1,nee) &
                          - rho_inv*tdcvdx(n,1,nec)*tusp(n,nee)*rho_inv
    return_value(n,2,nee) = rho_inv*tdcvdx(n,2,nee) &
                          - rho_inv*tdcvdx(n,2,nec)*tusp(n,nee)*rho_inv
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
  !$ACC PARALLEL LOOP PRESENT(tusp,return_value)
  do n = n1,n2
    !
    rho_inv = one/tusp(n,nec)
    !
    return_value(n,1,nee) = return_value(n,1,nee) &
                          - return_value(n,1,nmx)*tusp(n,nmx)*rho_inv &
                          - return_value(n,1,nmy)*tusp(n,nmy)*rho_inv
    !
    return_value(n,2,nee) = return_value(n,2,nee) &
                          - return_value(n,2,nmx)*tusp(n,nmx)*rho_inv &
                          - return_value(n,2,nmy)*tusp(n,nmy)*rho_inv
    !
  end do
  !
  ! Finish computing the temperature gradient by dividing through by
  ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
  !
  !$ACC PARALLEL LOOP PRESENT(return_value)
  do n = n1,n2
    return_value(n,1,nee) = return_value(n,1,nee) * cv_inv
    return_value(n,2,nee) = return_value(n,2,nee) * cv_inv
  end do
  !
end subroutine grad_cv_to_grad_pv_gpu_2d
!
!###############################################################################
!
pure subroutine grad_cv_to_grad_pv_gpu_nd(n1,n2,tusp,tdcvdx,return_value)
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: tusp(:,:)
  real(wp),    intent(in) :: tdcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: cv_inv,rho_inv
  !
continue
  !
  nr = size(tdcvdx,dim=2)
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  ! Copy the density gradient
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(tdcvdx,return_value)
    do n = n1,n2
      return_value(n,1,nec) = tdcvdx(n,1,nec)
      return_value(n,2,nec) = tdcvdx(n,2,nec)
      return_value(n,3,nec) = tdcvdx(n,3,nec)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(tdcvdx,return_value)
    do n = n1,n2
      return_value(n,1,nec) = tdcvdx(n,1,nec)
      return_value(n,2,nec) = tdcvdx(n,2,nec)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(tdcvdx,return_value)
    do n = n1,n2
      return_value(n,1,nec) = tdcvdx(n,1,nec)
    end do
  end if
  !
  ! Using the chain rule:
  !
  ! 1. Compute the velocity gradients (m = 2,nq-1)
  ! 2. Compute the total energy gradient and store in the
  !    location for the temperature gradient (m = nq)
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(tusp,tdcvdx,return_value)
    do n = n1,n2
      !
      rho_inv = one/tusp(n,nec)
      !
      return_value(n,1,nmx) = rho_inv*tdcvdx(n,1,nmx) &
                            - rho_inv*tdcvdx(n,1,nec)*tusp(n,nmx)*rho_inv
      return_value(n,2,nmx) = rho_inv*tdcvdx(n,2,nmx) &
                            - rho_inv*tdcvdx(n,2,nec)*tusp(n,nmx)*rho_inv
      return_value(n,3,nmx) = rho_inv*tdcvdx(n,3,nmx) &
                            - rho_inv*tdcvdx(n,3,nec)*tusp(n,nmx)*rho_inv
      !
      return_value(n,1,nmy) = rho_inv*tdcvdx(n,1,nmy) &
                            - rho_inv*tdcvdx(n,1,nec)*tusp(n,nmy)*rho_inv
      return_value(n,2,nmy) = rho_inv*tdcvdx(n,2,nmy) &
                            - rho_inv*tdcvdx(n,2,nec)*tusp(n,nmy)*rho_inv
      return_value(n,3,nmy) = rho_inv*tdcvdx(n,3,nmy) &
                            - rho_inv*tdcvdx(n,3,nec)*tusp(n,nmy)*rho_inv
      !
      return_value(n,1,nmz) = rho_inv*tdcvdx(n,1,nmz) &
                            - rho_inv*tdcvdx(n,1,nec)*tusp(n,nmz)*rho_inv
      return_value(n,2,nmz) = rho_inv*tdcvdx(n,2,nmz) &
                            - rho_inv*tdcvdx(n,2,nec)*tusp(n,nmz)*rho_inv
      return_value(n,3,nmz) = rho_inv*tdcvdx(n,3,nmz) &
                            - rho_inv*tdcvdx(n,3,nec)*tusp(n,nmz)*rho_inv
      !
      return_value(n,1,nee) = rho_inv*tdcvdx(n,1,nee) &
                            - rho_inv*tdcvdx(n,1,nec)*tusp(n,nee)*rho_inv
      return_value(n,2,nee) = rho_inv*tdcvdx(n,2,nee) &
                            - rho_inv*tdcvdx(n,2,nec)*tusp(n,nee)*rho_inv
      return_value(n,3,nee) = rho_inv*tdcvdx(n,3,nee) &
                            - rho_inv*tdcvdx(n,3,nec)*tusp(n,nee)*rho_inv
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(tusp,tdcvdx,return_value)
    do n = n1,n2
      !
      rho_inv = one/tusp(n,nec)
      !
      return_value(n,1,nmx) = rho_inv*tdcvdx(n,1,nmx) &
                            - rho_inv*tdcvdx(n,1,nec)*tusp(n,nmx)*rho_inv
      return_value(n,2,nmx) = rho_inv*tdcvdx(n,2,nmx) &
                            - rho_inv*tdcvdx(n,2,nec)*tusp(n,nmx)*rho_inv
      !
      return_value(n,1,nmy) = rho_inv*tdcvdx(n,1,nmy) &
                            - rho_inv*tdcvdx(n,1,nec)*tusp(n,nmy)*rho_inv
      return_value(n,2,nmy) = rho_inv*tdcvdx(n,2,nmy) &
                            - rho_inv*tdcvdx(n,2,nec)*tusp(n,nmy)*rho_inv
      !
      return_value(n,1,nee) = rho_inv*tdcvdx(n,1,nee) &
                            - rho_inv*tdcvdx(n,1,nec)*tusp(n,nee)*rho_inv
      return_value(n,2,nee) = rho_inv*tdcvdx(n,2,nee) &
                            - rho_inv*tdcvdx(n,2,nec)*tusp(n,nee)*rho_inv
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(tusp,tdcvdx,return_value)
    do n = n1,n2
      !
      rho_inv = one/tusp(n,nec)
      !
      return_value(n,1,nmx) = rho_inv*tdcvdx(n,1,nmx) &
                            - rho_inv*tdcvdx(n,1,nec)*tusp(n,nmx)*rho_inv
      !
      return_value(n,1,nee) = rho_inv*tdcvdx(n,1,nee) &
                            - rho_inv*tdcvdx(n,1,nec)*tusp(n,nee)*rho_inv
      !
    end do
  end if
  !
  ! The location for the temperature gradient currently contains
  ! the total energy gradient from the previous loop.
  ! Now use the velocity gradients to subtract the gradient of
  ! kinetic-energy-per-unit-mass from the total energy gradient.
  !
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(tusp,return_value)
    do n = n1,n2
      !
      rho_inv = one/tusp(n,nec)
      !
      return_value(n,1,nee) = return_value(n,1,nee) &
                            - return_value(n,1,nmx)*tusp(n,nmx)*rho_inv &
                            - return_value(n,1,nmy)*tusp(n,nmy)*rho_inv &
                            - return_value(n,1,nmz)*tusp(n,nmz)*rho_inv
      !
      return_value(n,2,nee) = return_value(n,2,nee) &
                            - return_value(n,2,nmx)*tusp(n,nmx)*rho_inv &
                            - return_value(n,2,nmy)*tusp(n,nmy)*rho_inv &
                            - return_value(n,2,nmz)*tusp(n,nmz)*rho_inv
      !
      return_value(n,3,nee) = return_value(n,3,nee) &
                            - return_value(n,3,nmx)*tusp(n,nmx)*rho_inv &
                            - return_value(n,3,nmy)*tusp(n,nmy)*rho_inv &
                            - return_value(n,3,nmz)*tusp(n,nmz)*rho_inv
      !
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(tusp,return_value)
    do n = n1,n2
      !
      rho_inv = one/tusp(n,nec)
      !
      return_value(n,1,nee) = return_value(n,1,nee) &
                            - return_value(n,1,nmx)*tusp(n,nmx)*rho_inv &
                            - return_value(n,1,nmy)*tusp(n,nmy)*rho_inv
      !
      return_value(n,2,nee) = return_value(n,2,nee) &
                            - return_value(n,2,nmx)*tusp(n,nmx)*rho_inv &
                            - return_value(n,2,nmy)*tusp(n,nmy)*rho_inv
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(tusp,return_value)
    do n = n1,n2
      return_value(n,1,nee) = return_value(n,1,nee) &
                            - return_value(n,1,nmx)*tusp(n,nmx)/tusp(n,nec)
    end do
  end if
  !
  ! Finish computing the temperature gradient by dividing through by
  ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(return_value)
    do n = n1,n2
      return_value(n,1,nee) = return_value(n,1,nee) * cv_inv
      return_value(n,2,nee) = return_value(n,2,nee) * cv_inv
      return_value(n,3,nee) = return_value(n,3,nee) * cv_inv
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(return_value)
    do n = n1,n2
      return_value(n,1,nee) = return_value(n,1,nee) * cv_inv
      return_value(n,2,nee) = return_value(n,2,nee) * cv_inv
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(return_value)
    do n = n1,n2
      return_value(n,1,nee) = return_value(n,1,nee) * cv_inv
    end do
  end if
  !
end subroutine grad_cv_to_grad_pv_gpu_nd
!
!###############################################################################
!
pure subroutine vortex_solution_gpu(n1,n2,xyz,return_value)
  !
  !.. Use Statements ..
  use ovar, only : velref,aref,uref,machref
  use ovar, only : gam,rhoref,tref,pref,itestcase
  use ovar, only : vortex_bigR,vortex_bigGam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: xyz(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: c0,c1,c2,c3,c4,c5,bigR,bigGam
  real(wp) :: gm1,gm1_inv,gam_over_gm1
  real(wp) :: rhoref_inv,aref_inv,rhoa2_inv,tref_inv
  real(wp) :: vx_ref,vy_ref
  real(wp) :: x,y,a,du,dv,dT,vx,vy,vz,temp,rho,pres
  !
  logical(lk) :: running_vortex_problem
  !
continue
  !
  nr = size(xyz,dim=1)
  !
  ! Some needed constants
  !
  running_vortex_problem = any(itestcase==Vortex_Problems)
  !
  c0 = deltafun( running_vortex_problem )
  !
  gm1 = gam - one
  gm1_inv = one/gm1
  ! NOTE: Multiply by c0 because -> If doing density/entropy transport, make the
  !       pressure constant everywhere and set it to the free-stream pressure
  gam_over_gm1 = gam*gm1_inv*c0
  !
  rhoref_inv = one/rhoref
  aref_inv = one/aref
  rhoa2_inv = rhoref_inv*aref_inv*aref_inv
  tref_inv = one/tref
  !
  vx_ref = velref(1)
  vy_ref = velref(2)
  !
  ! For density/entropy transport or Shu vortex -> bigR = 1;   bigGam = 5
  ! For vortex          transport               -> bigR = 1.5; bigGam = 13.5
  !
  bigR   = merge(vortex_bigR,  one, running_vortex_problem)
  bigGam = merge(vortex_bigGam,five,running_vortex_problem)
  !
  c1 = half/(bigR*bigR)
  c2 = half*bigGam/pi
  c3 = half*machref*machref*gm1*tref
  c4 = uref/bigR
  !
  ! If doing density transport, make c5 = 0 to zero out pressure part of energy
  !
  c5 = gm1_inv * (one - deltafun(itestcase==Density_Transport))
  !
  !$ACC PARALLEL LOOP PRESENT(xyz,return_value)
  do n = n1,n2
    !
    x = xyz(1,n)
    y = xyz(2,n)
    !
    a = c2*exp((one - x*x - y*y)*c1)
    !
    du = -c4*a*y
    dv =  c4*a*x
    dT = -c3*a*a
    !
    ! Only add on the velocity pertubations if doing vortex transport
    ! i.e., the velocity is constant if doing density/entropy transport
    !
    vx = vx_ref + du*c0
    vy = vy_ref + dv*c0
    temp = tref  + dT
    !
    rho  = rhoref  * (temp*tref_inv)**(gm1_inv)
    pres = c5*pref * (temp*tref_inv)**(gam_over_gm1)
    !
    ! Evaluate the dimensional initial conditions using conserved variables
    !
    return_value(n,nec) = rho * rhoref_inv
    return_value(n,nmx) = vx * aref_inv
    return_value(n,nmy) = vy * aref_inv
    return_value(n,nee) = pres * rhoa2_inv
    !
  end do
  !
  if (nr == 3) then
    vz = velref(3)*aref_inv
    !$ACC PARALLEL LOOP PRESENT(return_value)
    do n = n1,n2
      return_value(n,nmz) = vz
    end do
  end if
  !
  call tvsp2u_gpu(n1,n2,return_value)
  !
end subroutine vortex_solution_gpu
!
!###############################################################################
!
pure subroutine vortex_offset_gpu(xyz,return_value)
  !
  !.. Use Statements ..
  use ovar, only : itcur,velref,aref
  use ovar, only : Period_Timesteps,Period_Time
  use ovar, only : Grid_Length,Grid_Min
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: xyz(:,:)
  real(wp), intent(inout) :: return_value(:,:)
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
end subroutine vortex_offset_gpu
!
!###############################################################################
!
end module gpu_mod
