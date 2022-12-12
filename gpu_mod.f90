module gpu_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,ntk,ntl,nq
  use eqn_idx, only : ncs,npres,ntemp,ngamm,nenth,ncp
  !
  implicit none
  !
  private
  !
  public :: cv2pv_gpu
  public :: pv2cv_gpu
  public :: grad_cv_to_grad_pv_gpu
  public :: vortex_solution_gpu
  !
  interface grad_cv_to_grad_pv_gpu
    module procedure grad_cv_to_grad_pv_gpu, &
                     grad_cv_to_grad_pv_sp
  end interface grad_cv_to_grad_pv_gpu
  !
contains
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine cv2pv_gpu(n1,n2,cv,gp,return_value,swap_t_for_rho)
  !
  !.. Use Statements ..
  use ovar, only : turb_is_on
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: swap_t_for_rho
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: rho_inv
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call cv2pv_gpu_3d(n1,n2,cv,gp,return_value)
  else if (nr == 2) then
    call cv2pv_gpu_2d(n1,n2,cv,gp,return_value)
  else
    call cv2pv_gpu_1d(n1,n2,cv,gp,return_value)
  end if
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(cv,return_value)
    do m = 1,ncs
      do n = n1,n2
        return_value(n,m) = cv(n,m)/cv(n,nec)
      end do
    end do
  end if
  !
  if (turb_is_on) then
    !$ACC PARALLEL LOOP PRESENT(cv,return_value)
    do n = n1,n2
      rho_inv = one/cv(n,nec)
      return_value(n,ntk) = cv(n,ntk)*rho_inv
      return_value(n,ntl) = cv(n,ntl)*rho_inv
    end do
  end if
  !
  if (present(swap_t_for_rho)) then
    if (swap_t_for_rho) then
      !$ACC PARALLEL LOOP PRESENT(gp,return_value)
      do n = n1,n2
        return_value(n,nec) = gp(n,ntemp)
      end do
    end if
  end if
  !
end subroutine cv2pv_gpu
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine cv2pv_gpu_3d(n1,n2,cv,gp,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: rho_inv
  !
continue
  !
  !$ACC PARALLEL LOOP PRESENT(cv,gp,return_value)
  do n = n1,n2
    !
    rho_inv = one/cv(n,nec)
    !
    return_value(n,nec) = cv(n,nec)
    return_value(n,nmx) = cv(n,nmx)*rho_inv
    return_value(n,nmy) = cv(n,nmy)*rho_inv
    return_value(n,nmz) = cv(n,nmz)*rho_inv
    return_value(n,nee) = gp(n,npres)
    !
  end do
  !
end subroutine cv2pv_gpu_3d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine cv2pv_gpu_2d(n1,n2,cv,gp,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: rho_inv
  !
continue
  !
  !$ACC PARALLEL LOOP PRESENT(cv,gp,return_value)
  do n = n1,n2
    !
    rho_inv = one/cv(n,nec)
    !
    return_value(n,nec) = cv(n,nec)
    return_value(n,nmx) = cv(n,nmx)*rho_inv
    return_value(n,nmy) = cv(n,nmy)*rho_inv
    return_value(n,nee) = gp(n,npres)
    !
  end do
  !
end subroutine cv2pv_gpu_2d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine cv2pv_gpu_1d(n1,n2,cv,gp,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: rho_inv
  !
continue
  !
  !$ACC PARALLEL LOOP PRESENT(cv,gp,return_value)
  do n = n1,n2
    !
    rho_inv = one/cv(n,nec)
    !
    return_value(n,nec) = cv(n,nec)
    return_value(n,nmx) = cv(n,nmx)*rho_inv
    return_value(n,nee) = gp(n,npres)
    !
  end do
  !
end subroutine cv2pv_gpu_1d
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine pv2cv_gpu(n1,n2,pv,gp)
  !
  !.. Use Statements ..
  use ovar, only : include_tke_in_energy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: pv(:,:)
  real(wp), intent(inout) :: gp(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: rho,rhoe,tke,tls
  !
continue
  !
  if (ncs > 0) then
    !
    ! Convert the primitive variables for a thermally
    ! perfect gas mixture to conserved variables.
    !
    call pv2cv_gpu_tpg(n1,n2,pv,gp)
    !
  else
    !
    ! Convert the primitive variables for a calorically
    ! perfect gas to conserved variables.
    !
    call pv2cv_gpu_cpg(n1,n2,pv)
    !
  end if
  !
  ! Convert the turbulence primitive variables to conserved variables and
  ! add the turbulence kinetic energy contribution to the total energy.
  !
  if (include_tke_in_energy) then
    !$ACC PARALLEL LOOP PRESENT(pv)
    do n = n1,n2
      rho = pv(n,nec)
      rho = pv(n,nee)
      tke = pv(n,ntk)
      tls = pv(n,ntl)
      pv(n,nee) = rhoe + rho*tke
      pv(n,ntk) = rho*tke
      pv(n,ntl) = rho*tls
    end do
  end if
  !
end subroutine pv2cv_gpu
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine pv2cv_gpu_tpg(n1,n2,pv,gp)
  !
  !.. Use Statements ..
  use thermo_mod, only : compute_mixture_rgas
  use thermo_mod, only : compute_mixture_enthalpy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: pv(:,:)
  real(wp), intent(inout) :: gp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
continue
  !
  ! For a gas mixture, we have everything we need to compute the
  ! conserved variables except for total energy which requires
  ! the static enthalpy. To compute the static enthalpy, we need to:
  !   1) Compute mixture gas constant from species mass fractions.
  !   2) Compute static temperature from static density,
  !      static pressure, and mixture gas constant.
  !   3) Compute static enthalpy from static temperature
  !      and species mass fractions.
  !
  ! Compute the new mixture gas constant using
  ! the updated species mass fractions.
  !
  call compute_mixture_rgas(n1,n2,pv(:,1:ncs),gp(:,ngamm))
  !
  ! Compute the new static temperature using the updated
  ! static density, static pressure, and mixture gas constant.
  !
  !$ACC PARALLEL LOOP PRESENT(pv,gp)
  do n = n1,n2
    gp(n,ntemp) = pv(n,nee) / (pv(n,nec)*gp(n,ngamm))
  end do
  !
  ! Compute the new static enthalpy using the updated
  ! mass fractions and static temperature.
  !
  call compute_mixture_enthalpy(n1,n2,pv(:,1:ncs),gp(:,ntemp),gp(:,nenth))
  !
  ! Convert the updated mass fractions to species densities.
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(pv)
  do m = 1,ncs
    do n = n1,n2
      pv(n,m) = pv(n,m)*pv(n,nec)
    end do
  end do
  !
  ! Finally, convert the primitive variables to conserved variables,
  ! excluding the species equations which were converted to species
  ! densities in the above loop.
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call pv2cv_gpu_tpg_3d(n1,n2,pv,gp)
  else if (nr == 2) then
    call pv2cv_gpu_tpg_2d(n1,n2,pv,gp)
  else
    call pv2cv_gpu_tpg_1d(n1,n2,pv,gp)
  end if
  !
end subroutine pv2cv_gpu_tpg
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine pv2cv_gpu_tpg_3d(n1,n2,pv,gp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: pv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: rho,vx,vy,vz,q2,pres,enth
  !
continue
  !
  ! Convert the primitive variables to conserved variables.
  !
  !$ACC PARALLEL LOOP PRESENT(pv,gp)
  do n = n1,n2
    !
    rho  = pv(n,nec)
    vx   = pv(n,nmx)
    vy   = pv(n,nmy)
    vz   = pv(n,nmz)
    pres = pv(n,nee)
    !
    enth = gp(n,nenth)
    !
    q2 = vx*vx + vy*vy + vz*vz
    !
   !pv(n,nec) = rho ! NOT NEEDED, WASTE OF FLOATING POINT OPERATIONS
    pv(n,nmx) = rho*vx
    pv(n,nmy) = rho*vy
    pv(n,nmz) = rho*vz
    pv(n,nee) = rho*enth - pres + half*rho*q2
    !
  end do
  !
end subroutine pv2cv_gpu_tpg_3d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine pv2cv_gpu_tpg_2d(n1,n2,pv,gp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: pv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: rho,vx,vy,q2,pres,enth
  !
continue
  !
  ! Convert the primitive variables to conserved variables.
  !
  !$ACC PARALLEL LOOP PRESENT(pv,gp)
  do n = n1,n2
    !
    rho  = pv(n,nec)
    vx   = pv(n,nmx)
    vy   = pv(n,nmy)
    pres = pv(n,nee)
    !
    enth = gp(n,nenth)
    !
    q2 = vx*vx + vy*vy
    !
   !pv(n,nec) = rho ! NOT NEEDED, WASTE OF FLOATING POINT OPERATIONS
    pv(n,nmx) = rho*vx
    pv(n,nmy) = rho*vy
    pv(n,nee) = rho*enth - pres + half*rho*q2
    !
  end do
  !
end subroutine pv2cv_gpu_tpg_2d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine pv2cv_gpu_tpg_1d(n1,n2,pv,gp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: pv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: rho,vx,pres,enth
  !
continue
  !
  ! Convert the primitive variables to conserved variables.
  !
  !$ACC PARALLEL LOOP PRESENT(pv,gp)
  do n = n1,n2
    !
    rho  = pv(n,nec)
    vx   = pv(n,nmx)
    pres = pv(n,nee)
    !
    enth = gp(n,nenth)
    !
   !pv(n,nec) = rho ! NOT NEEDED, WASTE OF FLOATING POINT OPERATIONS
    pv(n,nmx) = rho*vx
    pv(n,nee) = rho*enth - pres + half*rho*vx*vx
    !
  end do
  !
end subroutine pv2cv_gpu_tpg_1d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine pv2cv_gpu_cpg(n1,n2,pv)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: pv(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: rho,rhoe,tke,tls
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call pv2cv_gpu_cpg_3d(n1,n2,pv)
  else if (nr == 2) then
    call pv2cv_gpu_cpg_2d(n1,n2,pv)
  else
    call pv2cv_gpu_cpg_1d(n1,n2,pv)
  end if
  !
end subroutine pv2cv_gpu_cpg
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine pv2cv_gpu_cpg_3d(n1,n2,pv)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: pv(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1_inv,q2,rho,vx,vy,vz,pres
  !
continue
  !
  gm1_inv = one/(gam-one)
  !
  !$ACC PARALLEL LOOP PRESENT(pv)
  do n = n1,n2
    !
    rho  = pv(n,nec)
    vx   = pv(n,nmx)
    vy   = pv(n,nmy)
    vz   = pv(n,nmz)
    pres = pv(n,nee)
    !
    q2 = vx*vx + vy*vy + vz*vz
    !
   !pv(n,nec) = rho ! NOT NEEDED, WASTE OF FLOATING POINT OPERATIONS
    pv(n,nmx) = rho*vx
    pv(n,nmy) = rho*vy
    pv(n,nmz) = rho*vz
    pv(n,nee) = pres*gm1_inv + half*rho*q2
    !
  end do
  !
end subroutine pv2cv_gpu_cpg_3d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine pv2cv_gpu_cpg_2d(n1,n2,pv)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: pv(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1_inv,q2,rho,vx,vy,pres
  !
continue
  !
  gm1_inv = one/(gam-one)
  !
  !$ACC PARALLEL LOOP PRESENT(pv)
  do n = n1,n2
    !
    rho  = pv(n,nec)
    vx   = pv(n,nmx)
    vy   = pv(n,nmy)
    pres = pv(n,nee)
    !
    q2 = vx*vx + vy*vy
    !
   !pv(n,nec) = rho ! NOT NEEDED, WASTE OF FLOATING POINT OPERATIONS
    pv(n,nmx) = rho*vx
    pv(n,nmy) = rho*vy
    pv(n,nee) = pres*gm1_inv + half*rho*q2
    !
  end do
  !
end subroutine pv2cv_gpu_cpg_2d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine pv2cv_gpu_cpg_1d(n1,n2,pv)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: pv(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1_inv,rho,vx,pres
  !
continue
  !
  gm1_inv = one/(gam-one)
  !
  !$ACC PARALLEL LOOP PRESENT(pv)
  do n = n1,n2
    !
    rho  = pv(n,nec)
    vx   = pv(n,nmx)
    pres = pv(n,nee)
    !
   !pv(n,nec) = rho ! NOT NEEDED, WASTE OF FLOATING POINT OPERATIONS
    pv(n,nmx) = rho*vx
    pv(n,nee) = pres*gm1_inv + half*rho*vx*vx
    !
  end do
  !
end subroutine pv2cv_gpu_cpg_1d
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine grad_cv_to_grad_pv_sp(uin,gpin,dcvdxin,return_value)
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: uin(:)
  real(wp),    intent(in) :: gpin(:)
  real(wp),    intent(in) :: dcvdxin(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: l,m
  !
  !.. Local Arrays ..
  real(wp) :: usp(1:1,1:size(uin))
  real(wp) :: gpsp(1:1,1:size(gpin))
  real(wp) :: dcvdx(1:1,1:size(dcvdxin,dim=1),1:size(dcvdxin,dim=2))
  real(wp) :: dpvdx(1:1,1:size(dcvdxin,dim=1),1:size(dcvdxin,dim=2))
  !
continue
  !
  do m = 1,size(uin)
    usp(1,m) = uin(m)
  end do
  !
  do m = 1,size(gpin)
    gpsp(1,m) = gpin(m)
  end do
  !
  do m = 1,size(dcvdxin,dim=2)
    do l = 1,size(dcvdxin,dim=1)
      dcvdx(1,l,m) = dcvdxin(l,m)
    end do
  end do
  !
  call grad_cv_to_grad_pv_gpu(1,1,usp,gpsp,dcvdx,dpvdx)
  !
  do m = 1,size(dpvdx,dim=2)
    do l = 1,size(dpvdx,dim=1)
      return_value(l,m) = dpvdx(1,l,m)
    end do
  end do
  !
end subroutine grad_cv_to_grad_pv_sp
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine grad_cv_to_grad_pv_gpu(n1,n2,usp,gpsp,dcvdx,return_value, &
                                  output_variables)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: gpsp(:,:)
  real(wp),    intent(in) :: dcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(inout) :: output_variables(:,:)
  !
  !.. Local Scalars ..
  integer  :: nr
  !
continue
  !
  nr = size(dcvdx,dim=2)
  !
  if (nr == 3) then
    call grad_cv_to_grad_pv_cpg_3d(n1,n2,usp,dcvdx,return_value)
  else if (nr == 2) then
    call grad_cv_to_grad_pv_cpg_2d(n1,n2,usp,dcvdx,return_value)
  else
    call grad_cv_to_grad_pv_cpg_1d(n1,n2,usp,dcvdx,return_value)
  end if
  !
  if (ncs > 0) then
    if (present(output_variables)) then
      if (nr == 3) then
        call adjust_temp_grad_for_tpg_3d(n1,n2,usp,gpsp,dcvdx,return_value, &
                                         output_variables)
      else if (nr == 2) then
        call adjust_temp_grad_for_tpg_2d(n1,n2,usp,gpsp,dcvdx,return_value, &
                                         output_variables)
      else
        call adjust_temp_grad_for_tpg_1d(n1,n2,usp,gpsp,dcvdx,return_value, &
                                         output_variables)
      end if
    else
      if (nr == 3) then
        call adjust_temp_gradient_for_tpg_3d(n1,n2,usp,gpsp,dcvdx,return_value)
      else if (nr == 2) then
        call adjust_temp_gradient_for_tpg_2d(n1,n2,usp,gpsp,dcvdx,return_value)
      else
        call adjust_temp_gradient_for_tpg_1d(n1,n2,usp,gpsp,dcvdx,return_value)
      end if
    end if
  end if
  !
end subroutine grad_cv_to_grad_pv_gpu
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine grad_cv_to_grad_pv_cpg_3d(n1,n2,usp,dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: dcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: cv_inv,rho_inv
  real(wp) :: decdx,decdy,decdz
  real(wp) :: vx,vy,vz,ee
  !
continue
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  ! Using the product rule:
  !
  ! 1. Copy the density gradient
  ! 2. Compute the velocity gradients (nmx,nmy,nmz)
  ! 3. Compute the total energy gradient and store in the
  !    location for the temperature gradient (nee)
  !
  !$ACC PARALLEL LOOP PRESENT(usp,dcvdx,return_value)
  do n = n1,n2
    !
    rho_inv = one/usp(n,nec)
    !
    decdx = dcvdx(n,1,nec)
    decdy = dcvdx(n,2,nec)
    decdz = dcvdx(n,3,nec)
    !
    return_value(n,1,nec) = decdx
    return_value(n,2,nec) = decdy
    return_value(n,3,nec) = decdz
    !
    decdx = decdx*rho_inv
    decdy = decdy*rho_inv
    decdz = decdz*rho_inv
    !
    vx = usp(n,nmx)*rho_inv
    vy = usp(n,nmy)*rho_inv
    vz = usp(n,nmz)*rho_inv
    ee = usp(n,nee)*rho_inv
    !
    return_value(n,1,nmx) = rho_inv*dcvdx(n,1,nmx) - decdx*vx
    return_value(n,2,nmx) = rho_inv*dcvdx(n,2,nmx) - decdy*vx
    return_value(n,3,nmx) = rho_inv*dcvdx(n,3,nmx) - decdz*vx
    !
    return_value(n,1,nmy) = rho_inv*dcvdx(n,1,nmy) - decdx*vy
    return_value(n,2,nmy) = rho_inv*dcvdx(n,2,nmy) - decdy*vy
    return_value(n,3,nmy) = rho_inv*dcvdx(n,3,nmy) - decdz*vy
    !
    return_value(n,1,nmz) = rho_inv*dcvdx(n,1,nmz) - decdx*vz
    return_value(n,2,nmz) = rho_inv*dcvdx(n,2,nmz) - decdy*vz
    return_value(n,3,nmz) = rho_inv*dcvdx(n,3,nmz) - decdz*vz
    !
    return_value(n,1,nee) = rho_inv*dcvdx(n,1,nee) - decdx*ee
    return_value(n,2,nee) = rho_inv*dcvdx(n,2,nee) - decdy*ee
    return_value(n,3,nee) = rho_inv*dcvdx(n,3,nee) - decdz*ee
    !
    ! The location for the temperature gradient currently contains
    ! the total energy gradient from the previous loop.
    ! Now use the velocity gradients to subtract the gradient of
    ! kinetic-energy-per-unit-mass from the total energy gradient.
    !
    return_value(n,1,nee) = return_value(n,1,nee) &
                          - return_value(n,1,nmx)*vx &
                          - return_value(n,1,nmy)*vy &
                          - return_value(n,1,nmz)*vz
    !
    return_value(n,2,nee) = return_value(n,2,nee) &
                          - return_value(n,2,nmx)*vx &
                          - return_value(n,2,nmy)*vy &
                          - return_value(n,2,nmz)*vz
    !
    return_value(n,3,nee) = return_value(n,3,nee) &
                          - return_value(n,3,nmx)*vx &
                          - return_value(n,3,nmy)*vy &
                          - return_value(n,3,nmz)*vz
    !
    ! Finish computing the temperature gradient by dividing through by
    ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
    !
    return_value(n,1,nee) = return_value(n,1,nee) * cv_inv
    return_value(n,2,nee) = return_value(n,2,nee) * cv_inv
    return_value(n,3,nee) = return_value(n,3,nee) * cv_inv
    !
  end do
  !
end subroutine grad_cv_to_grad_pv_cpg_3d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine grad_cv_to_grad_pv_cpg_3d_original(n1,n2,usp,dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: dcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: cv_inv,rho_inv
  real(wp) :: decdx,decdy,decdz
  !
continue
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  ! Using the product rule:
  !
  ! 1. Copy the density gradient
  ! 2. Compute the velocity gradients (nmx,nmy,nmz)
  ! 3. Compute the total energy gradient and store in the
  !    location for the temperature gradient (nee)
  !
  !$ACC PARALLEL LOOP PRESENT(usp,dcvdx,return_value)
  do n = n1,n2
    !
    rho_inv = one/usp(n,nec)
    !
    decdx = dcvdx(n,1,nec)
    decdy = dcvdx(n,2,nec)
    decdz = dcvdx(n,3,nec)
    !
    return_value(n,1,nec) = decdx
    return_value(n,2,nec) = decdy
    return_value(n,3,nec) = decdz
    !
    decdx = decdx*rho_inv*rho_inv
    decdy = decdy*rho_inv*rho_inv
    decdz = decdz*rho_inv*rho_inv
    !
    return_value(n,1,nmx) = rho_inv*dcvdx(n,1,nmx) - decdx*usp(n,nmx)
    return_value(n,2,nmx) = rho_inv*dcvdx(n,2,nmx) - decdy*usp(n,nmx)
    return_value(n,3,nmx) = rho_inv*dcvdx(n,3,nmx) - decdz*usp(n,nmx)
    !
    return_value(n,1,nmy) = rho_inv*dcvdx(n,1,nmy) - decdx*usp(n,nmy)
    return_value(n,2,nmy) = rho_inv*dcvdx(n,2,nmy) - decdy*usp(n,nmy)
    return_value(n,3,nmy) = rho_inv*dcvdx(n,3,nmy) - decdz*usp(n,nmy)
    !
    return_value(n,1,nmz) = rho_inv*dcvdx(n,1,nmz) - decdx*usp(n,nmz)
    return_value(n,2,nmz) = rho_inv*dcvdx(n,2,nmz) - decdy*usp(n,nmz)
    return_value(n,3,nmz) = rho_inv*dcvdx(n,3,nmz) - decdz*usp(n,nmz)
    !
    return_value(n,1,nee) = rho_inv*dcvdx(n,1,nee) - decdx*usp(n,nee)
    return_value(n,2,nee) = rho_inv*dcvdx(n,2,nee) - decdy*usp(n,nee)
    return_value(n,3,nee) = rho_inv*dcvdx(n,3,nee) - decdz*usp(n,nee)
    !
   !return_value(n,l,m) = rho_inv*dcvdx(n,l,m) - &
   !                      rho_inv*dcvdx(n,l,nec)*usp(n,m)
    !
  end do
  !
  ! The location for the temperature gradient currently contains
  ! the total energy gradient from the previous loop.
  ! Now use the velocity gradients to subtract the gradient of
  ! kinetic-energy-per-unit-mass from the total energy gradient.
  !
  !$ACC PARALLEL LOOP PRESENT(usp,return_value)
  do n = n1,n2
    !
    rho_inv = one/usp(n,nec)
    !
    return_value(n,1,nee) = return_value(n,1,nee) &
                          - return_value(n,1,nmx)*usp(n,nmx)*rho_inv &
                          - return_value(n,1,nmy)*usp(n,nmy)*rho_inv &
                          - return_value(n,1,nmz)*usp(n,nmz)*rho_inv
    !
    return_value(n,2,nee) = return_value(n,2,nee) &
                          - return_value(n,2,nmx)*usp(n,nmx)*rho_inv &
                          - return_value(n,2,nmy)*usp(n,nmy)*rho_inv &
                          - return_value(n,2,nmz)*usp(n,nmz)*rho_inv
    !
    return_value(n,3,nee) = return_value(n,3,nee) &
                          - return_value(n,3,nmx)*usp(n,nmx)*rho_inv &
                          - return_value(n,3,nmy)*usp(n,nmy)*rho_inv &
                          - return_value(n,3,nmz)*usp(n,nmz)*rho_inv
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
end subroutine grad_cv_to_grad_pv_cpg_3d_original
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine grad_cv_to_grad_pv_cpg_2d(n1,n2,usp,dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: dcvdx(:,:,:)
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
  !$ACC PARALLEL LOOP PRESENT(dcvdx,return_value)
  do n = n1,n2
    return_value(n,1,nec) = dcvdx(n,1,nec)
    return_value(n,2,nec) = dcvdx(n,2,nec)
  end do
  !
  ! Using the product rule:
  !
  ! 1. Compute the velocity gradients (nmx,nmy)
  ! 2. Compute the total energy gradient and store in the
  !    location for the temperature gradient (nee)
  !
  !$ACC PARALLEL LOOP PRESENT(usp,dcvdx,return_value)
  do n = n1,n2
    !
    rho_inv = one/usp(n,nec)
    !
    return_value(n,1,nmx) = rho_inv*dcvdx(n,1,nmx) &
                          - rho_inv*dcvdx(n,1,nec)*usp(n,nmx)*rho_inv
    return_value(n,2,nmx) = rho_inv*dcvdx(n,2,nmx) &
                          - rho_inv*dcvdx(n,2,nec)*usp(n,nmx)*rho_inv
    !
    return_value(n,1,nmy) = rho_inv*dcvdx(n,1,nmy) &
                          - rho_inv*dcvdx(n,1,nec)*usp(n,nmy)*rho_inv
    return_value(n,2,nmy) = rho_inv*dcvdx(n,2,nmy) &
                          - rho_inv*dcvdx(n,2,nec)*usp(n,nmy)*rho_inv
    !
    return_value(n,1,nee) = rho_inv*dcvdx(n,1,nee) &
                          - rho_inv*dcvdx(n,1,nec)*usp(n,nee)*rho_inv
    return_value(n,2,nee) = rho_inv*dcvdx(n,2,nee) &
                          - rho_inv*dcvdx(n,2,nec)*usp(n,nee)*rho_inv
    !
   !return_value(n,l,m) = rho_inv*dcvdx(n,l,m) - &
   !                      rho_inv*dcvdx(n,l,nec)*usp(n,m)
    !
  end do
  !
  ! The location for the temperature gradient currently contains
  ! the total energy gradient from the previous loop.
  ! Now use the velocity gradients to subtract the gradient of
  ! kinetic-energy-per-unit-mass from the total energy gradient.
  !
  !$ACC PARALLEL LOOP PRESENT(usp,return_value)
  do n = n1,n2
    !
    rho_inv = one/usp(n,nec)
    !
    return_value(n,1,nee) = return_value(n,1,nee) &
                          - return_value(n,1,nmx)*usp(n,nmx)*rho_inv &
                          - return_value(n,1,nmy)*usp(n,nmy)*rho_inv
    !
    return_value(n,2,nee) = return_value(n,2,nee) &
                          - return_value(n,2,nmx)*usp(n,nmx)*rho_inv &
                          - return_value(n,2,nmy)*usp(n,nmy)*rho_inv
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
end subroutine grad_cv_to_grad_pv_cpg_2d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine grad_cv_to_grad_pv_cpg_1d(n1,n2,usp,dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: dcvdx(:,:,:)
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
  !$ACC PARALLEL LOOP PRESENT(dcvdx,return_value)
  do n = n1,n2
    return_value(n,1,nec) = dcvdx(n,1,nec)
  end do
  !
  ! Using the product rule:
  !
  ! 1. Compute the velocity gradients (m = 2,nq-1)
  ! 2. Compute the total energy gradient and store in the
  !    location for the temperature gradient (m = nq)
  !
  !$ACC PARALLEL LOOP PRESENT(usp,dcvdx,return_value)
  do n = n1,n2
    !
    rho_inv = one/usp(n,nec)
    !
    return_value(n,1,nmx) = rho_inv*dcvdx(n,1,nmx) &
                          - rho_inv*dcvdx(n,1,nec)*usp(n,nmx)*rho_inv
    !
    return_value(n,1,nee) = rho_inv*dcvdx(n,1,nee) &
                          - rho_inv*dcvdx(n,1,nec)*usp(n,nee)*rho_inv
    !
  end do
  !
  ! The location for the temperature gradient currently contains
  ! the total energy gradient from the previous loop.
  ! Now use the velocity gradients to subtract the gradient of
  ! kinetic-energy-per-unit-mass from the total energy gradient.
  !
  !
  !$ACC PARALLEL LOOP PRESENT(usp,return_value)
  do n = n1,n2
    return_value(n,1,nee) = return_value(n,1,nee) &
                          - return_value(n,1,nmx)*usp(n,nmx)/usp(n,nec)
  end do
  !
  ! Finish computing the temperature gradient by dividing through by
  ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
  !
  !$ACC PARALLEL LOOP PRESENT(return_value)
  do n = n1,n2
    return_value(n,1,nee) = return_value(n,1,nee) * cv_inv
  end do
  !
end subroutine grad_cv_to_grad_pv_cpg_1d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine adjust_temp_gradient_for_tpg_3d(n1,n2,usp,gpsp,dcvdx, &
                                           return_value)
  !
  !.. Use Statements ..
  use ovar,       only : rgasref_nd,gam
  use thermo_mod, only : rspec
  use thermo_mod, only : compute_species_cp
  use thermo_mod, only : compute_mixture_dcpdt
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: gpsp(:,:)
  real(wp),    intent(in) :: dcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,np
  real(wp) :: rho_inv,mass_frac
  real(wp) :: decdx,decdy,decdz,dtdx,dtdy,dtdz
  real(wp) :: val,valx,valy,valz
  real(wp) :: temp,gamm,cp
  real(wp) :: constant_cv,coeff
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: dcpdt(:)
  real(wp), allocatable :: cpi(:,:)
  !
continue
  !
  np = size(usp,dim=1)
  !
  constant_cv = rgasref_nd / (gam - one)
  !
  allocate ( cpi(1:np,1:ncs) )
  allocate ( dcpdt(1:np) )
  !$ACC DATA CREATE(cpi,dcpdt)
  !
  ! The location for the temperature gradient currently contains the temperature
  ! gradient based on a constant 'specific heat at constant volume', i.e.,
  ! constant gamma and Cp. For a gas mixture, gamma and Cp vary w-r-t mass
  ! fractions and temperature so we need to correct the temperature gradient
  ! accordingly.
  !
  ! Using the product rule:
  !
  ! 1. Compute the species mass fraction gradients. Also compute the species
  !    mass fractions which will be needed to compute d(Cp)/dT.
  !    NOTE: Temporarily store the species mass fractions in cpi
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,cpi,dcvdx,return_value)
  do m = 1,ncs
    do n = n1,n2
      !
      rho_inv = one/usp(n,nec)
      !
      mass_frac = usp(n,m)*rho_inv
      !
      cpi(n,m) = mass_frac
      !
      decdx = dcvdx(n,1,nec)*rho_inv
      decdy = dcvdx(n,2,nec)*rho_inv
      decdz = dcvdx(n,3,nec)*rho_inv
      !
      return_value(n,1,m) = rho_inv*dcvdx(n,1,m) - decdx*mass_frac
      return_value(n,2,m) = rho_inv*dcvdx(n,2,m) - decdy*mass_frac
      return_value(n,3,m) = rho_inv*dcvdx(n,3,m) - decdz*mass_frac
      !
    end do
  end do
  !
  ! 2. Compute the derivative of the mixture Cp w-r-t temperature: d(Cp)/dT.
  !    NOTE: The species mass fractions are temporarily being stored in cpi.
  !
  call compute_mixture_dcpdt(n1,n2,cpi,gpsp(:,ntemp),dcpdt)
  !
  ! 3. Compute the species Cp.
  ! NOTE: cpi currently contains the species mass fractions. These are not
  !       needed in compute_species_cp, so the species mass fractions in
  !       cpi will simply be overwritten by the species Cp.
  !
  call compute_species_cp(n1,n2,gpsp(:,ntemp),cpi)
  !
  ! 4. Compute the temperature gradients of a thermally perfect gas mixture
  !
  !$ACC PARALLEL LOOP PRESENT(rspec,cpi,dcpdt,gpsp,return_value)
  do n = n1,n2
    !
    ! Compute the leading coefficient for the temperature gradient. For a
    ! calorically perfect gas, this coefficient is 1/Cv; for a thermally
    ! perfect gas mixture, this coefficient is:
    !          1                      1                       gamma
    !  ----------------- = ----------------------- = -----------------------
    !  Cv + T*[d(Cp)/dT]   Cp/gamma + T*[d(Cp)/dT]   Cp + gamma*T*[d(Cp)/dT]
    !
    temp = gpsp(n,ntemp)
    gamm = gpsp(n,ngamm)
    cp   = gpsp(n,ncp)
    !
    coeff = gamm / (cp + gamm*temp*dcpdt(n))
    !
    ! Compute the contribution to the temperature gradient
    ! from the species mass fractions gradients.
    !
    valx = zero
    valy = zero
    valz = zero
    !
    !$ACC LOOP SEQ REDUCTION(+:valx,valy,valz)
    do m = 1,ncs
      val = rspec(m) - cpi(n,m)
      valx = valx + val*return_value(n,1,m)
      valy = valy + val*return_value(n,2,m)
      valz = valz + val*return_value(n,3,m)
    end do
    !
    ! Multiply the previously computed temperature gradients by the constant
    ! Cv used in the grad_cv_to_grad_pv_cpg_[321]d subroutines.
    !
    dtdx = return_value(n,1,nee)*constant_cv
    dtdy = return_value(n,2,nee)*constant_cv
    dtdz = return_value(n,3,nee)*constant_cv
    !
    ! Add the contributions of the species mass fraction gradients to the
    ! temperature gradient and then scale it by the proper leading coefficient
    ! computed above.
    !
    return_value(n,1,nee) = coeff*(dtdx + temp*valx)
    return_value(n,2,nee) = coeff*(dtdy + temp*valy)
    return_value(n,3,nee) = coeff*(dtdz + temp*valz)
    !
  end do
  !
  !$ACC END DATA
  deallocate ( dcpdt )
  deallocate ( cpi )
  !
end subroutine adjust_temp_gradient_for_tpg_3d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine adjust_temp_gradient_for_tpg_2d(n1,n2,usp,gpsp,dcvdx, &
                                           return_value)
  !
  !.. Use Statements ..
  use ovar,       only : rgasref_nd,gam
  use thermo_mod, only : rspec
  use thermo_mod, only : compute_species_cp
  use thermo_mod, only : compute_mixture_dcpdt
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: gpsp(:,:)
  real(wp),    intent(in) :: dcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,np
  real(wp) :: rho_inv,mass_frac
  real(wp) :: decdx,decdy,dtdx,dtdy
  real(wp) :: val,valx,valy
  real(wp) :: temp,gamm,cp
  real(wp) :: constant_cv,coeff
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: dcpdt(:)
  real(wp), allocatable :: cpi(:,:)
  !
continue
  !
  np = size(usp,dim=1)
  !
  constant_cv = rgasref_nd / (gam - one)
  !
  allocate ( cpi(1:np,1:ncs) )
  allocate ( dcpdt(1:np) )
  !$ACC DATA CREATE(cpi,dcpdt)
  !
  ! The location for the temperature gradient currently contains the temperature
  ! gradient based on a constant 'specific heat at constant volume', i.e.,
  ! constant gamma and Cp. For a gas mixture, gamma and Cp vary w-r-t mass
  ! fractions and temperature so we need to correct the temperature gradient
  ! accordingly.
  !
  ! Using the product rule:
  !
  ! 1. Compute the species mass fraction gradients. Also compute the species
  !    mass fractions which will be needed to compute d(Cp)/dT.
  !    NOTE: Temporarily store the species mass fractions in cpi
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,cpi,dcvdx,return_value)
  do m = 1,ncs
    do n = n1,n2
      !
      rho_inv = one/usp(n,nec)
      !
      mass_frac = usp(n,m)*rho_inv
      !
      cpi(n,m) = mass_frac
      !
      decdx = dcvdx(n,1,nec)*rho_inv
      decdy = dcvdx(n,2,nec)*rho_inv
      !
      return_value(n,1,m) = rho_inv*dcvdx(n,1,m) - decdx*mass_frac
      return_value(n,2,m) = rho_inv*dcvdx(n,2,m) - decdy*mass_frac
      !
    end do
  end do
  !
  ! 2. Compute the derivative of the mixture Cp w-r-t temperature: d(Cp)/dT.
  !    NOTE: The species mass fractions are temporarily being stored in cpi.
  !
  call compute_mixture_dcpdt(n1,n2,cpi,gpsp(:,ntemp),dcpdt)
  !
  ! 3. Compute the species Cp.
  ! NOTE: cpi currently contains the species mass fractions. These are not
  !       needed in compute_species_cp, so the species mass fractions in
  !       cpi will simply be overwritten by the species Cp.
  !
  call compute_species_cp(n1,n2,gpsp(:,ntemp),cpi)
  !
  ! 4. Compute the temperature gradients of a thermally perfect gas mixture
  !
  !$ACC PARALLEL LOOP PRESENT(rspec,cpi,dcpdt,gpsp,return_value)
  do n = n1,n2
    !
    ! Compute the leading coefficient for the temperature gradient. For a
    ! calorically perfect gas, this coefficient is 1/Cv; for a thermally
    ! perfect gas mixture, this coefficient is:
    !          1                      1                       gamma
    !  ----------------- = ----------------------- = -----------------------
    !  Cv + T*[d(Cp)/dT]   Cp/gamma + T*[d(Cp)/dT]   Cp + gamma*T*[d(Cp)/dT]
    !
    temp = gpsp(n,ntemp)
    gamm = gpsp(n,ngamm)
    cp   = gpsp(n,ncp)
    !
    coeff = gamm / (cp + gamm*temp*dcpdt(n))
    !
    ! Compute the contribution to the temperature gradient
    ! from the species mass fractions gradients.
    !
    valx = zero
    valy = zero
    !
    !$ACC LOOP SEQ REDUCTION(+:valx,valy)
    do m = 1,ncs
      val = rspec(m) - cpi(n,m)
      valx = valx + val*return_value(n,1,m)
      valy = valy + val*return_value(n,2,m)
    end do
    !
    ! Multiply the previously computed temperature gradients by the constant
    ! Cv used in the grad_cv_to_grad_pv_cpg_[321]d subroutines.
    !
    dtdx = return_value(n,1,nee)*constant_cv
    dtdy = return_value(n,2,nee)*constant_cv
    !
    ! Add the contributions of the species mass fraction gradients to the
    ! temperature gradient and then scale it by the proper leading coefficient
    ! computed above.
    !
    return_value(n,1,nee) = coeff*(dtdx + temp*valx)
    return_value(n,2,nee) = coeff*(dtdy + temp*valy)
    !
  end do
  !
  !$ACC END DATA
  deallocate ( dcpdt )
  deallocate ( cpi )
  !
end subroutine adjust_temp_gradient_for_tpg_2d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine adjust_temp_gradient_for_tpg_1d(n1,n2,usp,gpsp,dcvdx, &
                                           return_value)
  !
  !.. Use Statements ..
  use ovar,       only : rgasref_nd,gam
  use thermo_mod, only : rspec
  use thermo_mod, only : compute_species_cp
  use thermo_mod, only : compute_mixture_dcpdt
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: gpsp(:,:)
  real(wp),    intent(in) :: dcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,np
  real(wp) :: rho_inv,mass_frac
  real(wp) :: decdx,dtdx
  real(wp) :: val,valx
  real(wp) :: temp,gamm,cp
  real(wp) :: constant_cv,coeff
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: dcpdt(:)
  real(wp), allocatable :: cpi(:,:)
  !
continue
  !
  np = size(usp,dim=1)
  !
  constant_cv = rgasref_nd / (gam - one)
  !
  allocate ( cpi(1:np,1:ncs) )
  allocate ( dcpdt(1:np) )
  !$ACC DATA CREATE(cpi,dcpdt)
  !
  ! The location for the temperature gradient currently contains the temperature
  ! gradient based on a constant 'specific heat at constant volume', i.e.,
  ! constant gamma and Cp. For a gas mixture, gamma and Cp vary w-r-t mass
  ! fractions and temperature so we need to correct the temperature gradient
  ! accordingly.
  !
  ! Using the product rule:
  !
  ! 1. Compute the species mass fraction gradients. Also compute the species
  !    mass fractions which will be needed to compute d(Cp)/dT.
  !    NOTE: Temporarily store the species mass fractions in cpi
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,cpi,dcvdx,return_value)
  do m = 1,ncs
    do n = n1,n2
      !
      rho_inv = one/usp(n,nec)
      !
      mass_frac = usp(n,m)*rho_inv
      !
      cpi(n,m) = mass_frac
      !
      decdx = dcvdx(n,1,nec)*rho_inv
      !
      return_value(n,1,m) = rho_inv*dcvdx(n,1,m) - decdx*mass_frac
      !
    end do
  end do
  !
  ! 2. Compute the derivative of the mixture Cp w-r-t temperature: d(Cp)/dT.
  !    NOTE: The species mass fractions are temporarily being stored in cpi.
  !
  call compute_mixture_dcpdt(n1,n2,cpi,gpsp(:,ntemp),dcpdt)
  !
  ! 3. Compute the species Cp.
  ! NOTE: cpi currently contains the species mass fractions. These are not
  !       needed in compute_species_cp, so the species mass fractions in
  !       cpi will simply be overwritten by the species Cp.
  !
  call compute_species_cp(n1,n2,gpsp(:,ntemp),cpi)
  !
  ! 4. Compute the temperature gradients of a thermally perfect gas mixture
  !
  !$ACC PARALLEL LOOP PRESENT(rspec,cpi,dcpdt,gpsp,return_value)
  do n = n1,n2
    !
    ! Compute the leading coefficient for the temperature gradient. For a
    ! calorically perfect gas, this coefficient is 1/Cv; for a thermally
    ! perfect gas mixture, this coefficient is:
    !          1                      1                       gamma
    !  ----------------- = ----------------------- = -----------------------
    !  Cv + T*[d(Cp)/dT]   Cp/gamma + T*[d(Cp)/dT]   Cp + gamma*T*[d(Cp)/dT]
    !
    temp = gpsp(n,ntemp)
    gamm = gpsp(n,ngamm)
    cp   = gpsp(n,ncp)
    !
    coeff = gamm / (cp + gamm*temp*dcpdt(n))
    !
    ! Compute the contribution to the temperature gradient
    ! from the species mass fractions gradients.
    !
    valx = zero
    !
    !$ACC LOOP SEQ REDUCTION(+:valx)
    do m = 1,ncs
      val = rspec(m) - cpi(n,m)
      valx = valx + val*return_value(n,1,m)
    end do
    !
    ! Multiply the previously computed temperature gradients by the constant
    ! Cv used in the grad_cv_to_grad_pv_cpg_[321]d subroutines.
    !
    dtdx = return_value(n,1,nee)*constant_cv
    !
    ! Add the contributions of the species mass fraction gradients to the
    ! temperature gradient and then scale it by the proper leading coefficient
    ! computed above.
    !
    return_value(n,1,nee) = coeff*(dtdx + temp*valx)
    !
  end do
  !
  !$ACC END DATA
  deallocate ( dcpdt )
  deallocate ( cpi )
  !
end subroutine adjust_temp_gradient_for_tpg_1d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine adjust_temp_grad_for_tpg_3d(n1,n2,usp,gpsp,dcvdx, &
                                       return_value,output_variables)
  !
  !.. Use Statements ..
  use ovar,       only : rgasref_nd,gam
  use thermo_mod, only : rspec
  use thermo_mod, only : compute_species_cp
  use thermo_mod, only : compute_species_dcpdt
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: gpsp(:,:)
  real(wp),    intent(in) :: dcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  real(wp), intent(inout) :: output_variables(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,np
  real(wp) :: rho_inv,ymf
  real(wp) :: decdx,decdy,decdz,dtdx,dtdy,dtdz
  real(wp) :: val,valx,valy,valz
  real(wp) :: temp,gamm,cp,dcpdt
  real(wp) :: constant_cv,coeff,cv_inv
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: dcpidt(:,:)
  real(wp), allocatable :: cpi(:,:)
  real(wp), allocatable :: mass_frac(:,:)
  !
continue
  !
  np = size(usp,dim=1)
  !
  cv_inv = (gam - one) / rgasref_nd
  constant_cv = rgasref_nd / (gam - one)
  !
  allocate ( mass_frac(1:np,1:ncs) )
  allocate ( cpi(1:np,1:ncs) )
  allocate ( dcpidt(1:np,1:ncs) )
  !$ACC DATA CREATE(cpi,dcpidt)
  !
  ! The location for the temperature gradient currently contains the temperature
  ! gradient based on a constant 'specific heat at constant volume', i.e.,
  ! constant gamma and Cp. For a gas mixture, gamma and Cp vary w-r-t mass
  ! fractions and temperature so we need to correct the temperature gradient
  ! accordingly.
  !
  ! Using the product rule:
  !
  ! 1. Compute the species mass fraction gradients. Also compute the species
  !    mass fractions which will be needed to compute d(Cp)/dT.
  !    NOTE: Temporarily store the species mass fractions in cpi
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,cpi,dcvdx,return_value)
  do m = 1,ncs
    do n = n1,n2
      !
      rho_inv = one/usp(n,nec)
      !
      ymf = usp(n,m)*rho_inv
      !
      mass_frac(n,m) = ymf
      !
      decdx = dcvdx(n,1,nec)*rho_inv
      decdy = dcvdx(n,2,nec)*rho_inv
      decdz = dcvdx(n,3,nec)*rho_inv
      !
      return_value(n,1,m) = rho_inv*dcvdx(n,1,m) - decdx*ymf
      return_value(n,2,m) = rho_inv*dcvdx(n,2,m) - decdy*ymf
      return_value(n,3,m) = rho_inv*dcvdx(n,3,m) - decdz*ymf
      !
    end do
  end do
  !
  ! 2. Compute the derivative of the mixture Cp w-r-t temperature: d(Cp)/dT.
  !    NOTE: The species mass fractions are temporarily being stored in cpi.
  !
  call compute_species_dcpdt(n1,n2,gpsp(:,ntemp),dcpidt)
  !
  ! 3. Compute the species Cp.
  ! NOTE: cpi currently contains the species mass fractions. These are not
  !       needed in compute_species_cp, so the species mass fractions in
  !       cpi will simply be overwritten by the species Cp.
  !
  call compute_species_cp(n1,n2,gpsp(:,ntemp),cpi)
  !
  ! 4. Compute the temperature gradients of a thermally perfect gas mixture
  !
  !$ACC PARALLEL LOOP PRESENT(rspec,cpi,dcpidt,gpsp,return_value)
  do n = n1,n2
    !
    ! Multiply the previously computed temperature gradients by the constant
    ! Cv used in the grad_cv_to_grad_pv_cpg_[321]d subroutines.
    !
    output_variables(n,1) = return_value(n,1,nee)
    output_variables(n,2) = return_value(n,2,nee)
    output_variables(n,3) = return_value(n,3,nee)
    !
    dtdx = return_value(n,1,nee)*constant_cv
    dtdy = return_value(n,2,nee)*constant_cv
    dtdz = return_value(n,3,nee)*constant_cv
    !
    ! Compute the leading coefficient for the temperature gradient. For a
    ! calorically perfect gas, this coefficient is 1/Cv; for a thermally
    ! perfect gas mixture, this coefficient is:
    !          1                      1                       gamma
    !  ----------------- = ----------------------- = -----------------------
    !  Cv + T*[d(Cp)/dT]   Cp/gamma + T*[d(Cp)/dT]   Cp + gamma*T*[d(Cp)/dT]
    !
    temp = gpsp(n,ntemp)
    gamm = gpsp(n,ngamm)
    cp   = gpsp(n,ncp)
    !
    dcpdt = zero
    do m = 1,ncs
      dcpdt = dcpdt + mass_frac(n,m)*dcpidt(n,m)
    end do
    !
    coeff = gamm / (cp + gamm*temp*dcpdt)
    !
    output_variables(n,3+2*ncs+1) = dcpdt
    output_variables(n,3+2*ncs+2) = temp
    output_variables(n,3+2*ncs+3) = cp
    output_variables(n,3+2*ncs+4) = gamm
    output_variables(n,3+2*ncs+5) = coeff
    output_variables(n,3+2*ncs+6) = constant_cv
    output_variables(n,3+2*ncs+7) = cv_inv
    !
    ! Compute the contribution to the temperature gradient
    ! from the species mass fractions gradients.
    !
    valx = zero
    valy = zero
    valz = zero
    !
    !$ACC LOOP SEQ REDUCTION(+:valx,valy,valz)
    do m = 1,ncs
      !
      val = rspec(m) - cpi(n,m)
      valx = valx + val*return_value(n,1,m)
      valy = valy + val*return_value(n,2,m)
      valz = valz + val*return_value(n,3,m)
      !
      output_variables(n,3+m) = -val
      output_variables(n,3+ncs+m) = dcpidt(n,m)
      !
    end do
    !
    ! Add the contributions of the species mass fraction gradients to the
    ! temperature gradient and then scale it by the proper leading coefficient
    ! computed above.
    !
    return_value(n,1,nee) = coeff*(dtdx + temp*valx)
    return_value(n,2,nee) = coeff*(dtdy + temp*valy)
    return_value(n,3,nee) = coeff*(dtdz + temp*valz)
    !
  end do
  !
  !$ACC END DATA
  deallocate ( dcpidt )
  deallocate ( cpi )
  deallocate ( mass_frac )
  !
end subroutine adjust_temp_grad_for_tpg_3d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine adjust_temp_grad_for_tpg_2d(n1,n2,usp,gpsp,dcvdx, &
                                       return_value,output_variables)
  !
  !.. Use Statements ..
  use ovar,       only : rgasref_nd,gam
  use thermo_mod, only : rspec
  use thermo_mod, only : compute_species_cp
  use thermo_mod, only : compute_species_dcpdt
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: gpsp(:,:)
  real(wp),    intent(in) :: dcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  real(wp), intent(inout) :: output_variables(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,np
  real(wp) :: rho_inv,ymf
  real(wp) :: decdx,decdy,dtdx,dtdy
  real(wp) :: val,valx,valy
  real(wp) :: temp,gamm,cp,dcpdt
  real(wp) :: constant_cv,coeff,cv_inv
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: dcpidt(:,:)
  real(wp), allocatable :: cpi(:,:)
  real(wp), allocatable :: mass_frac(:,:)
  !
continue
  !
  np = size(usp,dim=1)
  !
  cv_inv = (gam - one) / rgasref_nd
  constant_cv = rgasref_nd / (gam - one)
  !
  allocate ( mass_frac(1:np,1:ncs) )
  allocate ( cpi(1:np,1:ncs) )
  allocate ( dcpidt(1:np,1:ncs) )
  !$ACC DATA CREATE(cpi,dcpidt)
  !
  ! The location for the temperature gradient currently contains the temperature
  ! gradient based on a constant 'specific heat at constant volume', i.e.,
  ! constant gamma and Cp. For a gas mixture, gamma and Cp vary w-r-t mass
  ! fractions and temperature so we need to correct the temperature gradient
  ! accordingly.
  !
  ! Using the product rule:
  !
  ! 1. Compute the species mass fraction gradients. Also compute the species
  !    mass fractions which will be needed to compute d(Cp)/dT.
  !    NOTE: Temporarily store the species mass fractions in cpi
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,cpi,dcvdx,return_value)
  do m = 1,ncs
    do n = n1,n2
      !
      rho_inv = one/usp(n,nec)
      !
      ymf = usp(n,m)*rho_inv
      !
      mass_frac(n,m) = ymf
      !
      decdx = dcvdx(n,1,nec)*rho_inv
      decdy = dcvdx(n,2,nec)*rho_inv
      !
      return_value(n,1,m) = rho_inv*dcvdx(n,1,m) - decdx*ymf
      return_value(n,2,m) = rho_inv*dcvdx(n,2,m) - decdy*ymf
      !
    end do
  end do
  !
  ! 2. Compute the derivative of the mixture Cp w-r-t temperature: d(Cp)/dT.
  !    NOTE: The species mass fractions are temporarily being stored in cpi.
  !
  call compute_species_dcpdt(n1,n2,gpsp(:,ntemp),dcpidt)
  !
  ! 3. Compute the species Cp.
  ! NOTE: cpi currently contains the species mass fractions. These are not
  !       needed in compute_species_cp, so the species mass fractions in
  !       cpi will simply be overwritten by the species Cp.
  !
  call compute_species_cp(n1,n2,gpsp(:,ntemp),cpi)
  !
  ! 4. Compute the temperature gradients of a thermally perfect gas mixture
  !
  !$ACC PARALLEL LOOP PRESENT(rspec,cpi,dcpidt,gpsp,return_value)
  do n = n1,n2
    !
    ! Multiply the previously computed temperature gradients by the constant
    ! Cv used in the grad_cv_to_grad_pv_cpg_[321]d subroutines.
    !
    output_variables(n,1) = return_value(n,1,nee)
    output_variables(n,2) = return_value(n,2,nee)
    !
    dtdx = return_value(n,1,nee)*constant_cv
    dtdy = return_value(n,2,nee)*constant_cv
    !
    ! Compute the leading coefficient for the temperature gradient. For a
    ! calorically perfect gas, this coefficient is 1/Cv; for a thermally
    ! perfect gas mixture, this coefficient is:
    !          1                      1                       gamma
    !  ----------------- = ----------------------- = -----------------------
    !  Cv + T*[d(Cp)/dT]   Cp/gamma + T*[d(Cp)/dT]   Cp + gamma*T*[d(Cp)/dT]
    !
    temp = gpsp(n,ntemp)
    gamm = gpsp(n,ngamm)
    cp   = gpsp(n,ncp)
    !
    dcpdt = zero
    do m = 1,ncs
      dcpdt = dcpdt + mass_frac(n,m)*dcpidt(n,m)
    end do
    !
    coeff = gamm / (cp + gamm*temp*dcpdt)
    !
    output_variables(n,2+2*ncs+1) = dcpdt
    output_variables(n,2+2*ncs+2) = temp
    output_variables(n,2+2*ncs+3) = cp
    output_variables(n,2+2*ncs+4) = gamm
    output_variables(n,2+2*ncs+5) = coeff
    output_variables(n,2+2*ncs+6) = constant_cv
    output_variables(n,2+2*ncs+7) = cv_inv
    !
    ! Compute the contribution to the temperature gradient
    ! from the species mass fractions gradients.
    !
    valx = zero
    valy = zero
    !
    !$ACC LOOP SEQ REDUCTION(+:valx,valy,valz)
    do m = 1,ncs
      !
      val = rspec(m) - cpi(n,m)
      valx = valx + val*return_value(n,1,m)
      valy = valy + val*return_value(n,2,m)
      !
      output_variables(n,2+m) = -val
      output_variables(n,2+ncs+m) = dcpidt(n,m)
      !
    end do
    !
    ! Add the contributions of the species mass fraction gradients to the
    ! temperature gradient and then scale it by the proper leading coefficient
    ! computed above.
    !
    return_value(n,1,nee) = coeff*(dtdx + temp*valx)
    return_value(n,2,nee) = coeff*(dtdy + temp*valy)
    !
  end do
  !
  !$ACC END DATA
  deallocate ( dcpidt )
  deallocate ( cpi )
  deallocate ( mass_frac )
  !
end subroutine adjust_temp_grad_for_tpg_2d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine adjust_temp_grad_for_tpg_1d(n1,n2,usp,gpsp,dcvdx, &
                                       return_value,output_variables)
  !
  !.. Use Statements ..
  use ovar,       only : rgasref_nd,gam
  use thermo_mod, only : rspec
  use thermo_mod, only : compute_species_cp
  use thermo_mod, only : compute_species_dcpdt
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: gpsp(:,:)
  real(wp),    intent(in) :: dcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  real(wp), intent(inout) :: output_variables(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,np
  real(wp) :: rho_inv,ymf
  real(wp) :: decdx,dtdx
  real(wp) :: val,valx
  real(wp) :: temp,gamm,cp,dcpdt
  real(wp) :: constant_cv,coeff,cv_inv
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: dcpidt(:,:)
  real(wp), allocatable :: cpi(:,:)
  real(wp), allocatable :: mass_frac(:,:)
  !
continue
  !
  np = size(usp,dim=1)
  !
  cv_inv = (gam - one) / rgasref_nd
  constant_cv = rgasref_nd / (gam - one)
  !
  allocate ( mass_frac(1:np,1:ncs) )
  allocate ( cpi(1:np,1:ncs) )
  allocate ( dcpidt(1:np,1:ncs) )
  !$ACC DATA CREATE(cpi,dcpidt)
  !
  ! The location for the temperature gradient currently contains the temperature
  ! gradient based on a constant 'specific heat at constant volume', i.e.,
  ! constant gamma and Cp. For a gas mixture, gamma and Cp vary w-r-t mass
  ! fractions and temperature so we need to correct the temperature gradient
  ! accordingly.
  !
  ! Using the product rule:
  !
  ! 1. Compute the species mass fraction gradients. Also compute the species
  !    mass fractions which will be needed to compute d(Cp)/dT.
  !    NOTE: Temporarily store the species mass fractions in cpi
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,cpi,dcvdx,return_value)
  do m = 1,ncs
    do n = n1,n2
      !
      rho_inv = one/usp(n,nec)
      !
      ymf = usp(n,m)*rho_inv
      !
      mass_frac(n,m) = ymf
      !
      decdx = dcvdx(n,1,nec)*rho_inv
      !
      return_value(n,1,m) = rho_inv*dcvdx(n,1,m) - decdx*ymf
      !
    end do
  end do
  !
  ! 2. Compute the derivative of the mixture Cp w-r-t temperature: d(Cp)/dT.
  !    NOTE: The species mass fractions are temporarily being stored in cpi.
  !
  call compute_species_dcpdt(n1,n2,gpsp(:,ntemp),dcpidt)
  !
  ! 3. Compute the species Cp.
  ! NOTE: cpi currently contains the species mass fractions. These are not
  !       needed in compute_species_cp, so the species mass fractions in
  !       cpi will simply be overwritten by the species Cp.
  !
  call compute_species_cp(n1,n2,gpsp(:,ntemp),cpi)
  !
  ! 4. Compute the temperature gradients of a thermally perfect gas mixture
  !
  !$ACC PARALLEL LOOP PRESENT(rspec,cpi,dcpidt,gpsp,return_value)
  do n = n1,n2
    !
    ! Multiply the previously computed temperature gradients by the constant
    ! Cv used in the grad_cv_to_grad_pv_cpg_[321]d subroutines.
    !
    output_variables(n,1) = return_value(n,1,nee)
    !
    dtdx = return_value(n,1,nee)*constant_cv
    !
    ! Compute the leading coefficient for the temperature gradient. For a
    ! calorically perfect gas, this coefficient is 1/Cv; for a thermally
    ! perfect gas mixture, this coefficient is:
    !          1                      1                       gamma
    !  ----------------- = ----------------------- = -----------------------
    !  Cv + T*[d(Cp)/dT]   Cp/gamma + T*[d(Cp)/dT]   Cp + gamma*T*[d(Cp)/dT]
    !
    temp = gpsp(n,ntemp)
    gamm = gpsp(n,ngamm)
    cp   = gpsp(n,ncp)
    !
    dcpdt = zero
    do m = 1,ncs
      dcpdt = dcpdt + mass_frac(n,m)*dcpidt(n,m)
    end do
    !
    coeff = gamm / (cp + gamm*temp*dcpdt)
    !
    output_variables(n,1+2*ncs+1) = dcpdt
    output_variables(n,1+2*ncs+2) = temp
    output_variables(n,1+2*ncs+3) = cp
    output_variables(n,1+2*ncs+4) = gamm
    output_variables(n,1+2*ncs+5) = coeff
    output_variables(n,1+2*ncs+6) = constant_cv
    output_variables(n,1+2*ncs+7) = cv_inv
    !
    ! Compute the contribution to the temperature gradient
    ! from the species mass fractions gradients.
    !
    valx = zero
    !
    !$ACC LOOP SEQ REDUCTION(+:valx)
    do m = 1,ncs
      !
      val = rspec(m) - cpi(n,m)
      valx = valx + val*return_value(n,1,m)
      !
      output_variables(n,1+m) = -val
      output_variables(n,1+ncs+m) = dcpidt(n,m)
      !
    end do
    !
    ! Add the contributions of the species mass fraction gradients to the
    ! temperature gradient and then scale it by the proper leading coefficient
    ! computed above.
    !
    return_value(n,1,nee) = coeff*(dtdx + temp*valx)
    !
  end do
  !
  !$ACC END DATA
  deallocate ( dcpidt )
  deallocate ( cpi )
  deallocate ( mass_frac )
  !
end subroutine adjust_temp_grad_for_tpg_1d
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
!#define DEBUG_VORTEX_SOLUTION
#ifndef DISABLE_PURE
#ifndef DEBUG_VORTEX_SOLUTION
pure &
#endif
#endif
subroutine vortex_solution_gpu(n1,n2,xyz,return_value)
  !
  !.. Use Statements ..
  use ovar, only : velref,aref,uref,machref
  use ovar, only : gam,rhoref,tref,pref,itestcase
  use ovar, only : vortex_bigR,vortex_bigGam
  use ovar, only : itcur,Period_Timesteps,Period_Time
  use ovar, only : Grid_Length,Grid_min
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
  real(wp) :: vx_ref,vy_ref,delt_since_lastperiod
  real(wp) :: a,du,dv,dT,vx,vy,vz,q2,temp,rho,pres
  real(wp) :: x,y,dx_prop,dy_prop,xc,yc,dx,dy
  !
  logical(lk) :: running_vortex_problem
  !
continue
  !
  nr = size(xyz,dim=1)
  !
  ! Compute the time that has elapsed since the start of the current period
  !
  delt_since_lastperiod = real(mod(itcur,Period_Timesteps),wp) &
                        / real(Period_Timesteps,wp)
  delt_since_lastperiod = delt_since_lastperiod * Period_Time
  !
  ! Find the x and y distances the vortex has
  ! traveled/propagated since the start of the current period
  !
  dx_prop = deltafun(velref(1)>eps3) * (velref(1)/aref) * delt_since_lastperiod
  dy_prop = deltafun(velref(2)>eps3) * (velref(2)/aref) * delt_since_lastperiod
  !
  ! Now translate these into a periodic domain
  ! with the grid lengths for each direction
  !
  dx_prop = mod(dx_prop,Grid_Length(1))
  dy_prop = mod(dy_prop,Grid_Length(2))
  !
  ! Some needed constants
  !
  ! running_vortex_problem = fals -> advection only problem
  ! running_vortex_problem = true -> actual vortex problem
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
#ifdef DEBUG_VORTEX_SOLUTION
  write (iout,1) "dx_prop",dx_prop
  write (iout,1) "dy_prop",dy_prop
  write (iout,1) "delt",delt_since_lastperiod
  write (iout,1) "gam",gam
  write (iout,1) "rhoref",rhoref
  write (iout,1) "aref",aref
  write (iout,1) "vx_ref",vx_ref
  write (iout,1) "vy_ref",vy_ref
  write (iout,1) "tref",tref
  write (iout,1) "rhoref_inv",rhoref_inv
  write (iout,1) "aref_inv",aref_inv
  write (iout,1) "rhoa2_inv",rhoa2_inv
  write (iout,1) "tref_inv",tref_inv
  write (iout,1) "bigR",bigR
  write (iout,1) "bigGam",bigGam
  write (iout,1) "c0",c0
  write (iout,1) "c1",c1
  write (iout,1) "c2",c2
  write (iout,1) "c3",c3
  write (iout,1) "c4",c4
  write (iout,1) "c5",c5
  1 format (2x,a,t12," = ",es14.6)
#endif
  !
  !$ACC PARALLEL LOOP PRESENT(xyz,return_value)
  do n = n1,n2
    !
    xc = xyz(1,n)
    yc = xyz(2,n)
    !
    dx = dx_prop
    dy = dy_prop
    !
    x = xc - merge( dx , dx - Grid_Length(1) , xc - dx >= Grid_Min(1) )
    y = yc - merge( dy , dy - Grid_Length(2) , yc - dy >= Grid_Min(2) )
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
    rho  = rhoref * (temp*tref_inv)**(gm1_inv)
    pres = pref   * (temp*tref_inv)**(gam_over_gm1)
    !
    ! Store the nondimensional primitive variables in return_value
    !
    return_value(n,nec) = rho * rhoref_inv
    return_value(n,nmx) = vx * aref_inv
    return_value(n,nmy) = vy * aref_inv
    return_value(n,nee) = pres * rhoa2_inv
    !
  end do
  !
  ! Set the z-velocity to the reference z-velocity if this is 3D
  !
  if (nr == 3) then
    vz = velref(3)*aref_inv
    !$ACC PARALLEL LOOP PRESENT(return_value)
    do n = n1,n2
      return_value(n,nmz) = vz
    end do
  end if
  !
  ! Convert primitive variables to conserved variables
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(return_value)
    do n = n1,n2
      !
      rho  = return_value(n,nec)
      vx   = return_value(n,nmx)
      vy   = return_value(n,nmy)
      vz   = return_value(n,nmz)
      pres = return_value(n,nee)
      !
      q2 = vx*vx + vy*vy + vz*vz
      !
      return_value(n,nmx) = rho*vx
      return_value(n,nmy) = rho*vy
      return_value(n,nmz) = rho*vz
      return_value(n,nee) = pres*c5 + half*rho*q2
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(return_value)
    do n = n1,n2
      !
      rho  = return_value(n,nec)
      vx   = return_value(n,nmx)
      vy   = return_value(n,nmy)
      pres = return_value(n,nee)
      !
      q2 = vx*vx + vy*vy
      !
      return_value(n,nmx) = rho*vx
      return_value(n,nmy) = rho*vy
      return_value(n,nee) = pres*c5 + half*rho*q2
      !
    end do
  end if
  !
end subroutine vortex_solution_gpu
!
!###############################################################################
!
end module gpu_mod
