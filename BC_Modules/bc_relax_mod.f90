module bc_relax_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmb,nme,nee,nq
  use eqn_idx, only : ncs,npres,ntemp,ncp,ngamm,nenth,ng
  !
  implicit none
  !
  private
  !
  public :: get_relax_coefficient
  public :: get_relaxed_cv_bc
  public :: get_relaxed_pv_bc
  public :: get_relaxed_total_bc
  !
contains
!
!###############################################################################
!
#ifdef DISABLE_FOR_GPU
pure &
#endif
subroutine get_relax_coefficient(relax,rlx)
  !
  !.. Use Statements ..
  use ovar, only : relaxation_t
  use ovar, only : itcur,time
  !
  !.. Formal Arguments ..
  type(relaxation_t), intent(in) :: relax
  real(wp),        intent(inout) :: rlx
  !
  !.. Local Scalars ..
  real(wp) :: rlx_iter,rlx_time
  !
continue
  !
  ! Determine the relaxation coefficient
  !
  rlx = one
  rlx_iter = one
  rlx_time = one
  if (relax%iter > 0) then
    rlx_iter = min( one , real(itcur,kind=wp)/real(relax%iter,kind=wp) )
  end if
  if (relax%time > zero) then
    rlx_time = min( one , time/relax%time )
  end if
  rlx = max( zero , min(relax%value,rlx_iter,rlx_time,one) )
  !
end subroutine get_relax_coefficient
!
!###############################################################################
!
#ifdef DISABLE_FOR_GPU
!pure &
#endif
subroutine get_relaxed_cv_bc(relax,cv_fs,cv_bc,gp_bc, &
                             cv_rlx,pv_rlx,gp_rlx,return_rlx)
  !
  !.. Use Statements ..
  use ovar, only : relaxation_t
  !
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Formal Arguments ..
  type(relaxation_t), intent(in) :: relax
  real(wp),           intent(in) :: cv_fs(:)
  real(wp),           intent(in) :: cv_bc(:)
  real(wp),           intent(in) :: gp_bc(:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(inout) :: cv_rlx(:)
  real(wp), optional, intent(inout) :: pv_rlx(:)
  real(wp), optional, intent(inout) :: gp_rlx(:)
  real(wp), optional, intent(inout) :: return_rlx
  !
  !.. Local Scalars ..
  integer  :: m
  real(wp) :: rlx
  !
  !.. Local Arrays ..
  real(wp) :: cv(1:size(cv_bc))
  real(wp) :: pv(1:size(cv_bc))
  real(wp) :: gp(1:size(gp_bc))
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_relaxed_cv_bc"
  !
continue
  !
  ! Determine the relaxation coefficient
  !
  call get_relax_coefficient(relax,rlx)
  !
  do m = 1,size(cv_bc)
    cv(m) = (one-rlx)*cv_fs(m) + rlx*cv_bc(m)
  end do
  !
  gp(npres) = gp_bc(npres)
  gp(ntemp) = gp_bc(ntemp)
  !
  call update_thermo_properties(cv,gp,caller_name=pname)
  call update_transport_properties(cv,gp)
  !
  if (ncs > 0) then
    do m = 1,ncs
      pv(m) = cv(m)/cv(nec)
    end do
  end if
  !
  pv(nec) = cv(nec)
  do m = nmb,nme
    pv(m) = cv(m)/cv(nec)
  end do
  pv(nee) = gp(npres)
  !
  if (present(cv_rlx)) then
    cv_rlx(:) = cv(:)
  end if
  !
  if (present(pv_rlx)) then
    pv_rlx(:) = pv(:)
   !pv_rlx(nec) = rtot_rlx
   !pv_rlx(nee) = ptot_rlx
  end if
  !
  if (present(gp_rlx)) then
    gp_rlx(:) = gp(:)
  end if
  !
  if (present(return_rlx)) then
    return_rlx = rlx
  end if
  !
end subroutine get_relaxed_cv_bc
!
!###############################################################################
!
#ifdef DISABLE_FOR_GPU
!pure &
#endif
subroutine get_relaxed_pv_bc(relax,pv_fs,pv_bc, &
                             cv_rlx,pv_rlx,gp_rlx,return_rlx)
  !
  !.. Use Statements ..
  use ovar, only : relaxation_t
  use ovar, only : gam,rgasref_nd
  !
  use thermo_mod,    only : compute_mixture_rgas
  use thermo_mod,    only : compute_mixture_enthalpy
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Formal Arguments ..
  type(relaxation_t), intent(in) :: relax
  real(wp),           intent(in) :: pv_fs(:)
  real(wp),           intent(in) :: pv_bc(:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(inout) :: cv_rlx(:)
  real(wp), optional, intent(inout) :: pv_rlx(:)
  real(wp), optional, intent(inout) :: gp_rlx(:)
  real(wp), optional, intent(inout) :: return_rlx
  !
  !.. Local Scalars ..
  integer  :: m
  real(wp) :: rlx
  real(wp) :: rho,pres,temp,enth,q2,rgas
  !
  !.. Local Arrays ..
  real(wp) :: cv(1:size(pv_bc))
  real(wp) :: pv(1:size(pv_bc))
  real(wp) :: gp(1:ng)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_relaxed_pv_bc"
  !
continue
  !
  ! Determine the relaxation coefficient
  !
  call get_relax_coefficient(relax,rlx)
  !
  ! Compute the relaxation boundary static density and static pressure
  !
  rho  = (one-rlx)*pv_fs(nec) + rlx*pv_bc(nec)
  pres = (one-rlx)*pv_fs(nee) + rlx*pv_bc(nee)
  !
  cv(nec) = rho
  pv(nec) = rho
  pv(nee) = pres
  !
  ! Copy the boundary velocities, without any relaxation
  !
  do m = nmb,nme
    pv(m) = pv_bc(m)
    cv(m) = pv(m)*rho
  end do
  q2 = dot_product( pv(nmb:nme) , pv(nmb:nme) )
  !
  ! Compute the relaxation boundary species mass fractions
  !
  if (ncs > 0) then
    do m = 1,ncs
      pv(m) = (one-rlx)*pv_fs(m) + rlx*pv_bc(m)
      cv(m) = pv(m)*rho
    end do
  end if
  !
  if (ncs > 0) then
    !
    call compute_mixture_rgas(pv(1:ncs),rgas)
    !
    temp = pres / (rgas*rho)
    !
    call compute_mixture_enthalpy(pv(1:ncs),temp,enth)
    !
    cv(nee) = rho*enth - pres + half*rho*q2
    !
  else
    !
    temp = pres / (rgasref_nd*rho)
    !
    cv(nee) = pres/(gam-one) + half*rho*q2
    !
  end if
  !
  gp(npres) = pres
  gp(ntemp) = temp
  !
  call update_thermo_properties(cv,gp,caller_name=pname)
  call update_transport_properties(cv,gp)
  !
  if (present(cv_rlx)) then
    cv_rlx(:) = cv(:)
  end if
  !
  if (present(pv_rlx)) then
    pv_rlx(:) = pv(:)
   !pv_rlx(nec) = rtot_rlx
   !pv_rlx(nee) = ptot_rlx
  end if
  !
  if (present(gp_rlx)) then
    gp_rlx(:) = gp(:)
  end if
  !
  if (present(return_rlx)) then
    return_rlx = rlx
  end if
  !
end subroutine get_relaxed_pv_bc
!
!###############################################################################
!
#ifdef DISABLE_FOR_GPU
!#define DEBUG_RELAXED_TOTAL
pure &
#endif
subroutine get_relaxed_total_bc(relax,pv_fs,gp_fs,pv_bc,gp_bc, &
                                cv_rlx,pv_rlx,gp_rlx,return_rlx)
  !
  !.. Use Statements ..
  use ovar, only : relaxation_t
  use ovar, only : gam,rgasref_nd
  use ovar, only : tref,tmin,tmax,conv_temp
  !
  use thermo_mod, only : compute_mixture_rgas
  use thermo_mod, only : compute_mixture_gamma
  use thermo_mod, only : compute_mixture_enthalpy
  use thermo_mod, only : compute_mixture_cp
  use thermo_mod, only : compute_cp_over_temp_integral
  !
  !.. Formal Arguments ..
  type(relaxation_t), intent(in) :: relax
  real(wp),           intent(in) :: pv_fs(:)
  real(wp),           intent(in) :: gp_fs(:)
  real(wp),           intent(in) :: pv_bc(:)
  real(wp),           intent(in) :: gp_bc(:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(inout) :: cv_rlx(:)
  real(wp), optional, intent(inout) :: pv_rlx(:)
  real(wp), optional, intent(inout) :: gp_rlx(:)
  real(wp), optional, intent(inout) :: return_rlx
  !
  !.. Local Scalars ..
  integer  :: m,k,i
  real(wp) :: rlx,tmin_nd,tmax_nd
  real(wp) :: total_cond,rtot_rlx,ptot_rlx,ttot_rlx,htot_rlx
  real(wp) :: aspd_fs,mach_fs,rtot_fs,ptot_fs
  real(wp) :: aspd_bc,mach_bc,rtot_bc,ptot_bc
  real(wp) :: gamm,gm1,one_over_gm1,gam_over_gm1,one_over_1mg,gam_over_1mg
  real(wp) :: func,dfunc,dtemp,dterr,val,aspd,vmag
  real(wp) :: rgas,q2,rho,pres,temp,enth,cp,htot
  !
  !.. Local Arrays ..
  real(wp) :: cv_bc(1:size(pv_bc))
  real(wp) :: cv(1:size(pv_bc))
  real(wp) :: pv(1:size(pv_bc))
  real(wp) :: gp(1:size(gp_bc))
  real(wp) :: dir_cos(nmb:nme)
  !
  !.. Local Parameters ..
  integer,  parameter :: newton_iter = 100
 !real(wp), parameter :: conv_temp = eps6
  !
continue
  !
  ! Determine the relaxation coefficient
  !
  call get_relax_coefficient(relax,rlx)
  !
  q2 = dot_product( pv_bc(nmb:nme) , pv_bc(nmb:nme) )
  !
  cv_bc(nec) = pv_bc(nec)
  cv_bc(nmb:nme) = pv_bc(nmb:nme)*pv_bc(nec)
  if (ncs > 0) then
    cv_bc(1:ncs) = pv_bc(1:ncs)*pv_bc(nec)
    cv_bc(nee) = pv_bc(nec)*gp_bc(nenth) - gp_bc(npres) + half*pv_bc(nec)*q2
  else
    cv_bc(nee) = gp_bc(npres)/(gam-one) + half*pv_bc(nec)*q2
  end if
#ifdef DEBUG_RELAXED_TOTAL
  write (iout,2) (i,pv_bc(i),i=1,nq)
  write (iout,3) (i,cv_bc(i),i=1,nq)
  write (iout,4) (i,gp_bc(i),i=1,ng)
  write (iout,1) "rlx",rlx
  1 format (4x,a," = ",es14.6)
  2 format (20(4x,"pv_bc(",i0,") = ",es14.6,:,/))
  3 format (20(4x,"cv_bc(",i0,") = ",es14.6,:,/))
  4 format (20(4x,"gp_bc(",i0,") = ",es14.6,:,/))
#endif
  !
  ! Compute the inflow conditions
  !
  rho = pv_bc(nec)
  pres = gp_bc(npres)
  !
  if (ncs > 0) then
    gamm = gp_bc(ngamm)
  else
    gamm = gam
  end if
  !
  gm1 = gamm - one
  one_over_gm1 = one / gm1
  gam_over_gm1 = gamm * one_over_gm1
  !
  aspd_bc = sqrt( max( gamm*pres/rho , rndoff ) )
  !
  q2 = dot_product( pv_bc(nmb:nme) , pv_bc(nmb:nme) )
  !
  dir_cos(nmb:nme) = pv_bc(nmb:nme)/max(sqrt(q2),rndoff)
  !
  mach_bc = sqrt(q2) / aspd_bc
  !
  total_cond = one + half*gm1*mach_bc*mach_bc
  !
  rtot_bc = rho  * total_cond**one_over_gm1
  ptot_bc = pres * total_cond**gam_over_gm1
  !
#ifdef DEBUG_RELAXED_TOTAL
  write (iout,1) "gamm",gamm
  write (iout,1) "aspd_bc",aspd_bc
  write (iout,1) "pv_bc(nec)",pv_bc(nec)
  write (iout,1) "pv_bc(nee)",pv_bc(nee)
  write (iout,1) "gp_bc(npres)",gp_bc(npres)
  write (iout,1) "total_cond",total_cond
  write (iout,1) "rtot_bc",rtot_bc
  write (iout,1) "ptot_bc",ptot_bc
  5 format (20(4x,"   pv(",i0,") = ",es14.6,:,/))
#endif
  !
  if (abs(rlx-one) < eps6) then
    !
    cv(:) = cv_bc(:)
    pv(:) = pv_bc(:)
    gp(:) = gp_bc(:)
    !
    ptot_rlx = ptot_bc
    rtot_rlx = rtot_bc
    !
  else
    !
    ! Compute the free stream conditions
    !
    if (ncs > 0) then
      gamm = gp_fs(ngamm)
    else
      gamm = gam
      rgas = rgasref_nd
    end if
    !
    gm1 = gamm - one
    one_over_gm1 = one / gm1
    gam_over_gm1 = gamm * one_over_gm1
    one_over_1mg = -one_over_gm1
    gam_over_1mg = -gam_over_gm1
    !
    aspd_fs = sqrt( max( gamm*pv_fs(nee)/pv_fs(nec) , rndoff ) )
    !
    mach_fs = norm2( pv_fs(nmb:nme) ) / aspd_fs
    !
    total_cond = one + half*gm1*mach_fs*mach_fs
    !
    rtot_fs = pv_fs(nec) * total_cond**one_over_gm1
    ptot_fs = pv_fs(nee) * total_cond**gam_over_gm1
    !
#ifdef DEBUG_RELAXED_TOTAL
    write (iout,1) "gamm",gamm
    write (iout,1) "aspd_fs",aspd_fs
    write (iout,1) "pv_fs(nec)",pv_fs(nec)
    write (iout,1) "pv_fs(nee)",pv_fs(nee)
    write (iout,1) "gp_fs(npres)",gp_fs(npres)
    write (iout,1) "total_cond",total_cond
    write (iout,1) "rtot_fs",rtot_fs
    write (iout,1) "ptot_fs",ptot_fs
#endif
    !
    ! Compute the relaxation boundary total density and total pressure
    !
    rtot_rlx = (one-rlx)*rtot_fs + rlx*rtot_bc
    ptot_rlx = (one-rlx)*ptot_fs + rlx*ptot_bc
    !
    if (ncs > 0) then
      !
      tmin_nd = tmin/tref
      tmax_nd = tmax/tref
      !
      do m = 1,ncs
        pv(m) = (one-rlx)*pv_fs(m) + rlx*pv_bc(m)
      end do
      call compute_mixture_rgas(pv(1:ncs),rgas)
      !
      ttot_rlx = ptot_rlx/(rgas*rtot_rlx)
      !
      call compute_mixture_enthalpy(pv(1:ncs),ttot_rlx,htot_rlx)
      !
      temp = ttot_rlx/total_cond
      !
      newton_iteration: &
      do k = 1,newton_iter
        !
        call compute_mixture_enthalpy(pv(1:ncs),temp,enth)
        !
        call compute_mixture_cp(pv(1:ncs),temp,cp)
        !
        gamm = rgas
        !
        call compute_mixture_gamma(cp,gamm)
        !
        htot = enth + half * gamm*rgas*temp * mach_bc*mach_bc
        !
        func = htot - htot_rlx
        dfunc = cp + half * gamm*rgas * mach_bc*mach_bc
        !
        dtemp = func/dfunc
        !
        temp = min( tmax_nd , max(tmin_nd,temp-dtemp) )
        !
        dterr = abs(dtemp)/temp
        !
        if (dterr < conv_temp) then
          exit newton_iteration
        end if
        !
      end do newton_iteration
      !
      call compute_cp_over_temp_integral(pv(1:ncs),val,temp,ttot_rlx)
      !
      pres = ptot_rlx / exp(val/rgas)
      rho = pres / (rgas*temp)
      !
      call compute_mixture_enthalpy(pv(1:ncs),temp,enth)
      !
      call compute_mixture_cp(pv(1:ncs),temp,cp)
      !
      gamm = rgas
      call compute_mixture_gamma(cp,gamm)
      !
      aspd = sqrt(max(gamm*rgas*temp,rndoff))
      vmag = mach_bc*aspd
      !
      do m = 1,ncs
        cv(m) = pv(m)*rho
      end do
      !
      cv(nec) = rho
      pv(nec) = rho
      do m = nmb,nme
        pv(m) = vmag*dir_cos(m)
        cv(m) = rho*pv(m)
      end do
      pv(nee) = pres
      cv(nee) = rho*enth - pres + half*rho*vmag*vmag
      !
      gp(npres) = pres
      gp(ntemp) = temp
      gp(ngamm) = gamm
      gp(nenth) = enth
      gp(ncp)   = cp
      !
    else
      !
      ttot_rlx = ptot_rlx/(rgasref_nd*rtot_rlx)
      !
      rho  = rtot_rlx * total_cond**one_over_1mg
      pres = ptot_rlx * total_cond**gam_over_1mg
      temp = ttot_rlx/total_cond
      !
      aspd = sqrt(max(gam*rgasref_nd*temp,rndoff))
      !
      vmag = mach_bc*aspd
      !
      cv(nec) = rho
      pv(nec) = rho
      do m = nmb,nme
        pv(m) = vmag*dir_cos(m)
        cv(m) = rho*pv(m)
      end do
      pv(nee) = pres
      cv(nee) = pres*one_over_gm1 + half*rho*vmag*vmag
      !
      gp(npres) = pres
      gp(ntemp) = temp
      !
    end if
    !
  end if
  !
#ifdef DEBUG_RELAXED_TOTAL
  write (iout,5) (i,pv(i), i=1,nq)
  write (iout,1) "rtot_rlx",rtot_rlx
  write (iout,1) "ptot_rlx",ptot_rlx
 !write (iout,1) "ttot_rlx",ttot_rlx
#endif
  !
  if (present(cv_rlx)) then
    cv_rlx(:) = cv(:)
  end if
  !
  if (present(pv_rlx)) then
    pv_rlx(:) = pv(:)
    pv_rlx(nec) = rtot_rlx
    pv_rlx(nee) = ptot_rlx
  end if
  !
  if (present(gp_rlx)) then
    gp_rlx(:) = gp(:)
  end if
  !
  if (present(return_rlx)) then
    return_rlx = rlx
  end if
  !
end subroutine get_relaxed_total_bc
!
!###############################################################################
!
end module bc_relax_mod
