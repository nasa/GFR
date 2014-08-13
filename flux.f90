module flux_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use ovar, only : flux_method
  !
  implicit none
  !
  private
  !
 !include "kind_defs.h"
 !include "constants.h"
  !
  public :: correctd
  public :: correct_flux_project_cv
  public :: correct_flux_project_pv
  public :: correct_flux_project_flx
  public :: derflux1
  public :: derflux2
  public :: derflux3
  public :: dg_weak_form
  public :: dg_strong_form
  public :: interp_cv_to_flxpts
  public :: interp_dusp_to_flxpts
  public :: derflux1_qdtr
  public :: derflux1_qdtr_with_visc
  public :: correct_interface_flux_qdtr
  public :: correct_interface_flux_qdtr_with_visc
  public :: discontinuous_gradient_at_solpts
  public :: correct_gradients_at_solpts
  !
  real(wp), parameter :: amco2m = 1.0e-4_wp
  real(wp), parameter :: alphax = 0.1875_wp
  real(wp), parameter :: betax  = 0.125_wp
  !
  ! CL : coefficient for the length of the face. This is needed because of
  !      the difference between how the reference triangle and reference
  !      quad are defined. 
  !      = 1.0 for triangles
  !      = 0.5 for quads
  !
  real(wp), dimension(1:2), parameter :: cl = (/ one , half /)
  !
  ! Generalization parameters for the diffusion scheme
  !
  ! LDG_BETA : Constant that determines the directional bias of
  !            the diffusion scheme for each coordinate direction
  !        = [0,0,0] -> Results in BR2, a centered biased scheme
  !   
  real(wp), parameter :: LDG_beta(3) = (/ zero, zero, zero /)
  !
  ! LDG_TAU : Constant that controls the penalty term
  !           term of the common diffusive flux
  !        = 
  real(wp), parameter :: LDG_tau = zero

  !
  integer, parameter :: ROEFLUX = 1
  integer, parameter :: HLLC = 2
  integer, parameter :: LDFSS = 3
  !
  logical, parameter :: output_boundary_fluxes = fals
  logical, save :: first_output_done = fals
  !
contains
!
!###############################################################################
!
pure function upwind_flux_cv(unrm,cvl,cvr) result(return_value)
  !
  !... input unit normal unx,uny, edge length ds, rl,ul,vl,pl, rr,ur,vr,pr
  !... output: f1,f2,f3,f4 : fluxes.
  !... Roe's flux with Osher's simple wave model (Huynh, SIAM J. Num. 1995)
  !
  use geovar, only : neqns,nr
  use ovar, only : gam
  !
  real(wp), dimension(1:nr), intent(in) :: unrm
  real(wp), dimension(1:neqns), intent(in) :: cvl
  real(wp), dimension(1:neqns), intent(in) :: cvr
  !
  real(wp), dimension(1:neqns) :: return_value
  !
  real(wp) :: gm1,gp1,unx,uny
  real(wp) :: rl,ul,vl,pl,hl,ql,rql,betal
  real(wp) :: rr,ur,vr,pr,hr,qr,rqr,betar
  real(wp) :: dw1,hsr1,eta1,udw1
  real(wp) :: dw4,hsr4,eta4,udw4
  real(wp) :: rtd,utd,vtd,htd,atd2,atd,qtd
  real(wp) :: atd2nz,atd2rt,atdnz,atdrt
  !
continue
  !
  gm1 = gam - one
  gp1 = gam + one
  !
  unx = unrm(1)
  uny = unrm(2)
  !
  rl = cvl(1)
  ul = cvl(2)/cvl(1)
  vl = cvl(3)/cvl(1)
  pl = gm1*(cvl(4) - half*rl*(ul**2+vl**2))
  !
  rr = cvr(1)
  ur = cvr(2)/cvr(1)
  vr = cvr(3)/cvr(1)
  pr = gm1*(cvr(4) - half*rr*(ur**2+vr**2))
  !
  hl = half*(ul**2+vl**2) + gam/gm1 * pl/rl
  hr = half*(ur**2+vr**2) + gam/gm1 * pr/rr
  ql = ul*unx + vl*uny
  qr = ur*unx + vr*uny
  !
  !... square root averaging
  rtd   = sqrt(rr*rl)
  betal = rl / (rl+rtd)
  betar = one - betal
  utd   = betal*ul + betar*ur
  vtd   = betal*vl + betar*vr
  htd   = betal*hl + betar*hr
  atd2  = max(gm1*(htd - half*(utd**2+vtd**2)),zero)
  atd   = sqrt(atd2)
  qtd   = utd*unx + vtd*uny
  !
  atd2nz = max(atd2,rndoff)
  atdnz = max(atd,rndoff)
  atd2rt = atd2 / atd2nz
  atdrt = atd / atdnz
  !
  if (qtd >= zero) then
    !
    dw1  = half*((pr-pl)/atd2nz*atd2rt - (qr-ql)*rtd/atdnz*atdrt)
    hsr1 = -c0_25*gp1*dw1*atd/rtd
    eta1 = max(-abs(qtd-atd) + hsr1, zero)
    udw1 = dw1*(min(qtd-atd,zero) - half*eta1)
    rql  = rl*ql
    !
    return_value(1) =    rql + udw1
    return_value(2) = ul*rql + udw1*(utd-atd*unx) + pl*unx
    return_value(3) = vl*rql + udw1*(vtd-atd*uny) + pl*uny
    return_value(4) = hl*rql + udw1*(htd-atd*qtd)
    !
  else
    !
    dw4  = half*((pr-pl)/atd2nz*atd2rt + (qr-ql)*rtd/atdnz*atdrt)
    hsr4 = c0_25*gp1*dw4*atd/rtd
    eta4 = max(-abs(qtd+atd) + hsr4, zero)
    udw4 = dw4*(max(qtd+atd,zero) + half*eta4)
    rqr  = rr*qr
    !
    return_value(1) =    rqr - udw4
    return_value(2) = ur*rqr - udw4*(utd+atd*unx) + pr*unx
    return_value(3) = vr*rqr - udw4*(vtd+atd*uny) + pr*uny
    return_value(4) = hr*rqr - udw4*(htd+atd*qtd)
    !
  end if
  !
end function upwind_flux_cv
!
!###############################################################################
!
pure function upwind_flux_pv(unrm,pvl,pvr) result(return_value)
  !
  !... input unit normal unx,uny, edge length ds, rl,ul,vl,pl, rr,ur,vr,pr
  !... output: f1,f2,f3,f4 : fluxes.
  !... Roe's flux with Osher's simple wave model (Huynh, SIAM J. Num. 1995)
  !
  use geovar, only : neqns,nr
  use ovar, only : gam
  !
  real(wp), dimension(1:nr), intent(in) :: unrm
  real(wp), dimension(1:neqns), intent(in) :: pvl
  real(wp), dimension(1:neqns), intent(in) :: pvr
  !
  real(wp), dimension(1:neqns) :: return_value
  !
  real(wp) :: gm1,gp1,unx,uny
  real(wp) :: rl,ul,vl,pl,hl,ql,rql,betal
  real(wp) :: rr,ur,vr,pr,hr,qr,rqr,betar
  real(wp) :: dw1,hsr1,eta1,udw1
  real(wp) :: dw4,hsr4,eta4,udw4
  real(wp) :: rtd,utd,vtd,htd,atd2,atd,qtd
  !
continue
  !
  gm1 = gam - one
  gp1 = gam + one
  !
  unx = unrm(1)
  uny = unrm(2)
  !
  rl = pvl(1)
  ul = pvl(2)
  vl = pvl(3)
  pl = pvl(4)
  !
  rr = pvr(1)
  ur = pvr(2)
  vr = pvr(3)
  pr = pvr(4)
  !
  hl = half*(ul**2+vl**2) + gam/gm1 * pl/rl
  hr = half*(ur**2+vr**2) + gam/gm1 * pr/rr
  ql = ul*unx + vl*uny
  qr = ur*unx + vr*uny
  !
  !... square root averaging
  rtd   = sqrt(rr*rl)
  betal = rl / (rl+rtd)
  betar = one - betal
  utd   = betal*ul + betar*ur
  vtd   = betal*vl + betar*vr
  htd   = betal*hl + betar*hr
  atd2  = gm1*(htd - half*(utd**2+vtd**2))
  atd   = sqrt(atd2)
  qtd   = utd*unx + vtd*uny
  !
  if (qtd >= zero) then
    !
    dw1  = half*((pr-pl)/atd2 - (qr-ql)*rtd/atd)
    hsr1 = -c0_25*gp1*dw1*atd/rtd
    eta1 = max(-abs(qtd-atd) + hsr1, zero)
    udw1 = dw1*(min(qtd-atd,zero) - half*eta1)
    rql  = rl*ql
    !
    return_value(1) =    rql + udw1
    return_value(2) = ul*rql + udw1*(utd-atd*unx) + pl*unx
    return_value(3) = vl*rql + udw1*(vtd-atd*uny) + pl*uny
    return_value(4) = hl*rql + udw1*(htd-atd*qtd)
    !
  else
    !
    dw4  = half*((pr-pl)/atd2 + (qr-ql)*rtd/atd)
    hsr4 = c0_25*gp1*dw4*atd/rtd
    eta4 = max(-abs(qtd+atd) + hsr4, zero)
    udw4 = dw4*(max(qtd+atd,zero) + half*eta4)
    rqr  = rr*qr
    !
    return_value(1) =    rqr - udw4
    return_value(2) = ur*rqr - udw4*(utd+atd*unx) + pr*unx
    return_value(3) = vr*rqr - udw4*(vtd+atd*uny) + pr*uny
    return_value(4) = hr*rqr - udw4*(htd+atd*qtd)
    !
  end if
  !
end function upwind_flux_pv
!
!###############################################################################
!
pure function hllc_upwind_flux_cv(unrm,cvl,cvr) result(return_value)
  !
  !... input unit normal unx,uny, edge length ds, rl,ul,vl,pl, rr,ur,vr,pr
  !... output: f1,f2,f3,f4 : fluxes.
  !... Roe's flux with Osher's simple wave model (Huynh, SIAM J. Num. 1995)
  !
  use geovar, only : neqns,nr
  use ovar, only : gam
  !
  real(wp), dimension(1:nr), intent(in) :: unrm
  real(wp), dimension(1:neqns), intent(in) :: cvl
  real(wp), dimension(1:neqns), intent(in) :: cvr
  !
  !.. Function Return Value ..
  real(wp), dimension(1:neqns) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: gm1,unx,uny
  real(wp) :: rl,ul,vl,pl,hl,al,vnl,rhoel,wgl,sl,dsl
  real(wp) :: rr,ur,vr,pr,hr,ar,vnr,rhoer,wgr,sr,dsr
  real(wp) :: srt,ra,ua,va,ha,q2a,vna,aa
  real(wp) :: rho,rhou,rhov,rhoe,pstar,sm
  !
continue
  !
  gm1 = gam - one
  !
  unx = unrm(1)
  uny = unrm(2)
  !
  rl = cvl(1)
  ul = cvl(2)/cvl(1)
  vl = cvl(3)/cvl(1)
  pl = gm1*(cvl(4) - half*rl*(ul**2+vl**2))
  hl = half*(ul**2+vl**2) + gam/gm1 * pl/rl
  vnl = ul*unx + vl*uny
  rhoel = cvl(4)
  !
  rr = cvr(1)
  ur = cvr(2)/cvr(1)
  vr = cvr(3)/cvr(1)
  pr = gm1*(cvr(4) - half*rr*(ur**2+vr**2))
  hr = half*(ur**2+vr**2) + gam/gm1 * pr/rr
  vnr = ur*unx + vr*uny
  rhoer = cvr(4)
  !
  srt = sqrt(max(rr/rl,rndoff))
  wgl = one/(srt+one)
  wgr = srt/(srt+one)
  !
  ra = sqrt(max(rl*rr,rndoff))
  ua = ul*wgl + ur*wgr
  va = vl*wgl + vr*wgr
  ha = hl*wgl + hr*wgr
  !
  q2a = ua**2 + va**2
  vna = ua*unx + va*uny
  !
  aa = sqrt(max(gm1*(ha-half*q2a),rndoff))
  !
  sl = min( vnl-al , vna-aa )
  sr = max( vnr+ar , vna+aa )
  !
  if (sl > zero) then
    !
    return_value(1) = vnl*rl
    return_value(2) = vnl*rl*ul + pl*unx
    return_value(3) = vnl*rl*vl + pl*uny
    return_value(4) = vnl*rhoel + pl*vnl
    !
  else if (sr < zero) then
    !
    return_value(1) = vnr*rr
    return_value(2) = vnr*rr*ur + pr*unx
    return_value(3) = vnr*rr*vr + pr*uny
    return_value(4) = vnr*rhoer + pr*vnr
    !
  else
    !
    dsl = sl - vnl
    dsr = sr - vnr
    !
    sm = ((rr*vnr*dsr-pr) - (rl*vnl*dsl-pl)) / (rr*dsr-rl*dsl)
    !
    pstar = pl*dsl*(sm-vnl) + pl
    !
    if (sm > zero) then
      !
      rho = dsl*rl/(sl-sm)
      rhou = (dsl*rl*ul + (pstar-pl)*unx)/(sl-sm)
      rhov = (dsl*rl*vl + (pstar-pl)*uny)/(sl-sm)
      rhoe = (dsl*rhoel - pl*vnl + pstar*sm)/(sl-sm)
      !
      return_value(1) = sm*rho
      return_value(2) = sm*rhou + pstar*unx
      return_value(3) = sm*rhov + pstar*uny
      return_value(4) = sm*rhoe + pstar*sm
      !
    else
      !
      rho = dsr*rr/(sr-sm)
      rhou = (dsr*rr*ur + (pstar-pr)*unx)/(sr-sm)
      rhov = (dsr*rr*vr + (pstar-pr)*uny)/(sr-sm)
      rhoe = (dsr*rhoer - pr*vnr + pstar*sm)/(sr-sm)
      !
      return_value(1) = sm*rho
      return_value(2) = sm*rhou + pstar*unx
      return_value(3) = sm*rhov + pstar*uny
      return_value(4) = sm*rhoe + pstar*sm
      !
    end if
    !
  end if
  !
end function hllc_upwind_flux_cv
!
!###############################################################################
!
pure function ldfss_upwind_flux_cv(unrm,cvl,cvr) result(return_value)
  !
  !... input unit normal unx,uny, edge length ds, rl,ul,vl,pl, rr,ur,vr,pr
  !... output: f1,f2,f3,f4 : fluxes.
  !... Roe's flux with Osher's simple wave model (Huynh, SIAM J. Num. 1995)
  !
  use geovar, only : neqns,nr
  use ovar, only : gam
  !
  real(wp), dimension(1:nr), intent(in) :: unrm
  real(wp), dimension(1:neqns), intent(in) :: cvl
  real(wp), dimension(1:neqns), intent(in) :: cvr
  !
  real(wp), dimension(1:neqns) :: return_value
  !
  real(wp) :: gm1,unx,uny,ahalf
  real(wp) :: rl,ul,vl,pl,hl,al
  real(wp) :: rr,ur,vr,pr,hr,ar
  real(wp) :: alpha_l,beta_l,mach_l,mstar_l
  real(wp) :: alpha_r,beta_r,mach_r,mstar_r
  real(wp) :: machp_half,uc_l,flux_l,pplus_l
  real(wp) :: machm_half,uc_r,flux_r,pplus_r
  real(wp) :: m_half,mstar_half
  real(wp) :: del_p,sum_p,p_net
  !
continue
  !
  gm1 = gam - one
  !
  unx = unrm(1)
  uny = unrm(2)
  !
  rl = cvl(1)
  ul = cvl(2)/cvl(1)
  vl = cvl(3)/cvl(1)
  pl = gm1*(cvl(4) - half*rl*(ul**2+vl**2))
  hl = half*(ul**2+vl**2) + gam/gm1 * pl/rl
  al = sqrt(max(gam*pl/rl,rndoff))
  !
  rr = cvr(1)
  ur = cvr(2)/cvr(1)
  vr = cvr(3)/cvr(1)
  pr = gm1*(cvr(4) - half*rr*(ur**2+vr**2))
  hr = half*(ur**2+vr**2) + gam/gm1 * pr/rr
  ar = sqrt(max(gam*pr/rr,rndoff))
  !
  ahalf = half*(al+ar)
  !
  mach_l = (unx*ul + uny*vl)/ahalf
  mach_r = (unx*ur + uny*vr)/ahalf
  !
  alpha_l = half*(one + sign(one,mach_l))
  alpha_r = half*(one - sign(one,mach_r))
  !
  beta_l = -max( zero , one-real(int(abs(mach_l)),kind=wp) )
  beta_r = -max( zero , one-real(int(abs(mach_r)),kind=wp) )
  !
  mstar_l =  one4*(mach_l + one)**2
  mstar_r = -one4*(mach_r - one)**2
  !
  m_half = sqrt( half*(mstar_l**2 + mstar_r**2) )
  !
  mstar_half = one4*beta_l*beta_r*(m_half - one)**2
  !
  del_p = pl - pr
  sum_p = pl + pr
  !
  machp_half = mstar_half*(one - del_p/sum_p - two*abs(del_p)/pl)
  machm_half = mstar_half*(one + del_p/sum_p - two*abs(del_p)/pr)
  !
  uc_l = ahalf*(alpha_l*(one+beta_l)*mach_l - beta_l*mstar_l - machp_half)
  uc_r = ahalf*(alpha_r*(one+beta_r)*mach_r - beta_r*mstar_r + machm_half)
  !
  flux_l = uc_l*rl
  flux_r = uc_r*rr
  !
  pplus_l = one4 * (two-mach_l) * (mach_l+one)**2
  pplus_r = one4 * (two+mach_r) * (mach_r-one)**2
  !
  p_net = pl*(alpha_l*(one+beta_l) - beta_l*pplus_l) + &
          pr*(alpha_r*(one+beta_r) - beta_r*pplus_r)
  !
  return_value(1) = flux_l    + flux_r
  return_value(2) = flux_l*ul + flux_r*ur + p_net*unx
  return_value(3) = flux_l*vl + flux_r*vr + p_net*uny
  return_value(4) = flux_l*hl + flux_r*hr
  !
end function ldfss_upwind_flux_cv
!
!###############################################################################
!
pure function ausm_plus_upwind_flux(unrm,pvl,pvr) result(return_value)
  !
  !   Compute Euler fluxes using AUSM+
  !... input unit normal unx,uny, edge length ds, rl,ul,vl,pl, rr,ur,vr,pr
  !... output: f1,f2,f3,f4 : fluxes
  !
  use geovar, only : neqns,nr
  use ovar, only : gam
  !
  real(wp), dimension(1:nr), intent(in) :: unrm
  real(wp), dimension(1:neqns), intent(in) :: pvl
  real(wp), dimension(1:neqns), intent(in) :: pvr
  !
  real(wp), dimension(1:neqns) :: return_value
  !
  real(wp) :: gm1,unx,uny
  real(wp) :: rl,ul,vl,pl,hl,al,v2l,a2l,h0l,uul,f1l
  real(wp) :: rr,ur,vr,pr,hr,ar,v2r,a2r,h0r,uur,f1r
  real(wp) :: alr,alr2,rm2,rmref2,coef1,coef2,prec_f
  real(wp) :: aml,absml,am1p,aml1,addp,aml2,ampl,ppl
  real(wp) :: amr,absmr,am1m,amr1,addm,amr2,ammr,pmr
  real(wp) :: amhalf,denom,amhalfp,amhalfm,phalf,amco2
  !
continue
  !
  gm1 = gam - one
  !
  unx = unrm(1)
  uny = unrm(2)
  !
  rl = pvl(1)
  ul = pvl(2)
  vl = pvl(3)
  pl = pvl(4)
  !
  rr = pvr(1)
  ur = pvr(2)
  vr = pvr(3)
  pr = pvr(4)
  !
  hl  = gam*pl/(gm1*rl)
  v2l = half*(ul**2 + vl**2)
  h0l = hl + v2l
  a2l = gm1*hl
  al  = sqrt(a2l)
  uul = unx*ul + uny*vl
  !
  hr  = gam*pr/(gm1*rr)
  v2r = half*(ur**2 + vr**2)
  h0r = hr + v2r
  a2r = gm1*hr
  ar  = sqrt(a2r)
  uur = unx*ur + uny*vr
  !
  alr    = half*(al + ar)
  alr2   = alr**2
  rm2    = (v2r+v2l)/alr2
  amco2  = max(amco2m,c0_25*rm2)
  rmref2 = min( one,max( rm2,amco2 ) )
  coef1  = half * (one+rmref2)
  coef2  = one - coef1
  prec_f = one - (one-sqrt(rmref2))**2
  alr    = prec_f*alr
  !
  aml   = (coef1*uul + coef2*uur)/alr
  amr   = (coef1*uur + coef2*uul)/alr
  absml = abs(aml)
  absmr = abs(amr)
  am1p  = half*(aml + absml)
  am1m  = half*(amr - absmr)
  !
  if (absml < one) then
    aml1 = aml + one
    addp = (aml1*(aml - one))**2
    aml2 = c0_25*aml1**2
    ampl = aml2 + betax*addp
    ppl  = aml2*(two - aml) + alphax*aml*addp
  else
    ampl = am1p
    ppl  = am1p/aml
  end if
  !
  if (absmr < one) then
    amr1 = amr - one
    addm = (amr1*(amr + one))**2
    amr2 = c0_25*amr1**2
    ammr = -amr2 - betax*addm
    pmr  = amr2*(two + amr) - alphax*amr*addm
  else
    ammr = am1m
    pmr  = am1m/amr
  end if
  !
  amhalf  = ampl + ammr
  denom   = ( rmref2*(pl+pr)*(rl+rr) - zero*(pl-pr)*(rl-rr) )/sqrt(rl*rr)
  amhalf  = amhalf + half*(pl-pr)*((ampl-ammr) - (am1p-am1m))/denom
  amhalfp = half*(amhalf + abs(amhalf))
  amhalfm = amhalf - amhalfp
  !
  phalf = ppl*pl + pmr*pr - ppl*pmr*(rl+rr)/two*alr*(uur-uul)
  !
  f1l = rl*amhalfp*alr
  f1r = rr*amhalfm*alr
  !
  return_value(1) = f1l     + f1r
  return_value(2) = f1l*ul  + f1r*ur  + unx*phalf
  return_value(3) = f1l*vl  + f1r*vr  + uny*phalf
  return_value(4) = f1l*h0l + f1r*h0r
  !
end function ausm_plus_upwind_flux
!
!###############################################################################
!
subroutine interp_cv_to_flxpts()
  !
  use geovar,  only : neqns,ncell,n_sp2dt,n_sp2dq,pts_in_cell_ptr
  use geovar,  only : flxpts_of_cell_ptr,flxpts_of_cell
  use geovar,  only : interp_s2f,interp_s2g,n_sp1de,n_fp1de
  use flowvar, only : usp,faceusp
  !
  use parallel_mod, only : exchange_faceusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,nl,k,m
  integer :: nf,f1,f2,kf,lf,sf
  !
  !.. Local Arrays ..
  real(wp), dimension(1:n_sp2dq) :: lmat
  real(wp), dimension(1:n_sp2dq,1:neqns) :: tusp
  real(wp), dimension(1:n_sp2dq,1:4*n_fp1de,1:2) :: interp_mat
  !
continue
  !
  !######################################################################
  !######################################################################
  ! WE NEED TO DO SOMETHING FOR THE COMMUNICATION AND PERIODIC BOUNDARY
  ! CONDITIONS TO INTERPOLATE THE GHOST CELL SOLUTION POINTS TO THE FLUX
  ! POINTS ON THE BOUNDARY FACE.
  !
  ! AFTER CALLING THIS ROUTINE, WE NEED TO EXCHANGE THE FACE DATA
  ! FOR COMMUNICATION AND PERIODIC BOUNDARY CONDITIONS USING MPI
  ! NON-BLOCKING ROUTINES. THIS WAY WE CAN INITIALIZE THIS COMMUNICATION
  ! AND THEN CONTINUE ON TO COMPUTING THE DISCONTINUOUS FLUX DERIVATIVES
  ! AT THE INTERNAL SOLUTION POINTS WHILE THE COMMUNICATION FINISHES.
  ! THIS CAN BE DONE BECAUSE WE DONT NEED THE INTERFACE INFORMATION UNTIL
  ! AFTER THAT STEP.
  !######################################################################
  !######################################################################
  !
  if (n_fp1de /= n_sp1de) then
    interp_mat(:,:,:) = interp_s2g(:,:,:)
  else
    interp_mat(:,:,:) = interp_s2f(:,:,:)
  end if
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  do nc = 1,ncell
    !
    n1 = pts_in_cell_ptr(nc)+1
    n2 = pts_in_cell_ptr(nc+1)
    np = n2-n1+1
    nl = merge(1,2,np == n_sp2dt)
    !
    f1 = flxpts_of_cell_ptr(nc)+1
    f2 = flxpts_of_cell_ptr(nc+1)
    !
    ! Create transpose of usp array to reduce cache misses
    !
    do k = 1,np
      do m = 1,neqns
        tusp(k,m) = usp(m,k+n1-1)
      end do
    end do
    !
    ! Now Interpolate the interior solution gradients to the flux points
    !
    !
    do n = f1,f2
      !
      nf = flxpts_of_cell(1,n)
      kf = flxpts_of_cell(2,n)
      lf = flxpts_of_cell(3,n) ! ideally this should evaluate to 1:f2-f1+1 in this loop
      sf = flxpts_of_cell(4,n)
      !
      !######################################################################
      !######################################################################
      ! IF USING LOBATTO SOLUTION POINTS WE SHOULD CHANGE lf TO GIVE THE
      ! INTERIOR SOLUTION POINT THAT CORRESPONDS TO THE CURRENT FLUX POINT
      ! TO PREVENT HAVING TO DO THE ENTIRE DOT_PRODUCT COMPUTATION
      !######################################################################
      !######################################################################
      !
     !write (mypnum+110,'(*(i6))') nc,n,f1,f2,nf,kf,lf,sf
      !
      lmat(1:np) = interp_mat(1:np,lf,nl)
      !
      do m = 1,neqns
        faceusp(m,kf,sf,nf) = dot_product(lmat(1:np),tusp(1:np,m))
      end do
      !
    end do
    !
  end do
  !
  ! Exchange the left and right states for faces on communication boundaries
  !
  if (ncpu > 1) call exchange_faceusp(faceusp)
  !
end subroutine interp_cv_to_flxpts
!
!###############################################################################
!
subroutine interp_dusp_to_flxpts()
  !
  use geovar,  only : nr,neqns,ncell,n_sp2dt,n_sp2dq,pts_in_cell_ptr
  use geovar,  only : flxpts_of_cell_ptr,flxpts_of_cell
  use geovar,  only : interp_s2f,interp_s2g,n_sp1de,n_fp1de
  use flowvar, only : dusp,facedusp
  !
  use parallel_mod, only : exchange_facedusp
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,nl,k,l,m
  integer :: nf,f1,f2,kf,lf,sf
  !
  !.. Local Arrays ..
  real(wp), dimension(1:n_sp2dq) :: lmat
  real(wp), dimension(1:n_sp2dq,1:nr,1:neqns) :: tdusp
  real(wp), dimension(1:n_sp2dq,1:4*n_fp1de,1:2) :: interp_mat
  !
continue
  !
  !######################################################################
  !######################################################################
  ! WE NEED TO DO SOMETHING FOR THE COMMUNICATION AND PERIODIC BOUNDARY
  ! CONDITIONS TO INTERPOLATE THE GHOST CELL SOLUTION POINTS TO THE FLUX
  ! POINTS ON THE BOUNDARY FACE.
  !
  ! AFTER CALLING THIS ROUTINE, WE NEED TO EXCHANGE THE FACE DATA
  ! FOR COMMUNICATION AND PERIODIC BOUNDARY CONDITIONS USING MPI
  ! NON-BLOCKING ROUTINES. THIS WAY WE CAN INITIALIZE THIS COMMUNICATION
  ! AND THEN CONTINUE ON TO COMPUTING THE DISCONTINUOUS FLUX DERIVATIVES
  ! AT THE INTERNAL SOLUTION POINTS WHILE THE COMMUNICATION FINISHES.
  ! THIS CAN BE DONE BECAUSE WE DONT NEED THE INTERFACE INFORMATION UNTIL
  ! AFTER THAT STEP.
  !######################################################################
  !######################################################################
  !
  if (n_fp1de /= n_sp1de) then
    interp_mat(:,:,:) = interp_s2g(:,:,:)
  else
    interp_mat(:,:,:) = interp_s2f(:,:,:)
  end if
  !
  ! Loop through each interior cell and interpolate the conserved
  ! variables at the solution points to the flux points of each cell
  !
  do nc = 1,ncell
    !
    n1 = pts_in_cell_ptr(nc)+1
    n2 = pts_in_cell_ptr(nc+1)
    np = n2-n1+1
    nl = merge(1,2,np == n_sp2dt)
    !
    f1 = flxpts_of_cell_ptr(nc)+1
    f2 = flxpts_of_cell_ptr(nc+1)
    !
    ! Create transpose of dusp array to reduce cache misses
    !
    do k = 1,np
      do m = 1,neqns
        do l = 1,nr
          tdusp(k,l,m) = dusp(l,m,k+n1-1)
        end do
      end do
    end do
    !
    ! Now Interpolate the interior solution gradients to the flux points
    !
    do n = f1,f2
      !
      nf = flxpts_of_cell(1,n)
      kf = flxpts_of_cell(2,n)
      lf = flxpts_of_cell(3,n) ! ideally this should evaluate to 1:f2-f1+1 in this loop
      sf = flxpts_of_cell(4,n)
      !
      !######################################################################
      !######################################################################
      ! IF USING LOBATTO SOLUTION POINTS WE SHOULD CHANGE lf TO GIVE THE
      ! INTERIOR SOLUTION POINT THAT CORRESPONDS TO THE CURRENT FLUX POINT
      ! TO PREVENT HAVING TO DO THE ENTIRE DOT_PRODUCT COMPUTATION
      !######################################################################
      !######################################################################
      !
     !write (mypnum+110,'(*(i6))') nc,n,f1,f2,nf,kf,lf,sf
      !
      lmat(1:np) = interp_mat(1:np,lf,nl)
      !
      do m = 1,neqns
        do l = 1,nr
          facedusp(l,m,kf,sf,nf) = dot_product(lmat(1:np),tdusp(1:np,l,m))
        end do
      end do
      !
    end do
    !
  end do
  !
  ! Exchange the left and right states for faces on communication boundaries
  !
  if (ncpu > 1) call exchange_facedusp(facedusp)
  !
end subroutine interp_dusp_to_flxpts
!
!###############################################################################
!###############################################################################
!###############################################################################
!
! Procedures for using quadrature to compute
! the divergence of only inviscid fluxes
!
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine derflux1_qdtr
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : neqns,ncell,n_sp1de,n_fp1de,n_sp2dt,n_sp2dq,nr
  use geovar,  only : pts_in_cell_ptr,interp_s2q,n_qp2dq,n_qp2dt
  use geovar,  only : flxpts_of_cell,flxpts_of_cell_ptr,projct_q2s
  use geovar,  only : interp_s2f,interp_s2g,using_quadrature
  use metvar,  only : dx,dy,dr_mat,ds_mat,geofa
  use metvar,  only : dqr_LPmat,dqs_LPmat,jacobian
  use flowvar, only : usp,residual,faceflx,faceusp
  use ovar,    only : dflux_opt
  !
  !.. Local Scalars ..
  integer :: k,l,m,n,nw,n1,n2,nsp,nqp
  integer :: f1,f2,nf,sf,kf,lf
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)      :: unrm
  real(wp), dimension(1:neqns)   :: ff,gg
  real(wp), dimension(1:n_sp2dq) :: lmat
  real(wp), dimension(1:n_qp2dq) :: dxdr,dxds
  real(wp), dimension(1:n_qp2dq) :: dydr,dyds
  real(wp), dimension(1:n_qp2dq) :: qjac
  !
  real(wp), dimension(1:n_sp2dq,1:neqns) :: tusp
  real(wp), dimension(1:neqns,1:n_qp2dq) :: uqp
  !
  real(wp), dimension(1:n_sp2dq,1:neqns,1:nr) :: rsflx_sp
  real(wp), dimension(1:n_sp2dq,1:neqns,1:nr) :: xyflx_sp
  real(wp), dimension(1:n_sp2dq,1:neqns,1:nr) :: dflx_sp
  real(wp), dimension(1:n_qp2dq,1:neqns,1:nr) :: rsflx_qp
  real(wp), dimension(1:n_qp2dq,1:neqns,1:nr) :: xyflx_qp
  !
continue
  !
  ! Initialize faceflx to zero before we start
  !
  faceflx(:,:,:,:) = zero
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    nsp = n2-n1+1 ! number of solution points in this cell
    nqp = merge(n_qp2dt,n_qp2dq,nsp==n_sp2dt) ! no. of QP in this cell
    !
    nw = merge( 1 , 2 , nsp == n_sp2dt ) ! check if cell n is tria or quad
    !
    f1 = flxpts_of_cell_ptr(n)+1
    f2 = flxpts_of_cell_ptr(n+1)
    !
    ! Store the transpose of usp to reduce cache misses
    !
    tusp(1:nsp,1:neqns) = transpose(usp(1:neqns,n1:n2))
    !
    ! Get the information local to this cell
    !
    if (using_quadrature) then
      !
      do k = 1,nqp
        !
        ! Get the conservative solution at the quadrature points
        !
        do m = 1,neqns
          uqp(m,k) = dot_product( interp_s2q(1:nsp,k,nw) , tusp(1:nsp,m) )
        end do
        !
        ! Get the metrics for this cell
        !
        dxdr(k) = dot_product( interp_s2q(1:nsp,k,nw) , dx(n1:n2)%dr )
        dxds(k) = dot_product( interp_s2q(1:nsp,k,nw) , dx(n1:n2)%ds )
        dydr(k) = dot_product( interp_s2q(1:nsp,k,nw) , dy(n1:n2)%dr )
        dyds(k) = dot_product( interp_s2q(1:nsp,k,nw) , dy(n1:n2)%ds )
        qjac(k) = dot_product( interp_s2q(1:nsp,k,nw) , jacobian(n1:n2) )
        !
      end do
      !
    else
      !
      uqp(1:neqns,1:nsp) = usp(1:neqns,n1:n2)
      !
      dxdr(1:nsp) = dx(n1:n2)%dr
      dxds(1:nsp) = dx(n1:n2)%ds
      dydr(1:nsp) = dy(n1:n2)%dr
      dyds(1:nsp) = dy(n1:n2)%ds
      qjac(1:nsp) = jacobian(n1:n2)
      !
    end if
    !
    rsflx_qp(:,:,:) = zero ! zero for dot product later with hybrid tria/quad
    xyflx_qp(:,:,:) = zero ! zero for dot product later with hybrid tria/quad
    !
    ! Get the flux at each quadrature point in the current cell
    !
    do k = 1,nqp
      !
      xyflx_qp(k,:,1) = f_invs_flx_cv( uqp(:,k) ) ! x-directional flux
      xyflx_qp(k,:,2) = g_invs_flx_cv( uqp(:,k) ) ! y-directional flux
      !
      rsflx_qp(k,:,1) =  xyflx_qp(k,:,1)*dyds(k) - xyflx_qp(k,:,2)*dxds(k)
      rsflx_qp(k,:,2) = -xyflx_qp(k,:,1)*dydr(k) + xyflx_qp(k,:,2)*dxdr(k)
      !
    end do
    !
    ! Two different methods for getting the derivative at the solution points
    !
    if (.not.using_quadrature .or. dflux_opt == derivative_of_proj_fluxes) then
      !
      ! Project the fluxes back to the solution polynomial order
      !
      if (using_quadrature) then
        !
        do l = 1,nr
          do m = 1,neqns
            do k = 1,nsp
              xyflx_sp(k,m,l) = dot_product( projct_q2s(1:nqp,k,nw) , xyflx_qp(1:nqp,m,l) )
              rsflx_sp(k,m,l) = dot_product( projct_q2s(1:nqp,k,nw) , rsflx_qp(1:nqp,m,l) )
            end do
          end do
        end do
        !
      else
        !
        xyflx_sp(1:nsp,1:neqns,1:nr) = xyflx_qp(1:nsp,1:neqns,1:nr)
        rsflx_sp(1:nsp,1:neqns,1:nr) = rsflx_qp(1:nsp,1:neqns,1:nr)
        !
      end if
      !
      ! Take the derivative of the projected fluxes
      !
      do m = 1,neqns
        do k = 1,nsp
          dflx_sp(k,m,1) = dot_product( dr_mat(1:nsp,k,nw) , rsflx_sp(1:nsp,m,1) )
          dflx_sp(k,m,2) = dot_product( ds_mat(1:nsp,k,nw) , rsflx_sp(1:nsp,m,2) )
        end do
      end do
      !
      ! Interpolate the fluxes to the interface flux points and dot this
      ! interpolated vector with the face normal at each flux point
      ! NOTE: We dont need to exchange this information for communication and
      !       periodic boundary faces since this flux is only used within the
      !       current cell.
      !
      do k = f1,f2
        !
        nf = flxpts_of_cell(1,k)
        kf = flxpts_of_cell(2,k)
        lf = flxpts_of_cell(3,k)
        sf = flxpts_of_cell(4,k)
        !
        ! Get the unit normal for this face flux/solution point
        !
        unrm(:) = unit_vector( geofa(:,kf,nf) )
        !
        ! WE NEED TO ADD SOMETHING HERE IN CASE LOBATTO NODES ARE USED. LOOK AT 
        ! ROUTINE interp_cv_to_flxpts FOR MORE DETAILS.
        !
        lmat(1:nsp) = interp_s2f(1:nsp,lf,nw)
        !
        do l = 1,nr
          do m = 1,neqns
            faceflx(m,kf,sf,nf) = faceflx(m,kf,sf,nf) + &
                                  dot_product(lmat(1:nsp),xyflx_sp(1:nsp,m,l))*unrm(l)
          end do
        end do
        !
      end do
      !
    else if (dflux_opt == project_derivatives) then
      !
      ! Compute the derivative of the fluxes at the quadrature points
      ! and project these derivatives to the solution points
      !
      do k = 1,nsp
        do m = 1,neqns
          dflx_sp(k,m,1) = dot_product( dqr_LPmat(1:nqp,k,nw) , rsflx_qp(1:nqp,m,1) )
          dflx_sp(k,m,2) = dot_product( dqs_LPmat(1:nqp,k,nw) , rsflx_qp(1:nqp,m,2) )
        end do
      end do
      !
      ! Interpolate the fluxes to the interface flux points and dot this
      ! interpolated vector with the face normal at each flux point
      ! NOTE: We dont need to exchange this information for communication and
      !       periodic boundary faces since this flux is only used within the
      !       current cell.
      !
     !do k = f1,f2
     !  !
     !  nf = flxpts_of_cell(1,k)
     !  kf = flxpts_of_cell(2,k)
     !  lf = flxpts_of_cell(3,k)
     !  sf = flxpts_of_cell(4,k)
     !  !
     !  ! #################################################################
     !  ! #################################################################
     !  ! #################################################################
     !  ! ##  ERROR! THIS WILL NOT WORK BECAUSE WE NEED TO INTERPOLATE  ##
     !  ! ##         THE FACE NORMALS TO THE FACE QUADRATURE POINTS!    ##
     !  ! #################################################################
     !  ! #################################################################
     !  ! #################################################################
     !  !
     !  ! Get the unit normal for this face flux/solution point
     !  !
     !  unrm(:) = unit_vector( geofa(:,kf,nf) )
     !  !
     !  ! WE NEED TO ADD SOMETHING HERE IN CASE LOBATTO NODES ARE USED. LOOK AT 
     !  ! ROUTINE interp_cv_to_flxpts FOR MORE DETAILS.
     !  !
     !  lmat(1:nsp) = interp_s2g(1:nsp,lf,nw)
     !  !
     !  do l = 1,nr
     !    do m = 1,neqns
     !      faceflx(m,kf,sf,nf) = faceflx(m,kf,sf,nf) + &
     !                            dot_product(lmat(1:nsp),xyflx_qp(1:nsp,m,l))*unrm(l)
     !    end do
     !  end do
     !  !
     !end do
      !
    end if
    !
    ! Add the derivatives of the flux to the residual
    !
    do m = 1,neqns
      do k = 1,nsp
        residual(m,k+n1-1) = dflx_sp(k,m,1) + dflx_sp(k,m,2)
      end do
    end do
    !
   !do k = 1,nsp
   !  do m = 1,neqns
   !   !residual(m,k+n1-1) = dot_product(dr_mat(1:nsp,k,nw),ftdc(1:nsp,m)) - &
   !   !                     dot_product(ds_mat(1:nsp,k,nw),gtdc(1:nsp,m))
   !   !residual(m,k+n1-1) = dot_product(dqr_mat(:,k,nw),ftdc(:,m)) - &
   !   !                     dot_product(dqs_mat(:,k,nw),gtdc(:,m))
   !    residual(m,k+n1-1) = dot_product(dr_mat(1:nsp,k,nw),dfsp(1:nsp,m)) - &
   !                         dot_product(ds_mat(1:nsp,k,nw),dgsp(1:nsp,m))
   !  end do
   !end do
    !
   !do m = 1,neqns
   !  residual(m,n1:n2) = matmul(invM_mat(:,:,nw),residual(m,n1:n2))!/jacobian(n1:n2)
   !end do
    !
  end do
  !
end subroutine derflux1_qdtr
!
!###############################################################################
!
subroutine correct_interface_flux_qdtr()
  !
  use geovar,  only : neqns,nfbnd,nface,n_sp1de,n_fp1de,n_sp2dt,n_sp2dq,nr
  use geovar,  only : bface,cells_on_face,pts_in_cell_ptr,xyz,flxwts_1d
  use geovar,  only : interp_s2f,n_qp1de,interp_f2g,projct_g2f,using_quadrature
  use ovar,    only : flux_divergence,dflux_opt,itcur,nst
  use metvar,  only : cormat,geofa,volj
  use flowvar, only : usp,residual,faceflx,faceusp
  !
  use generic_mod,  only : pressure_cv_sp
  use parallel_mod, only : wait_faceusp
  !
  !.. Local Scalars ..
  integer   :: k,l,m,n,nf,r1,r2,l1,l2
  integer   :: jl,jr,wl,wr,nsp,nqp
  integer   :: kbl,kbr,kl,kr,iu
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)    :: unrm
  real(wp), dimension(1:neqns) :: ul,ur
  real(wp), dimension(1:neqns) :: dfl,dfr
  real(wp), dimension(1:neqns) :: dgl,dgr
  real(wp), dimension(1:neqns) :: fnl,fnr,fnc
  real(wp), dimension(1:neqns) :: srcl,srcr
  real(wp), dimension(1:nr,1:n_qp1de)    :: rnrm
  real(wp), dimension(1:neqns,1:n_sp1de) :: fjl,fjr
  real(wp), dimension(1:neqns,1:n_fp1de) :: lflx,rflx
  real(wp), dimension(1:neqns,1:n_qp1de) :: luqp,ruqp,fncq
  !
continue
  !
  nsp = n_sp1de
  nqp = n_qp1de
  !
  ! Make sure that the non-blocking exchange of faceusp has completed
  !
  if (ncpu > 1) call wait_faceusp(faceusp)
  !
  ! Loop over all the faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  if (output_boundary_fluxes) write (102,110) itcur,nst
  !
  do nf = 1,nface
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kbr = cells_on_face(6,nf) ! ending    flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    wl = merge( 1 , 2 , l2-l1+1 == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , r2-r1+1 == n_sp2dt ) ! check if jr is a tria or quad
    !
    ! Interpolate the conserved variables and face normals at
    ! the solution flux points to the quadrature flux points
    !
    if (using_quadrature) then
      !
      do k = 1,nqp
        !
        do m = 1,neqns
          luqp(m,k) = dot_product( interp_f2g(1:nsp,k) , faceusp(m,1:nsp,1,nf) )
          ruqp(m,k) = dot_product( interp_f2g(1:nsp,k) , faceusp(m,1:nsp,2,nf) )
        end do
        !
        do l = 1,nr
          rnrm(l,k) = dot_product( interp_f2g(1:nsp,k) , geofa(l,1:nsp,nf) )
        end do
        !
      end do
      !
    else
      !
      lflx(1:neqns,1:nsp) = faceflx(1:neqns,1:nsp,1,nf)
      rflx(1:neqns,1:nsp) = faceflx(1:neqns,1:nsp,2,nf)
      luqp(1:neqns,1:nsp) = faceusp(1:neqns,1:nsp,1,nf)
      ruqp(1:neqns,1:nsp) = faceusp(1:neqns,1:nsp,2,nf)
      rnrm(1:nr,1:nsp) = geofa(1:nr,1:nsp,nf)
      !
    end if
    !
    if (output_boundary_fluxes) then
      if (nf <= nfbnd) then
        write (102,111) nf,trim(bc_string(bface(1,nf)))
        do k = 1,nsp
          write (102,112) k,luqp(1:neqns,k),pressure_cv_sp(luqp(1:neqns,k))
          write (102,113) k,ruqp(1:neqns,k),pressure_cv_sp(ruqp(1:neqns,k))
        end do
      end if
    end if
    !
    ! Get the discontinuous fluxes at the interface
    !
    if (dflux_opt == project_derivatives) then
      !
      ! Get the discontinuous fluxes at the quadrature flux points
      !
      do k = 1,nqp
        lflx(1:neqns,k) = normal_invs_flx_cv( luqp(:,k) , rnrm(:,k) )
        rflx(1:neqns,k) = normal_invs_flx_cv( ruqp(:,k) , rnrm(:,k) )
      end do
      !
    else
      !
      ! Get the discontinuous fluxes at the solution flux points
      !
      lflx(1:neqns,1:nsp) = faceflx(1:neqns,1:nsp,1,nf)
      rflx(1:neqns,1:nsp) = faceflx(1:neqns,1:nsp,2,nf)
      !
    end if
    !
    ! Get the common flux at each quadrature flux point
    !
    do k = 1,nqp
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( rnrm(:,k) )
      !
      ! Get the numerical flux at the current flux point
      !
      if (flux_method == ROEFLUX) then
        fncq(:,k) = upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      else if (flux_method == HLLC) then
        fncq(:,k) = hllc_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      else if (flux_method == LDFSS) then
        fncq(:,k) = ldfss_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      end if
      !
    end do
    !
    ! Project the common flux to the solution flux points
    !
    if (using_quadrature) then
      !
      if (dflux_opt == project_derivatives) then
        !
        lflx(1:neqns,1:nqp) = fncq(1:neqns,1:nqp) - lflx(1:neqns,1:nqp)
        rflx(1:neqns,1:nqp) = fncq(1:neqns,1:nqp) - rflx(1:neqns,1:nqp)
        !
        do k = 1,nsp
          do m = 1,neqns
            fjl(m,k) = dot_product( projct_g2f(1:nqp,k) , lflx(m,1:nqp) )
            fjr(m,k) = dot_product( projct_g2f(1:nqp,k) , rflx(m,1:nqp) )
          end do
        end do
        !
      else
        !
        do k = 1,nsp
          !
          do m = 1,neqns
            fnc(m) = dot_product( projct_g2f(1:nqp,k) , fncq(m,1:nqp) )
          end do
          !
          fjl(:,k) = fnc(:) - lflx(:,k)
          fjr(:,k) = fnc(:) - rflx(:,k)
          !
        end do
        !
      end if
      !
    else
      !
      fjl(1:neqns,1:nsp) = fncq(1:neqns,1:nsp) - lflx(1:neqns,1:nsp)
      fjr(1:neqns,1:nsp) = fncq(1:neqns,1:nsp) - rflx(1:neqns,1:nsp)
      ! 
    end if
    !
    ! Compute the flux jump at each solution flux point
    !
   !do k = 1,nsp
   !  fjl(:,k) = fnc(:,k) - lflx(:,k)
   !  fjr(:,k) = fnc(:,k) - rflx(:,k)
   !end do
   !!
   !! Initialize the conservation source terms
   !!
   !srcl(:) = zero
   !srcr(:) = zero
   !!
   !1 format ("     lflx        rflx        fjl         fjr         fnc",/,*(i8))
   !4 format ("     luqp        ruqp",/,*(i8))
   !3 format (*(i6))
   !2 format (*(es12.4))
   !iu = 103 
   !write (iu,1) nf,jl,jr
   !write (iu+1,4) nf,jl,jr
    do k = 1,nsp
      !
      ! Get the index for the correction function 
      ! corresponding to the current interface flux point
      !
      kl = kbl-1 + k ! index for left  cell
      kr = kbr+1 - k ! index for right cell
      !
     !write (iu,3) k,kl,kr
     !write (iu+1,3) k,kl,kr
     !do m = 1,neqns
     !  write (iu,2) lflx(m,k),rflx(m,k),fjl(m,k),fjr(m,k),fnc(m,k)
     !  write (iu+1,2) luqp(m,k),ruqp(m,k)
     !end do
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( rnrm(:,k) )
      dsl = ds * cl(wl)
      dsr = ds * cl(wr)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
      do n = l1,l2
        do m = 1,neqns
          residual(m,n) = residual(m,n) + dsl*fjl(m,k)*cormat(n-l1+1,kl,wl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jr
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      if (nf > nfbnd) then
        do n = r1,r2
          do m = 1,neqns
            residual(m,n) = residual(m,n) - dsr*fjr(m,k)*cormat(n-r1+1,kr,wr)
          end do
        end do
      end if
      !
      ! From ZJ, to ensure conservation with chain-rule approach
      !
     !if (flux_divergence == ChainRuleCV) then
     !  do m = 1,neqns   ! ARE DSL AND DSR CORRECT HERE?
     !    srcl(m) = srcl(m) + dsl*lflx(m,k)*flxwts_1d(k)
     !    srcr(m) = srcr(m) + dsr*rflx(m,k)*flxwts_1d(k)
     !  end do
     !end if
      !
    end do
    !
    ! Now add this faces contribution to the surface conservation
    ! source term to each solution point in the left and right cells
    !
   !if (flux_divergence == ChainRuleCV) then
   !  !
   !  do n = l1,l2
   !    do m = 1,neqns
   !      residual(m,n) = residual(m,n) + srcl(m)/volj(n)
   !    end do
   !  end do
   !  !
   !  if (nf > nfbnd) then
   !    do n = r1,r2
   !      do m = 1,neqns
   !        residual(m,n) = residual(m,n) + srcr(m)/volj(n)
   !      end do
   !    end do
   !  end if
   !  !
   !end if
    !
  end do
  !
  110 format (//"itcur = ",i0,"; rkstep = ",i0)
  111 format ("Boundary face ",i0," with BC = ",a)
  112 format ("FP-",i0,"; Left : ",*(es14.6))
  113 format ("FP-",i0,"; Right: ",*(es14.6))
  !
end subroutine correct_interface_flux_qdtr
!
!###############################################################################
!###############################################################################
!###############################################################################
!
! Procedures for computing the divergence of inviscid and viscous fluxes
!
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine derflux1_qdtr_with_visc
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : neqns,ncell,n_sp1de,n_fp1de,n_sp2dt,n_sp2dq,nr
  use geovar,  only : pts_in_cell_ptr,interp_s2q,n_qp2dq,n_qp2dt
  use geovar,  only : flxpts_of_cell,flxpts_of_cell_ptr,projct_q2s
  use geovar,  only : interp_s2f,interp_s2g,using_quadrature
  use metvar,  only : dx,dy,dr_mat,ds_mat,geofa
  use metvar,  only : dqr_LPmat,dqs_LPmat,jacobian
  use flowvar, only : usp,dusp,residual,faceflx,faceusp
  use ovar,    only : dflux_opt
  !
  !.. Local Scalars ..
  integer :: k,l,m,n,nw,n1,n2,np,nq
  integer :: f1,f2,nf,sf,kf,lf
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)      :: unrm
  real(wp), dimension(1:neqns,1:nr)   :: vflx
  real(wp), dimension(1:n_sp2dq) :: lmat
  !
  real(wp), dimension(1:n_sp2dq,1:neqns,1:nr) :: rsflx
  real(wp), dimension(1:n_sp2dq,1:neqns,1:nr) :: xyflx
  !
continue
  !
  nq = neqns
  !
  ! Initialize faceflx to zero before we start
  !
  faceflx(:,:,:,:) = zero
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    np = n2-n1+1 ! number of solution points in this cell
    !
    nw = merge( 1 , 2 , np == n_sp2dt ) ! check if cell n is tria or quad
    !
    f1 = flxpts_of_cell_ptr(n)+1
    f2 = flxpts_of_cell_ptr(n+1)
    !
    if (.not. first_output_done .and. output_boundary_fluxes) then
      write (120,110) n
      do k = n1,n2
        write (120,111) k,dx(k)%dr,dx(k)%ds,dy(k)%dr,dy(k)%ds,jacobian(k) 
      end do
    end if
    110 format ("Cell ",i0)
    111 format ("SP = ",i0,"; dxdr=",es14.6,", dxds=",es14.6, &
                ", dydr=",es14.6,", dyds=",es14.6,", jacobian=",es14.6)
    !
    ! Get the flux at each quadrature point in the current cell
    !
    do k = 1,np
      !
      l = k + n1-1
      !
      ! Inviscid fluxes
      !
      xyflx(k,:,1) = f_invs_flx_cv( usp(:,l) ) ! x-directional flux
      xyflx(k,:,2) = g_invs_flx_cv( usp(:,l) ) ! y-directional flux
      !
      ! Viscous fluxes
      !
      vflx(1:nq,1:nr) = visc_flx_cv( usp(1:nq,l) , dusp(1:nr,1:nq,l) )
      !
      xyflx(k,1:nq,1:nr) = xyflx(k,1:nq,1:nr) - vflx(1:nq,1:nr)
      !
      ! Transform the fluxes in the x-y physical
      ! domain to the r-s computational domain
      !
      rsflx(k,:,1) =  xyflx(k,:,1)*dy(l)%ds - xyflx(k,:,2)*dx(l)%ds
      rsflx(k,:,2) = -xyflx(k,:,1)*dy(l)%dr + xyflx(k,:,2)*dx(l)%dr
      !
    end do
    !
    ! Compute the divergence of the discontinuous fluxes 
    !
    do m = 1,nq
      do k = 1,np
        residual(m,k+n1-1) = dot_product( dr_mat(1:np,k,nw) , rsflx(1:np,m,1) ) + &
                             dot_product( ds_mat(1:np,k,nw) , rsflx(1:np,m,2) )
      end do
    end do
    !
    ! Interpolate the fluxes to the interface flux points and dot this
    ! interpolated vector with the face normal at each flux point
    ! NOTE: We dont need to exchange this information for communication and
    !       periodic boundary faces since this flux is only used within the
    !       current cell.
    !
    do k = f1,f2
      !
      nf = flxpts_of_cell(1,k)
      kf = flxpts_of_cell(2,k)
      lf = flxpts_of_cell(3,k)
      sf = flxpts_of_cell(4,k)
      !
      ! Get the unit normal for this face flux/solution point
      !
      unrm(:) = unit_vector( geofa(:,kf,nf) )
      !
      ! WE NEED TO ADD SOMETHING HERE IN CASE LOBATTO NODES ARE USED. LOOK AT 
      ! ROUTINE interp_cv_to_flxpts FOR MORE DETAILS.
      !
      lmat(1:np) = interp_s2f(1:np,lf,nw)
      !
      do l = 1,nr
        do m = 1,nq
          faceflx(m,kf,sf,nf) = faceflx(m,kf,sf,nf) + &
                                dot_product(lmat(1:np),xyflx(1:np,m,l))*unrm(l)
        end do
      end do
      !
    end do
    !
   !do k = 1,np
   !  do m = 1,nq
   !   !residual(m,k+n1-1) = dot_product(dr_mat(1:np,k,nw),ftdc(1:np,m)) - &
   !   !                     dot_product(ds_mat(1:np,k,nw),gtdc(1:np,m))
   !   !residual(m,k+n1-1) = dot_product(dqr_mat(:,k,nw),ftdc(:,m)) - &
   !   !                     dot_product(dqs_mat(:,k,nw),gtdc(:,m))
   !    residual(m,k+n1-1) = dot_product(dr_mat(1:np,k,nw),dfsp(1:np,m)) - &
   !                         dot_product(ds_mat(1:np,k,nw),dgsp(1:np,m))
   !  end do
   !end do
    !
   !do m = 1,nq
   !  residual(m,n1:n2) = matmul(invM_mat(:,:,nw),residual(m,n1:n2))!/jacobian(n1:n2)
   !end do
    !
  end do
  !
  first_output_done = true
  !
end subroutine derflux1_qdtr_with_visc
!
!###############################################################################
!
subroutine correct_interface_flux_qdtr_with_visc()
  !
  use geovar,  only : neqns,nfbnd,nface,n_sp1de,n_fp1de,n_sp2dt,n_sp2dq,nr
  use geovar,  only : bface,cells_on_face,pts_in_cell_ptr,xyz,flxwts_1d
  use geovar,  only : interp_s2f,n_qp1de,interp_f2g,projct_g2f,using_quadrature
  use ovar,    only : flux_divergence,dflux_opt,itcur,nst
  use metvar,  only : cormat,geofa,volj
  use flowvar, only : usp,residual,faceflx,faceusp,facedusp
  !
  use generic_mod, only : pressure_cv_sp
  use parallel_mod, only : wait_facedusp
  !
  !.. Local Scalars ..
  integer   :: k,l,m,n,nf,r1,r2,l1,l2
  integer   :: jl,jr,wl,wr,nsp
  integer   :: kbl,kbr,kl,kr,iu
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)    :: unrm
  real(wp), dimension(1:neqns) :: ul,ur
  real(wp), dimension(1:neqns) :: dfl,dfr
  real(wp), dimension(1:neqns) :: dgl,dgr
  real(wp), dimension(1:neqns) :: fnl,fnr,fnc
  real(wp), dimension(1:neqns) :: srcl,srcr
  real(wp), dimension(1:neqns) :: vflx
  real(wp), dimension(1:nr,1:n_sp1de)    :: rnrm
  real(wp), dimension(1:neqns,1:n_sp1de) :: fjl,fjr
  real(wp), dimension(1:neqns,1:n_sp1de) :: lflx,rflx
  real(wp), dimension(1:neqns,1:n_sp1de) :: luqp,ruqp,fncq
  real(wp), dimension(1:nr,1:neqns,1:n_sp1de) :: lduqp,rduqp
  !
continue
  !
  nsp = n_sp1de
  !
  ! Make sure that the non-blocking exchange of faceusp has completed
  !
  if (ncpu > 1) call wait_facedusp(facedusp)
  !
  ! Loop over all the faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  if (output_boundary_fluxes) write (101,110) itcur,nst
  !
  do nf = 1,nface
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kbr = cells_on_face(6,nf) ! ending    flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    wl = merge( 1 , 2 , l2-l1+1 == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , r2-r1+1 == n_sp2dt ) ! check if jr is a tria or quad
    !
    lflx(1:neqns,1:nsp) = faceflx(1:neqns,1:nsp,1,nf)
    rflx(1:neqns,1:nsp) = faceflx(1:neqns,1:nsp,2,nf)
    !
    luqp(1:neqns,1:nsp) = faceusp(1:neqns,1:nsp,1,nf)
    ruqp(1:neqns,1:nsp) = faceusp(1:neqns,1:nsp,2,nf)
    !
    if (output_boundary_fluxes) then
      if (nf <= nfbnd) then
        write (101,111) nf,trim(bc_string(bface(1,nf)))
        do k = 1,nsp
          write (101,112) k,luqp(1:neqns,k),pressure_cv_sp(luqp(1:neqns,k))
          write (101,113) k,ruqp(1:neqns,k),pressure_cv_sp(ruqp(1:neqns,k))
        end do
      end if
    end if
    !
    lduqp(1:nr,1:neqns,1:nsp) = facedusp(1:nr,1:neqns,1:nsp,1,nf)
    if (nf <= nfbnd) then
      rduqp(1:nr,1:neqns,1:nsp) = lduqp(1:nr,1:neqns,1:nsp)
    else
      rduqp(1:nr,1:neqns,1:nsp) = facedusp(1:nr,1:neqns,1:nsp,2,nf)
    end if
    !
    rnrm(1:nr,1:nsp) = geofa(1:nr,1:nsp,nf)
    !
    ! Get the common flux at each quadrature flux point
    !
    do k = 1,nsp
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( rnrm(:,k) )
      !
      ! Get the inviscid upwind flux at the current flux point
      !
      if (flux_method == ROEFLUX) then
        fncq(:,k) = upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      else if (flux_method == HLLC) then
        fncq(:,k) = hllc_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      else if (flux_method == LDFSS) then
        fncq(:,k) = ldfss_upwind_flux_cv ( unrm(:) , luqp(:,k) , ruqp(:,k) )
      end if
      !
      ! Get the common viscous flux at the current flux point
      !
      vflx(:) = common_visc_flx_cv( unrm(:), luqp(:,k), lduqp(:,:,k), &
                                             ruqp(:,k), rduqp(:,:,k) )
      fncq(:,k) = fncq(:,k) - vflx(:)
      !
    end do
    !
    fjl(1:neqns,1:nsp) = fncq(1:neqns,1:nsp) - lflx(1:neqns,1:nsp)
    fjr(1:neqns,1:nsp) = fncq(1:neqns,1:nsp) - rflx(1:neqns,1:nsp)
    !
   !!
   !! Initialize the conservation source terms
   !!
   !srcl(:) = zero
   !srcr(:) = zero
   !!
    do k = 1,nsp
      !
      ! Get the index for the correction function 
      ! corresponding to the current interface flux point
      !
      kl = kbl-1 + k ! index for left  cell
      kr = kbr+1 - k ! index for right cell
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( rnrm(:,k) )
      dsl = ds * cl(wl)
      dsr = ds * cl(wr)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
      do n = l1,l2
        do m = 1,neqns
          residual(m,n) = residual(m,n) + dsl*fjl(m,k)*cormat(n-l1+1,kl,wl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jr
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      if (nf > nfbnd) then
        do n = r1,r2
          do m = 1,neqns
            residual(m,n) = residual(m,n) - dsr*fjr(m,k)*cormat(n-r1+1,kr,wr)
          end do
        end do
      end if
      !
      ! From ZJ, to ensure conservation with chain-rule approach
      !
     !if (flux_divergence == ChainRuleCV) then
     !  do m = 1,neqns   ! ARE DSL AND DSR CORRECT HERE?
     !    srcl(m) = srcl(m) + dsl*lflx(m,k)*flxwts_1d(k)
     !    srcr(m) = srcr(m) + dsr*rflx(m,k)*flxwts_1d(k)
     !  end do
     !end if
      !
    end do
    !
    ! Now add this faces contribution to the surface conservation
    ! source term to each solution point in the left and right cells
    !
   !if (flux_divergence == ChainRuleCV) then
   !  !
   !  do n = l1,l2
   !    do m = 1,neqns
   !      residual(m,n) = residual(m,n) + srcl(m)/volj(n)
   !    end do
   !  end do
   !  !
   !  if (nf > nfbnd) then
   !    do n = r1,r2
   !      do m = 1,neqns
   !        residual(m,n) = residual(m,n) + srcr(m)/volj(n)
   !      end do
   !    end do
   !  end if
   !  !
   !end if
    !
  end do
  !
  110 format (//"itcur = ",i0,"; rkstep = ",i0)
  111 format ("Boundary face ",i0," with BC = ",a)
  112 format ("FP-",i0,"; Left : ",*(es14.6))
  113 format ("FP-",i0,"; Right: ",*(es14.6))
  !
end subroutine correct_interface_flux_qdtr_with_visc
!
!###############################################################################
!
subroutine correct_flux_project_cv()
  !
  use geovar,  only : neqns,nfbnd,nface,n_sp1de,n_sp2dt,n_sp2dq
  use geovar,  only : bface,cells_on_face,pts_in_cell_ptr,xyz,flxwts_1d
  use ovar,    only : flux_divergence
  use metvar,  only : cormat,geofa,volj
  use geovar,  only : interp_s2f
  use flowvar, only : usp,residual
  !
  !.. Local Scalars ..
  logical   :: is_1d_bc
  integer   :: k,l,m,n,nf,r1,r2,l1,l2
  integer   :: jl,jr,nl,nr,wl,wr
  integer   :: kcorbl,kcorbr,kcorl,kcorr
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:2)  :: unrm
  real(wp), dimension(1:neqns) :: ul,ur
  real(wp), dimension(1:neqns) :: dfl,dfr
  real(wp), dimension(1:neqns) :: dgl,dgr
  real(wp), dimension(1:neqns) :: fnl,fnr
  real(wp), dimension(1:neqns) :: fjl,fjr
  real(wp), dimension(1:neqns) :: fnc
  real(wp), dimension(1:neqns) :: srcl,srcr
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lusp
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rusp
  !
continue
  !
  ! Loop over all the faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  do nf = 1,nface
    !
    ! Check if this is a boundary face needing a 1D treatment of the flux points
    !
    is_1d_bc = fals
    if (nf <= nfbnd) is_1d_bc = .not. any(bface(1,nf) == bc_comm)
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kcorbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kcorbr = cells_on_face(6,nf) ! beginning flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    nl = l2-l1+1 ! number of solution points in left  cell
    nr = r2-r1+1 ! number of solution points in right cell
    !
    wl = merge( 1 , 2 , nl == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , nr == n_sp2dt ) ! check if jr is a tria or quad
    !
    ! Store the transpose of usp for the left
    ! and right cells to reduce cache misses
    !
    do k = l1,l2
      lusp (k-l1+1,:) = usp(:,k)
    end do
    !
    do k = r1,r2
      rusp(k-r1+1,:) = usp(:,k)
    end do
   !!
   !! Initialize the conservation source terms
   !!
   !srcl(:) = zero
   !srcr(:) = zero
    !
    ! Loop through the flux points on the current face
    !
    do k = 1,n_sp1de
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( geofa(:,k,nf) )
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( geofa(:,k,nf) )
      dsl = ds * cl(wl)
      dsr = ds * cl(wr)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kcorl = kcorbl-1 + k ! index for left  cell
      kcorr = kcorbr+1 - k ! index for right cell
      !
      ! Use the Lagrange polynomial to project the conservative variables
      ! from the left cells to the current flux point on the interface
      !
      do m = 1,neqns
        ul(m) = dot_product(interp_s2f(1:nl,kcorl,wl),usp(m,l1:l2)) ! left  state
      end do
      !
      ! Now get the information at flux point on the right side
      !
      if (is_1d_bc) then
        !
        ! Copy the specified boundary condition in usp into ur and use
        ! this value to get the discontinuous fluxes at this flux point
        !
        ur(:) = rusp(kcorr,:)
        !
      else
        !
        ! Use the Lagrange polynomial to project the conserved solution
        ! variables and the discontinuous fluxes from the right (non-ghost)
        ! cell to the current flux point on the interface.
        !
        do m = 1,neqns
          ur(m)  = dot_product(interp_s2f(1:nr,kcorr,wr), rusp(1:nr,m)) ! right solution
        end do
        !
      end if
      !
      ! Get the numerical flux at the current flux point
      !
      if (flux_method == ROEFLUX) then
        fnc(:) = upwind_flux_cv ( unrm(:) , ul , ur )
      else if (flux_method == HLLC) then
        fnc(:) = hllc_upwind_flux_cv ( unrm(:) , ul , ur )
      else if (flux_method == LDFSS) then
        fnc(:) = ldfss_upwind_flux_cv ( unrm(:) , ul , ur )
      end if
      !
      ! Get left and right discontinuous fluxes at the current flux point
      !
      dfl(:) = f_invs_flx_cv(ul)
      dgl(:) = g_invs_flx_cv(ul)
      dfr(:) = f_invs_flx_cv(ur)
      dgr(:) = g_invs_flx_cv(ur)
      !
      ! Get the discontinuous flux function at the current
      ! interface point in the direction of the face normal
      !
      fnl(:) = dfl(:)*unrm(1) + dgl(:)*unrm(2)
      fnr(:) = dfr(:)*unrm(1) + dgr(:)*unrm(2)
      !
      ! The jump in the flux at the current interface point
      !
      fjl(:) = fnc(:) - fnl(:)
      fjr(:) = fnc(:) - fnr(:)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
      do n = l1,l2
        do m = 1,neqns
          residual(m,n) = residual(m,n) + dsl*fjl(m)*cormat(n-l1+1,kcorl,wl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jr
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      if (nf > nfbnd) then
        do n = r1,r2
          do m = 1,neqns
            residual(m,n) = residual(m,n) - dsr*fjr(m)*cormat(n-r1+1,kcorr,wr)
          end do
        end do
      end if
      !
      ! From ZJ, to ensure conservation with chain-rule approach
      !
     !if (flux_divergence == ChainRuleCV) then
     !  do m = 1,neqns   ! ARE DSL AND DSR CORRECT HERE?
     !    srcl(m) = srcl(m) + dsl*fnl(m)*flxwts_1d(k)
     !    srcr(m) = srcr(m) + dsr*fnr(m)*flxwts_1d(k)
     !  end do
     !end if
      !
    end do
    !
    ! Now add this faces contribution to the surface conservation
    ! source term to each solution point in the left and right cells
    !
   !if (flux_divergence == ChainRuleCV) then
   !  !
   !  do n = l1,l2
   !    do m = 1,neqns
   !      residual(m,n) = residual(m,n) + srcl(m)/volj(n)
   !    end do
   !  end do
   !  !
   !  if (nf > nfbnd) then
   !    do n = r1,r2
   !      do m = 1,neqns
   !        residual(m,n) = residual(m,n) + srcr(m)/volj(n)
   !      end do
   !    end do
   !  end if
   !  !
   !end if
    !
  end do
  !
end subroutine correct_flux_project_cv
!
!###############################################################################
!
subroutine correct_flux_project_pv()
  !
  use geovar,  only : neqns,nfbnd,nface,n_sp1de,n_sp2dt, &
                      cells_on_face,pts_in_cell_ptr
  use geovar,  only : interp_s2f
  use metvar,  only : cormat,geofa
  use flowvar, only : vsp,residual
  use generic_mod, only : vsp2u_sp
  !
  !.. Local Scalars ..
  integer   :: k,m,n,nf,r1,r2,l1,l2
  integer   :: jl,jr,nl,nr,wl,wr
  integer   :: kcorbl,kcorbr,kcorl,kcorr
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:2)  :: unrm
  real(wp), dimension(1:neqns) :: vl,vr
  real(wp), dimension(1:neqns) :: dfl,dfr
  real(wp), dimension(1:neqns) :: dgl,dgr
  real(wp), dimension(1:neqns) :: fnl,fnr
  real(wp), dimension(1:neqns) :: fjl,fjr
  real(wp), dimension(1:neqns) :: fnc
  !
continue
  !
  ! Loop over all the faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  do nf = 1,nface
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kcorbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kcorbr = cells_on_face(6,nf) ! beginning flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    nl = l2-l1+1 ! number of solution points in left  cell
    nr = r2-r1+1 ! number of solution points in right cell
    !
    wl = merge( 1 , 2 , nl == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , nr == n_sp2dt ) ! check if jr is a tria or quad
    !
    ! Loop through the flux points on the current face
    !
    do k = 1,n_sp1de
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( geofa(:,k,nf) )
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( geofa(:,k,nf) )
      dsl = ds * cl(wl)
      dsr = ds * cl(wr)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kcorl = kcorbl-1 + k ! index for left  cell
      kcorr = kcorbr+1 - k ! index for right cell
      !
      ! Use the Lagrange polynomial to project the primitive variables from
      ! the left and right cells to the current flux point on the interface
      !
      do m = 1,neqns
        vl(m) = dot_product(interp_s2f(1:nl,kcorl,wl),vsp(m,l1:l2)) ! left  state
        vr(m) = dot_product(interp_s2f(1:nr,kcorr,wr),vsp(m,r1:r2)) ! right state
      end do
      !
      ! Get the numerical flux at the current flux point
      !
      if (flux_method == ROEFLUX) then
        fnc(:) = upwind_flux_cv ( unrm(:) , vsp2u_sp(vl) , vsp2u_sp(vr) )
      else if (flux_method == HLLC) then
        fnc(:) = hllc_upwind_flux_cv ( unrm(:) , vsp2u_sp(vl) , vsp2u_sp(vr) )
      else if (flux_method == LDFSS) then
        fnc(:) = ldfss_upwind_flux_cv ( unrm(:) , vsp2u_sp(vl) , vsp2u_sp(vr) )
      end if
      !
      ! Get left and right discontinuous fluxes at the current flux point
      !
      dfl(:) = f_invs_flx_pv(vl)
      dgl(:) = g_invs_flx_pv(vl)
      dfr(:) = f_invs_flx_pv(vr)
      dgr(:) = g_invs_flx_pv(vr)
      !
      ! Get the discontinuous flux function at the current
      ! interface point in the direction of the face normal
      !
      fnl(:) = dfl(:)*unrm(1) + dgl(:)*unrm(2)
      fnr(:) = dfr(:)*unrm(1) + dgr(:)*unrm(2)
      !
      ! The jump in the flux at the current interface point
      !
      fjl(:) = fnc(:) - fnl(:)
      fjr(:) = fnc(:) - fnr(:)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
      do n = l1,l2
        do m = 1,neqns
          residual(m,n) = residual(m,n) + dsl*fjl(m)*cormat(n-l1+1,kcorl,wl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jr
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      if (nf > nfbnd) then
        do n = r1,r2
          do m = 1,neqns
            residual(m,n) = residual(m,n) - dsr*fjr(m)*cormat(n-r1+1,kcorr,wr)
          end do
        end do
      end if
      !
    end do
    !
  end do
  !
end subroutine correct_flux_project_pv
!
!###############################################################################
!
subroutine correct_flux_project_flx()
  !
  use geovar,  only : neqns,nfbnd,nface,n_sp1de,n_sp2dt,n_sp2dq, &
                      cells_on_face,pts_in_cell_ptr,bface,xyz
  use geovar,  only : interp_s2f
  use metvar,  only : cormat,geofa,dx,dy
  use flowvar, only : usp,residual
  !
  !.. Local Scalars ..
  logical   :: is_1d_bc
  integer   :: k,m,n,nf,r1,r2,l1,l2
  integer   :: jl,jr,nl,nr,wl,wr,iu
  integer   :: kcorbl,kcorbr,kcorl,kcorr
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:2)  :: unrm
  real(wp), dimension(1:neqns) :: ul,ur
  real(wp), dimension(1:neqns) :: dfl,dfr
  real(wp), dimension(1:neqns) :: dgl,dgr
  real(wp), dimension(1:neqns) :: fnl,fnr
  real(wp), dimension(1:neqns) :: fjl,fjr
  real(wp), dimension(1:neqns) :: fnc
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lusp
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rusp
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lflxf
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lflxg
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rflxf
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rflxg
  !
continue
  !
  ! Loop over all the faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  do nf = 1,nface
    !
    ! Check if this is a boundary face needing a 1D treatment of the flux points
    !
    is_1d_bc = fals
    if (nf <= nfbnd) is_1d_bc = .not. any(bface(1,nf) == bc_comm)
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kcorbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kcorbr = cells_on_face(6,nf) ! beginning flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    nl = l2-l1+1 ! number of solution points in left  cell
    nr = r2-r1+1 ! number of solution points in right cell
    !
    wl = merge( 1 , 2 , nl == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , nr == n_sp2dt ) ! check if jr is a tria or quad
    !
    ! Compute the fluxes at each solution point and store the
    ! transpose of usp for the left cell to reduce cache misses
    !
    do k = l1,l2
      lflxf(k-l1+1,:) = f_invs_flx_cv( usp(:,k) )
      lflxg(k-l1+1,:) = g_invs_flx_cv( usp(:,k) )
      lusp (k-l1+1,:) = usp(:,k)
    end do
    !
    ! Store the transpose of usp for the right cell to reduce cache misses
    !
    do k = r1,r2
      rusp(k-r1+1,:) = usp(:,k)
    end do
    !
    ! Compute the right cell fluxes at each solution point for only the
    ! interior cells and communication/periodic boundary conditions
    !
   !1 format ("     lflx        rflx        fjl         fjr         fnc",/,*(i8))
   !4 format ("     luqp        ruqp",/,*(i8))
   !3 format (*(i6))
   !2 format (i6,*(es12.4))
   !iu = 105+10*mypnum
   !write (iu,1) nf,jl,jr
   !write (iu+1,3) nf,jl,jr
   !if (nf <= nfbnd) write (iu+1,3) (bface(m,nf),m=1,size(bface,1))
   !do k = r1,r2
   !  write (iu+1,2) k,(usp(m,k),m=1,neqns),(xyz(m,k),m=1,2)
   !end do
    if (.not. is_1d_bc) then
      do k = r1,r2
        rflxf(k-r1+1,:) = f_invs_flx_cv( usp(:,k) )
        rflxg(k-r1+1,:) = g_invs_flx_cv( usp(:,k) )
      end do
    end if
    !
    ! Loop through the flux points on the current face
    !
    do k = 1,n_sp1de
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( geofa(:,k,nf) )
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( geofa(:,k,nf) )
      dsl = ds * cl(wl)
      dsr = ds * cl(wr)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kcorl = kcorbl-1 + k ! index for left  cell
      kcorr = kcorbr+1 - k ! index for right cell
      !
      ! Use the Lagrange polynomial to project the conserved
      ! solution variables and the discontinuous fluxes from the
      ! left cell to the current flux point on the interface
      !
      do m = 1,neqns
        ul(m)  = dot_product(interp_s2f(1:nl,kcorl,wl), lusp(1:nl,m)) ! left solution
        dfl(m) = dot_product(interp_s2f(1:nl,kcorl,wl),lflxf(1:nl,m)) ! left f flux
        dgl(m) = dot_product(interp_s2f(1:nl,kcorl,wl),lflxg(1:nl,m)) ! left g flux
      end do
      !
      ! Now get the information at flux point on the right side
      !
      if (is_1d_bc) then
        !
        ! Copy the specified boundary condition in usp into ur and use
        ! this value to get the discontinuous fluxes at this flux point
        !
        !########################################################
        !########################################################
        ! NEED TO CHECK IF kcorr IS THE RIGHT INDEX THAT NEEDS TO
        ! BE USED HERE. I THINK kcorr VARIES FROM 1:3*n_sp1de
        ! AND NOT FROM 1:n_sp2dt OR 1:n_sp2dq. 
        !########################################################
        !########################################################
        !
       !write (*,*) "kcorr = ",kcorr
        ur(:) = rusp(kcorr,:)
        dfr(:) = f_invs_flx_cv( ur )
        dgr(:) = g_invs_flx_cv( ur )
        !
      else
        !
        ! Use the Lagrange polynomial to project the conserved solution
        ! variables and the discontinuous fluxes from the right (non-ghost)
        ! cell to the current flux point on the interface.
        !
        do m = 1,neqns
          ur(m)  = dot_product(interp_s2f(1:nr,kcorr,wr), rusp(1:nr,m)) ! right solution
          dfr(m) = dot_product(interp_s2f(1:nr,kcorr,wr),rflxf(1:nr,m)) ! right f flux
          dgr(m) = dot_product(interp_s2f(1:nr,kcorr,wr),rflxg(1:nr,m)) ! right g flux
        end do
        !
      end if
      !
      ! Get the numerical flux at the current flux point
      !
      if (flux_method == ROEFLUX) then
        fnc(:) = upwind_flux_cv ( unrm(:) , ul , ur )
      else if (flux_method == HLLC) then
        fnc(:) = hllc_upwind_flux_cv ( unrm(:) , ul , ur )
      else if (flux_method == LDFSS) then
        fnc(:) = ldfss_upwind_flux_cv ( unrm(:) , ul , ur )
      end if
      !
      ! Get the discontinuous flux function at the current
      ! interface point in the direction of the face normal
      !
      fnl(:) = dfl(:)*unrm(1) + dgl(:)*unrm(2)
      fnr(:) = dfr(:)*unrm(1) + dgr(:)*unrm(2)
      !
      ! The jump in the flux at the current interface point
      !
      fjl(:) = fnc(:) - fnl(:)
      fjr(:) = fnc(:) - fnr(:)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
     !write (iu,3) k,kcorl,kcorr
     !write (iu+1,3) k,kcorl,kcorr
     !do m = 1,neqns
     !  write (iu,2) fnl(m),fnr(m),fjl(m),fjr(m),fnc(m)
     !  write (iu+1,2) ul(m),ur(m)
     !end do
      do n = l1,l2
        do m = 1,neqns
          residual(m,n) = residual(m,n) + dsl*fjl(m)*cormat(n-l1+1,kcorl,wl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jr
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      if (nf > nfbnd) then
        do n = r1,r2
          do m = 1,neqns
            residual(m,n) = residual(m,n) - dsr*fjr(m)*cormat(n-r1+1,kcorr,wr)
          end do
        end do
      end if
      !
    end do
    !
  end do
  !
end subroutine correct_flux_project_flx
!
!###############################################################################
!
subroutine derflux1
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : neqns,ncell,n_sp2dt,n_sp2dq, &
                      pts_in_cell_ptr,write_variable
  use metvar,  only : dx,dy,dr_mat,ds_mat
  use flowvar, only : usp,vsp,residual
  !
  !.. Local Scalars ..
  integer :: k,m,n,nw,n1,n2,np
  !
  !.. Local Arrays ..
  real(wp), dimension(1:neqns)           :: ff,gg
  real(wp), dimension(1:n_sp2dq,1:neqns) :: ftdc,gtdc
  !
continue
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    np = n2-n1+1 ! number of solution points in this cell
    !
    nw = merge( 1 , 2 , np == n_sp2dt ) ! check if cell n is tria or quad
    !
    ftdc(:,:) = zero ! zero for dot product later with hybrid tria/quad
    gtdc(:,:) = zero ! zero for dot product later with hybrid tria/quad
    !
    ! Get the flux at each solution point in the current cell
    !
    do k = n1,n2
      !
      ff(:) = f_invs_flx_cv( usp(:,k) ) ! x-directional flux
      gg(:) = g_invs_flx_cv( usp(:,k) ) ! y-directional flux
      !
      ftdc(k-n1+1,:) =  ff(:)*dy(k)%ds - gg(:)*dx(k)%ds
      gtdc(k-n1+1,:) =  ff(:)*dy(k)%dr - gg(:)*dx(k)%dr
      !
    end do
    !
    ! Compute the derivatives of the fluxes at each
    ! solution point and add to the residual
    !
    do k = 1,np
      do m = 1,neqns
        residual(m,k+n1-1) = dot_product(dr_mat(1:np,k,nw),ftdc(1:np,m)) - &
                             dot_product(ds_mat(1:np,k,nw),gtdc(1:np,m))
      end do
    end do
    !
  end do
  !
end subroutine derflux1
!
!###############################################################################
!
subroutine derflux2
  !
  !... this subroutine get derivative of fluxes (no interaction) by
  !... the chain rule (non-conservative formula,
  !...  eq. 95 p. 124 Thompson et al. Grid Generation book)
  !... for each cell, get uxi and uet at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : neqns,ncell,n_sp2dt,n_sp2dq, &
                      pts_in_cell_ptr
  use metvar,  only : dx,dy,dr_mat,ds_mat
  use flowvar, only : usp,residual
  !
  !.. Local Scalars ..
  integer :: k,m,n,n1,n2,nw,np,ip
  !
  !.. Local Arrays ..
  real(wp), dimension(1:neqns)           :: src
  real(wp), dimension(1:neqns,1:neqns)   :: dfdu,dgdu
  real(wp), dimension(1:neqns,1:neqns)   :: c_f,c_g
  real(wp), dimension(1:neqns,1:n_sp2dq) :: ucr,ucs
  real(wp), dimension(1:n_sp2dq,1:neqns) :: tusp
  !
continue
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    np = n2-n1+1 ! number of solution points in this cell
    !
    nw = merge(1,2,np == n_sp2dt) ! check if cell n is tria or quad
    !
    ! Store the transpose of usp for this cell to reduce cache misses
    !
    do k = 1,np
      do m = 1,neqns
        tusp(k,m) = usp(m,k+n1-1)
      end do
    end do
    !
    ! Compute the derivative of the conserved variables at each solution point
    !
    do k = 1,np
      do m = 1,neqns
        ucr(m,k) = dot_product( dr_mat(1:np,k,nw) , tusp(1:np,m) )
        ucs(m,k) = dot_product( ds_mat(1:np,k,nw) , tusp(1:np,m) )
      end do
    end do
   !!
   !! Initialize the conservation source term
   !!
   !src(:) = zero
    !
    ! Compute the derivatives of the fluxes at each solution point
    !
    do k = n1,n2
      !
      ip = k-n1+1 ! local index of current solution point
      !
      ! Compute the Jacobian of the fluxes at the current sol-pt
      !
      dfdu(:,:) = f_invs_flx_jacobian_cv( usp(:,k) )
      dgdu(:,:) = g_invs_flx_jacobian_cv( usp(:,k) )
      !
      ! Compute the Jacobian in the direction of the face normal
      !
      c_f(:,:) = dy(k)%ds*dfdu(:,:) - dx(k)%ds*dgdu(:,:)
      c_g(:,:) = dy(k)%dr*dfdu(:,:) - dx(k)%dr*dgdu(:,:)
      !
      ! Use the chain rule and add to the residual
      !
      do m = 1,neqns
        residual(m,k) = dot_product(ucr(:,ip),c_f(:,m)) - &
                        dot_product(ucs(:,ip),c_g(:,m))
      end do
   !  !
   !  ! From ZJ, to ensure conservation with chain-rule approach
   !  !
   !  do m = 1,neqns
   !    src(m) = src(m) + residual(m,k)*jacobian(k)*weights2d(ip,nw)
   !  end do
      !
    end do
   !!
   !! Now at the volume conservation source term
   !! to each solution point in this cell
   !!
   !do k = n1,n2
   !  do m = 1,neqns
   !    residual(m,k) = residual(m,k) - src(m)/volj(k)
   !  end do
   !end do
    !
  end do
  !
end subroutine derflux2
!
!###############################################################################
!
subroutine derflux3
  !
  !... this subroutine get derivative of fluxes (no interaction) by
  !... applying the chain rule on UTD
  !... for each cell, get uxi and uet at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : neqns,ncell,n_sp2dt,n_sp2dq, &
                      pts_in_cell_ptr
  use metvar,  only : dx,dy,dr_mat,ds_mat, &
                      jacobian,dJ
  use flowvar, only : usp,residual
  !
  !.. Local Scalars ..
  integer :: k,m,n,nw,n1,n2,np,ip
  !
  !.. Local Arrays ..
  real(wp), dimension(1:neqns,1:neqns)      :: dfdu,dgdu
  real(wp), dimension(1:neqns,1:neqns)      :: c_f,c_g
  real(wp), dimension(1:neqns,1:n_sp2dq) :: ucr,ucs
  real(wp), dimension(1:n_sp2dq,1:neqns) :: uja
  real(wp), dimension(1:neqns,1:n_sp2dq) :: ujar,ujas
  !
continue
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    np = n2-n1+1
    !
    nw = merge(1,2,np == n_sp2dt) ! check if cell n is tria or quad
    !
    ! Compute u_tilde (= conserved variables * jacobian) at each solution point
    !
    do k = n1,n2
      uja(k-n1+1,:) = usp(:,k)*jacobian(k)
    end do
    !
    ! Compute the derivative of u_tilde at each solution point
    !
    do k = 1,np
      do m = 1,neqns
        ujar(m,k) = dot_product( dr_mat(1:np,k,nw) , uja(1:np,m) )
        ujas(m,k) = dot_product( ds_mat(1:np,k,nw) , uja(1:np,m) )
      end do
    end do
    !
    ! Compute the derivatives of the fluxes at each solution point
    !
    do k = n1,n2
      !
      ip = k-n1+1 ! local index of current solution point
      !
      ! Compute the derivative of the conserved
      ! variables for the current sol-pt
      !
      ! = 1/J * { [d(u*J)/dr] - [u * d(J)/dr] }
      ! = 1/J * { [J * d(u)/dr] + [u * d(J)/dr] - [u * d(J)/dr] }
      ! = 1/J * { [J * d(u)/dr] }
      ! = d(u)/dr
      !
      ucr(:,ip) = (one/jacobian(k)) * (ujar(:,ip) - usp(:,k)*dJ(k)%dr)
      ucs(:,ip) = (one/jacobian(k)) * (ujas(:,ip) - usp(:,k)*dJ(k)%ds)
      !
      ! Compute the Jacobian of the fluxes at the current sol-pt
      !
      dfdu(:,:) = f_invs_flx_jacobian_cv( usp(:,k) )
      dgdu(:,:) = g_invs_flx_jacobian_cv( usp(:,k) )
      !
      ! Compute the Jacobian in the direction of the face normal
      !
      c_f(:,:) = dy(k)%ds*dfdu(:,:) - dx(k)%ds*dgdu(:,:)
      c_g(:,:) = dy(k)%dr*dfdu(:,:) - dx(k)%dr*dgdu(:,:)
      !
      ! Use the chain rule and add to the residual
      !
      do m = 1,neqns
        residual(m,k) = dot_product(ucr(:,ip), c_f(:,m)) - &
                        dot_product(ucs(:,ip), c_g(:,m))
      end do
      !
    end do
    !
  end do
  !
end subroutine derflux3
!
!###############################################################################
!
subroutine dg_strong_form_old()
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : n_sp1de,n_sp2dt,n_sp2dq, &
                      neqns,ncell,nface,nfbnd, &
                      pts_in_cell_ptr,cells_on_face
  use geovar,  only : interp_s2f
  use metvar,  only : cormat_dg,geofa, &
                      dx,dy,dr_mat,ds_mat
  use flowvar, only : usp,residual
  !
  !.. Local Scalars ..
  integer   :: k,m,n,nf,r1,r2,l1,l2
  integer   :: jl,jr,nl,nr,nw,n1,n2,np,wl,wr
  integer   :: kcorbl,kcorbr,kcorl,kcorr
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:2)  :: unrm
  real(wp), dimension(1:neqns) :: ul,ur
  real(wp), dimension(1:neqns) :: dfl,dfr
  real(wp), dimension(1:neqns) :: dgl,dgr
  real(wp), dimension(1:neqns) :: fnl,fnr
  real(wp), dimension(1:neqns) :: fjl,fjr
  real(wp), dimension(1:neqns) :: fnc
  real(wp), dimension(1:neqns) :: ff,gg
  real(wp), dimension(1:n_sp2dq,1:neqns) :: ftdc
  real(wp), dimension(1:n_sp2dq,1:neqns) :: gtdc
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lusp
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rusp
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lflxf
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lflxg
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rflxf
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rflxg
  !
continue
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    np = n2-n1+1 ! number of solution points in this cell
    !
    nw = merge( 1 , 2 , np == n_sp2dt ) ! check if cell n is tria or quad
    !
    ftdc(:,:) = zero ! zero for dot product later with hybrid tria/quad
    gtdc(:,:) = zero ! zero for dot product later with hybrid tria/quad
    !
    ! Get the flux at each solution point in the current cell
    !
    do k = n1,n2
      !
      ff(:) = f_invs_flx_cv( usp(:,k) ) ! x-directional flux
      gg(:) = g_invs_flx_cv( usp(:,k) ) ! y-directional flux
      !
      ftdc(k-n1+1,:) =  ff(:)*dy(k)%ds - gg(:)*dx(k)%ds
      gtdc(k-n1+1,:) =  ff(:)*dy(k)%dr - gg(:)*dx(k)%dr
      !
    end do
    !
    ! Compute the derivatives of the fluxes at each
    ! solution point and add to the residual
    !
    do k = 1,np
      do m = 1,neqns
        residual(m,k+n1-1) = dot_product(dr_mat(1:np,k,nw),ftdc(1:np,m)) - &
                             dot_product(ds_mat(1:np,k,nw),gtdc(1:np,m))
      end do
    end do
    !
  end do
  !
  ! Loop over all the faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  do nf = 1,nface
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kcorbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kcorbr = cells_on_face(6,nf) ! beginning flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    nl = l2-l1+1 ! number of solution points in left  cell
    nr = r2-r1+1 ! number of solution points in right cell
    !
    wl = merge( 1 , 2 , nl == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , nr == n_sp2dt ) ! check if jr is a tria or quad
    !
    ! Compute the fluxes at each solution point and store the
    ! transpose of usp for the left cell to reduce cache misses
    !
    do k = l1,l2
      lflxf(k-l1+1,:) = f_invs_flx_cv( usp(:,k) )
      lflxg(k-l1+1,:) = g_invs_flx_cv( usp(:,k) )
      lusp (k-l1+1,:) = usp(:,k)
    end do
    !
    ! Compute the fluxes at each solution point and store the
    ! transpose of usp for the right cell to reduce cache misses
    ! 
    do k = r1,r2
      rflxf(k-r1+1,:) = f_invs_flx_cv( usp(:,k) )
      rflxg(k-r1+1,:) = g_invs_flx_cv( usp(:,k) )
      rusp (k-r1+1,:) = usp(:,k)
    end do
    !
    ! Loop through the flux points on the current face
    !
    do k = 1,n_sp1de
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( geofa(:,k,nf) )
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( geofa(:,k,nf) )
      dsl = ds * cl(wl)
      dsr = ds * cl(wr)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kcorl = kcorbl-1 + k ! index for left  cell
      kcorr = kcorbr+1 - k ! index for right cell
      !
      ! Use the Lagrange polynomial to project the conservative variables from
      ! the left and right cells to the current flux point on the interface
      ! in order to compute the numerical flux at the current flux point
      !
      do m = 1,neqns
        ul(m) = dot_product(interp_s2f(1:nl,kcorl,wl),lusp(1:nl,m)) ! left  state
        ur(m) = dot_product(interp_s2f(1:nr,kcorr,wr),rusp(1:nr,m)) ! right state
      end do
      !
      ! Get the numerical flux at the current flux point
      !
      if (flux_method == ROEFLUX) then
        fnc(:) = upwind_flux_cv ( unrm(:) , ul , ur )
      else if (flux_method == HLLC) then
        fnc(:) = hllc_upwind_flux_cv ( unrm(:) , ul , ur )
      else if (flux_method == LDFSS) then
        fnc(:) = ldfss_upwind_flux_cv ( unrm(:) , ul , ur )
      end if
      !
      ! Get left and right discontinuous fluxes at the current flux point
      !
      do m = 1,neqns
        dfl(m) = dot_product(interp_s2f(1:nl,kcorl,wl),lflxf(1:nl,m))
        dgl(m) = dot_product(interp_s2f(1:nl,kcorl,wl),lflxg(1:nl,m))
        dfr(m) = dot_product(interp_s2f(1:nr,kcorr,wr),rflxf(1:nr,m))
        dgr(m) = dot_product(interp_s2f(1:nr,kcorr,wr),rflxg(1:nr,m))
      end do
      !
      ! Get the discontinuous flux function at the current
      ! interface point in the direction of the face normal
      !
      fnl(:) = dfl(:)*unrm(1) + dgl(:)*unrm(2)
      fnr(:) = dfr(:)*unrm(1) + dgr(:)*unrm(2)
      !
      ! The jump in the flux at the current interface point
      !
      fjl(:) = fnc(:) - fnl(:)
      fjr(:) = fnc(:) - fnr(:)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
      do n = l1,l2
        do m = 1,neqns
          residual(m,n) = residual(m,n) + dsl*fjl(m)*cormat_dg(n-l1+1,kcorl,wl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jr
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      if (nf > nfbnd) then
        do n = r1,r2
          do m = 1,neqns
            residual(m,n) = residual(m,n) - dsr*fjr(m)*cormat_dg(n-r1+1,kcorr,wr)
          end do
        end do
      end if
      !
    end do
    !
  end do
  !
end subroutine dg_strong_form_old
!
!###############################################################################
!
subroutine dg_strong_form()
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : n_sp1de,n_sp2dt,n_sp2dq, &
                      neqns,ncell,nface,nfbnd, &
                      pts_in_cell_ptr,cells_on_face
  use geovar,  only : n_qp1de,n_qp2dq,n_qp2dt
  use geovar,  only : interp_f2g,interp_s2q,qdtwts_2d
  use geovar,  only : interp_s2f,interp_s2g
  use metvar,  only : cormat_dg,geofa,dqr_mat,dqs_mat, &
                      dx,dy,dr_mat,ds_mat,jacobian
  use flowvar, only : usp,residual
  !
  !.. Local Scalars ..
  integer   :: k,m,n,nf,r1,r2,l1,l2,nsp,nqp
  integer   :: jl,jr,nl,nr,nw,n1,n2,np,wl,wr
  integer   :: kcorbl,kcorbr,kcorl,kcorr
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:2)  :: unrm
  real(wp), dimension(1:neqns) :: ul,ur
  real(wp), dimension(1:neqns) :: dfl,dfr
  real(wp), dimension(1:neqns) :: dgl,dgr
  real(wp), dimension(1:neqns) :: fnl,fnr
  real(wp), dimension(1:neqns) :: fjl,fjr
  real(wp), dimension(1:neqns) :: fnc
  real(wp), dimension(1:neqns) :: ff,gg
  real(wp), dimension(1:n_qp2dq) :: dxdr
  real(wp), dimension(1:n_qp2dq) :: dxds
  real(wp), dimension(1:n_qp2dq) :: dydr
  real(wp), dimension(1:n_qp2dq) :: dyds
  real(wp), dimension(1:neqns,1:n_qp2dq) :: uqp
  real(wp), dimension(1:n_qp2dq,1:neqns) :: ftdc
  real(wp), dimension(1:n_qp2dq,1:neqns) :: gtdc
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lusp
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rusp
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lflxf
  real(wp), dimension(1:n_sp2dq,1:neqns) :: lflxg
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rflxf
  real(wp), dimension(1:n_sp2dq,1:neqns) :: rflxg
  real(wp), dimension(1:n_qp2dq) :: qjac
  real(wp), dimension(1:2,1:n_qp1de) :: facenorms
  !
continue
  !
  nsp = n_sp1de
  nqp = n_qp1de
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    np = n2-n1+1 ! number of solution points in this cell
    !
    nw = merge( 1 , 2 , np == n_sp2dt ) ! check if cell n is tria or quad
    !
    do k = 1,n_qp2dq
      !
      ! Get the conservative solution at the quadrature points
      !
      do m = 1,neqns
        uqp(m,k) = dot_product( interp_s2q(1:np,k,nw) , usp(m,n1:n2) )
      end do
      !
      ! Get the metrics for this cell
      !
      dxdr(k) = dot_product( interp_s2q(1:np,k,nw) , dx(n1:n2)%dr )
      dxds(k) = dot_product( interp_s2q(1:np,k,nw) , dx(n1:n2)%ds )
      dydr(k) = dot_product( interp_s2q(1:np,k,nw) , dy(n1:n2)%dr )
      dyds(k) = dot_product( interp_s2q(1:np,k,nw) , dy(n1:n2)%ds )
      qjac(k) = dot_product( interp_s2q(1:np,k,nw) , jacobian(n1:n2) )
      !
    end do
    !
    !
    ftdc(:,:) = zero ! zero for dot product later with hybrid tria/quad
    gtdc(:,:) = zero ! zero for dot product later with hybrid tria/quad
    !
    ! Get the flux at each solution point in the current cell
    !
    do k = 1,n_qp2dq
      !
      ff(:) = f_invs_flx_cv( uqp(:,k) ) ! x-directional flux
      gg(:) = g_invs_flx_cv( uqp(:,k) ) ! y-directional flux
      !
      ftdc(k,:) =  ff(:)*dyds(k) - gg(:)*dxds(k)
      gtdc(k,:) =  ff(:)*dydr(k) - gg(:)*dxdr(k)
      !
      ftdc(k,:) = ftdc(k,:) * qdtwts_2d(k,nw)
      gtdc(k,:) = gtdc(k,:) * qdtwts_2d(k,nw)
      !
    end do
    !
    ! Compute the derivatives of the fluxes at each
    ! solution point and add to the residual
    !
    do k = 1,np
      do m = 1,neqns
        residual(m,k+n1-1) = dot_product(dqr_mat(:,k,nw),ftdc(:,m)) - &
                             dot_product(dqs_mat(:,k,nw),gtdc(:,m))
      end do
    end do
    !
  end do
  !
  ! Loop over all the faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  do nf = 1,nface
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kcorbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kcorbr = cells_on_face(6,nf) ! beginning flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    nl = l2-l1+1 ! number of solution points in left  cell
    nr = r2-r1+1 ! number of solution points in right cell
    !
    wl = merge( 1 , 2 , nl == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , nr == n_sp2dt ) ! check if jr is a tria or quad
    !
    ! Compute the fluxes at each solution point and store the
    ! transpose of usp for the left cell to reduce cache misses
    !
    do k = l1,l2
      lflxf(k-l1+1,:) = f_invs_flx_cv( usp(:,k) )
      lflxg(k-l1+1,:) = g_invs_flx_cv( usp(:,k) )
      lusp (k-l1+1,:) = usp(:,k)
    end do
    !
    ! Compute the fluxes at each solution point and store the
    ! transpose of usp for the right cell to reduce cache misses
    ! 
    do k = r1,r2
      rflxf(k-r1+1,:) = f_invs_flx_cv( usp(:,k) )
      rflxg(k-r1+1,:) = g_invs_flx_cv( usp(:,k) )
      rusp (k-r1+1,:) = usp(:,k)
    end do
    !
    ! Loop through the flux points on the current face
    !
    do k = 1,n_qp1de
      do n = 1,2
        facenorms(n,k) = dot_product( interp_f2g(1:nsp,k) , geofa(n,1:nsp,nf) )
      end do
    end do
    !
   !do k = 1,n_sp1de
    do k = 1,n_qp1de
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( facenorms(:,k) )
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( facenorms(:,k) )
      dsl = ds * cl(wl)
      dsr = ds * cl(wr)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kcorl = kcorbl-1 + k ! index for left  cell
      kcorr = kcorbr+1 - k ! index for right cell
      !
      ! Use the Lagrange polynomial to project the conservative variables from
      ! the left and right cells to the current flux point on the interface
      ! in order to compute the numerical flux at the current flux point
      !
      do m = 1,neqns
       !ul(m) = dot_product(interp_s2f(1:nl,kcorl,wl),lusp(1:nl,m)) ! left  state
        ul(m) = dot_product(interp_s2g(1:nl,kcorl,wl),lusp(1:nl,m)) ! left  state
       !ur(m) = dot_product(interp_s2f(1:nr,kcorr,wr),rusp(1:nr,m)) ! right state
        ur(m) = dot_product(interp_s2g(1:nr,kcorr,wr),rusp(1:nr,m)) ! right state
      end do
      !
      ! Get the numerical flux at the current flux point
      !
      if (flux_method == ROEFLUX) then
        fnc(:) = upwind_flux_cv ( unrm(:) , ul , ur )
      else if (flux_method == HLLC) then
        fnc(:) = hllc_upwind_flux_cv ( unrm(:) , ul , ur )
      else if (flux_method == LDFSS) then
        fnc(:) = ldfss_upwind_flux_cv ( unrm(:) , ul , ur )
      end if
      !
      ! Get left and right discontinuous fluxes at the current flux point
      !
      do m = 1,neqns
       !dfl(m) = dot_product(interp_s2f(1:nl,kcorl,wl),lflxf(1:nl,m))
        dfl(m) = dot_product(interp_s2g(1:nl,kcorl,wl),lflxf(1:nl,m))
       !dgl(m) = dot_product(interp_s2f(1:nl,kcorl,wl),lflxg(1:nl,m))
        dgl(m) = dot_product(interp_s2g(1:nl,kcorl,wl),lflxg(1:nl,m))
       !dfr(m) = dot_product(interp_s2f(1:nr,kcorr,wr),rflxf(1:nr,m))
        dfr(m) = dot_product(interp_s2g(1:nr,kcorr,wr),rflxf(1:nr,m))
       !dgr(m) = dot_product(interp_s2f(1:nr,kcorr,wr),rflxg(1:nr,m))
        dgr(m) = dot_product(interp_s2g(1:nr,kcorr,wr),rflxg(1:nr,m))
      end do
      !
      ! Get the discontinuous flux function at the current
      ! interface point in the direction of the face normal
      !
      fnl(:) = dfl(:)*unrm(1) + dgl(:)*unrm(2)
      fnr(:) = dfr(:)*unrm(1) + dgr(:)*unrm(2)
      !
      ! The jump in the flux at the current interface point
      !
      fjl(:) = fnc(:) - fnl(:)
      fjr(:) = fnc(:) - fnr(:)
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
      do n = l1,l2
        do m = 1,neqns
          residual(m,n) = residual(m,n) + dsl*fjl(m)*cormat_dg(n-l1+1,kcorl,wl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jr
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      if (nf > nfbnd) then
        do n = r1,r2
          do m = 1,neqns
            residual(m,n) = residual(m,n) - dsr*fjr(m)*cormat_dg(n-r1+1,kcorr,wr)
          end do
        end do
      end if
      !
    end do
    !
  end do
  !
end subroutine dg_strong_form
!
!###############################################################################
!
subroutine dg_weak_form()
  !
  !... this subroutine get derivative of fluxes (no interaction)
  !... for each cell, get f, g, ftd, gtd at the solution points
  !
  !.. Use Statements ..
  use geovar,  only : n_sp1de,n_sp2dt,n_sp2dq, &
                      neqns,ncell,nface,nfbnd,bface, &
                      cells_on_face,pts_in_cell_ptr
  use geovar,  only : interp_s2f
  use metvar,  only : cormat_dg,geofa, &
                      dx,dy,dr_hat,ds_hat
  use flowvar, only : usp,residual
  !
  !.. Local Scalars ..
  logical   :: is_1d_bc
  integer   :: k,l,m,n,nf,r1,r2,l1,l2
  integer   :: jl,jr,nl,nr,nw,n1,n2,np,wl,wr
  integer   :: kcorbl,kcorbr,kcorl,kcorr
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:2)            :: unrm
  real(wp), dimension(1:neqns)           :: fnc
  real(wp), dimension(1:neqns)           :: ul,ur
  real(wp), dimension(1:neqns)           :: ff,gg
  real(wp), dimension(1:neqns)           :: fjl,fjr
  real(wp), dimension(1:n_sp2dq,1:neqns) :: ftdc,gtdc
  !
continue
  !
  ! Loop over the cells and get the derivative of the
  ! fluxes at each solution point within each cell
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    np = n2-n1+1 ! number of solution points in this cell
    !
    nw = merge( 1 , 2 , np == n_sp2dt ) ! check if cell n is tria or quad
    !
    ftdc(:,:) = zero ! zero for dot product later with hybrid tria/quad
    gtdc(:,:) = zero ! zero for dot product later with hybrid tria/quad
    !
    ! Get the flux at each solution point in the current cell
    !
    do k = n1,n2
      !
      ff(:) = f_invs_flx_cv( usp(:,k) ) ! x-directional flux
      gg(:) = g_invs_flx_cv( usp(:,k) ) ! y-directional flux
      !
      ftdc(k-n1+1,:) =  ff(:)*dy(k)%ds - gg(:)*dx(k)%ds
      gtdc(k-n1+1,:) =  ff(:)*dy(k)%dr - gg(:)*dx(k)%dr
      !
    end do
    !
    ! Compute the derivatives of the fluxes at each
    ! solution point and add to the residual
    !
    do k = 1,np
      do m = 1,neqns
        residual(m,k+n1-1) = dot_product(dr_hat(1:np,k,nw),ftdc(1:np,m)) - &
                             dot_product(ds_hat(1:np,k,nw),gtdc(1:np,m))
      end do
    end do
    !
  end do
  !
  ! Loop over all the faces and reconstruct the
  ! interface fluxes at each cell solution point
  !
  do nf = 1,nface
    !
    ! Check if this is a boundary face needing a 1D treatment of the flux points
    !
    is_1d_bc = fals
    if (nf <= nfbnd) is_1d_bc = .not. any(bface(1,nf) == bc_comm)
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kcorbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kcorbr = cells_on_face(6,nf) ! beginning flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    nl = l2-l1+1 ! number of solution points in left  cell
    nr = r2-r1+1 ! number of solution points in right cell
    !
    wl = merge( 1 , 2 , nl == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , nr == n_sp2dt ) ! check if jr is a tria or quad
    !
    ! Loop through the flux points on the current face
    !
    do k = 1,n_sp1de
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( geofa(:,k,nf) )
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( geofa(:,k,nf) )
      dsl = ds * cl(wl)
      dsr = ds * cl(wr)
      !
      ! Get the index for the correction function
      ! corresponding to the current interface flux point
      !
      kcorl = kcorbl-1 + k ! index for left  cell
      kcorr = kcorbr+1 - k ! index for right cell
      !
      ! Use the Lagrange polynomial to project the conservative variables from
      ! the left and right cells to the current flux point on the interface
      !
      do m = 1,neqns
        ul(m) = dot_product(interp_s2f(1:nl,kcorl,wl),usp(m,l1:l2)) ! left  state
        ur(m) = dot_product(interp_s2f(1:nr,kcorr,wr),usp(m,r1:r2)) ! right state
      end do
      !
      ! Get the numerical flux at the current flux point
      !
      if (flux_method == ROEFLUX) then
        fnc(:) = upwind_flux_cv ( unrm(:) , ul , ur )
      else if (flux_method == HLLC) then
        fnc(:) = hllc_upwind_flux_cv ( unrm(:) , ul , ur )
      else if (flux_method == LDFSS) then
        fnc(:) = ldfss_upwind_flux_cv ( unrm(:) , ul , ur )
      end if
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
      do n = l1,l2
        do m = 1,neqns
          residual(m,n) = residual(m,n) + dsl*fnc(m)*cormat_dg(n-l1+1,kcorl,wl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jr
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      if (nf > nfbnd) then
        do n = r1,r2
          do m = 1,neqns
            residual(m,n) = residual(m,n) - dsr*fnc(m)*cormat_dg(n-r1+1,kcorr,wr)
          end do
        end do
      end if
      !
    end do
    !
  end do
  !
end subroutine dg_weak_form
!
!###############################################################################
!###############################################################################
!###############################################################################
!
! Procedures for inviscid fluxes
!
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function normal_invs_flx_cv(uu,rnrm) result(return_value)
  !
  use geovar, only : neqns,nr
  !
  real(wp), dimension(1:neqns), intent(in) :: uu
  real(wp), dimension(1:nr),    intent(in) :: rnrm
  !
  real(wp), dimension(1:neqns) :: return_value
  !
  real(wp), dimension(1:nr)    :: unrm
  real(wp), dimension(1:neqns) :: f_flx,g_flx
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
  use geovar, only : neqns
  use ovar, only : gam
  !
  real(wp), dimension(1:neqns), intent(in) :: uu
  !
  real(wp), dimension(1:neqns) :: return_value
  !
  real(wp) :: pres
  !
continue
  !
  pres = (gam-one)*(uu(neqns) - half*selfdot(uu(2:neqns-1))/uu(1))
  !
  return_value(1) = uu(2)
  return_value(2:neqns-1) = uu(2)*uu(2:neqns-1)/uu(1)
  return_value(2) = return_value(2) + pres
  return_value(neqns) = uu(2)/uu(1) * (uu(neqns) + pres)
  !
end function f_invs_flx_cv
!
!###############################################################################
!
pure function g_invs_flx_cv(uu) result(return_value)
  !
  use geovar, only : neqns
  use ovar, only : gam
  !
  real(wp), dimension(1:neqns), intent(in) :: uu
  !
  real(wp), dimension(1:neqns) :: return_value
  !
  real(wp) :: pres
  !
continue
  !
  pres = (gam-one)*(uu(neqns) - half*selfdot(uu(2:neqns-1))/uu(1))
  !
  return_value(1) = uu(3)
  return_value(2:neqns-1) = uu(3)*uu(2:neqns-1)/uu(1)
  return_value(3) = return_value(3) + pres
  return_value(neqns) = uu(3)/uu(1) * (uu(neqns) + pres)
  !
end function g_invs_flx_cv
!
!###############################################################################
!
pure function f_invs_flx_pv(vv) result(return_value) 
  !
  use geovar, only : neqns
  use ovar, only : gam
  !
  real(wp), dimension(1:neqns), intent(in) :: vv
  !
  real(wp), dimension(1:neqns) :: return_value
  !
  real(wp) :: ener
  !
continue
  !
  ener = vv(neqns)/(gam-one) + half*vv(1)*selfdot(vv(2:neqns-1))
  !
  return_value(1) = vv(1)*vv(2)
  return_value(2:neqns-1) = vv(1)*vv(2:neqns-1) * vv(2)
  return_value(2) = return_value(2) + vv(neqns)
  return_value(neqns) = (ener + vv(neqns)) * vv(2)
  !
end function f_invs_flx_pv
!
!###############################################################################
!
pure function g_invs_flx_pv(vv) result(return_value)
  !
  use geovar, only : neqns
  use ovar, only : gam
  !
  real(wp), dimension(1:neqns), intent(in) :: vv
  !
  real(wp), dimension(1:neqns) :: return_value
  !
  real(wp) :: ener
  !
continue
  !
  ener = vv(neqns)/(gam-one) + half*vv(1)*selfdot(vv(2:neqns-1))
  !
  return_value(1) = vv(1)*vv(3)
  return_value(2:neqns-1) = vv(1)*vv(2:neqns-1) * vv(3)
  return_value(3) = return_value(3) + vv(neqns)
  return_value(neqns) = (ener + vv(neqns)) * vv(3)
  !
end function g_invs_flx_pv
!
!###############################################################################
!
pure function f_invs_flx_jacobian_cv(uu) result(return_value)
  !
  use geovar, only : neqns
  use ovar, only : gam
  !
  real(wp), dimension(1:neqns), intent(in) :: uu
  !
  real(wp), dimension(1:neqns,1:neqns) :: return_value
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
  use geovar, only : neqns
  use ovar, only : gam
  !
  real(wp), dimension(1:neqns), intent(in) :: uu
  !
  real(wp), dimension(1:neqns,1:neqns) :: return_value
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
  use geovar, only : neqns
  use ovar, only : gam
  !
  real(wp), dimension(1:neqns), intent(in) :: vv
  !
  real(wp), dimension(1:neqns,1:neqns) :: return_value
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
  use geovar, only : neqns
  use ovar, only : gam
  !
  real(wp), dimension(1:neqns), intent(in) :: vv
  !
  real(wp), dimension(1:neqns,1:neqns) :: return_value
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
!###############################################################################
!###############################################################################
!
! Procedures for viscous fluxes
!
!###############################################################################
!###############################################################################
!###############################################################################
!
pure function visc_flx_cv(uu,dcvdx) result(return_value) 
  !
  !.. Use Statements ..
  use geovar, only : nr,neqns
  use ovar,   only : cpref_nd,Pr
  !
  use generic_mod, only : viscosity_cv_sp
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:neqns),      intent(in) :: uu
  real(wp), dimension(1:nr,1:neqns), intent(in) :: dcvdx
  !
  !.. Function Result ..
  real(wp), dimension(1:neqns,1:nr) :: return_value
  !
  !.. Local Scalars ..
  integer  :: nq,vb,ve
  real(wp) :: dudx,dudy,dvdx,dvdy
  real(wp) :: mvsc,mcdt,dTdx,dTdy
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: vel
  real(wp), dimension(1:nr,1:neqns) :: dpvdx
  !
continue
  !
  nq = neqns ! alias for neqns
  vb = 2     ! first momentum equation
  ve = nq-1  ! last  momentum equation 
  !
  ! Get the gradients of the primitive variables using the 
  ! gradients of the conservative variables and the product rule
  !
  dpvdx(1:nr,1:nq) = grad_cv_to_grad_pv_sp(uu,dcvdx)
  !
  ! Get the velocity vector
  !
  vel(:) = uu(vb:ve)/uu(1)
  !
  mvsc = viscosity_cv_sp(uu)
  mcdt = cpref_nd * mvsc / Pr
  !
  dudx = dpvdx(1,2)
  dudy = dpvdx(2,2)
  dvdx = dpvdx(1,3)
  dvdy = dpvdx(2,3)
  dTdx = dpvdx(1,nq)
  dTdy = dpvdx(2,nq)
  !
  ! x-directional viscous fluxes
  !
  return_value( 1,1) = zero
  return_value( 2,1) = mvsc * two3 * (two*dudx - dvdy)
  return_value( 3,1) = mvsc * (dudy + dvdx)
  return_value(nq,1) = mcdt*dTdx + dot_product( vel , return_value(vb:ve,1) )
  !
  ! y-directional viscous fluxes
  !
  return_value( 1,2) = zero
  return_value( 2,2) = mvsc * (dudy + dvdx)
  return_value( 3,2) = mvsc * two3 * (two*dvdy - dudx)
  return_value(nq,2) = mcdt*dTdy + dot_product( vel , return_value(vb:ve,2) )
  !
end function visc_flx_cv
!
!###############################################################################
!
pure function common_visc_flx_cv(unrm,ul,duldx,ur,durdx) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : nr,neqns
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),         intent(in) :: unrm
  real(wp), dimension(1:neqns),      intent(in) :: ul,ur
  real(wp), dimension(1:nr,1:neqns), intent(in) :: duldx,durdx
  !
  !.. Function Result ..
  real(wp), dimension(1:neqns) :: return_value
  !
  !.. Local Scalars ..
  integer  :: l,m
  real(wp) :: vjump
  !
  !.. Local Arrays ..
  real(wp), dimension(1:neqns,1:nr) :: lvflx
  real(wp), dimension(1:neqns,1:nr) :: rvflx
  real(wp), dimension(1:nr,1:neqns) :: vflux
  !
continue
  !
  ! Get the left and right viscous fluxes
  !
  lvflx(1:neqns,1:nr) = visc_flx_cv(ul,duldx)
  rvflx(1:neqns,1:nr) = visc_flx_cv(ur,durdx)
  !
  ! Use the LDG method to compute the common viscous flux vector
  !
  ! Start by applying the averaging operator {{.}} to the viscous flux vectors
  !
  do m = 1,neqns
    do l = 1,nr
      vflux(l,m) = half * (lvflx(m,l) + rvflx(m,l))
    end do
  end do
  !
  ! Apply the differencing (jump) operator [[.]] to the interface viscous
  ! flux vectors and then add it to the common viscous flux vector
  !
  do m = 1,neqns
    !
    vjump = dot_product( unrm(1:nr) , lvflx(m,1:nr)-rvflx(m,1:nr) )
    !
    do l = 1,nr
      vflux(l,m) = vflux(l,m) + LDG_beta(l)*vjump
    end do
    !
  end do
  !
  ! Finally, add the penalty term, (simply the differencing (jump) operator
  ! [[.]] applied to the interface solutions) to the common viscous flux vector
  !
  do m = 1,neqns
    do l = 1,nr
      vflux(l,m) = vflux(l,m) + LDG_tau * (ul(m)*unrm(l) - ur(m)*unrm(l))
    end do
  end do
  !
  ! Finally, dot the common viscous flux vector with the face unit normal
  !
  do m = 1,neqns
    return_value(m) = dot_product( vflux(1:nr,m) , unrm(1:nr) )
  end do
  !
end function common_visc_flx_cv
!
!###############################################################################
!
pure function grad_cv_to_grad_pv(rusp,dcdx,dcdy) result(return_value)
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use geovar, only : nr,neqns
  use ovar,   only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rusp
  real(wp), dimension(:,:), intent(in) :: dcdx
  real(wp), dimension(:,:), intent(in) :: dcdy
  !
  !.. Function Return Value ..
  real(wp), dimension(1:nr,1:neqns,1:size(rusp,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: k,m,n,np
  real(wp) :: cv_inv,rho_inv,rho_inv2
  !
  !.. Local Arrays ..
  real(wp), dimension(1:neqns,1:size(rusp,dim=2)) :: dpdx,dpdy
  !
continue
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  ! Get the number of solution points for the current cell
  !
  np = size(rusp,dim=2)
  !
  ! Compute the primitive variable gradients for the current cell
  !
  do k = 1,np
    !
    rho_inv = one/rusp(1,k)
    rho_inv2 = rho_inv * rho_inv
    !
    ! Copy the density gradient
    !
    dpdx(1,k) = dcdx(1,k)
    dpdy(1,k) = dcdy(1,k)
    !
    ! Using the chain rule:
    !
    ! 1. Compute the velocity gradients (m = 2,neqns-1)
    ! 2. Compute the total energy gradient and store in the
    !    location for the temperature gradient (m = neqns)
    !
    do m = 2,neqns
      dpdx(m,k) = rho_inv*dcdx(m,k) - rho_inv2*rusp(m,k)*dcdx(1,k)
      dpdy(m,k) = rho_inv*dcdy(m,k) - rho_inv2*rusp(m,k)*dcdy(1,k)
    end do
    !
    ! The location for the temperature gradient currently contains
    ! the total energy gradient from the previous loop.
    ! Now use the velocity gradients to subtract the gradient of
    ! kinetic-energy-per-unit-mass from the total energy gradient.
    ! 
    do m = 2,neqns-1
      dpdx(neqns,k) = dpdx(neqns,k) - rho_inv*rusp(m,k)*dpdx(m,k)
      dpdy(neqns,k) = dpdy(neqns,k) - rho_inv*rusp(m,k)*dpdy(m,k)
    end do
    !
    ! Finish computing the temperature gradient by dividing through by
    ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
    !
    dpdx(neqns,k) = dpdx(neqns,k) * cv_inv
    dpdy(neqns,k) = dpdy(neqns,k) * cv_inv
    !
  end do
  !
  ! Put primitive variable gradients in the dusp array
  !
  do k = 1,np
    do m = 1,neqns
      return_value(1,m,k) = dpdx(m,k)
      return_value(2,m,k) = dpdy(m,k)
    end do
  end do
  !
end function grad_cv_to_grad_pv
!
!###############################################################################
!
pure function grad_cv_to_grad_pv_sp(rusp,dcvdx) result(return_value)
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use geovar, only : nr,neqns
  use ovar,   only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:),   intent(in) :: rusp
  real(wp), dimension(:,:), intent(in) :: dcvdx
  !
  !.. Function Return Value ..
  real(wp), dimension(1:nr,1:neqns) :: return_value
  !
  !.. Local Scalars ..
  integer  :: l,m
  real(wp) :: cv_inv,rho_inv,rho_inv2
  !
continue
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  rho_inv = one/rusp(1)
  rho_inv2 = rho_inv * rho_inv
  !
  ! Copy the density gradient
  !
  return_value(1:nr,1) = dcvdx(1:nr,1)
  !
  ! Using the chain rule:
  !
  ! 1. Compute the velocity gradients (m = 2,neqns-1)
  ! 2. Compute the total energy gradient and store in the
  !    location for the temperature gradient (m = neqns)
  !
  do m = 2,neqns
    do l = 1,nr
      return_value(l,m) = rho_inv*dcvdx(l,m) - rho_inv2*rusp(m)*dcvdx(l,1)
    end do
  end do
  !
  ! The location for the temperature gradient currently contains
  ! the total energy gradient from the previous loop.
  ! Now use the velocity gradients to subtract the gradient of
  ! kinetic-energy-per-unit-mass from the total energy gradient.
  ! 
  do m = 2,neqns-1
    do l = 1,nr
      return_value(l,neqns) = return_value(l,neqns) &
                            - rho_inv*rusp(m)*return_value(l,m)
    end do
  end do
  !
  ! Finish computing the temperature gradient by dividing through by
  ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
  !
  return_value(1:nr,neqns) = return_value(1:nr,neqns) * cv_inv
  !
end function grad_cv_to_grad_pv_sp
!
!###############################################################################
!###############################################################################
!###############################################################################
!
! Procedures for getting the conservative variable gradients
!
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine discontinuous_gradient_at_solpts()
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use geovar,  only : neqns,ncell,n_sp2dt,n_sp2dq
  use geovar,  only : pts_in_cell_ptr
  use metvar,  only : dr_mat,ds_mat,dx,dy,jacobian
  use flowvar, only : usp,dusp
  !
  !.. Local Scalars ..
  integer  :: k,l,m,n,n1,n2,np,nw
  real(wp) :: drdx,drdy,dsdx,dsdy
  real(wp) :: dudr,duds,jac_inv
  !
  !.. Local Arrays ..
  real(wp), dimension(1:n_sp2dq,1:neqns) :: tusp
  !
continue
  !
  ! Loop over the cells and get the derivative of
  ! the primitive variables at each solution point
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for solpts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for solpts in cell n
    np = n2-n1+1 ! number of solution points in this cell
    !
    nw = merge(1,2,np == n_sp2dt) ! check if cell n is tria or quad
    !
    ! Store the transpose of usp for this cell to reduce cache misses
    !
    do k = 1,np
      do m = 1,neqns
        tusp(k,m) = usp(m,k+n1-1)
      end do
    end do
    !
    ! Compute the gradient of the conservative variables at each solution point
    !
    do k = 1,np
      !
      l = k+n1-1 ! total index of solution point on current partition
      !
      jac_inv = one/jacobian(l) ! inverse of det(J)
      !
      ! Compute the contravariant metrics
      !
      drdx =  jac_inv * dy(l)%ds
      drdy = -jac_inv * dx(l)%ds
      dsdx = -jac_inv * dy(l)%dr
      dsdy =  jac_inv * dx(l)%dr
      !
      do m = 1,neqns
        !
        dudr = dot_product( dr_mat(1:np,k,nw) , tusp(1:np,m) )
        duds = dot_product( ds_mat(1:np,k,nw) , tusp(1:np,m) )
        !
        dusp(1,m,l) = drdx*dudr + dsdx*duds
        dusp(2,m,l) = drdy*dudr + dsdy*duds
        !
      end do
      !
    end do
    !
  end do
  !
end subroutine discontinuous_gradient_at_solpts
!
!###############################################################################
!
subroutine correct_gradients_at_solpts()
  !
  !.. Use Statements ..
  use geovar,  only : neqns,nfbnd,nface,n_sp1de,n_sp2dt,nr
  use geovar,  only : cells_on_face,pts_in_cell_ptr
  use metvar,  only : cormat,geofa
  use flowvar, only : dusp,faceusp
  !
  use parallel_mod, only : wait_faceusp
  !
  !.. Local Scalars ..
  integer  :: k,l,m,n,nf,r1,r2,l1,l2
  integer  :: jl,jr,wl,wr,nsp
  integer  :: kbl,kbr,kl,kr
  real(wp) :: cbias,hfpb,hfmb
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)              :: unrm
  real(wp), dimension(1:neqns)           :: ljfp,rjfp
  real(wp), dimension(1:neqns,1:n_sp1de) :: lufp,rufp
  !
continue
  !
  nsp = n_sp1de
  !
  ! Make sure that the non-blocking exchange of faceusp has completed
  !
  if (ncpu > 1) call wait_faceusp(faceusp)
  !
  ! Loop over all the faces and add correct the gradient of the
  ! continuous solution polynomial in the cells adjacent to the face
  !
  do nf = 1,nface
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kbl = cells_on_face(5,nf) ! beginning flux point for left  cell
    kbr = cells_on_face(6,nf) ! ending    flux point for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    wl = merge( 1 , 2 , l2-l1+1 == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , r2-r1+1 == n_sp2dt ) ! check if jr is a tria or quad
    !
    ! For all the flux points on the current face, extract the conservative
    ! conservative variables for the left and right cells from faceusp
    !
    lufp(1:neqns,1:nsp) = faceusp(1:neqns,1:nsp,1,nf)
    rufp(1:neqns,1:nsp) = faceusp(1:neqns,1:nsp,2,nf)
    !
    ! Correct the gradients of the conservative variables using the
    ! solution jumps at the interface and the correction function
    !
    ! Use the LDG flux formulation for computing the common interface solutions
    !   u_com = half*(u_left + u_right) - &
    !                (u_left - u_right)*dot_product(LDG_beta,unit_normal)
    !
    do k = 1,nsp
      !
      ! Get the index for the correction function 
      ! corresponding to the current interface flux point
      !
      kl = kbl-1 + k ! index for left  cell
      kr = kbr+1 - k ! index for right cell
      !
      ! Get the unit normal at this flux point
      !
      unrm(1:nr) = unit_vector( geofa(1:nr,k,nf) )
      !
      ! Compute the LDG directional parameters
      !
      cbias = dot_product( LDG_beta(1:nr) , unrm(1:nr) )
      !
      hfpb = half + cbias
      hfmb = half - cbias
      !
      ! Compute the interface jumps for the left and right solutions
      ! NOTE: The common interface solution simplifies to
      !         u_com = (half - cbias)*u_left + (half + cbias)*u_right
      !       Using jump_left = u_com - u_left and jump_right = u_com - u_right
      !         jump_left  = (half + cbias)*(u_right-u_left)
      !         jump_right = (half - cbias)*(u_left-u_right)
      !
      do m = 1,neqns
        ljfp(m) = hfpb*(rufp(m,k) - lufp(m,k))
        rjfp(m) = hfmb*(lufp(m,k) - rufp(m,k))
      end do
      !
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the left cell
      !
      do n = l1,l2
        do m = 1,neqns
          do l = 1,nr
            dusp(l,m,n) = dusp(l,m,n) + unrm(l)*ljfp(m)*cormat(n-l1+1,kl,wl)
          end do
        end do
      end do
      !
      ! Now use the correction function with the interface jump to correct
      ! the gradients of the conservative variables in the right cell
      ! NOTE: only do this for interior cells, i.e. skip this for ghost cells
      !
      if (nf > nfbnd) then
        do n = r1,r2
          do m = 1,neqns
            do l = 1,nr
              dusp(l,m,n) = dusp(l,m,n) + unrm(l)*rjfp(m)*cormat(n-r1+1,kr,wr)
            end do
          end do
        end do
      end if
      !
    end do
    !
  end do
  !
end subroutine correct_gradients_at_solpts
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine correctd
  !
  !... this subroutine corrects the derivatives
  !... ftdxi_(m,k,j) and gtdet_(m,k,j)
  !
  !.. Use Statements ..
  use geovar,  only : neqns,nface,nfbnd,n_sp1de,n_sp2dt,n_sp2dq, &
                      cells_on_face, &
                      pts_on_face, &
                      pts_in_cell_ptr,write_variable
  use metvar,  only : cormat,geofa
  use flowvar, only : residual,vsp
  use generic_mod, only : vsp2u_sp
  !
  !.. Local Scalars ..
  integer   :: nf,k,i,jl,jr,kl,kr,l1,l2,r1,r2,wl,wr,m
  integer   :: kcorbl,kcorbr,kcorl,kcorr
  real(wp) :: ds,dsl,dsr
  !
  !.. Local Arrays ..
  real(wp), dimension(1:2)  :: unrm
  real(wp), dimension(1:neqns) :: vl,vr
  real(wp), dimension(1:neqns) :: fnl,fnr,fjl,fjr,fnc
  !
continue
  !
  ! Loop over all the faces and correct the fluxes at the solution points
  !
  do nf = 1,nface
    !
    jl = cells_on_face(1,nf) ! left  cell on face
    jr = cells_on_face(2,nf) ! right cell on face
    !
    kcorbl = cells_on_face(5,nf) ! starting face sol-pt for left  cell
    kcorbr = cells_on_face(6,nf) ! starting face sol-pt for right cell
    !
    l1 = pts_in_cell_ptr(jl)+1 ! beginning index for sol-pts in left cell
    l2 = pts_in_cell_ptr(jl+1) ! ending    index for sol-pts in left cell
    !
    r1 = pts_in_cell_ptr(jr)+1 ! beginning index for sol-pts in right cell
    r2 = pts_in_cell_ptr(jr+1) ! ending    index for sol-pts in right cell
    !
    wl = merge( 1 , 2 , l2-l1+1 == n_sp2dt ) ! check if jl is a tria or quad
    wr = merge( 1 , 2 , r2-r1+1 == n_sp2dt ) ! check if jr is a tria or quad
    !
    ! Loop through the solution points on the current face
    !
    do k = 1,n_sp1de
      !
      ! Get the unit normal at this flux point
      !
      unrm(:) = unit_vector( geofa(:,k,nf) )
      !
      ! Get the face length for the left and right cell
      !
      ds  = rnorm2( geofa(:,k,nf) )
      dsl = ds * cl(wl)
      dsr = ds * cl(wr)
      !
      kl = pts_on_face(1,k,nf) ! current face sol-pt of left  cell
      kr = pts_on_face(2,k,nf) ! current face sol-pt of right cell
      !
      vl = vsp(:,kl) ! primitive variables at kl
      vr = vsp(:,kr) ! primitive variables at kr
      !
      ! Get the numerical flux at the current face sol-pt
      !
      if (flux_method == ROEFLUX) then
        fnc(:) = upwind_flux_cv ( unrm(:) , vsp2u_sp(vl) , vsp2u_sp(vr) )
      else if (flux_method == HLLC) then
        fnc(:) = hllc_upwind_flux_cv ( unrm(:) , vsp2u_sp(vl) , vsp2u_sp(vr) )
      else if (flux_method == LDFSS) then
        fnc(:) = ldfss_upwind_flux_cv ( unrm(:) , vsp2u_sp(vl) , vsp2u_sp(vr) )
      end if
      !
      ! The flux at the sol-pts kl and kr in the direction of the face normal
      !
      fnl(:) = f_invs_flx_pv(vl)*unrm(1) + g_invs_flx_pv(vl)*unrm(2)
      fnr(:) = f_invs_flx_pv(vr)*unrm(1) + g_invs_flx_pv(vr)*unrm(2)
      !
      ! The jump in the flux at the sol-pts kl and kr
      !
      fjl(:) = fnc(:) - fnl(:)
      fjr(:) = fnc(:) - fnr(:)
      !
      kcorl = kcorbl-1 + k ! index for the correction function at sol-pt kl
      kcorr = kcorbr+1 - k ! index for the correction function at sol-pt kr
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jl
      !
      do i = l1,l2
        do m = 1,neqns
          residual(m,i) = residual(m,i) + dsl*fjl(m)*cormat(i-l1+1,kcorl,wl)
        end do
      end do
      !
      ! Add the flux correction to the residuals of the sol-pts in cell jr
      !
      if (nf > nfbnd) then
        do i = r1,r2
          do m = 1,neqns
            residual(m,i) = residual(m,i) - dsr*fjr(m)*cormat(i-r1+1,kcorr,wr)
          end do
        end do
      end if
      !
    end do
    !
  end do
  !
end subroutine correctd
!
!###############################################################################
!
end module flux_mod
