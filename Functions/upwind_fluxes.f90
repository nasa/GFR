pure function upwind_flux_cv(unrm,cvl,cvr) result(return_value)
  !
  !... input unit normal unx,uny, edge length ds, rl,ul,vl,pl, rr,ur,vr,pr
  !... output: f1,f2,f3,f4 : fluxes.
  !... Roe's flux with Osher's simple wave model (Huynh, SIAM J. Num. 1995)
  !
  use geovar, only : nr
  use ovar, only : gam
  !
  real(wp), dimension(1:nr), intent(in) :: unrm
  real(wp), dimension(1:nq), intent(in) :: cvl
  real(wp), dimension(1:nq), intent(in) :: cvr
  !
  real(wp), dimension(1:nq) :: return_value
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
  use geovar, only : nr
  use ovar, only : gam
  !
  real(wp), dimension(1:nr), intent(in) :: unrm
  real(wp), dimension(1:nq), intent(in) :: pvl
  real(wp), dimension(1:nq), intent(in) :: pvr
  !
  real(wp), dimension(1:nq) :: return_value
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
  use geovar, only : nr
  use ovar, only : gam
  !
  real(wp), dimension(1:nr), intent(in) :: unrm
  real(wp), dimension(1:nq), intent(in) :: cvl
  real(wp), dimension(1:nq), intent(in) :: cvr
  !
  !.. Function Return Value ..
  real(wp), dimension(1:nq) :: return_value
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
  use geovar, only : nr
  use ovar, only : gam
  !
  real(wp), dimension(1:nr), intent(in) :: unrm
  real(wp), dimension(1:nq), intent(in) :: cvl
  real(wp), dimension(1:nq), intent(in) :: cvr
  !
  real(wp), dimension(1:nq) :: return_value
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
  use geovar, only : nr
  use ovar, only : gam
  !
  real(wp), dimension(1:nr), intent(in) :: unrm
  real(wp), dimension(1:nq), intent(in) :: pvl
  real(wp), dimension(1:nq), intent(in) :: pvr
  !
  real(wp), dimension(1:nq) :: return_value
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
