pure subroutine roe_upwind_flux_gpu_3d(n1,n2,facenrm,face_cv,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gp1,ds
  real(wp) :: rl,presl,hl,q2l,ql,rql,betal
  real(wp) :: rr,presr,hr,q2r,qr,rqr,betar
  real(wp) :: rl_inv,vxl,vyl,vzl,el
  real(wp) :: rr_inv,vxr,vyr,vzr,er
  real(wp) :: vxtd,vytd,vztd
  real(wp) :: unx,uny,unz
  real(wp) :: rnx,rny,rnz
  real(wp) :: dw1,hsr1,eta1,udw1
  real(wp) :: dw4,hsr4,eta4,udw4
  real(wp) :: rtd,htd,atd2,atd,q2td,qtd
  real(wp) :: atd2nz,atd2rt,atdnz,atdrt
  !
continue
  !
  gm1 = gam - one
  gp1 = gam + one
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(facenrm,face_cv,return_value)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt(rnx*rnx + rny*rny + rnz*rnz)
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    rl = face_cv(n,nec,1)
    rl_inv = one/rl
    vxl = face_cv(n,nmx,1)*rl_inv
    vyl = face_cv(n,nmy,1)*rl_inv
    vzl = face_cv(n,nmz,1)*rl_inv
    el = face_cv(n,nee,1)
    q2l = vxl*vxl + vyl*vyl + vzl*vzl
    presl = gm1*(el - half*rl*q2l)
    !
    rr = face_cv(n,nec,2)
    rr_inv = one/rr
    vxr = face_cv(n,nmx,2)*rr_inv
    vyr = face_cv(n,nmy,2)*rr_inv
    vzr = face_cv(n,nmz,2)*rr_inv
    er = face_cv(n,nee,2)
    q2r = vxr*vxr + vyr*vyr + vzr*vzr
    presr = gm1*(er - half*rr*q2r)
    !
    hl = half*q2l + gam/gm1 * presl/rl
    hr = half*q2r + gam/gm1 * presr/rr
    ql = vxl*unx + vyl*uny + vzl*unz
    qr = vxr*unx + vyr*uny + vzr*unz
    !
    !... square root averaging
    rtd   = sqrt(rr*rl)
    betal = rl / (rl+rtd)
    betar = one - betal
    vxtd  = betal*vxl + betar*vxr
    vytd  = betal*vyl + betar*vyr
    vztd  = betal*vzl + betar*vzr
    htd   = betal*hl + betar*hr
    q2td  = vxtd*vxtd + vytd*vytd + vztd*vztd
    atd2  = max(gm1*(htd - half*q2td),zero)
    atd   = sqrt(atd2)
    qtd   = vxtd*unx + vytd*uny + vztd*unz
    !
    atd2nz = max(atd2,rndoff)
    atdnz = max(atd,rndoff)
    atd2rt = atd2 / atd2nz
    atdrt = atd / atdnz
    !
    if (qtd >= zero) then
      !
      dw1  = half*((presr-presl)/atd2nz*atd2rt - (qr-ql)*rtd/atdnz*atdrt)
      hsr1 = -one4*gp1*dw1*atd/rtd
      eta1 = max(-abs(qtd-atd) + hsr1, zero)
      udw1 = dw1*(min(qtd-atd,zero) - half*eta1)
      rql  = rl*ql
      !
      return_value(n,nec) = ds*(    rql + udw1                           )
      return_value(n,nmx) = ds*(vxl*rql + udw1*(vxtd-atd*unx) + presl*unx)
      return_value(n,nmy) = ds*(vyl*rql + udw1*(vytd-atd*uny) + presl*uny)
      return_value(n,nmz) = ds*(vzl*rql + udw1*(vztd-atd*unz) + presl*unz)
      return_value(n,nee) = ds*( hl*rql + udw1*(  htd-atd*qtd)           )
      !
    else
      !
      dw4  = half*((presr-presl)/atd2nz*atd2rt + (qr-ql)*rtd/atdnz*atdrt)
      hsr4 = one4*gp1*dw4*atd/rtd
      eta4 = max(-abs(qtd+atd) + hsr4, zero)
      udw4 = dw4*(max(qtd+atd,zero) + half*eta4)
      rqr  = rr*qr
      !
      return_value(n,nec) = ds*(    rqr - udw4                           )
      return_value(n,nmx) = ds*(vxr*rqr - udw4*(vxtd+atd*unx) + presr*unx)
      return_value(n,nmy) = ds*(vyr*rqr - udw4*(vytd+atd*uny) + presr*uny)
      return_value(n,nmz) = ds*(vzr*rqr - udw4*(vztd+atd*unz) + presr*unz)
      return_value(n,nee) = ds*( hr*rqr - udw4*(  htd+atd*qtd)           )
      !
    end if
    !
  end do
  !
end subroutine roe_upwind_flux_gpu_3d
!
!###############################################################################
!
pure subroutine roe_upwind_flux_gpu_2d(n1,n2,facenrm,face_cv,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gp1,ds
  real(wp) :: rl,presl,hl,q2l,ql,rql,betal
  real(wp) :: rr,presr,hr,q2r,qr,rqr,betar
  real(wp) :: rl_inv,vxl,vyl,el
  real(wp) :: rr_inv,vxr,vyr,er
  real(wp) :: vxtd,vytd
  real(wp) :: unx,uny
  real(wp) :: rnx,rny
  real(wp) :: dw1,hsr1,eta1,udw1
  real(wp) :: dw4,hsr4,eta4,udw4
  real(wp) :: rtd,htd,atd2,atd,q2td,qtd
  real(wp) :: atd2nz,atd2rt,atdnz,atdrt
  !
continue
  !
  gm1 = gam - one
  gp1 = gam + one
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(facenrm,face_cv,return_value)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt(rnx*rnx + rny*rny)
    unx = rnx/ds
    uny = rny/ds
    !
    rl = face_cv(n,nec,1)
    rl_inv = one/rl
    vxl = face_cv(n,nmx,1)*rl_inv
    vyl = face_cv(n,nmy,1)*rl_inv
    el = face_cv(n,nee,1)
    q2l = vxl*vxl + vyl*vyl
    presl = gm1*(el - half*rl*q2l)
    !
    rr = face_cv(n,nec,2)
    rr_inv = one/rr
    vxr = face_cv(n,nmx,2)*rr_inv
    vyr = face_cv(n,nmy,2)*rr_inv
    er = face_cv(n,nee,2)
    q2r = vxr*vxr + vyr*vyr
    presr = gm1*(er - half*rr*q2r)
    !
    hl = half*q2l + gam/gm1 * presl/rl
    hr = half*q2r + gam/gm1 * presr/rr
    ql = vxl*unx + vyl*uny
    qr = vxr*unx + vyr*uny
    !
    !... square root averaging
    rtd   = sqrt(rr*rl)
    betal = rl / (rl+rtd)
    betar = one - betal
    vxtd  = betal*vxl + betar*vxr
    vytd  = betal*vyl + betar*vyr
    htd   = betal*hl + betar*hr
    q2td  = vxtd*vxtd + vytd*vytd
    atd2  = max(gm1*(htd - half*q2td),zero)
    atd   = sqrt(atd2)
    qtd   = vxtd*unx + vytd*uny
    !
    atd2nz = max(atd2,rndoff)
    atdnz = max(atd,rndoff)
    atd2rt = atd2 / atd2nz
    atdrt = atd / atdnz
    !
    if (qtd >= zero) then
      !
      dw1  = half*((presr-presl)/atd2nz*atd2rt - (qr-ql)*rtd/atdnz*atdrt)
      hsr1 = -one4*gp1*dw1*atd/rtd
      eta1 = max(-abs(qtd-atd) + hsr1, zero)
      udw1 = dw1*(min(qtd-atd,zero) - half*eta1)
      rql  = rl*ql
      !
      return_value(n,nec) = ds*(    rql + udw1                           )
      return_value(n,nmx) = ds*(vxl*rql + udw1*(vxtd-atd*unx) + presl*unx)
      return_value(n,nmy) = ds*(vyl*rql + udw1*(vytd-atd*uny) + presl*uny)
      return_value(n,nee) = ds*( hl*rql + udw1*(  htd-atd*qtd)           )
      !
    else
      !
      dw4  = half*((presr-presl)/atd2nz*atd2rt + (qr-ql)*rtd/atdnz*atdrt)
      hsr4 = one4*gp1*dw4*atd/rtd
      eta4 = max(-abs(qtd+atd) + hsr4, zero)
      udw4 = dw4*(max(qtd+atd,zero) + half*eta4)
      rqr  = rr*qr
      !
      return_value(n,nec) = ds*(    rqr - udw4                           )
      return_value(n,nmx) = ds*(vxr*rqr - udw4*(vxtd+atd*unx) + presr*unx)
      return_value(n,nmy) = ds*(vyr*rqr - udw4*(vytd+atd*uny) + presr*uny)
      return_value(n,nee) = ds*( hr*rqr - udw4*(  htd+atd*qtd)           )
      !
    end if
    !
  end do
  !
end subroutine roe_upwind_flux_gpu_2d
!
!###############################################################################
!
pure subroutine roe_upwind_flux_gpu_nd(n1,n2,facenrm,face_cv,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: gm1,gp1,ds
  real(wp) :: rl,presl,hl,q2l,ql,rql,betal
  real(wp) :: rr,presr,hr,q2r,qr,rqr,betar
  real(wp) :: dw1,hsr1,eta1,udw1
  real(wp) :: dw4,hsr4,eta4,udw4
  real(wp) :: rtd,htd,atd2,atd,qtd
  real(wp) :: atd2nz,atd2rt,atdnz,atdrt
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(facenrm,dim=1)) :: unrm,rnrm
  real(wp), dimension(1:size(facenrm,dim=1)) :: vell,velr,veltd
  real(wp), dimension(1:nq) :: cvl,cvr,flx
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  gm1 = gam - one
  gp1 = gam + one
  !
  !$ACC PARALLEL LOOP &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRIVATE(unrm,rnrm,vell,velr,veltd,cvl,cvr,flx) &
  !$ACC&              PRESENT(facenrm,face_cv,return_value)
  do n = n1,n2
    !
    rnrm(1:nr) = facenrm(1:nr,n)
    ds = norm2(rnrm(1:nr))
    unrm(1:nr) = rnrm(1:nr) / ds
    cvl(1:nq) = face_cv(n,1:nq,1)
    cvr(1:nq) = face_cv(n,1:nq,2)
    !
    rl = cvl(nec)
    vell = cvl(nmb:nme)/cvl(nec)
    q2l = selfdot(vell)
    presl = gm1*(cvl(nee) - half*rl*q2l)
    !
    rr = cvr(nec)
    velr = cvr(nmb:nme)/cvr(nec)
    q2r = selfdot(velr)
    presr = gm1*(cvr(nee) - half*rr*q2r)
    !
    hl = half*q2l + gam/gm1 * presl/rl
    hr = half*q2r + gam/gm1 * presr/rr
    ql = dot_product(vell,unrm)
    qr = dot_product(velr,unrm)
    !
    !... square root averaging
    rtd   = sqrt(rr*rl)
    betal = rl / (rl+rtd)
    betar = one - betal
    veltd = betal*vell + betar*velr
    htd   = betal*hl + betar*hr
    atd2  = max(gm1*(htd - half*selfdot(veltd)),zero)
    atd   = sqrt(atd2)
    qtd   = dot_product(veltd,unrm)
    !
    atd2nz = max(atd2,rndoff)
    atdnz = max(atd,rndoff)
    atd2rt = atd2 / atd2nz
    atdrt = atd / atdnz
    !
    if (qtd >= zero) then
      !
      dw1  = half*((presr-presl)/atd2nz*atd2rt - (qr-ql)*rtd/atdnz*atdrt)
      hsr1 = -one4*gp1*dw1*atd/rtd
      eta1 = max(-abs(qtd-atd) + hsr1, zero)
      udw1 = dw1*(min(qtd-atd,zero) - half*eta1)
      rql  = rl*ql
      !
      flx(nec)     =      rql + udw1
      flx(nmb:nme) = vell*rql + udw1*(veltd-atd*unrm) + presl*unrm
      flx(nee)     =   hl*rql + udw1*(  htd-atd*qtd)
      !
    else
      !
      dw4  = half*((presr-presl)/atd2nz*atd2rt + (qr-ql)*rtd/atdnz*atdrt)
      hsr4 = one4*gp1*dw4*atd/rtd
      eta4 = max(-abs(qtd+atd) + hsr4, zero)
      udw4 = dw4*(max(qtd+atd,zero) + half*eta4)
      rqr  = rr*qr
      !
      flx(nec)     =      rqr - udw4
      flx(nmb:nme) = velr*rqr - udw4*(veltd+atd*unrm) + presr*unrm
      flx(nee)     =   hr*rqr - udw4*(  htd+atd*qtd)
      !
    end if
    !
    return_value(n,1:nq) = flx(1:nq)*ds
    !
  end do
  !
end subroutine roe_upwind_flux_gpu_nd
!
!###############################################################################
!
pure subroutine hllc_upwind_flux_gpu(n1,n2,facenrm,face_cv,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nfp,nr
  real(wp) :: ds
  !
  !.. Local Arrays ..
  real(wp) :: rnrm(1:size(facenrm,dim=1))
  real(wp) :: unrm(1:size(facenrm,dim=1))
  real(wp) :: cvl(1:nq)
  real(wp) :: cvr(1:nq)
  real(wp) :: flx(1:nq)
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  !$ACC PARALLEL LOOP &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRIVATE(unrm,rnrm,cvl,cvr,flx) &
  !$ACC&              PRESENT(facenrm,face_cv,return_value)
  do n = n1,n2
    rnrm(1:nr) = facenrm(1:nr,n)
    ds = norm2(rnrm(1:nr))
    unrm(1:nr) = rnrm(1:nr) / ds
    cvl(1:nq) = face_cv(n,1:nq,1)
    cvr(1:nq) = face_cv(n,1:nq,2)
    call hllc_upwind_flux_sub(unrm,cvl,cvr,flx)
    return_value(n,1:nq) = flx(1:nq)*ds
  end do
  !
end subroutine hllc_upwind_flux_gpu
!
!###############################################################################
!
pure subroutine ldfss_upwind_flux_gpu(n1,n2,facenrm,face_cv,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nfp,nr
  real(wp) :: ds
  !
  !.. Local Arrays ..
  real(wp) :: rnrm(1:size(facenrm,dim=1))
  real(wp) :: unrm(1:size(facenrm,dim=1))
  real(wp) :: cvl(1:nq)
  real(wp) :: cvr(1:nq)
  real(wp) :: flx(1:nq)
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  !$ACC PARALLEL LOOP &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRIVATE(unrm,rnrm,cvl,cvr,flx) &
  !$ACC&              PRESENT(facenrm,face_cv,return_value)
  do n = n1,n2
    rnrm(1:nr) = facenrm(1:nr,n)
    ds = norm2(rnrm(1:nr))
    unrm(1:nr) = rnrm(1:nr) / ds
    cvl(1:nq) = face_cv(n,1:nq,1)
    cvr(1:nq) = face_cv(n,1:nq,2)
    call ldfss_upwind_flux_sub(unrm,cvl,cvr,flx)
    return_value(n,1:nq) = flx(1:nq)*ds
  end do
  !
end subroutine ldfss_upwind_flux_gpu
!
!###############################################################################
!
pure subroutine ausm_plus_upwind_flux_gpu(n1,n2,facenrm,face_cv,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nfp,nr
  real(wp) :: ds
  !
  !.. Local Arrays ..
  real(wp) :: rnrm(1:size(facenrm,dim=1))
  real(wp) :: unrm(1:size(facenrm,dim=1))
  real(wp) :: cvl(1:nq)
  real(wp) :: cvr(1:nq)
  real(wp) :: flx(1:nq)
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  !$ACC PARALLEL LOOP &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRIVATE(unrm,rnrm,cvl,cvr,flx) &
  !$ACC&              PRESENT(facenrm,face_cv,return_value)
  do n = n1,n2
    rnrm(1:nr) = facenrm(1:nr,n)
    ds = norm2(rnrm(1:nr))
    unrm(1:nr) = rnrm(1:nr) / ds
    cvl(1:nq) = face_cv(n,1:nq,1)
    cvr(1:nq) = face_cv(n,1:nq,2)
    call ausm_plus_upwind_flux_sub(unrm,cvl,cvr,flx)
    return_value(n,1:nq) = flx(1:nq)*ds
  end do
  !
end subroutine ausm_plus_upwind_flux_gpu
!
!###############################################################################
!
pure subroutine rotated_rhll_upwind_flux_gpu(n1,n2,facenrm,face_cv,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nfp,nr
  real(wp) :: ds
  !
  !.. Local Arrays ..
  real(wp) :: rnrm(1:size(facenrm,dim=1))
  real(wp) :: unrm(1:size(facenrm,dim=1))
  real(wp) :: cvl(1:nq)
  real(wp) :: cvr(1:nq)
  real(wp) :: flx(1:nq)
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  if (nr == 3) then
    !
    !$ACC PARALLEL LOOP &
    !$ACC&              VECTOR_LENGTH(32) &
    !$ACC&              PRIVATE(unrm,rnrm,cvl,cvr,flx) &
    !$ACC&              PRESENT(facenrm,face_cv,return_value)
    do n = n1,n2
      rnrm(1:3) = facenrm(1:3,n)
      ds = norm2(rnrm(1:3))
      unrm(1:3) = rnrm(1:3) / ds
      cvl(1:nq) = face_cv(n,1:nq,1)
      cvr(1:nq) = face_cv(n,1:nq,2)
      call rotated_rhll_upwind_flux_sub_3d(unrm,cvl,cvr,flx)
      return_value(n,1:nq) = flx(1:nq)*ds
    end do
    !
  else if (nr == 2) then
    !
    !$ACC PARALLEL LOOP &
    !$ACC&              VECTOR_LENGTH(32) &
    !$ACC&              PRIVATE(unrm,rnrm,cvl,cvr,flx) &
    !$ACC&              PRESENT(facenrm,face_cv,return_value)
    do n = n1,n2
      rnrm(1:2) = facenrm(1:2,n)
      ds = norm2(rnrm(1:2))
      unrm(1:2) = rnrm(1:2) / ds
      cvl(1:nq) = face_cv(n,1:nq,1)
      cvr(1:nq) = face_cv(n,1:nq,2)
      call rotated_rhll_upwind_flux_sub_2d(unrm,cvl,cvr,flx)
      return_value(n,1:nq) = flx(1:nq)*ds
    end do
    !
  end if
  !
end subroutine rotated_rhll_upwind_flux_gpu
!
!###############################################################################
!
pure function upwind_flux_cv(unrm,cvl,cvr) result(return_value)
  !
  !... input unit normal unx,uny, edge length ds, rl,ul,vl,presl, rr,ur,vr,presr
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
  real(wp) :: gm1,gp1
  real(wp) :: rl,presl,hl,q2l,ql,rql,betal
  real(wp) :: rr,presr,hr,q2r,qr,rqr,betar
  real(wp) :: dw1,hsr1,eta1,udw1
  real(wp) :: dw4,hsr4,eta4,udw4
  real(wp) :: rtd,htd,atd2,atd,qtd
  real(wp) :: atd2nz,atd2rt,atdnz,atdrt
  !
  real(wp), dimension(1:nr) :: vell,velr,veltd
  !
continue
  !
  gm1 = gam - one
  gp1 = gam + one
  !
  rl = cvl(nec)
  vell = cvl(nmb:nme)/cvl(nec)
  q2l = selfdot(vell)
  presl = gm1*(cvl(nee) - half*rl*q2l)
  !
  rr = cvr(nec)
  velr = cvr(nmb:nme)/cvr(nec)
  q2r = selfdot(velr)
  presr = gm1*(cvr(nee) - half*rr*q2r)
  !
  hl = half*q2l + gam/gm1 * presl/rl
  hr = half*q2r + gam/gm1 * presr/rr
  ql = dot_product(vell,unrm)
  qr = dot_product(velr,unrm)
  !
  !... square root averaging
  rtd   = sqrt(rr*rl)
  betal = rl / (rl+rtd)
  betar = one - betal
  veltd = betal*vell + betar*velr
  htd   = betal*hl + betar*hr
  atd2  = max(gm1*(htd - half*selfdot(veltd)),zero)
  atd   = sqrt(atd2)
  qtd   = dot_product(veltd,unrm)
  !
  atd2nz = max(atd2,rndoff)
  atdnz = max(atd,rndoff)
  atd2rt = atd2 / atd2nz
  atdrt = atd / atdnz
  !
  if (qtd >= zero) then
    !
    dw1  = half*((presr-presl)/atd2nz*atd2rt - (qr-ql)*rtd/atdnz*atdrt)
    hsr1 = -one4*gp1*dw1*atd/rtd
    eta1 = max(-abs(qtd-atd) + hsr1, zero)
    udw1 = dw1*(min(qtd-atd,zero) - half*eta1)
    rql  = rl*ql
    !
    return_value(nec)     =      rql + udw1
    return_value(nmb:nme) = vell*rql + udw1*(veltd-atd*unrm) + presl*unrm
    return_value(nee)     =   hl*rql + udw1*(  htd-atd*qtd)
    !
  else
    !
    dw4  = half*((presr-presl)/atd2nz*atd2rt + (qr-ql)*rtd/atdnz*atdrt)
    hsr4 = one4*gp1*dw4*atd/rtd
    eta4 = max(-abs(qtd+atd) + hsr4, zero)
    udw4 = dw4*(max(qtd+atd,zero) + half*eta4)
    rqr  = rr*qr
    !
    return_value(nec)     =      rqr - udw4
    return_value(nmb:nme) = velr*rqr - udw4*(veltd+atd*unrm) + presr*unrm
    return_value(nee)     =   hr*rqr - udw4*(  htd+atd*qtd)
    !
  end if
  !
end function upwind_flux_cv
!
!###############################################################################
!
pure function upwind_flux_pv(unrm,pvl,pvr) result(return_value)
  !
  !... input unit normal unx,uny, edge length ds, rl,ul,vl,presl, rr,ur,vr,presr
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
  real(wp) :: gm1,gp1
  real(wp) :: rl,presl,hl,ql,rql,betal
  real(wp) :: rr,presr,hr,qr,rqr,betar
  real(wp) :: dw1,hsr1,eta1,udw1
  real(wp) :: dw4,hsr4,eta4,udw4
  real(wp) :: rtd,htd,atd2,atd,qtd
  !
  real(wp), dimension(1:nr) :: vell,velr,veltd
  !
continue
  !
  gm1 = gam - one
  gp1 = gam + one
  !
  rl = pvl(nec)
  vell = pvl(nmb:nme)
  presl = pvl(nee)
  !
  rr = pvr(nec)
  velr = pvr(nmb:nme)
  presr = pvr(nee)
  !
  hl = half*selfdot(vell) + gam/gm1 * presl/rl
  hr = half*selfdot(velr) + gam/gm1 * presr/rr
  ql = dot_product(vell,unrm)
  qr = dot_product(velr,unrm)
  !
  !... square root averaging
  rtd   = sqrt(rr*rl)
  betal = rl / (rl+rtd)
  betar = one - betal
  veltd = betal*vell + betar*velr
  htd   = betal*hl + betar*hr
  atd2  = gm1*(htd - half*selfdot(veltd))
  atd   = sqrt(atd2)
  qtd   = dot_product(veltd,unrm)
  !
  if (qtd >= zero) then
    !
    dw1  = half*((presr-presl)/atd2 - (qr-ql)*rtd/atd)
    hsr1 = -one4*gp1*dw1*atd/rtd
    eta1 = max(-abs(qtd-atd) + hsr1, zero)
    udw1 = dw1*(min(qtd-atd,zero) - half*eta1)
    rql  = rl*ql
    !
    return_value(nec)     =      rql + udw1
    return_value(nmb:nme) = vell*rql + udw1*(veltd-atd*unrm) + presl*unrm
    return_value(nee)     =   hl*rql + udw1*(  htd-atd*qtd)
    !
  else
    !
    dw4  = half*((presr-presl)/atd2 + (qr-ql)*rtd/atd)
    hsr4 = one4*gp1*dw4*atd/rtd
    eta4 = max(-abs(qtd+atd) + hsr4, zero)
    udw4 = dw4*(max(qtd+atd,zero) + half*eta4)
    rqr  = rr*qr
    !
    return_value(nec)     =      rqr - udw4
    return_value(nmb:nme) = velr*rqr - udw4*(veltd+atd*unrm) + presr*unrm
    return_value(nee)     =   hr*rqr - udw4*(  htd+atd*qtd)
    !
  end if
  !
end function upwind_flux_pv
!
!###############################################################################
!
pure function hllc_upwind_flux_cv(unrm,cvl,cvr) result(return_value)
  !
  !... input unit normal unx,uny, edge length ds, rl,ul,vl,presl, rr,ur,vr,presr
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
  real(wp) :: gm1
  real(wp) :: rl,presl,hl,al,q2l,vnl,rhoel,wgl,sl,dsl
  real(wp) :: rr,presr,hr,ar,q2r,vnr,rhoer,wgr,sr,dsr
  real(wp) :: srt,ra,ha,q2a,vna,aa
  real(wp) :: rho,rhou,rhov,rhoe,pstar,sm
  !
  real(wp), dimension(1:nr) :: vell,velr,vela,rvel
  !
continue
  !
  gm1 = gam - one
  !
  rl = cvl(nec)
  vell = cvl(nmb:nme)/cvl(nec)
  q2l = selfdot(vell)
  presl = gm1*(cvl(nee) - half*rl*q2l)
  hl = half*q2l + gam/gm1 * presl/rl
  vnl = dot_product(vell,unrm)
  rhoel = cvl(nee)
  !
  rr = cvr(nec)
  velr = cvr(nmb:nme)/cvr(nec)
  q2r = selfdot(velr)
  presr = gm1*(cvr(nee) - half*rr*q2r)
  hr = half*q2r + gam/gm1 * presr/rr
  vnr = dot_product(velr,unrm)
  rhoer = cvr(nee)
  !
  srt = sqrt(max(rr/rl,rndoff))
  wgl = one/(srt+one)
  wgr = srt/(srt+one)
  !
  ra = sqrt(max(rl*rr,rndoff))
  vela = vell*wgl + velr*wgr
  ha = hl*wgl + hr*wgr
  !
  q2a = selfdot(vela)
  vna = dot_product(vela,unrm)
  !
  aa = sqrt(max(gm1*(ha-half*q2a),rndoff))
  !
  sl = min( vnl-al , vna-aa )
  sr = max( vnr+ar , vna+aa )
  !
  if (sl > zero) then
    !
    return_value(nec)     = vnl*rl
    return_value(nmb:nme) = vnl*rl*vell + presl*unrm
    return_value(nee)     = vnl*rhoel   + presl*vnl
    !
  else if (sr < zero) then
    !
    return_value(nec)     = vnr*rr
    return_value(nmb:nme) = vnr*rr*velr + presr*unrm
    return_value(nee)     = vnr*rhoer   + presr*vnr
    !
  else
    !
    dsl = sl - vnl
    dsr = sr - vnr
    !
    sm = ((rr*vnr*dsr-presr) - (rl*vnl*dsl-presl)) / (rr*dsr-rl*dsl)
    !
    pstar = presl*dsl*(sm-vnl) + presl
    !
    if (sm > zero) then
      !
      rho  = dsl*rl/(sl-sm)
      rvel = (dsl*rl*vell + (pstar-presl)*unrm)/(sl-sm)
      rhoe = (dsl*rhoel - presl*vnl + pstar*sm)/(sl-sm)
      !
    else
      !
      rho  = dsr*rr/(sr-sm)
      rvel = (dsr*rr*velr + (pstar-presr)*unrm)/(sr-sm)
      rhoe = (dsr*rhoer - presr*vnr + pstar*sm)/(sr-sm)
      !
    end if
    !
    return_value(nec)     = sm*rho
    return_value(nmb:nme) = sm*rvel + pstar*unrm
    return_value(nee)     = sm*rhoe + pstar*sm
    !
  end if
  !
end function hllc_upwind_flux_cv
!
!###############################################################################
!
pure function ldfss_upwind_flux_cv(unrm,cvl,cvr) result(return_value)
  !
  !... input unit normal unx,uny, edge length ds, rl,ul,vl,presl, rr,ur,vr,presr
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
  real(wp) :: gm1,ahalf
  real(wp) :: rl,presl,hl,al,q2l
  real(wp) :: rr,presr,hr,ar,q2r
  real(wp) :: alpha_l,beta_l,mach_l,mstar_l
  real(wp) :: alpha_r,beta_r,mach_r,mstar_r
  real(wp) :: machp_half,uc_l,flux_l,pplus_l
  real(wp) :: machm_half,uc_r,flux_r,pplus_r
  real(wp) :: m_half,mstar_half
  real(wp) :: del_p,sum_p,p_net
  !
  real(wp), dimension(1:nr) :: vell,velr
  !
continue
  !
  gm1 = gam - one
  !
  rl = cvl(nec)
  vell = cvl(nmb:nme)/cvl(nec)
  q2l = selfdot(vell)
  presl = gm1*(cvl(nee) - half*rl*q2l)
  hl = half*q2l + gam/gm1 * presl/rl
  al = sqrt(max(gam*presl/rl,rndoff))
  !
  rr = cvr(nec)
  velr = cvr(nmb:nme)/cvr(nec)
  q2r = selfdot(velr)
  presr = gm1*(cvr(nee) - half*rr*q2r)
  hr = half*q2r + gam/gm1 * presr/rr
  ar = sqrt(max(gam*presr/rr,rndoff))
  !
  ahalf = half*(al+ar)
  !
  mach_l = dot_product(vell,unrm)/ahalf
  mach_r = dot_product(velr,unrm)/ahalf
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
  del_p = presl - presr
  sum_p = presl + presr
  !
  machp_half = mstar_half*(one - del_p/sum_p - two*abs(del_p)/presl)
  machm_half = mstar_half*(one + del_p/sum_p - two*abs(del_p)/presr)
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
  p_net = presl*(alpha_l*(one+beta_l) - beta_l*pplus_l) + &
          presr*(alpha_r*(one+beta_r) - beta_r*pplus_r)
  !
  return_value(nec)     = flux_l      + flux_r
  return_value(nmb:nme) = flux_l*vell + flux_r*velr + p_net*unrm
  return_value(nee)     = flux_l*hl   + flux_r*hr
  !
end function ldfss_upwind_flux_cv
!
!###############################################################################
!
pure function ausm_plus_upwind_flux(unrm,pvl,pvr) result(return_value)
  !
  !   Compute Euler fluxes using AUSM+
  !... input unit normal unx,uny, edge length ds, rl,ul,vl,presl, rr,ur,vr,presr
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
  real(wp) :: gm1
  real(wp) :: rl,presl,hl,al,v2l,a2l,h0l,uul,f1l
  real(wp) :: rr,presr,hr,ar,v2r,a2r,h0r,uur,f1r
  real(wp) :: alr,alr2,rm2,rmref2,coef1,coef2,prec_f
  real(wp) :: aml,absml,am1p,aml1,addp,aml2,ampl,ppl
  real(wp) :: amr,absmr,am1m,amr1,addm,amr2,ammr,pmr
  real(wp) :: amhalf,denom,amhalfp,amhalfm,phalf,amco2
  !
  real(wp), dimension(1:nr) :: vell,velr
  !
  !.. Local Parameters ..
  real(wp), parameter :: amco2m = 1.0e-4_wp         
  real(wp), parameter :: alphax = 0.1875_wp
  real(wp), parameter :: betax  = 0.125_wp
  !
continue
  !
  gm1 = gam - one
  !
  rl = pvl(nec)
  vell = pvl(nmb:nme)
  presl = pvl(nee)
  !
  rr = pvr(nec)
  velr = pvr(nmb:nme)
  presr = pvr(nee)
  !
  hl  = gam*presl/(gm1*rl)
  v2l = half * selfdot(vell)
  h0l = hl + v2l
  a2l = gm1*hl
  al  = sqrt(a2l)
  uul = dot_product(vell,unrm)
  !
  hr  = gam*presr/(gm1*rr)
  v2r = half * selfdot(velr)
  h0r = hr + v2r
  a2r = gm1*hr
  ar  = sqrt(a2r)
  uur = dot_product(velr,unrm)
  !
  alr    = half*(al + ar)
  alr2   = alr**2
  rm2    = (v2r+v2l)/alr2
  amco2  = max(amco2m,one4*rm2)
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
    aml2 = one4*aml1**2
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
    amr2 = one4*amr1**2
    ammr = -amr2 - betax*addm
    pmr  = amr2*(two + amr) - alphax*amr*addm
  else
    ammr = am1m
    pmr  = am1m/amr
  end if
  !
  amhalf  = ampl + ammr
  denom   = ( rmref2*(presl+presr)*(rl+rr) - &
                zero*(presl-presr)*(rl-rr) ) / sqrt(rl*rr)
  amhalf  = amhalf + half*(presl-presr)*((ampl-ammr) - (am1p-am1m))/denom
  amhalfp = half*(amhalf + abs(amhalf))
  amhalfm = amhalf - amhalfp
  !
  phalf = ppl*presl + pmr*presr - ppl*pmr*(rl+rr)/two*alr*(uur-uul)
  !
  f1l = rl*amhalfp*alr
  f1r = rr*amhalfm*alr
  !
  return_value(nec)     = f1l      + f1r
  return_value(nmb:nme) = f1l*vell + f1r*velr + unrm*phalf
  return_value(nee)     = f1l*h0l  + f1r*h0r
  !
end function ausm_plus_upwind_flux
!
!###############################################################################
!
pure function rotated_rhll_upwind_flux_cv_2d(unrm,cvl,cvr) result(return_value)
 !
 !.. Use Statements ..
 use geovar, only : nr
 use ovar,   only : gam
 !
 !.. Formal Arguments ..
 real(wp), dimension(1:nr), intent(in) :: unrm ! unit face normal
 real(wp), dimension(1:nq), intent(in) :: cvl
 real(wp), dimension(1:nq), intent(in) :: cvr
 !
 !.. Function Result ..
 real(wp), dimension(1:nq) :: return_value
 !
 real(wp) :: nx, ny
 real(wp) :: nx1, ny1, nx2, ny2      ! Rotated normals, n1 and n2
 real(wp) :: tx, ty                  ! Tangent vector (taken as n1)
 real(wp) :: alpha1, alpha2          ! Projections of the new normals
 real(wp) :: vxL, vxR, vyL, vyR      ! Velocity components.
 real(wp) :: rhoL, rhoR, presL, presR      ! Primitive variables.
 real(wp) :: vnL, vnR, vtL, vtR      ! Normal and tagent velocities
 real(wp) :: aL, aR, HL, HR          ! Speeds of sound and total enthalpy
 real(wp) :: RT,rho,vx,vy,H,a        ! Roe-averages
 real(wp) :: vn, vt                  ! Normal and tagent velocities (Roe-ave)
 real(wp) :: drho,dvx,dvy,dvn,dvt,dp ! Wave strenghs
 real(wp) :: abs_dq                  ! Magnitude of the velocity difference
 real(wp) :: SRp,SLm                 ! Wave speeds for the HLL part
 real(wp) :: temp
 !
 !.. Local Arrays ..
 real(wp) :: dV(1:4)              ! Wave strenghs
 real(wp) :: abs_ws(1:4), ws(1:4) ! Abs wave speeds
 real(wp) :: dws(1:4)             ! width of a parabolic fit for entropy fix
 real(wp) :: Rv(1:nq,1:4)         ! Right-eigevectors
 real(wp) :: fL(1:nq), fR(1:nq)   ! Fluxes
 real(wp) :: diss(1:nq)           ! Dissipation term
 !
continue
  !
  ! Face normal vector (unit vector)
  !
  nx = unrm(1)
  ny = unrm(2)
  !
  ! Primitive and other variables.
  !
  ! Left state
  !
  rhoL = cvl(1)
  vxL = cvl(2)/cvl(1)
  vyL = cvl(3)/cvl(1)
  presL = (gam-one)*( cvl(4) - half*rhoL*(vxL*vxL+vyL*vyL) )
  aL = sqrt(gam*presL/rhoL)
  HL = ( cvl(4) + presL ) / rhoL
  !
  ! Right state
  !
  rhoR = cvr(1)
  vxR = cvr(2)/cvr(1)
  vyR = cvr(3)/cvr(1)
  presR = (gam-one)*( cvr(4) - half*rhoR*(vxR*vxR+vyR*vyR) )
  aR = sqrt(gam*presR/rhoR)
  HR = ( cvr(4) + presR ) / rhoR
  !
  vnL = vxL*nx + vyL*ny
  vnR = vxR*nx + vyR*ny
  !
  ! Compute the flux.
  !
  fL(1) = rhoL*vnL
  fL(2) = rhoL*vnL * vxL + presL*nx
  fL(3) = rhoL*vnL * vyL + presL*ny
  fL(4) = rhoL*vnL *  HL
  !
  fR(1) = rhoR*vnR
  fR(2) = rhoR*vnR * vxR + presR*nx
  fR(3) = rhoR*vnR * vyR + presR*ny
  fR(4) = rhoR*vnR *  HR
  !
  ! Define n1 and n2, and compute alpha1 and alpha2: (4.2) in original paper.
  ! (NB: n1 and n2 may need to be frozen at some point during 
  !      a steady calculation to fully make it converge. For time-accurate 
  !      calculation, this is fine.)
  ! NB: For a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).
  !
  abs_dq = sqrt( (vxR-vxL)**2+(vyR-vyL)**2 )
  !
  if (abs_dq > eps12) then
    nx1 = (vxR-vxL)/abs_dq
    ny1 = (vyR-vyL)/abs_dq
  else
    nx1 = -ny 
    ny1 =  nx
  endif
  !
  alpha1 = nx * nx1 + ny * ny1 
  !
  ! To make alpha1 always positive.
  !
  temp = sign(one,alpha1)
  nx1 = temp * nx1
  ny1 = temp * ny1
  alpha1 = temp * alpha1
  !
  ! Take n2 as perpendicular to n1.
  !
  nx2 = -ny1
  ny2 =  nx1
  alpha2 = nx * nx2 + ny * ny2
  !
  ! To make alpha2 always positive.
  !
  temp = sign(one,alpha2)
  nx2 = temp * nx2
  ny2 = temp * ny2
  alpha2 = temp * alpha2
  !
  ! Now we are going to compute the Roe flux with n2 as the normal
  ! and n1 as the tagent vector, with modified wave speeds (5.12)
  !
  ! Compute the Roe Averages
  !
  RT = sqrt(rhoR/rhoL)
  rho = RT*rhoL
  vx = (vxL+RT*vxR)/(one+RT)
  vy = (vyL+RT*vyR)/(one+RT)
  H = ( HL+RT* HR)/(one+RT)
  a = sqrt( (gam-one)*(H-half*(vx*vx+vy*vy)) )
  vn = vx*nx2+vy*ny2
  vt = vx*nx1+vy*ny1
  !
  ! Wave Strengths (remember that n2 is the normal and n1 is the tangent.)
  !
  vnL = vxL*nx2 + vyL*ny2
  vnR = vxR*nx2 + vyR*ny2
  vtL = vxL*nx1 + vyL*ny1
  vtR = vxR*nx1 + vyR*ny1
  !
  drho = rhoR - rhoL 
  dp =   presR - presL
  dvn =  vnR - vnL
  dvt =  vtR - vtL
  !
  dV(1) = (dp - rho*a*dvn )/(two*a*a)
  dV(2) =  rho*dvt/a
  dV(3) =  drho - dp/(a*a)
  dV(4) = (dp + rho*a*dvn )/(two*a*a)
  !
  ! Wave Speeds for Roe flux part.
  !
  ws(1) = vn-a
  ws(2) = vn
  ws(3) = vn
  ws(4) = vn+a
  abs_ws  = abs(ws)
  !
  ! Harten's Entropy Fix JCP(1983), 49, pp357-393:
  ! only for the nonlinear fields.
  !
  dws(1) = one5
  if ( abs_ws(1) < dws(1) ) then
    abs_ws(1) = half*(abs_ws(1)*abs_ws(1)/dws(1)+dws(1))
  end if
  !
  dws(4) = one5
  if ( abs_ws(4) < dws(4) ) then
    abs_ws(4) = half*(abs_ws(4)*abs_ws(4)/dws(4)+dws(4))
  end if
  !
  ! HLL wave speeds, evaluated with [nx1,ny1] (=tangent wrt n2).
  !
  SRp = max( zero, vtR + aR, vt + a)
  SLm = min( zero, vtL - aL, vt - a)
  !
  ! Modified wave speeds for the Rotated-RHLL flux: (5.12) in original paper.
  !
  ws = alpha2*abs_ws - ( alpha2*(SRp+SLm)*ws + two*alpha1*SRp*SLm )/ (SRp-SLm)
  !
  ! Right Eigenvectors: with n2 as normal and n1 as tangent.
  !
  tx = nx1
  ty = ny1
  !
  Rv(1,1) = one    
  Rv(2,1) = vx - a*nx2
  Rv(3,1) = vy - a*ny2
  Rv(4,1) =  H - vn*a
  !
  Rv(1,2) = zero
  Rv(2,2) = a*tx
  Rv(3,2) = a*ty
  Rv(4,2) = a*vt
  !
  Rv(1,3) = one
  Rv(2,3) = vx
  Rv(3,3) = vy 
  Rv(4,3) = half*(vx*vx+vy*vy)
  !
  Rv(1,4) = one
  Rv(2,4) = vx + a*nx2
  Rv(3,4) = vy + a*ny2
  Rv(4,4) =  H + vn*a
  !
  ! Dissipation Term: Roe dissipation with the modified wave speeds.
  !
  diss(:) = ws(1)*dV(1)*Rv(:,1) + &
            ws(2)*dV(2)*Rv(:,2) + &
            ws(3)*dV(3)*Rv(:,3) + &
            ws(4)*dV(4)*Rv(:,4)
  !
  ! Compute the Rotated-RHLL flux.
  !
  return_value = (SRp*fL - SLm*fR)/(SRp-SLm) - half*diss
  !
end function rotated_rhll_upwind_flux_cv_2d
!
!###############################################################################
!
pure function rotated_rhll_upwind_flux_cv_3d(unrm,cvl,cvr) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr), intent(in) :: unrm ! unit face normal
  real(wp), dimension(1:nq), intent(in) :: cvl
  real(wp), dimension(1:nq), intent(in) :: cvr
  !
  !.. Function Result ..
  real(wp), dimension(1:nq) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: nx, ny, nz                ! Face normal vector
  real(wp) :: uL, uR, vL, vR, wL, wR    ! Velocity components.
  real(wp) :: rhoL, rhoR, presL, presR        ! Primitive variables.
  real(wp) :: qnL, qnR                  ! Normal velocity
  real(wp) :: aL, aR, HL, HR            ! Speed of sound, Total enthalpy
  real(wp) :: RT,rho,u,v,w,H,a,qn       ! Roe-averages
  real(wp) :: drho,dqn,dp               ! Wave strengths
  real(wp) :: du, dv, dw                ! Velocity conponent differences
  real(wp) :: SRp,SLm                   ! Wave speeds for the HLL part
  real(wp) :: nx1, ny1, nz1             ! Vector along which HLL is applied
  real(wp) :: nx2, ny2, nz2             ! Vector along which Roe is applied
  real(wp) :: alpha1, alpha2            ! Projections of the new normals
  real(wp) :: abs_dq                    ! Magnitude of the velocity difference
  real(wp) :: temp, tempx, tempy, tempz ! Temporary variables
  !
  !.. Local Arrays
  real(wp) :: LdU(1:4)            ! Wave strengths
  real(wp) :: eig(1:4)            ! Eigenvalues
  real(wp) :: ws(1:4)             ! Abs Wave speeds
  real(wp) :: R(1:nq,1:4)         ! Right-eigenvectors
  real(wp) :: dws(1:4)            ! Width of a parabolic fit for entropy fix
  real(wp) :: fL(1:nq), fR(1:nq)  ! Fluxes
  real(wp) :: diss(1:nq)          ! Dissipation term
  !
continue
  !
  ! Face normal vector (unit vector)
  !
  nx = unrm(1)
  ny = unrm(2)
  nz = unrm(3)
  !
  ! Primitive and other variables.
  !
  ! Left state
  rhoL = cvl(1)
  uL = cvl(2)/cvl(1)
  vL = cvl(3)/cvl(1)
  wL = cvl(4)/cvl(1)
  qnL = uL*nx + vL*ny + wL*nz
  presL = (gam-one)*( cvl(5) - half*rhoL*(uL*uL+vL*vL+wL*wL) )
  aL = sqrt(gam*presL/rhoL)
  HL = aL*aL/(gam-one) + half*(uL*uL+vL*vL+wL*wL)
  !
  ! Right state
  rhoR = cvr(1)
  uR = cvr(2)/cvr(1)
  vR = cvr(3)/cvr(1)
  wR = cvr(4)/cvr(1)
  qnR = uR*nx + vR*ny + wR*nz
  presR = (gam-one)*( cvr(5) - half*rhoR*(uR*uR+vR*vR+wR*wR) )
  aR = sqrt(gam*presR/rhoR)
  HR = aR*aR/(gam-one) + half*(uR*uR+vR*vR+wR*wR)
  !
  ! Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)
  !
  fL(1) = rhoL*qnL
  fL(2) = rhoL*qnL * uL + presL*nx
  fL(3) = rhoL*qnL * vL + presL*ny
  fL(4) = rhoL*qnL * wL + presL*nz
  fL(5) = rhoL*qnL * HL
  !
  fR(1) = rhoR*qnR
  fR(2) = rhoR*qnR * uR + presR*nx
  fR(3) = rhoR*qnR * vR + presR*ny
  fR(4) = rhoR*qnR * wR + presR*nz
  fR(5) = rhoR*qnR * HR
  !
  ! ----------------------------------------------------------------------------
  ! Define n1 and n2, and compute alpha1 and alpha2: (4.2) in original paper.
  ! NOTE: n1 and n2 may need to be frozen at some point during 
  !       a steady calculation to fully make it converge. For time-accurate 
  !       calculation, it is not necessary.
  ! NOTE: For a boundary face, you may want to set (nx2,ny2,nz2)=(nx,ny,nz),
  !       (nx1,ny1)=tangent vector, i.e., use the Roe flux.
  !
  abs_dq = sqrt( (uR-uL)**2 + (vR-vL)**2 + (wR-wL)**2 )
  !
  if (abs_dq > eps12) then
    !
    ! n1 = Velocity difference vector: normal to shock or tangent to shear
    !
    nx1 = (uR-uL)/abs_dq
    ny1 = (vR-vL)/abs_dq
    nz1 = (wR-wL)/abs_dq
    !
  else
    !
    ! n1 = Face tangent vector if abs_dq is too small.
    !
    ! NOTE: There are infinitely many choices for the tangent vector.
    !       The best choice may be discovered in future.
    !       Here, we choose a vector in a plane (essentially 2D vector).
    !
    ! NOTE: We must be careful and make sure the vector is not a zero vector.
    !
    ! E.g., if n = (1,0,0), then, (0,-nz,ny) = (0,0,0). This is a problem.
    !       To avoid such a situation, we choose the 2D tangential vector
    !       having the maximum length.
    !
    tempx = ny*ny + nz*nz
    tempy = nz*nz + nx*nx
    tempz = nx*nx + ny*ny
    !
    if ( tempx >= tempy .and. tempx >= tempz ) then
      nx1 =  zero
      ny1 = -nz
      nz1 =  ny
    else if ( tempy >= tempx .and. tempy >= tempz ) then
      nx1 = -nz
      ny1 =  zero
      nz1 =  nx
    else if ( tempz >= tempx .and. tempz >= tempy ) then
      nx1 = -ny
      ny1 =  nx
      nz1 =  zero
    end if
    !
    ! Make it the unit vector.
    temp = sqrt( nx1*nx1 + ny1*ny1 + nz1*nz1 )
    nx1 = nx1/temp
    ny1 = ny1/temp
    nz1 = nz1/temp
    !
  end if
  !
  alpha1 = nx*nx1 + ny*ny1 + nz*nz1
  !
  ! Make alpha1 always positive.
  temp = sign(one,alpha1)
  nx1 = temp * nx1
  ny1 = temp * ny1
  nz1 = temp * nz1
  alpha1 = temp * alpha1
  !
  ! n2 = direction perpendicular to n1.
  ! NOTE: There are infinitely many choices for this vector.
  !       The best choice may be discovered in future.
  ! Here, we employ the formula (4.4) in the paper:
  !     (nx2,ny2,nz2) = (n1xn)xn1 / |(n1xn)xn1|    ('x' is the vector product.)
  !
  ! (tempx,tempy,tempz) = n1xn
  tempx = ny1*nz - nz1*ny
  tempy = nz1*nx - nx1*nz
  tempz = nx1*ny - ny1*nx
  !
  ! (nx2,ny2,nz2) = (n1xn)xn1
  nx2 = tempy*nz1 - tempz*ny1
  ny2 = tempz*nx1 - tempx*nz1
  nz2 = tempx*ny1 - tempy*nx1
  !
  ! Make n2 the unit vector
  temp = sqrt( nx2*nx2 + ny2*ny2 + nz2*nz2 )
  nx2 = nx2/temp
  ny2 = ny2/temp
  nz2 = nz2/temp
  !
  alpha2 = nx*nx2 + ny*ny2 + nz*nz2
  !
  ! Make alpha2 always positive.
  temp = sign(one,alpha2)
  nx2 = temp * nx2
  ny2 = temp * ny2
  nz2 = temp * nz2
  alpha2 = temp * alpha2
  !
  ! ----------------------------------------------------------------------------
  ! Now we are going to compute the Roe flux with n2
  ! as the normal with modified wave speeds (5.12).
  !
  ! NOTE: The Roe flux here is computed without tangent vectors.
  !       See "I do like CFD, VOL.1, Second Edition" for details:
  !       page 80, Equation (3.6.28).
  !
  ! First compute the Roe-averaged quantities
  !
  ! NOTE: See http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
  !       the Roe-averaged density.
  !
  RT = sqrt(rhoR/rhoL)
  rho = RT*rhoL                                      ! Roe-ave density.
  u = (uL + RT*uR)/(one + RT)                        ! Roe-ave x-velocity
  v = (vL + RT*vR)/(one + RT)                        ! Roe-ave y-velocity
  w = (wL + RT*wR)/(one + RT)                        ! Roe-ave z-velocity
  H = (HL + RT*HR)/(one + RT)                        ! Roe-ave total enthalpy
  a = sqrt( (gam-one)*(H-half*(u*u + v*v + w*w)) ) ! Roe-ave speed of sound
  !
  ! ----------------------------------------------------
  ! Compute the wave speed estimates for the HLL part,
  ! following Einfeldt:
  !
  ! B. Einfeldt, On Godunov-type methods for gas dynamics,
  ! SIAM Journal on Numerical Analysis 25 (2) (1988) 294--318.
  !
  ! NOTE: HLL is actually applied to n1, but this is
  !       all we need to incorporate HLL. See JCP2008 paper.
  !
  qn  = u *nx1 + v *ny1 + w *nz1
  qnL = uL*nx1 + vL*ny1 + wL*nz1
  qnR = uR*nx1 + vR*ny1 + wR*nz1
  SLm = min( zero, qn - a, qnL - aL ) ! Minimum wave speed estimate
  SRp = max( zero, qn + a, qnR + aR ) ! Maximum wave speed estimate
  !
  ! This is the only place where n1=(nx1,ny1,nz1) is used.
  ! n1=(nx1,ny1,nz1) is never used below.
  ! ----------------------------------------------------
  !
  ! Wave Strengths
  !
  qn  = u *nx2 + v *ny2 + w *nz2
  qnL = uL*nx2 + vL*ny2 + wL*nz2
  qnR = uR*nx2 + vR*ny2 + wR*nz2
  !
  drho = rhoR - rhoL  ! Density difference
  dp =   presR - presL    ! Pressure difference
  dqn =  qnR - qnL   ! Normal velocity difference
  !
  LdU(1) = (dp - rho*a*dqn )/(two*a*a) ! Left-moving acoustic wave strength
  LdU(2) =  drho - dp/(a*a)            ! Entropy wave strength
  LdU(3) = (dp + rho*a*dqn )/(two*a*a) ! Right-moving acoustic wave strength
  LdU(4) = rho                 ! Shear wave strength (not really, just a factor)
  !
  ! Wave Speed (Eigenvalues)
  !
  eig(1) = qn-a ! Left-moving acoustic wave velocity
  eig(2) = qn   ! Entropy wave velocity
  eig(3) = qn+a ! Right-moving acoustic wave velocity
  eig(4) = qn   ! Shear wave velocity
  !
  ! Absolute values of the wave speeds (Eigenvalues)
  !
  ws(1) = abs(qn-a) ! Left-moving acoustic wave speed
  ws(2) = abs(qn)   ! Entropy wave speed
  ws(3) = abs(qn+a) ! Right-moving acoustic wave speed
  ws(4) = abs(qn)   ! Shear wave speed
  !
  ! Harten's Entropy Fix JCP(1983), 49, pp357-393: only for the
  ! nonlinear fields.
  ! NOTE: It avoids vanishing wave speeds by making a parabolic fit near ws = 0.
  !
  dws(1) = one5
  if ( ws(1) < dws(1) ) then
    ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  end if
  !
  dws(3) = one5
  if ( ws(3) < dws(3) ) then
    ws(3) = half * ( ws(3)*ws(3)/dws(3)+dws(3) )
  end if
  !
  ! Combine wave speeds for Rotated-RHLL: Eq.(5.12) in original JCP2008 paper.
  !
  ws = alpha2*ws - (alpha1*two*SRp*SLm + alpha2*(SRp+SLm)*eig)/(SRp-SLm)
  !
  ! Below, we compute the Roe dissipation term in the direction n2 with the
  ! above modified wave speeds. HLL wave speeds act something like the entropy
  ! fix or eigenvalue limiting; they contribute only by the amount given by
  ! the fraction, alpha1 (less than or equal to 1.0). See JCP2008 paper.
  !
  ! Right Eigenvectors:
  ! NOTE: Two shear wave components are combined into one, so that tangent
  !       vectors are not required. And thats why there are only 4 vectors here.
  !
  ! Left-moving acoustic wave
  R(1,1) = one    
  R(2,1) = u - a*nx2
  R(3,1) = v - a*ny2
  R(4,1) = w - a*nz2
  R(5,1) = H - a*qn
  !
  ! Entropy wave
  R(1,2) = one
  R(2,2) = u
  R(3,2) = v 
  R(4,2) = w
  R(5,2) = half*(u*u + v*v + w*w)
  !
  ! Right-moving acoustic wave
  R(1,3) = one
  R(2,3) = u + a*nx2
  R(3,3) = v + a*ny2
  R(4,3) = w + a*nz2
  R(5,3) = H + a*qn
  !
  ! Two shear wave components combined into one (wave strength incorporated).
  du = uR - uL
  dv = vR - vL
  dw = wR - wL
  !
  R(1,4) = zero
  R(2,4) = du - dqn*nx2
  R(3,4) = dv - dqn*ny2
  R(4,4) = dw - dqn*nz2
  R(5,4) = u*du + v*dv + w*dw - qn*dqn
  !
  ! Dissipation Term: Roe dissipation with the modified wave speeds.
  !    |An|dU = R|Lambda|L*dU
  !           = sum_k of [ ws(k) * R(:,k) * L*dU(k) ], where n=n2.
  !
  diss(:) = ws(1)*LdU(1)*R(:,1) + &
            ws(2)*LdU(2)*R(:,2) + &
            ws(3)*LdU(3)*R(:,3) + &
            ws(4)*LdU(4)*R(:,4)
  !
  ! Compute the Rotated-RHLL flux, it looks like HLL flux with Roe dissipation
  !
  return_value = (SRp*fL - SLm*fR)/(SRp-SLm) - half*diss
  !
  ! Normal max wave speed in the normal direction.
  !  wsn = abs(qn) + a
  !
end function rotated_rhll_upwind_flux_cv_3d
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
!###############################################################################
!###############################################################################
!
pure subroutine roe_upwind_flux_sub(unrm,cvl,cvr,return_value)
  !$ACC ROUTINE SEQ
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),    intent(in) :: unrm
  real(wp), dimension(1:nq),    intent(in) :: cvl
  real(wp), dimension(1:nq),    intent(in) :: cvr
  real(wp), dimension(1:nq), intent(inout) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: gm1,gp1
  real(wp) :: rl,presl,hl,q2l,ql,rql,betal
  real(wp) :: rr,presr,hr,q2r,qr,rqr,betar
  real(wp) :: dw1,hsr1,eta1,udw1
  real(wp) :: dw4,hsr4,eta4,udw4
  real(wp) :: rtd,htd,atd2,atd,qtd
  real(wp) :: atd2nz,atd2rt,atdnz,atdrt
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: vell,velr,veltd
  !
continue
  !
  gm1 = gam - one
  gp1 = gam + one
  !
  rl = cvl(nec)
  vell = cvl(nmb:nme)/cvl(nec)
  q2l = selfdot(vell)
  presl = gm1*(cvl(nee) - half*rl*q2l)
  !
  rr = cvr(nec)
  velr = cvr(nmb:nme)/cvr(nec)
  q2r = selfdot(velr)
  presr = gm1*(cvr(nee) - half*rr*q2r)
  !
  hl = half*q2l + gam/gm1 * presl/rl
  hr = half*q2r + gam/gm1 * presr/rr
  ql = dot_product(vell,unrm)
  qr = dot_product(velr,unrm)
  !
  !... square root averaging
  rtd   = sqrt(rr*rl)
  betal = rl / (rl+rtd)
  betar = one - betal
  veltd = betal*vell + betar*velr
  htd   = betal*hl + betar*hr
  atd2  = max(gm1*(htd - half*selfdot(veltd)),zero)
  atd   = sqrt(atd2)
  qtd   = dot_product(veltd,unrm)
  !
  atd2nz = max(atd2,rndoff)
  atdnz = max(atd,rndoff)
  atd2rt = atd2 / atd2nz
  atdrt = atd / atdnz
  !
  if (qtd >= zero) then
    !
    dw1  = half*((presr-presl)/atd2nz*atd2rt - (qr-ql)*rtd/atdnz*atdrt)
    hsr1 = -one4*gp1*dw1*atd/rtd
    eta1 = max(-abs(qtd-atd) + hsr1, zero)
    udw1 = dw1*(min(qtd-atd,zero) - half*eta1)
    rql  = rl*ql
    !
    return_value(nec)     =      rql + udw1
    return_value(nmb:nme) = vell*rql + udw1*(veltd-atd*unrm) + presl*unrm
    return_value(nee)     =   hl*rql + udw1*(  htd-atd*qtd)
    !
  else
    !
    dw4  = half*((presr-presl)/atd2nz*atd2rt + (qr-ql)*rtd/atdnz*atdrt)
    hsr4 = one4*gp1*dw4*atd/rtd
    eta4 = max(-abs(qtd+atd) + hsr4, zero)
    udw4 = dw4*(max(qtd+atd,zero) + half*eta4)
    rqr  = rr*qr
    !
    return_value(nec)     =      rqr - udw4
    return_value(nmb:nme) = velr*rqr - udw4*(veltd+atd*unrm) + presr*unrm
    return_value(nee)     =   hr*rqr - udw4*(  htd+atd*qtd)
    !
  end if
  !
end subroutine roe_upwind_flux_sub
!
!###############################################################################
!
pure subroutine hllc_upwind_flux_sub(unrm,cvl,cvr,return_value)
  !$ACC ROUTINE SEQ
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),    intent(in) :: unrm
  real(wp), dimension(1:nq),    intent(in) :: cvl
  real(wp), dimension(1:nq),    intent(in) :: cvr
  real(wp), dimension(1:nq), intent(inout) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: gm1
  real(wp) :: rl,presl,hl,al,q2l,vnl,rhoel,wgl,sl,dsl
  real(wp) :: rr,presr,hr,ar,q2r,vnr,rhoer,wgr,sr,dsr
  real(wp) :: srt,ra,ha,q2a,vna,aa
  real(wp) :: rho,rhou,rhov,rhoe,pstar,sm
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: vell,velr,vela,rvel
  !
continue
  !
  gm1 = gam - one
  !
  rl = cvl(nec)
  vell = cvl(nmb:nme)/cvl(nec)
  q2l = selfdot(vell)
  presl = gm1*(cvl(nee) - half*rl*q2l)
  hl = half*q2l + gam/gm1 * presl/rl
  vnl = dot_product(vell,unrm)
  rhoel = cvl(nee)
  !
  rr = cvr(nec)
  velr = cvr(nmb:nme)/cvr(nec)
  q2r = selfdot(velr)
  presr = gm1*(cvr(nee) - half*rr*q2r)
  hr = half*q2r + gam/gm1 * presr/rr
  vnr = dot_product(velr,unrm)
  rhoer = cvr(nee)
  !
  srt = sqrt(max(rr/rl,rndoff))
  wgl = one/(srt+one)
  wgr = srt/(srt+one)
  !
  ra = sqrt(max(rl*rr,rndoff))
  vela = vell*wgl + velr*wgr
  ha = hl*wgl + hr*wgr
  !
  q2a = selfdot(vela)
  vna = dot_product(vela,unrm)
  !
  aa = sqrt(max(gm1*(ha-half*q2a),rndoff))
  !
  sl = min( vnl-al , vna-aa )
  sr = max( vnr+ar , vna+aa )
  !
  if (sl > zero) then
    !
    return_value(nec)     = vnl*rl
    return_value(nmb:nme) = vnl*rl*vell + presl*unrm
    return_value(nee)     = vnl*rhoel   + presl*vnl
    !
  else if (sr < zero) then
    !
    return_value(nec)     = vnr*rr
    return_value(nmb:nme) = vnr*rr*velr + presr*unrm
    return_value(nee)     = vnr*rhoer   + presr*vnr
    !
  else
    !
    dsl = sl - vnl
    dsr = sr - vnr
    !
    sm = ((rr*vnr*dsr-presr) - (rl*vnl*dsl-presl)) / (rr*dsr-rl*dsl)
    !
    pstar = presl*dsl*(sm-vnl) + presl
    !
    if (sm > zero) then
      !
      rho  = dsl*rl/(sl-sm)
      rvel = (dsl*rl*vell + (pstar-presl)*unrm)/(sl-sm)
      rhoe = (dsl*rhoel - presl*vnl + pstar*sm)/(sl-sm)
      !
    else
      !
      rho  = dsr*rr/(sr-sm)
      rvel = (dsr*rr*velr + (pstar-presr)*unrm)/(sr-sm)
      rhoe = (dsr*rhoer - presr*vnr + pstar*sm)/(sr-sm)
      !
    end if
    !
    return_value(nec)     = sm*rho
    return_value(nmb:nme) = sm*rvel + pstar*unrm
    return_value(nee)     = sm*rhoe + pstar*sm
    !
  end if
  !
end subroutine hllc_upwind_flux_sub
!
!###############################################################################
!
pure subroutine ldfss_upwind_flux_sub(unrm,cvl,cvr,return_value)
  !$ACC ROUTINE SEQ
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),    intent(in) :: unrm
  real(wp), dimension(1:nq),    intent(in) :: cvl
  real(wp), dimension(1:nq),    intent(in) :: cvr
  real(wp), dimension(1:nq), intent(inout) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: gm1,ahalf
  real(wp) :: rl,presl,hl,al,q2l
  real(wp) :: rr,presr,hr,ar,q2r
  real(wp) :: alpha_l,beta_l,mach_l,mstar_l
  real(wp) :: alpha_r,beta_r,mach_r,mstar_r
  real(wp) :: machp_half,uc_l,flux_l,pplus_l
  real(wp) :: machm_half,uc_r,flux_r,pplus_r
  real(wp) :: m_half,mstar_half
  real(wp) :: del_p,sum_p,p_net
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: vell,velr
  !
continue
  !
  gm1 = gam - one
  !
  rl = cvl(nec)
  vell = cvl(nmb:nme)/cvl(nec)
  q2l = selfdot(vell)
  presl = gm1*(cvl(nee) - half*rl*q2l)
  hl = half*q2l + gam/gm1 * presl/rl
  al = sqrt(max(gam*presl/rl,rndoff))
  !
  rr = cvr(nec)
  velr = cvr(nmb:nme)/cvr(nec)
  q2r = selfdot(velr)
  presr = gm1*(cvr(nee) - half*rr*q2r)
  hr = half*q2r + gam/gm1 * presr/rr
  ar = sqrt(max(gam*presr/rr,rndoff))
  !
  ahalf = half*(al+ar)
  !
  mach_l = dot_product(vell,unrm)/ahalf
  mach_r = dot_product(velr,unrm)/ahalf
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
  del_p = presl - presr
  sum_p = presl + presr
  !
  machp_half = mstar_half*(one - del_p/sum_p - two*abs(del_p)/presl)
  machm_half = mstar_half*(one + del_p/sum_p - two*abs(del_p)/presr)
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
  p_net = presl*(alpha_l*(one+beta_l) - beta_l*pplus_l) + &
          presr*(alpha_r*(one+beta_r) - beta_r*pplus_r)
  !
  return_value(nec)     = flux_l      + flux_r
  return_value(nmb:nme) = flux_l*vell + flux_r*velr + p_net*unrm
  return_value(nee)     = flux_l*hl   + flux_r*hr
  !
end subroutine ldfss_upwind_flux_sub
!
!###############################################################################
!
pure subroutine ausm_plus_upwind_flux_sub(unrm,pvl,pvr,return_value)
  !$ACC ROUTINE SEQ
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),    intent(in) :: unrm
  real(wp), dimension(1:nq),    intent(in) :: pvl
  real(wp), dimension(1:nq),    intent(in) :: pvr
  real(wp), dimension(1:nq), intent(inout) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: gm1
  real(wp) :: rl,presl,hl,al,v2l,a2l,h0l,uul,f1l
  real(wp) :: rr,presr,hr,ar,v2r,a2r,h0r,uur,f1r
  real(wp) :: alr,alr2,rm2,rmref2,coef1,coef2,prec_f
  real(wp) :: aml,absml,am1p,aml1,addp,aml2,ampl,ppl
  real(wp) :: amr,absmr,am1m,amr1,addm,amr2,ammr,pmr
  real(wp) :: amhalf,denom,amhalfp,amhalfm,phalf,amco2
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: vell,velr
  !
  !.. Local Parameters ..
  real(wp), parameter :: amco2m = 1.0e-4_wp         
  real(wp), parameter :: alphax = 0.1875_wp
  real(wp), parameter :: betax  = 0.125_wp
  !
continue
  !
  gm1 = gam - one
  !
  rl = pvl(nec)
  vell = pvl(nmb:nme)
  presl = pvl(nee)
  !
  rr = pvr(nec)
  velr = pvr(nmb:nme)
  presr = pvr(nee)
  !
  hl  = gam*presl/(gm1*rl)
  v2l = half * selfdot(vell)
  h0l = hl + v2l
  a2l = gm1*hl
  al  = sqrt(a2l)
  uul = dot_product(vell,unrm)
  !
  hr  = gam*presr/(gm1*rr)
  v2r = half * selfdot(velr)
  h0r = hr + v2r
  a2r = gm1*hr
  ar  = sqrt(a2r)
  uur = dot_product(velr,unrm)
  !
  alr    = half*(al + ar)
  alr2   = alr**2
  rm2    = (v2r+v2l)/alr2
  amco2  = max(amco2m,one4*rm2)
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
    aml2 = one4*aml1**2
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
    amr2 = one4*amr1**2
    ammr = -amr2 - betax*addm
    pmr  = amr2*(two + amr) - alphax*amr*addm
  else
    ammr = am1m
    pmr  = am1m/amr
  end if
  !
  amhalf  = ampl + ammr
  denom   = ( rmref2*(presl+presr)*(rl+rr) - &
                zero*(presl-presr)*(rl-rr) ) / sqrt(rl*rr)
  amhalf  = amhalf + half*(presl-presr)*((ampl-ammr) - (am1p-am1m))/denom
  amhalfp = half*(amhalf + abs(amhalf))
  amhalfm = amhalf - amhalfp
  !
  phalf = ppl*presl + pmr*presr - ppl*pmr*(rl+rr)/two*alr*(uur-uul)
  !
  f1l = rl*amhalfp*alr
  f1r = rr*amhalfm*alr
  !
  return_value(nec)     = f1l      + f1r
  return_value(nmb:nme) = f1l*vell + f1r*velr + unrm*phalf
  return_value(nee)     = f1l*h0l  + f1r*h0r
  !
end subroutine ausm_plus_upwind_flux_sub
!
!###############################################################################
!
pure subroutine rotated_rhll_upwind_flux_sub_2d(unrm,cvl,cvr,return_value)
  !$ACC ROUTINE SEQ
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),    intent(in) :: unrm ! unit face normal
  real(wp), dimension(1:nq),    intent(in) :: cvl
  real(wp), dimension(1:nq),    intent(in) :: cvr
  real(wp), dimension(1:nq), intent(inout) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: nx, ny
  real(wp) :: nx1, ny1, nx2, ny2      ! Rotated normals, n1 and n2
  real(wp) :: tx, ty                  ! Tangent vector (taken as n1)
  real(wp) :: alpha1, alpha2          ! Projections of the new normals
  real(wp) :: vxL, vxR, vyL, vyR      ! Velocity components.
  real(wp) :: rhoL, rhoR, presL, presR      ! Primitive variables.
  real(wp) :: vnL, vnR, vtL, vtR      ! Normal and tagent velocities
  real(wp) :: aL, aR, HL, HR          ! Speeds of sound and total enthalpy
  real(wp) :: RT,rho,vx,vy,H,a        ! Roe-averages
  real(wp) :: vn, vt                  ! Normal and tagent velocities (Roe-ave)
  real(wp) :: drho,dvx,dvy,dvn,dvt,dp ! Wave strenghs
  real(wp) :: abs_dq                  ! Magnitude of the velocity difference
  real(wp) :: SRp,SLm                 ! Wave speeds for the HLL part
  real(wp) :: temp
  !
  !.. Local Arrays ..
  real(wp) :: dV(1:4)              ! Wave strenghs
  real(wp) :: abs_ws(1:4), ws(1:4) ! Abs wave speeds
  real(wp) :: dws(1:4)             ! width of a parabolic fit for entropy fix
  real(wp) :: Rv(1:nq,1:4)         ! Right-eigevectors
  real(wp) :: fL(1:nq), fR(1:nq)   ! Fluxes
  real(wp) :: diss(1:nq)           ! Dissipation term
  !
continue
  !
  ! Face normal vector (unit vector)
  !
  nx = unrm(1)
  ny = unrm(2)
  !
  ! Primitive and other variables.
  !
  ! Left state
  !
  rhoL = cvl(1)
  vxL = cvl(2)/cvl(1)
  vyL = cvl(3)/cvl(1)
  presL = (gam-one)*( cvl(4) - half*rhoL*(vxL*vxL+vyL*vyL) )
  aL = sqrt(gam*presL/rhoL)
  HL = ( cvl(4) + presL ) / rhoL
  !
  ! Right state
  !
  rhoR = cvr(1)
  vxR = cvr(2)/cvr(1)
  vyR = cvr(3)/cvr(1)
  presR = (gam-one)*( cvr(4) - half*rhoR*(vxR*vxR+vyR*vyR) )
  aR = sqrt(gam*presR/rhoR)
  HR = ( cvr(4) + presR ) / rhoR
  !
  vnL = vxL*nx + vyL*ny
  vnR = vxR*nx + vyR*ny
  !
  ! Compute the flux.
  !
  fL(1) = rhoL*vnL
  fL(2) = rhoL*vnL * vxL + presL*nx
  fL(3) = rhoL*vnL * vyL + presL*ny
  fL(4) = rhoL*vnL *  HL
  !
  fR(1) = rhoR*vnR
  fR(2) = rhoR*vnR * vxR + presR*nx
  fR(3) = rhoR*vnR * vyR + presR*ny
  fR(4) = rhoR*vnR *  HR
  !
  ! Define n1 and n2, and compute alpha1 and alpha2: (4.2) in original paper.
  ! (NB: n1 and n2 may need to be frozen at some point during 
  !      a steady calculation to fully make it converge. For time-accurate 
  !      calculation, this is fine.)
  ! NB: For a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).
  !
  abs_dq = sqrt( (vxR-vxL)**2+(vyR-vyL)**2 )
  !
  if (abs_dq > eps12) then
    nx1 = (vxR-vxL)/abs_dq
    ny1 = (vyR-vyL)/abs_dq
  else
    nx1 = -ny 
    ny1 =  nx
  endif
  !
  alpha1 = nx * nx1 + ny * ny1 
  !
  ! To make alpha1 always positive.
  !
  temp = sign(one,alpha1)
  nx1 = temp * nx1
  ny1 = temp * ny1
  alpha1 = temp * alpha1
  !
  ! Take n2 as perpendicular to n1.
  !
  nx2 = -ny1
  ny2 =  nx1
  alpha2 = nx * nx2 + ny * ny2
  !
  ! To make alpha2 always positive.
  !
  temp = sign(one,alpha2)
  nx2 = temp * nx2
  ny2 = temp * ny2
  alpha2 = temp * alpha2
  !
  ! Now we are going to compute the Roe flux with n2 as the normal
  ! and n1 as the tagent vector, with modified wave speeds (5.12)
  !
  ! Compute the Roe Averages
  !
  RT = sqrt(rhoR/rhoL)
  rho = RT*rhoL
  vx = (vxL+RT*vxR)/(one+RT)
  vy = (vyL+RT*vyR)/(one+RT)
  H = ( HL+RT* HR)/(one+RT)
  a = sqrt( (gam-one)*(H-half*(vx*vx+vy*vy)) )
  vn = vx*nx2+vy*ny2
  vt = vx*nx1+vy*ny1
  !
  ! Wave Strengths (remember that n2 is the normal and n1 is the tangent.)
  !
  vnL = vxL*nx2 + vyL*ny2
  vnR = vxR*nx2 + vyR*ny2
  vtL = vxL*nx1 + vyL*ny1
  vtR = vxR*nx1 + vyR*ny1
  !
  drho = rhoR - rhoL 
  dp =   presR - presL
  dvn =  vnR - vnL
  dvt =  vtR - vtL
  !
  dV(1) = (dp - rho*a*dvn )/(two*a*a)
  dV(2) =  rho*dvt/a
  dV(3) =  drho - dp/(a*a)
  dV(4) = (dp + rho*a*dvn )/(two*a*a)
  !
  ! Wave Speeds for Roe flux part.
  !
  ws(1) = vn-a
  ws(2) = vn
  ws(3) = vn
  ws(4) = vn+a
  abs_ws  = abs(ws)
  !
  ! Harten's Entropy Fix JCP(1983), 49, pp357-393:
  ! only for the nonlinear fields.
  !
  dws(1) = one5
  if ( abs_ws(1) < dws(1) ) then
    abs_ws(1) = half*(abs_ws(1)*abs_ws(1)/dws(1)+dws(1))
  end if
  !
  dws(4) = one5
  if ( abs_ws(4) < dws(4) ) then
    abs_ws(4) = half*(abs_ws(4)*abs_ws(4)/dws(4)+dws(4))
  end if
  !
  ! HLL wave speeds, evaluated with [nx1,ny1] (=tangent wrt n2).
  !
  SRp = max( zero, vtR + aR, vt + a)
  SLm = min( zero, vtL - aL, vt - a)
  !
  ! Modified wave speeds for the Rotated-RHLL flux: (5.12) in original paper.
  !
  ws = alpha2*abs_ws - ( alpha2*(SRp+SLm)*ws + two*alpha1*SRp*SLm )/ (SRp-SLm)
  !
  ! Right Eigenvectors: with n2 as normal and n1 as tangent.
  !
  tx = nx1
  ty = ny1
  !
  Rv(1,1) = one    
  Rv(2,1) = vx - a*nx2
  Rv(3,1) = vy - a*ny2
  Rv(4,1) =  H - vn*a
  !
  Rv(1,2) = zero
  Rv(2,2) = a*tx
  Rv(3,2) = a*ty
  Rv(4,2) = a*vt
  !
  Rv(1,3) = one
  Rv(2,3) = vx
  Rv(3,3) = vy 
  Rv(4,3) = half*(vx*vx+vy*vy)
  !
  Rv(1,4) = one
  Rv(2,4) = vx + a*nx2
  Rv(3,4) = vy + a*ny2
  Rv(4,4) =  H + vn*a
  !
  ! Dissipation Term: Roe dissipation with the modified wave speeds.
  !
  diss(:) = ws(1)*dV(1)*Rv(:,1) + &
            ws(2)*dV(2)*Rv(:,2) + &
            ws(3)*dV(3)*Rv(:,3) + &
            ws(4)*dV(4)*Rv(:,4)
  !
  ! Compute the Rotated-RHLL flux.
  !
  return_value = (SRp*fL - SLm*fR)/(SRp-SLm) - half*diss
  !
end subroutine rotated_rhll_upwind_flux_sub_2d
!
!###############################################################################
!
pure subroutine rotated_rhll_upwind_flux_sub_3d(unrm,cvl,cvr,return_value)
  !$ACC ROUTINE SEQ
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),    intent(in) :: unrm ! unit face normal
  real(wp), dimension(1:nq),    intent(in) :: cvl
  real(wp), dimension(1:nq),    intent(in) :: cvr
  real(wp), dimension(1:nq), intent(inout) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: nx, ny, nz                ! Face normal vector
  real(wp) :: uL, uR, vL, vR, wL, wR    ! Velocity components.
  real(wp) :: rhoL, rhoR, presL, presR        ! Primitive variables.
  real(wp) :: qnL, qnR                  ! Normal velocity
  real(wp) :: aL, aR, HL, HR            ! Speed of sound, Total enthalpy
  real(wp) :: RT,rho,u,v,w,H,a,qn       ! Roe-averages
  real(wp) :: drho,dqn,dp               ! Wave strengths
  real(wp) :: du, dv, dw                ! Velocity conponent differences
  real(wp) :: SRp,SLm                   ! Wave speeds for the HLL part
  real(wp) :: nx1, ny1, nz1             ! Vector along which HLL is applied
  real(wp) :: nx2, ny2, nz2             ! Vector along which Roe is applied
  real(wp) :: alpha1, alpha2            ! Projections of the new normals
  real(wp) :: abs_dq                    ! Magnitude of the velocity difference
  real(wp) :: temp, tempx, tempy, tempz ! Temporary variables
  !
  !.. Local Arrays
  real(wp) :: LdU(1:4)            ! Wave strengths
  real(wp) :: eig(1:4)            ! Eigenvalues
  real(wp) :: ws(1:4)             ! Abs Wave speeds
  real(wp) :: R(1:nq,1:4)         ! Right-eigenvectors
  real(wp) :: dws(1:4)            ! Width of a parabolic fit for entropy fix
  real(wp) :: fL(1:nq), fR(1:nq)  ! Fluxes
  real(wp) :: diss(1:nq)          ! Dissipation term
  !
continue
  !
  ! Face normal vector (unit vector)
  !
  nx = unrm(1)
  ny = unrm(2)
  nz = unrm(3)
  !
  ! Primitive and other variables.
  !
  ! Left state
  rhoL = cvl(1)
  uL = cvl(2)/cvl(1)
  vL = cvl(3)/cvl(1)
  wL = cvl(4)/cvl(1)
  qnL = uL*nx + vL*ny + wL*nz
  presL = (gam-one)*( cvl(5) - half*rhoL*(uL*uL+vL*vL+wL*wL) )
  aL = sqrt(gam*presL/rhoL)
  HL = aL*aL/(gam-one) + half*(uL*uL+vL*vL+wL*wL)
  !
  ! Right state
  rhoR = cvr(1)
  uR = cvr(2)/cvr(1)
  vR = cvr(3)/cvr(1)
  wR = cvr(4)/cvr(1)
  qnR = uR*nx + vR*ny + wR*nz
  presR = (gam-one)*( cvr(5) - half*rhoR*(uR*uR+vR*vR+wR*wR) )
  aR = sqrt(gam*presR/rhoR)
  HR = aR*aR/(gam-one) + half*(uR*uR+vR*vR+wR*wR)
  !
  ! Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)
  !
  fL(1) = rhoL*qnL
  fL(2) = rhoL*qnL * uL + presL*nx
  fL(3) = rhoL*qnL * vL + presL*ny
  fL(4) = rhoL*qnL * wL + presL*nz
  fL(5) = rhoL*qnL * HL
  !
  fR(1) = rhoR*qnR
  fR(2) = rhoR*qnR * uR + presR*nx
  fR(3) = rhoR*qnR * vR + presR*ny
  fR(4) = rhoR*qnR * wR + presR*nz
  fR(5) = rhoR*qnR * HR
  !
  ! ----------------------------------------------------------------------------
  ! Define n1 and n2, and compute alpha1 and alpha2: (4.2) in original paper.
  ! NOTE: n1 and n2 may need to be frozen at some point during 
  !       a steady calculation to fully make it converge. For time-accurate 
  !       calculation, it is not necessary.
  ! NOTE: For a boundary face, you may want to set (nx2,ny2,nz2)=(nx,ny,nz),
  !       (nx1,ny1)=tangent vector, i.e., use the Roe flux.
  !
  abs_dq = sqrt( (uR-uL)**2 + (vR-vL)**2 + (wR-wL)**2 )
  !
  if (abs_dq > eps12) then
    !
    ! n1 = Velocity difference vector: normal to shock or tangent to shear
    !
    nx1 = (uR-uL)/abs_dq
    ny1 = (vR-vL)/abs_dq
    nz1 = (wR-wL)/abs_dq
    !
  else
    !
    ! n1 = Face tangent vector if abs_dq is too small.
    !
    ! NOTE: There are infinitely many choices for the tangent vector.
    !       The best choice may be discovered in future.
    !       Here, we choose a vector in a plane (essentially 2D vector).
    !
    ! NOTE: We must be careful and make sure the vector is not a zero vector.
    !
    ! E.g., if n = (1,0,0), then, (0,-nz,ny) = (0,0,0). This is a problem.
    !       To avoid such a situation, we choose the 2D tangential vector
    !       having the maximum length.
    !
    tempx = ny*ny + nz*nz
    tempy = nz*nz + nx*nx
    tempz = nx*nx + ny*ny
    !
    if ( tempx >= tempy .and. tempx >= tempz ) then
      nx1 =  zero
      ny1 = -nz
      nz1 =  ny
    else if ( tempy >= tempx .and. tempy >= tempz ) then
      nx1 = -nz
      ny1 =  zero
      nz1 =  nx
    else if ( tempz >= tempx .and. tempz >= tempy ) then
      nx1 = -ny
      ny1 =  nx
      nz1 =  zero
    end if
    !
    ! Make it the unit vector.
    temp = sqrt( nx1*nx1 + ny1*ny1 + nz1*nz1 )
    nx1 = nx1/temp
    ny1 = ny1/temp
    nz1 = nz1/temp
    !
  end if
  !
  alpha1 = nx*nx1 + ny*ny1 + nz*nz1
  !
  ! Make alpha1 always positive.
  temp = sign(one,alpha1)
  nx1 = temp * nx1
  ny1 = temp * ny1
  nz1 = temp * nz1
  alpha1 = temp * alpha1
  !
  ! n2 = direction perpendicular to n1.
  ! NOTE: There are infinitely many choices for this vector.
  !       The best choice may be discovered in future.
  ! Here, we employ the formula (4.4) in the paper:
  !     (nx2,ny2,nz2) = (n1xn)xn1 / |(n1xn)xn1|    ('x' is the vector product.)
  !
  ! (tempx,tempy,tempz) = n1xn
  tempx = ny1*nz - nz1*ny
  tempy = nz1*nx - nx1*nz
  tempz = nx1*ny - ny1*nx
  !
  ! (nx2,ny2,nz2) = (n1xn)xn1
  nx2 = tempy*nz1 - tempz*ny1
  ny2 = tempz*nx1 - tempx*nz1
  nz2 = tempx*ny1 - tempy*nx1
  !
  ! Make n2 the unit vector
  temp = sqrt( nx2*nx2 + ny2*ny2 + nz2*nz2 )
  nx2 = nx2/temp
  ny2 = ny2/temp
  nz2 = nz2/temp
  !
  alpha2 = nx*nx2 + ny*ny2 + nz*nz2
  !
  ! Make alpha2 always positive.
  temp = sign(one,alpha2)
  nx2 = temp * nx2
  ny2 = temp * ny2
  nz2 = temp * nz2
  alpha2 = temp * alpha2
  !
  ! ----------------------------------------------------------------------------
  ! Now we are going to compute the Roe flux with n2
  ! as the normal with modified wave speeds (5.12).
  !
  ! NOTE: The Roe flux here is computed without tangent vectors.
  !       See "I do like CFD, VOL.1, Second Edition" for details:
  !       page 80, Equation (3.6.28).
  !
  ! First compute the Roe-averaged quantities
  !
  ! NOTE: See http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
  !       the Roe-averaged density.
  !
  RT = sqrt(rhoR/rhoL)
  rho = RT*rhoL                                      ! Roe-ave density.
  u = (uL + RT*uR)/(one + RT)                        ! Roe-ave x-velocity
  v = (vL + RT*vR)/(one + RT)                        ! Roe-ave y-velocity
  w = (wL + RT*wR)/(one + RT)                        ! Roe-ave z-velocity
  H = (HL + RT*HR)/(one + RT)                        ! Roe-ave total enthalpy
  a = sqrt( (gam-one)*(H-half*(u*u + v*v + w*w)) ) ! Roe-ave speed of sound
  !
  ! ----------------------------------------------------
  ! Compute the wave speed estimates for the HLL part,
  ! following Einfeldt:
  !
  ! B. Einfeldt, On Godunov-type methods for gas dynamics,
  ! SIAM Journal on Numerical Analysis 25 (2) (1988) 294--318.
  !
  ! NOTE: HLL is actually applied to n1, but this is
  !       all we need to incorporate HLL. See JCP2008 paper.
  !
  qn  = u *nx1 + v *ny1 + w *nz1
  qnL = uL*nx1 + vL*ny1 + wL*nz1
  qnR = uR*nx1 + vR*ny1 + wR*nz1
  SLm = min( zero, qn - a, qnL - aL ) ! Minimum wave speed estimate
  SRp = max( zero, qn + a, qnR + aR ) ! Maximum wave speed estimate
  !
  ! This is the only place where n1=(nx1,ny1,nz1) is used.
  ! n1=(nx1,ny1,nz1) is never used below.
  ! ----------------------------------------------------
  !
  ! Wave Strengths
  !
  qn  = u *nx2 + v *ny2 + w *nz2
  qnL = uL*nx2 + vL*ny2 + wL*nz2
  qnR = uR*nx2 + vR*ny2 + wR*nz2
  !
  drho = rhoR - rhoL  ! Density difference
  dp =   presR - presL    ! Pressure difference
  dqn =  qnR - qnL   ! Normal velocity difference
  !
  LdU(1) = (dp - rho*a*dqn )/(two*a*a) ! Left-moving acoustic wave strength
  LdU(2) =  drho - dp/(a*a)            ! Entropy wave strength
  LdU(3) = (dp + rho*a*dqn )/(two*a*a) ! Right-moving acoustic wave strength
  LdU(4) = rho                 ! Shear wave strength (not really, just a factor)
  !
  ! Wave Speed (Eigenvalues)
  !
  eig(1) = qn-a ! Left-moving acoustic wave velocity
  eig(2) = qn   ! Entropy wave velocity
  eig(3) = qn+a ! Right-moving acoustic wave velocity
  eig(4) = qn   ! Shear wave velocity
  !
  ! Absolute values of the wave speeds (Eigenvalues)
  !
  ws(1) = abs(qn-a) ! Left-moving acoustic wave speed
  ws(2) = abs(qn)   ! Entropy wave speed
  ws(3) = abs(qn+a) ! Right-moving acoustic wave speed
  ws(4) = abs(qn)   ! Shear wave speed
  !
  ! Harten's Entropy Fix JCP(1983), 49, pp357-393: only for the
  ! nonlinear fields.
  ! NOTE: It avoids vanishing wave speeds by making a parabolic fit near ws = 0.
  !
  dws(1) = one5
  if ( ws(1) < dws(1) ) then
    ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  end if
  !
  dws(3) = one5
  if ( ws(3) < dws(3) ) then
    ws(3) = half * ( ws(3)*ws(3)/dws(3)+dws(3) )
  end if
  !
  ! Combine wave speeds for Rotated-RHLL: Eq.(5.12) in original JCP2008 paper.
  !
  ws = alpha2*ws - (alpha1*two*SRp*SLm + alpha2*(SRp+SLm)*eig)/(SRp-SLm)
  !
  ! Below, we compute the Roe dissipation term in the direction n2 with the
  ! above modified wave speeds. HLL wave speeds act something like the entropy
  ! fix or eigenvalue limiting; they contribute only by the amount given by
  ! the fraction, alpha1 (less than or equal to 1.0). See JCP2008 paper.
  !
  ! Right Eigenvectors:
  ! NOTE: Two shear wave components are combined into one, so that tangent
  !       vectors are not required. And thats why there are only 4 vectors here.
  !
  ! Left-moving acoustic wave
  R(1,1) = one    
  R(2,1) = u - a*nx2
  R(3,1) = v - a*ny2
  R(4,1) = w - a*nz2
  R(5,1) = H - a*qn
  !
  ! Entropy wave
  R(1,2) = one
  R(2,2) = u
  R(3,2) = v 
  R(4,2) = w
  R(5,2) = half*(u*u + v*v + w*w)
  !
  ! Right-moving acoustic wave
  R(1,3) = one
  R(2,3) = u + a*nx2
  R(3,3) = v + a*ny2
  R(4,3) = w + a*nz2
  R(5,3) = H + a*qn
  !
  ! Two shear wave components combined into one (wave strength incorporated).
  du = uR - uL
  dv = vR - vL
  dw = wR - wL
  !
  R(1,4) = zero
  R(2,4) = du - dqn*nx2
  R(3,4) = dv - dqn*ny2
  R(4,4) = dw - dqn*nz2
  R(5,4) = u*du + v*dv + w*dw - qn*dqn
  !
  ! Dissipation Term: Roe dissipation with the modified wave speeds.
  !    |An|dU = R|Lambda|L*dU
  !           = sum_k of [ ws(k) * R(:,k) * L*dU(k) ], where n=n2.
  !
  diss(:) = ws(1)*LdU(1)*R(:,1) + &
            ws(2)*LdU(2)*R(:,2) + &
            ws(3)*LdU(3)*R(:,3) + &
            ws(4)*LdU(4)*R(:,4)
  !
  ! Compute the Rotated-RHLL flux, it looks like HLL flux with Roe dissipation
  !
  return_value = (SRp*fL - SLm*fR)/(SRp-SLm) - half*diss
  !
  ! Normal max wave speed in the normal direction.
  !  wsn = abs(qn) + a
  !
end subroutine rotated_rhll_upwind_flux_sub_3d
!
!###############################################################################
!
