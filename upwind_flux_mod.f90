module upwind_flux_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  use eqn_idx, only : ncs,npres,ngamm,nenth
  !
  implicit none
  !
  private
  !
  public :: upwind_flux_gpu
  !
  public :: upwind_flux_cv
  public :: hllc_upwind_flux_cv
  public :: ldfss_upwind_flux_cv
 !public :: ausm_plus_upwind_flux_cv
  public :: rotated_rhll_upwind_flux_cv_2d
  public :: rotated_rhll_upwind_flux_cv_3d
  !
  !.. Module Parameters ..
  real(wp), parameter :: upaefix = zero
 !real(wp), parameter :: upaefix = -half
  !
contains
!
!###############################################################################
!
!pure &
subroutine upwind_flux_gpu(n1,n2,facenrm,face_cv,face_gp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : invs_flux_method
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,n,np
  real(wp) :: mfl,mfr
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: lflux(:)
  real(wp), allocatable :: rflux(:)
  !
continue
  !
  np = size(face_cv,dim=1)
  !
  allocate ( lflux(1:np) , source=zero )
  allocate ( rflux(1:np) , source=zero )
  !$ACC DATA CREATE(lflux,rflux)
  !
  if (invs_flux_method == invs_flux_ROE) then
    !
    call roe_upwind_flux_gpu(n1,n2,facenrm,face_cv,face_gp, &
                             lflux,rflux,return_value)
    !
  else if (invs_flux_method == invs_flux_HLLC) then
    !
    call hllc_upwind_flux_gpu(n1,n2,facenrm,face_cv,face_gp, &
                              lflux,rflux,return_value)
    !
  else if (invs_flux_method == invs_flux_LDFSS) then
    !
    call ldfss_upwind_flux_gpu(n1,n2,facenrm,face_cv,face_gp, &
                               lflux,rflux,return_value)
    !
  else if (invs_flux_method == invs_flux_AUSM_Plus) then
    !
    call ausm_plus_upwind_flux_gpu(n1,n2,facenrm,face_cv,face_gp, &
                                   lflux,rflux,return_value)
    !
  else if (invs_flux_method == invs_flux_Rotated_RHLL) then
    !
    call rotated_rhll_upwind_flux_gpu(n1,n2,facenrm,face_cv,face_gp, &
                                      lflux,rflux,return_value)
    !
  end if
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv,lflux,rflux,return_value)
    do i = 1,ncs
      do n = n1,n2
        mfl = face_cv(n,i,1)/face_cv(n,nec,1)
        mfr = face_cv(n,i,2)/face_cv(n,nec,2)
        return_value(n,i) = mfl*lflux(n) + mfr*rflux(n)
      end do
    end do
  end if
  !
  !$ACC END DATA
  deallocate ( lflux )
  deallocate ( rflux )
  !
end subroutine upwind_flux_gpu
!
!###############################################################################
!
!pure &
subroutine roe_upwind_flux_gpu(n1,n2,facenrm,face_cv,face_gp, &
                                    lflux,rflux,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  use eqn_idx, only : aux_gamma,aux_enthalpy
  use eqn_idx, only : aux_var_is_not_gamma
  use eqn_idx, only : aux_var_is_not_enthalpy
  use eqn_idx, only : set_auxiliary_variable
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,n,nr,np
  real(wp) :: gam_over_gm1
  !
  !.. Local Parameters ..
 !logical(lk), parameter :: use_new_roe_subs = true
  logical(lk), parameter :: use_new_roe_subs = fals
  !
continue
  !
  nr = size(facenrm,dim=1)
  np = size(face_cv,dim=1)
  !
  ! Store the static enthalpy for the left and right states in the lflux and
  ! rflux arrays, respectively. After returning from roe_upwind_flux_tpg_*,
  ! the gammas stored in lflux and rflux will be replaced with the left and
  ! right fluxes for the species continuity equations.
  !
  if (use_new_roe_subs .or. ncs > 0) then
    !
    ! Move enthalpy into the auxiliary gas property variable spot
    !
    if ( aux_var_is_not_enthalpy() ) then
      gam_over_gm1 = gam / (gam - one)
      !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv,face_gp)
      do i = 1,2
        do n = n1,n2
          face_gp(n,nenth,i) = gam_over_gm1*face_gp(n,npres,i)/face_cv(n,nec,i)
        end do
      end do
      call set_auxiliary_variable(aux_enthalpy)
    end if
    !
    ! Copy the static enthalpy into lflux and rflux
    !
    !$ACC PARALLEL LOOP PRESENT(face_gp,lflux,rflux)
    do n = n1,n2
      lflux(n) = face_gp(n,nenth,1)
      rflux(n) = face_gp(n,nenth,2)
    end do
    !
    ! Move gamma into the auxiliary gas property variable spot
    !
    if ( aux_var_is_not_gamma() ) then
      !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_gp)
      do i = 1,2
        do n = n1,n2
          face_gp(n,ngamm,i) = gam
        end do
      end do
      call set_auxiliary_variable(aux_gamma)
    end if
    !
  end if
  !
  if (use_new_roe_subs) then
    !
    ! Use the new Roe upwinding subroutines regardless of whether or not
    ! this is a calorically perfect gas or a thermally perfect gas mixture.
    !
    if (nr == 3) then
      call roe_upwind_flux_tpg_3d(n1,n2,facenrm,face_cv,face_gp, &
                                  lflux,rflux,return_value)
    else if (nr == 2) then
      call roe_upwind_flux_tpg_2d(n1,n2,facenrm,face_cv,face_gp, &
                                  lflux,rflux,return_value)
    else
      call roe_upwind_flux_tpg_1d(n1,n2,facenrm,face_cv,face_gp, &
                                  lflux,rflux,return_value)
    end if
    !
  else
    !
    if (ncs > 0) then
      !
      ! Use the new Roe upwinding subroutines only
      ! if this is a thermally perfect gas mixture.
      !
      if (nr == 3) then
        call roe_upwind_flux_tpg_3d(n1,n2,facenrm,face_cv,face_gp, &
                                    lflux,rflux,return_value)
      else if (nr == 2) then
        call roe_upwind_flux_tpg_2d(n1,n2,facenrm,face_cv,face_gp, &
                                    lflux,rflux,return_value)
      else
        call roe_upwind_flux_tpg_1d(n1,n2,facenrm,face_cv,face_gp, &
                                    lflux,rflux,return_value)
      end if
      !
    else
      !
      ! Use the old Roe upwinding subroutines only
      ! if this is a calorically perfect gas.
      !
      if (nr == 3) then
        call roe_upwind_flux_cpg_3d(n1,n2,facenrm,face_cv,face_gp,return_value)
      else if (nr == 2) then
        call roe_upwind_flux_cpg_2d(n1,n2,facenrm,face_cv,face_gp,return_value)
      else
        call roe_upwind_flux_cpg_1d(n1,n2,facenrm,face_cv,face_gp,return_value)
      end if
      !
    end if
    !
  end if
  !
end subroutine roe_upwind_flux_gpu
!
!###############################################################################
!
!pure &
subroutine hllc_upwind_flux_gpu(n1,n2,facenrm,face_cv,face_gp, &
                                     lflux,rflux,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  use eqn_idx, only : aux_gamma
  use eqn_idx, only : aux_var_is_not_gamma
  use eqn_idx, only : set_auxiliary_variable
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,n,nr
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  if ( aux_var_is_not_gamma() ) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_gp)
    do i = 1,2
      do n = n1,n2
        face_gp(n,ngamm,i) = gam
      end do
    end do
    call set_auxiliary_variable(aux_gamma)
  end if
  !
  if (nr == 3) then
    call hllc_upwind_flux_gpu_3d(n1,n2,facenrm,face_cv,face_gp, &
                                 lflux,rflux,return_value)
  else if (nr == 2) then
    call hllc_upwind_flux_gpu_2d(n1,n2,facenrm,face_cv,face_gp, &
                                 lflux,rflux,return_value)
  else
    call hllc_upwind_flux_gpu_1d(n1,n2,facenrm,face_cv,face_gp, &
                                 lflux,rflux,return_value)
  end if
  !
end subroutine hllc_upwind_flux_gpu
!
!###############################################################################
!
!pure &
subroutine ldfss_upwind_flux_gpu(n1,n2,facenrm,face_cv,face_gp, &
                                      lflux,rflux,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  use eqn_idx, only : aux_gamma,aux_enthalpy
  use eqn_idx, only : aux_var_is_not_gamma
  use eqn_idx, only : aux_var_is_not_enthalpy
  use eqn_idx, only : set_auxiliary_variable
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp), intent(inout) :: face_gp(:,:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,n,np,nr
  real(wp) :: rho_inv,al,ar
  real(wp) :: gam_over_gm1,pres,gamm
  !
continue
  !
  nr = size(facenrm,dim=1)
  np = size(face_cv,dim=1)
  !
  ! Move gamma into the auxiliary gas property variable spot
  !
  if ( aux_var_is_not_gamma() ) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_gp)
    do i = 1,2
      do n = n1,n2
        face_gp(n,ngamm,i) = gam
      end do
    end do
    call set_auxiliary_variable(aux_gamma)
  end if
  !
  ! Compute the common speed of sound and store in lflux
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp,lflux)
  do n = n1,n2
    !
    rho_inv = one/face_cv(n,nec,1)
    pres = face_gp(n,npres,1)
    gamm = face_gp(n,ngamm,1)
    al = sqrt(max(gamm*pres*rho_inv,rndoff))
    !
    rho_inv = one/face_cv(n,nec,2)
    pres = face_gp(n,npres,2)
    gamm = face_gp(n,ngamm,2)
    ar = sqrt(max(gamm*pres*rho_inv,rndoff))
    !
    lflux(n) = half*(al + ar)
    !
  end do
  !
  ! Move enthalpy into the auxiliary gas property variable spot
  !
  if ( aux_var_is_not_enthalpy() ) then
    gam_over_gm1 = gam / (gam - one)
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_gp)
    do i = 1,2
      do n = n1,n2
        face_gp(n,nenth,i) = gam_over_gm1*face_gp(n,npres,i)/face_cv(n,nec,i)
      end do
    end do
    call set_auxiliary_variable(aux_enthalpy)
  end if
  !
  ! Compute the upwind flux using LDFSS
  !
  if (nr == 3) then
    call ldfss_upwind_flux_gpu_3d(n1,n2,facenrm,face_cv,face_gp, &
                                  lflux,rflux,return_value)
  else if (nr == 2) then
    call ldfss_upwind_flux_gpu_2d(n1,n2,facenrm,face_cv,face_gp, &
                                  lflux,rflux,return_value)
  else
    call ldfss_upwind_flux_gpu_1d(n1,n2,facenrm,face_cv,face_gp, &
                                  lflux,rflux,return_value)
  end if
  !
end subroutine ldfss_upwind_flux_gpu
!
!###############################################################################
!
pure subroutine ausm_plus_upwind_flux_gpu(n1,n2,facenrm,face_cv,face_gp, &
                                          lflux,rflux,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  if (nr == 3) then
    call ausm_plus_upwind_flux_gpu_3d(n1,n2,facenrm,face_cv,face_gp, &
                                      lflux,rflux,return_value)
  else if (nr == 2) then
    call ausm_plus_upwind_flux_gpu_2d(n1,n2,facenrm,face_cv,face_gp, &
                                      lflux,rflux,return_value)
  else
    call ausm_plus_upwind_flux_gpu_nd(n1,n2,facenrm,face_cv,face_gp, &
                                      lflux,rflux,return_value)
  end if
  !
end subroutine ausm_plus_upwind_flux_gpu
!
!###############################################################################
!
pure subroutine rotated_rhll_upwind_flux_gpu(n1,n2,facenrm,face_cv,face_gp, &
                                             lflux,rflux,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  if (nr == 3) then
    call rotated_rhll_upwind_flux_gpu_3d(n1,n2,facenrm,face_cv,face_gp, &
                                         lflux,rflux,return_value)
  else if (nr == 2) then
    call rotated_rhll_upwind_flux_gpu_2d(n1,n2,facenrm,face_cv,face_gp, &
                                         lflux,rflux,return_value)
  end if
  !
end subroutine rotated_rhll_upwind_flux_gpu
!
!###############################################################################
!
pure subroutine roe_upwind_flux_cpg_3d(n1,n2,facenrm,face_cv,face_gp, &
                                       return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gp1,gam_over_gm1,ds
  real(wp) :: rl,presl,h0l,q2l,vnl,betal
  real(wp) :: rr,presr,h0r,q2r,vnr,betar
  real(wp) :: rl_inv,vxl,vyl,vzl
  real(wp) :: rr_inv,vxr,vyr,vzr
  real(wp) :: ra_mx,ra_my,ra_mz
  real(wp) :: unx,uny,unz
  real(wp) :: rnx,rny,rnz
  real(wp) :: dw1,hsr1,eta1,udw1,rvnl
  real(wp) :: dw4,hsr4,eta4,udw4,rvnr
  real(wp) :: ra_ec,ra_h0,ra_a2,ra_a,ra_q2,ra_vn
!#define PROTECT_AGAINST_ATD_EQUAL_ZERO
#ifdef PROTECT_AGAINST_ATD_EQUAL_ZERO
  real(wp) :: atd2nz,atd2rt,atdnz,atdrt
#endif
  !
continue
  !
  gm1 = gam - one
  gp1 = gam + one
  !
  gam_over_gm1 = gam / gm1
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,return_value)
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
    q2l = vxl*vxl + vyl*vyl + vzl*vzl
    presl = face_gp(n,npres,1)
    h0l = half*q2l + gam_over_gm1*presl*rl_inv
    !
    rr = face_cv(n,nec,2)
    rr_inv = one/rr
    vxr = face_cv(n,nmx,2)*rr_inv
    vyr = face_cv(n,nmy,2)*rr_inv
    vzr = face_cv(n,nmz,2)*rr_inv
    q2r = vxr*vxr + vyr*vyr + vzr*vzr
    presr = face_gp(n,npres,2)
    h0r = half*q2r + gam_over_gm1*presr*rr_inv
    !
    vnl = vxl*unx + vyl*uny + vzl*unz
    vnr = vxr*unx + vyr*uny + vzr*unz
    !
    ! Roe averaged variables
    !
    ra_ec = sqrt(rr*rl)
    betal = rl / (rl+ra_ec)
    betar = one - betal
    ra_mx = betal*vxl + betar*vxr
    ra_my = betal*vyl + betar*vyr
    ra_mz = betal*vzl + betar*vzr
    ra_h0 = betal*h0l + betar*h0r
    ra_q2 = ra_mx*ra_mx + ra_my*ra_my + ra_mz*ra_mz
    ra_a2 = max(gm1*(ra_h0 - half*ra_q2),zero)
    ra_a  = sqrt(ra_a2)
    ra_vn = ra_mx*unx + ra_my*uny + ra_mz*unz
    !
#ifdef PROTECT_AGAINST_ATD_EQUAL_ZERO
    atd2nz = max(ra_a2,rndoff)
    atdnz = max(ra_a,rndoff)
    atd2rt = ra_a2 / atd2nz
    atdrt = ra_a / atdnz
#endif
    !
    if (ra_vn >= zero) then
      !
#ifdef PROTECT_AGAINST_ATD_EQUAL_ZERO
      dw1  = half*((presr-presl)/atd2nz*atd2rt - (vnr-vnl)*ra_ec/atdnz*atdrt)
#else
      dw1  = half*((presr-presl)/ra_a2 - (vnr-vnl)*ra_ec/ra_a)
#endif
      hsr1 = -one4*gp1*dw1*ra_a/ra_ec
      eta1 = max(-abs(ra_vn-ra_a) + hsr1, zero)
      udw1 = dw1*(min(ra_vn-ra_a,zero) - half*eta1)
      rvnl = rl*vnl
      !
      return_value(n,nec) = ds*(    rvnl + udw1                             )
      return_value(n,nmx) = ds*(vxl*rvnl + udw1*(ra_mx-ra_a*unx) + presl*unx)
      return_value(n,nmy) = ds*(vyl*rvnl + udw1*(ra_my-ra_a*uny) + presl*uny)
      return_value(n,nmz) = ds*(vzl*rvnl + udw1*(ra_mz-ra_a*unz) + presl*unz)
      return_value(n,nee) = ds*(h0l*rvnl + udw1*(ra_h0-ra_a*ra_vn)          )
      !
    else
      !
#ifdef PROTECT_AGAINST_ATD_EQUAL_ZERO
      dw4  = half*((presr-presl)/atd2nz*atd2rt + (vnr-vnl)*ra_ec/atdnz*atdrt)
#else
      dw4  = half*((presr-presl)/ra_a2 + (vnr-vnl)*ra_ec/ra_a)
#endif
      hsr4 = one4*gp1*dw4*ra_a/ra_ec
      eta4 = max(-abs(ra_vn+ra_a) + hsr4, zero)
      udw4 = dw4*(max(ra_vn+ra_a,zero) + half*eta4)
      rvnr = rr*vnr
      !
      return_value(n,nec) = ds*(    rvnr - udw4                             )
      return_value(n,nmx) = ds*(vxr*rvnr - udw4*(ra_mx+ra_a*unx) + presr*unx)
      return_value(n,nmy) = ds*(vyr*rvnr - udw4*(ra_my+ra_a*uny) + presr*uny)
      return_value(n,nmz) = ds*(vzr*rvnr - udw4*(ra_mz+ra_a*unz) + presr*unz)
      return_value(n,nee) = ds*(h0r*rvnr - udw4*(ra_h0+ra_a*ra_vn)          )
      !
    end if
    !
  end do
  !
end subroutine roe_upwind_flux_cpg_3d
!
!###############################################################################
!
!pure &
subroutine roe_upwind_flux_cpg_2d(n1,n2,facenrm,face_cv,face_gp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gp1,gam_over_gm1,ds
  real(wp) :: rl,presl,h0l,q2l,vnl,betal
  real(wp) :: rr,presr,h0r,q2r,vnr,betar
  real(wp) :: rl_inv,vxl,vyl
  real(wp) :: rr_inv,vxr,vyr
  real(wp) :: ra_mx,ra_my
  real(wp) :: unx,uny
  real(wp) :: rnx,rny
  real(wp) :: dw1,hsr1,eta1,udw1,rvnl
  real(wp) :: dw4,hsr4,eta4,udw4,rvnr
  real(wp) :: ra_ec,ra_h0,ra_a2,ra_a,ra_q2,ra_vn
#ifdef PROTECT_AGAINST_ATD_EQUAL_ZERO
  real(wp) :: atd2nz,atd2rt,atdnz,atdrt
#endif
  !
continue
  !
  gm1 = gam - one
  gp1 = gam + one
  !
  gam_over_gm1 = gam / gm1
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,return_value)
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
    q2l = vxl*vxl + vyl*vyl
    presl = face_gp(n,npres,1)
    h0l = half*q2l + gam_over_gm1*presl*rl_inv
    !
    rr = face_cv(n,nec,2)
    rr_inv = one/rr
    vxr = face_cv(n,nmx,2)*rr_inv
    vyr = face_cv(n,nmy,2)*rr_inv
    q2r = vxr*vxr + vyr*vyr
    presr = face_gp(n,npres,2)
    h0r = half*q2r + gam_over_gm1*presr*rr_inv
    !
    vnl = vxl*unx + vyl*uny
    vnr = vxr*unx + vyr*uny
    !
    ! Roe averaged variables
    !
    ra_ec = sqrt(rr*rl)
    betal = rl / (rl+ra_ec)
    betar = one - betal
    ra_mx = betal*vxl + betar*vxr
    ra_my = betal*vyl + betar*vyr
    ra_h0 = betal*h0l + betar*h0r
    ra_q2 = ra_mx*ra_mx + ra_my*ra_my
    ra_a2 = max(gm1*(ra_h0 - half*ra_q2),zero)
    ra_a  = sqrt(ra_a2)
    ra_vn = ra_mx*unx + ra_my*uny
    !
#ifdef PROTECT_AGAINST_ATD_EQUAL_ZERO
    atd2nz = max(ra_a2,rndoff)
    atdnz = max(ra_a,rndoff)
    atd2rt = ra_a2 / atd2nz
    atdrt = ra_a / atdnz
#endif
    !
    if (ra_vn >= zero) then
      !
#ifdef PROTECT_AGAINST_ATD_EQUAL_ZERO
      dw1  = half*((presr-presl)/atd2nz*atd2rt - (vnr-vnl)*ra_ec/atdnz*atdrt)
#else
      dw1  = half*((presr-presl)/ra_a2 - (vnr-vnl)*ra_ec/ra_a)
#endif
      hsr1 = -one4*gp1*dw1*ra_a/ra_ec
      eta1 = max(-abs(ra_vn-ra_a) + hsr1, zero)
      udw1 = dw1*(min(ra_vn-ra_a,zero) - half*eta1)
      rvnl = rl*vnl
      !
      return_value(n,nec) = ds*(    rvnl + udw1                             )
      return_value(n,nmx) = ds*(vxl*rvnl + udw1*(ra_mx-ra_a*unx) + presl*unx)
      return_value(n,nmy) = ds*(vyl*rvnl + udw1*(ra_my-ra_a*uny) + presl*uny)
      return_value(n,nee) = ds*(h0l*rvnl + udw1*(ra_h0-ra_a*ra_vn)          )
      !
    else
      !
#ifdef PROTECT_AGAINST_ATD_EQUAL_ZERO
      dw4  = half*((presr-presl)/atd2nz*atd2rt + (vnr-vnl)*ra_ec/atdnz*atdrt)
#else
      dw4  = half*((presr-presl)/ra_a2 + (vnr-vnl)*ra_ec/ra_a)
#endif
      hsr4 = one4*gp1*dw4*ra_a/ra_ec
      eta4 = max(-abs(ra_vn+ra_a) + hsr4, zero)
      udw4 = dw4*(max(ra_vn+ra_a,zero) + half*eta4)
      rvnr = rr*vnr
      !
      return_value(n,nec) = ds*(    rvnr - udw4                             )
      return_value(n,nmx) = ds*(vxr*rvnr - udw4*(ra_mx+ra_a*unx) + presr*unx)
      return_value(n,nmy) = ds*(vyr*rvnr - udw4*(ra_my+ra_a*uny) + presr*uny)
      return_value(n,nee) = ds*(h0r*rvnr - udw4*(ra_h0+ra_a*ra_vn)          )
      !
    end if
!#define DEBUG_UPWIND_FLUXES
#ifdef DEBUG_UPWIND_FLUXES
    block
      use ovar, only : itcur
      write (200+itcur,1) n,return_value(n,nec),return_value(n,nmx), &
                            return_value(n,nmy),return_value(n,nee)
      write (210+itcur,1) n,rl,rl_inv,vxl,vyl,q2l,presl,h0l,vnl
      write (220+itcur,1) n,rr,rr_inv,vxr,vyr,q2r,presr,h0r,vnr
      write (230+itcur,1) n,ra_ec,betal,betar,ra_mx,ra_my, &
                            ra_h0,ra_q2,ra_a2,ra_a,ra_vn
#ifdef PROTECT_AGAINST_ATD_EQUAL_ZERO
      write (240+itcur,1) n,atd2nz,atdnz,atd2rt,atdrt
#endif
      if (ra_vn >= zero) then
        write (250+itcur,2) n,"ra_vn >= 0",dw1,hsr1,eta1,udw1,rvnl
      else
        write (250+itcur,2) n,"ra_vn <  0",dw4,hsr4,eta4,udw4,rvnr
      end if
      write (260+itcur,3) n,presr,presl,presr-presl,vnr,vnl,vnr-vnl
      1 format (i6,2x,10es18.10)
      2 format (i6,2x,a,2x,5es18.10)
      3 format (i6,2x,6es23.15)
    end block
#endif
    !
  end do
  !
end subroutine roe_upwind_flux_cpg_2d
!
!###############################################################################
!
pure subroutine roe_upwind_flux_cpg_1d(n1,n2,facenrm,face_cv,face_gp, &
                                       return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,gp1,gam_over_gm1,ds
  real(wp) :: rl,presl,h0l,q2l,vnl,betal
  real(wp) :: rr,presr,h0r,q2r,vnr,betar
  real(wp) :: rl_inv,vxl
  real(wp) :: rr_inv,vxr
  real(wp) :: ra_mx
  real(wp) :: unx
  real(wp) :: rnx
  real(wp) :: dw1,hsr1,eta1,udw1,rvnl
  real(wp) :: dw4,hsr4,eta4,udw4,rvnr
  real(wp) :: ra_ec,ra_h0,ra_a2,ra_a,ra_q2,ra_vn
#ifdef PROTECT_AGAINST_ATD_EQUAL_ZERO
  real(wp) :: atd2nz,atd2rt,atdnz,atdrt
#endif
  !
continue
  !
  gm1 = gam - one
  gp1 = gam + one
  !
  gam_over_gm1 = gam / gm1
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,return_value)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    ds = sqrt( rnx*rnx )
    unx = rnx/ds
    !
    rl = face_cv(n,nec,1)
    rl_inv = one/rl
    vxl = face_cv(n,nmx,1)*rl_inv
    q2l = vxl*vxl
    presl = face_gp(n,npres,1)
    h0l = half*q2l + gam_over_gm1*presl*rl_inv
    !
    rr = face_cv(n,nec,2)
    rr_inv = one/rr
    vxr = face_cv(n,nmx,2)*rr_inv
    q2r = vxr*vxr
    presr = face_gp(n,npres,2)
    h0r = half*q2r + gam_over_gm1*presr*rr_inv
    !
    vnl = vxl*unx
    vnr = vxr*unx
    !
    ! Roe averaged variables
    !
    ra_ec = sqrt(rr*rl)
    betal = rl / (rl+ra_ec)
    betar = one - betal
    ra_mx = betal*vxl + betar*vxr
    ra_h0 = betal*h0l + betar*h0r
    ra_q2 = ra_mx*ra_mx
    ra_a2 = max(gm1*(ra_h0 - half*ra_q2),zero)
    ra_a  = sqrt(ra_a2)
    ra_vn = ra_mx*unx
    !
#ifdef PROTECT_AGAINST_ATD_EQUAL_ZERO
    atd2nz = max(ra_a2,rndoff)
    atdnz = max(ra_a,rndoff)
    atd2rt = ra_a2 / atd2nz
    atdrt = ra_a / atdnz
#endif
    !
    if (ra_vn >= zero) then
      !
#ifdef PROTECT_AGAINST_ATD_EQUAL_ZERO
      dw1  = half*((presr-presl)/atd2nz*atd2rt - (vnr-vnl)*ra_ec/atdnz*atdrt)
#else
      dw1  = half*((presr-presl)/ra_a2 - (vnr-vnl)*ra_ec/ra_a)
#endif
      hsr1 = -one4*gp1*dw1*ra_a/ra_ec
      eta1 = max(-abs(ra_vn-ra_a) + hsr1, zero)
      udw1 = dw1*(min(ra_vn-ra_a,zero) - half*eta1)
      rvnl = rl*vnl
      !
      return_value(n,nec) = ds*(    rvnl + udw1                             )
      return_value(n,nmx) = ds*(vxl*rvnl + udw1*(ra_mx-ra_a*unx) + presl*unx)
      return_value(n,nee) = ds*(h0l*rvnl + udw1*(ra_h0-ra_a*ra_vn)          )
      !
    else
      !
#ifdef PROTECT_AGAINST_ATD_EQUAL_ZERO
      dw4  = half*((presr-presl)/atd2nz*atd2rt + (vnr-vnl)*ra_ec/atdnz*atdrt)
#else
      dw4  = half*((presr-presl)/ra_a2 + (vnr-vnl)*ra_ec/ra_a)
#endif
      hsr4 = one4*gp1*dw4*ra_a/ra_ec
      eta4 = max(-abs(ra_vn+ra_a) + hsr4, zero)
      udw4 = dw4*(max(ra_vn+ra_a,zero) + half*eta4)
      rvnr = rr*vnr
      !
      return_value(n,nec) = ds*(    rvnr - udw4                             )
      return_value(n,nmx) = ds*(vxr*rvnr - udw4*(ra_mx+ra_a*unx) + presr*unx)
      return_value(n,nee) = ds*(h0r*rvnr - udw4*(ra_h0+ra_a*ra_vn)          )
      !
    end if
    !
  end do
  !
end subroutine roe_upwind_flux_cpg_1d
!
!###############################################################################
!
pure subroutine roe_upwind_flux_tpg_3d(n1,n2,facenrm,face_cv,face_gp, &
                                       lflux,rflux,return_value)
  !
  ! Routine to calculate the upwind flux difference splitting using the
  ! approximate Riemann solver of Roe for a mixture of gases in 3D
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rl,rl_inv,vxl,vyl,vzl,etotl,presl,hl,gaml,vnl,h0l,eil
  real(wp) :: rr,rr_inv,vxr,vyr,vzr,etotr,presr,hr,gamr,vnr,h0r,eir
  real(wp) :: drho,dvx,dvy,dvz,dvn,dpres,dq2
  real(wp) :: ra_weight_l,ra_weight_r
  real(wp) :: ra_ec,ra_mx,ra_my,ra_mz,ra_vn,ra_pres
  real(wp) :: ra_gam,ra_gm1,ra_a2,ra_a,ra_h0
  real(wp) :: alf1,alf2,alf3,bet1,bet2,bet3
  real(wp) :: b1,b2,bx,by,bz,b6,b7,b8
  real(wp) :: fec_l,fmx_l,fmy_l,fmz_l,fee_l
  real(wp) :: fec_r,fmx_r,fmy_r,fmz_r,fee_r
  real(wp) :: fec_d,fmx_d,fmy_d,fmz_d,fee_d
  !
continue
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp) &
  !$ACC&              PRESENT(lflux,rflux,return_value)
  do n = n1,n2
    !
    ! Cell face area
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    ! Left(-) state variables
    !
    rl = face_cv(n,nec,1)
    rl_inv = one/rl
    vxl = face_cv(n,nmx,1)*rl_inv
    vyl = face_cv(n,nmy,1)*rl_inv
    vzl = face_cv(n,nmz,1)*rl_inv
    etotl = face_cv(n,nee,1)
    presl = face_gp(n,npres,1)
    gaml = face_gp(n,ngamm,1)
    hl = lflux(n)
    !
    vnl = vxl*unx + vyl*uny + vzl*unz
    h0l = hl + half*(vxl*vxl+vyl*vyl+vzl*vzl)
    eil = hl - presl*rl_inv
    !
    ! Right(+) state variables
    !
    rr = face_cv(n,nec,2)
    rr_inv = one/rr
    vxr = face_cv(n,nmx,2)*rr_inv
    vyr = face_cv(n,nmy,2)*rr_inv
    vzr = face_cv(n,nmz,2)*rr_inv
    etotr = face_cv(n,nee,2)
    presr = face_gp(n,npres,2)
    gamr = face_gp(n,ngamm,2)
    hr = rflux(n)
    !
    vnr = vxr*unx + vyr*uny + vzr*unz
    h0r = hr + half*(vxr*vxr+vyr*vyr+vzr*vzr)
    eir = hr - presr*rr_inv
    !
    ! Delta state variables
    !
    drho = rr - rl
    dvx = vxr - vxl
    dvy = vyr - vyl
    dvz = vzr - vzl
    dvn = vnr - vnl
    dpres = presr - presl
    dq2 = dvx*dvx + dvy*dvy + dvz*dvz
    !
    ! Roe averaged density
    !
    ra_ec = sqrt(rr*rl)
    !
    ! Roe averaging weights
    !
    ra_weight_l = rl / (rl+ra_ec)
    ra_weight_r = one - ra_weight_l
    !
    ! Roe averaged velocities
    !
    ra_mx = ra_weight_r*vxr + vxl*ra_weight_l
    ra_my = ra_weight_r*vyr + vyl*ra_weight_l
    ra_mz = ra_weight_r*vzr + vzl*ra_weight_l
    !
    ! Roe averaged normal velocity
    !
    ra_vn = ra_mx*unx + ra_my*uny + ra_mz*unz
    !
    ! Roe averaged pressure (special averaging for pressure used here)
    !
    ra_pres = ra_weight_r*presl + presr*ra_weight_l
    !
    ! Roe averaged gamma
    !
    ra_gam = ra_weight_r*gamr + gaml*ra_weight_l
    ra_gm1 = ra_gam - one
    !
    ! Roe averaged speed of sound
    !
    ra_a2 = ra_gam*ra_pres/ra_ec + &
            half*ra_ec*ra_gm1*dq2 / ((sqrt(rl)+sqrt(rr))**2)
    ra_a = sqrt(ra_a2)
    !
    ! Roe averaged total enthalpy
    !
    ra_h0 = ra_weight_r*h0r + h0l*ra_weight_l
    !
    ! Compute alphas
    !
    alf1 = half * (dpres/ra_a2 + ra_ec*dvn/ra_a)
    alf2 = half * (dpres/ra_a2 - ra_ec*dvn/ra_a)
    alf3 = drho - dpres/ra_a2
    !
    ! Compute betas
    !
    bet1 = alf1 * abs(ra_vn + ra_a)
    bet2 = alf2 * abs(ra_vn - ra_a)
    bet3 = alf3 * abs(ra_vn)
    !
    ! Compute Bi
    !
    b1 = bet1 + bet2 + bet3
    b2 = ra_a * (bet1 - bet2)
    bx = ra_ec * (dvx - unx*dvn) * abs(ra_vn)
    by = ra_ec * (dvy - uny*dvn) * abs(ra_vn)
    bz = ra_ec * (dvz - unz*dvn) * abs(ra_vn)
    b6 = bet3 * ra_a2 / ra_gm1
    b7 = ra_mx*bx + ra_my*by + ra_mz*bz
    ! NOTE: For a calorically perfect gas, gaml=gamr=ra_gam=gam, and
    !       I did the math and b8 should equal zero (at least analytically)
    b8 = (ra_ec*(eir-eil) + (ra_pres/ra_ec*drho-dpres)/ra_gm1) * abs(ra_vn)
    !
    ! Compute the dissipative flux = 1/2*sum(ALFi*RALMi*right eigenvector)
    !
    fec_d = b1
    fmx_d = b1*ra_mx + b2*unx + bx
    fmy_d = b1*ra_my + b2*uny + by
    fmz_d = b1*ra_mz + b2*unz + bz
    fee_d = b1*ra_h0 + b2*ra_vn - b6 + b7 + b8
    !
    ! Compute the left(-) convective fluxes
    !
    fec_l = vnl*rl
    fmx_l = vxl*fec_l + presl*unx
    fmy_l = vyl*fec_l + presl*uny
    fmz_l = vzl*fec_l + presl*unz
    fee_l = vnl*etotl + presl*vnl
    !
    ! Compute the right(+) convective fluxes
    !
    fec_r = vnr*rr
    fmx_r = vxr*fec_r + presr*unx
    fmy_r = vyr*fec_r + presr*uny
    fmz_r = vzr*fec_r + presr*unz
    fee_r = vnr*etotr + presr*vnr
    !
    ! Compute the total upwind flux from the Roe scheme
    !
    return_value(n,nec) = half * ds * (fec_l + fec_r - fec_d)
    return_value(n,nmx) = half * ds * (fmx_l + fmx_r - fmx_d)
    return_value(n,nmy) = half * ds * (fmy_l + fmy_r - fmy_d)
    return_value(n,nmz) = half * ds * (fmz_l + fmz_r - fmz_d)
    return_value(n,nee) = half * ds * (fee_l + fee_r - fee_d)
    !
    ! Save the left and right flux for species continuity equations
    !
    lflux(n) = half*ds*(fec_l + ra_ec*abs(ra_vn) - b1*ra_weight_l)
    rflux(n) = half*ds*(fec_r - ra_ec*abs(ra_vn) - b1*ra_weight_r)
    !
  end do
  !
end subroutine roe_upwind_flux_tpg_3d
!
!###############################################################################
!
pure subroutine roe_upwind_flux_tpg_2d(n1,n2,facenrm,face_cv,face_gp, &
                                       lflux,rflux,return_value)
  !
  ! Routine to calculate the upwind flux difference splitting using the
  ! approximate Riemann solver of Roe for a mixture of gases in 2D
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: rl,rl_inv,vxl,vyl,etotl,presl,hl,gaml,vnl,h0l,eil
  real(wp) :: rr,rr_inv,vxr,vyr,etotr,presr,hr,gamr,vnr,h0r,eir
  real(wp) :: drho,dvx,dvy,dvn,dpres,dq2
  real(wp) :: ra_weight_l,ra_weight_r
  real(wp) :: ra_ec,ra_mx,ra_my,ra_vn,ra_pres
  real(wp) :: ra_gam,ra_gm1,ra_a2,ra_a,ra_h0
  real(wp) :: alf1,alf2,alf3,bet1,bet2,bet3
  real(wp) :: b1,b2,bx,by,b6,b7,b8
  real(wp) :: fec_l,fmx_l,fmy_l,fee_l
  real(wp) :: fec_r,fmx_r,fmy_r,fee_r
  real(wp) :: fec_d,fmx_d,fmy_d,fee_d
  !
continue
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp) &
  !$ACC&              PRESENT(lflux,rflux,return_value)
  do n = n1,n2
    !
    ! Cell face area
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    ! Left(-) state variables
    !
    rl = face_cv(n,nec,1)
    rl_inv = one/rl
    vxl = face_cv(n,nmx,1)*rl_inv
    vyl = face_cv(n,nmy,1)*rl_inv
    etotl = face_cv(n,nee,1)
    presl = face_gp(n,npres,1)
    gaml = face_gp(n,ngamm,1)
    hl = lflux(n)
    !
    vnl = vxl*unx + vyl*uny
    h0l = hl + half*(vxl*vxl+vyl*vyl)
    eil = hl - presl*rl_inv
    !
    ! Right(+) state variables
    !
    rr = face_cv(n,nec,2)
    rr_inv = one/rr
    vxr = face_cv(n,nmx,2)*rr_inv
    vyr = face_cv(n,nmy,2)*rr_inv
    etotr = face_cv(n,nee,2)
    presr = face_gp(n,npres,2)
    gamr = face_gp(n,ngamm,2)
    hr = rflux(n)
    !
    vnr = vxr*unx + vyr*uny
    h0r = hr + half*(vxr*vxr+vyr*vyr)
    eir = hr - presr*rr_inv
    !
    ! Delta state variables
    !
    drho = rr - rl
    dvx = vxr - vxl
    dvy = vyr - vyl
    dvn = vnr - vnl
    dpres = presr - presl
    dq2 = dvx*dvx + dvy*dvy
    !
    ! Roe averaged density
    !
    ra_ec = sqrt(rr*rl)
    !
    ! Roe averaging weights
    !
    ra_weight_l = rl / (rl+ra_ec)
    ra_weight_r = one - ra_weight_l
    !
    ! Roe averaged velocities
    !
    ra_mx = ra_weight_r*vxr + vxl*ra_weight_l
    ra_my = ra_weight_r*vyr + vyl*ra_weight_l
    !
    ! Roe averaged normal velocity
    !
    ra_vn = ra_mx*unx + ra_my*uny
    !
    ! Roe averaged pressure (special averaging for pressure used here)
    !
    ra_pres = ra_weight_r*presl + presr*ra_weight_l
    !
    ! Roe averaged gamma
    !
    ra_gam = ra_weight_r*gamr + gaml*ra_weight_l
    ra_gm1 = ra_gam - one
    !
    ! Roe averaged speed of sound
    !
    ra_a2 = ra_gam*ra_pres/ra_ec + &
            half*ra_ec*ra_gm1*dq2 / ((sqrt(rl)+sqrt(rr))**2)
    ra_a = sqrt(ra_a2)
    !
    ! Roe averaged total enthalpy
    !
    ra_h0 = ra_weight_r*h0r + h0l*ra_weight_l
    !
    ! Compute alphas
    !
    alf1 = half * (dpres/ra_a2 + ra_ec*dvn/ra_a)
    alf2 = half * (dpres/ra_a2 - ra_ec*dvn/ra_a)
    alf3 = drho - dpres/ra_a2
    !
    ! Compute betas
    !
    bet1 = alf1 * abs(ra_vn + ra_a)
    bet2 = alf2 * abs(ra_vn - ra_a)
    bet3 = alf3 * abs(ra_vn)
    !
    ! Compute Bi
    !
    b1 = bet1 + bet2 + bet3
    b2 = ra_a * (bet1 - bet2)
    bx = ra_ec * (dvx - unx*dvn) * abs(ra_vn)
    by = ra_ec * (dvy - uny*dvn) * abs(ra_vn)
    b6 = bet3 * ra_a2 / ra_gm1
    b7 = ra_mx*bx + ra_my*by
    ! NOTE: For a calorically perfect gas, gaml=gamr=ra_gam=gam, and
    !       I did the math and b8 should equal zero (at least analytically)
    b8 = (ra_ec*(eir-eil) + (ra_pres/ra_ec*drho-dpres)/ra_gm1) * abs(ra_vn)
    !
    ! Compute the dissipative flux = 1/2*sum(ALFi*RALMi*right eigenvector)
    !
    fec_d = b1
    fmx_d = b1*ra_mx + b2*unx + bx
    fmy_d = b1*ra_my + b2*uny + by
    fee_d = b1*ra_h0 + b2*ra_vn - b6 + b7 + b8
    !
    ! Compute the left(-) convective fluxes
    !
    fec_l = vnl*rl
    fmx_l = vxl*fec_l + presl*unx
    fmy_l = vyl*fec_l + presl*uny
    fee_l = vnl*etotl + presl*vnl
    !
    ! Compute the right(+) convective fluxes
    !
    fec_r = vnr*rr
    fmx_r = vxr*fec_r + presr*unx
    fmy_r = vyr*fec_r + presr*uny
    fee_r = vnr*etotr + presr*vnr
    !
    ! Compute the total upwind flux from the Roe scheme
    !
    return_value(n,nec) = half * ds * (fec_l + fec_r - fec_d)
    return_value(n,nmx) = half * ds * (fmx_l + fmx_r - fmx_d)
    return_value(n,nmy) = half * ds * (fmy_l + fmy_r - fmy_d)
    return_value(n,nee) = half * ds * (fee_l + fee_r - fee_d)
    !
    ! Save the left and right flux for species continuity equations
    !
    lflux(n) = half*ds*(fec_l + ra_ec*abs(ra_vn) - b1*ra_weight_l)
    rflux(n) = half*ds*(fec_r - ra_ec*abs(ra_vn) - b1*ra_weight_r)
    !
  end do
  !
end subroutine roe_upwind_flux_tpg_2d
!
!###############################################################################
!
pure subroutine roe_upwind_flux_tpg_1d(n1,n2,facenrm,face_cv,face_gp, &
                                       lflux,rflux,return_value)
  !
  ! Routine to calculate the upwind flux difference splitting using the
  ! approximate Riemann solver of Roe for a mixture of gases in 1D
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1
  integer,     intent(in) :: n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,unx
  real(wp) :: rl,rl_inv,vxl,etotl,presl,hl,gaml,vnl,h0l,eil
  real(wp) :: rr,rr_inv,vxr,etotr,presr,hr,gamr,vnr,h0r,eir
  real(wp) :: drho,dvx,dvn,dpres
  real(wp) :: ra_weight_l,ra_weight_r
  real(wp) :: ra_ec,ra_mx,ra_vn,ra_pres
  real(wp) :: ra_gam,ra_gm1,ra_a2,ra_a,ra_h0
  real(wp) :: alf1,alf2,alf3,bet1,bet2,bet3
  real(wp) :: b1,b2,bx,b6,b7,b8
  real(wp) :: fec_l,fmx_l,fee_l
  real(wp) :: fec_r,fmx_r,fee_r
  real(wp) :: fec_d,fmx_d,fee_d
  !
continue
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp) &
  !$ACC&              PRESENT(lflux,rflux,return_value)
  do n = n1,n2
    !
    ! Cell face area
    !
    rnx = facenrm(1,n)
    ds = sqrt( rnx*rnx )
    unx = rnx/ds
    !
    ! Left(-) state variables
    !
    rl = face_cv(n,nec,1)
    rl_inv = one/rl
    vxl = face_cv(n,nmx,1)*rl_inv
    etotl = face_cv(n,nee,1)
    presl = face_gp(n,npres,1)
    gaml = face_gp(n,ngamm,1)
    hl = lflux(n)
    !
    vnl = vxl*unx
    h0l = hl + half*vxl*vxl
    eil = hl - presl*rl_inv
    !
    ! Right(+) state variables
    !
    rr = face_cv(n,nec,2)
    rr_inv = one/rr
    vxr = face_cv(n,nmx,2)*rr_inv
    etotr = face_cv(n,nee,2)
    presr = face_gp(n,npres,2)
    gamr = face_gp(n,ngamm,2)
    hr = rflux(n)
    !
    vnr = vxr*unx
    h0r = hr + half*vxr*vxr
    eir = hr - presr*rr_inv
    !
    ! Delta state variables
    !
    drho = rr - rl
    dvx = vxr - vxl
    dvn = vnr - vnl
    dpres = presr - presl
    !
    ! Roe averaged density
    !
    ra_ec = sqrt(rr*rl)
    !
    ! Roe averaging weights
    !
    ra_weight_l = rl / (rl+ra_ec)
    ra_weight_r = one - ra_weight_l
    !
    ! Roe averaged velocities
    !
    ra_mx = ra_weight_r*vxr + vxl*ra_weight_l
    !
    ! Roe averaged normal velocity
    !
    ra_vn = ra_mx*unx
    !
    ! Roe averaged pressure (special averaging for pressure used here)
    !
    ra_pres = ra_weight_r*presl + presr*ra_weight_l
    !
    ! Roe averaged gamma
    !
    ra_gam = ra_weight_r*gamr + gaml*ra_weight_l
    ra_gm1 = ra_gam - one
    !
    ! Roe averaged speed of sound
    !
    ra_a2 = ra_gam*ra_pres/ra_ec + &
            half*ra_ec*ra_gm1*dvx*dvx / ((sqrt(rl)+sqrt(rr))**2)
    ra_a = sqrt(ra_a2)
    !
    ! Roe averaged total enthalpy
    !
    ra_h0 = ra_weight_r*h0r + h0l*ra_weight_l
    !
    ! Compute alphas
    !
    alf1 = half * (dpres/ra_a2 + ra_ec*dvn/ra_a)
    alf2 = half * (dpres/ra_a2 - ra_ec*dvn/ra_a)
    alf3 = drho - dpres/ra_a2
    !
    ! Compute betas
    !
    bet1 = alf1 * abs(ra_vn + ra_a)
    bet2 = alf2 * abs(ra_vn - ra_a)
    bet3 = alf3 * abs(ra_vn)
    !
    ! Compute Bi
    !
    b1 = bet1 + bet2 + bet3
    b2 = ra_a * (bet1 - bet2)
    bx = ra_ec * (dvx - unx*dvn) * abs(ra_vn)
    b6 = bet3 * ra_a2 / ra_gm1
    b7 = ra_mx*bx
    ! NOTE: For a calorically perfect gas, gaml=gamr=ra_gam=gam, and
    !       I did the math and b8 should equal zero (at least analytically)
    b8 = (ra_ec*(eir-eil) + (ra_pres/ra_ec*drho-dpres)/ra_gm1) * abs(ra_vn)
    !
    ! Compute the dissipative flux = 1/2*sum(ALFi*RALMi*right eigenvector)
    !
    fec_d = b1
    fmx_d = b1*ra_mx + b2*unx + bx
    fee_d = b1*ra_h0 + b2*ra_vn - b6 + b7 + b8
    !
    ! Compute the left(-) convective fluxes
    !
    fec_l = vnl*rl
    fmx_l = vxl*fec_l + presl*unx
    fee_l = vnl*etotl + presl*vnl
    !
    ! Compute the right(+) convective fluxes
    !
    fec_r = vnr*rr
    fmx_r = vxr*fec_r + presr*unx
    fee_r = vnr*etotr + presr*vnr
    !
    ! Compute the total upwind flux from the Roe scheme
    !
    return_value(n,nec) = half * ds * (fec_l + fec_r - fec_d)
    return_value(n,nmx) = half * ds * (fmx_l + fmx_r - fmx_d)
    return_value(n,nee) = half * ds * (fee_l + fee_r - fee_d)
    !
    ! Save the left and right flux for species continuity equations
    !
    lflux(n) = half*ds*(fec_l + ra_ec*abs(ra_vn) - b1*ra_weight_l)
    rflux(n) = half*ds*(fec_r - ra_ec*abs(ra_vn) - b1*ra_weight_r)
    !
  end do
  !
end subroutine roe_upwind_flux_tpg_1d
!
!###############################################################################
!
pure subroutine hllc_upwind_flux_gpu_3d(n1,n2,facenrm,face_cv,face_gp, &
                                        lflux,rflux,return_value)
  !
  ! Routine to calculate the convective flux in any given direction
  ! for a mixture of gases using the HLLC method
  !
  ! Based on the journal article:
  ! Batten, P., Clarke, N., Lambert, C. and Causon, D.M.,
  ! ON THE CHOICE OF WAVE SPEEDS FOR THE HLLC RIEMANN SOLVER
  ! SIAM, Journal for Scientific Computation,
  ! Vol. 18, No. 6, pp. 1553-1570, Nov. 1997.
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,rnz,unx,uny,unz
  real(wp) :: rl,vxl,vyl,vzl,presl,gaml,q2l,etotl,vnl,al,sl,rl_inv
  real(wp) :: rr,vxr,vyr,vzr,presr,gamr,q2r,etotr,vnr,ar,sr,rr_inv
  real(wp) :: fec_l,fmx_l,fmy_l,fmz_l,fee_l,ra_weight_l
  real(wp) :: fec_r,fmx_r,fmy_r,fmz_r,fee_r,ra_weight_r
  real(wp) :: ra_ec,ra_mx,ra_my,ra_mz,ra_pres
  real(wp) :: ra_a,ra_vn,ra_gam,ra_gm1
  real(wp) :: dvx,dvy,dvz,dq2,pstar,sm
  real(wp) :: slmsmi,dsl,pstrmpl,rhosl,rhousl,rhovsl,rhowsl,rhoesl,efix_l
  real(wp) :: srmsmi,dsr,pstrmpr,rhosr,rhousr,rhovsr,rhowsr,rhoesr,efix_r
  real(wp) :: dels
  real(wp) :: sdrostrl,sdrostrr
  !
continue
  !
  ! Compute the continuity, momentum and energy fluxes
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp,facenrm) &
  !$ACC&              PRESENT(lflux,rflux,return_value)
  do n = n1,n2
    !
    ! Cell face area
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    rl     = face_cv(n,nec,1)
    rl_inv = one/rl
    vxl    = face_cv(n,nmx,1)*rl_inv
    vyl    = face_cv(n,nmy,1)*rl_inv
    vzl    = face_cv(n,nmz,1)*rl_inv
    etotl  = face_cv(n,nee,1)
    presl  = face_gp(n,npres,1)
    gaml   = face_gp(n,ngamm,1)
    !
    rr     = face_cv(n,nec,2)
    rr_inv = one/rl
    vxr    = face_cv(n,nmx,2)*rr_inv
    vyr    = face_cv(n,nmy,2)*rr_inv
    vzr    = face_cv(n,nmz,2)*rr_inv
    etotr  = face_cv(n,nee,2)
    presr  = face_gp(n,npres,2)
    gamr   = face_gp(n,ngamm,2)
    !
    ! Left(-) state variables
    !
    q2l = vxl*vxl + vyl*vyl + vzl*vzl
    vnl = vxl*unx + vyl*uny + vzl*unz
    al  = sqrt(max(gaml*presl/rl,rndoff))
    !
    ! Right(+) state variables
    !
    q2r = vxr*vxr + vyr*vyr + vzr*vzr
    vnr = vxr*unx + vyr*uny + vzr*unz
    ar  = sqrt(max(gamr*presr/rr,rndoff))
    !
    ! Left(-) fluxes
    !
    fec_l = rl*vnl
    fmx_l = vxl*fec_l + presl*unx
    fmy_l = vyl*fec_l + presl*uny
    fmz_l = vzl*fec_l + presl*unz
    fee_l = (etotl + presl)*vnl
    !
    ! Right(+) fluxes
    !
    fec_r = rr*vnr
    fmx_r = vxr*fec_r + presr*unx
    fmy_r = vyr*fec_r + presr*uny
    fmz_r = vzr*fec_r + presr*unz
    fee_r = (etotr + presr)*vnr
    !
    ! Roe averaged density
    !
    ra_ec = sqrt(rl*rr)
    !
    ! Roe averaging weights
    !
    ra_weight_l = rl / (rl+ra_ec)
    ra_weight_r = one - ra_weight_l
    !
    ! Roe averaged velocities
    !
    ra_mx = ra_weight_r*vxr + vxl*ra_weight_l
    ra_my = ra_weight_r*vyr + vyl*ra_weight_l
    ra_mz = ra_weight_r*vzr + vzl*ra_weight_l
    !
    ! Roe averaged contravarient velocity
    !
    ra_vn = ra_mx*unx + ra_my*uny + ra_mz*unz
    !
    ! Roe averaged pressure (special averaging for pressure used here)
    !
    ra_pres = ra_weight_r*presl + presr*ra_weight_l
    !
    ! Roe averaged gamma
    !
    ra_gam = ra_weight_r*gamr + gaml*ra_weight_l
    ra_gm1 = ra_gam - one
    !
    ! Roe averaged speed of sound
    !
    !ra_a = sqrt(ra_gam*ra_pres/ra_ec)
    !
    dvx = vxr - vxl
    dvy = vyr - vyl
    dvz = vzr - vzl
    dq2 = dvx*dvx + dvy*dvy + dvz*dvz
    !
    ra_a = ra_gam*ra_pres/ra_ec + &
           half*ra_ec*ra_gm1*dq2 / ((sqrt(rl)+sqrt(rr))**2)
    ra_a = sqrt(ra_a)
    !
    ! Acoustic wave speed estimates
    !
    sl = min(vnl-al, ra_vn-ra_a)
    sr = max(vnr+ar, ra_vn+ra_a)
    !
    ! Compute local smallest allowable acoustic eigenvalue
    !
    dels = max(abs(upaefix)*(abs(ra_mx)+abs(ra_my)+abs(ra_mz)+ra_a), zero)
    !
    ! Entropy fix for acoustic eigenvalues
    !
    efix_l = dels + max(abs(sl)-dels, &
                        min(zero,half*(sl*sl-dels*dels)/(max(dels,rndoff))))
    efix_r = dels + max(abs(sr)-dels, &
                        min(zero,half*(sr*sr-dels*dels)/(max(dels,rndoff))))
    !
    sl = sign(efix_l, sl)
    sr = sign(efix_r, sr)
    !
    lflux(n) = zero
    rflux(n) = zero
    !
    ! Assign/compute fluxes depending on wave speeds
    !
    if (sl > zero) then
      !
      ! Left state fluxes
      !
      lflux(n) = fec_l
      !
      return_value(n,nec) = fec_l
      return_value(n,nmx) = fmx_l
      return_value(n,nmy) = fmy_l
      return_value(n,nmz) = fmz_l
      return_value(n,nee) = fee_l
      !
    else if (sr < zero) then
      !
      ! Right state fluxes
      !
      rflux(n) = fec_r
      !
      return_value(n,nec) = fec_r
      return_value(n,nmx) = fmx_r
      return_value(n,nmy) = fmy_r
      return_value(n,nmz) = fmz_r
      return_value(n,nee) = fee_r
      !
    else
      !
      ! HLLC interior state based fluxes
      !
      ! Compute the velocity of the contact discontinuity
      !
      dsl = sl - vnl
      dsr = sr - vnr
      !
      sm = ((rr*vnr*dsr-presr) - (rl*vnl*dsl-presl)) / (rr*dsr-rl*dsl)
      !
      ! Compute the interior pressure
      !
      pstar = rl*(vnl-sl)*(vnl-sm) + presl
      !
      if (sm > zero) then
        !
        ! Compute the HLLC fluxes using left interior states
        !
        slmsmi  = one/(sl-sm)
        pstrmpl = pstar-presl
        !
        rhosl  = dsl*rl*slmsmi
        rhousl = rhosl*vxl + pstrmpl*unx*slmsmi
        rhovsl = rhosl*vyl + pstrmpl*uny*slmsmi
        rhowsl = rhosl*vzl + pstrmpl*unz*slmsmi
        rhoesl = (dsl*etotl - presl*vnl + pstar*sm)*slmsmi
        !
        sdrostrl = sl*(rhosl-rl)
        !
        ! Mean flow equations
        !
        lflux(n) = fec_l + sdrostrl
        !
        return_value(n,nec) = fec_l + sdrostrl
        return_value(n,nmx) = fmx_l + sl*(rhousl-rl*vxl)
        return_value(n,nmy) = fmy_l + sl*(rhovsl-rl*vyl)
        return_value(n,nmz) = fmz_l + sl*(rhowsl-rl*vzl)
        return_value(n,nee) = fee_l + sl*(rhoesl-etotl)
        !
      else
        !
        ! Compute the HLLC fluxes using right interior states
        !
        srmsmi  = one/(sr-sm)
        pstrmpr = pstar-presr
        !
        rhosr  = rr*dsr*srmsmi
        rhousr = rhosr*vxr + pstrmpr*unx*srmsmi
        rhovsr = rhosr*vyr + pstrmpr*uny*srmsmi
        rhowsr = rhosr*vzr + pstrmpr*unz*srmsmi
        rhoesr = (dsr*etotr - presr*vnr + pstar*sm)*srmsmi
        !
        sdrostrr = sr*(rhosr-rr)
        !
        ! Mean flow equations
        !
        rflux(n) = fec_r + sdrostrr
        !
        return_value(n,nec) = fec_r + sdrostrr
        return_value(n,nmx) = fmx_r + sr*(rhousr-rr*vxr)
        return_value(n,nmy) = fmy_r + sr*(rhovsr-rr*vyr)
        return_value(n,nmz) = fmz_r + sr*(rhowsr-rr*vzr)
        return_value(n,nee) = fee_r + sr*(rhoesr-etotr)
        !
      end if
      !
    end if
    !
    ! Multiply by cell face area
    !
    lflux(n) = lflux(n) * ds
    rflux(n) = rflux(n) * ds
    !
    return_value(n,nec) = return_value(n,nec) * ds
    return_value(n,nmx) = return_value(n,nmx) * ds
    return_value(n,nmy) = return_value(n,nmy) * ds
    return_value(n,nmz) = return_value(n,nmz) * ds
    return_value(n,nee) = return_value(n,nee) * ds
    !
  end do
  !
  ! Finished computing the convective flux in any given direction
  ! for a mixture of gases using the HLLC method
  !
end subroutine hllc_upwind_flux_gpu_3d
!
!###############################################################################
!
pure subroutine hllc_upwind_flux_gpu_2d(n1,n2,facenrm,face_cv,face_gp, &
                                        lflux,rflux,return_value)
  !
  ! Routine to calculate the convective flux in any given direction
  ! for a mixture of gases using the HLLC method
  !
  ! Based on the journal article:
  ! Batten, P., Clarke, N., Lambert, C. and Causon, D.M.,
  ! ON THE CHOICE OF WAVE SPEEDS FOR THE HLLC RIEMANN SOLVER
  ! SIAM, Journal for Scientific Computation,
  ! Vol. 18, No. 6, pp. 1553-1570, Nov. 1997.
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,rny,unx,uny
  real(wp) :: rl,vxl,vyl,presl,gaml,q2l,etotl,vnl,al,sl,rl_inv
  real(wp) :: rr,vxr,vyr,presr,gamr,q2r,etotr,vnr,ar,sr,rr_inv
  real(wp) :: fec_l,fmx_l,fmy_l,fee_l,ra_weight_l
  real(wp) :: fec_r,fmx_r,fmy_r,fee_r,ra_weight_r
  real(wp) :: ra_ec,ra_mx,ra_my,ra_pres
  real(wp) :: ra_a,ra_vn,ra_gam,ra_gm1
  real(wp) :: dvx,dvy,dq2,pstar,sm
  real(wp) :: slmsmi,dsl,pstrmpl,rhosl,rhousl,rhovsl,rhoesl,efix_l
  real(wp) :: srmsmi,dsr,pstrmpr,rhosr,rhousr,rhovsr,rhoesr,efix_r
  real(wp) :: dels
  real(wp) :: sdrostrl,sdrostrr
  !
continue
  !
  ! Compute the continuity, momentum and energy fluxes
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp,facenrm) &
  !$ACC&              PRESENT(lflux,rflux,return_value)
  do n = n1,n2
    !
    ! Cell face area
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    rl     = face_cv(n,nec,1)
    rl_inv = one/rl
    vxl    = face_cv(n,nmx,1)*rl_inv
    vyl    = face_cv(n,nmy,1)*rl_inv
    etotl  = face_cv(n,nee,1)
    presl  = face_gp(n,npres,1)
    gaml   = face_gp(n,ngamm,1)
    !
    rr     = face_cv(n,nec,2)
    rr_inv = one/rl
    vxr    = face_cv(n,nmx,2)*rr_inv
    vyr    = face_cv(n,nmy,2)*rr_inv
    etotr  = face_cv(n,nee,2)
    presr  = face_gp(n,npres,2)
    gamr   = face_gp(n,ngamm,2)
    !
    ! Left(-) state variables
    !
    q2l = vxl*vxl + vyl*vyl
    vnl = vxl*unx + vyl*uny
    al  = sqrt(max(gaml*presl/rl,rndoff))
    !
    ! Right(+) state variables
    !
    q2r = vxr*vxr + vyr*vyr
    vnr = vxr*unx + vyr*uny
    ar  = sqrt(max(gamr*presr/rr,rndoff))
    !
    ! Left(-) fluxes
    !
    fec_l = rl*vnl
    fmx_l = vxl*fec_l + presl*unx
    fmy_l = vyl*fec_l + presl*uny
    fee_l = (etotl + presl)*vnl
    !
    ! Right(+) fluxes
    !
    fec_r = rr*vnr
    fmx_r = vxr*fec_r + presr*unx
    fmy_r = vyr*fec_r + presr*uny
    fee_r = (etotr + presr)*vnr
    !
    ! Roe averaged density
    !
    ra_ec = sqrt(rl*rr)
    !
    ! Roe averaging weights
    !
    ra_weight_l = rl / (rl+ra_ec)
    ra_weight_r = one - ra_weight_l
    !
    ! Roe averaged velocities
    !
    ra_mx = ra_weight_r*vxr + vxl*ra_weight_l
    ra_my = ra_weight_r*vyr + vyl*ra_weight_l
    !
    ! Roe averaged contravarient velocity
    !
    ra_vn = ra_mx*unx + ra_my*uny
    !
    ! Roe averaged pressure (special averaging for pressure used here)
    !
    ra_pres = ra_weight_r*presl + presr*ra_weight_l
    !
    ! Roe averaged gamma
    !
    ra_gam = ra_weight_r*gamr + gaml*ra_weight_l
    ra_gm1 = ra_gam - one
    !
    ! Roe averaged speed of sound
    !
    !ra_a = sqrt(ra_gam*ra_pres/ra_ec)
    !
    dvx = vxr - vxl
    dvy = vyr - vyl
    dq2 = dvx*dvx + dvy*dvy
    !
    ra_a = ra_gam*ra_pres/ra_ec + &
           half*ra_ec*ra_gm1*dq2 / ((sqrt(rl)+sqrt(rr))**2)
    ra_a = sqrt(ra_a)
    !
    ! Acoustic wave speed estimates
    !
    sl = min(vnl-al, ra_vn-ra_a)
    sr = max(vnr+ar, ra_vn+ra_a)
    !
    ! Compute local smallest allowable acoustic eigenvalue
    !
    dels = max(abs(upaefix)*(abs(ra_mx)+abs(ra_my)+ra_a), zero)
    !
    ! Entropy fix for acoustic eigenvalues
    !
    efix_l = dels + max(abs(sl)-dels, &
                        min(zero,half*(sl*sl-dels*dels)/(max(dels,rndoff))))
    efix_r = dels + max(abs(sr)-dels, &
                        min(zero,half*(sr*sr-dels*dels)/(max(dels,rndoff))))
    !
    sl = sign(efix_l, sl)
    sr = sign(efix_r, sr)
    !
    lflux(n) = zero
    rflux(n) = zero
    !
    ! Assign/compute fluxes depending on wave speeds
    !
    if (sl > zero) then
      !
      ! Left state fluxes
      !
      lflux(n) = fec_l
      !
      return_value(n,nec) = fec_l
      return_value(n,nmx) = fmx_l
      return_value(n,nmy) = fmy_l
      return_value(n,nee) = fee_l
      !
    else if (sr < zero) then
      !
      ! Right state fluxes
      !
      rflux(n) = fec_r
      !
      return_value(n,nec) = fec_r
      return_value(n,nmx) = fmx_r
      return_value(n,nmy) = fmy_r
      return_value(n,nee) = fee_r
      !
    else
      !
      ! HLLC interior state based fluxes
      !
      ! Compute the velocity of the contact discontinuity
      !
      dsl = sl - vnl
      dsr = sr - vnr
      !
      sm = ((rr*vnr*dsr-presr) - (rl*vnl*dsl-presl)) / (rr*dsr-rl*dsl)
      !
      ! Compute the interior pressure
      !
      pstar = rl*(vnl-sl)*(vnl-sm) + presl
      !
      if (sm > zero) then
        !
        ! Compute the HLLC fluxes using left interior states
        !
        slmsmi  = one/(sl-sm)
        pstrmpl = pstar-presl
        !
        rhosl  = dsl*rl*slmsmi
        rhousl = rhosl*vxl + pstrmpl*unx*slmsmi
        rhovsl = rhosl*vyl + pstrmpl*uny*slmsmi
        rhoesl = (dsl*etotl - presl*vnl + pstar*sm)*slmsmi
        !
        sdrostrl = sl*(rhosl-rl)
        !
        ! Mean flow equations
        !
        lflux(n) = fec_l + sdrostrl
        !
        return_value(n,nec) = fec_l + sdrostrl
        return_value(n,nmx) = fmx_l + sl*(rhousl-rl*vxl)
        return_value(n,nmy) = fmy_l + sl*(rhovsl-rl*vyl)
        return_value(n,nee) = fee_l + sl*(rhoesl-etotl)
        !
      else
        !
        ! Compute the HLLC fluxes using right interior states
        !
        srmsmi  = one/(sr-sm)
        pstrmpr = pstar-presr
        !
        rhosr  = rr*dsr*srmsmi
        rhousr = rhosr*vxr + pstrmpr*unx*srmsmi
        rhovsr = rhosr*vyr + pstrmpr*uny*srmsmi
        rhoesr = (dsr*etotr - presr*vnr + pstar*sm)*srmsmi
        !
        sdrostrr = sr*(rhosr-rr)
        !
        ! Mean flow equations
        !
        rflux(n) = fec_r + sdrostrr
        !
        return_value(n,nec) = fec_r + sdrostrr
        return_value(n,nmx) = fmx_r + sr*(rhousr-rr*vxr)
        return_value(n,nmy) = fmy_r + sr*(rhovsr-rr*vyr)
        return_value(n,nee) = fee_r + sr*(rhoesr-etotr)
        !
      end if
      !
    end if
    !
    ! Multiply by cell face area
    !
    lflux(n) = lflux(n) * ds
    rflux(n) = rflux(n) * ds
    !
    return_value(n,nec) = return_value(n,nec) * ds
    return_value(n,nmx) = return_value(n,nmx) * ds
    return_value(n,nmy) = return_value(n,nmy) * ds
    return_value(n,nee) = return_value(n,nee) * ds
    !
  end do
  !
  ! Finished computing the convective flux in any given direction
  ! for a mixture of gases using the HLLC method
  !
end subroutine hllc_upwind_flux_gpu_2d
!
!###############################################################################
!
pure subroutine hllc_upwind_flux_gpu_1d(n1,n2,facenrm,face_cv,face_gp, &
                                        lflux,rflux,return_value)
  !
  ! Routine to calculate the convective flux in any given direction
  ! for a mixture of gases using the HLLC method
  !
  ! Based on the journal article:
  ! Batten, P., Clarke, N., Lambert, C. and Causon, D.M.,
  ! ON THE CHOICE OF WAVE SPEEDS FOR THE HLLC RIEMANN SOLVER
  ! SIAM, Journal for Scientific Computation,
  ! Vol. 18, No. 6, pp. 1553-1570, Nov. 1997.
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: ds,rnx,unx
  real(wp) :: rl,vxl,presl,gaml,q2l,etotl,vnl,al,sl,rl_inv
  real(wp) :: rr,vxr,presr,gamr,q2r,etotr,vnr,ar,sr,rr_inv
  real(wp) :: fec_l,fmx_l,fee_l,ra_weight_l
  real(wp) :: fec_r,fmx_r,fee_r,ra_weight_r
  real(wp) :: ra_ec,ra_mx,ra_pres
  real(wp) :: ra_a,ra_vn,ra_gam,ra_gm1
  real(wp) :: dvx,pstar,sm
  real(wp) :: slmsmi,dsl,pstrmpl,rhosl,rhousl,rhoesl,efix_l
  real(wp) :: srmsmi,dsr,pstrmpr,rhosr,rhousr,rhoesr,efix_r
  real(wp) :: dels
  real(wp) :: sdrostrl,sdrostrr
  !
continue
  !
  ! Compute the continuity, momentum and energy fluxes
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,face_gp,facenrm) &
  !$ACC&              PRESENT(lflux,rflux,return_value)
  do n = n1,n2
    !
    ! Cell face area
    !
    rnx = facenrm(1,n)
    ds = sqrt( rnx*rnx )
    unx = rnx/ds
    !
    rl     = face_cv(n,nec,1)
    rl_inv = one/rl
    vxl    = face_cv(n,nmx,1)*rl_inv
    etotl  = face_cv(n,nee,1)
    presl  = face_gp(n,npres,1)
    gaml   = face_gp(n,ngamm,1)
    !
    rr     = face_cv(n,nec,2)
    rr_inv = one/rl
    vxr    = face_cv(n,nmx,2)*rr_inv
    etotr  = face_cv(n,nee,2)
    presr  = face_gp(n,npres,2)
    gamr   = face_gp(n,ngamm,2)
    !
    ! Left(-) state variables
    !
    q2l = vxl*vxl
    vnl = vxl*unx
    al  = sqrt(max(gaml*presl/rl,rndoff))
    !
    ! Right(+) state variables
    !
    q2r = vxr*vxr
    vnr = vxr*unx
    ar  = sqrt(max(gamr*presr/rr,rndoff))
    !
    ! Left(-) fluxes
    !
    fec_l = rl*vnl
    fmx_l = vxl*fec_l + presl*unx
    fee_l = (etotl + presl)*vnl
    !
    ! Right(+) fluxes
    !
    fec_r = rr*vnr
    fmx_r = vxr*fec_r + presr*unx
    fee_r = (etotr + presr)*vnr
    !
    ! Roe averaged density
    !
    ra_ec = sqrt(rl*rr)
    !
    ! Roe averaging weights
    !
    ra_weight_l = rl / (rl+ra_ec)
    ra_weight_r = one - ra_weight_l
    !
    ! Roe averaged velocities
    !
    ra_mx = ra_weight_r*vxr + vxl*ra_weight_l
    !
    ! Roe averaged contravarient velocity
    !
    ra_vn = ra_mx*unx
    !
    ! Roe averaged pressure (special averaging for pressure used here)
    !
    ra_pres = ra_weight_r*presl + presr*ra_weight_l
    !
    ! Roe averaged gamma
    !
    ra_gam = ra_weight_r*gamr + gaml*ra_weight_l
    ra_gm1 = ra_gam - one
    !
    ! Roe averaged speed of sound
    !
    !ra_a = sqrt(ra_gam*ra_pres/ra_ec)
    !
    dvx = vxr - vxl
    !
    ra_a = ra_gam*ra_pres/ra_ec + &
           half*ra_ec*ra_gm1*dvx*dvx / ((sqrt(rl)+sqrt(rr))**2)
    ra_a = sqrt(ra_a)
    !
    ! Acoustic wave speed estimates
    !
    sl = min(vnl-al, ra_vn-ra_a)
    sr = max(vnr+ar, ra_vn+ra_a)
    !
    ! Compute local smallest allowable acoustic eigenvalue
    !
    dels = max(abs(upaefix)*(abs(ra_mx)+ra_a), zero)
    !
    ! Entropy fix for acoustic eigenvalues
    !
    efix_l = dels + max(abs(sl)-dels, &
                        min(zero,half*(sl*sl-dels*dels)/(max(dels,rndoff))))
    efix_r = dels + max(abs(sr)-dels, &
                        min(zero,half*(sr*sr-dels*dels)/(max(dels,rndoff))))
    !
    sl = sign(efix_l, sl)
    sr = sign(efix_r, sr)
    !
    lflux(n) = zero
    rflux(n) = zero
    !
    ! Assign/compute fluxes depending on wave speeds
    !
    if (sl > zero) then
      !
      ! Left state fluxes
      !
      lflux(n) = fec_l
      !
      return_value(n,nec) = fec_l
      return_value(n,nmx) = fmx_l
      return_value(n,nee) = fee_l
      !
    else if (sr < zero) then
      !
      ! Right state fluxes
      !
      rflux(n) = fec_r
      !
      return_value(n,nec) = fec_r
      return_value(n,nmx) = fmx_r
      return_value(n,nee) = fee_r
      !
    else
      !
      ! HLLC interior state based fluxes
      !
      ! Compute the velocity of the contact discontinuity
      !
      dsl = sl - vnl
      dsr = sr - vnr
      !
      sm = ((rr*vnr*dsr-presr) - (rl*vnl*dsl-presl)) / (rr*dsr-rl*dsl)
      !
      ! Compute the interior pressure
      !
      pstar = rl*(vnl-sl)*(vnl-sm) + presl
      !
      if (sm > zero) then
        !
        ! Compute the HLLC fluxes using left interior states
        !
        slmsmi  = one/(sl-sm)
        pstrmpl = pstar-presl
        !
        rhosl  = dsl*rl*slmsmi
        rhousl = rhosl*vxl + pstrmpl*unx*slmsmi
        rhoesl = (dsl*etotl - presl*vnl + pstar*sm)*slmsmi
        !
        sdrostrl = sl*(rhosl-rl)
        !
        ! Mean flow equations
        !
        lflux(n) = fec_l + sdrostrl
        !
        return_value(n,nec) = fec_l + sdrostrl
        return_value(n,nmx) = fmx_l + sl*(rhousl-rl*vxl)
        return_value(n,nee) = fee_l + sl*(rhoesl-etotl)
        !
      else
        !
        ! Compute the HLLC fluxes using right interior states
        !
        srmsmi  = one/(sr-sm)
        pstrmpr = pstar-presr
        !
        rhosr  = rr*dsr*srmsmi
        rhousr = rhosr*vxr + pstrmpr*unx*srmsmi
        rhoesr = (dsr*etotr - presr*vnr + pstar*sm)*srmsmi
        !
        sdrostrr = sr*(rhosr-rr)
        !
        ! Mean flow equations
        !
        rflux(n) = fec_r + sdrostrr
        !
        return_value(n,nec) = fec_r + sdrostrr
        return_value(n,nmx) = fmx_r + sr*(rhousr-rr*vxr)
        return_value(n,nee) = fee_r + sr*(rhoesr-etotr)
        !
      end if
      !
    end if
    !
    ! Multiply by cell face area
    !
    lflux(n) = lflux(n) * ds
    rflux(n) = rflux(n) * ds
    !
    return_value(n,nec) = return_value(n,nec) * ds
    return_value(n,nmx) = return_value(n,nmx) * ds
    return_value(n,nee) = return_value(n,nee) * ds
    !
  end do
  !
  ! Finished computing the convective flux in any given direction
  ! for a mixture of gases using the HLLC method
  !
end subroutine hllc_upwind_flux_gpu_1d
!
!###############################################################################
!
pure subroutine ldfss_upwind_flux_gpu_3d(n1,n2,facenrm,face_cv,face_gp, &
                                         lflux,rflux,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,ds,unx,uny,unz,rnx,rny,rnz
  real(wp) :: rl,rl_inv,vxl,vyl,vzl,vnl,presl,h0l,al,q2l
  real(wp) :: rr,rr_inv,vxr,vyr,vzr,vnr,presr,h0r,ar,q2r
  real(wp) :: alpha_l,beta_l,mach_l,mstar_l
  real(wp) :: alpha_r,beta_r,mach_r,mstar_r
  real(wp) :: machp_half,uc_l,flux_l,pplus_l
  real(wp) :: machm_half,uc_r,flux_r,pplus_r
  real(wp) :: m_half,mstar_half,ahalf
  real(wp) :: del_p,sum_p,p_net
  !
continue
  !
  !$ACC PARALLEL LOOP GANG WORKER VECTOR &
  !$ACC&              PRESENT(facenrm,face_cv,face_gp) &
  !$ACC&              PRESENT(lflux,rflux,return_value)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    rl = face_cv(n,nec,1)
    rl_inv = one/rl
    vxl = face_cv(n,nmx,1)*rl_inv
    vyl = face_cv(n,nmy,1)*rl_inv
    vzl = face_cv(n,nmz,1)*rl_inv
    presl = face_gp(n,npres,1)
    h0l = face_gp(n,nenth,1) + half*(vxl*vxl+vyl*vyl+vzl*vzl)
    vnl = vxl*unx + vyl*uny + vzl*unz
    !
    rr = face_cv(n,nec,2)
    rr_inv = one/rr
    vxr = face_cv(n,nmx,2)*rr_inv
    vyr = face_cv(n,nmy,2)*rr_inv
    vzr = face_cv(n,nmz,2)*rr_inv
    presr = face_gp(n,npres,2)
    h0r = face_gp(n,nenth,2) + half*(vxr*vxr+vyr*vyr+vzr*vzr)
    vnr = vxr*unx + vyr*uny + vzr*unz
    !
    ahalf = lflux(n)
    !
    mach_l = vnl/ahalf
    mach_r = vnr/ahalf
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
    return_value(n,nec) = flux_l     + flux_r
    return_value(n,nmx) = flux_l*vxl + flux_r*vxr + p_net*unx
    return_value(n,nmy) = flux_l*vyl + flux_r*vyr + p_net*uny
    return_value(n,nmz) = flux_l*vzl + flux_r*vzr + p_net*unz
    return_value(n,nee) = flux_l*h0l + flux_r*h0r
    !
    lflux(n) = flux_l
    rflux(n) = flux_r
    !
  end do
  !
end subroutine ldfss_upwind_flux_gpu_3d
!
!###############################################################################
!
pure subroutine ldfss_upwind_flux_gpu_2d(n1,n2,facenrm,face_cv,face_gp, &
                                         lflux,rflux,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,ds,unx,uny,rnx,rny
  real(wp) :: rl,rl_inv,vxl,vyl,vnl,presl,h0l,al,q2l
  real(wp) :: rr,rr_inv,vxr,vyr,vnr,presr,h0r,ar,q2r
  real(wp) :: alpha_l,beta_l,mach_l,mstar_l
  real(wp) :: alpha_r,beta_r,mach_r,mstar_r
  real(wp) :: machp_half,uc_l,flux_l,pplus_l
  real(wp) :: machm_half,uc_r,flux_r,pplus_r
  real(wp) :: m_half,mstar_half,ahalf
  real(wp) :: del_p,sum_p,p_net
  !
continue
  !
  !$ACC PARALLEL LOOP GANG WORKER VECTOR &
  !$ACC&              PRESENT(facenrm,face_cv,face_gp) &
  !$ACC&              PRESENT(lflux,rflux,return_value)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    rl = face_cv(n,nec,1)
    rl_inv = one/rl
    vxl = face_cv(n,nmx,1)*rl_inv
    vyl = face_cv(n,nmy,1)*rl_inv
    presl = face_gp(n,npres,1)
    h0l = face_gp(n,nenth,1) + half*(vxl*vxl+vyl*vyl)
    vnl = vxl*unx + vyl*uny
    !
    rr = face_cv(n,nec,2)
    rr_inv = one/rr
    vxr = face_cv(n,nmx,2)*rr_inv
    vyr = face_cv(n,nmy,2)*rr_inv
    presr = face_gp(n,npres,2)
    h0r = face_gp(n,nenth,2) + half*(vxr*vxr+vyr*vyr)
    vnr = vxr*unx + vyr*uny
    !
    ahalf = lflux(n)
    !
    mach_l = vnl/ahalf
    mach_r = vnr/ahalf
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
    return_value(n,nec) = flux_l     + flux_r
    return_value(n,nmx) = flux_l*vxl + flux_r*vxr + p_net*unx
    return_value(n,nmy) = flux_l*vyl + flux_r*vyr + p_net*uny
    return_value(n,nee) = flux_l*h0l + flux_r*h0r
    !
    lflux(n) = flux_l
    rflux(n) = flux_r
    !
  end do
  !
end subroutine ldfss_upwind_flux_gpu_2d
!
!###############################################################################
!
pure subroutine ldfss_upwind_flux_gpu_1d(n1,n2,facenrm,face_cv,face_gp, &
                                         lflux,rflux,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,ds,unx,rnx
  real(wp) :: rl,vxl,presl,h0l,al,q2l
  real(wp) :: rr,vxr,presr,h0r,ar,q2r
  real(wp) :: alpha_l,beta_l,mach_l,mstar_l
  real(wp) :: alpha_r,beta_r,mach_r,mstar_r
  real(wp) :: machp_half,uc_l,flux_l,pplus_l
  real(wp) :: machm_half,uc_r,flux_r,pplus_r
  real(wp) :: m_half,mstar_half,ahalf
  real(wp) :: del_p,sum_p,p_net
  !
continue
  !
  !$ACC PARALLEL LOOP GANG WORKER VECTOR &
  !$ACC&              PRESENT(facenrm,face_cv,face_gp) &
  !$ACC&              PRESENT(lflux,rflux,return_value)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    ds = sqrt( rnx*rnx )
    unx = rnx/ds
    !
    rl = face_cv(n,nec,1)
    vxl = face_cv(n,nmx,1)/rl
    presl = face_gp(n,npres,1)
    h0l = face_gp(n,nenth,1) + half*vxl*vxl
    !
    rr = face_cv(n,nec,2)
    vxr = face_cv(n,nmx,2)/rr
    presr = face_gp(n,npres,2)
    h0r = face_gp(n,nenth,2) + half*vxr*vxr
    !
    ahalf = lflux(n)
    !
    mach_l = vxl*unx/ahalf
    mach_r = vxr*unx/ahalf
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
    return_value(n,nec) = flux_l     + flux_r
    return_value(n,nmx) = flux_l*vxl + flux_r*vxr + p_net*unx
    return_value(n,nee) = flux_l*h0l + flux_r*h0r
    !
    lflux(n) = flux_l
    rflux(n) = flux_r
    !
  end do
  !
end subroutine ldfss_upwind_flux_gpu_1d
!
!###############################################################################
!
pure subroutine ausm_plus_upwind_flux_gpu_3d(n1,n2,facenrm,face_cv,face_gp, &
                                             lflux,rflux,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: gm1,ds,unx,uny,unz,rnx,rny,rnz
  real(wp) :: rl,rl_inv,vxl,vyl,vzl,presl,al,h0l,vnl,q2l
  real(wp) :: rr,rr_inv,vxr,vyr,vzr,presr,ar,h0r,vnr,q2r
  real(wp) :: alr,rm2,rmref2,coef1,coef2,prec_f
  real(wp) :: aml,absml,am1p,aml2,ampl,ppl,f1l
  real(wp) :: amr,absmr,am1m,amr2,ammr,pmr,f1r
  real(wp) :: amhalf,denom,amhalfp,amhalfm,phalf
  !
  !.. Local Parameters ..
  real(wp), parameter :: a1 = one/eight
  real(wp), parameter :: a2 = half
  real(wp), parameter :: a3 = three/eight
  !
  real(wp), parameter :: b1 = three/sixteen
  real(wp), parameter :: b2 = five/eight
  real(wp), parameter :: b3 = three*five/sixteen
  real(wp), parameter :: b4 = half
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  gm1 = gam - one
  !
  !$ACC PARALLEL LOOP GANG WORKER VECTOR &
  !$ACC&              PRESENT(facenrm,face_cv,face_gp,aave) &
  !$ACC&              PRESENT(htotl,htotr,lflux,rflux,return_value)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    rl = face_cv(n,nec,1)
    rl_inv = one/rl
    vxl = face_cv(n,nmx,1)*rl_inv
    vyl = face_cv(n,nmy,1)*rl_inv
    vzl = face_cv(n,nmz,1)*rl_inv
    q2l = vxl*vxl + vyl*vyl + vzl*vzl
    presl = face_gp(n,npres,1)
    h0l = half*q2l + gam/gm1 * presl/rl
    vnl = vxl*unx + vyl*uny + vzl*unz
    al = sqrt(max(gam*presl/rl,rndoff))
    !
    rr = face_cv(n,nec,2)
    rr_inv = one/rr
    vxr = face_cv(n,nmx,2)*rr_inv
    vyr = face_cv(n,nmy,2)*rr_inv
    vzr = face_cv(n,nmz,2)*rr_inv
    q2r = vxr*vxr + vyr*vyr + vzr*vzr
    presr = face_gp(n,npres,2)
    h0r = half*q2r + gam/gm1 * presr/rr
    vnr = vxr*unx + vyr*uny + vzr*unz
    ar = sqrt(max(gam*presr/rr,rndoff))
    !
    alr    = half*(al + ar)
    rm2    = half*(q2l+q2r)/(alr*alr)
    rmref2 = min( one , max( rm2 , max(eps4,one4*rm2) ) )
    coef1  = half * (one+rmref2)
    coef2  = one - coef1
    prec_f = one - (one-sqrt(rmref2))**2
    alr    = prec_f*alr
    !
    aml   = (coef1*vnl + coef2*vnr)/alr
    amr   = (coef1*vnr + coef2*vnl)/alr
    absml = abs(aml)
    absmr = abs(amr)
    am1p  = half*(aml + absml)
    am1m  = half*(amr - absmr)
    !
    if (absml < one) then
     !aml1 = aml + one
     !addp = (aml1*(aml - one))**2
     !aml2 = one4*aml1**2
     !ampl = aml2 + one8*addp
     !ppl  = aml2*(two - aml) + (three/sixteen)*aml*addp
      aml2 = aml*aml
      ampl =  a1*aml2*aml2 + a2*aml + a3
      ppl  =  b1*aml2*aml2*aml - b2*aml2*aml + b3*aml + b4
    else
      ampl = am1p
      ppl  = am1p/aml
    end if
    !
    if (absmr < one) then
     !amr1 = amr - one
     !addm = (amr1*(amr + one))**2
     !amr2 = one4*amr1**2
     !ammr = -amr2 - one8*addm
     !pmr  = amr2*(two + amr) - (three/sixteen)*amr*addm
      amr2 = amr*amr
      ammr = -a1*amr2*amr2 + a2*amr - a3
      pmr  = -b1*amr2*amr2*amr + b2*amr2*amr - b3*amr + b4
    else
      ammr = am1m
      pmr  = am1m/amr
    end if
    !
    denom   = ( rmref2*(presl+presr)*(rl+rr) ) / sqrt(rl*rr)
    amhalf  = ampl + ammr + half*(presl-presr)*((ampl-ammr) - (am1p-am1m))/denom
    amhalfp = half*(amhalf + abs(amhalf))
    amhalfm = amhalf - amhalfp
    !
    phalf = ppl*presl + pmr*presr - ppl*pmr*(rl+rr)/two*alr*(vnr-vnl)
    !
    f1l = rl*amhalfp*alr
    f1r = rr*amhalfm*alr
    !
    return_value(n,nec) = f1l     + f1r
    return_value(n,nmx) = f1l*vxl + f1r*vxr + phalf*unx
    return_value(n,nmy) = f1l*vyl + f1r*vyr + phalf*uny
    return_value(n,nmz) = f1l*vzl + f1r*vzr + phalf*unz
    return_value(n,nee) = f1l*h0l + f1r*h0r
    !
    lflux(n) = f1l
    rflux(n) = f1r
    !
  end do
  !
end subroutine ausm_plus_upwind_flux_gpu_3d
!
!###############################################################################
!
pure subroutine ausm_plus_upwind_flux_gpu_2d(n1,n2,facenrm,face_cv,face_gp, &
                                             lflux,rflux,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: gm1,ds,unx,uny,rnx,rny
  real(wp) :: rl,rl_inv,vxl,vyl,presl,al,h0l,vnl,q2l
  real(wp) :: rr,rr_inv,vxr,vyr,presr,ar,h0r,vnr,q2r
  real(wp) :: alr,rm2,rmref2,coef1,coef2,prec_f
  real(wp) :: aml,absml,am1p,aml2,ampl,ppl,f1l
  real(wp) :: amr,absmr,am1m,amr2,ammr,pmr,f1r
  real(wp) :: amhalf,denom,amhalfp,amhalfm,phalf
  !
  !.. Local Parameters ..
  real(wp), parameter :: a1 = one/eight
  real(wp), parameter :: a2 = half
  real(wp), parameter :: a3 = three/eight
  !
  real(wp), parameter :: b1 = three/sixteen
  real(wp), parameter :: b2 = five/eight
  real(wp), parameter :: b3 = three*five/sixteen
  real(wp), parameter :: b4 = half
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  gm1 = gam - one
  !
  !$ACC PARALLEL LOOP GANG WORKER VECTOR &
  !$ACC&              PRESENT(facenrm,face_cv,face_gp,aave) &
  !$ACC&              PRESENT(htotl,htotr,lflux,rflux,return_value)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    rl = face_cv(n,nec,1)
    rl_inv = one/rl
    vxl = face_cv(n,nmx,1)*rl_inv
    vyl = face_cv(n,nmy,1)*rl_inv
    q2l = vxl*vxl + vyl*vyl
    presl = face_gp(n,npres,1)
    h0l = half*q2l + gam/gm1 * presl/rl
    vnl = vxl*unx + vyl*uny
    al = sqrt(max(gam*presl/rl,rndoff))
    !
    rr = face_cv(n,nec,2)
    rr_inv = one/rr
    vxr = face_cv(n,nmx,2)*rr_inv
    vyr = face_cv(n,nmy,2)*rr_inv
    coef2  = one - coef1
    prec_f = one - (one-sqrt(rmref2))**2
    alr    = prec_f*alr
    !
    aml   = (coef1*vnl + coef2*vnr)/alr
    amr   = (coef1*vnr + coef2*vnl)/alr
    absml = abs(aml)
    absmr = abs(amr)
    am1p  = half*(aml + absml)
    am1m  = half*(amr - absmr)
    !
    if (absml < one) then
     !aml1 = aml + one
     !addp = (aml1*(aml - one))**2
     !aml2 = one4*aml1**2
     !ampl = aml2 + one8*addp
     !ppl  = aml2*(two - aml) + (three/sixteen)*aml*addp
      aml2 = aml*aml
      ampl =  a1*aml2*aml2 + a2*aml + a3
      ppl  =  b1*aml2*aml2*aml - b2*aml2*aml + b3*aml + b4
    else
      ampl = am1p
      ppl  = am1p/aml
    end if
    !
    if (absmr < one) then
     !amr1 = amr - one
     !addm = (amr1*(amr + one))**2
     !amr2 = one4*amr1**2
     !ammr = -amr2 - one8*addm
     !pmr  = amr2*(two + amr) - (three/sixteen)*amr*addm
      amr2 = amr*amr
      ammr = -a1*amr2*amr2 + a2*amr - a3
      pmr  = -b1*amr2*amr2*amr + b2*amr2*amr - b3*amr + b4
    else
      ammr = am1m
      pmr  = am1m/amr
    end if
    !
    denom   = ( rmref2*(presl+presr)*(rl+rr) ) / sqrt(rl*rr)
    amhalf  = ampl + ammr + half*(presl-presr)*((ampl-ammr) - (am1p-am1m))/denom
    amhalfp = half*(amhalf + abs(amhalf))
    amhalfm = amhalf - amhalfp
    !
    phalf = ppl*presl + pmr*presr - ppl*pmr*(rl+rr)/two*alr*(vnr-vnl)
    !
    f1l = rl*amhalfp*alr
    f1r = rr*amhalfm*alr
    !
    return_value(n,nec) = f1l     + f1r
    return_value(n,nmx) = f1l*vxl + f1r*vxr + phalf*unx
    return_value(n,nmy) = f1l*vyl + f1r*vyr + phalf*uny
    return_value(n,nee) = f1l*h0l + f1r*h0r
    !
    lflux(n) = f1l
    rflux(n) = f1r
    !
  end do
  !
end subroutine ausm_plus_upwind_flux_gpu_2d
!
!###############################################################################
!
pure subroutine ausm_plus_upwind_flux_gpu_nd(n1,n2,facenrm,face_cv,face_gp, &
                                             lflux,rflux,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: gm1,ds
  real(wp) :: rl,presl,hl,al,v2l,a2l,h0l,uul,f1l,q2l
  real(wp) :: rr,presr,hr,ar,v2r,a2r,h0r,uur,f1r,q2r
  real(wp) :: alr,alr2,rm2,rmref2,coef1,coef2,prec_f
  real(wp) :: aml,absml,am1p,aml1,addp,aml2,ampl,ppl
  real(wp) :: amr,absmr,am1m,amr1,addm,amr2,ammr,pmr
  real(wp) :: amhalf,denom,amhalfp,amhalfm,phalf,amco2
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(facenrm,dim=1)) :: unrm,rnrm
  real(wp), dimension(1:size(facenrm,dim=1)) :: vell,velr
  !
  !.. Local Parameters ..
  real(wp), parameter :: amco2m = 1.0e-4_wp
  real(wp), parameter :: alphax = 0.1875_wp
  real(wp), parameter :: betax  = 0.125_wp
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  gm1 = gam - one
  !
  !$ACC PARALLEL LOOP GANG WORKER VECTOR &
  !$ACC&              PRIVATE(unrm,rnrm,vell,velr) &
  !$ACC&              PRESENT(facenrm,face_cv,face_gp,aave) &
  !$ACC&              PRESENT(htotl,htotr,lflux,rflux,return_value)
  do n = n1,n2
    !
    rnrm(1:nr) = facenrm(1:nr,n)
    ds = sqrt( dot_product( rnrm(1:nr) , rnrm(1:nr) ) )
    unrm(1:nr) = rnrm(1:nr) / ds
    !
    rl = face_cv(n,nec,1)
    vell = face_cv(n,nmb:nme,1)/rl
    q2l = dot_product(vell,vell)
    presl = face_gp(n,npres,1)
    !
    rr = face_cv(n,nec,2)
    velr = face_cv(n,nmb:nme,2)/rr
    q2r = dot_product(velr,velr)
    presr = face_gp(n,npres,2)
    !
    hl  = gam*presl/(gm1*rl)
    v2l = half * dot_product(vell,vell)
    h0l = hl + v2l
    a2l = gm1*hl
    al  = sqrt(a2l)
    uul = dot_product(vell,unrm)
    !
    hr  = gam*presr/(gm1*rr)
    v2r = half * dot_product(velr,velr)
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
    return_value(n,nec)     = f1l      + f1r
    return_value(n,nmb:nme) = f1l*vell + f1r*velr + unrm*phalf
    return_value(n,nee)     = f1l*h0l  + f1r*h0r
    !
    lflux(n) = f1l
    rflux(n) = f1r
    !
  end do
  !
end subroutine ausm_plus_upwind_flux_gpu_nd
!
!###############################################################################
!
pure subroutine rotated_rhll_upwind_flux_gpu_3d(n1,n2,facenrm,face_cv,face_gp, &
                                                lflux,rflux,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,ds,unx,uny,unz,rnx,rny,rnz
  real(wp) :: rl_inv,rr_inv,q2l,q2r
  real(wp) :: vxl, vxr, vyl, vyr, vzl, vzr    ! Velocity components.
  real(wp) :: rl, rr, presl, presr        ! Primitive variables.
  real(wp) :: vnl, vnr                  ! Normal velocity
  real(wp) :: al, ar, h0l, h0r            ! Speed of sound, Total enthalpy
  real(wp) :: RT,rho,u,v,w,H,a,qn       ! Roe-averages
  real(wp) :: drho,dqn,dp               ! Wave strengths
  real(wp) :: du, dv, dw                ! Velocity conponent differences
  real(wp) :: SRp,SLm                   ! Wave speeds for the HLL part
  real(wp) :: nx1, ny1, nz1             ! Vector along which HLL is applied
  real(wp) :: nx2, ny2, nz2             ! Vector along which Roe is applied
  real(wp) :: alpha1, alpha2            ! Projections of the new normals
  real(wp) :: abs_dq                    ! Magnitude of the velocity difference
  real(wp) :: temp, tempx, tempy, tempz ! Temporary variables
  real(wp) :: ws1,ws2,ws3,ws4
  real(wp) :: eig1,eig2,eig3,eig4
  real(wp) :: LdU1,LdU2,LdU3,LdU4
  real(wp) :: diss_ec,diss_mx,diss_my,diss_mz,diss_ee
  real(wp) :: fl_ec,fl_mx,fl_my,fl_mz,fl_ee
  real(wp) :: fr_ec,fr_mx,fr_my,fr_mz,fr_ee
  real(wp) :: r1_ec,r1_mx,r1_my,r1_mz,r1_ee
  real(wp) :: r2_ec,r2_mx,r2_my,r2_mz,r2_ee
  real(wp) :: r3_ec,r3_mx,r3_my,r3_mz,r3_ee
  real(wp) :: r4_ec,r4_mx,r4_my,r4_mz,r4_ee
  !
continue
  !
  gm1 = gam - one
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,aave) &
  !$ACC&              PRESENT(htotl,htotr,lflux,rflux,return_value)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt( rnx*rnx + rny*rny + rnz*rnz )
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
    !
    ! Primitive and other variables.
    !
    ! Left state
    rl = face_cv(n,nec,1)
    rl_inv = one/rl
    vxl = face_cv(n,nmx,1)*rl_inv
    vyl = face_cv(n,nmy,1)*rl_inv
    vzl = face_cv(n,nmz,1)*rl_inv
    q2l = vxl*vxl + vyl*vyl + vzl*vzl
    presl = face_gp(n,npres,1)
    al = sqrt(max(gam*presl/rl,rndoff))
    h0l = al*al/gm1 + half*q2l
    vnl = vxl*unx + vyl*uny + vzl*unz
    !
    ! Right state
    rr = face_cv(n,nec,2)
    rr_inv = one/rr
    vxr = face_cv(n,nmx,2)*rr_inv
    vyr = face_cv(n,nmy,2)*rr_inv
    vzr = face_cv(n,nmz,2)*rr_inv
    q2r = vxr*vxr + vyr*vyr + vzr*vzr
    presr = face_gp(n,npres,2)
    ar = sqrt(max(gam*presr/rr,rndoff))
    h0r = ar*ar/gm1 + half*q2r
    vnr = vxr*unx + vyr*uny + vzr*unz
    !
    ! --------------------------------------------------------------------------
    ! Define n1 and n2, and compute alpha1 and alpha2: (4.2) in original paper.
    ! NOTE: n1 and n2 may need to be frozen at some point during
    !       a steady calculation to fully make it converge. For time-accurate
    !       calculation, it is not necessary.
    ! NOTE: For a boundary face, you may want to set
    !         (nx2,ny2,nz2)=(unx,uny,unz),  (nx1,ny1)=tangent vector,
    !       i.e., use the Roe flux.
    !
    abs_dq = sqrt( (vxr-vxl)**2 + (vyr-vyl)**2 + (vzr-vzl)**2 )
    !
    if (abs_dq > eps12) then
      !
      ! n1 = Velocity difference vector: normal to shock or tangent to shear
      !
      nx1 = (vxr-vxl)/abs_dq
      ny1 = (vyr-vyl)/abs_dq
      nz1 = (vzr-vzl)/abs_dq
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
      ! E.g., if n = (1,0,0), then, (0,-unz,uny) = (0,0,0). This is a problem.
      !       To avoid such a situation, we choose the 2D tangential vector
      !       having the maximum length.
      !
      tempx = uny*uny + unz*unz
      tempy = unz*unz + unx*unx
      tempz = unx*unx + uny*uny
      !
      if ( tempx >= tempy .and. tempx >= tempz ) then
        nx1 =  zero
        ny1 = -unz
        nz1 =  uny
      else if ( tempy >= tempx .and. tempy >= tempz ) then
        nx1 = -unz
        ny1 =  zero
        nz1 =  unx
      else if ( tempz >= tempx .and. tempz >= tempy ) then
        nx1 = -uny
        ny1 =  unx
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
    alpha1 = unx*nx1 + uny*ny1 + unz*nz1
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
    !     (nx2,ny2,nz2) = (n1xn)xn1 / |(n1xn)xn1|   ('x' is the vector product.)
    !
    ! (tempx,tempy,tempz) = n1xn
    tempx = ny1*unz - nz1*uny
    tempy = nz1*unx - nx1*unz
    tempz = nx1*uny - ny1*unx
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
    alpha2 = unx*nx2 + uny*ny2 + unz*nz2
    !
    ! Make alpha2 always positive.
    temp = sign(one,alpha2)
    nx2 = temp * nx2
    ny2 = temp * ny2
    nz2 = temp * nz2
    alpha2 = temp * alpha2
    !
    ! --------------------------------------------------------------------------
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
    RT = sqrt(rr/rl)
    rho = RT*rl                                ! Roe-ave density.
    u = (vxl + RT*vxr)/(one + RT)              ! Roe-ave x-velocity
    v = (vyl + RT*vyr)/(one + RT)              ! Roe-ave y-velocity
    w = (vzl + RT*vzr)/(one + RT)              ! Roe-ave z-velocity
    H = (h0l + RT*h0r)/(one + RT)              ! Roe-ave total enthalpy
    a = sqrt( gm1*(H-half*(u*u + v*v + w*w)) ) ! Roe-ave speed of sound
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
    vnl = vxl*nx1 + vyl*ny1 + vzl*nz1
    vnr = vxr*nx1 + vyr*ny1 + vzr*nz1
    SLm = min( zero, qn - a, vnl - al ) ! Minimum wave speed estimate
    SRp = max( zero, qn + a, vnr + ar ) ! Maximum wave speed estimate
    !
    ! This is the only place where n1=(nx1,ny1,nz1) is used.
    ! n1=(nx1,ny1,nz1) is never used below.
    ! ----------------------------------------------------
    !
    ! Wave Strengths
    !
    qn  = u *nx2 + v *ny2 + w *nz2
    vnl = vxl*nx2 + vyl*ny2 + vzl*nz2
    vnr = vxr*nx2 + vyr*ny2 + vzr*nz2
    !
    drho = rr - rl  ! Density difference
    dp =   presr - presl    ! Pressure difference
    dqn =  vnr - vnl   ! Normal velocity difference
    !
    LdU1 = (dp - rho*a*dqn )/(two*a*a) ! Left-moving acoustic wave strength
    LdU2 =  drho - dp/(a*a)            ! Entropy wave strength
    LdU3 = (dp + rho*a*dqn )/(two*a*a) ! Right-moving acoustic wave strength
    LdU4 = rho                 ! Shear wave strength (not really, just a factor)
    !
    ! Wave Speed (Eigenvalues)
    !
    eig1 = qn-a ! Left-moving acoustic wave velocity
    eig2 = qn   ! Entropy wave velocity
    eig3 = qn+a ! Right-moving acoustic wave velocity
    eig4 = qn   ! Shear wave velocity
    !
    ! Absolute values of the wave speeds (Eigenvalues)
    !
    ws1 = abs(qn-a) ! Left-moving acoustic wave speed
    ws2 = abs(qn)   ! Entropy wave speed
    ws3 = abs(qn+a) ! Right-moving acoustic wave speed
    ws4 = abs(qn)   ! Shear wave speed
    !
    ! Harten's Entropy Fix JCP(1983), 49, pp357-393: only for the
    ! nonlinear fields.
    ! NOTE: Avoids vanishing wave speeds by making a parabolic fit near ws = 0.
    !
    if ( ws1 < one5 ) then
      ws1 = half * ( ws1*ws1/one5 + one5 )
    end if
    !
    if ( ws3 < one5 ) then
      ws3 = half * ( ws3*ws3/one5 + one5 )
    end if
    !
    ! Combine wave speeds for Rotated-RHLL: Eq.(5.12) in original JCP2008 paper.
    !
    ws1 = alpha2*ws1 - (alpha1*two*SRp*SLm + alpha2*(SRp+SLm)*eig1)/(SRp-SLm)
    ws2 = alpha2*ws2 - (alpha1*two*SRp*SLm + alpha2*(SRp+SLm)*eig2)/(SRp-SLm)
    ws3 = alpha2*ws3 - (alpha1*two*SRp*SLm + alpha2*(SRp+SLm)*eig3)/(SRp-SLm)
    ws4 = alpha2*ws4 - (alpha1*two*SRp*SLm + alpha2*(SRp+SLm)*eig4)/(SRp-SLm)
    !
    ! Below, we compute the Roe dissipation term in the direction n2 with the
    ! above modified wave speeds. HLL wave speeds act something like the entropy
    ! fix or eigenvalue limiting; they contribute only by the amount given by
    ! the fraction, alpha1 (less than or equal to 1.0). See JCP2008 paper.
    !
    ! Right Eigenvectors:
    ! NOTE: Two shear wave components are combined into one, so that tangent
    !       vectors are not required, thats why there are only 4 vectors here.
    !
    ! Left-moving acoustic wave
    r1_ec = one
    r1_mx = u - a*nx2
    r1_my = v - a*ny2
    r1_mz = w - a*nz2
    r1_ee = H - a*qn
    !
    ! Entropy wave
    r2_ec = one
    r2_mx = u
    r2_my = v
    r2_mz = w
    r2_ee = half*(u*u + v*v + w*w)
    !
    ! Right-moving acoustic wave
    r3_ec = one
    r3_mx = u + a*nx2
    r3_my = v + a*ny2
    r3_mz = w + a*nz2
    r3_ee = H + a*qn
    !
    ! Two shear wave components combined into one (wave strength incorporated).
    du = vxr - vxl
    dv = vyr - vyl
    dw = vzr - vzl
    !
    r4_ec = zero
    r4_mx = du - dqn*nx2
    r4_my = dv - dqn*ny2
    r4_mz = dw - dqn*nz2
    r4_ee = u*du + v*dv + w*dw - qn*dqn
    !
    ! Dissipation Term: Roe dissipation with the modified wave speeds.
    !    |An|dU = R|Lambda|L*dU
    !           = sum_k of [ ws(k) * R(:,k) * L*dU(k) ], where n=n2.
    !
    ws1 = ws1 * LdU1
    ws2 = ws2 * LdU2
    ws3 = ws3 * LdU3
    ws4 = ws4 * LdU4
    !
    diss_ec = ws1*r1_ec + ws2*r2_ec + ws3*r3_ec + ws4*r4_ec
    diss_mx = ws1*r1_mx + ws2*r2_mx + ws3*r3_mx + ws4*r4_mx
    diss_my = ws1*r1_my + ws2*r2_my + ws3*r3_my + ws4*r4_my
    diss_mz = ws1*r1_mz + ws2*r2_mz + ws3*r3_mz + ws4*r4_mz
    diss_ee = ws1*r1_ee + ws2*r2_ee + ws3*r3_ee + ws4*r4_ee
    !
    ! Compute the physical flux: fL = Fn(vxl) and fR = Fn(vxr)
    !
    fl_ec = rl*vnl
    fl_mx = rl*vnl * vxl + presl*unx
    fl_my = rl*vnl * vyl + presl*uny
    fl_mz = rl*vnl * vzl + presl*unz
    fl_ee = rl*vnl * h0l
    !
    fr_ec = rr*vnr
    fr_mx = rr*vnr * vxr + presr*unx
    fr_my = rr*vnr * vyr + presr*uny
    fr_mz = rr*vnr * vzr + presr*unz
    fr_ee = rr*vnr * h0r
    !
    ! Compute the Rotated-RHLL flux, it looks like HLL flux with Roe dissipation
    !
    return_value(n,nec) = (SRp*fl_ec - SLm*fr_ec)/(SRp-SLm) - half*diss_ec
    return_value(n,nmx) = (SRp*fl_mx - SLm*fr_mx)/(SRp-SLm) - half*diss_mx
    return_value(n,nmy) = (SRp*fl_my - SLm*fr_my)/(SRp-SLm) - half*diss_my
    return_value(n,nmz) = (SRp*fl_mz - SLm*fr_mz)/(SRp-SLm) - half*diss_mz
    return_value(n,nee) = (SRp*fl_ee - SLm*fr_ee)/(SRp-SLm) - half*diss_ee
    !
  end do
  !
  ! Normal max wave speed in the normal direction.
  !  wsn = abs(qn) + a
  !
end subroutine rotated_rhll_upwind_flux_gpu_3d
!
!###############################################################################
!
pure subroutine rotated_rhll_upwind_flux_gpu_2d(n1,n2,facenrm,face_cv,face_gp, &
                                                lflux,rflux,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp), intent(inout) :: lflux(:)
  real(wp), intent(inout) :: rflux(:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,ds,rnx,rny,unx,uny
  real(wp) :: rl_inv,rr_inv,q2l,q2r
  real(wp) :: nx, ny
  real(wp) :: nx1, ny1, nx2, ny2      ! Rotated normals, n1 and n2
  real(wp) :: tx, ty                  ! Tangent vector (taken as n1)
  real(wp) :: alpha1, alpha2          ! Projections of the new normals
  real(wp) :: vxl, vxr, vyl, vyr      ! Velocity components.
  real(wp) :: rl, rr, presl, presr      ! Primitive variables.
  real(wp) :: vnl, vnr, vtl, vtr      ! Normal and tagent velocities
  real(wp) :: al, ar, h0l, h0r          ! Speeds of sound and total enthalpy
  real(wp) :: RT,rho,vx,vy,H,a        ! Roe-averages
  real(wp) :: vn, vt                  ! Normal and tagent velocities (Roe-ave)
  real(wp) :: drho,dvx,dvy,dvn,dvt,dp ! Wave strenghs
  real(wp) :: abs_dq                  ! Magnitude of the velocity difference
  real(wp) :: SRp,SLm                 ! Wave speeds for the HLL part
  real(wp) :: temp
  real(wp) :: ws1,ws2,ws3,ws4
  real(wp) :: eig1,eig2,eig3,eig4
  real(wp) :: LdU1,LdU2,LdU3,LdU4
  real(wp) :: diss_ec,diss_mx,diss_my,diss_ee
  real(wp) :: fl_ec,fl_mx,fl_my,fl_ee
  real(wp) :: fr_ec,fr_mx,fr_my,fr_ee
  real(wp) :: r1_ec,r1_mx,r1_my,r1_ee
  real(wp) :: r2_ec,r2_mx,r2_my,r2_ee
  real(wp) :: r3_ec,r3_mx,r3_my,r3_ee
  real(wp) :: r4_ec,r4_mx,r4_my,r4_ee
  !
continue
  !
  gm1 = gam - one
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,aave) &
  !$ACC&              PRESENT(htotl,htotr,lflux,rflux,return_value)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt( rnx*rnx + rny*rny )
    unx = rnx/ds
    uny = rny/ds
    !
    ! Left state
    !
    rl = face_cv(n,nec,1)
    rl_inv = one/rl
    vxl = face_cv(n,nmx,1)*rl_inv
    vyl = face_cv(n,nmy,1)*rl_inv
    q2l = vxl*vxl + vyl*vyl
    presl = face_gp(n,npres,1)
    al = sqrt(max(gam*presl/rl,rndoff))
    h0l = al*al/gm1 + half*q2l
    vnl = vxl*unx + vyl*uny
    !
    ! Right state
    !
    rr = face_cv(n,nec,2)
    vxr = face_cv(n,nmx,2)*rr_inv
    vyr = face_cv(n,nmy,2)*rr_inv
    q2r = vxr*vxr + vyr*vyr
    presr = face_gp(n,npres,2)
    ar = sqrt(max(gam*presr/rr,rndoff))
    h0r = ar*ar/gm1 + half*q2r
    vnr = vxr*unx + vyr*uny
    !
    ! Define n1 and n2, and compute alpha1 and alpha2: (4.2) in original paper.
    ! (NB: n1 and n2 may need to be frozen at some point during
    !      a steady calculation to fully make it converge. For time-accurate
    !      calculation, this is fine.)
    ! NB: For a boundary face, set (nx2,ny2)=(unx,uny), (nx1,ny1)=(-uny,unx).
    !
    abs_dq = sqrt( (vxr-vxl)**2+(vyr-vyl)**2 )
    !
    if (abs_dq > eps12) then
      nx1 = (vxr-vxl)/abs_dq
      ny1 = (vyr-vyl)/abs_dq
    else
      nx1 = -uny
      ny1 =  unx
    end if
    !
    alpha1 = unx * nx1 + uny * ny1
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
    alpha2 = unx * nx2 + uny * ny2
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
    RT = sqrt(rr/rl)
    rho = RT*rl
    vx = (vxl+RT*vxr)/(one+RT)
    vy = (vyl+RT*vyr)/(one+RT)
    H = ( h0l+RT* h0r)/(one+RT)
    a = sqrt( gm1*(H-half*(vx*vx+vy*vy)) )
    vn = vx*nx2+vy*ny2
    vt = vx*nx1+vy*ny1
    !
    ! Wave Strengths (remember that n2 is the normal and n1 is the tangent.)
    !
    vnl = vxl*nx2 + vyl*ny2
    vnr = vxr*nx2 + vyr*ny2
    vtl = vxl*nx1 + vyl*ny1
    vtr = vxr*nx1 + vyr*ny1
    !
    drho = rr - rl
    dp =   presr - presl
    dvn =  vnr - vnl
    dvt =  vtr - vtl
    !
    LdU1 = (dp - rho*a*dvn )/(two*a*a) ! Left-moving acoustic wave strength
    LdU2 =  rho*dvt/a          ! Shear wave strength (not really, just a factor)
    LdU3 =  drho - dp/(a*a)            ! Entropy wave strength
    LdU4 = (dp + rho*a*dvn )/(two*a*a) ! Right-moving acoustic wave strength
    !
    ! Wave Speeds for Roe flux part.
    !
    eig1 = vn-a
    eig2 = vn
    eig3 = vn
    eig4 = vn+a
    !
    ws1 = abs(vn-a)
    ws2 = abs(vn)
    ws3 = abs(vn)
    ws4 = abs(vn+a)
    !
    ! Harten's Entropy Fix JCP(1983), 49, pp357-393:
    ! only for the nonlinear fields.
    !
    if ( ws1 < one5 ) then
      ws1 = half * ( ws1*ws1/one5 + one5 )
    end if
    !
    if ( ws4 < one5 ) then
      ws4 = half * ( ws4*ws4/one5 + one5 )
    end if
    !
    ! HLL wave speeds, evaluated with [nx1,ny1] (=tangent wrt n2).
    !
    SRp = max( zero, vtr + ar, vt + a)
    SLm = min( zero, vtl - al, vt - a)
    !
    ! Modified wave speeds for the Rotated-RHLL flux: (5.12) in original paper.
    !
    ws1 = alpha2*ws1 - (two*alpha1*SRp*SLm + alpha2*(SRp+SLm)*eig1)/(SRp-SLm)
    ws2 = alpha2*ws2 - (two*alpha1*SRp*SLm + alpha2*(SRp+SLm)*eig2)/(SRp-SLm)
    ws3 = alpha2*ws3 - (two*alpha1*SRp*SLm + alpha2*(SRp+SLm)*eig3)/(SRp-SLm)
    ws4 = alpha2*ws4 - (two*alpha1*SRp*SLm + alpha2*(SRp+SLm)*eig4)/(SRp-SLm)
    !
    ! Right Eigenvectors: with n2 as normal and n1 as tangent.
    !
    r1_ec = one
    r1_mx = vx - a*nx2
    r1_my = vy - a*ny2
    r1_ee =  H - vn*a
    !
    r2_ec = zero
    r2_mx = a*nx1
    r2_my = a*ny1
    r2_ee = a*vt
    !
    r3_ec = one
    r3_mx = vx
    r3_my = vy
    r3_ee = half*(vx*vx+vy*vy)
    !
    r4_ec = one
    r4_mx = vx + a*nx2
    r4_my = vy + a*ny2
    r4_ee =  H + vn*a
    !
    ! Dissipation Term: Roe dissipation with the modified wave speeds.
    !
    ws1 = ws1 * LdU1
    ws2 = ws2 * LdU2
    ws3 = ws3 * LdU3
    ws4 = ws4 * LdU4
    !
    diss_ec = ws1*r1_ec + ws2*r2_ec + ws3*r3_ec + ws4*r4_ec
    diss_mx = ws1*r1_mx + ws2*r2_mx + ws3*r3_mx + ws4*r4_mx
    diss_my = ws1*r1_my + ws2*r2_my + ws3*r3_my + ws4*r4_my
    diss_ee = ws1*r1_ee + ws2*r2_ee + ws3*r3_ee + ws4*r4_ee
    !
    ! Compute the flux.
    !
    fl_ec = rl*vnl
    fl_mx = rl*vnl * vxl + presl*unx
    fl_my = rl*vnl * vyl + presl*uny
    fl_ee = rl*vnl *  h0l
    !
    fr_ec = rr*vnr
    fr_mx = rr*vnr * vxr + presr*unx
    fr_my = rr*vnr * vyr + presr*uny
    fr_ee = rr*vnr *  h0r
    !
    ! Compute the Rotated-RHLL flux.
    !
    return_value(n,nec) = (SRp*fl_ec - SLm*fr_ec)/(SRp-SLm) - half*diss_ec
    return_value(n,nmx) = (SRp*fl_mx - SLm*fr_mx)/(SRp-SLm) - half*diss_mx
    return_value(n,nmy) = (SRp*fl_my - SLm*fr_my)/(SRp-SLm) - half*diss_my
    return_value(n,nee) = (SRp*fl_ee - SLm*fr_ee)/(SRp-SLm) - half*diss_ee
    !
  end do
  !
end subroutine rotated_rhll_upwind_flux_gpu_2d
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
  real(wp) :: rtd,htd,atd2,atd,vntd
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
  ! Roe averaged variables
  !
  rtd   = sqrt(rr*rl)
  betal = rl / (rl+rtd)
  betar = one - betal
  veltd = betal*vell + betar*velr
  htd   = betal*hl + betar*hr
  atd2  = max(gm1*(htd - half*selfdot(veltd)),zero)
  atd   = sqrt(atd2)
  vntd  = dot_product(veltd,unrm)
  !
  atd2nz = max(atd2,rndoff)
  atdnz = max(atd,rndoff)
  atd2rt = atd2 / atd2nz
  atdrt = atd / atdnz
  !
  if (vntd >= zero) then
    !
    dw1  = half*((presr-presl)/atd2nz*atd2rt - (qr-ql)*rtd/atdnz*atdrt)
    hsr1 = -one4*gp1*dw1*atd/rtd
    eta1 = max(-abs(vntd-atd) + hsr1, zero)
    udw1 = dw1*(min(vntd-atd,zero) - half*eta1)
    rql  = rl*ql
    !
    return_value(nec)     =      rql + udw1
    return_value(nmb:nme) = vell*rql + udw1*(veltd-atd*unrm) + presl*unrm
    return_value(nee)     =   hl*rql + udw1*(  htd-atd*vntd)
    !
  else
    !
    dw4  = half*((presr-presl)/atd2nz*atd2rt + (qr-ql)*rtd/atdnz*atdrt)
    hsr4 = one4*gp1*dw4*atd/rtd
    eta4 = max(-abs(vntd+atd) + hsr4, zero)
    udw4 = dw4*(max(vntd+atd,zero) + half*eta4)
    rqr  = rr*qr
    !
    return_value(nec)     =      rqr - udw4
    return_value(nmb:nme) = velr*rqr - udw4*(veltd+atd*unrm) + presr*unrm
    return_value(nee)     =   hr*rqr - udw4*(  htd+atd*vntd)
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
  ! Roe averaged variables
  !
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
    pstar = rl*dsl*(sm-vnl) + presl
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
  end if
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
!
end module upwind_flux_mod
