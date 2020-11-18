module visc_flux_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  use eqn_idx, only : ncs,ntemp,nvisc,ncond,ng
  !
  implicit none
  !
  private
  !
  public :: visc_flx_cv_gpu
  public :: visc_wall_flx_gpu
  public :: visc_flx_cv
  public :: visc_wall_flx_cv
  public :: common_vflx_ave_fluxes_gpu
  public :: common_vflx_ave_sol_and_grad_gpu
  public :: common_visc_flx_cv
  !
contains
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
#ifndef TIME_MINIAPPS_FOR_GPU_COMPARISON
pure &
#endif
subroutine visc_flx_cv_gpu(n1,n2,cv,gp,dcvdx,return_value,output_visc_fluxes)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  real(wp),    intent(in) :: dcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  logical(lk), optional, intent(in) :: output_visc_fluxes
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = size(return_value,dim=2)
  !
  if (nr == 3) then
    call visc_flx_cv_gpu_3d(n1,n2,cv,gp,dcvdx,return_value)
  else if (nr == 2) then
    call visc_flx_cv_gpu_2d(n1,n2,cv,gp,dcvdx,return_value,output_visc_fluxes)
  else
    call visc_flx_cv_gpu_1d(n1,n2,cv,gp,dcvdx,return_value)
  end if
  !
end subroutine visc_flx_cv_gpu
!
!###############################################################################
!
#ifndef TIME_MINIAPPS_FOR_GPU_COMPARISON
pure &
#endif
subroutine visc_flx_cv_gpu_3d(n1,n2,cv,gp,dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar,       only : schmidt_lam
  use gpu_mod,    only : grad_cv_to_grad_pv_gpu
  use thermo_mod, only : compute_species_enthalpy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  real(wp),    intent(in) :: dcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,np
  real(wp) :: mvsc,mcdt
  real(wp) :: rho_inv,ave_ext_rate
  real(wp) :: dudx,dudy,dudz
  real(wp) :: dvdx,dvdy,dvdz
  real(wp) :: dwdx,dwdy,dwdz
  real(wp) :: tauxx,tauyy,tauzz
  real(wp) :: tauxy,tauxz,tauyz
  real(wp) :: val,hi,dsdx,dsdy,dsdz
  real(wp) :: schmidt_lam_inv
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: hsi(:,:)
  real(wp), allocatable :: dpvdx(:,:,:)
  !
continue
  !
  np = size(return_value,dim=1)
  !
  allocate ( dpvdx(1:np,1:3,1:nq) )
  !$ACC ENTER DATA CREATE(dpvdx)
  !
  ! Compute the gradients of the primitive variables
  !
  call grad_cv_to_grad_pv_gpu(n1,n2,cv,gp,dcvdx,dpvdx)
  !
  !$ACC PARALLEL LOOP PRESENT(cv,gp,return_value,dpvdx)
  do n = n1,n2
    !
    mvsc = gp(n,nvisc)
    mcdt = gp(n,ncond)
    !
    dudx = mvsc * dpvdx(n,1,nmx)
    dudy = mvsc * dpvdx(n,2,nmx)
    dudz = mvsc * dpvdx(n,3,nmx)
    !
    dvdx = mvsc * dpvdx(n,1,nmy)
    dvdy = mvsc * dpvdx(n,2,nmy)
    dvdz = mvsc * dpvdx(n,3,nmy)
    !
    dwdx = mvsc * dpvdx(n,1,nmz)
    dwdy = mvsc * dpvdx(n,2,nmz)
    dwdz = mvsc * dpvdx(n,3,nmz)
    !
    ave_ext_rate = two3 * (dudx + dvdy + dwdz)
    !
    tauxx = dudx + dudx - ave_ext_rate
    tauyy = dvdy + dvdy - ave_ext_rate
    tauzz = dwdz + dwdz - ave_ext_rate
    !
    tauxy = dvdx + dudy
    tauxz = dwdx + dudz
    tauyz = dwdy + dvdz
    !
    ! Density flux is zero
    !
    return_value(n,1,nec) = zero
    return_value(n,2,nec) = zero
    return_value(n,3,nec) = zero
    !
    ! Momentum equations
    !
    return_value(n,1,nmx) = tauxx
    return_value(n,2,nmx) = tauxy
    return_value(n,3,nmx) = tauxz
    !
    return_value(n,1,nmy) = tauxy
    return_value(n,2,nmy) = tauyy
    return_value(n,3,nmy) = tauyz
    !
    return_value(n,1,nmz) = tauxz
    return_value(n,2,nmz) = tauyz
    return_value(n,3,nmz) = tauzz
    !
    ! Energy Equation
    !
    rho_inv = one/cv(n,nec)
    !
    return_value(n,1,nee) = mcdt*dpvdx(n,1,nee) + rho_inv*cv(n,nmx)*tauxx &
                                                + rho_inv*cv(n,nmy)*tauxy &
                                                + rho_inv*cv(n,nmz)*tauxz
    !
    return_value(n,2,nee) = mcdt*dpvdx(n,2,nee) + rho_inv*cv(n,nmx)*tauxy &
                                                + rho_inv*cv(n,nmy)*tauyy &
                                                + rho_inv*cv(n,nmz)*tauyz
    !
    return_value(n,3,nee) = mcdt*dpvdx(n,3,nee) + rho_inv*cv(n,nmx)*tauxz &
                                                + rho_inv*cv(n,nmy)*tauyz &
                                                + rho_inv*cv(n,nmz)*tauzz
    !
  end do
  !
  if (ncs > 0) then
    !
    allocate ( hsi(1:np,1:ncs) )
    !$ACC ENTER DATA CREATE(hsi)
    !
    call compute_species_enthalpy(n1,n2,gp(:,ntemp),hsi)
    !
    schmidt_lam_inv = one / schmidt_lam
    !
    do m = 1,ncs
      !
      !$ACC PARALLEL LOOP PRESENT(hsi,gp,dpvdx,return_value)
      do n = n1,n2
        !
        val = gp(n,nvisc) * schmidt_lam_inv
        hi = hsi(n,m)
        !
        dsdx = val*dpvdx(n,1,m)
        dsdy = val*dpvdx(n,2,m)
        dsdz = val*dpvdx(n,3,m)
        !
        return_value(n,1,m) = dsdx
        return_value(n,2,m) = dsdy
        return_value(n,3,m) = dsdz
        !
        return_value(n,1,nee) = return_value(n,1,nee) + hi*dsdx
        return_value(n,2,nee) = return_value(n,2,nee) + hi*dsdy
        return_value(n,3,nee) = return_value(n,3,nee) + hi*dsdz
        !
      end do
      !
    end do
    !
    !$ACC EXIT DATA DELETE(hsi)
    deallocate ( hsi )
    !
  end if
  !
  !$ACC EXIT DATA DELETE(dpvdx)
  deallocate ( dpvdx )
  !
end subroutine visc_flx_cv_gpu_3d
!
!###############################################################################
!
#ifndef TIME_MINIAPPS_FOR_GPU_COMPARISON
pure &
#endif
subroutine visc_flx_cv_gpu_2d(n1,n2,cv,gp,dcvdx,return_value,output_visc_fluxes)
  !
  !.. Use Statements ..
  use ovar,       only : schmidt_lam
  use gpu_mod,    only : grad_cv_to_grad_pv_gpu
  use thermo_mod, only : compute_species_enthalpy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  real(wp),    intent(in) :: dcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  logical(lk), optional, intent(in) :: output_visc_fluxes
  !
  !.. Local Scalars ..
  integer  :: m,n,np
  real(wp) :: mvsc,mcdt
  real(wp) :: rho_inv,ave_ext_rate
  real(wp) :: dudx,dudy
  real(wp) :: dvdx,dvdy
  real(wp) :: tauxx,tauyy,tauxy
  real(wp) :: val,hi,dsdx,dsdy
  real(wp) :: schmidt_lam_inv
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: hsi(:,:)
  real(wp), allocatable :: dpvdx(:,:,:)
  real(wp), allocatable :: mixture_flux(:,:)
  !
continue
  !
  np = size(return_value,dim=1)
  !
  allocate ( dpvdx(1:np,1:2,1:nq) )
  !$ACC ENTER DATA CREATE(dpvdx)
  !
  ! Compute the gradients of the primitive variables
  !
  call grad_cv_to_grad_pv_gpu(n1,n2,cv,gp,dcvdx,dpvdx)
  !
  !$ACC PARALLEL LOOP PRESENT(cv,gp,return_value,dpvdx)
  do n = n1,n2
    !
    mvsc = gp(n,nvisc)
    mcdt = gp(n,ncond)
    !
    dudx = mvsc * dpvdx(n,1,nmx)
    dudy = mvsc * dpvdx(n,2,nmx)
    !
    dvdx = mvsc * dpvdx(n,1,nmy)
    dvdy = mvsc * dpvdx(n,2,nmy)
    !
    ave_ext_rate = two3 * (dudx + dvdy)
    !
    tauxx = dudx + dudx - ave_ext_rate
    tauyy = dvdy + dvdy - ave_ext_rate
    tauxy = dvdx + dudy
    !
    ! Density flux is zero
    !
    return_value(n,1,nec) = zero
    return_value(n,2,nec) = zero
    !
    ! Momentum equations
    !
    return_value(n,1,nmx) = tauxx
    return_value(n,2,nmx) = tauxy
    !
    return_value(n,1,nmy) = tauxy
    return_value(n,2,nmy) = tauyy
    !
    ! Energy Equation
    !
    rho_inv = one/cv(n,nec)
    !
    return_value(n,1,nee) = mcdt*dpvdx(n,1,nee) + rho_inv*cv(n,nmx)*tauxx &
                                                + rho_inv*cv(n,nmy)*tauxy
    !
    return_value(n,2,nee) = mcdt*dpvdx(n,2,nee) + rho_inv*cv(n,nmx)*tauxy &
                                                + rho_inv*cv(n,nmy)*tauyy
    !
  end do
  !
  if (present(output_visc_fluxes)) then
    if (output_visc_fluxes) then
      block
        !
        use ovar, only : itcur
        use geovar, only : xyz
        integer :: io,l,m,n
        character(len=50) :: fname
        !
        write (fname,1) itcur
        open (newunit=io,file=fname,status="replace",action="write")
        do n = 1,size(return_value,dim=1)
          write (io,2) n, (xyz(l,n), l=1,2), &
                         ((return_value(n,l,m), l=1,2), m=nec,nee)
        end do
        close (io)
        !
        1 format ("base_visc_fluxes_at_sp.",i0,".dat")
        2 format (i7,2x,10(2es14.6,:,","))
        !
      end block
    end if
  end if
  !
  if (ncs > 0) then
    !
    allocate ( mixture_flux(1:np,1:2) , source=zero )
    !$ACC ENTER DATA CREATE(mixture_flux)
    allocate ( hsi(1:np,1:ncs) )
    !$ACC ENTER DATA CREATE(hsi)
    !
    call compute_species_enthalpy(n1,n2,gp(:,ntemp),hsi)
    !
    schmidt_lam_inv = one / schmidt_lam
    !
    do m = 1,ncs
      !
      !$ACC PARALLEL LOOP PRESENT(hsi,gp,dpvdx,return_value)
      do n = n1,n2
        !
        val = gp(n,nvisc) * schmidt_lam_inv
        hi = hsi(n,m)
        !
        dsdx = val*dpvdx(n,1,m)
        dsdy = val*dpvdx(n,2,m)
        !
        return_value(n,1,m) = dsdx
        return_value(n,2,m) = dsdy
        !
        mixture_flux(n,1) = mixture_flux(n,1) + hi*dsdx
        mixture_flux(n,2) = mixture_flux(n,2) + hi*dsdx
        !
       !return_value(n,1,nee) = return_value(n,1,nee) + hi*dsdx
       !return_value(n,2,nee) = return_value(n,2,nee) + hi*dsdy
        !
      end do
      !
    end do
    !
    !$ACC EXIT DATA DELETE(hsi)
    deallocate ( hsi )
    !
    if (present(output_visc_fluxes)) then
      if (output_visc_fluxes) then
        block
          !
          use ovar, only : itcur
          use geovar, only : xyz
          integer :: io,l,m,n
          character(len=50) :: fname
          !
          write (fname,3) itcur
          open (newunit=io,file=fname,status="replace",action="write")
          do n = 1,size(return_value,dim=1)
            write (io,4) n, (xyz(l,n), l=1,2), &
                           ((return_value(n,l,m), l=1,2), m=1,ncs), &
                            (mixture_flux(n,l), l=1,2)
          end do
          close (io)
          3 format ("species_visc_fluxes_at_sp.",i0,".dat")
          4 format (i7,2x,10(2es14.6,:,","))
          !
        end block
      end if
    end if
    !
    !$ACC PARALLEL LOOP PRESENT(return_value,mixture_flux)
    do n = n1,n2
      return_value(n,1,nee) = return_value(n,1,nee) + mixture_flux(n,1)
      return_value(n,2,nee) = return_value(n,2,nee) + mixture_flux(n,2)
    end do
    !
    !$ACC EXIT DATA DELETE(mixture_flux)
    deallocate ( mixture_flux )
    !
  end if
  !
  !$ACC EXIT DATA DELETE(dpvdx)
  deallocate ( dpvdx )
  !
end subroutine visc_flx_cv_gpu_2d
!
!###############################################################################
!
#ifndef TIME_MINIAPPS_FOR_GPU_COMPARISON
pure &
#endif
subroutine visc_flx_cv_gpu_1d(n1,n2,cv,gp,dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar,       only : schmidt_lam
  use gpu_mod,    only : grad_cv_to_grad_pv_gpu
  use thermo_mod, only : compute_species_enthalpy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv(:,:)
  real(wp),    intent(in) :: gp(:,:)
  real(wp),    intent(in) :: dcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,np
  real(wp) :: tauxx,dsdx
  real(wp) :: schmidt_lam_inv
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: hsi(:,:)
  real(wp), allocatable :: dpvdx(:,:,:)
  !
continue
  !
  np = size(return_value,dim=1)
  !
  allocate ( dpvdx(1:np,1:1,1:nq) )
  !$ACC ENTER DATA CREATE(dpvdx)
  !
  ! Get the gradients of the primitive variables using the
  ! gradients of the conservative variables and the product rule
  !
  call grad_cv_to_grad_pv_gpu(n1,n2,cv,gp,dcvdx,dpvdx)
  !
  ! Density flux is zero
  !
  !$ACC PARALLEL LOOP PRESENT(cv,gp,return_value,dpvdx)
  do n = n1,n2
    !
    tauxx = two * two3 * gp(n,nvisc) * dpvdx(n,1,nmx)
    !
    ! Density flux is zero
    !
    return_value(n,1,nec) = zero
    !
    ! Momentum equation
    !
    return_value(n,1,nmx) = tauxx
    !
    ! Energy Equation
    !
    return_value(n,1,nee) = gp(n,ncond)*dpvdx(n,1,nee) &
                          + tauxx*cv(n,nmx)/cv(n,nec)
    !
  end do
  !
  if (ncs > 0) then
    !
    allocate ( hsi(1:np,1:ncs) )
    !$ACC ENTER DATA CREATE(hsi)
    !
    call compute_species_enthalpy(n1,n2,gp(:,ntemp),hsi)
    !
    schmidt_lam_inv = one / schmidt_lam
    !
    do m = 1,ncs
      !
      !$ACC PARALLEL LOOP PRESENT(hsi,gp,dpvdx,return_value)
      do n = n1,n2
        !
        dsdx = gp(n,nvisc)*schmidt_lam_inv*dpvdx(n,1,m)
        !
        return_value(n,1,m) = dsdx
        !
        return_value(n,1,nee) = return_value(n,1,nee) + hsi(n,m)*dsdx
        !
      end do
      !
    end do
    !
    !$ACC EXIT DATA DELETE(hsi)
    deallocate ( hsi )
    !
  end if
  !
  !$ACC EXIT DATA DELETE(dpvdx)
  deallocate ( dpvdx )
  !
end subroutine visc_flx_cv_gpu_1d
!
!###############################################################################
!
#ifndef TIME_MINIAPPS_FOR_GPU_COMPARISON
pure &
#endif
function visc_flx_cv1(uu,gp,dcvdx) result(return_value)
  !
  !.. Use Statements ..
  use geovar,  only : nr
  use gpu_mod, only : grad_cv_to_grad_pv_gpu
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nq),      intent(in) :: uu
  real(wp), dimension(1:ng),      intent(in) :: gp
  real(wp), dimension(1:nr,1:nq), intent(in) :: dcvdx
  !
  !.. Function Result ..
  real(wp), dimension(1:nr,1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,l
  real(wp) :: mvsc,mcdt
  real(wp) :: ave_ext_rate
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: vel
  real(wp), dimension(1:nr,1:nq) :: dpvdx
  !
continue
  !
  ! Initialize the viscous fluxes to zero.
  ! NOTE: This takes care of the continuity equation
  !
  return_value(1:nr,1:nq) = zero
  !
  ! Get the gradients of the primitive variables using the
  ! gradients of the conservative variables and the product rule
  !
  call grad_cv_to_grad_pv_gpu(uu,gp,dcvdx,dpvdx)
  !
  ! Get the velocity vector
  !
  vel(:) = uu(nmb:nme)/uu(nec)
  !
  ! Compute the molecular viscosity and molecular conductivity
  !
  mvsc = gp(nvisc)
  mcdt = gp(ncond)
  !
  ! Compute the average extension rate
  !
  ave_ext_rate = zero
  do l = 1,nr
    ave_ext_rate = ave_ext_rate + dpvdx(l,nmb-1+l)
  end do
  ave_ext_rate = two3 * ave_ext_rate
  !
  ! Compute the viscous fluxes for the momentum equations which are just
  ! the components of the viscous stress tensor.
  ! NOTE: We are initially excluding the molecular viscosity term which
  !       will be multiplied to the whole tau array after these loops.
  !
  do i = 1,nr
    do j = 1,nr
        return_value(j,nmb-1+i) = dpvdx(i,nmb-1+j) + dpvdx(j,nmb-1+i)
    end do
  end do
  !
  ! Subtract the average extension rate from the normal stresses
  !
  do l = 1,nr
    return_value(l,nmb-1+l) = return_value(l,nmb-1+l) - ave_ext_rate
  end do
  !
  ! Multiply the viscous stresses by the molecular viscosity
  !
  return_value(1:nr,nmb:nme) = mvsc * return_value(1:nr,nmb:nme)
  !
  ! Compute the viscous fluxes for the energy equation
  !
  do l = 1,nr
    return_value(l,nq) = mcdt*dpvdx(l,nq) + dot( vel , return_value(l,nmb:nme) )
  end do
  !
end function visc_flx_cv1
!
!###############################################################################
!
#ifndef TIME_MINIAPPS_FOR_GPU_COMPARISON
pure &
#endif
subroutine visc_flx_cv(uu,gp,dcvdx,return_value)
! !$ACC ROUTINE SEQ
  !
  !.. Use Statements ..
  use geovar,  only : nr
  use gpu_mod, only : grad_cv_to_grad_pv_gpu
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nq),      intent(in) :: uu
  real(wp), dimension(1:ng),      intent(in) :: gp
  real(wp), dimension(1:nr,1:nq), intent(in) :: dcvdx
  !
  !.. Function Result ..
  real(wp), dimension(1:nq,1:nr), intent(inout) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,l
  real(wp) :: mvsc,mcdt
  real(wp) :: ave_ext_rate
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: vel
  real(wp), dimension(1:nr,1:nq) :: dpvdx
  !
continue
  !
  ! Initialize the viscous fluxes to zero.
  ! NOTE: This takes care of the continuity equation
  !
 !i = 0
  return_value(1:nq,1:nr) = zero
  !
  ! Get the gradients of the primitive variables using the
  ! gradients of the conservative variables and the product rule
  !
  call grad_cv_to_grad_pv_gpu(uu,gp,dcvdx,dpvdx)
  !
  ! Get the velocity vector
  !
  vel(:) = uu(nmb:nme)/uu(nec)
  !
  ! Compute the molecular viscosity and molecular conductivity
  !
  mvsc = gp(nvisc)
  mcdt = gp(ncond)
  !
  ! Compute the average extension rate
  !
  ave_ext_rate = zero
  do l = 1,nr
    ave_ext_rate = ave_ext_rate + dpvdx(l,nmb-1+l)
  end do
  ave_ext_rate = two3 * ave_ext_rate
  !
  ! Compute the viscous fluxes for the momentum equations which are just
  ! the components of the viscous stress tensor.
  ! NOTE: We are initially excluding the molecular viscosity term which
  !       will be multiplied to the whole tau array after these loops.
  !
  do j = 1,nr
    do i = 1,nr
        return_value(nmb-1+i,j) = dpvdx(i,nmb-1+j) + dpvdx(j,nmb-1+i)
    end do
  end do
  !
  ! Subtract the average extension rate from the normal stresses
  !
  do l = 1,nr
    return_value(nmb-1+l,l) = return_value(nmb-1+l,l) - ave_ext_rate
  end do
  !
  ! Multiply the viscous stresses by the molecular viscosity
  !
  return_value(nmb:nme,1:nr) = mvsc * return_value(nmb:nme,1:nr)
  !
  ! Compute the viscous fluxes for the energy equation
  !
  do l = 1,nr
    return_value(nq,l) = mcdt*dpvdx(l,nq) + &
                         dot( vel(1:nr) , return_value(nmb:nme,l) )
  end do
  !
!end function visc_flx_cv
end subroutine visc_flx_cv
!
!###############################################################################
!
#ifndef TIME_MINIAPPS_FOR_GPU_COMPARISON
pure &
#endif
function normal_visc_flx_cv(rnrm,uu,gp,dcvdx) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nq),      intent(in) :: uu
  real(wp), dimension(1:ng),      intent(in) :: gp
  real(wp), dimension(1:nr),      intent(in) :: rnrm
  real(wp), dimension(1:nr,1:nq), intent(in) :: dcvdx
  !
  !.. Function Result ..
  real(wp), dimension(1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer :: l
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)      :: unrm
  real(wp), dimension(1:nq,1:nr) :: vflx
  !
continue
  !
  unrm(:) = unit_vector( rnrm(:) )
  !
 !vflx(1:nq,1:nr) = visc_flx_cv(uu,gp,dcvdx)
  call visc_flx_cv(uu,gp,dcvdx,vflx)
  !
  return_value(1:nq) = unrm(1)*vflx(1:nq,1)
  do l = 2,nr
    return_value(1:nq) = return_value(1:nq) + unrm(l)*vflx(1:nq,l)
  end do
  !
end function normal_visc_flx_cv
!
!###############################################################################
!
pure &
subroutine visc_wall_flx_gpu(n1,n2,facenrm,face_cv,face_gp,face_dcvdx, &
                             return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: face_dcvdx(:,:,:,:)
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
    call visc_wall_flx_gpu_3d(n1,n2,facenrm,face_cv,face_gp,face_dcvdx, &
                              return_value)
  else if (nr == 2) then
    call visc_wall_flx_gpu_2d(n1,n2,facenrm,face_cv,face_gp,face_dcvdx, &
                              return_value)
  else
    call visc_wall_flx_gpu_1d(n1,n2,facenrm,face_cv,face_gp,face_dcvdx, &
                              return_value)
  end if
  !
end subroutine visc_wall_flx_gpu
!
!###############################################################################
!
pure subroutine visc_wall_flx_gpu_3d(n1,n2,facenrm,face_cv,face_gp,face_dcvdx, &
                                     return_value)
  !
  !.. Use Statements ..
  use ovar,       only : schmidt_lam
  use thermo_mod, only : compute_species_enthalpy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: face_dcvdx(:,:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,np
  real(wp) :: rho_inv,mvsc,mcdt,ave_ext_rate,tflux
  real(wp) :: tauxx,tauyy,tauzz,tauxy,tauxz,tauyz
  real(wp) :: rnx,velx,dudx,dvdx,dwdx,dtdx,tflx_x,uflux
  real(wp) :: rny,vely,dudy,dvdy,dwdy,dtdy,tflx_y,vflux
  real(wp) :: rnz,velz,dudz,dvdz,dwdz,dtdz,tflx_z,wflux
  real(wp) :: val,dsdx,dsdy,dsdz,sflux,schmidt_lam_inv
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: hsi(:,:)
  !
continue
  !
  np = size(return_value,dim=1)
  !
  ! For slip wall: all gradients /= zero
  !                uwall /= zero
  !
  ! For adiabatic wall: Temperature gradient == zero
  !                     Velocity gradients /= zero
  !                     uwall(nmb:nme) [i.e., wall velocities] == zero
  !
  ! For isothermal wall: all gradients /= zero
  !                      uwall(nmb:nme) [i.e., wall velocities] == zero
  !
  ! Initialize the viscous fluxes to zero.
  ! NOTE: This takes care of the continuity equation
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,face_dcvdx,return_value)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    !
    rho_inv = one/face_cv(n,nec,2)
    !
    velx = face_cv(n,nmx,2)*rho_inv
    vely = face_cv(n,nmy,2)*rho_inv
    velz = face_cv(n,nmz,2)*rho_inv
    !
    mvsc = face_gp(n,nvisc,2)
    mcdt = face_gp(n,ncond,2)
    !
    ! NOTE: FOR WALL BCs, THE GRADIENTS STORED IN FACE_DCVDX(:,:,2)
    !       ARE ACTUALLY THE PRIMITVE VARIABLE GRADIENTS
    !
    dudx = mvsc * face_dcvdx(n,1,nmx,2)
    dudy = mvsc * face_dcvdx(n,2,nmx,2)
    dudz = mvsc * face_dcvdx(n,3,nmx,2)
    !
    dvdx = mvsc * face_dcvdx(n,1,nmy,2)
    dvdy = mvsc * face_dcvdx(n,2,nmy,2)
    dvdz = mvsc * face_dcvdx(n,3,nmy,2)
    !
    dwdx = mvsc * face_dcvdx(n,1,nmz,2)
    dwdy = mvsc * face_dcvdx(n,2,nmz,2)
    dwdz = mvsc * face_dcvdx(n,3,nmz,2)
    !
    dtdx = mcdt * face_dcvdx(n,1,nee,2)
    dtdy = mcdt * face_dcvdx(n,2,nee,2)
    dtdz = mcdt * face_dcvdx(n,3,nee,2)
    !
    ave_ext_rate = two3 * (dudx + dvdy + dwdz)
    !
    tauxx = dudx + dudx - ave_ext_rate
    tauyy = dvdy + dvdy - ave_ext_rate
    tauzz = dwdz + dwdz - ave_ext_rate
    !
    tauxy = dvdx + dudy
    tauxz = dwdx + dudz
    tauyz = dwdy + dvdz
    !
    tflx_x = dtdx + velx*tauxx + vely*tauxy + velz*tauxz
    tflx_y = dtdy + velx*tauxy + vely*tauyy + velz*tauyz
    tflx_z = dtdz + velx*tauxz + vely*tauyz + velz*tauzz
    !
    uflux = rnx*tauxx  + rny*tauxy  + rnz*tauxz
    vflux = rnx*tauxy  + rny*tauyy  + rnz*tauyz
    wflux = rnx*tauxz  + rny*tauyz  + rnz*tauzz
    tflux = rnx*tflx_x + rny*tflx_y + rnz*tflx_z
    !
    ! Subtract the normal viscous flux from the
    ! inviscid flux already stored in return_value
    !
    return_value(n,nmx) = return_value(n,nmx) - uflux
    return_value(n,nmy) = return_value(n,nmy) - vflux
    return_value(n,nmz) = return_value(n,nmz) - wflux
    return_value(n,nee) = return_value(n,nee) - tflux
    !
  end do
  !
  ! Compute the viscous fluxes for the species continuity equations
  ! and include their contributions to the energy equation.
  !
  if (ncs > 0) then
    !
    allocate ( hsi(1:np,1:ncs) )
    !$ACC ENTER DATA CREATE(hsi)
    !
    call compute_species_enthalpy(n1,n2,face_gp(:,ntemp,2),hsi)
    !
    schmidt_lam_inv = one / schmidt_lam
    !
    do m = 1,ncs
      !
      !$ACC PARALLEL LOOP PRESENT(facenrm,hsi,face_gp,face_dcvdx,return_value)
      do n = n1,n2
        !
        val = face_gp(n,nvisc,2) * schmidt_lam_inv
        !
        dsdx = val*face_dcvdx(n,1,m,2)*facenrm(1,n)
        dsdy = val*face_dcvdx(n,2,m,2)*facenrm(2,n)
        dsdz = val*face_dcvdx(n,3,m,2)*facenrm(3,n)
        !
        sflux = dsdx + dsdy + dsdz
        !
        return_value(n,m)   = return_value(n,m)   - sflux
        return_value(n,nee) = return_value(n,nee) - sflux*hsi(n,m)
        !
      end do
      !
    end do
    !
    !$ACC EXIT DATA DELETE(hsi)
    deallocate ( hsi )
    !
  end if
  !
end subroutine visc_wall_flx_gpu_3d
!
!###############################################################################
!
pure &
subroutine visc_wall_flx_gpu_2d(n1,n2,facenrm,face_cv,face_gp,face_dcvdx, &
                                return_value)
  !
  !.. Use Statements ..
  use ovar,       only : schmidt_lam
  use thermo_mod, only : compute_species_enthalpy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: face_dcvdx(:,:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,np
  real(wp) :: rho_inv,mvsc,mcdt,ave_ext_rate,tflux
  real(wp) :: tauxx,tauyy,tauxy
  real(wp) :: rnx,velx,dudx,dvdx,dwdx,dtdx,tflx_x,uflux
  real(wp) :: rny,vely,dudy,dvdy,dwdy,dtdy,tflx_y,vflux
  real(wp) :: val,dsdx,dsdy,sflux,schmidt_lam_inv
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: hsi(:,:)
  !
continue
  !
  np = size(return_value,dim=1)
  !
  ! For slip wall: all gradients /= zero
  !                uwall /= zero
  !
  ! For adiabatic wall: Temperature gradient == zero
  !                     Velocity gradients /= zero
  !                     uwall(nmb:nme) [i.e., wall velocities] == zero
  !
  ! For isothermal wall: all gradients /= zero
  !                      uwall(nmb:nme) [i.e., wall velocities] == zero
  !
  ! Initialize the viscous fluxes to zero.
  ! NOTE: This takes care of the continuity equation
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,face_dcvdx,return_value)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    !
    rho_inv = one/face_cv(n,nec,2)
    !
    velx = face_cv(n,nmx,2)*rho_inv
    vely = face_cv(n,nmy,2)*rho_inv
    !
    mvsc = face_gp(n,nvisc,2)
    mcdt = face_gp(n,ncond,2)
    !
    ! NOTE: FOR WALL BCs, THE GRADIENTS STORED IN FACE_DCVDX(:,:,2)
    !       ARE ACTUALLY THE PRIMITVE VARIABLE GRADIENTS
    !
    dudx = mvsc * face_dcvdx(n,1,nmx,2)
    dudy = mvsc * face_dcvdx(n,2,nmx,2)
    !
    dvdx = mvsc * face_dcvdx(n,1,nmy,2)
    dvdy = mvsc * face_dcvdx(n,2,nmy,2)
    !
    dtdx = mcdt * face_dcvdx(n,1,nee,2)
    dtdy = mcdt * face_dcvdx(n,2,nee,2)
    !
    ave_ext_rate = two3 * (dudx + dvdy)
    !
    tauxx = dudx + dudx - ave_ext_rate
    tauyy = dvdy + dvdy - ave_ext_rate
    tauxy = dvdx + dudy
    !
    tflx_x = dtdx + velx*tauxx + vely*tauxy
    tflx_y = dtdy + velx*tauxy + vely*tauyy
    !
    uflux = rnx*tauxx  + rny*tauxy
    vflux = rnx*tauxy  + rny*tauyy
    tflux = rnx*tflx_x + rny*tflx_y
    !
    ! Subtract the normal viscous flux from the
    ! inviscid flux already stored in return_value
    !
    return_value(n,nmx) = return_value(n,nmx) - uflux
    return_value(n,nmy) = return_value(n,nmy) - vflux
    return_value(n,nee) = return_value(n,nee) - tflux
    !
  end do
  !
  ! Compute the viscous fluxes for the species continuity equations
  ! and include their contributions to the energy equation.
  !
  if (ncs > 0) then
    !
    allocate ( hsi(1:np,1:ncs) )
    !$ACC ENTER DATA CREATE(hsi)
    !
    call compute_species_enthalpy(n1,n2,face_gp(:,ntemp,2),hsi)
    !
    schmidt_lam_inv = one / schmidt_lam
    !
    do m = 1,ncs
      !
      !$ACC PARALLEL LOOP PRESENT(facenrm,hsi,face_gp,face_dcvdx,return_value)
      do n = n1,n2
        !
        val = face_gp(n,nvisc,2) * schmidt_lam_inv
        !
        dsdx = val*face_dcvdx(n,1,m,2)*facenrm(1,n)
        dsdy = val*face_dcvdx(n,2,m,2)*facenrm(2,n)
        !
        sflux = dsdx + dsdy
        !
        return_value(n,m)   = return_value(n,m)   - sflux
        return_value(n,nee) = return_value(n,nee) - sflux*hsi(n,m)
        !
      end do
      !
    end do
    !
    !$ACC EXIT DATA DELETE(hsi)
    deallocate ( hsi )
    !
  end if
  !
end subroutine visc_wall_flx_gpu_2d
!
!###############################################################################
!
pure subroutine visc_wall_flx_gpu_1d(n1,n2,facenrm,face_cv,face_gp,face_dcvdx, &
                                     return_value)
  !
  !.. Use Statements ..
  use ovar,       only : schmidt_lam
  use thermo_mod, only : compute_species_enthalpy
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: face_dcvdx(:,:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,np
  real(wp) :: mvsc,mcdt,tflux,tauxx
  real(wp) :: rnx,velx,dudx,dtdx,tflx_x,uflux
  real(wp) :: val,sflux,schmidt_lam_inv
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: hsi(:,:)
  !
continue
  !
  np = size(return_value,dim=1)
  !
  ! For slip wall: all gradients /= zero
  !                uwall /= zero
  !
  ! For adiabatic wall: Temperature gradient == zero
  !                     Velocity gradients /= zero
  !                     uwall(nmb:nme) [i.e., wall velocities] == zero
  !
  ! For isothermal wall: all gradients /= zero
  !                      uwall(nmb:nme) [i.e., wall velocities] == zero
  !
  ! Initialize the viscous fluxes to zero.
  ! NOTE: This takes care of the continuity equation
  !
  !$ACC PARALLEL LOOP PRESENT(facenrm,face_cv,face_gp,face_dcvdx,return_value)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    !
    velx = face_cv(n,nmx,2)/face_cv(n,nec,2)
    !
    mvsc = face_gp(n,nvisc,2)
    mcdt = face_gp(n,ncond,2)
    !
    ! NOTE: FOR WALL BCs, THE GRADIENTS STORED IN FACE_DCVDX(:,:,2)
    !       ARE ACTUALLY THE PRIMITVE VARIABLE GRADIENTS
    !
    dudx = mvsc * face_dcvdx(n,1,nmx,2)
    dtdx = mcdt * face_dcvdx(n,1,nee,2)
    !
    tauxx = two * two3 * dudx
    !
    tflx_x = dtdx + velx*tauxx
    !
    uflux = rnx*tauxx
    tflux = rnx*tflx_x
    !
    ! Subtract the normal viscous flux from the
    ! inviscid flux already stored in return_value
    !
    return_value(n,nmx) = return_value(n,nmx) - uflux
    return_value(n,nee) = return_value(n,nee) - tflux
    !
  end do
  !
  ! Compute the viscous fluxes for the species continuity equations
  ! and include their contributions to the energy equation.
  !
  if (ncs > 0) then
    !
    allocate ( hsi(1:np,1:ncs) )
    !$ACC ENTER DATA CREATE(hsi)
    !
    call compute_species_enthalpy(n1,n2,face_gp(:,ntemp,2),hsi)
    !
    schmidt_lam_inv = one / schmidt_lam
    !
    do m = 1,ncs
      !
      !$ACC PARALLEL LOOP PRESENT(facenrm,hsi,face_gp,face_dcvdx,return_value)
      do n = n1,n2
        !
        val = face_gp(n,nvisc,2) * schmidt_lam_inv
        !
        sflux = val*face_dcvdx(n,1,m,2)*facenrm(1,n)
        !
        return_value(n,m)   = return_value(n,m)   - sflux
        return_value(n,nee) = return_value(n,nee) - sflux*hsi(n,m)
        !
      end do
      !
    end do
    !
    !$ACC EXIT DATA DELETE(hsi)
    deallocate ( hsi )
    !
  end if
  !
end subroutine visc_wall_flx_gpu_1d
!
!###############################################################################
!
pure subroutine visc_wall_flx_sub(unrm,uwall,gpwall,dpvdx,return_value)
! !$ACC ROUTINE SEQ
  !
  ! UNRM  : Outward pointing unit normal to the wall
  !         boundary face for current flux point
  ! UWALL : Conservative variables exact solution at the wall
  ! DPVDX : Prescribed wall gradients of the primitive variables at the
  !         current flux point on the wall boundary face. Most gradients
  !         are extrapolated from the interior cell flux point on the
  !         boundary except for certain special cases, e.g., an adiabatic
  !         wall where dT/dx = dT/dy = 0.
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),      intent(in) :: unrm
  real(wp), dimension(1:nq),      intent(in) :: uwall
  real(wp), dimension(1:ng),      intent(in) :: gpwall
  real(wp), dimension(1:nr,1:nq), intent(in) :: dpvdx
  real(wp), dimension(1:nq),   intent(inout) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,l,m
  real(wp) :: mvsc,mcdt
  real(wp) :: ave_ext_rate
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)      :: vel
  real(wp), dimension(1:nr,1:nq) :: vflx
  !
continue
  !
  ! For slip wall: all gradients /= zero
  !                uwall /= zero
  !
  ! For adiabatic wall: Temperature gradient == zero
  !                     Velocity gradients /= zero
  !                     uwall(nmb:nme) [i.e., wall velocities] == zero
  !
  ! For isothermal wall: all gradients /= zero
  !                      uwall(nmb:nme) [i.e., wall velocities] == zero
  !
  ! Initialize the viscous fluxes to zero.
  ! NOTE: This takes care of the continuity equation
  !
  vflx(:,:) = zero
  !
  ! Get the velocity vector
  !
  vel(:) = uwall(nmb:nme)/uwall(nec)
  !
  ! Compute the molecular viscosity and molecular conductivity
  !
  mvsc = gpwall(nvisc)
  mcdt = gpwall(ncond)
  !
  ! Compute the average extension rate
  !
  ave_ext_rate = zero
  do l = 1,nr
    ave_ext_rate = ave_ext_rate + dpvdx(l,nmb-1+l)
  end do
  ave_ext_rate = two3 * ave_ext_rate
  !
  ! Compute the viscous fluxes for the momentum equations which are just
  ! the components of the viscous stress tensor.
  ! NOTE: We are initially excluding the molecular viscosity term which
  !       will be multiplied to the whole tau array after these loops.
  !
  do j = 1,nr
    do i = 1,nr
        vflx(i,nmb-1+j) = dpvdx(i,nmb-1+j) + dpvdx(j,nmb-1+i)
    end do
  end do
  !
  ! Subtract the average extension rate from the normal stresses
  !
  do l = 1,nr
    vflx(l,nmb-1+l) = vflx(l,nmb-1+l) - ave_ext_rate
  end do
  !
  ! Multiply the viscous stresses by the molecular viscosity
  !
  vflx(1:nr,nmb:nme) = mvsc * vflx(1:nr,nmb:nme)
  !
  ! Compute the viscous fluxes for the energy equation
  !
  do l = 1,nr
    vflx(l,nq) = mcdt*dpvdx(l,nq) + dot( vel , vflx(l,nmb:nme) )
  end do
  !
  ! Compute the normal viscous flux
  !
  do m = 1,nq
    return_value(m) = dot( vflx(1:nr,m) , unrm(1:nr) )
  end do
  !
end subroutine visc_wall_flx_sub
!
!###############################################################################
!
pure function visc_wall_flx_cv(unrm,uwall,gpwall,dpvdx) result(return_value)
  !
  ! UNRM  : Outward pointing unit normal to the wall
  !         boundary face for current flux point
  ! UWALL : Conservative variables exact solution at the wall
  ! DPVDX : Prescribed wall gradients of the primitive variables at the
  !         current flux point on the wall boundary face. Most gradients
  !         are extrapolated from the interior cell flux point on the
  !         boundary except for certain special cases, e.g., an adiabatic
  !         wall where dT/dx = dT/dy = 0.
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),      intent(in) :: unrm
  real(wp), dimension(1:nq),      intent(in) :: uwall
  real(wp), dimension(1:ng),      intent(in) :: gpwall
  real(wp), dimension(1:nr,1:nq), intent(in) :: dpvdx
  !
  !.. Function Result ..
  real(wp), dimension(1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,l,m
  real(wp) :: mvsc,mcdt
  real(wp) :: ave_ext_rate
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)      :: vel
  real(wp), dimension(1:nr,1:nq) :: vflx
  !
continue
  !
  ! For slip wall: all gradients /= zero
  !                uwall /= zero
  !
  ! For adiabatic wall: Temperature gradient == zero
  !                     Velocity gradients /= zero
  !                     uwall(nmb:nme) [i.e., wall velocities] == zero
  !
  ! For isothermal wall: all gradients /= zero
  !                      uwall(nmb:nme) [i.e., wall velocities] == zero
  !
  ! Initialize the viscous fluxes to zero.
  ! NOTE: This takes care of the continuity equation
  !
  vflx(:,:) = zero
  !
  ! Get the velocity vector
  !
  vel(:) = uwall(nmb:nme)/uwall(nec)
  !
  ! Compute the molecular viscosity and molecular conductivity
  !
  mvsc = gpwall(nvisc)
  mcdt = gpwall(ncond)
  !
  ! Compute the average extension rate
  !
  ave_ext_rate = zero
  do l = 1,nr
    ave_ext_rate = ave_ext_rate + dpvdx(l,nmb-1+l)
  end do
  ave_ext_rate = two3 * ave_ext_rate
  !
  ! Compute the viscous fluxes for the momentum equations which are just
  ! the components of the viscous stress tensor.
  ! NOTE: We are initially excluding the molecular viscosity term which
  !       will be multiplied to the whole tau array after these loops.
  !
  do j = 1,nr
    do i = 1,nr
        vflx(i,nmb-1+j) = dpvdx(i,nmb-1+j) + dpvdx(j,nmb-1+i)
    end do
  end do
  !
  ! Subtract the average extension rate from the normal stresses
  !
  do l = 1,nr
    vflx(l,nmb-1+l) = vflx(l,nmb-1+l) - ave_ext_rate
  end do
  !
  ! Multiply the viscous stresses by the molecular viscosity
  !
  vflx(1:nr,nmb:nme) = mvsc * vflx(1:nr,nmb:nme)
  !
  ! Compute the viscous fluxes for the energy equation
  !
  do l = 1,nr
    vflx(l,nq) = mcdt*dpvdx(l,nq) + dot( vel , vflx(l,nmb:nme) )
  end do
  !
  ! Compute the normal viscous flux
  !
  do m = 1,nq
    return_value(m) = dot( vflx(1:nr,m) , unrm(1:nr) )
  end do
  !
end function visc_wall_flx_cv
!
!###############################################################################
!
!pure &
     subroutine common_vflx_ave_fluxes_gpu(n1,n2,facenrm,face_cv,face_gp, &
                                           face_dcvdx,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: face_dcvdx(:,:,:,:)
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
    call common_vflx_ave_fluxes_gpu_3d(n1,n2,facenrm,face_cv,face_gp, &
                                       face_dcvdx,return_value)
  else if (nr == 2) then
    call common_vflx_ave_fluxes_gpu_2d(n1,n2,facenrm,face_cv,face_gp, &
                                       face_dcvdx,return_value)
  else
    call common_vflx_ave_fluxes_gpu_1d(n1,n2,facenrm,face_cv,face_gp, &
                                       face_dcvdx,return_value)
  end if
  !
end subroutine common_vflx_ave_fluxes_gpu
!
!###############################################################################
!
#ifndef TIME_MINIAPPS_FOR_GPU_COMPARISON
pure &
#endif
subroutine common_vflx_ave_fluxes_gpu_3d(n1,n2,facenrm,face_cv,face_gp, &
                                              face_dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar, only : LDG_tau, LDG_beta
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: face_dcvdx(:,:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,m,np
  real(wp) :: rnx,rny,rnz
  real(wp) :: unx,uny,unz
  real(wp) :: lcv,lfx,lfy,lfz
  real(wp) :: rcv,rfx,rfy,rfz
  real(wp) :: vavex,vavey,vavez
  real(wp) :: ds,vjump,penalty,vflux
  real(wp) :: vflx_x,vflx_y,vflx_z
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: lvflx(:,:,:)
  real(wp), allocatable :: rvflx(:,:,:)
  !
continue
  !
  np = size(face_dcvdx,dim=1)
  !
  allocate ( lvflx(1:np,1:3,1:nq) )
  allocate ( rvflx(1:np,1:3,1:nq) )
  !
  !$ACC DATA CREATE(lvflx,rvflx)
  !
  ! Get the left viscous fluxes
  !
  call visc_flx_cv_gpu_3d(n1,n2,face_cv(:,:,1),face_gp(:,:,1), &
                          face_dcvdx(:,:,:,1),lvflx)
  !
  ! Get the right viscous fluxes
  !
  call visc_flx_cv_gpu_3d(n1,n2,face_cv(:,:,2),face_gp(:,:,2), &
                          face_dcvdx(:,:,:,2),rvflx)
  !
  ! Use the LDG method to compute the common viscous flux vector
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(facenrm,face_cv,return_value,lvflx,rvflx)
  do m = 1,nq
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
      lcv = face_cv(n,m,1)
      rcv = face_cv(n,m,2)
      !
      lfx = lvflx(n,1,m)
      lfy = lvflx(n,2,m)
      lfz = lvflx(n,3,m)
      !
      rfx = rvflx(n,1,m)
      rfy = rvflx(n,2,m)
      rfz = rvflx(n,3,m)
      !
      ! Compute the common averaging terms by applying the
      ! averaging operator {{.}} to the viscous flux vectors
      !
      vavex = half*(lfx + rfx)
      vavey = half*(lfy + rfy)
      vavez = half*(lfz + rfz)
      !
      ! Compute the common jump term by applying the differencing (jump)
      ! operator [[.]] to the interface viscous flux vectors
      !
      vjump = unx*(lfx-rfx) + uny*(lfy-rfy) + unz*(lfz-rfz)
      !
      ! Compute the penalty term (simply the differencing (jump)
      ! operator [[.]] applied to the interface solutions)
      !
      penalty = LDG_tau*(rcv-lcv)
      !
      ! Compute the directional common interface viscous fluxes
      !
      vflx_x = vavex + LDG_beta(1)*vjump + penalty*unx
      vflx_y = vavey + LDG_beta(2)*vjump + penalty*uny
      vflx_z = vavez + LDG_beta(3)*vjump + penalty*unz
      !
      ! Use area face normals to compute the final common
      ! viscous flux normal to the face at this flux point
      !
      vflux = rnx*vflx_x + rny*vflx_y + rnz*vflx_z
      !
      ! Subtract the common viscous flux from the common
      ! inviscid flux already stored within return_value
      !
      return_value(n,m) = return_value(n,m) - vflux
      !
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( rvflx )
  deallocate ( lvflx )
  !
end subroutine common_vflx_ave_fluxes_gpu_3d
!
!###############################################################################
!
!pure &
subroutine common_vflx_ave_fluxes_gpu_2d(n1,n2,facenrm,face_cv,face_gp, &
                                              face_dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar, only : LDG_tau, LDG_beta
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: face_dcvdx(:,:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,m,np
  real(wp) :: rnx,rny
  real(wp) :: unx,uny
  real(wp) :: lcv,lfx,lfy
  real(wp) :: rcv,rfx,rfy
  real(wp) :: vavex,vavey
  real(wp) :: ds,vjump,penalty,vflux
  real(wp) :: vflx_x,vflx_y
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: lvflx(:,:,:)
  real(wp), allocatable :: rvflx(:,:,:)
  !
continue
  !
  np = size(face_dcvdx,dim=1)
  !
  allocate ( lvflx(1:np,1:2,1:nq) )
  allocate ( rvflx(1:np,1:2,1:nq) )
  !
  !$ACC DATA CREATE(lvflx,rvflx)
  !
  ! Get the left viscous fluxes
  !
  call visc_flx_cv_gpu_2d(n1,n2,face_cv(:,:,1),face_gp(:,:,1), &
                          face_dcvdx(:,:,:,1),lvflx)
  !
  ! Get the right viscous fluxes
  !
  call visc_flx_cv_gpu_2d(n1,n2,face_cv(:,:,2),face_gp(:,:,2), &
                          face_dcvdx(:,:,:,2),rvflx)
  !
  ! Use the LDG method to compute the common viscous flux vector
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(facenrm,face_cv,return_value,lvflx,rvflx)
  do m = 1,nq
    do n = n1,n2
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      ds = sqrt(rnx*rnx + rny*rny)
      unx = rnx/ds
      uny = rny/ds
      !
      lcv = face_cv(n,m,1)
      rcv = face_cv(n,m,2)
      !
      lfx = lvflx(n,1,m)
      lfy = lvflx(n,2,m)
      !
      rfx = rvflx(n,1,m)
      rfy = rvflx(n,2,m)
      !
      ! Compute the common averaging terms by applying the
      ! averaging operator {{.}} to the viscous flux vectors
      !
      vavex = half*(lfx + rfx)
      vavey = half*(lfy + rfy)
      !
      ! Compute the common jump term by applying the differencing (jump)
      ! operator [[.]] to the interface viscous flux vectors
      !
      vjump = unx*(lfx-rfx) + uny*(lfy-rfy)
      !
      ! Compute the penalty term (simply the differencing (jump)
      ! operator [[.]] applied to the interface solutions)
      !
      penalty = LDG_tau*(rcv-lcv)
      !
      ! Compute the directional common interface viscous fluxes
      !
      vflx_x = vavex + LDG_beta(1)*vjump + penalty*unx
      vflx_y = vavey + LDG_beta(2)*vjump + penalty*uny
      !
      ! Use area face normals to compute the final common
      ! viscous flux normal to the face at this flux point
      !
      vflux = rnx*vflx_x + rny*vflx_y
      !
      ! Subtract the common viscous flux from the common
      ! inviscid flux already stored within return_value
      !
      return_value(n,m) = return_value(n,m) - vflux
      !
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( rvflx )
  deallocate ( lvflx )
  !
end subroutine common_vflx_ave_fluxes_gpu_2d
!
!###############################################################################
!
#ifndef TIME_MINIAPPS_FOR_GPU_COMPARISON
pure &
#endif
subroutine common_vflx_ave_fluxes_gpu_1d(n1,n2,facenrm,face_cv,face_gp, &
                                              face_dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar, only : LDG_tau, LDG_beta
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: face_dcvdx(:,:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,m,np
  real(wp) :: rnx,lcv,lfx,rcv,rfx,vavex
  real(wp) :: vjump,penalty,vflux,vflx_x
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: lvflx(:,:,:)
  real(wp), allocatable :: rvflx(:,:,:)
  !
continue
  !
  np = size(face_dcvdx,dim=1)
  !
  allocate ( lvflx(1:np,1:1,1:nq) )
  allocate ( rvflx(1:np,1:1,1:nq) )
  !
  !$ACC DATA CREATE(lvflx,rvflx)
  !
  ! Get the left viscous fluxes
  !
  call visc_flx_cv_gpu_1d(n1,n2,face_cv(:,:,1),face_gp(:,:,1), &
                          face_dcvdx(:,:,:,1),lvflx)
  !
  ! Get the right viscous fluxes
  !
  call visc_flx_cv_gpu_1d(n1,n2,face_cv(:,:,2),face_gp(:,:,2), &
                          face_dcvdx(:,:,:,2),rvflx)
  !
  ! Use the LDG method to compute the common viscous flux vector
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(facenrm,face_cv,return_value,lvflx,rvflx)
  do m = 1,nq
    do n = n1,n2
      !
      rnx = facenrm(1,n)
      !
      lcv = face_cv(n,m,1)
      rcv = face_cv(n,m,2)
      !
      lfx = lvflx(n,1,m)
      !
      rfx = rvflx(n,1,m)
      !
      ! Compute the common averaging terms by applying the
      ! averaging operator {{.}} to the viscous flux vectors
      !
      vavex = half*(lfx + rfx)
      !
      ! Compute the common jump term by applying the differencing (jump)
      ! operator [[.]] to the interface viscous flux vectors
      !
      vjump = rnx*(lfx-rfx)
      !
      ! Compute the penalty term (simply the differencing (jump)
      ! operator [[.]] applied to the interface solutions)
      !
      penalty = LDG_tau*(rcv-lcv)
      !
      ! Compute the directional common interface viscous fluxes
      !
      vflx_x = vavex + LDG_beta(1)*vjump + penalty*rnx
      !
      ! Use area face normals to compute the final common
      ! viscous flux normal to the face at this flux point
      !
      vflux = rnx*vflx_x
      !
      ! Subtract the common viscous flux from the common
      ! inviscid flux already stored within return_value
      !
      return_value(n,m) = return_value(n,m) - vflux
      !
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( rvflx )
  deallocate ( lvflx )
  !
end subroutine common_vflx_ave_fluxes_gpu_1d
!
!###############################################################################
!
!pure &
     subroutine common_vflx_ave_sol_and_grad_gpu(n1,n2,facenrm,face_cv, &
                                                 face_dcvdx,return_value)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_dcvdx(:,:,:,:)
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
    call common_vflx_ave_sol_and_grad_gpu_3d(n1,n2,facenrm,face_cv, &
                                             face_dcvdx,return_value)
  else if (nr == 2) then
    call common_vflx_ave_sol_and_grad_gpu_2d(n1,n2,facenrm,face_cv, &
                                             face_dcvdx,return_value)
  else
    call common_vflx_ave_sol_and_grad_gpu_1d(n1,n2,facenrm,face_cv, &
                                             face_dcvdx,return_value)
  end if
  !
end subroutine common_vflx_ave_sol_and_grad_gpu
!
!###############################################################################
!
!pure &
    subroutine common_vflx_ave_sol_and_grad_gpu_3d(n1,n2,facenrm,face_cv, &
                                                    face_dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar,          only : LDG_tau, LDG_beta
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_dcvdx(:,:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,m,np,nr
  real(wp) :: vflux,ds,lcv,rcv,vjump,penalty
  real(wp) :: rnx,unx,ldx,rdx,vavex
  real(wp) :: rny,uny,ldy,rdy,vavey
  real(wp) :: rnz,unz,ldz,rdz,vavez
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: com_cv(:,:)
  real(wp), allocatable :: com_gp(:,:)
  real(wp), allocatable :: com_dcvdx(:,:,:)
  real(wp), allocatable :: com_vflx(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "common_vflx_ave_sol_and_grad_gpu_3d"
  !
continue
  !
  !
  np = size(face_dcvdx,dim=1)
  !
  allocate ( com_cv(1:np,1:nq) )
  allocate ( com_gp(1:np,1:ng) )
  allocate ( com_dcvdx(1:np,1:3,1:nq) )
  allocate ( com_vflx(1:np,1:3,1:nq) )
  !
  !$ACC DATA CREATE(com_cv,com_gp,com_dcvdx,com_vflx)
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(facenrm,face_cv,face_dcvdx,com_cv,com_dcvdx)
  do m = 1,nq
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
      lcv = face_cv(n,m,1)
      rcv = face_cv(n,m,2)
      !
      ldx = face_dcvdx(n,1,m,1)
      ldy = face_dcvdx(n,2,m,1)
      ldz = face_dcvdx(n,3,m,1)
      !
      rdx = face_dcvdx(n,1,m,2)
      rdy = face_dcvdx(n,2,m,2)
      rdz = face_dcvdx(n,3,m,2)
      !
      ! Get the common solution at the interface using a simple average
      !
      com_cv(n,m) = half*(lcv + rcv)
      !
      ! Apply the averaging operator {{.}} to the gradients
      !
      vavex = half*(ldx + rdx)
      vavey = half*(ldy + rdy)
      vavez = half*(ldz + rdz)
      !
      ! Compute the common jump term by applying the
      ! differencing (jump) operator [[.]] to the gradients
      !
      vjump = unx*(ldx-rdx) + uny*(ldy-rdy) + unz*(ldz-rdz)
      !
      ! Compute the penalty term (simply the differencing (jump)
      ! operator [[.]] applied to the interface solutions)
      !
      penalty = LDG_tau*(rcv-lcv)
      !
      ! Finally, compute the directional common interface gradients
      !
      com_dcvdx(n,1,m) = vavex + LDG_beta(1)*vjump + penalty*unx
      com_dcvdx(n,2,m) = vavey + LDG_beta(2)*vjump + penalty*uny
      com_dcvdx(n,3,m) = vavez + LDG_beta(3)*vjump + penalty*unz
      !
    end do
  end do
  !
  ! Update the common thermodynamic and transport properties from
  ! the common solution
  !
  call update_thermo_properties(n1,n2,com_cv,com_gp, &
                                loc="face",caller_name=pname)
  call update_transport_properties(n1,n2,com_cv,com_gp)
  !
  ! Compute the common viscous flux vector using
  ! the common solution and common gradient
  !
  call visc_flx_cv_gpu_3d(n1,n2,com_cv,com_gp,com_dcvdx,com_vflx)
  !
  ! Subtract the common viscous fluxes from the common
  ! inviscid fluxes already stored within return_value
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(facenrm,return_value,com_vflx)
  do m = 1,nq
    do n = n1,n2
      !
      ! Use area face normals to compute the final common
      ! viscous flux normal to the face at this flux point
      !
      vflux = facenrm(1,n)*com_vflx(n,1,m) &
            + facenrm(2,n)*com_vflx(n,2,m) &
            + facenrm(3,n)*com_vflx(n,3,m)
      !
      ! Subtract the common viscous flux from the common
      ! inviscid flux already stored within return_value
      !
      return_value(n,m) = return_value(n,m) - vflux
      !
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( com_vflx )
  deallocate ( com_dcvdx )
  deallocate ( com_gp )
  deallocate ( com_cv )
  !
end subroutine common_vflx_ave_sol_and_grad_gpu_3d
!
!###############################################################################
!
!pure &
subroutine common_vflx_ave_sol_and_grad_gpu_2d(n1,n2,facenrm,face_cv, &
                                                    face_dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar,          only : LDG_tau, LDG_beta
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_dcvdx(:,:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,m,np,nr
  real(wp) :: vflux,ds,lcv,rcv,vjump,penalty
  real(wp) :: rnx,unx,ldx,rdx,vavex
  real(wp) :: rny,uny,ldy,rdy,vavey
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: com_cv(:,:)
  real(wp), allocatable :: com_gp(:,:)
  real(wp), allocatable :: com_dcvdx(:,:,:)
  real(wp), allocatable :: com_vflx(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "common_vflx_ave_sol_and_grad_gpu_2d"
  !
continue
  !
  np = size(face_dcvdx,dim=1)
  !
  allocate ( com_cv(1:np,1:nq) )
  allocate ( com_gp(1:np,1:ng) )
  allocate ( com_dcvdx(1:np,1:2,1:nq) )
  allocate ( com_vflx(1:np,1:2,1:nq) )
  !
  !$ACC DATA CREATE(com_cv,com_gp,com_dcvdx,com_vflx)
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(facenrm,face_cv,face_dcvdx,com_cv,com_dcvdx)
  do m = 1,nq
    do n = n1,n2
      !
      rnx = facenrm(1,n)
      rny = facenrm(2,n)
      ds = sqrt(rnx*rnx + rny*rny)
      unx = rnx/ds
      uny = rny/ds
      !
      lcv = face_cv(n,m,1)
      rcv = face_cv(n,m,2)
      !
      ldx = face_dcvdx(n,1,m,1)
      ldy = face_dcvdx(n,2,m,1)
      !
      rdx = face_dcvdx(n,1,m,2)
      rdy = face_dcvdx(n,2,m,2)
      !
      ! Get the common solution at the interface using a simple average
      !
      com_cv(n,m) = half*(lcv + rcv)
      !
      ! Apply the averaging operator {{.}} to the gradients
      !
      vavex = half*(ldx + rdx)
      vavey = half*(ldy + rdy)
      !
      ! Compute the common jump term by applying the
      ! differencing (jump) operator [[.]] to the gradients
      !
      vjump = unx*(ldx-rdx) + uny*(ldy-rdy)
      !
      ! Compute the penalty term (simply the differencing (jump)
      ! operator [[.]] applied to the interface solutions)
      !
      penalty = LDG_tau*(rcv-lcv)
      !
      ! Finally, compute the directional common interface gradients
      !
      com_dcvdx(n,1,m) = vavex + LDG_beta(1)*vjump + penalty*unx
      com_dcvdx(n,2,m) = vavey + LDG_beta(2)*vjump + penalty*uny
      !
    end do
  end do
  !
  ! Update the common thermodynamic and transport properties from
  ! the common solution
  !
  call update_thermo_properties(n1,n2,com_cv,com_gp, &
                                loc="face",caller_name=pname)
  call update_transport_properties(n1,n2,com_cv,com_gp)
  !
  ! Compute the common viscous flux vector using
  ! the common solution and common gradient
  !
  call visc_flx_cv_gpu_2d(n1,n2,com_cv,com_gp,com_dcvdx,com_vflx)
  !
  ! Subtract the common viscous fluxes from the common
  ! inviscid fluxes already stored within return_value
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(facenrm,return_value,com_vflx)
  do m = 1,nq
    do n = n1,n2
      !
      ! Use area face normals to compute the final common
      ! viscous flux normal to the face at this flux point
      !
      vflux = facenrm(1,n)*com_vflx(n,1,m) + facenrm(2,n)*com_vflx(n,2,m)
      !
      ! Subtract the common viscous flux from the common
      ! inviscid flux already stored within return_value
      !
      return_value(n,m) = return_value(n,m) - vflux
      !
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( com_vflx )
  deallocate ( com_dcvdx )
  deallocate ( com_gp )
  deallocate ( com_cv )
  !
end subroutine common_vflx_ave_sol_and_grad_gpu_2d
!
!###############################################################################
!
!pure &
    subroutine common_vflx_ave_sol_and_grad_gpu_1d(n1,n2,facenrm,face_cv, &
                                                    face_dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar,          only : LDG_tau, LDG_beta
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_dcvdx(:,:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,m,np
  real(wp) :: lcv,rcv,vjump,penalty
  real(wp) :: rnx,ldx,rdx,vavex
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: com_cv(:,:)
  real(wp), allocatable :: com_gp(:,:)
  real(wp), allocatable :: com_dcvdx(:,:,:)
  real(wp), allocatable :: com_vflx(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "common_vflx_ave_sol_and_grad_gpu_1d"
  !
continue
  !
  np = size(face_dcvdx,dim=1)
  !
  allocate ( com_cv(1:np,1:nq) )
  allocate ( com_gp(1:np,1:ng) )
  allocate ( com_dcvdx(1:np,1:1,1:nq) )
  allocate ( com_vflx(1:np,1:1,1:nq) )
  !
  !$ACC DATA CREATE(com_cv,com_gp,com_dcvdx,com_vflx)
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(facenrm,face_cv,face_dcvdx,com_cv,com_dcvdx)
  do m = 1,nq
    do n = n1,n2
      !
      rnx = facenrm(1,n)
      !
      lcv = face_cv(n,m,1)
      rcv = face_cv(n,m,2)
      !
      ldx = face_dcvdx(n,1,m,1)
      rdx = face_dcvdx(n,1,m,2)
      !
      ! Get the common solution at the interface using a simple average
      !
      com_cv(n,m) = half*(lcv + rcv)
      !
      ! Apply the averaging operator {{.}} to the gradients
      !
      vavex = half*(ldx + rdx)
      !
      ! Compute the common jump term by applying the
      ! differencing (jump) operator [[.]] to the gradients
      !
      vjump = rnx*(ldx-rdx)
      !
      ! Compute the penalty term (simply the differencing (jump)
      ! operator [[.]] applied to the interface solutions)
      !
      penalty = LDG_tau*(rcv-lcv)
      !
      ! Finally, compute the directional common interface gradients
      !
      com_dcvdx(n,1,m) = vavex + LDG_beta(1)*vjump + penalty*rnx
      !
    end do
  end do
  !
  ! Update the common thermodynamic and transport properties from
  ! the common solution
  !
  call update_thermo_properties(n1,n2,com_cv,com_gp, &
                                loc="face",caller_name=pname)
  call update_transport_properties(n1,n2,com_cv,com_gp)
  !
  ! Compute the common viscous flux vector using
  ! the common solution and common gradient
  !
  call visc_flx_cv_gpu_1d(n1,n2,com_cv,com_gp,com_dcvdx,com_vflx)
  !
  ! Subtract the common viscous fluxes from the common
  ! inviscid fluxes already stored within return_value
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(facenrm,return_value,com_vflx)
  do m = 1,nq
    do n = n1,n2
      !
      ! Use area face normals to compute the final common
      ! viscous flux normal to the face at this flux point
      ! and subtract it from the common inviscid flux
      ! already stored within return_value
      !
      return_value(n,m) = return_value(n,m) - facenrm(1,n)*com_vflx(n,1,m)
      !
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( com_vflx )
  deallocate ( com_dcvdx )
  deallocate ( com_gp )
  deallocate ( com_cv )
  !
end subroutine common_vflx_ave_sol_and_grad_gpu_1d
!
!###############################################################################
!
!pure &
    subroutine common_visc_flx_sub(unrm,ul,gpl,duldx,ur,gpr,durdx,return_value)
! !$ACC ROUTINE SEQ
  !
  !.. Use Statements ..
  use geovar,        only : nr
  use ovar,          only : visc_variables_to_ave
  use ovar,          only : LDG_tau, LDG_beta
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),      intent(in) :: unrm
  real(wp), dimension(1:nq),      intent(in) :: ul,ur
  real(wp), dimension(1:ng),      intent(in) :: gpl,gpr
  real(wp), dimension(1:nr,1:nq), intent(in) :: duldx,durdx
  real(wp), dimension(1:nq),   intent(inout) :: return_value
  !
  !.. Local Scalars ..
  integer  :: l,m
  real(wp) :: vjump,qjump,penalty
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq)      :: ucom
  real(wp), dimension(1:ng)      :: gpcom
  real(wp), dimension(1:nr,1:nq) :: qcom
  real(wp), dimension(1:nq,1:nr) :: lvflx
  real(wp), dimension(1:nq,1:nr) :: rvflx
  real(wp), dimension(1:nr,1:nq) :: vflux
  real(wp), dimension(1:nq,1:nr) :: fvcom
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "common_visc_flx_sub"
  !
continue
  !
  if (visc_variables_to_ave == visc_interface_fluxes) then
    !
    ! Get the left and right viscous fluxes
    !
    call visc_flx_cv(ul,gpl,duldx,lvflx)
   !lvflx(1:nq,1:nr) = visc_flx_cv(ul,gpl,duldx)
    call visc_flx_cv(ur,gpr,durdx,rvflx)
   !rvflx(1:nq,1:nr) = visc_flx_cv(ur,gpr,durdx)
    !
    ! Use the LDG method to compute the common viscous flux vector
    !
    ! Start by applying the averaging operator {{.}} to the viscous flux vectors
    !
    do m = 1,nq
      do l = 1,nr
        vflux(l,m) = half * (lvflx(m,l) + rvflx(m,l))
      end do
    end do
    !
    ! Apply the differencing (jump) operator [[.]] to the interface viscous
    ! flux vectors and then add it to the common viscous flux vector
    !
    do m = 1,nq
      !
      vjump = dot( unrm(1:nr) , lvflx(m,1:nr)-rvflx(m,1:nr) )
      !
      do l = 1,nr
        vflux(l,m) = vflux(l,m) + LDG_beta(l)*vjump
      end do
      !
    end do
    !
    ! Finally, add the penalty term, (simply the differencing (jump)
    ! operator [[.]] applied to the interface solutions) to the
    ! common viscous flux vector
    !
    do m = 1,nq
      penalty = LDG_tau * ( ur(m) - ul(m) )
      do l = 1,nr
        vflux(l,m) = vflux(l,m) + penalty * unrm(l)
      end do
    end do
    !
    ! Finally, dot the common viscous flux vector with the face unit normal
    !
    do m = 1,nq
      return_value(m) = dot( vflux(1:nr,m) , unrm(1:nr) )
    end do
    !
  else if (visc_variables_to_ave == visc_solution_and_gradients) then
    !
    ! Get the common solution at the interface using a simple average
    !
    do m = 1,nq
      ucom(m) = half * ( ul(m) + ur(m) )
    end do
    !
    ! Use the LDG method to compute the common gradient
    !
    ! Start by applying the averaging operator {{.}} to the gradients
    !
    do m = 1,nq
      do l = 1,nr
        qcom(l,m) = half * (duldx(l,m) + durdx(l,m))
      end do
    end do
    !
    ! Apply the differencing (jump) operator [[.]] to the
    ! gradients and then add it to the common gradient
    !
    do m = 1,nq
      !
      qjump = dot( unrm(1:nr) , duldx(1:nr,m)-durdx(1:nr,m) )
      !
      do l = 1,nr
        qcom(l,m) = qcom(l,m) + LDG_beta(l)*qjump
      end do
      !
    end do
    !
    ! Finally, add the penalty term, (simply the differencing (jump) operator
    ! [[.]] applied to the interface solutions) to the common gradient
    !
    do m = 1,nq
      penalty = LDG_tau * ( ur(m) - ul(m) )
      do l = 1,nr
        qcom(l,m) = qcom(l,m) + penalty * unrm(l)
      end do
    end do
    !
    call update_thermo_properties(ucom,gpcom,caller_name=pname)
    call update_transport_properties(ucom,gpcom)
    !
    ! Compute the common viscous flux vector using
    ! the common solution and common gradient
    !
    call visc_flx_cv(ucom,gpcom,qcom,fvcom)
   !fvcom(1:nq,1:nr) = visc_flx_cv(ucom,gpcom,qcom)
    !
    ! Finally, dot the common viscous flux vector with the face unit normal
    !
    do m = 1,nq
      return_value(m) = dot( fvcom(m,1:nr) , unrm(1:nr) )
    end do
    !
  end if
  !
end subroutine common_visc_flx_sub
!
!###############################################################################
!
!pure &
    function common_visc_flx_cv(unrm,ul,gpl,duldx,ur,gpr,durdx) &
                          result(return_value)
  !
  !.. Use Statements ..
  use geovar,        only : nr
  use ovar,          only : visc_variables_to_ave
  use ovar,          only : LDG_tau, LDG_beta
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),      intent(in) :: unrm
  real(wp), dimension(1:nq),      intent(in) :: ul,ur
  real(wp), dimension(1:ng),      intent(in) :: gpl,gpr
  real(wp), dimension(1:nr,1:nq), intent(in) :: duldx,durdx
  !
  !.. Function Result ..
  real(wp), dimension(1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer  :: l,m
  real(wp) :: vjump,qjump,penalty
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq)      :: ucom
  real(wp), dimension(1:ng)      :: gpcom
  real(wp), dimension(1:nr,1:nq) :: qcom
  real(wp), dimension(1:nq,1:nr) :: lvflx
  real(wp), dimension(1:nq,1:nr) :: rvflx
  real(wp), dimension(1:nr,1:nq) :: vflux
  real(wp), dimension(1:nq,1:nr) :: fvcom
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "common_visc_flx_cv"
  !
continue
  !
  if (visc_variables_to_ave == visc_interface_fluxes) then
    !
    ! Get the left and right viscous fluxes
    !
    call visc_flx_cv(ul,gpl,duldx,lvflx)
   !lvflx(1:nq,1:nr) = visc_flx_cv(ul,gpl,duldx)
    call visc_flx_cv(ur,gpr,durdx,rvflx)
   !rvflx(1:nq,1:nr) = visc_flx_cv(ur,gpr,durdx)
    !
    ! Use the LDG method to compute the common viscous flux vector
    !
    ! Start by applying the averaging operator {{.}} to the viscous flux vectors
    !
    do m = 1,nq
      do l = 1,nr
        vflux(l,m) = half * (lvflx(m,l) + rvflx(m,l))
      end do
    end do
    !
    ! Apply the differencing (jump) operator [[.]] to the interface viscous
    ! flux vectors and then add it to the common viscous flux vector
    !
    do m = 1,nq
      !
      vjump = dot( unrm(1:nr) , lvflx(m,1:nr)-rvflx(m,1:nr) )
      !
      do l = 1,nr
        vflux(l,m) = vflux(l,m) + LDG_beta(l)*vjump
      end do
      !
    end do
    !
    ! Finally, add the penalty term, (simply the differencing (jump)
    ! operator [[.]] applied to the interface solutions) to the
    ! common viscous flux vector
    !
    do m = 1,nq
      penalty = LDG_tau * ( ur(m) - ul(m) )
      do l = 1,nr
        vflux(l,m) = vflux(l,m) + penalty * unrm(l)
      end do
    end do
    !
    ! Finally, dot the common viscous flux vector with the face unit normal
    !
    do m = 1,nq
      return_value(m) = dot( vflux(1:nr,m) , unrm(1:nr) )
    end do
    !
  else if (visc_variables_to_ave == visc_solution_and_gradients) then
    !
    ! Get the common solution at the interface using a simple average
    !
    do m = 1,nq
      ucom(m) = half * ( ul(m) + ur(m) )
    end do
    !
    ! Use the LDG method to compute the common gradient
    !
    ! Start by applying the averaging operator {{.}} to the gradients
    !
    do m = 1,nq
      do l = 1,nr
        qcom(l,m) = half * (duldx(l,m) + durdx(l,m))
      end do
    end do
    !
    ! Apply the differencing (jump) operator [[.]] to the
    ! gradients and then add it to the common gradient
    !
    do m = 1,nq
      !
      qjump = dot( unrm(1:nr) , duldx(1:nr,m)-durdx(1:nr,m) )
      !
      do l = 1,nr
        qcom(l,m) = qcom(l,m) + LDG_beta(l)*qjump
      end do
      !
    end do
    !
    ! Finally, add the penalty term, (simply the differencing (jump) operator
    ! [[.]] applied to the interface solutions) to the common gradient
    !
    do m = 1,nq
      penalty = LDG_tau * ( ur(m) - ul(m) )
      do l = 1,nr
        qcom(l,m) = qcom(l,m) + penalty * unrm(l)
      end do
    end do
    !
    call update_thermo_properties(ucom,gpcom,caller_name=pname)
    call update_transport_properties(ucom,gpcom)
    !
    ! Compute the common viscous flux vector using
    ! the common solution and common gradient
    !
    call visc_flx_cv(ucom,gpcom,qcom,fvcom)
   !fvcom(1:nq,1:nr) = visc_flx_cv(ucom,gpcom,qcom)
    !
    ! Finally, dot the common viscous flux vector with the face unit normal
    !
    do m = 1,nq
      return_value(m) = dot( fvcom(m,1:nr) , unrm(1:nr) )
    end do
    !
  end if
  !
end function common_visc_flx_cv
!
!###############################################################################
!
end module visc_flux_mod
