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
pure subroutine visc_flx_cv_gpu_3d(usp,dusp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd,cpref_nd,Pr
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: dusp(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer :: l,m,n,np
  real(wp) :: cv_inv,cpPr,mvsc,mcdt
  real(wp) :: rho_inv,ave_ext_rate
  real(wp) :: dudx,dudy,dudz
  real(wp) :: dvdx,dvdy,dvdz
  real(wp) :: dwdx,dwdy,dwdz
  real(wp) :: velx,vely,velz
  real(wp) :: dtdx,dtdy,dtdz
  real(wp) :: tauxx,tauyy,tauzz
  real(wp) :: tauxy,tauxz,tauyz
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: tusp(:,:)
  real(wp), allocatable :: tdpvdx(:,:,:)
  !
continue
  !
  np = size(return_value,dim=1)
  !
  cv_inv = (gam - one) / rgasref_nd
  cpPr = cpref_nd/Pr
  !
  allocate ( tusp(1:np,1:nq) )
  allocate ( tdpvdx(1:np,1:3,1:nq) )
  !$ACC DATA CREATE(tusp,tdpvdx) &
  !$ACC&     PRESENT(return_value) &
  !$ACC&     COPYIN(usp,dusp)
  !
  ! Create a transposed copy of usp
  ! NOTE: EMBED 1/RHO IN ALL VARIABLES EXCEPT RHO ITSELF
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = 1,np
    !
    rho_inv = one/usp(nec,n)
    !
    tusp(n,nec) = usp(nec,n)
    tusp(n,nmx) = usp(nmx,n)*rho_inv
    tusp(n,nmy) = usp(nmy,n)*rho_inv
    tusp(n,nmz) = usp(nmz,n)*rho_inv
    tusp(n,nee) = usp(nee,n)*rho_inv
    !
  end do
  !
  ! Temporarily store the transposed conserved variable
  ! gradients in the return_value array
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = 1,np
    !
    return_value(n,1,nec) = dusp(1,nec,n)
    return_value(n,2,nec) = dusp(2,nec,n)
    return_value(n,3,nec) = dusp(3,nec,n)
    !
    return_value(n,1,nmx) = dusp(1,nmx,n)
    return_value(n,2,nmx) = dusp(2,nmx,n)
    return_value(n,3,nmx) = dusp(3,nmx,n)
    !
    return_value(n,1,nmy) = dusp(1,nmy,n)
    return_value(n,2,nmy) = dusp(2,nmy,n)
    return_value(n,3,nmy) = dusp(3,nmy,n)
    !
    return_value(n,1,nmz) = dusp(1,nmz,n)
    return_value(n,2,nmz) = dusp(2,nmz,n)
    return_value(n,3,nmz) = dusp(3,nmz,n)
    !
    return_value(n,1,nee) = dusp(1,nee,n)
    return_value(n,2,nee) = dusp(2,nee,n)
    return_value(n,3,nee) = dusp(3,nee,n)
    !
  end do
  !
  ! Get the gradients of the primitive variables using the 
  ! gradients of the conservative variables and the product rule
  !
  call grad_cv_to_grad_pv_gpu_3d(tusp,return_value,tdpvdx, &
                                 tusp_includes_rho_inv=true)
  ! 
  ! Replace energy variable with molecular viscosity
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = 1,np
    tusp(n,nee) = viscosity_cv_sp( usp(:,n) )
  end do
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = 1,np
    !
    velx = tusp(n,nmx)
    vely = tusp(n,nmy)
    velz = tusp(n,nmz)
    !
    mvsc = tusp(n,nee)
    mcdt = mvsc*cpPr
    !
    dudx = mvsc * tdpvdx(n,1,nmx)
    dudy = mvsc * tdpvdx(n,2,nmx)
    dudz = mvsc * tdpvdx(n,3,nmx)
    !
    dvdx = mvsc * tdpvdx(n,1,nmy)
    dvdy = mvsc * tdpvdx(n,2,nmy)
    dvdz = mvsc * tdpvdx(n,3,nmy)
    !
    dwdx = mvsc * tdpvdx(n,1,nmz)
    dwdy = mvsc * tdpvdx(n,2,nmz)
    dwdz = mvsc * tdpvdx(n,3,nmz)
    !
    dtdx = mcdt * tdpvdx(n,1,nee)
    dtdy = mcdt * tdpvdx(n,2,nee)
    dtdz = mcdt * tdpvdx(n,3,nee)
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
    ! Energy equation
    !
    return_value(n,1,nee) = dtdx + velx*tauxx + vely*tauxy + velz*tauxz
    return_value(n,2,nee) = dtdy + velx*tauxy + vely*tauyy + velz*tauyz
    return_value(n,3,nee) = dtdz + velx*tauxz + vely*tauyz + velz*tauzz
    !
  end do
  !
  !$ACC END DATA
  deallocate ( tusp )
  deallocate ( tdpvdx )
  !
end subroutine visc_flx_cv_gpu_3d
!
!###############################################################################
!
pure subroutine visc_flx_cv_gpu_2d(usp,dusp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd,cpref_nd,Pr
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: dusp(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: l,m,n,np
  real(wp) :: cv_inv,cpPr,mvsc
  real(wp) :: rho_inv,ave_ext_rate
  real(wp) :: dudx,dudy,dudz
  real(wp) :: dvdx,dvdy,dvdz
  real(wp) :: dwdx,dwdy,dwdz
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: tusp(:,:)
  real(wp), allocatable :: tdpvdx(:,:,:)
  !
continue
  !
  np = size(return_value,dim=1)
  !
  cv_inv = (gam - one) / rgasref_nd
  cpPr = cpref_nd/Pr
  !
  allocate ( tusp(1:np,1:nq) )
  allocate ( tdpvdx(1:np,1:2,1:nq) )
  !$ACC DATA CREATE(tusp,tdpvdx) &
  !$ACC&     PRESENT(return_value) &
  !$ACC&     COPYIN(usp,dusp)
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = 1,np
    tusp(n,nec) = usp(nec,n)
  end do
  !
  ! NOTE: EMBED 1/RHO IN ALL OTHER VARIABLES
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do m = nmb,nq
    do n = 1,np
      tusp(n,m) = usp(m,n)/tusp(n,nec)
    end do
  end do
  !
  ! Temporarily store the transposed conserved variable
  ! gradients in the return_value array
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(3)
  do n = 1,np
    do m = 1,nq
      do l = 1,2
        return_value(n,l,m) = dusp(l,m,n)
      end do
    end do
  end do
  !
  ! Get the gradients of the primitive variables using the 
  ! gradients of the conservative variables and the product rule
  !
  call grad_cv_to_grad_pv_gpu_2d(tusp,return_value,tdpvdx, &
                                 tusp_includes_rho_inv=true)
  ! 
  ! Replace energy variable with molecular viscosity
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = 1,np
    tusp(n,nee) = viscosity_cv_sp( usp(:,n) )
  end do
  !
  ! Density flux is zero
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do l = 1,2
    do n = 1,np
      return_value(n,l,nec) = zero
    end do
  end do
  !
  ! Momentum equations
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = 1,np
    !
    mvsc = tusp(n,nee)
    !
    dudx = mvsc * tdpvdx(n,1,nmx)
    dudy = mvsc * tdpvdx(n,2,nmx)
    !
    dvdx = mvsc * tdpvdx(n,1,nmy)
    dvdy = mvsc * tdpvdx(n,2,nmy)
    !
    ave_ext_rate = two3 * (dudx + dvdy)
    !
    return_value(n,1,nmx) = dudx + dudx - ave_ext_rate
    return_value(n,2,nmx) = dvdx + dudy
    !
    return_value(n,1,nmy) = dudy + dvdx
    return_value(n,2,nmy) = dvdy + dvdy - ave_ext_rate
    !
  end do
  !
  ! Energy Equation
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do l = 1,2
    do n = 1,np
      return_value(n,l,nee) = cpPr*tusp(n,nee)*tdpvdx(n,l,nee) &
                            + tusp(n,nmx)*return_value(n,l,nmx) &
                            + tusp(n,nmy)*return_value(n,l,nmy)
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( tusp )
  deallocate ( tdpvdx )
  !
end subroutine visc_flx_cv_gpu_2d
!
!###############################################################################
!
pure subroutine visc_flx_cv_gpu_nd(usp,dusp,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd,cpref_nd,Pr
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: dusp(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer :: l,m,n,np,nr
  real(wp) :: cv_inv,cpPr,mvsc
  real(wp) :: rho_inv,ave_ext_rate
  real(wp) :: dudx,dudy,dudz
  real(wp) :: dvdx,dvdy,dvdz
  real(wp) :: dwdx,dwdy,dwdz
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: tusp(:,:)
  real(wp), allocatable :: tdpvdx(:,:,:)
  !
continue
  !
  np = size(return_value,dim=1)
  nr = size(return_value,dim=2)
  !
  cv_inv = (gam - one) / rgasref_nd
  cpPr = cpref_nd/Pr
  !
  allocate ( tusp(1:np,1:nq) )
  allocate ( tdpvdx(1:np,1:nr,1:nq) )
  !$ACC DATA CREATE(tusp,tdpvdx) &
  !$ACC&     PRESENT(return_value) &
  !$ACC&     COPYIN(usp,dusp)
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = 1,np
    tusp(n,nec) = usp(nec,n)
  end do
  !
  ! NOTE: EMBED 1/RHO IN ALL OTHER VARIABLES
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do m = nmb,nq
    do n = 1,np
      tusp(n,m) = usp(m,n)/tusp(n,nec)
    end do
  end do
  !
  ! Temporarily store the transposed conserved variable
  ! gradients in the return_value array
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(3)
  do n = 1,np
    do m = 1,nq
      do l = 1,nr
        return_value(n,l,m) = dusp(l,m,n)
      end do
    end do
  end do
  !
  ! Get the gradients of the primitive variables using the 
  ! gradients of the conservative variables and the product rule
  !
 !call grad_cv_to_grad_pv_sp_sub(uu,dcvdx,dpvdx)
  call grad_cv_to_grad_pv_gpu_nd(tusp,return_value,tdpvdx, &
                                 tusp_includes_rho_inv=true)
  ! 
  ! Replace energy variable with molecular viscosity
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = 1,np
    tusp(n,nee) = viscosity_cv_sp( usp(:,n) )
  end do
  !
  ! Density flux is zero
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do l = 1,nr
    do n = 1,np
      return_value(n,l,nec) = zero
    end do
  end do
  !
  ! Momentum equations
  !
  if (nr == 3) then
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
    do n = 1,np
      !
      mvsc = tusp(n,nee)
      !
      dudx = mvsc * tdpvdx(n,1,nmx)
      dudy = mvsc * tdpvdx(n,2,nmx)
      dudz = mvsc * tdpvdx(n,3,nmx)
      !
      dvdx = mvsc * tdpvdx(n,1,nmy)
      dvdy = mvsc * tdpvdx(n,2,nmy)
      dvdz = mvsc * tdpvdx(n,3,nmy)
      !
      dwdx = mvsc * tdpvdx(n,1,nmz)
      dwdy = mvsc * tdpvdx(n,2,nmz)
      dwdz = mvsc * tdpvdx(n,3,nmz)
      !
      ave_ext_rate = two3 * (dudx + dvdy + dwdz)
      !
      return_value(n,1,nmx) = dudx + dudx - ave_ext_rate
      return_value(n,2,nmx) = dvdx + dudy
      return_value(n,3,nmx) = dwdx + dudz
      !
      return_value(n,1,nmy) = dudy + dvdx
      return_value(n,2,nmy) = dvdy + dvdy - ave_ext_rate
      return_value(n,3,nmy) = dwdy + dvdz
      !
      return_value(n,1,nmz) = dudz + dwdx
      return_value(n,2,nmz) = dvdz + dwdy
      return_value(n,3,nmz) = dwdz + dwdz - ave_ext_rate
      !
    end do
    !
  else if (nr == 2) then
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
    do n = 1,np
      !
      mvsc = tusp(n,nee)
      !
      dudx = mvsc * tdpvdx(n,1,nmx)
      dudy = mvsc * tdpvdx(n,2,nmx)
      !
      dvdx = mvsc * tdpvdx(n,1,nmy)
      dvdy = mvsc * tdpvdx(n,2,nmy)
      !
      ave_ext_rate = two3 * (dudx + dvdy)
      !
      return_value(n,1,nmx) = dudx + dudx - ave_ext_rate
      return_value(n,2,nmx) = dvdx + dudy
      !
      return_value(n,1,nmy) = dudy + dvdx
      return_value(n,2,nmy) = dvdy + dvdy - ave_ext_rate
      !
    end do
    !
  else
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
    do n = 1,np
      return_value(n,1,nmx) = two * two3 * tusp(n,nee) * tdpvdx(n,1,nmx)
    end do
    !
  end if
  !
  ! Energy Equation
  !
  if (nr == 3) then
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do l = 1,nr
      do n = 1,np
        return_value(n,l,nee) = cpPr*tusp(n,nee)*tdpvdx(n,l,nee) &
                              + tusp(n,nmx)*return_value(n,l,nmx) &
                              + tusp(n,nmy)*return_value(n,l,nmy) &
                              + tusp(n,nmz)*return_value(n,l,nmz)
      end do
    end do
    !
  else if (nr == 2) then
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do l = 1,nr
      do n = 1,np
        return_value(n,l,nee) = cpPr*tusp(n,nee)*tdpvdx(n,l,nee) &
                              + tusp(n,nmx)*return_value(n,l,nmx) &
                              + tusp(n,nmy)*return_value(n,l,nmy)
      end do
    end do
    !
  else
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do l = 1,nr
      do n = 1,np
        return_value(n,l,nee) = cpPr*tusp(n,nee)*tdpvdx(n,l,nee) &
                              + tusp(n,nmx)*return_value(n,l,nmx)
      end do
    end do
    !
  end if
  !
  !$ACC END DATA
  deallocate ( tusp )
  deallocate ( tdpvdx )
  !
end subroutine visc_flx_cv_gpu_nd
!
!###############################################################################
!
pure subroutine grad_cv_to_grad_pv_gpu_3d(tusp,tdcvdx,return_value, &
                                          tusp_includes_rho_inv)
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: tusp(:,:)
  real(wp),    intent(in) :: tdcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: tusp_includes_rho_inv
  !
  !.. Local Scalars ..
  integer  :: n,np
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
  np = size(tdcvdx,dim=1)
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  ! Embed 1/rho in the all variables except density
  !
  if (tusp_needs_rho_inv) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
    do n = 1,np
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = 1,np
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = 1,np
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = 1,np
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = 1,np
    return_value(n,1,nee) = return_value(n,1,nee) * cv_inv
    return_value(n,2,nee) = return_value(n,2,nee) * cv_inv
    return_value(n,3,nee) = return_value(n,3,nee) * cv_inv
  end do
  !
  ! Remove the embedded 1/rho in the all variables except density
  !
  if (tusp_needs_rho_inv) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
    do n = 1,np
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
pure subroutine grad_cv_to_grad_pv_gpu_2d(tusp,tdcvdx,return_value, &
                                          tusp_includes_rho_inv)
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: tusp(:,:)
  real(wp),    intent(in) :: tdcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: tusp_includes_rho_inv
  !
  !.. Local Scalars ..
  integer  :: l,m,n,np
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
  np = size(tdcvdx,dim=1)
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  ! Embed 1/rho in the all variables except density
  !
  if (tusp_needs_rho_inv) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do m = nmb,nq
      do n = 1,np
        tusp(n,m) = tusp(n,m)/tusp(n,nec)
      end do
    end do
  end if
  !
  ! Copy the density gradient
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do l = 1,2
    do n = 1,np
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(3)
  do m = nmb,nq
    do l = 1,2
      do n = 1,np
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
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do l = 1,2
      do n = 1,np
        return_value(n,l,nee) = return_value(n,l,nee) &
                              - return_value(n,l,m)*tusp(n,m)
      end do
    end do
  end do
#else
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do l = 1,2
    do n = 1,np
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do l = 1,2
    do n = 1,np
      return_value(n,l,nee) = return_value(n,l,nee) * cv_inv
    end do
  end do
  !
  ! Remove the embedded 1/rho in the all variables except density
  !
  if (tusp_needs_rho_inv) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do m = nmb,nq
      do n = 1,np
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
pure subroutine grad_cv_to_grad_pv_gpu_nd(tusp,tdcvdx,return_value, &
                                          tusp_includes_rho_inv)
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: tusp(:,:)
  real(wp),    intent(in) :: tdcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: tusp_includes_rho_inv
  !
  !.. Local Scalars ..
  integer  :: l,m,n,np,nr
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
  np = size(tdcvdx,dim=1)
  nr = size(tdcvdx,dim=2)
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  ! Embed 1/rho in the all variables except density
  !
  if (tusp_needs_rho_inv) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do m = nmb,nq
      do n = 1,np
        tusp(n,m) = tusp(n,m)/tusp(n,nec)
      end do
    end do
  end if
  !
  ! Copy the density gradient
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do l = 1,nr
    do n = 1,np
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(3)
  do m = nmb,nq
    do l = 1,nr
      do n = 1,np
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
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do l = 1,nr
      do n = 1,np
        return_value(n,l,nee) = return_value(n,l,nee) &
                              - return_value(n,l,m)*tusp(n,m)
      end do
    end do
  end do
#else
  if (nr == 3) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do l = 1,3
      do n = 1,np
        return_value(n,l,nee) = return_value(n,l,nee) &
                              - return_value(n,l,nmx)*tusp(n,nmx) &
                              - return_value(n,l,nmy)*tusp(n,nmy) &
                              - return_value(n,l,nmz)*tusp(n,nmz)
      end do
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do l = 1,2
      do n = 1,np
        return_value(n,l,nee) = return_value(n,l,nee) &
                              - return_value(n,l,nmx)*tusp(n,nmx) &
                              - return_value(n,l,nmy)*tusp(n,nmy)
      end do
    end do
  else
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
    do n = 1,np
      return_value(n,1,nee) = return_value(n,1,nee) &
                            - return_value(n,1,nmx)*tusp(n,nmx)
    end do
  end if
#endif
  !
  ! Finish computing the temperature gradient by dividing through by
  ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do l = 1,nr
    do n = 1,np
      return_value(n,l,nee) = return_value(n,l,nee) * cv_inv
    end do
  end do
  !
  ! Remove the embedded 1/rho in the all variables except density
  !
  if (tusp_needs_rho_inv) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do m = nmb,nq
      do n = 1,np
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
pure subroutine visc_flx_cv_pts_gpu_3d(n1,n2,tusp,tdcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd,cpref_nd,Pr
  use gpu_mod, only : viscosity_cv_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: tusp(:,:)
  real(wp),    intent(in) :: tdcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: l,m,n,np
  real(wp) :: cv_inv,cpPr,mvsc,mcdt
  real(wp) :: rho_inv,ave_ext_rate
  real(wp) :: dudx,dudy,dudz
  real(wp) :: dvdx,dvdy,dvdz
  real(wp) :: dwdx,dwdy,dwdz
  real(wp) :: velx,vely,velz
  real(wp) :: tauxx,tauyy,tauzz
  real(wp) :: tauxy,tauxz,tauyz
  real(wp) :: dtdx,dtdy,dtdz
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: visc(:)
  real(wp), allocatable :: tdpvdx(:,:,:)
  !
continue
  !
  np = size(return_value,dim=1)
  !
  cv_inv = (gam - one) / rgasref_nd
  cpPr = cpref_nd/Pr
  !
  allocate ( visc(1:np) )
  allocate ( tdpvdx(1:np,1:3,1:nq) )
  !$ACC DATA CREATE(visc,tdpvdx) &
  !$ACC&     PRESENT(tusp,tdcvdx,return_value)
  ! 
  ! Compute the molecular viscosity
  !
  call viscosity_cv_gpu(n1,n2,tusp,visc)
  !
  ! NOTE: EMBED 1/RHO IN ALL OTHER VARIABLES
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = n1,n2
    rho_inv = one/tusp(n,nec)
    tusp(n,nmx) = tusp(n,nmx)*rho_inv
    tusp(n,nmy) = tusp(n,nmy)*rho_inv
    tusp(n,nmz) = tusp(n,nmz)*rho_inv
    tusp(n,nee) = tusp(n,nee)*rho_inv
  end do
  !
  ! Get the gradients of the primitive variables using the 
  ! gradients of the conservative variables and the product rule
  !
  call grad_cv_to_grad_pv_pts_gpu_3d(n1,n2,tusp,tdcvdx,tdpvdx, &
                                     tusp_includes_rho_inv=true)
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = n1,n2
    !
    velx = tusp(n,nmx)
    vely = tusp(n,nmy)
    velz = tusp(n,nmz)
    !
    mvsc = visc(n)
    mcdt = mvsc*cpPr
    !
    dudx = mvsc * tdpvdx(n,1,nmx)
    dudy = mvsc * tdpvdx(n,2,nmx)
    dudz = mvsc * tdpvdx(n,3,nmx)
    !
    dvdx = mvsc * tdpvdx(n,1,nmy)
    dvdy = mvsc * tdpvdx(n,2,nmy)
    dvdz = mvsc * tdpvdx(n,3,nmy)
    !
    dwdx = mvsc * tdpvdx(n,1,nmz)
    dwdy = mvsc * tdpvdx(n,2,nmz)
    dwdz = mvsc * tdpvdx(n,3,nmz)
    !
    dtdx = mcdt * tdpvdx(n,1,nee)
    dtdy = mcdt * tdpvdx(n,2,nee)
    dtdz = mcdt * tdpvdx(n,3,nee)
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
    return_value(n,1,nee) = dtdx + velx*tauxx + vely*tauxy + velz*tauxz
    return_value(n,2,nee) = dtdy + velx*tauxy + vely*tauyy + velz*tauyz
    return_value(n,3,nee) = dtdz + velx*tauxz + vely*tauyz + velz*tauzz
    !
  end do
  !
  !$ACC END DATA
  deallocate ( tdpvdx )
  deallocate ( visc )
  !
end subroutine visc_flx_cv_pts_gpu_3d
!
!###############################################################################
!
pure subroutine visc_flx_cv_pts_gpu_2d(n1,n2,tusp,tdcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd,cpref_nd,Pr
  use gpu_mod, only : viscosity_cv_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: tusp(:,:)
  real(wp),    intent(in) :: tdcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: l,m,n,np
  real(wp) :: cv_inv,cpPr,mvsc
  real(wp) :: rho_inv,ave_ext_rate
  real(wp) :: dudx,dudy
  real(wp) :: dvdx,dvdy
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: visc(:)
  real(wp), allocatable :: tdpvdx(:,:,:)
  !
continue
  !
  np = size(return_value,dim=1)
  !
  cv_inv = (gam - one) / rgasref_nd
  cpPr = cpref_nd/Pr
  !
  allocate ( visc(1:np) )
  allocate ( tdpvdx(1:np,1:2,1:nq) )
  !$ACC DATA CREATE(visc,tdpvdx) &
  !$ACC&     PRESENT(tusp,tdcvdx,return_value)
  ! 
  ! Compute the molecular viscosity
  !
  call viscosity_cv_gpu(n1,n2,tusp,visc)
  !
  ! NOTE: EMBED 1/RHO IN ALL OTHER VARIABLES
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do m = nmb,nq
    do n = n1,n2
      tusp(n,m) = tusp(n,m)/tusp(n,nec)
    end do
  end do
  !
  ! Get the gradients of the primitive variables using the 
  ! gradients of the conservative variables and the product rule
  !
  call grad_cv_to_grad_pv_pts_gpu_2d(n1,n2,tusp,tdcvdx,tdpvdx, &
                                     tusp_includes_rho_inv=true)
  !
  ! Density flux is zero
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do l = 1,2
    do n = n1,n2
      return_value(n,l,nec) = zero
    end do
  end do
  !
  ! Momentum equations
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = n1,n2
    !
    mvsc = visc(n)
    !
    dudx = mvsc * tdpvdx(n,1,nmx)
    dudy = mvsc * tdpvdx(n,2,nmx)
    !
    dvdx = mvsc * tdpvdx(n,1,nmy)
    dvdy = mvsc * tdpvdx(n,2,nmy)
    !
    ave_ext_rate = two3 * (dudx + dvdy)
    !
    return_value(n,1,nmx) = dudx + dudx - ave_ext_rate
    return_value(n,2,nmx) = dvdx + dudy
    !
    return_value(n,1,nmy) = dudy + dvdx
    return_value(n,2,nmy) = dvdy + dvdy - ave_ext_rate
    !
  end do
  !
  ! Energy Equation
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do l = 1,2
    do n = n1,n2
      return_value(n,l,nee) = cpPr*visc(n)*tdpvdx(n,l,nee) &
                            + tusp(n,nmx)*return_value(n,l,nmx) &
                            + tusp(n,nmy)*return_value(n,l,nmy)
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( tdpvdx )
  deallocate ( visc )
  !
end subroutine visc_flx_cv_pts_gpu_2d
!
!###############################################################################
!
pure subroutine visc_flx_cv_pts_gpu_nd(n1,n2,tusp,tdcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd,cpref_nd,Pr
  use gpu_mod, only : viscosity_cv_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: tusp(:,:)
  real(wp),    intent(in) :: tdcvdx(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: l,m,n,np,nr
  real(wp) :: cv_inv,cpPr,mvsc
  real(wp) :: rho_inv,ave_ext_rate
  real(wp) :: dudx,dudy,dudz
  real(wp) :: dvdx,dvdy,dvdz
  real(wp) :: dwdx,dwdy,dwdz
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: visc(:)
  real(wp), allocatable :: tdpvdx(:,:,:)
  !
continue
  !
  np = size(return_value,dim=1)
  nr = size(return_value,dim=2)
  !
  cv_inv = (gam - one) / rgasref_nd
  cpPr = cpref_nd/Pr
  !
  allocate ( visc(1:np) )
  allocate ( tdpvdx(1:np,1:nr,1:nq) )
  !$ACC DATA CREATE(visc,tdpvdx) &
  !$ACC&     PRESENT(tusp,tdcvdx,return_value)
  ! 
  ! Compute the molecular viscosity
  !
  call viscosity_cv_gpu(n1,n2,tusp,visc)
  !
  ! NOTE: EMBED 1/RHO IN ALL OTHER VARIABLES
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do m = nmb,nq
    do n = n1,n2
      tusp(n,m) = tusp(n,m)/tusp(n,nec)
    end do
  end do
  !
  ! Get the gradients of the primitive variables using the 
  ! gradients of the conservative variables and the product rule
  !
 !call grad_cv_to_grad_pv_sp_sub(uu,dcvdx,dpvdx)
  call grad_cv_to_grad_pv_pts_gpu_nd(n1,n2,tusp,tdcvdx,tdpvdx, &
                                     tusp_includes_rho_inv=true)
  !
  ! Density flux is zero
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do l = 1,nr
    do n = n1,n2
      return_value(n,l,nec) = zero
    end do
  end do
  !
  ! Momentum equations
  !
  if (nr == 3) then
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
    do n = n1,n2
      !
      mvsc = visc(n)
      !
      dudx = mvsc * tdpvdx(n,1,nmx)
      dudy = mvsc * tdpvdx(n,2,nmx)
      dudz = mvsc * tdpvdx(n,3,nmx)
      !
      dvdx = mvsc * tdpvdx(n,1,nmy)
      dvdy = mvsc * tdpvdx(n,2,nmy)
      dvdz = mvsc * tdpvdx(n,3,nmy)
      !
      dwdx = mvsc * tdpvdx(n,1,nmz)
      dwdy = mvsc * tdpvdx(n,2,nmz)
      dwdz = mvsc * tdpvdx(n,3,nmz)
      !
      ave_ext_rate = two3 * (dudx + dvdy + dwdz)
      !
      return_value(n,1,nmx) = dudx + dudx - ave_ext_rate
      return_value(n,2,nmx) = dvdx + dudy
      return_value(n,3,nmx) = dwdx + dudz
      !
      return_value(n,1,nmy) = dudy + dvdx
      return_value(n,2,nmy) = dvdy + dvdy - ave_ext_rate
      return_value(n,3,nmy) = dwdy + dvdz
      !
      return_value(n,1,nmz) = dudz + dwdx
      return_value(n,2,nmz) = dvdz + dwdy
      return_value(n,3,nmz) = dwdz + dwdz - ave_ext_rate
      !
    end do
    !
  else if (nr == 2) then
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
    do n = n1,n2
      !
      mvsc = visc(n)
      !
      dudx = mvsc * tdpvdx(n,1,nmx)
      dudy = mvsc * tdpvdx(n,2,nmx)
      !
      dvdx = mvsc * tdpvdx(n,1,nmy)
      dvdy = mvsc * tdpvdx(n,2,nmy)
      !
      ave_ext_rate = two3 * (dudx + dvdy)
      !
      return_value(n,1,nmx) = dudx + dudx - ave_ext_rate
      return_value(n,2,nmx) = dvdx + dudy
      !
      return_value(n,1,nmy) = dudy + dvdx
      return_value(n,2,nmy) = dvdy + dvdy - ave_ext_rate
      !
    end do
    !
  else
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
    do n = n1,n2
      return_value(n,1,nmx) = two * two3 * visc(n) * tdpvdx(n,1,nmx)
    end do
    !
  end if
  !
  ! Energy Equation
  !
  if (nr == 3) then
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do l = 1,3
      do n = n1,n2
        return_value(n,l,nee) = cpPr*visc(n)*tdpvdx(n,l,nee) &
                              + tusp(n,nmx)*return_value(n,l,nmx) &
                              + tusp(n,nmy)*return_value(n,l,nmy) &
                              + tusp(n,nmz)*return_value(n,l,nmz)
      end do
    end do
    !
  else if (nr == 2) then
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do l = 1,2
      do n = n1,n2
        return_value(n,l,nee) = cpPr*visc(n)*tdpvdx(n,l,nee) &
                              + tusp(n,nmx)*return_value(n,l,nmx) &
                              + tusp(n,nmy)*return_value(n,l,nmy)
      end do
    end do
    !
  else
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do l = 1,nr
      do n = n1,n2
        return_value(n,l,nee) = cpPr*visc(n)*tdpvdx(n,l,nee) &
                              + tusp(n,nmx)*return_value(n,l,nmx)
      end do
    end do
    !
  end if
  !
  !$ACC END DATA
  deallocate ( tdpvdx )
  deallocate ( visc )
  !
end subroutine visc_flx_cv_pts_gpu_nd
!
!###############################################################################
!
pure subroutine grad_cv_to_grad_pv_pts_gpu_3d(n1,n2,tusp,tdcvdx,return_value, &
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
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
  do n = n1,n2
    return_value(n,1,nee) = return_value(n,1,nee) * cv_inv
    return_value(n,2,nee) = return_value(n,2,nee) * cv_inv
    return_value(n,3,nee) = return_value(n,3,nee) * cv_inv
  end do
  !
  ! Remove the embedded 1/rho in the all variables except density
  !
  if (tusp_needs_rho_inv) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
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
end subroutine grad_cv_to_grad_pv_pts_gpu_3d
!
!###############################################################################
!
pure subroutine grad_cv_to_grad_pv_pts_gpu_2d(n1,n2,tusp,tdcvdx,return_value, &
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
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do m = nmb,nq
      do n = n1,n2
        tusp(n,m) = tusp(n,m)/tusp(n,nec)
      end do
    end do
  end if
  !
  ! Copy the density gradient
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(3)
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
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do l = 1,2
      do n = n1,n2
        return_value(n,l,nee) = return_value(n,l,nee) &
                              - return_value(n,l,m)*tusp(n,m)
      end do
    end do
  end do
#else
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do l = 1,2
    do n = n1,n2
      return_value(n,l,nee) = return_value(n,l,nee) * cv_inv
    end do
  end do
  !
  ! Remove the embedded 1/rho in the all variables except density
  !
  if (tusp_needs_rho_inv) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do m = nmb,nq
      do n = n1,n2
        tusp(n,m) = tusp(n,m)*tusp(n,nec)
      end do
    end do
  end if
  !$ACC END DATA
  !
end subroutine grad_cv_to_grad_pv_pts_gpu_2d
!
!###############################################################################
!
pure subroutine grad_cv_to_grad_pv_pts_gpu_nd(n1,n2,tusp,tdcvdx,return_value, &
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
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do m = nmb,nq
      do n = n1,n2
        tusp(n,m) = tusp(n,m)/tusp(n,nec)
      end do
    end do
  end if
  !
  ! Copy the density gradient
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(3)
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
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do l = 1,nr
      do n = n1,n2
        return_value(n,l,nee) = return_value(n,l,nee) &
                              - return_value(n,l,m)*tusp(n,m)
      end do
    end do
  end do
#else
  if (nr == 3) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do l = 1,3
      do n = n1,n2
        return_value(n,l,nee) = return_value(n,l,nee) &
                              - return_value(n,l,nmx)*tusp(n,nmx) &
                              - return_value(n,l,nmy)*tusp(n,nmy) &
                              - return_value(n,l,nmz)*tusp(n,nmz)
      end do
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do l = 1,2
      do n = n1,n2
        return_value(n,l,nee) = return_value(n,l,nee) &
                              - return_value(n,l,nmx)*tusp(n,nmx) &
                              - return_value(n,l,nmy)*tusp(n,nmy)
      end do
    end do
  else
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
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
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do l = 1,nr
    do n = n1,n2
      return_value(n,l,nee) = return_value(n,l,nee) * cv_inv
    end do
  end do
  !
  ! Remove the embedded 1/rho in the all variables except density
  !
  if (tusp_needs_rho_inv) then
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do m = nmb,nq
      do n = n1,n2
        tusp(n,m) = tusp(n,m)*tusp(n,nec)
      end do
    end do
  end if
  !$ACC END DATA
  !
end subroutine grad_cv_to_grad_pv_pts_gpu_nd
!
!###############################################################################
!
pure function visc_flx_cv1(uu,dcvdx) result(return_value) 
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : cpref_nd,Pr
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nq),      intent(in) :: uu
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
  dpvdx(1:nr,1:nq) = grad_cv_to_grad_pv_sp(uu,dcvdx)
  !
  ! Get the velocity vector
  !
  vel(:) = uu(nmb:nme)/uu(nec)
  !
  ! Compute the molecular viscosity and molecular conductivity
  !
  mvsc = viscosity_cv_sp(uu)
  mcdt = cpref_nd * mvsc / Pr
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
!pure function visc_flx_cv(uu,dcvdx) result(return_value) 
pure subroutine visc_flx_cv(uu,dcvdx,return_value)
  !$ACC ROUTINE SEQ
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : cpref_nd,Pr
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nq),      intent(in) :: uu
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
 !dpvdx(1:nr,1:nq) = grad_cv_to_grad_pv_sp_sub(uu,dcvdx)
  call grad_cv_to_grad_pv_sp_sub(uu,dcvdx,dpvdx)
  !
  ! Get the velocity vector
  !
  vel(:) = uu(nmb:nme)/uu(nec)
  !
  ! Compute the molecular viscosity and molecular conductivity
  !
  mvsc = viscosity_cv_sp(uu)
  mcdt = cpref_nd * mvsc / Pr
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
pure function normal_visc_flx_cv(rnrm,uu,dcvdx) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nq),      intent(in) :: uu
  real(wp), dimension(1:nr),         intent(in) :: rnrm
  real(wp), dimension(1:nr,1:nq), intent(in) :: dcvdx
  !
  !.. Function Result ..
  real(wp), dimension(1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer :: l
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr)         :: unrm
  real(wp), dimension(1:nq,1:nr) :: vflx
  !
continue
  !
  unrm(:) = unit_vector( rnrm(:) )
  !
 !vflx(1:nq,1:nr) = visc_flx_cv(uu,dcvdx)
  call visc_flx_cv(uu,dcvdx,vflx)
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
pure subroutine visc_wall_flx_gpu_3d(n1,n2,facenrm,face_cv, &
                                     face_dcvdx,return_value)
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
  use ovar, only : cpref_nd,Pr
  use gpu_mod, only : viscosity_cv_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_dcvdx(:,:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,np
  real(wp) :: rho_inv,mvsc,mcdt,ave_ext_rate,tflux,cpPr
  real(wp) :: tauxx,tauyy,tauzz,tauxy,tauxz,tauyz
  real(wp) :: rnx,velx,dudx,dvdx,dwdx,dtdx,tflx_x,uflux
  real(wp) :: rny,vely,dudy,dvdy,dwdy,dtdy,tflx_y,vflux
  real(wp) :: rnz,velz,dudz,dvdz,dwdz,dtdz,tflx_z,wflux
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: visc(:)
  !
continue
  !
  np = size(face_dcvdx,dim=1)
  !
  cpPr = cpref_nd/Pr
  !
  allocate ( visc(1:np) )
  !
  !$ACC DATA CREATE(visc) &
  !$ACC&     PRESENT(facenrm,face_cv,face_dcvdx,return_value)
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
  call viscosity_cv_gpu(n1,n2,face_cv(:,:,2),visc)
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
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
    mvsc = visc(n)
    mcdt = mvsc*cpPr
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
  !$ACC END DATA
  deallocate ( visc )
  !
end subroutine visc_wall_flx_gpu_3d
!
!###############################################################################
!
pure subroutine visc_wall_flx_gpu_2d(n1,n2,facenrm,face_cv, &
                                     face_dcvdx,return_value)
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
  use ovar, only : cpref_nd,Pr
  use gpu_mod, only : viscosity_cv_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_dcvdx(:,:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,np
  real(wp) :: rho_inv,mvsc,mcdt,ave_ext_rate,tflux,cpPr
  real(wp) :: tauxx,tauyy,tauxy
  real(wp) :: rnx,velx,dudx,dvdx,dwdx,dtdx,tflx_x,uflux
  real(wp) :: rny,vely,dudy,dvdy,dwdy,dtdy,tflx_y,vflux
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: visc(:)
  !
continue
  !
  np = size(face_dcvdx,dim=1)
  !
  cpPr = cpref_nd/Pr
  !
  allocate ( visc(1:np) )
  !
  !$ACC DATA CREATE(visc) &
  !$ACC&     PRESENT(facenrm,face_cv,face_dcvdx,return_value)
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
  call viscosity_cv_gpu(n1,n2,face_cv(:,:,2),visc)
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
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
    mvsc = visc(n)
    mcdt = mvsc*cpPr
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
  !$ACC END DATA
  deallocate ( visc )
  !
end subroutine visc_wall_flx_gpu_2d
!
!###############################################################################
!
pure subroutine visc_wall_flx_gpu_nd(n1,n2,facenrm,face_cv, &
                                     face_dcvdx,return_value)
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
  use ovar, only : cpref_nd,Pr
  use gpu_mod, only : viscosity_cv_gpu
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_dcvdx(:,:,:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,np,nr
  real(wp) :: rho_inv,mvsc,mcdt,ave_ext_rate,tflux,cpPr
  real(wp) :: tauxx,tauyy,tauzz,tauxy,tauxz,tauyz
  real(wp) :: rnx,velx,dudx,dvdx,dwdx,dtdx,tflx_x,uflux
  real(wp) :: rny,vely,dudy,dvdy,dwdy,dtdy,tflx_y,vflux
  real(wp) :: rnz,velz,dudz,dvdz,dwdz,dtdz,tflx_z,wflux
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: visc(:)
  !
continue
  !
  np = size(face_dcvdx,dim=1)
  nr = size(face_dcvdx,dim=2)
  !
  cpPr = cpref_nd/Pr
  !
  allocate ( visc(1:np) )
  !
  !$ACC DATA CREATE(visc) &
  !$ACC&     PRESENT(facenrm,face_cv,face_dcvdx,return_value)
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
  call viscosity_cv_gpu(n1,n2,face_cv(:,:,2),visc)
  !
  if (nr == 3) then
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
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
      mvsc = visc(n)
      mcdt = mvsc*cpPr
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
  else if (nr == 2) then
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
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
      mvsc = visc(n)
      mcdt = mvsc*cpPr
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
  else
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32)
    do n = n1,n2
      !
      rnx = facenrm(1,n)
      !
      velx = face_cv(n,nmx,2)/face_cv(n,nec,2)
      !
      mvsc = visc(n)
      mcdt = mvsc*cpPr
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
  end if
  !
  !$ACC END DATA
  deallocate ( visc )
  !
end subroutine visc_wall_flx_gpu_nd
!
!###############################################################################
!
pure subroutine visc_wall_flx_sub(unrm,uwall,dpvdx,return_value) 
  !$ACC ROUTINE SEQ
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
  use ovar,   only : cpref_nd,Pr
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),      intent(in) :: unrm
  real(wp), dimension(1:nq),      intent(in) :: uwall
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
  mvsc = viscosity_cv_sp(uwall)
  mcdt = cpref_nd * mvsc / Pr
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
pure function visc_wall_flx_cv(unrm,uwall,dpvdx) result(return_value) 
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
  use ovar,   only : cpref_nd,Pr
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),      intent(in) :: unrm
  real(wp), dimension(1:nq),      intent(in) :: uwall
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
  mvsc = viscosity_cv_sp(uwall)
  mcdt = cpref_nd * mvsc / Pr
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
pure subroutine common_vflx_ave_fluxes_gpu_3d(n1,n2,facenrm,face_cv, &
                                              face_dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar, only : LDG_tau, LDG_beta
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
  real(wp) :: rnx,rny,rnz
  real(wp) :: unx,uny,unz
  real(wp) :: lcv,lfx,lfy,lfz
  real(wp) :: rcv,rfx,rfy,rfz
  real(wp) :: vavex,vavey,vavez
  real(wp) :: ds,vjump,penalty,vflux
  real(wp) :: vflx_x,vflx_y,vflx_z
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: tusp(:,:)
  real(wp), allocatable :: lvflx(:,:,:)
  real(wp), allocatable :: rvflx(:,:,:)
  !
continue
  !
  np = size(face_dcvdx,dim=1)
  !
  allocate ( tusp(1:np,1:nq) )
  allocate ( lvflx(1:np,1:3,1:nq) )
  allocate ( rvflx(1:np,1:3,1:nq) )
  !
  !$ACC DATA CREATE(tusp,lvflx,rvflx) &
  !$ACC&     PRESENT(facenrm,face_cv,face_dcvdx,return_value)
  !
  ! Get the left viscous fluxes
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do m = 1,nq
    do n = n1,n2
      tusp(n,m) = face_cv(n,m,1)
    end do
  end do
  call visc_flx_cv_pts_gpu_3d(n1,n2,tusp,face_dcvdx(:,:,:,1),lvflx)
  !
  ! Get the right viscous fluxes
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do m = 1,nq
    do n = n1,n2
      tusp(n,m) = face_cv(n,m,2)
    end do
  end do
  call visc_flx_cv_pts_gpu_3d(n1,n2,tusp,face_dcvdx(:,:,:,2),rvflx)
  !
  ! Use the LDG method to compute the common viscous flux vector
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
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
  deallocate ( tusp )
  !
end subroutine common_vflx_ave_fluxes_gpu_3d
!
!###############################################################################
!
pure subroutine common_vflx_ave_fluxes_gpu_2d(n1,n2,facenrm,face_cv, &
                                              face_dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar, only : LDG_tau, LDG_beta
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
  real(wp) :: rnx,rny
  real(wp) :: unx,uny
  real(wp) :: lcv,lfx,lfy
  real(wp) :: rcv,rfx,rfy
  real(wp) :: vavex,vavey
  real(wp) :: ds,vjump,penalty,vflux
  real(wp) :: vflx_x,vflx_y
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: tusp(:,:)
  real(wp), allocatable :: lvflx(:,:,:)
  real(wp), allocatable :: rvflx(:,:,:)
  !
continue
  !
  np = size(face_dcvdx,dim=1)
  !
  allocate ( tusp(1:np,1:nq) )
  allocate ( lvflx(1:np,1:2,1:nq) )
  allocate ( rvflx(1:np,1:2,1:nq) )
  !
  !$ACC DATA CREATE(tusp,lvflx,rvflx) &
  !$ACC&     PRESENT(facenrm,face_cv,face_dcvdx,return_value)
  !
  ! Get the left viscous fluxes
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do m = 1,nq
    do n = n1,n2
      tusp(n,m) = face_cv(n,m,1)
    end do
  end do
  call visc_flx_cv_pts_gpu_2d(n1,n2,tusp,face_dcvdx(:,:,:,1),lvflx)
  !
  ! Get the right viscous fluxes
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do m = 1,nq
    do n = n1,n2
      tusp(n,m) = face_cv(n,m,2)
    end do
  end do
  call visc_flx_cv_pts_gpu_2d(n1,n2,tusp,face_dcvdx(:,:,:,2),rvflx)
  !
  ! Use the LDG method to compute the common viscous flux vector
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
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
  deallocate ( tusp )
  !
end subroutine common_vflx_ave_fluxes_gpu_2d
!
!###############################################################################
!
pure subroutine common_vflx_ave_fluxes_gpu_nd(n1,n2,facenrm,face_cv, &
                                              face_dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar, only : LDG_tau, LDG_beta
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
  real(wp) :: rnx,rny,rnz
  real(wp) :: unx,uny,unz
  real(wp) :: lcv,lfx,lfy,lfz
  real(wp) :: rcv,rfx,rfy,rfz
  real(wp) :: vavex,vavey,vavez
  real(wp) :: ds,vjump,penalty,vflux
  real(wp) :: vflx_x,vflx_y,vflx_z
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: tusp(:,:)
  real(wp), allocatable :: lvflx(:,:,:)
  real(wp), allocatable :: rvflx(:,:,:)
  !
continue
  !
  np = size(face_dcvdx,dim=1)
  nr = size(face_dcvdx,dim=2)
  !
  allocate ( tusp(1:np,1:nq) )
  allocate ( lvflx(1:np,1:nr,1:nq) )
  allocate ( rvflx(1:np,1:nr,1:nq) )
  !
  !$ACC DATA CREATE(tusp,lvflx,rvflx) &
  !$ACC&     PRESENT(facenrm,face_cv,face_dcvdx,return_value)
  !
  ! Get the left viscous fluxes
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do m = 1,nq
    do n = n1,n2
      tusp(n,m) = face_cv(n,m,1)
    end do
  end do
  call visc_flx_cv_pts_gpu_nd(n1,n2,tusp,face_dcvdx(:,:,:,1),lvflx)
  !
  ! Get the right viscous fluxes
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
  do m = 1,nq
    do n = n1,n2
      tusp(n,m) = face_cv(n,m,2)
    end do
  end do
  call visc_flx_cv_pts_gpu_nd(n1,n2,tusp,face_dcvdx(:,:,:,2),rvflx)
  !
  ! Use the LDG method to compute the common viscous flux vector
  !
  if (nr == 3) then
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
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
  else if (nr == 2) then
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
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
  else
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
    do m = 1,nq
      do n = n1,n2
        !
        rnx = facenrm(1,n)
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
  end if
  !
  !$ACC END DATA
  deallocate ( rvflx )
  deallocate ( lvflx )
  deallocate ( tusp )
  !
end subroutine common_vflx_ave_fluxes_gpu_nd
!
!###############################################################################
!
pure subroutine common_vflx_ave_sol_and_grad_gpu_3d(n1,n2,facenrm,face_cv, &
                                                    face_dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar, only : LDG_tau, LDG_beta
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
  real(wp), allocatable :: com_dcvdx(:,:,:)
  real(wp), allocatable :: com_vflx(:,:,:)
  !
continue
  !
  np = size(face_dcvdx,dim=1)
  !
  allocate ( com_cv(1:np,1:nq) )
  allocate ( com_dcvdx(1:np,1:3,1:nq) )
  allocate ( com_vflx(1:np,1:3,1:nq) )
  !
  !$ACC DATA CREATE(com_cv,com_dcvdx,com_vflx) &
  !$ACC&     PRESENT(facenrm,face_cv,face_dcvdx,return_value)
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
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
  ! Compute the common viscous flux vector using
  ! the common solution and common gradient
  !
  call visc_flx_cv_pts_gpu_3d(n1,n2,com_cv,com_dcvdx,com_vflx)
  !
  ! Subtract the common viscous fluxes from the common
  ! inviscid fluxes already stored within return_value
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
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
  deallocate ( com_cv )
  !
end subroutine common_vflx_ave_sol_and_grad_gpu_3d
!
!###############################################################################
!
pure subroutine common_vflx_ave_sol_and_grad_gpu_2d(n1,n2,facenrm,face_cv, &
                                                    face_dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar, only : LDG_tau, LDG_beta
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
  real(wp), allocatable :: com_dcvdx(:,:,:)
  real(wp), allocatable :: com_vflx(:,:,:)
  !
continue
  !
  np = size(face_dcvdx,dim=1)
  !
  allocate ( com_cv(1:np,1:nq) )
  allocate ( com_dcvdx(1:np,1:2,1:nq) )
  allocate ( com_vflx(1:np,1:2,1:nq) )
  !
  !$ACC DATA CREATE(com_cv,com_dcvdx,com_vflx) &
  !$ACC&     PRESENT(facenrm,face_cv,face_dcvdx,return_value)
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
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
  ! Compute the common viscous flux vector using
  ! the common solution and common gradient
  !
  call visc_flx_cv_pts_gpu_2d(n1,n2,com_cv,com_dcvdx,com_vflx)
  !
  ! Subtract the common viscous fluxes from the common
  ! inviscid fluxes already stored within return_value
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
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
  deallocate ( com_cv )
  !
end subroutine common_vflx_ave_sol_and_grad_gpu_2d
!
!###############################################################################
!
pure subroutine common_vflx_ave_sol_and_grad_gpu_nd(n1,n2,facenrm,face_cv, &
                                                    face_dcvdx,return_value)
  !
  !.. Use Statements ..
  use ovar, only : LDG_tau, LDG_beta
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
  real(wp), allocatable :: com_dcvdx(:,:,:)
  real(wp), allocatable :: com_vflx(:,:,:)
  !
continue
  !
  np = size(face_dcvdx,dim=1)
  nr = size(face_dcvdx,dim=2)
  !
  allocate ( com_cv(1:np,1:nq) )
  allocate ( com_dcvdx(1:np,1:nr,1:nq) )
  allocate ( com_vflx(1:np,1:nr,1:nq) )
  !
  !$ACC DATA CREATE(com_cv,com_dcvdx,com_vflx) &
  !$ACC&     PRESENT(facenrm,face_cv,face_dcvdx,return_value)
  !
  if (nr == 3) then
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
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
  else if (nr == 2) then
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
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
  else
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
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
  end if
  !
  ! Compute the common viscous flux vector using
  ! the common solution and common gradient
  !
  call visc_flx_cv_pts_gpu_nd(n1,n2,com_cv,com_dcvdx,com_vflx)
  !
  ! Subtract the common viscous fluxes from the common
  ! inviscid fluxes already stored within return_value
  !
  if (nr == 3) then
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
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
  else if (nr == 2) then
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
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
  else
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) COLLAPSE(2)
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
  end if
  !
  !$ACC END DATA
  deallocate ( com_vflx )
  deallocate ( com_dcvdx )
  deallocate ( com_cv )
  !
end subroutine common_vflx_ave_sol_and_grad_gpu_nd
!
!###############################################################################
!
pure subroutine common_visc_flx_sub(unrm,ul,duldx,ur,durdx,return_value)
  !$ACC ROUTINE SEQ
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : visc_variables_to_ave
  use ovar,   only : LDG_tau, LDG_beta
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),      intent(in) :: unrm
  real(wp), dimension(1:nq),      intent(in) :: ul,ur
  real(wp), dimension(1:nr,1:nq), intent(in) :: duldx,durdx
  real(wp), dimension(1:nq),   intent(inout) :: return_value
  !
  !.. Local Scalars ..
  integer  :: l,m
  real(wp) :: vjump,qjump,penalty
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq)      :: ucom
  real(wp), dimension(1:nr,1:nq) :: qcom
  real(wp), dimension(1:nq,1:nr) :: lvflx
  real(wp), dimension(1:nq,1:nr) :: rvflx
  real(wp), dimension(1:nr,1:nq) :: vflux
  real(wp), dimension(1:nq,1:nr) :: fvcom
  !
continue
  !
  if (visc_variables_to_ave == visc_interface_fluxes) then
    !
    ! Get the left and right viscous fluxes
    !
    call visc_flx_cv(ul,duldx,lvflx)
   !lvflx(1:nq,1:nr) = visc_flx_cv(ul,duldx)
    call visc_flx_cv(ur,durdx,rvflx)
   !rvflx(1:nq,1:nr) = visc_flx_cv(ur,durdx)
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
    ! Compute the common viscous flux vector using
    ! the common solution and common gradient
    !
    call visc_flx_cv(ucom,qcom,fvcom)
   !fvcom(1:nq,1:nr) = visc_flx_cv(ucom,qcom)
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
pure function common_visc_flx_cv(unrm,ul,duldx,ur,durdx) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : visc_variables_to_ave
  use ovar,   only : LDG_tau, LDG_beta
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nr),         intent(in) :: unrm
  real(wp), dimension(1:nq),      intent(in) :: ul,ur
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
  real(wp), dimension(1:nr,1:nq) :: qcom
  real(wp), dimension(1:nq,1:nr) :: lvflx
  real(wp), dimension(1:nq,1:nr) :: rvflx
  real(wp), dimension(1:nr,1:nq) :: vflux
  real(wp), dimension(1:nq,1:nr) :: fvcom
  !
continue
  !
  if (visc_variables_to_ave == visc_interface_fluxes) then
    !
    ! Get the left and right viscous fluxes
    !
    call visc_flx_cv(ul,duldx,lvflx)
   !lvflx(1:nq,1:nr) = visc_flx_cv(ul,duldx)
    call visc_flx_cv(ur,durdx,rvflx)
   !rvflx(1:nq,1:nr) = visc_flx_cv(ur,durdx)
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
    ! Compute the common viscous flux vector using
    ! the common solution and common gradient
    !
    call visc_flx_cv(ucom,qcom,fvcom)
   !fvcom(1:nq,1:nr) = visc_flx_cv(ucom,qcom)
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
