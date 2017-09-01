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
pure subroutine visc_flx_cv_gpu(usp,dusp,return_value)
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
  !$ACC&     COPY(return_value) &
  !$ACC&     COPYIN(usp,dusp)
  !
  !$ACC PARALLEL LOOP
  do n = 1,np
    tusp(n,nec) = usp(nec,n)
  end do
  !
  ! NOTE: EMBED 1/RHO IN ALL OTHER VARIABLES
  !
  !$ACC PARALLEL LOOP COLLAPSE(2)
  do m = nmb,nq
    do n =1,np
      tusp(n,m) = usp(m,n)/tusp(n,nec)
    end do
  end do
  !
  ! Temporarily store the transposed conserved variable
  ! gradients in the return_value array
  !
  !$ACC PARALLEL LOOP COLLAPSE(3)
  do n = 1,np
    do m = 1,nq
      do l = 1,nr
        return_value(n,l,m) = dusp(l,m,n)
      end do
    end do
  end do
  !
 !!$ACC PARALLEL LOOP COLLAPSE(3)
 !do m = 1,nq
 !  do l = 1,nr
 !    do n = 1,np
 !      return_value(n,l,m) = return_value(n,l,m)/tusp(n,nec)
 !    end do
 !  end do
 !end do
  !
  ! Get the gradients of the primitive variables using the 
  ! gradients of the conservative variables and the product rule
  !
 !call grad_cv_to_grad_pv_sp_sub(uu,dcvdx,dpvdx)
  call grad_cv_to_grad_pv_gpu(tusp,return_value,tdpvdx, &
                              tusp_includes_rho_inv=true)
  ! 
  ! Replace energy variable with molecular viscosity
  !
  !$ACC PARALLEL LOOP
  do n = 1,np
    tusp(n,nee) = viscosity_cv_sp( usp(:,n) )
  end do
  !
  ! Density flux is zero
  !
  !$ACC PARALLEL LOOP COLLAPSE(2)
  do l = 1,nr
    do n = 1,np
      return_value(n,l,nec) = zero
    end do
  end do
  !
  ! Momentum equations
 !!
 !! Compute the average extension rate
 !!
 !ave_ext_rate = zero
 !do l = 1,nr
 !  ave_ext_rate = ave_ext_rate + dpvdx(l,nmb-1+l)
 !end do
 !ave_ext_rate = two3 * ave_ext_rate
 !!
 !! Compute the viscous fluxes for the momentum equations which are just
 !! the components of the viscous stress tensor.
 !! NOTE: We are initially excluding the molecular viscosity term which
 !!       will be multiplied to the whole tau array after these loops.
 !!
 !do j = 1,nr
 !  do i = 1,nr
 !      return_value(nmb-1+i,j) = dpvdx(i,nmb-1+j) + dpvdx(j,nmb-1+i)
 !  end do
 !end do
 !!
 !! Get the velocity vector
 !!
 !vel(:) = uu(nmb:nme)/uu(nec)
 !!
 !! Compute the molecular viscosity and molecular conductivity
 !!
 !mvsc = viscosity_cv_sp(uu)
 !mcdt = cpref_nd * mvsc / Pr
 !!
 !! Subtract the average extension rate from the normal stresses
 !!
 !do l = 1,nr
 !  return_value(nmb-1+l,l) = return_value(nmb-1+l,l) - ave_ext_rate
 !end do
 !!
 !! Multiply the viscous stresses by the molecular viscosity
 !!
 !return_value(nmb:nme,1:nr) = mvsc * return_value(nmb:nme,1:nr)
  !
  if (nr == 3) then
    !
    !$ACC PARALLEL LOOP
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
    !$ACC PARALLEL LOOP
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
    !$ACC PARALLEL LOOP
    do n = 1,np
      return_value(n,1,nmx) = two * two3 * tusp(n,nee) * tdpvdx(n,1,nmx)
    end do
    !
  end if
  !
  ! Energy Equation
  !
 !do l = 1,nr
 !  return_value(nq,l) = mcdt*dpvdx(l,nq) + &
 !                       dot( vel(1:nr) , return_value(nmb:nme,l) )
 !end do
  if (nr == 3) then
    !
    !$ACC PARALLEL LOOP COLLAPSE(2)
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
    !$ACC PARALLEL LOOP COLLAPSE(2)
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
    !$ACC PARALLEL LOOP COLLAPSE(2)
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
end subroutine visc_flx_cv_gpu
!
!###############################################################################
!
pure subroutine grad_cv_to_grad_pv_gpu(tusp,tdcvdx,return_value, &
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
    !$ACC PARALLEL LOOP COLLAPSE(2)
    do m = nmb,nq
      do n = 1,np
        tusp(n,m) = tusp(n,m)/tusp(n,nec)
      end do
    end do
  end if
  !
  ! Copy the density gradient
  !
  !$ACC PARALLEL LOOP COLLAPSE(2)
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
  !$ACC PARALLEL LOOP COLLAPSE(3)
  do m = nmb,nq
    do l = 1,nr
     !dpvdx(l,m) = rho_inv*dcvdx(l,m) - rho_inv*var(m)*dcvdx(l,nec)
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
    !$ACC PARALLEL LOOP COLLAPSE(2)
    do l = 1,nr
     !dpvdx(l,nee) = dpvdx(l,nee) - var(m)*dpvdx(l,m)
      do n = 1,np
        return_value(n,l,nee) = return_value(n,l,nee) &
                              - return_value(n,l,m)*tusp(n,m)
      end do
    end do
  end do
#else
  if (nr == 3) then
    !$ACC PARALLEL LOOP COLLAPSE(2)
    do l = 1,3
     !dpvdx(l,nee) = dpvdx(l,nee) - var(m)*dpvdx(l,m)
      do n = 1,np
        return_value(n,l,nee) = return_value(n,l,nee) &
                              - return_value(n,l,nmx)*tusp(n,nmx) &
                              - return_value(n,l,nmy)*tusp(n,nmy) &
                              - return_value(n,l,nmz)*tusp(n,nmz)
      end do
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP COLLAPSE(2)
    do l = 1,2
     !dpvdx(l,nee) = dpvdx(l,nee) - var(m)*dpvdx(l,m)
      do n = 1,np
        return_value(n,l,nee) = return_value(n,l,nee) &
                              - return_value(n,l,nmx)*tusp(n,nmx) &
                              - return_value(n,l,nmy)*tusp(n,nmy)
      end do
    end do
  else
    !$ACC PARALLEL LOOP
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
 !dpvdx(1:nr,nee) = dpvdx(1:nr,nee) * cv_inv
  !$ACC PARALLEL LOOP COLLAPSE(2)
  do l = 1,nr
    do n = 1,np
      return_value(n,l,nee) = return_value(n,l,nee) * cv_inv
    end do
  end do
  !
  ! Remove the embedded 1/rho in the all variables except density
  !
  if (tusp_needs_rho_inv) then
    !$ACC PARALLEL LOOP COLLAPSE(2)
    do m = nmb,nq
      do n = 1,np
        tusp(n,m) = tusp(n,m)*tusp(n,nec)
      end do
    end do
  end if
  !$ACC END DATA
  !
! rho_inv = one/rusp(nec)
! !
! var(1:nq) = rho_inv*rusp(1:nq)
! !
! ! Copy the density gradient
! !
! return_value(1:nr,nec) = dcvdx(1:nr,nec)
! !
! ! Using the chain rule:
! !
! ! 1. Compute the velocity gradients (m = 2,nq-1)
! ! 2. Compute the total energy gradient and store in the
! !    location for the temperature gradient (m = nq)
! !
! do m = nmb,nee
!   do l = 1,nr
!     return_value(l,m) = rho_inv*dcvdx(l,m) - rho_inv*var(m)*dcvdx(l,nec)
!   end do
! end do
! !
! ! The location for the temperature gradient currently contains
! ! the total energy gradient from the previous loop.
! ! Now use the velocity gradients to subtract the gradient of
! ! kinetic-energy-per-unit-mass from the total energy gradient.
! ! 
! do m = nmb,nme
!   do l = 1,nr
!     return_value(l,nee) = return_value(l,nee) - var(m)*return_value(l,m)
!   end do
! end do
! !
! ! Finish computing the temperature gradient by dividing through by
! ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
! !
! return_value(1:nr,nee) = return_value(1:nr,nee) * cv_inv
! !
end subroutine grad_cv_to_grad_pv_gpu
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
  return_value = zero
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
  vflx = zero
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
