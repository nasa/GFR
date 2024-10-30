!pure subroutine grad_cv_to_grad_pv_gpu(tusp,tdcvdx,return_value, &
!                                       tusp_includes_rho_inv)
!  !
!  ! This routine computes the discontinuous gradient at each solution point
!  !
!  !.. Use Statements ..
!  use ovar, only : rgasref_nd,gam
!  !
!  !.. Formal Arguments ..
!  real(wp), intent(inout) :: tusp(:,:)
!  real(wp),    intent(in) :: tdcvdx(:,:,:)
!  real(wp), intent(inout) :: return_value(:,:,:)
!  !
!  !.. Optional Arguments ..
!  logical(lk), optional, intent(in) :: tusp_includes_rho_inv
!  !
!  !.. Local Scalars ..
!  integer  :: l,m,n,np,nr
!  real(wp) :: cv_inv,rho_inv
!  !
!  logical(lk) :: tusp_needs_rho_inv
!  !
!continue
!  !
!  tusp_needs_rho_inv = true
!  if (present(tusp_includes_rho_inv)) then
!    tusp_needs_rho_inv = .not. tusp_includes_rho_inv
!  end if
!  !
!  np = size(tdcvdx,dim=1)
!  nr = size(tdcvdx,dim=2)
!  !
!  cv_inv = (gam - one) / rgasref_nd
!  !
!  ! Embed 1/rho in the all variables except density
!  !
!  if (tusp_needs_rho_inv) then
!    !$ACC PARALLEL LOOP COLLAPSE(2) &
!    !$ACC&              PRESENT(tusp)
!    do m = nmb,nq
!      do n = 1,np
!        tusp(n,m) = tusp(n,m)/tusp(n,nec)
!      end do
!    end do
!  end if
!  !
!  ! Copy the density gradient
!  !
!  !$ACC PARALLEL LOOP COLLAPSE(2) &
!  !$ACC&              PRESENT(tdcvdx,return_value)
!  do l = 1,nr
!    do n = 1,np
!      return_value(n,l,nec) = tdcvdx(n,l,nec)
!    end do
!  end do
!  !
!  ! Using the chain rule:
!  !
!  ! 1. Compute the velocity gradients (m = 2,nq-1)
!  ! 2. Compute the total energy gradient and store in the
!  !    location for the temperature gradient (m = nq)
!  !
!  !$ACC PARALLEL LOOP COLLAPSE(3) &
!  !$ACC&              PRESENT(tusp,tdcvdx,return_value)
!  do m = nmb,nq
!    do l = 1,nr
!     !dpvdx(l,m) = rho_inv*dcvdx(l,m) - rho_inv*var(m)*dcvdx(l,nec)
!      do n = 1,np
!        rho_inv = one/tusp(n,nec)
!        return_value(n,l,m) = rho_inv*tdcvdx(n,l,m) - &
!                              rho_inv*tdcvdx(n,l,nec)*tusp(n,m)
!      end do
!    end do
!  end do
!  !
!  ! The location for the temperature gradient currently contains
!  ! the total energy gradient from the previous loop.
!  ! Now use the velocity gradients to subtract the gradient of
!  ! kinetic-energy-per-unit-mass from the total energy gradient.
!  !
!  !
!#ifdef NOT_SURE_ABOUT_LOOP
!  do m = nmb,nme
!    !$ACC PARALLEL LOOP COLLAPSE(2) &
!    !$ACC&              PRESENT(tusp,return_value)
!    do l = 1,nr
!     !dpvdx(l,nee) = dpvdx(l,nee) - var(m)*dpvdx(l,m)
!      do n = 1,np
!        return_value(n,l,nee) = return_value(n,l,nee) &
!                              - return_value(n,l,m)*tusp(n,m)
!      end do
!    end do
!  end do
!#else
!  if (nr == 3) then
!    !$ACC PARALLEL LOOP COLLAPSE(2) &
!    !$ACC&              PRESENT(tusp,return_value)
!    do l = 1,3
!     !dpvdx(l,nee) = dpvdx(l,nee) - var(m)*dpvdx(l,m)
!      do n = 1,np
!        return_value(n,l,nee) = return_value(n,l,nee) &
!                              - return_value(n,l,nmx)*tusp(n,nmx) &
!                              - return_value(n,l,nmy)*tusp(n,nmy) &
!                              - return_value(n,l,nmz)*tusp(n,nmz)
!      end do
!    end do
!  else if (nr == 2) then
!    !$ACC PARALLEL LOOP COLLAPSE(2) &
!    !$ACC&              PRESENT(tusp,return_value)
!    do l = 1,2
!     !dpvdx(l,nee) = dpvdx(l,nee) - var(m)*dpvdx(l,m)
!      do n = 1,np
!        return_value(n,l,nee) = return_value(n,l,nee) &
!                              - return_value(n,l,nmx)*tusp(n,nmx) &
!                              - return_value(n,l,nmy)*tusp(n,nmy)
!      end do
!    end do
!  else
!    !$ACC PARALLEL LOOP PRESENT(tusp,return_value)
!    do n = 1,np
!      return_value(n,1,nee) = return_value(n,1,nee) &
!                            - return_value(n,1,nmx)*tusp(n,nmx)
!    end do
!  end if
!#endif
!  !
!  ! Finish computing the temperature gradient by dividing through by
!  ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
!  !
! !dpvdx(1:nr,nee) = dpvdx(1:nr,nee) * cv_inv
!  !$ACC PARALLEL LOOP COLLAPSE(2) &
!  !$ACC&              PRESENT(return_value)
!  do l = 1,nr
!    do n = 1,np
!      return_value(n,l,nee) = return_value(n,l,nee) * cv_inv
!    end do
!  end do
!  !
!  ! Remove the embedded 1/rho in the all variables except density
!  !
!  if (tusp_needs_rho_inv) then
!    !$ACC PARALLEL LOOP COLLAPSE(2) &
!    !$ACC&              PRESENT(tusp)
!    do m = nmb,nq
!      do n = 1,np
!        tusp(n,m) = tusp(n,m)*tusp(n,nec)
!      end do
!    end do
!  end if
!  !
!! rho_inv = one/rusp(nec)
!! !
!! var(1:nq) = rho_inv*rusp(1:nq)
!! !
!! ! Copy the density gradient
!! !
!! return_value(1:nr,nec) = dcvdx(1:nr,nec)
!! !
!! ! Using the chain rule:
!! !
!! ! 1. Compute the velocity gradients (m = 2,nq-1)
!! ! 2. Compute the total energy gradient and store in the
!! !    location for the temperature gradient (m = nq)
!! !
!! do m = nmb,nee
!!   do l = 1,nr
!!     return_value(l,m) = rho_inv*dcvdx(l,m) - rho_inv*var(m)*dcvdx(l,nec)
!!   end do
!! end do
!! !
!! ! The location for the temperature gradient currently contains
!! ! the total energy gradient from the previous loop.
!! ! Now use the velocity gradients to subtract the gradient of
!! ! kinetic-energy-per-unit-mass from the total energy gradient.
!! !
!! do m = nmb,nme
!!   do l = 1,nr
!!     return_value(l,nee) = return_value(l,nee) - var(m)*return_value(l,m)
!!   end do
!! end do
!! !
!! ! Finish computing the temperature gradient by dividing through by
!! ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
!! !
!! return_value(1:nr,nee) = return_value(1:nr,nee) * cv_inv
!! !
!end subroutine grad_cv_to_grad_pv_gpu
!
!###############################################################################
!
pure function grad_cv_to_grad_pv(rusp,dcdx,dcdy) result(return_value)
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  use ovar,        only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: rusp
  real(wp), dimension(:,:), intent(in) :: dcdx
  real(wp), dimension(:,:), intent(in) :: dcdy
  !
  !.. Function Return Value ..
  real(wp), dimension(1:nr,1:nq,1:size(rusp,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: k,m,n,np
  real(wp) :: cv_inv,rho_inv,rho_inv2
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq,1:size(rusp,dim=2)) :: dpdx,dpdy
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
    rho_inv = one/rusp(nec,k)
    rho_inv2 = rho_inv * rho_inv
    !
    ! Copy the density gradient
    !
    dpdx(nec,k) = dcdx(nec,k)
    dpdy(nec,k) = dcdy(nec,k)
    !
    ! Using the chain rule:
    !
    ! 1. Compute the velocity gradients (m = 2,nq-1)
    ! 2. Compute the total energy gradient and store in the
    !    location for the temperature gradient (m = nq)
    !
    do m = nmb,nee
      dpdx(m,k) = rho_inv*dcdx(m,k) - rho_inv2*rusp(m,k)*dcdx(nec,k)
      dpdy(m,k) = rho_inv*dcdy(m,k) - rho_inv2*rusp(m,k)*dcdy(nec,k)
    end do
    !
    ! The location for the temperature gradient currently contains
    ! the total energy gradient from the previous loop.
    ! Now use the velocity gradients to subtract the gradient of
    ! kinetic-energy-per-unit-mass from the total energy gradient.
    !
    do m = nmb,nme
      dpdx(nee,k) = dpdx(nee,k) - rho_inv*rusp(m,k)*dpdx(m,k)
      dpdy(nee,k) = dpdy(nee,k) - rho_inv*rusp(m,k)*dpdy(m,k)
    end do
    !
    ! Finish computing the temperature gradient by dividing through by
    ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
    !
    dpdx(nee,k) = dpdx(nee,k) * cv_inv
    dpdy(nee,k) = dpdy(nee,k) * cv_inv
    !
  end do
  !
  ! Put primitive variable gradients in the dusp array
  !
  do k = 1,np
    do m = 1,nq
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
  use geomvar_mod, only : nr
  use ovar,        only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(:),   intent(in) :: rusp
  real(wp), dimension(:,:), intent(in) :: dcvdx
  !
  !.. Function Return Value ..
  real(wp), dimension(1:nr,1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer  :: l,m
  real(wp) :: cv_inv,rho_inv
  !
  !.. Local Arrays ..
  real(wp) :: var(1:nq)
  !
continue
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  rho_inv = one/rusp(nec)
  !
  var(1:nq) = rho_inv*rusp(1:nq)
  !
  ! Copy the density gradient
  !
  return_value(1:nr,nec) = dcvdx(1:nr,nec)
  !
  ! Using the chain rule:
  !
  ! 1. Compute the velocity gradients (m = 2,nq-1)
  ! 2. Compute the total energy gradient and store in the
  !    location for the temperature gradient (m = nq)
  !
  do m = nmb,nee
    do l = 1,nr
      return_value(l,m) = rho_inv*dcvdx(l,m) - rho_inv*var(m)*dcvdx(l,nec)
    end do
  end do
  !
  ! The location for the temperature gradient currently contains
  ! the total energy gradient from the previous loop.
  ! Now use the velocity gradients to subtract the gradient of
  ! kinetic-energy-per-unit-mass from the total energy gradient.
  !
  do m = nmb,nme
    do l = 1,nr
      return_value(l,nee) = return_value(l,nee) - var(m)*return_value(l,m)
    end do
  end do
  !
  ! Finish computing the temperature gradient by dividing through by
  ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
  !
  return_value(1:nr,nee) = return_value(1:nr,nee) * cv_inv
  !
end function grad_cv_to_grad_pv_sp
!
!###############################################################################
!
pure subroutine grad_cv_to_grad_pv_sp_sub(rusp,dcvdx,return_value)
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  use ovar,        only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  real(wp), dimension(1:nq),   intent(in) :: rusp
  real(wp), dimension(1:nr,1:nq), intent(in) :: dcvdx
  !
  !.. Function Return Value ..
  real(wp), dimension(1:nr,1:nq), intent(inout) :: return_value
  !
  !.. Local Scalars ..
  integer  :: l,m
  real(wp) :: cv_inv,rho_inv
  !
  !.. Local Arrays ..
  real(wp) :: var(1:nq)
  !
continue
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  rho_inv = one/rusp(nec)
  !
  var(1:nq) = rho_inv*rusp(1:nq)
  !
  ! Copy the density gradient
  !
  return_value(1:nr,nec) = dcvdx(1:nr,nec)
  !
  ! Using the chain rule:
  !
  ! 1. Compute the velocity gradients (m = 2,nq-1)
  ! 2. Compute the total energy gradient and store in the
  !    location for the temperature gradient (m = nq)
  !
  do m = nmb,nee
    do l = 1,nr
      return_value(l,m) = rho_inv*dcvdx(l,m) - rho_inv*var(m)*dcvdx(l,nec)
    end do
  end do
  !
  ! The location for the temperature gradient currently contains
  ! the total energy gradient from the previous loop.
  ! Now use the velocity gradients to subtract the gradient of
  ! kinetic-energy-per-unit-mass from the total energy gradient.
  !
  do m = nmb,nme
    do l = 1,nr
      return_value(l,nee) = return_value(l,nee) - var(m)*return_value(l,m)
    end do
  end do
  !
  ! Finish computing the temperature gradient by dividing through by
  ! the specific-heat-at-constant-volume, i.e. cv = Rgas/(gamma-1)
  !
  return_value(1:nr,nee) = return_value(1:nr,nee) * cv_inv
  !
end subroutine grad_cv_to_grad_pv_sp_sub
!
!###############################################################################
!
