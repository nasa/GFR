pure function grad_cv_to_grad_pv(rusp,dcdx,dcdy) result(return_value)
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : rgasref_nd,gam
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
  use geovar, only : nr
  use ovar,   only : rgasref_nd,gam
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
  real(wp) :: cv_inv,rho_inv,rho_inv2,temp
  !
  !.. Local Arrays ..
  real(wp) :: var(1:nq)
  !
continue
  !
  cv_inv = (gam - one) / rgasref_nd
  !
  rho_inv = one/rusp(nec)
  rho_inv2 = rho_inv * rho_inv
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
