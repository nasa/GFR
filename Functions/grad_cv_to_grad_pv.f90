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
  real(wp) :: cv_inv,rho_inv,rho_inv2,temp
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
      return_value(l,neqns) = return_value(l,neqns) - rho_inv*rusp(m)*return_value(l,m)
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
!
