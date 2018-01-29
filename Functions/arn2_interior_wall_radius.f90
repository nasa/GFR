pure function arn2_interior_wall_radius(x,eps_1) result(return_value)
  !
  ! NOTE: This function assumes the incoming xyz is in meters
  !       since all the parameter coordinates are in meters.
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: x
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: eps_1
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: e1
  !
  !.. Local Parameters ..
  real(wp), parameter :: pa(2) = [-0.196596_wp,0.0762_wp]
  real(wp), parameter :: pb(2) = [-0.006350_wp,0.0254_wp]
  !
  real(wp), parameter :: c0 =  0.0256084_wp
  real(wp), parameter :: c1 =  0.0667531_wp
  real(wp), parameter :: c2 =  5.4427200_wp
  real(wp), parameter :: c3 = 19.5023000_wp
  !
continue
  !
  e1 = zero
  if (present(eps_1)) e1 = eps_1
  !
  if (x <= pa(1)) then
    return_value = pa(2)
  else if (x >= pb(1)) then
    return_value = pb(2)
  else
    return_value = c0 + c1*x + c2*x*x + c3*x*x*x + e1
  end if
  !
end function arn2_interior_wall_radius
!
!###############################################################################
!
pure function point_is_in_arn2_nozzle(xyz) result(return_value)
  !
  ! NOTE: This function assumes the incoming xyz is in meters
  !       since all the parameter coordinates are in meters.
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:)
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: x,r,radius
  !
  !.. Local Parameters ..
  real(wp), parameter :: x_in   = -0.196596_wp
  real(wp), parameter :: y_in   =  0.0762_wp
  real(wp), parameter :: x_exit =  0.0_wp
  !
  real(wp), parameter :: eps_1 = eps3
  !
continue
  !
  return_value = fals
  !
  x = xyz(1)
  !
  if (x >= x_in-eps_1 .and. x < x_exit-eps_1) then
    !
    r = norm2(xyz(2:))
    !
    if (r <= y_in+eps_1) then
      radius = arn2_interior_wall_radius(x,eps_1)
      return_value = ( r <= radius )
    end if
    !
  end if
  !
end function point_is_in_arn2_nozzle
!
!###############################################################################
!
