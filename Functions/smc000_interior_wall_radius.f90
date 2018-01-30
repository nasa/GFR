pure function smc000_interior_wall_radius(x,eps_1,eps_2) result(return_value)
  !
  ! NOTE: This function assumes the incoming xyz is in meters
  !       since all the parameter coordinates are in meters.
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: x
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: eps_1
  real(wp), optional, intent(in) :: eps_2
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: e1,e2
  !
  !.. Local Parameters ..
 !real(wp), parameter :: paa(2) = [-0.250825_wp,          0.0762_wp]
  real(wp), parameter :: pa(2) = [-0.1957832_wp,         0.0762_wp]
  real(wp), parameter :: pf(2) = [-0.199486039157406_wp, 0.0510713513659322_wp]
  real(wp), parameter :: pb(2) = [-0.187954680464026_wp, 0.07370291708031857_wp]
  real(wp), parameter :: pc(2) = [-0.138839103315909_wp, 0.04867728055674608_wp]
  real(wp), parameter :: pg(2) = [-0.069650951155604_wp, 0.1844666748430524_wp]
  real(wp), parameter :: pd(2) = [-0.0830072_wp,         0.03265306999550994_wp]
  real(wp), parameter :: pe(2) = [ 0.0_wp,               0.0254_wp ]
  !
 !real(wp), parameter :: lg2 = dot(pf-pa,pf-pa)
  real(wp), parameter :: lg2 = (pf(1)-pa(1))*(pf(1)-pa(1)) + &
                               (pf(2)-pa(2))*(pf(2)-pa(2))
 !real(wp), parameter :: lh2 = dot(pg-pc,pg-pc)
  real(wp), parameter :: lh2 = (pg(1)-pc(1))*(pg(1)-pc(1)) + &
                               (pg(2)-pc(2))*(pg(2)-pc(2))
  !
  real(wp), parameter :: mcd = (pb(2)-pc(2)) / (pb(1)-pc(1))
  real(wp), parameter :: bcd = (pb(1)*pc(2)-pb(2)*pc(1)) / (pb(1)-pc(1))
  !
  real(wp), parameter :: mef = (pd(2)-pe(2)) / (pd(1)-pe(1))
  real(wp), parameter :: bef = (pd(1)*pe(2)-pd(2)*pe(1)) / (pd(1)-pe(1))
  !
continue
  !
  e1 = zero
  if (present(eps_1)) e1 = eps_1
  !
  e2 = zero
  if (present(eps_2)) e2 = eps_2
  !
  if (x <= pa(1)) then
    return_value = pa(2)
  else if (x > pa(1) .and. x < pb(1)) then
    return_value = pf(2) + sqrt( lg2 - (x-pf(1))**2 ) &
                 + e1 ! add small value if optional argument is present
  else if (x >= pb(1) .and. x <= pc(1)) then
    return_value = mcd*x + bcd &
                 + e1 ! add small value if optional argument is present
  else if (x > pc(1) .and. x < pd(1)) then
    return_value = pg(2) - sqrt( lh2 - (x-pg(1))**2 ) &
                 + e2 ! add small value if optional argument is present
  else if (x >= pd(1) .and. x < pe(1)) then
    return_value = mef*x + bef
  else if (x >= pe(1)) then
    return_value = pe(2)
  end if
  !
end function smc000_interior_wall_radius
!
!###############################################################################
!
pure function point_is_in_smc000_nozzle(xyz) result(return_value)
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
  real(wp), parameter :: x_in   = -0.250825_wp
  real(wp), parameter :: y_in   =  0.0762_wp
  real(wp), parameter :: x_exit =  0.0_wp
  !
  real(wp), parameter :: eps_1 = eps6
  real(wp), parameter :: eps_2 = eps3
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
      radius = smc000_interior_wall_radius(x,eps_1,eps_2)
      return_value = ( r <= radius )
    end if
    !
  end if
  !
end function point_is_in_smc000_nozzle
!
!###############################################################################
!
