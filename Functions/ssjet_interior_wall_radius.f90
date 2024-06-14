pure function ssjet_interior_wall_radius(xin,opt_eps) result(return_value)
  !
  ! NOTE: This function assumes the incoming xyz is in meters
  !       since all the parameter coordinates are in meters.
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xin
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: opt_eps
  !
  !.. Function Results ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: x,radius
  !
  !.. Local Parameters ..
  ! Ending locations for each section
  real(wp), parameter :: pa = -12.91_wp
  real(wp), parameter :: pb = -12.46_wp
  real(wp), parameter :: pc = -11.5648_wp
  real(wp), parameter :: pd =  -6.01697_wp
  real(wp), parameter :: pe =  -3.02368_wp
  real(wp), parameter :: pf =  -1.40164_wp
  real(wp), parameter :: pg =  zero
  ! Curvefit coefficients for each section
  ! curve fit of geometry found by Jim DeBonis using excel spreadsheet
  real(wp), parameter :: ca(1) = [three]
  real(wp), parameter :: cb(2) = [2.024545_wp,0.951512_wp]
  real(wp), parameter :: cc(2) = [-3.4586_wp,-0.50952_wp]
  real(wp), parameter :: cd(7) = [ 56.760_wp, &
                                   36.830_wp, &
                                    9.6905_wp, &
                                    1.2833_wp, &
                                    0.088276_wp, &
                                    0.0028479_wp, &
                                    2.9499E-05_wp]
  real(wp), parameter :: ce(7) = [ 1.0413_wp, &
                                   0.14554_wp, &
                                   0.036498_wp, &
                                   0.0017165_wp, &
                                  -3.3744E-05_wp, &
                                  -1.5264E-06_wp, &
                                  -5.9181E-09_wp]
  real(wp), parameter :: cf(7) = [  9.6209_wp, &
                                   24.930_wp, &
                                   29.434_wp, &
                                   18.170_wp, &
                                    6.1608_wp, &
                                    1.0881_wp, &
                                    0.078353_wp]
  real(wp), parameter :: cg(7) = [ 0.99999_wp, &
                                  -0.00034635_wp, &
                                  -0.010837_wp, &
                                   0.016500_wp, &
                                   0.0076309_wp, &
                                   0.0019854_wp, &
                                   0.00030604_wp]
  !
continue
  !
  x = xin * meters_to_inches
  !
  if (x <= pa) then                   ! section a
    radius = three
  else if (x > pa .and. x <= pb) then ! section b
    radius = cb(1) + sqrt(cb(2) - (x-pa)**2)
  else if (x > pb .and. x <= pc) then ! section c
    radius = cc(1) + cc(2)*x
  else if (x > pc .and. x <= pd) then ! section d
    radius = cd(1) + cd(2)*x + cd(3)*x**2 + cd(4)*x**3 + cd(5)*x**4 &
                                          + cd(6)*x**5 + cd(7)*x**6
  else if (x > pd .and. x <= pe) then ! section e
    radius = ce(1) + ce(2)*x + ce(3)*x**2 + ce(4)*x**3 + ce(5)*x**4 &
                                          + ce(6)*x**5 + ce(7)*x**6
  else if (x > pe .and. x <= pf) then ! section f
    radius = cf(1) + cf(2)*x + cf(3)*x**2 + cf(4)*x**3 + cf(5)*x**4 &
                                          + cf(6)*x**5 + cf(7)*x**6
  else if (x > pf .and. x < pg) then  ! section g
    radius = cg(1) + cg(2)*x + cg(3)*x**2 + cg(4)*x**3 + cg(5)*x**4 &
                                          + cg(6)*x**5 + cg(7)*x**6
  else if (x >= pg) then
    radius = one
  end if
  !
  if (present(opt_eps)) then
    ! Scale the radius by a small value if optional argument is present
    radius = radius*(one+opt_eps)
  end if
  !
  return_value = radius * inches_to_meters
  !
end function ssjet_interior_wall_radius
!
!###############################################################################
!
pure function point_is_in_ssjet_nozzle(xyz) result(return_value)
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
  real(wp), parameter :: y_in   = 0.0762_wp ! 3 inches in meters
  real(wp), parameter :: x_exit = zero
  !
  real(wp), parameter :: eps_1 = eps6
  real(wp), parameter :: eps_2 = eps2
  !
continue
  !
  return_value = fals
  !
  x = xyz(1)
  !
  if (x < x_exit+eps_1) then
    !
    r = norm2(xyz(2:))
    !
    ! Only check if this point is inside the nozzle if the radius for the
    ! current point is less than the largest radius of the nozzle (3 inches)
    if (r <= y_in+eps_1) then
      radius = ssjet_interior_wall_radius(x,opt_eps=eps_2)
      return_value = ( r <= radius )
    end if
    !
  end if
  !
end function point_is_in_ssjet_nozzle
!
!###############################################################################
!
