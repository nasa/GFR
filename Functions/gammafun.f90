!
!###############################################################################
!
elemental function gammafun(x) result(return_value)
  !
  real(wp), intent(in) :: x
  !
  real(wp) :: return_value
  !
  integer   :: ind,k
  real(wp) :: xx,pr,s,g
  !
  real(wp), parameter :: s_start = 0.426401432711220868_wp
  real(wp), parameter :: gamma_constant(1:11) = &
             [ -0.524741987629368444e+00_wp,  0.116154405493589130e+00_wp, &
               -0.765978624506602380e-02_wp,  0.899719449391378898e-04_wp, &
               -0.194536980009534621e-07_wp,  0.199382839513630987e-10_wp, &
               -0.204209590209541319e-11_wp,  0.863896817907000175e-13_wp, &
                0.152237501608472336e-13_wp, -0.825725175277719950e-14_wp, &
                0.299734782205224610e-14_wp ]
  !
continue
  !
  xx = x
  return_value = one
  !
  do
    if (abs(xx-one) < eps14) return
    if (xx < one) exit
    xx = xx - one
    return_value = return_value*xx
  end do
  !
  ind = 0
  !
  if (xx < half) then
    ind = 1
    return_value = return_value*pi/sin(pi*xx)
    xx = one - xx
  end if
  !
  pr = one
  s = s_start
  !
  do k = 1,11
    pr = pr * (xx-real(k,kind=wp)) / (xx+real(k-1,kind=wp))
    s  = s + gamma_constant(k)*pr
  end do
  !
  g = s * exp(one-xx) * (xx+(nine/two))**(xx-half)
  !
  if (ind == 1) then
    return_value = return_value/g
  else
    return_value = return_value*g
  end if
  !
end function gammafun
!
!###############################################################################
!
elemental function log_gammafun(x) result(return_value)
  real(wp), intent(in) :: x
  real(wp) :: return_value
  return_value = log( gammafun(x) )
end function log_gammafun
!
!###############################################################################
!
