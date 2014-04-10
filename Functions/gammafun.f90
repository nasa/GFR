!
!###############################################################################
!
elemental function gammafun(x) result(return_value)
  !
  real(vbp), intent(in) :: x
  !
  real(vbp) :: return_value
  !
  integer   :: ind,k
  real(vbp) :: xx,pr,s,g
  !
  real(vbp), parameter :: s_start = 0.426401432711220868_vbp
  real(vbp), parameter :: gamma_constant(1:11) = &
             [ -0.524741987629368444e+00_vbp,  0.116154405493589130e+00_vbp, &
               -0.765978624506602380e-02_vbp,  0.899719449391378898e-04_vbp, &
               -0.194536980009534621e-07_vbp,  0.199382839513630987e-10_vbp, &
               -0.204209590209541319e-11_vbp,  0.863896817907000175e-13_vbp, &
                0.152237501608472336e-13_vbp, -0.825725175277719950e-14_vbp, &
                0.299734782205224610e-14_vbp ]
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
    pr = pr * (xx-real(k,vbp)) / (xx+real(k-1,vbp))
    s  = s + gamma_constant(k)*pr
  end do
  !
  g = s * exp(one-xx) * (xx+c4_50)**(xx-half)
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
  real(vbp), intent(in) :: x
  real(vbp) :: return_value
  return_value = log( gammafun(x) )
end function log_gammafun
!
!###############################################################################
!
elemental function binomial(power,term) result(return_value)
  integer, intent(in) :: power,term
  real(vbp) :: return_value
continue
  return_value = log_gammafun( real(power     +1,vbp) ) - &
                 log_gammafun( real(      term+1,vbp) ) - &
                 log_gammafun( real(power-term+1,vbp) )
  return_value = anint( exp(return_value) )
end function binomial
!
!###############################################################################
!
