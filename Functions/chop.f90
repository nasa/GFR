elemental function chop(a) result(return_value)
  real(qp), intent(in) :: a
  real(wp) :: return_value
  real(wp) :: a_wp
  real(qp), parameter :: epswp = real( rndoff , kind=qp )
  a_wp = real( a , kind=wp )
  return_value = a_wp * merge( one , zero , abs(a)>epswp )
 !return_value = a_wp!* merge( one , zero , abs(a)>epswp )
end function chop
elemental function chop_qp(a) result(return_value)
  real(qp), intent(in) :: a
  real(qp) :: return_value
  real(qp), parameter :: epswp = real( rndoff , kind=qp )
  return_value = a * merge( qone , qzero , abs(a)>epswp )
 !return_value = a!* merge( qone , qzero , abs(a)>epswp )
end function chop_qp
