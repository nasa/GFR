elemental function chop(a) result(return_value)
  real(qp), intent(in) :: a
  real(wp) :: return_value
  real(wp) :: a_wp
  real(wp), parameter :: epswp = epsilon(0.0_wp)
  real(qp), parameter :: epswp_qp = real( epswp , kind=qp )
  a_wp = real( a , kind=wp )
  return_value = a_wp * merge( one , zero , abs(a)>epswp_qp )
end function chop
