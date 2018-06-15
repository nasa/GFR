elemental function chop_WP(a) result(return_value)
  integer, parameter :: lp = WP
  real(QP), intent(in) :: a
  real(lp) :: return_value
  real(QP), parameter :: eps = real( rndoff , kind=QP )
  return_value = real( merge( a , 0.0_QP , abs(a)>eps ) , kind=lp )
 !return_value = real( a , kind=lp )
end function chop_WP
!
elemental function chop_SP(a,b) result(return_value)
  integer, parameter :: lp = SP
  real(QP), intent(in) :: a
  real(lp), intent(in) :: b
  real(lp) :: return_value
  real(QP), parameter :: eps = real( rndoff , kind=QP )
  return_value = real( merge( a , 0.0_QP , abs(a)>eps ) , kind=lp )
 !return_value = real( a , kind=lp )
end function chop_SP
!
elemental function chop_DP(a,b) result(return_value)
  integer, parameter :: lp = DP
  real(QP), intent(in) :: a
  real(lp), intent(in) :: b
  real(lp) :: return_value
  real(QP), parameter :: eps = real( rndoff , kind=QP )
  return_value = real( merge( a , 0.0_QP , abs(a)>eps ) , kind=lp )
 !return_value = real( a , kind=lp )
end function chop_DP
!
elemental function chop_QP(a,b) result(return_value)
  integer, parameter :: lp = QP
  real(QP), intent(in) :: a
  real(lp), intent(in) :: b
  real(lp) :: return_value
  real(QP), parameter :: eps = real( rndoff , kind=QP )
  return_value = real( merge( a , 0.0_QP , abs(a)>eps ) , kind=lp )
 !return_value = real( a , kind=lp )
end function chop_QP
!
