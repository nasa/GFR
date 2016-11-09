pure function cross_product_SP(a,b,is_2D) result(return_value)
integer, parameter :: lp = SP
real(lp), parameter :: zero = 0.0_lp
real(lp), dimension(3), intent(in) :: a,b
logical(lk),  optional,     intent(in) :: is_2D
real(lp), dimension(3) :: return_value
logical(lk) :: use_3D
use_3D = .true.
if (present(is_2D)) then
  use_3D = .not. is_2D
end if
if (use_3D) then
  return_value(1) = a(2)*b(3) - a(3)*b(2)
  return_value(2) = a(3)*b(1) - a(1)*b(3)
  return_value(3) = a(1)*b(2) - a(2)*b(1)
else
  return_value(1) = a(2) - b(2)
  return_value(2) = b(1) - a(1)
  return_value(3) = zero
end if
end function cross_product_SP
!
pure function cross_product_DP(a,b,is_2D) result(return_value)
integer, parameter :: lp = DP
real(lp), parameter :: zero = 0.0_lp
real(lp), dimension(3), intent(in) :: a,b
logical(lk),  optional,     intent(in) :: is_2D
real(lp), dimension(3) :: return_value
logical(lk) :: use_3D
use_3D = .true.
if (present(is_2D)) then
  use_3D = .not. is_2D
end if
if (use_3D) then
  return_value(1) = a(2)*b(3) - a(3)*b(2)
  return_value(2) = a(3)*b(1) - a(1)*b(3)
  return_value(3) = a(1)*b(2) - a(2)*b(1)
else
  return_value(1) = a(2) - b(2)
  return_value(2) = b(1) - a(1)
  return_value(3) = zero
end if
end function cross_product_DP
!
pure function cross_product_QP(a,b,is_2D) result(return_value)
integer, parameter :: lp = QP
real(lp), parameter :: zero = 0.0_lp
real(lp), dimension(3), intent(in) :: a,b
logical(lk),  optional,     intent(in) :: is_2D
real(lp), dimension(3) :: return_value
logical(lk) :: use_3D
use_3D = .true.
if (present(is_2D)) then
  use_3D = .not. is_2D
end if
if (use_3D) then
  return_value(1) = a(2)*b(3) - a(3)*b(2)
  return_value(2) = a(3)*b(1) - a(1)*b(3)
  return_value(3) = a(1)*b(2) - a(2)*b(1)
else
  return_value(1) = a(2) - b(2)
  return_value(2) = b(1) - a(1)
  return_value(3) = zero
end if
end function cross_product_QP
