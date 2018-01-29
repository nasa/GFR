elemental function is_nan_SP(a) result(return_value)
  real(sp), intent(in) :: a
  logical(lk) :: return_value
  return_value = (a /= a)
end function is_nan_SP
elemental function is_nan_DP(a) result(return_value)
  real(dp), intent(in) :: a
  logical(lk) :: return_value
  return_value = (a /= a)
end function is_nan_DP
elemental function is_nan_QP(a) result(return_value)
  real(qp), intent(in) :: a
  logical(lk) :: return_value
  return_value = (a /= a)
end function is_nan_QP
!
!###############################################################################
!
elemental function is_finite_SP(a) result(return_value)
  real(sp), intent(in) :: a
  logical(lk) :: return_value
  return_value = ( a <= huge(a) .and. a >= -huge(a) )
end function is_finite_SP
elemental function is_finite_DP(a) result(return_value)
  real(dp), intent(in) :: a
  logical(lk) :: return_value
  return_value = ( a <= huge(a) .and. a >= -huge(a) )
end function is_finite_DP
elemental function is_finite_QP(a) result(return_value)
  real(qp), intent(in) :: a
  logical(lk) :: return_value
  return_value = ( a <= huge(a) .and. a >= -huge(a) )
end function is_finite_QP
!
!###############################################################################
!
elemental function is_negative_SP(a) result(return_value)
  real(sp), intent(in) :: a
  logical(lk) :: return_value
  return_value = ( a < real(0,kind=kind(a)) )
end function is_negative_SP
elemental function is_negative_DP(a) result(return_value)
  real(dp), intent(in) :: a
  logical(lk) :: return_value
  return_value = ( a < real(0,kind=kind(a)) )
end function is_negative_DP
elemental function is_negative_QP(a) result(return_value)
  real(qp), intent(in) :: a
  logical(lk) :: return_value
  return_value = ( a < real(0,kind=kind(a)) )
end function is_negative_QP
!
!###############################################################################
!
