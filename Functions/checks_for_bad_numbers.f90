elemental function is_nan_SP(a) result(return_value)
  use, intrinsic :: iso_fortran_env, only : sp => real32
  real(sp), intent(in) :: a
  logical :: return_value
  return_value = (a /= a)
end function is_nan_SP
elemental function is_nan_DP(a) result(return_value)
  use, intrinsic :: iso_fortran_env, only : dp => real64
  real(dp), intent(in) :: a
  logical :: return_value
  return_value = (a /= a)
end function is_nan_DP
elemental function is_nan_QP(a) result(return_value)
  use, intrinsic :: iso_fortran_env, only : qp => real128
  real(qp), intent(in) :: a
  logical :: return_value
  return_value = (a /= a)
end function is_nan_QP
!
!###############################################################################
!
elemental function is_finite_SP(a) result(return_value)
  use, intrinsic :: iso_fortran_env, only : sp => real32
  real(sp), intent(in) :: a
  logical :: return_value
  return_value = ( a <= huge(a) .and. a >= -huge(a) )
end function is_finite_SP
elemental function is_finite_DP(a) result(return_value)
  use, intrinsic :: iso_fortran_env, only : dp => real64
  real(dp), intent(in) :: a
  logical :: return_value
  return_value = ( a <= huge(a) .and. a >= -huge(a) )
end function is_finite_DP
elemental function is_finite_QP(a) result(return_value)
  use, intrinsic :: iso_fortran_env, only : qp => real128
  real(qp), intent(in) :: a
  logical :: return_value
  return_value = ( a <= huge(a) .and. a >= -huge(a) )
end function is_finite_QP
!
!###############################################################################
!
elemental function is_negative_SP(a) result(return_value)
  use, intrinsic :: iso_fortran_env, only : sp => real32
  real(sp), intent(in) :: a
  logical :: return_value
  return_value = ( a < real(0,kind=kind(a)) )
end function is_negative_SP
elemental function is_negative_DP(a) result(return_value)
  use, intrinsic :: iso_fortran_env, only : dp => real64
  real(dp), intent(in) :: a
  logical :: return_value
  return_value = ( a < real(0,kind=kind(a)) )
end function is_negative_DP
elemental function is_negative_QP(a) result(return_value)
  use, intrinsic :: iso_fortran_env, only : qp => real128
  real(qp), intent(in) :: a
  logical :: return_value
  return_value = ( a < real(0,kind=kind(a)) )
end function is_negative_QP
!
!###############################################################################
!
