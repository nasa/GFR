elemental function deltafun_DP_LDK(logic) result(return_value)
  logical, intent(in) :: logic
  real(dp) :: return_value
  return_value = merge(1.0_dp,0.0_dp,logic)
end function deltafun_DP_LDK
!
elemental function deltafun_DP_LK(logic) result(return_value)
  use, intrinsic :: iso_fortran_env, only : logical_kinds
  integer, parameter :: lp = logical_kinds(1)
  logical(lp), intent(in) :: logic
  real(dp) :: return_value
  return_value = merge(1.0_dp,0.0_dp,logic)
end function deltafun_DP_LK
!
elemental function deltafun_QP_LDK(logic) result(return_value)
  logical, intent(in) :: logic
  real(qp) :: return_value
  return_value = merge(1.0_qp,0.0_qp,logic)
end function deltafun_QP_LDK
!
elemental function deltafun_QP_LK(logic) result(return_value)
  use, intrinsic :: iso_fortran_env, only : logical_kinds
  integer, parameter :: lp = logical_kinds(1)
  logical(lp), intent(in) :: logic
  real(qp) :: return_value
  return_value = merge(1.0_qp,0.0_qp,logic)
end function deltafun_QP_LK
