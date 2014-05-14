elemental function deltafun_DP(logic) result(return_value)
  logical, intent(in) :: logic
  real(dp) :: return_value
  return_value = merge(1.0_dp,0.0_dp,logic)
end function deltafun_DP
!
elemental function deltafun_QP(logic) result(return_value)
  logical, intent(in) :: logic
  real(qp) :: return_value
  return_value = merge(1.0_qp,0.0_qp,logic)
end function deltafun_QP
