pure function reciprocal_DP(input_value) result(return_value)
  real(dp), intent(in) :: input_value
  real(dp) :: return_value
  real(dp) :: input_nonzero,input_ratio
  input_nonzero = max(abs(input_value),real(epsilon(0.0_dp),kind=dp))
  input_ratio = input_value / input_nonzero
  return_value = 1.0_dp / input_nonzero*input_ratio
end function reciprocal_DP
!
pure function reciprocal_QP(input_value) result(return_value)
  real(qp), intent(in) :: input_value
  real(qp) :: return_value
  real(qp) :: input_nonzero,input_ratio
  input_nonzero = max(abs(input_value),real(epsilon(0.0_dp),kind=qp))
  input_ratio = input_value / input_nonzero
  return_value = 1.0_qp / input_nonzero*input_ratio
end function reciprocal_QP
