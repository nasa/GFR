elemental function reciprocal_SP(input_value) result(return_value)
  real(SP), intent(in) :: input_value
  real(SP) :: return_value
  real(SP) :: input_nonzero,input_ratio
  input_nonzero = max(abs(input_value),real(epsilon(0.0_SP),kind=SP))
  input_ratio = input_value / input_nonzero
  return_value = 1.0_SP / input_nonzero*input_ratio
end function reciprocal_SP
!
elemental function reciprocal_DP(input_value) result(return_value)
  real(DP), intent(in) :: input_value
  real(DP) :: return_value
  real(DP) :: input_nonzero,input_ratio
  input_nonzero = max(abs(input_value),real(epsilon(0.0_DP),kind=DP))
  input_ratio = input_value / input_nonzero
  !
  ! return_value takes the following values based on input_value
  !
  !   iv = input_value
  !   eps = epsilon
  !
  !   if (iv > 0) then
  !                                   iv          1       1
  !                 return_value = --------- = ------ = ----
  !                                |iv|*|iv|    |iv|     iv
  !   else if (iv < 0) then
  !                                  -|iv|       -1       1
  !                 return_value = --------- = ------ = ----
  !                                |iv|*|iv|    |iv|     iv
  !   else if (iv == 0) then
  !                                    iv             0
  !                 return_value = ----------- = ----------- = 0
  !                                |eps|*|eps|   |eps|*|eps|
  !   end if
  !
  ! so
  !
  ! if (input_value /= 0) =>     return_value = 1/input_value
  ! if (input_value == 0) =>     return_value = 0
  !
  return_value = 1.0_DP / input_nonzero*input_ratio
end function reciprocal_DP
!
elemental function reciprocal_QP(input_value) result(return_value)
  real(QP), intent(in) :: input_value
  real(QP) :: return_value
  real(QP) :: input_nonzero,input_ratio
  input_nonzero = max(abs(input_value),real(epsilon(0.0_DP),kind=QP))
  input_ratio = input_value / input_nonzero
  return_value = 1.0_QP / input_nonzero*input_ratio
end function reciprocal_QP
