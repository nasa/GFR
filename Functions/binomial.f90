elemental function binomial(power,term) result(return_value)
  integer, intent(in) :: power,term
  real(wp) :: return_value
continue
  return_value = log_gamma( real(power     +1,kind=wp) ) - &
                 log_gamma( real(      term+1,kind=wp) ) - &
                 log_gamma( real(power-term+1,kind=wp) )
  return_value = anint( exp(return_value) )
end function binomial
!
!###############################################################################
!
