function random_i4(a) result(return_value)
  integer(i4), intent(in) :: a
  integer(i4) :: return_value
  real(sp) :: rnum,b
  if (.not. random_is_initialized) then
    call random_seed()
    random_is_initialized = true
  end if
  call random_number(rnum)
  b = real( sign( max(abs(a),10) , a ) , kind=sp )
  return_value = nint( b * rnum )
end function random_i4
function random_sp(a) result(return_value)
  real(sp), intent(in) :: a
  real(sp) :: return_value
  real(sp) :: rnum
  if (.not. random_is_initialized) then
    call random_seed()
    random_is_initialized = true
  end if
  call random_number(rnum)
  return_value = a * rnum
end function random_sp
function random_dp(a) result(return_value)
  real(dp), intent(in) :: a
  real(dp) :: return_value
  real(dp) :: rnum
  if (.not. random_is_initialized) then
    call random_seed()
    random_is_initialized = true
  end if
  call random_number(rnum)
  return_value = a * rnum
end function random_dp
function random_qp(a) result(return_value)
  real(qp), intent(in) :: a
  real(qp) :: return_value
  real(qp) :: rnum
  if (.not. random_is_initialized) then
    call random_seed()
    random_is_initialized = true
  end if
  call random_number(rnum)
  return_value = a * rnum
end function random_qp
