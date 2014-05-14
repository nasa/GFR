function outerprod_DP(a,b) result(return_value)
  real(dp), dimension(:), intent(in) :: a
  real(dp), dimension(:), intent(in) :: b
  real(dp), dimension(size(a),size(b)) :: return_value
  return_value = spread(a,dim=2,ncopies=size(b)) * &
                 spread(b,dim=1,ncopies=size(a))
end function outerprod_DP
!
function outerprod_QP(a,b) result(return_value)
  real(qp), dimension(:), intent(in) :: a
  real(qp), dimension(:), intent(in) :: b
  real(qp), dimension(size(a),size(b)) :: return_value
  return_value = spread(a,dim=2,ncopies=size(b)) * &
                 spread(b,dim=1,ncopies=size(a))
end function outerprod_QP
