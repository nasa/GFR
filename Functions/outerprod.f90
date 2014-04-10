function outerprod(a,b) result(return_value)
  !
  real(vbp), dimension(:), intent(in) :: a
  real(vbp), dimension(:), intent(in) :: b
  !
  real(vbp), dimension(size(a),size(b)) :: return_value
  !
  return_value = spread(a,dim=2,ncopies=size(b)) * &
                 spread(b,dim=1,ncopies=size(a))
  !
end function outerprod
