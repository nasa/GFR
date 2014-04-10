pure function cross_product(a,b) result(return_value)
real(vbp), dimension(3), intent(in) :: a,b
real(vbp), dimension(3) :: return_value
return_value(1) = a(2)*b(3) - a(3)*b(2)
return_value(2) = a(3)*b(1) - a(1)*b(3)
return_value(3) = a(1)*b(2) - a(2)*b(1)
end function cross_product
