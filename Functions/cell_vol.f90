pure function area_tria(xyz_pts) result(return_value)
  !
  ! This function computes the area of a triangular element
  ! from the Jacobian of the coordinates
  !
  !.. Formal Arguments ..
  real(wp), dimension(2,3), intent(in) :: xyz_pts
  !
  !.. Function output ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: x31,y31,x32,y32
  !
continue
  !
  x31 = xyz_pts(1,3) - xyz_pts(1,1)
  y31 = xyz_pts(2,3) - xyz_pts(2,1)
  !
  x32 = xyz_pts(1,3) - xyz_pts(1,2)
  y32 = xyz_pts(2,3) - xyz_pts(2,2)
  !
  return_value = half * (x31*y32 - x32*y31)
  !
end function area_tria
!
!###############################################################################
!
pure function area_quad(xyz_pts) result(return_value)
  !
  ! This function computes the area of a quadrilateral element
  ! from the Jacobian of the coordinates
  !
  !.. Formal Arguments ..
  real(wp), dimension(2,4), intent(in) :: xyz_pts
  !
  !.. Function output ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: x31,y31,x42,y42
  !
continue
  !
  x31 = xyz_pts(1,3) - xyz_pts(1,1)
  y31 = xyz_pts(2,3) - xyz_pts(2,1)
  !
  x42 = xyz_pts(1,4) - xyz_pts(1,2)
  y42 = xyz_pts(2,4) - xyz_pts(2,2)
  !
  return_value = x31*y42 - x42*y31
  !
end function area_quad
!
!###############################################################################
!
pure function vol_tet(xyz_pts) result(return_value)
  !
  ! This function computes the volume of a tetrahedral element
  ! from the Jacobian of the coordinates
  !
  !.. Formal Arguments ..
  real(wp), dimension(3,4), intent(in) :: xyz_pts
  !
  !.. Function output ..
  real(wp) :: return_value
  !
  !.. Local Arrays ..
  real(wp), dimension(3) :: vec21,vec31,vec41
  !
continue
  !
  vec21(:) = xyz_pts(:,2) - xyz_pts(:,1)
  vec31(:) = xyz_pts(:,3) - xyz_pts(:,1)
  vec41(:) = xyz_pts(:,4) - xyz_pts(:,1)
  !
  return_value = dot(vec21,cross_product(vec31,vec41))/six
  !
end function vol_tet
!
!###############################################################################
!
pure function vol_pyr(xyz_pts) result(return_value)
  !
  ! This function computes the volume of a pyramidal element
  ! by dividing the pyramid into two tetrahedrals
  !
  !.. Formal Arguments ..
  real(wp), dimension(3,5), intent(in) :: xyz_pts
  !
  !.. Function output ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: vol1,vol2
  !
continue
  !
  vol1 = vol_tet( xyz_pts(:,[1,2,3,5]) )
  !
  vol2 = vol_tet( xyz_pts(:,[1,3,4,5]) )
  !
  return_value = vol1 + vol2
  !
end function vol_pyr
!
!###############################################################################
!
pure function vol_pri(xyz_pts) result(return_value)
  !
  ! This function computes the volume of a prismatic element
  ! by dividing the prism into one pyramid and one tetrahedral
  !
  !.. Formal Arguments ..
  real(wp), dimension(3,6), intent(in) :: xyz_pts
  !
  !.. Function output ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: vol1,vol2
  !
continue
  !
  vol1 = vol_tet( xyz_pts(:,[1,6,4,5]) )
  !
  vol2 = vol_pyr( xyz_pts(:,[2,5,6,3,1]) )
  !
  return_value = vol1 + vol2
  !
end function vol_pri
!
!###############################################################################
!
pure function vol_hex(xyz_pts) result(return_value)
  !
  ! This function computes the volume of a hexahedral element
  !
  !.. Formal Arguments ..
  real(wp), dimension(3,8), intent(in) :: xyz_pts
  !
  !.. Function output ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: det1,det2,det3
  !
  !.. Local Arrays ..
  real(wp), dimension(3) :: vec11,vec22,vec33
  real(wp), dimension(3) :: vec21,vec31,vec12
  real(wp), dimension(3) :: vec32,vec13,vec23
  !
continue
  !
  vec11(:) = xyz_pts(:,7) - xyz_pts(:,2) + xyz_pts(:,8) - xyz_pts(:,1)
  vec22(:) = xyz_pts(:,7) - xyz_pts(:,4) + xyz_pts(:,6) - xyz_pts(:,1)
  vec33(:) = xyz_pts(:,7) - xyz_pts(:,5) + xyz_pts(:,3) - xyz_pts(:,1)
  !
  vec21(:) = xyz_pts(:,7) - xyz_pts(:,4)
  vec31(:) = xyz_pts(:,3) - xyz_pts(:,1)
  vec12(:) = xyz_pts(:,8) - xyz_pts(:,1)
  !
  vec32(:) = xyz_pts(:,7) - xyz_pts(:,5)
  vec13(:) = xyz_pts(:,7) - xyz_pts(:,2)
  vec23(:) = xyz_pts(:,6) - xyz_pts(:,1)
  !
  det1 = dot(vec11,cross_product(vec21,vec31))
  det2 = dot(vec12,cross_product(vec22,vec32))
  det3 = dot(vec13,cross_product(vec23,vec33))
  !
  return_value = (det1 + det2 + det3) / twelve
  !
end function vol_hex
!
!###############################################################################
!
