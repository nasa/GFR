pure function face_nodes_quad(face_loc,cell_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer,               intent(in) :: face_loc
  integer, dimension(:), intent(in) :: cell_nodes
  !
  !.. Function Return Value
  integer, dimension(1:2) :: return_value
  !
  !.. Local Parameters ..
  integer, parameter :: fnodes(1:2,1:4) = reshape( [ 1, 2, &
                                                     2, 3, &
                                                     3, 4, &
                                                     4, 1 ] , [2,4] )
  !
continue
  !
  return_value = cell_nodes(fnodes(:,face_loc))
  !
end function face_nodes_quad
!
!###############################################################################
!
pure function face_nodes_tria(face_loc,cell_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer,               intent(in) :: face_loc
  integer, dimension(:), intent(in) :: cell_nodes
  !
  !.. Function Return Value
  integer, dimension(1:2) :: return_value
  !
  !.. Local Parameters ..
  integer, parameter :: fnodes(1:2,1:3) = reshape( [ 1, 2, &
                                                     2, 3, &
                                                     3, 1 ] , [2,3] )
  !
continue
  !
  return_value = cell_nodes(fnodes(:,face_loc))
  !
end function face_nodes_tria
!
!###############################################################################
!
pure function face_nodes_tetr(face_loc,cell_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer,               intent(in) :: face_loc
  integer, dimension(:), intent(in) :: cell_nodes
  !
  !.. Function Return Value
  integer, dimension(1:3) :: return_value
  !
  !.. Local Parameters ..
  integer, parameter :: fnodes(1:3,1:4) = reshape( [ 2, 3, 4, &
                                                     3, 1, 4, &
                                                     4, 1, 2, &
                                                     1, 3, 2 ] , [3,4] )
  !
continue
  !
  return_value = cell_nodes(fnodes(:,face_loc))
  !
  ! Tetrahedron elements are not currently set up so set
  ! return_value to the currently not set up error value of -4
  !
  return_value = -4
  !
end function face_nodes_tetr
!
!###############################################################################
!
pure function face_nodes_pyra(face_loc,cell_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer,               intent(in) :: face_loc
  integer, dimension(:), intent(in) :: cell_nodes
  !
  !.. Function Return Value
  integer, dimension(1:4) :: return_value
  !
  !.. Local Parameters ..
  integer, parameter :: fnodes(1:4,1:5) = reshape( [ 2, 3, 5, 0, &
                                                     3, 4, 5, 0, &
                                                     4, 1, 5, 0, &
                                                     1, 2, 5, 0, &
                                                     1, 4, 3, 2 ] , [4,5] )
  !
continue
  !
  return_value = cell_nodes(fnodes(:,face_loc))
  !
  ! Pyramid elements are not currently set up so set
  ! return_value to the currently not set up error value of -5
  !
  return_value = -5
  !
end function face_nodes_pyra
!
!###############################################################################
!
pure function face_nodes_pris(face_loc,cell_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer,               intent(in) :: face_loc
  integer, dimension(:), intent(in) :: cell_nodes
  !
  !.. Function Return Value
  integer, dimension(1:4) :: return_value
  !
  !.. Local Parameters ..
  integer, parameter :: fnodes(1:4,1:5) = reshape( [ 1, 3, 2, 0, &
                                                     4, 5, 6, 0, &
                                                     1, 2, 5, 4, &
                                                     2, 3, 6, 5, &
                                                     4, 6, 3, 1 ] , [4,5] )
  !
continue
  !
  return_value = cell_nodes(fnodes(:,face_loc))
  !
  ! Prism elements are not currently set up so set
  ! return_value to the currently not set up error value of -6
  !
  return_value = -6
  !
end function face_nodes_pris
!
!###############################################################################
!
pure function face_nodes_hexa(face_loc,cell_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer,               intent(in) :: face_loc
  integer, dimension(:), intent(in) :: cell_nodes
  !
  !.. Function Return Value
  integer, dimension(1:4) :: return_value
  !
  !.. Local Parameters ..
  integer, parameter :: fnodes(1:4,1:6) = reshape( [ 1, 4, 3, 2, &
                                                     5, 6, 7, 8, &
                                                     2, 3, 7, 6, &
                                                     5, 8, 4, 1, &
                                                     1, 2, 6, 5, &
                                                     3, 4, 8, 7 ] , [4,6] )
  !
continue
  !
  return_value = cell_nodes(fnodes(:,face_loc))
  !
  ! Hexahedron elements are not currently set up so set
  ! return_value to the currently not set up error value of -8
  !
  return_value = -8
  !
end function face_nodes_hexa
!
!###############################################################################
!
