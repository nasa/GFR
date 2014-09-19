pure function host_face_quad(int_nodes,bnd_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: int_nodes
  integer, dimension(:), intent(in) :: bnd_nodes
  !
  !.. Function Return Value
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Arrays ..
  integer, dimension(1:4) :: marker
  !
  !.. Local Parameters ..
  integer, parameter :: num_face_nodes = 2
  integer, parameter :: fnodes(1:2,1:4) = reshape( [ 1, 2, &
                                                     2, 3, &
                                                     3, 4, &
                                                     4, 1 ] , [2,4] )
  !
continue
  !
  ! Initialize return_value with the face not found value of -99
  !
  return_value = -99
  !
  ! Initialize the marker array to zero
  !
  marker(:) = 0
  !
  ! Mark each host cell node that is on the boundary
  !
  do i = 1,num_face_nodes
    where (int_nodes(1:4) == bnd_nodes(i)) marker = 1
  end do
  !
  ! Now find which face corresponds to the marked nodes
  !
  do i = 1,size(fnodes,dim=2)
    if ( sum( marker(fnodes(:,i)) ) == num_face_nodes ) then
      return_value = i
      exit
    end if
  end do
  !
end function host_face_quad
!
!###############################################################################
!
pure function host_face_tria(int_nodes,bnd_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: int_nodes
  integer, dimension(:), intent(in) :: bnd_nodes
  !
  !.. Function Return Value
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Arrays ..
  integer, dimension(1:3) :: marker
  !
  !.. Local Parameters ..
  integer, parameter :: num_face_nodes = 2
  integer, parameter :: fnodes(1:2,1:3) = reshape( [ 1, 2, &
                                                     2, 3, &
                                                     3, 1 ] , [2,3] )
  !
continue
  !
  ! Initialize return_value with the face not found value of -99
  !
  return_value = -99
  !
  ! Initialize the marker array to zero
  !
  marker(:) = 0
  !
  ! Mark each host cell node that is on the boundary
  !
  do i = 1,num_face_nodes
    where (int_nodes(1:3) == bnd_nodes(i)) marker = 1
  end do
  !
  ! Now find which face corresponds to the marked nodes
  !
  do i = 1,size(fnodes,dim=2)
    if ( sum( marker(fnodes(:,i)) ) == num_face_nodes ) then
      return_value = i
      exit
    end if
  end do
  !
end function host_face_tria
!
!###############################################################################
!
pure function host_face_tetr(int_nodes,bnd_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: int_nodes
  integer, dimension(:), intent(in) :: bnd_nodes
  !
  !.. Function Return Value
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Arrays ..
  integer, dimension(1:4) :: marker
  !
  !.. Local Parameters ..
  integer, parameter :: num_face_nodes = 3
  integer, parameter :: fnodes(1:3,1:4) = reshape( [ 2, 3, 4, &
                                                     3, 1, 4, &
                                                     4, 1, 2, &
                                                     1, 3, 2 ] , [3,4] )
  !
continue
  !
  ! Initialize return_value with the face not found value of -99
  !
  return_value = -99
  !
  ! Initialize the marker array to zero
  !
  marker(:) = 0
  !
  ! Mark each host cell node that is on the boundary
  !
  do i = 1,num_face_nodes
    where (int_nodes(1:4) == bnd_nodes(i)) marker = 1
  end do
  !
  ! Now find which face corresponds to the marked nodes
  !
  do i = 1,size(fnodes,dim=2)
    if ( sum( marker(fnodes(:,i)) ) == num_face_nodes ) then
      return_value = i
      exit
    end if
  end do
  !
  ! Tetrahedron elements are not currently set up so set
  ! return_value to the currently not set up error value of -4
  !
  return_value = -4
  !
end function host_face_tetr
!
!###############################################################################
!
pure function host_face_pyra(int_nodes,bnd_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: int_nodes
  integer, dimension(:), intent(in) :: bnd_nodes
  !
  !.. Function Return Value
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Arrays ..
  integer, dimension(1:5) :: marker
  !
  !.. Local Parameters ..
  integer, parameter :: num_face_nodes = 3
  integer, parameter :: fnodes(1:4,1:5) = reshape( [ 2, 3, 5, 0, &
                                                     3, 4, 5, 0, &
                                                     4, 1, 5, 0, &
                                                     1, 2, 5, 0, &
                                                     1, 4, 3, 2 ] , [4,5] )
  !
continue
  !
  ! Initialize return_value with the face not found value of -99
  !
  return_value = -99
  !
  ! Initialize the marker array to zero
  !
  marker(:) = 0
  !
  ! Mark each host cell node that is on the boundary
  ! NOTE: We only check the first three boundary nodes because those three nodes
  !       are enough to identify which of the five faces is on the boundary.
  !       This way we dont have to differentiate between the quad and triangle
  !       faces of the pyramid element.
  !
  do i = 1,num_face_nodes
    where (int_nodes(1:5) == bnd_nodes(i)) marker = 1
  end do
  !
  ! Now find which face corresponds to the marked nodes
  !
  do i = 1,size(fnodes,dim=2)
    if ( sum( marker(fnodes(:,i)) ) == num_face_nodes ) then
      return_value = i
      exit
    end if
  end do
  !
  ! Pyramid elements are not currently set up so set
  ! return_value to the currently not set up error value of -5
  !
  return_value = -5
  !
end function host_face_pyra
!
!###############################################################################
!
pure function host_face_pris(int_nodes,bnd_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: int_nodes
  integer, dimension(:), intent(in) :: bnd_nodes
  !
  !.. Function Return Value
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Arrays ..
  integer, dimension(1:6) :: marker
  !
  !.. Local Parameters ..
  integer, parameter :: num_face_nodes = 3
  integer, parameter :: fnodes(1:4,1:5) = reshape( [ 1, 3, 2, 0, &
                                                     4, 5, 6, 0, &
                                                     1, 2, 5, 4, &
                                                     2, 3, 6, 5, &
                                                     4, 6, 3, 1 ] , [4,5] )
  !
continue
  !
  ! Initialize return_value with the face not found value of -99
  !
  return_value = -99
  !
  ! Initialize the marker array to zero
  !
  marker(:) = 0
  !
  ! Mark each host cell node that is on the boundary
  ! NOTE: We only check the first three boundary nodes because those three nodes
  !       are enough to identify which of the five faces is on the boundary.
  !       This way we dont have to differentiate between the quad and triangle
  !       faces of the prism element.
  !
  do i = 1,num_face_nodes
    where (int_nodes(1:6) == bnd_nodes(i)) marker = 1
  end do
  !
  ! Now find which face corresponds to the marked nodes
  !
  do i = 1,size(fnodes,dim=2)
    if ( sum( marker(fnodes(:,i)) ) == num_face_nodes ) then
      return_value = i
      exit
    end if
  end do
  !
  ! Prism elements are not currently set up so set
  ! return_value to the currently not set up error value of -6
  !
  return_value = -6
  !
end function host_face_pris
!
!###############################################################################
!
pure function host_face_hexa(int_nodes,bnd_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: int_nodes
  integer, dimension(:), intent(in) :: bnd_nodes
  !
  !.. Function Return Value
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Arrays ..
  integer, dimension(1:8) :: marker
  !
  !.. Local Parameters ..
  integer, parameter :: num_face_nodes = 4
  integer, parameter :: fnodes(1:4,1:6) = reshape( [ 1, 4, 3, 2, &
                                                     5, 6, 7, 8, &
                                                     2, 3, 7, 6, &
                                                     5, 8, 4, 1, &
                                                     1, 2, 6, 5, &
                                                     3, 4, 8, 7 ] , [4,6] )
  !
continue
  !
  ! Initialize return_value with the face not found value of -99
  !
  return_value = -99
  !
  ! Initialize the marker array to zero
  !
  marker(:) = 0
  !
  ! Mark each host cell node that is on the boundary
  !
  do i = 1,num_face_nodes
    where (int_nodes(1:8) == bnd_nodes(i)) marker = 1
  end do
  !
  ! Now find which face corresponds to the marked nodes
  !
  do i = 1,size(fnodes,dim=2)
    if ( sum( marker(fnodes(:,i)) ) == num_face_nodes ) then
      return_value = i
      exit
    end if
  end do
  !
  ! Hexahedron elements are not currently set up so set
  ! return_value to the currently not set up error value of -8
  !
  return_value = -8
  !
end function host_face_hexa
!
!###############################################################################
!
