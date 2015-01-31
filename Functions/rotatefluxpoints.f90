pure function RotateFluxPoints_Quad_ab(al,ar,bl,br,face_order) result(return_value)
  !
  !.. Formal Arguments ..
  integer,                intent(in) :: face_order
  real(wp), dimension(:), intent(in) :: al
  real(wp), dimension(:), intent(in) :: ar
  real(wp), dimension(:), intent(in) :: bl
  real(wp), dimension(:), intent(in) :: br
  !
  !.. Function Result ..
  integer :: return_value(1:(face_order+1)**2)
  !
  !.. Local Scalars ..
  integer :: nep,nfp
  !
  !.. Local Arrays ..
  integer, dimension(1:face_order+1,1:face_order+1) :: fp_idx
  !
continue
  !
  nep = face_order+1
  nfp = nep*nep
  !
  ! Create the offset for the flux point indices for the left face
  !
  fp_idx = reshape( intseq(0,nfp-1) , [nep,nep] )
  !
  return_value = 0
  !
  if ( abs(dot(al,ar)-norm2(al)*norm2(ar)) <= eps6 ) then
    !
    if ( abs(dot(bl,br)-norm2(bl)*norm2(br)) <= eps6 ) then
      !
      ! Rotation 1
      !
      return_value = reshape( fp_idx , [nfp] )
      !
    else if ( abs(dot(bl,br)+norm2(bl)*norm2(br)) <= eps6 ) then
      !
      ! Rotation 2
      !
      fp_idx = fp_idx(:,nep:1:-1)
      !
      return_value = reshape( fp_idx , [nfp] )
      !
    end if
    !
  else if ( abs(dot(al,ar)+norm2(al)*norm2(ar)) <= eps6 ) then
    !
    if ( abs(dot(bl,br)-norm2(bl)*norm2(br)) <= eps6 ) then
      !
      ! Rotation 3
      !
      fp_idx = fp_idx(nep:1:-1,:)
      !
      return_value = reshape( fp_idx , [nfp] )
      !
    else if ( abs(dot(bl,br)+norm2(bl)*norm2(br)) <= eps6 ) then
      !
      ! Rotation 4
      !
      fp_idx = fp_idx(nep:1:-1,nep:1:-1)
      !
      return_value = reshape( fp_idx , [nfp] )
      !
    end if
    !
  else if ( abs(dot(al,br)-norm2(al)*norm2(br)) <= eps6 ) then
    !
    if ( abs(dot(bl,ar)-norm2(bl)*norm2(ar)) <= eps6 ) then
      !
      ! Rotation 5
      !
      fp_idx = transpose(fp_idx)
      !
      return_value = reshape( fp_idx , [nfp] )
      !
    else if ( abs(dot(bl,ar)+norm2(bl)*norm2(ar)) <= eps6 ) then
      !
      ! Rotation 6
      !
      fp_idx = transpose(fp_idx(nep:1:-1,:))
      !
      return_value = reshape( fp_idx , [nfp] )
      !
    end if
    !
  else if ( abs(dot(al,br)+norm2(al)*norm2(br)) <= eps6 ) then
    !
    if ( abs(dot(bl,ar)-norm2(bl)*norm2(ar)) <= eps6 ) then
      !
      ! Rotation 7
      !
      fp_idx = transpose(fp_idx(:,nep:1:-1))
      !
      return_value = reshape( fp_idx , [nfp] )
      !
    else if ( abs(dot(bl,ar)+norm2(bl)*norm2(ar)) <= eps6 ) then
      !
      ! Rotation 8
      !
      fp_idx = transpose(fp_idx(nep:1:-1,nep:1:-1))
      !
      return_value = reshape( fp_idx , [nfp] )
      !
    end if
    !
  end if
  !
end function RotateFluxPoints_Quad_ab
!
!###############################################################################
!
pure function RotateFluxPoints_Tria_ab(al,ar,bl,br,nfp) result(return_value)
  !
  !.. Formal Arguments ..
  integer,                intent(in) :: nfp
  real(wp), dimension(:), intent(in) :: al
  real(wp), dimension(:), intent(in) :: ar
  real(wp), dimension(:), intent(in) :: bl
  real(wp), dimension(:), intent(in) :: br
  !
  !.. Function Result ..
  integer :: return_value(1:nfp)
  !
  !.. Local Scalars ..
  integer :: nep
  !
continue
  !
  return_value = 0
  !
end function RotateFluxPoints_Tria_ab
!
!###############################################################################
!
pure function RotateFluxPoints_Quad(left_nodes,right_nodes,face_order) &
                             result(return_value)
  !
  !.. Formal Arguments ..
  integer,               intent(in) :: face_order
  integer, dimension(:), intent(in) :: left_nodes
  integer, dimension(:), intent(in) :: right_nodes
  !
  !.. Function Result ..
  integer, dimension(1:(face_order+1)**2) :: return_value
  !
  !.. Local Scalars ..
  integer :: nfnods,nfp,nep,n,m
  integer :: i1,i2,iint,j1,j2,jint
  !
  !.. Local Arrays ..
  integer, dimension(1:3,1:4) :: idx
  integer, dimension(1:size(right_nodes)) :: rotated_nodes
  integer, dimension(1:face_order+1,1:face_order+1) :: fp_idx
  !
continue
  !
  nep = face_order+1
  nfp = nep*nep
  !
  ! Create the offset for the flux point indices for the left face
  !
  fp_idx = reshape( intseq(0,nfp-1) , [nep,nep] )
  !
  nfnods = size(left_nodes)
  !
  rotated_nodes = right_nodes
  !
  return_value = 0
  !
  transpose_loop: do m = 1,2
    !
    idx(1,1) = 1   ! starting index  for first dimension of fp_idx
    idx(2,1) = nep ! ending   index  for first dimension of fp_idx
    idx(3,1) = 1   ! index increment for first dimension of fp_idx
    !
    idx(1,2) = 1   ! starting index  for second dimension of fp_idx
    idx(2,2) = nep ! ending   index  for second dimension of fp_idx
    idx(3,2) = 1   ! index increment for second dimension of fp_idx
    !
    rotation_loop: do n = 1,nfnods
      !
      if (all(rotated_nodes == left_nodes)) then
        !
        return_value = reshape( fp_idx(idx(1,1):idx(2,1):idx(3,1), &
                                idx(1,2):idx(2,2):idx(3,2)) , [nfp] )
        !
        exit transpose_loop
        !
      end if
      !
      ! Reverse the array ordering of dimension 'n' of fp_idx for the next
      ! iteration of rotation_loop.
      ! NOTE: If n is 3 or 4, we dont care about those values in the idx array.
      !       If n is 3, we just need both dimensions of fp_idx to be reversed
      !       and this has already happened when n was 1 and 2. If n is 4, we
      !       will be leaving this loop anyways and the idx array will be reset
      !       if rotation_loop needs to be done a second time.
      !
      idx(1,n) = nep
      idx(2,n) = 1
      idx(3,n) = -1
      !
      ! Perform a circular shift on the right nodes of the face to see if this
      ! rotation will match the left nodes.
      !
      rotated_nodes = cshift(rotated_nodes,shift=1)
      !
    end do rotation_loop
    !
    ! If the above loop didnt work the first time, reverse the node orderings
    ! and perform the rotation loop again.
    !
    rotated_nodes = right_nodes(nfnods:1:-1)
    !
    ! The flux point indices require a transpose operation since we are having
    ! to reverse the node orderings.
    !
    fp_idx = transpose(fp_idx)
    !
  end do transpose_loop
  !
end function RotateFluxPoints_Quad
!
!###############################################################################
!
subroutine RotateFluxPoints_Quad_sub(left_nodes,left_geom,left_side, &
                                     right_nodes,right_geom,right_side, &
                                     face_order,return_value,rotated)
  !
  !.. Formal Arguments ..
  integer, dimension(:),    intent(in) :: right_nodes
  integer,                  intent(in) :: left_geom
  integer,                  intent(in) :: left_side
  integer, dimension(:),    intent(in) :: left_nodes
  integer,                  intent(in) :: right_geom
  integer,                  intent(in) :: right_side
  integer,                  intent(in) :: face_order
  integer, dimension(:), intent(inout) :: return_value
  integer,               intent(inout) :: rotated
  !
  !.. Local Scalars ..
  integer :: nfnods,nfp,nep,n,m
  integer :: i1,i2,iint,j1,j2,jint
  !
  !.. Local Arrays ..
  integer, dimension(1:3,1:4) :: idx
  integer, dimension(1:4) :: rotated_nodes
  integer, dimension(1:4) :: left_face
  integer, dimension(1:4) :: right_face
  integer, dimension(1:4) :: abcd_left
  integer, dimension(1:4) :: abcd_right
  integer, dimension(1:face_order+1,1:face_order+1) :: fp_idx
  !
  !.. Local Parameters ..
  integer, parameter :: Min_Faces(1:3) = [Xmin3D_Face,Ymin3D_Face,Zmin3D_Face]
  integer, parameter :: Max_Faces(1:3) = [Xmax3D_Face,Ymax3D_Face,Zmax3D_Face]
  !
  integer, parameter :: Xmin_ABCD(1:4) = [4,5,1,8]
  integer, parameter :: Xmax_ABCD(1:4) = [3,6,2,7]
  integer, parameter :: Ymin_ABCD(1:4) = [2,5,1,6]
  integer, parameter :: Ymax_ABCD(1:4) = [3,8,4,7]
  integer, parameter :: Zmin_ABCD(1:4) = [2,4,1,3]
  integer, parameter :: Zmax_ABCD(1:4) = [6,8,5,7]
  !
  integer, parameter :: rotation(1:4,1:8) = reshape( [1,2,3,4, 4,3,2,1, &
                                                      3,4,1,2, 2,1,4,3, &
                                                      2,1,3,4, 4,3,1,2, &
                                                      3,4,2,1, 1,2,4,3] , [4,8] )
  character(len=*), parameter :: cface(1:6) = ["z-","z+","x+","x-","y-","y+"]
  !
continue
  !
  ! ########################################################################
  ! ########################################################################
  ! ########################################################################
  ! ########################################################################
  ! ########################################################################
  ! #####                                                              #####
  ! #####    WARNING!!! WARNING!!! WARNING!!! WARNING!!! WARNING!!!    #####
  ! #####    WARNING!!! WARNING!!! WARNING!!! WARNING!!! WARNING!!!    #####
  ! #####    WARNING!!! WARNING!!! WARNING!!! WARNING!!! WARNING!!!    #####
  ! #####                                                              #####
  ! #####  CURRENTLY THIS SUB ASSUMES THAT THE CELLS ON THE QUAD FACE  #####
  ! #####  ARE BOTH HEXAHEDRAS! THIS IS AN INVALID ASSUMPTION AND      #####
  ! #####  THIS SUB NEEDS TO BE FURTHER GENERALIZED SO THAT ANY CELL   #####
  ! #####  TYPE WITH A QUAD FACE CAN BE USED!                          #####
  ! #####                                                              #####
  ! #####    WARNING!!! WARNING!!! WARNING!!! WARNING!!! WARNING!!!    #####
  ! #####    WARNING!!! WARNING!!! WARNING!!! WARNING!!! WARNING!!!    #####
  ! #####    WARNING!!! WARNING!!! WARNING!!! WARNING!!! WARNING!!!    #####
  ! #####                                                              #####
  ! ########################################################################
  ! ########################################################################
  ! ########################################################################
  ! ########################################################################
  ! ########################################################################
  !
  !
  nep = face_order+1
  nfp = nep*nep
  !
  ! Create the offset for the flux point indices for the left face
  !
  fp_idx = reshape( intseq(0,nfp-1) , [nep,nep] )
  !
  nfnods = size(left_nodes)
  !
  rotated_nodes = right_nodes
  !
  rotated = 0
  !
 !write (161,2) left_nodes
  select case (left_side)
    case (Xmin3D_Face)
      abcd_left = left_nodes(Xmin_ABCD)
   !  write (161,4) Xmin_ABCD
    case (Xmax3D_Face)
      abcd_left = left_nodes(Xmax_ABCD)
   !  write (161,4) Xmax_ABCD
    case (Ymin3D_Face)
      abcd_left = left_nodes(Ymin_ABCD)
   !  write (161,4) Ymin_ABCD
    case (Ymax3D_Face)
      abcd_left = left_nodes(Ymax_ABCD)
   !  write (161,4) Ymax_ABCD
    case (Zmin3D_Face)
      abcd_left = left_nodes(Zmin_ABCD)
   !  write (161,4) Zmin_ABCD
    case (Zmax3D_Face)
      abcd_left = left_nodes(Zmax_ABCD)
   !  write (161,4) Zmax_ABCD
  end select
  !
  if (right_geom == Geom_Quad) then
    !
    abcd_right = right_nodes(Zmin_ABCD)
    !
  else
    !
   !write (161,3) right_nodes
    select case (right_side)
      case (Xmin3D_Face)
      abcd_right = right_nodes(Xmin_ABCD)
   !  write (161,4) Xmin_ABCD
    case (Xmax3D_Face)
      abcd_right = right_nodes(Xmax_ABCD)
   !  write (161,4) Xmax_ABCD
    case (Ymin3D_Face)
      abcd_right = right_nodes(Ymin_ABCD)
   !  write (161,4) Ymin_ABCD
    case (Ymax3D_Face)
      abcd_right = right_nodes(Ymax_ABCD)
   !  write (161,4) Ymax_ABCD
    case (Zmin3D_Face)
      abcd_right = right_nodes(Zmin_ABCD)
   !  write (161,4) Zmin_ABCD
    case (Zmax3D_Face)
      abcd_right = right_nodes(Zmax_ABCD)
   !  write (161,4) Zmax_ABCD
    end select
    !
  end if
  !
  if (all(abcd_left == abcd_right(rotation(:,1)))) then
    !
    return_value = return_value + reshape( fp_idx , [nfp] )
    rotated = 1
    !
  else if (all(abcd_left == abcd_right(rotation(:,2)))) then
    !
    return_value = return_value + reshape( fp_idx(:,nep:1:-1) , [nfp] )
    rotated = 2
    !
  else if (all(abcd_left == abcd_right(rotation(:,3)))) then
    !
    return_value = return_value + reshape( fp_idx(nep:1:-1,:) , [nfp] )
    rotated = 3
    !
  else if (all(abcd_left == abcd_right(rotation(:,4)))) then
    !
   !fp_idx = transpose(fp_idx)
    return_value = return_value + reshape( fp_idx(nep:1:-1,nep:1:-1) , [nfp] )
    rotated = 4
    !
  else if (all(abcd_left == abcd_right(rotation(:,5)))) then
    !
    fp_idx = transpose(fp_idx)
    return_value = return_value + reshape( fp_idx , [nfp] )
    rotated = 5
    !
  else if (all(abcd_left == abcd_right(rotation(:,6)))) then
    !
    fp_idx = transpose(fp_idx)
    return_value = return_value + reshape( fp_idx(:,nep:1:-1) , [nfp] )
    rotated = 6
    !
  else if (all(abcd_left == abcd_right(rotation(:,7)))) then
    !
    fp_idx = transpose(fp_idx)
    return_value = return_value + reshape( fp_idx(nep:1:-1,:) , [nfp] )
    rotated = 7
    !
  else if (all(abcd_left == abcd_right(rotation(:,8)))) then
    !
    fp_idx = transpose(fp_idx)
    return_value = return_value + reshape( fp_idx(nep:1:-1,nep:1:-1) , [nfp] )
    rotated = 8
    !
  end if
 !if (size(right_nodes) == 8) then
 !  write (161,5) cface(left_side),cface(right_side),rotated
 ! !write (161,1) abcd_left
 ! !if (rotated >= 1) then
 ! !  write (161,1) abcd_right
 ! !  write (161,1) (abcd_right(rotation(m,rotated)),m=1,4)
 ! !else
 ! !  write (161,1) abcd_right
 ! !end if
 ! !write (161,*)
 !end if
 !1 format (*(i6))
 !2 format ("Left indices : ",*(i6))
 !3 format ("Right indices: ",*(i6))
 !4 format ("               ",*(i6))
 !5 format (4x,a,4x,a,*(i6))
 !!
 !transpose_loop: do m = 1,2
 !  !
 ! !idx(1,1) = nep ! starting index  for first dimension of fp_idx
 ! !idx(2,1) = 1   ! ending   index  for first dimension of fp_idx
 ! !idx(3,1) = -1  ! index increment for first dimension of fp_idx
 ! !!
 ! !idx(1,2) = nep ! starting index  for second dimension of fp_idx
 ! !idx(2,2) = 1   ! ending   index  for second dimension of fp_idx
 ! !idx(3,2) = -1  ! index increment for second dimension of fp_idx
 ! !!
 !  idx(1,1) = 1   ! starting index  for first dimension of fp_idx
 !  idx(2,1) = nep ! ending   index  for first dimension of fp_idx
 !  idx(3,1) = 1   ! index increment for first dimension of fp_idx
 !  !
 !  idx(1,2) = 1   ! starting index  for second dimension of fp_idx
 !  idx(2,2) = nep ! ending   index  for second dimension of fp_idx
 !  idx(3,2) = 1   ! index increment for second dimension of fp_idx
 !  !
 !  rotation_loop: do n = 1,nfnods
 !    !
 !    if (all(rotated_nodes == left_nodes)) then
 !      !
 !      rotation = (m-1)*nfnods + n
 !      return_value = return_value + reshape( fp_idx(idx(1,1):idx(2,1):idx(3,1), &
 !                                                    idx(1,2):idx(2,2):idx(3,2)) , [nfp] )
 !      !
 !      exit transpose_loop
 !      !
 !    end if
 !    !
 !    ! Reverse the array ordering of dimension 'n' of fp_idx for the next
 !    ! iteration of rotation_loop.
 !    ! NOTE: If n is 3 or 4, we dont care about those values in the idx array.
 !    !       If n is 3, we just need both dimensions of fp_idx to be reversed
 !    !       and this has already happened when n was 1 and 2. If n is 4, we
 !    !       will be leaving this loop anyways and the idx array will be reset
 !    !       if rotation_loop needs to be done a second time.
 !    !
 !   !idx(1,n) = 1
 !   !idx(2,n) = nep
 !   !idx(3,n) = 1
 !   !!
 !    idx(1,n) = nep
 !    idx(2,n) = 1
 !    idx(3,n) = -1
 !    !
 !    ! Perform a circular shift on the right nodes of the face to see if this
 !    ! rotation will match the left nodes.
 !    !
 !    rotated_nodes = cshift(rotated_nodes,shift=1)
 !    !
 !  end do rotation_loop
 !  !
 !  ! If the above loop didnt work the first time, reverse the node orderings
 !  ! and perform the rotation loop again.
 !  !
 !  rotated_nodes = right_nodes(nfnods:1:-1)
 !  !
 !  ! The flux point indices require a transpose operation since we are having
 !  ! to reverse the node orderings.
 !  !
 !  fp_idx = transpose(fp_idx)
 !  !
 !end do transpose_loop
  !
end subroutine RotateFluxPoints_Quad_sub
!
!###############################################################################
!
subroutine RotateFluxPoints_Quad_sub_old(left_nodes,right_nodes,face_order, &
                                     return_value,rotation)
  !
  !.. Formal Arguments ..
  integer,               intent(in) :: face_order
  integer, dimension(:), intent(in) :: left_nodes
  integer, dimension(:), intent(in) :: right_nodes
  integer, dimension(:), intent(inout) :: return_value
  integer,               intent(inout) :: rotation
  !
  !.. Local Scalars ..
  integer :: nfnods,nfp,nep,n,m
  integer :: i1,i2,iint,j1,j2,jint
  !
  !.. Local Arrays ..
  integer, dimension(1:3,1:4) :: idx
  integer, dimension(1:size(right_nodes)) :: rotated_nodes
  integer, dimension(1:face_order+1,1:face_order+1) :: fp_idx
  !
continue
  !
  nep = face_order+1
  nfp = nep*nep
  !
  ! Create the offset for the flux point indices for the left face
  !
  fp_idx = reshape( intseq(0,nfp-1) , [nep,nep] )
  !
  nfnods = size(left_nodes)
  !
  rotated_nodes = right_nodes
  !
  rotation = 0
  !
  transpose_loop: do m = 1,2
    !
   !idx(1,1) = nep ! starting index  for first dimension of fp_idx
   !idx(2,1) = 1   ! ending   index  for first dimension of fp_idx
   !idx(3,1) = -1  ! index increment for first dimension of fp_idx
   !!
   !idx(1,2) = nep ! starting index  for second dimension of fp_idx
   !idx(2,2) = 1   ! ending   index  for second dimension of fp_idx
   !idx(3,2) = -1  ! index increment for second dimension of fp_idx
   !!
    idx(1,1) = 1   ! starting index  for first dimension of fp_idx
    idx(2,1) = nep ! ending   index  for first dimension of fp_idx
    idx(3,1) = 1   ! index increment for first dimension of fp_idx
    !
    idx(1,2) = 1   ! starting index  for second dimension of fp_idx
    idx(2,2) = nep ! ending   index  for second dimension of fp_idx
    idx(3,2) = 1   ! index increment for second dimension of fp_idx
    !
    rotation_loop: do n = 1,nfnods
      !
      if (all(rotated_nodes == left_nodes)) then
        !
        rotation = (m-1)*nfnods + n
        return_value = return_value + reshape( fp_idx(idx(1,1):idx(2,1):idx(3,1), &
                                                      idx(1,2):idx(2,2):idx(3,2)) , [nfp] )
        !
        exit transpose_loop
        !
      end if
      !
      ! Reverse the array ordering of dimension 'n' of fp_idx for the next
      ! iteration of rotation_loop.
      ! NOTE: If n is 3 or 4, we dont care about those values in the idx array.
      !       If n is 3, we just need both dimensions of fp_idx to be reversed
      !       and this has already happened when n was 1 and 2. If n is 4, we
      !       will be leaving this loop anyways and the idx array will be reset
      !       if rotation_loop needs to be done a second time.
      !
     !idx(1,n) = 1
     !idx(2,n) = nep
     !idx(3,n) = 1
     !!
      idx(1,n) = nep
      idx(2,n) = 1
      idx(3,n) = -1
      !
      ! Perform a circular shift on the right nodes of the face to see if this
      ! rotation will match the left nodes.
      !
      rotated_nodes = cshift(rotated_nodes,shift=1)
      !
    end do rotation_loop
    !
    ! If the above loop didnt work the first time, reverse the node orderings
    ! and perform the rotation loop again.
    !
    rotated_nodes = right_nodes(nfnods:1:-1)
    !
    ! The flux point indices require a transpose operation since we are having
    ! to reverse the node orderings.
    !
    fp_idx = transpose(fp_idx)
    !
  end do transpose_loop
  !
end subroutine RotateFluxPoints_Quad_sub_old
!
!###############################################################################
!
pure function RotateFluxPoints_Tria(left_nodes,right_nodes,nfp) &
                             result(return_value)
  !
  !.. Formal Arguments ..
  integer,               intent(in) :: nfp
  integer, dimension(:), intent(in) :: left_nodes
  integer, dimension(:), intent(in) :: right_nodes
  !
  !.. Function Result ..
  integer, dimension(1:nfp) :: return_value
  !
continue
  !
  return_value = 0
  !
end function RotateFluxPoints_Tria
