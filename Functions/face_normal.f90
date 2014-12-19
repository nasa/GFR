pure function face_normal_SP(face_nodes) result(return_value)
  !.. Local Parameters..
  integer,  parameter :: lp = SP
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: half = 0.5_lp
  !.. Formal Arguments ..
  real(lp), dimension(:,:), intent(in) :: face_nodes
  !.. Function Return Value ..
  real(lp), dimension(1:size(face_nodes,dim=1)) :: return_value
  !.. Local Scalars ..
  integer  :: ndim,npts
  real(lp) :: coef
  !.. Local Arrays ..
  real(lp), dimension(1:3) :: d1,d2,d3
  !
continue
  !
  ndim = size(face_nodes,dim=1)
  npts = size(face_nodes,dim=2)
  !
  d1 = one ! Initialize d1 and d2 to 1 in case face_nodes
  d2 = one ! only contains x,y coordinates
  !
  select case (npts)
    case (2)
      coef = one
      ! I think d1 needs to be point 2 and d2 needs to be point 1
      d1(1:ndim) = face_nodes(1:ndim,2)
      d2(1:ndim) = face_nodes(1:ndim,1)
      ! If that above order is incorrect, try this order
     !d1(1:ndim) = face_nodes(1:ndim,1)
     !d2(1:ndim) = face_nodes(1:ndim,2)
    case (3,4)
      coef = half
      d1 = face_nodes(:,npts-1) - face_nodes(:,     1)
      d2 = face_nodes(:,npts  ) - face_nodes(:,npts-2)
    case default
      return_value = one
      return
  end select
  !
  d3 = cross_product( d1 , d2 )
  !
  return_value(1:ndim) = coef * d3(1:ndim)
  !
end function face_normal_SP
!
!###############################################################################
!
pure function face_normal_DP(face_nodes) result(return_value)
  !.. Local Parameters..
  integer,  parameter :: lp = DP
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: half = 0.5_lp
  !.. Formal Arguments ..
  real(lp), dimension(:,:), intent(in) :: face_nodes
  !.. Function Return Value ..
  real(lp), dimension(1:size(face_nodes,dim=1)) :: return_value
  !.. Local Scalars ..
  integer  :: ndim,npts
  real(lp) :: coef
  !.. Local Arrays ..
  real(lp), dimension(1:3) :: d1,d2,d3
  !
continue
  !
  ndim = size(face_nodes,dim=1)
  npts = size(face_nodes,dim=2)
  !
  d1 = one ! Initialize d1 and d2 to 1 in case face_nodes
  d2 = one ! only contains x,y coordinates
  !
  select case (npts)
    case (2)
      coef = one
      ! I think d1 needs to be point 2 and d2 needs to be point 1
      d1(1:ndim) = face_nodes(1:ndim,2)
      d2(1:ndim) = face_nodes(1:ndim,1)
      ! If that above order is incorrect, try this order
     !d1(1:ndim) = face_nodes(1:ndim,1)
     !d2(1:ndim) = face_nodes(1:ndim,2)
    case (3,4)
      coef = half
      d1 = face_nodes(:,npts-1) - face_nodes(:,     1)
      d2 = face_nodes(:,npts  ) - face_nodes(:,npts-2)
    case default
      return_value = one
      return
  end select
  !
  d3 = cross_product( d1 , d2 )
  !
  return_value(1:ndim) = coef * d3(1:ndim)
  !
end function face_normal_DP
!
!###############################################################################
!
pure function face_normal_QP(face_nodes) result(return_value)
  !.. Local Parameters..
  integer,  parameter :: lp = QP
  real(lp), parameter :: one  = 1.0_lp
  real(lp), parameter :: half = 0.5_lp
  !.. Formal Arguments ..
  real(lp), dimension(:,:), intent(in) :: face_nodes
  !.. Function Return Value ..
  real(lp), dimension(1:size(face_nodes,dim=1)) :: return_value
  !.. Local Scalars ..
  integer  :: ndim,npts
  real(lp) :: coef
  !.. Local Arrays ..
  real(lp), dimension(1:3) :: d1,d2,d3
  !
continue
  !
  ndim = size(face_nodes,dim=1)
  npts = size(face_nodes,dim=2)
  !
  d1 = one ! Initialize d1 and d2 to 1 in case face_nodes
  d2 = one ! only contains x,y coordinates
  !
  select case (npts)
    case (2)
      coef = one
      ! I think d1 needs to be point 2 and d2 needs to be point 1
      d1(1:ndim) = face_nodes(1:ndim,2)
      d2(1:ndim) = face_nodes(1:ndim,1)
      ! If that above order is incorrect, try this order
     !d1(1:ndim) = face_nodes(1:ndim,1)
     !d2(1:ndim) = face_nodes(1:ndim,2)
    case (3,4)
      coef = half
      d1 = face_nodes(:,npts-1) - face_nodes(:,     1)
      d2 = face_nodes(:,npts  ) - face_nodes(:,npts-2)
    case default
      return_value = one
      return
  end select
  !
  d3 = cross_product( d1 , d2 )
  !
  return_value(1:ndim) = coef * d3(1:ndim)
  !
end function face_normal_QP
!
!###############################################################################
!
