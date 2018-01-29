pure function compute_jacobian_SP(metrics) result(return_value)
  !.. Local Parameters ..
  integer,  parameter :: lp = SP
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !.. Formal Arguments ..
  real(lp), dimension(:,:), intent(in) :: metrics
  !.. Function Result ..
  real(lp) :: return_value
  !.. Local Scalars ..
  integer :: n1,n2
  !.. Local Arrays ..
  real(lp), dimension(1:3) :: d1,d2,d3
  !
continue
  !
  n1 = size(metrics,dim=1)
  n2 = size(metrics,dim=2)
  !
  ! Initialize d1 and d2 to zero in case the metrics array only
  ! contains x,y information. If there is no z information, the
  ! z component of these arrays will not be overwritten and will
  ! remain zero. These z=zero components will make the vector
  ! returned by the cross_product function below to be in only
  ! the k direction.
  !
  d1 = zero
  d2 = zero
  !
  ! Load the metrics into the local vector arrays
  !
  d1(1:n1) = metrics(1:n1,1)
  d2(1:n1) = metrics(1:n1,2)
  !
  if (n2 == 2) then
    !
    ! Compute the cross product of the (xi,eta) or (r,s) metrics
    !
    d3 = cross_product( d1 , d2 )
    !
    ! For a pure 2D flow, there shouldnt be any z information so
    ! d3 should only contain a k component. If this is the case,
    ! the following norm2 computation will just extract this
    ! k component which is the value of the 2D jacobian. If there
    ! is z information, all of the i,j,k components will contain
    ! information and the jacobian is the magnitude of this vector.
    !
    return_value = norm2( d3 )
    !
  else if (n2 >= 3) then
    !
    ! Load the zeta metrics for a 3D problem
    !
    d3 = metrics(1:n1,3)
    !
    ! Compute the jacobian using the equation
    !   jacobian = d3 .dot. (d1 .cross. d2)
    !
    return_value = dot( d3 , cross_product( d1 , d2 ) )
    !
  end if
  !
end function compute_jacobian_SP
!
!###############################################################################
!
pure function compute_jacobian_DP(metrics) result(return_value)
  !.. Local Parameters ..
  integer,  parameter :: lp = DP
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !.. Formal Arguments ..
  real(lp), dimension(:,:), intent(in) :: metrics
  !.. Function Result ..
  real(lp) :: return_value
  !.. Local Scalars ..
  integer :: n1,n2
  !.. Local Arrays ..
  real(lp), dimension(1:3) :: d1,d2,d3
  !
continue
  !
  n1 = size(metrics,dim=1)
  n2 = size(metrics,dim=2)
  !
  ! Initialize d1 and d2 to zero in case the metrics array only
  ! contains x,y information. If there is no z information, the
  ! z component of these arrays will not be overwritten and will
  ! remain zero. These z=zero components will make the vector
  ! returned by the cross_product function below to be in only
  ! the k direction.
  !
  d1 = zero
  d2 = zero
  !
  ! Load the metrics into the local vector arrays
  !
  d1(1:n1) = metrics(1:n1,1)
  d2(1:n1) = metrics(1:n1,2)
  !
  if (n2 == 2) then
    !
    ! Compute the cross product of the (xi,eta) or (r,s) metrics
    !
    d3 = cross_product( d1 , d2 )
    !
    ! For a pure 2D flow, there shouldnt be any z information so
    ! d3 should only contain a k component. If this is the case,
    ! the following norm2 computation will just extract this
    ! k component which is the value of the 2D jacobian. If there
    ! is z information, all of the i,j,k components will contain
    ! information and the jacobian is the magnitude of this vector.
    !
    return_value = norm2( d3 )
    !
  else if (n2 >= 3) then
    !
    ! Load the zeta metrics for a 3D problem
    !
    d3 = metrics(1:n1,3)
    !
    ! Compute the jacobian using the equation
    !   jacobian = d3 .dot. (d1 .cross. d2)
    !
    return_value = dot( d3 , cross_product( d1 , d2 ) )
    !
  end if
  !
end function compute_jacobian_DP
!
!###############################################################################
!
pure function compute_jacobian_QP(metrics) result(return_value)
  !.. Local Parameters ..
  integer,  parameter :: lp = QP
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !.. Formal Arguments ..
  real(lp), dimension(:,:), intent(in) :: metrics
  !.. Function Result ..
  real(lp) :: return_value
  !.. Local Scalars ..
  integer :: n1,n2
  !.. Local Arrays ..
  real(lp), dimension(1:3) :: d1,d2,d3
  !
continue
  !
  n1 = size(metrics,dim=1)
  n2 = size(metrics,dim=2)
  !
  ! Initialize d1 and d2 to zero in case the metrics array only
  ! contains x,y information. If there is no z information, the
  ! z component of these arrays will not be overwritten and will
  ! remain zero. These z=zero components will make the vector
  ! returned by the cross_product function below to be in only
  ! the k direction.
  !
  d1 = zero
  d2 = zero
  !
  ! Load the metrics into the local vector arrays
  !
  d1(1:n1) = metrics(1:n1,1)
  d2(1:n1) = metrics(1:n1,2)
  !
  if (n2 == 2) then
    !
    ! Compute the cross product of the (xi,eta) or (r,s) metrics
    !
    d3 = cross_product( d1 , d2 )
    !
    ! For a pure 2D flow, there shouldnt be any z information so
    ! d3 should only contain a k component. If this is the case,
    ! the following norm2 computation will just extract this
    ! k component which is the value of the 2D jacobian. If there
    ! is z information, all of the i,j,k components will contain
    ! information and the jacobian is the magnitude of this vector.
    !
    return_value = norm2( d3 )
    !
  else if (n2 >= 3) then
    !
    ! Load the zeta metrics for a 3D problem
    !
    d3 = metrics(1:n1,3)
    !
    ! Compute the jacobian using the equation
    !   jacobian = d3 .dot. (d1 .cross. d2)
    ! NOTE: Must use dot_product here because there is no quad precision
    !       BLAS dot function
    !
    return_value = dot_product( d3 , cross_product( d1 , d2 ) )
    !
  end if
  !
end function compute_jacobian_QP
!
!###############################################################################
!
