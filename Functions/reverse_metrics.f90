pure function reverse_metrics_SP(metrics) result(return_value)
  !.. Local Parameters ..
  integer,  parameter :: lp = SP
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !.. Formal Arguments ..
  real(lp), dimension(:,:), intent(in) :: metrics
  !.. Function Result ..
  real(lp) :: return_value(1:size(metrics,dim=1),1:size(metrics,dim=2))
  !.. Local Scalars ..
  integer :: n1,n2
  !.. Local Arrays ..
  real(lp), dimension(1:3,1:3) :: met
  real(lp), dimension(1:3,1:3) :: rev_met
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
  met = zero
  met(:,3) = one
  !
  met(1:n1,1:n2) = metrics
  !
  ! Load the metrics into the local vector arrays
  !
  rev_met(:,1) = cross_product( met(:,2) , met(:,3) )
  rev_met(:,2) = cross_product( met(:,3) , met(:,1) )
  rev_met(:,3) = cross_product( met(:,1) , met(:,2) )
  !
  return_value = rev_met(1:n1,1:n2)
  !
end function reverse_metrics_SP
!
!###############################################################################
!
pure function reverse_metrics_DP(metrics) result(return_value)
  !.. Local Parameters ..
  integer,  parameter :: lp = DP
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !.. Formal Arguments ..
  real(lp), dimension(:,:), intent(in) :: metrics
  !.. Function Result ..
  real(lp) :: return_value(1:size(metrics,dim=1),1:size(metrics,dim=2))
  !.. Local Scalars ..
  integer :: n1,n2
  !.. Local Arrays ..
  real(lp), dimension(1:3,1:3) :: met
  real(lp), dimension(1:3,1:3) :: rev_met
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
  met = zero
  met(:,3) = one
  !
  met(1:n1,1:n2) = metrics
  !
  ! Load the metrics into the local vector arrays
  !
  rev_met(:,1) = cross_product( met(:,2) , met(:,3) )
  rev_met(:,2) = cross_product( met(:,3) , met(:,1) )
  rev_met(:,3) = cross_product( met(:,1) , met(:,2) )
  !
  return_value = rev_met(1:n1,1:n2)
  !
end function reverse_metrics_DP
!
!###############################################################################
!
pure function reverse_metrics_QP(metrics) result(return_value)
  !.. Local Parameters ..
  integer,  parameter :: lp = QP
  real(lp), parameter :: zero = 0.0_lp
  real(lp), parameter :: one  = 1.0_lp
  !.. Formal Arguments ..
  real(lp), dimension(:,:), intent(in) :: metrics
  !.. Function Result ..
  real(lp) :: return_value(1:size(metrics,dim=1),1:size(metrics,dim=2))
  !.. Local Scalars ..
  integer :: n1,n2
  !.. Local Arrays ..
  real(lp), dimension(1:3,1:3) :: met
  real(lp), dimension(1:3,1:3) :: rev_met
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
  met = zero
  met(:,3) = one
  !
  met(1:n1,1:n2) = metrics
  !
  ! Load the metrics into the local vector arrays
  !
  rev_met(:,1) = cross_product( met(:,2) , met(:,3) )
  rev_met(:,2) = cross_product( met(:,3) , met(:,1) )
  rev_met(:,3) = cross_product( met(:,1) , met(:,2) )
  !
  return_value = rev_met(1:n1,1:n2)
  !
end function reverse_metrics_QP
!
!###############################################################################
!
