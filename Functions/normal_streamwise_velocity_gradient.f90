pure function normal_streamwise_velocity_gradient_cv(cv,dcvdx,rnrm) &
                                              result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:),   intent(in) :: cv
  real(wp), dimension(:,:), intent(in) :: dcvdx
  real(wp), dimension(:),   intent(in) :: rnrm
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,m,nr
  real(wp) :: rho_inv,vnrm,vtan_inv
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(rnrm)) :: unrm
  real(wp), dimension(1:size(rnrm)) :: vel
  real(wp), dimension(1:size(rnrm)) :: vcoef
  real(wp), dimension(1:size(rnrm)) :: dvtdx
  real(wp), dimension(1:size(rnrm),1:size(rnrm)) :: dpvdx
  !
continue
  !
  ! Get some quick problem size parameters
  !
  nr = size(rnrm) ! total number of dimensions
  !
  ! Compute the inverse of density
  !
  rho_inv = one / cv(nec)
  !
  ! First make sure that the normal is the unit normal
  !
  unrm(1:nr) = unit_vector( rnrm(1:nr) )
  !
  ! Get the velocities from the conserved variables
  !
  vel(1:nr) = cv(nmb:nme) * rho_inv
  !
  ! Compute the normal velocity
  !
  vnrm = dot( vel , unrm )
  !
  ! Compute the reciprocal of the streamwise/tangent velocity
  ! NOTE: For a given vector incident to a face, a right triangle can
  !       be formed by 1) the part of the vector normal to the face,
  !       2) the part of the vector tangent to the face, 3) and the
  !       hypotenuese is the vector itself. Thus, using the pythagorean
  !       theorem of vector^2 = normal^2 + tangent^2, the streamwise/
  !       tangent velocity can be computed as
  !     streamwise/tangent velocity = sqrt( (magnitude of vector)^2 +
  !                                         (normal component of vector)^2 )
  !
  vtan_inv = reciprocal( sqrt( selfdot(vel) - vnrm*vnrm ) )
  !
  ! Compute the directional components of the tangential velocity
  !
  vcoef(1:nr) = vel(1:nr) - vnrm*unrm(1:nr)
  !
  ! Get the gradients of the velocities from the gradients of
  ! the momentum and the conservative variables
  !
  do m = 1,nr
    do l = 1,nr
      dpvdx(l,m) = rho_inv*dcvdx(l,nec+m) - rho_inv*vel(m)*dcvdx(l,nec)
    end do
  end do
  !
  ! Compute the directional derivatives of the streamwise velocity
  !
  do l = 1,nr
    dvtdx(l) = vtan_inv * dot( dpvdx(l,1:nr) , vcoef(1:nr) )
  end do
  !
  ! Finally get the normal of the streamwise velocity gradient
  !
  return_value = dot( dvtdx , unrm )
  !
end function normal_streamwise_velocity_gradient_cv
!
!###############################################################################
!
pure function normal_streamwise_velocity_gradient_pv(pv,dpvdx,rnrm) &
                                              result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), dimension(:),   intent(in) :: pv
  real(wp), dimension(:,:), intent(in) :: dpvdx
  real(wp), dimension(:),   intent(in) :: rnrm
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  integer :: l,m,nr
  real(wp) :: vnrm,vtan_inv
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(rnrm)) :: unrm
  real(wp), dimension(1:size(rnrm)) :: vcoef
  real(wp), dimension(1:size(rnrm)) :: dvtdx
  !
continue
  !
  ! Get some quick problem size parameters
  !
  nr = size(rnrm) ! total number of dimensions
  !
  ! First make sure that the normal is the unit normal
  !
  unrm(1:nr) = unit_vector( rnrm(1:nr) )
  !
  ! Compute the normal velocity
  !
  vnrm = dot( pv(nmb:nme) , unrm )
  !
  ! Compute the reciprocal of the streamwise/tangent velocity
  ! NOTE: For a given vector incident to a face, a right triangle can
  !       be formed by 1) the part of the vector normal to the face,
  !       2) the part of the vector tangent to the face, 3) and the
  !       hypotenuese is the vector itself. Thus, using the pythagorean
  !       theorem of vector^2 = normal^2 + tangent^2, the streamwise/
  !       tangent velocity can be computed as
  !     streamwise/tangent velocity = sqrt( (magnitude of vector)^2 +
  !                                         (normal component of vector)^2 )
  !
  vtan_inv = reciprocal( sqrt( selfdot(pv(nmb:nme)) - vnrm*vnrm ) )
  !
  ! Compute the directional components of the tangential velocity
  !
  vcoef(1:nr) = pv(nmb:nme) - vnrm*unrm(1:nr)
  !
  ! Compute the directional derivatives of the streamwise velocity
  !
  do l = 1,nr
    dvtdx(l) = vtan_inv * dot( dpvdx(l,nmb:nme) , vcoef(1:nr) )
  end do
  !
  ! Finally get the normal of the streamwise velocity gradient
  !
  return_value = dot( dvtdx , unrm )
  !
end function normal_streamwise_velocity_gradient_pv
!
!###############################################################################
!
