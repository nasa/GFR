module bc_periodic_mod
  !
  !.. Use Statements
  use module_kind_types
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nmb,nme,nee
  use eqn_idx_mod, only : ncs
  !
  implicit none
  !
  private
  !
  public :: get_bc_periodic_solution
  public :: get_bc_periodic_gradient
  !
contains
!
!###############################################################################
!
!pure &
subroutine get_bc_periodic_solution(n1,n2,face_cv,periodic_list, &
                                    rotation_matrix)
  !
  ! #############################################################
  ! #############################################################
  ! ##  Periodic Boundary Condition
  ! #############################################################
  ! #############################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved variables at the flux points
  ! PERIODIC_LIST : For each face point n from n1 to n2, periodic_list(n) gives
  !                 the face point p that matches/corresponds to n based on
  !                 the periodic connectivity. Once this subroutine is finished,
  !                 the following should be true:
  !                   face_cv(n,:,2) = face_cv(periodic_list(n-n1+1),:,1)
  !                   face_cv(periodic_list(n-n1+1),:,2) = face_cv(n,:,1)
  !
  !.. Formal Arguments ..
  integer,               intent(in) :: n1,n2
  integer,               intent(in) :: periodic_list(:)
  real(wp), allocatable, intent(in) :: rotation_matrix(:,:)
  real(wp),           intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (allocated(rotation_matrix)) then
    !
    ! Compute the rotated periodic boundary condition solution
    ! NOTE: A rotated periodic BC can only be 2D or 3D
    !
    if (nr == 3) then
      call get_bc_periodic_rotated_sol_3d(n1,n2,face_cv,periodic_list, &
                                          rotation_matrix)
    else if (nr == 2) then
      call get_bc_periodic_rotated_sol_2d(n1,n2,face_cv,periodic_list, &
                                          rotation_matrix)
    end if
    !
  else
    !
    ! Copy over the periodic boundary condition solution
    !
    if (nr == 3) then
      call get_bc_periodic_solution_3d(n1,n2,face_cv,periodic_list)
    else if (nr == 2) then
      call get_bc_periodic_solution_2d(n1,n2,face_cv,periodic_list)
    else
      call get_bc_periodic_solution_nd(n1,n2,face_cv,periodic_list)
    end if
    !
  end if
  !
end subroutine get_bc_periodic_solution
!
!###############################################################################
!
!pure &
subroutine get_bc_periodic_rotated_sol_3d(n1,n2,face_cv,periodic_list, &
                                          rotation_matrix)
  !
  ! #############################################################
  ! #############################################################
  ! ##  Rotated Periodic Boundary Condition
  ! #############################################################
  ! #############################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved variables at the flux points
  ! PERIODIC_LIST : For each face point n from n1 to n2, periodic_list(n) gives
  !                 the face point p that matches/corresponds to n based on
  !                 the periodic connectivity. Once this subroutine is finished,
  !                 the following should be true:
  !                   face_cv(n,:,2) = face_cv(periodic_list(n-n1+1),:,1)
  !                   face_cv(periodic_list(n-n1+1),:,2) = face_cv(n,:,1)
  ! ROTATION_MATRIX : coordinate transformation matrix used to rotate the
  !                   velocity vector for a rotated periodic BC
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp),    intent(in) :: rotation_matrix(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: rho,rhou,rhov,rhow,rhoe
  real(wp) :: rxx,rxy,rxz,ryx,ryy,ryz,rzx,rzy,rzz
  !
continue
  !
  ! Extract the components of the rotation matrix
  !
  rxx = rotation_matrix(1,1)
  rxy = rotation_matrix(1,2)
  rxz = rotation_matrix(1,3)
  !
  ryx = rotation_matrix(2,1)
  ryy = rotation_matrix(2,2)
  ryz = rotation_matrix(2,3)
  !
  rzx = rotation_matrix(3,1)
  rzy = rotation_matrix(3,2)
  rzz = rotation_matrix(3,3)
  !
#ifdef DEBUG_ON
!#define DEBUG_PERIODIC_ROTATIONS
#endif
#ifdef DEBUG_PERIODIC_ROTATIONS
  write (mypnum+7000,'(80("="))')
#endif
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,periodic_list)
  do n = n1,n2
    !
    rho  = face_cv(periodic_list(n-n1+1),nec,1)
    rhou = face_cv(periodic_list(n-n1+1),nmx,1)
    rhov = face_cv(periodic_list(n-n1+1),nmy,1)
    rhow = face_cv(periodic_list(n-n1+1),nmz,1)
    rhoe = face_cv(periodic_list(n-n1+1),nee,1)
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rhou*rxx + rhov*rxy + rhow*rxz
    face_cv(n,nmy,2) = rhou*ryx + rhov*ryy + rhow*ryz
    face_cv(n,nmz,2) = rhou*rzx + rhov*rzy + rhow*rzz
    face_cv(n,nee,2) = rhoe
#ifdef DEBUG_PERIODIC_ROTATIONS
    write (mypnum+7000,1) &
                          n, (face_cv(n,m,1),m=nmx,nmz), &
      periodic_list(n-n1+1), rhou,rhov,rhow, &
                             (face_cv(n,m,2),m=nmx,nmz)
#endif
    !
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv,periodic_list)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = face_cv(periodic_list(n-n1+1),m,1)
      end do
    end do
  end if
  !
  1 format (1x,"n=",i5," : ",3es12.4,"   ;", &
            1x,"k=",i5," : ",3es12.4," -> ",3es12.4)
  !
end subroutine get_bc_periodic_rotated_sol_3d
!
!###############################################################################
!
pure subroutine get_bc_periodic_rotated_sol_2d(n1,n2,face_cv,periodic_list, &
                                               rotation_matrix)
  !
  ! #############################################################
  ! #############################################################
  ! ##  Rotated Periodic Boundary Condition
  ! #############################################################
  ! #############################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved variables at the flux points
  ! PERIODIC_LIST : For each face point n from n1 to n2, periodic_list(n) gives
  !                 the face point p that matches/corresponds to n based on
  !                 the periodic connectivity. Once this subroutine is finished,
  !                 the following should be true:
  !                   face_cv(n,:,2) = face_cv(periodic_list(n-n1+1),:,1)
  !                   face_cv(periodic_list(n-n1+1),:,2) = face_cv(n,:,1)
  ! ROTATION_MATRIX : coordinate transformation matrix used to rotate the
  !                   velocity vector for a rotated periodic BC
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp),    intent(in) :: rotation_matrix(:,:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n
  real(wp) :: rho,rhou,rhov,rhoe
  real(wp) :: rxx,rxy,ryx,ryy
  !
continue
  !
  ! Extract the components of the rotation matrix
  !
  rxx = rotation_matrix(1,1)
  rxy = rotation_matrix(1,2)
  !
  ryx = rotation_matrix(2,1)
  ryy = rotation_matrix(2,2)
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,periodic_list)
  do n = n1,n2
    !
    rho  = face_cv(periodic_list(n-n1+1),nec,1)
    rhou = face_cv(periodic_list(n-n1+1),nmx,1)
    rhov = face_cv(periodic_list(n-n1+1),nmy,1)
    rhoe = face_cv(periodic_list(n-n1+1),nee,1)
    !
    face_cv(n,nec,2) = rho
    face_cv(n,nmx,2) = rhou*rxx + rhov*rxy
    face_cv(n,nmy,2) = rhou*ryx + rhov*ryy
    face_cv(n,nee,2) = rhoe
    !
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv,periodic_list)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = face_cv(periodic_list(n-n1+1),m,1)
      end do
    end do
  end if
  !
end subroutine get_bc_periodic_rotated_sol_2d
!
!###############################################################################
!
pure subroutine get_bc_periodic_solution_3d(n1,n2,face_cv,periodic_list)
  !
  ! #############################################################
  ! #############################################################
  ! ##  Periodic Boundary Condition
  ! #############################################################
  ! #############################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved variables at the flux points
  ! PERIODIC_LIST : For each face point n from n1 to n2, periodic_list(n) gives
  !                 the face point p that matches/corresponds to n based on
  !                 the periodic connectivity. Once this subroutine is finished,
  !                 the following should be true:
  !                   face_cv(n,:,2) = face_cv(periodic_list(n-n1+1),:,1)
  !                   face_cv(periodic_list(n-n1+1),:,2) = face_cv(n,:,1)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: m,n
  !
continue
  !
  ! NOTE: We may need to reverse the copy order, i.e.,
  !          gufp(:,k) = pufp(:,nfp+1-k)
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,periodic_list)
  do n = n1,n2
    face_cv(n,nec,2) =  face_cv(periodic_list(n-n1+1),nec,1)
    face_cv(n,nmx,2) =  face_cv(periodic_list(n-n1+1),nmx,1)
    face_cv(n,nmy,2) =  face_cv(periodic_list(n-n1+1),nmy,1)
    face_cv(n,nmz,2) =  face_cv(periodic_list(n-n1+1),nmz,1)
    face_cv(n,nee,2) =  face_cv(periodic_list(n-n1+1),nee,1)
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv,periodic_list)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = face_cv(periodic_list(n-n1+1),m,1)
      end do
    end do
  end if
  !
end subroutine get_bc_periodic_solution_3d
!
!###############################################################################
!
pure subroutine get_bc_periodic_solution_2d(n1,n2,face_cv,periodic_list)
  !
  ! #############################################################
  ! #############################################################
  ! ##  Periodic Boundary Condition
  ! #############################################################
  ! #############################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved variables at the flux points
  ! PERIODIC_LIST : For each face point n from n1 to n2, periodic_list(n) gives
  !                 the face point p that matches/corresponds to n based on
  !                 the periodic connectivity. Once this subroutine is finished,
  !                 the following should be true:
  !                   face_cv(n,:,2) = face_cv(periodic_list(n-n1+1),:,1)
  !                   face_cv(periodic_list(n-n1+1),:,2) = face_cv(n,:,1)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: m,n
  !
continue
  !
  ! NOTE: We may need to reverse the copy order, i.e.,
  !          gufp(:,k) = pufp(:,nfp+1-k)
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,periodic_list)
  do n = n1,n2
    face_cv(n,nec,2) =  face_cv(periodic_list(n-n1+1),nec,1)
    face_cv(n,nmx,2) =  face_cv(periodic_list(n-n1+1),nmx,1)
    face_cv(n,nmy,2) =  face_cv(periodic_list(n-n1+1),nmy,1)
    face_cv(n,nee,2) =  face_cv(periodic_list(n-n1+1),nee,1)
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv,periodic_list)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = face_cv(periodic_list(n-n1+1),m,1)
      end do
    end do
  end if
  !
end subroutine get_bc_periodic_solution_2d
!
!###############################################################################
!
pure subroutine get_bc_periodic_solution_nd(n1,n2,face_cv,periodic_list)
  !
  ! #############################################################
  ! #############################################################
  ! ##  Periodic Boundary Condition
  ! #############################################################
  ! #############################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved variables at the flux points
  ! PERIODIC_LIST : For each face point n from n1 to n2, periodic_list(n) gives
  !                 the face point p that matches/corresponds to n based on
  !                 the periodic connectivity. Once this subroutine is finished,
  !                 the following should be true:
  !                   face_cv(n,:,2) = face_cv(periodic_list(n-n1+1),:,1)
  !                   face_cv(periodic_list(n-n1+1),:,2) = face_cv(n,:,1)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: m,n
  !
continue
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(face_cv,periodic_list)
  do m = 1,size(face_cv,dim=2)
    do n = n1,n2
      face_cv(n,m,2) = face_cv(periodic_list(n-n1+1),m,1)
    end do
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv,periodic_list)
    do m = 1,ncs
      do n = n1,n2
        face_cv(n,m,2) = face_cv(periodic_list(n-n1+1),m,1)
      end do
    end do
  end if
  !
end subroutine get_bc_periodic_solution_nd
!
!###############################################################################
!
!pure &
subroutine get_bc_periodic_gradient(n1,n2,face_dcvdx,periodic_list, &
                                    rotation_matrix)
  !
  ! #############################################################
  ! #############################################################
  ! ##  Periodic Boundary Condition for gradients
  ! #############################################################
  ! #############################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
  ! PERIODIC_LIST : For each face point n from n1 to n2, periodic_list(n) gives
  !                 the face point p that matches/corresponds to n based on
  !                 the periodic connectivity. Once this subroutine is finished,
  !                 the following should be true:
  !           face_dcvdx(n,:,:,2) = face_dcvdx(periodic_list(n-n1+1),:,:,1)
  !           face_dcvdx(periodic_list(n-n1+1),:,:,2) = face_dcvdx(n,:,:,1)
  !
  !.. Formal Arguments ..
  integer,               intent(in) :: n1,n2
  integer,               intent(in) :: periodic_list(:)
  real(wp), allocatable, intent(in) :: rotation_matrix(:,:)
  real(wp),           intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = size(face_dcvdx,dim=2)
  !
  if (allocated(rotation_matrix)) then
    !
    ! Compute the rotated periodic boundary condition gradients
    ! NOTE: A rotated periodic BC can only be 2D or 3D
    !
    if (nr == 3) then
      call get_bc_periodic_rotated_grad_3d(n1,n2,face_dcvdx,periodic_list, &
                                           rotation_matrix)
    else if (nr == 2) then
      call get_bc_periodic_rotated_grad_2d(n1,n2,face_dcvdx,periodic_list, &
                                           rotation_matrix)
    end if
    !
  else
    !
    ! Copy over the periodic boundary condition gradients
    !
    if (nr == 3) then
      call get_bc_periodic_gradient_3d(n1,n2,face_dcvdx,periodic_list)
    else if (nr == 2) then
      call get_bc_periodic_gradient_2d(n1,n2,face_dcvdx,periodic_list)
    else
      call get_bc_periodic_gradient_nd(n1,n2,face_dcvdx,periodic_list)
    end if
    !
  end if
  !
end subroutine get_bc_periodic_gradient
!
!###############################################################################
!
!pure &
subroutine get_bc_periodic_rotated_grad_3d(n1,n2,face_dcvdx, &
                                           periodic_list,rotation_matrix)
  !
  ! #############################################################
  ! #############################################################
  ! ##  Periodic Boundary Condition for gradients
  ! #############################################################
  ! #############################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
  ! PERIODIC_LIST : For each face point n from n1 to n2, periodic_list(n) gives
  !                 the face point p that matches/corresponds to n based on
  !                 the periodic connectivity. Once this subroutine is finished,
  !                 the following should be true:
  !           face_dcvdx(n,:,:,2) = face_dcvdx(periodic_list(n-n1+1),:,:,1)
  !           face_dcvdx(periodic_list(n-n1+1),:,:,2) = face_dcvdx(n,:,:,1)
  ! ROTATION_MATRIX : coordinate transformation matrix used to rotate the
  !                   velocity vector for a rotated periodic BC
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp),    intent(in) :: rotation_matrix(:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: m,n
  real(wp) :: decdx,dmxdx,dmydx,dmzdx,deedx,rdmxdx,rdmydx,rdmzdx
  real(wp) :: decdy,dmxdy,dmydy,dmzdy,deedy,rdmxdy,rdmydy,rdmzdy
  real(wp) :: decdz,dmxdz,dmydz,dmzdz,deedz,rdmxdz,rdmydz,rdmzdz
  real(wp) :: rxx,rxy,rxz,ryx,ryy,ryz,rzx,rzy,rzz
  real(wp) :: axx,axy,axz,ayx,ayy,ayz,azx,azy,azz
  !
  !.. Local Arrays ..
  real(wp) :: rmat(1:3,1:3)
  real(wp) :: tnsr(1:3,1:3)
  real(wp) :: rslt(1:3,1:3)
  !
  !.. Local Parameters ..
  logical(lk), parameter :: use_matmul = true
 !logical(lk), parameter :: use_matmul = fals
  !
continue
  !
  ! Extract the components of the rotation matrix
  !
  rxx = rotation_matrix(1,1)
  rxy = rotation_matrix(1,2)
  rxz = rotation_matrix(1,3)
  !
  ryx = rotation_matrix(2,1)
  ryy = rotation_matrix(2,2)
  ryz = rotation_matrix(2,3)
  !
  rzx = rotation_matrix(3,1)
  rzy = rotation_matrix(3,2)
  rzz = rotation_matrix(3,3)
  !
 !rmat(1:3,1:3) = rotation_matrix(1:3,1:3)
  rmat(1:3,1:3) = transpose(rotation_matrix(1:3,1:3))
#ifdef DEBUG_PERIODIC_ROTATIONS
  write (mypnum+8000,'(80("#"))')
#endif
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx,periodic_list,rmat,rtrn,tnsr,rslt)
  do n = n1,n2
    !
    ! Extract the gradients from the matching periodic face point
    !
    decdx = face_dcvdx(periodic_list(n-n1+1),1,nec,1)
    decdy = face_dcvdx(periodic_list(n-n1+1),2,nec,1)
    decdz = face_dcvdx(periodic_list(n-n1+1),3,nec,1)
    !
    dmxdx = face_dcvdx(periodic_list(n-n1+1),1,nmx,1)
    dmxdy = face_dcvdx(periodic_list(n-n1+1),2,nmx,1)
    dmxdz = face_dcvdx(periodic_list(n-n1+1),3,nmx,1)
    !
    dmydx = face_dcvdx(periodic_list(n-n1+1),1,nmy,1)
    dmydy = face_dcvdx(periodic_list(n-n1+1),2,nmy,1)
    dmydz = face_dcvdx(periodic_list(n-n1+1),3,nmy,1)
    !
    dmzdx = face_dcvdx(periodic_list(n-n1+1),1,nmz,1)
    dmzdy = face_dcvdx(periodic_list(n-n1+1),2,nmz,1)
    dmzdz = face_dcvdx(periodic_list(n-n1+1),3,nmz,1)
    !
    deedx = face_dcvdx(periodic_list(n-n1+1),1,nee,1)
    deedy = face_dcvdx(periodic_list(n-n1+1),2,nee,1)
    deedz = face_dcvdx(periodic_list(n-n1+1),3,nee,1)
    !
    ! tnsr(1:3,1) = [dmxdx,dmxdy,dmxdz]
    ! tnsr(1:3,2) = [dmydx,dmydy,dmydz]
    ! tnsr(1:3,3) = [dmzdx,dmzdy,dmzdz]
    tnsr(1:3,1:3) = face_dcvdx(periodic_list(n-n1+1),1:3,nmx:nmz,1)
   !!
   !! tnsr(1:3,1) = [dmxdx,dmydx,dmzdx]
   !! tnsr(1:3,2) = [dmxdy,dmydy,dmzdy]
   !! tnsr(1:3,3) = [dmxdz,dmydz,dmzdz]
   !tnsr(:,:) = transpose(tnsr(:,:))
    !
    rslt(1:3,1:3) = matmul( rmat , matmul( tnsr , transpose(rmat) ) )
    !
    ! Compute matrix components resulting from matrix-matrix product of
    ! [rotation-matrix] x [velocity-gradient-tensor]
    !
    axx = dmxdx*rxx + dmxdy*rxy + dmxdz*rxz
    axy = dmydx*rxx + dmydy*rxy + dmydz*rxz
    axz = dmzdx*rxx + dmzdy*rxy + dmzdz*rxz
    !
    ayx = dmxdx*ryx + dmxdy*ryy + dmxdz*ryz
    ayy = dmydx*ryx + dmydy*ryy + dmydz*ryz
    ayz = dmzdx*ryx + dmzdy*ryy + dmzdz*ryz
    !
    azx = dmxdx*rzx + dmxdy*rzy + dmxdz*rzz
    azy = dmydx*rzx + dmydy*rzy + dmydz*rzz
    azz = dmzdx*rzx + dmzdy*rzy + dmzdz*rzz
    !
    rdmxdx = axx*rxx + axy*rxy + axz*rxz
    rdmxdy = ayx*rxx + ayy*rxy + ayz*rxz
    rdmxdz = azx*rxx + azy*rxy + azz*rxz
    !
    rdmydx = axx*ryx + axy*ryy + axz*ryz
    rdmydy = ayx*ryx + ayy*ryy + ayz*ryz
    rdmydz = azx*ryx + azy*ryy + azz*ryz
    !
    rdmzdx = axx*rzx + axy*rzy + axz*rzz
    rdmzdy = ayx*rzx + ayy*rzy + ayz*rzz
    rdmzdz = azx*rzx + azy*rzy + azz*rzz
    !
    ! Rotate the gradients of the matching periodic face point and save
    ! to face_dcvdx
    ! NOTE: The density and energy gradients are vector rotations which
    !       is just a matrix-vector product of
    !       [rotation-matrix] x [gradient-vector].
    !       The velocity gradient is a tensor rotation which is a
    !       matrix-matrix-matrix product of
    !       [rotation-matrix] x [grad-vel-tensor] x [Transpose(rotation-matrix)]
    !
    face_dcvdx(n,1,nec,2) = decdx*rxx + decdy*rxy + decdz*rxz
    face_dcvdx(n,2,nec,2) = decdx*ryx + decdy*ryy + decdz*ryz
    face_dcvdx(n,3,nec,2) = decdx*rzx + decdy*rzy + decdz*rzz
    !
    face_dcvdx(n,1,nmx,2) = merge( rslt(1,1) , rdmxdx , use_matmul )
    face_dcvdx(n,2,nmx,2) = merge( rslt(2,1) , rdmxdy , use_matmul )
    face_dcvdx(n,3,nmx,2) = merge( rslt(3,1) , rdmxdz , use_matmul )
    !
    face_dcvdx(n,1,nmy,2) = merge( rslt(1,2) , rdmydx , use_matmul )
    face_dcvdx(n,2,nmy,2) = merge( rslt(2,2) , rdmydy , use_matmul )
    face_dcvdx(n,3,nmy,2) = merge( rslt(3,2) , rdmydz , use_matmul )
    !
    face_dcvdx(n,1,nmz,2) = merge( rslt(1,3) , rdmzdx , use_matmul )
    face_dcvdx(n,2,nmz,2) = merge( rslt(2,3) , rdmzdy , use_matmul )
    face_dcvdx(n,3,nmz,2) = merge( rslt(3,3) , rdmzdz , use_matmul )
    !
    face_dcvdx(n,1,nee,2) = deedx*rxx + deedy*rxy + deedz*rxz
    face_dcvdx(n,2,nee,2) = deedx*ryx + deedy*ryy + deedz*ryz
    face_dcvdx(n,3,nee,2) = deedx*rzx + deedy*rzy + deedz*rzz
#ifdef DEBUG_PERIODIC_ROTATIONS
block
  character(len=*), parameter :: ceq(1:5) = ["ec","mx","my","mz","ee"]
    write (mypnum+8000,1) n,periodic_list(n-n1+1), &
         ceq(nec),(face_dcvdx(n,m,nec,1),m=1,3), &
                  (face_dcvdx(n,m,nec,2),m=1,3), &
                  decdx,decdy,decdz, &
         ceq(nmx),(face_dcvdx(n,m,nmx,1),m=1,3), &
                  (face_dcvdx(n,m,nmx,2),m=1,3), &
                  dmxdx,dmxdy,dmxdz, &
         ceq(nmy),(face_dcvdx(n,m,nmy,1),m=1,3), &
                  (face_dcvdx(n,m,nmy,2),m=1,3), &
                  dmydx,dmydy,dmydz, &
         ceq(nmz),(face_dcvdx(n,m,nmz,1),m=1,3), &
                  (face_dcvdx(n,m,nmz,2),m=1,3), &
                  dmzdx,dmzdy,dmzdz, &
         ceq(nee),(face_dcvdx(n,m,nee,1),m=1,3), &
                  (face_dcvdx(n,m,nee,2),m=1,3), &
                  deedx,deedy,deedz
end block
  1 format (80("="),/, &
            21x,"n = ",i6,14x,";",16x,"k = ",i6,/, &
            80("-"),/, &
            5(" d",a2,3x,3es12.4,"  ;  ",3es12.4," <- ",3es12.4,:,/))
#endif
    !
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_dcvdx,periodic_list)
    do m = 1,ncs
      do n = n1,n2
        !
        ! Extract the gradients from the matching periodic face point
        !
        decdx = face_dcvdx(periodic_list(n-n1+1),1,m,1)
        decdy = face_dcvdx(periodic_list(n-n1+1),2,m,1)
        decdz = face_dcvdx(periodic_list(n-n1+1),3,m,1)
        !
        face_dcvdx(n,1,m,2) = decdx*rxx + decdy*rxy + decdz*rxz
        face_dcvdx(n,2,m,2) = decdx*ryx + decdy*ryy + decdz*ryz
        face_dcvdx(n,3,m,2) = decdx*rzx + decdy*rzy + decdz*rzz
        !
      end do
    end do
  end if
  !
end subroutine get_bc_periodic_rotated_grad_3d
!
!###############################################################################
!
pure subroutine get_bc_periodic_rotated_grad_2d(n1,n2,face_dcvdx, &
                                                periodic_list,rotation_matrix)
  !
  ! #############################################################
  ! #############################################################
  ! ##  Periodic Boundary Condition for gradients
  ! #############################################################
  ! #############################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
  ! PERIODIC_LIST : For each face point n from n1 to n2, periodic_list(n) gives
  !                 the face point p that matches/corresponds to n based on
  !                 the periodic connectivity. Once this subroutine is finished,
  !                 the following should be true:
  !           face_dcvdx(n,:,:,2) = face_dcvdx(periodic_list(n-n1+1),:,:,1)
  !           face_dcvdx(periodic_list(n-n1+1),:,:,2) = face_dcvdx(n,:,:,1)
  ! ROTATION_MATRIX : coordinate transformation matrix used to rotate the
  !                   velocity vector for a rotated periodic BC
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp),    intent(in) :: rotation_matrix(:,:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: m,n
  real(wp) :: decdx,dmxdx,dmydx,deedx,rdmxdx,rdmydx
  real(wp) :: decdy,dmxdy,dmydy,deedy,rdmxdy,rdmydy
  real(wp) :: rxx,rxy,ryx,ryy
  real(wp) :: axx,axy,ayx,ayy
  !
  !.. Local Arrays ..
  real(wp) :: rmat(1:2,1:2)
  real(wp) :: tnsr(1:2,1:2)
  real(wp) :: rslt(1:2,1:2)
  !
  !.. Local Parameters ..
  logical(lk), parameter :: use_matmul = true
 !logical(lk), parameter :: use_matmul = fals
  !
continue
  !
  ! Extract the components of the rotation matrix
  !
  rxx = rotation_matrix(1,1)
  rxy = rotation_matrix(1,2)
  !
  ryx = rotation_matrix(2,1)
  ryy = rotation_matrix(2,2)
  !
 !rmat(1:2,1:2) = rotation_matrix(1:2,1:2)
  rmat(1:2,1:2) = transpose(rotation_matrix(1:2,1:2))
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx,periodic_list,rmat,rtrn,tnsr,rslt)
  do n = n1,n2
    !
    ! Extract the gradients from the matching periodic face point
    !
    decdx = face_dcvdx(periodic_list(n-n1+1),1,nec,1)
    decdy = face_dcvdx(periodic_list(n-n1+1),2,nec,1)
    !
    dmxdx = face_dcvdx(periodic_list(n-n1+1),1,nmx,1)
    dmxdy = face_dcvdx(periodic_list(n-n1+1),2,nmx,1)
    !
    dmydx = face_dcvdx(periodic_list(n-n1+1),1,nmy,1)
    dmydy = face_dcvdx(periodic_list(n-n1+1),2,nmy,1)
    !
    deedx = face_dcvdx(periodic_list(n-n1+1),1,nee,1)
    deedy = face_dcvdx(periodic_list(n-n1+1),2,nee,1)
    !
    ! tnsr(1:2,1) = [dmxdx,dmxdy]
    ! tnsr(1:2,2) = [dmydx,dmydy]
    tnsr(1:2,1:2) = face_dcvdx(periodic_list(n-n1+1),1:2,nmx:nmy,1)
   !!
   !! tnsr(1:2,1) = [dmxdx,dmydx]
   !! tnsr(1:2,2) = [dmxdy,dmydy]
   !tnsr(:,:) = transpose(tnsr(:,:))
    !
    rslt(1:2,1:2) = matmul( rmat , matmul( tnsr , transpose(rmat) ) )
    !
    ! Compute matrix components resulting from matrix-matrix product of
    ! [rotation-matrix] x [velocity-gradient-tensor]
    !
    axx = dmxdx*rxx + dmxdy*rxy
    axy = dmydx*rxx + dmydy*rxy
    !
    ayx = dmxdx*ryx + dmxdy*ryy
    ayy = dmydx*ryx + dmydy*ryy
    !
    rdmxdx = axx*rxx + axy*rxy
    rdmxdy = ayx*rxx + ayy*rxy
    !
    rdmydx = axx*ryx + axy*ryy
    rdmydy = ayx*ryx + ayy*ryy
    !
    ! Rotate the gradients of the matching periodic face point and save
    ! to face_dcvdx
    ! NOTE: The density and energy gradients are vector rotations which
    !       is just a matrix-vector product of
    !       [rotation-matrix] x [gradient-vector].
    !       The velocity gradient is a tensor rotation which is a
    !       matrix-matrix-matrix product of
    !       [rotation-matrix] x [grad-vel-tensor] x [Transpose(rotation-matrix)]
    !
    face_dcvdx(n,1,nec,2) = decdx*rxx + decdy*rxy
    face_dcvdx(n,2,nec,2) = decdx*ryx + decdy*ryy
    !
    face_dcvdx(n,1,nmx,2) = merge( rslt(1,1) , rdmxdx , use_matmul )
    face_dcvdx(n,2,nmx,2) = merge( rslt(2,1) , rdmxdy , use_matmul )
    !
    face_dcvdx(n,1,nmy,2) = merge( rslt(1,2) , rdmydx , use_matmul )
    face_dcvdx(n,2,nmy,2) = merge( rslt(2,2) , rdmydy , use_matmul )
    !
    face_dcvdx(n,1,nee,2) = deedx*rxx + deedy*rxy
    face_dcvdx(n,2,nee,2) = deedx*ryx + deedy*ryy
    !
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_dcvdx,periodic_list)
    do m = 1,ncs
      do n = n1,n2
        !
        ! Extract the gradients from the matching periodic face point
        !
        decdx = face_dcvdx(periodic_list(n-n1+1),1,m,1)
        decdy = face_dcvdx(periodic_list(n-n1+1),2,m,1)
        !
        face_dcvdx(n,1,m,2) = decdx*rxx + decdy*rxy
        face_dcvdx(n,2,m,2) = decdx*ryx + decdy*ryy
        !
      end do
    end do
  end if
  !
end subroutine get_bc_periodic_rotated_grad_2d
!
!###############################################################################
!
pure subroutine get_bc_periodic_gradient_3d(n1,n2,face_dcvdx,periodic_list)
  !
  ! #############################################################
  ! #############################################################
  ! ##  Periodic Boundary Condition for gradients
  ! #############################################################
  ! #############################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
  ! PERIODIC_LIST : For each face point n from n1 to n2, periodic_list(n) gives
  !                 the face point p that matches/corresponds to n based on
  !                 the periodic connectivity. Once this subroutine is finished,
  !                 the following should be true:
  !           face_dcvdx(n,:,:,2) = face_dcvdx(periodic_list(n-n1+1),:,:,1)
  !           face_dcvdx(periodic_list(n-n1+1),:,:,2) = face_dcvdx(n,:,:,1)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: m,n
  !
continue
  !
  ! Simply copy the interior gradients to the boundary flux points
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx,periodic_list)
  do n = n1,n2
    !
    face_dcvdx(n,1,nec,2) = face_dcvdx(periodic_list(n-n1+1),1,nec,1)
    face_dcvdx(n,2,nec,2) = face_dcvdx(periodic_list(n-n1+1),2,nec,1)
    face_dcvdx(n,3,nec,2) = face_dcvdx(periodic_list(n-n1+1),3,nec,1)
    !
    face_dcvdx(n,1,nmx,2) = face_dcvdx(periodic_list(n-n1+1),1,nmx,1)
    face_dcvdx(n,2,nmx,2) = face_dcvdx(periodic_list(n-n1+1),2,nmx,1)
    face_dcvdx(n,3,nmx,2) = face_dcvdx(periodic_list(n-n1+1),3,nmx,1)
    !
    face_dcvdx(n,1,nmy,2) = face_dcvdx(periodic_list(n-n1+1),1,nmy,1)
    face_dcvdx(n,2,nmy,2) = face_dcvdx(periodic_list(n-n1+1),2,nmy,1)
    face_dcvdx(n,3,nmy,2) = face_dcvdx(periodic_list(n-n1+1),3,nmy,1)
    !
    face_dcvdx(n,1,nmz,2) = face_dcvdx(periodic_list(n-n1+1),1,nmz,1)
    face_dcvdx(n,2,nmz,2) = face_dcvdx(periodic_list(n-n1+1),2,nmz,1)
    face_dcvdx(n,3,nmz,2) = face_dcvdx(periodic_list(n-n1+1),3,nmz,1)
    !
    face_dcvdx(n,1,nee,2) = face_dcvdx(periodic_list(n-n1+1),1,nee,1)
    face_dcvdx(n,2,nee,2) = face_dcvdx(periodic_list(n-n1+1),2,nee,1)
    face_dcvdx(n,3,nee,2) = face_dcvdx(periodic_list(n-n1+1),3,nee,1)
    !
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_dcvdx,periodic_list)
    do m = 1,ncs
      do n = n1,n2
        face_dcvdx(n,1,m,2) = face_dcvdx(periodic_list(n-n1+1),1,m,1)
        face_dcvdx(n,2,m,2) = face_dcvdx(periodic_list(n-n1+1),2,m,1)
        face_dcvdx(n,3,m,2) = face_dcvdx(periodic_list(n-n1+1),3,m,1)
      end do
    end do
  end if
  !
end subroutine get_bc_periodic_gradient_3d
!
!###############################################################################
!
pure subroutine get_bc_periodic_gradient_2d(n1,n2,face_dcvdx,periodic_list)
  !
  ! #############################################################
  ! #############################################################
  ! ##  Periodic Boundary Condition for gradients
  ! #############################################################
  ! #############################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
  ! PERIODIC_LIST : For each face point n from n1 to n2, periodic_list(n) gives
  !                 the face point p that matches/corresponds to n based on
  !                 the periodic connectivity. Once this subroutine is finished,
  !                 the following should be true:
  !           face_dcvdx(n,:,:,2) = face_dcvdx(periodic_list(n-n1+1),:,:,1)
  !           face_dcvdx(periodic_list(n-n1+1),:,:,2) = face_dcvdx(n,:,:,1)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: m,n
  !
continue
  !
  ! Simply copy the interior gradients to the boundary flux points
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx,periodic_list)
  do n = n1,n2
    !
    face_dcvdx(n,1,nec,2) = face_dcvdx(periodic_list(n-n1+1),1,nec,1)
    face_dcvdx(n,2,nec,2) = face_dcvdx(periodic_list(n-n1+1),2,nec,1)
    !
    face_dcvdx(n,1,nmx,2) = face_dcvdx(periodic_list(n-n1+1),1,nmx,1)
    face_dcvdx(n,2,nmx,2) = face_dcvdx(periodic_list(n-n1+1),2,nmx,1)
    !
    face_dcvdx(n,1,nmy,2) = face_dcvdx(periodic_list(n-n1+1),1,nmy,1)
    face_dcvdx(n,2,nmy,2) = face_dcvdx(periodic_list(n-n1+1),2,nmy,1)
    !
    face_dcvdx(n,1,nee,2) = face_dcvdx(periodic_list(n-n1+1),1,nee,1)
    face_dcvdx(n,2,nee,2) = face_dcvdx(periodic_list(n-n1+1),2,nee,1)
    !
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_dcvdx,periodic_list)
    do m = 1,ncs
      do n = n1,n2
        face_dcvdx(n,1,m,2) = face_dcvdx(periodic_list(n-n1+1),1,m,1)
        face_dcvdx(n,2,m,2) = face_dcvdx(periodic_list(n-n1+1),2,m,1)
      end do
    end do
  end if
  !
end subroutine get_bc_periodic_gradient_2d
!
!###############################################################################
!
pure subroutine get_bc_periodic_gradient_nd(n1,n2,face_dcvdx,periodic_list)
  !
  ! #############################################################
  ! #############################################################
  ! ##  Periodic Boundary Condition for gradients
  ! #############################################################
  ! #############################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
  ! PERIODIC_LIST : For each face point n from n1 to n2, periodic_list(n) gives
  !                 the face point p that matches/corresponds to n based on
  !                 the periodic connectivity. Once this subroutine is finished,
  !                 the following should be true:
  !           face_dcvdx(n,:,:,2) = face_dcvdx(periodic_list(n-n1+1),:,:,1)
  !           face_dcvdx(periodic_list(n-n1+1),:,:,2) = face_dcvdx(n,:,:,1)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: l,m,n
  !
continue
  !
  ! Simply copy the interior gradients to the boundary flux points
  !
  !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(face_dcvdx,periodic_list)
  do m = 1,size(face_dcvdx,dim=3)
    do l = 1,size(face_dcvdx,dim=2)
      do n = n1,n2
        face_dcvdx(n,l,m,2) = face_dcvdx(periodic_list(n-n1+1),l,m,1)
      end do
    end do
  end do
  !
end subroutine get_bc_periodic_gradient_nd
!
!###############################################################################
!
end module bc_periodic_mod
