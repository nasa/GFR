module bc_periodic_mod
  !
  !.. Use Statements
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee
  use eqn_idx, only : ncs
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
pure subroutine get_bc_periodic_solution(n1,n2,face_cv,periodic_list)
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
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    call get_bc_periodic_solution_3d(n1,n2,face_cv,periodic_list)
  else if (nr == 2) then
    call get_bc_periodic_solution_2d(n1,n2,face_cv,periodic_list)
  else
    call get_bc_periodic_solution_nd(n1,n2,face_cv,periodic_list)
  end if
  !
end subroutine get_bc_periodic_solution
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
#ifdef REVERSE_PERIODIC
    face_cv(n,nmz,2) = -face_cv(periodic_list(n-n1+1),nmz,1)
#else
    face_cv(n,nmz,2) =  face_cv(periodic_list(n-n1+1),nmz,1)
#endif
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
#ifdef REVERSE_PERIODIC
    face_cv(n,nmy,2) = -face_cv(periodic_list(n-n1+1),nmy,1)
#else
    face_cv(n,nmy,2) =  face_cv(periodic_list(n-n1+1),nmy,1)
#endif
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
#ifdef REVERSE_PERIODIC
  !
  ! NOTE: We may need to reverse the copy order, i.e.,
  !          gufp(:,k) = pufp(:,nfp+1-k)
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,periodic_list)
  do n = n1,n2
    face_cv(n,nme,2) = -face_cv(periodic_list(n-n1+1),nme,1)
  end do
#endif
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
pure subroutine get_bc_periodic_gradient(n1,n2,face_dcvdx,periodic_list)
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
  integer :: nr
  !
continue
  !
  nr = size(face_dcvdx,dim=2)
  !
  if (nr == 3) then
    call get_bc_periodic_gradient_3d(n1,n2,face_dcvdx,periodic_list)
  else if (nr == 2) then
    call get_bc_periodic_gradient_2d(n1,n2,face_dcvdx,periodic_list)
  else
    call get_bc_periodic_gradient_nd(n1,n2,face_dcvdx,periodic_list)
  end if
  !
end subroutine get_bc_periodic_gradient
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
