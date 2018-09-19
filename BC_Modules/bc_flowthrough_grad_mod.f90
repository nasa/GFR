module bc_flowthrough_grad_mod
  !
  !.. Use Statements
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nee
  use eqn_idx, only : ncs
  !
  implicit none
  !
  private
  !
  public :: get_bc_flowthrough_gradient
  !
contains
!
!###############################################################################
!
pure subroutine get_bc_flowthrough_gradient(n1,n2,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Flow-through boundary condition for gradients
  ! #####################################################
  ! #####################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
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
    call get_bc_flowthrough_gradient_3d(n1,n2,face_dcvdx)
  else if (nr == 2) then
    call get_bc_flowthrough_gradient_2d(n1,n2,face_dcvdx)
  else
    call get_bc_flowthrough_gradient_nd(n1,n2,face_dcvdx)
  end if
  !
end subroutine get_bc_flowthrough_gradient
!
!###############################################################################
!
pure subroutine get_bc_flowthrough_gradient_3d(n1,n2,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Flow-through boundary condition for gradients
  ! #####################################################
  ! #####################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: m,n
  !
continue
  !
  ! Simply copy the interior gradients to the boundary flux points
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
  do n = n1,n2
    !
    face_dcvdx(n,1,nec,2) = face_dcvdx(n,1,nec,1)
    face_dcvdx(n,2,nec,2) = face_dcvdx(n,2,nec,1)
    face_dcvdx(n,3,nec,2) = face_dcvdx(n,3,nec,1)
    !
    face_dcvdx(n,1,nmx,2) = face_dcvdx(n,1,nmx,1)
    face_dcvdx(n,2,nmx,2) = face_dcvdx(n,2,nmx,1)
    face_dcvdx(n,3,nmx,2) = face_dcvdx(n,3,nmx,1)
    !
    face_dcvdx(n,1,nmy,2) = face_dcvdx(n,1,nmy,1)
    face_dcvdx(n,2,nmy,2) = face_dcvdx(n,2,nmy,1)
    face_dcvdx(n,3,nmy,2) = face_dcvdx(n,3,nmy,1)
    !
    face_dcvdx(n,1,nmz,2) = face_dcvdx(n,1,nmz,1)
    face_dcvdx(n,2,nmz,2) = face_dcvdx(n,2,nmz,1)
    face_dcvdx(n,3,nmz,2) = face_dcvdx(n,3,nmz,1)
    !
    face_dcvdx(n,1,nee,2) = face_dcvdx(n,1,nee,1)
    face_dcvdx(n,2,nee,2) = face_dcvdx(n,2,nee,1)
    face_dcvdx(n,3,nee,2) = face_dcvdx(n,3,nee,1)
    !
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_dcvdx)
    do m = 1,ncs
      do n = n1,n2
        face_dcvdx(n,1,m,2) = face_dcvdx(n,1,m,1)
        face_dcvdx(n,2,m,2) = face_dcvdx(n,2,m,1)
        face_dcvdx(n,3,m,2) = face_dcvdx(n,3,m,1)
      end do
    end do
  end if
  !
end subroutine get_bc_flowthrough_gradient_3d
!
!###############################################################################
!
pure subroutine get_bc_flowthrough_gradient_2d(n1,n2,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Flow-through boundary condition for gradients
  ! #####################################################
  ! #####################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: m,n
  !
continue
  !
  ! Simply copy the interior gradients to the boundary flux points
  !
  !$ACC PARALLEL LOOP PRESENT(face_dcvdx)
  do n = n1,n2
    !
    face_dcvdx(n,1,nec,2) = face_dcvdx(n,1,nec,1)
    face_dcvdx(n,2,nec,2) = face_dcvdx(n,2,nec,1)
    !
    face_dcvdx(n,1,nmx,2) = face_dcvdx(n,1,nmx,1)
    face_dcvdx(n,2,nmx,2) = face_dcvdx(n,2,nmx,1)
    !
    face_dcvdx(n,1,nmy,2) = face_dcvdx(n,1,nmy,1)
    face_dcvdx(n,2,nmy,2) = face_dcvdx(n,2,nmy,1)
    !
    face_dcvdx(n,1,nee,2) = face_dcvdx(n,1,nee,1)
    face_dcvdx(n,2,nee,2) = face_dcvdx(n,2,nee,1)
    !
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_dcvdx)
    do m = 1,ncs
      do n = n1,n2
        face_dcvdx(n,1,m,2) = face_dcvdx(n,1,m,1)
        face_dcvdx(n,2,m,2) = face_dcvdx(n,2,m,1)
      end do
    end do
  end if
  !
end subroutine get_bc_flowthrough_gradient_2d
!
!###############################################################################
!
pure subroutine get_bc_flowthrough_gradient_nd(n1,n2,face_dcvdx)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Flow-through boundary condition for gradients
  ! #####################################################
  ! #####################################################
  !
  ! N1         : Beginning face point for this boundary condition
  ! N2         : Ending    face point for this boundary condition
  ! FACE_DCVDX : Gradients of the conserved variables at the flux points
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_dcvdx(:,:,:,:)
  !
  !.. Local Scalars ..
  integer :: n,m,l
  !
continue
  !
  ! Simply copy the interior gradients to the boundary flux points
  !
  !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(face_dcvdx)
  do m = 1,size(face_dcvdx,dim=3)
    do l = 1,size(face_dcvdx,dim=2)
      do n = n1,n2
        face_dcvdx(n,l,m,2) = face_dcvdx(n,l,m,1)
      end do
    end do
  end do
  !
end subroutine get_bc_flowthrough_gradient_nd
!
!###############################################################################
!
end module bc_flowthrough_grad_mod
