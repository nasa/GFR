module bc_sup_out_mod
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
  public :: get_bc_sup_outflow_solution
  !
contains
!
!###############################################################################
!
pure subroutine get_bc_sup_outflow_solution(n1,n2,face_cv)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic outflow boundary condition
  ! ##    - Extrapolate all variables
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved variables at the flux points
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0) then
    call get_bc_sup_outflow_solution_3d(n1,n2,face_cv)
  else if (nr == 2 .and. ncs == 0) then
    call get_bc_sup_outflow_solution_2d(n1,n2,face_cv)
  else
    call get_bc_sup_outflow_solution_nd(n1,n2,face_cv)
  end if
  !
end subroutine get_bc_sup_outflow_solution
!
!###############################################################################
!
pure subroutine get_bc_sup_outflow_solution_3d(n1,n2,face_cv)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic outflow boundary condition
  ! ##    - Extrapolate all variables
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved variables at the flux points
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv)
  do n = n1,n2
    face_cv(n,nec,2) = face_cv(n,nec,1)
    face_cv(n,nmx,2) = face_cv(n,nmx,1)
    face_cv(n,nmy,2) = face_cv(n,nmy,1)
    face_cv(n,nmz,2) = face_cv(n,nmz,1)
    face_cv(n,nee,2) = face_cv(n,nee,1)
  end do
  !
end subroutine get_bc_sup_outflow_solution_3d
!
!###############################################################################
!
pure subroutine get_bc_sup_outflow_solution_2d(n1,n2,face_cv)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic outflow boundary condition
  ! ##    - Extrapolate all variables
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved variables at the flux points
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv)
  do n = n1,n2
    face_cv(n,nec,2) = face_cv(n,nec,1)
    face_cv(n,nmx,2) = face_cv(n,nmx,1)
    face_cv(n,nmy,2) = face_cv(n,nmy,1)
    face_cv(n,nee,2) = face_cv(n,nee,1)
  end do
  !
end subroutine get_bc_sup_outflow_solution_2d
!
!###############################################################################
!
pure subroutine get_bc_sup_outflow_solution_nd(n1,n2,face_cv)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Supersonic outflow boundary condition
  ! ##    - Extrapolate all variables
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved variables at the flux points
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  !.. Local Scalars ..
  integer :: m,n
  !
continue
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv)
  do m = 1,size(face_cv,dim=2)
    do n = n1,n2
      face_cv(n,m,2) = face_cv(n,m,1)
    end do
  end do
  !
end subroutine get_bc_sup_outflow_solution_nd
!
!###############################################################################
!
end module bc_sup_out_mod
