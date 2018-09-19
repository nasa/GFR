module bc_fixed_mod
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
  public :: get_bc_fixed_solution
  !
contains
!
!###############################################################################
!
pure subroutine get_bc_fixed_solution(n1,n2,face_cv,cv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Fixed boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved variables at the flux points
  ! CV_IN   : Conserved variables specified for this BC
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0) then
    call get_bc_fixed_solution_3d(n1,n2,face_cv,cv_in)
  else if (nr == 2 .and. ncs == 0) then
    call get_bc_fixed_solution_2d(n1,n2,face_cv,cv_in)
  else
    call get_bc_fixed_solution_nd(n1,n2,face_cv,cv_in)
  end if
  !
end subroutine get_bc_fixed_solution
!
!###############################################################################
!
pure subroutine get_bc_fixed_solution_3d(n1,n2,face_cv,cv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Fixed boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved variables at the flux points
  ! CV_IN   : Conserved variables specified for this BC
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:)
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! Fix all variables to the specified flow conditions for this boundary
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,cv_in)
  do n = n1,n2
    face_cv(n,nec) = cv_in(nec)
    face_cv(n,nmx) = cv_in(nmx)
    face_cv(n,nmy) = cv_in(nmy)
    face_cv(n,nmz) = cv_in(nmz)
    face_cv(n,nee) = cv_in(nee)
  end do
  !
end subroutine get_bc_fixed_solution_3d
!
!###############################################################################
!
pure subroutine get_bc_fixed_solution_2d(n1,n2,face_cv,cv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Fixed boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved variables at the flux points
  ! CV_IN   : Conserved variables specified for this BC
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:)
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! Fix all variables to the specified flow conditions for this boundary
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,cv_in)
  do n = n1,n2
    face_cv(n,nec) = cv_in(nec)
    face_cv(n,nmx) = cv_in(nmx)
    face_cv(n,nmy) = cv_in(nmy)
    face_cv(n,nee) = cv_in(nee)
  end do
  !
end subroutine get_bc_fixed_solution_2d
!
!###############################################################################
!
pure subroutine get_bc_fixed_solution_nd(n1,n2,face_cv,cv_in)
  !
  ! #####################################################
  ! #####################################################
  ! ##  Fixed boundary condition
  ! #####################################################
  ! #####################################################
  !
  ! N1      : Beginning face point for this boundary condition
  ! N2      : Ending    face point for this boundary condition
  ! FACE_CV : Conserved variables at the flux points
  ! CV_IN   : Conserved variables specified for this BC
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: cv_in(:)
  real(wp), intent(inout) :: face_cv(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n
  !
continue
  !
  ! Fix all variables to the specified flow conditions for this boundary
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv,cv_in)
  do m = 1,size(cv_in)
    do n = n1,n2
      face_cv(n,m) = cv_in(m)
    end do
  end do
  !
end subroutine get_bc_fixed_solution_nd
!
!###############################################################################
!
end module bc_fixed_mod
