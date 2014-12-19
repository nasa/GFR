module derivatives_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use generic_types_mod, only : matrix
  !
  implicit none
  !
  private
  !
  ! derivative_matrices : derived type that stores the derivative matrices of
  !                       each local coordinate direction for a given
  !                       combination of cell geometry and cell order
  !
  type, public :: derivative_matrices
    type(matrix), allocatable :: d(:)
  end type derivative_matrices
  !
  ! dmat : derivative matrix for all combinations of cell geometry and 
  !        and cell order that are possible in the current simulation
  !
  !     For cell 'n' :
  !
  !         The derivative matrix (giving the contribution of each solution
  !         point within cell 'n' to the solution point 'k' that is currently
  !         of interest) is given by the expression
  !
  !         dmat(this_geom,this_order)%d(l)%mat(:,k)
  !
  !     Where :
  !
  !        this_geom = cell_geom(n)
  !        this_order = cell_order(n)
  !        l = the local coordinate direction within the reference cell,
  !            i.e., (1,2,3) for (xi,eta,zeta) or (r,s,t)
  !        k = the solution point within cell 'n' that is currently of interest
  !
  type(derivative_matrices), public, save, target, allocatable :: dmat(:,:)
  !
  ! dqmat : same as dmat but gives the derivative matrices for the quadrature
  !         points
  !
  type(derivative_matrices), public, save, target, allocatable :: dqmat(:,:)
  !
end module derivatives_mod
!
!###############################################################################
!
