module projection_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use generic_types_mod, only : matrix
  !
  implicit none
  !
  private
  !
  ! projection_matrix : derived type that stores the matrices for
  !                     projecting a polynomial quantity within a
  !                     cell to a face point or cell point for a given
  !                     combination of cell geometry and cell order
  !
  type, public :: projection_matrix
    type(matrix), allocatable :: toFace(:) ! project to a face point
    type(matrix), allocatable :: toCell(:) ! project to a cell point
  end type projection_matrix
  !
  ! projct : projection matrix for all combinations of cell geometry and 
  !          and cell order that are possible in the current simulation
  !
  !     For cell 'n' :
  !
  !         The projection matrix (giving the contribution of each solution
  !         point within cell 'n' to the solution point 'k' that is currently
  !         of interest) is given by the expression
  !
  !         projct(this_geom,this_order)%toFace(to_order)%mat(:,k)
  !         projct(this_geom,this_order)%toCell(to_order)%mat(:,k)
  !
  !     Where :
  !
  !        this_geom = cell_geom(n)
  !        this_order = cell_order(n)
  !        to_order = order of face or cell to which we are projecting
  !        k = if using toFace, the face solution point we are projecting to
  !            if using toCell, the cell solution point we are projecting to
  !
  type(projection_matrix), public, save, target, allocatable :: projct(:,:)
  !
end module projection_mod
!
!###############################################################################
!
