module interpolation_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use generic_types_mod, only : matrix
  !
  implicit none
  !
  private
  !
  ! interpolation_matrix : derived type that stores the matrices for
  !                        interpolating a polynomial quantity within a
  !                        cell to a face point or cell point for a given
  !                        combination of cell geometry and cell order
  !
  type, public :: interpolation_matrix
    type(matrix), allocatable :: toFace(:) ! interpolate to a face point
    type(matrix), allocatable :: toCell(:) ! interpolate to a cell point
  end type interpolation_matrix
  !
  ! interp : interpolation matrix for all combinations of cell geometry
  !          and cell order that are possible in the current simulation
  !
  !     For cell 'n' :
  !
  !         The interpolation matrix (giving the contribution of each solution
  !         point within cell 'n' to the solution point 'k' that is currently
  !         of interest) is given by the expression
  !
  !         interp(this_geom,this_order)%toFace(to_order)%mat(:,k)
  !         interp(this_geom,this_order)%toCell(to_order)%mat(:,k)
  !
  !     Where :
  !
  !        this_geom = cell_geom(n)
  !        this_order = cell_order(n)
  !        to_order = order of face or cell to which we are interpolating
  !        k = if using toFace, the face solution point we are interpolating to
  !            if using toCell, the cell solution point we are interpolating to
  !
  type(interpolation_matrix), public, save, target, allocatable :: interp(:,:)
  !
  ! outerp : "outerpolation" matrix for all combinations of cell geometry
  !          and cell order that are possible in the current simulation.
  !          These matrices interpolate the solution from the cell solution
  !          points to a set of nodal points only meant for output and to
  !          visualize the solution. Get it?! Outerpolation! Interpolate to
  !          the output points! Outerpolation! I should probably go home and
  !          get some sleep...
  !
  type(interpolation_matrix), public, save, target, allocatable :: outerp(:,:)
  !
end module interpolation_mod
!
!###############################################################################
!
