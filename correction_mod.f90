module correction_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use generic_types_mod, only : matrix
  !
  implicit none
  !
  private
  !
  ! correct : correction matrix for all combinations of cell geometry and 
  !           and cell order that are possible in the current simulation
  !
  !     For cell 'n' :
  !
  !         The correction matrix (giving the contribution to each solution
  !         point within cell 'n' from the interface flux point 'k' that is
  !         currently of interest) is given by the expression
  !
  !         correct(this_geom,this_order)%mat(:,k)
  !
  !     Where :
  !
  !        this_geom = cell_geom(n)
  !        this_order = cell_order(n)
  !        k = the flux point on a face of cell 'n' that is currently of interest
  !
  type(matrix), public, save, target, allocatable :: correct
  !
end module correction_mod
!
!###############################################################################
!
