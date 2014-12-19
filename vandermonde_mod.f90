module vandermonde_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use generic_types_mod, only : matrix
  !
  implicit none
  !
  private
  !
  ! vandermonde_matrix : derived type that stores the matrices for converting
  !                      the polynomial solution within a cell between the nodal
  !                      and modal forms (two matrices, one for converting in
  !                      each of the two directions) for a given combination of
  !                      cell geometry and cell order.
  !
  ! vandermonde_matrix%nodal2modal : converts a nodal solution to modal form
  ! vandermonde_matrix%modal2nodal : converts a modal solution to nodal form
  !
  ! NOTE: the modal2nodal matrix is the 'actual' generalized Vandermonde matrix
  !       whereas the nodal2modal matrix is the inverse of the generalized
  !       Vandermonde matrix
  !
  type, public :: vandermonde_matrix
    type(matrix), allocatable :: nodal2modal
    type(matrix), allocatable :: modal2nodal
  end type vandermonde_matrix
  !
  ! vand : projection matrix for all combinations of cell geometry and 
  !        and cell order that are possible in the current simulation
  !
  !     For cell 'n' :
  !
  !         The expression for converting the nodal solution in cell 'n' to
  !         its modal form is
  !
  !     matmul( vand(this_geom,this_order)%nodal2modal%mat , nodal_solution )
  !
  !         The expression for converting the modal solution in cell 'n' to
  !         its nodal form is
  !
  !     matmul( vand(this_geom,this_order)%modal2nodal%mat , modal_solution )
  !
  !     Where :
  !        nsp = # of solution points within this cell geometry determined
  !              by the order of the cell
  !        this_geom = cell_geom(n)
  !        this_order = cell_order(n)
  !        mat has dimensions and extents mat(1:nsp,1:nsp)
  !        nodal_solution(1:nsp) : nodal solution vector for cell 'n'
  !        modal_solution(1:nsp) : modal solution vector for cell 'n'
  !
  type(vandermonde_matrix), public, save, target, allocatable :: vand(:,:)
  !
end module vandermonde_mod
!
!###############################################################################
!
