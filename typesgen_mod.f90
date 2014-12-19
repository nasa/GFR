module generic_types_mod
  !
  use module_kind_types, only : wp
  !
  implicit none
  !
  private
  !
  ! matrix : derived type that stores a matrix quantity
  !            Examples : interpolation matrix, derivative matrix, etc.
  type, public :: matrix
    real(wp), allocatable :: mat(:,:)
  end type matrix
  !
end module generic_types_mod
!
!###############################################################################
!
