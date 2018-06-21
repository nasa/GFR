module octstat_mod
!
! Status codes for performing an octtree based search through xyz data.
!
!   OCT_STATUS_OK              Everything worked fine.
!   OCT_STATUS_MEMORY_ERROR    Error alloc/deallocating memory.
!   OCT_STATUS_NOT_FOUND       Point not found in octtree or octnode.
!   OCT_STATUS_DATA_ERROR      Error evaluating xyz data at a given index.
!
! -------------------------------------------------------------------------

  implicit none
  public

  integer ,parameter :: OCT_STATUS_OK           = 0
  integer ,parameter :: OCT_STATUS_MEMORY_ERROR = 1
  integer ,parameter :: OCT_STATUS_NOT_FOUND    = 2
  integer ,parameter :: OCT_STATUS_DATA_ERROR   = 3

!##########################################################################
end module octstat_mod
