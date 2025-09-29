module openmp_mod
  !
#ifdef _OPENMP
  use omp_lib, only : openmp_get_num_threads => omp_get_num_threads, &
                      openmp_get_thread_num  => omp_get_thread_num
#endif
  !
  implicit none
  !
  private
  !
  public :: omp_get_num_threads
  public :: omp_get_thread_num
  !
contains
!
!###############################################################################
!
function omp_get_num_threads() result(return_value)
  !
  !.. Function Result ..
  integer :: return_value
  !
continue
  !
#ifdef _OPENMP
  return_value = openmp_get_num_threads()
#else
  return_value = 1
#endif
  !
end function omp_get_num_threads
!
!###############################################################################
!
function omp_get_thread_num() result(return_value)
  !
  !.. Function Result ..
  integer :: return_value
  !
continue
  !
#ifdef _OPENMP
  return_value = openmp_get_thread_num()
#else
  return_value = 0
#endif
  !
end function omp_get_thread_num
!
!###############################################################################
!
end module openmp_mod
