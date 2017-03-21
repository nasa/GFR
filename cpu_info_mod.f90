module cpu_info_mod
  !
  use module_kind_types
  !
  implicit none
  !
  private
  !
  public :: get_cpu_locality
  !
contains
!
!###############################################################################
!
subroutine get_cpu_locality
  !
  !.. Use Statements
  use iso_c_binding, only : c_int
  use ovar, only : all_host_num,all_cpu_num
  !
  !.. Local Scalars ..
  integer :: n,i,ierr
  integer :: num_hosts
  integer(INT_MPI) :: name_len
  character(len=MPI_MAX_PROCESSOR_NAME) :: host_name
  integer(c_int) :: cpu_num
  !
  !.. Automatic Arrays ..
  integer, dimension(1:3) :: host_int
  integer, dimension(1:3,1:ncpu) :: all_host_int
  logical(lk), dimension(1:ncpu) :: host_mask
  !
  character(len=MPI_MAX_PROCESSOR_NAME) :: all_host_name(1:ncpu)
  !
  !.. Allocatable Arrays ..
  character(len=MPI_MAX_PROCESSOR_NAME), allocatable :: unique_hosts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_cpu_locality"
  !
  !.. C Function Interfaces ..
  interface
    function findmycpu() bind(c,name="findmycpu_")
      use, intrinsic :: iso_c_binding, only : c_int
      integer(c_int) :: findmycpu
    end function findmycpu
  end interface
  !
continue
  !
  allocate ( all_host_num(1:ncpu) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"all_host_num",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( all_cpu_num(1:ncpu) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"all_cpu_num",1,__LINE__,__FILE__,ierr,error_message)
  !
  call mpi_get_processor_name(host_name,name_len,mpierr)
  !
  ! Convert my host name to 3 - 32bit integers
  !
  host_int(1) = transfer(host_name(1:4),host_int(1))
  host_int(2) = transfer(host_name(5:8),host_int(2))
  host_int(3) = transfer(host_name(9:12),host_int(3))
  !
  call mpi_allgather(host_int,    3_int_mpi,mpi_inttyp, &
                     all_host_int,3_int_mpi,mpi_inttyp,MPI_COMM_WORLD,mpierr)
  !
  do n = 1,ncpu
    all_host_name(n)(1:4)  = transfer(all_host_int(1,n),host_name(1:4))
    all_host_name(n)(5:8)  = transfer(all_host_int(2,n),host_name(5:8))
    all_host_name(n)(9:12) = transfer(all_host_int(3,n),host_name(9:12))
  end do
  !
  host_mask = true
  !
  do n = 1,ncpu
    !
    if (.not. host_mask(n)) cycle
    !
    do i = n+1,ncpu
      if (host_mask(i)) then
        host_mask(i) = all_host_name(n) /= all_host_name(i)
      end if
    end do
    !
  end do
  !
  num_hosts = count(host_mask)
  !
  allocate ( unique_hosts(1:num_hosts) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"unique_hosts",1,__LINE__,__FILE__,ierr,error_message)
  !
  unique_hosts = pack(all_host_name,host_mask)
  !
  do n = 1,num_hosts
    where (all_host_name(:) == unique_hosts(n)) all_host_num(:) = n
  end do
  !
  cpu_num = findmycpu()
  !
  call mpi_allgather(cpu_num,    1_int_mpi,mpi_inttyp, &
                     all_cpu_num,1_int_mpi,mpi_inttyp,MPI_COMM_WORLD,mpierr)
  !
end subroutine get_cpu_locality
!
!###############################################################################
!
end module cpu_info_mod
