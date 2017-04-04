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
  integer(INT_MPI), parameter :: max_int = 1 + MPI_MAX_PROCESSOR_NAME/4
  !
!!!#define DEBUG_MPI_GROUPS
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
  integer :: n,i,i1,i2,ierr,num_hosts
  integer :: cores_per_cpu,host_cpu_num
  integer(INT_MPI) :: name_len
  integer(INT_MPI) :: color,key
  character(len=MPI_MAX_PROCESSOR_NAME) :: host_name
  integer(c_int) :: cpu_num
  integer(INT_MPI) :: host_size,cpu_size
  integer(INT_MPI) :: host_roots_size,cpu_roots_size
  integer(INT_MPI) :: host_roots_rank,cpu_roots_rank
  !
  !.. Automatic Arrays ..
  integer     :: host_int(1:max_int)
  integer     :: all_host_int(1:max_int,1:ncpu)
  logical(lk) :: host_mask(1:ncpu)
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
  call alloc_error(pname,"all_host_num",1,__LINE__,__FILE__, &
                   ierr,error_message,skip_alloc_pause)
  !
  allocate ( all_cpu_num(1:ncpu) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"all_cpu_num",1,__LINE__,__FILE__, &
                   ierr,error_message,skip_alloc_pause)
  !
  call mpi_get_processor_name(host_name,name_len,mpierr)
  !
  my_host_name = trim(adjustl(host_name))
  !
  ! Convert my host name to a set of 32bit integers
  !
  do i = 1,max_int
    i1 = min( (i-1)*4 + 1 , MPI_MAX_PROCESSOR_NAME )
    i2 = min( i*4         , MPI_MAX_PROCESSOR_NAME )
    host_int(i) = transfer(host_name(i1:i2),host_int(i))
  end do
  !
  ! Collect the integer-converted host names of all mpi processes
  !
  call mpi_allgather(    host_int,size(host_int,kind=int_mpi),mpi_inttyp, &
                     all_host_int,size(host_int,kind=int_mpi),mpi_inttyp, &
                     MPI_COMM_WORLD,mpierr)
  !
  ! Convert the integer host names of all the mpi processes back to characters
  !
  do n = 1,ncpu
    do i = 1,max_int
      i1 = min( (i-1)*4 + 1 , MPI_MAX_PROCESSOR_NAME )
      i2 = min( i*4         , MPI_MAX_PROCESSOR_NAME )
      all_host_name(n)(i1:i2)  = transfer(all_host_int(i,n),host_name(i1:i2))
    end do
  end do
  !
  ! Create a mask for all the unique host names
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
  ! Count the number of unique host names
  !
  num_hosts = count(host_mask)
  !
  ! Create a list of the unique host names
  !
  allocate ( unique_hosts(1:num_hosts) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"unique_hosts",1,__LINE__,__FILE__, &
                   ierr,error_message,skip_alloc_pause)
  !
  unique_hosts(:) = pack(all_host_name,host_mask)
  !
  ! Give each unique host a unique integer ID and give
  ! this ID to all processes on the respective host
  !
  do n = 1,num_hosts
    where (all_host_name(:) == unique_hosts(n)) all_host_num(:) = n
    if (host_name == unique_hosts(n)) my_host_num = n
  end do
  !
  ! Get the CPU number of each process local to each host
  !
  cpu_num = findmycpu()
  my_cpu_num = int(cpu_num,kind=kind(my_cpu_num))
  !
  ! Collect the local CPU numbers for all MPI processes
  !
  call mpi_allgather( my_cpu_num,1_int_mpi,mpi_inttyp, &
                     all_cpu_num,1_int_mpi,mpi_inttyp, &
                     MPI_COMM_WORLD,mpierr)
  !
  cores_per_cpu = (maxval(all_cpu_num)+1) / 2
  !
 !n = count(all_host_num == my_host_num)
 !!
 !allocate ( my_host_grp(1:n) , source=0_int_mpi , &
 !           stat=ierr , errmsg=error_message )
 !call alloc_error(pname,"my_host_grp",1,__LINE__,__FILE__,ierr,error_message)
 !!
 !my_host_grp(:) = int( pack(intseq(0,ncpu-1),all_host_num==my_host_num) , &
 !                        kind=int_mpi )
 !!
 !allocate ( my_cpu_grp(1:cores_per_cpu) , source=0_int_mpi , &
 !           stat=ierr , errmsg=error_message )
 !call alloc_error(pname,"my_cpu_grp",1,__LINE__,__FILE__,ierr,error_message)
 !!
 !i = 0
 !do n = 1,size(my_host_grp)
 !  host_cpu_num = all_cpu_num(my_host_grp(n)+1)
 !  if (host_cpu_num/cores_per_cpu == my_cpu_num/cores_per_cpu) then
 !    i = i + 1
 !    my_cpu_grp(i) = my_host_grp(n)
 !  end if
 !end do
  !
  ! Create new communicators for all processes on each unique host
  !
  color = int(my_host_num,kind=INT_MPI)
  key = int(my_cpu_num,kind=INT_MPI)
  call mpi_comm_split(MPI_COMM_WORLD,color,key,my_host_comm,mpierr)
  !
  ! Get the total number of processes and the local rank
  ! of this process within the my_host_comm communicator
  !
  call mpi_comm_size(my_host_comm,host_size,mpierr)
  call mpi_comm_rank(my_host_comm,my_host_rank,mpierr)
  !
  allocate ( my_host_grp(1:host_size) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"my_host_grp",1,__LINE__,__FILE__, &
                   ierr,error_message,skip_alloc_pause)
  !
  ! Collect the global ranks of all the MPI processes included within
  ! the my_host_comm communicator that the current process is a part of
  !
  call mpi_allgather(mypnum,1_int_mpi,MPI_INTEGER, &
                     my_host_grp,1_int_mpi,MPI_INTEGER, &
                     my_host_comm,mpierr)
  !
  ! Create new communicators for all processes on each physical cpu
  !
  color = int( 2*(my_host_num-1) + my_cpu_num/cores_per_cpu , kind=INT_MPI )
  key = int(my_cpu_num,kind=INT_MPI)
  call mpi_comm_split(MPI_COMM_WORLD,color,key,my_cpu_comm,mpierr)
  !
  ! Get the total number of processes and the local rank
  ! of this process within the my_cpu_comm communicator
  !
  call mpi_comm_size(my_cpu_comm,cpu_size,mpierr)
  call mpi_comm_rank(my_cpu_comm,my_cpu_rank,mpierr)
  !
  allocate ( my_cpu_grp(1:cpu_size) , source=0_int_mpi , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"my_cpu_grp",1,__LINE__,__FILE__, &
                   ierr,error_message,skip_alloc_pause)
  !
  ! Collect the global ranks of all the MPI processes included within
  ! the my_cpu_comm communicator that the current process is a part of
  !
  call mpi_allgather(mypnum,1_int_mpi,MPI_INTEGER, &
                     my_cpu_grp,1_int_mpi,MPI_INTEGER, &
                     my_cpu_comm,mpierr)
  !
  ! Create a new communicator for all the root processes on each unique host
  !
  if (my_cpu_num == 0) then
    color = int(my_cpu_num,kind=INT_MPI)
    i_am_host_root = true
  else
    color = MPI_UNDEFINED
    i_am_host_root = fals
  end if
  key = mypnum
  call mpi_comm_split(MPI_COMM_WORLD,color,key,host_roots_comm,mpierr)
  !
  ! If the current process is the root process for a host (e.g. NAS Node),
  ! get the total number of processes and the local rank
  ! of this process within the host_roots_comm communicator
  !
  if (i_am_host_root) then
    !
    my_host_root = mypnum
    !
    call mpi_comm_size(host_roots_comm,host_roots_size,mpierr)
    call mpi_comm_rank(host_roots_comm,host_roots_rank,mpierr)
    !
    allocate ( host_roots_grp(1:host_roots_size) , source=0_int_mpi , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"host_roots_grp",1,__LINE__,__FILE__, &
                     ierr,error_message,skip_alloc_pause)
    !
    call mpi_gather(mypnum,1_int_mpi,MPI_INTEGER, &
                    host_roots_grp,1_int_mpi,MPI_INTEGER, &
                    glb_root,host_roots_comm,mpierr)
    !
  end if
  !
  ! Broadcast the host_roots information to all processors
  !
  call mpi_bcast(host_roots_size,1_int_mpi,MPI_INTEGER, &
                 glb_root,MPI_COMM_WORLD,mpierr)
  !
  if (.not.allocated(host_roots_grp)) then
    allocate ( host_roots_grp(1:host_roots_size) , source=0_int_mpi , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"host_roots_grp",1,__LINE__,__FILE__, &
                     ierr,error_message,skip_alloc_pause)
  end if
  !
  call mpi_bcast(host_roots_grp,host_roots_size,MPI_INTEGER, &
                 glb_root,MPI_COMM_WORLD,mpierr)
  !
  ! Create a new communicator for all the root processes on each physical CPU
  !
  if (any(my_cpu_num == [0,cores_per_cpu])) then
    color = 1_int_mpi
    i_am_cpu_root = true
  else
    color = MPI_UNDEFINED
    i_am_cpu_root = fals
  end if
  key = mypnum
  call mpi_comm_split(MPI_COMM_WORLD,color,key,cpu_roots_comm,mpierr)
  !
  ! If the current process is the root process for a physical
  ! cpu, get the total number of processes and the local rank
  ! of this process within the cpu_roots_comm communicator
  !
  if (i_am_cpu_root) then
    !
    my_cpu_root = mypnum
    !
    call mpi_comm_size(cpu_roots_comm,cpu_roots_size,mpierr)
    call mpi_comm_rank(cpu_roots_comm,cpu_roots_rank,mpierr)
    !
    allocate ( cpu_roots_grp(1:cpu_roots_size) , source=0_int_mpi , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cpu_roots_grp",1,__LINE__,__FILE__, &
                     ierr,error_message,skip_alloc_pause)
    !
    call mpi_gather(mypnum,1_int_mpi,MPI_INTEGER, &
                    cpu_roots_grp,1_int_mpi,MPI_INTEGER, &
                    glb_root,cpu_roots_comm,mpierr)
    !
  end if
  !
  ! Broadcast the cpu_roots information to all processors
  !
  call mpi_bcast(cpu_roots_size,1_int_mpi,MPI_INTEGER, &
                 glb_root,MPI_COMM_WORLD,mpierr)
  !
  if (.not.allocated(cpu_roots_grp)) then
    allocate ( cpu_roots_grp(1:cpu_roots_size) , source=0_int_mpi , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cpu_roots_grp",1,__LINE__,__FILE__, &
                     ierr,error_message,skip_alloc_pause)
  end if
  !
  call mpi_bcast(cpu_roots_grp,cpu_roots_size,MPI_INTEGER, &
                 glb_root,MPI_COMM_WORLD,mpierr)
  !
  ! Broadcast the respective host root to all processors in each host group
  !
  call mpi_bcast(my_host_root,1_int_mpi,MPI_INTEGER, &
                 host_root,my_host_comm,mpierr)
  !
  ! Broadcast the respective cpu root to all processors in each cpu group
  !
  call mpi_bcast(my_cpu_root,1_int_mpi,MPI_INTEGER, &
                 cpu_root,my_cpu_comm,mpierr)
  !
#ifdef DEBUG_MPI_GROUPS
 !if (mypnum == 0) then
 !  do n = 1,size(all_cpu_num)
 !    write (150,1) n,all_cpu_num(n)
 !  end do
 !  flush (150)
 !  do n = 1,size(unique_hosts)
 !    write (151,2) n,trim(adjustl(unique_hosts(n)))
 !  end do
 !  flush (151)
 !  do n = 1,size(all_host_num)
 !    write (152,3) n,all_host_num(n)
 !  end do
 !  flush (152)
 !  do n = 1,size(all_host_name)
 !    write (153,4) n,trim(adjustl(all_host_name(n)))
 !  end do
 !  flush (153)
 !end if
 !!
  1 format ("all_cpu_num(",i5,") = ",i0)
  2 format ("unique_hosts(",i5,") = '",a,"'")
  3 format ("all_host_num(",i5,") = ",i0)
  4 format ("all_host_name(",i5,") = '",a,"'")
  5 format ("mypnum = ",i5,"    my_cpu_num = ",i5,"    host_name = '",a,"'")
  6 format ("GLBCPU = ",i3,"; ",a," = ",i2,"  :  ",a," =",*(1x,i3,:))
  7 format (31x,a," =",*(1x,i3,:))
  8 format (14x,a," = ",i2,"     ",a," =",*(1x,i3,:))
  9 format ("GLBCPU = ",i3,"; ",12x,"  :  ",a," =",*(1x,i3,:))
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  do n = 1,ncpu
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == n-1) then
     !write (iout,5) mypnum,my_cpu_num,trim(adjustl(host_name))
      if (i_am_host_root) then
        write (iout,6) mypnum,"HR-RANK",host_roots_rank, &
                   "host_roots_grp",(host_roots_grp(i),i=1,size(host_roots_grp))
        if (i_am_cpu_root) then
          write (iout,8) "CR-RANK",cpu_roots_rank, & 
                   " cpu_roots_grp",(cpu_roots_grp(i),i=1,size(cpu_roots_grp))
        end if
        write (iout,8) "LH-RANK",my_host_rank, &
                   "   my_host_grp",(my_host_grp(i),i=1,size(my_host_grp))
        write (iout,8) "LC-RANK",my_cpu_rank, &
                   "    my_cpu_grp",(my_cpu_grp(i),i=1,size(my_cpu_grp))
      else if (i_am_cpu_root) then
        write (iout,9) mypnum, &
                   "host_roots_grp",(host_roots_grp(i),i=1,size(host_roots_grp))
        write (iout,8) "CR-RANK",cpu_roots_rank, &
                   " cpu_roots_grp",(cpu_roots_grp(i),i=1,size(cpu_roots_grp))
        write (iout,8) "LH-RANK",my_host_rank, &
                   "   my_host_grp",(my_host_grp(i),i=1,size(my_host_grp))
        write (iout,8) "LC-RANK",my_cpu_rank, &
                   "    my_cpu_grp",(my_cpu_grp(i),i=1,size(my_cpu_grp))
      else
        write (iout,9) mypnum, &
                   "host_roots_grp",(host_roots_grp(i),i=1,size(host_roots_grp))
        write (iout,7) &
                   " cpu_roots_grp",(cpu_roots_grp(i),i=1,size(cpu_roots_grp))
        write (iout,8) "LH-RANK",my_host_rank, &
                   "   my_host_grp",(my_host_grp(i),i=1,size(my_host_grp))
        write (iout,8) "LC-RANK",my_cpu_rank, &
                   "    my_cpu_grp",(my_cpu_grp(i),i=1,size(my_cpu_grp))
      end if
      flush (iout)
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end do
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  call stop_recon(stop_mpi,pname,__LINE__,__FILE__, &
                  "AFTER DUMPING CPU LOCALITY!!")
#endif
  !
end subroutine get_cpu_locality
!
!###############################################################################
!
end module cpu_info_mod
