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
!!!#define DEBUG_COMMUNICATORS
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
  integer :: l,n,i,i1,i2,ierr,num_hosts
  integer :: cores_per_cpu,host_cpu_num
  integer :: num_cpus,max_cpus
  integer :: min_cpu_num,my_new_cpu_num
  integer(INT_MPI) :: name_len
  integer(INT_MPI) :: color,key
  character(len=MPI_MAX_PROCESSOR_NAME) :: host_name
  integer(c_int) :: cpu_num,num_cpus_c,max_cpus_c
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
  type :: host_t
    integer, allocatable :: mpi_num(:)
    integer, allocatable :: cpu_num(:)
    integer, allocatable :: new_cpu_num(:)
  end type host_t
  !
  !.. Allocatable Arrays ..
  logical(lk), allocatable :: msk(:)
 !logical(lk), allocatable :: mpi_num_mask(:)
  !
  character(len=MPI_MAX_PROCESSOR_NAME), allocatable :: unique_hosts(:)
  !
  type(host_t), allocatable :: host(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_cpu_locality"
  !
  !.. C Function Interfaces ..
  interface
    !
    function findmycpu() bind(c,name="findmycpu_")
      use, intrinsic :: iso_c_binding, only : c_int
      integer(c_int) :: findmycpu
    end function findmycpu
    !
   !function get_num_cpus_per_host() bind(c,name="get_num_cpus_per_host_")
   !  use, intrinsic :: iso_c_binding, only : c_int
   !  integer(c_int) :: get_num_cpus_per_host
   !end function get_num_cpus_per_host
   !!
   !function get_max_cpus_per_host() bind(c,name="get_max_cpus_per_host_")
   !  use, intrinsic :: iso_c_binding, only : c_int
   !  integer(c_int) :: get_max_cpus_per_host
   !end function get_max_cpus_per_host
   !!
  end interface
  !
continue
  !
  ! Initialize communication variables
  !
  my_host_comm = MPI_COMM_WORLD
  my_cpu_comm = MPI_COMM_WORLD
  host_roots_comm = MPI_COMM_WORLD
  cpu_roots_comm = MPI_COMM_WORLD
  !
  if (ncpu == 1) return
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
  100 format (5x,a," = '",a,"'",:,t60,a,":",i0)
  101 format (5x,a," = ",i0,:,t60,a,":",i0)
  102 format (5x,a,"(",i0,") = '",a,"'",t60,a,":",i0)
  103 format (5x,a,"(",i0,") = ",i0,t60,a,":",i0)
  104 format (5x,a," = ",l1,:,t60,a,":",i0)
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
    if (my_host_name == unique_hosts(n)) my_host_num = n
  end do
  !
  ! Get the CPU number of each process local to each host
  !
 !write (*,*) "b"
#ifdef PBS_ENV
  cpu_num = findmycpu()
  my_cpu_num = int(cpu_num,kind=kind(my_cpu_num))
#else
  my_cpu_num = int(mypnum,kind=kind(my_cpu_num))
#endif
#ifdef DEBUG_COMMUNICATORS
  write (iout,100) "my_host_name",trim(adjustl(my_host_name)),__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
 !!
 !num_cpus_c = get_num_cpus_per_host()
 !num_cpus = int(num_cpus_c,kind=kind(num_cpus))
#ifdef DEBUG_COMMUNICATORS
 !write (iout,101) "num_cpus",num_cpus,__FILE__,__LINE__
 !flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
 !!
 !max_cpus_c = get_max_cpus_per_host()
 !max_cpus = int(max_cpus_c,kind=kind(max_cpus))
#ifdef DEBUG_COMMUNICATORS
 !write (iout,101) "max_cpus",max_cpus,__FILE__,__LINE__
 !flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
  ! Collect the local CPU numbers for all MPI processes
  !
  call mpi_allgather( my_cpu_num,1_int_mpi,mpi_inttyp, &
                     all_cpu_num,1_int_mpi,mpi_inttyp, &
                     MPI_COMM_WORLD,mpierr)
  !
#ifdef DEBUG_COMMUNICATORS
  write (iout,103) ("all_cpu_num",n,all_cpu_num(n), &
                    __FILE__,__LINE__,n=1,size(all_cpu_num))
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
  ! Create a list of the MPI ranks and CPU numbers on each host
  !
 !allocate ( mpi_num_mask(1:ncpu) , source=fals , &
 !           stat=ierr , errmsg=error_message )
 !call alloc_error(pname,"mpi_num_mask",1,__LINE__,__FILE__, &
 !                 ierr,error_message,skip_alloc_pause)
  !
  allocate ( host(1:num_hosts) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"host",1,__LINE__,__FILE__,ierr,error_message)
  !
  do n = 1,num_hosts
   !!
   !!######################################################################
   !!                              METHOD 1
   !!----------------------------------------------------------------------
   !!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   !!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
   !!
   !mpi_num_mask(:) = fals
   !do i = 1,size(all_host_num)
   !  if (all_host_num(i) == n) then
   !    mpi_num_mask(i) = true
   !  end if
   !end do
   !!
   !! F2003 AUTO-REALLOCATION
   !host(n)%mpi_num = pack( intseq(1,ncpu) , mask=mpi_num_mask )
   !!
   !!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
   !!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   !!----------------------------------------------------------------------
   !!######################################################################
   !!
   !!######################################################################
   !!                              METHOD 2
   !!----------------------------------------------------------------------
   !!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   !!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
   !!
   !mpi_num_mask(:) = ( all_host_num(:) == n )
   !!
   !! F2003 AUTO-REALLOCATION
   !host(n)%mpi_num = pack( intseq(1,ncpu) , mask=mpi_num_mask )
   !!
   !!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
   !!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   !!----------------------------------------------------------------------
   !!######################################################################
    !
    !######################################################################
    !                              METHOD 3
    !----------------------------------------------------------------------
    !||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    !VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
    !
    ! F2003 AUTO-REALLOCATION
    ! NOTE: This will make the numbering in mpi_num go from 1 to ncpu
    host(n)%mpi_num = pack( intseq(1,ncpu) , mask=(all_host_num(:)==n) )
    !
    !AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    !||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    !----------------------------------------------------------------------
    !######################################################################
    !
    ! F2003 AUTO-REALLOCATION
    host(n)%cpu_num = all_cpu_num( host(n)%mpi_num )
    !
    ! Now subtract 1 from host(:)%mpi_num so the numbering
    ! goes from 0 to ncpu-1 and matches the values in mypnum
    !
    host(n)%mpi_num(:) = host(n)%mpi_num(:) - 1
    !
    ! Now create new CPU numbering that goes from 0 to size(host(:)%cpu_num)-1
    !
    ! F2003 AUTO-REALLOCATION
    host(n)%new_cpu_num = host(n)%cpu_num
    !
    allocate ( msk(1:size(host(n)%cpu_num)) , source=true , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"msk",1,__LINE__,__FILE__,ierr,error_message)
    !
    l = 0
    do i = 1,size(host(n)%cpu_num)
      min_cpu_num = minloc(host(n)%cpu_num,mask=msk,dim=1)
      host(n)%new_cpu_num(min_cpu_num) = l
      l = l + 1
      msk(min_cpu_num) = fals
    end do
    !
    deallocate ( msk , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
    !
  end do
  !
  ! Find the new CPU number for this MPI process
  !
  do i = 1,size(host(my_host_num)%mpi_num)
    if (mypnum == host(my_host_num)%mpi_num(i)) then
      my_new_cpu_num = host(my_host_num)%new_cpu_num(i)
      exit
    end if
  end do
  !
  ! Create new communicators for all processes on each unique host
  !
  color = int(my_host_num,kind=INT_MPI)
#ifdef DEBUG_COMMUNICATORS
  write (iout,101) "mpi_comm_split(host,color)",color,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  key = int(my_new_cpu_num,kind=INT_MPI)
 !key = int(my_cpu_num,kind=INT_MPI)
#ifdef DEBUG_COMMUNICATORS
  write (iout,101) "mpi_comm_split(host,key)",key,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  call mpi_comm_split(MPI_COMM_WORLD,color,key,my_host_comm,mpierr)
#ifdef DEBUG_COMMUNICATORS
  write (iout,101) "mpi_comm_split(host,mpierr)",mpierr,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
  ! Get the total number of processes and the local rank
  ! of this process within the my_host_comm communicator
  !
  call mpi_comm_size(my_host_comm,host_size,mpierr)
#ifdef DEBUG_COMMUNICATORS
  write (iout,101) "host_size",host_size,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  call mpi_comm_rank(my_host_comm,my_host_rank,mpierr)
#ifdef DEBUG_COMMUNICATORS
  write (iout,101) "my_host_rank",my_host_rank,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
#ifdef PBS_ENV
 !cores_per_cpu = max( 1 , (maxval(all_cpu_num)+1) / 2 )
 !cores_per_cpu = max( 1 , size(host(my_host_num)%cpu_num) / 2 )
  cores_per_cpu = max( 1 , host_size/2 )
#else
 !cores_per_cpu = maxval(all_cpu_num)+1
 !cores_per_cpu = max( 1 , size(host(my_host_num)%cpu_num) )
  cores_per_cpu = max( 1 , host_size )
#endif
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
#ifdef DEBUG_COMMUNICATORS
  write (iout,103) ("my_host_grp",n,my_host_grp(n), &
                    __FILE__,__LINE__,n=1,size(my_host_grp))
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
  ! Create new communicators for all processes on each physical cpu
  !
#ifdef DEBUG_COMMUNICATORS
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  do n = 1,ncpu
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == n-1) then
      write (iout,101) "my_host_num",my_host_num,__FILE__,__LINE__
      write (iout,101) "my_cpu_num",my_cpu_num,__FILE__,__LINE__
      write (iout,101) "my_new_cpu_num",my_new_cpu_num,__FILE__,__LINE__
      write (iout,101) "cores_per_cpu",cores_per_cpu,__FILE__,__LINE__
      flush (iout)
    end if
    flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end do
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
 !write (*,*) "c"
  color = int( 2*(my_host_num-1) + my_host_rank/cores_per_cpu , kind=INT_MPI )
 !color = int( 2*(my_host_num-1) + my_new_cpu_num/cores_per_cpu , kind=INT_MPI )
#ifdef DEBUG_COMMUNICATORS
  write (iout,101) "mpi_comm_split(cpu,color)",color,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  key = int(my_new_cpu_num,kind=INT_MPI)
 !key = int(my_cpu_num,kind=INT_MPI)
#ifdef DEBUG_COMMUNICATORS
  write (iout,101) "mpi_comm_split(cpu,key)",key,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  call mpi_comm_split(MPI_COMM_WORLD,color,key,my_cpu_comm,mpierr)
#ifdef DEBUG_COMMUNICATORS
  write (iout,101) "mpi_comm_split(cpu,mpierr)",mpierr,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
  ! Get the total number of processes and the local rank
  ! of this process within the my_cpu_comm communicator
  !
  call mpi_comm_size(my_cpu_comm,cpu_size,mpierr)
#ifdef DEBUG_COMMUNICATORS
  write (iout,101) "cpu_size",cpu_size,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  call mpi_comm_rank(my_cpu_comm,my_cpu_rank,mpierr)
#ifdef DEBUG_COMMUNICATORS
  write (iout,101) "my_cpu_rank",my_cpu_rank,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
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
#ifdef DEBUG_COMMUNICATORS
  write (iout,103) ("my_cpu_grp",n,my_cpu_grp(n), &
                    __FILE__,__LINE__,n=1,size(my_cpu_grp))
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
  ! Create a new communicator for all the root processes on each unique host
  !
 !if (my_cpu_num == 0) then
 !  color = int(my_cpu_num,kind=INT_MPI)
  if (my_host_rank == 0) then
    color = int(my_host_rank,kind=INT_MPI)
    i_am_host_root = true
  else
    color = MPI_UNDEFINED
    i_am_host_root = fals
  end if
#ifdef DEBUG_COMMUNICATORS
  write (iout,104) "i_am_host_root",i_am_host_root,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
#ifdef DEBUG_COMMUNICATORS
  write (iout,101) "mpi_comm_split(host_roots,color)",color,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  key = mypnum
#ifdef DEBUG_COMMUNICATORS
  write (iout,101) "mpi_comm_split(host_roots,key)",key,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  call mpi_comm_split(MPI_COMM_WORLD,color,key,host_roots_comm,mpierr)
#ifdef DEBUG_COMMUNICATORS
  write (iout,101) "mpi_comm_split(host_roots,mpierr)",mpierr,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
  ! If the current process is the root process for a host (e.g. NAS Node),
  ! get the total number of processes and the local rank
  ! of this process within the host_roots_comm communicator
  !
 !write (*,*) "d"
  if (i_am_host_root) then
    !
    my_host_root = mypnum
#ifdef DEBUG_COMMUNICATORS
    write (iout,101) "my_host_root",my_host_root,__FILE__,__LINE__
    flush (iout) ; call mpi_barrier(host_roots_comm,mpierr)
#endif
    !
    call mpi_comm_size(host_roots_comm,host_roots_size,mpierr)
#ifdef DEBUG_COMMUNICATORS
    write (iout,101) "host_roots_size",host_roots_size,__FILE__,__LINE__
    flush (iout) ; call mpi_barrier(host_roots_comm,mpierr)
#endif
    call mpi_comm_rank(host_roots_comm,host_roots_rank,mpierr)
#ifdef DEBUG_COMMUNICATORS
    write (iout,101) "host_roots_rank",host_roots_rank,__FILE__,__LINE__
    flush (iout) ; call mpi_barrier(host_roots_comm,mpierr)
#endif
    !
    allocate ( host_roots_grp(1:host_roots_size) , source=0_int_mpi , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"host_roots_grp",1,__LINE__,__FILE__, &
                     ierr,error_message,skip_alloc_pause)
    !
    call mpi_gather(mypnum,1_int_mpi,MPI_INTEGER, &
                    host_roots_grp,1_int_mpi,MPI_INTEGER, &
                    glb_root,host_roots_comm,mpierr)
#ifdef DEBUG_COMMUNICATORS
    write (iout,103) ("host_roots_grp",n,host_roots_grp(n),__FILE__,__LINE__, &
                      n=1,size(host_roots_grp))
    flush (iout) ; call mpi_barrier(host_roots_comm,mpierr)
#endif
    !
  end if
  !
  ! Broadcast the host_roots information to all processors
  !
  call mpi_bcast(host_roots_size,1_int_mpi,MPI_INTEGER, &
                 glb_root,MPI_COMM_WORLD,mpierr)
#ifdef DEBUG_COMMUNICATORS
  write (iout,101) "host_roots_size",host_roots_size,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
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
#ifdef DEBUG_COMMUNICATORS
  write (iout,103) ("host_roots_grp",n,host_roots_grp(n),__FILE__,__LINE__, &
                    n=1,size(host_roots_grp))
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
  ! Create a new communicator for all the root processes on each physical CPU
  !
 !if (any(my_cpu_num == [0,cores_per_cpu])) then
  if (any(my_cpu_rank == [0,cores_per_cpu])) then
    color = 1_int_mpi
    i_am_cpu_root = true
  else
    color = MPI_UNDEFINED
    i_am_cpu_root = fals
  end if
#ifdef DEBUG_COMMUNICATORS
  write (iout,104) "i_am_cpu_root",i_am_cpu_root,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
#ifdef DEBUG_COMMUNICATORS
  write (iout,101) "mpi_comm_split(cpu_roots,color)",color,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  key = mypnum
#ifdef DEBUG_COMMUNICATORS
  write (iout,101) "mpi_comm_split(cpu_roots,key)",key,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  call mpi_comm_split(MPI_COMM_WORLD,color,key,cpu_roots_comm,mpierr)
#ifdef DEBUG_COMMUNICATORS
  write (iout,101) "mpi_comm_split(cpu_roots,mpierr)",mpierr,__FILE__,__LINE__
  flush (iout) ; call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
  ! If the current process is the root process for a physical
  ! cpu, get the total number of processes and the local rank
  ! of this process within the cpu_roots_comm communicator
  !
 !write (*,*) "e"
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
 !write (*,*) "f"
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
  6 format ("GLBCPU = ",i3,"; ",a," = ",i2,"  :  ",a," =",100(1x,i3,:))
 !6 format ("GLBCPU = ",i3,"; ",a," = ",i2,"  :  ",a," =",*(1x,i3,:))
  7 format (31x,a," =",100(1x,i3,:))
 !7 format (31x,a," =",*(1x,i3,:))
  8 format (14x,a," = ",i2,"     ",a," =",100(1x,i3,:))
 !8 format (14x,a," = ",i2,"     ",a," =",*(1x,i3,:))
  9 format ("GLBCPU = ",i3,"; ",12x,"  :  ",a," =",100(1x,i3,:))
 !9 format ("GLBCPU = ",i3,"; ",12x,"  :  ",a," =",*(1x,i3,:))
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
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__, &
                  "AFTER DUMPING CPU LOCALITY!!")
#endif
  !
end subroutine get_cpu_locality
!
!###############################################################################
!
end module cpu_info_mod
