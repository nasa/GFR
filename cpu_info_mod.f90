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
  integer :: num_cpus,max_cpus,num_cores
  integer :: min_cpu_num,my_new_cpu_num
  integer(INT_MPI) :: name_len
  integer(INT_MPI) :: color,key
  character(len=MPI_MAX_PROCESSOR_NAME) :: host_name
  integer(c_int) :: cpu_num,num_cores_c,max_cores_c
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
    integer              :: num_cpus
    integer              :: num_cores
    integer, allocatable :: mpi_num(:)
    integer, allocatable :: cpu_num(:)
    integer, allocatable :: new_cpu_num(:)
  end type host_t
  !
  !.. Allocatable Arrays ..
  logical(lk), allocatable :: msk(:)
  !
  character(len=MPI_MAX_PROCESSOR_NAME), allocatable :: unique_hosts(:)
  !
  type(host_t), allocatable :: host(:)
  !
  integer, allocatable :: host_cores(:)
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
    function get_num_cpus_per_host() bind(c,name="get_num_cpus_per_host_")
      use, intrinsic :: iso_c_binding, only : c_int
      integer(c_int) :: get_num_cpus_per_host
    end function get_num_cpus_per_host
    !
    function get_max_cpus_per_host() bind(c,name="get_max_cpus_per_host_")
      use, intrinsic :: iso_c_binding, only : c_int
      integer(c_int) :: get_max_cpus_per_host
    end function get_max_cpus_per_host
    !
  end interface
  !
continue
  !
  call debug_timer(entering_procedure,pname)
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
  ! Get the host name for each processor
  !
  all_host_num = init_array(ncpu,0) ! F2003 AUTO-REALLOC
  all_cpu_num  = init_array(ncpu,0) ! F2003 AUTO-REALLOC
  !
  ! Get the host name for the current processor
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
 !!
 !! Count the number of unique host names
 !!
 !num_hosts = count(host_mask)
 !!
 !! Create a list of the unique host names
 !!
 !allocate ( unique_hosts(1:num_hosts) , stat=ierr , errmsg=error_message )
 !call alloc_error(pname,"unique_hosts",1,__LINE__,__FILE__, &
 !                 ierr,error_message,skip_alloc_pause)
 !!
 !unique_hosts(:) = pack(all_host_name,host_mask)
  !
  ! Create a list of the unique host names
  !
  unique_hosts = pack( all_host_name(:) , host_mask(:) ) ! F2003 AUTO-REALLOC
  num_hosts = size(unique_hosts)
  !
  ! Give each unique host a unique integer ID and give
  ! this ID to all processes on the respective host
  !
  do n = 1,num_hosts
    where (all_host_name(:) == unique_hosts(n)) all_host_num(:) = n
    if (my_host_name == unique_hosts(n)) my_host_num = n
  end do
  !
  !#####################################################################
  !#####################################################################
  !
  ! Get the CPU number of each process local to each host
  !
  !#####################################################################
  !#####################################################################
  !
#ifdef PBS_ENV
  cpu_num = findmycpu()
  my_cpu_num = int(cpu_num,kind=kind(my_cpu_num))
#else
  my_cpu_num = int(mypnum,kind=kind(my_cpu_num))
#endif
  !
  ! Collect the local CPU numbers for all MPI processes
  !
  call mpi_allgather( my_cpu_num,1_int_mpi,mpi_inttyp, &
                     all_cpu_num,1_int_mpi,mpi_inttyp, &
                     MPI_COMM_WORLD,mpierr)
  !
  !
  !#####################################################################
  !#####################################################################
  !
  ! Create a list of the MPI ranks and CPU numbers on each host
  !
  !#####################################################################
  !#####################################################################
  !
  allocate ( host(1:num_hosts) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"host",1,__LINE__,__FILE__,ierr,error_message)
  !
  do n = 1,num_hosts
    !
#ifdef PBS_ENV
    host(n)%num_cpus = 2
#else
    host(n)%num_cpus = 1
#endif
    !
    ! F2003 AUTO-REALLOCATION
    ! NOTE: This will make the numbering in mpi_num go from 1 to ncpu
    host(n)%mpi_num = pack( intseq(1,ncpu) , mask=(all_host_num(:)==n) )
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
    msk = init_array(size(host(n)%cpu_num),true) ! F2003 AUTO-REALLOC
    !
    l = 0
    do i = 1,size(host(n)%cpu_num)
      min_cpu_num = minloc(host(n)%cpu_num,mask=msk,dim=1)
      host(n)%new_cpu_num(min_cpu_num) = l
      l = l + 1
      msk(min_cpu_num) = fals
    end do
    !
  end do
  !
  if (allocated(msk)) then
    deallocate ( msk , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
  end if
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
  !#####################################################################
  !#####################################################################
  !
  ! Create new communicators for all processes on each unique host
  !
  !#####################################################################
  !#####################################################################
  !
  color = int(my_host_num,kind=INT_MPI)
  key = int(my_new_cpu_num,kind=INT_MPI)
 !key = int(my_cpu_num,kind=INT_MPI)
  call mpi_comm_split(MPI_COMM_WORLD,color,key,my_host_comm,mpierr)
  !
  ! Get the total number of processes and the local rank
  ! of this process within the my_host_comm communicator
  !
  call mpi_comm_size(my_host_comm,host_size,mpierr)
  call mpi_comm_rank(my_host_comm,my_host_rank,mpierr)
  !
  my_host_grp = init_array(host_size,0_int_mpi) ! F2003 AUTO-REALLOC
  !
  ! Collect the global ranks of all the MPI processes included within
  ! the my_host_comm communicator that the current process is a part of
  !
  call mpi_allgather(mypnum,1_int_mpi,MPI_INTEGER, &
                     my_host_grp,1_int_mpi,MPI_INTEGER, &
                     my_host_comm,mpierr)
  !
  !#####################################################################
  !#####################################################################
  !
  ! Create a new communicator for all the root processes on each unique host
  !
  !#####################################################################
  !#####################################################################
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
    host_roots_grp = init_array(host_roots_size,0_int_mpi) ! F2003 AUTO-REALLOC
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
    host_roots_grp = init_array(host_roots_size,0_int_mpi) ! F2003 AUTO-REALLOC
  end if
  !
  call mpi_bcast(host_roots_grp,host_roots_size,MPI_INTEGER, &
                 glb_root,MPI_COMM_WORLD,mpierr)
  !
  !#####################################################################
  !#####################################################################
  !
  ! Get the number of processor cores on each host
  !
  !#####################################################################
  !#####################################################################
  !
  host_cores = init_array(num_hosts,0) ! F2003 AUTO-REALLOC
  !
  do n = 1,num_hosts
    if (n == my_host_num) then
      !
      num_cores_c = get_num_cpus_per_host()
      ! get_num_cpus_per_host returns number of hyper-threading cores,
      ! so divide this by 2 since there are generally 2 threads per core.
      num_cores = int(num_cores_c,kind=kind(num_cores)) / 2
      !
    ! num_cores_c = get_max_cpus_per_host()
    ! num_cores = int(max_cpus_c,kind=kind(num_cores)) / 2
      !
      if (num_cores <= 0) then
       !num_cores = max( 1 , (maxval(all_cpu_num)+1) )
       !num_cores = max( 1 , size(host(my_host_num)%cpu_num) )
        num_cores = max( 1 , host_size )
      end if
      !
      host_cores(n) = num_cores
      !
    end if
  end do
  !
  call mpi_allreduce(MPI_IN_PLACE,host_cores,size(host_cores,kind=int_mpi), &
                     mpi_inttyp,MPI_MAX,MPI_COMM_WORLD,mpierr)
  !
  !
  do n = 1,num_hosts
    host(n)%num_cores = host_cores(n)
  end do
  !
  if (allocated(host_cores)) then
    deallocate ( host_cores , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"host_cores",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  !#####################################################################
  !#####################################################################
  !
  ! Create new communicators for all processes on each physical cpu
  !
  !#####################################################################
  !#####################################################################
  !
  cores_per_cpu = host(my_host_num)%num_cores / host(my_host_num)%num_cpus
  cores_per_cpu = max(1,cores_per_cpu)
  color = int( 2*(my_host_num-1) + my_host_rank/cores_per_cpu , kind=INT_MPI )
  key = int(my_new_cpu_num,kind=INT_MPI)
 !key = int(my_cpu_num,kind=INT_MPI)
  call mpi_comm_split(MPI_COMM_WORLD,color,key,my_cpu_comm,mpierr)
  !
  ! Get the total number of processes and the local rank
  ! of this process within the my_cpu_comm communicator
  !
  call mpi_comm_size(my_cpu_comm,cpu_size,mpierr)
  call mpi_comm_rank(my_cpu_comm,my_cpu_rank,mpierr)
  !
  my_cpu_grp = init_array(cpu_size,0_int_mpi) ! F2003 AUTO-REALLOC
  !
  ! Collect the global ranks of all the MPI processes included within
  ! the my_cpu_comm communicator that the current process is a part of
  !
  call mpi_allgather(mypnum,1_int_mpi,MPI_INTEGER, &
                     my_cpu_grp,1_int_mpi,MPI_INTEGER, &
                     my_cpu_comm,mpierr)
  !
  !#####################################################################
  !#####################################################################
  !
  ! Create a new communicator for all the root processes on each physical CPU
  !
  !#####################################################################
  !#####################################################################
  !
  color = MPI_UNDEFINED
  i_am_cpu_root = fals
  !
  min_cpu_num = minloc(host(my_host_num)%cpu_num,dim=1)
  if (mypnum == host(my_host_num)%mpi_num(min_cpu_num)) then
    color = 1_int_mpi
    i_am_cpu_root = true
  end if
  !
  min_cpu_num = minloc( host(my_host_num)%cpu_num , &
                        mask=host(my_host_num)%cpu_num >= cores_per_cpu , &
                        dim=1 )
  if (mypnum == host(my_host_num)%mpi_num(min_cpu_num)) then
    color = 1_int_mpi
    i_am_cpu_root = true
  end if
  !
!!if (any(my_cpu_num == [0,cores_per_cpu])) then
! if (any(my_cpu_rank == [0,cores_per_cpu])) then
!   color = 1_int_mpi
!   i_am_cpu_root = true
! else
!   color = MPI_UNDEFINED
!   i_am_cpu_root = fals
! end if
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
    cpu_roots_grp = init_array(cpu_roots_size,0_int_mpi) ! F2003 AUTO-REALLOC
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
    cpu_roots_grp = init_array(cpu_roots_size,0_int_mpi) ! F2003 AUTO-REALLOC
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
  !#####################################################################
  !#####################################################################
  !
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
 !write (iout,723) mypnum,my_host_name
 !flush (iout)
 !723 format ("CPU-",i0," : my_host_name = ",a)
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
  6 format ("GLBCPU = ",i5,"; ",a," = ",i3," :  ",a," =",1000(1x,i0,:))
 !6 format ("GLBCPU = ",i5,"; ",a," = ",i2,"  :  ",a," =",*(1x,i5,:))
  7 format (33x,a," =",1000(1x,i0,:))
 !7 format (31x,a," =",*(1x,i5,:))
  8 format (14x,a," = ",i5,"    ",a," =",1000(1x,i0,:))
 !8 format (14x,a," = ",i2,"     ",a," =",*(1x,i5,:))
  9 format ("GLBCPU = ",i5,"; ",12x,"  :  ",a," =",1000(1x,i0,:))
 !9 format ("GLBCPU = ",i5,"; ",12x,"  :  ",a," =",*(1x,i5,:))
  11 format (200("#"))
  12 format (200("-"))
  17 format (2x,a10,21x,a," =",1000(1x,i0,:))
  18 format (2x,a10,2x,a," = ",i5,"    ",a," =",1000(1x,i0,:))
  19 format (23x,"MPI_COMM_WORLD  = ",i0,/, &
             23x,"MPI_COMM_NULL   = ",i0,/, &
             23x,"MPI_COMM_SELF   = ",i0,/, &
             23x,"my_host_comm    = ",i0,/, &
             23x,"my_cpu_comm     = ",i0,/, &
             23x,"host_roots_comm = ",i0,/, &
             23x,"cpu_roots_comm  = ",i0)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  do n = 1,ncpu
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == n-1) then
     !write (iout,5) mypnum,my_cpu_num,trim(adjustl(host_name))
      if (i_am_host_root) then
        write (iout,11)
        write (iout,6) mypnum,"HR-RANK",host_roots_rank, &
                   "host_roots_grp",(host_roots_grp(i),i=1,size(host_roots_grp))
        if (i_am_cpu_root) then
          write (iout,18) trim(adjustl(my_host_name)),"CR-RANK",cpu_roots_rank, &
                   " cpu_roots_grp",(cpu_roots_grp(i),i=1,size(cpu_roots_grp))
        end if
        write (iout,8) "LH-RANK",my_host_rank, &
                   "   my_host_grp",(my_host_grp(i),i=1,size(my_host_grp))
        write (iout,8) "LC-RANK",my_cpu_rank, &
                   "    my_cpu_grp",(my_cpu_grp(i),i=1,size(my_cpu_grp))
        write (iout,19) MPI_COMM_WORLD%MPI_VAL, &
                        MPI_COMM_NULL%MPI_VAL, &
                        MPI_COMM_SELF%MPI_VAL, &
                        my_host_comm%MPI_VAL, &
                        my_cpu_comm%MPI_VAL, &
                        host_roots_comm%MPI_VAL, &
                        cpu_roots_comm%MPI_VAL
        write (iout,11)
      else if (i_am_cpu_root) then
        write (iout,12)
        write (iout,9) mypnum, &
                   "host_roots_grp",(host_roots_grp(i),i=1,size(host_roots_grp))
        write (iout,18) trim(adjustl(my_host_name)),"CR-RANK",cpu_roots_rank, &
                   " cpu_roots_grp",(cpu_roots_grp(i),i=1,size(cpu_roots_grp))
        write (iout,8) "LH-RANK",my_host_rank, &
                   "   my_host_grp",(my_host_grp(i),i=1,size(my_host_grp))
        write (iout,8) "LC-RANK",my_cpu_rank, &
                   "    my_cpu_grp",(my_cpu_grp(i),i=1,size(my_cpu_grp))
        write (iout,19) MPI_COMM_WORLD%MPI_VAL, &
                        MPI_COMM_NULL%MPI_VAL, &
                        MPI_COMM_SELF%MPI_VAL, &
                        my_host_comm%MPI_VAL, &
                        my_cpu_comm%MPI_VAL, &
                        host_roots_comm%MPI_VAL, &
                        cpu_roots_comm%MPI_VAL
        write (iout,12)
      else
        write (iout,9) mypnum, &
                   "host_roots_grp",(host_roots_grp(i),i=1,size(host_roots_grp))
        write (iout,17) trim(adjustl(my_host_name)), &
                   " cpu_roots_grp",(cpu_roots_grp(i),i=1,size(cpu_roots_grp))
        write (iout,8) "LH-RANK",my_host_rank, &
                   "   my_host_grp",(my_host_grp(i),i=1,size(my_host_grp))
        write (iout,8) "LC-RANK",my_cpu_rank, &
                   "    my_cpu_grp",(my_cpu_grp(i),i=1,size(my_cpu_grp))
        write (iout,19) MPI_COMM_WORLD%MPI_VAL, &
                        MPI_COMM_NULL%MPI_VAL, &
                        MPI_COMM_SELF%MPI_VAL, &
                        my_host_comm%MPI_VAL, &
                        my_cpu_comm%MPI_VAL, &
                        host_roots_comm%MPI_VAL, &
                        cpu_roots_comm%MPI_VAL
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
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_cpu_locality
!
!###############################################################################
!
end module cpu_info_mod
