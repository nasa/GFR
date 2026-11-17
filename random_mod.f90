module random_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  implicit none
  !
  private
  !
  public :: set_random_seed
  public :: get_random_seed
  !
  integer, public, allocatable, save :: iseed(:)
  !
contains
!
!###############################################################################
!
subroutine set_random_seed(seed)
  !
  !***********************************************************************
  ! Changes the seed in the random number generator to produce
  ! a different set of numbers
  !***********************************************************************
  !
  !.. Use Statements ..
  use ovar, only : SEM_use_random_seed_from_restart
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: seed
  !
  !.. Local Scalars ..
  integer :: isize,iclk,i,ierr
  !
  !.. Local Arrays ..
  integer :: date(1:8)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "set_random_seed"
  !
continue
  !
  ! If iseed is allocated, that means the random seed was read in from the
  ! restart file. If the input variable SEM_use_random_seed_from_restart is
  ! true, keep iseed allocated and we just need to call random_seed with the
  ! restart seed. If SEM_use_random_seed_from_restart is false, deallocate
  ! iseed so that a new random seed is created.
  !
  if (allocated(iseed)) then
    if (.not. SEM_use_random_seed_from_restart) then
      deallocate ( iseed , stat=ierr, errmsg=error_message )
      call alloc_error(pname,"iseed",2,__LINE__,__FILE__,ierr,error_message)
    end if
  end if
  !
  ! Only allocate and create the new random seed if the iseed array is
  ! not yet allocated.
  !
  if (.not. allocated(iseed)) then
    !
    call random_seed( size=isize )
    ! isize = 2 (ifort, Athena, 11/2026)
   !if (mypnum == glb_root) write (iout,*) 'isize = ',isize
    iseed = init_array(isize,0) ! F2003 AUTO-REALLOC
    !
    if (present(seed)) then
      !
      do i = 1,isize
        if (any(seed == [0,1])) then
          iseed(i) = seed + i - 1
        else
          iseed(i) = seed**i
        end if
      end do
      !
    else
      !
      if (mypnum == glb_root) then
        !
        call date_and_time( values=date )
        iclk = date(8)
       !write (iout,*) 'iclk  = ',iclk
        !
        do i = 1,isize
          iseed(i) = iclk**i
        end do
        !
      end if
      !
      call mpi_bcast(iseed,size(iseed,kind=int_mpi),mpi_inttyp, &
                     glb_root,MPI_COMM_WORLD,mpierr)
      !
    end if
    !
  end if
  !
  call random_seed( put=iseed )
  !
end subroutine set_random_seed
!
!###############################################################################
!
subroutine get_random_seed()
  !
  !***********************************************************************
  ! Extracts the current seed in the random number generator
  !***********************************************************************
  !
  !.. Local Scalars ..
  integer :: isize
  !
continue
  !
  call random_seed( size=isize )
  !if (mypnum == glb_root) write (iout,*) 'isize = ',isize
  iseed = init_array(isize,0) ! F2003 AUTO-REALLOC
  !
  call random_seed( get=iseed )
  !
end subroutine get_random_seed
!
!###############################################################################
!
end module random_mod
!
!###############################################################################
!
