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
subroutine set_random_seed()
  !
  !***********************************************************************
  ! Changes the seed in the random number generator to produce
  ! a different set of numbers
  !***********************************************************************
  !
  !.. Local Scalars ..
  integer :: isize,iclk,i
  !
  !.. Local Arrays ..
  integer :: date(1:8)
  !
continue
  !
  ! Only allocate and create the new random seed if the iseed array is
  ! not yet allocated. If iseed is allocated, that means the random seed
  ! was read in from the restart file. In that case, we just need to
  ! call random_seed with the restart seed.
  !
  if (.not. allocated(iseed)) then
    !
    call random_seed( size=isize )
   !if (mypnum == glb_root) write (iout,*) 'isize = ',isize
    iseed = spread(0,dim=1,ncopies=isize) ! F2003 AUTO-REALLOC
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
  iseed = spread(0,dim=1,ncopies=isize) ! F2003 AUTO-REALLOC
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
