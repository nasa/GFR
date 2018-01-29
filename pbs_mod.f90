module pbs_mod
  !
  use module_kind_types
  use, intrinsic :: iso_c_binding, only : pbs_int => c_long
  !
  implicit none
  !
  private
  !
  public :: pbs_remaining_walltime
  !
#ifdef PBS_ENV
  interface
    subroutine pbs_time_left(seconds_left) bind(c,name="pbs_time_left_")
      import :: pbs_int
      integer(pbs_int) :: seconds_left
    end subroutine pbs_time_left
  end interface
#endif
  !
  integer(pbs_int), save :: initial_wtl
  integer(pbs_int), save :: previous_wtl
  real(wp),         save :: initial_wtime
  real(wp),         save :: previous_wtime
  logical(lk),      save :: pbs_mod_needs_initializing = true
  !
contains
!
!###############################################################################
!
subroutine pbs_remaining_walltime(walltime_expiring,wall_time_remaining)
  !
  !.. Formal Arguments ..
  logical(lk), intent(out) :: walltime_expiring
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(out) :: wall_time_remaining
  !
  !.. Local Scalars ..
  integer(pbs_int) :: current_wtl
  real(wp)         :: current_wtime
  real(wp)         :: time_elapsed
  real(wp)         :: next_time_est
  real(wp)         :: wall_time_left
  !
  !.. Local Arrays ..
  real(wp) :: time_values(1:2)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "pbs_remaining_walltime"
  real(wp), parameter :: max_sec = real( 1000*3600-1 , kind=wp )
  real(wp), parameter :: five_min = real( 300 , kind=wp )
  !
continue
  !
  current_wtl = 0_pbs_int
  current_wtime = zero
  !
  if (mypnum == 0) then
    !
    if (pbs_mod_needs_initializing) then
      call initialize_pbs_mod(walltime_expiring)
    end if
    !
#ifdef PBS_ENV
   !call pbs_time_left(current_wtl)
    ! pbs_time_left seems to not be reliable enough to use it each time
    current_wtl = previous_wtl
    current_wtime = mpi_wtime()
   !write (iout,1) current_wtl,current_wtime
    !
    if (previous_wtl-current_wtl == 0_pbs_int) then
      !
      ! If the time between calls to pbs_time_left was too small such that
      ! the returned integers are the same, use the mpi_wtime() function to
      ! compute the time elapsed between the current time and the previous
      ! time this subroutine was called.
      !
      ! NOTE : previous_wtime should be less than current_wtime since
      !        mpi_wtime counts up as time progresses
      !
      time_elapsed = current_wtime - previous_wtime
      time_elapsed = min( time_elapsed , max_sec )
      wall_time_left = real(initial_wtl,kind=kind(wall_time_left)) - &
                       (current_wtime-initial_wtime)
      !
    else
      !
      ! Use integers returned from pbs_time_left to compute the time
      ! elapsed since the previous time this subroutine was called
      !
      ! NOTE : current_wtl should be less than previous_wtl since pbs_time_left
      !        counts down the remaining time for the job as time progresses
      !
      time_elapsed = real(previous_wtl-current_wtl,kind=kind(time_elapsed))
      time_elapsed = min( time_elapsed , max_sec )
      wall_time_left = real(current_wtl,kind=kind(wall_time_left))
      !
    end if
    !
    next_time_est = time_elapsed * (two+half)
    !
    walltime_expiring = (wall_time_left-next_time_est <= zero) .or. &
                        (wall_time_left-five_min <= zero)
    !
    previous_wtl = current_wtl
    previous_wtime = current_wtime
    !
#else
    wall_time_left = huge(zero)
    walltime_expiring = .false.
#endif
    !
   !write (iout,1) current_wtl,wall_time_left
    time_values(1) = merge(-one,one,walltime_expiring)
    time_values(2) = wall_time_left
    !
  end if
  !
  if (ncpu > 1) then
    call mpi_bcast(time_values,size(time_values,kind=int_mpi), &
                   mpi_flttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
  end if
  !
  walltime_expiring = (time_values(1) < zero)
  if (present(wall_time_remaining)) then
    wall_time_remaining = time_values(2)
   !if (mypnum == 0) write (iout,2) wall_time_remaining
  end if
  !
  1 format (/,2x,"        current_wtl = ",i0,/, &
              2x,"     wall_time_left = ",es25.17)
  2 format (2x,  "wall_time_remaining = ",es25.17,/)
  !
end subroutine pbs_remaining_walltime
!
!###############################################################################
!
subroutine initialize_pbs_mod(walltime_expiring)
  !
  !.. Formal Arguments ..
  logical(lk), intent(out) :: walltime_expiring
  !
  !.. Local Scalars ..
  integer(pbs_int) :: current_wtl
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_pbs_mod"
  !
continue
  !
  walltime_expiring = fals
  pbs_mod_needs_initializing = fals
  !
#ifdef PBS_ENV
  call pbs_time_left(current_wtl)
#else
  current_wtl = huge(int(0,kind=kind(current_wtl)))
#endif
  write (iout,1) current_wtl,ncpu
  !
  ! Initialize the module variables storing the previous wall time left
  ! and previous mpi_wtime
  !
  initial_wtl = current_wtl
  initial_wtime = mpi_wtime()
  !
  previous_wtl = initial_wtl
  previous_wtime = initial_wtime
  !
  ! Format Statements
  !
  1 format (/," Allocated Wall Time = ",i0," s. running on ",i0," CPUs",/)
  !
end subroutine initialize_pbs_mod
!
!###############################################################################
!
end module pbs_mod
