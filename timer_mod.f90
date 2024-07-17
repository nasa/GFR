module timer_mod
  !
  use module_kind_types
  !
  implicit none
  !
  private
  !
  public :: solver_timer
  public :: check_solver_timer
  !
  type :: procedure_times_t
    character(len=:), allocatable :: pname
    integer  :: n = 0
    real(wp) :: times(1:100)
    integer  :: num = 0
    real(wp) :: ave = zero
    logical(lk) :: stall_detected = fals
  contains
    procedure :: update_counter => update_times_counter
    procedure :: update_ave => update_timer_average
  end type procedure_times_t
  !
  type(procedure_times_t), allocatable :: solver(:)
  !
contains
!
!###############################################################################
!
subroutine check_solver_timer()
  !
  !.. Use Statements ..
  use ovar, only : current_rk_stage
  !
  !.. Local Scalars ..
  integer  :: n,idx
  real(wp) :: max_time,ave_time
  logical(lk) :: stall_detected
  !
continue
#ifdef ENABLE_SOLVER_TIMER
  !
  do idx = 1,size(solver)
    !
    if (solver(idx)%stall_detected) then
      !
      n = solver(idx)%n
      !
      max_time = maxval( solver(idx)%times(1:n) )
      ave_time = solver(idx)%ave
      !
      write (iout,1) mypnum,current_rk_stage,n, &
                     solver(idx)%pname,max_time,ave_time
      flush (iout)
      !
    end if
    !
    ! Reset the times and stall indicator
    !
    solver(idx)%n = 0
    solver(idx)%times(:) = zero
    solver(idx)%stall_detected = fals
    !
  end do
  !
  1 format ("CPU # ",i0," RK stage # ",i0," stalled at call # ",i0, &
            " to procedure ",a," [",f12.6," sec. vs ",f12.6," ave.]")
  !
#endif
end subroutine check_solver_timer
!
!###############################################################################
!
subroutine solver_timer(timer_flag,pname,idx)
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: timer_flag
  character(len=*), intent(in) :: pname
  integer,       intent(inout) :: idx
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: time,end_time
  !
continue
#ifdef ENABLE_SOLVER_TIMER
  !
  if (timer_flag == entering_procedure) then
    !
    call find_procedure_index(pname,idx)
    !
    n = solver(idx)%n
    !
    call gfr_timer( solver(idx)%times(n) )
    !
  else if (timer_flag == leaving_procedure) then
    !
    call gfr_timer( end_time )
    !
    call find_procedure_index(pname,idx,dont_update_counter=true)
    !
    n = solver(idx)%n
    !
    time = end_time - solver(idx)%times(n)
    !
    call solver(idx)%update_ave(time)
    !
    solver(idx)%times(n) = time
    !
  end if
  !
#endif
end subroutine solver_timer
!
!###############################################################################
!
subroutine update_timer_average(this,time)
  !
  !.. Passed Argument from Invoking Object ..
  class(procedure_times_t), intent(inout) :: this
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: time
  !
  !.. Local Scalars ..
  integer  :: num
  real(wp) :: rnum,rnp1,ave
  !
continue
  !
  if (this%num == 0) then
    !
    this%ave = time
    this%num = 1
    !
  else
    !
    num = this%num
    ave = this%ave
    !
    rnum = real(num,kind=wp)
    rnp1 = real(num+1,kind=wp)
    !
    if (num < 1000) then
      if (time < two*ten*ave) then
        this%ave = (ave*rnum + time)/rnp1
        this%num = num + 1
      else
        if (num > 100 .and. time > eps1) then
          this%stall_detected = true
        end if 
      end if
    else if (time < two*ten*ave) then
      this%ave = (ave*rnum + time)/rnp1
      this%num = num + 1
    else
      if (time > eps1) then
        this%stall_detected = true
      end if
    end if
    !
  end if
  !
end subroutine update_timer_average
!
!###############################################################################
!
subroutine update_times_counter(this)
  !
  !.. Passed Argument from Invoking Object ..
  class(procedure_times_t), intent(inout) :: this
  !
continue
  !
  this%n = merge( this%n+1 , 1 , this%n+1 <= ubound(this%times,dim=1) )
  !
end subroutine update_times_counter
!
!###############################################################################
!
subroutine find_procedure_index(pname,idx,dont_update_counter)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: pname
  integer,       intent(inout) :: idx
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: dont_update_counter
  !
  !.. Local Scalars ..
  integer :: n,len_pname,ierr
  logical(lk) :: update_counter
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: upname
  !
  !.. Local Allocatable Arrays ..
  type(procedure_times_t), allocatable :: tmp_solver(:)
  !
continue
  !
  update_counter = true
  if (present(dont_update_counter)) then
    update_counter = (.not. dont_update_counter)
  end if
  !
  if (.not. allocated(solver)) then
    !
    allocate ( solver(1:1) , stat=ierr , errmsg=error_message )
    !
    idx = 1
    !
    solver(idx)%pname = pname ! F2003 AUTO-REALLOC
    solver(idx)%times(:) = zero
    solver(idx)%n = 1
    !
    return
    !
  end if
  !
  upname = uppercase(pname) ! F2003 AUTO-REALLOC
  len_pname = len_trim(upname)
  !
  ! First check to see if idx already contains the correct value
  !
  if (idx >= lbound(solver,dim=1) .and. idx <= ubound(solver,dim=1)) then
    if (upname == uppercase(solver(idx)%pname)) then
      ! idx is already good, so dont change it
      if (update_counter) then
        call solver(idx)%update_counter
      end if
      return
    end if
  end if
  !
  ! idx doesnt have the correct value, so loop through solver to find
  ! the element that has the same pname value.
  !
  do n = 1,size(solver)
    if (len_pname == len(solver(n)%pname)) then
      if (upname == uppercase(solver(n)%pname)) then
        idx = n
        if (update_counter) call solver(idx)%update_counter
        return
      end if
    end if
  end do
  !
  ! If we made it here, the input pname doesnt exist yet in solver, so we
  ! need to expand it to add pname
  !
  n = size(solver)
  allocate ( tmp_solver(1:n+1) , stat=ierr , errmsg=error_message )
  !
  tmp_solver(1:n) = solver(1:n)
  !
  call move_alloc( from=tmp_solver , to=solver )
  !
  idx = n + 1
  !
  solver(idx)%pname = pname ! F2003 AUTO-REALLOC
  solver(idx)%times(:) = zero
  solver(idx)%n = 1

end subroutine find_procedure_index
!
!###############################################################################
!
end module timer_mod
