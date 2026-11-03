subroutine debug_timer(time_flag,message,processor,part_num,acc,pid,return_time)
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: time_flag
  !
  character(len=*), intent(in) :: message
  !
  !.. Optional Arguments ..
  integer,  optional, intent(in) :: processor
  integer,  optional, intent(in) :: part_num
  integer,  optional, intent(in) :: acc
  integer,  optional, intent(in) :: pid
  real(wp), optional, intent(out) :: return_time
  !
#ifdef PROFILE_ON
  !.. Local Scalars ..
  real(wp), save :: beg_wall_time = zero
  real(wp)       :: end_wall_time
  real(wp)       :: wall_time
  real(wp)       :: temp_time
  real(wp)       :: report_time
  !
  integer, save :: procedure_counter = 0
  integer, save :: spot_counter = 0
  integer       :: i,output_cpu,output_id
  integer       :: my_spot
  !
  character(len=200) :: msg1,msg2
  !
  character(len=300) :: out_message
  !
  integer, parameter :: max_timers = 50
  integer, parameter :: max_spot_timers = 100
  !
  !.. Local Arrays ..
  real(wp), save :: beg_procedure_time(max_timers)
  real(wp), save :: accumulated_time(max_timers)
  real(wp), save :: accumulation_beg_time(max_timers)
  !
  real(wp),           save :: spot_beg_wall_time(max_spot_timers)
  character(len=200), save :: spot_messages(max_spot_timers)
  !
  !.. Local Parameters ..
  character(len=200) :: leading_space
  !
continue
  !
  if (debug_timer_disabled) return
  !
  write (leading_space,'(200(" "))')
  !
  report_time = zero
  output_cpu = glb_root
 !output_cpu = mypnum
  if (present(processor)) output_cpu = processor
  output_id = mypnum
  if (present(pid)) output_id = pid
  !
  if (time_flag == start_timer) then
    !
    if (present(part_num)) then
      write (out_message,6) message,part_num
    else
      write (out_message,5) message
    end if
    !
    spot_counter = spot_counter + 1
    !
    i = 3 + procedure_counter*2
    !
    call gfr_timer( beg_wall_time )
    !
    spot_messages(spot_counter) = trim(adjustl(message))
    spot_beg_wall_time(spot_counter) = beg_wall_time
    !
    if (output_id == output_cpu) then
      write (iout,1) leading_space(1:i),trim(adjustl(out_message))
    end if
    !
  else if (time_flag == stop_timer) then
    !
    call gfr_timer( end_wall_time )
    !
    my_spot = 0
    msg1 = uppercase(trim(adjustl(message)))
    do i = 1,spot_counter
      msg2 = uppercase(trim(adjustl(spot_messages(i))))
      if (msg1 == msg2) then
        my_spot = i
        exit
      end if
    end do
    !
    if (my_spot > 0) then
      wall_time = end_wall_time - spot_beg_wall_time(my_spot)
      spot_counter = max(0,spot_counter - 1)
      do i = my_spot,spot_counter
        spot_messages(i) = spot_messages(i+1)
        spot_beg_wall_time(i) = spot_beg_wall_time(i+1)
      end do
      spot_messages(spot_counter+1) = ""
      spot_beg_wall_time(spot_counter+1) = zero
    else
      wall_time = end_wall_time - beg_wall_time
    end if
    !
    if (present(part_num)) then
      write (out_message,6) message,part_num
    else
      write (out_message,5) message
    end if
    !
    i = 3 + procedure_counter*2
    !
    if (output_id == output_cpu) then
      write (iout,2) leading_space(1:i),wall_time,trim(adjustl(out_message))
    end if
    !
    report_time = wall_time
    !
    ! Reset beg_wall_time in case this procedure is called twice
    ! in a row with time_flag set to stop_timer
    !
    call gfr_timer( beg_wall_time )
    !
  else if (time_flag == reset_accumulated_time) then
    !
    if (present(acc)) then
      accumulated_time(acc) = zero
    end if
    !
  else if (time_flag == report_accumulated_time) then
    !
    if (present(part_num)) then
      write (out_message,6) message,part_num
    else
      write (out_message,5) message
    end if
    !
    i = 3 + procedure_counter*2
    !
    if (output_id == output_cpu) then
      write (iout,2) leading_space(1:i),accumulated_time(acc), &
                     trim(adjustl(out_message))
    end if
    report_time = accumulated_time(acc)
    !
  else if (time_flag == start_accumulation) then
    !
    if (present(acc)) then
      call gfr_timer( accumulation_beg_time(acc) )
    end if
    !
  else if (time_flag == stop_accumulation) then
    !
    if (present(acc)) then
      !
      call gfr_timer( temp_time )
      !
      wall_time = temp_time - accumulation_beg_time(acc)
      !
      accumulated_time(acc) = accumulated_time(acc) + wall_time
      !
      report_time = wall_time
      !
    end if
    !
  else if (time_flag == entering_procedure) then
    !
    i = 1 + procedure_counter*2
    !
    if (output_id == output_cpu) then
      write (iout,3) leading_space(1:i),trim(adjustl(message))
    end if
    !
    procedure_counter = procedure_counter + 1
    !
    call gfr_timer( beg_procedure_time(procedure_counter) )
    !
  else if (time_flag == leaving_procedure) then
    !
    call gfr_timer( end_wall_time )
    !
    wall_time = end_wall_time - beg_procedure_time(procedure_counter)
    !
    procedure_counter = max( 0 , procedure_counter-1 )
    !
    i = 1 + procedure_counter*2
    !
    if (output_id == output_cpu) then
      write (iout,4) leading_space(1:i),trim(adjustl(message)),wall_time
    end if
    !
    report_time = wall_time
    !
  end if
  !
  if (present(return_time)) then
    return_time = report_time
  end if
  !
  flush (iout)
  !
  ! Format Statements
  !
  1 format (a,"## Starting timer for: ",a)
  2 format (a,"     -- TIME WAS",f12.6," seconds for ",a)
  3 format (a,"=> Entering procedure : ",a)
  4 format (a," <= Leaving procedure : ",a,"   after",f12.6," seconds")
  5 format (a)
  6 format (a," for partition # ",i0)
  !
#endif
  !
end subroutine debug_timer
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
subroutine pure_debug_timer(time_flag,message,processor, &
                            part_num,acc,pid,return_time)
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: time_flag
  !
  character(len=*), intent(in) :: message
  !
  !.. Optional Arguments ..
  integer,  optional, intent(in) :: processor
  integer,  optional, intent(in) :: part_num
  integer,  optional, intent(in) :: acc
  integer,  optional, intent(in) :: pid
  real(wp), optional, intent(out) :: return_time
  !
continue
  !
#ifdef DISABLE_PURE
  call debug_timer(time_flag,message, &
                   processor=processor,part_num=part_num, &
                   acc=acc,pid=pid,return_time=return_time)
#endif
  !
end subroutine pure_debug_timer
!
!###############################################################################
!
