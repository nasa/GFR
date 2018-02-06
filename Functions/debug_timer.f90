subroutine debug_timer(time_flag,message,processor,part_num,acc)
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: time_flag
  !
  character(len=*), intent(in) :: message
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: processor
  integer, optional, intent(in) :: part_num
  integer, optional, intent(in) :: acc
  !
#ifdef PROFILE_ON
  !.. Local Scalars ..
  real(wp), save :: beg_wall_time
  real(wp)       :: end_wall_time
  real(wp)       :: wall_time
  !
  integer, save :: procedure_counter = 0
  integer       :: i,hours,minutes,seconds,total_sec
  integer       :: output_cpu
  !
  character(len=300) :: out_message
  !
  !.. Local Arrays ..
  real(wp), save :: beg_procedure_time(10)
  real(wp), save :: accumulated_time(10)
  real(wp), save :: accumulation_beg_time(10)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: leading_space = "                    "
  !
continue
  !
  output_cpu = 0
 !output_cpu = mypnum
  if (present(processor)) output_cpu = processor
  !
  if (time_flag == start_timer) then
    !
    if (present(part_num)) then
      write (out_message,6) message,part_num
    else
      write (out_message,5) message
    end if
    !
    i = 3 + procedure_counter*2
    !
    beg_wall_time = mpi_wtime()
   !beg_wall_time = dclock()
   !call cpu_time( beg_wall_time )
    !
    if (mypnum == output_cpu) then
      write (iout,1) leading_space(1:i),trim(adjustl(out_message))
    end if
    !
  else if (time_flag == stop_timer) then
    !
    end_wall_time = mpi_wtime()
   !end_wall_time = dclock()
   !call cpu_time( end_wall_time )
    wall_time = end_wall_time - beg_wall_time
    !
    total_sec = nint(wall_time)
    !
    hours   = total_sec / 3600
    minutes = (total_sec - hours*3600) / 60
    seconds = total_sec - hours*3600 - minutes*60
    !
    if (present(part_num)) then
      write (out_message,6) message,part_num
    else
      write (out_message,5) message
    end if
    !
    i = 3 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,2) leading_space(1:i),wall_time,trim(adjustl(out_message))
    end if
    !
    ! Reset beg_wall_time in case this procedure is called twice
    ! in a row with time_flag set to stop_timer
    !
    beg_wall_time = mpi_wtime()
   !beg_wall_time = dclock()
   !call cpu_time( beg_wall_time )
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
    if (mypnum == output_cpu) then
      write (iout,2) leading_space(1:i),accumulated_time(acc), &
                     trim(adjustl(out_message))
    end if
    !
  else if (time_flag == start_accumulation) then
    !
    if (present(acc)) then
      accumulation_beg_time(acc) = mpi_wtime()
    end if
    !
  else if (time_flag == stop_accumulation) then
    !
    if (present(acc)) then
      wall_time = mpi_wtime() - accumulation_beg_time(acc)
      accumulated_time(acc) = accumulated_time(acc) + wall_time
    end if
    !
  else if (time_flag == entering_procedure) then
    !
    i = 1 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,3) leading_space(1:i),trim(adjustl(message))
    end if
    !
    procedure_counter = procedure_counter + 1
    !
    beg_procedure_time(procedure_counter) = mpi_wtime()
   !beg_procedure_time(procedure_counter) = dclock()
   !call cpu_time( beg_procedure_time(procedure_counter) )
    !
  else if (time_flag == leaving_procedure) then
    !
    end_wall_time = mpi_wtime()
   !end_wall_time = dclock()
   !call cpu_time( end_wall_time )
    wall_time = end_wall_time - beg_procedure_time(procedure_counter)
    !
    procedure_counter = max( 0 , procedure_counter-1 )
    !
    i = 1 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,4) leading_space(1:i),trim(adjustl(message)),wall_time
    end if
    !
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
