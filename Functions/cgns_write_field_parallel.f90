  subroutine cgns_write_field_parallel()
    !
    ! This routine is placed here simply to prevent from having to place this
    ! whole chunk of code after each variable that needs to be written to the
    ! CGNS file. Everything used here is inherited from the parent routine.
    !
  continue
    !
    if (interpolate_before_output_variable) then
      var = real( vartmp , kind(var) )
    else
      var(:) = real( interpolate_for_output(vartmp) , kind(var) )
    end if
    !
    ! Perform some quick modifications to make sure all values in var
    ! are valid before writing to the CGNS file
    !
    where (is_nan(var)) var = zero
    where (.not. is_finite(var))
      var = merge(-1.0e+10_CRT,1.0e+10_CRT,is_negative(var))
     !where (is_negative(var))
     !  var = -1.0e+10_CRT
     !else where
     !  var =  1.0e+10_CRT
     !end where
    end where
    !
    ! Reset the queue.
    !
    call set_cgns_queue(BEGIN_CGNS_QUEUE,pname,cgns_sol_file)
    !
    ! Add the data for current variable to the output queue.
    !
    call cgp_field_write_data_f(ifile_s,ibase_s,izone_s,isol_s,isvar_s(n), &
                                my_pts_beg,my_pts_end,var,cgierr)
    call cgns_error(pname,cgns_sol_file,"cgp_field_write_data_f", &
                    cgierr,__LINE__,__FILE__,"Write data: Solution variable",n)
    !
    ! Finally, flush the queue containing
    ! the current variable to the CGNS file.
    !
    call flush_cgns_queue(pname,cgns_sol_file)
    !
  end subroutine cgns_write_field_parallel
