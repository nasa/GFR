  subroutine cgns_write_field_data()
    !
    ! This routine is placed here simply to prevent from having to place this
    ! whole chunk of code after each variable that needs to be written to the
    ! CGNS file. Everything used here is inherited from the parent routine.
    !
    !.. Local Scalars ..
    integer :: n
    !
    !.. Local Parameters ..
    ! ZERO_CRT : The value 0.0 in CRT precision, i.e. the precision
    !            of real variables written to CGNS files.
    real(CRT), parameter :: zero_CRT = 0.0_CRT
    ! EPS_CRT : CRT machine precision/epsilon, which is the value closest to
    !           zero that can be represented by a real in CRT precision.
    real(CRT), parameter :: eps_CRT = epsilon(zero_CRT)
    ! FTZ_CRT : This parameter contains the same value as eps_CRT but
    !           represented in the regular real working precision (WP). All
    !           working precision values less than or equal to the value of
    !           this parameter will be flushed to zero (FTZ) to prevent any
    !           floating point underflow exceptions from occurring.
    real(wp), parameter :: ftz_CRT = real( eps_CRT , kind=wp )
    !
  continue
    !
    ! Flush-to-zero all elements of the vartmp array which are smaller than
    ! machine precision for the real precision used for writing to CGNS files.
    !
    if (CRT /= wp) then
      where (abs(vartmp(:)) <= ftz_CRT) vartmp = zero
    end if
    !
    ! Copy the vartmp array over to the var array, which
    ! will convert the working precision data to the
    ! real precision used for writing to CGNS files.
    !
    var(:) = real( vartmp , kind(var) )
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
    ! Write the current variable data to the CGNS file
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_field_write_data_f")
    call cgp_field_write_data_f(ifile_s,ibase_s,izone_s,isol_s,ivar_s, &
                                my_pts_beg,my_pts_end,var,cgierr)
    call cgns_error(pname,cgns_sol_file,"cgp_field_write_data_f", &
                    cgierr,__LINE__,__FILE__,"Write data: Solution variable",nv)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_field_write_data_f")
    !
  end subroutine cgns_write_field_data
