  subroutine cgns_write_field_header(n,vname,ifld,line_number,file_name)
    !
    ! This routine is placed here simply to prevent from having to place this
    ! whole chunk of code after each variable that needs to be written to the
    ! CGNS file. Aside from the formal arguments, everything used here is
    ! inherited from the parent routine.
    !
    !.. Formal Arguments ..
    integer,          intent(in) :: n
    character(len=*), intent(in) :: vname
    integer(CBT),  intent(inout) :: ifld
    integer,          intent(in) :: line_number
    character(len=*), intent(in) :: file_name
    !
  continue
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_field_write_f")
    call cgp_field_write_f(ifile,ibase,izone,isol, &
                           CGNS_KIND_REAL,trim(vname),ifld,cgierr)
    call cgns_error(pname,cgns_sol_file,"cgp_field_write_f", &
                    cgierr,line_number,file_name, &
                    "Write node: Solution variable",n)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_field_write_f")
    !
  end subroutine cgns_write_field_header
