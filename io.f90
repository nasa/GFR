module io_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nq
  use eqn_idx, only : ncs
  !
  implicit none
  !
  private
  !
  ! ##################################
  ! ###  Module Public Procedures  ###
  ! ##################################
  !
  public :: read_input_file
  public :: read_grid_file
  public :: write_final_error
  public :: write_residual_file
  public :: write_summary_file
  public :: check_face_flux_point_orderings
  public :: check_cell_flux_point_orderings
  public :: test_communication
  !
  ! ################################
  ! ###  Module Local Variables  ###
  ! ################################
  !
  ! Module variables used for debugging
  !
  integer, save :: current_cell
  !
  interface g2c
    module procedure g2c_i4, g2c_i8, g2c_r4, g2c_r8, g2c_r16
  end interface g2c
  !
contains
!
!###############################################################################
!
subroutine read_input_file(input_mode)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: input_mode
  !
  !.. Local Scalars ..
  integer :: i,j,ierr,input_status
  integer :: command_count,command_length
  character(len=1000) :: command_value
  character(len=1000) :: input_fname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_input_file"
  character(len=*), parameter :: default_fname = "gfr.in"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Check for additional command line arguments
  !
  command_count = command_argument_count()
  !
  one_pass_loop: do i = 1,1
    !
    if (command_count == 1) then
      !
      call get_command_argument(1,input_fname,command_length,ierr)
      !
#define ONLY_ROOT_CHECKS_INPUT_FILE
#ifdef ONLY_ROOT_CHECKS_INPUT_FILE
      call check_input_file_root_only(input_fname,input_status)
#else
      call check_input_file(input_fname,input_status)
#endif
      !
      if (input_status == 0) exit one_pass_loop
      !
    else if (command_count > 1) then
      !
      do j = 1,command_count
        !
        call get_command_argument(j,command_value,command_length,ierr)
        !
        if (trim(adjustl(command_value)) == "-f") then
          !
          call get_command_argument(j+1,input_fname,command_length,ierr)
          !
#ifdef ONLY_ROOT_CHECKS_INPUT_FILE
          call check_input_file_root_only(input_fname,input_status)
#else
          call check_input_file(input_fname,input_status)
#endif
          !
          if (input_status == 0) exit one_pass_loop
          !
        end if
        !
      end do
      !
    end if
    !
    input_fname = default_fname
    !
#ifdef ONLY_ROOT_CHECKS_INPUT_FILE
    call check_input_file_root_only(input_fname,input_status)
#else
    call check_input_file(input_fname,input_status)
#endif
    !
    if (input_status /= 0) then
      call write_nml_input_file('gfr.in.default')
      write (error_message,11)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
    !
  end do one_pass_loop
  call check_for_bad_number(__LINE__,__FILE__)
  !
  call read_nml_input_file(input_fname,input_mode)
  call check_for_bad_number(__LINE__,__FILE__)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format statements
  !
  11 format ("ERROR FINDING A VALID INPUT FILE!")
  !
end subroutine read_input_file
!
!###############################################################################
!
subroutine check_input_file(input_fname,input_status)
  !
  !.. Formal Arguments ..
  character(len=*),  intent(in) :: input_fname
  integer,          intent(out) :: input_status
  !
  !.. Local Scalars ..
  integer :: ierr,inpt
  logical(ldk) :: file_exists
  character(len=80) :: text
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_input_file"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  input_status = -1
  !
  ! Inquire whether the input file name exists
  !
  inquire (file=input_fname,exist=file_exists)
  !
  if (file_exists) then
    !
    ! Try opening the input file name
    !
    open (newunit=inpt,file=input_fname,status="old",action="read", &
                       form="formatted",iostat=ierr,iomsg=error_message)
    call io_error(pname,input_fname,1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Read the first line of the input file
    !
    read (inpt,*,iostat=ierr,err=100,end=100) text
    !
    close (inpt,iostat=ierr,iomsg=error_message)
    call io_error(pname,input_fname,2,__LINE__,__FILE__,ierr,error_message)
    !
    ! Check the first line for neccessary string
    !
    text = trim(adjustl(text))
    text = uppercase( text )
    !
    ! Change input_status to 0 if the first line
    ! matches the required namelist value
    !
    if (text(1:6) == "&INPUT") input_status = 0
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  return
  !
  ! Read error continuation statements
  !
  100 continue
  write (error_message,10) trim(adjustl(input_fname)), ierr
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  !
  ! Format statements
  !
  10 format (" Error checking the format of the input file '",a,"'!", &
             " IO Error = ",i0)
  !
end subroutine check_input_file
!
!###############################################################################
!
subroutine check_input_file_root_only(input_fname,input_status)
  !
  !.. Formal Arguments ..
  character(len=*),  intent(in) :: input_fname
  integer,          intent(out) :: input_status
  !
  !.. Local Scalars ..
  integer :: ierr,inpt
  logical(ldk) :: file_exists
  character(len=80) :: text
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_input_file_root_only"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  input_status = -1
  !
  ! Inquire whether the input file name exists
  !
  if (mypnum == glb_root) then
    !
    inquire (file=input_fname,exist=file_exists)
    !
    if (file_exists) then
      !
      ! Try opening the input file name
      !
      open (newunit=inpt,file=input_fname,status="old",action="read", &
                         form="formatted",iostat=ierr,iomsg=error_message)
      call io_error(pname,input_fname,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Read the first line of the input file
      !
      read (inpt,*,iostat=ierr,err=100,end=100) text
      !
      close (inpt,iostat=ierr,iomsg=error_message)
      call io_error(pname,input_fname,2,__LINE__,__FILE__,ierr,error_message)
      !
      ! Check the first line for neccessary string
      !
      text = trim(adjustl(text))
      text = uppercase( text )
      !
      ! Change input_status to 0 if the first line
      ! matches the required namelist value
      !
      if (text(1:6) == "&INPUT") input_status = 0
      !
    end if
    !
  end if
  !
  if (ncpu > 1) then
    call mpi_bcast(input_status,1_int_mpi,mpi_inttyp, &
                   glb_root,MPI_COMM_WORLD,mpierr)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  return
  !
  ! Read error continuation statements
  !
  100 continue
  write (error_message,10) trim(adjustl(input_fname)), ierr
  call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
  !
  ! Format statements
  !
  10 format (" Error checking the format of the input file '",a,"'!", &
             " IO Error = ",i0)
  !
end subroutine check_input_file_root_only
!
!###############################################################################
!
subroutine read_nml_input_file(input_fname,input_mode)
  !
  !.. Use Statements ..
  use order_mod, only : initialize_geom_solpts
  use order_mod, only : n_order,p_order,q_order,o_order,e_order,t_order
  !
  use geovar, only : nr,using_quadrature,xyz_nodes
  !
  use input_namelist_mod ! use all the variables in the input namelist
  !
  use ovar,   only : Total_Periods,dflux_opt
  use ovar,   only : velref,aref,uref,Lref,time_ref
  use ovar,   only : Period_Timesteps,Period_Time
  use ovar,   only : Grid_Length
  use ovar,   only : LDG_beta,LDG_tau
  use ovar,   only : marker_c1,marker_c2,marker_dpsi
  use ovar,   only : Default_Marker_c1,Default_Marker_c1_JI
  use ovar,   only : Default_Marker_c2,Default_Marker_c2_JI
  use ovar,   only : disable_cgns_output
  use ovar,   only : output_KEDR
  !
  use quadrature_mod, only : init_geom_quadrature_rules
  use quadrature_mod, only : init_face_quadrature_rules
  !
  use restart_mod, only : check_restart_dir
  !
  use oerk_mod,  only : init_optimized_erk_coefficients
  !
  use eperk_mod,  only : init_eperk_coefficients
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: input_fname
  integer,          intent(in) :: input_mode
  !
  !.. Local Scalars ..
  integer  :: i,it1,it2,ierr,ival,inpt
  integer  :: Period_Timesteps_X,Period_Timesteps_Y
  logical(lk)  :: use_qdtr_face_pts,mms_error,pulse_test_error
  real(wp) :: Period_Time_X,Period_Time_Y
  real(wp) :: real_Total_Periods
  real(wp) :: aoa_radians
  real(wp) :: cfl_sum,cfl_time
  real(wp) :: xmin,xmax,psi0,l10p
  !
  !.. Local Arrays ..
  character(len=80), dimension(12) :: input_text
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_nml_input_file"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate all the variables within the input namelist
  ! to their default sizes and values
  !
  call allocate_namelist_variables
  !
  ! Open the input file
  !
  call debug_timer(start_timer,"actually reading the input file")
  open (newunit=inpt,file=input_fname,status="old",action="read", &
                     form="formatted",iostat=ierr,iomsg=error_message)
  call io_error(pname,input_fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read in the input file using the input namelist
  !
  read (inpt,nml=input)
  !
  ! Close the input file now that we are finished reading it
  !
  close (inpt,iostat=ierr,iomsg=error_message)
  call io_error(pname,input_fname,2,__LINE__,__FILE__,ierr,error_message)
  call debug_timer(stop_timer,"actually reading the input file")
  !
  ! Resize all the allocatable variables within the input
  ! namelist to sizes more appropriate to their input values
  !
  call resize_namelist_allocatables
  !
  ! For interpolation pass, shuffle some variables around to trick
  ! the lower level readers to read interpto_gridfile which will be
  ! interpolated onto.  If no name provided for interpto_gridfile
  ! then skip interpolation and treat as normal solver read.
  !
  if (input_mode == 1 .and. len_trim(adjustl(interpto_gridfile)) > 0) then
   !if (ncpu > 1) then
   !  write (error_message,107)
   !  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
   !end if
    grid_format              = interpto_grid_format
    gridfile                 = interpto_gridfile
    plot3d_cutfile           = interpto_plot3d_cutfile
    plot3d_bcsfile           = interpto_plot3d_bcsfile
    grid_scaling_factor      = interpto_grid_scaling_factor
    plot3d_i_stride          = interpto_plot3d_i_stride
    plot3d_j_stride          = interpto_plot3d_j_stride
    plot3d_k_stride          = interpto_plot3d_k_stride
    plot3d_all_strides       = interpto_plot3d_all_strides
    plot3d_agglomerate_order = interpto_plot3d_agglomerate_order
    solution_file            = ""
    !
    solution_order   = interpto_solution_order
    projection_order = interpto_projection_order
    quadrature_order = interpto_quadrature_order
    output_order     = interpto_output_order
    error_order      = interpto_error_order
    all_orders       = interpto_all_orders
    !
  end if
  !
  ! Check the lustre striping parameters
  !
  call check_lustre_striping
  call check_for_bad_number(__LINE__,__FILE__)
  !
  ! Check the output directory
  !
  output_dir = trim(adjustl(output_dir))
  call check_output_dir
  call check_for_bad_number(__LINE__,__FILE__)
  !
  ! Check for the CGNS output directory
  ! NOTE: This will be done later in the post-processor
  !       module if continuous_output is true, so only
  !       do this if continuous_output is false.
  !
  if (.not. continuous_output) then
    call check_cgns_dir
    call check_for_bad_number(__LINE__,__FILE__)
  end if
  !
  ! Check the restart output directory
  !
  call check_restart_dir
  call check_for_bad_number(__LINE__,__FILE__)
  !
  call debug_timer(start_timer,"read_nml_input_file pre read_grid_file")
  !
  ! Adjust the name of the grid and solution files
  !
  gridfile = trim(adjustl(gridfile))
  solution_file = trim(adjustl(solution_file))
  !
  ! Additional checks on solution file.
  !
  if (len_trim(adjustl(solution_file)) > 0) then
    !
    ! Make sure we have a reader for that solution format.
    !
    if (.not. any(solution_format == [Plot3D_Format,Restart_Format])) then
      write (error_message,108)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
    !
    ! Make sure the grid and solution file formats are compatible.
    !
    if ((solution_format == Plot3D_Format .and. &
             grid_format /= Plot3D_Format) .or. &
        (solution_format == CGNS_Format .and. &
             grid_format /= CGNS_Format)) then
      write (error_message,109)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
  end if
  !
  ! Make sure both profile_io_restart and profile_io_cgns
  ! are true if profile_io_all is true
  !
  if (profile_io_all) then
    profile_io_restart = true
    profile_io_cgns = true
  end if
  !
  ! Allow for the use of the variable all_orders to make
  ! solution_order (n_order), projection_order (p_order),
  ! and quadrature_order (q_order) all the same value without
  ! having to explicitly specify each one individually.
  !
  if (all_orders >= 0) then
    solution_order = all_orders
    projection_order = all_orders
    quadrature_order = all_orders
  end if
  !
  ! Set the order variables using the input parameters
  !
  n_order = solution_order
  p_order = merge(projection_order,solution_order,projection_order > 0)
  q_order = merge(quadrature_order,solution_order,quadrature_order > 0)
  o_order = merge(    output_order,solution_order,    output_order > 0)
  e_order = merge(     error_order,solution_order,     error_order > 0)
  t_order = max(n_order,dmat_grad_temp_order)
  !
  ! Allow for the use of the variable all_points to make
  ! loc_solution_pts, loc_flux_pts, and loc_quadrature_pts all the same value
  ! without having to explicitly specify each one individually.
  !
  if ( any(all_points == [1,2]) ) then
    loc_solution_pts = all_points
    loc_flux_pts = all_points
    loc_quadrature_pts = all_points
  end if
  !
  if (loc_solution_pts /= loc_flux_pts) then
    if (mypnum == glb_root) write (iout,101)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  !
  if (loc_solution_pts == Legendre_Gauss_Lobatto) then
    if (continuous_output) then
      write (error_message,106)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
  end if
  !
  ! Disable use_bc_conflicts if Lobatto points are not being used.
  !
  one_pass: do i = 1,1
    !
    if (any(loc_solution_pts   == Gauss_Lobatto)) exit one_pass
    if (any(loc_flux_pts       == Gauss_Lobatto)) exit one_pass
    if (any(loc_quadrature_pts == Gauss_Lobatto)) exit one_pass
    !
    use_bc_conflicts = fals
    !
  end do one_pass
  !
  ! Adjust the solution order and compute the number
  ! of solution points for the different cell shapes
  !
  n_order = max(p_order,n_order)
  !
  ! Make sure the quadrature order is at least the same as the
  ! solution order and other quadrature variables are valid
  !
  q_order = max(n_order,q_order)
  !
  if (q_order > n_order) then
    using_quadrature = true
  end if
  if (loc_solution_pts /= loc_quadrature_pts) then
    using_quadrature = true
  end if
  !
  ! Get the number of solution points for each of the different geometry types
  ! This includes quadrature and face values
  !
  use_qdtr_face_pts = using_quadrature .and. dflux_opt==project_derivatives
  call debug_timer(stop_timer,"read_nml_input_file pre read_grid_file")
  !
  call initialize_geom_solpts(use_qdtr_face_pts)
  call check_for_bad_number(__LINE__,__FILE__)
  !
  ! Now that we have the name of the grid file and the location
  ! of the flux and solution points, read in the grid file
  !
  call read_grid_file
  call check_for_bad_number(__LINE__,__FILE__)
  !
  call init_face_quadrature_rules
  call check_for_bad_number(__LINE__,__FILE__)
  call init_geom_quadrature_rules
  call check_for_bad_number(__LINE__,__FILE__)
  !
  call debug_timer(start_timer,"read_nml_input_file post read_grid_file")
  !
  ! If using a gas mixture and the requested upwind method isnt
  ! compatible, change the upwind method to the Roe solver.
  !
  if (ncs > 0) then
    if (invs_flux_method /= invs_flux_ROE  .and. &
        invs_flux_method /= invs_flux_Rusanov .and. &
        invs_flux_method /= invs_flux_HLLC .and. &
        invs_flux_method /= invs_flux_LDFSS) then
      invs_flux_method = invs_flux_ROE
      if (mypnum == glb_root) then
        write (iout,102)
        flush (iout)
      end if
    end if
  end if
  !
  ! If using a moving reference frame and the requested upwind method
  ! isnt compatible, change the upwind method to the Roe solver.
  !
  if ( moving_ref_frame%is_active() ) then
    if (invs_flux_method /= invs_flux_ROE  .and. &
        invs_flux_method /= invs_flux_Rusanov) then
      invs_flux_method = invs_flux_ROE
      if (mypnum == glb_root) then
        write (iout,110)
        flush (iout)
      end if
    end if
  end if
  !
  ! Make sure the variable index for the resolution marker is valid
  !
  ! Make sure it is non-negative
  Resolution_Marker_Variable = max(0,Resolution_Marker_Variable)
  ! Make sure it is one of the mean flow variables if greater than 0
  if (Resolution_Marker_Variable > 0) then
    if (Resolution_Marker_Variable > nr+2) then
      ! If it is greater than the energy equation (nr+2), reset it to density
      Resolution_Marker_Variable = nec
    else
      ! Otherwise, offset the variable by the number of species equations
      Resolution_Marker_Variable = Resolution_Marker_Variable + ncs
    end if
  end if
 !!
  ! Reset default marker controls when using jump indicator.
  !
  if (Limiter_Option == 2) then
    Default_Marker_c1 = Default_Marker_c1_JI
    Default_Marker_c2 = Default_Marker_c2_JI
    marker_c1 = Default_Marker_c1
    marker_c2 = Default_Marker_c2
  end if
  !
#if DUMP_MARKER_INPUTS == 2
  ! Dump marker inputs before they are processed.
  if (mypnum == glb_root) then
    psi0 = -one*(marker_c1+marker_c2*log10(real(n_order,kind=wp)))
    write (iout,582) "BEFORE", &
                     Limiter_Marker_c1, &
                     Limiter_Marker_c2, &
                     Limiter_Marker_dpsi, &
                     Limiter_Marker_Min, &
                     Limiter_Marker_Max, &
                     marker_c1, &
                     marker_c2, &
                     marker_dpsi, &
                     psi0, &
                     psi0-marker_dpsi, &
                     psi0+marker_dpsi
  end if
  flush (iout)
  582 format (//,2x,a," modifying marker parameters:",//, &
                 4x,"Limiter_Marker_c1   = ",f10.6,/, &
                 4x,"Limiter_Marker_c2   = ",f10.6,/, &
                 4x,"Limiter_Marker_dpsi = ",f10.6,/, &
                 4x,"Limiter_Marker_Min  = ",f10.6,/, &
                 4x,"Limiter_Marker_Max  = ",f10.6,//, &
                 4x,"        marker_c1   = ",f10.6,/, &
                 4x,"        marker_c2   = ",f10.6,/, &
                 4x,"        marker_dpsi = ",f10.6,/, &
                 4x,"        marker_psi0 = ",f10.6,/, &
                 4x,"        marker_min  = ",f10.6,/, &
                 4x,"        marker_max  = ",f10.6,//)
#endif
  !
  if (Limiter_Marker_c1   > zero) marker_c1   = Limiter_Marker_c1
  if (Limiter_Marker_c2   > zero) marker_c2   = Limiter_Marker_c2
  if (Limiter_Marker_dpsi > zero) marker_dpsi = Limiter_Marker_dpsi
  !
  if (Limiter_Marker_Max < zero .and. Limiter_Marker_Min > zero) then
    Limiter_Marker_Min = Limiter_Marker_Max - two*marker_dpsi
  else if (Limiter_Marker_Min < zero .and. Limiter_Marker_Max > zero) then
    Limiter_Marker_Max = Limiter_Marker_Min + two*marker_dpsi
  end if
  !
  l10p = log10(real(n_order,kind=wp))
  !
  if (all([Limiter_Marker_Min,Limiter_Marker_Max] < zero)) then
    !
    marker_dpsi = half*abs(Limiter_Marker_Min-Limiter_Marker_Max)
    !
    ! Confirmed via Mathematica that the following options:
    !   1) holding marker_c1 constant and using it plus Limiter_Marker_min
    !      and marker_dpsi to compute a new marker_c2
    !   2) holding marker_c2 constant and using it plus Limiter_Marker_min
    !      and marker_dpsi to compute a new marker_c1
    ! result in the SAME psi0, which makes them equivalent.
    !
    if (marker_c1 /= Default_Marker_c1) then
      marker_c2 = -one*(Limiter_Marker_Min + marker_dpsi + marker_c1) / l10p
    else
      marker_c1 = -one*(Limiter_Marker_Min + marker_dpsi + marker_c2*l10p)
    end if
    !
  end if
  !
  ! Copy the computed values back into the input parameters,
  ! so that the echoed user inputs reflect what the code really uses.
  !
  Limiter_Marker_c1 = marker_c1
  Limiter_Marker_c2 = marker_c2
  Limiter_Marker_dpsi = marker_dpsi
  !
  Limiter_Marker_Min = -one*(marker_c1 + marker_c2*l10p) - marker_dpsi
  Limiter_Marker_Max = -one*(marker_c1 + marker_c2*l10p) + marker_dpsi
 !!
#if DUMP_MARKER_INPUTS == 2
  ! Dump marker inputs after they are processed. Then stop.
  if (mypnum == glb_root) then
    psi0 = -one*(marker_c1+marker_c2*log10(real(n_order,kind=wp)))
    write (iout,582) " AFTER", &
                     Limiter_Marker_c1, &
                     Limiter_Marker_c2, &
                     Limiter_Marker_dpsi, &
                     Limiter_Marker_Min, &
                     Limiter_Marker_Max, &
                     marker_c1, &
                     marker_c2, &
                     marker_dpsi, &
                     psi0, &
                     psi0-marker_dpsi, &
                     psi0+marker_dpsi
  end if
#endif
  !
#if DUMP_MARKER_INPUTS == 1
  if (mypnum == glb_root) then
    psi0 = -one*(marker_c1 + marker_c2*log10(real(n_order,kind=wp)))
    write (iout,30) marker_c1,marker_c2,marker_dpsi, &
                    psi0-marker_dpsi,psi0,psi0+marker_dpsi
  end if
  30 format (/,2x,"marker_c1   = ",f11.6,/, &
               2x,"marker_c2   = ",f11.6,/, &
               2x,"marker_dpsi = ",f11.6,//, &
               2x,"psi0 - dpsi = ",f11.6,/, &
               2x,"psi0        = ",f11.6,/, &
               2x,"psi0 + dpsi = ",f11.6,/)
#endif
  !
#if DUMP_MARKER_INPUTS && DUMP_MARKER_INPUTS_THEN_STOP
  flush (iout)
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,"Dumping marker parameters!")
#endif
  !
  ! Make sure the number of RK stages is correct for some RK types
  !
  if (.not. any(abs(Runge_Kutta_Scheme) == Valid_RK_Schemes(:))) then
    write (error_message,31)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  else if (Runge_Kutta_Scheme == TVD2_RK) then
    num_rk_stages = 2
  else if (Runge_Kutta_Scheme == TVD3_RK) then
    num_rk_stages = 3
  else if (Runge_Kutta_Scheme == TVD4_RK) then
    num_rk_stages = 5
  else if (Runge_Kutta_Scheme == CK4_RK) then
    num_rk_stages = 5
  else if (any(abs(Runge_Kutta_Scheme) == OERK_Schemes)) then
    call init_optimized_erk_coefficients
    Runge_Kutta_Scheme = OERK
  else if (any(abs(Runge_Kutta_Scheme) == EPERK_Schemes)) then
    call init_eperk_coefficients
    if (any(abs(Runge_Kutta_Scheme) == OEPERK_Schemes)) then
      Runge_Kutta_Scheme = EPERK
    else
      Runge_Kutta_Scheme = -EPERK
    end if
  end if
  !
  ! Make sure the initialization and error locations are valid
  !
  if (all(loc_init_pts /= Legendre_Nodes)) loc_init_pts = loc_solution_pts
  if (all(loc_error_pts /= Legendre_Nodes)) loc_error_pts = loc_solution_pts
  !
  ! Make sure the output order is greater than or
  ! equal to the solution order if it is positive
  !
  if (o_order > 0) o_order = max(o_order,n_order+0)
  !
  ! Make sure the error order is greater than or
  ! equal to the solution order if it is positive
  !
  if (e_order > 0) e_order = max(e_order,n_order+0)
  !
  ! Check that the Method of Manufactured Solutions (MMS) option is valid
  !
  if (mms_opt /= 0) then
    !
    mms_error = fals
    if (any(mms_opt == [1,2]) .and. nr < 2) then ! 1 and 2 are for 2D or 3D
      write (error_message,21)
      mms_error = true
    else if (mms_opt == 3 .and. nr > 2) then ! 3 is for 1D or 2D
      write (error_message,22)
      mms_error = true
    else if (mms_opt == 4 .and. nr /= 2) then ! 4 is for 2D only
      write (error_message,23)
      mms_error = true
    else if (all(mms_opt /= [1,2,3,4])) then
      write (error_message,24)
      mms_error = true
    end if
    if (mms_error) then
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
    !
    ! If mms_output_solution is true, set the MMS initialization factor
    ! to one and set the number of time steps and final time to zero
    !
    if (mms_output_solution) then
      mms_init = one
      num_timesteps = 0
      Final_Time = zero
    end if
    !
  end if
  !
 !if (itestcase == Channel_Flow) then
 !  call read_channel_init_file
 !  call check_for_bad_number(__LINE__,__FILE__)
 !end if
  !
  ! Check that input for density pulse test case is valid.
  !
  if (itestcase == Density_Pulse_Test) then
    pulse_test_error = true
    if      ( .not. any( pulse_test_form == [0,1,2,3,4] ) ) then
      write (error_message,41)
    else if ( .not. any( pulse_test_form_dim == [1,2] ) ) then
      write (error_message,42)
    else if ( pulse_test_length .le. 0._wp ) then
      write (error_message,43)
    else if ( pulse_test_width  .le. 0._wp ) then
      write (error_message,44)
    else if ( pulse_test_gap    .lt. 0._wp ) then
      write (error_message,45)
    else if ( pulse_test_repeat .lt. 0 ) then
      write (error_message,46)
    else
      pulse_test_error = fals
    endif
    if (pulse_test_error) then
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
  end if
  !
  ! Compute periodic information for the transport problems
  !
  if (any(itestcase == Transport_Problems)) then
    !
    cfl = one
    cfl_beg = one
    cfl_end = one
    cfl_cycles = 1
    !
    if (abs(velref(1)) < eps3) then
      !
      Period_Time_Y = Grid_Length(2) * reciprocal(velref(2)/aref)
      Period_Time_Y = merge( Period_Time_Y, &
                             Grid_Length(2), &
                             Period_Time_Y > eps3 )
      Period_Timesteps_Y = nint(Period_Time_Y/constant_dt)
      Period_Timesteps = Period_Timesteps_Y
      !
     !if (mypnum == glb_root) then
     !  write (*,*) Period_Time_Y,Period_Timesteps_Y,Period_Timesteps
     !end if
      !
    else if (abs(velref(2)) < eps3) then
      !
      Period_Time_X = Grid_Length(1) * reciprocal(velref(1)/aref)
      Period_Time_X = merge( Period_Time_X, &
                             Grid_Length(1), &
                             Period_Time_X > eps3 )
      Period_Timesteps_X = nint(Period_Time_X/constant_dt)
      Period_Timesteps = Period_Timesteps_X
      !
     !if (mypnum == glb_root) then
     !  write (*,*) Period_Time_Y,Period_Timesteps_Y,Period_Timesteps
     !end if
      !
    else
      !
      Period_Time_X = Grid_Length(1) / (velref(1)/aref)
      Period_Time_Y = Grid_Length(2) / (velref(2)/aref)
      !
      Period_Timesteps_X = nint(Period_Time_X/constant_dt)
      Period_Timesteps_Y = nint(Period_Time_Y/constant_dt)
      !
      if (Period_Timesteps_X == Period_Timesteps_Y) then
        Period_Timesteps = Period_Timesteps_X
      else
        Period_Timesteps = Period_Timesteps_X * Period_Timesteps_Y &
                         / mod( max(Period_Timesteps_X,Period_Timesteps_Y), &
                                min(Period_Timesteps_X,Period_Timesteps_Y) )
      end if
      !
     !if (mypnum == glb_root) then
     !write (*,*) Period_Time_X,Period_Time_Y
     !write (*,*) Period_Timesteps_X,Period_Timesteps_Y
     !write (*,*) Period_Timesteps
     !end if
      !
    end if
    !
    Period_Time = constant_dt*real(Period_Timesteps,kind=wp)
    !
    ! If Final_Time is negative, it is giving the total
    ! number of periods that the simulation will run
    !
    if (Final_Time < zero) then
      Total_Periods = nint( abs(Final_Time) )
      Final_Time = Period_Time * real(Total_Periods,kind=wp)
    else
      real_Total_Periods = Final_Time*reciprocal(Period_Time)
      Total_Periods = nint(real_Total_Periods)
      if (abs(real_Total_Periods-real(Total_Periods,kind=wp)) > ten*eps6) then
        if (mypnum == glb_root) write (iout,103) real_Total_Periods
      end if
    end if
   !if (mypnum == glb_root) write (*,*) Period_Time,Total_Periods,Final_Time
    !
    num_timesteps = nint( Final_Time / constant_dt )
    !
  else if (any(itestcase == Taylor_Problems)) then
    !
    output_KEDR = true
    !
   !cfl_beg = CFL
   !cfl_end = CFL
   !cfl_cycles = 1
    !
    if (Final_Time > zero) then
      !
      ! Need to adjust the nondimensional final time which uses uref as the
      ! velocity for nondimensionalization and we use aref within the code
      !
     !num_timesteps = nint( Final_Time / &
     !                      (constant_dt*time_ref*time_scaling_factor) )
      !
     !Final_Time = Final_Time * time_ref
      !
    end if
    !
    ! Make sure the timestep type is global (positive)
    !
    Timestep_Type = abs(Timestep_Type)
    !
  else
    !
    ! time_ref is now set in initialize_reference
    !
    if (convert_restart_to_cgns) then
      !
      num_timesteps = 0
      !
    else if (Timestep_Type == Constant_Timestep) then
      !
      if (cfl_cycles < 2) then
        cfl_beg = CFL
        cfl_end = CFL
        cfl_cycles = 1
        if (Final_Time > zero) then
          num_timesteps = nint( Final_Time / (CFL*constant_dt*time_ref) )
        end if
      else
        cfl_end = max(cfl_end,eps6)
        cfl_beg = min(cfl_beg,cfl_end)
        cfl_sum = real(cfl_cycles,kind=wp)*(cfl_beg+cfl_end)/two
        cfl_time = cfl_sum*constant_dt*time_ref
        if (Final_Time > zero) then
          num_timesteps = cfl_cycles
          if (cfl_time < Final_Time) then
            num_timesteps = num_timesteps + &
                     nint( (Final_Time-cfl_time) / (CFL*constant_dt*time_ref) )
          end if
        end if
      end if
      !
      !
    else
      !
      if (cfl_cycles < 2) then
        cfl_beg = CFL
        cfl_end = CFL
        cfl_cycles = 0
      else
        cfl_end = max(cfl_end,eps6)    ! cfl_end must be >= 1.0e-6
        cfl_beg = min(cfl_beg,cfl_end) ! cfl_beg must be <= cfl_end
      end if
      !
    end if
    !
  end if
  !
  ! Set the time scaling factor to the flow-through-time for the channel flow
  !
  if (abs(itestcase) == Channel_Flow) then
    !
    xmin = minval(xyz_nodes(1,:))
    xmax = maxval(xyz_nodes(1,:))
    !
    time_scaling_factor = (machref*aref) / abs(xmax-xmin)
    !
  end if
  !
  ! Determine the frequency of dumping the results and solution files
  !
  if (any(itestcase == Transport_Problems)) then
    !
    if (results_interval == 0) then
      results_interval = max(Period_Timesteps/500,1)
    else if (results_interval < 0) then
      results_interval = max(Period_Timesteps/abs(results_interval),1)
    end if
    !
    if (output_interval == 0) then
      output_interval = Period_Timesteps
    else if (output_interval < 0) then
      output_interval = abs(output_interval)*Period_Timesteps
    end if
    !
    if (iter_out_interval == 0) then
      iter_out_interval = Period_Timesteps
    else if (iter_out_interval < 0) then
      iter_out_interval = max(Period_Timesteps/abs(iter_out_interval),1)
    end if
    !
  else
    !
    if (results_interval == 0) then
      results_interval = num_timesteps/500
    else if (results_interval < 0) then
      results_interval = num_timesteps/abs(results_interval)
    end if
    results_interval = max(results_interval,1)
    !
    if (any(itestcase == Taylor_Problems)) then
      if (output_interval == 0) then
        disable_cgns_output = true
      end if
    end if
    !
    if (output_frequency > zero) then
      output_interval = num_timesteps+1
    else
      if (output_interval == 0) then
        output_interval = num_timesteps/50
      else if (output_interval < 0) then
        output_interval = num_timesteps/abs(output_interval)
      end if
      output_interval = max(output_interval,1)
    end if
    !
    if (iter_out_interval == 0) then
      iter_out_interval = 1
    else if (iter_out_interval < 0) then
      iter_out_interval = num_timesteps/abs(iter_out_interval)
    end if
    iter_out_interval = max(iter_out_interval,1)
    !
  end if
  !
  if (visc_flux_method /= visc_flux_BR2) then
    if (visc_flux_method == visc_flux_LDG) then
      LDG_beta = half
    end if
    if (nr == 2) then
      LDG_tau = one/ten
     !LDG_tau = zero
    else if (nr == 3) then
      LDG_tau = one
    end if
  end if
  !
  ! Now write out a copy of the input file to standard output
  !
  if (mypnum == glb_root) then
    !
    ! If interpolating, make sure to only write once with the
    ! non-manipulated namelist input values.
    !
    if (input_mode == 2 .or. len_trim(adjustl(interpto_gridfile)) == 0) then
      !
      write (iout,104)
      write (iout,nml=input,delim='apostrophe')
      write (iout,105)
      !
    end if
    !
  end if
  !
  !$ACC UPDATE DEVICE(nr,pr,walls_are_exact,sub_inflow_method, &
  !$ACC&              LDG_beta,LDG_tau,visc_variables_to_ave, &
  !$ACC&              invs_flux_method,Resolution_Marker_Variable)
  !
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  call debug_timer(stop_timer,"read_nml_input_file post read_grid_file")
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  101 format (/, &
     "***************************** ERROR *****************************",/, &
     "   THE REQUESTED LOCATIONS FOR SOLUTION POINTS AND FLUX POINTS",/, &
     "   DO NOT MATCH. CURRENTLY, THE PROCEDURES USED TO CREATE THE",/, &
     "   CORRECTION MATRICES ASSUMES THESE TWO POINTS LOCATIONS ARE",/, &
     "   THE SAME. THESE PROCEDURES WILL NEED TO BE UPDATED IN ORDER",/, &
     "   TO ALLOW FOR THIS CAPABILITY.",/, &
     "***************************** ERROR *****************************",/)
  102 format (/, &
      " ************************** WARNING!!!! **************************",//, &
      "  AN INVALID INVISCID FLUX SCHEME FOR A GAS MIXTURE WAS REQUESTED",/, &
      "  BY THE INPUT FILE. THE VALID UPWINDING INVISCID FLUX SCHEMES",/, &
      "  FOR A GAS MIXTURE ARE THE ROE, HLLC, OR LDFSS SCHEMES.",//, &
      "    THE INVISCID FLUX SCHEME HAS BEEN RESET TO THE ROE SCHEME.",//, &
      " ************************** WARNING!!!! **************************",/)
  103 format (/, &
     " ********************* WARNING!!! *********************",/, &
     "  IT APPEARS THE FINAL TIME DOES NOT CORRESPOND TO THE",/, &
     "  EXACT END OF A PERIOD FOR THE VORTEX PROBLEM. THE",/, &
     "  FINAL TIME CORRESPONDS TO",f12.5," TOTAL PERIODS.",/, &
     " ********************* WARNING!!! *********************",/)
  104 format (/,"#######################################",/, &
                "#####      INPUT PARAMETERS       #####",/, &
                "#######################################",/)
  105 format (/,"#######################################",/, &
                "#####   END OF INPUT PARAMETERS   #####",/, &
                "#######################################",/)
  106 format ("Lobatto solutions points are not yet compatible with the ", &
              "option to correct the output solution so that it is continuous!")
  107 format ("Interpolation cannot be done in parallel.")
  108 format (/, &
     "***************************** ERROR **********************",/, &
     "   NO READER AVAILABLE FOR SPECIFIED SOLUTION FORMAT.",/, &
     "***************************** ERROR **********************",/)
  109 format (/, &
     "***************************** ERROR **********************",/, &
     "   THE SPECIFIED GRID AND SOLUTION FILE FORMATS ARE NOT",/, &
     "   COMPATIBLE.",/, &
     "***************************** ERROR **********************",/)
  110 format (/, &
      " ************************** WARNING!!!! **************************",//, &
      "  AN INVALID INVISCID FLUX SCHEME FOR A MOVING REFERENCE FRAME",/, &
      "  WAS REQUESTED BY THE INPUT FILE. THE VALID UPWINDING INVISCID",/, &
      "  FLUX SCHEMES FOR A MOVING REFERENCE FRAME ARE THE ROE OR",/, &
      "  RUSANOV SCHEMES. THE INVISCID FLUX SCHEME HAS BEEN RESET TO",/, &
      "  THE ROE SCHEME.",//, &
      " ************************** WARNING!!!! **************************",/)
  21 format ("MMS ERROR: MMS solutions 1 and 2 are allowed for 2D and 3D only!")
  22 format ("MMS ERROR: MMS solution 3 is allowed for 1D and 2D only!")
  23 format ("MMS ERROR: MMS solution 4 is allowed for 2D only!")
  24 format ("MMS ERROR: The MMS solution chosen is invalid!")
  !
  31 format ("RUNGE-KUTTA ERROR: Requested Runge-Kutta scheme is invalid!")
  !
  41 format('ERROR: Density pulse test requres  pulse_test_form between 0:4.')
  42 format('ERROR: Density pulse test requres  pulse_test_form_dim between 1:2.')
  43 format('ERROR: Density pulse test requres  pulse_test_length >0.')
  44 format('ERROR: Density pulse test requres  pulse_test_width >0.')
  45 format('ERROR: Density pulse test requres  pulse_test_gap >=0.')
  46 format('ERROR: Density pulse test requres  pulse_test_repeat >=0.')
  !
end subroutine read_nml_input_file
!
!###############################################################################
!
subroutine write_nml_input_file(output_fname)
  !
  !.. Use Statements ..
  use input_namelist_mod ! use all the variables in the input namelist
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: output_fname
  !
  !.. Local Scalars ..
  integer ionml,ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_nml_input_file"
  !
continue
  !
  ! Allocate all the variables within the input namelist
  ! to their default sizes and values
  !
  call allocate_namelist_variables(opt_n_bc_input=4, &
                                   opt_n_species=3, &
                                   opt_char_len=0)
  !
  if (mypnum == glb_root) then
    open (newunit=ionml,file=output_fname,status="new",action="write", &
                        delim='apostrophe', &
                        position="rewind",iostat=ierr,iomsg=error_message)
    call io_error(pname,trim(adjustl(gridfile)),1,__LINE__,__FILE__, &
                  ierr,error_message)
    write (ionml,nml=input,delim='apostrophe')
    write (iout,101) trim(output_fname)
    close(ionml)
  end if
  !
  101 format(/,"Default namelist input file written to:  ",A)
  !
end subroutine write_nml_input_file
!
!###############################################################################
!
subroutine check_lustre_striping
  !
  !.. Use Statements ..
  use ovar, only : lustre_stripe_count
  use ovar, only : lustre_stripe_size
  use ovar, only : lustre_stripe_unit
  use ovar, only : lustre_stripe_factor
  !
  !.. Local Scalars ..
  integer :: i,j,ip,im,ierr,stripe_num
  character(len=1) :: stripe_unit
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_lustre_striping"
  character(len=*), parameter :: units = "KMG"
  character(len=*), parameter :: numbers(1:10) = ["0","1","2","3","4", &
                                                  "5","6","7","8","9"]
  !
  integer, parameter :: max_lustre_count = 120
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Check the Lustre stripe count
  !
  if (mod(ncpu,lustre_stripe_count) /= 0) then
    !
    count_loop: do i = 1,lustre_stripe_count
      ip = min( lustre_stripe_count + i , max_lustre_count )
      im = max( lustre_stripe_count - i , 0 )
      if (mod(ncpu,ip) == 0) then
        if (ncpu/ip >= 4) then
          lustre_stripe_count = ip
          exit count_loop
        end if
      else if (mod(ncpu,im) == 0) then
        lustre_stripe_count = im
        exit count_loop
      end if
    end do count_loop
    !
  end if
  !
  ! Save the stripe count for MPI
  !
  write (lustre_stripe_factor,4) lustre_stripe_count
  !
  ! Check the Lustre stripe size
  !
  stripe_num = 0
  !
  lustre_stripe_size = trim(adjustl(lustre_stripe_size))
  !
  size_loop: do i = 1,len_trim(lustre_stripe_size)
    !
    if (all(lustre_stripe_size(i:i) /= numbers)) then
      !
      read (lustre_stripe_size(1:i-1),*,iostat=ierr) stripe_num
      !
      if (ierr /= 0) then
        stripe_num = 0
        exit size_loop
      end if
      !
      j = i-1 + scan( uppercase(lustre_stripe_size(i:)) , units )
      !
      if (j /= 0) then
        !
        select case (uppercase(lustre_stripe_size(j:j)))
          case ("K")
            stripe_num = stripe_num * 1024
          case ("M")
            stripe_num = stripe_num * 1024 * 1024
          case ("G")
            stripe_num = stripe_num * 1024 * 1024 * 1024
        end select
        !
      end if
      !
      exit size_loop
      !
    end if
    !
  end do size_loop
  !
  stripe_num = max(0,stripe_num)
  ! Make sure the stripe size is a multiple of 64 KB (65536 Bytes)
  stripe_num = nint( real(stripe_num,kind=wp) / 65536.0_wp ) * 65536
  ! Save the stripe number in bytes for MPI
  write (lustre_stripe_unit,4) stripe_num
  ! Convert the stripe size to kilobytes
  stripe_num  = stripe_num / 1024
  stripe_unit = "k"
  !
  if (stripe_num == 0) then
    !
    if (mypnum == glb_root) write (iout,1)
    lustre_stripe_size = "0"
    !
  else
    !
    if (stripe_num >= 1024) then
      stripe_num = stripe_num / 1024
      stripe_unit = "m"
    end if
    !
    if (stripe_num >= 1024) then
      stripe_num = stripe_num / 1024
      stripe_unit = "g"
    end if
    !
    write (lustre_stripe_size,2) stripe_num,stripe_unit
    !
  end if
  !
#ifdef PBS_ENV
  if (mypnum == glb_root) then
    write (iout,3) lustre_stripe_count, &
                   trim(adjustl(lustre_stripe_size))
  end if
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/,20("!"),5x,"WARNING",5x,20("!"),/, &
              "!!!! THE LUSTRE STRIPE SIZE WAS NOT VALID AND  !!!!",/, &
              "!!!! SO IT HAS BEEN RESET TO THE DEFAULT VALUE !!!!",/, &
              20("!"),5x,"WARNING",5x,20("!"),/)
  2 format (i0,a1)
  3 format (/,"The Lustre striping command to be used is:",/, &
              "lfs setstripe -c ",i0," -S ",a,1x,"[file/directory]",/)
  4 format (i0)
  !
end subroutine check_lustre_striping
!
!###############################################################################
!
subroutine check_output_dir
  !
  !.. Use Statements ..
  use ovar, only : dump_fluxes
  use ovar, only : output_dir
  !
  !.. Local Scalars ..
  integer :: clen,ierr,iodir,ii
  integer :: exitstat,cmndstat
  character(len=200) :: test_file
  character(len=200) :: test_dir
  character(len=300) :: command
  character(len=300) :: cmndmsg
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_output_dir"
  !
  character(len=*), parameter :: nobackup = "/nobackup/sspiegel"
  character(len=*), parameter :: homedir = "/u/sspiegel"
  !
  logical(ldk), parameter :: wait_for_completion = .true.
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize test_file and test_dir
  !
  test_file = " "
  test_dir = trim(adjustl(output_dir))
  !
  ! Only perform the following tests using the root processor
  ! to minimize contention and prevent each MPI process from
  ! executing the same command.
  !
  root_proc_only: if (mypnum == glb_root) then
    !
    one_pass: do ii = 1,1
      !
      ! Check for a directory "/" at the end of the input output_dir
      !
      test_dir = trim(adjustl(output_dir))
      !
      clen = len_trim(test_dir)
      if (test_dir(clen:clen) /= "/") then
        test_dir = trim(test_dir)//"/"
      end if
      !
      ! Try to open a file in the given directory
      !
      test_file = trim(test_dir)//"output_dir.test_file.dat"
      !
      open (newunit=iodir,file=trim(test_file),status="replace", &
                          action="write",iostat=ierr)
      !
      ! If the open was successful, close and delete the test file and return
      !
      if (ierr == 0) then
        close (iodir,status="delete",iostat=ierr)
        output_dir = trim(adjustl(test_dir))
        exit one_pass
      end if
      !
      ! Reset test_file and test_dir
      !
      test_file = " "
      test_dir = trim(adjustl(output_dir))
      !
      ! Try to append nobackup to the beginning of output_dir
      ! and try a test file in that directory
      !
      if (test_dir(1:len(nobackup)) /= nobackup) then
        !
        if (test_dir(1:1) /= "/") then
          test_dir = nobackup//"/"//trim(test_dir)
        else
          test_dir = nobackup//trim(test_dir)
        end if
        !
        ! Try to open a file in the new output_dir
        !
        test_file = trim(test_dir)//"output_dir.test_file.dat"
        !
        open (newunit=iodir,file=trim(test_file),status="replace", &
                            action="write",iostat=ierr)
        !
        ! If the open was successfull, close and delete the test file,
        ! set output_dir to test_dir, and exit the one pass loop
        !
        if (ierr == 0) then
          close (iodir,status="delete",iostat=ierr)
          output_dir = trim(adjustl(test_dir))
          exit one_pass
        end if
        !
      end if
      !
      ! Reset test_file and test_dir
      !
      test_file = " "
      test_dir = trim(adjustl(output_dir))
      !
      ! Finally, try to append homedir to the beginning of output_dir
      ! and try a test file in that directory
      !
      if (test_dir(1:len(homedir)) /= homedir) then
        !
        if (test_dir(1:1) /= "/") then
          test_dir = homedir//"/"//trim(test_dir)
        else
          test_dir = homedir//trim(test_dir)
        end if
        !
        ! Try to open a file in the new output_dir
        !
        test_file = trim(test_dir)//"output_dir.test_file.dat"
        !
        open (newunit=iodir,file=trim(test_file),status="replace", &
                            action="write",iostat=ierr)
        !
        ! If the open was successfull, close and delete the test file,
        ! set output_dir to test_dir, and exit the one pass loop
        !
        if (ierr == 0) then
          close (iodir,status="delete",iostat=ierr)
          output_dir = trim(adjustl(test_dir))
          exit one_pass
        end if
        !
      end if
      !
      ! If we made it here, we cant seem to find the output directory that
      ! was given in the input file so preprend "BAD_DIR" to output_dir
      !
      output_dir = "BAD_DIR" // trim(output_dir)
      !
    end do one_pass
    !
  end if root_proc_only
  !
  ! Broadcast the value of output_dir to all the non-root MPI processes
  ! NOTE: Use test_dir instead of output_dir for the MPI broadcasts since
  !       some versions of MPI have problems with allocatable characters
  !
  test_dir = trim(adjustl(output_dir))
  clen = len_trim(test_dir)
  !
  call mpi_bcast(clen,1_int_mpi,mpi_inttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(test_dir(1:clen),int(clen,kind=int_mpi),MPI_CHARACTER, &
                 0_int_mpi,MPI_COMM_WORLD,mpierr)
  !
  ! If output_dir starts with "BAD_DIR", we cant seem to find the
  ! output directory that was given in the input file so exit
  !
  if (test_dir(1:7) == "BAD_DIR") then
    !
    output_dir = test_dir(8:clen)
    if (mypnum == glb_root) write (iout,1) trim(adjustl(output_dir))
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
    !
  else
    !
    output_dir = trim(adjustl(test_dir))
    !
    ! Create the flux directory if the dump_flux flag was given
    ! in the input file
    !
    if (dump_fluxes) then
      if (mypnum == glb_root) then
        !
        command = "mkdir -p "//trim(output_dir)//"flux_files"
        !
        call execute_command(command,wait_for_completion,exitstat, &
                             cmndstat,cmndmsg,pname)
       !call execute_command_line( trim(command) , wait_for_completion , &
       !                           exitstat , cmndstat , cmndmsg )
       !call write_command_results(pname,command,exitstat,cmndstat,cmndmsg)
        !
        if (exitstat /= 0 .or. cmndstat /= 0) then
          write (error_message,2) trim(output_dir)//"flux_files"
          call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
        end if
        !
      end if
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  1 format (/,"******************** ERROR! *********************",/, &
              "     CANT SEEM TO FIND THE OUTPUT DIRECTORY",/, &
              "       THAT WAS GIVEN IN THE INPUT FILE.",/, &
              "      '",a,"'",/, &
              "     PLEASE CHECK THE INPUT FILE TO CORRECT",/, &
              "        ANY MISTAKES BEFORE CONTINUING!",/, &
              "******************** ERROR! *********************",/)
  2 format ("ERROR CREATING THE DIRECTORY: ",a)
  !
end subroutine check_output_dir
!
!###############################################################################
!
subroutine check_cgns_dir
  !
  !.. Use Statements ..
  use ovar, only : output_dir,cgns_dir
  use ovar, only : lustre_stripe_count
  use ovar, only : lustre_stripe_size
  use ovar, only : matrix_output_unit
  !
  !.. Local Scalars ..
  integer :: clen,ierr,iodir
  integer :: exitstat,cmndstat
  logical(lk) :: using_cgns_dir
  character(len=300) :: command
  character(len=300) :: cmndmsg
  character(len=200) :: test_file
  character(len=200) :: test_dir
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_cgns_dir"
  !
  logical(ldk), parameter :: wait_for_completion = .true.
  !
continue
  !
  if (matrix_output_unit /= iout) return
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize exitstat and cmndmsg
  !
  exitstat = 0
  cmndmsg = ""
  !
  ! Try to open a file in the directory named CGNS_files
  !
  test_dir = trim(output_dir) // "CGNS_files/"
  test_file = trim(test_dir) // "cgns_dir.test_file.dat"
  !
  ! Only perform the following tests using the root processor
  ! to minimize contention and prevent each MPI process from
  ! executing the same command.
  !
  if (mypnum == glb_root) then
    !
    open (newunit=iodir,file=trim(test_file),status="replace", &
                        action="write",iostat=ierr)
    !
    if (ierr == 0) then
      !
      ! The open was successful, so close and delete the test file
      !
      close (iodir,status="delete",iostat=ierr)
      !
      ! Set the logical saying that we are using a separate CGNS directory
      !
      using_cgns_dir = true
      cgns_dir = trim(test_dir)
      !
    else
      !
      ! The open failed, so more than likely the directory CGNS_files
      ! doesnt exist. Try to use the intrinsic execute_command_line
      ! to create the directory.
      !
      command = "mkdir -p "//trim(test_dir)
      !
      call execute_command(command,wait_for_completion,exitstat, &
                           cmndstat,cmndmsg,pname)
     !call execute_command_line( trim(command) , wait_for_completion , &
     !                           exitstat , cmndstat , cmndmsg )
     !call write_command_results(pname,command,exitstat,cmndstat,cmndmsg)
      !
      if (exitstat /= 0 .or. cmndstat /= 0) then
        using_cgns_dir = fals
        cgns_dir = trim(output_dir)
      else
        using_cgns_dir = true
        cgns_dir = trim(test_dir)
      end if
      !
    end if
    !
#ifdef PBS_ENV
    !
    ! If we are using a separate CGNS directory and we are running
    ! this through PBS on NAS, try to modify the lfs striping for
    ! the CGNS directory for better I/O performance.
    !
    if (using_cgns_dir) then
      !
     !command = "lfs setstripe -c 64 -s 4m "//trim(cgns_dir)
      write (command,1) max(0,lustre_stripe_count), &
                        trim(adjustl(lustre_stripe_size)), &
                        trim(adjustl(cgns_dir))
      write (iout,2) trim(adjustl(command))
      !
      call execute_command(command,wait_for_completion,exitstat, &
                           cmndstat,cmndmsg,pname)
     !call execute_command_line( trim(command) , wait_for_completion , &
     !                           exitstat , cmndstat , cmndmsg )
     !call write_command_results(pname,command,exitstat,cmndstat,cmndmsg)
      !
      ! We arent going to test the status of executing this command because
      ! it worked or it didnt. Either way, we will continue on and output
      ! the CGNS solution to the CGNS directory.
      !
    end if
    !
#endif
    !
  end if
  !
  ! Broadcast the value of cgns_dir to all the non-root MPI processes
  !
  clen = len_trim(cgns_dir)
  call mpi_bcast(clen,1_int_mpi,mpi_inttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
  call mpi_bcast(cgns_dir(1:clen),int(clen,kind=int_mpi),MPI_CHARACTER, &
                 0_int_mpi,MPI_COMM_WORLD,mpierr)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("lfs setstripe -c ",i0," -S ",a,1x,a)
  2 format ("Lustre striping command = '",a,"'")
  !
end subroutine check_cgns_dir
!
!###############################################################################
!
subroutine read_grid_file()
  !
  !.. Use Statements ..
  use geovar,        only : nr,grid
  use geovar,        only : output_bface
  use module_cgns,   only : read_cgns_gridfile
  use module_gmsh,   only : read_gmsh_gridfile
  use module_plot3d, only : read_plot3d_gridfile
  use order_mod,     only : find_valid_geometries
  use ovar,          only : use_bc_conflicts
  use ovar,          only : itestcase
  use ovar,          only : postproc_debug_grid
  use ovar,          only : grid_format
  use ovar,          only : gridfile
  use ovar,          only : plot3d_cutfile
  use ovar,          only : plot3d_bcsfile
  use ovar,          only : solution_file
  use ovar,          only : correct_cgns_grid_points
  !
  use periodic_mod, only : match_periodic_faces
  !
  use initialization_mod, only : initialize_reference
  use initialization_mod, only : create_bc_in
  !
  use internal_grid_mod, only : create_internal_grid
  use internal_grid_mod, only : write_internal_grid
  !
  use eqn_idx, only : create_equation_indices
  !
  use synthetic_eddy_mod, only : separate_sem_boundary_faces
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_grid_file"
  !
  ! Routine that uses the variable grid_format to choose
  ! which routine to use in order to read the grid file.
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  call memory_pause("reading input file","reading grid file")
  if (grid_format == Internally_Generated_Grid .or. postproc_debug_grid) then
    call create_internal_grid
  else if (grid_format == Gambit_Neutral) then
    call read_gambit_neutral_gridfile
  else if (grid_format == Gmsh_Format) then
    call read_gmsh_gridfile( gridfile )
  else if (grid_format == CGNS_Format) then
    call read_cgns_gridfile( gridfile )
  else if (grid_format == Plot3D_Format) then
    call read_plot3d_gridfile( gridfile , solution_file , &
                               plot3d_cutfile , plot3d_bcsfile )
    call check_for_bad_number(__LINE__,__FILE__)
  end if
  !
  call output_bface(strng="after.grid_read.",only_glb_root=true)
  !
  ! Alter bface due to the input boundary conditions
  ! NOTE: This will apply the correct BC to bface(1,:) and fill bface(6,:)
  !
  call apply_bc_input_to_bface
  call output_bface(strng="after.apply_bc_input.",only_glb_root=true)
  !
  ! Create the grid derived type if it hasnt already been created.
  ! This means that the grid is linear.
  !
  if (.not. allocated(grid)) then
    if (i_am_host_root) then
      call create_linear_grid_dt
    end if
  end if
  !
  if (correct_cgns_grid_points) then
    call check_ho_points_are_equal_distance
  end if
  !
block
  use geovar, only : ncell,nfbnd
 !call write_internal_grid(ncell,nfbnd)
 !call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
end block
  !
  ! Match periodic faces specified using bc_input
  !
  call write_namelist
  call match_periodic_faces
  call output_bface(strng="after.match_periodic.",only_glb_root=true)
  !
  ! If running certain special cases and the grid was not internally generated,
  ! make sure the grid extents are valid for the current test case
  !
 !if (grid_format /= Internally_Generated_Grid) then
    if (any(abs(itestcase) == Transport_Problems) .or. &
        itestcase == Taylor_Green_Vortex) then
      call rescale_grid_for_vortex_problems
    end if
 !end if
  !
  ! Create the equation indices now that we know the dimension of the grid
  !
  call create_equation_indices(nr,itestcase)
  call check_for_bad_number(__LINE__,__FILE__)
  !
  ! Initialize the reference conditions
  !
  call initialize_reference
  call check_for_bad_number(__LINE__,__FILE__)
  !
  ! Before doing anything with the boundary conditions, check to see if
  ! there are any SEM boundary conditions. If there are any, read in the
  ! SEM profiles and separate any boundary faces outside of the profile
  ! height into new subsonic inflow boundaries.
  !
  call separate_sem_boundary_faces
  call output_bface(strng="after.separate_sem.",only_glb_root=true)
  !
  ! First, initialize the bc_in array
  ! NOTE: This will also fill bface(6,:) with the correct value
  !
  call create_bc_in
  call check_for_bad_number(__LINE__,__FILE__)
  call output_bface(strng="after.create_bc_in.",only_glb_root=true)
  !
  ! Reorder the boundary faces
  !
  call reorder_boundary_faces
  call check_for_bad_number(__LINE__,__FILE__)
  call output_bface(strng="after.reorder_bface.",only_glb_root=true)
  !
  ! Check that all the nodes are in counter-clockwise
  ! ordering for each cell and adjust if needed
  !
  call memory_pause("reading grid file","checking the node orderings")
  call check_node_ordering
  call check_for_bad_number(__LINE__,__FILE__)
  !
  ! Adjust the grid coordinates and nodes_of_cell
  ! connectivity to account for ghost cells
  !
  call memory_pause("checking the node orderings","adjusting for ghost cells")
  call adjust_for_ghost_cells
  call check_for_bad_number(__LINE__,__FILE__)
  !
  ! Find any grid nodes on the boundaries that separate two boundary
  ! faces with non-matching boundary conditions
  !
  if (use_bc_conflicts) then
    call memory_pause("adjusting for ghost cells", &
                      "finding boundary condition conflicts")
    call find_bc_conflicts
    call check_for_bad_number(__LINE__,__FILE__)
    call memory_pause("finding boundary condition conflicts", &
                      "finding valid geometries")
  else
    call memory_pause("adjusting for ghost cells","finding valid geometries")
  end if
  !
  ! Go through the grid and find which geometries are valid for
  ! the current simulation
  !
  call find_valid_geometries
  call check_for_bad_number(__LINE__,__FILE__)
  call memory_pause("finding valid geometries", &
                    "finish processing the remainder of the input file")
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_grid_file
!
!###############################################################################
!
subroutine write_namelist()
  use input_namelist_mod
#ifdef DEBUG_ON
!#define DUMP_NAMELIST_BEFORE_PERIODIC
#ifdef DUMP_NAMELIST_BEFORE_PERIODIC
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  if (mypnum == glb_root) then
    write (iout,104)
    write (iout,nml=input,delim='apostrophe')
    write (iout,105)
    flush (iout)
  end if
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  104 format (/,"#######################################",/, &
                "#####      INPUT PARAMETERS       #####",/, &
                "#######################################",/)
  105 format (/,"#######################################",/, &
                "#####   END OF INPUT PARAMETERS   #####",/, &
                "#######################################",/)
#endif
#endif
end subroutine write_namelist
!
!###############################################################################
!
subroutine check_ho_points_are_equal_distance()
  !
  !.. Use Statements ..
  use geovar, only : grid,ncell,nfbnd
  use module_cgns, only : reorder_to_cgns_pts
  use ovar,   only : cgns_grid_needs_correction
  !
  !.. Local Scalars ..
  integer :: nc,nf,i,j,k,nd,ncount
  integer :: this_geom,this_order
  !
  !.. Local Allocatable Arrays ..
  integer,     allocatable :: pts(:)
  integer,     allocatable :: pts1(:)
  integer,     allocatable :: pts2(:)
  integer,     allocatable :: pts3(:)
  real(wp),    allocatable :: xyz_pts(:,:)
  real(wp),    allocatable :: dsq(:,:,:)
  real(wp),    allocatable :: dsh(:,:,:,:)
  logical(lk), allocatable :: is_eq_dsq(:,:)
  logical(lk), allocatable :: is_eq_dsh(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_ho_points_are_equal_distance"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (mypnum == glb_root) then
    ncount = 0
    do nc = 1,ncell
      !
      this_geom = grid%elem(nc)%prop%geom
      this_order = grid%elem(nc)%prop%order
      !
      if (all(this_geom /= [Geom_Quad,Geom_Hexa])) cycle
      !
      pts = grid%elem(nc)%pts(:) ! F2003 AUTO-REALLOC
      ! F2003 AUTO-REALLOC
      pts1 = reorder_to_cgns_pts(grid%elem(nc)%prop%elem_type,pts)
      pts2 = reorder_to_cgns_pts(grid%elem(nc)%prop%elem_type, &
                                 intseq(1,size(pts)),from_cgns=true)
      !
      nd = geom_dimen(this_geom)
      !
      xyz_pts = grid%xyz(1:nd,pts(:)) ! F2003 AUTO-REALLOC
      !
      if (this_geom == Geom_Quad) then
        call check_ho_point_spacings_quad(this_order,xyz_pts,dsq,is_eq_dsq)
        if (.not. all(is_eq_dsq(:,:))) then
          ncount = ncount+1
#ifndef PBS_ENV
#ifdef DEBUG_ON
          write (iout,1) nc
#endif
#endif
        end if
      else
        call check_ho_point_spacings_hexa(this_order,xyz_pts,dsh, &
                                          is_eq_dsh,nc,pts2)
        if (.not. all(is_eq_dsh(:,:,:))) then
          ncount = ncount+1
#ifdef PBS_ENV
          if (ncount > 0) exit
#endif
#ifndef PBS_ENV
#ifdef DEBUG_ON
          write (iout,2) nc
#endif
#endif
        end if
      end if
      !
    end do
#ifdef PBS_ENV
    if (ncount > 0) write (iout,4)
#else
    if (ncount > 0) write (iout,3) ncount
#endif
  end if
  !
  if (ncpu > 1) then
    call mpi_bcast(ncount,1_int_mpi,mpi_inttyp,glb_root,MPI_COMM_WORLD,mpierr)
  end if
  !
  if (ncount > 0) then
    cgns_grid_needs_correction = true
  end if
  !
 !call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("The high-order points for quad cell #",i0, &
            " do not appear to be equally spaced!")
  2 format ("The high-order points for hex cell #",i0, &
            " do not appear to be equally spaced!")
  3 format (/,80("#"),/, &
              2x,"ERROR: Found ",i0," cells in which the high-order points", &
                 " do not appear",/, &
              9x,"to be equally spaced! This has the potential to introduce", &
                 " additional",/, &
              9x,"errors into the simulation because the mapping of the cell", &
                 " coordinates",/, &
              9x,"from physical space to computational space is based on the", &
                 " assumption",/, &
              9x,"that the input grid points in each cell are equally", &
                 " spaced.",/, &
              9x,"YOU HAVE BEEN WARNED!",/, &
              80("#"),/)
  !
  4 format (/,80("#"),/, &
    2x,"ERROR: At least once cell was found to have high-order points that",/, &
    9x,"do not appear to be equally spaced! This has the potential",/, &
    9x,"to introduce additional errors into the simulation because",/, &
    9x,"the mapping of the cell coordinates from physical space to",/, &
    9x,"computational space is based on the assumption that the",/, &
    9x,"input grid points in each cell are equally spaced. We will",/, &
    9x,"attempt a method to correct the grid points in each cell",/, &
    9x,"so that they will end up equally spaced, but this method",/, &
    9x,"is very experimental and unproven.",/, &
    9x,"YOU HAVE BEEN WARNED!",/, &
              80("#"),/)
  !
end subroutine check_ho_points_are_equal_distance
!
!###############################################################################
!
subroutine check_ho_point_spacings_quad(p,xyz,dsq,is_eq_dsq)
  !
  use ovar, only : epsval => tvrref
  !
  !.. Formal Arguments ..
  integer,                   intent(in) :: p
  real(wp),                  intent(in) :: xyz(:,:)
  real(wp),    allocatable, intent(out) :: dsq(:,:,:)
  logical(lk), allocatable, intent(out) :: is_eq_dsq(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,i1,i2,j1,j2,np
  real(wp) :: eps
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: ds(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_ho_point_spacings_hexa"
  !
continue
  !
  eps = epsval
  !
  np = p+1
  !
  dsq = init_array(p,np,2,zero) ! F2003 AUTO-REALLOC
  is_eq_dsq = init_array(np,2,fals) ! F2003 AUTO-REALLOC
  ds = init_array(p,zero) ! F2003 AUTO-REALLOC
  !
  ! Distances between points in i direction
  !
  do j = 1,np
    !
    do i = 1,p
      i1 = map2dto1d(i  ,j,np)
      i2 = map2dto1d(i+1,j,np)
      ds(i) = norm2( xyz(:,i2) - xyz(:,i1) )
      dsq(i,j,1) = ds(i)
    end do
    !
    is_eq_dsq(j,1) = all_are_equal( ds(:) , acceptable_error=eps )
    !
  end do
  !
  ! Distances between points in j direction
  !
  do i = 1,np
    !
    do j = 1,p
      j1 = map2dto1d(i,j  ,np)
      j2 = map2dto1d(i,j+1,np)
      ds(j) = norm2( xyz(:,j2) - xyz(:,j1) )
      dsq(j,i,2) = ds(j)
    end do
    !
    is_eq_dsq(i,2) = all_are_equal( ds(:) , acceptable_error=eps )
    !
  end do
  !
  contains
  !
  pure function map2dto1d(i,j,n) result(return_value)
    integer, intent(in) :: i,j,n
    integer :: return_value
    return_value = i + n*(j-1)
  end function map2dto1d
  !
end subroutine check_ho_point_spacings_quad
!
!###############################################################################
!
subroutine check_ho_point_spacings_hexa(p,xyz,dsh,is_eq_dsh,nc,pts)
  !
  use ovar,   only : epsval => tvrref
  use ovar,   only : out_eps => tinref
  use ovar,   only : gsf => grid_scaling_factor
  use geovar, only : bface
  !
  !.. Formal Arguments ..
  integer,                   intent(in) :: p
  real(wp),                  intent(in) :: xyz(:,:)
  real(wp),    allocatable, intent(out) :: dsh(:,:,:,:)
  logical(lk), allocatable, intent(out) :: is_eq_dsh(:,:,:)
  integer,                   intent(in) :: nc
  integer,                   intent(in) :: pts(:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,i1,i2,j1,j2,k1,k2,np
  real(wp) :: eps,dl,da,xmin,xmax,rmin,rmax
  !
  integer  :: l,m,n
  real(wp) :: x,y,z,r,t
  !
  real(wp) :: vdsn,vdsh,vdsr
  character(len=14) :: cdsn,cdsh,cdsr
  !
  character(len=15) :: cell_loc
  !
  real(wp) :: output_tolerance
  real(wp) :: dshn_tolerance
  real(wp) :: dsr_tolerance
  !
  logical(lk), save :: first_output_i = true
  logical(lk), save :: first_output_j = true
  logical(lk), save :: first_output_k = true
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: ds(:)
  real(wp), allocatable :: dsn(:,:,:,:)
  real(wp), allocatable :: dxh(:,:,:)
  real(wp), allocatable :: dxn(:,:,:)
  real(wp), allocatable :: dxr(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_ho_point_spacings_hexa"
  !
  logical(lk) :: debug_output_ho_point_spacings
#ifdef PBS_ENV
  logical(lk), parameter :: turn_off_output_override = true
#else
  logical(lk), parameter :: turn_off_output_override = true
 !logical(lk), parameter :: turn_off_output_override = fals
#endif
  !
continue
  !
  eps = epsval
 !write (iout,*) "eps = ",eps
  !
  dshn_tolerance = eps9
  dsr_tolerance = eps4
  !
  np = p+1
  !
  dsh = reshape( [zero] , shape=[p,np,np,3] , pad=[zero] ) ! F2003 AUTO-REALLOC
  dsn = reshape( [zero] , shape=[p,np,np,3] , pad=[zero] ) ! F2003 AUTO-REALLOC
  is_eq_dsh = init_array(np,np,3,fals) ! F2003 AUTO-REALLOC
  ds = init_array(p,zero) ! F2003 AUTO-REALLOC
  !
  dxh = init_array(np,np,3,zero) ! F2003 AUTO-REALLOC
  dxn = init_array(np,np,3,zero) ! F2003 AUTO-REALLOC
  dxr = init_array(np,np,3,zero) ! F2003 AUTO-REALLOC
  !
  debug_output_ho_point_spacings = fals
  xmin = minval(xyz(1,:))/gsf
  xmax = maxval(xyz(1,:))/gsf
  rmin = minval( sqrt( (xyz(2,:)**2) + (xyz(3,:)**2) ) )/gsf
  rmax = maxval( sqrt( (xyz(2,:)**2) + (xyz(3,:)**2) ) )/gsf
 !if (xmin > -3.45_wp .and. xmax < -3.25_wp .and. &
 !    rmin > 0.88_wp .and. rmax < one) then
  if (xmin > -3.55_wp .and. xmax < -3.15_wp .and. &
      rmax < one .and. any(bface(2,:)==nc)) then
      debug_output_ho_point_spacings = true
  end if
  !
  if (turn_off_output_override) then
    debug_output_ho_point_spacings = fals
  end if
  !
  ! Check all of the i-direction segments to see if the points
  ! along each segment are equally spaced from each other.
  !
  do k = 1,np
    do j = 1,np
      !
      dl = zero
      do i = 1,p
        i1 = map3dto1d(i  ,j,k,np)
        i2 = map3dto1d(i+1,j,k,np)
        ds(i) = norm2( xyz(:,i2) - xyz(:,i1) )
        dl = dl + ds(i)
      end do
      !
      da = dl / real(p,kind=wp)
      !
      dsh(1:p,j,k,1) = ds(1:p)
      !
      ds(1:p) = ds(1:p) / da
      dsn(1:p,j,k,1) = ds(1:p)
      !
      is_eq_dsh(j,k,1) = all_are_equal( ds(:) , acceptable_error=eps )
     !write (681,5) minval(ds(:)), maxval(ds(:)), dl, &
     !    minval(dsh(:,j,k,1)), maxval(dsh(:,j,k,1)), da, is_eq_dsh(j,k,1)
     !5 format (2(3(2x,es14.6),2x),2x,l1)
     !write (681,6) abs( maxval(ds(:)) - minval(ds(:)) ), &
     !    abs( maxval(dsh(:,j,k,1)) - minval(dsh(:,j,k,1)) ), is_eq_dsh(j,k,1)
     !6 format (2(2x,es14.6),2x,l1)
      !
      vdsn = abs( maxval(ds(:)) - minval(ds(:)) )
      vdsh = abs( maxval(dsh(:,j,k,1)) - minval(dsh(:,j,k,1)) )
      vdsr = abs( maxval(ds(:)) / minval(ds(:)) )
      !
      dxh(j,k,1) = vdsh
      dxn(j,k,1) = vdsn
      dxr(j,k,1) = vdsr
      !
      cdsn = "  ------------"
      cdsh = "  ------------"
      cdsr = "  ------------"
      !
      if (vdsn > dshn_tolerance)         write (cdsn,16) vdsn
      if (vdsh > dshn_tolerance)         write (cdsh,16) vdsh
      if (abs(vdsr-one) > dsr_tolerance) write (cdsr,16) vdsr
      !
      if (debug_output_ho_point_spacings) then
        write (681,17) adjustr(cdsn),adjustr(cdsh),adjustr(cdsr), &
                       is_eq_dsh(j,k,1)
      end if
      16 format (es14.6)
      17 format (3(2x,a14),2x,l1)
      !
    end do
  end do
  if (debug_output_ho_point_spacings) then
    write (681,*)
  end if
  !
  ! Check all of the j-direction segments to see if the points
  ! along each segment are equally spaced from each other.
  !
  do i = 1,np
    do k = 1,np
      !
      dl = zero
      do j = 1,p
        j1 = map3dto1d(i,j  ,k,np)
        j2 = map3dto1d(i,j+1,k,np)
        ds(j) = norm2( xyz(:,j2) - xyz(:,j1) )
        dl = dl + ds(j)
      end do
      !
      da = dl / real(p,kind=wp)
      !
      dsh(1:p,k,i,2) = ds(1:p)
      !
      ds(1:p) = ds(1:p) / da
      dsn(1:p,k,i,2) = ds(1:p)
      !
      is_eq_dsh(k,i,2) = all_are_equal( ds(:) , acceptable_error=eps )
     !write (682,5) minval(ds(:)), maxval(ds(:)), dl, &
     !    minval(dsh(:,k,i,2)), maxval(dsh(:,k,i,2)), da, is_eq_dsh(k,i,2)
     !write (682,6) abs( maxval(ds(:)) - minval(ds(:)) ), &
     !    abs( maxval(dsh(:,k,i,2)) - minval(dsh(:,k,i,2)) ), is_eq_dsh(k,i,2)
      !
      vdsn = abs( maxval(ds(:)) - minval(ds(:)) )
      vdsh = abs( maxval(dsh(:,k,i,2)) - minval(dsh(:,k,i,2)) )
      vdsr = abs( maxval(ds(:)) / minval(ds(:)) )
      !
      dxh(k,i,2) = vdsh
      dxn(k,i,2) = vdsn
      dxr(k,i,2) = vdsr
      !
      cdsn = "  ------------"
      cdsh = "  ------------"
      cdsr = "  ------------"
      !
      if (vdsn > dshn_tolerance)         write (cdsn,16) vdsn
      if (vdsh > dshn_tolerance)         write (cdsh,16) vdsh
      if (abs(vdsr-one) > dsr_tolerance) write (cdsr,16) vdsr
      !
      if (debug_output_ho_point_spacings) then
        write (682,17) adjustr(cdsn),adjustr(cdsh),adjustr(cdsr), &
                       is_eq_dsh(k,i,2)
      end if
      !
    end do
  end do
  if (debug_output_ho_point_spacings) then
    write (682,*)
  end if
  !
  ! Check all of the k-direction segments to see if the points
  ! along each segment are equally spaced from each other.
  !
  do j = 1,np
    do i = 1,np
      !
      dl = zero
      do k = 1,p
        k1 = map3dto1d(i,j,k  ,np)
        k2 = map3dto1d(i,j,k+1,np)
        ds(k) = norm2( xyz(:,k2) - xyz(:,k1) )
        dl = dl + ds(k)
      end do
      !
      da = dl / real(p,kind=wp)
      !
      dsh(1:p,i,j,3) = ds(1:p)
      !
      ds(1:p) = ds(1:p) / da
      dsn(1:p,i,j,3) = ds(1:p)
      !
      is_eq_dsh(i,j,3) = all_are_equal( ds(:) , acceptable_error=eps )
     !write (683,5) minval(ds(:)), maxval(ds(:)), dl, &
     !    minval(dsh(:,i,j,3)), maxval(dsh(:,i,j,3)), da, is_eq_dsh(i,j,3)
     !write (683,6) abs( maxval(ds(:)) - minval(ds(:)) ), &
     !    abs( maxval(dsh(:,i,j,3)) - minval(dsh(:,i,j,3)) ), is_eq_dsh(i,j,3)
      !
      vdsn = abs( maxval(ds(:)) - minval(ds(:)) )
      vdsh = abs( maxval(dsh(:,i,j,3)) - minval(dsh(:,i,j,3)) )
      vdsr = abs( maxval(ds(:)) / minval(ds(:)) )
      !
      dxh(i,j,3) = vdsh
      dxn(i,j,3) = vdsn
      dxr(i,j,3) = vdsr
      !
      cdsn = "  ------------"
      cdsh = "  ------------"
      cdsr = "  ------------"
      !
      if (vdsn > dshn_tolerance)         write (cdsn,16) vdsn
      if (vdsh > dshn_tolerance)         write (cdsh,16) vdsh
      if (abs(vdsr-one) > dsr_tolerance) write (cdsr,16) vdsr
      !
      if (debug_output_ho_point_spacings) then
        write (683,17) adjustr(cdsn),adjustr(cdsh),adjustr(cdsr), &
                       is_eq_dsh(i,j,3)
      end if
      !
    end do
  end do
  if (debug_output_ho_point_spacings) then
    write (683,*)
  end if
  !
  !=============================================================================
  !
  if (debug_output_ho_point_spacings) then
    !
    output_tolerance = out_eps
    !
    if (any(bface(2,:) == nc)) then
      cell_loc = "(Boundary Cell)"
    else
      cell_loc = "(Interior Cell)"
    end if
    !
    if (any(dxn(:,:,1) > output_tolerance)) then
      if (first_output_i) then
        write (691,11) nc, cell_loc
        first_output_i = fals
      else
        write (691,7) nc, cell_loc
      end if
      do k = 1,np
        do j = 1,np
          if (dxn(j,k,1) > output_tolerance) then
            i = 1
            n = map3dto1d(i,j,k,np)
            x = xyz(1,n)
            y = xyz(2,n)
            z = xyz(3,n)
            r = sqrt(y*y+z*z)
            t = atan2(z,y)
            t = (t + merge(two*pi,zero,t<zero)) * radians2degrees
            write (691,8) pts(n),x/gsf,y/gsf,z/gsf,r/gsf,t
            do l = 1,p
              i = l+1
              n = map3dto1d(i,j,k,np)
              x = xyz(1,n)
              y = xyz(2,n)
              z = xyz(3,n)
              r = sqrt(y*y+z*z)
              t = atan2(z,y)
              t = (t + merge(two*pi,zero,t<zero)) * radians2degrees
              write (691,9) pts(n),x/gsf,y/gsf,z/gsf,r/gsf,t, &
                            dsh(l,j,k,1),dsn(l,j,k,1)
            end do
          end if
        end do
      end do
      write (691,10)
    end if
    !
    if (any(dxn(:,:,2) > output_tolerance)) then
      if (first_output_j) then
        write (692,11) nc, cell_loc
        first_output_j = fals
      else
        write (692,7) nc, cell_loc
      end if
      do i = 1,np
        do k = 1,np
          if (dxn(k,i,2) > output_tolerance) then
            j = 1
            n = map3dto1d(i,j,k,np)
            x = xyz(1,n)
            y = xyz(2,n)
            z = xyz(3,n)
            r = sqrt(y*y+z*z)
            t = atan2(z,y)
            t = (t + merge(two*pi,zero,t<zero)) * radians2degrees
            write (692,8) pts(n),x/gsf,y/gsf,z/gsf,r/gsf,t
            do l = 1,p
              j = l+1
              n = map3dto1d(i,j,k,np)
              x = xyz(1,n)
              y = xyz(2,n)
              z = xyz(3,n)
              r = sqrt(y*y+z*z)
              t = atan2(z,y)
              t = (t + merge(two*pi,zero,t<zero)) * radians2degrees
              write (692,9) pts(n),x/gsf,y/gsf,z/gsf,r/gsf,t, &
                            dsh(l,k,i,2),dsn(l,k,i,2)
            end do
          end if
        end do
      end do
      write (692,10)
    end if
    !
    if (any(dxn(:,:,3) > output_tolerance)) then
      if (first_output_k) then
        write (693,11) nc, cell_loc
        first_output_k = fals
      else
        write (693,7) nc, cell_loc
      end if
      do j = 1,np
        do i = 1,np
          if (dxn(i,j,3) > output_tolerance) then
            k = 1
            n = map3dto1d(i,j,k,np)
            x = xyz(1,n)
            y = xyz(2,n)
            z = xyz(3,n)
            r = sqrt(y*y+z*z)
            t = atan2(z,y)
            t = (t + merge(two*pi,zero,t<zero)) * radians2degrees
            write (693,8) pts(n),x/gsf,y/gsf,z/gsf,r/gsf,t
            do l = 1,p
              k = l+1
              n = map3dto1d(i,j,k,np)
              x = xyz(1,n)
              y = xyz(2,n)
              z = xyz(3,n)
              r = sqrt(y*y+z*z)
              t = atan2(z,y)
              t = (t + merge(two*pi,zero,t<zero)) * radians2degrees
              write (693,9) pts(n),x/gsf,y/gsf,z/gsf,r/gsf,t, &
                            dsh(l,i,j,3),dsn(l,i,j,3)
            end do
          end if
        end do
      end do
      write (693,10)
    end if
    !
  end if
  !
  7 format (/,92("#"),/,2x,"Cell # ",i0,2x,a)
  8 format (92("-"),/,2x,i3,2x,5f11.6)
  9 format (2x,i3,2x,5f11.6,2(3x,f12.8))
  10 format (92("-"))
  11 format (56x,"Theta",/,2x,"Pt#",8x,"X",10x,"Y",10x,"Z",10x,"R", &
             8x,"Degrees",8x,"|ds|",9x,"|ds|/|da|",/, &
             92("#"),/,2x,"Cell # ",i0,2x,a)
  contains
  !
  pure function map3dto1d(i,j,k,n) result(return_value)
    integer, intent(in) :: i,j,k,n
    integer :: return_value
    return_value = i + n*(j-1) + n*n*(k-1)
  end function map3dto1d
  !
end subroutine check_ho_point_spacings_hexa
!
!###############################################################################
!
subroutine get_ho_segments()
  !
  !=============================================================================
  !
  integer, parameter :: p2_i(3,9) = reshape( [ 1, 9, 2, &
                                              12,21,10, &
                                               4,11, 3, &
                                              13,22,14, &
                                              25,27,23, &
                                              16,24,15, &
                                               5,17, 6, &
                                              20,26,18, &
                                               8,19, 7], shape(p2_i) )
  !
  integer, parameter :: p2_j(3,9) = reshape( transpose(p2_i) , shape(p2_j) )
  !
  integer, parameter :: p2_k(3,9) = reshape( transpose(p2_j) , shape(p2_k) )
  !
  !=============================================================================
  !
  integer, parameter :: p3_i(4,16) = reshape( [ 1, 9,10, 2, &
                                               16,33,34,11, &
                                               15,36,35,12, &
                                                4,14,13, 3, &
                                               17,37,38,19, &
                                               50,57,58,41, &
                                               49,60,59,42, &
                                               23,46,45,21, &
                                               18,40,39,20, &
                                               51,61,62,44, &
                                               52,64,63,43, &
                                               24,47,48,22, &
                                                5,25,26, 6, &
                                               32,53,54,27, &
                                               31,56,55,28, &
                                                8,30,29, 7], shape(p3_i) )
  !
  integer, parameter :: p3_j(4,16) = reshape( transpose(p3_i) , shape(p3_j) )
  !
  integer, parameter :: p3_k(4,16) = reshape( transpose(p3_j) , shape(p3_k) )
  !
  !=============================================================================
  !
  integer, parameter :: p4_i(5,25) = reshape( [  1,  9, 10, 11,  2, &
                                                20, 45, 46, 47, 12, &
                                                19, 52, 53, 48, 13, &
                                                18, 51, 50, 49, 14, &
                                                 4, 17, 16, 15,  3, &
                                                21, 54, 55, 56, 24, &
                                                83, 99,100,101, 63, &
                                                82,106,107,102, 64, &
                                                81,105,104,103, 65, &
                                                30, 74, 73, 72, 27, &
                                                22, 61, 62, 57, 25, &
                                                84,108,109,110, 70, &
                                                89,115,116,111, 71, &
                                                88,114,113,112, 66, &
                                                31, 75, 80, 79, 28, &
                                                23, 60, 59, 58, 26, &
                                                85,117,118,119, 69, &
                                                86,124,125,120, 68, &
                                                87,123,122,121, 67, &
                                                32, 76, 77, 78, 29, &
                                                 5, 33, 34, 35,  6, &
                                                44, 90, 91, 92, 36, &
                                                43, 97, 98, 93, 37, &
                                                42, 96, 95, 94, 38, &
                                                 8, 41, 40, 39,  7],shape(p4_i))
  !
  integer, parameter :: p4_j(5,25) = reshape( transpose(p4_i) , shape(p4_j) )
  !
  integer, parameter :: p4_k(5,25) = reshape( transpose(p4_j) , shape(p4_k) )
  !
end subroutine get_ho_segments
!
!###############################################################################
!
subroutine apply_bc_input_to_bface()
  !
  !.. Use Statements ..
  use geovar, only : bface
  use ovar,   only : bc_input,bc_names
  !
  !.. Local Scalars ..
  integer :: n,nf
  !
  character(len=32) :: bcnam,bcinp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "apply_bc_input_to_bface"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  do n = 1,size(bc_input)
    !
    bface_loop: do nf = 1,size(bface,dim=2)
      !
      ! NOTE: This check on bc_comm will only cycle periodic boundary faces
      !       that have already had bface(1,nf) set to bc_periodic. This will
      !       only happen for
      !         1) most of the internal grid cases (TGV, TCF, vortex, etc.)
      !         2) when the string giving the BC name in either the CGNS grid
      !            file or the Plot3D BCS file contains the substring
      !            'periodic', e.g., 'x-periodic-1'. (This was the old way of
      !            specifying periodic BCs, and the new way of using
      !            bc_input to specify periodic BCs is preferred.)
      !
      ! NOTE: The new way of specifying periodic BCs using bc_input(:)%periodic
      !       will not yet have set bface(1,nf)=bc_periodic, so this cycle
      !       statement will get ignored.
      !
      if (any(bface(1,nf) == bc_comm) .or. bface(5,nf) == 0) cycle bface_loop
      !
      bcnam = trim(adjustl( uppercase( bc_names(bface(5,nf)) ) ))
      bcinp = trim(adjustl( uppercase( bc_input(n)%name      ) ))
      !
      if (bcnam == bcinp) then
        !
        ! Save the index for this boundary condition
        ! within the bc_in array to bface(6,:)
        !
        bface(6,nf) = n
        !
        ! Also save this index to bface(11,:) if the first dimension of
        ! bface is large enough. Its unsure if the first dimension of bface
        ! will be permanently increased to 11, but if it is, bface(11,:)
        ! is meant to store this value.
        !
        if (size(bface,dim=1) >= 11) then
          bface(11,nf) = n
        end if
        !
        ! Force the boundary condition to that specified in the input file
        ! if the value of set_bc_value is not equal to bc_unknown
        !
        if (bc_input(n)%set_bc_value /= bc_unknown) then
          bface(1,nf) = bc_input(n)%set_bc_value
        end if
        !
      end if
      !
    end do bface_loop
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine apply_bc_input_to_bface
!
!###############################################################################
!
subroutine reorder_boundary_faces()
  !
  !.. Use Statements ..
  use geovar, only : ncell,nfbnd,bface
  use geovar, only : num_bc_grps
  use geovar, only : grid,grid_t
  !
  !.. Local Scalars ..
  integer :: i,j,k,n,nc,nf,ngrps,ierr
  integer :: min_grp,max_grp,min_bcin,max_bcin
  !
  logical(lk) :: contains_periodic_faces
  logical(lk) :: errors_found
  !
  !.. Local Allocatable Arrays ..
  logical(lk),  allocatable :: msk(:)
  integer,      allocatable :: face_grps(:,:,:)
  integer,      allocatable :: sorted_grps(:,:)
  integer,      allocatable :: tmp_bface(:,:)
  integer,      allocatable :: new2old(:)
  integer,      allocatable :: old2new(:)
  type(grid_t), allocatable :: tmp_grid
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reorder_boundary_faces"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  contains_periodic_faces = fals
  errors_found = fals
  !
  ! Find all the physical boundary faces
  !
  allocate ( msk(1:nfbnd) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",1,__LINE__,__FILE__,ierr,error_message)
  !
  do nf = 1,nfbnd
    if (bface(1,nf) /= bc_periodic .and. bface(1,nf) /= bc_cpu_bnd) then
      msk(nf) = true
    end if
    if (bface(1,nf) == bc_periodic) then
      contains_periodic_faces = true
    end if
  end do
  !
  ! Find the maximum BC group number
  !
  if (count(msk) == 0) then
    min_grp = 1
    max_grp = 0
    min_bcin = 1
    max_bcin = 0
  else
    min_grp = minval( bface(5,1:nfbnd) , mask=msk )
    max_grp = maxval( bface(5,1:nfbnd) , mask=msk )
    min_bcin = minval( bface(6,1:nfbnd) , mask=msk )
    max_bcin = maxval( bface(6,1:nfbnd) , mask=msk )
  end if
  !
  allocate ( face_grps(1:2,min_bcin:max_bcin,min_grp:max_grp) , source=0 , &
             stat=ierr , errmsg=error_message )
  !
  face_grps(1,:,:) = 0
  face_grps(2,:,:) = bc_unknown
  !
  do nf = 1,nfbnd
    !
    if (.not. msk(nf)) cycle
    !
    j = bface(5,nf)
    i = bface(6,nf)
    !
    face_grps(1,i,j) = face_grps(1,i,j) + 1
    !
    ! Make sure each BC group / BC flow conditions combo has the same BC type
    !
    if (face_grps(2,i,j) == bc_unknown) then
      face_grps(2,i,j) = bface(1,nf)
    else if (face_grps(2,i,j) /= bface(1,nf)) then
      errors_found = true
      if (mypnum == glb_root) then
        write (iout,1) nf,bface(1,nf),i,j
      end if
    end if
    !
  end do
  !
  deallocate ( msk , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
  !
  if (errors_found) then
    write (error_message,2)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ngrps = count( face_grps(1,:,:) > 0 )
  if (contains_periodic_faces) ngrps = ngrps + 1
  !
  num_bc_grps = ngrps
  !
  flush (iout)
  if (mypnum == glb_root) then
    write (iout,100) pname,min_grp,max_grp,min_bcin,max_bcin,num_bc_grps
    100 format (/,2x,"In procedure ",a," : ",/, &
                  6x,"   min_grp  = ",i0,/, &
                  6x,"   max_grp  = ",i0,/, &
                  6x,"   min_bcin = ",i0,/, &
                  6x,"   max_bcin = ",i0,/, &
                  6x,"num_bc_grps = ",i0,/)
  end if
  flush (iout)
  !
  allocate ( sorted_grps(1:2,1:ngrps) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sorted_grps",1,__LINE__,__FILE__,ierr,error_message)
  !
  n = 0
  !
  ! First, find all the wall boundary groups
  !
  do k = 1,size(bc_walls)
    do j = lbound(face_grps,dim=3),ubound(face_grps,dim=3)
      do i = lbound(face_grps,dim=2),ubound(face_grps,dim=2)
        if (face_grps(2,i,j) == bc_walls(k)) then
          n = n + 1
          sorted_grps(1,n) = j
          sorted_grps(2,n) = i
          face_grps(2,i,j) = bc_unknown
        end if
      end do
    end do
  end do
  !
  ! Then find all the now-wall boundary groups
  !
  do j = lbound(face_grps,dim=3),ubound(face_grps,dim=3)
    do i = lbound(face_grps,dim=2),ubound(face_grps,dim=2)
      if (face_grps(2,i,j) /= bc_unknown) then
        n = n + 1
        sorted_grps(1,n) = j
        sorted_grps(2,n) = i
        face_grps(2,i,j) = bc_unknown
      end if
    end do
  end do
  !
  deallocate ( face_grps , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_grps",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! new2old(nf) = the old face number for a  new face nf
  ! old2new(nf) = the new face number for an old face nf
  !
  allocate ( new2old(1:nfbnd) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"new2old",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( old2new(1:nfbnd) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"old2new",1,__LINE__,__FILE__,ierr,error_message)
  !
  i = 0
  do n = 1,ngrps
    !
    if (n == ngrps .and. contains_periodic_faces) then
      !
      do nf = 1,nfbnd
        if (bface(1,nf) == bc_periodic) then
          i = i + 1
          new2old(i) = nf
          old2new(nf) = i
        end if
      end do
      !
    else
      !
      do nf = 1,nfbnd
        if (bface(5,nf) == sorted_grps(1,n) .and. &
            bface(6,nf) == sorted_grps(2,n)) then
          i = i + 1
          new2old(i) = nf
          old2new(nf) = i
        end if
      end do
      !
    end if
    !
  end do
  !
  deallocate ( sorted_grps , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"sorted_grps",2,__LINE__,__FILE__,ierr,error_message)
  !
  tmp_bface = bface ! F2003 AUTO-REALLOCATION
  !
  do nf = 1,nfbnd
    !
    tmp_bface(:,nf) = bface(:,new2old(nf))
    !
    if (tmp_bface(1,nf) == bc_periodic) then
      tmp_bface(5,nf) = old2new( bface(5,new2old(nf)) )
    end if
    !
  end do
  !
  call move_alloc( from=tmp_bface , to=bface )
  !
  deallocate ( new2old , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"new2old",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( old2new , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"old2new",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Copy all of grid into tmp_grid
  !
  if (i_am_host_root) then
    !
    tmp_grid = grid ! F2003 AUTO-REALLOCATION
    !
    nc = ncell
    !
    do nf = 1,nfbnd
      !
      nc = nc + 1
      !
      ! F2003 AUTO-REALLOCATION OF COMPONENT ARRAYS
      tmp_grid%elem(nc) = grid%elem( bface(10,nf) )
      !
      bface(10,nf) = nc
      !
    end do
    !
    call move_alloc( from=tmp_grid , to=grid )
    !
  end if
  !
  ! Broadcast the new bface to all the non-host-root MPI processes
  !
  if (ncpu > 1) then
    call mpi_bcast(bface,size(bface,kind=int_mpi),mpi_inttyp, &
                   host_root,my_host_comm,mpierr)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("bface(1,",i0,") = ",i0,"; bc_in(",i0,"); group # = ",i0)
  2 format ("ERRORS DETECTED WHILE TRYING TO RESORT THE BOUNDARY FACES")
  !
end subroutine reorder_boundary_faces
!
!###############################################################################
!
subroutine read_gambit_neutral_gridfile()
  !
  ! Routine to read in a grid file using Gambit Neutral formatting.
  !
  use geovar, only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geovar, only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geovar, only : cell_geom,cell_order
  use ovar,   only : loc_solution_pts,loc_flux_pts
  use ovar,   only : gridfile,grid_scaling_factor
  use order_mod, only : n_order
  !
  !.. Local Scalars ..
  integer :: i,j,k,n,nf,n1,n2,i1,nc,pc,pf
  integer :: icell,inode,iside,ierr,itype
  integer :: ngrps,nbsets,ndfvl,nentry,nvalues
  integer :: iogrd
  !
  character(len=80) :: text
  !
  !.. Local Allocatable Arrays ..
  logical(lk), allocatable, dimension(:) :: bc_done
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_gambit_neutral_gridfile"
  integer, dimension(1:7), parameter :: convert_geom = [Geom_Edge,Geom_Quad, &
                                                        Geom_Tria,Geom_Hexa, &
                                                        Geom_Pris,Geom_Tetr, &
                                                        Geom_Pyra]
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Open the grid file
  !
  open (newunit=iogrd,file=gridfile,status="old",action="read", &
                      position="rewind",iostat=ierr,iomsg=error_message)
  call io_error(pname,trim(adjustl(gridfile)),1,__LINE__,__FILE__, &
                ierr,error_message)
  !
  ! Skip the file header
  !
  do i = 1,6
    read (iogrd,1) text
  end do
  !
  ! Read the grid size and dimensions
  !
  read (iogrd,*) nnode,ncell,ngrps,nbsets,nr,ndfvl
  !
  allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  allocate ( cell_geom(1:ncell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_geom",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_order(1:ncell) , source=n_order , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_order",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Skip some lines
  !
  read (iogrd,1) text
  read (iogrd,1) text
  !
  ! Skip over the node coordinates
  !
  do i = 1,nnode
    read (iogrd,1) text
  end do
  !
  ! Skip some more lines
  !
  read (iogrd,1) text
  read (iogrd,1) text
  !
  ! Read through the cell connectivities to
  ! find the number of nodes defining each cell
  !
  do i = 1,ncell
    read (iogrd,*) n,cell_geom(i),nodes_of_cell_ptr(n+1)
  end do
  !
  ! Convert the cell geometry values from the Gambit neutral grid file
  !
  !   If 'igeom' is the value given in the grid file
  !       igeom = 1  -> edge
  !       igeom = 2  -> quadrilateral
  !       igeom = 3  -> triangle
  !       igeom = 4  -> brick/hexahedron
  !       igeom = 5  -> wedge/prism
  !       igeom = 6  -> tetrahedron
  !       igeom = 7  -> pyramid
  !
  !   We need them stored as
  !       cell_geom(i) = 1  -> edge
  !       cell_geom(i) = 2  -> quadrilateral
  !       cell_geom(i) = 3  -> triangle
  !       cell_geom(i) = 4  -> tetrahedron
  !       cell_geom(i) = 5  -> pyramid
  !       cell_geom(i) = 6  -> wedge/prism
  !       cell_geom(i) = 8  -> brick/hexahedron
  !
  do i = 1,ncell
    cell_geom(i) = convert_geom( cell_geom(i) )
  end do
  !
  ! Reshuffle nodes_of_cell_ptr
  !
  do i = 2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(i) = nodes_of_cell_ptr(i) + nodes_of_cell_ptr(i-1)
  end do
  !
  ! Skip the group sections
  !
  read (iogrd,1) text
  do i = 1,ngrps
    read (iogrd,1) text
    do
      read (iogrd,1) text
      if (trim(adjustl(text)) == "ENDOFSECTION") exit
    end do
  end do
  !
  ! Cycle through the different boundary conditions
  ! to find the total number of boundary faces
  !
  nfbnd = 0
  do i = 1,nbsets
    read (iogrd,1) text
    read (iogrd,*) text,itype,nentry,nvalues
    do j = 1,nentry
      read (iogrd,1) text
    end do
    read (iogrd,1) text
    nfbnd = nfbnd + nentry
  end do
  !
  ! Now go back to the beginning and read in the grid data
  !
  rewind (iogrd)
  !
  ! Allocate the arrays needed to store this data
  !
  allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_nodes",1,__LINE__,__FILE__,ierr,error_message)
  allocate ( nodes_of_cell(1:last(nodes_of_cell_ptr)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr,error_message)
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
  allocate ( bc_done(1:nfbnd) , source=true , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bc_done",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Skip the header and other data that is already known
  !
  do i = 1,9
    read (iogrd,1) text
  end do
  !
  ! Read in the nodal coordinates
  !
  do i = 1,nnode
    read (iogrd,*) n,(xyz_nodes(k,n),k=1,nr)
  end do
  !
  ! Apply the grid scaling factor to the coordinates of the grid nodes
  !
  do n = 1,size(xyz_nodes,dim=2)
    do i = 1,size(xyz_nodes,dim=1)
      xyz_nodes(i,n) = grid_scaling_factor*xyz_nodes(i,n)
    end do
  end do
  !
  ! Skip a couple lines
  !
  read (iogrd,1) text
  read (iogrd,1) text
  !
  ! Read the nodes defining each cell into nodes_of_cell
  !
  do i = 1,ncell
    read (iogrd,*) n,itype,j,(nodes_of_cell(k),k=nodes_of_cell_ptr(n)+1, &
                                                 nodes_of_cell_ptr(n+1))
  end do
  !
  ! Skip the group sections
  !
  read (iogrd,1) text
  do i = 1,ngrps
    read (iogrd,1) text
    do
      read (iogrd,1) text
      if (trim(adjustl(text)) == "ENDOFSECTION") exit
    end do
  end do
  !
  ! Read in all the boundary condition sets
  !
  n = 0
  do i = 1,nbsets
    !
    read (iogrd,1) text
    read (iogrd,*) text,itype,nentry,nvalues
    !
    ! Convert the BC string to itts respective integer value
    !
    itype = bc_string_to_integer(text)
    !
    ! Mark these boundaries as not finished if this is a periodic boundary
    !
    if (itype == bc_periodic) then
      bc_done(n+1:n+nentry) = fals
    end if
    !
    do j = 1,nentry
      n = n + 1
      bface(1,n) = itype
      read (iogrd,*) (bface(k,n),k=2,4),(bface(4+k,n),k=1,nvalues)
      if (itype /= bc_periodic) then
        bface(5,n) = i ! boundary condition group for this boundary face
      end if
    end do
    !
    read (iogrd,1) text
    !
  end do
  !
  close (iogrd,iostat=ierr,iomsg=error_message)
  call io_error(pname,trim(adjustl(gridfile)),2,__LINE__,__FILE__, &
                ierr,error_message)
  !
  ! Convert the cell geometry values from the Gambit neutral grid file
  !
  do nf = 1,nfbnd
    bface(3,nf) = convert_geom( bface(3,nf) )
  end do
  !
  ! Save the index for the ghost cell on the face in bface(10,:)
  !
  bface(10,:) = ncell + intseq(1,nfbnd)
  !
  ! Lets do a quick modification to the periodic boundary conditions
  ! For each periodic boundary condition, we will find the opposite
  ! connecting periodic boundary condition and link the two together
  !
  do nf = 1,nfbnd
    !
    ! Skip the boundary faces that are not periodic or have been completed
    !
    if (bc_done(nf)) cycle
    !
    ! Cells that are included in this periodic bc
    !
    nc = bface(2,nf)
    pc = bface(5,nf)
    !
    ! Find the other periodic boundary face with these two cells
    !
    search_loop: do n = 1,nfbnd
      if (bc_done(n)) cycle
      if (bface(2,n)==pc .and. bface(5,n)==nc) then
        pf = n
        exit search_loop
      end if
    end do search_loop
    !
    ! Link the two boundary faces to each other
    !
    bface(5,nf) = pf          ! boundary face of connecting cell
    bface(6,nf) = pc          ! connecting cell
    bface(8,nf) = bface(3,pf) ! connecting cell geometry type
    bface(9,nf) = bface(4,pf) ! face of connecting cell on boundary
    !
    bface(5,pf) = nf          ! boundary face of connecting cell
    bface(6,pf) = nc          ! connecting cell
    bface(8,pf) = bface(3,nf) ! connecting cell geometry type
    bface(9,pf) = bface(4,nf) ! face of connecting cell on boundary
    !
    ! Mark these two periodic boundary faces as finished
    !
    bc_done(nf) = true
    bc_done(pf) = true
    !
  end do
  !
  deallocate ( bc_done , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bc_done",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Make sure that Lobatto nodes are requested if the grid contains triangle
  ! elements since Gauss nodes is not currently set up for triangles.
  !
  if (any(cell_geom == Geom_Tria)) then
    write (*,3)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  !
  if (loc_solution_pts /= Legendre_Gauss_Lobatto .or. &
      loc_flux_pts /= Legendre_Gauss_Lobatto) then
    if (any(cell_geom == Geom_Tria)) then
      write (*,2)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
    end if
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format statements
  !
  1 format (a)
  2 format (/," Quadrature points besides Gauss-Lobatto nodes are not",/, &
              " currently supported for a grid containing triangle cells.",//)
  3 format (///," Triangles are temporarily disabled due to significant",/, &
                " changes in the quadrilateral part of the code. These ",/, &
                " changes have yet to be applied to triangles.",///)
  !
end subroutine read_gambit_neutral_gridfile
!
!###############################################################################
!
subroutine rescale_grid_for_vortex_problems()
  !
  !.. Use Statements ..
  use ovar,   only : itestcase
  use ovar,   only : Lref,VortexGrid_Lref
  use geovar, only : nr,grid,xyz_nodes,xyz
  !
  !.. Local Scalars ..
  integer  :: l
  real(wp) :: llenmax
  !
  logical(lk) :: grid_xyz_was_allocated
  logical(lk) :: xyz_was_allocated
  !
  !.. Local Arrays ..
  real(wp) :: lmin(1:nr)
  real(wp) :: lmax(1:nr)
  real(wp) :: llen(1:nr)
  real(wp) :: lave(1:nr)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: tmp_grid_xyz(:,:)
  real(wp), allocatable :: tmp_xyz(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "rescale_grid_for_vortex_problems"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  Lref = one
  !
  ! Make sure VortexGrid_Lref for the Taylor Green vortex problem is pi
  !
  if (itestcase == Taylor_Green_Vortex) then
    VortexGrid_Lref = PI
  end if
  !
  ! Get the min/max grid extents in each direction
  ! NOTE: If grid%xyz exists, xyz_nodes is a subset of grid%xyz and (in nearly
  !       all cases) will include the min/max grid coordinates. The only case
  !       where this might not be the case is with a highly curved grid, and
  !       this shouldnt be the case for either the isentropic vortex or Taylor
  !       Green vortex cases.
  !
  do l = 1,nr
    lmin(l) = minval(xyz_nodes(l,:))
    lmax(l) = maxval(xyz_nodes(l,:))
    llen(l) = half*(lmax(l)-lmin(l)) ! half length in each direction
    lave(l) = half*(lmax(l)+lmin(l))
  end do
  !
  ! Get the half length for the coordinate with the largest extent
  !
  llenmax = maxval(llen)
  !
  ! Create temporary arrays for grid%xyz and xyz so we dont have to check
  ! their allocation status for every adjustment for the grid
  ! NOTE: xyz shouldnt exist yet, but we will scale it just in case
  !
  tmp_grid_xyz = init_array(nr,1,zero) ! F2003 AUTO-REALLOC
  tmp_xyz = init_array(nr,1,zero)      ! F2003 AUTO-REALLOC
  !
  grid_xyz_was_allocated = fals
  if (allocated(grid)) then
    if (allocated(grid%xyz)) then
      grid_xyz_was_allocated = true
      call move_alloc( from=grid%xyz , to=tmp_grid_xyz )
    end if
  end if
  !
  xyz_was_allocated = fals
  if (allocated(xyz)) then
    xyz_was_allocated = true
    call move_alloc( from=xyz , to=tmp_xyz )
  end if
  !
  ! Center the grid coordinates, and scale them so the coordinate
  ! direction with the largest extent is scaled to (-1,1)
  !
  do l = 1,nr
    xyz_nodes(l,:) = (xyz_nodes(l,:) - lave(l)) / llenmax
    tmp_grid_xyz(l,:) = (tmp_grid_xyz(l,:) - lave(l)) / llenmax
    tmp_xyz(l,:) = (tmp_xyz(l,:) - lave(l)) / llenmax
  end do
  !
  ! Finally scale the grid so the coordinate direction with the
  ! largest extent is scaled to (-VortexGrid_Lref,VortexGrid_Lref)
  !
  xyz_nodes(:,:) = xyz_nodes(:,:) * abs(VortexGrid_Lref)
  tmp_grid_xyz(:,:) = tmp_grid_xyz(:,:) * abs(VortexGrid_Lref)
  tmp_xyz(:,:) = tmp_xyz(:,:) * abs(VortexGrid_Lref)
  !
  ! Move the temporary arrays for grid%xyz and xyz back into the actual
  ! arrays. We only want to do this if the arrays were already allocated
  ! before making the adjustments to the grid coordinates. Otherwise,
  ! we dont want to copy the (nr,1) sized temporary arrays over because
  ! we want to leave these array unallocated as they were when entering
  ! this subroutine.
  !
  if (grid_xyz_was_allocated) then
    call move_alloc( from=tmp_grid_xyz , to=grid%xyz )
  end if
  !
  if (xyz_was_allocated) then
    call move_alloc( from=tmp_xyz , to=xyz )
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine rescale_grid_for_vortex_problems
!
!###############################################################################
!
subroutine create_linear_grid_dt()
  !
  !.. Use Statements ..
  use geovar, only : ncell,nfbnd,grid,xyz_nodes,cell_geom
  use geovar, only : bface,nodes_of_cell_ptr,nodes_of_cell
  !
  !.. Local Scalars ..
  integer :: n,nc,nf,np,n1,n2,ierr
  integer :: host_cell,host_geom,host_side
  integer :: nfnods,face_geom
  !
  character(len=100) :: array_name
  !
  !.. Local Arrays ..
  integer :: ip(0:8)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_linear_grid_dt"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  allocate ( grid , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid",1,__LINE__,__FILE__,ierr,error_message)
  !
  grid%xyz = xyz_nodes(:,:) ! F2003 AUTO-REALLOC
  !
  allocate ( grid%elem(1:ncell+nfbnd) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid%elem",1,__LINE__,__FILE__,ierr,error_message)
  !
  do nc = 1,ncell
    !
    current_cell = nc
    !
    n1 = nodes_of_cell_ptr(nc)+1
    n2 = nodes_of_cell_ptr(nc+1)
    np = n2-n1+1
    !
    grid%elem(nc)%prop%npts = np
    grid%elem(nc)%prop%order = 1
    grid%elem(nc)%prop%geom = cell_geom(nc)
    grid%elem(nc)%prop%elem_type = cell_geom(nc)
    !
    grid%elem(nc)%pts = nodes_of_cell(n1:n2)   ! USING F2003 AUTO-REALLOCATION
    grid%elem(nc)%nodes = nodes_of_cell(n1:n2) ! USING F2003 AUTO-REALLOCATION
    !
    grid%elem(nc)%tags = [0,0] ! F2003 AUTO-REALLOCATION
    !
    grid%elem(nc) = reorder_linear_pts(grid%elem(nc))
    !
  end do
  !
  do nf = 1,nfbnd
    !
    nc = ncell + nf
    !
    current_cell = nc
    !
    host_cell = bface(2,nf)
    host_geom = bface(3,nf)
    host_side = bface(4,nf)
    nfnods = num_face_nodes(host_side,host_geom)
    !
    n1 = nodes_of_cell_ptr(host_cell)+1
    n2 = nodes_of_cell_ptr(host_cell+1)
    !
    ip(1:nfnods) = side_nodes(1:nfnods,host_side,host_geom)
    !
    face_geom = cell_face_geom(host_side,host_geom)
    !
    grid%elem(nc)%prop%npts = nfnods
    grid%elem(nc)%prop%order = 1
    grid%elem(nc)%prop%geom = face_geom
    grid%elem(nc)%prop%elem_type = face_geom
    !
    ! USING F2003 AUTO-REALLOCATION
    grid%elem(nc)%pts = nodes_of_cell( n1 + ip(1:nfnods) )
    grid%elem(nc)%nodes = nodes_of_cell( n1 + ip(1:nfnods) )
    !
    grid%elem(nc)%host_cell = host_cell
    !
    grid%elem(nc)%tags = [bface(1,nf),0] ! F2003 AUTO-REALLOCATION
    !
    grid%elem(nc) = reorder_linear_pts(grid%elem(nc))
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  1 format ("grid%elem(",i0,")%",a)
  !
end subroutine create_linear_grid_dt
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function reorder_linear_pts(elem,reverse) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : elem_t
  !
  !.. Formal Arguments ..
  type(elem_t), intent(in) :: elem
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: reverse
  !
  !.. Function Result ..
  type(elem_t) :: return_value
  !
  !.. Local Scalars ..
  integer :: n,n1
  logical(lk) :: reverse_reorder
  !
  !.. Local Parameters ..
  integer, parameter :: quad(4) = [1,2,4,3]
  integer, parameter :: pyra(5) = [1,2,4,3,5]
  integer, parameter :: hexa(8) = [1,2,4,3,5,6,8,7]
  !
continue
  !
  n1 = lbound(elem%pts,dim=1)
  !
  return_value = elem
  !
  reverse_reorder = fals
  if (present(reverse)) then
    reverse_reorder = reverse
  end if
  !
  if (reverse_reorder) then
    select case (elem%prop%geom)
      case (Geom_Quad)
        return_value%pts( quad(:) ) = elem%pts(:)
      case (Geom_Pyra)
        return_value%pts( pyra(:) ) = elem%pts(:)
      case (Geom_Hexa)
        return_value%pts( hexa(:) ) = elem%pts(:)
    end select
  else
    select case (elem%prop%geom)
      case (Geom_Quad)
        return_value%pts(:) = elem%pts( n1-1 + quad(:) )
      case (Geom_Pyra)
        return_value%pts(:) = elem%pts( n1-1 + pyra(:) )
      case (Geom_Hexa)
        return_value%pts(:) = elem%pts( n1-1 + hexa(:) )
    end select
  end if
  !
end function reorder_linear_pts
!
!###############################################################################
!
subroutine tecplot_rotations(node_order)
  !
  !.. Formal Arguments ..
  integer, dimension(:,:), intent(in) :: node_order
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,n,nn,m,ierr,iotec
  !
  !.. Local Arrays ..
  integer, dimension(3,8) :: coord
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "tecplot_rotations"
  character(len=*), parameter :: fname = "rotations.tec.dat"
  !
  integer, parameter :: base(3,8) = reshape( [0,0,0, 1,0,0, 1,1,0, 0,1,0, &
                                              0,0,1, 1,0,1, 1,1,1, 0,1,1] , &
                                             shape(base) )
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  open (newunit=iotec,file=fname,status="replace",action="write", &
                      form="formatted",position="rewind", &
                      iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (iotec,1)
  !
  do nn = 1,size(node_order,dim=2)
    !
    j = (nn-1)/4
    i = mod(nn+3,4)
    !
    coord = base
    coord(1,:) = coord(1,:) + 2*i
    coord(3,:) = coord(3,:) + 2*j
    !
    write (iotec,2) nn
    do k = 1,size(coord,dim=2)
      write (iotec,4) (base(l,node_order(k,nn)), l=1,size(coord,dim=1)), &
                      node_order(k,nn), k
    end do
   !write (iotec,4) (node_order(k,nn), k=1,size(node_order,dim=1))
    write (iotec,4) (k, k=1,size(node_order,dim=1))
    !
    write (iotec,3) "r",nn,4*(nn-1)+1
    write (iotec,4) node_order(1,nn),node_order(2,nn)
    write (iotec,3) "s",nn,4*(nn-1)+1
    write (iotec,4) node_order(1,nn),node_order(4,nn)
    write (iotec,3) "t",nn,4*(nn-1)+1
    write (iotec,4) node_order(1,nn),node_order(5,nn)
    !
  end do
  !
  close (iotec,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format('VARIABLES = "X", "Y", "Z", "Node", "Index"')
  2 format('ZONE T="Rotation ',i0,'", N=8, E=1, ', &
           'DATAPACKING=POINT, ZONETYPE=FEBRICK')
  3 format('ZONE T="',a,' - R(',i0,')", N=8, E=1, ', &
           'DATAPACKING=POINT, ZONETYPE=FELINESEG',/, &
           'VARSHARELIST=([1-5]=',i0,')')
  4 format(100(1x,i0))
 !4 format(*(1x,i0))
  !
end subroutine tecplot_rotations
!
!###############################################################################
!
subroutine tecplot_cell_rotations(nodes_of_cell_ptr,nodes_of_cell, &
                                  old_noc,xyz_nodes,rotation,node_order)
  !
  !.. Formal Arguments ..
  integer, dimension(:),    intent(in) :: nodes_of_cell_ptr
  integer, dimension(:),    intent(in) :: nodes_of_cell
  integer, dimension(:),    intent(in) :: old_noc
  integer, dimension(:,:),  intent(in) :: node_order
  real(wp), dimension(:,:), intent(in) :: xyz_nodes
  integer, dimension(:),    intent(in) :: rotation
  !
  !.. Local Scalars ..
  integer :: l,m,n,nc,ierr,iotec,n1,n2
  !
  !.. Local Arrays ..
  integer, dimension(8) :: ip
  integer, dimension(8) :: op
  integer, dimension(8) :: rev_rot
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "tecplot_cell_rotations"
  character(len=*), parameter :: fname = "cell_rotations.tec.dat"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  open (newunit=iotec,file=fname,status="replace",action="write", &
                      form="formatted",position="rewind", &
                      iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (iotec,1)
  !
  write (iotec,6) size(xyz_nodes,dim=2),size(rotation)
  !
  do n = 1,size(xyz_nodes,dim=2)
    write (iotec,4) (nint(xyz_nodes(l,n)), l=1,size(xyz_nodes,dim=1)), &
                    n,n,n
  end do
  !
  do nc = 1,size(nodes_of_cell_ptr)-1
    n1 = nodes_of_cell_ptr(nc)+1
    n2 = nodes_of_cell_ptr(nc+1)
    write (iotec,5) (old_noc(n), n=n1,n2)
  end do
  !
  do nc = 1,size(nodes_of_cell_ptr)-1
    !
    n1 = nodes_of_cell_ptr(nc)+1
    n2 = nodes_of_cell_ptr(nc+1)
    !
    ip = nodes_of_cell(n1:n2)
    op = old_noc(n1:n2)
    !
    write (iotec,2) nc,rotation(nc)
    !
    rev_rot(node_order(:,rotation(nc))+1) = intseq(1,8)
    !
    do n = 1,8
      !
     !write (iotec,4) (nint(xyz_nodes(l,op(n))), l=1,size(xyz_nodes,dim=1)), &
     !                ip(n),op(n),node_order(n,rotation(nc))+1
      write (iotec,4) (nint(xyz_nodes(l,ip(n))), l=1,size(xyz_nodes,dim=1)), &
                      n, ip(n), rev_rot(n)
      !
    end do
    !
   !write (iotec,5) (node_order(n,rotation(nc))+1, n=1,8)
    write (iotec,5) (n, n=1,8)
    !
    write (iotec,3) "r",nc,4*(nc-1)+2
    write (iotec,5) 1,2
   !write (iotec,5) rev_rot(1),rev_rot(2)
   !write (iotec,5) node_order(1,rotation(nc))+1,node_order(2,rotation(nc))+1
    write (iotec,3) "s",nc,4*(nc-1)+2
    write (iotec,5) 1,4
   !write (iotec,5) rev_rot(1),rev_rot(4)
   !write (iotec,5) node_order(1,rotation(nc))+1,node_order(4,rotation(nc))+1
    write (iotec,3) "t",nc,4*(nc-1)+2
    write (iotec,5) 1,5
   !write (iotec,5) rev_rot(1),rev_rot(5)
   !write (iotec,5) node_order(1,rotation(nc))+1,node_order(5,rotation(nc))+1
    !
  end do
  !
  !
  close (iotec,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format('VARIABLES = "X", "Y", "Z", "NoC", "Original NoC", "Local Index"')
  2 format('ZONE T="Cell ',i0,' - R(',i0,')", N=8, E=1, ', &
           'DATAPACKING=POINT, ZONETYPE=FEBRICK')
  3 format('ZONE T="',a,' - Cell ',i0,'", N=8, E=1, ', &
           'DATAPACKING=POINT, ZONETYPE=FELINESEG',/, &
           'VARSHARELIST=([1-6]=',i0,')')
  4 format(3(1x,i0),2(1x,i3),1x,i0)
  5 format(100(1x,i0))
 !5 format(*(1x,i0))
  6 format('ZONE T="Full Grid", N=',i0,', E=',i0,', ', &
           'DATAPACKING=POINT, ZONETYPE=FEBRICK')
  !
end subroutine tecplot_cell_rotations
!
!###############################################################################
!
subroutine check_node_ordering()
  !
  !.. Use Statements ..
  use geovar, only : nr,ncell,nfbnd,nodes_of_cell_ptr,cell_geom
  use geovar, only : xyz_nodes,nodes_of_cell,bface
  !
  !.. Local Scalars ..
  integer  :: i,n,nf,np,n1,n2,ip2,ip4
  real(wp) :: cell_vol
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: cn_xyz(:,:)
  real(wp), allocatable :: cell_nodes(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_node_ordering"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Lets do a quick test on the nodes defining each cell to
  ! make sure they are all in counter-clockwise order. If they
  ! arent, switch the ordering of the nodes so that they are.
  !
  ! First lets look at the cells on the grid boundaries since we
  ! also need to change bface(4) if the ordering is wrong
  !
  call debug_timer(start_timer,"checking node order of boundary faces")
  do nf = 1,nfbnd
    !
    n = bface(2,nf)
    !
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    np = n2-n1+1
    !
    ! Compute the signed area/volume of the current host cell
    !
    cell_nodes = xyz_nodes(1:nr,nodes_of_cell(n1:n2)) ! F2003 AUTO-REALLOC
    !
    select case (cell_geom(n))
      case (Geom_Edge)
        cell_vol = length_edge( cell_nodes(:,:) )
      case (Geom_Quad)
        cell_vol = area_quad( cell_nodes(:,:) )
      case (Geom_Tria)
        cell_vol = area_tria( cell_nodes(:,:) )
      case (Geom_Tetr)
        cell_vol = vol_tet( cell_nodes(:,:) )
      case (Geom_Pyra)
        cell_vol = vol_pyr( cell_nodes(:,:) )
      case (Geom_Pris)
        cell_vol = vol_pri( cell_nodes(:,:) )
      case (Geom_Hexa)
        cell_vol = vol_hex( cell_nodes(:,:) )
      case default
        write (iout,1)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
    end select
    !
    ! If the cell volume is negative, this means the cell is in clockwise
    ! ordering.  Switch two nodes (2 and 3 for triangles; 2 and 4 for quads)
    ! to change the cell to counter-clockwise ordering.
    !
    if (cell_vol < zero) then
      !
      if (any(cell_geom(n) == Geom_2D)) then
        !
        ip2 = nodes_of_cell(n1+1)
        ip4 = nodes_of_cell(n2  ) ! for triangles, ip4 = ip3 because n2 = n1+2
        !
        nodes_of_cell(n1+1) = ip4 ! for triangles, ip4 = ip3
        nodes_of_cell(n2  ) = ip2
        !
        ! We also need to adjust the value in bface(4,nf) depending on cell type
        ! For triangles, faces need to go from 1,2,3 to 3,2,1 ordering.
        !                This is done using:  bface(4) = 4 - bface(4)
        ! For quad, faces need to go from 1,2,3,4 to 4,3,2,1 ordering.
        !                This is done using:  bface(4) = 5 - bface(4)
        !
        bface(4,nf) = merge(4,5,cell_geom(n)==Geom_Tria) - bface(4,nf)
        !
        ! Now we have to loop through the remaining boundary faces
        ! and apply this correction to bface(4) for any other boundary
        ! faces that have this same host cell. Otherwise, when it gets
        ! to the other boundary face > nf with this same host cell, the
        ! node ordering will already be fixed and it will not make it into
        ! this cell_vol<zero conditional to apply a correction to bface(4).
        !
        do i = nf+1,nfbnd
          if (bface(2,i) == n) then
            bface(4,i) = merge(4,5,cell_geom(n)==Geom_Tria) - bface(4,i)
          end if
        end do
        !
      else
        !
        write (iout,2)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
        !
      end if
      !
    end if
    !
  end do
  call debug_timer(stop_timer,"checking node order of boundary faces")
  !
  ! Now lets do it over all the cells in the grid. This includes checking
  ! the cells on the grid boundaries a second time but they shouldnt be
  ! adjusted this time and the additional cost shouldnt hurt too much.
  !
  !
  ! Edge cells
  !
  call debug_timer(start_timer,"checking node order of edge cells")
  if (any(cell_geom(:) == Geom_Edge)) then
    cn_xyz = init_array(nr,geom_nodes(Geom_Edge),zero) ! F2003 AUTO-REALLOC
    do n = 1,ncell
      if (cell_geom(n) == Geom_Edge) then
        n1 = nodes_of_cell_ptr(n)+1
        n2 = nodes_of_cell_ptr(n+1)
        cn_xyz(:,:) = xyz_nodes(1:nr,nodes_of_cell(n1:n2))
        cell_vol = length_edge( cn_xyz(:,:) )
        if (cell_vol < zero) then
          write (iout,3) geom_name(Geom_Edge)
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
        end if
      end if
    end do
  end if
  call debug_timer(stop_timer,"checking node order of edge cells")
  !
  ! Quad cells
  !
  call debug_timer(start_timer,"checking node order of quad cells")
  if (any(cell_geom(:) == Geom_Quad)) then
    cn_xyz = init_array(nr,geom_nodes(Geom_Quad),zero) ! F2003 AUTO-REALLOC
    do n = 1,ncell
      if (cell_geom(n) == Geom_Quad) then
        n1 = nodes_of_cell_ptr(n)+1
        n2 = nodes_of_cell_ptr(n+1)
        cn_xyz(:,:) = xyz_nodes(1:nr,nodes_of_cell(n1:n2))
        cell_vol = area_quad( cn_xyz(:,:) )
        if (cell_vol < zero) then
          ip2 = nodes_of_cell(n1+1)
          ip4 = nodes_of_cell(n2  )
          nodes_of_cell(n1+1) = ip4
          nodes_of_cell(n2  ) = ip2
        end if
      end if
    end do
  end if
  call debug_timer(stop_timer,"checking node order of quad cells")
  !
  ! Tria cells
  !
  call debug_timer(start_timer,"checking node order of tria cells")
  if (any(cell_geom(:) == Geom_Tria)) then
    cn_xyz = init_array(nr,geom_nodes(Geom_Tria),zero) ! F2003 AUTO-REALLOC
    do n = 1,ncell
      if (cell_geom(n) == Geom_Tria) then
        n1 = nodes_of_cell_ptr(n)+1
        n2 = nodes_of_cell_ptr(n+1)
        cn_xyz(:,:) = xyz_nodes(1:nr,nodes_of_cell(n1:n2))
        cell_vol = area_tria( cn_xyz(:,:) )
        if (cell_vol < zero) then
          ip2 = nodes_of_cell(n1+1)
          ip4 = nodes_of_cell(n2  )
          nodes_of_cell(n1+1) = ip4
          nodes_of_cell(n2  ) = ip2
        end if
      end if
    end do
  end if
  call debug_timer(stop_timer,"checking node order of tria cells")
  !
  ! Tetr cells
  !
  call debug_timer(start_timer,"checking node order of tetr cells")
  if (any(cell_geom(:) == Geom_Tetr)) then
    cn_xyz = init_array(nr,geom_nodes(Geom_Tetr),zero) ! F2003 AUTO-REALLOC
    do n = 1,ncell
      if (cell_geom(n) == Geom_Tetr) then
        n1 = nodes_of_cell_ptr(n)+1
        n2 = nodes_of_cell_ptr(n+1)
        cn_xyz(:,:) = xyz_nodes(1:nr,nodes_of_cell(n1:n2))
        cell_vol = vol_tet( cn_xyz(:,:) )
        if (cell_vol < zero) then
          write (iout,3) geom_name(Geom_Tetr)
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
        end if
      end if
    end do
  end if
  call debug_timer(stop_timer,"checking node order of tetr cells")
  !
  ! Pyra cells
  !
  call debug_timer(start_timer,"checking node order of pyra cells")
  if (any(cell_geom(:) == Geom_Pyra)) then
    cn_xyz = init_array(nr,geom_nodes(Geom_Pyra),zero) ! F2003 AUTO-REALLOC
    do n = 1,ncell
      if (cell_geom(n) == Geom_Pyra) then
        n1 = nodes_of_cell_ptr(n)+1
        n2 = nodes_of_cell_ptr(n+1)
        cn_xyz(:,:) = xyz_nodes(1:nr,nodes_of_cell(n1:n2))
        cell_vol = vol_pyr( cn_xyz(:,:) )
        if (cell_vol < zero) then
          write (iout,3) geom_name(Geom_Pyra)
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
        end if
      end if
    end do
  end if
  call debug_timer(stop_timer,"checking node order of pyra cells")
  !
  ! Pris cells
  !
  call debug_timer(start_timer,"checking node order of pris cells")
  if (any(cell_geom(:) == Geom_Pris)) then
    cn_xyz = init_array(nr,geom_nodes(Geom_Pris),zero) ! F2003 AUTO-REALLOC
    do n = 1,ncell
      if (cell_geom(n) == Geom_Pris) then
        n1 = nodes_of_cell_ptr(n)+1
        n2 = nodes_of_cell_ptr(n+1)
        cn_xyz(:,:) = xyz_nodes(1:nr,nodes_of_cell(n1:n2))
        cell_vol = vol_pri( cn_xyz(:,:) )
        if (cell_vol < zero) then
          write (iout,3) geom_name(Geom_Pris)
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
        end if
      end if
    end do
  end if
  call debug_timer(stop_timer,"checking node order of pris cells")
  !
  ! Hexa cells
  !
  call debug_timer(start_timer,"checking node order of hexa cells")
  if (any(cell_geom(:) == Geom_Hexa)) then
    cn_xyz = init_array(nr,geom_nodes(Geom_Hexa),zero) ! F2003 AUTO-REALLOC
    do n = 1,ncell
      if (cell_geom(n) == Geom_Hexa) then
        n1 = nodes_of_cell_ptr(n)+1
        n2 = nodes_of_cell_ptr(n+1)
        cn_xyz(:,:) = xyz_nodes(1:nr,nodes_of_cell(n1:n2))
        cell_vol = vol_hex( cn_xyz(:,:) )
        if (cell_vol < zero) then
          write (iout,3) geom_name(Geom_Hexa)
          call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
        end if
      end if
    end do
  end if
  call debug_timer(stop_timer,"checking node order of hexa cells")
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/, &
    " ############################## ERROR! ##############################",/, &
    "     INVALID CELL GEOMETRY WHEN CHECKING THE CELL NODE ORDERINGS!",/, &
    " ############################## ERROR! ##############################",/)
  2 format (/, &
    " ############################## ERROR! ##############################",/, &
    "   NEGATIVE VOLUME DETECTED WHEN CHECKING THE CELL NODE ORDERINGS!",/, &
    "         UNFORTUNATELY, THE NODES CANNOT BE REARRANGED",/, &
    "         TO FIX THIS FOR 3D SIMULATIONS!",/, &
    " ############################## ERROR! ##############################",/)
  3 format (/, &
    " ############################## ERROR! ##############################",/, &
    "   NEGATIVE VOLUME DETECTED WHEN CHECKING THE ",a, &
                                     " CELL NODE ORDERINGS!",/, &
    "         UNFORTUNATELY, THE NODES CANNOT BE REARRANGED",/, &
    "         TO FIX THIS FOR 3D SIMULATIONS!",/, &
    " ############################## ERROR! ##############################",/)
  !
end subroutine check_node_ordering
!
!###############################################################################
!
subroutine adjust_for_ghost_cells()
  !
  ! Routine to make some adjustments to the data from the grid file
  ! to account for ghost cells and new nodes associated with them
  !
  !.. Use Statements ..
  use geovar, only : nr,ncell,nfbnd,nnode,n_totcel,n_totnod
  use geovar, only : bface,nodes_of_cell_ptr,nodes_of_cell,xyz_nodes
  use geovar, only : cell_geom,cell_order
  use ovar,   only : Grid_Min,Grid_Max,Grid_Length
  !
  !.. Local Scalars ..
  integer :: g1,g2,h1,nf,nfnods,switch_index
  integer :: ghost_cell,host_cell,host_geom,host_side
  !
  !.. Local Arrays ..
  integer, dimension(1:4) :: ip
  integer, dimension(1:4) :: face_nodes
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "adjust_for_ghost_cells"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Before we start adding ghost cells, lets get the extents of the grid domain
  ! NOTE: ASSUME GRID BOUNDARIES OF CONSTANT X AND Y
  !
  Grid_Min = zero
  Grid_Max = zero
  !
  Grid_Min(1:nr) = minval(xyz_nodes(1:nr,:),dim=2)
  Grid_Max(1:nr) = maxval(xyz_nodes(1:nr,:),dim=2)
  Grid_Length = Grid_Max - Grid_Min
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Total number of interior plus ghost cells
  !
  n_totcel = ncell + nfbnd
  !
  ! Reallocate nodes_of_cell_ptr to account for ghost cells
  !
  call reallocate( nodes_of_cell_ptr , n_totcel+1 )
  !
  ! Set the number of nodes for each ghost cell.
  ! NOTE: The geometry of the ghost cells is one dimension less than the host
  !       cells. For example, the ghost cell geometry will be an edge for a 2D
  !       grid, and it will be either a triangle or a quadrilateral for a 3D
  !       grid.
  !
  do nf = 1,nfbnd
    !
    host_geom = bface(3,nf)
    host_side = bface(4,nf)
    !
    ! Make sure the host geometry is valid
    !
    if (all(host_geom /= [Geom_1D,Geom_2D,Geom_3D])) then
      write (error_message,100)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
    !
    nodes_of_cell_ptr(ncell+nf+1) = num_face_nodes(host_side,host_geom)
    !
  end do
  !
  ! Reshuffle the ghost cell portion of nodes_of_cell_ptr
  !
  do nf = ncell+2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(nf) = nodes_of_cell_ptr(nf) + nodes_of_cell_ptr(nf-1)
  end do
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Reallocate nodes_of_cell to account for ghost cells
  !
  call reallocate( nodes_of_cell , last(nodes_of_cell_ptr) )
  !
  ! Reallocate cell_geom to account for ghost cells
  !
  call reallocate( cell_geom , n_totcel )
  !
  ! Reallocate cell_order to account for ghost cells
  !
  call reallocate( cell_order , n_totcel )
  !
  ! Loop over the boundary faces and find the nodes of each ghost cell.
  ! The nodes are in counter-clockwise ordering. Also, while we have the
  ! information available, fill in the ghost cell data for the arrays
  ! cell_geom and cell_order.
  !
  do nf = 1,nfbnd
    !
    ghost_cell = nf + ncell
    !
    ! Node indices for the ghost cell
    !
    g1 = nodes_of_cell_ptr(ghost_cell)+1
    g2 = nodes_of_cell_ptr(ghost_cell+1)
    !
    ! Information about boundary face
    !
    host_cell = bface(2,nf)
    host_geom = bface(3,nf)
    host_side = bface(4,nf)
    !
    ! Starting node index for the interior cell
    !
    h1 = nodes_of_cell_ptr(host_cell)+1
    !
    ! Get the number of nodes defining this face
    !
    nfnods = num_face_nodes(host_side,host_geom)
    !
    ! Check that the face geometry matches the number of ghost nodes
    !
    if (g2-g1+1 /= nfnods) then
      write (error_message,102)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
    !
    ! Node offsets for this face of the host cell
    !
    ip(1:nfnods) = side_nodes(1:nfnods,host_side,host_geom)
    !
    ! Use the node offsets to get the actual indices for these nodes
    !
    face_nodes(1:nfnods) = nodes_of_cell( h1 + ip(1:nfnods) )
    !
    ! Assign nodes defining boundary face from the host cell to the ghost cell
    ! NOTE: We have to reverse the rotation/ordering of the face nodes so that
    !       the nodes of the ghost cell are defined such that the normal of
    !       each of its faces points outwards (this is the same as the interior
    !       cells). For a 2D host cell, there are only two nodes so we just
    !       need to switch them, i.e. node1 <=> node2.  For a 3D host cell,
    !       switching either nodes 1 & 3 or nodes 3 & 4 will do the trick.
    !       Since triangle faces only have three nodes, we will switch nodes
    !       1 & 3 for both triangle and quad faces so that we dont have to
    !       distinguish between the two.
    !
    ! Copy over the face nodes to the ghost cell,
    ! keeping the same ordering as the host cell.
    !
    nodes_of_cell(g1:g2) = face_nodes(1:nfnods)
    !
    ! Get the index of the node to switch with node 1.
    ! NOTE: switch_index = nfnods (i.e. 2) for a 2D host cell
    !       switch_index = nfnods (i.e. 3) for a triangular boundary face
    !       switch_index = 3               for a quadrilateral boundary face
    ! ADDITIONAL NOTE: switch_index also evaluates to nfnods (i.e. 1) for
    !                  a 1D host cell which is the correct value if the
    !                  rest of the code is ever set up to run 1D.
    !
    switch_index = min( nfnods , 3 )
    !
    ! Switch the two face nodes required to reverse the face node ordering
    !
    nodes_of_cell(g1) = face_nodes(switch_index)
    nodes_of_cell(g1-1+switch_index) = face_nodes(1)
    !
    ! Get the geometry of the ghost cell (i.e. face geometry) and set the order
    ! of the ghost cell to the same as the host cell.
    !
    cell_geom(ghost_cell) = geom_of_face(nfnods)
    cell_order(ghost_cell) = cell_order(host_cell)
    !
  end do
  !
  n_totnod = nnode
  !
  ! If using multiple processors, save the global grid arrays
  !
  call save_global_grid
  !
  call debug_timer(leaving_procedure,pname)
  !
  1 format ("HOSTCELL = ",i0,"  PCELL = ",i0,"   (delta x,delta y) = ",2es24.11)
  2 format (/," ********************* WARNING!!! *********************",/, &
              "  IT APPEARS THE FINAL TIME DOES NOT CORRESPOND TO THE",/, &
              "  EXACT END OF A PERIOD FOR THE VORTEX PROBLEM. THE",/, &
              "  FINAL TIME CORRESPONDS TO",f12.5," TOTAL PERIODS.",/, &
              " ********************* WARNING!!! *********************",/)
  100 format ("Not sure what nodes to use for the ghost cell based on ", &
              "the host cell geometry encountered.")
  101 format ("The face geometry across a periodic boundary does not match!")
  102 format ("The number of ghost nodes does not match this face geometry!")
  !
end subroutine adjust_for_ghost_cells
!
!###############################################################################
!
subroutine save_global_grid
  !
  !.. Use Statements ..
  use geovar,    only : nr,nnode,ncell,nfbnd,n_totnod,n_totcel
  use geovar,    only : n_global_node,n_global_cell,n_global_fbnd
  use geovar,    only : n_global_totnod,n_global_totcel
  use geovar,    only : n_global_solpts,n_global_totpts
  use geovar,    only : xyz_nodes,global_xyz_nodes
  use geovar,    only : nodes_of_cell_ptr,global_nofc_ptr
  use geovar,    only : nodes_of_cell,global_nofc
  use geovar,    only : global_pinc_ptr
  use geovar,    only : cell_geom,global_cell_geom
  use geovar,    only : cell_order,global_cell_order
  use order_mod, only : geom_solpts
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,this_geom,this_order,ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "save_global_grid"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  n_global_node = nnode
  n_global_cell = ncell
  n_global_fbnd = nfbnd
  n_global_totnod = n_totnod
  n_global_totcel = n_totcel
  !
 !n1 = size(xyz_nodes,dim=1)
 !n2 = size(xyz_nodes,dim=2)
 !allocate ( global_xyz_nodes(1:n1,1:n2) , source=xyz_nodes , &
 !           stat=ierr , errmsg=error_message )
 !call alloc_error(pname,"global_xyz_nodes",1,__LINE__,__FILE__,ierr, &
 !                 error_message)
 !!
 !n1 = size(nodes_of_cell_ptr)
 !allocate ( global_nofc_ptr(1:n1) , source=nodes_of_cell_ptr , &
 !           stat=ierr , errmsg=error_message )
 !call alloc_error(pname,"global_nofc_ptr",1,__LINE__,__FILE__,ierr, &
 !                 error_message)
 !!
 !n1 = size(nodes_of_cell)
 !allocate ( global_nofc(1:n1) , source=nodes_of_cell , &
 !           stat=ierr , errmsg=error_message )
 !call alloc_error(pname,"global_nofc",1,__LINE__,__FILE__,ierr,error_message)
  !
  n1 = size(cell_geom)
  allocate ( global_cell_geom(1:n1) , source=cell_geom , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"global_cell_geom",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  n1 = size(cell_order)
  allocate ( global_cell_order(1:n1) , source=cell_order , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"global_cell_order",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! Create the global integer pointer for the beginning and ending
  ! indices for the solution points in each cell. This will be done
  ! later but a global version is never made before the global grid
  ! is partitioned and each processor creates its own separate
  ! local grid.
  !
  allocate ( global_pinc_ptr(1:n_global_totcel+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"global_pinc_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  do n = 1,n_global_totcel
    !
    this_geom = global_cell_geom(n)
    this_order = global_cell_order(n)
    !
    global_pinc_ptr(n+1) = global_pinc_ptr(n) &
                         + geom_solpts(this_geom,this_order)
    !
  end do
  !
  n_global_solpts = global_pinc_ptr(ncell+1)
  n_global_totpts = last(global_pinc_ptr)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine save_global_grid
!
!###############################################################################
!
subroutine find_bc_conflicts
  !
  !.. Use Statements ..
  use geovar, only : nfbnd,bface,nnode
  use order_mod, only : n_fp1de
  use geovar, only : nodes_of_cell_ptr,nodes_of_cell
  use geovar, only : bc_conflict
  use geovar, only : cell_geom
  !
  !.. Local Scalars ..
  integer :: n,nf,i,ip,is,j,n1,n2,nf1,nf2,ierr
  integer :: host_cell,host_geom,host_side
  integer :: n_2bc_nodes,nfm,nfnods
  character(len=100) :: array_name
  !
  !.. Local Arrays ..
  integer, dimension(1:4) :: node_offsets
  integer, dimension(Geom_Min:Geom_Max) :: error_count
  !
  !.. Local Allocatable Arrays ..
  logical(lk), allocatable, dimension(:) :: mskp
  integer, allocatable, dimension(:) :: faces_with_node_ptr
  integer, allocatable, dimension(:) :: faces_with_node
  integer, allocatable, dimension(:,:) :: bface_nodes
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_bc_conflicts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  error_count = 0
  !
  allocate ( bface_nodes(1:4,1:nfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface_nodes",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Find the nodes that define each boundary face
  !
  do nf = 1,nfbnd
    !
    host_cell = bface(2,nf) ! host cell of boundary face
    host_geom = bface(3,nf) ! geometry type of host cell
    host_side = bface(4,nf) ! face of host cell on boundary
    !
    n1 = nodes_of_cell_ptr(host_cell)+1 ! beginning index for host cell nodes
    n2 = nodes_of_cell_ptr(host_cell+1) ! ending    index for host cell nodes
    !
    ! Get the number of nodes defining this face
    !
    nfnods = num_face_nodes(host_side,host_geom)
    !
    ! Node offsets for this face of the host cell
    !
    node_offsets(1:nfnods) = side_nodes(1:nfnods,host_side,host_geom)
    !
    ! Use the node offsets to get the actual indices for these nodes
    !
    bface_nodes(1:nfnods,nf) = nodes_of_cell( n1 + node_offsets(1:nfnods) )
    !
  end do
  !
  ! Check for any bad nodes
  !
  if (any(bface_nodes < 0)) then
    !
    error_count(Geom_Node) = count(bface_nodes(1,:) == -Geom_Node)
    error_count(Geom_Edge) = count(bface_nodes(1,:) == -Geom_Edge)
    error_count(Geom_Tria) = count(bface_nodes(1,:) == -Geom_Tria)
    error_count(Geom_Quad) = count(bface_nodes(1,:) == -Geom_Quad)
    error_count(Geom_Tetr) = count(bface_nodes(1,:) == -Geom_Tetr)
    error_count(Geom_Pyra) = count(bface_nodes(1,:) == -Geom_Pyra)
    error_count(Geom_Pris) = count(bface_nodes(1,:) == -Geom_Pris)
    error_count(Geom_Hexa) = count(bface_nodes(1,:) == -Geom_Hexa)
    error_count(Geom_Unknown) = count(bface_nodes(1,:) == -Geom_Unknown)
    !
    write (error_message,4)
    do n = Geom_Min,Geom_Max
      if (error_count(n) > 0) then
        write (error_message,5) trim(error_message)//new_line('a'), &
                                error_count(n),trim(Geom_String(n))
      end if
    end do
    write (error_message,6) trim(error_message)//new_line('a'),sum(error_count)
    !
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    !
  end if
  !
  ! Mark the nodes defining each boundary face, storing ahead
  ! NOTE: Skip all communication type boundary faces
  !
  allocate ( faces_with_node_ptr(1:nnode+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faces_with_node_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  do nf = 1,nfbnd
    if (any(bface(1,nf) == bc_comm)) cycle ! skip all communication boundaries
    do i = 1,count(bface_nodes(:,nf) > 0)
      ip = bface_nodes(i,nf) + 1
      faces_with_node_ptr(ip) = faces_with_node_ptr(ip) + 1
    end do
  end do
  !
  ! Reshuffle faces_with_node_ptr
  !
  do n = 2,size(faces_with_node_ptr)
    faces_with_node_ptr(n) = faces_with_node_ptr(n) + faces_with_node_ptr(n-1)
  end do
  !
  allocate ( faces_with_node(1:last(faces_with_node_ptr)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faces_with_node",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! Now store the faces surrounding each node in faces_with_node
  ! NOTE: Skip all communication type boundary faces
  !
  do nf = 1,nfbnd
    if (any(bface(1,nf) == bc_comm)) cycle ! skip all communication boundaries
    do i = 1,count(bface_nodes(:,nf) > 0)
      ip = bface_nodes(i,nf)
      is = faces_with_node_ptr(ip) + 1
      faces_with_node_ptr(ip) = is
      faces_with_node(is) = nf
    end do
  end do
  !
  ! Finally, reorder faces_with_node_ptr
  !
  do n = size(faces_with_node_ptr),2,-1
    faces_with_node_ptr(n) = faces_with_node_ptr(n-1)
  end do
  faces_with_node_ptr(1) = 0
  !
  ! Now go through each node and see if any of the boundary faces
  ! containing the node have non-matching boundary conditions.
  ! Mark any nodes that are found.
  !
  allocate ( mskp(1:nnode) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"mskp",1,__LINE__,__FILE__,ierr,error_message)
  !
  do n = 1,size(faces_with_node_ptr)-1
    nf1 = faces_with_node_ptr(n)+1
    nf2 = faces_with_node_ptr(n+1)
   !write (iout,'(1x,i0,":",i0,3x,l1,3x,*(1x,i0,:))') &
   !         nf1,nf2,all_not_equal( bface(1,faces_with_node(nf1:nf2)) ), &
   !         (faces_with_node(i),i=nf1,nf2)
    mskp(n) = all_not_equal( bface(1,faces_with_node(nf1:nf2)) )
  end do
  !
  allocate ( bc_conflict(1:count(mskp)) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bc_conflict",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Use the mask and faces_with_node arrays to extract the nodes
  ! and corresponding boundary faces that have non-matching
  ! boundary conditions.
  !
  i = 0
  do n = 1,size(faces_with_node_ptr)-1
    !
    if (.not. mskp(n)) cycle
    !
    i = i + 1
    !
    nf1 = faces_with_node_ptr(n)+1
    nf2 = faces_with_node_ptr(n+1)
    !
    ! Find the "master" boundary face with this node, i.e., the face
    ! that will contain the flux point that overrides the other flux
    ! point solutions at this node.
    !
    nfm = bc_priority(bface(1,faces_with_node(nf1:nf2)))
    if (nfm == 0) then
      write (iout,*) "Error: Function bc_priority returned value 0!"
    end if
    !
    allocate ( bc_conflict(i)%slave(1:nf2-nf1) , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "bc_conflict(",i,")%slave"
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Save the master boundary face information and set the
    ! partition number to the root partition number (zero).
    !
    bc_conflict(i)%master%face = faces_with_node(nf1-1+nfm)
    bc_conflict(i)%master%cpu = 0
    !
    ! Automatically mark that this partition is part of the communication
    ! for this conflict and is only local to this processor. This is the
    ! required default for the serial case and will be altered as necessary
    ! if running in parallel.
    !
    bc_conflict(i)%not_my_conflict = fals
    bc_conflict(i)%is_local = true
    !
    ! Store all faces for the "slave" boundary faces.
    ! NOTE: When nf >= nfm, omitting the "-1" in the index to faces_with_node
    !       is necessary to skip the index "nf1-1+nfm" within the array.
    !
    do nf = 1,nf2-nf1
      if (nf < nfm) then
        bc_conflict(i)%slave(nf)%face = faces_with_node(nf1-1+nf)
      else if (nf >= nfm) then
        bc_conflict(i)%slave(nf)%face = faces_with_node(nf1  +nf)
      end if
      bc_conflict(i)%slave(nf)%cpu = 0
    end do
    !
    ! Using the node for this conflict, determine which flux point
    ! of the boundary face will need to be exchanged.
    ! NOTE: This will need to be generalized for 3D grids where
    !       the conflicts will be along a cell edge and not just
    !       at a single node.
    !
    ! Start with the master face
    !
    nf = bc_conflict(i)%master%face
    !
    ! ##########################################################################
    ! ##########################################################################
    ! ###################### !!!! START OF TEMPORARY !!!! ######################
    ! ##########################################################################
    ! ##########################################################################
    !
    ! Quick check to see if the host cell geometry is not 2D. The following
    ! conditional statements need to be altered before they are valid
    ! for non-2D cells so kill the program if the host cell is not 2D.
    !
    if (all(cell_geom(bface(2,nf)) /= Geom_2D)) then
      write (error_message,7)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
    !
    ! ##########################################################################
    ! ##########################################################################
    ! ####################### !!!! END OF TEMPORARY !!!! #######################
    ! ##########################################################################
    ! ##########################################################################
    !
    if (n == bface_nodes(1,nf)) then
      bc_conflict(i)%master%flux_point = 1
    else if (n == bface_nodes(2,nf)) then
      bc_conflict(i)%master%flux_point = n_fp1de
    else
      write (iout,8)
    end if
    !
    ! Now do the same for the slave faces
    !
    do j = 1,size(bc_conflict(i)%slave)
      !
      nf = bc_conflict(i)%slave(j)%face
      !
      ! ########################################################################
      ! ########################################################################
      ! ##################### !!!! START OF TEMPORARY !!!! #####################
      ! ########################################################################
      ! ########################################################################
      !
      ! Quick check to see if the host cell geometry is not 2D. The following
      ! conditional statements need to be altered before they are valid
      ! for non-2D cells so kill the program if the host cell is not 2D.
      !
      if (all(cell_geom(bface(2,nf)) /= Geom_2D)) then
        write (error_message,7)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      !
      ! ########################################################################
      ! ########################################################################
      ! ###################### !!!! END OF TEMPORARY !!!! ######################
      ! ########################################################################
      ! ########################################################################
      !
      if (n == bface_nodes(1,nf)) then
        bc_conflict(i)%slave(j)%flux_point = 1
      else if (n == bface_nodes(2,nf)) then
        bc_conflict(i)%slave(j)%flux_point = n_fp1de
      else
        write (iout,8)
      end if
      !
    end do
    !
  end do
  !
  ! Deallocate the local arrays before leaving
  !
  deallocate ( bface_nodes , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface_nodes",2,__LINE__,__FILE__,ierr,error_message)
  deallocate ( faces_with_node_ptr , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faces_with_node_ptr",2,__LINE__,__FILE__,ierr, &
                   error_message)
  deallocate ( faces_with_node , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faces_with_node",2,__LINE__,__FILE__,ierr, &
                   error_message)
  deallocate ( mskp , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"mskp",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (a,i0,a)
  2 format (//," Number of BC Conflicts = ",i0,/)
  3 format (" Conflict #",i0,": Node = ",i0,"; Faces = ",100(i0,:,", "))
 !3 format (" Conflict #",i0,": Node = ",i0,"; Faces = ",*(i0,:,", "))
  4 format ("ERROR FINDING THE NODES OF")
  5 format (a,20x,i10,' - ',a)
  6 format (a,23x,"HOST CELLS ON ",i0," BOUNDARY FACES!")
  7 format ("The following conditionals need updating for non-2D Geometries.")
  8 format (" ERROR: Conflict node doesn't seem to match a boundary node!")
  !
end subroutine find_bc_conflicts
!
!###############################################################################
!
subroutine check_face_flux_point_orderings()
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  use geovar, only : ncell,cell
  use geovar, only : nfbnd,nface,face,fp
  use geovar, only : nr,xyz,bface,bc_grps
  use interpolation_mod, only : interp
  use flowvar, only : face_var,faceusp
  use flowvar, only : facexyz
  use parallel_mod, only : exchange_faceusp,wait_faceusp
  use parallel_mod, only : exchange_faceusp_dt,wait_faceusp_dt
  use geovar, only : xyz_nodes
  use ovar, only : VortexGrid_Lref
  use geovar, only : nodes_of_cell,nodes_of_cell_ptr
  !
  !.. Local Scalars ..
  integer :: i,j,nf,nfp,kf,kl,kr,l,ierr
  integer :: l1,l2,r1,r2,npl,npr,pf,f0,ntfp
  integer :: n1,n2,n0,np,bc_type
  integer :: face_geom,face_order
  integer :: rotation,fp_offset
  integer :: left_cell,left_geom,left_order
  integer :: right_cell,right_geom,right_order
  integer :: l1file,r1file,d1file
  integer :: l2file,r2file,d2file
  real(wp) :: a,val
  real(wp) :: ds1,ds2,ds3,ds4,ds5,ds6
  character(len=25) :: face_string
  character(len=100) :: array_name
  integer :: il1,il2,ir1,ir2
  !
  !.. Local Arrays ..
  type(face_var), allocatable :: face_xyz(:)
  integer, allocatable :: left_nodes(:)
  integer, allocatable :: right_nodes(:)
  integer, allocatable :: node_order(:)
  logical(lk), allocatable :: msk(:)
  real(wp) :: left_xyz(1:nr)
  real(wp) :: right_xyz(1:nr)
  !
  integer, allocatable :: abc_left(:)
  integer, allocatable :: abc_right(:)
  integer, allocatable :: new_right_nodes(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_face_flux_point_orderings"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  call output_element_solpts_and_flxpts
  !
  call check_periodic_coordinates
  !
  do i = 1,ncell
    n1 = cell(i)%beg_sp
    n2 = cell(i)%end_sp
    write (102,55) i
    do j = n1,n2
      write (102,56) (xyz(l,j), l=1,nr)
    end do
  end do
  55 format ("Cell # ",i0)
  56 format (3es16.7)
  !
  face_geom = Geom_Tria
  do j = 0,1
    do face_order = lbound(fp%geord,dim=2),ubound(fp%geord,dim=2)
      write (103+j,553) face_order
      do i = lbound(fp%geord(face_geom,face_order)%rot,dim=1), &
             ubound(fp%geord(face_geom,face_order)%rot,dim=1)
        write (103+j,554) i
        do kf = lbound(fp%geord(face_geom,face_order)%rot(i)%pts,dim=1), &
                ubound(fp%geord(face_geom,face_order)%rot(i)%pts,dim=1)
          write (103+j,555) kf,fp%geord(face_geom,face_order)%rot(i)%pts(kf)+1
        end do
      end do
    end do
    face_geom = Geom_Quad
  end do
  553 format (/,40("#"),/,"Face order = ",i0)
  554 format (2x,"Rotation # ",i0)
  555 format (4x,"kf = ",i3,",  pts(kf) = ",i3)
  !
  l1file = 1101 + mypnum
  r1file = 1201 + mypnum
  d1file = 1301 + mypnum
  l2file = 2101 + mypnum
  r2file = 2201 + mypnum
  d2file = 2301 + mypnum
  !
  write (l1file,*)
  write (r1file,*)
  write (d1file,*)
  write (l2file,*)
  write (r2file,*)
  write (d2file,*)
  !
  allocate ( face_xyz(1:nface) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_xyz",1,__LINE__,__FILE__,ierr,error_message)
  !
  do nf = 1,nface
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    nfp = geom_solpts(face_geom,face_order)
    allocate ( face_xyz(nf)%v(1:nq,1:nfp,1:2) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,'("face_xyz(",i0,")%v")') nf
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  end do
  !
  ntfp = face(nface)%end_fp
  allocate ( faceusp(1:ntfp,1:nq,1:2) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faceusp",1,__LINE__,__FILE__,ierr,error_message)
  !
  do nf = 1,nfbnd
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    left_cell = face(nf)%left%cell
    !
    left_geom = cell(left_cell)%geom
    left_order = cell(left_cell)%order
    !
    l1 = cell(left_cell)%beg_sp
    l2 = cell(left_cell)%end_sp
    npl = l2-l1+1
    !
    nfp = geom_solpts(face_geom,face_order)
    !
    do kf = 1,nfp
      !
      kl = face(nf)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                        rot(face(nf)%left%rotation)% &
                                     pts(kf)
      !
      ! Interpolate the coordinates for this flux point for the left cell
      !
      do l = 1,nr
        val = zero
        do i = interp(left_geom,left_order)%toFace(face_order)%row_ptr(kl)+1, &
               interp(left_geom,left_order)%toFace(face_order)%row_ptr(kl+1)
          !
          j = interp(left_geom,left_order)%toFace(face_order)%row_idx(i)
          a = interp(left_geom,left_order)%toFace(face_order)%row_val(i)
          !
          val = val + a * xyz(l,l1-1+j)
          !
        end do
        face_xyz(nf)%v(l,kf,1) = val
      end do
      !
    end do
    !
    if (bface(1,nf) == bc_periodic) then
      !
      pf = bface(5,nf)
      !
      right_cell = face(pf)%left%cell
      !
      right_geom = cell(right_cell)%geom
      right_order = cell(right_cell)%order
      !
      r1 = cell(right_cell)%beg_sp
      r2 = cell(right_cell)%end_sp
      npr = r2-r1+1
      !
      rotation = face(nf)%right%rotation
      fp_offset = face(pf)%left%fp_offset
      !
      do kf = 1,nfp
        !
        kr = fp_offset + fp%geord(face_geom,face_order)% &
                                          rot(rotation)% &
                         pts(kf)
        !
        ! Interpolate the coordinates for this flux point for the right cell
        !
        do l = 1,nr
          val = zero
          do i = interp(right_geom,right_order)% &
                             toFace(face_order)%row_ptr(kr)+1, &
                 interp(right_geom,right_order)% &
                             toFace(face_order)%row_ptr(kr+1)
            !
            j = interp(right_geom,right_order)%toFace(face_order)%row_idx(i)
            a = interp(right_geom,right_order)%toFace(face_order)%row_val(i)
            !
            val = val + a * xyz(l,r1-1+j)
            !
          end do
          face_xyz(nf)%v(l,kf,2) = val
        end do
        !
      end do
      !
    else if (bface(1,nf) /= bc_cpu_bnd) then
      !
      do kf = 1,nfp
        do l = 1,nr
          face_xyz(nf)%v(l,kf,2) = face_xyz(nf)%v(l,kf,1)
        end do
      end do
      !
    end if
    !
  end do
  !
  do nf = nfbnd+1,nface
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    left_cell = face(nf)%left%cell
    right_cell = face(nf)%right%cell
    !
    left_geom = cell(left_cell)%geom
    left_order = cell(left_cell)%order
    !
    right_geom = cell(right_cell)%geom
    right_order = cell(right_cell)%order
    !
    l1 = cell(left_cell)%beg_sp
    l2 = cell(left_cell)%end_sp
    npl = l2-l1+1
    !
    r1 = cell(right_cell)%beg_sp
    r2 = cell(right_cell)%end_sp
    npr = r2-r1+1
    !
    nfp = geom_solpts(face_geom,face_order)
    !
    do kf = 1,nfp
      !
      kl = face(nf)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                        rot(face(nf)%left%rotation)% &
                                     pts(kf)
      kr = face(nf)%right%fp_offset + fp%geord(face_geom,face_order)% &
                                        rot(face(nf)%right%rotation)% &
                                      pts(kf)
      !
      ! Interpolate the coordinates for this flux point for the left cell
      !
      do l = 1,nr
        val = zero
        do i = interp(left_geom,left_order)%toFace(face_order)%row_ptr(kl)+1, &
               interp(left_geom,left_order)%toFace(face_order)%row_ptr(kl+1)
          !
          j = interp(left_geom,left_order)%toFace(face_order)%row_idx(i)
          a = interp(left_geom,left_order)%toFace(face_order)%row_val(i)
          !
          val = val + a * xyz(l,l1-1+j)
          !
        end do
        face_xyz(nf)%v(l,kf,1) = val
      end do
      !
      ! Interpolate the coordinates for this flux point for the right cell
      !
      do l = 1,nr
        val = zero
        do i = interp(right_geom,right_order)%toFace(face_order)%row_ptr(kr)+1,&
               interp(right_geom,right_order)%toFace(face_order)%row_ptr(kr+1)
          !
          j = interp(right_geom,right_order)%toFace(face_order)%row_idx(i)
          a = interp(right_geom,right_order)%toFace(face_order)%row_val(i)
          !
          val = val + a * xyz(l,r1-1+j)
          !
        end do
        face_xyz(nf)%v(l,kf,2) = val
      end do
      !
    end do
    !
  end do
  !
  faceusp(:,:,:) = zero
  !
  do nf = 1,nface
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
   !do i = 1,2
    do kf = 1,nfp
      do l = 1,nr
   !    write (iout,100) mypnum,nf,f0,kf,l,i
   !    100 format ("CPU # ",i0," : ",i5,i5,i3,i2,i2)
   !    flush (iout)
        faceusp(f0+kf,l,1) = face_xyz(nf)%v(l,kf,1)
      end do
    end do
   !end do
   !!
    if (nf > nfbnd) then
      do kf = 1,nfp
        do l = 1,nr
          faceusp(f0+kf,l,2) = face_xyz(nf)%v(l,kf,2)
        end do
      end do
    end if
    !
  end do
  !
  do i = 1,size(bc_grps)
    !
    n1 = bc_grps(i)%beg_fp
    n2 = bc_grps(i)%end_fp
    !
    n0 = bc_grps(i)%beg_fp - 1
    np = bc_grps(i)%end_fp - n0
    !
    bc_type = bc_grps(i)%bc_type
    !
    if (bc_type == bc_periodic) then
      !
      call transfer_periodic(n1,n2,faceusp,bc_grps(i)%periodic_list)
      !
    else
      !
      do kf = 1,np
        do l = 1,nr
          faceusp(n0+kf,l,2) = faceusp(n0+kf,l,1)
        end do
      end do
      !
    end if
    !
  end do
  !
  if (ncpu > 1) then
    !
    call exchange_faceusp_dt(face_xyz)
    call wait_faceusp_dt(face_xyz)
    !
    call exchange_faceusp
    call wait_faceusp
    !
  end if
  !
  do nf = 1,nface
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    rotation = 0
    !
    face_string = ""
    if (nf <= nfbnd) then
      if (bface(1,nf) == bc_periodic) then
        pf = bface(5,nf)
        face_string = "(periodic face)"
        rotation = bface(11,nf)/10000
        write (mypnum+4201,924) &
          nf-nfbnd,nf,face(nf)%left%cell_face,face(pf)%left%cell_face, &
          rotation
  924 format (i6,1x,i5,6x,i1,i2,2x,i0)
      else if (bface(1,nf) == bc_cpu_bnd) then
        face_string = "(cpu boundary)"
      else
        face_string = "(non-communication BC)"
      end if
    end if
    !
    if (rotation /= 0) then
      write (l1file,12) nf-nfbnd,trim(adjustl(face_string)),rotation
      write (r1file,12) nf-nfbnd,trim(adjustl(face_string)),rotation
      write (d1file,12) nf-nfbnd,trim(adjustl(face_string)),rotation
      write (l2file,12) nf-nfbnd,trim(adjustl(face_string)),rotation
      write (r2file,12) nf-nfbnd,trim(adjustl(face_string)),rotation
      write (d2file,12) nf-nfbnd,trim(adjustl(face_string)),rotation
    else
      write (l1file,12) nf-nfbnd,trim(adjustl(face_string))
      write (r1file,12) nf-nfbnd,trim(adjustl(face_string))
      write (d1file,12) nf-nfbnd,trim(adjustl(face_string))
      write (l2file,12) nf-nfbnd,trim(adjustl(face_string))
      write (r2file,12) nf-nfbnd,trim(adjustl(face_string))
      write (d2file,12) nf-nfbnd,trim(adjustl(face_string))
    end if
    !
    do kf = 1,nfp
      !
      write (l1file,11) (merge(zero,face_xyz(nf)%v(l,kf,1), &
                                abs(face_xyz(nf)%v(l,kf,1))<=eps14),l=1,nr)
      write (l2file,11) (merge(zero,faceusp(f0+kf,l,1), &
                                abs(faceusp(f0+kf,l,1))<=eps14),l=1,nr)
      !
      write (r1file,11) (merge(zero,face_xyz(nf)%v(l,kf,2), &
                                abs(face_xyz(nf)%v(l,kf,2))<=eps14),l=1,nr)
      write (r2file,11) (merge(zero,faceusp(f0+kf,l,2), &
                                abs(faceusp(f0+kf,l,2))<=eps14),l=1,nr)
      !
      write (d1file,11) (abs( face_xyz(nf)%v(l,kf,2) - &
                              face_xyz(nf)%v(l,kf,1)) , l=1,nr), &
                        norm2( face_xyz(nf)%v(1:nr,kf,2) - &
                               face_xyz(nf)%v(1:nr,kf,1) )
      write (d2file,11) (abs( faceusp(f0+kf,l,2) - &
                              faceusp(f0+kf,l,1)) , l=1,nr), &
                        norm2( faceusp(f0+kf,1:nr,2) - &
                               faceusp(f0+kf,1:nr,1) )
      !
    end do
    !
  end do
  !
  do kf = 0,8
    if (any(face(:)%left%rotation == kf)) then
      write (iout,*) "there exists  left%rotation = ",kf
    end if
    if (any(face(:)%right%rotation == kf)) then
      write (iout,*) "there exists right%rotation = ",kf
    end if
  end do
  !
  deallocate ( face_xyz , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_xyz",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( faceusp , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faceusp",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
 !call stop_gfr(stop_mpi,pname,__LINE__,__FILE__, &
 !              "After dumping flux point coordinates")
  !
  ! Format Statements
  !
 !1 format (i4,*(es23.15))
  1 format (i4,100(es12.4))
 !1 format (i4,*(es12.4))
  2 format (/,"#######################################################",/, &
              "#################  FACE ",i0,"  ###########################",/, &
              "#######################################################",/)
 11 format (1x,100(es12.4))
!11 format (1x,*(es12.4))
 12 format ("FACE ",i0,4x,a,:,4x,"rotation = ",i0)
  !
end subroutine check_face_flux_point_orderings
!
!###############################################################################
!
subroutine check_periodic_coordinates()
  !
  !.. Use Statements ..
  use geovar,  only : nr,bc_grps
  use flowvar, only : facexyz
  !
  !.. Local Scalars ..
  integer :: i,k,l,n,n1,n2
  integer :: lpfile,rpfile
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_periodic_coordinates"
  !
continue
  !
  lpfile = 3101 + mypnum
  rpfile = 3201 + mypnum
  !
  do i = 1,size(bc_grps)
    !
    if (bc_grps(i)%bc_type == bc_periodic) then
      !
      write (lpfile,10) i
      write (rpfile,10) i
      !
      n1 = bc_grps(i)%beg_fp
      n2 = bc_grps(i)%end_fp
      !
      do n = n1,n2
        !
        k = bc_grps(i)%periodic_list(n-n1+1)
        !
        write (lpfile,11) n,(facexyz(l,n), l=1,nr)
        write (rpfile,11) k,(facexyz(l,k), l=1,nr)
        !
        if (mod(n-n1+1,10) == 0) then
          write (lpfile,*)
          write (rpfile,*)
        end if
        !
      end do
      !
    end if
    !
  end do
  !
  10 format ("BC_GRP(",i0,")",/,60("="))
  11 format (i8,100(es12.4))
  !
end subroutine check_periodic_coordinates
!
!###############################################################################
!
subroutine output_element_solpts_and_flxpts()
  !
  !.. Use Statements ..
  use element_type_mod,  only : std_elem
  use order_mod,         only : n_min_geom,n_min_order
  use order_mod,         only : n_max_geom,n_max_order
  use order_mod,         only : geom_is_used,geom_solpts
  use geovar,            only : fp
  use interpolation_mod, only : interp
  !
  !.. Local Scalars ..
  integer :: this_geom,this_order,face_geom
  integer :: ndim,io1,ierr,n,l,k0,nfp,kf,k,lf,i,j
  real(wp) :: a,val
  !
  character(len=100) :: f1
  !
  !.. Local Arrays ..
  integer :: ip(1:8)
  real(wp), allocatable :: face_val(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_element_solpts_and_flxpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  geom_loop: &
  do this_geom = n_min_geom,n_max_geom
    !
    if (.not. geom_is_used(this_geom)) cycle geom_loop
    if (all(this_geom /= Geom_3D)) cycle geom_loop
    !
    ndim = geom_dimen(this_geom)
    !
    order_loop: &
    do this_order = n_min_order,n_max_order
      !
      write (f1,1) Geom_Short_Name(this_geom),this_order
      !
      open (newunit=io1,file=trim(adjustl(f1)),status="replace",action="write",&
                        iostat=ierr,iomsg=error_message)
      call io_error(pname,f1,1,__LINE__,__FILE__,ierr,error_message)
      !
      ! Write out the nodes defining this standard element
      !
      write (io1,10) Geom_Short_Name(this_geom),geom_nodes(this_geom)
      !
      if (this_geom == Geom_Tetr) then
        ip(:) = [1,2,3,3,4,4,4,4]
      else if (this_geom == Geom_Pyra) then
        ip(:) = [1,2,3,4,5,5,5,5]
      else if (this_geom == Geom_Pris) then
        ip(:) = [1,2,3,3,4,5,6,6]
      else if (this_geom == Geom_Hexa) then
        ip(:) = [1,2,3,4,5,6,7,8]
      end if
      do n = 1,geom_nodes(this_geom)
        write (io1,11) (ref_elem_nodes_xyz(l,n,this_geom), l=1,ndim),n,n
      end do
      write (io1,12) (ip(n), n=1,8)
      !
      ! Write out the solution points for this geometry
      !
      nfp = size(std_elem(this_geom,this_order)%pts,dim=2)
      !
      write (io1,20) Geom_Short_Name(this_geom),nfp
      !
      do n = 1,nfp
        write (io1,11) (std_elem(this_geom,this_order)%pts(l,n), l=1,ndim),n,n
      end do
      !
      ! Write out the face points stored in std_elem
      !
      k0 = 0 ! offset for start of face points for current face
      do n = 1,geom_faces(this_geom)
        !
        face_geom = cell_face_geom(n,this_geom) ! geometry of this face
        !
        nfp = geom_solpts(face_geom,this_order) ! # of face points for this face
        !
        write (io1,30) Geom_Short_Name(this_geom),n, nfp
        !
        do kf = 1,nfp
          k = k0 + kf
          write (io1,11) &
            (std_elem(this_geom,this_order)%face_pts(l,k), l=1,ndim), kf, k
        end do
        !
        k0 = k0 + nfp ! update offset for start of face points for next face
        !
      end do
     !!
     !! Write out the face points found by interpolating the solution points
     !!
     !face_val = init_array(ndim,zero) ! F2003 AUTO-REALLOC
     !!
     !k0 = 1
     !do n = 1,geom_faces(this_geom)
     !  !
     !  face_geom = cell_face_geom(n,this_geom) ! geometry of this face
     !  !
     !  nfp = geom_solpts(face_geom,this_order) ! # of face points for this face
     !  !
     !  write (io1,40) Geom_Short_Name(this_geom),n, nfp
     !  !
     !  do kf = 1,nfp
     !    !
     !    lf = k0 + fp%geord(face_geom,this_order)%rot(1)%pts(kf)
     !    !
     !    do l = 1,ndim
     !      val = zero
     !      do i = interp(this_geom,this_order)% &
     !                       toFace(this_order)%row_ptr(lf)+1, &
     !             interp(this_geom,this_order)% &
     !                       toFace(this_order)%row_ptr(lf+1)
     !        !
     !        j = interp(this_geom,this_order)%toFace(this_order)%row_idx(i)
     !        a = interp(this_geom,this_order)%toFace(this_order)%row_val(i)
     !        !
     !        val = val + a * std_elem(this_geom,this_order)%pts(l,j)
     !        !
     !      end do
     !      face_val(l) = val
     !    end do
     !    !
     !    write (io1,11) (face_val(l), l=1,ndim), kf, lf
     !    !
     !  end do
     !  !
     !  k0 = k0 + nfp
     !  !
     !end do
      !
      close (io1,iostat=ierr,iomsg=error_message)
      call io_error(pname,f1,2,__LINE__,__FILE__,ierr,error_message)
      !
    end do order_loop
    !
  end do geom_loop
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1  format (a4,".p",i0,".tec.dat")
  10 format ('VARIABLES = "X", "Y", "Z", "local idx", "global idx"',/, &
             'ZONE T="',a,' standard element", N=',i0,', E=1, ', &
             'ZONETYPE=FEBRICK, DATAPACKING=POINT')
  11 format (3es16.8,2x,i4,2x,i4)
  12 format (10(1x,i0))
  20 format ('ZONE T="',a,' solution points", I=',i0,',J=1,K=1, ', &
             'ZONETYPE=ORDERED, DATAPACKING=POINT')
  30 format ('ZONE T="',a,' face ',i1,' points", I=',i0,',J=1,K=1, ', &
             'ZONETYPE=ORDERED, DATAPACKING=POINT')
  40 format ('ZONE T="',a,' face ',i1,' interpolated points", I=',i0, &
             ',J=1,K=1, ZONETYPE=ORDERED, DATAPACKING=POINT')
  !
end subroutine output_element_solpts_and_flxpts
!
!###############################################################################
!
pure subroutine transfer_periodic(n1,n2,face_cv,periodic_list)
  !
  integer,     intent(in) :: n1,n2
  integer,     intent(in) :: periodic_list(:)
  real(wp), intent(inout) :: face_cv(:,:,:)
  !
  integer :: n
  !
continue
  !
  do n = n1,n2
    face_cv(n,1,2) = face_cv(periodic_list(n-n1+1),1,1)
    face_cv(n,2,2) = face_cv(periodic_list(n-n1+1),2,1)
    face_cv(n,3,2) = face_cv(periodic_list(n-n1+1),3,1)
  end do
  !
end subroutine transfer_periodic
!
!###############################################################################
!
subroutine check_cell_flux_point_orderings()
  !
  !.. Use Statements ..
  use order_mod, only : maxSP,geom_solpts
  use geovar, only : ncell,cell
  use geovar, only : nfbnd,nface,face,fp
  use geovar, only : nr,xyz,bface
  use interpolation_mod, only : interp
  use flowvar, only : face_var,faceusp
  use parallel_mod, only : exchange_faceusp,wait_faceusp
  use parallel_mod, only : exchange_faceusp_dt,wait_faceusp_dt
  !
  !.. Local Scalars ..
  integer :: i,j,nf,nc,np,n1,n2,l,n,sf,kf,lf,pf,nfp,ierr
  integer :: this_geom,this_order
  integer :: face_geom,face_order,f0,ntfp
  integer :: rotation,fp_offset
  integer :: l1file,r1file,d1file
  integer :: l2file,r2file,d2file
  real(wp) :: a
  character(len=25) :: face_string
  character(len=100) :: array_name
  !
  !.. Local Arrays ..
  real(wp) :: txyz(1:maxSP,1:nr)
  real(wp) :: xyz_val(1:nr)
  !
  !.. Local Allocatable Arrays ..
  type(face_var), allocatable :: face_xyz(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_cell_flux_point_orderings"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  l1file = 1401 + mypnum
  r1file = 1501 + mypnum
  d1file = 1601 + mypnum
  l2file = 2401 + mypnum
  r2file = 2501 + mypnum
  d2file = 2601 + mypnum
  !
  write (l1file,*)
  write (r1file,*)
  write (d1file,*)
  write (l2file,*)
  write (r2file,*)
  write (d2file,*)
  !
  allocate ( face_xyz(1:nface) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_xyz",1,__LINE__,__FILE__,ierr,error_message)
  !
  do nf = 1,nface
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    nfp = geom_solpts(face_geom,face_order)
    allocate ( face_xyz(nf)%v(1:nq,1:nfp,1:2) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,'("face_xyz(",i0,")%v")') nf
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  end do
  !
  ntfp = face(nface)%end_fp
  allocate ( faceusp(1:ntfp,1:nq,1:2) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faceusp",1,__LINE__,__FILE__,ierr,error_message)
  !
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    np = n2-n1+1
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    do kf = 1,np
      do l = 1,nr
        txyz(kf,l) = xyz(l,kf+n1-1)
      end do
    end do
    !
    do n = 1,size(cell(nc)%face)
      !
      nf = cell(nc)%face(n)%idx
      sf = cell(nc)%face(n)%side
      !
      face_geom = cell(nc)%face(n)%geom
      face_order = cell(nc)%face(n)%order
      !
      do kf = 1,geom_solpts(face_geom,face_order)
        !
        lf = cell(nc)%face(n)%fp_offset + fp%geord(face_geom,face_order)% &
                                          rot(cell(nc)%face(n)%rotation)% &
                                          pts(kf)
        !
        xyz_val(1:nr) = zero
        do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1, &
               interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
          !
          j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          do l = 1,nr
            xyz_val(l) = xyz_val(l) + a * txyz(j,l)
          end do
          !
        end do
        !
        face_xyz(nf)%v(1:nr,kf,sf) = xyz_val(1:nr)
        !
      end do
      !
    end do
    !
  end do
  !
  do nf = 1,nfbnd
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    nfp = geom_solpts(face_geom,face_order)
    !
    if (bface(1,nf) == bc_periodic) then
      !
      pf = bface(5,nf)
      !
      rotation = face(nf)%right%rotation
      fp_offset = face(nf)%right%fp_offset
      !
      do kf = 1,nfp
        do l = 1,nr
          lf = fp_offset + fp%geord(face_geom,face_order)% &
                                            rot(rotation)% &
                           pts(kf)
          face_xyz(nf)%v(l,kf,2) = face_xyz(pf)%v(l,lf,1)
        end do
      end do
      !
    else if (bface(1,nf) /= bc_cpu_bnd) then
      !
      do kf = 1,nfp
        do l = 1,nr
          face_xyz(nf)%v(l,kf,2) = face_xyz(nf)%v(l,kf,1)
        end do
      end do
      !
    end if
    !
  end do
  !
  faceusp(:,:,:) = zero
  !
  do nf = 1,nface
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    do i = 1,2
      do kf = 1,nfp
        do l = 1,nr
          faceusp(f0+kf,l,i) = face_xyz(nf)%v(l,kf,i)
        end do
      end do
    end do
    !
  end do
  !
  if (ncpu > 1) then
    !
    call exchange_faceusp_dt(face_xyz)
    call wait_faceusp_dt(face_xyz)
    !
    call exchange_faceusp
    call wait_faceusp
    !
  end if
  !
  do nf = 1,nface
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    face_string = ""
    if (nf <= nfbnd) then
      if (bface(1,nf) == bc_periodic) then
        face_string = "(periodic face)"
      else if (bface(1,nf) == bc_cpu_bnd) then
        face_string = "(cpu boundary)"
      else
        face_string = "(non-communication BC)"
      end if
    end if
    !
    write (l1file,12) nf-nfbnd,trim(adjustl(face_string))
    write (r1file,12) nf-nfbnd,trim(adjustl(face_string))
    write (d1file,12) nf-nfbnd,trim(adjustl(face_string))
    write (l2file,12) nf-nfbnd,trim(adjustl(face_string))
    write (r2file,12) nf-nfbnd,trim(adjustl(face_string))
    write (d2file,12) nf-nfbnd,trim(adjustl(face_string))
    !
    do kf = 1,nfp
      !
      write (l1file,11) (merge(zero,face_xyz(nf)%v(l,kf,1), &
                                abs(face_xyz(nf)%v(l,kf,1))<=eps14),l=1,nr)
      write (l2file,11) (merge(zero,faceusp(f0+kf,l,1), &
                                abs(faceusp(f0+kf,l,1))<=eps14),l=1,nr)
      !
      write (r1file,11) (merge(zero,face_xyz(nf)%v(l,kf,2), &
                                abs(face_xyz(nf)%v(l,kf,2))<=eps14),l=1,nr)
      write (r2file,11) (merge(zero,faceusp(f0+kf,l,2), &
                                abs(faceusp(f0+kf,l,2))<=eps14),l=1,nr)
      !
      write (d1file,11) (abs( face_xyz(nf)%v(l,kf,2) - &
                              face_xyz(nf)%v(l,kf,1) ) , l=1,nr), &
                       norm2( face_xyz(nf)%v(1:nr,kf,2) - &
                              face_xyz(nf)%v(1:nr,kf,1) )
      write (d2file,11) (abs( faceusp(f0+kf,l,2) - &
                              faceusp(f0+kf,l,1) ) , l=1,nr), &
                       norm2( faceusp(f0+kf,1:nr,2) - &
                              faceusp(f0+kf,1:nr,1) )
      !
    end do
    !
  end do
  !
  deallocate ( face_xyz , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_xyz",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( faceusp , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faceusp",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__, &
                "After dumping flux point coordinates")
  !
  ! Format Statements
  !
 !1 format (i4,*(es23.15))
  1 format (i4,100(es12.4))
 !1 format (i4,*(es12.4))
  2 format (/,"#######################################################",/, &
              "#################  CELL ",i0,"  ###########################",/, &
              "#######################################################",/)
 11 format (1x,100(es12.4))
!11 format (1x,*(es12.4))
 12 format ("FACE ",i0,4x,a)
  !
end subroutine check_cell_flux_point_orderings
!
!###############################################################################
!
subroutine test_communication()
  !
  !.. Use Statements ..
  use geovar, only : nfbnd,nface,face,bface
  use flowvar, only : face_var,faceusp
  use parallel_mod, only : exchange_faceusp,wait_faceusp
  use parallel_mod, only : exchange_faceusp_dt,wait_faceusp_dt
  !
  !.. Local Scalars ..
  integer :: i,nf,l,kf,pf,nfp,f0,ntfp,ierr
  integer :: l1file,r1file,d1file
  integer :: l2file,r2file,d2file
  character(len=25) :: face_string
  character(len=100) :: array_name
  !
  !.. Local Allocatable Arrays ..
  type(face_var), allocatable :: face_xyz(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "test_communication"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  l1file = 1701 + mypnum
  r1file = 1801 + mypnum
  d1file = 1901 + mypnum
  l2file = 2701 + mypnum
  r2file = 2801 + mypnum
  d2file = 2901 + mypnum
  !
  write (l1file,*)
  write (r1file,*)
  write (d1file,*)
  write (l2file,*)
  write (r2file,*)
  write (d2file,*)
  !
  allocate ( face_xyz(1:nface) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_xyz",1,__LINE__,__FILE__,ierr,error_message)
  !
  do nf = 1,nface
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    allocate ( face_xyz(nf)%v(1:nq,1:nfp,1:2) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,'("face_xyz(",i0,")%v")') nf
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  end do
  !
  ntfp = face(nface)%end_fp
  allocate ( faceusp(1:ntfp,1:nq,1:2) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faceusp",1,__LINE__,__FILE__,ierr,error_message)
  !
  i = 0
  !
  do nf = 1,nfbnd
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    do kf = 1,nfp
      do l = 1,nq
        i = i + 1
        face_xyz(nf)%v(l,kf,1) = real(i,kind=wp)
      end do
    end do
    !
  end do
  !
  do nf = 1,nfbnd
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    if (bface(1,nf) == bc_periodic) then
      !
      pf = bface(5,nf)
      !
      do kf = 1,nfp
        do l = 1,nq
          face_xyz(nf)%v(l,kf,2) = face_xyz(pf)%v(l,kf,1)
        end do
      end do
      !
    else if (bface(1,nf) /= bc_cpu_bnd) then
      !
      do kf = 1,nfp
        do l = 1,nq
          face_xyz(nf)%v(l,kf,2) = face_xyz(nf)%v(l,kf,1)
        end do
      end do
      !
    end if
    !
  end do
  !
  faceusp(:,:,:) = zero
  !
  do nf = 1,nfbnd
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    do i = 1,2
      do kf = 1,nfp
        do l = 1,nq
          faceusp(f0+kf,l,i) = face_xyz(nf)%v(l,kf,i)
        end do
      end do
    end do
    !
  end do
  !
  if (ncpu > 1) then
    !
    call exchange_faceusp_dt(face_xyz)
    call wait_faceusp_dt(face_xyz)
    !
    call exchange_faceusp
    call wait_faceusp
    !
  end if
  !
  do nf = 1,nfbnd
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    if (bface(1,nf) == bc_periodic) then
      face_string = "(periodic face)"
    else if (bface(1,nf) == bc_cpu_bnd) then
      face_string = "(cpu boundary)"
    else
      face_string = "(non-communication BC)"
    end if
    !
    write (l1file,12) nf-nfbnd,trim(adjustl(face_string))
    write (r1file,12) nf-nfbnd,trim(adjustl(face_string))
    write (d1file,12) nf-nfbnd,trim(adjustl(face_string))
    write (l2file,12) nf-nfbnd,trim(adjustl(face_string))
    write (r2file,12) nf-nfbnd,trim(adjustl(face_string))
    write (d2file,12) nf-nfbnd,trim(adjustl(face_string))
    !
    do kf = 1,nfp
      !
      write (l1file,13) (nint(face_xyz(nf)%v(l,kf,1)),l=1,nq)
      write (l2file,13) (nint(faceusp(f0+kf,l,1)),l=1,nq)
      !
      write (r1file,13) (nint(face_xyz(nf)%v(l,kf,2)),l=1,nq)
      write (r2file,13) (nint(faceusp(f0+kf,l,2)),l=1,nq)
      !
      write (d1file,11) (abs( face_xyz(nf)%v(l,kf,2) - &
                              face_xyz(nf)%v(l,kf,1) ) , l=1,nq), &
                       norm2( face_xyz(nf)%v(1:nq,kf,2) - &
                              face_xyz(nf)%v(1:nq,kf,1) )
      write (d2file,11) (abs( faceusp(f0+kf,l,2) - &
                              faceusp(f0+kf,l,1) ) , l=1,nq), &
                       norm2( faceusp(f0+kf,1:nq,2) - &
                              faceusp(f0+kf,1:nq,1) )
      !
    end do
    !
  end do
  !
  deallocate ( face_xyz , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_xyz",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( faceusp , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faceusp",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
 !1 format (i4,*(es23.15))
  1 format (i4,100(es12.4))
 !1 format (i4,*(es12.4))
  2 format (/,"#######################################################",/, &
              "#################  CELL ",i0,"  ###########################",/, &
              "#######################################################",/)
 11 format (1x,100(es12.4))
!11 format (1x,*(es12.4))
 12 format ("FACE ",i0,4x,a)
 13 format (1x,100(i9))
  !
end subroutine test_communication
!
!###############################################################################
!
subroutine write_final_error
  !
  use geovar, only : nr,n_global_cell
  use geovar, only : n_global_solpts,global_pinc_ptr
  use ovar,   only : itestcase,flux_divergence
  use ovar,   only : uverr
  use order_mod, only : n_order
  !
  !.. Local Scalars ..
  integer :: n,k
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_final_error"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  n = nint(uverr(1,4,2))
  k = nint(uverr(1,5,2))
  write (iout,31) itestcase,flux_divergence, &
                  n_order,n_global_cell,n_global_solpts, &
                  n,global_pinc_ptr(n)+k
  !
  write (iout,32) uverr(nec,1,2),  uverr(nec,3,2),  uverr(nec,2,2)
  write (iout,33) uverr(nmx,1,2),  uverr(nmx,3,2),  uverr(nmx,2,2)
  if (nr > 1) &
  write (iout,34) uverr(nmy,1,2),  uverr(nmy,3,2),  uverr(nmy,2,2)
  if (nr > 2) &
  write (iout,35) uverr(nmz,1,2),  uverr(nmz,3,2),  uverr(nmz,2,2)
  write (iout,36) uverr(nq+1,1,1), uverr(nq+1,3,1), uverr(nq+1,2,1)
  write (iout,37) uverr(nq+1,1,2), uverr(nq+1,3,2), uverr(nq+1,2,2)
  !
  call debug_timer(leaving_procedure,pname)
  !
  31 format (/,'itestcase = ',i0,/, &
               'flux_divergence  = ',i0,//, &
               'Order of solution polynomial = ',i10,/, &
               'Number of grid cells         = ',i10,/, &
               'Number of solution points    = ',i10,//, &
               'Max density error in cell = ',i10,/, &
               '        at solution point = ',i10,/)
  32 format (/,'  Error  -      L1           Max            L2',/, &
               ' --------------------------------------------------',/, &
               ' Density ',3es14.6)
  33 format (  '   Vx    ',3es14.6)
  34 format (  '   Vy    ',3es14.6)
  35 format (  '   Vz    ',3es14.6)
  36 format (  ' Entropy ',5es14.6)
  37 format (  ' Vel Mag ',5es14.6,/, &
        ' --------------------------------------------------',/)
  !
end subroutine write_final_error
!
!###############################################################################
!
subroutine write_residual_file(iflag,icount)
  !
  use geovar,    only : nr,n_global_solpts
  use ovar,      only : elapsed_wall_time,work_units
  use ovar,      only : uverr,total_volume
  use ovar,      only : time,time_ref,time_scaling_factor
  use ovar,      only : output_time_averaging,ave_start_time
  use ovar,      only : itestcase
  use order_mod, only : p_order
  use thermo_mod, only : species_name
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iflag
  integer, intent(in) :: icount
  !
  !.. Local Scalars ..
  integer  :: i,ierr
  real(wp) :: comp_cost
  character(len=25) :: suffix
  !
  !.. Local Saved Scalars ..
  integer,     save :: iores
  logical(lk), save :: file_is_opened = fals
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq+1) :: u2,v2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_residual_file"
  character(len=*), parameter :: fname = "residual.dat"
  character(len=*), parameter :: uvw(1:3) = ["u","v","w"]
 !logical(lk),      parameter :: write_work_units = true
  logical(lk),      parameter :: write_work_units = fals
  !
!#define OUTPUT_LINF_AND_L1_RESIDUALS
#ifdef OUTPUT_LINF_AND_L1_RESIDUALS
  integer, save :: io1,ioinf
  character(len=*), parameter :: fname_linf = "residual.linf.dat"
  character(len=*), parameter :: fname_l1 = "residual.l1.dat"
#endif
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Open the results file only if iflag is +1
  !
  if (.not. file_is_opened) then
    !
    open (newunit=iores,file=fname,status="replace",action="write", &
                        iostat=ierr,iomsg=error_message)
    call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
    !
    write (iores,1)
    !
#ifdef OUTPUT_LINF_AND_L1_RESIDUALS
    open (newunit=ioinf,file=fname_linf,status="replace",action="write", &
                        iostat=ierr,iomsg=error_message)
    call io_error(pname,fname_linf,1,__LINE__,__FILE__,ierr,error_message)
    !
    write (ioinf,1)
    !
    open (newunit=io1,file=fname_l1,status="replace",action="write", &
                        iostat=ierr,iomsg=error_message)
    call io_error(pname,fname_l1,1,__LINE__,__FILE__,ierr,error_message)
    !
    write (io1,1)
#endif
    !
    suffix = '",'
    if (output_time_averaging) suffix = '<sub>ave</sub>",'
   !if (output_time_averaging) suffix = "<sub>ave</sub>" // trim(suffix)
    !
    if (ncs > 0) then
      do i = 1,ncs
        write (iores,2) "<greek>r</greek><sub>" // &
                        trim(adjustl(species_name(i))) // &
                        "</sub>" // trim(suffix)
      end do
    end if
    !
    write (iores,2) "<greek>r</greek>" // trim(suffix)
    do i = 1,nr
      write (iores,3) "<greek>r</greek>" // uvw(i) // trim(suffix)
    end do
    write (iores,2) "Energy" // trim(suffix)
    write (iores,2) "Entropy" // trim(suffix)
    !
    if (ncs > 0) then
      do i = 1,ncs
        write (iores,2) "Y<sub>" // &
                        trim(adjustl(species_name(i))) // &
                        "</sub>" // trim(suffix)
      end do
    end if
    !
    write (iores,2) "T" // trim(suffix)
    if (abs(itestcase) == THX_Phase_3) then
      write (iores,2) "Nozzle mdot" // '",'
      write (iores,2) "Injector mdot" // '",'
      write (iores,2) "Blowing Ratio" // '",'
    else
      do i = 1,nr
        write (iores,2) uvw(i) // trim(suffix)
      end do
    end if
    write (iores,2) "p" // trim(suffix)
    write (iores,2) "||V||" // trim(suffix)
    if (write_work_units) then
      write (iores,2) "Work units" // '",'
    else
      write (iores,2) "Time" // '",'
    end if
    write (iores,2) "Wall Time" // '",'
    if     ( nr == 1 ) then
      write (iores,2) "xloc" // '"'
    elseif ( nr == 2 ) then
      write (iores,2) "xloc" // '",'
      write (iores,2) "yloc" // '"'
    else
      write (iores,2) "xloc" // '",'
      write (iores,2) "yloc" // '",'
      write (iores,2) "zloc" // '"'
    endif
    !
    write (iores,4) p_order, one/sqrt(real(n_global_solpts,kind=wp))
    !
#ifdef OUTPUT_LINF_AND_L1_RESIDUALS
    if (ncs > 0) then
      do i = 1,ncs
        write (ioinf,2) "<greek>r</greek><sub>" // &
                        trim(adjustl(species_name(i))) // &
                        "</sub>" // trim(suffix)
      end do
    end if
    !
    write (ioinf,2) "<greek>r</greek>" // trim(suffix)
    do i = 1,nr
      write (ioinf,3) "<greek>r</greek>" // uvw(i) // trim(suffix)
    end do
    write (ioinf,2) "Energy" // trim(suffix)
    write (ioinf,2) "Entropy" // trim(suffix)
    !
    if (ncs > 0) then
      do i = 1,ncs
        write (ioinf,2) "Y<sub>" // &
                        trim(adjustl(species_name(i))) // &
                        "</sub>" // trim(suffix)
      end do
    end if
    !
    write (ioinf,2) "T" // trim(suffix)
    do i = 1,nr
      write (ioinf,2) uvw(i) // trim(suffix)
    end do
    write (ioinf,2) "p" // trim(suffix)
    write (ioinf,2) "||V||" // trim(suffix)
    if (write_work_units) then
      write (ioinf,2) "Work units",'"'
    else
      write (ioinf,2) "Time",'"'
    end if
    !
    write (ioinf,4) p_order, one/sqrt(real(n_global_solpts,kind=wp))
    !
    if (ncs > 0) then
      do i = 1,ncs
        write (io1,2) "<greek>r</greek><sub>" // &
                        trim(adjustl(species_name(i))) // &
                        "</sub>" // trim(suffix)
      end do
    end if
    !
    write (io1,2) "<greek>r</greek>" // trim(suffix)
    do i = 1,nr
      write (io1,3) "<greek>r</greek>" // uvw(i) // trim(suffix)
    end do
    write (io1,2) "Energy" // trim(suffix)
    write (io1,2) "Entropy" // trim(suffix)
    !
    if (ncs > 0) then
      do i = 1,ncs
        write (io1,2) "Y<sub>" // &
                        trim(adjustl(species_name(i))) // &
                        "</sub>" // trim(suffix)
      end do
    end if
    !
    write (io1,2) "T" // trim(suffix)
    do i = 1,nr
      write (io1,2) uvw(i) // trim(suffix)
    end do
    write (io1,2) "p" // trim(suffix)
    write (io1,2) "||V||" // trim(suffix)
    if (write_work_units) then
      write (io1,2) "Work units",'"'
    else
      write (io1,2) "Time",'"'
    end if
    !
    write (io1,4) p_order, one/sqrt(real(n_global_solpts,kind=wp))
    !
#endif
    !
    file_is_opened = true
    !
  end if
  !
  if (write_work_units) then
    comp_cost = log10(merge(work_units,one,work_units > zero))
  else
    comp_cost = time*time_ref*time_scaling_factor
    if (output_time_averaging) then
      comp_cost = comp_cost - ave_start_time*time_ref*time_scaling_factor
    end if
  end if
  !
  ! L2 Norm
  u2(:) = merge( uverr(:,2,1) , epsilon(zero)/total_volume , uverr(:,2,1)>zero )
  v2(:) = merge( uverr(:,2,2) , epsilon(zero)/total_volume , uverr(:,2,2)>zero )
  !
  do i = 1,nq+1
    u2(i) = log10(u2(i))
    v2(i) = log10(v2(i))
  end do
  !
  if (abs(itestcase) == THX_Phase_3) then
    v2(nmx) = uverr(nmx,2,2)
    v2(nmy) = uverr(nmy,2,2)
    v2(nmz) = uverr(nmz,2,2)
  end if
  !
  write (iores,5) icount,(u2(i),i=1,nq+1), (v2(i),i=1,nq+1), &
                  comp_cost, elapsed_wall_time, &
                  real(uverr(1,6:6+nr-1,1),sp)
  flush (iores)
  !
#ifdef OUTPUT_LINF_AND_L1_RESIDUALS
  ! Infinity Norm
  u2(:) = merge( uverr(:,3,1) , epsilon(zero)/total_volume , uverr(:,3,1)>zero )
  v2(:) = merge( uverr(:,3,2) , epsilon(zero)/total_volume , uverr(:,3,2)>zero )
  !
  write (ioinf,5) icount,(log10(u2(i)),i=1,nq+1), &
                         (log10(v2(i)),i=1,nq+1), &
                         comp_cost
  flush (ioinf)
  !
  ! L1 Norm
  u2(:) = merge( uverr(:,1,1) , epsilon(zero)/total_volume , uverr(:,1,1)>zero )
  v2(:) = merge( uverr(:,1,2) , epsilon(zero)/total_volume , uverr(:,1,2)>zero )
  !
  write (io1,5) icount,(log10(u2(i)),i=1,nq+1), &
                       (log10(v2(i)),i=1,nq+1), &
                       comp_cost
  flush (io1)
#endif
  !
  ! Close the results file only if iflag is -1
  !
  if (iflag == -1) then
    close (iores,iostat=ierr,iomsg=error_message)
    call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
#ifdef OUTPUT_LINF_AND_L1_RESIDUALS
    close (ioinf,iostat=ierr,iomsg=error_message)
    call io_error(pname,fname_linf,2,__LINE__,__FILE__,ierr,error_message)
    close (io1,iostat=ierr,iomsg=error_message)
    call io_error(pname,fname_l1,2,__LINE__,__FILE__,ierr,error_message)
#endif
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ('VARIABLES= "Timestep",')
  2 format ('"',a,a)
  3 format ('"',a,a,a)
  4 format ('Zone T="Residual Error"',/, &
            'AUXDATA order="',i0,'"',/, &
            'AUXDATA h="',f10.6,'"')
  5 format (i10,100(es17.9))
 !5 format (i10,*(es17.9))
  !
end subroutine write_residual_file
!
!###############################################################################
!
subroutine write_summary_file(iter)
  !
  use geovar, only : n_global_solpts,n_global_cell
  !
  use order_mod, only : n_order
  !
  use ovar, only : output_summary_file
  use ovar, only : itestcase,grid_format,itcur,tmrst
  use ovar, only : time,time_ref,time_scaling_factor
  use ovar, only : output_time_averaging,ave_start_time
  use ovar, only : resource_list
  use ovar, only : gfr_cpu_time
  !
  use internal_grid_mod, only : write_channel_grid_info
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iter
  !
  !.. Local Scalars ..
  integer  :: i,io,ierr
  real(wp) :: out_tstar_ave,out_tstar_total,out_tstar_job,ave_tstar_job
  real(wp) :: out_time_ave,out_time_total,out_time_job,ave_time_job
  !
  character(len=14) :: c_out_tstar_ave,c_ave_tstar_job
  character(len=14) :: c_out_tstar_total,c_out_tstar_job
  character(len=14) :: c_out_time_ave,c_ave_time_job
  character(len=14) :: c_out_time_total,c_out_time_job
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_summary_file"
  character(len=*), parameter :: fname = "summary.log"
  !
continue
  !
  if (.not. output_summary_file) return
  !
  if (mypnum == glb_root) then
    !
    ! Final physical time completed across all jobs including restarts
    out_time_total = time*time_ref
    ! Final tstar completed across all jobs including restarts
    out_tstar_total = out_time_total*time_scaling_factor
    !
    ! Physical time completed during this job
    out_time_job = out_time_total - tmrst*time_ref
    ! Tstar completed during this job
    out_tstar_job = out_time_job*time_scaling_factor
    !
    if (output_time_averaging) then
      ! Total physical time averaged across all jobs including restarts
      out_time_ave = out_time_total - ave_start_time*time_ref
      ! Physical time averaged during this job
      ave_time_job = out_time_job - merge(ave_start_time*time_ref,zero, &
                                          ave_start_time>tmrst)
    else
      out_time_ave = zero
      ave_time_job = zero
    end if
    ! Total tstar time-averaged across all jobs including restarts
    out_tstar_ave = out_time_ave*time_scaling_factor
    ! Tstar time-averaged during this job
    ave_tstar_job = ave_time_job*time_scaling_factor
    !
    c_out_tstar_ave = g2c(out_tstar_ave)
    c_ave_tstar_job = g2c(ave_tstar_job)
    c_out_tstar_total = g2c(out_tstar_total)
    c_out_tstar_job = g2c(out_tstar_job)
    c_out_time_ave = g2c(out_time_ave)
    c_ave_time_job = g2c(ave_time_job)
    c_out_time_total = g2c(out_time_total)
    c_out_time_job = g2c(out_time_job)
    !
    open (newunit=io,file=fname,status="replace",action="write", &
                     iostat=ierr,iomsg=error_message)
    call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
    !
    write (io,2) g2c(n_global_cell),"global cells"
    write (io,2) g2c(n_global_solpts),"global DoF"
    write (io,1) n_order,"solution order"
    write (io,2) g2c(iter),"iterations completed with this job"
    write (io,2) g2c(itcur), &
               "total iterations completed across all jobs (including restarts)"
    !
    !
    !
    if (abs(time_scaling_factor-one) > eps6) then
      write (io,2) c_out_tstar_total,"final t-star reached across all jobs"
      write (io,2) c_out_tstar_job,"t-star elapsed during this job"
    end if
    write (io,2) c_out_time_total,"final physical time reached across all jobs"
    write (io,2) c_out_time_job, &
                              "simulation physical time elapsed during this job"
    !
    !
    !
    if (output_time_averaging) then
      if (abs(time_scaling_factor-one) > eps6) then
        write (io,2) c_out_tstar_ave, &
                      "total t-star that has been time-averaged across all jobs"
        write (io,2) c_ave_tstar_job, &
                      "total t-star time-averaged during this job"
      end if
      write (io,2) c_out_time_ave, &
               "total physical time that has been time-averaged across all jobs"
      write (io,2) c_ave_time_job, &
               "total physical time time-averaged during this job"
    end if
    !
    !
    !
    if (allocated(resource_list)) then
      write (io,3) (resource_list(i)%n, &
                    resource_list(i)%model, &
                    trim(merge("s"," ",resource_list(i)%n>1)), &
                    i=1,size(resource_list))
    end if
    !
    write (io,1) ncpu,"MPI process",trim(merge("es","  ",ncpu>1))
    !
    call gfr_cpu_time(-1,io_unit=io)
    write (io,*)
    !
    if (itestcase == Channel_Flow) then
      if (grid_format == Internally_Generated_Grid) then
        call write_channel_grid_info(iunit=io)
      end if
    end if
    !
    close (io,iostat=ierr,iomsg=error_message)
    call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
    !
  end if
  !
  ! Format Statements
  !
  1 format (2x,i14," : ",a,:,a)
  2 format (2x,a14," : ",a)
  3 format (/,6x,"Job run on : ",100(i0,1x,a," node",a,:,", "))
  !
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
end subroutine write_summary_file
!
!###############################################################################
!
pure &
function g2c_i4(val,width,digit,include_comma) result(cval)
  !
  use, intrinsic :: iso_fortran_env, only : lp => int32
  !
  !.. Formal Arguments ..
  integer(lp), intent(in) :: val
  !
  !.. Optional Arguments ..
  integer,     optional, intent(in) :: width
  integer,     optional, intent(in) :: digit
  logical(lk), optional, intent(in) :: include_comma
  !
  !.. Function Result ..
  character(len=:), allocatable :: cval
  !
  !.. Local Scalars ..
  integer :: w,d
  !
  logical(lk) :: add_commas
  !
  character(len=20)  :: cfmt
  character(len=100) :: ctmp
  !
  !.. Local Parameters ..
  logical(lk), parameter :: add_commas_default = true
  !
continue
  !
  w = 14
  d = -1
  if (present(width)) w = merge( width , w , width >= 0 )
  if (present(digit)) d = merge( digit , d , digit >= 0 )
  !
  add_commas = add_commas_default
  if (present(include_comma)) add_commas = include_comma
  !
  if (w == 0) then
    d = -1
  end if
  !
  ! Create the format string for reading val
  !
  write (cfmt,1) w
  !
  if (d >= 0) then
    write (cfmt,2) trim(adjustl(cfmt)),d
  end if
  !
  write (cfmt,3) trim(adjustl(cfmt))
  !
  ! Write val into ctmp
  write (ctmp,cfmt) val
  !
  ! Add thousands separators if add_commas is true
  if (add_commas) ctmp = add_thousands_separators(trim(ctmp))
  !
  ! Copy only the trimmed value of ctmp into cval as the function result.
  ! NOTE: Using only trim(ctmp) instead of trim(adjustl(ctmp))
  !       preserves the width value used to create ctmp.
  cval = trim(ctmp) ! F2003 AUTO-REALLOC
  !
  ! Format Statements
  !
  1 format ("(i",i0)
  2 format (a,".",i0)
  3 format (a,")")
  !
end function g2c_i4
!
!###############################################################################
!
pure &
function g2c_i8(val,width,digit,include_comma) result(cval)
  !
  use, intrinsic :: iso_fortran_env, only : lp => int64
  !
  !.. Formal Arguments ..
  integer(lp), intent(in) :: val
  !
  !.. Optional Arguments ..
  integer,     optional, intent(in) :: width
  integer,     optional, intent(in) :: digit
  logical(lk), optional, intent(in) :: include_comma
  !
  !.. Function Result ..
  character(len=:), allocatable :: cval
  !
  !.. Local Scalars ..
  integer :: w,d
  !
  logical(lk) :: add_commas
  !
  character(len=20)  :: cfmt
  character(len=100) :: ctmp
  !
  !.. Local Parameters ..
  logical(lk), parameter :: add_commas_default = true
  !
continue
  !
  w = 14
  d = -1
  if (present(width)) w = merge( width , w , width >= 0 )
  if (present(digit)) d = merge( digit , d , digit >= 0 )
  !
  add_commas = add_commas_default
  if (present(include_comma)) add_commas = include_comma
  !
  if (w == 0) then
    d = -1
  end if
  !
  ! Create the format string for reading val
  !
  write (cfmt,1) w
  !
  if (d >= 0) then
    write (cfmt,2) trim(adjustl(cfmt)),d
  end if
  !
  write (cfmt,3) trim(adjustl(cfmt))
  !
  ! Write val into ctmp
  write (ctmp,cfmt) val
  !
  ! Add thousands separators if add_commas is true
  if (add_commas) ctmp = add_thousands_separators(trim(ctmp))
  !
  ! Copy only the trimmed value of ctmp into cval as the function result.
  ! NOTE: Using only trim(ctmp) instead of trim(adjustl(ctmp))
  !       preserves the width value used to create ctmp.
  cval = trim(ctmp) ! F2003 AUTO-REALLOC
  !
  ! Format Statements
  !
  1 format ("(i",i0)
  2 format (a,".",i0)
  3 format (a,")")
  !
end function g2c_i8
!
!###############################################################################
!
pure &
function g2c_r4(val,width,digit,include_comma) result(cval)
  !
  use, intrinsic :: iso_fortran_env, only : lp => real32
  !
  !.. Formal Arguments ..
  real(lp), intent(in) :: val
  !
  !.. Optional Arguments ..
  integer,     optional, intent(in) :: width
  integer,     optional, intent(in) :: digit
  logical(lk), optional, intent(in) :: include_comma
  !
  !.. Function Result ..
  character(len=:), allocatable :: cval
  !
  !.. Local Scalars ..
  integer :: i,n,w,d
  !
  logical(lk) :: add_commas
  !
  character(len=20)  :: cfmt,cfmt1
  character(len=100) :: ctmp,cdec
  !
  character(len=:), allocatable :: csig
  !
  !.. Local Parameters ..
  logical(lk), parameter :: add_commas_default = true
  !
continue
  !
  w = 14
  d = 7
  if (present(width)) w = merge( width , w , width >= 0 )
  if (present(digit)) d = merge( digit , d , digit >= 0 )
  !
  add_commas = add_commas_default
  if (present(include_comma)) add_commas = include_comma
  !
  ! Extract the significand of the number as a character including
  ! the minus sign if negative, and thousand separators if desired
  csig = get_significand(val) ! F2003 AUTO-REALLOC
  !
  ! Get the length of the significand + 1 for the decimal point
  n = len_trim(adjustl(csig)) + 1
  !
  if (w-n >= d) then
    !
    ! This number will fit the specified width using F formatting
    !
    if (d == 0) then
      !
      write (cfmt,1) w-1
      write (ctmp,cfmt) trim(adjustl(csig))
      !
    else
      !
      ! Get the non-integer part of val and convert this to an integer value
      ! with d-digits, and then write that integer to a character string
      !
      i = nint( (val - real(int(val),kind=lp)) * (real(10,kind=lp)**d) )
      !
      write (cfmt1,2) d,d
      write (cdec,cfmt1) i
      !
      ! Write the significand and decimal part of val to ctmp
      !
      write (cfmt,3) w-d-1,d
      write (ctmp,cfmt) trim(adjustl(csig)),trim(adjustl(cdec))
      !
    end if
    !
  else
    !
    ! This number will not fit the specified width using F formatting,
    ! so use ES formatting instead.
    !
    ! Adjust the number of digits behind the decimal point if the specified
    ! number of digits will not be able to fit the specified width.
    !
    if (w-d < 7) d = w-7
    !
    ! Create the write format for val using ES formatting with w and d
    write (cfmt,4) w,d
    !
    ! Write val into ctmp
    write (ctmp,cfmt) val
    !
  end if
  !
  cval = trim(ctmp) ! F2003 AUTO-REALLOC
  !
  ! Format Statements
  !
  1 format ("(a",i0,",'.')")
  2 format ("(i",i0,".",i0,")")
  3 format ("(a",i0,",'.',a",i0,")")
  4 format ("(es",i0,".",i0,")")
  !
  contains
  !=============================================================================
  pure &
  function get_significand(a) result(return_value)
    real(lp), intent(in) :: a
    character(len=:), allocatable :: return_value
    character(len=100) :: c
    write (c,'(i0)') int(abs(a))
    if (add_commas) c = add_thousands_separators(trim(c))
    return_value = trim(adjustl(c)) ! F2003 AUTO-REALLOC
    if (a < real(0,kind=lp)) then
      return_value = '-' // return_value ! F2003 AUTO-REALLOC
    end if
  end function get_significand
  !=============================================================================
end function g2c_r4
!
!###############################################################################
!
pure &
function g2c_r8(val,width,digit,include_comma) result(cval)
  !
  use, intrinsic :: iso_fortran_env, only : lp => real64
  !
  !.. Formal Arguments ..
  real(lp), intent(in) :: val
  !
  !.. Optional Arguments ..
  integer,     optional, intent(in) :: width
  integer,     optional, intent(in) :: digit
  logical(lk), optional, intent(in) :: include_comma
  !
  !.. Function Result ..
  character(len=:), allocatable :: cval
  !
  !.. Local Scalars ..
  integer :: i,n,w,d
  !
  logical(lk) :: add_commas
  !
  character(len=20)  :: cfmt,cfmt1
  character(len=100) :: ctmp,cdec
  !
  character(len=:), allocatable :: csig
  !
  !.. Local Parameters ..
  logical(lk), parameter :: add_commas_default = true
  !
continue
  !
  w = 14
  d = 7
  if (present(width)) w = merge( width , w , width >= 0 )
  if (present(digit)) d = merge( digit , d , digit >= 0 )
  !
  add_commas = add_commas_default
  if (present(include_comma)) add_commas = include_comma
  !
  ! Extract the significand of the number as a character including
  ! the minus sign if negative, and thousand separators if desired
  csig = get_significand(val) ! F2003 AUTO-REALLOC
  !
  ! Get the length of the significand + 1 for the decimal point
  n = len_trim(adjustl(csig)) + 1
  !
  if (w-n >= d) then
    !
    ! This number will fit the specified width using F formatting
    !
    if (d == 0) then
      !
      write (cfmt,1) w-1
      write (ctmp,cfmt) trim(adjustl(csig))
      !
    else
      !
      ! Get the non-integer part of val and convert this to an integer value
      ! with d-digits, and then write that integer to a character string
      !
      i = nint( (val - real(int(val),kind=lp)) * (real(10,kind=lp)**d) )
      !
      write (cfmt1,2) d,d
      write (cdec,cfmt1) i
      !
      ! Write the significand and decimal part of val to ctmp
      !
      write (cfmt,3) w-d-1,d
      write (ctmp,cfmt) trim(adjustl(csig)),trim(adjustl(cdec))
      !
    end if
    !
  else
    !
    ! This number will not fit the specified width using F formatting,
    ! so use ES formatting instead.
    !
    ! Adjust the number of digits behind the decimal point if the specified
    ! number of digits will not be able to fit the specified width.
    !
    if (w-d < 7) d = w-7
    !
    ! Create the write format for val using ES formatting with w and d
    write (cfmt,4) w,d
    !
    ! Write val into ctmp
    write (ctmp,cfmt) val
    !
  end if
  !
  cval = trim(ctmp) ! F2003 AUTO-REALLOC
  !
  ! Format Statements
  !
  1 format ("(a",i0,",'.')")
  2 format ("(i",i0,".",i0,")")
  3 format ("(a",i0,",'.',a",i0,")")
  4 format ("(es",i0,".",i0,")")
  !
  contains
  !=============================================================================
  pure &
  function get_significand(a) result(return_value)
    real(lp), intent(in) :: a
    character(len=:), allocatable :: return_value
    character(len=100) :: c
    write (c,'(i0)') int(abs(a))
    if (add_commas) c = add_thousands_separators(trim(c))
    return_value = trim(adjustl(c)) ! F2003 AUTO-REALLOC
    if (a < real(0,kind=lp)) then
      return_value = '-' // return_value ! F2003 AUTO-REALLOC
    end if
  end function get_significand
  !=============================================================================
end function g2c_r8
!
!###############################################################################
!
pure &
function g2c_r16(val,width,digit,include_comma) result(cval)
  !
  use, intrinsic :: iso_fortran_env, only : lp => real128
  !
  !.. Formal Arguments ..
  real(lp), intent(in) :: val
  !
  !.. Optional Arguments ..
  integer,     optional, intent(in) :: width
  integer,     optional, intent(in) :: digit
  logical(lk), optional, intent(in) :: include_comma
  !
  !.. Function Result ..
  character(len=:), allocatable :: cval
  !
  !.. Local Scalars ..
  integer :: i,n,w,d
  !
  logical(lk) :: add_commas
  !
  character(len=20)  :: cfmt,cfmt1
  character(len=100) :: ctmp,cdec
  !
  character(len=:), allocatable :: csig
  !
  !.. Local Parameters ..
  logical(lk), parameter :: add_commas_default = true
  !
continue
  !
  w = 14
  d = 7
  if (present(width)) w = merge( width , w , width >= 0 )
  if (present(digit)) d = merge( digit , d , digit >= 0 )
  !
  add_commas = add_commas_default
  if (present(include_comma)) add_commas = include_comma
  !
  ! Extract the significand of the number as a character including
  ! the minus sign if negative, and thousand separators if desired
  csig = get_significand(val) ! F2003 AUTO-REALLOC
  !
  ! Get the length of the significand + 1 for the decimal point
  n = len_trim(adjustl(csig)) + 1
  !
  if (w-n >= d) then
    !
    ! This number will fit the specified width using F formatting
    !
    if (d == 0) then
      !
      write (cfmt,1) w-1
      write (ctmp,cfmt) trim(adjustl(csig))
      !
    else
      !
      ! Get the non-integer part of val and convert this to an integer value
      ! with d-digits, and then write that integer to a character string
      !
      i = nint( (val - real(int(val),kind=lp)) * (real(10,kind=lp)**d) )
      !
      write (cfmt1,2) d,d
      write (cdec,cfmt1) i
      !
      ! Write the significand and decimal part of val to ctmp
      !
      write (cfmt,3) w-d-1,d
      write (ctmp,cfmt) trim(adjustl(csig)),trim(adjustl(cdec))
      !
    end if
    !
  else
    !
    ! This number will not fit the specified width using F formatting,
    ! so use ES formatting instead.
    !
    ! Adjust the number of digits behind the decimal point if the specified
    ! number of digits will not be able to fit the specified width.
    !
    if (w-d < 7) d = w-7
    !
    ! Create the write format for val using ES formatting with w and d
    write (cfmt,4) w,d
    !
    ! Write val into ctmp
    write (ctmp,cfmt) val
    !
  end if
  !
  cval = trim(ctmp) ! F2003 AUTO-REALLOC
  !
  ! Format Statements
  !
  1 format ("(a",i0,",'.')")
  2 format ("(i",i0,".",i0,")")
  3 format ("(a",i0,",'.',a",i0,")")
  4 format ("(es",i0,".",i0,")")
  !
  contains
  !=============================================================================
  pure &
  function get_significand(a) result(return_value)
    real(lp), intent(in) :: a
    character(len=:), allocatable :: return_value
    character(len=100) :: c
    write (c,'(i0)') int(abs(a))
    if (add_commas) c = add_thousands_separators(trim(c))
    return_value = trim(adjustl(c)) ! F2003 AUTO-REALLOC
    if (a < real(0,kind=lp)) then
      return_value = '-' // return_value ! F2003 AUTO-REALLOC
    end if
  end function get_significand
  !=============================================================================
end function g2c_r16
!
!###############################################################################
!
pure &
function add_thousands_separators(cnum) result(return_value)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: cnum
  !
  !.. Function Result ..
  character(len=:), allocatable :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,n
  !
  character(len=100) :: cr1,cr2
  character(len=1)   :: c
  character(len=len_trim(cnum)) :: csep
  !
continue
  !
  cr1 = ''
  cr2 = ''
  csep = ''
  !
  ! Create a reversed version of the input cnum
  !
  i = 0
  do n = len_trim(cnum),1,-1
    i = i + 1
    cr1(i:i) = cnum(n:n)
  end do
  !
  ! Use cr1, the reversed version of cnum, to add a comma every
  ! third digit and save the result in cr2. If a space or minus
  ! sign is encountered at any time, pad the remaining part of
  ! cr2 with spaces and exit the loop.
  !
  i = 0
  j = 0
  create_cr2: &
  do n = 1,len_trim(cnum)
    i = i + 1
    j = j + 1
    c = cr1(n:n)
    cr2(i:i) = c
    if (c == '-' .or. c == ' ') then
      cr2(i+1:) = ' '
      exit create_cr2
    end if
    if (j == 3) then
      j = 0
      i = i + 1
      c = cr1(n+1:n+1)
      if (c /= '-' .and. c /= ' ') then
        cr2(i:i) = ','
      else
        cr2(i:i) = c
        cr2(i+1:) = ' '
        exit create_cr2
      end if
    end if
  end do create_cr2
  !
  ! Reverse the reversed number containing commas for thousands separators
  ! into cr1 to create a regular/forward version of the number including
  ! the thousands separators.
  !
  i = 0
  do n = len_trim(cr2),1,-1
    i = i + 1
    cr1(i:i) = cr2(n:n)
  end do
  !
  ! Copy the comma separated version of the number into the function result,
  ! and make sure to include left space padding if the cnum input is longer
  ! than the thousands separated result.
  !
  if (len_trim(cnum) > len_trim(cr1)) then
    csep = trim(cr1)
    csep = adjustr(csep)
    return_value = csep ! F2003 AUTO-REALLOC
  else
    return_value = trim(adjustl(cr1)) ! F2003 AUTO-REALLOC
  end if
  !
end function add_thousands_separators
!
!###############################################################################
!
include "Functions/last.f90"
#include "Functions/bc_priority.f90"
include "Functions/cell_vol.f90"
include "Functions/bc_string_to_integer.f90"
#include "Functions/check_for_bad_number.f90"
!
!###############################################################################
!
end module io_mod
!
!###############################################################################
!
