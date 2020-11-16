module io_mod
  !
#include <mpi_defs.h>
  !
  !.. Use Statements ..
  use module_kind_types
  use module_cgns_types
  use eqn_idx, only : nec,nmb,nme,nmx,nmy,nmz,nee,nq
  use eqn_idx, only : ncs,npres,ntemp,ngamm,nenth,ncp,nvisc,ncond,ng
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
  public :: write_parallel_cgns
  public :: write_residual_file
  public :: check_face_flux_point_orderings
  public :: check_cell_flux_point_orderings
  public :: test_communication
  public :: TaylorGreen_8s_vorticity_solution
  !
  ! ###########################
  ! ###  Module Parameters  ###
  ! ###########################
  !
  logical(lk), parameter :: check_gradients = fals
 !logical(lk), parameter :: check_gradients = true
  !
  ! ################################
  ! ###  Module Local Variables  ###
  ! ################################
  !
  !.. Output File Names ..
  character(len=*),   parameter :: cgns_base_fname = "solution.grid.cgns"
  character(len=*),   parameter :: cgns_grid_fname = "solution.grid.cgns"
  character(len=*),   parameter :: cgns_time_fname = "solution.time.cgns"
  character(len=170),      save :: cgns_dir = ""
  character(len=170),      save :: cgns_grid_file
  character(len=170),      save :: cgns_time_file
  !
  character(len=170),      save :: cgnsfile_fmt
  !
  !.. CGNS Path Variables ..
  integer(CBT),      save :: ifile_n
  integer(CBT),      save :: ibase_n
  integer(CBT),      save :: izone_n
  integer(CBT),      save :: isect_n
  integer(CBT), parameter :: iiter_n = 1_CBT
  !
  !.. CGNS Time Dependent Variables ..
  integer, save :: cgns_zone_counter
  !
  real(sp),             save, allocatable :: soltime(:)
  character(len=CGLEN), save, allocatable :: solname(:)
  character(len=CGLEN), save, allocatable :: solfile(:)
  !
  character(len=*), parameter :: base_name = "gfr_base"
  character(len=*), parameter :: zone_name = "Solution"
  !
  !.. CGNS Parallel Point Indices ..
  integer(CST), save :: n_my_output_pts
  integer(CST), save :: my_pts_beg
  integer(CST), save :: my_pts_end
  !
  integer(CST), dimension(1:3) :: zone_size
  !
  ! DISABLE_CGNS_OUTPUT: Private module variable to disable CGNS output. It
  !                      is currently only relevant to the Taylor-Green Vortex
  !                      problem and is more of a temporary addition until the
  !                      output logic/code becomes more robust.
  !
  logical(lk), save :: disable_cgns_output = fals
  !
  ! Generic interface for interpolate_for_output
  !
  interface interpolate_for_output
    module procedure interpolate_for_output_1d, &
                     interpolate_for_output_2d, &
                     interpolate_for_output_3d
  end interface interpolate_for_output
  !
  ! Module variables used for debugging
  !
  integer, save :: current_cell
  integer, save :: save_n1
  integer, save :: save_n2
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
      call check_input_file(input_fname,input_status)
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
          call check_input_file(input_fname,input_status)
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
    call check_input_file(input_fname,input_status)
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
  use ovar,   only : Default_Marker_c1
  !
  use quadrature_mod, only : init_geom_quadrature_rules
  use quadrature_mod, only : init_face_quadrature_rules
  !
  use restart_mod, only : check_restart_dir
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: input_fname
  integer,          intent(in) :: input_mode
  !
  !.. Local Scalars ..
  integer  :: i,it1,it2,ierr,ival,inpt
  integer  :: Period_Timesteps_X,Period_Timesteps_Y
  logical(lk)  :: use_qdtr_face_pts,mms_error
  real(wp) :: Period_Time_X,Period_Time_Y
  real(wp) :: real_Total_Periods
  real(wp) :: aoa_radians
  real(wp) :: cfl_sum,cfl_time
  real(wp) :: xmin,xmax,psi0
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
    solution_order           = interpto_solution_order
    all_orders               = interpto_all_orders
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
  p_order = projection_order
  q_order = quadrature_order
  o_order = output_order
  e_order = error_order
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
    if (mypnum == 0) write (iout,101)
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
  call init_geom_quadrature_rules
  call check_for_bad_number(__LINE__,__FILE__)
  call init_face_quadrature_rules
  call check_for_bad_number(__LINE__,__FILE__)
  !
  ! Change the upwind method to LDFSS if using gas mixture
  !
  if (ncs > 0) then
    if (invs_flux_method /= invs_flux_ROE  .and. &
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
 !if (mypnum == glb_root) then
 !  psi0 = -one*(marker_c1+marker_c2*log10(real(n_order,kind=wp)))
 !  write (iout,582) "BEFORE", &
 !                   Limiter_Marker_c1, &
 !                   Limiter_Marker_c2, &
 !                   Limiter_Marker_dpsi, &
 !                   Limiter_Marker_Min, &
 !                   Limiter_Marker_Max, &
 !                   marker_c1, &
 !                   marker_c2, &
 !                   marker_dpsi, &
 !                   psi0, &
 !                   psi0-marker_dpsi, &
 !                   psi0+marker_dpsi
 !end if
 !flush (iout)
 !582 format (//,2x,a," modifying marker parameters:",//, &
 !               4x,"Limiter_Marker_c1   = ",f10.6,/, &
 !               4x,"Limiter_Marker_c2   = ",f10.6,/, &
 !               4x,"Limiter_Marker_dpsi = ",f10.6,/, &
 !               4x,"Limiter_Marker_Min  = ",f10.6,/, &
 !               4x,"Limiter_Marker_Max  = ",f10.6,//, &
 !               4x,"        marker_c1   = ",f10.6,/, &
 !               4x,"        marker_c2   = ",f10.6,/, &
 !               4x,"        marker_dpsi = ",f10.6,/, &
 !               4x,"        marker_psi0 = ",f10.6,/, &
 !               4x,"        marker_min  = ",f10.6,/, &
 !               4x,"        marker_max  = ",f10.6,//)
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
      marker_c2 = -one*(Limiter_Marker_Min + marker_dpsi + marker_c1) / &
                   log10(real(n_order,kind=wp))
    else
      marker_c1 = -one*(Limiter_Marker_Min + &
                        marker_dpsi + marker_c2*log10(real(n_order,kind=wp)))
    end if
    !
  end if
 !!
 !if (mypnum == glb_root) then
 !  psi0 = -one*(marker_c1+marker_c2*log10(real(n_order,kind=wp)))
 !  write (iout,582) " AFTER", &
 !                   Limiter_Marker_c1, &
 !                   Limiter_Marker_c2, &
 !                   Limiter_Marker_dpsi, &
 !                   Limiter_Marker_Min, &
 !                   Limiter_Marker_Max, &
 !                   marker_c1, &
 !                   marker_c2, &
 !                   marker_dpsi, &
 !                   psi0, &
 !                   psi0-marker_dpsi, &
 !                   psi0+marker_dpsi
 !end if
 !flush (iout)
 !call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,"Dumping marker parameters!")
  !
!!if (mypnum == glb_root) then
!!  psi0 = -one*(marker_c1 + marker_c2*log10(real(n_order,kind=wp)))
!!  write (iout,30) marker_c1,marker_c2,marker_dpsi, &
!!                  psi0-marker_dpsi,psi0,psi0+marker_dpsi
!!end if
!!30 format (/,2x,"marker_c1   = ",f11.6,/, &
!!             2x,"marker_c2   = ",f11.6,/, &
!!             2x,"marker_dpsi = ",f11.6,//, &
!!             2x,"psi0 - dpsi = ",f11.6,/, &
!!             2x,"psi0        = ",f11.6,/, &
!!             2x,"psi0 + dpsi = ",f11.6,/)
  !
  ! Make sure the number of RK stages is correct for some RK types
  !
  if (Runge_Kutta_Scheme == TVD2_RK) then
    num_rk_stages = 2
  else if (Runge_Kutta_Scheme == TVD3_RK) then
    num_rk_stages = 3
  else if (Runge_Kutta_Scheme == TVD4_RK) then
    num_rk_stages = 5
  else if (Runge_Kutta_Scheme == CK4_RK) then
    num_rk_stages = 5
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
     !if (mypnum == 0) then
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
     !if (mypnum == 0) then
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
     !if (mypnum == 0) then
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
        if (mypnum == 0) write (iout,103) real_Total_Periods
      end if
    end if
   !if (mypnum == 0) write (*,*) Period_Time,Total_Periods,Final_Time
    !
    num_timesteps = nint( Final_Time / constant_dt )
    !
  else if (itestcase == Taylor_Green_Vortex) then
    !
   !cfl_beg = CFL
   !cfl_end = CFL
   !cfl_cycles = 1
    !
    ! Set the variable time_ref to adjust our nondimension time to match
    ! the nondimensional time used for the Taylor-Green Vortex problem
    !
    time_ref = machref
    !
    if (Final_Time > zero) then
      !
      ! Need to adjust the nondimensional final time which uses uref as the
      ! velocity for nondimensionalization and we use aref within the code
      !
      num_timesteps = nint( Final_Time / constant_dt )
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
    ! Set time_ref to Lref/aref for nondimensionalizing time
    ! NOTE: Lref is always just one
    !
    if (itestcase == Infinite_Cylinder) then
      time_ref = machref
    else
      time_ref = one / aref
    end if
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
    if (itestcase == Taylor_Green_Vortex) then
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
      " ************************** WARNING!!!! *************************",/)
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
  21 format ("MMS ERROR: MMS solutions 1 and 2 are allowed for 2D and 3D only!")
  22 format ("MMS ERROR: MMS solution 3 is allowed for 1D and 2D only!")
  23 format ("MMS ERROR: MMS solution 4 is allowed for 2D only!")
  24 format ("MMS ERROR: The MMS solution chosen is invalid!")
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
  call allocate_namelist_variables
  !
  if (mypnum == 0) then
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
    if (mypnum == 0) write (iout,1)
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
  if (mypnum == 0) then
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
  root_proc_only: if (mypnum == 0) then
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
    if (mypnum == 0) write (iout,1) trim(adjustl(output_dir))
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
      if (mypnum == 0) then
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
  use ovar, only : output_dir
  use ovar, only : lustre_stripe_count
  use ovar, only : lustre_stripe_size
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
  if (mypnum == 0) then
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
  use geovar,         only : nr,grid
  use module_cgns,    only : read_cgns_gridfile
  use module_gmsh,    only : read_gmsh_gridfile
  use module_plot3d,  only : read_plot3d_gridfile
  use quadrature_mod, only : find_valid_geometries
  use ovar,           only : use_bc_conflicts
  use ovar,           only : itestcase
  use ovar,           only : postproc_debug_grid
  use ovar,           only : grid_format
  use ovar,           only : gridfile
  use ovar,           only : plot3d_cutfile
  use ovar,           only : plot3d_bcsfile
  use ovar,           only : solution_file
  use ovar,           only : flame_speed_test_dof
  !
  use initialization_mod, only : initialize_reference
  use initialization_mod, only : create_bc_in
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
  if (postproc_debug_grid) then
    call create_postproc_debug_grid
  else if (itestcase == Taylor_Green_Vortex .and. &
           grid_format == Internally_Generated_Grid) then
    call create_TaylorGreenVortex_grid
  else if (any(abs(itestcase) == Transport_Problems) .and. &
           grid_format == Internally_Generated_Grid) then
    call create_IsentropicVortex_grid
  else if (abs(itestcase) == CH4_Flame_Speed_Test .and. &
           flame_speed_test_dof > 0) then
    call create_flame_speed_test_grid
  else if (abs(itestcase) == Channel_Flow .and. &
           grid_format == Internally_Generated_Grid) then
    call create_channel_flow_grid
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
  if (.not. allocated(grid)) then
    if (i_am_host_root) then
      call create_linear_grid_dt
    end if
  end if
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
  ! Alter bface due to the input boundary conditions
  ! NOTE: This will apply the correct BC to bface(1,:) and fill bface(6,:)
  !
  call apply_bc_input_to_bface
  !
  ! Before doing anything with the boundary conditions, check to see if
  ! there are any SEM boundary conditions. If there are any, read in the
  ! SEM profiles and separate any boundary faces outside of the profile
  ! height into new subsonic inflow boundaries.
  !
  call separate_sem_boundary_faces
  !
  ! First, initialize the bc_in array
  ! NOTE: This will also fill bface(6,:) with the correct value
  !
  call create_bc_in
  !
  ! Reorder the boundary faces
  !
  call reorder_boundary_faces
  call check_for_bad_number(__LINE__,__FILE__)
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
subroutine create_flame_speed_test_grid()
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use ovar,      only : Lref,itestcase
  use ovar,      only : flame_speed_test_dof
  use geovar,    only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geovar,    only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geovar,    only : cell_geom,cell_order
  use ovar,      only : grid_format,bc_names
  !
  !.. Local Scalars ..
  integer  :: n,nc,nf,np,ierr
  real(wp) :: dl
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_flame_speed_test_grid"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Reset the grid format input parameter to identify
  ! the grid as internally generated
  !
  grid_format = Internally_Generated_Grid
  !
  ! Make the reference length one centimeter for the flame speed test case
  !
  Lref = one/(ten*ten)
  !
  ! The flame speed test case is always 1D
  !
  nr = 1
  !
  ! Estimate the number of cells required to get the requested DOF
  !
  nc = nint( real(flame_speed_test_dof,kind=wp) / real(n_order+1,kind=wp) )
  !
  np = nc + 1 ! number of grid nodes in each direction
  !
  ncell = nc    ! total number of grid cells
  nnode = np    ! total number of grid nodes
  nfbnd = 2     ! total number of boundary faces
  !
  ! Allocate all the arrays needed to define the grid
  !
  allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_nodes",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  allocate ( nodes_of_cell(1:geom_nodes(Geom_Edge)*ncell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_geom(1:ncell) , source=Geom_Edge , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_geom",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_order(1:ncell) , source=n_order , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_order",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! #######
  ! #######
  ! ## 1 ##
  ! #######
  ! #######
  !
  ! Create the nodes_of_cell_ptr array
  !
  nodes_of_cell_ptr(1) = 0
  do n = 2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(n) = nodes_of_cell_ptr(n-1) + geom_nodes(Geom_Edge)
  end do
  !
  ! Create the nodes_of_cell array
  !
  do n = 1,ncell
    nodes_of_cell( nodes_of_cell_ptr(n)+1 ) = n
    nodes_of_cell( nodes_of_cell_ptr(n)+2 ) = n+1
  end do
  !
  ! #######
  ! #######
  ! ## 2 ##
  ! #######
  ! #######
  !
  ! Compute the cell size for a Cartesian/uniform grid of total length Lref
  !
  dl = Lref / real(ncell,kind=wp)
  !
  ! Create the grid node coordinates as a uniform 1D grid from -Lref/2 to Lref/2
  !
  do n = 1,nnode
    xyz_nodes(1,n) = -half*Lref + dl*real(n-1,kind=wp)
  end do
  !
  ! #######
  ! #######
  ! ## 3 ##
  ! #######
  ! #######
  !
  ! Create the bface array
  !
  ! Since this is a 1D grid, there are only 2 boundary faces: i-min and i-max
  !
  ! Create the i-min boundary face
  !
  nf = 1
  !
  bface( 1,nf) = bc_sub_inflow
  bface( 2,nf) = 1      ! host cell index
  bface( 3,nf) = Geom_Edge
  bface( 4,nf) = 1      ! x-min face
  bface( 5,nf) = 1      ! this face belongs to boundary group 1
  bface( 6,nf) = 1      ! index of array bc_in(:) to get BC flow conditions
  bface( 7,nf) = 0      ! not used for physical boundary
  bface( 8,nf) = 0      ! not used for physical boundary
  bface( 9,nf) = 0      ! not used for physical boundary
  bface(10,nf) = ncell + nf ! index of ghost cell/face
  !
  ! Create the i-max boundary face
  !
  nf = 2
  !
  bface( 1,nf) = bc_sub_outflow
  bface( 2,nf) = ncell  ! host cell index
  bface( 3,nf) = Geom_Edge
  bface( 4,nf) = 2      ! x-max face
  bface( 5,nf) = 2      ! this face belongs to boundary group 2
  bface( 6,nf) = 2      ! index of array bc_in(:) to get BC flow conditions
  bface( 7,nf) = 0      ! not used for physical boundary
  bface( 8,nf) = 0      ! not used for physical boundary
  bface( 9,nf) = 0      ! not used for physical boundary
  bface(10,nf) = ncell + nf ! index of ghost cell/face
  !
  ! #######
  ! #######
  ! ## 3 ##
  ! #######
  ! #######
  !
  ! Create the bc_names array
  !
  allocate ( bc_names(1:2) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bc_names",1,__LINE__,__FILE__,ierr,error_message)
  !
  bc_names(1) = "INFLOW"
  bc_names(2) = "OUTFLOW"
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine create_flame_speed_test_grid
!
!###############################################################################
!
subroutine create_IsentropicVortex_grid()
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use ovar,      only : Lref,itestcase
  use ovar,      only : VortexGrid_DOF,VortexGrid_Lref
  use ovar,      only : VortexGrid_random_perturb
  use ovar,      only : VortexGrid_perturb_percent
  use geovar,    only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geovar,    only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geovar,    only : cell_geom,cell_order
  use ovar,      only : grid_format
  !
  !.. Local Scalars ..
  integer :: nc,np,ierr,i,j,l,n
  integer :: nf_min,nf_max
  real(r4) :: x,y,dl
  real(r4) :: rbas,radj,rval
  !
  !.. Local Allocatable Arrays ..
  real(r4), allocatable :: perturb(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_IsentropicVortex_grid"
 !logical(lk), parameter :: test_random_bc = true
  logical(lk), parameter :: test_random_bc = fals
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Reset the grid format input parameter to identify
  ! the grid as internally generated
  !
  grid_format = Internally_Generated_Grid
  !
  ! If we are generating the grid on the fly, make the reference length 1
  !
  Lref = one
  !
  ! The isentropic Euler vortex problem is always 2D
  !
  nr = 2
  !
  ! Estimate the number of cells required to get the requested DOF
  !
  nc = nint( real(VortexGrid_DOF,kind=wp) / real(n_order+1,kind=wp) )
  !
  np = nc + 1 ! number of grid nodes in each direction
  !
  ncell = nc*nc ! total number of grid cells
  nnode = np*np ! total number of grid nodes
  nfbnd = 4*nc  ! total number of boundary faces
  !
  ! Allocate all the arrays needed to define the grid
  !
  allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_nodes",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  allocate ( nodes_of_cell(1:geom_nodes(Geom_Quad)*ncell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_geom(1:ncell) , source=Geom_Quad , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_geom",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_order(1:ncell) , source=n_order , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_order",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( perturb(1:nr,1:nnode) , source=real(0,kind=r4) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"perturb",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! #######
  ! #######
  ! ## 1 ##
  ! #######
  ! #######
  !
  ! Create the nodes_of_cell_ptr array
  !
  nodes_of_cell_ptr(1) = 0
  do n = 2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(n) = nodes_of_cell_ptr(n-1) + geom_nodes(Geom_Quad)
  end do
  !
  ! Create the nodes_of_cell array
  !
  do j = 1,nc
    do i = 1,nc
      !
      n = nodes_of_cell_ptr( map2dto1d(i,j,nc) )
      !
      nodes_of_cell(n+1) = map2dto1d(i  ,j  ,np)
      nodes_of_cell(n+2) = map2dto1d(i+1,j  ,np)
      nodes_of_cell(n+3) = map2dto1d(i+1,j+1,np)
      nodes_of_cell(n+4) = map2dto1d(i  ,j+1,np)
      !
    end do
  end do
  !
  ! #######
  ! #######
  ! ## 2 ##
  ! #######
  ! #######
  !
  ! Compute the cell length for a Cartesian/uniform grid
  !
  dl = real(1,kind=r4) / real(np-1,kind=r4)
  !
  ! Compute the random perturbations if the flag was set on input
  !
  if (VortexGrid_random_perturb) then
    !
    rbas = dl * real(abs(VortexGrid_perturb_percent),kind=r4) &
              / real(100,kind=r4)
    radj = rbas / real(2,kind=r4)
    !
    ! Randomly perturb the crap out of everything
    !
    do j = 2,np-1
      do i = 2,np-1
        !
        n = map2dto1d(i,j,np)
        !
        do l = 1,nr
          rval = random(rbas) - radj
          perturb(l,n) = rval
        end do
        !
      end do
    end do
    !
    ! Have the root processor broadcast its perturbations to all the
    ! other processors so that the perturbations will be consistent
    ! on all processors.
    !
    call mpi_bcast(perturb,size(perturb,kind=int_mpi), &
                   MPI_REAL,0_int_mpi,MPI_COMM_WORLD,mpierr)
    !
  end if
  !
  ! Create the grid node coordinates as a square from (0,0) to (1,1)
  !
  do j = 1,np
    !
    y = dl * real(j-1,kind=r4)
    !
    do i = 1,np
      !
      x = dl * real(i-1,kind=r4)
      !
      n = map2dto1d(i,j,np)
      !
      xyz_nodes(1,n) = real( x + perturb(1,n) , kind=wp )
      xyz_nodes(2,n) = real( y + perturb(2,n) , kind=wp )
      !
    end do
  end do
  !
  ! Deallocate perturb since it is no longer needed
  !
  deallocate ( perturb , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"perturb",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Scale the grid node coordinates to the square (-1,-1) to (1,1)
  !
  xyz_nodes = two*xyz_nodes - one
  !
  ! Scale the grid node coordinates to the square (-Lref,-Lref) to (+Lref,+Lref)
  !
  xyz_nodes = abs(VortexGrid_Lref) * xyz_nodes
  !
  ! #######
  ! #######
  ! ## 3 ##
  ! #######
  ! #######
  !
  ! Create the bface array
  !
  only_root_creates_bface: if (mypnum == 0) then
    !
    nf_max = 0
    !
    ! Loop through the 4 grid boundaries
    ! NOTE: This corresponds to the grid boundaries in
    !       the order [imin,imax,jmin,jmax]
    !
    nf_min = nf_max
    !
    do j = 1,nc
      !
      nf_min = nf_min + 1
      nf_max = nf_min + nc
      !
      bface( 1,nf_min) = bc_periodic
      bface( 2,nf_min) = map2dto1d(1 ,j,nc)
      bface( 3,nf_min) = Geom_Quad
      bface( 4,nf_min) = HostFace2D(1)
      bface( 5,nf_min) = nf_max
      bface( 6,nf_min) = map2dto1d(nc,j,nc)
      bface( 8,nf_min) = Geom_Quad
      bface( 9,nf_min) = HostFace2D(2)
      bface(10,nf_min) = ncell + nf_min
      !
      bface( 1,nf_max) = bc_periodic
      bface( 2,nf_max) = bface(6,nf_min)
      bface( 3,nf_max) = bface(8,nf_min)
      bface( 4,nf_max) = bface(9,nf_min)
      bface( 5,nf_max) = nf_min
      bface( 6,nf_max) = bface(2,nf_min)
      bface( 8,nf_max) = bface(3,nf_min)
      bface( 9,nf_max) = bface(4,nf_min)
      bface(10,nf_max) = ncell + nf_max
      !
      ! Adjust boundary conditions for stationary vortex or
      ! vortex propagation in only the vertical direction
      !
      if (any(itestcase == [-Shu_Vortex,-Vortex_Transport])) then
        !
        if (test_random_bc) then
          l = size(bc_heirarchy)
          l = max( 1 , min(l,random(l)) )
          bface(1,nf_min) = bc_heirarchy(l)
        else
          bface(1,nf_min) = bc_characteristic
        end if
        bface(  5,nf_min) = 1
        bface(6:9,nf_min) = 0
        !
        if (test_random_bc) then
          l = size(bc_heirarchy)
          l = max( 1 , min(l,random(l)) )
          bface(1,nf_max) = bc_heirarchy(l)
        else
          bface(1,nf_max) = bc_characteristic
        end if
        bface(  5,nf_max) = 1
        bface(6:9,nf_max) = 0
        !
      end if
      !
    end do
    !
    nf_min = nf_max
    !
    do i = 1,nc
      !
      nf_min = nf_min + 1
      nf_max = nf_min + nc
      !
      bface( 1,nf_min) = bc_periodic
      bface( 2,nf_min) = map2dto1d(i,1 ,nc)
      bface( 3,nf_min) = Geom_Quad
      bface( 4,nf_min) = HostFace2D(3)
      bface( 5,nf_min) = nf_max
      bface( 6,nf_min) = map2dto1d(i,nc,nc)
      bface( 8,nf_min) = Geom_Quad
      bface( 9,nf_min) = HostFace2D(4)
      bface(10,nf_min) = ncell + nf_min
      !
      bface( 1,nf_max) = bc_periodic
      bface( 2,nf_max) = bface(6,nf_min)
      bface( 3,nf_max) = bface(8,nf_min)
      bface( 4,nf_max) = bface(9,nf_min)
      bface( 5,nf_max) = nf_min
      bface( 6,nf_max) = bface(2,nf_min)
      bface( 8,nf_max) = bface(3,nf_min)
      bface( 9,nf_max) = bface(4,nf_min)
      bface(10,nf_max) = ncell + nf_max
      !
      ! Adjust boundary conditions for stationary vortex
      !
      if (itestcase == -Shu_Vortex) then
        !
        if (test_random_bc) then
          l = size(bc_heirarchy)
          l = max( 1 , min(l,random(l)) )
          bface(1,nf_min) = bc_heirarchy(l)
        else
          bface(1,nf_min) = bc_characteristic
        end if
        bface(  5,nf_min) = 1
        bface(6:9,nf_min) = 0
        !
        if (test_random_bc) then
          l = size(bc_heirarchy)
          l = max( 1 , min(l,random(l)) )
          bface(1,nf_max) = bc_heirarchy(l)
        else
          bface(1,nf_max) = bc_characteristic
        end if
        bface(  5,nf_max) = 1
        bface(6:9,nf_max) = 0
        !
      end if
      !
    end do
    !
  end if only_root_creates_bface
  !
  if (ncpu > 1) then
    call mpi_bcast(bface,size(bface,kind=int_mpi), &
                   mpi_inttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  contains
  !
  pure function map2dto1d(i,j,n) result(return_value)
    integer, intent(in) :: i,j,n
    integer :: return_value
    return_value = i + n*(j-1)
  end function map2dto1d
  !
end subroutine create_IsentropicVortex_grid
!
!###############################################################################
!
subroutine create_TaylorGreenVortex_grid()
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use ovar,      only : VortexGrid_DOF,Lref
  use ovar,      only : VortexGrid_random_perturb
  use ovar,      only : VortexGrid_perturb_percent
  use geovar,    only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geovar,    only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geovar,    only : cell_geom,cell_order
  use ovar,      only : grid_format
  !
  !.. Local Scalars ..
  integer :: nc,np,ierr,i,j,k,l,n,n1,n2
  integer :: nf_min,nf_max
  real(r4) :: x,y,z,dl
  real(r4) :: rbas,radj,rval
  !
  !.. Local Allocatable Arrays ..
  real(r4), allocatable :: perturb(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_TaylorGreenVortex_grid"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Reset the grid format input parameter to identify
  ! the grid as internally generated
  !
  grid_format = Internally_Generated_Grid
  !
  ! If we are generating the grid on the fly, make the reference length 1
  !
  Lref = one
  !
  ! The Taylor-Green vortex problem is always 3D
  !
  nr = 3
  !
  ! Estimate the number of cells required to get the requested DOF
  !
  nc = nint( real(VortexGrid_DOF,kind=wp) / real(n_order+1,kind=wp) )
  !
  np = nc + 1 ! number of grid nodes in each direction
  !
  ncell = nc*nc*nc ! total number of grid cells
  nnode = np*np*np ! total number of grid nodes
  nfbnd = 6*nc*nc  ! total number of boundary faces
  !
  ! Allocate all the arrays needed to define the grid
  !
  allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_nodes",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  allocate ( nodes_of_cell(1:geom_nodes(Geom_Hexa)*ncell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_geom(1:ncell) , source=Geom_Hexa , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_geom",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_order(1:ncell) , source=n_order , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_order",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( perturb(1:nr,1:nnode) , source=real(0,kind=r4) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"perturb",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! #######
  ! #######
  ! ## 1 ##
  ! #######
  ! #######
  !
  ! Create the nodes_of_cell_ptr array
  !
  nodes_of_cell_ptr(1) = 0
  do n = 2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(n) = nodes_of_cell_ptr(n-1) + geom_nodes(Geom_Hexa)
  end do
  !
  ! Create the nodes_of_cell array
  !
  do k = 1,nc
    do j = 1,nc
      do i = 1,nc
        !
        n = nodes_of_cell_ptr( map3dto1d(i,j,k,nc) )
        !
        nodes_of_cell(n+1) = map3dto1d(i  ,j  ,k  ,np)
        nodes_of_cell(n+2) = map3dto1d(i+1,j  ,k  ,np)
        nodes_of_cell(n+3) = map3dto1d(i+1,j+1,k  ,np)
        nodes_of_cell(n+4) = map3dto1d(i  ,j+1,k  ,np)
        nodes_of_cell(n+5) = map3dto1d(i  ,j  ,k+1,np)
        nodes_of_cell(n+6) = map3dto1d(i+1,j  ,k+1,np)
        nodes_of_cell(n+7) = map3dto1d(i+1,j+1,k+1,np)
        nodes_of_cell(n+8) = map3dto1d(i  ,j+1,k+1,np)
        !
      end do
    end do
  end do
  !
  ! #######
  ! #######
  ! ## 2 ##
  ! #######
  ! #######
  !
  ! Compute the cell length for a Cartesian/uniform grid
  !
  dl = real(1,kind=r4) / real(np-1,kind=r4)
  !
  ! Compute the random perturbations if the flag was set on input
  !
  if (VortexGrid_random_perturb) then
    !
    rbas = dl * real(abs(VortexGrid_perturb_percent),kind=r4) &
              / real(100,kind=r4)
    radj = rbas / real(2,kind=r4)
    !
    if (VortexGrid_perturb_percent < 0) then
      !
      ! Randomly perturb the crap out of everything
      !
      do k = 2,np-1
        do j = 2,np-1
          do i = 2,np-1
            !
            n = map3dto1d(i,j,k,np)
            !
            do l = 1,nr
              rval = random(rbas) - radj
              perturb(l,n) = rval
            end do
            !
          end do
        end do
      end do
      !
    else
      !
      ! Get the random perturbation in along i for y
      !
      do i = 2,np-1
        rval = random(rbas) - radj
        do k = 2,np-1
          do j = 2,np-1
            n = map3dto1d(i,j,k,np)
            perturb(2,n) = rval
          end do
        end do
      end do
      !
      ! Get the random perturbation along j for z
      !
      do j = 2,np-1
        rval = random(rbas) - radj
        do k = 2,np-1
          do i = 2,np-1
            n = map3dto1d(i,j,k,np)
            perturb(3,n) = rval
          end do
        end do
      end do
      !
      ! Get the random perturbation along k for x
      !
      do k = 2,np-1
        rval = random(rbas) - radj
        do j = 2,np-1
          do i = 2,np-1
            n = map3dto1d(i,j,k,np)
            perturb(1,n) = rval
          end do
        end do
      end do
      !
    end if
    !
    ! Have the root processor broadcast its perturbations to all the
    ! other processors so that the perturbations will be consistent
    ! on all processors.
    !
    call mpi_bcast(perturb,size(perturb,kind=int_mpi), &
                   MPI_REAL,0_int_mpi,MPI_COMM_WORLD,mpierr)
    !
  end if
  !
  ! Create the grid node coordinates as a cube from (0,0,0) to (1,1,1)
  !
  do k = 1,np
    !
    z = dl * real(k-1,kind=r4)
    !
    do j = 1,np
      !
      y = dl * real(j-1,kind=r4)
      !
      do i = 1,np
        !
        x = dl * real(i-1,kind=r4)
        !
        n = map3dto1d(i,j,k,np)
        !
        xyz_nodes(1,n) = real( x + perturb(1,n) , kind=wp )
        xyz_nodes(2,n) = real( y + perturb(2,n) , kind=wp )
        xyz_nodes(3,n) = real( z + perturb(3,n) , kind=wp )
        !
      end do
    end do
  end do
  !
  ! Deallocate perturb since it is no longer needed
  !
  deallocate ( perturb , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"perturb",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Scale the grid node coordinates to the cube (-1,-1,-1) to (1,1,1)
  !
  xyz_nodes = two*xyz_nodes - one
  !
  ! Scale the grid node coordinates to the cube (-pi,-pi,-pi) to (pi,pi,pi)
  !
  xyz_nodes = PI * xyz_nodes
  !
  ! #######
  ! #######
  ! ## 3 ##
  ! #######
  ! #######
  !
  ! Create the bface array
  !
  nf_max = 0
  !
  ! Loop through the 6 grid boundaries
  ! NOTE: This corresponds to the grid boundaries in
  !       the order [imin,imax,jmin,jmax,kmin,kmax]
  !
  nf_min = nf_max
  !
  do k = 1,nc
    do j = 1,nc
      !
      n1 = nf_min + 1
      n2 = nf_min + 1 + nc*nc
      !
      bface( 1,n1) = bc_periodic
      bface( 2,n1) = map3dto1d(1 ,j,k,nc)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(1)
      bface( 5,n1) = n2
      bface( 6,n1) = map3dto1d(nc,j,k,nc)
      bface( 8,n1) = Geom_Hexa
      bface( 9,n1) = HostFace3D(2)
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_periodic
      bface( 2,n2) = bface(6,n1)
      bface( 3,n2) = bface(8,n1)
      bface( 4,n2) = bface(9,n1)
      bface( 5,n2) = n1
      bface( 6,n2) = bface(2,n1)
      bface( 8,n2) = bface(3,n1)
      bface( 9,n2) = bface(4,n1)
      bface(10,n2) = ncell + n2
      !
      nf_min = n1
      nf_max = n2
      !
    end do
  end do
  !
  nf_min = nf_max
  !
  do k = 1,nc
    do i = 1,nc
      !
      n1 = nf_min + 1
      n2 = nf_min + 1 + nc*nc
      !
      bface( 1,n1) = bc_periodic
      bface( 2,n1) = map3dto1d(i,1 ,k,nc)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(3)
      bface( 5,n1) = n2
      bface( 6,n1) = map3dto1d(i,nc,k,nc)
      bface( 8,n1) = Geom_Hexa
      bface( 9,n1) = HostFace3D(4)
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_periodic
      bface( 2,n2) = bface(6,n1)
      bface( 3,n2) = bface(8,n1)
      bface( 4,n2) = bface(9,n1)
      bface( 5,n2) = n1
      bface( 6,n2) = bface(2,n1)
      bface( 8,n2) = bface(3,n1)
      bface( 9,n2) = bface(4,n1)
      bface(10,n2) = ncell + n2
      !
      nf_min = n1
      nf_max = n2
      !
    end do
  end do
  !
  nf_min = nf_max
  !
  do j = 1,nc
    do i = 1,nc
      !
      n1 = nf_min + 1
      n2 = nf_min + 1 + nc*nc
      !
      bface( 1,n1) = bc_periodic
      bface( 2,n1) = map3dto1d(i,j,1 ,nc)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(5)
      bface( 5,n1) = n2
      bface( 6,n1) = map3dto1d(i,j,nc,nc)
      bface( 8,n1) = Geom_Hexa
      bface( 9,n1) = HostFace3D(6)
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_periodic
      bface( 2,n2) = bface(6,n1)
      bface( 3,n2) = bface(8,n1)
      bface( 4,n2) = bface(9,n1)
      bface( 5,n2) = n1
      bface( 6,n2) = bface(2,n1)
      bface( 8,n2) = bface(3,n1)
      bface( 9,n2) = bface(4,n1)
      bface(10,n2) = ncell + n2
      !
      nf_min = n1
      nf_max = n2
      !
    end do
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  contains
  !
  pure function map3dto1d(i,j,k,n) result(return_value)
    integer, intent(in) :: i,j,k,n
    integer :: return_value
    return_value = i + n*(j-1) + n*n*(k-1)
  end function map3dto1d
  !
end subroutine create_TaylorGreenVortex_grid
!
!###############################################################################
!
subroutine create_channel_flow_grid()
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use ovar,      only : Lref,loc_solution_pts,itestcase
  use ovar,      only : channel_grid_ave_xplus
  use ovar,      only : channel_grid_ave_yplus
  use ovar,      only : channel_grid_ave_zplus
  use ovar,      only : channel_grid_yplus_wall
  use ovar,      only : channel_yplus_of_one
  use geovar,    only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geovar,    only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geovar,    only : cell_geom,cell_order
  use ovar,      only : grid_format
  !
  use polynomial_mod, only : legendre_gauss_quadrature
  use polynomial_mod, only : legendre_gauss_lobatto_quadrature
  !
  !.. Local Scalars ..
  integer  :: ierr,i,j,k,l,n,n1,n2
  integer  :: nx,ny,nz,nf_min,nf_max
  real(wp) :: xplus_ave,yplus_ave,zplus_ave,dx,dy,dz
  real(wp) :: ds1,ds2,dl1,dl2,a,b,sf,rnum,rden,sferr,val,ue,se
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: xi(:)
  real(wp), allocatable :: x(:)
  real(wp), allocatable :: y(:)
  real(wp), allocatable :: z(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_channel_flow_grid"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Reset the grid format input parameter to identify
  ! the grid as internally generated
  !
  grid_format = Internally_Generated_Grid
  !
  ! If we are generating the grid on the fly, make the reference length 1
  !
  Lref = one
  !
  ! The channel flow problem is always 3D
  !
  nr = 3
  !
  ! #######
  ! #######
  ! ## 1 ##
  ! #######
  ! #######
  !
  ! Reset any of the grid parameters to their default values if they are zero
  !
  if (channel_yplus_of_one    == zero) channel_yplus_of_one = 8.42134e-6_wp
  if (channel_grid_ave_xplus  == zero) channel_grid_ave_xplus = five*ten
  if (channel_grid_ave_yplus  == zero) channel_grid_ave_yplus = five*five
  if (channel_grid_ave_zplus  == zero) channel_grid_ave_zplus = five*five
  if (channel_grid_yplus_wall == zero) channel_grid_yplus_wall = -one
  !
  ! Compute the cell average x+,y+,z+ for the grid
  ! For a given coordinate direction (say [i]), the cell average [i]+ is
  !   if (channel_grid_ave_[i]plus > zero) then
  !     [i]plus_ave = channel_grid_ave_[i]plus
  !   else
  !     [i]plus_ave = abs(channel_grid_ave_[i]plus)*real(n_order+2,kind=wp)
  !   end if
  ! NOTE: We use (P+2) here because the delta plus value is the
  !       average distance between the DoF within a cell (or the
  !       subcells), and there are P+2 subcells in each cell.
  !
  xplus_ave = abs(channel_grid_ave_xplus)*merge(one, &
                                                real(n_order+2,kind=wp), &
                                                channel_grid_ave_xplus > zero)
  yplus_ave = abs(channel_grid_ave_yplus)*merge(one, &
                                                real(n_order+2,kind=wp), &
                                                channel_grid_ave_yplus > zero)
  zplus_ave = abs(channel_grid_ave_zplus)*merge(one, &
                                                real(n_order+2,kind=wp), &
                                                channel_grid_ave_zplus > zero)
  !
  ! Get the number of grid points in each direction
  !
  nx = nint( channel_grid_Lx / (channel_yplus_of_one*xplus_ave) )
  ny = nint( channel_grid_Ly / (channel_yplus_of_one*yplus_ave) )
  nz = nint( channel_grid_Lz / (channel_yplus_of_one*zplus_ave) )
  !
  ! Compute the x and z cell spacings
  !
  dx = channel_grid_Lx / real(nx,kind=wp)
  dz = channel_grid_Lz / real(nz,kind=wp)
  !
  ! Initialize the x,y,z grid spacing arrays
  !
  x = init_array(nx+1,zero) ! F2003 AUTO-REALLOC
  y = init_array(ny+1,zero) ! F2003 AUTO-REALLOC
  z = init_array(nz+1,zero) ! F2003 AUTO-REALLOC
  !
  ! #######
  ! #######
  ! ## 2 ##
  ! #######
  ! #######
  !
  ! Get the grid cell spacings for x,y,z
  !
  ! ######
  ! # 2A #
  ! ######
  !
  ! Get the grid point locations for x
  !
  do n = 2,nx
    x(n) = x(n-1) + dx
  end do
  !
  ! ######
  ! # 2B #
  ! ######
  !
  ! Get the grid point locations for z
  !
  do n = 2,nz
    z(n) = z(n-1) + dz
  end do
  !
  ! ######
  ! # 2C #
  ! ######
  !
  ! 2C.1 : Get the beginning and ending spacings needed
  !        to generate the tanh spacing
  !
  ! Need to get dl1/ds1 and dl2/ds2 to create tanh spacing for y spacing
  !
  if (channel_grid_yplus_wall < zero) then
    !
    ! The i-th point of the first cell is at y+=1
    !
    i = min( n_order+1 , max( 1 , nint(abs(channel_grid_yplus_wall)) ) )
    !
    ! Get the location of the 1D solution points in reference space
    !
    xi = init_array(n_order+1,zero) ! F2003 AUTO-REALLOC
    !
    if (loc_solution_pts == Legendre_Gauss_Lobatto) then
      i = min(i+1,n_order+1) ! Offset for 1st Lobatto point being at the wall
      call legendre_gauss_lobatto_quadrature(xi)
    else
      call legendre_gauss_quadrature(xi)
    end if
    !
    xi(:) = half*(xi(:) + one)
    !
    dl1 = channel_yplus_of_one / xi(i)
    !
  else
    !
    ! The height of the first cell off the wall is channel_grid_yplus_wall
    !
    dl1 = channel_yplus_of_one * channel_grid_yplus_wall
    !
  end if
  !
  dl2 = dl1
  !
  ! Map the beginning and ending grid cell
  ! spacings to the reference interval [0,1].
  !
  ds1 = dl1 / channel_grid_Ly
  ds2 = dl2 / channel_grid_Ly
  !
  ! 2C.2 : Get the a and b parameters to begin solving for the tanh y-spacing
  !
  a = sqrt(ds2/ds1)
  b = one / (real(ny,kind=wp)*sqrt(ds1*ds2))
  !
  ! 2C.3 : Use b to solve the stretching factor function, then use the
  !        stretching factor to get the y-spacing.
  !
  if (b >= 1.1_wp) then
    !
    ! Use a Newton iteration to solve the stretching factor function
    !    sinh(sf)/sf = b
    ! for the stretching factor, sf.
    !
    sf = two ! initial guess for sf
    !
    newton_tanh_loop: &
    do i = 1,100
      !
      rnum = sinh(sf)/sf - b
      rden = cosh(sf)/sf - sinh(sf)/(sf*sf)
      !
      sferr = rnum / rden
      !
      if (abs(sferr) < eps6) then
        exit newton_tanh_loop
      end if
      !
      sf = sf - sferr
      !
      ! Make sure sf doesnt become 0
      !
      if (abs(sf) < 0.01_wp) then
        sf = 0.01_wp
      end if
      !
    end do newton_tanh_loop
    !
    ! The stretching factor function has two solutions for sf given any
    ! value of b (ignoring sf=0). Since the stretching factor function is
    ! an even function, the two possible solutions are +sf and -sf. We only
    ! care about +sf so take the absolute value of whatever result we get
    ! from the Newton iteration.
    !
    sf = abs(sf)
    !
    ! Get the distribution y(n)
    !
    do n = 0,ny
      val = real(n,kind=wp)/real(ny,kind=wp)
      ue = half*(one + tanh(sf*(val-half))/tanh(sf*half))
      se = ue / (a + (one-a)*ue)
      y(n) = se * channel_grid_Ly
    end do
    !
  else if (b <= 0.9_wp) then
    !
    ! Use a Newton iteration to solve the stretching factor function
    !    sin(sf)/sf = b
    ! for the stretching factor, sf.
    !
    sf = two ! initial guess for sf
    !
    newton_tan_loop: &
    do i = 1,100
      !
      rnum = sin(sf)/sf - b
      rden = cos(sf)/sf - sin(sf)/(sf*sf)
      !
      sferr = rnum / rden
      !
      if (abs(sferr) < eps6) then
        exit newton_tan_loop
      end if
      !
      sf = sf - sferr
      !
      ! Make sure sf doesnt become 0
      !
      if (abs(sf) < 0.01_wp) then
        sf = 0.01_wp
      end if
      !
    end do newton_tan_loop
    !
    ! The stretching factor function has two solutions for sf given any
    ! value of b (ignoring sf=0). Since the stretching factor function is
    ! an even function, the two possible solutions are +sf and -sf. We only
    ! care about +sf so take the absolute value of whatever result we get
    ! from the Newton iteration.
    !
    sf = abs(sf)
    !
    ! Get the distribution y(n)
    !
    do n = 0,ny
      val = real(n,kind=wp)/real(ny,kind=wp)
      ue = half*(one + tan(sf*(val-half))/tan(sf*half))
      se = ue / (a + (one-a)*ue)
      y(n) = se * channel_grid_Ly
    end do
    !
  else
    !
    ! Both of the above Newton iterations can have numerical difficulties
    ! if the stretching factor is close to 1, therefore use this simpler
    ! stretching function that doesnt require a Newton iteration.
    !
    do n = 0,ny
      val = real(n,kind=wp)/real(ny,kind=wp)
      ue = val*(one + two*(b-one)*(val-half)*(one-val))
      se = ue / (a + (one-a)*ue)
      y(n) = se * channel_grid_Ly
    end do
    !
  end if
  !
  ! #######
  ! #######
  ! ## 3 ##
  ! #######
  ! #######
  !
  ! We can finish creating the grid now that we have the x,y,z grid spacings
  !
  ncell = nx*ny*nz                    ! total number of grid cells
  nnode = (nx+1)*(ny+1)*(nz+1)        ! total number of grid nodes
  nfbnd = 2*nx*ny + 2*nx*nz + 2*ny*nz ! total number of boundary faces
  !
  ! Allocate all the arrays needed to define the grid
  !
  allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_nodes",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  allocate ( nodes_of_cell(1:geom_nodes(Geom_Hexa)*ncell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_geom(1:ncell) , source=Geom_Hexa , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_geom",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_order(1:ncell) , source=n_order , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_order",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! #######
  ! #######
  ! ## 4 ##
  ! #######
  ! #######
  !
  ! Create the nodes_of_cell_ptr array
  !
  nodes_of_cell_ptr(1) = 0
  do n = 2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(n) = nodes_of_cell_ptr(n-1) + geom_nodes(Geom_Hexa)
  end do
  !
  ! Create the nodes_of_cell array
  !
  do k = 1,nz
    do j = 1,ny
      do i = 1,nx
        !
        n = nodes_of_cell_ptr( map3dto1d(i,j,k,nx,ny) )
        !
        nodes_of_cell(n+1) = map3dto1d(i  ,j  ,k  ,nx+1,ny+1)
        nodes_of_cell(n+2) = map3dto1d(i+1,j  ,k  ,nx+1,ny+1)
        nodes_of_cell(n+3) = map3dto1d(i+1,j+1,k  ,nx+1,ny+1)
        nodes_of_cell(n+4) = map3dto1d(i  ,j+1,k  ,nx+1,ny+1)
        nodes_of_cell(n+5) = map3dto1d(i  ,j  ,k+1,nx+1,ny+1)
        nodes_of_cell(n+6) = map3dto1d(i+1,j  ,k+1,nx+1,ny+1)
        nodes_of_cell(n+7) = map3dto1d(i+1,j+1,k+1,nx+1,ny+1)
        nodes_of_cell(n+8) = map3dto1d(i  ,j+1,k+1,nx+1,ny+1)
        !
      end do
    end do
  end do
  !
  ! #######
  ! #######
  ! ## 5 ##
  ! #######
  ! #######
  !
  ! Create the grid node coordinates
  !
  do k = 1,nz+1
    do j = 1,ny+1
      do i = 1,nx+1
        !
        n = map3dto1d(i,j,k,nx+1,ny+1)
        !
        xyz_nodes(1,n) = x(i)
        xyz_nodes(2,n) = y(j)
        xyz_nodes(3,n) = z(k)
        !
      end do
    end do
  end do
  !
  ! #######
  ! #######
  ! ## 6 ##
  ! #######
  ! #######
  !
  ! Create the bface array
  !
  nf_max = 0
  !
  ! Loop through the 6 grid boundaries
  ! NOTE: This corresponds to the grid boundaries in
  !       the order [imin,imax,jmin,jmax,kmin,kmax]
  !
  nf_min = nf_max
  !
  do k = 1,nz
    do j = 1,ny
      !
      n1 = nf_min + 1
      n2 = nf_min + 1 + ny*nz
      !
      if (itestcase == -Channel_Flow) then
        !
        bface( 1,n1) = bc_sem
        bface( 2,n1) = map3dto1d(1 ,j,k,nx,ny)
        bface( 3,n1) = Geom_Hexa
        bface( 4,n1) = HostFace3D(1)
        bface( 5,n1) = 1
        bface( 6,n1) = 1
        bface(10,n1) = ncell + n1
        !
        bface( 1,n2) = bc_sub_outflow
        bface( 2,n2) = map3dto1d(nx,j,k,nx,ny)
        bface( 3,n2) = Geom_Hexa
        bface( 4,n2) = HostFace3D(2)
        bface( 5,n2) = 2
        bface( 6,n2) = 2
        bface(10,n2) = ncell + n2
        !
      else
        !
        bface( 1,n1) = bc_periodic
        bface( 2,n1) = map3dto1d(1 ,j,k,nx,ny)
        bface( 3,n1) = Geom_Hexa
        bface( 4,n1) = HostFace3D(1)
        bface( 5,n1) = n2
        bface( 6,n1) = map3dto1d(nx,j,k,nx,ny)
        bface( 8,n1) = Geom_Hexa
        bface( 9,n1) = HostFace3D(2)
        bface(10,n1) = ncell + n1
        !
        bface( 1,n2) = bc_periodic
        bface( 2,n2) = bface(6,n1)
        bface( 3,n2) = bface(8,n1)
        bface( 4,n2) = bface(9,n1)
        bface( 5,n2) = n1
        bface( 6,n2) = bface(2,n1)
        bface( 8,n2) = bface(3,n1)
        bface( 9,n2) = bface(4,n1)
        bface(10,n2) = ncell + n2
        !
      end if
      !
      nf_min = n1
      nf_max = n2
      !
    end do
  end do
  !
  nf_min = nf_max
  !
  do k = 1,nz
    do i = 1,nx
      !
      n1 = nf_min + 1
      n2 = nf_min + 1 + nx*nz
      !
      bface( 1,n1) = bc_periodic
      bface( 2,n1) = map3dto1d(i,1 ,k,nx,ny)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(3)
      bface( 5,n1) = n2
      bface( 6,n1) = map3dto1d(i,ny,k,nx,ny)
      bface( 8,n1) = Geom_Hexa
      bface( 9,n1) = HostFace3D(4)
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_periodic
      bface( 2,n2) = bface(6,n1)
      bface( 3,n2) = bface(8,n1)
      bface( 4,n2) = bface(9,n1)
      bface( 5,n2) = n1
      bface( 6,n2) = bface(2,n1)
      bface( 8,n2) = bface(3,n1)
      bface( 9,n2) = bface(4,n1)
      bface(10,n2) = ncell + n2
      !
      nf_min = n1
      nf_max = n2
      !
    end do
  end do
  !
  nf_min = nf_max
  !
  do j = 1,ny
    do i = 1,nx
      !
      n1 = nf_min + 1
      n2 = nf_min + 1 + nx*ny
      !
      bface( 1,n1) = bc_periodic
      bface( 2,n1) = map3dto1d(i,j,1 ,nx,ny)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(5)
      bface( 5,n1) = n2
      bface( 6,n1) = map3dto1d(i,j,nz,nx,ny)
      bface( 8,n1) = Geom_Hexa
      bface( 9,n1) = HostFace3D(6)
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_periodic
      bface( 2,n2) = bface(6,n1)
      bface( 3,n2) = bface(8,n1)
      bface( 4,n2) = bface(9,n1)
      bface( 5,n2) = n1
      bface( 6,n2) = bface(2,n1)
      bface( 8,n2) = bface(3,n1)
      bface( 9,n2) = bface(4,n1)
      bface(10,n2) = ncell + n2
      !
      nf_min = n1
      nf_max = n2
      !
    end do
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  contains
  !
  pure function map3dto1d(i,j,k,nx,ny) result(return_value)
    integer, intent(in) :: i,j,k,nx,ny
    integer :: return_value
    return_value = i + nx*(j-1) + nx*ny*(k-1)
  end function map3dto1d
  !
end subroutine create_channel_flow_grid
!
!###############################################################################
!
subroutine create_postproc_debug_grid
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use ovar,      only : Lref
  use geovar,    only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geovar,    only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geovar,    only : cell_geom,cell_order
  use ovar,      only : VortexGrid_random_perturb
  use ovar,      only : VortexGrid_perturb_percent
  use ovar,      only : grid_format
  !
  !.. Local Scalars ..
  integer :: nc,np,ierr,i,j,k,l,n
  integer :: n1,n2,nf_min,nf_max
  integer :: nc_min,nc_max
  integer :: no_min,no_max
  logical(ldk) :: bnd_face_error
  real(r4) :: x,y,z,dl
  real(r4) :: rbas,radj,rval
  !
  !.. Local Arrays ..
  integer, dimension(1:4)      :: top,bot,ip
  integer, dimension(1:8)      :: cell_nodes
  integer, dimension(1:8,1:24) :: node_order
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: lpoin(:)
  integer,  allocatable :: original_noc(:)
  real(r4), allocatable :: perturb(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_postproc_debug_grid"
  !
  integer, parameter :: cells_per_cpu = 3**3
  integer, parameter :: basic_order(8,3) = reshape( [1,2,3,4,5,6,7,8, &
                                                     2,6,7,3,1,5,8,4, &
                                                     1,5,6,2,4,8,7,3] , &
                                                    shape(basic_order) )
  integer, parameter :: flip(4) = [1,4,3,2]
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Reset the grid format input parameter to identify
  ! the grid as internally generated
  !
  grid_format = Internally_Generated_Grid
  !
  bnd_face_error = fals
  !
  ! Create the node_order array to prepare for randomizing the node orderings
  ! NOTE: We subtract one so node_order can be used as an offset
  !
  n = 0
  do j = 1,3
    bot = basic_order(1:4,j)
    top = basic_order(5:8,j)
    do i = 0,3
      n = n + 1
      node_order(:,n) = [cshift(bot,shift=i),cshift(top,shift=i)] - 1
    end do
    bot = top(flip)
    top = basic_order(flip,j)
    do i = 0,3
      n = n + 1
      node_order(:,n) = [cshift(bot,shift=i),cshift(top,shift=i)] - 1
    end do
  end do
  !
  ! If we are generating the grid on the fly, make the reference length 1
  !
  Lref = one
  !
  nr = 3
  !
  ! Estimate the total number of cells
  !
  nc = nint( real(ncpu*cells_per_cpu,kind=wp)**one3 )
  !
  np = nc + 1 ! number of grid nodes in each direction
  !
  ncell = nc*nc*nc ! total number of grid cells
  nnode = np*np*np ! total number of grid nodes
  nfbnd = 6*nc*nc  ! total number of boundary faces
  !
  ! Compute the cell length for a Cartesian/uniform grid
  !
  dl = real(1,kind=r4) / real(np-1,kind=r4)
  !
  ! Allocate all the arrays needed to define the grid
  !
  allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_nodes",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  allocate ( nodes_of_cell(1:geom_nodes(Geom_Hexa)*ncell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_geom(1:ncell) , source=Geom_Hexa , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_geom",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_order(1:ncell) , source=n_order , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_order",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! ONLY HAVE THE ROOT PROCESSOR CREATE THE GRID SINCE THERE IS A LOT OF
  ! RANDOMIZATION INVOLVED. WE WILL BROADCAST THE RESULTS TO ALL THE OTHER
  ! PROCESSORS AFTERWARDS SO THAT ALL PROCESSORS HAVE THE SAME GRID DATA.
  !
  only_root_creates_grid: if (mypnum == 0) then
    !
    ! #######
    ! #######
    ! ## 1 ##
    ! #######
    ! #######
    !
    ! Create the nodes_of_cell_ptr array
    !
    nodes_of_cell_ptr(1) = 0
    do n = 2,size(nodes_of_cell_ptr)
      nodes_of_cell_ptr(n) = nodes_of_cell_ptr(n-1) + geom_nodes(Geom_Hexa)
    end do
    !
    ! Create the nodes_of_cell array
    !
    do k = 1,nc
      do j = 1,nc
        do i = 1,nc
          !
          n = nodes_of_cell_ptr( map3dto1d(i,j,k,nc) )
          !
          nodes_of_cell(n+1) = map3dto1d(i  ,j  ,k  ,np)
          nodes_of_cell(n+2) = map3dto1d(i+1,j  ,k  ,np)
          nodes_of_cell(n+3) = map3dto1d(i+1,j+1,k  ,np)
          nodes_of_cell(n+4) = map3dto1d(i  ,j+1,k  ,np)
          nodes_of_cell(n+5) = map3dto1d(i  ,j  ,k+1,np)
          nodes_of_cell(n+6) = map3dto1d(i+1,j  ,k+1,np)
          nodes_of_cell(n+7) = map3dto1d(i+1,j+1,k+1,np)
          nodes_of_cell(n+8) = map3dto1d(i  ,j+1,k+1,np)
          !
        end do
      end do
    end do
    !
    ! Randomize the node ordering for all cells
    !
    if (VortexGrid_random_perturb) then
      !
      ! Create a copy of nodes_of_cell before it is randomized
      !
      allocate ( original_noc(1:size(nodes_of_cell)) , source=nodes_of_cell , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"original_noc",1,__LINE__,__FILE__,ierr, &
                       error_message)
      !
      do k = 1,nc
        do j = 1,nc
          do i = 1,nc
            !
            n = map3dto1d(i,j,k,nc)
            !
            l = random(size(node_order,dim=2))
            l = min( max(1,l) , size(node_order,dim=2) )
            n1 = nodes_of_cell_ptr(n)+1
            n2 = nodes_of_cell_ptr(n+1)
            nodes_of_cell(n1:n2) = nodes_of_cell( n1 + node_order(:,l) )
            !
          end do
        end do
      end do
      !
    end if
    !
    ! #######
    ! #######
    ! ## 2 ##
    ! #######
    ! #######
    !
    allocate ( perturb(1:nr,1:nnode) , source=real(0,kind=r4) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"perturb",1,__LINE__,__FILE__,ierr,error_message)
    !
    if (VortexGrid_random_perturb) then
      !
      ! Compute the random perturbations of the interior node coordinates
      !
      rbas = real(abs(VortexGrid_perturb_percent),kind=r4) / real(100,kind=r4)
      radj = rbas / real(2,kind=r4)
      !
      ! Randomly perturb the crap out of everything
      !
      do k = 2,np-1
        do j = 2,np-1
          do i = 2,np-1
            !
            n = map3dto1d(i,j,k,np)
            !
            do l = 1,nr
              rval = random(rbas) - radj
              perturb(l,n) = rval
            end do
            !
          end do
        end do
      end do
      !
    end if
    !
    ! Create the grid node coordinates as a cube from (0,0,0) to (1,1,1)
    !
    do k = 1,np
      !
      z = real(k-1,kind=r4)
      !
      do j = 1,np
        !
        y = real(j-1,kind=r4)
        !
        do i = 1,np
          !
          x = real(i-1,kind=r4)
          !
          n = map3dto1d(i,j,k,np)
          !
          xyz_nodes(1,n) = real( x + perturb(1,n) , kind=wp )
          xyz_nodes(2,n) = real( y + perturb(2,n) , kind=wp )
          xyz_nodes(3,n) = real( z + perturb(3,n) , kind=wp )
          !
        end do
      end do
    end do
    !
    ! Deallocate perturb since it is no longer needed
    !
    deallocate ( perturb , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"perturb",2,__LINE__,__FILE__,ierr,error_message)
    !
    ! #######
    ! #######
    ! ## 3 ##
    ! #######
    ! #######
    !
    allocate ( lpoin(1:nnode) , source=0 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"lpoin",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Create the bface array
    !
    nf_max = 0
    !
    ! Loop through the 6 grid boundaries
    ! NOTE: This corresponds to the grid boundaries in
    !       the order [imin,imax,jmin,jmax,kmin,kmax]
    !
    nf_min = nf_max
    !
    do k = 1,nc
      do j = 1,nc
        !
        nf_min = nf_min + 1
        nf_max = nf_min + nc*nc
        !
        nc_min = map3dto1d(1 ,j,k,nc)
        nc_max = map3dto1d(nc,j,k,nc)
        !
        no_min = 0
        no_max = 0
        !
        if (VortexGrid_random_perturb) then
          !
          ! ******************
          ! *** X-Min Face ***
          ! ******************
          !
          n1 = nodes_of_cell_ptr(nc_min)+1
          n2 = nodes_of_cell_ptr(nc_min+1)
          !
          ! Mark the nodes on this face
          !
          lpoin( original_noc(n1 + Xmin3D) ) = 1
          !
          ! Extract the nodes of this cell in the new randomized order
          !
          cell_nodes(1:8) = nodes_of_cell(n1:n2)
          !
          ! Find which side of the host cell is now on the boundary
          !
          find_xmin_face: do n = 1,6
            ! Node offsets for this face of the cell
            ! NOTE: The +1 is because the indexing in side_nodes starts at 0.
            ip(1:4) = side_nodes(1:4,n,Geom_Hexa) + 1
            ! Check if all nodes for this face have been marked
            if (sum(lpoin(cell_nodes(ip))) == 4) then
              no_min = n                ! Save the cell face on the boundary
              lpoin(cell_nodes(ip)) = 0 ! Reset the marked nodes to zero
              exit find_xmin_face
            end if
          end do find_xmin_face
          !
          ! Make sure that no_min was set
          !
          if (no_min == 0) then
            bnd_face_error = true
            write (iout,1) "X-Min"
          end if
          !
          ! ******************
          ! *** X-Max Face ***
          ! ******************
          !
          n1 = nodes_of_cell_ptr(nc_max)+1
          n2 = nodes_of_cell_ptr(nc_max+1)
          !
          ! Mark the nodes on this face
          !
          lpoin( original_noc(n1 + Xmax3D) ) = 1
          !
          ! Extract the nodes of this cell in the new randomized order
          !
          cell_nodes(1:8) = nodes_of_cell(n1:n2)
          !
          ! Find which side of the host cell is now on the boundary
          !
          find_xmax_face: do n = 1,6
            ! Node offsets for this face of the cell
            ! NOTE: The +1 is because the indexing in side_nodes starts at 0.
            ip(1:4) = side_nodes(1:4,n,Geom_Hexa) + 1
            ! Check if all nodes for this face have been marked
            if (sum(lpoin(cell_nodes(ip))) == 4) then
              no_max = n                ! Save the cell face on the boundary
              lpoin(cell_nodes(ip)) = 0 ! Reset the marked nodes to zero
              exit find_xmax_face
            end if
          end do find_xmax_face
          !
          ! Make sure that no_max was set
          !
          if (no_max == 0) then
            bnd_face_error = true
            write (iout,1) "X-Max"
          end if
          !
        else
          !
          no_min = HostFace3D(1)
          no_max = HostFace3D(2)
          !
        end if
        !
        bface( 1,nf_min) = bc_periodic
        bface( 2,nf_min) = nc_min
        bface( 3,nf_min) = Geom_Hexa
        bface( 4,nf_min) = no_min
        bface( 5,nf_min) = nf_max
        bface( 6,nf_min) = nc_max
        bface( 8,nf_min) = Geom_Hexa
        bface( 9,nf_min) = no_max
        bface(10,nf_min) = ncell + nf_min
        !
        bface( 1,nf_max) = bc_periodic
        bface( 2,nf_max) = bface(6,nf_min)
        bface( 3,nf_max) = bface(8,nf_min)
        bface( 4,nf_max) = bface(9,nf_min)
        bface( 5,nf_max) = nf_min
        bface( 6,nf_max) = bface(2,nf_min)
        bface( 8,nf_max) = bface(3,nf_min)
        bface( 9,nf_max) = bface(4,nf_min)
        bface(10,nf_max) = ncell + nf_max
        !
      end do
    end do
    !
    nf_min = nf_max
    !
    do k = 1,nc
      do i = 1,nc
        !
        nf_min = nf_min + 1
        nf_max = nf_min + nc*nc
        !
        nc_min = map3dto1d(i,1 ,k,nc)
        nc_max = map3dto1d(i,nc,k,nc)
        !
        no_min = 0
        no_max = 0
        !
        if (VortexGrid_random_perturb) then
          !
          ! ******************
          ! *** Y-Min Face ***
          ! ******************
          !
          n1 = nodes_of_cell_ptr(nc_min)+1
          n2 = nodes_of_cell_ptr(nc_min+1)
          !
          ! Mark the nodes on this face
          !
          lpoin( original_noc(n1 + Ymin3D) ) = 1
          !
          ! Extract the nodes of this cell in the new randomized order
          !
          cell_nodes(1:8) = nodes_of_cell(n1:n2)
          !
          ! Find which side of the host cell is now on the boundary
          !
          find_ymin_face: do n = 1,6
            ! Node offsets for this face of the cell
            ! NOTE: The +1 is because the indexing in side_nodes starts at 0.
            ip(1:4) = side_nodes(1:4,n,Geom_Hexa) + 1
            ! Check if all nodes for this face have been marked
            if (sum(lpoin(cell_nodes(ip))) == 4) then
              no_min = n                ! Save the cell face on the boundary
              lpoin(cell_nodes(ip)) = 0 ! Reset the marked nodes to zero
              exit find_ymin_face
            end if
          end do find_ymin_face
          !
          ! Make sure that no_min was set
          !
          if (no_min == 0) then
            bnd_face_error = true
            write (iout,1) "Y-Min"
          end if
          !
          ! ******************
          ! *** Y-Max Face ***
          ! ******************
          !
          n1 = nodes_of_cell_ptr(nc_max)+1
          n2 = nodes_of_cell_ptr(nc_max+1)
          !
          ! Mark the nodes on this face
          !
          lpoin( original_noc(n1 + Ymax3D) ) = 1
          !
          ! Extract the nodes of this cell in the new randomized order
          !
          cell_nodes(1:8) = nodes_of_cell(n1:n2)
          !
          ! Find which side of the host cell is now on the boundary
          !
          find_ymax_face: do n = 1,6
            ! Node offsets for this face of the cell
            ! NOTE: The +1 is because the indexing in side_nodes starts at 0.
            ip(1:4) = side_nodes(1:4,n,Geom_Hexa) + 1
            ! Check if all nodes for this face have been marked
            if (sum(lpoin(cell_nodes(ip))) == 4) then
              no_max = n                ! Save the cell face on the boundary
              lpoin(cell_nodes(ip)) = 0 ! Reset the marked nodes to zero
              exit find_ymax_face
            end if
          end do find_ymax_face
          !
          ! Make sure that no_max was set
          !
          if (no_max == 0) then
            bnd_face_error = true
            write (iout,1) "Y-Max"
          end if
          !
        else
          !
          no_min = HostFace3D(3)
          no_max = HostFace3D(4)
          !
        end if
        !
        bface( 1,nf_min) = bc_periodic
        bface( 2,nf_min) = nc_min
        bface( 3,nf_min) = Geom_Hexa
        bface( 4,nf_min) = no_min
        bface( 5,nf_min) = nf_max
        bface( 6,nf_min) = nc_max
        bface( 8,nf_min) = Geom_Hexa
        bface( 9,nf_min) = no_max
        bface(10,nf_min) = ncell + nf_min
        !
        bface( 1,nf_max) = bc_periodic
        bface( 2,nf_max) = bface(6,nf_min)
        bface( 3,nf_max) = bface(8,nf_min)
        bface( 4,nf_max) = bface(9,nf_min)
        bface( 5,nf_max) = nf_min
        bface( 6,nf_max) = bface(2,nf_min)
        bface( 8,nf_max) = bface(3,nf_min)
        bface( 9,nf_max) = bface(4,nf_min)
        bface(10,nf_max) = ncell + nf_max
        !
      end do
    end do
    !
    nf_min = nf_max
    !
    do j = 1,nc
      do i = 1,nc
        !
        nf_min = nf_min + 1
        nf_max = nf_min + nc*nc
        !
        nc_min = map3dto1d(i,j,1 ,nc)
        nc_max = map3dto1d(i,j,nc,nc)
        !
        no_min = 0
        no_max = 0
        !
        if (VortexGrid_random_perturb) then
          !
          ! ******************
          ! *** Z-Min Face ***
          ! ******************
          !
          n1 = nodes_of_cell_ptr(nc_min)+1
          n2 = nodes_of_cell_ptr(nc_min+1)
          !
          ! Mark the nodes on this face
          !
          lpoin( original_noc(n1 + Zmin3D) ) = 1
          !
          ! Extract the nodes of this cell in the new randomized order
          !
          cell_nodes(1:8) = nodes_of_cell(n1:n2)
          !
          ! Find which side of the host cell is now on the boundary
          !
          find_zmin_face: do n = 1,6
            ! Node offsets for this face of the cell
            ! NOTE: The +1 is because the indexing in side_nodes starts at 0.
            ip(1:4) = side_nodes(1:4,n,Geom_Hexa) + 1
            ! Check if all nodes for this face have been marked
            if (sum(lpoin(cell_nodes(ip))) == 4) then
              no_min = n                ! Save the cell face on the boundary
              lpoin(cell_nodes(ip)) = 0 ! Reset the marked nodes to zero
              exit find_zmin_face
            end if
          end do find_zmin_face
          !
          ! Make sure that no_min was set
          !
          if (no_min == 0) then
            bnd_face_error = true
            write (iout,1) "Z-Min"
          end if
          !
          ! ******************
          ! *** Z-Max Face ***
          ! ******************
          !
          n1 = nodes_of_cell_ptr(nc_max)+1
          n2 = nodes_of_cell_ptr(nc_max+1)
          !
          ! Mark the nodes on this face
          !
          lpoin( original_noc(n1 + Zmax3D) ) = 1
          !
          ! Extract the nodes of this cell in the new randomized order
          !
          cell_nodes(1:8) = nodes_of_cell(n1:n2)
          !
          ! Find which side of the host cell is now on the boundary
          !
          find_zmax_face: do n = 1,6
            ! Node offsets for this face of the cell
            ! NOTE: The +1 is because the indexing in side_nodes starts at 0.
            ip(1:4) = side_nodes(1:4,n,Geom_Hexa) + 1
            ! Check if all nodes for this face have been marked
            if (sum(lpoin(cell_nodes(ip))) == 4) then
              no_max = n                ! Save the cell face on the boundary
              lpoin(cell_nodes(ip)) = 0 ! Reset the marked nodes to zero
              exit find_zmax_face
            end if
          end do find_zmax_face
          !
          ! Make sure that no_max was set
          !
          if (no_max == 0) then
            bnd_face_error = true
            write (iout,1) "Z-Max"
          end if
          !
        else
          !
          no_min = HostFace3D(5)
          no_max = HostFace3D(6)
          !
        end if
        !
        bface( 1,nf_min) = bc_periodic
        bface( 2,nf_min) = nc_min
        bface( 3,nf_min) = Geom_Hexa
        bface( 4,nf_min) = no_min
        bface( 5,nf_min) = nf_max
        bface( 6,nf_min) = nc_max
        bface( 8,nf_min) = Geom_Hexa
        bface( 9,nf_min) = no_max
        bface(10,nf_min) = ncell + nf_min
        !
        bface( 1,nf_max) = bc_periodic
        bface( 2,nf_max) = bface(6,nf_min)
        bface( 3,nf_max) = bface(8,nf_min)
        bface( 4,nf_max) = bface(9,nf_min)
        bface( 5,nf_max) = nf_min
        bface( 6,nf_max) = bface(2,nf_min)
        bface( 8,nf_max) = bface(3,nf_min)
        bface( 9,nf_max) = bface(4,nf_min)
        bface(10,nf_max) = ncell + nf_max
        !
      end do
    end do
    !
    if (allocated(lpoin)) then
      deallocate ( lpoin , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"lpoin",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    if (allocated(original_noc))  then
      deallocate ( original_noc , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"original_noc",2,__LINE__,__FILE__,ierr, &
                       error_message)
    end if
    !
  end if only_root_creates_grid
  !
  ! Have the root processor broadcast the value of bnd_face_error so that
  ! the simulation can be aborted before continuing if there was an error
  !
  if (ncpu > 1) then
    call mpi_bcast(bnd_face_error,1_int_mpi, &
                   MPI_LOGICAL,0_int_mpi,MPI_COMM_WORLD,mpierr)
  end if
  !
  if (bnd_face_error) then
    write (error_message,2)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Have the root processor broadcast all the grid
  ! information to all the other processors.
  !
  if (ncpu > 1) then
    !
    ! Call mpi_barrier just to make sure the processor are syncronized
    !
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    !
    call mpi_bcast(xyz_nodes,size(xyz_nodes,kind=int_mpi), &
                   mpi_flttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
    !
    call mpi_bcast(nodes_of_cell_ptr,size(nodes_of_cell_ptr,kind=int_mpi), &
                   mpi_inttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
    !
    call mpi_bcast(nodes_of_cell,size(nodes_of_cell,kind=int_mpi), &
                   mpi_inttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
    !
    call mpi_bcast(bface,size(bface,kind=int_mpi), &
                   mpi_inttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
    !
    ! Call mpi_barrier just to make sure the processor are syncronized
    !
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    !
  end if
  !
  ! Scale the grid node coordinates to the cube (0,0,0) to (1,1,1)
  !
 !xyz_nodes = dl * xyz_nodes
  !
  ! Scale the grid node coordinates to the cube (-1,-1,-1) to (1,1,1)
  !
 !xyz_nodes = two*xyz_nodes - one
  !
  ! Scale the grid node coordinates to the cube (-pi,-pi,-pi) to (pi,pi,pi)
  !
 !xyz_nodes = PI * xyz_nodes
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (" Error finding new side of host cell on ",a," boundary face!")
  2 format ("An error was detected finding the side of one or more host ", &
            "cells that is on the boundary after randomizing the order of ", &
            "the nodes defining each cell!")
  !
  contains
  !
  pure function map3dto1d(i,j,k,n) result(return_value)
    integer, intent(in) :: i,j,k,n
    integer :: return_value
    return_value = i + n*(j-1) + n*n*(k-1)
  end function map3dto1d
  !
end subroutine create_postproc_debug_grid
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
  grid%xyz = xyz_nodes(:,:) ! USING F2003 AUTO-REALLOCATION
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
 pure &
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
      case (Geom_Hexa)
        return_value%pts( hexa(:) ) = elem%pts(:)
    end select
  else
    select case (elem%prop%geom)
      case (Geom_Quad)
        return_value%pts(:) = elem%pts( n1-1 + quad(:) )
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
  !.. Local Arrays ..
  real(wp), dimension(1:nr,1:8) :: cell_nodes
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
    cell_nodes(1:nr,1:np) = xyz_nodes(1:nr,nodes_of_cell(n1:n2))
    !
    select case (cell_geom(n))
      case (Geom_Edge)
        cell_vol = length_edge( cell_nodes(1:nr,1:np) )
      case (Geom_Quad)
        cell_vol = area_quad( cell_nodes(1:nr,1:np) )
      case (Geom_Tria)
        cell_vol = area_tria( cell_nodes(1:nr,1:np) )
      case (Geom_Tetr)
        cell_vol = vol_tet( cell_nodes(1:nr,1:np) )
      case (Geom_Pyra)
        cell_vol = vol_pyr( cell_nodes(1:nr,1:np) )
      case (Geom_Pris)
        cell_vol = vol_pri( cell_nodes(1:nr,1:np) )
      case (Geom_Hexa)
        cell_vol = vol_hex( cell_nodes(1:nr,1:np) )
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
  !
  ! Now lets do it over all the cells in the grid. This includes checking
  ! the cells on the grid boundaries a second time but they shouldnt be
  ! adjusted this time and the additional cost shouldnt hurt too much.
  !
  do n = 1,ncell
    !
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    !
    ! Compute the signed area/volume of the current host cell
    !
    cell_nodes(1:nr,1:np) = xyz_nodes(1:nr,nodes_of_cell(n1:n2))
    !
    select case (cell_geom(n))
      case (Geom_Edge)
        cell_vol = length_edge( cell_nodes(1:nr,1:np) )
      case (Geom_Quad)
        cell_vol = area_quad( cell_nodes(1:nr,1:np) )
      case (Geom_Tria)
        cell_vol = area_tria( cell_nodes(1:nr,1:np) )
      case (Geom_Tetr)
        cell_vol = vol_tet( cell_nodes(1:nr,1:np) )
      case (Geom_Pyra)
        cell_vol = vol_pyr( cell_nodes(1:nr,1:np) )
      case (Geom_Pris)
        cell_vol = vol_pri( cell_nodes(1:nr,1:np) )
      case (Geom_Hexa)
        cell_vol = vol_hex( cell_nodes(1:nr,1:np) )
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
      else
        !
        write (iout,2)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
        !
      end if
    end if
    !
  end do
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
subroutine adjust_for_same_dimension_ghost_cells()
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
  integer  :: g1,g2,i,h1,h2,nf,nn,ni,np
  integer  :: ghost_cell,host_cell,host_geom,inode,host_side
  integer  :: pf,p1,p2
  integer  :: ni1,ni2,ierr
  integer  :: periodic_cell,periodic_geom,periodic_side,pcnods
  real(wp) :: theta
  integer :: nfnods
  !
  !.. Local Arrays ..
  integer,  dimension(1:4)      :: ip,periodic_ip
  integer,  dimension(1:4)      :: face_nodes
  integer,  dimension(1:4)      :: pface_nodes
  real(wp), dimension(1:nr)     :: ptf
  real(wp), dimension(1:nr)     :: xyz_h1
  real(wp), dimension(1:nr)     :: xyz_h2
  real(wp), dimension(1:nr)     :: xyz_ni
  real(wp), dimension(1:nr)     :: ds,ds1,ds2
  real(wp), dimension(1:nr)     :: rnrm,pnrm
  real(wp), dimension(1:nr,1:8) :: periodic_nodes
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "adjust_for_same_dimension_ghost_cells"
  !
  integer, parameter :: rotation = 0
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
  ! Set the number of nodes for each ghost cell
  !
  do nf = 1,nfbnd
    !
    if (bface(1,nf) /= bc_periodic) then
      !
      ! Not a periodic boundary, number of nodes is the same as interior
      !
      inode = geom_nodes( bface(3,nf) )
      !
    else
      !
      ! Periodic boundary, get the number of
      ! nodes for the connected periodic cell
      !
      pf = bface(5,nf)
      !
      inode = geom_nodes( bface(3,pf) )
      !
    end if
    !
    nodes_of_cell_ptr(ncell+nf+1) = inode
    !
  end do
  !
  ! Reshuffle the ghost cell portion of nodes_of_cell_ptr
  !
  do i = ncell+2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(i) = nodes_of_cell_ptr(i) + nodes_of_cell_ptr(i-1)
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
  ! Loop over the boundary faces and find the nodes of each ghost cell.
  ! The first two nodes are on the boundary face while the other nodes
  ! are newly created and to be defined in step 3. The nodes are in
  ! counter-clockwise ordering.
  !
  inode = nnode
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
    ! Node offsets for this face of the host cell
    !
    ip(1:nfnods) = side_nodes(1:nfnods,host_side,host_geom)
    !
    ! Use the node offsets to get the actual indices for these nodes
    !
    face_nodes(1:nfnods) = nodes_of_cell( h1 + ip(1:nfnods) )
    !
    ! Copy the face node indices into their correct locations for the ghost cell
    !
    nodes_of_cell(g1:g2) = ghost_nodes(face_nodes(1:nfnods),host_geom,host_side)
   !!
   !! Assign nodes defining boundary face from the host cell to the ghost cell
   !! NOTE: We have to reverse the rotation/ordering of the face nodes so that
   !!       the nodes of the ghost cell are defined such that the normal of
   !!       each of its faces points outwards (this is the same as the interior
   !!       cells). For a 2D host cell, there are only two nodes so we just
   !!       need to switch them, i.e. node1 <=> node2.  For a 3D host cell,
   !!       switching either nodes 1 & 3 or nodes 3 & 4 will do the trick.
   !!       Since triangle faces only have three nodes, we will switch nodes
   !!       1 & 3 for both triangle and quad faces so that we dont have to
   !!       distinguish between the two.
   !!
   !! Copy over the face nodes to the ghost cell,
   !! keeping the same ordering as the host cell.
   !!
   !nodes_of_cell(g1:g1-1+nfnods) = face_nodes(1:nfnods)
   !!
   !! Get the index of the node to switch with node 1.
   !! NOTE: switch_index = nfnods (i.e. 2) for a 2D host cell
   !!       switch_index = nfnods (i.e. 3) for a triangular boundary face
   !!       switch_index = 3               for a quadrilateral boundary face
   !! ADDITIONAL NOTE: switch_index also evaluates to nfnods (i.e. 1) for
   !!                  a 1D host cell which is the correct value if the
   !!                  rest of the code is ever set up to run 1D.
   !!
   !switch_index = min( nfnods , 3 )
   !!
   !! Switch the two face nodes required to reverse the face node ordering
   !!
   !nodes_of_cell(g1) = face_nodes(switch_index)
   !nodes_of_cell(g1-1+switch_index) = face_nodes(1)
    !
    if (all(host_geom /= [Geom_1D,Geom_2D,Geom_3D])) then
      write (error_message,100)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
    !
    ! Create new nodes for the remaining ghost cell nodes
    !
    do i = g1,g2
      if (nodes_of_cell(i) < 0) then
        inode = inode + 1
        nodes_of_cell(i) = inode
      end if
    end do
    !
  end do
  !
  ! Total number of nodes after the addition of
  ! the newly created nodes for the ghost cells
  !
  n_totnod = inode
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Reallocate cell_geom to account for ghost cells
  !
  call reallocate( cell_geom , n_totcel )
  !
  ! Reallocate cell_order to account for ghost cells
  !
  call reallocate( cell_order , n_totcel )
  !
  ! Fill in the information for the ghost cells
  !
  do nf = 1,nfbnd
    cell_geom(nf+ncell) = cell_geom(bface(2,nf))
    cell_order(nf+ncell) = cell_order(bface(2,nf))
  end do
  !
  ! #######
  ! ## 4 ##
  ! #######
  !
  ! Reallocate xyz_nodes to account for ghost nodes
  !
  call reallocate( xyz_nodes , new_size2=n_totnod )
  !
  ! Loop over boundary faces to compute the
  ! coordinates of the newly created ghost nodes
  !
  do nf = 1,nfbnd
    !
    ! Ghost cell and the corresponding node indices
    !
    ghost_cell = nf + ncell
    !
    g1 = nodes_of_cell_ptr(ghost_cell)+1
    g2 = nodes_of_cell_ptr(ghost_cell+1)
    !
    ! Interior cell and the corresponding node indices
    !
    host_cell = bface(2,nf)
    !
    h1 = nodes_of_cell_ptr(host_cell)+1
    h2 = nodes_of_cell_ptr(host_cell+1)
    !
    ! Information about boundary face
    !
    host_geom = bface(3,nf)
    host_side = bface(4,nf)
    !
    ! Get the number of nodes defining this face
    !
    nfnods = num_face_nodes(host_side,host_geom)
    !
    ! Node offsets for this face of the host cell
    !
    ip(1:nfnods) = side_nodes(1:nfnods,host_side,host_geom)
    !
    ! Use the node offsets to get the actual indices for these nodes
    !
    face_nodes(1:nfnods) = nodes_of_cell( h1 + ip(1:nfnods) )
    !
    ! Outward unit normal to boundary face.
    !
    rnrm = face_normal( xyz_nodes(1:nr,face_nodes(1:nfnods)) )
    rnrm = unit_vector( rnrm )
    !
    ! Store the coordinates for the first and last node of the cell
    !
    xyz_h1(:) = xyz_nodes(:,h1)
    xyz_h2(:) = xyz_nodes(:,h2)
    !
    ! Conditional to treat periodic boundary conditions separately
    !
    if (bface(1,nf) /= bc_periodic) then
      !
      ! This is not a periodic boundary face
      !
      ! Center of boundary face by averaging the
      ! coordinates of the two boundary nodes
      !
      ptf(1:nr) = node_center( xyz_nodes(1:nr,face_nodes(1:nfnods)) )
      !
      ! Coordinates for nodes g1 to g1-1+nfnods of the ghost cell have already
      ! been defined so we need to define nodes g1+nfnods to g2.
      !
      if (host_geom == 3) then
        !
        ! Ghost node 3 for triangle
        !
        ni = nodes_of_cell( h1 + mod(host_side+1,host_geom) )
        nn = nodes_of_cell( g1 + 2 )
        !
        ! Store the coordinates for node ni
        !
        xyz_ni(:) = xyz_nodes(:,ni)
        !
        ! Get the coordinates of ghost node 3
        !
        xyz_nodes(:,nn) = merge(rotate(xyz_h1(:),xyz_h2(:),xyz_ni(:)), &
                                reflect(xyz_ni(:),ptf,rnrm), &
                                rotation == 1)
        !
      else if (host_geom == 4) then
        !
        ! Ghost node 3 for quadrilateral
        !
        nn = nodes_of_cell( g1 + 2 )
        ni = nodes_of_cell( h1 + &
                            mod(host_side+2-max(0,min(1,rotation)),host_geom) )
        !
        ! Store the coordinates for node ni
        !
        xyz_ni(:) = xyz_nodes(:,ni)
        !
        ! Get the coordinates of ghost node 3
        !
        xyz_nodes(:,nn) = merge(rotate(xyz_h1(:),xyz_h2(:),xyz_ni(:)), &
                                reflect(xyz_ni(:),ptf,rnrm), &
                                rotation == 1)
        !
        ! Ghost node 4 for quadrilateral
        !
        nn = nodes_of_cell( g1 + 3 )
        ni = nodes_of_cell( h1 + &
                            mod(host_side+1+max(0,min(1,rotation)),host_geom) )
        !
        ! Store the coordinates for node ni
        !
        xyz_ni(:) = xyz_nodes(:,ni)
        !
        ! Get the coordinates of ghost node 4
        !
        xyz_nodes(:,nn) = merge(rotate(xyz_h1(:),xyz_h2(:),xyz_ni(:)), &
                                reflect(xyz_ni(:),ptf,rnrm), &
                                rotation == 1)
        !
      end if
      !
    else
      !
      ! This is a periodic boundary face
      !
      ! Boundary face linked to this periodic boundary face
      !
      pf = bface(5,nf)
      !
      ! Linked cell and the corresponding node indices
      !
      periodic_cell = bface(2,pf)
      !
      p1 = nodes_of_cell_ptr(periodic_cell)+1
      p2 = nodes_of_cell_ptr(periodic_cell+1)
      !
      ! Information about the linked face
      !
      periodic_geom = bface(3,pf)
      periodic_side = bface(4,pf)
      !
      ! Check that the face geometry matches across this periodic boundary
      !
      if (num_face_nodes(periodic_side,periodic_geom) /= nfnods) then
        write (error_message,101)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      !
      ! Get the number of nodes defining the periodic cell
      !
      pcnods = geom_nodes(periodic_geom)
      !
      ! Node coordinates for periodic cell
      !
      periodic_nodes(1:nr,1:pcnods) = xyz_nodes(1:nr,nodes_of_cell(p1:p2))
      !
      ! Node offsets for this face of the periodic cell
      !
      periodic_ip(1:nfnods) = side_nodes(1:nfnods,periodic_side,periodic_geom)
      !
      ! Use the node offsets to get the actual indices for these nodes
      !
      pface_nodes(1:nfnods) = nodes_of_cell( p1 + periodic_ip(1:nfnods) )
      !
      ! Outward unit normal to boundary face.
      !
      pnrm = face_normal( xyz_nodes(1:nr,pface_nodes(1:nfnods)) )
      pnrm = unit_vector( pnrm )
      !
      ! Center of boundary face for periodic_cell
      !
      ptf(1:nr) = node_center( xyz_nodes(1:nr,pface_nodes(1:nfnods)) )
      !
      ! Nodes on the boundary face for periodic_cell
      ! NOTE: NOW STORED IN pface_nodes(1:nfnods)
      !
      ! Nodes on the boundary face for host_cell
      ! NOTE: NOW STORED IN face_nodes(1:nfnods)
      !
      ! ########################################################################
      ! ########################################################################
      ! ########################################################################
      ! ##############  ||||||||||                   ||||||||||  ###############
      ! ##############  VVVVVVVVVV  NEED TO WORK ON  VVVVVVVVVV  ###############
      ! ########################################################################
      ! ########################################################################
      ! ########################################################################
      !
      ! Get the angle between the two outward unit normals
      !
      theta = acos(dot(pnrm,rnrm))
      ! NEED ANOTHER ANGLE FOR 3D !!!!!!
      !
      ! Translation of the periodic_cell nodes so that ptf is at the origin
      !
      do i = 1,pcnods
        periodic_nodes(1:nr,i) = periodic_nodes(1:nr,i) - ptf(1:nr)
      end do
      !
      ! Rotate the nodes of periodic_cell about ptf theta+180 degrees
      !
      do i = 1,pcnods
        periodic_nodes(1:nr,i) = rotate_2d( periodic_nodes(1:nr,i) , theta+pi )
      end do
      !
      ! Translation of the periodic_cell nodes back to original ptf
      !
      do i = 1,pcnods
        periodic_nodes(1:nr,i) = periodic_nodes(1:nr,i) + ptf(1:nr)
      end do
      !
      ! Find the distances between the two corresponding nodes of
      ! the two boundary faces. The two boundary faces should now be
      ! parallel to each other but with opposite outward unit normals
      ! and therefore the distances should be equal.
      !
     !ds1(:) = xyz_nodes(:,ni1) - periodic_nodes(:,1+mod(periodic_side,pcnods))
     !ds2(:) = xyz_nodes(:,ni2) - periodic_nodes(:,periodic_side)
     !!
     !if (abs(ds1(1)-ds2(1))>eps6 .or. abs(ds1(2)-ds2(2))>eps6) then
     !  write (*,1) host_cell,periodic_cell,abs(ds1-ds2)
     !end if
     !!
     !! Move the interior nodes of periodic_cell to the ghost nodes of host_cell
     !!
     !do i = nfnods+1,pcnods
     !  !
     !  nn = nodes_of_cell( g1 + i-1 )
     !  np = nodes_of_cell( p1 + mod(periodic_side+i-2,pcnods) )
     !  !
     !  xyz_nodes(:,nn) = xyz_nodes(:,np) + ds1(:)
     !  !
     !end do
      !
      ! ########################################################################
      ! ########################################################################
      ! ########################################################################
      ! ##############  AAAAAAAAAA  NEED TO WORK ON  AAAAAAAAAA  ###############
      ! ##############  ||||||||||                   ||||||||||  ###############
      ! ########################################################################
      ! ########################################################################
      ! ########################################################################
      !
    end if
    !
  end do
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
  !
end subroutine adjust_for_same_dimension_ghost_cells
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
  use flowvar, only : face_var,faceusp_array
  use parallel_mod, only : exchange_faceusp,wait_faceusp
  use parallel_mod, only : exchange_faceusp_dt,wait_faceusp_dt
  !
  !.. Local Scalars ..
  integer :: i,j,nf,nfp,kf,kl,kr,l,ierr
  integer :: l1,l2,r1,r2,npl,npr,pf,f0,ntfp
  integer :: n1,n2,n0,np,bc_type
  integer :: face_geom,face_order
  integer :: left_cell,left_geom,left_order
  integer :: right_cell,right_geom,right_order
  integer :: l1file,r1file,d1file
  integer :: l2file,r2file,d2file
  real(wp) :: a,val
  character(len=25) :: face_string
  character(len=100) :: array_name
  !
  !.. Local Arrays ..
  type(face_var), allocatable :: face_xyz(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_face_flux_point_orderings"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
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
  allocate ( faceusp_array(1:ntfp,1:nq,1:2) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faceusp_array",1,__LINE__,__FILE__,ierr,error_message)
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
      do kf = 1,nfp
        !
        kr = face(pf)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                          rot(face(pf)%left%rotation)% &
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
         !face_xyz(nf)%v(l,kf,2) = val
          face_xyz(nf)%v(l,nfp-kf+1,2) = val
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
  faceusp_array(:,:,:) = zero
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
        faceusp_array(f0+kf,l,1) = face_xyz(nf)%v(l,kf,1)
      end do
    end do
   !end do
   !!
    if (nf > nfbnd) then
      do kf = 1,nfp
        do l = 1,nr
          faceusp_array(f0+kf,l,2) = face_xyz(nf)%v(l,kf,2)
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
      call transfer_periodic(n1,n2,faceusp_array,bc_grps(i)%periodic_list)
      !
    else
      !
      do kf = 1,np
        do l = 1,nr
          faceusp_array(n0+kf,l,2) = faceusp_array(n0+kf,l,1)
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
    face_string = ""
    if (nf < nfbnd) then
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
      write (l2file,11) (merge(zero,faceusp_array(f0+kf,l,1), &
                                abs(faceusp_array(f0+kf,l,1))<=eps14),l=1,nr)
      !
      write (r1file,11) (merge(zero,face_xyz(nf)%v(l,kf,2), &
                                abs(face_xyz(nf)%v(l,kf,2))<=eps14),l=1,nr)
      write (r2file,11) (merge(zero,faceusp_array(f0+kf,l,2), &
                                abs(faceusp_array(f0+kf,l,2))<=eps14),l=1,nr)
      !
      write (d1file,11) (abs( face_xyz(nf)%v(l,kf,2) - &
                              face_xyz(nf)%v(l,kf,1)) , l=1,nr), &
                        norm2( face_xyz(nf)%v(1:nr,kf,2) - &
                               face_xyz(nf)%v(1:nr,kf,1) )
      write (d2file,11) (abs( faceusp_array(f0+kf,l,2) - &
                              faceusp_array(f0+kf,l,1)) , l=1,nr), &
                        norm2( faceusp_array(f0+kf,1:nr,2) - &
                               faceusp_array(f0+kf,1:nr,1) )
      !
    end do
    !
  end do
  !
  deallocate ( face_xyz , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_xyz",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( faceusp_array , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faceusp_array",2,__LINE__,__FILE__,ierr,error_message)
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
 12 format ("FACE ",i0,4x,a)
  !
end subroutine check_face_flux_point_orderings
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
  use flowvar, only : face_var,faceusp_array
  use parallel_mod, only : exchange_faceusp,wait_faceusp
  use parallel_mod, only : exchange_faceusp_dt,wait_faceusp_dt
  !
  !.. Local Scalars ..
  integer :: i,j,nf,nc,np,n1,n2,l,n,sf,kf,lf,pf,nfp,ierr
  integer :: this_geom,this_order
  integer :: face_geom,face_order,f0,ntfp
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
  allocate ( faceusp_array(1:ntfp,1:nq,1:2) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faceusp_array",1,__LINE__,__FILE__,ierr,error_message)
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
      do kf = 1,nfp
        do l = 1,nr
         !face_xyz(nf)%v(l,kf,2) = face_xyz(pf)%v(l,kf,1)
          face_xyz(nf)%v(l,kf,2) = face_xyz(pf)%v(l,nfp-kf+1,1)
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
  faceusp_array(:,:,:) = zero
  !
  do nf = 1,nface
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    do i = 1,2
      do kf = 1,nfp
        do l = 1,nr
          faceusp_array(f0+kf,l,i) = face_xyz(nf)%v(l,kf,i)
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
      write (l2file,11) (merge(zero,faceusp_array(f0+kf,l,1), &
                                abs(faceusp_array(f0+kf,l,1))<=eps14),l=1,nr)
      !
      write (r1file,11) (merge(zero,face_xyz(nf)%v(l,kf,2), &
                                abs(face_xyz(nf)%v(l,kf,2))<=eps14),l=1,nr)
      write (r2file,11) (merge(zero,faceusp_array(f0+kf,l,2), &
                                abs(faceusp_array(f0+kf,l,2))<=eps14),l=1,nr)
      !
      write (d1file,11) (abs( face_xyz(nf)%v(l,kf,2) - &
                              face_xyz(nf)%v(l,kf,1) ) , l=1,nr), &
                       norm2( face_xyz(nf)%v(1:nr,kf,2) - &
                              face_xyz(nf)%v(1:nr,kf,1) )
      write (d2file,11) (abs( faceusp_array(f0+kf,l,2) - &
                              faceusp_array(f0+kf,l,1) ) , l=1,nr), &
                       norm2( faceusp_array(f0+kf,1:nr,2) - &
                              faceusp_array(f0+kf,1:nr,1) )
      !
    end do
    !
  end do
  !
  deallocate ( face_xyz , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_xyz",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( faceusp_array , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faceusp_array",2,__LINE__,__FILE__,ierr,error_message)
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
  use flowvar, only : face_var,faceusp_array
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
  allocate ( faceusp_array(1:ntfp,1:nq,1:2) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faceusp_array",1,__LINE__,__FILE__,ierr,error_message)
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
  faceusp_array(:,:,:) = zero
  !
  do nf = 1,nfbnd
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    do i = 1,2
      do kf = 1,nfp
        do l = 1,nq
          faceusp_array(f0+kf,l,i) = face_xyz(nf)%v(l,kf,i)
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
      write (l2file,13) (nint(faceusp_array(f0+kf,l,1)),l=1,nq)
      !
      write (r1file,13) (nint(face_xyz(nf)%v(l,kf,2)),l=1,nq)
      write (r2file,13) (nint(faceusp_array(f0+kf,l,2)),l=1,nq)
      !
      write (d1file,11) (abs( face_xyz(nf)%v(l,kf,2) - &
                              face_xyz(nf)%v(l,kf,1) ) , l=1,nq), &
                       norm2( face_xyz(nf)%v(1:nq,kf,2) - &
                              face_xyz(nf)%v(1:nq,kf,1) )
      write (d2file,11) (abs( faceusp_array(f0+kf,l,2) - &
                              faceusp_array(f0+kf,l,1) ) , l=1,nq), &
                       norm2( faceusp_array(f0+kf,1:nq,2) - &
                              faceusp_array(f0+kf,1:nq,1) )
      !
    end do
    !
  end do
  !
  deallocate ( face_xyz , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_xyz",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( faceusp_array , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faceusp_array",2,__LINE__,__FILE__,ierr,error_message)
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
  use ovar,      only : work_units,uverr,total_volume
  use ovar,      only : time,time_ref,time_scaling_factor
  use ovar,      only : output_time_averaging,ave_start_time
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
    do i = 1,nr
      write (iores,2) uvw(i) // trim(suffix)
    end do
    write (iores,2) "p" // trim(suffix)
    write (iores,2) "||V||" // trim(suffix)
    if (write_work_units) then
      write (iores,2) "Work units",'"'
    else
      write (iores,2) "Time",'"'
    end if
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
  write (iores,5) icount,(log10(u2(i)),i=1,nq+1), &
                         (log10(v2(i)),i=1,nq+1), &
                         comp_cost
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
subroutine write_parallel_cgns(iflag)
  !
  !.. Use Statements ..
  use ovar, only : itcur,itrst,num_timesteps,time,time_ref
  use ovar, only : write_TaylorGreen_full_solution
  use ovar, only : completely_disable_cgns
  use ovar, only : output_time_averaging
  use ovar, only : time_scaling_factor
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iflag
  !
  !.. Local Scalars ..
  integer :: n,ierr,final_iter_digits
  !
  character(len=CGLEN) :: tname
  character(len=100)   :: char_num
  character(len=170)   :: sol_file
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_parallel_cgns"
  !
continue
  !
  if (completely_disable_cgns) return
  !
  if (disable_cgns_output) then
    if (iflag == add_to_CGNS) then
      if (.not. write_TaylorGreen_full_solution) then
        return
      end if
    end if
  end if
  !
  call debug_timer(entering_procedure,pname)
  !
  if (iflag == open_CGNS) then
    !
    ! Create the format string for creating the CGNS file names
    !
    write (char_num,4) itrst + num_timesteps
    final_iter_digits = max(4,len_trim(adjustl(char_num)))
    !
    write (cgnsfile_fmt,3) final_iter_digits
    !
    ! Create the CGNS grid file
    !
    call write_parallel_cgns_grid_file
    !
    ! Initialize the zone counter
    !
    cgns_zone_counter = 0
    !
  end if
  !
  ! Get the name for the current solution node
  !
  write (tname,'("FlowSolution",i0)') cgns_zone_counter
 !write (sol_file,'("solution.",i0,".cgns")') cgns_zone_counter
  if (iflag == -3) then
    write (sol_file,cgnsfile_fmt) max(0,itcur-1),".NaN"
  else
    write (sol_file,cgnsfile_fmt) itcur,""
  end if
  !
  ! Reallocate the soltime and solname arrays so we
  ! can concatenate this solution onto those arrays
  !
  call reallocate(soltime,cgns_zone_counter+1)
  call reallocate(solname,cgns_zone_counter+1)
  call reallocate(solfile,cgns_zone_counter+1)
  soltime(size(soltime)) = real(time*time_ref*time_scaling_factor,kind(soltime))
  solname(size(solname)) = tname
  solfile(size(solfile)) = sol_file
  !
  ! Write the solution file
  !
  call write_parallel_cgns_solution_file(sol_file,tname)
  !
  ! Write the solution file
  !
  if (output_time_averaging) then
    call write_parallel_cgns_time_ave_file
  end if
  !
  ! ################################################
  ! #####   FINISHED WRITING FIELD VARIABLES   #####
  ! ################################################
  !
  call set_cgns_queue(STOP_CGNS_QUEUE,pname,cgns_grid_file)
  !
  ! Write the time iterative solution file
  !
  if (any(iflag == [close_CGNS,-3])) then
    call write_parallel_cgns_time_file
  end if
  !
  ! Update the zone counter
  !
  cgns_zone_counter = cgns_zone_counter + 1
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("#############################################",/, &
            "    FOR PROCESSOR #",i0)
  2 format (a," = ",i0)
  3 format ("('solution.',i0.",i0,",a,'.cgns')")
  4 format (i0)
  !
end subroutine write_parallel_cgns
!
!###############################################################################
!
subroutine write_parallel_cgns_time_file()
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : profile_io_cgns,profile_io_cgns_subs
  !
  !.. Local Scalars ..
  integer :: n
  !
  integer(CBT) :: ifile_t,ibase_t,izone_t,iiter_t
  !
  real(wp) :: beg_wall_time
  real(wp) :: end_wall_time
  real(wp) :: cgns_wall_time
  !
  character(len=CGLEN) :: sect_name
  character(len=300)   :: link_path
  !
  !.. Local Arrays ..
  integer(CST) :: idata(1:2)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_parallel_cgns_time_file"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
#ifndef SPECIAL_FOR_PGI
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      write (iout,1) pname
      flush (iout)
      beg_wall_time = mpi_wtime()
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  iiter_t = 1_CBT
  !
  call cgp_mpi_comm_f(mpi_comm_self_val,cgierr)
  !
  if (mypnum == glb_root) then
    !
    cgns_time_file = trim(adjustl(cgns_dir)) // cgns_time_fname
    !
    call cgns_timer(start_cgns_file_timer,profile_io_cgns_subs, &
                    cgns_time_file,override_barrier=true)
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_open_f",override_barrier=true)
    call cgp_open_f(cgns_time_file,CG_MODE_WRITE,ifile_t,cgierr)
    call cgns_error(pname,cgns_time_file,"cgp_open_f", &
                    cgierr,__LINE__,__FILE__, &
                    "Mode = CG_MODE_WRITE",not_parallel=true)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_open_f",override_barrier=true)
    !
    ! Set the CGNS base (who knows what this is for)
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_base_write_f",override_barrier=true)
    call cg_base_write_f(ifile_t,base_name,int(nr,kind=CBT), &
                         int(nr,kind=CBT),ibase_t,cgierr)
    call cgns_error(pname,cgns_time_file,"cg_base_write_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_base_write_f",override_barrier=true)
    !
    ! Write the title for this solution
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_zone_write_f",override_barrier=true)
    call cg_zone_write_f(ifile_t,ibase_t,zone_name,zone_size, &
                         UNSTRUCTURED,izone_t,cgierr)
    call cgns_error(pname,cgns_time_file,"cg_zone_write_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_zone_write_f",override_barrier=true)
    !
    ! Create a link to the grid coordinates in the grid file
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_link_write_f",override_barrier=true)
    call cg_goto_f(ifile_t,ibase_t,cgierr, &
                   "Zone_t",izone_t, &
                   CGNS_GOTO_END)
    call cgns_error(pname,cgns_time_file,"cg_goto_f",cgierr, &
                    __LINE__,__FILE__,not_parallel=true)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_link_write_f",override_barrier=true)
    !
    link_path = base_name//"/"//zone_name//"/GridCoordinates"
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_link_write_f",override_barrier=true)
    call cg_link_write_f("GridCoordinates",trim(cgns_base_fname), &
                         trim(link_path),cgierr)
    call cgns_error(pname,cgns_time_file,"cg_link_write_f(GridCoordinates)", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_link_write_f",override_barrier=true)
    !
    ! Create a link to the element connectivity in the grid file
    !
    if (nr == 1) then
      sect_name = "Edges"
    else if (nr == 2) then
      sect_name = "Quadrilaterals"
    else if (nr == 3) then
      sect_name = "Hexahedrals"
    end if
    link_path = base_name//"/"//zone_name//"/"//trim(sect_name)
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_link_write_f",override_barrier=true)
    call cg_link_write_f(trim(sect_name),trim(cgns_base_fname), &
                         trim(link_path),cgierr)
    call cgns_error(pname,cgns_time_file, &
                    "cg_link_write_f("//trim(sect_name)//")", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_link_write_f",override_barrier=true)
    !
    !
    ! Create a links to the solution node in each of the solution files
    !
   !call cg_goto_f(ifile_t,ibase_t,cgierr, &
   !               "Zone_t",izone_t, &
   !               CGNS_GOTO_END)
   !call cgns_error(pname,cgns_time_file,"cg_goto_f", &
   !                cgierr,__LINE__,__FILE__,not_parallel=true)
    !
    do n = 1,size(solname)
      !
      link_path = base_name//"/"//zone_name//"/"//trim(solname(n))
      !
      call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                      "cg_link_write_f",override_barrier=true)
      call cg_link_write_f(trim(solname(n)),trim(solfile(n)), &
                           trim(link_path),cgierr)
      call cgns_error(pname,cgns_time_file,"cg_link_write_f", &
                      cgierr,__LINE__,__FILE__,not_parallel=true)
      call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                      "cg_link_write_f",override_barrier=true)
      !
    end do
    !
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_biter_write_f",override_barrier=true)
    call cg_biter_write_f(ifile_t,ibase_t,'TimeIterValues', &
                          size(soltime,kind=CBT),cgierr)
    call cgns_error(pname,cgns_time_file,"cg_biter_write_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_biter_write_f",override_barrier=true)
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_goto_f",override_barrier=true)
    call cg_goto_f(ifile_t,ibase_t,cgierr, &
                   "BaseIterativeData_t",iiter_t, &
                   CGNS_GOTO_END)
    call cgns_error(pname,cgns_time_file,"cg_goto_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_goto_f",override_barrier=true)
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_array_write_f",override_barrier=true)
    call cg_array_write_f("TimeValues",CGNS_KIND_REAL,1_CBT, &
                          size(soltime,kind=CST),soltime,cgierr)
    call cgns_error(pname,cgns_time_file,"cg_array_write_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_array_write_f",override_barrier=true)
    !
    ! Update the zone iterative node
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_ziter_write_f",override_barrier=true)
    call cg_ziter_write_f(ifile_t,ibase_t,izone_t,"ZoneIterativeData",cgierr)
    call cgns_error(pname,cgns_time_file,"cg_ziter_write_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_ziter_write_f",override_barrier=true)
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_goto_f",override_barrier=true)
    call cg_goto_f(ifile_t,ibase_t,cgierr, &
                   "Zone_t",izone_t, &
                   "ZoneIterativeData_t",iiter_t, &
                   CGNS_GOTO_END)
    call cgns_error(pname,cgns_time_file,"cg_goto_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_goto_f",override_barrier=true)
    !
    idata(:) = int( [32,size(solname)] , kind=kind(idata) )
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_array_write_f",override_barrier=true)
    call cg_array_write_f("FlowSolutionPointers",CHARACTER,2_CBT, &
                          idata,solname,cgierr)
    call cgns_error(pname,cgns_time_file,"cg_array_write_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_array_write_f",override_barrier=true)
    !
    ! Identify the flow solution within the CGNS file as time dependent
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_simulation_type_write_f",override_barrier=true)
    call cg_simulation_type_write_f(ifile_t,ibase_t,TIMEACCURATE,cgierr)
    call cgns_error(pname,cgns_time_file,"cg_simulation_type_write_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cg_simulation_type_write_f",override_barrier=true)
    !
    ! Close the CGNS FILE
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_close_f",override_barrier=true)
    call cgp_close_f(ifile_t,cgierr)
    call cgns_error(pname,cgns_time_file,"cgp_close_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_close_f",override_barrier=true)
    !
    call cgns_timer(stop_cgns_file_timer,profile_io_cgns_subs, &
                    cgns_time_file,override_barrier=true)
    !
  end if
  !
  call cgp_mpi_comm_f(mpi_comm_world_val,cgierr)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      end_wall_time = mpi_wtime()
      cgns_wall_time = end_wall_time - beg_wall_time
      write (iout,2) pname,trim(adjustl(cgns_time_file)),cgns_wall_time
      flush (iout)
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," Entering ",a)
  2 format (1x,"Leaving ",a,/,:, &
            5x,"Wall time taken to write the file '",a, &
            "' = ",f12.6," seconds",/)
  !
end subroutine write_parallel_cgns_time_file
!
!###############################################################################
!
subroutine write_parallel_cgns_grid_file()
  !
  !.. Use Statements ..
  use geovar, only : nr,n_solpts,xyz
  use ovar,   only : profile_io_cgns,profile_io_cgns_subs
  use ovar,   only : grid_scaling_factor
  use ovar,   only : channel_nondim_output_coord
  !
  !.. Local Scalars ..
  integer :: n,ierr
  !
  integer(CST) :: my_cells_beg,my_cells_end
  integer(CST) :: n_total_output_cells
  !
  real(wp) :: beg_wall_time
  real(wp) :: end_wall_time
  real(wp) :: cgns_wall_time
  !
  character(len=CGLEN) :: varname
  !
  !.. Local Arrays ..
  integer(CBT) :: ixyz_n(1:nr)
  !
  !.. Local Allocatable Arrays ..
  integer(CST), allocatable :: ipelem(:,:)
  !
  real(CRT), allocatable :: var(:)
  real(wp),  allocatable :: vartmp(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_parallel_cgns_grid_file"
  character(len=*), parameter :: char_xyz(1:3) = ["X","Y","Z"]
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
#ifndef SPECIAL_FOR_PGI
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      write (iout,1) pname
      flush (iout)
      beg_wall_time = mpi_wtime()
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  cgns_grid_file = trim(adjustl(cgns_dir)) // cgns_grid_fname
  !
  ! Initialize the variables and data arrays needed for CGNS output
  ! NOTE: This call to routine get_initial_data_for_cgns will initialize
  !       the module variables my_pts_beg, my_pts_end, and n_my_output_pts.
  !
  call get_initial_data_for_cgns(my_cells_beg,my_cells_end,ipelem)
  !
  ! Extract the total number of cells from the zone_size array
  !
  n_total_output_cells = zone_size(2)
  !
  ! Open the cgns grid file
  !
  call cgns_timer(start_cgns_file_timer,profile_io_cgns_subs,cgns_grid_file)
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cgp_open_f")
  call cgp_open_f(cgns_grid_file,CG_MODE_WRITE,ifile_n,cgierr)
  call cgns_error(pname,cgns_grid_file,"cgp_open_f", &
                  cgierr,__LINE__,__FILE__, &
                  "Mode = CG_MODE_WRITE")
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cgp_open_f")
  !
  ! Set the CGNS base (who knows what this is for)
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_base_write_f")
  call cg_base_write_f(ifile_n,base_name,int(nr,kind=CBT), &
                       int(nr,kind=CBT),ibase_n,cgierr)
  call cgns_error(pname,cgns_grid_file,"cg_base_write_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_base_write_f")
  !
  ! Write the title for this solution
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_zone_write_f")
  call cg_zone_write_f(ifile_n,ibase_n,zone_name,zone_size, &
                       UNSTRUCTURED,izone_n,cgierr)
  call cgns_error(pname,cgns_grid_file,"cg_zone_write_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_zone_write_f")
  !
  ! Create data nodes for the coordinates
  !
  do n = 1,nr
    varname = "Coordinate"//char_xyz(n)
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_coord_write_f")
    call cgp_coord_write_f(ifile_n,ibase_n,izone_n,CGNS_KIND_REAL, &
                           trim(varname),ixyz_n(n),cgierr)
    call cgns_error(pname,cgns_grid_file,"cgp_coord_write_f", &
                    cgierr,__LINE__,__FILE__,"Write node: Coordinate",n)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_coord_write_f")
  end do
  !
  ! Now create and write the coordinate data in parallel
  !
  ! Allocate the temporary variable arrays
  !
  allocate ( vartmp(1:n_my_output_pts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"vartmp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( var(1:n_my_output_pts) , source=real(0,CRT) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"var",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Stop/Disable the CGNS queue if it is not to be used
  !
  call disable_cgns_queue(pname,cgns_grid_file)
  !
  do n = 1,nr
    !
    ! Interpolate vartmp to the desired output locations and
    ! convert to the real type required for the CGNS file.
    !
    call interpolate_for_output( xyz(n,1:n_solpts) , vartmp )
    !
    if (channel_nondim_output_coord) then
      vartmp(:) = vartmp(:) / 0.003048_wp
    end if
    !
    var(:) = real( vartmp/grid_scaling_factor , kind(var) )
    !
    ! Reset the queue.
    !
    call set_cgns_queue(BEGIN_CGNS_QUEUE,pname,cgns_grid_file)
    !
    ! Add the data for current coordinate variable to the output queue.
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_coord_write_data_f")
    call cgp_coord_write_data_f(ifile_n,ibase_n,izone_n,ixyz_n(n), &
                                my_pts_beg,my_pts_end,var,cgierr)
    call cgns_error(pname,cgns_grid_file,"cgp_coord_write_data_f", &
                    cgierr,__LINE__,__FILE__,"Write data: Coordinate",n)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_coord_write_data_f")
    !
    ! Finally, flush the queue containing the
    ! current coordinate variable to the CGNS file.
    !
    call flush_cgns_queue(pname,cgns_grid_file)
    !
  end do
  !
  ! Deallocate the var and vartmp arrays
  !
  if (allocated(var)) then
    deallocate ( var , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"var",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(vartmp)) then
    deallocate ( vartmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vartmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Create the data node for the grid connectivity
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cgp_section_write_f")
  if (nr == 1) then
    call cgp_section_write_f(ifile_n,ibase_n,izone_n,"Edges", &
                             BAR_2,1_CST,n_total_output_cells, &
                             0_CBT,isect_n,cgierr)
  else if (nr == 2) then
    call cgp_section_write_f(ifile_n,ibase_n,izone_n,"Quadrilaterals", &
                             QUAD_4,1_CST,n_total_output_cells, &
                             0_CBT,isect_n,cgierr)
  else if (nr == 3) then
    call cgp_section_write_f(ifile_n,ibase_n,izone_n,"Hexahedrals", &
                             HEXA_8,1_CST,n_total_output_cells, &
                             0_CBT,isect_n,cgierr)
  end if
  call cgns_error(pname,cgns_grid_file,"cgp_section_write_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cgp_section_write_f")
  !
  ! Write the element connectivity in parallel
  !
  ! Reset the queue.
  !
  call set_cgns_queue(BEGIN_CGNS_QUEUE,pname,cgns_grid_file)
  !
  ! Add the element connectivity data to the output queue.
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cgp_elements_write_data_f")
  call cgp_elements_write_data_f(ifile_n,ibase_n,izone_n,isect_n, &
                                 my_cells_beg,my_cells_end,ipelem,cgierr)
  call cgns_error(pname,cgns_grid_file,"cgp_elements_write_data_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cgp_elements_write_data_f")
  !
  ! Finally, flush the queue containing the
  ! element connectivity to the CGNS file.
  !
  call flush_cgns_queue(pname,cgns_grid_file)
  !
  ! Deallocate ipelem
  !
  if (allocated(ipelem)) then
    deallocate ( ipelem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ipelem",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Close the CGNS FILE
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cgp_close_f")
  call cgp_close_f(ifile_n,cgierr)
  call cgns_error(pname,cgns_grid_file,"cgp_close_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cgp_close_f")
  !
  call cgns_timer(stop_cgns_file_timer,profile_io_cgns_subs,cgns_grid_file)
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      end_wall_time = mpi_wtime()
      cgns_wall_time = end_wall_time - beg_wall_time
      write (iout,2) pname,trim(adjustl(cgns_grid_file)),cgns_wall_time
      flush (iout)
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," Entering ",a)
  2 format (1x,"Leaving ",a,/,:, &
            5x,"Wall time taken to write the file '",a, &
            "' = ",f12.6," seconds",/)
  !
end subroutine write_parallel_cgns_grid_file
!
!###############################################################################
!
subroutine write_parallel_cgns_solution_file(sol_file,tname)
  !
  !.. Use Statements ..
  use geovar,  only : nr,n_solpts
  use ovar,    only : governing_equations
  use ovar,    only : profile_io_cgns,profile_io_cgns_subs
  use ovar,    only : max_wall_time => max_cgns_wall_time
  use flowvar, only : usp,dusp,gpsp
  !
  use metrics_mod, only : time_metrics
  !
  use thermo_mod,    only : species_name
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  use flux_mod, only : discontinuous_gradient_at_solpts
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: sol_file
  character(len=*), intent(in) :: tname
  !
  !.. Local Scalars ..
  integer :: i,k,l,m,n,nc,npts,ierr
  !
  real(wp) :: beg_wall_time
  real(wp) :: end_wall_time
  real(wp) :: cgns_wall_time
  !
  character(len=CGLEN) :: varname
  character(len=CGLEN) :: sect_name
  character(len=170)   :: cgns_sol_file
  character(len=300)   :: link_path
  !
  !.. Local Allocatable Arrays ..
#ifdef SPECIAL_FOR_INTEL
  character(len=:), allocatable :: char_sname(:)
#else
  character(len=CGLEN), allocatable :: char_sname(:)
#endif
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_parallel_cgns_solution_file"
  character(len=*), parameter :: char_xyz(1:3) = ["X","Y","Z"]
  character(len=*), parameter :: char_ijk(1:3) = ["Xi  ","Eta ","Zeta"]
  !
  !.. Local Scalars ..
  integer(CBT) :: ifile_s,ibase_s,izone_s,isol_s
  !
  !.. Local Allocatable Arrays ..
  integer(CBT), allocatable, dimension(:) :: isvar_s
  !
  real(CRT), allocatable :: var(:)
  real(wp),  allocatable :: vartmp(:)
  real(wp),  allocatable :: interp_usp(:,:)
  real(wp),  allocatable :: interp_gpsp(:,:)
  real(wp),  allocatable :: interp_dusp(:,:,:)
  real(wp),  allocatable :: dusp_copy(:,:,:)
  !
 !logical(lk), parameter :: output_discontinuous_gradients = true
  logical(lk), parameter :: output_discontinuous_gradients = fals
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
#ifndef SPECIAL_FOR_PGI
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      write (iout,1) pname
      flush (iout)
      beg_wall_time = mpi_wtime()
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Get the list of variables that will be written to the CGNS solution file
  !
  call get_cgns_sol_variables(char_sname)
  !
  ! Get the name of the current solution file
  !
  cgns_sol_file = trim(adjustl(cgns_dir)) // trim(adjustl(sol_file))
  !
  ! Open the CGNS solution file if iflag is 1
  !
  call cgns_timer(start_cgns_file_timer,profile_io_cgns_subs,cgns_sol_file)
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cgp_open_f")
  call cgp_open_f(cgns_sol_file,CG_MODE_WRITE,ifile_s,cgierr)
  call io_error(pname,cgns_sol_file,1,__LINE__,__FILE__,cgierr)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cgp_open_f")
  !
  ! Set the CGNS base (who knows what this is for)
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_base_write_f")
  call cg_base_write_f(ifile_s,base_name,int(nr,kind=CBT), &
                       int(nr,kind=CBT),ibase_s,cgierr)
  call cgns_error(pname,cgns_sol_file,"cg_base_write_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_base_write_f")
  !
  ! Write the title for this solution
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_zone_write_f")
  call cg_zone_write_f(ifile_s,ibase_s,zone_name,zone_size, &
                       UNSTRUCTURED,izone_s,cgierr)
  call cgns_error(pname,cgns_sol_file,"cg_zone_write_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_zone_write_f")
  !
  ! Create a link to the grid coordinates in the grid file
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cg_goto_f")
  call cg_goto_f(ifile_s,ibase_s,cgierr, &
                 "Zone_t",izone_s, &
                 CGNS_GOTO_END)
  call cgns_error(pname,cgns_sol_file,"cg_goto_f",cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cg_goto_f")
  !
  link_path = base_name//"/"//zone_name//"/GridCoordinates"
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  call cg_link_write_f("GridCoordinates",trim(cgns_base_fname), &
                       trim(link_path),cgierr)
  call cgns_error(pname,cgns_sol_file,"cg_link_write_f(GridCoordinates)", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  !
  ! Create a link to the element connectivity in the grid file
  !
  if (nr == 1) then
    sect_name = "Edges"
  else if (nr == 2) then
    sect_name = "Quadrilaterals"
  else if (nr == 3) then
    sect_name = "Hexahedrals"
  end if
  link_path = base_name//"/"//zone_name//"/"//trim(sect_name)
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  call cg_link_write_f(trim(sect_name),trim(cgns_base_fname), &
                       trim(link_path),cgierr)
  call cgns_error(pname,cgns_sol_file, &
                  "cg_link_write_f("//trim(sect_name)//")", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  !
  ! Allocate var and vartmp if they arent already allocated
  !
  if (allocated(vartmp)) then
    deallocate ( vartmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vartmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(interp_usp)) then
    deallocate ( interp_usp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interp_usp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(interp_gpsp)) then
    deallocate ( interp_gpsp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interp_gpsp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(interp_dusp)) then
    deallocate ( interp_dusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interp_dusp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (.not. allocated(var)) then
    allocate ( var(1:n_my_output_pts) , source=real(0,CRT) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"var",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Write the flow variables
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_sol_write_f")
  call cg_sol_write_f(ifile_s,ibase_s,izone_s,tname,VERTEX,isol_s,cgierr)
  call cgns_error(pname,cgns_sol_file,"cg_sol_write_f",cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_sol_write_f")
  !
  ! #######################################
  ! #####   WRITING FIELD VARIABLES   #####
  ! #######################################
  !
  npts = n_my_output_pts
  !
  allocate ( vartmp(1:npts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"vartmp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( interp_usp(1:npts,1:nq) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"interp_usp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( interp_dusp(1:npts,1:nr,1:nq) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"interp_dusp",1,__LINE__,__FILE__,ierr,error_message)
  !
  !
  !
  call interpolate_for_output( usp , interp_usp )
  !
  if (allocated(dusp)) then
    !
    if (output_discontinuous_gradients) then
      !
      ! Calling discontinuous_gradient_at_solpts will replace the
      ! contents of dusp with the discontinuous gradient, so save
      ! the contents of dusp in dusp_copy.
      !
      dusp_copy = dusp ! F2003 AUTOMATIC-REALLOCATION
      !
      call discontinuous_gradient_at_solpts()
      !
    end if
    !
    call interpolate_for_output( dusp , interp_dusp )
    !
    if (output_discontinuous_gradients) then
      !
      ! Now that we have the discontinuous gradients in interp_dusp,
      ! copy the original contents of dusp that were saved in dusp_copy
      ! back into dusp, then deallocate dusp_copy.
      !
      dusp(:,:,:) = dusp_copy(:,:,:)
      !
      deallocate ( dusp_copy , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"dusp_copy",2,__LINE__,__FILE__,ierr,error_message)
      !
    end if
    !
  end if
  !
  if (check_gradients) then
    !
    ! ################################################
    ! #####   WRITING GRADIENT FIELD VARIABLES   #####
    ! ################################################
    !
    ! Count the total number of solution variables that will be written
    ! and allocate isvar_s accordingly
    !
    n = nr * (nme - nec)
    allocate ( isvar_s(1:n) , source=0 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"isvar_s",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Create data nodes for the solution variables
    !
    ! Density Gradients
    !
    do l = 1,nr
      varname = "Density_d"//char_xyz(l)
      call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                      "cgp_field_write_f")
      call cgp_field_write_f(ifile_s,ibase_s,izone_s,isol_s,CGNS_KIND_REAL, &
                             trim(varname),isvar_s(l),cgierr)
      call cgns_error(pname,cgns_sol_file,"cgp_field_write_f", &
                      cgierr,__LINE__,__FILE__, &
                      "Write node: Gradient variable",l)
      call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                      "cgp_field_write_f")
    end do
    !
    ! Velocity Gradients
    !
    do m = 1,nr
      do l = 1,nr
        n = m*nr + l
        varname = "Velocity"//char_xyz(m)//"_d"//char_xyz(l)
        call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                        "cgp_field_write_f")
        call cgp_field_write_f(ifile_s,ibase_s,izone_s,isol_s,CGNS_KIND_REAL, &
                               trim(varname),isvar_s(n),cgierr)
        call cgns_error(pname,cgns_sol_file,"cgp_field_write_f", &
                        cgierr,__LINE__,__FILE__, &
                        "Write node: Gradient variable",n)
        call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                        "cgp_field_write_f")
      end do
    end do
    !
    ! Now create the output data and write it to the CGNS file in parallel
    !
    grad_eqn: do m = nec,nme
      !
      grad_dir: do l = 1,nr
        !
        n = (m-nec)*nr + l
        !
        vartmp(1:npts) = grad_of_selected_pv(l,m,interp_usp(1:npts,:), &
                                                interp_dusp(1:npts,:,:))
        !
        var(:) = real( vartmp , kind(var) )
        !
        ! Reset the queue.
        !
        call set_cgns_queue(BEGIN_CGNS_QUEUE,pname,cgns_sol_file)
        !
        ! Add the data for current variable to the output queue.
        !
        call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                        "cgp_field_write_data_f")
        call cgp_field_write_data_f(ifile_s,ibase_s,izone_s,isol_s,isvar_s(n), &
                                    my_pts_beg,my_pts_end,var,cgierr)
        call cgns_error(pname,cgns_sol_file,"cgp_field_write_data_f", &
                        cgierr,__LINE__,__FILE__, &
                        "Write data: Gradient variable",n)
        call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                        "cgp_field_write_data_f")
        !
        ! Finally, flush the queue containing
        ! the current variable to the CGNS file.
        !
        call flush_cgns_queue(pname,cgns_sol_file)
        !
      end do grad_dir
      !
    end do grad_eqn
    !
  else
    !
    ! ################################################
    ! #####   WRITING SOLUTION FIELD VARIABLES   #####
    ! ################################################
    !
    ! Count the total number of solution variables that will be written
    ! and allocate isvar_s accordingly
    !
    n = 0
    do l = 1,size(char_sname)
      if (trim(char_sname(l)) == "Velocity") then
        n = n + nr
      else if (trim(char_sname(l)) == "Vorticity") then
        n = n + merge(nr,1,nr==3)
      else if (trim(char_sname(l)) == "Momentum") then
        n = n + nr
      else if (trim(char_sname(l)) == "GradDensity") then
        n = n + nr
      else if (trim(char_sname(l)) == "GradMomentum") then
        n = n + nr*nr
      else if (trim(char_sname(l)) == "GradVelocity") then
        n = n + nr*nr
      else if (trim(char_sname(l)) == "GradEnergy") then
        n = n + nr
      else if (trim(char_sname(l)) == "GradTemperature") then
        n = n + nr
      else if (trim(char_sname(l)) == "DmatTemperature") then
        n = n + nr
      else if (trim(char_sname(l)) == "GradMassFractions") then
        n = n + ncs*nr
      else if (trim(char_sname(l)) == "GradSpeciesDensities") then
        n = n + ncs*nr
      else if (trim(char_sname(l)) == "Viscosity") then
        if (governing_equations == NavierStokes_Eqns) then
          n = n + 1
        end if
      else if (trim(char_sname(l)) == "ResolutionMarker") then
        n = n + 2
      else if (trim(char_sname(l)) == "JumpMarker") then
        n = n + 2
      else if (trim(char_sname(l)) == "TimeMetrics") then
        if (allocated(time_metrics)) then
          n = n + nr + 1
        end if
      else if (trim(char_sname(l)) == "CoMetrics") then
     !  n = n + nr*nr
      else if (trim(char_sname(l)) == "ContraMetrics") then
        n = n + nr*nr
      else if (trim(char_sname(l)) == "MassFractions") then
        if (ncs > 0) then
          n = n + ncs
        end if
      else if (trim(char_sname(l)) == "MoleFractions") then
        if (ncs > 0) then
          n = n + ncs
        end if
      else if (trim(char_sname(l)) == "SpeciesDensities") then
        if (ncs > 0) then
          n = n + ncs
        end if
      else if (trim(char_sname(l)) == "SpeciesEnthalpy") then
        if (ncs > 0) then
          n = n + ncs
        end if
      else if (trim(char_sname(l)) == "SpeciesCp") then
        if (ncs > 0) then
          n = n + ncs
        end if
      else if (trim(char_sname(l)) == "cpi/R") then
        if (ncs > 0) then
          n = n + ncs
        end if
      else if (trim(char_sname(l)) == "hi/RT") then
        if (ncs > 0) then
          n = n + ncs
        end if
      else if (trim(char_sname(l)) == "hi/R") then
        if (ncs > 0) then
          n = n + ncs
        end if
      else if (trim(char_sname(l)) == "dTdx_Breakdown") then
        if (ncs > 0) then
          n = n + (2*nr + 2*ncs + 7 + ncs*nr)
        end if
      else if (trim(char_sname(l)) == "QCriterion") then
        n = n + 2
      else
        n = n + 1
      end if
    end do
    !
    allocate ( isvar_s(1:n) , source=0 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"isvar_s",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Create data nodes for the solution variables
    !
    n = 0
    do k = 1,size(char_sname)
      if (trim(char_sname(k)) == "Velocity") then
        do m = 1,nr
          n = n + 1
          varname = "Velocity"//char_xyz(m)
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      else if (trim(char_sname(k)) == "Vorticity") then
        do m = 1,nr
          n = n + 1
          varname = "Vorticity"//char_xyz(merge(m,3,nr==3))
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          if (nr < 3) exit
        end do
      else if (trim(char_sname(k)) == "Momentum") then
        do m = 1,nr
          n = n + 1
          varname = "Momentum"//char_xyz(m)
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      else if (trim(char_sname(k)) == "GradDensity") then
        do l = 1,nr
          n = n + 1
          varname = "DensityGradient"//char_xyz(l)
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      else if (trim(char_sname(k)) == "GradMomentum") then
        do m = 1,nr
          do l = 1,nr
            n = n + 1
            varname = "Momentum"//char_xyz(m)//"Gradient"//char_xyz(l)
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
        end do
      else if (trim(char_sname(k)) == "GradVelocity") then
        do m = 1,nr
          do l = 1,nr
            n = n + 1
            varname = "Velocity"//char_xyz(m)//"Gradient"//char_xyz(l)
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
        end do
      else if (trim(char_sname(k)) == "GradEnergy") then
        do l = 1,nr
          n = n + 1
          varname = "EnergyGradient"//char_xyz(l)
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      else if (trim(char_sname(k)) == "GradTemperature") then
        do l = 1,nr
          n = n + 1
          varname = "TemperatureGradient"//char_xyz(l)
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      else if (trim(char_sname(k)) == "DmatTemperature") then
        do l = 1,nr
          n = n + 1
          varname = "TemperatureDMatrix"//char_xyz(l)
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end do
      else if (trim(char_sname(k)) == "GradMassFractions") then
        do m = 1,ncs
          do l = 1,nr
            n = n + 1
            varname = trim(adjustl(species_name(m))) // &
                      "MassFractionGradient"//char_xyz(l)
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
        end do
      else if (trim(char_sname(k)) == "GradSpeciesDensities") then
        do m = 1,ncs
          do l = 1,nr
            n = n + 1
            varname = trim(adjustl(species_name(m))) // &
                      "DensityGradient"//char_xyz(l)
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
        end do
      else if (trim(char_sname(k)) == "Viscosity") then
        if (governing_equations == NavierStokes_Eqns) then
          n = n + 1
          varname = trim(char_sname(k))
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end if
      else if (trim(char_sname(k)) == "ResolutionMarker") then
        !
        n = n + 1
        varname = trim(char_sname(k))
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        !
        n = n + 1
        varname = "Scaled"//trim(char_sname(k))
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        !
      else if (trim(char_sname(k)) == "JumpMarker") then
        !
        n = n + 1
        varname = trim(char_sname(k))
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        !
        n = n + 1
        varname = "Scaled"//trim(char_sname(k))
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        !
      else if (trim(char_sname(k)) == "CoMetrics") then
     !  do m = 1,nr
     !    do i = 1,nr
     !      n = n + 1
     !      varname = "Covariant Metric d" // &
     !                trim(adjustl(char_xyz(m))) // "d" // char_ijk(i)
     !      call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
     !    end do
     !  end do
      else if (trim(char_sname(k)) == "ContraMetrics") then
        do m = 1,nr
          do i = 1,nr
            n = n + 1
            varname = "Contravariant Metric d" // &
                      trim(adjustl(char_ijk(m))) // "d" // char_xyz(i)
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
        end do
      else if (trim(char_sname(k)) == "TimeMetrics") then
        if (allocated(time_metrics)) then
          do m = 1,nr
            n = n + 1
            varname = "Time Metric d" // char_xyz(m)
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
          n = n + 1
          varname = "Time Metric dA"
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        end if
      else if (trim(char_sname(k)) == "MassFractions") then
        if (ncs > 0) then
          do m = 1,ncs
            n = n + 1
            varname = "MassFraction"//trim(adjustl(species_name(m)))
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
        end if
      else if (trim(char_sname(k)) == "MoleFractions") then
        if (ncs > 0) then
          do m = 1,ncs
            n = n + 1
            varname = "MoleFraction"//trim(adjustl(species_name(m)))
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
        end if
      else if (trim(char_sname(k)) == "SpeciesDensities") then
        if (ncs > 0) then
          do m = 1,ncs
            n = n + 1
            varname = "SpeciesDensity"//trim(adjustl(species_name(m)))
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
        end if
      else if (trim(char_sname(k)) == "SpeciesEnthalpy") then
        if (ncs > 0) then
          do m = 1,ncs
            n = n + 1
            varname = "SpeciesEnthalpy"//trim(adjustl(species_name(m)))
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
        end if
      else if (trim(char_sname(k)) == "SpeciesCp") then
        if (ncs > 0) then
          do m = 1,ncs
            n = n + 1
            varname = "SpeciesCp"//trim(adjustl(species_name(m)))
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
        end if
      else if (trim(char_sname(k)) == "cpi/R") then
        if (ncs > 0) then
          do m = 1,ncs
            n = n + 1
            varname = "cpi_over_R_"//trim(adjustl(species_name(m)))
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
        end if
      else if (trim(char_sname(k)) == "hi/R") then
        if (ncs > 0) then
          do m = 1,ncs
            n = n + 1
            varname = "hi_over_R_"//trim(adjustl(species_name(m)))
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
        end if
      else if (trim(char_sname(k)) == "hi/RT") then
        if (ncs > 0) then
          do m = 1,ncs
            n = n + 1
            varname = "hi_over_RT_"//trim(adjustl(species_name(m)))
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
        end if
      else if (trim(char_sname(k)) == "dTdx_Breakdown") then
        if (ncs > 0) then
          do l = 1,nr
            n = n + 1
            varname = "CPG_dTd"//char_xyz(l)
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
          do m = 1,ncs
            n = n + 1
            varname = "Cv("//trim(adjustl(species_name(m)))//")"
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
          do m = 1,ncs
            n = n + 1
            varname = "dCp("//trim(adjustl(species_name(m)))//")dT"
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
          !
          n = n + 1
          varname = "dCp(mix)dT"
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          !
          n = n + 1
          varname = "T_from_gradcv2gradpv"
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          !
          n = n + 1
          varname = "Cp(mix)"
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          !
          n = n + 1
          varname = "Gamma(Mix)"
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          !
          n = n + 1
          varname = "dTdx_lead_coefficient"
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          !
          n = n + 1
          varname = "CPG_constant_Cv"
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          !
          n = n + 1
          varname = "CPG_1_over_constant_Cv"
          call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          !
          do l = 1,nr
            n = n + 1
            varname = "TPG_dTd"//char_xyz(l)
            call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
          end do
          !
          do m = 1,ncs
            do l = 1,nr
              n = n + 1
              varname = "dY("//trim(adjustl(species_name(m)))// &
                        ")d"//char_xyz(l)//"_from_dTd"//char_xyz(l)
              call cgns_write_field_header(n,varname,isvar_s(n), &
                                           __LINE__,__FILE__)
            end do
          end do
        end if
      else if (trim(char_sname(k)) == "QCriterion") then
        !
        n = n + 1
        varname = trim(char_sname(k))
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        !
        n = n + 1
        varname = trim(char_sname(k))//"_Positive"
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
        !
      else
        n = n + 1
        if (trim(char_sname(k)) == "Energy") then
          varname = "EnergyStagnationDensity"
        else if (trim(char_sname(k)) == "T0") then
          varname = "TemperatureStagnation"
        else if (trim(char_sname(k)) == "P0") then
          varname = "PressureStagnation"
        else if (trim(char_sname(k)) == "H0") then
          varname = "EnthalpyStagnation"
        else if (trim(char_sname(k)) == "Viscosity") then
          varname = "ViscosityMolecular"
        else if (trim(char_sname(k)) == "Conductivity") then
          varname = "ThermalConductivity"
        else if (trim(char_sname(k)) == "Cp") then
          varname = "SpecificHeatPressure"
        else if (trim(char_sname(k)) == "Gamma") then
          varname = "SpecificHeatRatio"
        else if (trim(char_sname(k)) == "Rgas") then
          varname = "IdealGasConstant"
        else
          varname = trim(char_sname(k))
        end if
        call cgns_write_field_header(n,varname,isvar_s(n),__LINE__,__FILE__)
      end if
    end do
    !
    !
    !
    allocate ( interp_gpsp(1:npts,1:ng) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interp_gpsp",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Initialize interp_gpsp to the interpolated values of gpsp
    ! Interpolate the values in gpsp to initialize interp_gpsp. We need to
    ! initialize the temperature in interp_gpsp in case this is for a
    ! thermally perfect gas mixture which needs an initial guess for
    ! temperature when performing the Newton solve for temperature.
    !
    call interpolate_for_output( gpsp , interp_gpsp )
    !
    ! Now actually compute the updated thermodynamic and transport properties
    ! based on the interpolated values of the conserved variables.
    !
    call update_thermo_properties(1,npts,interp_usp,interp_gpsp, &
                                  caller_name=pname)
    call update_transport_properties(1,npts,interp_usp,interp_gpsp)
    !
    ! Now create the output data and write it to the CGNS file in parallel
    !
    n = 0
    do l = 1,size(char_sname)
      !
      call output_variable(ifile_s,ibase_s,izone_s,isol_s,isvar_s, &
                           cgns_sol_file,char_sname(l),n, &
                           interp_usp,interp_gpsp,interp_dusp,vartmp,var)
      !
    end do
    !
  end if
  !
  ! ################################################
  ! #####   FINISHED WRITING FIELD VARIABLES   #####
  ! ################################################
  !
  call set_cgns_queue(STOP_CGNS_QUEUE,pname,cgns_sol_file)
  !
  ! Close the CGNS FILE
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cgp_close_f")
  call cgp_close_f(ifile_s,cgierr)
  call io_error(pname,cgns_sol_file,2,__LINE__,__FILE__,cgierr)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cgp_close_f")
  !
  call cgns_timer(stop_cgns_file_timer,profile_io_cgns_subs,cgns_sol_file)
  !
  ! Deallocate the arrays var and vartmp
  !
  if (allocated(var)) then
    deallocate ( var , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"var",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(vartmp)) then
    deallocate ( vartmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vartmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(interp_usp)) then
    deallocate ( interp_usp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interp_usp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(interp_dusp)) then
    deallocate ( interp_dusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interp_dusp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(isvar_s)) then
    deallocate ( isvar_s , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"isvar_s",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      end_wall_time = mpi_wtime()
      cgns_wall_time = end_wall_time - beg_wall_time
      max_wall_time = max(max_wall_time,cgns_wall_time)
      write (iout,2) pname,trim(adjustl(cgns_sol_file)),cgns_wall_time
      flush (iout)
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," Entering ",a)
  2 format (1x,"Leaving ",a,/,:, &
            5x,"Wall time taken to write the file '",a, &
            "' = ",f12.6," seconds",/)
  !
  !-----------------------------------------------------------------------------
  !
  contains
  !
#include "Functions/cgns_write_field_header.f90"
  !
  !-----------------------------------------------------------------------------
  !
end subroutine write_parallel_cgns_solution_file
!
!###############################################################################
!
subroutine output_variable(ifile_s,ibase_s,izone_s,isol_s,isvar_s, &
                           cgns_sol_file,vname,nv, &
                           cv,gp,dcvdx,vartmp,var)
  !
  !.. Use Statements ..
  use geovar,  only : nr,n_solpts,ncell,cell
  use ovar,    only : aref,rhoref,time_ref,tref
  use ovar,    only : gam,pref,pref_nd,machref
  use ovar,    only : rgasref,cpref,Pr
  use ovar,    only : itestcase,governing_equations
  use ovar,    only : profile_io_cgns_subs
  use ovar,    only : Limiter_Option
  use flowvar, only : dtsp
  !
  use metrics_mod, only : time_metrics,metrics_array,jac
  !
  use module_limiters, only : marker,scaled_marker
  use module_limiters, only : resolution_indicator
  use module_limiters, only : jump_indicator
  !
  use thermo_mod, only : rspec
  use thermo_mod, only : compute_mixture_rgas
  use thermo_mod, only : compute_species_enthalpy
  use thermo_mod, only : compute_species_cp
  !
  use transport_mod, only : compute_mole_fractions
  !
  use gpu_mod, only : grad_cv_to_grad_pv_gpu
  !
  !.. Formal Arguments ..
  integer(CBT),     intent(in) :: ifile_s
  integer(CBT),     intent(in) :: ibase_s
  integer(CBT),     intent(in) :: izone_s
  integer(CBT),     intent(in) :: isol_s
  integer(CBT),     intent(in) :: isvar_s(:)
  character(len=*), intent(in) :: cgns_sol_file
  character(len=*), intent(in) :: vname
  integer,       intent(inout) :: nv
  real(wp),      intent(inout) :: cv(:,:)
  real(wp),      intent(inout) :: gp(:,:)
  real(wp),         intent(in) :: dcvdx(:,:,:)
  real(wp),      intent(inout) :: vartmp(:)
  real(CRT),     intent(inout) :: var(:)
  !
  !.. Local Scalars ..
  integer  :: i,j,l,m,n,nc,np,n0,npts,ierr
  integer  :: nmi,nmj,l_comp,l_phys,nvar
  integer  :: this_order,limit_order
  real(wp) :: freestream_ref,rhoa,rhoa2,time_ref_inv
  real(wp) :: cp_dim,cond_dim,dim_val
  real(wp) :: rho_inv,vx,vy,vz,a2,q2,s
  real(wp) :: gamm,gm1,gam_over_gm1,mach,pres,temp
  real(wp) :: mi,mj,drhodi,drhodj,dmjdi,dmidj
  real(wp) :: w12,w13,w23,wijwij
  real(wp) :: s12,s13,s23,sijsij,s11,s22,s33,skk
  real(wp) :: dudx,dudy,dudz
  real(wp) :: dvdx,dvdy,dvdz
  real(wp) :: dwdx,dwdy,dwdz
  real(wp) :: qmin,qmax
  !
  integer(CBT) :: ivar_s
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: vartmp_at_sp(:)
  real(wp), allocatable :: spcvar(:,:)
  real(wp), allocatable :: dtdx(:,:)
  real(wp), allocatable :: dpvdx(:,:,:)
  real(wp), allocatable :: other_variables(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "output_variable"
  !
continue
  !
  npts = size(vartmp)
  !
  freestream_ref = zero
  rhoa = rhoref*aref
  rhoa2 = rhoa*aref
  cp_dim = aref*aref/tref
  cond_dim = rhoa*cp_dim
  !
  ! Pre-load vartmp with Mach number if the current variable
  ! is Mach number, total pressure, or total temperature.
  !
  select case (trim(adjustl(vname)))
    !
    case ("Mach","P0","T0")
      !
      if (ncs == 0) then
        do n = 1,npts
          gp(n,ngamm) = gam
        end do
      end if
      !
      if (nr == 3) then
        do n = 1,npts
          rho_inv = one/cv(n,nec)
          vx = cv(n,nmx)*rho_inv
          vy = cv(n,nmy)*rho_inv
          vz = cv(n,nmz)*rho_inv
          a2 = max( gp(n,ngamm)*gp(n,npres)*rho_inv , rndoff )
          q2 = vx*vx + vy*vy + vz*vz
          vartmp(n) = sqrt( q2 / a2 ) - freestream_ref
        end do
      else if (nr == 2) then
        do n = 1,npts
          rho_inv = one/cv(n,nec)
          vx = cv(n,nmx)*rho_inv
          vy = cv(n,nmy)*rho_inv
          a2 = max( gp(n,ngamm)*gp(n,npres)*rho_inv , rndoff )
          q2 = vx*vx + vy*vy
          vartmp(n) = sqrt( q2 / a2 ) - freestream_ref
        end do
      else
        do n = 1,npts
          rho_inv = one/cv(n,nec)
          vx = cv(n,nmx)*rho_inv
          a2 = max( gp(n,ngamm)*gp(n,npres)*rho_inv , rndoff )
          vartmp(n) = sqrt( vx*vx / a2 ) - freestream_ref
        end do
      end if
      !
    case ("GradVelocity","GradTemperature","GradMassFractions", &
          "GradTempMagnitude","VorticityMagnitude","StrainRateMagnitude", &
          "QCriterion","DivergenceOfVelocity")
      !
      allocate ( dpvdx(1:npts,1:nr,1:nq) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"dpvdx",1,__LINE__,__FILE__,ierr,error_message)
      !
      call grad_cv_to_grad_pv_gpu(1,npts,cv,gp,dcvdx,dpvdx)
      !
    case ("dTdx_Breakdown")
      !
      if (ncs > 0) then
        !
        allocate ( dpvdx(1:npts,1:nr,1:nq) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"dpvdx",1,__LINE__,__FILE__,ierr,error_message)
        !
        nvar = nr + 2*ncs + 7
        allocate ( other_variables(1:npts,1:nvar) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"other_variables",1,__LINE__,__FILE__, &
                         ierr,error_message)
        !
        call grad_cv_to_grad_pv_gpu(1,npts,cv,gp,dcvdx,dpvdx, &
                                    output_variables=other_variables)
        !
      end if
      !
    case ("SpeciesEnthalpy","SpeciesCp","hi/RT","hi/R","cpi/R")
      !
      if (ncs > 0) then
        !
        allocate ( spcvar(1:npts,1:ncs) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"spcvar",1,__LINE__,__FILE__,ierr,error_message)
        !
        select case (trim(adjustl(vname)))
          case ("SpeciesCp","cpi/R")
            call compute_species_cp(1,npts,gp(1:npts,ntemp),spcvar, &
                                    dimensional=true)
          case default
            call compute_species_enthalpy(1,npts,gp(1:npts,ntemp),spcvar, &
                                          dimensional=true)
        end select
        !
      end if
      !
    case ("MoleFractions")
      !
      if (ncs > 0) then
        !
        allocate ( spcvar(1:npts,1:ncs) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"spcvar",1,__LINE__,__FILE__,ierr,error_message)
        !
        call compute_mole_fractions(1,npts,cv,spcvar)
        !
      end if
      !
  end select
  !
  ! Compute the current output variable and write it to the CGNS file
  !
  select case (trim(adjustl(vname)))
    !
    ! #####################
    ! ###  NEED TO ADD  ###
    ! #####################
    ! ---------------------
    !   FOR GAS MIXTURES
    ! ---------------------
    ! sensible enthalpy
    ! ---------------------
    !   FOR CPG AND TPG
    ! ---------------------
    ! pitot pressure
    ! strain rate
    ! pressure gradient
    ! magnitude of pressure gradient
    !
    case ("Density") ! Density
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (itestcase == Freestream_Preservation) freestream_ref = rhoref
      !
      do n = 1,npts
        vartmp(n) = cv(n,nec) * rhoref - freestream_ref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("Momentum") ! Momentum
      !!{{{
      do m = nmb,nme
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,npts
          vartmp(n) = cv(n,m) * rhoa
        end do
        !
        call cgns_write_field_data
        !
      end do
      !!}}}
    case ("Energy") ! Energy
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,npts
        vartmp(n) = cv(n,nee) * rhoa2
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("SpeciesDensities") ! Species densities
      !!{{{
      do m = 1,ncs
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,npts
          vartmp(n) = cv(n,m)
        end do
        !
        call cgns_write_field_data
        !
      end do
      !!}}}
    case ("MassFractions") ! Species mass fractions
      !!{{{
      do m = 1,ncs
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,npts
          vartmp(n) = cv(n,m)/cv(n,nec)
        end do
        !
        call cgns_write_field_data
        !
      end do
      !!}}}
    case ("MoleFractions") ! Species mole fractions
      !!{{{
      do m = 1,ncs
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        vartmp(1:npts) = spcvar(1:npts,m)
        !
        call cgns_write_field_data
        !
      end do
      !!}}}
    case ("GradSpeciesDensities") ! Species densities gradients
      !!{{{
      do m = 1,ncs
        do l = 1,nr
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          do n = 1,npts
            vartmp(n) = dcvdx(n,l,m) * rhoref
          end do
          !
          call cgns_write_field_data
          !
        end do
      end do
      !!}}}
    case ("GradDensity") ! Density gradient
      !!{{{
      do l = 1,nr
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,npts
          vartmp(n) = dcvdx(n,l,nec) * rhoref
        end do
        !
        call cgns_write_field_data
        !
      end do
      !!}}}
    case ("GradDensityMagnitude") ! Magnitude of density gradient
      !!{{{
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,npts
        vartmp(n) = norm2( dcvdx(n,1:nr,nec) ) * rhoref
      end do
      !
      call cgns_write_field_data
      !
      !!}}}
    case ("GradMomentum") ! Momentum gradients
      !!{{{
      do m = nmb,nme
        do l = 1,nr
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          do n = 1,npts
            vartmp(n) = dcvdx(n,l,m) * rhoa
          end do
          !
          call cgns_write_field_data
          !
        end do
      end do
      !!}}}
    case ("GradEnergy") ! Total energy gradient
      !!{{{
      do l = 1,nr
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,npts
          vartmp(n) = dcvdx(n,l,nee) * rhoa2
        end do
        !
        call cgns_write_field_data
        !
      end do
      !!}}}
    case ("GradMassFractions") ! Species mass fraction gradients
      !!{{{
      do m = 1,ncs
        do l = 1,nr
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          do n = 1,npts
            vartmp(n) = dpvdx(n,l,m)
          end do
          !
          call cgns_write_field_data
          !
        end do
      end do
      !!}}}
    case ("GradVelocity") ! Velocity gradients
      !!{{{
      do m = nmb,nme
        do l = 1,nr
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          do n = 1,npts
            vartmp(n) = dpvdx(n,l,m) * aref
          end do
          !
          call cgns_write_field_data
          !
        end do
      end do
      !!}}}
    case ("GradTemperature") ! Temperature gradient
      !!{{{
      do l = 1,nr
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,npts
          vartmp(n) = dpvdx(n,l,nee) * tref
        end do
        !
        call cgns_write_field_data
        !
      end do
      !!}}}
    case ("GradTempMagnitude") ! Magnitude of temperature gradient
      !!{{{
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,npts
        vartmp(n) = norm2( dpvdx(n,1:nr,nee) ) * tref
      end do
      !
      call cgns_write_field_data
      !
      !!}}}
    case ("dTdx_Breakdown") ! Breakdown of terms to compute dT/dx
      !!{{{
      if (ncs > 0) then
        do m = 1,nvar
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          do n = 1,npts
            vartmp(n) = other_variables(n,m)
          end do
          !
          call cgns_write_field_data
          !
        end do
        !
        do l = 1,nr
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          do n = 1,npts
            vartmp(n) = dpvdx(n,l,nee)
          end do
          !
          call cgns_write_field_data
          !
        end do
        !
        do m = 1,ncs
          do l = 1,nr
            !
            nv = nv + 1
            ivar_s = isvar_s(nv)
            !
            do n = 1,npts
              vartmp(n) = dpvdx(n,l,m)
            end do
            !
            call cgns_write_field_data
            !
          end do
        end do
      end if
      !!}}}
    case ("DmatTemperature") ! Temperature gradient using derivative matrices
      !!{{{
      !
      allocate ( dtdx(1:n_solpts,1:nr) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"dtdx",1,__LINE__,__FILE__,ierr,error_message)
      !
      call compute_dmatrix_grad_temperature(dtdx)
      !
      if (.not. allocated(vartmp_at_sp)) then
        allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
      !
      do l = 1,nr
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,n_solpts
          vartmp_at_sp(n) = dtdx(n,l) * tref
        end do
        !
        call interpolate_for_output( vartmp_at_sp , vartmp )
        !
        call cgns_write_field_data
        !
      end do
      !
      if (allocated(dtdx)) then
        deallocate ( dtdx , stat=ierr , errmsg=error_message )
        call alloc_error(pname,"dtdx",2,__LINE__,__FILE__,ierr,error_message)
      end if
      !
      !!}}}
    case ("Velocity") ! Velocity
      !!{{{
      do m = nmb,nme
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        if (itestcase == Freestream_Preservation) then
          freestream_ref = merge( aref*machref , zero , m==nmb )
        end if
        !
        do n = 1,npts
          vartmp(n) = cv(n,m)/cv(n,nec) * aref - freestream_ref
        end do
        !
        call cgns_write_field_data
        !
      end do
      !!}}}
    case ("Pressure") ! Pressure
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (itestcase == Freestream_Preservation) freestream_ref = pref
      !
      do n = 1,npts
        vartmp(n) = gp(n,npres) * rhoa2 - freestream_ref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("P0") ! Total Pressure
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      ! NOTE: vartmp has been pre-loaded with the Mach number, and the auxiliary
      !       variable in gp has already been switched to gamma
      !
      do n = 1,npts
        mach = vartmp(n)
        pres = gp(n,npres) * rhoa2
        gamm = gp(n,ngamm)
        gm1 = gamm - one
        gam_over_gm1 = gamm / gm1
        vartmp(n) = pres*( (one + half*gm1*mach*mach)**gam_over_gm1 )
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("Mach") ! Mach number
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      ! NOTE: vartmp has been pre-loaded with the Mach number
      !
      if (itestcase == Freestream_Preservation) then
        do n = 1,npts
          vartmp(n) = vartmp(n) - machref
        end do
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("Entropy") ! Entropy
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (ncs == 0) then
        do n = 1,npts
          gp(n,ngamm) = gam
        end do
      end if
      !
      if (itestcase == Inviscid_Gaussian_Bump) then
        do n = 1,npts
          s = (gp(n,npres)/pref_nd) / (cv(n,nec)**gp(n,ngamm))
          vartmp(n) = abs( one - s )
        end do
      else
        do n = 1,npts
          vartmp(n) = log(gp(n,npres)) - gp(n,ngamm)*log(cv(n,nec))
        end do
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("Temperature") ! Temperature
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (itestcase == Freestream_Preservation) freestream_ref = tref
      !
      do n = 1,npts
        vartmp(n) = gp(n,ntemp) * tref - freestream_ref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("T0") ! Total Temperature
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      ! NOTE: vartmp has been pre-loaded with the Mach number, and the auxiliary
      !       variable in gp has already been switched to gamma
      !
      do n = 1,npts
        mach = vartmp(n)
        temp = gp(n,ntemp) * tref
        gm1 = gp(n,ngamm) - one
        vartmp(n) = temp*(one + half*gm1*mach*mach)
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("Enthalpy") ! Enthalpy
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (ncs > 0) then
        do n = 1,npts
          vartmp(n) = gp(n,nenth) * aref*aref
        end do
      else
        gam_over_gm1 = gam / (gam-one)
        do n = 1,npts
          vartmp(n) = gam_over_gm1*gp(n,npres)/cv(n,nec) * aref*aref
        end do
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("H0") ! Total Enthalpy
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (ncs > 0) then
        do n = 1,npts
          vartmp(n) = gp(n,nenth) * aref*aref
        end do
      else
        gam_over_gm1 = gam / (gam-one)
        do n = 1,npts
          vartmp(n) = gam_over_gm1*gp(n,npres)/cv(n,nec) * aref*aref
        end do
      end if
      !
      if (nr == 3) then
        do n = 1,npts
          rho_inv = one/cv(n,nec)
          vx = cv(n,nmx)*rho_inv * aref
          vy = cv(n,nmy)*rho_inv * aref
          vz = cv(n,nmz)*rho_inv * aref
          q2 = vx*vx + vy*vy + vz*vz
          vartmp(n) = vartmp(n) + half*q2
        end do
      else if (nr == 2) then
        do n = 1,npts
          rho_inv = one/cv(n,nec)
          vx = cv(n,nmx)*rho_inv * aref
          vy = cv(n,nmy)*rho_inv * aref
          q2 = vx*vx + vy*vy
          vartmp(n) = vartmp(n) + half*q2
        end do
      else
        do n = 1,npts
          vx = cv(n,nmx)/cv(n,nec) * aref
          vartmp(n) = vartmp(n) + half*vx*vx
        end do
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("Rgas") ! Mixture gas constant
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (ncs > 0) then
        do m = 1,ncs
          do n = 1,npts
            cv(n,m) = cv(n,m)/cv(n,nec)
          end do
        end do
        call compute_mixture_rgas(1,npts,cv(1:npts,1:ncs),vartmp, &
                                  dimensional=true)
        do m = 1,ncs
          do n = 1,npts
            cv(n,m) = cv(n,m)*cv(n,nec)
          end do
        end do
      else
        do n = 1,npts
          vartmp(n) = rgasref
        end do
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("Gamma") ! Ratio of specific heats
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (ncs > 0) then
        do n = 1,npts
          vartmp(n) = gp(n,ngamm)
        end do
      else
        do n = 1,npts
          vartmp(n) = gam
        end do
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("Cp") ! Specific heat at constant pressure
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (ncs > 0) then
        do n = 1,npts
          vartmp(n) = gp(n,ncp) * cp_dim
        end do
      else
        do n = 1,npts
          vartmp(n) = cpref
        end do
      end if
      !
      call cgns_write_field_data
      !!}}}
    case ("Viscosity") ! Viscosity
      !!{{{
      if (governing_equations == NavierStokes_Eqns) then
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,npts
          vartmp(n) = gp(n,nvisc) * rhoa
        end do
        !
        call cgns_write_field_data
        !
      end if
      !!}}}
    case ("Conductivity") ! Conductivity
      !!{{{
      if (governing_equations == NavierStokes_Eqns) then
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,npts
          vartmp(n) = gp(n,ncond) * cond_dim
        end do
        !
        call cgns_write_field_data
        !
      end if
      !!}}}
    case ("Prandtl") ! Mixture Prandtl number
      !!{{{
      if (governing_equations == NavierStokes_Eqns) then
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        if (ncs > 0) then
          do n = 1,npts
            vartmp(n) = gp(n,nvisc) * gp(n,ncp) / gp(n,ncond)
          end do
        else
          do n = 1,npts
            vartmp(n) = Pr
          end do
        end if
        !
        call cgns_write_field_data
        !
      end if
      !!}}}
    case ("SpeciesEnthalpy") ! Species enthalpies
      !!{{{
      if (ncs > 0) then
        do m = 1,ncs
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          vartmp(1:npts) = spcvar(1:npts,m)
          !
          call cgns_write_field_data
          !
        end do
      end if
      !!}}}
    case ("SpeciesCp") ! Species specific heats at constant pressure
      !!{{{
      if (ncs > 0) then
        do m = 1,ncs
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          vartmp(1:npts) = spcvar(1:npts,m)
          !
          call cgns_write_field_data
          !
        end do
      end if
      !!}}}
    case ("cpi/R") ! Species Cp / Species gas constant
      !!{{{
      if (ncs > 0) then
        do m = 1,ncs
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          dim_val = one/(rspec(m)*cp_dim)
          vartmp(1:npts) = spcvar(1:npts,m)*dim_val
          !
          call cgns_write_field_data
          !
        end do
      end if
      !!}}}
    case ("hi/RT") ! Species enthalpies / (Species gas constant * temperature)
      !!{{{
      if (ncs > 0) then
        do m = 1,ncs
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          dim_val = rspec(m)*cp_dim*tref
          do n = 1,npts
            vartmp(n) = spcvar(n,m)/(gp(n,ntemp)*dim_val)
          end do
          !
          call cgns_write_field_data
          !
        end do
      end if
      !!}}}
    case ("hi/R") ! Species enthalpies / Species gas constant
      !!{{{
      if (ncs > 0) then
        do m = 1,ncs
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          dim_val = one/(rspec(m)*cp_dim)
          vartmp(1:npts) = spcvar(1:npts,m)*dim_val
          !
          call cgns_write_field_data
          !
        end do
      end if
      !!}}}
    case ("Vorticity") ! Vorticity
      !!{{{
      time_ref_inv = one/time_ref
      !
      vorticity_loop: &
      do l = 1,nr
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        i = (l+1) - ( (l+1)/4 )*3
        j = (l+2) - ( (l+2)/4 )*3
        !
        nmi = nmb-1+i
        nmj = nmb-1+j
        !
        do n = 1,npts
          !
          rho_inv = one/cv(n,nec)
          !
          mi = cv(n,nmi)
          mj = cv(n,nmj)
          !
          drhodi = dcvdx(n,i,nec) * time_ref_inv
          drhodj = dcvdx(n,j,nec) * time_ref_inv
          dmjdi = dcvdx(n,i,nmj) * time_ref_inv
          dmidj = dcvdx(n,j,nmi) * time_ref_inv
          !
          vartmp(n) = rho_inv*( (dmjdi-dmidj) - rho_inv*(mj*drhodi-mi*drhodj) )
          !
        end do
        !
        call cgns_write_field_data
        !
        if (nr < 3) exit vorticity_loop
        !
      end do vorticity_loop
      !!}}}
    case ("VorticityMagnitude") ! Vorticity magnitude
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,npts
        !
        w12 = dpvdx(n,2,nmx) - dpvdx(n,1,nmy)
        w13 = dpvdx(n,3,nmx) - dpvdx(n,1,nmz)
        w23 = dpvdx(n,3,nmy) - dpvdx(n,2,nmz)
        !
        vartmp(n) = sqrt( w12*w12 + w13*w13 + w23*w23 ) * aref
        !
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("StrainRateMagnitude") ! Strain rate magnitude
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,npts
        !
        dudx = dpvdx(n,1,nmx)
        dudy = dpvdx(n,2,nmx)
        dudz = dpvdx(n,3,nmx)
        !
        dvdx = dpvdx(n,1,nmy)
        dvdy = dpvdx(n,2,nmy)
        dvdz = dpvdx(n,3,nmy)
        !
        dwdx = dpvdx(n,1,nmz)
        dwdy = dpvdx(n,2,nmz)
        dwdz = dpvdx(n,3,nmz)
        !
        skk = one3 * (dudx + dvdy + dwdz)
        !
        s11 = dudx - skk
        s22 = dvdy - skk
        s33 = dwdz - skk
        !
        s12 = half * (dudy + dvdx)
        s13 = half * (dudz + dwdx)
        s23 = half * (dvdz + dwdy)
        !
        sijsij = s11**2 + s22**2 + s33**2 + two*(s12**2 + s13**2 + s23**2)
        !
        vartmp(n) = sqrt( two * sijsij ) * aref
        !
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("QCriterion") ! Q-criterion
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,npts
        !
        dudx = dpvdx(n,1,nmx)
        dudy = dpvdx(n,2,nmx)
        dudz = dpvdx(n,3,nmx)
        !
        dvdx = dpvdx(n,1,nmy)
        dvdy = dpvdx(n,2,nmy)
        dvdz = dpvdx(n,3,nmy)
        !
        dwdx = dpvdx(n,1,nmz)
        dwdy = dpvdx(n,2,nmz)
        dwdz = dpvdx(n,3,nmz)
        !
        w12 = half * (dudy - dvdx)
        w13 = half * (dudz - dwdx)
        w23 = half * (dvdz - dwdy)
        !
        wijwij = two * (w12**2 + w13**2 + w23**2)
        !
        skk = one3 * (dudx + dvdy + dwdz)
        !
        s11 = dudx - skk
        s22 = dvdy - skk
        s33 = dwdz - skk
        !
        s12 = half * (dudy + dvdx)
        s13 = half * (dudz + dwdx)
        s23 = half * (dvdz + dwdy)
        !
        sijsij = s11**2 + s22**2 + s33**2 + two*(s12**2 + s13**2 + s23**2)
        !
        vartmp(n) = half*(wijwij - sijsij)
       !vartmp(n) = vartmp(n) * aref*aref
        !
      end do
      !
      ! Save the original Q-criterion
      !
      vartmp_at_sp = vartmp(:) ! F2003 AUTO-REALLOC
      !
      ! Normalize the whole Q-criterion to the interval 0->1
      !
      qmin = minval(vartmp(:))
      qmax = maxval(vartmp(:))
      !
      vartmp(:) = (vartmp(:) - qmin)/(qmax-qmin)
      !
      call cgns_write_field_data
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      ! Copy the original positive values of Q-Criterion back to vartmp
      ! and make the negative values zero
      !
      vartmp(:) = merge( vartmp_at_sp(:) , zero , vartmp_at_sp(:) >= zero )
      !
      call cgns_write_field_data
      !
      !!}}}
    case ("DivergenceOfVelocity") ! Divergence of velocity
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,npts
        vartmp(n) = (dpvdx(n,1,nmx) + dpvdx(n,2,nmy) + dpvdx(n,3,nmz)) * aref
      end do
      !
      call cgns_write_field_data
      !!}}}
    case ("ResolutionMarker")
      !!{{{
      if (.not. allocated(vartmp_at_sp)) then
        allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
                         ierr,error_message)
      end if
      !
      ! Compute the resolution marker if it has yet to be computed
      !
      if (Limiter_Option >= 0 .or. &
          .not. allocated(marker) .or. &
          .not. allocated(scaled_marker)) then
        call resolution_indicator(apply_projection=fals)
      end if
      !
      ! Actual resolution marker value
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do nc = 1,ncell
        n0 = cell(nc)%beg_sp - 1
        np = cell(nc)%end_sp - n0
        do n = 1,np
          vartmp_at_sp(n0+n) = marker(nc)
        end do
      end do
      !
      call interpolate_for_output( vartmp_at_sp , vartmp )
      !
      call cgns_write_field_data
      !
      ! Scaled resolution marker value
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      if (Limiter_Option <= -1) then
        do nc = 1,ncell
          n0 = cell(nc)%beg_sp - 1
          np = cell(nc)%end_sp - n0
          this_order = cell(nc)%order
          limit_order = int( real(this_order,kind=wp) * &
                             (one-scaled_marker(nc)) )
          limit_order = max(1,min(this_order,limit_order))
          do n = 1,np
            vartmp_at_sp(n0+n) = limit_order
          end do
        end do
      else
        do nc = 1,ncell
          n0 = cell(nc)%beg_sp - 1
          np = cell(nc)%end_sp - n0
          do n = 1,np
            vartmp_at_sp(n0+n) = scaled_marker(nc)
          end do
        end do
      end if
      !
      call interpolate_for_output( vartmp_at_sp , vartmp )
      !
      call cgns_write_field_data
      !!}}}
    case ("JumpMarker")
      !!{{{
      allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      ! Compute the pressure jump marker if it has yet to be computed
      !
      if (Limiter_Option <= 1 .or. &
          .not. allocated(marker) .or. &
          .not. allocated(scaled_marker)) then
        call jump_indicator
      end if
      !
      ! Actual pressure jump marker value
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do nc = 1,ncell
        n0 = cell(nc)%beg_sp - 1
        np = cell(nc)%end_sp - n0
        do n = 1,np
          vartmp_at_sp(n0+n) = marker(nc)
        end do
      end do
      !
      call interpolate_for_output( vartmp_at_sp , vartmp )
      !
      call cgns_write_field_data
      !
      ! Scaled pressure jump marker value
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do nc = 1,ncell
        n0 = cell(nc)%beg_sp - 1
        np = cell(nc)%end_sp - n0
        do n = 1,np
          vartmp_at_sp(n0+n) = scaled_marker(nc)
        end do
      end do
      !
      call interpolate_for_output( vartmp_at_sp , vartmp )
      !
      call cgns_write_field_data
      !!}}}
    case ("Partition")
      !!{{{
      allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,n_solpts
        vartmp_at_sp(n) = mypnum+1
      end do
      !
      call interpolate_for_output( vartmp_at_sp , vartmp )
      !
      call cgns_write_field_data
      !!}}}
    case ("Timestep")
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      call interpolate_for_output( dtsp(1:n_solpts) , vartmp , nondim=time_ref )
      !
      call cgns_write_field_data
      !!}}}
    case ("CoMetrics")
     !!!{{{
     !allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
     !           stat=ierr , errmsg=error_message )
     !call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
     !                 ierr,error_message)
     !!
     !do l_phys = 1,nr
     !  do l_comp = 1,nr
     !    !
     !    nv = nv + 1
     !    ivar_s = isvar_s(nv)
     !    !
     !    do n = 1,n_solpts
     !      vartmp_at_sp(n) = metrics_array(n,l_phys,l_comp)
     !    end do
     !    !
     !    call interpolate_for_output( vartmp_at_sp , vartmp )
     !    !
     !    call cgns_write_field_data
     !    !
     !  end do
     !end do
     !!!}}}
    case ("ContraMetrics")
      !!{{{
      allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      do l_comp = 1,nr
        do l_phys = 1,nr
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          do n = 1,n_solpts
            vartmp_at_sp(n) = jac(n)*metrics_array(n,l_comp,l_phys)
          end do
          !
          call interpolate_for_output( vartmp_at_sp , vartmp )
          !
          call cgns_write_field_data
          !
        end do
      end do
      !!}}}
    case ("TimeMetrics")
      !!{{{
      allocate ( vartmp_at_sp(1:n_solpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"vartmp_at_sp",1,__LINE__,__FILE__, &
                       ierr,error_message)
      !
      if (allocated(time_metrics)) then
        !
        do l = 1,nr
          !
          nv = nv + 1
          ivar_s = isvar_s(nv)
          !
          do n = 1,n_solpts
            vartmp_at_sp(n) = one/time_metrics(l,n)
          end do
          !
          call interpolate_for_output( vartmp_at_sp , vartmp )
          !
          call cgns_write_field_data
          !
        end do
        !
        nv = nv + 1
        ivar_s = isvar_s(nv)
        !
        do n = 1,n_solpts
          vartmp_at_sp(n) = time_metrics(nr+1,n) * time_ref
        end do
        !
        call interpolate_for_output( vartmp_at_sp , vartmp )
        !
        call cgns_write_field_data
        !
      end if
      !!}}}
    case default
      !!{{{
      nv = nv + 1
      ivar_s = isvar_s(nv)
      !
      do n = 1,npts
        vartmp(n) = cv(n,nec) * rhoref
      end do
      !
      call cgns_write_field_data
      !!}}}
  end select
  !
  if (allocated(vartmp_at_sp)) then
    deallocate ( vartmp_at_sp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vartmp_at_sp",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  if (allocated(dpvdx)) then
    deallocate ( dpvdx , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"dpvdx",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(spcvar)) then
    deallocate ( spcvar , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"spcvar",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(other_variables)) then
    deallocate ( other_variables , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"other_variables",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  !-----------------------------------------------------------------------------
  !
  contains
  !
#include "Functions/cgns_write_field_data.f90"
  !
  !-----------------------------------------------------------------------------
  !
end subroutine output_variable
!
!###############################################################################
!
subroutine get_cgns_sol_variables(variable_list)
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : itestcase,governing_equations
  use ovar,   only : output_vorticity
  use ovar,   only : output_qcriterion
  use ovar,   only : output_conserved_variables
  use ovar,   only : output_temperature
  use ovar,   only : output_thermo_variables
  use ovar,   only : output_gradients
  use ovar,   only : output_transport_variables
  use ovar,   only : output_stagnation_temperature
  use ovar,   only : output_stagnation_pressure
  use ovar,   only : output_stagnation_variables
  use ovar,   only : output_marker_variables
  use ovar,   only : Resolution_Marker_Variable
  use ovar,   only : output_species_mass_fractions
  use ovar,   only : output_species_mole_fractions
  use ovar,   only : output_extra_species_variables
  use ovar,   only : output_debug_species_variables
  use ovar,   only : output_grdnt_species_variables
  !
  !.. Formal Arguments ..
#ifdef SPECIAL_FOR_INTEL
  character(len=:), allocatable, intent(out) :: variable_list(:)
#else
  character(len=CGLEN), allocatable, intent(out) :: variable_list(:)
#endif
  !
  !.. Local Scalars ..
  integer :: i,j,n,m
  !
  character(len=CGLEN) :: varname
  character(len=CGLEN) :: varcheck
  !
  !.. Local Allocatable Arrays ..
#ifdef SPECIAL_FOR_INTEL
  character(len=:), allocatable :: temp_variables(:)
#else
  character(len=CGLEN), allocatable :: temp_variables(:)
#endif
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_cgns_sol_variables"
  !
  character(len=*), parameter :: default_list(1:4) = ["Density ", &
                                                      "Velocity", &
                                                      "Pressure", &
                                                      "Mach    "]
  !
  ! AVAILABLE QUANTITIES TO OUTPUT AND ORDER IN WHICH THEY WILL BE OUTPUT
  !
  character(len=*), parameter :: var_order(1:49) = ["Density             ", &
                                                    "Momentum            ", &
                                                    "Energy              ", &
                                                    "Velocity            ", &
                                                    "Pressure            ", &
                                                    "Temperature         ", &
                                                    "Mach                ", &
                                                    "MassFractions       ", &
                                                    "MoleFractions       ", &
                                                    "SpeciesDensities    ", &
                                                    "Enthalpy            ", &
                                                    "Cp                  ", &
                                                    "Entropy             ", &
                                                    "P0                  ", &
                                                    "T0                  ", &
                                                    "H0                  ", &
                                                    "Viscosity           ", &
                                                    "Conductivity        ", &
                                                    "Prandtl             ", &
                                                    "Gamma               ", &
                                                    "Rgas                ", &
                                                    "SpeciesEnthalpy     ", &
                                                    "SpeciesCp           ", &
                                                    "cpi/R               ", &
                                                    "hi/RT               ", &
                                                    "hi/R                ", &
                                                    "GradSpeciesDensities", &
                                                    "GradDensity         ", &
                                                    "GradMomentum        ", &
                                                    "GradEnergy          ", &
                                                    "GradMassFractions   ", &
                                                    "GradVelocity        ", &
                                                    "Vorticity           ", &
                                                    "GradTemperature     ", &
                                                    "ResolutionMarker    ", &
                                                    "JumpMarker          ", &
                                                    "Timestep            ", &
                                                    "Partition           ", &
                                                    "ContraMetrics       ", &
                                                    "TimeMetrics         ", &
                                                    "dTdx_Breakdown      ", &
                                                    "DmatTemperature     ", &
                                                    "CoMetrics           ", &
                                                    "GradTempMagnitude   ", &
                                                    "GradDensityMagnitude", &
                                                    "VorticityMagnitude  ", &
                                                    "StrainRateMagnitude ", &
                                                    "QCriterion          ", &
                                                    "DivergenceOfVelocity"]
  !
  ! Conserved Variables
  character(len=*), parameter :: consv_var(1:3) = ["Density ", &
                                                   "Momentum", &
                                                   "Energy  "]
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate variable_list to the default_list
  !
  call add_sol_variables(variable_list,default_list)
  !
  ! Make sure vorticity is one of the variables that will
  ! be output if running the Taylor-Green vortex problem.
  !
  if (itestcase == Taylor_Green_Vortex) then
    call add_sol_variable(variable_list,"Vorticity")
  end if
  !
  ! Output vorticity and other variables useful for identifying
  ! vortical structures
  !
  if (output_vorticity) then
    call add_sol_variable(variable_list,"Vorticity")
    call add_sol_variable(variable_list,"DivergenceOfVelocity")
    if (nr > 2) then
      call add_sol_variable(variable_list,"VorticityMagnitude")
      call add_sol_variable(variable_list,"StrainRateMagnitude")
      call add_sol_variable(variable_list,"QCriterion")
    end if
  end if
  !
  ! Output Q-criterion without all the other vorticity variables
  !
  if (output_qcriterion) then
    if (nr > 2) then
      call add_sol_variable(variable_list,"QCriterion")
    end if
  end if
  !
  ! Make sure entropy is one of the variables that will be
  ! output if running any of the isentropic vortex problems.
  !
  if (any(itestcase == Transport_Problems)) then
    call add_sol_variable(variable_list,"Entropy")
  end if
  !
  if (output_conserved_variables) then
    call add_sol_variables(variable_list,consv_var)
    if (ncs > 1) then
      call add_sol_variable(variable_list,"SpeciesDensities")
    end if
  end if
  !
  if (output_temperature) then
    call add_sol_variable(variable_list,"Temperature")
  end if
  !
  if (output_thermo_variables) then
    call add_sol_variable(variable_list,"Pressure")
    call add_sol_variable(variable_list,"Temperature")
    call add_sol_variable(variable_list,"Enthalpy")
    call add_sol_variable(variable_list,"Cp")
    call add_sol_variable(variable_list,"Entropy")
  end if
  !
  if (output_stagnation_temperature) then
    call add_sol_variable(variable_list,"T0")
  end if
  !
  if (output_stagnation_pressure) then
    call add_sol_variable(variable_list,"P0")
  end if
  !
  if (output_stagnation_variables) then
    call add_sol_variable(variable_list,"P0")
    call add_sol_variable(variable_list,"T0")
    call add_sol_variable(variable_list,"H0")
  end if
  !
  ! Output the marker quantities if requested
  !
  if (output_marker_variables) then
    call add_sol_variable(variable_list,"ResolutionMarker")
    if (Resolution_Marker_Variable > 0) then
      if (Resolution_Marker_Variable == nec .or. ncs > 1) then
        m = 1
      else if (Resolution_Marker_Variable == nee) then
        m = 3
      else
        m = 2
      end if
      call add_sol_variable(variable_list,trim(consv_var(m)))
    else
      call add_sol_variables(variable_list,consv_var)
    end if
  end if
  !
  ! Add species mass fractions to the variable list if using a gas mixture.
  ! Also output any other requested variables for gas mixtures.
  !
  if (ncs > 0) then
    !
    if (output_species_mass_fractions) then
      call add_sol_variable(variable_list,"MassFractions")
    end if
    !
    if (output_species_mole_fractions) then
      call add_sol_variable(variable_list,"MoleFractions")
    end if
    !
    if (output_extra_species_variables) then
      call add_sol_variable(variable_list,"SpeciesDensities")
      call add_sol_variable(variable_list,"Enthalpy")
      call add_sol_variable(variable_list,"H0")
      call add_sol_variable(variable_list,"Gamma")
      call add_sol_variable(variable_list,"Cp")
      call add_sol_variable(variable_list,"Rgas")
      call add_sol_variable(variable_list,"Momentum")
      call add_sol_variable(variable_list,"Energy")
      call add_sol_variable(variable_list,"P0")
      call add_sol_variable(variable_list,"T0")
    end if
    !
    if (output_grdnt_species_variables) then
      call add_sol_variable(variable_list,"GradMassFractions")
      call add_sol_variable(variable_list,"GradSpeciesDensities")
      call add_sol_variable(variable_list,"GradDensity")
      call add_sol_variable(variable_list,"GradMomentum")
      call add_sol_variable(variable_list,"GradEnergy")
      call add_sol_variable(variable_list,"GradVelocity")
      call add_sol_variable(variable_list,"GradTemperature")
    end if
    !
    if (output_debug_species_variables) then
      call add_sol_variable(variable_list,"SpeciesEnthalpy")
      call add_sol_variable(variable_list,"SpeciesCp")
      call add_sol_variable(variable_list,"hi/RT")
      call add_sol_variable(variable_list,"hi/R")
      call add_sol_variable(variable_list,"cpi/R")
      call add_sol_variable(variable_list,"dTdx_Breakdown")
    end if
    !
  end if
  !
  ! Output transport properties if solving the Navier-Stokes equations
  !
  if (output_transport_variables) then
    if (governing_equations == NavierStokes_Eqns) then
      !
      call add_sol_variable(variable_list,"Viscosity")
      call add_sol_variable(variable_list,"Conductivity")
      !
      if (ncs > 0) then
        if (output_extra_species_variables) then
          call add_sol_variable(variable_list,"Prandtl")
        end if
      end if
      !
    end if
  end if
  !
  ! If outputing gradients was requested, add the gradients of whatever
  ! variables are in the variable list.
  !
  if (output_gradients) then
    !
    allocate ( temp_variables , source=variable_list )
    !
    do i = 1,size(temp_variables)
      select case (uppercase(trim(adjustl(temp_variables(i)))))
        case ("DENSITY")
          call add_sol_variable(variable_list,"GradDensity")
        case ("VELOCITY")
          call add_sol_variable(variable_list,"GradVelocity")
        case ("ENERGY")
          call add_sol_variable(variable_list,"GradEnergy")
        case ("TEMPERATURE")
          call add_sol_variable(variable_list,"GradTemperature")
        case ("MOMENTUM")
          call add_sol_variable(variable_list,"GradMomentum")
        case ("MASSFRACTIONS")
          call add_sol_variable(variable_list,"GradMassFractions")
        case ("SPECIESDENSITIES")
          call add_sol_variable(variable_list,"GradSpeciesDensities")
      end select
    end do
    !
    deallocate ( temp_variables )
    !
  end if
  !
  ! Make sure there are no duplicate variables being output
  !
  allocate ( temp_variables , source=variable_list )
  !
  n = 0
  copy_loop: &
  do i = 1,size(variable_list)
    varcheck = uppercase( trim(adjustl(variable_list(i))) )
    do j = 1,n
      varname = uppercase( trim(adjustl(temp_variables(j))) )
      if (varname == varcheck) then
        ! This variable already exists in temp_variables, so dont add this
        ! to temp_variables and go on to the next variable in variable_list.
        cycle copy_loop
      end if
    end do
    n = n + 1
    temp_variables(n) = trim(adjustl(variable_list(i)))
  end do copy_loop
  !
  ! Use F2003 AUTO-REALLOCATION to cutoff the duplicate variables
  ! still contained in temp_variables(n+1:size(temp_variables)).
  !
  temp_variables = temp_variables(1:n)
  !
  ! Use move_alloc to over-write variable_list with temp_variables
  !
  call move_alloc( from=temp_variables , to=variable_list )
  !
  ! Put the output variables in order
  !
  allocate ( temp_variables , source=variable_list )
  !
  n = 0
  do i = 1,size(var_order)
    varname = uppercase( trim(adjustl(var_order(i))) )
    reorder_loop: &
    do j = 1,size(variable_list)
      varcheck = uppercase( trim(adjustl(variable_list(j))) )
      if (varname == varcheck) then
        n = n + 1
        temp_variables(n) = trim(adjustl(variable_list(j)))
        exit reorder_loop
      end if
    end do reorder_loop
  end do
  !
  ! Use move_alloc to over-write variable_list with temp_variables
  !
  call move_alloc( from=temp_variables , to=variable_list )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_cgns_sol_variables
!
!###############################################################################
!
subroutine add_sol_variable(variable_list,new_variable)
  !
  !.. Formal Arguments ..
  character(len=*),                     intent(in) :: new_variable
#ifdef SPECIAL_FOR_INTEL
  character(len=:),     allocatable, intent(inout) :: variable_list(:)
#else
  character(len=CGLEN), allocatable, intent(inout) :: variable_list(:)
#endif
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "add_sol_variable"
  !
continue
  !
  if (.not. allocated(variable_list)) then
    !
#ifdef SPECIAL_FOR_INTEL
    n = len_trim(adjustl(new_variable))
    ! USING F2003 AUTO-REALLOCATION
    variable_list = [ character(len=n) :: trim(adjustl(new_variable)) ]
#else
    ! USING F2003 AUTO-REALLOCATION
    variable_list = [ adjustl(new_variable) ]
#endif
    !
  else
    !
#ifdef SPECIAL_FOR_INTEL
    n = max( len_trim(adjustl(new_variable)) , len(variable_list) )
    ! USING F2003 AUTO-REALLOCATION
    variable_list = [ character(len=n) :: variable_list , &
                                          trim(adjustl(new_variable)) ]
#else
    ! USING F2003 AUTO-REALLOCATION
    variable_list = [ variable_list , adjustl(new_variable) ]
#endif
    !
  end if
  !
end subroutine add_sol_variable
!
!###############################################################################
!
subroutine add_sol_variables(variable_list,new_variables)
  !
  !.. Formal Arguments ..
  character(len=*),                     intent(in) :: new_variables(:)
#ifdef SPECIAL_FOR_INTEL
  character(len=:),     allocatable, intent(inout) :: variable_list(:)
#else
  character(len=CGLEN), allocatable, intent(inout) :: variable_list(:)
#endif
  !
  !.. Local Scalars ..
  integer :: i,ibeg,n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "add_sol_variables"
  !
continue
  !
  if (.not. allocated(variable_list)) then
    !
#ifdef SPECIAL_FOR_INTEL
    n = len_trim(adjustl(new_variables(1)))
    ! USING F2003 AUTO-REALLOCATION
    variable_list = [ character(len=n) :: trim(adjustl(new_variables(1))) ]
#else
    ! USING F2003 AUTO-REALLOCATION
    variable_list = [ adjustl(new_variables(1)) ]
#endif
    !
    ibeg = 2
    !
  else
    !
    ibeg = 1
    !
  end if
  !
  do i = ibeg,size(new_variables)
#ifdef SPECIAL_FOR_INTEL
    n = max( len_trim(adjustl(new_variables(i))) , len(variable_list) )
    ! USING F2003 AUTO-REALLOCATION
    variable_list = [ character(len=n) :: variable_list , &
                                          trim(adjustl(new_variables(i))) ]
#else
    ! USING F2003 AUTO-REALLOCATION
    variable_list = [ variable_list , adjustl(new_variables(i)) ]
#endif
  end do
  !
end subroutine add_sol_variables
!
!###############################################################################
!
subroutine write_parallel_cgns_time_ave_file()
  !
  !.. Use Statements ..
  use geovar,  only : nr,n_solpts
  use ovar,    only : aref,rhoref,tref,tkeref,tlsref
  use ovar,    only : profile_io_cgns,profile_io_cgns_subs
  use ovar,    only : max_wall_time => max_time_ave_cgns_wall_time
  use flowvar, only : uavesp,time_ave_variables
  !
  !.. Local Scalars ..
  integer :: l,nv,naq,ierr
  !
  real(wp) :: beg_wall_time
  real(wp) :: end_wall_time
  real(wp) :: cgns_wall_time
  !
  character(len=1)     :: ta_char
  character(len=1)     :: vmult
  character(len=CGLEN) :: varname
  character(len=CGLEN) :: sect_name
  character(len=170)   :: cgns_sol_file
  character(len=300)   :: link_path
  !
  integer(CBT) :: ifile_s,ibase_s,izone_s,isol_s,ivar_s
  !
  !.. Local Allocatable Arrays ..
  integer(CBT), allocatable :: isvar_s(:)
  !
  real(CRT), allocatable :: var(:)
  real(wp),  allocatable :: vartmp(:)
  real(wp),  allocatable :: var_dim(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_parallel_cgns_time_ave_file"
  character(len=*), parameter :: tname = "TimeAveragedSolution"
  character(len=*), parameter :: fname = "solution.time_ave.cgns"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
#ifndef SPECIAL_FOR_PGI
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      write (iout,1) pname
      flush (iout)
      beg_wall_time = mpi_wtime()
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  naq = size(uavesp,dim=2)
  !
  cgns_sol_file = trim(adjustl(cgns_dir)) // fname
  !
  ! Open the CGNS solution file if iflag is 1
  !
  call cgns_timer(start_cgns_file_timer,profile_io_cgns_subs,cgns_sol_file)
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cgp_open_f")
  call cgp_open_f(cgns_sol_file,CG_MODE_WRITE,ifile_s,cgierr)
  call io_error(pname,cgns_sol_file,1,__LINE__,__FILE__,cgierr)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cgp_open_f")
  !
  ! Set the CGNS base (who knows what this is for)
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_base_write_f")
  call cg_base_write_f(ifile_s,base_name,int(nr,kind=CBT), &
                       int(nr,kind=CBT),ibase_s,cgierr)
  call cgns_error(pname,cgns_sol_file,"cg_base_write_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_base_write_f")
  !
  ! Write the title for this solution
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_zone_write_f")
  call cg_zone_write_f(ifile_s,ibase_s,zone_name,zone_size, &
                       UNSTRUCTURED,izone_s,cgierr)
  call cgns_error(pname,cgns_sol_file,"cg_zone_write_f", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_zone_write_f")
  !
  ! Create a link to the grid coordinates in the grid file
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cg_goto_f")
  call cg_goto_f(ifile_s,ibase_s,cgierr, &
                 "Zone_t",izone_s, &
                 CGNS_GOTO_END)
  call cgns_error(pname,cgns_sol_file,"cg_goto_f",cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cg_goto_f")
  !
  link_path = base_name//"/"//zone_name//"/GridCoordinates"
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  call cg_link_write_f("GridCoordinates",trim(cgns_base_fname), &
                       trim(link_path),cgierr)
  call cgns_error(pname,cgns_sol_file,"cg_link_write_f(GridCoordinates)", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  !
  ! Create a link to the element connectivity in the grid file
  !
  if (nr == 2) then
    sect_name = "Quadrilaterals"
  else if (nr == 3) then
    sect_name = "Hexahedrals"
  end if
  link_path = base_name//"/"//zone_name//"/"//trim(sect_name)
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  call cg_link_write_f(trim(sect_name),trim(cgns_base_fname), &
                       trim(link_path),cgierr)
  call cgns_error(pname,cgns_sol_file, &
                  "cg_link_write_f("//trim(sect_name)//")", &
                  cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_link_write_f")
  !
  ! Allocate isvar_s and var_dim to the number of time averaged variables
  !
  allocate ( isvar_s(1:naq) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"isvar_s",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( var_dim(1:naq) , source=one , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"var_dim",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Write the flow variables
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_sol_write_f")
  call cg_sol_write_f(ifile_s,ibase_s,izone_s,tname,VERTEX,isol_s,cgierr)
  call cgns_error(pname,cgns_sol_file,"cg_sol_write_f",cgierr,__LINE__,__FILE__)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                  "cg_sol_write_f")
  !
  !
  ! ##############################################################
  ! #####   WRITING TIME AVERAGED SOLUTION FIELD VARIABLES   #####
  ! ##############################################################
  !
  ! Create data nodes for the time-averaged solution variables
  !
  do nv = 1,naq
    !
    vmult = ""
    varname = ""
    do l = 1,len_trim(time_ave_variables(nv))
      ta_char = uppercase( time_ave_variables(nv)(l:l) )
      select case (ta_char)
        case ("D")
          varname = trim(adjustl(varname)) // trim(vmult) // "Rho"
          var_dim(nv) = var_dim(nv) * rhoref
        case ("X")
          varname = trim(adjustl(varname)) // trim(vmult) // "Vx"
          var_dim(nv) = var_dim(nv) * aref
        case ("Y")
          varname = trim(adjustl(varname)) // trim(vmult) // "Vy"
          var_dim(nv) = var_dim(nv) * aref
        case ("Z")
          varname = trim(adjustl(varname)) // trim(vmult) // "Vz"
          var_dim(nv) = var_dim(nv) * aref
        case ("E")
          varname = trim(adjustl(varname)) // trim(vmult) // "E"
          var_dim(nv) = var_dim(nv) * rhoref * aref * aref
        case ("P")
          varname = trim(adjustl(varname)) // trim(vmult) // "P"
          var_dim(nv) = var_dim(nv) * rhoref * aref * aref
        case ("T")
          varname = trim(adjustl(varname)) // trim(vmult) // "T"
          var_dim(nv) = var_dim(nv) * tref
        case ("R")
          varname = trim(adjustl(varname)) // trim(vmult) // "Vr"
          var_dim(nv) = var_dim(nv) * aref
        case ("F")
          varname = trim(adjustl(varname)) // trim(vmult) // "Vf"
          var_dim(nv) = var_dim(nv) * aref
        case ("K")
          varname = trim(adjustl(varname)) // trim(vmult) // "TKE"
          var_dim(nv) = var_dim(nv) * tkeref
        case ("W")
          varname = trim(adjustl(varname)) // trim(vmult) // "TLS"
          var_dim(nv) = var_dim(nv) * tlsref
      end select
      vmult = "*"
    end do
    !
    call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_field_write_f")
    call cgp_field_write_f(ifile_s,ibase_s,izone_s,isol_s,CGNS_KIND_REAL, &
                           trim(adjustl(varname)),isvar_s(nv),cgierr)
    call cgns_error(pname,cgns_sol_file,"cgp_field_write_f", &
                    cgierr,__LINE__,__FILE__, &
                    "Write node: Time-averaged variable",nv)
    call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs, &
                    "cgp_field_write_f")
    !
  end do
  !
  ! Preallocate var for writing variables in the correct precision
  ! to the CGNS file
  !
  if (allocated(vartmp)) then
    deallocate ( vartmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vartmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(var)) then
    deallocate ( var , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"var",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( vartmp(1:n_my_output_pts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"vartmp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( var(1:n_my_output_pts) , source=real(0,CRT) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"var",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Now create the output data and write it to the CGNS file in parallel
  !
  do nv = 1,naq
    !
    call interpolate_for_output( uavesp(1:n_solpts,nv) , &
                                 vartmp , nondim=var_dim(nv) )
   !vartmp = uavesp(1:n_solpts,nv) * var_dim(nv) ! F2003 AUTO-REALLOCATION
    !
    ivar_s = isvar_s(nv)
    !
    call cgns_write_field_data
    !
  end do
  !
  ! ################################################
  ! #####   FINISHED WRITING FIELD VARIABLES   #####
  ! ################################################
  !
  call set_cgns_queue(STOP_CGNS_QUEUE,pname,cgns_sol_file)
  !
  ! Close the CGNS FILE
  !
  call cgns_timer(start_cgns_procedure_timer,profile_io_cgns_subs,"cgp_close_f")
  call cgp_close_f(ifile_s,cgierr)
  call io_error(pname,cgns_sol_file,2,__LINE__,__FILE__,cgierr)
  call cgns_timer(stop_cgns_procedure_timer,profile_io_cgns_subs,"cgp_close_f")
  !
  call cgns_timer(stop_cgns_file_timer,profile_io_cgns_subs,cgns_sol_file)
  !
  ! Deallocate the arrays var and vartmp
  !
  if (allocated(var)) then
    deallocate ( var , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"var",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(vartmp)) then
    deallocate ( vartmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vartmp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(isvar_s)) then
    deallocate ( isvar_s , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"isvar_s",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(var_dim)) then
    deallocate ( var_dim , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"var_dim",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      end_wall_time = mpi_wtime()
      cgns_wall_time = end_wall_time - beg_wall_time
      max_wall_time = max(max_wall_time,cgns_wall_time)
      write (iout,2) pname,trim(adjustl(cgns_sol_file)),cgns_wall_time
      flush (iout)
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," Entering ",a)
  2 format (1x,"Leaving ",a,/,:, &
            5x,"Wall time taken to write the file '",a, &
            "' = ",f12.6," seconds",/)
  !
  !-----------------------------------------------------------------------------
  !
  contains
  !
#include "Functions/cgns_write_field_data.f90"
  !
  !-----------------------------------------------------------------------------
  !
end subroutine write_parallel_cgns_time_ave_file
!
!###############################################################################
!
subroutine get_initial_data_for_cgns(my_cells_beg,my_cells_end,ipelem)
  !
  !.. Use Statements ..
  use geovar,    only : ncell,n_global_cell,cell
  use geovar,    only : pst_geom
  use geovar,    only : global_cell_geom
  use geovar,    only : global_cell_order
  use geovar,    only : cell_map
  use order_mod, only : o_order,geom_solpts
  !
  !.. Formal Arguments ..
  integer(CST),              intent(inout) :: my_cells_beg
  integer(CST),              intent(inout) :: my_cells_end
  integer(CST), allocatable, intent(inout) :: ipelem(:,:)
  !
  !.. Local Scalars ..
  integer :: n,nc,np,ierr,icpu
  integer :: gmin,gmax,this_subcells
  integer :: this_cell,this_geom
  integer :: this_order,output_order
  !
  integer(CST) :: p1,this_npts
  integer(CST) :: n_my_output_cells
  integer(CST) :: n_total_output_cells
  integer(CST) :: n_total_output_pts
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_initial_data_for_cgns"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize n_total_output_cells and n_total_output_pts which give:
  !
  ! n_total_output_cells: total number of output cells for the global grid
  ! n_total_output_pts: total number of output points for the global grid
  !
  n_total_output_pts = 0_CST
  n_total_output_cells = 0_CST
  !
  do icpu = 1,ncpu
    !
    if (icpu-1 == mypnum) then
      my_pts_beg = n_total_output_pts + 1_CST
      my_cells_beg = n_total_output_cells + 1_CST
    end if
    !
    do nc = 1,size(cell_map(icpu)%loc_to_glb)
      !
      ! Index for grid cell on the current processor
      !
      this_cell = cell_map(icpu)%loc_to_glb(nc)
      !
      ! Geometry of this grid cell
      !
      this_geom = global_cell_geom(this_cell)
      !
      ! Polynomial order for this grid cell
      !
      this_order = global_cell_order(this_cell)
      !
      ! Number of solution points for this cell geometry at this order
      !
      output_order = max(this_order,o_order)
      !
      this_npts = int(geom_solpts(this_geom,output_order),kind=kind(this_npts))
      !
      ! Get the number of subcells for this cell geometry at this order.
      ! The number of subcells is basically:
      !     this_order**(spatial dimension of this_geom) or
      ! NOTE: Make sure that the number of subcells is at least 1,
      !       particularly for the case of this_order=0
      !
      this_subcells = max( 1 , output_order**geom_dimen(this_geom) )
      !
      ! Add the number of solution points for this cell geometry at
      ! this order to the running count of the global total number of
      ! output solution points.
      !
      n_total_output_pts = n_total_output_pts + this_npts
      !
      ! Add the number of subcells for this cell geometry at this order
      ! to the running count of the global total number of output
      ! subcells.
      !
      n_total_output_cells = n_total_output_cells + &
                             int(this_subcells,kind=kind(n_total_output_cells))
      !
    end do
    !
    if (icpu-1 == mypnum) then
      my_pts_end = n_total_output_pts
      my_cells_end = n_total_output_cells
    end if
    !
  end do
  !
  ! Save n_total_output_pts and n_total_output_cells to zone_size
  !
  zone_size = [n_total_output_pts,n_total_output_cells,0_CST]
  !
  ! n_my_output_pts: number of output points for this processor
  !
  n_my_output_pts = my_pts_end - my_pts_beg + 1_CST
  !
  ! n_my_output_cells: number of output cell for this processor
  !
  n_my_output_cells = my_cells_end - my_cells_beg + 1_CST
  !
  ! Some quick stuff needed before creating the connectivity array
  !
  gmin = minval( global_cell_geom(1:n_global_cell) )
  gmax = maxval( global_cell_geom(1:n_global_cell) )
  !
  if (any(Geom_Unknown == [gmin,gmax])) then
    write (error_message,1)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  else if (geom_dimen(gmin) /= geom_dimen(gmax)) then
    write (error_message,2)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Get the number of nodes needed to define grid cells of this dimension
  ! NOTE: We will be creating degenerate grid cells for unstructured type
  !       cells such as Geom_Tria,Geom_Tetr,Geom_Pyra,Geom_Pris.
  ! Therefore, the number of nodes needed are:
  !         Geom_0D  => 2**0 = 1
  !         Geom_1D  => 2**1 = 2
  !         Geom_2D  => 2**2 = 4
  !         Geom_3D  => 2**3 = 8
  !
  np = 2**geom_dimen(gmin)
  !
  ! Create the connectivity array
  !
  if (allocated(ipelem)) then
    deallocate ( ipelem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ipelem",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( ipelem(1:np,1:n_my_output_cells) , source=0_CST , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ipelem",1,__LINE__,__FILE__,ierr,error_message)
  !
  nc = 0
  p1 = my_pts_beg
  !
  do this_cell = 1,ncell
    !
    ! Geometry of this grid cell
    !
    this_geom = cell(this_cell)%geom
    !
    ! Polynomial order for this grid cell
    !
    this_order = cell(this_cell)%order
    !
    ! Number of solution points for this cell geometry at this order
    !
    output_order = max(this_order,o_order)
    !
    this_npts = int(geom_solpts(this_geom,output_order),kind=kind(this_npts))
    !
    ! Get the number of subcells for this cell geometry at this order.
    ! The number of subcells is basically:
    !     this_order**(spatial dimension of this_geom) or
    ! NOTE: Make sure that the number of subcells is at least 1,
    !       particularly for the case of this_order=0
    !
    this_subcells = max( 1 , output_order**geom_dimen(this_geom) )
    !
    ! Loop through the number of subcells and add the subconnectivity
    ! to the first solution point of the current grid cell to get the
    ! connectivity for the current subcell.
    !
    do n = 1,this_subcells
      !
      nc = nc + 1
      !
      ipelem(:,nc) = p1 + int( pst_geom(this_geom,output_order)% &
                                              connectivity(:,n) , &
                               kind=kind(ipelem) )
      !
    end do
    !
    p1 = p1 + this_npts
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (" A grid cell of an unknown geometry type was found!")
  2 format (" The minimum and maximum dimensions of the grid cells do not", &
            " match. For example, this would happen if there is a 2D cell", &
            " within a 3D grid. This configuration is currently not supported.")
  !
end subroutine get_initial_data_for_cgns
!
!###############################################################################
!
subroutine TaylorGreen_8s_vorticity_solution()
  !
  !.. Use Statements ..
  use ovar, only : write_TaylorGreen_8s_solution
  use ovar, only : profile_io_cgns
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  integer(CBT) :: ifile_v,ibase_v,izone_v
  integer(CBT) :: isect_v,isol_v,isvar_v
  !
  integer(CST) :: my_beg,my_end,my_pts
  !
  real(wp) :: beg_wall_time
  real(wp) :: end_wall_time
  real(wp) :: cgns_wall_time
  !
  character(len=170) :: cgns_vort_file
  !
  integer(INT_MPI) :: color,key
  _MPI_COMM_TYPE_ :: vort_comm
  logical(lk) :: i_am_vort_cpu
  !
  !.. Local Arrays ..
  integer(CST), dimension(1:3) :: my_vort_pts
  integer(CST), dimension(1:3) :: my_vort_cells
  !
  integer(CBT), dimension(1:2) :: ixyz_v
  integer(CST), dimension(1:3) :: vort_zone_size
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:) :: my_vort_faces
  !
  integer(CST), allocatable, dimension(:,:) :: ipelem
  !
  real(CRT), allocatable, dimension(:) :: vary
  real(CRT), allocatable, dimension(:) :: varz
  !
  !.. Local Parameters ..
  character(len=*), parameter :: bname = "TGV_vort_8s_base"
  character(len=*), parameter :: zname = "TGV_vort_8s_zone"
  character(len=*), parameter :: fname = "TGV_vorticity_8s.cgns"
  character(len=*), parameter :: pname = "TaylorGreen_8s_vorticity_solution"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
#ifndef SPECIAL_FOR_PGI
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      write (iout,1) pname
      flush (iout)
      beg_wall_time = mpi_wtime()
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  cgns_vort_file = trim(adjustl(cgns_dir)) // fname
  !
  ! Get/collect all the data necessary to output the 8s vorticity solution
  !
  call get_data_for_8s_vorticity_solution(my_vort_pts,my_vort_cells, &
                                          vort_zone_size,my_vort_faces,ipelem)
  !
  my_pts = my_vort_pts(1)
  my_beg = my_vort_pts(2)
  my_end = my_vort_pts(3)
  !
  i_am_vort_cpu = (my_pts > 0)
  color = merge( 1_INT_MPI , MPI_UNDEFINED , i_am_vort_cpu )
  key = 0_INT_MPI
  !
  call mpi_comm_split(MPI_COMM_WORLD,color,key,vort_comm,mpierr)
  !
  if (i_am_vort_cpu) then
    !
#ifdef USE_MPI_F08
    call cgp_mpi_comm_f(vort_comm%mpi_val,cgierr)
#else
    call cgp_mpi_comm_f(vort_comm,cgierr)
#endif
    !
    ! Open the CGNS file
    !
    call cgp_open_f(cgns_vort_file,CG_MODE_WRITE,ifile_v,cgierr)
    call io_error(pname,cgns_vort_file,1,__LINE__,__FILE__,cgierr)
    !
    ! Set the CGNS base (who knows what this is for)
    ! NOTE: We are creating this base 2-dimensional with surface (2D) cells
    !       since this file will contain a 2D solution along a boundary face
    !
    call cg_base_write_f(ifile_v,bname,2_CBT,2_CBT,ibase_v,cgierr)
    call cgns_error(pname,cgns_vort_file,"cg_base_write_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    !
    ! Write the title for this solution
    !
    call cg_zone_write_f(ifile_v,ibase_v,zname,vort_zone_size, &
                         UNSTRUCTURED,izone_v,cgierr)
    call cgns_error(pname,cgns_vort_file,"cg_zone_write_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    !
    ! Create data nodes for the coordinates
    !
    call cgp_coord_write_f(ifile_v,ibase_v,izone_v,CGNS_KIND_REAL, &
                           "CoordinateY",ixyz_v(1),cgierr)
    call cgns_error(pname,cgns_vort_file,"cgp_coord_write_f", &
                    cgierr,__LINE__,__FILE__,"Write node: Coordinate = 1 (Y)", &
                    not_parallel=true)
    !
    call cgp_coord_write_f(ifile_v,ibase_v,izone_v,CGNS_KIND_REAL, &
                           "CoordinateZ",ixyz_v(2),cgierr)
    call cgns_error(pname,cgns_vort_file,"cgp_coord_write_f", &
                    cgierr,__LINE__,__FILE__,"Write node: Coordinate = 2 (Z)", &
                    not_parallel=true)
    !
    ! Now create and write the coordinate data in parallel
    !
    ! Stop/Disable the CGNS queue if it is not to be used
    !
    call disable_cgns_queue(pname,cgns_vort_file)
    !
    ! For each boundary face within the grid region of interest, interpolate the
    ! coordinates of the solution points for the host cell to the boundary face
    ! and this interpolate those face points to the output face points.
    !
    allocate ( vary(1:my_pts) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vary",1,__LINE__,__FILE__,ierr,error_message)
    vary = real( zero , kind=kind(vary) )
    !
    allocate ( varz(1:my_pts) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"varz",1,__LINE__,__FILE__,ierr,error_message)
    varz = real( zero , kind=kind(varz) )
    !
    ! Interpolate the coordinates to the output points
    !
    call get_vort_face_coordinates(my_vort_faces,vary,varz)
    !
    ! Add the data for the y-coordinate to the output queue
    ! NOTE: It is an error for a processor to try an write data if it doesnt
    !       have any data to write, so make sure that only those processors
    !       that have data to write call the CGNS write subroutine.
    !
    cgierr = 0
    call cgp_coord_write_data_f(ifile_v,ibase_v,izone_v,ixyz_v(1), &
                                my_beg,my_end,vary,cgierr)
    call cgns_error(pname,cgns_vort_file,"cgp_coord_write_data_f", &
                    cgierr,__LINE__,__FILE__,"Write data: Coordinate = 1 (Y)", &
                    not_parallel=true)
    !
    ! Add the data for the z-coordinate to the output queue
    ! NOTE: It is an error for a processor to try an write data if it doesnt
    !       have any data to write, so make sure that only those processors
    !       that have data to write call the CGNS write subroutine.
    !
    cgierr = 0
    call cgp_coord_write_data_f(ifile_v,ibase_v,izone_v,ixyz_v(2), &
                                my_beg,my_end,varz,cgierr)
    call cgns_error(pname,cgns_vort_file,"cgp_coord_write_data_f", &
                    cgierr,__LINE__,__FILE__,"Write data: Coordinate = 2 (Z)", &
                    not_parallel=true)
    !
    ! Finally, flush the queue containing the
    ! coordinate variables to the CGNS file
    !
    call flush_cgns_queue(pname,cgns_vort_file)
    !
    ! Create the data node for the grid connectivity
    !
    call cgp_section_write_f(ifile_v,ibase_v,izone_v,"Quadrilaterals", &
                             QUAD_4,1_CST,vort_zone_size(2), &
                             0_CBT,isect_v,cgierr)
    call cgns_error(pname,cgns_vort_file,"cgp_section_write_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    !
    ! Write the element connectivity in parallel
    !
    ! Reset the queue.
    !
    call set_cgns_queue(BEGIN_CGNS_QUEUE,pname,cgns_vort_file)
    !
    my_beg = my_vort_cells(2)
    my_end = my_vort_cells(3)
    !
    ! Add the element connectivity data to the output queue.
    ! NOTE: It is an error for a processor to try an write data if it doesnt
    !       have any data to write, so make sure that only those processors
    !       that have data to write call the CGNS write subroutine.
    !
    cgierr = 0
    call cgp_elements_write_data_f(ifile_v,ibase_v,izone_v,isect_v, &
                                   my_beg,my_end,ipelem,cgierr)
    call cgns_error(pname,cgns_vort_file,"cgp_elements_write_data_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    !
    ! Finally, flush the queue containing the
    ! element connectivity to the CGNS file.
    !
    call flush_cgns_queue(pname,cgns_vort_file)
    !
    ! Deallocate the arrays varz, and ipelem
    ! NOTE: Keep vary allocated and use it to store the vorticity variable
    !
    if (allocated(varz)) then
      deallocate ( varz , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"varz",2,__LINE__,__FILE__,ierr,error_message)
    end if
    if (allocated(ipelem)) then
      deallocate ( ipelem , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"ipelem",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    ! Evaluate and write the vorticity solution
    !
    my_beg = my_vort_pts(2)
    my_end = my_vort_pts(3)
    !
    ! Create the solution node
    !
    call cg_sol_write_f(ifile_v,ibase_v,izone_v,"VorticitySolution", &
                        VERTEX,isol_v,cgierr)
    call cgns_error(pname,cgns_vort_file,"cg_sol_write_f", &
                    cgierr,__LINE__,__FILE__,not_parallel=true)
    !
    ! Create a data node for the vorticity solution
    !
    call cgp_field_write_f(ifile_v,ibase_v,izone_v,isol_v,CGNS_KIND_REAL, &
                           "VorticityMagnitude",isvar_v,cgierr)
    call cgns_error(pname,cgns_vort_file,"cgp_field_write_f", &
                    cgierr,__LINE__,__FILE__, &
                    "Write node: Solution variable = 1 (VorticityMagnitude)", &
                    not_parallel=true)
    !
    ! Compute the vorticity magnitude at each face point
    !
    call get_vort_face_solution(my_vort_faces,vary)
    !
    ! Reset the queue.
    !
    call set_cgns_queue(BEGIN_CGNS_QUEUE,pname,cgns_vort_file)
    !
    ! Add the data for current variable to the output queue.
    ! NOTE: It is an error for a processor to try an write data if it doesnt
    !       have any data to write, so make sure that only those processors
    !       that have data to write call the CGNS write subroutine.
    !
    cgierr = 0
    call cgp_field_write_data_f(ifile_v,ibase_v,izone_v,isol_v,isvar_v, &
                                my_beg,my_end,vary,cgierr)
    call cgns_error(pname,cgns_vort_file,"cgp_field_write_data_f", &
                    cgierr,__LINE__,__FILE__, &
                    "Write data: Solution variable = 1 (VorticityMagnitude)", &
                    not_parallel=true)
    !
    ! Finally, flush the queue containing
    ! the current variable to the CGNS file.
    !
    call flush_cgns_queue(pname,cgns_vort_file)
    !
    ! Deallocate vary now that we are finished with it
    !
    if (allocated(vary)) then
      deallocate ( vary , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"vary",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    ! Close the CGNS file
    !
    call cgp_close_f(ifile_v,cgierr)
    call io_error(pname,cgns_vort_file,2,__LINE__,__FILE__,cgierr)
    !
  end if
  !
  call cgp_mpi_comm_f(mpi_comm_world_val,cgierr)
  !
  if (profile_io_cgns) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    if (mypnum == glb_root) then
      end_wall_time = mpi_wtime()
      cgns_wall_time = end_wall_time - beg_wall_time
      write (iout,2) pname,trim(adjustl(cgns_vort_file)),cgns_wall_time
      flush (iout)
    end if
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
#endif
  !
  write_TaylorGreen_8s_solution = fals
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," Entering ",a)
  2 format (1x,"Leaving ",a,/,:, &
            5x,"Wall time taken to write the file '",a, &
            "' = ",f12.6," seconds",/)
  !
end subroutine TaylorGreen_8s_vorticity_solution
!
!###############################################################################
!
subroutine get_data_for_8s_vorticity_solution(my_vort_pts, &
                                              my_vort_cells, &
                                              vort_zone_size, &
                                              my_vort_faces, &
                                              ipelem)
  !
  !.. Use Statements ..
  use ovar, only : Lref,VortexGrid_DOF
  use geovar, only : nfbnd,face
  use geovar, only : xyz_nodes,pst_geom
  use order_mod, only : geom_solpts,n_order
  !
  !.. Formal Arguments ..
  integer(CST),              intent(inout) :: my_vort_pts(1:3)
  integer(CST),              intent(inout) :: my_vort_cells(1:3)
  integer(CST),              intent(inout) :: vort_zone_size(1:3)
  integer,      allocatable, intent(inout) :: my_vort_faces(:)
  integer(CST), allocatable, intent(inout) :: ipelem(:,:)
  !
  !.. Local Scalars ..
  integer :: nc,nf,n,ierr,ip
  integer :: this_face,face_geom,face_order
  integer :: face_subcells
  !
  integer(INT_MPI) :: int_size,cvd_size
  _MPI_DATA_TYPE_ :: mpi_vorttyp
  !
  integer(CST) :: p1,face_npts
  integer(CST) :: n_my_vort_pts,n_my_vort_cells
  integer(CST) :: my_vort_pts_beg,my_vort_cells_beg
  integer(CST) :: my_vort_pts_end,my_vort_cells_end
  integer(CST) :: n_total_vort_cells,n_total_vort_pts
  !
  real(wp) :: LPI,half_LPI,dx_eps,cell_dx,x,y,z
  !
  !.. Local Allocatable Arrays ..
  logical(lk),      allocatable, dimension(:)   :: maskf
  integer(CST), allocatable, dimension(:,:) :: cpu_vort_data
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_data_for_8s_vorticity"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Compute a few initial constants relating to the grid length
  !
  LPI = Lref * PI
  half_LPI = half * LPI
  !
  ! Compute the average cell height and then use that to get a small
  ! percentage of the cell height to use a limit to identify boundary
  ! faces within the grid region of interest
  !
  cell_dx = two * LPI / anint( real(VortexGrid_DOF,kind=wp) / &
                               real(n_order+1,kind=wp) , kind=wp )
  dx_eps = cell_dx * eps3 / ten
  !
  ! Initialize the variables giving the total number of vorticity points
  ! and subcells as well as the variables giving the same information for
  ! the respective information local to the current processor
  !
  n_total_vort_pts = 0_CST
  n_total_vort_cells = 0_CST
  !
  n_my_vort_pts = 0_CST
  n_my_vort_cells = 0_CST
  !
  allocate ( maskf(1:nfbnd) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"maskf",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Find the boundary faces that lie on the plane x = -PI*Lref within the region
  ! y = [0,0.5*PI*Lref] and z = [0.5*PI*Lref,PI*Lref]
  !
  maskf_loop: do nf = 1,nfbnd
    !
    bface_node_loop: do n = 1,size(face(nf)%nodes)
      !
      face_geom = face(nf)%geom
      face_order = face(nf)%order ! order of boundary face
      !
      ! NOTE: We dont care here if any of the nodes stored in face(nf)%nodes
      !       for this face are duplicates. The if statements will just produce
      !       a false value for the second time.  We only care if the value
      !       ip that is returned is zero because this will cause the xyz_nodes
      !       array to go out-of-bounds so cycle the inner loop if this
      !       situation arrises.
      !
      ip = face(nf)%nodes(n)
      !
      if (ip == 0) cycle bface_node_loop
      !
      x = xyz_nodes(1,ip)
      y = xyz_nodes(2,ip)
      z = xyz_nodes(3,ip)
      !
      if (abs(x + LPI) <= eps6) then
        if (z >= half_LPI + dx_eps) then
          if ( (y >= dx_eps) .and. (y <= half_LPI - dx_eps) ) then
            !
            ! Mark this face for quick identification later
            !
            maskf(nf) = true
            !
            ! Number of solution points for this face geometry at this order
            !
            face_npts = int(geom_solpts(face_geom,face_order), &
                            kind=kind(face_npts))
            !
            ! Get the number of subcells for this face geometry at this order.
            ! The number of subcells is basically:
            !     face_order**(spatial dimension of face_geom) or
            ! NOTE: Make sure that the number of subcells is at least 1,
            !       particularly for the case of face_order=0
            !
            face_subcells = max( 1 , face_order**geom_dimen(face_geom) )
            !
            ! Add the number of solution points and subcells to the total
            ! number local to this processor
            !
            n_my_vort_pts   = n_my_vort_pts   + face_npts
            n_my_vort_cells = n_my_vort_cells + &
                              int(face_subcells,kind=kind(n_my_vort_cells))
            !
            cycle maskf_loop
            !
          end if
        end if
      end if
      !
    end do bface_node_loop
    !
  end do maskf_loop
  !
  ! Create the array that stores the index for all boundary faces on this
  ! processor that lie within the grid region of interest
  !
  if (allocated(my_vort_faces)) then
    deallocate ( my_vort_faces , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"my_vort_faces",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  allocate ( my_vort_faces(1:count(maskf)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"my_vort_faces",1,__LINE__,__FILE__,ierr,error_message)
  !
  my_vort_faces = pack( intseq(1,nfbnd) , mask=maskf )
  !
  ! Deallocate the maskf array since it is no longer needed
  !
  deallocate ( maskf , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"maskf",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! We need to collect this information for all processors
  !
  allocate ( cpu_vort_data(1:2,1:ncpu) , source=0_CST , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cpu_vort_data",1,__LINE__,__FILE__,ierr,error_message)
  !
  cpu_vort_data(1,mypnum+1) = n_my_vort_pts
  cpu_vort_data(2,mypnum+1) = n_my_vort_cells
  !
  if (ncpu > 1) then
    !
#ifdef USE_MPI_F08
    call mpi_sizeof(reshape(cpu_vort_data,[size(cpu_vort_data)]), &
                    int_size,mpierr)
    call mpi_type_match_size(MPI_TYPECLASS_INTEGER,int_size,mpi_vorttyp,mpierr)
#else
    if (CST == I8) then
      mpi_vorttyp = MPI_INTEGER8
    else
      mpi_vorttyp = MPI_INTEGER
    end if
#endif
    !
    cvd_size = size(cpu_vort_data,dim=1,kind=int_mpi)
    call mpi_allgather(MPI_IN_PLACE, cvd_size,mpi_vorttyp, &
                       cpu_vort_data,cvd_size,mpi_vorttyp, &
                       MPI_COMM_WORLD,mpierr)
    !
  end if
  !
  ! Sum up all the vorticity points and subcells from all processors
  !
  n_total_vort_pts   = sum( cpu_vort_data(1,:) )
  n_total_vort_cells = sum( cpu_vort_data(2,:) )
  !
  ! Save the total number of vorticity points and subcells to vort_zone_size
  !
  vort_zone_size = [n_total_vort_pts,n_total_vort_cells,0_CST]
  !
  ! Get the beginning index for the vorticity points and subcells local to this
  ! processor by summing up all the vorticity points and subcells from all
  ! processors with a smaller rank than the this processor.
  ! NOTE: The second index of cpu_vort_data goes from 1:ncpu whereas the value
  !       given by mypnum is from 0:ncpu-1. Therefore, summing from 1:mypnum
  !       gives the range of all processors with a smaller rank.
  !
  my_vort_pts_beg   = sum( cpu_vort_data(1,1:mypnum) ) + 1_CST
  my_vort_cells_beg = sum( cpu_vort_data(2,1:mypnum) ) + 1_CST
  !
  ! Get the ending index for the vorticity points and subcells local to this
  ! processor by adding the numbers found in fmask_loop to the beginning indices
  !
  my_vort_pts_end   = (my_vort_pts_beg-1)   + n_my_vort_pts
  my_vort_cells_end = (my_vort_cells_beg-1) + n_my_vort_cells
  !
  ! Save these values in the my_vort_pts and my_vort_cells arrays
  !
  my_vort_pts   = [n_my_vort_pts,  my_vort_pts_beg,  my_vort_pts_end]
  my_vort_cells = [n_my_vort_cells,my_vort_cells_beg,my_vort_cells_end]
  !
  ! Deallocate the cpu_vort_data array since it is no longer needed
  !
  deallocate ( cpu_vort_data , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cpu_vort_data",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the connectivity array
  ! NOTE: The Taylor-Green vortex case is 3D so we know the boundary faces
  !       will always be 2D and they will require 4 nodes to define each
  !       subcell within the CGNS file
  !
  if (allocated(ipelem)) then
    deallocate ( ipelem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ipelem",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( ipelem(1:4,1:n_my_vort_cells) , source=0_CST , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ipelem",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Initialize the counters nc and p1 for this processor
  !
  nc = 0
  p1 = my_vort_pts_beg
  !
  ! Loop through the boundary faces on this processor that are within the
  ! grid region of interest and get the subcell connectivity for each face
  !
  do nf = 1,size(my_vort_faces)
    !
    this_face = my_vort_faces(nf)
    !
    face_geom = face(this_face)%geom   ! geometry of boundary face
    face_order = face(this_face)%order ! order of boundary face
    !
    ! Number of solution points for this face geometry at this order
    !
    face_npts = int(geom_solpts(face_geom,face_order),kind=kind(face_npts))
    !
    ! Get the number of subcells for this face geometry at this order.
    ! The number of subcells is basically:
    !     face_order**(spatial dimension of face_geom) or
    ! NOTE: Make sure that the number of subcells is at least 1,
    !       particularly for the case of face_order=0
    !
    face_subcells = max( 1 , face_order**geom_dimen(face_geom) )
    !
    ! Loop through the number of subcells and add the subconnectivity
    ! to the first solution point of the current grid cell to get the
    ! connectivity for the current subcell.
    !
    do n = 1,face_subcells
      !
      nc = nc + 1
      !
      ipelem(:,nc) = p1 + int( pst_geom(face_geom,face_order)% &
                                            connectivity(:,n) , &
                               kind=kind(ipelem) )
      !
    end do
    !
    p1 = p1 + face_npts
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_data_for_8s_vorticity_solution
!
!###############################################################################
!
subroutine get_vort_face_coordinates(my_vort_faces,y,z)
  !
  !.. Use Statements ..
  use order_mod,         only : maxSP,maxFP
  use order_mod,         only : geom_solpts
  use geovar,            only : face,cell,xyz,fp
  use interpolation_mod, only : interp
  use interpolation_mod, only : outerp,outerp_mv_mult
  !
  !.. Formal Arguments ..
  integer,   dimension(:),    intent(in) :: my_vort_faces
  real(CRT), dimension(:), intent(inout) :: y
  real(CRT), dimension(:), intent(inout) :: z
  !
  !.. Local Scalars ..
  integer :: i,j,nf,np,n1,n2,nfp,f1,f2,k,kl
  integer :: this_face,face_geom,face_order
  integer :: host_cell,host_geom,host_order
  real(wp) :: a
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxFP) :: face_y
  real(wp), dimension(1:maxFP) :: face_z
  real(wp), dimension(1:maxSP,1:2) :: txyz
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_vort_face_coordinates"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  f2 = 0
  !
  do nf = 1,size(my_vort_faces)
    !
    this_face = my_vort_faces(nf) ! boundary face of current profile
    !
    face_geom = face(this_face)%geom   ! geometry of boundary face
    face_order = face(this_face)%order ! geometry of boundary face
    !
    host_cell = face(this_face)%left%cell ! host cell on current boundary face
    !
    host_geom = cell(host_cell)%geom   ! host cell geometry type
    host_order = cell(host_cell)%order ! host cell solution order
    !
    n1 = cell(host_cell)%beg_sp ! beginning index for sol-pts in host cell
    n2 = cell(host_cell)%end_sp ! ending    index for sol-pts in host cell
    np = n2-n1+1                ! number of sol-pts in host cell
    !
    nfp = geom_solpts(face_geom,face_order) ! number of sol-pts on boundary face
    f1 = f2 + 1                  ! beginning index for sol-pts on boundary face
    f2 = f2 + nfp                 ! ending    index for sol-pts on boundary face
    !
    ! Create transpose of the xyz(2:3,n1:n2) to reduce cache misses
    !
    txyz(1:np,1:2) = transpose( xyz(2:3,n1:n2) )
    !
    ! Compute the profile quantities at each data point for this face
    !
    do k = 1,nfp
      !
      ! Index of flx-pt local to the current host cell
      !
      kl = face(this_face)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                        rot(face(this_face)%left%rotation)% &
                                            pts(k)
      !
      ! First interpolate to the face flux points
      !
      face_y(k) = zero
      face_z(k) = zero
      do i = interp(host_geom,host_order)%toFace(face_order)%row_ptr(kl)+1, &
             interp(host_geom,host_order)%toFace(face_order)%row_ptr(kl+1)
        !
        j = interp(host_geom,host_order)%toFace(face_order)%row_idx(i)
        a = interp(host_geom,host_order)%toFace(face_order)%row_val(i)
        !
        face_y(k) = face_y(k) + a * txyz(j,1)
        face_z(k) = face_z(k) + a * txyz(j,2)
        !
      end do
      !
    end do
    !
    ! Finally, interpolate the face flux points
    ! to the output face solution points
    !
    if (allocated(outerp(face_geom,face_order)%output)) then
      !
      call outerp_mv_mult( outerp(face_geom,face_order)%output , face_y(1:nfp) )
      call outerp_mv_mult( outerp(face_geom,face_order)%output , face_z(1:nfp) )
      !
    end if
    !
    y(f1:f2) = real( face_y(1:nfp) , kind=kind(y) )
    z(f1:f2) = real( face_z(1:nfp) , kind=kind(z) )
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_vort_face_coordinates
!
!###############################################################################
!
subroutine get_vort_face_solution(my_vort_faces,vorticity)
  !
  !.. Use Statements ..
  use order_mod, only : maxSP,maxFP
  use order_mod, only : geom_solpts
  use geovar, only : nr,face,cell,fp
  use ovar, only : time_ref
  use interpolation_mod, only : interp
  use interpolation_mod, only : outerp,outerp_mv_mult
  use flowvar, only : usp,dusp
  !
  !.. Formal Arguments ..
  integer,   dimension(:),    intent(in) :: my_vort_faces
  real(CRT), dimension(:), intent(inout) :: vorticity
  !
  !.. Local Scalars ..
  integer :: i,j,k,kl,l,m,nf,n0,nfp,f1,f2
  integer :: this_face,face_geom,face_order
  integer :: host_cell,host_geom,host_order
  real(wp) :: a
  real(wp) :: rho_inv,vx,vy,vz
  real(wp) :: decdx,dmydx,dmzdx
  real(wp) :: decdy,dmxdy,dmzdy
  real(wp) :: decdz,dmxdz,dmydz
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: vort_vec
  real(wp), dimension(1:nq,1:maxFP) :: cvfp
  real(wp), dimension(1:maxFP,1:nq) :: cv_tmp
  real(wp), dimension(1:nr,1:nq,1:maxFP) :: dcvdx
  real(wp), dimension(1:maxFP,1:nr,1:nq) :: dcvdx_tmp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_vort_face_solution"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  f2 = 0
  !
  do nf = 1,size(my_vort_faces)
    !
    this_face = my_vort_faces(nf) ! boundary face of current profile
    !
    face_geom = face(this_face)%geom   ! geometry of boundary face
    face_order = face(this_face)%order ! geometry of boundary face
    !
    host_cell = face(this_face)%left%cell ! host cell on current boundary face
    !
    host_geom = cell(host_cell)%geom   ! host cell geometry type
    host_order = cell(host_cell)%order ! host cell solution order
    !
    n0 = cell(host_cell)%beg_sp-1 ! beginning index for sol-pts in host cell
    !
    nfp = geom_solpts(face_geom,face_order) ! number of sol-pts on boundary face
    f1 = f2 + 1                  ! beginning index for sol-pts on boundary face
    f2 = f2 + nfp                 ! ending    index for sol-pts on boundary face
    !
    ! Compute the profile quantities at each data point for this face
    !
    do k = 1,nfp
      !
      ! Index of flx-pt local to the current host cell
      !
      kl = face(this_face)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                        rot(face(this_face)%left%rotation)% &
                                            pts(k)
      !
      ! First interpolate to the face flux points
      !
      do m = 1,nq
        cv_tmp(k,m) = zero
        dcvdx_tmp(k,1:nr,m) = zero
        do i = interp(host_geom,host_order)%toFace(face_order)%row_ptr(kl)+1, &
               interp(host_geom,host_order)%toFace(face_order)%row_ptr(kl+1)
          !
          j = interp(host_geom,host_order)%toFace(face_order)%row_idx(i)
          a = interp(host_geom,host_order)%toFace(face_order)%row_val(i)
          !
          cv_tmp(k,m) = cv_tmp(k,m) + a * usp(n0+j,m)
          do l = 1,nr
            dcvdx_tmp(k,l,m) = dcvdx_tmp(k,l,m) + a * dusp(n0+j,l,m)
          end do
          !
        end do
      end do
      !
    end do
    !
    ! Now, interpolate the vorticity at the face flux points
    ! to the output face solution points
    !
    if (allocated(outerp(face_geom,face_order)%output)) then
      !
      do m = 1,nq
        call outerp_mv_mult( outerp(face_geom,face_order)%output , &
                             cv_tmp(1:nfp,m) , cvfp(m,1:nfp) )
        do l = 1,nr
          call outerp_mv_mult( outerp(face_geom,face_order)%output , &
                               dcvdx_tmp(1:nfp,l,m) , dcvdx(l,m,1:nfp) )
        end do
      end do
      !
    else
      !
      cvfp(1:nq,1:nfp) = transpose( cv_tmp(1:nfp,1:nq) )
      !
      do k = 1,nfp
        do m = 1,nq
          do l = 1,nr
            dcvdx(l,m,k) = dcvdx_tmp(k,l,m)
          end do
        end do
      end do
      !
    end if
    !
    ! Finally, compute the magnitude of the vorticity vector at each face point
    !
    do k = 1,nfp
      !
      ! Compute the vorticity vector at this face point using the
      ! gradients of the momentums and the product rule
      !
      rho_inv = one/cvfp(nec,k)
      vx = cvfp(nmx,k)*rho_inv
      vy = cvfp(nmy,k)*rho_inv
      vz = cvfp(nmz,k)*rho_inv
      !
      decdx = dcvdx(1,nec,k)
      decdy = dcvdx(2,nec,k)
      decdz = dcvdx(3,nec,k)
      !
      dmxdy = dcvdx(2,nmx,k)
      dmxdz = dcvdx(3,nmx,k)
      !
      dmydx = dcvdx(1,nmy,k)
      dmydz = dcvdx(3,nmy,k)
      !
      dmzdx = dcvdx(1,nmz,k)
      dmzdy = dcvdx(2,nmz,k)
      !
      vort_vec(1) = rho_inv*( (dmzdy-dmydz) - (vz*decdy-vy*decdz) )
      vort_vec(2) = rho_inv*( (dmxdz-dmzdx) - (vx*decdz-vz*decdx) )
      vort_vec(3) = rho_inv*( (dmydx-dmxdy) - (vy*decdx-vx*decdy) )
      !
      ! Compute the magnitude of the vorticity vector
      ! NOTE: Dont forget to adjust the nondimensionalization of the
      !       vorticity by dividing through by time_ref to get vorticity*
      !
      vorticity(f1-1+k) = real( norm2(vort_vec) / time_ref , &
                                kind=kind(vorticity) )
      !
    end do
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_vort_face_solution
!
!###############################################################################
!
pure function grad_of_selected_pv(ndr,npv,uu,dcvdx) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  integer,                    intent(in) :: ndr ! index for derivative direction
  integer,                    intent(in) :: npv ! index for primitive variable
  real(wp), dimension(:,:),   intent(in) :: uu  ! conservative solution
  real(wp), dimension(:,:,:), intent(in) :: dcvdx ! conservative var gradients
  !
  !.. Function Results ..
  real(wp), dimension(1:size(dcvdx,dim=1)) :: return_value
  !
  !.. Local Scalars ..
  integer :: nr,np,n
  real(wp) :: cv_inv,rho_inv,rho_inv2,q2
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(dcvdx,dim=3)) :: coef
  !
continue
  !
  nr = size(dcvdx,dim=2)
  np = size(dcvdx,dim=1)
  !
  if (npv == nec) then
    !
    do n = 1,np
      return_value(n) = dcvdx(n,ndr,npv)
    end do
    !
  else if (npv == nee) then
    !
    cv_inv = (gam - one) / rgasref_nd
    !
    do n = 1,np
      !
      rho_inv = one/uu(n,nec)
      rho_inv2 = rho_inv*rho_inv
      !
      q2 = dot( uu(n,nmb:nme) , uu(n,nmb:nme) )
      !
      coef(nec)     = -rho_inv2 * (uu(n,nee) - rho_inv*q2)
      coef(nmb:nme) = -rho_inv2 * uu(n,nmb:nme)
      coef(nee)     =  rho_inv
      !
      return_value(n) = cv_inv * dot( coef , dcvdx(n,ndr,1:nq) )
      !
    end do
    !
  else
    !
    do n = 1,np
      !
      rho_inv = one/uu(n,nec)
      rho_inv2 = rho_inv*rho_inv
      !
      return_value(n) = rho_inv  * dcvdx(n,ndr,npv) - &
                        rho_inv2 * dcvdx(n,ndr,nec) * uu(n,npv)
      !
    end do
    !
  end if
  !
end function grad_of_selected_pv
!
!###############################################################################
!
pure subroutine interpolate_for_output_1d(array,return_value,nondim)
  !
  !.. Use Statements ..
  use geovar,            only : ncell,cell
  use interpolation_mod, only : outerp
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: array(:)
  real(wp), intent(inout) :: return_value(:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: nondim
  !
  !.. Local Scalars ..
  integer  :: i,k,n,k0,kp,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,val
  !
continue
  !
  ! Interpolate the data
  !
  k0 = 0
  kp = 0
  !
  do nc = 1,ncell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    k0 = k0 + kp
    !
    if (allocated(outerp(this_geom,this_order)%output)) then
      !
      kp = size(outerp(this_geom,this_order)%output%mat,dim=1)
      !
      do k = 1,kp
        !
        val = zero
        !
        do i = outerp(this_geom,this_order)%output%col_ptr(k)+1, &
               outerp(this_geom,this_order)%output%col_ptr(k+1)
          !
          n = outerp(this_geom,this_order)%output%col_idx(i)
          a = outerp(this_geom,this_order)%output%col_val(i)
          !
          val = val + a * array(n0+n)
          !
        end do
        !
        return_value(k0+k) = val
        !
      end do
      !
    else
      !
      kp = np
      !
      do n = 1,np
        return_value(n0+n) = array(n0+n)
      end do
      !
    end if
    !
  end do
  !
  ! Dimensionalize the variable if the nondimensionalization factor is present
  !
  if (present(nondim)) then
    do n = 1,size(return_value)
      return_value(n) = return_value(n) * nondim
    end do
  end if
  !
end subroutine interpolate_for_output_1d
!
!###############################################################################
!
pure subroutine interpolate_for_output_2d(array,return_value,nondim)
  !
  !.. Use Statements ..
  use geovar,            only : ncell,cell
  use interpolation_mod, only : outerp
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: array(:,:)
  real(wp), intent(inout) :: return_value(:,:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: nondim
  !
  !.. Local Scalars ..
  integer  :: i,k,m,n,k0,kp,nc,np,n0,nvar
  integer  :: this_geom,this_order
  real(wp) :: a,val
  !
continue
  !
  nvar = size(array,dim=2)
  !
  ! Interpolate the data
  !
  k0 = 0
  kp = 0
  !
  do nc = 1,ncell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    k0 = k0 + kp
    !
    if (allocated(outerp(this_geom,this_order)%output)) then
      !
      kp = size(outerp(this_geom,this_order)%output%mat,dim=1)
      !
      do m = 1,nvar
        do k = 1,kp
          !
          val = zero
          !
          do i = outerp(this_geom,this_order)%output%col_ptr(k)+1, &
                 outerp(this_geom,this_order)%output%col_ptr(k+1)
            !
            n = outerp(this_geom,this_order)%output%col_idx(i)
            a = outerp(this_geom,this_order)%output%col_val(i)
            !
            val = val + a * array(n0+n,m)
            !
          end do
          !
          return_value(k0+k,m) = val
          !
        end do
      end do
      !
    else
      !
      kp = np
      !
      do m = 1,nvar
        do n = 1,np
          return_value(k0+n,m) = array(n0+n,m)
        end do
      end do
      !
    end if
    !
  end do
  !
  ! Dimensionalize the variables if the nondimensionalization factor is present
  !
  if (present(nondim)) then
    do m = 1,size(return_value,dim=2)
      do n = 1,size(return_value,dim=1)
        return_value(n,m) = return_value(n,m) * nondim
      end do
    end do
  end if
  !
end subroutine interpolate_for_output_2d
!
!###############################################################################
!
pure subroutine interpolate_for_output_3d(array,return_value,nondim)
  !
  !.. Use Statements ..
  use geovar,            only : ncell,cell
  use interpolation_mod, only : outerp
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: array(:,:,:)
  real(wp), intent(inout) :: return_value(:,:,:)
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: nondim
  !
  !.. Local Scalars ..
  integer  :: i,l,k,m,n,k0,kp,nc,np,n0,nvar,ndim
  integer  :: this_geom,this_order
  real(wp) :: a,val
  !
continue
  !
  ndim = size(array,dim=2)
  nvar = size(array,dim=3)
  !
  ! Interpolate the data
  !
  k0 = 0
  kp = 0
  !
  do nc = 1,ncell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    k0 = k0 + kp
    !
    if (allocated(outerp(this_geom,this_order)%output)) then
      !
      kp = size(outerp(this_geom,this_order)%output%mat,dim=1)
      !
      do m = 1,nvar
        do l = 1,ndim
          do k = 1,kp
            !
            val = zero
            !
            do i = outerp(this_geom,this_order)%output%col_ptr(k)+1, &
                   outerp(this_geom,this_order)%output%col_ptr(k+1)
              !
              n = outerp(this_geom,this_order)%output%col_idx(i)
              a = outerp(this_geom,this_order)%output%col_val(i)
              !
              val = val + a * array(n0+n,l,m)
              !
            end do
            !
            return_value(k0+k,l,m) = val
            !
          end do
        end do
      end do
      !
    else
      !
      kp = np
      !
      do m = 1,nvar
        do l = 1,ndim
          do n = 1,np
            return_value(k0+n,l,m) = array(n0+n,l,m)
          end do
        end do
      end do
      !
    end if
    !
  end do
  !
  ! Dimensionalize the variables if the nondimensionalization factor is present
  !
  if (present(nondim)) then
    do m = 1,size(return_value,dim=3)
      do l = 1,size(return_value,dim=2)
        do n = 1,size(return_value,dim=1)
          return_value(n,l,m) = return_value(n,l,m) * nondim
        end do
      end do
    end do
  end if
  !
end subroutine interpolate_for_output_3d
!
!###############################################################################
!
subroutine compute_dmatrix_grad_temperature(dtdx)
  !
  !.. Use Statements
  use geovar,            only : nr,n_solpts,ncell,cell
  use order_mod,         only : n_order,t_order,geom_solpts
  use metrics_mod,       only : metrics_array
  use interpolation_mod, only : interp
  use derivatives_mod,   only : deriv
  use projection_mod,    only : projct
  use flowvar,           only : usp,gpsp
  !
  use thermo_mod, only : update_thermo_properties
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: dtdx(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,l,m,n,nc,np,n0
  integer  :: this_geom,this_order
  integer  :: nqp,num_qpts,q0,ierr
  real(wp) :: a,t,val
  !
  !.. Local Arrays ..
  real(wp) :: dtdi(1:3)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: cv_qp(:,:)
  real(wp), allocatable :: gp_qp(:,:)
  real(wp), allocatable :: temp(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_dmatrix_grad_temperature"
  !
continue
  !
  allocate ( temp(1:n_solpts) , source=zero , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp",1,__LINE__,__FILE__,ierr,error_message)
  !
  if (t_order <= n_order) then
    !
    do n = 1,n_solpts
      temp(n) = gpsp(n,ntemp)
    end do
    !
  else
    !
    ! 1. Interpolate conserved variables to quadrature points and compute
    !    temperature at quadrature points
    !
    nqp = geom_solpts(cell(1)%geom,t_order)
    num_qpts = ncell*nqp
    !
    allocate ( cv_qp(1:num_qpts,1:nq) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cv_qp",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( gp_qp(1:num_qpts,1:ng) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"gp_qp",1,__LINE__,__FILE__,ierr,error_message)
    !
    !
    !
    q0 = 0
    !
    do nc = 1,ncell
      !
      n0 = cell(nc)%beg_sp - 1
      !
      this_geom = cell(nc)%geom   ! geometry type of this cell
      this_order = cell(nc)%order ! solution order of this cell
      !
      ! Now interpolate the conserved variables from the solution points
      ! to the quadrature points
      !
      do n = 1,nqp
        do m = 1,nq
          !
          val = zero
          !
          do i = interp(this_geom,this_order)%toCell(t_order)%col_ptr(n)+1, &
                 interp(this_geom,this_order)%toCell(t_order)%col_ptr(n+1)
            !
            j = interp(this_geom,this_order)%toCell(t_order)%col_idx(i)
            a = interp(this_geom,this_order)%toCell(t_order)%col_val(i)
            !
            val = val + a * usp(n0+j,m)
            !
          end do
          !
          cv_qp(q0+n,m) = val
          !
        end do
        !
      end do
      !
      q0 = q0 + nqp
      !
    end do
    !
    call update_thermo_properties(1,num_qpts,cv_qp,gp_qp,caller_name=pname)
    !
    if (allocated(cv_qp)) then
      deallocate ( cv_qp , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"cv_qp",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    ! 2. Project temperature at quadrature points to solution order polynomial
    !    at solution points
    !
    q0 = 0
    !
    do nc = 1,ncell
      !
      n0 = cell(nc)%beg_sp - 1
      np = cell(nc)%end_sp - n0
      !
      this_geom = cell(nc)%geom   ! geometry type of this cell
      this_order = cell(nc)%order ! solution order of this cell
      !
      ! Now interpolate the derivatives of the conserved variables
      ! from the solution points to the quadrature points
      !
      do n = 1,np
        !
        t = zero
        !
        do i = projct(this_geom,this_order)%CellOI(t_order)%row_ptr(n)+1, &
               projct(this_geom,this_order)%CellOI(t_order)%row_ptr(n+1)
          !
          j = projct(this_geom,this_order)%CellOI(t_order)%row_idx(i)
          a = projct(this_geom,this_order)%CellOI(t_order)%row_val(i)
          !
          t = t + a * gp_qp(q0+j,ntemp)
          !
        end do
        !
        temp(n0+n) = t
        !
      end do
      !
      q0 = q0 + nqp
      !
    end do
    !
    if (allocated(gp_qp)) then
      deallocate ( gp_qp , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"gp_qp",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
  end if
  !
  ! 4. Use derivative matrices to compute gradient of temperature 
  !    in computation space
  !
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1  ! beginning index for solpts in cell nc
    np = cell(nc)%end_sp-n0 ! number of solution points in this cell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Compute the gradient of the conservative variables at each solution point
    !
    do n = 1,np
      !
      do l = 1,nr
        !
        dtdi(l) = zero
        !
        do i = deriv(this_geom,this_order)%d(l)%row_ptr(n)+1, &
               deriv(this_geom,this_order)%d(l)%row_ptr(n+1)
          !
          j = deriv(this_geom,this_order)%d(l)%row_idx(i)
          a = deriv(this_geom,this_order)%d(l)%row_val(i)
          !
          dtdi(l) = dtdi(l) + a * temp(n0+j)
          !
        end do
        !
      end do
      !
      do l = 1,nr
        dtdx(n0+n,l) = dot( dtdi(1:nr) , metrics_array(n0+n,1:nr,l) )
      end do
      !
    end do
    !
  end do
  !
  if (allocated(temp)) then
    deallocate ( temp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"temp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
end subroutine compute_dmatrix_grad_temperature
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
