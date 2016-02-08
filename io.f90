module io_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmb,nme,nmx,nmy,nmz,nee,ntk,ntl,nq
  use ovar,    only : grid_format
  use ovar,    only : gridfile
  use ovar,    only : output_dir
  use ovar,    only : plot3d_cutfile
  use ovar,    only : plot3d_bcsfile
  !
  implicit none
  !
  private
  !
  include "cgnslib_f95.h"
  !
  ! ##################################
  ! ###  Module Public Procedures  ###
  ! ##################################
  !
  public :: read_input_file
  public :: read_grid_file
  public :: write_final_error
  public :: write_parallel_cgns
  public :: write_results_file
  public :: check_flux_point_orderings
  public :: TaylorGreen_8s_vorticity_solution
  !
  ! ###########################
  ! ###  Module Parameters  ###
  ! ###########################
  !
  ! Plot3D grid formats
  !
  character(len=11), parameter :: clas_rk(1:4) = ["RK1-Classic","RK2-Classic", &
                                                  "RK3-Classic","RK4-Classic"]
  character(len=11), parameter :: char_rk(1:5) = ["           ", &
                                                  "RK2-SSP    ","RK3-SSP    ", &
                                                  "RK4-SSP    ","RK4-CK     "]
  character(len=3), parameter :: char_init(-1:2) = ["iRa","iSP","iGa","iLo"]
  !
#ifdef CGNS_ADF
  integer, parameter :: file_type = CG_FILE_ADF
#elif CGNS_ADF2
  integer, parameter :: file_type = CG_FILE_ADF2
#else
  integer, parameter :: file_type = CG_FILE_HDF5
#endif
  !
  logical, parameter :: check_gradients = fals
 !logical, parameter :: check_gradients = true
  !
  ! ################################
  ! ###  Module Local Variables  ###
  ! ################################
  !
  !.. Output File Names ..
  character(len=*),      parameter :: cgns_base_fname = "solution.grid.cgns"
  character(len=170), public, save :: fname_nodal = "solution.nodal.cgns"
  character(len=170), public, save :: fname_modal = "solution.modal.cgns"
  character(len=170), public, save :: error_file = "error.dat"
  character(len=170),         save :: cgns_dir = ""
  character(len=170),         save :: cgns_grid_file
  !
  !.. CGNS Path Variables ..
  integer, save :: ifile_n
  integer, save :: ibase_n
  integer, save :: izone_n
  integer, save :: isect_n
  integer, save :: icoord_n
  integer, save :: isol_n
  integer, save :: ifld_n
  !
  integer, save :: ifile_m
  integer, save :: ibase_m
  integer, save :: izone_m
  integer, save :: isect_m
  integer, save :: icoord_m
  integer, save :: isol_m
  integer, save :: ifld_m
  !
  !.. CGNS Time Dependent Variables ..
  integer, save :: cgnsnodalzone_counter
  integer, save :: cgnsmodalzone_counter
  integer, save :: cgns_sol_counter
  !
  real(sp),          save, allocatable, dimension(:) :: soltime
  real(sp),          save, allocatable, dimension(:) :: modtime
  character(len=32), save, allocatable, dimension(:) :: solname
  character(len=32), save, allocatable, dimension(:) :: modname
  !
  character(len=*), parameter :: base_name = "recon_base"
  character(len=*), parameter :: zone_name = "Solution"
  !
  !.. CGNS Kind Types ..
  integer, save :: cgns_kind_real = RealSingle
  !
  !.. CGNS Parallel Point Indices ..
  integer, save :: n_my_output_pts
  integer, save :: my_pts_beg
  integer, save :: my_pts_end
  !
  integer, dimension(1:3) :: zone_size
  !
  ! DISABLE_CGNS_OUTPUT: Private module variable to disable CGNS output. It
  !                      is currently only relevant to the Taylor-Green Vortex
  !                      problem and is more of a temporary addition until the
  !                      output logic/code becomes more robust.
  !
  logical, save :: disable_cgns_output = fals
  !
  !.. CGNS Data Array ..
  !
  ! NOTE: This is to ensure coordinate and solution variables
  !       are the same real type
  !
  integer, parameter :: cgns_flt = r4
  !
  real(cgns_flt), allocatable, dimension(:) :: var
  !
  interface interpolate_for_output
    module procedure interpolate_for_output_1d, &
                     interpolate_for_output_2d
  end interface interpolate_for_output
  !
contains
!
!###############################################################################
!
subroutine read_input_file
  !
  !.. Local Scalars ..
  integer :: ierr
  character(len=80) :: text
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_input_file"
  character(len=*), parameter :: fname = "recon3.in"
  !
continue
  !
  open (inpt,file=fname,status="old",action="read", &
             form="formatted",iostat=ierr)
  call io_error(pname,fname,1,ierr)
  !
  read (inpt,*,iostat=ierr,err=100,end=100) text
  !
  close (inpt,iostat=ierr)
  call io_error(pname,fname,2,ierr)
  !
  text = trim(adjustl(text))
  text = uppercase( text )
  !
  if (text(1:6) == "&INPUT") then
    call read_nml_input_file
  else
    goto 100
  end if
  !
  return
  !
  ! Read error continuation statements
  !
  100 write (iout,10) ierr ; call stop_recon(stop_mpi,pname)
  !
  ! Format statements
  !
  10 format (/," Error checking the format of the input file 'recon3.in'!",/, &
               " IO Error = ",i0,/)
  !
end subroutine read_input_file
!
!###############################################################################
!
subroutine read_nml_input_file
  !
  !.. Use Statements ..
  use order_mod, only : initialize_geom_solpts
  !
  use geovar, only : nr,using_quadrature
  !
  use input_namelist_mod ! use all the variables in the input namelist
  !
  use ovar,   only : Total_Periods
  use ovar,   only : velref,aref,uref,Lref,time_ref
  use ovar,   only : Period_Timesteps,Period_Time
  use ovar,   only : Grid_Length
  use ovar,   only : LDG_beta,LDG_tau
  !
  use quadrature_mod, only : init_geom_quadrature_rules
  use quadrature_mod, only : init_face_quadrature_rules
  !
  use eqn_idx, only : create_equation_indices
  !
  use initialization_mod, only : initialize_reference
  !
  use channel_mod, only : read_channel_init_file
  !
  !.. Local Scalars ..
  integer  :: it1,it2,ierr,ival
  integer  :: Period_Timesteps_X,Period_Timesteps_Y
  logical  :: use_qdtr_face_pts,mms_error
  real(wp) :: Period_Time_X,Period_Time_Y
  real(wp) :: real_Total_Periods
  real(wp) :: aoa_radians
  !
  !.. Local Arrays ..
  character(len=80), dimension(12) :: input_text
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_nml_input_file"
  character(len=*), parameter :: fname = "recon3.in"
  !
continue
  !
  ! Open the input file
  !
  open (inpt,file=fname,status="old",action="read", &
             form="formatted",iostat=ierr)
  call io_error(pname,fname,1,ierr)
  !
  ! Read in the input file using the input namelist
  !
  read (inpt,nml=input)
  !
  ! Close the input file now that we are finished reading it
  !
  close (inpt,iostat=ierr)
  call io_error(pname,fname,2,ierr)
  !
  ! Check the output directory
  !
  output_dir = trim(adjustl(output_dir))
  call check_output_dir
  !
  ! Check for the CGNS output directory
  !
  call check_cgns_dir
  !
  ! Append the name of the solution and error files to output_dir
  !
  fname_nodal = trim(adjustl(cgns_dir)) // trim(adjustl(fname_nodal))
  fname_modal = trim(adjustl(cgns_dir)) // trim(adjustl(fname_modal))
  error_file = trim(output_dir) // trim(adjustl(error_file))
  !
  ! Adjust the name of the grid file
  !
  gridfile = trim(adjustl(gridfile))
  !
  ! Allow for the use of the variable all_orders to make
  ! n_order, p_order, and q_order all the same value without
  ! having to explicitly specify each one individually.
  !
  if (all_orders >= 0) then
    n_order = all_orders
    p_order = all_orders
    q_order = all_orders
  end if
  !
  ! Allow for the use of the variable all_points to make
  ! iop_solpts, iop_flxpts, and iop_qdtpts all the same value
  ! without having to explicitly specify each one individually.
  !
  if ( any(all_points == [1,2]) ) then
    iop_solpts = all_points
    iop_flxpts = all_points
    iop_qdtpts = all_points
  end if
  !
  if (iop_solpts /= iop_flxpts) then
    if (mypnum == 0) write (iout,101)
    call stop_recon(stop_mpi,pname)
  end if
  !
  if (iop_solpts == Legendre_Gauss_Lobatto) then
    if (continuous_output) then
      write (error_message,106)
      call stop_recon(stop_mpi,pname,error_message)
    end if
  end if
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
  using_quadrature = .not. (q_order==n_order .and. iop_solpts==iop_qdtpts)
  !
  ! Get the number of solution points for each of the different geometry types
  ! This includes quadrature and face values
  !
  use_qdtr_face_pts = using_quadrature .and. dflux_opt==project_derivatives
  !
  call initialize_geom_solpts(use_qdtr_face_pts)
  !
  ! Now that we have the name of the grid file and the location
  ! of the flux and solution points, read in the grid file
  !
  call read_grid_file
  !
  call init_geom_quadrature_rules
  call init_face_quadrature_rules
  !
  ! Create the equation indices now that we know the dimension of the grid
  !
  call create_equation_indices(nr)
  !
  ! Make sure the number of RK stages is correct for some RK types
  !
  if (iRKtype == TVD2_RK) then
    num_stages = 2
  else if (iRKtype == TVD3_RK) then
    num_stages = 3
  else if (iRKtype == TVD4_RK) then
    num_stages = 5
  else if (iRKtype == CK4_RK) then
    num_stages = 5
  end if
  !
  ! Make sure the initialization and error locations are valid
  !
  if (all(loc_init_pts /= Legendre_Nodes)) loc_init_pts = iop_solpts
  if (all(loc_error_pts /= Legendre_Nodes)) loc_error_pts = iop_solpts
  !
  ! Make sure the output order is greater than or
  ! equal to the solution order if it is positive
  !
  if (o_order > 0) o_order = max(o_order,n_order+0)
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
    if (mms_error) call stop_recon(stop_mpi,pname,error_message)
    !
    ! If mms_output_solution is true, set the MMS initialization factor
    ! to one and set the number of time steps and final time to zero 
    !
    if (mms_output_solution) then
      mms_init = one
      num_ts = 0
      Final_Time = zero
    end if
    !
  end if
  !
  if (itestcase == Channel_Flow) then
    call read_channel_init_file
  end if
  !
  ! Initialize the reference conditions
  !
  call initialize_reference
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
      Period_Timesteps_Y = nint(Period_Time_Y/dt_in)
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
      Period_Timesteps_X = nint(Period_Time_X/dt_in)
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
      Period_Timesteps_X = nint(Period_Time_X/dt_in)
      Period_Timesteps_Y = nint(Period_Time_Y/dt_in)
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
    Period_Time = dt_in*real(Period_Timesteps,kind=wp)
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
    num_ts = nint( Final_Time / dt_in )
    !
  else if (itestcase == Taylor_Green_Vortex) then
    !
    cfl_beg = CFL
    cfl_end = CFL
    cfl_cycles = 1
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
      Final_Time = Final_Time / time_ref
      !
      num_ts = nint( Final_Time / dt_in )
      !
    end if
    !
    ! Make sure the timestep type is global (positive)
    !
    Timestep_Type = abs(Timestep_Type)
    !
  else
    !
    if (Timestep_Type == Constant_Timestep) then
      !
      cfl_beg = CFL
      cfl_end = CFL
      cfl_cycles = 1
      !
      if (Final_Time > zero) then
        num_ts = nint( Final_Time / dt_in )
      end if
      !
    else
      !
      cfl_end = max(cfl_end,eps6)    ! cfl_end must be >= 1.0e-6
      cfl_cycles = max(cfl_cycles,1) ! cfl_cycles must be at least 1
      cfl_beg = min(cfl_beg,cfl_end) ! cfl_beg must be <= cfl_end
      !
    end if
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
      results_interval = num_ts/500
    else if (results_interval < 0) then
      results_interval = num_ts/abs(results_interval)
    end if
    results_interval = max(results_interval,1)
    !
    if (itestcase == Taylor_Green_Vortex) then
      if (output_interval == 0) then
        disable_cgns_output = true
      end if
    end if
    !
    if (output_interval == 0) then
      output_interval = num_ts/50
    else if (output_interval < 0) then
      output_interval = num_ts/abs(output_interval)
    end if
    output_interval = max(output_interval,1)
    !
    if (iter_out_interval == 0) then
      iter_out_interval = 1
    else if (iter_out_interval < 0) then
      iter_out_interval = num_ts/abs(iter_out_interval)
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
  if (mypnum == 0) then
    !
    write (iout,104)
    write (iout,nml=input)
    write (iout,105)
    !
  end if
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
  21 format ("MMS ERROR: MMS solutions 1 and 2 are allowed for 2D and 3D only!")
  22 format ("MMS ERROR: MMS solution 3 is allowed for 1D and 2D only!")
  23 format ("MMS ERROR: MMS solution 4 is allowed for 2D only!")
  24 format ("MMS ERROR: The MMS solution chosen is invalid!")
  !
end subroutine read_nml_input_file
!
!###############################################################################
!
subroutine check_output_dir
  !
  integer :: clen,ierr,iodir,ii
  character(len=200) :: test_file
  character(len=200) :: test_dir
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_output_dir"
  !
  character(len=*), parameter :: nobackup = "/nobackup/sspiegel"
  character(len=*), parameter :: homedir = "/u/sspiegel"
  !
continue
  !
  ! Initialize test_file and test_dir
  !
  test_file = " "
  test_dir = " "
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
      clen = len_trim(output_dir)
      if (output_dir(clen:clen) /= "/") then
        output_dir = trim(output_dir)//"/"
      end if
      !
      ! Try to open a file in the given directory
      !
      test_file = trim(output_dir)//"output_dir.test_file.dat"
      !
      open (newunit=iodir,file=trim(test_file),status="replace", &
                          action="write",iostat=ierr)
      !
      ! If the open was successful, close and delete the test file and return
      !
      if (ierr == 0) then
        close (iodir,status="delete",iostat=ierr)
        exit one_pass
      end if
      !
      ! Reset test_file
      !
      test_file = " "
      !
      ! Try to append nobackup to the beginning of output_dir
      ! and try a test file in that directory
      !
      if (output_dir(1:len(nobackup)) /= nobackup) then
        !
        if (output_dir(1:1) /= "/") then
          test_dir = nobackup//"/"//trim(output_dir)
        else
          test_dir = nobackup//trim(output_dir)
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
          output_dir = trim(test_dir)
          exit one_pass
        end if
        !
      end if
      !
      ! Reset test_file and test_dir
      !
      test_file = " "
      test_dir = " "
      !
      ! Finally, try to append homedir to the beginning of output_dir
      ! and try a test file in that directory
      !
      if (output_dir(1:len(homedir)) /= homedir) then
        !
        if (output_dir(1:1) /= "/") then
          test_dir = homedir//"/"//trim(output_dir)
        else
          test_dir = homedir//trim(output_dir)
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
          output_dir = trim(test_dir)
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
  !
  clen = len_trim(output_dir)
  call mpi_bcast(clen,1,mpi_inttyp,0,mpi_comm_world,mpierr)
  call mpi_bcast(output_dir(1:clen),clen,MPI_CHARACTER,0,mpi_comm_world,mpierr)
  !
  ! If output_dir starts with "BAD_DIR", we cant seem to find the
  ! output directory that was given in the input file so exit
  !
  if (output_dir(1:7) == "BAD_DIR") then
    output_dir = output_dir(8:clen)
    if (mypnum == 0) write (iout,1) trim(output_dir)
    call stop_recon(stop_mpi,pname)
  end if
  !
  1 format (/,"******************** ERROR! *********************",/, &
              "     CANT SEEM TO FIND THE OUTPUT DIRECTORY",/, &
              "       THAT WAS GIVEN IN THE INPUT FILE.",/, &
              "      '",a,"'",/, &
              "     PLEASE CHECK THE INPUT FILE TO CORRECT",/, &
              "        ANY MISTAKES BEFORE CONTINUING!",/, &
              "******************** ERROR! *********************",/)
  !
end subroutine check_output_dir
!
!###############################################################################
!
subroutine check_cgns_dir
  !
  integer :: clen,ierr,iodir
  integer :: exitstat,cmndstat
  logical :: using_cgns_dir
  character(len=300) :: command
  character(len=300) :: cmndmsg
  character(len=200) :: test_file
  character(len=200) :: test_dir
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_cgns_dir"
  !
  logical, parameter :: wait_for_completion = true
  !
continue
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
      call execute_command_line( trim(command) , wait_for_completion , &
                                 exitstat , cmndstat , cmndmsg )
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
      command = "lfs setstripe -c 64 -s 4m "//trim(cgns_dir)
      !
      call execute_command_line( trim(command) , wait_for_completion , &
                                 exitstat , cmndstat , cmndmsg )
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
  call mpi_bcast(clen,1,mpi_inttyp,0,mpi_comm_world,mpierr)
  call mpi_bcast(cgns_dir(1:clen),clen,MPI_CHARACTER,0,mpi_comm_world,mpierr)
  !
end subroutine check_cgns_dir
!
!###############################################################################
!
subroutine read_grid_file()
  !
  !.. Use Statements ..
  use module_gmsh,    only : read_gmsh_gridfile
  use module_plot3d,  only : read_plot3d_gridfile
  use quadrature_mod, only : find_valid_geometries
  use ovar,           only : use_bc_conflicts
  use ovar,           only : itestcase,TaylorGreen_DOF
  !
  ! Routine that uses the variable grid_format to choose
  ! which routine to use in order to read the grid file.
  !        
continue
  !
  call memory_pause("reading input file","reading grid file")
  if (itestcase == Taylor_Green_Vortex .and. TaylorGreen_DOF > 0) then
    call create_TaylorGreenVortex_grid
  else if (grid_format == Gambit_Neutral) then
    call read_gambit_neutral_gridfile
  else if (grid_format == Gmsh_Format) then
    call read_gmsh_gridfile( gridfile )
  else if (grid_format == Plot3D_Format) then
    call read_plot3d_gridfile( gridfile , plot3d_cutfile , plot3d_bcsfile )
  end if
  !
  ! Check that all the nodes are in counter-clockwise
  ! ordering for each cell and adjust if needed
  !
  call memory_pause("reading grid file","checking the node orderings")
  call check_node_ordering
  !
  ! Adjust the grid coordinates and nodes_of_cell
  ! connectivity to account for ghost cells
  !
  call memory_pause("checking the node orderings","adjusting for ghost cells")
  call adjust_for_ghost_cells
  !
  ! Find any grid nodes on the boundaries that separate two boundary
  ! faces with non-matching boundary conditions
  !
  if (use_bc_conflicts) then
    call memory_pause("adjusting for ghost cells", &
                      "finding boundary condition conflicts")
    call find_bc_conflicts
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
  call memory_pause("finding valid geometries", &
                    "finish processing the remainder of the input file")
  !
end subroutine read_grid_file
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
  use ovar,   only : iop_solpts,iop_flxpts
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
  logical, allocatable, dimension(:) :: bc_done
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
  ! Open the grid file
  !
  open (newunit=iogrd,file=gridfile,status="old",action="read", &
                      position="rewind",iostat=ierr)
  call io_error(pname,trim(adjustl(gridfile)),1,ierr)
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
  allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , stat=ierr )
  call alloc_error(pname,"nodes_of_cell_ptr",1,ierr)
  !
  allocate ( cell_geom(1:ncell) , source=0 , stat=ierr )
  call alloc_error(pname,"cell_geom",1,ierr)
  !
  allocate ( cell_order(1:ncell) , source=n_order , stat=ierr )
  call alloc_error(pname,"cell_order",1,ierr)
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
  allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , stat=ierr )
  call alloc_error(pname,"xyz_nodes",1,ierr)
  allocate ( nodes_of_cell(1:last(nodes_of_cell_ptr)) , source=0 , stat=ierr )
  call alloc_error(pname,"nodes_of_cell",1,ierr)
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , stat=ierr )
  call alloc_error(pname,"bface",1,ierr)
  allocate ( bc_done(1:nfbnd) , source=true , stat=ierr )
  call alloc_error(pname,"bc_done",1,ierr)
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
  close (iogrd,iostat=ierr)
  call io_error(pname,trim(adjustl(gridfile)),2,ierr)
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
  deallocate ( bc_done , stat=ierr )
  call alloc_error(pname,"bc_done",2,ierr)
  !
  ! Make sure that Lobatto nodes are requested if the grid contains triangle
  ! elements since Gauss nodes is not currently set up for triangles.
  !
  if (any(cell_geom == Geom_Tria)) then
    write (*,3)
    call stop_recon(stop_mpi,pname)
  end if
  !
  if (iop_solpts /= Legendre_Gauss_Lobatto .or. &
      iop_flxpts /= Legendre_Gauss_Lobatto) then
    if (any(cell_geom == Geom_Tria)) then
      write (*,2)
      call stop_recon(stop_mpi,pname)
    end if
  end if
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
subroutine create_TaylorGreenVortex_grid()
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use ovar,      only : TaylorGreen_DOF,Lref
  use ovar,      only : TaylorGreen_random_perturb
  use ovar,      only : TaylorGreen_perturb_percent
  use geovar,    only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geovar,    only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geovar,    only : cell_geom,cell_order
  !
  !.. Local Scalars ..
  integer :: nc,np,ierr,i,j,k,n
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
  nc = nint( real(TaylorGreen_DOF,kind=wp) / real(n_order+1,kind=wp) )
  !
  np = nc + 1 ! number of grid nodes in each direction
  !
  ncell = nc*nc*nc ! total number of grid cells
  nnode = np*np*np ! total number of grid nodes
  nfbnd = 6*nc*nc  ! total number of boundary faces
  !
  ! Allocate all the arrays needed to define the grid
  !
  allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , stat=ierr )
  call alloc_error(pname,"xyz_nodes",1,ierr)
  !
  allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , stat=ierr )
  call alloc_error(pname,"nodes_of_cell_ptr",1,ierr)
  !
  allocate ( nodes_of_cell(1:8*ncell) , source=0 , stat=ierr )
  call alloc_error(pname,"nodes_of_cell",1,ierr)
  !
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , stat=ierr )
  call alloc_error(pname,"bface",1,ierr)
  !
  allocate ( cell_geom(1:ncell) , source=Geom_Hexa , stat=ierr )
  call alloc_error(pname,"cell_geom",1,ierr)
  !
  allocate ( cell_order(1:ncell) , source=n_order , stat=ierr )
  call alloc_error(pname,"cell_order",1,ierr)
  !
  allocate ( perturb(1:nr,1:nnode) , source=real(0,kind=r4) , stat=ierr )
  call alloc_error(pname,"perturb",1,ierr)
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
    nodes_of_cell_ptr(n) = nodes_of_cell_ptr(n-1) + 8
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
  if (TaylorGreen_random_perturb) then
    !
    rbas = dl * real(abs(TaylorGreen_perturb_percent),kind=r4) &
              / real(100,kind=r4)
    radj = rbas / real(2,kind=r4)
    !
    if (TaylorGreen_perturb_percent < 0) then
      !
      ! Randomly perturb the crap out of everything
      !
      do k = 2,np-1
        do j = 2,np-1
          do i = 2,np-1
            !
            n = map3dto1d(i,j,k,np)
            !
            rval = random(rbas) - radj
            perturb(1,n) = rval
            !
            rval = random(rbas) - radj
            perturb(2,n) = rval
            !
            rval = random(rbas) - radj
            perturb(3,n) = rval
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
  deallocate ( perturb , stat=ierr )
  call alloc_error(pname,"perturb",2,ierr)
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
      nf_min = nf_min + 1
      nf_max = nf_min + nc*nc
      !
      bface( 1,nf_min) = bc_periodic
      bface( 2,nf_min) = map3dto1d(1 ,j,k,nc)
      bface( 3,nf_min) = Geom_Hexa
      bface( 4,nf_min) = HostFace3D(1)
      bface( 5,nf_min) = nf_max
      bface( 6,nf_min) = map3dto1d(nc,j,k,nc)
      bface( 8,nf_min) = Geom_Hexa
      bface( 9,nf_min) = HostFace3D(2)
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
      bface( 1,nf_min) = bc_periodic
      bface( 2,nf_min) = map3dto1d(i,1 ,k,nc)
      bface( 3,nf_min) = Geom_Hexa
      bface( 4,nf_min) = HostFace3D(3)
      bface( 5,nf_min) = nf_max
      bface( 6,nf_min) = map3dto1d(i,nc,k,nc)
      bface( 8,nf_min) = Geom_Hexa
      bface( 9,nf_min) = HostFace3D(4)
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
      bface( 1,nf_min) = bc_periodic
      bface( 2,nf_min) = map3dto1d(i,j,1 ,nc)
      bface( 3,nf_min) = Geom_Hexa
      bface( 4,nf_min) = HostFace3D(5)
      bface( 5,nf_min) = nf_max
      bface( 6,nf_min) = map3dto1d(i,j,nc,nc)
      bface( 8,nf_min) = Geom_Hexa
      bface( 9,nf_min) = HostFace3D(6)
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
        call stop_recon(stop_mpi,pname)
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
        call stop_recon(stop_mpi,pname)
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
        call stop_recon(stop_mpi,pname)
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
        call stop_recon(stop_mpi,pname)
        !
      end if
    end if
    !
  end do
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
      call stop_recon(stop_mpi,pname,error_message)
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
      call stop_recon(stop_mpi,pname,error_message)
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
  n_global_node = nnode
  n_global_cell = ncell
  n_global_fbnd = nfbnd
  n_global_totnod = n_totnod
  n_global_totcel = n_totcel
  !
  n1 = size(xyz_nodes,dim=1)
  n2 = size(xyz_nodes,dim=2)
  allocate ( global_xyz_nodes(1:n1,1:n2) , source=xyz_nodes , stat=ierr )
  call alloc_error(pname,"global_xyz_nodes",1,ierr)
  !
  n1 = size(nodes_of_cell_ptr)
  allocate ( global_nofc_ptr(1:n1) , source=nodes_of_cell_ptr , stat=ierr )
  call alloc_error(pname,"global_nofc_ptr",1,ierr)
  !
  n1 = size(nodes_of_cell)
  allocate ( global_nofc(1:n1) , source=nodes_of_cell , stat=ierr )
  call alloc_error(pname,"global_nofc",1,ierr)
  !
  n1 = size(cell_geom)
  allocate ( global_cell_geom(1:n1) , source=cell_geom , stat=ierr )
  call alloc_error(pname,"global_cell_geom",1,ierr)
  !
  n1 = size(cell_order)
  allocate ( global_cell_order(1:n1) , source=cell_order , stat=ierr )
  call alloc_error(pname,"global_cell_order",1,ierr)
  !
  ! Create the global integer pointer for the beginning and ending
  ! indices for the solution points in each cell. This will be done
  ! later but a global version is never made before the global grid
  ! is partitioned and each processor creates its own separate
  ! local grid.
  !
  allocate ( global_pinc_ptr(1:n_global_totcel+1) , source=0 , stat=ierr )
  call alloc_error(pname,"global_pinc_ptr",1,ierr)
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
  integer  :: g1,g2,i,h1,h2,nf,ng,ni,np
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
      call stop_recon(stop_mpi,pname,error_message)
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
        ng = nodes_of_cell( g1 + 2 )
        !
        ! Store the coordinates for node ni
        !
        xyz_ni(:) = xyz_nodes(:,ni)
        !
        ! Get the coordinates of ghost node 3
        !
        xyz_nodes(:,ng) = merge(rotate(xyz_h1(:),xyz_h2(:),xyz_ni(:)), &
                                reflect(xyz_ni(:),ptf,rnrm), &
                                rotation == 1)
        !
      else if (host_geom == 4) then
        !
        ! Ghost node 3 for quadrilateral
        !
        ng = nodes_of_cell( g1 + 2 )
        ni = nodes_of_cell( h1 + &
                            mod(host_side+2-max(0,min(1,rotation)),host_geom) )
        !
        ! Store the coordinates for node ni
        !
        xyz_ni(:) = xyz_nodes(:,ni)
        !
        ! Get the coordinates of ghost node 3
        !
        xyz_nodes(:,ng) = merge(rotate(xyz_h1(:),xyz_h2(:),xyz_ni(:)), &
                                reflect(xyz_ni(:),ptf,rnrm), &
                                rotation == 1)
        !
        ! Ghost node 4 for quadrilateral
        !
        ng = nodes_of_cell( g1 + 3 )
        ni = nodes_of_cell( h1 + &
                            mod(host_side+1+max(0,min(1,rotation)),host_geom) )
        !
        ! Store the coordinates for node ni
        !
        xyz_ni(:) = xyz_nodes(:,ni)
        !
        ! Get the coordinates of ghost node 4
        !
        xyz_nodes(:,ng) = merge(rotate(xyz_h1(:),xyz_h2(:),xyz_ni(:)), &
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
        call stop_recon(stop_mpi,pname,error_message)
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
     !  ng = nodes_of_cell( g1 + i-1 )
     !  np = nodes_of_cell( p1 + mod(periodic_side+i-2,pcnods) )
     !  !
     !  xyz_nodes(:,ng) = xyz_nodes(:,np) + ds1(:)
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
  logical, allocatable, dimension(:) :: mskp
  integer, allocatable, dimension(:) :: faces_with_node_ptr
  integer, allocatable, dimension(:) :: faces_with_node
  integer, allocatable, dimension(:,:) :: bface_nodes
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_bc_conflicts"
  !
continue
  !
  error_count = 0
  !
  allocate ( bface_nodes(1:4,1:nfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface_nodes",1,ierr,error_message)
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
    call stop_recon(stop_mpi,pname,error_message)
    !
  end if
  !
  ! Mark the nodes defining each boundary face, storing ahead
  ! NOTE: Skip all communication type boundary faces
  !
  allocate ( faces_with_node_ptr(1:nnode+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faces_with_node_ptr",1,ierr,error_message)
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
  call alloc_error(pname,"faces_with_node",1,ierr,error_message)
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
  call alloc_error(pname,"mskp",1,ierr,error_message)
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
  call alloc_error(pname,"bc_conflict",1,ierr,error_message)
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
    call alloc_error(pname,array_name,1,ierr,error_message)
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
      call stop_recon(stop_mpi,pname,error_message)
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
        call stop_recon(stop_mpi,pname,error_message)
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
  call alloc_error(pname,"bface_nodes",2,ierr,error_message)
  deallocate ( faces_with_node_ptr , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faces_with_node_ptr",2,ierr,error_message)
  deallocate ( faces_with_node , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faces_with_node",2,ierr,error_message)
  deallocate ( mskp , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"mskp",2,ierr,error_message)
  !
  ! Format Statements
  !
  1 format (a,i0,a)
  2 format (//," Number of BC Conflicts = ",i0,/)
  3 format (" Conflict #",i0,": Node = ",i0,"; Faces = ",*(i0,:,", "))
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
subroutine check_flux_point_orderings()
  !
  !.. Use Statements ..
  use order_mod, only : maxSP,geom_solpts
  use geovar, only : ncell,cell
  use geovar, only : nfbnd,nface,face
  use geovar, only : nr,xyz
  use interpolation_mod, only : interp
  !
  !.. Local Scalars ..
  integer :: nf,nfp,k,kl,kr,l
  integer :: l1,l2,r1,r2,npl,npr
  integer :: face_geom,face_order
  integer :: left_cell,left_geom,left_order
  integer :: right_cell,right_geom,right_order
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: xyz_fp
#ifndef USE_MATRIX_PROCEDURES
  real(wp), dimension(1:maxSP) :: lmat
#endif
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_flux_point_orderings"
  !
continue
  !
  write (151,*)
  write (152,*)
  !
  do nf = nfbnd+1,nface
    !
    write (151,2) nf-nfbnd
    write (152,2) nf-nfbnd
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
    do k = 1,nfp
      !
      kl = face(nf)%left%flxpts(k)
      kr = face(nf)%right%flxpts(k)
      !
#ifdef USE_MATRIX_PROCEDURES
      do l = 1,nr
        xyz_fp(l) = &
        interp(left_geom,left_order)%toFace(face_order)%dot( kl , xyz(l,l1:l2) )
      end do
#else
      lmat(1:npl) = &
        interp(left_geom,left_order)%toFace(face_order)%mat(1:npl,kl)
      do l = 1,nr
        xyz_fp(l) = dot( lmat(1:npl) , xyz(l,l1:l2) )
      end do
#endif
      !
      write (151,1) kl,xyz_fp(1:nr)
      !
#ifdef USE_MATRIX_PROCEDURES
      do l = 1,nr
        xyz_fp(l) = &
    interp(right_geom,right_order)%toFace(face_order)%dot( kr , xyz(l,r1:r2) )
      end do
#else
      lmat(1:npr) = &
        interp(right_geom,right_order)%toFace(face_order)%mat(1:npr,kr)
      do l = 1,nr
        xyz_fp(l) = dot( lmat(1:npr) , xyz(l,r1:r2) )
      end do
#endif
      !
      write (152,1) kr,xyz_fp(1:nr)
      !
    end do
    !
  end do
  !
 !call stop_recon(stop_mpi,pname)
  !
  ! Format Statements
  !
 !1 format (i4,*(es23.15))
  1 format (i4,*(es12.4))
  2 format (/,"#######################################################",/, &
              "#################  FACE ",i0,"  ###########################",/, &
              "#######################################################",/)
  !
end subroutine check_flux_point_orderings
!
!###############################################################################
!
subroutine write_final_error
  !
  use geovar, only : nr,n_global_cell
  use geovar, only : n_global_solpts,global_pinc_ptr
  use ovar,   only : itestcase,flux_divergence
  use ovar,   only : vnloc,vkloc,vlone,vlinf,vltwo
  use ovar,   only : unloc,ukloc,ulone,ulinf,ultwo
  use order_mod, only : n_order
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_final_error"
  !
continue
  !
  write (iout,31) itestcase,flux_divergence, &
                  n_order,n_global_cell,n_global_solpts, &
                  vnloc(1),global_pinc_ptr(vnloc(1))+vkloc(1)
  !
  write (iout,32) vlone(nec),  vlinf(nec),  vltwo(nec)
  write (iout,33) vlone(nmx),  vlinf(nmx),  vltwo(nmx)
  if (nr > 1) &
  write (iout,34) vlone(nmy),  vlinf(nmy),  vltwo(nmy)
  if (nr > 2) &
  write (iout,35) vlone(nmz),  vlinf(nmz),  vltwo(nmz)
  write (iout,36) ulone(nq+1), ulinf(nq+1), ultwo(nq+1)
  write (iout,37) vlone(nq+1), vlinf(nq+1), vltwo(nq+1)
  !
  31 format (/,'itestcase =',i2,/, &
               'flux_divergence  =',i2,//, &
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
subroutine write_results_file(iflag,icount)
  !
  use geovar,  only : nr,n_global_solpts
  use ovar,    only : ultwo,vltwo,work_units
  use order_mod, only : p_order
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iflag
  integer, intent(in) :: icount
  !
  !.. Local Scalars ..
  integer  :: i,ierr
  real(wp) :: comp_cost
  !
  !.. Local Saved Scalars ..
  integer, save :: iores
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nq+1) :: u2,v2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_results_file"
  character(len=*), parameter :: fname = "results.dat"
  !
continue
  !
  ! Open the results file only if iflag is +1
  !
  if (iflag == 1) then
    !
    open (newunit=iores,file=fname,status="replace",action="write",iostat=ierr)
    call io_error(pname,fname,1,ierr)
    !
    if (nr == 1) then
      write (iores,1) p_order, one/sqrt(real(n_global_solpts,kind=wp))
    else if (nr == 2) then
      write (iores,2) p_order, one/sqrt(real(n_global_solpts,kind=wp))
    else if (nr == 3) then
      write (iores,3) p_order, one/sqrt(real(n_global_solpts,kind=wp))
    end if
    !
    return
    !
  end if
  !
  comp_cost = log10(merge(work_units,one,work_units > zero))
  !
  u2(:) = ultwo(:)
  v2(:) = vltwo(:)
  !
  where (u2 <= zero)
    u2 = one
  end where
  where (v2 <= zero)
    v2 = one
  end where
  !
  write (iores,4) icount,(log10(u2(i)),i=1,nq+1), &
                         (log10(v2(i)),i=2,nq+1), &
                         comp_cost
  !
  ! Close the results file only if iflag is -1
  !
  if (iflag == -1) then
    close (iores,iostat=ierr)
    call io_error(pname,fname,2,ierr)
  end if
  !
  ! Format Statements
  !
  1 format ('VARIABLES= "Timestep" "rho" "rhou" "Energy" &
           &"Entropy" "u" "Pressure" "||V||" "Work units"',/, &
            'Zone T="Residual Error"',/, &
            'AUXDATA order="',i0,'"',/, &
            'AUXDATA h="',f10.6,'"')
  2 format ('VARIABLES= "Timestep" "rho" "rhou" "rhov" "Energy" &
           &"Entropy" "u" "v" "Pressure" "||V||" "Work units"',/, &
            'Zone T="Residual Error"',/, &
            'AUXDATA order="',i0,'"',/, &
            'AUXDATA h="',f10.6,'"')
  3 format ('VARIABLES= "Timestep" "rho" "rhou" "rhov" "rhow" "Energy" &
           &"Entropy" "u" "v" "w" "Pressure" "||V||" "Work units"',/, &
            'Zone T="Residual Error"',/, &
            'AUXDATA order="',i0,'"',/, &
            'AUXDATA h="',f10.6,'"')
  4 format (i10,*(es17.9))
  !
end subroutine write_results_file
!
!###############################################################################
!
subroutine write_parallel_cgns(iflag)
  !
  !.. Use Statements ..
  use geovar,  only : nr,n_solpts,xyz
  use ovar,    only : time,aref,uref,rhoref,time_ref
  use ovar,    only : itestcase,cgns_use_queue
  use ovar,    only : interpolate_before_output_variable
  use ovar,    only : write_TaylorGreen_full_solution
  use ovar,    only : completely_disable_cgns
  use flowvar, only : usp,dusp
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iflag
  !
  !.. Local Scalars ..
  integer :: n,ierr
  integer :: my_cells_beg,my_cells_end
  integer :: n_total_output_pts
  integer :: n_total_output_cells
  !
  character(len=32) :: tname
  character(len=32) :: varname
  character(len=170) :: sol_file
  character(len=300) :: link_path
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_parallel_cgns"
  character(len=*), parameter :: char_xyz(1:3) = ["X","Y","Z"]
  !
  !.. Local Arrays ..
  integer, dimension(1:2) :: idata
  integer, dimension(1:nr) :: ixyz_n
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable, dimension(:,:) :: ipelem
  real(wp), allocatable, dimension(:)   :: vartmp
  !
continue
  !
  if (completely_disable_cgns) return
  !
  if (disable_cgns_output) then
    if (iflag == 0) then
      if (.not. write_TaylorGreen_full_solution) then
        return
      end if
    end if
  end if
  !
  if (iflag == 1) then
    !
    cgns_grid_file = trim(adjustl(cgns_dir)) // cgns_base_fname
    !
    cgns_kind_real = merge(RealDouble,RealSingle,kind(var)==wp)
    !
    ! Initialize the variables and data arrays needed for CGNS output
    ! NOTE: This call to routine get_initial_data_for_cgns will initialize
    !       the module variables my_pts_beg, my_pts_end, and n_my_output_pts.
    !
    call get_initial_data_for_cgns(my_cells_beg,my_cells_end,ipelem)
    !
    ! Extract the total number of pts and cells from the zone_size array
    !
    n_total_output_pts = zone_size(1)
    n_total_output_cells = zone_size(2)
    !
    ! Allocate the temporary variable arrays
    !
    allocate ( vartmp(1:n_solpts) , source=zero , stat=ierr )
    call alloc_error(pname,"vartmp",1,ierr)
    allocate ( var(1:n_my_output_pts) , source=real(0,cgns_flt) , stat=ierr )
    call alloc_error(pname,"var",1,ierr)
    !
    ! Open the cgns output file if iflag is 1
    !
    ! Initialize the zone counter
    !
    cgnsnodalzone_counter = 0
    !
    call cgp_open_f(cgns_grid_file,CG_MODE_WRITE,ifile_n,ierr)
    call io_error(pname,cgns_grid_file,1,ierr)
    !
    ! Set the CGNS base (who knows what this is for)
    !
    call cg_base_write_f(ifile_n,base_name,nr,nr,ibase_n,ierr)
    call cgns_error(pname,cgns_grid_file,"cg_base_write_f",ierr)
    !
    ! Write the title for this solution
    !
    call cg_zone_write_f(ifile_n,ibase_n,zone_name,zone_size, &
                         Unstructured,izone_n,ierr)
    call cgns_error(pname,cgns_grid_file,"cg_zone_write_f",ierr)
    !
    ! Create data nodes for the coordinates
    !
    do n = 1,nr
      varname = "Coordinate"//char_xyz(n)
      call cgp_coord_write_f(ifile_n,ibase_n,izone_n,cgns_kind_real, &
                             trim(varname),ixyz_n(n),ierr)
      call cgns_error(pname,cgns_grid_file,"cgp_coord_write_f",ierr, &
                      "Write node: Coordinate",n)
    end do
    !
    ! Now create and write the coordinate data in parallel
    !
    if (.not. cgns_use_queue) then
      call cgp_queue_set_f(0,ierr)
      call cgns_error(pname,cgns_grid_file,"cgp_queue_set_f",ierr)
    end if
    !
    do n = 1,nr
      !
      ! Save the current coordinate to the local array vartmp
      !
      vartmp(1:n_solpts) = xyz(n,1:n_solpts)
      !
      ! Interpolate vartmp to the desired output locations and
      ! convert to the real type required for the CGNS file.
      !
      var(:) = real( interpolate_for_output(vartmp) , kind(var) )
      !
      ! Reset the queue.
      !
      if (cgns_use_queue) then
        call cgp_queue_set_f(1,ierr)
        call cgns_error(pname,cgns_grid_file,"cgp_queue_set_f",ierr)
      end if
      !
      ! Add the data for current coordinate variable to the output queue.
      !
      call cgp_coord_write_data_f(ifile_n,ibase_n,izone_n,ixyz_n(n), &
                                  my_pts_beg,my_pts_end,var,ierr)
      call cgns_error(pname,cgns_grid_file,"cgp_coord_write_data_f",ierr, &
                      "Write data: Coordinate",n)
      !
      ! Finally, flush the queue containing the
      ! current coordinate variable to the CGNS file.
      !
      if (cgns_use_queue) then
        call cgp_queue_flush_f(ierr)
        call cgns_error(pname,cgns_grid_file,"cgp_queue_flush_f",ierr)
      end if
      !
    end do
    !
    ! Create the data node for the grid connectivity
    !
    if (nr == 2) then
      call cgp_section_write_f(ifile_n,ibase_n,izone_n,"Quadrilaterals", &
                               QUAD_4,1,n_total_output_cells,0,isect_n,ierr)
    else if (nr == 3) then
      call cgp_section_write_f(ifile_n,ibase_n,izone_n,"Hexahedrals", &
                               HEXA_8,1,n_total_output_cells,0,isect_n,ierr)
    end if
    call cgns_error(pname,cgns_grid_file,"cgp_section_write_f",ierr)
    !
    ! Write the element connectivity in parallel
    !
    ! Reset the queue.
    !
    if (cgns_use_queue) then
      call cgp_queue_set_f(1,ierr)
      call cgns_error(pname,cgns_grid_file,"cgp_queue_set_f",ierr)
    end if
    !
    ! Add the element connectivity data to the output queue.
    !
    call cgp_elements_write_data_f(ifile_n,ibase_n,izone_n,isect_n, &
                                   my_cells_beg,my_cells_end,ipelem,ierr)
    call cgns_error(pname,cgns_grid_file,"cgp_elements_write_data_f",ierr)
    !
    ! Finally, flush the queue containing the
    ! element connectivity to the CGNS file.
    !
    if (cgns_use_queue) then
      call cgp_queue_flush_f(ierr)
      call cgns_error(pname,cgns_grid_file,"cgp_queue_flush_f",ierr)
    end if
    !
    ! Deallocate the arrays var, vartmp, and ipelem
    !
    if (allocated(var)) then
      deallocate ( var , stat=ierr )
      call alloc_error(pname,"var",2,ierr)
    end if
    if (allocated(vartmp)) then
      deallocate ( vartmp , stat=ierr )
      call alloc_error(pname,"vartmp",2,ierr)
    end if
    if (allocated(ipelem)) then
      deallocate ( ipelem , stat=ierr )
      call alloc_error(pname,"ipelem",2,ierr)
    end if
    !
  else!if (iflag == 2) then
    !
    ! The CGNS file was dumped last time this subroutine was called
    ! so the CGNS file needs to be reopened in modification mode.
    !
    call cgp_open_f(cgns_grid_file,CG_MODE_MODIFY,ifile_n,ierr)
    call io_error(pname,cgns_grid_file,3,ierr)
    !
  end if
  !
  ! Get the name for the current solution node
  !
  write (tname,'("FlowSolution",i0)') cgnsnodalzone_counter
  write (sol_file,'("solution.",i0,".cgns")') cgnsnodalzone_counter
  !
  ! Reallocate the soltime and solname arrays so we
  ! can concatenate this solution onto those arrays
  !
  call reallocate(soltime,cgnsnodalzone_counter+1)
  call reallocate(solname,cgnsnodalzone_counter+1)
  soltime(size(soltime)) = real(time*time_ref,kind(soltime))
  solname(size(solname)) = tname
  !
  ! Create a link to the solution node in the current solution file
  !
  call cg_goto_f(ifile_n,ibase_n,ierr,"Zone_t",izone_n,"end")
  call cgns_error(pname,cgns_grid_file,"cg_goto_f",ierr)
  !
  link_path = base_name//"/"//zone_name//"/"//trim(tname)
  !
  call cg_link_write_f(trim(tname),trim(sol_file),trim(link_path),ierr)
  call cgns_error(pname,cgns_grid_file,"cg_link_write_f",ierr)
  !
  ! Write the solution file
  !
  call write_solution_parallel_cgns(sol_file,tname)
  !
  ! ################################################
  ! #####   FINISHED WRITING FIELD VARIABLES   #####
  ! ################################################
  !
  if (cgns_use_queue) then
    call cgp_queue_set_f(0,ierr)
    call cgns_error(pname,cgns_grid_file,"cgp_queue_set_f",ierr)
  end if
  !
  ! If last iteration or a solution dump was requested, close CGNS file
  !
 !if (any(iflag == [-1,-2])) then
    !
    ! Update the base iterative node
    !
    call cg_biter_write_f(ifile_n,ibase_n,'TimeIterValues',size(soltime),ierr)
    call cgns_error(pname,cgns_grid_file,"cg_biter_write_f",ierr)
    !
    call cg_goto_f(ifile_n,ibase_n,ierr,"BaseIterativeData_t",1,"end")
    call cgns_error(pname,cgns_grid_file,"cg_goto_f",ierr)
    !
    call cg_array_write_f("TimeValues",cgns_kind_real,1, &
                          size(soltime),soltime,ierr)
    call cgns_error(pname,cgns_grid_file,"cg_array_write_f",ierr)
    !
    ! Update the zone iterative node
    !
    call cg_ziter_write_f(ifile_n,ibase_n,izone_n,"ZoneIterativeData",ierr)
    call cgns_error(pname,cgns_grid_file,"cg_ziter_write_f",ierr)
    !
    call cg_goto_f(ifile_n,ibase_n,ierr,"Zone_t",izone_n, &
                                        "ZoneIterativeData_t",1,"end")
    call cgns_error(pname,cgns_grid_file,"cg_goto_f",ierr)
    !
    idata(:) = [32,size(solname)]
    call cg_array_write_f("FlowSolutionPointers",Character,2,idata,solname,ierr)
    call cgns_error(pname,cgns_grid_file,"cg_array_write_f",ierr)
    !
    ! Identify the flow solution within the CGNS file as time dependent
    !
    call cg_simulation_type_write_f(ifile_n,ibase_n,TimeAccurate,ierr)
    call cgns_error(pname,cgns_grid_file,"cg_simulation_type_write_f",ierr)
    !
    ! Close the CGNS FILE
    !
    call cgp_close_f(ifile_n,ierr)
    call io_error(pname,cgns_grid_file,2,ierr)
    !
 !end if
  !
  ! Update the zone counter
  !
  cgnsnodalzone_counter = cgnsnodalzone_counter + 1
  !
  ! Format Statements
  !
  1 format ("#############################################",/, &
            "    FOR PROCESSOR #",i0)
  2 format (a," = ",i0)
  !
end subroutine write_parallel_cgns
!
!###############################################################################
!
subroutine write_solution_parallel_cgns(sol_file,tname)
  !
  !.. Use Statements ..
  use geovar,  only : nr,n_solpts
  use ovar,    only : aref,rhoref,time_ref
  use ovar,    only : cgns_use_queue
  use ovar,    only : interpolate_before_output_variable
  use flowvar, only : usp,dusp
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: sol_file
  character(len=*), intent(in) :: tname
  !
  !.. Local Scalars ..
  integer :: l,m,n,ierr
  integer :: n_total_output_pts
  integer :: n_total_output_cells
  !
  character(len=32) :: varname
  character(len=32) :: sect_name
  character(len=170) :: cgns_sol_file
  character(len=300) :: link_path
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_solution_parallel_cgns"
  character(len=*), parameter :: char_xyz(1:3) = ["X","Y","Z"]
  character(len=*), parameter :: char_sname(1:4) = ["Density  ", &
                                                    "Velocity ", &
                                                    "Pressure ", &
                                                    "Vorticity" ]
  !
  !.. Local Arrays ..
  integer :: ifile_s,ibase_s,izone_s,isol_s
  integer, dimension(1:(nme-nec)*nr) :: idvar_s
  integer, dimension(1:size(char_sname)+2*(nr-1)) :: isvar_s
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable, dimension(:)     :: vartmp
  real(wp), allocatable, dimension(:,:)   :: interp_usp
  real(wp), allocatable, dimension(:,:,:) :: interp_dusp
  !
continue
  !
  cgns_sol_file = trim(adjustl(cgns_dir)) // trim(adjustl(sol_file))
  !
  ! Open the cgns output file if iflag is 1
  !
  ! Initialize the zone counter
  !
  call cgp_open_f(cgns_sol_file,CG_MODE_WRITE,ifile_s,ierr)
  call io_error(pname,cgns_sol_file,1,ierr)
  !
  ! Set the CGNS base (who knows what this is for)
  !
  call cg_base_write_f(ifile_s,base_name,nr,nr,ibase_s,ierr)
  call cgns_error(pname,cgns_sol_file,"cg_base_write_f",ierr)
  !
  ! Write the title for this solution
  !
  call cg_zone_write_f(ifile_s,ibase_s,zone_name,zone_size, &
                       Unstructured,izone_s,ierr)
  call cgns_error(pname,cgns_sol_file,"cg_zone_write_f",ierr)
  !
  ! Create a link to the grid coordinates in the grid file
  !
  call cg_goto_f(ifile_s,ibase_s,ierr,"Zone_t",izone_s,"end")
  call cgns_error(pname,cgns_sol_file,"cg_goto_f",ierr)
  !
  link_path = base_name//"/"//zone_name//"/GridCoordinates"
  !
  call cg_link_write_f("GridCoordinates",trim(cgns_base_fname), &
                       trim(link_path),ierr)
  call cgns_error(pname,cgns_sol_file,"cg_link_write_f(GridCoordinates)",ierr)
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
  call cg_link_write_f(trim(sect_name),trim(cgns_base_fname), &
                       trim(link_path),ierr)
  call cgns_error(pname,cgns_sol_file, &
                  "cg_link_write_f("//trim(sect_name)//")",ierr)
  !
  ! Allocate var and vartmp if they arent already allocated
  !
  if (.not. allocated(vartmp)) then
    allocate ( vartmp(1:n_solpts) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vartmp",1,ierr,error_message)
  end if
  if (.not. allocated(var)) then
    allocate ( var(1:n_my_output_pts) , source=real(0,cgns_flt) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"var",1,ierr,error_message)
  end if
  if (.not. allocated(interp_usp)) then
    allocate ( interp_usp(1:nq,1:n_solpts) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interp_usp",1,ierr,error_message)
  end if
  if (.not. allocated(interp_dusp)) then
    allocate ( interp_dusp(1:nr,1:nq,1:n_solpts) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interp_dusp",1,ierr,error_message)
  end if
  !
  ! Write the flow variables
  !
  call cg_sol_write_f(ifile_s,ibase_s,izone_s,tname,Vertex,isol_s,ierr)
  call cgns_error(pname,cgns_sol_file,"cg_sol_write_f",ierr)
  !
  ! #######################################
  ! #####   WRITING FIELD VARIABLES   #####
  ! #######################################
  !
  if (check_gradients) then
    !
    ! ################################################
    ! #####   WRITING GRADIENT FIELD VARIABLES   #####
    ! ################################################
    !
    ! Create data nodes for the solution variables
    !
    ! Density Gradients
    !
    do l = 1,nr
      varname = "Density_d"//char_xyz(l)
      call cgp_field_write_f(ifile_s,ibase_s,izone_s,isol_s,cgns_kind_real, &
                             trim(varname),idvar_s(l),ierr)
      call cgns_error(pname,cgns_sol_file,"cgp_field_write_f",ierr, &
                      "Write node: Gradient variable",l)
    end do
    !
    ! Velocity Gradients
    !
    do m = 1,nr
      do l = 1,nr
        n = m*nr + l
        varname = "Velocity"//char_xyz(m)//"_d"//char_xyz(l)
        call cgp_field_write_f(ifile_s,ibase_s,izone_s,isol_s,cgns_kind_real, &
                               trim(varname),idvar_s(n),ierr)
        call cgns_error(pname,cgns_sol_file,"cgp_field_write_f",ierr, &
                        "Write node: Gradient variable",n)
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
        vartmp(1:n_solpts) = grad_of_selected_pv(l,m,usp(:,1:n_solpts), &
                                                  dusp(:,:,1:n_solpts))
        !
        var(:) = real( interpolate_for_output(vartmp) , kind(var) )
        !
        ! Reset the queue.
        !
        if (cgns_use_queue) then
          call cgp_queue_set_f(1,ierr)
          call cgns_error(pname,cgns_sol_file,"cgp_queue_set_f",ierr)
        end if
        !
        ! Add the data for current variable to the output queue.
        !
        call cgp_field_write_data_f(ifile_s,ibase_s,izone_s,isol_s,idvar_s(n), &
                                    my_pts_beg,my_pts_end,var,ierr)
        call cgns_error(pname,cgns_sol_file,"cgp_field_write_data_f",ierr, &
                        "Write data: Gradient variable",n)
        !
        ! Finally, flush the queue containing
        ! the current variable to the CGNS file.
        !
        if (cgns_use_queue) then
          call cgp_queue_flush_f(ierr)
          call cgns_error(pname,cgns_sol_file,"cgp_queue_flush_f",ierr)
        end if
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
    ! Create data nodes for the solution variables
    !
    n = 0
    do l = 1,size(char_sname)
      if (trim(char_sname(l)) == "Velocity") then
        do m = 1,nr
          n = n + 1
          varname = "Velocity"//char_xyz(m)
          call cgp_field_write_f(ifile_s,ibase_s,izone_s,isol_s, &
                                 cgns_kind_real,trim(varname),isvar_s(n),ierr)
          call cgns_error(pname,cgns_sol_file,"cgp_field_write_f",ierr, &
                          "Write node: Solution variable",n)
        end do
      else if (trim(char_sname(l)) == "Vorticity") then
        do m = 1,nr
          n = n + 1
          varname = "Vorticity"//char_xyz(m)
          call cgp_field_write_f(ifile_s,ibase_s,izone_s,isol_s, &
                                 cgns_kind_real,trim(varname),isvar_s(n),ierr)
          call cgns_error(pname,cgns_sol_file,"cgp_field_write_f",ierr, &
                          "Write node: Solution variable",n)
        end do
      else
        n = n + 1
        call cgp_field_write_f(ifile_s,ibase_s,izone_s,isol_s,cgns_kind_real, &
                               trim(char_sname(l)),isvar_s(n),ierr)
        call cgns_error(pname,cgns_sol_file,"cgp_field_write_f",ierr, &
                        "Write node: Solution variable",n)
      end if
    end do
    !
    ! Interpolate the conservative variables for output
    !
    if (interpolate_before_output_variable) then
      interp_usp(1:nq,1:n_solpts) = &
            interpolate_for_output( interp_usp(1:nq,1:n_solpts) )
      do l = 1,nq
        interp_dusp(1:nr,l,1:n_solpts) = &
                 interpolate_for_output( interp_dusp(1:nr,l,1:n_solpts) )
      end do
    else
      interp_usp(1:nq,1:n_solpts) = usp(1:nq,1:n_solpts)
      interp_dusp(1:nr,1:nq,1:n_solpts) = dusp(1:nr,1:nq,1:n_solpts)
    end if
    !
    ! Now create the output data and write it to the CGNS file in parallel
    !
    n = 0
    do l = 1,size(char_sname)
      !
      select case (trim(char_sname(l)))
        case ("Density") ! Density
          n = n + 1
          vartmp(1:n_solpts) = interp_usp(nec,1:n_solpts)
          vartmp = vartmp * rhoref
          call cgns_write_field_parallel
        case ("Velocity") ! Velocity
          do m = nmb,nme
            n = n + 1
            vartmp(1:n_solpts) = interp_usp(m,1:n_solpts) / &
                                 interp_usp(nec,1:n_solpts)
            vartmp = vartmp * aref
            call cgns_write_field_parallel
          end do
        case ("Pressure") ! Pressure
          n = n + 1
          vartmp = pressure_cv( interp_usp )
          vartmp = vartmp * rhoref * aref * aref
          call cgns_write_field_parallel
        case ("Mach") ! Mach number
          n = n + 1
          vartmp = mach_number_cv( interp_usp )
          call cgns_write_field_parallel
        case ("Entropy") ! Entropy
          n = n + 1
          vartmp = entropy_cv( interp_usp )
          call cgns_write_field_parallel
        case ("Vorticity") ! Vorticity
         !vartmp = entropy_cv( interp_usp )
         !vartmp = vorticity_cv( interp_usp )
          do m = 1,nr
            n = n + 1
            vartmp = vorticity_dusp( interp_usp , interp_dusp , m )
            vartmp = vartmp / time_ref
            call cgns_write_field_parallel
          end do
        case default
          n = n + 1
          vartmp(1:n_solpts) = interp_usp(1,1:n_solpts)
          call cgns_write_field_parallel
      end select
      !
    end do
    !
  end if
  !
  ! ################################################
  ! #####   FINISHED WRITING FIELD VARIABLES   #####
  ! ################################################
  !
  if (cgns_use_queue) then
    call cgp_queue_set_f(0,ierr)
    call cgns_error(pname,cgns_sol_file,"cgp_queue_set_f",ierr)
  end if
  !
  ! Close the CGNS FILE
  !
  call cgp_close_f(ifile_s,ierr)
  call io_error(pname,cgns_grid_file,2,ierr)
  !
  ! Deallocate the arrays var and vartmp
  !
  if (allocated(var)) then
    deallocate ( var , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"var",2,ierr,error_message)
  end if
  if (allocated(vartmp)) then
    deallocate ( vartmp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vartmp",2,ierr,error_message)
  end if
  if (allocated(interp_usp)) then
    deallocate ( interp_usp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interp_usp",2,ierr,error_message)
  end if
  if (allocated(interp_dusp)) then
    deallocate ( interp_dusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"interp_dusp",2,ierr,error_message)
  end if
  !
  ! Format Statements
  !
  1 format ("#############################################",/, &
            "    FOR PROCESSOR #",i0)
  2 format (a," = ",i0)
  !
  !-----------------------------------------------------------------------------
  !
  contains
  !
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
    ! Reset the queue.
    !
    if (cgns_use_queue) then
      call cgp_queue_set_f(1,ierr)
      call cgns_error(pname,cgns_sol_file,"cgp_queue_set_f",ierr)
    end if
    !
    ! Add the data for current variable to the output queue.
    !
    call cgp_field_write_data_f(ifile_s,ibase_s,izone_s,isol_s,isvar_s(n), &
                                my_pts_beg,my_pts_end,var,ierr)
    call cgns_error(pname,cgns_sol_file,"cgp_field_write_data_f",ierr, &
                    "Write data: Solution variable",n)
    !
    ! Finally, flush the queue containing
    ! the current variable to the CGNS file.
    !
    if (cgns_use_queue) then
      call cgp_queue_flush_f(ierr)
      call cgns_error(pname,cgns_sol_file,"cgp_queue_flush_f",ierr)
    end if
    !
  end subroutine cgns_write_field_parallel
  !
  !-----------------------------------------------------------------------------
  !
end subroutine write_solution_parallel_cgns
!
!###############################################################################
!
subroutine get_initial_data_for_cgns(my_cells_beg,my_cells_end,ipelem)
  !
  !.. Use Statements ..
  use geovar,       only : ncell,n_global_cell,cell
  use geovar,       only : pst_geom
  use geovar,       only : global_cell_geom
  use geovar,       only : global_cell_order
  use parallel_mod, only : cell_map
  use order_mod,    only : geom_solpts
  !
  !.. Formal Arguments ..
  integer,                              intent(inout) :: my_cells_beg
  integer,                              intent(inout) :: my_cells_end
  integer, allocatable, dimension(:,:), intent(inout) :: ipelem
  !
  !.. Local Scalars ..
  integer :: n,nc,np,ierr,p1,icpu
  integer :: gmin,gmax,this_npts,this_subcells
  integer :: this_cell,this_geom,this_order
  integer :: n_total_output_cells
  integer :: n_total_output_pts
  integer :: n_my_output_cells
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_initial_data_for_cgns"
  !
continue
  !
  ! Initialize n_total_output_cells and n_total_output_pts which give:
  !
  ! n_total_output_cells: total number of output cells for the global grid
  ! n_total_output_pts: total number of output points for the global grid
  !
  n_total_output_pts = 0
  n_total_output_cells = 0
  !
  do icpu = 1,ncpu
    !
    if (icpu-1 == mypnum) then
      my_pts_beg = n_total_output_pts + 1
      my_cells_beg = n_total_output_cells + 1
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
      this_npts = geom_solpts(this_geom,this_order)
      !
      ! Get the number of subcells for this cell geometry at this order.
      ! The number of subcells is basically:
      !     this_order**(spatial dimension of this_geom) or
      ! NOTE: Make sure that the number of subcells is at least 1,
      !       particularly for the case of this_order=0
      !
      this_subcells = max( 1 , this_order**geom_dimen(this_geom) )
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
      n_total_output_cells = n_total_output_cells + this_subcells
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
  zone_size = [n_total_output_pts,n_total_output_cells,0]
  !
  ! n_my_output_pts: number of output points for this processor
  !
  n_my_output_pts = my_pts_end - my_pts_beg + 1
  !
  ! n_my_output_cells: number of output cell for this processor
  !
  n_my_output_cells = my_cells_end - my_cells_beg + 1
  !
  ! Some quick stuff needed before creating the connectivity array
  !
  gmin = minval( global_cell_geom(1:n_global_cell) )
  gmax = maxval( global_cell_geom(1:n_global_cell) )
  !
  if (any(Geom_Unknown == [gmin,gmax])) then
    write (error_message,1)
    call stop_recon(stop_mpi,pname,error_message)
  else if (geom_dimen(gmin) /= geom_dimen(gmax)) then
    write (error_message,2)
    call stop_recon(stop_mpi,pname,error_message)
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
    call alloc_error(pname,"ipelem",2,ierr,error_message)
  end if
  !
  allocate ( ipelem(1:np,1:n_my_output_cells) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ipelem",1,ierr,error_message)
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
    this_npts = geom_solpts(this_geom,this_order)
    !
    ! Get the number of subcells for this cell geometry at this order.
    ! The number of subcells is basically:
    !     this_order**(spatial dimension of this_geom) or
    ! NOTE: Make sure that the number of subcells is at least 1,
    !       particularly for the case of this_order=0
    !
    this_subcells = max( 1 , this_order**geom_dimen(this_geom) )
    !
    ! Loop through the number of subcells and add the subconnectivity
    ! to the first solution point of the current grid cell to get the
    ! connectivity for the current subcell.
    !
    do n = 1,this_subcells
      !
      nc = nc + 1
      !
      ipelem(:,nc) = p1 + pst_geom(this_geom,this_order)%connectivity(:,n)
      !
    end do
    !
    p1 = p1 + this_npts
    !
  end do
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
  use ovar, only : cgns_use_queue
  use ovar, only : write_TaylorGreen_8s_solution
  !
  !.. Local Scalars ..
  integer :: my_beg,my_end,my_pts,ierr
  integer :: ifile_v,ibase_v,izone_v
  integer :: isect_v,isol_v,isvar_v
  !
  character(len=170) :: cgns_vort_file
  !
  !.. Local Arrays ..
  integer,  dimension(1:2) :: ixyz_v
  integer,  dimension(1:3) :: my_vort_pts
  integer,  dimension(1:3) :: my_vort_cells
  integer,  dimension(1:3) :: vort_zone_size
  !
  !.. Local Allocatable Arrays ..
  integer,        allocatable, dimension(:)   :: my_vort_faces
  integer,        allocatable, dimension(:,:) :: ipelem
  real(cgns_flt), allocatable, dimension(:)   :: vary
  real(cgns_flt), allocatable, dimension(:)   :: varz
  !
  !.. Local Parameters ..
  character(len=*), parameter :: bname = "TGV_vort_8s_base"
  character(len=*), parameter :: zname = "TGV_vort_8s_zone"
  character(len=*), parameter :: fname = "TGV_vorticity_8s.cgns"
  character(len=*), parameter :: pname = "TaylorGreen_8s_vorticity_solution"
  !
continue
  !
  cgns_vort_file = trim(adjustl(cgns_dir)) // fname
  !
  ! Get/collect all the data necessary to output the 8s vorticity solution
  !
  call get_data_for_8s_vorticity_solution(my_vort_pts,my_vort_cells, &
                                          vort_zone_size,my_vort_faces,ipelem)
  !
  ! Open the CGNS file
  !
  call cgp_open_f(cgns_vort_file,CG_MODE_WRITE,ifile_v,ierr)
  call io_error(pname,cgns_vort_file,1,ierr)
  !
  ! Set the CGNS base (who knows what this is for)
  ! NOTE: We are creating this base 2-dimensional with surface (2D) cells
  !       since this file will contain a 2D solution along a boundary face
  !
  call cg_base_write_f(ifile_v,bname,2,2,ibase_v,ierr)
  call cgns_error(pname,cgns_vort_file,"cg_base_write_f",ierr)
  !
  ! Write the title for this solution
  !
  call cg_zone_write_f(ifile_v,ibase_v,zname,vort_zone_size, &
                       Unstructured,izone_v,ierr)
  call cgns_error(pname,cgns_vort_file,"cg_zone_write_f",ierr)
  !
  ! Create data nodes for the coordinates
  !
  call cgp_coord_write_f(ifile_v,ibase_v,izone_v,cgns_kind_real, &
                         "CoordinateY",ixyz_v(1),ierr)
  call cgns_error(pname,cgns_vort_file,"cgp_coord_write_f",ierr, &
                  "Write node: Coordinate = 1 (Y)")
  !
  call cgp_coord_write_f(ifile_v,ibase_v,izone_v,cgns_kind_real, &
                         "CoordinateZ",ixyz_v(2),ierr)
  call cgns_error(pname,cgns_vort_file,"cgp_coord_write_f",ierr, &
                  "Write node: Coordinate = 2 (Z)")
  !
  ! Now create and write the coordinate data in parallel
  !
  if (.not. cgns_use_queue) then
    call cgp_queue_set_f(0,ierr)
    call cgns_error(pname,cgns_vort_file,"cgp_queue_set_f",ierr)
  end if
  !
  ! For each boundary face within the grid region of interest, interpolate the
  ! coordinates of the solution points for the host cell to the boundary face
  ! and this interpolate those face points to the output face points.
  !
  my_pts = my_vort_pts(1)
  my_beg = my_vort_pts(2)
  my_end = my_vort_pts(3)
  !
  allocate ( vary(1:my_pts) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"vary",1,ierr,error_message)
  vary = real( zero , kind=kind(vary) )
  !
  allocate ( varz(1:my_pts) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"varz",1,ierr,error_message)
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
  ierr = 0
  if (my_end-my_beg+1 > 0) then
    call cgp_coord_write_data_f(ifile_v,ibase_v,izone_v,ixyz_v(1), &
                                my_beg,my_end,vary,ierr)
  end if
  call cgns_error(pname,cgns_vort_file,"cgp_coord_write_data_f",ierr, &
                  "Write data: Coordinate = 1 (Y)")
  !
  ! Add the data for the z-coordinate to the output queue
  ! NOTE: It is an error for a processor to try an write data if it doesnt
  !       have any data to write, so make sure that only those processors
  !       that have data to write call the CGNS write subroutine.
  !
  ierr = 0
  if (my_end-my_beg+1 > 0) then
    call cgp_coord_write_data_f(ifile_v,ibase_v,izone_v,ixyz_v(2), &
                                my_beg,my_end,varz,ierr)
  end if
  call cgns_error(pname,cgns_vort_file,"cgp_coord_write_data_f",ierr, &
                  "Write data: Coordinate = 2 (Z)")
  !
  ! Finally, flush the queue containing the
  ! coordinate variables to the CGNS file
  !
  if (cgns_use_queue) then
    call cgp_queue_flush_f(ierr)
    call cgns_error(pname,cgns_vort_file,"cgp_queue_flush_f",ierr)
  end if
  !
  ! Create the data node for the grid connectivity
  !
  call cgp_section_write_f(ifile_v,ibase_v,izone_v,"Quadrilaterals", &
                           QUAD_4,1,vort_zone_size(2),0,isect_v,ierr)
  call cgns_error(pname,cgns_vort_file,"cgp_section_write_f",ierr)
  !
  ! Write the element connectivity in parallel
  !
  ! Reset the queue.
  !
  if (cgns_use_queue) then
    call cgp_queue_set_f(1,ierr)
    call cgns_error(pname,cgns_vort_file,"cgp_queue_set_f",ierr)
  end if
  !
  my_beg = my_vort_cells(2)
  my_end = my_vort_cells(3)
  !
  ! Add the element connectivity data to the output queue.
  ! NOTE: It is an error for a processor to try an write data if it doesnt
  !       have any data to write, so make sure that only those processors
  !       that have data to write call the CGNS write subroutine.
  !
  ierr = 0
  if (my_end-my_beg+1 > 0) then
    call cgp_elements_write_data_f(ifile_v,ibase_v,izone_v,isect_v, &
                                   my_beg,my_end,ipelem,ierr)
  end if
  call cgns_error(pname,cgns_vort_file,"cgp_elements_write_data_f",ierr)
  !
  ! Finally, flush the queue containing the
  ! element connectivity to the CGNS file.
  !
  if (cgns_use_queue) then
    call cgp_queue_flush_f(ierr)
    call cgns_error(pname,cgns_vort_file,"cgp_queue_flush_f",ierr)
  end if
  !
  ! Deallocate the arrays varz, and ipelem
  ! NOTE: Keep vary allocated and use it to store the vorticity variable
  !
  if (allocated(varz)) then
    deallocate ( varz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"varz",2,ierr,error_message)
  end if
  if (allocated(ipelem)) then
    deallocate ( ipelem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ipelem",2,ierr,error_message)
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
                      Vertex,isol_v,ierr)
  call cgns_error(pname,cgns_vort_file,"cg_sol_write_f",ierr)
  !
  ! Create a data node for the vorticity solution
  !
  call cgp_field_write_f(ifile_v,ibase_v,izone_v,isol_v,cgns_kind_real, &
                         "VorticityMagnitude",isvar_v,ierr)
  call cgns_error(pname,cgns_vort_file,"cgp_field_write_f",ierr, &
                  "Write node: Solution variable = 1 (VorticityMagnitude)")
  !
  ! Compute the vorticity magnitude at each face point
  !
  call get_vort_face_solution(my_vort_faces,vary)
  !
  ! Reset the queue.
  !
  if (cgns_use_queue) then
    call cgp_queue_set_f(1,ierr)
    call cgns_error(pname,cgns_vort_file,"cgp_queue_set_f",ierr)
  end if
  !
  ! Add the data for current variable to the output queue.
  ! NOTE: It is an error for a processor to try an write data if it doesnt
  !       have any data to write, so make sure that only those processors
  !       that have data to write call the CGNS write subroutine.
  !
  ierr = 0
  if (my_end-my_beg+1 > 0) then
    call cgp_field_write_data_f(ifile_v,ibase_v,izone_v,isol_v,isvar_v, &
                                my_beg,my_end,vary,ierr)
  end if
  call cgns_error(pname,cgns_vort_file,"cgp_field_write_data_f",ierr, &
                  "Write data: Solution variable = 1 (VorticityMagnitude)")
  !
  ! Finally, flush the queue containing
  ! the current variable to the CGNS file.
  !
  if (cgns_use_queue) then
    call cgp_queue_flush_f(ierr)
    call cgns_error(pname,cgns_vort_file,"cgp_queue_flush_f",ierr)
  end if
  !
  ! Deallocate vary now that we are finished with it
  !
  if (allocated(vary)) then
    deallocate ( vary , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vary",2,ierr,error_message)
  end if
  !
  ! Close the CGNS file
  !
  call cgp_close_f(ifile_v,ierr)
  call io_error(pname,cgns_vort_file,2,ierr)
  !
  write_TaylorGreen_8s_solution = fals
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
  use ovar, only : Lref,TaylorGreen_DOF
  use geovar, only : nfbnd,face,nodes_on_face
  use geovar, only : xyz_nodes,pst_geom
  use order_mod, only : geom_solpts,n_order
  !
  !.. Formal Arguments ..
  integer,              dimension(1:3), intent(inout) :: my_vort_pts
  integer,              dimension(1:3), intent(inout) :: my_vort_cells
  integer,              dimension(1:3), intent(inout) :: vort_zone_size
  integer, allocatable, dimension(:),   intent(inout) :: my_vort_faces
  integer, allocatable, dimension(:,:), intent(inout) :: ipelem
  !
  !.. Local Scalars ..
  integer :: nc,nf,p1,n,ierr,ip
  integer :: this_face,face_geom,face_order
  integer :: face_npts,face_subcells
  integer :: my_vort_cells_beg,my_vort_pts_beg
  integer :: my_vort_cells_end,my_vort_pts_end
  integer :: n_my_vort_cells,n_my_vort_pts
  integer :: n_total_vort_cells,n_total_vort_pts
  !
  real(wp) :: LPI,half_LPI,dx_eps,cell_dx,x,y,z
  !
  !.. Local Allocatable Arrays ..
  logical, allocatable, dimension(:)   :: maskf
  integer, allocatable, dimension(:,:) :: cpu_vort_data
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_data_for_8s_vorticity"
  !
continue
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
  cell_dx = two * LPI / anint( real(TaylorGreen_DOF,kind=wp) / &
                               real(n_order+1,kind=wp) , kind=wp )
  dx_eps = cell_dx * eps3 / ten
  !
  ! Initialize the variables giving the total number of vorticity points
  ! and subcells as well as the variables giving the same information for
  ! the respective information local to the current processor
  !
  n_total_vort_pts = 0
  n_total_vort_cells = 0
  !
  n_my_vort_pts = 0
  n_my_vort_cells = 0
  !
  allocate ( maskf(1:nfbnd) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"maskf",1,ierr,error_message)
  !
  ! Find the boundary faces that lie on the plane x = -PI*Lref within the region
  ! y = [0,0.5*PI*Lref] and z = [0.5*PI*Lref,PI*Lref]
  !
  maskf_loop: do nf = 1,nfbnd
    !
    bface_node_loop: do n = 1,size(nodes_on_face,dim=1)
      !
      face_geom = face(nf)%geom
      face_order = face(nf)%order ! order of boundary face
      !
      ! NOTE: We dont care here if any of the nodes stored in nodes_on_face
      !       for this face are duplicates. The if statements will just produce
      !       a false value for the second time.  We only care if the value
      !       ip that is returned is zero because this will cause the xyz_nodes
      !       array to go out-of-bounds so cycle the inner loop if this
      !       situation arrises.
      !       
      ip = nodes_on_face(n,nf)
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
            face_npts = geom_solpts(face_geom,face_order)
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
            n_my_vort_cells = n_my_vort_cells + face_subcells
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
    call alloc_error(pname,"my_vort_faces",2,ierr,error_message)
  end if
  !
  allocate ( my_vort_faces(1:count(maskf)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"my_vort_faces",1,ierr,error_message)
  !
  my_vort_faces = pack( intseq(1,nfbnd) , mask=maskf )
  !
  ! Deallocate the maskf array since it is no longer needed
  !
  deallocate ( maskf , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"maskf",2,ierr,error_message)
  !
  ! We need to collect this information for all processors
  !
  allocate ( cpu_vort_data(1:2,1:ncpu) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cpu_vort_data",1,ierr,error_message)
  !
  cpu_vort_data(1,mypnum+1) = n_my_vort_pts
  cpu_vort_data(2,mypnum+1) = n_my_vort_cells
  !
  if (ncpu > 1) then
    call mpi_allgather(MPI_IN_PLACE, size(cpu_vort_data,dim=1),mpi_inttyp, &
                       cpu_vort_data,size(cpu_vort_data,dim=1),mpi_inttyp, &
                       mpi_comm_world,mpierr)
  end if
  !
  ! Sum up all the vorticity points and subcells from all processors
  !
  n_total_vort_pts   = sum( cpu_vort_data(1,:) )
  n_total_vort_cells = sum( cpu_vort_data(2,:) )
  !
  ! Save the total number of vorticity points and subcells to vort_zone_size
  !
  vort_zone_size = [n_total_vort_pts,n_total_vort_cells,0]
  !
  ! Get the beginning index for the vorticity points and subcells local to this
  ! processor by summing up all the vorticity points and subcells from all
  ! processors with a smaller rank than the this processor.
  ! NOTE: The second index of cpu_vort_data goes from 1:ncpu whereas the value
  !       given by mypnum is from 0:ncpu-1. Therefore, summing from 1:mypnum
  !       gives the range of all processors with a smaller rank.
  !
  my_vort_pts_beg   = sum( cpu_vort_data(1,1:mypnum) ) + 1
  my_vort_cells_beg = sum( cpu_vort_data(2,1:mypnum) ) + 1
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
  call alloc_error(pname,"cpu_vort_data",2,ierr,error_message)
  !
  ! Create the connectivity array
  ! NOTE: The Taylor-Green vortex case is 3D so we know the boundary faces
  !       will always be 2D and they will require 4 nodes to define each
  !       subcell within the CGNS file
  !
  if (allocated(ipelem)) then
    deallocate ( ipelem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ipelem",2,ierr,error_message)
  end if
  !
  allocate ( ipelem(1:4,1:n_my_vort_cells) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ipelem",1,ierr,error_message)
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
    face_npts = geom_solpts(face_geom,face_order)
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
      ipelem(:,nc) = p1 + pst_geom(face_geom,face_order)%connectivity(:,n)
      !
    end do
    !
    p1 = p1 + face_npts
    !
  end do
  !
end subroutine get_data_for_8s_vorticity_solution
!
!###############################################################################
!
subroutine get_vort_face_coordinates(my_vort_faces,y,z)
  !
  !.. Use Statements ..
  use order_mod, only : maxSP,maxFP
  use order_mod, only : geom_solpts
  use geovar, only : face,cell,xyz
  use interpolation_mod, only : interp
  use interpolation_mod, only : outerp
  !
  !.. Formal Arguments ..
  integer,        dimension(:),    intent(in) :: my_vort_faces
  real(cgns_flt), dimension(:), intent(inout) :: y
  real(cgns_flt), dimension(:), intent(inout) :: z
  !
  !.. Local Scalars ..
  integer :: nf,np,n1,n2,fp,f1,f2,k,kl
  integer :: this_face,face_geom,face_order
  integer :: host_cell,host_geom,host_order
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxFP) :: face_y
  real(wp), dimension(1:maxFP) :: face_z
  real(wp), dimension(1:maxSP,1:2) :: txyz
#ifndef USE_MATRIX_PROCEDURES
  real(wp), dimension(1:maxSP) :: lmat
#endif
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_vort_face_coordinates"
  !
continue
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
    fp = geom_solpts(face_geom,face_order) ! number of sol-pts on boundary face
    f1 = f2 + 1                  ! beginning index for sol-pts on boundary face
    f2 = f2 + fp                 ! ending    index for sol-pts on boundary face
    !
    ! Create transpose of the xyz(2:3,n1:n2) to reduce cache misses
    !
    txyz(1:np,1:2) = transpose( xyz(2:3,n1:n2) )
    !
    ! Compute the profile quantities at each data point for this face
    !
    do k = 1,size(face(this_face)%left%flxpts)
      !
      ! Index of flx-pt local to the current host cell
      !
      kl = face(this_face)%left%flxpts(k)
      !
      ! First interpolate to the face flux points
      !
#ifdef USE_MATRIX_PROCEDURES
      face_y(k) = &
       interp(host_geom,host_order)%toFace(face_order)%dot( kl , txyz(1:np,1) )
      face_z(k) = &
       interp(host_geom,host_order)%toFace(face_order)%dot( kl , txyz(1:np,2) )
#else
      lmat(1:np) = interp(host_geom,host_order)%toFace(face_order)%mat(1:np,kl)
      !
      face_y(k) = dot( lmat(1:np) , txyz(1:np,1) )
      face_z(k) = dot( lmat(1:np) , txyz(1:np,2) )
#endif
      !
    end do
    !
    ! Finally, interpolate the face flux points
    ! to the output face solution points
    !
    if (associated(outerp(face_geom,face_order)%output)) then
      !
#ifdef USE_MATRIX_PROCEDURES
      face_y(1:fp) = outerp(face_geom,face_order)%output%mv_mult( face_y(1:fp) )
      face_z(1:fp) = outerp(face_geom,face_order)%output%mv_mult( face_z(1:fp) )
#else
      face_y(1:fp) = matmul( outerp(face_geom,face_order)%output%mat, &
                             face_y(1:fp) )
      face_z(1:fp) = matmul( outerp(face_geom,face_order)%output%mat, &
                             face_z(1:fp) )
#endif
      !
    end if
    !
    y(f1:f2) = real( face_y(1:fp) , kind=kind(y) )
    z(f1:f2) = real( face_z(1:fp) , kind=kind(z) )
    !
  end do
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
  use geovar, only : nr,face,cell
  use ovar, only : time_ref
  use interpolation_mod, only : interp
  use interpolation_mod, only : outerp
  use flowvar, only : usp,dusp
  !
  !.. Formal Arguments ..
  integer,        dimension(:),    intent(in) :: my_vort_faces
  real(cgns_flt), dimension(:), intent(inout) :: vorticity
  !
  !.. Local Scalars ..
  integer :: k,kl,l,m,nf,np,n1,n2,fp,f1,f2
  integer :: this_face,face_geom,face_order
  integer :: host_cell,host_geom,host_order
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: vort_vec
  real(wp), dimension(1:nr,1:nq) :: dpvdx
  real(wp), dimension(1:nq,1:maxFP) :: cvfp
  real(wp), dimension(1:maxFP,1:nq) :: cv_tmp
  real(wp), dimension(1:maxSP,1:nq) :: tusp
  real(wp), dimension(1:nr,1:nq,1:maxFP) :: dcvdx
  real(wp), dimension(1:maxFP,1:nr,1:nq) :: dcvdx_tmp
  real(wp), dimension(1:maxSP,1:nr,1:nq) :: tdusp
#ifndef USE_MATRIX_PROCEDURES
  real(wp), dimension(1:maxSP) :: lmat
#endif
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_vort_face_solution"
  !
continue
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
    fp = geom_solpts(face_geom,face_order) ! number of sol-pts on boundary face
    f1 = f2 + 1                  ! beginning index for sol-pts on boundary face
    f2 = f2 + fp                 ! ending    index for sol-pts on boundary face
    !
    ! Create transpose of the xyz(2:3,n1:n2) to reduce cache misses
    !
    tusp(1:np,1:nq) = transpose( usp(1:nq,n1:n2) )
    do k = 1,np
      tdusp(k,1:nr,1:nq) = dusp(1:nr,1:nq,n1-1+k)
    end do
    !
    ! Compute the profile quantities at each data point for this face
    !
    do k = 1,size(face(this_face)%left%flxpts)
      !
      ! Index of flx-pt local to the current host cell
      !
      kl = face(this_face)%left%flxpts(k)
      !
      ! First interpolate to the face flux points
      !
#ifdef USE_MATRIX_PROCEDURES
      do m = 1,nq
        cv_tmp(k,m) = &
    interp(host_geom,host_order)%toFace(face_order)%dot( kl , tusp(1:np,m) )
        do l = 1,nr
          dcvdx_tmp(k,l,m) = &
    interp(host_geom,host_order)%toFace(face_order)%dot( kl , tdusp(1:np,l,m) )
        end do
      end do
#else
      lmat(1:np) = interp(host_geom,host_order)%toFace(face_order)%mat(1:np,kl)
      !
      do m = 1,nq
        cv_tmp(k,m) = dot( lmat(1:np) , tusp(1:np,m) )
        do l = 1,nr
          dcvdx_tmp(k,l,m) =  dot( lmat(1:np) , tdusp(1:np,l,m) )
        end do
      end do
#endif
      !
    end do
    !
    ! Now, interpolate the vorticity at the face flux points
    ! to the output face solution points
    !
    if (associated(outerp(face_geom,face_order)%output)) then
      !
#ifdef USE_MATRIX_PROCEDURES
      do m = 1,nq
        cvfp(m,1:fp) = &
            outerp(face_geom,face_order)%output%mv_mult( cv_tmp(1:fp,m) )
        do l = 1,nr
          dcvdx(l,m,1:fp) = &
            outerp(face_geom,face_order)%output%mv_mult( dcvdx_tmp(1:fp,l,m) )
        end do
      end do
#else
      do m = 1,nq
        cvfp(m,1:fp) = &
            matmul( outerp(face_geom,face_order)%output%mat , cv_tmp(1:fp,m) )
        do l = 1,nr
          dcvdx(l,m,1:fp) = &
            matmul( outerp(face_geom,face_order)%output%mat , &
                    dcvdx_tmp(1:fp,l,m) )
        end do
      end do
#endif
      !
    else
      !
      cvfp(1:nq,1:fp) = transpose( cv_tmp(1:fp,1:nq) )
      !
      do k = 1,fp
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
    do k = 1,fp
      !
      ! Compute the derivatives of the primitive variables
      !
      dpvdx(1:nr,1:nq) = grad_cv_to_grad_pv_sp( cvfp(1:nq,k) , &
                                                dcvdx(1:nr,1:nq,k) )
      !
      ! Compute the vorticity vector at this face point
      !
      vort_vec(1) = dpvdx(2,nmz) - dpvdx(3,nmy)
      vort_vec(2) = dpvdx(3,nmx) - dpvdx(1,nmz)
      vort_vec(3) = dpvdx(1,nmy) - dpvdx(2,nmx)
      !
      ! Compute the magnitude of the vorticity vector
      ! NOTE: Dont forget to adjust the nondimensionalization of the
      !       vorticity by dividing through by time_ref to get vorticity*
      !        
      vorticity(f1-1+k) = real( rnorm2(vort_vec) / time_ref , &
                                kind=kind(vorticity) )
      !
    end do
    !
  end do
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
  real(wp), dimension(1:size(dcvdx,dim=3)) :: return_value
  !
  !.. Local Scalars ..
  integer :: nr,np,n
  real(wp) :: cv_inv,rho_inv,rho_inv2,q2
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(dcvdx,dim=2)) :: coef
  !
continue
  !
  nr = size(dcvdx,dim=1)
  np = size(dcvdx,dim=3)
  !
  if (npv == nec) then
    !
    do n = 1,np
      return_value(n) = dcvdx(ndr,npv,n)
    end do
    !
  else if (npv == nee) then
    !
    cv_inv = (gam - one) / rgasref_nd
    !
    do n = 1,np
      !
      rho_inv = one/uu(nec,n)
      rho_inv2 = rho_inv*rho_inv
      !
      q2 = dot( uu(nmb:nme,n) , uu(nmb:nme,n) )
      !
      coef(nec)     = -rho_inv2 * (uu(nee,n) - rho_inv*q2)
      coef(nmb:nme) = -rho_inv2 * uu(nmb:nme,n)
      coef(nee)     =  rho_inv
      !
      return_value(n) = cv_inv * dot( coef , dcvdx(ndr,1:nq,n) )
      !
    end do
    !
  else
    !
    do n = 1,np
      !
      rho_inv = one/uu(nec,n)
      rho_inv2 = rho_inv*rho_inv
      !
      return_value(n) = rho_inv  * dcvdx(ndr,npv,n) - &
                        rho_inv2 * dcvdx(ndr,nec,n) * uu(npv,n)
      !
    end do
    !
  end if
  !
end function grad_of_selected_pv
!
!###############################################################################
!
pure function interpolate_for_output_1d(array) result(return_value)
  !
  !... calculates d(utd)/dt
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts,n_global_solpts
  use geovar,  only : ncell,n_global_cell,cell
  use geovar,  only : global_pinc_ptr
  use geovar,  only : global_cell_geom
  use geovar,  only : global_cell_order
  !
  use interpolation_mod, only : outerp
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: array
  !
  !.. Function Return Value ..
  real(wp), dimension(1:size(array)) :: return_value
  !
  !.. Local Scalars ..
  integer :: nc,n1,n2,np,n
  integer :: this_geom,this_order
  !
continue 
  !
  ! Interpolate the data depending on if we are working on
  ! a local or global data array.
  !
  if (size(array) == n_solpts) then
    !
    do nc = 1,ncell
      !
      this_geom = cell(nc)%geom
      this_order = cell(nc)%order
      !
      n1 = cell(nc)%beg_sp
      n2 = cell(nc)%end_sp
      np = n2-n1+1
      !
      if (associated(outerp(this_geom,this_order)%output)) then
#ifdef USE_MATRIX_PROCEDURES
        return_value(n1:n2) = &
          outerp(this_geom,this_order)%output%mv_mult( array(n1:n2) )
#else
        return_value(n1:n2) = &
          matmul( outerp(this_geom,this_order)%output%mat , array(n1:n2) )
#endif
      else
        return_value(n1:n2) = array(n1:n2)
      end if
      !
    end do
    !
  else if (size(array) == n_global_solpts) then
    !
    do nc = 1,n_global_cell
      !
      this_geom = global_cell_geom(nc)
      this_order = global_cell_order(nc)
      !
      n1 = global_pinc_ptr(nc)+1
      n2 = global_pinc_ptr(nc+1)
      np = n2-n1+1
      !
      if (associated(outerp(this_geom,this_order)%output)) then
#ifdef USE_MATRIX_PROCEDURES
        return_value(n1:n2) = &
          outerp(this_geom,this_order)%output%mv_mult( array(n1:n2) )
#else
        return_value(n1:n2) = &
          matmul( outerp(this_geom,this_order)%output%mat , array(n1:n2) )
#endif
      else
        return_value(n1:n2) = array(n1:n2)
      end if
      !
    end do
    !
  else
    !
    do n = 1,size(array)
      return_value(n) = array(n)
    end do
    !
  end if
  ! 
end function interpolate_for_output_1d
!
!###############################################################################
!
pure function interpolate_for_output_2d(array) result(return_value)
  !
  !... calculates d(utd)/dt
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts,n_global_solpts
  use geovar,  only : ncell,n_global_cell,cell
  use geovar,  only : global_pinc_ptr
  use geovar,  only : global_cell_geom
  use geovar,  only : global_cell_order
  !
  use interpolation_mod, only : outerp
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: array
  !
  !.. Function Return Value ..
  real(wp) :: return_value( 1:size(array,dim=1) , 1:size(array,dim=2) )
  !
  !.. Local Scalars ..
  integer :: nc,n1,n2,np,n,npts,nvar,m
  integer :: this_geom,this_order
  !
continue 
  !
  nvar = size(array,dim=1)
  npts = size(array,dim=2)
  !
  ! Interpolate the data depending on if we are working on
  ! a local or global data array.
  !
  if (npts == n_solpts) then
    !
    do nc = 1,ncell
      !
      this_geom = cell(nc)%geom
      this_order = cell(nc)%order
      !
      n1 = cell(nc)%beg_sp
      n2 = cell(nc)%end_sp
      np = n2-n1+1
      !
      if (associated(outerp(this_geom,this_order)%output)) then
        do m = 1,nvar
#ifdef USE_MATRIX_PROCEDURES
          return_value(m,n1:n2) = &
            outerp(this_geom,this_order)%output%mv_mult( array(m,n1:n2) )
#else
          return_value(m,n1:n2) = &
            matmul( outerp(this_geom,this_order)%output%mat , array(m,n1:n2) )
#endif
        end do
      else
        do n = n1,n2
          do m = 1,nvar
            return_value(m,n) = array(m,n)
          end do
        end do
      end if
      !
    end do
    !
  else if (npts == n_global_solpts) then
    !
    do nc = 1,n_global_cell
      !
      this_geom = global_cell_geom(nc)
      this_order = global_cell_order(nc)
      !
      n1 = global_pinc_ptr(nc)+1
      n2 = global_pinc_ptr(nc+1)
      np = n2-n1+1
      !
      if (associated(outerp(this_geom,this_order)%output)) then
        do m = 1,nvar
#ifdef USE_MATRIX_PROCEDURES
          return_value(m,n1:n2) = &
            outerp(this_geom,this_order)%output%mv_mult( array(m,n1:n2) )
#else
          return_value(m,n1:n2) = &
            matmul( outerp(this_geom,this_order)%output%mat , array(m,n1:n2) )
#endif
        end do
      else
        do n = n1,n2
          do m = 1,nvar
            return_value(m,n) = array(m,n)
          end do
        end do
      end if
      !
    end do
    !
  else
    !
    ! No interpolation, just copy the array
    !
    do n = 1,npts
      do m = 1,nvar
        return_value(m,n) = array(m,n)
      end do
    end do
    !
  end if
  ! 
end function interpolate_for_output_2d
!
!###############################################################################
!
pure function vorticity_cv(uu) result(return_value)
  !
  !.. Use Statements ..
  use order_mod, only : maxpts
  use geovar, only : nr,ncell,cell
  use metrics_mod, only : jacobian,metrics
  use derivatives_mod, only : deriv
  !
  real(wp), dimension(:,:), intent(in) :: uu
  !
  real(wp), dimension(1:size(uu,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,np,k
  integer :: this_geom,this_order
  real(wp) :: dudr,duds,dvdr,dvds
  real(wp) :: drdx,dsdx,drdy,dsdy
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr,1:maxpts) :: vel
  !
continue
  !
  do n = 1,ncell
    !
    n1 = cell(n)%beg_sp ! beginning index for sol-pts in cell n
    n2 = cell(n)%end_sp ! ending    index for sol-pts in cell n
    np = n2-n1+1
    !
    this_geom = cell(n)%geom
    this_order = cell(n)%order
    !
    do k = n1,n2
      vel(1:nr,k-n1+1) = uu(nmb:nme,k)/uu(1,k)
    end do
    !
    do k = n1,n2
      !
#ifdef USE_MATRIX_PROCEDURES
      dudr = deriv(this_geom,this_order)%d(1)%dot( k-n1+1 , vel(1,1:np) )
      duds = deriv(this_geom,this_order)%d(2)%dot( k-n1+1 , vel(1,1:np) )
      dvdr = deriv(this_geom,this_order)%d(1)%dot( k-n1+1 , vel(2,1:np) )
      dvds = deriv(this_geom,this_order)%d(2)%dot( k-n1+1 , vel(2,1:np) )
#else
      dudr = dot(deriv(this_geom,this_order)%d(1)%mat(1:np,k-n1+1),vel(1,1:np))
      duds = dot(deriv(this_geom,this_order)%d(2)%mat(1:np,k-n1+1),vel(1,1:np))
      dvdr = dot(deriv(this_geom,this_order)%d(1)%mat(1:np,k-n1+1),vel(2,1:np))
      dvds = dot(deriv(this_geom,this_order)%d(2)%mat(1:np,k-n1+1),vel(2,1:np))
#endif
      !
      drdx = metrics(1,1,k,2)/jacobian(k) ! =  dyds/J
      dsdx = metrics(2,1,k,2)/jacobian(k) ! = -dydr/J
      drdy = metrics(1,2,k,2)/jacobian(k) ! = -dxds/J
      dsdy = metrics(2,2,k,2)/jacobian(k) ! =  dxdr/J
      !
      ! 2D vorticity = d/dx[v(r,s)] - d/dy[u(r,s)]
      !
      !         and using chain rule gives
      !
      ! 2D vorticity = ( (dr/dx)*(dv/dr) + (ds/dx)*(dv/ds) ) - 
      !                ( (dr/dy)*(du/dr) + (ds/dy)*(du/ds) )
      !
      return_value(k) = (drdx*dvdr + dsdx*dvds) - &
                        (drdy*dudr + dsdy*duds)
      !                   
    end do
    !
  end do
  !
end function vorticity_cv
!
!###############################################################################
!
pure function vorticity_pv(vv) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : ncell,cell
  use metrics_mod, only : jacobian,metrics
  use derivatives_mod, only : deriv
  !
  real(wp), dimension(:,:), intent(in) :: vv
  !
  real(wp), dimension(1:size(vv,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,np,k
  integer :: this_geom,this_order
  real(wp) :: dudr,duds,dvdr,dvds
  real(wp) :: drdx,dsdx,drdy,dsdy
  !
continue
  !
  do n = 1,ncell
    !
    n1 = cell(n)%beg_sp ! beginning index for sol-pts in cell n
    n2 = cell(n)%end_sp ! ending    index for sol-pts in cell n
    np = n2-n1+1
    !
    this_geom = cell(n)%geom
    this_order = cell(n)%order
    !
    do k = n1,n2
      !
#ifdef USE_MATRIX_PROCEDURES
      dudr = deriv(this_geom,this_order)%d(1)%dot( k-n1+1 , vv(2,n1:n2) )
      duds = deriv(this_geom,this_order)%d(2)%dot( k-n1+1 , vv(2,n1:n2) )
      dvdr = deriv(this_geom,this_order)%d(1)%dot( k-n1+1 , vv(3,n1:n2) )
      dvds = deriv(this_geom,this_order)%d(2)%dot( k-n1+1 , vv(3,n1:n2) )
#else
      dudr = dot(deriv(this_geom,this_order)%d(1)%mat(1:np,k-n1+1),vv(2,n1:n2))
      duds = dot(deriv(this_geom,this_order)%d(2)%mat(1:np,k-n1+1),vv(2,n1:n2))
      dvdr = dot(deriv(this_geom,this_order)%d(1)%mat(1:np,k-n1+1),vv(3,n1:n2))
      dvds = dot(deriv(this_geom,this_order)%d(2)%mat(1:np,k-n1+1),vv(3,n1:n2))
#endif
      !
      drdx = metrics(1,1,k,2)/jacobian(k) ! =  dyds/J
      dsdx = metrics(2,1,k,2)/jacobian(k) ! = -dydr/J
      drdy = metrics(1,2,k,2)/jacobian(k) ! = -dxds/J
      dsdy = metrics(2,2,k,2)/jacobian(k) ! =  dxdr/J
      !
      ! 2D vorticity = d/dx[v(r,s)] - d/dy[u(r,s)]
      !
      !         and using chain rule gives
      !
      ! 2D vorticity = ( (dr/dx)*(dv/dr) + (ds/dx)*(dv/ds) ) - 
      !                ( (dr/dy)*(du/dr) + (ds/dy)*(du/ds) )
      !
      return_value(k) = (drdx*dvdr + dsdx*dvds) - &
                        (drdy*dudr + dsdy*duds)
      !                   
    end do                
    !
  end do
  !
end function vorticity_pv
!
!###############################################################################
!
pure function vorticity_dusp(usp,dusp,dir) result(return_value)
  !
  integer,                    intent(in) :: dir
  real(wp), dimension(:,:),   intent(in) :: usp
  real(wp), dimension(:,:,:), intent(in) :: dusp
  !
  real(wp), dimension(1:size(dusp,dim=3)) :: return_value
  !
  !.. Local Scalars ..
  integer :: i,j,k,n1,n2
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(dusp,dim=1),1:size(dusp,dim=2)) :: dpvdx
  !
  !.. Local Parameters ..
  integer, parameter, dimension(3) :: idx = [2,3,1]
  integer, parameter, dimension(3) :: jdx = [3,1,2]
  !
continue
  !
  n1 = size(dusp,dim=1)
  n2 = size(dusp,dim=2)
  !
  i = idx(dir)
  j = jdx(dir)
  !
  do k = 1,size(dusp,dim=3)
    dpvdx(1:n1,1:n2) =  grad_cv_to_grad_pv_sp( usp(:,k) , dusp(:,:,k) )
    return_value(k) = dpvdx(i,nmb-1+j) - dpvdx(j,nmb-1+i)
  end do
  !
end function vorticity_dusp
!
!###############################################################################
!
include "Functions/last.f90"
include "Functions/vortex.f90"
include "Functions/bc_priority.f90"
include "Functions/entropy.f90"
include "Functions/pressure.f90"
include "Functions/mach_number.f90"
include "Functions/cell_vol.f90"
include "Functions/bc_string_to_integer.f90"
include "Functions/grad_cv_to_grad_pv.f90"
!
!###############################################################################
!
end module io_mod
!
!###############################################################################
!
