module io_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmb,nme,nee,ntk,ntl,nq
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
  public :: read_restart_file
  public :: write_restart_file
  public :: write_final_error
  public :: write_nodal_parallel_cgns
  public :: write_results_file
  public :: check_flux_point_orderings
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
  character(len=170), public, save :: fname_nodal = "solution.nodal.cgns"
  character(len=170), public, save :: fname_modal = "solution.modal.cgns"
  character(len=170), public, save :: error_file = "error.dat"
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
  !
  real(sp),          save, allocatable, dimension(:) :: soltime
  real(sp),          save, allocatable, dimension(:) :: modtime
  character(len=32), save, allocatable, dimension(:) :: solname
  character(len=32), save, allocatable, dimension(:) :: modname
  !
  !.. CGNS Kind Types ..
  integer, save :: cgns_kind_real
  !
  !.. CGNS Parallel Point Indices ..
  integer, save :: n_my_output_pts
  integer, save :: my_pts_beg
  integer, save :: my_pts_end
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
  use ovar,   only : velref,aref,uref
  use ovar,   only : Period_Timesteps,Period_Time
  use ovar,   only : Grid_Length
  use ovar,   only : LDG_beta,LDG_tau
  !
  use quadrature_mod, only : init_geom_quadrature_rules
  use quadrature_mod, only : init_face_quadrature_rules
  !
  use eqn_idx, only : create_equation_indices
  !
  !.. Local Scalars ..
  integer  :: it1,it2,ierr,ival
  integer  :: Period_Timesteps_X,Period_Timesteps_Y
  logical  :: use_qdtr_face_pts,mms_error
  real(wp) :: Period_Time_X,Period_Time_Y
  real(wp) :: real_Total_Periods
  real(wp) :: aoa_radians,pres_check
  character(len=1000) :: error_message
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
  ! Append the name of the solution and error files to output_dir
  !
  fname_nodal = trim(output_dir)//trim(adjustl(fname_nodal))
  fname_modal = trim(output_dir)//trim(adjustl(fname_modal))
  error_file = trim(output_dir)//trim(adjustl(error_file))
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
      write (error_message,201)
      mms_error = true
    else if (mms_opt == 3 .and. nr > 2) then ! 3 is for 1D or 2D
      write (error_message,202)
      mms_error = true
    else if (mms_opt == 4 .and. nr /= 2) then ! 4 is for 2D only
      write (error_message,203)
      mms_error = true
    else if (all(mms_opt /= [1,2,3,4])) then
      write (error_message,204)
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
  ! Check the reference conditions for consistency
  !
  if (rhoref < zero) then
    rhoref = pref/(rgasref*tref)
  else if (pref < zero) then
    pref = rhoref*rgasref*tref
  else if (tref < zero) then
    tref = pref/(rgasref*rhoref)
  else
    pres_check = rhoref*rgasref*tref
    if (abs(pres_check-pref) > eps6) then
     !write (iout,102) pref,pres_check
      write (error_message,102) pref,pres_check
      call stop_recon(stop_mpi,pname,error_message)
    end if
  end if
  !
  ! Compute some of the remaining reference conditions
  !
  aref = sqrt( gam * pref / rhoref )
  uref = machref * aref
  velref(1) = uref * cos(beta_aoaref  * degrees2radians) * &
                     cos(alpha_aoaref * degrees2radians)
  velref(2) = uref * cos(beta_aoaref  * degrees2radians) * &
                     sin(alpha_aoaref * degrees2radians)
  velref(3) = uref * sin(beta_aoaref  * degrees2radians)
  !
  ! Redefine the reference conditions for the linear advection problems
  ! and the original Shu vortex problem because using the dimensional speed
  ! of sound for nondimensionalizing velocity will cause issues
  !
  if (any(itestcase == Advec_Problems) .or. &
      abs(itestcase) == Shu_Vortex .or. &
      itestcase == -Vortex_Transport) then
    !
    ! Alteration to make the Shu case propagate in the
    ! y-direction only instead of in both x and y directions
    !
    if (itestcase == -Shu_Vortex) then
      alpha_aoaref = 90.0_wp
      beta_aoaref = 0.0_wp
      velref(1) = zero
      velref(2) = zero
      velref(3) = zero
      itestcase = abs(itestcase)
    else if (itestcase == -Vortex_Transport) then
      alpha_aoaref = 90.0_wp
      beta_aoaref = 0.0_wp
      velref(1) = zero
      velref(2) = one
      velref(3) = zero
      itestcase = Shu_Vortex
    else
      alpha_aoaref = 45.0_wp
      beta_aoaref = 0.0_wp
      velref(1) = one
      velref(2) = one
      velref(3) = zero
    end if
    !
    machref = one/sqrt(gam)
    aref = one
    uref = one
    rhoref = one
    pref = one
    tref = one
    rgasref = one
    !
  end if
  !
  ! Initialize the remaining reference conditions
  !
  call initialize_reference
  !
  ! Initialize the viscosity model
  !
  call initialize_viscosity
  !
  ! Initialize the turbulence equations
  !
  call initialize_turbulence
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
      Period_Time_Y = merge( Period_Time_Y , Grid_Length(2) , Period_Time_Y > eps3 )
      Period_Timesteps_Y = nint(Period_Time_Y/dt_in)
      Period_Timesteps = Period_Timesteps_Y
      !
     !if (mypnum == 0) write (*,*) Period_Time_Y,Period_Timesteps_Y,Period_Timesteps
      !
    else if (abs(velref(2)) < eps3) then
      !
      Period_Time_X = Grid_Length(1) * reciprocal(velref(1)/aref)
      Period_Time_X = merge( Period_Time_X , Grid_Length(1) , Period_Time_X > eps3 )
      Period_Timesteps_X = nint(Period_Time_X/dt_in)
      Period_Timesteps = Period_Timesteps_X
      !
     !if (mypnum == 0) write (*,*) Period_Time_Y,Period_Timesteps_Y,Period_Timesteps
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
  else
    !
    cfl_end = max(cfl_end,eps6)    ! cfl_end must be >= 1.0e-6
    cfl_cycles = max(cfl_cycles,1) ! cfl_cycles must be at least 1
    cfl_beg = min(cfl_beg,cfl_end) ! cfl_beg must be <= cfl_end
    !
    if (Final_Time > zero) then
      num_ts = nint( Final_Time / dt_in )
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
  102 format ( &
     " Free-stream conditions in the input file are not consistent!",/, &
     " Reference pressure given    : ",es14.6,/, &
     " Reference pressure computed : ",es14.6,/, &
     " Please correct the input file or set one of density, pressure,", &
     " and temperature to a negative value so it is explicitly computed.")
  !
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
  201 format ("MMS ERROR: MMS solutions 1 and 2 are allowed for 2D and 3D only!")
  202 format ("MMS ERROR: MMS solution 3 is allowed for 1D and 2D only!")
  203 format ("MMS ERROR: MMS solution 4 is allowed for 2D only!")
  204 format ("MMS ERROR: The MMS solution chosen is invalid!")
  !
end subroutine read_nml_input_file
!
!###############################################################################
!
subroutine check_output_dir
  !
  integer :: clen,ierr,iodir
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
 !write (iout,*)
 !write (iout,*) "Checking output directory"
 !write (iout,*)
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
    return
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
    ! set output_dir to test_dir, and return
    !
    if (ierr == 0) then
      close (iodir,status="delete",iostat=ierr)
      output_dir = trim(test_dir)
      return
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
    ! set output_dir to test_dir, and return
    !
    if (ierr == 0) then
      close (iodir,status="delete",iostat=ierr)
      output_dir = trim(test_dir)
      return
    end if
    !
  end if
  !
  ! If we have made it here, we cant seem to find the output directory
  ! that was given in the input file so exit
  !
  write (iout,1) trim(output_dir)
  call stop_recon(stop_mpi,pname)
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
subroutine initialize_viscosity
  !
  !.. Use Statements ..
  use ovar, only : muref, muref_nd, reyref
  use ovar, only : pv_in, rhoref, aref, machref, tref
  use ovar, only : suth_muref, suth_c1
  use ovar, only : suth_Tref, suth_Sref
  !
  !.. Local Scalars ..
  real(wp) :: muref_scale,tref_scale
  !
continue
  !
  muref_scale = (one + suth_Sref/suth_Tref) / sqrt(suth_Tref)
  !
  if (reyref <= zero) then
    !
    suth_c1 = suth_muref * muref_scale
    !
    muref_nd = viscosity_pv_sp(pv_in)
    muref    = muref_nd * rhoref * aref
    !
    reyref = rhoref * machref * aref * reciprocal( muref )
    !
  else
    !
    muref = rhoref * machref * aref / reyref
    !
    muref_nd = muref / (rhoref * aref)
    !
    tref_scale = tref * sqrt(tref) / (tref + suth_Sref)
    !
    suth_muref = muref / muref_scale / tref_scale
    !
    suth_c1 = suth_muref * muref_scale
    !
  end if
  !
end subroutine initialize_viscosity
!
!###############################################################################
!
subroutine initialize_reference
  !
  !.. Use Statements ..
  use ovar,   only : rhoref,aref,tref,gam
  use ovar,   only : velref,velref_nd
  use ovar,   only : pref,pref_nd,tmin,tmax
  use ovar,   only : rgasref,rgasref_nd
  use ovar,   only : cpref,cpref_nd
  use ovar,   only : cv_in,pv_in,machref
  use ovar,   only : ttotref,ttotref_nd
  use ovar,   only : ptotref,ptotref_nd
  use ovar,   only : ptot2p_ratio
  !
  !.. Local Scalars ..
  integer  :: ierr
  real(wp) :: rhoa,rhoa2,tda2
  real(wp) :: gam_over_gm1
  real(wp) :: gm1_over_gam
  real(wp) :: total_cond
  !
  character(len=1000) :: error_message
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_reference"
  !
continue
  !
  rhoa = rhoref * aref
  rhoa2 = rhoa * aref
  tda2 = tref / aref / aref
  gam_over_gm1 = gam / (gam-one)
  gm1_over_gam = one / gam_over_gm1
  !
  ! Allocate the free stream inflow arrays
  !
  allocate ( cv_in(1:nq) , source=zero , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cv_in",1,ierr,error_message)
  allocate ( pv_in(1:nq) , source=zero , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"pv_in",1,ierr,error_message)
  !
  ! Initialize the inflow boundary conditions to the reference conditions
  ! NOTE: This is hard-wired for 2D and needs to be generalized for 3D
  !
  pv_in(nec) = rhoref
  pv_in(nmb:nme) = velref
  pv_in(nee) = pref
  !
  cv_in(nec) = pv_in(nec)
  cv_in(nmb:nme) = pv_in(nec) * pv_in(nmb:nme)
  cv_in(nee) = pv_in(nee)/(gam-one) + half*pv_in(nec)*selfdot(pv_in(nmb:nme))
  !
  ! Nondimensional reference velocities
  velref_nd = velref / aref
  !
  ! Nondimensional reference pressure
  pref_nd = pref / rhoa2
  !
  ! Nondimensional reference gas constant
  rgasref_nd = rgasref * tda2
  !
  ! Dimensional reference specific heat at constant pressure
  cpref    = rgasref * gam / (gam - one)
  ! Nondimensional reference specific heat at constant pressure
  cpref_nd = cpref * tda2
  !
  ! Finally, nondimensionalize the inflow boundary conditions
  !
  cv_in(nec)     = cv_in(nec) / rhoref
  cv_in(nmb:nme) = cv_in(nmb:nme) / rhoa
  cv_in(nee)     = cv_in(nee) / rhoa2
  !
  pv_in(nec)     = pv_in(nec) / rhoref
  pv_in(nmb:nme) = pv_in(nmb:nme) / aref
  pv_in(nee)     = pv_in(nee) / rhoa2
  !
  ! Total pressure and temperature
  !
  if (ptot2p_ratio < zero) then
    !
    total_cond = one + half*(gam-one)*machref*machref
    !
    ttotref = tref * total_cond
    ptotref = pref * (total_cond**gam_over_gm1)
    !
    ptot2p_ratio = ptotref/pref
    !
  else
    !
    ttotref = tref * (ptot2p_ratio**gm1_over_gam)
    ptotref = pref * ptot2p_ratio
    !
  end if
  !
  ttotref_nd = ttotref / tref
  ptotref_nd = ptotref / rhoa2
  !
  tmin = 10.0_wp
  tmax = 15000.0_wp
  !
end subroutine initialize_reference
!
!###############################################################################
!
subroutine initialize_turbulence()
  !
  !.. Use Statements ..
  use ovar, only : tkeref, tkeref_nd
  use ovar, only : tlsref, tlsref_nd
  use ovar, only : tvrref, tinref
  use ovar, only : rhoref, muref
  use ovar, only : aref, machref
  use ovar, only : cv_in, pv_in
  use ovar, only : turbulence_model
  use ovar, only : turb_is_on
  use ovar, only : include_tke_in_energy
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_turbulence"
  !
continue
  !
  ! Turbulence equations
  !
  if (turbulence_model /= laminar_flow) then
    !
    ! Set the flag indicating that turbulence is active
    !
    turb_is_on = true
    !
    ! Set flag that says whether TKE is included in the total energy
    !
    if (abs(turbulence_model) == k_omega) then
      include_tke_in_energy = true
    end if
    !
    ! Dimensional reference conditions
    !
    tkeref = 1.5_wp * (tinref*machref*aref)**2
    tlsref = rhoref * tkeref / (tvrref*muref)
    !
    ! Nondimensional reference conditions
    !
    tkeref_nd = tkeref / (rhoref * aref * aref)
    tlsref_nd = tlsref / (rhoref * aref)
    !
    ! Set the turbulence conserved variables inflow conditions
    !
    cv_in(ntk) = cv_in(nec) * tkeref_nd
    cv_in(ntl) = cv_in(nec) * tlsref_nd
    !
    ! Add the tke contribution to the energy equation
    !
    cv_in(nee) = cv_in(nee) + cv_in(ntk)
    !
    ! Recompute the primitive variable inflow conditions
    !
    pv_in(1:nq) =  usp2v_sp( cv_in(1:nq) )
    !
  end if
  !
end subroutine initialize_turbulence
!
!###############################################################################
!
subroutine read_grid_file()
  !
  !.. Use Statements ..
  use module_gmsh,   only : read_gmsh_gridfile
  use module_plot3d, only : read_plot3d_gridfile
  use quadrature_mod, only : find_valid_geometries
  use ovar,           only : use_bc_conflicts
  !
  ! Routine that uses the variable grid_format to choose
  ! which routine to use in order to read the grid file.
  !        
continue
  !
  if (grid_format == Gambit_Neutral) then
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
  call check_node_ordering
  !
  ! Adjust the grid coordinates and nodes_of_cell
  ! connectivity to account for ghost cells
  !
  call adjust_for_ghost_cells
  !
  ! Find any grid nodes on the boundaries that separate two boundary
  ! faces with non-matching boundary conditions
  !
  if (use_bc_conflicts) call find_bc_conflicts
  !
  ! Go through the grid and find which geometries are valid for
  ! the current simulation
  !
  call find_valid_geometries
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
  character(len=1000) :: error_message
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
    ! NOTE: switch_index evaluates to nfnods (i.e. 2) for a 2D host cell
    !       switch_index evaluates to nfnods (i.e. 3) for a triangular boundary face
    !       switch_index evaluates to 3 for a quadrilateral boundary face
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
  use geovar,    only : global_pinc_ptr,global_xyz
  use geovar,    only : cell_geom,global_cell_geom
  use geovar,    only : cell_order,global_cell_order
  use flowvar,   only : global_usp
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
  ! Create the global pts_in_cell_ptr. This will be done later
  ! but a global version is never made before the global grid
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
    global_pinc_ptr(n+1) = global_pinc_ptr(n) + geom_solpts(this_geom,this_order)
    !
  end do
  !
  n_global_solpts = global_pinc_ptr(ncell+1)
  n_global_totpts = last(global_pinc_ptr)
  !
  ! Now the only remaining global grid arrays to create are global_xyz
  ! and global_usp. We will have to do a collective communication of the
  ! local arrays to accomplish this after they are created. We will only
  ! allocate the array here since we have the parameters needed to do so.
  !
  allocate ( global_xyz(1:nr,1:n_global_totpts) , source=zero , stat=ierr )
  call alloc_error(pname,"global_xyz",1,ierr)
  !
 !allocate ( global_usp(1:nq,1:n_global_totpts) , source=zero , stat=ierr )
 !call alloc_error(pname,"global_usp",1,ierr)
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
  character(len=1000) :: error_message
  integer :: nfnods
  !
  !.. Local Arrays ..
  integer,  dimension(1:4)      :: ip,periodic_ip
  integer,  dimension(1:4)      :: face_nodes
  integer,  dimension(1:4)      :: pface_nodes
  real(wp), dimension(1:nr)     :: ptf
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
   !! NOTE: switch_index evaluates to nfnods (i.e. 2) for a 2D host cell
   !!       switch_index evaluates to nfnods (i.e. 3) for a triangular boundary face
   !!       switch_index evaluates to 3 for a quadrilateral boundary face
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
    ! Conditional to treat periodic boundary conditions separately
    !
    if (bface(1,nf) /= bc_periodic) then
      !
      ! This is not a periodic boundary face
      !
      ! Center of boundary face by averaging the coordinates of the two boundary nodes
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
        ! Get the coordinates of ghost node 3
        !
        xyz_nodes(:,ng) = merge(rotate(xyz_nodes(:,h1),xyz_nodes(:,h2),xyz_nodes(:,ni)), &
                                reflect(xyz_nodes(:,ni),ptf,rnrm), &
                                rotation == 1)
        !
      else if (host_geom == 4) then
        !
        ! Ghost node 3 for quadrilateral
        !
        ng = nodes_of_cell( g1 + 2 )
        ni = nodes_of_cell( h1 + mod(host_side+2-max(0,min(1,rotation)),host_geom) )
        !
        ! Get the coordinates of ghost node 3
        !
        xyz_nodes(:,ng) = merge(rotate(xyz_nodes(:,h1),xyz_nodes(:,h2),xyz_nodes(:,ni)), &
                                reflect(xyz_nodes(:,ni),ptf,rnrm), &
                                rotation == 1)
        !
        ! Ghost node 4 for quadrilateral
        !
        ng = nodes_of_cell( g1 + 3 )
        ni = nodes_of_cell( h1 + mod(host_side+1+max(0,min(1,rotation)),host_geom) )
        !
        ! Get the coordinates of ghost node 4
        !
        xyz_nodes(:,ng) = merge(rotate(xyz_nodes(:,h1),xyz_nodes(:,h2),xyz_nodes(:,ni)), &
                                reflect(xyz_nodes(:,ni),ptf,rnrm), &
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
      ! Nodes on the boundary face for periodic_cell ! NOW STORED IN pface_nodes(1:nfnods)
      !
      ! Nodes on the boundary face for host_cell ! NOW STORED IN face_nodes(1:nfnods)
      !
      ! ##########################################################################
      ! ##########################################################################
      ! ##########################################################################
      ! ###############  ||||||||||                   ||||||||||  ################
      ! ###############  VVVVVVVVVV  NEED TO WORK ON  VVVVVVVVVV  ################
      ! ##########################################################################
      ! ##########################################################################
      ! ##########################################################################
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
      ! ##########################################################################
      ! ##########################################################################
      ! ##########################################################################
      ! ###############  AAAAAAAAAA  NEED TO WORK ON  AAAAAAAAAA  ################
      ! ###############  ||||||||||                   ||||||||||  ################
      ! ##########################################################################
      ! ##########################################################################
      ! ##########################################################################
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
  character(len=1000) :: error_message
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
  allocate ( bface_nodes(1:4,1:nfbnd) , source=0 , stat=ierr )
  call alloc_error(pname,"bface_nodes",1,ierr)
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
  allocate ( faces_with_node_ptr(1:nnode+1) , source=0 , stat=ierr )
  call alloc_error(pname,"faces_with_node_ptr",1,ierr)
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
  allocate ( faces_with_node(1:last(faces_with_node_ptr)) , source=0 , stat=ierr )
  call alloc_error(pname,"faces_with_node",1,ierr)
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
  allocate ( mskp(1:nnode) , source=fals , stat=ierr )
  call alloc_error(pname,"mskp",1,ierr)
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
  allocate ( bc_conflict(1:count(mskp)) , stat=ierr )
  call alloc_error(pname,"bc_conflict",1,ierr)
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
    allocate ( bc_conflict(i)%slave(1:nf2-nf1) , stat=ierr )
    write (array_name,1) "bc_conflict(",i,")%slave"
    call alloc_error(pname,array_name,1,ierr)
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
      write (iout,*) "ERROR: Conflict node doesn't seem to match a boundary node!"
    end if
    !
    ! Now do the same for the slave faces
    !
    do j = 1,size(bc_conflict(i)%slave)
      !
      nf = bc_conflict(i)%slave(j)%face
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
        bc_conflict(i)%slave(j)%flux_point = 1
      else if (n == bface_nodes(2,nf)) then
        bc_conflict(i)%slave(j)%flux_point = n_fp1de
      else
        write (iout,*) "ERROR: Conflict node doesn't seem to match a boundary node!"
      end if
      !
    end do
    !
  end do
  !
  ! Deallocate the local arrays before leaving
  !
  deallocate ( bface_nodes , stat=ierr )
  call alloc_error(pname,"bface_nodes",2,ierr)
  deallocate ( faces_with_node_ptr , stat=ierr )
  call alloc_error(pname,"faces_with_node_ptr",2,ierr)
  deallocate ( faces_with_node , stat=ierr )
  call alloc_error(pname,"faces_with_node",2,ierr)
  deallocate ( mskp , stat=ierr )
  call alloc_error(pname,"mskp",2,ierr)
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
  !
end subroutine find_bc_conflicts
!
!###############################################################################
!
subroutine read_restart_file
  !
  !.. Use Statements ..
  use geovar,       only : ncell,pts_in_cell_ptr
  use geovar,       only : n_global_cell,global_nofc_ptr,global_pinc_ptr
  use ovar,         only : itrst,tmrst,dtrst
  use flowvar,      only : usp,global_usp
  use parallel_mod, only : cell_map
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,gp,ierr
  integer :: iorst
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_restart_file"
  character(len=*), parameter :: fname = "out.rst"
  !
continue
  !
  open(newunit=iorst,file=fname,status="old",action="read", &
                     form="unformatted",iostat=ierr)
  call io_error(pname,fname,1,ierr)
  !
  read (iorst) itrst, tmrst, dtrst
  read (iorst) global_usp
  !
  close (iorst,iostat=ierr)
  call io_error(pname,fname,2,ierr)
  !
  ! If using multiple processors, copy the solution data from 
  ! only those cells that are local to this processor into usp
  !
  if (ncpu > 1) then
    !
    do nc = 1,ncell
      !
      ! Indices of the datat for this cell in the usp array
      !
      n1 = pts_in_cell_ptr(nc)+1
      n2 = pts_in_cell_ptr(nc+1)
      np = n2-n1+1
      !
      ! Starting point of the data for this cell in the global_usp array
      !
      gp = global_pinc_ptr( cell_map(mypnum+1)%loc_to_glb(nc) )
      !
      ! Copy the solution data from the global_usp array for this cell
      !
      usp(:,n1:n2) = global_usp(:,gp+1:gp+np)
      !
    end do
    !
  else
    !
    ! For a serial job, just copy the global_usp array into usp
    !
    usp(:,:) = global_usp(:,:)
    !
  end if
  !
end subroutine read_restart_file
!
!###############################################################################
!
subroutine write_restart_file(itnum)
  !
  !.. Use Statements ..
  use ovar,    only : itcur,time,d_t
  use flowvar, only : global_usp
  use parallel_mod, only : collect_global_solution
  !
  !.. Formal Arguments ..
  integer, intent(in) :: itnum
  !
  !.. Local Scalars ..
  integer :: ierr,iorst
  character(len=40) :: fname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_restart_file"
  !
continue
  !
  ! Collect the global solution
  !
  call collect_global_solution
  !
  ! Only use the root partition to output the restart file
  !
  if (mypnum == 0) then
    !
    if (itnum > 0) then
      write (fname,'("out.",i0.8,".rst")') itcur
    else
      write (fname,'("out.",i0.8,".done.rst")') itcur
    end if
    !
    open (newunit=iorst,file=trim(adjustl(fname)),status="replace", &
                        action="write",form='unformatted',iostat=ierr)
    call io_error(pname,trim(adjustl(fname)),1,ierr)
    !
    write (iorst) itcur, time, d_t
    write (iorst) global_usp
    !
    close (iorst,iostat=ierr)
    call io_error(pname,trim(adjustl(fname)),2,ierr)
    !
  end if
  !
  call mpi_barrier(mpi_comm_world,mpierr)
  !
end subroutine write_restart_file
!
!###############################################################################
!
subroutine write_final_error
  !
  use geovar, only : n_global_cell
  use geovar, only : n_global_solpts,global_pinc_ptr
  use ovar,   only : itestcase,flux_divergence
  use ovar,   only : vnloc,vkloc,vlone,vlinf,vltwo
  use ovar,   only : unloc,ukloc,ulone,ulinf,ultwo
  use order_mod, only : n_order
  !
  !.. Aliased Use Statements ..
  use geovar, only : xyz => global_xyz
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_final_error"
  !
continue
  !
  write (iout,*)
  !
  write (iout,34) itestcase,flux_divergence, &
                  n_order,n_global_cell,n_global_solpts, &
                  vnloc(1),global_pinc_ptr(vnloc(1))+vkloc(1), &
                  xyz(1,global_pinc_ptr(vnloc(1))+vkloc(1)), &
                  xyz(2,global_pinc_ptr(vnloc(1))+vkloc(1))
  !
  write (iout,35) (vlone(i), &
                   vlinf(i), &
                   vltwo(i), &
                   xyz(1,global_pinc_ptr(vnloc(i))+vkloc(i)), &
                   xyz(2,global_pinc_ptr(vnloc(i))+vkloc(i)), i=1,nq-1)
  write (iout,36) (ulone(i), &
                   ulinf(i), &
                   ultwo(i), &
                   xyz(1,global_pinc_ptr(unloc(i))+ukloc(i)), &
                   xyz(2,global_pinc_ptr(unloc(i))+ukloc(i)), i=nq+1,nq+1)
  write (iout,37) (vlone(i), &
                   vlinf(i), &
                   vltwo(i), &
                   xyz(1,global_pinc_ptr(vnloc(i))+vkloc(i)), &
                   xyz(2,global_pinc_ptr(vnloc(i))+vkloc(i)), i=nq+1,nq+1)
  !
  34 format ('itestcase =',i2,/, &
             'flux_divergence  =',i2,//, &
             'Order of solution polynomial = ',i10,/, &
             'Number of grid cells         = ',i10,/, &
             'Number of solution points    = ',i10,//, &
             'Max density error in cell = ',i10,/, &
             '        at solution point = ',i10,/, &
             '         with coordinates = ',2es14.6,/)
  35 format (/,'  Error  -      L1           Max            L2     @      x             y',/,&
        ' ------------------------------------------------------------------------------',/,&
        ' Density ',5es14.6,/,&
        '   Vx    ',5es14.6,/,&
        '   Vy    ',5es14.6)
  36 format (' Entropy ',5es14.6)
  37 format (' Vel Mag ',5es14.6,/,&
        ' ------------------------------------------------------------------------------',/)
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
subroutine write_nodal_parallel_cgns(iflag)
  !
  !.. Use Statements ..
  use geovar,  only : nr,n_solpts,xyz
  use ovar,    only : time,aref,uref,rhoref
  use ovar,    only : itestcase,cgns_use_queue
  use ovar,    only : interpolate_before_output_variable
  use flowvar, only : usp,dusp
  !
  use geovar,  only : geovar_alloc_status
  use ovar,    only : ovar_alloc_status
  use flowvar, only : flowvar_alloc_status
  use metrics_mod,  only : metrics_alloc_status
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iflag
  !
  !.. Local Scalars ..
  integer :: l,m,n,ierr
  integer :: my_cells_beg,my_cells_end
  integer :: n_total_output_pts
  integer :: n_total_output_cells
  !
  character(len=32) :: tname
  character(len=32) :: varname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_nodal_parallel_cgns"
  character(len=*), parameter :: char_xyz(1:3) = ["X","Y","Z"]
  character(len=*), parameter :: char_sname(1:6) = ["Density  ", &
                                                    "Velocity ", &
                                                    "Pressure ", &
                                                    "Mach     ", &
                                                    "Entropy  ", &
                                                    "Vorticity" ]
  !
  !.. Local Arrays ..
  integer, dimension(1:2) :: idata
  integer, dimension(1:3) :: idim1
  integer, dimension(1:nr) :: ixyz_n
  integer, dimension(1:(nme-nec)*nr) :: idvar_n
  integer, dimension(1:size(char_sname)+nr-1) :: isvar_n
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable, dimension(:,:) :: ipelem
  real(wp), allocatable, dimension(:)   :: vartmp
  real(r4), allocatable, dimension(:)   :: var
  real(wp), allocatable, dimension(:,:) :: interp_usp
  !
continue
  !
  if (iflag == 1) then
    !
    cgns_kind_real = merge(RealDouble,RealSingle,kind(var)==wp)
    !
    ! Initialize the variables and data arrays needed for CGNS output
    ! NOTE: This call to routine get_initial_data_for_cgns will initialize
    !       the module variables my_pts_beg, my_pts_end, and n_my_output_pts.
    !
    call get_initial_data_for_cgns(my_cells_beg,my_cells_end,idim1,ipelem)
    !
    ! Extract the total number of pts and cells from the idim1 array
    !
    n_total_output_pts = idim1(1)
    n_total_output_cells = idim1(2)
    !
    ! Allocate the temporary variable arrays
    !
    allocate ( vartmp(1:n_solpts) , source=zero , stat=ierr )
    call alloc_error(pname,"vartmp",1,ierr)
    allocate ( var(1:n_my_output_pts) , source=real(0,r4) , stat=ierr )
    call alloc_error(pname,"var",1,ierr)
    !
    ! Open the cgns output file if iflag is 1
    !
    ! Initialize the zone counter
    !
    cgnsnodalzone_counter = 0
    !
    call cgp_open_f(fname_nodal,CG_MODE_WRITE,ifile_n,ierr)
    call io_error(pname,trim(fname_nodal),1,ierr)
    !
    ! Set the CGNS base (who knows what this is for)
    !
    call cg_base_write_f(ifile_n,"recon_base",nr,nr,ibase_n,ierr)
    call cgns_error(pname,fname_nodal,"cg_base_write_f",ierr)
    !
    ! Write the title for this solution
    !
    call cg_zone_write_f(ifile_n,ibase_n,"Solution",idim1, &
                         Unstructured,izone_n,ierr)
    call cgns_error(pname,fname_nodal,"cg_zone_write_f",ierr)
    !
    ! Create data nodes for the coordinates
    !
    do l = 1,nr
      varname = "Coordinate"//char_xyz(l)
      call cgp_coord_write_f(ifile_n,ibase_n,izone_n,cgns_kind_real, &
                             trim(varname),ixyz_n(l),ierr)
      call cgns_error(pname,fname_nodal,"cgp_coord_write_f",ierr, &
                      "Write node: Coordinate",n)
    end do
    !
    ! Now create and write the coordinate data in parallel
    !
    if (.not. cgns_use_queue) then
      call cgp_queue_set_f(0,ierr)
      call cgns_error(pname,fname_nodal,"cgp_queue_set_f",ierr)
    end if
    !
    do l = 1,nr
      !
      ! Save the current coordinate to the local array vartmp
      !
      vartmp(1:n_solpts) = xyz(l,1:n_solpts)
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
        call cgns_error(pname,fname_nodal,"cgp_queue_set_f",ierr)
      end if
      !
      ! Add the data for current coordinate variable to the output queue.
      !
      call cgp_coord_write_data_f(ifile_n,ibase_n,izone_n,ixyz_n(l), &
                                  my_pts_beg,my_pts_end,var,ierr)
      call cgns_error(pname,fname_nodal,"cgp_coord_write_data_f",ierr, &
                      "Write data: Coordinate",n)
      !
      ! Finally, flush the queue containing the
      ! current coordinate variable to the CGNS file.
      !
      if (cgns_use_queue) then
        call cgp_queue_flush_f(ierr)
        call cgns_error(pname,fname_nodal,"cgp_queue_flush_f",ierr)
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
    call cgns_error(pname,fname_nodal,"cgp_section_write_f",ierr)
    !
    ! Write the element connectivity in parallel
    !
    ! Reset the queue.
    !
    if (cgns_use_queue) then
      call cgp_queue_set_f(1,ierr)
      call cgns_error(pname,fname_nodal,"cgp_queue_set_f",ierr)
    end if
    !
    ! Add the element connectivity data to the output queue.
    !
    call cgp_elements_write_data_f(ifile_n,ibase_n,izone_n,isect_n, &
                                   my_cells_beg,my_cells_end,ipelem,ierr)
    call cgns_error(pname,fname_nodal,"cgp_elements_write_data_f",ierr)
    !
    ! Finally, flush the queue containing the
    ! element connectivity to the CGNS file.
    !
    if (cgns_use_queue) then
      call cgp_queue_flush_f(ierr)
      call cgns_error(pname,fname_nodal,"cgp_queue_flush_f",ierr)
    end if
    !
    ! Deallocate ipelem since we are done with it
    !
    deallocate ( ipelem , stat=ierr )
    call alloc_error(pname,"ipelem",2,ierr)
    !
  else if (iflag == 2) then
    !
    ! The CGNS file was dumped last time this subroutine was called
    ! so the CGNS file needs to be reopened in modification mode.
    !
    call cgp_open_f(fname_nodal,CG_MODE_MODIFY,ifile_n,ierr)
    call io_error(pname,trim(fname_nodal),3,ierr)
    !
  end if
  !
  ! Allocate var and vartmp if they arent already allocated
  !
  if (.not. allocated(vartmp)) then
    allocate ( vartmp(1:n_solpts) , source=zero , stat=ierr )
    call alloc_error(pname,"vartmp",1,ierr)
  end if
  if (.not. allocated(var)) then
    allocate ( var(1:n_my_output_pts) , source=real(0,r4) , stat=ierr )
    call alloc_error(pname,"var",1,ierr)
  end if
  if (.not. allocated(interp_usp)) then
    allocate ( interp_usp(1:nq,1:n_solpts) , source=zero , stat=ierr )
    call alloc_error(pname,"interp_usp",1,ierr)
  end if
  !
  ! Get the name for the current solution node
  !
  write (tname,'("FlowSolution",i0)') cgnsnodalzone_counter
  !
  ! Reallocate the soltime and solname arrays so we
  ! can concatenate this solution onto those arrays
  !
  call reallocate(soltime,cgnsnodalzone_counter+1)
  call reallocate(solname,cgnsnodalzone_counter+1)
  soltime(size(soltime)) = real(time,kind(soltime))
  solname(size(solname)) = tname
  !
  ! Write the flow variables
  !
  call cg_sol_write_f(ifile_n,ibase_n,izone_n,tname,Vertex,isol_n,ierr)
  call cgns_error(pname,fname_nodal,"cg_sol_write_f",ierr)
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
      call cgp_field_write_f(ifile_n,ibase_n,izone_n,isol_n,cgns_kind_real, &
                             trim(varname),idvar_n(l),ierr)
      call cgns_error(pname,fname_nodal,"cgp_field_write_f",ierr, &
                      "Write node: Gradient variable",l)
    end do
    !
    ! Velocity Gradients
    !
    do m = 1,nr
      do l = 1,nr
        n = m*nr + l
        varname = "Velocity"//char_xyz(m)//"_d"//char_xyz(l)
        call cgp_field_write_f(ifile_n,ibase_n,izone_n,isol_n,cgns_kind_real, &
                               trim(varname),idvar_n(n),ierr)
        call cgns_error(pname,fname_nodal,"cgp_field_write_f",ierr, &
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
          call cgns_error(pname,fname_nodal,"cgp_queue_set_f",ierr)
        end if
        !
        ! Add the data for current variable to the output queue.
        !
        call cgp_field_write_data_f(ifile_n,ibase_n,izone_n,isol_n,idvar_n(n), &
                                    my_pts_beg,my_pts_end,var,ierr)
        call cgns_error(pname,fname_nodal,"cgp_field_write_data_f",ierr, &
                        "Write data: Gradient variable",n)
        !
        ! Finally, flush the queue containing
        ! the current variable to the CGNS file.
        !
        if (cgns_use_queue) then
          call cgp_queue_flush_f(ierr)
          call cgns_error(pname,fname_nodal,"cgp_queue_flush_f",ierr)
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
          call cgp_field_write_f(ifile_n,ibase_n,izone_n,isol_n,cgns_kind_real, &
                                 trim(varname),isvar_n(n),ierr)
          call cgns_error(pname,fname_nodal,"cgp_field_write_f",ierr, &
                          "Write node: Solution variable",n)
        end do
      else
        n = n + 1
        call cgp_field_write_f(ifile_n,ibase_n,izone_n,isol_n,cgns_kind_real, &
                               trim(char_sname(l)),isvar_n(n),ierr)
        call cgns_error(pname,fname_nodal,"cgp_field_write_f",ierr, &
                        "Write node: Solution variable",n)
      end if
    end do
    !
    ! Interpolate the conservative variables for output
    !
    if (interpolate_before_output_variable) then
      interp_usp(1:nq,1:n_solpts) = &
                           interpolate_for_output( usp(1:nq,1:n_solpts) )
    else
      interp_usp(1:nq,1:n_solpts) = usp(1:nq,1:n_solpts)
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
            vartmp(1:n_solpts) = interp_usp(m,1:n_solpts)/interp_usp(nec,1:n_solpts)
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
          n = n + 1
         !vartmp = entropy_cv( interp_usp )
          vartmp = vorticity_cv( interp_usp )
          call cgns_write_field_parallel
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
    call cgns_error(pname,fname_nodal,"cgp_queue_set_f",ierr)
  end if
  !
  ! Deallocate the arrays var and vartmp
  !
  if (allocated(var)) then
    deallocate ( var , stat=ierr )
    call alloc_error(pname,"var",2,ierr)
  end if
  if (allocated(vartmp)) then
    deallocate ( vartmp , stat=ierr )
    call alloc_error(pname,"vartmp",2,ierr)
  end if
  if (allocated(interp_usp)) then
    deallocate ( interp_usp , stat=ierr )
    call alloc_error(pname,"interp_usp",2,ierr)
  end if
  !
  ! If last iteration or a solution dump was requested, close CGNS file
  !
  if (any(iflag == [-1,-2])) then
    !
    ! Update the base iterative node
    !
    call cg_biter_write_f(ifile_n,ibase_n,'TimeIterValues',size(soltime),ierr)
    call cgns_error(pname,fname_nodal,"cg_biter_write_f",ierr)
    !
    call cg_goto_f(ifile_n,ibase_n,ierr,"BaseIterativeData_t",1,"end")
    call cgns_error(pname,fname_nodal,"cg_goto_f",ierr)
    !
    call cg_array_write_f("TimeValues",cgns_kind_real,1,size(soltime),soltime,ierr)
    call cgns_error(pname,fname_nodal,"cg_array_write_f",ierr)
    !
    ! Update the zone iterative node
    !
    call cg_ziter_write_f(ifile_n,ibase_n,izone_n,"ZoneIterativeData",ierr)
    call cgns_error(pname,fname_nodal,"cg_ziter_write_f",ierr)
    !
    call cg_goto_f(ifile_n,ibase_n,ierr,"Zone_t",izone_n, &
                                        "ZoneIterativeData_t",1,"end")
    call cgns_error(pname,fname_nodal,"cg_goto_f",ierr)
    !
    idata(:) = [32,size(solname)]
    call cg_array_write_f("FlowSolutionPointers",Character,2,idata,solname,ierr)
    call cgns_error(pname,fname_nodal,"cg_array_write_f",ierr)
    !
    ! Identify the flow solution within the CGNS file as time dependent
    !
    call cg_simulation_type_write_f(ifile_n,ibase_n,TimeAccurate,ierr)
    call cgns_error(pname,fname_nodal,"cg_simulation_type_write_f",ierr)
    !
    ! Close the CGNS FILE
    !
    call cgp_close_f(ifile_n,ierr)
    call io_error(pname,trim(fname_nodal),2,ierr)
    !
  end if
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
      call cgns_error(pname,fname_nodal,"cgp_queue_set_f",ierr)
    end if
    !
    ! Add the data for current variable to the output queue.
    !
    call cgp_field_write_data_f(ifile_n,ibase_n,izone_n,isol_n,isvar_n(n), &
                                my_pts_beg,my_pts_end,var,ierr)
    call cgns_error(pname,fname_nodal,"cgp_field_write_data_f",ierr, &
                    "Write data: Solution variable",n)
    !
    ! Finally, flush the queue containing
    ! the current variable to the CGNS file.
    !
    if (cgns_use_queue) then
      call cgp_queue_flush_f(ierr)
      call cgns_error(pname,fname_nodal,"cgp_queue_flush_f",ierr)
    end if
    !
  end subroutine cgns_write_field_parallel
  !
  !-----------------------------------------------------------------------------
  !
end subroutine write_nodal_parallel_cgns
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
  real(wp), dimension(1:maxSP) :: lmat
  real(wp), dimension(1:nr) :: xyz_fp
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
      lmat(1:npl) = interp(left_geom,left_order)%toFace(face_order)%mat(1:npl,kl)
      do l = 1,nr
        xyz_fp(l) = dot( lmat(1:npl) , xyz(l,l1:l2) )
      end do
      !
      write (151,1) kl,xyz_fp(1:nr)
      !
      lmat(1:npr) = interp(right_geom,right_order)%toFace(face_order)%mat(1:npr,kr)
      do l = 1,nr
        xyz_fp(l) = dot( lmat(1:npr) , xyz(l,r1:r2) )
      end do
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
subroutine get_initial_data_for_cgns(my_cells_beg,my_cells_end,idim1,ipelem)
  !
  !.. Use Statements ..
  use geovar,       only : ncell,nodes_of_cell_ptr
  use geovar,       only : n_global_cell,global_nofc_ptr
  use geovar,       only : pst_geom
 !use geovar,       only : geom_conn
  use geovar,       only : cell_geom,global_cell_geom
  use geovar,       only : cell_order,global_cell_order
  use parallel_mod, only : cell_map
  use order_mod,    only : geom_solpts
  !
  !.. Formal Arguments ..
  integer,                              intent(inout) :: my_cells_beg
  integer,                              intent(inout) :: my_cells_end
  integer,              dimension(:),   intent(inout) :: idim1
  integer, allocatable, dimension(:,:), intent(inout) :: ipelem
  !
  !.. Local Scalars ..
  integer :: n,nc,np,ierr,p1,icpu
  integer :: gmin,gmax,this_npts,this_subcells
  integer :: this_cell,this_geom,this_order
  integer :: n_total_output_cells
  integer :: n_total_output_pts
  integer :: n_my_output_cells
  character(len=200) :: error_message
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_initial_data_for_cgns"
  !
continue
  !
  ! Initialize n_total output_cells and n_total_output_pts which give:
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
  ! Save n_total_output_pts and n_total_output_cells to idim1
  !
  idim1 = [n_total_output_pts,n_total_output_cells,0]
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
    call alloc_error(pname,"ipelem",1,ierr,error_message)
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
    this_geom = cell_geom(this_cell)
    !
    ! Polynomial order for this grid cell
    !
    this_order = cell_order(this_cell)
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
pure function grad_of_selected_pv(ndr,npv,uu,dcvdx) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : rgasref_nd,gam
  !
  !.. Formal Arguments ..
  integer,                    intent(in) :: ndr ! index for derivative direction
  integer,                    intent(in) :: npv ! index for primitive variable
  real(wp), dimension(:,:),   intent(in) :: uu  ! conservative solution
  real(wp), dimension(:,:,:), intent(in) :: dcvdx ! conservative variable gradients
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
  use geovar,  only : ncell,n_global_cell
  use geovar,  only : pts_in_cell_ptr,global_pinc_ptr
  use geovar,  only : cell_geom,global_cell_geom
  use geovar,  only : cell_order,global_cell_order
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
      this_geom = cell_geom(nc)
      this_order = cell_order(nc)
      !
      n1 = pts_in_cell_ptr(nc)+1
      n2 = pts_in_cell_ptr(nc+1)
      np = n2-n1+1
      !
      if (associated(outerp(this_geom,this_order)%output)) then
        return_value(n1:n2) = matmul( outerp(this_geom,this_order)%output , &
                                      array(n1:n2) )
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
        return_value(n1:n2) = matmul( outerp(this_geom,this_order)%output , &
                                      array(n1:n2) )
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
  use geovar,  only : ncell,n_global_cell
  use geovar,  only : pts_in_cell_ptr,global_pinc_ptr
  use geovar,  only : cell_geom,global_cell_geom
  use geovar,  only : cell_order,global_cell_order
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
      this_geom = cell_geom(nc)
      this_order = cell_order(nc)
      !
      n1 = pts_in_cell_ptr(nc)+1
      n2 = pts_in_cell_ptr(nc+1)
      np = n2-n1+1
      !
      if (associated(outerp(this_geom,this_order)%output)) then
        do m = 1,nvar
          return_value(m,n1:n2) = matmul( outerp(this_geom,this_order)%output , &
                                          array(m,n1:n2) )
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
          return_value(m,n1:n2) = matmul( outerp(this_geom,this_order)%output , &
                                          array(m,n1:n2) )
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
  use geovar, only : nr,ncell,pts_in_cell_ptr
  use geovar, only : cell_geom,cell_order
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
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    np = n2-n1+1
    !
    this_geom = cell_geom(n)
    this_order = cell_order(n)
    !
    do k = n1,n2
      vel(1:nr,k-n1+1) = uu(nmb:nme,k)/uu(1,k)
    end do
    !
    do k = n1,n2
      !
      dudr = dot(deriv(this_geom,this_order)%d(1)%mat(1:np,k-n1+1),vel(1,1:np))
      duds = dot(deriv(this_geom,this_order)%d(2)%mat(1:np,k-n1+1),vel(1,1:np))
      dvdr = dot(deriv(this_geom,this_order)%d(1)%mat(1:np,k-n1+1),vel(2,1:np))
      dvds = dot(deriv(this_geom,this_order)%d(2)%mat(1:np,k-n1+1),vel(2,1:np))
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
  use geovar, only : ncell,pts_in_cell_ptr
  use geovar, only : cell_geom,cell_order
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
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    np = n2-n1+1
    !
    this_geom = cell_geom(n)
    this_order = cell_order(n)
    !
    do k = n1,n2
      !
      dudr = dot(deriv(this_geom,this_order)%d(1)%mat(1:np,k-n1+1),vv(2,n1:n2))
      duds = dot(deriv(this_geom,this_order)%d(2)%mat(1:np,k-n1+1),vv(2,n1:n2))
      dvdr = dot(deriv(this_geom,this_order)%d(1)%mat(1:np,k-n1+1),vv(3,n1:n2))
      dvds = dot(deriv(this_geom,this_order)%d(2)%mat(1:np,k-n1+1),vv(3,n1:n2))
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
include "Functions/last.f90"
include "Functions/vortex.f90"
include "Functions/bc_priority.f90"
include "Functions/viscosity.f90"
include "Functions/entropy.f90"
include "Functions/pressure.f90"
include "Functions/usp2v.f90"
include "Functions/vsp2u.f90"
include "Functions/mach_number.f90"
include "Functions/cell_vol.f90"
include "Functions/bc_string_to_integer.f90"
!
!###############################################################################
!
end module io_mod
!
!###############################################################################
!
