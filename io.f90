module io_mod
  !
  !.. Use Statements ..
  use module_kind_types
  use ovar, only : grid_format
  use ovar, only : gridfile
  use ovar, only : output_dir
  use ovar, only : plot3d_cutfile
  use ovar, only : plot3d_bcsfile
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
 !public :: read_input_file_test
  public :: read_grid_file
  public :: read_restart_file
  public :: write_restart_file
  public :: write_tecplot_file
  public :: write_final_error
  public :: write_final_advec_error
  public :: write_final_vortex_error
  public :: write_entropy_error
  public :: write_cgns_solution
  public :: write_results_file
  public :: write_advec_results_file
  public :: write_vortex_results_file
  public :: write_tecplot_partitions
  public :: create_output_arrays
  !
  ! ###########################
  ! ###  Module Parameters  ###
  ! ###########################
  !
  integer, parameter :: iorst = 10
  integer, parameter :: iotec = 11
  integer, parameter :: ioxyd = 12
  integer, parameter :: iocnd = 13
  integer, parameter :: iobcd = 14
  integer, parameter :: iores = 15
  integer, parameter :: iodir = 16
  !
  integer, parameter :: patran_neutral = 1
  integer, parameter :: gambit_neutral = 2
  integer, parameter :: gmsh_file = 3
  !
  ! Plot3D grid formats
  !
  integer, parameter :: plot2d_file = 4
  integer, parameter :: plot3d_file = 5
  integer, parameter :: plot3d_formats(1:2) = [ plot2d_file, plot3d_file ]
  !
  character(len=11), parameter :: clas_rk(1:4) = ["RK1-Classic","RK2-Classic", &
                                                  "RK3-Classic","RK4-Classic"]
  character(len=11), parameter :: char_rk(1:5) = ["           ", &
                                                  "RK2-SSP    ","RK3-SSP    ", &
                                                  "RK4-SSP    ","RK4-CK     "]
  character(len=3), parameter :: char_init(-1:2) = ["iRa","iSP","iGa","iLo"]
  !
#ifdef DCGNS_ADF
  integer, parameter :: file_type = CG_FILE_ADF
#elif DCGNS_ADF2
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
  !.. Input Variables ..
 !integer,            public, save :: grid_format = gambit_neutral
 !character(len=150), public, save :: gridfile = "grid_files/meilin/quad.40x40.2X.neu"
 !character(len=150), public, save :: output_dir = "."
 !character(len=150), public, save :: plot3d_cutfile = ""
 !character(len=150), public, save :: plot3d_bcsfile = ""
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
  real(vhp),         save, allocatable, dimension(:) :: soltime
  real(vhp),         save, allocatable, dimension(:) :: modtime
  character(len=32), save, allocatable, dimension(:) :: solname
  character(len=32), save, allocatable, dimension(:) :: modname
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
continue
  !
  open (inpt,file="recon3.in",status="old",action="read", &
             form="formatted",iostat=ierr)
  call io_error("read_input_file","recon3.in",1,ierr)
  !
  read (inpt,*,iostat=ierr,err=100,end=100) text
  !
  close (inpt,iostat=ierr)
  call io_error("read_input_file","recon3.in",2,ierr)
  !
  text = trim(adjustl(text))
  text = uppercase( text )
  !
  if (text(1:6) == "&INPUT") then
    call read_nml_input_file
  else
    call read_old_input_file
  end if
  !
  return
  !
  ! Read error continuation statements
  !
  100 write (iout,10) ierr ; call stop_recon(stop_mpi,"read_input_file")
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
  !... kCE: number of points on each Curved Edge
  !...(linear edge = 2; parabolic edge = 3)
  !
  !... itestcase = 3 for Shu's vortex problem;
  !...           = 0 or 1 for bump problem
  !...           = 2 Ringleb
  !... iopSP = 1 Gauss points for SP; = 2 Lobatto points for SP;
  !...       = 3 equi distant including two boundaries for SP
  !... iopCE = 1 Gauss points for SP; = 2 Lobatto points for SP;
  !...       = 3 equi distant including two boundaries for SP
  !... we only use Lobatto and equidistant points here
  !... irestart = 1 for restart run; = 0 for fresh start
  !... flux_divergence = 1 calculate derivative (no interaction) using
  !...               Lagrange polynomials
  !...          = 2 calculate derivative (no interaction) using
  !...               the chain rule for U
  !...          = 3 calculate derivative (no interaction) using
  !...               the chain rule for Uja
  !... isupersonic = 1 for supersonic flow (bump problem),
  !...             = 0 for subsonic
  !...
  !... it1 was iopSP and it2, kSP1d (but we do not need these anymore)
  !
  use geovar, only : p_order,n_order,n_sp1de,n_sp2dt,n_sp2dq,n_fp1de
  use geovar, only : o_order,all_orders
  use geovar, only : q_order,n_qp1de,n_qp2dt,n_qp2dq,using_quadrature
  use ovar,   only : itestcase,num_stages,irestart, &
                     flux_divergence,iop_correction, &
                     num_ts,cfl,dt_in,gam,output_solution, &
                     output_interval,restart_interval,dflux_opt, &
                     results_interval,invs_flux_method,visc_flux_method, &
                     iRKtype,CK4_RK,TVD2_RK,TVD3_RK,TVD4_RK, &
                     Final_Time,Total_Periods,all_points, &
                     iop_solpts,iop_flxpts,iop_tripts,iop_qdtpts,iop_outpts, &
                     Filter_Option,Filter_Order, &
                     Filter_Cutoff,Filter_Alpha,Filter_etac, &
                     rhoref,tref,pref,machref,aoaref, &
                     vxref,vyref,aref,uref,rgasref,reyref, &
                     metis_option,init_loc,error_loc, &
                     Period_Timesteps,Period_Time,iter_out_interval, &
                     suth_muref,suth_Tref,suth_Sref,Pr,Prt
  use ovar,   only : Grid_Length_X,Grid_Length_Y
  use ovar,   only : bound_error,error_bound_limit
  use ovar,   only : vortex_bigR,vortex_bigGam
  use ovar,   only : Timestep_Type
  use ovar,   only : output_plot3d_to_tecplot
  use ovar,   only : walls_are_exact
  use ovar,   only : input
  !
  !.. Local Scalars ..
  integer   :: it1,it2,ierr,ival
  integer   :: Period_Timesteps_X,Period_Timesteps_Y
  real(vbp) :: Period_Time_X,Period_Time_Y
  real(vbp) :: real_Total_Periods
  real(vbp) :: aoa_radians,pres_check
  !
  !.. Local Arrays ..
  character(len=80), dimension(12) :: input_text
  !
 !!.. Input Namelist ..
 !! GRID INFORMATION
 !namelist / input / grid_format, gridfile, plot3d_cutfile, plot3d_bcsfile
 !! BASE FR OPTIONS
 !namelist / input / n_order, iop_solpts, iop_flxpts, iop_tripts, iop_outpts, &
 !                   flux_divergence, iop_correction
 !! SOLVER OPTIONS
 !namelist / input / iRKtype, num_stages, invs_flux_method, visc_flux_method
 !! TIME STEPPING INFORMATION
 !namelist / input / num_ts, Final_Time, cfl, dt_in, Timestep_Type
 !! OVER-INTEGRATION OPTIONS
 !namelist / input / p_order, q_order, iop_qdtpts, dflux_opt
 !! FORCE GLOBAL ORDER AND SOLUTION POINTS
 !namelist / input / all_orders, all_points
 !! PARALLEL OPTIONS
 !namelist / input / metis_option
 !! INITIALIZATION / TEST CASE
 !namelist / input / itestcase
 !! RESTART OPTIONS
 !namelist / input / irestart, restart_interval
 !! OUTPUT OPTIONS
 !namelist / input / output_dir, iter_out_interval, output_solution, &
 !                   output_interval, results_interval, o_order
 !! NONDIMENSIONAL REFERENCE CONDITIONS
 !namelist / input / machref, reyref, gam
 !! REFERENCE CONDITIONS
 !namelist / input / rhoref, tref, pref, aoaref, rgasref
 !! LAMINAR AND TURBULENT PRANDTL NUMBERS
 !namelist / input / Pr, Prt
 !! SUTHERLANDS LAW COEFFICIENTS
 !namelist / input / suth_muref, suth_Tref, suth_Sref
 !! FILTER OPTIONS
 !namelist / input / Filter_Option, Filter_etac, Filter_Order, Filter_Alpha
 !! BOUNDING ERROR OPTIONS
 !namelist / input / bound_error, error_bound_limit
 !! INITIALIZATION/ERROR LOCATIONS
 !namelist / input / init_loc, error_loc
 !! VORTEX CONSTANTS
 !namelist / input / vortex_bigR, vortex_bigGam
 !! OPTION TO OUTPUT PLOT3D GRID TO TECPLOT IN UNSTRUCTURED FORMAT
 !namelist / input / output_plot3d_to_tecplot
 !! LOGICAL IDENTIFYING HOW WALL SOLUTIONS ARE TREATED
 !namelist / input / walls_are_exact
  !
continue
  !
  ! Open the input file
  !
  open (inpt,file="recon3.in",status="old",action="read", &
             form="formatted",iostat=ierr)
  call io_error("read_nml_input_file","recon3.in",1,ierr)
  !
  ! Read in the input file using the input namelist
  !
  read (inpt,nml=input)
  !
  ! Close the input file now that we are finished reading it
  !
  close (inpt,iostat=ierr)
  call io_error("read_nml_input_file","recon3.in",2,ierr)
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
  if (all_orders > 0) then
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
    call stop_recon(stop_mpi,"read_nml_input_file")
  end if
  !
  ! Adjust the solution order and compute the number
  ! of solution points for the different cell shapes
  !
  n_order = max(p_order,n_order)
  n_sp1de = n_order + 1
  n_sp2dt = n_sp1de*(n_sp1de+1)/2
  n_sp2dq = n_sp1de*n_sp1de
  !
  ! Now that we have the name of the grid file and the location
  ! of the flux and solution points, read in the grid file
  !
  call read_grid_file
  !
  ! Make sure the quadrature order is at least the same as the
  ! solution order and other quadrature variables are valid
  !
  q_order = max(n_order,q_order)
  !
  n_qp1de = q_order + 1
  n_qp2dt = n_qp1de*(n_qp1de+1)/2
  n_qp2dq = n_qp1de*n_qp1de
  !
  using_quadrature = .not. (q_order==n_order .and. iop_solpts==iop_qdtpts)
  !
  if (using_quadrature .and. dflux_opt == project_derivatives) then
    n_fp1de = n_qp1de
  else
    n_fp1de = n_sp1de
  end if
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
  if (all(abs(init_loc) /= [1,2])) init_loc = 0
  if (abs(error_loc) /= 1) error_loc = 0
  !
  ! Compute the filter cutoff level
  !
  Filter_Cutoff = nint( max(zero, min(one,Filter_etac) ) * real(n_order,vbp) )
  !
  ! Make sure the output order is greater than or
  ! equal to the solution order if it is positive
  !
  if (o_order > 0) o_order = max(o_order,n_order+0)
  !
  ! Make any necessary adjustments to the input parameters
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
      write (iout,99) pref,pres_check
      call stop_recon(stop_mpi,"read_nml_input_file")
    end if
  end if
  !
  ! Compute some of the remaining reference conditions
  !
  aref = sqrt( gam * pref / rhoref )
  uref = machref * aref
  vxref = uref * cos(aoaref * degrees2radians)
  vyref = uref * sin(aoaref * degrees2radians)
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
      aoaref = 90.0_vbp
      vxref = zero
      vyref = zero
      itestcase = abs(itestcase)
    else if (itestcase == -Vortex_Transport) then
      aoaref = 90.0_vbp
      vxref = zero
      vyref = one
      itestcase = Shu_Vortex
    else
      aoaref = 45.0_vbp
      vxref = one
      vyref = one
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
  ! Compute periodic information for the transport problems
  !
  if (any(itestcase == Transport_Problems)) then
    !
    cfl = one
    !
    if (abs(vxref) < eps3) then
      !
      Period_Time_Y = Grid_Length_Y * reciprocal(vyref/aref)
      Period_Time_Y = merge( Period_Time_Y , Grid_Length_Y , Period_Time_Y > eps3 )
      Period_Timesteps_Y = nint(Period_Time_Y/dt_in)
      Period_Timesteps = Period_Timesteps_Y
      !
     !if (mypnum == 0) write (*,*) Period_Time_Y,Period_Timesteps_Y,Period_Timesteps
      !
    else if (abs(vyref) < eps3) then
      !
      Period_Time_X = Grid_Length_X * reciprocal(vxref/aref)
      Period_Time_X = merge( Period_Time_X , Grid_Length_X , Period_Time_X > eps3 )
      Period_Timesteps_X = nint(Period_Time_X/dt_in)
      Period_Timesteps = Period_Timesteps_X
      !
     !if (mypnum == 0) write (*,*) Period_Time_Y,Period_Timesteps_Y,Period_Timesteps
      !
    else
      !
      Period_Time_X = Grid_Length_X / (vxref/aref)
      Period_Time_Y = Grid_Length_Y / (vyref/aref)
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
    Period_Time = dt_in*real(Period_Timesteps,vbp)
    !
    ! If Final_Time is negative, it is giving the total
    ! number of periods that the simulation will run
    !
    if (Final_Time < zero) then
      Total_Periods = nint( abs(Final_Time) )
      Final_Time = Period_Time * real(Total_Periods,vbp)
    else
      real_Total_Periods = Final_Time*reciprocal(Period_Time)
      Total_Periods = nint(real_Total_Periods)
      if (abs(real_Total_Periods-real(Total_Periods,vbp)) > ten*eps6) then
        if (mypnum == 0) write (iout,2) real_Total_Periods
      end if
    end if
   !if (mypnum == 0) write (*,*) Period_Time,Total_Periods,Final_Time
    !
    num_ts = nint( Final_Time / dt_in )
    !
  else
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
  ! Now write out a copy of the input file to standard output
  !
  if (mypnum == 0) then
    !
    write (iout,5)
    write (iout,nml=input)
   !write (iout,1) trim(input_text(1))
   !write (iout,1) trim(output_dir)
   !write (iout,1) trim(input_text(2))
   !write (iout,50) grid_format, trim(gridfile)
   !write (iout,1) trim(input_text(3))
   !write (iout,51) p_order, n_order, iop_solpts, iop_flxpts, iop_tripts, metis_option
   !write (iout,1) trim(input_text(4))
   !write (iout,52) flux_divergence, iop_correction, &
   !                q_order, iop_qdtpts, dflux_opt
   !write (iout,1) trim(input_text(5))
   !write (iout,53) itestcase, iRKtype, num_stages, irestart, init_loc, error_loc
   !write (iout,1) trim(input_text(6))
   !write (iout,54) num_ts, Final_Time, cfl, dt_in, invs_flux_method, iter_out_interval
   !write (iout,1) trim(input_text(7))
   !write (iout,55) Filter_Option, Filter_etac, Filter_Cutoff, &
   !                Filter_Order, Filter_Alpha
   !write (iout,1) trim(input_text(8))
   !write (iout,56) output_solution, output_interval, &
   !                restart_interval, results_interval
   !write (iout,1) trim(input_text(11))
   !write (iout,57) bound_error,error_bound_limit
   !write (iout,1) trim(input_text(12))
   !write (iout,58) o_order
    write (iout,6)
    !
  end if
  !
  1  format (a)
  5  format (/,"#######################################",/, &
               "#####      INPUT PARAMETERS       #####",/, &
               "#######################################",/)
  6  format (/,"#######################################",/, &
               "#####   END OF INPUT PARAMETERS   #####",/, &
               "#######################################",/)
  50 format (i8,8x,a)
  51 format (3x,i2,7x,i2,8x,i2,10x,i2,10x,i2,11x,i2)
  52 format (7x,i2,15x,i2,2x,3(10x,i2))
  53 format (3x,i2,7x,i2,8x,i2,8x,i2,9x,i2,8x,i2)
  54 format (i10,2x,f10.3,3x,f8.4,3x,f9.6,9x,i0,8x,i10)
  55 format (i8,2x,f7.3,3x,i8,3x,i8,3x,f7.3)
  56 format (1x,i8,9x,i8,9x,i8,9x,i8)
  57 format (6x,l1,12x,f8.3)
  58 format (i7)
  78 format (" Dimensional Free-stream conditions",/, &
             " Density         = ",es14.6," kg/m^3",/, &
             " Pressure        = ",es14.6," N/m^2",/, &
             " Temperature     = ",es14.6," K",/, &
             " Mach Number     = ",es14.6,/, &
             " Speed of Sound  = ",es14.6," m/s",/, &
             " Velocity        = ",es14.6," m/s",/, &
             " Angle of Attack = ",es14.6," degrees"/, &
             " Gamma           = ",es14.6,/, &
             " Gas Constant    = ",es14.6," J/(kg*K)",/)
  99 format (/, &
     "****************************** ERROR ******************************",/, &
     " Free-stream conditions in the input file are not consistent! ",/, &
     " Reference pressure given    : ",es14.6,/, &
     " Reference pressure computed : ",es14.6,/, &
     " Please correct the input file or set one of density, pressure, ",/, &
     " and temperature to a negative value so it is explicitly computed. ",/, &
     "****************************** ERROR ******************************",/)
  2 format (/," ********************* WARNING!!! *********************",/, &
              "  IT APPEARS THE FINAL TIME DOES NOT CORRESPOND TO THE",/, &
              "  EXACT END OF A PERIOD FOR THE VORTEX PROBLEM. THE",/, &
              "  FINAL TIME CORRESPONDS TO",f12.5," TOTAL PERIODS.",/, &
              " ********************* WARNING!!! *********************",/)
  101 format (/, &
     "***************************** ERROR *****************************",/, &
     "   THE REQUESTED LOCATIONS FOR SOLUTION POINTS AND FLUX POINTS",/, &
     "   DO NOT MATCH. CURRENTLY, THE PROCEDURES USED TO CREATE THE",/, &
     "   CORRECTION MATRICES ASSUMES THESE TWO POINTS LOCATIONS ARE",/, &
     "   THE SAME. THESE PROCEDURES WILL NEED TO BE UPDATED IN ORDER",/, &
     "   TO ALLOW FOR THIS CAPABILITY.",/, &
     "***************************** ERROR *****************************",/)
  !
end subroutine read_nml_input_file
!
!###############################################################################
!
subroutine read_old_input_file
  !
  !... kCE: number of points on each Curved Edge
  !...(linear edge = 2; parabolic edge = 3)
  !
  !... itestcase = 3 for Shu's vortex problem;
  !...           = 0 or 1 for bump problem
  !...           = 2 Ringleb
  !... iopSP = 1 Gauss points for SP; = 2 Lobatto points for SP;
  !...       = 3 equi distant including two boundaries for SP
  !... iopCE = 1 Gauss points for SP; = 2 Lobatto points for SP;
  !...       = 3 equi distant including two boundaries for SP
  !... we only use Lobatto and equidistant points here
  !... irestart = 1 for restart run; = 0 for fresh start
  !... flux_divergence = 1 calculate derivative (no interaction) using
  !...               Lagrange polynomials
  !...          = 2 calculate derivative (no interaction) using
  !...               the chain rule for U
  !...          = 3 calculate derivative (no interaction) using
  !...               the chain rule for Uja
  !... isupersonic = 1 for supersonic flow (bump problem),
  !...             = 0 for subsonic
  !...
  !... it1 was iopSP and it2, kSP1d (but we do not need these anymore)
  !
  use geovar, only : p_order,n_order,n_sp1de,n_sp2dt,n_sp2dq,n_fp1de
  use geovar, only : o_order
  use geovar, only : q_order,n_qp1de,n_qp2dt,n_qp2dq,using_quadrature
  use ovar,   only : itestcase,num_stages,irestart, &
                     flux_divergence,iop_correction, &
                     num_ts,cfl,dt_in,gam,output_solution, &
                     output_interval,restart_interval, &
                     results_interval,invs_flux_method,dflux_opt, &
                     iRKtype,CK4_RK,TVD2_RK,TVD3_RK,TVD4_RK, &
                     Final_Time,Total_Periods, &
                     iop_solpts,iop_flxpts,iop_tripts,iop_qdtpts, &
                     Filter_Option,Filter_Order, &
                     Filter_Cutoff,Filter_Alpha,Filter_etac, &
                     rhoref,tref,pref,machref,aoaref, &
                     vxref,vyref,aref,uref,rgasref, &
                     metis_option,init_loc,error_loc, &
                     Period_Timesteps,Period_Time,iter_out_interval
  use ovar,   only : Grid_Length_X,Grid_Length_Y
  use ovar,   only : bound_error,error_bound_limit
 !use ovar,   only : input
  !
  !.. Local Scalars ..
  integer   :: it1,it2,ierr,ival
  integer   :: Period_Timesteps_X,Period_Timesteps_Y
  real(vbp) :: Period_Time_X,Period_Time_Y
  real(vbp) :: real_Total_Periods
  real(vbp) :: aoa_radians,pres_check
  !
  !.. Local Arrays ..
  character(len=80), dimension(12) :: input_text
  !
  !.. Input Namelist ..
  ! GRID INFORMATION
  namelist / input / grid_format, gridfile
  ! BASE FR OPTIONS
  namelist / input / n_order, iop_solpts, iop_flxpts, iop_tripts, &
                     flux_divergence, iop_correction
  ! SOLVER OPTIONS
  namelist / input / iRKtype, num_stages, invs_flux_method
  ! TIME STEPPING INFORMATION
  namelist / input / num_ts, Final_Time, cfl, dt_in
  ! OVER-INTEGRATION OPTIONS
  namelist / input / p_order, q_order, iop_qdtpts, dflux_opt
  ! PARALLEL OPTIONS
  namelist / input / metis_option
  ! INITIALIZATION / TEST CASE
  namelist / input / itestcase
  ! RESTART OPTIONS
  namelist / input / irestart, restart_interval
  ! OUTPUT OPTIONS
  namelist / input / output_dir, iter_out_interval, output_solution, &
                     output_interval, results_interval, o_order
  ! REFERENCE CONDITIONS
  namelist / input / rhoref, tref, pref, machref, aoaref, gam, rgasref
  ! FILTER OPTIONS
  namelist / input / Filter_Option, Filter_etac, Filter_Order, Filter_Alpha
  ! BOUNDING ERROR OPTIONS
  namelist / input / bound_error, error_bound_limit
  ! INITIALIZATION/ERROR LOCATIONS
  namelist / input / init_loc, error_loc
  !
continue
  !
  ! Open the input file
  !
  open (inpt,file="recon3.in",status="old",action="read", &
             form="formatted",iostat=ierr)
  call io_error("read_old_input_file","recon3.in",1,ierr)
  !
  ! Read the title line of the input file and the output directory
  !
  read (inpt,1) input_text(1)
  read (inpt,*) output_dir
  !
  output_dir = trim(adjustl(output_dir))
  call check_output_dir
  !
  ! Append the name of the solution and error files to output_dir
  !
  fname_nodal = trim(output_dir)//"solution.nodal.cgns"
  fname_modal = trim(output_dir)//"solution.modal.cgns"
  error_file = trim(output_dir)//"error.dat"
  !
  ! Read the name of the grid file and its format
  !
  read (inpt,1) input_text(2)
  read (inpt,*) grid_format, gridfile
  !
  gridfile = trim(adjustl(gridfile))
  !
  ! Read the polynomial order, node locations, and MeTiS option
  !
  read (inpt,1) input_text(3)
  read (inpt,*) p_order, n_order, iop_solpts, iop_flxpts, iop_tripts, metis_option
  !
  if (iop_solpts /= iop_flxpts) then
    if (mypnum == 0) write (iout,101)
    call stop_recon(stop_mpi,"read_old_input_file")
  end if
  !
  n_order = max(p_order,n_order)
  n_sp1de = n_order + 1
  n_sp2dt = n_sp1de*(n_sp1de+1)/2
  n_sp2dq = n_sp1de*n_sp1de
  !
  ! Now that we have the name of the grid file and the location
  ! of the flux and solution points, read in the grid file
  !
  call read_grid_file
  !
  ! Read the option for computing the discontinuous flux derivatives
  ! and the option of what variables to use when computing the flux corrections
  !
  read (inpt,1) input_text(4)
  read (inpt,*) flux_divergence, iop_correction, &
                q_order, iop_qdtpts, dflux_opt
  !
  ! Make sure the quadrature order is at least the same as the solution order
  !
  q_order = max(n_order,q_order)
  !
  n_qp1de = q_order + 1
  n_qp2dt = n_qp1de*(n_qp1de+1)/2
  n_qp2dq = n_qp1de*n_qp1de
  !
  using_quadrature = .not. (q_order==n_order .and. iop_solpts==iop_qdtpts)
  !
  if (using_quadrature .and. dflux_opt == project_derivatives) then
    n_fp1de = n_qp1de
  else
    n_fp1de = n_sp1de
  end if
  !
  ! Read the testcase, Runge-Kutta options, restart, etc
  !
  read (inpt,1) input_text(5)
  read (inpt,*) itestcase, iRKtype, num_stages, irestart, init_loc, error_loc
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
  if (all(abs(init_loc) /= [1,2])) init_loc = 0
  if (abs(error_loc) /= 1) error_loc = 0
  !
  ! Read the time step information and the method for computing the interface fluxes
  !
  read (inpt,1) input_text(6)
  read (inpt,*) num_ts, Final_Time, cfl, dt_in, invs_flux_method, iter_out_interval
  !
  ! Read the filter options
  !
  read (inpt,1) input_text(7)
  read (inpt,*) Filter_Option, Filter_etac, Filter_Order, Filter_Alpha
  !
  Filter_Cutoff = nint( max(zero, min(one,Filter_etac) ) * real(n_order,vbp) )
  !
  ! Read the output options and intervals
  !
  read (inpt,1) input_text(8)
  read (inpt,*) output_solution, output_interval, &
                restart_interval, results_interval
  !
  ! Read the reference conditions for the flow variables
  !
  read (inpt,1) input_text(9)
  read (inpt,1) input_text(10)
  read (inpt,*) rhoref, tref, pref, machref, aoaref, gam, rgasref
  !
  ! Read the information on bounding the region to compute error
  !
  read (inpt,1) input_text(11)
  read (inpt,*) ival,error_bound_limit
  bound_error = (ival /= 0)
  !
  ! Read the output order
  !
  read (inpt,1) input_text(12)
  if (input_text(12)(1:1) /= "!") then
    read (inpt,*) o_order
    if (o_order > 0) o_order = max(o_order,n_order+0)
  else
    input_text(12) = " Output_Order"
    o_order = -n_order
  end if
  !
  ! Close the input file now that we are finished reading it
  !
  close (inpt,iostat=ierr)
  call io_error("read_old_input_file","recon3.in",2,ierr)
  !
  ! Make any necessary adjustments to the input parameters
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
      write (iout,99) pref,pres_check
      call stop_recon(stop_mpi,"read_old_input_file")
    end if
  end if
  !
  ! Compute some of the remaining reference conditions
  !
  aref = sqrt( gam * pref / rhoref )
  uref = machref * aref
  vxref = uref * cos(aoaref * degrees2radians)
  vyref = uref * sin(aoaref * degrees2radians)
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
      aoaref = 90.0_vbp
      vxref = zero
      vyref = zero
      itestcase = abs(itestcase)
    else if (itestcase == -Vortex_Transport) then
      aoaref = 90.0_vbp
      vxref = zero
      vyref = one
      itestcase = Shu_Vortex
    else
      aoaref = 45.0_vbp
      vxref = one
      vyref = one
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
  ! Compute periodic information for the transport problems
  !
  if (any(itestcase == Transport_Problems)) then
    !
    cfl = one
    !
    if (abs(vxref) < eps3) then
      !
      Period_Time_Y = Grid_Length_Y * reciprocal(vyref/aref)
      Period_Time_Y = merge( Period_Time_Y , Grid_Length_Y , Period_Time_Y > eps3 )
      Period_Timesteps_Y = nint(Period_Time_Y/dt_in)
      Period_Timesteps = Period_Timesteps_Y
      !
     !if (mypnum == 0) write (*,*) Period_Time_Y,Period_Timesteps_Y,Period_Timesteps
      !
    else if (abs(vyref) < eps3) then
      !
      Period_Time_X = Grid_Length_X * reciprocal(vxref/aref)
      Period_Time_X = merge( Period_Time_X , Grid_Length_X , Period_Time_X > eps3 )
      Period_Timesteps_X = nint(Period_Time_X/dt_in)
      Period_Timesteps = Period_Timesteps_X
      !
     !if (mypnum == 0) write (*,*) Period_Time_Y,Period_Timesteps_Y,Period_Timesteps
      !
    else
      !
      Period_Time_X = Grid_Length_X / (vxref/aref)
      Period_Time_Y = Grid_Length_Y / (vyref/aref)
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
    Period_Time = dt_in*real(Period_Timesteps,vbp)
    !
    ! If Final_Time is negative, it is giving the total
    ! number of periods that the simulation will run
    !
    if (Final_Time < zero) then
      Total_Periods = nint( abs(Final_Time) )
      Final_Time = Period_Time * real(Total_Periods,vbp)
    else
      real_Total_Periods = Final_Time*reciprocal(Period_Time)
      Total_Periods = nint(real_Total_Periods)
      if (abs(real_Total_Periods-real(Total_Periods,vbp)) > ten*eps6) then
        if (mypnum == 0) write (iout,2) real_Total_Periods
      end if
    end if
   !if (mypnum == 0) write (*,*) Period_Time,Total_Periods,Final_Time
    !
    num_ts = nint( Final_Time / dt_in )
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
  end if
  !
  ! Now write out a copy of the input file to standard output
  !
  if (mypnum == 0) then
    !
    write (iout,*)
    write (iout,1) trim(input_text(1))
    write (iout,1) trim(output_dir)
    write (iout,1) trim(input_text(2))
    write (iout,50) grid_format, trim(gridfile)
    write (iout,1) trim(input_text(3))
    write (iout,51) p_order, n_order, iop_solpts, iop_flxpts, iop_tripts, metis_option
    write (iout,1) trim(input_text(4))
    write (iout,52) flux_divergence, iop_correction, &
                    q_order, iop_qdtpts, dflux_opt
    write (iout,1) trim(input_text(5))
    write (iout,53) itestcase, iRKtype, num_stages, irestart, init_loc, error_loc
    write (iout,1) trim(input_text(6))
    write (iout,54) num_ts, Final_Time, cfl, dt_in, &
                    invs_flux_method, iter_out_interval
    write (iout,1) trim(input_text(7))
    write (iout,55) Filter_Option, Filter_etac, Filter_Cutoff, &
                    Filter_Order, Filter_Alpha
    write (iout,1) trim(input_text(8))
    write (iout,56) output_solution, output_interval, &
                    restart_interval, results_interval
    write (iout,1) trim(input_text(11))
    write (iout,57) bound_error,error_bound_limit
    write (iout,1) trim(input_text(12))
    write (iout,58) o_order
    write (iout,*)
    !
  end if
  !
  1  format (a)
  50 format (i8,8x,a)
  51 format (3x,i2,7x,i2,8x,i2,10x,i2,10x,i2,11x,i2)
  52 format (7x,i2,15x,i2,2x,3(10x,i2))
  53 format (3x,i2,7x,i2,8x,i2,8x,i2,9x,i2,8x,i2)
  54 format (i10,2x,f10.3,3x,f8.4,3x,f9.6,9x,i0,8x,i10)
  55 format (i8,2x,f7.3,3x,i8,3x,i8,3x,f7.3)
  56 format (1x,i8,9x,i8,9x,i8,9x,i8)
  57 format (6x,l1,12x,f8.3)
  58 format (i7)
  78 format (" Dimensional Free-stream conditions",/, &
             " Density         = ",es14.6," kg/m^3",/, &
             " Pressure        = ",es14.6," N/m^2",/, &
             " Temperature     = ",es14.6," K",/, &
             " Mach Number     = ",es14.6,/, &
             " Speed of Sound  = ",es14.6," m/s",/, &
             " Velocity        = ",es14.6," m/s",/, &
             " Angle of Attack = ",es14.6," degrees"/, &
             " Gamma           = ",es14.6,/, &
             " Gas Constant    = ",es14.6," J/(kg*K)",/)
  99 format (/, &
     "****************************** ERROR ******************************",/, &
     " Free-stream conditions in the input file are not consistent! ",/, &
     " Reference pressure given    : ",es14.6,/, &
     " Reference pressure computed : ",es14.6,/, &
     " Please correct the input file or set one of density, pressure, ",/, &
     " and temperature to a negative value so it is explicitly computed. ",/, &
     "****************************** ERROR ******************************",/)
  2 format (/," ********************* WARNING!!! *********************",/, &
              "  IT APPEARS THE FINAL TIME DOES NOT CORRESPOND TO THE",/, &
              "  EXACT END OF A PERIOD FOR THE VORTEX PROBLEM. THE",/, &
              "  FINAL TIME CORRESPONDS TO",f12.5," TOTAL PERIODS.",/, &
              " ********************* WARNING!!! *********************",/)
  101 format (/, &
     "***************************** ERROR *****************************",/, &
     "   THE REQUESTED LOCATIONS FOR SOLUTION POINTS AND FLUX POINTS",/, &
     "   DO NOT MATCH. CURRENTLY, THE PROCEDURES USED TO CREATE THE",/, &
     "   CORRECTION MATRICES ASSUMES THESE TWO POINTS LOCATIONS ARE",/, &
     "   THE SAME. THESE PROCEDURES WILL NEED TO BE UPDATED IN ORDER",/, &
     "   TO ALLOW FOR THIS CAPABILITY.",/, &
     "***************************** ERROR *****************************",/)
  !
end subroutine read_old_input_file
!
!###############################################################################
!
subroutine check_output_dir
  !
  integer :: clen,ierr
  character(len=200) :: test_file
  character(len=200) :: test_dir
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
  open (iodir,file=trim(test_file),status="replace",action="write",iostat=ierr)
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
    open (iodir,file=trim(test_file),status="replace",action="write",iostat=ierr)
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
    open (iodir,file=trim(test_file),status="replace",action="write",iostat=ierr)
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
  call stop_recon(stop_mpi)
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
  use ovar, only : muref, muref_nd, reyref
  use ovar, only : pv_in, rhoref, aref, machref, tref
  use ovar, only : suth_muref, suth_c1
  use ovar, only : suth_Tref, suth_Sref
  !
  use generic_mod, only : viscosity_pv_sp
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
    reyref = rhoref * machref * aref / muref
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
  use geovar, only : nq => neqns
  use ovar,   only : rhoref,aref,tref,gam
  use ovar,   only : vxref,vxref_nd
  use ovar,   only : vyref,vyref_nd
  use ovar,   only : pref,pref_nd
  use ovar,   only : rgasref,rgasref_nd
  use ovar,   only : cpref,cpref_nd
  use ovar,   only : cv_in,pv_in
  !
  !.. Local Scalars ..
  integer  :: ierr
  real(wp) :: rhoa,rhoa2,tda2
  !
continue
  !
  rhoa = rhoref * aref
  rhoa2 = rhoa * aref
  tda2 = tref / aref / aref
  !
  ! Allocate the free stream inflow arrays
  !
  allocate ( cv_in(1:nq) , source=zero , stat=ierr )
  call alloc_error("initialize_reference","cv_in",1,ierr)
  allocate ( pv_in(1:nq) , source=zero , stat=ierr )
  call alloc_error("initialize_reference","pv_in",1,ierr)
  !
  ! Initialize the inflow boundary conditions to the reference conditions
  ! NOTE: This is hard-wired for 2D and needs to be generalized for 3D
  !
  pv_in(1) = rhoref
  pv_in(2) = vxref
  pv_in(3) = vyref
  pv_in(4) = pref
  !
  cv_in(1) = pv_in(1)
  cv_in(2:nq-1) = pv_in(1) * pv_in(2:nq-1)
  cv_in(nq) = pv_in(nq)/(gam-one) + half*pv_in(1)*selfdot(pv_in(2:nq-1))
  !
  ! Nondimensional reference velocities
  vxref_nd = vxref / aref
  vyref_nd = vyref / aref
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
  cv_in(1) = cv_in(1) / rhoref
  cv_in(2:nq-1) = cv_in(2:nq-1) / rhoa
  cv_in(nq) = cv_in(nq) / rhoa2
  !
  pv_in(1) = pv_in(1) / rhoref
  pv_in(2:nq-1) = pv_in(2:nq-1) / aref
  pv_in(nq) = pv_in(nq) / rhoa2
  !
end subroutine initialize_reference
!
!###############################################################################
!
subroutine read_grid_file()
  !
  !.. Use Statements ..
  use module_gmsh,   only : read_gmsh_gridfile
  use module_plot3d, only : read_plot3d_gridfile
  !
  ! Routine that uses the variable grid_format to choose
  ! which routine to use in order to read the grid file.
  !        
continue
  !
  if (grid_format == patran_neutral) then
    call read_patran_neutral_gridfile
  else if (grid_format == gambit_neutral) then
    call read_gambit_neutral_gridfile
  else if (grid_format == gmsh_file) then
    call read_gmsh_gridfile( gridfile )
  else if (any(grid_format == plot3d_formats)) then
    call read_plot3d_gridfile( grid_format , gridfile , plot3d_cutfile , &
                                                        plot3d_bcsfile )
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
  use ovar,   only : iop_solpts,iop_flxpts
  !
  !.. Local Scalars ..
  integer   :: i,j,k,n,nf,n1,n2,i1,nc,pc,pf
  integer   :: icell,inode,iside,ierr,itype
  integer   :: ngrps,nbsets,ndfvl,nentry,nvalues
  !
  character(len=80) :: text
  !
  !.. Local Allocatable Arrays ..
  logical, allocatable, dimension(:) :: bc_done
  !
continue
  !
  ! Open the grid file
  !
  open (iogrd,file=gridfile,status="old",action="read", &
        position="rewind",iostat=ierr)
  call io_error("read_gambit_neutral_gridfile",trim(adjustl(gridfile)),1,ierr)
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
  call alloc_error("read_gambit_neutral_gridfile","nodes_of_cell_ptr",1,ierr)
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
    read (iogrd,*) n,itype,nodes_of_cell_ptr(n+1)
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
  call alloc_error("read_gambit_neutral_gridfile","xyz_nodes",1,ierr)
  allocate ( nodes_of_cell(1:last(nodes_of_cell_ptr)) , source=0 , stat=ierr )
  call alloc_error("read_gambit_neutral_gridfile","nodes_of_cell",1,ierr)
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , stat=ierr )
  call alloc_error("read_gambit_neutral_gridfile","bface",1,ierr)
  allocate ( bc_done(1:nfbnd) , source=true , stat=ierr )
  call alloc_error("read_gambit_neutral_gridfile","bc_done",1,ierr)
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
    if (trim(adjustl(text)) == "SubInflow") then
      itype = bc_sub_inflow
    else if (trim(adjustl(text)) == "SubOutflow") then
      itype = bc_sub_outflow
    else if (trim(adjustl(text)) == "SupInflow") then
      itype = bc_sup_inflow
    else if (trim(adjustl(text)) == "SupOutflow") then
      itype = bc_sup_outflow
    else if (trim(adjustl(text)) == "Characteristic") then
      itype = bc_characteristic
    else if (trim(adjustl(text)) == "Wall") then
      itype = bc_slip_wall
    else if (trim(adjustl(text)) == "FreeStream") then
      itype = bc_freestream
    else if (trim(adjustl(text)) == "Periodic") then
      itype = bc_periodic
      bc_done(n+1:n+nentry) = fals
    else
      itype = bc_unknown
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
  call io_error("read_gambit_neutral_gridfile",trim(adjustl(gridfile)),2,ierr)
  !
  ! Save the index for the ghost cell on the face in bface(10,:)
  !
  bface(10,:) = ncell + [1:nfbnd]
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
  call alloc_error("read_gambit_neutral_gridfile","bc_done",2,ierr)
  !
  ! Make sure that Lobatto nodes are requested if the grid contains triangle
  ! elements since Gauss nodes is not currently set up for triangles.
  !
  do n = 1,ncell
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    if (n2-n1+1 == 3) then
      write (*,3)
      call stop_recon(stop_mpi,"read_gambit_neutral_gridfile")
    end if
  end do
  if (iop_solpts /= Legendre_Gauss_Lobatto .or. &
      iop_flxpts /= Legendre_Gauss_Lobatto) then
    do n = 1,ncell
      n1 = nodes_of_cell_ptr(n)+1
      n2 = nodes_of_cell_ptr(n+1)
      if (n2-n1+1 == 3) then
        write (*,2)
        call stop_recon(stop_mpi,"read_gambit_neutral_gridfile")
      end if
    end do
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
subroutine read_patran_neutral_gridfile
  !
  use geovar
  !
  !  converts a patran neutral file generated by Gridgen using the
  !  NCC Analysis SoftWare selection to geometry output for our purpose
  !
  !.. Local Allocatable Arrays ..
  integer,   allocatable, dimension(:,:) :: edge,c_vtemp
  real(vbp), allocatable, dimension(:)   :: x,y,area_
  !
  !.. Local Scalars ..
  integer   :: ndum,i,k,ierr
  real(vbp) :: rdum
  !
continue
  !
  open(iogrd,file=trim(adjustl(gridfile)),status="old",action="read", &
             position="rewind",iostat=ierr)
  call io_error("read_patran_neutral_gridfile",trim(adjustl(gridfile)),1,ierr)
  !
 !open(ioxyd,file="x_y.data",      status="replace",action="write")
 !open(iocnd,file="cell_node.data",status="replace",action="write")
 !open(iobcd,file="b_c.data",      status="replace",action="write")
 !!
 !!  read patran neutral file to get nnode, ncell, and nfbnd
 !!
 !read (iogrd,*)
 !read (iogrd,*)
 !read (iogrd,*) ndum,ndum,ndum,ndum, nnode, ncell
 !read (iogrd,*)
 !!
 !!   allocate memory for x and y and c_vtemp
 !allocate (x(nnode))
 !allocate (y(nnode))
 !allocate (c_vtemp(3,ncell))
 !!
 !do i = 1,nnode
 !  read (iogrd,*)
 !  read (iogrd,11) x(i),y(i)
 !  read (iogrd,*)
 !end do
 !!
 !do i = 1,ncell
 !  read (iogrd,*)
 !  read (iogrd,*)
 !  read (iogrd,*) c_vtemp(1,i),c_vtemp(2,i),c_vtemp(3,i)
 !end do
 !!
 !nfbnd = 0
 !do
 !  read (iogrd,*,end = 20)
 !  read (iogrd,*)
 !  read (iogrd,*)
 !  nfbnd = nfbnd+1
 !end do
 !20 continue
 !!
 !deallocate (x, y, c_vtemp)
 !!
 !!.. with nnode, ncell, and nfbnd known, now allocate memory
 !!... for geovar
 !allocate (xyz_nodes(2,-nfbnd:nnode))
 !allocate (c_v(3,-nfbnd:ncell), c_e(3,-nfbnd:ncell),c_nbr(3,ncell))
 !allocate (be_v(2,nfbnd), be_ic(nfbnd), be_e(nfbnd), &
 !          be_pe(nfbnd), be_pc(nfbnd), id_be(nfbnd))
 !!
 !!now local memory
 !allocate (area_(-nfbnd:ncell))
 !!...
 !rewind(iogrd)
 !read (iogrd,*)
 !read (iogrd,*)
 !read (iogrd,*)
 !read (iogrd,*)
 !!
 !do i = 1,nnode
 !  read (iogrd,*)
 !  read (iogrd,11)xyz_nodes(1,i),xyz_nodes(2,i)
 !  read (iogrd,*)
 !end do
 !!
 !do i = 1,ncell
 !  read (iogrd,*)
 !  read (iogrd,*)
 !  read (iogrd,*) c_v(1,i),c_v(2,i),c_v(3,i)
 !end do
 !!
 !!... we do not use the value i = 0, but initialze anyway
 !xyz_nodes(1,0) = zero
 !xyz_nodes(2,0) = zero
 !!
 !!   allocate memory for edge
 !allocate (edge(3,nfbnd))
 !!
 !do i = 1,nfbnd
 !  read (iogrd,*) ndum,be_ic(i)
 !  read (iogrd,10) ndum,edge(1,i),edge(2,i),edge(3,i)
 !  read (iogrd,*) rdum
 !  id_be(i) = int(rdum)
 !end do
 !!
 !!   determine which edge is on boundary
 !do i = 1,nfbnd
 !  !
 !  if (edge(3,i) == 0) then
 !    be_v(1,i) = c_v(1,be_ic(i))
 !    be_v(2,i) = c_v(2,be_ic(i))
 !  else
 !    if (edge(2,i) == 0) then
 !      be_v(1,i) = c_v(3,be_ic(i))
 !      be_v(2,i) = c_v(1,be_ic(i))
 !    else
 !      if (edge(1,i) == 0) then
 !        be_v(1,i) = c_v(2,be_ic(i))
 !        be_v(2,i) = c_v(3,be_ic(i))
 !      end if
 !    end if
 !  end if
 !  !
 !end do
 !!
 !!  write output for x and y coordinate data; cell to vertex data;
 !!  and boundary condition data.
 !!
 !write (ioxyd,*) nnode
 !do i = 1,nnode
 !  write (ioxyd,15) xyz_nodes(1,i),xyz_nodes(2,i)
 !end do
 !!
 !write (iocnd,*) ncell
 !do i = 1,ncell
 !  write (iocnd,211) c_v(1,i),c_v(2,i),c_v(3,i)
 !end do
 !!
 !! reset id_be to match cese convention
 !do i = 1,nfbnd
 !  !
 !  if (id_be(i) == 1) then ! inflow
 !    id_be(i) = 1
 !  else if (id_be(i) == 10) then ! inside wall
 !    id_be(i) = 0
 !  else if (id_be(i) == 11) then ! outside wall
 !    id_be(i) = 3
 !  else if (id_be(i) == 12) then ! outflow
 !    id_be(i) = 2
 !  end if
 !  !
 !end do
 !!
 !write (iobcd,*) nfbnd,0
 !write (iobcd,222) (id_be(i),i = 1,nfbnd)
 !write (iobcd,222) (be_ic(i),i = 1,nfbnd)
 !write (iobcd,222) ((be_v(k,i),k = 1,2),i = 1,nfbnd)
 !!
 !deallocate (edge)
 !!
 !close (ioxyd)
 !close (iocnd)
 !close (iobcd)
  close (iogrd,iostat=ierr)
  call io_error("read_patran_neutral_gridfile",trim(adjustl(gridfile)),2,ierr)
  !
  11 format (2e16.8)
  10 format (i9,i1,i1,i1)
  15 format (2e16.8)
  211 format (1x,3i10)
  222 format (1x,10i8)
  !
end subroutine read_patran_neutral_gridfile
!
!###############################################################################
!
subroutine check_node_ordering()
  !
  !.. Use Statements ..
  use geovar, only : ncell,nfbnd,nodes_of_cell_ptr, &
                     xyz_nodes,nodes_of_cell,bface
  !
  !.. Local Scalars ..
  integer   :: i,n,nf,n1,n2,ip1,ip2,ip3,ip4
  real(vbp) :: x31,y31,x42,y42,ss
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
    !
    ip1 = nodes_of_cell(n1+0)
    ip2 = nodes_of_cell(n1+1)
    ip3 = nodes_of_cell(n1+2)
    ip4 = nodes_of_cell(n2  ) ! for triangles, ip4 = ip3 because n2 = n1+2
    !
    x31 = xyz_nodes(1,ip3) - xyz_nodes(1,ip1)
    y31 = xyz_nodes(2,ip3) - xyz_nodes(2,ip1)
    !
    x42 = xyz_nodes(1,ip4) - xyz_nodes(1,ip2) ! for triangles x42 = x3-x2
    y42 = xyz_nodes(2,ip4) - xyz_nodes(2,ip2) ! for triangles y42 = y3-y2
    !
    ! Take the cross product - x31 X x42
    !
    ss = x31*y42 - x42*y31
    !
    ! If the cross product is negative, this means the cell is in clockwise
    ! ordering.  Switch two nodes (2 and 3 for triangles; 2 and 4 for quads)
    ! to change the cell to counter-clockwise ordering.
    !
    if (ss < zero) then
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
      bface(4,nf) = merge(4,5,bface(3,nf)==3) - bface(4,nf)
      !
      ! Now we have to loop through the remaining boundary faces
      ! and apply this correction to bface(4) for any other boundary
      ! faces that have this same host cell. Otherwise, when it gets
      ! to the other boundary face > nf with this same host cell, the
      ! node ordering will already be fixed and it will not make it
      ! into this ss<zero conditional to apply a correction to bface(4).
      !
      do i = nf+1,nfbnd
        if (bface(2,i) == n) then
          bface(4,i) = merge(4,5,bface(3,i)==3) - bface(4,i)
        end if
      end do
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
    ip1 = nodes_of_cell(n1+0)
    ip2 = nodes_of_cell(n1+1)
    ip3 = nodes_of_cell(n1+2)
    ip4 = nodes_of_cell(n2  ) ! for triangles, ip4 = ip3 because n2 = n1+2
    !
    x31 = xyz_nodes(1,ip3) - xyz_nodes(1,ip1)
    y31 = xyz_nodes(2,ip3) - xyz_nodes(2,ip1)
    !
    x42 = xyz_nodes(1,ip4) - xyz_nodes(1,ip2) ! for triangles x42 = x3-x2
    y42 = xyz_nodes(2,ip4) - xyz_nodes(2,ip2) ! for triangles y42 = y3-y2
    !
    ! Take the cross product - x31 X x42
    !
    ss = x31*y42 - x42*y31
    !
    ! If the cross product is negative, this means the cell is in clockwise
    ! ordering.  Switch two nodes (2 and 3 for triangles; 2 and 4 for quads)
    ! to change the cell have counter-clockwise ordering.
    !
    if (ss < zero) then
      !
      nodes_of_cell(n1+1) = ip4 ! for triangles, ip4 = ip3
      nodes_of_cell(n2  ) = ip2
      !
    end if
    !
  end do
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
  use geovar, only : nr,ncell,nfbnd,nnode,n_totcel,n_totnod, &
                     bface,nodes_of_cell_ptr,nodes_of_cell,xyz_nodes
 !use ovar,   only : Period_Time,Period_Timesteps
 !use ovar,   only : Final_Time,Total_Periods
 !use ovar,   only : vxref,vyref,aref,dt_in,machref
  use ovar,   only : Grid_Min_X,Grid_Max_X,Grid_Length_X
  use ovar,   only : Grid_Min_Y,Grid_Max_Y,Grid_Length_Y
  !
  !.. Local Scalars ..
  integer   :: g1,g2,i,i1,i2,n1,n2,nf,ng,ni,np
  integer   :: gcell,icell,igeom,inode,iside
  integer   :: pf,pcell,p1,p2,pgeom,pside
  integer   :: np1,np2,ni1,ni2,ierr
  integer   :: Period_Timesteps_X,Period_Timesteps_Y
  real(vbp) :: Period_Time_X,Period_Time_Y
  real(vbp) :: theta,real_Total_Periods
  !
  !.. Local Arrays ..
  real(vbp), dimension(1:nr)     :: ptf,rnrm,ds
  real(vbp), dimension(1:nr)     :: ds1,ds2,pnrm,inrm
  real(vbp), dimension(1:nr,1:4) :: xyp
  !
  !.. Local Allocatable Arrays ..
  integer,   allocatable, dimension(:)   :: tmp_int
  real(vbp), allocatable, dimension(:,:) :: tmp_real
  !
  integer, parameter :: rotation = 0
  !
continue
  !
  ! Before we start adding ghost cells, lets get the extents of the grid domain
  ! NOTE: ASSUME GRID BOUNDARIES OF CONSTANT X AND Y
  !
  Grid_Min_X = minval(xyz_nodes(1,:))
  Grid_Max_X = maxval(xyz_nodes(1,:))
  Grid_Min_Y = minval(xyz_nodes(2,:))
  Grid_Max_Y = maxval(xyz_nodes(2,:))
  Grid_Length_X = Grid_Max_X - Grid_Min_X
  Grid_Length_Y = Grid_Max_Y - Grid_Min_Y
  !
 !if (vxref < eps3) then
 !  !
 !  Period_Time_Y = Grid_Length_Y / machref
 !  !
 !  Period_Timesteps_Y = nint(Period_Time_Y/dt_in)
 !  !
 !  Period_Timesteps = Period_Timesteps_Y
 !  !
 !else if (vyref < eps3) then
 !  !
 !  Period_Time_X = Grid_Length_X / machref
 !  !
 !  Period_Timesteps_X = nint(Period_Time_X/dt_in)
 !  !
 !  Period_Timesteps = Period_Timesteps_X
 !  !
 !else
 !  !
 !  Period_Time_X = Grid_Length_X / (vxref/aref)
 !  Period_Time_Y = Grid_Length_Y / (vyref/aref)
 !  !
 !  Period_Timesteps_X = nint(Period_Time_X/dt_in)
 !  Period_Timesteps_Y = nint(Period_Time_Y/dt_in)
 !  !
 !  if (Period_Timesteps_X == Period_Timesteps_Y) then
 !    Period_Timesteps = Period_Timesteps_X
 !  else
 !    Period_Timesteps = Period_Timesteps_X * Period_Timesteps_Y &
 !                     / mod( max(Period_Timesteps_X,Period_Timesteps_Y), &
 !                            min(Period_Timesteps_X,Period_Timesteps_Y) )
 !  end if
 !  !
 !end if
 !!
 !Period_Time = dt_in*real(Period_Timesteps,vbp)
 !!
 !if (Total_Periods /= 0) then
 !  Final_Time = Period_Time*real(Total_Periods,vbp)
 !else
 !  real_Total_Periods = Final_Time/Period_Time
 !  Total_Periods = nint(real_Total_Periods)
 !  if (abs(real_Total_Periods-real(Total_Periods,vbp)) > ten*eps6) then
 !    write (iout,2) real_Total_Periods
 !  end if
 !end if
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
  allocate ( tmp_int(1:n_totcel+1) , source=0 , stat=ierr )
  call alloc_error("adjust_for_ghost_cells","tmp_int",1,ierr)
  ! 
  tmp_int(1:size(nodes_of_cell_ptr)) = nodes_of_cell_ptr
  call move_alloc( from=tmp_int , to=nodes_of_cell_ptr )
  !
  ! Set the number of nodes for each ghost cell
  !
  do nf = 1,nfbnd
    !
    if (bface(1,nf) /= bc_periodic) then
      !
      ! Not a periodic boundary, number of nodes is the same as interior
      !
      inode = merge(3,4,bface(3,nf)==3)
      !
    else
      !
      pf = bface(5,nf)
      pcell = bface(2,pf)
      !
      n1 = nodes_of_cell_ptr(pcell)+1
      n2 = nodes_of_cell_ptr(pcell+1)
      !
      inode = merge(3,4,n2-n1+1==3)
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
  allocate ( tmp_int(1:last(nodes_of_cell_ptr)) , source=0 , stat=ierr )
  call alloc_error("adjust_for_ghost_cells","tmp_int",1,ierr)
  !
  tmp_int(1:size(nodes_of_cell)) = nodes_of_cell
  call move_alloc( from=tmp_int , to=nodes_of_cell )
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
    gcell = nf + ncell
    !
    ! Node indices for the ghost cell
    !
    g1 = nodes_of_cell_ptr(gcell)+1
    g2 = nodes_of_cell_ptr(gcell+1)
    !
    ! Information about boundary face
    !
    icell = bface(2,nf)
    igeom = merge(3,4,bface(3,nf)==3)
    iside = bface(4,nf)
    !
    ! Starting node index for the interior cell
    !
    i1 = nodes_of_cell_ptr(icell)+1
    !
    ! Assign nodes on boundary from interior cell
    !
    nodes_of_cell(g1+0) = nodes_of_cell(i1+mod(iside,igeom))
    nodes_of_cell(g1+1) = nodes_of_cell(i1+iside-1)
    !
    ! Create new nodes for the remaining ghost cell nodes
    !
    do i = g1+2,g2
      inode = inode + 1
      nodes_of_cell(i) = inode
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
  ! Reallocate xyz_nodes to account for ghost nodes
  !
  allocate ( tmp_real(1:nr,1:n_totnod) , source=zero , stat=ierr )
  call alloc_error("adjust_for_ghost_cells","tmp_real",1,ierr)
  !
  tmp_real(1:nr,1:nnode) = xyz_nodes
  call move_alloc( from=tmp_real , to=xyz_nodes )
  !
  ! Loop over boundary faces to compute the
  ! coordinates of the newly created ghost nodes
  !
  do nf = 1,nfbnd
    !
    ! Ghost cell and the corresponding node indices
    !
    gcell = nf + ncell
    !
    g1 = nodes_of_cell_ptr(gcell)+1
    g2 = nodes_of_cell_ptr(gcell+1)
    !
    ! Interior cell and the corresponding node indices
    !
    icell = bface(2,nf)
    !
    i1 = nodes_of_cell_ptr(icell)+1
    i2 = nodes_of_cell_ptr(icell+1)
    !
    ! Information about boundary face
    !
    igeom = merge(3,4,bface(3,nf)==3)
    iside = bface(4,nf)
    !
    ! Conditional to treat periodic boundary conditions separately
    !
    if (bface(1,nf) /= bc_periodic) then
      !
      ! This is not a periodic boundary face
      !
      ! Nodes on face
      !
      n1 = nodes_of_cell(g1)
      n2 = nodes_of_cell(g1+1)
      !
      ! Center of boundary face by averaging the coordinates of the two boundary nodes
      !
      ptf(:) = half*(xyz_nodes(:,n1)+xyz_nodes(:,n2))
      !
      ! Outward unit normal to boundary face. NOTE, we are using the
      ! ordering of the boundary nodes for the ghost cell which is
      ! opposite that of the interior cell so the normal is the negative
      ! of the regular normal equation. (We are essentially finding the
      ! inward unit normal of the boundary face within the ghost cell)
      !
      rnrm(1) = xyz_nodes(2,n1) - xyz_nodes(2,n2)
      rnrm(2) = xyz_nodes(1,n2) - xyz_nodes(1,n1)
      rnrm = unit_vector(rnrm)
      !
      ! Coordinates for nodes g1 and g1+1 of the ghost cell have already been defined
      ! so we need to define nodes g1+2 and g2 (g1+2==g2 for triangular cells).
      !
      if (igeom == 3) then
        !
        ! Ghost node 3 for triangle
        !
        ni = nodes_of_cell( i1 + mod(iside+1,igeom) )
        ng = nodes_of_cell( g1 + 2 )
        !
        ! Get the coordinates of ghost node 3
        !
        xyz_nodes(:,ng) = merge(rotate(xyz_nodes(:,n1),xyz_nodes(:,n2),xyz_nodes(:,ni)), &
                                reflect(xyz_nodes(:,ni),ptf,rnrm), &
                                rotation == 1)
        !
      else if (igeom == 4) then
        !
        ! Ghost node 3 for quadrilateral
        !
        ng = nodes_of_cell( g1 + 2 )
        ni = nodes_of_cell( i1 + mod(iside+2-max(0,min(1,rotation)),igeom) )
        !
        ! Get the coordinates of ghost node 3
        !
        xyz_nodes(:,ng) = merge(rotate(xyz_nodes(:,n1),xyz_nodes(:,n2),xyz_nodes(:,ni)), &
                                reflect(xyz_nodes(:,ni),ptf,rnrm), &
                                rotation == 1)
        !
        ! Ghost node 4 for quadrilateral
        !
        ng = nodes_of_cell( g1 + 3 )
        ni = nodes_of_cell( i1 + mod(iside+1+max(0,min(1,rotation)),igeom) )
        !
        ! Get the coordinates of ghost node 4
        !
        xyz_nodes(:,ng) = merge(rotate(xyz_nodes(:,n1),xyz_nodes(:,n2),xyz_nodes(:,ni)), &
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
      pcell = bface(2,pf)
      !
      p1 = nodes_of_cell_ptr(pcell)+1
      p2 = nodes_of_cell_ptr(pcell+1)
      !
      ! Information about the linked face
      !
      pgeom = merge(3,4,bface(3,pf)==3)
      pside = bface(4,pf)
      !
      ! Node coordinates for periodic cell
      !
      xyp(:,1:pgeom) = xyz_nodes(:,nodes_of_cell(p1:p2))
      !
      ! Nodes on the boundary face for pcell
      !
      np1 = nodes_of_cell( p1 + pside - 1 )
      np2 = nodes_of_cell( p1 + mod(pside,pgeom) )
      !
      ! Nodes on the boundary face for icell
      !
      ni1 = nodes_of_cell( i1 + iside - 1)
      ni2 = nodes_of_cell( i1 + mod(iside,igeom) )
      !
      ! Center of boundary face for pcell
      !
      ptf(:) = half * (xyz_nodes(:,np1) + xyz_nodes(:,np2))
      !
      ! Outward unit normal for boundary face of pcell
      !
      pnrm(1) = xyz_nodes(2,np2) - xyz_nodes(2,np1)
      pnrm(2) = xyz_nodes(1,np1) - xyz_nodes(1,np2)
      pnrm = unit_vector(pnrm)
      !
      ! Outward unit normal for boundary face of icell
      !
      inrm(1) = xyz_nodes(2,ni2) - xyz_nodes(2,ni1)
      inrm(2) = xyz_nodes(1,ni1) - xyz_nodes(1,ni2)
      inrm = unit_vector(inrm)
      !
      ! Get the angle between the two outward unit normals
      !
      theta = acos(dot_product(pnrm,inrm))
      !
      ! Translation of the pcell nodes so that ptf is at the origin
      !
      do i = 1,pgeom
        xyp(:,i) = xyp(:,i) - ptf(:)
      end do
      !
      ! Rotate xyp about ptf theta+180 degrees
      !
      do i = 1,pgeom
        xyp(:,i) = rotate_2d( xyp(:,i) , theta+pi )
      end do
      !
      ! Translation of the pcell nodes back to original ptf
      !
      do i = 1,pgeom
        xyp(:,i) = xyp(:,i) + ptf(:)
      end do
      !
      ! Find the distances between the two corresponding nodes of
      ! the two boundary faces. The two boundary faces should now be
      ! parallel to each other but with opposite outward unit normals
      ! and therefore the distances should be equal. 
      !
      ds1(:) = xyz_nodes(:,ni1) - xyp(:,1+mod(pside,pgeom))
      ds2(:) = xyz_nodes(:,ni2) - xyp(:,pside)
      !
      if (abs(ds1(1)-ds2(1))>eps6 .or. abs(ds1(2)-ds2(2))>eps6) then
        write (*,1) icell,pcell,abs(ds1-ds2)
      end if
      !
      ! Move the interior nodes of pcell to the ghost nodes of icell
      !
      do i = 3,pgeom
        !
        ng = nodes_of_cell( g1 + i-1 )
        np = nodes_of_cell( p1 + mod(pside+i-2,pgeom) )
        !
        xyz_nodes(:,ng) = xyz_nodes(:,np) + ds1(:)
        !
      end do
      !
    end if
    !
  end do
  !
  ! If using multiple processors, save the global grid arrays
  !
  call save_global_grid
  !
  1 format ("ICELL = ",i0,"  PCELL = ",i0,"   (delta x,delta y) = ",2es24.11)
  2 format (/," ********************* WARNING!!! *********************",/, &
              "  IT APPEARS THE FINAL TIME DOES NOT CORRESPOND TO THE",/, &
              "  EXACT END OF A PERIOD FOR THE VORTEX PROBLEM. THE",/, &
              "  FINAL TIME CORRESPONDS TO",f12.5," TOTAL PERIODS.",/, &
              " ********************* WARNING!!! *********************",/)
  !
end subroutine adjust_for_ghost_cells
!
!###############################################################################
!
subroutine save_global_grid
  !
  !.. Use Statements ..
  use geovar
  use flowvar, only : global_usp,global_vsp
  !
  !.. Local Scalars ..
  integer :: n,np,n1,n2,ierr
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
  call alloc_error("save_global_grid","global_xyz_nodes",1,ierr)
  !
  n1 = size(nodes_of_cell_ptr)
  allocate ( global_nofc_ptr(1:n1) , source=nodes_of_cell_ptr , stat=ierr )
  call alloc_error("save_global_grid","global_nofc_ptr",1,ierr)
  !
  n1 = size(nodes_of_cell)
  allocate ( global_nofc(1:n1) , source=nodes_of_cell , stat=ierr )
  call alloc_error("save_global_grid","global_nofc",1,ierr)
  !
  ! Create the global pts_in_cell_ptr. This will be done later
  ! but a global version is never made before the global grid
  ! is partitioned and each processor creates its own separate
  ! local grid.
  !
  allocate ( global_pinc_ptr(1:n_totcel+1) , source=0 , stat=ierr )
  call alloc_error("save_global_grid","global_pinc_ptr",1,ierr)
  !
  do n = 1,n_totcel
    np = nodes_of_cell_ptr(n+1) - nodes_of_cell_ptr(n)
    global_pinc_ptr(n+1) = global_pinc_ptr(n) + merge(n_sp2dt,n_sp2dq,np==3)
  end do
  !
  n_global_solpts = global_pinc_ptr(ncell+1)
  n_global_totpts = last(global_pinc_ptr)
  !
  ! Now the only remaining global grid array to create is global_xyz
  ! We will have to do a collective communication of the local xyz
  ! arrays to accomplish this after they are created. We will allocate
  ! the array here though since we have the parameters needed to do so.
  !
  allocate ( global_xyz(1:nr,1:n_global_totpts) , source=zero , stat=ierr )
  call alloc_error("save_global_grid","global_xyz",1,ierr)
  !
  ! Also allocate the global solution arrays for use later
  !
  allocate ( global_usp(1:neqns,1:n_global_totpts) , source=zero , stat=ierr )
  call alloc_error("save_global_grid","global_usp",1,ierr)
  allocate ( global_vsp(1:neqns,1:n_global_totpts) , source=zero , stat=ierr )
  call alloc_error("save_global_grid","global_vsp",1,ierr)
  !
  ! Create the global cell offset array
  !
  allocate ( global_cell_offset(1:n_global_totcel) , source=0 , stat=ierr )
  call alloc_error("save_global_grid","global_cell_offset",1,ierr)
  !
  do n = 1,n_global_totcel-1
    np = global_nofc_ptr(n+1) - global_nofc_ptr(n)
    global_cell_offset(n+1) = global_cell_offset(n) + &
                              merge(n_sp2dt,n_sp2dq,np==3)
  end do
  !
end subroutine save_global_grid
!
!###############################################################################
!
subroutine read_restart_file
  !
  !.. Use Statements ..
  use geovar,       only : ncell,n_sp2dt,n_sp2dq,pts_in_cell_ptr, &
                           n_global_cell,global_nofc_ptr,global_cell_offset
  use ovar,         only : itrst,tmrst,dtrst
  use flowvar,      only : usp,global_usp
  use parallel_mod, only : cell_map,cell_map_ptr
  !
  !.. Local Scalars ..
  integer :: n,nc,nc1,n1,n2,np,gp,ierr
  !
continue
  !
  open(iorst,file="out.rst",status="old",action="read", &
             form="unformatted",iostat=ierr)
  call io_error("read_restart_file","out.rst",1,ierr)
  !
  read (iorst) itrst, tmrst, dtrst
  read (iorst) global_usp
  !
  close (iorst,iostat=ierr)
  call io_error("read_restart_file","out.rst",2,ierr)
  !
  ! If using multiple processors, copy the solution data from 
  ! only those cells that are local to this processor into usp
  !
  if (ncpu > 1) then
    !
    nc1 = cell_map_ptr(mypnum+1)
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
      gp = global_cell_offset( cell_map(nc1+nc) )
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
  integer :: ierr
  character(len=40) :: fname
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
    open (iorst,file=trim(adjustl(fname)),status="replace",action="write", &
                form='unformatted',iostat=ierr)
    call io_error("write_restart_file",trim(adjustl(fname)),1,ierr)
    !
    write (iorst) itcur, time, d_t
    write (iorst) global_usp
    !
    close (iorst,iostat=ierr)
    call io_error("write_restart_file",trim(adjustl(fname)),2,ierr)
    !
  end if
  !
  call mpi_barrier(mpi_comm_world,mpierr)
  !
end subroutine write_restart_file
!
!###############################################################################
!
subroutine write_tecplot_file(itnum)
  !
  ! Output the solution in tecplot format
  !
  !.. Use Statements ..
  use geovar,  only : ncell,nr,neqns,n_solpts,n_order,nodes_of_cell_ptr,xyz, &
                      pts_in_cell_ptr,tria_connectivity,quad_connectivity
  use ovar,    only : gam,pv_in
  use flowvar, only : vsp
  !
  !.. Formal Arguments ..
  integer, intent(in) :: itnum
  !
  !.. Local Scalars ..
  integer :: i,k,n,n1,n2,p1,ierr
  character(len=40) :: fname
  !
continue
  !
  if (itnum >= 0) then
    write (fname,'("out.tp.",i0.8,".dat")') itnum
  else
    fname = "out.tp.final.dat"
  end if
  !
  open (iotec,file=trim(adjustl(fname)),status="replace", &
              action="write",form="formatted",iostat=ierr)
  call io_error("write_tecplot_file",trim(adjustl(fname)),1,ierr)
  !
  write (iotec,1)
  write (iotec,2) n_solpts,ncell*n_order*n_order
  !
  ! Write the coordinates of the solution points
  !
  do n = 1,nr
    write (iotec,3) (xyz(n,i), i=1,n_solpts)
  end do
  !
  ! Write the flow variables
  !
  write (iotec,3) (vsp(1,i), i=1,n_solpts)
  do n = 2,neqns-1
    write (iotec,3) (vsp(n,i), i=1,n_solpts)
  end do
  write (iotec,3) (vsp(neqns,i), i=1,n_solpts)
  write (iotec,3) ( sqrt( selfdot(vsp(2:neqns-1,i))*vsp(1,i) / &
                         (gam*vsp(neqns,i))                  ), i=1,n_solpts)
  write (iotec,3) ( log( (vsp(neqns,i)/pv_in(neqns)) / &
                         ((vsp(1,i)/pv_in(1))**gam)  ), i=1,n_solpts)
  !
  ! Write out the connectivity matrix
  !
  do n = 1,ncell
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    p1 = pts_in_cell_ptr(n)+1
    do i = 1,n_order*n_order
      write (iotec,4) (p1+merge(tria_connectivity(k,i), &
                                quad_connectivity(k,i), &
                                n2-n1+1 == 3), k=1,4)
    end do
  end do
  !
  ! Finished writing tecplot solution file
  !
  close (iotec,iostat=ierr)
  call io_error("write_tecplot_file",trim(adjustl(fname)),2,ierr)
  !
  ! Format Statements
  !
  1 format ('VARIABLES = "X", "Y", "RHO", "u", "v", ', &
            '"Pressure", "Mach", "Entropy"')
  2 format ('ZONE T="Solution", N=',i8,', E=',i8, &
            ', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL' )
  3 format (5es18.10)
  4 format (10i8)
  !
end subroutine write_tecplot_file
!
!###############################################################################
!
subroutine write_tecplot_exact_sol()
  !
  ! Output the solution in tecplot format
  !
  !.. Use Statements ..
  use geovar,  only : nr,neqns,n_order,tria_connectivity,quad_connectivity
  use ovar,    only : itcur,Period_Timesteps,vxref,vyref,aref,time, &
                      Grid_Length_X,Grid_Length_Y,Grid_Min_X,Grid_Min_Y, &
                      Period_Time,itestcase
  !
  !.. Aliased Use Statements ..
  use geovar,  only : ncell => n_global_cell
  use geovar,  only : n_solpts => n_global_solpts
  use geovar,  only : nodes_of_cell_ptr => global_nofc_ptr
  use geovar,  only : pts_in_cell_ptr => global_pinc_ptr
  use geovar,  only : xyz => global_xyz
  !
  !.. Local Scalars ..
  integer :: i,k,n,n1,n2,p1,ierr
  real(vbp) :: delx,dely,xoffset,yoffset
  logical, save :: is_first_write = true
  !
  !.. Local Allocatable Arrays ..
  real(vbp), allocatable, dimension(:,:) :: exact
  !
  !.. Local Parameters ..
  integer, parameter :: io = 21
  !
continue
  !
  if (is_first_write) then
    open (io,file="out.exact.tec.dat",status="replace", &
             action="write",iostat=ierr)
    call io_error("write_tecplot_exact_sol","out.exact.tec.dat",1,ierr)
    write (io,1)
  end if
  !
  allocate ( exact(1:neqns,1:n_solpts) , source=zero , stat=ierr )
  call alloc_error("write_tecplot_exact_sol","exact",1,ierr)
  !
  ! If running one of the vortex transport problems, compute the change
  ! in position of the vortex since the start of the current period.
  !
  if (any(itestcase == Transport_Problems)) call vortex_position(delx,dely)
  !
  ! Get the exact solution at the current time
  !
  do n = 1,n_solpts
    !
    xoffset = merge( delx , delx-Grid_Length_X , xyz(1,n)-delx >= Grid_Min_X )
    yoffset = merge( dely , dely-Grid_Length_Y , xyz(2,n)-dely >= Grid_Min_Y )
    !
    exact(:,n) = vortex2D( xyz(1,n)-xoffset , xyz(2,n)-yoffset )
    !
  end do
  !
  if (is_first_write) then
    write (io,2) time,n_solpts,ncell*n_order*n_order
  else
    write (io,3) time,n_solpts,ncell*n_order*n_order
  end if
  !
  ! Write the coordinates of the solution points
  !
  if (is_first_write) then
    do n = 1,nr
      write (io,4) (xyz(n,i), i=1,n_solpts)
    end do
  end if
  !
  ! Write the density
  !
  write (io,4) (exact(1,i), i=1,n_solpts)
  !
  ! Write out the connectivity matrix
  !
  if (is_first_write) then
    do n = 1,ncell
      n1 = nodes_of_cell_ptr(n)+1
      n2 = nodes_of_cell_ptr(n+1)
      p1 = pts_in_cell_ptr(n)+1
      do i = 1,n_order*n_order
        write (io,5) (p1+merge(tria_connectivity(k,i), &
                               quad_connectivity(k,i), &
                               n2-n1+1 == 3), k=1,4)
      end do
    end do
  end if
  !
  ! Finished writing tecplot solution file
  !
  deallocate ( exact , stat=ierr )
  call alloc_error("write_tecplot_exact_sol","exact",2,ierr)
  !
  is_first_write = fals
  !
  ! Format Statements
  !
  1 format ('VARIABLES = "X", "Y", "RHO"')
  2 format ('ZONE T="t=',f6.2,'", N=',i0,', E=',i0, &
            ', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL' )
  3 format ('ZONE T="t=',f6.2,'", N=',i0,', E=',i0, &
            ', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL', &
            ', VARSHARELIST = ([1,2]=1), CONNECTIVITYSHAREZONE = 1')
  4 format (5es18.10)
  5 format (10i8)
  !
end subroutine write_tecplot_exact_sol
!
!###############################################################################
!
subroutine write_final_error
  !
  use geovar, only : neqns,n_order,n_global_cell, &
                     n_global_solpts,global_cell_offset
  use ovar,   only : itestcase,flux_divergence, &
                     vnloc,vkloc,vlone,vlinf,vltwo, &
                     unloc,ukloc,ulone,ulinf,ultwo
  !
  !.. Aliased Use Statements ..
  use geovar, only : xyz => global_xyz
  !
  !.. Local Scalars ..
  integer :: i
  !
continue
  !
  write (iout,*)
  !
 !write (iout,33) itestcase,iopSP,iopCE,kSP,kCE,ibdvu,flux_divergence,ncell, &
 !                itcur,time,nloc,kloc,xyz(1,kloc),xyz(2,kloc)
  write (iout,34) itestcase,flux_divergence, &
                  n_order,n_global_cell,n_global_solpts, &
                  vnloc(1),global_cell_offset(vnloc(1))+vkloc(1), &
                  xyz(1,global_cell_offset(vnloc(1))+vkloc(1)), &
                  xyz(2,global_cell_offset(vnloc(1))+vkloc(1))
  !
  write (iout,35) (vlone(i), &
                   vlinf(i), &
                   vltwo(i), &
                   xyz(1,global_cell_offset(vnloc(i))+vkloc(i)), &
                   xyz(2,global_cell_offset(vnloc(i))+vkloc(i)), i=1,neqns-1)
  write (iout,36) (ulone(i), &
                   ulinf(i), &
                   ultwo(i), &
                   xyz(1,global_cell_offset(unloc(i))+ukloc(i)), &
                   xyz(2,global_cell_offset(unloc(i))+ukloc(i)), i=neqns+1,neqns+1)
  write (iout,37) (vlone(i), &
                   vlinf(i), &
                   vltwo(i), &
                   xyz(1,global_cell_offset(vnloc(i))+vkloc(i)), &
                   xyz(2,global_cell_offset(vnloc(i))+vkloc(i)), i=neqns+1,neqns+1)
  !
 !33 format (/'itestcase =',i2,', iopSP,CE =',2i2,', kSP,CE =',2i2, &
 !           ', ibdvu =',i2,', flux_divergence =',i2,/,'ncell ',i7,', itcur ', i7, &
 !           ', time',f8.2, /,'max density err at k,j = ', i7, i7,'  at x, y ', 2e10.2)
 !34 format (/'itestcase =',i2,', kSP =',i2,', flux_divergence =',i2,/, &
 !            'ncell ',i7,', itcur ',i7,', time',f8.2, /, &
 !            'max density err at k,j = ', i7, i4,'  at x, y ', 2e10.2)
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
subroutine write_final_vortex_error
  !
  use geovar, only : neqns,n_order,n_global_cell, &
                     n_global_solpts,global_cell_offset
  use ovar,   only : itestcase,flux_divergence,vortex_error
  !
  !.. Aliased Use Statements ..
  use geovar, only : xyz => global_xyz
  !
  !.. Local Scalars ..
  integer :: i,n,m
  integer :: nloc,kloc,gloc
  !
  !.. Local Arrays ..
  integer, dimension(1:neqns+1,1:2) :: loc
  !
continue
  !
  write (iout,*)
  !
  do n = 1,2
    do m = 1,neqns+1
      loc(m,n) = global_cell_offset(nint(vortex_error(m,n,4))) + vortex_error(m,n,5)
    end do
  end do
  !
  nloc = nint(vortex_error(1,1,4))
  kloc = nint(vortex_error(1,1,5))
  gloc = global_cell_offset(nloc) + kloc
  !
  write (iout,34) itestcase,flux_divergence, &
                  n_order,n_global_cell,n_global_solpts, &
                  nloc,gloc, xyz(1,gloc), xyz(2,gloc)
  !
  ! Write out Gauss quadrature error
  !
  write (iout,35) "Gauss  ", &
                  (vortex_error(m,1,1),vortex_error(m,1,3),vortex_error(m,1,2), &
                  xyz(1,loc(m,1)),xyz(2,loc(m,1)), m=1,neqns+1)
  !
  ! Write out Lobatto quadrature error
  !
  write (iout,35) "Lobatto", &
                  (vortex_error(m,2,1),vortex_error(m,2,3),vortex_error(m,2,2), &
                  xyz(1,loc(m,2)),xyz(2,loc(m,2)), m=1,neqns+1)
  !
  34 format ('itestcase =',i2,/, &
             'flux_divergence  =',i2,//, &
             'Order of solution polynomial = ',i10,/, &
             'Number of grid cells         = ',i10,/, &
             'Number of solution points    = ',i10,//, &
             'Max density error in cell = ',i10,/, &
             '        at solution point = ',i10,/, &
             '         with coordinates = ',2es14.6,/)
  35 format (/,' ',a,'  -       L1           Max            L2     @      x             y',/,&
  ' ------------------------------------------------------------------------------',/,&
  ' Density  ',5es14.6,/,&
  ' Vx       ',5es14.6,/,&
  ' Vy       ',5es14.6,/,&
  ' Pressure ',5es14.6,/,&
  ' Entropy  ',5es14.6,/,&
  ' ------------------------------------------------------------------------------',/)
  !
end subroutine write_final_vortex_error
!
!###############################################################################
!
subroutine write_final_advec_error
  !
  use geovar, only : neqns,n_order,n_global_cell, &
                     n_global_solpts,global_cell_offset
  use ovar,   only : itestcase,flux_divergence,advec_error
  !
  !.. Aliased Use Statements ..
  use geovar, only : xyz => global_xyz
  !
  !.. Local Scalars ..
  integer :: i,n,m
  integer :: nloc,kloc,gloc
  !
  !.. Local Arrays ..
  integer, dimension(1:3,1:3) :: loc
  !
continue
  !
  write (iout,*)
  !
  do n = 1,3
    do m = 1,3
      loc(m,n) = global_cell_offset(nint(advec_error(m,n,4))) + advec_error(m,n,5)
    end do
  end do
  !
  nloc = nint(advec_error(1,1,4))
  kloc = nint(advec_error(1,1,5))
  gloc = global_cell_offset(nloc) + kloc
  !
  write (iout,34) itestcase,flux_divergence, &
                  n_order,n_global_cell,n_global_solpts, &
                  nloc,gloc, xyz(1,gloc), xyz(2,gloc)
  !
  ! Write out Gauss quadrature error
  !
  write (iout,35) "Gauss  ", &
                  (advec_error(m,1,1),advec_error(m,1,3),advec_error(m,1,2), &
                  xyz(1,loc(m,1)),xyz(2,loc(m,1)), m=1,3)
  !
  ! Write out Radau quadrature error
  !
  write (iout,35) "Radau  ", &
                  (advec_error(m,2,1),advec_error(m,2,3),advec_error(m,2,2), &
                  xyz(1,loc(m,2)),xyz(2,loc(m,2)), m=1,3)
  !
  ! Write out Lobatto quadrature error
  !
  write (iout,35) "Lobatto", &
                  (advec_error(m,3,1),advec_error(m,3,3),advec_error(m,3,2), &
                  xyz(1,loc(m,3)),xyz(2,loc(m,3)), m=1,3)
  !
  34 format ('itestcase =',i2,/, &
             'flux_divergence  =',i2,//, &
             'Order of solution polynomial = ',i10,/, &
             'Number of grid cells         = ',i10,/, &
             'Number of solution points    = ',i10,//, &
             'Max density error in cell = ',i10,/, &
             '        at solution point = ',i10,/, &
             '         with coordinates = ',2es14.6,/)
  35 format (/,' ',a,'  -       L1           Max            L2     @      x             y',/,&
  ' ------------------------------------------------------------------------------',/,&
  ' Density  ',5es14.6,/,&
  ' Pressure ',5es14.6,/,&
  ' Entropy  ',5es14.6,/,&
  ' ------------------------------------------------------------------------------',/)
  !
end subroutine write_final_advec_error
!
!###############################################################################
!
subroutine write_entropy_error
  !
  use geovar, only : ncell,n_sp1de
  use ovar
  !
  !
continue
  !
 !write (iout,34) itestcase,iopSP,iopCE,ioparcL,kSP,kCE,ibdvu,flux_divergence
  write (iout,35) itestcase,n_sp1de,flux_divergence
  !
  write (iout,38) ncell,nentmax,kentmax,xentmax,yentmax
  !
  write (iout,40) entmax,entL2
  !
  34 format (/'case',i2,', op SP CE',2i2,', arcL',i2, &
              ', kSP CE',2i2, ', ibdvu',i2,' flux_divergence',i2)
  35 format (/'case',i2,', kSP ',i2,' flux_divergence',i2)
  38 format ('ncell',i7, ' max En at j,k ',2i8, ', x,y ', 2f8.3)
  40 format ('max en', e10.3, ', L2 en', e10.3/)
  !
end subroutine write_entropy_error
!
!###############################################################################
!
subroutine write_results_file(iflag,icount)
  !
  use geovar,  only : neqns,n_global_solpts,p_order
  use ovar,    only : ultwo,vltwo,work_units
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iflag
  integer, intent(in) :: icount
  !
  !.. Local Scalars ..
  integer  :: i,ierr
  real(wp) :: comp_cost
  !
  !.. Local Arrays ..
  real(wp), dimension(1:neqns+1) :: u2,v2
  !
continue
  !
  ! Open the results file only if iflag is +1
  !
  if (iflag == 1) then
    !
    open (iores,file="results.dat",status="replace",action="write",iostat=ierr)
    call io_error("write_results_file","results.dat",1,ierr)
    !
    write (iores,1) p_order, one/sqrt(real(n_global_solpts,kind=wp))
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
  write (iores,2) icount,(log10(u2(i)),i=1,neqns+1), &
                         (log10(v2(i)),i=2,neqns+1), &
                         comp_cost
  !
  ! Close the results file only if iflag is -1
  !
  if (iflag == -1) then
    close (iores,iostat=ierr)
    call io_error("write_results_file","results.dat",2,ierr)
  end if
  !
 !call write_tecplot_exact_sol
  !
  1 format ('VARIABLES= "Timestep" "rho" "rhou" "rhov" "Energy" &
           &"Entropy" "u" "v" "Pressure" "||V||" "Work units"',/, &
            'Zone T="Residual Error"',/, &
            'AUXDATA order="',i0,'"',/, &
            'AUXDATA h="',f10.6,'"')
  2 format (i10,*(es17.9))
  !
end subroutine write_results_file
!
!###############################################################################
!
subroutine write_vortex_results_file(iflag,icount)
  !
  use geovar,  only : neqns,n_global_solpts,n_order,n_global_cell
  use ovar,    only : vortex_error,work_units,Period_Timesteps
  use ovar,    only : iop_solpts,dt_in,iRKtype,num_stages
  use ovar,    only : init_loc
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iflag
  integer, intent(in) :: icount
  !
  !.. Local Scalars ..
  integer   :: i,n,m,ierr
  real(vbp) :: h,comp_cost,period
 !character(len=11) :: rk_string
 !character(len=15) :: zonetitle
  !
  !.. Local Arrays ..
  real(vbp), dimension(1:neqns+1,1:2,1:3) :: e2
  !
continue
  !
  h = log10(one/sqrt(real(n_global_solpts,vbp)))
  !
 !write (zonetitle,2) n_order, merge("Gauss","Lobatto",iop_solpts==Gauss_Legendre)
 !rk_string = merge(char_rk(iRKtype),clas_rk(num_stages),iRKtype==1)
  !
  if (iflag == 1) then
    open (iores,file=trim(error_file),status="replace", &
                action="write",iostat=ierr)
    call io_error("write_vortex_results_file",trim(error_file),1,ierr)
    write (iores,4)
    write (iores,3) n_order, &
                    trim(merge("Gauss  ","Lobatto",iop_solpts==Legendre_Gauss)), &
                    nint(sqrt(real(n_global_cell,vbp))),char_init(init_loc), &
                    n_order, h, dt_in, &
                    trim(merge(clas_rk(num_stages),char_rk(iRKtype),iRKtype==1))
    return
  end if
  !
  comp_cost = work_units
  period = real(icount,vbp)/real(Period_Timesteps,vbp)
  !
  e2(:,:,:) = vortex_error(:,:,1:3)
  !
  where (e2 <= zero) e2 = one
  !
  write (iores,1) period,((log10(e2(m,n,2)), m=1,neqns+1), n=1,2), &
                         ((log10(e2(m,n,1)), m=1,neqns+1), n=1,2), &
                         ((log10(e2(m,n,3)), m=1,neqns+1), n=1,2), &
                         comp_cost
  !
  if (iflag == -1) then
    close (iores,iostat=ierr)
    call io_error("write_vortex_results_file",trim(error_file),2,ierr)
  end if
  !
 !call write_tecplot_exact_sol
  !
  1 format (es14.6,30es17.9,es14.6)
  2 format ("P",i0," - ",a)
  3 format ('ZONE',/, &
            'T="P',i0,' - ',a,' ',i0,' ',a,'"',/, &
            'AUXDATA p="',i0,'"',/, &
            'AUXDATA h="',es14.7,'"',/, &
            'AUXDATA dt="',es14.6,'"',/, &
            'AUXDATA TimeMethod="',a,'"')
  4 format ('VARIABLES = "Period"',/, &
            '"L<sub>2</sub> Density Gauss"',/, &
            '"L<sub>2</sub> V<sub>x</sub> Gauss"',/, &
            '"L<sub>2</sub> V<sub>y</sub> Gauss"',/, &
            '"L<sub>2</sub> Pressure Gauss"',/, &
            '"L<sub>2</sub> Entropy Gauss"',/, &
            '"L<sub>2</sub> Density Lobatto"',/, &
            '"L<sub>2</sub> V<sub>x</sub> Lobatto"',/, &
            '"L<sub>2</sub> V<sub>y</sub> Lobatto"',/, &
            '"L<sub>2</sub> Pressure Lobatto"',/, &
            '"L<sub>2</sub> Entropy Lobatto"',/, &
            '"L<sub>1</sub> Density Gauss"',/, &
            '"L<sub>1</sub> V<sub>x</sub> Gauss"',/, &
            '"L<sub>1</sub> V<sub>y</sub> Gauss"',/, &
            '"L<sub>1</sub> Pressure Gauss"',/, &
            '"L<sub>1</sub> Entropy Gauss"',/, &
            '"L<sub>1</sub> Density Lobatto"',/, &
            '"L<sub>1</sub> V<sub>x</sub> Lobatto"',/, &
            '"L<sub>1</sub> V<sub>y</sub> Lobatto"',/, &
            '"L<sub>1</sub> Pressure Lobatto"',/, &
            '"L<sub>1</sub> Entropy Lobatto"',/, &
            '"L<math>%</math> Density Gauss"',/, &
            '"L<math>%</math> V<sub>x</sub> Gauss"',/, &
            '"L<math>%</math> V<sub>y</sub> Gauss"',/, &
            '"L<math>%</math> Pressure Gauss"',/, &
            '"L<math>%</math> Entropy Gauss"',/, &
            '"L<math>%</math> Density Lobatto"',/, &
            '"L<math>%</math> V<sub>x</sub> Lobatto"',/, &
            '"L<math>%</math> V<sub>y</sub> Lobatto"',/, &
            '"L<math>%</math> Pressure Lobatto"',/, &
            '"L<math>%</math> Entropy Lobatto"',/, &
            '"Work Units"')
  !
end subroutine write_vortex_results_file
!
!###############################################################################
!
subroutine write_advec_results_file(iflag,icount)
  !
  use geovar,  only : n_global_solpts,n_order,n_global_cell
  use ovar,    only : advec_error,work_units,Period_Timesteps
  use ovar,    only : iop_solpts,dt_in,iRKtype,num_stages
  use ovar,    only : init_loc
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iflag
  integer, intent(in) :: icount
  !
  !.. Local Scalars ..
  integer   :: i,n,m,ierr
  real(vbp) :: h,comp_cost,period
 !character(len=11) :: rk_string
 !character(len=15) :: zonetitle
  !
  !.. Local Arrays ..
  real(vbp), dimension(1:3,1:3,1:3) :: e2
  !
continue
  !
  h = log10(one/sqrt(real(n_global_solpts,vbp)))
  !
 !write (zonetitle,2) n_order, merge("Gauss","Lobatto",iop_solpts==Gauss_Legendre)
 !rk_string = merge(char_rk(iRKtype),clas_rk(num_stages),iRKtype==1)
  !
  if (iflag == 1) then
    open (iores,file=trim(error_file),status="replace", &
                action="write",iostat=ierr)
    call io_error("write_advec_results_file",trim(error_file),1,ierr)
    write (iores,4)
    write (iores,3) n_order, &
                    trim(merge("Gauss  ","Lobatto",iop_solpts==Legendre_Gauss)), &
                    nint(sqrt(real(n_global_cell,vbp))),char_init(init_loc), &
                    n_order, h, dt_in, &
                    trim(merge(clas_rk(num_stages),char_rk(iRKtype),iRKtype==1))
    return
  end if
  !
  comp_cost = work_units
  period = real(icount,vbp)/real(Period_Timesteps,vbp)
  !
  e2(:,:,:) = advec_error(:,:,1:3)
  !
  where (e2 <= zero) e2 = one
  !
  write (iores,1) period,((e2(m,n,2), n=1,3), m=1,3), &
                         ((e2(m,n,1), n=1,3), m=1,3), &
                         ((e2(m,n,3), n=1,3), m=1,3), &
                         comp_cost
  !
  if (iflag == -1) then
    close (iores,iostat=ierr)
    call io_error("write_advec_results_file",trim(error_file),2,ierr)
  end if
  !
 !call write_tecplot_exact_sol
  !
  1 format (es14.6,27es17.9,es14.6)
  2 format ("P",i0," - ",a)
  3 format ('ZONE',/, &
            'T="P',i0,' - ',a,' ',i0,' ',a,'"',/, &
            'AUXDATA p="',i0,'"',/, &
            'AUXDATA h="',es14.7,'"',/, &
            'AUXDATA dt="',es14.6,'"',/, &
            'AUXDATA TimeMethod="',a,'"')
  4 format ('VARIABLES = "Period"',/, &
            '"L<sub>2</sub> Density Gauss"',/, &
            '"L<sub>2</sub> Density Radau"',/, &
            '"L<sub>2</sub> Density Lobatto"',/, &
            '"L<sub>2</sub> Pressure Gauss"',/, &
            '"L<sub>2</sub> Pressure Radau"',/, &
            '"L<sub>2</sub> Pressure Lobatto"',/, &
            '"L<sub>2</sub> Entropy Gauss"',/, &
            '"L<sub>2</sub> Entropy Radau"',/, &
            '"L<sub>2</sub> Entropy Lobatto"',/, &
            '"L<sub>1</sub> Density Gauss"',/, &
            '"L<sub>1</sub> Density Radau"',/, &
            '"L<sub>1</sub> Density Lobatto"',/, &
            '"L<sub>1</sub> Pressure Gauss"',/, &
            '"L<sub>1</sub> Pressure Radau"',/, &
            '"L<sub>1</sub> Pressure Lobatto"',/, &
            '"L<sub>1</sub> Entropy Gauss"',/, &
            '"L<sub>1</sub> Entropy Radau"',/, &
            '"L<sub>1</sub> Entropy Lobatto"',/, &
            '"L<math>%</math> Density Gauss"',/, &
            '"L<math>%</math> Density Radau"',/, &
            '"L<math>%</math> Density Lobatto"',/, &
            '"L<math>%</math> Pressure Gauss"',/, &
            '"L<math>%</math> Pressure Radau"',/, &
            '"L<math>%</math> Pressure Lobatto"',/, &
            '"L<math>%</math> Entropy Gauss"',/, &
            '"L<math>%</math> Entropy Radau"',/, &
            '"L<math>%</math> Entropy Lobatto"',/, &
            '"Work Units"')
  !
end subroutine write_advec_results_file
!
!###############################################################################
!
subroutine write_cgns_solution(iflag)
  !
  !.. Use Statements ..
  use geovar,       only : neqns,n_sp2dt,ncell,n_solpts,pts_in_cell_ptr
  use ovar,         only : gam,pv_in,itestcase,output_solution
  use metvar,       only : dx,dy,dr_mat,ds_mat,jacobian
  use flowvar,      only : usp,vsp,uoldsp,dusp
  use parallel_mod, only : collect_global_solution
  use generic_mod,  only : vsp2u,usp2v,entropy_pv_sp
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iflag
  !
  !.. Local Scalars ..
  integer   :: n,n1,n2,np,nw,k
  real(vbp) :: dudr,duds,dvdr,dvds
  !
continue
  !
  ! First we need to change the conserved solution variables for output
  ! NOTE: Primitive solution variables will stay the same
  !
  ! usp(1,:) = density    -> Mach number
  ! usp(2,:) = x-momentum -> Vorticity
  ! usp(3,:) = y-momentum -> ?
  ! usp(4,:) = energy     -> Entropy
  !
  ! Mach number
  !
  if (check_gradients) then
    !
    do k = 1,n_solpts
      do n = 1,neqns
        usp(n,k) = dusp(1,n,k)
        vsp(n,k) = dusp(2,n,k)
      end do
    end do
    !
  else
    !
    if (itestcase /= Density_Transport) then
      do k = 1,n_solpts
        usp(1,k) = sqrt( selfdot(vsp(2:neqns-1,k)) * &
                         vsp(1,k)/(gam*vsp(neqns,k)) )
      end do
    else
      usp(1,1:n_solpts) = one
    end if
    !
    ! Entropy
    !
    if (itestcase /= Density_Transport) then
      do k = 1,n_solpts
        usp(4,k) = entropy_pv_sp( vsp(:,k) )
      end do
    else
      usp(4,1:n_solpts) = one
    end if
    !
    ! Vorticity
    !
    do n = 1,ncell
      !
      n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
      n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
      np = n2-n1+1
      !
      nw = merge(1,2,np==n_sp2dt) ! check if cell n is tria or quad
      !
      do k = n1,n2
        !
        dudr = dot_product(dr_mat(1:np,k-n1+1,nw),vsp(2,n1:n2))
        duds = dot_product(ds_mat(1:np,k-n1+1,nw),vsp(2,n1:n2))
        dvdr = dot_product(dr_mat(1:np,k-n1+1,nw),vsp(3,n1:n2))
        dvds = dot_product(ds_mat(1:np,k-n1+1,nw),vsp(3,n1:n2))
        !
        usp(2,k) = ( ( dy(k)%ds*dvdr - dy(k)%dr*dvds ) - &
                     ( dx(k)%dr*duds - dx(k)%ds*dudr ) ) / jacobian(k)
        !
      end do
      !
    end do
    !
  end if
  !
  ! Collect the global solution arrays into global_usp and global_vsp
  !
  call collect_global_solution
  !
  ! Output the nodal and/or modal solution to the cgns file(s)
  !
  if (mypnum == 0) then
   !do n = 1,size(pv_in)
   !  write (iout,'(" pv_in(",i0,") = ",g12.4)') n,pv_in(n)
   !end do
    if (output_solution == 1) then
      call write_nodal_cgns(iflag)
    else if (output_solution == 2) then
      call write_modal_cgns(iflag)
    else
      call write_modal_cgns(iflag)
      call write_nodal_cgns(iflag)
    end if
   !call write_tecplot_exact_sol
  end if
  !
  call mpi_barrier(mpi_comm_world,mpierr)
  !
  ! Copy the saved conserved variables back into usp
  !
  usp(:,1:n_solpts) = uoldsp(:,1:n_solpts)
  vsp(:,1:n_solpts) = usp2v(usp(:,1:n_solpts))
  !
end subroutine write_cgns_solution
!
!###############################################################################
!
subroutine write_nodal_cgns(iflag)
  !
  !.. Use Statements ..
  use geovar,  only : p_order,n_order,n_sp2dt,n_sp2dq,neqns
  use geovar,  only : tria_connectivity,quad_connectivity
  use ovar,    only : time,aref,uref,vxref,vyref,gam,rhoref
  use ovar,    only : itcur,Period_Timesteps,Period_Time,itestcase
  use ovar,    only : Grid_Length_X,Grid_Length_Y,Grid_Min_X,Grid_Min_Y
  use generic_mod, only : usp2v,project_to_porder,entropy_pv_sp
  !
  !.. Aliased Use Statements ..
  use geovar,  only : ncell => n_global_cell
  use geovar,  only : n_solpts => n_global_solpts
  use geovar,  only : nodes_of_cell_ptr => global_nofc_ptr
  use geovar,  only : pts_in_cell_ptr => global_pinc_ptr
  use geovar,  only : xyz => global_xyz
  use flowvar, only : usp => global_usp
  use flowvar, only : vsp => global_vsp
  use flowvar, only : uoldsp
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iflag
  !
  !.. Local Scalars ..
  integer   :: i,j,n,n1,n2,ierr,p1,k,m
  integer   :: nelem,kind_real,nsp,nw,np
  integer   :: n_outpts,no1d,noc,nopt,nopq
  real(vbp) :: cvel,sn,se
  real(vbp) :: delx,dely,xoffset,yoffset
  !
  character(len=32) :: tname
  !
  !.. Local Arrays ..
  integer, dimension(1:2) :: idata
  integer, dimension(1:3) :: idim1
  real(wp), dimension(1:n_sp2dq,1:neqns) :: tusp
  !
  !.. Local Allocatable Arrays ..
  integer,   allocatable, dimension(:,:) :: ipelem
  real(wp),  allocatable, dimension(:)   :: vartmp
  real(vhp), allocatable, dimension(:)   :: var
  !
continue
  !
  noc = size(tria_connectivity,2)
  !
  no1d = nint(sqrt(real(noc,kind=wp)))
 !no1d = n_order
  !
 !noc = no1d*no1d
  nopt = (no1d+1)*(no1d+2)/2
  nopq = (no1d+1)*(no1d+1)
  nelem = ncell*noc
  !
  n_outpts = 0
  do n = 1,ncell
    np = nodes_of_cell_ptr(n+1) - nodes_of_cell_ptr(n)
    n_outpts = n_outpts + merge(nopt,nopq,np==3)
  end do
  !
  ! Variables specific to the Vortex transport cases
  !
  if ( any(itestcase == Transport_Problems) ) then
    !
    ! If running one of the vortex transport problems, compute the change
    ! in position of the vortex since the start of the current period.
    !
    if (any(itestcase == Transport_Problems)) call vortex_position(delx,dely)
    !
    cvel = aref*reciprocal(uref)
    !
  else
    !
    cvel = one
    !
  end if
  !
  kind_real = merge(RealDouble,RealSingle,kind(var)==vbp)
  !
  allocate ( vartmp(1:n_solpts) , source=zero , stat=ierr )
  call alloc_error("write_nodal_cgns","vartmp",1,ierr)
  allocate ( var(1:n_outpts) , source=real(0,vhp) , stat=ierr )
  call alloc_error("write_nodal_cgns","var",1,ierr)
  !
  ! Open the cgns output file if iflag is 1
  !
  if (iflag == 1) then
    !
    ! Initialize the zone counter
    !
    cgnsnodalzone_counter = 0
    !
    call cg_set_file_type_f(file_type,ierr)
    call cg_open_f(fname_nodal,MODE_WRITE,ifile_n,ierr)
    call io_error("write_nodal_cgns",trim(fname_nodal),1,ierr)
    !
    ! Set the CGNS base (who knows what this is for)
    !
    call cg_base_write_f(ifile_n,"recon_base",2,2,ibase_n,ierr)
    !
    ! Create the connectivity array
    !
    allocate ( ipelem(1:4,1:nelem) , source=0 , stat=ierr )
    call alloc_error("write_nodal_cgns","ipelem",1,ierr)
    !
    j = 0
    p1 = 1
    do n = 1,ncell
      n1 = nodes_of_cell_ptr(n)+1
      n2 = nodes_of_cell_ptr(n+1)
     !p1 = pts_in_cell_ptr(n)+1
     !do i = 1,n_order*n_order
      do i = 1,noc
        j = j + 1
        ipelem(:,j) = p1 + merge(tria_connectivity(:,i), &
                                 quad_connectivity(:,i), &
                                 n2-n1+1 == 3)
      end do
      p1 = p1 + merge(nopt,nopq,n2-n1+1 == 3)
    end do
    !
    ! Write the title for this solution
    !
    idim1(:) = [n_outpts,nelem,0]
    call cg_zone_write_f(ifile_n,ibase_n,"Solution",idim1, &
                         Unstructured,izone_n,ierr)
    !
    ! Write the x-coordinates
    !
    do i = 1,n_solpts
      vartmp(i) = xyz(1,i)
    end do
    var(:) = real( interpolate_for_output(vartmp) , kind(var) )
    call cg_coord_write_f(ifile_n,ibase_n,izone_n,kind_real,"CoordinateX", &
                          var(1:n_outpts),icoord_n,ierr)
    !
    ! Write the y-coordinates
    !
    do i = 1,n_solpts
      vartmp(i) = xyz(2,i)
    end do
    var(:) = real( interpolate_for_output(vartmp) , kind(var) )
    call cg_coord_write_f(ifile_n,ibase_n,izone_n,kind_real,"CoordinateY", &
                          var(1:n_outpts),icoord_n,ierr)
    !
    ! Write the grid connectivity
    !
    call cg_section_write_f(ifile_n,ibase_n,izone_n,"Quadrilaterals", &
                            QUAD_4,1,nelem,0,ipelem(1:4,1:nelem),isect_n,ierr)
    !
    ! Deallocate ipelem since we are done with it
    !
    deallocate ( ipelem , stat=ierr )
    call alloc_error("write_nodal_cgns","ipelem",2,ierr)
    !
#ifndef IO_FOR_NAS
    !
    ! If not using NAS, we want to possibly look at the solution as the simulation
    ! is still running. Using CGNS with HDF5, Tecplot wont be able to open the
    ! solution file if we keep it open during the whole simulation so we need to
    ! close it each time after we are done writing to it.
    !
    ! Close the CGNS file temporarily
    !
    call cg_close_f(ifile_n,ierr)
    call io_error("write_nodal_cgns",trim(fname_nodal),2,ierr)
    !
#endif
    !
  end if
  !
 !write (tname,'(i0.4," Time = ",es9.2)') cgnsnodalzone_counter,time
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
#ifndef IO_FOR_NAS
  !
  ! If not using NAS, we want to possibly look at the solution as the simulation
  ! is still running. Using CGNS with HDF5, Tecplot wont be able to open the
  ! solution file if we keep it open during the whole simulation so we need to
  ! close it each time after we are done writing to it.
  !
  ! Reopen the CGNS file in modification mode
  !
  call cg_open_f(fname_nodal,CG_MODE_MODIFY,ifile_n,ierr)
  call io_error("write_nodal_cgns",trim(fname_nodal),3,ierr)
  !
#endif
  !
  ! Write the flow variables
  !
  call cg_sol_write_f(ifile_n,ibase_n,izone_n,tname,Vertex,isol_n,ierr)
  !
  if (check_gradients) then
    !
    ! Density
    !
    do i = 1,n_solpts
      vartmp(i) = usp(1,i)
    end do
    var(:) = real( interpolate_for_output(vartmp) , kind(var) )
    call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"Density dx", &
                          var(1:n_outpts),ifld_n,ierr)
    !
    ! x-velocity
    !
    do i = 1,n_solpts
      vartmp(i) = (usp(2,i) - uoldsp(2,i)/uoldsp(1,i)*usp(1,i)) / uoldsp(1,i)
    end do
    var(:) = real( interpolate_for_output(vartmp) , kind(var) )
    call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"VelocityX dx", &
                          var(1:n_outpts),ifld_n,ierr)
    !
    ! y-velocity
    !
    do i = 1,n_solpts
      vartmp(i) = (usp(3,i) - uoldsp(3,i)/uoldsp(1,i)*usp(1,i)) / uoldsp(1,i)
    end do
    var(:) = real( interpolate_for_output(vartmp) , kind(var) )
    call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"VelocityY dx", &
                          var(1:n_outpts),ifld_n,ierr)
    !
    ! Density
    !
    do i = 1,n_solpts
      vartmp(i) = vsp(1,i)
    end do
    var(:) = real( interpolate_for_output(vartmp) , kind(var) )
    call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"Density dy", &
                          var(1:n_outpts),ifld_n,ierr)
    !
    ! x-velocity
    !
    do i = 1,n_solpts
      vartmp(i) = (vsp(2,i) - uoldsp(2,i)/uoldsp(1,i)*vsp(1,i)) / uoldsp(1,i)
    end do
    var(:) = real( interpolate_for_output(vartmp) , kind(var) )
    call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"VelocityX dy", &
                          var(1:n_outpts),ifld_n,ierr)
    !
    ! y-velocity
    !
    do i = 1,n_solpts
      vartmp(i) = (vsp(3,i) - uoldsp(3,i)/uoldsp(1,i)*vsp(1,i)) / uoldsp(1,i)
    end do
    var(:) = real( interpolate_for_output(vartmp) , kind(var) )
    call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"VelocityY dy", &
                          var(1:n_outpts),ifld_n,ierr)
    !
  else
    !
    ! Density
    !
    do i = 1,n_solpts
      vartmp(i) = vsp(1,i)
    end do
    var(:) = real( interpolate_for_output(vartmp) , kind(var) )
    call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"Density", &
                          var(1:n_outpts),ifld_n,ierr)
    !
    ! x-velocity
    !
    do i = 1,n_solpts
      vartmp(i) = vsp(2,i)
    end do
    var(:) = real( interpolate_for_output(vartmp) , kind(var) )
    call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"VelocityX", &
                          var(1:n_outpts),ifld_n,ierr)
    !
    ! y-velocity
    !
    do i = 1,n_solpts
      vartmp(i) = vsp(3,i)
    end do
    var(:) = real( interpolate_for_output(vartmp) , kind(var) )
    call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"VelocityY", &
                          var(1:n_outpts),ifld_n,ierr)
    !
    ! Pressure
    !
    do i = 1,n_solpts
      vartmp(i) = vsp(4,i)*cvel**2
    end do
    var(:) = real( interpolate_for_output(vartmp) , kind(var) )
    call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"Pressure", &
                          var(1:n_outpts),ifld_n,ierr)
    !
    ! Mach number
    !
    do i = 1,n_solpts
      vartmp(i) = usp(1,i)
    end do
    var(:) = real( interpolate_for_output(vartmp) , kind(var) )
    call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"Mach", &
                          var(1:n_outpts),ifld_n,ierr)
    !
    ! Entropy
    !
    do i = 1,n_solpts
      vartmp(i) = usp(4,i)
    end do
    var(:) = real( interpolate_for_output(vartmp) , kind(var) )
    call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"Entropy", &
                          var(1:n_outpts),ifld_n,ierr)
    !
    ! Vorticity
    !
    do i = 1,n_solpts
      vartmp(i) = usp(2,i)*cvel
    end do
    var(:) = real( interpolate_for_output(vartmp) , kind(var) )
    call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"Vorticity", &
                          var(1:n_outpts),ifld_n,ierr)
    !
  end if
! !
! ! Create the exact solution at each solution point
! !
! do i = 1,n_solpts
!   xoffset = merge(delx,delx-Grid_Length_X,xyz(1,i)-delx >= Grid_Min_X)
!   yoffset = merge(dely,dely-Grid_Length_Y,xyz(2,i)-dely >= Grid_Min_Y)
!   usp(:,i) = vortex2d( xyz(1,i)-xoffset , xyz(2,i)-yoffset )
! end do
! !
! call project_to_porder(ncell,pts_in_cell_ptr,usp)
! !
! ! Convert these to primitive variables
! !
! usp(1:neqns,1:n_solpts) = usp2v(usp(1:neqns,1:n_solpts))
! !
! ! Density error
! !
! do i = 1,n_solpts
!   vartmp(i) = usp(1,i)-vsp(1,i)
! end do
! var(:) = real( interpolate_for_output(vartmp) , kind(var) )
! call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"Error Density", &
!                       var(1:n_outpts),ifld_n,ierr)
  !
  ! x-velocity error
  !
! do i = 1,n_solpts
!   vartmp(i) = usp(2,i)-vsp(2,i)
! end do
! var(:) = real( interpolate_for_output(vartmp) , kind(var) )
! call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"Error VelocityX", &
!                       var(1:n_outpts),ifld_n,ierr)
  !
  ! y-velocity Error
  !
! do i = 1,n_solpts
!   vartmp(i) = usp(3,i)-vsp(3,i)
! end do
! var(:) = real( interpolate_for_output(vartmp) , kind(var) )
! call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"Error VelocityY", &
!                       var(1:n_outpts),ifld_n,ierr)
  !
  ! Entropy error
  !
 !do i = 1,n_solpts
 !  sn = entropy_pv_sp( vsp(:,i) )
 !  se = entropy_pv_sp( usp(:,i) )
 !  vartmp(i) = sn-se
 !end do
 !var(:) = real( interpolate_for_output(vartmp) , kind(var) )
 !call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"Error Entropy", &
 !                      var(1:n_outpts),ifld_n,ierr)
  !
  ! Pressure error
  !
 !do i = 1,n_solpts
 !  vartmp(i) = usp(4,i)-vsp(4,i)
 !end do
 !var(:) = real( interpolate_for_output(vartmp) , kind(var) )
 !call cg_field_write_f(ifile_n,ibase_n,izone_n,isol_n,kind_real,"Error Pressure", &
 !                      var(1:n_outpts),ifld_n,ierr)
  !
  ! Deallocate var array
  !
  deallocate ( var , stat=ierr )
  call alloc_error("write_nodal_cgns","var",2,ierr)
  !
#ifdef IO_FOR_NAS
  !
  ! If last iteration, close CGNS file
  !
  if (iflag == -1) then
    !
    ! Update the base iterative node
    !
    call cg_biter_write_f(ifile_n,ibase_n,'TimeIterValues',size(soltime),ierr)
    call cg_goto_f(ifile_n,ibase_n,ierr,"BaseIterativeData_t",1,"end")
    call cg_array_write_f("TimeValues",kind_real,1,size(soltime),soltime,ierr)
    !
    ! Update the zone iterative node
    !
    call cg_ziter_write_f(ifile_n,ibase_n,izone_n,"ZoneIterativeData",ierr)
    call cg_goto_f(ifile_n,ibase_n,ierr,"Zone_t",izone_n, &
                                        "ZoneIterativeData_t",1,"end")
    idata(:) = [32,size(solname)]
    call cg_array_write_f("FlowSolutionPointers",Character,2,idata,solname,ierr)
    !
    ! Identify the flow solution within the CGNS file as time dependent
    !
    call cg_simulation_type_write_f(ifile_n,ibase_n,TimeAccurate,ierr)
    !
    ! Close the CGNS FILE
    !
    call cg_close_f(ifile_n,ierr)
    call io_error("write_nodal_cgns",trim(fname_nodal),2,ierr)
    !
  end if
  !
#else
  !
  ! Update the base iterative node
  !
  call cg_biter_write_f(ifile_n,ibase_n,'TimeIterValues',size(soltime),ierr)
  call cg_goto_f(ifile_n,ibase_n,ierr,"BaseIterativeData_t",1,"end")
  call cg_array_write_f("TimeValues",kind_real,1,size(soltime),soltime,ierr)
  !
  ! Update the zone iterative node
  !
  call cg_ziter_write_f(ifile_n,ibase_n,izone_n,"ZoneIterativeData",ierr)
  call cg_goto_f(ifile_n,ibase_n,ierr,"Zone_t",izone_n, &
                                      "ZoneIterativeData_t",1,"end")
  idata(:) = [32,size(solname)]
  call cg_array_write_f("FlowSolutionPointers",Character,2,idata,solname,ierr)
  !
  ! Identify the flow solution within the CGNS file as time dependent
  !
  call cg_simulation_type_write_f(ifile_n,ibase_n,TimeAccurate,ierr)
  !
  ! If not using NAS, we want to possibly look at the solution as the simulation
  ! is still running. Using CGNS with HDF5, Tecplot wont be able to open the
  ! solution file if we keep it open during the whole simulation so we need to
  ! close it each time after we are done writing to it.
  !
  ! Close the CGNS FILE
  !
  call cg_close_f(ifile_n,ierr)
  call io_error("write_nodal_cgns",trim(fname_nodal),2,ierr)
  !
#endif
  !
  ! Update the zone counter
  !
  cgnsnodalzone_counter = cgnsnodalzone_counter + 1
  !
end subroutine write_nodal_cgns
!
!###############################################################################
!
subroutine write_modal_cgns(iflag)
  !
  !.. Use Statements ..
  use geovar,  only : nr,n_sp1de,n_sp2dq,vand_quad_inv,vand_quad
  use ovar,    only : time
  !
  !.. Aliased Use Statements ..
  use geovar,  only : ncell => n_global_cell
  use geovar,  only : n_solpts => n_global_solpts
  use geovar,  only : nodes_of_cell => global_nofc
  use geovar,  only : nodes_of_cell_ptr => global_nofc_ptr
  use geovar,  only : pts_in_cell_ptr => global_pinc_ptr
  use geovar,  only : xyz_nodes => global_xyz_nodes
  use flowvar, only : usp => global_usp
  use flowvar, only : vsp => global_vsp
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iflag
  !
  !.. Local Scalars ..
  integer   :: i,j,k,n,n1,n2,ierr,p1
  integer   :: npoin,kind_real
  !
  character(len=32) :: tname
  !
  !.. Local Arrays ..
  integer,   dimension(1:2) :: idata
  integer,   dimension(1:3) :: idim1
  real(vbp), dimension(1:n_sp1de+1) :: x1d
  real(vbp), dimension(1:nr,1:4) :: xy
  real(vhp), dimension(1:n_sp2dq) :: var2
  !
  !.. Local Allocatable Arrays ..
  integer,   allocatable, dimension(:,:) :: ipelem
  real(vhp), allocatable, dimension(:)   :: var
  real(vbp), allocatable, dimension(:,:) :: coord
  !
continue
  !
  npoin = ncell*(n_sp1de+1)*(n_sp1de+1)
  kind_real = merge(RealDouble,RealSingle,kind(var)==vbp)
  !
  ! If triangle cells exist, return
  !
  if (n_solpts /= ncell*n_sp2dq) return
  !
  ! Allocate var for cell-centered data
  !
  allocate ( var(1:npoin) , source=real(0,vhp) , stat=ierr )
  call alloc_error("write_modal_cgns","var",1,ierr)
  !
  ! Open the cgns output file
  !
  if (iflag == 1) then
    !
    ! Initialize the zone counter
    !
    cgnsmodalzone_counter = 0
    !
    call cg_set_file_type_f(file_type,ierr)
    call cg_open_f(fname_modal,MODE_WRITE,ifile_m,ierr)
    call io_error("write_modal_cgns",trim(fname_modal),1,ierr)
    !
    ! Set the CGNS base (who knows what this is for)
    !
    call cg_base_write_f(ifile_m,"recon_base",2,2,ibase_m,ierr)
    !
    ! Create the connectivity array
    !
    allocate ( ipelem(1:4,1:n_solpts) , source=0 , stat=ierr )
    call alloc_error("write_modal_cgns","ipelem",1,ierr)
    allocate ( coord(1:nr,1:npoin) , source=zero , stat=ierr )
    call alloc_error("write_modal_cgns","coord",1,ierr)
    !
    ! Create the connectivity array
    !
    k = 1
    do n = 1,ncell
      do j = 1,n_sp1de
        do i = 1,n_sp1de
          ipelem(1,k) = (n_sp1de+1)*(n_sp1de+1)*(n-1) + (n_sp1de+1)*(j-1) + i
          ipelem(2,k) = (n_sp1de+1)*(n_sp1de+1)*(n-1) + (n_sp1de+1)*(j-1) + i + 1
          ipelem(3,k) = (n_sp1de+1)*(n_sp1de+1)*(n-1) + (n_sp1de+1)*j + i + 1
          ipelem(4,k) = (n_sp1de+1)*(n_sp1de+1)*(n-1) + (n_sp1de+1)*j + i
          k = k + 1
        end do
      end do
    end do
    !
    ! Create the coordinates
    !
    do i = 1,n_sp1de+1
      x1d(i) = -one + two*real(i-1,vbp)/real(n_sp1de,vbp)
    end do
    !
    k = 1
    do n = 1,ncell
      !
      n1 = nodes_of_cell_ptr(n)+1
      n2 = nodes_of_cell_ptr(n+1)
      !
      do j = n1,n2
        xy(:,j-n1+1) = xyz_nodes(:,nodes_of_cell(j))
      end do
      !
      do j = 1,n_sp1de+1
        do i = 1,n_sp1de+1
          coord(:,k) = c0_25 * ( (one-x1d(j)) * (one-x1d(i)) * xy(:,1) + &
                                 (one-x1d(j)) * (one+x1d(i)) * xy(:,2) + &
                                 (one+x1d(j)) * (one+x1d(i)) * xy(:,3) + &
                                 (one+x1d(j)) * (one-x1d(i)) * xy(:,4) )
          k = k + 1
        end do
      end do
      !
    end do
    !
    ! Write the title for this solution
    !
    idim1(:) = [npoin,n_solpts,0]
    call cg_zone_write_f(ifile_m,ibase_m,"Solution",idim1,Unstructured,izone_m,ierr)
    !
    ! Write the x-coordinates
    !
    do i = 1,npoin
      var(i) = real(coord(1,i),kind(var))
    end do
    call cg_coord_write_f(ifile_m,ibase_m,izone_m,kind_real,"CoordinateX", &
                          var(1:npoin),icoord_m,ierr)
    !
    ! Write the y-coordinates
    !
    do i = 1,npoin
      var(i) = real(coord(2,i),kind(var))
    end do
    call cg_coord_write_f(ifile_m,ibase_m,izone_m,kind_real,"CoordinateY", &
                          var(1:npoin),icoord_m,ierr)
    !
    ! Write the grid connectivity
    !
    call cg_section_write_f(ifile_m,ibase_m,izone_m,"Quadrilaterals",QUAD_4,1,n_solpts,0, &
                            ipelem(1:4,1:n_solpts),isect_m,ierr)
    !
    deallocate ( ipelem , stat=ierr )
    call alloc_error("write_modal_cgns","ipelem",2,ierr)
    deallocate ( coord , stat=ierr )
    call alloc_error("write_modal_cgns","coord",2,ierr)
    !
#ifndef IO_FOR_NAS
    !
    ! If not using NAS, we want to possibly look at the solution as the simulation
    ! is still running. Using CGNS with HDF5, Tecplot wont be able to open the
    ! solution file if we keep it open during the whole simulation so we need to
    ! close it each time after we are done writing to it.
    !
    ! Close the CGNS file temporarily
    !
    call cg_close_f(ifile_m,ierr)
    call io_error("write_modal_cgns",trim(fname_modal),2,ierr)
    !
#endif
    !
  end if
  !
 !write (tname,'(i0.4," Time = ",es9.2)') cgnsmodalzone_counter,time
  write (tname,'("FlowSolution",i0)') cgnsmodalzone_counter
  !
  ! Reallocate the modtime and modname arrays so we
  ! can concatenate this solution onto those arrays
  !
  call reallocate(modtime,cgnsmodalzone_counter+1)
  call reallocate(modname,cgnsmodalzone_counter+1)
  modtime(size(modtime)) = real(time,kind(modtime))
  modname(size(modname)) = tname
  !
#ifndef IO_FOR_NAS
  !
  ! If not using NAS, we want to possibly look at the solution as the simulation
  ! is still running. Using CGNS with HDF5, Tecplot wont be able to open the
  ! solution file if we keep it open during the whole simulation so we need to
  ! close it each time after we are done writing to it.
  !
  ! Reopen the CGNS file in modification mode
  !
  call cg_open_f(fname_modal,CG_MODE_MODIFY,ifile_m,ierr)
  call io_error("write_modal_cgns",trim(fname_modal),3,ierr)
  !
#endif
  !
  ! Write the flow variables
  !
  call cg_sol_write_f(ifile_m,ibase_m,izone_m,tname,CellCenter,isol_m,ierr)
  !
  ! Convert solution variables to modal form by multiplying it
  ! with the inverse of the Vandermonde matrix and then output
  ! the modal form to the CGNS file
  !
  ! Density
  !
  do n = 1,ncell
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    var(n1:n2) = real( matmul( vand_quad_inv , vsp(1,n1:n2) ) , kind(var) )
  end do
  call cg_field_write_f(ifile_m,ibase_m,izone_m,isol_m,kind_real,"Density", &
                        var(1:n_solpts),ifld_m,ierr)
  !
  ! x-velocity
  !
  do n = 1,ncell
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    var(n1:n2) = real( matmul( vand_quad_inv , vsp(2,n1:n2) ) , kind(var) )
  end do
  call cg_field_write_f(ifile_m,ibase_m,izone_m,isol_m,kind_real,"VelocityX", &
                        var(1:n_solpts),ifld_m,ierr)
  !
  ! Entropy
  !
  do n = 1,ncell
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    var(n1:n2) = real( matmul( vand_quad_inv , usp(4,n1:n2) ) , kind(var) )
  end do
  call cg_field_write_f(ifile_m,ibase_m,izone_m,isol_m,kind_real,"Entropy", &
                        var(1:n_solpts),ifld_m,ierr)
  !
  ! Vorticity
  !
  do n = 1,ncell
    n1 = pts_in_cell_ptr(n)+1 ! beginning index for sol-pts in cell n
    n2 = pts_in_cell_ptr(n+1) ! ending    index for sol-pts in cell n
    var(n1:n2) = real( matmul( vand_quad_inv , usp(2,n1:n2) ) , kind(var) )
  end do
  call cg_field_write_f(ifile_m,ibase_m,izone_m,isol_m,kind_real,"Vorticity", &
                        var(1:n_solpts),ifld_m,ierr)
  !
  ! Deallocate var array
  !
  deallocate ( var , stat=ierr )
  call alloc_error("write_modal_cgns","var",2,ierr)
  !
#ifdef IO_FOR_NAS
  !
  ! If last iteration, close CGNS file
  !
  if (iflag == -1) then
    !
    ! Update the base iterative node
    !
    call cg_biter_write_f(ifile_m,ibase_m,"TimeIterValues",size(modtime),ierr)
    call cg_goto_f(ifile_m,ibase_m,ierr,"BaseIterativeData_t",1,"end")
    call cg_array_write_f("TimeValues",kind_real,1,size(modtime),modtime,ierr)
    !
    ! Update the zone iterative node
    !
    call cg_ziter_write_f(ifile_m,ibase_m,izone_m,"ZoneIterativeData",ierr)
    call cg_goto_f(ifile_m,ibase_m,ierr,"Zone_t",izone_m, &
                                        "ZoneIterativeData_t",1,"end")
    idata(:) = [32,size(modname)]
    call cg_array_write_f("FlowSolutionPointers",Character,2,idata,modname,ierr)
    !
    ! Identify the flow solution within the CGNS file as time dependent
    !
    call cg_simulation_type_write_f(ifile_m,ibase_m,TimeAccurate,ierr)
    !
    ! Close the CGNS file
    !
    call cg_close_f(ifile_m,ierr)
    call io_error("write_modal_cgns",trim(fname_modal),2,ierr)
    !
  end if
  !
#else
  !
  ! Update the base iterative node
  !
  call cg_biter_write_f(ifile_m,ibase_m,"TimeIterValues",size(modtime),ierr)
  call cg_goto_f(ifile_m,ibase_m,ierr,"BaseIterativeData_t",1,"end")
  call cg_array_write_f("TimeValues",kind_real,1,size(modtime),modtime,ierr)
  !
  ! Update the zone iterative node
  !
  call cg_ziter_write_f(ifile_m,ibase_m,izone_m,"ZoneIterativeData",ierr)
  call cg_goto_f(ifile_m,ibase_m,ierr,"Zone_t",izone_m, &
                                      "ZoneIterativeData_t",1,"end")
  idata(:) = [32,size(modname)]
  call cg_array_write_f("FlowSolutionPointers",Character,2,idata,modname,ierr)
  !
  ! Identify the flow solution within the CGNS file as time dependent
  !
  call cg_simulation_type_write_f(ifile_m,ibase_m,TimeAccurate,ierr)
  !
  ! If not using NAS, we want to possibly look at the solution as the simulation
  ! is still running. Using CGNS with HDF5, Tecplot wont be able to open the
  ! solution file if we keep it open during the whole simulation so we need to
  ! close it each time after we are done writing to it.
  !
  ! Close the CGNS file
  !
  call cg_close_f(ifile_m,ierr)
  call io_error("write_modal_cgns",trim(fname_modal),2,ierr)
  !
#endif
  !
  ! Update the zone counter
  !
  cgnsmodalzone_counter = cgnsmodalzone_counter + 1
  !
end subroutine write_modal_cgns
!
!###############################################################################
!
subroutine write_tecplot_partitions()
  !
  !.. Use Statements ..
  use geovar, only : ncell,nr,n_solpts,n_order,nodes_of_cell_ptr,xyz, &
                     pts_in_cell_ptr,tria_connectivity,quad_connectivity, &
                     n_totpts,n_totcel
  use ovar,    only : aref,uref
  use flowvar, only : vsp
  use parallel_mod, only : nonblocking_exchange
  !
  !.. Local Scalars ..
  integer   :: i,k,n,n1,n2,p1,p2,ierr
  real(vbp) :: cvel
  character(len=40) :: fname
  !
  !.. Local Allocatable Arrays ..
  real(vbp), allocatable, dimension(:) :: test_array
  !
continue
  !
  cvel = aref/uref
  !
  write (fname,6) mypnum
  !
  open (iotec,file=trim(adjustl(fname)),status="replace", &
              action="write",form="formatted",iostat=ierr)
  call io_error("write_tecplot_partitions",trim(adjustl(fname)),1,ierr)
  !
  write (iotec,9)
  !
  write (iotec,2) mypnum,n_totpts,n_totcel*n_order*n_order
  !
  do n = 1,nr
    write (iotec,4) (xyz(n,i),i=1,n_totpts)
  end do
  !
  write (iotec,4) (vsp(1,i),i=1,n_totpts)
  write (iotec,4) (vsp(2,i)*cvel,i=1,n_totpts)
  write (iotec,4) (vsp(3,i)*cvel,i=1,n_totpts)
  write (iotec,4) (vsp(4,i)*cvel**2,i=1,n_totpts)
  !
  do n = 1,n_totcel
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    p1 = pts_in_cell_ptr(n)+1
    do i = 1,n_order*n_order
      write (iotec,5) (p1+merge(tria_connectivity(k,i), &
                                quad_connectivity(k,i), &
                                n2-n1+1 == 3), k=1,4)
    end do
  end do
  !
  call nonblocking_exchange(vsp)
  !
  write (iotec,10) mypnum,n_totpts,n_totcel*n_order*n_order
  !
  write (iotec,4) (vsp(1,i),i=1,n_totpts)
  write (iotec,4) (vsp(2,i)*cvel,i=1,n_totpts)
  write (iotec,4) (vsp(3,i)*cvel,i=1,n_totpts)
  write (iotec,4) (vsp(4,i)*cvel**2,i=1,n_totpts)
  !
  close (iotec,iostat=ierr)
  call io_error("write_tecplot_partitions",trim(adjustl(fname)),2,ierr)
 !!
 !! Write individual cells as zones
 !!
 !write (fname,7) mypnum
 !!
 !open (iotec,file=trim(adjustl(fname)),status="replace", &
 !            action="write",form="formatted",iostat=ierr)
 !!
 !write (iotec,1)
 !write (iotec,9)
 !!
 !do n = 1,n_totcel
 !  !
 !  write (iotec,8) mypnum,n,n_totpts,n_order*n_order
 !  !
 !  do i = 1,nr
 !    write (iotec,4) (xyz(i,k),k=1,n_totpts)
 !  end do
 ! !write (iotec,4) (test_array(k),k=1,n_totpts)
 !  !
 !  n1 = nodes_of_cell_ptr(n)+1
 !  n2 = nodes_of_cell_ptr(n+1)
 !  p1 = pts_in_cell_ptr(n)+1
 !  do i = 1,n_order*n_order
 !    write (iotec,5) (p1+merge(tria_connectivity(k,i), &
 !                              quad_connectivity(k,i), &
 !                              n2-n1+1 == 3), k=1,4)
 !  end do
 !  !
 !end do
 !!
 !close (iotec,iostat=ierr)
  !
 !deallocate ( test_array , stat=ierr )
  !
  1 format ('VARIABLES = "X", "Y"')
  2  format ('ZONE T="Block ',i0,'", N=',i0,', E=',i0, &
             ', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL')
  10 format ('ZONE T="Block ',i0,'", N=',i0,', E=',i0, &
             ', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL', &
             ', VARSHARELIST=([1,2]=1), CONNECTIVITYSHAREZONE=1')
  3 format ('ZONE T="Ghost ',i0,'", N=',i0,', E=',i0, &
            ', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL')
  4 format (5es18.10)
  5 format (10i8)
  6 format ("partition_sp.",i0,".tec")
  7 format ("partition_cell_sp.",i0,".tec")
  8  format ('ZONE T="Block ',i0,', Cell ',i0,'", N=',i0,', E=',i0, &
             ', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL')
  9 format ('VARIABLES = "X", "Y", "Density", "Vx", "Vy", "Pressure"')
  !
end subroutine write_tecplot_partitions
!
!###############################################################################
!
subroutine create_output_arrays
  !
continue
  !
end subroutine create_output_arrays
!
!###############################################################################
!
!pure function interpolate_for_output(array,npts,ityp) result(return_value)
pure function interpolate_for_output(array) result(return_value)
  !
  !... calculates d(utd)/dt
  !
  !.. Use Statements ..
  use geovar,  only : n_sp2dq
  use geovar,  only : ncell => n_global_cell
  use geovar,  only : pts_in_cell_ptr => global_pinc_ptr
  use geovar,  only : numer2gauss,numer2lobatto
  use ovar,    only : iop_outpts,iop_solpts
 !use geovar,  only : n_sp2dt,o_order
 !use geovar,  only : interp_s2o
  !
  !.. Formal Arguments ..
  real(wp), dimension(:), intent(in) :: array
  !
  !.. Function Return Value ..
  real(wp), dimension(1:size(array)) :: return_value
  !
  !.. Local Scalars ..
  integer :: nc,n1,n2,np,n
 !integer :: k,k1,n,n1,n2,nsp,nw,nop,nopq,nopt,iloc
  !
  !.. Local Arrays ..
  real(wp), dimension(1:n_sp2dq,1:n_sp2dq) :: interp_mat
  !
continue 
  !
  !
  !
  if (iop_outpts == iop_solpts) then
    !
    !
    !
    do n = 1,size(array)
      return_value(n) = array(n)
    end do
    !
  else
    !
    !
    !
    if (iop_outpts == Legendre_Gauss_Lobatto) then
      !
      interp_mat(:,:) = numer2lobatto(:,:)
      !
    else if (iop_outpts == Legendre_Gauss) then
      !
      interp_mat(:,:) = numer2gauss(:,:)
      !
    else
      !
      do n = 1,size(array)
        return_value(n) = array(n)
      end do
      !
      return
      !
    end if
    !
    !
    !
    do nc = 1,ncell
      !
      n1 = pts_in_cell_ptr(nc)+1
      n2 = pts_in_cell_ptr(nc+1)
      np = n2-n1+1
      !
      if (np == n_sp2dq) then
        return_value(n1:n2) = matmul( interp_mat , array(n1:n2) )
      else
        return_value(n1:n2) = array(n1:n2)
      end if
      !
    end do
    !
  end if
  !
 !k1 = 0
 !!
 !if (npts == size(array)) then
 !  !
 !  do n = 1,npts
 !    return_value(n) = array(n)
 !  end do
 !  !
 !else
 !  !
 !  if (present(ityp)) then
 !    iloc = max(1,min(ityp,2))
 !  else
 !    iloc = 1
 !  end if
 !  !
 !  nopt = (o_order+1)*(o_order+2)/2
 !  nopq = (o_order+1)*(o_order+1)
 !  !
 !  do n = 1,ncell
 !    !
 !    n1 = pts_in_cell_ptr(n)+1
 !    n2 = pts_in_cell_ptr(n+1)
 !    nsp = n2-n1+1
 !    nw = merge(1,2,nsp==n_sp2dt)
 !    !
 !    nop = merge(nopt,nopq,nsp==n_sp2dt)
 !    !
 !    do k = 1,nop
 !      return_value(k1+k) = dot_product( interp_s2o(1:nsp,k,nw,iloc) , array(n1:n2) )
 !    end do
 !    !
 !    k1 = k1 + nop
 !    !
 !  end do
 !  !
 !end if
  !
end function interpolate_for_output
!
!###############################################################################
!
subroutine reallocate1_i4(array,new_size)
  integer, intent(in) :: new_size
  integer(i4), allocatable, dimension(:), intent(inout) :: array
  integer(i4), allocatable, dimension(:) :: temp_array
  integer :: ierr
continue
  allocate ( temp_array(1:new_size) , source=0_i4 , stat=ierr )
  call alloc_error("reallocate1_i4","temp_array",1,ierr)
  if (allocated(array)) temp_array(1:size(array)) = array
  call move_alloc( from=temp_array , to=array )
end subroutine reallocate1_i4
!
!###############################################################################
!
subroutine reallocate1_i8(array,new_size)
  integer, intent(in) :: new_size
  integer(i8), allocatable, dimension(:), intent(inout) :: array
  integer(i8), allocatable, dimension(:) :: temp_array
  integer :: ierr
continue
  allocate ( temp_array(1:new_size) , source=0_i8 , stat=ierr )
  call alloc_error("reallocate1_i8","temp_array",1,ierr)
  if (allocated(array)) temp_array(1:size(array)) = array
  call move_alloc( from=temp_array , to=array )
end subroutine reallocate1_i8
!
!###############################################################################
!
subroutine reallocate1_r4(array,new_size)
  integer, intent(in) :: new_size
  real(vhp), allocatable, dimension(:), intent(inout) :: array
  real(vhp), allocatable, dimension(:) :: temp_array
  integer :: ierr
continue
  allocate ( temp_array(1:new_size) , source=0.0_vhp , stat=ierr )
  call alloc_error("reallocate1_r4","temp_array",1,ierr)
  if (allocated(array)) temp_array(1:size(array)) = array
  call move_alloc( from=temp_array , to=array )
end subroutine reallocate1_r4
!
!###############################################################################
!
subroutine reallocate1_r8(array,new_size)
  integer, intent(in) :: new_size
  real(vbp), allocatable, dimension(:), intent(inout) :: array
  real(vbp), allocatable, dimension(:) :: temp_array
  integer :: ierr
continue
  allocate ( temp_array(1:new_size) , source=0.0_vbp , stat=ierr )
  call alloc_error("reallocate1_r8","temp_array",1,ierr)
  if (allocated(array)) temp_array(1:size(array)) = array
  call move_alloc( from=temp_array , to=array )
end subroutine reallocate1_r8
!
!###############################################################################
!
subroutine reallocate1_c(array,new_size)
  integer, intent(in) :: new_size
  character(len=*), allocatable, dimension(:), intent(inout) :: array
  character(len=len(array)), allocatable, dimension(:) :: temp_array
  integer :: ierr
continue
  allocate ( temp_array(1:new_size) , stat=ierr )
  call alloc_error("reallocate1_c","temp_array",1,ierr)
  if (allocated(array)) temp_array(1:size(array)) = array
  call move_alloc( from=temp_array , to=array )
end subroutine reallocate1_c
!
!###############################################################################
!
include "Functions/last.f90"
include "Functions/vortex.f90"
!
!###############################################################################
!
end module io_mod
!
!###############################################################################
!
