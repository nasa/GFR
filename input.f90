module input_mod
  !
  use module_kind_types
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
  !
  ! ###########################
  ! ###  Module Parameters  ###
  ! ###########################
  !
  !
  ! ################################
  ! ###  Module Local Variables  ###
  ! ################################
  !
  !
  ! ########################
  ! ###  Input Namelist  ###
  ! ########################
  !
  !.. Input Namelist ..
  ! GRID INFORMATION
  namelist / inpt_nml / grid_format, gridfile
  ! BASE FR OPTIONS
  namelist / inpt_nml / n_order, iop_solpts, iop_flxpts, iop_tripts, &
                     flux_divergence, iop_correction
  ! SOLVER OPTIONS
  namelist / inpt_nml / iRKtype, num_stages, flux_method
  ! TIME STEPPING INFORMATION
  namelist / inpt_nml / num_ts, Final_Time, cfl, dt_in
  ! OVER-INTEGRATION OPTIONS
  namelist / inpt_nml / p_order, q_order, iop_qdtpts, dflux_opt
  ! PARALLEL OPTIONS
  namelist / inpt_nml / metis_option
  ! INITIALIZATION / TEST CASE
  namelist / inpt_nml / itestcase
  ! RESTART OPTIONS
  namelist / inpt_nml / irestart, restart_interval
  ! OUTPUT OPTIONS
  namelist / inpt_nml / output_dir, iter_out_interval, output_solution, &
                     output_interval, results_interval, o_order
  ! REFERENCE CONDITIONS
  namelist / inpt_nml / rhoref, tref, pref, machref, aoaref, gam, rgasref
  ! FILTER OPTIONS
  namelist / inpt_nml / Filter_Option, Filter_etac, Filter_Order, Filter_Alpha
  ! BOUNDING ERROR OPTIONS
  namelist / inpt_nml / bound_error, error_bound_limit
  ! INITIALIZATION/ERROR LOCATIONS
  namelist / inpt_nml / init_loc, error_loc
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
  use geovar, only : o_order
  use geovar, only : q_order,n_qp1de,n_qp2dt,n_qp2dq,using_quadrature
  use ovar,   only : itestcase,num_stages,irestart, &
                     flux_divergence,iop_correction, &
                     num_ts,cfl,dt_in,gam,output_solution, &
                     output_interval,restart_interval, &
                     results_interval,flux_method,dflux_opt, &
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
  namelist / input / iRKtype, num_stages, flux_method
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
  Final_Time = zero
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
   !write (iout,54) num_ts, Final_Time, cfl, dt_in, flux_method, iter_out_interval
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
                     results_interval,flux_method,dflux_opt, &
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
  namelist / input / iRKtype, num_stages, flux_method
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
  Final_Time = zero
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
  read (inpt,*) num_ts, Final_Time, cfl, dt_in, flux_method, iter_out_interval
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
    write (iout,54) num_ts, Final_Time, cfl, dt_in, flux_method, iter_out_interval
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
end module input_mod
