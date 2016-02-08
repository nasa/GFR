program Euler_tri
  !
  use module_kind_types
  !
  use geovar, only : ncell,n_global_cell
  use geovar, only : nnode,n_global_node
  use geovar, only : nfbnd,n_global_fbnd
  use geovar, only : n_solpts,n_global_solpts
  use geovar, only : bface,xyz_nodes,nface
  use geovar, only : cell_geom,cell_order,global_cell_order
  use geovar, only : nodes_of_cell_ptr,nodes_of_cell
  use geovar, only : geovar_memory_usage
  !
  use ovar, only : max_wall_time,wall_time_left
  use ovar, only : walltime_expiring
  use ovar, only : metis_option,convergence_reached
  use ovar, only : total_volume,local_volume,itestcase
  use ovar, only : nst,num_stages,iRKtype,time,itcur,num_ts
  use ovar, only : governing_equations,visc_flux_method
  use ovar, only : results_interval,iter_out_interval
  use ovar, only : restart_interval,output_interval
  use ovar, only : dump_memory_usage,mms_opt
  use ovar, only : recon_cpu_time,this_is_final_timestep
  use ovar, only : ovar_memory_usage
  use ovar, only : write_TaylorGreen_8s_solution
  use ovar, only : write_TaylorGreen_full_solution
  !
  use flowvar, only : usp,uoldsp
  use flowvar, only : flowvar_memory_usage
  !
  use quadrature_mod, only : solpts_edge
  use quadrature_mod, only : quadrature_memory_usage
  !
  use io_mod, only : read_input_file
  use io_mod, only : write_final_error
  use io_mod, only : write_parallel_cgns
  use io_mod, only : write_results_file
  use io_mod, only : check_flux_point_orderings
  use io_mod, only : TaylorGreen_8s_vorticity_solution
  !
  use restart_mod, only : write_restart_file
  !
  use mappings_mod, only : get_solution_points
  !
  use connectivity_mod, only : get_grid_connectivity
  !
  use metrics_mod, only : get_metrics
  use metrics_mod, only : metrics_memory_usage
  !
  use vandermonde_mod, only : init_vandermonde_matrices
  use vandermonde_mod, only : vandermonde_memory_usage
  !
  use derivatives_mod, only : init_derivative_matrices
  use derivatives_mod, only : derivatives_memory_usage
  !
  use correction_mod, only : init_correction_matrices
  use correction_mod, only : correction_memory_usage
  !
  use interpolation_mod, only : init_interpolation_matrices
  use interpolation_mod, only : init_outerpolation_matrices
  use interpolation_mod, only : interpolation_memory_usage
  !
  use projection_mod, only : init_projection_matrices
  use projection_mod, only : projection_memory_usage
  !
  use filter_mod, only : init_filter_matrices
  use filter_mod, only : filter_memory_usage
  !
  use initialization_mod, only : initialize_flow
  !
  use module_bc,   only : get_faceusp_bc
  use module_bc,   only : get_facedusp_bc
!!  use module_bc,   only : compute_bc_profiles
  !
  use generic_mod, only : init_generic_mod
  use generic_mod, only : rgeterr
  use generic_mod, only : getmaxent
  use generic_mod, only : compute_final_error
  use generic_mod, only : compute_TaylorGreen_KEDR
  use generic_mod, only : generic_memory_usage
  !
  use time_mod, only : get_dtsp
  use time_mod, only : rresi
  use time_mod, only : rupdate_Classic_RK
  use time_mod, only : rupdate_CK4_RK
  use time_mod, only : rupdate_TVD2_RK
  use time_mod, only : rupdate_TVD3_RK
  use time_mod, only : rupdate_TVD4_RK
  !
  use flux_mod, only : interp_cv_to_flxpts
  use flux_mod, only : interp_dusp_to_flxpts
  use flux_mod, only : flux_gradient_at_solpts_with_qdtr
  use flux_mod, only : flux_gradient_at_solpts_visc
  use flux_mod, only : correct_interface_flux_with_qdtr
  use flux_mod, only : correct_interface_flux_visc
  use flux_mod, only : discontinuous_gradient_at_solpts
  use flux_mod, only : correct_gradients_at_solpts
  use flux_mod, only : correct_gradients_br2
  !
  use parallel_mod, only : partition_grid
  use parallel_mod, only : create_serial_cell_map
  use parallel_mod, only : parallel_memory_usage
  !
  use module_limiters, only : initialize_limiters
  use module_limiters, only : compute_cell_ave
  use module_limiters, only : limiters_memory_usage
  !
  use mms_mod, only : mms_src
  use mms_mod, only : initialize_mms
  !
  use channel_mod, only : channel_source
  !
  use bnd_profiles_mod, only : init_bnd_profile_structures
  use bnd_profiles_mod, only : compute_bnd_profiles
  use bnd_profiles_mod, only : bnd_profiles_memory_usage
  !
  use module_plot3d, only : plot3d_memory_usage
  !
  use module_gmsh, only : gmsh_memory_usage
  !
  use iso_c_binding, only : c_long
  !
  implicit none
  !
  character(len=*), parameter :: pname = "Euler_tri(main)"
  !
  integer  :: iter,i,ierr,k,npart
  integer  :: completed_timesteps
  real(wp) :: dtmin
  !
  character(len=50) :: memfile
  !
  integer, parameter :: iostp = 2683
  !
  integer, parameter :: reopen_CGNS = 2
  integer, parameter :: open_CGNS = 1
  integer, parameter :: add_to_CGNS = 0
  integer, parameter :: close_CGNS = -1
  integer, parameter :: dump_CGNS = -2
  !
  integer :: CGNS_action
  !
  logical :: request_stop
  logical :: request_solution
  !
 !logical, parameter :: ignore_errors = true
  logical, parameter :: ignore_errors = fals
  !
 !integer, parameter :: test_param = merge( 1 , 2 , ignore_errors )
  !
  integer(c_long) :: current_wtl
  !
  ! Interface for pbs_time_left procedure
  !
#ifdef PBS_ENV
  interface
    subroutine pbs_time_left(seconds_left) bind(c,name="pbs_time_left_")
      use, intrinsic :: iso_c_binding, only : c_long
      integer(c_long) :: seconds_left
    end subroutine pbs_time_left
  end interface
#endif
  !
continue
  !
  ! Initialize the overall MPI session
  !
  call initialize_mpi_environment
  !
  ! Initialize the wall time
  !
  call recon_cpu_time(0)
  !
  ! Get the number of seconds allowed for this simulation
  !
  walltime_expiring = fals
  if (mypnum == 0) then
#ifdef PBS_ENV
    call pbs_time_left(current_wtl)
#else
    current_wtl = huge(int(0,kind=c_long))
#endif
    write (iout,16) current_wtl,ncpu
    max_wall_time = real(current_wtl,kind=wp)
    wall_time_left = max_wall_time
  end if
  !
  ! Broadcast the maximum wall time to all processors
  !
  call mpi_bcast(max_wall_time,1,mpi_flttyp,0,mpi_comm_world,mpierr)
  !
  ! Initialize all the random constants that must be initialized at run time
  !
  call initialize_runtime_parameters
  call memory_pause("base initializations","reading input and grid files")
  !
  ! Read in the input file
  !
  call read_input_file
  call memory_pause("reading input and grid files","partitioning the grid")
  !
  ! Partition the grid if using more than one processor. 
  ! Upon entry:  These arrays define the global grid
  ! Upon return: These arrays define the partition local
  !              to the respective processor along with
  !              connectivities for the partition boundaries
  npart = ncpu
 !npart = 1200
  if (npart > 1) then
    call partition_grid(npart,metis_option,bface,nodes_of_cell, &
                        nodes_of_cell_ptr,xyz_nodes,cell_geom,cell_order)
    call memory_pause("partitioning the grid","getting the solution points")
  else
    call create_serial_cell_map(ncell)
    call memory_pause("getting serial cell map","getting the solution points")
  end if
  !
  ! Get the coordinates of the solution points within each grid cell
  !
  call get_solution_points
  call memory_pause("getting the solution points","getting grid connectivity")
  !
  ! Get the grid connectivity and metrics
  !
  call get_grid_connectivity(bface,cell_geom,nodes_of_cell,nodes_of_cell_ptr)
  call memory_pause("getting grid connectivity","initializing vandermonde matrices")
  !
  ! Compute the remaining quadrature information
  ! (Vandermonde matrix, stiffness matrix, LIFT / correction function, etc.)
  !
  call init_vandermonde_matrices
  call memory_pause("initializing vandermonde matrices","initializing derivative matrices")
  call init_derivative_matrices
  call memory_pause("initializing derivative matrices","initializing correction matrices")
  call init_correction_matrices
  call memory_pause("initializing correction matrices","initializing interpolation matrices")
  call init_interpolation_matrices
  call memory_pause("initializing interpolation matrices","initializing outerpolation matrices")
  call init_outerpolation_matrices
  call memory_pause("initializing outerpolation matrices","initializing projection matrices")
  call init_projection_matrices
  call memory_pause("initializing projection matrices","initializing filter matrices")
  call init_filter_matrices
  call memory_pause("initializing filter matrices","initializing boundary profile stuff")
  call init_bnd_profile_structures
  call memory_pause("initializing boundary profile stuff","initializing metrics")
  !
  ! Get the grid metrics
  !
  call get_metrics
  call memory_pause("initializing metrics","initializing solution")
  !
  ! Output information about each partition
  !
  do i = 0,ncpu-1
    if (mypnum == i) then
      write (iout,10) mypnum, ncell, nnode, nfbnd, nface, n_solpts, local_volume
    end if
    call mpi_barrier(mpi_comm_world,mpierr)
  end do
  call mpi_barrier(mpi_comm_world,mpierr)
  !
  ! Output information about the global grid
  !
  if (mypnum == 0) then
    write (iout,11) n_global_cell, n_global_node, n_global_fbnd, &
                    n_global_solpts, total_volume, &
                    (solpts_edge(k), k=1,size(solpts_edge))
  end if
  call mpi_barrier(mpi_comm_world,mpierr)
  !
  ! Initialize the MMS module if it is being used. 
  !
  if (mms_opt /= 0) call initialize_mms
  !
  ! Initialize the flow variables for all solution points or
  ! read in the restart values if restarting a simulation
  !
  call initialize_flow
  call memory_pause("initializing solution","writing initial CGNS file")
  !
  ! Output the initialized solution
  !
  call write_parallel_cgns(open_CGNS)
  CGNS_action = add_to_CGNS
  call memory_pause("writing initial CGNS file","initializing results file")
  !
  if (mypnum == 0) call write_results_file(1,itcur)
  call memory_pause("initializing results file","initializing generic module")
  !
  if (mypnum == 0) then
    write (iout,19) minval(global_cell_order), maxval(global_cell_order)
  end if
  !
 !call check_flux_point_orderings
  !
  ! Initialize any necessary information for use in the module generic_mod
  !
  call init_generic_mod
  if (num_ts > 0) then
    call memory_pause("initializing generic module","initial interpolation of C.V. to flux points")
  else
    if (dump_memory_usage) then
      call memory_pause("initializing generic module","dumping memory usage")
    else
      call memory_pause("initializing generic module","main iteration loop")
    end if
  end if
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%                                    %%%
  !%%%  Beginning of main iteration loop  %%%
  !%%%                                    %%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  ! Initialize the wall time
  !
  call recon_cpu_time(0)
  !
 !call initialize_limiters
 !call compute_cell_ave
  !
  !
  if (num_ts > 0) then
    !
    call interp_cv_to_flxpts
    call memory_pause("initial interpolation of C.V. to flux points", &
                      "initial evaluation of boundary conditions")
    call get_faceusp_bc
    call memory_pause("initial evaluation of boundary conditions","initial TGV KEDR evaluation")
    !
    if (itestcase == Taylor_Green_Vortex) then
      call compute_TaylorGreen_KEDR(0,zero)
      if (dump_memory_usage) then
        call memory_pause("initial TGV KEDR evaluation","dumping memory usage")
      else
        call memory_pause("initial TGV KEDR evaluation","main iteration loop")
      end if
    end if
    !
  end if
  !
#ifdef PAUSE_BETWEEN_ITERATIONS
  if (mypnum == 0) then
    write (iout,8260)
  end if
  8260 format (//," Pausing before starting the main iteration loop!",/, &
                    " Press any key to continue!",//)
  read (*,*)
#endif
  !
  ! If requested, dump the memory usage of the advanced data types from some of
  ! the modules
  !
  if (dump_memory_usage) then
    !
    if (ncpu > 1) then
      write (memfile,'("memory_usage.",i0,".dat")') mypnum
    else
      memfile = "memory_usage.dat"
    end if
    open (newunit=i,file=trim(adjustl(memfile)),status="replace", &
          action="write",form="formatted",position="rewind",iostat=ierr)
    call io_error(pname,memfile,1,ierr)
    !
    call geovar_memory_usage(i)
    call vandermonde_memory_usage(i)
    call derivatives_memory_usage(i)
    call correction_memory_usage(i)
    call interpolation_memory_usage(i)
    call quadrature_memory_usage(i)
    call projection_memory_usage(i)
    call filter_memory_usage(i)
    call ovar_memory_usage(i)
    call limiters_memory_usage(i)
    call generic_memory_usage(i)
    call plot3d_memory_usage(i)
    call bnd_profiles_memory_usage(i)
    call gmsh_memory_usage(i)
    call metrics_memory_usage(i)
    call flowvar_memory_usage(i)
    call parallel_memory_usage(i)
    !
    close (i,iostat=ierr)
    call io_error(pname,memfile,2,ierr)
    !
    call stop_recon(stop_mpi,pname,"After dumping memory usage")
    !
  end if
  !
  main_loop: do iter = 1, num_ts
    !
    ! Current iteration count including restart iterations
    !
    itcur = itcur + 1
    !
    ! Store a copy of the solution from the previous time step
    !
    uoldsp(:,1:n_solpts) = usp(:,1:n_solpts)
    !
    ! Update the global time step
    !
   !call get_dtsp(dtmin)
   !!
   !! Update elapsed time
   !! NOTE: This is only time accurate when using a global timestep
   !!
   !time = time + dtmin
   !d_t = dtmin
    !
    ! Loop over the Runge-Kutta stages
    !
    rk_loop: do nst = 1,num_stages
      !
      call interp_cv_to_flxpts
      call get_faceusp_bc
      if (nst == 1) call get_dtsp(dtmin)
      !
      ! ######################################################################
      ! ######################################################################
      ! ######################################################################
      ! ######################################################################
      ! ######################################################################
      ! #####                                                            #####
      ! #####  CURRENTLY THE INVISCID OVER-INTEGRATION PATH IS DISABLED  #####
      ! #####  BECAUSE THE ARRAYS interp(:,:)%toCell(:) and projct(:,:)  #####
      ! #####  STILL NEED TO BE SET UP PROPERLY FOR THIS PATH TO WORK!!  #####
      ! #####                                                            #####
      ! #####       ||||||||||||||||||||||||||||||||||||||||||||||       #####
      ! #####       VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV       #####
      ! #####                                                            #####
      ! ######################################################################
      ! ######################################################################
      ! ######################################################################
      ! ######################################################################
      ! ######################################################################
      !
     !if (governing_equations == Euler_Eqns) then
     !  !
     !  call flux_gradient_at_solpts_with_qdtr
     !  call correct_interface_flux_with_qdtr
     !  !
     !else if (governing_equations == NavierStokes_Eqns) then
        !
      if (governing_equations == NavierStokes_Eqns) then
        call discontinuous_gradient_at_solpts
        if (visc_flux_method == visc_flux_BR2) then
          call correct_gradients_br2
        else if (visc_flux_method /= visc_flux_BR2) then
          call correct_gradients_at_solpts
          call interp_dusp_to_flxpts
        end if
        call get_facedusp_bc
      end if
        call flux_gradient_at_solpts_visc
        call correct_interface_flux_visc
        !
     !end if
      !
      ! ######################################################################
      ! ######################################################################
      ! ######################################################################
      ! ######################################################################
      ! ######################################################################
      ! #####                                                            #####
      ! #####       AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA       #####
      ! #####       ||||||||||||||||||||||||||||||||||||||||||||||       #####
      ! #####                                                            #####
      ! #####  CURRENTLY THE INVISCID OVER-INTEGRATION PATH IS DISABLED  #####
      ! #####  BECAUSE THE ARRAYS interp(:,:)%toCell(:) and projct(:,:)  #####
      ! #####  STILL NEED TO BE SET UP PROPERLY FOR THIS PATH TO WORK!!  #####
      ! #####                                                            #####
      ! ######################################################################
      ! ######################################################################
      ! ######################################################################
      ! ######################################################################
      ! ######################################################################
      !
      ! Compute the source terms for the method of manufactured solutions (MMS)
      !
      if (mms_opt /= 0) call mms_src
      !
      ! Compute the source terms for the channel flow
      !
      if (itestcase == Channel_Flow) call channel_source
      !
      ! Compute the final residual for each solution point
      !
      call rresi
      !
      ! Store residual in rssprk(:,:,nst) and update the solutions in usp
      !
      if (iRKtype == Classic_RK) then
        call rupdate_Classic_RK
      else if (iRKtype == TVD2_RK) then
        call rupdate_TVD2_RK
      else if (iRKtype == TVD3_RK) then
        call rupdate_TVD3_RK
      else if (iRKtype == TVD4_RK) then
        call rupdate_TVD4_RK
      else if (iRKtype == CK4_RK) then
        call rupdate_CK4_RK
      end if
      !
    end do rk_loop
    !
    completed_timesteps = iter
    !
    ! Exit the main loop if we have reached the requested final simulation
    ! time when running a time accurate simulation with a global time step.
    !
    if (this_is_final_timestep .and. iter/=num_ts) then
      !
      call rgeterr
      !
      ! Exit the main loop since we have reached
      ! the requested final simulation time
      !
      exit main_loop
      !
    end if
    !
    ! Compute the residual error from the density and output to terminal
    !
    if (mod(iter,iter_out_interval) == 0) then
      !
      call rgeterr
      !
      ! Exit the main loop if the convergence criteria were met
      !
      if (convergence_reached) exit main_loop
      !
      if (iter /= num_ts) then
        if (walltime_expiring) exit main_loop
      end if
      !
    end if
    !
    ! Output the error at the specified interval
    !
    if (mod(iter,results_interval)==0 .and. iter/=num_ts) then
      !
      call recon_cpu_time(1)
      !
      call compute_final_error
      if (mypnum == 0) call write_results_file(0,itcur)
      !
    end if
    !
    ! Every 50 (1000 for NAS) iterations
    !    1. Compute the residual error for this time step and exit
    !       the main loop if the remaining wall time is almost zero
    !    2. Check for a file requesting a premature stop
    !
#ifdef PBS_ENV
    if (mod(iter,1000) == 0) then
#else
    if (mod(iter,50) == 0) then
#endif
      !
      ! Only do the rest if this is not the last time step
      !
      if (iter /= num_ts) then
        !
        ! Have the root processor check for various temporary files
        !
        if (mypnum == 0) then
          !
          ! Check for a stop request
          !
          inquire (file="STOP_RECON",exist=request_stop)
          !
          ! If stop was requested, delete the
          ! STOP_RECON file and exit the main loop
          !
          if (request_stop) then
            !
            write (iout,99)
            !
            ! We dont care about the exit status here since we are
            ! going to be stopping execution after this anyways
            !
            open (iostp,file="STOP_RECON",iostat=ierr)
            close (iostp,status="delete",iostat=ierr)
            !
          end if
          !
          ! Check for a dump of the solution file
          !
          inquire (file="DUMP_SOLUTION",exist=request_solution)
          !
          ! If solution dump was requested, delete the DUMP_SOLUTION file
          !
          if (request_solution) then
            !
            write (iout,97)
            !
            open (iostp,file="DUMP_SOLUTION",iostat=ierr)
            close (iostp,status="delete",iostat=ierr)
            !
          end if
          !
        end if
        !
        ! Have the root processor broadcast request_stop to
        ! the other processors so all processors can exit the
        ! main loop together if the STOP_RECON file was found
        !
        call mpi_bcast(request_stop,1,mpi_logical,0,mpi_comm_world,mpierr)
        if (request_stop) exit main_loop
        !
        ! Have the root processor broadcast request_solution to
        ! the other processors
        !
        call mpi_bcast(request_solution,1,mpi_logical,0,mpi_comm_world,mpierr)
        !
        ! If a solution dump was requested, dump the solution and then let
        ! the CGNS subroutine know that the file needs to be reopened.
        !
        if (request_solution) then
          call write_parallel_cgns(dump_CGNS)
          CGNS_action = reopen_CGNS
        end if
        !
      end if
      !
    end if
    !
    ! Output the restart file at specified step interval
    !
    if (restart_interval > 0) then
      if (mod(iter,restart_interval)==0 .and. iter/=num_ts) then
        call write_restart_file(itcur)
      end if
    end if
    !
    ! Output the solution at the specified interval
    !
    if (output_interval > 0) then
      if (mod(iter,output_interval)==0 .and. iter/=num_ts) then
        call write_parallel_cgns(CGNS_action)
      end if
    end if
    !
    if (itestcase == Taylor_Green_Vortex) then
      if (write_TaylorGreen_8s_solution) then
        call TaylorGreen_8s_vorticity_solution
      end if
      if (write_TaylorGreen_full_solution) then
        call write_parallel_cgns(CGNS_action)
        write_TaylorGreen_full_solution = fals
      end if
    end if
    !
    ! If solving the Taylor-Green Vortex problem, compute the various integrated
    ! quantities needed to get the different the kinetic energy dissipation rates
    !
    if (itestcase == Taylor_Green_Vortex .and. iter/=num_ts) then
      call compute_TaylorGreen_KEDR(iter,dtmin)
    end if
    !
#ifdef PAUSE_BETWEEN_ITERATIONS
    if (iter /= num_ts) then
      if (mypnum == 0) then
        write (iout,8261) iter,num_ts
      end if
      8261 format (//," Pausing after completing ",i0," of ",i0, &
                        " trips through the main iteration loop!",/, &
                        " Press any key to continue!",//)
      read (*,*)
    end if
#endif
    !
  end do main_loop
  !
#ifdef PAUSE_BETWEEN_ITERATIONS
  if (mypnum == 0) then
    write (iout,8262)
  end if
  8262 format (//," Pausing after finishing the main iteration loop!",/, &
                    " Press any key to continue!",//)
  read (*,*)
#endif
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%                              %%%
  !%%%  End of main iteration loop  %%%
  !%%%                              %%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  ! Get the elapsed wall time
  !
  call recon_cpu_time(1)
  !
  ! If solving the Taylor-Green Vortex problem,
  ! finalize the kinetic energy dissipation rates
  !
  if (num_ts > 0) then
    if (itestcase == Taylor_Green_Vortex) then
      call compute_TaylorGreen_KEDR(-completed_timesteps,dtmin)
    end if
  end if
  !
  !.. get min and max entropy value and location
  if (.not. ignore_errors) call getmaxent
  !
  call write_parallel_cgns(close_CGNS)
  !
  call write_restart_file(-1)
  !
  if (.not. ignore_errors) then
    !
    call compute_final_error
    !
    call compute_bnd_profiles
    !
    if (mypnum == 0) then
      !
      write (iout,40) itcur,time
      call recon_cpu_time(-1)
      !
      call write_final_error
      call write_results_file(-1,itcur)
      !
    end if
    !
  end if
  !
  flush (iout)
  call mpi_barrier(mpi_comm_world,mpierr)
  !
  if (mypnum == 0) write (iout,50)
  !
  ! Finalize the overall MPI session and stop execution
  !
  call finalize_mpi_environment
  !
  !
  !
  10 format (/,' For Processor ',i0,/,'------------------',/, &
               '  ncell = ',i0,'  nnode = ',i0,'  nfbnd = ',i0, &
               '  nface = ',i0,'  nDOF = ',i0,/, &
               '  Partition Volume = ',f13.6)
  11 format (/,' For the Global Grid ',/,'---------------------',/, &
               '  ncell = ',i0,'  nnode = ',i0,'  nfbnd = ',i0,'  nDOF = ',i0,/, &
               '  Total Grid Volume = ',f13.6,/, &
               '  1D Solution Points = ',6f9.5,:,*(/,23x,6f9.5))
  16 format (/," Allocated Wall Time = ",i0," s. running on ",i0," CPUs",/)
  !
  19 format (/,' Lowest  Cell Order = ',i4,/, &
               ' Highest Cell Order = ',i4,/)
  !
  40 format (//' Final iteration count   = ',i12,/, &
               ' Simulation elapsed time = ',f12.3)
  41 format (/,"TOTAL SIMULATION WALL TIME WAS ",f18.6," s.",/,&
               "TOTAL COMPUTATIONAL COST WAS   ",f18.6," work units.")
  !
  50 format (/'Done'/)
  !
  97 format (/,"*******************************************************",/, &
               "  FOUND FILE REQUESTING A DUMP OF THE SOLUTION FILE! ",/, &
               "*******************************************************",/)
  98 format (/,"*******************************************************",/, &
               " VORTEX HAS APPEARED TO COLLAPSE!",/, &
               " THE SIMULATION WAS RUN ANOTHER 25 PERIODS", &
               " PAST THE POINT OF COLLAPSE!",/, &
               " RECON IS NOW STOPPING EXECUTION AND", &
               " DUMPING THE FINAL RESULTS!",/, &
               "*******************************************************",/)
  99 format (/,"*******************************************************",/, &
               " FOUND FILE REQUESTING TO STOP THE EXECUTION OF RECON! ",/, &
               "   STOPPING EXECUTION AND DUMPING THE FINAL RESULTS!",/, &
               "*******************************************************",/)
  !
end program Euler_tri
!
!###############################################################################
!
