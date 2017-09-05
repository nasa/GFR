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
  !
  use ovar, only : metis_option,convergence_reached
  use ovar, only : total_volume,local_volume,itestcase
  use ovar, only : num_rk_stages,Runge_Kutta_Scheme
  use ovar, only : time,time_ref,itcur,num_timesteps
  use ovar, only : governing_equations,visc_flux_method
  use ovar, only : results_interval,iter_out_interval
  use ovar, only : restart_interval,output_interval
  use ovar, only : dump_memory_usage,mms_opt
  use ovar, only : recon_cpu_time,this_is_final_timestep
  use ovar, only : write_TaylorGreen_8s_solution
  use ovar, only : write_TaylorGreen_full_solution
  use ovar, only : continuous_output,completely_disable_cgns
  use ovar, only : check_flux_point_coordinates
  use ovar, only : cfl,cfl_beg,cfl_end,cfl_cycles
  use ovar, only : adjust_cfl_cycles
  use ovar, only : output_time_averaging,time_scaling_factor
  use ovar, only : time_average_restart_files,ave_start_time
  use ovar, only : Filtering_Interval,Limiting_Interval
  !
  use input_namelist_mod, only : cfl_adjustment
  use input_namelist_mod, only : restart_interupt
  !
  use flowvar, only : usp,uoldsp
  use flowvar, only : faceusp,facedusp
  !
  use quadrature_mod, only : solpts_edge
  !
  use io_mod, only : read_input_file
  use io_mod, only : write_final_error
  use io_mod, only : write_parallel_cgns
  use io_mod, only : write_residual_file
  use io_mod, only : check_face_flux_point_orderings
  use io_mod, only : check_cell_flux_point_orderings
  use io_mod, only : TaylorGreen_8s_vorticity_solution
  !
  use restart_mod, only : write_restart_file
  !
  use mappings_mod, only : get_solution_points
  !
  use connectivity_mod, only : get_grid_connectivity
  !
  use metrics_mod, only : get_metrics
  !
  use vandermonde_mod, only : init_vandermonde_matrices
  !
  use derivatives_mod, only : init_derivative_matrices
  !
  use correction_mod, only : init_correction_matrices
  !
  use interpolation_mod, only : init_interpolation_matrices
  use interpolation_mod, only : init_outerpolation_matrices
  !
  use projection_mod, only : init_projection_matrices
  !
  use filter_mod, only : init_filter_matrices
  !
  use initialization_mod, only : initialize_flow
  !
  use module_bc, only : get_faceusp_bc
  use module_bc, only : get_facedusp_bc
  !
  use generic_mod, only : init_generic_mod
  use generic_mod, only : write_iter_stats
  use generic_mod, only : getmaxent
  use generic_mod, only : compute_residual_error
  use generic_mod, only : compute_TaylorGreen_KEDR
  !
  use time_mod, only : get_dtsp
  use time_mod, only : rresi
  use time_mod, only : rupdate_Classic_RK
  use time_mod, only : rupdate_CK4_RK
  use time_mod, only : rupdate_TVD2_RK
  use time_mod, only : rupdate_TVD3_RK
  use time_mod, only : rupdate_TVD4_RK
  use time_mod, only : update_time_ave
  !
  use flux_mod, only : interp_cv_to_flxpts
  use flux_mod, only : interp_dusp_to_flxpts
  use flux_mod, only : flux_gradient_at_solpts_visc
  use flux_mod, only : flux_gradient_at_solpts_visc_gpu
  use flux_mod, only : flux_gradient_at_solpts_visc_gpu1
  use flux_mod, only : correct_interface_flux_visc
  use flux_mod, only : discontinuous_gradient_at_solpts
  use flux_mod, only : correct_gradients_at_solpts
  use flux_mod, only : correct_gradients_br2
  use flux_mod, only : new_flux_gradient_at_solpts
  use flux_mod, only : get_invs_common_fluxes_interior
  use flux_mod, only : get_visc_common_fluxes_interior
  use flux_mod, only : get_invs_common_fluxes_bc
  use flux_mod, only : get_visc_common_fluxes_bc
  !
  use parallel_mod, only : partition_grid
  use parallel_mod, only : create_serial_cell_map
  use parallel_mod, only : exchange_faceusp
  use parallel_mod, only : exchange_facedusp
  use parallel_mod, only : wait_faceusp
  use parallel_mod, only : wait_facedusp
  !
  use module_limiters, only : initialize_limiters
  use module_limiters, only : compute_cell_ave
  !
  use mms_mod, only : mms_src
  use mms_mod, only : initialize_mms
  !
  use channel_mod, only : channel_source
  !
  use bnd_profiles_mod, only : init_bnd_profile_structures
  use bnd_profiles_mod, only : compute_bnd_profiles
  !
  use postprocessor_mod, only : new_write_parallel_cgns
  !
  use pbs_mod, only : pbs_remaining_walltime
  !
  use cpu_info_mod, only : get_cpu_locality
  !
  use averaging_mod, only : average_restart_files
  !
  use, intrinsic :: iso_c_binding,   only : c_long
  !
  implicit none
  !
  character(len=*), parameter :: pname = "Euler_tri(main)"
  !
  integer  :: iter,nst,i,ierr,k,npart
  integer  :: completed_timesteps
  real(wp) :: dtmin,out_time
  !
  character(len=50) :: memfile
  !
  integer :: CGNS_action
  !
  logical(lk)  :: exit_main_loop
  logical(lk)  :: walltime_expiring
  logical(ldk) :: request_stop
  logical(ldk) :: request_restart
  logical(ldk) :: request_solution
  logical(ldk) :: request_cfl_adjustment
  logical(ldk) :: request_restart_interupt
  !
  logical(lk)  :: Apply_Filter
  logical(lk)  :: Apply_Limiter
  !
 !logical(lk), parameter :: ignore_errors = true
  logical(lk), parameter :: ignore_errors = fals
  !
  character(len=300) :: overflow_error = ""
  character(len=300) :: divide_by_zero_error = ""
  character(len=300) :: invalid_error = ""
  !
  real(wp) :: cfl_tmp(1:4)
  !
continue
  !
  request_restart = .false.
  request_solution = .false.
  !
  ! Initialize the overall MPI session
  !
  call initialize_mpi_environment
  !
  call get_cpu_locality
  !
  ! Initialize the wall time
  !
  call recon_cpu_time(0)
  !
  ! Get the number of seconds allowed for this simulation
  !
  call pbs_remaining_walltime(walltime_expiring)
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
  call get_grid_connectivity(bface,cell_geom,cell_order,xyz_nodes, &
                             nodes_of_cell,nodes_of_cell_ptr)
  call memory_pause("getting grid connectivity", &
                    "initializing vandermonde matrices")
  !
  ! Compute the remaining quadrature information
  ! (Vandermonde matrix, stiffness matrix, LIFT / correction function, etc.)
  !
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  call init_vandermonde_matrices
  call memory_pause("initializing vandermonde matrices", &
                    "initializing derivative matrices")
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  call init_derivative_matrices
  call memory_pause("initializing derivative matrices", &
                    "initializing correction matrices")
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  call init_correction_matrices
  call memory_pause("initializing correction matrices", &
                    "initializing interpolation matrices")
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  call init_interpolation_matrices
  call memory_pause("initializing interpolation matrices", &
                    "initializing outerpolation matrices")
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  call init_outerpolation_matrices
  call memory_pause("initializing outerpolation matrices", &
                    "initializing projection matrices")
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  call init_projection_matrices
  call memory_pause("initializing projection matrices", &
                    "initializing filter matrices")
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  call init_filter_matrices
  call memory_pause("initializing filter matrices", &
                    "initializing boundary profile stuff")
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  call init_bnd_profile_structures
  call memory_pause("initializing boundary profile stuff", &
                    "initializing metrics")
  !
  ! Get the grid metrics
  !
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  call get_metrics
  call memory_pause("initializing metrics","initializing solution")
  !
  ! Uncomment the call to check_flux_point_orderings if you
  ! suspect that there is an error between the orderings of
  ! the flux points from the two cells on an interface
  !
  if (check_flux_point_coordinates) then
    call check_face_flux_point_orderings
    call check_cell_flux_point_orderings
  end if
  !
  ! Output information about the global grid
  !
  if (mypnum == glb_root) then
    write (iout,11) n_global_cell, n_global_node, n_global_fbnd, &
                    n_global_solpts, total_volume, &
                    (solpts_edge(k), k=1,size(solpts_edge))
  end if
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  ! If requested, dump the memory usage of the advanced data types from some of
  ! the modules
  !
  if (dump_memory_usage) then
    call report_memory_usage(n=1,request_stop=fals)
  end if
  !
  ! Initialize the MMS module if it is being used. 
  !
  if (mms_opt /= 0) call initialize_mms
  !
  ! Initialize the flow variables for all solution points or
  ! read in the restart values if restarting a simulation
  !
  call initialize_flow
  call memory_pause("initializing solution","initializing generic module")
  !
  ! Time average the list of restart files
  !
  if (time_average_restart_files) then
    call average_restart_files
  end if
  !
  ! If requested, dump the memory usage of the advanced data types from some of
  ! the modules
  !
  if (dump_memory_usage) then
    call report_memory_usage(n=2,request_stop=true)
  end if
  !
  if (mypnum == glb_root) then
    write (iout,19) minval(global_cell_order), maxval(global_cell_order)
  end if
  !
  ! Initialize any necessary information for use in the module generic_mod
  !
  call init_generic_mod
  if (num_timesteps > 0) then
    call memory_pause("initializing generic module", &
                      "initial interpolation of C.V. to flux points")
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
  call initialize_limiters
  !
  if (num_timesteps > 0) then
    !
    call interp_cv_to_flxpts
    call memory_pause("initial interpolation of C.V. to flux points", &
                      "exchanging C.V. between communication faces")
    if (ncpu > 1) call exchange_faceusp(faceusp)
    call memory_pause("exchanging C.V. between communication faces", &
                      "initial evaluation of boundary conditions")
    call get_faceusp_bc
    call memory_pause("initial evaluation of boundary conditions", &
                      "waiting for C.V. exchanges to complete")
    if (ncpu > 1) call wait_faceusp(faceusp)
    call memory_pause("waiting for C.V. exchanges to complete", &
                      "initial TGV KEDR evaluation")
    !
    if (itestcase == Taylor_Green_Vortex) then
      call compute_TaylorGreen_KEDR(0,zero)
      if (dump_memory_usage) then
        call memory_pause("initial TGV KEDR evaluation","dumping memory usage")
      else
        call memory_pause("initial TGV KEDR evaluation", &
                          "writing initial CGNS file")
      end if
    end if
    !
  end if
  !
#ifdef PAUSE_BETWEEN_ITERATIONS
  if (mypnum == glb_root) then
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
    call report_memory_usage(n=3,request_stop=true)
  end if
  !
  ! Output the initialized solution
  !
  if (continuous_output) then
    call new_write_parallel_cgns
  else
    call write_parallel_cgns(open_CGNS)
  end if
  CGNS_action = add_to_CGNS
  call memory_pause("writing initial CGNS file","main iteration loop")
  !
  main_loop: do iter = 1, num_timesteps
    !
    ! Current iteration count including restart iterations
    !
    itcur = itcur + 1
    !$ACC UPDATE DEVICE(itcur)
    !
    ! Store a copy of the solution from the previous time step
    !
    uoldsp(:,1:n_solpts) = usp(:,1:n_solpts)
    !
    ! Loop over the Runge-Kutta stages
    !
    rk_loop: do nst = 1,num_rk_stages
      !
      Apply_Filter = (Filtering_Interval == 0)
      Apply_Limiter = (Limiting_Interval == 0)
      if (nst == num_rk_stages) then
        if (Filtering_Interval > 0) then
          Apply_Filter = (mod(iter,Filtering_Interval) == 0)
        end if
        if (Limiting_Interval > 0) then
          Apply_Limiter = (mod(iter,Limiting_Interval) == 0)
        end if
      end if
#define NEW_FLUX_REARRANGEMENT
#ifdef NEW_FLUX_REARRANGEMENT
      call interp_cv_to_flxpts
      if (ncpu > 1) call exchange_faceusp(faceusp)
      if (governing_equations == NavierStokes_Eqns) then
        call discontinuous_gradient_at_solpts
      end if
      call get_faceusp_bc
      call get_invs_common_fluxes_interior
      !
      if (ncpu > 1) call wait_faceusp(faceusp)
      !
      if (governing_equations == NavierStokes_Eqns) then
        if (visc_flux_method == visc_flux_BR2) then
          call correct_gradients_br2
        else
          call correct_gradients_at_solpts
          call interp_dusp_to_flxpts
        end if
        if (ncpu > 1) call exchange_facedusp(facedusp)
        call get_facedusp_bc
        call get_visc_common_fluxes_interior
      end if
      !
      call get_invs_common_fluxes_bc
      !
      if (governing_equations == NavierStokes_Eqns) then
        if (ncpu > 1) call wait_facedusp(facedusp)
        call get_visc_common_fluxes_bc
      end if
      !
      call new_flux_gradient_at_solpts
#else
      !
      call interp_cv_to_flxpts
      !
      if (ncpu > 1) call exchange_faceusp(faceusp)
      !
      call get_faceusp_bc
      !
      if (ncpu > 1) call wait_faceusp(faceusp)
      !
     !if (nst == 1) call get_dtsp(dtmin)
      !
      if (governing_equations == NavierStokes_Eqns) then
        !
        call discontinuous_gradient_at_solpts
        !
        if (visc_flux_method == visc_flux_BR2) then
          call correct_gradients_br2
        else if (visc_flux_method /= visc_flux_BR2) then
          call correct_gradients_at_solpts
          call interp_dusp_to_flxpts
        end if
        !
        if (ncpu > 1) call exchange_facedusp(facedusp)
        !
        call get_facedusp_bc
        !
        if (ncpu > 1) call wait_facedusp(facedusp)
        !
      end if
      !
     !call flux_gradient_at_solpts_visc(iter,nst)
     !call flux_gradient_at_solpts_visc_gpu(iter,nst)
      call flux_gradient_at_solpts_visc_gpu1(iter,nst)
      call correct_interface_flux_visc(iter,nst)
#endif
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
      if (nst == 1) call get_dtsp(dtmin)
      !
      ! Store residual in rssprk(:,:,nst) and update the solutions in usp
      !
      if (Runge_Kutta_Scheme == Classic_RK) then
        call rupdate_Classic_RK(nst,Apply_Filter,Apply_Limiter)
      else if (Runge_Kutta_Scheme == TVD2_RK) then
        call rupdate_TVD2_RK(nst,Apply_Filter,Apply_Limiter)
      else if (Runge_Kutta_Scheme == TVD3_RK) then
        call rupdate_TVD3_RK(nst,Apply_Filter,Apply_Limiter)
      else if (Runge_Kutta_Scheme == TVD4_RK) then
        call rupdate_TVD4_RK(nst,Apply_Filter,Apply_Limiter)
      else if (Runge_Kutta_Scheme == CK4_RK) then
        call rupdate_CK4_RK(nst,Apply_Filter,Apply_Limiter)
      end if
      !
    end do rk_loop
    !
    completed_timesteps = iter
    !
    ! Update the current time now that the Runge-Kutta loop is finished
    !
    time = time + dtmin
    !$ACC UPDATE DEVICE(time)
    !
   !call check_for_bad_number(__LINE__)
    !
    ! Update the time averaged variables
    !
    if (output_time_averaging) then
      call update_time_ave
    end if
    !
    ! Compute and output the residual error at the specified interval
    !
    if ( (mod(iter,iter_out_interval) == 0 .or. this_is_final_timestep) .or. &
         (mod(itcur,results_interval) == 0 .and. iter/=num_timesteps) ) then
      !
      ! Compute the residual error
      !
      call compute_residual_error
      !
      if (mod(iter,iter_out_interval) == 0 .or. this_is_final_timestep) then
        !
        ! Output the density residual error to stdout
        !
        call write_iter_stats(walltime_expiring,exit_main_loop)
        !
        ! Exit the main loop if the convergence criteria were met
        !
        if (exit_main_loop) exit main_loop
        !
      end if
      !
      ! Output the residual error at the specified interval
      !
      if (mod(itcur,results_interval)==0 .and. iter/=num_timesteps) then
        !
        call recon_cpu_time(1)
        !
        if (mypnum == glb_root) call write_residual_file(0,itcur)
        !
      end if
      !
    end if
    !
    call check_for_bad_number(__LINE__)
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
      if (iter /= num_timesteps) then
        !
        ! Have the root processor check for various temporary files
        !
        if (mypnum == glb_root) then
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
            open (newunit=i,file="STOP_RECON",iostat=ierr)
            close (i,status="delete",iostat=ierr)
            !
          end if
          !
          ! Check for a dump of a restart file
          !
          inquire (file="DUMP_RESTART",exist=request_restart)
          !
          ! If restart dump was requested, delete the DUMP_RESTART file
          !
          if (request_restart) then
            !
            write (iout,96)
            !
            open (newunit=i,file="DUMP_RESTART",iostat=ierr)
            close (i,status="delete",iostat=ierr)
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
            open (newunit=i,file="DUMP_SOLUTION",iostat=ierr)
            close (i,status="delete",iostat=ierr)
            !
          end if
          !
          ! Check for a CFL adjustment file
          !
          inquire (file="ADJUST_CFL",exist=request_cfl_adjustment)
          !
          ! If a CFL adjustment was requested, delete the ADJUST_CFL file
          !
          if (request_cfl_adjustment) then
            !
            write (iout,95)
            !
            open (newunit=i,file="ADJUST_CFL",iostat=ierr)
            read (i,nml=cfl_adjustment)
            close (i,status="delete",iostat=ierr)
            !
            cfl_tmp = [cfl,cfl_beg,cfl_end,real(cfl_cycles,kind=wp)]
            !
          end if
         !!
         !! Check for a file requesting an interuption of the current
         !! simulation to restart it with a new restart file
         !!
         !inquire (file="RESTART_INTERUPT",exist=request_restart_interupt)
         !!
         !!
         !! Delete the RESTART_INTERUPT file if restart interupt was requested
         !!
         !if (request_restart_interupt) then
         !  !
         !  write (iout,96)
         !  !
         !  open (newunit=i,file="RESTART_INTERUPT",iostat=ierr)
         !  close (i,status="delete",iostat=ierr)
         !  !
         !end if
          !
        end if
        !
        ! Have the root processor broadcast request_stop to
        ! the other processors so all processors can exit the
        ! main loop together if the STOP_RECON file was found
        !
        call mpi_bcast(request_stop,1_int_mpi,MPI_LOGICAL, &
                       0_int_mpi,MPI_COMM_WORLD,mpierr)
        if (request_stop) exit main_loop
        !
        ! Have the root processor broadcast request_restart to
        ! the other processors
        !
        call mpi_bcast(request_restart,1_int_mpi,MPI_LOGICAL, &
                       0_int_mpi,MPI_COMM_WORLD,mpierr)
        !
        ! If a restart dump was requested, write a restart file
        !
        if (request_restart) then
          call write_restart_file(itcur)
        end if
        !
        ! Have the root processor broadcast request_solution to
        ! the other processors
        !
        call mpi_bcast(request_solution,1_int_mpi,MPI_LOGICAL, &
                       0_int_mpi,MPI_COMM_WORLD,mpierr)
        !
        ! If a solution dump was requested, dump the solution and then let
        ! the CGNS subroutine know that the file needs to be reopened.
        !
        if (request_solution) then
          if (continuous_output) then
            call new_write_parallel_cgns
          else
            call write_parallel_cgns(dump_CGNS)
          end if
          CGNS_action = reopen_CGNS
        end if
        !
        call mpi_bcast(request_cfl_adjustment,1_int_mpi,MPI_LOGICAL, &
                       0_int_mpi,MPI_COMM_WORLD,mpierr)
        !
        if (request_cfl_adjustment) then
          !
          call mpi_bcast(cfl_tmp,size(cfl_tmp,kind=int_mpi),mpi_flttyp, &
                         0_int_mpi,MPI_COMM_WORLD,mpierr)
          !
          cfl = cfl_tmp(1)
          cfl_beg = cfl_tmp(2)
          cfl_end = cfl_tmp(3)
          cfl_cycles = nint(cfl_tmp(4))
          !
          call adjust_cfl_cycles
          !
        end if
        !
       !call mpi_bcast(request_restart_interupt,1_int_mpi,MPI_LOGICAL, &
       !               0_int_mpi,MPI_COMM_WORLD,mpierr)
        !
      end if
      !
    end if
    !
    ! Output the restart file at specified step interval
    !
    if (restart_interval > 0 .and. .not.request_restart) then
      if (mod(itcur,restart_interval)==0 .and. iter/=num_timesteps) then
        call write_restart_file(itcur)
      end if
    end if
    !
    ! Output the solution at the specified interval
    !
    if (output_interval > 0 .and. .not.request_solution) then
      if (mod(itcur,output_interval)==0 .and. iter/=num_timesteps) then
        if (continuous_output) then
          call new_write_parallel_cgns
        else
          call write_parallel_cgns(CGNS_action)
        end if
      end if
    end if
    !
    if (itestcase == Taylor_Green_Vortex) then
      if (write_TaylorGreen_8s_solution) then
        if (completely_disable_cgns) then
          write_TaylorGreen_8s_solution = fals
        else
          call TaylorGreen_8s_vorticity_solution
        end if
      end if
      if (write_TaylorGreen_full_solution) then
        if (continuous_output) then
          call new_write_parallel_cgns
        else
          call write_parallel_cgns(CGNS_action)
        end if
        write_TaylorGreen_full_solution = fals
      end if
    end if
    !
    ! If solving the Taylor-Green Vortex problem, compute
    ! the various integrated quantities needed to get the
    ! different the kinetic energy dissipation rates
    !
    if (itestcase == Taylor_Green_Vortex .and. iter/=num_timesteps) then
      call compute_TaylorGreen_KEDR(iter,dtmin)
    end if
    !
#ifdef PAUSE_BETWEEN_ITERATIONS
    if (iter /= num_timesteps) then
      if (mypnum == glb_root) then
        write (iout,8261) iter,num_timesteps
      end if
      8261 format (//," Pausing after completing ",i0," of ",i0, &
                        " trips through the main iteration loop!",/, &
                        " Press any key to continue!",//)
      read (*,*)
    end if
#endif
    !
    ! Reset the request options
    !
    request_restart = .false.
    request_solution = .false.
    !
  end do main_loop
  !
#ifdef PAUSE_BETWEEN_ITERATIONS
  if (mypnum == glb_root) then
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
  if (num_timesteps > 0) then
    if (itestcase == Taylor_Green_Vortex) then
      call compute_TaylorGreen_KEDR(-completed_timesteps,dtmin)
    end if
  end if
  !
  call write_restart_file(-1)
  !
  if (continuous_output) then
    call new_write_parallel_cgns
  else
    call write_parallel_cgns(close_CGNS)
  end if
  !
  if (.not. ignore_errors) then
    !
    call compute_residual_error
    !
    ! Get min and max entropy value and location
    call getmaxent
    !
    call compute_bnd_profiles
    !
    if (mypnum == glb_root) then
      !
      out_time = time*time_ref*time_scaling_factor
      if (output_time_averaging) then
        out_time = out_time - ave_start_time*time_ref*time_scaling_factor
      end if
      write (iout,40) itcur,out_time
      call recon_cpu_time(-1)
      !
      call write_final_error
      call write_residual_file(-1,itcur)
      !
    end if
    !
  end if
  !
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  if (mypnum == glb_root) write (iout,50)
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
               '  ncell = ',i0,'  nnode = ',i0, &
               '  nfbnd = ',i0,'  nDOF = ',i0,/, &
               '  Total Grid Volume = ',f13.6,/, &
               '  1D Solution Points = ',100(6f9.5,:,/,23x))
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
  95 format (/,"*******************************************************",/, &
               "    FOUND FILE REQUESTING AN ADJUSTMENT TO THE CFL! ",/, &
               "*******************************************************",/)
  96 format (/,"*******************************************************",/, &
               "  FOUND FILE REQUESTING A DUMP OF A RESTART FILE! ",/, &
               "*******************************************************",/)
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
contains
!
!###############################################################################
!
subroutine check_for_bad_number(line_num)
  !
#ifndef SPECIAL_FOR_GCC
  use, intrinsic :: ieee_exceptions, only : ieee_usual
  use, intrinsic :: ieee_exceptions, only : ieee_get_flag
#endif
  !
  !.. Formal Arguments ..
  integer, intent(in) :: line_num
  !
#ifndef SPECIAL_FOR_GCC
  !.. Local Arrays ..
  integer      :: bad_number_found(1:size(ieee_usual))
  logical(ldk) :: signaling_exceptions(1:size(ieee_usual))
  !
  character(len=*), parameter :: pname = "check_for_bad_number"
  !
continue
  !
  ! First check to make sure that no floating-point exceptions are signaling
  !
  call ieee_get_flag(ieee_usual,signaling_exceptions)
  !
  bad_number_found = merge(1,0,signaling_exceptions)
  !
  call mpi_allreduce(MPI_IN_PLACE,bad_number_found, &
                     size(bad_number_found,kind=int_mpi), &
                     mpi_inttyp,MPI_MAX,MPI_COMM_WORLD,mpierr)
  !
  if (sum(bad_number_found) > 0) then
    !
    ! A floating-point exception was detected. Dump solution and restart
    ! files for the solution at the previous time step and abort the
    ! simulation
    !
    ! Copy the solution from the previous time step back into
    ! usp so we can dump the solution and restart files before
    ! aborting the simulation
    !
    usp(:,1:n_solpts) = uoldsp(:,1:n_solpts)
    !
    if (continuous_output) then
      call new_write_parallel_cgns
    else
      call write_parallel_cgns(-3)
    end if
    !
    call write_restart_file(-2)
    !
    if (bad_number_found(1) > 0) write (overflow_error,101)
    if (bad_number_found(2) > 0) write (divide_by_zero_error,102)
    if (bad_number_found(3) > 0) write (invalid_error,103)
    write (error_message,100) itcur, &
                              trim(overflow_error), &
                              trim(divide_by_zero_error), &
                              trim(invalid_error)
    !
    call stop_recon(stop_mpi,pname,line_num,__FILE__,error_message)
    !
  end if
  !
  100 format ("A SIGNALING FLOATING POINT EXCEPTION WAS ", &
              "DETECTED AFTER UPDATING THE SOLUTION AFTER ", &
              "TIME STEP # ",i0," !!! DUMPING SOLUTION AND ", &
              "RESTART FILES FOR THE PREVIOUS TIME STEP!!!", &
              3(2x,a,:))
  !
  101 format ("OVERFLOW EXCEPTION DETECTED!!!")
  102 format ("DIVIDE-BY-ZERO EXCEPTION DETECTED!!!")
  103 format ("INVALID NUMBER (NaN) EXCEPTION DETECTED!!!")
#endif
  !
end subroutine check_for_bad_number
!
!###############################################################################
!
subroutine report_memory_usage(n,request_stop)
  !
  !.. Use Statements ..
  use geovar,            only : geovar_memory_usage
  use ovar,              only : ovar_memory_usage
  use flowvar,           only : flowvar_memory_usage
  use quadrature_mod,    only : quadrature_memory_usage
  use metrics_mod,       only : metrics_memory_usage
  use vandermonde_mod,   only : vandermonde_memory_usage
  use derivatives_mod,   only : derivatives_memory_usage
  use correction_mod,    only : correction_memory_usage
  use interpolation_mod, only : interpolation_memory_usage
  use projection_mod,    only : projection_memory_usage
  use filter_mod,        only : filter_memory_usage
  use generic_mod,       only : generic_memory_usage
  use parallel_mod,      only : parallel_memory_usage
  use module_limiters,   only : limiters_memory_usage
  use bnd_profiles_mod,  only : bnd_profiles_memory_usage
  use module_plot3d,     only : plot3d_memory_usage
  use module_gmsh,       only : gmsh_memory_usage
  use module_cgns,       only : cgns_memory_usage
  !
  !.. Optional Arguments ..
  integer,     optional, intent(in) :: n
  logical(lk), optional, intent(in) :: request_stop
  !
  !.. Local Scalars ..
  integer :: i,ierr,ncpu_digits
  !
  character(len=100) :: memfile
  character(len=100) :: frmt
  character(len=100) :: char_num
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "report_memory_usage"
  logical(lk),      parameter :: skip_output = true
  !
continue
  !
  frmt = "('memory_usage'"
  !
  if (ncpu > 1) then
    !
    write (char_num,1) ncpu
    ncpu_digits = max(1,len_trim(adjustl(char_num)))
    !
    write (frmt,10) trim(adjustl(frmt)),ncpu_digits,ncpu_digits
    !
  end if
  !
  if (present(n)) then
    write (frmt,11) trim(adjustl(frmt)),n
  end if
  !
  frmt = trim(adjustl(frmt)) // ",'.dat')"
  !
  if (ncpu > 1) then
    write (memfile,frmt) mypnum
  else
    write (memfile,frmt)
  end if
  !
  open (newunit=i,file=trim(adjustl(memfile)),status="replace", &
        action="write",form="formatted",position="rewind", &
        iostat=ierr,iomsg=error_message)
  call io_error(pname,memfile,1,__LINE__,__FILE__,ierr,error_message)
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
  call cgns_memory_usage(i)
  !
  close (i,iostat=ierr,iomsg=error_message)
  call io_error(pname,memfile,2,__LINE__,__FILE__,ierr,error_message)
  !
  if (present(request_stop)) then
    if (request_stop) then
      write (error_message,20)
      call stop_recon(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
  end if
  !
  ! Format Statements
  !
  1  format (i0)
  10 format (a,",'.CPU-',i",i0,".",i0)
  11 format (a,",'.CP-",i0,"'")
  20 format ("After dumping memory usage")
  !
end subroutine report_memory_usage
!
end program Euler_tri
!
!###############################################################################
!
