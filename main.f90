program GFR_main
  !
  use module_kind_types
  !
  use geovar, only : ncell,n_global_cell
  use geovar, only : nnode,n_global_node
  use geovar, only : nfbnd,n_global_fbnd
  use geovar, only : nr,n_solpts
  use geovar, only : n_global_solpts
  use geovar, only : bface,xyz_nodes,nface
  use geovar, only : cell_geom,cell_order,global_cell_order
  use geovar, only : nodes_of_cell_ptr,nodes_of_cell
  use geovar, only : move_face_to_device
  use geovar, only : move_cell_to_device
  use geovar, only : move_fp_to_device
  use geovar, only : move_bcgrps_to_device
  use geovar, only : move_cellmap_to_device
  use geovar, only : cell
  use geovar, only : n_flxpts,n_bnd_flxpts,n_wall_flxpts
  use geovar, only : n_qdtpts,n_flx_qdtpts,n_bnd_qdtpts,n_wall_qdtpts
  !
  use eqn_idx, only : nq
  !
  use ovar, only : metis_option,convergence_reached
  use ovar, only : total_volume,local_volume,itestcase
  use ovar, only : num_rk_stages,Runge_Kutta_Scheme
  use ovar, only : time,time_ref,itcur,num_timesteps,rl2mx
  use ovar, only : governing_equations,visc_flux_method
  use ovar, only : results_interval,iter_out_interval
  use ovar, only : restart_interval,output_interval
  use ovar, only : is_at_output_frequency
  use ovar, only : dump_memory_usage,mms_opt
  use ovar, only : gfr_cpu_time,this_is_final_timestep
  use ovar, only : write_TaylorGreen_8s_solution
  use ovar, only : write_TaylorGreen_full_solution
  use ovar, only : continuous_output,completely_disable_cgns
  use ovar, only : check_flux_point_coordinates
  use ovar, only : cfl,cfl_beg,cfl_end,cfl_cycles,cfl_start_iter
  use ovar, only : adjust_cfl_cycles
  use ovar, only : output_time_averaging,time_scaling_factor
  use ovar, only : time_average_restart_files,ave_start_time
  use ovar, only : Filtering_Interval,Limiting_Interval
  use ovar, only : move_bcin_to_device
  use ovar, only : channel_spatial_ave_restart
  use ovar, only : if_interpto_newgrid
  use ovar, only : convert_restart_to_cgns
  use ovar, only : start_averaging_timestep
  use ovar, only : start_averaging_tstar,init_tstar
  use ovar, only : current_rk_stage
  use ovar, only : Filter_Option
  !
  use input_namelist_mod, only : cfl_adjustment
  use input_namelist_mod, only : restart_interupt
  !
  use flowvar, only : usp,dusp
  !
  use order_mod, only : n_order,p_order,q_order
  !
  use quadrature_mod, only : solpts_edge
  !
  use io_mod, only : read_input_file
  use io_mod, only : write_final_error
  use io_mod, only : write_parallel_cgns
  use io_mod, only : write_residual_file
  use io_mod, only : test_communication
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
  use vandermonde_mod, only : move_vand_to_device
  use vandermonde_mod, only : move_monomials_to_device
  !
  use derivatives_mod, only : init_derivative_matrices
  use derivatives_mod, only : move_deriv_to_device
  !
  use correction_mod, only : init_correction_matrices
  use correction_mod, only : move_correct_to_device
  !
  use interpolation_mod, only : init_interpolation_matrices
  use interpolation_mod, only : init_outerpolation_matrices
  use interpolation_mod, only : move_interp_to_device
  use interpolation_mod, only : move_outerp_to_device
  !
  use projection_mod, only : init_projection_matrices
  use projection_mod, only : move_projct_to_device
  use projection_mod, only : project_to_porder
  !
  use filter_mod, only : init_filter_matrices
  use filter_mod, only : move_filter_to_device
  use filter_mod, only : filter_usp
  !
  use initialization_mod, only : initialize_flow
  use initialization_mod, only : initialize_dusp
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
  use time_mod, only : update_Classic_RK
  use time_mod, only : update_CK_RK4
  use time_mod, only : update_TVD_RK2
  use time_mod, only : update_TVD_RK3
  use time_mod, only : update_TVD_RK4
  use time_mod, only : save_old_solution
  !
  use flux_mod, only : interp_cv_to_flxpts
  use flux_mod, only : interp_dusp_to_flxpts
!!use flux_mod, only : flux_gradient_at_solpts_visc
!!use flux_mod, only : flux_gradient_at_solpts_visc_gpu
!!use flux_mod, only : correct_interface_flux_visc
  use flux_mod, only : discontinuous_gradient_at_solpts
!!use flux_mod, only : correct_gradients_at_solpts
  use flux_mod, only : correct_gradients_at_solpts_gpu
!!use flux_mod, only : correct_gradients_br2
  use flux_mod, only : correct_gradients_br2_gpu
  use flux_mod, only : flux_gradient_at_solpts
  use flux_mod, only : get_invs_common_fluxes_interior
  use flux_mod, only : get_visc_common_fluxes_interior
  use flux_mod, only : get_invs_common_fluxes_bc
  use flux_mod, only : get_visc_common_fluxes_bc
  use flux_mod, only : init_cbias
  !
  use reactions_mod, only : compute_wdot
  !
  use parallel_mod, only : partition_grid
  use parallel_mod, only : create_serial_cell_map
  use parallel_mod, only : exchange_faceusp
  use parallel_mod, only : exchange_facedusp
  use parallel_mod, only : wait_faceusp
  use parallel_mod, only : wait_facedusp
  !
  use module_limiters, only : initialize_limiters
  !
  use mms_mod, only : mms_src
  use mms_mod, only : initialize_mms
  !
  use channel_mod, only : channel_source
  use channel_mod, only : spatial_ave_channel_flow
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
  use averaging_mod, only : initialize_time_ave
  use averaging_mod, only : update_time_ave
  use averaging_mod, only : average_restart_files
  !
  use synthetic_eddy_mod, only : output_sem_averages
  !
  use interpsol_mod, only : mv_global_to_grid2
  use interpsol_mod, only : mv_grid2_to_global
  use interpsol_mod, only : interpsol
!-DEBUG_BEG:
  use interpsol_mod, only : interpsol_debug
  use mydebug_mod,   only : mydebug
!-DEBUG_END:
  !
  use, intrinsic :: iso_c_binding,   only : c_long
  !
  implicit none
  !
  character(len=*), parameter :: pname = "GFR_main"
  !
  integer  :: iter,nst,i,ierr,k,npart
  integer  :: completed_timesteps
  integer  :: newton_failures
  real(wp) :: dtmin,out_time
  real(wp) :: current_tstar
  !
  character(len=50) :: memfile
  !
  integer :: CGNS_action
  integer :: input_mode
  !
  logical(lk)  :: completed_initial_cgns_output
  logical(lk)  :: exit_main_loop
  logical(lk)  :: walltime_expiring
  logical(lk)  :: output_solution
  logical(lk)  :: output_restart
  !
  logical(ldk) :: request_stop
  logical(ldk) :: request_restart
  logical(ldk) :: request_solution
  logical(ldk) :: request_cfl_adjustment
  logical(ldk) :: request_restart_interupt
  logical(ldk) :: request_time_ave
  logical(ldk) :: request_reactions
  logical(ldk) :: request_array(1:7)
  !
  logical(lk)  :: Apply_Filter
  logical(lk)  :: Apply_Limiter
  !
 !logical(lk), parameter :: ignore_errors = true
  logical(lk), parameter :: ignore_errors = fals
  !
  real(wp) :: cfl_tmp(1:5)
  !
  integer :: final_timestep_from_stop_request
  !
continue
  !
  call check_for_bad_number(__LINE__)
  request_restart = .false.
  request_solution = .false.
  request_time_ave = .false.
  request_restart_interupt = .false.
  request_cfl_adjustment = .false.
  request_reactions = .false.
  request_stop = .false.
  completed_initial_cgns_output = fals
  final_timestep_from_stop_request = 0
  !
  ! Initialize the overall MPI session
  !
  call initialize_mpi_environment
  call check_for_bad_number(__LINE__)
  !
  call get_cpu_locality
  call check_for_bad_number(__LINE__)
  !
  ! Initialize the wall time
  !
  call gfr_cpu_time(0)
  call check_for_bad_number(__LINE__)
  !
  ! Get the number of seconds allowed for this simulation
  !
  call pbs_remaining_walltime(walltime_expiring)
  call check_for_bad_number(__LINE__)
  !
  ! Initialize all the random constants that must be initialized at run time
  !
  call initialize_runtime_parameters
  call memory_pause("base initializations","reading input and grid files")
  call check_for_bad_number(__LINE__)
  !
  ! Perform preprocessing loop twice, if necessary.
  !   input_mode = 1  Interpolation Mode
  !   input_mode = 2  Normal Solver Mode
  ! During mode 1, try to load data into the interpolation destination set.
  ! However, if interpto_gridfile name is empty then drop into mode 2.
  ! During mode 2, load data into the normal solver set.
  preprocessing_loop: &
  do input_mode = 1,2
    !
    ! Read in the input file
    !
    call read_input_file(input_mode)
    call memory_pause("reading input and grid files","partitioning the grid")
    call check_for_bad_number(__LINE__)
    !
    ! Partition the grid if using more than one processor.
    ! Upon entry:  These arrays define the global grid
    ! Upon return: These arrays define the partition local
    !              to the respective processor along with
    !              connectivities for the partition boundaries
    npart = ncpu
    if (input_mode == 2) npart = 1
    if (npart > 1) then
      call partition_grid(npart,metis_option,bface,nodes_of_cell, &
                          nodes_of_cell_ptr,xyz_nodes,cell_geom,cell_order)
      call memory_pause("partitioning the grid","getting the solution points")
    else
      call create_serial_cell_map(ncell)
      call memory_pause("getting serial cell map","getting the solution points")
    end if
    call check_for_bad_number(__LINE__)
    !
    ! Get the coordinates of the solution points within each grid cell
    !
!-DEBUG_BEG:
!!!!call mydebug(pname,__LINE__,__FILE__)
!!!!write(iout,*) 'DEBUG[Euler]:  Calling get_solution_points next.'
!!!!flush(iout)
!-DEBUG_END:
    call get_solution_points
!-DEBUG_BEG:
!!!!write(iout,*) 'DEBUG[Euler]:  Back from get_solution_points.'
!!!!flush(iout)
!!!!call mydebug(pname,__LINE__,__FILE__)
!-DEBUG_END:
    call memory_pause("getting the solution points","getting grid connectivity")
    call check_for_bad_number(__LINE__)
    !
    ! Get the grid connectivity and metrics
    !
    call get_grid_connectivity(bface,cell_geom,cell_order,xyz_nodes, &
                               nodes_of_cell,nodes_of_cell_ptr)
    !
    call memory_pause("getting grid connectivity", &
                      "initializing vandermonde matrices")
    call check_for_bad_number(__LINE__)
    !
    ! Compute the remaining quadrature information
    ! (Vandermonde matrix, stiffness matrix, LIFT / correction function, etc.)
    !
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    call init_vandermonde_matrices
    call memory_pause("initializing vandermonde matrices", &
                      "initializing derivative matrices")
    call check_for_bad_number(__LINE__)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    call init_derivative_matrices
    call memory_pause("initializing derivative matrices", &
                      "initializing correction matrices")
    call check_for_bad_number(__LINE__)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    call init_correction_matrices
    call memory_pause("initializing correction matrices", &
                      "initializing interpolation matrices")
    call check_for_bad_number(__LINE__)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    call init_interpolation_matrices
    call memory_pause("initializing interpolation matrices", &
                      "initializing outerpolation matrices")
    call check_for_bad_number(__LINE__)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    call init_outerpolation_matrices
    call memory_pause("initializing outerpolation matrices", &
                      "initializing projection matrices")
    call check_for_bad_number(__LINE__)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    call init_projection_matrices
    call memory_pause("initializing projection matrices", &
                      "initializing filter matrices")
    call check_for_bad_number(__LINE__)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    call init_filter_matrices
    call memory_pause("initializing filter matrices", &
                      "initializing boundary profile stuff")
    call check_for_bad_number(__LINE__)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    call init_bnd_profile_structures
    call memory_pause("initializing boundary profile stuff", &
                      "initializing metrics")
    call check_for_bad_number(__LINE__)
    !
    ! Get the grid metrics
    !
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    call get_metrics
    !
    ! Special handling for interpolation mode (input_mode=1)
    !
    if (input_mode == 1) then
      if ( if_interpto_newgrid() ) then
        ! We will be interpolating the sol from one grid to another.
        ! Data for the destination grid to be used with interpolation
        ! has just been loaded into the normal global arrays.
        !
        ! Allocate the solution array, which would normally be done
        ! during initialization.
        allocate ( usp(1:n_solpts,1:nq) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"usp",1,__LINE__,__FILE__,ierr,error_message)
        !
        ! Move this data to another set of variables, so that during the
        ! second pass of this preprocessing loop we can read the
        ! interpolation source grid/solution into these same global variables.
        call mv_global_to_grid2()
!-DEBUG_BEG:
!-DEBUG_CMT: Print debug info to check if arrays moved ok.
        if (mypnum == glb_root) then
          write(iout,*) 'Data has been moved to the g2_* variables.'
        end if
        call mydebug(pname,__LINE__,__FILE__)
!-DEBUG_END:
        call memory_pause("initializing metrics", &
                          "rereading input and grid files")
        !
      else
        ! We will not be interpolating the sol from one grid to another.
        ! Data was loaded into the normal solver set instead.
        ! No need to proceed with second input_mode pass.
        exit preprocessing_loop
      end if
    end if
    !
  end do preprocessing_loop
  !
  call memory_pause("initializing metrics","initializing solution")
  call check_for_bad_number(__LINE__)
  !
  ! Uncomment the call to check_flux_point_orderings if you
  ! suspect that there is an error between the orderings of
  ! the flux points from the two cells on an interface
  !
  if (check_flux_point_coordinates) then
    call test_communication
    call check_face_flux_point_orderings
    call check_cell_flux_point_orderings
  end if
  !
  ! Output information about the global grid
  !
  if (mypnum == glb_root) then
    write (iout,11) n_global_cell, n_global_node, n_global_fbnd, &
                    n_global_solpts, &
                    n_solpts,n_flxpts,n_bnd_flxpts,n_wall_flxpts, &
                    n_qdtpts,n_flx_qdtpts,n_bnd_qdtpts,n_wall_qdtpts, &
                    total_volume, &
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
  if (output_time_averaging) then
    call initialize_time_ave
  end if
  call check_for_bad_number(__LINE__)
  call memory_pause("initializing solution","initializing generic module")
  call init_cbias
  !
  ! Interpolate solution onto a new grid.
  !
  if ( if_interpto_newgrid() ) then
    !
    ! Start by projecting the source solution to p_order
    !
    if (p_order < n_order) then
      call project_to_porder(usp)
    end if
    !
!-DEBUG_BEG:
!-DEBUG_CMT:  Write source grid/sol to cgns to examine, then stop
!-DEBUG_CMT:  before interpolating.
    if (.false.) then
   !if (.true.) then
      if (mypnum == glb_root) then
        write (iout,*) 'DEBUG[Euler] Write source grid/solution to cgns.'
        flush (iout)
      end if
      call write_parallel_cgns(open_CGNS)
      write (error_message,210)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
!-DEBUG_END:
    !
    ! Interpolate current global solution to grid 2.
    !
    if (mypnum == glb_root) write (iout,201)
    !
!-DEBUG_BEG:
!-DEBUG_CMT: Debug interpolation at a given point and stop.
!-DEBUG_CMT: Much faster than interpolating a whole grid, and
!-DEBUG_CMT: could test a location that is not a grid point.
    if (.false.) then
      call gfr_cpu_time(0)
      call interpsol_debug()
      call gfr_cpu_time(1)
      call gfr_cpu_time(-1)
      write (error_message,211)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
!-DEBUG_END:
    !
    call gfr_cpu_time(0)
    call interpsol()
    if (mypnum == glb_root) then
      write (iout,202)
      flush (iout)
    end if
    call gfr_cpu_time(1)
    call gfr_cpu_time(-1)
    !
    ! Move grid 2 data to current global arrays.
    ! This puts the destination grid/solution back into the global arrays.
    !
    call mv_grid2_to_global()
    !
    ! Finish by projecting or filtering the destination solution
    !
    if (p_order < n_order) then
      call project_to_porder(usp)
    else
      if (Filter_Option > 0) then
        call filter_usp
      end if
    end if
    !
    ! Allocate the derivative array, which would normally be done
    ! during initialization.  This is only needed for meaningful
    ! derivative quantities in CGNS output, else they will be zero.
    !
    if (.not. allocated(dusp)) then
      allocate ( dusp(1:n_solpts,1:nr,1:nq) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"dusp",1,__LINE__,__FILE__,ierr,error_message)
     !call initialize_dusp
    end if
    !
    ! Write restart file of the interpolated solution.
    !
    if ( mypnum == glb_root ) write (iout,207)
    call write_restart_file(itcur)
    !
    ! Write CGNS file of the interpolated solution.
    !
    if ( mypnum == glb_root ) write (iout,208)
    call write_parallel_cgns(open_CGNS)
    !
    ! Force termination with a message.
    !
    write (error_message,210)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    !
  end if
  !
  ! Time average the list of restart files
  !
  if (time_average_restart_files) then
    call average_restart_files
  end if
  !
  call move_face_to_device
  call move_cell_to_device
  call move_fp_to_device
  call move_bcgrps_to_device
  call move_correct_to_device
  call move_deriv_to_device
  call move_interp_to_device
  call move_outerp_to_device
  call move_filter_to_device
  call move_projct_to_device
  call move_vand_to_device
  call move_monomials_to_device
  call move_cellmap_to_device
  call move_bcin_to_device
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
  call check_for_bad_number(__LINE__)
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%                                    %%%
  !%%%  Beginning of main iteration loop  %%%
  !%%%                                    %%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  call initialize_limiters
  call check_for_bad_number(__LINE__)
  !
  ! Initialize the wall time
  !
  call gfr_cpu_time(0)
  !
  call reset_fp_exceptions
  !
 !if (num_timesteps > 0) then
    !
    call interp_cv_to_flxpts(0,0) ! ; call check_for_bad_number(__LINE__)
    call memory_pause("initial interpolation of C.V. to flux points", &
                      "exchanging C.V. between communication faces")
    call check_for_bad_number(__LINE__)
    if (ncpu > 1) then
      call exchange_faceusp       ! ; call check_for_bad_number(__LINE__)
!#define IMMEDIATE_WAIT_FOR_MPI
#ifdef IMMEDIATE_WAIT_FOR_MPI
      call wait_faceusp           ! ; call check_for_bad_number(__LINE__)
#endif
    end if
    call memory_pause("exchanging C.V. between communication faces", &
                      "initial evaluation of boundary conditions")
    call check_for_bad_number(__LINE__)
    call get_faceusp_bc(newton_failures) ! ; call check_for_bad_number(__LINE__)
    call memory_pause("initial evaluation of boundary conditions", &
                      "waiting for C.V. exchanges to complete")
    call check_for_bad_number(__LINE__,newton_failures)
    if (ncpu > 1) then
#ifndef IMMEDIATE_WAIT_FOR_MPI
      call wait_faceusp           ! ; call check_for_bad_number(__LINE__)
#endif
    end if
    call memory_pause("waiting for C.V. exchanges to complete", &
                      "initial TGV KEDR evaluation")
    call check_for_bad_number(__LINE__)
    !
    ! Get initial time step values in case they are needed for CGNS output
    !
    call get_dtsp(dtmin,final_timestep_from_stop_request)
    call check_for_bad_number(__LINE__)
    !
    if (itestcase == Taylor_Green_Vortex) then
      call compute_TaylorGreen_KEDR(0,zero)
      if (dump_memory_usage) then
        call memory_pause("initial TGV KEDR evaluation","dumping memory usage")
      else
        call memory_pause("initial TGV KEDR evaluation", &
                          "writing initial CGNS file")
      end if
      call check_for_bad_number(__LINE__)
    end if
    !
 !end if
  !
  ! Spatial average a channel flow restart file
  !
  if (channel_spatial_ave_restart) then
    call spatial_ave_channel_flow
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
    call new_write_parallel_cgns(open_CGNS)
  else
    call write_parallel_cgns(open_CGNS)
  end if
  if (convert_restart_to_cgns) then
    if (.not. continuous_output) then
      call write_parallel_cgns(close_CGNS)
    end if
    call write_restart_file(-1)
    write (error_message,209)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  completed_initial_cgns_output = true
  CGNS_action = add_to_CGNS
  call memory_pause("writing initial CGNS file","main iteration loop")
  call check_for_bad_number(__LINE__)
  !
  ! Initialize uoldsp and gpoldsp if num_timesteps is zero
  !
  if (num_timesteps == 0) then
    call save_old_solution
    call check_for_bad_number(__LINE__)
  end if
  !
  ! Initialize the wall time
  !
  call gfr_cpu_time(0)
  !
  restart_loop: &
  do
    !
    main_loop: &
    do iter = 1, num_timesteps
      !
      ! Current iteration count including restart iterations
      !
      itcur = itcur + 1
      !$ACC UPDATE DEVICE(itcur)
      !
      ! Store a copy of the solution from the previous time step
      !
      call save_old_solution              !; call check_for_bad_number(__LINE__)
      !
      ! Loop over the Runge-Kutta stages
      !
      rk_loop: &
      do nst = 1,num_rk_stages
        !
        current_rk_stage = nst
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
        !
        call interp_cv_to_flxpts(iter,nst)  
                                           ! call check_for_bad_number(__LINE__)
        !
        if (ncpu > 1) then
          call exchange_faceusp           !; call check_for_bad_number(__LINE__)
#ifdef IMMEDIATE_WAIT_FOR_MPI
          call wait_faceusp               !; call check_for_bad_number(__LINE__)
#endif
        end if
        !
        if (governing_equations == NavierStokes_Eqns) then
          call discontinuous_gradient_at_solpts
                                          !  call check_for_bad_number(__LINE__)
        end if
        !
        call get_faceusp_bc(newton_failures)
        call check_for_bad_number(__LINE__,newton_failures)
        !
        call get_invs_common_fluxes_interior
                                          !  call check_for_bad_number(__LINE__)
        !
        if (ncpu > 1) then
#ifndef IMMEDIATE_WAIT_FOR_MPI
          call wait_faceusp               !; call check_for_bad_number(__LINE__)
#endif
        end if
        !
        if (governing_equations == NavierStokes_Eqns) then
          if (visc_flux_method == visc_flux_BR2) then
           !call correct_gradients_br2
            call correct_gradients_br2_gpu  
                                           ! call check_for_bad_number(__LINE__)
          else
            call correct_gradients_at_solpts_gpu
                                          !  call check_for_bad_number(__LINE__)
           !call correct_gradients_at_solpts 
                                           ! call check_for_bad_number(__LINE__)
            call interp_dusp_to_flxpts    !; call check_for_bad_number(__LINE__)
          end if
          if (ncpu > 1) then
            call exchange_facedusp        !; call check_for_bad_number(__LINE__)
#ifdef IMMEDIATE_WAIT_FOR_MPI
            call wait_facedusp            !; call check_for_bad_number(__LINE__)
#endif
          end if
          call get_facedusp_bc            !; call check_for_bad_number(__LINE__)
          call get_visc_common_fluxes_interior
                                          !  call check_for_bad_number(__LINE__)
        end if
        !
        call get_invs_common_fluxes_bc    !; call check_for_bad_number(__LINE__)
        !
        if (governing_equations == NavierStokes_Eqns) then
          if (ncpu > 1) then
#ifndef IMMEDIATE_WAIT_FOR_MPI
            call wait_facedusp            !; call check_for_bad_number(__LINE__)
#endif
          end if
          call get_visc_common_fluxes_bc  !; call check_for_bad_number(__LINE__)
        end if
        !
        call flux_gradient_at_solpts(iter,nst)
                                          !  call check_for_bad_number(__LINE__)
        !
        ! Compute the source terms for method of manufactured solutions (MMS)
        !
        if (mms_opt /= 0) then
          call mms_src                    !; call check_for_bad_number(__LINE__)
        end if
        !
        ! Compute the source terms for the channel flow
        !
        if (itestcase == Channel_Flow) then
          call channel_source             !; call check_for_bad_number(__LINE__)
        end if
        !
        ! Compute the species source terms for chemically reacting flows
        !
        call compute_wdot(iter)           !; call check_for_bad_number(__LINE__)
        !
        ! Compute the final residual for each solution point
        !
        call rresi(iter,nst)              !; call check_for_bad_number(__LINE__)
        !
        if (nst == 1) then
          call get_dtsp(dtmin,final_timestep_from_stop_request)
                                          !  call check_for_bad_number(__LINE__)
        end if
        !
        ! Store residual in rssprk(:,:,nst) and update the solutions in usp
        !
        if (Runge_Kutta_Scheme == Classic_RK) then
          call update_Classic_RK(nst,Apply_Filter,Apply_Limiter)
        else if (Runge_Kutta_Scheme == TVD2_RK) then
          call update_TVD_RK2(nst,Apply_Filter,Apply_Limiter)
        else if (Runge_Kutta_Scheme == TVD3_RK) then
          call update_TVD_RK3(nst,Apply_Filter,Apply_Limiter,iter)
        else if (Runge_Kutta_Scheme == TVD4_RK) then
          call update_TVD_RK4(nst,Apply_Filter,Apply_Limiter)
        else if (Runge_Kutta_Scheme == CK4_RK) then
          call update_CK_RK4(nst,Apply_Filter,Apply_Limiter)
        end if
       !call check_for_bad_number(__LINE__)
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
          call gfr_cpu_time(1)
          !
          if (mypnum == glb_root) call write_residual_file(0,itcur)
          !
        end if
        !
      end if
      !
      call check_for_bad_number(__LINE__)
      !
      ! Initialize time averaging if start_averaging_timestep
      ! is equal to the current time step
      !
      if (.not. output_time_averaging) then
        if (start_averaging_timestep > 0) then
          if (itcur == start_averaging_timestep) then
            !
            if (mypnum == glb_root) then
              write (iout,80)
              flush (iout)
            end if
            !
            output_time_averaging = true
            rl2mx = -huge(zero)
            !
            call initialize_time_ave
            !
          end if
        end if
      end if
      !
      ! Initialize time averaging if the current non-dimensional elapsed
      ! simulation time (t* or t_star) since the beginning of the simulation
      ! is greater than start_averaging_tstar
      !
      if (.not. output_time_averaging) then
        if (start_averaging_tstar > zero) then
          !
         !current_tstar = time*time_ref*time_scaling_factor - init_tstar
          current_tstar = time*time_ref*time_scaling_factor
          !
          if (current_tstar > start_averaging_tstar) then
            !
            if (mypnum == glb_root) then
              write (iout,81)
              flush (iout)
            end if
            !
            output_time_averaging = true
            rl2mx = -huge(zero)
            !
            call initialize_time_ave
            !
          end if
          !
        end if
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
        if (iter /= num_timesteps) then
          !
          ! Have the root processor check for various temporary files
          !
          if (mypnum == glb_root) then
            !
            ! Check for a stop request
            !
            inquire (file="STOP_GFR",exist=request_stop)
            !
            ! If stop was requested, delete the
            ! STOP_GFR file and exit the main loop
            !
            if (request_stop) then
              !
              write (iout,99)
              !
              ! We dont care about the exit status here since we are
              ! going to be stopping execution after this anyways
              !
              open (newunit=i,file="STOP_GFR",iostat=ierr)
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
              cfl_tmp = [cfl,cfl_beg,cfl_end, &
                         real(cfl_cycles,kind=wp), &
                         real(cfl_start_iter,kind=wp)]
              !
              ! If we are adjusting the CFL cycling, automatically
              ! write out a restart file of the solution before the
              ! new CFL cycling takes effect.
              !
              if (cfl_start_iter < itcur) then
                request_restart = .true.
              end if
              !
            end if
            !
            ! Check for a start time-averaging request
            !
            inquire (file="START_TIME_AVERAGING",exist=request_time_ave)
            !
            ! If stop was requested, delete the
            ! STOP_GFR file and exit the main loop
            !
            if (request_time_ave) then
              !
              if (output_time_averaging) then
                request_time_ave = .false.
              else
                write (iout,94)
              end if
              !
              ! We dont care about the exit status here
              !
              open (newunit=i,file="START_TIME_AVERAGING",iostat=ierr)
              close (i,status="delete",iostat=ierr)
              !
            end if
           !!
           !! Check for a start reactions request
           !!
           !inquire (file="START_REACTIONS",exist=request_reactions)
           !!
           !! If reactions were requested, read the reactions namelist
           !! to get the iteration to begin reactions
           !!
           !if (request_reactions) then
           !  !
           !  write (iout,96)
           !  !
           !  open (newunit=i,file="START_REACTIONS",iostat=ierr)
           !  read (i,nml=start_reactions)
           !  close (i,status="delete",iostat=ierr)
           !  !
           !  ! If we are starting chemical reactions, automatically
           !  ! write out a restart file of the solution before the
           !  ! reactions begin.
           !  !
           !  if (reactions_start_iter < itcur) then
           !    request_restart = .true.
           !  end if
           !  !
           !end if
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
            request_array(:) = [request_stop,request_restart,request_solution, &
                                request_cfl_adjustment,request_restart_interupt, &
                                request_time_ave,request_reactions]
            !
          end if
          !
          ! Have the root processor broadcast all the request logicals to
          ! the other processors so they can proceed as needed
          !
          call mpi_bcast(request_array,size(request_array,kind=int_mpi), &
                         MPI_LOGICAL,0_int_mpi,MPI_COMM_WORLD,mpierr)
          !
          request_stop = request_array(1)
          request_restart = request_array(2)
          request_solution = request_array(3)
          request_cfl_adjustment = request_array(4)
          request_restart_interupt = request_array(5)
          request_time_ave = request_array(6)
          request_reactions = request_array(7)
          !
          ! Exit if the STOP_GFR file was found
          !
          if (request_stop) exit main_loop
          !
          ! If a restart dump was requested, write a restart file
          !
          if (request_restart) then
            call write_restart_file(itcur)
          end if
          !
          ! If a solution dump was requested, dump the solution and then let
          ! the CGNS subroutine know that the file needs to be reopened.
          !
          if (request_solution) then
            if (continuous_output) then
              call new_write_parallel_cgns(dump_CGNS)
            else
              call write_parallel_cgns(dump_CGNS)
            end if
            CGNS_action = reopen_CGNS
          end if
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
            cfl_start_iter = nint(cfl_tmp(5))
            !
            if (cfl_start_iter < itcur) then
              call adjust_cfl_cycles
            end if
            !
          end if
          !
          if (request_time_ave) then
            !
            output_time_averaging = true
            rl2mx = -huge(zero)
            !
            call initialize_time_ave
            !
          end if
          !
        end if
        !
      end if
      !
      if (iter /= num_timesteps) then
        !
        ! If the current time step is equal to cfl_start_iter,
        ! begin the adjustment of the CFL cycling
        !
        output_restart = fals
        if (cfl_start_iter == itcur) then
          call adjust_cfl_cycles
          output_restart = true
        end if
        !
        ! Output the restart file at specified step interval
        !
        if (restart_interval > 0 .and. .not.request_restart) then
          if (mod(itcur,restart_interval) == 0) then
            output_restart = true
          end if
        end if
        !
        if (output_restart .and. .not.request_restart) then
          call write_restart_file(itcur)
        end if
        !
        ! Output the solution at the specified interval
        !
        output_solution = fals
        if (output_interval > 0 .and. .not.request_solution) then
          if (mod(itcur,output_interval) == 0) then
            output_solution = true
          end if
        end if
        !
        if (is_at_output_frequency) then
          output_solution = true
          is_at_output_frequency = fals
        end if
        !
        if (output_solution .and. .not.request_solution) then
          if (continuous_output) then
            call new_write_parallel_cgns(CGNS_action)
          else
            call write_parallel_cgns(CGNS_action)
          end if
        end if
        !
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
            call new_write_parallel_cgns(CGNS_action)
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
    exit restart_loop
    !
  end do restart_loop
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
  call gfr_cpu_time(1)
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
    call new_write_parallel_cgns(close_CGNS)
  else
    call write_parallel_cgns(close_CGNS)
  end if
  !
  call output_sem_averages()
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
      write (iout,40) ncpu,itcur,out_time
      call gfr_cpu_time(-1)
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
  11 format (/,'  For the Global Grid ',/, &
               ' ---------------------',/, &
               '  ncell = ',i0,'  nnode = ',i0, &
               '  nfbnd = ',i0,'  nDOF = ',i0,//, &
               '  n_solpts      = ',i0,/, &
               '  n_flxpts      = ',i0,/, &
               '  n_bnd_flxpts  = ',i0,/, &
               '  n_wall_flxpts = ',i0,//, &
               '  n_qdtpts      = ',i0,/, &
               '  n_flx_qdtpts  = ',i0,/, &
               '  n_bnd_qdtpts  = ',i0,/, &
               '  n_wall_qdtpts = ',i0,//, &
               '  Total Grid Volume = ',f13.6,/, &
               '  1D Solution Points = ',100(6f9.5,:,/,23x))
  !
  19 format (/,' Lowest  Cell Order = ',i4,/, &
               ' Highest Cell Order = ',i4,/)
  !
  40 format (//,' Number of CPU cores used = ',i12,//, &
                ' Final iteration count    = ',i12,/, &
                ' Simulation elapsed time  = ',f12.3)
  41 format (/,"TOTAL SIMULATION WALL TIME WAS ",f18.6," s.",/,&
               "TOTAL COMPUTATIONAL COST WAS   ",f18.6," work units.")
  !
  50 format (/'Done'/)
  !
  80 format (/, &
    1x,72("*"),//, &
    25x,"START OF TIME AVERAGING!",//, &
    4x,"THE CURRENT ITERATION IS EQUAL TO THE VALUE OF THE INPUT PARAMETER",/, &
    4x,"START_AVERAGING_TIMESTEP THAT WAS GIVEN IN THE INPUT FILE!",//, &
    1x,72("*"),/)
  81 format (/, &
       1x,69("*"),//, &
       23x,"START OF TIME AVERAGING!",//, &
       4x,"THE CURRENT ELAPSED SIMULATION TIME SINCE THE BEGINNING OF THIS",/, &
       4x,"SIMULATION IS GREATER THAN THE VALUE OF THE INPUT PARAMETER",/, &
       4x,"START_AVERAGING_TSTAR THAT WAS GIVEN IN THE INPUT FILE!",//, &
       1x,69("*"),/)
  !
  94 format (/,"*******************************************************",/, &
               "    FOUND FILE REQUESTING TO START TIME AVERAGING! ",/, &
               "*******************************************************",/)
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
               " GFR IS NOW STOPPING EXECUTION AND", &
               " DUMPING THE FINAL RESULTS!",/, &
               "*******************************************************",/)
  99 format (/,"*******************************************************",/, &
               " FOUND FILE REQUESTING TO STOP THE EXECUTION OF GFR! ",/, &
               "   STOPPING EXECUTION AND DUMPING THE FINAL RESULTS!",/, &
               "*******************************************************",/)
  !
  201 format (/,"INTERPOLATING SOLUTION TO NEW GRID.")
  202 format (/,"INTERPOLATION COMPLETE.")
  207 format (/,"WRITING RESTART FILE.")
  208 format (/,"WRITING CGNS FILE.")
  209 format ("CONVERTING RESTART FILE TO CGNS SOLUTION FILE.")
  210 format ("STOPPING EXECUTION SO YOU MAY EXAMINE INTERPOLATED SOLUTION.")
  211 format ("STOPPING EXECUTION SO YOU MAY EXAMINE interpsol_debug.tec")
  !
contains
!
!###############################################################################
!
subroutine check_for_bad_number(line_num,newton_failures,force_signaling)
  !
  !.. Use Statements ..
  use eqn_idx, only : nq
  use geovar,  only : n_solpts
  use flowvar, only : usp,uoldsp
#ifndef SPECIAL_FOR_GCC
  use, intrinsic :: ieee_exceptions, only : ieee_usual
  use, intrinsic :: ieee_exceptions, only : ieee_get_flag
#endif
  !
  !.. Formal Arguments ..
  integer, intent(in) :: line_num
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: newton_failures
  integer, optional, intent(in) :: force_signaling
  !
#ifndef SPECIAL_FOR_GCC
  !
  !.. Local Scalars ..
  integer  :: n,m
  real(wp) :: out_time
  !
  logical(lk) :: error_found
  !
  character(len=300) :: overflow_error
  character(len=300) :: divide_by_zero_error
  character(len=300) :: invalid_error
  character(len=300) :: newton_error
  !
  character(len=10000) :: output_message
  !
  !.. Local Arrays ..
  integer      :: bad_number_found(1:size(ieee_usual))
  logical(ldk) :: signaling_exceptions(1:size(ieee_usual))
  !
  character(len=*), parameter :: pname = "check_for_bad_number"
  !
continue
  !
  error_found = fals
  !
  overflow_error = ""
  divide_by_zero_error = ""
  invalid_error = ""
  newton_error = ""
  output_message = ""
  !
  ! First check to make sure that no floating-point exceptions are signaling
  !
  call ieee_get_flag(ieee_usual,signaling_exceptions)
  !
  bad_number_found = merge(1,0,signaling_exceptions)
  !
  if (present(force_signaling)) then
    if (force_signaling > 0 .and. force_signaling <= size(ieee_usual)) then
      bad_number_found(force_signaling) = 1
    else if (force_signaling < 0) then
      bad_number_found(:) = 1
    end if
  end if
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,bad_number_found, &
                       size(bad_number_found,kind=int_mpi), &
                       mpi_inttyp,MPI_MAX,MPI_COMM_WORLD,mpierr)
  end if
  !
  if (present(newton_failures)) then
    if (newton_failures > 0) then
      error_found = true
      write (newton_error,104) newton_failures
    end if
  end if
  !
  if (sum(bad_number_found) > 0) then
    error_found = true
    if (bad_number_found(1) > 0) write (overflow_error,101)
    if (bad_number_found(2) > 0) write (divide_by_zero_error,102)
    if (bad_number_found(3) > 0) write (invalid_error,103)
  end if
  !
  if (error_found) then
    !
    ! A floating-point exception was detected or the subsonic inflow
    ! Newton iteration failed. Dump solution and restart files for the
    ! solution at the previous time step and abort the simulation.
    !
    ! Copy the solution from the previous time step back into
    ! usp so we can dump the solution and restart files before
    ! aborting the simulation
    !
    if (mypnum == glb_root) then
      !
      write (error_message,200) itcur
      output_message = reformat_string(error_message, &
                                       line_length=70, &
                                       leading_blanks=5)
      !
      write (iout,1) trim(output_message)
      flush (iout)
      !
      write (error_message,201) trim(adjustl(overflow_error)), &
                                trim(adjustl(divide_by_zero_error)), &
                                trim(adjustl(invalid_error)), &
                                trim(adjustl(newton_error))
      output_message = reformat_string(error_message, &
                                       line_length=53, &
                                       leading_blanks=21)
      !
      write (iout,2) trim(adjustl(output_message))
      flush (iout)
      !
    end if
    !
    if (completed_initial_cgns_output) then
      if (allocated(usp) .and. allocated(uoldsp)) then
        !
        !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(uoldsp,usp)
        do m = 1,nq
          do n = 1,n_solpts
            usp(n,m) = uoldsp(n,m)
          end do
        end do
        !
        if (continuous_output) then
          call new_write_parallel_cgns(-3)
        else
          call write_parallel_cgns(-3)
        end if
        !
        call write_restart_file(-2)
        !
      end if
    end if
    !
    out_time = time*time_ref*time_scaling_factor
    write (error_message,100) itcur,out_time, &
                              trim(adjustl(overflow_error)), &
                              trim(adjustl(divide_by_zero_error)), &
                              trim(adjustl(invalid_error)), &
                              trim(adjustl(newton_error))
    !
    call stop_gfr(stop_mpi,pname,line_num,__FILE__,error_message)
    !
  end if
  !
  1 format (//,1x,79("#"),//,5x,a,/)
  2 format (5x,"ERROR MESSAGE: '",a,"'",//,1x,79("#"),//)
  !
  100 format ("A SIGNALING FLOATING POINT EXCEPTION WAS ", &
              "DETECTED AFTER UPDATING THE SOLUTION AFTER ", &
              "TIME STEP # ",i0," AT SOLUTION TIME ",es18.10, &
              "!!! DUMPING SOLUTION AND ", &
              "RESTART FILES FOR THE PREVIOUS TIME STEP!!!", &
              4(2x,a,:))
  !
  101 format ("OVERFLOW EXCEPTION DETECTED!!!")
  102 format ("DIVIDE-BY-ZERO EXCEPTION DETECTED!!!")
  103 format ("INVALID NUMBER (NaN) EXCEPTION DETECTED!!!")
  104 format (" FAILURES WERE DETECTED USING NEWTON'S METHOD", &
              " TO COMPUTE THE SUBSONIC INFLOW TEMPERATURE!", &
              "       newton_failures = ",i0)
  !
  200 format ("A SIGNALING FLOATING POINT EXCEPTION WAS ", &
              "DETECTED AFTER UPDATING THE SOLUTION AFTER ", &
              "TIME STEP # ",i0," !!! DUMPING SOLUTION AND ", &
              "RESTART FILES FOR THE PREVIOUS TIME STEP!!!")
  201 format (4(2x,a,:))
#endif
  !
end subroutine check_for_bad_number
!
!###############################################################################
!
subroutine reset_fp_exceptions
  !
#ifndef SPECIAL_FOR_GCC
  use, intrinsic :: ieee_exceptions, only : ieee_usual
  use, intrinsic :: ieee_exceptions, only : ieee_set_flag
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Arrays ..
  logical(ldk) :: reset_exceptions(1:size(ieee_usual))
  !
  character(len=*), parameter :: pname = "reset_fp_exceptions"
  !
continue
  !
  do n = 1,size(reset_exceptions)
    reset_exceptions(n) = .false.
  end do
  !
  ! First check to make sure that no floating-point exceptions are signaling
  !
  call ieee_set_flag(ieee_usual,reset_exceptions)
#endif
  !
end subroutine reset_fp_exceptions
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
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
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
end program GFR_main
!
!###############################################################################
!
