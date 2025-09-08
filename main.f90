program GFR_main
  !
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  !
  use geomvar_mod, only : ncell,n_global_cell
  use geomvar_mod, only : nnode,n_global_node
  use geomvar_mod, only : nfbnd,n_global_fbnd,n_global_face
  use geomvar_mod, only : nr,n_solpts
  use geomvar_mod, only : n_global_solpts
  use geomvar_mod, only : bface,xyz_nodes,nface
  use geomvar_mod, only : cell_geom,cell_order,global_cell_order
  use geomvar_mod, only : nodes_of_cell_ptr,nodes_of_cell
  use geomvar_mod, only : move_face_to_device
  use geomvar_mod, only : move_cell_to_device
  use geomvar_mod, only : move_fp_to_device
  use geomvar_mod, only : move_bcgrps_to_device
  use geomvar_mod, only : move_cellmap_to_device
  use geomvar_mod, only : cell
  use geomvar_mod, only : n_flxpts,n_bnd_flxpts,n_wall_flxpts
  use geomvar_mod, only : n_qdtpts,n_flx_qdtpts,n_bnd_qdtpts,n_wall_qdtpts
  use geomvar_mod, only : output_bface
  !
  use eqn_idx_mod, only : nq,ng
  !
  use ovar, only : metis_option,convergence_reached
  use ovar, only : total_volume,local_volume,itestcase
  use ovar, only : num_rk_stages,Runge_Kutta_Scheme
  use ovar, only : time,time_ref,itcur,num_timesteps,rl2mx
  use ovar, only : results_interval,iter_out_interval,output_eperk_cfl
  use ovar, only : restart_interval,output_interval
  use ovar, only : is_at_output_frequency
  use ovar, only : is_at_restart_frequency
  use ovar, only : dump_memory_usage,mms_opt
  use ovar, only : gfr_cpu_time,this_is_final_timestep
  use ovar, only : write_TaylorGreen_8s_solution
  use ovar, only : write_TaylorGreen_full_solution
  use ovar, only : output_KEDR
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
  use ovar, only : interpto_input_mode
  use ovar, only : spatial_ave_restart_z_planes
  use ovar, only : spatial_ave_restart_x_and_z_planes
  use ovar, only : CFL_reduction_for_NaN_restart
  use ovar, only : matrix_output_unit
  !
  use input_namelist_mod, only : cfl_adjustment
  use input_namelist_mod, only : restart_interupt
  !
  use flowvar_mod, only : usp,gpsp,dusp
  use flowvar_mod, only : faceusp,facegpsp
  !
  use order_mod, only : n_order,p_order,q_order
  !
  use element_type_mod, only : solpts_edge
  !
  use io_mod, only : read_input_file
  use io_mod, only : write_final_error
  use io_mod, only : write_residual_file
  use io_mod, only : write_summary_file
  use io_mod, only : test_communication
  use io_mod, only : check_face_flux_point_orderings
  use io_mod, only : check_cell_flux_point_orderings
  !
  use write_cgns_mod, only : write_parallel_cgns
  use write_cgns_mod, only : TaylorGreen_8s_vorticity_solution
  use write_cgns_mod, only : convert_restart_files_to_cgns
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
  use initialization_mod, only : get_face_and_edge_point_coordinates
  !
  use bc_mod, only : get_faceusp_bc
  !
  use generic_mod, only : init_generic_mod
  use generic_mod, only : write_iter_header
  use generic_mod, only : write_iter_stats
  use generic_mod, only : getmaxent
  use generic_mod, only : compute_residual_error
  use generic_mod, only : compute_KEDR
  !
  use dt_mod, only : get_dtsp
  use erk_mod, only : update_Classic_RK
  use erk_mod, only : update_CK_RK4
  use erk_mod, only : update_TVD_RK2
  use erk_mod, only : update_TVD_RK3
  use erk_mod, only : update_TVD_RK4
  use rhs_mod, only : save_old_solution
  !
  use oerk_mod, only : update_oerk
  !
  use eperk_mod, only : update_eperk
  use eperk_mod, only : get_eperk_adaptive_dtsp
  !
  use flux_mod, only : interp_cv_to_flxpts
  use flux_mod, only : init_cbias
  !
  use parallel_mod, only : partition_grid
  use parallel_mod, only : create_serial_cell_map
  use parallel_mod, only : exchange_faceusp
  use parallel_mod, only : wait_faceusp
  !
  use limiters_mod, only : initialize_limiters
  !
  use mms_mod, only : initialize_mms
  !
  use channel_mod, only : spatial_ave_channel_flow
  use channel_mod, only : spatial_ave_z_planes
  use channel_mod, only : spatial_ave_x_and_z_planes
  !
  use bnd_profiles_mod, only : init_bnd_profile_structures
  use bnd_profiles_mod, only : compute_bnd_profiles
  !
  use postprocessor_mod, only : new_write_parallel_cgns
  use postprocessor_mod, only : convert_restart_files_to_continuous_cgns
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
  use thx_mod, only : compute_blowing_ratio
  !
  use interpsol_mod, only : mv_global_to_grid2
  use interpsol_mod, only : mv_grid2_to_global
  use interpsol_mod, only : interpsol
!-DEBUG_BEG:
  use interpsol_mod, only : interpsol_debug
  use mydebug_mod,   only : mydebug
!-DEBUG_END:
  !
  use cgns_mod, only : root_broadcasts_grid_dt
  !
  use periodic_mod, only : check_periodic_list_indices
  !
  use thermo_mod, only : update_thermo_properties
  !
  use transport_mod, only : update_transport_properties
  !
  use timer_mod, only : solver_timer
  use timer_mod, only : check_solver_timer
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
  real(wp) :: next_dtmin
  integer  :: n_rejected_timesteps
  logical(lk) :: previous_timestep_rejected
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
  !
  integer :: request_array(1:7)
  !
  logical(lk)  :: Apply_Filter
  logical(lk)  :: Apply_Limiter
  !
 !logical(lk), parameter :: ignore_errors = true
  logical(lk), parameter :: ignore_errors = fals
  !
  real(wp) :: cfl_tmp(1:5)
  real(wp) :: cfl_save(1:5)
  !
  integer :: n_par_fbnd,nfint
  integer :: n_sum_fbnd,n_sum_fint
  !
  integer :: final_timestep_from_stop_request
  integer :: previous_stop_request_timestep
  !
  character(len=10)  :: text
  character(len=100) :: dbgfile,dbgfrmt
  !
  integer            :: num_CFL_restarts
  integer, parameter :: max_CFL_restarts = 0
  !
  real(wp) :: br,mdot_inj,mdot_noz
  !
  integer, save :: idx1 = -1
  !
continue
  !
  ! Save the default IEEE floating point status
  !
  call get_default_ieee_status
  !
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
  num_CFL_restarts = 0
  !
  n_rejected_timesteps = 0
  previous_timestep_rejected = fals
  !
  ! Initialize the overall MPI session
  !
  call initialize_mpi_environment
  call timed_check_for_bad_number(__LINE__)
  !
#ifdef DEBUG_VERBOSE_MPI_CALLS
  write (text,'(i0)') ncpu
  write (dbgfrmt,850) len_trim(adjustl(text)),len_trim(adjustl(text))
  850 format ('("mpi_calls.CPU-",i',i0,'.',i0,',".log")')
  write (dbgfile,dbgfrmt) mypnum
  open (newunit=my_dbg_out,file=dbgfile,status="replace",action="write", &
                           iostat=ierr,iomsg=error_message)
  call io_error(pname,dbgfile,1,__LINE__,__FILE__,ierr,error_message)
#endif
  !
  call get_cpu_locality
  call timed_check_for_bad_number(__LINE__)
  !
  ! Initialize the wall time
  !
  call gfr_cpu_time(0)
  call timed_check_for_bad_number(__LINE__)
  !
  ! Get the number of seconds allowed for this simulation
  !
  call pbs_remaining_walltime(walltime_expiring)
  call timed_check_for_bad_number(__LINE__)
  !
  ! Initialize all the random constants that must be initialized at run time
  !
  call initialize_runtime_parameters
  call memory_pause("base initializations","reading input and grid files")
  call timed_check_for_bad_number(__LINE__)
  !
  if (mypnum == glb_root) then
    if (matrix_output_unit /= iout) then
      open (newunit=matrix_output_unit,file="matrix_output.dat", &
            status="replace",action="write", &
            iostat=ierr,iomsg=error_message)
      call io_error(pname,"matrix_output.dat",1,__LINE__,__FILE__, &
                    ierr,error_message)
    end if
  end if
  if (ncpu > 1) then
    call mpi_bcast(matrix_output_unit,1_int_mpi,mpi_inttyp, &
                   glb_root,MPI_COMM_WORLD,mpierr)
  end if
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
    if ( if_interpto_newgrid() ) then
      interpto_input_mode = input_mode
    end if
    !
    call output_bface(strng="after.read_input.",only_glb_root=true)
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
    if (input_mode == 2) then
      call root_broadcasts_grid_dt
    end if
    !
    call output_bface(strng="after.partitioning.",only_glb_root=true)
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
    call output_bface(strng="after.get_sp.",only_glb_root=true)
    !
    ! Get the grid connectivity and metrics
    !
    call get_grid_connectivity(bface,cell_geom,cell_order,xyz_nodes, &
                               nodes_of_cell,nodes_of_cell_ptr,input_mode)
    !
    call output_bface(strng="after.grid_connectivity.",only_glb_root=true)
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
    call check_for_bad_number(__LINE__)
    !
   !error_message = "STOP_FOR MATRIX OUTPUT"
   !call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
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
        ! Allocate the gas properties array, which would normally be done
        ! during initialization.
        allocate ( gpsp(1:n_solpts,1:ng) , source=zero , &
                   stat=ierr , errmsg=error_message )
        call alloc_error(pname,"gpsp",1,__LINE__,__FILE__,ierr,error_message)
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
  call check_periodic_list_indices
  if (check_flux_point_coordinates) then
    call get_face_and_edge_point_coordinates(force_facexyz=true)
    call test_communication
    call check_face_flux_point_orderings
    call check_cell_flux_point_orderings
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  !
  ! Output information about the global grid
  !
  n_sum_fbnd = 0
  n_par_fbnd = 0
  n_sum_fint = 0
  nfint = nface-nfbnd
  call mpi_allreduce(nfbnd,n_sum_fbnd,1_int_mpi,mpi_inttyp, &
                     MPI_SUM,MPI_COMM_WORLD,mpierr)
  call mpi_allreduce(nfint,n_sum_fint,1_int_mpi,mpi_inttyp, &
                     MPI_SUM,MPI_COMM_WORLD,mpierr)
  n_par_fbnd = (n_sum_fbnd - n_global_fbnd)/2
  !
  if (mypnum == glb_root) then
    write (iout,11) n_global_cell, n_global_node, n_global_fbnd, &
                    n_global_face-n_global_fbnd, n_global_face, &
                    n_sum_fbnd,n_par_fbnd,n_sum_fint, &
                    n_global_solpts, &
                    n_solpts,n_flxpts,n_bnd_flxpts,n_wall_flxpts, &
                    n_qdtpts,n_flx_qdtpts,n_bnd_qdtpts,n_wall_qdtpts, &
                    total_volume,nr, &
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
      call project_to_porder(usp,update_gpsp=true)
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
      call project_to_porder(usp,update_gpsp=true)
    else
      if (Filter_Option > 0) then
        call filter_usp(update_gpsp=true)
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
    call update_thermo_properties(1,n_bnd_flxpts,faceusp,facegpsp, &
                                  caller_name=pname)
    call check_for_bad_number(__LINE__)
    call update_transport_properties(1,n_bnd_flxpts,faceusp,facegpsp)
    call check_for_bad_number(__LINE__)
    !
    ! Get initial time step values in case they are needed for CGNS output
    !
    call get_dtsp(dtmin,final_timestep_from_stop_request)
    call check_for_bad_number(__LINE__)
    !
    if (output_KEDR) then
      call compute_KEDR(0,zero)
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
  else
    if (spatial_ave_restart_z_planes .or. &
        spatial_ave_restart_x_and_z_planes) then
      if (spatial_ave_restart_z_planes) then
        call spatial_ave_z_planes
      end if
      if (spatial_ave_restart_x_and_z_planes) then
        call spatial_ave_x_and_z_planes
      end if
      write (error_message,'(a)') "After spatially averaging a restart file."
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
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
  if (abs(itestcase) == THX_Phase_3) then
    call compute_blowing_ratio(br,mdot_noz,mdot_inj)
  end if
  !
  ! If requested, dump the memory usage of the advanced data types from some of
  ! the modules
  !
  if (dump_memory_usage) then
    call report_memory_usage(n=3,request_stop=true)
  end if
  !
  ! Reset the IEEE floating point status to the default values
  !
  call reset_ieee_status
  !
  if (convert_restart_to_cgns) then
    !
    ! Convert the current restart file or list
    ! of restart files into CGNS solution files.
    !
    if (.not. continuous_output) then
      call convert_restart_files_to_cgns
    else
      call convert_restart_files_to_continuous_cgns
    end if
    !
  else
    !
    ! Output the initialized solution
    !
    if (num_timesteps > 0) then
      if (continuous_output) then
        call new_write_parallel_cgns(open_CGNS)
      else
        call write_parallel_cgns(open_CGNS)
      end if
    end if
    !
  end if
  !
  completed_initial_cgns_output = true
  CGNS_action = add_to_CGNS
  call memory_pause("writing initial CGNS file","main iteration loop")
 !call check_for_bad_number(__LINE__)
  !
  ! Reset the IEEE floating point status to the default values
  !
  call reset_ieee_status
  !
  ! Initialize uoldsp and gpoldsp if num_timesteps is zero
  !
  if (num_timesteps == 0) then
    call save_old_solution
    call check_for_bad_number(__LINE__)
  end if
  !
  ! Write out the initial solution to a restart file before starting the solver.
  ! This is necessary for using the capability to automatically restart with
  ! a reduced CFL when a NaN is detected during the solver.
  !
  if (itcur == 0) then
    call write_restart_file(itcur)
  end if
  !
  call write_iter_header(output_header=true)
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
      call solver_timer(entering_procedure,pname,idx1)
      !
      ! Current iteration count including restart iterations
      !
      itcur = itcur + 1
      !$ACC UPDATE DEVICE(itcur)
      !
      ! Store a copy of the solution from the previous time step
      !
      call save_old_solution
      call debugging_checks_for_bad_number(__LINE__)
      !
      ! Loop over the Runge-Kutta stages
      !
      rk_loop: &
      do nst = 1,num_rk_stages
        !
        current_rk_stage = nst
        !
        Apply_Filter = (Filtering_Interval <= 0)
        Apply_Limiter = (Limiting_Interval <= 0)
        if (nst == num_rk_stages) then
          if (Filtering_Interval > 0) then
            Apply_Filter = (mod(iter,Filtering_Interval) == 0)
          end if
          if (Limiting_Interval > 0) then
            Apply_Limiter = (mod(iter,Limiting_Interval) == 0)
          end if
        end if
        !
        ! Compute the residual/RHS for this Runge-Kutta stage
        ! NOTE: We want to skip this if
        !         1) we are using an embedded-pair Runge-Kutta method,
        !         2) this is the first Runge-Kutta stage, and
        !         3) this is not the first time step
        !       because this was already done after computing the final stage
        !       of the previous time step.
        !
        if (abs(Runge_Kutta_Scheme) /= EPERK) then
          call compute_rhs(iter,nst,newton_failures)
        else ! Runge_Kutta_Scheme == EPERK
          if (iter == 1 .or. nst > 1) then
            call compute_rhs(iter,nst,newton_failures)
          end if
        end if
        !
        ! Compute the CFL-based time step size if this is the first
        ! Runge-Kutta stage and not using an embedded-pair RK method
        ! NOTE: We also need to do this if this is the first time step and
        !       we are using an embedded-pair Runge-Kutta method.
        !
        if (nst == 1) then
          if (Runge_Kutta_Scheme /= EPERK) then
            call get_dtsp(dtmin,final_timestep_from_stop_request)
            call debugging_checks_for_bad_number(__LINE__)
          else ! Runge_Kutta_Scheme == EPERK
            if (iter == 1) then
              call get_dtsp(dtmin,final_timestep_from_stop_request)
              call debugging_checks_for_bad_number(__LINE__)
            else if (output_eperk_cfl) then
              call get_dtsp(dtmin,final_timestep_from_stop_request, &
                            skip_updates=true)
            end if
          end if
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
        else if (Runge_Kutta_Scheme == OERK) then
          call update_OERK(nst,Apply_Filter,Apply_Limiter)
        else if (abs(Runge_Kutta_Scheme) == EPERK) then
          call update_EPERK(nst,iter,Apply_Filter,Apply_Limiter)!, &
                           !dtmin,final_timestep_from_stop_request)
        end if
        call check_for_bad_number(__LINE__)
        !
      end do rk_loop
      !
      ! If using an embedded-pair Runge-Kutta method, compute the time error
      ! using the embedded-pair stage and use that to either reject the
      ! previous time step or compute a new dtmin for the next time step.
      !
      if (abs(Runge_Kutta_Scheme) == EPERK) then
        !
        current_rk_stage = 1 ! this needs to be 1 for SEM BCs
        call compute_rhs(iter,current_rk_stage,newton_failures)
        !
        call get_EPERK_adaptive_dtsp(iter,dtmin,next_dtmin, &
                                     final_timestep_from_stop_request, &
                                     previous_timestep_rejected)
        !
        ! If the previous time step was rejected, we need to do some stuff
        ! before redoing the previous time step with a new time step size.
        !
        if (previous_timestep_rejected) then
         !!
         !! With the previous time step rejected, the subroutine
         !! get_EPERK_adaptive_dtsp will reset usp back to uoldsp. We
         !! need to compute residual with the reset usp, because
         !! the Runge-Kutta loop wont compute this for the first stage.
         !! NOTE: This might technically make the SEM inconsistent because
         !!       calling this again will propagate the eddies twice in one
         !!       time step. However, I dont think this will cause any
         !!       significant issues.
         !!
         !call compute_rhs(iter,current_rk_stage,newton_failures)
          !
          ! Finally, adjust some of the counters and skip to the start of the
          ! main loop so redo the previous time step with a new time step size.
          !
          n_rejected_timesteps = n_rejected_timesteps + 1
          itcur = itcur - 1
          if (mypnum == glb_root) then
            write (iout,*) "TIME STEP REJECTED"
            flush (iout)
          end if
          call mpi_barrier(MPI_COMM_WORLD,mpierr)
          !
          dtmin = next_dtmin
          !
          cycle main_loop
          !
        else
          time = time + dtmin
          dtmin = next_dtmin
        end if
        !
      else
        !
        ! Update the current time now that the Runge-Kutta loop is finished
        !
        time = time + dtmin
        !
      end if
      !$ACC UPDATE DEVICE(time)
      !
      completed_timesteps = iter - n_rejected_timesteps
      !
      call debugging_checks_for_bad_number(__LINE__)
      !
      ! Update the time averaged variables
      !
      if (output_time_averaging) then
        call update_time_ave
      end if
      call debugging_checks_for_bad_number(__LINE__)
      !
      ! Compute and output the residual error at the specified interval
      !
      if ((mod(itcur,iter_out_interval) == 0 .or. this_is_final_timestep) .or. &
          (mod(itcur,results_interval) == 0 .and. iter/=num_timesteps)) then
        !
        ! Compute the residual error
        !
        call compute_residual_error
        !
        if (mod(itcur,iter_out_interval) == 0 .or. this_is_final_timestep) then
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
      call debugging_checks_for_bad_number(__LINE__)
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
      check_for_interupt_files: &
#ifdef DEBUG_ON
      if (mod(itcur,1) == 0) then
#else
#ifdef PBS_ENV
      if (mod(itcur,100) == 0) then
#else
      if (mod(itcur,50) == 0) then
#endif
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
            ! Save the previous time step from the last stop request
            ! NOTE: It probably isnt needed to save this to a new variable,
            !       but this just makes it clear that if request_stop is true
            !       the value read in from the STOP_GFR file isnt what is
            !       tested for the 'else if' part of the conditional.
            !
            previous_stop_request_timestep = final_timestep_from_stop_request
            !
            ! If stop was requested, delete the
            ! STOP_GFR file and exit the main loop
            !
            if (request_stop) then
              !
              ! We dont care about the exit status here since we are
              ! going to be stopping execution after this anyways
              !
              open (newunit=i,file="STOP_GFR",iostat=ierr)
              !
              read (i,*,iostat=ierr) final_timestep_from_stop_request
              !
              if (ierr /= 0) then
                final_timestep_from_stop_request = 1
              else if (final_timestep_from_stop_request <= 0) then
                final_timestep_from_stop_request = 1
              end if
              !
              close (i,status="delete",iostat=ierr)
              !
              if (final_timestep_from_stop_request <= itcur) then
                write (iout,99)
              else
                write (iout,199) final_timestep_from_stop_request
              end if
              !
            else if (previous_stop_request_timestep > 0) then
              !
              ! Manually set request_stop. If a previous stop request provided
              ! a final timestep and there was no STOP_GFR file when checking 
              ! this time, creating the request_array below will zero out the
              ! previous final timestep if request_stop is false. We want to
              ! keep the previous value of final_timestep_from_stop_request so
              ! the simulation stops as was requested.
              !
              request_stop = .true.
              final_timestep_from_stop_request = previous_stop_request_timestep
              !
            end if
            !
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
              ! Save the current CFL settings for the global root process
              !
              cfl_save = [cfl,cfl_beg,cfl_end, &
                          real(cfl_cycles,kind=wp), &
                          real(cfl_start_iter,kind=wp)]
              !
              open (newunit=i,file="ADJUST_CFL",iostat=ierr)
              read (i,nml=cfl_adjustment)
              close (i,status="delete",iostat=ierr)
              !
              ! Copy the CFL adjustment inputs into cfl_tmp so the global
              ! root processes can broadcast them to all the other processes.
              !
              cfl_tmp = [cfl,cfl_beg,cfl_end, &
                         real(cfl_cycles,kind=wp), &
                         real(cfl_start_iter,kind=wp)]
              !
              ! Restore the saved CFL values for the global root process
              ! now that the inputs from the CFL adjustment file have been
              ! copied into cfl_tmp
              !
              cfl = cfl_save(1)
              cfl_beg = cfl_save(2)
              cfl_end = cfl_save(3)
              cfl_cycles = nint(cfl_save(4))
              cfl_start_iter = nint(cfl_save(5))
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
            request_array(:) = [merge(1,0,request_stop), &
                                merge(1,0,request_restart), &
                                merge(1,0,request_solution), &
                                merge(1,0,request_cfl_adjustment), &
                                merge(1,0,request_restart_interupt), &
                                merge(1,0,request_time_ave), &
                                merge(1,0,request_reactions)]
            !
            if (request_stop) then
              request_array(1) = final_timestep_from_stop_request
            end if
            !
          end if
          !
          ! Have the root processor broadcast all the request logicals to
          ! the other processors so they can proceed as needed
          !
          call mpi_bcast(request_array,size(request_array,kind=int_mpi), &
                         mpi_inttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
          !
          final_timestep_from_stop_request = request_array(1)
          !
          request_stop             = (request_array(1) /= 0)
          request_restart          = (request_array(2) /= 0)
          request_solution         = (request_array(3) /= 0)
          request_cfl_adjustment   = (request_array(4) /= 0)
          request_restart_interupt = (request_array(5) /= 0)
          request_time_ave         = (request_array(6) /= 0)
          request_reactions        = (request_array(7) /= 0)
          !
          ! Exit if the STOP_GFR file was found
          !
          if (request_stop) then
            if (final_timestep_from_stop_request <= itcur) exit main_loop
          end if
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
            ! Reset the IEEE floating point status to the default values
            call reset_ieee_status
            CGNS_action = reopen_CGNS
          end if
          !
          if (request_cfl_adjustment) then
            !
            call mpi_bcast(cfl_tmp,size(cfl_tmp,kind=int_mpi),mpi_flttyp, &
                           0_int_mpi,MPI_COMM_WORLD,mpierr)
            !
           !cfl = cfl_tmp(1)
           !cfl_beg = cfl_tmp(2)
           !cfl_end = cfl_tmp(3)
           !cfl_cycles = nint(cfl_tmp(4))
           !cfl_start_iter = nint(cfl_tmp(5))
           !!
           !if (cfl_start_iter < itcur) then
           !  call adjust_cfl_cycles
           !end if
            call adjust_cfl_cycles(cfl_tmp)
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
      end if check_for_interupt_files
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
        if (is_at_restart_frequency) then
          output_restart = true
        end if
        !
        if (output_restart .and. .not.request_restart) then
          call write_restart_file(itcur)
          is_at_restart_frequency = fals
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
          ! Reset the IEEE floating point status to the default values
          call reset_ieee_status
        end if
        !
      end if
      !
      if (any(itestcase == Taylor_Problems)) then
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
          call write_restart_file(itcur)
          ! Reset the IEEE floating point status to the default values
          call reset_ieee_status
          write_TaylorGreen_full_solution = fals
        end if
      end if
      !
      ! If solving the Taylor-Green Vortex problem, compute
      ! the various integrated quantities needed to get the
      ! different the kinetic energy dissipation rates
      !
      if (output_KEDR .and. iter/=num_timesteps) then
        call compute_KEDR(itcur,dtmin)
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
      call solver_timer(leaving_procedure,pname,idx1)
      !
      call check_solver_timer()
      !
    end do main_loop
    !
    ! For now, always exit the restart loop
    !
    if (abs(CFL_reduction_for_NaN_restart) > one) then
      !
      exit restart_loop
      !
    else
      !
      if (num_CFL_restarts == max_CFL_restarts) then
        exit restart_loop
      end if
      !
      num_CFL_restarts = num_CFL_restarts + 1
      !
      ! Reset the IEEE floating point status to the default values
      !
      call reset_ieee_status
      !
    end if
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
  ! Write out a file providing a summary of the simulation
  !
  call write_summary_file(iter)
  !
  ! If solving the Taylor-Green Vortex problem,
  ! finalize the kinetic energy dissipation rates
  !
  if (num_timesteps > 0) then
    if (output_KEDR) then
      call compute_KEDR(-completed_timesteps,dtmin)
    end if
  end if
  !
  if (itcur /= 0) then
    call write_restart_file(-1)
  end if
  !
  if (continuous_output) then
    call new_write_parallel_cgns(close_CGNS)
  else
    call write_parallel_cgns(close_CGNS)
  end if
  ! Reset the IEEE floating point status to the default values
  call reset_ieee_status
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
               '  ncell = ',i0,t21,'  nnode = ',i0,/, &
               '  nfbnd = ',i0,t21,'  nfint = ',i0,t41,'  nface = ',i0,/, &
               '  nsbnd = ',i0,t21,'  npbnd = ',i0,t41,'  nsint = ',i0,//, &
               '  n_global_DOF  = ',i0,//, &
               '  n_solpts      = ',i0,/, &
               '  n_flxpts      = ',i0,/, &
               '  n_bnd_flxpts  = ',i0,/, &
               '  n_wall_flxpts = ',i0,//, &
               '  n_qdtpts      = ',i0,/, &
               '  n_flx_qdtpts  = ',i0,/, &
               '  n_bnd_qdtpts  = ',i0,/, &
               '  n_wall_qdtpts = ',i0,//, &
               '  Total Grid Volume = ',g16.8," [m^",i1,"]",/, &
               '  1D Solution Points = ',100(9f9.5,:,/,23x))
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
  199 format (/,"*******************************************************",/, &
                "  FOUND FILE REQUESTING TO STOP THE EXECUTION OF GFR!",/, &
                "  STOPPING EXECUTION AND DUMPING THE FINAL RESULTS",/, &
                "  AFTER ITERATION # ",i0," !",/, &
                "*******************************************************",/)
  !
  201 format (/,"INTERPOLATING SOLUTION TO NEW GRID.")
  202 format (/,"INTERPOLATION COMPLETE.")
  207 format (/,"WRITING RESTART FILE.")
  208 format (/,"WRITING CGNS FILE.")
  210 format ("STOPPING EXECUTION SO YOU MAY EXAMINE INTERPOLATED SOLUTION.")
  211 format ("STOPPING EXECUTION SO YOU MAY EXAMINE interpsol_debug.tec")
  !
contains
!
!###############################################################################
!
subroutine compute_rhs(iter,nst,newton_failures)
  !
  !.. Use Statements ..
  use ovar, only : governing_equations
  use ovar, only : visc_flux_method
  use ovar, only : mms_opt
  use ovar, only : itestcase
  !
  use bc_mod, only : get_faceusp_bc
  use bc_mod, only : get_facedusp_bc
  !
  use parallel_mod, only : exchange_faceusp
  use parallel_mod, only : wait_faceusp
  use parallel_mod, only : exchange_facedusp
  use parallel_mod, only : wait_facedusp
  !
  use mms_mod, only : mms_src
  !
  use channel_mod, only : channel_source
  !
  use reactions_mod, only : compute_wdot
  !
  use mrf_mod, only : mrf_src_cv
  !
  use rhs_mod, only : rresi
  !
  use flux_mod, only : interp_cv_to_flxpts
  use flux_mod, only : discontinuous_gradient_at_solpts
  use flux_mod, only : get_invs_common_fluxes_interior
  use flux_mod, only : correct_gradients_br2_gpu
  use flux_mod, only : correct_gradients_at_solpts_gpu
  use flux_mod, only : interp_dusp_to_flxpts
  use flux_mod, only : get_visc_common_fluxes_interior
  use flux_mod, only : get_invs_common_fluxes_bc
  use flux_mod, only : get_visc_common_fluxes_bc
  use flux_mod, only : flux_gradient_at_solpts
!!use flux_mod, only : correct_gradients_at_solpts
!!use flux_mod, only : correct_gradients_br2
  !
  !.. Formal Arguments ..
  integer,    intent(in) :: iter
  integer,    intent(in) :: nst
  integer, intent(inout) :: newton_failures
  !
  !.. Local Saved Scalars ..
  integer, save :: idx = -1
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_rhs"
  !
continue
  !
  call solver_timer(entering_procedure,pname,idx)
  !
  call global_debug_timer(entering_procedure,pname)
  !
  ! Interpolate the conserved variables from
  ! the solution points to the face/flux points
  !
 !call interp_cv_to_flxpts(iter,nst)
 !call debugging_checks_for_bad_number(__LINE__)
  !
  ! Initiate the exchange of face/flux point solutions
  ! between adjacent partitions
  !
  if (ncpu > 1) then
    call exchange_faceusp
    call debugging_checks_for_bad_number(__LINE__)
#ifdef IMMEDIATE_WAIT_FOR_MPI
    call wait_faceusp
    call debugging_checks_for_bad_number(__LINE__)
#endif
  end if
  !
  ! Compute the discontinuous gradient of the solution at the solution points
  !
 !if (governing_equations == NavierStokes_Eqns) then
 !  call discontinuous_gradient_at_solpts
 !  call debugging_checks_for_bad_number(__LINE__)
 !end if
  !
  ! Get the boundary conditions for the boundary face/flux points
  !
  call get_faceusp_bc(newton_failures)
  call debugging_checks_for_bad_number(__LINE__,newton_failures)
  !
  ! Compute the inviscid common fluxes at all the face/flux points that
  ! are not on a boundary face
  !
  call get_invs_common_fluxes_interior
  call debugging_checks_for_bad_number(__LINE__)
  !
  ! Ensure that the exchange of face/flux point solutions between adjacent
  ! partitions has been completed, because that information will be needed
  ! in the next step.
  !
  if (ncpu > 1) then
#ifndef IMMEDIATE_WAIT_FOR_MPI
    call wait_faceusp
    call debugging_checks_for_bad_number(__LINE__)
#endif
  end if
  !
  if (governing_equations == NavierStokes_Eqns) then
    !
    ! Compute the corrected gradients at the solution points and interpolate
    ! the corrected gradients to all the face/flux points
    !
    if (visc_flux_method == visc_flux_BR2) then
     !call correct_gradients_br2
      call correct_gradients_br2_gpu
      call debugging_checks_for_bad_number(__LINE__)
    else
      call correct_gradients_at_solpts_gpu
      call debugging_checks_for_bad_number(__LINE__)
     !call correct_gradients_at_solpts
     !call debugging_checks_for_bad_number(__LINE__)
      call interp_dusp_to_flxpts
      call debugging_checks_for_bad_number(__LINE__)
    end if
    !
    ! Initiate the exchange of face/flux point corrected
    ! gradients between adjacent partitions
    !
    if (ncpu > 1) then
      call exchange_facedusp
      call debugging_checks_for_bad_number(__LINE__)
#ifdef IMMEDIATE_WAIT_FOR_MPI
      call wait_facedusp
      call debugging_checks_for_bad_number(__LINE__)
#endif
    end if
    !
    ! Compute the gradient boundary conditions for the boundary face/flux points
    !
    call get_facedusp_bc
    call debugging_checks_for_bad_number(__LINE__)
    !
    ! Compute the viscous common fluxes at all the face/flux points that
    ! are not on a boundary face
    !
    call get_visc_common_fluxes_interior
    call debugging_checks_for_bad_number(__LINE__)
  end if
  !
  ! Compute the inviscid common fluxes for all boundary face/flux points
  !
  call get_invs_common_fluxes_bc
  call debugging_checks_for_bad_number(__LINE__)
  !
  if (governing_equations == NavierStokes_Eqns) then
    !
    ! Ensure that the exchange of face/flux point corrected gradients
    ! between adjacent partitions has been completed, because that
    ! information will be needed in the next step.
    !
    if (ncpu > 1) then
#ifndef IMMEDIATE_WAIT_FOR_MPI
      call wait_facedusp
      call debugging_checks_for_bad_number(__LINE__)
#endif
    end if
    !
    ! Compute the viscous common fluxes for all the boundary face/flux points
    !
    call get_visc_common_fluxes_bc
    call debugging_checks_for_bad_number(__LINE__)
    !
  end if
  !
  ! Compute the corrected flux divergence at the solution points
  ! NOTE: The name of this subroutine is a holdover from an older version
  !       and misleading, so it probably needs to be updated to what it
  !       actually does.
  !
  call flux_gradient_at_solpts(iter,nst)
  call debugging_checks_for_bad_number(__LINE__)
  !
  ! Compute the source terms for a moving reference frame
  !
  call mrf_src_cv
  call debugging_checks_for_bad_number(__LINE__)
  !
  ! Compute the source terms for method of manufactured solutions (MMS)
  !
  if (mms_opt /= 0) then
    call mms_src
    call debugging_checks_for_bad_number(__LINE__)
  end if
  !
  ! Compute the source terms for the channel flow
  !
  if (any(itestcase == [Channel_Flow,Pipe_Channel_Flow])) then
    call channel_source
    call debugging_checks_for_bad_number(__LINE__)
  end if
  !
  ! Compute the species source terms for chemically reacting flows
  !
  call compute_wdot(iter)
  call debugging_checks_for_bad_number(__LINE__)
  !
  ! Compute the final residual for each solution point
  !
  call rresi(iter,nst)
  call debugging_checks_for_bad_number(__LINE__)
  !
  call global_debug_timer(leaving_procedure,pname)
  !
  call solver_timer(leaving_procedure,pname,idx)
  !
end subroutine compute_rhs
!
!###############################################################################
!
#define DISABLE_DEBUGGING_CHECKS_FOR_BAD_NUMBER
subroutine debugging_checks_for_bad_number(line_num,newton_failures, &
                                           force_signaling)
  !
  ! This subroutine is a wrapper to quickly enable a high-frequency check for
  ! floating point exceptions within the solver. This will add significant
  ! overhead to the computational cost per iteration due to frequent use of
  ! MPI reductions and interpolations to find bad points.
  ! NOTE: If this is "disabled", it will still be enabled if newton_failures
  !       is present and non-zero to stop execution when a Newton iteration
  !       fails to converge.
  !
  !.. Formal Arguments ..
  integer, intent(in) :: line_num
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: newton_failures
  integer, optional, intent(in) :: force_signaling
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "debugging_checks_for_bad_number"
  !
continue
#ifndef DISABLE_DEBUGGING_CHECKS_FOR_BAD_NUMBER
  call check_for_bad_number(line_num,newton_failures,force_signaling)
#else
  if (present(newton_failures)) then
    if (newton_failures > 0) then
      call check_for_bad_number(line_num,newton_failures,force_signaling)
    end if
  end if
#endif
end subroutine debugging_checks_for_bad_number
!
!###############################################################################
!
subroutine check_for_bad_number(line_num,newton_failures,force_signaling)
  !
  !.. Use Statements ..
  use eqn_idx_mod, only : nq
  use geomvar_mod, only : n_solpts
  use flowvar_mod, only : usp,uoldsp
!#ifndef SPECIAL_FOR_GCC
!  use, intrinsic :: ieee_exceptions, only : ieee_usual
!  use, intrinsic :: ieee_exceptions, only : ieee_get_flag
!#endif
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
  integer  :: n,m,n_bad_pts,num_newton_failures
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
  logical(lk), parameter :: report_fpe_after_bad_point_confirmation = true
 !logical(lk), parameter :: report_fpe_after_bad_point_confirmation = fals
  !
continue
#ifndef DISABLE_CHECK_FOR_BAD_NUMBER
  !
  error_found = fals
  !
  overflow_error = ""
  divide_by_zero_error = ""
  invalid_error = ""
  newton_error = ""
  output_message = ""
  !
  num_newton_failures = 0
  if (present(newton_failures)) then
    num_newton_failures = newton_failures
  end if
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
  if (num_newton_failures > 0) then
    error_found = true
    write (newton_error,104) num_newton_failures
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
    if (.not. report_fpe_after_bad_point_confirmation) then
      call report_signaling_exception_found(line_num,overflow_error, &
                                            divide_by_zero_error, &
                                            invalid_error, &
                                            newton_error)
    end if
    !
    call find_bad_number(n_bad_pts)
    !
    if (ncpu > 1) then
      call mpi_allreduce(MPI_IN_PLACE,n_bad_pts,1_int_mpi, &
                         mpi_inttyp,MPI_MAX,MPI_COMM_WORLD,mpierr)
    end if
    !
    if (n_bad_pts > 0 .or. num_newton_failures > 0) then
      !
      if (report_fpe_after_bad_point_confirmation) then
        call report_signaling_exception_found(line_num,overflow_error, &
                                              divide_by_zero_error, &
                                              invalid_error, &
                                              newton_error)
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
      write (error_message,100) itcur,current_rk_stage,out_time, &
                                trim(adjustl(overflow_error)), &
                                trim(adjustl(divide_by_zero_error)), &
                                trim(adjustl(invalid_error)), &
                                trim(adjustl(newton_error))
      !
      call stop_gfr(stop_mpi,pname,line_num_with_bad_number, &
                    src_file_with_bad_number,error_message)
      !
    else
      !
      ! If are here, a floating point exception seemed to be found, but the
      ! density, pressure, and temperature all seem to be fine at both the
      ! solution points and face points. Reset the IEEE floating point status
      ! to the default values in case this was a false positive.
      !
     !call reset_ieee_status
      call reset_fp_exceptions
      !
    end if
    !
  end if
  !
  100 format ("A SIGNALING FLOATING POINT EXCEPTION WAS ", &
              "DETECTED AFTER UPDATING THE SOLUTION AFTER ", &
              "TIME STEP # ",i0,", RK STAGE # ",i0, &
              ", AT SOLUTION TIME ",es18.10, &
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
#endif
#endif
  !
end subroutine check_for_bad_number
!
!###############################################################################
!
subroutine timed_check_for_bad_number(line_num,newton_failures,force_signaling)
  !
  !.. Use Statements ..
  use eqn_idx_mod, only : nq
  use geomvar_mod, only : n_solpts
  use flowvar_mod, only : usp,uoldsp
!#ifndef SPECIAL_FOR_GCC
!  use, intrinsic :: ieee_exceptions, only : ieee_usual
!  use, intrinsic :: ieee_exceptions, only : ieee_get_flag
!#endif
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
  integer  :: n,m,n_bad_pts,num_newton_failures
  real(wp) :: out_time
  !
  logical(lk) :: error_found
  !
  character(len=300) :: timer_message
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
  character(len=*), parameter :: pname = "timed_check_for_bad_number"
  !
  logical(lk), parameter :: report_fpe_after_bad_point_confirmation = true
 !logical(lk), parameter :: report_fpe_after_bad_point_confirmation = fals
  !
continue
#ifndef DISABLE_CHECK_FOR_BAD_NUMBER
  !
  write (timer_message,10) pname,line_num,__FILE__
  10 format (a," called from line ",i0," of ",a)
  !
  call debug_timer(entering_procedure,timer_message)
  !
  error_found = fals
  !
  overflow_error = ""
  divide_by_zero_error = ""
  invalid_error = ""
  newton_error = ""
  output_message = ""
  !
  num_newton_failures = 0
  if (present(newton_failures)) then
    num_newton_failures = newton_failures
  end if
  !
  ! First check to make sure that no floating-point exceptions are signaling
  !
  call debug_timer(start_timer, &
                   "ieee_usual in timed_check_for_bad_number(main)")
  call ieee_get_flag(ieee_usual,signaling_exceptions)
  call debug_timer(stop_timer, &
                   "ieee_usual in timed_check_for_bad_number(main)")
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
    call debug_timer(start_timer, &
                     "mpi_allreduce in timed_check_for_bad_number(main)")
    call mpi_allreduce(MPI_IN_PLACE,bad_number_found, &
                       size(bad_number_found,kind=int_mpi), &
                       mpi_inttyp,MPI_MAX,MPI_COMM_WORLD,mpierr)
    call debug_timer(stop_timer, &
                     "mpi_allreduce in timed_check_for_bad_number(main)")
  end if
  !
  if (num_newton_failures > 0) then
    error_found = true
    write (newton_error,104) num_newton_failures
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
    if (.not. report_fpe_after_bad_point_confirmation) then
      call report_signaling_exception_found(line_num,overflow_error, &
                                            divide_by_zero_error, &
                                            invalid_error, &
                                            newton_error)
    end if
    !
    call find_bad_number(n_bad_pts)
    !
    if (ncpu > 1) then
      call mpi_allreduce(MPI_IN_PLACE,n_bad_pts,1_int_mpi, &
                         mpi_inttyp,MPI_MAX,MPI_COMM_WORLD,mpierr)
    end if
    !
    if (n_bad_pts > 0 .or. num_newton_failures > 0) then
      !
      if (report_fpe_after_bad_point_confirmation) then
        call report_signaling_exception_found(line_num,overflow_error, &
                                              divide_by_zero_error, &
                                              invalid_error, &
                                              newton_error)
      end if
      !
      flush (iout)
      call mpi_barrier(MPI_COMM_WORLD,mpierr)
      flush (iout)
      call mpi_barrier(MPI_COMM_WORLD,mpierr)
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
      write (error_message,100) itcur,current_rk_stage,out_time, &
                                trim(adjustl(overflow_error)), &
                                trim(adjustl(divide_by_zero_error)), &
                                trim(adjustl(invalid_error)), &
                                trim(adjustl(newton_error))
      !
      call stop_gfr(stop_mpi,pname,line_num_with_bad_number, &
                    src_file_with_bad_number,error_message)
      !
    else
      !
      ! If are here, a floating point exception seemed to be found, but the
      ! density, pressure, and temperature all seem to be fine at both the
      ! solution points and face points. Reset the IEEE floating point status
      ! to the default values in case this was a false positive.
      !
     !call reset_ieee_status
      call reset_fp_exceptions
      !
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,timer_message)
  !
  100 format ("A SIGNALING FLOATING POINT EXCEPTION WAS ", &
              "DETECTED AFTER UPDATING THE SOLUTION AFTER ", &
              "TIME STEP # ",i0,", RK STAGE # ",i0, &
              ", AT SOLUTION TIME ",es18.10, &
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
#endif
#endif
  !
end subroutine timed_check_for_bad_number
!
!###############################################################################
!
subroutine report_signaling_exception_found(line_num,overflow_error, &
                                            divide_by_zero_error, &
                                            invalid_error,newton_error)
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: line_num
  character(len=*), intent(in) :: overflow_error
  character(len=*), intent(in) :: divide_by_zero_error
  character(len=*), intent(in) :: invalid_error
  character(len=*), intent(in) :: newton_error
  !
  !.. Local Scalars ..
  character(len=10000) :: output_message
  !
continue
  !
  if (len_trim(adjustl(src_file_with_bad_number)) == 0) then
    src_file_with_bad_number = __FILE__
    line_num_with_bad_number = line_num
  end if
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
    write (error_message,200) itcur,current_rk_stage
    output_message = reformat_string(error_message, &
                                     line_length=70, &
                                     leading_blanks=5)
    !
    write (iout,1) trim(output_message)
    flush (iout)
    !
    write (iout,202) line_num_with_bad_number, &
                     trim(adjustl(src_file_with_bad_number))
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
  ! Format Statements
  !
  1 format (//,1x,79("#"),//,5x,a,/)
  2 format (5x,"ERROR MESSAGE: '",a,"'",//,1x,79("#"),//)
  !
  200 format ("A SIGNALING FLOATING POINT EXCEPTION WAS ", &
              "DETECTED AFTER UPDATING THE SOLUTION AFTER ", &
              "TIME STEP # ",i0,", RK STAGE # ",i0," !!! ", &
              "DUMPING SOLUTION AND RESTART FILES FOR THE ", &
              "PREVIOUS TIME STEP!!!")
  201 format (4(2x,a,:))
  202 format (5x,"FLOATING POINT CHECK WAS CALLED ON LINE #  ",i0,/,&
              35x,"OF THE FILE '",a,"'",/)
  !
end subroutine report_signaling_exception_found
!
!###############################################################################
!
subroutine find_bad_number(n_bad_pts)
  !
  !.. Use Statements ..
  use eqn_idx_mod,       only : nq,nec,ntemp,npres
  use ovar,              only : grid_scaling_factor,rhoref,aref,tref
  use ovar,              only : itestcase
  use geomvar_mod,       only : nr,cell,face,fp,nfbnd,nface
  use geomvar_mod,       only : n_solpts,n_bnd_flxpts,n_flxpts,xyz
  use flowvar_mod,       only : usp,gpsp,faceusp,facegpsp
  use flowvar_mod,       only : facexyz
  use interpolation_mod, only : interp
  !
 !use, intrinsic :: ieee_arithmetic, only : is_finite => ieee_is_finite
 !use, intrinsic :: ieee_arithmetic, only : is_negative => ieee_is_negative
  !
  !.. Formal Arguments ..
  integer, intent(inout) :: n_bad_pts
  !
  !.. Local Scalars ..
  integer  :: n,nc,nf,nfp,n0,sf,f0,i,i1,i2,j,kf,lf,l,ic
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  integer  :: one_pass
  real(wp) :: a
  !
  logical(lk) :: use_interp
  !
  character(len=15)  :: cxyz_out
  character(len=200) :: inv_fmt,neg_fmt
  !
  !.. Local Arrays ..
  real(wp) :: vars(1:3)
  real(wp) :: xyz_fp(1:nr)
  real(wp) :: xyz_out(1:nr)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: xyz_of_bad_pts(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_bad_number"
  !
continue
  !
  n_bad_pts = 0
 !xyz_of_bad_pts = init_array(5,1000,zero) ! F2003 AUTO-REALLOC
  !
  one_pass_solpt_loop: &
  do one_pass = 1,1
    !
    ! Skip this loop if any of the needed arrays have yet to be allocated
    !
    if (.not. allocated(xyz) .or. &
        .not. allocated(usp) .or. &
        .not. allocated(gpsp)) exit one_pass_solpt_loop
    !
    ! Find a negative density, temperature, or pressure at a solution point
    !
    write (inv_fmt,1) "Invalid",nr
    write (neg_fmt,1) "Negative",nr
    !
   !do n = 1,n_solpts
    do nc = 1,ncell
      do n = cell(nc)%beg_sp,cell(nc)%end_sp
        vars(1) = usp(n,nec)
        vars(2) = gpsp(n,npres)
        vars(3) = gpsp(n,ntemp)
        if (.not. all(is_finite(vars))) then
          call add_bad_pt(n_bad_pts,1,nc,xyz(:,n),xyz_of_bad_pts)
          call get_xyz_out(xyz(:,n),xyz_out,cxyz_out)
         !write (iout,inv_fmt) (xyz(l,n)/grid_scaling_factor, l=1,nr), &
          write (iout,inv_fmt) trim(adjustl(cxyz_out)), (xyz_out(l), l=1,nr), &
                               vars(1)*rhoref, &
                               vars(2)*rhoref*aref*aref, &
                               vars(3)*tref
          flush (iout)
        else if (any(is_negative(vars))) then
          call add_bad_pt(n_bad_pts,1,nc,xyz(:,n),xyz_of_bad_pts)
          call get_xyz_out(xyz(:,n),xyz_out,cxyz_out)
         !write (iout,neg_fmt) (xyz(l,n)/grid_scaling_factor, l=1,nr), &
          write (iout,neg_fmt) trim(adjustl(cxyz_out)), (xyz_out(l), l=1,nr), &
                               vars(1)*rhoref, &
                               vars(2)*rhoref*aref*aref, &
                               vars(3)*tref
          flush (iout)
        end if
      end do
    end do
    !
  end do one_pass_solpt_loop
  !
  ! Determine if we need to use the interpolation matrix
  ! to get the coordinates of the face points
  !
  use_interp = allocated(interp)
  !
  if (allocated(facexyz)) then
    if (size(facexyz,dim=2) == n_flxpts) then
      use_interp = fals
    end if
  end if
  !
  ! Find a negative density, temperature, or pressure at a face point
  !
  if (use_interp) then
    !
    one_pass_interp_face_loop: &
    do one_pass = 1,1
      !
      ! Skip this loop if any of the needed arrays have yet to be allocated
      !
      if (.not. allocated(xyz) .or. &
          .not. allocated(cell) .or. &
          .not. allocated(face) .or. &
          .not. allocated(fp) .or. &
          .not. allocated(faceusp) .or. &
          .not. allocated(facegpsp)) exit one_pass_interp_face_loop
      !
      ! Find a negative density, temperature,
      ! or pressure at a boundary face point
      !
      write (inv_fmt,2) "Invalid boundary",nr
      write (neg_fmt,2) "Negative boundary",nr
      !
      bnd_nf_loop: &
      do nf = 1,nfbnd
        !
        face_geom = face(nf)%geom
        face_order = face(nf)%order
        !
        f0 = face(nf)%beg_fp - 1
        nfp = face(nf)%end_fp - f0
        !
        bnd_kf_loop: &
        do kf = 1,nfp
          bnd_sf_loop: &
          do sf = 1,2
            !
            vars(1) = faceusp(f0+kf,nec,sf)
            vars(2) = facegpsp(f0+kf,npres,sf)
            vars(3) = facegpsp(f0+kf,ntemp,sf)
            !
            if (.not. all(is_finite(vars)) .or. &
                any(is_negative(vars))) then
              !
              nc = face(nf)%left%cell
              !
              this_geom = cell(nc)%geom
              this_order = cell(nc)%order
              !
              n0 = cell(nc)%beg_sp - 1
              !
              lf = face(nf)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                                rot(face(nf)%left%rotation)% &
                                             pts(kf)
              !
              i1 = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1
              i2 = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
              !
              xyz_fp(1:nr) = zero
              do i = i1,i2
                j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
                a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
                xyz_fp(1:nr) = xyz_fp(1:nr) + a * xyz(1:nr,n0+j)
              end do
              !
              if (.not. all(is_finite(vars))) then
                call add_bad_pt(n_bad_pts,2,nc,xyz_fp(:),xyz_of_bad_pts)
                call get_xyz_out(xyz_fp,xyz_out,cxyz_out)
                write (iout,inv_fmt) sf, &
                                    !(xyz_fp(l)/grid_scaling_factor, l=1,nr), &
                                     trim(adjustl(cxyz_out)), &
                                     (xyz_out(l), l=1,nr), &
                                     vars(1)*rhoref, &
                                     vars(2)*rhoref*aref*aref, &
                                     vars(3)*tref
                flush (iout)
              else if (any(is_negative(vars))) then
                call add_bad_pt(n_bad_pts,2,nc,xyz_fp(:),xyz_of_bad_pts)
                call get_xyz_out(xyz_fp,xyz_out,cxyz_out)
                write (iout,neg_fmt) sf, &
                                    !(xyz_fp(l)/grid_scaling_factor, l=1,nr), &
                                     trim(adjustl(cxyz_out)), &
                                     (xyz_out(l), l=1,nr), &
                                     vars(1)*rhoref, &
                                     vars(2)*rhoref*aref*aref, &
                                     vars(3)*tref
                flush (iout)
              end if
              !
            end if
          end do bnd_sf_loop
        end do bnd_kf_loop
      end do bnd_nf_loop
      !
      ! Find a negative density, temperature,
      ! or pressure at an interior face point
      !
      write (inv_fmt,2) "Invalid interior",nr
      write (neg_fmt,2) "Negative interior",nr
      !
      int_nf_loop: &
      do nf = nfbnd+1,nface
        !
        face_geom = face(nf)%geom
        face_order = face(nf)%order
        !
        f0 = face(nf)%beg_fp - 1
        nfp = face(nf)%end_fp - f0
        !
        int_kf_loop: &
        do kf = 1,nfp
          int_sf_loop: &
          do sf = 1,2
            !
            vars(1) = faceusp(f0+kf,nec,sf)
            vars(2) = facegpsp(f0+kf,npres,sf)
            vars(3) = facegpsp(f0+kf,ntemp,sf)
            !
            if (.not. all(is_finite(vars)) .or. &
                any(is_negative(vars))) then
              !
              nc = face(nf)%left%cell
              !
              this_geom = cell(nc)%geom
              this_order = cell(nc)%order
              !
              n0 = cell(nc)%beg_sp - 1
              !
              lf = face(nf)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                                rot(face(nf)%left%rotation)% &
                                             pts(kf)
              !
              i1 = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1
              i2 = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
              !
              xyz_fp(1:nr) = zero
              do i = i1,i2
                j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
                a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
                xyz_fp(1:nr) = xyz_fp(1:nr) + a * xyz(1:nr,n0+j)
              end do
              !
              if (.not. all(is_finite(vars))) then
                ic = merge(face(nf)%left%cell,face(nf)%right%cell,sf==1)
                call add_bad_pt(n_bad_pts,3,ic,xyz_fp(:),xyz_of_bad_pts)
                call get_xyz_out(xyz_fp,xyz_out,cxyz_out)
                write (iout,inv_fmt) sf, &
                                    !(xyz_fp(l)/grid_scaling_factor, l=1,nr), &
                                     trim(adjustl(cxyz_out)), &
                                     (xyz_out(l), l=1,nr), &
                                     vars(1)*rhoref, &
                                     vars(2)*rhoref*aref*aref, &
                                     vars(3)*tref
                flush (iout)
              else if (any(is_negative(vars))) then
                ic = merge(face(nf)%left%cell,face(nf)%right%cell,sf==1)
                call add_bad_pt(n_bad_pts,3,ic,xyz_fp(:),xyz_of_bad_pts)
                call get_xyz_out(xyz_fp,xyz_out,cxyz_out)
                write (iout,neg_fmt) sf, &
                                    !(xyz_fp(l)/grid_scaling_factor, l=1,nr), &
                                     trim(adjustl(cxyz_out)), &
                                     (xyz_out(l), l=1,nr), &
                                     vars(1)*rhoref, &
                                     vars(2)*rhoref*aref*aref, &
                                     vars(3)*tref
                flush (iout)
              end if
              !
            end if
          end do int_sf_loop
        end do int_kf_loop
      end do int_nf_loop
      !
    end do one_pass_interp_face_loop
    !
  else
    !
    one_pass_noninterp_face_loop: &
    do one_pass = 1,1
      !
      ! Skip this loop if any of the needed arrays have yet to be allocated
      !
      if (.not. allocated(facexyz) .or. &
          .not. allocated(faceusp) .or. &
          .not. allocated(facegpsp)) exit one_pass_noninterp_face_loop
      !
      ! Find a negative density, temperature,
      ! or pressure at a boundary face point
      !
      write (inv_fmt,2) "Invalid boundary",nr
      write (neg_fmt,2) "Negative boundary",nr
      !
     !do n = 1,n_bnd_flxpts
      do nf = 1,nfbnd
        !
        f0 = face(nf)%beg_fp - 1
        nfp = face(nf)%end_fp - f0
        !
        nc = face(nf)%left%cell
        !
        do kf = 1,nfp
          n = f0 + kf
          do sf = 1,2
            vars(1) = faceusp(n,nec,sf)
            vars(2) = facegpsp(n,npres,sf)
            vars(3) = facegpsp(n,ntemp,sf)
            if (.not. all(is_finite(vars))) then
              call add_bad_pt(n_bad_pts,2,nc,facexyz(:,n),xyz_of_bad_pts)
              call get_xyz_out(facexyz(:,n),xyz_out,cxyz_out)
              write (iout,inv_fmt) sf, &
                              !(facexyz(l,n)/grid_scaling_factor, l=1,nr), &
                               trim(adjustl(cxyz_out)), (xyz_out(l), l=1,nr), &
                               vars(1)*rhoref, &
                               vars(2)*rhoref*aref*aref, &
                               vars(3)*tref
              flush (iout)
            else if (any(is_negative(vars))) then
              call add_bad_pt(n_bad_pts,2,nc,facexyz(:,n),xyz_of_bad_pts)
              call get_xyz_out(facexyz(:,n),xyz_out,cxyz_out)
              write (iout,neg_fmt) sf, &
                              !(facexyz(l,n)/grid_scaling_factor, l=1,nr), &
                               trim(adjustl(cxyz_out)), (xyz_out(l), l=1,nr), &
                               vars(1)*rhoref, &
                               vars(2)*rhoref*aref*aref, &
                               vars(3)*tref
              flush (iout)
            end if
          end do
        end do
      end do
      !
      ! Find a negative density, temperature,
      ! or pressure at an interior face point
      !
      write (inv_fmt,2) "Invalid interior",nr
      write (neg_fmt,2) "Negative interior",nr
      !
     !do n = n_bnd_flxpts+1,n_flxpts
      do nf = nfbnd+1,nface
        !
        f0 = face(nf)%beg_fp - 1
        nfp = face(nf)%end_fp - f0
        !
        do kf = 1,nfp
          n = f0 + kf
          do sf = 1,2
            nc = merge(face(nf)%left%cell,face(nf)%right%cell,sf==1)
            vars(1) = faceusp(n,nec,sf)
            vars(2) = facegpsp(n,npres,sf)
            vars(3) = facegpsp(n,ntemp,sf)
            if (.not. all(is_finite(vars))) then
              call add_bad_pt(n_bad_pts,3,nc,facexyz(:,n),xyz_of_bad_pts)
              call get_xyz_out(facexyz(:,n),xyz_out,cxyz_out)
              write (iout,inv_fmt) sf, &
                              !(facexyz(l,n)/grid_scaling_factor, l=1,nr), &
                               trim(adjustl(cxyz_out)), (xyz_out(l), l=1,nr), &
                               vars(1)*rhoref, &
                               vars(2)*rhoref*aref*aref, &
                               vars(3)*tref
              flush (iout)
            else if (any(is_negative(vars))) then
              call add_bad_pt(n_bad_pts,3,nc,facexyz(:,n),xyz_of_bad_pts)
              call get_xyz_out(facexyz(:,n),xyz_out,cxyz_out)
              write (iout,neg_fmt) sf, &
                              !(facexyz(l,n)/grid_scaling_factor, l=1,nr), &
                               trim(adjustl(cxyz_out)), (xyz_out(l), l=1,nr), &
                               vars(1)*rhoref, &
                               vars(2)*rhoref*aref*aref, &
                               vars(3)*tref
              flush (iout)
            end if
          end do
        end do
      end do
      !
    end do one_pass_noninterp_face_loop
    !
  end if
  !
  if (any(abs(itestcase) == [Channel_Flow,Pipe_Channel_Flow])) then
    call report_bad_pts(n_bad_pts,xyz_of_bad_pts)
  end if
  !
  ! Format Statements
  !
 !1 format ('(2x,"',a,' SP value found at [x,y,z] = ",', &
  1 format ('(2x,"',a,' SP value found at [",a,"] = ",', &
            i0,'es14.6,4x,"r,p,t = ",3es14.6)')
 !2 format ('(2x,"',a,' FP value found at side ",i1," of [x,y,z] = ",', &
  2 format ('(2x,"',a,' FP value found at side ",i1," of [",a,"] = ",', &
            i0,'es14.6,4x,"r,p,t = ",3es14.6)')
  !
end subroutine find_bad_number
!
!###############################################################################
!
subroutine get_xyz_out(xyz,xyz_out,cxyz_out)
  !
  !.. Use Statements ..
  use eqn_idx_mod, only : NCD_Axial,NCD_Vert,NCD_Span
  use geomvar_mod, only : nr
  use ovar,        only : grid_scaling_factor,itestcase
  !
  !.. Formal Arguments ..
  real(wp),            intent(in) :: xyz(:)
  real(wp),         intent(inout) :: xyz_out(:)
  character(len=*), intent(inout) :: cxyz_out
  !
  !.. Local Scalars ..
  real(wp)    :: x,y,z,r,theta
  logical(lk) :: report_axisymmetric_xyz
  !
continue
  !
  report_axisymmetric_xyz = fals
  !
  if (nr < 3) then
    report_axisymmetric_xyz = fals
  else if (any(abs(itestcase) == [Subsonic_Jet, &
                                  Supersonic_Jet, &
                                  Taylor_Couette_Flow, &
                                  Pipe_Channel_Flow, &
                                  UIUC_Scramjet])) then
    report_axisymmetric_xyz = true
  end if
  !
  if (report_axisymmetric_xyz) then
    !
    x = xyz(NCD_Axial)/grid_scaling_factor
    y = xyz(NCD_Vert )/grid_scaling_factor
    z = xyz(NCD_Span )/grid_scaling_factor
    !
    r = sqrt(y*y + z*z)
    theta = atan2(z,y)
    if (theta < zero) theta = theta + two*pi
    theta = theta*radians2degrees
    !
    xyz_out(1:3) = [x,r,theta]
    cxyz_out = "Axial,R,Theta"
    !
  else
    !
    xyz_out(1:nr) = xyz(1:nr)/grid_scaling_factor
    cxyz_out = "x,y,z"
    !
  end if
  !
end subroutine get_xyz_out
!
!###############################################################################
!
subroutine add_bad_pt(n,loc,nc,xyz,xyz_bad_pts)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  !
  !.. Formal Arguments ..
  integer,               intent(inout) :: n
  integer,                  intent(in) :: loc
  integer,                  intent(in) :: nc
  real(wp),                 intent(in) :: xyz(:)
  real(wp), allocatable, intent(inout) :: xyz_bad_pts(:,:)
  !
  !.. Local Scalars ..
  integer :: l
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: tmp(:,:)
  !
continue
  !
  if (.not. allocated(xyz_bad_pts)) then
    xyz_bad_pts = init_array(5,1,zero) ! F2003 AUTO-REALLOC
  end if
  !
  n = n + 1
  !
  l = size(xyz_bad_pts,dim=2)
  !
  if (n > l) then
    tmp = init_array(5,n,zero) ! F2003 AUTO-REALLOC
    tmp(1:5,1:l) = xyz_bad_pts(1:5,1:l)
    call move_alloc( from=tmp , to=xyz_bad_pts )
  end if
  !
  xyz_bad_pts(1:nr,n) = xyz(1:nr)
  xyz_bad_pts(4,n) = real(loc,kind=wp)
  xyz_bad_pts(5,n) = real(nc,kind=wp)
  !
end subroutine add_bad_pt
!
!###############################################################################
!
subroutine report_bad_pts(npts,xyz)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr,xyz_nodes,nodes_of_cell_ptr,nodes_of_cell
  use ovar, only : yplus1 => channel_yplus_of_one
  use ovar, only : itestcase,grid_scaling_factor
  use internal_grid_mod, only : xplus_cell_ave
  use internal_grid_mod, only : zplus_cell_ave
  use internal_grid_mod, only : yplus_cell_wall
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: npts
  real(wp), intent(in) :: xyz(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,loc,nc,k1,k2,kp,k,l,i
  real(wp) :: xp_frst_cell,zp_frst_cell,yp_bot_wall
  real(wp) :: xp_last_cell,zp_last_cell,yp_top_wall
  !
  !.. Local Arrays ..
  real(wp) :: xyz_cell(1:nr)
  real(wp) :: xp_cell(1:nr)
  logical(lk) :: bnd_cell(1:2,1:nr)
  !
continue
  !
  if (any(abs(itestcase) == [Channel_Flow,Pipe_Channel_Flow])) then
    !
    xp_frst_cell = xplus_cell_ave
    xp_last_cell = (channel_grid_Lx/yplus1) - xplus_cell_ave
    !
    zp_frst_cell = zplus_cell_ave
    zp_last_cell = (channel_grid_Lz/yplus1) - zplus_cell_ave
    !
    yp_bot_wall = yplus_cell_wall
    yp_top_wall = (channel_grid_Ly/yplus1) - yplus_cell_wall
    !
  end if
  !
  do n = 1,npts
    !
    loc = nint(xyz(4,n))
    nc = nint(xyz(5,n))
    !
    xyz_cell(:) = zero
    !
    k1 = nodes_of_cell_ptr(nc)+1
    k2 = nodes_of_cell_ptr(nc+1)
    kp = k2-k1+1
    do k = k1,k2
      i = nodes_of_cell(k)
      do l = 1,nr
        xyz_cell(l) = xyz_cell(l) + xyz_nodes(l,i)/real(kp,kind=wp)
      end do
    end do
    !
    xyz_cell(:) = xyz_cell(:)/grid_scaling_factor
    !
    if (any(abs(itestcase) == [Channel_Flow,Pipe_Channel_Flow])) then
      !
      xp_cell(:) = xyz_cell(:)/yplus1
      !
      bnd_cell(1,1) = ( xp_cell(1) < xp_frst_cell )
      bnd_cell(2,1) = ( xp_cell(1) > xp_last_cell )
      !
      bnd_cell(1,2) = ( xp_cell(2) < yp_bot_wall )
      bnd_cell(2,2) = ( xp_cell(2) > yp_top_wall )
      !
      bnd_cell(1,3) = ( xp_cell(3) < zp_frst_cell )
      bnd_cell(2,3) = ( xp_cell(3) > zp_last_cell )
      !
      xyz_cell(1) = ten*ten*xyz_cell(1)/channel_grid_Lx
      xyz_cell(2) = ten*ten*xyz_cell(2)/channel_grid_Ly
      xyz_cell(3) = ten*ten*xyz_cell(3)/channel_grid_Lz
      !
      write (iout,1) (xp_cell(l), l=1,nr), (xyz_cell(l), l=1,nr)
      !
      if (bnd_cell(1,1)) write (iout,2) "first cell","x"
      if (bnd_cell(2,1)) write (iout,2) "last  cell","x"
      if (bnd_cell(1,2)) write (iout,2) "first cell on bottom wall","y"
      if (bnd_cell(2,2)) write (iout,2) "first cell on top    wall","y"
      if (bnd_cell(1,3)) write (iout,2) "first cell","z"
      if (bnd_cell(2,3)) write (iout,2) "last  cell","z"
      !
    end if
    !
   !write (iout,*) (xyz(l,n), l=1,nr)
    !
  end do
  !
  1 format (/,2x,"For cell that encountered a bad number:",/, &
              4x,"Cell center coordinates in       +-coordinates: ", &
                 "[x+,y+,z+] = [",f9.3,",",f9.3,",",f9.3,"]",/, &
              4x,"Cell center coordinates in   % of grid extents: ", &
                 "[x%,y%,z%] = [",f8.2,"%,",f8.2,"%,",f8.2,"%]")
  2 format (4x,"Bad number was found in ",a," of ",a," direction.")
  !
end subroutine report_bad_pts
!
!###############################################################################
!
!subroutine reset_fp_exceptions
!  !
!#ifndef SPECIAL_FOR_GCC
!  use, intrinsic :: ieee_exceptions, only : ieee_usual
!  use, intrinsic :: ieee_exceptions, only : ieee_set_flag
!  !
!  !.. Local Scalars ..
!  integer :: n
!  !
!  !.. Local Arrays ..
!  logical(ldk) :: reset_exceptions(1:size(ieee_usual))
!  !
!  character(len=*), parameter :: pname = "reset_fp_exceptions"
!  !
!continue
!  !
!  do n = 1,size(reset_exceptions)
!    reset_exceptions(n) = .false.
!  end do
!  !
!  ! First check to make sure that no floating-point exceptions are signaling
!  !
!  call ieee_set_flag(ieee_usual,reset_exceptions)
!#endif
!  !
!end subroutine reset_fp_exceptions
!
!###############################################################################
!
subroutine report_memory_usage(n,request_stop)
  !
  !.. Use Statements ..
  use geomvar_mod,       only : geomvar_memory_usage
  use ovar,              only : ovar_memory_usage
  use flowvar_mod,       only : flowvar_memory_usage
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
  use limiters_mod,      only : limiters_memory_usage
  use bnd_profiles_mod,  only : bnd_profiles_memory_usage
  use plot3d_mod,        only : plot3d_memory_usage
  use gmsh_mod,          only : gmsh_memory_usage
  use cgns_mod,          only : cgns_memory_usage
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
  call geomvar_memory_usage(i)
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
!###############################################################################
!
!!#define PROFILE_ON
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end program GFR_main
!
!###############################################################################
!
