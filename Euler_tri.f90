program Euler_tri
  !
  use module_kind_types
  !
  use geovar, only : ncell,n_global_cell
  use geovar, only : nnode,n_global_node
  use geovar, only : nfbnd,n_global_fbnd
  use geovar, only : n_solpts,n_global_solpts
  use geovar, only : bface,xyz_nodes
  use geovar, only : cell_geom,cell_order,global_cell_order
  use geovar, only : nodes_of_cell_ptr,nodes_of_cell
  use geovar, only : geovar_memory_usage
  !
  use ovar, only : beg_wall_time,end_wall_time,max_wall_time
  use ovar, only : work_units,wall_time,walltime_expiring
  use ovar, only : vortex_has_collapsed,Density_Maximum_Residual
  use ovar, only : final_countdown,Period_Timesteps
  use ovar, only : metis_option,convergence_reached
  use ovar, only : total_volume,local_volume,itestcase
  use ovar, only : nst,num_stages,iRKtype,time,itcur,num_ts
  use ovar, only : governing_equations,visc_flux_method
  use ovar, only : results_interval,iter_out_interval
  use ovar, only : restart_interval,output_interval
  use ovar, only : dump_memory_usage
  !
  use flowvar, only : usp,uoldsp
  !
  use quadrature_mod, only : solpts_edge
  use quadrature_mod, only : quadrature_memory_usage
  !
  use io_mod, only : read_input_file
  use io_mod, only : write_restart_file
  use io_mod, only : write_final_error
  use io_mod, only : write_nodal_parallel_cgns
  use io_mod, only : write_results_file
  use io_mod, only : check_flux_point_orderings
  !
  use mappings_mod, only : get_solution_points
  !
  use connectivity_mod, only : get_grid_connectivity
  !
  use metrics_mod, only : get_metrics
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
  !
  use module_limiters, only : initialize_limiters
  use module_limiters, only : compute_cell_ave
  !
  use iso_c_binding, only : c_long
  !
  implicit none
  !
  character(len=*), parameter :: pname = "Euler_tri(main)"
  !
  integer  :: iter,i,ierr,k
  real(wp) :: dtmin
  !
  integer, parameter :: iostp = 2683
  !
  logical :: request_stop
  !
 !logical, parameter :: ignore_errors = true
  logical, parameter :: ignore_errors = fals
  !
  integer, parameter :: test_param = merge( 1 , 2 , ignore_errors )
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
  beg_wall_time = mpi_wtime()
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
    write (iout,16) current_wtl
    max_wall_time = real(current_wtl,kind=wp)
  end if
  !
  ! Broadcast the maximum wall time to all processors
  !
  call mpi_bcast(max_wall_time,1,mpi_flttyp,0,mpi_comm_world,mpierr)
  !
  ! Initialize all the random constants that must be initialized at run time
  !
  call initialize_runtime_parameters
  !
  ! Initialize the Vortex collapsing stuff
  !
  vortex_has_collapsed = fals
  Density_Maximum_Residual = -huge(zero)
  final_countdown = 0
  !
  ! Read in the input file
  !
  call read_input_file
  !
  ! Partition the grid if using more than one processor. 
  ! Upon entry:  These arrays define the global grid
  ! Upon return: These arrays define the partition local
  !              to the respective processor along with
  !              connectivities for the partition boundaries
  !
  if (ncpu > 1) then
    call partition_grid(ncpu,metis_option,bface,nodes_of_cell, &
                        nodes_of_cell_ptr,xyz_nodes,cell_geom,cell_order)
  else
    call create_serial_cell_map(ncell)
  end if
  !
  ! Get the coordinates of the solution points within each grid cell
  !
  call get_solution_points
  !
  ! Get the grid connectivity and metrics
  !
  call get_grid_connectivity
  !
  ! Compute the remaining quadrature information
  ! (Vandermonde matrix, stiffness matrix, LIFT / correction function, etc.)
  !
  call init_vandermonde_matrices
  call init_derivative_matrices
  call init_correction_matrices
  call init_interpolation_matrices
  call init_outerpolation_matrices
  call init_projection_matrices
  call init_filter_matrices
  !
  ! If requested, dump the memory usage of the advanced data types from some of
  ! the modules
  !
  if (dump_memory_usage) then
    !
    open (newunit=i,file="memory_usage.dat",status="replace", &
          action="write",form="formatted",position="rewind",iostat=ierr)
    call io_error(pname,"memory_usage.dat",1,ierr)
    !
    call geovar_memory_usage(i)
    call vandermonde_memory_usage(i)
    call derivatives_memory_usage(i)
    call correction_memory_usage(i)
    call interpolation_memory_usage(i)
    call quadrature_memory_usage(i)
    call projection_memory_usage(i)
    call filter_memory_usage(i)
    !
    close (i,iostat=ierr)
    call io_error(pname,"memory_usage.dat",2,ierr)
    !
  end if
  !
  ! Get the grid metrics
  !
  call get_metrics
  !
  ! Output information about each partition
  !
  do i = 0,ncpu-1
    if (mypnum == i) then
      write (iout,10) mypnum
      write (iout,12) ncell, nnode, nfbnd, n_solpts
      write (iout,13) local_volume
    end if
    call mpi_barrier(mpi_comm_world,mpierr)
  end do
  call mpi_barrier(mpi_comm_world,mpierr)
  !
  ! Output information about the global grid
  !
  if (mypnum == 0) then
    write (iout,11)
    write (iout,12) n_global_cell,n_global_node,n_global_fbnd,n_global_solpts
    write (iout,14) total_volume
    write (iout,15) (solpts_edge(k), k=1,size(solpts_edge))
  end if
  call mpi_barrier(mpi_comm_world,mpierr)
  !
  ! Initialize the flow variables for all solution points or
  ! read in the restart values if restarting a simulation
  !
  call initialize_flow
  !
  ! Output the initialized solution
  !
  call write_nodal_parallel_cgns(1)
  !
  work_units = zero
  if (mypnum == 0) call write_results_file(1,itcur)
  !
  if (mypnum == 0) then
    write (iout,19) minval(global_cell_order), maxval(global_cell_order)
  end if
  !
  call check_flux_point_orderings
  !
  ! Initialize any necessary information for use in the module generic_mod
  !
  call init_generic_mod
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%                                    %%%
  !%%%  Beginning of main iteration loop  %%%
  !%%%                                    %%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
 !call initialize_limiters
 !call compute_cell_ave
  !
  !
  if (num_ts > 0) then
    call interp_cv_to_flxpts
    call get_faceusp_bc
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
        else if (visc_flux_method == visc_flux_BR1) then
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
    end if
    !
    ! Output the error at the specified interval
    !
    if (mod(iter,results_interval)==0 .and. iter/=num_ts) then
      !
      end_wall_time = mpi_wtime()
      wall_time = end_wall_time-beg_wall_time
      work_units = wall_time/one_work_unit
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
        if (walltime_expiring) exit main_loop
        !
        ! Have the root processor check for the file
        !
        if (mypnum == 0) then
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
        end if
        !
        ! Have the root processor broadcast request_stop to
        ! the other processors so all processors can exit the
        ! main loop together if the STOP_RECON file was found
        !
        call mpi_bcast(request_stop,1,mpi_logical,0,mpi_comm_world,mpierr)
        if (request_stop) exit main_loop
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
        call write_nodal_parallel_cgns(0)
      end if
    end if
    !
    ! Continue the final countdown if the vortex has collapsed
    !
#ifndef PBS_ENV
    if (vortex_has_collapsed) then
      final_countdown = final_countdown + 1
      if (mod(final_countdown,25*Period_Timesteps) == 0) then
        if (mypnum == 0) write (iout,98)
        exit main_loop
      end if
    end if
#endif
    !
  end do main_loop
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
  end_wall_time = mpi_wtime()
  wall_time = end_wall_time-beg_wall_time
  work_units = wall_time/one_work_unit
  !
  !.. get min and max entropy value and location
  if (.not. ignore_errors) call getmaxent
  !
  call write_nodal_parallel_cgns(-1)
  !
  call write_restart_file(-1)
  !
  if (.not. ignore_errors) then
    !
    call compute_final_error
    !
   !if (ncpu == 1 .and. fals) call compute_bc_profiles
    !
    if (mypnum == 0) then
      !
      write (iout,40) itcur,time
      write (iout,41) wall_time,work_units
      !
      call write_final_error
      call write_results_file(-1,itcur)
      !
    end if
    !
  end if
  !
  write (iout,50)
  !
  ! Finalize the overall MPI session and stop execution
  !
  call finalize_mpi_environment
  !
  !
  !
  10 format (/,' For Processor ',i0,/,'------------------')
  11 format (/,' For the Global Grid ',/,'---------------------')
  12 format ('  ncell = ',i0,'  nnode = ',i0,'  nfbnd = ',i0,'  nDOF = ',i0)
  13 format ('  Partition Volume = ',f13.6)
  14 format ('  Total Grid Volume = ',f13.6)
  15 format ('  1D Solution Points = ',8f9.5)
  16 format (/," Allocated Wall Time = ",i0," s.",/)
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
