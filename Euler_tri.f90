!include "geovar_mod.f90"
!include "metvar_mod.f90"
!include "ovar_mod.f90"
!include "flowvar_mod.f90"
!include "mpivar_mod.f90"
!include "parallel.f90"
!include "poly.f90"
!include "generic.f90"
!include "io.f90"
!include "initialize.f90"
!include "startup.f90"
!include "triangle.f90"
!include "metrics.f90"
!include "flux.f90"
!
program Euler_tri
  !
  use module_kind_types
  use geovar
  use ovar
  use metvar
  use flowvar
 !use mpivar
  use io_mod, only : read_input_file, &
                    !read_input_file_test, &
                     read_grid_file, &
                     write_tecplot_file, &
                     write_restart_file, &
                     write_final_error, &
                     write_final_advec_error, &
                     write_final_vortex_error, &
                     write_entropy_error, &
                     write_cgns_solution, &
                     write_nodal_parallel_cgns, &
                     write_results_file, &
                     write_advec_results_file, &
                     write_vortex_results_file, &
                     write_tecplot_partitions
  use metrics_mod, only : get_solution_points_2d, &
                          get_grid_connectivity, &
                          get_metrics
  use startup_mod, only : get_quadrature_data
  use startup_mod, only : old_get_quadrature_data
  use initialization_mod, only : initialize_flow
  use module_bc,   only : get_faceusp_bc
  use module_bc,   only : get_facedusp_bc
  use module_bc,   only : compute_bc_profiles
  use generic_mod, only : rgetdtmax,get_dtsp, &
                          getbc,rresi, &
                          rgeterr,getmaxent, &
                          compute_final_error, &
                          compute_final_advec_error, &
                          compute_final_vortex_error, &
                          rupdate_Classic_RK, &
                          rupdate_CK4_RK, &
                          rupdate_TVD2_RK, &
                          rupdate_TVD3_RK, &
                          rupdate_TVD4_RK, &
                          test_error
  use flux_mod, only : flux_divergence_lp, &
                       flux_divergence_cr, &
                       flux_divergence_jac_cr, &
                       correct_flux_project_cv, &
                       correct_flux_project_pv, &
                       correct_flux_project_flx, &
                       dg_weak_form, &
                       dg_strong_form, &
                       interp_cv_to_flxpts, &
                       interp_dusp_to_flxpts, &
                       derflux1_qdtr, &
                       derflux1_qdtr_with_visc, &
                       correct_interface_flux_qdtr, &
                       correct_interface_flux_qdtr_with_visc, &
                       discontinuous_gradient_at_solpts, &
                       correct_gradients_at_solpts, &
                       correct_gradients_br2
  use parallel_mod, only : partition_grid, &
                           create_serial_cell_map, &
                           test_nonblocking_fexchange
  !
  use iso_c_binding, only : c_long
  !
  implicit none
  !
 !include "kind_defs.h"
 !include "constants.h"
  !
  integer  :: iter,i,ierr,k
  real(wp) :: dtmin
  !
  integer, parameter :: iostp = 2683
  !
 !integer, parameter :: LagrangePolynomial = 1
 !integer, parameter :: ChainRuleCV = 2
 !integer, parameter :: ChainRuleUtd = 3
 !integer, parameter :: DGWeakForm = 4
 !integer, parameter :: DGStrongForm = 5
 !integer, parameter :: OnPrimitiveVariables = 1
 !integer, parameter :: OnConservativeVariables = 2
 !integer, parameter :: OnFluxes = 3
  !
  logical :: request_stop
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
    max_wall_time = real(current_wtl,kind=vbp)
  end if
  !
  ! Broadcast the maximum wall time to all processors
  !
  call mpi_bcast(max_wall_time,1,mpi_flttyp,0,mpi_comm_world,mpierr)
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
    call partition_grid(ncpu,n_order,metis_option,bface,nodes_of_cell, &
                        nodes_of_cell_ptr,xyz_nodes)
  else
    call create_serial_cell_map(ncell)
  end if
  !
  ! Get the coordinates of the solution points within each grid cell
  !
  call get_solution_points_2d
  !
  ! Get the grid connectivity and metrics
  !
  call get_grid_connectivity
  !
  ! Compute the remaining quadrature information
  ! (Vandermonde matrix, stiffness matrix, LIFT / correction function, etc.)
  !
  call get_quadrature_data
 !call old_get_quadrature_data
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
    write (iout,15) (solpts_1d(k), k=1,n_sp1de)
  end if
  call mpi_barrier(mpi_comm_world,mpierr)
 !call stop_recon(stop_mpi,"main: before initialize_flow")
  !
  ! Initialize the flow variables for all solution points or
  ! read in the restart values if restarting a simulation
  !
  call initialize_flow
  if (any(flux_divergence == [TestQuadrature,TestViscous])) then
    call interp_cv_to_flxpts
    call get_faceusp_bc
  else
    call getbc
  end if
  !
  ! Output the initialized solution
  !
  if (cgns_is_parallel) then
    call write_nodal_parallel_cgns(1)
  else
    call write_cgns_solution(1)
  end if
  !
  work_units = zero
  if (any(itestcase == Advec_Problems)) then
    if (mypnum == 0) call write_advec_results_file(1,itcur)
  else if (any(itestcase == Vortex_Problems)) then
    if (mypnum == 0) call write_vortex_results_file(1,itcur)
  else
    if (mypnum == 0) call write_results_file(1,itcur)
  end if
  !
  ! If running the vortex or density problem,
  ! reset results_interval if it is set to zero
  !
 !if (any(itestcase == Transport_Problems)) then
 !  if (results_interval == 0) results_interval = max(Period_Timesteps/500,1)
 !  if (output_interval == 0) output_interval = Period_Timesteps
 !end if
  !
  if (mypnum == 0) then
    write (iout,19) p_order,n_order,q_order,n_sp1de,n_qp1de, &
                    n_sp2dt,n_qp2dt,n_sp2dq,n_qp2dq
  end if
  !
  ! Write the header for the terminal output 
  !
 !call test_nonblocking_fexchange
 !call stop_recon(stop_mpi,"main: before starting main iteration loop")
  if (mypnum == 0) write (iout,20)
  !
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%                                    %%%
  !%%%  Beginning of main iteration loop  %%%
  !%%%                                    %%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  main_loop: do iter = 1, num_ts
    !
    ! Current iteration count including restart iterations
    !
    itcur = itcur + 1
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
      ! Apply a correction to the flux derivatives by
      ! interpolating the interior fluxes to the interface
      !
      ! Compute the derivative of the discontinuous flux function
      !
      if (flux_divergence == LagrangePolynomial) then
        !
        call getbc
        !
        ! Use the derivative of the Lagrange polynomial
        !
        if (q_order > n_order) then
          call derflux1_qdtr
        else
          call flux_divergence_lp
        end if
        if (iop_correction == OnPrimitiveVariables) then
          call correct_flux_project_pv
        else if (iop_correction == OnConservativeVariables) then
          call correct_flux_project_cv
        else if (iop_correction == OnFluxes) then
          call correct_flux_project_flx
        end if
        !
      else if (flux_divergence == ChainRuleCV) then
        !
        call getbc
        !
        ! Use the chain rule on the conserved variables
        !
        call flux_divergence_cr
        if (iop_correction == OnPrimitiveVariables) then
          call correct_flux_project_pv
        else if (iop_correction == OnConservativeVariables) then
          call correct_flux_project_cv
        else if (iop_correction == OnFluxes) then
          call correct_flux_project_flx
        end if
        !
      else if (flux_divergence == ChainRuleUtd) then
        !
        call getbc
        !
        ! Use the chain rule on u*J
        !
        call flux_divergence_jac_cr
        if (iop_correction == OnPrimitiveVariables) then
          call correct_flux_project_pv
        else if (iop_correction == OnConservativeVariables) then
          call correct_flux_project_cv
        else if (iop_correction == OnFluxes) then
          call correct_flux_project_flx
        end if
        !
      else if (flux_divergence == DGWeakForm) then
        !
        call getbc
        !
        ! Use the weak form of the DG method
        !
        call dg_weak_form
        !
      else if (flux_divergence == DGStrongForm) then
        !
        call getbc
        !
        ! Use the weak form of the DG method
        !
        call dg_strong_form
        !
      else if (flux_divergence == TestQuadrature) then
        !
        call interp_cv_to_flxpts
        call get_faceusp_bc
        if (nst == 1) call get_dtsp(dtmin)
        call derflux1_qdtr
        call correct_interface_flux_qdtr
        !
      else if (flux_divergence == TestViscous) then
        !
        call interp_cv_to_flxpts
        call get_faceusp_bc
        if (nst == 1) call get_dtsp(dtmin)
        call discontinuous_gradient_at_solpts
        if (visc_flux_method == visc_flux_BR2) then
          call correct_gradients_br2
        else if (visc_flux_method == visc_flux_BR1) then
          call correct_gradients_at_solpts
          call interp_dusp_to_flxpts
        end if
        call get_facedusp_bc
        call derflux1_qdtr_with_visc
        call correct_interface_flux_qdtr_with_visc
       !call correct_flux_project_flx
        !
      end if
      !
      ! Compute the final residual for each solution point
      !
      call rresi
      !
      ! Store residual in rssprk(:,:,nst) and
      ! update the solutions in usp and vsp
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
    if (mod(iter,iter_out_interval) == 0) call rgeterr
    !
    ! Output the error at the specified interval
    !
    if (mod(iter,results_interval)==0 .and. iter/=num_ts) then
      !
      end_wall_time = mpi_wtime()
      wall_time = end_wall_time-beg_wall_time
      work_units = wall_time/one_work_unit
      !
      if (any(itestcase == Advec_Problems)) then
        call compute_final_advec_error
        if (mypnum == 0) call write_advec_results_file(0,itcur)
      else if (any(itestcase == Vortex_Problems)) then
        call compute_final_vortex_error
        if (mypnum == 0) call write_vortex_results_file(0,itcur)
      else
        call compute_final_error
        if (mypnum == 0) call write_results_file(0,itcur)
      end if
      !
    end if
    !
    ! Store the new solution as the backup for the next iteration
    !
    uoldsp(:,1:n_solpts) = usp(:,1:n_solpts)
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
        if (cgns_is_parallel) then
          call write_nodal_parallel_cgns(0)
        else
          call write_cgns_solution(0)
        end if
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
  if (all(flux_divergence /= [TestQuadrature,TestViscous])) then
    call getbc
  end if
  !
  !.. get min and max entropy value and location
  call getmaxent
  !
 !call write_tecplot_file(-1)
  if (cgns_is_parallel) then
    call write_nodal_parallel_cgns(-1)
  else
    call write_cgns_solution(-1)
  end if
  !
  call write_restart_file(-1)
  !
  if (any(itestcase == Advec_Problems)) then
    call compute_final_advec_error
  else if (any(itestcase == Vortex_Problems)) then
    call compute_final_vortex_error
  else
    call compute_final_error
  end if
  !
  if (ncpu == 1) call compute_bc_profiles
  if (mypnum == 0) then
    !
    write (iout,40) itcur,time
    write (iout,41) wall_time,work_units
    !
    if (any(itestcase == Advec_Problems)) then
      call write_final_advec_error
      call write_advec_results_file(-1,itcur)
    else if (any(itestcase == Vortex_Problems)) then
      call write_final_vortex_error
      call write_vortex_results_file(-1,itcur)
    else
      call write_final_error
      call write_results_file(-1,itcur)
    end if
    !
   !if (itestcase <= 1 .or. itestcase == 6) then
   !  call write_entropy_error
   !end if
    !
    write (iout,50)
    !
  end if
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
  19 format (/,' p_order = ',i4,/, &
               ' n_order = ',i4,/, &
               ' q_order = ',i4,/, &
               ' n_sp1de = ',i4,/, &
               ' n_qp1de = ',i4,/, &
               ' n_sp2dt = ',i4,/, &
               ' n_qp2dt = ',i4,/, &
               ' n_sp2dq = ',i4,/, &
               ' n_qp2dq = ',i4,/)
  20 format (/,5x,'iter',5x,'d_t',8x,'Time',6x,'L2-Entropy',2x,'L2-Density', &
               4x,'Cell',3x,'Solpt',6x,'x',10x,'y',6x,'Time Left')
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
