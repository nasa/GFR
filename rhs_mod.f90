module rhs_mod
  !
!#define PROFILE_ON
  !
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  use eqn_idx_mod, only : ncs,npres,ntemp,ngamm
  !
  use ovar, only : force_nd_subroutines
  use timer_mod, only : solver_timer
  !
  implicit none
  !
  private
  !
 !public :: compute_rhs
  public :: rresi
  public :: save_old_solution
  public :: stabilize_usp
  public :: update_faceusp
  public :: update_facegpsp
  !
contains
!
!###############################################################################
!
!subroutine compute_rhs(iter,nst,newton_failures)
!  !
!  ! NEED TO FIND WAY TO MOVE THIS HERE FROM main.f90
!  !
!end subroutine compute_rhs
!
!###############################################################################
!
subroutine save_old_solution()
  !
  !.. Use Statements ..
  use ovar,        only : gam
  use flowvar_mod, only : usp,uoldsp
  use flowvar_mod, only : gpsp,gpoldsp
  !
  !.. Local Scalars ..
  integer :: m,n,np,nr
  !
  !. Local Saved Scalars ..
  integer, save :: idx = -1
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "save_old_solution"
  !
continue
  !
  call solver_timer(entering_procedure,pname,idx)
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  np = size(uoldsp,dim=1)
  !
 !!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(uoldsp,usp)
 !do m = 1,nq
 !  do n = 1,np
 !    uoldsp(n,m) = usp(n,m)
 !  end do
 !end do
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(uoldsp,usp)
    do n = 1,np
      uoldsp(n,nec) = usp(n,nec)
      uoldsp(n,nmx) = usp(n,nmx)
      uoldsp(n,nmy) = usp(n,nmy)
      uoldsp(n,nmz) = usp(n,nmz)
      uoldsp(n,nee) = usp(n,nee)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(uoldsp,usp)
    do n = 1,np
      uoldsp(n,nec) = usp(n,nec)
      uoldsp(n,nmx) = usp(n,nmx)
      uoldsp(n,nmy) = usp(n,nmy)
      uoldsp(n,nee) = usp(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(uoldsp,usp)
    do m = 1,nq
      do n = 1,np
        uoldsp(n,m) = usp(n,m)
      end do
    end do
  end if
  !
  !$ACC PARALLEL LOOP PRESENT(gpoldsp,gpsp)
  do n = 1,np
    gpoldsp(n,npres) = gpsp(n,npres)
    gpoldsp(n,ntemp) = gpsp(n,ntemp)
  end do
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP PRESENT(gpoldsp,gpsp)
    do n = 1,np
      gpoldsp(n,ngamm) = gpsp(n,ngamm)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(gpoldsp,gpsp)
    do n = 1,np
      gpoldsp(n,ngamm) = gam
    end do
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  call solver_timer(leaving_procedure,pname,idx)
  !
end subroutine save_old_solution
!
!###############################################################################
!
subroutine rresi(iter,nst)
  !
  !... calculates d(utd)/dt
  !
  !.. Use Statements ..
  use order_mod,   only : p_order,n_order
  use geomvar_mod, only : n_solpts
  use ovar,        only : itestcase
  use metrics_mod, only : jac
  use flowvar_mod, only : residual
  !
  use projection_mod, only : project_to_porder
  !
  !.. Formal Arguments ..
  integer, intent(in) :: iter
  integer, intent(in) :: nst
  !
  !.. Local Scalars ..
  integer  :: io,m,n,nr
  real(wp) :: jac_inv
  !
  !. Local Saved Scalars ..
  integer, save :: idx = -1
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "rresi"
  !
continue
  !
  call solver_timer(entering_procedure,pname,idx)
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(jac,residual)
    do n = 1,n_solpts
      jac_inv = one/jac(n)
      residual(n,nec) = -residual(n,nec)*jac_inv
      residual(n,nmx) = -residual(n,nmx)*jac_inv
      residual(n,nmy) = -residual(n,nmy)*jac_inv
      residual(n,nmz) = -residual(n,nmz)*jac_inv
      residual(n,nee) = -residual(n,nee)*jac_inv
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(jac,residual)
    do n = 1,n_solpts
      jac_inv = one/jac(n)
      residual(n,nec) = -residual(n,nec)*jac_inv
      residual(n,nmx) = -residual(n,nmx)*jac_inv
      residual(n,nmy) = -residual(n,nmy)*jac_inv
      residual(n,nee) = -residual(n,nee)*jac_inv
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(jac,residual)
    do m = 1,nq
      do n = 1,n_solpts
        residual(n,m) = -residual(n,m)/jac(n)
      end do
    end do
  end if
  !
 !io = 100 + 10*iter + nst
 !!$ACC UPDATE SELF(residual)
 !do n = 1,n_solpts
 !  write (io,100) (residual(n,m), m=1,nq)
 !end do
 !100 format (5es16.8)
  !
  ! Zero out the residuals for all but the density
  ! equation if running the density transport test case
  !
  if (itestcase == Density_Transport) then
    if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
      !$ACC PARALLEL LOOP PRESENT(residual)
      do n = 1,n_solpts
        residual(n,nmx) = residual(n,nec)
        residual(n,nmy) = residual(n,nec)
        residual(n,nmz) = residual(n,nec)
        residual(n,nee) = residual(n,nec)
      end do
    else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
      !$ACC PARALLEL LOOP PRESENT(residual)
      do n = 1,n_solpts
        residual(n,nmx) = residual(n,nec)
        residual(n,nmy) = residual(n,nec)
        residual(n,nee) = residual(n,nec)
      end do
    else
      !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(residual)
      do m = 1,nec-1
        do n = 1,n_solpts
          residual(n,m) = residual(n,nec)
        end do
      end do
      !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(residual)
      do m = nec+1,nq
        do n = 1,n_solpts
          residual(n,m) = residual(n,nec)
        end do
      end do
    end if
  end if
  !
#ifdef DEBUG_ON
!#define DEBUG_MRF_FREESTREAM_TCF
#ifdef DEBUG_MRF_FREESTREAM_TCF
block
  use geomvar_mod,only : xyz
  use flowvar_mod, only : dtsp
  real(wp) :: x,y,z,r,t
  if (abs(itestcase) == Freestream_TCF) then
    write (371,1) iter,nst
    do n = 1,n_solpts
      if (any(abs(dtsp(n)*residual(n,:)) > rndoff)) then
        x = xyz(1,n)
        y = xyz(2,n)
        z = xyz(3,n)
        r = sqrt(y*y+z*z)
        t = atan2(z,y)
        t = (t + merge(two*pi,zero,t<zero))*radians2degrees
        write (371,2) n,x,r,t,(dtsp(n)*residual(n,m), m=1,nq)
      end if
    end do
  end if
  1 format (/,2x,"Iteration = ",i4,5x,"RK stage = ",i1,/,80("-"))
  2 format (i8,3x,3(2x,f12.6),3x,5es20.10)
end block
#endif
#endif
  !
  if (p_order < n_order) then
    call project_to_porder(residual)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  call solver_timer(leaving_procedure,pname,idx)
  !
end subroutine rresi
!
!###############################################################################
!
subroutine check_total_temperature(np,usp,gpsp,nst)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd,tref,tmin,tmax
  use ovar, only : ideal_total_temperature
  use ovar, only : limit_total_temperature
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  integer,     intent(in) :: nst
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: gpsp(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr,loc_ttmin,loc_ttmax
  integer  :: num_min_limits
  integer  :: num_max_limits
  real(wp) :: gm1,rgas_inv,gm1_over_2gmR
  real(wp) :: rho,rho_inv,pres,temp,total_temp
  real(wp) :: velx,vely,velz,q2,a2
  real(wp) :: tmin_lim,tmax_lim,ttmin,ttmax
  !
  !.. Local Arrays ..
  integer :: num_temp_limits(1:2)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_total_temperature"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (ncs > 0) return
  !
  gm1 = gam - one
  rgas_inv = one / rgasref_nd
  gm1_over_2gmR = half*gm1*rgas_inv/gam
  !
  num_min_limits = 0
  num_max_limits = 0
  !
  ttmin = tmax/tref ! initialize ttmin to the largest  possible temperature
  ttmax = tmin/tref ! initialize ttmax to the smallest possible temperature
  !
  tmin_lim = ideal_total_temperature*(one - limit_total_temperature)
  tmax_lim = ideal_total_temperature*(one + limit_total_temperature)
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    !
    !$ACC PARALLEL LOOP PRESENT(usp,gpsp) &
    !$ACC&              REDUCTION(+:num_min_limits,num_max_limits) &
    !$ACC&              REDUCTION(max:ttmax) &
    !$ACC&              REDUCTION(min:ttmin)
    do n = 1,np
      !
      rho = usp(n,nec)
      !
      rho_inv = one/rho
      !
      velx = usp(n,nmx)*rho_inv
      vely = usp(n,nmy)*rho_inv
      velz = usp(n,nmz)*rho_inv
      !
      q2 = velx*velx + vely*vely + velz*velz
      !
      pres = gpsp(n,npres)
      temp = gpsp(n,ntemp)
      !
      a2 = max( gam*pres*rho_inv , rndoff )
      !
      total_temp = temp*(one + half*gm1*q2/a2)
      !
      if (total_temp > ttmax) then
        ttmax = total_temp
        loc_ttmax = n
      else if (total_temp < ttmin) then
        ttmin = total_temp
        loc_ttmin = n
      end if
      !
      if (total_temp > tmax_lim) then
        total_temp = tmax_lim
        num_max_limits = num_max_limits + 1
      else if (total_temp < tmin_lim) then
        total_temp = tmin_lim
        num_min_limits = num_min_limits + 1
      end if
      !
      temp = total_temp - q2*gm1_over_2gmR
      gpsp(n,ntemp) = temp
      !
!#define APPLY_TTOTAL_LIMIT_TO_PRESSURE
#ifdef APPLY_TTOTAL_LIMIT_TO_PRESSURE
      pres = rho*rgasref_nd*temp
      usp(n,nee) = pres/gm1 + half*rho*q2
      gpsp(n,npres) = pres
#else
      usp(n,nec) = rgas_inv*pres/temp
#endif
      !
    end do
    !
  else if (nr == 2) then
    !
    !$ACC PARALLEL LOOP PRESENT(usp,gpsp) &
    !$ACC&              REDUCTION(+:num_min_limits,num_max_limits) &
    !$ACC&              REDUCTION(max:ttmax) &
    !$ACC&              REDUCTION(min:ttmin)
    do n = 1,np
      !
      rho = usp(n,nec)
      !
      rho_inv = one/rho
      !
      velx = usp(n,nmx)*rho_inv
      vely = usp(n,nmy)*rho_inv
      !
      q2 = velx*velx + vely*vely
      !
      pres = gpsp(n,npres)
      temp = gpsp(n,ntemp)
      !
      a2 = max( gam*pres*rho_inv , rndoff )
      !
      total_temp = temp*(one + half*gm1*q2/a2)
      !
      if (total_temp > ttmax) then
        ttmax = total_temp
        loc_ttmax = n
      else if (total_temp < ttmin) then
        ttmin = total_temp
        loc_ttmin = n
      end if
      !
      if (total_temp > tmax_lim) then
        total_temp = tmax_lim
        num_max_limits = num_max_limits + 1
      else if (total_temp < tmin_lim) then
        total_temp = tmin_lim
        num_min_limits = num_min_limits + 1
      end if
      !
      temp = total_temp - q2*gm1_over_2gmR
      gpsp(n,ntemp) = temp
      !
#ifdef APPLY_TTOTAL_LIMIT_TO_PRESSURE
      pres = rho*rgasref_nd*temp
      usp(n,nee) = pres/gm1 + half*rho*q2
      gpsp(n,npres) = pres
#else
      usp(n,nec) = rgas_inv*pres/temp
#endif
      !
    end do
    !
  else
    !
    !$ACC PARALLEL LOOP PRESENT(usp,gpsp) &
    !$ACC&              REDUCTION(+:num_min_limits,num_max_limits) &
    !$ACC&              REDUCTION(max:ttmax) &
    !$ACC&              REDUCTION(min:ttmin)
    do n = 1,np
      !
      rho = usp(n,nec)
      !
      rho_inv = one/rho
      !
      velx = usp(n,nmx)*rho_inv
      !
      q2 = velx*velx
      !
      pres = gpsp(n,npres)
      temp = gpsp(n,ntemp)
      !
      a2 = max( gam*pres*rho_inv , rndoff )
      !
      total_temp = temp*(one + half*gm1*q2/a2)
      !
      if (total_temp > ttmax) then
        ttmax = total_temp
        loc_ttmax = n
      else if (total_temp < ttmin) then
        ttmin = total_temp
        loc_ttmin = n
      end if
      !
      if (total_temp > tmax_lim) then
        total_temp = tmax_lim
        num_max_limits = num_max_limits + 1
      else if (total_temp < tmin_lim) then
        total_temp = tmin_lim
        num_min_limits = num_min_limits + 1
      end if
      !
      temp = total_temp - q2*gm1_over_2gmR
      gpsp(n,ntemp) = temp
      !
#ifdef APPLY_TTOTAL_LIMIT_TO_PRESSURE
      pres = rho*rgasref_nd*temp
      usp(n,nee) = pres/gm1 + half*rho*q2
      gpsp(n,npres) = pres
#else
      usp(n,nec) = rgas_inv*pres/temp
#endif
      !
    end do
    !
  end if
  !
!#define OUTPUT_TOTAL_TEMP_LIMITS
#ifdef OUTPUT_TOTAL_TEMP_LIMITS
  num_temp_limits(1:2) = [num_min_limits,num_max_limits]
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,num_temp_limits, &
                       size(num_temp_limits,kind=int_mpi),mpi_inttyp, &
                       MPI_SUM,MPI_COMM_WORLD,mpierr)
    call mpi_allreduce(MPI_IN_PLACE,ttmin,1_int_mpi,mpi_flttyp, &
                       MPI_MIN,MPI_COMM_WORLD,mpierr)
    call mpi_allreduce(MPI_IN_PLACE,ttmax,1_int_mpi,mpi_flttyp, &
                       MPI_MAX,MPI_COMM_WORLD,mpierr)
  end if
  if (mypnum == glb_root) then
    if (any(num_temp_limits > 0)) then
      write (iout,100) nst,num_temp_limits(1),ttmin, &
                           num_temp_limits(2),ttmax
      flush (iout)
    end if
  end if
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  100 format (6x,"RK Stage #",i0,t20,":  # min T0 limits = ",i0,",", &
                                 t55,"min T0 = ",f8.2,/, &
                                 t23,"# max T0 limits = ",i0,",", &
                                 t55,"max T0 = ",f8.2)
  !
end subroutine check_total_temperature
!
!###############################################################################
!
subroutine stabilize_usp(nst,Apply_Filter,Apply_Sponge,Apply_Limiter, &
                         caller_name)
  !
  !.. Use Statements ..
  use geomvar_mod,         only : n_solpts,n_flxpts
  use flowvar_mod,         only : usp,gpsp,faceusp,facegpsp
  use ovar,                only : Limiter_Option,Filter_Option
  use ovar,                only : Filter_Alpha
  use ovar,                only : ideal_total_temperature
  use ovar,                only : Sponge_Filter_Option
  use ovar,                only : thermofix
  use ovar,                only : governing_equations
  use ovar,                only : positivity, entropy_limiter
  use filter_mod,          only : filter_usp
  use limiters_mod,        only : resolution_indicator
  use limiters_mod,        only : limit_by_projecting
  use limiters_mod,        only : apply_sponge_filter
  use entropy_limiter_mod, only : perform_entropy_limiting
  use entropy_limiter_mod, only : update_parallel_entropy_bounds
  use mlp_mod,             only : limit_by_mlp
  use flux_mod,            only : interp_cv_to_flxpts
  use flux_mod,            only : discontinuous_gradient_at_solpts
  !
  use thermofix_mod,       only : check_thermo_properties
  use thermofix_mod,       only : check_cell_thermo_properties
  use thermo_mod,          only : update_mixture_density
  use thermo_mod,          only : update_thermo_properties
  use transport_mod,       only : update_transport_properties
  !
  use positivity_mod,      only : positivity_preserving_limiter
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: nst           !< Current Runge-Kutta stage.
  logical(lk),      intent(in) :: Apply_Filter  !< Apply filter  this stage.
  logical(lk),      intent(in) :: Apply_Sponge  !< Apply sponge  this stage.
  logical(lk),      intent(in) :: Apply_Limiter !< Apply limiter this stage.
  character(len=*), intent(in) :: caller_name   !< Name of calling routine.
  !
  !. Local Saved Scalars ..
  integer, save :: idx = -1
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "stabilize_usp"
  !
continue
  !
  call solver_timer(entering_procedure,pname,idx)
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Update the mixture density if solving a gas mixture
  ! NOTE: This will also be called first thing upon entering subroutine
  !       update_thermo_properties. However, we need to make sure this
  !       is done before calling stabilize_usp. It is a relatively simple
  !       operation so the double work shouldnt amount to much additional
  !       computational cost.
  !
  call update_mixture_density(1,n_solpts,usp)
  call check_for_bad_number(__LINE__,__FILE__)
  !
  ! Apply the sponge filter to cells in sponge regions
  !
  if (Sponge_Filter_Option == 1) then ! The sponge filter is active.
    if (Apply_Sponge) then  ! Apply sponge filter during this stage.
#ifdef PROFILE_ON
      if (mypnum == glb_root) write (iout,1) nst
      1 format (8x,"Applying Sponge Filter : RK Stage = ",i0)
#endif
      call apply_sponge_filter(usp)
      call check_for_bad_number(__LINE__,__FILE__)
    end if
  end if
  !
  ! Apply the global/background exponential filter. This filter is applied
  ! to every cell regardless of whether or not the solution in the cell
  ! is stable/smooth/well-resolved.
  ! NOTE: Apply this global filter before the local limiter filter only if
  !       Filter_Alpha <= 1.0. This corresponds to a very weak filter that
  !       will damp the highest mode by at most 64 percent (Alpha=1.0), which
  !       should still allow the resolution marker to find troubled cells.
  !
  if (Filter_Option /= Filter_Option__Off .and. Filter_Alpha <= one) then
    ! A weak global filter is active
    if (Apply_Filter) then ! Apply the global filter during this stage
#ifdef PROFILE_ON
      if (mypnum == glb_root) write (iout,1) nst
#endif
      call filter_usp
      call check_for_bad_number(__LINE__,__FILE__)
    end if
  end if
  !
  ! Ensures the positivity of density and pressure
  if (Positivity%Option /= POSLIM_OPT__OFF) then
    call positivity_preserving_limiter(nst)
  end if
  !
  if (entropy_limiter%Option /= ENTLIM_OPT__OFF) then
    call update_parallel_entropy_bounds()
    call perform_entropy_limiting(nst)
  end if
  !
  ! Apply realizable total temperature limiter
  !
  if (ideal_total_temperature > zero) then
    call update_thermo_properties(1,n_solpts,usp,gpsp, &
                                  loc="cell",caller_name=caller_name)
    call check_total_temperature(n_solpts,usp,gpsp,nst)
    !
    call update_thermo_properties(1,n_flxpts,faceusp,facegpsp,"stabilize_usp")
  end if
  !
  ! Apply the shock/troubled-cell limiter. The cells will first be checked
  ! to find any with signs of instability or shocks, and these cells will
  ! be marked as 'troubled'. The limiting method (filtering, projection, etc.)
  ! will then be applied to only the troubled cells and all other cells
  ! will remain unchanged.
  !
  if (Limiter_Option /= Limiter_Option__Off) then  ! The limiter is active.
    if (Apply_Limiter) then  ! Apply limiter during this stage.
#ifdef PROFILE_ON
      if (mypnum == glb_root) write (iout,2) nst
      2 format (8x,"Applying Limiter : RK Stage = ",i0)
#endif
      select case (Limiter_Option)
        case (Limiter_Option__ProjectionScaled)
          call resolution_indicator(nst=nst,apply_projection=true)
        case (Limiter_Option__Filter, &
              Limiter_Option__FilterScaled, &
              Limiter_Option__ProjectionRecursive)
          call limit_by_projecting(nst=nst)
        case (Limiter_Option__MLP)
          call limit_by_mlp(nst=nst,relax=one)
        case default
          ! Should not get here.
      end select
      call check_for_bad_number(__LINE__,__FILE__)
    end if
  end if
  !
  ! Apply the global/background exponential filter. This filter is applied
  ! to every cell regardless of whether or not the solution in the cell
  ! is stable/smooth/well-resolved.
  ! NOTE: Apply the global filter after the local limiter filter if
  !       Filter_Alpha > 1, because this global filter will damp the
  !       highest mode too much preventing the resolution marker from
  !       finding any troubled cells.
  !
  if (Filter_Option /= Filter_Option__Off .and. Filter_Alpha > one) then
    ! A strong global filter is active
    if (Apply_Filter) then ! Apply the global filter during this stage
#ifdef PROFILE_ON
      if (mypnum == glb_root) write (iout,3) nst
      3 format (8x,"Applying Standard Filter : RK Stage = ",i0)
#endif
      call filter_usp
      call check_for_bad_number(__LINE__,__FILE__)
    end if
  end if
  !
  ! Update the thermodynamic and transport properties at the solution points
  !
  call update_thermo_properties(1,n_solpts,usp,gpsp, &
                                loc="cell",caller_name=caller_name)
  call check_for_bad_number(__LINE__,__FILE__)
  !
  ! Check for any negative thermodynamic properties at the solution points
  !
 !call check_cell_thermo_properties()  ! obsolete
  !
  ! Update the thermodynamic and transport properties at the solution points
  !
  call update_transport_properties(1,n_solpts,usp,gpsp)
  call check_for_bad_number(__LINE__,__FILE__)
  !
  ! Interpolate the new updated solution to the face points
  !
  call interp_cv_to_flxpts
  call check_for_bad_number(__LINE__,__FILE__)
  !
  ! Check the cell and face thermo properties
  !
  call check_thermo_properties( thermofix%check_cells, thermofix%check_faces, &
                                thermofix%check_nodes, update_faces=true, &
                                verbose=thermofix%verbose )
  call check_for_bad_number(__LINE__,__FILE__)
  !
  ! Compute the discontinuous gradient at the solution points
  !
  if (governing_equations == NavierStokes_Eqns) then
    call discontinuous_gradient_at_solpts
    call check_for_bad_number(__LINE__,__FILE__)
  end if
  !
  !
  call debug_timer(leaving_procedure,pname)
  !
  call solver_timer(leaving_procedure,pname,idx)
  !
end subroutine stabilize_usp
!
!###############################################################################
!
subroutine update_gp_face_dusp_arrays(caller_name,skip_thermo_check)
  !
  !.. Use Statements ..
  use geomvar_mod,     only : n_solpts
  use flowvar_mod,     only : usp,gpsp
  use ovar,            only : thermofix
  use ovar,            only : governing_equations
  use flux_mod,        only : interp_cv_to_flxpts
  use flux_mod,        only : discontinuous_gradient_at_solpts
  !
  use thermofix_mod,   only : check_thermo_properties
  use thermo_mod,      only : update_mixture_density
  use thermo_mod,      only : update_thermo_properties
  use transport_mod,   only : update_transport_properties
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: caller_name
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: skip_thermo_check
  !
  !.. Local Scalars ..
  logical(lk) :: perform_thermo_check
  !
  !. Local Saved Scalars ..
  integer, save :: idx = -1
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_gp_face_dusp_arrays"
  !
continue
  !
  call solver_timer(entering_procedure,pname,idx)
  !
  call debug_timer(entering_procedure,pname)
  !
  perform_thermo_check = true
  if (present(skip_thermo_check)) then
    perform_thermo_check = (.not. skip_thermo_check)
  end if
  !
  ! Update the thermodynamic and transport properties at the solution points
  !
  call update_thermo_properties(1,n_solpts,usp,gpsp, &
                                loc="cell",caller_name=caller_name)
  !
  ! Update the thermodynamic and transport properties at the solution points
  !
  call update_transport_properties(1,n_solpts,usp,gpsp)
  !
  ! Interpolate the new updated solution to the face points
  !
  call interp_cv_to_flxpts
  !
  ! Check for any negative thermodynamic properties
  ! at the solution points or flux points
  !
  if (perform_thermo_check) then
    call check_thermo_properties(thermofix%check_cells, thermofix%check_faces, &
                                 thermofix%check_nodes, update_faces=true, &
                                 verbose=thermofix%verbose )
  end if
  !
  ! Compute the discontinuous gradient at the solution points
  !
  if (governing_equations == NavierStokes_Eqns) then
    call discontinuous_gradient_at_solpts
  end if
  !
  !
  call debug_timer(leaving_procedure,pname)
  !
  call solver_timer(leaving_procedure,pname,idx)
  !
end subroutine update_gp_face_dusp_arrays
!
!###############################################################################
!
subroutine update_faceusp(iter,nst,newton_failures)
  !
  !.. Use Statements ..
  use bc_mod, only : get_faceusp_bc
  !
  use parallel_mod, only : exchange_faceusp
  use parallel_mod, only : wait_faceusp
  !
  use flux_mod, only : interp_cv_to_flxpts
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
  character(len=*), parameter :: pname = "update_faceusp"
  !
continue
  !
  call solver_timer(entering_procedure,pname,idx)
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Interpolate the conserved variables from
  ! the solution points to the face/flux points
  ! NOTE: Skip this if this is the first Runge-Kutta stage, since this is
  !       done immediately after the last RK stage.
  !
  call interp_cv_to_flxpts(iter,nst)
                                    !  call check_for_bad_number(__LINE__)
  !
  ! Initiate the exchange of face/flux point solutions
  ! between adjacent partitions
  !
  if (ncpu > 1) then
    call exchange_faceusp           !; call check_for_bad_number(__LINE__)
#ifdef IMMEDIATE_WAIT_FOR_MPI
    call wait_faceusp               !; call check_for_bad_number(__LINE__)
#endif
  end if
  !
  ! Get the boundary conditions for the boundary face/flux points
  !
  call get_faceusp_bc(newton_failures)
 !call check_for_bad_number(__LINE__,newton_failures)
  !
  ! Ensure that the exchange of face/flux point solutions between adjacent
  ! partitions has been completed, because that information will be needed
  ! in the next step.
  !
  if (ncpu > 1) then
#ifndef IMMEDIATE_WAIT_FOR_MPI
    call wait_faceusp               !; call check_for_bad_number(__LINE__)
#endif
  end if
  !
  ! Finally, update facegpsp using the updated faceusp
  !
  call update_facegpsp
  !
  call debug_timer(leaving_procedure,pname)
  !
  call solver_timer(leaving_procedure,pname,idx)
  !
end subroutine update_faceusp
!
!###############################################################################
!
subroutine update_facegpsp()
  !
  !.. Use Statements ..
  use geomvar_mod,   only : n_bnd_flxpts
  use flowvar_mod,   only : faceusp,facegpsp
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Local Scalars ..
  integer :: n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_facegpsp"
  !
continue
  !
  ! Update the thermo and transport properties for all the
  ! exterior flux points of all the boundary faces.
  !
  n1 = 1
  n2 = n_bnd_flxpts
  !
  call update_thermo_properties(n1,n2,faceusp(:,:,2),facegpsp(:,:,2), &
                                      loc="face",caller_name=pname)
  call update_transport_properties(n1,n2,faceusp(:,:,2),facegpsp(:,:,2))
  !
end subroutine update_facegpsp
!
!###############################################################################
!
#include "Functions/debug_timer.f90"
! The following FPP definition disables the body of check_for_bad_number
#define DISABLE_CHECK_FOR_BAD_NUMBER
#include "Functions/check_for_bad_number.f90"
!
!###############################################################################
!
end module rhs_mod
