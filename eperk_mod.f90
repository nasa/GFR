module eperk_mod
  !
  use module_kind_types
  !
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq,ncs
  !
  use ovar, only : force_nd_subroutines
  use timer_mod, only : solver_timer
  !
  implicit none
  !
  private
  !
  public :: update_EPERK
  public :: get_EPERK_adaptive_dtsp
  public :: init_eperk_coefficients
  !
  !
  integer, save :: this_iter
  !
  logical(lk), save :: dtmin_was_modified = fals
  real(wp),    save :: previous_dtmin = zero
  !
  real(wp), save, allocatable :: error_dt(:,:)
  real(wp), save, allocatable :: my_error(:)
  real(wp), save, allocatable :: all_error(:,:)
  !
  real(wp), save, allocatable :: uepsp(:,:)
  real(wp), save, allocatable :: rhs_save(:,:)
  !
  real(wp), save, allocatable :: oeperk_coef_c(:)
  real(wp), save, allocatable :: oeperk_coef_b(:)
  real(wp), save, allocatable :: oeperk_coef_bhat(:)
  real(wp), save, allocatable :: oeperk_coef_g1(:)
  real(wp), save, allocatable :: oeperk_coef_g2(:)
  real(wp), save, allocatable :: oeperk_coef_g3(:)
  real(wp), save, allocatable :: oeperk_coef_d(:)
  !
  real(wp), save :: pseudo_time
  !
  real(wp), save :: order_exp
  !
  real(wp), save :: err_dt_atol
  real(wp), save :: err_dt_rtol
  !
  real(wp), save :: safety_factor
  real(wp), save :: safety_min
  real(wp), save :: safety_max
  !
  integer, save :: n_controllers
  !
  real(wp), save, allocatable :: pid_beta(:)
  !
  logical(lk), save :: EPERK_IS_FSAL = fals
  !
contains
!
!###############################################################################
!
subroutine update_EPERK(nst,iter,Apply_Filter,Apply_Sponge,Apply_Limiter)!, &
                       !dtmin,final_timestep_from_stop_request)
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_solpts
  use flowvar_mod, only : residual,usp,uoldsp,rssprk,dtsp,rk_sum
  use ovar,        only : time,num_rk_stages
  use rhs_mod,     only : stabilize_usp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  integer,     intent(in) :: iter
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Sponge
  logical(lk), intent(in) :: Apply_Limiter
 !real(wp),   intent(out) :: dtmin
 !integer,     intent(in) :: final_timestep_from_stop_request
  !
  !. Local Saved Scalars ..
  integer, save :: idx = -1
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK"
  !
continue
  !
  call solver_timer(entering_procedure,pname,idx)
  !
  call debug_timer(entering_procedure,pname)
  !
  this_iter = iter
  !
  if (iter == 1) then
    if (nst == 1) then
      if (allocated(rssprk)) then
        rssprk(:,:,:) = zero
      else
        if (.not.allocated(uepsp)) then
          uepsp = init_array(n_solpts,nq,zero) ! F2003 AUTO-REALLOC
        end if
      end if
      if (.not. allocated(rhs_save)) then
        rhs_save = init_array(n_solpts,nq,zero) ! F2003 AUTO-REALLOC
      end if
    end if
  end if
  !
  ! Update the solution using the current Runge-Kutta stage
  !
  if (num_rk_stages == 12) then
    !
    if (nst == 1) then
      pseudo_time = time
      call update_EPERK_DP853_S1(n_solpts,dtsp,residual,uoldsp, &
                                 usp,rssprk,pseudo_time,rhs_save)
    else if (nst == 2) then
      call update_EPERK_DP853_S2(n_solpts,dtsp,residual,uoldsp, &
                                 usp,rssprk,pseudo_time)
    else if (nst == 3) then
      call update_EPERK_DP853_S3(n_solpts,dtsp,residual,uoldsp, &
                                 usp,rssprk,pseudo_time)
    else if (nst == 4) then
      call update_EPERK_DP853_S4(n_solpts,dtsp,residual,uoldsp, &
                                 usp,rssprk,pseudo_time)
    else if (nst == 5) then
      call update_EPERK_DP853_S5(n_solpts,dtsp,residual,uoldsp, &
                                 usp,rssprk,pseudo_time)
    else if (nst == 6) then
      call update_EPERK_DP853_S6(n_solpts,dtsp,residual,uoldsp, &
                                 usp,rssprk,pseudo_time)
    else if (nst == 7) then
      call update_EPERK_DP853_S7(n_solpts,dtsp,residual,uoldsp, &
                                 usp,rssprk,pseudo_time)
    else if (nst == 8) then
      call update_EPERK_DP853_S8(n_solpts,dtsp,residual,uoldsp, &
                                 usp,rssprk,pseudo_time)
    else if (nst == 9) then
      call update_EPERK_DP853_S9(n_solpts,dtsp,residual,uoldsp, &
                                 usp,rssprk,pseudo_time)
    else if (nst == 10) then
      call update_EPERK_DP853_S10(n_solpts,dtsp,residual,uoldsp, &
                                  usp,rssprk,pseudo_time)
    else if (nst == 11) then
      call update_EPERK_DP853_S11(n_solpts,dtsp,residual,uoldsp, &
                                  usp,rssprk,pseudo_time)
    else if (nst == 12) then
      call update_EPERK_DP853_S12(n_solpts,dtsp,residual,uoldsp,usp,rssprk)
    end if
    !
  else if (num_rk_stages == 6) then
    !
    if (nst == 1) then
      pseudo_time = time
      call update_EPERK_DP54_S1(n_solpts,dtsp,residual,uoldsp, &
                                usp,rssprk,pseudo_time,rhs_save)
    else if (nst == 2) then
      call update_EPERK_DP54_S2(n_solpts,dtsp,residual,uoldsp, &
                                usp,rssprk,pseudo_time)
    else if (nst == 3) then
      call update_EPERK_DP54_S3(n_solpts,dtsp,residual,uoldsp, &
                                usp,rssprk,pseudo_time)
    else if (nst == 4) then
      call update_EPERK_DP54_S4(n_solpts,dtsp,residual,uoldsp, &
                                usp,rssprk,pseudo_time)
    else if (nst == 5) then
      call update_EPERK_DP54_S5(n_solpts,dtsp,residual,uoldsp, &
                                usp,rssprk,pseudo_time)
    else if (nst == 6) then
      call update_EPERK_DP54_S6(n_solpts,dtsp,residual,uoldsp,usp,rssprk)
    end if
    !
  else
    !
    call update_Optimized_EPERK(nst,n_solpts,dtsp,residual,uoldsp,usp, &
                                rk_sum,uepsp,rhs_save)
    !
  end if
  !
  ! Apply filtering or limiting to stabilize the solution, and then update
  ! the thermodynamic and transport properties at the solution points
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Sponge,Apply_Limiter, &
                     caller_name=pname)
  !
  call debug_timer(leaving_procedure,pname)
  !
  call solver_timer(leaving_procedure,pname,idx)
  !
end subroutine update_EPERK
!
!###############################################################################
!
subroutine get_EPERK_adaptive_dtsp(iter,dtmin,next_dtmin, &
                                   final_timestep_from_stop_request, &
                                   previous_timestep_rejected)
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_solpts
  use flowvar_mod, only : residual,usp,uoldsp,rssprk,dtsp
  use ovar,        only : d_t,time,num_rk_stages
  use dt_mod,      only : get_dtsp
  use dt_mod,      only : checks_on_dtmin
  !
  !.. Formal Arguments ..
  integer,      intent(in) :: iter
  real(wp),     intent(in) :: dtmin
  real(wp),    intent(out) :: next_dtmin
  integer,      intent(in) :: final_timestep_from_stop_request
  logical(lk), intent(out) :: previous_timestep_rejected
  !
  !.. Local Scalars ..
  integer  :: n,m
  real(wp) :: dtmin_save
  !
  !. Local Saved Scalars ..
  integer, save :: idx = -1
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_EPERK_adaptive_dtsp"
  !
 !logical(lk), parameter :: use_CFL_dtsp = true
  logical(lk), parameter :: use_CFL_dtsp = fals
  !
continue
  !
  call solver_timer(entering_procedure,pname,idx)
  !
  call debug_timer(entering_procedure,pname)
  !
  if (num_rk_stages /= 12) then
    if (.not.allocated(uepsp)) then
      uepsp = init_array(n_solpts,nq,zero) ! F2003 AUTO-REALLOC
    end if
  end if
  !
  previous_timestep_rejected = fals
  !
  next_dtmin = dtmin
  !
  if (use_CFL_dtsp) then
    !
    call get_dtsp(next_dtmin,final_timestep_from_stop_request)
    !
  else
    !
    ! If dtmin was modified by checks_on_dtmin last time step, we need to
    ! reset it to the value before it was modified so the adaptive controller
    ! doesnt use the modified value when computing the next time step size.
    !
    if (dtmin_was_modified) then
      next_dtmin = eps1*previous_dtmin
      dtmin_was_modified = fals
      error_dt(:,:) = one
    end if
    !
    ! Shift the time step error values before computing
    ! the time step error from the previous time step.
    ! NOTE: The use of eoshift here is equivalent to the following loops
    !           do n = n_controllers,2,-1
    !             do m = 1,nq
    !               error_dt(m,n) = error_dt(m,n-1)
    !             end do
    !           end do
    !           error_dt(:,1) = zero
    !
    error_dt(:,:) = eoshift(error_dt(:,:),shift=to_end,dim=2)
    !
    if (num_rk_stages == 12) then
      !
      call update_EPERK_DP853_S13(n_solpts,rssprk)
      !
      call get_embedded_pair_time_error_dop853(n_solpts,uoldsp,usp, &
                                               rssprk,error_dt(1:nq,1))
      !
    else if (num_rk_stages == 6) then
      !
      ! Compute the embedded pair solution using the RHS computed with
      ! the solution from the last time step.
      !
      call update_EPERK_DP54_S7(n_solpts,dtsp,residual,uoldsp,rssprk,uepsp)
      !
      ! Compute the time error using the embedded pair solutions
      !
      call get_embedded_pair_time_error_wtd_nrm(n_solpts,usp,uepsp,uoldsp, &
                                                error_dt(1:nq,1))
     !call get_embedded_pair_time_error_reg_nrm(n_solpts,usp,uepsp, &
     !                                          error_dt(1:nq,1))
      !
    else
      !
      if (EPERK_is_FSAL) then
        !
        call update_Optimized_EPERK_FSAL(n_solpts,dtsp,residual,uepsp)
        !
      end if
      !
      call get_embedded_pair_time_error_wtd_nrm(n_solpts,usp,uepsp,uoldsp, &
                                                error_dt(1:nq,1))
     !call get_embedded_pair_time_error_reg_nrm(n_solpts,usp,uepsp, &
     !                                          error_dt(1:nq,1))
      !
    end if
    !
   !call compute_adaptive_time_step(n_solpts,error_dt,dtsp,next_dtmin)
    call compute_adaptive_time_step_with_rejection(n_solpts,error_dt, &
                                                   dtsp,next_dtmin, &
                                                   previous_timestep_rejected)
    !
    ! Perform some checks on next_dtmin to make sure the time of the new
    ! solution doesnt skip past a milestone for output or final solution time
    !
    if (.not. previous_timestep_rejected) then
      !
      dtmin_save = next_dtmin
      !
      ! The checks on dtmin are dependent on time including the dtmin from
      ! the time step that was just completed. However, we need to subtract
      ! dtmin from time after performing these checks since time will be
      ! updated at a later point.
      !
      time = time + dtmin
      call checks_on_dtmin(next_dtmin,dtsp,final_timestep_from_stop_request)
      time = time - dtmin
      !
      ! Check to see if dtmin was modified in checks_on_dtmin. If it
      ! was modified, save the unmodified value to the module variable
      ! previous_dtmin to use when computing the next time step size.
      !
      if (abs(dtmin_save-next_dtmin) > eps6) then
        dtmin_was_modified = true
        previous_dtmin = dtmin_save
      end if
      !
    else
      !
      ! If the previous time step was rejected, reset usp back to uoldsp and
      ! residual back to its value based on uoldsp.
      !
      do m = 1,nq
        do n = 1,n_solpts
          usp(n,m) = uoldsp(n,m)
          residual(n,m) = rhs_save(n,m)
        end do
      end do
      !
    end if
    !
  end if
  !
  ! Update d_t with the new dtmin
  !
  d_t = next_dtmin
  !
  call debug_timer(leaving_procedure,pname)
  !
  call solver_timer(leaving_procedure,pname,idx)
  !
end subroutine get_EPERK_adaptive_dtsp
!
!###############################################################################
!
subroutine get_embedded_pair_time_error_dop853(np,uoldsp,usp,rssprk,dt_error)
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_global_solpts
 !use ovar,        only : err_dt_atol,err_dt_rtol
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp), intent(inout) :: dt_error(:)
  !
  !.. Local Scalars ..
  integer  :: n,m
  real(wp) :: lsum1,lerr1,u0,lden
  real(wp) :: lsum2,lerr2,u1
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_embedded_pair_time_error_dop853"
  !
continue
  !
  ! Compute the error for this partition
  !
  do m = 1,nq
    !
    lsum1 = zero
    lsum2 = zero
    !
    do n = 1,np
      !
      u1 = usp(n,m)
      u0 = uoldsp(n,m)
      !
      lden = err_dt_atol + err_dt_rtol*max( abs(u1) , abs(u0) )
      !
      lerr2 = rssprk(n,m,2)/lden
      lerr1 = rssprk(n,m,3)/lden
      !
      lsum1 = lsum1 + lerr1*lerr1
      lsum2 = lsum2 + lerr2*lerr2
      !
    end do
    !
    my_error(   m) = lsum1
    my_error(nq+m) = lsum2
    !
  end do
  !
  ! Collect the computed error sum across all MPI processes
  !
  if (ncpu > 1) then
    call mpi_allgather(my_error ,size(my_error,kind=int_mpi),mpi_flttyp, &
                       all_error,size(my_error,kind=int_mpi),mpi_flttyp, &
                       MPI_COMM_WORLD,mpierr)
  else
    all_error(1:2*nq,1) = my_error(1:2*nq)
  end if
  !
  ! Compute the final embedded pair time error
  !
  do m = 1,nq
    !
    lsum1 = sum(all_error(   m,:))
    lsum2 = sum(all_error(nq+m,:))
    !
    lden = lsum1 + eps2*lsum2
    if (lden <= zero) lden = one
    !
    lden = lden*real(n_global_solpts,kind=wp)
    !
    dt_error(m) = lsum1 / sqrt(lden)
    !
  end do
  !
end subroutine get_embedded_pair_time_error_dop853
!
!###############################################################################
!
subroutine get_embedded_pair_time_error_wtd_nrm(np,usp,uepsp,uoldsp,dt_error)
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_global_solpts
 !use ovar,        only : err_dt_atol,err_dt_rtol
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: usp(:,:)
  real(wp),    intent(in) :: uepsp(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: dt_error(:)
  !
  !.. Local Scalars ..
  integer  :: n,m
  real(wp) :: lsum,lnum,lden,lerr,up,um,u0
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_embedded_pair_time_error_wtd_nrm"
  !
continue
  !
  ! Compute the error for this partition
  !
  do m = 1,nq
    !
    lsum = zero
    !
    do n = 1,np
      !
      up = usp(n,m)
      um = uepsp(n,m)
      u0 = uoldsp(n,m)
      !
      lnum = up-um
     !lden = err_dt_atol + err_dt_rtol*max( abs(up) , abs(um) )
      lden = err_dt_atol + err_dt_rtol*max( abs(up) , abs(u0) )
      !
      lerr = lnum/lden
      !
      lsum = lsum + lerr*lerr
      !
    end do
    !
    my_error(m) = lsum
    !
  end do
  !
  ! Collect the computed error sum across all MPI processes
  !
  if (ncpu > 1) then
    call mpi_allgather(my_error ,size(my_error,kind=int_mpi),mpi_flttyp, &
                       all_error,size(my_error,kind=int_mpi),mpi_flttyp, &
                       MPI_COMM_WORLD,mpierr)
  else
    all_error(1:nq,1) = my_error(1:nq)
  end if
  !
  ! Compute the final embedded pair time error
  !
  do m = 1,nq
    !
    lsum = sum(all_error(m,:)) / real(n_global_solpts,kind=wp)
    !
    dt_error(m) = sqrt(lsum)
    !
  end do
  !
end subroutine get_embedded_pair_time_error_wtd_nrm
!
!###############################################################################
!
subroutine get_embedded_pair_time_error_reg_nrm(np,usp,uepsp,dt_error)
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_global_solpts
 !use ovar,        only : err_dt_atol
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: np
  real(wp), intent(in) :: usp(:,:)
  real(wp), intent(in) :: uepsp(:,:)
  real(wp), intent(inout) :: dt_error(:)
  !
  !.. Local Scalars ..
  integer  :: n,m
  real(wp) :: lsum,lerr,up,um
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_embedded_pair_time_error_reg_nrm"
  !
continue
  !
  ! Compute the error for this partition
  !
  do m = 1,nq
    !
    lsum = zero
    !
    do n = 1,np
      !
      up = usp(n,m)
      um = uepsp(n,m)
      !
      lerr = up-um
      !
      lsum = lsum + lerr*lerr
      !
    end do
    !
    my_error(m) = lsum
    !
  end do
  !
  ! Collect the computed error sum across all MPI processes
  !
  if (ncpu > 1) then
    call mpi_allgather(my_error ,size(my_error,kind=int_mpi),mpi_flttyp, &
                       all_error,size(my_error,kind=int_mpi),mpi_flttyp, &
                       MPI_COMM_WORLD,mpierr)
  else
    all_error(1:nq,1) = my_error(1:nq)
  end if
  !
  ! Compute the final embedded pair time error
  !
  do m = 1,nq
    !
    lsum = sum(all_error(m,:))
    lsum = lsum / real(n_global_solpts,kind=wp)
    !
    lerr = sqrt(lsum)
    !
   !dt_error(m) = lerr/eps9
    dt_error(m) = lerr/err_dt_atol
    !
  end do
  !
end subroutine get_embedded_pair_time_error_reg_nrm
!
!###############################################################################
!
subroutine compute_adaptive_time_step(np,dt_error,dtsp,dtmin)
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_global_solpts
 !use ovar,        only : err_dt_atol
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dt_error(:,:)
  real(wp), intent(inout) :: dtsp(:)
  real(wp), intent(inout) :: dtmin
  !
  !.. Local Scalars ..
  integer  :: i,m,n
  real(wp) :: ltot,dtnew
  !
  !.. Local Arrays ..
  real(wp) :: delta_t(1:np)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_adaptive_time_step"
  !
continue
  !
  do m = 1,nq
    ltot = one
    do i = 1,n_controllers
      ltot = ltot * (dt_error(m,i)**(pid_beta(i)*order_exp))
    end do
    delta_t(m) = ltot
  end do
  !
  dtnew = safety_factor * dtmin / maxval(delta_t(1:nq))
  !
  dtmin = dtnew
  !
  do n = 1,np
    dtsp(n) = dtmin
  end do
  !
end subroutine compute_adaptive_time_step
!
!###############################################################################
!
subroutine compute_adaptive_time_step_with_rejection(np,dt_error,dtsp,dtmin, &
                                                     previous_timestep_rejected)
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_global_solpts
 !use ovar,        only : err_dt_atol
  use ovar,        only : itcur
  !
  !.. Formal Arguments ..
  integer,      intent(in) :: np
  real(wp),     intent(in) :: dt_error(:,:)
  real(wp),  intent(inout) :: dtsp(:)
  real(wp),  intent(inout) :: dtmin
  logical(lk), intent(out) :: previous_timestep_rejected
  !
  !.. Local Scalars ..
  integer  :: i,m,n
  real(wp) :: ltot,dtnew,dt_m,errdt,errdt_prev,dt
  !
  !.. Local Arrays ..
  real(wp) :: delta_t(1:np)
  !
  !.. Local Saved Scalars ..
  real(wp), save :: previous_dt_m = eps4
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = &
                                     "compute_adaptive_time_step_with_rejection"
  !
continue
  !
  do m = 1,nq
    ltot = one
    do i = 1,n_controllers
      ltot = ltot * (dt_error(m,i)**(pid_beta(i)*order_exp))
    end do
    delta_t(m) = ltot
  end do
  !
  ! NOTE: dtnew = dtmin/delta_t, we want to find the max delta_t to get
  !       the minimum time step size required by the equations.
  !
  m = maxloc( delta_t(1:nq) , dim=1 )
 !if (mypnum == glb_root) write (iout,'(2x,i1)') m
  dt_m = dt_error(m,1) ! equivalent to err   from dopri5
  errdt = delta_t(m)   ! equivalent to fac11 from dopri5
  !
 !if (mypnum == glb_root) then
 !  write (9300,1) this_iter,itcur,(delta_t(m), m=1,nq)
 !  write (9301,1) this_iter,itcur,(dt_error(m,1), m=1,nq)
 !end if
 !1 format (i6,i6,2x,5es14.6)
  !
  if (dt_m <= one) then
    !
    ! This new time step size is accepted
    !
    previous_timestep_rejected = fals
    !
    errdt_prev = previous_dt_m**(0.2_wp*order_exp) ! unsure where 0.2 comes from
    !
    ! The following two are equivalent to one/fac from dopri5
    ! NOTE: The one5 and ten come from trying to limit the new time step
    !       size adjustment to
    !            0.2 <= 0.9*dt <= 10.0
    !       where 0.9 is some safety factor
    !
    ltot = errdt_prev / errdt
   !dt = min( ten , max(one5,safety_factor*ltot) )
    dt = min( safety_max , max(safety_min,safety_factor*ltot) )
    !
    ! Save the previous dt_m (err from dopri5) for the next time step
    !
   !if (mypnum == glb_root) then
   !  write (9302,1) this_iter,itcur,previous_dt_m,errdt_prev,errdt,ltot,dt
   !end if
    previous_dt_m = max(dt_m,eps4)
    !
    dtnew = dt * dtmin
    !
  else
    !
    ! Reject the previous time step
    !
    previous_timestep_rejected = true
    !
    dtnew = dtmin / min(five,errdt/safety_factor)
    !
  end if
  !
  dtmin = dtnew
  !
  do n = 1,np
    dtsp(n) = dtmin
  end do
  !
end subroutine compute_adaptive_time_step_with_rejection
!
!###############################################################################
!
subroutine update_EPERK_DP54_S1(np,dtsp,residual,uoldsp,usp,rssprk, &
                                pseudo_time,rhs_save)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  real(wp), intent(inout) :: pseudo_time
  real(wp), intent(inout) :: rhs_save(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP54_S1"
  !
  real(wp), parameter :: a1 = 1.0_wp/5.0_wp
  real(wp), parameter :: c1 = 1.0_wp/5.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Update pseudo_time
  !
  pseudo_time = pseudo_time + c1
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp,rhs_save)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = dtsp(n)*residual(n,nec)
      k1mx = dtsp(n)*residual(n,nmx)
      k1my = dtsp(n)*residual(n,nmy)
      k1mz = dtsp(n)*residual(n,nmz)
      k1ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + a1*k1ec
      usp(n,nmx) = u0mx + a1*k1mx
      usp(n,nmy) = u0my + a1*k1my
      usp(n,nmz) = u0mz + a1*k1mz
      usp(n,nee) = u0ee + a1*k1ee
      !
      rssprk(n,nec,1) = k1ec
      rssprk(n,nmx,1) = k1mx
      rssprk(n,nmy,1) = k1my
      rssprk(n,nmz,1) = k1mz
      rssprk(n,nee,1) = k1ee
      !
      rhs_save(n,nec) = residual(n,nec)
      rhs_save(n,nmx) = residual(n,nmx)
      rhs_save(n,nmy) = residual(n,nmy)
      rhs_save(n,nmz) = residual(n,nmz)
      rhs_save(n,nee) = residual(n,nee)
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp,rhs_save)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = dtsp(n)*residual(n,nec)
      k1mx = dtsp(n)*residual(n,nmx)
      k1my = dtsp(n)*residual(n,nmy)
      k1ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + a1*k1ec
      usp(n,nmx) = u0mx + a1*k1mx
      usp(n,nmy) = u0my + a1*k1my
      usp(n,nee) = u0ee + a1*k1ee
      !
      rssprk(n,nec,1) = k1ec
      rssprk(n,nmx,1) = k1mx
      rssprk(n,nmy,1) = k1my
      rssprk(n,nee,1) = k1ee
      !
      rhs_save(n,nec) = residual(n,nec)
      rhs_save(n,nmx) = residual(n,nmx)
      rhs_save(n,nmy) = residual(n,nmy)
      rhs_save(n,nee) = residual(n,nee)
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    !$ACC&                          PRESENT(residual,rhs_save)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,1) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) + a1*rssprk(n,m,1)
        !
        rhs_save(n,m) = residual(n,m)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP54_S1
!
!###############################################################################
!
subroutine update_EPERK_DP54_S2(np,dtsp,residual,uoldsp,usp,rssprk,pseudo_time)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  real(wp), intent(inout) :: pseudo_time
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k2ec,k2mx,k2my,k2mz,k2ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP54_S2"
  !
  real(wp), parameter :: a1 = 3.0_wp/40.0_wp
  real(wp), parameter :: a2 = 9.0_wp/40.0_wp
  real(wp), parameter :: c2 = 1.0_wp/10.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Update pseudo_time
  !
  pseudo_time = pseudo_time + c2
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1mz = a1*rssprk(n,nmz,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k2ec = dtsp(n)*residual(n,nec)
      k2mx = dtsp(n)*residual(n,nmx)
      k2my = dtsp(n)*residual(n,nmy)
      k2mz = dtsp(n)*residual(n,nmz)
      k2ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + a2*k2ec
      usp(n,nmx) = u0mx + k1mx + a2*k2mx
      usp(n,nmy) = u0my + k1my + a2*k2my
      usp(n,nmz) = u0mz + k1mz + a2*k2mz
      usp(n,nee) = u0ee + k1ee + a2*k2ee
      !
      rssprk(n,nec,2) = k2ec
      rssprk(n,nmx,2) = k2mx
      rssprk(n,nmy,2) = k2my
      rssprk(n,nmz,2) = k2mz
      rssprk(n,nee,2) = k2ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k2ec = dtsp(n)*residual(n,nec)
      k2mx = dtsp(n)*residual(n,nmx)
      k2my = dtsp(n)*residual(n,nmy)
      k2ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + a2*k2ec
      usp(n,nmx) = u0mx + k1mx + a2*k2mx
      usp(n,nmy) = u0my + k1my + a2*k2my
      usp(n,nee) = u0ee + k1ee + a2*k2ee
      !
      rssprk(n,nec,2) = k2ec
      rssprk(n,nmx,2) = k2mx
      rssprk(n,nmy,2) = k2my
      rssprk(n,nee,2) = k2ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,2) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) + a1*rssprk(n,m,1) &
                               + a2*rssprk(n,m,2)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP54_S2
!
!###############################################################################
!
subroutine update_EPERK_DP54_S3(np,dtsp,residual,uoldsp,usp,rssprk,pseudo_time)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  real(wp), intent(inout) :: pseudo_time
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k2ec,k2mx,k2my,k2mz,k2ee
  real(wp) :: k3ec,k3mx,k3my,k3mz,k3ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP54_S3"
  !
  real(wp), parameter :: a1 = 44.0_wp/45.0_wp
  real(wp), parameter :: a2 = -56.0_wp/15.0_wp
  real(wp), parameter :: a3 = 32.0_wp/9.0_wp
  real(wp), parameter :: c3 = half
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Update pseudo_time
  !
  pseudo_time = pseudo_time + c3
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1mz = a1*rssprk(n,nmz,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k2ec = a2*rssprk(n,nec,2)
      k2mx = a2*rssprk(n,nmx,2)
      k2my = a2*rssprk(n,nmy,2)
      k2mz = a2*rssprk(n,nmz,2)
      k2ee = a2*rssprk(n,nee,2)
      !
      k3ec = dtsp(n)*residual(n,nec)
      k3mx = dtsp(n)*residual(n,nmx)
      k3my = dtsp(n)*residual(n,nmy)
      k3mz = dtsp(n)*residual(n,nmz)
      k3ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k2ec + a3*k3ec
      usp(n,nmx) = u0mx + k1mx + k2mx + a3*k3mx
      usp(n,nmy) = u0my + k1my + k2my + a3*k3my
      usp(n,nmz) = u0mz + k1mz + k2mz + a3*k3mz
      usp(n,nee) = u0ee + k1ee + k2ee + a3*k3ee
      !
      rssprk(n,nec,3) = k3ec
      rssprk(n,nmx,3) = k3mx
      rssprk(n,nmy,3) = k3my
      rssprk(n,nmz,3) = k3mz
      rssprk(n,nee,3) = k3ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k2ec = a2*rssprk(n,nec,2)
      k2mx = a2*rssprk(n,nmx,2)
      k2my = a2*rssprk(n,nmy,2)
      k2ee = a2*rssprk(n,nee,2)
      !
      k3ec = dtsp(n)*residual(n,nec)
      k3mx = dtsp(n)*residual(n,nmx)
      k3my = dtsp(n)*residual(n,nmy)
      k3ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k2ec + a3*k3ec
      usp(n,nmx) = u0mx + k1mx + k2mx + a3*k3mx
      usp(n,nmy) = u0my + k1my + k2my + a3*k3my
      usp(n,nee) = u0ee + k1ee + k2ee + a3*k3ee
      !
      rssprk(n,nec,3) = k3ec
      rssprk(n,nmx,3) = k3mx
      rssprk(n,nmy,3) = k3my
      rssprk(n,nee,3) = k3ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,3) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) + a1*rssprk(n,m,1) &
                               + a2*rssprk(n,m,2) &
                               + a3*rssprk(n,m,3)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP54_S3
!
!###############################################################################
!
subroutine update_EPERK_DP54_S4(np,dtsp,residual,uoldsp,usp,rssprk,pseudo_time)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  real(wp), intent(inout) :: pseudo_time
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k2ec,k2mx,k2my,k2mz,k2ee
  real(wp) :: k3ec,k3mx,k3my,k3mz,k3ee
  real(wp) :: k4ec,k4mx,k4my,k4mz,k4ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP54_S4"
  !
  real(wp), parameter :: a1 = 19372.0_wp/6561.0_wp
  real(wp), parameter :: a2 = -25360.0_wp/2187.0_wp
  real(wp), parameter :: a3 = 64448.0_wp/6561.0_wp
  real(wp), parameter :: a4 = -212.0_wp/729.0_wp
  real(wp), parameter :: c4 = 4.0_wp/45.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Update pseudo_time
  !
  pseudo_time = pseudo_time + c4
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1mz = a1*rssprk(n,nmz,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k2ec = a2*rssprk(n,nec,2)
      k2mx = a2*rssprk(n,nmx,2)
      k2my = a2*rssprk(n,nmy,2)
      k2mz = a2*rssprk(n,nmz,2)
      k2ee = a2*rssprk(n,nee,2)
      !
      k3ec = a3*rssprk(n,nec,3)
      k3mx = a3*rssprk(n,nmx,3)
      k3my = a3*rssprk(n,nmy,3)
      k3mz = a3*rssprk(n,nmz,3)
      k3ee = a3*rssprk(n,nee,3)
      !
      k4ec = dtsp(n)*residual(n,nec)
      k4mx = dtsp(n)*residual(n,nmx)
      k4my = dtsp(n)*residual(n,nmy)
      k4mz = dtsp(n)*residual(n,nmz)
      k4ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k2ec + k3ec + a4*k4ec
      usp(n,nmx) = u0mx + k1mx + k2mx + k3mx + a4*k4mx
      usp(n,nmy) = u0my + k1my + k2my + k3my + a4*k4my
      usp(n,nmz) = u0mz + k1mz + k2mz + k3mz + a4*k4mz
      usp(n,nee) = u0ee + k1ee + k2ee + k3ee + a4*k4ee
      !
      rssprk(n,nec,4) = k4ec
      rssprk(n,nmx,4) = k4mx
      rssprk(n,nmy,4) = k4my
      rssprk(n,nmz,4) = k4mz
      rssprk(n,nee,4) = k4ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k2ec = a2*rssprk(n,nec,2)
      k2mx = a2*rssprk(n,nmx,2)
      k2my = a2*rssprk(n,nmy,2)
      k2ee = a2*rssprk(n,nee,2)
      !
      k3ec = a3*rssprk(n,nec,3)
      k3mx = a3*rssprk(n,nmx,3)
      k3my = a3*rssprk(n,nmy,3)
      k3ee = a3*rssprk(n,nee,3)
      !
      k4ec = dtsp(n)*residual(n,nec)
      k4mx = dtsp(n)*residual(n,nmx)
      k4my = dtsp(n)*residual(n,nmy)
      k4ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k2ec + k3ec + a4*k4ec
      usp(n,nmx) = u0mx + k1mx + k2mx + k3mx + a4*k4mx
      usp(n,nmy) = u0my + k1my + k2my + k3my + a4*k4my
      usp(n,nee) = u0ee + k1ee + k2ee + k3ee + a4*k4ee
      !
      rssprk(n,nec,4) = k4ec
      rssprk(n,nmx,4) = k4mx
      rssprk(n,nmy,4) = k4my
      rssprk(n,nee,4) = k4ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,4) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) + a1*rssprk(n,m,1) &
                               + a2*rssprk(n,m,2) &
                               + a3*rssprk(n,m,3) &
                               + a4*rssprk(n,m,4)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP54_S4
!
!###############################################################################
!
subroutine update_EPERK_DP54_S5(np,dtsp,residual,uoldsp,usp,rssprk,pseudo_time)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  real(wp), intent(inout) :: pseudo_time
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k2ec,k2mx,k2my,k2mz,k2ee
  real(wp) :: k3ec,k3mx,k3my,k3mz,k3ee
  real(wp) :: k4ec,k4mx,k4my,k4mz,k4ee
  real(wp) :: k5ec,k5mx,k5my,k5mz,k5ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP54_S5"
  !
  real(wp), parameter :: a1 = 9017.0_wp/3168.0_wp
  real(wp), parameter :: a2 = -355.0_wp/33.0_wp
  real(wp), parameter :: a3 = 46732.0_wp/5247.0_wp
  real(wp), parameter :: a4 = 49.0_wp/176.0_wp
  real(wp), parameter :: a5 = -5103.0_wp/18656.0_wp
  real(wp), parameter :: c5 = 1.0_wp/9.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Update pseudo_time
  !
  pseudo_time = pseudo_time + c5
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1mz = a1*rssprk(n,nmz,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k2ec = a2*rssprk(n,nec,2)
      k2mx = a2*rssprk(n,nmx,2)
      k2my = a2*rssprk(n,nmy,2)
      k2mz = a2*rssprk(n,nmz,2)
      k2ee = a2*rssprk(n,nee,2)
      !
      k3ec = a3*rssprk(n,nec,3)
      k3mx = a3*rssprk(n,nmx,3)
      k3my = a3*rssprk(n,nmy,3)
      k3mz = a3*rssprk(n,nmz,3)
      k3ee = a3*rssprk(n,nee,3)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4mz = a4*rssprk(n,nmz,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = dtsp(n)*residual(n,nec)
      k5mx = dtsp(n)*residual(n,nmx)
      k5my = dtsp(n)*residual(n,nmy)
      k5mz = dtsp(n)*residual(n,nmz)
      k5ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k2ec + k3ec + k4ec + a5*k5ec
      usp(n,nmx) = u0mx + k1mx + k2mx + k3mx + k4mx + a5*k5mx
      usp(n,nmy) = u0my + k1my + k2my + k3my + k4my + a5*k5my
      usp(n,nmz) = u0mz + k1mz + k2mz + k3mz + k4mz + a5*k5mz
      usp(n,nee) = u0ee + k1ee + k2ee + k3ee + k4ee + a5*k5ee
      !
      rssprk(n,nec,5) = k5ec
      rssprk(n,nmx,5) = k5mx
      rssprk(n,nmy,5) = k5my
      rssprk(n,nmz,5) = k5mz
      rssprk(n,nee,5) = k5ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k2ec = a2*rssprk(n,nec,2)
      k2mx = a2*rssprk(n,nmx,2)
      k2my = a2*rssprk(n,nmy,2)
      k2ee = a2*rssprk(n,nee,2)
      !
      k3ec = a3*rssprk(n,nec,3)
      k3mx = a3*rssprk(n,nmx,3)
      k3my = a3*rssprk(n,nmy,3)
      k3ee = a3*rssprk(n,nee,3)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = dtsp(n)*residual(n,nec)
      k5mx = dtsp(n)*residual(n,nmx)
      k5my = dtsp(n)*residual(n,nmy)
      k5ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k2ec + k3ec + k4ec + a5*k5ec
      usp(n,nmx) = u0mx + k1mx + k2mx + k3mx + k4mx + a5*k5mx
      usp(n,nmy) = u0my + k1my + k2my + k3my + k4my + a5*k5my
      usp(n,nee) = u0ee + k1ee + k2ee + k3ee + k4ee + a5*k5ee
      !
      rssprk(n,nec,5) = k5ec
      rssprk(n,nmx,5) = k5mx
      rssprk(n,nmy,5) = k5my
      rssprk(n,nee,5) = k5ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,5) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) + a1*rssprk(n,m,1) &
                               + a2*rssprk(n,m,2) &
                               + a3*rssprk(n,m,3) &
                               + a4*rssprk(n,m,4) &
                               + a5*rssprk(n,m,5)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP54_S5
!
!###############################################################################
!
subroutine update_EPERK_DP54_S6(np,dtsp,residual,uoldsp,usp,rssprk)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k3ec,k3mx,k3my,k3mz,k3ee
  real(wp) :: k4ec,k4mx,k4my,k4mz,k4ee
  real(wp) :: k5ec,k5mx,k5my,k5mz,k5ee
  real(wp) :: k6ec,k6mx,k6my,k6mz,k6ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP54_S6"
  !
  real(wp), parameter :: a1 = 35.0_wp/384.0_wp
  real(wp), parameter :: a3 = 500.0_wp/1113.0_wp
  real(wp), parameter :: a4 = 125.0_wp/192.0_wp
  real(wp), parameter :: a5 = -2187.0_wp/6784.0_wp
  real(wp), parameter :: a6 = 11.0_wp/84.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1mz = a1*rssprk(n,nmz,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k3ec = a3*rssprk(n,nec,3)
      k3mx = a3*rssprk(n,nmx,3)
      k3my = a3*rssprk(n,nmy,3)
      k3mz = a3*rssprk(n,nmz,3)
      k3ee = a3*rssprk(n,nee,3)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4mz = a4*rssprk(n,nmz,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = a5*rssprk(n,nec,5)
      k5mx = a5*rssprk(n,nmx,5)
      k5my = a5*rssprk(n,nmy,5)
      k5mz = a5*rssprk(n,nmz,5)
      k5ee = a5*rssprk(n,nee,5)
      !
      k6ec = dtsp(n)*residual(n,nec)
      k6mx = dtsp(n)*residual(n,nmx)
      k6my = dtsp(n)*residual(n,nmy)
      k6mz = dtsp(n)*residual(n,nmz)
      k6ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k3ec + k4ec + k5ec + a6*k6ec
      usp(n,nmx) = u0mx + k1mx + k3mx + k4mx + k5mx + a6*k6mx
      usp(n,nmy) = u0my + k1my + k3my + k4my + k5my + a6*k6my
      usp(n,nmz) = u0mz + k1mz + k3mz + k4mz + k5mz + a6*k6mz
      usp(n,nee) = u0ee + k1ee + k3ee + k4ee + k5ee + a6*k6ee
      !
      rssprk(n,nec,6) = k6ec
      rssprk(n,nmx,6) = k6mx
      rssprk(n,nmy,6) = k6my
      rssprk(n,nmz,6) = k6mz
      rssprk(n,nee,6) = k6ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k3ec = a3*rssprk(n,nec,3)
      k3mx = a3*rssprk(n,nmx,3)
      k3my = a3*rssprk(n,nmy,3)
      k3ee = a3*rssprk(n,nee,3)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = a5*rssprk(n,nec,5)
      k5mx = a5*rssprk(n,nmx,5)
      k5my = a5*rssprk(n,nmy,5)
      k5ee = a5*rssprk(n,nee,5)
      !
      k6ec = dtsp(n)*residual(n,nec)
      k6mx = dtsp(n)*residual(n,nmx)
      k6my = dtsp(n)*residual(n,nmy)
      k6ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k3ec + k4ec + k5ec + a6*k6ec
      usp(n,nmx) = u0mx + k1mx + k3mx + k4mx + k5mx + a6*k6mx
      usp(n,nmy) = u0my + k1my + k3my + k4my + k5my + a6*k6my
      usp(n,nee) = u0ee + k1ee + k3ee + k4ee + k5ee + a6*k6ee
      !
      rssprk(n,nec,6) = k6ec
      rssprk(n,nmx,6) = k6mx
      rssprk(n,nmy,6) = k6my
      rssprk(n,nee,6) = k6ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,6) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) + a1*rssprk(n,m,1) &
                               + a3*rssprk(n,m,3) &
                               + a4*rssprk(n,m,4) &
                               + a5*rssprk(n,m,5) &
                               + a6*rssprk(n,m,6)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP54_S6
!
!###############################################################################
!
subroutine update_EPERK_DP54_S7(np,dtsp,residual,uoldsp,rssprk,uepsp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp), intent(inout) :: uepsp(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k3ec,k3mx,k3my,k3mz,k3ee
  real(wp) :: k4ec,k4mx,k4my,k4mz,k4ee
  real(wp) :: k5ec,k5mx,k5my,k5mz,k5ee
  real(wp) :: k6ec,k6mx,k6my,k6mz,k6ee
  real(wp) :: k7ec,k7mx,k7my,k7mz,k7ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP54_S7"
  !
  real(wp), parameter :: a1 = 5179.0_wp/57600.0_wp
  real(wp), parameter :: a3 = 7571.0_wp/16695.0_wp
  real(wp), parameter :: a4 = 393.0_wp/640.0_wp
  real(wp), parameter :: a5 = -92097.0_wp/339200.0_wp
  real(wp), parameter :: a6 = 187.0_wp/2100.0_wp
  real(wp), parameter :: a7 = 1.0_wp/40.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(uepsp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1mz = a1*rssprk(n,nmz,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k3ec = a3*rssprk(n,nec,3)
      k3mx = a3*rssprk(n,nmx,3)
      k3my = a3*rssprk(n,nmy,3)
      k3mz = a3*rssprk(n,nmz,3)
      k3ee = a3*rssprk(n,nee,3)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4mz = a4*rssprk(n,nmz,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = a5*rssprk(n,nec,5)
      k5mx = a5*rssprk(n,nmx,5)
      k5my = a5*rssprk(n,nmy,5)
      k5mz = a5*rssprk(n,nmz,5)
      k5ee = a5*rssprk(n,nee,5)
      !
      k6ec = a6*rssprk(n,nec,6)
      k6mx = a6*rssprk(n,nmx,6)
      k6my = a6*rssprk(n,nmy,6)
      k6mz = a6*rssprk(n,nmz,6)
      k6ee = a6*rssprk(n,nee,6)
      !
      k7ec = dtsp(n)*residual(n,nec)
      k7mx = dtsp(n)*residual(n,nmx)
      k7my = dtsp(n)*residual(n,nmy)
      k7mz = dtsp(n)*residual(n,nmz)
      k7ee = dtsp(n)*residual(n,nee)
      !
      uepsp(n,nec) = u0ec + k1ec + k3ec + k4ec + k5ec + k6ec + a7*k7ec
      uepsp(n,nmx) = u0mx + k1mx + k3mx + k4mx + k5mx + k6mx + a7*k7mx
      uepsp(n,nmy) = u0my + k1my + k3my + k4my + k5my + k6my + a7*k7my
      uepsp(n,nmz) = u0mz + k1mz + k3mz + k4mz + k5mz + k6mz + a7*k7mz
      uepsp(n,nee) = u0ee + k1ee + k3ee + k4ee + k5ee + k6ee + a7*k7ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(uepsp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k3ec = a3*rssprk(n,nec,3)
      k3mx = a3*rssprk(n,nmx,3)
      k3my = a3*rssprk(n,nmy,3)
      k3ee = a3*rssprk(n,nee,3)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = a5*rssprk(n,nec,5)
      k5mx = a5*rssprk(n,nmx,5)
      k5my = a5*rssprk(n,nmy,5)
      k5ee = a5*rssprk(n,nee,5)
      !
      k6ec = a6*rssprk(n,nec,6)
      k6mx = a6*rssprk(n,nmx,6)
      k6my = a6*rssprk(n,nmy,6)
      k6ee = a6*rssprk(n,nee,6)
      !
      k7ec = dtsp(n)*residual(n,nec)
      k7mx = dtsp(n)*residual(n,nmx)
      k7my = dtsp(n)*residual(n,nmy)
      k7ee = dtsp(n)*residual(n,nee)
      !
      uepsp(n,nec) = u0ec + k1ec + k3ec + k4ec + k5ec + k6ec + a7*k7ec
      uepsp(n,nmx) = u0mx + k1mx + k3mx + k4mx + k5mx + k6mx + a7*k7mx
      uepsp(n,nmy) = u0my + k1my + k3my + k4my + k5my + k6my + a7*k7my
      uepsp(n,nee) = u0ee + k1ee + k3ee + k4ee + k5ee + k6ee + a7*k7ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(uepsp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        uepsp(n,m) = uoldsp(n,m) + a1*rssprk(n,m,1) &
                                 + a3*rssprk(n,m,3) &
                                 + a4*rssprk(n,m,4) &
                                 + a5*rssprk(n,m,5) &
                                 + a6*rssprk(n,m,6) &
                                 + a7*dtsp(n)*residual(n,m)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP54_S7
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine update_EPERK_DP853_S1(np,dtsp,residual,uoldsp,usp,rssprk, &
                                 pseudo_time,rhs_save)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  real(wp), intent(inout) :: pseudo_time
  real(wp), intent(inout) :: rhs_save(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP853_S1"
  !
  real(wp), parameter :: a1 = 0.0526001519587677318785587544488_wp
  real(wp), parameter :: c1 = 0.0526001519587677318785587544488_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Update pseudo_time
  !
  pseudo_time = pseudo_time + c1
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp,rhs_save)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = dtsp(n)*residual(n,nec)
      k1mx = dtsp(n)*residual(n,nmx)
      k1my = dtsp(n)*residual(n,nmy)
      k1mz = dtsp(n)*residual(n,nmz)
      k1ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + a1*k1ec
      usp(n,nmx) = u0mx + a1*k1mx
      usp(n,nmy) = u0my + a1*k1my
      usp(n,nmz) = u0mz + a1*k1mz
      usp(n,nee) = u0ee + a1*k1ee
      !
      rssprk(n,nec,1) = k1ec
      rssprk(n,nmx,1) = k1mx
      rssprk(n,nmy,1) = k1my
      rssprk(n,nmz,1) = k1mz
      rssprk(n,nee,1) = k1ee
      !
      rhs_save(n,nec) = residual(n,nec)
      rhs_save(n,nmx) = residual(n,nmx)
      rhs_save(n,nmy) = residual(n,nmy)
      rhs_save(n,nmz) = residual(n,nmz)
      rhs_save(n,nee) = residual(n,nee)
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp,rhs_save)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = dtsp(n)*residual(n,nec)
      k1mx = dtsp(n)*residual(n,nmx)
      k1my = dtsp(n)*residual(n,nmy)
      k1ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + a1*k1ec
      usp(n,nmx) = u0mx + a1*k1mx
      usp(n,nmy) = u0my + a1*k1my
      usp(n,nee) = u0ee + a1*k1ee
      !
      rssprk(n,nec,1) = k1ec
      rssprk(n,nmx,1) = k1mx
      rssprk(n,nmy,1) = k1my
      rssprk(n,nee,1) = k1ee
      !
      rhs_save(n,nec) = residual(n,nec)
      rhs_save(n,nmx) = residual(n,nmx)
      rhs_save(n,nmy) = residual(n,nmy)
      rhs_save(n,nee) = residual(n,nee)
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp) &
    !$ACC&                          PRESENT(residual,rhs_save)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,1) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) + a1*rssprk(n,m,1)
        !
        rhs_save(n,m) = residual(n,m)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP853_S1
!
!###############################################################################
!
subroutine update_EPERK_DP853_S2(np,dtsp,residual,uoldsp,usp,rssprk,pseudo_time)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  real(wp), intent(inout) :: pseudo_time
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k2ec,k2mx,k2my,k2mz,k2ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP853_S2"
  !
  real(wp), parameter :: a1 = 0.0197250569845378994544595329183_wp
  real(wp), parameter :: a2 = 0.0591751709536136983633785987549_wp
  real(wp), parameter :: c2 = 0.0789002279381515978178381316732_wp - &
                              0.0526001519587677318785587544488_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Update pseudo_time
  !
  pseudo_time = pseudo_time + c2
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1mz = a1*rssprk(n,nmz,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k2ec = dtsp(n)*residual(n,nec)
      k2mx = dtsp(n)*residual(n,nmx)
      k2my = dtsp(n)*residual(n,nmy)
      k2mz = dtsp(n)*residual(n,nmz)
      k2ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + a2*k2ec
      usp(n,nmx) = u0mx + k1mx + a2*k2mx
      usp(n,nmy) = u0my + k1my + a2*k2my
      usp(n,nmz) = u0mz + k1mz + a2*k2mz
      usp(n,nee) = u0ee + k1ee + a2*k2ee
      !
      rssprk(n,nec,2) = k2ec
      rssprk(n,nmx,2) = k2mx
      rssprk(n,nmy,2) = k2my
      rssprk(n,nmz,2) = k2mz
      rssprk(n,nee,2) = k2ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k2ec = dtsp(n)*residual(n,nec)
      k2mx = dtsp(n)*residual(n,nmx)
      k2my = dtsp(n)*residual(n,nmy)
      k2ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + a2*k2ec
      usp(n,nmx) = u0mx + k1mx + a2*k2mx
      usp(n,nmy) = u0my + k1my + a2*k2my
      usp(n,nee) = u0ee + k1ee + a2*k2ee
      !
      rssprk(n,nec,2) = k2ec
      rssprk(n,nmx,2) = k2mx
      rssprk(n,nmy,2) = k2my
      rssprk(n,nee,2) = k2ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,2) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) + a1*rssprk(n,m,1) &
                               + a2*rssprk(n,m,2)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP853_S2
!
!###############################################################################
!
subroutine update_EPERK_DP853_S3(np,dtsp,residual,uoldsp,usp,rssprk,pseudo_time)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  real(wp), intent(inout) :: pseudo_time
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k3ec,k3mx,k3my,k3mz,k3ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP853_S3"
  !
  real(wp), parameter :: a1 = 0.0295875854768068491816892993775_wp
  real(wp), parameter :: a3 = 0.0887627564304205475450678981324_wp
  real(wp), parameter :: c3 = 0.118350341907227396726757197510_wp - &
                              0.0789002279381515978178381316732_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Update pseudo_time
  !
  pseudo_time = pseudo_time + c3
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1mz = a1*rssprk(n,nmz,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k3ec = dtsp(n)*residual(n,nec)
      k3mx = dtsp(n)*residual(n,nmx)
      k3my = dtsp(n)*residual(n,nmy)
      k3mz = dtsp(n)*residual(n,nmz)
      k3ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + a3*k3ec
      usp(n,nmx) = u0mx + k1mx + a3*k3mx
      usp(n,nmy) = u0my + k1my + a3*k3my
      usp(n,nmz) = u0mz + k1mz + a3*k3mz
      usp(n,nee) = u0ee + k1ee + a3*k3ee
      !
      rssprk(n,nec,3) = k3ec
      rssprk(n,nmx,3) = k3mx
      rssprk(n,nmy,3) = k3my
      rssprk(n,nmz,3) = k3mz
      rssprk(n,nee,3) = k3ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k3ec = dtsp(n)*residual(n,nec)
      k3mx = dtsp(n)*residual(n,nmx)
      k3my = dtsp(n)*residual(n,nmy)
      k3ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + a3*k3ec
      usp(n,nmx) = u0mx + k1mx + a3*k3mx
      usp(n,nmy) = u0my + k1my + a3*k3my
      usp(n,nee) = u0ee + k1ee + a3*k3ee
      !
      rssprk(n,nec,3) = k3ec
      rssprk(n,nmx,3) = k3mx
      rssprk(n,nmy,3) = k3my
      rssprk(n,nee,3) = k3ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,3) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) + a1*rssprk(n,m,1) &
                               + a3*rssprk(n,m,3)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP853_S3
!
!###############################################################################
!
subroutine update_EPERK_DP853_S4(np,dtsp,residual,uoldsp,usp,rssprk,pseudo_time)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  real(wp), intent(inout) :: pseudo_time
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k3ec,k3mx,k3my,k3mz,k3ee
  real(wp) :: k4ec,k4mx,k4my,k4mz,k4ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP853_S4"
  !
  real(wp), parameter :: a1 =  0.241365134159266685502369798665_wp
  real(wp), parameter :: a3 = -0.884549479328286085344864962717_wp
  real(wp), parameter :: a4 =  0.924834003261792003115737966543_wp
  real(wp), parameter :: c4 =  0.281649658092772603273242802490_wp - &
                               0.118350341907227396726757197510_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Update pseudo_time
  !
  pseudo_time = pseudo_time + c4
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1mz = a1*rssprk(n,nmz,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k3ec = a3*rssprk(n,nec,3)
      k3mx = a3*rssprk(n,nmx,3)
      k3my = a3*rssprk(n,nmy,3)
      k3mz = a3*rssprk(n,nmz,3)
      k3ee = a3*rssprk(n,nee,3)
      !
      k4ec = dtsp(n)*residual(n,nec)
      k4mx = dtsp(n)*residual(n,nmx)
      k4my = dtsp(n)*residual(n,nmy)
      k4mz = dtsp(n)*residual(n,nmz)
      k4ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k3ec + a4*k4ec
      usp(n,nmx) = u0mx + k1mx + k3mx + a4*k4mx
      usp(n,nmy) = u0my + k1my + k3my + a4*k4my
      usp(n,nmz) = u0mz + k1mz + k3mz + a4*k4mz
      usp(n,nee) = u0ee + k1ee + k3ee + a4*k4ee
      !
      rssprk(n,nec,4) = k4ec
      rssprk(n,nmx,4) = k4mx
      rssprk(n,nmy,4) = k4my
      rssprk(n,nmz,4) = k4mz
      rssprk(n,nee,4) = k4ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k3ec = a3*rssprk(n,nec,3)
      k3mx = a3*rssprk(n,nmx,3)
      k3my = a3*rssprk(n,nmy,3)
      k3ee = a3*rssprk(n,nee,3)
      !
      k4ec = dtsp(n)*residual(n,nec)
      k4mx = dtsp(n)*residual(n,nmx)
      k4my = dtsp(n)*residual(n,nmy)
      k4ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k3ec + a4*k4ec
      usp(n,nmx) = u0mx + k1mx + k3mx + a4*k4mx
      usp(n,nmy) = u0my + k1my + k3my + a4*k4my
      usp(n,nee) = u0ee + k1ee + k3ee + a4*k4ee
      !
      rssprk(n,nec,4) = k4ec
      rssprk(n,nmx,4) = k4mx
      rssprk(n,nmy,4) = k4my
      rssprk(n,nee,4) = k4ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,4) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) + a1*rssprk(n,m,1) &
                               + a3*rssprk(n,m,3) &
                               + a4*rssprk(n,m,4)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP853_S4
!
!###############################################################################
!
subroutine update_EPERK_DP853_S5(np,dtsp,residual,uoldsp,usp,rssprk,pseudo_time)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  real(wp), intent(inout) :: pseudo_time
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k4ec,k4mx,k4my,k4mz,k4ee
  real(wp) :: k5ec,k5mx,k5my,k5mz,k5ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP853_S5"
  !
  real(wp), parameter :: a1 = 0.0370370370370370370370370370370_wp
  real(wp), parameter :: a4 = 0.170828608729473871279604482173_wp
  real(wp), parameter :: a5 = 0.125467687566822425016691814123_wp
  real(wp), parameter :: c5 = one/three - 0.281649658092772603273242802490_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Update pseudo_time
  !
  pseudo_time = pseudo_time + c5
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1mz = a1*rssprk(n,nmz,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4mz = a4*rssprk(n,nmz,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = dtsp(n)*residual(n,nec)
      k5mx = dtsp(n)*residual(n,nmx)
      k5my = dtsp(n)*residual(n,nmy)
      k5mz = dtsp(n)*residual(n,nmz)
      k5ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k4ec + a5*k5ec
      usp(n,nmx) = u0mx + k1mx + k4mx + a5*k5mx
      usp(n,nmy) = u0my + k1my + k4my + a5*k5my
      usp(n,nmz) = u0mz + k1mz + k4mz + a5*k5mz
      usp(n,nee) = u0ee + k1ee + k4ee + a5*k5ee
      !
      rssprk(n,nec,5) = k5ec
      rssprk(n,nmx,5) = k5mx
      rssprk(n,nmy,5) = k5my
      rssprk(n,nmz,5) = k5mz
      rssprk(n,nee,5) = k5ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = dtsp(n)*residual(n,nec)
      k5mx = dtsp(n)*residual(n,nmx)
      k5my = dtsp(n)*residual(n,nmy)
      k5ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k4ec + a5*k5ec
      usp(n,nmx) = u0mx + k1mx + k4mx + a5*k5mx
      usp(n,nmy) = u0my + k1my + k4my + a5*k5my
      usp(n,nee) = u0ee + k1ee + k4ee + a5*k5ee
      !
      rssprk(n,nec,5) = k5ec
      rssprk(n,nmx,5) = k5mx
      rssprk(n,nmy,5) = k5my
      rssprk(n,nee,5) = k5ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,5) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) + a1*rssprk(n,m,1) &
                               + a4*rssprk(n,m,4) &
                               + a5*rssprk(n,m,5)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP853_S5
!
!###############################################################################
!
subroutine update_EPERK_DP853_S6(np,dtsp,residual,uoldsp,usp,rssprk,pseudo_time)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  real(wp), intent(inout) :: pseudo_time
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k4ec,k4mx,k4my,k4mz,k4ee
  real(wp) :: k5ec,k5mx,k5my,k5mz,k5ee
  real(wp) :: k6ec,k6mx,k6my,k6mz,k6ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP853_S6"
  !
  real(wp), parameter :: a1 =  0.037109375_wp
  real(wp), parameter :: a4 =  0.170252211019544039314978060272_wp
  real(wp), parameter :: a5 =  0.0602165389804559606850219397283_wp
  real(wp), parameter :: a6 = -0.017578125_wp
  real(wp), parameter :: c6 =  0.25_wp - one/three
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Update pseudo_time
  !
  pseudo_time = pseudo_time + c6
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1mz = a1*rssprk(n,nmz,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4mz = a4*rssprk(n,nmz,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = a5*rssprk(n,nec,5)
      k5mx = a5*rssprk(n,nmx,5)
      k5my = a5*rssprk(n,nmy,5)
      k5mz = a5*rssprk(n,nmz,5)
      k5ee = a5*rssprk(n,nee,5)
      !
      k6ec = dtsp(n)*residual(n,nec)
      k6mx = dtsp(n)*residual(n,nmx)
      k6my = dtsp(n)*residual(n,nmy)
      k6mz = dtsp(n)*residual(n,nmz)
      k6ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k4ec + k5ec + a6*k6ec
      usp(n,nmx) = u0mx + k1mx + k4mx + k5mx + a6*k6mx
      usp(n,nmy) = u0my + k1my + k4my + k5my + a6*k6my
      usp(n,nmz) = u0mz + k1mz + k4mz + k5mz + a6*k6mz
      usp(n,nee) = u0ee + k1ee + k4ee + k5ee + a6*k6ee
      !
      rssprk(n,nec,6) = k6ec
      rssprk(n,nmx,6) = k6mx
      rssprk(n,nmy,6) = k6my
      rssprk(n,nmz,6) = k6mz
      rssprk(n,nee,6) = k6ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = a5*rssprk(n,nec,5)
      k5mx = a5*rssprk(n,nmx,5)
      k5my = a5*rssprk(n,nmy,5)
      k5ee = a5*rssprk(n,nee,5)
      !
      k6ec = dtsp(n)*residual(n,nec)
      k6mx = dtsp(n)*residual(n,nmx)
      k6my = dtsp(n)*residual(n,nmy)
      k6ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k4ec + k5ec + a6*k6ec
      usp(n,nmx) = u0mx + k1mx + k4mx + k5mx + a6*k6mx
      usp(n,nmy) = u0my + k1my + k4my + k5my + a6*k6my
      usp(n,nee) = u0ee + k1ee + k4ee + k5ee + a6*k6ee
      !
      rssprk(n,nec,6) = k6ec
      rssprk(n,nmx,6) = k6mx
      rssprk(n,nmy,6) = k6my
      rssprk(n,nee,6) = k6ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,6) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) + a1*rssprk(n,m,1) &
                               + a4*rssprk(n,m,4) &
                               + a5*rssprk(n,m,5) &
                               + a6*rssprk(n,m,6)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP853_S6
!
!###############################################################################
!
subroutine update_EPERK_DP853_S7(np,dtsp,residual,uoldsp,usp,rssprk,pseudo_time)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  real(wp), intent(inout) :: pseudo_time
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k4ec,k4mx,k4my,k4mz,k4ee
  real(wp) :: k5ec,k5mx,k5my,k5mz,k5ee
  real(wp) :: k6ec,k6mx,k6my,k6mz,k6ee
  real(wp) :: k7ec,k7mx,k7my,k7mz,k7ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP853_S7"
  !
  real(wp), parameter :: a1 =  0.0370920001185047927108779319836_wp
  real(wp), parameter :: a4 =  0.170383925712239993810214054705_wp
  real(wp), parameter :: a5 =  0.107262030446373284651809199168_wp
  real(wp), parameter :: a6 = -0.0153194377486244017527936158236_wp
  real(wp), parameter :: a7 =  0.00827378916381402288758473766002_wp
  real(wp), parameter :: c7 =  0.307692307692307692307692307692_wp - 0.25_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Update pseudo_time
  !
  pseudo_time = pseudo_time + c7
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1mz = a1*rssprk(n,nmz,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4mz = a4*rssprk(n,nmz,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = a5*rssprk(n,nec,5)
      k5mx = a5*rssprk(n,nmx,5)
      k5my = a5*rssprk(n,nmy,5)
      k5mz = a5*rssprk(n,nmz,5)
      k5ee = a5*rssprk(n,nee,5)
      !
      k6ec = a6*rssprk(n,nec,6)
      k6mx = a6*rssprk(n,nmx,6)
      k6my = a6*rssprk(n,nmy,6)
      k6mz = a6*rssprk(n,nmz,6)
      k6ee = a6*rssprk(n,nee,6)
      !
      k7ec = dtsp(n)*residual(n,nec)
      k7mx = dtsp(n)*residual(n,nmx)
      k7my = dtsp(n)*residual(n,nmy)
      k7mz = dtsp(n)*residual(n,nmz)
      k7ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k4ec + k5ec + k6ec + a7*k7ec
      usp(n,nmx) = u0mx + k1mx + k4mx + k5mx + k6mx + a7*k7mx
      usp(n,nmy) = u0my + k1my + k4my + k5my + k6my + a7*k7my
      usp(n,nmz) = u0mz + k1mz + k4mz + k5mz + k6mz + a7*k7mz
      usp(n,nee) = u0ee + k1ee + k4ee + k5ee + k6ee + a7*k7ee
      !
      rssprk(n,nec,7) = k7ec
      rssprk(n,nmx,7) = k7mx
      rssprk(n,nmy,7) = k7my
      rssprk(n,nmz,7) = k7mz
      rssprk(n,nee,7) = k7ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = a5*rssprk(n,nec,5)
      k5mx = a5*rssprk(n,nmx,5)
      k5my = a5*rssprk(n,nmy,5)
      k5ee = a5*rssprk(n,nee,5)
      !
      k6ec = a6*rssprk(n,nec,6)
      k6mx = a6*rssprk(n,nmx,6)
      k6my = a6*rssprk(n,nmy,6)
      k6ee = a6*rssprk(n,nee,6)
      !
      k7ec = dtsp(n)*residual(n,nec)
      k7mx = dtsp(n)*residual(n,nmx)
      k7my = dtsp(n)*residual(n,nmy)
      k7ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k4ec + k5ec + k6ec + a7*k7ec
      usp(n,nmx) = u0mx + k1mx + k4mx + k5mx + k6mx + a7*k7mx
      usp(n,nmy) = u0my + k1my + k4my + k5my + k6my + a7*k7my
      usp(n,nee) = u0ee + k1ee + k4ee + k5ee + k6ee + a7*k7ee
      !
      rssprk(n,nec,7) = k7ec
      rssprk(n,nmx,7) = k7mx
      rssprk(n,nmy,7) = k7my
      rssprk(n,nee,7) = k7ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,7) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) + a1*rssprk(n,m,1) &
                               + a4*rssprk(n,m,4) &
                               + a5*rssprk(n,m,5) &
                               + a6*rssprk(n,m,6) &
                               + a7*rssprk(n,m,7)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP853_S7
!
!###############################################################################
!
subroutine update_EPERK_DP853_S8(np,dtsp,residual,uoldsp,usp,rssprk,pseudo_time)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  real(wp), intent(inout) :: pseudo_time
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k4ec,k4mx,k4my,k4mz,k4ee
  real(wp) :: k5ec,k5mx,k5my,k5mz,k5ee
  real(wp) :: k6ec,k6mx,k6my,k6mz,k6ee
  real(wp) :: k7ec,k7mx,k7my,k7mz,k7ee
  real(wp) :: k8ec,k8mx,k8my,k8mz,k8ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP853_S8"
  !
  real(wp), parameter :: a1 =  0.624110958716075717114429577812_wp
  real(wp), parameter :: a4 = -3.36089262944694129406857109825_wp
  real(wp), parameter :: a5 = -0.868219346841726006818189891453_wp
  real(wp), parameter :: a6 =  27.5920996994467083049415600797_wp
  real(wp), parameter :: a7 =  20.1540675504778934086186788979_wp
  real(wp), parameter :: a8 = -43.4898841810699588477366255144_wp
  real(wp), parameter :: c8 =  0.651282051282051282051282051282_wp - &
                               0.307692307692307692307692307692_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Update pseudo_time
  !
  pseudo_time = pseudo_time + c8
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1mz = a1*rssprk(n,nmz,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4mz = a4*rssprk(n,nmz,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = a5*rssprk(n,nec,5)
      k5mx = a5*rssprk(n,nmx,5)
      k5my = a5*rssprk(n,nmy,5)
      k5mz = a5*rssprk(n,nmz,5)
      k5ee = a5*rssprk(n,nee,5)
      !
      k6ec = a6*rssprk(n,nec,6)
      k6mx = a6*rssprk(n,nmx,6)
      k6my = a6*rssprk(n,nmy,6)
      k6mz = a6*rssprk(n,nmz,6)
      k6ee = a6*rssprk(n,nee,6)
      !
      k7ec = a7*rssprk(n,nec,7)
      k7mx = a7*rssprk(n,nmx,7)
      k7my = a7*rssprk(n,nmy,7)
      k7mz = a7*rssprk(n,nmz,7)
      k7ee = a7*rssprk(n,nee,7)
      !
      k8ec = dtsp(n)*residual(n,nec)
      k8mx = dtsp(n)*residual(n,nmx)
      k8my = dtsp(n)*residual(n,nmy)
      k8mz = dtsp(n)*residual(n,nmz)
      k8ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k4ec + k5ec + k6ec + k7ec + a8*k8ec
      usp(n,nmx) = u0mx + k1mx + k4mx + k5mx + k6mx + k7mx + a8*k8mx
      usp(n,nmy) = u0my + k1my + k4my + k5my + k6my + k7my + a8*k8my
      usp(n,nmz) = u0mz + k1mz + k4mz + k5mz + k6mz + k7mz + a8*k8mz
      usp(n,nee) = u0ee + k1ee + k4ee + k5ee + k6ee + k7ee + a8*k8ee
      !
      rssprk(n,nec,8) = k8ec
      rssprk(n,nmx,8) = k8mx
      rssprk(n,nmy,8) = k8my
      rssprk(n,nmz,8) = k8mz
      rssprk(n,nee,8) = k8ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = a5*rssprk(n,nec,5)
      k5mx = a5*rssprk(n,nmx,5)
      k5my = a5*rssprk(n,nmy,5)
      k5ee = a5*rssprk(n,nee,5)
      !
      k6ec = a6*rssprk(n,nec,6)
      k6mx = a6*rssprk(n,nmx,6)
      k6my = a6*rssprk(n,nmy,6)
      k6ee = a6*rssprk(n,nee,6)
      !
      k7ec = a7*rssprk(n,nec,7)
      k7mx = a7*rssprk(n,nmx,7)
      k7my = a7*rssprk(n,nmy,7)
      k7ee = a7*rssprk(n,nee,7)
      !
      k8ec = dtsp(n)*residual(n,nec)
      k8mx = dtsp(n)*residual(n,nmx)
      k8my = dtsp(n)*residual(n,nmy)
      k8ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k4ec + k5ec + k6ec + k7ec + a8*k8ec
      usp(n,nmx) = u0mx + k1mx + k4mx + k5mx + k6mx + k7mx + a8*k8mx
      usp(n,nmy) = u0my + k1my + k4my + k5my + k6my + k7my + a8*k8my
      usp(n,nee) = u0ee + k1ee + k4ee + k5ee + k6ee + k7ee + a8*k8ee
      !
      rssprk(n,nec,8) = k8ec
      rssprk(n,nmx,8) = k8mx
      rssprk(n,nmy,8) = k8my
      rssprk(n,nee,8) = k8ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,8) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) + a1*rssprk(n,m,1) &
                               + a4*rssprk(n,m,4) &
                               + a5*rssprk(n,m,5) &
                               + a6*rssprk(n,m,6) &
                               + a7*rssprk(n,m,7) &
                               + a8*rssprk(n,m,8)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP853_S8
!
!###############################################################################
!
subroutine update_EPERK_DP853_S9(np,dtsp,residual,uoldsp,usp,rssprk,pseudo_time)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  real(wp), intent(inout) :: pseudo_time
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k4ec,k4mx,k4my,k4mz,k4ee
  real(wp) :: k5ec,k5mx,k5my,k5mz,k5ee
  real(wp) :: k6ec,k6mx,k6my,k6mz,k6ee
  real(wp) :: k7ec,k7mx,k7my,k7mz,k7ee
  real(wp) :: k8ec,k8mx,k8my,k8mz,k8ee
  real(wp) :: k9ec,k9mx,k9my,k9mz,k9ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP853_S9"
  !
  real(wp), parameter :: a1 =  0.477662536438264365890433908527_wp
  real(wp), parameter :: a4 = -2.48811461997166764192642586468_wp
  real(wp), parameter :: a5 = -0.590290826836842996371446475743_wp
  real(wp), parameter :: a6 =  21.2300514481811942347288949897_wp
  real(wp), parameter :: a7 =  15.2792336328824235832596922938_wp
  real(wp), parameter :: a8 = -33.2882109689848629194453265587_wp
  real(wp), parameter :: a9 = -0.0203312017085086261358222928593_wp
  real(wp), parameter :: c9 =  0.6_wp - 0.651282051282051282051282051282_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Update pseudo_time
  !
  pseudo_time = pseudo_time + c9
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1mz = a1*rssprk(n,nmz,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4mz = a4*rssprk(n,nmz,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = a5*rssprk(n,nec,5)
      k5mx = a5*rssprk(n,nmx,5)
      k5my = a5*rssprk(n,nmy,5)
      k5mz = a5*rssprk(n,nmz,5)
      k5ee = a5*rssprk(n,nee,5)
      !
      k6ec = a6*rssprk(n,nec,6)
      k6mx = a6*rssprk(n,nmx,6)
      k6my = a6*rssprk(n,nmy,6)
      k6mz = a6*rssprk(n,nmz,6)
      k6ee = a6*rssprk(n,nee,6)
      !
      k7ec = a7*rssprk(n,nec,7)
      k7mx = a7*rssprk(n,nmx,7)
      k7my = a7*rssprk(n,nmy,7)
      k7mz = a7*rssprk(n,nmz,7)
      k7ee = a7*rssprk(n,nee,7)
      !
      k8ec = a8*rssprk(n,nec,8)
      k8mx = a8*rssprk(n,nmx,8)
      k8my = a8*rssprk(n,nmy,8)
      k8mz = a8*rssprk(n,nmz,8)
      k8ee = a8*rssprk(n,nee,8)
      !
      k9ec = dtsp(n)*residual(n,nec)
      k9mx = dtsp(n)*residual(n,nmx)
      k9my = dtsp(n)*residual(n,nmy)
      k9mz = dtsp(n)*residual(n,nmz)
      k9ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k4ec + k5ec + k6ec + k7ec + k8ec + a9*k9ec
      usp(n,nmx) = u0mx + k1mx + k4mx + k5mx + k6mx + k7mx + k8mx + a9*k9mx
      usp(n,nmy) = u0my + k1my + k4my + k5my + k6my + k7my + k8my + a9*k9my
      usp(n,nmz) = u0mz + k1mz + k4mz + k5mz + k6mz + k7mz + k8mz + a9*k9mz
      usp(n,nee) = u0ee + k1ee + k4ee + k5ee + k6ee + k7ee + k8ee + a9*k9ee
      !
      rssprk(n,nec,9) = k9ec
      rssprk(n,nmx,9) = k9mx
      rssprk(n,nmy,9) = k9my
      rssprk(n,nmz,9) = k9mz
      rssprk(n,nee,9) = k9ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = a5*rssprk(n,nec,5)
      k5mx = a5*rssprk(n,nmx,5)
      k5my = a5*rssprk(n,nmy,5)
      k5ee = a5*rssprk(n,nee,5)
      !
      k6ec = a6*rssprk(n,nec,6)
      k6mx = a6*rssprk(n,nmx,6)
      k6my = a6*rssprk(n,nmy,6)
      k6ee = a6*rssprk(n,nee,6)
      !
      k7ec = a7*rssprk(n,nec,7)
      k7mx = a7*rssprk(n,nmx,7)
      k7my = a7*rssprk(n,nmy,7)
      k7ee = a7*rssprk(n,nee,7)
      !
      k8ec = a8*rssprk(n,nec,8)
      k8mx = a8*rssprk(n,nmx,8)
      k8my = a8*rssprk(n,nmy,8)
      k8ee = a8*rssprk(n,nee,8)
      !
      k9ec = dtsp(n)*residual(n,nec)
      k9mx = dtsp(n)*residual(n,nmx)
      k9my = dtsp(n)*residual(n,nmy)
      k9ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec + k1ec + k4ec + k5ec + k6ec + k7ec + k8ec + a9*k9ec
      usp(n,nmx) = u0mx + k1mx + k4mx + k5mx + k6mx + k7mx + k8mx + a9*k9mx
      usp(n,nmy) = u0my + k1my + k4my + k5my + k6my + k7my + k8my + a9*k9my
      usp(n,nee) = u0ee + k1ee + k4ee + k5ee + k6ee + k7ee + k8ee + a9*k9ee
      !
      rssprk(n,nec,9) = k9ec
      rssprk(n,nmx,9) = k9mx
      rssprk(n,nmy,9) = k9my
      rssprk(n,nee,9) = k9ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,9) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) + a1*rssprk(n,m,1) &
                               + a4*rssprk(n,m,4) &
                               + a5*rssprk(n,m,5) &
                               + a6*rssprk(n,m,6) &
                               + a7*rssprk(n,m,7) &
                               + a8*rssprk(n,m,8) &
                               + a9*rssprk(n,m,9)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP853_S9
!
!###############################################################################
!
subroutine update_EPERK_DP853_S10(np,dtsp,residual,uoldsp,usp,rssprk, &
                                  pseudo_time)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  real(wp), intent(inout) :: pseudo_time
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k4ec,k4mx,k4my,k4mz,k4ee
  real(wp) :: k5ec,k5mx,k5my,k5mz,k5ee
  real(wp) :: k6ec,k6mx,k6my,k6mz,k6ee
  real(wp) :: k7ec,k7mx,k7my,k7mz,k7ee
  real(wp) :: k8ec,k8mx,k8my,k8mz,k8ee
  real(wp) :: k9ec,k9mx,k9my,k9mz,k9ee
  real(wp) :: k10ec,k10mx,k10my,k10mz,k10ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP853_S10"
  !
  real(wp), parameter :: a1  = -0.93714243008598732571704021658_wp
  real(wp), parameter :: a4  =  5.18637242884406370830023853209_wp
  real(wp), parameter :: a5  =  1.09143734899672957818500254654_wp
  real(wp), parameter :: a6  = -8.14978701074692612513997267357_wp
  real(wp), parameter :: a7  = -18.5200656599969598641566180701_wp
  real(wp), parameter :: a8  =  22.7394870993505042818970056734_wp
  real(wp), parameter :: a9  =  2.49360555267965238987089396762_wp
  real(wp), parameter :: a10 = -3.0467644718982195003823669022_wp
  real(wp), parameter :: c10 =  0.857142857142857142857142857142_wp - 0.6_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Update pseudo_time
  !
  pseudo_time = pseudo_time + c10
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1mz = a1*rssprk(n,nmz,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4mz = a4*rssprk(n,nmz,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = a5*rssprk(n,nec,5)
      k5mx = a5*rssprk(n,nmx,5)
      k5my = a5*rssprk(n,nmy,5)
      k5mz = a5*rssprk(n,nmz,5)
      k5ee = a5*rssprk(n,nee,5)
      !
      k6ec = a6*rssprk(n,nec,6)
      k6mx = a6*rssprk(n,nmx,6)
      k6my = a6*rssprk(n,nmy,6)
      k6mz = a6*rssprk(n,nmz,6)
      k6ee = a6*rssprk(n,nee,6)
      !
      k7ec = a7*rssprk(n,nec,7)
      k7mx = a7*rssprk(n,nmx,7)
      k7my = a7*rssprk(n,nmy,7)
      k7mz = a7*rssprk(n,nmz,7)
      k7ee = a7*rssprk(n,nee,7)
      !
      k8ec = a8*rssprk(n,nec,8)
      k8mx = a8*rssprk(n,nmx,8)
      k8my = a8*rssprk(n,nmy,8)
      k8mz = a8*rssprk(n,nmz,8)
      k8ee = a8*rssprk(n,nee,8)
      !
      k9ec = a9*rssprk(n,nec,9)
      k9mx = a9*rssprk(n,nmx,9)
      k9my = a9*rssprk(n,nmy,9)
      k9mz = a9*rssprk(n,nmz,9)
      k9ee = a9*rssprk(n,nee,9)
      !
      k10ec = dtsp(n)*residual(n,nec)
      k10mx = dtsp(n)*residual(n,nmx)
      k10my = dtsp(n)*residual(n,nmy)
      k10mz = dtsp(n)*residual(n,nmz)
      k10ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec+k1ec+k4ec+k5ec+k6ec+k7ec+k8ec+k9ec+a10*k10ec
      usp(n,nmx) = u0mx+k1mx+k4mx+k5mx+k6mx+k7mx+k8mx+k9mx+a10*k10mx
      usp(n,nmy) = u0my+k1my+k4my+k5my+k6my+k7my+k8my+k9my+a10*k10my
      usp(n,nmz) = u0mz+k1mz+k4mz+k5mz+k6mz+k7mz+k8mz+k9mz+a10*k10mz
      usp(n,nee) = u0ee+k1ee+k4ee+k5ee+k6ee+k7ee+k8ee+k9ee+a10*k10ee
      !
      rssprk(n,nec,10) = k10ec
      rssprk(n,nmx,10) = k10mx
      rssprk(n,nmy,10) = k10my
      rssprk(n,nmz,10) = k10mz
      rssprk(n,nee,10) = k10ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = a5*rssprk(n,nec,5)
      k5mx = a5*rssprk(n,nmx,5)
      k5my = a5*rssprk(n,nmy,5)
      k5ee = a5*rssprk(n,nee,5)
      !
      k6ec = a6*rssprk(n,nec,6)
      k6mx = a6*rssprk(n,nmx,6)
      k6my = a6*rssprk(n,nmy,6)
      k6ee = a6*rssprk(n,nee,6)
      !
      k7ec = a7*rssprk(n,nec,7)
      k7mx = a7*rssprk(n,nmx,7)
      k7my = a7*rssprk(n,nmy,7)
      k7ee = a7*rssprk(n,nee,7)
      !
      k8ec = a8*rssprk(n,nec,8)
      k8mx = a8*rssprk(n,nmx,8)
      k8my = a8*rssprk(n,nmy,8)
      k8ee = a8*rssprk(n,nee,8)
      !
      k9ec = a9*rssprk(n,nec,9)
      k9mx = a9*rssprk(n,nmx,9)
      k9my = a9*rssprk(n,nmy,9)
      k9ee = a9*rssprk(n,nee,9)
      !
      k10ec = dtsp(n)*residual(n,nec)
      k10mx = dtsp(n)*residual(n,nmx)
      k10my = dtsp(n)*residual(n,nmy)
      k10ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec+k1ec+k4ec+k5ec+k6ec+k7ec+k8ec+k9ec+a10*k10ec
      usp(n,nmx) = u0mx+k1mx+k4mx+k5mx+k6mx+k7mx+k8mx+k9mx+a10*k10mx
      usp(n,nmy) = u0my+k1my+k4my+k5my+k6my+k7my+k8my+k9my+a10*k10my
      usp(n,nee) = u0ee+k1ee+k4ee+k5ee+k6ee+k7ee+k8ee+k9ee+a10*k10ee
      !
      rssprk(n,nec,10) = k10ec
      rssprk(n,nmx,10) = k10mx
      rssprk(n,nmy,10) = k10my
      rssprk(n,nee,10) = k10ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,10) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) +  a1*rssprk(n,m, 1) &
                               +  a4*rssprk(n,m, 4) &
                               +  a5*rssprk(n,m, 5) &
                               +  a6*rssprk(n,m, 6) &
                               +  a7*rssprk(n,m, 7) &
                               +  a8*rssprk(n,m, 8) &
                               +  a9*rssprk(n,m, 9) &
                               + a10*rssprk(n,m,10)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP853_S10
!
!###############################################################################
!
subroutine update_EPERK_DP853_S11(np,dtsp,residual,uoldsp,usp,rssprk, &
                                  pseudo_time)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  real(wp), intent(inout) :: pseudo_time
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k4ec,k4mx,k4my,k4mz,k4ee
  real(wp) :: k5ec,k5mx,k5my,k5mz,k5ee
  real(wp) :: k6ec,k6mx,k6my,k6mz,k6ee
  real(wp) :: k7ec,k7mx,k7my,k7mz,k7ee
  real(wp) :: k8ec,k8mx,k8my,k8mz,k8ee
  real(wp) :: k9ec,k9mx,k9my,k9mz,k9ee
  real(wp) :: k10ec,k10mx,k10my,k10mz,k10ee
  real(wp) :: k11ec,k11mx,k11my,k11mz,k11ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP853_S11"
  !
  real(wp), parameter :: a1  =  2.27331014751653820792359768449_wp
  real(wp), parameter :: a4  = -10.5344954667372501984066689879_wp
  real(wp), parameter :: a5  = -2.00087205822486249909675718444_wp
  real(wp), parameter :: a6  = -17.9589318631187989172765950534_wp
  real(wp), parameter :: a7  =  27.9488845294199600508499808837_wp
  real(wp), parameter :: a8  = -2.85899827713502369474065508674_wp
  real(wp), parameter :: a9  = -8.87285693353062954433549289258_wp
  real(wp), parameter :: a10 =  12.3605671757943030647266201528_wp
  real(wp), parameter :: a11 =  0.643392746015763530355970484046_wp
  real(wp), parameter :: c11 =  one - 0.857142857142857142857142857142_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  ! Update pseudo_time
  !
  pseudo_time = pseudo_time + c11
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1mz = a1*rssprk(n,nmz,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4mz = a4*rssprk(n,nmz,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = a5*rssprk(n,nec,5)
      k5mx = a5*rssprk(n,nmx,5)
      k5my = a5*rssprk(n,nmy,5)
      k5mz = a5*rssprk(n,nmz,5)
      k5ee = a5*rssprk(n,nee,5)
      !
      k6ec = a6*rssprk(n,nec,6)
      k6mx = a6*rssprk(n,nmx,6)
      k6my = a6*rssprk(n,nmy,6)
      k6mz = a6*rssprk(n,nmz,6)
      k6ee = a6*rssprk(n,nee,6)
      !
      k7ec = a7*rssprk(n,nec,7)
      k7mx = a7*rssprk(n,nmx,7)
      k7my = a7*rssprk(n,nmy,7)
      k7mz = a7*rssprk(n,nmz,7)
      k7ee = a7*rssprk(n,nee,7)
      !
      k8ec = a8*rssprk(n,nec,8)
      k8mx = a8*rssprk(n,nmx,8)
      k8my = a8*rssprk(n,nmy,8)
      k8mz = a8*rssprk(n,nmz,8)
      k8ee = a8*rssprk(n,nee,8)
      !
      k9ec = a9*rssprk(n,nec,9)
      k9mx = a9*rssprk(n,nmx,9)
      k9my = a9*rssprk(n,nmy,9)
      k9mz = a9*rssprk(n,nmz,9)
      k9ee = a9*rssprk(n,nee,9)
      !
      k10ec = a10*rssprk(n,nec,10)
      k10mx = a10*rssprk(n,nmx,10)
      k10my = a10*rssprk(n,nmy,10)
      k10mz = a10*rssprk(n,nmz,10)
      k10ee = a10*rssprk(n,nee,10)
      !
      k11ec = dtsp(n)*residual(n,nec)
      k11mx = dtsp(n)*residual(n,nmx)
      k11my = dtsp(n)*residual(n,nmy)
      k11mz = dtsp(n)*residual(n,nmz)
      k11ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec+k1ec+k4ec+k5ec+k6ec+k7ec+k8ec+k9ec+k10ec+a11*k11ec
      usp(n,nmx) = u0mx+k1mx+k4mx+k5mx+k6mx+k7mx+k8mx+k9mx+k10mx+a11*k11mx
      usp(n,nmy) = u0my+k1my+k4my+k5my+k6my+k7my+k8my+k9my+k10my+a11*k11my
      usp(n,nmz) = u0mz+k1mz+k4mz+k5mz+k6mz+k7mz+k8mz+k9mz+k10mz+a11*k11mz
      usp(n,nee) = u0ee+k1ee+k4ee+k5ee+k6ee+k7ee+k8ee+k9ee+k10ee+a11*k11ee
      !
      rssprk(n,nec,11) = k11ec
      rssprk(n,nmx,11) = k11mx
      rssprk(n,nmy,11) = k11my
      rssprk(n,nmz,11) = k11mz
      rssprk(n,nee,11) = k11ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k4ec = a4*rssprk(n,nec,4)
      k4mx = a4*rssprk(n,nmx,4)
      k4my = a4*rssprk(n,nmy,4)
      k4ee = a4*rssprk(n,nee,4)
      !
      k5ec = a5*rssprk(n,nec,5)
      k5mx = a5*rssprk(n,nmx,5)
      k5my = a5*rssprk(n,nmy,5)
      k5ee = a5*rssprk(n,nee,5)
      !
      k6ec = a6*rssprk(n,nec,6)
      k6mx = a6*rssprk(n,nmx,6)
      k6my = a6*rssprk(n,nmy,6)
      k6ee = a6*rssprk(n,nee,6)
      !
      k7ec = a7*rssprk(n,nec,7)
      k7mx = a7*rssprk(n,nmx,7)
      k7my = a7*rssprk(n,nmy,7)
      k7ee = a7*rssprk(n,nee,7)
      !
      k8ec = a8*rssprk(n,nec,8)
      k8mx = a8*rssprk(n,nmx,8)
      k8my = a8*rssprk(n,nmy,8)
      k8ee = a8*rssprk(n,nee,8)
      !
      k9ec = a9*rssprk(n,nec,9)
      k9mx = a9*rssprk(n,nmx,9)
      k9my = a9*rssprk(n,nmy,9)
      k9ee = a9*rssprk(n,nee,9)
      !
      k10ec = a10*rssprk(n,nec,10)
      k10mx = a10*rssprk(n,nmx,10)
      k10my = a10*rssprk(n,nmy,10)
      k10ee = a10*rssprk(n,nee,10)
      !
      k11ec = dtsp(n)*residual(n,nec)
      k11mx = dtsp(n)*residual(n,nmx)
      k11my = dtsp(n)*residual(n,nmy)
      k11ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec+k1ec+k4ec+k5ec+k6ec+k7ec+k8ec+k9ec+k10ec+a11*k11ec
      usp(n,nmx) = u0mx+k1mx+k4mx+k5mx+k6mx+k7mx+k8mx+k9mx+k10mx+a11*k11mx
      usp(n,nmy) = u0my+k1my+k4my+k5my+k6my+k7my+k8my+k9my+k10my+a11*k11my
      usp(n,nee) = u0ee+k1ee+k4ee+k5ee+k6ee+k7ee+k8ee+k9ee+k10ee+a11*k11ee
      !
      rssprk(n,nec,11) = k11ec
      rssprk(n,nmx,11) = k11mx
      rssprk(n,nmy,11) = k11my
      rssprk(n,nee,11) = k11ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,11) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) +  a1*rssprk(n,m, 1) &
                               +  a4*rssprk(n,m, 4) &
                               +  a5*rssprk(n,m, 5) &
                               +  a6*rssprk(n,m, 6) &
                               +  a7*rssprk(n,m, 7) &
                               +  a8*rssprk(n,m, 8) &
                               +  a9*rssprk(n,m, 9) &
                               + a10*rssprk(n,m,10) &
                               + a11*rssprk(n,m,11)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP853_S11
!
!###############################################################################
!
subroutine update_EPERK_DP853_S12(np,dtsp,residual,uoldsp,usp,rssprk)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rssprk(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: u0ec,u0mx,u0my,u0mz,u0ee
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k6ec,k6mx,k6my,k6mz,k6ee
  real(wp) :: k7ec,k7mx,k7my,k7mz,k7ee
  real(wp) :: k8ec,k8mx,k8my,k8mz,k8ee
  real(wp) :: k9ec,k9mx,k9my,k9mz,k9ee
  real(wp) :: k10ec,k10mx,k10my,k10mz,k10ee
  real(wp) :: k11ec,k11mx,k11my,k11mz,k11ee
  real(wp) :: k12ec,k12mx,k12my,k12mz,k12ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP853_S12"
  !
  real(wp), parameter :: a1  =  0.0542937341165687622380535766363_wp
  real(wp), parameter :: a6  =  4.45031289275240888144113950566_wp
  real(wp), parameter :: a7  =  1.89151789931450038304281599044_wp
  real(wp), parameter :: a8  = -5.8012039600105847814672114227_wp
  real(wp), parameter :: a9  =  0.311164366957819894408916062370_wp
  real(wp), parameter :: a10 = -0.152160949662516078556178806805_wp
  real(wp), parameter :: a11 =  0.201365400804030348374776537501_wp
  real(wp), parameter :: a12 =  0.0447106157277725905176885569043_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0mz = uoldsp(n,nmz)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1mz = a1*rssprk(n,nmz,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k6ec = a6*rssprk(n,nec,6)
      k6mx = a6*rssprk(n,nmx,6)
      k6my = a6*rssprk(n,nmy,6)
      k6mz = a6*rssprk(n,nmz,6)
      k6ee = a6*rssprk(n,nee,6)
      !
      k7ec = a7*rssprk(n,nec,7)
      k7mx = a7*rssprk(n,nmx,7)
      k7my = a7*rssprk(n,nmy,7)
      k7mz = a7*rssprk(n,nmz,7)
      k7ee = a7*rssprk(n,nee,7)
      !
      k8ec = a8*rssprk(n,nec,8)
      k8mx = a8*rssprk(n,nmx,8)
      k8my = a8*rssprk(n,nmy,8)
      k8mz = a8*rssprk(n,nmz,8)
      k8ee = a8*rssprk(n,nee,8)
      !
      k9ec = a9*rssprk(n,nec,9)
      k9mx = a9*rssprk(n,nmx,9)
      k9my = a9*rssprk(n,nmy,9)
      k9mz = a9*rssprk(n,nmz,9)
      k9ee = a9*rssprk(n,nee,9)
      !
      k10ec = a10*rssprk(n,nec,10)
      k10mx = a10*rssprk(n,nmx,10)
      k10my = a10*rssprk(n,nmy,10)
      k10mz = a10*rssprk(n,nmz,10)
      k10ee = a10*rssprk(n,nee,10)
      !
      k11ec = a11*rssprk(n,nec,11)
      k11mx = a11*rssprk(n,nmx,11)
      k11my = a11*rssprk(n,nmy,11)
      k11mz = a11*rssprk(n,nmz,11)
      k11ee = a11*rssprk(n,nee,11)
      !
      k12ec = dtsp(n)*residual(n,nec)
      k12mx = dtsp(n)*residual(n,nmx)
      k12my = dtsp(n)*residual(n,nmy)
      k12mz = dtsp(n)*residual(n,nmz)
      k12ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec+k1ec+k6ec+k7ec+k8ec+k9ec+k10ec+k11ec+a12*k12ec
      usp(n,nmx) = u0mx+k1mx+k6mx+k7mx+k8mx+k9mx+k10mx+k11mx+a12*k12mx
      usp(n,nmy) = u0my+k1my+k6my+k7my+k8my+k9my+k10my+k11my+a12*k12my
      usp(n,nmz) = u0mz+k1mz+k6mz+k7mz+k8mz+k9mz+k10mz+k11mz+a12*k12mz
      usp(n,nee) = u0ee+k1ee+k6ee+k7ee+k8ee+k9ee+k10ee+k11ee+a12*k12ee
      !
      rssprk(n,nec,12) = k12ec
      rssprk(n,nmx,12) = k12mx
      rssprk(n,nmy,12) = k12my
      rssprk(n,nmz,12) = k12mz
      rssprk(n,nee,12) = k12ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do n = 1,np
      !
      u0ec = uoldsp(n,nec)
      u0mx = uoldsp(n,nmx)
      u0my = uoldsp(n,nmy)
      u0ee = uoldsp(n,nee)
      !
      k1ec = a1*rssprk(n,nec,1)
      k1mx = a1*rssprk(n,nmx,1)
      k1my = a1*rssprk(n,nmy,1)
      k1ee = a1*rssprk(n,nee,1)
      !
      k6ec = a6*rssprk(n,nec,6)
      k6mx = a6*rssprk(n,nmx,6)
      k6my = a6*rssprk(n,nmy,6)
      k6ee = a6*rssprk(n,nee,6)
      !
      k7ec = a7*rssprk(n,nec,7)
      k7mx = a7*rssprk(n,nmx,7)
      k7my = a7*rssprk(n,nmy,7)
      k7ee = a7*rssprk(n,nee,7)
      !
      k8ec = a8*rssprk(n,nec,8)
      k8mx = a8*rssprk(n,nmx,8)
      k8my = a8*rssprk(n,nmy,8)
      k8ee = a8*rssprk(n,nee,8)
      !
      k9ec = a9*rssprk(n,nec,9)
      k9mx = a9*rssprk(n,nmx,9)
      k9my = a9*rssprk(n,nmy,9)
      k9ee = a9*rssprk(n,nee,9)
      !
      k10ec = a10*rssprk(n,nec,10)
      k10mx = a10*rssprk(n,nmx,10)
      k10my = a10*rssprk(n,nmy,10)
      k10ee = a10*rssprk(n,nee,10)
      !
      k11ec = a11*rssprk(n,nec,11)
      k11mx = a11*rssprk(n,nmx,11)
      k11my = a11*rssprk(n,nmy,11)
      k11ee = a11*rssprk(n,nee,11)
      !
      k12ec = dtsp(n)*residual(n,nec)
      k12mx = dtsp(n)*residual(n,nmx)
      k12my = dtsp(n)*residual(n,nmy)
      k12ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = u0ec+k1ec+k6ec+k7ec+k8ec+k9ec+k10ec+k11ec+a12*k12ec
      usp(n,nmx) = u0mx+k1mx+k6mx+k7mx+k8mx+k9mx+k10mx+k11mx+a12*k12mx
      usp(n,nmy) = u0my+k1my+k6my+k7my+k8my+k9my+k10my+k11my+a12*k12my
      usp(n,nee) = u0ee+k1ee+k6ee+k7ee+k8ee+k9ee+k10ee+k11ee+a12*k12ee
      !
      rssprk(n,nec,12) = k12ec
      rssprk(n,nmx,12) = k12mx
      rssprk(n,nmy,12) = k12my
      rssprk(n,nee,12) = k12ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,12) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) +  a1*rssprk(n,m, 1) &
                               +  a6*rssprk(n,m, 6) &
                               +  a7*rssprk(n,m, 7) &
                               +  a8*rssprk(n,m, 8) &
                               +  a9*rssprk(n,m, 9) &
                               + a10*rssprk(n,m,10) &
                               + a11*rssprk(n,m,11) &
                               + a12*rssprk(n,m,12)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP853_S12
!
!###############################################################################
!
subroutine update_EPERK_DP853_S13(np,rssprk)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp), intent(inout) :: rssprk(:,:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr,m
  real(wp) :: k1ec,k1mx,k1my,k1mz,k1ee
  real(wp) :: k6ec,k6mx,k6my,k6mz,k6ee
  real(wp) :: k7ec,k7mx,k7my,k7mz,k7ee
  real(wp) :: k8ec,k8mx,k8my,k8mz,k8ee
  real(wp) :: k9ec,k9mx,k9my,k9mz,k9ee
  real(wp) :: k10ec,k10mx,k10my,k10mz,k10ee
  real(wp) :: k11ec,k11mx,k11my,k11mz,k11ee
  real(wp) :: k12ec,k12mx,k12my,k12mz,k12ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK_DP853_S13"
  !
  real(wp), parameter :: b1  =  0.0542937341165687622380535766363_wp - &
                                0.244094488188976377952755905512_wp
  real(wp), parameter :: b6  =  4.45031289275240888144113950566_wp
  real(wp), parameter :: b7  =  1.89151789931450038304281599044_wp
  real(wp), parameter :: b8  = -5.8012039600105847814672114227_wp
  real(wp), parameter :: b9  =  0.311164366957819894408916062370_wp - &
                                0.733846688281611857341361741547_wp
  real(wp), parameter :: b10 = -0.152160949662516078556178806805_wp
  real(wp), parameter :: b11 =  0.201365400804030348374776537501_wp
  real(wp), parameter :: b12 =  0.0447106157277725905176885569043_wp - &
                                0.0220588235294117647058823529412_wp
  !
  !
  real(wp), parameter :: e1  =  0.01312004499419488073250102996_wp
  real(wp), parameter :: e6  = -1.225156446376204440720569753_wp
  real(wp), parameter :: e7  = -0.4957589496572501915214079952_wp
  real(wp), parameter :: e8  =  1.664377182454986536961530415_wp
  real(wp), parameter :: e9  = -0.3503288487499736816886487290_wp
  real(wp), parameter :: e10 =  0.3341791187130174790297318841_wp
  real(wp), parameter :: e11 =  0.08192320648511571246570742613_wp
  real(wp), parameter :: e12 = -0.02235530786388629525884427845_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(rssprk)
    do n = 1,np
      !
      k1ec = rssprk(n,nec,1)
      k1mx = rssprk(n,nmx,1)
      k1my = rssprk(n,nmy,1)
      k1mz = rssprk(n,nmz,1)
      k1ee = rssprk(n,nee,1)
      !
      k6ec = rssprk(n,nec,6)
      k6mx = rssprk(n,nmx,6)
      k6my = rssprk(n,nmy,6)
      k6mz = rssprk(n,nmz,6)
      k6ee = rssprk(n,nee,6)
      !
      k7ec = rssprk(n,nec,7)
      k7mx = rssprk(n,nmx,7)
      k7my = rssprk(n,nmy,7)
      k7mz = rssprk(n,nmz,7)
      k7ee = rssprk(n,nee,7)
      !
      k8ec = rssprk(n,nec,8)
      k8mx = rssprk(n,nmx,8)
      k8my = rssprk(n,nmy,8)
      k8mz = rssprk(n,nmz,8)
      k8ee = rssprk(n,nee,8)
      !
      k9ec = rssprk(n,nec,9)
      k9mx = rssprk(n,nmx,9)
      k9my = rssprk(n,nmy,9)
      k9mz = rssprk(n,nmz,9)
      k9ee = rssprk(n,nee,9)
      !
      k10ec = rssprk(n,nec,10)
      k10mx = rssprk(n,nmx,10)
      k10my = rssprk(n,nmy,10)
      k10mz = rssprk(n,nmz,10)
      k10ee = rssprk(n,nee,10)
      !
      k11ec = rssprk(n,nec,11)
      k11mx = rssprk(n,nmx,11)
      k11my = rssprk(n,nmy,11)
      k11mz = rssprk(n,nmz,11)
      k11ee = rssprk(n,nee,11)
      !
      k12ec = rssprk(n,nec,12)
      k12mx = rssprk(n,nmx,12)
      k12my = rssprk(n,nmy,12)
      k12mz = rssprk(n,nmz,12)
      k12ee = rssprk(n,nee,12)
      !
      ! Use rssprk(:,:,2) and rssprk(:,:,3) to store these values
      ! since we dont need those residuals any more
      !
      rssprk(n,nec,2) = b1*k1ec + b6*k6ec   + b7*k7ec   + b8*k8ec + &
                        b9*k9ec + b10*k10ec + b11*k11ec + b12*k12ec
      rssprk(n,nmx,2) = b1*k1mx + b6*k6mx   + b7*k7mx   + b8*k8mx + &
                        b9*k9mx + b10*k10mx + b11*k11mx + b12*k12mx
      rssprk(n,nmy,2) = b1*k1my + b6*k6my   + b7*k7my   + b8*k8my + &
                        b9*k9my + b10*k10my + b11*k11my + b12*k12my
      rssprk(n,nmz,2) = b1*k1mz + b6*k6mz   + b7*k7mz   + b8*k8mz + &
                        b9*k9mz + b10*k10mz + b11*k11mz + b12*k12mz
      rssprk(n,nee,2) = b1*k1ee + b6*k6ee   + b7*k7ee   + b8*k8ee + &
                        b9*k9ee + b10*k10ee + b11*k11ee + b12*k12ee
      !
      rssprk(n,nec,3) = e1*k1ec + e6*k6ec   + e7*k7ec   + e8*k8ec + &
                        e9*k9ec + e10*k10ec + e11*k11ec + e12*k12ec
      rssprk(n,nmx,3) = e1*k1mx + e6*k6mx   + e7*k7mx   + e8*k8mx + &
                        e9*k9mx + e10*k10mx + e11*k11mx + e12*k12mx
      rssprk(n,nmy,3) = e1*k1my + e6*k6my   + e7*k7my   + e8*k8my + &
                        e9*k9my + e10*k10my + e11*k11my + e12*k12my
      rssprk(n,nmz,3) = e1*k1mz + e6*k6mz   + e7*k7mz   + e8*k8mz + &
                        e9*k9mz + e10*k10mz + e11*k11mz + e12*k12mz
      rssprk(n,nee,3) = e1*k1ee + e6*k6ee   + e7*k7ee   + e8*k8ee + &
                        e9*k9ee + e10*k10ee + e11*k11ee + e12*k12ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not. force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(rssprk)
    do n = 1,np
      !
      k1ec = rssprk(n,nec,1)
      k1mx = rssprk(n,nmx,1)
      k1my = rssprk(n,nmy,1)
      k1ee = rssprk(n,nee,1)
      !
      k6ec = rssprk(n,nec,6)
      k6mx = rssprk(n,nmx,6)
      k6my = rssprk(n,nmy,6)
      k6ee = rssprk(n,nee,6)
      !
      k7ec = rssprk(n,nec,7)
      k7mx = rssprk(n,nmx,7)
      k7my = rssprk(n,nmy,7)
      k7ee = rssprk(n,nee,7)
      !
      k8ec = rssprk(n,nec,8)
      k8mx = rssprk(n,nmx,8)
      k8my = rssprk(n,nmy,8)
      k8ee = rssprk(n,nee,8)
      !
      k9ec = rssprk(n,nec,9)
      k9mx = rssprk(n,nmx,9)
      k9my = rssprk(n,nmy,9)
      k9ee = rssprk(n,nee,9)
      !
      k10ec = rssprk(n,nec,10)
      k10mx = rssprk(n,nmx,10)
      k10my = rssprk(n,nmy,10)
      k10ee = rssprk(n,nee,10)
      !
      k11ec = rssprk(n,nec,11)
      k11mx = rssprk(n,nmx,11)
      k11my = rssprk(n,nmy,11)
      k11ee = rssprk(n,nee,11)
      !
      k12ec = rssprk(n,nec,12)
      k12mx = rssprk(n,nmx,12)
      k12my = rssprk(n,nmy,12)
      k12ee = rssprk(n,nee,12)
      !
      ! Use rssprk(:,:,2) and rssprk(:,:,3) to store these values
      ! since we dont need those residuals any more
      !
      rssprk(n,nec,2) = b1*k1ec + b6*k6ec   + b7*k7ec   + b8*k8ec + &
                        b9*k9ec + b10*k10ec + b11*k11ec + b12*k12ec
      rssprk(n,nmx,2) = b1*k1mx + b6*k6mx   + b7*k7mx   + b8*k8mx + &
                        b9*k9mx + b10*k10mx + b11*k11mx + b12*k12mx
      rssprk(n,nmy,2) = b1*k1my + b6*k6my   + b7*k7my   + b8*k8my + &
                        b9*k9my + b10*k10my + b11*k11my + b12*k12my
      rssprk(n,nee,2) = b1*k1ee + b6*k6ee   + b7*k7ee   + b8*k8ee + &
                        b9*k9ee + b10*k10ee + b11*k11ee + b12*k12ee
      !
      rssprk(n,nec,3) = e1*k1ec + e6*k6ec   + e7*k7ec   + e8*k8ec + &
                        e9*k9ec + e10*k10ec + e11*k11ec + e12*k12ec
      rssprk(n,nmx,3) = e1*k1mx + e6*k6mx   + e7*k7mx   + e8*k8mx + &
                        e9*k9mx + e10*k10mx + e11*k11mx + e12*k12mx
      rssprk(n,nmy,3) = e1*k1my + e6*k6my   + e7*k7my   + e8*k8my + &
                        e9*k9my + e10*k10my + e11*k11my + e12*k12my
      rssprk(n,nee,3) = e1*k1ee + e6*k6ee   + e7*k7ee   + e8*k8ee + &
                        e9*k9ee + e10*k10ee + e11*k11ee + e12*k12ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(rssprk)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,2) =  b1*rssprk(n,m, 1) &
                      +  b6*rssprk(n,m, 6) &
                      +  b7*rssprk(n,m, 7) &
                      +  b8*rssprk(n,m, 8) &
                      +  b9*rssprk(n,m, 9) &
                      + b10*rssprk(n,m,10) &
                      + b11*rssprk(n,m,11) &
                      + b12*rssprk(n,m,12)
        !
        rssprk(n,m,3) =  e1*rssprk(n,m, 1) &
                      +  e6*rssprk(n,m, 6) &
                      +  e7*rssprk(n,m, 7) &
                      +  e8*rssprk(n,m, 8) &
                      +  e9*rssprk(n,m, 9) &
                      + e10*rssprk(n,m,10) &
                      + e11*rssprk(n,m,11) &
                      + e12*rssprk(n,m,12)
        !
      end do
    end do
  end if
  !
end subroutine update_EPERK_DP853_S13
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine init_eperk_coefficients()
  !
  !.. Use Statements ..
  use ovar, only : Runge_Kutta_Scheme,num_rk_stages
  !
  !#############################################################################
  !#############################################################################
  ! 5-Stage/3rd-order with embedded 2nd-order
  ! NOTE: Some of these coefficients come from Table B1 on page 31
  !       but all of them can be found in the file
  !       https://github.com/ranocha/Optimized-RK-CFD/blob/main/3Sstarp35.txt
  !#############################################################################
  !#############################################################################
  !
  ! NOTE: This is the diagonal component of the beta matrix from the file. The
  !       algorithm in the paper refers to using the beta_i coefficients
  !       provided in the paper, but these coefficients are being used instead
  !       for some reason that I cant figure out or remember why. Some quick
  !       tests show that this beta array works better than the b array below,
  !       so this may be correct.
  real(wp), parameter :: OEPERK_5_3_Coef_beta(1:5) = &
                              [ 0.2300285062878154351930669430512780706_wp, &
                                0.3021457892454169700189445968126242994_wp, &
                                0.8025601039472704213300183888573974531_wp, &
                                0.4362158997637629844305216319994356355_wp, &
                                0.1129268494470295369172265188216779157_wp]
  !
  ! NOTE: This is beta_i from the table, but vector b from the file
  real(wp), parameter :: OEPERK_5_3_Coef_b(1:6) = &
                              [ 0.1147931563369900682037379182772608287_wp, &
                                0.08933559295232859013880114997436974196_wp, &
                                0.4355858717379231779899161991033964256_wp, &
                                0.2473585295257286267503182138232950881_wp, &
                                0.1129268494470295369172265188216779157_wp, &
                                zero]
  !
  real(wp), parameter :: OEPERK_5_3_Coef_bhat(1:6) = &
                              [ 0.1046363371354093758897668305991705199_wp, &
                                0.09520431574956758809511173383346476348_wp, &
                                0.4482446645568668405072421350300379357_wp, &
                                0.2449030295461310135957132640369862245_wp, &
                                0.1070116530120251819121660365003405564_wp, &
                                zero]
  !
  ! NOTE: This is only found in the file, and can be computed by summing each
  !       row of the A matrix of the Butcher coefficients.
  real(wp), parameter :: OEPERK_5_3_Coef_c(1:5) = &
                              [ zero, &
                                0.2300285062878154351930669430512780706_wp, &
                                0.4050049049262914975700372321130661410_wp, &
                                0.8947823877926760224705450466361360720_wp, &
                                0.7235108137218888081489570284485201518_wp]
  !
  real(wp), parameter :: OEPERK_5_3_Coef_g1(1:5) = &
                              [ zero, &
                                0.2587669070352079020144955303389306026_wp, &
                               -0.1324366873994502973977035353758550057_wp, &
                                0.05055601231460399101814291350373559483_wp, &
                                0.5670552807902877312521811889846000976_wp]
  !
  real(wp), parameter :: OEPERK_5_3_Coef_g2(1:5) = &
                              [ one, &
                                0.5528418745102160639901976698795928733_wp, &
                                0.6731844400389673824374042790213570079_wp, &
                                0.2803103804507635075215805236096803381_wp, &
                                0.5521508873507393276457754945308880998_wp]
  !
  real(wp), parameter :: OEPERK_5_3_Coef_g3(1:5) = &
                              [ zero, &
                                zero, &
                                zero, &
                                0.2752585813446636957256614568573008811_wp, &
                               -0.8950548709279785077579454232514633376_wp]
  !
  real(wp), parameter :: OEPERK_5_3_Coef_d(1:5) = &
                              [ one, &
                                0.3407687209321455242558804921815861422_wp, &
                                0.3414399280584625023244387687873774697_wp, &
                                0.7229302732875589702087936723400941329_wp, &
                                zero]
  !
  !#############################################################################
  !#############################################################################
  ! FSAL 5-Stage/3rd-order with embedded 2nd-order
  ! NOTE: Some of these coefficients come from Table B2 on page 31
  !       but all of them can be found in the file
  !      https://github.com/ranocha/Optimized-RK-CFD/blob/main/3SstarpFSAL35.txt
  !#############################################################################
  !#############################################################################
  !
  ! NOTE: This is the diagonal component of the beta matrix from the file. The
  !       algorithm in the paper refers to using the beta_i coefficients
  !       provided in the paper, but these coefficients are being used instead
  !       for some reason that I cant figure out or remember why. Some quick
  !       tests show that this beta array works better than the b array below,
  !       so this may be correct.
  real(wp), parameter :: OEPERK_5_3_FSAL_Coef_beta(1:5) = &
                              [ 0.2300298624518076223899418286314123354_wp, &
                                0.3021434166948288809034402119555380003_wp, &
                                0.8025606185416310937583009085873554681_wp, &
                                0.4362158943603440930655148245148766471_wp, &
                                0.1129272530455059129782111662594436580_wp]
  !
  ! NOTE: This is beta_i from the table, but vector b from the file
  real(wp), parameter :: OEPERK_5_3_FSAL_Coef_b(1:6) = &
                              [ 0.1147935971023541171733601324486904546_wp, &
                                0.08933442853113315592708384523126474636_wp, &
                                0.4355871025008616992483722693795608738_wp, &
                                0.2473576188201451146729725866810402672_wp, &
                                0.1129272530455059129782111662594436580_wp, &
                                zero]
  !
  real(wp), parameter :: OEPERK_5_3_FSAL_Coef_bhat(1:6) = &
                              [ 0.09484166705035703392326247283838082847_wp, &
                                0.1726371339430353766966762629176676070_wp, &
                                0.3998243189084371024483169698618455770_wp, &
                                0.1718016807580178450618829007973835152_wp, &
                                0.05881914422155740300718268359027168467_wp, &
                                0.1020760551185952388626787099944507877_wp]
  !
  ! NOTE: This is only found in the file, and can be computed by summing each
  !       row of the A matrix of the Butcher coefficients.
  real(wp), parameter :: OEPERK_5_3_FSAL_Coef_c(1:5) = &
                              [ zero, &
                                0.2300298624518076223899418286314123354_wp, &
                                0.4050046072094990912268498160116125481_wp, &
                                0.8947822893693433545220710894560512805_wp, &
                                0.7235136928826589010272834603680114769_wp]
  !
  real(wp), parameter :: OEPERK_5_3_FSAL_Coef_g1(1:5) = &
                              [ zero, &
                                0.2587771979725733308135192812685323706_wp, &
                               -0.1324380360140723382965420909764953437_wp, &
                                0.05056033948190826045833606441415585735_wp, &
                                0.5670532000739313812633197158607642990_wp]
  !
  real(wp), parameter :: OEPERK_5_3_FSAL_Coef_g2(1:5) = &
                              [ one, &
                                0.5528354909301389892439698870483746541_wp, &
                                0.6731871608203061824849561782794643600_wp, &
                                0.2803103963297672407841316576323901761_wp, &
                                0.5521525447020610386070346724931300367_wp]
  !
  real(wp), parameter :: OEPERK_5_3_FSAL_Coef_g3(1:5) = &
                              [ zero, &
                                zero, &
                                zero, &
                                0.2752563273304676380891217287572780582_wp, &
                               -0.8950526174674033822276061734289327568_wp]
  !
  real(wp), parameter :: OEPERK_5_3_FSAL_Coef_d(1:5) = &
                              [ one, &
                                0.3407655879334525365094815965895763636_wp, &
                                0.3414382655003386206551709871126405331_wp, &
                                0.7229275366787987419692007421895451953_wp, &
                                zero]
 !!
 !real(wp), parameter :: OERK_VO_5_3_Coef_c(1:5) =  &
 !                               [ zero, &
 !                                 0.37469540309464722_wp, &
 !                                 0.47077676669889185_wp, &
 !                                 0.97644725493183659_wp, &
 !                                 0.87244510145462129_wp]
 !real(wp), parameter :: OERK_VO_5_3_Coef_b(1:5) =  &
 !                               [ 0.37469540309464722_wp, &
 !                                 0.43191672433390471_wp, &
 !                                 0.70039384794436799_wp, &
 !                                 0.12238246651043405_wp, &
 !                                 0.20492309083147120_wp]
 !real(wp), parameter :: OERK_VO_5_3_Coef_bhat(1:6) =  &
 !                               [ 0.37469540309464722_wp, &
 !                                 0.43191672433390471_wp, &
 !                                 0.70039384794436799_wp, &
 !                                 0.12238246651043405_wp, &
 !                                 0.20492309083147120_wp, &
 !                                 zero]
 !real(wp), parameter :: OERK_VO_5_3_Coef_g1(1:5) =  &
 !                               [ zero, &
 !                                -0.40733194644117998_wp, &
 !                                 0.11227306818680416_wp, &
 !                                 0.41979779027566899_wp, &
 !                                 0.56736097340799185_wp]
 !real(wp), parameter :: OERK_VO_5_3_Coef_g2(1:5) =  &
 !                               [ one, &
 !                                 0.89628898021155823_wp, &
 !                                 0.37050635458955933_wp, &
 !                                 0.25928503140256132_wp, &
 !                                 0.22874377431566192_wp]
 !real(wp), parameter :: OERK_VO_5_3_Coef_g3(1:5) =  &
 !                               [ zero, &
 !                                 zero, &
 !                                 zero, &
 !                                -0.22943304087221603_wp, &
 !                                -0.28162905914613162_wp]
 !real(wp), parameter :: OERK_VO_5_3_Coef_d(1:5) =  &
 !                               [ one, &
 !                                 0.57017655857935046_wp, &
 !                                 0.82580645444239575_wp, &
 !                                 0.72658540605340605_wp, &
 !                                 zero]
  !
  !#############################################################################
  !#############################################################################
  ! 9-Stage/4th-order with embedded 3rd-order
  ! NOTE: Some of these coefficients come from Table B3 on page 32
  !       but all of them can be found in the file
  !       https://github.com/ranocha/Optimized-RK-CFD/blob/main/3Sstarp49.txt
  !#############################################################################
  !#############################################################################
  !
  ! NOTE: This is the diagonal component of the beta matrix from the file. The
  !       algorithm in the paper refers to using the beta_i coefficients
  !       provided in the paper, but these coefficients are being used instead
  !       for some reason that I cant figure out or remember why. Some quick
  !       tests show that this beta array works better than the b array below,
  !       so this may be correct.
  real(wp), parameter :: OEPERK_9_4_Coef_beta(1:9) = &
                              [ 0.2836343531977826022543660465926414772_wp, &
                                0.9736497978646965372894268287659773644_wp, &
                                0.3382358566377620380505126936670933370_wp, &
                               -0.3584937820217850715182820651063453804_wp, &
                               -0.004113955814725134294322006403954822487_wp, &
                                1.427968962196019024010757034274849198_wp, &
                                0.01808467712038743032991177525728915926_wp, &
                                0.1605771316794521018947553625079465692_wp, &
                                0.2952226811394310028003810072027839487_wp]
  !
  ! NOTE: This is beta_i from the table, but vector b from the file
  real(wp), parameter :: OEPERK_9_4_Coef_b(1:10) = &
                              [ 0.04503731969165884304041981629148469971_wp, &
                                0.1859217322011968812563859888433403777_wp, &
                                0.03329727509207630932171676116314110008_wp, &
                               -0.004784222621050198909820741390895649698_wp, &
                                0.004055848062637567925908043629915811671_wp, &
                                0.4185027999682794463309031355073933444_wp, &
                               -0.004381894507474277848407591859322000026_wp, &
                                0.02712846097324442608251358061215836749_wp, &
                                0.2952226811394310028003810072027839487_wp, &
                                zero]
  !
  real(wp), parameter :: OEPERK_9_4_Coef_bhat(1:10) = &
                              [ 0.04550655927970944948340364817140593012_wp, &
                                0.1175968310492638562142460384341959193_wp, &
                                0.03658257330515213200375475084421083608_wp, &
                               -0.005311555834355629559010061596928357525_wp, &
                                0.005178250012713127329531367677410650996_wp, &
                                0.4954639022118682638697706200022961443_wp, &
                               -0.005999303132737865921441409466809521699_wp, &
                                0.09405093434568315929035250835218733824_wp, &
                                0.2169318087627035072893925375820310602_wp, &
                                zero]
  !
  ! NOTE: This is only found in the file, and can be computed by summing each
  !       row of the A matrix of the Butcher coefficients.
  real(wp), parameter :: OEPERK_9_4_Coef_c(1:9) = &
                              [ zero, &
                                0.2836343531977826022543660465926414772_wp, &
                                0.5484073767552486705240014599676811834_wp, &
                                0.3687229456675706936558667052479014150_wp, &
                               -0.6806119916032093175251948474173648331_wp, &
                                0.3518526451892056368706593492732753284_wp, &
                                1.665941920204672094647868254892387293_wp, &
                                0.9715276989307335935187466054546761665_wp, &
                                0.9051569554420045339601721625247585643_wp]
  !
  real(wp), parameter :: OEPERK_9_4_Coef_g1(1:9) = &
                              [ zero, &
                               -4.655641301259180308677051498071354582_wp, &
                               -0.7720264924836063859141482018013692338_wp, &
                               -4.024423213419724605695005429153112050_wp, &
                               -0.02129685246739018613087466942802498152_wp, &
                               -2.435022519234470128602335652131234586_wp, &
                                0.01985627480986167686791439120784668251_wp, &
                               -0.2810790112885283952929218377438668784_wp, &
                                0.1689434895835535695524003319503844110_wp]
  !
  real(wp), parameter :: OEPERK_9_4_Coef_g2(1:9) = &
                              [ one, &
                                2.499262752607825957145627300817258023_wp, &
                                0.5866820365436136799319929406678132638_wp, &
                                1.205141365412670762568835277881144391_wp, &
                                0.3474793796700868848597960521248007941_wp, &
                                1.321346140128723105871355808477092220_wp, &
                                0.3119636324379370564023292317172847140_wp, &
                                0.4351419055894087609560896967082486864_wp, &
                                0.2359698299440788299161958168555704234_wp]
  !
  real(wp), parameter :: OEPERK_9_4_Coef_g3(1:9) = &
                              [ zero, &
                                zero, &
                                zero, &
                                0.7621037111138170045618771082985664430_wp, &
                               -0.1981182159087218433914909510116664154_wp, &
                               -0.6228960706317566993192689455719570179_wp, &
                               -0.3752246993432626328289874575355102038_wp, &
                               -0.3355436539000946543242869676125143358_wp, &
                               -0.04560963110717484359015342341157302403_wp]
  !
  real(wp), parameter :: OEPERK_9_4_Coef_d(1:9) = &
                              [ one, &
                                1.262923854387806460989545005598562667_wp, &
                                0.7574967177560872438940839460448329992_wp, &
                                0.5163591158111222863455531895152351544_wp, &
                               -0.02746333792042827389548936599648122146_wp, &
                               -0.4382674653941770848797864513655752318_wp, &
                                1.273587103668392811985704533534301656_wp, &
                               -0.6294740045442794829622796613103492913_wp, &
                                zero]
  !
  !#############################################################################
  !#############################################################################
  ! FSAL 9-Stage/4th-order with embedded 3rd-order
  ! NOTE: Some of these coefficients come from Table B4 on page 32
  !       but all of them can be found in the file
  !      https://github.com/ranocha/Optimized-RK-CFD/blob/main/3SstarpFSAL49.txt
  !#############################################################################
  !#############################################################################
  !
  ! NOTE: This is the diagonal component of the beta matrix from the file. The
  !       algorithm in the paper refers to using the beta_i coefficients
  !       provided in the paper, but these coefficients are being used instead
  !       for some reason that I cant figure out or remember why. Some quick
  !       tests show that this beta array works better than the b array below,
  !       so this may be correct.
  real(wp), parameter :: OEPERK_9_4_FSAL_Coef_beta(1:9) = &
                              [ 0.2836343005184365275160654678626695428_wp, &
                                0.9736500104654741223716056170419660217_wp, &
                                0.3382359225242515288768487569778320563_wp, &
                               -0.3584943611106183357043212309791897386_wp, &
                               -0.004113944068471528211627210454497620358_wp, &
                                1.427968894048586363415504654313371031_wp, &
                                0.01808470948394314017665968411915568633_wp, &
                                0.1605770645946802213926893453819236685_wp, &
                                0.2952227015964591648775833803635147962_wp]
  !
  ! NOTE: This is beta_i from the table, but vector b from the file
  real(wp), parameter :: OEPERK_9_4_FSAL_Coef_b(1:10) = &
                              [ 0.04503732627263753698356970706617404465_wp, &
                                0.1859217303699847950262276860012454333_wp, &
                                0.03329729672569717599759560403851202805_wp, &
                               -0.004784204180958975587114459316829942677_wp, &
                                0.004055835961031310727671557609188874328_wp, &
                                0.4185027772596074197662616795629003544_wp, &
                               -0.004381901968919326084347037216500072323_wp, &
                                0.02712843796446089829255188189179448399_wp, &
                                0.2952227015964591648775833803635147962_wp, &
                                zero]
  !
  real(wp), parameter :: OEPERK_9_4_FSAL_Coef_bhat(1:10) = &
                              [ 0.02483675912451591196775756814283216443_wp, &
                                0.1866327774562103796990092260942180726_wp, &
                                0.05671080795936984495604436622517631183_wp, &
                               -0.003447695439149287702616943808570747099_wp, &
                                0.003602245056516636472203469198006404016_wp, &
                                0.4545570622145088936800484247980581766_wp, &
                               -0.0002434665289427612407531544765622888855_wp, &
                                0.06642755361103549971517945063138312147_wp, &
                                0.1613697079523505006226025497715177578_wp, &
                                0.04955424859358438183052504342394102722_wp]
  !
  ! NOTE: This is only found in the file, and can be computed by summing each
  !       row of the A matrix of the Butcher coefficients.
  real(wp), parameter :: OEPERK_9_4_FSAL_Coef_c(1:9) = &
                              [ zero, &
                                0.2836343005184365275160654678626695428_wp, &
                                0.5484076570002894365286665352032296535_wp, &
                                0.3687228761669438493478872632332010073_wp, &
                               -0.6806126440140844191258463830024463902_wp, &
                                0.3518526124230705801739919476290327750_wp, &
                                1.665941994879593315477304663913129942_wp, &
                                0.9715279295934715835299192116436237065_wp, &
                                0.9051569840159589594903399929316959062_wp]
  !
  real(wp), parameter :: OEPERK_9_4_FSAL_Coef_g1(1:9) = &
                              [ zero, &
                               -4.655641447335068552684422206224169103_wp, &
                               -0.7720265099645871829248487209517314217_wp, &
                               -4.024436690519806086742256154738379161_wp, &
                               -0.02129676284018530966221583708648634733_wp, &
                               -2.435022509790109546199372365866450709_wp, &
                                0.01985627297131987000579523283542615256_wp, &
                               -0.2810791146791038566946663374735713961_wp, &
                                0.1689434168754859644351230590422137972_wp]
  !
  real(wp), parameter :: OEPERK_9_4_FSAL_Coef_g2(1:9) = &
                              [ one, &
                                2.499262792574495009336242992898153462_wp, &
                                0.5866820377718875577451517985847920081_wp, &
                                1.205146086523094569925592464380295241_wp, &
                                0.3474793722186732780030762737753849272_wp, &
                                1.321346060965113109321230804210670518_wp, &
                                0.3119636464694193615946633676950358444_wp, &
                                0.4351419539684379261368971206040518552_wp, &
                                0.2359698130028753572503744518147537768_wp]
  !
  real(wp), parameter :: OEPERK_9_4_FSAL_Coef_g3(1:9) = &
                              [ zero, &
                                zero, &
                                zero, &
                                0.7621006678721315291614677352949377871_wp, &
                               -0.1981182504339400567765766904309673119_wp, &
                               -0.6228959218699007450469629366684127462_wp, &
                               -0.3752248380775956442989480369774937099_wp, &
                               -0.3355438309135169811915662336248989661_wp, &
                               -0.04560955005031121479972862973705108039_wp]
  !
  real(wp), parameter :: OEPERK_9_4_FSAL_Coef_d(1:9) = &
                              [ one, &
                                1.262923876648114432874834923838556100_wp, &
                                0.7574967189685911558308119415539596711_wp, &
                                0.5163589453140728104667573195005629833_wp, &
                               -0.02746327421802609557034437892013640319_wp, &
                               -0.4382673178127944142238606608356542890_wp, &
                                1.273587294602656522645691372699677063_wp, &
                               -0.6294740283927400326554066998751383342_wp, &
                                zero]
  !
  !#############################################################################
  !#############################################################################
  ! FSAL 10-Stage/5th-order with embedded 4th-order
  ! NOTE: Some of these coefficients come from Table B5 on page 33
  !       but all of them can be found in the file
  !       https://github.com/ranocha/Optimized-RK-CFD/blob/main/3Sstarp510.txt
  !#############################################################################
  !#############################################################################
  !
  ! NOTE: This is the diagonal component of the beta matrix from the file. The
  !       algorithm in the paper refers to using the beta_i coefficients
  !       provided in the paper, but these coefficients are being used instead
  !       for some reason that I cant figure out or remember why. Some quick
  !       tests show that this beta array works better than the b array below,
  !       so this may be correct.
  real(wp), parameter :: OEPERK_10_5_Coef_beta(1:10) = &
                              [ 0.2597883575710995826783320802193635406_wp, &
                                0.01777008800169541694837687556103565007_wp, &
                                0.2481636637328140606807905234325691851_wp, &
                                0.7941736827560429420202759490815682546_wp, &
                                0.3885391296871822541486945325814526190_wp, &
                                0.1455051664264339366757555740296587660_wp, &
                                0.1587517379462528932413419955691782412_wp, &
                                0.1650605631567659573994022720500446501_wp, &
                                0.2118093299943235065178000892467421832_wp, &
                                0.1559392340339606299335442956580114440_wp]
  !
  ! NOTE: This is beta_i from the table, but vector b from the file
  real(wp), parameter :: OEPERK_10_5_Coef_b(1:11) = &
                              [-0.002280102305596364773323878383881954511_wp, &
                                0.01407393020823230537861040991952849386_wp, &
                                0.2332691794172822486743039657924919496_wp, &
                                0.04808266700465181307162297999657715930_wp, &
                                0.4119003221139622842134291677033040683_wp, &
                               -0.1291461071364752805327361051196128312_wp, &
                                0.1220746011038579789984601943748468541_wp, &
                                0.04357858803113387764356338334851554715_wp, &
                                0.1025076875289905073925255867102192694_wp, &
                                0.1559392340339606299335442956580114440_wp, &
                                zero]
  !
  real(wp), parameter :: OEPERK_10_5_Coef_bhat(1:11) = &
                              [ 0.05734588484676193812418453938089759359_wp, &
                                0.01971447518039733870541652912891291496_wp, &
                                0.07215296605683716720707226840456658773_wp, &
                                0.1739659489807939956977075317768151880_wp, &
                                0.3703693600445487815015171515640585668_wp, &
                               -0.1215599039055065009827765147821222534_wp, &
                                0.1180372945491121604465067725859678821_wp, &
                                0.04155688823364870056536983972605056553_wp, &
                                0.1227886627910379901351569893551486490_wp, &
                                0.1456284232223684285998448928597043056_wp, &
                                zero]
  !
  ! NOTE: This is only found in the file, and can be computed by summing each
  !       row of the A matrix of the Butcher coefficients.
  real(wp), parameter :: OEPERK_10_5_Coef_c(1:10) = &
                              [ zero, &
                                0.2597883575710995826783320802193635406_wp, &
                                0.09904573115730917688557891428202061598_wp, &
                                0.2155511882303785204133426661931565216_wp, &
                                0.5007950078421880417512789524851012021_wp, &
                                0.5592251914858131230054392022144328176_wp, &
                                0.5449986973408778242805929551952000165_wp, &
                                0.7615224662599497796472095353126697300_wp, &
                                0.8427062083059167761623893618875787414_wp, &
                                0.9152209807185253394871325258038753352_wp]
  !
  real(wp), parameter :: OEPERK_10_5_Coef_g1(1:10) = &
                              [ zero, &
                                0.4043660078504695837542588769963326988_wp, &
                               -0.8503427464263185087039788184485627962_wp, &
                               -6.950894167072419998080989313353063399_wp, &
                                0.9238765225328278557805080247596562995_wp, &
                               -2.563178039957404359875124580586147888_wp, &
                                0.2545744869966347362604059848503340890_wp, &
                                0.3125831733863168874151935287174374515_wp, &
                               -0.7007114800567584871263283872289072079_wp, &
                                0.4839620970980726631935174740648996010_wp]
  !
  real(wp), parameter :: OEPERK_10_5_Coef_g2(1:10) = &
                              [ one, &
                                0.6871467069752345566001768382316915820_wp, &
                                1.093024760468898686510433898645775908_wp, &
                                3.225975382330161123625348062949430509_wp, &
                                1.041153700841396427100436517666787823_wp, &
                                1.292821488864702752767390075072674807_wp, &
                                0.7391462769297006312785029455392854586_wp, &
                                0.1239129257039300081860496157739352186_wp, &
                                0.1842753479366766790220633908793933781_wp, &
                                0.05712788942697077644959290025755003720_wp]
  !
  real(wp), parameter :: OEPERK_10_5_Coef_g3(1:10) = &
                              [ zero, &
                                zero, &
                                zero, &
                               -2.393405159342139386425044844626597490_wp, &
                               -1.902854422095986544338294743445530533_wp, &
                               -2.820042210583207174321941694153843259_wp, &
                               -1.832698464130564949123807896975136336_wp, &
                               -0.2199094510750697865007677774395365522_wp, &
                               -0.4082430660384876496971887725512427800_wp, &
                               -0.1377669791121207993339861855818881150_wp]
  !
  real(wp), parameter :: OEPERK_10_5_Coef_d(1:10) = &
                              [ one, &
                               -0.1331778409133849616712007380176762548_wp, &
                                0.8260422785246030254485064732649153253_wp, &
                                1.513700430513332405798616943654007796_wp, &
                               -1.305810063177048110528482211982726539_wp, &
                                3.036678789342507704281817524408221954_wp, &
                               -1.449458267074592489788800461540171106_wp, &
                                3.834313873320957483471400258279635203_wp, &
                                4.122293971923324492772059928094971199_wp, &
                                zero]
  !
  !#############################################################################
  !#############################################################################
  ! FSAL 10-Stage/5th-order with embedded 4th-order
  ! NOTE: Some of these coefficients come from Table B6 on page 34
  !       but all of them can be found in the file
  !     https://github.com/ranocha/Optimized-RK-CFD/blob/main/3SstarpFSAL510.txt
  !#############################################################################
  !#############################################################################
  !
  ! NOTE: This is the diagonal component of the beta matrix from the file. The
  !       algorithm in the paper refers to using the beta_i coefficients
  !       provided in the paper, but these coefficients are being used instead
  !       for some reason that I cant figure out or remember why. Some quick
  !       tests show that this beta array works better than the b array below,
  !       so this may be correct.
  real(wp), parameter :: OEPERK_10_5_FSAL_Coef_beta(1:10) = &
                              [ 0.2597883554788674084039539165398464630_wp, &
                                0.01777008889438867858759149597539211023_wp, &
                                0.2481636629715501931294746189266601496_wp, &
                                0.7941736871152005775821844297293296135_wp, &
                                0.3885391285642019129575902994397298066_wp, &
                                0.1455051657916305055730603387469193768_wp, &
                                0.1587517385964749337690916959584348979_wp, &
                                0.1650605617880053419242434594242509601_wp, &
                                0.2118093284937153836908655490906875007_wp, &
                                0.1559392342362059886106995325687547506_wp]
  !
  ! NOTE: This is beta_i from the table, but vector b from the file
  real(wp), parameter :: OEPERK_10_5_FSAL_Coef_b(1:11) = &
                              [-0.002280100321836980811830528665041532799_wp, &
                                0.01407393115790186300730580636032878435_wp, &
                                0.2332691775508456597719992034291118324_wp, &
                                0.04808266741353862546318531020856621860_wp, &
                                0.4119003217706951892385733111000873172_wp, &
                               -0.1291461067807736321056740833501596735_wp, &
                                0.1220746013848710098878384114422516148_wp, &
                                0.04357858583174420432201228508067333299_wp, &
                                0.1025076877568080726158907518254273554_wp, &
                                0.1559392342362059886106995325687547506_wp, &
                                zero]
  !
  real(wp), parameter :: OEPERK_10_5_FSAL_Coef_bhat(1:11) = &
                              [-0.02019255440012066080909442770590267512_wp, &
                                0.02737903480959184339932730854141598275_wp, &
                                0.3028818636145965534365173822296811090_wp, &
                               -0.03656843880622222190071445247906780540_wp, &
                                0.3982664774676767729863101188528827405_wp, &
                               -0.05715959421140685436681459970502471634_wp, &
                                0.09849855103848558320961101178888983150_wp, &
                                0.06654601552456084978615342374581437947_wp, &
                                0.09073479542748112726465375642050504556_wp, &
                                0.08432289325330803924891866923939606351_wp, &
                                0.04529095628204896774513180907141004447_wp]
  !
  ! NOTE: This is only found in the file, and can be computed by summing each
  !       row of the A matrix of the Butcher coefficients.
  real(wp), parameter :: OEPERK_10_5_FSAL_Coef_c(1:10) = &
                              [ zero, &
                                0.2597883554788674084039539165398464630_wp, &
                                0.09904573247592460887087003212056568980_wp, &
                                0.2155511890524058691860390281856497503_wp, &
                                0.5007950088969676776844289399972611534_wp, &
                                0.5592251911688643533787800688765883636_wp, &
                                0.5449986978853637084972622392134732553_wp, &
                                0.7615224694532590139829150720490417596_wp, &
                                0.8427062083267360939805493320684741215_wp, &
                                0.9152209805057669959657927210873423883_wp]
  !
  real(wp), parameter :: OEPERK_10_5_FSAL_Coef_g1(1:10) = &
                              [ zero, &
                                0.4043660121685749695640462197806189975_wp, &
                               -0.8503427289575839690883191973980814832_wp, &
                               -6.950894175262117526410215315179482885_wp, &
                                0.9238765192731084931855438934978371889_wp, &
                               -2.563178056509891340215942413817786020_wp, &
                                0.2545744879365226143946122067064118430_wp, &
                                0.3125831707411998258746812355492206137_wp, &
                               -0.7007114414440507927791249989236719346_wp, &
                                0.4839621016023833375810172323297465039_wp]
  !
  real(wp), parameter :: OEPERK_10_5_FSAL_Coef_g2(1:10) = &
                              [ one, &
                                0.6871467028161416909922221357014564412_wp, &
                                1.093024748914750833700799552463885117_wp, &
                                3.225975379607193001678365742708874597_wp, &
                                1.041153702510101386914019859778740444_wp, &
                                1.292821487912164945157744726076279306_wp, &
                                0.7391462755788122847651304143259254381_wp, &
                                0.1239129251371800313941948224441873274_wp, &
                                0.1842753472370123193132193302369345580_wp, &
                                0.05712788998796583446479387686662738843_wp]
  !
  real(wp), parameter :: OEPERK_10_5_FSAL_Coef_g3(1:10) = &
                              [ zero, &
                                zero, &
                                zero, &
                               -2.393405133244194727221124311276648940_wp, &
                               -1.902854422421760920850597670305403139_wp, &
                               -2.820042207399977261483046412236557428_wp, &
                               -1.832698465277380999601896111079977378_wp, &
                               -0.2199094483084671192328083958346519535_wp, &
                               -0.4082430635847870963724591602173546218_wp, &
                               -0.1377669797880289713535665985132703979_wp]
  !
  real(wp), parameter :: OEPERK_10_5_FSAL_Coef_d(1:10) = &
                              [ one, &
                               -0.1331778419508803397033287009506932673_wp, &
                                0.8260422814750207498262063505871077303_wp, &
                                1.513700425755728332485300719652378197_wp, &
                               -1.305810059935023735972298885749903694_wp, &
                                3.036678802924163246003321318996156380_wp, &
                               -1.449458274398895177922690618003584514_wp, &
                                3.834313899176362315089976408899373409_wp, &
                                4.122293760012985409330881631526514714_wp, &
                                zero]
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  err_dt_atol = 5.0e-5_wp
  err_dt_rtol = 5.0e-5_wp
  !
  safety_factor = 0.99_wp
  !
  safety_min = 0.2_wp
  safety_max = 10.0_wp
  !
  n = 1
  !
  if (abs(Runge_Kutta_Scheme) == EPERK_DormandPrince_5thOrder) then
    !
    EPERK_is_FSAL = true
    !
    num_rk_stages = 6
    !
   !pid_beta = [one]      ! F2003 AUTO-REALLOC
   !pid_beta = [0.85_wp]  ! F2003 AUTO-REALLOC
    pid_beta = [0.7_wp,-0.2_wp] ! F2003 AUTO-REALLOC
    !
    order_exp = one5
    !
  else if (abs(Runge_Kutta_Scheme) == EPERK_DormandPrince_8thOrder) then
    !
    EPERK_is_FSAL = fals
    !
    num_rk_stages = 12
    !
   !pid_beta = [one]      ! F2003 AUTO-REALLOC
   !pid_beta = [0.85_wp]  ! F2003 AUTO-REALLOC
    pid_beta = [0.7_wp,-0.2_wp] ! F2003 AUTO-REALLOC
    !
    order_exp = one8
    !
    err_dt_atol = 1.0e-5_wp
    err_dt_rtol = 1.0e-5_wp
    !
    n = 2
    !
  else if (Runge_Kutta_Scheme < 0) then
    !
    ! Non-FSAL optimized embedded-pair ERK methods
    !
    Runge_Kutta_Scheme = abs(Runge_Kutta_Scheme)
    !
    EPERK_is_FSAL = fals
    !
    if (Runge_Kutta_Scheme == OEPERK_3rdOrder) then
      !
      num_rk_stages = 5
      !
      oeperk_coef_c    = OEPERK_5_3_Coef_c(:)    ! F2003 AUTO-REALLOC
      oeperk_coef_b    = OEPERK_5_3_Coef_beta(:) ! F2003 AUTO-REALLOC
      oeperk_coef_bhat = OEPERK_5_3_Coef_bhat(:) ! F2003 AUTO-REALLOC
      oeperk_coef_g1   = OEPERK_5_3_Coef_g1(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_g2   = OEPERK_5_3_Coef_g2(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_g3   = OEPERK_5_3_Coef_g3(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_d    = OEPERK_5_3_Coef_d(:)    ! F2003 AUTO-REALLOC
      !
      pid_beta = [0.64_wp,-0.31_wp,0.04_wp] ! F2003 AUTO-REALLOC
      !
      order_exp = one3
      !
    else if (Runge_Kutta_Scheme == OEPERK_4thOrder) then
      !
      num_rk_stages = 9
      !
      oeperk_coef_c    = OEPERK_9_4_Coef_c(:)    ! F2003 AUTO-REALLOC
      oeperk_coef_b    = OEPERK_9_4_Coef_beta(:) ! F2003 AUTO-REALLOC
      oeperk_coef_bhat = OEPERK_9_4_Coef_bhat(:) ! F2003 AUTO-REALLOC
      oeperk_coef_g1   = OEPERK_9_4_Coef_g1(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_g2   = OEPERK_9_4_Coef_g2(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_g3   = OEPERK_9_4_Coef_g3(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_d    = OEPERK_9_4_Coef_d(:)    ! F2003 AUTO-REALLOC
      !
      pid_beta = [0.25_wp,-0.12_wp] ! F2003 AUTO-REALLOC
      !
      order_exp = one4
      !
     !err_dt_atol = 1.0e-5_wp
     !err_dt_rtol = 1.0e-5_wp
     !!
     !safety_factor = 0.90_wp
      !
      safety_min = 0.5_wp
      safety_max = 1.5_wp
      !
    else if (Runge_Kutta_Scheme == OEPERK_5thOrder) then
      !
      num_rk_stages = 10
      !
      oeperk_coef_c    = OEPERK_10_5_Coef_c(:)    ! F2003 AUTO-REALLOC
      oeperk_coef_b    = OEPERK_10_5_Coef_beta(:) ! F2003 AUTO-REALLOC
      oeperk_coef_bhat = OEPERK_10_5_Coef_bhat(:) ! F2003 AUTO-REALLOC
      oeperk_coef_g1   = OEPERK_10_5_Coef_g1(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_g2   = OEPERK_10_5_Coef_g2(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_g3   = OEPERK_10_5_Coef_g3(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_d    = OEPERK_10_5_Coef_d(:)    ! F2003 AUTO-REALLOC
      !
      pid_beta = [0.47_wp,-0.20_wp,0.06_wp] ! F2003 AUTO-REALLOC
      !
      order_exp = one5
      !
      err_dt_atol = 1.0e-7_wp
      err_dt_rtol = 1.0e-7_wp
      !
      safety_factor = 0.90_wp
      !
      safety_min = 0.50_wp
      safety_max = 1.50_wp
      !
    end if
    !
  else
    !
    EPERK_is_FSAL = true
    !
    ! FSAL optimized embedded-pair ERK methods
    !
    if (Runge_Kutta_Scheme == OEPERK_3rdOrder) then
      !
      num_rk_stages = 5
      !
      oeperk_coef_c    = OEPERK_5_3_FSAL_Coef_c(:)    ! F2003 AUTO-REALLOC
      oeperk_coef_b    = OEPERK_5_3_FSAL_Coef_beta(:) ! F2003 AUTO-REALLOC
      oeperk_coef_bhat = OEPERK_5_3_FSAL_Coef_bhat(:) ! F2003 AUTO-REALLOC
      oeperk_coef_g1   = OEPERK_5_3_FSAL_Coef_g1(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_g2   = OEPERK_5_3_FSAL_Coef_g2(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_g3   = OEPERK_5_3_FSAL_Coef_g3(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_d    = OEPERK_5_3_FSAL_Coef_d(:)    ! F2003 AUTO-REALLOC
      !
      pid_beta = [0.70_wp,-0.23_wp] ! F2003 AUTO-REALLOC
      !
      order_exp = one3
      !
    else if (Runge_Kutta_Scheme == OEPERK_4thOrder) then
      !
      num_rk_stages = 9
      !
      oeperk_coef_c    = OEPERK_9_4_FSAL_Coef_c(:)    ! F2003 AUTO-REALLOC
      oeperk_coef_b    = OEPERK_9_4_FSAL_Coef_beta(:) ! F2003 AUTO-REALLOC
      oeperk_coef_bhat = OEPERK_9_4_FSAL_Coef_bhat(:) ! F2003 AUTO-REALLOC
      oeperk_coef_g1   = OEPERK_9_4_FSAL_Coef_g1(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_g2   = OEPERK_9_4_FSAL_Coef_g2(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_g3   = OEPERK_9_4_FSAL_Coef_g3(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_d    = OEPERK_9_4_FSAL_Coef_d(:)    ! F2003 AUTO-REALLOC
      !
      pid_beta = [0.38_wp,-0.18_wp,0.01_wp] ! F2003 AUTO-REALLOC
      !
      order_exp = one4
      !
     !err_dt_atol = 1.0e-6_wp
     !err_dt_rtol = 1.0e-6_wp
     !!
     !safety_factor = 0.90_wp
      !
      safety_min = 0.50_wp
      safety_max = 1.50_wp
      !
    else !if (Runge_Kutta_Scheme == OEPERK_5thOrder) then
      !
      num_rk_stages = 10
      !
      oeperk_coef_c    = OEPERK_10_5_FSAL_Coef_c(:)    ! F2003 AUTO-REALLOC
      oeperk_coef_b    = OEPERK_10_5_FSAL_Coef_beta(:) ! F2003 AUTO-REALLOC
      oeperk_coef_bhat = OEPERK_10_5_FSAL_Coef_bhat(:) ! F2003 AUTO-REALLOC
      oeperk_coef_g1   = OEPERK_10_5_FSAL_Coef_g1(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_g2   = OEPERK_10_5_FSAL_Coef_g2(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_g3   = OEPERK_10_5_FSAL_Coef_g3(:)   ! F2003 AUTO-REALLOC
      oeperk_coef_d    = OEPERK_10_5_FSAL_Coef_d(:)    ! F2003 AUTO-REALLOC
      !
      pid_beta = [0.45_wp,-0.13_wp] ! F2003 AUTO-REALLOC
      !
      order_exp = one5
      !
      err_dt_atol = 1.0e-7_wp
      err_dt_rtol = 1.0e-7_wp
      !
      safety_factor = 0.90_wp
      !
      safety_min = 0.50_wp
      safety_max = 1.50_wp
      !
    end if
    !
  end if
  !
  n_controllers = size(pid_beta)
  ! !
  ! NOTE: error_dt needs to be initialized to one because it is used
  !       in a product function/operator type loop. For the cases of:
  !         1) n_controllers >= 2, and computing the new time step size
  !            before the 2nd time step, OR
  !         2) n_controllers == 3, and computing the new time step size
  !            before the 3rd time step
  !       the previous time step errors have yet to be defined and so their
  !       contribution to the product needs to be one.
  !
  error_dt  = init_array(nq,n_controllers,one) ! F2003 AUTO-REALLOC
  !
  my_error  = init_array(n*nq,zero)      ! F2003 AUTO-REALLOC
  all_error = init_array(n*nq,ncpu,zero) ! F2003 AUTO-REALLOC
  !
end subroutine init_eperk_coefficients
!
!###############################################################################
!
subroutine update_Optimized_EPERK(nst,np,dtsp,residual,uoldsp,usp,rk_sum, &
                                  uepsp,rhs_save)
  !
  !.. Use Statements ..
  use ovar, only : time
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rk_sum(:,:)
  real(wp), intent(inout) :: uepsp(:,:)
  real(wp), intent(inout) :: rhs_save(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,m,nr
  real(wp) :: c,b,bh,g1,g2,g3,d,dt
  real(wp) :: u1ec,u1mx,u1my,u1mz,u1ee,u1m
  real(wp) :: u2ec,u2mx,u2my,u2mz,u2ee,u2m
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_Optimized_EPERK"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nst == 1) then
    pseudo_time = time
    rk_sum(:,:) = zero
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uepsp)
    do m = 1,nq
      do n = 1,np
        uepsp(n,m) = usp(n,m)
      end do
    end do
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(residual,rhs_save)
    do m = 1,nq
      do n = 1,np
        rhs_save(n,m) = residual(n,m)
      end do
    end do
  end if
  !
  c  = oeperk_coef_c(nst)
  b  = oeperk_coef_b(nst)
  bh = oeperk_coef_bhat(nst)
  g1 = oeperk_coef_g1(nst)
  g2 = oeperk_coef_g2(nst)
  g3 = oeperk_coef_g3(nst)
  d  = oeperk_coef_d(nst)
  !
 !if (mypnum == glb_root) &
 !write (iout,1) nst,c,b,bh,g1,g2,g3,d
 !1 format (/,"nst = ",i1," :    c = ",f15.12,/, &
 !                       10x,"   b = ",f15.12,/, &
 !                       10x,"bhat = ",f15.12,/, &
 !                       10x,"  g1 = ",f15.12,/, &
 !                       10x,"  g2 = ",f15.12,/, &
 !                       10x,"  g3 = ",f15.12,/, &
 !                       10x,"   d = ",f15.12)
  !
  pseudo_time = pseudo_time + c*dtsp(1)
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual,rk_sum,uepsp)
    do n = 1,np
      !
      dt = dtsp(n)
      !
      u1ec = usp(n,nec)
      u1mx = usp(n,nmx)
      u1my = usp(n,nmy)
      u1mz = usp(n,nmz)
      u1ee = usp(n,nee)
      !
      u2ec = rk_sum(n,nec) + d*u1ec
      u2mx = rk_sum(n,nmx) + d*u1mx
      u2my = rk_sum(n,nmy) + d*u1my
      u2mz = rk_sum(n,nmz) + d*u1mz
      u2ee = rk_sum(n,nee) + d*u1ee
      !
      usp(n,nec) = g2*u2ec + g1*u1ec + g3*uoldsp(n,nec) + b*dt*residual(n,nec)
      usp(n,nmx) = g2*u2mx + g1*u1mx + g3*uoldsp(n,nmx) + b*dt*residual(n,nmx)
      usp(n,nmy) = g2*u2my + g1*u1my + g3*uoldsp(n,nmy) + b*dt*residual(n,nmy)
      usp(n,nmz) = g2*u2mz + g1*u1mz + g3*uoldsp(n,nmz) + b*dt*residual(n,nmz)
      usp(n,nee) = g2*u2ee + g1*u1ee + g3*uoldsp(n,nee) + b*dt*residual(n,nee)
      !
      rk_sum(n,nec) = u2ec
      rk_sum(n,nmx) = u2mx
      rk_sum(n,nmy) = u2my
      rk_sum(n,nmz) = u2mz
      rk_sum(n,nee) = u2ee
      !
      uepsp(n,nec) = uepsp(n,nec) + bh*dt*residual(n,nec)
      uepsp(n,nmx) = uepsp(n,nmx) + bh*dt*residual(n,nmx)
      uepsp(n,nmy) = uepsp(n,nmy) + bh*dt*residual(n,nmy)
      uepsp(n,nmz) = uepsp(n,nmz) + bh*dt*residual(n,nmz)
      uepsp(n,nee) = uepsp(n,nee) + bh*dt*residual(n,nee)
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual,rk_sum,uepsp)
    do n = 1,np
      !
      dt = dtsp(n)
      !
      u1ec = usp(n,nec)
      u1mx = usp(n,nmx)
      u1my = usp(n,nmy)
      u1ee = usp(n,nee)
      !
      u2ec = rk_sum(n,nec) + d*u1ec
      u2mx = rk_sum(n,nmx) + d*u1mx
      u2my = rk_sum(n,nmy) + d*u1my
      u2ee = rk_sum(n,nee) + d*u1ee
      !
      usp(n,nec) = g2*u2ec + g1*u1ec + g3*uoldsp(n,nec) + b*dt*residual(n,nec)
      usp(n,nmx) = g2*u2mx + g1*u1mx + g3*uoldsp(n,nmx) + b*dt*residual(n,nmx)
      usp(n,nmy) = g2*u2my + g1*u1my + g3*uoldsp(n,nmy) + b*dt*residual(n,nmy)
      usp(n,nee) = g2*u2ee + g1*u1ee + g3*uoldsp(n,nee) + b*dt*residual(n,nee)
      !
      rk_sum(n,nec) = u2ec
      rk_sum(n,nmx) = u2mx
      rk_sum(n,nmy) = u2my
      rk_sum(n,nee) = u2ee
      !
      uepsp(n,nec) = uepsp(n,nec) + bh*dt*residual(n,nec)
      uepsp(n,nmx) = uepsp(n,nmx) + bh*dt*residual(n,nmx)
      uepsp(n,nmy) = uepsp(n,nmy) + bh*dt*residual(n,nmy)
      uepsp(n,nee) = uepsp(n,nee) + bh*dt*residual(n,nee)
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,dtsp,residual) &
    !$ACC&                          PRESENT(rk_sum,uepsp)
    do m = 1,nq
      do n = 1,np
        !
        u1m = usp(n,m)
        !
        u2m = rk_sum(n,m) + d*u1m
        !
        usp(n,m) = g2*u2m + g1*u1m + g3*uoldsp(n,m) + b*dtsp(n)*residual(n,m)
        !
        rk_sum(n,m) = u2m
        !
        uepsp(n,m) = uepsp(n,m) + bh*dtsp(n)*residual(n,m)
        !
      end do
    end do
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_Optimized_EPERK
!
!###############################################################################
!
subroutine update_Optimized_EPERK_FSAL(np,dtsp,residual,uepsp)
  !
  !.. Use Statements ..
  use ovar, only : time,num_rk_stages
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: uepsp(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,m,nr
  real(wp) :: bhat
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_Optimized_EPERK_FSAL"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  !
  !
  bhat = oeperk_coef_bhat(num_rk_stages+1)
  !
  !
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(dtsp,residual,uepsp)
    do n = 1,np
      !
      uepsp(n,nec) = uepsp(n,nec) + bhat*dtsp(n)*residual(n,nec)
      uepsp(n,nmx) = uepsp(n,nmx) + bhat*dtsp(n)*residual(n,nmx)
      uepsp(n,nmy) = uepsp(n,nmy) + bhat*dtsp(n)*residual(n,nmy)
      uepsp(n,nmz) = uepsp(n,nmz) + bhat*dtsp(n)*residual(n,nmz)
      uepsp(n,nee) = uepsp(n,nee) + bhat*dtsp(n)*residual(n,nee)
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(dtsp,residual,uepsp)
    do n = 1,np
      !
      uepsp(n,nec) = uepsp(n,nec) + bhat*dtsp(n)*residual(n,nec)
      uepsp(n,nmx) = uepsp(n,nmx) + bhat*dtsp(n)*residual(n,nmx)
      uepsp(n,nmy) = uepsp(n,nmy) + bhat*dtsp(n)*residual(n,nmy)
      uepsp(n,nee) = uepsp(n,nee) + bhat*dtsp(n)*residual(n,nee)
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(dtsp,residual,uepsp)
    do m = 1,nq
      do n = 1,np
        !
        uepsp(n,m) = uepsp(n,m) + bhat*dtsp(n)*residual(n,m)
        !
      end do
    end do
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_Optimized_EPERK_FSAL
!
!###############################################################################
!
end module eperk_mod
