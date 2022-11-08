module eperk_mod
  !
  use module_kind_types
  !
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq,ncs
  !
  use ovar, only : force_nd_subroutines
  !
  implicit none
  !
  private
  !
  public :: update_EPERK
  !
  !
  !
  logical(lk), save :: dtmin_was_modified = fals
  real(wp),    save :: previous_dtmin = zero
  !
  real(wp), save, allocatable :: error_dt(:,:)
  real(wp), save, allocatable :: my_error(:)
  real(wp), save, allocatable :: all_error(:,:)
  !
  real(wp), save, allocatable :: uepsp(:,:)
  !
  real(wp), save :: pseudo_time
  !
  real(wp), parameter :: err_dt_atol = 1.0e-4_wp
  real(wp), parameter :: err_dt_rtol = 1.0e-4_wp
  !
 !integer, parameter :: n_controllers = 1
 !real(wp), parameter :: pid_beta(1:n_controllers) = [one]
  !
  integer, parameter :: n_controllers = 2
  real(wp), parameter :: pid_beta(1:n_controllers) = [-0.2_wp,0.7_wp]
  !
 !integer, parameter :: n_controllers = 3
 !real(wp), parameter :: pid_beta(1:n_controllers) = [0.08_wp,-0.28_wp,0.41_wp]
 !real(wp), parameter :: pid_beta(1:n_controllers) = [0.1_wp,-0.34_wp,0.49_wp]
  !
contains
!
!###############################################################################
!
subroutine update_EPERK(nst,iter,Apply_Filter,Apply_Limiter, &
                        dtmin,final_timestep_from_stop_request)
  !
  !.. Use Statements ..
  use geovar,   only : n_solpts
  use flowvar,  only : residual,usp,uoldsp,rssprk,dtsp
  use ovar,     only : time
  use time_mod, only : stabilize_usp,get_dtsp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  integer,     intent(in) :: iter
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  real(wp),   intent(out) :: dtmin
  integer,     intent(in) :: final_timestep_from_stop_request
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_EPERK"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nst == 1) then
    !
    ! Get the time step size if this is the first Runge-Kutta stage
    !
    if (iter == 1) then
      !
      ! If this is the first time step of the solver, compute
      ! the time step size using the traditional method
      !
      rssprk(:,:,:) = zero
      call get_dtsp(dtmin,final_timestep_from_stop_request)
      !
    else
      !
      ! If this is not the first time step, use the time error from the
      ! embedded pair Runge-Kutta method to adaptively adjust the time
      ! step size
      !
     !call get_dtsp(dtmin,final_timestep_from_stop_request)
      call get_EPERK_adaptive_dtsp(iter,dtmin,final_timestep_from_stop_request)
      !
    end if
  end if
  !
  ! Update the solution using the current Runge-Kutta stage
  !
  if (nst == 1) then
    pseudo_time = time
    call update_EPERK_DP54_S1(n_solpts,dtsp,residual,uoldsp, &
                              usp,rssprk,pseudo_time)
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
  ! Apply filtering or limiting to stabilize the solution, and then update
  ! the thermodynamic and transport properties at the solution points
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter,caller_name=pname)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_EPERK
!
!###############################################################################
!
subroutine get_EPERK_adaptive_dtsp(iter,dtmin,final_timestep_from_stop_request)
  !
  !.. Use Statements ..
  use geovar,   only : n_solpts
  use flowvar,  only : residual,usp,uoldsp,rssprk,dtsp
  use ovar,     only : d_t
  use time_mod, only : checks_on_dtmin
  !
  !.. Formal Arguments ..
  integer,   intent(in) :: iter
  real(wp), intent(out) :: dtmin
  integer,   intent(in) :: final_timestep_from_stop_request
  !
  !.. Local Scalars ..
  integer  :: n,m
  real(wp) :: dtmin_save
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_EPERK_adaptive_dtsp"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (.not.allocated(error_dt)) then
   !error_dt = init_array(nq,n_controllers,zero) ! F2003 AUTO-REALLOC
    error_dt = init_array(nq,n_controllers,one) ! F2003 AUTO-REALLOC
  end if
  !
  if (.not.allocated(uepsp)) then
    uepsp = init_array(n_solpts,nq,zero) ! F2003 AUTO-REALLOC
  end if
  !
  if (.not.allocated(my_error)) then
    my_error = init_array(nq,zero) ! F2003 AUTO-REALLOC
  end if
  !
  if (.not.allocated(all_error)) then
    all_error = init_array(nq,ncpu,zero) ! F2003 AUTO-REALLOC
  end if
  !
  ! If dtmin was modified by checks_on_dtmin last time step, we need to
  ! reset it to the value before it was modified so the adaptive controller
  ! doesnt use the modified value when computing the next time step size.
  !
  if (dtmin_was_modified) then
    dtmin = previous_dtmin
    dtmin_was_modified = fals
  end if
  !
  ! Shift the time step error values before computing
  ! the time step error from the previous time step.
  !
  do n = 2,n_controllers
    do m = 1,nq
      error_dt(m,n-1) = error_dt(m,n)
    end do
  end do
  error_dt(:,n_controllers) = zero
  !
  ! Compute the embedded pair solution using the RHS computed with
  ! the solution from the last time step.
  !
  call update_EPERK_DP54_S7(n_solpts,dtsp,residual,uoldsp,usp,rssprk,uepsp)
  !
  ! Compute the time error using the embedded pair solutions
  !
 !call get_embedded_pair_time_error_wtd_nrm(n_solpts,usp,uepsp)
  call get_embedded_pair_time_error_reg_nrm(n_solpts,usp,uepsp)
  !
  call compute_adaptive_time_step(n_solpts,dtsp,dtmin)
  !
  ! Perform some checks on dtmin to make sure the time of the new solution
  ! doesnt skip past a milestone for output or final solution time
  !
  dtmin_save = dtmin
  !
  call checks_on_dtmin(dtmin,final_timestep_from_stop_request)
  !
  ! Check to see if dtmin was modified in checks_on_dtmin. If it was modified,
  ! save the unmodified value to the module variable previous_dtmin to use
  ! when computing the next time step size.
  !
  if (abs(dtmin_save-dtmin) > eps6) then
    dtmin_was_modified = true
    previous_dtmin = dtmin_save
  end if
  !
  ! Update d_t with the new dtmin
  !
  d_t = dtmin
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_EPERK_adaptive_dtsp
!
!###############################################################################
!
subroutine get_embedded_pair_time_error_wtd_nrm(np,usp,uepsp)
  !
  !.. Use Statements ..
  use geovar, only : n_global_solpts
 !use ovar,   only : err_dt_atol,err_dt_rtol
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: np
  real(wp), intent(in) :: usp(:,:)
  real(wp), intent(in) :: uepsp(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,m
  real(wp) :: lsum,lnum,lden,lerr,up,um
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
      !
      lnum = up-um
      lden = err_dt_atol + err_dt_rtol*max( abs(up) , abs(um) )
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
    lerr = sqrt(lsum)
    !
    error_dt(m,n_controllers) = one/lerr
    !
  end do
  !
end subroutine get_embedded_pair_time_error_wtd_nrm
!
!###############################################################################
!
subroutine get_embedded_pair_time_error_reg_nrm(np,usp,uepsp)
  !
  !.. Use Statements ..
  use geovar, only : n_global_solpts
 !use ovar,   only : err_dt_atol
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: np
  real(wp), intent(in) :: usp(:,:)
  real(wp), intent(in) :: uepsp(:,:)
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
   !error_dt(m,n_controllers) = eps9/lerr
    error_dt(m,n_controllers) = err_dt_atol/lerr
    !
  end do
  !
end subroutine get_embedded_pair_time_error_reg_nrm
!
!###############################################################################
!
subroutine compute_adaptive_time_step(np,dtsp,dtmin)
  !
  !.. Use Statements ..
  use geovar, only : n_global_solpts
 !use ovar,   only : err_dt_atol
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
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
  real(wp), parameter :: one5 = one/five
  !
continue
  !
  do m = 1,nq
    ltot = one
    do i = 1,n_controllers
      ltot = ltot * (error_dt(m,i)**(pid_beta(i)*one5))
    end do
    delta_t(m) = ltot
  end do
  !
  dtnew = 0.9_wp * minval(delta_t(1:nq)) * dtmin
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
subroutine update_EPERK_DP54_S1(np,dtsp,residual,uoldsp,usp,rssprk,pseudo_time)
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
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,residual,dtsp)
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
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        rssprk(n,m,1) = dtsp(n)*residual(n,m)
        !
        usp(n,m) = uoldsp(n,m) + a1*rssprk(n,m,1)
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
subroutine update_EPERK_DP54_S7(np,dtsp,residual,uoldsp,usp,rssprk,uepsp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: usp(:,:)
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
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,residual,dtsp)
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
!
end module eperk_mod
!!! !
!!! !
!!!   !
!!!   dt = dtsp(n)
!!!   !
!!!   t0 = time
!!!   t1 = t0+dt*c1
!!!   t2 = t0+dt*c2
!!!   t3 = t0+dt*c3
!!!   t4 = t0+dt*c4
!!!   t5 = t0+dt*c5
!!!   t6 = t0+dt*c6
!!!   !
!!!   ! INPUT: U0,USP from previous time step
!!!   !        dt = time step size from previous time step
!!!   !        t
!!!   !        k0-k5 = values from previous time steps
!!!   !
!!!   ! First compute new time step size
!!!   !
!!!   k6 = dt * odefun(time,usp)
!!!   u7 = u0 + a1s7*k0 + a3s7*k2 + a4s7*k3 + a5s7*k4 + a6s7*k5 + a7s7*k6
!!!   !
!!!   E = norm(usp-u7)
!!!   dt = 0.9 * dt * (TOL/E)**(1/5)
!!!   !
!!!   ! Next update u0 with unchanged usp
!!!   !
!!!   u0 = usp
!!!   !
!!!   ! Stage 1
!!!   !
!!!   k0 = dt * odefun(t0, u0)
!!!   u1 = u0 + k0*a1s1
!!!   !
!!!   k1 = dt * odefun(t1,u1)
!!!   u2 = u0 + k0*a1s2 + k1*a2s2
!!!   !
!!!   k2 = dt * odefun(t2,u2)
!!!   u3 = u0 + k0*a1s3 + k1*a2s3 + k2*a3s3
!!!   !
!!!   k3 = dt * odefun(t3,u3)
!!!   u4 = u0 + k0*a1s4 + k1*a2s4 + k2*a3s4 + k3*a4s4
!!!   !
!!!   k4 = dt * odefun(t4,u4)
!!!   u5 = u0 + k0*a1s5 + k1*a2s5 + k2*a3s5 + k3*a4s5 + k4*a5s5
!!!   !
!!!   k5 = dt * odefun(t5,u5)
!!!   u6 = u0 + k0*a1s6 + k2*a3s6 + k3*a4s6 + k4*a5s6 + k5*a6s6
!!!   !
!!!   time = t6
!!!   usp = u6
!!!   !
!!!   unp1 = u6
!!!   
!!!   E = norm(unp1-usm1)
!!!   dt = 0.9 * dt * (TOL/E)^(1/5)   
!!!   
!!!   y(:,k) = unp1
!!!   t(k) = t(k-1) + dt
!!!   k = k + 1
!!! !
!!! !
!!! !
!!! !
!!! !
!!! !
!!! !
!!! !
!!! !
!!! !
!!!   !
!!!   dt = dtsp(n)
!!!   !
!!!   u0 = uoldsp(n)
!!!   t0 = time
!!!   !
!!!   !
!!!   k1 = dt * odefun(t0, u0)
!!!   !
!!!   t1 = t0+dt*(c1-c0)
!!!   t2 = t1+dt*(c2-c1)
!!!   t3 = t2+dt*(c3-c2)
!!!   t4 = t3+dt*(c4-c3)
!!!   t5 = t4+dt*(c5-c4)
!!!   t6 = t5+dt*(c6-c5)
!!!   !
!!!   !
!!!   u1   = u0 + a1s1*dt*f(u0)
!!!   k2   = dt * odefun(u1)
!!!   !
!!!   u2   = u0 + a1s2*dt*f(u0) + a2s2*dt*f(u1)
!!!   k3   = dt * odefun(u2)
!!!   !
!!!   u3   = u0 + a1s3*dt*f(u0) + a2s3*dt*f(u1) + a3s3*dt*f(u2)
!!!   k4   = dt * odefun(u3)
!!!   !
!!!   u4   = u0 + a1s4*dt*f(u0) + a2s4*dt*f(u1) + &
!!!               a3s4*dt*f(u2) + a4s4*dt*f(u3)
!!!   k5   = dt * odefun(u4)
!!!   !
!!!   u5   = u0 + a1s5*dt*f(u0) + a2s5*dt*f(u1) + &
!!!               a3s5*dt*f(u2) + a4s5*dt*f(u3) + a5s5*dt*f(u4)
!!!   k6   = dt * odefun(u5)
!!!   !
!!!   u6   = u0 + a1s6*dt*f(u0)                 + &
!!!               a3s6*dt*f(u2) + a4s6*dt*f(u3) + &
!!!               a5s6*dt*f(u4) + a6s6*dt*f(u5)
!!!   k7   = dt * odefun(u6)
!!!   !
!!!   !
!!!   unp1 = u6
!!!   usm1 = u0 + a1s7*dt*f(u0)                 + &
!!!               a3s7*dt*f(u2) + a4s7*dt*f(u3) + &
!!!               a5s7*dt*f(u4) + a6s7*dt*f(u5) + a7s7*dt*f(u6)
!!! 
!!! 
!!! 
!!! 
!!! 
!!! 
!!! 
!!! 
!!!   dt = dtsp(n)
!!!   !
!!!   u0 = uoldsp(n)
!!!   t0 = time
!!!   !
!!!   !
!!!   k1 = dt * odefun(t0, u0)
!!!   !
!!!   t1 = t0+dt*c1
!!!   u1 = u0+k1*a1s1
!!!   k2 = dt * odefun(t1,u1)
!!!   !
!!!   t2 = t0+dt*c2
!!!   u2 = u0+k1*a1s2 + k2*a2s2
!!!   k3 = dt * odefun(t2,u2)
!!!   !
!!!   t3 = t0+dt*c3
!!!   u3 = u0 + k1*a1s3 + k2*a2s3 + k3*a3s3
!!!   k4 = dt * odefun(t3,u3)
!!!   !
!!!   t4 = t0+dt*c4
!!!   u4 = u0 + k1*a1s4 + k2*a2s4 + k3*a3s4 + k4*a4s4
!!!   k5 = dt * odefun(t4,u4)
!!!   !
!!!   t5 = t0+dt*c5
!!!   u5 = u0 + k1*a1s5 + k2*a2s5 + k3*a3s5 + k4*a4s5 + k5*a5s5
!!!   k6 = dt * odefun(t5,u5)
!!!   !
!!!   t6 = t0+dt*c6
!!!   u6 = u0 + k1*a1s6 + k3*a3s6 + k4*a4s6 + k5*a5s6 + k6*a6s6
!!!   k7 = dt * odefun(t6,u6)
!!!   !
!!!   !
!!!   unp1 = u6
!!!   usm1 = u0 + a1s7*k1 + a3s7*k3 + a4s7*k4 + a5s7*k5 + a6s7*k6 + a7s7*k7
!!!   
!!!   E = norm(unp1-usm1)
!!!   dt = 0.9 * dt * (TOL/E)^(1/5)   
!!!   
!!!   y(:,k) = unp1
!!!   t(k) = t(k-1) + dt
!!!   k = k + 1
!!! 
!!! 
!!! 
!!!   ! stage 1
!!!   !
