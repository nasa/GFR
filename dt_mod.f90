module dt_mod
  !
!#define PROFILE_ON
  !
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nmb,nme
  use eqn_idx_mod, only : ncs,npres,nvisc,ngamm
  use timer_mod, only : solver_timer
  !
  implicit none
  !
  private
  !
  public :: get_dtsp
  public :: checks_on_dtmin
  !
contains
!
!###############################################################################
!
subroutine get_dtsp(dtmin,final_timestep_from_stop_request,skip_updates)
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_solpts
  use ovar,        only : Timestep_Type
  use ovar,        only : d_t,cfl ! modified
  use flowvar_mod, only : dtsp ! modified
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: dtmin
  integer,     intent(in) :: final_timestep_from_stop_request
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: skip_updates
  !
  !.. Local Scalars ..
  integer     :: n
  real(wp)    :: dtmin_new
  logical(lk) :: update_dt_values
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: dt_sp(:)
  !
  !. Local Saved Scalars ..
  integer, save :: idx = -1
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_dtsp"
  !
continue
  !
  call solver_timer(entering_procedure,pname,idx)
  !
  call debug_timer(entering_procedure,pname)
  !
  update_dt_values = true
  if (present(skip_updates)) then
    update_dt_values = (.not. skip_updates)
  end if
  !
  allocate ( dt_sp(1:n_solpts) )
  !$ACC DATA CREATE(dt_sp)
  !
  if (Timestep_Type == Constant_Timestep) then
    !
    call get_current_cfl(cfl)
    !
    call get_constant_timestep(dtmin_new,dt_sp)
    !
  else if (abs(Timestep_Type) == Cell_Timestep) then
    !
    call get_current_cfl(cfl)
    !
    if (Timestep_Type == Global_Cell_Timestep) then
      call get_cell_global_timestep(dtmin_new,dt_sp)
    else if (Timestep_Type == Local_Cell_Timestep) then
      call get_cell_local_timestep(dtmin_new,dt_sp)
    end if
    !
  else if (abs(Timestep_Type) == Point_Timestep) then
    !
    call get_current_cfl(cfl)
    !
    if (Timestep_Type == Global_Point_Timestep) then
      call get_point_global_timestep(dtmin_new,dt_sp)
    else if (Timestep_Type == Local_Point_Timestep) then
      call get_point_local_timestep(dtmin_new,dt_sp)
    end if
    !
  end if
  !
  ! Perform some checks on dtmin_new to make sure the time of the new solution
  ! doesnt skip past a milestone for output or final solution time
  !
  call checks_on_dtmin(dtmin_new,dt_sp,final_timestep_from_stop_request, &
                       skip_updates)
  !
  ! If the optional input skip_updates is not present or is present and false,
  ! update dtmin, dtsp, and d_t
  !
  if (update_dt_values) then
    !
    !$ACC PARALLEL LOOP PRESENT(dt_sp,dtsp)
    do n = 1,n_solpts
      dtsp(n) = dt_sp(n)
    end do
    !
    ! Update dtmin and d_t with the new dtmin
    !
    dtmin = dtmin_new
    d_t = dtmin
    !
  else
    !
    ! This means that skip_updates is present and true, which should only
    ! happen if using an EPERK method and this subroutine is being used just
    ! get an estimated CFL number for the adaptive time step size computed
    ! by the EPERK method.
    !
    ! The dtmin_new computed from the above subroutines was multiplied by the
    ! current CFL value, so we need to divide that out in order to compare this
    ! dtmin_new with the input argument dtmin that was computed by the EPERK
    ! method.
    !
    dtmin_new = dtmin_new / CFL
    !
    ! Now get the approximate CFL by dividing the input argument
    ! dtmin (which was computed by the EPERK adaptive time step
    ! controller) by dtmin_new (which is the CFL-based dtmin with
    ! an assumed CFL of 1).
    !
    CFL = dtmin / dtmin_new
    !
  end if
  !
  !$ACC END DATA
  deallocate ( dt_sp )
  !
  call debug_timer(leaving_procedure,pname)
  !
  call solver_timer(leaving_procedure,pname,idx)
  !
end subroutine get_dtsp
!
!###############################################################################
!
subroutine checks_on_dtmin(dtmin,dt_sp,final_timestep_from_stop_request, &
                           skip_updates)
  !
  !.. Use Statements ..
  use ovar, only : itestcase,Final_Time
  use ovar, only : output_frequency
  use ovar, only : restart_frequency
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: dtmin
  real(wp), intent(inout) :: dt_sp(:)
  integer,     intent(in) :: final_timestep_from_stop_request
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: skip_updates
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "checks_on_dtmin"
  !
continue
  !
  if (itestcase == Taylor_Green_Vortex) then
    call check_tstar_for_TaylorGreen(dtmin,dt_sp,skip_updates)
  else if (abs(itestcase) == Taylor_Couette_Flow) then
    call check_tstar_for_TaylorCouette(dtmin,dt_sp,skip_updates)
  end if
  !
  ! Check to see if we need to output a new solution or restart
  ! file based on output_frequency or restart_frequency
  !
  if (output_frequency > zero .or. restart_frequency > zero) then
    call check_output_frequency(dtmin,dt_sp,skip_updates)
  end if
  !
  ! Make sure the timestep size doesnt make the final time exceed Final_Time
  !
  if (Final_time > zero .or. final_timestep_from_stop_request /= 0) then
    call check_final_time(dtmin,dt_sp,final_timestep_from_stop_request, &
                          skip_updates)
  end if
  !
end subroutine checks_on_dtmin
!
!###############################################################################
!
subroutine get_current_cfl(cfl)
  !
  !.. Use Statements ..
  use ovar, only : cfl_beg,cfl_end,cfl_cycles,itcur
  use ovar, only : cfl_cycle_start,cfl_cycle_end
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: cfl
  !
  !.. Local Scalars ..
  real(wp) :: cycle_ratio
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_current_cfl"
  !
continue
  !
  if (itcur < cfl_cycle_end) then
    if (cfl_cycles < 2) then
      cycle_ratio = one
    else
      cycle_ratio = real(itcur-cfl_cycle_start-1,kind=wp) / &
                    real(cfl_cycles-1,kind=wp)
    end if
    !
    cfl = cfl_beg + (cfl_end - cfl_beg) * cycle_ratio
  else
    cfl = cfl_end
  end if
  !
end subroutine get_current_cfl
!
!###############################################################################
!
subroutine get_constant_timestep(dtmin,dt_sp)
  !
  !.. Use Statements ..
  use ovar,        only : constant_dt,CFL
  use geomvar_mod, only : n_solpts
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  real(wp), intent(out) :: dt_sp(:)
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_constant_timestep"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  dtmin = CFL * constant_dt
  !
  !$ACC PARALLEL LOOP PRESENT(dt_sp)
  do n = 1,n_solpts
    dt_sp(n) = dtmin
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_constant_timestep
!
!###############################################################################
!
subroutine get_cell_global_timestep(dtmin,dt_sp)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nfbnd,nface,n_solpts,ncell
  use geomvar_mod, only : face,cell,n_bnd_flxpts,n_flxpts
  !
  use ovar, only : CFL,governing_equations
  !
  use metrics_mod, only : facenrm,facejac
  !
  use flowvar_mod, only : faceusp
  use flowvar_mod, only : facegpsp
  use flowvar_mod, only : facexyz
  use flowvar_mod, only : dt_save
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  real(wp), intent(out) :: dt_sp(:)
  !
  !.. Local Scalars ..
  integer  :: n,kf,nc,nf,f1,f2,n1,n2
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: eigni_fp,eignv_fp
  real(wp) :: dtl,hl,lcl
  real(wp) :: dtr,hr,rcl
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: viscl(:)
  real(wp), allocatable :: viscr(:)
  real(wp), allocatable :: eigni(:,:)
  real(wp), allocatable :: psq(:,:)
  real(wp), allocatable :: dt_fp(:,:)
  real(wp), allocatable :: dt_cell(:)
  !
  integer, parameter :: make_p_plus_one = 0
 !integer, parameter :: make_p_plus_one = 1
  !
  character(len=*), parameter :: pname = "get_cell_global_timestep"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  allocate ( viscl(1:n_flxpts) )
  allocate ( viscr(1:n_flxpts) )
  allocate ( psq(1:n_flxpts,1:2) )
  allocate ( eigni(1:n_flxpts,1:2) )
  allocate ( dt_fp(1:n_flxpts,1:2) )
  !
  !$ACC DATA CREATE(viscl,viscr,eigni,psq,dt_fp)
  !
  n1 = 1
  n2 = n_flxpts
  !
  ! First, compute the left and right inviscid eigenvalues at each flux point
  !
  call invs_eigenvalue(n1,n2,faceusp,facegpsp,facenrm,facexyz,eigni)
  !
  ! Copy the left and right viscosities into the local arrays.
  ! NOTE: This is done because facegpsp(:,nvisc,:) doesnt exist for an
  !       inviscid simulation, so using a local array with viscosity set to
  !       zero prevents us from having to have either inviscid and viscous
  !       versions of this subroutine or if statements at the lowest level
  !       of the loop to compute the final time step size.
  !
  if (governing_equations == NavierStokes_Eqns) then
    !$ACC PARALLEL LOOP PRESENT(facegpsp,viscl,viscr)
    do n = n1,n2
      viscl(n) = facegpsp(n,nvisc,1)
      viscr(n) = facegpsp(n,nvisc,2)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(viscl,viscr)
    do n = n1,n2
      viscl(n) = zero
      viscr(n) = zero
    end do
  end if
  !
  dtmin = huge(zero)
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(face,cell,psq)
  do nf = 1,nfbnd
    !
    f1 = face(nf)%beg_fp
    f2 = face(nf)%end_fp
    !
    left_cell  = face(nf)%left%cell
    left_geom  = cell(left_cell)%geom
    left_order = cell(left_cell)%order
    !
    lcl = cl(left_geom)
    !
    !$ACC LOOP VECTOR
    do kf = f1,f2
      psq(kf,1) = lcl * real( (left_order+make_p_plus_one)**2 , kind=wp )
      psq(kf,2) = lcl * real( (left_order+make_p_plus_one)**2 , kind=wp )
    end do
    !
  end do
  !
  !$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) &
  !$ACC&                   PRESENT(face,cell,psq)
  do nf = nfbnd+1,nface
    !
    f1 = face(nf)%beg_fp
    f2 = face(nf)%end_fp
    !
    left_cell  = face(nf)%left%cell
    right_cell = face(nf)%right%cell
    !
    left_geom  = cell(left_cell)%geom
    left_order = cell(left_cell)%order
    !
    right_geom  = cell(right_cell)%geom
    right_order = cell(right_cell)%order
    !
    lcl = cl(left_geom)
    rcl = cl(right_geom)
    !
    !$ACC LOOP VECTOR
    do kf = f1,f2
      psq(kf,1) = lcl * real( (left_order +make_p_plus_one)**2 , kind=wp )
      psq(kf,2) = rcl * real( (right_order+make_p_plus_one)**2 , kind=wp )
    end do
    !
  end do
  !
  !$ACC PARALLEL LOOP REDUCTION(min:dtmin) &
  !$ACC&              PRESENT(facejac,eigni,viscl,psq,dt_fp)
  do kf = 1,n_bnd_flxpts
    !
    ! NOTE: facejac is actually area/facejac
    !
    hl = psq(kf,1) * half * facejac(kf,1)
    !
    eigni_fp = eigni(kf,1)
    eignv_fp = hl * viscl(kf)
    !
    dtl = CFL / (hl * (eigni_fp + eignv_fp))
    !
#define IGNORE_NEGATIVE_TIMESTEPS
#ifdef IGNORE_NEGATIVE_TIMESTEPS
    if (dtl <= zero) then
      dtl = huge(zero)
    end if
#endif
    !
    dtmin = min(dtmin,dtl)
    !
    dt_fp(kf,1) = dtl
    dt_fp(kf,2) = dtl
    !
  end do
  !
  !$ACC PARALLEL LOOP REDUCTION(min:dtmin) &
  !$ACC&              PRESENT(facejac,eigni,viscl,viscr,psq,dt_fp)
  do kf = n_bnd_flxpts+1,n_flxpts
    !
    ! NOTE: facejac is actually area/facejac
    !
    hl = half * psq(kf,1) * facejac(kf,1)
    hr = half * psq(kf,2) * facejac(kf,2)
    !
    eigni_fp = max( eigni(kf,1) , eigni(kf,2) )
    eignv_fp = max( hl * viscl(kf) , hr * viscr(kf) )
    !
    dtl = CFL / (max(hl,hr) * (eigni_fp + eignv_fp))
#ifdef IGNORE_NEGATIVE_TIMESTEPS
    if (dtl <= zero) then
      dtl = huge(zero)
    end if
#endif
    dtmin = min(dtmin,dtl)
   !dtl = CFL / (hl * (eigni_fp + eignv_fp))
   !dtr = CFL / (hr * (eigni_fp + eignv_fp))
   !dtmin = min(dtmin,dtl,dtr)
    !
    dt_fp(kf,1) = CFL / (hl * (eigni(kf,1) + hl*viscl(kf)) )
    dt_fp(kf,2) = CFL / (hr * (eigni(kf,2) + hr*viscr(kf)) )
    !
  end do
  !
  if (allocated(dt_save)) then
    dt_cell = init_array(ncell,huge(zero)) ! F2003 AUTO-REALLOC
    do nf = 1,nfbnd
      f1 = face(nf)%beg_fp
      f2 = face(nf)%end_fp
      left_cell = face(nf)%left%cell
      dtl = minval(dt_fp(f1:f2,1))
      dt_cell(left_cell) = min(dt_cell(left_cell),dtl)
    end do
    do nf = nfbnd+1,nface
      f1 = face(nf)%beg_fp
      f2 = face(nf)%end_fp
      left_cell = face(nf)%left%cell
      right_cell = face(nf)%right%cell
      dtl = minval(dt_fp(f1:f2,1))
      dtr = minval(dt_fp(f1:f2,2))
      dt_cell(left_cell) = min(dt_cell(left_cell),dtl)
      dt_cell(right_cell) = min(dt_cell(right_cell),dtr)
    end do
    do nc = 1,ncell
      n1 = cell(nc)%beg_sp
      n2 = cell(nc)%end_sp
      dt_save(n1:n2) = dt_cell(nc)
    end do
    deallocate ( dt_cell )
  end if
  !
  ! Get the minimum timestep across all processors
  !
  if (ncpu > 1) then
    ! Update the value of dtmin on the CPU
    !$ACC UPDATE SELF(dtmin)
    call mpi_allreduce(MPI_IN_PLACE,dtmin,1_int_mpi,mpi_flttyp, &
                       MPI_MIN,MPI_COMM_WORLD,mpierr)
    ! Update the value of dtmin on the GPU
    !$ACC UPDATE DEVICE(dtmin)
  end if
  !
  !$ACC PARALLEL LOOP PRESENT(dt_sp)
  do n = 1,n_solpts
    dt_sp(n) = dtmin
  end do
  !
  !$ACC END DATA
  deallocate ( psq )
  deallocate ( eigni )
  deallocate ( viscr )
  deallocate ( viscl )
  deallocate ( dt_fp )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_cell_global_timestep
!
!###############################################################################
!
subroutine get_cell_local_timestep(dtmin,dt_sp)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nfbnd,nface,ncell
  use geomvar_mod, only : face,cell,n_flxpts
  !
  use ovar, only : CFL,governing_equations
  !
  use metrics_mod, only : facenrm,facejac
  !
  use flowvar_mod, only : faceusp
  use flowvar_mod, only : facegpsp
  use flowvar_mod, only : facexyz
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  real(wp), intent(out) :: dt_sp(:)
  !
  !.. Local Scalars ..
  integer  :: k,n,kf,nc,nf,f0,nfp,n1,n2
  integer  :: host_cell,host_geom,host_order
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: dtc,eigni_fp,eignv_fp
  real(wp) :: dtl,hl,psql
  real(wp) :: dtr,hr,psqr
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: viscl(:)
  real(wp), allocatable :: viscr(:)
  real(wp), allocatable :: dt_cell(:)
  real(wp), allocatable :: eigni(:,:)
  !
  integer, parameter :: make_p_plus_one = 0
 !integer, parameter :: make_p_plus_one = 1
  !
  character(len=*), parameter :: pname = "get_cell_local_timestep"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  allocate ( viscl(1:n_flxpts) )
  allocate ( viscr(1:n_flxpts) )
  allocate ( dt_cell(1:ncell) )
  allocate ( eigni(1:n_flxpts,1:2) )
  !
  !$ACC DATA CREATE(viscl,viscr,eigni,dt_cell)
  !
  n1 = 1
  n2 = n_flxpts
  !
  ! First, compute the left and right inviscid eigenvalues at each flux point
  !
  call invs_eigenvalue(n1,n2,faceusp,facegpsp,facenrm,facexyz,eigni)
  !
  ! Copy the left and right viscosities into the local arrays.
  ! NOTE: This is done because facegpsp(:,nvisc,:) doesnt exist for an
  !       inviscid simulation, so using a local array with viscosity set to
  !       zero prevents us from having to have either inviscid and viscous
  !       versions of this subroutine or if statements at the lowest level
  !       of the loop to compute the final time step size.
  !
  if (governing_equations == NavierStokes_Eqns) then
    !$ACC PARALLEL LOOP PRESENT(facegpsp,viscl,viscr)
    do n = n1,n2
      viscl(n) = facegpsp(n,nvisc,1)
      viscr(n) = facegpsp(n,nvisc,2)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(viscl,viscr)
    do n = n1,n2
      viscl(n) = zero
      viscr(n) = zero
    end do
  end if
  !
  dtmin = huge(zero)
  !$ACC PARALLEL LOOP
  do nc = 1,ncell
    dt_cell(nc) = dtmin
  end do
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(face,cell,facejac) &
  !$ACC&              PRESENT(dt_cell,eigni,viscl)
  do nf = 1,nfbnd
    !
    left_cell = face(nf)%left%cell  ! left  cell on face
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    left_geom  = cell(left_cell)%geom  ! left cell geometry
    left_order = cell(left_cell)%order ! left cell order
    !
    psql = half*cl(left_geom)*real( (left_order+make_p_plus_one)**2 , kind=wp )
    !
    dtl = dt_cell(left_cell)
    !
    !$ACC LOOP VECTOR REDUCTION(min:dtl)
    do k = 1,nfp
      !
      kf = f0 + k
      !
      ! NOTE: facejac is actually area/facejac
      !
      hl = psql * facejac(kf,1)
      !
      eigni_fp = eigni(kf,1)
      eignv_fp = hl * viscl(kf)
      !
      ! Compute the timestep for this flux point
      !
      dtl = min( dtl , CFL / (hl * (eigni_fp + eignv_fp) ) )
      !
    end do
    !
    dt_cell(left_cell) = dtl
    !
  end do
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(face,cell,facejac) &
  !$ACC&              PRESENT(dt_cell,eigni,viscl,viscr)
  do nf = nfbnd+1,nface
    !
    left_cell  = face(nf)%left%cell  ! left  cell on face
    right_cell = face(nf)%right%cell ! right cell on face
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    left_geom  = cell(left_cell)%geom  ! left cell geometry
    left_order = cell(left_cell)%order ! left cell order
    !
    right_geom  = cell(right_cell)%geom  ! right cell geometry
    right_order = cell(right_cell)%order ! right cell order
    !
    psql = half * cl(left_geom) &
                * real( (left_order +make_p_plus_one)**2 , kind=wp )
    psqr = half * cl(right_geom) &
                * real( (right_order+make_p_plus_one)**2 , kind=wp )
    !
    dtl = dt_cell(left_cell)
    dtr = dt_cell(right_cell)
    !
    !$ACC LOOP VECTOR REDUCTION(min:dtl,dtr)
    do k = 1,nfp
      !
      kf = f0 + k
      !
      ! Left and right cell jacobians at the face point
      ! NOTE: facejac is actually area/facejac
      !
      hl = psql * facejac(kf,1)
      hr = psqr * facejac(kf,2)
      !
      eigni_fp = max( eigni(kf,1) , eigni(kf,2) )
      eignv_fp = max( hl * viscl(kf) , hr * viscr(kf) )
      !
      ! Compute the timestep for this flux point
      !
      dtl = min( dtl , CFL / (hl * (eigni_fp + eignv_fp)) )
      dtr = min( dtr , CFL / (hr * (eigni_fp + eignv_fp)) )
      !
    end do
    !
    dt_cell(left_cell)  = dtl
    dt_cell(right_cell) = dtr
    !
  end do
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) &
  !$ACC&              REDUCTION(min:dtmin) &
  !$ACC&              PRESENT(cell,dt_cell,dt_sp)
  do nc = 1,ncell
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    dtc = dt_cell(nc)
    dtmin = min(dtmin,dtc)
    do n = n1,n2
      dt_sp(n) = dtc
    end do
  end do
  !
  !$ACC END DATA
  deallocate ( eigni )
  deallocate ( dt_cell )
  deallocate ( viscr )
  deallocate ( viscl )
  !
  ! Get the minimum timestep across all processors
  !
  if (ncpu > 1) then
    ! Update the value of dtmin on the CPU
    !$ACC UPDATE SELF(dtmin)
    call mpi_allreduce(MPI_IN_PLACE,dtmin,1_int_mpi,mpi_flttyp, &
                       MPI_MIN,MPI_COMM_WORLD,mpierr)
    ! Update the value of dtmin on the GPU
    !$ACC UPDATE DEVICE(dtmin)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_cell_local_timestep
!
!###############################################################################
!
!pure &
subroutine invs_eigenvalue(n1,n2,face_cv,face_gp,facenrm,facexyz,eigni)
  !
  !.. Use Statements ..
  use mrf_mod, only : is_in_mrf
  use mrf_mod, only : abs_to_mrf,mrf_to_abs
  use mrf_mod, only : switch_vel_ref_frame
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp), intent(inout) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp),    intent(in) :: facexyz(:,:)
  real(wp), intent(inout) :: eigni(:,:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "invs_eigenvalue"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = size(facenrm,dim=1)
  !
  ! First, compute the speed of sound component of the eigenvalue and
  ! store it in eigni
  !
  call invs_eigenvalue_aspd(n1,n2,face_cv,face_gp,eigni)
  !
  ! If a rotating reference frame is being used, convert the face
  ! velocity/momentum from the absolute/inertial frame to the rotating frame.
  !
  if ( is_in_mrf() ) then
    call switch_vel_ref_frame(n1,n2,face_cv,facexyz,abs_to_mrf)
  end if
  !
  ! Now, add the velocity component of the eigenvalue to the eigni array
  !
  if (nr == 3) then      ! this should also work for ncs > 0
    call invs_eigenvalue_3d(n1,n2,face_cv,facenrm,eigni)
  else if (nr == 2) then ! this should also work for ncs > 0
    call invs_eigenvalue_2d(n1,n2,face_cv,facenrm,eigni)
  else
    call invs_eigenvalue_nd(n1,n2,face_cv,facenrm,eigni)
  end if
  !
  ! If a rotating reference frame is being used, convert the face
  ! velocity/momentum from the rotating from back to the
  ! absolute/inertial frame.
  !
  if ( is_in_mrf() ) then
    call switch_vel_ref_frame(n1,n2,face_cv,facexyz,mrf_to_abs)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine invs_eigenvalue
!
!###############################################################################
!
pure subroutine invs_eigenvalue_aspd(n1,n2,face_cv,face_gp,aspd)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp), intent(inout) :: aspd(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,n
  real(wp) :: rho_inv,pres,gamm
  !
continue
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv,face_gp,aspd)
    do i = 1,2
      do n = n1,n2
        !
        rho_inv = one/face_cv(n,nec,i)
        pres = face_gp(n,npres,i)
        gamm = face_gp(n,ngamm,i)
        !
        aspd(n,i) = sqrt( max( gamm*pres*rho_inv , rndoff ) )
        !
      end do
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(face_cv,face_gp,aspd)
    do i = 1,2
      do n = n1,n2
        !
        rho_inv = one/face_cv(n,nec,i)
        pres = face_gp(n,npres,i)
        !
        aspd(n,i) = sqrt( max( gam*pres*rho_inv , rndoff ) )
        !
      end do
    end do
  end if
  !
end subroutine invs_eigenvalue_aspd
!
!###############################################################################
!
pure subroutine invs_eigenvalue_3d(n1,n2,face_cv,facenrm,eigni)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: eigni(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: rho_inv
  real(wp) :: velx,vely,velz,vface
  real(wp) :: rnx,rny,rnz,unx,uny,unz,ds
  !
continue
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,facenrm,eigni)
  do n = n1,n2
    !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds = sqrt(rnx*rnx + rny*rny + rnz*rnz)
    unx = rnx/ds
    uny = rny/ds
    unz = rnz/ds
#endif
    !
    ! Left side of face
    !
    rho_inv = one/face_cv(n,nec,1)
    !
    velx = face_cv(n,nmx,1)*rho_inv
    vely = face_cv(n,nmy,1)*rho_inv
    velz = face_cv(n,nmz,1)*rho_inv
    !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    vface = unx*velx + uny*vely + unz*velz
#else
    vface = sqrt(velx*velx + vely*vely + velz*velz)
#endif
    !
    eigni(n,1) = eigni(n,1) + abs(vface)
    !
    ! Right side of face
    !
    rho_inv = one/face_cv(n,nec,2)
    !
    velx = face_cv(n,nmx,2)*rho_inv
    vely = face_cv(n,nmy,2)*rho_inv
    velz = face_cv(n,nmz,2)*rho_inv
    !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    vface = unx*velx + uny*vely + unz*velz
#else
    vface = sqrt(velx*velx + vely*vely + velz*velz)
#endif
    !
    eigni(n,2) = eigni(n,2) + abs(vface)
    !
  end do
  !
end subroutine invs_eigenvalue_3d
!
!###############################################################################
!
pure subroutine invs_eigenvalue_2d(n1,n2,face_cv,facenrm,eigni)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: eigni(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: rho_inv
  real(wp) :: velx,vely,vface
  real(wp) :: rnx,rny,unx,uny,ds
  !
continue
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,facenrm,eigni)
  do n = n1,n2
    !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds = sqrt(rnx*rnx + rny*rny)
    unx = rnx/ds
    uny = rny/ds
#endif
    !
    ! Left side of face
    !
    rho_inv = one/face_cv(n,nec,1)
    !
    velx = face_cv(n,nmx,1)*rho_inv
    vely = face_cv(n,nmy,1)*rho_inv
    !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    vface = unx*velx + uny*vely
#else
    vface = sqrt(velx*velx + vely*vely)
#endif
    !
    eigni(n,1) = eigni(n,1) + abs(vface)
    !
    ! Right side of face
    !
    rho_inv = one/face_cv(n,nec,2)
    !
    velx = face_cv(n,nmx,2)*rho_inv
    vely = face_cv(n,nmy,2)*rho_inv
    !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    vface = unx*velx + uny*vely
#else
    vface = sqrt(velx*velx + vely*vely)
#endif
    !
    eigni(n,2) = eigni(n,2) + abs(vface)
    !
  end do
  !
end subroutine invs_eigenvalue_2d
!
!###############################################################################
!
pure subroutine invs_eigenvalue_nd(n1,n2,face_cv,facenrm,eigni)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: eigni(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nd
  real(wp) :: rho_inv,vface
  !
  !.. Local Arrays ..
  real(wp) :: unrm(1:size(facenrm,dim=1))
  real(wp) :: vel(nmb:nme)
  !
continue
  !
  nd = size(facenrm,dim=1)
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,facenrm,eigni,unrm,vel)
  do n = n1,n2
    !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    unrm(1:nd) = facenrm(1:nd,n)/norm2(facenrm(1:nd,n))
#endif
    !
    ! Left side of face
    !
    vel(nmb:nme) = face_cv(n,nmb:nme,1)/face_cv(n,nec,1)
    !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    vface = dot_product( unrm(1:nd) , vel(nmb:nme) )
#else
    vface = norm2(vel(nmb:nme))
#endif
    !
    eigni(n,1) = eigni(n,1) + abs(vface)
    !
    ! Right side of face
    !
    vel(nmb:nme) = face_cv(n,nmb:nme,2)/face_cv(n,nec,2)
    !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    vface = dot_product( unrm(1:nd) , vel(nmb:nme) )
#else
    vface = norm2(vel(nmb:nme))
#endif
    !
    eigni(n,2) = eigni(n,2) + abs(vface)
    !
  end do
  !
end subroutine invs_eigenvalue_nd
!
!###############################################################################
!
subroutine get_point_global_timestep(dtmin,dt_sp)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr,n_solpts,xyz
  use ovar,        only : gam,CFL
  use ovar,        only : moving_ref_frame
  use metrics_mod, only : time_metrics
  use flowvar_mod, only : usp
  use flowvar_mod, only : gpsp
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  real(wp), intent(out) :: dt_sp(:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: da,dx,dy,dz,dtpt,dtc
  real(wp) :: velx,vely,velz
  real(wp) :: rho_inv,pres,gamm
  real(wp) :: x,y,z,x_r,y_r,z_r
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: aspd(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_point_global_timestep"
  !
 !integer, parameter :: make_p_plus_one = 0
  integer, parameter :: make_p_plus_one = 1
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  allocate ( aspd(1:n_solpts) )
  !$ACC DATA CREATE(aspd)
  !
  call moving_ref_frame%get_rotation(x_r,y_r,z_r)
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP PRESENT(usp,gpsp,aspd)
    do n = 1,n_solpts
      !
      rho_inv = one/usp(n,nec)
      pres = gpsp(n,npres)
      gamm = gpsp(n,ngamm)
      !
      aspd(n) = sqrt( max( gamm*pres*rho_inv , rndoff ) )
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(usp,gpsp,aspd)
    do n = 1,n_solpts
      !
      rho_inv = one/usp(n,nec)
      pres = gpsp(n,npres)
      !
      aspd(n) = sqrt( max( gam*pres*rho_inv , rndoff ) )
      !
    end do
  end if
  !
  !
  !
  dtmin = huge(zero)
  !
  if (nr == 3) then ! this should also work for ncs > 0
    !
    !$ACC PARALLEL LOOP REDUCTION(min:dtmin) &
    !$ACC&              PRESENT(usp,aspd,time_metrics,dt_sp,xyz)
    do n = 1,n_solpts
      !
      x = xyz(1,n)
      y = xyz(2,n)
      z = xyz(3,n)
      !
      dx = time_metrics(1,n)
      dy = time_metrics(2,n)
      dz = time_metrics(3,n)
      da = time_metrics(4,n)
      !
      rho_inv = one/usp(n,nec)
      !
      velx = usp(n,nmx)*rho_inv - (y_r*z - z_r*y)
      vely = usp(n,nmy)*rho_inv - (z_r*x - x_r*z)
      velz = usp(n,nmz)*rho_inv - (x_r*y - y_r*x)
      !
      dtc = dx*abs(velx) + dy*abs(vely) + dz*abs(velz) + da*aspd(n)
      !
      dtpt = CFL/dtc
      !
      dtmin = min( dtpt , dtmin )
      !
      dt_sp(n) = dtpt
      !
    end do
    !
  else if (nr == 2) then ! this should also work for ncs > 0
    !
    !$ACC PARALLEL LOOP REDUCTION(min:dtmin) &
    !$ACC&              PRESENT(usp,aspd,time_metrics,dt_sp)
    do n = 1,n_solpts
      !
      dx = time_metrics(1,n)
      dy = time_metrics(2,n)
      da = time_metrics(3,n)
      !
      rho_inv = one/usp(n,nec)
      !
      velx = usp(n,nmx)*rho_inv
      vely = usp(n,nmy)*rho_inv
      !
      dtc = dx*abs(velx) + dy*abs(vely) + da*aspd(n)
      !
      dtpt = CFL/dtc
      !
      dtmin = min( dtpt , dtmin )
      !
      dt_sp(n) = dtpt
      !
    end do
    !
  else
    !
    !$ACC PARALLEL LOOP REDUCTION(min:dtmin) &
    !$ACC&              PRESENT(usp,aspd,time_metrics,dt_sp)
    do n = 1,n_solpts
      !
      dx = time_metrics(1,n)
      da = time_metrics(2,n)
      !
      dtc = dx*abs(usp(n,nmx))/usp(n,nec) + da*aspd(n)
      !
      dtpt = CFL/dtc
      !
      dtmin = min( dtpt , dtmin )
      !
      dt_sp(n) = dtpt
      !
    end do
    !
  end if
  !
  !$ACC END DATA
  deallocate ( aspd )
  !
  ! Get the minimum timestep across all processors
  !
  if (ncpu > 1) then
    ! Update the value of dtmin on the CPU
    !$ACC UPDATE SELF(dtmin)
    call mpi_allreduce(MPI_IN_PLACE,dtmin,1_int_mpi,mpi_flttyp, &
                       MPI_MIN,MPI_COMM_WORLD,mpierr)
    ! Update the value of dtmin on the GPU
    !$ACC UPDATE DEVICE(dtmin)
  end if
  !
  !$ACC PARALLEL LOOP PRESENT(dt_sp)
  do n = 1,n_solpts
    dt_sp(n) = dtmin
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_point_global_timestep
!
!###############################################################################
!
subroutine get_point_local_timestep(dtmin,dt_sp)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr,n_solpts,xyz
  use ovar,        only : gam,CFL
  use ovar,        only : moving_ref_frame
  use metrics_mod, only : time_metrics
  use flowvar_mod, only : usp
  use flowvar_mod, only : gpsp
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  real(wp), intent(out) :: dt_sp(:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: da,dx,dy,dz,dtpt,dtc
  real(wp) :: velx,vely,velz
  real(wp) :: rho_inv,pres,gamm
  real(wp) :: x,y,z,x_r,y_r,z_r
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: aspd(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_point_local_timestep"
  !
 !integer, parameter :: make_p_plus_one = 0
  integer, parameter :: make_p_plus_one = 1
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  allocate ( aspd(1:n_solpts) )
  !$ACC DATA CREATE(aspd)
  !
  call moving_ref_frame%get_rotation(x_r,y_r,z_r)
  !
  if (ncs > 0) then
    !$ACC PARALLEL LOOP PRESENT(usp,gpsp,aspd)
    do n = 1,n_solpts
      !
      rho_inv = one/usp(n,nec)
      pres = gpsp(n,npres)
      gamm = gpsp(n,ngamm)
      !
      aspd(n) = sqrt( max( gamm*pres*rho_inv , rndoff ) )
      !
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(usp,gpsp,aspd)
    do n = 1,n_solpts
      !
      rho_inv = one/usp(n,nec)
      pres = gpsp(n,npres)
      !
      aspd(n) = sqrt( max( gam*pres*rho_inv , rndoff ) )
      !
    end do
  end if
  !
  !
  !
  dtmin = huge(zero)
  !
  if (nr == 3) then ! this should also work for ncs > 0
    !
    !$ACC PARALLEL LOOP REDUCTION(min:dtmin) &
    !$ACC&              PRESENT(usp,aspd,time_metrics,dt_sp,xyz)
    do n = 1,n_solpts
      !
      x = xyz(1,n)
      y = xyz(2,n)
      z = xyz(3,n)
      !
      dx = time_metrics(1,n)
      dy = time_metrics(2,n)
      dz = time_metrics(3,n)
      da = time_metrics(4,n)
      !
      rho_inv = one/usp(n,nec)
      !
      velx = usp(n,nmx)*rho_inv - (y_r*z - z_r*y)
      vely = usp(n,nmy)*rho_inv - (z_r*x - x_r*z)
      velz = usp(n,nmz)*rho_inv - (x_r*y - y_r*x)
      !
      dtc = dx*abs(velx) + dy*abs(vely) + dz*abs(velz) + da*aspd(n)
      !
      dtpt = CFL/dtc
      !
      dtmin = min( dtpt , dtmin )
      !
      dt_sp(n) = dtpt
      !
    end do
    !
  else if (nr == 2) then ! this should also work for ncs > 0
    !
    !$ACC PARALLEL LOOP REDUCTION(min:dtmin) &
    !$ACC&              PRESENT(usp,aspd,time_metrics,dt_sp)
    do n = 1,n_solpts
      !
      dx = time_metrics(1,n)
      dy = time_metrics(2,n)
      da = time_metrics(3,n)
      !
      rho_inv = one/usp(n,nec)
      !
      velx = usp(n,nmx)*rho_inv
      vely = usp(n,nmy)*rho_inv
      !
      dtc = dx*abs(velx) + dy*abs(vely) + da*aspd(n)
      !
      dtpt = CFL/dtc
      !
      dtmin = min( dtpt , dtmin )
      !
      dt_sp(n) = dtpt
      !
    end do
    !
  else
    !
    !$ACC PARALLEL LOOP REDUCTION(min:dtmin) &
    !$ACC&              PRESENT(usp,aspd,time_metrics,dt_sp)
    do n = 1,n_solpts
      !
      dx = time_metrics(1,n)
      da = time_metrics(2,n)
      !
      dtc = dx*abs(usp(n,nmx))/usp(n,nec) + da*aspd(n)
      !
      dtpt = CFL/dtc
      !
      dtmin = min( dtpt , dtmin )
      !
      dt_sp(n) = dtpt
      !
    end do
    !
  end if
  !
  !$ACC END DATA
  deallocate ( aspd )
  !
  ! Get the minimum timestep across all processors
  !
  if (ncpu > 1) then
    ! Update the value of dtmin on the CPU
    !$ACC UPDATE SELF(dtmin)
    call mpi_allreduce(MPI_IN_PLACE,dtmin,1_int_mpi,mpi_flttyp, &
                       MPI_MIN,MPI_COMM_WORLD,mpierr)
    ! Update the value of dtmin on the GPU
    !$ACC UPDATE DEVICE(dtmin)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_point_local_timestep
!
!###############################################################################
!
subroutine check_tstar_for_TaylorGreen(dtmin,dt_sp,skip_updates)
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_solpts
  use ovar,        only : time,time_ref
  use ovar,        only : write_TaylorGreen_8s_solution ! modified
  use ovar,        only : write_TaylorGreen_full_solution ! modified
  !
  !.. Formal Arguments ..
  real(wp),    intent(inout) :: dtmin
  real(wp),    intent(inout) :: dt_sp(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: skip_updates
  !
  !.. Local Scalars ..
  integer     :: n,ichk
  real(wp)    :: tstar,dtstar,ideal_tstar
  logical(lk) :: reset_dtsp,checking_for_eight
  logical(lk) :: write_tgv_8s_solution
  logical(lk) :: write_tgv_full_solution
  logical(lk) :: perform_updates
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_tstar_for_TaylorGreen"
  !
continue
  !
  perform_updates = true
  if (present(skip_updates)) then
    perform_updates = (.not. skip_updates)
  end if
  !
  reset_dtsp = fals
  !
  tstar = time*time_ref
  dtstar = dtmin*time_ref
  !
  checking_for_eight = fals
  if (tstar < three) then
    ideal_tstar = three
  else if (tstar < five) then
    ideal_tstar = five
  else if (tstar < seven) then
    ideal_tstar = seven
  else if (tstar < eight) then
    ideal_tstar = eight
    checking_for_eight = true
  else if (tstar < nine) then
    ideal_tstar = nine
  else if (tstar < 11.0_wp) then
    ideal_tstar = 11.0_wp
  else if (tstar < 15.0_wp) then
    ideal_tstar = 15.0_wp
  else
   !ideal_tstar = 20.0_wp
    ideal_tstar = 25.0_wp
  end if
  !
  ! Use the function check_against_ideal_time to see if the next 1 or 2 time
  ! steps will exceed the ideal output time.
  !
  ichk = check_against_ideal_time(tstar,dtstar,ideal_tstar)
  !
  write_tgv_8s_solution = write_TaylorGreen_8s_solution
  write_tgv_full_solution = write_TaylorGreen_full_solution
  !
  if (ichk == -1) then
    !
    ! If the simulation time after the upcoming time step (using the current
    ! time step size) is within 1e-6 of the ideal time for output, just
    ! continue on using the current time step size (i.e. do not reset dtmin),
    ! and set the necessary flags to create the desired output at this
    ! simulation time.
    !
    if (checking_for_eight) then
      write_tgv_8s_solution = true
      write_tgv_full_solution = true
    else
      write_tgv_full_solution = true
    end if
    !
  else if (ichk == 1) then
    !
    ! If the simulation time after the upcoming time step (using the current
    ! time step size) will exceed the ideal time for output, adjust dtmin
    ! so that the simulation time after this time step is equal to the
    ! ideal time. Also set the necessary flags to create the desired output
    ! at this simulation time.
    ! NOTE: We are adding 1e-14 to the ideal time to make sure that the
    !       simulation time after this time step will exceed the ideal
    !       time, and there isnt some issue with round off errors making
    !       the updated simulation time smaller than the ideal time.
    !
    if (checking_for_eight) then
      write_tgv_8s_solution = true
      write_tgv_full_solution = true
    else
      write_tgv_full_solution = true
    end if
    !
    reset_dtsp = true
    !
    dtmin = ((ideal_tstar+eps14) - tstar) / time_ref
    !
  else if (ichk == 2) then
    !
    ! It looks like the next time step will exceed the ideal time for output.
    ! To smooth the modification of the time step size, adjust dtmin to
    ! half the remaining difference between the current simulation time and
    ! the ideal time.
    ! NOTE: We dont want to set the flags to create the desired output
    !       because the simulation time after this time step wont be the
    !       ideal time. This should happen during the next time step, and
    !       at that point the previous conditional above this one should
    !       be met.
    !
    reset_dtsp = true
    !
    dtmin = (half * (ideal_tstar - tstar)) / time_ref
    !
  end if
  !
  if (perform_updates) then
    write_TaylorGreen_8s_solution = write_tgv_8s_solution
    write_TaylorGreen_full_solution = write_tgv_full_solution
  end if
  !
 !!$ACC UPDATE DEVICE(dtmin)
  !
  if (reset_dtsp) then
    !$ACC PARALLEL LOOP PRESENT(dt_sp)
    do n = 1,n_solpts
      dt_sp(n) = dtmin
    end do
  end if
  !
end subroutine check_tstar_for_TaylorGreen
!
!###############################################################################
!
subroutine check_tstar_for_TaylorCouette(dtmin,dt_sp,skip_updates)
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_solpts
  use ovar,        only : time,time_ref
  use ovar,        only : write_TaylorGreen_full_solution ! modified
  !
  !.. Formal Arguments ..
  real(wp),    intent(inout) :: dtmin
  real(wp),    intent(inout) :: dt_sp(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: skip_updates
  !
  !.. Local Scalars ..
  integer     :: n,ichk
  real(wp)    :: tstar,dtstar,ideal_tstar
  logical(lk) :: reset_dtsp,checking_for_eight
  logical(lk) :: write_tcf_full_solution
  logical(lk) :: perform_updates
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_tstar_for_TaylorCouette"
  !
  real(wp), parameter :: t1 = 22.058_wp
  real(wp), parameter :: t2 = 24.016_wp
  real(wp), parameter :: t3 = 27.200_wp
  real(wp), parameter :: t4 = 33.500_wp
  !
continue
  !
  perform_updates = true
  if (present(skip_updates)) then
    perform_updates = (.not. skip_updates)
  end if
  !
  reset_dtsp = fals
  !
  tstar = time*time_ref
  dtstar = dtmin*time_ref
  !
  if (tstar < t1) then
    ideal_tstar = t1
  else if (tstar < t2) then
    ideal_tstar = t2
  else if (tstar < t3) then
    ideal_tstar = t3
  else if (tstar < t4) then
    ideal_tstar = t4
  else
    ideal_tstar = 80.0_wp
  end if
  !
  ! Use the function check_against_ideal_time to see if the next 1 or 2 time
  ! steps will exceed the ideal output time.
  !
  ichk = check_against_ideal_time(tstar,dtstar,ideal_tstar)
  !
  write_tcf_full_solution = write_TaylorGreen_full_solution
  !
  if (ichk == -1) then
    !
    ! If the simulation time after the upcoming time step (using the current
    ! time step size) is within 1e-6 of the ideal time for output, just
    ! continue on using the current time step size (i.e. do not reset dtmin),
    ! and set the necessary flags to create the desired output at this
    ! simulation time.
    !
    write_tcf_full_solution = true
    !
  else if (ichk == 1) then
    !
    ! If the simulation time after the upcoming time step (using the current
    ! time step size) will exceed the ideal time for output, adjust dtmin
    ! so that the simulation time after this time step is equal to the
    ! ideal time. Also set the necessary flags to create the desired output
    ! at this simulation time.
    ! NOTE: We are adding 1e-14 to the ideal time to make sure that the
    !       simulation time after this time step will exceed the ideal
    !       time, and there isnt some issue with round off errors making
    !       the updated simulation time smaller than the ideal time.
    !
    write_tcf_full_solution = true
    !
    reset_dtsp = true
    !
    dtmin = ((ideal_tstar+eps14) - tstar) / time_ref
    !
  else if (ichk == 2) then
    !
    ! It looks like the next time step will exceed the ideal time for output.
    ! To smooth the modification of the time step size, adjust dtmin to
    ! half the remaining difference between the current simulation time and
    ! the ideal time.
    ! NOTE: We dont want to set the flags to create the desired output
    !       because the simulation time after this time step wont be the
    !       ideal time. This should happen during the next time step, and
    !       at that point the previous conditional above this one should
    !       be met.
    !
    reset_dtsp = true
    !
    dtmin = (half * (ideal_tstar - tstar)) / time_ref
    !
  end if
  !
  if (perform_updates) then
    write_TaylorGreen_full_solution = write_tcf_full_solution
  end if
  !
 !!$ACC UPDATE DEVICE(dtmin)
  !
  if (reset_dtsp) then
    !$ACC PARALLEL LOOP PRESENT(dt_sp)
    do n = 1,n_solpts
      dt_sp(n) = dtmin
    end do
  end if
  !
end subroutine check_tstar_for_TaylorCouette
!
!###############################################################################
!
function check_against_ideal_time(current_time,current_dt,ideal_time) &
                           result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: current_time
  real(wp), intent(in) :: current_dt
  real(wp), intent(in) :: ideal_time
  !
  !.. Function Results ..
  integer :: return_value
  !
  !.. Local Parameters ..
  real(wp), parameter :: scale_factor = 1.75_wp
  real(wp), parameter :: eps_time = eps6
  !
continue
  !
  ! Default function result is 0, which means that the
  ! next 1 or 2 time steps will NOT exceed the ideal time.
  !
  return_value = 0
  !
  ! Perform some checks to see if the next 1 or 2 time steps will reach
  ! or exceed the ideal simulation time.
  ! NOTE: The parameter scale_factor is used to see if we will hit a
  !       simulation time within the next two time steps. We want to
  !       use a value somewhere between 1 and 2 so that we can adjust the
  !       next two time steps to reach that ideal simulation time smoothly,
  !       instead of only adjusting 1 time step size which can result in a
  !       really tiny time step size for the next time step if the difference
  !       between the ideal time and the current time happens to be just
  !       larger than the current time step size.
  !
  !
  if ((current_time+scale_factor*current_dt) > (ideal_time-eps_time)) then
    ! 
    if (abs(ideal_time - (current_time+current_dt)) <= eps_time) then
      !
      ! The simulation time after the upcoming time step (using the current
      ! time step size) is within 1e-6 of the ideal time. Return a value of
      ! -1 to indicate this condition.
      !
      return_value = -1
      !
    else if ((current_time+current_dt) > (ideal_time-eps_time)) then
      !
      ! The simulation time after the upcoming time step (using the current
      ! time step size) will exceed the ideal time, so return a value of
      ! 1 to indicate the ideal time will be exceeded in 1 time step.
      !
      return_value = 1
      !
    else
      !
      ! It looks like the next time step will exceed the ideal time, so
      ! return a value of 2 to indicate the ideal time step will be exceeded
      ! in 2 time steps.
      !
      return_value = 2
      !
    end if
    !
  end if
  !
end function check_against_ideal_time
!
!###############################################################################
!
subroutine check_output_frequency(dtmin,dt_sp,skip_updates)
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_solpts
  use ovar,        only : time,time_ref,time_scaling_factor
  use ovar,        only : output_frequency
  use ovar,        only : is_at_output_frequency ! modified
  use ovar,        only : restart_frequency
  use ovar,        only : is_at_restart_frequency ! modified
  use ovar,        only : output_time_averaging,ave_start_time
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: dtmin
  real(wp), intent(inout) :: dt_sp(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: skip_updates
  !
  !.. Local Scalars ..
  integer     :: n
  integer     :: old_out_int,new_out_int
  integer     :: old_rst_int,new_rst_int
  real(wp)    :: old_time,new_time
  real(wp)    :: new_out_time,new_rst_time
  logical(lk) :: is_at_output_freq
  logical(lk) :: is_at_restart_freq
  logical(lk) :: perform_updates
  !
  !.. Local Parameters ..
  logical(lk), parameter :: update_dtmin_for_output_frequency = fals
 !logical(lk), parameter :: update_dtmin_for_output_frequency = true
  !
continue
  !
  perform_updates = true
  if (present(skip_updates)) then
    perform_updates = (.not. skip_updates)
  end if
  !
  is_at_output_freq = is_at_output_frequency
  is_at_restart_freq = is_at_restart_frequency
  !
  old_time = time*time_ref*time_scaling_factor
  new_time = old_time + dtmin*time_ref*time_scaling_factor
  !
  if (output_time_averaging) then
    old_time = old_time - ave_start_time*time_ref*time_scaling_factor
    new_time = new_time - ave_start_time*time_ref*time_scaling_factor
  end if
  !
  ! Number of output frequency intervals that have
  ! been completed at the end of the previous timestep.
  !
  old_out_int = floor( old_time/output_frequency )
  !
  ! Number of output frequency intervals that will have
  ! been completed after finishing the current timestep.
  !
  new_out_int = floor( new_time/output_frequency )
  !
  ! If the number of output frequency intervals completed increases
  ! after completing the current timestep, set the logical flag that
  ! tells the main loop to output a new solution file.
  !
  if (new_out_int > old_out_int) then
    is_at_output_freq = true
  end if
  !
  ! Number of restart frequency intervals that have
  ! been completed at the end of the previous timestep.
  !
  old_rst_int = floor( old_time/restart_frequency )
  !
  ! Number of restart frequency intervals that will have
  ! been completed after finishing the current timestep.
  !
  new_rst_int = floor( new_time/restart_frequency )
  !
  ! If the number of restart frequency intervals completed increases
  ! after completing the current timestep, set the logical flag that
  ! tells the main loop to output a new restart file.
  !
  if (new_rst_int > old_rst_int) then
    is_at_restart_freq = true
  end if
  !
  ! NOTE: WITHOUT ALTERING DTMIN, WE ARE ASSUMING THAT THE TIMESTEP SIZE IS
  !       SUFFICIENTLY SMALL ENOUGH SUCH THAT THE SOLUTION TIME AT THE END OF
  !       THE CURRENT TIMESTEP IS VERY NEARLY CLOSE (OR CLOSE ENOUGH) TO A
  !       MULTPLE OF THE OUTPUT FREQUENCY.
  !
  !       SHOULD WE ALTER DTMIN (AND THUS THE DTSP ARRAY) SO THAT THE SOLUTION
  !       TIME AT THE END OF THE CURRENT TIMESTEP IS ACTUALLY AT A MULTIPLE OF
  !       THE OUTPUT FREQUENCY?
  !
  if ( (is_at_output_freq .or. is_at_restart_freq) .and. &
       (update_dtmin_for_output_frequency) ) then
    !
    new_out_time = new_time
    new_rst_time = new_time
    !
    if (is_at_output_freq) then
      new_out_time = output_frequency*real(new_out_int,kind=wp) / &
                     (time_ref*time_scaling_factor)
    end if
    !
    if (is_at_restart_freq) then
      new_rst_time = restart_frequency*real(new_rst_int,kind=wp) / &
                     (time_ref*time_scaling_factor)
    end if
    !
    ! Compute the new dtmin
    !
    if (output_time_averaging) then
      dtmin = min(new_out_time,new_rst_time) + ave_start_time - time
    else
      dtmin = min(new_out_time,new_rst_time) - time
    end if
    !$ACC UPDATE DEVICE(dtmin)
    !
    !$ACC PARALLEL LOOP PRESENT(dt_sp)
    do n = 1,n_solpts
      dt_sp(n) = dtmin
    end do
    !
  end if
  !
  if (perform_updates) then
    is_at_output_frequency = is_at_output_freq
    is_at_restart_frequency = is_at_restart_freq
  end if
  !
end subroutine check_output_frequency
!
!###############################################################################
!
subroutine check_final_time(dtmin,dt_sp,final_timestep_from_stop_request, &
                            skip_updates)
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_solpts
  use ovar,        only : itcur,time,Final_Time,time_ref,time_scaling_factor
  use ovar,        only : this_is_final_timestep ! modified
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: dtmin
  real(wp), intent(inout) :: dt_sp(:)
  integer,     intent(in) :: final_timestep_from_stop_request
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: skip_updates
  !
  !.. Local Scalars ..
  integer     :: n,ichk
  real(wp)    :: chk_time
  logical(lk) :: reset_dtsp
  logical(lk) :: this_is_last_timestep
  logical(lk) :: perform_updates
  !
continue
  !
  perform_updates = true
  if (present(skip_updates)) then
    perform_updates = (.not. skip_updates)
  end if
  !
  this_is_last_timestep = this_is_final_timestep
  !
  reset_dtsp = fals
  !
  if (Final_Time > zero) then
    !
    chk_time = (Final_Time/time_scaling_factor)/time_ref + eps9
    !
    ichk = check_against_ideal_time(time,dtmin,chk_time)
    !
    if (ichk == -1) then
      !
      this_is_last_timestep = true
      !
    else if (ichk == 1) then
      !
      this_is_last_timestep = true
      !
      reset_dtsp = true
      !
      dtmin = (Final_Time/time_scaling_factor)/time_ref - time
      !
    else if (ichk == 2) then
      !
      reset_dtsp = true
      !
      dtmin = half * ((Final_Time/time_scaling_factor)/time_ref - time)
      !
    end if
    !
  end if
  !
  ! Check to see if a final time step number was provided by a temp file
  !
  if (final_timestep_from_stop_request /= 0) then
    if (final_timestep_from_stop_request <= itcur) then
      this_is_last_timestep = true
    end if
  end if
  !
 !!$ACC UPDATE DEVICE(dtmin)
  !
  if (reset_dtsp) then
    !$ACC PARALLEL LOOP PRESENT(dt_sp)
    do n = 1,n_solpts
      dt_sp(n) = dtmin
    end do
  end if
  !
  if (perform_updates) then
    this_is_final_timestep = this_is_last_timestep
  end if
  !
end subroutine check_final_time
!
!###############################################################################
!
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module dt_mod
