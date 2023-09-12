module time_mod
  !
!#define PROFILE_ON
  !
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  use eqn_idx, only : ncs,npres,nvisc,ngamm
  !
  use ovar, only : force_nd_subroutines
  !
  implicit none
  !
  private
  !
  public :: get_dtsp
  public :: checks_on_dtmin
  public :: update_Classic_RK
  public :: update_CK_RK4
  public :: update_TVD_RK2
  public :: update_TVD_RK3
  public :: update_TVD_RK4
  !
  interface write_array
    module procedure write_1d_array, write_2d_array, write_3d_array
  end interface write_array
  !
contains
!
!###############################################################################
!
subroutine get_dtsp(dtmin,final_timestep_from_stop_request,skip_updates)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : Timestep_Type
  use ovar,    only : d_t,cfl ! modified
  use flowvar, only : dtsp ! modified
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
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_dtsp"
  !
continue
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
  ! Check to see if we need to output a new
  ! solution file based on output_frequency
  !
  if (output_frequency > zero) then
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
  use ovar,    only : constant_dt,CFL
  use geovar,  only : n_solpts
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
  use geovar, only : nfbnd,nface,n_solpts
  use geovar, only : face,cell,n_bnd_flxpts,n_flxpts
  !
  use ovar, only : CFL,governing_equations
  !
  use metrics_mod, only : facenrm,facejac
  !
  use flowvar, only : faceusp
  use flowvar, only : facegpsp
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  real(wp), intent(out) :: dt_sp(:)
  !
  !.. Local Scalars ..
  integer  :: n,kf,nf,f1,f2,n1,n2
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: ds,eigni_fp,eignv_fp
  real(wp) :: dtl,hl,lcl
  real(wp) :: dtr,hr,rcl
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: viscl(:)
  real(wp), allocatable :: viscr(:)
  real(wp), allocatable :: eigni(:,:)
  real(wp), allocatable :: dsf(:)
  real(wp), allocatable :: psq(:,:)
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
  allocate ( eigni(1:n_flxpts,1:2) )
  allocate ( dsf(1:n_flxpts) )
  allocate ( psq(1:n_flxpts,1:2) )
  !
  !$ACC DATA CREATE(viscl,viscr,eigni,dsf,psq)
  !
  n1 = 1
  n2 = n_flxpts
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
  call invs_eigenvalue(n1,n2,faceusp,facegpsp,facenrm,eigni,dsf)
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
  !$ACC&              PRESENT(facejac,dsf,eigni,viscl,psq)
  do kf = 1,n_bnd_flxpts
    !
    hl = psq(kf,1) * half * dsf(kf) / facejac(kf,1)
    !
    eigni_fp = eigni(kf,1)
    eignv_fp = hl * viscl(kf)
    !
    dtl = CFL / (hl * (eigni_fp + eignv_fp))
    !
    dtmin = min(dtmin,dtl)
    !
  end do
  !
  !$ACC PARALLEL LOOP REDUCTION(min:dtmin) &
  !$ACC&              PRESENT(facejac,dsf,eigni,viscl,viscr,psq)
  do kf = n_bnd_flxpts+1,n_flxpts
    !
    ds = half * dsf(kf)
    !
    hl = ds * psq(kf,1) / facejac(kf,1)
    hr = ds * psq(kf,2) / facejac(kf,2)
    !
    eigni_fp = max( eigni(kf,1) , eigni(kf,2) )
    eignv_fp = max( hl * viscl(kf) , hr * viscr(kf) )
    !
    dtl = CFL / (max(hl,hr) * (eigni_fp + eignv_fp))
    dtmin = min(dtmin,dtl)
   !dtl = CFL / (hl * (eigni_fp + eignv_fp))
   !dtr = CFL / (hr * (eigni_fp + eignv_fp))
   !dtmin = min(dtmin,dtl,dtr)
    !
  end do
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
  deallocate ( dsf )
  deallocate ( eigni )
  deallocate ( viscr )
  deallocate ( viscl )
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
  use geovar, only : nfbnd,nface,ncell
  use geovar, only : face,cell,n_flxpts
  !
  use ovar, only : CFL,governing_equations
  !
  use metrics_mod, only : facenrm,facejac
  !
  use flowvar, only : faceusp
  use flowvar, only : facegpsp
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
  real(wp) :: ds,dtc,eigni_fp,eignv_fp
  real(wp) :: dtl,hl,psql
  real(wp) :: dtr,hr,psqr
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: viscl(:)
  real(wp), allocatable :: viscr(:)
  real(wp), allocatable :: dt_cell(:)
  real(wp), allocatable :: eigni(:,:)
  real(wp), allocatable :: dsf(:)
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
  allocate ( dsf(1:n_flxpts) )
  !
  !$ACC DATA CREATE(viscl,viscr,eigni,dt_cell,dsf)
  !
  n1 = 1
  n2 = n_flxpts
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
  call invs_eigenvalue(n1,n2,faceusp,facegpsp,facenrm,eigni,dsf)
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
  !$ACC&              PRESENT(dt_cell,dsf,eigni,viscl)
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
      hl = dsf(kf) * psql / facejac(kf,1)
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
  !$ACC&              PRESENT(dt_cell,dsf,eigni,viscl,viscr)
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
      ds = dsf(kf)
      !
      ! Left and right cell jacobians at the face point
      !
      hl = ds * psql / facejac(kf,1)
      hr = ds * psqr / facejac(kf,2)
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
  deallocate ( dsf )
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
subroutine invs_eigenvalue(n1,n2,face_cv,face_gp,facenrm,eigni,ds)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: eigni(:,:)
  real(wp), intent(inout) :: ds(:)
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
  if (nr == 3) then      ! this should also work for ncs > 0
    call invs_eigenvalue_3d(n1,n2,face_cv,face_gp,facenrm,eigni,ds)
  else if (nr == 2) then ! this should also work for ncs > 0
    call invs_eigenvalue_2d(n1,n2,face_cv,face_gp,facenrm,eigni,ds)
  else
    call invs_eigenvalue_nd(n1,n2,face_cv,face_gp,facenrm,eigni,ds)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine invs_eigenvalue
!
!###############################################################################
!
pure subroutine invs_eigenvalue_3d(n1,n2,face_cv,face_gp,facenrm,eigni,ds)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: eigni(:,:)
  real(wp), intent(inout) :: ds(:)
  !
  !.. Local Scalars ..
  integer  :: i,n
  real(wp) :: rho_inv,pres,gamm
  real(wp) :: velx,vely,velz,vface
  real(wp) :: rnx,rny,rnz,unx,uny,unz,dl
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: aspd(:,:)
  !
continue
  !
  allocate ( aspd(n1:n2,1:2) )
  !$ACC DATA CREATE(aspd)
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
  !
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,facenrm,aspd,eigni,ds)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    rnz = facenrm(3,n)
    ds(n) = sqrt(rnx*rnx + rny*rny + rnz*rnz)
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    dl = one/ds(n)
    unx = rnx*dl
    uny = rny*dl
    unz = rnz*dl
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
    eigni(n,1) = abs(vface) + aspd(n,1)
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
    eigni(n,2) = abs(vface) + aspd(n,2)
    !
  end do
  !
  !$ACC END DATA
  deallocate ( aspd )
  !
end subroutine invs_eigenvalue_3d
!
!###############################################################################
!
pure subroutine invs_eigenvalue_2d(n1,n2,face_cv,face_gp,facenrm,eigni,ds)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: eigni(:,:)
  real(wp), intent(inout) :: ds(:)
  !
  !.. Local Scalars ..
  integer  :: i,n
  real(wp) :: rho_inv,pres,gamm
  real(wp) :: velx,vely,vface
  real(wp) :: rnx,rny,unx,uny,dl
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: aspd(:,:)
  !
continue
  !
  allocate ( aspd(n1:n2,1:2) )
  !$ACC DATA CREATE(aspd)
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
  !$ACC PARALLEL LOOP PRESENT(face_cv,facenrm,aspd,eigni,ds)
  do n = n1,n2
    !
    rnx = facenrm(1,n)
    rny = facenrm(2,n)
    ds(n) = sqrt(rnx*rnx + rny*rny)
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    dl = one/ds(n)
    unx = rnx*dl
    uny = rny*dl
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
    eigni(n,1) = abs(vface) + aspd(n,1)
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
    eigni(n,2) = abs(vface) + aspd(n,2)
    !
  end do
  !
  !$ACC END DATA
  deallocate ( aspd )
  !
end subroutine invs_eigenvalue_2d
!
!###############################################################################
!
pure subroutine invs_eigenvalue_nd(n1,n2,face_cv,face_gp,facenrm,eigni,ds)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: face_gp(:,:,:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: eigni(:,:)
  real(wp), intent(inout) :: ds(:)
  !
  !.. Local Scalars ..
  integer  :: i,n
  real(wp) :: rho_inv,pres,gamm,vface
  !
  !.. Local Arrays ..
  real(wp) :: unrm(1:size(facenrm,dim=1))
  real(wp) :: vel(nmb:nme)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: aspd(:,:)
  !
continue
  !
  allocate ( aspd(n1:n2,1:2) )
  !$ACC DATA CREATE(aspd,unrm,vel)
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
  !$ACC PARALLEL LOOP PRESENT(face_cv,facenrm,aspd,eigni,ds,unrm,vel)
  do n = n1,n2
    !
    ds(n) = norm2(facenrm(:,n))
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    unrm(:) = facenrm(:,n)/ds(n)
#endif
    !
    ! Left side of face
    !
    vel(nmb:nme) = face_cv(n,nmb:nme,1)/face_cv(n,nec,1)
    !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    vface = dot_product( unrm(:) , vel(nmb:nme) )
#else
    vface = norm2(vel(nmb:nme))
#endif
    !
    eigni(n,1) = abs(vface) + aspd(n,1)
    !
    ! Right side of face
    !
    vel(nmb:nme) = face_cv(n,nmb:nme,2)/face_cv(n,nec,2)
    !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    vface = dot_product( facenrm(:,n) , vel(nmb:nme) )
#else
    vface = norm2(vel(nmb:nme))
#endif
    !
    eigni(n,2) = abs(vface) + aspd(n,2)
    !
  end do
  !
  !$ACC END DATA
  deallocate ( aspd )
  !
end subroutine invs_eigenvalue_nd
!
!###############################################################################
!
subroutine get_point_global_timestep(dtmin,dt_sp)
  !
  !.. Use Statements ..
  use geovar,      only : nr,n_solpts
  use ovar,        only : gam,CFL
  use metrics_mod, only : time_metrics
  use flowvar,     only : usp
  use flowvar,     only : gpsp
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
    !$ACC&              PRESENT(usp,aspd,time_metrics,dt_sp)
    do n = 1,n_solpts
      !
      dx = time_metrics(1,n)
      dy = time_metrics(2,n)
      dz = time_metrics(3,n)
      da = time_metrics(4,n)
      !
      rho_inv = one/usp(n,nec)
      !
      velx = usp(n,nmx)*rho_inv
      vely = usp(n,nmy)*rho_inv
      velz = usp(n,nmz)*rho_inv
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
  use geovar,      only : nr,n_solpts
  use ovar,        only : gam,CFL
  use metrics_mod, only : time_metrics
  use flowvar,     only : usp
  use flowvar,     only : gpsp
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
    !$ACC&              PRESENT(usp,aspd,time_metrics,dt_sp)
    do n = 1,n_solpts
      !
      dx = time_metrics(1,n)
      dy = time_metrics(2,n)
      dz = time_metrics(3,n)
      da = time_metrics(4,n)
      !
      rho_inv = one/usp(n,nec)
      !
      velx = usp(n,nmx)*rho_inv
      vely = usp(n,nmy)*rho_inv
      velz = usp(n,nmz)*rho_inv
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
  use geovar,  only : n_solpts
  use ovar,    only : time,time_ref
  use ovar,    only : write_TaylorGreen_8s_solution ! modified
  use ovar,    only : write_TaylorGreen_full_solution ! modified
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
  use geovar,  only : n_solpts
  use ovar,    only : time,time_ref
  use ovar,    only : write_TaylorGreen_full_solution ! modified
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
  use geovar,  only : n_solpts
  use ovar,    only : time,time_ref,time_scaling_factor
  use ovar,    only : output_frequency
  use ovar,    only : is_at_output_frequency ! modified
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: dtmin
  real(wp), intent(inout) :: dt_sp(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: skip_updates
  !
  !.. Local Scalars ..
  integer     :: n,old_int,new_int
  real(wp)    :: old_time,new_time
  logical(lk) :: is_at_output_freq
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
  !
  old_time = time*time_ref*time_scaling_factor
  new_time = old_time + dtmin*time_ref*time_scaling_factor
  !
  ! Number of output frequency intervals that have
  ! been completed at the end of the previous timestep.
  !
  old_int = floor( old_time/output_frequency )
  !
  ! Number of output frequency intervals that will have
  ! been completed after finishing the current timestep.
  !
  new_int = floor( new_time/output_frequency )
  !
  ! If the number of output frequency intervals completed increases
  ! after completing the current timestep, set the logical flag that
  ! tells the main loop to output a new solution file.
  !
  if (new_int > old_int) then
    is_at_output_freq = true
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
  if (update_dtmin_for_output_frequency) then
    !
    if (is_at_output_freq) then
      !
      new_time = output_frequency*real(new_int,kind=wp)
      !
      ! Nondimensionalize new_time to match time and dtmin
      !
      dtmin = time - new_time/time_ref/time_scaling_factor
      !$ACC UPDATE DEVICE(dtmin)
      !
      !$ACC PARALLEL LOOP PRESENT(dt_sp)
      do n = 1,n_solpts
        dt_sp(n) = dtmin
      end do
      !
    end if
    !
  end if
  !
  if (perform_updates) then
    is_at_output_frequency = is_at_output_freq
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
  use geovar,  only : n_solpts
  use ovar,    only : itcur,time,Final_Time,time_ref
  use ovar,    only : this_is_final_timestep ! modified
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
    chk_time = Final_Time/time_ref + eps9
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
      dtmin = Final_Time/time_ref - time
      !
    else if (ichk == 2) then
      !
      reset_dtsp = true
      !
      dtmin = half * (Final_Time/time_ref - time)
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
subroutine update_CK_RK4(nst,Apply_Filter,Apply_Limiter)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : rk_sum,residual,usp,dtsp
  use rhs_mod, only : stabilize_usp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_CK_RK4"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nst == 1) then
    call update_CK_RK4_S1(n_solpts,dtsp,residual,usp,rk_sum)
  else if (nst == 2) then
    call update_CK_RK4_S2(n_solpts,dtsp,residual,usp,rk_sum)
  else if (nst == 3) then
    call update_CK_RK4_S3(n_solpts,dtsp,residual,usp,rk_sum)
  else if (nst == 4) then
    call update_CK_RK4_S4(n_solpts,dtsp,residual,usp,rk_sum)
  else if (nst == 5) then
    call update_CK_RK4_S5(n_solpts,dtsp,residual,usp,rk_sum)
  end if
  !
  ! Apply filtering or limiting to stabilize the solution, and then update
  ! the thermodynamic and transport properties at the solution points
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter,caller_name=pname)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_CK_RK4
!
!###############################################################################
!
pure subroutine update_CK_RK4_S1(np,dtsp,residual,usp,rk_sum)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rk_sum(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: sum_ec,sum_mx,sum_my,sum_mz,sum_ee,sum_eq
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_CK_RK4_S1"
  !
  real(wp), parameter :: b = 1432997174477.0_wp / 9575080441755.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = dtsp(n)*residual(n,nec)
      sum_mx = dtsp(n)*residual(n,nmx)
      sum_my = dtsp(n)*residual(n,nmy)
      sum_mz = dtsp(n)*residual(n,nmz)
      sum_ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nmz) = usp(n,nmz) + b*sum_mz
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nmz) = sum_mz
      rk_sum(n,nee) = sum_ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = dtsp(n)*residual(n,nec)
      sum_mx = dtsp(n)*residual(n,nmx)
      sum_my = dtsp(n)*residual(n,nmy)
      sum_ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nee) = sum_ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,rk_sum,dtsp,residual)
    do m = 1,nq
      do n = 1,np
        !
        sum_eq = dtsp(n)*residual(n,m)
        !
        usp(n,m) = usp(n,m) + b*sum_eq
        !
        rk_sum(n,m) = sum_eq
        !
      end do
    end do
  end if
  !
end subroutine update_CK_RK4_S1
!
!###############################################################################
!
pure subroutine update_CK_RK4_S2(np,dtsp,residual,usp,rk_sum)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rk_sum(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: sum_ec,sum_mx,sum_my,sum_mz,sum_ee,sum_eq
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_CK_RK4_S2"
  !
  real(wp), parameter :: a = -567301805773.0_wp / 1357537059087.0_wp
  real(wp), parameter :: b = 5161836677717.0_wp / 13612068292357.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_my = a*rk_sum(n,nmy) + dtsp(n)*residual(n,nmy)
      sum_mz = a*rk_sum(n,nmz) + dtsp(n)*residual(n,nmz)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nmz) = usp(n,nmz) + b*sum_mz
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nmz) = sum_mz
      rk_sum(n,nee) = sum_ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_my = a*rk_sum(n,nmy) + dtsp(n)*residual(n,nmy)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nee) = sum_ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,rk_sum,dtsp,residual)
    do m = 1,nq
      do n = 1,np
        !
        sum_eq = a*rk_sum(n,m) + dtsp(n)*residual(n,m)
        !
        usp(n,m) = usp(n,m) + b*sum_eq
        !
        rk_sum(n,m) = sum_eq
        !
      end do
    end do
  end if
  !
end subroutine update_CK_RK4_S2
!
!###############################################################################
!
pure subroutine update_CK_RK4_S3(np,dtsp,residual,usp,rk_sum)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rk_sum(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: sum_ec,sum_mx,sum_my,sum_mz,sum_ee,sum_eq
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_CK_RK4_S3"
  !
  real(wp), parameter :: a = -2404267990393.0_wp / 2016746695238.0_wp
  real(wp), parameter :: b = 1720146321549.0_wp / 2090206949498.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_my = a*rk_sum(n,nmy) + dtsp(n)*residual(n,nmy)
      sum_mz = a*rk_sum(n,nmz) + dtsp(n)*residual(n,nmz)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nmz) = usp(n,nmz) + b*sum_mz
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nmz) = sum_mz
      rk_sum(n,nee) = sum_ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_my = a*rk_sum(n,nmy) + dtsp(n)*residual(n,nmy)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nee) = sum_ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,rk_sum,dtsp,residual)
    do m = 1,nq
      do n = 1,np
        !
        sum_eq = a*rk_sum(n,m) + dtsp(n)*residual(n,m)
        !
        usp(n,m) = usp(n,m) + b*sum_eq
        !
        rk_sum(n,m) = sum_eq
        !
      end do
    end do
  end if
  !
end subroutine update_CK_RK4_S3
!
!###############################################################################
!
pure subroutine update_CK_RK4_S4(np,dtsp,residual,usp,rk_sum)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rk_sum(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: sum_ec,sum_mx,sum_my,sum_mz,sum_ee,sum_eq
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_CK_RK4_S4"
  !
  real(wp), parameter :: a = -3550918686646.0_wp / 2091501179385.0_wp
  real(wp), parameter :: b = 3134564353537.0_wp / 4481467310338.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_my = a*rk_sum(n,nmy) + dtsp(n)*residual(n,nmy)
      sum_mz = a*rk_sum(n,nmz) + dtsp(n)*residual(n,nmz)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nmz) = usp(n,nmz) + b*sum_mz
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nmz) = sum_mz
      rk_sum(n,nee) = sum_ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_my = a*rk_sum(n,nmy) + dtsp(n)*residual(n,nmy)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nee) = sum_ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,rk_sum,dtsp,residual)
    do m = 1,nq
      do n = 1,np
        !
        sum_eq = a*rk_sum(n,m) + dtsp(n)*residual(n,m)
        !
        usp(n,m) = usp(n,m) + b*sum_eq
        !
        rk_sum(n,m) = sum_eq
        !
      end do
    end do
  end if
  !
end subroutine update_CK_RK4_S4
!
!###############################################################################
!
pure subroutine update_CK_RK4_S5(np,dtsp,residual,usp,rk_sum)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rk_sum(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: sum_ec,sum_mx,sum_my,sum_mz,sum_ee,sum_eq
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_CK_RK4_S5"
  !
  real(wp), parameter :: a = -1275806237668.0_wp / 842570457699.0_wp
  real(wp), parameter :: b = 2277821191437.0_wp / 14882151754819.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_my = a*rk_sum(n,nmy) + dtsp(n)*residual(n,nmy)
      sum_mz = a*rk_sum(n,nmz) + dtsp(n)*residual(n,nmz)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nmz) = usp(n,nmz) + b*sum_mz
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nmz) = sum_mz
      rk_sum(n,nee) = sum_ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_my = a*rk_sum(n,nmy) + dtsp(n)*residual(n,nmy)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nee) = sum_ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,rk_sum,dtsp,residual)
    do m = 1,nq
      do n = 1,np
        !
        sum_eq = a*rk_sum(n,m) + dtsp(n)*residual(n,m)
        !
        usp(n,m) = usp(n,m) + b*sum_eq
        !
        rk_sum(n,m) = sum_eq
        !
      end do
    end do
  end if
  !
end subroutine update_CK_RK4_S5
!
!###############################################################################
!
subroutine update_TVD_RK2(nst,Apply_Filter,Apply_Limiter)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : residual,usp,uoldsp,dtsp
  use rhs_mod, only : stabilize_usp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK2"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nst == 1) then
    call update_TVD_RK2_S1(n_solpts,dtsp,residual,usp)
  else
    call update_TVD_RK2_S2(n_solpts,dtsp,uoldsp,residual,usp)
  end if
  !
  ! Apply filtering or limiting to stabilize the solution, and then update
  ! the thermodynamic and transport properties at the solution points
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter,caller_name=pname)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_TVD_RK2
!
!###############################################################################
!
pure subroutine update_TVD_RK2_S1(np,dtsp,residual,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK2_S1"
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,usp)
    do n = 1,np
      usp(n,nec) = usp(n,nec) + dtsp(n)*residual(n,nec)
      usp(n,nmx) = usp(n,nmx) + dtsp(n)*residual(n,nmx)
      usp(n,nmy) = usp(n,nmy) + dtsp(n)*residual(n,nmy)
      usp(n,nmz) = usp(n,nmz) + dtsp(n)*residual(n,nmz)
      usp(n,nee) = usp(n,nee) + dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,usp)
    do n = 1,np
      usp(n,nec) = usp(n,nec) + dtsp(n)*residual(n,nec)
      usp(n,nmx) = usp(n,nmx) + dtsp(n)*residual(n,nmx)
      usp(n,nmy) = usp(n,nmy) + dtsp(n)*residual(n,nmy)
      usp(n,nee) = usp(n,nee) + dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(residual,dtsp,usp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = usp(n,m) + dtsp(n)*residual(n,m)
      end do
    end do
  end if
  !
end subroutine update_TVD_RK2_S1
!
!###############################################################################
!
pure subroutine update_TVD_RK2_S2(np,dtsp,uoldsp,residual,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK2_S2"
  !
  real(wp), parameter :: a = half
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a * ( uoldsp(n,nec) + usp(n,nec) + dtsp(n)*residual(n,nec) )
      usp(n,nmx) = a * ( uoldsp(n,nmx) + usp(n,nmx) + dtsp(n)*residual(n,nmx) )
      usp(n,nmy) = a * ( uoldsp(n,nmy) + usp(n,nmy) + dtsp(n)*residual(n,nmy) )
      usp(n,nmz) = a * ( uoldsp(n,nmz) + usp(n,nmz) + dtsp(n)*residual(n,nmz) )
      usp(n,nee) = a * ( uoldsp(n,nee) + usp(n,nee) + dtsp(n)*residual(n,nee) )
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a * ( uoldsp(n,nec) + usp(n,nec) + dtsp(n)*residual(n,nec) )
      usp(n,nmx) = a * ( uoldsp(n,nmx) + usp(n,nmx) + dtsp(n)*residual(n,nmx) )
      usp(n,nmy) = a * ( uoldsp(n,nmy) + usp(n,nmy) + dtsp(n)*residual(n,nmy) )
      usp(n,nee) = a * ( uoldsp(n,nee) + usp(n,nee) + dtsp(n)*residual(n,nee) )
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(residual,dtsp,uoldsp,usp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = a * ( uoldsp(n,m) + usp(n,m) + dtsp(n)*residual(n,m) )
      end do
    end do
  end if
  !
end subroutine update_TVD_RK2_S2
!
!###############################################################################
!
subroutine update_TVD_RK3(nst,Apply_Filter,Apply_Limiter,iter)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : residual,usp,uoldsp,dtsp
  use rhs_mod, only : stabilize_usp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst,iter
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK3"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nst == 1) then
    call update_TVD_RK3_S1(n_solpts,dtsp,residual,usp)
  else if (nst == 2) then
    call update_TVD_RK3_S2(n_solpts,dtsp,uoldsp,residual,usp)
  else if (nst == 3) then
    call update_TVD_RK3_S3(n_solpts,dtsp,uoldsp,residual,usp)
  end if
  !
  ! Apply filtering or limiting to stabilize the solution, and then update
  ! the thermodynamic and transport properties at the solution points
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter,caller_name=pname)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_TVD_RK3
!
!###############################################################################
!
pure subroutine update_TVD_RK3_S1(np,dtsp,residual,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK3_S1"
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,usp)
    do n = 1,np
      usp(n,nec) = usp(n,nec) + dtsp(n)*residual(n,nec)
      usp(n,nmx) = usp(n,nmx) + dtsp(n)*residual(n,nmx)
      usp(n,nmy) = usp(n,nmy) + dtsp(n)*residual(n,nmy)
      usp(n,nmz) = usp(n,nmz) + dtsp(n)*residual(n,nmz)
      usp(n,nee) = usp(n,nee) + dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,usp)
    do n = 1,np
      usp(n,nec) = usp(n,nec) + dtsp(n)*residual(n,nec)
      usp(n,nmx) = usp(n,nmx) + dtsp(n)*residual(n,nmx)
      usp(n,nmy) = usp(n,nmy) + dtsp(n)*residual(n,nmy)
      usp(n,nee) = usp(n,nee) + dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(residual,dtsp,usp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = usp(n,m) + dtsp(n)*residual(n,m)
      end do
    end do
  end if
  !
end subroutine update_TVD_RK3_S1
!
!###############################################################################
!
pure subroutine update_TVD_RK3_S2(np,dtsp,uoldsp,residual,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK3_S2"
  !
  real(wp), parameter :: a = three4
  real(wp), parameter :: b = one4
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + b*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + b*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + b*dtsp(n)*residual(n,nmy)
      usp(n,nmz) = a*uoldsp(n,nmz) + b*usp(n,nmz) + b*dtsp(n)*residual(n,nmz)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + b*dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + b*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + b*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + b*dtsp(n)*residual(n,nmy)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + b*dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(residual,dtsp,uoldsp,usp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = a*uoldsp(n,m) + b*usp(n,m) + b*dtsp(n)*residual(n,m)
      end do
    end do
  end if
  !
end subroutine update_TVD_RK3_S2
!
!###############################################################################
!
pure subroutine update_TVD_RK3_S3(np,dtsp,uoldsp,residual,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK3_S3"
  !
  real(wp), parameter :: a = one3
  real(wp), parameter :: b = two3
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + b*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + b*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + b*dtsp(n)*residual(n,nmy)
      usp(n,nmz) = a*uoldsp(n,nmz) + b*usp(n,nmz) + b*dtsp(n)*residual(n,nmz)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + b*dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + b*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + b*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + b*dtsp(n)*residual(n,nmy)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + b*dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(residual,dtsp,uoldsp,usp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = a*uoldsp(n,m) + b*usp(n,m) + b*dtsp(n)*residual(n,m)
      end do
    end do
  end if
  !
end subroutine update_TVD_RK3_S3
!
!###############################################################################
!
subroutine update_TVD_RK4(nst,Apply_Filter,Apply_Limiter)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : residual,usp,uoldsp,rk_sum,dtsp
  use rhs_mod, only : stabilize_usp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK4"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nst == 1) then
    call update_TVD_RK4_S1(n_solpts,dtsp,uoldsp,residual,usp)
  else if (nst == 2) then
    call update_TVD_RK4_S2(n_solpts,dtsp,uoldsp,residual,usp)
  else if (nst == 3) then
    call update_TVD_RK4_S3(n_solpts,dtsp,uoldsp,residual,usp,rk_sum)
  else if (nst == 4) then
    call update_TVD_RK4_S4(n_solpts,dtsp,uoldsp,residual,usp,rk_sum)
  else if (nst == 5) then
    call update_TVD_RK4_S5(n_solpts,dtsp,uoldsp,residual,usp,rk_sum)
  end if
  !
  ! Apply filtering or limiting to stabilize the solution, and then update
  ! the thermodynamic and transport properties at the solution points
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter,caller_name=pname)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_TVD_RK4
!
!###############################################################################
!
pure subroutine update_TVD_RK4_S1(np,dtsp,uoldsp,residual,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK4_S1"
  !
  real(wp), parameter :: c = 0.39175222700395_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = uoldsp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = uoldsp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nmz) = uoldsp(n,nmz) + c*dtsp(n)*residual(n,nmz)
      usp(n,nee) = uoldsp(n,nee) + c*dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = uoldsp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = uoldsp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nee) = uoldsp(n,nee) + c*dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,dtsp,residual)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + c*dtsp(n)*residual(n,m)
      end do
    end do
  end if
  !
end subroutine update_TVD_RK4_S1
!
!###############################################################################
!
pure subroutine update_TVD_RK4_S2(np,dtsp,uoldsp,residual,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK4_S2"
  !
  real(wp), parameter :: a = 0.44437049406734_wp
  real(wp), parameter :: b = 0.55562950593266_wp
  real(wp), parameter :: c  = 0.36841059262959_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nmz) = a*uoldsp(n,nmz) + b*usp(n,nmz) + c*dtsp(n)*residual(n,nmz)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,dtsp,residual)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = a*uoldsp(n,m) + b*usp(n,m) + c*dtsp(n)*residual(n,m)
      end do
    end do
  end if
  !
end subroutine update_TVD_RK4_S2
!
!###############################################################################
!
pure subroutine update_TVD_RK4_S3(np,dtsp,uoldsp,residual,usp,rk_sum)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rk_sum(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: stage2_ec,stage2_mx,stage2_my,stage2_mz,stage2_ee,stage2_eq
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK4_S3"
  !
  real(wp), parameter :: a = 0.62010185138540_wp
  real(wp), parameter :: b = 0.37989814861460_wp
  real(wp), parameter :: c = 0.25189177424738_wp
  real(wp), parameter :: d = 0.51723167208978_wp
  real(wp), parameter :: e = 0.12759831133288_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual,rk_sum)
    do n = 1,np
      !
      stage2_ec = usp(n,nec)
      stage2_mx = usp(n,nmx)
      stage2_my = usp(n,nmy)
      stage2_mz = usp(n,nmz)
      stage2_ee = usp(n,nee)
      !
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nmz) = a*uoldsp(n,nmz) + b*usp(n,nmz) + c*dtsp(n)*residual(n,nmz)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*dtsp(n)*residual(n,nee)
      !
      rk_sum(n,nec) = d*stage2_ec + e*usp(n,nec)
      rk_sum(n,nmx) = d*stage2_mx + e*usp(n,nmx)
      rk_sum(n,nmy) = d*stage2_my + e*usp(n,nmy)
      rk_sum(n,nmz) = d*stage2_mz + e*usp(n,nmz)
      rk_sum(n,nee) = d*stage2_ee + e*usp(n,nee)
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual,rk_sum)
    do n = 1,np
      !
      stage2_ec = usp(n,nec)
      stage2_mx = usp(n,nmx)
      stage2_my = usp(n,nmy)
      stage2_ee = usp(n,nee)
      !
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*dtsp(n)*residual(n,nee)
      !
      rk_sum(n,nec) = d*stage2_ec + e*usp(n,nec)
      rk_sum(n,nmx) = d*stage2_mx + e*usp(n,nmx)
      rk_sum(n,nmy) = d*stage2_my + e*usp(n,nmy)
      rk_sum(n,nee) = d*stage2_ee + e*usp(n,nee)
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,dtsp,residual,rk_sum)
    do m = 1,nq
      do n = 1,np
        !
        stage2_eq = usp(n,m)
        !
        usp(n,m) = a*uoldsp(n,m) + b*usp(n,m) + c*dtsp(n)*residual(n,m)
        !
        rk_sum(n,m) = d*stage2_eq + e*usp(n,m)
        !
      end do
    end do
  end if
  !
end subroutine update_TVD_RK4_S3
!
!###############################################################################
!
pure subroutine update_TVD_RK4_S4(np,dtsp,uoldsp,residual,usp,rk_sum)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rk_sum(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK4_S4"
  !
  real(wp), parameter :: a = 0.17807995410773_wp
  real(wp), parameter :: b = 0.82192004589227_wp
  real(wp), parameter :: c = 0.54497475021237_wp
  real(wp), parameter :: d = 0.08460416338212_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual,rk_sum)
    do n = 1,np
      !
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nmz) = a*uoldsp(n,nmz) + b*usp(n,nmz) + c*dtsp(n)*residual(n,nmz)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*dtsp(n)*residual(n,nee)
      !
      rk_sum(n,nec) = rk_sum(n,nec) + d*dtsp(n)*residual(n,nec)
      rk_sum(n,nmx) = rk_sum(n,nmx) + d*dtsp(n)*residual(n,nmx)
      rk_sum(n,nmy) = rk_sum(n,nmy) + d*dtsp(n)*residual(n,nmy)
      rk_sum(n,nmz) = rk_sum(n,nmz) + d*dtsp(n)*residual(n,nmz)
      rk_sum(n,nee) = rk_sum(n,nee) + d*dtsp(n)*residual(n,nee)
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual,rk_sum)
    do n = 1,np
      !
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*dtsp(n)*residual(n,nee)
      !
      rk_sum(n,nec) = rk_sum(n,nec) + d*dtsp(n)*residual(n,nec)
      rk_sum(n,nmx) = rk_sum(n,nmx) + d*dtsp(n)*residual(n,nmx)
      rk_sum(n,nmy) = rk_sum(n,nmy) + d*dtsp(n)*residual(n,nmy)
      rk_sum(n,nee) = rk_sum(n,nee) + d*dtsp(n)*residual(n,nee)
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,dtsp,residual,rk_sum)
    do m = 1,nq
      do n = 1,np
        !
        usp(n,m) = a*uoldsp(n,m) + b*usp(n,m) + c*dtsp(n)*residual(n,m)
        !
        rk_sum(n,m) = rk_sum(n,m) + d*dtsp(n)*residual(n,m)
        !
      end do
    end do
  end if
  !
end subroutine update_TVD_RK4_S4
!
!###############################################################################
!
pure subroutine update_TVD_RK4_S5(np,dtsp,uoldsp,res,usp,rk_sum)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: res(:,:)
  real(wp),    intent(in) :: rk_sum(:,:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: c
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK4_S5"
  !
  real(wp), parameter :: a = 0.00683325884039_wp
  real(wp), parameter :: b = 0.34833675773694_wp
  real(wp), parameter :: d = 0.22600748319395_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,res,rk_sum)
    do n = 1,np
      !
      c = d*dtsp(n)
      !
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*res(n,nec) + rk_sum(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*res(n,nmx) + rk_sum(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*res(n,nmy) + rk_sum(n,nmy)
      usp(n,nmz) = a*uoldsp(n,nmz) + b*usp(n,nmz) + c*res(n,nmz) + rk_sum(n,nmz)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*res(n,nee) + rk_sum(n,nee)
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,res,rk_sum)
    do n = 1,np
      !
      c = d*dtsp(n)
      !
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*res(n,nec) + rk_sum(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*res(n,nmx) + rk_sum(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*res(n,nmy) + rk_sum(n,nmy)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*res(n,nee) + rk_sum(n,nee)
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,dtsp,res,rk_sum)
    do m = 1,nq
      do n = 1,np
        !
        c = d*dtsp(n)
        !
        usp(n,m) = a*uoldsp(n,m) + b*usp(n,m) + c*res(n,m) + rk_sum(n,m)
        !
      end do
    end do
  end if
  !
end subroutine update_TVD_RK4_S5
!
!###############################################################################
!
subroutine update_Classic_RK(nst,Apply_Filter,Apply_Limiter)
  !
  ! First, store d(utd)/dt for runge kutta time
  ! stepping then update usp
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : num_rk_stages
  use flowvar, only : rssprk,residual,usp,uoldsp,dtsp
  use rhs_mod, only : stabilize_usp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_Classic_RK"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  ! Store the residual for this Runge-Kutta step in rssprk
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,rssprk)
    do n = 1,n_solpts
      rssprk(n,nec,nst) = residual(n,nec)
      rssprk(n,nmx,nst) = residual(n,nmx)
      rssprk(n,nmy,nst) = residual(n,nmy)
      rssprk(n,nmz,nst) = residual(n,nmz)
      rssprk(n,nee,nst) = residual(n,nee)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,rssprk)
    do n = 1,n_solpts
      rssprk(n,nec,nst) = residual(n,nec)
      rssprk(n,nmx,nst) = residual(n,nmx)
      rssprk(n,nmy,nst) = residual(n,nmy)
      rssprk(n,nee,nst) = residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(residual,rssprk)
    do m = 1,nq
      do n = 1,n_solpts
        rssprk(n,m,nst) = residual(n,m)
      end do
    end do
  end if
  !
  ! Update usp
  !
  if (num_rk_stages == 1) then
    call update_Classic_RK1_S1(n_solpts,uoldsp,rssprk,dtsp,usp)
  else if (num_rk_stages == 2) then
    if (nst == 1) then
      call update_Classic_RK2_S1(n_solpts,uoldsp,rssprk,dtsp,usp)
    else if (nst == 2) then
      call update_Classic_RK2_S2(n_solpts,uoldsp,rssprk,dtsp,usp)
    end if
  else if (num_rk_stages == 3) then
    if (nst == 1) then
      call update_Classic_RK3_S1(n_solpts,uoldsp,rssprk,dtsp,usp)
    else if (nst == 2) then
      call update_Classic_RK3_S2(n_solpts,uoldsp,rssprk,dtsp,usp)
    else if (nst == 3) then
      call update_Classic_RK3_S3(n_solpts,uoldsp,rssprk,dtsp,usp)
    end if
  else if (num_rk_stages == 3) then
    if (nst == 1) then
      call update_Classic_RK4_S1(n_solpts,uoldsp,rssprk,dtsp,usp)
    else if (nst == 2) then
      call update_Classic_RK4_S2(n_solpts,uoldsp,rssprk,dtsp,usp)
    else if (nst == 3) then
      call update_Classic_RK4_S3(n_solpts,uoldsp,rssprk,dtsp,usp)
    else if (nst == 4) then
      call update_Classic_RK4_S4(n_solpts,uoldsp,rssprk,dtsp,usp)
    end if
  end if
  !
  ! Apply filtering or limiting to stabilize the solution, and then update
  ! the thermodynamic and transport properties at the solution points
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter,caller_name=pname)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_Classic_RK
!
!###############################################################################
!
pure subroutine update_Classic_RK1_S1(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,1)
      usp(n,nmz) = uoldsp(n,nmz) + dtsp(n)*rssprk(n,nmz,1)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,1)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,1)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,1)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + dtsp(n)*rssprk(n,m,1)
      end do
    end do
  end if
  !
end subroutine update_Classic_RK1_S1
!
!###############################################################################
!
pure subroutine update_Classic_RK2_S1(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nmz) = uoldsp(n,nmz) + half*dtsp(n)*rssprk(n,nmz,1)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,1)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,1)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + half*dtsp(n)*rssprk(n,m,1)
      end do
    end do
  end if
  !
end subroutine update_Classic_RK2_S1
!
!###############################################################################
!
pure subroutine update_Classic_RK2_S2(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,2)
      usp(n,nmz) = uoldsp(n,nmz) + dtsp(n)*rssprk(n,nmz,2)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,2)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,2)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,2)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + dtsp(n)*rssprk(n,m,2)
      end do
    end do
  end if
  !
end subroutine update_Classic_RK2_S2
!
!###############################################################################
!
pure subroutine update_Classic_RK3_S1(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + one3*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + one3*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + one3*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nmz) = uoldsp(n,nmz) + one3*dtsp(n)*rssprk(n,nmz,1)
      usp(n,nee) = uoldsp(n,nee) + one3*dtsp(n)*rssprk(n,nee,1)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + one3*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + one3*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + one3*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nee) = uoldsp(n,nee) + one3*dtsp(n)*rssprk(n,nee,1)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + one3*dtsp(n)*rssprk(n,m,1)
      end do
    end do
  end if
  !
end subroutine update_Classic_RK3_S1
!
!###############################################################################
!
pure subroutine update_Classic_RK3_S2(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + two3*dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + two3*dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + two3*dtsp(n)*rssprk(n,nmy,2)
      usp(n,nmz) = uoldsp(n,nmz) + two3*dtsp(n)*rssprk(n,nmz,2)
      usp(n,nee) = uoldsp(n,nee) + two3*dtsp(n)*rssprk(n,nee,2)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + two3*dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + two3*dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + two3*dtsp(n)*rssprk(n,nmy,2)
      usp(n,nee) = uoldsp(n,nee) + two3*dtsp(n)*rssprk(n,nee,2)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + two3*dtsp(n)*rssprk(n,m,2)
      end do
    end do
  end if
  !
end subroutine update_Classic_RK3_S2
!
!###############################################################################
!
pure subroutine update_Classic_RK3_S3(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: c1,c2
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      !
      c1 = one4*dtsp(n)
      c2 = three4*dtsp(n)
      !
      usp(n,nec) = uoldsp(n,nec) + c1*rssprk(n,nec,1) + c2*rssprk(n,nec,3)
      usp(n,nmx) = uoldsp(n,nmx) + c1*rssprk(n,nmx,1) + c2*rssprk(n,nmx,3)
      usp(n,nmy) = uoldsp(n,nmy) + c1*rssprk(n,nmy,1) + c2*rssprk(n,nmy,3)
      usp(n,nmz) = uoldsp(n,nmz) + c1*rssprk(n,nmz,1) + c2*rssprk(n,nmz,3)
      usp(n,nee) = uoldsp(n,nee) + c1*rssprk(n,nee,1) + c2*rssprk(n,nee,3)
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      !
      c1 = one4*dtsp(n)
      c2 = three4*dtsp(n)
      !
      usp(n,nec) = uoldsp(n,nec) + c1*rssprk(n,nec,1) + c2*rssprk(n,nec,3)
      usp(n,nmx) = uoldsp(n,nmx) + c1*rssprk(n,nmx,1) + c2*rssprk(n,nmx,3)
      usp(n,nmy) = uoldsp(n,nmy) + c1*rssprk(n,nmy,1) + c2*rssprk(n,nmy,3)
      usp(n,nee) = uoldsp(n,nee) + c1*rssprk(n,nee,1) + c2*rssprk(n,nee,3)
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        c1 = one4*dtsp(n)
        c2 = three4*dtsp(n)
        !
        usp(n,m) = uoldsp(n,m) + c1*rssprk(n,m,1) + c2*rssprk(n,m,3)
        !
      end do
    end do
  end if
  !
end subroutine update_Classic_RK3_S3
!
!###############################################################################
!
pure subroutine update_Classic_RK4_S1(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nmz) = uoldsp(n,nmz) + half*dtsp(n)*rssprk(n,nmz,1)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,1)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,1)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + half*dtsp(n)*rssprk(n,m,1)
      end do
    end do
  end if
  !
end subroutine update_Classic_RK4_S1
!
!###############################################################################
!
pure subroutine update_Classic_RK4_S2(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,2)
      usp(n,nmz) = uoldsp(n,nmz) + half*dtsp(n)*rssprk(n,nmz,2)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,2)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,2)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,2)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + half*dtsp(n)*rssprk(n,m,2)
      end do
    end do
  end if
  !
end subroutine update_Classic_RK4_S2
!
!###############################################################################
!
pure subroutine update_Classic_RK4_S3(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,3)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,3)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,3)
      usp(n,nmz) = uoldsp(n,nmz) + dtsp(n)*rssprk(n,nmz,3)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,3)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,3)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,3)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,3)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,3)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + dtsp(n)*rssprk(n,m,3)
      end do
    end do
  end if
  !
end subroutine update_Classic_RK4_S3
!
!###############################################################################
!
pure subroutine update_Classic_RK4_S4(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: c1,c2
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      !
      c1 = one6*dtsp(n)
      c2 = one3*dtsp(n)
      !
      usp(n,nec) = uoldsp(n,nec) + c1*rssprk(n,nec,1) + c2*rssprk(n,nec,2) &
                                 + c2*rssprk(n,nec,3) + c1*rssprk(n,nec,4)
      usp(n,nmx) = uoldsp(n,nmx) + c1*rssprk(n,nmx,1) + c2*rssprk(n,nmx,2) &
                                 + c2*rssprk(n,nmx,3) + c1*rssprk(n,nmx,4)
      usp(n,nmy) = uoldsp(n,nmy) + c1*rssprk(n,nmy,1) + c2*rssprk(n,nmy,2) &
                                 + c2*rssprk(n,nmy,3) + c1*rssprk(n,nmy,4)
      usp(n,nmz) = uoldsp(n,nmz) + c1*rssprk(n,nmz,1) + c2*rssprk(n,nmz,2) &
                                 + c2*rssprk(n,nmz,3) + c1*rssprk(n,nmz,4)
      usp(n,nee) = uoldsp(n,nee) + c1*rssprk(n,nee,1) + c2*rssprk(n,nee,2) &
                                 + c2*rssprk(n,nee,3) + c1*rssprk(n,nee,4)
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      !
      c1 = one6*dtsp(n)
      c2 = one3*dtsp(n)
      !
      usp(n,nec) = uoldsp(n,nec) + c1*rssprk(n,nec,1) + c2*rssprk(n,nec,2) &
                                 + c2*rssprk(n,nec,3) + c1*rssprk(n,nec,4)
      usp(n,nmx) = uoldsp(n,nmx) + c1*rssprk(n,nmx,1) + c2*rssprk(n,nmx,2) &
                                 + c2*rssprk(n,nmx,3) + c1*rssprk(n,nmx,4)
      usp(n,nmy) = uoldsp(n,nmy) + c1*rssprk(n,nmy,1) + c2*rssprk(n,nmy,2) &
                                 + c2*rssprk(n,nmy,3) + c1*rssprk(n,nmy,4)
      usp(n,nee) = uoldsp(n,nee) + c1*rssprk(n,nee,1) + c2*rssprk(n,nee,2) &
                                 + c2*rssprk(n,nee,3) + c1*rssprk(n,nee,4)
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        c1 = one6*dtsp(n)
        c2 = one3*dtsp(n)
        !
        usp(n,m) = uoldsp(n,m) + c1*rssprk(n,m,1) + c2*rssprk(n,m,2) &
                               + c2*rssprk(n,m,3) + c1*rssprk(n,m,4)
        !
      end do
    end do
  end if
  !
end subroutine update_Classic_RK4_S4
!
!###############################################################################
!
subroutine write_1d_array(io_start,iter,nst,array)
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: io_start,iter,nst
  real(wp), intent(in) :: array(:)
  !
  !.. Local Scalars ..
  integer :: io,n
  !
continue
  !
  !$ACC UPDATE SELF(array)
  io = io_start + 10*(iter-1) + nst
  do n = 1,size(array,dim=1)
    write (io,100) array(n)
  end do
  100 format (es16.8)
  !$ACC UPDATE DEVICE(array)
  !
end subroutine write_1d_array
!
!###############################################################################
!
subroutine write_2d_array(io_start,iter,nst,array)
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: io_start,iter,nst
  real(wp), intent(in) :: array(:,:)
  !
  !.. Local Scalars ..
  integer :: io,n,m
  !
continue
  !
  !$ACC UPDATE SELF(array)
  io = io_start + 10*(iter-1) + nst
  do n = 1,size(array,dim=1)
    write (io,100) (array(n,m), m=1,size(array,dim=2))
  end do
  100 format (5es16.8)
  !$ACC UPDATE DEVICE(array)
  !
end subroutine write_2d_array
!
!###############################################################################
!
subroutine write_3d_array(io_start,iter,nst,array)
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: io_start,iter,nst
  real(wp), intent(in) :: array(:,:,:)
  !
  !.. Local Scalars ..
  integer :: io,n,m,l
  !
continue
  !
  !$ACC UPDATE SELF(array)
  io = io_start + 10*(iter-1) + nst
  do n = 1,size(array,dim=1)
    write (io,100) ((array(n,l,m), l=1,size(array,dim=2)), &
                                   m=1,size(array,dim=3))
  end do
  100 format (15es11.3)
  !$ACC UPDATE DEVICE(array)
  !
end subroutine write_3d_array
!
!###############################################################################
!
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module time_mod
