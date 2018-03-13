module time_mod
  !
!#define PROFILE_ON
  !
  use module_kind_types, global_debug_timer => debug_timer
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq,ntk,ntl
  !
  implicit none
  !
  private
  !
  public :: get_dtsp
  public :: rresi
  public :: update_Classic_RK
  public :: update_CK_RK4
  public :: update_TVD_RK2
  public :: update_TVD_RK3
  public :: update_TVD_RK4
  public :: save_old_solution
  !
 !real(wp), parameter :: o = zero
 !!
 !! Classic N stage Runge-Kutta coefficients
 !!
 !real(wp), parameter, dimension(4,4,4) :: &
 !alphaRK = reshape( &
 !                      ! 1-Stage Coefficients
 !     (/ONE ,o,o,o,o, o  ,o,o, o  ,o,  o   ,o, o  , o  , o  , o  ,   &
 !                      ! 2-Stage Coefficients
 !       HALF,o,o,o,o,ONE ,o,o, o  ,o,  o   ,o, o  , o  , o  , o  ,   &
 !                      ! 3-Stage Coefficients
 !       ONE3,o,o,o,o,TWO3,o,o,ONE4,o,THREE4,o, o  , o  , o  , o  ,   &
 !                      ! 4-Stage Coefficients
 !       HALF,o,o,o,o,HALF,o,o, o  ,o, ONE  ,o,ONE6,ONE3,ONE3,ONE6/), &
 !     (/4,4,4/) )
 !!
 !! TVD 4th-order/5-stage Runge-Kutta coefficients
 !!
 !real(wp), parameter, dimension(1:2,1:5) :: &
 !  betaRK_TVD = reshape( (/ 0.00000000000000_wp,0.39175222700395_wp, &
 !                           0.00000000000000_wp,0.36841059262959_wp, &
 !                           0.00000000000000_wp,0.25189177424738_wp, &
 !                           0.00000000000000_wp,0.54497475021237_wp, &
 !                           0.08460416338212_wp,0.22600748319395_wp  /), &
 !                        (/2,5/) )
 !!
 !real(wp), parameter, dimension(1:3,1:5) :: &
 !  alphaRK_TVD = reshape( &
 !       (/ 1.00000000000000_wp,0.00000000000000_wp,0.0_wp, &
 !          0.44437049406734_wp,0.55562950593266_wp,0.0_wp, &
 !          0.62010185138540_wp,0.37989814861460_wp,0.0_wp, &
 !          0.17807995410773_wp,0.82192004589227_wp,0.0_wp, &
 !          0.00683325884039_wp,0.34833675773694_wp,1.0_wp  /), &
 !       (/3,5/) )
 !!
 !real(wp), parameter, dimension(1:5) :: gammaRK_TVD = (/ 0.00000000000000_wp, &
 !                                                        0.51723167208978_wp, &
 !                                                        0.12759831133288_wp, &
 !                                                        0.00000000000000_wp, &
 !                                                        0.00000000000000_wp /)
 !!
 !character(len=100), save :: iter_out_fmt
 !!
 !real(wp), save, allocatable :: ave(:)
 !!
 !type :: time_ave_var_t
 !  integer, allocatable :: idx(:)
 !end type time_ave_var_t
 !!
 !type(time_ave_var_t), save, allocatable :: time_ave_var(:)
  !
  interface write_array
    module procedure write_1d_array, write_2d_array, write_3d_array
  end interface write_array
  !
contains
!
!###############################################################################
!
subroutine get_dtsp(dtmin)
  !
  !.. Use Statements ..
  use ovar, only : Timestep_Type,d_t,itestcase,Final_Time
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  !
continue
  !
  if (Timestep_Type == Constant_Timestep) then
    !
    call get_current_cfl
    !
    call get_constant_timestep(dtmin)
    !
  else if (abs(Timestep_Type) == Cell_Timestep) then
    !
    call get_current_cfl
    !
    if (Timestep_Type == Global_Cell_Timestep) then
      call get_cell_global_timestep_gpu(dtmin)
     !call get_cell_timestep(dtmin)
    else if (Timestep_Type == Local_Cell_Timestep) then
      call get_cell_local_timestep_gpu(dtmin)
     !call get_cell_timestep(dtmin)
    end if
    !
  else if (abs(Timestep_Type) == Point_Timestep) then
    !
    call get_current_cfl
    !
    call get_point_timestep(dtmin)
    !
  end if
  !
  if (itestcase == Taylor_Green_Vortex) then
    call check_tstar_for_TaylorGreen(dtmin)
  end if
  !
  ! Make sure the timestep size doesnt make the final time exceed Final_Time
  !
  if (Final_time > zero) then
    call check_final_time(dtmin)
  end if
  !
  !
  d_t = dtmin
  !
end subroutine get_dtsp
!
!###############################################################################
!
subroutine get_current_cfl
  !
  !.. Use Statements ..
  use ovar, only : cfl,cfl_beg,cfl_end,cfl_cycles,itcur
  use ovar, only : cfl_cycle_start,cfl_cycle_end
  !
  !.. Local Scalars ..
  real(wp) :: cycle_ratio
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
subroutine get_constant_timestep(dtmin)
  !
  !.. Use Statements ..
  use ovar,    only : constant_dt,time_ref,Final_Time,time,CFL
  use ovar,    only : this_is_final_timestep
  use flowvar, only : dtsp
  use geovar,  only : n_solpts
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  dtmin = CFL * constant_dt
  !
  if (Final_Time > zero) then
    if (time + dtmin > Final_Time/time_ref + eps9) then
      dtmin = Final_Time/time_ref - time
      this_is_final_timestep = true
    end if
  end if
  !
  !$ACC PARALLEL LOOP PRESENT(dtsp)
  do n = 1,n_solpts
    dtsp(n) = dtmin
  end do
  !
end subroutine get_constant_timestep
!
!###############################################################################
!
subroutine get_cell_global_timestep_gpu(dtmin)
  !
  !.. Use Statements ..
  use geovar, only : nfbnd,nface,n_solpts
  use geovar, only : face,cell,n_bnd_flxpts,n_flxpts
  !
  use ovar, only : CFL,governing_equations
  !
  use metrics_mod, only : facenrm,facejac
  !
  use flowvar, only : dtsp
  use flowvar, only : faceusp_array
  !
  use gpu_mod, only : viscosity_cv_gpu
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
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
  character(len=*), parameter :: pname = "get_cell_global_timestep_gpu"
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
    call viscosity_cv_gpu(n1,n2,faceusp_array(:,:,1),viscl)
    call viscosity_cv_gpu(n1,n2,faceusp_array(:,:,2),viscr)
  else
    !$ACC PARALLEL LOOP PRESENT(viscl,viscr)
    do n = n1,n2
      viscl(n) = zero
      viscr(n) = zero
    end do
  end if
  !
  call invs_eigenvalue_gpu(n1,n2,faceusp_array,facenrm,eigni,dsf)
 !call invs_eigenvalue_gpu(n1,n2,faceusp_array(:,:,1),facenrm,eigni(:,1),dsf)
 !call invs_eigenvalue_gpu(n1,n2,faceusp_array(:,:,2),facenrm,eigni(:,2),dsf)
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
  !$ACC PARALLEL LOOP PRESENT(dtsp)
  do n = 1,n_solpts
    dtsp(n) = dtmin
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
end subroutine get_cell_global_timestep_gpu
!
!###############################################################################
!
subroutine get_cell_local_timestep_gpu(dtmin)
  !
  !.. Use Statements ..
  use geovar, only : nfbnd,nface,ncell
  use geovar, only : face,cell,n_flxpts
  !
  use ovar, only : CFL,governing_equations
  !
  use metrics_mod, only : facenrm,facejac
  !
  use flowvar, only : dtsp
  use flowvar, only : faceusp_array
  !
  use gpu_mod, only : viscosity_cv_gpu
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
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
  character(len=*), parameter :: pname = "get_cell_local_timestep_gpu"
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
    call viscosity_cv_gpu(n1,n2,faceusp_array(:,:,1),viscl)
    call viscosity_cv_gpu(n1,n2,faceusp_array(:,:,2),viscr)
  else
    !$ACC PARALLEL LOOP PRESENT(viscl,viscr)
    do n = n1,n2
      viscl(n) = zero
      viscr(n) = zero
    end do
  end if
  !
  call invs_eigenvalue_gpu(n1,n2,faceusp_array,facenrm,eigni,dsf)
 !call invs_eigenvalue_gpu(n1,n2,faceusp_array(:,:,1),facenrm,eigni(:,1),dsf)
 !call invs_eigenvalue_gpu(n1,n2,faceusp_array(:,:,2),facenrm,eigni(:,2),dsf)
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
  !$ACC&              PRESENT(cell,dt_cell,dtsp)
  do nc = 1,ncell
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    dtc = dt_cell(nc)
    dtmin = min(dtmin,dtc)
    do n = n1,n2
      dtsp(n) = dtc
    end do
  end do
  !
 !!$ACC UPDATE SELF(dtmin)
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
    call mpi_allreduce(MPI_IN_PLACE,dtmin,1_int_mpi,mpi_flttyp, &
                       MPI_MIN,MPI_COMM_WORLD,mpierr)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_cell_local_timestep_gpu
!
!###############################################################################
!
pure subroutine invs_eigenvalue_gpu(n1,n2,face_cv,facenrm,eigni,ds)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: eigni(:,:)
  real(wp), intent(inout) :: ds(:)
  !
  !.. Local Scalars ..
  integer :: nr
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  if (nr == 3) then
    call invs_eigenvalue_gpu_3d(n1,n2,face_cv,facenrm,eigni,ds)
  else if (nr == 2) then
    call invs_eigenvalue_gpu_2d(n1,n2,face_cv,facenrm,eigni,ds)
  else
    call invs_eigenvalue_gpu_nd(n1,n2,face_cv,facenrm,eigni,ds)
  end if
  !
end subroutine invs_eigenvalue_gpu
!
!###############################################################################
!
pure subroutine invs_eigenvalue_gpu_3d(n1,n2,face_cv,facenrm,eigni,ds)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: eigni(:,:)
  real(wp), intent(inout) :: ds(:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,rho,rho_inv,rhoe,pres
  real(wp) :: velx,vely,velz,q2,vface,aspd
  real(wp) :: rnx,rny,rnz,unx,uny,unz,dl
  !
continue
  !
  gm1 = gam - one
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,facenrm,eigni,ds)
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
    rho = face_cv(n,nec,1)
    !
    rho_inv = one/rho
    !
    velx = face_cv(n,nmx,1)*rho_inv
    vely = face_cv(n,nmy,1)*rho_inv
    velz = face_cv(n,nmz,1)*rho_inv
    rhoe = face_cv(n,nee,1)
    !
    q2 = velx*velx + vely*vely + velz*velz
    !
    pres = gm1*(rhoe - half*rho*q2)
    !
    aspd = sqrt( max( gam*pres*rho_inv , rndoff ) )
    !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    vface = unx*velx + uny*vely + unz*velz
#else
    vface = sqrt(q2)
#endif
    !
    eigni(n,1) = abs(vface) + aspd
    !
    ! Right side of face
    !
    rho = face_cv(n,nec,2)
    !
    rho_inv = one/rho
    !
    velx = face_cv(n,nmx,2)*rho_inv
    vely = face_cv(n,nmy,2)*rho_inv
    velz = face_cv(n,nmz,2)*rho_inv
    rhoe = face_cv(n,nee,2)
    !
    q2 = velx*velx + vely*vely + velz*velz
    !
    pres = gm1*(rhoe - half*rho*q2)
    !
    aspd = sqrt( max( gam*pres*rho_inv , rndoff ) )
    !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    vface = unx*velx + uny*vely + unz*velz
#else
    vface = sqrt(q2)
#endif
    !
    eigni(n,2) = abs(vface) + aspd
    !
  end do
  !
end subroutine invs_eigenvalue_gpu_3d
!
!###############################################################################
!
pure subroutine invs_eigenvalue_gpu_2d(n1,n2,face_cv,facenrm,eigni,ds)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: eigni(:,:)
  real(wp), intent(inout) :: ds(:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: gm1,rho,rho_inv,rhoe,pres
  real(wp) :: velx,vely,q2,vface,aspd
  real(wp) :: rnx,rny,unx,uny,dl
  !
continue
  !
  gm1 = gam - one
  !
  !$ACC PARALLEL LOOP PRESENT(face_cv,facenrm,eigni,ds)
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
    rho = face_cv(n,nec,1)
    !
    rho_inv = one/rho
    !
    velx = face_cv(n,nmx,1)*rho_inv
    vely = face_cv(n,nmy,1)*rho_inv
    rhoe = face_cv(n,nee,1)
    !
    q2 = velx*velx + vely*vely
    !
    pres = gm1*(rhoe - half*rho*q2)
    !
    aspd = sqrt( max( gam*pres*rho_inv , rndoff ) )
    !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    vface = unx*velx + uny*vely
#else
    vface = sqrt(q2)
#endif
    !
    eigni(n,1) = abs(vface) + aspd
    !
    ! Right side of face
    !
    rho = face_cv(n,nec,2)
    !
    rho_inv = one/rho
    !
    velx = face_cv(n,nmx,2)*rho_inv
    vely = face_cv(n,nmy,2)*rho_inv
    rhoe = face_cv(n,nee,2)
    !
    q2 = velx*velx + vely*vely
    !
    pres = gm1*(rhoe - half*rho*q2)
    !
    aspd = sqrt( max( gam*pres*rho_inv , rndoff ) )
    !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
    vface = unx*velx + uny*vely
#else
    vface = sqrt(q2)
#endif
    !
    eigni(n,2) = abs(vface) + aspd
    !
  end do
  !
end subroutine invs_eigenvalue_gpu_2d
!
!###############################################################################
!
pure subroutine invs_eigenvalue_gpu_nd(n1,n2,face_cv,facenrm,eigni,ds)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: n1,n2
  real(wp),    intent(in) :: face_cv(:,:,:)
  real(wp),    intent(in) :: facenrm(:,:)
  real(wp), intent(inout) :: eigni(:,:)
  real(wp), intent(inout) :: ds(:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: gm1,rho,rho_inv,rhoe,pres
  real(wp) :: velx,vely,velz,q2,vface,aspd
  real(wp) :: rnx,rny,rnz,unx,uny,unz,dl
  !
continue
  !
  nr = size(facenrm,dim=1)
  !
  gm1 = gam - one
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(face_cv,facenrm,eigni,ds)
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
      rho = face_cv(n,nec,1)
      !
      rho_inv = one/rho
      !
      velx = face_cv(n,nmx,1)*rho_inv
      vely = face_cv(n,nmy,1)*rho_inv
      velz = face_cv(n,nmz,1)*rho_inv
      rhoe = face_cv(n,nee,1)
      !
      q2 = velx*velx + vely*vely + velz*velz
      !
      pres = gm1*(rhoe - half*rho*q2)
      !
      aspd = sqrt( max( gam*pres*rho_inv , rndoff ) )
      !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
      vface = unx*velx + uny*vely + unz*velz
#else
      vface = sqrt(q2)
#endif
      !
      eigni(n,1) = abs(vface) + aspd
      !
      ! Right side of face
      !
      rho = face_cv(n,nec,2)
      !
      rho_inv = one/rho
      !
      velx = face_cv(n,nmx,2)*rho_inv
      vely = face_cv(n,nmy,2)*rho_inv
      velz = face_cv(n,nmz,2)*rho_inv
      rhoe = face_cv(n,nee,2)
      !
      q2 = velx*velx + vely*vely + velz*velz
      !
      pres = gm1*(rhoe - half*rho*q2)
      !
      aspd = sqrt( max( gam*pres*rho_inv , rndoff ) )
      !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
      vface = unx*velx + uny*vely + unz*velz
#else
      vface = sqrt(q2)
#endif
      !
      eigni(n,2) = abs(vface) + aspd
      !
    end do
    !
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(face_cv,facenrm,eigni,ds)
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
      rho = face_cv(n,nec,1)
      !
      rho_inv = one/rho
      !
      velx = face_cv(n,nmx,1)*rho_inv
      vely = face_cv(n,nmy,1)*rho_inv
      rhoe = face_cv(n,nee,1)
      !
      q2 = velx*velx + vely*vely
      !
      pres = gm1*(rhoe - half*rho*q2)
      !
      aspd = sqrt( max( gam*pres*rho_inv , rndoff ) )
      !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
      vface = unx*velx + uny*vely
#else
      vface = sqrt(q2)
#endif
      !
      eigni(n,1) = abs(vface) + aspd
      !
      ! Right side of face
      !
      rho = face_cv(n,nec,2)
      !
      rho_inv = one/rho
      !
      velx = face_cv(n,nmx,2)*rho_inv
      vely = face_cv(n,nmy,2)*rho_inv
      rhoe = face_cv(n,nee,2)
      !
      q2 = velx*velx + vely*vely
      !
      pres = gm1*(rhoe - half*rho*q2)
      !
      aspd = sqrt( max( gam*pres*rho_inv , rndoff ) )
      !
#ifdef USE_VN_FOR_INVS_EIGENVALUES
      vface = unx*velx + uny*vely
#else
      vface = sqrt(q2)
#endif
      !
      eigni(n,2) = abs(vface) + aspd
      !
    end do
    !
  else
    !$ACC PARALLEL LOOP PRESENT(face_cv,facenrm,eigni,ds)
    do n = n1,n2
      !
      ds(n) = abs(facenrm(1,n))
      !
      ! Left side of face
      !
      rho = face_cv(n,nec,1)
      velx = face_cv(n,nmx,1)/rho
      rhoe = face_cv(n,nee,1)
      !
      q2 = velx*velx
      !
      pres = gm1*(rhoe - half*rho*velx*velx)
      !
      aspd = sqrt( max( gam*pres/rho , rndoff ) )
      !
      eigni(n,1) = abs(velx) + aspd
      !
      ! Right side of face
      !
      rho = face_cv(n,nec,2)
      velx = face_cv(n,nmx,2)/rho
      rhoe = face_cv(n,nee,2)
      !
      q2 = velx*velx
      !
      pres = gm1*(rhoe - half*rho*q2)
      !
      aspd = sqrt( max( gam*pres/rho , rndoff ) )
      !
      eigni(n,2) = abs(velx) + aspd
      !
    end do
    !
  end if
  !
end subroutine invs_eigenvalue_gpu_nd
!
!###############################################################################
!
subroutine get_cell_timestep(dtmin)
  !
  !.. Use Statements ..
  use order_mod, only : n_order,maxSP,geom_solpts
  !
  use geovar, only : nr,nfbnd,nface
  use geovar, only : face,cell,fp
  !
  use ovar, only : gam,Timestep_Type,time_ref,Final_Time,time,CFL
  use ovar, only : governing_equations
  use ovar, only : this_is_final_timestep
  !
  use metrics_mod, only : geofa,metrics_dt,facejac
  !
  use flowvar, only : dtsp
  use flowvar, only : faceusp_array
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  !
  !.. Local Scalars ..
  integer  :: nf,kf,sf,m
  integer  :: l1,l2,npl
  integer  :: r1,r2,npr
  integer  :: face_geom,face_order,nfp,f0
  integer  :: host_cell,host_geom,host_order
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  real(wp) :: ds,eigni,eignv,pp1sq
  real(wp) :: dsl,vell,al,jacl,dtl,mul,hl
  real(wp) :: dsr,velr,ar,jacr,dtr,mur,hr
  real(wp) :: psql,psqr
  !
  !.. Local Arrays ..
  real(wp) :: vl(1:nq)
  real(wp) :: vr(1:nq)
  !
  integer, parameter :: make_p_plus_one = 0
 !integer, parameter :: make_p_plus_one = 1
  !
  logical(lk), parameter :: use_sqrt = true
 !logical(lk), parameter :: use_sqrt = fals
  !
  integer, parameter :: lsf = 1
  integer, parameter :: rsf = 2
  !
continue
  !
  mul = zero
  mur = zero
  !
  dtmin = huge(zero)
  dtsp(:) = dtmin
  !
  do nf = 1,nfbnd
    !
    host_cell = face(nf)%left%cell ! left cell on face
    !
    host_geom  = cell(host_cell)%geom  ! host cell geometry
    host_order = cell(host_cell)%order ! host cell order
    !
    l1 = cell(host_cell)%beg_sp ! beginning index for sol-pts in left cell
    l2 = cell(host_cell)%end_sp ! ending    index for sol-pts in left cell
    npl = l2-l1+1               ! number of solution points in left cell
    !
    face_geom = face(nf)%geom   ! geometry of the current face
    face_order = face(nf)%order ! order    of the current face
    !
    psql = real( (host_order+make_p_plus_one)**2 , kind=wp )
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    do kf = 1,nfp
      !
      ! Get the primitive variables at this face point
      !
      vl(1:nq) = usp2v_sp( faceusp_array(f0+kf,1:nq,lsf) )
      !
      if (governing_equations == NavierStokes_Eqns) then
        mul = viscosity_pv_sp( vl(1:nq) )
      end if
      !
      ds = norm2( geofa(nf)%v(1:nr,kf) )
      dsl = ds * cl(host_geom) ! face length/area for left cell
      !
      ! Normal velocity at this flux point
      !
     !vell = dot( vl(nmb:nme) , unit_vector(geofa(nf)%v(1:nr,kf)) )
      vell = selfdot( vl(nmb:nme) )
      if (use_sqrt) vell = sqrt(vell)
      !
      ! Speed of sound at this face
      !
      al = max( gam*vl(nee)/vl(nec) , rndoff )
      if (use_sqrt) al = sqrt(al)
      !
      hl = half * dsl / facejac(f0+kf,1)
      !
      eigni = abs(vell) + al
      eignv = psql * hl * mul
      !
      ! Compute the timestep for this flux point
      !
      dtl = min( dtsp(l1) , CFL / (psql * hl * (eigni + eignv)) )
      !
      ! Store the minimum timestep for this cell in the location
      ! of the first solution point for this cell
      !
      dtsp(l1:l2) = dtl
      !
      dtmin = min( dtmin , dtl )
      !
    end do
    !
  end do
  !
  do nf = nfbnd+1,nface
    !
    left_cell  = face(nf)%left%cell  ! left  cell on face
    right_cell = face(nf)%right%cell ! right cell on face
    !
    left_geom  = cell(left_cell)%geom  ! left cell geometry
    left_order = cell(left_cell)%order ! left cell order
    !
    l1 = cell(left_cell)%beg_sp ! beginning index for sol-pts in left  cell
    l2 = cell(left_cell)%end_sp ! ending    index for sol-pts in left cell
    npl = l2-l1+1               ! number of solution points in left  cell
    !
    right_geom  = cell(right_cell)%geom  ! right cell geometry
    right_order = cell(right_cell)%order ! right cell order
    !
    r1 = cell(right_cell)%beg_sp ! beginning index for sol-pts in right cell
    r2 = cell(right_cell)%end_sp ! ending    index for sol-pts in right cell
    npr = r2-r1+1                ! number of solution points in right cell
    !
    face_geom = face(nf)%geom   ! geometry of the current face
    face_order = face(nf)%order ! order    of the current face
    !
    psql = real( (left_order +make_p_plus_one)**2 , kind=wp )
    psqr = real( (right_order+make_p_plus_one)**2 , kind=wp )
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    do kf = 1,nfp
      !
      ! Get the primitive variables at this face point
      !
      vl(1:nq) = usp2v_sp( faceusp_array(f0+kf,1:nq,lsf) )
      vr(1:nq) = usp2v_sp( faceusp_array(f0+kf,1:nq,rsf) )
      !
      if (governing_equations == NavierStokes_Eqns) then
        mul = viscosity_pv_sp( vl(1:nq) )
        mur = viscosity_pv_sp( vr(1:nq) )
      end if
      !
      ds = norm2( geofa(nf)%v(1:nr,kf) )
      dsl = ds * cl(left_geom) ! face length/area for left  cell
      dsr = ds * cl(right_geom) ! face length/area for right cell
      !
      ! Normal velocity at this flux point
      !
     !vell = dot( vl(nmb:nme) , unit_vector(geofa(nf)%v(1:nr,kf)) )
     !velr = dot( vr(nmb:nme) , unit_vector(geofa(nf)%v(1:nr,kf)) )
      vell = selfdot( vl(nmb:nme) )
      if (use_sqrt) vell = sqrt(vell)
      velr = selfdot( vr(nmb:nme) )
      if (use_sqrt) velr = sqrt(velr)
      !
      ! Speed of sound at this face
      !
      al = max( gam*vl(nee)/vl(nec) , rndoff )
      if (use_sqrt) al = sqrt(al)
      ar = max( gam*vr(nee)/vr(nec) , rndoff )
      if (use_sqrt) ar = sqrt(ar)
      !
      hl = half * dsl / facejac(f0+kf,1)
      hr = half * dsr / facejac(f0+kf,2)
      !
      eigni = max( abs(vell) + al , abs(velr) + ar )
      eignv = max( psql * hl * mul , psqr * hr * mur )
      !
      ! Compute the timestep for this flux point
      !
      dtl = min( dtsp(l1) , CFL / (psql * hl * (eigni + eignv) ) )
      dtr = min( dtsp(r1) , CFL / (psqr * hr * (eigni + eignv) ) )
      !
      ! Store the minimum timestep for all solution points in both cells
      !
      dtsp(l1:l2) = dtl
      dtsp(r1:r2) = dtr
      !
      dtmin = min( dtmin , dtl , dtr )
      !
    end do
    !
  end do
  !
  ! Get the minimum timestep across all processors
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,dtmin,1_int_mpi,mpi_flttyp, &
                       MPI_MIN,MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Limit the local timesteps if using a global timestepping
  !
  if (Timestep_Type == Global_Cell_Timestep) then
    !
    ! Set all time step sizes to the global minimum timestep
    !
    dtsp(:) = dtmin
    !
  end if
  !
  ! Make sure the timestep size doesnt make the final time exceed Final_Time
  !
  if (Final_Time > zero) then
    if (time + dtmin > Final_Time/time_ref + eps9) then
      dtsp(:) = Final_Time/time_ref - time
      this_is_final_timestep = true
    end if
  end if
  !
end subroutine get_cell_timestep
!
!###############################################################################
!
subroutine get_point_timestep(dtmin)
  !
  !.. Use Statements ..
  use geovar,  only : nr,ncell,cell
  use ovar,    only : gam,time_ref,Final_Time,time,CFL,Timestep_Type
  use ovar,    only : this_is_final_timestep
  use ovar,    only : governing_equations
  use metrics_mod,  only : metrics_dt
  use flowvar, only : usp,dtsp
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  !
  !.. Local Scalars ..
  integer  :: n,nc,k,n1,n2,np
  integer  :: this_geom,this_order
  real(wp) :: dtpt,dtc,eigni,eignv
  real(wp) :: vmag,aspd,psq,visc
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: vel
  real(wp), dimension(1:nq) :: pv
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_point_timestep"
  !
 !integer, parameter :: make_p_plus_one = 0
  integer, parameter :: make_p_plus_one = 1
  !
continue
  !
  dtmin = huge(zero)
  !
  visc = zero
  !
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp ! beginning index of sol-pts for this cell
    n2 = cell(nc)%end_sp ! ending    index of sol-pts for this cell
    np = n2-n1+1         ! number of sol-pts in this cell
    !
    this_geom = cell(nc)%geom   ! geometry type of this cell
    this_order = cell(nc)%order ! solution order of this cell
    !
    psq = real( (this_order+make_p_plus_one)**2 , kind=wp )
    !
    dtc = CFL * two / psq
    !
    do k = 1,np
      !
      n = k + n1-1
      !
      pv(1:nq) = usp2v_sp( usp(n,1:nq) )
      !
      if (governing_equations == NavierStokes_Eqns) then
        visc = viscosity_cv_sp( usp(n,1:nq) )
      end if
      !
      vel(1:nr) = matmul( metrics_dt(nc)%met(k,1:nr,1:nr) , pv(nmb:nme) )
      !
      vmag = norm2(vel)
      !
      aspd = sqrt( max( gam*pv(nee)/pv(nec) , rndoff ) )
      !
      eigni = vmag + aspd
      !
      eignv = psq * visc
      !
      dtpt = dtc / (eigni + eignv)
      !
      dtmin = min( dtpt , dtmin )
      !
      dtsp(n) = dtpt
      !
    end do
    !
  end do
  !
  ! Get the minimum timestep across all processors
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,dtmin,1_int_mpi,mpi_flttyp, &
                       MPI_MIN,MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Limit the local timesteps if using a global timestepping
  !
  if (Timestep_Type == Global_Point_Timestep) then
    !
    ! Set all time step sizes to the minimum timestep
    !
    dtsp(:) = dtmin
    !
    ! Make sure the timestep size doesnt make the final time exceed Final_Time
    !
    if (Final_Time > zero) then
      if (time + dtmin > Final_Time/time_ref + eps9) then
        dtsp(:) = Final_Time/time_ref - time
        this_is_final_timestep = true
      end if
    end if
    !
  end if
  !
end subroutine get_point_timestep
!
!###############################################################################
!
subroutine old_get_point_timestep(dtmin)
  !
  !.. Use Statements ..
  use order_mod,  only : n_order
  use geovar,  only : nr,ncell,cell
  use ovar,    only : gam,time_ref,Final_Time,time,CFL,Timestep_Type
  use ovar,    only : this_is_final_timestep
  use metrics_mod,  only : metrics_dt
  use flowvar, only : usp,dtsp
  !
  !.. Formal Arguments ..
  real(wp), intent(out) :: dtmin
  !
  !.. Local Scalars ..
  integer  :: k,n,nc,n1,n2
  real(wp) :: dtpt,gm1,dtc,eign
  real(wp) :: dxdr,dydr,dxds,dyds
  real(wp) :: rho,vx,vy,pres,aspd
  real(wp) :: alfa,beta,gama
  !
  real(wp) :: vel(1:nr)
  real(wp) :: met(1:nr,1:nr)
  !
continue
  !
  gm1 = gam - one
  !
  dtmin = huge(zero)
  !
  dtc = CFL * two / real( (n_order+1)**2 , kind=wp )
  !
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    !
    do n = n1,n2
      !
      k = n-n1+1
      !
      met(1:nr,1:nr) = invert_matrix(metrics_dt(nc)%met(k,1:nr,1:nr))
      dxdr = met(1,1)
      dydr = met(2,1)
      dxds = met(1,2)
      dyds = met(2,1)
      !
      rho  = usp(n,nec)
      vel  = usp(n,nmb:nme)/rho
      pres = gm1 * (usp(n,nee) - half*rho*selfdot(vel))
      !
      aspd = sqrt( max( gam*pres/rho , rndoff ) )
      !
      alfa = abs( vel(1)*dydr - vel(2)*dxdr )
      beta = abs( vel(1)*dyds - vel(2)*dxds )
      !
      gama = aspd*sqrt(dxdr**2 + dydr**2 + dxds**2 + dyds**2)
      !
      eign = alfa + beta + gama
      !
      dtpt = dtc * metrics_dt(nc)%jac(k) / eign
      !
      dtmin = min( dtpt , dtmin )
      !
      dtsp(n) = dtpt
      !
    end do
    !
  end do
  !
  ! Get the minimum timestep across all processors
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,dtmin,1_int_mpi,mpi_flttyp, &
                       MPI_MIN,MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Limit the local timesteps if using a global timestepping
  !
  if (Timestep_Type == Global_Point_Timestep) then
    !
    ! Set all time step sizes to the minimum timestep
    !
    dtsp(:) = dtmin
    !
    ! Make sure the timestep size doesnt make the final time exceed Final_Time
    !
    if (Final_Time > zero) then
      if (time + dtmin > Final_Time/time_ref + eps9) then
        dtsp(:) = Final_Time/time_ref - time
        this_is_final_timestep = true
      end if
    end if
    !
  end if
  !
end subroutine old_get_point_timestep
!
!###############################################################################
!
subroutine check_tstar_for_TaylorGreen(dtmin)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : dtsp
  use ovar,    only : time,time_ref
  use ovar,    only : write_TaylorGreen_8s_solution
  use ovar,    only : write_TaylorGreen_full_solution
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: dtmin
  !
  !.. Local Scalars ..
  integer     :: n
  real(wp)    :: tstar,new_tstar
  logical(lk) :: reset_dtsp
  !
continue
  !
  reset_dtsp = fals
  !
  tstar = time*time_ref
  new_tstar = tstar + dtmin*time_ref
  !
  if (tstar < eight .and. new_tstar >= eight) then
    !
    dtmin = (eight - tstar) / time_ref
    write_TaylorGreen_8s_solution = true
    reset_dtsp = true
    !
  else if (tstar < three .and. new_tstar >= three) then
    !
    dtmin = (three - tstar) / time_ref
    write_TaylorGreen_full_solution = true
    reset_dtsp = true
    !
  else if (tstar < five .and. new_tstar >= five) then
    !
    dtmin = (five - tstar) / time_ref
    write_TaylorGreen_full_solution = true
    reset_dtsp = true
    !
  else if (tstar < seven .and. new_tstar >= seven) then
    !
    dtmin = (seven - tstar) / time_ref
    write_TaylorGreen_full_solution = true
    reset_dtsp = true
    !
  else if (tstar < nine .and. new_tstar >= nine) then
    !
    dtmin = (nine - tstar) / time_ref
    write_TaylorGreen_full_solution = true
    reset_dtsp = true
    !
  else if (tstar < 11.0_wp .and. new_tstar >= 11.0_wp) then
    !
    dtmin = (11.0_wp - tstar) / time_ref
    write_TaylorGreen_full_solution = true
    reset_dtsp = true
    !
  else if (tstar < 15.0_wp .and. new_tstar >= 15.0_wp) then
    !
    dtmin = (15.0_wp - tstar) / time_ref
    write_TaylorGreen_full_solution = true
    reset_dtsp = true
    !
  end if
  !
 !!$ACC UPDATE DEVICE(dtmin)
  !
  if (reset_dtsp) then
    !$ACC PARALLEL LOOP PRESENT(dtsp)
    do n = 1,n_solpts
      dtsp(n) = dtmin
    end do
  end if
  !
end subroutine check_tstar_for_TaylorGreen
!
!###############################################################################
!
subroutine check_final_time(dtmin)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : dtsp
  use ovar,    only : time,Final_Time,time_ref
  use ovar,    only : this_is_final_timestep
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: dtmin
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: new_time,chk_time
  !
continue
  !
  new_time = time + dtmin
  chk_time = Final_Time/time_ref + eps9
  !
  if (new_time > chk_time) then
    dtmin = Final_Time/time_ref - time
    this_is_final_timestep = true
  end if
  !
 !!$ACC UPDATE DEVICE(dtmin)
  !
  if (this_is_final_timestep) then
    !$ACC PARALLEL LOOP PRESENT(dtsp)
    do n = 1,n_solpts
      dtsp(n) = dtmin
    end do
  end if
  !
end subroutine check_final_time
!
!###############################################################################
!
subroutine save_old_solution()
  !
  !.. Use Statements ..
  use flowvar, only : usp,uoldsp
  !
  !.. Local Scalars ..
  integer :: n,np,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "save_old_solution"
  !
continue
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
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(uoldsp,usp)
    do n = 1,np
      uoldsp(n,nec) = usp(n,nec)
      uoldsp(n,nmx) = usp(n,nmx)
      uoldsp(n,nmy) = usp(n,nmy)
      uoldsp(n,nmz) = usp(n,nmz)
      uoldsp(n,nee) = usp(n,nee)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(uoldsp,usp)
    do n = 1,np
      uoldsp(n,nec) = usp(n,nec)
      uoldsp(n,nmx) = usp(n,nmx)
      uoldsp(n,nmy) = usp(n,nmy)
      uoldsp(n,nee) = usp(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(uoldsp,usp)
    do n = 1,np
      uoldsp(n,nec) = usp(n,nec)
      uoldsp(n,nmx) = usp(n,nmx)
      uoldsp(n,nee) = usp(n,nee)
    end do
  end if
  !
  call debug_timer(leaving_procedure,pname)
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
  use geovar,      only : n_solpts,cell
  use ovar,        only : itestcase
  use metrics_mod, only : jac
  use flowvar,     only : residual
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
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "rresi"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(jac,residual)
    do n = 1,n_solpts
      jac_inv = one/jac(n)
      residual(n,nec) = -residual(n,nec)*jac_inv
      residual(n,nmx) = -residual(n,nmx)*jac_inv
      residual(n,nmy) = -residual(n,nmy)*jac_inv
      residual(n,nmz) = -residual(n,nmz)*jac_inv
      residual(n,nee) = -residual(n,nee)*jac_inv
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(jac,residual)
    do n = 1,n_solpts
      jac_inv = one/jac(n)
      residual(n,nec) = -residual(n,nec)*jac_inv
      residual(n,nmx) = -residual(n,nmx)*jac_inv
      residual(n,nmy) = -residual(n,nmy)*jac_inv
      residual(n,nee) = -residual(n,nee)*jac_inv
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(jac,residual)
    do n = 1,n_solpts
      jac_inv = one/jac(n)
      residual(n,nec) = -residual(n,nec)*jac_inv
      residual(n,nmx) = -residual(n,nmx)*jac_inv
      residual(n,nee) = -residual(n,nee)*jac_inv
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
    if (nr == 3) then
      !$ACC PARALLEL LOOP PRESENT(residual)
      do n = 1,n_solpts
        residual(n,nmx) = residual(n,nec)
        residual(n,nmy) = residual(n,nec)
        residual(n,nmz) = residual(n,nec)
        residual(n,nee) = residual(n,nec)
      end do
    else if (nr == 2) then
      !$ACC PARALLEL LOOP PRESENT(residual)
      do n = 1,n_solpts
        residual(n,nmx) = residual(n,nec)
        residual(n,nmy) = residual(n,nec)
        residual(n,nee) = residual(n,nec)
      end do
    else
      !$ACC PARALLEL LOOP PRESENT(residual)
      do n = 1,n_solpts
        residual(n,nmx) = residual(n,nec)
        residual(n,nee) = residual(n,nec)
      end do
    end if
  end if
  !
  if (p_order < n_order) then
    call project_to_porder(residual)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine rresi
!
!###############################################################################
!
subroutine stabilize_usp(nst,Apply_Filter,Apply_Limiter)
  !
  !.. Use Statements ..
  use ovar,            only : Limiter_Option,Filter_Option
  use filter_mod,      only : filter_usp
  use module_limiters, only : resolution_indicator
  use module_limiters, only : limit_by_projecting
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "stabilize_usp"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Apply the filter
  !
  if (Filter_Option > 0) then
    if (Apply_Filter) then
#ifdef PROFILE_ON
      if (mypnum == glb_root) write (iout,1) nst
      1 format (8x,"Applying Filter : nst = ",i0)
#endif
      call filter_usp
    end if
  end if
  !
  ! Apply the limiter
  !
  if (Limiter_Option /= 0) then
    if (Apply_Limiter) then
#ifdef PROFILE_ON
      if (mypnum == glb_root) write (iout,2) nst
      2 format (8x,"Applying Limiter : nst = ",i0)
#endif
      if (Limiter_Option == -1) then
        call resolution_indicator(apply_projection=true)
      else if (Limiter_Option >= 1) then
        call limit_by_projecting
      end if
    end if
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine stabilize_usp
!
!###############################################################################
!
subroutine update_CK_RK4(nst,Apply_Filter,Apply_Limiter)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : rk_sum,residual,usp,dtsp
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
  else if (nst == 2) then
    call update_CK_RK4_S3(n_solpts,dtsp,residual,usp,rk_sum)
  else if (nst == 2) then
    call update_CK_RK4_S4(n_solpts,dtsp,residual,usp,rk_sum)
  else if (nst == 2) then
    call update_CK_RK4_S5(n_solpts,dtsp,residual,usp,rk_sum)
  end if
  !
  ! Apply filtering or limiting to stabilize the solution
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter)
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
  integer  :: n,nr
  real(wp) :: sum_ec,sum_mx,sum_my,sum_mz,sum_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK4_S1"
  !
  real(wp), parameter :: b = 1432997174477.0_wp / 9575080441755.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
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
  else if (nr == 2) then
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
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = dtsp(n)*residual(n,nec)
      sum_mx = dtsp(n)*residual(n,nmx)
      sum_ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nee) = sum_ee
      !
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
  integer  :: n,nr
  real(wp) :: sum_ec,sum_mx,sum_my,sum_mz,sum_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK4_S2"
  !
  real(wp), parameter :: a = -567301805773.0_wp / 1357537059087.0_wp
  real(wp), parameter :: b = 5161836677717.0_wp / 13612068292357.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
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
  else if (nr == 2) then
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
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nee) = sum_ee
      !
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
  integer  :: n,nr
  real(wp) :: sum_ec,sum_mx,sum_my,sum_mz,sum_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK4_S3"
  !
  real(wp), parameter :: a = -2404267990393.0_wp / 2016746695238.0_wp
  real(wp), parameter :: b = 1720146321549.0_wp / 2090206949498.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
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
  else if (nr == 2) then
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
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nee) = sum_ee
      !
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
  integer  :: n,nr
  real(wp) :: sum_ec,sum_mx,sum_my,sum_mz,sum_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK4_S4"
  !
  real(wp), parameter :: a = -3550918686646.0_wp / 2091501179385.0_wp
  real(wp), parameter :: b = 3134564353537.0_wp / 4481467310338.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
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
  else if (nr == 2) then
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
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nee) = sum_ee
      !
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
  integer  :: n,nr
  real(wp) :: sum_ec,sum_mx,sum_my,sum_mz,sum_ee
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK4_S5"
  !
  real(wp), parameter :: a = -1275806237668.0_wp / 842570457699.0_wp
  real(wp), parameter :: b = 2277821191437.0_wp / 14882151754819.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
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
  else if (nr == 2) then
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
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nee) = sum_ee
      !
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
  ! Apply filtering or limiting to stabilize the solution
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter)
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
  integer :: n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK2_S1"
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,usp)
    do n = 1,np
      usp(n,nec) = usp(n,nec) + dtsp(n)*residual(n,nec)
      usp(n,nmx) = usp(n,nmx) + dtsp(n)*residual(n,nmx)
      usp(n,nmy) = usp(n,nmy) + dtsp(n)*residual(n,nmy)
      usp(n,nmz) = usp(n,nmz) + dtsp(n)*residual(n,nmz)
      usp(n,nee) = usp(n,nee) + dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,usp)
    do n = 1,np
      usp(n,nec) = usp(n,nec) + dtsp(n)*residual(n,nec)
      usp(n,nmx) = usp(n,nmx) + dtsp(n)*residual(n,nmx)
      usp(n,nmy) = usp(n,nmy) + dtsp(n)*residual(n,nmy)
      usp(n,nee) = usp(n,nee) + dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,usp)
    do n = 1,np
      usp(n,nec) = usp(n,nec) + dtsp(n)*residual(n,nec)
      usp(n,nmx) = usp(n,nmx) + dtsp(n)*residual(n,nmx)
      usp(n,nee) = usp(n,nee) + dtsp(n)*residual(n,nee)
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
  integer :: n,nr
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
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a * ( uoldsp(n,nec) + usp(n,nec) + dtsp(n)*residual(n,nec) )
      usp(n,nmx) = a * ( uoldsp(n,nmx) + usp(n,nmx) + dtsp(n)*residual(n,nmx) )
      usp(n,nmy) = a * ( uoldsp(n,nmy) + usp(n,nmy) + dtsp(n)*residual(n,nmy) )
      usp(n,nmz) = a * ( uoldsp(n,nmz) + usp(n,nmz) + dtsp(n)*residual(n,nmz) )
      usp(n,nee) = a * ( uoldsp(n,nee) + usp(n,nee) + dtsp(n)*residual(n,nee) )
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a * ( uoldsp(n,nec) + usp(n,nec) + dtsp(n)*residual(n,nec) )
      usp(n,nmx) = a * ( uoldsp(n,nmx) + usp(n,nmx) + dtsp(n)*residual(n,nmx) )
      usp(n,nmy) = a * ( uoldsp(n,nmy) + usp(n,nmy) + dtsp(n)*residual(n,nmy) )
      usp(n,nee) = a * ( uoldsp(n,nee) + usp(n,nee) + dtsp(n)*residual(n,nee) )
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a * ( uoldsp(n,nec) + usp(n,nec) + dtsp(n)*residual(n,nec) )
      usp(n,nmx) = a * ( uoldsp(n,nmx) + usp(n,nmx) + dtsp(n)*residual(n,nmx) )
      usp(n,nee) = a * ( uoldsp(n,nee) + usp(n,nee) + dtsp(n)*residual(n,nee) )
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
  ! Apply filtering or limiting to stabilize the solution
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter)
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
  integer :: n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK3_S1"
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,usp)
    do n = 1,np
      usp(n,nec) = usp(n,nec) + dtsp(n)*residual(n,nec)
      usp(n,nmx) = usp(n,nmx) + dtsp(n)*residual(n,nmx)
      usp(n,nmy) = usp(n,nmy) + dtsp(n)*residual(n,nmy)
      usp(n,nmz) = usp(n,nmz) + dtsp(n)*residual(n,nmz)
      usp(n,nee) = usp(n,nee) + dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,usp)
    do n = 1,np
      usp(n,nec) = usp(n,nec) + dtsp(n)*residual(n,nec)
      usp(n,nmx) = usp(n,nmx) + dtsp(n)*residual(n,nmx)
      usp(n,nmy) = usp(n,nmy) + dtsp(n)*residual(n,nmy)
      usp(n,nee) = usp(n,nee) + dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,usp)
    do n = 1,np
      usp(n,nec) = usp(n,nec) + dtsp(n)*residual(n,nec)
      usp(n,nmx) = usp(n,nmx) + dtsp(n)*residual(n,nmx)
      usp(n,nee) = usp(n,nee) + dtsp(n)*residual(n,nee)
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
  integer :: n,nr
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
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + b*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + b*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + b*dtsp(n)*residual(n,nmy)
      usp(n,nmz) = a*uoldsp(n,nmz) + b*usp(n,nmz) + b*dtsp(n)*residual(n,nmz)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + b*dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + b*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + b*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + b*dtsp(n)*residual(n,nmy)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + b*dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + b*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + b*dtsp(n)*residual(n,nmx)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + b*dtsp(n)*residual(n,nee)
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
  integer :: n,nr
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
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + b*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + b*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + b*dtsp(n)*residual(n,nmy)
      usp(n,nmz) = a*uoldsp(n,nmz) + b*usp(n,nmz) + b*dtsp(n)*residual(n,nmz)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + b*dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + b*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + b*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + b*dtsp(n)*residual(n,nmy)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + b*dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + b*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + b*dtsp(n)*residual(n,nmx)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + b*dtsp(n)*residual(n,nee)
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
  use flowvar, only : rssprk,residual,usp,uoldsp,rk_sum,dtsp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Scalars ..
  integer :: n,m
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
  ! Apply filtering or limiting to stabilize the solution
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter)
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
  integer :: n,nr
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
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = uoldsp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = uoldsp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nmz) = uoldsp(n,nmz) + c*dtsp(n)*residual(n,nmz)
      usp(n,nee) = uoldsp(n,nee) + c*dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = uoldsp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = uoldsp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nee) = uoldsp(n,nee) + c*dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = uoldsp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nee) = uoldsp(n,nee) + c*dtsp(n)*residual(n,nee)
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
  integer :: n,nr
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
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nmz) = a*uoldsp(n,nmz) + b*usp(n,nmz) + c*dtsp(n)*residual(n,nmz)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*dtsp(n)*residual(n,nee)
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
  integer  :: n,nr
  real(wp) :: stage2_ec,stage2_mx,stage2_my,stage2_mz,stage2_ee
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
  if (nr == 3) then
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
  else if (nr == 2) then
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
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual,rk_sum)
    do n = 1,np
      !
      stage2_ec = usp(n,nec)
      stage2_mx = usp(n,nmx)
      stage2_ee = usp(n,nee)
      !
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*dtsp(n)*residual(n,nee)
      !
      rk_sum(n,nec) = d*stage2_ec + e*usp(n,nec)
      rk_sum(n,nmx) = d*stage2_mx + e*usp(n,nmx)
      rk_sum(n,nee) = d*stage2_ee + e*usp(n,nee)
      !
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
  integer :: n,nr
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
  if (nr == 3) then
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
  else if (nr == 2) then
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
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual,rk_sum)
    do n = 1,np
      !
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*dtsp(n)*residual(n,nee)
      !
      rk_sum(n,nec) = rk_sum(n,nec) + d*dtsp(n)*residual(n,nec)
      rk_sum(n,nmx) = rk_sum(n,nmx) + d*dtsp(n)*residual(n,nmx)
      rk_sum(n,nee) = rk_sum(n,nee) + d*dtsp(n)*residual(n,nee)
      !
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
  integer  :: n,nr
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
  if (nr == 3) then
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
  else if (nr == 2) then
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
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,res,rk_sum)
    do n = 1,np
      !
      c = d*dtsp(n)
      !
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*res(n,nec) + rk_sum(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*res(n,nmx) + rk_sum(n,nmx)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*res(n,nee) + rk_sum(n,nee)
      !
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
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Scalars ..
  integer :: n,nr
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
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(residual,rssprk)
    do n = 1,n_solpts
      rssprk(n,nec,nst) = residual(n,nec)
      rssprk(n,nmx,nst) = residual(n,nmx)
      rssprk(n,nmy,nst) = residual(n,nmy)
      rssprk(n,nmz,nst) = residual(n,nmz)
      rssprk(n,nee,nst) = residual(n,nee)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(residual,rssprk)
    do n = 1,n_solpts
      rssprk(n,nec,nst) = residual(n,nec)
      rssprk(n,nmx,nst) = residual(n,nmx)
      rssprk(n,nmy,nst) = residual(n,nmy)
      rssprk(n,nee,nst) = residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(residual,rssprk)
    do n = 1,n_solpts
      rssprk(n,nec,nst) = residual(n,nec)
      rssprk(n,nmx,nst) = residual(n,nmx)
      rssprk(n,nee,nst) = residual(n,nee)
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
  ! Apply filtering or limiting to stabilize the solution
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter)
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
  integer :: n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,1)
      usp(n,nmz) = uoldsp(n,nmz) + dtsp(n)*rssprk(n,nmz,1)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,1)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,1)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,1)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,1)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,1)
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
  integer  :: n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nmz) = uoldsp(n,nmz) + half*dtsp(n)*rssprk(n,nmz,1)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,1)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,1)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,1)
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
  integer :: n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,2)
      usp(n,nmz) = uoldsp(n,nmz) + dtsp(n)*rssprk(n,nmz,2)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,2)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,2)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,2)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,2)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,2)
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
  integer  :: n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + one3*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + one3*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + one3*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nmz) = uoldsp(n,nmz) + one3*dtsp(n)*rssprk(n,nmz,1)
      usp(n,nee) = uoldsp(n,nee) + one3*dtsp(n)*rssprk(n,nee,1)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + one3*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + one3*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + one3*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nee) = uoldsp(n,nee) + one3*dtsp(n)*rssprk(n,nee,1)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + one3*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + one3*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nee) = uoldsp(n,nee) + one3*dtsp(n)*rssprk(n,nee,1)
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
  integer  :: n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + two3*dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + two3*dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + two3*dtsp(n)*rssprk(n,nmy,2)
      usp(n,nmz) = uoldsp(n,nmz) + two3*dtsp(n)*rssprk(n,nmz,2)
      usp(n,nee) = uoldsp(n,nee) + two3*dtsp(n)*rssprk(n,nee,2)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + two3*dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + two3*dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + two3*dtsp(n)*rssprk(n,nmy,2)
      usp(n,nee) = uoldsp(n,nee) + two3*dtsp(n)*rssprk(n,nee,2)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + two3*dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + two3*dtsp(n)*rssprk(n,nmx,2)
      usp(n,nee) = uoldsp(n,nee) + two3*dtsp(n)*rssprk(n,nee,2)
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
  integer  :: n,nr
  real(wp) :: c1,c2
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
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
  else if (nr == 2) then
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
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      !
      c1 = one4*dtsp(n)
      c2 = three4*dtsp(n)
      !
      usp(n,nec) = uoldsp(n,nec) + c1*rssprk(n,nec,1) + c2*rssprk(n,nec,3)
      usp(n,nmx) = uoldsp(n,nmx) + c1*rssprk(n,nmx,1) + c2*rssprk(n,nmx,3)
      usp(n,nee) = uoldsp(n,nee) + c1*rssprk(n,nee,1) + c2*rssprk(n,nee,3)
      !
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
  integer  :: n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nmz) = uoldsp(n,nmz) + half*dtsp(n)*rssprk(n,nmz,1)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,1)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,1)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,1)
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
  integer  :: n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,2)
      usp(n,nmz) = uoldsp(n,nmz) + half*dtsp(n)*rssprk(n,nmz,2)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,2)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,2)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,2)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,2)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,2)
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
  integer :: n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,3)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,3)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,3)
      usp(n,nmz) = uoldsp(n,nmz) + dtsp(n)*rssprk(n,nmz,3)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,3)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,3)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,3)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,3)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,3)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,3)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,3)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,3)
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
  integer  :: n,nr
  real(wp) :: c1,c2
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3) then
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
  else if (nr == 2) then
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
      usp(n,nee) = uoldsp(n,nee) + c1*rssprk(n,nee,1) + c2*rssprk(n,nee,2) &
                                 + c2*rssprk(n,nee,3) + c1*rssprk(n,nee,4)
      !
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
include "Functions/usp2v.f90"
include "Functions/viscosity.f90"
!include "Functions/temperature.f90"
!
!###############################################################################
!
end module time_mod
