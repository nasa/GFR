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
  public :: rupdate_Classic_RK
  public :: rupdate_CK4_RK
  public :: rupdate_TVD2_RK
  public :: rupdate_TVD3_RK
  public :: rupdate_TVD4_RK
 !public :: update_time_ave
 !public :: initialize_time_ave
  public :: save_old_solution
  !
  real(wp), parameter :: o = zero
  !
  ! Classic N stage Runge-Kutta coefficients
  !
  real(wp), parameter, dimension(4,4,4) :: &
  alphaRK = reshape( &
                        ! 1-Stage Coefficients
       (/ONE ,o,o,o,o, o  ,o,o, o  ,o,  o   ,o, o  , o  , o  , o  ,   &
                        ! 2-Stage Coefficients
         HALF,o,o,o,o,ONE ,o,o, o  ,o,  o   ,o, o  , o  , o  , o  ,   &
                        ! 3-Stage Coefficients
         ONE3,o,o,o,o,TWO3,o,o,ONE4,o,THREE4,o, o  , o  , o  , o  ,   &
                        ! 4-Stage Coefficients
         HALF,o,o,o,o,HALF,o,o, o  ,o, ONE  ,o,ONE6,ONE3,ONE3,ONE6/), &
       (/4,4,4/) )
  !
  ! TVD 4th-order/5-stage Runge-Kutta coefficients
  !
  real(wp), parameter, dimension(1:2,1:5) :: &
    betaRK_TVD = reshape( (/ 0.00000000000000_wp,0.39175222700395_wp, &
                             0.00000000000000_wp,0.36841059262959_wp, &
                             0.00000000000000_wp,0.25189177424738_wp, &
                             0.00000000000000_wp,0.54497475021237_wp, &
                             0.08460416338212_wp,0.22600748319395_wp  /), &
                          (/2,5/) )
  !
  real(wp), parameter, dimension(1:3,1:5) :: &
    alphaRK_TVD = reshape( &
         (/ 1.00000000000000_wp,0.00000000000000_wp,0.0_wp, &
            0.44437049406734_wp,0.55562950593266_wp,0.0_wp, &
            0.62010185138540_wp,0.37989814861460_wp,0.0_wp, &
            0.17807995410773_wp,0.82192004589227_wp,0.0_wp, &
            0.00683325884039_wp,0.34833675773694_wp,1.0_wp  /), &
         (/3,5/) )
  !
  real(wp), parameter, dimension(1:5) :: gammaRK_TVD = (/ 0.00000000000000_wp, &
                                                          0.51723167208978_wp, &
                                                          0.12759831133288_wp, &
                                                          0.00000000000000_wp, &
                                                          0.00000000000000_wp /)
  !
  character(len=100), save :: iter_out_fmt
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
 !!$ACC UPDATE SELF(dtmin)
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
  ! Get the minimum timestep across all processors
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,dtmin,1_int_mpi,mpi_flttyp, &
                       MPI_MIN,MPI_COMM_WORLD,mpierr)
  end if
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
#ifdef USE_FACE_ARRAYS
  use flowvar, only : faceusp_array
#else
  use flowvar, only : faceusp
#endif
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
#ifdef USE_FACE_ARRAYS
    nfp = face(nf)%end_fp - f0
#else
    nfp = geom_solpts(face_geom,face_order)
#endif
    !
    do kf = 1,nfp
      !
      ! Get the primitive variables at this face point
      !
#ifdef USE_FACE_ARRAYS
      vl(1:nq) = usp2v_sp( faceusp_array(f0+kf,1:nq,lsf) )
#else
      vl(1:nq) = usp2v_sp( faceusp(nf)%v(1:nq,kf,lsf) )
#endif
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
#ifdef USE_FACE_ARRAYS
    nfp = face(nf)%end_fp - f0
#else
    nfp = geom_solpts(face_geom,face_order)
#endif
    !
    do kf = 1,nfp
      !
      ! Get the primitive variables at this face point
      !
#ifdef USE_FACE_ARRAYS
      vl(1:nq) = usp2v_sp( faceusp_array(f0+kf,1:nq,lsf) )
      vr(1:nq) = usp2v_sp( faceusp_array(f0+kf,1:nq,rsf) )
#else
      vl(1:nq) = usp2v_sp( faceusp(nf)%v(1:nq,kf,lsf) )
      vr(1:nq) = usp2v_sp( faceusp(nf)%v(1:nq,kf,rsf) )
#endif
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
continue
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
    call project_to_porder(cell,residual)
  end if
  !
end subroutine rresi
!
!###############################################################################
!
subroutine rupdate_Classic_RK(nst,Apply_Filter,Apply_Limiter)
  !
  ! First, store d(utd)/dt for runge kutta time
  ! stepping then update usp
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : num_rk_stages
  use ovar,    only : Limiter_Option
  use ovar,    only : Filter_Option
  use flowvar, only : rssprk,residual,usp,uoldsp,dtsp
  use filter_mod, only : filter_usp
  use module_limiters, only : resolution_indicator
  use module_limiters, only : limit_by_projecting
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Scalars ..
  integer   :: n,m,nt
  real(wp) :: temp
  !
continue
  !
  ! Store the residual for this Runge-Kutta step in rssprk
  !
  do m = 1,nq
    do n = 1,n_solpts
      rssprk(n,m,nst) = residual(n,m)
    end do
  end do
  !
  ! Update usp
  !
  do m = 1,nq
    do n = 1,n_solpts
      temp = zero
      do nt = 1,nst
        temp = temp + alphark(nt,nst,num_rk_stages)*rssprk(n,m,nt)
      end do
      usp(n,m) = uoldsp(n,m) + dtsp(n)*temp
    end do
  end do
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
end subroutine rupdate_Classic_RK
!
!###############################################################################
!
subroutine rupdate_CK4_RK(nst,Apply_Filter,Apply_Limiter)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : Limiter_Option,Filter_Option
  use flowvar, only : rssprk,residual,usp,dtsp
  use filter_mod, only : filter_usp
  use module_limiters, only : resolution_indicator
  use module_limiters, only : limit_by_projecting
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Scalars ..
  integer :: n,m
  !
  ! Carpenter-Kennedy 4th-order/5-stage Runge-Kutta coefficients
  !
  real(wp), parameter, dimension(1:5) :: &
             alfa = (/                                   0.0_wp, &
                        -567301805773.0_wp/ 1357537059087.0_wp, &
                       -2404267990393.0_wp/ 2016746695238.0_wp, &
                       -3550918686646.0_wp/ 2091501179385.0_wp, &
                       -1275806237668.0_wp/  842570457699.0_wp  /)

  real(wp), parameter, dimension(1:5) :: &
             beta = (/  1432997174477.0_wp/ 9575080441755.0_wp, &
                        5161836677717.0_wp/13612068292357.0_wp, &
                        1720146321549.0_wp/ 2090206949498.0_wp, &
                        3134564353537.0_wp/ 4481467310338.0_wp, &
                        2277821191437.0_wp/14882151754819.0_wp  /)
  !
continue
  !
  do m = 1,nq
    do n = 1,n_solpts
      rssprk(n,m,1) = alfa(nst)*rssprk(n,m,1) + dtsp(n)*residual(n,m)
      usp(n,m) = usp(n,m) + beta(nst)*rssprk(n,m,1)
    end do
  end do
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
end subroutine rupdate_CK4_RK
!
!###############################################################################
!
subroutine rupdate_TVD2_RK(nst,Apply_Filter,Apply_Limiter)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : Limiter_Option,Filter_Option
  use flowvar, only : residual,usp,uoldsp,dtsp
  use filter_mod, only : filter_usp
  use module_limiters, only : resolution_indicator
  use module_limiters, only : limit_by_projecting
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Scalars ..
  integer   :: n,m
  !
  ! TVD 2nd-order/2-stage Runge-Kutta Coefficients
  !
  real(wp), parameter, dimension(1:2) :: alfa = (/ zero,one  /)
  real(wp), parameter, dimension(1:2) :: beta = (/ one ,half /)
  !
continue
  !
  do m = 1,nq
    do n = 1,n_solpts
      usp(n,m) = beta(nst) * (alfa(nst)*uoldsp(n,m) + &
                              usp(n,m) + &
                              dtsp(n)*residual(n,m))
    end do
  end do
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
end subroutine rupdate_TVD2_RK
!
!###############################################################################
!
subroutine rupdate_TVD3_RK(nst,Apply_Filter,Apply_Limiter,iter)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : Limiter_Option,Filter_Option
  use flowvar, only : residual,usp,uoldsp,dtsp
  use filter_mod, only : filter_usp
  use module_limiters, only : resolution_indicator
  use module_limiters, only : limit_by_projecting
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst,iter
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: ac,bc
  !
  ! TVD 2nd-order/2-stage Runge-Kutta Coefficients
  !
  real(wp), parameter, dimension(1:3) :: alfa = (/ zero,three,half /)
  real(wp), parameter, dimension(1:3) :: beta = (/ one ,one4 ,two3 /)
  !
continue
  !
  nr = nme-nmb+1
  !
  bc = beta(nst)
  ac = alfa(nst)*bc
  !
 !call write_array(100,iter,nst,usp)
 !call write_array(300,iter,nst,uoldsp)
 !call write_array(500,iter,nst,residual)
 !call write_array(700,iter,nst,dtsp)
  !
  if (nr == 3) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,n_solpts
      usp(n,nec) = ac*uoldsp(n,nec) + bc*usp(n,nec) + bc*dtsp(n)*residual(n,nec)
      usp(n,nmx) = ac*uoldsp(n,nmx) + bc*usp(n,nmx) + bc*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = ac*uoldsp(n,nmy) + bc*usp(n,nmy) + bc*dtsp(n)*residual(n,nmy)
      usp(n,nmz) = ac*uoldsp(n,nmz) + bc*usp(n,nmz) + bc*dtsp(n)*residual(n,nmz)
      usp(n,nee) = ac*uoldsp(n,nee) + bc*usp(n,nee) + bc*dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,n_solpts
      usp(n,nec) = ac*uoldsp(n,nec) + bc*usp(n,nec) + bc*dtsp(n)*residual(n,nec)
      usp(n,nmx) = ac*uoldsp(n,nmx) + bc*usp(n,nmx) + bc*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = ac*uoldsp(n,nmy) + bc*usp(n,nmy) + bc*dtsp(n)*residual(n,nmy)
      usp(n,nee) = ac*uoldsp(n,nee) + bc*usp(n,nee) + bc*dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,n_solpts
      usp(n,nec) = ac*uoldsp(n,nec) + bc*usp(n,nec) + bc*dtsp(n)*residual(n,nec)
      usp(n,nmx) = ac*uoldsp(n,nmx) + bc*usp(n,nmx) + bc*dtsp(n)*residual(n,nmx)
      usp(n,nee) = ac*uoldsp(n,nee) + bc*usp(n,nee) + bc*dtsp(n)*residual(n,nee)
    end do
  end if
 !!$ACC PARALLEL LOOP GANG WORKER VECTOR COLLAPSE(2) &
 !!$ACC&              PRESENT(residual,dtsp,uoldsp,usp)
 !do m = 1,nq
 !  do n = 1,n_solpts
 !    usp(n,m) = beta(nst) * (alfa(nst)*uoldsp(n,m) + &
 !                            usp(n,m) + &
 !                            dtsp(n)*residual(n,m))
 !  end do
 !end do
  !
 !call write_array(200,iter,nst,usp)
 !call write_array(400,iter,nst,uoldsp)
 !call write_array(600,iter,nst,residual)
 !call write_array(800,iter,nst,dtsp)
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
end subroutine rupdate_TVD3_RK
!
!###############################################################################
!
subroutine rupdate_TVD4_RK(nst,Apply_Filter,Apply_Limiter)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : Limiter_Option,Filter_Option
  use flowvar, only : rssprk,residual,usp,uoldsp,usp_sum,dtsp
  use filter_mod, only : filter_usp
  use module_limiters, only : resolution_indicator
  use module_limiters, only : limit_by_projecting
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Scalars ..
  integer :: n,m
  !
  ! TVD 4th-order/5-stage Runge-Kutta Coefficients
  !
  real(wp), parameter, dimension(1:2,1:5) :: &
    beta = reshape( (/ 0.00000000000000_wp,0.39175222700395_wp, &
                       0.00000000000000_wp,0.36841059262959_wp, &
                       0.00000000000000_wp,0.25189177424738_wp, &
                       0.00000000000000_wp,0.54497475021237_wp, &
                       0.08460416338212_wp,0.22600748319395_wp  /), &
                    (/2,5/) )
  !
  real(wp), parameter, dimension(1:3,1:5) :: &
    alfa = reshape( (/ 1.00000000000000_wp,0.00000000000000_wp,0.0_wp, &
                       0.44437049406734_wp,0.55562950593266_wp,0.0_wp, &
                       0.62010185138540_wp,0.37989814861460_wp,0.0_wp, &
                       0.17807995410773_wp,0.82192004589227_wp,0.0_wp, &
                       0.00683325884039_wp,0.34833675773694_wp,1.0_wp  /), &
                    (/3,5/) )
  !
  real(wp), parameter, dimension(1:5) :: gama = (/ 0.00000000000000_wp, &
                                                    0.51723167208978_wp, &
                                                    0.12759831133288_wp, &
                                                    0.00000000000000_wp, &
                                                    0.00000000000000_wp  /)
  !
continue
  !
  if (nst == 1) usp_sum(:,:) = zero
  !
  do m = 1,nq
    do n = 1,n_solpts
      !
      rssprk(n,m,1) = beta(1,nst)*rssprk(n,m,1) + beta(2,nst)*residual(n,m)
      !
      usp(n,m) = alfa(1,nst)* uoldsp(n,m) + &
                 alfa(2,nst)*    usp(n,m) + &
                 alfa(3,nst)*usp_sum(n,m) + dtsp(n)*rssprk(n,m,1)
      !
      usp_sum(n,m) = usp_sum(n,m) + gama(nst)*usp(n,m)
      !
    end do
  end do
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
end subroutine rupdate_TVD4_RK
!
!###############################################################################
!
!subroutine initialize_time_ave()
!  !
!  !.. Use Statements ..
!  use geovar,  only : nr,n_solpts
!  use ovar,    only : time,ave_start_time
!  use ovar,    only : time_ave_vel_is_axisymm
!  use flowvar, only : time_ave_variables
!  use flowvar, only : uavesp
!  !
!  !.. Local Scalars ..
!  integer :: l,m,nv,nvar,ierr
!  character(len=1)  :: ta_char,v1,v2,v3
!  character(len=30) :: array_name
!  !
!  !.. Local Parameters ..
!  character(len=*), parameter :: pname = "initialize_time_ave"
!  !
!continue
!  !
!  call debug_timer(entering_procedure,pname)
!  !
!  ! Allocate and create the time_ave_variables array if it hasnt already
!  ! been allocated; otherwise, count the number of time-averaged variables it
!  ! contains. This character array identifies which time-averaged variables
!  ! are computed and written to the time-averaged solution/restart files.
!  !
!  if (.not. allocated(time_ave_variables)) then
!    !
!    v1 = "x"
!    v2 = merge( "r" , "y" , time_ave_vel_is_axisymm )
!    v3 = merge( "f" , "z" , time_ave_vel_is_axisymm )
!    !
!    ! Since the time_ave_variables array isnt allocated, this means that a
!    ! time-averaged restart file wasnt read and the start time of the
!    ! time-averaging needs to be initialized to the current time
!    !
!    ave_start_time = time
!    !
!   !nvar = merge(23,11,nr==3)
!   !allocate ( time_ave_variables(1:nvar) , stat=ierr , errmsg=error_message )
!   !call alloc_error(pname,"time_ave_variables",1,__LINE__,__FILE__, &
!   !                 ierr,error_message)
!    !
!    if (nr == 3) then
!      !
!      ! First 5 are the conserved variables
!      !
!      ! USING F2003 AUTO-REALLOCATION
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             "d" , "dx" , "dy" , "dz" , "e" ]
!      !
!      ! Next 5 are the primitive variables
!      ! (with temperature instead of density)
!      !
!      ! USING F2003 AUTO-REALLOCATION
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , &
!                             "t" , "x" , "y" , "z" , "p" ]
!      !
!      ! Remaining are derived variables
!      !
!      ! ALL OF THE BELOW ASSIGNMENTS USE F2003 AUTO-REALLOCATION TO
!      ! APPEND A NEW ARRAY ELEMENT ONTO THE time_ave_variables ARRAY
!      ! FOR EACH ADDITIONAL TIME AVERAGED VARIABLE
!      !
!      if (time_ave_vel_is_axisymm) then
!        time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                               time_ave_variables(:) , "d" // v2 ]
!        time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                               time_ave_variables(:) , "d" // v3 ]
!        time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                               time_ave_variables(:) , v2 ]
!        time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                               time_ave_variables(:) , v3 ]
!      end if
!      !
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , v1 // v1 ]
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , v2 // v2 ]
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , v3 // v3 ]
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , v1 // v2 ]
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , v1 // v3 ]
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , v2 // v3 ]
!      !
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , v1 // "t" ]
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , v2 // "t" ]
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , v3 // "t" ]
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , "dt" ]
!      !
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , "d" // v1 // v1 ]
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , "d" // v2 // v2 ]
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , "d" // v3 // v3 ]
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , "d" // v1 // v2 ]
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , "d" // v1 // v3 ]
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , "d" // v2 // v3 ]
!      !
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , "d" // v1 // "t" ]
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , "d" // v2 // "t" ]
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , "d" // v3 // "t" ]
!      !
!      nvar = size(time_ave_variables)
!      !
!    else
!      !
!      ! First 4 are the conserved variables
!      !
!      ! USING F2003 AUTO-REALLOCATION
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             "d" , "dx" , "dy" , "e" ]
!      !
!      ! Next 4 are the primitive variables
!      ! (with temperature instead of density)
!      !
!      ! USING F2003 AUTO-REALLOCATION
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , "t" , "x" , "y" , "p" ]
!      !
!      ! Remaining are derived variables
!      !
!      ! USING F2003 AUTO-REALLOCATION
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , "xx" ]
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , "yy" ]
!      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
!                             time_ave_variables(:) , "xy" ]
!      !
!      nvar = size(time_ave_variables)
!      !
!    end if
!    !
!  else
!    !
!    nvar = size(time_ave_variables)
!    !
!  end if
!  !
!  if (allocated(time_ave_var)) then
!    deallocate ( time_ave_var , stat=ierr , errmsg=error_message )
!    call alloc_error(pname,"time_ave_var",2,__LINE__,__FILE__, &
!                     ierr,error_message)
!  end if
!  !
!  allocate ( time_ave_var(1:nvar) , stat=ierr , errmsg=error_message )
!  call alloc_error(pname,"time_ave_var",1,__LINE__,__FILE__,ierr,error_message)
!  !
!  do m = 1,size(time_ave_variables)
!    !
!    time_ave_variables(m) = trim(adjustl(time_ave_variables(m)))
!    !
!    nv = len_trim(time_ave_variables(m))
!    !
!    allocate ( time_ave_var(m)%idx(1:nv) , source=0 , &
!               stat=ierr , errmsg=error_message )
!    write (array_name,1) m
!    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
!    !
!    do l = 1,nv
!      !
!      ta_char = uppercase( time_ave_variables(m)(l:l) )
!      !
!      select case (ta_char)
!        case ("D")
!          time_ave_var(m)%idx(l) = nec
!        case ("X")
!          time_ave_var(m)%idx(l) = nmx
!        case ("Y")
!          time_ave_var(m)%idx(l) = nmy
!        case ("Z")
!          time_ave_var(m)%idx(l) = nmz
!        case ("E")
!          time_ave_var(m)%idx(l) = nee
!        case ("P")
!          time_ave_var(m)%idx(l) = nee
!        case ("T")
!          time_ave_var(m)%idx(l) = nq+1
!        case ("R")
!          time_ave_var(m)%idx(l) = nq+2
!        case ("F")
!          time_ave_var(m)%idx(l) = nq+3
!        case ("K")
!          time_ave_var(m)%idx(l) = ntk
!        case ("W")
!          time_ave_var(m)%idx(l) = ntl
!        case default
!          ! ??? Should probably produce an error message and abort
!          time_ave_var(m)%idx(l) = nec
!      end select
!    end do
!  end do
!  !
!  ! Allocate the local module array used for computing
!  ! the current update to the time-averaged variables.
!  !
!  if (allocated(ave)) then
!    deallocate ( ave , stat=ierr , errmsg=error_message )
!    call alloc_error(pname,"ave",2,__LINE__,__FILE__,ierr,error_message)
!  end if
!  !
!  allocate ( ave(1:nvar) , source=zero , stat=ierr , errmsg=error_message )
!  call alloc_error(pname,"ave",1,__LINE__,__FILE__,ierr,error_message)
!  !
!  ! Allocate and initialize the time-averaged variables
!  ! solution array if it hasnt already been allocated.
!  !
!  if (.not. allocated(uavesp)) then
!    !
!    allocate ( uavesp(1:n_solpts,1:nvar) , source=zero , &
!               stat=ierr , errmsg=error_message )
!    call alloc_error(pname,"uavesp",1,__LINE__,__FILE__,ierr,error_message)
!    !$ACC UPDATE DEVICE(uavesp)
!    !
!    ! Initialize the time-averaged variables to the initial/restart
!    ! solution since they werent read from a time-averaged restart file.
!    !
!    call update_time_ave
!    !
!  end if
!  !
!  call debug_timer(leaving_procedure,pname)
!  !
!  ! Format Statements
!  !
!  1 format ("time_ave_var(",i0,")%idx")
!  !
!end subroutine initialize_time_ave
!!
!!###############################################################################
!!
!subroutine update_time_ave()
!  !
!  !.. Use Statements ..
!  use geovar,  only : xyz,n_solpts
!  use flowvar, only : usp,uavesp
!  use ovar,    only : time,prev_ave_time
!  use ovar,    only : ave_start_time,saved_ave_time
!  !
!  !.. Local Scalars ..
!  integer  :: m,n,nvar
!  integer  :: nc1,nc2,nm1,nm2
!  real(wp) :: a,b,y,z,st,ct,theta
!  !
!  !.. Local Arrays ..
!  real(wp) :: time_val(1:2)
!  real(wp) :: ave_var(1:nq+3)
!  !
!  !.. Local Parameters ..
!  character(len=*), parameter :: pname = "update_time_ave"
!  !
!continue
!  !
!  call debug_timer(entering_procedure,pname)
!  !
!  if (time == zero .or. time == ave_start_time) then
!    time_val(1) = zero
!    time_val(2) = one
!  else
!    time_val(1) = (prev_ave_time-ave_start_time) / (time-ave_start_time)
!    time_val(2) = one - time_val(1) ! to try and ensure that sum(time_val) = 1
!  end if
!  !
!  saved_ave_time = prev_ave_time
!  !
!  ! Broadcast the time_val array so hopefully all MPI processes
!  ! have the exact same values for the two averaging coefficients
!  ! (in case there might be roundoff differences between processes).
!  !
!  call mpi_bcast(time_val,size(time_val,kind=int_mpi),mpi_flttyp, &
!                 glb_root,MPI_COMM_WORLD,mpierr)
!  !
!  a = time_val(1)
!  b = time_val(2)
!  !
!  nvar = size(time_ave_var)
!  !
!  ! Adjust the coordinate and momentum indices for conversion to cylindrical
!  ! coordinates depending whether this is a 2D or 3D simulation
!  !
!  if (size(xyz,dim=1) == 2) then
!    nc1 = 1; nc2 = 2; nm1 = nmx; nm2 = nmy
!  else
!    nc1 = 2; nc2 = 3; nm1 = nmy; nm2 = nmz
!  end if
!  !
!  do n = 1,n_solpts
!    !
!    y = xyz(nc1,n)
!    z = xyz(nc2,n)
!    !
!    if (all([y,z] == zero)) then
!      !
!      st = zero
!      ct = zero
!      !
!    else
!      !
!      theta = atan2(z,y)
!      if (theta < zero) theta = theta + two*pi
!      !
!      st = sin(theta)
!      ct = cos(theta)
!      !
!    end if
!    !
!    ! Primitive variables
!    ave_var(1:nq) = usp2v_sp( usp(n,1:nq) )
!    !
!    ! Temperature
!    ave_var(nq+1) = temperature_pv_sp( ave_var(1:nq) )
!    !
!    ! Radial    instantaneous velocity
!    ave_var(nq+2) = st*ave_var(nm2) + ct*ave_var(nm1)
!    ! Azimuthal instantaneous velocity
!    ave_var(nq+3) = ct*ave_var(nm2) - st*ave_var(nm1)
!    !
!    ! First, update the time-averaged conservative variables
!    !
!    do m = 1,nq
!      uavesp(n,m) = a*uavesp(n,m) + b*usp(n,m)
!    end do
!    !
!    ! Next, update the time-averaged primitive variables
!    ! (with temperature instead of density)
!    !
!    uavesp(n,nq+1) = a*uavesp(n,nq+1) + b*ave_var(nq+1)
!    do m = 2,nq
!      uavesp(n,nq+m) = a*uavesp(n,nq+m) + b*ave_var(m)
!    end do
!    !
!    ! Update the remaining time-averaged variables
!    !
!    do m = 2*nq+1,nvar
!      uavesp(n,m) = a*uavesp(n,m) + b*product( ave_var(time_ave_var(m)%idx) )
!    end do
!    !
!  end do
!  !
!  ! Update prev_ave_time to the current time
!  !
!  prev_ave_time = time
!  !
!  call debug_timer(leaving_procedure,pname)
!  !
!end subroutine update_time_ave
!!
!!###############################################################################
!!
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
subroutine debug_timer(time_flag,message,processor,part_num,acc)
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: time_flag
  !
  character(len=*), intent(in) :: message
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: processor
  integer, optional, intent(in) :: part_num
  integer, optional, intent(in) :: acc
  !
#ifdef PROFILE_ON
  !.. Local Scalars ..
  real(wp), save :: beg_wall_time
  real(wp)       :: end_wall_time
  real(wp)       :: wall_time
  !
  integer, save :: procedure_counter = 0
  integer       :: i,hours,minutes,seconds,total_sec
  integer       :: output_cpu
  !
  character(len=300) :: out_message
  !
  !.. Local Arrays ..
  real(wp), save :: beg_procedure_time(10)
  real(wp), save :: accumulated_time(10)
  real(wp), save :: accumulation_beg_time(10)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: leading_space = "                    "
  !
continue
  !
  output_cpu = 0
  if (present(processor)) output_cpu = processor
  !
  if (time_flag == start_timer) then
    !
    if (present(part_num)) then
      write (out_message,6) message,part_num
    else
      write (out_message,5) message
    end if
    !
    i = 3 + procedure_counter*2
    !
    beg_wall_time = mpi_wtime()
   !beg_wall_time = dclock()
   !call cpu_time( beg_wall_time )
    !
    if (mypnum == output_cpu) then
      write (iout,1) leading_space(1:i),trim(adjustl(out_message))
    end if
    !
  else if (time_flag == stop_timer) then
    !
    end_wall_time = mpi_wtime()
   !end_wall_time = dclock()
   !call cpu_time( end_wall_time )
    wall_time = end_wall_time - beg_wall_time
    !
    total_sec = nint(wall_time)
    !
    hours   = total_sec / 3600
    minutes = (total_sec - hours*3600) / 60
    seconds = total_sec - hours*3600 - minutes*60
    !
    if (present(part_num)) then
      write (out_message,6) message,part_num
    else
      write (out_message,5) message
    end if
    !
    i = 3 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,2) leading_space(1:i),wall_time,trim(adjustl(out_message))
    end if
    !
    ! Reset beg_wall_time in case this procedure is called twice
    ! in a row with time_flag set to stop_timer
    !
    beg_wall_time = mpi_wtime()
   !beg_wall_time = dclock()
   !call cpu_time( beg_wall_time )
    !
  else if (time_flag == reset_accumulated_time) then
    !
    if (present(acc)) then
      accumulated_time(acc) = zero
    end if
    !
  else if (time_flag == report_accumulated_time) then
    !
    if (present(part_num)) then
      write (out_message,6) message,part_num
    else
      write (out_message,5) message
    end if
    !
    i = 3 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,2) leading_space(1:i),accumulated_time(acc), &
                     trim(adjustl(out_message))
    end if
    !
  else if (time_flag == start_accumulation) then
    !
    if (present(acc)) then
      accumulation_beg_time(acc) = mpi_wtime()
    end if
    !
  else if (time_flag == stop_accumulation) then
    !
    if (present(acc)) then
      wall_time = mpi_wtime() - accumulation_beg_time(acc)
      accumulated_time(acc) = accumulated_time(acc) + wall_time
    end if
    !
  else if (time_flag == entering_procedure) then
    !
    i = 1 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,3) leading_space(1:i),trim(adjustl(message))
    end if
    !
    procedure_counter = procedure_counter + 1
    !
    beg_procedure_time(procedure_counter) = mpi_wtime()
   !beg_procedure_time(procedure_counter) = dclock()
   !call cpu_time( beg_procedure_time(procedure_counter) )
    !
  else if (time_flag == leaving_procedure) then
    !
    end_wall_time = mpi_wtime()
   !end_wall_time = dclock()
   !call cpu_time( end_wall_time )
    wall_time = end_wall_time - beg_procedure_time(procedure_counter)
    !
    procedure_counter = max( 0 , procedure_counter-1 )
    !
    i = 1 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,4) leading_space(1:i),trim(adjustl(message)),wall_time
    end if
    !
  end if
  !
  flush (iout)
  !
  ! Format Statements
  !
  1 format (a,"## Starting timer for: ",a)
  2 format (a,"     -- TIME WAS",f12.6," seconds for ",a)
  3 format (a,"=> Entering procedure : ",a)
  4 format (a," <= Leaving procedure : ",a,"   after",f12.6," seconds")
  5 format (a)
  6 format (a," for partition # ",i0)
  !
#endif
  !
end subroutine debug_timer
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
