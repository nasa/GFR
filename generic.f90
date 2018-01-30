module generic_mod
  !
!#define PROFILE_ON
  !
  use module_kind_types, global_debug_timer => debug_timer
  use eqn_idx, only : nec,nmb,nme,nmx,nmy,nmz,nee,nq,ntk,ntl
  !
  implicit none
  !
  private
  !
  public :: init_generic_mod
  public :: write_iter_stats
  public :: write_iter_stats_gpu
  public :: getmaxent
  public :: compute_residual_error
  public :: compute_residual_error_gpu
  public :: compute_TaylorGreen_KEDR
  public :: compute_TaylorGreen_KEDR_gpu
  public :: generic_memory_usage
  !
  character(len=100), save :: iter_out_fmt
  !
  real(wp),    save :: Density_Maximum_Residual = -huge(zero)
  logical(lk), save :: vortex_has_collapsed = fals
  integer,     save :: final_countdown = 0
  !
  real(wp), save :: rl2mx = -huge(zero)
  real(wp), save :: rl2mx1 = -huge(zero)
  !
  ! TGV_KEDR : array used for the Taylor-Green vortex problem that
  !            contains the temporal evolution of various integrated
  !            quantities relating to the kinetic energy
  !
  integer,  save              :: iotg
  real(wp), save, allocatable :: tgv_kedr(:,:)
  real(wp), save, allocatable :: sum_kedr(:,:)
  real(wp), save              :: tgv_dimen(1:7) = zero
  !
  !$ACC DECLARE CREATE(tgv_kedr,sum_kedr,tgv_dimen)
  !
  integer,  save              :: iotg1
  real(wp), save, allocatable :: tgv_kedr1(:,:)
  real(wp), save, allocatable :: sum_kedr1(:,:)
  real(wp), save              :: tgv_dimen1(1:7) = zero
  !
contains
!
!###############################################################################
!
subroutine compute_residual_error_gpu()
  !
  !.. Use Statements ..
  use geovar,       only : nr,n_solpts,ncell,cell
  use ovar,         only : total_volume,uverr
  use ovar,         only : use_async_for_error
  use parallel_mod, only : cell_map
  !
  !.. Local Scalars ..
  integer  :: k,l,m,n,nc,n1,n2,l1,l2,l3,e1,e2,e3
  real(wp) :: vol_inv
  !
  !.. Local Arrays ..
 !real(wp), dimension(1:nq+1,1:5+nr,1:2)        :: error
  real(wp), dimension(1:nq+1,1:5+nr,1:2,1:ncpu) :: all_error
  !
  ! error(1:nq+1,:,:) = error value for each flow equation + extra variable
  ! error(:,1,:) = error for conserved variables [rho,rho(u,v,..),rhoe,entropy]
  ! error(:,2,:) = error for primitive variables [rho,u,v,..,pressure,vel mag]
  ! error(:,:,1) = L1    norm of total error for all solution points
  ! error(:,:,2) = L2%   norm of total error for all solution points
  ! error(:,:,3) = L-inf norm (max value) of error for all solution points
  ! error(:,:,4) = ID of grid cell containing maximum value of error
  ! error(:,:,5) = ID of solution point with maximum value of error
  ! error(:,:,6) = x-coordinate of sol-pt with max value of error
  ! error(:,:,7) = y-coordinate of sol-pt with max value of error (if 2D or 3D)
  ! error(:,:,8) = z-coordinate of sol-pt with max value of error (if 3D
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: cell_for_sp(:)
  real(wp), allocatable :: tunew(:,:)
  real(wp), allocatable :: tuold(:,:)
  real(wp), allocatable :: error(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_residual_error_gpu"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  l1 = 1
  l2 = min(2,nr)
  l3 = min(3,nr)
  !
  e1 = size(uverr,dim=1)
  e2 = size(uverr,dim=2)
  e3 = size(uverr,dim=3)
  !
  allocate ( error(1:e1,1:e2,1:e3) , source=zero )
  allocate ( tunew(1:n_solpts,1:nq+1) , source=zero )
  allocate ( tuold(1:n_solpts,1:nq+1) , source=zero )
  allocate ( cell_for_sp(1:n_solpts) , source=0 )
  !
  !$ACC DATA CREATE(error,tunew,tuold,cell_for_sp)
  !
  !$ACC PARALLEL LOOP GANG &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,cell_for_sp)
  do nc = 1,ncell
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    !$ACC LOOP VECTOR
    do n = n1,n2
      cell_for_sp(n) = nc
    end do
  end do
  !
  call get_unew_solpts_gpu(tunew)
  call get_uold_solpts_gpu(tuold)
  !
  if (use_async_for_error) then
    call compute_error_gpu_async(tunew,tuold,error(:,:,1))
  else
    call compute_error_gpu(tunew,tuold,error(:,:,1))
  end if
  !
  call get_vnew_solpts_gpu(tunew)
  call get_vold_solpts_gpu(tuold)
  !
  if (use_async_for_error) then
    call compute_error_gpu_async(tunew,tuold,error(:,:,2))
  else
    call compute_error_gpu(tunew,tuold,error(:,:,2))
  end if
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              VECTOR_LENGTH(32) &
  !$ACC&              PRESENT(cell,cell_map,error,cell_for_sp)
  do l = 1,2
    do m = 1,nq+1
      n = nint(error(m,5,l))
      nc = cell_for_sp(n)
      k = n - cell(nc)%beg_sp + 1
      error(m,4,l) = real(cell_map(mypnum+1)%loc_to_glb(nc),kind=wp)
      error(m,5,l) = real(k,kind=wp)
    end do
  end do
  !
  !$ACC UPDATE SELF(error)
  !
!#define DEBUG_WRITE_ITER_OUT
#ifdef DEBUG_WRITE_ITER_OUT
  do l = 1,size(error,dim=3)
    do n = 1,size(error,dim=2)
      do m = 1,size(error,dim=1)
        write (*,101) m,n,l,error(m,n,l)
      end do
    end do
  end do
  write (*,*)
#endif
  !
  ! Finish computing the various global norms for all variables
  !
  vol_inv = one/total_volume
  !
  if (ncpu > 1) then
    !
    ! Collect the error information from all processors
    !
    call mpi_allgather(error,    size(error,kind=int_mpi),mpi_flttyp, &
                       all_error,size(error,kind=int_mpi),mpi_flttyp, &
                       MPI_COMM_WORLD,mpierr)
    !
    ! Finish computing the various error norms across all processors
    !
    do l = 1,2
      ! l : 1 => conserved variables + entropy
      !     2 => primitive variables (temp,vel,pres) + velocity magnitude
      do m = 1,nq+1
        !
        ! L1 norm for the current variable
        uverr(m,1,l) =       sum(all_error(m,1,l,:)) * vol_inv
        ! L2 norm for the current variable
        uverr(m,2,l) = sqrt( sum(all_error(m,2,l,:)) * vol_inv )
        ! Infinity norm for the current variable
        uverr(m,3,l) = maxval(all_error(m,3,l,:))
        !
        ! Processor with max inf norm for the current variable
        k = maxloc(all_error(m,3,l,:),dim=1)
        ! ID of grid cell with max inf norm for the current variable
        uverr(m,4,l) = all_error(m,4,l,k)
        ! ID of sol point with max inf norm for the current variable
        uverr(m,5,l) = all_error(m,5,l,k)
        ! Coordinates of sol point with max inf norm for the current variable
        uverr(m,5+l1,l) = all_error(m,5+l1,l,k)
        uverr(m,5+l2,l) = all_error(m,5+l2,l,k)
        uverr(m,5+l3,l) = all_error(m,5+l3,l,k)
        !
      end do
      !
    end do
    !
  else
    !
    ! Copy all the error values into the uverr array
    !
    !$ACC PARALLEL LOOP COLLAPSE(3) &
    !$ACC&              VECTOR_LENGTH(32) &
    !$ACC&              PRESENT(uverr,error)
    do l = 1,e3
   !do l = 1,size(error,dim=3) ! DOESNT WORK
      do n = 1,e2
   !  do n = 1,size(error,dim=2) ! DOESNT WORK
        do m = 1,e1
   !    do m = 1,size(error,dim=1) ! DOESNT WORK
          uverr(m,n,l) = error(m,n,l)
        end do
      end do
    end do
    !
    ! Finalize the L1 and L2 norms
    !
    !$ACC PARALLEL LOOP COLLAPSE(2) &
    !$ACC&              VECTOR_LENGTH(32) &
    !$ACC&              PRESENT(uverr)
    do l = 1,e3
   !do l = 1,size(uverr,dim=3) ! DOESNT WORK
      do m = 1,e1
     !  do m = 1,size(uverr,dim=1) ! DOESNT WORK
        uverr(m,1,l) =       uverr(m,1,l) * vol_inv
        uverr(m,2,l) = sqrt( uverr(m,2,l) * vol_inv )
      end do
    end do
    !
  end if
  !$ACC UPDATE SELF(uverr)
#ifdef DEBUG_WRITE_ITER_OUT
  !$ACC PARALLEL LOOP SEQ PRESENT(uverr)
  do n = 1,size(uverr,dim=2)
    !$ACC LOOP SEQ
    do m = 1,size(uverr,dim=1)
      write (*,*) m,n,uverr(m,n,1)
    end do
  end do
  write (*,*)
  !$ACC PARALLEL LOOP SEQ PRESENT(uverr)
  do n = 1,size(uverr,dim=2)
    !$ACC LOOP SEQ
    do m = 1,size(uverr,dim=1)
      write (*,*) m,n,uverr(m,n,2)
    end do
  end do
#endif
  !
  !$ACC END DATA
  deallocate ( cell_for_sp )
  deallocate ( tuold )
  deallocate ( tunew )
  deallocate ( error )
  !
#ifdef DEBUG_WRITE_ITER_OUT
  write (*,*)
  do l = 1,size(uverr,dim=3)
    do n = 1,size(uverr,dim=2)
      do m = 1,size(uverr,dim=1)
        write (*,100) m,n,l,uverr(m,n,l)
      end do
    end do
  end do
  write (*,*)
  100 format ("uverr(",i0,",",i0,",",i0,") = ",es23.15)
  101 format ("error(",i0,",",i0,",",i0,") = ",es23.15)
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_residual_error_gpu
!
!###############################################################################
!
subroutine compute_error_gpu(tunew,tuold,error)
  !
  !.. Use Statements ..
  use geovar,      only : nr,n_solpts,xyz
  use metrics_mod, only : wtdjac
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: tunew(:,:)
  real(wp),    intent(in) :: tuold(:,:)
  real(wp), intent(inout) :: error(:,:)
  !
  !.. Local Scalars ..
  integer  :: k,m,n,l1,l2,l3
  real(wp) :: linf,lone,ltwo,new_val,old_val,linf_max
  real(wp) :: lone_sum,ltwo_sum
  !
  !.. Local Allocatble Arrays ..
  real(wp), allocatable :: linf_error(:)
  !
  ! error(1:nq+1,:,:) = error value for each flow equation + extra variable
  ! error(:,1) = L1    norm of total error for all solution points
  ! error(:,2) = L2%   norm of total error for all solution points
  ! error(:,3) = L-inf norm (max value) of error for all solution points
  ! error(:,4) = ID of grid cell containing maximum value of error
  ! error(:,5) = ID of solution point with maximum value of error
  ! error(:,6) = x-coordinate of sol-pt with max value of error
  ! error(:,7) = y-coordinate of sol-pt with max value of error (if 2D or 3D)
  ! error(:,8) = z-coordinate of sol-pt with max value of error (if 3D
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_error_gpu"
  !
continue
  !
  l1 = 1
  l2 = min(2,nr)
  l3 = min(3,nr)
  !
  call debug_timer(entering_procedure,pname)
  !
  allocate ( linf_error(1:n_solpts) , source=zero )
  !$ACC DATA CREATE(linf_error)
  !
  do m = 1,nq+1
    !
    lone_sum = zero
    ltwo_sum = zero
    !
    !$ACC PARALLEL LOOP REDUCTION(+:lone_sum,ltwo_sum) &
    !$ACC&              PRESENT(tunew,tuold,wtdjac,linf_error)
    do n = 1,n_solpts
      !
      new_val = tunew(n,m)
      old_val = tuold(n,m)
      !
      linf = abs( new_val - old_val )
      lone = linf*wtdjac(n)
      ltwo = linf*lone
      !
      lone_sum = lone_sum + lone
      ltwo_sum = ltwo_sum + ltwo
      !
      linf_error(n) = linf
      !
    end do
    !
    !$ACC UPDATE SELF(linf_error)
    !
    linf_max = -huge(zero)
    do n = 1,n_solpts
      if (linf_error(n) > linf_max) then
        linf_max = linf_error(n)
        k = n
      end if
    end do
   !k = maxloc(linf_error(1:n_solpts),dim=1)
    !
#ifdef DEBUG_WRITE_ITER_OUT
    write (iout,100) "lone_sum = ",lone_sum
    write (iout,100) "ltwo_sum = ",ltwo_sum
    100 format (a,es16.8)
#endif
    error(m,1) = lone_sum
    error(m,2) = ltwo_sum
    error(m,3) = linf_error(k)
    error(m,5) = real(k,kind=wp)
    error(m,5+l1) = xyz(l1,k)
    error(m,5+l2) = xyz(l2,k)
    error(m,5+l3) = xyz(l3,k)
    !
  end do
  !
  !$ACC UPDATE DEVICE(error)
  !
  !$ACC END DATA
  deallocate ( linf_error )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_error_gpu
!
!###############################################################################
!
subroutine compute_error_gpu_async(tunew,tuold,error)
  !
  !.. Use Statements ..
  use geovar,      only : nr,n_solpts,xyz
  use metrics_mod, only : wtdjac
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: tunew(:,:)
  real(wp),    intent(in) :: tuold(:,:)
  real(wp), intent(inout) :: error(:,:)
  !
  !.. Local Scalars ..
  integer  :: k,m,n,l1,l2,l3
  real(wp) :: lone,ltwo,linf_max
  !
  !.. Local Arrays ..
  real(wp) :: lone_sum(1:nq+1)
  real(wp) :: ltwo_sum(1:nq+1)
  !
  !.. Local Allocatble Arrays ..
  real(wp), allocatable :: linf_error(:,:)
  !
  ! error(1:nq+1,:,:) = error value for each flow equation + extra variable
  ! error(:,1) = L1    norm of total error for all solution points
  ! error(:,2) = L2%   norm of total error for all solution points
  ! error(:,3) = L-inf norm (max value) of error for all solution points
  ! error(:,4) = ID of grid cell containing maximum value of error
  ! error(:,5) = ID of solution point with maximum value of error
  ! error(:,6) = x-coordinate of sol-pt with max value of error
  ! error(:,7) = y-coordinate of sol-pt with max value of error (if 2D or 3D)
  ! error(:,8) = z-coordinate of sol-pt with max value of error (if 3D
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_error_gpu_async"
  !
continue
  !
  l1 = 1
  l2 = min(2,nr)
  l3 = min(3,nr)
  !
  call debug_timer(entering_procedure,pname)
  !
  allocate ( linf_error(1:n_solpts,1:nq+1) )
  !$ACC DATA CREATE(linf_error,lone_sum,ltwo_sum)
  !
  !$ACC PARALLEL LOOP COLLAPSE(2) &
  !$ACC&              PRESENT(tunew,tuold,linf_error)
  do m = 1,nq+1
    do n = 1,n_solpts
      linf_error(n,m) = abs( tunew(n,m) - tuold(n,m) )
    end do
  end do
  !
  !$ACC UPDATE SELF(linf_error)
  !
  !$ACC PARALLEL LOOP GANG ASYNC &
  !$ACC&              NUM_GANGS(6) &
  !$ACC&              VECTOR_LENGTH(512) &
  !$ACC&              PRESENT(wtdjac,linf_error,lone_sum,ltwo_sum)
  do m = 1,nq+1
    !
    lone = zero
    ltwo = zero
    !
    !$ACC LOOP REDUCTION(+:lone,ltwo)
    do n = 1,n_solpts
     !lone = linf_error(n,m)*wtdjac(n)
     !ltwo = linf_error(n,m)*lone
     !lone_sum = lone_sum + lone
     !ltwo_sum = ltwo_sum + ltwo
      lone = lone + wtdjac(n)*linf_error(n,m)
      ltwo = ltwo + wtdjac(n)*linf_error(n,m)*linf_error(n,m)
    end do
    !
    lone_sum(m) = lone
    ltwo_sum(m) = ltwo
    !
  end do
  !
  ! PERFORM THE FOLLOWING LOOP ON THE CPU WHILE THE
  ! ABOVE LOOP IS DONE ON THE GPU
  !
  do m = 1,nq+1
    !
    linf_max = -huge(zero)
    !
    do n = 1,n_solpts
      if (linf_error(n,m) > linf_max) then
        linf_max = linf_error(n,m)
        k = n
      end if
    end do
   !k = maxloc(linf_error(:,m),dim=1)
    !
   !error(m,3) = linf_error(k,m)
    error(m,3) = linf_max
    error(m,5) = real(k,kind=wp)
    error(m,5+l1) = xyz(l1,k)
    error(m,5+l2) = xyz(l2,k)
    error(m,5+l3) = xyz(l3,k)
    !
  end do
  !
  ! IF IT HASNT YET, WAIT FOR THE PREVIOUS GPU LOOP TO FINISH AND
  ! COPY THE lone_sum AND ltwo_sum ARRAYS FROM THE GPU TO THE CPU
  !
  !$ACC UPDATE SELF(lone_sum,ltwo_sum) WAIT
  !
  ! SAVE THE L1 AND L2 SUMS INTO THE ERROR ARRAY
  !
  do m = 1,nq+1
    error(m,1) = lone_sum(m)
    error(m,2) = ltwo_sum(m)
  end do
  !
  ! COPY THE ERROR ARRAY BACK TO THE GPU
  !
  !$ACC UPDATE DEVICE(error)
  !
  !$ACC END DATA
  deallocate ( linf_error )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_error_gpu_async
!
!###############################################################################
!
subroutine get_unew_solpts_gpu(tunew)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : usp,uavesp
  use ovar,    only : itestcase,pref_nd,gam
  use ovar,    only : output_time_averaging
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: tunew(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: rho,rhou,rhov,rhow,rhoe,rhoq2
  real(wp) :: gm1,gm1_over_pref_nd,pres
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_unew_solpts_gpu"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  gm1 = gam - one
  !
  if (output_time_averaging) then
    !
    if (nr == 3) then
      !$ACC PARALLEL LOOP PRESENT(uavesp,tunew)
      do n = 1,n_solpts
        tunew(n,nec) = uavesp(n,nec)
        tunew(n,nmx) = uavesp(n,nmx)
        tunew(n,nmy) = uavesp(n,nmy)
        tunew(n,nmz) = uavesp(n,nmz)
        tunew(n,nee) = uavesp(n,nee)
      end do
    else if (nr == 2) then
      !$ACC PARALLEL LOOP PRESENT(uavesp,tunew)
      do n = 1,n_solpts
        tunew(n,nec) = uavesp(n,nec)
        tunew(n,nmx) = uavesp(n,nmx)
        tunew(n,nmy) = uavesp(n,nmy)
        tunew(n,nee) = uavesp(n,nee)
      end do
    else
      !$ACC PARALLEL LOOP PRESENT(uavesp,tunew)
      do n = 1,n_solpts
        tunew(n,nec) = uavesp(n,nec)
        tunew(n,nmx) = uavesp(n,nmx)
        tunew(n,nee) = uavesp(n,nee)
      end do
    end if
    !
  else
    !
    if (nr == 3) then
      !$ACC PARALLEL LOOP PRESENT(usp,tunew)
      do n = 1,n_solpts
        tunew(n,nec) = usp(n,nec)
        tunew(n,nmx) = usp(n,nmx)
        tunew(n,nmy) = usp(n,nmy)
        tunew(n,nmz) = usp(n,nmz)
        tunew(n,nee) = usp(n,nee)
      end do
    else if (nr == 2) then
      !$ACC PARALLEL LOOP PRESENT(usp,tunew)
      do n = 1,n_solpts
        tunew(n,nec) = usp(n,nec)
        tunew(n,nmx) = usp(n,nmx)
        tunew(n,nmy) = usp(n,nmy)
        tunew(n,nee) = usp(n,nee)
      end do
    else
      !$ACC PARALLEL LOOP PRESENT(usp,tunew)
      do n = 1,n_solpts
        tunew(n,nec) = usp(n,nec)
        tunew(n,nmx) = usp(n,nmx)
        tunew(n,nee) = usp(n,nee)
      end do
    end if
    !
  end if
  !
  if (itestcase == Inviscid_Gaussian_Bump) then
    !
    gm1_over_pref_nd = gm1 / pref_nd
    !
    if (nr == 3) then
      !$ACC PARALLEL LOOP PRESENT(tunew)
      do n = 1,n_solpts
        !
        rho  = tunew(n,nec)
        rhou = tunew(n,nmx)
        rhov = tunew(n,nmy)
        rhow = tunew(n,nmz)
        rhoe = tunew(n,nee)
        !
        rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
        !
        pres = gm1_over_pref_nd * (rhoe - half*rhoq2/rho)
        !
        tunew(n,nq+1) = pres / (rho**gam)
        !
      end do
    else if (nr == 2) then
      !$ACC PARALLEL LOOP PRESENT(tunew)
      do n = 1,n_solpts
        !
        rho  = tunew(n,nec)
        rhou = tunew(n,nmx)
        rhov = tunew(n,nmy)
        rhoe = tunew(n,nee)
        !
        rhoq2 = rhou*rhou + rhov*rhov
        !
        pres = gm1_over_pref_nd * (rhoe - half*rhoq2/rho)
        !
        tunew(n,nq+1) = pres / (rho**gam)
        !
      end do
    else
      !$ACC PARALLEL LOOP PRESENT(tunew)
      do n = 1,n_solpts
        !
        rho  = tunew(n,nec)
        rhou = tunew(n,nmx)
        rhoe = tunew(n,nee)
        !
        rhoq2 = rhou*rhou
        !
        pres = gm1_over_pref_nd * (rhoe - half*rhoq2/rho)
        !
        tunew(n,nq+1) = pres / (rho**gam)
        !
      end do
    end if
    !
  else
    !
    if (nr == 3) then
      !$ACC PARALLEL LOOP PRESENT(tunew)
      do n = 1,n_solpts
        !
        rho  = tunew(n,nec)
        rhou = tunew(n,nmx)
        rhov = tunew(n,nmy)
        rhow = tunew(n,nmz)
        rhoe = tunew(n,nee)
        !
        rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
        !
        pres = gm1 * (rhoe - half*rhoq2/rho)
        !
        tunew(n,nq+1) = log(pres) - gam*log(rho)
        !
      end do
    else if (nr == 2) then
      !$ACC PARALLEL LOOP PRESENT(tunew)
      do n = 1,n_solpts
        !
        rho  = tunew(n,nec)
        rhou = tunew(n,nmx)
        rhov = tunew(n,nmy)
        rhoe = tunew(n,nee)
        !
        rhoq2 = rhou*rhou + rhov*rhov
        !
        pres = gm1 * (rhoe - half*rhoq2/rho)
        !
        tunew(n,nq+1) = log(pres) - gam*log(rho)
        !
      end do
    else
      !$ACC PARALLEL LOOP PRESENT(tunew)
      do n = 1,n_solpts
        !
        rho  = tunew(n,nec)
        rhou = tunew(n,nmx)
        rhoe = tunew(n,nee)
        !
        rhoq2 = rhou*rhou
        !
        pres = gm1 * (rhoe - half*rhoq2/rho)
        !
        tunew(n,nq+1) = log(pres) - gam*log(rho)
        !
      end do
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_unew_solpts_gpu
!
!###############################################################################
!
subroutine get_uold_solpts_gpu(tuold)
  !
  !.. Use Statements ..
  use ovar,    only : prev_ave_time,saved_ave_time,ave_start_time
  use ovar,    only : output_time_averaging,read_time_ave_restart
  use ovar,    only : itcur,itrst,bc_in,itestcase,gam,mms_opt
  use flowvar, only : usp,uoldsp,uavesp
  use geovar,  only : n_solpts,xyz
  use mms_mod, only : mms_solution_at_xyz_gpu
  use gpu_mod, only : vortex_solution_gpu
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: tuold(:,:)
  !
  !.. Local Scalars ..
  integer     :: n,nr
  real(wp)    :: a,b,gm1,pres,rhoq2
  real(wp)    :: rho,rhou,rhov,rhow,rhoe
  logical(lk) :: use_uavesp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_uold_solpts_gpu"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  gm1 = gam - one
  !
  if (mms_opt /= 0) then
    !
    ! Evaluate the exact analytical MMS solution at the error points
    !
    call mms_solution_at_xyz_gpu(1,n_solpts,xyz,tuold, &
                                 return_nondim=true, &
                                 return_cv=true)
    !
  else if (any(itestcase == Transport_Problems)) then
    !
    ! Evaluate the exact analytical solution to the vortex transport problem
    !
    call vortex_solution_gpu(1,n_solpts,xyz,tuold)
    !
  else if (itestcase == Freestream_Preservation) then
    !
    if (nr == 3) then
      !$ACC PARALLEL LOOP PRESENT(bc_in,tuold)
      do n = 1,n_solpts
        tuold(n,nec) = bc_in(0)%cv(nec)
        tuold(n,nmx) = bc_in(0)%cv(nmx)
        tuold(n,nmy) = bc_in(0)%cv(nmy)
        tuold(n,nmz) = bc_in(0)%cv(nmz)
        tuold(n,nee) = bc_in(0)%cv(nee)
      end do
    else if (nr == 2) then
      !$ACC PARALLEL LOOP PRESENT(bc_in,tuold)
      do n = 1,n_solpts
        tuold(n,nec) = bc_in(0)%cv(nec)
        tuold(n,nmx) = bc_in(0)%cv(nmx)
        tuold(n,nmy) = bc_in(0)%cv(nmy)
        tuold(n,nee) = bc_in(0)%cv(nee)
      end do
    else
      !$ACC PARALLEL LOOP PRESENT(bc_in,tuold)
      do n = 1,n_solpts
        tuold(n,nec) = bc_in(0)%cv(nec)
        tuold(n,nmx) = bc_in(0)%cv(nmx)
        tuold(n,nee) = bc_in(0)%cv(nee)
      end do
    end if
    !
  else
    !
    use_uavesp = fals
    !
    if (output_time_averaging) then
      use_uavesp = true
      if (any(saved_ave_time == [zero,ave_start_time]) .or. &
          (.not. read_time_ave_restart .and. itcur-itrst < 100)) then
        use_uavesp = fals
      end if
    end if
    !
    if (use_uavesp) then
      !
      ! Compute uavesp at the previous averaging time
      !
      a = (prev_ave_time-saved_ave_time) / (saved_ave_time-ave_start_time)
      b = one - a
      !
      if (nr == 3) then
        !$ACC PARALLEL LOOP PRESENT(usp,uavesp,tuold)
        do n = 1,n_solpts
          tuold(n,nec) = b * uavesp(n,nec) - a * usp(n,nec)
          tuold(n,nmx) = b * uavesp(n,nmx) - a * usp(n,nmx)
          tuold(n,nmy) = b * uavesp(n,nmy) - a * usp(n,nmy)
          tuold(n,nmz) = b * uavesp(n,nmz) - a * usp(n,nmz)
          tuold(n,nee) = b * uavesp(n,nee) - a * usp(n,nee)
        end do
      else if (nr == 2) then
        !$ACC PARALLEL LOOP PRESENT(usp,uavesp,tuold)
        do n = 1,n_solpts
          tuold(n,nec) = b * uavesp(n,nec) - a * usp(n,nec)
          tuold(n,nmx) = b * uavesp(n,nmx) - a * usp(n,nmx)
          tuold(n,nmy) = b * uavesp(n,nmy) - a * usp(n,nmy)
          tuold(n,nee) = b * uavesp(n,nee) - a * usp(n,nee)
        end do
      else
        !$ACC PARALLEL LOOP PRESENT(usp,uavesp,tuold)
        do n = 1,n_solpts
          tuold(n,nec) = b * uavesp(n,nec) - a * usp(n,nec)
          tuold(n,nmx) = b * uavesp(n,nmx) - a * usp(n,nmx)
          tuold(n,nee) = b * uavesp(n,nee) - a * usp(n,nee)
        end do
      end if
      !
    else
      !
      if (nr == 3) then
        !$ACC PARALLEL LOOP PRESENT(uoldsp,tuold)
        do n = 1,n_solpts
          tuold(n,nec) = uoldsp(n,nec)
          tuold(n,nmx) = uoldsp(n,nmx)
          tuold(n,nmy) = uoldsp(n,nmy)
          tuold(n,nmz) = uoldsp(n,nmz)
          tuold(n,nee) = uoldsp(n,nee)
        end do
      else if (nr == 2) then
        !$ACC PARALLEL LOOP PRESENT(uoldsp,tuold)
        do n = 1,n_solpts
          tuold(n,nec) = uoldsp(n,nec)
          tuold(n,nmx) = uoldsp(n,nmx)
          tuold(n,nmy) = uoldsp(n,nmy)
          tuold(n,nee) = uoldsp(n,nee)
        end do
      else
        !$ACC PARALLEL LOOP PRESENT(uoldsp,tuold)
        do n = 1,n_solpts
          tuold(n,nec) = uoldsp(n,nec)
          tuold(n,nmx) = uoldsp(n,nmx)
          tuold(n,nee) = uoldsp(n,nee)
        end do
      end if
      !
    end if
    !
  end if
  !
  if (itestcase == Inviscid_Gaussian_Bump) then
    !
    !$ACC PARALLEL LOOP PRESENT(tuold)
    do n = 1,n_solpts
      tuold(n,nq+1) = one
    end do
    !
  else
    !
    if (nr == 3) then
      !$ACC PARALLEL LOOP PRESENT(tuold)
      do n = 1,n_solpts
        !
        rho  = tuold(n,nec)
        rhou = tuold(n,nmx)
        rhov = tuold(n,nmy)
        rhow = tuold(n,nmz)
        rhoe = tuold(n,nee)
        !
        rhoq2 = rhou*rhou + rhov*rhov + rhow*rhow
        !
        pres = gm1 * (rhoe - half*rhoq2/rho)
        !
        tuold(n,nq+1) = log(pres) - gam*log(rho)
        !
      end do
    else if (nr == 2) then
      !$ACC PARALLEL LOOP PRESENT(tuold)
      do n = 1,n_solpts
        !
        rho  = tuold(n,nec)
        rhou = tuold(n,nmx)
        rhov = tuold(n,nmy)
        rhoe = tuold(n,nee)
        !
        rhoq2 = rhou*rhou + rhov*rhov
        !
        pres = gm1 * (rhoe - half*rhoq2/rho)
        !
        tuold(n,nq+1) = log(pres) - gam*log(rho)
        !
      end do
    else
      !$ACC PARALLEL LOOP PRESENT(tuold)
      do n = 1,n_solpts
        !
        rho  = tuold(n,nec)
        rhou = tuold(n,nmx)
        rhoe = tuold(n,nee)
        !
        rhoq2 = rhou*rhou
        !
        pres = gm1 * (rhoe - half*rhoq2/rho)
        !
        tuold(n,nq+1) = log(pres) - gam*log(rho)
        !
      end do
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_uold_solpts_gpu
!
!###############################################################################
!
subroutine get_vnew_solpts_gpu(tvnew)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : uavesp
  use ovar,    only : gam,rgasref_nd
  use ovar,    only : output_time_averaging
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: tvnew(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nr
  real(wp) :: velx,vely,velz,q2,rho,rho_inv
  real(wp) :: rhoe,pres,temp,rgas_inv,gm1
  !
  real(wp) :: rhou,rhov,rhow,rhoq2,temp1,temp2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_vnew_solpts_gpu"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  gm1 = gam - one
  !
  ! Get the primitive variables of the current solution at the error points
  !
  if (output_time_averaging) then
    if (nr == 3) then
      !$ACC PARALLEL LOOP PRESENT(uavesp,tvnew)
      do n = 1,n_solpts
        !
        velx = uavesp(n,nq+nmx)
        vely = uavesp(n,nq+nmy)
        velz = uavesp(n,nq+nmz)
        !
        tvnew(n,nec) = uavesp(n,nq+nec)
        tvnew(n,nmx) = velx
        tvnew(n,nmy) = vely
        tvnew(n,nmz) = velz
        tvnew(n,nee) = uavesp(n,nq+nee)
        tvnew(n,nq+1) = sqrt(velx*velx + vely*vely + velz*velz)
        !
      end do
    else if (nr == 2) then
      !$ACC PARALLEL LOOP PRESENT(uavesp,tvnew)
      do n = 1,n_solpts
        !
        velx = uavesp(n,nq+nmx)
        vely = uavesp(n,nq+nmy)
        !
        tvnew(n,nec) = uavesp(n,nq+nec)
        tvnew(n,nmx) = velx
        tvnew(n,nmy) = vely
        tvnew(n,nee) = uavesp(n,nq+nee)
        tvnew(n,nq+1) = sqrt(velx*velx + vely*vely)
        !
      end do
    else
      !$ACC PARALLEL LOOP PRESENT(uavesp,tvnew)
      do n = 1,n_solpts
        !
        velx = uavesp(n,nq+nmx)
        !
        tvnew(n,nec) = uavesp(n,nq+nec)
        tvnew(n,nmx) = velx
        tvnew(n,nee) = uavesp(n,nq+nee)
        tvnew(n,nq+1) = sqrt(velx*velx)
        !
      end do
    end if
  else
    !
    rgas_inv = one/rgasref_nd
    !
    if (nr == 3) then
      !$ACC PARALLEL LOOP PRESENT(tvnew)
      do n = 1,n_solpts
        !
        rho = tvnew(n,nec)
        !
        rho_inv = one/rho
        !
        velx = tvnew(n,nmx)*rho_inv
        vely = tvnew(n,nmy)*rho_inv
        velz = tvnew(n,nmz)*rho_inv
        rhoe = tvnew(n,nee)
        !
        q2 = velx*velx + vely*vely + velz*velz
        !
        pres = gm1 * (rhoe - half*rho*q2)
        temp = pres*rho_inv*rgas_inv
        !
        tvnew(n,nec) = temp
        tvnew(n,nmx) = velx
        tvnew(n,nmy) = vely
        tvnew(n,nmz) = velz
        tvnew(n,nee) = pres
        tvnew(n,nq+1) = sqrt(q2)
        !
      end do
    else if (nr == 2) then
      !$ACC PARALLEL LOOP PRESENT(tvnew)
      do n = 1,n_solpts
        !
        rho = tvnew(n,nec)
        !
        rho_inv = one/rho
        !
        velx = tvnew(n,nmx)*rho_inv
        vely = tvnew(n,nmy)*rho_inv
        rhoe = tvnew(n,nee)
        !
        q2 = velx*velx + vely*vely
        !
        pres = gm1 * (rhoe - half*rho*q2)
        temp = pres*rho_inv*rgas_inv
        !
        tvnew(n,nec) = temp
        tvnew(n,nmx) = velx
        tvnew(n,nmy) = vely
        tvnew(n,nee) = pres
        tvnew(n,nq+1) = sqrt(q2)
        !
      end do
    else
      !$ACC PARALLEL LOOP PRESENT(tvnew)
      do n = 1,n_solpts
        !
        rho = tvnew(n,nec)
        !
        rho_inv = one/rho
        !
        velx = tvnew(n,nmx)*rho_inv
        rhoe = tvnew(n,nee)
        !
        q2 = velx*velx
        !
        pres = gm1 * (rhoe - half*rho*q2)
        temp = pres*rho_inv*rgas_inv
        !
        tvnew(n,nec) = temp
        tvnew(n,nmx) = velx
        tvnew(n,nee) = pres
        tvnew(n,nq+1) = sqrt(q2)
        !
      end do
    end if
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_vnew_solpts_gpu
!
!###############################################################################
!
subroutine get_vold_solpts_gpu(tvold)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : usp,uavesp
  use ovar,    only : prev_ave_time,saved_ave_time,ave_start_time
  use ovar,    only : output_time_averaging,read_time_ave_restart
  use ovar,    only : itcur,itrst,gam,rgasref_nd
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: tvold(:,:)
  !
  !.. Local Scalars ..
  integer     :: n,nr
  real(wp)    :: a,b,velx,vely,velz,q2,rho,rho_inv
  real(wp)    :: rhoe,pres,temp,rgas_inv,gm1
  logical(lk) :: use_uavesp
  !
  real(wp) :: rhou,rhov,rhow,rhoq2,temp1,temp2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_vold_solpts_gpu"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  gm1 = gam - one
  rgas_inv = one/rgasref_nd
  !
  ! Get the primitive variables of the current solution at the error points
  !
  use_uavesp = fals
  !
  if (output_time_averaging) then
    use_uavesp = true
    if (any(saved_ave_time == [zero,ave_start_time]) .or. &
        (.not. read_time_ave_restart .and. itcur-itrst < 100)) then
      use_uavesp = fals
    end if
  end if
  !
  if (use_uavesp) then
    !
    ! Compute uavesp at the previous averaging time
    !
    a = (prev_ave_time-saved_ave_time) / (saved_ave_time-ave_start_time)
    b = one - a
    !
    if (nr == 3) then
      !$ACC PARALLEL LOOP PRESENT(usp,uavesp,tvold)
      do n = 1,n_solpts
        !
        rho = usp(n,nec)
        !
        rho_inv = one/rho
        !
        velx = usp(n,nmx)*rho_inv
        vely = usp(n,nmy)*rho_inv
        velz = usp(n,nmz)*rho_inv
        rhoe = usp(n,nee)
        !
        q2 = velx*velx + vely*vely + velz*velz
        !
        pres = gm1 * (rhoe - half*rho*q2)
        temp = pres*rho_inv*rgas_inv
        !
        velx = b * uavesp(n,nq+nmx) - a * velx
        vely = b * uavesp(n,nq+nmy) - a * vely
        velz = b * uavesp(n,nq+nmz) - a * velz
        !
        q2 = b * sqrt(velx*velx + vely*vely + velz*velz) - a * sqrt(q2)
        !
        tvold(n,nec) = b * uavesp(n,nq+nec) - a * temp
        tvold(n,nmx) = velx
        tvold(n,nmy) = vely
        tvold(n,nmz) = velz
        tvold(n,nee) = b * uavesp(n,nq+nee) - a * pres
        tvold(n,nq+1) = q2
        !
      end do
    else if (nr == 2) then
      !$ACC PARALLEL LOOP PRESENT(usp,uavesp,tvold)
      do n = 1,n_solpts
        !
        rho = usp(n,nec)
        !
        rho_inv = one/rho
        !
        velx = usp(n,nmx)*rho_inv
        vely = usp(n,nmy)*rho_inv
        rhoe = usp(n,nee)
        !
        q2 = velx*velx + vely*vely
        !
        pres = gm1 * (rhoe - half*rho*q2)
        temp = pres*rho_inv*rgas_inv
        !
        velx = b * uavesp(n,nq+nmx) - a * velx
        vely = b * uavesp(n,nq+nmy) - a * vely
        !
        q2 = b * sqrt(velx*velx + vely*vely) - a * sqrt(q2)
        !
        tvold(n,nec) = b * uavesp(n,nq+nec) - a * temp
        tvold(n,nmx) = velx
        tvold(n,nmy) = vely
        tvold(n,nee) = b * uavesp(n,nq+nee) - a * pres
        tvold(n,nq+1) = q2
        !
      end do
    else
      !$ACC PARALLEL LOOP PRESENT(usp,uavesp,tvold)
      do n = 1,n_solpts
        !
        rho = usp(n,nec)
        !
        rho_inv = one/rho
        !
        velx = usp(n,nmx)*rho_inv
        rhoe = usp(n,nee)
        !
        q2 = velx*velx
        !
        pres = gm1 * (rhoe - half*rho*q2)
        temp = pres*rho_inv*rgas_inv
        !
        velx = b * uavesp(n,nq+nmx) - a * velx
        !
        q2 = b * sqrt(velx*velx) - a * sqrt(q2)
        !
        tvold(n,nec) = b * uavesp(n,nq+nec) - a * temp
        tvold(n,nmx) = velx
        tvold(n,nee) = b * uavesp(n,nq+nee) - a * pres
        tvold(n,nq+1) = q2
        !
      end do
    end if
    !
  else
    !
    if (nr == 3) then
      !$ACC PARALLEL LOOP PRESENT(tvold)
      do n = 1,n_solpts
        !
        rho = tvold(n,nec)
        !
        rho_inv = one/rho
        !
        velx = tvold(n,nmx)*rho_inv
        vely = tvold(n,nmy)*rho_inv
        velz = tvold(n,nmz)*rho_inv
        rhoe = tvold(n,nee)
        !
        q2 = velx*velx + vely*vely + velz*velz
        !
        pres = gm1 * (rhoe - half*rho*q2)
        temp = pres*rho_inv*rgas_inv
        !
        tvold(n,nec) = temp
        tvold(n,nmx) = velx
        tvold(n,nmy) = vely
        tvold(n,nmz) = velz
        tvold(n,nee) = pres
        tvold(n,nq+1) = sqrt(q2)
        !
      end do
    else if (nr == 2) then
      !$ACC PARALLEL LOOP PRESENT(tvold)
      do n = 1,n_solpts
        !
        rho = tvold(n,nec)
        !
        rho_inv = one/rho
        !
        velx = tvold(n,nmx)*rho_inv
        vely = tvold(n,nmy)*rho_inv
        rhoe = tvold(n,nee)
        !
        q2 = velx*velx + vely*vely
        !
        pres = gm1 * (rhoe - half*rho*q2)
        temp = pres*rho_inv*rgas_inv
        !
        tvold(n,nec) = temp
        tvold(n,nmx) = velx
        tvold(n,nmy) = vely
        tvold(n,nee) = pres
        tvold(n,nq+1) = sqrt(q2)
        !
      end do
    else
      !$ACC PARALLEL LOOP PRESENT(tvold)
      do n = 1,n_solpts
        !
        rho = tvold(n,nec)
        !
        rho_inv = one/rho
        !
        velx = tvold(n,nmx)*rho_inv
        rhoe = tvold(n,nee)
        !
        q2 = velx*velx
        !
        pres = gm1 * (rhoe - half*rho*q2)
        temp = pres*rho_inv*rgas_inv
        !
        tvold(n,nec) = temp
        tvold(n,nmx) = velx
        tvold(n,nee) = pres
        tvold(n,nq+1) = sqrt(q2)
        !
      end do
    end if
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_vold_solpts_gpu
!
!###############################################################################
!
subroutine compute_residual_error()
  !
  !.. Use Statements ..
  use order_mod, only : maxEP
  use geovar,  only : nr,ncell,cell
  use ovar,    only : total_volume,itestcase
  use ovar,    only : uverr1
  use parallel_mod, only : cell_map
  use interpolation_mod, only : outerp
  !
  !.. Local Scalars ..
  integer  :: this_cell,n1,n2,np,k,m,n,l
  integer  :: this_geom,this_order,l1,l2,l3
  real(wp) :: sn,se
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxEP) :: wtdjac
  real(wp), dimension(1:nq,1:maxEP) :: uold,vold
  real(wp), dimension(1:nq,1:maxEP) :: unew,vnew
  real(wp), dimension(1:nr,1:maxEP) :: xyz_err
  real(wp), dimension(1:nq+1,1:2) :: ett
  real(wp), dimension(1:nq+1,1:2,1:5+nr) :: error
  real(wp), dimension(1:nq+1,1:2,1:5+nr,1:ncpu) :: all_error
  !
  ! error(1:nq+1,:,:) = error value for each flow equation + extra variable
  ! error(:,1,:) = error for conserved variables [rho,rho(u,v,..),rhoe,entropy]
  ! error(:,2,:) = error for primitive variables [rho,u,v,..,pressure,vel mag]
  ! error(:,:,1) = L1    norm of total error for all solution points
  ! error(:,:,2) = L2%   norm of total error for all solution points
  ! error(:,:,3) = L-inf norm (max value) of error for all solution points
  ! error(:,:,4) = ID of grid cell containing maximum value of error
  ! error(:,:,5) = ID of solution point with maximum value of error
  ! error(:,:,6) = x-coordinate of sol-pt with max value of error
  ! error(:,:,7) = y-coordinate of sol-pt with max value of error (if 2D or 3D)
  ! error(:,:,8) = z-coordinate of sol-pt with max value of error (if 3D
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_residual_error"
  !
continue
  !
  l1 = 1
  l2 = min(2,nr)
  l3 = min(3,nr)
  !
  call debug_timer(entering_procedure,pname)
  !
  error(:,:,:) = zero  ! intialize everything to zero
  error(:,:,3) = -eps6 ! make this slightly negative in case all error is zero
  !
  do this_cell = 1,ncell
    !
    n1 = cell(this_cell)%beg_sp
    n2 = cell(this_cell)%end_sp
    np = n2-n1+1
    !
    this_geom = cell(this_cell)%geom
    this_order = cell(this_cell)%order
    !
    if (allocated(outerp(this_geom,this_order)%error)) then
      if (allocated(outerp(this_geom,this_order)%error%mat)) then
        np = size(outerp(this_geom,this_order)%error%mat,dim=1)
      end if
    end if
    !
    call collect_error_arrays(this_cell,unew(1:nq,1:np),uold(1:nq,1:np), &
                              wtdjac(1:np),xyz_err(1:nr,1:np), &
                              vnew(1:nq,1:np),vold(1:nq,1:np))
    !
    ! Now loop through the solution points of this cell to compute the
    ! error between the solutions of the current and previous time steps
    !
    do k = 1,np
      !
      ! Write out an error to stdout if density or pressure is negative
      !
      if (any([unew(nec,k),vnew(nee,k)] <= zero)) then
        write (iout,1) mypnum,cell_map(mypnum+1)%loc_to_glb(this_cell),k, &
                       unew(nec,k),vnew(nee,k),vnew(nec,k)
      end if
      !
      ! Error for conserved variables at the current solution point
      !
      if (itestcase == Inviscid_Gaussian_Bump) then
        sn = entropy_cv_sp( unew(:,k) , log_opt=fals , use_pref_nd=true )
        se = one
      else
        sn = entropy_cv_sp( unew(:,k) )
        se = entropy_cv_sp( uold(:,k) )
      end if
      !
      ett(1:nq,1) = abs( unew(:,k) - uold(:,k) )
      ett(nq+1,1) = abs(sn - se)
      !
      ! Error for primitive variables at the current solution point
      !
      ett(1:nq,2) = abs( vnew(:,k) - vold(:,k) )
     !ett(nq+1,2) = norm2( ett(nmb:nme,2) )
      ett(nq+1,2) = abs( norm2(vnew(nmb:nme,k)) - norm2(vold(nmb:nme,k)) )
      !
      ! Now add the error for all variables at the current solution
      ! point to the running integral sum for the L1 and L2 error norms
      !
      error(:,:,1) = error(:,:,1) + wtdjac(k)*ett(:,:)
      error(:,:,2) = error(:,:,2) + wtdjac(k)*ett(:,:)*ett(:,:)
      !
      ! If any of these errors are greater than the previous maximum
      ! values, save the error value and record this location
      !
      where (abs(ett(:,:)) > error(:,:,3))
        error(:,:,3) = abs(ett(:,:))
        error(:,:,4) = real(cell_map(mypnum+1)%loc_to_glb(this_cell),kind=wp)
        error(:,:,5) = real(k,kind=wp)
        error(:,:,5+l1) = xyz_err(l1,k)
        error(:,:,5+l2) = xyz_err(l2,k)
        error(:,:,5+l3) = xyz_err(l3,k)
      end where
      !
    end do
    !
  end do
  !
  ! Collect the error information from all processors
  !
  if (ncpu > 1) then
    call mpi_allgather(error,    size(error,kind=int_mpi),mpi_flttyp, &
                       all_error,size(error,kind=int_mpi),mpi_flttyp, &
                       MPI_COMM_WORLD,mpierr)
  else
    all_error(:,:,:,1) = error(:,:,:)
  end if
  !
  ! Finish computing the various global norms for all variables
  !
  do n = 1,2
    !
    ! n : 1 => conserved variables + entropy
    !     2 => primitive variables (temp,vel,pres) + velocity magnitude
    !
   !! L1 norm for the current set of variables
   !uverr1(:,1,n) =       sum(all_error(:,n,1,:),dim=2) / total_volume
   !! L2 norm for the current set of variables
   !uverr1(:,2,n) = sqrt( sum(all_error(:,n,2,:),dim=2) / total_volume )
   !! Infinity norm for the current set of variables
   !uverr1(:,3,n) = maxval(all_error(:,n,3,:),dim=2)
    !
    do m = 1,nq+1
      !
      ! L1 norm for the current variable
      uverr1(m,1,n) =       sum(all_error(m,n,1,:)) / total_volume
      ! L2 norm for the current variable
      uverr1(m,2,n) = sqrt( sum(all_error(m,n,2,:)) / total_volume )
      ! Infinity norm for the current variable
      uverr1(m,3,n) = maxval(all_error(m,n,3,:))
      !
      ! Processor with max inf norm for the current variable
      k = maxloc(all_error(m,n,3,:),dim=1)
      ! ID of grid cell with max inf norm for the current variable
      uverr1(m,4,n) = all_error(m,n,4,k)
      ! ID of sol point with max inf norm for the current variable
      uverr1(m,5,n) = all_error(m,n,5,k)
      ! Coordinates of sol point with max inf norm for the current variable
      uverr1(m,6:5+nr,n) = all_error(m,n,6:5+nr,k)
      !
    end do
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("ERROR DETECTED! CPU #",i0," =>  cell=",i0,", SP=",i0," :",/, &
            8x,"rho=",es13.6,", pres=",es13.6,", temp=",es13.6)
  !
end subroutine compute_residual_error
!
!###############################################################################
!
subroutine write_iter_stats_gpu(walltime_expiring,exit_main_loop)
  !
  !.. Use Statements ..
  use geovar,  only : global_pinc_ptr
  use ovar,    only : itcur,d_t,cfl,time,minimum_timesteps,time_ref
  use ovar,    only : convergence_reached,itestcase,Period_Timesteps
  use ovar,    only : convergence_order_max_res,convergence_order_abs_res
  use ovar,    only : uverr,this_is_final_timestep,iter_out_interval
  use ovar,    only : output_time_averaging,time_scaling_factor,ave_start_time
  !
  use pbs_mod, only : pbs_remaining_walltime
  !
  !.. Formal Arguments ..
  logical(lk), intent(inout) :: walltime_expiring
  logical(lk), intent(  out) :: exit_main_loop
  !
  !.. Local Scalars ..
  integer  :: hrs_left,min_left,sec_left
  integer  :: nloc,kloc
  real(wp) :: xentmax,yentmax
  real(wp) :: entL2
  real(wp) :: wall_time_left,out_time
  real(wp) :: orders_converged,output_L2val
  real(wp) :: rltwo
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_iter_stats_gpu"
  !
continue
  !
  exit_main_loop = fals
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Get the density and entropy residual stats from uverr
  !
  rltwo = uverr(nec,2,1)
  entL2 = uverr(nq+1,2,1)
  !
  nloc = nint( uverr(nec,4,1) )
  kloc = nint( uverr(nec,5,1) )
  !
  ! Coordinates Location of the max inf norm for entropy
  xentmax = uverr(nq+1,6,1)
  yentmax = uverr(nq+1,7,1)
  !
  ! Maximum L2 norm of density for the current simulation
  rl2mx = max( rltwo , rl2mx )
  !
  if (any(itestcase == Transport_Problems)) then
    !
    output_L2val = entL2
    convergence_reached = fals
    !
    ! Check for a collapse of the vortex, or update the final
    ! countdown if a collapse has already been detected
    !
    if (vortex_has_collapsed) then
      !
      ! Need to subtract from final_countdown the number of time steps
      ! that have been completed since the last call to write_iter_stats_gpu
      !
      final_countdown = max(final_countdown - iter_out_interval,0)
      !
    else if (rltwo > zero) then
      !
      if (log10(rltwo) > Density_Maximum_Residual) then
        Density_Maximum_Residual = log10(rltwo)
      end if
      if (Density_Maximum_Residual > -two) then
        vortex_has_collapsed = true
        final_countdown = 26*Period_Timesteps - mod(itcur,Period_Timesteps)
      end if
      !
    end if
    !
  else
    !
    orders_converged = log10(max(rl2mx,eps17)) - log10(max(rltwo,eps17))
    output_L2val = orders_converged
    !
    ! Check for convergence
    !
    if (itcur > minimum_timesteps) then
      convergence_reached = (-orders_converged <= convergence_order_max_res) &
                            .or. (log10(rltwo) <= convergence_order_abs_res)
      if (itestcase == Taylor_Green_Vortex) then
        convergence_reached = fals
      end if
    end if
    !
  end if
  !
  ! Output the results for the current time step
  !
  call pbs_remaining_walltime(walltime_expiring,wall_time_left)
  !
  if (mypnum == glb_root) then
    !
    out_time = time*time_ref*time_scaling_factor
    if (output_time_averaging) then
      out_time = out_time - ave_start_time*time_ref*time_scaling_factor
    end if
#ifdef PBS_ENV
    !
    sec_left = abs(nint(wall_time_left))
    !
    hrs_left = sec_left/3600
    min_left = (sec_left - hrs_left*3600)/60
    sec_left = sec_left - hrs_left*3600 - min_left*60
    ! Make hrs_left negative if wall_time_left happens to be negative
    hrs_left = sign(hrs_left,int(wall_time_left))
    !
#ifdef DEBUG_WRITE_ITER_OUT
    100 format (a,"'",a,"'")
    101 format (a,i0)
    102 format (a,es23.15)
    write (*,100) "iter_out_fmt = ",trim(adjustl(iter_out_fmt))
    write (*,101) "itcur = ",itcur
    write (*,102) "d_t = ",d_t
    write (*,102) "out_time = ",out_time
    write (*,102) "cfl = ",cfl
    write (*,102) "rltwo = ",rltwo
    write (*,102) "output_L2val = ",output_L2val
    write (*,101) "nloc = ",nloc
    write (*,101) "kloc = ",kloc
    write (*,101) "global_pinc_ptr(nloc) = ",global_pinc_ptr(nloc)
    write (*,102) "xentmax = ",xentmax
    write (*,102) "yentmax = ",yentmax
    write (*,101) "hrs_left = ",hrs_left
    write (*,101) "min_left = ",min_left
    write (*,101) "sec_left = ",sec_left
#endif
    !
   !write (iout,iter_out_fmt) itcur,d_t*time_ref,time*time_ref, &
    write (iout,iter_out_fmt) itcur,d_t,out_time, &
                              cfl,rltwo,output_L2val,nloc, &
                              global_pinc_ptr(nloc)+kloc, &
                              xentmax,yentmax, &
                              hrs_left,min_left,sec_left
    !
#else
    !
   !write (iout,iter_out_fmt) itcur,d_t*time_ref,time*time_ref, &
    write (iout,iter_out_fmt) itcur,d_t,out_time, &
                              cfl,rltwo,output_L2val,nloc, &
                              global_pinc_ptr(nloc)+kloc, &
                              xentmax,yentmax
    !
    !
#endif
    !
    if (convergence_reached) then
      if (-orders_converged <= convergence_order_max_res) write (iout,3)
      if (log10(rltwo) <= convergence_order_abs_res) write (iout,4)
    end if
    !
    if (walltime_expiring) write (iout,2)
    !
    if (this_is_final_timestep) write (iout,5)
    !
  end if
  !
  ! If a vortex collapse has been detected and final_countdown has reached 0,
  ! activate convergence_reached so that the main_loop can be exited when
  ! we return from this subroutine to the main program.
  !
  if (vortex_has_collapsed .and. final_countdown==0) then
    if (mypnum == glb_root) write (iout,6)
    convergence_reached = true
  end if
  !
  if (convergence_reached) exit_main_loop = true
  if (this_is_final_timestep) exit_main_loop = true
  if (walltime_expiring) exit_main_loop = true
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  2 format (/,"*******************************************************",/, &
              " THE SIMULATION WALL TIME REQUESTED IS EXPIRING SOON!",/, &
              " STOPPING EXECUTION AND DUMPING THE FINAL RESULTS!",/, &
              "*******************************************************",/)
  3 format (/,"*******************************************************",/, &
              " THE CONVERGENCE CRITERION FOR THE ORDERS OF",/, &
              " MAGNITUDE REDUCTION OF THE RESIDUAL ERROR",/, &
              " RELATIVE TO THE MAXIMUM RESIDUAL ERROR HAS BEEN MET!",/, &
              "*******************************************************",/)
  4 format (/,"*******************************************************",/, &
              " THE CONVERGENCE CRITERION FOR THE ORDER OF MAGNITUDE",/, &
              " REDUCTION OF THE ABSOLUTE RESIDUAL ERROR HAS BEEN MET!",/, &
              "*******************************************************",/)
  5 format (/,"*******************************************************",/, &
              " THE SIMULATION HAS REACHED THE REQUESTED FINAL TIME!",/, &
              "*******************************************************",/)
  6 format (/,"******************************************************",/, &
              " VORTEX HAS APPEARED TO COLLAPSE! THE SIMULATION WAS",/, &
              " RUN ANOTHER 25 PERIODS PAST THE POINT OF COLLAPSE!",/, &
              " STOPPING EXECUTION AND DUMPING THE FINAL RESULTS!",/, &
              "******************************************************",/)
  !
end subroutine write_iter_stats_gpu
!
!###############################################################################
!
subroutine write_iter_stats(walltime_expiring,exit_main_loop)
  !
  !.. Use Statements ..
  use geovar,  only : global_pinc_ptr
  use ovar,    only : itcur,d_t,cfl,time,minimum_timesteps,time_ref
  use ovar,    only : convergence_reached,itestcase,Period_Timesteps
  use ovar,    only : convergence_order_max_res,convergence_order_abs_res
  use ovar,    only : uverr1,this_is_final_timestep,iter_out_interval
  use ovar,    only : output_time_averaging,time_scaling_factor,ave_start_time
  !
  use pbs_mod, only : pbs_remaining_walltime
  !
  !.. Formal Arguments ..
  logical(lk), intent(inout) :: walltime_expiring
  logical(lk), intent(  out) :: exit_main_loop
  !
  !.. Local Scalars ..
  integer  :: hrs_left,min_left,sec_left
  integer  :: nloc,kloc
  real(wp) :: xentmax,yentmax
  real(wp) :: entL2
  real(wp) :: wall_time_left,out_time
  real(wp) :: orders_converged,output_L2val
  real(wp) :: rltwo
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_iter_stats"
  !
continue
  !
  exit_main_loop = fals
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Get the density and entropy residual stats from uverr1
  !
  rltwo = uverr1(nec,2,1)
  entL2 = uverr1(nq+1,2,1)
  !
  nloc = nint( uverr1(nec,4,1) )
  kloc = nint( uverr1(nec,5,1) )
  !
  ! Coordinates Location of the max inf norm for entropy
  xentmax = uverr1(nq+1,6,1)
  yentmax = uverr1(nq+1,7,1)
  !
  ! Maximum L2 norm of density for the current simulation
  rl2mx1 = max( rltwo , rl2mx1 )
  !
  if (any(itestcase == Transport_Problems)) then
    !
    output_L2val = entL2
    convergence_reached = fals
    !
    ! Check for a collapse of the vortex, or update the final
    ! countdown if a collapse has already been detected
    !
    if (vortex_has_collapsed) then
      !
      ! Need to subtract from final_countdown the number of time steps
      ! that have been completed since the last call to write_iter_stats
      !
      final_countdown = max(final_countdown - iter_out_interval,0)
      !
    else if (rltwo > zero) then
      !
      if (log10(rltwo) > Density_Maximum_Residual) then
        Density_Maximum_Residual = log10(rltwo)
      end if
      if (Density_Maximum_Residual > -two) then
        vortex_has_collapsed = true
        final_countdown = 26*Period_Timesteps - mod(itcur,Period_Timesteps)
      end if
      !
    end if
    !
  else
    !
    orders_converged = log10(max(rl2mx1,eps17)) - log10(max(rltwo,eps17))
    output_L2val = orders_converged
    !
    ! Check for convergence
    !
    if (itcur > minimum_timesteps) then
      convergence_reached = (-orders_converged <= convergence_order_max_res) &
                            .or. (log10(rltwo) <= convergence_order_abs_res)
      if (itestcase == Taylor_Green_Vortex) then
        convergence_reached = fals
      end if
    end if
    !
  end if
  !
  ! Output the results for the current time step
  !
  call pbs_remaining_walltime(walltime_expiring,wall_time_left)
  !
  if (mypnum == glb_root) then
    !
    out_time = time*time_ref*time_scaling_factor
    if (output_time_averaging) then
      out_time = out_time - ave_start_time*time_ref*time_scaling_factor
    end if
#ifdef PBS_ENV
    !
    sec_left = abs(nint(wall_time_left))
    !
    hrs_left = sec_left/3600
    min_left = (sec_left - hrs_left*3600)/60
    sec_left = sec_left - hrs_left*3600 - min_left*60
    ! Make hrs_left negative if wall_time_left happens to be negative
    hrs_left = sign(hrs_left,int(wall_time_left))
    !
   !write (iout,iter_out_fmt) itcur,d_t*time_ref,time*time_ref, &
    write (iout,iter_out_fmt) itcur,d_t,out_time, &
                              cfl,rltwo,output_L2val,nloc, &
                              global_pinc_ptr(nloc)+kloc, &
                              xentmax,yentmax, &
                              hrs_left,min_left,sec_left
    !
#else
    !
   !write (iout,iter_out_fmt) itcur,d_t*time_ref,time*time_ref, &
    write (iout,iter_out_fmt) itcur,d_t,out_time, &
                              cfl,rltwo,output_L2val,nloc, &
                              global_pinc_ptr(nloc)+kloc, &
                              xentmax,yentmax
    !
    !
#endif
    !
    if (convergence_reached) then
      if (-orders_converged <= convergence_order_max_res) write (iout,3)
      if (log10(rltwo) <= convergence_order_abs_res) write (iout,4)
    end if
    !
    if (walltime_expiring) write (iout,2)
    !
    if (this_is_final_timestep) write (iout,5)
    !
  end if
  !
  ! If a vortex collapse has been detected and final_countdown has reached 0,
  ! activate convergence_reached so that the main_loop can be exited when
  ! we return from this subroutine to the main program.
  !
  if (vortex_has_collapsed .and. final_countdown==0) then
    if (mypnum == glb_root) write (iout,6)
    convergence_reached = true
  end if
  !
  if (convergence_reached) exit_main_loop = true
  if (this_is_final_timestep) exit_main_loop = true
  if (walltime_expiring) exit_main_loop = true
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  2 format (/,"*******************************************************",/, &
              " THE SIMULATION WALL TIME REQUESTED IS EXPIRING SOON!",/, &
              " STOPPING EXECUTION AND DUMPING THE FINAL RESULTS!",/, &
              "*******************************************************",/)
  3 format (/,"*******************************************************",/, &
              " THE CONVERGENCE CRITERION FOR THE ORDERS OF",/, &
              " MAGNITUDE REDUCTION OF THE RESIDUAL ERROR",/, &
              " RELATIVE TO THE MAXIMUM RESIDUAL ERROR HAS BEEN MET!",/, &
              "*******************************************************",/)
  4 format (/,"*******************************************************",/, &
              " THE CONVERGENCE CRITERION FOR THE ORDER OF MAGNITUDE",/, &
              " REDUCTION OF THE ABSOLUTE RESIDUAL ERROR HAS BEEN MET!",/, &
              "*******************************************************",/)
  5 format (/,"*******************************************************",/, &
              " THE SIMULATION HAS REACHED THE REQUESTED FINAL TIME!",/, &
              "*******************************************************",/)
  6 format (/,"******************************************************",/, &
              " VORTEX HAS APPEARED TO COLLAPSE! THE SIMULATION WAS",/, &
              " RUN ANOTHER 25 PERIODS PAST THE POINT OF COLLAPSE!",/, &
              " STOPPING EXECUTION AND DUMPING THE FINAL RESULTS!",/, &
              "******************************************************",/)
  !
end subroutine write_iter_stats
!
!###############################################################################
!
subroutine init_generic_mod
  !
  !.. Use Statements ..
  use geovar, only : nr,n_global_cell,n_global_solpts
  use ovar,   only : itcur,num_timesteps,itestcase
  use ovar,   only : uverr,uverr1,output_time_averaging
  !
  use pbs_mod, only : pbs_remaining_walltime
  !
  !.. Local Scalars ..
  integer     :: ierr
  integer     :: final_iter_digits
  integer     :: max_cell_digits
  integer     :: max_solpt_digits
  integer     :: hrs_left,min_left,sec_left
  logical(lk) :: walltime_expiring
  real(wp)    :: wall_time_left
  character(len=100) :: char_num
  character(len=1000) :: iter_header
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_generic_mod"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate the error array
  !
  if (.not. allocated(uverr)) then
    allocate ( uverr(1:nq+1,1:5+nr,1:2) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"uverr",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
  end if
  !$ACC UPDATE DEVICE(uverr)
  if (.not. allocated(uverr1)) then
    allocate ( uverr1(1:nq+1,1:5+nr,1:2) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"uverr1",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
  end if
  !
  write (char_num,1) itcur + num_timesteps
  final_iter_digits = max(4,len_trim(adjustl(char_num)))
  !
  write (char_num,1) n_global_cell
  max_cell_digits = max(4,len_trim(adjustl(char_num)))
  !
  write (char_num,1) n_global_solpts
  max_solpt_digits = max(5,len_trim(adjustl(char_num)))
  !
 !if (output_time_averaging) then
 !  write (iter_out_fmt,101) final_iter_digits
 !else
    write (iter_out_fmt,100) final_iter_digits,max_cell_digits,max_solpt_digits
 !end if
  !
  ! Add on the remaining simulation time to the output if
  ! running in the PBS environment
  !
#ifdef PBS_ENV
  iter_out_fmt = iter_out_fmt(1:len_trim(iter_out_fmt)) // ',1x,i0,2(":",i2.2)'
#endif
  !
  ! Put the last ')' onto iter_out_fmt
  !
  iter_out_fmt = iter_out_fmt(1:len_trim(iter_out_fmt)) // ")"
  !
  ! Now get the format for the header
  !
  if (any(itestcase == Transport_Problems)) then
    write (iter_header,200) final_iter_digits-3, &
                            max_cell_digits-2, &
                            max_solpt_digits-3
  else if (itestcase == Laminar_BndLyr) then
    write (iter_header,201) final_iter_digits-3, &
                            max_cell_digits-2, &
                            max_solpt_digits-3
  else if (itestcase == Taylor_Green_Vortex) then
   !write (iter_header,203) final_iter_digits-3, &
    write (iter_header,202) final_iter_digits-3, &
                            max_cell_digits-2, &
                            max_solpt_digits-3
 !else if (output_time_averaging) then
 !  write (iter_header,204) final_iter_digits-3
  else
    write (iter_header,202) final_iter_digits-3, &
                            max_cell_digits-2, &
                            max_solpt_digits-3
  end if
  !
#ifdef PBS_ENV
  iter_header = iter_header(1:len_trim(iter_header)) // ",4x,'Time Left'"
  !
  ! Get the remaining wall time before the solver begins
  !
  call pbs_remaining_walltime(walltime_expiring,wall_time_left)
  !
  sec_left = abs(nint(wall_time_left))
  !
  hrs_left = sec_left/3600
  min_left = (sec_left - hrs_left*3600)/60
  sec_left = sec_left - hrs_left*3600 - min_left*60
  ! Make hrs_left negative if wall_time_left happens to be negative
  hrs_left = sign(hrs_left,int(wall_time_left))
  !
  if (mypnum == 0) write (iout,2) hrs_left,min_left,sec_left
  !
#endif
  !
  ! Put the last ')' onto iter_header
  !
  iter_header = iter_header(1:len_trim(iter_header)) // ")"
  !
  ! Write the header for the terminal output
  !
  if (mypnum == 0) write (iout,iter_header)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (i0)
  2 format (/,"  After Pre-Processor and Initializations =>", &
              "  Wall Time Remaining = ",i0,2(":",i2.2))
  100 format ("(1x,i",i0,",5es12.4,2x,i",i0,",2x,i",i0,",2es9.1")
  101 format ("(1x,i",i0,",8es12.4")
  200 format ("(/,",i0,"x,'iter',5x,'d_t',8x,'Time',6x,'L2-Density',2x,", &
              "'L2-Entropy',",i0,"x,'Cell',",i0,"x,'Solpt',5x,'x',8x,'y'")
  201 format ("(/,",i0,"x,'iter',5x,'d_t',8x,'Time',6x,'L2-X-Mmntm',2x,", &
              "'OrdersDrop',2x,'Max-DensL2',",i0,"x,'Cell',",i0, &
              "x,'Solpt',5x,'x',8x,'y'")
 !202 format ("(/,",i0,"x,'iter',5x,'d_t',8x,'Time',10x,'CFL',5x,", &
  202 format ("(/,",i0,"x,'iter',3x,'d_t* (ND)',3x,'Time (s)',7x,'CFL',5x,", &
              "'L2-Density',2x,'Max-DensL2',",i0,"x,'Cell',",i0, &
              "x,'Solpt',5x,'x',8x,'y'")
 !202 format ("(/,",i0,"x,'iter',5x,'d_t',8x,'Time',6x,'L2-Density',2x,", &
 !            "'OrdersDrop',2x,'Max-DensL2',",i0,"x,'Cell',",i0, &
 !            "x,'Solpt',5x,'x',8x,'y'")
  203 format ("(/,",i0,"x,'iter',5x,'dt*',8x,'Time*',5x,'L2-Density',2x,", &
              "'OrdersDrop',2x,'Max-DensL2',",i0,"x,'Cell',",i0, &
              "x,'Solpt',5x,'x',8x,'y'")
  204 format ("(/,",i0,"x,'iter',5x,'dt*',8x,'Time*',5x,'L2-Density',4x,", &
              "'L2-Vx',7x,'L2-Vy',7x,'L2-Vz',7x,'L2-Pres'")
  !
end subroutine init_generic_mod
!
!###############################################################################
!
subroutine getmaxent
  !
  !.. Use Statements ..
  use ovar, only : uverr
  !
  !.. Local Scalars ..
  integer  :: nentmax,kentmax
  real(wp) :: xentmax,yentmax,zentmax
  real(wp) :: entmax,entL1,entL2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "getmaxent"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  entL1   =       uverr(nq+1,1,1)
  entL2   =       uverr(nq+1,2,1)
  entmax  =       uverr(nq+1,3,1)
  nentmax = nint( uverr(nq+1,4,1) )
  kentmax = nint( uverr(nq+1,5,1) )
  xentmax =       uverr(nq+1,6,1)
  if (size(uverr,dim=2) >= 7) &
  yentmax =       uverr(nq+1,7,1)
  if (size(uverr,dim=2) == 8) &
  zentmax =       uverr(nq+1,8,1)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine getmaxent
!
!###############################################################################
!
subroutine collect_error_arrays(this_cell,unew,uold,wtdjac,xyz_err,vnew,vold)
  !
  !.. Use Statements ..
  use order_mod,         only : maxSP
  use geovar,            only : nr,cell,xyz
  use ovar,              only : mms_opt,itestcase,bc_in
  use ovar,              only : itcur,itrst,read_time_ave_restart
  use ovar,              only : prev_ave_time,saved_ave_time,ave_start_time
  use ovar,              only : output_time_averaging
  use flowvar,           only : usp,uoldsp,uavesp
  use interpolation_mod, only : outerp,outerp_mv_mult
  use quadrature_mod,    only : std_elem
  use metrics_mod,       only : jac
  use mms_mod,           only : mms_solution_at_xyz
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: this_cell
  real(wp), dimension(:,:), intent(inout) :: unew
  real(wp), dimension(:,:), intent(inout) :: uold
  real(wp), dimension(:,:), intent(inout) :: xyz_err
  real(wp), dimension(:),   intent(inout) :: wtdjac
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(inout) :: vnew(:,:)
  real(wp), optional, intent(inout) :: vold(:,:)
  !
  !.. Local Scalars ..
  integer :: k,l,m,n,n0,n1,n2,np,nep
  integer :: this_geom,this_order
  real(wp) :: a,b
  logical(lk) :: use_error_matrix
  !
  !.. Local Arrays ..
  real(wp) :: xyz_offset(1:nr,1:maxSP)
  real(wp) :: uprev(1:nq,1:maxSP)
  real(wp) :: vprev(1:nq,1:maxSP)
  real(wp) :: ucurr(1:nq,1:maxSP)
  real(wp) :: vcurr(1:nq,1:maxSP)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "collect_error_arrays"
  !
continue
  !
 !call debug_timer(entering_procedure,pname)
  !
  this_geom = cell(this_cell)%geom
  this_order = cell(this_cell)%order
  !
  n1 = cell(this_cell)%beg_sp
  n2 = cell(this_cell)%end_sp
  np = n2-n1+1
  n0 = n1-1
  !
  ! Determine if we will be using the error matrix or not
  !
  use_error_matrix = fals
  !
  if (allocated(outerp(this_geom,this_order)%error)) then
    if (allocated(outerp(this_geom,this_order)%error%mat)) then
      use_error_matrix = true
    end if
  end if
  !
  ! Get the number of error points
  !
  if (use_error_matrix) then
    nep = size(outerp(this_geom,this_order)%error%mat,dim=1)
  else
    nep = np
  end if
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Get the conserved variables of the current solution at the error points
  !
  if (output_time_averaging) then
    do n = 1,np
      do m = 1,nq
        ucurr(m,n) = uavesp(n0+n,m)
      end do
    end do
   !ucurr(1:nq,1:np) = uavesp(1:nq,n1:n2)
  else
    do n = 1,np
      do m = 1,nq
        ucurr(m,n) = usp(n0+n,m)
      end do
    end do
   !ucurr(1:nq,1:np) = usp(1:nq,n1:n2)
  end if
  !
  if (use_error_matrix) then
    !
    do m = 1,nq
      call outerp_mv_mult( outerp(this_geom,this_order)%error , &
                           ucurr(m,1:np) , unew(m,1:nep) )
    end do
    !
  else
    !
    unew(1:nq,1:nep) = ucurr(1:nq,1:np)
    !
  end if
  !
  ! Get the primitive variables of the current solution at the error points
  !
  if (present(vnew)) then
    !
    if (output_time_averaging) then
      do n = 1,np
        do m = 1,nq
          vcurr(m,n) = uavesp(n0+n,nq+m)
        end do
      end do
     !vcurr(1:nq,1:np) = uavesp(nq+1:2*nq,n1:n2)
    else
      do n = 1,np
        vcurr(1:nq,n) = usp2v_sp( usp(n0+n,1:nq) , swap_t_for_rho=true )
      end do
     !vcurr(1:nq,1:np) = usp2v( usp(1:nq,n1:n2) , swap_t_for_rho=true )
    end if
    !
    if (use_error_matrix) then
      !
      do m = 1,nq
        call outerp_mv_mult( outerp(this_geom,this_order)%error , &
                             vcurr(m,1:np) , vnew(m,1:nep) )
      end do
      !
    else
      !
      vnew(1:nq,1:nep) = vcurr(1:nq,1:np)
      !
    end if
    !
  end if
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Get the value of the metric jacobians multiplied
  ! by the quadrature weights at the error points
  !
  if (use_error_matrix) then
    !
    call outerp_mv_mult( outerp(this_geom,this_order)%error , &
                         jac(n1:n2) , wtdjac(1:nep) )
    !
    wtdjac(1:nep) = wtdjac(1:nep) * outerp(this_geom,this_order)%err_wts
    !
  else
    !
    wtdjac(1:nep) = jac(n1:n2) * std_elem(this_geom,this_order)%wts
    !
  end if
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Get the coordinates of the error points
  !
  if (use_error_matrix) then
    !
    do l = 1,nr
      call outerp_mv_mult( outerp(this_geom,this_order)%error , &
                           xyz(l,n1:n2) , xyz_err(l,1:nep) )
    end do
    !
  else
    !
    xyz_err(1:nr,1:nep) = xyz(1:nr,n1:n2)
    !
  end if
  !
  ! #######
  ! ## 4 ##
  ! #######
  !
  ! Get either
  !   (a) the exact analytical solution or
  !   (b) solution from the previous time step
  ! at the error points.
  !
  if (mms_opt /= 0) then
    !
    ! Evaluate the exact analytical MMS solution at the error points
    !
    do k = 1,nep
      uold(1:nq,k) = mms_solution_at_xyz( xyz_err(1:nr,k), &
                                          return_nondim=true, &
                                          return_cv=true )
    end do
    !
    if (present(vold)) then
      vold(1:nq,1:nep) = usp2v( uold(1:nq,1:nep) , swap_t_for_rho=true )
    end if
    !
  else if (any(itestcase == Transport_Problems)) then
    !
    ! Evaluate the exact analytical solution to the vortex transport problem
    !
    xyz_offset(1:nr,1:nep) = vortex_offset( xyz_err(1:nr,1:nep) )
    !
    do k = 1,nep
      uold(1:nq,k) = vortex_solution( xyz_offset(1:nr,k) )
    end do
    !
    if (present(vold)) then
      vold(1:nq,1:nep) = usp2v( uold(1:nq,1:nep) , swap_t_for_rho=true )
    end if
    !
  else if (itestcase == Freestream_Preservation) then
    !
    do k = 1,nep
      uold(1:nq,k) = bc_in(0)%cv(1:nq)
    end do
    !
    if (present(vold)) then
      do k = 1,nep
        vold(1:nq,k) = bc_in(0)%pv(1:nq)
      end do
    end if
    !
  else
    !
    if (output_time_averaging) then
      ! Compute uavesp at the previous averaging time
      if (any(saved_ave_time == [zero,ave_start_time]) .or. &
          (.not. read_time_ave_restart .and. itcur-itrst < 100)) then
        do n = 1,np
          do m = 1,nq
            uprev(m,n) = uoldsp(n0+n,m)
          end do
        end do
       !uprev(1:nq,1:np) = uoldsp(1:nq,n1:n2)
      else
        a = (prev_ave_time-saved_ave_time) / (saved_ave_time-ave_start_time)
        b = one - a
        do n = 1,np
          do m = 1,nq
            uprev(m,n) = b * uavesp(n0+n,m) - a * usp(n0+n,m)
          end do
        end do
       !uprev(1:nq,1:np) = b * uavesp(1:nq,n1:n2) - a * usp(1:nq,n1:n2)
       !if (any(uprev(nec,1:np) <= zero)) then
       !  uprev(1:nq,1:np) = uoldsp(1:nq,n1:n2)
       !else if (any(pressure_cv(uprev(1:nq,1:np)) <= zero)) then
       !  uprev(1:nq,1:np) = uoldsp(1:nq,n1:n2)
       !end if
      end if
    else
      do n = 1,np
        do m = 1,nq
          uprev(m,n) = uoldsp(n0+n,m)
        end do
      end do
     !uprev(1:nq,1:np) = uoldsp(1:nq,n1:n2)
    end if
    !
    if (use_error_matrix) then
      !
      ! Get the solution from the previous time step at the error points
      !
      do m = 1,nq
        call outerp_mv_mult( outerp(this_geom,this_order)%error , &
                             uprev(m,1:np) , uold(m,1:nep) )
      end do
      !
    else
      !
      ! Use the solution from the previous time step
      !
      uold(1:nq,1:nep) = uprev(1:nq,1:np)
      !
    end if
    !
    if (present(vold)) then
      !
      if (output_time_averaging) then
        ! Compute uavesp at the previous averaging time
        if (any(saved_ave_time == [zero,ave_start_time]) .or. &
            (.not. read_time_ave_restart .and. itcur-itrst < 100)) then
          do n = 1,np
            vprev(1:nq,n) = usp2v_sp( uoldsp(n0+n,1:nq) , swap_t_for_rho=true )
          end do
         !vprev(1:nq,1:np) = usp2v( uoldsp(1:nq,n1:n2) , swap_t_for_rho=true )
        else
          ! a and b were already computed above, no need to do it again
          do n = 1,np
            vprev(1:nq,n) = b * uavesp(n0+n,nq+1:2*nq) - &
                            a * usp2v_sp( usp(n0+n,1:nq) , swap_t_for_rho=true )
          end do
         !vprev(1:nq,1:np) = b * uavesp(nq+1:2*nq,n1:n2) - &
         !                   a * usp2v( usp(1:nq,n1:n2) , swap_t_for_rho=true )
        end if
      else
        do n = 1,np
          vprev(1:nq,n) = usp2v_sp( uoldsp(n0+n,1:nq) , swap_t_for_rho=true )
        end do
       !vprev(1:nq,1:np) = usp2v( uoldsp(1:nq,n1:n2) , swap_t_for_rho=true )
      end if
      !
      if (use_error_matrix) then
        !
        ! Get the solution from the previous time step at the error points
        !
        do m = 1,nq
          call outerp_mv_mult( outerp(this_geom,this_order)%error , &
                               vprev(m,1:np) , vold(m,1:nep) )
        end do
        !
      else
        !
        ! Use the solution from the previous time step
        !
        vold(1:nq,1:nep) = vprev(1:nq,1:np)
        !
      end if
      !
    end if
    !
  end if
  !
 !call debug_timer(leaving_procedure,pname)
  !
end subroutine collect_error_arrays
!
!###############################################################################
!
subroutine compute_TaylorGreen_KEDR_gpu(iter,dtmin)
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use geovar, only : nr,ncell,n_solpts
  use ovar, only : muref_nd,total_volume
  use ovar, only : time,time_ref
  use ovar, only : VortexGrid_DOF
  use ovar, only : iter_out_interval
  use flowvar, only : usp,dusp
  use quadrature_mod, only : std_elem
  use metrics_mod, only : wtdjac
  use gpu_mod, only : grad_cv_to_grad_pv_gpu
  use gpu_mod, only : tusp2v_gpu
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: iter
  real(wp), intent(in) :: dtmin
  !
  !.. Local Scalars ..
  integer :: k,n,nc,np,n0,tgv_idx
  integer :: ierr,nmin,nmax,nbeg,nend
  integer :: this_cell,this_geom,this_order
  !
  logical(lk) :: first_output,final_output
  logical(lk) :: this_is_output_iter
  !
  real(wp) :: kem,kec,kep,dtm,dtp,vol_inv
  real(wp) :: solution_time,muref_star
  real(wp) :: sum_kedr_1,sum_kedr_2
  real(wp) :: sum_kedr_3,sum_kedr_4
  real(wp) :: rho,rho_inv,pres,q2,vort,div_vel
  real(wp) :: a1,a2,a3,a4,a5,a6,strain
  real(wp) :: velx,dudx,dvdx,dwdx,vortx
  real(wp) :: vely,dudy,dvdy,dwdy,vorty
  real(wp) :: velz,dudz,dvdz,dwdz,vortz
  !
  !.. Local Saved Scalars ..
  logical(lk), save :: completed_first_output = fals
  logical(lk), save :: needs_initialization = true
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: tvsp(:,:)
  real(wp), allocatable :: tdpvdx(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_TaylorGreen_KEDR_gpu"
  character(len=*), parameter :: fname = "TaylorGreen_KEDR.dat"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nmin = 1
  nmax = iter_out_interval
  !
  ! Get the current index value for the kedr array
  ! NOTE: tgv_idx will be zero if iter is 0, otherwise it will
  !       take a value between 1 and iter_out_interval
  !
  tgv_idx = mod( abs(iter)-1 , nmax ) + 1
  !
  ! Create some logical variables to help readability
  !
  final_output = (iter < 0)
  first_output = (iter == nmax) .or. &
                 (final_output .and. .not. completed_first_output)
  this_is_output_iter = ((tgv_idx == nmax .or. final_output) .and. (iter /= 0))
  !
  !####################################################
  !####################################################
  !####################################################
  !##########                                ##########
  !##########     INITIALIZATION SECTION     ##########
  !##########                                ##########
  !####################################################
  !####################################################
  !####################################################
  !
  ! If this is before the first iteration, initialize everything needed
  ! for this subroutine
  !
  if (needs_initialization) then
    !
    ! Reset tgv_idx to 0 for the first time through this subroutine
    !
    tgv_idx = 0
    !
    ! NOTE: tgv_kedr(7,-1) will store the solution time for the last
    !       output of the TGV results so initializing tgv_kedr to zero
    !       is exactly what we need to start
    !
    allocate ( tgv_kedr(1:7,-1:nmax) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"tgv_kedr",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !$ACC UPDATE DEVICE(tgv_kedr)
    !
    allocate ( sum_kedr(1:4,0:nmax) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"sum_kedr",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !$ACC UPDATE DEVICE(sum_kedr)
    !
    ! Initialize the array to adjust the nondimensionalization of the
    ! integrated results
    ! NOTE: The comments above each value indicate the reason for
    !       the number of time_refs
    !
    ! Adjust the nondimensionalizatin of muref
    !
    muref_star = muref_nd / time_ref
    !
    ! two for square of velocity derivatives (from vorticity**2)
    tgv_dimen(2) = muref_star / time_ref / time_ref
    ! two for square of velocity derivatives (from strain**2)
    tgv_dimen(3) = two * muref_star / time_ref / time_ref
    ! two for pressure, one for velocity derivative
    tgv_dimen(4) = -one / time_ref / time_ref / time_ref
    ! two for velocity squared in kinetic energy term
    tgv_dimen(5) = half / time_ref / time_ref
    ! two for square of velocity derivatives (from vorticity**2)
    tgv_dimen(6) = half / time_ref / time_ref
    !
    ! Set tgv_kedr(7,-1) to time the first time through
    ! time will have the value zero for a new simulation or the
    ! restart time if restarting a simulation
    !
    tgv_kedr(7,-1) = time * time_ref
    !$ACC UPDATE DEVICE(tgv_kedr,sum_kedr,tgv_dimen)
    !
    ! Open the output file for the root process
    !
    if (mypnum == 0) then
      !
      nc = nint( real(VortexGrid_DOF,kind=wp) / real(n_order+1,kind=wp) )
      !
      open (newunit=iotg,file=fname,status="replace",action="write", &
            iostat=ierr,iomsg=error_message)
      call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
      !
      write (iotg,1) n_order,nc*(n_order+1),n_order, &
                     VortexGrid_DOF,nc*(n_order+1)
      !
    end if
    !
    ! Unset needs_initialization now that everything has been initialized
    !
    needs_initialization = fals
    !
  end if
  !
  allocate ( tvsp(1:n_solpts,1:nq) )
  allocate ( tdpvdx(1:n_solpts,1:nr,1:nq) )
  !
  !$ACC DATA CREATE(tvsp,tdpvdx)
  !
  call grad_cv_to_grad_pv_gpu(1,n_solpts,usp,dusp,tdpvdx)
  !
  call tusp2v_gpu(1,n_solpts,usp,tvsp)
  !
  !######################################################
  !######################################################
  !######################################################
  !##########                                  ##########
  !##########     MAIN COMPUTATION SECTION     ##########
  !##########                                  ##########
  !######################################################
  !######################################################
  !######################################################
  !
  ! Loop over all the solution points to compute the
  ! integrated quantities local to each processor
  !
  sum_kedr_1 = zero
  sum_kedr_2 = zero
  sum_kedr_3 = zero
  sum_kedr_4 = zero
  !
  !$ACC PARALLEL LOOP REDUCTION(+:sum_kedr_1,sum_kedr_2) &
  !$ACC&              PRESENT(wtdjac,tvsp,tdpvdx)
  do n = 1,n_solpts
    !
    rho = tvsp(n,nec)
    !
    dudx = tdpvdx(n,1,nmx)
    dudy = tdpvdx(n,2,nmx)
    dudz = tdpvdx(n,3,nmx)
    !
    dvdx = tdpvdx(n,1,nmy)
    dvdy = tdpvdx(n,2,nmy)
    dvdz = tdpvdx(n,3,nmy)
    !
    dwdx = tdpvdx(n,1,nmz)
    dwdy = tdpvdx(n,2,nmz)
    dwdz = tdpvdx(n,3,nmz)
    !
    vortx = dwdy - dvdz
    vorty = dudz - dwdx
    vortz = dvdx - dudy
    vort  = vortx*vortx + vorty*vorty + vortz*vortz
    !
    a1 = dudy + dvdx
    a2 = dudz + dwdx
    a3 = dvdz + dwdy
    !
    a4 = dudx * (dudx - dvdy)
    a5 = dvdy * (dvdy - dwdz)
    a6 = dwdz * (dwdz - dudx)
    !
    strain = half*(a1*a1 + a2*a2 + a3*a3) + two3*(a4 + a5 + a6)
    !
    sum_kedr_1 = sum_kedr_1 + wtdjac(n)*rho*vort
    sum_kedr_2 = sum_kedr_2 + wtdjac(n)*strain
    !
  end do
  !
  !$ACC PARALLEL LOOP REDUCTION(+:sum_kedr_3,sum_kedr_4) &
  !$ACC&              PRESENT(wtdjac,tvsp,tdpvdx)
  do n = 1,n_solpts
    !
    rho  = tvsp(n,nec)
    velx = tvsp(n,nmx)
    vely = tvsp(n,nmy)
    velz = tvsp(n,nmz)
    pres = tvsp(n,nee)
    !
    q2 = velx*velx + vely*vely + velz*velz
    !
    dudx = tdpvdx(n,1,nmx)
    dvdy = tdpvdx(n,2,nmy)
    dwdz = tdpvdx(n,3,nmz)
    !
    div_vel = dudx + dvdy + dwdz
    !
    sum_kedr_3 = sum_kedr_3 + wtdjac(n)*pres*div_vel
    sum_kedr_4 = sum_kedr_4 + wtdjac(n)*rho*q2
    !
  end do
  !
  sum_kedr(1,tgv_idx) = sum_kedr_1
  sum_kedr(2,tgv_idx) = sum_kedr_2
  sum_kedr(3,tgv_idx) = sum_kedr_3
  sum_kedr(4,tgv_idx) = sum_kedr_4
  !$ACC UPDATE DEVICE(sum_kedr)
  !
  !
  !$ACC END DATA
  deallocate ( tdpvdx )
  deallocate ( tvsp )
  !
  ! Store the current nondimensional time step
  !
  tgv_kedr(7,tgv_idx) = dtmin * time_ref
  !$ACC UPDATE DEVICE(tgv_kedr)
  !
  !#####################################################
  !#####################################################
  !#####################################################
  !##########                                 ##########
  !##########     INTERVAL OUTPUT SECTION     ##########
  !##########                                 ##########
  !#####################################################
  !#####################################################
  !#####################################################
  !
  ! If this is an iteration to output the results, collect the KEDR values
  ! since the last output and dump the results to the Taylor-Green results file.
  !
  if (this_is_output_iter) then
    !
    ! Sum up the KEDR values across all the processors
    !
    if (ncpu > 1) then
      !
      call mpi_allreduce(MPI_IN_PLACE,sum_kedr,size(sum_kedr,kind=int_mpi), &
                         mpi_flttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
      !
      !$ACC UPDATE DEVICE(sum_kedr)
    end if
    !
    ! If this this the first output interval, change nmin to 0 so that the
    ! KEDR values for the initial conditions are also modified
    !
    if (first_output) nmin = 0
    !
    ! Finish computing the different KEDR
    !
    vol_inv = one/total_volume
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) PRESENT(tgv_kedr,sum_kedr)
    do n = nmin,nmax
      tgv_kedr(2,n) = sum_kedr(1,n) * tgv_dimen(2) * vol_inv
      tgv_kedr(3,n) = sum_kedr(2,n) * tgv_dimen(3) * vol_inv
      tgv_kedr(4,n) = sum_kedr(3,n) * tgv_dimen(4) * vol_inv
      tgv_kedr(5,n) = sum_kedr(4,n) * tgv_dimen(5) * vol_inv
      tgv_kedr(6,n) = sum_kedr(1,n) * tgv_dimen(6) * vol_inv
    end do
    !
    ! Compute the central difference of the KEDR to get its rate of change
    ! over time
    ! NOTE: The KEDR found by differencing the integrated kinetic energy is
    !       actually given as KEDR(Ek) = -d/dt (Ek) so the differencing
    !       equation below has been multiplied through by a -1.
    !
    nbeg = merge(     1     ,   0    , first_output )
    nend = merge( tgv_idx-1 , nmax-1 , final_output )
    !
    !$ACC PARALLEL LOOP VECTOR_LENGTH(32) PRESENT(tgv_kedr)
    do n = nbeg,nend
      !
      kem = tgv_kedr(5,n-1)
      kec = tgv_kedr(5,n)
      kep = tgv_kedr(5,n+1)
      !
      dtm = tgv_kedr(7,n)
      dtp = tgv_kedr(7,n+1)
      !
      tgv_kedr(1,n) =  ( (dtm/dtp)*(kec-kep) + (dtp/dtm)*(kem-kec) ) / (dtm+dtp)
      !
    end do
    !$ACC UPDATE SELF(tgv_kedr)
    !
    ! Do a one sided difference on the initial conditions if this is the
    ! first output interval
    !
    if (first_output) then
      tgv_kedr(1,0) = (tgv_kedr(5,0) - tgv_kedr(5,1)) / tgv_kedr(7,1)
    end if
    !
    ! If this is the last output interval, do a one sided difference
    ! on the final solution
    !
    if (final_output) then
      n = tgv_idx
      tgv_kedr(1,n) = (tgv_kedr(5,n-1) - tgv_kedr(5,n)) / tgv_kedr(7,n)
    end if
    !
    ! Have the root processor output the Taylor-Green KEDR solutions
    !
    if (mypnum == 0) then
      !
      solution_time = tgv_kedr(7,-1)
      do n = 0,nend
        solution_time = solution_time + tgv_kedr(7,n)
        write (iotg,2) solution_time, (tgv_kedr(k,n), k=1,size(tgv_kedr,dim=1))
      end do
      !
      ! Perform a few final tasks if this is the final output
      !
      if (final_output) then
        !
        ! Output the last value of tgv_kedr
        !
        n = tgv_idx
        solution_time = solution_time + tgv_kedr(7,n)
        write (iotg,2) solution_time, (tgv_kedr(k,n), k=1,size(tgv_kedr,dim=1))
        !
        ! Close the Taylor-Green KEDR solution file
        !
        close (iotg,iostat=ierr,iomsg=error_message)
        call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
        !
      end if
      !
    else
      !
      solution_time = sum( tgv_kedr(7,-1:nend) )
      !
    end if
    !
    ! Peform some maintence on the tgv_kedr array before continuing
    !
    if (final_output) then
      !
      ! Deallocate tgv_kedr since it is no longer needed
      !
      deallocate ( tgv_kedr , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"tgv_kedr",2,__LINE__,__FILE__,ierr, &
                       error_message,skip_alloc_pause)
      !
      deallocate ( sum_kedr , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"sum_kedr",2,__LINE__,__FILE__,ierr, &
                       error_message,skip_alloc_pause)
      !
    else
      !
      ! Copy tgv_kedr(:,nmax-1:nmax) to tgv_kedr(:,-1:0) so they can
      ! be used for the next output interval. Set tgv_kedr(7,-1) to
      ! the last solution time that was output. Initialize all the other
      ! values to zero.
      !
      tgv_kedr(:,-1) = [tgv_kedr(1:6,nmax-1),solution_time]
      tgv_kedr(:,0)  = tgv_kedr(:,nmax)
      tgv_kedr(:,1:nmax) = zero
      sum_kedr = zero
      !$ACC UPDATE DEVICE(tgv_kedr,sum_kedr)
      !
    end if
    !
    ! Set completed_first_output to true since
    ! we have now completed at least one output
    !
    completed_first_output = true
    !
  end if
  !
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ('VARIABLES = "t*"',/, &
            '"<greek>e</greek>(E<sub>k</sub>)"',/, &
            '"<greek>e</greek>(<greek>z</greek>)"',/, &
            '"<greek>e</greek><sub>1</sub> (Deviatoric Strain Rate)"',/, &
            '"<greek>e</greek><sub>3</sub> (Pressure Dialation)"',/, &
            '"E<sub>k</sub>"',/, &
            '"<greek>z</greek>"',/, &
            '"<greek>D</greek>t*"',/, &
            'ZONE T="P',i0,' - ',i0,'<sup>3</sup>"',/, &
            'AUXDATA order="',i0,'"',/, &
            'AUXDATA RDOF="',i0,'<sup>3</sup>"',/, &
            'AUXDATA DOF="',i0,'<sup>3</sup>"')
  2 format (100(es17.9))
 !2 format (*(es17.9))
  !
end subroutine compute_TaylorGreen_KEDR_gpu
!
!###############################################################################
!
subroutine compute_TaylorGreen_KEDR(iter,dtmin)
  !
  !.. Use Statements ..
  use order_mod, only : maxpts,n_order
  use geovar, only : nr,ncell,cell
  use ovar, only : muref_nd,total_volume
  use ovar, only : time,time_ref
  use ovar, only : VortexGrid_DOF
  use ovar, only : iter_out_interval
  use flowvar, only : usp,dusp
  use interpolation_mod, only : outerp,outerp_mv_mult
  use quadrature_mod, only : std_elem
  use metrics_mod, only : jac
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: iter
  real(wp), intent(in) :: dtmin
  !
  !.. Local Scalars ..
  integer :: k,l,m,n,nc,nep,n1,n2,tgv_idx
  integer :: ierr,nke,ndt,nmin,nmax,nbeg,nend
  integer :: this_cell,this_geom,this_order
  !
  logical(lk) :: first_output,final_output
  logical(lk) :: this_is_output_iter
  !
  real(wp) :: kem,kec,kep,dtm,dtp
  real(wp) :: solution_time
  real(wp) :: div_vel,muref_star
  !
  !.. Local Saved Scalars ..
  logical(lk), save :: completed_first_output = fals
  logical(lk), save :: needs_initialization = true
  !
  !.. Local Arrays ..
  real(wp) :: vort(1:nr)
  real(wp) :: wtdjac(1:maxpts)
  real(wp) :: strain(1:nr,1:nr)
  real(wp) :: vsp(1:nq,1:maxpts)
  real(wp) :: dcvdx(1:nr,1:nq,1:maxpts)
  real(wp) :: dpvdx(1:nr,1:nq,1:maxpts)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_TaylorGreen_KEDR"
  character(len=*), parameter :: fname = "TaylorGreen_KEDR.dat"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nke = 5
  ndt = 7
  !
  nmin = 1
  nmax = iter_out_interval
  !
  ! Get the current index value for the kedr array
  ! NOTE: tgv_idx will be zero if iter is 0, otherwise it will
  !       take a value between 1 and iter_out_interval
  !
  tgv_idx = mod( abs(iter)-1 , nmax ) + 1
  !
  ! Create some logical variables to help readability
  !
  final_output = (iter < 0)
  first_output = (iter == nmax) .or. &
                 (final_output .and. .not. completed_first_output)
  this_is_output_iter = ((tgv_idx == nmax .or. final_output) .and. (iter /= 0))
  !
  !####################################################
  !####################################################
  !####################################################
  !##########                                ##########
  !##########     INITIALIZATION SECTION     ##########
  !##########                                ##########
  !####################################################
  !####################################################
  !####################################################
  !
  ! If this is before the first iteration, initialize everything needed
  ! for this subroutine
  !
  if (needs_initialization) then
    !
    ! Reset tgv_idx to 0 for the first time through this subroutine
    !
    tgv_idx = 0
    !
    ! NOTE: tgv_kedr(ndt,-1) will store the solution time for the last
    !       output of the TGV results so initializing tgv_kedr to zero
    !       is exactly what we need to start
    !
    allocate ( tgv_kedr1(1:7,-1:nmax) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"tgv_kedr1",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !
    allocate ( sum_kedr1(1:4,0:nmax) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"sum_kedr1",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !
    ! Initialize the array to adjust the nondimensionalization of the
    ! integrated results
    ! NOTE: The comments above each value indicate the reason for
    !       the number of time_refs
    !
    ! Adjust the nondimensionalizatin of muref
    !
    muref_star = muref_nd / time_ref
    !
    ! two for square of velocity derivatives (from vorticity**2)
    tgv_dimen1(2) = muref_star / time_ref / time_ref
    ! two for square of velocity derivatives (from strain**2)
    tgv_dimen1(3) = two * muref_star / time_ref / time_ref
    ! two for pressure, one for velocity derivative
    tgv_dimen1(4) = -one / time_ref / time_ref / time_ref
    ! two for velocity squared in kinetic energy term
    tgv_dimen1(5) = half / time_ref / time_ref
    ! two for square of velocity derivatives (from vorticity**2)
    tgv_dimen1(6) = half / time_ref / time_ref
    !
    ! Set tgv_kedr(ndt,-1) to time the first time through
    ! time will have the value zero for a new simulation or the
    ! restart time if restarting a simulation
    !
    tgv_kedr1(ndt,-1) = time * time_ref
    !
    ! Open the output file for the root process
    !
    if (mypnum == 0) then
      !
      nc = nint( real(VortexGrid_DOF,kind=wp) / real(n_order+1,kind=wp) )
      !
      open (newunit=iotg1,file=fname,status="replace",action="write", &
            iostat=ierr,iomsg=error_message)
      call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
      !
      write (iotg1,1) n_order,nc*(n_order+1),n_order, &
                     VortexGrid_DOF,nc*(n_order+1)
      !
    end if
    !
    ! Unset needs_initialization now that everything has been initialized
    !
    needs_initialization = fals
    !
  end if
  !
  !######################################################
  !######################################################
  !######################################################
  !##########                                  ##########
  !##########     MAIN COMPUTATION SECTION     ##########
  !##########                                  ##########
  !######################################################
  !######################################################
  !######################################################
  !
  ! Store the current nondimensional time step
  !
  tgv_kedr1(ndt,tgv_idx) = dtmin * time_ref
  !
  ! Loop through the cells to compute the integrated
  ! quantities local to each processor
  !
  do this_cell = 1,ncell
    !
    n1 = cell(this_cell)%beg_sp
    n2 = cell(this_cell)%end_sp
    !
    this_geom = cell(this_cell)%geom
    this_order = cell(this_cell)%order
    !
    ! Collect the primitive variables, velocity gradients, jacobian, and
    ! quadrature weights at the error quadrature points for the current cell
    !
    if (allocated(outerp(this_geom,this_order)%error)) then
      !
      nep = size(outerp(this_geom,this_order)%error%mat,dim=1)
      !
      ! The error quadrature points are different from the solution points so
      ! we must interpolate these values from the solution points to the
      ! error quadrature points
      !
      ! NOTE: The vsp array will temporarily hold the conservative variables
      do m = 1,nq
        call outerp_mv_mult( outerp(this_geom,this_order)%error , &
                             usp(n1:n2,m) , vsp(m,1:nep) )
      end do
      !
      do m = 1,nq
        do l = 1,nr
          call outerp_mv_mult( outerp(this_geom,this_order)%error , &
                               dusp(n1:n2,l,m) , dcvdx(l,m,1:nep) )
        end do
      end do
      !
      call outerp_mv_mult( outerp(this_geom,this_order)%error , &
                           jac(n1:n2) , wtdjac(1:nep) )
      !
      wtdjac(1:nep) = wtdjac(1:nep) * outerp(this_geom,this_order)%err_wts
      !
    else
      !
      nep = n2-n1+1
      !
      ! The error quadrature points and solution points are co-located so
      ! just copy these values into the local arrays.
      ! NOTE: The vsp array will temporarily hold the conservative variables
      !
      do n = 1,nep
        do m = 1,nq
          vsp(m,n) = usp(n1-1+n,m)
          do l = 1,nr
            dcvdx(l,m,n) = dusp(n1-1+n,l,m)
          end do
        end do
      end do
     !vsp(1:nq,1:nep) = usp(1:nq,n1:n2)
     !dcvdx(1:nr,1:nq,1:nep) = dusp(1:nr,1:nq,n1:n2)
      wtdjac(1:nep) = jac(n1:n2) * std_elem(this_geom,this_order)%wts
      !
    end if
    !
    ! Convert the conservative variable derivatives
    ! to primitive variable derivatives
    ! NOTE: The vsp array currently contains the conservative variables
    !
    do k = 1,nep
      dpvdx(1:nr,1:nq,k) = grad_cv_to_grad_pv_sp( vsp(1:nq,k) , &
                                                  dcvdx(1:nr,1:nq,k) )
    end do
    !
    ! Convert the conservative variables in vsp to primitive variables
    !
    vsp(1:nq,1:nep) = usp2v( vsp(1:nq,1:nep) )
    !
    ! Loop through the solution points of the current
    !
    do k = 1,nep
      !
      ! Compute the divergence of the velocity
      !
      div_vel = dpvdx(1,nmx,k) + dpvdx(2,nmy,k) + dpvdx(3,nmz,k)
      !
      ! Compute the deviatoric strain rate tensor
      !
      do m = 1,nr
        do l = 1,nr
          strain(l,m) = half * (dpvdx(l,nmb-1+m,k) + dpvdx(m,nmb-1+l,k))
        end do
        strain(m,m) = strain(m,m) - one3*div_vel
      end do
      !
      ! Compute the vorticity vector
      !
      vort(1) = dpvdx(2,nmz,k) - dpvdx(3,nmy,k) ! dw/dy - dv/dz
      vort(2) = dpvdx(3,nmx,k) - dpvdx(1,nmz,k) ! du/dz - dw/dx
      vort(3) = dpvdx(1,nmy,k) - dpvdx(2,nmx,k) ! dv/dx - du/dy
      !
      ! Compute the different kinetic energy dissipation rates (KEDR)
      !
      ! Theoretical KEDR based on the integrated enstrophy
      !
      sum_kedr1(1,tgv_idx) = sum_kedr1(1,tgv_idx) + &
                            wtdjac(k)*vsp(nec,k)*selfdot( vort )
      !
      ! Contribution to theoretical KEDR from the deviatoric strain rate tensor
      !
      sum_kedr1(2,tgv_idx) = sum_kedr1(2,tgv_idx) + &
                            wtdjac(k)*sum(strain*strain)
      !
      ! Contribution to theoretical KEDR from the pressure dialation term
      !
      sum_kedr1(3,tgv_idx) = sum_kedr1(3,tgv_idx) + &
                            wtdjac(k)*vsp(nee,k)*div_vel
      !
      ! Integrated kinetic energy (later, it will
      ! be differenced in time to get KEDR)
      !
      sum_kedr1(4,tgv_idx) = sum_kedr1(4,tgv_idx) + &
                            wtdjac(k)*vsp(nec,k)*selfdot( vsp(nmb:nme,k) )
      !
    end do
    !
  end do
  !
  !#####################################################
  !#####################################################
  !#####################################################
  !##########                                 ##########
  !##########     INTERVAL OUTPUT SECTION     ##########
  !##########                                 ##########
  !#####################################################
  !#####################################################
  !#####################################################
  !
  ! If this is an iteration to output the results, collect the KEDR values
  ! since the last output and dump the results to the Taylor-Green results file.
  !
  if (this_is_output_iter) then
    !
    ! Sum up the KEDR values across all the processors
    !
    if (ncpu > 1) then
      !
      call mpi_allreduce(MPI_IN_PLACE,sum_kedr1,size(sum_kedr1,kind=int_mpi), &
                         mpi_flttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
      !
    end if
    !
    ! If this this the first output interval, change nmin to 0 so that the
    ! KEDR values for the initial conditions are also modified
    !
    if (first_output) nmin = 0
    !
    ! Finish computing the different KEDR
    !
    tgv_kedr1(2:5,nmin:nmax) = sum_kedr1(1:4,nmin:nmax) / total_volume
    !
    ! Copy tgv_kedr(2) into tgv_kedr(6)
    !
    tgv_kedr1(6,nmin:nmax) = tgv_kedr1(2,nmin:nmax)
    !
    ! Adjust the nondimensionalization for some of these KEDR values
    !
    do n = nmin,nmax
      tgv_kedr1(2,n) = tgv_kedr1(2,n) * tgv_dimen1(2)
      tgv_kedr1(3,n) = tgv_kedr1(3,n) * tgv_dimen1(3)
      tgv_kedr1(4,n) = tgv_kedr1(4,n) * tgv_dimen1(4)
      tgv_kedr1(5,n) = tgv_kedr1(5,n) * tgv_dimen1(5)
      tgv_kedr1(6,n) = tgv_kedr1(6,n) * tgv_dimen1(6)
    end do
    !
    ! Compute the central difference of the KEDR to get its rate of change
    ! over time
    ! NOTE: The KEDR found by differencing the integrated kinetic energy is
    !       actually given as KEDR(Ek) = -d/dt (Ek) so the differencing
    !       equation below has been multiplied through by a -1.
    !
    nbeg = merge(     1     ,   0    , first_output )
    nend = merge( tgv_idx-1 , nmax-1 , final_output )
    !
    do n = nbeg,nend
      !
      kem = tgv_kedr1(nke,n-1)
      kec = tgv_kedr1(nke,n)
      kep = tgv_kedr1(nke,n+1)
      !
      dtm = tgv_kedr1(ndt,n)
      dtp = tgv_kedr1(ndt,n+1)
      !
      tgv_kedr1(1,n) =  ( (dtm/dtp)*(kec-kep) + &
                          (dtp/dtm)*(kem-kec) ) / (dtm+dtp)
      !
    end do
    !
    ! Do a one sided difference on the initial conditions if this is the
    ! first output interval
    !
    if (first_output) then
      n = 0
      tgv_kedr1(1,n) = (tgv_kedr1(nke,n)-tgv_kedr1(nke,n+1))/tgv_kedr1(ndt,n+1)
    end if
    !
    ! If this is the last output interval, do a one sided difference
    ! on the final solution
    !
    if (final_output) then
      n = tgv_idx
      tgv_kedr1(1,n) = (tgv_kedr1(nke,n-1)-tgv_kedr1(nke,n)) / tgv_kedr1(ndt,n)
    end if
    !
    ! Have the root processor output the Taylor-Green KEDR solutions
    !
    if (mypnum == 0) then
      !
      solution_time = tgv_kedr1(ndt,-1)
      do n = 0,nend
        solution_time = solution_time + tgv_kedr1(ndt,n)
        write (iotg1,2) solution_time, &
                        (tgv_kedr1(k,n), k=1,size(tgv_kedr1,dim=1))
      end do
      !
      ! Perform a few final tasks if this is the final output
      !
      if (final_output) then
        !
        ! Output the last value of tgv_kedr
        !
        n = tgv_idx
        solution_time = solution_time + tgv_kedr1(ndt,n)
        write (iotg1,2) solution_time, &
                        (tgv_kedr1(k,n), k=1,size(tgv_kedr1,dim=1))
        !
        ! Close the Taylor-Green KEDR solution file
        !
        close (iotg1,iostat=ierr,iomsg=error_message)
        call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
        !
      end if
      !
    else
      !
      solution_time = sum( tgv_kedr1(ndt,-1:nend) )
      !
    end if
    !
    ! Peform some maintence on the tgv_kedr array before continuing
    !
    if (final_output) then
      !
      ! Deallocate tgv_kedr since it is no longer needed
      !
      deallocate ( tgv_kedr1 , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"tgv_kedr1",2,__LINE__,__FILE__,ierr, &
                       error_message,skip_alloc_pause)
      !
      deallocate ( sum_kedr1 , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"sum_kedr1",2,__LINE__,__FILE__,ierr, &
                       error_message,skip_alloc_pause)
      !
    else
      !
      ! Copy tgv_kedr(:,nmax-1:nmax) to tgv_kedr(:,-1:0) so they can
      ! be used for the next output interval. Set tgv_kedr(ndt,-1) to
      ! the last solution time that was output. Initialize all the other
      ! values to zero.
      !
      tgv_kedr1(:,-1) = [tgv_kedr1(1:6,nmax-1),solution_time]
      tgv_kedr1(:,0)  = tgv_kedr1(:,nmax)
      tgv_kedr1(:,1:nmax) = zero
      sum_kedr1 = zero
      !
    end if
    !
    ! Set completed_first_output to true since
    ! we have now completed at least one output
    !
    completed_first_output = true
    !
  end if
  !
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ('VARIABLES = "t*"',/, &
            '"<greek>e</greek>(E<sub>k</sub>)"',/, &
            '"<greek>e</greek>(<greek>z</greek>)"',/, &
            '"<greek>e</greek><sub>1</sub> (Deviatoric Strain Rate)"',/, &
            '"<greek>e</greek><sub>3</sub> (Pressure Dialation)"',/, &
            '"E<sub>k</sub>"',/, &
            '"<greek>z</greek>"',/, &
            '"<greek>D</greek>t*"',/, &
            'ZONE T="P',i0,' - ',i0,'<sup>3</sup>"',/, &
            'AUXDATA order="',i0,'"',/, &
            'AUXDATA RDOF="',i0,'<sup>3</sup>"',/, &
            'AUXDATA DOF="',i0,'<sup>3</sup>"')
  2 format (100(es17.9))
  !
end subroutine compute_TaylorGreen_KEDR
!
!###############################################################################
!
subroutine generic_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
  !
  !.. Formal Arguments ..
  integer, optional, intent(in) :: iunit
  !
  !.. Local Scalars ..
  integer :: iou
  character(len=36) :: array_name
  type(memory) :: array
  type(memory) :: total
  !
#ifndef DISABLE_DTIO
continue
  !
  iou = iout
  if (present(iunit)) iou = iunit
  !
  call total%reset
  !
  write (iou,1)
  !
  if (allocated(tgv_kedr)) then
    call array%set( storage_size(tgv_kedr,kind=inttype) , &
                            size(tgv_kedr,kind=inttype) )
    call total%add(array)
    write (array_name,5) "tgv_kedr"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "tgv_kedr"
  end if
  if (allocated(sum_kedr)) then
    call array%set( storage_size(sum_kedr,kind=inttype) , &
                            size(sum_kedr,kind=inttype) )
    call total%add(array)
    write (array_name,5) "sum_kedr"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "sum_kedr"
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module generic_mod")
  2 format (" Array: ",a7,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine generic_memory_usage
!
!###############################################################################
!
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
include "Functions/usp2v.f90"
include "Functions/entropy.f90"
include "Functions/pressure.f90"
include "Functions/vortex.f90"
include "Functions/grad_cv_to_grad_pv.f90"
!
!###############################################################################
!
end module generic_mod
