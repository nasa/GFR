module generic_mod
  !
!#define PROFILE_ON
  !
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  use eqn_idx, only : nec,nmb,nme,nmx,nmy,nmz,nee,nq
  use eqn_idx, only : ncs,ncb,nce,npres,ntemp,ngamm
  !
  implicit none
  !
  private
  !
  public :: init_generic_mod
  public :: write_iter_stats
  public :: getmaxent
  public :: compute_residual_error
  public :: compute_KEDR
  public :: generic_memory_usage
  !
  character(len=100), save :: iter_out_fmt
  !
  real(wp),    save :: Density_Maximum_Residual = -huge(zero)
  logical(lk), save :: vortex_has_collapsed = fals
  integer,     save :: final_countdown = 0
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
subroutine compute_residual_error()
  !
  !.. Use Statements ..
  use geovar,       only : nr,n_solpts,ncell,cell,cell_map
  use ovar,         only : total_volume,uverr
  use ovar,         only : use_async_for_error
  use ovar,         only : itcur
  use flowvar,      only : usp,uoldsp
  !
  !.. Local Scalars ..
  integer  :: k,l,m,n,nc,n1,n2,l1,l2,l3,e1,e2,e3
  real(wp) :: vol_inv
  !
  integer :: io
  character(len=20) :: fname
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
  real(wp), allocatable :: unew(:,:)
  real(wp), allocatable :: uold(:,:)
  real(wp), allocatable :: error(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_residual_error"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  l1 = 1
  l2 = min(2,nr)
  l3 = min(3,nr)
  !
  e1 = size(uverr,dim=1) ! nq+1
  e2 = size(uverr,dim=2) ! nr+5
  e3 = size(uverr,dim=3) ! 2
  !
  allocate ( error(1:e1,1:e2,1:e3) , source=zero )
  allocate ( unew(1:n_solpts,1:nq+1) , source=zero )
  allocate ( uold(1:n_solpts,1:nq+1) , source=zero )
  allocate ( cell_for_sp(1:n_solpts) , source=0 )
  !
  !$ACC DATA CREATE(error,unew,uold,cell_for_sp)
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
  call get_unew_solpts_gpu(unew)
  call get_uold_solpts_gpu(uold)
  !
!#define OUTPUT_U_ARRAYS
#ifdef OUTPUT_U_ARRAYS
  200 format (i7,2x,10es18.10)
  201 format (a,".",i0,".dat")
  !
  if (itcur == 1) then
    write (fname,201) "usp",0
    open (newunit=io,file=fname,status="replace",action="write")
    do n = 1,n_solpts
      write (io,200) n,(uoldsp(n,m),m=1,size(uoldsp,dim=2))
    end do
    close (io)
  end if
  !
  write (fname,201) "usp",itcur
  open (newunit=io,file=fname,status="replace",action="write")
  do n = 1,n_solpts
    write (io,200) n,(usp(n,m),m=1,size(usp,dim=2))
  end do
  close (io)
  !
  write (fname,201) "unew",itcur
  open (newunit=io,file=fname,status="replace",action="write")
  do n = 1,n_solpts
    write (io,200) n,(unew(n,m),m=1,size(unew,dim=2))
  end do
  close (io)
  !
  write (fname,201) "uold",itcur
  open (newunit=io,file=fname,status="replace",action="write")
  do n = 1,n_solpts
    write (io,200) n,(uold(n,m),m=1,size(uold,dim=2))
  end do
  close (io)
#endif
  !
  if (use_async_for_error) then
    call compute_error_gpu_async(unew,uold,error(:,:,1))
  else
    call compute_error_gpu(unew,uold,error(:,:,1),1)
  end if
  !
  call get_vnew_solpts_gpu(unew)
  call get_vold_solpts_gpu(uold)
  !
#ifdef OUTPUT_U_ARRAYS
  write (fname,201) "vnew",itcur
  open (newunit=io,file=fname,status="replace",action="write")
  do n = 1,n_solpts
    write (io,200) n,(unew(n,m),m=1,size(unew,dim=2))
  end do
  close (io)
  !
  write (fname,201) "vold",itcur
  open (newunit=io,file=fname,status="replace",action="write")
  do n = 1,n_solpts
    write (io,200) n,(uold(n,m),m=1,size(uold,dim=2))
  end do
  close (io)
#endif
  !
  if (use_async_for_error) then
    call compute_error_gpu_async(unew,uold,error(:,:,2))
  else
    call compute_error_gpu(unew,uold,error(:,:,2),2)
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
   !if (mypnum == glb_root) then
   !  do l = 1,2
   !    do m = 1,nq+1
   !      write (150,150) (all_error(m,1,l,n), n=1,ncpu), &
   !                      sum(all_error(m,2,l,:)) * vol_inv
   !      write (151,150) (all_error(m,2,l,n), n=1,ncpu), &
   !                      sqrt( sum(all_error(m,2,l,:)) * vol_inv )
   !    end do
   !    write (150,*)
   !    write (151,*)
   !  end do
   !  write (150,*)
   !  write (151,*)
   !  150 format (20es14.6)
   !  flush (150)
   !  flush (151)
   !end if
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
  deallocate ( uold )
  deallocate ( unew )
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
end subroutine compute_residual_error
!
!###############################################################################
!
subroutine compute_error_gpu(unew,uold,error,cv_pv_switch)
  !
  !.. Use Statements ..
  use geovar,      only : nr,n_solpts,xyz
  use metrics_mod, only : wtdjac
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: unew(:,:)
  real(wp),    intent(in) :: uold(:,:)
  real(wp), intent(inout) :: error(:,:)
  integer,     intent(in) :: cv_pv_switch
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
  character(len=*), parameter :: cvar(6,3) = &
      reshape( ["rho ","rhou","rhoE"," s  ","    ","    ", &
                "rho ","rhou","rhov","rhoE"," s  ","    ", &
                "rho ","rhou","rhov","rhow","rhoE"," s  "], shape(cvar) )
  character(len=*), parameter :: pvar(6,3) = &
      reshape( ["rho "," vx ","pres","vmag","    ","    ", &
                "rho "," vx "," vy ","pres","vmag","    ", &
                "rho "," vx "," vy "," vz ","pres","vmag"], shape(pvar) )
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
    !$ACC&              PRESENT(unew,uold,wtdjac,linf_error)
    do n = 1,n_solpts
      !
      new_val = unew(n,m)
      old_val = uold(n,m)
      !
     !linf = abs( new_val - old_val )
      if (m == nq+1 .and. cv_pv_switch == 2) then
        linf = norm2( abs(unew(n,nmb:nme) - uold(n,nmb:nme)) )
      else
        linf = abs( unew(n,m) - uold(n,m) )
      end if
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
    if (cv_pv_switch == 1) then
      write (iout,100) "lone_sum",cvar(m,nr),lone_sum
      write (iout,100) "ltwo_sum",cvar(m,nr),ltwo_sum
    else
      write (iout,100) "lone_sum",pvar(m,nr),lone_sum
      write (iout,100) "ltwo_sum",pvar(m,nr),ltwo_sum
    end if
    100 format (a,"(",a,") = ",es16.8)
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
subroutine compute_error_gpu_async(unew,uold,error)
  !
  !.. Use Statements ..
  use geovar,      only : nr,n_solpts,xyz
  use metrics_mod, only : wtdjac
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: unew(:,:)
  real(wp),    intent(in) :: uold(:,:)
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
  !$ACC&              PRESENT(unew,uold,linf_error)
  do m = 1,nq+1
    do n = 1,n_solpts
      linf_error(n,m) = abs( unew(n,m) - uold(n,m) )
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
subroutine get_unew_solpts_gpu(unew)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : usp,uavesp,gpsp
  use ovar,    only : itestcase,pref_nd,gam
  use ovar,    only : output_time_averaging
  use ovar,    only : force_nd_subroutines
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: unew(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: rho,pres,gamm
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
  if (output_time_averaging) then
    !
    if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
      !$ACC PARALLEL LOOP PRESENT(uavesp,unew)
      do n = 1,n_solpts
        unew(n,nec) = uavesp(n,nec)
        unew(n,nmx) = uavesp(n,nmx)
        unew(n,nmy) = uavesp(n,nmy)
        unew(n,nmz) = uavesp(n,nmz)
        unew(n,nee) = uavesp(n,nee)
      end do
    else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
      !$ACC PARALLEL LOOP PRESENT(uavesp,unew)
      do n = 1,n_solpts
        unew(n,nec) = uavesp(n,nec)
        unew(n,nmx) = uavesp(n,nmx)
        unew(n,nmy) = uavesp(n,nmy)
        unew(n,nee) = uavesp(n,nee)
      end do
    else
      !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(uavesp,unew)
      do m = 1,nq
        do n = 1,n_solpts
          unew(n,m) = uavesp(n,m)
        end do
      end do
    end if
    !
  else
    !
    if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
      !$ACC PARALLEL LOOP PRESENT(usp,unew)
      do n = 1,n_solpts
        unew(n,nec) = usp(n,nec)
        unew(n,nmx) = usp(n,nmx)
        unew(n,nmy) = usp(n,nmy)
        unew(n,nmz) = usp(n,nmz)
        unew(n,nee) = usp(n,nee)
      end do
    else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
      !$ACC PARALLEL LOOP PRESENT(usp,unew)
      do n = 1,n_solpts
        unew(n,nec) = usp(n,nec)
        unew(n,nmx) = usp(n,nmx)
        unew(n,nmy) = usp(n,nmy)
        unew(n,nee) = usp(n,nee)
      end do
    else
      !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,unew)
      do m = 1,nq
        do n = 1,n_solpts
          unew(n,m) = usp(n,m)
        end do
      end do
    end if
    !
  end if
  !
  if (ncs == 0) then
    !$ACC PARALLEL LOOP PRESENT(gpsp)
    do n = 1,n_solpts
      gpsp(n,ngamm) = gam
    end do
  end if
  !
  if (itestcase == Inviscid_Gaussian_Bump) then
    !
    !$ACC PARALLEL LOOP PRESENT(usp,gpsp,unew)
    do n = 1,n_solpts
      !
      rho  = usp(n,nec)
      pres = gpsp(n,npres)/pref_nd
      gamm = gpsp(n,ngamm)
      !
      unew(n,nq+1) = pres / (rho**gamm)
      !
    end do
    !
  else
    !
    !$ACC PARALLEL LOOP PRESENT(usp,gpsp,unew)
    do n = 1,n_solpts
      !
      rho = usp(n,nec)
      pres = gpsp(n,npres)
      gamm = gpsp(n,ngamm)
      !
      unew(n,nq+1) = log(pres) - gamm*log(rho)
      !
    end do
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_unew_solpts_gpu
!
!###############################################################################
!
subroutine get_uold_solpts_gpu(uold)
  !
  !.. Use Statements ..
  use ovar,    only : prev_ave_time,saved_ave_time,ave_start_time
  use ovar,    only : output_time_averaging,read_time_ave_restart
  use ovar,    only : itcur,itrst,bc_in,itestcase,mms_opt
  use flowvar, only : usp,uoldsp,uavesp,gpoldsp
  use geovar,  only : n_solpts,xyz
  use mms_mod, only : mms_solution_at_xyz_gpu
  use gpu_mod, only : vortex_solution_gpu
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: uold(:,:)
  !
  !.. Local Scalars ..
  integer     :: m,n,nr
  real(wp)    :: a,b,rho,pres,gamm
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
  if (mms_opt /= 0) then
    !
    ! Evaluate the exact analytical MMS solution at the error points
    !
    call mms_solution_at_xyz_gpu(1,n_solpts,xyz,uold, &
                                 return_nondim=true, &
                                 return_cv=true)
    !
  else if (any(itestcase == Transport_Problems)) then
    !
    ! Evaluate the exact analytical solution to the vortex transport problem
    !
    call vortex_solution_gpu(1,n_solpts,xyz,uold)
    !
  else if (itestcase == Freestream_Preservation) then
    !
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(bc_in,uold)
    do m = 1,nq
      do n = 1,n_solpts
        uold(n,m) = bc_in(0)%cv(m)
      end do
    end do
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
      !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uavesp,uold)
      do m = 1,nq
        do n = 1,n_solpts
          uold(n,m) = b * uavesp(n,m) - a * usp(n,m)
        end do
      end do
      !
    else
      !
      !$ACC PARALLEL LOOP PRESENT(uoldsp,uold)
      do m = 1,nq
        do n = 1,n_solpts
          uold(n,m) = uoldsp(n,m)
        end do
      end do
      !
    end if
    !
  end if
  !
  if (itestcase == Inviscid_Gaussian_Bump) then
    !
    !$ACC PARALLEL LOOP PRESENT(uold)
    do n = 1,n_solpts
      uold(n,nq+1) = one
    end do
    !
  else
    !
    !$ACC PARALLEL LOOP PRESENT(uoldsp,gpoldsp,uold)
    do n = 1,n_solpts
      !
      rho  = uoldsp(n,nec)
      pres = gpoldsp(n,npres)
      gamm = gpoldsp(n,ngamm)
      !
      uold(n,nq+1) = log(pres) - gamm*log(rho)
      !
    end do
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_uold_solpts_gpu
!
!###############################################################################
!
subroutine get_vnew_solpts_gpu(vnew)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : uavesp,gpsp
  use ovar,    only : output_time_averaging
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: vnew(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: velx,vely,velz,rho_inv,pres,temp
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
  ! Get the primitive variables of the current solution at the error points
  !
  if (output_time_averaging) then
    if (nr == 3) then
      !$ACC PARALLEL LOOP PRESENT(uavesp,vnew)
      do n = 1,n_solpts
        !
        velx = uavesp(n,nq+nmx)
        vely = uavesp(n,nq+nmy)
        velz = uavesp(n,nq+nmz)
        !
        vnew(n,nec) = uavesp(n,nq+nec)
        vnew(n,nmx) = velx
        vnew(n,nmy) = vely
        vnew(n,nmz) = velz
        vnew(n,nee) = uavesp(n,nq+nee)
        vnew(n,nq+1) = sqrt(velx*velx + vely*vely + velz*velz)
        !
      end do
    else if (nr == 2) then
      !$ACC PARALLEL LOOP PRESENT(uavesp,vnew)
      do n = 1,n_solpts
        !
        velx = uavesp(n,nq+nmx)
        vely = uavesp(n,nq+nmy)
        !
        vnew(n,nec) = uavesp(n,nq+nec)
        vnew(n,nmx) = velx
        vnew(n,nmy) = vely
        vnew(n,nee) = uavesp(n,nq+nee)
        vnew(n,nq+1) = sqrt(velx*velx + vely*vely)
        !
      end do
    else
      !$ACC PARALLEL LOOP PRESENT(uavesp,vnew)
      do n = 1,n_solpts
        !
        velx = uavesp(n,nq+nmx)
        !
        vnew(n,nec) = uavesp(n,nq+nec)
        vnew(n,nmx) = velx
        vnew(n,nee) = uavesp(n,nq+nee)
        vnew(n,nq+1) = sqrt(velx*velx)
        !
      end do
    end if
    !
    if (ncs > 0) then
      !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(uavesp,vnew)
      do m = 1,ncs
        do n = 1,n_solpts
          vnew(n,m) = uavesp(n,nq+m)
        end do
      end do
    end if
    !
  else
    !
    ! Need to first convert species densities to mass fractions
    ! before temperature is saved in the density location of vnew
    !
    if (ncs > 0) then
      !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(vnew)
      do m = 1,ncs
        do n = 1,n_solpts
          vnew(n,m) = vnew(n,m)/vnew(n,nec)
        end do
      end do
    end if
    !
    if (nr == 3) then
      !$ACC PARALLEL LOOP PRESENT(gpsp,vnew)
      do n = 1,n_solpts
        !
        rho_inv = one/vnew(n,nec)
        !
        velx = vnew(n,nmx)*rho_inv
        vely = vnew(n,nmy)*rho_inv
        velz = vnew(n,nmz)*rho_inv
        !
        vnew(n,nec) = gpsp(n,ntemp)
        vnew(n,nmx) = velx
        vnew(n,nmy) = vely
        vnew(n,nmz) = velz
        vnew(n,nee) = gpsp(n,npres)
        vnew(n,nq+1) = sqrt(velx*velx + vely*vely + velz*velz)
        !
      end do
    else if (nr == 2) then
      !$ACC PARALLEL LOOP PRESENT(gpsp,vnew)
      do n = 1,n_solpts
        !
        rho_inv = one/vnew(n,nec)
        !
        velx = vnew(n,nmx)*rho_inv
        vely = vnew(n,nmy)*rho_inv
        !
        vnew(n,nec) = gpsp(n,ntemp)
        vnew(n,nmx) = velx
        vnew(n,nmy) = vely
        vnew(n,nee) = gpsp(n,npres)
        vnew(n,nq+1) = sqrt(velx*velx + vely*vely)
        !
      end do
    else
      !$ACC PARALLEL LOOP PRESENT(gpsp,vnew)
      do n = 1,n_solpts
        !
        velx = vnew(n,nmx)/vnew(n,nec)
        !
        vnew(n,nec) = gpsp(n,ntemp)
        vnew(n,nmx) = velx
        vnew(n,nee) = gpsp(n,npres)
        vnew(n,nq+1) = abs(velx)
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
subroutine get_vold_solpts_gpu(vold)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : usp,uavesp,gpsp,gpoldsp
  use ovar,    only : prev_ave_time,saved_ave_time,ave_start_time
  use ovar,    only : output_time_averaging,read_time_ave_restart
  use ovar,    only : itcur,itrst
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: vold(:,:)
  !
  !.. Local Scalars ..
  integer     :: m,n,nr
  real(wp)    :: a,b,velx,vely,velz,q2,rho_inv
  logical(lk) :: use_uavesp
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
    ! Need to first convert species densities to mass fractions
    ! before temperature is saved in the density location of vold
    !
    if (ncs > 0) then
      !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uavesp,vold)
      do m = ncb,nce
        do n = 1,n_solpts
          vold(n,m) = b * uavesp(n,nq+m) - a * usp(n,m)/usp(n,nec)
        end do
      end do
    end if
    !
    if (nr == 3) then
      !$ACC PARALLEL LOOP PRESENT(usp,uavesp,gpsp,vold)
      do n = 1,n_solpts
        !
        rho_inv = one/usp(n,nec)
        !
        velx = usp(n,nmx)*rho_inv
        vely = usp(n,nmy)*rho_inv
        velz = usp(n,nmz)*rho_inv
        !
        q2 = velx*velx + vely*vely + velz*velz
        !
        velx = b * uavesp(n,nq+nmx) - a * velx
        vely = b * uavesp(n,nq+nmy) - a * vely
        velz = b * uavesp(n,nq+nmz) - a * velz
        !
        q2 = b * sqrt(velx*velx + vely*vely + velz*velz) - a * sqrt(q2)
        !
        vold(n,nec) = b * uavesp(n,nq+nec) - a * gpsp(n,ntemp)
        vold(n,nmx) = velx
        vold(n,nmy) = vely
        vold(n,nmz) = velz
        vold(n,nee) = b * uavesp(n,nq+nee) - a * gpsp(n,npres)
        vold(n,nq+1) = q2
        !
      end do
    else if (nr == 2) then
      !$ACC PARALLEL LOOP PRESENT(usp,uavesp,gpsp,vold)
      do n = 1,n_solpts
        !
        rho_inv = one/usp(n,nec)
        !
        velx = usp(n,nmx)*rho_inv
        vely = usp(n,nmy)*rho_inv
        !
        q2 = velx*velx + vely*vely
        !
        velx = b * uavesp(n,nq+nmx) - a * velx
        vely = b * uavesp(n,nq+nmy) - a * vely
        !
        q2 = b * sqrt(velx*velx + vely*vely) - a * sqrt(q2)
        !
        vold(n,nec) = b * uavesp(n,nq+nec) - a * gpsp(n,ntemp)
        vold(n,nmx) = velx
        vold(n,nmy) = vely
        vold(n,nee) = b * uavesp(n,nq+nee) - a * gpsp(n,npres)
        vold(n,nq+1) = q2
        !
      end do
    else
      !$ACC PARALLEL LOOP PRESENT(usp,uavesp,gpsp,vold)
      do n = 1,n_solpts
        !
        velx = usp(n,nmx)/usp(n,nec)
        !
        q2 = velx*velx
        !
        velx = b * uavesp(n,nq+nmx) - a * velx
        !
        q2 = b * sqrt(velx*velx) - a * sqrt(q2)
        !
        vold(n,nec) = b * uavesp(n,nq+nec) - a * gpsp(n,ntemp)
        vold(n,nmx) = velx
        vold(n,nee) = b * uavesp(n,nq+nee) - a * gpsp(n,npres)
        vold(n,nq+1) = q2
        !
      end do
    end if
    !
  else
    !
    ! Need to first convert species densities to mass fractions
    ! before temperature is saved in the density location of vold
    !
    if (ncs > 0) then
      !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(vold)
      do m = ncb,nce
        do n = 1,n_solpts
          vold(n,m) = vold(n,m) / vold(n,nec)
        end do
      end do
    end if
    !
    if (nr == 3) then
      !$ACC PARALLEL LOOP PRESENT(vold,gpoldsp)
      do n = 1,n_solpts
        !
        rho_inv = one/vold(n,nec)
        !
        velx = vold(n,nmx)*rho_inv
        vely = vold(n,nmy)*rho_inv
        velz = vold(n,nmz)*rho_inv
        !
        vold(n,nec) = gpoldsp(n,ntemp)
        vold(n,nmx) = velx
        vold(n,nmy) = vely
        vold(n,nmz) = velz
        vold(n,nee) = gpoldsp(n,npres)
        vold(n,nq+1) = sqrt(velx*velx + vely*vely + velz*velz)
        !
      end do
    else if (nr == 2) then
      !$ACC PARALLEL LOOP PRESENT(vold,gpoldsp)
      do n = 1,n_solpts
        !
        rho_inv = one/vold(n,nec)
        !
        velx = vold(n,nmx)*rho_inv
        vely = vold(n,nmy)*rho_inv
        !
        vold(n,nec) = gpoldsp(n,ntemp)
        vold(n,nmx) = velx
        vold(n,nmy) = vely
        vold(n,nee) = gpoldsp(n,npres)
        vold(n,nq+1) = sqrt(velx*velx + vely*vely)
        !
      end do
    else
      !$ACC PARALLEL LOOP PRESENT(vold,gpoldsp)
      do n = 1,n_solpts
        !
        velx = vold(n,nmx)/vold(n,nec)
        !
        vold(n,nec) = gpoldsp(n,ntemp)
        vold(n,nmx) = velx
        vold(n,nee) = gpoldsp(n,npres)
        vold(n,nq+1) = abs(velx)
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
subroutine write_iter_stats(walltime_expiring,exit_main_loop)
  !
  !.. Use Statements ..
  use geovar,  only : global_pinc_ptr,n_global_cell
  use ovar,    only : itcur,d_t,cfl,time,minimum_timesteps,time_ref
  use ovar,    only : convergence_reached,itestcase,Period_Timesteps
  use ovar,    only : output_limiter_statistics,nlim
  use ovar,    only : convergence_order_max_res,convergence_order_abs_res
  use ovar,    only : rl2mx,uverr,this_is_final_timestep,iter_out_interval
  use ovar,    only : output_time_averaging,time_scaling_factor,ave_start_time
  use ovar,    only : max_restart_wall_time,max_time_ave_rst_wall_time
  use ovar,    only : max_cgns_wall_time,max_time_ave_cgns_wall_time
  !
  use order_mod, only : n_order
  !
  use pbs_mod, only : pbs_remaining_walltime
  !
  use thx_mod, only : compute_blowing_ratio
  !
  !.. Formal Arguments ..
  logical(lk), intent(inout) :: walltime_expiring
  logical(lk), intent(  out) :: exit_main_loop
  !
  !.. Local Scalars ..
  integer  :: hrs_left,min_left,sec_left
  integer  :: n,nloc,kloc,ierr
  real(wp) :: xentmax,yentmax
  real(wp) :: entL2
  real(wp) :: wall_time_left,out_time
  real(wp) :: orders_converged,output_L2val
  real(wp) :: rltwo,rlim
  real(wp) :: blowing_ratio,mdot_nozzle,mdot_injector
  !
  !.. Local Arrays ..
  real(wp) :: max_wall_times(1:4)
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
  if (abs(itestcase) == THX_Phase_3) then
    call compute_blowing_ratio(blowing_ratio,mdot_nozzle,mdot_injector)
    uverr(nmx,2,2) = mdot_nozzle
    uverr(nmy,2,2) = mdot_injector
    uverr(nmz,2,2) = blowing_ratio
  end if
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
  if (size(uverr,dim=2) >= 7) then
    yentmax = uverr(nq+1,7,1)
  else
    yentmax = zero
  end if
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
    orders_converged = log10(max(rl2mx,eps17)) - log10(max(rltwo,eps17))
    output_L2val = orders_converged
    !
    ! Check for convergence
    !
    if (itcur > minimum_timesteps) then
      convergence_reached = (-orders_converged <= convergence_order_max_res) &
                            .or. (log10(rltwo) <= convergence_order_abs_res)
      if (any(itestcase == Taylor_Problems)) then
        convergence_reached = fals
      end if
    end if
    !
  end if
  !
  if (output_limiter_statistics) then
    !
    if (.not. allocated(nlim)) then
      allocate ( nlim(0:n_order) , source=0 , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"nlim",1,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    if (ncpu > 1) then
      call mpi_allreduce(MPI_IN_PLACE,nlim,size(nlim,kind=int_mpi), &
                         mpi_inttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
    end if
    !
  end if
  !
  ! Output the results for the current time step
  !
  max_wall_times(1:4) = [ max_restart_wall_time, &
                          max_time_ave_rst_wall_time, &
                          max_cgns_wall_time, &
                          max_time_ave_cgns_wall_time ]
  !
  if (.not. output_time_averaging) then
    max_wall_times(2) = -one
    max_wall_times(4) = -one
  end if
  !
  call pbs_remaining_walltime(walltime_expiring, &
                              wall_time_remaining=wall_time_left, &
                              io_time_needed=max_wall_times)
  !
  if (mypnum == glb_root) then
    !
    if (output_limiter_statistics) then
      do n = lbound(nlim,dim=1),ubound(nlim,dim=1)
        if (nlim(n) > 0) then
          rlim = ten*ten * real(nlim(n),kind=wp) / real(n_global_cell,kind=wp)
          write (iout,10) nlim(n),rlim,n
        end if
      end do
    end if
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
    if (abs(itestcase) == THX_Phase_3) then
     !write (iout,iter_out_fmt) itcur,d_t*time_ref,time*time_ref, &
      write (iout,iter_out_fmt) itcur,d_t,out_time, &
                                cfl,output_L2val, &
                                mdot_nozzle,mdot_injector,blowing_ratio, &
                                nloc,global_pinc_ptr(nloc)+kloc, &
                                hrs_left,min_left,sec_left
    else
     !write (iout,iter_out_fmt) itcur,d_t*time_ref,time*time_ref, &
      write (iout,iter_out_fmt) itcur,d_t,out_time, &
                                cfl,rltwo,output_L2val,nloc, &
                                global_pinc_ptr(nloc)+kloc, &
                                xentmax,yentmax, &
                                hrs_left,min_left,sec_left
    end if
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
  10 format (8x,i8," Cells (",f7.3," %) limited to P",i0)
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
  use ovar,   only : uverr,output_time_averaging
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
  if (abs(itestcase) == THX_Phase_3) then
    write (iter_out_fmt,102) final_iter_digits,max_cell_digits,max_solpt_digits
  else
 !if (output_time_averaging) then
 !  write (iter_out_fmt,101) final_iter_digits
 !else
    write (iter_out_fmt,100) final_iter_digits,max_cell_digits,max_solpt_digits
  end if
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
  else if (any(itestcase == Taylor_Problems)) then
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
  102 format ("(1x,i",i0,",7es12.4,2x,i",i0,",2x,i",i0)
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
subroutine compute_KEDR(iter,dtmin)
  !
  !.. Use Statements ..
  use order_mod,        only : n_order
  use geovar,           only : nr,ncell,n_solpts,n_global_solpts
  use geovar,           only : n_global_cell,global_cell_geom
  use ovar,             only : muref_nd,total_volume
  use ovar,             only : time,time_ref,machref
  use ovar,             only : VortexGrid_DOF
  use ovar,             only : iter_out_interval,itrst
  use ovar,             only : itestcase
  use flowvar,          only : usp,dusp,gpsp
  use element_type_mod, only : std_elem
  use metrics_mod,      only : wtdjac
  use gpu_mod,          only : grad_cv_to_grad_pv_gpu
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: iter
  real(wp), intent(in) :: dtmin
  !
  !.. Local Scalars ..
  integer :: k,n,np,n0,ndof_1d,tgv_idx
  integer :: ierr,nmin,nmax,nbeg,nend
  integer :: this_cell,this_geom,this_order
  !
  logical(lk) :: first_output,final_output
  logical(lk) :: this_is_output_iter
  !
  real(wp) :: kem,kec,kep,dtm,dtp,vol_inv
  real(wp) :: solution_time,muref_star
  real(wp) :: ke_coef,kedr_coef
  real(wp) :: sum_kedr_1,sum_kedr_2
  real(wp) :: sum_kedr_3,sum_kedr_4
  real(wp) :: rho,rho_inv,pres,q2,vort,div_vel
  real(wp) :: a1,a2,a3,a4,a5,a6,strain
  real(wp) :: velx,dudx,dvdx,dwdx,vortx
  real(wp) :: vely,dudy,dvdy,dwdy,vorty
  real(wp) :: velz,dudz,dvdz,dwdz,vortz
  !
  character(len=5) :: grid_geom
  !
  character(len=50) :: fname
  !
  !.. Local Saved Scalars ..
  logical(lk), save :: completed_first_output = fals
  logical(lk), save :: needs_initialization = true
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: tdpvdx(:,:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_KEDR"
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
  first_output = (iter == nmax+itrst) .or. &
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
    ! NOTE: This is now being used to compute the KEDR for both the
    !       Taylor-Green vortex (TGV) and Taylor-Couette (TCF) flows.
    !       Both cases nondimensionalize the different KEDR values
    !       by uref, but TCF also nondimensionalizes any velocity
    !       derivative terms by d which is the difference in radius
    !       between the two cylinders. For TGV, time_ref=machref, and
    !       for TCF, time_ref=machref/d. Dividing by machref basically
    !       switches the nondimensionalization of a term from aref to
    !       uref. To keep things the same between TGV and TCF, we
    !       use machref for re-nondimensionalizing all velocity terms
    !       and use time_ref (machref for TGV, and machref/d for TCF)
    !       to re-nondimensionalize all velocity derivative terms.
    ! NOTE: The comments above each value indicate the reason for
    !       the number of machrefs and time_refs
    !
    ! For some reason, ZJs TCF paper doesnt include the 1/2 term
    ! when computing the integrated KE or enstrophy, so we have to 
    ! also create a coefficient to adjust these terms accordingly.
    !
    ke_coef = merge(one,half,abs(itestcase) == Taylor_Couette_Flow)
    !
    ! Create the adjustment coefficient for the finite-difference KEDR.
    ! This is needed for the TCF since the 1/2 coefficient in the KEDR
    ! is implicitly included when tgv_dimen(5) multiplies the KE by 1/2
    ! for the TGV case. Since ZJ doesnt include the 1/2 for the KE for
    ! the TCF case, we need to multiply that in after computing the KEDR.
    !
    kedr_coef = half/ke_coef
    !
    ! Adjust the nondimensionalization of muref so it is
    ! nondimensionalized by uref instead of aref
    muref_star = muref_nd / machref
    !
    tgv_dimen(1) = kedr_coef
    ! Two time_refs for square of velocity derivatives (from vorticity**2)
    tgv_dimen(2) = muref_star / time_ref / time_ref
    ! Two time_refs for square of velocity derivatives (from strain**2)
    tgv_dimen(3) = two * muref_star / time_ref / time_ref
    ! Two machrefs for pressure, one time_ref for velocity derivative
    tgv_dimen(4) = -one / machref / machref / time_ref
    ! Two machrefs for velocity squared in kinetic energy term
    tgv_dimen(5) = ke_coef / machref / machref
    ! Two time_refs for square of velocity derivatives (from vorticity**2)
    tgv_dimen(6) = ke_coef / time_ref / time_ref
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
    if (mypnum == glb_root) then
      !
      if (all(global_cell_geom(1:n_global_cell) == Geom_Tetr)) then
        grid_geom = "Tet"
      else if (all(global_cell_geom(1:n_global_cell) == Geom_Pyra)) then
        grid_geom = "Pyr"
      else if (all(global_cell_geom(1:n_global_cell) == Geom_Pris)) then
        grid_geom = "Pri"
      else if (all(global_cell_geom(1:n_global_cell) == Geom_Hexa)) then
        grid_geom = "Hex"
      else
        grid_geom = "Mixed"
      end if
      !
      ndof_1d = nint( real(n_global_solpts,kind=wp)**(one3) )
      !
      if (abs(itestcase) == Taylor_Couette_Flow) then
        fname = "TaylorCouette_KEDR.dat"
      else
        fname = "TaylorGreen_KEDR.dat"
      end if
      !
      open (newunit=iotg,file=fname,status="replace",action="write", &
            iostat=ierr,iomsg=error_message)
      call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
      !
      write (iotg,1) n_order,ndof_1d,trim(adjustl(grid_geom)), &
                     n_order,ndof_1d,trim(adjustl(grid_geom))
      !
    end if
    !
    ! Unset needs_initialization now that everything has been initialized
    !
    needs_initialization = fals
    !
  end if
  !
  allocate ( tdpvdx(1:n_solpts,1:nr,1:nq) )
  !
  !$ACC DATA CREATE(tdpvdx)
  !
  call grad_cv_to_grad_pv_gpu(1,n_solpts,usp,gpsp,dusp,tdpvdx)
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
  !$ACC&              PRESENT(wtdjac,usp,tdpvdx)
  do n = 1,n_solpts
    !
    rho = usp(n,nec)
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
  !$ACC&              PRESENT(wtdjac,usp,gpsp,tdpvdx)
  do n = 1,n_solpts
    !
    rho = usp(n,nec)
    rho_inv = one/rho
    !
    velx = usp(n,nmx)*rho_inv
    vely = usp(n,nmy)*rho_inv
    velz = usp(n,nmz)*rho_inv
    !
    pres = gpsp(n,npres)
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
     !tgv_kedr(1,n) = tgv_kedr(5,n)-tgv_kedr(5,n-1)          ! e(Ek) KEDR
      tgv_kedr(2,n) = sum_kedr(1,n) * tgv_dimen(2) * vol_inv ! e(Z) Z-KEDR
      tgv_kedr(3,n) = sum_kedr(2,n) * tgv_dimen(3) * vol_inv ! e1 (strain KEDR)
      tgv_kedr(4,n) = sum_kedr(3,n) * tgv_dimen(4) * vol_inv ! e3 (pres KEDR)
      tgv_kedr(5,n) = sum_kedr(4,n) * tgv_dimen(5) * vol_inv ! Ek (KE)
      tgv_kedr(6,n) = sum_kedr(1,n) * tgv_dimen(6) * vol_inv ! Z (Enstrophy)
    end do
    !
    ! Compute the central difference of the KEDR to get its rate of change
    ! over time
    ! NOTE: The KEDR found by differencing the integrated kinetic energy is
    !       actually given as KEDR(Ek) = -d/dt (Ek) so the differencing
    !       equation below has been multiplied through by a -1.
    !
    call compute_dEkdt_central_diff_O2(first_output,final_output,tgv_idx, &
                                       tgv_kedr(5,:),tgv_kedr(7,:), &
                                       tgv_kedr(1,:),nend)
   !call compute_dEkdt_backward_diff_O1(first_output,final_output,tgv_idx, &
   !                                    tgv_kedr(5,:),tgv_kedr(7,:), &
   !                                    tgv_kedr(1,:),nend)
    !
    ! Multiply the finite-difference KEDR for the Taylor-Couette flow by 1/2
    ! since tgv_dimen(5) for KE doesnt include the 1/2 like it does for TGV.
    ! NOTE: This uses kedr_coef to adjust the KEDR for both the TGV and TCF
    !       cases, but kedr_coef=1 for TGV and kedr_coef=1/2 for TCF.
    !
    !$ACC PARALLEL LOOP PRESENT(tgv_kedr)
    do n = 0,nend
      tgv_kedr(1,n) = tgv_kedr(1,n) * tgv_dimen(1)
    end do
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
            '"<greek>e</greek>(E<sub>k</sub>) (Num+Phys KEDR)"',/, &
            '"<greek>e</greek>(<greek>z</greek>) (Physical KEDR)"',/, &
            '"<greek>e</greek><sub>1</sub> (Deviatoric Strain Rate)"',/, &
            '"<greek>e</greek><sub>3</sub> (Pressure Dialation)"',/, &
            '"E<sub>k</sub> (Kinetic Energy)"',/, &
            '"<greek>z</greek> (Enstrophy)"',/, &
            '"<greek>D</greek>t*"',/, &
            'ZONE T="P',i0,' - ',i0,'<sup>3</sup> ',a,' Grid"',/, &
            'AUXDATA order="',i0,'"',/, &
            'AUXDATA DOF="',i0,'<sup>3</sup>"',/, &
            'AUXDATA geom="',a,'"')
  2 format (100(es17.9))
 !2 format (*(es17.9))
  !
end subroutine compute_KEDR
!
!###############################################################################
!
subroutine compute_dEkdt_central_diff_O2(first_output,final_output, &
                                         tgv_idx,ek,dt,dekdt,nend)
  !
  !.. Formal Arugments ..
  logical(lk), intent(in) :: first_output
  logical(lk), intent(in) :: final_output
  integer,     intent(in) :: tgv_idx
  real(wp),    intent(in) :: ek(-1:)
  real(wp),    intent(in) :: dt(-1:)
  real(wp), intent(inout) :: dekdt(-1:)
  integer,  intent(inout) :: nend
  !
  !.. Local Scalars ..
  integer  :: n,np,nbeg
  real(wp) :: kem,kec,kep,dtm,dtp
  !
continue
  !
  ! Compute the central difference of the KEDR to get its rate of change
  ! over time
  ! NOTE: The KEDR found by differencing the integrated kinetic energy is
  !       actually given as KEDR(Ek) = -d/dt (Ek) so the differencing
  !       equation below has been multiplied through by a -1.
  !
  np = ubound(ek,dim=1)
  !
  nbeg = merge(     1     ,  0   , first_output )
  nend = merge( tgv_idx-1 , np-1 , final_output )
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) PRESENT(tgv_kedr)
  do n = nbeg,nend
    !
    kem = ek(n-1)
    kec = ek(n)
    kep = ek(n+1)
    !
    dtm = dt(n)
    dtp = dt(n+1)
    !
    dekdt(n) =  ( (dtm/dtp)*(kec-kep) + (dtp/dtm)*(kem-kec) ) / (dtm+dtp)
    !
  end do
  !$ACC UPDATE SELF(tgv_kedr)
  !
  ! Do a one sided difference on the initial conditions if this is the
  ! first output interval
  !
  if (first_output) then
    dekdt(0) = (ek(0) - ek(1)) / dt(1)
  end if
  !
  ! If this is the last output interval, do a one sided difference
  ! on the final solution
  !
  if (final_output) then
   !n = min(nend+1,np)
    n = nend+1 ! nend should be tgv_idx-1 which should be <= np-1
    dekdt(n) = (ek(n-1) - ek(n)) / dt(n)
  end if
  !
end subroutine compute_dEkdt_central_diff_O2
!
!###############################################################################
!
subroutine compute_dEkdt_backward_diff_O1(first_output,final_output, &
                                         tgv_idx,ek,dt,dekdt,nend)
  !
  !.. Formal Arugments ..
  logical(lk), intent(in) :: first_output
  logical(lk), intent(in) :: final_output
  integer,     intent(in) :: tgv_idx
  real(wp),    intent(in) :: ek(-1:)
  real(wp),    intent(in) :: dt(-1:)
  real(wp), intent(inout) :: dekdt(-1:)
  integer,  intent(inout) :: nend
  !
  !.. Local Scalars ..
  integer  :: n,np,nbeg
  real(wp) :: kem,kec,dtm
  !
continue
  !
  ! Compute the central difference of the KEDR to get its rate of change
  ! over time
  ! NOTE: The KEDR found by differencing the integrated kinetic energy is
  !       actually given as KEDR(Ek) = -d/dt (Ek) so the differencing
  !       equation below has been multiplied through by a -1.
  !
  np = ubound(ek,dim=1)
  !
  nbeg = merge(     1     ,  0   , first_output )
  nend = merge( tgv_idx-1 , np-1 , final_output )
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) PRESENT(tgv_kedr)
  do n = nbeg,nend
    !
    kem = ek(n-1)
    kec = ek(n)
    !
    dtm = dt(n)
    !
    dekdt(n) =  (kem - kec) / dtm
    !
  end do
  !$ACC UPDATE SELF(tgv_kedr)
  !
  ! Do a one sided difference on the initial conditions if this is the
  ! first output interval
  !
  if (first_output) then
    dekdt(0) = (ek(0) - ek(1)) / dt(1)
  end if
  !
  ! If this is the last output interval, do a one sided difference
  ! on the final solution
  !
  if (final_output) then
   !n = min(nend+1,np)
    n = nend+1 ! nend should be tgv_idx-1 which should be <= np-1
    dekdt(n) = (ek(n-1) - ek(n)) / dt(n)
  end if
  !
end subroutine compute_dEkdt_backward_diff_O1
!
!###############################################################################
!
subroutine compute_dEkdt_backward_diff_O2(first_output,final_output, &
                                         tgv_idx,ek,dt,dekdt,nend)
  !
  !.. Formal Arugments ..
  logical(lk), intent(in) :: first_output
  logical(lk), intent(in) :: final_output
  integer,     intent(in) :: tgv_idx
  real(wp),    intent(in) :: ek(-1:)
  real(wp),    intent(in) :: dt(-1:)
  real(wp), intent(inout) :: dekdt(-1:)
  integer,  intent(inout) :: nend
  !
  !.. Local Scalars ..
  integer  :: n,np,nbeg
  real(wp) :: kem,kec,kep,dtm,dtp
  !
continue
  !
  ! Compute the central difference of the KEDR to get its rate of change
  ! over time
  ! NOTE: The KEDR found by differencing the integrated kinetic energy is
  !       actually given as KEDR(Ek) = -d/dt (Ek) so the differencing
  !       equation below has been multiplied through by a -1.
  !
  np = ubound(ek,dim=1)
  !
  nbeg = merge(     1     ,  0   , first_output )
  nend = merge( tgv_idx-1 , np-1 , final_output )
  !
  !$ACC PARALLEL LOOP VECTOR_LENGTH(32) PRESENT(tgv_kedr)
  do n = nbeg,nend
    !
    kem = ek(n-1)
    kec = ek(n)
    kep = ek(n+1)
    !
    dtm = dt(n)
    dtp = dt(n+1)
    !
    dekdt(n) =  ( (dtm/dtp)*(kec-kep) + (dtp/dtm)*(kem-kec) ) / (dtm+dtp)
    !
  end do
  !$ACC UPDATE SELF(tgv_kedr)
  !
  ! Do a one sided difference on the initial conditions if this is the
  ! first output interval
  !
  if (first_output) then
    dekdt(0) = (ek(0) - ek(1)) / dt(1)
  end if
  !
  ! If this is the last output interval, do a one sided difference
  ! on the final solution
  !
  if (final_output) then
   !n = min(nend+1,np)
    n = nend+1 ! nend should be tgv_idx-1 which should be <= np-1
    dekdt(n) = (ek(n-1) - ek(n)) / dt(n)
  end if
  !
end subroutine compute_dEkdt_backward_diff_O2
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
end module generic_mod
