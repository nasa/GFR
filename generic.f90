module generic_mod
  !
  use module_kind_types
  use eqn_idx, only : nec,nmb,nme,nmx,nmy,nmz,nee,nq,ntk,ntl
  !
  implicit none
  !
  private
  !
  public :: init_generic_mod
  public :: rgeterr
  public :: getmaxent
  public :: compute_final_error
  public :: compute_TaylorGreen_KEDR
  public :: generic_memory_usage
  !
  character(len=100), save :: iter_out_fmt
  !
  real(wp), save :: Density_Maximum_Residual = -huge(zero)
  logical(lk),  save :: vortex_has_collapsed = fals
  integer,  save :: final_countdown = 0
  !
  ! TGV_KEDR : array used for the Taylor-Green vortex problem that
  !            contains the temporal evolution of various integrated
  !            quantities relating to the kinetic energy
  !
  integer,  save                              :: iotg
  real(wp), save, allocatable, dimension(:,:) :: tgv_kedr
  real(wp), save, allocatable, dimension(:,:) :: sum_kedr
  real(wp), save,              dimension(1:7) :: tgv_dimen = zero
  !
contains
!
!###############################################################################
!
subroutine compute_final_error()
  !
  ! Compute the final error. The error is
  ! stored in the arrays ultwo and vltwo as
  !
  ! ultwo(1)      = density
  ! ultwo(nmb:nme) = directional momentums
  ! ultwo(nq)     = energy
  ! ultwo(nq+1)   = entropy
  !
  ! vltwo(1)      = density
  ! vltwo(nmb:nme) = directional velocities
  ! vltwo(nq)     = pressure
  ! vltwo(nq+1)   = magnitude of velocity
  !
  !.. Use Statements ..
  use order_mod, only : maxEP
  use geovar,  only : nr,ncell,cell
  use ovar,    only : vlinf,vlone,vltwo,vnloc,vkloc,total_volume
  use ovar,    only : ulinf,ulone,ultwo,unloc,ukloc,itestcase
  use parallel_mod, only : cell_map
  use interpolation_mod, only : outerp
  !
  !.. Local Scalars ..
  integer  :: this_cell,n1,n2,np,k,m
  integer  :: this_geom,this_order
  real(wp) :: sn,se
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxEP) :: jac
  real(wp), dimension(1:maxEP) :: weights
  real(wp), dimension(1:nq,1:maxEP) :: uold
  real(wp), dimension(1:nq,1:maxEP) :: unew
  real(wp), dimension(1:nr,1:maxEP) :: xyz_err
  real(wp), dimension(1:nq+1,1:2) :: ett
  real(wp), dimension(1:nq+1,1:2,1:5) :: error
  real(wp), dimension(1:nq+1,1:2,1:5,1:ncpu) :: all_error
  !
  ! error(1:nq+1,:,:) = error value for each flow equation + extra variable
  ! error(:,1,:) = error for conserved variables [rho,rho(u,v,..),rhoe,entropy]
  ! error(:,2,:) = error for primitive variables [rho,u,v,..,pressure,vel mag]
  ! error(:,:,1) = L1 norm of total error for all solution points
  ! error(:,:,2) = L2 norm of total error for all solution points
  ! error(:,:,3) = maximum value of error for all solution points
  ! error(:,:,4) = ID of grid cell containing maximum value of error
  ! error(:,:,5) = ID of solution point with maximum value of error
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_final_error"
  !
continue
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
                              jac(1:np),weights(1:np),xyz_err(1:nr,1:np))
    !
    ! Now loop through the solution points of this cell to compute the
    ! error between the solutions of the current and previous time steps
    !
    do k = 1,np
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
      ett(nq+1,1) = sn - se
      !
      ! Error for primitive variables at the current solution point
      !
      ett(1:nq,2) = abs( usp2v_sp(unew(:,k)) - usp2v_sp(uold(:,k)) )
      ett(nq+1,2) = rnorm2( ett(nmb:nme,2) )
      !
      ! Now add the error for all variables at the current solution
      ! point to the running integral sum for the L1 and L2 error norms
      !
      error(:,:,1) = error(:,:,1) + jac(k)*weights(k)*ett(:,:)
      error(:,:,2) = error(:,:,2) + jac(k)*weights(k)*ett(:,:)*ett(:,:)
      !
      ! If any of these errors are greater than the previous maximum
      ! values, save the error value and record this location
      !
      where (abs(ett(:,:)) > error(:,:,3))
        error(:,:,3) = abs(ett(:,:))
        error(:,:,4) = real(cell_map(mypnum+1)%loc_to_glb(this_cell),kind=wp)
        error(:,:,5) = real(k,kind=wp)
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
  ! Find the max infinity norm for each variable
  !
  ulinf(:) = maxval(all_error(:,1,3,:),dim=2)
  vlinf(:) = maxval(all_error(:,2,3,:),dim=2)
  !
  ! Finish computing the L1 and L2 norms for each variable
  !
  ulone(:) =       sum(all_error(:,1,1,:),dim=2) / total_volume
  ultwo(:) = sqrt( sum(all_error(:,1,2,:),dim=2) / total_volume )
  vlone(:) =       sum(all_error(:,2,1,:),dim=2) / total_volume
  vltwo(:) = sqrt( sum(all_error(:,2,2,:),dim=2) / total_volume )
  !
  ! Get the location of the max infinity norm for each variable
  !
  do m = 1,nq+1
    unloc(m) = nint( all_error(m,1,4,maxloc(all_error(m,1,3,:),dim=1)) )
    ukloc(m) = nint( all_error(m,1,5,maxloc(all_error(m,1,3,:),dim=1)) )
    vnloc(m) = nint( all_error(m,2,4,maxloc(all_error(m,2,3,:),dim=1)) )
    vkloc(m) = nint( all_error(m,2,5,maxloc(all_error(m,2,3,:),dim=1)) )
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_final_error
!
!###############################################################################
!
subroutine rgeterr(walltime_expiring)
  !
  !.. Use Statements ..
  use, intrinsic :: ieee_exceptions, only : ieee_usual,ieee_get_flag
  !
  use order_mod,  only : maxEP
  use geovar,  only : nr,ncell,cell
  use geovar,  only : global_pinc_ptr
  use ovar,    only : rlinf,rlone,rltwo,nloc,kloc,total_volume,gam,bc_in
  use ovar,    only : itcur,d_t,time,rl2mx,minimum_timesteps,time_ref
  use ovar,    only : convergence_reached,cfl,itestcase
  use ovar,    only : convergence_order_max_res,convergence_order_abs_res
  use ovar,    only : this_is_final_timestep
  use ovar,    only : iter_out_interval
  use ovar,    only : Period_Timesteps
  use parallel_mod, only : cell_map
  use interpolation_mod, only : outerp
  !
  use pbs_mod, only : pbs_remaining_walltime
  !
  !.. Formal Arguments ..
  logical(lk), intent(inout) :: walltime_expiring
  !
  !.. Local Scalars ..
  integer  :: this_cell,n1,n2,np,m,l,k
  integer  :: hrs_left,min_left,sec_left
  integer  :: this_geom,this_order
  integer  :: nentmax,kentmax
  real(wp) :: xentmax,yentmax
  real(wp) :: entmax,entL2
  real(wp) :: rho,pres,sn,se,delx,dely
  real(wp) :: wall_time_left
  real(wp) :: orders_converged,output_L2val
  !
  character(len=300) :: overflow_error
  character(len=300) :: divide_by_zero_error
  character(len=300) :: invalid_error
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxEP) :: jac
  real(wp), dimension(1:maxEP) :: weights
  real(wp), dimension(1:nr,1:maxEP) :: xyz_err
  real(wp), dimension(1:nq,1:maxEP) :: uold
  real(wp), dimension(1:nq,1:maxEP) :: unew
  real(wp), dimension(1:2) :: ett
  real(wp), dimension(1:2,1:7) :: error
  real(wp), dimension(1:2,1:7,1:ncpu) :: all_error
  !
  integer,      dimension(1:size(ieee_usual)) :: bad_number_found
  logical(ldk), dimension(1:size(ieee_usual)) :: signaling_exceptions
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "rgeterr"
  !
  real(wp), parameter :: max_sec = real( 1000*3600-1 , kind=wp )
  !
  ! error(:,1) = L1 norm of total error for all solution points
  ! error(:,2) = L2 norm of total error for all solution points
  ! error(:,3) = maximum value of error for all solution points
  ! error(:,4) = ID of grid cell containing maximum value of error
  ! error(:,5) = ID of solution point with maximum value of error
  ! error(:,6) = x-coordinate of solution point with maximum value of error
  ! error(:,7) = y-coordinate of solution point with maximum value of error
  ! error(1,:) = error values for density
  ! error(2,:) = error values for entropy
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  overflow_error = ""
  divide_by_zero_error = ""
  invalid_error = ""
  !
  error(:,:) = zero  ! intialize everything to zero
  error(:,3) = -eps6 ! make this slightly negative in case all error is zero
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
                              jac(1:np),weights(1:np),xyz_err(1:nr,1:np))
    !
    ! Now loop through the solution points of this cell to compute the
    ! error between the solutions of the current and previous time steps
    !
    do k = 1,np
      !
      ! Get the change in density for this solution point
      !
      if (itestcase == Laminar_BndLyr) then
        ett(1) = abs( unew(nmb,k) - uold(nmb,k) )
      else
        ett(1) = abs( unew(nec,k) - uold(nec,k) )
      end if
      !
      ! Compute the entropy at this solution point
      !
      rho = unew(nec,k)
      pres = pressure_cv_sp( unew(:,k) )
      !
      if (any([rho,pres] <= zero)) write (iout,111) this_cell,k,rho,pres
      111 format ("cell = ",i0,"; SP = ",i0,"; ",*(1x,es14.6))
      !
      if (itestcase == Inviscid_Gaussian_Bump) then
        sn = entropy_cv_sp( unew(:,k) , log_opt=fals , use_pref_nd=true )
        se = one
      else
        sn = entropy_cv_sp( unew(:,k) )
        se = entropy_cv_sp( bc_in(0)%cv(:) )
      end if
      !
      ett(2) = sn - se
      !
      ! Get the norms of density and entropy
      !
      error(:,1) = error(:,1) + jac(k)*weights(k)*ett(:)
      error(:,2) = error(:,2) + jac(k)*weights(k)*ett(:)*ett(:)
      !
      where (abs(ett(:)) > error(:,3))
        error(:,3) = abs(ett)
        error(:,4) = real(cell_map(mypnum+1)%loc_to_glb(this_cell),kind=wp)
        error(:,5) = real(k,kind=wp)
        error(:,6) = xyz_err(1,k)
        error(:,7) = xyz_err(2,k)
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
    all_error(:,:,1) = error(:,:)
  end if
  !
  rlinf = maxval(all_error(1,3,:))
  nloc  = nint( all_error(1,4,maxloc(all_error(1,3,:),dim=1)) )
  kloc  = nint( all_error(1,5,maxloc(all_error(1,3,:),dim=1)) )
  !
  entmax  = maxval(all_error(2,3,:))
  nentmax = nint( all_error(2,4,maxloc(all_error(2,3,:),dim=1)) )
  kentmax = nint( all_error(2,5,maxloc(all_error(2,3,:),dim=1)) )
  xentmax =       all_error(2,6,maxloc(all_error(2,3,:),dim=1))
  yentmax =       all_error(2,7,maxloc(all_error(2,3,:),dim=1))
  !
  rlone =       sum(all_error(1,1,:)) / total_volume
  rltwo = sqrt( sum(all_error(1,2,:)) / total_volume )
  entL2 = sqrt( sum(all_error(2,2,:)) / total_volume )
  !
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
      ! that have been completed since the last call to rgeterr
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
  if (mypnum == 0) then
    !
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
    write (iout,iter_out_fmt) itcur,d_t*time_ref,time*time_ref, &
                              cfl,rltwo,output_L2val,nloc, &
                              global_pinc_ptr(nloc)+kloc, &
                              xentmax,yentmax, &
                              hrs_left,min_left,sec_left
    !
#else
    !
    write (iout,iter_out_fmt) itcur,d_t*time_ref,time*time_ref, &
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
    if (mypnum == 0) write (iout,6)
    convergence_reached = true
  end if
  !
  ! Check to make sure that no floating-point exceptions are signaling
  !
  call ieee_get_flag(ieee_usual,signaling_exceptions)
  !
  bad_number_found = merge(1,0,signaling_exceptions)
  !
  call mpi_allreduce(MPI_IN_PLACE,bad_number_found, &
                     size(bad_number_found,kind=int_mpi), &
                     mpi_inttyp,MPI_MAX,MPI_COMM_WORLD,mpierr)
  !
  if (sum(bad_number_found) > 0) then
    if (bad_number_found(1) > 0) write (overflow_error,8)
    if (bad_number_found(2) > 0) write (divide_by_zero_error,9)
    if (bad_number_found(3) > 0) write (invalid_error,10)
    write (error_message,7) trim(overflow_error), &
                            trim(divide_by_zero_error), &
                            trim(invalid_error)
    call stop_recon(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (i9,4(1x,es11.4),2i8,2es9.1,1x,i3,2(":",i2.2))
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
  7 format ("A SIGNALING FLOATING POINT EXCEPTION WAS ", &
            "DETECTED WHILE COMPUTING THE RESIDUAL NORM!!!", &
            3(2x,a,:))
  !
  8  format ("OVERFLOW EXCEPTION DETECTED!!!")
  9  format ("DIVIDE-BY-ZERO EXCEPTION DETECTED!!!")
  10 format ("INVALID NUMBER (NaN) EXCEPTION DETECTED!!!")
  !
end subroutine rgeterr
!
!###############################################################################
!
subroutine init_generic_mod
  !
  !.. Use Statements ..
  use geovar, only : n_global_cell,n_global_solpts
  use ovar,   only : itcur,num_ts,rl2mx,itestcase
  use ovar,   only : ulinf,ulone,ultwo,unloc,ukloc
  use ovar,   only : vlinf,vlone,vltwo,vnloc,vkloc
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
  ! Allocate the error arrays
  !
  if (.not. allocated(unloc)) then
    allocate ( unloc(1:nq+1) , source=0 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"unloc",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
  end if
  if (.not. allocated(ukloc)) then
    allocate ( ukloc(1:nq+1) , source=0 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ukloc",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
  end if
  if (.not. allocated(ulinf)) then
    allocate ( ulinf(1:nq+1) , source=zero , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ulinf",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
  end if
  if (.not. allocated(ulone)) then
    allocate ( ulone(1:nq+1) , source=zero , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ulone",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
  end if
  if (.not. allocated(ultwo)) then
    allocate ( ultwo(1:nq+1) , source=zero , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ultwo",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
  end if
  if (.not. allocated(vnloc)) then
    allocate ( vnloc(1:nq+1) , source=0 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vnloc",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
  end if
  if (.not. allocated(vkloc)) then
    allocate ( vkloc(1:nq+1) , source=0 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vkloc",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
  end if
  if (.not. allocated(vlinf)) then
    allocate ( vlinf(1:nq+1) , source=zero , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vlinf",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
  end if
  if (.not. allocated(vlone)) then
    allocate ( vlone(1:nq+1) , source=zero , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vlone",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
  end if
  if (.not. allocated(vltwo)) then
    allocate ( vltwo(1:nq+1) , source=zero , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"vltwo",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
  end if
  !
  rl2mx = -huge(zero)
  !
  write (char_num,1) itcur + num_ts
  final_iter_digits = max(4,len_trim(adjustl(char_num)))
  !
  write (char_num,1) n_global_cell
  max_cell_digits = max(4,len_trim(adjustl(char_num)))
  !
  write (char_num,1) n_global_solpts
  max_solpt_digits = max(5,len_trim(adjustl(char_num)))
  !
  write (iter_out_fmt,100) final_iter_digits,max_cell_digits,max_solpt_digits
  !
  ! Add on the remaining simulation time to the output if
  ! running in the PBS environment
  !
#ifdef PBS_ENV
  iter_out_fmt = iter_out_fmt(1:len_trim(iter_out_fmt)) // ',1x,i3,2(":",i2.2)'
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
    write (iter_header,203) final_iter_digits-3, &
                            max_cell_digits-2, &
                            max_solpt_digits-3
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
              "  Wall Time Remaining = ",i3,2(":",i2.2))
  100 format ("(1x,i",i0,",5es12.4,2x,i",i0,",2x,i",i0,",2es9.1")
  200 format ("(/,",i0,"x,'iter',5x,'d_t',8x,'Time',6x,'L2-Density',2x,", &
              "'L2-Entropy',",i0,"x,'Cell',",i0,"x,'Solpt',5x,'x',8x,'y'")
  201 format ("(/,",i0,"x,'iter',5x,'d_t',8x,'Time',6x,'L2-X-Mmntm',2x,", &
              "'OrdersDrop',2x,'Max-DensL2',",i0,"x,'Cell',",i0, &
              "x,'Solpt',5x,'x',8x,'y'")
  202 format ("(/,",i0,"x,'iter',5x,'d_t',8x,'Time',10x,'CFL',5x,", &
              "'L2-Density',2x,'Max-DensL2',",i0,"x,'Cell',",i0, &
              "x,'Solpt',5x,'x',8x,'y'")
 !202 format ("(/,",i0,"x,'iter',5x,'d_t',8x,'Time',6x,'L2-Density',2x,", &
 !            "'OrdersDrop',2x,'Max-DensL2',",i0,"x,'Cell',",i0, &
 !            "x,'Solpt',5x,'x',8x,'y'")
  203 format ("(/,",i0,"x,'iter',5x,'dt*',8x,'Time*',5x,'L2-Density',2x,", &
              "'OrdersDrop',2x,'Max-DensL2',",i0,"x,'Cell',",i0, &
              "x,'Solpt',5x,'x',8x,'y'")
  !
end subroutine init_generic_mod
!
!###############################################################################
!
subroutine getmaxent
  !
  !... get max entropy
  !
  !.. Use Statements ..
  use order_mod,         only : maxEP,e_order
  use geovar,            only : nr,ncell,cell,xyz
  use ovar,              only : gam,bc_in,total_volume,itestcase
  use metrics_mod,       only : jacobian
  use flowvar,           only : usp
  use parallel_mod,      only : cell_map
  use quadrature_mod,    only : std_elem
  use interpolation_mod, only : outerp
  !
  !.. Local Scalars ..
  integer  :: this_cell,n1,n2,np,m,l,k
  integer  :: this_geom,this_order
  integer  :: nentmax,kentmax
  real(wp) :: xentmax,yentmax
  real(wp) :: entmax,entL2
  real(wp) :: ett,sn,se
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxEP) :: jac
  real(wp), dimension(1:maxEP) :: weights
  real(wp), dimension(1:nr,1:maxEP) :: xyz_err
  real(wp), dimension(1:nq,1:maxEP) :: unew
  real(wp), dimension(1:nq,1:maxEP) :: uold
  real(wp), dimension(1:6) :: error
  real(wp), dimension(1:6,1:ncpu) :: all_error
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "getmaxent"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  error(:) = zero  ! intialize everything to zero
  error(2) = -eps6 ! make this slightly negative in case all error is zero
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
    ! Interpolate the new and old solution to the error points
    !
    call collect_error_arrays(this_cell,unew(1:nq,1:np),uold(1:nq,1:np), &
                              jac(1:np),weights(1:np),xyz_err(1:nr,1:np))
    !
   !if (allocated(outerp(this_geom,this_order)%error)) then
   !  np = size(outerp(this_geom,this_order)%error%mat,dim=1)
   !  do l = 1,nr
   !    xyz_err(l,1:np) = outerp(this_geom,this_order)% &
   !                                             error% &
   !                      mv_mult( xyz(l,n1:n2) )
   !  end do
   !  do m = 1,nq
   !    unew(m,1:np) = outerp(this_geom,this_order)% &
   !                                          error% &
   !                   mv_mult( usp(m,n1:n2) )
   !  end do
   !  jac(1:np) = outerp(this_geom,this_order)% &
   !                                     error% &
   !              mv_mult( jacobian(n1:n2) )
   !  weights(1:np) = outerp(this_geom,this_order)%err_wts
   !else
   !  xyz_err(1:nr,1:np) = xyz(1:nr,n1:n2)
   !  unew(1:nq,1:np) = usp(1:nq,n1:n2)
   !  jac(1:np) = jacobian(n1:n2)
   !  weights(1:np) = std_elem(this_geom,this_order)%wts
   !end if
    !     
    ! Now loop through the solution points of this cell to compute the
    ! error between the solutions of the current and previous time steps
    !
    do k = 1,np
      !
      if (itestcase == Inviscid_Gaussian_Bump) then
        sn = entropy_cv_sp( unew(:,k) , log_opt=fals , use_pref_nd=true )
        se = one
      else
        sn = entropy_cv_sp( unew(:,k) )
        se = entropy_cv_sp( bc_in(0)%cv(:) )
      end if
      !
      ett = sn - se
      !
      ! L2 norm of entropy
      !
      error(1) = error(1) + jac(k)*weights(k)*ett*ett
      !
      if (abs(ett) > error(2)) then
        error(2) = abs(ett)
        error(3) = real(cell_map(mypnum+1)%loc_to_glb(this_cell),kind=wp)
        error(4) = real(k,kind=wp)
        error(5) = xyz_err(1,k)
        error(6) = xyz_err(2,k)
      end if
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
    all_error(:,1) = error(:)
  end if
  !
  entmax  = maxval(all_error(2,:))
  nentmax = nint( all_error(3,maxloc(all_error(2,:),dim=1)) )
  kentmax = nint( all_error(4,maxloc(all_error(2,:),dim=1)) )
  xentmax =       all_error(5,maxloc(all_error(2,:),dim=1))
  yentmax =       all_error(6,maxloc(all_error(2,:),dim=1))
  !
  entL2 = sqrt( sum(all_error(1,:)) / total_volume )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine getmaxent
!
!###############################################################################
!
subroutine collect_error_arrays(this_cell,unew,uold,jac,weights,xyz_err)
  !
  !.. Use Statements ..
  use order_mod,         only : maxSP
  use geovar,            only : nr,cell,xyz
  use ovar,              only : mms_opt,itestcase
  use flowvar,           only : usp,uoldsp
  use interpolation_mod, only : outerp
  use quadrature_mod,    only : std_elem
  use metrics_mod,       only : jacobian
  use mms_mod,           only : mms_solution_at_xyz
  !
  !.. Formal Arguments ..
  integer,                     intent(in) :: this_cell
  real(wp), dimension(:,:), intent(inout) :: unew
  real(wp), dimension(:,:), intent(inout) :: uold
  real(wp), dimension(:,:), intent(inout) :: xyz_err
  real(wp), dimension(:),   intent(inout) :: jac
  real(wp), dimension(:),   intent(inout) :: weights
  !
  !.. Local Scalars ..
  integer :: m,l,k,n1,n2,np
  integer :: this_geom,this_order
  logical(lk) :: use_error_matrix
  !
  !.. Local Arrays ..
  real(wp) :: xyz_offset(1:nr,1:maxSP)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "collect_error_arrays"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  this_geom = cell(this_cell)%geom
  this_order = cell(this_cell)%order
  !
  n1 = cell(this_cell)%beg_sp
  n2 = cell(this_cell)%end_sp
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
    np = size(outerp(this_geom,this_order)%error%mat,dim=1)
  else
    np = n2-n1+1
  end if
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Get the current solution at the error points
  !
  if (use_error_matrix) then
    !
    do m = 1,nq
      unew(m,1:np) = outerp(this_geom,this_order)% &
                                            error% &
                     mv_mult( usp(m,n1:n2) )
    end do
    !
  else
    !
    unew(1:nq,1:np) = usp(1:nq,n1:n2)
    !
  end if
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Get the local jacobians at the error points
  !
  if (use_error_matrix) then
    !
    jac(1:np) = outerp(this_geom,this_order)% &
                                       error% &
                mv_mult( jacobian(n1:n2) )
    !
  else
    !
    jac(1:np) = jacobian(n1:n2)
    !
  end if
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Get the quadrature weights for the error points
  !
  if (use_error_matrix) then
    !
    weights(1:np) = outerp(this_geom,this_order)%err_wts
    !
  else
    !
    weights(1:np) = std_elem(this_geom,this_order)%wts
    !
  end if
  !
  ! #######
  ! ## 4 ##
  ! #######
  !
  ! Get the coordinates of the error points
  !
  if (use_error_matrix) then
    !
    do l = 1,nr
      xyz_err(l,1:np) = outerp(this_geom,this_order)% &
                                               error% &
                        mv_mult( xyz(l,n1:n2) )
    end do
    !
  else
    !
    xyz_err(1:nr,1:np) = xyz(1:nr,n1:n2)
    !
  end if
  !
  ! #######
  ! ## 5 ##
  ! #######
  !
  ! Get either (a) the exact analytical solution or (b) solution from
  ! the previous time step at the error points.
  !
  if (mms_opt /= 0) then
    !
    ! Evaluate the exact analytical MMS solution at the error points
    !
    do k = 1,np
      uold(1:nq,k) = mms_solution_at_xyz( xyz_err(1:nr,k), &
                                          return_nondim=true, &
                                          return_cv=true )
    end do
    !
  else if (any(itestcase == Transport_Problems)) then
    !
    ! Evaluate the exact analytical solution to the vortex transport problem
    !
    xyz_offset(1:nr,1:np) = vortex_offset( xyz_err(1:nr,1:np) )
    !
    do k = 1,np
      uold(1:nq,k) = vortex_solution( xyz_offset(1:nr,k) )
    end do
    !
  else
    !
    if (use_error_matrix) then
      !
      ! Get the solution from the previous time step at the error points
      !
      do m = 1,nq
        uold(m,1:np) = outerp(this_geom,this_order)% &
                                              error% &
                       mv_mult( uoldsp(m,n1:n2) )
      end do
      !
    else
      !
      ! Use the solution from the previous time step
      !
      uold(1:nq,1:np) = uoldsp(1:nq,n1:n2)
      !
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine collect_error_arrays
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
  use interpolation_mod, only : outerp
  use quadrature_mod, only : std_elem
  use metrics_mod, only : jacobian
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: iter
  real(wp), intent(in) :: dtmin
  !
  !.. Local Scalars ..
  integer :: k,l,m,n,nc,np,n1,n2,tgv_idx
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
#ifdef DEBUG_ON
 !if (mypnum == 0) then
 !  if (first_output) then
 !    write (iout,13) "FIRST",iter,tgv_idx
 !  end if
 !  if (final_output) then
 !    write (iout,13) "FINAL",iter,tgv_idx
 !  end if
 !end if
 !13 format (/," THIS IS THE ",a," TGV-KEDR OUTPUT!",/, &
 !             "   niter=",i0,", tgv_idx=",i0,/)
#endif
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
    allocate ( tgv_kedr(1:7,-1:nmax) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"tgv_kedr",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !
    allocate ( sum_kedr(1:4,0:nmax) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"sum_kedr",1,__LINE__,__FILE__,ierr, &
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
    ! Set tgv_kedr(ndt,-1) to time the first time through
    ! time will have the value zero for a new simulation or the
    ! restart time if restarting a simulation
    !
    tgv_kedr(ndt,-1) = time * time_ref
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
  tgv_kedr(ndt,tgv_idx) = dtmin * time_ref
  !
  ! Loop through the cells to compute the integrated
  ! quantities local to each processor
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
    ! Collect the primitive variables, velocity gradients, jacobian, and
    ! quadrature weights at the error quadrature points for the current cell
    !
    if (allocated(outerp(this_geom,this_order)%error)) then
      !
      ! The error quadrature points are different from the solution points so
      ! we must interpolate these values from the solution points to the
      ! error quadrature points
      !
      ! NOTE: The vsp array will temporarily hold the conservative variables
      do m = 1,nq
        vsp(m,1:np) = outerp(this_geom,this_order)% &
                                             error% &
                      mv_mult( usp(m,n1:n2) )
      end do
      !
      do m = 1,nq
        do l = 1,nr
          dcvdx(l,m,1:np) = outerp(this_geom,this_order)% &
                                                   error% &
                            mv_mult( dusp(l,m,n1:n2) )
        end do
      end do
      !
      wtdjac(1:np) = outerp(this_geom,this_order)% &
                                            error% &
                     mv_mult( jacobian(n1:n2) )
      !
      wtdjac(1:np) = wtdjac(1:np) * outerp(this_geom,this_order)%err_wts
      !
    else
      !
      ! The error quadrature points and solution points are co-located so
      ! just copy these values into the local arrays.
      ! NOTE: The vsp array will temporarily hold the conservative variables
      !
      vsp(1:nq,1:np) = usp(1:nq,n1:n2)
      dcvdx(1:nr,1:nq,1:np) = dusp(1:nr,1:nq,n1:n2)
      wtdjac(1:np) = jacobian(n1:n2) * std_elem(this_geom,this_order)%wts
      !
    end if
    !
    ! Convert the conservative variable derivatives
    ! to primitive variable derivatives
    ! NOTE: The vsp array currently contains the conservative variables
    !
    do k = 1,np
      dpvdx(1:nr,1:nq,k) = grad_cv_to_grad_pv_sp( vsp(1:nq,k) , &
                                                  dcvdx(1:nr,1:nq,k) )
    end do
    !
    ! Convert the conservative variables in vsp to primitive variables
    !
    vsp(1:nq,1:np) = usp2v( vsp(1:nq,1:np) )
    !
    ! Loop through the solution points of the current 
    !
    do k = 1,np
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
      sum_kedr(1,tgv_idx) = sum_kedr(1,tgv_idx) + &
                            wtdjac(k)*vsp(nec,k)*selfdot( vort )
      !
      ! Contribution to theoretical KEDR from the deviatoric strain rate tensor
      !
      sum_kedr(2,tgv_idx) = sum_kedr(2,tgv_idx) + &
                            wtdjac(k)*sum(strain*strain)
      !
      ! Contribution to theoretical KEDR from the pressure dialation term
      !
      sum_kedr(3,tgv_idx) = sum_kedr(3,tgv_idx) + &
                            wtdjac(k)*vsp(nee,k)*div_vel
      !
      ! Integrated kinetic energy (later, it will
      ! be differenced in time to get KEDR)
      !
      sum_kedr(4,tgv_idx) = sum_kedr(4,tgv_idx) + &
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
#ifdef DEBUG_ON
   !write (mypnum+20,*)
   !write (mypnum+20,'("  iter = ",i0)') iter
   !do n = lbound(sum_kedr,dim=2),ubound(sum_kedr,dim=2)
   !  write (mypnum+20,15) n,(sum_kedr(k,n), &
   !      k=lbound(sum_kedr,dim=1),ubound(sum_kedr,dim=1))
   !end do
   !15 format (2x,i3,*(2x,es14.6))
#endif
    !
    if (ncpu > 1) then
      !
      call mpi_allreduce(MPI_IN_PLACE,sum_kedr,size(sum_kedr,kind=int_mpi), &
                         mpi_flttyp,MPI_SUM,MPI_COMM_WORLD,mpierr)
      !
    end if
    !
    ! If this this the first output interval, change nmin to 0 so that the
    ! KEDR values for the initial conditions are also modified
    !
    if (first_output) nmin = 0
    !
#ifdef DEBUG_ON
   !write (mypnum+20,*)
   !write (mypnum+20,14) nmin,nmax
   !do n = lbound(sum_kedr,dim=2),ubound(sum_kedr,dim=2)
   !  write (mypnum+20,15) n,(sum_kedr(k,n), &
   !      k=lbound(sum_kedr,dim=1),ubound(sum_kedr,dim=1))
   !end do
   !14 format ("  nmin = ",i0,",  nmax = ",i0)
#endif
    !
    ! Finish computing the different KEDR
    !
    tgv_kedr(2:5,nmin:nmax) = sum_kedr(1:4,nmin:nmax) / total_volume
    !
    ! Copy tgv_kedr(2) into tgv_kedr(6)
    !
    tgv_kedr(6,nmin:nmax) = tgv_kedr(2,nmin:nmax)
    !
    ! Adjust the nondimensionalization for some of these KEDR values
    !
    do n = nmin,nmax
      tgv_kedr(2,n) = tgv_kedr(2,n) * tgv_dimen(2)
      tgv_kedr(3,n) = tgv_kedr(3,n) * tgv_dimen(3)
      tgv_kedr(4,n) = tgv_kedr(4,n) * tgv_dimen(4)
      tgv_kedr(5,n) = tgv_kedr(5,n) * tgv_dimen(5)
      tgv_kedr(6,n) = tgv_kedr(6,n) * tgv_dimen(6)
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
      kem = tgv_kedr(nke,n-1)
      kec = tgv_kedr(nke,n)
      kep = tgv_kedr(nke,n+1)
      !
      dtm = tgv_kedr(ndt,n)
      dtp = tgv_kedr(ndt,n+1)
      !
#ifdef DEBUG_ON
     !if (mypnum == 0) then
     !  write (50,16) n,nbeg,nend
     !  write (50,17) n-1,n,n+1,n,n-1,n+1,n
     !  write (50,18) kem,kec,kep,dtm,dtp
     !end if
     !16 format ("  n=",i0,", nbeg=",i0,", nend=",i0)
     !17 format (6x,"kedr(",i0")",9x,"kedr(",i0")",9x,"kedr(",i0")",9x, &
     !           "dt(",i0,"-",i0,")",9x,"dt(",i0,"-",i0,")")
     !18 format (*(2x,es14.6))
     !write (iout,12) mypnum,n,nbeg,nend,kem,kec,kep,dtm,dtp
     !12 format ("  cpu ",i0,": n=",i0,", nbeg=",i0,", nend=",i0,*(2x,es14.6))
#endif
      tgv_kedr(1,n) =  ( (dtm/dtp)*(kec-kep) + (dtp/dtm)*(kem-kec) ) / (dtm+dtp)
      !
    end do
    !
    ! Do a one sided difference on the initial conditions if this is the
    ! first output interval
    !
    if (first_output) then
      n = 0
      tgv_kedr(1,n) = (tgv_kedr(nke,n) - tgv_kedr(nke,n+1)) / tgv_kedr(ndt,n+1)
    end if
    !
    ! If this is the last output interval, do a one sided difference
    ! on the final solution
    !
    if (final_output) then
      n = tgv_idx
      tgv_kedr(1,n) = (tgv_kedr(nke,n-1) - tgv_kedr(nke,n)) / tgv_kedr(ndt,n)
    end if
    !
    ! Have the root processor output the Taylor-Green KEDR solutions
    !
    if (mypnum == 0) then
      !
      solution_time = tgv_kedr(ndt,-1)
      do n = 0,nend
        solution_time = solution_time + tgv_kedr(ndt,n)
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
        solution_time = solution_time + tgv_kedr(ndt,n)
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
      solution_time = sum( tgv_kedr(ndt,-1:nend) )
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
      ! be used for the next output interval. Set tgv_kedr(ndt,-1) to
      ! the last solution time that was output. Initialize all the other
      ! values to zero.
      !
#ifdef DEBUG_ON
     !write (iout,11) mypnum,solution_time
     !11 format (2x,"cpu ",i0,": solution_time = ",es23.15)
#endif
      tgv_kedr(:,-1) = [tgv_kedr(1:6,nmax-1),solution_time]
      tgv_kedr(:,0)  = tgv_kedr(:,nmax)
      tgv_kedr(:,1:nmax) = zero
      sum_kedr = zero
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
  2 format (*(es17.9))
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
end subroutine generic_memory_usage
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
