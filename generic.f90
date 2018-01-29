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
  public :: write_iter_stats
  public :: getmaxent
  public :: compute_residual_error
  public :: compute_TaylorGreen_KEDR
  public :: generic_memory_usage
  !
  character(len=100), save :: iter_out_fmt
  !
  real(wp),    save :: Density_Maximum_Residual = -huge(zero)
  logical(lk), save :: vortex_has_collapsed = fals
  integer,     save :: final_countdown = 0
 !!
 !real(wp), save :: rl2mx = -huge(zero)
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
contains
!
!###############################################################################
!
subroutine compute_residual_error()
  !
  !.. Use Statements ..
  use order_mod, only : maxEP
  use geovar,  only : nr,ncell,cell
  use ovar,    only : total_volume,itestcase
  use ovar,    only : uverr
  use parallel_mod, only : cell_map
  use interpolation_mod, only : outerp
  !
  !.. Local Scalars ..
  integer  :: this_cell,n1,n2,np,k,m,n
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
#ifdef DEBUG_ON
  call debug_timer(entering_procedure,pname)
#endif
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
      ett(nq+1,1) = sn - se
      !
      ! Error for primitive variables at the current solution point
      !
      ett(1:nq,2) = abs( vnew(:,k) - vold(:,k) )
      ett(nq+1,2) = norm2( ett(nmb:nme,2) )
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
   !uverr(:,1,n) =       sum(all_error(:,n,1,:),dim=2) / total_volume
   !! L2 norm for the current set of variables
   !uverr(:,2,n) = sqrt( sum(all_error(:,n,2,:),dim=2) / total_volume )
   !! Infinity norm for the current set of variables
   !uverr(:,3,n) = maxval(all_error(:,n,3,:),dim=2)
    !
    do m = 1,nq+1
      !
      ! L1 norm for the current variable
      uverr(m,1,n) =       sum(all_error(m,n,1,:)) / total_volume
      ! L2 norm for the current variable
      uverr(m,2,n) = sqrt( sum(all_error(m,n,2,:)) / total_volume )
      ! Infinity norm for the current variable
      uverr(m,3,n) = maxval(all_error(m,n,3,:))
      !
      ! Processor with max inf norm for the current variable
      k = maxloc(all_error(m,n,3,:),dim=1)
      ! ID of grid cell with max inf norm for the current variable
      uverr(m,4,n) = all_error(m,n,4,k)
      ! ID of sol point with max inf norm for the current variable
      uverr(m,5,n) = all_error(m,n,5,k)
      ! Coordinates of sol point with max inf norm for the current variable
      uverr(m,6:5+nr,n) = all_error(m,n,6:5+nr,k)
      !
    end do
    !
  end do
  !
#ifdef DEBUG_ON
  call debug_timer(leaving_procedure,pname)
#endif
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
subroutine write_iter_stats(walltime_expiring,exit_main_loop)
  !
  !.. Use Statements ..
  use geovar,  only : global_pinc_ptr
  use ovar,    only : itcur,d_t,cfl,time,minimum_timesteps,time_ref
  use ovar,    only : convergence_reached,itestcase,Period_Timesteps
  use ovar,    only : convergence_order_max_res,convergence_order_abs_res
  use ovar,    only : rl2mx,uverr,this_is_final_timestep,iter_out_interval
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
#ifdef DEBUG_ON
  call debug_timer(entering_procedure,pname)
#endif
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
#ifdef DEBUG_ON
  call debug_timer(leaving_procedure,pname)
#endif
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
              "  Wall Time Remaining = ",i3,2(":",i2.2))
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
  use interpolation_mod, only : outerp
  use quadrature_mod,    only : std_elem
  use metrics_mod,       only : metrics_dt
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
  integer :: m,l,k,n1,n2,np,nep
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
    ucurr(1:nq,1:np) = uavesp(1:nq,n1:n2)
  else
    ucurr(1:nq,1:np) = usp(1:nq,n1:n2)
  end if
  !
  if (use_error_matrix) then
    !
    do m = 1,nq
      unew(m,1:nep) = outerp(this_geom,this_order)% &
                                             error% &
                      mv_mult( ucurr(m,1:np) )
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
      vcurr(1:nq,1:np) = uavesp(nq+1:2*nq,n1:n2)
    else
      vcurr(1:nq,1:np) = usp2v( usp(1:nq,n1:n2) , swap_t_for_rho=true )
    end if
    !
    if (use_error_matrix) then
      !
      do m = 1,nq
        vnew(m,1:nep) = outerp(this_geom,this_order)% &
                                               error% &
                        mv_mult( vcurr(m,1:np) )
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
    wtdjac(1:nep) = outerp(this_geom,this_order)% &
                                           error% &
                    mv_mult( metrics_dt(this_cell)%jac(:) )
    !
    wtdjac(1:nep) = wtdjac(1:nep) * outerp(this_geom,this_order)%err_wts
    !
  else
    !
    wtdjac(1:nep) = metrics_dt(this_cell)%jac(:) * &
                    std_elem(this_geom,this_order)%wts
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
      xyz_err(l,1:nep) = outerp(this_geom,this_order)% &
                                                error% &
                         mv_mult( xyz(l,n1:n2) )
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
        uprev(1:nq,1:np) = uoldsp(1:nq,n1:n2)
      else
        a = (prev_ave_time-saved_ave_time) / (saved_ave_time-ave_start_time)
        b = one - a
        uprev(1:nq,1:np) = b * uavesp(1:nq,n1:n2) - a * usp(1:nq,n1:n2)
       !if (any(uprev(nec,1:np) <= zero)) then
       !  uprev(1:nq,1:np) = uoldsp(1:nq,n1:n2)
       !else if (any(pressure_cv(uprev(1:nq,1:np)) <= zero)) then
       !  uprev(1:nq,1:np) = uoldsp(1:nq,n1:n2)
       !end if
      end if
    else
      uprev(1:nq,1:np) = uoldsp(1:nq,n1:n2)
    end if
    !
    if (use_error_matrix) then
      !
      ! Get the solution from the previous time step at the error points
      !
      do m = 1,nq
        uold(m,1:nep) = outerp(this_geom,this_order)% &
                                               error% &
                        mv_mult( uprev(m,1:np) )
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
          vprev(1:nq,1:np) = usp2v( uoldsp(1:nq,n1:n2) , swap_t_for_rho=true )
        else
          ! a and b were already computed above, no need to do it again
          vprev(1:nq,1:np) = b * uavesp(nq+1:2*nq,n1:n2) - &
                             a * usp2v( usp(1:nq,n1:n2) , swap_t_for_rho=true )
        end if
      else
        vprev(1:nq,1:np) = usp2v( uoldsp(1:nq,n1:n2) , swap_t_for_rho=true )
      end if
      !
      if (use_error_matrix) then
        !
        ! Get the solution from the previous time step at the error points
        !
        do m = 1,nq
          vold(m,1:nep) = outerp(this_geom,this_order)% &
                                                 error% &
                          mv_mult( vprev(m,1:np) )
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
  use metrics_mod, only : metrics_dt
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
        vsp(m,1:nep) = outerp(this_geom,this_order)% &
                                              error% &
                       mv_mult( usp(m,n1:n2) )
      end do
      !
      do m = 1,nq
        do l = 1,nr
          dcvdx(l,m,1:nep) = outerp(this_geom,this_order)% &
                                                    error% &
                             mv_mult( dusp(l,m,n1:n2) )
        end do
      end do
      !
      wtdjac(1:nep) = outerp(this_geom,this_order)% &
                                             error% &
                      mv_mult( metrics_dt(this_cell)%jac(:) )
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
      vsp(1:nq,1:nep) = usp(1:nq,n1:n2)
      dcvdx(1:nr,1:nq,1:nep) = dusp(1:nr,1:nq,n1:n2)
      wtdjac(1:nep) = metrics_dt(this_cell)%jac * &
                      std_elem(this_geom,this_order)%wts
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
  2 format (100(es17.9))
 !2 format (*(es17.9))
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
include "Functions/usp2v.f90"
include "Functions/entropy.f90"
include "Functions/pressure.f90"
include "Functions/vortex.f90"
include "Functions/grad_cv_to_grad_pv.f90"
!
!###############################################################################
!
end module generic_mod
