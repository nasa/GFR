module generic_mod
  !
  use module_kind_types
  use eqn_idx, only : nec,nmb,nme,nee,nq
  !
  implicit none
  !
  private
  !
  public :: init_generic_mod
  public :: rgeterr
  public :: getmaxent
  public :: compute_final_error
  !
  character(len=100), save :: iter_out_fmt
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
  use order_mod, only : maxpts
  use geovar,  only : ncell,pts_in_cell_ptr
  use geovar,  only : cell_geom,cell_order
  use ovar,    only : vlinf,vlone,vltwo,vnloc,vkloc,total_volume
  use ovar,    only : ulinf,ulone,ultwo,unloc,ukloc
  use metrics_mod,  only : jacobian
  use flowvar, only : usp,uoldsp
  use parallel_mod, only : cell_map,cell_map_ptr
  use quadrature_mod, only : std_elem
  use interpolation_mod, only : outerp
  !
  !.. Local Scalars ..
  integer  :: n,n1,n2,np,k,m,nc1
  integer  :: this_geom,this_order
  real(wp) :: sn,se
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxpts) :: jac
  real(wp), dimension(1:maxpts) :: weights
  real(wp), dimension(1:nq,1:maxpts) :: uold
  real(wp), dimension(1:nq,1:maxpts) :: unew
  real(wp), dimension(1:nq+1,1:2) :: ett
  real(wp), dimension(1:nq+1,1:2,1:5) :: error
  real(wp), dimension(1:nq+1,1:2,1:5,1:ncpu) :: all_error
  !
  ! error(1:nq+1,:,:) = error value for each flow equation + extra variable
  ! error(:,1,:) = error for conserved variables [rho,rho(u,v,..),rhoe,entropy]
  ! error(:,2,:) = error for primitive variables [rho,u,v,..,pressure,mag of vel]
  ! error(:,:,1) = L1 norm of total error for all solution points
  ! error(:,:,2) = L2 norm of total error for all solution points
  ! error(:,:,3) = maximum value of error for all solution points
  ! error(:,:,4) = ID of grid cell containing maximum value of error
  ! error(:,:,5) = ID of solution point with maximum value of error
  !
continue
  !
  error(:,:,:) = zero  ! intialize everything to zero
  error(:,:,3) = -eps6 ! make this slightly negative in case all error is zero
  !
  nc1 = cell_map_ptr(mypnum+1)
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1
    n2 = pts_in_cell_ptr(n+1)
    np = n2-n1+1
    !
    this_geom = cell_geom(n)
    this_order = cell_order(n)
    !
    ! Interpolate the new and old solution to the error points
    !
    if (associated(outerp(this_geom,this_order)%error)) then
      do m = 1,nq
        unew(m,1:np) = matmul( outerp(this_geom,this_order)%error , usp(m,n1:n2) )
        uold(m,1:np) = matmul( outerp(this_geom,this_order)%error , uoldsp(m,n1:n2) )
      end do
      jac(1:np) = matmul( outerp(this_geom,this_order)%error , jacobian(n1:n2) )
      weights(1:np) = outerp(this_geom,this_order)%err_wts
    else
      unew(1:nq,1:np) = usp(1:nq,n1:n2)
      uold(1:nq,1:np) = uoldsp(1:nq,n1:n2)
      jac(1:np) = jacobian(n1:n2)
      weights(1:np) = std_elem(this_geom,this_order)%wts
    end if
    !
    ! Now loop through the solution points of this cell to compute the
    ! error between the solutions of the current and previous time steps
    !
    do k = 1,np
      !
      ! Error for conserved variables at the current solution point
      !
      sn = entropy_cv_sp( unew(:,k) )
      se = entropy_cv_sp( uold(:,k) )
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
        error(:,:,4) = real(cell_map(nc1+n),kind=wp)
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
    call mpi_allgather(error,    size(error),mpi_flttyp, &
                       all_error,size(error),mpi_flttyp,mpi_comm_world,mpierr)
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
end subroutine compute_final_error
!
!###############################################################################
!
subroutine rgeterr()
  !
  !.. Use Statements ..
  use order_mod,  only : maxpts
  use geovar,  only : nr,ncell,pts_in_cell_ptr,xyz
  use geovar,  only : cell_geom,cell_order
  use geovar,  only : global_pinc_ptr,global_xyz
  use ovar,    only : rlinf,rlone,rltwo,nloc,kloc,total_volume,gam,cv_in
  use ovar,    only : nentmax,kentmax,xentmax,yentmax,entmax,entL2,itestcase
  use ovar,    only : itcur,d_t,time,rl2mx
  use ovar,    only : max_wall_time,wall_time_left
  use ovar,    only : beg_wall_time,wall_time
  use ovar,    only : walltime_expiring,convergence_reached
  use ovar,    only : convergence_order_max_res,convergence_order_abs_res
  use metrics_mod,  only : jacobian
  use flowvar, only : usp,uoldsp
  use parallel_mod, only : cell_map,cell_map_ptr
  use quadrature_mod, only : std_elem
  use interpolation_mod, only : outerp
  !
  !.. Local Scalars ..
  integer  :: n,n1,n2,np,m,l,k,nc1
  integer  :: hrs_left,min_left,sec_left
  integer  :: this_geom,this_order
  real(wp) :: rho,pres,sn,se
  real(wp) :: next_time_est
  real(wp) :: wall_time_last
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxpts) :: jac
  real(wp), dimension(1:maxpts) :: weights
  real(wp), dimension(1:nr,1:maxpts) :: xyz_err
  real(wp), dimension(1:nq,1:maxpts) :: uold
  real(wp), dimension(1:nq,1:maxpts) :: unew
  real(wp), dimension(1:2) :: ett
  real(wp), dimension(1:2,1:7) :: error
  real(wp), dimension(1:2,1:7,1:ncpu) :: all_error
  !
  !.. Local Parameters ..
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
  error(:,:) = zero  ! intialize everything to zero
  error(:,3) = -eps6 ! make this slightly negative in case all error is zero
  !
  nc1 = cell_map_ptr(mypnum+1)
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1
    n2 = pts_in_cell_ptr(n+1)
    np = n2-n1+1
    !
    this_geom = cell_geom(n)
    this_order = cell_order(n)
    !
    ! Interpolate the new and old solution to the error points
    !
    if (associated(outerp(this_geom,this_order)%error)) then
      do l = 1,nr
        xyz_err(l,1:np) = matmul( outerp(this_geom,this_order)%error , xyz(l,n1:n2) )
      end do
      do m = 1,nq
        unew(m,1:np) = matmul( outerp(this_geom,this_order)%error , usp(m,n1:n2) )
        uold(m,1:np) = matmul( outerp(this_geom,this_order)%error , uoldsp(m,n1:n2) )
      end do
      jac(1:np) = matmul( outerp(this_geom,this_order)%error , jacobian(n1:n2) )
      weights(1:np) = outerp(this_geom,this_order)%err_wts
    else
      xyz_err(1:nr,1:np) = xyz(1:nr,n1:n2)
      unew(1:nq,1:np) = usp(1:nq,n1:n2)
      uold(1:nq,1:np) = uoldsp(1:nq,n1:n2)
      jac(1:np) = jacobian(n1:n2)
      weights(1:np) = std_elem(this_geom,this_order)%wts
    end if
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
      if (any([rho,pres] <= zero)) write (iout,111) n,k,rho,pres
      111 format ("cell = ",i0,"; SP = ",i0,"; ",*(1x,es14.6))
      !
      sn = entropy_cv_sp( unew(:,k) )
      se = entropy_cv_sp( cv_in(:) )
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
        error(:,4) = real(cell_map(nc1+n),kind=wp)
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
    call mpi_allgather(error,    size(error),mpi_flttyp, &
                       all_error,size(error),mpi_flttyp,mpi_comm_world,mpierr)
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
  !
  rl2mx = max( rltwo , rl2mx )
  !
  if (any(itestcase == Transport_Problems)) then
    entL2 = sqrt( sum(all_error(2,2,:)) / total_volume )
  else
    entL2 = log10(rl2mx) - log10(rltwo)
    convergence_reached = (entL2 <= convergence_order_max_res) .or. &
                          (log10(rltwo) <= convergence_order_abs_res)
  end if
  !
  ! Output the results for the current time step
  !
  if (mypnum == 0) then
    !
#ifdef PBS_ENV
    !
    wall_time_last = wall_time_left
    wall_time = mpi_wtime() - beg_wall_time
    wall_time_left = min( max_wall_time - wall_time , max_sec )
    !
    next_time_est = (wall_time_last - wall_time_left) * (one+eps1)
    walltime_expiring = (wall_time_left-next_time_est <= zero)
    !
    sec_left = abs(nint(wall_time_left))
    !
    hrs_left = sec_left/3600
    min_left = (sec_left - hrs_left*3600)/60
    sec_left = sec_left - hrs_left*3600 - min_left*60
    ! Make hrs_left negative if wall_time_left happens to be negative
    hrs_left = sign(hrs_left,int(wall_time_left))
    !
    write (iout,iter_out_fmt) itcur,d_t,time,rltwo,entL2,nloc, &
                              global_pinc_ptr(nloc)+kloc, &
                              global_xyz(1,global_pinc_ptr(nloc)+kloc), &
                              global_xyz(2,global_pinc_ptr(nloc)+kloc), &
                              hrs_left,min_left,sec_left
    !
#else
    !
    walltime_expiring = fals
    !
    write (iout,iter_out_fmt) itcur,d_t,time,rltwo,entL2,nloc, &
                              global_pinc_ptr(nloc)+kloc, &
                              global_xyz(1,global_pinc_ptr(nloc)+kloc), &
                              global_xyz(2,global_pinc_ptr(nloc)+kloc)
    !
    !
#endif
    !
    if (convergence_reached) then
      if (entL2 <= convergence_order_max_res) write (iout,3)
      if (log10(rltwo) <= convergence_order_abs_res) write (iout,4)
    end if
    !
    if (walltime_expiring) write (iout,2)
    !
  end if
  !
  ! Broadcast walltime_expiring to all processors
  !
  call mpi_bcast(walltime_expiring,1,mpi_logical,0,mpi_comm_world,mpierr)
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
  !.. Local Scalars ..
  integer :: ierr
  integer :: final_iter_digits
  integer :: max_cell_digits
  integer :: max_solpt_digits
  character(len=100) :: char_num
  character(len=1000) :: iter_header
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "init_generic_mod"
  !
continue
  !
  ! Allocate the error arrays
  !
  if (.not. allocated(unloc)) then
    allocate ( unloc(1:nq+1) , source=0 , stat=ierr )
    call alloc_error(pname,"unloc",1,ierr)
  end if
  if (.not. allocated(ukloc)) then
    allocate ( ukloc(1:nq+1) , source=0 , stat=ierr )
    call alloc_error(pname,"ukloc",1,ierr)
  end if
  if (.not. allocated(ulinf)) then
    allocate ( ulinf(1:nq+1) , source=zero , stat=ierr )
    call alloc_error(pname,"ulinf",1,ierr)
  end if
  if (.not. allocated(ulone)) then
    allocate ( ulone(1:nq+1) , source=zero , stat=ierr )
    call alloc_error(pname,"ulone",1,ierr)
  end if
  if (.not. allocated(ultwo)) then
    allocate ( ultwo(1:nq+1) , source=zero , stat=ierr )
    call alloc_error(pname,"ultwo",1,ierr)
  end if
  if (.not. allocated(vnloc)) then
    allocate ( vnloc(1:nq+1) , source=0 , stat=ierr )
    call alloc_error(pname,"vnloc",1,ierr)
  end if
  if (.not. allocated(vkloc)) then
    allocate ( vkloc(1:nq+1) , source=0 , stat=ierr )
    call alloc_error(pname,"vkloc",1,ierr)
  end if
  if (.not. allocated(vlinf)) then
    allocate ( vlinf(1:nq+1) , source=zero , stat=ierr )
    call alloc_error(pname,"vlinf",1,ierr)
  end if
  if (.not. allocated(vlone)) then
    allocate ( vlone(1:nq+1) , source=zero , stat=ierr )
    call alloc_error(pname,"vlone",1,ierr)
  end if
  if (.not. allocated(vltwo)) then
    allocate ( vltwo(1:nq+1) , source=zero , stat=ierr )
    call alloc_error(pname,"vltwo",1,ierr)
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
  else
    write (iter_header,202) final_iter_digits-3, &
                            max_cell_digits-2, &
                            max_solpt_digits-3
  end if
  !
#ifdef PBS_ENV
  iter_header = iter_header(1:len_trim(iter_header)) // ",4x,'Time Left'"
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
  ! Format Statements
  !
  1 format (i0)
  100 format ("(1x,i",i0,",4es12.4,2x,i",i0,",2x,i",i0,",2es9.1")
  200 format ("(/,",i0,"x,'iter',5x,'d_t',8x,'Time',6x,'L2-Density',2x,", &
              "'L2-Entropy',",i0,"x,'Cell',",i0,"x,'Solpt',5x,'x',8x,'y'")
  201 format ("(/,",i0,"x,'iter',5x,'d_t',8x,'Time',6x,'L2-X-Mmntm',2x,", &
              "'OrdersDrop',",i0,"x,'Cell',",i0,"x,'Solpt',5x,'x',8x,'y'")
  202 format ("(/,",i0,"x,'iter',5x,'d_t',8x,'Time',6x,'L2-Density',2x,", &
              "'OrdersDrop',",i0,"x,'Cell',",i0,"x,'Solpt',5x,'x',8x,'y'")
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
  use order_mod, only : maxpts
  use geovar,  only : nr,ncell,pts_in_cell_ptr,xyz
  use geovar,  only : cell_geom,cell_order
  use ovar,    only : nentmax,kentmax,entmax,gam,cv_in
  use ovar,    only : xentmax,yentmax,entL2,total_volume,itestcase
  use metrics_mod,  only : jacobian
  use flowvar, only : usp
  use parallel_mod, only : cell_map,cell_map_ptr
  use quadrature_mod, only : std_elem
  use interpolation_mod, only : outerp
  !
  !.. Local Scalars ..
  integer  :: n,n1,n2,np,m,l,k,nc1
  integer  :: this_geom,this_order
  real(wp) :: ett,sn,se
  !
  !.. Local Arrays ..
  real(wp), dimension(1:maxpts) :: jac
  real(wp), dimension(1:maxpts) :: weights
  real(wp), dimension(1:nr,1:maxpts) :: xyz_err
  real(wp), dimension(1:nq,1:maxpts) :: unew
  real(wp), dimension(1:6) :: error
  real(wp), dimension(1:6,1:ncpu) :: all_error
  !
continue
  !
  error(:) = zero  ! intialize everything to zero
  error(2) = -eps6 ! make this slightly negative in case all error is zero
  !
  nc1 = cell_map_ptr(mypnum+1)
  !
  do n = 1,ncell
    !
    n1 = pts_in_cell_ptr(n)+1
    n2 = pts_in_cell_ptr(n+1)
    np = n2-n1+1
    !
    this_geom = cell_geom(n)
    this_order = cell_order(n)
    !
    ! Interpolate the new and old solution to the error points
    !
    if (associated(outerp(this_geom,this_order)%error)) then
      do l = 1,nr
        xyz_err(l,1:np) = matmul( outerp(this_geom,this_order)%error , xyz(l,n1:n2) )
      end do
      do m = 1,nq
        unew(m,1:np) = matmul( outerp(this_geom,this_order)%error , usp(m,n1:n2) )
      end do
      jac(1:np) = matmul( outerp(this_geom,this_order)%error , jacobian(n1:n2) )
      weights(1:np) = outerp(this_geom,this_order)%err_wts
    else
      xyz_err(1:nr,1:np) = xyz(1:nr,n1:n2)
      unew(1:nq,1:np) = usp(1:nq,n1:n2)
      jac(1:np) = jacobian(n1:n2)
      weights(1:np) = std_elem(this_geom,this_order)%wts
    end if
    !     
    ! Now loop through the solution points of this cell to compute the
    ! error between the solutions of the current and previous time steps
    !
    do k = 1,np
      !
      sn = entropy_cv_sp( unew(:,k) )
      se = entropy_cv_sp( cv_in(:) )
      !
      ett = sn - se
      !
      ! L2 norm of entropy
      !
      error(1) = error(1) + jac(k)*weights(k)*ett*ett
      !
      if (abs(ett) > error(2)) then
        error(2) = abs(ett)
        error(3) = real(cell_map(nc1+n),kind=wp)
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
    call mpi_allgather(error,    size(error),mpi_flttyp, &
                       all_error,size(error),mpi_flttyp,mpi_comm_world,mpierr)
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
end subroutine getmaxent
!
!###############################################################################
!
include "Functions/usp2v.f90"
include "Functions/entropy.f90"
include "Functions/pressure.f90"
!
!###############################################################################
!
end module generic_mod
