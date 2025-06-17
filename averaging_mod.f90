module averaging_mod
  !
  use module_kind_types
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq,ntk,ntl
  use eqn_idx_mod, only : npres,ntemp
  !
  implicit none
  !
  private
  !
  public :: initialize_time_ave
  public :: update_time_ave
  public :: average_restart_files
  !
  character(len=100), save, allocatable :: restart_file_list(:)
  real(wp), save, allocatable :: restart_times(:)
  !
  real(wp), save, allocatable :: ave(:)
  !
  type :: time_ave_var_t
    integer, allocatable :: idx(:)
  end type time_ave_var_t
  !
  type(time_ave_var_t), save, allocatable :: time_ave_var(:)
  !
contains
!
!###############################################################################
!
subroutine initialize_time_ave()
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr,n_solpts
  use ovar,        only : time,ave_start_time
  use ovar,        only : time_ave_vel_is_axisymm
  use flowvar_mod, only : time_ave_variables
  use flowvar_mod, only : uavesp
  use mrf_mod,     only : is_in_mrf
  !
  !.. Local Scalars ..
  integer :: l,m,nv,nvar,ierr
  character(len=1)  :: ta_char,v1,v2,v3
  character(len=30) :: array_name
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_time_ave"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate and create the time_ave_variables array if it hasnt already
  ! been allocated; otherwise, count the number of time-averaged variables it
  ! contains. This character array identifies which time-averaged variables
  ! are computed and written to the time-averaged solution/restart files.
  !
  if (.not. allocated(time_ave_variables)) then
    !
    v1 = "x"
    v2 = merge( "r" , "y" , time_ave_vel_is_axisymm )
    v3 = merge( "f" , "z" , time_ave_vel_is_axisymm )
    !
    ! Since the time_ave_variables array isnt allocated, this means that a
    ! time-averaged restart file wasnt read and the start time of the
    ! time-averaging needs to be initialized to the current time
    !
    ave_start_time = time
    !
   !nvar = merge(23,11,nr==3)
   !allocate ( time_ave_variables(1:nvar) , stat=ierr , errmsg=error_message )
   !call alloc_error(pname,"time_ave_variables",1,__LINE__,__FILE__, &
   !                 ierr,error_message)
    !
    if (nr == 3) then
      !
      ! First 5 are the conserved variables
      !
      ! USING F2003 AUTO-REALLOCATION
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             "d" , "dx" , "dy" , "dz" , "e" ]
      !
      ! Next 5 are the primitive variables
      ! (with temperature instead of density)
      !
      ! USING F2003 AUTO-REALLOCATION
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , &
                             "t" , "x" , "y" , "z" , "p" ]
      !
      ! Remaining are derived variables
      !
      ! ALL OF THE BELOW ASSIGNMENTS USE F2003 AUTO-REALLOCATION TO
      ! APPEND A NEW ARRAY ELEMENT ONTO THE time_ave_variables ARRAY
      ! FOR EACH ADDITIONAL TIME AVERAGED VARIABLE
      !
      if (time_ave_vel_is_axisymm) then
        time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                               time_ave_variables(:) , "d" // v2 ]
        time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                               time_ave_variables(:) , "d" // v3 ]
        time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                               time_ave_variables(:) , v2 ]
        time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                               time_ave_variables(:) , v3 ]
      end if
      !
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v1 // v1 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v2 // v2 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v3 // v3 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v1 // v2 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v1 // v3 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v2 // v3 ]
      !
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "dd" ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "tt" ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "pp" ]
      !
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v1 // "t" ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v2 // "t" ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , v3 // "t" ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "dt" ]
      !
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v1 // v1 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v2 // v2 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v3 // v3 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v1 // v2 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v1 // v3 ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v2 // v3 ]
      !
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v1 // "t" ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v2 // "t" ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "d" // v3 // "t" ]
      !
    else
      !
      ! First 4 are the conserved variables
      !
      ! USING F2003 AUTO-REALLOCATION
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             "d" , "dx" , "dy" , "e" ]
      !
      ! Next 4 are the primitive variables
      ! (with temperature instead of density)
      !
      ! USING F2003 AUTO-REALLOCATION
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "t" , "x" , "y" , "p" ]
      !
      ! Remaining are derived variables
      !
      ! USING F2003 AUTO-REALLOCATION
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "xx" ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "yy" ]
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "xy" ]
      !
    end if
    !
    ! Add total temperature and total pressure to the list of variables
    ! to average if using a moving reference frame
    !
    if ( is_in_mrf() ) then
      !
      ! Total temperature
      !
      ! USING F2003 AUTO-REALLOCATION
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "T" ]
      !
      ! Total pressure
      !
      ! USING F2003 AUTO-REALLOCATION
      time_ave_variables = [ character(len=len(time_ave_variables)) :: &
                             time_ave_variables(:) , "P" ]
      !
    end if
    !
    nvar = size(time_ave_variables)
    !
  else
    !
    nvar = size(time_ave_variables)
    !
  end if
  !
  if (allocated(time_ave_var)) then
    deallocate ( time_ave_var , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"time_ave_var",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  allocate ( time_ave_var(1:nvar) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"time_ave_var",1,__LINE__,__FILE__,ierr,error_message)
  !
  do m = 1,size(time_ave_variables)
    !
    time_ave_variables(m) = trim(adjustl(time_ave_variables(m)))
    !
    nv = len_trim(time_ave_variables(m))
    !
    allocate ( time_ave_var(m)%idx(1:nv) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) m
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    do l = 1,nv
      !
      ta_char = time_ave_variables(m)(l:l)
      !
      select case (ta_char)
        case ("d") ! density
          time_ave_var(m)%idx(l) = nec
        case ("x") ! axial velocity
          time_ave_var(m)%idx(l) = nmx
        case ("y") ! y velocity
          time_ave_var(m)%idx(l) = nmy
        case ("z") ! spanwise velocity
          time_ave_var(m)%idx(l) = nmz
        case ("e") ! total energy
          time_ave_var(m)%idx(l) = nee
        case ("p") ! static pressure
          time_ave_var(m)%idx(l) = nee
        case ("t") ! static temperature
          time_ave_var(m)%idx(l) = nq+1
        case ("r") ! radial velocity
          time_ave_var(m)%idx(l) = nq+2
        case ("f") ! theta velocity
          time_ave_var(m)%idx(l) = nq+3
        case ("P") ! total pressure
          time_ave_var(m)%idx(l) = nq+4
        case ("T") ! total temperature
          time_ave_var(m)%idx(l) = nq+5
        case ("k") ! turbulent kinetic energy
          time_ave_var(m)%idx(l) = ntk
        case ("w") ! turbulence frequency/dissipation
          time_ave_var(m)%idx(l) = ntl
        case default
          ! ??? Should probably produce an error message and abort
          time_ave_var(m)%idx(l) = nec
      end select
    end do
  end do
  !
  ! Allocate the local module array used for computing
  ! the current update to the time-averaged variables.
  !
  if (allocated(ave)) then
    deallocate ( ave , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ave",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( ave(1:nvar) , source=zero , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ave",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate and initialize the time-averaged variables
  ! solution array if it hasnt already been allocated.
  !
  if (.not. allocated(uavesp)) then
    !
    allocate ( uavesp(1:n_solpts,1:nvar) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"uavesp",1,__LINE__,__FILE__,ierr,error_message)
    !$ACC UPDATE DEVICE(uavesp)
    !
    ! Initialize the time-averaged variables to the initial/restart
    ! solution since they werent read from a time-averaged restart file.
    !
    call update_time_ave
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("time_ave_var(",i0,")%idx")
  !
end subroutine initialize_time_ave
!
!###############################################################################
!
subroutine update_time_ave()
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr
  use ovar,        only : time,prev_ave_time
  use ovar,        only : ave_start_time,saved_ave_time
  !
  !.. Local Arrays ..
  real(wp) :: time_val(1:2)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_time_ave"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (time == zero .or. time == ave_start_time) then
    time_val(1) = zero
    time_val(2) = one
  else
    time_val(1) = (prev_ave_time-ave_start_time) / (time-ave_start_time)
    time_val(2) = one - time_val(1) ! to try and ensure that sum(time_val) = 1
  end if
  !
  saved_ave_time = prev_ave_time
  !
  ! Broadcast the time_val array so hopefully all MPI processes
  ! have the exact same values for the two averaging coefficients
  ! (in case there might be roundoff differences between processes).
  !
  call mpi_bcast(time_val,size(time_val,kind=int_mpi),mpi_flttyp, &
                 glb_root,MPI_COMM_WORLD,mpierr)
  !
  if (nr == 3) then
    call update_time_ave_3d(time_val(1),time_val(2))
  else if (nr == 2) then
    call update_time_ave_2d(time_val(1),time_val(2))
  end if
  !
  call debug_timer(entering_procedure,pname)
  !
end subroutine update_time_ave
!
!###############################################################################
!
subroutine update_time_ave_3d(a,b)
  !
  !.. Use Statements ..
  use geomvar_mod, only : xyz,n_solpts
  use flowvar_mod, only : usp,gpsp,uavesp
  use ovar,        only : gam,time,prev_ave_time
  use ovar,        only : ave_start_time,saved_ave_time
  use ovar,        only : moving_ref_frame
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: a
  real(wp), intent(in) :: b
  !
  !.. Local Scalars ..
  integer  :: m,n,nvar
  real(wp) :: x,x_r,vx_r,vx
  real(wp) :: y,y_r,vy_r,vy
  real(wp) :: z,z_r,vz_r,vz
  real(wp) :: st,ct,theta,gm1,gam_over_gm1
  real(wp) :: r,p,t,q2,a2,m2,total_cond
  !
  !.. Local Arrays ..
  real(wp) :: ave_var(1:nq+5)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_time_ave_3d"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  gm1 = gam - one
  gam_over_gm1 = gam/gm1
  !
  nvar = size(time_ave_var)
  !
  call moving_ref_frame%get_rotation(x_r,y_r,z_r)
  !
  do n = 1,n_solpts
    !
    x = xyz(1,n)
    y = xyz(2,n)
    z = xyz(3,n)
    !
    if (all([y,z] == zero)) then
      !
      st = zero
      ct = zero
      !
    else
      !
      theta = atan2(z,y)
      if (theta < zero) theta = theta + two*pi
      !
      st = sin(theta)
      ct = cos(theta)
      !
    end if
    !
    r = usp(n,nec)
    p = gpsp(n,npres)
    t = gpsp(n,ntemp)
    !
    ! Compute the absolute velocity, aspd, Mach, and total conditions from usp
    !
    vx = usp(n,nmx)/r
    vy = usp(n,nmy)/r
    vz = usp(n,nmz)/r
    !
    q2 = vx*vx + vy*vy + vz*vz
    a2 = max( gam*p/r , rndoff )
    !
    m2 = q2/a2
    !
    total_cond = one + half*gm1*m2
    !
    ! Compute the reference frame velocity
    !
    vx_r = y_r*z - z_r*y
    vy_r = z_r*x - x_r*z
    vz_r = x_r*y - y_r*x
    !
    ! Compute the relative velocity
    !
    vx = vx - vx_r
    vy = vy - vy_r
    vz = vz - vz_r
    !
    ! Primitive variables
    ave_var(nec) = r
    ave_var(nmx) = vx
    ave_var(nmy) = vy
    ave_var(nmz) = vz
    ave_var(nee) = p
    !
    ! Temperature
    ave_var(nq+1) = t
    !
    ! Radial    instantaneous velocity
    ave_var(nq+2) = st*vz + ct*vy
    ! Azimuthal instantaneous velocity
    ave_var(nq+3) = ct*vz - st*vy
    !
    ! Total pressure
    ave_var(nq+4) = p*(total_cond**gam_over_gm1)
    ! Total temperature
    ave_var(nq+5) = t*total_cond
    !
    ! First, update the time-averaged conservative variables
    !
    do m = 1,nq
      uavesp(n,m) = a*uavesp(n,m) + b*usp(n,m)
    end do
    !
    ! Next, update the time-averaged primitive variables
    ! (with temperature instead of density)
    !
    uavesp(n,nq+1) = a*uavesp(n,nq+1) + b*ave_var(nq+1)
    do m = 2,nq
      uavesp(n,nq+m) = a*uavesp(n,nq+m) + b*ave_var(m)
    end do
    !
    ! Update the remaining time-averaged variables
    !
    do m = 2*nq+1,nvar
      uavesp(n,m) = a*uavesp(n,m) + b*product( ave_var(time_ave_var(m)%idx) )
    end do
    !
  end do
  !
  ! Update prev_ave_time to the current time
  !
  prev_ave_time = time
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_time_ave_3d
!
!###############################################################################
!
subroutine update_time_ave_2d(a,b)
  !
  !.. Use Statements ..
  use geomvar_mod, only : xyz,n_solpts
  use flowvar_mod, only : usp,gpsp,uavesp
  use ovar,        only : gam,time,prev_ave_time
  use ovar,        only : ave_start_time,saved_ave_time
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: a
  real(wp), intent(in) :: b
  !
  !.. Local Scalars ..
  integer  :: m,n,nvar
  real(wp) :: x,x_r,vx
  real(wp) :: y,y_r,vy
  real(wp) :: st,ct,theta,gm1,gam_over_gm1
  real(wp) :: r,p,t,q2,a2,m2,total_cond
  !
  !.. Local Arrays ..
  real(wp) :: ave_var(1:nq+5)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_time_ave_2d"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  gm1 = gam - one
  gam_over_gm1 = gam/gm1
  !
  nvar = size(time_ave_var)
  !
  do n = 1,n_solpts
    !
    x = xyz(1,n)
    y = xyz(2,n)
    !
    if (all([x,y] == zero)) then
      !
      st = zero
      ct = zero
      !
    else
      !
      theta = atan2(y,x)
      if (theta < zero) theta = theta + two*pi
      !
      st = sin(theta)
      ct = cos(theta)
      !
    end if
    !
    r = usp(n,nec)
    p = gpsp(n,npres)
    t = gpsp(n,ntemp)
    !
    vx = usp(n,nmx)/r
    vy = usp(n,nmy)/r
    !
    q2 = vx*vx + vy*vy
    a2 = max( gam*p/r , rndoff )
    !
    m2 = q2/a2
    !
    total_cond = one + half*gm1*m2
    !
    ! Primitive variables
    ave_var(nec) = r
    ave_var(nmx) = vx
    ave_var(nmy) = vy
    ave_var(nee) = p
    !
    ! Temperature
    ave_var(nq+1) = t
    !
    ! Radial    instantaneous velocity
    ave_var(nq+2) = st*vy + ct*vx
    ! Azimuthal instantaneous velocity
    ave_var(nq+3) = ct*vy - st*vx
    !
    ! Total temperature
    ave_var(nq+4) = t*total_cond
    ! Total pressure
    ave_var(nq+5) = p*(total_cond**gam_over_gm1)
    !
    ! First, update the time-averaged conservative variables
    !
    do m = 1,nq
      uavesp(n,m) = a*uavesp(n,m) + b*usp(n,m)
    end do
    !
    ! Next, update the time-averaged primitive variables
    ! (with temperature instead of density)
    !
    uavesp(n,nq+1) = a*uavesp(n,nq+1) + b*ave_var(nq+1)
    do m = 2,nq
      uavesp(n,nq+m) = a*uavesp(n,nq+m) + b*ave_var(m)
    end do
    !
    ! Update the remaining time-averaged variables
    !
    do m = 2*nq+1,nvar
      uavesp(n,m) = a*uavesp(n,m) + b*product( ave_var(time_ave_var(m)%idx) )
    end do
    !
  end do
  !
  ! Update prev_ave_time to the current time
  !
  prev_ave_time = time
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_time_ave_2d
!
!###############################################################################
!
subroutine average_restart_files()
  !
  !.. Use Statements ..
  use ovar,           only : restart_list_file,restart_file
  use ovar,           only : output_time_averaging,read_time_ave_restart
  use ovar,           only : ave_start_time,prev_ave_time
  use ovar,           only : time,itcur,tmrst,itrst
  use flowvar_mod,    only : uavesp
  use restart_mod,    only : read_restart_file,write_restart_file
  use write_cgns_mod, only : write_parallel_cgns
  !
  !.. Local Scalars ..
  integer :: ierr,io,m,n,np,nq
  !
  character(len=150) :: fname
  character(len=150) :: text
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "average_restart_files"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! First deallocate unnecessary arrays that were allocated for the solver but
  ! wont be needed since execution will abort at the end of this subroutine
  !
  call deallocate_solver_arrays
  !
  ! Read in the list of restart files that will be averaged
  !
  fname = trim(adjustl(restart_list_file))
  !
  open (newunit=io,file=fname,status="old",action="read", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  n = 0
  do
    read (io,1,iostat=ierr,iomsg=error_message) text
    if (ierr < 0) then
      exit
    else if (ierr > 0) then
      call io_error(pname,fname,-1,__LINE__,__FILE__,ierr,error_message)
    else
      n = n + 1
    end if
  end do
  !
  rewind (io)
  !
  allocate ( restart_times(1:n) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"restart_times",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( restart_file_list(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"restart_file_list",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  do n = 1,size(restart_file_list)
    read (io,1,iostat=ierr,iomsg=error_message) restart_file_list(n)
    call io_error(pname,fname,-1,__LINE__,__FILE__,ierr,error_message)
  end do
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  read_time_ave_restart = fals
  !
  restart_file = trim(adjustl(restart_file_list(1)))
  !
  if (mypnum == glb_root) then
    write (iout,2) trim(adjustl(restart_file))
  end if
  !
  call read_restart_file
  !
  time = tmrst
  itcur = itrst
  !
  restart_times(1) = time
  ave_start_time = time
  prev_ave_time = time
  !
  if (mypnum == glb_root) write (iout,3)
  !
  call initialize_time_ave
  !
  ! Initialize uavesp to zero
  !
  np = size(uavesp,dim=1)
  nq = size(uavesp,dim=2)
  !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(uavesp)
  do m = 1,nq
    do n = 1,np
      uavesp(n,m) = zero
    end do
  end do
  !
  do n = 2,size(restart_file_list)
    !
    restart_file = trim(adjustl(restart_file_list(n)))
    !
    if (mypnum == glb_root) then
      write (iout,4) trim(adjustl(restart_file))
    end if
    !
    call read_restart_file
    !
    time = tmrst
    itcur = itrst
    !
    restart_times(n) = time
    !
    if (mypnum == glb_root) write (iout,5)
    call update_time_ave
    !
  end do
  !
  output_time_averaging = true
  !
  if (mypnum == glb_root) write (iout,7) pname
  call write_restart_file(itcur)
  !
  if (mypnum == glb_root) write (iout,8) pname
  call write_parallel_cgns(open_CGNS)
  !
  call debug_timer(leaving_procedure,pname)
  !
  write (error_message,6)
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  !
  ! Format Statements
  !
  1 format (a)
  2 format (/,4x,"Reading initial restart file '",a,"'")
  3 format (/,4x,"Initializing the time-averaged variables")
  4 format (/,4x,"Working on restart file '",a,"'",/, &
              8x,"Reading the restart file")
  5 format (8x,"Updating the time-averaged variables")
  6 format ("STOPPING EXECUTION AFTER TIME-AVERAGING ", &
            "THE LIST OF RESTART FILES!")
  7 format (/,4x,"Calling write_restart_file from ",a)
  8 format (/,4x,"Calling write_parallel_cgns from ",a)
  !
end subroutine average_restart_files
!
!###############################################################################
!
subroutine deallocate_solver_arrays()
  !
  !.. Use Statements ..
  use geomvar_mod, only : face,grid
  use flowvar_mod, only : faceusp_dt,facedusp_dt,faceflx,facexyz_dt
  use flowvar_mod, only : residual,dtsp,uoldsp,rk_sum,rssprk
  use flowvar_mod, only : nodeusp_dt,nodedusp_dt
  use flowvar_mod, only : edgeusp_dt,edgedusp_dt,edgexyz_dt
  use flowvar_mod, only : global_usp,global_var
  use metrics_mod, only : metrics_dt,geofa,facenrm
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "deallocate_solver_arrays"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (allocated(face)) then
    deallocate ( face , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"face",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(grid)) then
    deallocate ( grid , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(faceusp_dt)) then
    deallocate ( faceusp_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"faceusp_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(facedusp_dt)) then
    deallocate ( facedusp_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facedusp_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(faceflx)) then
    deallocate ( faceflx , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"faceflx",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(facexyz_dt)) then
    deallocate ( facexyz_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facexyz_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(residual)) then
    deallocate ( residual , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"residual",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(dtsp)) then
    deallocate ( dtsp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"dtsp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(uoldsp)) then
    deallocate ( uoldsp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"uoldsp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(rk_sum)) then
    deallocate ( rk_sum , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"rk_sum",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(rssprk)) then
    deallocate ( rssprk , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"rssprk",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(nodeusp_dt)) then
    deallocate ( nodeusp_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodeusp_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(nodedusp_dt)) then
    deallocate ( nodedusp_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodedusp_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edgeusp_dt)) then
    deallocate ( edgeusp_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edgeusp_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edgedusp_dt)) then
    deallocate ( edgedusp_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edgedusp_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edgexyz_dt)) then
    deallocate ( edgexyz_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edgexyz_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(global_usp)) then
    deallocate ( global_usp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_usp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(global_var)) then
    deallocate ( global_var , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_var",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(metrics_dt)) then
    deallocate ( metrics_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"metrics_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(geofa)) then
    deallocate ( geofa , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"geofa",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(facenrm)) then
    deallocate ( facenrm , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facenrm",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine deallocate_solver_arrays
!
!###############################################################################
!
end module averaging_mod
