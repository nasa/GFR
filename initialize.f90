module initialization_mod
  !
  ! Module containing the functions and routines to
  ! initialize the solution
  !
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,ntk,ntl,nq
  use eqn_idx, only : ncs,npres,ntemp,ngamm,nenth,ncp,nvisc,ncond,ng
  !
  implicit none
  !
  private
  !
  public :: initialize_dusp
  public :: initialize_flow
  public :: initialize_reference
  public :: create_bc_in
  public :: get_face_and_edge_point_coordinates
  !
contains
!
!###############################################################################
!
subroutine initialize_flow()
  !
  ! Routine to initialize the flow
  !
  !.. Use Statements ..
  use matrix_types_mod, only : matrix
  !
  use geovar, only : nr,n_solpts,ncell,cell,xyz,grid
  !
  use ovar, only : bc_in
  use ovar, only : rhoref,aref,itestcase,jet_init_func
  use ovar, only : load_restart_file,d_t,constant_dt,itcur,time
  use ovar, only : init_tstar,time_ref,time_scaling_factor
  use ovar, only : dtrst,itrst,tmrst
  use ovar, only : adjust_cfl_cycles
  use ovar, only : solution_file, solution_format
  use ovar, only : use_SEM_profile_for_initialization
  !
  use flowvar, only : usp,uoldsp,gpsp
  use flowvar, only : allocate_flowvar_arrays
  !
  use order_mod, only : maxpts
  !
  use restart_mod, only : read_restart_file
  use restart_mod, only : create_restart_fname_formats
  !
  use projection_mod, only : project_to_porder
  !
  use interpolation_mod, only : outerp
  !
  use synthetic_eddy_mod, only : initialize_from_sem_profile
  !
  use thermo_mod,    only : update_thermo_properties
  use transport_mod, only : update_transport_properties
  !
  !.. Local Scalars ..
  integer :: n,n0,n1,n2,np,l,l1,l2,m
  integer :: gmin,gmax,omin,omax,ierr
  integer :: this_geom,this_order
  integer :: init_cond,this_cell
  integer :: mesh_order
  real(wp) :: rhoa,rhoa2,rhoref_inv,rhoa_inv,rhoa2_inv
  logical(lk) :: grid_sol_is_available
  character(len=200) :: array_name
  !
  !.. Local Allocatable Arrays ..
  real(wp),     allocatable :: xyz_init(:,:)
  type(matrix), allocatable :: outerp_inv(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_flow"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  !
  ! Initialization of the bc_in array has been
  ! moved to subroutine read_grid_file.
  !
  ! Some constants from the free-stream conditions
  !
  rhoa = rhoref * aref
  rhoa2 = rhoa * aref
  !
  rhoref_inv = one / rhoref
  rhoa_inv = one / rhoa
  rhoa2_inv = one / rhoa2
  !
  ! Allocate the arrays within the module flowvar
  !
  call allocate_flowvar_arrays
  !
  ! Initialize the iteration counter, time step, and elapsed time
  ! Overwrite with restart file values later, if available.
  !
  d_t   = constant_dt
  itcur = 0
  time  = zero
  !
  ! Check whether grid%sol has been read in.
  !
  grid_sol_is_available = fals
  if (allocated(grid)) then
    if (allocated(grid%sol)) grid_sol_is_available = true
  end if
  !
  ! Check for using a restart file to initialize the solution
  !
  if (len_trim(adjustl(solution_file)) > 0 .and. &
      solution_format == Restart_Format) then
    !
    ! Read in usp from the restart file
    !
    call read_restart_file(solution_file)
   !call read_restart_file(solution_file,force_serial_read=true)
    !
    ! Do not initialize the iteration counter, time step,
    ! and elapsed time to that in the restart file.
    ! Use the default values from above.
    !
  else if (load_restart_file) then
    !
    ! Read in usp and the current iteration, time step,
    ! and elapsed time from the restart file
    !
    call read_restart_file
    !
    ! Initialize the iteration counter, time step,
    ! and elapsed time to that in the restart file.
    !
    d_t   = dtrst
    itcur = itrst
    time  = tmrst
    !
    if (abs(itestcase) == CH4_3D_Diffusion_Flame .and. &
        jet_init_func == -1*(CH4_3D_Diffusion_Flame*100+1)) then
      if (mypnum == glb_root) then
        write (iout,*) "Calling detach_lip_flame"
        flush (iout)
      end if
      !
      ! Update the thermodynamics and transport properties
      ! at the solution points
      !
      call update_thermo_properties(1,n_solpts,usp,gpsp, &
                                    loc="cell",caller_name=pname)
      call update_transport_properties(1,n_solpts,usp,gpsp)
      call detach_lip_flame(xyz,usp,gpsp)
    end if
    !
  else if (grid_sol_is_available) then
    !
    ! THIS SECTION OF CODE WAS MOVED INTO MAPPINGS_MOD SO DO NOTHING!
    !
  else if (use_SEM_profile_for_initialization .and. &
           any(bc_in(:)%bc_type == bc_sem)) then
    !
    ! Initialize the flow based on the 1D inflow profile
    ! used for the SEM turbulence inflow generator.
    !
    call initialize_from_sem_profile(xyz,usp)
    !
  else
    !
    ! Initialize the solution based on the input parameter itestcase
    !
    if (itestcase == Subsonic_Jet) then
      call find_elem_for_initial_jet_profile(point_is_in_smc000_nozzle)
    else if (itestcase == -Subsonic_Jet) then
      call find_elem_for_initial_jet_profile(point_is_in_arn2_nozzle)
    else if (abs(itestcase) == Supersonic_Jet) then
      call find_elem_for_initial_jet_profile(point_is_in_ssjet_nozzle)
    end if
    !
    allocate ( xyz_init(1:nr,1:maxpts) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_init",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Allocate the local array that stores the inverse of the matrix used for
    ! interpolating from the solution points to the initialization points.
    ! NOTE: This array is used so that we only invert the matrix once for a
    !       given cell geometry/order combination.  We dont want to keep
    !       inverting the same matrix over and over for cells with the same
    !       geometry/order combination.
    !
    gmin = lbound(outerp,dim=1)
    gmax = ubound(outerp,dim=1)
    omin = lbound(outerp,dim=2)
    omax = ubound(outerp,dim=2)
    !
    allocate ( outerp_inv(gmin:gmax,omin:omax) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"outerp_inv",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Initialize the solution in each cell
    !
    do n = 1,ncell
      !
      this_geom = cell(n)%geom
      this_order = cell(n)%order
      !
      n1 = cell(n)%beg_sp
      n2 = cell(n)%end_sp
      np = n2-n1+1
      !
      init_cond = 0
      if (allocated(grid)) then
        if (allocated(grid%elem)) then
          if (allocated(grid%elem(n)%tags)) then
            if (size(grid%elem(n)%tags) >= 2) then
              init_cond = grid%elem(n)%tags(2)
            end if
          end if
        end if
      end if
      !
      ! Check whether the init component of outerp is associated with
      ! an interpolation matrix to determine where we need to compute
      ! the initial conditions.
      !
      if (.not. allocated(outerp(this_geom,this_order)%init)) then
        !
        ! The init component of outerp is not associated so just compute
        ! the initial conditions at the solution points
        !
        usp(n1:n2,1:nq) = initial_conditions( xyz(1:nr,n1:n2) , init_cond )
        !
      else
        !
        ! The init component of outerp is associated so we need to compute
        ! the initial conditions at the initialization points before
        ! interpolating it back to the solution points.
        !
        ! Check if we have already created an inverse of this interpolation
        ! matrix. If we havent, allocate its location in outerp_inv and
        ! compute this inverse matrix.
        !
        if (.not. allocated(outerp_inv(this_geom,this_order)%mat)) then
          !
          l1 = size(outerp(this_geom,this_order)%init%mat,dim=1)
          l2 = size(outerp(this_geom,this_order)%init%mat,dim=2)
          !
          allocate ( outerp_inv(this_geom,this_order)%mat(1:l1,1:l2) , &
                     source=zero , stat=ierr , errmsg = error_message )
          write (array_name,1) Geom_Name(this_geom),this_order
          call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                           error_message)
          !
#ifdef SPECIAL_FOR_INTEL
          outerp_inv(this_geom,this_order)%mat = &
            invert_matrix_lapack95( outerp(this_geom,this_order)%init%mat )
#else
          outerp_inv(this_geom,this_order)%mat = &
            invert_matrix( outerp(this_geom,this_order)%init%mat )
#endif
          !
        end if
        !
        ! Interpolate the cell solution points to the initialization points
        !
        do l = 1,nr
          xyz_init(l,1:np) = matmul( outerp(this_geom,this_order)%init%mat , &
                                     xyz(l,n1:n2) )
        end do
        !
        ! Evaluate the initial conditions at the initialization points
        !
        usp(n1:n2,1:nq) = initial_conditions( xyz_init(1:nr,1:np) , init_cond )
        !
        ! Interpolate the initial conditions from the initialization points
        ! to the solution points
        !
        do m = 1,nq
          usp(n1:n2,m) = matmul( outerp_inv(this_geom,this_order)%mat , &
                                 usp(n1:n2,m) )
        end do
        !
      end if
      !
    end do
    !
    ! Project the initial solution from the
    ! initialization points to the solution points
    !
   !call project_to_porder(usp)
    !
  end if
  !
  ! Update the thermodynamics and transport properties at the solution points
  !
  call update_thermo_properties(1,n_solpts,usp,gpsp, &
                                loc="cell",caller_name=pname)
  call update_transport_properties(1,n_solpts,usp,gpsp)
  !
  init_tstar = time*time_ref*time_scaling_factor
  !
  ! Create the character format strings for the
  ! restart file names and the linking commands
  !
  call create_restart_fname_formats(itcur)
  !
  ! Initialization of time-averaging has been moved to the main program
  !
  ! Initialize the CFL cycles
  !
  call adjust_cfl_cycles
  !
  ! Deallocate the grid derived type
  !
  if (allocated(grid)) then
    deallocate ( grid , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Initialize the solution at previous timestep array
  !
  !$ACC UPDATE DEVICE(usp)
  if (allocated(uoldsp)) then
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(uoldsp,usp)
    do m = 1,nq
      do n = 1,n_solpts
        uoldsp(n,m) = usp(n,m)
      end do
    end do
  end if
  !
  ! Initialize the solution gradient array
  !
  call initialize_dusp
  !
  ! Output the boundary conditions
  !
  call write_bc_conditions(ibc_beg=1,ibc_end=ubound(bc_in,dim=1), &
                           dimensional=true)
  call write_bc_conditions(ibc_beg=1,ibc_end=ubound(bc_in,dim=1), &
                           dimensional=fals)
  !
  ! Output the free-stream conditions
  !
  if (any(itestcase == Advec_Problems)) then
    call write_linear_advec_ref_conditions
  else
    call write_bc_conditions(ibc_end=0,dimensional=true)
    call write_bc_conditions(ibc_end=0,dimensional=fals)
  end if
  !
  ! Output the free-stream conditions
  !
  call old_output_freestream_conditions
  !
!!block
!!  !
!!  use gpu_mod, only : cv2pv_gpu
!!  integer :: io,n,m
!!  real(wp), allocatable :: vsp(:,:)
!!  !
!!  open (newunit=io,file="usp.init.dat",status="replace",action="write")
!!  do n = 1,n_solpts
!!    write (io,81) n,(usp(n,m),m=1,nq)
!!  end do
!!  close (io)
!!  !
!!  allocate ( vsp(1:n_solpts,1:nq) , source=zero )
!!  call cv2pv_gpu(1,n_solpts,usp(1:n_solpts,1:nq),gpsp(1:n_solpts,1:nq), &
!!                            vsp(1:n_solpts,1:nq))
!!  open (newunit=io,file="vsp.init.dat",status="replace",action="write")
!!  do n = 1,n_solpts
!!    write (io,81) n,(vsp(n,m),m=1,nq)
!!  end do
!!  close (io)
!!  deallocate ( vsp )
!!  !
!!  open (newunit=io,file="gpsp.init.dat",status="replace",action="write")
!!  do n = 1,n_solpts
!!    write (io,81) n,gpsp(n,npres),gpsp(n,ntemp),gpsp(n,nvisc),gpsp(n,ncond)
!!  end do
!!  close (io)
!!  !
!!  81 format (i5,2x,10es18.10)
!!  !
!!  !
!!  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,"After dumping usp.init.dat")
!!  !
!!end block
  !
  ! Deallocate local arrays before leaving
  !
  if (allocated(outerp_inv)) then
    deallocate ( outerp_inv , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"outerp_inv",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(xyz_init)) then
    deallocate ( xyz_init , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_init",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Create the facexyz and edgexyz arrays if
  ! using the continuous output option or MMS
  !
  call get_face_and_edge_point_coordinates
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("outerp_inv(",a,",",i0,")%mat")
  !
end subroutine initialize_flow
!
!###############################################################################
!
subroutine old_output_freestream_conditions
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : itestcase
  use ovar,   only : rhoref,aref
  use ovar,   only : pref,tref,gam,rgasref,cpref
  use ovar,   only : machref,velref,uref,muref,reyref
  use ovar,   only : alpha_aoaref,beta_aoaref
  use ovar,   only : pref_nd,velref_nd,muref_nd
  use ovar,   only : rgasref_nd,cpref_nd
  !
continue
  !
  if (mypnum == glb_root) then
    if (all(itestcase /= [Density_Transport,Entropy_Transport])) then
      !
      write (iout,10) rhoref,pref,tref,machref,aref,velref(1)
      if (nr >= 2) write (iout,11) velref(2)
      if (nr == 3) write (iout,12) velref(3)
      write (iout,13) alpha_aoaref,beta_aoaref,gam,rgasref,cpref,muref,reyref
      !
      write (iout,20) rhoref/rhoref,pref_nd,tref/tref,aref/aref,velref_nd(1)
      if (nr >= 2) write (iout,21) velref_nd(2)
      if (nr == 3) write (iout,22) velref_nd(3)
      write (iout,23) rgasref_nd,cpref_nd,muref_nd
      !
    else
      !
      write (iout,30) rhoref,pref,tref,machref,uref/machref, &
                      velref(1),velref(2),alpha_aoaref,beta_aoaref,gam, &
                      rgasref,muref,reyref
      !
    end if
    !
  end if
  !
  10 format (/, &
             " Dimensional Free-stream conditions",/, &
             " Density          = ",es14.6," kg/m^3",/, &
             " Pressure         = ",es14.6," N/m^2",/, &
             " Temperature      = ",es14.6," K",/, &
             " Mach Number      = ",es14.6,/, &
             " Speed of Sound   = ",es14.6," m/s",/, &
             " X-Velocity       = ",es14.6," m/s")
  11 format (" Y-Velocity       = ",es14.6," m/s")
  12 format (" Z-Velocity       = ",es14.6," m/s")
  13 format (" Alpha A.o.Attack = ",es14.6," degrees"/, &
             " Beta  A.o.Attack = ",es14.6," degrees"/, &
             " Gamma            = ",es14.6,/, &
             " Gas Constant     = ",es14.6," J/(kg*K)",/, &
             " Cp               = ",es14.6," J/(kg*K)",/, &
             " Viscosity        = ",es14.6," kg/(m*s)",/, &
             " Reynolds #       = ",es14.6,/)
  !
  20 format (/, &
             " Non-dimensional Free-stream conditions",/, &
             " Density         = ",es14.6,/, &
             " Pressure        = ",es14.6,/, &
             " Temperature     = ",es14.6,/, &
             " Speed of Sound  = ",es14.6,/, &
             " X-Velocity      = ",es14.6)
  21 format (" Y-Velocity      = ",es14.6)
  22 format (" Z-Velocity      = ",es14.6)
  23 format (" Gas Constant    = ",es14.6,/, &
             " Cp              = ",es14.6,/, &
             " Viscosity       = ",es14.6,/)
  !
  30 format (/, &
             " The density vortex transport problem is being",/, &
             " used to model the linear 2D advection equation.",//, &
             " The Non-dimensional Free-stream conditions are",/, &
             " Density          = ",es14.6,/, &
             " Pressure         = ",es14.6,/, &
             " Temperature      = ",es14.6,/, &
             " Mach Number      = ",es14.6,/, &
             " Speed of Sound   = ",es14.6,/, &
             " X-Velocity       = ",es14.6,/, &
             " Y-Velocity       = ",es14.6,/, &
             " Alpha A.o.Attack = ",es14.6," degrees"/, &
             " Beta  A.o.Attack = ",es14.6," degrees"/, &
             " Gamma            = ",es14.6,/, &
             " Gas Constant     = ",es14.6,/)
  !
end subroutine old_output_freestream_conditions
!
!###############################################################################
!
function initial_conditions(xyz,init_cond) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : ncell
  use ovar, only : bc_in
  use ovar, only : itestcase
  use ovar, only : mms_opt,mms_init
  use ovar, only : postproc_debug_grid
  use ovar, only : jet_init_func
  use mms_mod, only : mms_solution_at_xyz
  use channel_mod, only : channel_initialization
  use channel_mod, only : pipe_initialization
  use gpu_mod, only : vortex_solution_gpu
  !
  use taylor_couette_mod, only : laminar_taylor_couette_init
  use taylor_couette_mod, only : turbulent_taylor_couette_init
  use taylor_couette_mod, only : freestream_taylor_couette_init
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:,:)
  integer,  intent(in) :: init_cond
  !
  !.. Function Result ..
  real(wp), dimension(1:size(xyz,dim=2),1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer  :: ibc,k,m,np,ierr
  real(wp) :: div,dl,r,r0,nozzle_len
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: adj_xyz(:,:)
  !
continue
  !
  np = size(xyz,dim=2)
  !
  if (any(itestcase == Transport_Problems)) then
    !
    !###########################################################################
    !##### 2D ISENTROPIC VORTEX PROBLEM
    !###########################################################################
    !
    call vortex_solution_gpu(1,np,xyz,return_value)
    !
  else if (itestcase == Taylor_Green_Vortex) then
    !
    !###########################################################################
    !##### TAYLOR-GREEN VORTEX PROBLEM
    !###########################################################################
    !
    adj_xyz = xyz ! F2003 AUTO-REALLOCATION
    !
    if (postproc_debug_grid) then
      dl = real(ncell,kind=wp)**(-one3)
      adj_xyz = PI * (two*dl*xyz - one)
    end if
    !
    do k = 1,np
      return_value(k,:) = taylor_green_initialization( adj_xyz(:,k) )
    end do
    !
  else if (itestcase == Channel_Flow) then
    !
    !###########################################################################
    !##### CHANNEL FLOW
    !###########################################################################
    !
    return_value(:,:) = channel_initialization( xyz(:,:) )
    !
  else if (itestcase == Pipe_Channel_Flow) then
    !
    !###########################################################################
    !##### PIPE CHANNEL FLOW
    !###########################################################################
    !
    return_value(:,:) = channel_initialization( xyz(:,:) )
   !return_value(:,:) = pipe_initialization( xyz(:,:) )
    !
  else if (abs(itestcase) == Sod_Shock_Tube ) then
    !
    !###########################################################################
    !##### 1D SOD SHOCK TUBE PROBLEM
    !###########################################################################
    !
    do k = 1,np
      return_value(k,:) = sod_init( xyz(:,k) )
    end do
    !
  else if (abs(itestcase) == Density_Pulse_Test ) then
    !
    !###########################################################################
    !##### !D PERIODIC DENSITY PULSE
    !###########################################################################
    !
    do k = 1,np
      return_value(k,:) = pulse_test_init( xyz(:,k) )
    end do
    !
  else if (abs(itestcase) == Taylor_Couette_Flow) then
    !
    !###########################################################################
    !##### TAYLOR-COUETTE FLOW TEST CASE
    !###########################################################################
    !
    ! These initializations for the Taylor-Couette flow come from
    ! the 2021 Advances in Aerodynamics journal article by ZJ Wang:
    !      Benchmark for scale-resolving simulation
    !      with curved walls: the Taylor Couette flow
    !
    ! If itestcase == -Taylor_Couette_Flow
    !  - initialize to the laminar solution
    ! Else if itestcase == +Taylor_Couette_Flow
    !  - initialize to the Taylor-Green-inspired initial solution
    !
    if (itestcase < 0) then
      do k = 1,np
        return_value(k,:) = laminar_taylor_couette_init( xyz(:,k) )
      end do
    else
      do k = 1,np
        return_value(k,:) = turbulent_taylor_couette_init( xyz(:,k) )
      end do
    end if
    !
  else if (abs(itestcase) == Freestream_TCF) then
    !
    !###########################################################################
    !##### TAYLOR-COUETTE FLOW FREESTREAM PRESERVATION
    !###########################################################################
    !
    ! Initialize for freestream preservation using Taylor Couette flow geometry
    ! 
    do k = 1,np
      return_value(k,:) = freestream_taylor_couette_init( xyz(:,k) )
    end do
    !
  else if (mms_opt /= 0) then
    !
    !###########################################################################
    !##### METHOD OF MANUFACTURED SOLUTIONS
    !###########################################################################
    !
    do k = 1,np
      return_value(k,:) = mms_init * mms_solution_at_xyz( xyz(:,k) , &
                                                          return_nondim=true, &
                                                          return_cv=true )
    end do
    !
  else if (abs(itestcase) == Subsonic_Jet) then
    !
    !###########################################################################
    !##### SUBSONIC JETS
    !###########################################################################
    !
    if (init_cond /= 0) then
      !
      do k = 1,np
        return_value(k,:) = subsonic_jet_profile_using_astar(xyz(:,k), &
                                                         bc_in(0)%pv(:), &
                                                         bc_in(init_cond)%pv(:))
       !return_value(k,:) = subsonic_jet_profile(xyz(:,k),bc_in(0)%pv(:), &
       !                                         bc_in(init_cond)%pv(:))
      end do
      !
    else
      !
      do k = 1,np
        return_value(k,:) = bc_in(init_cond)%cv(:)
      end do
      !
    end if
    !
  else if (abs(itestcase) == Supersonic_Jet) then
    !
    !###########################################################################
    !##### SUPERSONIC JET
    !###########################################################################
    !
    if (init_cond /= 0) then
      !
      do k = 1,np
        return_value(k,:) = ssjet_profile(xyz(:,k),bc_in(0)%pv(:), &
                                          bc_in(init_cond)%pv(:))
      end do
      !
    else
      !
      do k = 1,np
        return_value(k,:) = bc_in(init_cond)%cv(:)
      end do
      !
    end if
    !
  else if (itestcase == Infinite_Cylinder) then
    !
    !###########################################################################
    !##### INFINITE CYLINDER TEST CASE
    !###########################################################################
    !
    do k = 1,np
      return_value(k,:) = initialize_cylinder(xyz(:,k),bc_in(0)%pv(:))
    end do
    !
  else if (any(abs(itestcase) == Dual_Stream_Mixing_Flows)) then
    !
    !###########################################################################
    !##### DUAL STREAM MIXING FLOW
    !###########################################################################
    !
    do k = 1,np
      ibc = merge( 2 , 1 , xyz(2,k) > zero )
      do m = 1,nq
        return_value(k,m) = bc_in(ibc)%cv(m)
      end do
    end do
    !
  else if (abs(itestcase) == CH4_3D_Diffusion_Flame) then
    !
    !###########################################################################
    !##### METHANE DIFFUSION FLAME TEST CASE
    !###########################################################################
    !
    if (jet_init_func == 0) then
      do k = 1,np
        return_value(k,:) = bc_in(1)%cv(:)
      end do
    else if (jet_init_func == 1) then
      do k = 1,np
        if (xyz(1,k) >= zero) then
          return_value(k,:) = bc_in(1)%cv(:)
        else
          ibc = merge( 2 , 1 , norm2(xyz(2:,k)) <= 0.03157_wp )
          do m = 1,nq
            return_value(k,m) = bc_in(ibc)%cv(m)
          end do
        end if
      end do
    else if (jet_init_func == 2) then
      do k = 1,np
        ibc = merge( 2 , 1 , norm2(xyz(2:,k)) <= 0.03157_wp )
        do m = 1,nq
          return_value(k,m) = bc_in(ibc)%cv(m)
        end do
      end do
    end if
    !
  else if (abs(itestcase) == CH4_Flame_Speed_Test) then
    !
    !###########################################################################
    !##### METHANE FLAME SPEED TEST
    !###########################################################################
    !
    do k = 1,np
      ibc = merge( 2 , 1 , xyz(1,k) >= zero )
      do m = 1,nq
        return_value(k,m) = bc_in(ibc)%cv(m)
      end do
    end do
    !
  else
    !
    !###########################################################################
    !##### GENERIC INITIALIZATION USING FREESTREAM OR A BC
    !###########################################################################
    !
    ibc = init_cond
    !
#ifndef PBS_ENV
    if (jet_init_func >= lbound(bc_in,dim=1)) then
      if (jet_init_func <= ubound(bc_in,dim=1)) then
        if (allocated(bc_in(jet_init_func)%cv)) then
          ibc = jet_init_func
        end if
      end if
    end if
#endif
    !
    do k = 1,np
      return_value(k,:) = bc_in(init_cond)%cv(:)
    end do
    !
  end if
  !
end function initial_conditions
!
!###############################################################################
!
pure function pv_to_cv(pv) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: pv(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:nq)
  !
  !.. Local Scalars ..
  real(wp) :: r,p,q2
  !
continue
  !
  r = pv(nec)
  p = pv(nee)
  !
  q2 = dot_product( pv(nmb:nme) , pv(nmb:nme) )
  !
  return_value(nec) = r
  return_value(nmb:nme) = r*pv(nmb:nme)
  return_value(nee) = p/(gam-one) + half*r*q2
  !
end function pv_to_cv
!
!###############################################################################
!
pure &
function ssjet_profile(xyz,fs_pv,bc_pv) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd,aref
  use ovar, only : jet_init_func
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:)
  real(wp), intent(in) :: fs_pv(:)
  real(wp), intent(in) :: bc_pv(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(fs_pv))
  !
  !.. Local Scalars ..
  integer  :: ierr
  real(wp) :: x,radius,nozzle_radius
  real(wp) :: gam_over_gm1,total_cond
  real(wp) :: p0,t0,asoa,mach0,mach1
  real(wp) :: p_bc,t_bc,aspd_bc,mach_bc
  real(wp) :: p,t,aspd,mach,blend
  !
  !.. Local Arrays ..
  real(wp) :: pv(1:size(fs_pv))
  real(wp) :: pv_init(1:size(fs_pv))
  !
  !.. Local Parameters ..
  ! This is the axial location that the axial blending begins which
  ! jet_init_func = 2 or 3.
  real(wp), parameter :: start_xblend  = three*ten*inches_to_meters
  ! Axial location of the nozzle throat
  real(wp), parameter :: throat_xloc   = -2.54585414585415_wp*inches_to_meters
  ! Radius of the nozzle at the throat
  real(wp), parameter :: throat_radius =  0.87762738764425_wp*inches_to_meters
  !
continue
  !
  gam_over_gm1 = gam/(gam-one)
  !
  ! Get the nozzle radius at this axial location and the radius at this point
  !
  x = xyz(1)
  radius = norm2(xyz(2:))
  nozzle_radius = ssjet_interior_wall_radius(x)
  !
  ! Compute the total pressure and total temperature for the BC
  !
  p_bc = bc_pv(nee)
  t_bc = p_bc / (rgasref_nd*bc_pv(nec))
  !
  aspd_bc = sqrt(gam*rgasref_nd*t_bc)
  mach_bc = norm2(bc_pv(nmb:nme))/aspd_bc
  !
  total_cond = one + half*(gam-one)*mach_bc*mach_bc
  !
  p0 = p_bc * (total_cond**gam_over_gm1)
  t0 = t_bc * total_cond
  !
  ! Compute A*/A
  !
  asoa = throat_radius/nozzle_radius
  !
  ! Initial guess for computing mach number based on location up or
  ! downstream of the nozzle throat (subsonic or supersonic)
  !
  if (x <= throat_xloc) then
    mach0 = eps15
    mach1 = one
  else
    mach0 = one
    mach1 = ten
  end if
  !
  ! Get Mach number from A*/A using
  !
  call get_mach_from_astar_over_a(mach0,mach1,asoa,mach,ierr)
  !
  ! Compute the primitive variables based on the A*/A mach number
  ! and BC total conditions
  !
  total_cond = one + half*(gam-one)*mach*mach
  !
  t = t0 / total_cond
  p = p0 / (total_cond**gam_over_gm1)
  !
  aspd = sqrt(gam*rgasref_nd*t)
  !
  pv(nec) = p/(rgasref_nd*t)
  pv(nmx) = mach*aspd
  pv(nmx+1:nme) = zero
  pv(nee) = p
  !
  ! Compute the blending coefficient between the freestream and BC states
  !
  blend = ssjet_blending(jet_init_func,x,start_xblend,radius,nozzle_radius)
  !
  ! Compute the primitive variables for the blended initial conditions
  ! 
  pv_init(nec) = blend*pv(nec) + (one-blend)*fs_pv(nec)
  pv_init(nmx) = blend*pv(nmx) + (one-blend)*fs_pv(nmx)
  pv_init(nmx+1:nme) = zero
 !pv_init(nmy) = blend/aref
  pv_init(nee) = fs_pv(nee)
  !
  ! Convert the primitive variables to conserved variables
  !
  return_value(:) = pv_to_cv(pv_init)
  !
end function ssjet_profile
!
!###############################################################################
!
pure subroutine get_mach_from_astar_over_a(mach0,mach1,asoa,machc,ierr)
  !
  ! Solves for mach number given A*/A using the bisection method
  ! initial guess is either subsonic (M0=0/M1=1) or supersonic (M0=1/M1=10)
  !
  ! bisection method
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: mach0
  real(wp),    intent(in) :: mach1
  real(wp),    intent(in) :: asoa
  real(wp), intent(inout) :: machc
  integer,  intent(inout) :: ierr
  !
  !.. Local Scalars ..
  integer  :: n,siga,sigb,sigc
  real(wp) :: fa,fb,fc,res,macha,machb
  !
  !.. Local Parameters ..
  integer,  parameter :: maxcount = 1000    ! max iterations
  real(wp), parameter :: c0 = (six/five)**3 ! constant for computing fa,fb,fc
  real(wp), parameter :: res_limit = eps8   ! iteration goal
  !
continue
  !
  macha = mach0
  machb = mach1
  !
  ! initial guess
  !
  res = huge(zero)
  !
  ierr = 1
  bisec: &
  do n = 1,maxcount
    !
    machc = (macha+machb)*half
    !
    fa = c0*macha*( (one+macha*macha/five)**(-3) ) - asoa
    fb = c0*machb*( (one+machb*machb/five)**(-3) ) - asoa
    fc = c0*machc*( (one+machc*machc/five)**(-3) ) - asoa
    !
    siga = int((fa)/abs(fa))
    sigb = int((fb)/abs(fb))
    sigc = int((fc)/abs(fc))
    !
    if (siga /= sigc) then
      res = abs(machc-macha)
      machb = machc
    else if (sigb /= sigc) then
      res = abs(machc-machb)
      macha = machc
    else
     !write(*,*) 'oops',siga,sigb,sigc
    end if
    !
    if (res <= res_limit) then
      ierr = 0
      exit bisec
    end if
    !
  end do bisec
  !
end subroutine get_mach_from_astar_over_a
!
!###############################################################################
!
pure function ssjet_blending(init_func,x,x0,r,r0) result(return_value)
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: init_func
  real(wp), intent(in) :: x,x0
  real(wp), intent(in) :: r,r0
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: xblend,xref,xdiv
  real(wp) :: rblend,rref,rdiv
  real(wp) :: val,div_exp,a,b,val_exp,val_spread,val_ref,ref_exp
  !
continue
  !
  xblend = one
  rblend = one
  !
  if (x < zero) then
    !
    ! This is upstream of the nozzle exit, so set the blending to 1
    ! so that it uses the flow conditions computed from the nozzle
    ! inflow BC.
    !
    xblend = one
    rblend = one
    !
  else
    !
    ! This is downstream of the nozzle exit
    !
    if (init_func == 1) then
      !
      ! Use a constant radial blending only
      !
      xblend = one
      !
      rref = (1.02_wp)*r0
      val_exp = one/four
      val = max( zero , min( x/(x0/(ten)) , one ) )**val_exp
      div_exp = three/two
      a = (one-val)**div_exp
      b = val**div_exp
      rdiv = rref/(ten*five*a + four*b)
      rblend = erf_blend(half,r,rref,rdiv)
      !
    else if (init_func == 2) then
      !
      ! Use a blending in both x and r with a constant radial spreading
      !
      xref = x0
      xdiv = half*xref
      xblend = erf_blend(half,x,xref,xdiv)
      !
      rref = (1.02_wp)*r0
      val_exp = one/four
      val = max( zero , min( x/(x0/(ten)) , one ) )**val_exp
      div_exp = three/two
      a = (one-val)**div_exp
      b = val**div_exp
      rdiv = rref/(ten*five*a + four*b)
      rblend = erf_blend(half,r,rref,rdiv)
      !
    else if (init_func == 3) then
      !
      ! Use a blending in both x and r with a downstream radial spreading,
      ! i.e., variable rref and rdiv
      !
      xref = x0
      xdiv = half*xref
      xblend = erf_blend(half,x,xref,xdiv)
      !
     !rref = r0  + x/xref
     !rdiv = half*rref + (2.5_wp*x)/xref
      !
      val_ref = one/(two*ten)
      ref_exp = one
      val_spread = 1.5_wp*0.0254_wp/one
     !val_spread = zero
      rref = (1.02_wp)*r0 + val_ref*(x/xref)**ref_exp
      val_exp = one/four
      val = max( zero , min( x/(x0/(ten)) , one ) )**val_exp
      div_exp = three/two
      a = (one-val)**div_exp
      b = val**div_exp
      rdiv = rref/(ten*five*a + four*b) + (val_spread*x)/xref
      rblend = erf_blend(half,r,rref,rdiv)
      !
    else
      !
      ! This will set the initial conditions to the free stream conditions
      !
      xblend = zero
      rblend = zero
      !
    end if
    !
  end if
  !
  return_value = xblend*rblend
  !
end function ssjet_blending
!
!###############################################################################
!
pure function erf_blend(value_at_zero,x,xref,xdiv) result(return_value)
  !
  ! Function to compute a blending coefficient using the error function
  !
  real(wp), intent(in) :: value_at_zero ! value to take for erf(0)
  real(wp), intent(in) :: x
  real(wp), intent(in) :: xref
  real(wp), intent(in) :: xdiv
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
continue
  !
  return_value = value_at_zero*( one - erf((x-xref)/xdiv) )
  !
end function erf_blend
!
!###############################################################################
!
pure &
function compute_astar_radius(radius,mach) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  real(wp), intent(in) :: radius
  real(wp), intent(in) :: mach
  !
  real(wp) :: return_value
  !
  real(wp) :: area,gp1,gm1,gexp,b,c,astar
  !
continue
  !
  gp1 = gam+one
  gm1 = gam-one
  gexp = half*gp1/gm1
  !
  area = radius*radius*pi
  !
  b = (half*gp1)**(-gexp)
  c = (one + half*gm1*mach*mach)**(gexp)
  !
  astar = area*mach/(b*c)
  return_value = sqrt(astar/pi)
  !
end function compute_astar_radius
!
!###############################################################################
!
pure &
function subsonic_jet_profile_using_astar(xyz,fs_pv,bc_pv) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd,aref
  use ovar, only : jet_init_func,itestcase
  use ovar, only : grid_scaling_factor
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:)
  real(wp), intent(in) :: fs_pv(:)
  real(wp), intent(in) :: bc_pv(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(fs_pv))
  !
  !.. Local Scalars ..
  integer  :: ierr
  real(wp) :: x,radius,nozzle_radius,xbeg,rbeg
  real(wp) :: gam_over_gm1,two_over_gm1,gm1_over_gam
  real(wp) :: p,t,aspd,mach,p0,t0,total_cond
  real(wp) :: p_bc,t_bc,aspd_bc,mach_bc,p0_bc,t0_bc
  real(wp) :: p_fs,t_fs,aspd_fs,mach_fs,p0_fs,t0_fs
  real(wp) :: mach_jet,p0_jet,t0_jet,astar,asoa
  real(wp) :: theta,aoa,vx,vy,vz,vr,vmag
  real(wp) :: div,rfunc,xfunc,blend
  !
  !.. Local Arrays ..
  real(wp) :: pv_init(1:size(fs_pv))
  !
  !.. Local Parameters ..
  ! This is the axial location that the axial blending begins which
  ! jet_init_func = 2 or 3.
  real(wp), parameter :: start_xblend  = three*ten*inches_to_meters
  !
  ! This is a fudge factor for making sure a point is inside the nozzle
  ! interior. The lip is 0.03 inches thick, so any radius greater than the
  ! (nozzle radius at a given axial location + 0.025 inches) is definitely
  ! outside the nozzle.
  real(wp), parameter :: radius_eps = 0.025_wp*inches_to_meters
  !
  ! This controls the radial blending of the mach number inside the nozzle.
  ! It is used to approximate a boundary layer near the surface, with a
  ! smaller number creating a thinner boundary layer.
  real(wp), parameter :: radius_bl = 0.05_wp
  !
  real(wp), parameter :: vr_div = three/five
  !
continue
  !
  gam_over_gm1 = gam/(gam-one)
  two_over_gm1 = two/(gam-one)
  gm1_over_gam = one/gam_over_gm1
  !
  ! Get the nozzle radius at this axial location and the radius at this point
  !
  x = xyz(1)
  radius = norm2(xyz(2:))
  if (itestcase < 0) then
    nozzle_radius = arn2_interior_wall_radius(x)
    xbeg = 7.74_wp*grid_scaling_factor
    rbeg = arn2_interior_wall_radius(xbeg)
  else
    nozzle_radius = smc000_interior_wall_radius(x)
    xbeg = 9.875_wp*grid_scaling_factor
    rbeg = smc000_interior_wall_radius(xbeg)
  end if
  !
  theta = atan2(xyz(3),xyz(2))
  !
  ! Compute the Mach number, P0, and T0 for the BC
  !
  p_bc = bc_pv(nee)
  t_bc = p_bc / (rgasref_nd*bc_pv(nec))
  !
  aspd_bc = sqrt(gam*rgasref_nd*t_bc)
  mach_bc = norm2(bc_pv(nmb:nme))/aspd_bc
  !
  total_cond = one + half*(gam-one)*mach_bc*mach_bc
  !
  p0_bc = p_bc * (total_cond**gam_over_gm1)
  t0_bc = t_bc * total_cond
  !
  ! Compute the freestream Mach number, P0, and T0
  !
  p_fs = fs_pv(nee)
  t_fs = p_fs / (rgasref_nd*fs_pv(nec))
  !
  aspd_fs = sqrt(gam*rgasref_nd*t_fs)
  mach_fs = norm2(fs_pv(nmb:nme))/aspd_bc
  !
  total_cond = one + half*(gam-one)*mach_fs*mach_fs
  !
  p0_fs = p_fs * (total_cond**gam_over_gm1)
  t0_fs = t_fs * total_cond
  !
  ! Get the Mach number, P0, and T0 for the jet
  ! NOTE: P0 and T0 are assumed constant from the nozzle inflow to the exit
  !
  p0_jet = p0_bc
  t0_jet = t0_bc
  ! Use the freestream pressure to get the ideal jet Mach no. at the nozzle exit
  mach_jet = sqrt( two_over_gm1*( ((p0_jet/p_fs)**gm1_over_gam) - one ) )
  !
  ! Compute A* from the BC mach number and inflow radius
  !
 !astar = compute_astar_radius(rbeg,mach_bc)
  astar = compute_astar_radius(0.0254_wp,mach_jet)
  !
  ! Compute A*/A at this location
  !
  asoa = astar/(nozzle_radius)
 !asoa = astar*astar/(nozzle_radius*nozzle_radius)
  !
  ! Compute the primitive variables for the blended initial conditions
  ! 
  if (x <= zero) then
    if (radius > nozzle_radius+radius_eps) then
      !
      ! I dont think it will get here because the subroutine
      ! find_elem_for_initial_jet_profile should leave this section of the
      ! domain with an initialization cell tag of 0 indicating that it
      ! should not use the subsonic jet initialization and just initialize
      ! these cells to the freestream conditions.
      !
      p0 = p0_fs
      t0 = t0_fs
      mach = mach_fs
      !
      aoa = zero
      !
    else
      !
      call get_mach_from_astar_over_a(eps15,one,asoa,mach,ierr)
     !call get_mach_from_astar_over_a(mach_bc,one,asoa,mach,ierr)
      !
      div = radius_bl*nozzle_radius
      blend = half * (one + erf((four/xbeg)*abs(x)-two))
      !
      t0   = blend*t0_bc   + (one-blend)*t0_jet
      p0   = blend*p0_bc   + (one-blend)*p0_jet
     !mach = blend*mach_bc + (one-blend)*mach_jet
      !
      ! Add a radial blending on the Mach number to approximate a boundary layer
      !
      rfunc = half * (one - erf((radius-nozzle_radius)/div))
      mach = mach*rfunc
      !
      div = vr_div*nozzle_radius
      rfunc = half * (one - erf((radius-nozzle_radius)/div))
      aoa = atan2( arn2_interior_wall_slope(x) , one ) * (one-rfunc)
      !
    end if
  else
    !
    blend = ssjet_blending(jet_init_func,x,start_xblend,radius,nozzle_radius)
    !
    t0   = blend*t0_jet   + (one-blend)*t0_fs
    p0   = blend*p0_jet   + (one-blend)*p0_fs
    mach = blend*mach_jet + (one-blend)*mach_fs
    !
    aoa = zero
    !
  end if
  !
  total_cond = one + half*(gam-one)*mach*mach
  !
  t = t0 / total_cond
  p = p0 / (total_cond**gam_over_gm1)
  !
  vmag = mach*sqrt(gam*rgasref_nd*t)
  !
  vx = vmag*cos(aoa)
  vr = vmag*sin(aoa)
  !
  vy = vr*cos(theta)
  vz = vr*sin(theta)
  !
  pv_init(nec) = p / (rgasref_nd*t)
  pv_init(nmx) = vx
 !pv_init(nmy) = astar/aref
  pv_init(nmy) = vy
 !pv_init(nmy) = real(ierr,kind=wp)/aref
  pv_init(nmz) = vz
 !pv_init(nmz) = asoa/aref
  pv_init(nee) = p
  !
  ! Convert the primitive variables to conserved variables
  !
  return_value(:) = pv_to_cv(pv_init)
  !
end function subsonic_jet_profile_using_astar
!
!###############################################################################
!
pure &
function subsonic_jet_profile(xyz,fs_pv,bc_pv) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam,rgasref_nd
  use ovar, only : jet_init_func,itestcase
  use ovar, only : grid_scaling_factor
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:)
  real(wp), intent(in) :: fs_pv(:)
  real(wp), intent(in) :: bc_pv(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(fs_pv))
  !
  !.. Local Scalars ..
  integer  :: ierr
  real(wp) :: x,radius,nozzle_radius,xbeg
  real(wp) :: gam_over_gm1,two_over_gm1,gm1_over_gam
  real(wp) :: p,t,aspd,mach,p0,t0,total_cond
  real(wp) :: p_bc,t_bc,aspd_bc,mach_bc,p0_bc,t0_bc
  real(wp) :: p_fs,t_fs,aspd_fs,mach_fs,p0_fs,t0_fs
  real(wp) :: mach_jet,p0_jet,t0_jet
  real(wp) :: div,rfunc,xfunc,blend
  !
  !.. Local Arrays ..
  real(wp) :: pv_init(1:size(fs_pv))
  !
  !.. Local Parameters ..
  ! This is the axial location that the axial blending begins which
  ! jet_init_func = 2 or 3.
  real(wp), parameter :: start_xblend  = three*ten*inches_to_meters
  !
  ! This is a fudge factor for making sure a point is inside the nozzle
  ! interior. The lip is 0.03 inches thick, so any radius greater than the
  ! (nozzle radius at a given axial location + 0.025 inches) is definitely
  ! outside the nozzle.
  real(wp), parameter :: radius_eps = 0.025_wp*inches_to_meters
  !
  ! This controls the radial blending of the mach number inside the nozzle.
  ! It is used to approximate a boundary layer near the surface, with a
  ! smaller number creating a thinner boundary layer.
  real(wp), parameter :: radius_bl = 0.05_wp
  !
continue
  !
  gam_over_gm1 = gam/(gam-one)
  two_over_gm1 = two/(gam-one)
  gm1_over_gam = one/gam_over_gm1
  !
  ! Get the nozzle radius at this axial location and the radius at this point
  !
  x = xyz(1)
  radius = norm2(xyz(2:))
  if (itestcase < 0) then
    nozzle_radius = arn2_interior_wall_radius(x)
    xbeg = 7.74_wp*grid_scaling_factor
  else
    nozzle_radius = smc000_interior_wall_radius(x)
    xbeg = 9.875_wp*grid_scaling_factor
  end if
  !
  ! Compute the Mach number, P0, and T0 for the BC
  !
  p_bc = bc_pv(nee)
  t_bc = p_bc / (rgasref_nd*bc_pv(nec))
  !
  aspd_bc = sqrt(gam*rgasref_nd*t_bc)
  mach_bc = norm2(bc_pv(nmb:nme))/aspd_bc
  !
  total_cond = one + half*(gam-one)*mach_bc*mach_bc
  !
  p0_bc = p_bc * (total_cond**gam_over_gm1)
  t0_bc = t_bc * total_cond
  !
  ! Compute the freestream Mach number, P0, and T0
  !
  p_fs = fs_pv(nee)
  t_fs = p_fs / (rgasref_nd*fs_pv(nec))
  !
  aspd_fs = sqrt(gam*rgasref_nd*t_fs)
  mach_fs = norm2(fs_pv(nmb:nme))/aspd_bc
  !
  total_cond = one + half*(gam-one)*mach_fs*mach_fs
  !
  p0_fs = p_fs * (total_cond**gam_over_gm1)
  t0_fs = t_fs * total_cond
  !
  ! Get the Mach number, P0, and T0 for the jet
  ! NOTE: P0 and T0 are assumed constant from the nozzle inflow to the exit
  !
  p0_jet = p0_bc
  t0_jet = t0_bc
  ! Use the freestream pressure to get the ideal jet Mach no. at the nozzle exit
  mach_jet = sqrt( two_over_gm1*( ((p0_jet/p_fs)**gm1_over_gam) - one ) )
  !
  ! Compute the primitive variables for the blended initial conditions
  ! 
  if (x <= zero) then
    if (radius > nozzle_radius+radius_eps) then
      !
      ! I dont think it will get here because the subroutine
      ! find_elem_for_initial_jet_profile should leave this section of the
      ! domain with an initialization cell tag of 0 indicating that it
      ! should not use the subsonic jet initialization and just initialize
      ! these cells to the freestream conditions.
      !
      p0 = p0_fs
      t0 = t0_fs
      mach = mach_fs
      !
    else
      !
      div = radius_bl*nozzle_radius
      blend = half * (one + erf((four/xbeg)*abs(x)-two))
      !
      t0   = blend*t0_bc   + (one-blend)*t0_jet
      p0   = blend*p0_bc   + (one-blend)*p0_jet
      mach = blend*mach_bc + (one-blend)*mach_jet
      !
      ! Add a radial blending on the Mach number to approximate a boundary layer
      !
      rfunc = half * (one - erf((radius-nozzle_radius)/div))
      mach = mach*rfunc
      !
    end if
  else
    !
    blend = ssjet_blending(jet_init_func,x,start_xblend,radius,nozzle_radius)
    !
    t0   = blend*t0_jet   + (one-blend)*t0_fs
    p0   = blend*p0_jet   + (one-blend)*p0_fs
    mach = blend*mach_jet + (one-blend)*mach_fs
    !
  end if
  !
  total_cond = one + half*(gam-one)*mach*mach
  !
  t = t0 / total_cond
  p = p0 / (total_cond**gam_over_gm1)
  !
  pv_init(nec) = p / (rgasref_nd*t)
  pv_init(nmx) = mach*sqrt(gam*rgasref_nd*t)
  pv_init(nmx+1:nme) = zero
  pv_init(nee) = p
  !
  ! Convert the primitive variables to conserved variables
  !
  return_value(:) = pv_to_cv(pv_init)
  !
end function subsonic_jet_profile
!
!###############################################################################
!
subroutine find_elem_for_initial_jet_profile(point_is_in_nozzle_interior)
  !
  !.. Use Statements ..
  use geovar, only : grid
  use ovar,   only : itestcase
  use ovar,   only : jet_init_func
  use ovar,   only : bc_in,bc_input
  !
  !.. Dummy Function Interface ..
  interface
    function point_is_in_nozzle_interior(xyz) result(return_value)
      import :: wp,lk
      real(wp), intent(in) :: xyz(:)
      logical(lk) :: return_value
    end function point_is_in_nozzle_interior
  end interface
  !
  !.. Local Scalars ..
  integer :: l,n,init,bc_idx
  integer :: nozzle_inflow_bc_idx
  !
  !.. Local Arrays ..
  real(wp) :: xyz_ave(1:3)
  !
  !.. Local Allocatable Arrays ..
  integer,          allocatable :: pts(:)
  character(len=:), allocatable :: bc_input_name
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_elem_for_initial_jet_profile"
  !
 !integer, parameter :: init_val = 1 ! radial error function profile
 !integer, parameter :: init_val = 2 ! inflow mach number limiter
 !integer, parameter :: init_val = 3 ! options 1 and 2 combined
  !
continue
  !
  ! Find the index of bc_in for the nozzle inflow boundary condition
  !
  nozzle_inflow_bc_idx = 1
  l = 0
  !
  do n = lbound(bc_in,dim=1),ubound(bc_in,dim=1)
    bc_idx = bc_in(n)%bc_input_idx
    if (bc_idx >= lbound(bc_input,dim=1) .and. &
        bc_idx <= ubound(bc_input,dim=1)) then
      ! F2003 AUTO-REALLOC
      bc_input_name = uppercase( trim(adjustl(bc_input(bc_idx)%name)) )
      if (bc_input_name == "NOZZLEINFLOW") l = n
      if (bc_input_name == "NOZZLE_INFLOW") l = n
      if (bc_input_name == "NOZZLE INFLOW") l = n
      if (bc_input_name == "INFLOWNOZZLE") l = n
      if (bc_input_name == "INFLOW_NOZZLE") l = n
      if (bc_input_name == "INFLOW NOZZLE") l = n
      if (bc_input_name == "INFLOWSUBSONIC") l = n
      if (l > 0) then
        nozzle_inflow_bc_idx = l
        exit
      end if
    end if
  end do
  !
  xyz_ave(:) = zero
  !
  ! Loop through all the grid cells and update the initial condition tag
  ! depending on whether or not the cell is inside or outside the nozzle
  !
  do n = 1,size(grid%elem)
    !
    if (grid%elem(n)%tags(1) == 0) then
      !
      pts = grid%elem(n)%pts(:) ! USING F2003 AUTO-REALLOCATION
      !
      ! Get the coordinates for the cell center by averaging
      ! all the coordinates of the cells grid points
      !
      do l = 1,size(grid%xyz,dim=1)
        xyz_ave(l) = sum(grid%xyz(l,pts(:))) / real(size(pts),kind=wp)
      end do
      !
      if (xyz_ave(1) > zero) then
        !
        ! If the axial location of the cell center is downstream of the
        ! nozzle exit (x > 0), mark the initial condition depending on
        ! the test case being run and the jet initialization function
        ! being used
        !
        init = nozzle_inflow_bc_idx
        !
      else
        !
        ! If the axial location of the cell center is upstream of the
        ! nozzle exit, use the point_is_in_nozzle_interior to determine
        ! if this cell is inside or outside the nozzle, and set the initial
        ! condition if this cell is inside the nozzle.
        !
       !init = merge(jet_init_func,0,point_is_in_nozzle_interior(xyz_ave))
        init = merge(nozzle_inflow_bc_idx,0, &
                     point_is_in_nozzle_interior(xyz_ave))
        !
      end if
      !
      ! Update the cell tags with the new initial condition
      !
      if (size(grid%elem(n)%tags) < 2) then
        grid%elem(n)%tags = [grid%elem(n)%tags(1),init] ! F2003 AUTO-REALLOC
      else
        grid%elem(n)%tags(2) = init
      end if
      !
    end if
    !
  end do
  !
end subroutine find_elem_for_initial_jet_profile
!
!###############################################################################
!
pure &
function initialize_cylinder(xyz,fs_pv) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : gam
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:)
  real(wp), intent(in) :: fs_pv(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(fs_pv))
  !
  !.. Local Scalars ..
  real(wp) :: a,b,c
  real(wp) :: rho,vx,vy,vz,q2,pres
  real(wp) :: r,theta,uinf
  !
  !.. Local Parameters ..
  real(wp), parameter :: sig = one/ten
  real(wp), parameter :: three2 = three/two
  real(wp), parameter :: five4 = five/four
  real(wp), parameter :: halfpi = half*pi
  !
continue
  !
  r = norm2(xyz(1:2))
  theta = acos(xyz(1)/r)
  !
  rho  = fs_pv(nec)
  uinf = fs_pv(nmx)
  vx   = uinf
  vy   = fs_pv(nmy)
  vz   = fs_pv(nmz)
  pres = fs_pv(nee)
  !
  if (r < three2) then
    !
    vx = uinf * cylinder_function(r,one,half)
    !
    a = sin( halfpi * xyz(3) )
    !
    if (r > half .and. r <= one) then
      b = cylinder_function(r,three4,one4)
    else
      b = one - cylinder_function(r,five4,one4)
    end if
    !
    if (theta > zero .and. theta < pi) then
      c = cylinder_function(theta,halfpi,halfpi)
    else
      c = one - cylinder_function(theta,three*halfpi,halfpi)
    end if
    !
    vy = sig * uinf * b * c
    vz = sig * uinf * b * a*a
    !
  end if
  !
  q2 = vx*vx + vy*vy + vz*vz
  !
  return_value(nec) = rho
  return_value(nmx) = rho*vx
  return_value(nmy) = rho*vy
  return_value(nmz) = rho*vz
  return_value(nee) = pres/(gam-one) + half*rho*q2
  !
end function initialize_cylinder
!
!###############################################################################
!
pure function cylinder_function(r,rm,ra) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: r,rm,ra
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  real(wp) :: r_min_rm,numer,denom
  !
continue
  !
  r_min_rm = r - rm
  !
  numer = two * ra * r_min_rm
  denom = ra*ra - r_min_rm*r_min_rm
  !
  return_value = half * (one + tanh(numer/denom))
  !
end function cylinder_function
!
!###############################################################################
!
subroutine detach_lip_flame(xyz,cv,gp)
  !
  !.. Use Statements ..
  use ovar,       only : tref
  use thermo_mod, only : compute_mixture_enthalpy
  use thermo_mod, only : compute_mixture_rgas
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: xyz(:,:)
  real(wp), intent(inout) :: cv(:,:)
  real(wp), intent(inout) :: gp(:,:)
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: x,y,z,r,r1,r2
  real(wp) :: rho,u,v,w,t,p,rhoe,h,rmix
  real(wp) :: t_limit
  !
  !.. Local Arrays ..
  real(wp) :: mass_frac(1:ncs)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "detach_lip_flame"
  !
continue
  !
  t_limit = 750.0_wp/tref
  !
  do n = 1,size(xyz,dim=2)
    !
    x = xyz(1,n)
    y = xyz(2,n)
    z = xyz(3,n)
    !
  ! if (x < -0.0003_wp .or. x > 0.017_wp) cycle
    if (x < -0.0003_wp .or. x > 0.000983_wp) cycle
    !
    r = sqrt(y*y + z*z)
    !
    r1 = 0.031253_wp
    r2 = 0.0323973_wp
  ! r1 = 0.029_wp
  ! r2 = 0.034_wp
  !!r1 = -0.1512_wp*x + 0.03157_wp
  !!r2 = one4*x + 0.03175_wp
    !
    if (r >= r1 .and. r <= r2) then
      !
      if (gp(n,ntemp) > t_limit) then
        !
        rho = cv(n,nec)
        u = cv(n,nmx)/rho
        v = cv(n,nmy)/rho
        w = cv(n,nmz)/rho
        t = t_limit
        p = gp(n,npres)
        !
        mass_frac(1:ncs) = cv(n,1:ncs)/rho
        !
        call compute_mixture_enthalpy(mass_frac,t,h)
        !
        call compute_mixture_rgas(mass_frac,rmix)
        !
        rho = p / (rmix*t)
        !
        rhoe = rho*h - p + half*rho*(u*u + v*v + w*w)
        !
        cv(n,1:ncs) = rho*mass_frac(1:ncs)
        !
        cv(n,nec) = rho
        cv(n,nmx) = rho*u
        cv(n,nmy) = rho*v
        cv(n,nmz) = rho*w
        cv(n,nee) = rhoe
        !
        gp(n,ntemp) = t
        gp(n,nenth) = h
        !
      end if
      !
    end if
    !
  end do
  !
end subroutine detach_lip_flame
!
!###############################################################################
!
subroutine initialize_dusp()
  !
  !.. Use Statements ..
  use geovar,  only : nr
  use flowvar, only : dusp
  use ovar,    only : force_nd_subroutines
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_dusp"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (allocated(dusp)) then
    !
    if (force_nd_subroutines) then
      !
      call initialize_dusp_nd(1,nq)
      !
    else
      !
      if (nr == 3) then
        call initialize_dusp_3d
      else if (nr == 2) then
        call initialize_dusp_2d
      else
        call initialize_dusp_1d
      end if
      !
      if (ncs > 0) then
        call initialize_dusp_nd(1,ncs)
      end if
      !
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine initialize_dusp
!
!###############################################################################
!
subroutine initialize_dusp_3d()
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use geovar,          only : ncell,cell
  use metrics_mod,     only : metrics_array
  use derivatives_mod, only : deriv
  use flowvar,         only : usp,dusp
  !
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  use omp_lib, only : omp_get_thread_num
  use omp_lib, only : omp_get_num_threads
#endif
#endif
  !
  !.. Local Scalars ..
  integer  :: i,j,n,nc,n0,np
  integer  :: this_geom,this_order
  real(wp) :: a
  real(wp) :: decdi,dmxdi,dmydi,dmzdi,deedi,didx,didy,didz
  real(wp) :: decdj,dmxdj,dmydj,dmzdj,deedj,djdx,djdy,djdz
  real(wp) :: decdk,dmxdk,dmydk,dmzdk,deedk,dkdx,dkdy,dkdz
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_dusp_3d"
  !
continue
  !
  ! Loop over the cells and get the derivative of
  ! the primitive variables at each solution point
  !
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP PARALLEL &
  !$OMP&         DEFAULT(SHARED) &
  !$OMP&         PRIVATE(tid,nc,n0,np,this_geom,this_order,n,i,j,a, &
  !$OMP&                 decdi,dmxdi,dmydi,dmzdi,deedi, &
  !$OMP&                 decdj,dmxdj,dmydj,dmzdj,deedj, &
  !$OMP&                 decdk,dmxdk,dmydk,dmzdk,deedk, &
  !$OMP&                 didx,djdx,dkdx,didy,djdy,dkdy,didz,djdz,dkdz)
  !
  tid = omp_get_thread_num()
  nthreads = omp_get_num_threads()
  !
  if (mypnum == glb_root) then
    if (tid == glb_root) then
      write (iout,thread_report) pname,nthreads
    end if
  end if
  !
  chunk_size = ceiling( real(ncell,kind=wp) / real(nthreads,kind=wp) )
  chunk_size = chunk_size / 10
  !
  !$OMP DO SCHEDULE(RUNTIME)
 !!$OMP DO SCHEDULE(GUIDED,chunk_size)
 !!$OMP DO SCHEDULE(STATIC,chunk_size)
#endif
#endif
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1  ! beginning index for solpts in cell nc
    np = cell(nc)%end_sp-n0 ! number of solution points in this cell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Compute the gradient of the conservative variables at each solution point
    !
    do n = 1,np
      !
      decdi = zero
      dmxdi = zero
      dmydi = zero
      dmzdi = zero
      deedi = zero
      !
      do i = deriv(this_geom,this_order)%d(1)%row_ptr(n)+1, &
             deriv(this_geom,this_order)%d(1)%row_ptr(n+1)
        !
        j = deriv(this_geom,this_order)%d(1)%row_idx(i)
        a = deriv(this_geom,this_order)%d(1)%row_val(i)
        !
        decdi = decdi + a * usp(n0+j,nec)
        dmxdi = dmxdi + a * usp(n0+j,nmx)
        dmydi = dmydi + a * usp(n0+j,nmy)
        dmzdi = dmzdi + a * usp(n0+j,nmz)
        deedi = deedi + a * usp(n0+j,nee)
        !
      end do
      !
      decdj = zero
      dmxdj = zero
      dmydj = zero
      dmzdj = zero
      deedj = zero
      !
      do i = deriv(this_geom,this_order)%d(2)%row_ptr(n)+1, &
             deriv(this_geom,this_order)%d(2)%row_ptr(n+1)
        !
        j = deriv(this_geom,this_order)%d(2)%row_idx(i)
        a = deriv(this_geom,this_order)%d(2)%row_val(i)
        !
        decdj = decdj + a * usp(n0+j,nec)
        dmxdj = dmxdj + a * usp(n0+j,nmx)
        dmydj = dmydj + a * usp(n0+j,nmy)
        dmzdj = dmzdj + a * usp(n0+j,nmz)
        deedj = deedj + a * usp(n0+j,nee)
        !
      end do
      !
      decdk = zero
      dmxdk = zero
      dmydk = zero
      dmzdk = zero
      deedk = zero
      !
      do i = deriv(this_geom,this_order)%d(3)%row_ptr(n)+1, &
             deriv(this_geom,this_order)%d(3)%row_ptr(n+1)
        !
        j = deriv(this_geom,this_order)%d(3)%row_idx(i)
        a = deriv(this_geom,this_order)%d(3)%row_val(i)
        !
        decdk = decdk + a * usp(n0+j,nec)
        dmxdk = dmxdk + a * usp(n0+j,nmx)
        dmydk = dmydk + a * usp(n0+j,nmy)
        dmzdk = dmzdk + a * usp(n0+j,nmz)
        deedk = deedk + a * usp(n0+j,nee)
        !
      end do
      !
      didx = metrics_array(n0+n,1,1)
      djdx = metrics_array(n0+n,2,1)
      dkdx = metrics_array(n0+n,3,1)
      !
      didy = metrics_array(n0+n,1,2)
      djdy = metrics_array(n0+n,2,2)
      dkdy = metrics_array(n0+n,3,2)
      !
      didz = metrics_array(n0+n,1,3)
      djdz = metrics_array(n0+n,2,3)
      dkdz = metrics_array(n0+n,3,3)
      !
      dusp(n0+n,1,nec) = decdi*didx + decdj*djdx + decdk*dkdx
      dusp(n0+n,2,nec) = decdi*didy + decdj*djdy + decdk*dkdy
      dusp(n0+n,3,nec) = decdi*didz + decdj*djdz + decdk*dkdz
      !
      dusp(n0+n,1,nmx) = dmxdi*didx + dmxdj*djdx + dmxdk*dkdx
      dusp(n0+n,2,nmx) = dmxdi*didy + dmxdj*djdy + dmxdk*dkdy
      dusp(n0+n,3,nmx) = dmxdi*didz + dmxdj*djdz + dmxdk*dkdz
      !
      dusp(n0+n,1,nmy) = dmydi*didx + dmydj*djdx + dmydk*dkdx
      dusp(n0+n,2,nmy) = dmydi*didy + dmydj*djdy + dmydk*dkdy
      dusp(n0+n,3,nmy) = dmydi*didz + dmydj*djdz + dmydk*dkdz
      !
      dusp(n0+n,1,nmz) = dmzdi*didx + dmzdj*djdx + dmzdk*dkdx
      dusp(n0+n,2,nmz) = dmzdi*didy + dmzdj*djdy + dmzdk*dkdy
      dusp(n0+n,3,nmz) = dmzdi*didz + dmzdj*djdz + dmzdk*dkdz
      !
      dusp(n0+n,1,nee) = deedi*didx + deedj*djdx + deedk*dkdx
      dusp(n0+n,2,nee) = deedi*didy + deedj*djdy + deedk*dkdy
      dusp(n0+n,3,nee) = deedi*didz + deedj*djdz + deedk*dkdz
      !
    end do
    !
  end do
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP END DO
  !$OMP END PARALLEL
#endif
#endif
  !
end subroutine initialize_dusp_3d
!
!###############################################################################
!
subroutine initialize_dusp_2d()
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use geovar,          only : ncell,cell
  use metrics_mod,     only : metrics_array
  use derivatives_mod, only : deriv
  use flowvar,         only : usp,dusp
  !
  !.. Local Scalars ..
  integer  :: i,j,n,nc,n0,np
  integer  :: this_geom,this_order
  real(wp) :: a
  real(wp) :: decdi,dmxdi,dmydi,deedi,didx,didy
  real(wp) :: decdj,dmxdj,dmydj,deedj,djdx,djdy
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_dusp_2d"
  !
continue
  !
  ! Loop over the cells and get the derivative of
  ! the primitive variables at each solution point
  !
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1  ! beginning index for solpts in cell nc
    np = cell(nc)%end_sp-n0 ! number of solution points in this cell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Compute the gradient of the conservative variables at each solution point
    !
    do n = 1,np
      !
      decdi = zero
      dmxdi = zero
      dmydi = zero
      deedi = zero
      !
      do i = deriv(this_geom,this_order)%d(1)%row_ptr(n)+1, &
             deriv(this_geom,this_order)%d(1)%row_ptr(n+1)
        !
        j = deriv(this_geom,this_order)%d(1)%row_idx(i)
        a = deriv(this_geom,this_order)%d(1)%row_val(i)
        !
        decdi = decdi + a * usp(n0+j,nec)
        dmxdi = dmxdi + a * usp(n0+j,nmx)
        dmydi = dmydi + a * usp(n0+j,nmy)
        deedi = deedi + a * usp(n0+j,nee)
        !
      end do
      !
      decdj = zero
      dmxdj = zero
      dmydj = zero
      deedj = zero
      !
      do i = deriv(this_geom,this_order)%d(2)%row_ptr(n)+1, &
             deriv(this_geom,this_order)%d(2)%row_ptr(n+1)
        !
        j = deriv(this_geom,this_order)%d(2)%row_idx(i)
        a = deriv(this_geom,this_order)%d(2)%row_val(i)
        !
        decdj = decdj + a * usp(n0+j,nec)
        dmxdj = dmxdj + a * usp(n0+j,nmx)
        dmydj = dmydj + a * usp(n0+j,nmy)
        deedj = deedj + a * usp(n0+j,nee)
        !
      end do
      !
      didx = metrics_array(n0+n,1,1)
      djdx = metrics_array(n0+n,2,1)
      !
      didy = metrics_array(n0+n,1,2)
      djdy = metrics_array(n0+n,2,2)
      !
      dusp(n0+n,1,nec) = decdi*didx + decdj*djdx
      dusp(n0+n,2,nec) = decdi*didy + decdj*djdy
      !
      dusp(n0+n,1,nmx) = dmxdi*didx + dmxdj*djdx
      dusp(n0+n,2,nmx) = dmxdi*didy + dmxdj*djdy
      !
      dusp(n0+n,1,nmy) = dmydi*didx + dmydj*djdx
      dusp(n0+n,2,nmy) = dmydi*didy + dmydj*djdy
      !
      dusp(n0+n,1,nee) = deedi*didx + deedj*djdx
      dusp(n0+n,2,nee) = deedi*didy + deedj*djdy
      !
    end do
    !
  end do
  !
end subroutine initialize_dusp_2d
!
!###############################################################################
!
subroutine initialize_dusp_1d()
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use geovar,          only : ncell,cell
  use metrics_mod,     only : metrics_array
  use derivatives_mod, only : deriv
  use flowvar,         only : usp,dusp
  !
  !.. Local Scalars ..
  integer  :: i,j,n,nc,n0,np
  integer  :: this_geom,this_order
  real(wp) :: a
  real(wp) :: decdi,dmxdi,deedi,didx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_dusp_1d"
  !
continue
  !
  ! Loop over the cells and get the derivative of
  ! the primitive variables at each solution point
  !
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1  ! beginning index for solpts in cell nc
    np = cell(nc)%end_sp-n0 ! number of solution points in this cell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Compute the gradient of the conservative variables at each solution point
    !
    do n = 1,np
      !
      decdi = zero
      dmxdi = zero
      deedi = zero
      !
      do i = deriv(this_geom,this_order)%d(1)%row_ptr(n)+1, &
             deriv(this_geom,this_order)%d(1)%row_ptr(n+1)
        !
        j = deriv(this_geom,this_order)%d(1)%row_idx(i)
        a = deriv(this_geom,this_order)%d(1)%row_val(i)
        !
        decdi = decdi + a * usp(n0+j,nec)
        dmxdi = dmxdi + a * usp(n0+j,nmx)
        deedi = deedi + a * usp(n0+j,nee)
        !
      end do
      !
      didx = metrics_array(n0+n,1,1)
      !
      dusp(n0+n,1,nec) = decdi*didx
      dusp(n0+n,1,nmx) = dmxdi*didx
      dusp(n0+n,1,nee) = deedi*didx
      !
    end do
    !
  end do
  !
end subroutine initialize_dusp_1d
!
!###############################################################################
!
subroutine initialize_dusp_nd(mbeg,mend)
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use order_mod,       only : maxSP
  use geovar,          only : nr,ncell,cell
  use metrics_mod,     only : metrics_array
  use derivatives_mod, only : deriv
  use flowvar,         only : usp,dusp
  !
  !.. Formal Arguments ..
  integer, intent(in) :: mbeg,mend
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,ll,m,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,val
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: dudr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_dusp_nd"
  !
 !logical(lk), parameter :: use_conservative_form = true
  logical(lk), parameter :: use_conservative_form = fals
  !
continue
  !
  ! Loop over the cells and get the derivative of
  ! the primitive variables at each solution point
  !
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Compute the gradient of the conservative variables at each solution point
    !
    do k = 1,np
      !
      n = n0+k ! total index of solution point on current partition
      !
      do m = mbeg,mend
        !
        if (use_conservative_form) then
          !
          do l = 1,nr
            !
            val = zero
            do ll = 1,nr
              do i = deriv(this_geom,this_order)%d(ll)%row_ptr(k)+1, &
                     deriv(this_geom,this_order)%d(ll)%row_ptr(k+1)
                !
                j = deriv(this_geom,this_order)%d(ll)%row_idx(i)
                a = deriv(this_geom,this_order)%d(ll)%row_val(i)
                !
                val = val + a*usp(n0+j,m)*metrics_array(n0+j,ll,l)
                !
              end do
            end do
            !
            dusp(n,l,m) = val
            !
           !dudr(1) = deriv(this_geom,this_order)%d(1)% &
           !          dot( k , tusp(1:np,m)*metrics_dt(nc)%met(:,1,l) )
           !dudr(2) = deriv(this_geom,this_order)%d(2)% &
           !          dot( k , tusp(1:np,m)*metrics_dt(nc)%met(:,2,l) )
           !dudr(3) = deriv(this_geom,this_order)%d(3)% &
           !          dot( k , tusp(1:np,m)*metrics_dt(nc)%met(:,3,l) )
           !!
           !dusp(n,l,m) = sum(dudr)!* jacobian(n)
            !
          end do
          !
        else
          !
          ! Compute the derivatives in r-s space
          !
          do l = 1,nr
            val = zero
            do i = deriv(this_geom,this_order)%d(l)%row_ptr(k)+1, &
                   deriv(this_geom,this_order)%d(l)%row_ptr(k+1)
              !
              j = deriv(this_geom,this_order)%d(l)%row_idx(i)
              a = deriv(this_geom,this_order)%d(l)%row_val(i)
              !
              val = val + a*usp(n0+j,m)
              !
            end do
            dudr(l) = val
          end do
         !dudr(1:nr) = deriv(this_geom,this_order)%grad( k , tusp(1:np,m) )
          !
          ! Transform the derivatives back to x-y space
          !
          do l = 1,nr
            dusp(n,l,m) = dot( dudr(1:nr) , metrics_array(n,1:nr,l) )
          end do
          !
        end if
        !
      end do
      !
    end do
    !
  end do
  !
end subroutine initialize_dusp_nd
!
!###############################################################################
!
subroutine get_face_and_edge_point_coordinates(force_facexyz)
  !
  use order_mod,         only : maxSP,geom_solpts
  use geovar,            only : nr,cell,face,xyz,fp,nface
  use flowvar,           only : edgexyz_dt,facexyz_dt
  use flowvar,           only : facexyz
  use interpolation_mod, only : interp
  use ovar,              only : continuous_output
  !
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  use omp_lib, only : omp_get_thread_num
  use omp_lib, only : omp_get_num_threads
#endif
#endif
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: force_facexyz
  !
  !.. Local Scalars ..
  integer  :: n,nc,ne,nf,np,n1,n2,f0
  integer  :: i,j,k,kf,l,lf,ierr
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  integer  :: edge_order
  real(wp) :: val,a
  !
  !.. Local Arrays ..
  real(wp) :: txyz(1:maxSP,1:nr)
  !
  !.. Local Allocatable Arrays ..
  logical(lk), allocatable :: msk(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_face_and_edge_point_coordinates"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (present(force_facexyz)) then
    if (force_facexyz) then
      if (.not. allocated(facexyz_dt)) then
        allocate ( facexyz_dt(1:nface) )
        do nf = 1,nface
          np = geom_solpts(face(nf)%geom,face(nf)%order)
          facexyz_dt(nf)%v = init_array(nr,np,zero) ! F2003 AUTO-REALLOC
        end do
      end if
      np = face(nface)%end_fp
      facexyz = init_array(nr,np,zero) ! F2003 AUTO-REALLOC
    end if
  end if
  !
  ! Only do this if facexyz is allocated which means that either
  ! the continuous output option or MMS is being used
  !
  if (allocated(facexyz_dt)) then
    !
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
    !$OMP PARALLEL &
    !$OMP&         DEFAULT(SHARED) &
    !$OMP&         PRIVATE(tid,nf,face_geom,face_order,f0,nc,n1,n2,np, &
    !$OMP&                 this_geom,this_order,k,l,txyz,kf,lf,val,i,j,a)
    !
    tid = omp_get_thread_num()
    nthreads = omp_get_num_threads()
    !
    if (mypnum == glb_root) then
      if (tid == glb_root) then
        write (iout,thread_report) pname,nthreads
      end if
    end if
    !
    chunk_size = ceiling( real(size(facexyz_dt),kind=wp) / &
                          real(nthreads,kind=wp) )
    chunk_size = chunk_size / 10
    !
    !$OMP DO SCHEDULE(RUNTIME)
   !!$OMP DO SCHEDULE(GUIDED,chunk_size)
   !!$OMP DO SCHEDULE(STATIC,chunk_size)
#endif
#endif
    do nf = 1,size(facexyz_dt)
      !
      face_geom = face(nf)%geom
      face_order = face(nf)%order
      f0 = face(nf)%beg_fp-1
      !
      nc = face(nf)%left%cell
      !
      n1 = cell(nc)%beg_sp
      n2 = cell(nc)%end_sp
      np = n2-n1+1
      !
      this_geom = cell(nc)%geom
      this_order = cell(nc)%order
      !
      ! Create transpose of xyz array to reduce cache misses
      !
      do k = 1,np
        do l = 1,nr
          txyz(k,l) = xyz(l,k+n1-1)
        end do
      end do
      !
      ! Interpolate the solution point coordinates to each flux point
      ! on the current face
      !
      do kf = 1,geom_solpts(face_geom,face_order)
        !
        lf = face(nf)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                          rot(face(nf)%left%rotation)% &
                                       pts(kf)
        !
        do l = 1,nr
          val = zero
          do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1,&
                 interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
            !
            j = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
            a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
            !
            val = val + a * txyz(j,l)
            !
          end do
          facexyz_dt(nf)%v(l,kf) = val
          facexyz(l,f0+kf) = val
        end do
        !
      end do
      !
    end do
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
    !$OMP END DO
    !$OMP END PARALLEL
#endif
#endif
    !
    !$ACC UPDATE DEVICE(facexyz)
    !
    if (.not. continuous_output) then
      deallocate ( facexyz_dt , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"facexyz_dt",2,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
    !
  end if
  !
  ! Only do this if edgexyz is allocated which means
  ! that the continuous output option is being used
  !
  if (allocated(edgexyz_dt)) then
    !
    ! Initialize edgexyz to zero
    !
#ifdef SPECIAL_FOR_PGI
    call edgexyz_dt%zero_ne2d_t
#else
    call edgexyz_dt%zero
#endif
    !
    ! Allocate the mask array to mark which edges have been finished
    !
    allocate ( msk(1:size(edgexyz_dt)) , source=fals , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"msk",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Loop through the cells and compute the coordinates
    ! of the edge points on each grid edge by interpolating
    ! the coordinates of the cell solution points
    !
    do nc = 1,size(cell)
      !
      n1 = cell(nc)%beg_sp
      n2 = cell(nc)%end_sp
      np = n2-n1+1
      !
      this_geom = cell(nc)%geom
      this_order = cell(nc)%order
      !
      ! Create transpose of xyz array to reduce cache misses
      !
      do k = 1,np
        do l = 1,nr
          txyz(k,l) = xyz(l,k+n1-1)
        end do
      end do
      !
      ! Interpolate the solution point coordinates to each edge point
      ! on the current edge
      !
      edges_of_cell: do n = 1,size(cell(nc)%edge)
        !
        ne = cell(nc)%edge(n)%idx
        !
        ! Skip to the next edge of this cell if this edge has already been done
        !
        if (msk(ne)) cycle edges_of_cell
        !
        edge_order = cell(nc)%edge(n)%order
        !
        ! Loop over the edge points on this edge and interpolate the
        ! coordinates of the cell solution points to each edge point
        !
        do kf = 1,edge_order+1
          !
          lf = cell(nc)%edge(n)%ep_offset + merge( edge_order+2-kf , kf , &
                                                  cell(nc)%edge(n)%rev_dir)
          !
          do l = 1,nr
            val = zero
            do i = interp(this_geom,this_order)% &
                             toEdge(edge_order)%row_ptr(lf)+1, &
                   interp(this_geom,this_order)% &
                             toEdge(edge_order)%row_ptr(lf+1)
              !
              j = interp(this_geom,this_order)%toEdge(edge_order)%row_idx(i)
              a = interp(this_geom,this_order)%toEdge(edge_order)%row_val(i)
              !
              val = val + a * txyz(j,l)
              !
            end do
            edgexyz_dt(ne)%v(l,kf) = val
          end do
          !
        end do
        !
        ! Mark this edge as finished
        !
        msk(ne) = true
        !
      end do edges_of_cell
      !
    end do
    !
    ! Deallocate the msk array before leaving
    !
    deallocate ( msk , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_face_and_edge_point_coordinates
!
!###############################################################################
!
subroutine initialize_reference
  !
  !.. Use Statements ..
  use ovar, only : gam,tmin,tmax
  use ovar, only : rhoref,aref,tref,pref
  use ovar, only : machref,reyref,time_ref
  use ovar, only : rgasref_nd,cpref_nd
  use ovar, only : cv_ref,pv_ref,gp_ref
  use ovar, only : suth_Sref,suth_Tref
  use ovar, only : alpha_aoaref,beta_aoaref
  use ovar, only : itestcase,Lref
  use ovar, only : ideal_total_temperature
  use ovar, only : limit_total_temperature
  use ovar, only : moving_ref_frame
  !
  use thermo_mod,    only : init_thermo_mod
  use transport_mod, only : init_transport_mod
  use reactions_mod, only : init_reactions_mod
  !
  use taylor_couette_mod, only : initialize_taylor_couette_problem
  !
  !.. Local Scalars ..
  integer  :: ierr
  real(wp) :: rhoa,rhoa2,tda2,q2
  real(wp) :: gam_over_gm1
  real(wp) :: gm1_over_gam
  real(wp) :: pres_check
  real(wp) :: ref_length
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_reference"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  tmin = 10.0_wp
  tmax = 15000.0_wp
  !
  ref_length = one
  !
  ! Allocate the free stream inflow arrays
  !
  if (allocated(cv_ref)) then
    deallocate ( cv_ref , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cv_ref",2,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
  end if
  if (allocated(pv_ref)) then
    deallocate ( pv_ref , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"pv_ref",2,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
  end if
  if (allocated(gp_ref)) then
    deallocate ( gp_ref , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"gp_ref",2,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
  end if
  !
  allocate ( cv_ref(1:nq) , source=zero , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cv_ref",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  allocate ( pv_ref(1:nq) , source=zero , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"pv_ref",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  allocate ( gp_ref(1:ng) , source=zero , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"gp_ref",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  !
  ! Set up the base reference conditions
  !
  if (any(itestcase == Advec_Problems) .or. &
      abs(itestcase) == Shu_Vortex .or. &
      abs(itestcase) == Shu_MultiVortex .or. &
      itestcase == -Vortex_Transport) then
    !
    ! Set up the base reference conditions for different
    ! versions of the isentropic Euler vortex problem
    ! NOTE: This also sets the values of the following variables in ovar:
    !             rgasref, gam, rhoref, pref, tref, cpref, aref,
    !             uref, machref, velref, alpha_aoaref, beta_aoaref,
    !             velref_nd, pref_nd, rgasref_nd, and cpref_nd
    !
    call initialize_vortex_problem(cv_ref,pv_ref,gp_ref)
    !
  else
    !
    if (itestcase == Taylor_Green_Vortex) then
      !
      ! If Lref is negative, the grid was not created on the fly and
      ! we need to add code to compute Lref and make sure the input
      ! grid is valid for the Taylor-Green vortex problem
      !
      if (Lref <= zero) then
        write (error_message,101)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      !
      ! Set up the base reference conditions for the Taylor-Green Vortex problem
      !
      alpha_aoaref = zero
      beta_aoaref = zero
      !
      tref = 300.0_wp ! (K)
      pref = 101325.0_wp ! (Pa)
      rhoref = -one
      !
      machref = 0.1_wp
      !
      if (reyref <= zero) then
        reyref = 1600.0_wp
      end if
      !
      suth_Sref = 110.5556_wp
      suth_Tref = 273.111_wp
      !
    else if (abs(itestcase) == Taylor_Couette_Flow) then
      !
      ! Set the reference conditions for the Taylor-Couette flow problem
      ! NOTE: After calling this subroutine, subroutines init_thermo_mod and
      !       init_transport_mod will be called to set the thermo and
      !       transport refrerence conditions. Some of these will have already
      !       been set in this initialize_taylor_couette_problem, but these
      !       values should stay the same even if they are recomputed. These
      !       two subroutines (thermo/transport) still need to be called
      !       because they set a lot of other reference variables that are
      !       not evaluated in this subroutine.
      !
      call initialize_taylor_couette_problem(ref_length)
      !
    end if
    !
    ! Initialize the thermodynamic gas properties module and
    ! compute the reference conditions for the gas properties.
    ! NOTE: This also sets the values of the following variables in ovar:
    !          rgasref, gam, rhoref, pref, tref, cpref, aref, uref,
    !          velref, velref_nd, pref_nd, rgasref_nd, and cpref_nd
    !
    call init_thermo_mod(cv_ref,pv_ref,gp_ref)
    !
    call check_for_bad_number(__LINE__,__FILE__)
    !
    ! Initialize the reactions module
    !
    call init_reactions_mod
    !
    call check_for_bad_number(__LINE__,__FILE__)
    !
  end if
  !
  ! Set the reference time
  !
  if (any(itestcase == [Taylor_Problems,Infinite_Cylinder])) then
    ! The time used within the solver is nondimensionalized by aref, so
    ! set the reference time so that it re-nondimensionalizes the
    ! time using uref instead of aref.
    ! NOTE: For the Taylor-Couette flow problem, time_ref needs to include
    !       the difference in radii between the two cylinders which was
    !       returned as ref_length from initialize_taylor_couette_problem.
    !       For all other problems, ref_length should just be 1.0.
    time_ref = machref/ref_length
  else
    ! Set the reference time to one/aref which converts the
    ! nondimensional time used in the solver to the physical time.
    ! NOTE: For the various (vortex/advection) transport problems,
    !       aref=1, so time_ref should just be 1.
    time_ref = one / aref
  end if
  !
  ! Nondimensionalize the ideal total temperature
  !
  if (ideal_total_temperature > zero) then
    ideal_total_temperature = ideal_total_temperature / tref
  end if
  !
  ! Initialize the transport module and compute
  ! the transport reference conditions
  ! NOTE: For the Taylor-Couette flow problem, the length scale used for the
  !       reference Reynolds number, reyref, is the difference in radii between
  !       the two cylinders. This distance was returned as ref_length from
  !       initialize_taylor_couette_problem. Since reyref has now been set, we
  !       will pass ref_length to init_transport_mod to scale muref by this
  !       reference length. For all other problems, ref_length should just be 1.
  !
  call init_transport_mod(cv_ref,gp_ref,ref_length)
  !
  call check_for_bad_number(__LINE__,__FILE__)
  !
  ! Initialize the turbulence reference conditions
  !
  call initialize_turbulence(cv_ref,pv_ref,gp_ref)
  !
  call check_for_bad_number(__LINE__,__FILE__)
  !
  ! Nondimensionalize the moving reference frame rotation rate
  ! and translational velocity
  !
  call moving_ref_frame%make_nondim()
  !
  !$ACC UPDATE DEVICE(aref,rhoref,tref,tmin,tmax, &
  !$ACC&              gam,rgasref_nd,cpref_nd, &
  !$ACC&              ideal_total_temperature, &
  !$ACC&              limit_total_temperature)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  101 format ( &
     " The grid being used for the Taylor-Green vortex problem was not", &
     " generated on the fly and seems to be input from a grid file. Code", &
     " still needs to be added to compute Lref and make sure the input", &
     " grid is valid for the Taylor-Green vortex problem.")
  !
end subroutine initialize_reference
!
!###############################################################################
!
subroutine initialize_vortex_problem(cv,pv,gp)
  !
  !.. Use Statements ..
  use ovar, only : rhoref,aref,tref,gam
  use ovar, only : velref,velref_nd,uref
  use ovar, only : pref,pref_nd,machref
  use ovar, only : rgasref,rgasref_nd
  use ovar, only : cpref,cpref_nd
  use ovar, only : alpha_aoaref,beta_aoaref
  use ovar, only : itestcase
  use ovar, only : tmin,tmax
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: cv(:)
  real(wp), intent(inout) :: pv(:)
  real(wp), intent(inout) :: gp(:)
  !
  !.. Local Scalars ..
  real(wp) :: q2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_vortex_problem"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Set up the base reference conditions for different
  ! versions of the isentropic Euler vortex problem
  !
  ! Set some reference conditions common to all of these vortex problems
  !
  beta_aoaref = zero
  velref(:) = zero
  !
  if (any(itestcase == [-Shu_Vortex,-Shu_MultiVortex])) then
    !
    ! Stationary vortex problem
    !
    alpha_aoaref = 90.0_wp
    itestcase = abs(itestcase)
    !
  else if (itestcase == -Vortex_Transport) then
    !
    ! Alteration to make the Shu case propagate in the
    ! y-direction only instead of in both x and y directions
    !
    alpha_aoaref = 90.0_wp
    velref(2) = one
    itestcase = Shu_Vortex
    !
  else
    !
    ! Original Shu vortex problem with diagonal propagation
    !
    alpha_aoaref = 45.0_wp
    velref(1) = one
    velref(2) = one
    !
  end if
  !
  gam = 1.4_wp
  machref = one/sqrt(gam)
  aref = one
  uref = one
  rhoref = one
  pref = one
  tref = one
  rgasref = one
  !
  pv(nec)     = rhoref / rhoref
  pv(nmb:nme) = velref / aref
  pv(nee)     = pref   / (rhoref*aref*aref)
  !
  q2 = dot_product( pv(nmb:nme) , pv(nmb:nme) )
  !
  cv(nec)     = pv(nec)
  cv(nmb:nme) = pv(nec)*pv(nmb:nme)
  cv(nee)     = pv(nee)/(gam-one) + half*pv(nec)*q2
  !
  gp(npres) = pv(nee)
  gp(ntemp) = tref / tref
  gp(ngamm) = gam
  !
  ! Nondimensional reference velocities
  velref_nd = velref / aref
  !
  ! Nondimensional reference pressure
  pref_nd = pref / (rhoref*aref*aref)
  !
  ! Nondimensional reference gas constant
  rgasref_nd = rgasref * tref / (aref*aref)
  !
  ! Dimensional reference specific heat at constant pressure
  cpref = rgasref * gam / (gam - one)
  ! Nondimensional reference specific heat at constant pressure
  cpref_nd = cpref * tref / (aref*aref)
  !
  tmin = 0.025_wp ! roughly default tmin/tref
  tmax = 50.0_wp  ! roughly default tmax/tref
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine initialize_vortex_problem
!
!###############################################################################
!
subroutine initialize_turbulence(cv,pv,gp)
  !
  !.. Use Statements ..
  use ovar, only : tkeref, tkeref_nd
  use ovar, only : tlsref, tlsref_nd
  use ovar, only : tvrref, tinref
  use ovar, only : rhoref, muref
  use ovar, only : aref, machref
  use ovar, only : turbulence_model
  use ovar, only : turb_is_on
  use ovar, only : include_tke_in_energy
  !
  use thermo_mod, only : update_thermo_properties
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: cv(:)
  real(wp), intent(inout) :: pv(:)
  real(wp), intent(inout) :: gp(:)
  !
  !.. Local Scalars ..
  integer :: m
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_turbulence"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Turbulence equations
  !
  if (turbulence_model /= laminar_flow) then
    !
    ! Set the flag indicating that turbulence is active
    !
    turb_is_on = true
    !
    ! Set flag that says whether TKE is included in the total energy
    !
    if (abs(turbulence_model) == k_omega) then
      include_tke_in_energy = true
    end if
    !
    ! Dimensional reference conditions
    !
    tkeref = 1.5_wp * (tinref*machref*aref)**2
    tlsref = rhoref * tkeref / (tvrref*muref)
    !
    ! Nondimensional reference conditions
    !
    tkeref_nd = tkeref / (rhoref * aref * aref)
    tlsref_nd = tlsref / (rhoref * aref)
    !
    ! Set the turbulence conserved variables inflow conditions
    !
    cv(ntk) = cv(nec) * tkeref_nd
    cv(ntl) = cv(nec) * tlsref_nd
    !
    ! Add the tke contribution to the energy equation
    !
    cv(nee) = cv(nee) + cv(ntk)
    !
    gp(:) = zero ! initialize all of gp to zero
    !
    ! Create the gas properties array from the conserved variables
    ! NOTE: Initialize the temperature location of the gas properties array
    !       to nondimensional tref (i.e., tref_nd = tref/tref = 1). This is
    !       needed if using a thermally perfect gas mixture, because we need
    !       an initial guess for temperature when using Newtons method to
    !       solve for temperature.
    !
    gp(ntemp) = one ! tref/tref
    !
    call update_thermo_properties(cv,gp,caller_name=pname)
    !
    ! Recompute the primitive variable inflow conditions
    !
    if (ncs > 0) then
      do m = 1,ncs
        pv(m) = cv(m)/cv(nec)  ! Species continuity equations
      end do
    end if
    !
    pv(nec) = cv(nec)          ! Copy density from cv
    !
    do m = nmb,nme
      pv(m) = cv(m)/cv(nec)    ! Momentum equations
    end do
    !
    pv(nee) = gp(npres)        ! Copy pressure from gp array
    !
    do m = nee+1,nq
      pv(m) = cv(m)/cv(nec)    ! Turbulence equations
    end do
    !
  end if
  !
  !$ACC UPDATE DEVICE(turb_is_on,include_tke_in_energy)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine initialize_turbulence
!
!###############################################################################
!
subroutine create_bc_in
  !
  ! This subroutine takes the bc_input derived type (which is used in the
  ! namelist input file to specify the flow parameters for each boundary
  ! condition) and creates the bc_in derived type (which, for each boundary
  ! condition, contains the actual boundary states and other boundary
  ! information [SEM, MRF, periodic rotations, etc.] that are needed/used
  ! in the solver.
  !
  ! This should be called before partitioning, but after 1) initializing the
  ! reference conditions (rhoref, aref, tref, etc.), 2) matching all periodic
  ! boundary faces, and 3) separating out the boundary faces on an SEM
  ! boundary condition that are outside the SEM profile height and moving
  ! those faces to a separate traditional (i.e., steady) inflow boundary
  ! condition.
  !
  !.. Use Statements ..
  use geovar, only : nr,ncell,grid,bface
  use ovar,   only : bc_input,bc_in,bc_names
  use ovar,   only : tmin,tmax,tref,rgasref
  use ovar,   only : rhoref,aref,gam
  use ovar,   only : cv_ref,pv_ref,gp_ref
  use ovar,   only : governing_equations
  use ovar,   only : iso_wall_temp
  !
  use periodic_mod, only : get_pbc_matrix
  !
  use thermo_mod, only : compute_mixture_rgas
  use thermo_mod, only : compute_mixture_gamma
  !
  !.. Local Scalars ..
  integer     :: i,ii,m,n,nc,nf,nbc,ierr,bc_idx
  real(wp)    :: r,p,t,tnd,r0,p0,t0,mach,aspd,visc,reyn,rgas,gamm
  real(wp)    :: alpha,beta,total_cond
  real(wp)    :: gm1,gam_over_gm1,one_over_gm1
  logical(lk) :: static_is_set
  logical(lk) :: total_is_set
  logical(lk) :: mach_is_set
  logical(lk) :: vel_is_set
  logical(lk) :: face_was_reset_to_slip_wall
  logical(lk) :: this_is_a_wall_bc
  !
  character(len=32) :: bcnam,bcinp
  !
  !.. Local Arrays ..
  real(wp) :: vel(1:nr)
  real(wp) :: cv(1:nq)
  real(wp) :: pv(1:nq)
  real(wp) :: gp(1:ng)
  real(wp) :: mass_frac(1:ncs)
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: bc_list(:)
  integer, allocatable :: rev_bc_list(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_bc_in"
  !
#ifdef DEBUG_ON
!#define USE_DUMP_BC_CONDITIONS
#endif
#ifdef USE_DUMP_BC_CONDITIONS
  logical(lk), parameter :: debug_dump_bc = true
#else
  logical(lk), parameter :: debug_dump_bc = fals
#endif
  !
continue
  !
  if (allocated(bc_in)) then
    return
  end if
  !
  call debug_timer(entering_procedure,pname)
  !
#ifdef DEBUG_ON
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  if (mypnum == glb_root) then
    write (iout,260)
    do n = lbound(bc_input,dim=1),ubound(bc_input,dim=1)
      write (iout,261) "bc_input",n,trim(adjustl(bc_input(n)%name)), &
                       bc_integer_to_string(bc_input(n)%set_bc_value)
    end do
    write (iout,260)
    flush (iout)
    write (iout,260)
    do n = lbound(bc_names,dim=1),ubound(bc_names,dim=1)
      write (iout,261) "bc_names",n,trim(adjustl(bc_names(n)))
    end do
    write (iout,260)
    flush (iout)
  end if
  260 format (/,2(2x,60("#"),/))
  261 format (2x,a,"(",i0.2,") = '",a,"'",:," ;  BC_Type = ",a)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  flush (iout)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
#endif
  !
  ! Change all no-slip wall boundary conditions to
  ! slip walls if we are solving the Euler equations.
  !
  if (governing_equations == Euler_Eqns) then
    !
    face_was_reset_to_slip_wall = fals
    !
    if (allocated(grid)) then
      if (allocated(grid%elem)) then
        do nc = ncell+1,ubound(grid%elem,dim=1)
          if (allocated(grid%elem(nc)%tags)) then
            if (any(grid%elem(nc)%tags(1) == bc_noslip_walls)) then
              grid%elem(nc)%tags(1) = bc_slip_wall
              face_was_reset_to_slip_wall = true
            end if
          end if
        end do
      end if
    end if
    !
    do nf = 1,size(bface,dim=2)
      if (any(bface(1,nf) == bc_noslip_walls)) then
        bface(1,nf) = bc_slip_wall
        face_was_reset_to_slip_wall = true
      end if
    end do
    !
    do n = 1,size(bc_input)
      if (any(bc_input(n)%set_bc_value == bc_noslip_walls)) then
        bc_input(n)%set_bc_value = bc_slip_wall
        face_was_reset_to_slip_wall = true
      end if
    end do
    !
    if (face_was_reset_to_slip_wall) then
      if (mypnum == glb_root) then
        write (iout,3)
        flush (iout)
      end if
    end if
    !
  end if
  !
  ! Find any boundary conditions that were set in the
  ! input namelist file using the bc_input array
  !
  do n = 1,size(bc_input)
    !
    bc_input(n)%name = adjustl(bc_input(n)%name)
    !
   !if (len_trim(bc_input(n)%name) > 0) then
    if ( bc_input(n)%is_not_default() ) then
      if (.not. allocated(bc_list)) then
        bc_list = [ n ] ! Using F2003 auto-reallocation
       !allocate ( bc_list(1:1) , source=n , stat=ierr , errmsg=error_message )
       !call alloc_error(pname,"bc_list",1,__LINE__,__FILE__,ierr,error_message)
      else
        bc_list = [ bc_list , n ] ! Using F2003 auto-reallocation
      end if
    end if
    !
  end do
  !
  nbc = 0
  if (allocated(bc_list)) nbc = size(bc_list)
  !
  allocate ( bc_in(0:nbc) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bc_in",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Set the reference boundary conditions
  !
  bc_in(0)%cv = cv_ref(:) ! Using F2003 auto-reallocation
  bc_in(0)%pv = pv_ref(:) ! Using F2003 auto-reallocation
  bc_in(0)%gp = gp_ref(:) ! Using F2003 auto-reallocation
  !
  ! Set the other boundary conditions
  !
  do i = 1,nbc
    !
    n = bc_list(i)
    !
    ! Save the index within bc_input for this bc_in
    !
    bc_in(i)%bc_input_idx = n
    !
    ! Set a logical that identifies if the
    ! current boundary condition is a wall
    !
    this_is_a_wall_bc = ( any(bc_input(n)%set_bc_value == bc_walls) )
    !
    if ( bc_input(n)%has_default_flow_conditions() .or. &
         bc_input(n)%set_bc_value == bc_periodic ) then
      !
      ! This boundary condition has either:
      !   1) all the default flow conditions but it contains a 
      !      non-zero-length name or specified BC type, or
      !   2) it is a periodic BC.
      !      NOTE: We include this exception in case any of the flow conditions
      !            in bc_input have been altered from the default values for a
      !            periodic BC, in which case the BC flow conditions dont
      !            matter since the left and right states will just be
      !            communicated between matching periodic faces, so just set
      !            the BC flow conditions to the reference flow conditions.
      !
      ! Set the flow conditions for this boundary to the reference
      ! flow conditions.
      !
#ifdef DEBUG_ON
      if (mypnum == glb_root) then
        write (iout,1) i,trim(adjustl(bc_input(n)%name)), &
          trim(adjustl(bc_integer_to_string(bc_input(n)%set_bc_value)))
      end if
#endif
      !
      if (this_is_a_wall_bc) then
        !
        ! Set the boundary conditions for this wall BC to the reference
        ! density, pressure, temperature, and zero velocity
        !
        pv(1:nec) = bc_in(0)%pv(1:nec)
        pv(nmb:nme) = zero
        pv(nee) = bc_in(0)%pv(nee)
        !
        gp(npres) = bc_in(0)%gp(npres)
        !
        if (bc_input(n)%set_bc_value == bc_isothermal_wall) then
          bc_in(i)%wall_temperature = iso_wall_temp/tref
          gp(ntemp) = bc_in(i)%wall_temperature
        else
          gp(ntemp) = bc_in(0)%gp(ntemp)
        end if
        !
        call compute_flow_conditions_from_pvt(cv,pv,gp)
        !
      else
        !
        ! Set the boundary conditions for this BC to the reference conditions
        !
        cv(:) = bc_in(0)%cv(:)
        pv(:) = bc_in(0)%pv(:)
        gp(:) = bc_in(0)%gp(:)
        !
      end if
      !
      bc_in(i)%cv = cv(:) ! Using F2003 auto-reallocation
      bc_in(i)%pv = pv(:) ! Using F2003 auto-reallocation
      bc_in(i)%gp = gp(:) ! Using F2003 auto-reallocation
      !
      ! If this is a periodic BC, we need to create the rotation matrix
      ! if this is a rotated periodic BC.
      ! NOTE: Only allocate pbc_rotation%mat if this is a rotated periodic BC
      !       and do NOT allocate if this is a translation-only periodic BC.
      !
      if ((any(abs(bc_input(n)%periodic%RotationAxis(1:3)) > eps6) .and. &
           abs(bc_input(n)%periodic%RotationAxis(4)) > eps6) .or. &
          (any(abs(bc_input(n)%periodic%RotationAngle(1:3)) > eps6))) then
        ! F2003 AUTO-REALLOC
        bc_in(i)%pbc_rotation%mat = get_pbc_matrix(bc_input(n)%periodic, &
                                                   ignore_translation=true)
      end if
      !
    else
      !
      static_is_set = fals
      total_is_set = fals
      mach_is_set = fals
      vel_is_set = fals
      !
      ! Copy over the relaxation information
      !
      bc_in(i)%relax%value = min(one,max(bc_input(n)%relax%value,zero))
      bc_in(i)%relax%time = max(bc_input(n)%relax%time,zero)
      bc_in(i)%relax%iter = max(bc_input(n)%relax%iter,0)
      !
      ! Copy over the SEM information
      !
      bc_in(i)%sem = bc_input(n)%sem
      !
      ! Make sure certain SEM information is realizable
      !
      bc_in(i)%sem%length_scale = max(zero,bc_in(i)%sem%length_scale)
      bc_in(i)%sem%profile_height = max(zero,bc_in(i)%sem%profile_height)
      !
      ! Trim the profile file character string
      !
      bc_in(i)%sem%profile_file = trim(adjustl(bc_in(i)%sem%profile_file))
      !
      ! Nondimensionalize the SEM convective velocity
      !
      bc_in(i)%sem%convective_velocity = bc_in(i)%sem%convective_velocity/aref
      !
      ! Set the flow conditions for this BC depending on whether
      ! or not a wall temperature has been specified for this BC.
      !
      if (this_is_a_wall_bc) then
        !
        if (bc_input(n)%set_bc_value == bc_isothermal_wall) then
          !
          ! Set the wall temperature for an isothermal wall
          !
          if ( (bc_input(n)%wall_temp >= tmin) .and. &
               (bc_input(n)%wall_temp <= tmax) ) then
            !
            ! A wall temperature has been specified for this BC
            ! NOTE: The only values likely to be used for an isothermal wall
            !       boundary condition are:
            !              bc_in(i)%wall_temperature
            !              bc_in(i)%pv(nmb:nme)
            !       All other quantities will come from the interior.
            !
            bc_in(i)%wall_temperature = bc_input(n)%wall_temp/tref
            !
          else
            !
            bc_in(i)%wall_temperature = iso_wall_temp/tref
            !
          end if
          !
          t = bc_in(i)%wall_temperature
          !
        else
          !
          ! Set temperature to the reference temperature for an adiabatic wall
          !
          t = bc_in(0)%gp(ntemp)
          !
        end if
        !
        ! Get the wall rotation speed from bc_input and save it
        ! to bc_in in the form that is needed for the solver.
        !
        bc_in(i)%rotation_rate(:) = get_bc_rotation( bc_input(n) )
        !
        ! Get the wall translational velocity from bc_input
        ! and save it to the temporary array vel(:)
        !
        vel(:) = get_bc_velocity( bc_input(n) , dimensional=fals )
        !
        ! Compute the flow conditions for this BC using the reference static
        ! density, the reference static pressure, the specified BC velocity,
        ! and the specified BC wall temperature.
        !
        r = bc_in(0)%pv(nec)
        p = bc_in(0)%gp(npres)
        !
        if (ncs > 0) then
          do m = 1,ncs
            pv(m) = bc_in(0)%pv(m)
          end do
        end if
        !
        pv(nec) = r
        pv(nmb:nme) = vel(1:nr)
        pv(nee) = p
        !
        gp(npres) = p
        gp(ntemp) = t
        !
        call compute_flow_conditions_from_pvt(cv,pv,gp)
        !
        bc_in(i)%cv = cv(:) ! USING F2003 AUTO-REALLOCATION
        bc_in(i)%pv = pv(:) ! USING F2003 AUTO-REALLOCATION
        bc_in(i)%gp = gp(:) ! USING F2003 AUTO-REALLOCATION
        !
      else
        !
        if (ncs > 0) mass_frac(1:ncs) = bc_input(n)%mass_frac(1:ncs)
        !
        r = bc_input(n)%rho_static
        p = bc_input(n)%p_static
        t = bc_input(n)%t_static
        !
        r0 = bc_input(n)%rho_total
        p0 = bc_input(n)%p_total
        t0 = bc_input(n)%t_total
        !
        mach  = bc_input(n)%mach
        alpha = bc_input(n)%alpha_aoa
        beta  = bc_input(n)%beta_aoa
        !
        if (ncs > 0) then
          call compute_mixture_rgas(mass_frac,rgas,dimensional=true)
          gamm = -one
        else
          gamm = gam
          rgas = rgasref
        end if
        !
        ! Get the BC velocity from bc_input
        !
        vel(:) = get_bc_velocity( bc_input(n) , dimensional=true )
        !
        iterate_bc_conditions: do ii = 1,5
          !
          ! The Mach number has been set if:  0 < mach < 25
          mach_is_set = ( abs(mach) < five*five .and. mach > zero )
          ! The velocity has been set if:  |vel| > 0
          vel_is_set  = ( norm2(vel) > zero )
          !
          call igl(rgas,r,p,t,static_is_set,n)
          call dump_bc_conditions(ii,gamm,rgas,r,p,t,r0,p0,t0, &
                                  vel,mach,alpha,beta,mass_frac, &
                                  static_is_set,total_is_set, &
                                  vel_is_set,mach_is_set, &
                                  "after igl(r,p,t)", &
                                  use_dump_bc=debug_dump_bc)
          !
          if (gamm <= zero) then
            if (t > tmin) then
              tnd = t/tref
              call compute_mixture_gamma(mass_frac,tnd,gamm,dimensional=true)
              call dump_bc_conditions(ii,gamm,rgas,r,p,t,r0,p0,t0, &
                                      vel,mach,alpha,beta,mass_frac, &
                                      static_is_set,total_is_set, &
                                      vel_is_set,mach_is_set, &
                                      "after compute_mixture_gamma(1)", &
                                      use_dump_bc=debug_dump_bc)
            end if
          end if
          !
          call igl(rgas,r0,p0,t0,total_is_set,n)
          call dump_bc_conditions(ii,gamm,rgas,r,p,t,r0,p0,t0, &
                                  vel,mach,alpha,beta,mass_frac, &
                                  static_is_set,total_is_set, &
                                  vel_is_set,mach_is_set, &
                                  "after igl(r0,p0,t0)", &
                                  use_dump_bc=debug_dump_bc)
          !
          if (gamm > zero) then
            !
            call mach_total_pres(gamm,p,p0,mach,mach_is_set,n)
            call dump_bc_conditions(ii,gamm,rgas,r,p,t,r0,p0,t0, &
                                    vel,mach,alpha,beta,mass_frac, &
                                    static_is_set,total_is_set, &
                                    vel_is_set,mach_is_set, &
                                    "after mach_total_pres", &
                                    use_dump_bc=debug_dump_bc)
            !
            call mach_total_temp(gamm,t,t0,mach,mach_is_set,n)
            call dump_bc_conditions(ii,gamm,rgas,r,p,t,r0,p0,t0, &
                                    vel,mach,alpha,beta,mass_frac, &
                                    static_is_set,total_is_set, &
                                    vel_is_set,mach_is_set, &
                                    "after mach_total_temp", &
                                    use_dump_bc=debug_dump_bc)
            !
            call mach_total_dens(gamm,r,r0,mach,mach_is_set,n)
            call dump_bc_conditions(ii,gamm,rgas,r,p,t,r0,p0,t0, &
                                    vel,mach,alpha,beta,mass_frac, &
                                    static_is_set,total_is_set, &
                                    vel_is_set,mach_is_set, &
                                    "after mach_total_dens", &
                                    use_dump_bc=debug_dump_bc)
            !
          else if (t0 > zero .and. mach_is_set) then
            !
            call mach_total_temp_tpg(mass_frac,rgas,gamm,t,t0,mach)
            call dump_bc_conditions(ii,gamm,rgas,r,p,t,r0,p0,t0, &
                                    vel,mach,alpha,beta,mass_frac, &
                                    static_is_set,total_is_set, &
                                    vel_is_set,mach_is_set, &
                                    "after mach_total_temp_tpg", &
                                    use_dump_bc=debug_dump_bc)
            !
          end if
          !
          if (gamm <= zero) then
            if (t > tmin) then
              tnd = t/tref
              call compute_mixture_gamma(mass_frac,tnd,gamm,dimensional=true)
              call dump_bc_conditions(ii,gamm,rgas,r,p,t,r0,p0,t0, &
                                      vel,mach,alpha,beta,mass_frac, &
                                      static_is_set,total_is_set, &
                                      vel_is_set,mach_is_set, &
                                      "after compute_mixture_gamma(2)", &
                                      use_dump_bc=debug_dump_bc)
            end if
          end if
          !
          call vel_from_static(gamm,rgas,r,p,t,vel,mach,alpha,beta, &
                               vel_is_set,mach_is_set)
          call dump_bc_conditions(ii,gamm,rgas,r,p,t,r0,p0,t0, &
                                  vel,mach,alpha,beta,mass_frac, &
                                  static_is_set,total_is_set, &
                                  vel_is_set,mach_is_set, &
                                  "after vel_from_static", &
                                  use_dump_bc=debug_dump_bc)
          !
          if (gamm <= zero) then
            if (t > tmin) then
              tnd = t/tref
              call compute_mixture_gamma(mass_frac,tnd,gamm,dimensional=true)
            end if
          end if
          !
          if (static_is_set .and. vel_is_set) exit iterate_bc_conditions
          !
          if (ii == 5) then
            write (error_message,2) n
            call dump_bc_conditions(ii,gamm,rgas,r,p,t,r0,p0,t0, &
                                    vel,mach,alpha,beta,mass_frac, &
                                    static_is_set,total_is_set, &
                                    vel_is_set,mach_is_set, &
                                    "Reached max BC iterations!", &
                                    use_dump_bc=true)
            call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
          end if
          !
        end do iterate_bc_conditions
        !
        if (ncs > 0) then
          do m = 1,ncs
            pv(m) = mass_frac(m)
          end do
        end if
        !
        pv(nec) = r
        pv(nmb:nme) = vel(1:nr)
        pv(nee) = p
        !
        gp(npres) = p
        gp(ntemp) = t
        !
        call compute_flow_conditions_from_pvt(cv,pv,gp,dimensional=true)
        !
        call check_bc_input(bc_input(n),pv,gp,n)
        !
        pv(nec) = pv(nec) / rhoref
        pv(nmb:nme) = pv(nmb:nme) / aref
        pv(nee) = pv(nee) / (rhoref*aref*aref)
        !
        gp(npres) = gp(npres) / (rhoref*aref*aref)
        gp(ntemp) = gp(ntemp) / tref
        !
        call compute_flow_conditions_from_pvt(cv,pv,gp)
        !
        bc_in(i)%cv = cv(:) ! Using F2003 auto-reallocation
        bc_in(i)%pv = pv(:) ! Using F2003 auto-reallocation
        bc_in(i)%gp = gp(:) ! Using F2003 auto-reallocation
        !
      end if
      !
    end if
    !
    ii = 0
    bface_loop: do nf = 1,size(bface,dim=2)
      !
     !if (any(bface(1,nf) == bc_comm) .or. bface(5,nf) == 0) cycle bface_loop
      if (bface(1,nf) == bc_cpu_bnd .or. bface(5,nf) == 0) cycle bface_loop
      !
     !if (mypnum == glb_root) then
     !  write (630,173) n,nf,bface(5,nf),mod(bface(11,nf),10000), &
     !                  trim(adjustl(bc_integer_to_string(bface(1,nf))))
     !  173 format (i3,i8,i8,i8,4x,a)
     !end if
      !
      if (bface(1,nf) == bc_periodic) then
        !
        ! If this is a periodic BC, we need to extract bc_idx from bface(11)
        ! to see if this is the same as the current bc_input index, meaning
        ! this face is part of the current bc_input and add it to the count
        ! of total boundary faces for the current bc_input. If it is:
        !   1) add this face to the count of total boundary faces for the
        !      current bc_input
        !   2) save the BC type to bc_in if it still bc_unknown
        !
        bc_idx = mod(bface(11,nf),10000)
        !
        if (bc_idx == n) then
          !
          ii = ii + 1
          !
          ! Save the boundary condition type to bc_in
          !
          if (bc_in(i)%bc_type == bc_unknown) then
            bc_in(i)%bc_type = bface(1,nf)
          end if
          !
        end if
        !
      else
        !
        ! If this is anything but a periodic BC, we need to use bface(5) to
        ! extract the name of the BC for this face from bc_names to see if
        ! it is the same as the name of the current bc_input, meaning this
        ! face is part of the current bc_input. If it is:
        !   1) add this face to the count of total boundary faces for the
        !      current bc_input
        !   2) save the index for this boundary condition within bc_in to
        !      bface(6) for this face
        !   3) save the index for this boundary condition within bc_in to
        !      bface(11) for this face
        !      NOTE: this is only done because periodic BCs use bface(6)
        !      to store periodic connectivity information and instead
        !      use bface(11) to store their bc_in index, and we just
        !      want to make it so all faces can access their bc_in index
        !      from bface(11).
        !   4) if the BC type for this BC has been set in the input file,
        !      i.e., bc_input(n)%set_bc_value /= bc_unknown,
        !      reset bface(1) to this BC type
        !   5) save the BC type to bc_in if it still bc_unknown
        !      
        bcnam = trim(adjustl( uppercase( bc_names(bface(5,nf)) ) ))
        bcinp = trim(adjustl( uppercase( bc_input(n)%name      ) ))
        !
        if (bcnam == bcinp) then
          !
          ii = ii + 1
          !
          ! Save the index for this boundary condition
          ! within the bc_in array to bface(6,:)
          !
          bface(6,nf) = i
          !
          bface(11,nf) = i
          !
          ! Force the boundary condition to that specified in the input file
          ! if the value of set_bc_value is not equal to bc_unknown
          !
          if (bc_input(n)%set_bc_value /= bc_unknown) then
            bface(1,nf) = bc_input(n)%set_bc_value
          end if
          !
          ! Save the boundary condition type to bc_in
          !
          if (bc_in(i)%bc_type == bc_unknown) then
            bc_in(i)%bc_type = bface(1,nf)
          end if
          !
        end if
        !
      end if
      !
    end do bface_loop
    !
#ifdef DEBUG_ON
    call write_bc_conditions(ibc_end=i,dimensional=true,num_bc_faces=ii)
#endif
    !
  end do
  !
  ! For periodic boundary faces, mod(bface(11,:),10000) still gives the
  ! corresponding index for bc_input. We need to change this to the
  ! corresponding index for bc_in.
  !
  ! The array bc_list maps the bc_in index, i, to the corresponding bc_input
  ! index, n, i.e., bc_input(n) == bc_input(bc_list(i)). We need to create
  ! a reverse mapping array that maps the bc_input index, n, to the
  ! corresponding bc_in index, i, i.e., bc_in(i) == bc_in(rev_bc_list(n)).
  !
  rev_bc_list = init_array(size(bc_input),0) ! F2003 AUTO-REALLOC
  do i = 1,size(bc_list)
    rev_bc_list(bc_list(i)) = i
  end do
  !
  ! Use the reverse bc_list array to convert the bc_input index saved in
  ! bface(11) for periodic boundary faces, to the corresponding bc_in index.
  !
  do nf = 1,size(bface,dim=2)
    !
    if (bface(1,nf) == bc_periodic) then
      !
      bc_idx = mod(bface(11,nf),10000)
      !
      if (rev_bc_list(bc_idx) > 0) then
        bface(11,nf) = bface(11,nf) - bc_idx + rev_bc_list(bc_idx)
      end if
      !
    end if
    !
  end do
  !
  !
  !
  if (ncpu > 1) then
    if (mypnum == glb_root) flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," Setting boundary condition for BC # ",i0," - '",a,"'",/, &
              "    to ",a)
  2 format ("Error iterating on BC # ",i0," !!!")
  3 format (/, &
   2x,67("#"),/, &
   2x,"###",26x,"WARNING!!",26x,"###",/, &
   2x,"###  THE SIMULATION APPEARS TO BE SOLVING THE EULER EQUATIONS,  ###",/, &
   2x,"###  BUT SOME BOUNDARY FACES WERE GIVEN NO-SLIP WALL BOUNDARY   ###",/, &
   2x,"###  CONDITIONS. THE BOUNDARY CONDITIONS FOR THESE FACES HAVE   ###",/, &
   2x,"###  BEEN RESET TO A SLIP WALL BOUNDARY CONDITION.              ###",/, &
   2x,"###",26x,"WARNING!!",26x,"###",/, &
   2x,67("#"),/)
  !
end subroutine create_bc_in
!
!###############################################################################
!
pure &
function get_bc_rotation(bc) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : bc_input_t,aref
  use ovar, only : moving_ref_frame
  !
  !.. Formal Arguments ..
  type(bc_input_t), intent(in) :: bc
  !
  !.. Function Result ..
  real(wp) :: return_value(1:3)
  !
continue
  !
  return_value(:) = zero
  !
  ! The values in bc%rotation_rate(:) should be in RPMs. Only copy the
  ! rotation speeds into return_value if greater than 1 RPM.
  !
  if (any(abs(bc%rotation_rate(:)) > one)) then
    !
    ! When saving the wall rotation speed into return_value:
    !   1. multiply by rpm_to_mps to convert from RPM to rad/s, and
    !   2. divide through by aref to normalize by the speed of sound
    ! NOTE: We are normalizing by aref now so that we dont have to do this
    !       later when the rad/s values in the bc_in%rotation_rate vector
    !       is converted to m/s by multiplying by the position vector.
    !
    return_value(:) = bc%rotation_rate(:) * (rpm_to_mps/aref)
    !
  end if
  !
  ! If we are using a moving reference frame and the rotation_rate in bc
  ! is specified in the moving reference frame, we need to convert the
  ! rotation in return_value to the absolute/inertial reference frame.
  !
  if ( bc%is_specified_in_mrf ) then
    if ( moving_ref_frame%is_active() ) then
      !
      ! The rotation in return_value is currently in the MRF, so we need to
      ! add the MRF rotation to it to convert to the absolute/inertial RF.
      ! NOTE: At this point, moving_ref_frame%rotation(:) should have
      !       already been nondimensionalized from RPM to (rad/s)/aref.
      !
      return_value(:) = return_value(:) + moving_ref_frame%rotation(:)
      !
    end if
  end if
  !
end function get_bc_rotation
!
!###############################################################################
!
pure &
function get_bc_velocity(bc,dimensional) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : bc_input_t,aref
  use ovar,   only : moving_ref_frame
  !
  !.. Formal Arguments ..
  type(bc_input_t), intent(in) :: bc
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: dimensional
  !
  !.. Function Result ..
  real(wp) :: return_value(1:nr)
  !
  !.. Local Scalars ..
  real(wp) :: bc_nondim,mrf_nondim
  !
  !.. Local Arrays ..
  real(wp) :: vel(1:3)
  !
continue
  !
  return_value(:) = zero
  !
  ! Get the nondimensionalization coefficients depending on if the
  ! optional argument dimensional is present with the value true.
  ! NOTE: This logic makes a nondimensional result the default, and it is
  !       assumed that the values in bc%vx, bc%vy, and bc%vz are
  !       dimensional, and the values in moving_ref_frame%translation
  !       are nondimensional.
  !
  bc_nondim = one/aref
  mrf_nondim = one
  if (present(dimensional)) then
    if (dimensional) then
      bc_nondim = one
      mrf_nondim = aref
    end if
  end if
  !
  ! Copy the values in bc%vx, bc%vy, and bc%vz to return_value.
  !
  vel(1:3) = [ bc%vx , bc%vy , bc%vz ]
  !
  return_value(1:nr) = vel(1:nr) * bc_nondim
  !
  ! If we are using a moving reference frame and the velocity in bc
  ! is specified in the moving reference frame, we need to convert the
  ! velocity in return_value to the absolute/inertial reference frame.
  !
  if ( bc%is_specified_in_mrf ) then
    if ( moving_ref_frame%is_active() ) then
      !
      ! The velocity in return_value is currently in the MRF,
      ! so we need to add the MRF translational velocity to it
      ! to convert to the absolute/inertial RF.
      ! NOTE: At this point, moving_ref_frame%translation(:) should have
      !       already been nondimensionalized.
      !
      return_value(:) = return_value(:) + &
                        moving_ref_frame%translation(:) * mrf_nondim
      !
    end if
  end if
  !
end function get_bc_velocity
!
!###############################################################################
!
subroutine write_bc_conditions(ibc_beg,ibc_end,dimensional,num_bc_faces)
  !
  !.. Use Statements ..
  use geovar,           only : nr,bface,face
  use ovar,             only : bc_input,bc_in
  use ovar,             only : tmin,tmax,tref
  use ovar,             only : rhoref,aref,cpref_nd
  use ovar,             only : governing_equations
  use ovar,             only : moving_ref_frame
  use metrics_mod,      only : facenrm
  use element_type_mod, only : std_elem
  use thermo_mod,       only : species_name
  !
  !.. Optional Arguments ..
  integer,     optional, intent(in) :: ibc_beg
  integer,     optional, intent(in) :: ibc_end
  logical(lk), optional, intent(in) :: dimensional
  integer,     optional, intent(in) :: num_bc_faces
  !
  !.. Local Scalars ..
  integer  :: ii,k,m,nf,nfp,f0,idx
  integer  :: this_bc,bc_type,bc_beg,bc_end,bc_idx
  real(wp) :: r,p,t,h,r0,p0,t0,h0,mach,aspd,visc,cond,reyn
  real(wp) :: vmag,alpha,beta,total_cond,gamm,rgas,cp
  real(wp) :: gm1,gam_over_gm1,one_over_gm1
  real(wp) :: rdim,adim,tdim,face_area,rpm_dim,one_rpm_nd
  !
  character(len=5)  :: cdim
  character(len=10) :: r_units,t_units,p_units,v_units,h_units
  character(len=10) :: cp_units,visc_units,cond_units,q_units,rpm_units
  !
  character(len=2000) :: isowall_fmt,adiwall_fmt,rhovel_fmt
  character(len=2000) :: isorotwall_fmt,adirotwall_fmt,mrfrotwall_fmt
  character(len=2000) :: thermo_fmt,transport_fmt
  character(len=2000) :: species_fmt,tmp_fmt,energy_fmt
  !
  !.. Local Arrays ..
  real(wp) :: vel(1:nr)
  real(wp) :: mass_frac(1:ncs)
  real(wp) :: rotation(1:3)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_bc_conditions"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  face_area = zero
  !
  ! Determine whether to output dimensional or nondimensional quantities for
  ! each boundary condition based on the optional argument 'dimensional'.
  !
  rdim = rhoref
  adim = aref
  tdim = tref
  rpm_dim = adim/rpm_to_mps
  !
  one_rpm_nd = one/rpm_dim
  !
  cdim = "D"
  t_units = " K"
  v_units = " m/s"
  r_units = " kg/m^3"
  p_units = " N/m^2"
  visc_units = " kg/(m*s)"
  cond_units = " W/(m*K)"
  cp_units = " J/(kg*K)"
  h_units = " J/kg"
  q_units = " m^2/s^2"
  rpm_units = " RPM"
  !
  if (present(dimensional)) then
    if (.not. dimensional) then
      !
      rdim = one
      adim = one
      tdim = one
      rpm_dim = one
      !
      cdim = "Non-d"
      t_units = ""
      v_units = ""
      r_units = ""
      p_units = ""
      visc_units = ""
      cond_units = ""
      cp_units = ""
      h_units = ""
      q_units = ""
      rpm_units = " rad/s"
      !
    end if
  end if
  !
  ! Create the format strings
  !
  write (isowall_fmt,102) t_units,v_units,v_units,v_units
  write (adiwall_fmt,103) v_units,v_units,v_units
  write (isorotwall_fmt,202) t_units,rpm_units,rpm_units,rpm_units
  write (adirotwall_fmt,203) rpm_units,rpm_units,rpm_units
  write (mrfrotwall_fmt,204) rpm_units,rpm_units,rpm_units
  write (rhovel_fmt,104) r_units,v_units,v_units,v_units
  write (thermo_fmt,105) p_units,t_units,h_units,v_units,r_units, &
                         p_units,t_units,h_units,cp_units,cp_units
  write (transport_fmt,106) visc_units,cond_units
  write (energy_fmt,109) p_units,p_units,p_units,p_units,q_units
  !
  if (ncs > 0) then
    !
    write (tmp_fmt,107) (trim(adjustl(species_name(m))), m=1,ncs)
    species_fmt = trim(adjustl(tmp_fmt))
    !
    write (tmp_fmt,108) trim(adjustl(species_fmt))
    species_fmt = trim(adjustl(tmp_fmt))
    !
  end if
  !
  ! Determine what boundary conditions from bc_in to output
  ! NOTE: If only one of ibc_beg and ibc_end is present:
  !            bc_beg = bc_end = value of whichever one is present
  !       If both ibc_beg and ibc_end are present:
  !            bc_beg = ibc_beg
  !            bc_end = ibc_end
  !       If neither ibc_beg or ibc_end are present:
  !            bc_beg = lbound(bc_in,dim=1)
  !            bc_end = ubound(bc_in,dim=1)
  !
  bc_beg = -1
  bc_end = -1
  !
  if (present(ibc_beg)) then
    bc_beg = ibc_beg
  end if
  !
  if (present(ibc_end)) then
    bc_end = ibc_end
  end if
  !
  if (all([bc_beg,bc_end] < 0)) then
    bc_beg = lbound(bc_in,dim=1)
    bc_end = ubound(bc_in,dim=1)
  else if (bc_beg < 0) then
    bc_beg = bc_end
  else if (bc_end < 0) then
    bc_end = bc_beg
  end if
  !
  ! Output the boundary condition(s)
  !
  do this_bc = bc_beg,bc_end
    !
    if (mypnum == glb_root) then
      !
      ! Write out the header for this boundary condition
      !
      if (this_bc == 0) then
        !
        write (iout,1) trim(adjustl(cdim))
        bc_type = bc_unknown
        !
      else
       !!
       !! FOR MULTIPLE MPI PROCESSES, THERE IS NO GUARANTEE THAT THE
       !! THE GLOBAL ROOT PROCESS WILL HAVE A BOUNDARY FACE SUCH THAT
       !!    bface(1,nf) == this_bc
       !!
       !bc_type_loop: do nf = 1,size(bface,dim=2)
       !  if (any(bface(1,nf) == bc_comm) .or. bface(5,nf) == 0) then
       !    cycle bc_type_loop
       !  end if
       !  if (bface(6,nf) == this_bc) then
       !    bc_type = bface(1,nf)
       !    exit bc_type_loop
       !  end if
       !end do bc_type_loop
        !
        bc_idx = bc_in(this_bc)%bc_input_idx
        bc_type = bc_in(this_bc)%bc_type
        !
        write (iout,2) trim(adjustl(cdim)), this_bc, &
                       trim(adjustl(bc_input(bc_idx)%name)), &
                       trim(adjustl(bc_integer_to_string(bc_type)))
        !
      end if
      !
      ! Write out the flow conditions for this boundary condition
      !
      if (bc_type == bc_isothermal_wall) then
        !
        t = bc_in(this_bc)%wall_temperature * tdim
        !
        if (any(abs(bc_in(this_bc)%rotation_rate(:)) >= one_rpm_nd)) then
          rotation(:) = bc_in(this_bc)%rotation_rate(:)
          write (iout,isorotwall_fmt) t, (rotation(ii)*rpm_dim,ii=1,nr)
          if ( moving_ref_frame%is_active() ) then
            rotation(:) = rotation(:) - moving_ref_frame%rotation(:)
            write (iout,mrfrotwall_fmt) (rotation(ii)*rpm_dim,ii=1,nr)
          end if
        else
          vel(1:nr) = bc_in(this_bc)%pv(nmb:nme) * adim
          write (iout,isowall_fmt) t,(vel(ii),ii=1,nr)
        end if
        !
      else if (bc_type == bc_adiabatic_wall) then
        !
        if (any(abs(bc_in(this_bc)%rotation_rate(:)) >= one_rpm_nd)) then
          rotation(:) = bc_in(this_bc)%rotation_rate(:)
          write (iout,adirotwall_fmt) (rotation(ii)*rpm_dim,ii=1,nr)
          if ( moving_ref_frame%is_active() ) then
            rotation(:) = rotation(:) - moving_ref_frame%rotation(:)
            write (iout,mrfrotwall_fmt) (rotation(ii)*rpm_dim,ii=1,nr)
          end if
        else
          vel(1:nr) = bc_in(this_bc)%pv(nmb:nme) * adim
          write (iout,adiwall_fmt) (vel(ii),ii=1,nr)
        end if
        !
      else
        !
        r = bc_in(this_bc)%cv(nec) * rdim
        p = bc_in(this_bc)%gp(npres) * rdim*adim*adim
        t = bc_in(this_bc)%gp(ntemp) * tdim
        gamm = bc_in(this_bc)%gp(ngamm)
        !
        gm1 = gamm - one
        one_over_gm1 = one / gm1
        gam_over_gm1 = gamm / gm1
        !
        if (ncs > 0) then
          !
          do m = 1,ncs
            mass_frac(m) = bc_in(this_bc)%pv(m)
          end do
          !
          h = bc_in(this_bc)%gp(nenth) * adim*adim
          cp = bc_in(this_bc)%gp(ncp) * adim*adim/tdim
          !
        else
          !
          cp = cpref_nd * adim*adim/tdim
          h = cp*t * adim*adim
          !
        end if
        !
        rgas = cp*gm1/gamm
        !
        vel(1:nr) = bc_in(this_bc)%pv(nmb:nme) * adim
        vmag = norm2(vel(1:nr))
        !
        aspd = sqrt(gamm*p/r)
        mach = vmag/aspd
        !
        total_cond = one + half*gm1*mach*mach
        !
        r0 = r * total_cond**one_over_gm1
        p0 = p * total_cond**gam_over_gm1
        t0 = t * total_cond
        h0 = h + half*vmag*vmag
        !
        alpha = zero
        beta = zero
        if (nr >= 2) then
          alpha = atan2(bc_in(this_bc)%pv(nmy),bc_in(this_bc)%pv(nmx))
          if (nr == 3) then
            if (vmag > eps6) then
              beta = asin(bc_in(this_bc)%pv(nmz)/(vmag/adim))
            end if
          end if
        end if
        alpha = alpha * radians2degrees
        beta = beta * radians2degrees
        !
        if (ncs > 0) then
          write (iout,species_fmt) (mass_frac(m),m=1,ncs)
        end if
        !
        write (iout,rhovel_fmt) r,(vel(ii),ii=1,nr)
        write (iout,thermo_fmt) p,t,h,gamm,mach,aspd, &
                                r0,p0,t0,h0,cp,rgas,alpha,beta
        write (iout,energy_fmt) r*h,p,half*r*vmag*vmag, &
                                r*h-p+half*r*vmag*vmag, &
                                vmag*vmag
        !
        if (governing_equations == NavierStokes_Eqns) then
          !
          visc = bc_in(this_bc)%gp(nvisc) * rdim*adim
          cond = bc_in(this_bc)%gp(ncond) * rdim*adim*adim*adim/tdim
          reyn = vmag*r/visc
          !
          write (iout,transport_fmt) visc,cond,reyn
          !
        end if
        !
      end if
      !
      flush (iout)
      !
    end if
    !
    if (present(num_bc_faces)) then
      !
      ii = num_bc_faces
      !
    else
      !
      ii = 0
      face_area = zero
      bface_loop: do nf = 1,size(bface,dim=2)
        !
        ! If one of the conditionals hit, prematurely exit the one-pass loop
        ! which will continue on with adding the current face to the boundary
        ! faces for this boundary condition. If none of the conditions hit,
        ! the end of the one-pass loop will be reached which contains a
        ! cycle statement for the bface loop, so it will continue on to the
        ! next boundary face without adding the current face.
        !
        one_pass_loop: &
        do idx = 1,1
          !
         !if (any(bface(1,nf) == bc_comm) .or. bface(5,nf)==0) cycle bface_loop
          if (bface(5,nf) == 0) cycle bface_loop
          !
          ! Check for non-communication BC
          !
          if (all(bface(1,nf) /= bc_comm)) then
            if (bface(6,nf) == this_bc) then
              exit one_pass_loop
            end if
          end if
          !
          ! Check for local periodic BC [bface(1)==bc_periodic] or
          ! MPI communication periodic BC [bface(1)=bc_cpu_bnd.and.bface(11)>0]
          !
          if ( (bface(1,nf) == bc_periodic) .or. &
               (bface(1,nf) == bc_cpu_bnd .and. bface(11,nf) > 9999) ) then
            if (mod(bface(11,nf),10000) == this_bc) then
              exit one_pass_loop
            end if
          end if
          !
          cycle bface_loop
          !
        end do one_pass_loop
        !
        ii = ii + 1
        if (allocated(face) .and. allocated(facenrm)) then
          f0 = face(nf)%beg_fp - 1
          nfp = face(nf)%end_fp - f0
          do k = 1,nfp
            face_area = face_area + norm2(facenrm(:,f0+k))* &
                        std_elem( face(nf)%geom , face(nf)%order )%wts(k)
          end do
        end if
        !
      end do bface_loop
      !
      call mpi_allreduce(MPI_IN_PLACE,ii,1_int_mpi,mpi_inttyp, &
                         MPI_SUM,MPI_COMM_WORLD,mpierr)
      !
      call mpi_allreduce(MPI_IN_PLACE,face_area,1_int_mpi,mpi_flttyp, &
                         MPI_SUM,MPI_COMM_WORLD,mpierr)
      !
    end if
    !
    if (mypnum == glb_root) then
      write (iout,7) ii
      write (iout,8) face_area
      flush (iout)
    end if
    !
  end do
  !
  if (ncpu > 1) then
    if (mypnum == glb_root) flush (iout)
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  102 format ('(4x,"Wall Temperature",t27,"= ",es14.6,"',a,'",/,', &
               '4x,"Wall X-velocity",t27,"= ",es14.6,"',a,'",/,:,', &
               '4x,"Wall Y-velocity",t27,"= ",es14.6,"',a,'",/,:,', &
               '4x,"Wall Z-velocity",t27,"= ",es14.6,"',a,'",/)')
  103 format ('(4x,"Wall X-velocity",t27,"= ",es14.6,"',a,'",/,:,', &
               '4x,"Wall Y-velocity",t27,"= ",es14.6,"',a,'",/,:,', &
               '4x,"Wall Z-velocity",t27,"= ",es14.6,"',a,'",/)')
  202 format ('(4x,"Wall Temperature",t27,"= ",es14.6,"',a,'",/,', &
               '4x,"Wall X-Rotation IRF",t27,"= ",es14.6,"',a,'",/,:,', &
               '4x,"Wall Y-Rotation IRF",t27,"= ",es14.6,"',a,'",/,:,', &
               '4x,"Wall Z-Rotation IRF",t27,"= ",es14.6,"',a,'",/)')
  203 format ('(4x,"Wall X-Rotation IRF",t27,"= ",es14.6,"',a,'",/,:,', &
               '4x,"Wall Y-Rotation IRF",t27,"= ",es14.6,"',a,'",/,:,', &
               '4x,"Wall Z-Rotation IRF",t27,"= ",es14.6,"',a,'",/)')
  204 format ('(4x,"Wall X-Rotation MRF",t27,"= ",es14.6,"',a,'",/,:,', &
               '4x,"Wall Y-Rotation MRF",t27,"= ",es14.6,"',a,'",/,:,', &
               '4x,"Wall Z-Rotation MRF",t27,"= ",es14.6,"',a,'",/)')
  104 format ('(4x,"Static Density",t27,"= ",es14.6,"',a,'",/,' &
               '4x,"X-velocity",t27,"= ",es14.6,"',a,'",:,/,' &
               '4x,"Y-velocity",t27,"= ",es14.6,"',a,'",:,/,' &
               '4x,"Z-velocity",t27,"= ",es14.6,"',a,'")')
  105 format ('(4x,"Static Pressure",t27,"= ",es14.6,"',a,'",/,' &
               '4x,"Static Temperature",t27,"= ",es14.6,"',a,'",/,' &
               '4x,"Static Enthalpy",t27,"= ",es14.6,"',a,'",/,' &
               '4x,"Gamma",t27,"= ",es14.6,/,' &
               '4x,"Mach Number",t27,"= ",es14.6,/,' &
               '4x,"Speed of Sound",t27,"= ",es14.6,"',a,'",/,' &
               '4x,"Total Density",t27,"= ",es14.6,"',a,'",/,' &
               '4x,"Total Pressure",t27,"= ",es14.6,"',a,'",/,' &
               '4x,"Total Temperature",t27,"= ",es14.6,"',a,'",/,' &
               '4x,"Total Enthalpy",t27,"= ",es14.6,"',a,'",/,' &
               '4x,"Cp",t27,"= ",es14.6,"',a,'",/,' &
               '4x,"Gas Constant",t27,"= ",es14.6,"',a,'",/,' &
               '4x,"Alpha Angle of Attack",t27,"= ",es14.6," degrees",/,' &
               '4x,"Beta  Angle of Attack",t27,"= ",es14.6," degrees")')
  106 format ('(4x,"Molecular Viscosity",t27,"= ",es14.6,"',a,'",/,' &
               '4x,"Thermal Conductivity",t27,"= ",es14.6,"',a,'",/,' &
               '4x,"Reynolds #",t27,"= ",es14.6,/)')
  107 format (100('4x,"',a,' Mass Fraction",t27,"= ",es14.6',:,',/,'))
  108 format ('(',a,')')
  109 format ('(4x,"rho*h",t27,"= ",es14.6,"',a,'",/,' &
               '4x,"Static Pressure",t27,"= ",es14.6,"',a,'",/,' &
               '4x,"0.5*rho*q^2",t27,"= ",es14.6,"',a,'",/,' &
               '4x,"rhoE",t27,"= ",es14.6,"',a,'",/,' &
               '4x,"q^2",t27,"= ",es14.6,"',a,'")')
  !
  1 format (/,2x,a,"imensional flow conditions for BC # 0 - ", &
                   "Reference/Freestream BCs")
  2 format (/,2x,a,"imensional flow conditions for BC # ",i0, &
                   " - '",a,"'",/, &
              4x,"BC Type = ",a)
 !3 format (4x,"Wall Temperature      = ",es14.6," K",/, &
 !          4x,"X-velocity            = ",es14.6," m/s",/,:, &
 !          4x,"Y-velocity            = ",es14.6," m/s",/,:, &
 !          4x,"Z-velocity            = ",es14.6," m/s",/)
 !4 format (4x,"Static Density        = ",es14.6," kg/m^3",/, &
 !          4x,"X-velocity            = ",es14.6," m/s",:,/, &
 !          4x,"Y-velocity            = ",es14.6," m/s",:,/, &
 !          4x,"Z-velocity            = ",es14.6," m/s")
 !5 format (4x,"Static Pressure       = ",es14.6," N/m^2",/, &
 !          4x,"Static Temperature    = ",es14.6," K",/, &
 !          4x,"Mach Number           = ",es14.6,/, &
 !          4x,"Speed of Sound        = ",es14.6," m/s",/, &
 !          4x,"Total Density         = ",es14.6," kg/m^3",/, &
 !          4x,"Total Pressure        = ",es14.6," N/m^2",/, &
 !          4x,"Total Temperature     = ",es14.6," K",/, &
 !          4x,"Alpha Angle of Attack = ",es14.6," degrees",/, &
 !          4x,"Beta  Angle of Attack = ",es14.6," degrees",/, &
 !          4x,"Viscosity             = ",es14.6," kg/(m*s)",/, &
 !          4x,"Reynolds #            = ",es14.6,/)
  7 format (4x,"Number of boundary faces for this BC = ",i0)
  8 format (4x,"               Face area for this BC = ",es18.10)
  !
end subroutine write_bc_conditions
!
!###############################################################################
!
subroutine write_linear_advec_ref_conditions
  !
  !.. Use Statements ..
  use ovar, only : bc_in
  use ovar, only : gam,rgasref,alpha_aoaref
  !
  !.. Local Scalars ..
  real(wp) :: r,p,t,mach,aspd,alpha
  !
  !.. Local Arrays ..
  real(wp) :: vel(1:2)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_linear_advec_ref_conditions"
  !
continue
  !
  r = bc_in(0)%cv(nec)
  p = bc_in(0)%gp(npres)
  t = bc_in(0)%gp(ntemp)
  !
  vel(1) = bc_in(0)%pv(nmx)
  vel(2) = bc_in(0)%pv(nmy)
  !
  aspd = sqrt(gam*p/r)
  mach = norm2(vel)/aspd
  !
  alpha = alpha_aoaref
  !
  if (mypnum == glb_root) then
    !
    write (iout,1) r,p,t,mach,aspd,vel(1),vel(2),alpha,zero,gam,rgasref
    flush (iout)
    !
  end if
  !
  if (ncpu > 1) then
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
  end if
  !
  ! Format Statements
  !
  1 format (/, &
            " The density vortex transport problem is being",/, &
            " used to model the linear 2D advection equation.",//, &
            " The Non-dimensional Free-stream conditions are",/, &
            " Density          = ",es14.6,/, &
            " Pressure         = ",es14.6,/, &
            " Temperature      = ",es14.6,/, &
            " Mach Number      = ",es14.6,/, &
            " Speed of Sound   = ",es14.6,/, &
            " X-Velocity       = ",es14.6,/, &
            " Y-Velocity       = ",es14.6,/, &
            " Alpha A.o.Attack = ",es14.6," degrees"/, &
            " Beta  A.o.Attack = ",es14.6," degrees"/, &
            " Gamma            = ",es14.6,/, &
            " Gas Constant     = ",es14.6,/)
  !
end subroutine write_linear_advec_ref_conditions
!
!###############################################################################
!
subroutine dump_bc_conditions(ii,gamm,rgas,r,p,t,r0,p0,t0, &
                              vel,mach,alpha,beta,mass_frac, &
                              static_is_set,total_is_set, &
                              vel_is_set,mach_is_set,string, &
                              use_dump_bc)
  !
  !.. Use Statements ..
  use thermo_mod, only : species_name
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: ii
  real(wp),         intent(in) :: gamm,rgas,r,r0,p,p0,t,t0
  real(wp),         intent(in) :: mach,alpha,beta
  real(wp),         intent(in) :: vel(:)
  real(wp),         intent(in) :: mass_frac(:)
  logical(lk),      intent(in) :: static_is_set,total_is_set
  logical(lk),      intent(in) :: vel_is_set,mach_is_set
  character(len=*), intent(in) :: string
  logical(lk),      intent(in) :: use_dump_bc
  !
  !.. Local Scalars ..
  integer :: i,io
  character(len=10) :: fname
  integer, save :: output = -1
  !
continue
  !
  if (use_dump_bc) then
    !
    if (mypnum == glb_root) then
      !
     !output = output + 1
     !!
     !write (fname,"(i0,'.log')") output
     !!
     !open (newunit=io,file=fname,status="replace",action="write")
     !write (io,1) trim(adjustl(string)),ii, &
     !               gamm,rgas,r,p,t,static_is_set, &
     !               r0,p0,t0,total_is_set, &
     !               sqrt(merge(gamm*p/r,zero,gamm*p/r>zero)), &
     !               mach,mach_is_set
     !write (io,2) (vel(i),i=1,size(vel))
     !write (io,3) alpha,beta,vel_is_set
     !if (ncs > 0) then
     !  do i = 1,ncs
     !    write (io,4) uppercase( trim(adjustl(species_name(i))) ), mass_frac(i)
     !  end do
     !else
     !  write (io,*)
     !end if
     !close (io)
      !
      write (iout,1) trim(adjustl(string)),ii, &
                     gamm,rgas,r,p,t,static_is_set, &
                     r0,p0,t0,total_is_set, &
                     sqrt(merge(gamm*p/r,zero,gamm*p/r>zero)), &
                     mach,mach_is_set
      write (iout,2) (vel(i),i=1,size(vel))
      write (iout,3) alpha,beta,vel_is_set
      if (ncs > 0) then
        do i = 1,ncs
          write (iout,4) uppercase( trim(adjustl(species_name(i))) ), &
                         mass_frac(i)
        end do
      else
        write (iout,*)
      end if
      !
    end if
    !
    flush (iout)
    !
  end if
  !
  ! Format Statements
  !
  1 format (/, &
            3x,"Dumping variables at location '",a, &
               "' for BC iteration # ",i0,/, &
            5x,"Gamma                 = ",es14.6,/, &
            5x,"Gas Constant          = ",es14.6,/, &
            5x,"Static Density        = ",es14.6,/, &
            5x,"Static Pressure       = ",es14.6,/, &
            5x,"Static Temperature    = ",es14.6,/, &
            5x,"static_is_set         = ",l1,/, &
            5x,"Total Density         = ",es14.6,/, &
            5x,"Total Pressure        = ",es14.6,/, &
            5x,"Total Temperature     = ",es14.6,/, &
            5x,"total_is_set          = ",l1,/, &
            5x,"Speed of Sound        = ",es14.6,/, &
            5x,"Mach Number           = ",es14.6,/, &
            5x,"mach_is_set           = ",l1)
  2 format (5x,"X-velocity            = ",es14.6,:,/, &
            5x,"Y-velocity            = ",es14.6,:,/, &
            5x,"Z-velocity            = ",es14.6)
  3 format (5x,"Alpha Angle of Attack = ",es14.6,/, &
            5x,"Beta  Angle of Attack = ",es14.6,/, &
            5x,"vel_is_set            = ",l1)
  4 format (5x,a," Mass Fraction",t28,"= ",es14.6)
  !
end subroutine dump_bc_conditions
!
!###############################################################################
!
subroutine check_bc_input(bc,pv,gp,n)
  !
  !.. Use Statements ..
  use geovar, only : nr
  use ovar,   only : bc_input_t
  !
  !.. Formal Arguments ..
  type(bc_input_t), intent(in) :: bc
  real(wp),         intent(in) :: pv(:)
  real(wp),         intent(in) :: gp(:)
  integer,          intent(in) :: n
  !
  !.. Local Scalars ..
  real(wp)    :: r,p,t,r0,p0,t0,vmag,bc_vel,mach,a,alpha,beta
  real(wp)    :: gamm,gm1,one_over_gm1,gam_over_gm1,err_tol
  real(wp)    :: total_cond
  logical(lk) :: errors_found
  !
  !.. Local Arrays ..
  real(wp) :: vel(1:nr)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_bc_input"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  errors_found = fals
  err_tol = eps6
  !
  !
  r = pv(nec)
  vel(:) = pv(nmb:nme)
  vmag = norm2(vel)
  !
  p = gp(npres)
  t = gp(ntemp)
  gamm = gp(ngamm)
  !
  a = sqrt(gamm*p/r)
  mach = vmag/a
  !
  gm1 = gamm - one
  one_over_gm1 = one / gm1
  gam_over_gm1 = gamm / gm1
  !
  total_cond = one + half*gm1*mach*mach
  r0 = r * (total_cond**one_over_gm1)
  p0 = p * (total_cond**gam_over_gm1)
  t0 = t *  total_cond
  !
  alpha = zero
  beta = zero
  if (size(vel) == 3) beta = asin(vel(3)/vmag) * radians2degrees
  if (size(vel) >= 2) alpha = atan2(vel(2),vel(1)) * radians2degrees
  !
  ! Check static temperature
  !
  if (bc%t_static > zero) then
    if (abs(bc%t_static - t) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,1) bc%t_static,t
      end if
      errors_found = true
    end if
  end if
  !
  ! Check static pressure
  !
  if (bc%p_static > zero) then
    if (abs(bc%p_static - p) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,2) bc%p_static,p
      end if
      errors_found = true
    end if
  end if
  !
  ! Check static density
  !
  if (bc%rho_static > zero) then
    if (abs(bc%rho_static - r) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,3) bc%rho_static,r
      end if
      errors_found = true
    end if
  end if
  !
  ! Check total temperature
  !
  if (bc%t_total > zero) then
    if (abs(bc%t_total - t0) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,4) bc%t_total,t0
      end if
      errors_found = true
    end if
  end if
  !
  ! Check total pressure
  !
  if (bc%p_total > zero) then
    if (abs(bc%p_total - p0) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,5) bc%p_total,p0
      end if
      errors_found = true
    end if
  end if
  !
  ! Check total density
  !
  if (bc%rho_total > zero) then
    if (abs(bc%rho_total - r0) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,6) bc%rho_total,r0
      end if
      errors_found = true
    end if
  end if
  !
  ! Check mach number
  !
  if (bc%mach > zero) then
    if (abs(bc%mach - mach) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,7) bc%mach,mach
      end if
      errors_found = true
    end if
  end if
  !
  ! Check alpha angle of attack
  !
  if (bc%alpha_aoa /= zero) then
    if (abs(bc%alpha_aoa - alpha) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,8) bc%alpha_aoa,alpha
      end if
      errors_found = true
    end if
  end if
  !
  ! Check beta angle of attack
  !
  if (bc%beta_aoa /= zero) then
    if (abs(bc%beta_aoa - beta) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,9) bc%beta_aoa,beta
      end if
      errors_found = true
    end if
  end if
  !
  ! Check velocity
  !
  bc_vel = norm2( [ bc%vx , bc%vy , bc%vz ] )
  !
  if (bc_vel > zero) then
    if (abs(bc_vel - vmag) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,10) bc_vel,vmag
      end if
      errors_found = true
    end if
    if (abs(bc%vx - vel(1)) > err_tol) then
      if (mypnum == glb_root) then
        if (.not.errors_found) write (iout,100) n
        write (iout,11) bc%vx,vel(1)
      end if
      errors_found = true
    end if
    if (size(vel) >= 2) then
      if (abs(bc%vy - vel(2)) > err_tol) then
        if (mypnum == glb_root) then
          if (.not.errors_found) write (iout,100) n
          write (iout,12) bc%vy,vel(2)
        end if
        errors_found = true
      end if
    end if
    if (size(vel) == 3) then
      if (abs(bc%vz - vel(3)) > err_tol) then
        if (mypnum == glb_root) then
          if (.not.errors_found) write (iout,100) n
          write (iout,13) bc%vz,vel(3)
        end if
        errors_found = true
      end if
    end if
  end if
  !
  if (errors_found) then
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1   format (5x,"Input       static temperature = ",es14.6,/, &
              5x,"Computed    static temperature = ",es14.6)
  2   format (5x,"Input          static pressure = ",es14.6,/, &
              5x,"Computed       static pressure = ",es14.6)
  3   format (5x,"Input           static density = ",es14.6,/, &
              5x,"Computed        static density = ",es14.6)
  4   format (5x,"Input        total temperature = ",es14.6,/, &
              5x,"Computed     total temperature = ",es14.6)
  5   format (5x,"Input           total pressure = ",es14.6,/, &
              5x,"Computed        total pressure = ",es14.6)
  6   format (5x,"Input            total density = ",es14.6,/, &
              5x,"Computed         total density = ",es14.6)
  7   format (5x,"Input              Mach number = ",es14.6,/, &
              5x,"Computed           Mach number = ",es14.6)
  8   format (5x,"Input    alpha angle of attack = ",es14.6,/, &
              5x,"Computed alpha angle of attack = ",es14.6)
  9   format (5x,"Input    beta  angle of attack = ",es14.6,/, &
              5x,"Computed beta  angle of attack = ",es14.6)
  10  format (5x,"Input       velocity magnitude = ",es14.6,/, &
              5x,"Computed    velocity magnitude = ",es14.6)
  11  format (5x,"Input               x-velocity = ",es14.6,/, &
              5x,"Computed            x-velocity = ",es14.6)
  12  format (5x,"Input               y-velocity = ",es14.6,/, &
              5x,"Computed            y-velocity = ",es14.6)
  13  format (5x,"Input               z-velocity = ",es14.6,/, &
              5x,"Computed            z-velocity = ",es14.6)
  100 format (/,3x,"An error was found checking the final ", &
                   "flow conditions for BC # ",i0," !!!",/)
  !
end subroutine check_bc_input
!
!###############################################################################
!
subroutine vel_from_static(gamm,rgas,r,p,t,vel,mach,alpha,beta, &
                           vel_is_set,mach_is_set)
  !
  !.. Formal Arguments ..
  real(wp),       intent(in) :: gamm,rgas,r,p,t
  real(wp),    intent(inout) :: vel(:)
  real(wp),    intent(inout) :: mach
  real(wp),    intent(inout) :: alpha
  real(wp),    intent(inout) :: beta
  logical(lk), intent(inout) :: vel_is_set
  logical(lk), intent(inout) :: mach_is_set
  !
  !.. Local Scalars ..
  real(wp) :: aspd
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "vel_from_static"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  aspd = -huge(zero)
  if (t > zero .and. gamm > zero) then
    aspd = sqrt(gamm*rgas*t)
  else if (r > zero .and. p > zero .and. gamm > zero) then
    aspd = sqrt(gamm*p/r)
  else if (vel_is_set .and. mach_is_set) then
    aspd = norm2(vel) / mach
  end if
  !
  if (vel_is_set) then
    if (size(vel) == 3) beta = asin(vel(3)/norm2(vel)) * radians2degrees
    if (size(vel) >= 2) alpha = atan2(vel(2),vel(1)) * radians2degrees
    if (aspd > zero .and. .not. mach_is_set) then
      mach = norm2(vel)/aspd
      mach_is_set = true
    end if
  else
    if (aspd > zero .and. mach_is_set) then
      vel(:) = mach*aspd
      if (size(vel) >= 2) then
        vel(1) = vel(1)*cos(alpha*degrees2radians)
        vel(2) = vel(2)*sin(alpha*degrees2radians)
      end if
      if (size(vel) == 3) then
        vel(1) = vel(1)*cos(beta *degrees2radians)
        vel(2) = vel(2)*cos(beta *degrees2radians)
        vel(3) = vel(3)*sin(beta *degrees2radians)
      end if
      vel_is_set = true
    end if
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine vel_from_static
!
!###############################################################################
!
subroutine mach_total_temp(gamm,t,t0,mach,mach_is_set,n)
  !
  !.. Formal Arguments ..
  real(wp),       intent(in) :: gamm
  real(wp),    intent(inout) :: t,t0
  real(wp),    intent(inout) :: mach
  logical(lk), intent(inout) :: mach_is_set
  integer,        intent(in) :: n
  !
  !.. Local Scalars ..
  real(wp) :: tmach,total_cond
  real(wp) :: gm1,two_over_gm1
  real(wp) :: tratio
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "mach_total_temp"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  gm1 = gamm - one
  !
  if (t > zero .and. t0 > zero) then
    !
    two_over_gm1 = two / gm1
    !
    tratio = t0/t
    !
    tmach = sqrt( two_over_gm1 * (tratio - one) )
    !
    if (mach_is_set) then
      if (abs(mach-tmach) > eps6) then
        write (error_message,1) n,mach,tmach
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
    else
      mach = tmach
      mach_is_set = true
    end if
    !
  else
    !
    if (mach_is_set) then
      total_cond = one + half*gm1*mach*mach
      if (t > zero) then
        t0 = t * total_cond
      else if (t0 > zero) then
        t = t0 / total_cond
      end if
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ( &
     " Error with BC # ",i0," !", &
     " Mach number specified for boundary condition = ",es14.6, &
     " Mach number computed from the stagnation temperature = ",es14.6, &
     " Please correct the input values for this BC!")
  !
end subroutine mach_total_temp
!
!###############################################################################
!
subroutine mach_total_temp_tpg(mass_frac,rgas,gamm,t,t0,mach)
  !
  !.. Use Statements ..
  use ovar, only : tref,tmin,tmax
  use ovar, only : conv_temp
  !
  use thermo_mod, only : compute_mixture_enthalpy
  use thermo_mod, only : compute_mixture_cp
  use thermo_mod, only : compute_mixture_gamma
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: mass_frac(:)
  real(wp),    intent(in) :: rgas,t0,mach
  real(wp), intent(inout) :: gamm
  real(wp), intent(inout) :: t
  !
  !.. Local Scalars ..
  integer  :: k
  real(wp) :: h0,h,cp
  real(wp) :: t0_nd,t_nd
  real(wp) :: func,dfunc,dtemp,dterr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "mach_total_temp_tpg"
  !
  integer,  parameter :: newton_iter = 100
 !real(wp), parameter :: conv_temp = eps6
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Compute total enthalpy
  !
  t0_nd = t0/tref
  call compute_mixture_enthalpy(mass_frac,t0_nd,h0,dimensional=true)
  !
  ! Use the total temperature to get an initial guess for gamma
  !
  call compute_mixture_gamma(mass_frac,t0_nd,gamm,dimensional=true)
  !
  t = t0 / (one + half*(gamm-one)*mach*mach)
  !
  ! Compute the static temperature using Newtons method
  !
  newton_iteration: &
  do k = 1,newton_iter
    !
    t_nd = t/tref
    call compute_mixture_enthalpy(mass_frac,t_nd,h,dimensional=true)
    !
    call compute_mixture_cp(mass_frac,t_nd,cp,dimensional=true)
    !
    gamm = rgas
    call compute_mixture_gamma(cp,gamm)
    !
    func = h0 - h - half*gamm*rgas*t*mach*mach
    dfunc = cp + half*gamm*rgas*mach*mach
    !
    dtemp = func/dfunc
    !
    t = min( tmax , max( t-dtemp , tmin ) )
    !
    dterr = abs(dtemp) / t
    !
    if (dterr < conv_temp) then
      exit newton_iteration
    end if
    !
  end do newton_iteration
  !
  ! Update the mixture gamma using the final temperature
  !
  t_nd = t/tref
  call compute_mixture_gamma(mass_frac,t_nd,gamm,dimensional=true)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine mach_total_temp_tpg
!
!###############################################################################
!
subroutine mach_total_pres(gamm,p,p0,mach,mach_is_set,n)
  !
  !.. Formal Arguments ..
  real(wp),       intent(in) :: gamm
  real(wp),    intent(inout) :: p,p0
  real(wp),    intent(inout) :: mach
  logical(lk), intent(inout) :: mach_is_set
  integer,        intent(in) :: n
  !
  !.. Local Scalars ..
  real(wp) :: pmach,total_cond
  real(wp) :: two_over_gm1,gm1
  real(wp) :: gam_over_gm1,pratio
  real(wp) :: gm1_over_gam
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "mach_total_pres"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  gm1 = gamm - one
  !
  if (p > zero .and. p0 > zero) then
    !
    gm1_over_gam = gm1 / gamm
    two_over_gm1 = two / gm1
    !
    pratio = (p0/p)**gm1_over_gam
    !
    pmach = sqrt( two_over_gm1 * (pratio - one) )
    !
    if (mach_is_set) then
      if (abs(mach-pmach) > eps6) then
        write (error_message,1) n,mach,pmach
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
    else
      mach = pmach
      mach_is_set = true
    end if
    !
  else
    !
    if (mach_is_set) then
      gam_over_gm1 = gamm / gm1
      total_cond = (one + half*gm1*mach*mach)**gam_over_gm1
      if (p > zero) then
        p0 = p * total_cond
      else if (p0 > zero) then
        p = p0 / total_cond
      end if
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ( &
     " Error with BC # ",i0," !", &
     " Mach number specified for boundary condition = ",es14.6, &
     " Mach number computed from the stagnation pressure = ",es14.6, &
     " Please correct the input values for this BC!")
  !
end subroutine mach_total_pres
!
!###############################################################################
!
subroutine mach_total_dens(gamm,r,r0,mach,mach_is_set,n)
  !
  !.. Formal Arguments ..
  real(wp),       intent(in) :: gamm
  real(wp),    intent(inout) :: r,r0
  real(wp),    intent(inout) :: mach
  logical(lk), intent(inout) :: mach_is_set
  integer,        intent(in) :: n
  !
  !.. Local Scalars ..
  real(wp) :: rmach,total_cond
  real(wp) :: two_over_gm1,gm1
  real(wp) :: one_over_gm1,rratio
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "mach_total_dens"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  gm1 = gamm - one
  !
  if (r > zero .and. r0 > zero) then
    !
    two_over_gm1 = two / gm1
    rratio = (r0/r)**gm1
    !
    rmach = sqrt( two_over_gm1 * (rratio - one) )
    !
    if (mach_is_set) then
      if (abs(mach-rmach) > eps6) then
        write (error_message,1) n,mach,rmach
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
    else
      mach = rmach
      mach_is_set = true
    end if
    !
  else
    !
    if (mach_is_set) then
      one_over_gm1 = one / gm1
      total_cond = (one + half*gm1*mach*mach)**one_over_gm1
      if (r > zero) then
        r0 = r * total_cond
      else if (r0 > zero) then
        r = r0 / total_cond
      end if
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ( &
     " Error with BC # ",i0," !", &
     " Mach number specified for boundary condition = ",es14.6, &
     " Mach number computed from the stagnation density = ",es14.6, &
     " Please correct the input values for this BC!")
  !
end subroutine mach_total_dens
!
!###############################################################################
!
subroutine igl(rgas,rho,pres,temp,success,n)
  !
  !.. Formal Arguments ..
  real(wp),       intent(in) :: rgas
  real(wp),    intent(inout) :: rho
  real(wp),    intent(inout) :: pres
  real(wp),    intent(inout) :: temp
  logical(lk), intent(inout) :: success
  integer,        intent(in) :: n
  !
  !.. Local Scalars ..
  real(wp) :: pres_check
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "igl"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  success = fals
  !
  if (count([rho,pres,temp] > zero) >= 2) then
    !
    if ( rho > zero .and. pres > zero .and. temp > zero ) then
      !
      pres_check = rho * rgas * temp
      !
      if (abs(pres_check-pres) > eps6) then
        write (error_message,1) n,pres,pres_check
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      !
    else if ( rho < zero .and. pres > zero .and. temp > zero ) then
      rho = pres / (rgas * temp)
    else if ( pres < zero .and. rho > zero .and. temp > zero ) then
      pres = rho * rgas * temp
    else if ( temp < zero .and. rho > zero .and. pres > zero ) then
      temp = pres / (rgas * rho)
    end if
    !
    success = true
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (" Error with BC # ",i0," !", &
            " Pressure specified for boundary condition = ",es14.6, &
            " Pressure computed from other specified conditions = ",es14.6, &
            " Please correct the flow conditions for this BC in the input", &
            " file or set one of density, pressure, and temperature to a", &
            " negative value so it is explicitly computed.")
  !
end subroutine igl
!
!###############################################################################
!
subroutine compute_flow_conditions_from_pvt(cv,pv,gp,dimensional)
  !
  ! Given pressure, velocity, and temperature (and species mass fractions if
  ! this is a thermally perfect gas mixture), compute the remaining flow
  ! variables in the cv, pv, and gp arrays.
  !
  !.. Use Statements ..
  use ovar, only : gam,tref,cpref,cpref_nd
  use ovar, only : rgasref,rgasref_nd
  !
  use thermo_mod,    only : compute_mixture_rgas
  use thermo_mod,    only : compute_mixture_gamma
  use thermo_mod,    only : compute_mixture_enthalpy
  use thermo_mod,    only : compute_mixture_cp
  use transport_mod, only : update_transport_properties
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: cv(:)
  real(wp), intent(inout) :: pv(:)
  real(wp), intent(inout) :: gp(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: dimensional
  !
  !.. Local Scalars ..
  integer  :: m
  real(wp) :: rho,q2,pres,temp,temp_nd,enth,cp,rgas
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_flow_conditions_from_pvt"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  pres = gp(npres)
  temp = gp(ntemp)
  !
  q2 = dot_product( pv(nmb:nme) , pv(nmb:nme) )
  !
  if (ncs > 0) then
    !
    call compute_mixture_rgas(pv(1:ncs),rgas,dimensional)
    !
    rho = pres / (rgas*temp)
    !
    do m = 1,ncs
      cv(m) = rho*pv(m)
    end do
    !
    temp_nd = temp
    if (present(dimensional)) then
      if (dimensional) then
        temp_nd = temp/tref
      end if
    end if
    !
    call compute_mixture_enthalpy(pv(1:ncs),temp_nd,gp(nenth),dimensional)
    call compute_mixture_cp(pv(1:ncs),temp_nd,gp(ncp),dimensional)
    !
    gp(ngamm) = rgas
    call compute_mixture_gamma(gp(ncp),gp(ngamm))
    !
    enth = gp(nenth)
    !
  else
    !
    rgas = rgasref_nd
    cp = cpref_nd
    !
    if (present(dimensional)) then
      if (dimensional) then
        rgas = rgasref
        cp = cpref
      end if
    end if
    !
    rho = pres / (rgas*temp)
    !
    enth = cp*temp
    !
  end if
  !
  cv(nec) = rho
  cv(nmb:nme) = rho*pv(nmb:nme)
  cv(nee) = rho*enth - pres + half*rho*q2
  !
  pv(nec) = rho
  pv(nee) = pres
  !
  call update_transport_properties(cv,gp)
  !
  if (ncs == 0) then
    gp(ngamm) = gam
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_flow_conditions_from_pvt
!
!###############################################################################
!
include "Functions/taylor_green.f90"
include "Functions/sod.f90"
#include "Functions/pulse.f90"
include "Functions/smc000_interior_wall_radius.f90"
include "Functions/arn2_interior_wall_radius.f90"
include "Functions/ssjet_interior_wall_radius.f90"
#include "Functions/check_for_bad_number.f90"
!
!###############################################################################
!
end module initialization_mod
