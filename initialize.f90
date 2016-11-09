module initialization_mod
  !
  ! Module containing the functions and routines to
  ! initialize the solution
  !
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,ntk,ntl,nq
  !
  implicit none
  !
  private
  !
  public :: initialize_flow
  public :: initialize_reference
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
  use generic_types_mod, only : matrix
  !
  use geovar, only : nr,n_solpts,ncell,cell,xyz
  !
  use ovar, only : rhoref,aref,itestcase
  use ovar, only : cv_in,iop_solpts
  use ovar, only : irestart,d_t,dt_in,itcur,time
  use ovar, only : dtrst,itrst,tmrst
  use ovar, only : rhoref,pref,tref,gam,rgasref,cpref
  use ovar, only : aref,machref,velref,uref,muref,reyref
  use ovar, only : alpha_aoaref,beta_aoaref
  use ovar, only : pref_nd,velref_nd,muref_nd
  use ovar, only : rgasref_nd,cpref_nd
  !
  use flowvar, only : usp,uoldsp
  use flowvar, only : allocate_flowvar_arrays
  !
  use order_mod, only : maxpts
  !
  use restart_mod, only : read_restart_file
  !
  use projection_mod, only : project_to_porder
  !
  use interpolation_mod, only : outerp
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,np,l,l1,l2,m
  integer :: gmin,gmax,omin,omax,ierr
  integer :: this_geom,this_order
  real(wp) :: rhoa2
  character(len=200) :: array_name
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: xyz_init(:,:)
  real(wp), allocatable :: uspstart(:,:)
  type(matrix), allocatable, target :: outerp_inv(:,:)
  !
  !.. Local Pointers ..
  real(wp), pointer, contiguous :: numer2init(:,:)
  real(wp), pointer, contiguous :: init2numer(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_flow"
  !
continue
  !
  ! Initialize all local pointers to disassociated
  !
  numer2init => null()
  init2numer => null()
  !
  ! Some constants from the free-stream conditions
  !
  rhoa2 = rhoref * aref * aref
  !
  ! Allocate the arrays within the module flowvar
  !
  call allocate_flowvar_arrays
  !
  ! Allocate the local arrays
  !
  allocate ( uspstart(1:nq,1:n_solpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"uspstart",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( xyz_init(1:nr,1:maxpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_init",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate the local array that stores the inverse of the matrix used for
  ! interpolating from the solution points to the initialization points.
  ! NOTE: This is array is used so that we only invert the matrix once for a
  !       given cell geometry/order combination.  We dont want to keep inverting
  !       the same matrix over and over for cells with the same geometry/order
  !       combination.
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
    ! Check whether the init component of outerp is associated with
    ! an interpolation matrix to determine where we need to compute
    ! the initial conditions.
    !
    if (.not. allocated(outerp(this_geom,this_order)%init)) then
      !
      ! The init component of outerp is not associated so just compute
      ! the initial conditions at the solution points
      !
      uspstart(1:nq,n1:n2) = initial_conditions( xyz(1:nr,n1:n2) )
      !
    else
      !
      ! The init component of outerp is associated so we need to compute
      ! the initial conditions at the initialization points before
      ! interpolating it back to the solution points.
      !
      ! Assign local pointer to this matrix as an alias to simplify the code
      !
      numer2init => outerp(this_geom,this_order)%init%mat
      !
      ! Check if we have already created an inverse of this interpolation
      ! matrix. If we havent, allocate its location in outerp_inv and
      ! compute this inverse matrix.
      !
      if (.not. allocated(outerp_inv(this_geom,this_order)%mat)) then
        !
        l1 = size(numer2init,dim=1)
        l2 = size(numer2init,dim=2)
        !
        allocate ( outerp_inv(this_geom,this_order)%mat(1:l1,1:l2) , &
                   source=zero , stat=ierr , errmsg = error_message )
        write (array_name,10) Geom_Name(this_geom),this_order
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message)
        !
        outerp_inv(this_geom,this_order)%mat = invert_matrix(numer2init)
        !
      end if
      !
      ! Assign local pointer to this matrix as an alias to simplify the code
      !
      init2numer => outerp_inv(this_geom,this_order)%mat
      !
      ! Interpolate the cell solution points to the initialization points
      !
      do l = 1,nr
        xyz_init(l,1:np) = matmul( numer2init , xyz(l,n1:n2) )
      end do
      !
      ! Evaluate the initial conditions at the initialization points
      !
      uspstart(1:nq,n1:n2) = initial_conditions( xyz_init(1:nr,1:np) )
      !
      ! Interpolate the initial conditions from the initialization points
      ! to the solution points
      !
      do m = 1,nq
        uspstart(m,n1:n2) = matmul( init2numer , uspstart(m,n1:n2) )
      end do
      !
      ! Disassociate all local pointers before continuing
      !
      if (associated(numer2init)) numer2init => null()
      if (associated(init2numer)) init2numer => null()
      !
    end if
    !
  end do
  !
  ! Check for a restart job, otherwise initialize usp to uspstart
  !
  if (irestart == 1) then
    !
    ! Read in usp and the current iteration, time step,
    ! and elapsed time from the restart file
    !
    call read_restart_file
    !
    ! Initialize the iteration counter, time step,
    ! and elapsed time to that in the restart file
    !
    d_t   = dtrst
    itcur = itrst
    time  = tmrst
    !
  else
    !
    ! Project the initial solution from the
    ! initialization points to the solution points
    !
   !call project_to_porder(cell,uspstart)
    !
    ! Initialize usp to uspstart
    !
    usp(:,1:n_solpts) = uspstart(:,1:n_solpts)
    !
    ! Initialize the iteration counter, time step, and elapsed time
    !
    d_t   = dt_in
    itcur = 0
    time  = zero
    !
  end if
  !
  ! Initialize the solution at previous timestep array
  !
  uoldsp(1:nq,1:n_solpts) = usp(1:nq,1:n_solpts)
  !
  ! Initialize the solution gradient array
  !
  call initialize_dusp
  !
  ! Output the free-stream conditions
  !
  if (mypnum == 0) then
    if (all(itestcase /= [Density_Transport,Entropy_Transport])) then
      !
      write (iout,1) rhoref,pref,tref,machref,aref,velref(1),velref(2), &
                     alpha_aoaref,beta_aoaref,gam,rgasref,cpref,muref,reyref
      write (iout,2) rhoref/rhoref,pref_nd,tref/tref,aref/aref, &
                     velref_nd(1),velref_nd(2),rgasref_nd,cpref_nd,muref_nd
      !
    else
      !
      write (iout,3) rhoref,pref,tref,machref,uref/machref, &
                     velref(1),velref(2),alpha_aoaref,beta_aoaref,gam, &
                     rgasref,muref,reyref
      !
    end if
    !
  end if
  !
  ! Deallocate local arrays and disassociate local pointers before leaving
  !
  if (allocated(uspstart)) then
    deallocate ( uspstart , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"uspstart",2,__LINE__,__FILE__,ierr,error_message)
  end if
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
  if (associated(numer2init)) numer2init => null()
  if (associated(init2numer)) init2numer => null()
  !
  ! Create the facexyz and edgexyz arrays if
  ! using the continuous output option or MMS
  !
  call get_face_and_edge_point_coordinates
  !
  ! Format Statements
  !
  1 format (/," Dimensional Free-stream conditions",/, &
              " Density          = ",es14.6," kg/m^3",/, &
              " Pressure         = ",es14.6," N/m^2",/, &
              " Temperature      = ",es14.6," K",/, &
              " Mach Number      = ",es14.6,/, &
              " Speed of Sound   = ",es14.6," m/s",/, &
              " X-Velocity       = ",es14.6," m/s",/, &
              " Y-Velocity       = ",es14.6," m/s",/, &
              " Alpha A.o.Attack = ",es14.6," degrees"/, &
              " Beta  A.o.Attack = ",es14.6," degrees"/, &
              " Gamma            = ",es14.6,/, &
              " Gas Constant     = ",es14.6," J/(kg*K)",/, &
              " Cp               = ",es14.6," J/(kg*K)",/, &
              " Viscosity        = ",es14.6," kg/(m*s)",/, &
              " Reynolds #       = ",es14.6,/)
  2 format (/," Non-dimensional Free-stream conditions",/, &
              " Density         = ",es14.6,/, &
              " Pressure        = ",es14.6,/, &
              " Temperature     = ",es14.6,/, &
              " Speed of Sound  = ",es14.6,/, &
              " X-Velocity      = ",es14.6,/, &
              " Y-Velocity      = ",es14.6,/, &
              " Gas Constant    = ",es14.6,/, &
              " Cp              = ",es14.6,/, &
              " Viscosity       = ",es14.6,/)
  3 format (/," The density vortex transport problem is being",/, &
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
  4 format (/," Computing the initial solution at the ",a," points.",/)
  5 format (/," ***************** ERROR! ****************",/, &
              "  THIS TEST CASE IS CURRENTLY NOT SET UP!",/, &
              " ***************** ERROR! ****************",/)
  10 format ("outerp_inv(",a,",",i0,")%mat")
  !
end subroutine initialize_flow
!
!###############################################################################
!
pure &
function initial_conditions(xyz) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : ncell
  use ovar, only : cv_in
  use ovar, only : itestcase
  use ovar, only : mms_opt,mms_init
  use ovar, only : postproc_debug_grid
  use mms_mod, only : mms_solution_at_xyz
  use channel_mod, only : channel_initialization
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: xyz
  !
  !.. Function Result ..
  real(wp), dimension(1:nq,1:size(xyz,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer  :: k,np,ierr
  real(wp) :: dl
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
    do k = 1,np
      return_value(:,k) = vortex_solution( xyz(:,k) )
    end do
    !
  else if (itestcase == Taylor_Green_Vortex) then
    !
    allocate ( adj_xyz , source=xyz , stat=ierr )
    !
    if (postproc_debug_grid) then
      dl = real(ncell,kind=wp)**(-one3)
      adj_xyz = PI * (two*dl*xyz - one)
    end if
    !
    do k = 1,np
      return_value(:,k) = taylor_green_initialization( adj_xyz(:,k) )
    end do
    !
  else if (itestcase == Channel_Flow) then
    !
    return_value(:,:) = channel_initialization( xyz(:,:) )
    !
  else if (mms_opt /= 0) then
    !
    do k = 1,np
      return_value(:,k) = mms_init * mms_solution_at_xyz( xyz(:,k) , &
                                                          return_nondim=true, &
                                                          return_cv=true )
    end do
    !
  else
    !
    do k = 1,np
      return_value(:,k) = cv_in(:)
    end do
    !
  end if
  !
end function initial_conditions
!
!###############################################################################
!
subroutine initialize_dusp()
  !
  ! This routine computes the discontinuous gradient at each solution point
  !
  !.. Use Statements ..
  use order_mod,  only : maxSP
  use geovar,  only : nr,ncell,cell
  use metrics_mod,  only : jacobian
  use metrics_mod,  only : metrics
  use derivatives_mod, only : deriv
  use flowvar, only : usp,dusp
  !
  !.. Local Scalars ..
  integer :: k,l,m,n,nc,n1,n2,np
  integer :: this_geom,this_order
  integer :: i,i1,i2,i3
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr) :: dudr
  real(wp), dimension(1:maxSP,1:nq) :: tusp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_dusp"
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
    n1 = cell(nc)%beg_sp ! beginning index for solpts in cell nc
    n2 = cell(nc)%end_sp ! ending    index for solpts in cell nc
    np = n2-n1+1 ! number of solution points in this cell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Store the transpose of usp for this cell to reduce cache misses
    !
    do k = 1,np
      do m = 1,nq
        tusp(k,m) = usp(m,k+n1-1)
      end do
    end do
    !
    ! Compute the gradient of the conservative variables at each solution point
    !
    do k = 1,np
      !
      n = k+n1-1 ! total index of solution point on current partition
      !
      i = this_order+1
      i3 = (k-1)/(i*i) + 1
      i2 = (k-1)/i - (i3-1)*i + 1
      i1 = k - (i2-1)*i - (i3-1)*i*i
      !
      do m = 1,nq
        !
        if (use_conservative_form) then
          !
          do l = 1,nr
            !
            dudr(1) = deriv(this_geom,this_order)%d(1)% &
              dot( k , tusp(1:np,m)*metrics(1,l,n1:n2,2)/jacobian(n1:n2) )
            !
            dudr(2) = deriv(this_geom,this_order)%d(2)% &
              dot( k , tusp(1:np,m)*metrics(2,l,n1:n2,2)/jacobian(n1:n2) )
            !
            dudr(3) = deriv(this_geom,this_order)%d(3)% &
              dot( k , tusp(1:np,m)*metrics(3,l,n1:n2,2)/jacobian(n1:n2) )
            !
            dusp(l,m,n) = sum(dudr)!* jacobian(n)
            !
          end do
          !
        else
          !
          ! Compute the derivatives in r-s space
          !
          do l = 1,nr
            dudr(l) = deriv(this_geom,this_order)%d(l)%dot( k , tusp(1:np,m) )
          end do
          !
          ! Transform the derivatives back to x-y space
          !
          do l = 1,nr
            dusp(l,m,n) = dot( dudr(1:nr) , metrics(1:nr,l,n,2) ) / jacobian(n)
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
end subroutine initialize_dusp
!
!###############################################################################
!
subroutine get_face_and_edge_point_coordinates
  !
  use order_mod,         only : maxSP
  use geovar,            only : nr,cell,face,xyz
  use flowvar,           only : edgexyz,facexyz
  use interpolation_mod, only : interp
  !
  !.. Local Scalars ..
  integer :: n,nc,ne,nf,np,n1,n2
  integer :: k,kf,l,lf,ierr
  integer :: face_order,edge_order
  integer :: this_order,this_geom
  !
  !.. Local Arrays ..
  real(wp) :: txyz(1:maxSP,1:nr)
  !
  !.. Local Allocatable Arrays ..
  logical(lk), allocatable, dimension(:) :: msk
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_face_and_edge_point_coordinates"
  !
continue
  !
  ! Only do this if facexyz is allocated which means that either
  ! the continuous output option or MMS is being used
  !
  if (allocated(facexyz)) then
    !
    do nf = 1,size(facexyz)
      !
      face_order = face(nf)%order
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
      do kf = 1,size(face(nf)%left%flxpts)
        !
        lf = face(nf)%left%flxpts(kf)
        !
        do l = 1,nr
          facexyz(nf)%v(l,kf) = interp(this_geom,this_order)% &
                                          toFace(face_order)% &
                                dot( lf , txyz(1:np,l) )
        end do
        !
      end do
      !
    end do
    !
  end if
  !
  ! Only do this if edgexyz is allocated which means
  ! that the continuous output option is being used
  !
  if (allocated(edgexyz)) then
    !
    ! Initialize edgexyz to zero
    !
    call edgexyz%zero
    !
    ! Allocate the mask array to mark which edges have been finished
    !
    allocate ( msk(1:size(edgexyz)) , source=fals , &
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
        do kf = 1,size(cell(nc)%edge(n)%cell_idx)
          !
          lf = cell(nc)%edge(n)%cell_idx(kf)
          !
          do l = 1,nr
            edgexyz(ne)%v(l,kf) = interp(this_geom,this_order)% &
                                            toEdge(edge_order)% &
                                  dot(lf,txyz(1:np,l))
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
end subroutine get_face_and_edge_point_coordinates
!
!###############################################################################
!
subroutine initialize_reference
  !
  !.. Use Statements ..
  use ovar,   only : rhoref,aref,tref,gam
  use ovar,   only : velref,velref_nd,uref
  use ovar,   only : pref,pref_nd,tmin,tmax
  use ovar,   only : rgasref,rgasref_nd
  use ovar,   only : cpref,cpref_nd,reyref
  use ovar,   only : cv_in,pv_in,machref
  use ovar,   only : ttotref,ttotref_nd
  use ovar,   only : ptotref,ptotref_nd
  use ovar,   only : ptot2p_ratio
  use ovar,   only : suth_Sref,suth_Tref
  use ovar,   only : alpha_aoaref,beta_aoaref
  use ovar,   only : itestcase,Lref
  !
  !.. Local Scalars ..
  integer  :: ierr
  real(wp) :: rhoa,rhoa2,tda2
  real(wp) :: gam_over_gm1
  real(wp) :: gm1_over_gam
  real(wp) :: total_cond
  real(wp) :: pres_check
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_reference"
  !
continue
  !
  ! Set up the base reference conditions
  !
  if (any(itestcase == Advec_Problems) .or. &
      abs(itestcase) == Shu_Vortex .or. &
      itestcase == -Vortex_Transport) then
    !
    ! Set up the base reference conditions for different
    ! versions of the isentropic Euler vortex problem
    !
    ! Set some reference conditions common to all of these vortex problems
    !
    beta_aoaref = zero
    velref(:) = zero
    !
    if (itestcase == -Shu_Vortex) then
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
    machref = one/sqrt(gam)
    aref = one
    uref = one
    rhoref = one
    pref = one
    tref = one
    rgasref = one
    !
  else if (itestcase == Taylor_Green_Vortex) then
    !
    ! If Lref is negative, the grid was not created on the fly and
    ! we need to add code to compute Lref and make sure the input
    ! grid is valid for the Taylor-Green vortex problem
    !
    if (Lref <= zero) then
      write (error_message,101)
      call stop_recon(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
    !
    ! Set up the base reference conditions for the Taylor-Green Vortex problem
    !
    alpha_aoaref = zero
    beta_aoaref = zero
    !
    gam = 1.4_wp
    tref = 300.0_wp ! (K)
    pref = 101325.0_wp ! (Pa)
    ! ( universal gas constant / molecular weight of air )
    rgasref = 8314.34_wp / 28.9651159_wp ! (J)/(kg*K)
    rhoref = pref / rgasref / tref ! (kg/m^3) ~ 1.17663794
    !
    aref = sqrt( gam*pref/rhoref ) ! (m/s) ~ 347.2169357
    machref = 0.1_wp
    uref = aref * machref ! (m/s) ~ 34.72169357
    !
    if (reyref <= zero) then
      reyref = 1600.0_wp
    end if
    !
    suth_Sref = 110.5556_wp
    suth_Tref = 273.111_wp
    !
  else
    !
    ! Set up the base reference conditions based
    ! on the input flow reference conditions
    !
    ! Check the reference conditions for consistency
    !
    if (rhoref < zero) then
      rhoref = pref/(rgasref*tref)
    else if (pref < zero) then
      pref = rhoref*rgasref*tref
    else if (tref < zero) then
      tref = pref/(rgasref*rhoref)
    else
      pres_check = rhoref*rgasref*tref
      if (abs(pres_check-pref) > eps6) then
       !write (iout,102) pref,pres_check
        write (error_message,102) pref,pres_check
        call stop_recon(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
    end if
    !
    ! Compute some of the remaining reference conditions
    !
    aref = sqrt( gam * pref / rhoref )
    uref = machref * aref
    velref(1) = uref * cos(beta_aoaref  * degrees2radians) * &
                       cos(alpha_aoaref * degrees2radians)
    velref(2) = uref * cos(beta_aoaref  * degrees2radians) * &
                       sin(alpha_aoaref * degrees2radians)
    velref(3) = uref * sin(beta_aoaref  * degrees2radians)
    !
  end if
  !
  ! Pre-compute some reference constants that will be used
  !
  rhoa = rhoref * aref
  rhoa2 = rhoa * aref
  tda2 = tref / aref / aref
  gam_over_gm1 = gam / (gam-one)
  gm1_over_gam = one / gam_over_gm1
  !
  ! Allocate the free stream inflow arrays
  !
  allocate ( cv_in(1:nq) , source=zero , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cv_in",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  allocate ( pv_in(1:nq) , source=zero , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"pv_in",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  !
  ! Initialize the inflow boundary conditions to the reference conditions
  ! NOTE: This is hard-wired for 2D and needs to be generalized for 3D
  !
  pv_in(nec) = rhoref
  pv_in(nmb:nme) = velref
  pv_in(nee) = pref
  !
  cv_in(nec) = pv_in(nec)
  cv_in(nmb:nme) = pv_in(nec) * pv_in(nmb:nme)
  cv_in(nee) = pv_in(nee)/(gam-one) + half*pv_in(nec)*selfdot(pv_in(nmb:nme))
  !
  ! Nondimensional reference velocities
  velref_nd = velref / aref
  !
  ! Nondimensional reference pressure
  pref_nd = pref / rhoa2
  !
  ! Nondimensional reference gas constant
  rgasref_nd = rgasref * tda2
  !
  ! Dimensional reference specific heat at constant pressure
  cpref    = rgasref * gam / (gam - one)
  ! Nondimensional reference specific heat at constant pressure
  cpref_nd = cpref * tda2
  !
  ! Finally, nondimensionalize the inflow boundary conditions
  !
  cv_in(nec)     = cv_in(nec) / rhoref
  cv_in(nmb:nme) = cv_in(nmb:nme) / rhoa
  cv_in(nee)     = cv_in(nee) / rhoa2
  !
  pv_in(nec)     = pv_in(nec) / rhoref
  pv_in(nmb:nme) = pv_in(nmb:nme) / aref
  pv_in(nee)     = pv_in(nee) / rhoa2
  !
  ! Total pressure and temperature
  !
  if (ptot2p_ratio < zero) then
    !
    total_cond = one + half*(gam-one)*machref*machref
    !
    ttotref = tref * total_cond
    ptotref = pref * (total_cond**gam_over_gm1)
    !
    ptot2p_ratio = ptotref/pref
    !
  else
    !
    ttotref = tref * (ptot2p_ratio**gm1_over_gam)
    ptotref = pref * ptot2p_ratio
    !
  end if
  !
  ttotref_nd = ttotref / tref
  ptotref_nd = ptotref / rhoa2
  !
  tmin = 10.0_wp
  tmax = 15000.0_wp
  !
  ! Initialize the viscosity reference conditions
  !
  call initialize_viscosity
  !
  ! Initialize the turbulence reference conditions
  !
  call initialize_turbulence
  !
  ! Format Statements
  !
  101 format ( &
     " The grid being used for the Taylor-Green vortex problem was not", &
     " generated on the fly and seems to be input from a grid file. Code", &
     " still needs to be added to compute Lref and make sure the input", &
     " grid is valid for the Taylor-Green vortex problem.")
  102 format ( &
     " Free-stream conditions in the input file are not consistent!",/, &
     " Reference pressure given    : ",es14.6,/, &
     " Reference pressure computed : ",es14.6,/, &
     " Please correct the input file or set one of density, pressure,", &
     " and temperature to a negative value so it is explicitly computed.")
  !
end subroutine initialize_reference
!
!###############################################################################
!
subroutine initialize_viscosity
  !
  !.. Use Statements ..
  use ovar, only : muref, muref_nd, reyref
  use ovar, only : pv_in, rhoref, aref, machref, tref
  use ovar, only : suth_muref, suth_c1
  use ovar, only : suth_Tref, suth_Sref
  !
  !.. Local Scalars ..
  real(wp) :: muref_scale,tref_scale
  !
continue
  !
  muref_scale = (one + suth_Sref/suth_Tref) / sqrt(suth_Tref)
  !
  if (reyref <= zero) then
    !
    suth_c1 = suth_muref * muref_scale
    !
    muref_nd = viscosity_pv_sp(pv_in)
    muref    = muref_nd * rhoref * aref
    !
    reyref = rhoref * machref * aref * reciprocal( muref )
    !
  else
    !
    muref = rhoref * machref * aref / reyref
    !
    muref_nd = muref / (rhoref * aref)
    !
    tref_scale = tref * sqrt(tref) / (tref + suth_Sref)
    !
    suth_muref = muref / muref_scale / tref_scale
    !
    suth_c1 = suth_muref * muref_scale
    !
  end if
  !
end subroutine initialize_viscosity
!
!###############################################################################
!
subroutine initialize_turbulence()
  !
  !.. Use Statements ..
  use ovar, only : tkeref, tkeref_nd
  use ovar, only : tlsref, tlsref_nd
  use ovar, only : tvrref, tinref
  use ovar, only : rhoref, muref
  use ovar, only : aref, machref
  use ovar, only : cv_in, pv_in
  use ovar, only : turbulence_model
  use ovar, only : turb_is_on
  use ovar, only : include_tke_in_energy
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_turbulence"
  !
continue
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
    cv_in(ntk) = cv_in(nec) * tkeref_nd
    cv_in(ntl) = cv_in(nec) * tlsref_nd
    !
    ! Add the tke contribution to the energy equation
    !
    cv_in(nee) = cv_in(nee) + cv_in(ntk)
    !
    ! Recompute the primitive variable inflow conditions
    !
    pv_in(1:nq) =  usp2v_sp( cv_in(1:nq) )
    !
  end if
  !
end subroutine initialize_turbulence
!
!###############################################################################
!
include "Functions/vortex.f90"
include "Functions/taylor_green.f90"
include "Functions/usp2v.f90"
include "Functions/viscosity.f90"
!
!###############################################################################
!
end module initialization_mod
