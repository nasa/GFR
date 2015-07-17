module initialization_mod
  !
  ! Module containing the functions and routines to
  ! initialize the solution
  !
  use module_kind_types
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  !
  implicit none
  !
  private
  !
  public :: initialize_flow
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
  use geovar, only : nr,n_solpts,ncell
  use geovar, only : pts_in_cell_ptr,xyz
  use geovar, only : cell_geom,cell_order
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
  use io_mod, only : read_restart_file
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
  character(len=200) :: error_message
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
  call alloc_error(pname,"uspstart",1,ierr,error_message)
  !
  allocate ( xyz_init(1:nr,1:maxpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_init",1,ierr,error_message)
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
  call alloc_error(pname,"outerp_inv",1,ierr,error_message)
  !
  ! Initialize the solution in each cell
  !
  do n = 1,ncell
    !
    this_geom = cell_geom(n)
    this_order = cell_order(n)
    !
    n1 = pts_in_cell_ptr(n)+1
    n2 = pts_in_cell_ptr(n+1)
    np = n2-n1+1
    !
    ! Check whether the init component of outerp is associated with
    ! an interpolation matrix to determine where we need to compute
    ! the initial conditions.
    !
    if (.not. associated(outerp(this_geom,this_order)%init)) then
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
      numer2init => outerp(this_geom,this_order)%init
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
        call alloc_error(pname,array_name,1,ierr,error_message)
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
   !call project_to_porder(ncell,cell_geom,cell_order,pts_in_cell_ptr,uspstart)
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
  ! Output the free-stream conditions
  !
  if (mypnum == 0) then
    if (itestcase /= Density_Transport .and. itestcase /= Entropy_Transport) then
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
    call alloc_error(pname,"uspstart",2,ierr,error_message)
  end if
  !
  if (allocated(outerp_inv)) then
    deallocate ( outerp_inv , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"outerp_inv",2,ierr,error_message)
  end if
  !
  if (allocated(xyz_init)) then
    deallocate ( xyz_init , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_init",2,ierr,error_message)
  end if
  !
  if (associated(numer2init)) numer2init => null()
  if (associated(init2numer)) init2numer => null()
  !
  ! If using MMS, we need to create the facexyz array
  !
  call get_boundary_fp_coordinates
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
pure function initial_conditions(xyz) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : cv_in
  use ovar, only : itestcase
  use ovar, only : mms_opt,mms_init
  use mms_mod, only : mms_solution_at_xyz
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: xyz
  !
  !.. Function Result ..
  real(wp), dimension(1:nq,1:size(xyz,dim=2)) :: return_value
  !
  !.. Local Scalars ..
  integer :: k,np
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
subroutine get_boundary_fp_coordinates
  !
  use order_mod, only : maxSP
  use geovar, only : nr,nfbnd,face,cell,xyz
  use flowvar, only : facexyz
  use interpolation_mod, only : interp
  !
  !.. Local Scalars ..
  integer :: nf,n1,n2,np
  integer :: kf,lf,l
  integer :: host_cell,host_geom
  integer :: host_order,face_order
  !
  !.. Local Arrays ..
  real(wp) :: lmat(1:maxSP)
  !
continue
  !
  ! Only do this if facexyz is allocated which means that we are using MMS
  !
  if (allocated(facexyz)) then
    !
    do nf = 1,nfbnd
      !
      face_order = face(nf)%order
      !
      host_cell = face(nf)%left%cell
      host_geom = cell(host_cell)%geom
      host_order = cell(host_cell)%order
      !
      n1 = cell(host_cell)%beg_sp
      n2 = cell(host_cell)%end_sp
      np = n2-n1+1
      !
      do kf = 1,size(face(nf)%left%flxpts)
        !
        lf = face(nf)%left%flxpts(kf)
        !
        lmat(1:np) = interp(host_geom,host_order)%toFace(face_order)%mat(1:np,lf)
        !
        do l = 1,nr
          facexyz(nf)%v(l,kf) = dot( lmat(1:np) , xyz(l,n1:n2) )
        end do
        !
      end do
      !
    end do
    !
  end if
  !
end subroutine get_boundary_fp_coordinates
!
!###############################################################################
!
include "Functions/vortex.f90"
!
!###############################################################################
!
end module initialization_mod
