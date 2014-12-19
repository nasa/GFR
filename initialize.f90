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
  use geovar,      only : nr,n_solpts,n_totpts,nface,ncell
  use geovar,      only : numer2gauss,numer2lobatto
  use geovar,      only : pts_in_cell_ptr,xyz
  !
 !use metvar,      only : 
  !
  use ovar,        only : rhoref,aref,iRKtype,num_stages,itestcase
  use ovar,        only : cv_in,init_loc
  use ovar,        only : irestart,d_t,dt_in,itcur,time
  use ovar,        only : dtrst,itrst,tmrst
  use ovar,        only : rhoref,pref,tref,gam,rgasref,cpref
  use ovar,        only : aref,machref,velref,uref,muref,reyref
  use ovar,        only : alpha_aoaref,beta_aoaref
  use ovar,        only : pref_nd,velref_nd,muref_nd
  use ovar,        only : rgasref_nd,cpref_nd
  !
  use flowvar,     only : usp,dusp,residual,dtsp,uoldsp
  use flowvar,     only : usp_sum,rssprk
  use flowvar,     only : faceusp,facedusp,faceflx
  !
  use order_mod,   only : n_fp1de,maxpts
  !
  use io_mod,      only : read_restart_file
  !
  use generic_mod, only : project_to_porder
  !
  !.. Local Scalars ..
  integer  :: n1,n2,n,np,k,l,m,ierr
  real(wp) :: rhoa2
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: xyz_init(:,:)
  real(wp), allocatable :: interp_mat(:,:)
  real(wp), allocatable, dimension(:,:) :: uspstart
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_flow"
  !
continue
  !
  ! Some constants from the free-stream conditions
  !
  rhoa2 = rhoref * aref * aref
  !
  ! Allocate the solution arrays
  !
  allocate ( usp(1:nq,1:n_totpts) , source=zero , stat=ierr )
  call alloc_error(pname,"usp",1,ierr)
  allocate ( dusp(1:nr,1:nq,1:n_totpts) , source=zero , stat=ierr )
  call alloc_error(pname,"dusp",1,ierr)
  allocate ( residual(1:nq,1:n_solpts) , source=zero , stat=ierr )
  call alloc_error(pname,"residual",1,ierr)
  allocate ( dtsp(1:n_solpts) , source=zero , stat=ierr )
  call alloc_error(pname,"dtsp",1,ierr)
  allocate ( uoldsp(1:nq,1:n_solpts) , source=zero , stat=ierr )
  call alloc_error(pname,"uoldsp",1,ierr)
  allocate ( faceusp(1:nq,1:n_fp1de,1:2,1:nface) , source=zero , stat=ierr )
  call alloc_error(pname,"faceusp",1,ierr)
  allocate ( facedusp(1:2,1:nq,1:n_fp1de,1:2,1:nface) , source=zero , stat=ierr )
  call alloc_error(pname,"facedusp",1,ierr)
  allocate ( faceflx(1:nq,1:n_fp1de,1:2,1:nface) , source=zero , stat=ierr )
  call alloc_error(pname,"faceflx",1,ierr)
  !
  ! Local allocatable array uspstart
  !
  allocate ( uspstart(1:nq,1:n_solpts) , source=zero , stat=ierr )
  call alloc_error(pname,"uspstart",1,ierr)
  !
  ! Allocate the storage for the Runge-Kutta methods
  !
  n = 0
  if (iRKtype == Classic_RK) then
    n = num_stages
  else if (iRKtype == CK4_RK .or. iRKtype == TVD4_RK) then
    n = 1
    if (iRKtype == TVD4_RK) then
      allocate ( usp_sum(1:nq,1:n_solpts) , source=zero , stat=ierr )
      call alloc_error(pname,"usp_sum",1,ierr)
    end if
  end if
  !
  if (n /= 0) then
    allocate ( rssprk(1:nq,1:n_solpts,1:n) , source=zero , stat=ierr )
    call alloc_error(pname,"rssprk",1,ierr)
  end if
  !
  ! Initialize uspstart depending on the test case requested
  !
  if (any(itestcase == [Generic_MeanFlow,Laminar_BndLyr])) then
    !
    ! Initialize the solution to free-stream conditions
    !
    do n = 1,n_solpts
      uspstart(:,n) = cv_in(:)
    end do
    !
  else if (any(itestcase == Transport_Problems)) then
    !
    ! Compute the initial solution for the vortex problem
    !
    if (any(init_loc == [1,2,-1])) then
      !
      allocate ( xyz_init(1:nr,1:maxpts) , source=zero , stat=ierr )
      call alloc_error(pname,"xyz_init",1,ierr)
      !
      if (abs(init_loc) == 1) then
        !
        write (*,4) "Gauss"
        n1 = size(numer2gauss,dim=1)
        n2 = size(numer2gauss,dim=2)
        allocate ( interp_mat(1:n1,1:n2) , source=numer2gauss , stat=ierr )
        call alloc_error(pname,"interp_mat",1,ierr)
        !
      else if (init_loc == 2) then
        !
        write (*,4) "Lobatto"
        n1 = size(numer2lobatto,dim=1)
        n2 = size(numer2lobatto,dim=2)
        allocate ( interp_mat(1:n1,1:n2) , source=numer2lobatto , stat=ierr )
        call alloc_error(pname,"interp_mat",1,ierr)
        !
      end if
      !
      do n = 1,ncell
        !
        n1 = pts_in_cell_ptr(n)+1
        n2 = pts_in_cell_ptr(n+1)
        np = n2-n1+1
        !
        do l = 1,nr
          xyz_init(l,1:np) = matmul( interp_mat , xyz(l,n1:n2) )
        end do
        !
        do k = 1,np
          uspstart(:,n1+k-1) = vortex2D( xyz_init(1,k) , xyz_init(2,k) )
        end do
        !
        do m = 1,nq
          uspstart(m,n1:n2) = matmul( interp_mat , uspstart(m,n1:n2) )
        end do
        !
      end do
      !
    else
      !
      do np = 1,n_solpts
        uspstart(:,np) = vortex2D( xyz(1,np) , xyz(2,np) )
      end do
      !
    end if
    !
  else
    !
    write (iout,5)
    call stop_recon(8,pname)
    !
  end if
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
    call project_to_porder(ncell,pts_in_cell_ptr,uspstart)
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
  ! Deallocate uspstart
  !
  deallocate ( uspstart , stat=ierr )
  call alloc_error(pname,"uspstart",2,ierr)
  if (allocated(interp_mat)) then
    deallocate ( interp_mat , stat=ierr )
    call alloc_error(pname,"interp_mat",2,ierr)
  end if
  if (allocated(xyz_init)) then
    deallocate ( xyz_init , stat=ierr )
    call alloc_error(pname,"xyz_init",2,ierr)
  end if
  !
  ! Create the profile structures
  !
 !if (ncpu == 1) call create_profile_structures
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
  !
end subroutine initialize_flow
!
!###############################################################################
!
!!subroutine create_profile_structures
!!  !
!!  !.. Use Statements ..
!!  use order_mod, only : n_fp1de,n_sp2dt,maxpts
!!  use geovar,    only : nfbnd,nbfai,nr,bface,xyz
!!  use geovar,    only : interp_s2f,pts_in_cell_ptr
!!  use flowvar,   only : bnd_prof
!!  !
!!  !.. Local Scalars ..
!!  integer :: host_cell,cell_face,ndata,ierr
!!  integer :: begfp,begbp,flxpt,bndpt,nqty
!!  integer :: min_grp,max_grp,nbcgrps,nfbnd_nc
!!  integer :: i,k,l,n,nf,np,wlocl,n1,n2,gid,ngf
!!  integer :: wall_qtys,flow_qtys
!!  character(len=100) :: array_name
!!  !
!!  !.. Local Arrays ..
!!  real(wp), dimension(1:maxpts) :: lmat
!!  !
!!  !.. Local Allocatable Arrays ..
!!  logical, allocatable, dimension(:) :: msk
!!  integer, allocatable, dimension(:) :: grp_ids
!!  integer, allocatable, dimension(:) :: bface_map
!!  integer, allocatable, dimension(:,:) :: bface_nc
!!  character(len=100), allocatable, dimension(:) :: flow_var
!!  !
!!  !.. Local Parameters ..
!!  character(len=*), parameter :: char_xyz(1:3) = ["X","Y","Z"]
!!  !
!!continue
!!  !
!!  ! Get the number of flow quantities and their names
!!  !
!!  flow_qtys = nr + nq + 2
!!  wall_qtys = flow_qtys + 6
!!  !
!!  allocate ( flow_var(1:wall_qtys) , stat=ierr )
!!  call alloc_error("create_profile_structures","flow_var",1,ierr)
!!  !
!!  do n = 1,nr
!!    flow_var(n) = "Coordinate"//char_xyz(n)
!!  end do
!!  flow_var(nr+1) = "Density"
!!  do n = 1,nr
!!    flow_var(nr+n+1) = "Velocity"//char_xyz(n)
!!  end do
!!  flow_var(nr+nq+0) = "Pressure"
!!  flow_var(nr+nq+1) = "Temperature"
!!  flow_var(nr+nq+2) = "Viscosity"
!!  flow_var(flow_qtys+1) = "SkinFrictionMagnitude"
!!  flow_var(flow_qtys+2) = "ReynoldsX"
!!  flow_var(flow_qtys+3) = "YPlus"
!!  flow_var(flow_qtys+4) = "WallStress"
!!  flow_var(flow_qtys+5) = "VelocityFriction"
!!  flow_var(flow_qtys+6) = "WallHeatFlux"
!!  !
!!  ! Create the non-communication boundary face msk array
!!  !
!!  allocate ( msk(1:nfbnd) , source=fals , stat=ierr )
!!  call alloc_error("create_profile_structures","msk",1,ierr)
!!  !
!!  msk(:) = (bface(1,:)/=bc_periodic .and. bface(1,:)/=bc_cpu_bnd)
!!  !
!!  ! Count the number of non-communication boundary faces
!!  !
!!  nfbnd_nc = count( msk )
!!  !
!!  ! Create the bface array for only non-communication type boundaries
!!  !
!!  allocate ( bface_nc(1:nbfai,1:nfbnd_nc) , source=0 , stat=ierr )
!!  call alloc_error("create_profile_structures","bface_nc",1,ierr)
!!  do n = 1,nbfai
!!    bface_nc(n,:) = pack( bface(n,:) , mask = msk(:) )
!!  end do
!!  !
!!  ! Create the non-communication boundary face mapping array. We need this
!!  ! to give us the original index for the boundary face (i.e., from 1:nfbnd)
!!  ! when working on the bface_nc array since the face index within bface_nc
!!  ! is not necessarily the same as its corresponding index in bface.
!!  !
!!  allocate ( bface_map(1:nfbnd_nc) , source=0 , stat=ierr )
!!  call alloc_error("create_profile_structures","bface_map",1,ierr)
!!  bface_map = pack( [1:nfbnd] , mask = msk )
!!  !
!!  ! Find the number of different boundary condition groups on this processor.
!!  !
!!  ! Get minimum and maximum group IDs in bface(5,:)
!!  ! for the non-communication boundaries.
!!  !
!!  min_grp = minval( bface_nc(5,:) )
!!  max_grp = maxval( bface_nc(5,:) )
!!  !
!!  ! Loop through the possible group IDs to see how many BC groups are
!!  ! present on this processor. At the same time, store any group IDs
!!  ! that are present.
!!  !
!!  nbcgrps = 0
!!  do n = min_grp,max_grp
!!    if (any(bface(5,:) == n)) then
!!      nbcgrps = nbcgrps + 1
!!      call reallocate( grp_ids , nbcgrps )
!!      grp_ids(nbcgrps) = n
!!    end if
!!  end do
!!  !
!!  ! Allocate the top level of the bnd_prof array now that we
!!  ! know how many BC groups are on this processor
!!  !
!!  allocate ( bnd_prof(1:nbcgrps) , stat=ierr )
!!  call alloc_error("create_profile_structures","bnd_prof",1,ierr)
!!  !
!!  ! Loop back through the boundary groups and get the indices
!!  ! of the boundary faces belonging to each group ID.
!!  !
!!  do n = 1,size(grp_ids)
!!    !
!!    ! Store this group ID in the bnd_prof array
!!    !
!!    gid = grp_ids(n)
!!    bnd_prof(n)%group_num = gid
!!    !
!!    ! Count the number of boundary faces with this group ID. Also get
!!    ! the number of data points that will be needed for this group ID.
!!    ! #########################################################
!!    ! #########################################################
!!    ! #####################  WARNING!!!!  #####################
!!    ! #########################################################
!!    ! #########################################################
!!    ! ###  THE VALUE ndata IS ONLY VALID FOR 2D WHERE THE   ###
!!    ! ###  NUMBER OF FACE POINTS IS THE NUMBER OF EDGE      ###
!!    ! ###  POINTS REGARDLESS OF THE HOST CELL GEOMETRY.     ###
!!    ! ###  FOR 3D, WE WILL HAVE TO WORRY ABOUT DIFFERING    ###
!!    ! ###  NUMBERS OF FACE POINTS DUE TO TRIANGLE AND       ###
!!    ! ###  QUADRILATERAL BOUNDARY FACES. WE WILL ALSO HAVE  ###
!!    ! ###  TO WORRY ABOUT PROVIDING THE MESH CONNECTIVITY   ###
!!    ! ###  OF THE BOUNDARY FACES FOR 3D.                    ###
!!    ! #########################################################
!!    ! #########################################################
!!    ! #####################  WARNING!!!!  #####################
!!    ! #########################################################
!!    ! #########################################################
!!    !
!!    !
!!    ngf = count(bface_nc(5,:) == gid)
!!    ndata = ngf * n_fp1de
!!    !
!!    ! Allocate the faces part of the bnd_prof array for this group ID
!!    !
!!    allocate ( bnd_prof(n)%faces(1:ngf) , source=0 , stat=ierr )
!!    write (array_name,1) "bnd_prof(",n,")%faces"
!!    call alloc_error("create_profile_structures",array_name,1,ierr)
!!    !
!!    ! Store the indices for the boundary faces that have this group ID
!!    ! NOTE: We need to use bface_map as the packing array to map the face
!!    !       index within bface_nc back to the index for bface.
!!    !
!!    bnd_prof(n)%faces = pack( bface_map , mask = (bface_nc(5,:) == gid) )
!!    !
!!    ! Get the name of the BC for this group ID using the first face in the group.
!!    !
!!    nf = bnd_prof(n)%faces(1)
!!    bnd_prof(n)%prof_name = trim( bc_integer_to_string(bface(1,nf)) )
!!    !
!!    ! Allocate the number of flow quantities
!!    !
!!    nqty = wall_qtys
!!   !nqty = merge( wall_qtys , flow_qtys , any(bface(1,nf) == bc_walls) )
!!    !
!!    allocate ( bnd_prof(n)%qty(1:nqty) , stat=ierr )
!!    write (array_name,1) "bnd_prof(",n,")%qty"
!!    call alloc_error("create_profile_structures",array_name,1,ierr)
!!    !
!!    ! Loop over the number of flow quantities and
!!    ! allocate the data locations for each quantity.
!!    !
!!    do i = 1,nqty
!!      !
!!      ! Store the name of the flow quantity
!!      !
!!      bnd_prof(n)%qty(i)%qty_name = flow_var(i)
!!      !
!!      ! Allocate the arrays containing the actual data for each flow quantity
!!      !
!!      allocate ( bnd_prof(n)%qty(i)%qty_data(1:ndata) , &
!!                 source=zero , stat=ierr )
!!      write (array_name,1) "bnd_prof(",n,")%qty(",i,")%qty_data"
!!      call alloc_error("create_profile_structures",array_name,1,ierr)
!!      !
!!    end do
!!    !
!!  end do
!!  !
!!  ! Fill the coordinate quantity sections of boundary profile
!!  !
!!  do n = 1,size(bnd_prof)
!!    !
!!    do i = 1,size(bnd_prof(n)%faces)
!!      !
!!      nf = bnd_prof(n)%faces(i) ! boundary face of current profile
!!      !
!!      host_cell = bface(2,nf) ! host cell on current boundary face
!!      cell_face = bface(4,nf) ! face of host cell corresponding to boundary face nf
!!      !
!!      n1 = pts_in_cell_ptr(host_cell)+1 ! beginning index for sol-pts in host cell
!!      n2 = pts_in_cell_ptr(host_cell+1) ! ending    index for sol-pts in host cell
!!      np = n2-n1+1                      ! number of sol-pts in host cell
!!      wlocl = merge(1,2,np == n_sp2dt)     ! check if host cell is tria or quad
!!      !
!!      begfp = (cell_face-1)*n_fp1de ! Beginning flx-pt index on host cell face
!!      begbp = (i-1)*n_fp1de         ! Beginning flx-pt index within qty_data
!!      !
!!      do k = 1,n_fp1de
!!        !
!!        flxpt = begfp + k ! Index of flx-pt local to the current host cell
!!        bndpt = begbp + k ! Index within qty_data for current flx-pt
!!        !
!!        lmat(1:np) = interp_s2f(1:np,flxpt,wlocl)
!!        !
!!        ! Interpolate the coordinates of the interior sol-pts to the
!!        ! current flux point on face nf
!!        !
!!        do l = 1,nr
!!          bnd_prof(n)%qty(l)%qty_data(bndpt) = dot( lmat(1:np) , xyz(l,n1:n2) )
!!        end do
!!        !
!!      end do
!!      !
!!    end do
!!    !
!!  end do
!! !!
!! !! Output the interpolated coordinates to Tecplot for debugging
!! !!
!! !open (iogrd,file="bndprof.tec.dat",status="replace",action="write", &
!! !            form="formatted",position="rewind",iostat=ierr)
!! !call io_error("create_profile_structures","bndprof.tec.dat",1,ierr)
!! !!
!! !write (iogrd,100) (trim(bnd_prof(1)%qty(l)%qty_name), l=1,nr)
!! !do n = 1,size(bnd_prof)
!! !  !
!! !  np = size(bnd_prof(n)%qty(1)%qty_data)
!! !  !
!! !  write (iogrd,101) trim(bnd_prof(n)%prof_name),np
!! !  !
!! !  do l = 1,nr
!! !    write (iogrd,102) (bnd_prof(n)%qty(l)%qty_data(i), i=1,np)
!! !  end do
!! !  !
!! !end do
!! !!
!! !close (iogrd,iostat=ierr)
!! !call io_error("create_profile_structures","bndprof.tec.dat",2,ierr)
!!  !
!!  ! Format Statements
!!  !
!!  1 format (a,*(i0,a,:))
!!  2 format (2x,i0,5x,"'",a,"'")
!!  3 format ("BC Group# = ",i0,": ",a)
!!  100 format ('VARIABLES = "',a,'" "',a,'"')
!!  101 format ('ZONE T="',a,'", I=',i0,', DATAPACKING=BLOCK')
!!  102 format (5es18.10)
!!  !
!!end subroutine create_profile_structures
!
!###############################################################################
!
include "Functions/vortex.f90"
include "Functions/usp2v.f90"
!
!###############################################################################
!
end module initialization_mod
