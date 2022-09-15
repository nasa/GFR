module internal_grid_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  private
  !
  public :: create_internal_grid
  !
contains
!
!###############################################################################
!
subroutine create_internal_grid()
  !
  !.. Use Statements ..
  use ovar, only : itestcase,postproc_debug_grid
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_internal_grid"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (postproc_debug_grid) then
    !
    call create_postproc_debug_grid
    !
  else if (itestcase == Taylor_Green_Vortex) then
    !
    call create_TaylorGreenVortex_grid
    !
  else if (any(abs(itestcase) == Transport_Problems)) then
    !
    call create_IsentropicVortex_grid
    !
  else if (abs(itestcase) == CH4_Flame_Speed_Test) then
    !
    if (flame_speed_test_dof <= 0) then
      flame_speed_test_dof = 512
    end if
    !
    call create_flame_speed_test_grid
    !
  else if (abs(itestcase) == Channel_Flow) then
    !
    call create_channel_flow_grid
    !
  else
    !
    if (mypnum == glb_root) write (iout,1)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/,"************************* ERROR! *************************",/, &
              "  THE GRID TYPE SEEMS TO INDICATE CREATING AN INTERNALLY",/, &
              "  GENERATED GRID, BUT THE TEST CASE REQUESTED DOES NOT",/, &
              "  SUPPORT INTERNALLY GENERATING A GRID!",/, &
              "************************* ERROR! *************************",/)
  !
end subroutine create_internal_grid
!
!###############################################################################
!
subroutine create_IsentropicVortex_grid()
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use ovar,      only : Lref,itestcase
  use ovar,      only : VortexGrid_DOF,VortexGrid_Lref
  use ovar,      only : VortexGrid_random_perturb
  use ovar,      only : VortexGrid_perturb_percent
  use geovar,    only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geovar,    only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geovar,    only : cell_geom,cell_order
  use ovar,      only : grid_format
  !
  !.. Local Scalars ..
  integer :: nc,np,ierr,i,j,l,n
  integer :: nf_min,nf_max
  real(r4) :: x,y,dl
  real(r4) :: rbas,radj,rval
  !
  !.. Local Allocatable Arrays ..
  real(r4), allocatable :: perturb(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_IsentropicVortex_grid"
 !logical(lk), parameter :: test_random_bc = true
  logical(lk), parameter :: test_random_bc = fals
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Reset the grid format input parameter to identify
  ! the grid as internally generated
  !
  grid_format = Internally_Generated_Grid
  !
  ! If we are generating the grid on the fly, make the reference length 1
  !
  Lref = one
  !
  ! The isentropic Euler vortex problem is always 2D
  !
  nr = 2
  !
  ! Estimate the number of cells required to get the requested DOF
  !
  nc = nint( real(VortexGrid_DOF,kind=wp) / real(n_order+1,kind=wp) )
  !
  np = nc + 1 ! number of grid nodes in each direction
  !
  ncell = nc*nc ! total number of grid cells
  nnode = np*np ! total number of grid nodes
  nfbnd = 4*nc  ! total number of boundary faces
  !
  ! Allocate all the arrays needed to define the grid
  !
  allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_nodes",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  allocate ( nodes_of_cell(1:geom_nodes(Geom_Quad)*ncell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_geom(1:ncell) , source=Geom_Quad , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_geom",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_order(1:ncell) , source=n_order , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_order",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( perturb(1:nr,1:nnode) , source=real(0,kind=r4) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"perturb",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! #######
  ! #######
  ! ## 1 ##
  ! #######
  ! #######
  !
  ! Create the nodes_of_cell_ptr array
  !
  nodes_of_cell_ptr(1) = 0
  do n = 2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(n) = nodes_of_cell_ptr(n-1) + geom_nodes(Geom_Quad)
  end do
  !
  ! Create the nodes_of_cell array
  !
  do j = 1,nc
    do i = 1,nc
      !
      n = nodes_of_cell_ptr( map2dto1d(i,j,nc) )
      !
      nodes_of_cell(n+1) = map2dto1d(i  ,j  ,np)
      nodes_of_cell(n+2) = map2dto1d(i+1,j  ,np)
      nodes_of_cell(n+3) = map2dto1d(i+1,j+1,np)
      nodes_of_cell(n+4) = map2dto1d(i  ,j+1,np)
      !
    end do
  end do
  !
  ! #######
  ! #######
  ! ## 2 ##
  ! #######
  ! #######
  !
  ! Compute the cell length for a Cartesian/uniform grid
  !
  dl = real(1,kind=r4) / real(np-1,kind=r4)
  !
  ! Compute the random perturbations if the flag was set on input
  !
  if (VortexGrid_random_perturb) then
    !
    rbas = dl * real(abs(VortexGrid_perturb_percent),kind=r4) &
              / real(100,kind=r4)
    radj = rbas / real(2,kind=r4)
    !
    ! Randomly perturb the crap out of everything
    !
    do j = 2,np-1
      do i = 2,np-1
        !
        n = map2dto1d(i,j,np)
        !
        do l = 1,nr
          rval = random(rbas) - radj
          perturb(l,n) = rval
        end do
        !
      end do
    end do
    !
    ! Have the root processor broadcast its perturbations to all the
    ! other processors so that the perturbations will be consistent
    ! on all processors.
    !
    call mpi_bcast(perturb,size(perturb,kind=int_mpi), &
                   MPI_REAL,0_int_mpi,MPI_COMM_WORLD,mpierr)
    !
  end if
  !
  ! Create the grid node coordinates as a square from (0,0) to (1,1)
  !
  do j = 1,np
    !
    y = dl * real(j-1,kind=r4)
    !
    do i = 1,np
      !
      x = dl * real(i-1,kind=r4)
      !
      n = map2dto1d(i,j,np)
      !
      xyz_nodes(1,n) = real( x + perturb(1,n) , kind=wp )
      xyz_nodes(2,n) = real( y + perturb(2,n) , kind=wp )
      !
    end do
  end do
  !
  ! Deallocate perturb since it is no longer needed
  !
  deallocate ( perturb , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"perturb",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Scale the grid node coordinates to the square (-1,-1) to (1,1)
  !
  xyz_nodes = two*xyz_nodes - one
  !
  ! Scale the grid node coordinates to the square (-Lref,-Lref) to (+Lref,+Lref)
  !
  xyz_nodes = abs(VortexGrid_Lref) * xyz_nodes
  !
  ! #######
  ! #######
  ! ## 3 ##
  ! #######
  ! #######
  !
  ! Create the bface array
  !
  only_root_creates_bface: if (mypnum == glb_root) then
    !
    nf_max = 0
    !
    ! Loop through the 4 grid boundaries
    ! NOTE: This corresponds to the grid boundaries in
    !       the order [imin,imax,jmin,jmax]
    !
    nf_min = nf_max
    !
    do j = 1,nc
      !
      nf_min = nf_min + 1
      nf_max = nf_min + nc
      !
      bface( 1,nf_min) = bc_periodic
      bface( 2,nf_min) = map2dto1d(1 ,j,nc)
      bface( 3,nf_min) = Geom_Quad
      bface( 4,nf_min) = HostFace2D(1)
      bface( 5,nf_min) = nf_max
      bface( 6,nf_min) = map2dto1d(nc,j,nc)
      bface( 8,nf_min) = Geom_Quad
      bface( 9,nf_min) = HostFace2D(2)
      bface(10,nf_min) = ncell + nf_min
      !
      bface( 1,nf_max) = bc_periodic
      bface( 2,nf_max) = bface(6,nf_min)
      bface( 3,nf_max) = bface(8,nf_min)
      bface( 4,nf_max) = bface(9,nf_min)
      bface( 5,nf_max) = nf_min
      bface( 6,nf_max) = bface(2,nf_min)
      bface( 8,nf_max) = bface(3,nf_min)
      bface( 9,nf_max) = bface(4,nf_min)
      bface(10,nf_max) = ncell + nf_max
      !
      ! Adjust boundary conditions for stationary vortex or
      ! vortex propagation in only the vertical direction
      !
      if (any(itestcase == [-Shu_Vortex,-Vortex_Transport])) then
        !
        if (test_random_bc) then
          l = size(bc_heirarchy)
          l = max( 1 , min(l,random(l)) )
          bface(1,nf_min) = bc_heirarchy(l)
        else
          bface(1,nf_min) = bc_characteristic
        end if
        bface(  5,nf_min) = 1
        bface(6:9,nf_min) = 0
        !
        if (test_random_bc) then
          l = size(bc_heirarchy)
          l = max( 1 , min(l,random(l)) )
          bface(1,nf_max) = bc_heirarchy(l)
        else
          bface(1,nf_max) = bc_characteristic
        end if
        bface(  5,nf_max) = 1
        bface(6:9,nf_max) = 0
        !
      end if
      !
    end do
    !
    nf_min = nf_max
    !
    do i = 1,nc
      !
      nf_min = nf_min + 1
      nf_max = nf_min + nc
      !
      bface( 1,nf_min) = bc_periodic
      bface( 2,nf_min) = map2dto1d(i,1 ,nc)
      bface( 3,nf_min) = Geom_Quad
      bface( 4,nf_min) = HostFace2D(3)
      bface( 5,nf_min) = nf_max
      bface( 6,nf_min) = map2dto1d(i,nc,nc)
      bface( 8,nf_min) = Geom_Quad
      bface( 9,nf_min) = HostFace2D(4)
      bface(10,nf_min) = ncell + nf_min
      !
      bface( 1,nf_max) = bc_periodic
      bface( 2,nf_max) = bface(6,nf_min)
      bface( 3,nf_max) = bface(8,nf_min)
      bface( 4,nf_max) = bface(9,nf_min)
      bface( 5,nf_max) = nf_min
      bface( 6,nf_max) = bface(2,nf_min)
      bface( 8,nf_max) = bface(3,nf_min)
      bface( 9,nf_max) = bface(4,nf_min)
      bface(10,nf_max) = ncell + nf_max
      !
      ! Adjust boundary conditions for stationary vortex
      !
      if (itestcase == -Shu_Vortex) then
        !
        if (test_random_bc) then
          l = size(bc_heirarchy)
          l = max( 1 , min(l,random(l)) )
          bface(1,nf_min) = bc_heirarchy(l)
        else
          bface(1,nf_min) = bc_characteristic
        end if
        bface(  5,nf_min) = 1
        bface(6:9,nf_min) = 0
        !
        if (test_random_bc) then
          l = size(bc_heirarchy)
          l = max( 1 , min(l,random(l)) )
          bface(1,nf_max) = bc_heirarchy(l)
        else
          bface(1,nf_max) = bc_characteristic
        end if
        bface(  5,nf_max) = 1
        bface(6:9,nf_max) = 0
        !
      end if
      !
    end do
    !
  end if only_root_creates_bface
  !
  if (ncpu > 1) then
    call mpi_bcast(bface,size(bface,kind=int_mpi), &
                   mpi_inttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  contains
  !
  pure function map2dto1d(i,j,n) result(return_value)
    integer, intent(in) :: i,j,n
    integer :: return_value
    return_value = i + n*(j-1)
  end function map2dto1d
  !
end subroutine create_IsentropicVortex_grid
!
!###############################################################################
!
subroutine create_TaylorGreenVortex_grid()
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use ovar,      only : VortexGrid_DOF,Lref
  use ovar,      only : VortexGrid_random_perturb
  use ovar,      only : VortexGrid_perturb_percent
  use geovar,    only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geovar,    only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geovar,    only : cell_geom,cell_order
  use ovar,      only : grid_format
  !
  !.. Local Scalars ..
  integer :: nc,np,ierr,i,j,k,l,n,n1,n2
  integer :: nf_min,nf_max
  real(r4) :: x,y,z,dl
  real(r4) :: rbas,radj,rval
  !
  !.. Local Allocatable Arrays ..
  real(r4), allocatable :: perturb(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_TaylorGreenVortex_grid"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Reset the grid format input parameter to identify
  ! the grid as internally generated
  !
  grid_format = Internally_Generated_Grid
  !
  ! If we are generating the grid on the fly, make the reference length 1
  !
  Lref = one
  !
  ! The Taylor-Green vortex problem is always 3D
  !
  nr = 3
  !
  ! Estimate the number of cells required to get the requested DOF
  !
  nc = nint( real(VortexGrid_DOF,kind=wp) / real(n_order+1,kind=wp) )
  !
  np = nc + 1 ! number of grid nodes in each direction
  !
  ncell = nc*nc*nc ! total number of grid cells
  nnode = np*np*np ! total number of grid nodes
  nfbnd = 6*nc*nc  ! total number of boundary faces
  !
  ! Allocate all the arrays needed to define the grid
  !
  allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_nodes",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  allocate ( nodes_of_cell(1:geom_nodes(Geom_Hexa)*ncell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_geom(1:ncell) , source=Geom_Hexa , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_geom",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_order(1:ncell) , source=n_order , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_order",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( perturb(1:nr,1:nnode) , source=real(0,kind=r4) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"perturb",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! #######
  ! #######
  ! ## 1 ##
  ! #######
  ! #######
  !
  ! Create the nodes_of_cell_ptr array
  !
  nodes_of_cell_ptr(1) = 0
  do n = 2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(n) = nodes_of_cell_ptr(n-1) + geom_nodes(Geom_Hexa)
  end do
  !
  ! Create the nodes_of_cell array
  !
  do k = 1,nc
    do j = 1,nc
      do i = 1,nc
        !
        n = nodes_of_cell_ptr( map3dto1d(i,j,k,nc) )
        !
        nodes_of_cell(n+1) = map3dto1d(i  ,j  ,k  ,np)
        nodes_of_cell(n+2) = map3dto1d(i+1,j  ,k  ,np)
        nodes_of_cell(n+3) = map3dto1d(i+1,j+1,k  ,np)
        nodes_of_cell(n+4) = map3dto1d(i  ,j+1,k  ,np)
        nodes_of_cell(n+5) = map3dto1d(i  ,j  ,k+1,np)
        nodes_of_cell(n+6) = map3dto1d(i+1,j  ,k+1,np)
        nodes_of_cell(n+7) = map3dto1d(i+1,j+1,k+1,np)
        nodes_of_cell(n+8) = map3dto1d(i  ,j+1,k+1,np)
        !
      end do
    end do
  end do
  !
  ! #######
  ! #######
  ! ## 2 ##
  ! #######
  ! #######
  !
  ! Compute the cell length for a Cartesian/uniform grid
  !
  dl = real(1,kind=r4) / real(np-1,kind=r4)
  !
  ! Compute the random perturbations if the flag was set on input
  !
  if (VortexGrid_random_perturb) then
    !
    rbas = dl * real(abs(VortexGrid_perturb_percent),kind=r4) &
              / real(100,kind=r4)
    radj = rbas / real(2,kind=r4)
    !
    if (VortexGrid_perturb_percent < 0) then
      !
      ! Randomly perturb the crap out of everything
      !
      do k = 2,np-1
        do j = 2,np-1
          do i = 2,np-1
            !
            n = map3dto1d(i,j,k,np)
            !
            do l = 1,nr
              rval = random(rbas) - radj
              perturb(l,n) = rval
            end do
            !
          end do
        end do
      end do
      !
    else
      !
      ! Get the random perturbation in along i for y
      !
      do i = 2,np-1
        rval = random(rbas) - radj
        do k = 2,np-1
          do j = 2,np-1
            n = map3dto1d(i,j,k,np)
            perturb(2,n) = rval
          end do
        end do
      end do
      !
      ! Get the random perturbation along j for z
      !
      do j = 2,np-1
        rval = random(rbas) - radj
        do k = 2,np-1
          do i = 2,np-1
            n = map3dto1d(i,j,k,np)
            perturb(3,n) = rval
          end do
        end do
      end do
      !
      ! Get the random perturbation along k for x
      !
      do k = 2,np-1
        rval = random(rbas) - radj
        do j = 2,np-1
          do i = 2,np-1
            n = map3dto1d(i,j,k,np)
            perturb(1,n) = rval
          end do
        end do
      end do
      !
    end if
    !
    ! Have the root processor broadcast its perturbations to all the
    ! other processors so that the perturbations will be consistent
    ! on all processors.
    !
    call mpi_bcast(perturb,size(perturb,kind=int_mpi), &
                   MPI_REAL,0_int_mpi,MPI_COMM_WORLD,mpierr)
    !
  end if
  !
  ! Create the grid node coordinates as a cube from (0,0,0) to (1,1,1)
  !
  do k = 1,np
    !
    z = dl * real(k-1,kind=r4)
    !
    do j = 1,np
      !
      y = dl * real(j-1,kind=r4)
      !
      do i = 1,np
        !
        x = dl * real(i-1,kind=r4)
        !
        n = map3dto1d(i,j,k,np)
        !
        xyz_nodes(1,n) = real( x + perturb(1,n) , kind=wp )
        xyz_nodes(2,n) = real( y + perturb(2,n) , kind=wp )
        xyz_nodes(3,n) = real( z + perturb(3,n) , kind=wp )
        !
      end do
    end do
  end do
  !
  ! Deallocate perturb since it is no longer needed
  !
  deallocate ( perturb , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"perturb",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Scale the grid node coordinates to the cube (-1,-1,-1) to (1,1,1)
  !
  xyz_nodes(:,:) = two*xyz_nodes(:,:) - one
  !
  ! Scale the grid node coordinates to the cube (-pi,-pi,-pi) to (pi,pi,pi)
  !
  xyz_nodes(:,:) = PI * xyz_nodes(:,:)
  !
  ! #######
  ! #######
  ! ## 3 ##
  ! #######
  ! #######
  !
  ! Create the bface array
  !
  nf_max = 0
  !
  ! Loop through the 6 grid boundaries
  ! NOTE: This corresponds to the grid boundaries in
  !       the order [imin,imax,jmin,jmax,kmin,kmax]
  !
  ! Min-X and Max-X boundary faces
  !
  n1 = 0          ! min-x counter
  n2 = n1 + nc*nc ! max-x counter
  do k = 1,nc
    do j = 1,nc
      !
      n1 = n1 + 1
      n2 = n2 + 1
      !
      bface( 1,n1) = bc_periodic
      bface( 2,n1) = map3dto1d(1 ,j,k,nc)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(1)
      bface( 5,n1) = n2
      bface( 6,n1) = map3dto1d(nc,j,k,nc)
      bface( 8,n1) = Geom_Hexa
      bface( 9,n1) = HostFace3D(2)
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_periodic
      bface( 2,n2) = bface(6,n1)
      bface( 3,n2) = bface(8,n1)
      bface( 4,n2) = bface(9,n1)
      bface( 5,n2) = n1
      bface( 6,n2) = bface(2,n1)
      bface( 8,n2) = bface(3,n1)
      bface( 9,n2) = bface(4,n1)
      bface(10,n2) = ncell + n2
      !
    end do
  end do
  !
  ! Min-Y and Max-Y boundary faces
  !
  n1 = 2*nc*nc    ! min-y counter
  n2 = n1 + nc*nc ! max-y counter
  do k = 1,nc
    do i = 1,nc
      !
      n1 = n1 + 1
      n2 = n2 + 1
      !
      bface( 1,n1) = bc_periodic
      bface( 2,n1) = map3dto1d(i,1 ,k,nc)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(3)
      bface( 5,n1) = n2
      bface( 6,n1) = map3dto1d(i,nc,k,nc)
      bface( 8,n1) = Geom_Hexa
      bface( 9,n1) = HostFace3D(4)
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_periodic
      bface( 2,n2) = bface(6,n1)
      bface( 3,n2) = bface(8,n1)
      bface( 4,n2) = bface(9,n1)
      bface( 5,n2) = n1
      bface( 6,n2) = bface(2,n1)
      bface( 8,n2) = bface(3,n1)
      bface( 9,n2) = bface(4,n1)
      bface(10,n2) = ncell + n2
      !
    end do
  end do
  !
  ! Min-Z and Max-Z boundary faces
  !
  n1 = 4*nc*nc    ! min-z counter
  n2 = n1 + nc*nc ! max-z counter
  do j = 1,nc
    do i = 1,nc
      !
      n1 = n1 + 1
      n2 = n2 + 1
      !
      bface( 1,n1) = bc_periodic
      bface( 2,n1) = map3dto1d(i,j,1 ,nc)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(5)
      bface( 5,n1) = n2
      bface( 6,n1) = map3dto1d(i,j,nc,nc)
      bface( 8,n1) = Geom_Hexa
      bface( 9,n1) = HostFace3D(6)
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_periodic
      bface( 2,n2) = bface(6,n1)
      bface( 3,n2) = bface(8,n1)
      bface( 4,n2) = bface(9,n1)
      bface( 5,n2) = n1
      bface( 6,n2) = bface(2,n1)
      bface( 8,n2) = bface(3,n1)
      bface( 9,n2) = bface(4,n1)
      bface(10,n2) = ncell + n2
      !
    end do
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  contains
  !
  pure function map3dto1d(i,j,k,n) result(return_value)
    integer, intent(in) :: i,j,k,n
    integer :: return_value
    return_value = i + n*(j-1) + n*n*(k-1)
  end function map3dto1d
  !
end subroutine create_TaylorGreenVortex_grid
!
!###############################################################################
!
subroutine create_channel_flow_grid()
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use ovar,      only : Lref,loc_solution_pts,itestcase
  use ovar,      only : channel_grid_ave_xplus
  use ovar,      only : channel_grid_ave_yplus
  use ovar,      only : channel_grid_ave_zplus
  use ovar,      only : channel_grid_yplus_wall
  use ovar,      only : channel_yplus_of_one
  use geovar,    only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geovar,    only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geovar,    only : cell_geom,cell_order
  use ovar,      only : grid_format
  !
  use polynomial_mod, only : legendre_gauss_quadrature
  use polynomial_mod, only : legendre_gauss_lobatto_quadrature
  !
  !.. Local Scalars ..
  integer  :: ierr,i,j,k,l,n,n1,n2
  integer  :: nx,ny,nz,nf_min,nf_max
  real(wp) :: xplus_ave,yplus_ave,zplus_ave,dx,dy,dz
  real(wp) :: ds1,ds2,dl1,dl2,a,b,sf,rnum,rden,sferr,val,ue,se
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: xi(:)
  real(wp), allocatable :: x(:)
  real(wp), allocatable :: y(:)
  real(wp), allocatable :: z(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_channel_flow_grid"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Reset the grid format input parameter to identify
  ! the grid as internally generated
  !
  grid_format = Internally_Generated_Grid
  !
  ! If we are generating the grid on the fly, make the reference length 1
  !
  Lref = one
  !
  ! The channel flow problem is always 3D
  !
  nr = 3
  !
  ! #######
  ! #######
  ! ## 1 ##
  ! #######
  ! #######
  !
  ! Reset any of the grid parameters to their default values if they are zero
  !
  if (channel_yplus_of_one    == zero) channel_yplus_of_one = 8.42134e-6_wp
  if (channel_grid_ave_xplus  == zero) channel_grid_ave_xplus = five*ten
  if (channel_grid_ave_yplus  == zero) channel_grid_ave_yplus = five*five
  if (channel_grid_ave_zplus  == zero) channel_grid_ave_zplus = five*five
  if (channel_grid_yplus_wall == zero) channel_grid_yplus_wall = -one
  !
  ! Compute the cell average x+,y+,z+ for the grid
  ! For a given coordinate direction (say [i]), the cell average [i]+ is
  !   if (channel_grid_ave_[i]plus > zero) then
  !     [i]plus_ave = channel_grid_ave_[i]plus
  !   else
  !     [i]plus_ave = abs(channel_grid_ave_[i]plus)*real(n_order+2,kind=wp)
  !   end if
  ! NOTE: We use (P+2) here because the delta plus value is the
  !       average distance between the DoF within a cell (or the
  !       subcells), and there are P+2 subcells in each cell.
  !
  xplus_ave = abs(channel_grid_ave_xplus)*merge(one, &
                                                real(n_order+2,kind=wp), &
                                                channel_grid_ave_xplus > zero)
  yplus_ave = abs(channel_grid_ave_yplus)*merge(one, &
                                                real(n_order+2,kind=wp), &
                                                channel_grid_ave_yplus > zero)
  zplus_ave = abs(channel_grid_ave_zplus)*merge(one, &
                                                real(n_order+2,kind=wp), &
                                                channel_grid_ave_zplus > zero)
  !
  ! Get the number of grid points in each direction
  !
  nx = nint( channel_grid_Lx / (channel_yplus_of_one*xplus_ave) )
  ny = nint( channel_grid_Ly / (channel_yplus_of_one*yplus_ave) )
  nz = nint( channel_grid_Lz / (channel_yplus_of_one*zplus_ave) )
  !
  ! Compute the x and z cell spacings
  !
  dx = channel_grid_Lx / real(nx,kind=wp)
  dz = channel_grid_Lz / real(nz,kind=wp)
  !
  ! Initialize the x,y,z grid spacing arrays
  !
  x = init_array(nx+1,zero) ! F2003 AUTO-REALLOC
  y = init_array(ny+1,zero) ! F2003 AUTO-REALLOC
  z = init_array(nz+1,zero) ! F2003 AUTO-REALLOC
  !
  ! #######
  ! #######
  ! ## 2 ##
  ! #######
  ! #######
  !
  ! Get the grid cell spacings for x,y,z
  !
  ! ######
  ! # 2A #
  ! ######
  !
  ! Get the grid point locations for x
  !
  do n = 2,nx
    x(n) = x(n-1) + dx
  end do
  !
  ! ######
  ! # 2B #
  ! ######
  !
  ! Get the grid point locations for z
  !
  do n = 2,nz
    z(n) = z(n-1) + dz
  end do
  !
  ! ######
  ! # 2C #
  ! ######
  !
  ! 2C.1 : Get the beginning and ending spacings needed
  !        to generate the tanh spacing
  !
  ! Need to get dl1/ds1 and dl2/ds2 to create tanh spacing for y spacing
  !
  if (channel_grid_yplus_wall < zero) then
    !
    ! The i-th point of the first cell is at y+=1
    !
    i = min( n_order+1 , max( 1 , nint(abs(channel_grid_yplus_wall)) ) )
    !
    ! Get the location of the 1D solution points in reference space
    !
    xi = init_array(n_order+1,zero) ! F2003 AUTO-REALLOC
    !
    if (loc_solution_pts == Legendre_Gauss_Lobatto) then
      i = min(i+1,n_order+1) ! Offset for 1st Lobatto point being at the wall
      call legendre_gauss_lobatto_quadrature(xi)
    else
      call legendre_gauss_quadrature(xi)
    end if
    !
    xi(:) = half*(xi(:) + one)
    !
    dl1 = channel_yplus_of_one / xi(i)
    !
  else
    !
    ! The height of the first cell off the wall is channel_grid_yplus_wall
    !
    dl1 = channel_yplus_of_one * channel_grid_yplus_wall
    !
  end if
  !
  dl2 = dl1
  !
  ! Map the beginning and ending grid cell
  ! spacings to the reference interval [0,1].
  !
  ds1 = dl1 / channel_grid_Ly
  ds2 = dl2 / channel_grid_Ly
  !
  ! 2C.2 : Get the a and b parameters to begin solving for the tanh y-spacing
  !
  a = sqrt(ds2/ds1)
  b = one / (real(ny,kind=wp)*sqrt(ds1*ds2))
  !
  ! 2C.3 : Use b to solve the stretching factor function, then use the
  !        stretching factor to get the y-spacing.
  !
  if (b >= 1.1_wp) then
    !
    ! Use a Newton iteration to solve the stretching factor function
    !    sinh(sf)/sf = b
    ! for the stretching factor, sf.
    !
    sf = two ! initial guess for sf
    !
    newton_tanh_loop: &
    do i = 1,100
      !
      rnum = sinh(sf)/sf - b
      rden = cosh(sf)/sf - sinh(sf)/(sf*sf)
      !
      sferr = rnum / rden
      !
      if (abs(sferr) < eps6) then
        exit newton_tanh_loop
      end if
      !
      sf = sf - sferr
      !
      ! Make sure sf doesnt become 0
      !
      if (abs(sf) < 0.01_wp) then
        sf = 0.01_wp
      end if
      !
    end do newton_tanh_loop
    !
    ! The stretching factor function has two solutions for sf given any
    ! value of b (ignoring sf=0). Since the stretching factor function is
    ! an even function, the two possible solutions are +sf and -sf. We only
    ! care about +sf so take the absolute value of whatever result we get
    ! from the Newton iteration.
    !
    sf = abs(sf)
    !
    ! Get the distribution y(n)
    !
    do n = 0,ny
      val = real(n,kind=wp)/real(ny,kind=wp)
      ue = half*(one + tanh(sf*(val-half))/tanh(sf*half))
      se = ue / (a + (one-a)*ue)
      y(n) = se * channel_grid_Ly
    end do
    !
  else if (b <= 0.9_wp) then
    !
    ! Use a Newton iteration to solve the stretching factor function
    !    sin(sf)/sf = b
    ! for the stretching factor, sf.
    !
    sf = two ! initial guess for sf
    !
    newton_tan_loop: &
    do i = 1,100
      !
      rnum = sin(sf)/sf - b
      rden = cos(sf)/sf - sin(sf)/(sf*sf)
      !
      sferr = rnum / rden
      !
      if (abs(sferr) < eps6) then
        exit newton_tan_loop
      end if
      !
      sf = sf - sferr
      !
      ! Make sure sf doesnt become 0
      !
      if (abs(sf) < 0.01_wp) then
        sf = 0.01_wp
      end if
      !
    end do newton_tan_loop
    !
    ! The stretching factor function has two solutions for sf given any
    ! value of b (ignoring sf=0). Since the stretching factor function is
    ! an even function, the two possible solutions are +sf and -sf. We only
    ! care about +sf so take the absolute value of whatever result we get
    ! from the Newton iteration.
    !
    sf = abs(sf)
    !
    ! Get the distribution y(n)
    !
    do n = 0,ny
      val = real(n,kind=wp)/real(ny,kind=wp)
      ue = half*(one + tan(sf*(val-half))/tan(sf*half))
      se = ue / (a + (one-a)*ue)
      y(n) = se * channel_grid_Ly
    end do
    !
  else
    !
    ! Both of the above Newton iterations can have numerical difficulties
    ! if the stretching factor is close to 1, therefore use this simpler
    ! stretching function that doesnt require a Newton iteration.
    !
    do n = 0,ny
      val = real(n,kind=wp)/real(ny,kind=wp)
      ue = val*(one + two*(b-one)*(val-half)*(one-val))
      se = ue / (a + (one-a)*ue)
      y(n) = se * channel_grid_Ly
    end do
    !
  end if
  !
  ! #######
  ! #######
  ! ## 3 ##
  ! #######
  ! #######
  !
  ! We can finish creating the grid now that we have the x,y,z grid spacings
  !
  ncell = nx*ny*nz                    ! total number of grid cells
  nnode = (nx+1)*(ny+1)*(nz+1)        ! total number of grid nodes
  nfbnd = 2*nx*ny + 2*nx*nz + 2*ny*nz ! total number of boundary faces
  !
  ! Allocate all the arrays needed to define the grid
  !
  allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_nodes",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  allocate ( nodes_of_cell(1:geom_nodes(Geom_Hexa)*ncell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_geom(1:ncell) , source=Geom_Hexa , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_geom",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_order(1:ncell) , source=n_order , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_order",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! #######
  ! #######
  ! ## 4 ##
  ! #######
  ! #######
  !
  ! Create the nodes_of_cell_ptr array
  !
  nodes_of_cell_ptr(1) = 0
  do n = 2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(n) = nodes_of_cell_ptr(n-1) + geom_nodes(Geom_Hexa)
  end do
  !
  ! Create the nodes_of_cell array
  !
  do k = 1,nz
    do j = 1,ny
      do i = 1,nx
        !
        n = nodes_of_cell_ptr( map3dto1d(i,j,k,nx,ny) )
        !
        nodes_of_cell(n+1) = map3dto1d(i  ,j  ,k  ,nx+1,ny+1)
        nodes_of_cell(n+2) = map3dto1d(i+1,j  ,k  ,nx+1,ny+1)
        nodes_of_cell(n+3) = map3dto1d(i+1,j+1,k  ,nx+1,ny+1)
        nodes_of_cell(n+4) = map3dto1d(i  ,j+1,k  ,nx+1,ny+1)
        nodes_of_cell(n+5) = map3dto1d(i  ,j  ,k+1,nx+1,ny+1)
        nodes_of_cell(n+6) = map3dto1d(i+1,j  ,k+1,nx+1,ny+1)
        nodes_of_cell(n+7) = map3dto1d(i+1,j+1,k+1,nx+1,ny+1)
        nodes_of_cell(n+8) = map3dto1d(i  ,j+1,k+1,nx+1,ny+1)
        !
      end do
    end do
  end do
  !
  ! #######
  ! #######
  ! ## 5 ##
  ! #######
  ! #######
  !
  ! Create the grid node coordinates
  !
  do k = 1,nz+1
    do j = 1,ny+1
      do i = 1,nx+1
        !
        n = map3dto1d(i,j,k,nx+1,ny+1)
        !
        xyz_nodes(1,n) = x(i)
        xyz_nodes(2,n) = y(j)
        xyz_nodes(3,n) = z(k)
        !
      end do
    end do
  end do
  !
  ! #######
  ! #######
  ! ## 6 ##
  ! #######
  ! #######
  !
  ! Create the bface array
  !
  nf_max = 0
  !
  ! Loop through the 6 grid boundaries
  ! NOTE: This corresponds to the grid boundaries in
  !       the order [imin,imax,jmin,jmax,kmin,kmax]
  !
  nf_min = nf_max
  !
  do k = 1,nz
    do j = 1,ny
      !
      n1 = nf_min + 1
      n2 = nf_min + 1 + ny*nz
      !
      if (itestcase == -Channel_Flow) then
        !
        bface( 1,n1) = bc_sem
        bface( 2,n1) = map3dto1d(1 ,j,k,nx,ny)
        bface( 3,n1) = Geom_Hexa
        bface( 4,n1) = HostFace3D(1)
        bface( 5,n1) = 1
        bface( 6,n1) = 1
        bface(10,n1) = ncell + n1
        !
        bface( 1,n2) = bc_sub_outflow
        bface( 2,n2) = map3dto1d(nx,j,k,nx,ny)
        bface( 3,n2) = Geom_Hexa
        bface( 4,n2) = HostFace3D(2)
        bface( 5,n2) = 2
        bface( 6,n2) = 2
        bface(10,n2) = ncell + n2
        !
      else
        !
        bface( 1,n1) = bc_periodic
        bface( 2,n1) = map3dto1d(1 ,j,k,nx,ny)
        bface( 3,n1) = Geom_Hexa
        bface( 4,n1) = HostFace3D(1)
        bface( 5,n1) = n2
        bface( 6,n1) = map3dto1d(nx,j,k,nx,ny)
        bface( 8,n1) = Geom_Hexa
        bface( 9,n1) = HostFace3D(2)
        bface(10,n1) = ncell + n1
        !
        bface( 1,n2) = bc_periodic
        bface( 2,n2) = bface(6,n1)
        bface( 3,n2) = bface(8,n1)
        bface( 4,n2) = bface(9,n1)
        bface( 5,n2) = n1
        bface( 6,n2) = bface(2,n1)
        bface( 8,n2) = bface(3,n1)
        bface( 9,n2) = bface(4,n1)
        bface(10,n2) = ncell + n2
        !
      end if
      !
      nf_min = n1
      nf_max = n2
      !
    end do
  end do
  !
  nf_min = nf_max
  !
  do k = 1,nz
    do i = 1,nx
      !
      n1 = nf_min + 1
      n2 = nf_min + 1 + nx*nz
      !
      bface( 1,n1) = bc_periodic
      bface( 2,n1) = map3dto1d(i,1 ,k,nx,ny)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(3)
      bface( 5,n1) = n2
      bface( 6,n1) = map3dto1d(i,ny,k,nx,ny)
      bface( 8,n1) = Geom_Hexa
      bface( 9,n1) = HostFace3D(4)
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_periodic
      bface( 2,n2) = bface(6,n1)
      bface( 3,n2) = bface(8,n1)
      bface( 4,n2) = bface(9,n1)
      bface( 5,n2) = n1
      bface( 6,n2) = bface(2,n1)
      bface( 8,n2) = bface(3,n1)
      bface( 9,n2) = bface(4,n1)
      bface(10,n2) = ncell + n2
      !
      nf_min = n1
      nf_max = n2
      !
    end do
  end do
  !
  nf_min = nf_max
  !
  do j = 1,ny
    do i = 1,nx
      !
      n1 = nf_min + 1
      n2 = nf_min + 1 + nx*ny
      !
      bface( 1,n1) = bc_periodic
      bface( 2,n1) = map3dto1d(i,j,1 ,nx,ny)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(5)
      bface( 5,n1) = n2
      bface( 6,n1) = map3dto1d(i,j,nz,nx,ny)
      bface( 8,n1) = Geom_Hexa
      bface( 9,n1) = HostFace3D(6)
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_periodic
      bface( 2,n2) = bface(6,n1)
      bface( 3,n2) = bface(8,n1)
      bface( 4,n2) = bface(9,n1)
      bface( 5,n2) = n1
      bface( 6,n2) = bface(2,n1)
      bface( 8,n2) = bface(3,n1)
      bface( 9,n2) = bface(4,n1)
      bface(10,n2) = ncell + n2
      !
      nf_min = n1
      nf_max = n2
      !
    end do
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  contains
  !
  pure function map3dto1d(i,j,k,nx,ny) result(return_value)
    integer, intent(in) :: i,j,k,nx,ny
    integer :: return_value
    return_value = i + nx*(j-1) + nx*ny*(k-1)
  end function map3dto1d
  !
end subroutine create_channel_flow_grid
!
!###############################################################################
!
subroutine create_flame_speed_test_grid()
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use ovar,      only : Lref,itestcase
  use ovar,      only : flame_speed_test_dof
  use geovar,    only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geovar,    only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geovar,    only : cell_geom,cell_order
  use ovar,      only : grid_format,bc_names
  !
  !.. Local Scalars ..
  integer  :: n,nc,nf,np,ierr
  real(wp) :: dl
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_flame_speed_test_grid"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Reset the grid format input parameter to identify
  ! the grid as internally generated
  !
  grid_format = Internally_Generated_Grid
  !
  ! Make the reference length one centimeter for the flame speed test case
  !
  Lref = one/(ten*ten)
  !
  ! The flame speed test case is always 1D
  !
  nr = 1
  !
  ! Estimate the number of cells required to get the requested DOF
  !
  nc = nint( real(flame_speed_test_dof,kind=wp) / real(n_order+1,kind=wp) )
  !
  np = nc + 1 ! number of grid nodes in each direction
  !
  ncell = nc    ! total number of grid cells
  nnode = np    ! total number of grid nodes
  nfbnd = 2     ! total number of boundary faces
  !
  ! Allocate all the arrays needed to define the grid
  !
  allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_nodes",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  allocate ( nodes_of_cell(1:geom_nodes(Geom_Edge)*ncell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_geom(1:ncell) , source=Geom_Edge , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_geom",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_order(1:ncell) , source=n_order , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_order",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! #######
  ! #######
  ! ## 1 ##
  ! #######
  ! #######
  !
  ! Create the nodes_of_cell_ptr array
  !
  nodes_of_cell_ptr(1) = 0
  do n = 2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(n) = nodes_of_cell_ptr(n-1) + geom_nodes(Geom_Edge)
  end do
  !
  ! Create the nodes_of_cell array
  !
  do n = 1,ncell
    nodes_of_cell( nodes_of_cell_ptr(n)+1 ) = n
    nodes_of_cell( nodes_of_cell_ptr(n)+2 ) = n+1
  end do
  !
  ! #######
  ! #######
  ! ## 2 ##
  ! #######
  ! #######
  !
  ! Compute the cell size for a Cartesian/uniform grid of total length Lref
  !
  dl = Lref / real(ncell,kind=wp)
  !
  ! Create the grid node coordinates as a uniform 1D grid from -Lref/2 to Lref/2
  !
  do n = 1,nnode
    xyz_nodes(1,n) = -half*Lref + dl*real(n-1,kind=wp)
  end do
  !
  ! #######
  ! #######
  ! ## 3 ##
  ! #######
  ! #######
  !
  ! Create the bface array
  !
  ! Since this is a 1D grid, there are only 2 boundary faces: i-min and i-max
  !
  ! Create the i-min boundary face
  !
  nf = 1
  !
  bface( 1,nf) = bc_sub_inflow
  bface( 2,nf) = 1      ! host cell index
  bface( 3,nf) = Geom_Edge
  bface( 4,nf) = 1      ! x-min face
  bface( 5,nf) = 1      ! this face belongs to boundary group 1
  bface( 6,nf) = 1      ! index of array bc_in(:) to get BC flow conditions
  bface( 7,nf) = 0      ! not used for physical boundary
  bface( 8,nf) = 0      ! not used for physical boundary
  bface( 9,nf) = 0      ! not used for physical boundary
  bface(10,nf) = ncell + nf ! index of ghost cell/face
  !
  ! Create the i-max boundary face
  !
  nf = 2
  !
  bface( 1,nf) = bc_sub_outflow
  bface( 2,nf) = ncell  ! host cell index
  bface( 3,nf) = Geom_Edge
  bface( 4,nf) = 2      ! x-max face
  bface( 5,nf) = 2      ! this face belongs to boundary group 2
  bface( 6,nf) = 2      ! index of array bc_in(:) to get BC flow conditions
  bface( 7,nf) = 0      ! not used for physical boundary
  bface( 8,nf) = 0      ! not used for physical boundary
  bface( 9,nf) = 0      ! not used for physical boundary
  bface(10,nf) = ncell + nf ! index of ghost cell/face
  !
  ! #######
  ! #######
  ! ## 3 ##
  ! #######
  ! #######
  !
  ! Create the bc_names array
  !
  allocate ( bc_names(1:2) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bc_names",1,__LINE__,__FILE__,ierr,error_message)
  !
  bc_names(1) = "INFLOW"
  bc_names(2) = "OUTFLOW"
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine create_flame_speed_test_grid
!
!###############################################################################
!
subroutine create_postproc_debug_grid
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use ovar,      only : Lref
  use geovar,    only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geovar,    only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geovar,    only : cell_geom,cell_order
  use ovar,      only : VortexGrid_random_perturb
  use ovar,      only : VortexGrid_perturb_percent
  use ovar,      only : grid_format
  !
  !.. Local Scalars ..
  integer :: nc,np,ierr,i,j,k,l,n
  integer :: n1,n2,nf_min,nf_max
  integer :: nc_min,nc_max
  integer :: no_min,no_max
  logical(ldk) :: bnd_face_error
  real(r4) :: x,y,z,dl
  real(r4) :: rbas,radj,rval
  !
  !.. Local Arrays ..
  integer, dimension(1:4)      :: top,bot,ip
  integer, dimension(1:8)      :: cell_nodes
  integer, dimension(1:8,1:24) :: node_order
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: lpoin(:)
  integer,  allocatable :: original_noc(:)
  real(r4), allocatable :: perturb(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_postproc_debug_grid"
  !
  integer, parameter :: cells_per_cpu = 3**3
  integer, parameter :: basic_order(8,3) = reshape( [1,2,3,4,5,6,7,8, &
                                                     2,6,7,3,1,5,8,4, &
                                                     1,5,6,2,4,8,7,3] , &
                                                    shape(basic_order) )
  integer, parameter :: flip(4) = [1,4,3,2]
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Reset the grid format input parameter to identify
  ! the grid as internally generated
  !
  grid_format = Internally_Generated_Grid
  !
  bnd_face_error = fals
  !
  ! Create the node_order array to prepare for randomizing the node orderings
  ! NOTE: We subtract one so node_order can be used as an offset
  !
  n = 0
  do j = 1,3
    bot = basic_order(1:4,j)
    top = basic_order(5:8,j)
    do i = 0,3
      n = n + 1
      node_order(:,n) = [cshift(bot,shift=i),cshift(top,shift=i)] - 1
    end do
    bot = top(flip)
    top = basic_order(flip,j)
    do i = 0,3
      n = n + 1
      node_order(:,n) = [cshift(bot,shift=i),cshift(top,shift=i)] - 1
    end do
  end do
  !
  ! If we are generating the grid on the fly, make the reference length 1
  !
  Lref = one
  !
  nr = 3
  !
  ! Estimate the total number of cells
  !
  nc = nint( real(ncpu*cells_per_cpu,kind=wp)**one3 )
  !
  np = nc + 1 ! number of grid nodes in each direction
  !
  ncell = nc*nc*nc ! total number of grid cells
  nnode = np*np*np ! total number of grid nodes
  nfbnd = 6*nc*nc  ! total number of boundary faces
  !
  ! Compute the cell length for a Cartesian/uniform grid
  !
  dl = real(1,kind=r4) / real(np-1,kind=r4)
  !
  ! Allocate all the arrays needed to define the grid
  !
  allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_nodes",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  allocate ( nodes_of_cell(1:geom_nodes(Geom_Hexa)*ncell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_geom(1:ncell) , source=Geom_Hexa , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_geom",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_order(1:ncell) , source=n_order , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_order",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! ONLY HAVE THE ROOT PROCESSOR CREATE THE GRID SINCE THERE IS A LOT OF
  ! RANDOMIZATION INVOLVED. WE WILL BROADCAST THE RESULTS TO ALL THE OTHER
  ! PROCESSORS AFTERWARDS SO THAT ALL PROCESSORS HAVE THE SAME GRID DATA.
  !
  only_root_creates_grid: if (mypnum == glb_root) then
    !
    ! #######
    ! #######
    ! ## 1 ##
    ! #######
    ! #######
    !
    ! Create the nodes_of_cell_ptr array
    !
    nodes_of_cell_ptr(1) = 0
    do n = 2,size(nodes_of_cell_ptr)
      nodes_of_cell_ptr(n) = nodes_of_cell_ptr(n-1) + geom_nodes(Geom_Hexa)
    end do
    !
    ! Create the nodes_of_cell array
    !
    do k = 1,nc
      do j = 1,nc
        do i = 1,nc
          !
          n = nodes_of_cell_ptr( map3dto1d(i,j,k,nc) )
          !
          nodes_of_cell(n+1) = map3dto1d(i  ,j  ,k  ,np)
          nodes_of_cell(n+2) = map3dto1d(i+1,j  ,k  ,np)
          nodes_of_cell(n+3) = map3dto1d(i+1,j+1,k  ,np)
          nodes_of_cell(n+4) = map3dto1d(i  ,j+1,k  ,np)
          nodes_of_cell(n+5) = map3dto1d(i  ,j  ,k+1,np)
          nodes_of_cell(n+6) = map3dto1d(i+1,j  ,k+1,np)
          nodes_of_cell(n+7) = map3dto1d(i+1,j+1,k+1,np)
          nodes_of_cell(n+8) = map3dto1d(i  ,j+1,k+1,np)
          !
        end do
      end do
    end do
    !
    ! Randomize the node ordering for all cells
    !
    if (VortexGrid_random_perturb) then
      !
      ! Create a copy of nodes_of_cell before it is randomized
      !
      allocate ( original_noc(1:size(nodes_of_cell)) , source=nodes_of_cell , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"original_noc",1,__LINE__,__FILE__,ierr, &
                       error_message)
      !
      do k = 1,nc
        do j = 1,nc
          do i = 1,nc
            !
            n = map3dto1d(i,j,k,nc)
            !
            l = random(size(node_order,dim=2))
            l = min( max(1,l) , size(node_order,dim=2) )
            n1 = nodes_of_cell_ptr(n)+1
            n2 = nodes_of_cell_ptr(n+1)
            nodes_of_cell(n1:n2) = nodes_of_cell( n1 + node_order(:,l) )
            !
          end do
        end do
      end do
      !
    end if
    !
    ! #######
    ! #######
    ! ## 2 ##
    ! #######
    ! #######
    !
    allocate ( perturb(1:nr,1:nnode) , source=real(0,kind=r4) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"perturb",1,__LINE__,__FILE__,ierr,error_message)
    !
    if (VortexGrid_random_perturb) then
      !
      ! Compute the random perturbations of the interior node coordinates
      !
      rbas = real(abs(VortexGrid_perturb_percent),kind=r4) / real(100,kind=r4)
      radj = rbas / real(2,kind=r4)
      !
      ! Randomly perturb the crap out of everything
      !
      do k = 2,np-1
        do j = 2,np-1
          do i = 2,np-1
            !
            n = map3dto1d(i,j,k,np)
            !
            do l = 1,nr
              rval = random(rbas) - radj
              perturb(l,n) = rval
            end do
            !
          end do
        end do
      end do
      !
    end if
    !
    ! Create the grid node coordinates as a cube from (0,0,0) to (1,1,1)
    !
    do k = 1,np
      !
      z = real(k-1,kind=r4)
      !
      do j = 1,np
        !
        y = real(j-1,kind=r4)
        !
        do i = 1,np
          !
          x = real(i-1,kind=r4)
          !
          n = map3dto1d(i,j,k,np)
          !
          xyz_nodes(1,n) = real( x + perturb(1,n) , kind=wp )
          xyz_nodes(2,n) = real( y + perturb(2,n) , kind=wp )
          xyz_nodes(3,n) = real( z + perturb(3,n) , kind=wp )
          !
        end do
      end do
    end do
    !
    ! Deallocate perturb since it is no longer needed
    !
    deallocate ( perturb , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"perturb",2,__LINE__,__FILE__,ierr,error_message)
    !
    ! #######
    ! #######
    ! ## 3 ##
    ! #######
    ! #######
    !
    allocate ( lpoin(1:nnode) , source=0 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"lpoin",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Create the bface array
    !
    nf_max = 0
    !
    ! Loop through the 6 grid boundaries
    ! NOTE: This corresponds to the grid boundaries in
    !       the order [imin,imax,jmin,jmax,kmin,kmax]
    !
    nf_min = nf_max
    !
    do k = 1,nc
      do j = 1,nc
        !
        nf_min = nf_min + 1
        nf_max = nf_min + nc*nc
        !
        nc_min = map3dto1d(1 ,j,k,nc)
        nc_max = map3dto1d(nc,j,k,nc)
        !
        no_min = 0
        no_max = 0
        !
        if (VortexGrid_random_perturb) then
          !
          ! ******************
          ! *** X-Min Face ***
          ! ******************
          !
          n1 = nodes_of_cell_ptr(nc_min)+1
          n2 = nodes_of_cell_ptr(nc_min+1)
          !
          ! Mark the nodes on this face
          !
          lpoin( original_noc(n1 + Xmin3D) ) = 1
          !
          ! Extract the nodes of this cell in the new randomized order
          !
          cell_nodes(1:8) = nodes_of_cell(n1:n2)
          !
          ! Find which side of the host cell is now on the boundary
          !
          find_xmin_face: do n = 1,6
            ! Node offsets for this face of the cell
            ! NOTE: The +1 is because the indexing in side_nodes starts at 0.
            ip(1:4) = side_nodes(1:4,n,Geom_Hexa) + 1
            ! Check if all nodes for this face have been marked
            if (sum(lpoin(cell_nodes(ip))) == 4) then
              no_min = n                ! Save the cell face on the boundary
              lpoin(cell_nodes(ip)) = 0 ! Reset the marked nodes to zero
              exit find_xmin_face
            end if
          end do find_xmin_face
          !
          ! Make sure that no_min was set
          !
          if (no_min == 0) then
            bnd_face_error = true
            write (iout,1) "X-Min"
          end if
          !
          ! ******************
          ! *** X-Max Face ***
          ! ******************
          !
          n1 = nodes_of_cell_ptr(nc_max)+1
          n2 = nodes_of_cell_ptr(nc_max+1)
          !
          ! Mark the nodes on this face
          !
          lpoin( original_noc(n1 + Xmax3D) ) = 1
          !
          ! Extract the nodes of this cell in the new randomized order
          !
          cell_nodes(1:8) = nodes_of_cell(n1:n2)
          !
          ! Find which side of the host cell is now on the boundary
          !
          find_xmax_face: do n = 1,6
            ! Node offsets for this face of the cell
            ! NOTE: The +1 is because the indexing in side_nodes starts at 0.
            ip(1:4) = side_nodes(1:4,n,Geom_Hexa) + 1
            ! Check if all nodes for this face have been marked
            if (sum(lpoin(cell_nodes(ip))) == 4) then
              no_max = n                ! Save the cell face on the boundary
              lpoin(cell_nodes(ip)) = 0 ! Reset the marked nodes to zero
              exit find_xmax_face
            end if
          end do find_xmax_face
          !
          ! Make sure that no_max was set
          !
          if (no_max == 0) then
            bnd_face_error = true
            write (iout,1) "X-Max"
          end if
          !
        else
          !
          no_min = HostFace3D(1)
          no_max = HostFace3D(2)
          !
        end if
        !
        bface( 1,nf_min) = bc_periodic
        bface( 2,nf_min) = nc_min
        bface( 3,nf_min) = Geom_Hexa
        bface( 4,nf_min) = no_min
        bface( 5,nf_min) = nf_max
        bface( 6,nf_min) = nc_max
        bface( 8,nf_min) = Geom_Hexa
        bface( 9,nf_min) = no_max
        bface(10,nf_min) = ncell + nf_min
        !
        bface( 1,nf_max) = bc_periodic
        bface( 2,nf_max) = bface(6,nf_min)
        bface( 3,nf_max) = bface(8,nf_min)
        bface( 4,nf_max) = bface(9,nf_min)
        bface( 5,nf_max) = nf_min
        bface( 6,nf_max) = bface(2,nf_min)
        bface( 8,nf_max) = bface(3,nf_min)
        bface( 9,nf_max) = bface(4,nf_min)
        bface(10,nf_max) = ncell + nf_max
        !
      end do
    end do
    !
    nf_min = nf_max
    !
    do k = 1,nc
      do i = 1,nc
        !
        nf_min = nf_min + 1
        nf_max = nf_min + nc*nc
        !
        nc_min = map3dto1d(i,1 ,k,nc)
        nc_max = map3dto1d(i,nc,k,nc)
        !
        no_min = 0
        no_max = 0
        !
        if (VortexGrid_random_perturb) then
          !
          ! ******************
          ! *** Y-Min Face ***
          ! ******************
          !
          n1 = nodes_of_cell_ptr(nc_min)+1
          n2 = nodes_of_cell_ptr(nc_min+1)
          !
          ! Mark the nodes on this face
          !
          lpoin( original_noc(n1 + Ymin3D) ) = 1
          !
          ! Extract the nodes of this cell in the new randomized order
          !
          cell_nodes(1:8) = nodes_of_cell(n1:n2)
          !
          ! Find which side of the host cell is now on the boundary
          !
          find_ymin_face: do n = 1,6
            ! Node offsets for this face of the cell
            ! NOTE: The +1 is because the indexing in side_nodes starts at 0.
            ip(1:4) = side_nodes(1:4,n,Geom_Hexa) + 1
            ! Check if all nodes for this face have been marked
            if (sum(lpoin(cell_nodes(ip))) == 4) then
              no_min = n                ! Save the cell face on the boundary
              lpoin(cell_nodes(ip)) = 0 ! Reset the marked nodes to zero
              exit find_ymin_face
            end if
          end do find_ymin_face
          !
          ! Make sure that no_min was set
          !
          if (no_min == 0) then
            bnd_face_error = true
            write (iout,1) "Y-Min"
          end if
          !
          ! ******************
          ! *** Y-Max Face ***
          ! ******************
          !
          n1 = nodes_of_cell_ptr(nc_max)+1
          n2 = nodes_of_cell_ptr(nc_max+1)
          !
          ! Mark the nodes on this face
          !
          lpoin( original_noc(n1 + Ymax3D) ) = 1
          !
          ! Extract the nodes of this cell in the new randomized order
          !
          cell_nodes(1:8) = nodes_of_cell(n1:n2)
          !
          ! Find which side of the host cell is now on the boundary
          !
          find_ymax_face: do n = 1,6
            ! Node offsets for this face of the cell
            ! NOTE: The +1 is because the indexing in side_nodes starts at 0.
            ip(1:4) = side_nodes(1:4,n,Geom_Hexa) + 1
            ! Check if all nodes for this face have been marked
            if (sum(lpoin(cell_nodes(ip))) == 4) then
              no_max = n                ! Save the cell face on the boundary
              lpoin(cell_nodes(ip)) = 0 ! Reset the marked nodes to zero
              exit find_ymax_face
            end if
          end do find_ymax_face
          !
          ! Make sure that no_max was set
          !
          if (no_max == 0) then
            bnd_face_error = true
            write (iout,1) "Y-Max"
          end if
          !
        else
          !
          no_min = HostFace3D(3)
          no_max = HostFace3D(4)
          !
        end if
        !
        bface( 1,nf_min) = bc_periodic
        bface( 2,nf_min) = nc_min
        bface( 3,nf_min) = Geom_Hexa
        bface( 4,nf_min) = no_min
        bface( 5,nf_min) = nf_max
        bface( 6,nf_min) = nc_max
        bface( 8,nf_min) = Geom_Hexa
        bface( 9,nf_min) = no_max
        bface(10,nf_min) = ncell + nf_min
        !
        bface( 1,nf_max) = bc_periodic
        bface( 2,nf_max) = bface(6,nf_min)
        bface( 3,nf_max) = bface(8,nf_min)
        bface( 4,nf_max) = bface(9,nf_min)
        bface( 5,nf_max) = nf_min
        bface( 6,nf_max) = bface(2,nf_min)
        bface( 8,nf_max) = bface(3,nf_min)
        bface( 9,nf_max) = bface(4,nf_min)
        bface(10,nf_max) = ncell + nf_max
        !
      end do
    end do
    !
    nf_min = nf_max
    !
    do j = 1,nc
      do i = 1,nc
        !
        nf_min = nf_min + 1
        nf_max = nf_min + nc*nc
        !
        nc_min = map3dto1d(i,j,1 ,nc)
        nc_max = map3dto1d(i,j,nc,nc)
        !
        no_min = 0
        no_max = 0
        !
        if (VortexGrid_random_perturb) then
          !
          ! ******************
          ! *** Z-Min Face ***
          ! ******************
          !
          n1 = nodes_of_cell_ptr(nc_min)+1
          n2 = nodes_of_cell_ptr(nc_min+1)
          !
          ! Mark the nodes on this face
          !
          lpoin( original_noc(n1 + Zmin3D) ) = 1
          !
          ! Extract the nodes of this cell in the new randomized order
          !
          cell_nodes(1:8) = nodes_of_cell(n1:n2)
          !
          ! Find which side of the host cell is now on the boundary
          !
          find_zmin_face: do n = 1,6
            ! Node offsets for this face of the cell
            ! NOTE: The +1 is because the indexing in side_nodes starts at 0.
            ip(1:4) = side_nodes(1:4,n,Geom_Hexa) + 1
            ! Check if all nodes for this face have been marked
            if (sum(lpoin(cell_nodes(ip))) == 4) then
              no_min = n                ! Save the cell face on the boundary
              lpoin(cell_nodes(ip)) = 0 ! Reset the marked nodes to zero
              exit find_zmin_face
            end if
          end do find_zmin_face
          !
          ! Make sure that no_min was set
          !
          if (no_min == 0) then
            bnd_face_error = true
            write (iout,1) "Z-Min"
          end if
          !
          ! ******************
          ! *** Z-Max Face ***
          ! ******************
          !
          n1 = nodes_of_cell_ptr(nc_max)+1
          n2 = nodes_of_cell_ptr(nc_max+1)
          !
          ! Mark the nodes on this face
          !
          lpoin( original_noc(n1 + Zmax3D) ) = 1
          !
          ! Extract the nodes of this cell in the new randomized order
          !
          cell_nodes(1:8) = nodes_of_cell(n1:n2)
          !
          ! Find which side of the host cell is now on the boundary
          !
          find_zmax_face: do n = 1,6
            ! Node offsets for this face of the cell
            ! NOTE: The +1 is because the indexing in side_nodes starts at 0.
            ip(1:4) = side_nodes(1:4,n,Geom_Hexa) + 1
            ! Check if all nodes for this face have been marked
            if (sum(lpoin(cell_nodes(ip))) == 4) then
              no_max = n                ! Save the cell face on the boundary
              lpoin(cell_nodes(ip)) = 0 ! Reset the marked nodes to zero
              exit find_zmax_face
            end if
          end do find_zmax_face
          !
          ! Make sure that no_max was set
          !
          if (no_max == 0) then
            bnd_face_error = true
            write (iout,1) "Z-Max"
          end if
          !
        else
          !
          no_min = HostFace3D(5)
          no_max = HostFace3D(6)
          !
        end if
        !
        bface( 1,nf_min) = bc_periodic
        bface( 2,nf_min) = nc_min
        bface( 3,nf_min) = Geom_Hexa
        bface( 4,nf_min) = no_min
        bface( 5,nf_min) = nf_max
        bface( 6,nf_min) = nc_max
        bface( 8,nf_min) = Geom_Hexa
        bface( 9,nf_min) = no_max
        bface(10,nf_min) = ncell + nf_min
        !
        bface( 1,nf_max) = bc_periodic
        bface( 2,nf_max) = bface(6,nf_min)
        bface( 3,nf_max) = bface(8,nf_min)
        bface( 4,nf_max) = bface(9,nf_min)
        bface( 5,nf_max) = nf_min
        bface( 6,nf_max) = bface(2,nf_min)
        bface( 8,nf_max) = bface(3,nf_min)
        bface( 9,nf_max) = bface(4,nf_min)
        bface(10,nf_max) = ncell + nf_max
        !
      end do
    end do
    !
    if (allocated(lpoin)) then
      deallocate ( lpoin , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"lpoin",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    if (allocated(original_noc))  then
      deallocate ( original_noc , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"original_noc",2,__LINE__,__FILE__,ierr, &
                       error_message)
    end if
    !
  end if only_root_creates_grid
  !
  ! Have the root processor broadcast the value of bnd_face_error so that
  ! the simulation can be aborted before continuing if there was an error
  !
  if (ncpu > 1) then
    call mpi_bcast(bnd_face_error,1_int_mpi, &
                   MPI_LOGICAL,0_int_mpi,MPI_COMM_WORLD,mpierr)
  end if
  !
  if (bnd_face_error) then
    write (error_message,2)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Have the root processor broadcast all the grid
  ! information to all the other processors.
  !
  if (ncpu > 1) then
    !
    ! Call mpi_barrier just to make sure the processor are syncronized
    !
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    !
    call mpi_bcast(xyz_nodes,size(xyz_nodes,kind=int_mpi), &
                   mpi_flttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
    !
    call mpi_bcast(nodes_of_cell_ptr,size(nodes_of_cell_ptr,kind=int_mpi), &
                   mpi_inttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
    !
    call mpi_bcast(nodes_of_cell,size(nodes_of_cell,kind=int_mpi), &
                   mpi_inttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
    !
    call mpi_bcast(bface,size(bface,kind=int_mpi), &
                   mpi_inttyp,0_int_mpi,MPI_COMM_WORLD,mpierr)
    !
    ! Call mpi_barrier just to make sure the processor are syncronized
    !
    call mpi_barrier(MPI_COMM_WORLD,mpierr)
    !
  end if
  !
  ! Scale the grid node coordinates to the cube (0,0,0) to (1,1,1)
  !
 !xyz_nodes = dl * xyz_nodes
  !
  ! Scale the grid node coordinates to the cube (-1,-1,-1) to (1,1,1)
  !
 !xyz_nodes = two*xyz_nodes - one
  !
  ! Scale the grid node coordinates to the cube (-pi,-pi,-pi) to (pi,pi,pi)
  !
 !xyz_nodes = PI * xyz_nodes
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (" Error finding new side of host cell on ",a," boundary face!")
  2 format ("An error was detected finding the side of one or more host ", &
            "cells that is on the boundary after randomizing the order of ", &
            "the nodes defining each cell!")
  !
  contains
  !
  pure function map3dto1d(i,j,k,n) result(return_value)
    integer, intent(in) :: i,j,k,n
    integer :: return_value
    return_value = i + n*(j-1) + n*n*(k-1)
  end function map3dto1d
  !
end subroutine create_postproc_debug_grid
!
!###############################################################################
!
end module internal_grid_mod
