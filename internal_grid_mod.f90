module internal_grid_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  implicit none
  !
  private
  !
  public :: create_internal_grid
  !
  logical(lk), parameter :: randomly_rotate_cells_debug = fals
 !logical(lk), parameter :: randomly_rotate_cells_debug = true
  !
  integer, parameter :: tet_random_rotation(1:4,1:12) = &
                                           reshape( [1,2,3,4, &
                                                     2,3,1,4, &
                                                     3,1,2,4, &
                                                     4,3,2,1, &
                                                     3,2,4,1, &
                                                     2,4,3,1, &
                                                     1,4,2,3, &
                                                     4,2,1,3, &
                                                     2,1,4,3, &
                                                     1,3,4,2, &
                                                     3,4,1,2, &
                                                     4,1,3,2] , [4,12] )
  !
  integer, parameter :: pyr_random_rotation(1:5,1:4) = &
                                           reshape( [1,2,3,4,5, &
                                                     2,3,4,1,5, &
                                                     3,4,1,2,5, &
                                                     4,1,2,3,5] , [5,4] )
  !
  integer, parameter :: pri_random_rotation(1:6,1:6) = &
                                           reshape( [1,2,3,4,5,6, &
                                                     2,3,1,5,6,4, &
                                                     3,1,2,6,4,5, &
                                                     6,5,4,3,2,1, &
                                                     5,4,6,2,1,3, &
                                                     4,6,5,1,3,2] , [6,6] )
  !
  integer, parameter :: hex_random_rotation(1:8,1:24) = &
                                           reshape( [1,2,3,4,5,6,7,8, &
                                                     2,3,4,1,6,7,8,5, &
                                                     3,4,1,2,7,8,5,6, &
                                                     4,1,2,3,8,5,6,7, &
                                                     8,7,6,5,4,3,2,1, &
                                                     7,6,5,8,3,2,1,4, &
                                                     6,5,8,7,2,1,4,3, &
                                                     5,8,7,6,1,4,3,2, &
                                                     1,4,8,5,2,3,7,6, &
                                                     4,8,5,1,3,7,6,2, &
                                                     8,5,1,4,7,6,2,3, &
                                                     5,1,4,8,6,2,3,7, &
                                                     6,7,3,2,5,8,4,1, &
                                                     7,3,2,6,8,4,1,5, &
                                                     3,2,6,7,4,1,5,8, &
                                                     2,6,7,3,1,5,8,4, &
                                                     1,5,6,2,4,8,7,3, &
                                                     5,6,2,1,8,7,3,4, &
                                                     6,2,1,5,7,3,4,8, &
                                                     2,1,5,6,3,4,8,7, &
                                                     3,7,8,4,2,6,5,1, &
                                                     7,8,4,3,6,5,1,2, &
                                                     8,4,3,7,5,1,2,6, &
                                                     4,3,7,8,1,2,6,5] , [8,24] )
  !
  interface tanh_spacing
    module procedure tanh_spacing_1d, tanh_spacing_md
  end interface tanh_spacing
  !
  interface average_1p_spacing
    module procedure average_1p_spacing_1d, average_1p_spacing_md
  end interface average_1p_spacing
  !
contains
!
!###############################################################################
!
subroutine create_internal_grid()
  !
  !.. Use Statements ..
  use ovar, only : itestcase,postproc_debug_grid
  use ovar, only : VortexGrid_is_3D,VortexGrid_Hex_to_Pyramids
  use ovar, only : flame_speed_test_dof
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
  else if (any(abs(itestcase) == Transport_Problems) .or. &
           itestcase == Taylor_Green_Vortex) then
    !
    if (itestcase == Taylor_Green_Vortex .or. VortexGrid_is_3D) then
      !
      if (VortexGrid_Hex_to_Pyramids) then
        call create_periodic_cube_pyramid_grid
      else
        call create_periodic_cube_hex_grid
      end if
      !
    else
      !
      call create_periodic_square_quad_grid
      !
    end if
    !
  else if (abs(itestcase) == Taylor_Couette_Flow) then
    !
    call create_taylor_couette_grid
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
subroutine create_taylor_couette_grid()
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use ovar,      only : Lref,bc_input
  use geovar,    only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geovar,    only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geovar,    only : grid,cell_geom,cell_order
  use ovar,      only : grid_format,bc_names
  !
  use ovar, only : TCF_Inner_Radius,TCF_Inner_Wall_Spacing
  use ovar, only : TCF_Outer_Radius,TCF_Outer_Wall_Spacing
  use ovar, only : TCF_Aspect_Ratio,TCF_Axial_Cells
  use ovar, only : TCF_Radial_Cells,TCF_Azimuthal_Cells
  !
  use module_cgns_types, only : HEXA_8,HEXA_27,HEXA_64,HEXA_125
  !
  use module_cgns, only : get_face_elem
  !
  !.. Local Scalars ..
  integer  :: ierr,i,j,k,imin,imax,jmin,jmax
  integer  :: nc,nf,n1,n2,nc_min,nc_max
  integer  :: ncx,ncr,nca,npx,npr,npa,npts
  integer  :: host_cell,host_side,elem_type
  integer  :: inner_wall_idx,outer_wall_idx
  integer  :: inflow_idx,outflow_idx
  real(wp) :: lx,special_b_coeff
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: elem_nodes(:,:)
  integer,  allocatable :: elem_pts(:,:)
  real(wp), allocatable :: radial_pts(:)
  real(wp), allocatable :: xyz_pts(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_taylor_couette_grid"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
 !bc_names = ["Inflow","InnerWall","OuterWall","Outflow"]
 !bc_names = ["InnerWall","OuterWall","Inflow","Outflow"]
  !
  inner_wall_idx = -1
  outer_wall_idx = -1
  inflow_idx = -1
  outflow_idx = -1
  !
  bc_names = bc_input(1:4)%name ! F2003 AUTO-REALLOC
  !
  do i = 1,size(bc_names)
    if (trim(adjustl(uppercase(bc_names(i)))) == "INNERWALL") then
      inner_wall_idx = i
    else if (trim(adjustl(uppercase(bc_names(i)))) == "OUTERWALL") then
      outer_wall_idx = i
    else if (trim(adjustl(uppercase(bc_names(i)))) == "INFLOW") then
      inflow_idx = i
    else if (trim(adjustl(uppercase(bc_names(i)))) == "OUTFLOW") then
      outflow_idx = i
    end if
  end do
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
  ! Get the length of the cylinders
  !
  lx = TCF_Aspect_Ratio * abs(TCF_Outer_Radius-TCF_Inner_Radius)
  !
  ! The Taylor-Couette flow problem is always 3D
  !
  nr = 3
  !
  ! Get the number of grid cells in each direction from the TCF inputs
  !
  ncr = TCF_Radial_Cells
  nca = TCF_Azimuthal_Cells
  ncx = TCF_Axial_Cells
  !
  ! Get the number of grid nodes in each direction
  !
  npr = ncr+1 ! radial nodes
  npa = nca   ! azimuthal nodes NOTE: self-connecting so no +1
  npx = ncx+1 ! axial nodes
  !
  ncell = ncr*nca*ncx ! total number of grid cells
  ! NOTE: The following is still true, but we have changed npa=nca (instead of
  !       npa=nca+1) which simplifies the logic for a bunch of the following
  !       code. This means anywhere it says to use nca and NOT npa, we can now
  !       use npa in those places since npa==nca.
  ! We have to be careful when computing the number of grid nodes for this grid.
  ! Because the grid in the azimuthal direction reconnects with itself,
  ! the number of grid nodes isnt just npr*npa*npx like for other simple
  ! structured grids. To get the total number of grid nodes, we need to multiply
  ! the number of nodes on an axial plane is (npr*nca) by the number of axial
  ! planes (number axial nodes).
  nnode = npr*npa*npx ! total number of grid nodes (NOTE use nca NOT npa)
  ! total number of boundary faces
  nfbnd = ncx*nca & ! number of boundary faces for inner cylinder
        + ncx*nca & ! number of boundary faces for outer cylinder
        + ncr*nca & ! number of boundary faces for  inflow/periodic-A
        + ncr*nca   ! number of boundary faces for outflow/periodic-B
  !
  ! Allocate all the arrays needed to define the grid
  ! NOTE: xyz_nodes and elem_nodes are allocated in high_order_TCF_grid_points
  !
  nodes_of_cell_ptr = init_array(ncell+1,0) ! F2003 AUTO-REALLOC
  !
  nodes_of_cell = init_array(geom_nodes(Geom_Hexa)*ncell,0) ! F03 AUTO-REALLOC
  !
  bface = init_array(nbfai,nfbnd,0) ! F2003 AUTO-REALLOC
  !
  cell_geom = init_array(ncell,Geom_Hexa) ! F2003 AUTO-REALLOC
  !
  cell_order = init_array(ncell,n_order) ! F2003 AUTO-REALLOC
  !
  ! #######
  ! #######
  ! ## 1 ##
  ! #######
  ! #######
  !
  ! Compute the radial cell spacings using a tanh distribution
  !
  special_b_coeff = -one
  check_recreate_loop: &
  do i = 1,1
    !
    ! This loop is to check if the inputs used for creating the grid for
    ! the Taylor-Couette flow problem match the specs for the initial
    ! grid created by Pointwise to test/debug this problem. If the inputs
    ! match the Pointwise grid specs, it is assumed we are trying to recreate
    ! that grid for testing the use of an internally generated grid versus
    ! the result of using the Pointwise grid. The radial tanh spacing
    ! generated by Pointwise doesnt match the spacings generated by the
    ! tanh_spacing function (with default settings) using the same input
    ! values. We can recreate the Pointwise tanh spacing by providing a
    ! special coefficient used in the computation of the b parameter that is
    ! needed to compute the tanh spacing. We can specify this special
    ! coefficient using the optional argument special_b_coeff to the function
    ! tanh_spacing. Using a value of +3.75 for special_b_coeff seems to do
    ! a good enough job of recreating the Pointwise tanh spacings.
    !
    ! If any of the TCF input values do not match those used for the initial
    ! Pointwise grid, exit this single-pass loop and leave special_b_coeff to
    ! its default value of -1.
    !
    if (TCF_Radial_Cells    /= 20) exit check_recreate_loop
    if (TCF_Azimuthal_Cells /= 48) exit check_recreate_loop
    if (TCF_Axial_Cells     /= 24) exit check_recreate_loop
    if (abs(TCF_Inner_Radius       - half)      > eps6) exit check_recreate_loop
    if (abs(TCF_Inner_Wall_Spacing - 0.002_wp)  > eps6) exit check_recreate_loop
    if (abs(TCF_Outer_Radius       - one)       > eps6) exit check_recreate_loop
    if (abs(TCF_Outer_Wall_Spacing - 0.0025_wp) > eps6) exit check_recreate_loop
    !
    ! It appears we are trying to recreate the initial Pointwise grid used
    ! for testing/debugging the Taylor-Couette flow problem, so set
    ! special_b_coeff to +3.75 which should do a decent job of recreating the
    ! radial tanh spacings generated by Pointwise.
    !
    special_b_coeff = 3.75_wp
    !
  end do check_recreate_loop
  !
  ! Compute the actual radial point coordinates using the tanh_spacing function
  !
  radial_pts = tanh_spacing(TCF_Inner_Radius, & ! F2003 AUTO-REALLOC
                            TCF_Outer_Radius, &
                            TCF_Inner_Wall_Spacing, &
                            TCF_Outer_Wall_Spacing, &
                            npr,special_b_coeff=special_b_coeff)
  !
  ! #######
  ! #######
  ! ## 2 ##
  ! #######
  ! #######
  !
  ! Create the linear grid cell connectivities (elem_nodes), the high-order
  ! grid cell connectivities (elem_pts), the linear grid node coordinates
  ! (xyz_nodes), and the high-order grid point coordinates (xyz_pts).
  !
  call high_order_TCF_grid_points(ncx,ncr,nca,lx,radial_pts,n_order, &
                                  elem_nodes,elem_pts,xyz_nodes,xyz_pts)
  !
  ! Deallocate radial_pts since it is no longer needed
  !
  if (allocated(radial_pts)) then
    deallocate ( radial_pts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"radial_pts",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! #######
  ! #######
  ! ## 3 ##
  ! #######
  ! #######
  !
  ! Create the nodes_of_cell and nodes_of_cell_ptr arrays. nodes_of_cell is
  ! just a 1D version of elem_nodes, and nodes_of_cell_ptr is an indexing
  ! array that gives the beginning/ending index in nodes_of_cell for each
  ! cell.
  !
  ! Create the nodes_of_cell_ptr array
  !
  nodes_of_cell_ptr = intseq(0,ncell)*geom_nodes(Geom_Hexa) ! F2003 AUTO-REALLOC
 !nodes_of_cell_ptr(1) = 0
 !do n = 2,size(nodes_of_cell_ptr)
 !  nodes_of_cell_ptr(n) = nodes_of_cell_ptr(n-1) + geom_nodes(Geom_Hexa)
 !end do
  !
  ! Create the nodes_of_cell array
  !
  ! F2003 AUTO-REALLOC
  nodes_of_cell = reshape( elem_nodes(:,:) , shape=[size(elem_nodes)] )
 !do n = 1,ncell
 !  !
 !  n0 = nodes_of_cell_ptr(n)
 !  !
 !  nodes_of_cell(n0+1) = elem_nodes(1,n)
 !  nodes_of_cell(n0+2) = elem_nodes(2,n)
 !  nodes_of_cell(n0+3) = elem_nodes(3,n)
 !  nodes_of_cell(n0+4) = elem_nodes(4,n)
 !  nodes_of_cell(n0+5) = elem_nodes(5,n)
 !  nodes_of_cell(n0+6) = elem_nodes(6,n)
 !  nodes_of_cell(n0+7) = elem_nodes(7,n)
 !  nodes_of_cell(n0+8) = elem_nodes(8,n)
 !  !
 !end do
  !
  ! #######
  ! #######
  ! ## 4 ##
  ! #######
  ! #######
  !
  ! Create the bface array
  !
  ! Loop through the 4 grid boundaries
  ! NOTE: This corresponds to the grid boundaries in
  !       the order [InnerWall,OuterWall,Inflow,Outflow] which
  !       correspond to the [jmin,jmax,imin,imax] boundaries.
  ! NOTE: There are no kmin and kmax boundaries because
  !       the k-direction self-connects with itself.
  !
  ! InnerWall/j-min and OuterWall/j-max boundaries
  !
  n1 = 0            ! j-min counter
  n2 = n1 + ncx*nca ! j-max counter
  do k = 1,nca
    do i = 1,ncx
      !
      jmin = 1
      jmax = ncr
      !
      nc_min = map3dto1d(i,jmin,k,ncx,ncr,nca)
      nc_max = map3dto1d(i,jmax,k,ncx,ncr,nca)
      !
      n1 = n1 + 1 ! j-min face
      n2 = n2 + 1 ! j-max face
      !
      bface( 1,n1) = bc_isothermal_wall
      bface( 2,n1) = nc_min
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(3)
      bface( 5,n1) = inner_wall_idx ! group for boundary face
      bface( 6,n1) = inner_wall_idx ! index of bc_input for this BC
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_adiabatic_wall
      bface( 2,n2) = nc_max
      bface( 3,n2) = Geom_Hexa
      bface( 4,n2) = HostFace3D(4)
      bface( 5,n2) = outer_wall_idx ! group for boundary face
      bface( 6,n2) = outer_wall_idx ! index of bc_input for this BC
      bface(10,n2) = ncell + n2
      !
    end do
  end do
  !
  ! Min-X and Max-X boundary faces
  !
  n1 = 2*ncx*nca    ! min-x counter
  n2 = n1 + ncr*nca ! max-x counter
  do k = 1,nca
    do j = 1,ncr
      !
      imin = 1
      imax = ncx
      !
      nc_min = map3dto1d(imin,j,k,ncx,ncr,nca)
      nc_max = map3dto1d(imax,j,k,ncx,ncr,nca)
      !
      n1 = n1 + 1 ! x-min face
      n2 = n2 + 1 ! x-max face
      !
      ! NOTE: The subroutine match_periodic_faces will set the values of
      !       bface([5,6,8,9],:), so we dont need to do that here.
      !
      bface( 1,n1) = bc_periodic
      bface( 2,n1) = nc_min
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(1)
      bface( 5,n1) = inflow_idx ! group for boundary face
      bface( 6,n1) = inflow_idx ! idx needed by match_periodic_faces
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_periodic
      bface( 2,n2) = nc_max
      bface( 3,n2) = Geom_Hexa
      bface( 4,n2) = HostFace3D(2)
      bface( 5,n2) = outflow_idx ! group for boundary face
      bface( 6,n2) = outflow_idx ! idx needed by match_periodic_faces
      bface(10,n2) = ncell + n2
      !
    end do
  end do
  !
  ! #######
  ! #######
  ! ## 5 ##
  ! #######
  ! #######
  !
  ! We need to create the grid derived type for the Taylor-Couette grids
  ! because we need to make the high-order grid points for the curved
  ! boundaries
  !
  if (allocated(grid)) then
    deallocate ( grid , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( grid , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( grid%elem(1:ncell+nfbnd) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid%elem",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create grid for the grid cells
  !
  if (n_order == 1) then
    elem_type = int(HEXA_8)
    npts = 8
  else if (n_order == 2) then
    elem_type = int(HEXA_27)
    npts = 27
  else if (n_order == 3) then
    elem_type = int(HEXA_64)
    npts = 64
  else ! n_order >= 4
    elem_type = int(HEXA_125)
    npts = 125
  end if
  !
  do nc = 1,ncell
    !
    grid%elem(nc)%prop%npts = npts
    grid%elem(nc)%prop%order = min(n_order,4)
    grid%elem(nc)%prop%geom = Geom_Hexa
    grid%elem(nc)%prop%elem_type = elem_type
    !
    grid%elem(nc)%nodes = elem_nodes(1:8,nc) ! F2003 AUTO-REALLOC
    grid%elem(nc)%pts = elem_pts(1:npts,nc) ! F2003 AUTO-REALLOC
    !
    grid%elem(nc)%tags = [0,0] ! F2003 AUTO-REALLOC
    !
  end do
  !
  ! Create grid for the boundary faces
  !
  do nf = 1,nfbnd
    !
    nc = ncell + nf
    !
    host_cell = bface(2,nf)
    host_side = bface(4,nf)
    !
    ! Use get_face_elem to get the face connectivity
    ! NOTE: This will assign all the %prop values, create the %nodes
    !       and %pts arrays, and copy the host cells %tags array.
    ! NOTE: The get_face_elem function works correctly whether the host
    !       cells point connectivity is in CGNS or GFR ordering. It will check
    !       to see if
    !         1) %pts(1:size(%nodes)) == %nodes (CGNS ordering), and then
    !            return the face elements %pts array in GFR ordering, or
    !         2) otherwise, it will reverse the GFR ordering of the host
    !            cells %pts to put it back in CGNS order, the return the
    !            face elements %pts array in GFR ordering
    !
    grid%elem(nc) = get_face_elem(grid%elem(host_cell),host_side,host_cell)
    !
    ! The get_face_elem function copies the host cells tags array to the
    ! face element, so we need to add the boundary condition to tags.
    !
    grid%elem(nc)%tags = [bface(1,nf),0] ! F2003 AUTO-REALLOC
    !
  end do
  !
  ! Copy the xyz_pts array into grid%xyz using the move_alloc subroutine
  ! NOTE: This will allocate grid%xyz with the same size and shape as xyz_pts,
  !       copy over the values, then deallocate xyz_pts.
  !
  call move_alloc( from=xyz_pts , to=grid%xyz )
  !
 !call write_internal_grid(ncell,nfbnd)
  !
  ! Deallocate elem_nodes and elem_pts before leaving
  !
  if (allocated(elem_nodes)) then
    deallocate ( elem_nodes , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"elem_nodes",2,__LINE__,__FILE__,ierr,error_message)
  end if
  if (allocated(elem_pts)) then
    deallocate ( elem_pts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"elem_pts",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  contains
  !
  pure function map3dto1d(i,j,k,ni,nj,nk) result(return_value)
    integer, intent(in) :: i,j,k,ni,nj,nk
    integer :: return_value
    integer :: kk
    ! If k>nk, this should mean that we are looking at a cell on the
    ! self-connecting azimuthal plane. For these cells, the k+1 nodes are
    ! actually the k=1 nodes so zero out the k part of the mapping.
    kk = merge( 0 , ni*nj*(k-1) , k > nk )
    return_value = i + ni*(j-1) + kk
  end function map3dto1d
  !
end subroutine create_taylor_couette_grid
!
!###############################################################################
!
subroutine high_order_TCF_grid_points(ncx,ncr,nca,length,radial_pts,order, &
                                      elem_nodes,elem_pts,xyz_nodes,xyz_pts)
  !
  !.. Formal Arguments ..
  integer,                intent(in) :: ncx
  integer,                intent(in) :: ncr
  integer,                intent(in) :: nca
  real(wp),               intent(in) :: length
  real(wp),               intent(in) :: radial_pts(:)
  integer,                intent(in) :: order
  integer,  allocatable, intent(out) :: elem_nodes(:,:)
  integer,  allocatable, intent(out) :: elem_pts(:,:)
  real(wp), allocatable, intent(out) :: xyz_nodes(:,:)
  real(wp), allocatable, intent(out) :: xyz_pts(:,:)
  !
  !.. Local Scalars ..
  integer  :: nnx,nnr,nna,npx,npr,npa,p
  integer  :: i,j,k,ic,jc,kc,ip,jp,kp
  integer  :: n,nc,lcell,lnode,cpts,npts
  real(wp) :: y,z,theta,divp,dr,dlpa,dlpx
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: ijkmap(:,:,:)
  real(wp), allocatable :: rpts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "high_order_TCF_grid_points"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Number of grid nodes in each direction
  !
  nnx = ncx + 1 ! axial     grid nodes
  nnr = ncr + 1 ! radial    grid nodes
  nna = nca     ! azimuthal grid nodes
  !
  ! Get the order of the high order grid
  ! NOTE: For now we are limiting this to P4 to match the maximum allowed
  !       by CGNS. One reason to do this is because we can then use the
  !       functions in cgns_mod to extract the high-order connectivity for
  !       each of the boundary face elements, which are needed for the
  !       grid derived type.
  !
  p = min(order,4)
  !
  ! Number of high order points in each direction
  !
  npx = ncx*p + 1 ! axial     high order points
  npr = ncr*p + 1 ! radial    high order points
  npa = nca*p     ! azimuthal high order points
  !
  !################################################################
  !################################################################
  !#####  First, we need to take the radial cell spacings in  #####
  !#####  radial_pts and add in the uniform spacings for the  #####
  !#####  high-order points between each of the grid nodes.   #####
  !################################################################
  !################################################################
  !
  rpts = init_array(npr,zero) ! F2003 AUTO-REALLOC
  !
  jp = 1
  rpts(jp) = radial_pts(1)
  divp = one/real(p,kind=wp)
  do jc = 2,ncr+1
    dr = divp * (radial_pts(jc)-radial_pts(jc-1))
    do n = 1,p
      jp = jp + 1
      rpts(jp) = rpts(jp-1) + dr
    end do
  end do
  !
  ! Also, compute the uniform high-order points spacings
  !
  ! Compute the uniform axial high-order point spacing
  dlpx = length / real(npx-1,kind=wp)
  ! Compute the uniform azimuthal high-order point spacing
  dlpa = two*pi / real(npa,kind=wp)
  !
  !###################################################################
  !###################################################################
  !#####  Next, we need to create the ijkmap array which allows  #####
  !#####  us to map between the high-order i,j,k point indices   #####
  !#####  to the respective 1D index. Normally, we would just    #####
  !#####  use a map3dto1d function, but the points need to be    #####
  !#####  ordered such that the all the linear grid nodes come   #####
  !#####  first followed by all the high-order points (well I    #####
  !#####  am not sure the NEED to be ordered this way, but this  #####
  !#####  seems to be the ordering convention used by others).   #####
  !#####  While we are in the loop creating ijkmap, we will      #####
  !#####  also create xyz_nodes since we have all the necessary  #####
  !#####  information available and it uses the exact same loop  #####
  !#####  structures.                                            #####
  !###################################################################
  !###################################################################
  !
  ! Compute the total number of grid nodes and allocate xyz_nodes
  !
  lnode = nnx*nnr*nna
  !
  xyz_nodes = init_array(3,lnode,zero) ! F2003 AUTO-REALLOC
  !
  ! Allocate ijkmap
  ! NOTE: The ijkmap(:,:,npa+1) == ijkmap(:,:,1) because
  !       the min/max k boundaries self-connect
  ! NOTE: We initialize this to -1 so that we will know which i,j,k
  !       points have already been given a 1D index
  !
  ijkmap = init_array(npx,npr,npa+1,-1)
  !
  ! First create the indices for the linear grid nodes, and also
  ! fill in xyz_nodes for these nodes
  !
  n = 0
  do k = 1,nna
    !
    kp = (k-1)*p + 1
    !
    theta = dlpa * real(kp-1,kind=wp)
    !
    do j = 1,nnr
      !
      jp = (j-1)*p + 1
      !
      y = rpts(jp)*cos(theta)
      z = rpts(jp)*sin(theta)
      !
      do i = 1,nnx
        !
        ip = (i-1)*p + 1
        !
        ! Save the 1D index for this grid node to ijkmap
        !
        n = n + 1
        !
        ijkmap(ip,jp,kp) = n
        !
        ! Compute the coordinates for this point and save it to xyz_nodes
        !
        xyz_nodes(1,n) = dlpx * real(ip-1,kind=wp)
        xyz_nodes(2,n) = y
        xyz_nodes(3,n) = z
        !
      end do
    end do
  end do
  !
  ! Compute the total number of high-order points and allocate xyz_pts
  !
  npts = npx*npr*npa
  !
  xyz_pts = init_array(3,npts,zero) ! F2003 AUTO-REALLOC
  !
  ! Copy xyz_nodes over to xyz_pts
  ! NOTE: n should equal lnode coming out of the above loops
  !
  xyz_pts(1:3,1:lnode) = xyz_nodes(1:3,1:lnode)
  !
  ! Now create the indices for the high-order grid points, and also
  ! fill in xyz_pts for these points
  !
  n = lnode ! n should already be equal to lnode, but just in case
  do kp = 1,npa
    theta = dlpa * real(kp-1,kind=wp)
    do jp = 1,npr
      do ip = 1,npx
        !
        ! If the value in ijkmap for the current point is not the
        ! initialized value of -1, then this point is a linear grid
        ! node and has already been given a 1D index, so skip to the
        ! next high-order point.
        !
        if (ijkmap(ip,jp,kp) /= -1) cycle
        !
        ! This high-order point is not a linear grid nodes, so give it a
        ! 1D index in ijkmap
        !
        n = n + 1
        !
        ijkmap(ip,jp,kp) = n
        !
        ! Compute the coordinates for this point and save it to xyz_pts
        !
        xyz_pts(1,n) = dlpx * real(ip-1,kind=wp)
        xyz_pts(2,n) = rpts(jp)*cos(theta)
        xyz_pts(3,n) = rpts(jp)*sin(theta)
        !
      end do
    end do
  end do
  !
  ! Finally fill the self-connected kmax plane in ijkmap
  ! by copying over the indices from the kmin plane
  !
  ijkmap(:,:,npa+1) = ijkmap(:,:,1)
  !
  !####################################################################
  !####################################################################
  !#####  Finally, create the elem_nodes and elem_pts arrays.     #####
  !#####  The elem_nodes array contains the node connectivity     #####
  !#####  for each of the linear grid cells, and the elem_pts     #####
  !#####  array contains the high-order point connectivity for    #####
  !#####  each of the high-order grid cells. The connectivity     #####
  !#####  in elem_nodes is just a 2D array version of the same    #####
  !#####  connectivity stored in nodes_of_cell. The connectivity  #####
  !#####  in elem_pts uses the GFR point ordering, (i.e., i,j,k   #####
  !#####  ordering).                                              #####
  !####################################################################
  !####################################################################
  !
  ! Compute the total number of grid cells and allocate elem_nodes
  !
  lcell = ncx*ncr*nca
  !
  elem_nodes = init_array(8,lcell,0) ! F2003 AUTO-REALLOC
  !
  ! Compute the number of high order points in each grid cell
  ! and allocate elem_pts.
  !
  cpts = (p+1)**3
  !
  elem_pts = init_array(cpts,lcell,0) ! F2003 AUTO-REALLOC
  !
  ! Loop over the cells and create the elem_nodes
  ! and elem_pts connectivities for each cell.
  !
  nc = 0
  k_cells: &
  do kc = 1,nca
    !
    k = (kc-1)*p + 1
    !
    j_cells: &
    do jc = 1,ncr
      !
      j = (jc-1)*p + 1
      !
      i_cells: &
      do ic = 1,ncx
        !
       !nc = map3dto1d(ic,jc,kc,ncx,ncr,nca)
        nc = nc + 1
        !
        i = (ic-1)*p + 1
        !
        ! Fill elem_nodes
        !
        elem_nodes(1,nc) = ijkmap(i  ,j  ,k  )
        elem_nodes(2,nc) = ijkmap(i+p,j  ,k  )
        elem_nodes(3,nc) = ijkmap(i+p,j+p,k  )
        elem_nodes(4,nc) = ijkmap(i  ,j+p,k  )
        elem_nodes(5,nc) = ijkmap(i  ,j  ,k+p)
        elem_nodes(6,nc) = ijkmap(i+p,j  ,k+p)
        elem_nodes(7,nc) = ijkmap(i+p,j+p,k+p)
        elem_nodes(8,nc) = ijkmap(i  ,j+p,k+p)
        !
        ! Fill elem_pts
        !
        n = 0
        do kp = k,k+p
          do jp = j,j+p
            do ip = i,i+p
              n = n + 1
              elem_pts(n,nc) = ijkmap(ip,jp,kp)
            end do
          end do
        end do
        !
      end do i_cells
    end do j_cells
  end do k_cells
  !
  call debug_timer(leaving_procedure,pname)
  !
!!contains
!!!
!!pure function map3dto1d(i,j,k,ni,nj,nk) result(return_value)
!!  integer, intent(in) :: i,j,k,ni,nj,nk
!!  integer :: return_value
!!  integer :: kk
!!  ! If k>nk, this should mean that we are looking at a cell on the
!!  ! self-connecting azimuthal plane. For these cells, the k+1 nodes are
!!  ! actually the k=1 nodes so zero out the k part of the mapping.
!!  kk = merge( 0 , ni*nj*(k-1) , k > nk )
!!  return_value = i + ni*(j-1) + kk
!!end function map3dto1d
  !
end subroutine high_order_TCF_grid_points
!
!###############################################################################
!
subroutine write_internal_grid(lcell,lfbnd)
  !
  !.. Use Statements ..
  use geovar, only : nr,grid,bface
  !
  !.. Formal Arguments ..
  integer, intent(in) :: lcell
  integer, intent(in) :: lfbnd
  !
  !.. Local Scalars ..
  integer :: nnod,npts,fnod,fpts
  integer :: nc,nf,i,l
  real(wp) :: x,y,z,r,t
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: nodes(:)
  integer, allocatable :: pts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_internal_grid"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nnod = 0
  npts = 0
  fnod = 0
  fpts = 0
  !
  do nc = 1,lcell
    nnod = max(nnod,size(grid%elem(nc)%nodes))
    npts = max(npts,size(grid%elem(nc)%pts))
  end do
  !
  do nf = 1,lfbnd
    nc = lcell + nf
    fnod = max(fnod,size(grid%elem(nc)%nodes))
    fpts = max(fpts,size(grid%elem(nc)%pts))
  end do
  !
  1 format (i6,1x,"|",1x,100i7)
  2 format (i6,1x,"|",1x,i3,i2,1x,"|",1x,100i7)
  3 format (i7,2(5x,3f12.6))
  4 format (89("="),/,"Cell ",i0, &
            t21,"x",t33,"y",t45,"z",t62,"x",t74,"r",t86,"t",/, &
            89("-"))
  5 format (89("="),/,"Boundary Face ",i0, &
            t21,"x",t33,"y",t45,"z",t62,"x",t74,"r",t86,"t",/, &
            89("-"))
  6 format (1000(a1))
  if (mypnum == glb_root) then
    !
    write (4101,1) 0,    (i, i=1,nnod)
    write (4102,2) 0,0,0,(i, i=1,fnod)
    write (5101,1) 0,    (i, i=1,npts)
    write (5102,2) 0,0,0,(i, i=1,fpts)
    !
    write (5101,6) ("-",i=1,npts*7+9)
    write (5102,6) ("-",i=1,fpts*7+17)
    !
    do nc = 1,lcell
      !
      nodes = grid%elem(nc)%nodes(:) ! F2003 AUTO-REALLOC
      pts   = grid%elem(nc)%pts(:)   ! F2003 AUTO-REALLOC
      !
      write (4101,1) nc,(nodes(i), i=1,size(nodes))
      write (5101,1) nc,(pts(i), i=1,size(pts))
      !
      write (4103,4) nc
      do i = 1,size(nodes)
        l = nodes(i)
        x = grid%xyz(1,l)
        y = grid%xyz(2,l)
        z = grid%xyz(3,l)
        r = sqrt(y*y+z*z)
        t = atan2(z,y)
        t = t + merge(two*pi,zero,t<zero)
        write (4103,3) l,x,y,z,x,r,t*radians2degrees
      end do
      !
      write (5103,4) nc
      do i = 1,size(pts)
        l = pts(i)
        x = grid%xyz(1,l)
        y = grid%xyz(2,l)
        z = grid%xyz(3,l)
        r = sqrt(y*y+z*z)
        t = atan2(z,y)
        t = t + merge(two*pi,zero,t<zero)
        write (5103,3) l,x,y,z,x,r,t*radians2degrees
      end do
      !
    end do
    !
    do nf = 1,lfbnd
      !
      nc = lcell + nf
      !
      nodes = grid%elem(nc)%nodes(:) ! F2003 AUTO-REALLOC
      pts   = grid%elem(nc)%pts(:)   ! F2003 AUTO-REALLOC
      !
      write (4102,2) nf,bface(1,nf),bface(4,nf),(nodes(i), i=1,size(nodes))
      write (5102,2) nf,bface(1,nf),bface(4,nf),(pts(i), i=1,size(pts))
      !
      write (4104,5) nf
      do i = 1,size(nodes)
        l = nodes(i)
        x = grid%xyz(1,l)
        y = grid%xyz(2,l)
        z = grid%xyz(3,l)
        r = sqrt(y*y+z*z)
        t = atan2(z,y)
        t = t + merge(two*pi,zero,t<zero)
        write (4104,3) l,x,y,z,x,r,t*radians2degrees
      end do
      !
      write (5104,5) nf
      do i = 1,size(pts)
        l = pts(i)
        x = grid%xyz(1,l)
        y = grid%xyz(2,l)
        z = grid%xyz(3,l)
        r = sqrt(y*y+z*z)
        t = atan2(z,y)
        t = t + merge(two*pi,zero,t<zero)
        write (5104,3) l,x,y,z,x,r,t*radians2degrees
      end do
      !
    end do
    !
  end if
  !
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine write_internal_grid
!
!###############################################################################
!
subroutine create_periodic_square_quad_grid()
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
  character(len=*), parameter :: pname = "create_periodic_square_quad_grid"
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
  ! Reset bc_input for this fully periodic grid so that the subroutine
  ! match_periodic_faces can correctly match all the periodic faces.
  !
  call reset_bc_input_for_internal_grid
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
  xyz_nodes = init_array(nr,nnode,zero) ! F2003 AUTO-REALLOC
  !
  nodes_of_cell_ptr = init_array(ncell+1,0) ! F2003 AUTO-REALLOC
  !
  nodes_of_cell = init_array(geom_nodes(Geom_Quad)*ncell,0) ! F03 AUTO-REALLOC
  !
  bface = init_array(nbfai,nfbnd,0) ! F2003 AUTO-REALLOC
  !
  cell_geom = init_array(ncell,Geom_Quad) ! F2003 AUTO-REALLOC
  !
  cell_order = init_array(ncell,n_order) ! F2003 AUTO-REALLOC
  !
  perturb = init_array(nr,nnode,real(0,kind=r4)) ! F2003 AUTO-REALLOC
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
  ! NOTE: Scaling the grid to (-Lref,-Lref) to (+Lref,+Lref) is now done
  !       in subroutine rescale_grid_for_vortex_problems.
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
      ! NOTE: The subroutine match_periodic_faces will set the values of
      !       bface([5,6,8,9],:), so we dont need to do that here.
      !
      ! face on x-min boundary
      !
      bface( 1,nf_min) = bc_periodic
      bface( 2,nf_min) = map2dto1d(1 ,j,nc)
      bface( 3,nf_min) = Geom_Quad
      bface( 4,nf_min) = HostFace2D(1)
      bface( 6,nf_min) = 1 ! index for bc_input needed by match_periodic_faces
      bface(10,nf_min) = ncell + nf_min
      !
      ! face on x-max boundary
      !
      bface( 1,nf_max) = bc_periodic
      bface( 2,nf_max) = map2dto1d(nc,j,nc)
      bface( 3,nf_max) = Geom_Quad
      bface( 4,nf_max) = HostFace2D(2)
      bface( 6,nf_max) = 2 ! index for bc_input needed by match_periodic_faces
      bface(10,nf_max) = ncell + nf_max
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
      ! NOTE: The subroutine match_periodic_faces will set the values of
      !       bface([5,6,8,9],:), so we dont need to do that here.
      !
      ! face on y-min boundary
      !
      bface( 1,nf_min) = bc_periodic
      bface( 2,nf_min) = map2dto1d(i,1 ,nc)
      bface( 3,nf_min) = Geom_Quad
      bface( 4,nf_min) = HostFace2D(3)
      bface( 6,nf_min) = 3 ! index for bc_input needed by match_periodic_faces
      bface(10,nf_min) = ncell + nf_min
      !
      ! face on y-max boundary
      !
      bface( 1,nf_max) = bc_periodic
      bface( 2,nf_max) = map2dto1d(i,nc,nc)
      bface( 3,nf_max) = Geom_Quad
      bface( 4,nf_max) = HostFace2D(4)
      bface( 6,nf_max) = 4 ! index for bc_input needed by match_periodic_faces
      bface(10,nf_max) = ncell + nf_max
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
end subroutine create_periodic_square_quad_grid
!
!###############################################################################
!
subroutine create_periodic_cube_hex_grid()
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
  integer :: nf_min,nf_max,n0,nc_min,nc_max
  integer :: imin,imax,imin_side,imax_side
  integer :: jmin,jmax,jmin_side,jmax_side
  integer :: kmin,kmax,kmin_side,kmax_side
  real(r4) :: x,y,z,dl
  real(r4) :: rbas,radj,rval
  !
  !.. Local Arrays ..
  integer :: ip(1:8)
  integer :: ifp(1:4)
  integer :: icp(1:8)
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: ipelem(:,:)
  real(r4), allocatable :: perturb(:,:)
  !
  logical(lk), allocatable :: msk(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_periodic_cube_hex_grid"
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
  ! Reset bc_input for this fully periodic grid so that the subroutine
  ! match_periodic_faces can correctly match all the periodic faces.
  !
  call reset_bc_input_for_internal_grid
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
  xyz_nodes = init_array(nr,nnode,zero) ! F2003 AUTO-REALLOC
  !
  ipelem = init_array(geom_nodes(Geom_Hexa),ncell,0) ! F2003 AUTO-REALLOC
  !
  nodes_of_cell_ptr = init_array(ncell+1,0) ! F2003 AUTO-REALLOC
  !
  nodes_of_cell = init_array(geom_nodes(Geom_Hexa)*ncell,0) ! F03 AUTO-REALLOC
  !
  bface = init_array(nbfai,nfbnd,0) ! F2003 AUTO-REALLOC
  !
  cell_geom = init_array(ncell,Geom_Hexa) ! F2003 AUTO-REALLOC
  !
  cell_order = init_array(ncell,n_order) ! F2003 AUTO-REALLOC
  !
  perturb = init_array(nr,nnode,real(0,kind=r4)) ! F2003 AUTO-REALLOC
  !
  ! #######
  ! #######
  ! ## 1 ##
  ! #######
  ! #######
  !
  ! Create ipelem
  !
  do k = 1,nc
    do j = 1,nc
      do i = 1,nc
        !
        n = map3dto1d(i,j,k,nc)
        !
        ip(1) = map3dto1d(i  ,j  ,k  ,np)
        ip(2) = map3dto1d(i+1,j  ,k  ,np)
        ip(3) = map3dto1d(i+1,j+1,k  ,np)
        ip(4) = map3dto1d(i  ,j+1,k  ,np)
        ip(5) = map3dto1d(i  ,j  ,k+1,np)
        ip(6) = map3dto1d(i+1,j  ,k+1,np)
        ip(7) = map3dto1d(i+1,j+1,k+1,np)
        ip(8) = map3dto1d(i  ,j+1,k+1,np)
        !
        if (randomly_rotate_cells_debug) then
          !
          ! Randomly rotate the cell connectivity of all the hex cells
          ! NOTE: This leaves the grid nodes in a structured ijk order.
          !       Only the cell connectivities are randomized.
          !       The cells will also still be in a structured ijk order,
          !       but what faces of each cell face in each ijk direction
          !       will be randomized.
          !
          l = min( 24 , max(1,random(25)) )
          !
          ipelem(1,n) = ip( hex_random_rotation(1,l) )
          ipelem(2,n) = ip( hex_random_rotation(2,l) )
          ipelem(3,n) = ip( hex_random_rotation(3,l) )
          ipelem(4,n) = ip( hex_random_rotation(4,l) )
          ipelem(5,n) = ip( hex_random_rotation(5,l) )
          ipelem(6,n) = ip( hex_random_rotation(6,l) )
          ipelem(7,n) = ip( hex_random_rotation(7,l) )
          ipelem(8,n) = ip( hex_random_rotation(8,l) )
          !
        else
          !
          ipelem(1,n) = ip(1)
          ipelem(2,n) = ip(2)
          ipelem(3,n) = ip(3)
          ipelem(4,n) = ip(4)
          ipelem(5,n) = ip(5)
          ipelem(6,n) = ip(6)
          ipelem(7,n) = ip(7)
          ipelem(8,n) = ip(8)
          !
        end if
        !
      end do
    end do
  end do
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
  do n = 1,ncell
    !
    n0 = nodes_of_cell_ptr(n)
    !
    nodes_of_cell(n0+1) = ipelem(1,n)
    nodes_of_cell(n0+2) = ipelem(2,n)
    nodes_of_cell(n0+3) = ipelem(3,n)
    nodes_of_cell(n0+4) = ipelem(4,n)
    nodes_of_cell(n0+5) = ipelem(5,n)
    nodes_of_cell(n0+6) = ipelem(6,n)
    nodes_of_cell(n0+7) = ipelem(7,n)
    nodes_of_cell(n0+8) = ipelem(8,n)
    !
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
  ! NOTE: Scaling the grid to (-Lref,-Lref,-Lref) to (+Lref,+Lref,+Lref)
  !       is now done in subroutine rescale_grid_for_vortex_problems.
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
  ! #######
  ! #######
  ! ## 3 ##
  ! #######
  ! #######
  !
  msk = init_array(nnode,fals) ! F2003 AUTO-REALLOC
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
      imin = 1
      imax = nc
      !
      nc_min = map3dto1d(imin,j,k,nc)
      nc_max = map3dto1d(imax,j,k,nc)
      !
      n1 = n1 + 1 ! x-min face
      n2 = n2 + 1 ! x-max face
      !
      if (.not. randomly_rotate_cells_debug) then
        !
        imin_side = HostFace3D(1)
        imax_side = HostFace3D(2)
        !
      else
        !
        ! get the structured ijk nodes of the x-min face
        !
        ip(1) = map3dto1d(imin  ,j  ,k  ,np)
        ip(2) = map3dto1d(imin+1,j  ,k  ,np)
        ip(3) = map3dto1d(imin+1,j+1,k  ,np)
        ip(4) = map3dto1d(imin  ,j+1,k  ,np)
        ip(5) = map3dto1d(imin  ,j  ,k+1,np)
        ip(6) = map3dto1d(imin+1,j  ,k+1,np)
        ip(7) = map3dto1d(imin+1,j+1,k+1,np)
        ip(8) = map3dto1d(imin  ,j+1,k+1,np)
        !
        ! get the face nodes that are on this x-min boundary face
        !
        ifp(1:4) = ip( side_nodes(1:4,HostFace3D(1),Geom_Hexa) + 1 )
        !
        ! mark these face nodes
        !
        msk(ifp(1:4)) = true
        !
        ! get the cell nodes for this cell in the randomized order
        !
        icp(1:8) = ipelem(1:8,nc_min)
        !
        ! loop through the different cell faces and find which face
        ! has the same nodes as the x-min boundary face
        !
        imin_side = -99
        imin_side_loop: &
        do l = 1,6
          if (count(msk(icp(side_nodes(1:4,l,Geom_Hexa)+1))) == 4) then
            imin_side = l
            exit imin_side_loop
          end if
        end do imin_side_loop
        !
        ! reset the marker for these face nodes
        !
        msk(ifp(1:4)) = fals
        !
        ! get the structured ijk nodes of the x-max face
        !
        ip(1) = map3dto1d(imax  ,j  ,k  ,np)
        ip(2) = map3dto1d(imax+1,j  ,k  ,np)
        ip(3) = map3dto1d(imax+1,j+1,k  ,np)
        ip(4) = map3dto1d(imax  ,j+1,k  ,np)
        ip(5) = map3dto1d(imax  ,j  ,k+1,np)
        ip(6) = map3dto1d(imax+1,j  ,k+1,np)
        ip(7) = map3dto1d(imax+1,j+1,k+1,np)
        ip(8) = map3dto1d(imax  ,j+1,k+1,np)
        !
        ifp(1:4) = ip( side_nodes(1:4,HostFace3D(2),Geom_Hexa) + 1 )
        !
        msk(ifp(1:4)) = true
        !
        icp(1:8) = ipelem(1:8,nc_max)
        imax_side = -99
        imax_side_loop: &
        do l = 1,6
          if (count(msk(icp(side_nodes(1:4,l,Geom_Hexa)+1))) == 4) then
            imax_side = l
            exit imax_side_loop
          end if
        end do imax_side_loop
        !
        msk(ifp(1:4)) = fals
        !
      end if
      !
      ! NOTE: The subroutine match_periodic_faces will set the values of
      !       bface([5,6,8,9],:), so we dont need to do that here.
      !
      bface( 1,n1) = bc_periodic
      bface( 2,n1) = nc_min
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = imin_side
      bface( 6,n1) = 1 ! index for bc_input needed by match_periodic_faces
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_periodic
      bface( 2,n2) = nc_max
      bface( 3,n2) = Geom_Hexa
      bface( 4,n2) = imax_side
      bface( 6,n2) = 2 ! index for bc_input needed by match_periodic_faces
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
      jmin = 1
      jmax = nc
      !
      nc_min = map3dto1d(i,jmin,k,nc)
      nc_max = map3dto1d(i,jmax,k,nc)
      !
      n1 = n1 + 1 ! y-min face
      n2 = n2 + 1 ! y-max face
      !
      if (.not. randomly_rotate_cells_debug) then
        !
        jmin_side = HostFace3D(3)
        jmax_side = HostFace3D(4)
        !
      else
        !
        ! get the structured ijk nodes of the y-min face
        !
        ip(1) = map3dto1d(i  ,jmin  ,k  ,np)
        ip(2) = map3dto1d(i+1,jmin  ,k  ,np)
        ip(3) = map3dto1d(i+1,jmin+1,k  ,np)
        ip(4) = map3dto1d(i  ,jmin+1,k  ,np)
        ip(5) = map3dto1d(i  ,jmin  ,k+1,np)
        ip(6) = map3dto1d(i+1,jmin  ,k+1,np)
        ip(7) = map3dto1d(i+1,jmin+1,k+1,np)
        ip(8) = map3dto1d(i  ,jmin+1,k+1,np)
        !
        ifp(1:4) = ip( side_nodes(1:4,HostFace3D(3),Geom_Hexa) + 1 )
        !
        msk(ifp(1:4)) = true
        !
        icp(1:8) = ipelem(1:8,nc_min)
        jmin_side = -99
        jmin_side_loop: &
        do l = 1,6
          if (count(msk(icp(side_nodes(1:4,l,Geom_Hexa)+1))) == 4) then
            jmin_side = l
            exit jmin_side_loop
          end if
        end do jmin_side_loop
        !
        msk(ifp(1:4)) = fals
        !
        ! get the structured ijk nodes of the y-max face
        !
        ip(1) = map3dto1d(i  ,jmax  ,k  ,np)
        ip(2) = map3dto1d(i+1,jmax  ,k  ,np)
        ip(3) = map3dto1d(i+1,jmax+1,k  ,np)
        ip(4) = map3dto1d(i  ,jmax+1,k  ,np)
        ip(5) = map3dto1d(i  ,jmax  ,k+1,np)
        ip(6) = map3dto1d(i+1,jmax  ,k+1,np)
        ip(7) = map3dto1d(i+1,jmax+1,k+1,np)
        ip(8) = map3dto1d(i  ,jmax+1,k+1,np)
        !
        ifp(1:4) = ip( side_nodes(1:4,HostFace3D(4),Geom_Hexa) + 1 )
        !
        msk(ifp(1:4)) = true
        !
        icp(1:8) = ipelem(1:8,nc_max)
        jmax_side = -99
        jmax_side_loop: &
        do l = 1,6
          if (count(msk(icp(side_nodes(1:4,l,Geom_Hexa)+1))) == 4) then
            jmax_side = l
            exit jmax_side_loop
          end if
        end do jmax_side_loop
        !
        msk(ifp(1:4)) = fals
        !
      end if
      !
      ! NOTE: The subroutine match_periodic_faces will set the values of
      !       bface([5,6,8,9],:), so we dont need to do that here.
      !
      bface( 1,n1) = bc_periodic
      bface( 2,n1) = nc_min
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = jmin_side
      bface( 6,n1) = 3 ! index for bc_input needed by match_periodic_faces
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_periodic
      bface( 2,n2) = nc_max
      bface( 3,n2) = Geom_Hexa
      bface( 4,n2) = jmax_side
      bface( 6,n2) = 4 ! index for bc_input needed by match_periodic_faces
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
      kmin = 1
      kmax = nc
      !
      nc_min = map3dto1d(i,j,kmin,nc)
      nc_max = map3dto1d(i,j,kmax,nc)
      !
      n1 = n1 + 1 ! z-min face
      n2 = n2 + 1 ! z-max face
      !
      if (.not. randomly_rotate_cells_debug) then
        !
        kmin_side = HostFace3D(5)
        kmax_side = HostFace3D(6)
        !
      else
        !
        ! get the structured ijk nodes of the y-min face
        !
        ip(1) = map3dto1d(i  ,j  ,kmin  ,np)
        ip(2) = map3dto1d(i+1,j  ,kmin  ,np)
        ip(3) = map3dto1d(i+1,j+1,kmin  ,np)
        ip(4) = map3dto1d(i  ,j+1,kmin  ,np)
        ip(5) = map3dto1d(i  ,j  ,kmin+1,np)
        ip(6) = map3dto1d(i+1,j  ,kmin+1,np)
        ip(7) = map3dto1d(i+1,j+1,kmin+1,np)
        ip(8) = map3dto1d(i  ,j+1,kmin+1,np)
        !
        ifp(1:4) = ip( side_nodes(1:4,HostFace3D(5),Geom_Hexa) + 1 )
        !
        msk(ifp(1:4)) = true
        !
        icp(1:8) = ipelem(1:8,nc_min)
        kmin_side = -99
        kmin_side_loop: &
        do l = 1,6
          if (count(msk(icp(side_nodes(1:4,l,Geom_Hexa)+1))) == 4) then
            kmin_side = l
            exit kmin_side_loop
          end if
        end do kmin_side_loop
        !
        msk(ifp(1:4)) = fals
        !
        ! get the structured ijk nodes of the y-max face
        !
        ip(1) = map3dto1d(i  ,j  ,kmax  ,np)
        ip(2) = map3dto1d(i+1,j  ,kmax  ,np)
        ip(3) = map3dto1d(i+1,j+1,kmax  ,np)
        ip(4) = map3dto1d(i  ,j+1,kmax  ,np)
        ip(5) = map3dto1d(i  ,j  ,kmax+1,np)
        ip(6) = map3dto1d(i+1,j  ,kmax+1,np)
        ip(7) = map3dto1d(i+1,j+1,kmax+1,np)
        ip(8) = map3dto1d(i  ,j+1,kmax+1,np)
        !
        ifp(1:4) = ip( side_nodes(1:4,HostFace3D(6),Geom_Hexa) + 1 )
        !
        msk(ifp(1:4)) = true
        !
        icp(1:8) = ipelem(1:8,nc_max)
        kmax_side = -99
        kmax_side_loop: &
        do l = 1,6
          if (count(msk(icp(side_nodes(1:4,l,Geom_Hexa)+1))) == 4) then
            kmax_side = l
            exit kmax_side_loop
          end if
        end do kmax_side_loop
        !
        msk(ifp(1:4)) = fals
        !
      end if
      !
      ! NOTE: The subroutine match_periodic_faces will set the values of
      !       bface([5,6,8,9],:), so we dont need to do that here.
      !
      bface( 1,n1) = bc_periodic
      bface( 2,n1) = nc_min
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = kmin_side
      bface( 6,n1) = 5 ! index for bc_input needed by match_periodic_faces
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_periodic
      bface( 2,n2) = nc_max
      bface( 3,n2) = Geom_Hexa
      bface( 4,n2) = kmax_side
      bface( 6,n2) = 6 ! index for bc_input needed by match_periodic_faces
      bface(10,n2) = ncell + n2
      !
    end do
  end do
  !
  ! Deallocate ipelem and msk since they are no longer needed
  !
  deallocate ( ipelem , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ipelem",2,__LINE__,__FILE__,ierr,error_message)
  !
  deallocate ( msk , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
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
end subroutine create_periodic_cube_hex_grid
!
!###############################################################################
!
subroutine create_periodic_cube_pyramid_grid()
  !
  !.. Use Statements ..
  use order_mod, only : n_order,geom_solpts
  use ovar,      only : VortexGrid_DOF,Lref
  use ovar,      only : VortexGrid_random_perturb
  use ovar,      only : VortexGrid_perturb_percent
  use geovar,    only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geovar,    only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geovar,    only : cell_geom,cell_order
  use ovar,      only : grid_format
  !
  !.. Local Scalars ..
  integer :: nc,np,ierr,i,j,k,l,n,n1,n2,n0
  integer :: ncf,ncp,nc_min,nc_max,nc_hex
  integer :: imin,imax,jmin,jmax,kmin,kmax,ip_mid
  integer :: npyr_dof_per_hex,n_hex_cells,n_hex_nodes
  integer :: nc_hex_min,nc_hex_max
  real(wp) :: x,y,z,dl
  real(r4) :: rbas,radj,rval
  !
  !.. Local Arrays ..
  integer :: ip_hex(1:8)
  integer :: ip_pyr(1:5,1:6)
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: ipelem(:,:)
  integer,  allocatable :: ipehex(:,:)
  real(r4), allocatable :: perturb(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_periodic_cube_pyramid_grid"
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
  ! Reset bc_input for this fully periodic grid so that the subroutine
  ! match_periodic_faces can correctly match all the periodic faces.
  !
  call reset_bc_input_for_internal_grid
  !
  ! Estimate the number of cells required to get the requested DOF. We are
  ! converting each hex cell into 6 pyramid cells by adding a center node
  ! to the hex cell; therefore, the total number of pyramid DoF per hex cell
  ! is 6*(# DoF per pyramid cell). We need the number of pyramid DoF per hex
  ! cell in each direction, so we need to take the cube root of the total
  ! number of pyramid DoF per hex cell.
  !
  ncp = geom_solpts(Geom_Pyra,n_order)
  npyr_dof_per_hex = nint( real(max(1,6*ncp),kind=wp)**one3 )
  !
  ! Get the number of cells in each direction
  !
  nc = nint( real(VortexGrid_DOF,kind=wp) / real(npyr_dof_per_hex,kind=wp) )
  !
  np = nc + 1 ! number of grid nodes in each direction
  !
  n_hex_cells = nc*nc*nc ! total number of hex grid cells
  n_hex_nodes = np*np*np ! total number of hex grid nodes
  !
  ncell = 6*n_hex_cells ! total number of pyramid grid cells
  nnode = n_hex_nodes + n_hex_cells ! total number of pyramid grid nodes
  nfbnd = 6*nc*nc  ! total number of boundary faces (same as for hex grid)
  !
  ! Allocate all the arrays needed to define the grid
  !
  xyz_nodes = init_array(nr,nnode,zero) ! F2003 AUTO-REALLOC
  !
  ipelem = init_array(geom_nodes(Geom_Pyra),ncell,0) ! F2003 AUTO-REALLOC
  ipehex = init_array(geom_nodes(Geom_Hexa),n_hex_cells,0) ! F2003 AUTO-REALLOC
  !
  nodes_of_cell_ptr = init_array(ncell+1,0) ! F2003 AUTO-REALLOC
  !
  nodes_of_cell = init_array(geom_nodes(Geom_Pyra)*ncell,0) ! F03 AUTO-REALLOC
  !
  bface = init_array(nbfai,nfbnd,0) ! F2003 AUTO-REALLOC
  !
  cell_geom = init_array(ncell,Geom_Pyra) ! F2003 AUTO-REALLOC
  !
  cell_order = init_array(ncell,n_order) ! F2003 AUTO-REALLOC
  !
  perturb = init_array(nr,nnode,real(0,kind=r4)) ! F2003 AUTO-REALLOC
  !
  ! #######
  ! #######
  ! ## 1 ##
  ! #######
  ! #######
  !
  ! Create ipelem
  !
  do k = 1,nc
    do j = 1,nc
      do i = 1,nc
        !
        nc_hex = map3dto1d(i,j,k,nc)
        !
        n0 = 6*(nc_hex-1)
        !
        ip_hex(1) = map3dto1d(i  ,j  ,k  ,np)
        ip_hex(2) = map3dto1d(i+1,j  ,k  ,np)
        ip_hex(3) = map3dto1d(i+1,j+1,k  ,np)
        ip_hex(4) = map3dto1d(i  ,j+1,k  ,np)
        ip_hex(5) = map3dto1d(i  ,j  ,k+1,np)
        ip_hex(6) = map3dto1d(i+1,j  ,k+1,np)
        ip_hex(7) = map3dto1d(i+1,j+1,k+1,np)
        ip_hex(8) = map3dto1d(i  ,j+1,k+1,np)
        !
        ipehex(1:8,nc_hex) = ip_hex(1:8)
        !
        ip_mid = n_hex_nodes + nc_hex
        !
        ! Get the cell connectivities of the 6 pyramids in this hex cel
        !
        ip_pyr(1:5,1) = [ ip_hex([1,4,8,5]) , ip_mid ] ! x-min hex face as base
        ip_pyr(1:5,2) = [ ip_hex([2,6,7,3]) , ip_mid ] ! x-max hex face as base
        ip_pyr(1:5,3) = [ ip_hex([1,5,6,2]) , ip_mid ] ! y-min hex face as base
        ip_pyr(1:5,4) = [ ip_hex([4,3,7,8]) , ip_mid ] ! y-max hex face as base
        ip_pyr(1:5,5) = [ ip_hex([1,2,3,4]) , ip_mid ] ! z-min hex face as base
        ip_pyr(1:5,6) = [ ip_hex([5,8,7,6]) , ip_mid ] ! z-max hex face as base
        !
        ! Randomly rotate the pyramid cell connectivities
        !
        if (randomly_rotate_cells_debug) then
          do ncf = 1,6
            l = min( 4 , max(1,random(5)) )
            ip_pyr(1:5,ncf) = ip_pyr( pyr_random_rotation(1:5,l) , ncf )
          end do
        end if
        !
        ! Copy the pyramid connectivities into ipelem
        !
       !do ncf = 1,6
       !  do l = 1,5
       !    ipelem(l,n0+ncf) = ip_pyr(l,ncf)
       !  end do
       !end do
        !
        ipelem(1:5,n0+1) = ip_pyr(1:5,1)
        ipelem(1:5,n0+2) = ip_pyr(1:5,2)
        ipelem(1:5,n0+3) = ip_pyr(1:5,3)
        ipelem(1:5,n0+4) = ip_pyr(1:5,4)
        ipelem(1:5,n0+5) = ip_pyr(1:5,5)
        ipelem(1:5,n0+6) = ip_pyr(1:5,6)
        !
      end do
    end do
  end do
  !
 !do n = 1,n_hex_cells
 !  write (8502,851) n,(ipehex(i,n), i=1,8)
 !end do
 !851 format (i5,3x,8i4)
 !do n = 1,ncell
 !  write (8501,851) n,(ipelem(i,n), i=1,5)
 !end do
  !
  ! Create the nodes_of_cell_ptr array
  !
  nodes_of_cell_ptr(1) = 0
  do n = 2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(n) = nodes_of_cell_ptr(n-1) + geom_nodes(Geom_Pyra)
  end do
  !
  ! Create the nodes_of_cell array
  !
  do n = 1,ncell
    !
    n0 = nodes_of_cell_ptr(n)
    !
    nodes_of_cell(n0+1) = ipelem(1,n)
    nodes_of_cell(n0+2) = ipelem(2,n)
    nodes_of_cell(n0+3) = ipelem(3,n)
    nodes_of_cell(n0+4) = ipelem(4,n)
    nodes_of_cell(n0+5) = ipelem(5,n)
    !
  end do
  !
  ! Deallocate ipelem since it is no longer needed
  !
  deallocate ( ipelem , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ipelem",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! #######
  ! #######
  ! ## 2 ##
  ! #######
  ! #######
  !
  ! Compute the cell length for a Cartesian/uniform grid
  !
  dl = real(1,kind=wp) / real(np-1,kind=wp)
  !
  ! Compute the random perturbations if the flag was set on input
  !
  if (VortexGrid_random_perturb) then
    !
    rbas = real(dl,kind=r4) * real(abs(VortexGrid_perturb_percent),kind=r4) &
              / real(100,kind=r4)
    radj = rbas / real(2,kind=r4)
    !
    ! Randomly perturb the crap out of the hex nodes
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
    ! Randomly perturb the crap out of the center nodes in each hex
    ! cell that defines the apex of the 6 pyramids
    !
    do k = 2,nc-1
      do j = 2,nc-1
        do i = 2,nc-1
          !
          n = n_hex_nodes + map3dto1d(i,j,k,nc)
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
  ! NOTE: Scaling the grid to (-Lref,-Lref,-Lref) to (+Lref,+Lref,+Lref)
  !       is now done in subroutine rescale_grid_for_vortex_problems.
  !
  do k = 1,np
    !
    z = dl * real(k-1,kind=wp)
    !
    do j = 1,np
      !
      y = dl * real(j-1,kind=wp)
      !
      do i = 1,np
        !
        x = dl * real(i-1,kind=wp)
        !
        n = map3dto1d(i,j,k,np)
        !
        xyz_nodes(1,n) = x
        xyz_nodes(2,n) = y
        xyz_nodes(3,n) = z
       !write (8503,853) n,(xyz_nodes(l,n),l=1,3)
       !853 format (i4,3x,3f8.4,:,3x,i4)
        !
      end do
    end do
  end do
  !
  ! Create the grid node coordinates for the center nodes in each hex
  ! cell that defines the apex of the 6 pyramids
  !
  do k = 1,nc
    !
    z = dl * (real(k,kind=wp) - half)
    !
    do j = 1,nc
      !
      y = dl * (real(j,kind=wp) - half)
      !
      do i = 1,nc
        !
        x = dl * (real(i,kind=wp) - half)
        !
        nc_hex = map3dto1d(i,j,k,nc)
        n = n_hex_nodes + nc_hex
        !
        xyz_nodes(1,n) = x
        xyz_nodes(2,n) = y
        xyz_nodes(3,n) = z
       !write (8504,853) n,(xyz_nodes(l,n),l=1,3),nc_hex
        !
      end do
      !
    end do
    !
  end do
  !
  ! Now randomly perturb the grid nodes
  !
  if (VortexGrid_random_perturb) then
    !
    do n = 1,nnode
      xyz_nodes(1,n) = xyz_nodes(1,n) + real(perturb(1,n),kind=wp)
      xyz_nodes(2,n) = xyz_nodes(2,n) + real(perturb(2,n),kind=wp)
      xyz_nodes(3,n) = xyz_nodes(3,n) + real(perturb(3,n),kind=wp)
    end do
    !
  end if
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
  ! Create the bface array
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
      imin = 1
      imax = nc
      !
      nc_hex_min = map3dto1d(imin,j,k,nc)
      nc_hex_max = map3dto1d(imax,j,k,nc)
      !
      nc_min = 6*(nc_hex_min-1) + 1 ! pyramid with x-min hex face as base
      nc_max = 6*(nc_hex_max-1) + 2 ! pyramid with x-max hex face as base
      !
      n1 = n1 + 1 ! x-min face
      n2 = n2 + 1 ! x-max face
      !
      ! NOTE: The subroutine match_periodic_faces will set the values of
      !       bface([5,6,8,9],:), so we dont need to do that here.
      !
      bface( 1,n1) = bc_periodic
      bface( 2,n1) = nc_min
      bface( 3,n1) = Geom_Pyra
      bface( 4,n1) = 5 ! always 5 since quad faces on periodic boundaries
      bface( 6,n1) = 1 ! index for bc_input needed by match_periodic_faces
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_periodic
      bface( 2,n2) = nc_max
      bface( 3,n2) = Geom_Pyra
      bface( 4,n2) = 5 ! always 5 since quad faces on periodic boundaries
      bface( 6,n2) = 2 ! index for bc_input needed by match_periodic_faces
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
      jmin = 1
      jmax = nc
      !
      nc_hex_min = map3dto1d(i,jmin,k,nc)
      nc_hex_max = map3dto1d(i,jmax,k,nc)
      !
      nc_min = 6*(nc_hex_min-1) + 3 ! pyramid with y-min hex face as base
      nc_max = 6*(nc_hex_max-1) + 4 ! pyramid with y-max hex face as base
      !
      n1 = n1 + 1 ! y-min face
      n2 = n2 + 1 ! y-max face
      !
      ! NOTE: The subroutine match_periodic_faces will set the values of
      !       bface([5,6,8,9],:), so we dont need to do that here.
      !
      bface( 1,n1) = bc_periodic
      bface( 2,n1) = nc_min
      bface( 3,n1) = Geom_Pyra
      bface( 4,n1) = 5 ! always 5 since quad faces on periodic boundaries
      bface( 6,n1) = 3 ! index for bc_input needed by match_periodic_faces
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_periodic
      bface( 2,n2) = nc_max
      bface( 3,n2) = Geom_Pyra
      bface( 4,n2) = 5 ! always 5 since quad faces on periodic boundaries
      bface( 6,n2) = 4 ! index for bc_input needed by match_periodic_faces
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
      kmin = 1
      kmax = nc
      !
      nc_hex_min = map3dto1d(i,j,kmin,nc)
      nc_hex_max = map3dto1d(i,j,kmax,nc)
      !
      nc_min = 6*(nc_hex_min-1) + 5 ! pyramid with z-min hex face as base
      nc_max = 6*(nc_hex_max-1) + 6 ! pyramid with z-max hex face as base
      !
      n1 = n1 + 1 ! z-min face
      n2 = n2 + 1 ! z-max face
      !
      ! NOTE: The subroutine match_periodic_faces will set the values of
      !       bface([5,6,8,9],:), so we dont need to do that here.
      !
      bface( 1,n1) = bc_periodic
      bface( 2,n1) = nc_min
      bface( 3,n1) = Geom_Pyra
      bface( 4,n1) = 5 ! always 5 since quad faces on periodic boundaries
      bface( 6,n1) = 5 ! index for bc_input needed by match_periodic_faces
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_periodic
      bface( 2,n2) = nc_max
      bface( 3,n2) = Geom_Pyra
      bface( 4,n2) = 5 ! always 5 since quad faces on periodic boundaries
      bface( 6,n2) = 6 ! index for bc_input needed by match_periodic_faces
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
end subroutine create_periodic_cube_pyramid_grid
!
!###############################################################################
!
subroutine reset_bc_input_for_internal_grid()
  !
  !.. Use Statements ..
  use ovar,   only : bc_input
  use geovar, only : nr
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reset_bc_input_for_internal_grid"
  !
continue
  !
  if (allocated(bc_input)) then
    deallocate ( bc_input , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bc_input",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( bc_input(1:2*nr) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bc_input",1,__LINE__,__FILE__,ierr,error_message)
  !
  bc_input(:)%set_bc_value = bc_periodic
  !
  bc_input(1)%name = "xmin"
  bc_input(1)%periodic%matching_bc_name = "xmax"
  bc_input(1)%periodic%matching_bc_idx  = 2
  bc_input(1)%periodic%translation(1)   = one
  !
  bc_input(2)%name = "xmax"
  bc_input(2)%periodic%matching_bc_name = "xmin"
  bc_input(2)%periodic%matching_bc_idx  = 1
  bc_input(2)%periodic%translation(1)   = -one
  !
  if (nr > 1) then
    bc_input(3)%name = "ymin"
    bc_input(3)%periodic%matching_bc_name = "ymax"
    bc_input(3)%periodic%matching_bc_idx  = 4
    bc_input(3)%periodic%translation(2)   = one
    !
    bc_input(4)%name = "ymax"
    bc_input(4)%periodic%matching_bc_name = "ymin"
    bc_input(4)%periodic%matching_bc_idx  = 3
    bc_input(4)%periodic%translation(2)   = -one
  end if
  !
  if (nr > 2) then
    bc_input(5)%name = "zmin"
    bc_input(5)%periodic%matching_bc_name = "zmax"
    bc_input(5)%periodic%matching_bc_idx  = 6
    bc_input(5)%periodic%translation(3)   = one
    !
    bc_input(6)%name = "zmax"
    bc_input(6)%periodic%matching_bc_name = "zmin"
    bc_input(6)%periodic%matching_bc_idx  = 5
    bc_input(6)%periodic%translation(3)   = -one
  end if
  !
end subroutine reset_bc_input_for_internal_grid
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
#ifndef DISABLE_PURE
pure &
#endif
function tanh_spacing_1d(x1in,x2in,ds1in,ds2in,npts,special_b_coeff) &
                  result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: x1in
  real(wp), intent(in) :: x2in
  real(wp), intent(in) :: ds1in
  real(wp), intent(in) :: ds2in
  integer,  intent(in) :: npts
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: special_b_coeff
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer  :: nl,np,l
  real(wp) :: dl,x
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: x1(:)
  real(wp), allocatable :: x2(:)
  real(wp), allocatable :: x1i(:)
  real(wp), allocatable :: x2i(:)
  real(wp), allocatable :: xp(:,:)
  !
continue
  !
  ! We are going to enforce the end spacings ds1in and ds2in,
  ! so we only want to solve for the npts-3 interior spacings
  ! (meaning npts-2 points on the line).
  !
  np = npts-2
  !
  ! Allocate the function result
  !
  return_value = init_array(npts,zero) ! F2003 AUTO-REALLOC
  !
  ! Set the end points of the function result
  ! to the end points input to the function
  !
  return_value(   1) = x1in
  return_value(npts) = x2in
  !
  ! Special cases for npts=2 (a line with just two end points)
  ! or npts=3 (a line with just one interior point)
  !
  if (np == 0) then
    return
  else if (np == 1) then
    return_value(2) = average_1p_spacing(x1in,x2in,ds1in,ds2in)
    return
  end if
  !
  ! nl is always 1 for this version of tanh_spacing
  !
  nl = 1
  !
  x1 = [x1in] ! F2003 AUTO-REALLOC
  x2 = [x2in] ! F2003 AUTO-REALLOC
  !
  ! Get the distance between the two points
  !
  dl = norm2( x1(:) - x2(:) )
  !
  ! Initialize the x1i and x2i arrays which will hold the coordinates
  ! for the first points off the ends, i.e., [x1,x1i,...,x2i,x2].
  !
  x1i = init_array(nl,zero) ! F2003 AUTO-REALLOC
  x2i = init_array(nl,zero) ! F2003 AUTO-REALLOC
  !
  ! Modify x1i and x2i to be the location for the 2nd and npts-1 points
  !
  do l = 1,nl
    x1i(l) = x1(l) + (x2(l)-x1(l))*ds1in/dl
    x2i(l) = x2(l) - (x2(l)-x1(l))*ds2in/dl
  end do
  !
  ! F2003 AUTO-REALLOC
  xp = get_tanh_spacing(x1i,x2i,ds1in,ds2in,np,special_b_coeff)
  !
  ! Copy the interior point spacings to the function result
  !
  return_value(2:np+1) = xp(1,1:np)
  !
end function tanh_spacing_1d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function average_1p_spacing_1d(x1,x2,ds1,ds2) result(return_value)
  !.. Formal Arguments ..
  real(wp), intent(in) :: x1
  real(wp), intent(in) :: x2
  real(wp), intent(in) :: ds1
  real(wp), intent(in) :: ds2
  !.. Function Result ..
  real(wp) :: return_value
continue
  return_value = x1 + (x2-x1)*ds1/(ds1+ds2)
end function average_1p_spacing_1d
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function average_1p_spacing_md(x1,x2,ds1,ds2) result(return_value)
  !.. Formal Arguments ..
  real(wp), intent(in) :: x1(:)
  real(wp), intent(in) :: x2(:)
  real(wp), intent(in) :: ds1
  real(wp), intent(in) :: ds2
  !.. Function Result ..
  real(wp) :: return_value(1:size(x1))
continue
  return_value(:) = x1(:) + (x2(:)-x1(:))*ds1/(ds1+ds2)
end function average_1p_spacing_md
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function tanh_spacing_md(x1in,x2in,ds1in,ds2in,npts,special_b_coeff) &
                  result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: x1in(:)
  real(wp), intent(in) :: x2in(:)
  real(wp), intent(in) :: ds1in
  real(wp), intent(in) :: ds2in
  integer,  intent(in) :: npts
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: special_b_coeff
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: nl,np,l
  real(wp) :: dl
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: x1(:)
  real(wp), allocatable :: x2(:)
  real(wp), allocatable :: x1i(:)
  real(wp), allocatable :: x2i(:)
  real(wp), allocatable :: xp(:,:)
  !
continue
  !
  ! We are going to enforce the end spacings ds1in and ds2in,
  ! so we only want to solve for the npts-3 interior spacings
  ! (meaning npts-2 points on the line).
  !
  np = npts-2
  !
  ! Make sure x1 and x2 have the same dimensions
  !
  nl = min(size(x1in),size(x2in))
  !
  x1 = x1in(1:nl) ! F2003 AUTO-REALLOC
  x2 = x2in(1:nl) ! F2003 AUTO-REALLOC
  !
  ! Allocate the function result
  !
  return_value = init_array(nl,npts,zero) ! F2003 AUTO-REALLOC
  !
  ! Set the end points of the function result
  ! to the end points input to the function
  !
  return_value(1:nl,   1) = x1(1:nl)
  return_value(1:nl,npts) = x2(1:nl)
  !
  ! Special cases for npts=2 (a line with just two end points)
  ! or npts=3 (a line with just one interior point)
  !
  if (np == 0) then
    return
  else if (np == 1) then
    x1i = average_1p_spacing(x1,x2,ds1in,ds2in) ! F2003 AUTO-REALLOC
    return_value(1:nl,2) = x1i(1:nl)
    return
  end if
  !
  ! Get the distance between the two points
  !
  dl = norm2( x1(:) - x2(:) )
  !
  ! Initialize the x1i and x2i arrays which will hold the coordinates
  ! for the first points off the ends, i.e., [x1,x1i,...,x2i,x2].
  !
  x1i = init_array(nl,zero) ! F2003 AUTO-REALLOC
  x2i = init_array(nl,zero) ! F2003 AUTO-REALLOC
  !
  ! Modify x1i and x2i to be the location for the 2nd and npts-1 points
  !
  do l = 1,nl
    x1i(l) = x1(l) + (x2(l)-x1(l))*ds1in/dl
    x2i(l) = x2(l) - (x2(l)-x1(l))*ds2in/dl
  end do
  !
  ! F2003 AUTO-REALLOC
  xp = get_tanh_spacing(x1i,x2i,ds1in,ds2in,np,special_b_coeff)
  !
  ! Copy the interior point spacings to the function result
  !
  return_value(1:nl,2:np+1) = xp(1:nl,1:np)
  !
end function tanh_spacing_md
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_tanh_spacing(x1in,x2in,ds1in,ds2in,np,special_b_coeff) &
                   result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: x1in(:)
  real(wp), intent(in) :: x2in(:)
  real(wp), intent(in) :: ds1in
  real(wp), intent(in) :: ds2in
  integer,  intent(in) :: np
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: special_b_coeff
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: n,nl,k,l
  real(wp) :: ds1,ds2,dl,dsl
  real(wp) :: a,b,sf,rnum,rden,sferr
  real(wp) :: val,ue,bval
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: x1(:)
  real(wp), allocatable :: x2(:)
  real(wp), allocatable :: ds(:)
  !
continue
  !
  ! bval is used as a coefficient when computing the b parameter for the
  ! tanh spacing. Using a coefficient of -1 seems to give the smoothest
  ! spacing with less dramatic stretching. The special value for b of +3.75
  ! is only for recreating the radial tanh spacings for a specific grid
  ! generated by Pointwise for the Taylor-Couette flow problem. This specific
  ! grid contained 20 radial cells, 48 azimuthal cells, and 24 axial cells,
  ! with an inner cylinder radius of 0.5 with wall spacing set to 0.002, and
  ! an outer cylinder radius of 1.0 with wall spacing set to 0.0025.
  !
  bval = -one
  if (present(special_b_coeff)) then
    bval = max( -one , min(special_b_coeff,3.75_wp) )
  end if
  !
  ! Make sure x1 and x2 have the same dimensions
  !
  nl = min(size(x1in),size(x2in))
  !
  x1 = x1in(1:nl) ! F2003 AUTO-REALLOC
  x2 = x2in(1:nl) ! F2003 AUTO-REALLOC
  !
  ! Initialize ds
  !
  ds = init_array(np,zero) ! F2003 AUTO-REALLOC
  !
  ! Get the distance between the two points
  !
  dl = norm2( x1(:) - x2(:) )
  !
  ! Map the beginning and ending grid cell
  ! spacings to the reference interval [0,1].
  !
  ds1 = ds1in/dl
  ds2 = ds2in/dl
  !
  ! 2C.2 : Get the a and b parameters to begin solving for the tanh y-spacing
  !
  a = sqrt(ds2/ds1)
 !b = one / (real(np-1,kind=wp)*sqrt(ds1*ds2))
  ! When enforcing the input spacings at the ends, this version for b
  ! is the closest match for the tanh spacing produced by Pointwise
  b = one / ( (real(np,kind=wp)+bval) * sqrt(ds1*ds2) )
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
    do k = 1,100
      !
      rnum = sinh(sf)/sf - b
      rden = cosh(sf)/sf - sinh(sf)/(sf*sf)
      !
      sferr = rnum / rden
      !
      if (abs(sferr) < eps10) then
        exit newton_tanh_loop
      end if
      !
      sf = sf - sferr
      !
      ! Make sure sf doesnt become 0
      !
      if (abs(sf) < eps2) then
        sf = eps2
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
    do n = 0,np-1
      val = real(n,kind=wp)/real(np-1,kind=wp)
      ue = half*(one + tanh(sf*(val-half))/tanh(sf*half))
      ds(n+1) = ue / (a + (one-a)*ue)
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
    do k = 1,100
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
      if (abs(sf) < eps2) then
        sf = eps2
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
    do n = 0,np-1
      val = real(n,kind=wp)/real(np-1,kind=wp)
      ue = half*(one + tan(sf*(val-half))/tan(sf*half))
      ds(n+1) = ue / (a + (one-a)*ue)
    end do
    !
  else
    !
    ! Both of the above Newton iterations can have numerical difficulties
    ! if the stretching factor is close to 1, therefore use this simpler
    ! stretching function that doesnt require a Newton iteration.
    !
    do n = 0,np-1
      val = real(n,kind=wp)/real(np-1,kind=wp)
      ue = val*(one + two*(b-one)*(val-half)*(one-val))
      ds(n+1) = ue / (a + (one-a)*ue)
    end do
    !
  end if
  !
  return_value = init_array(nl,np,zero) ! F2003 AUTO-REALLOC
  !
  return_value(1:nl, 1) = x1(1:nl)
  return_value(1:nl,np) = x2(1:nl)
  !
  do l = 1,nl
    !
    ! Get the length for this direction which is needed
    ! to unscale the stretching values from [0,1]
    !
    dsl = x2(l)-x1(l)
    !
    ! Finally loop through and use the scaled stretching values to get the
    ! coordinate for the points of this direction
    !
    do n = 2,np-1
      return_value(l,n) = x1(l) + ds(n) * dsl
    end do
    !
  end do
  !
end function get_tanh_spacing
!
!###############################################################################
!
end module internal_grid_mod
