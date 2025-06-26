module internal_grid_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  implicit none
  !
  private
  !
  !.. Module public procedures ..
  public :: create_internal_grid
  public :: write_internal_grid
  public :: write_channel_grid_info
  !
  !.. Module public scalars ..
  real(wp), public, save :: xplus_cell_ave
  real(wp), public, save :: yplus_cell_ave
  real(wp), public, save :: yplus_cell_wall
  real(wp), public, save :: yplus_cell_centerline
  real(wp), public, save :: zplus_cell_ave
  !
  !.. Module public allocatable arrays ..
  real(wp), public, save, allocatable :: yplus_dof_wall(:)
  !
  !.. Module Parameters ..
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
  use ovar, only : flame_speed_test_dof, sod_test_dof, pulse_test_dof, &
                   shu_osher_test_dof
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
  else if (any(abs(itestcase) == [Taylor_Couette_Flow,Freestream_TCF])) then
    !
    call create_taylor_couette_grid
    !
  else if (abs(itestcase) == Pipe_Channel_Flow) then
    !
    call create_pipe_channel_grid
    !
  else if (abs(itestcase) == CH4_Flame_Speed_Test) then
    !
    if (flame_speed_test_dof <= 0) then
      flame_speed_test_dof = 512
    end if
    !
    call create_flame_speed_test_grid
    !
  else if (abs(itestcase) == Sod_Shock_Tube) then
    !
    if (sod_test_dof == 0) sod_test_dof = 512
    !
    call create_sod_test_grid
    !
  else if (abs(itestcase) == Density_Pulse_Test) then
    !
    if (pulse_test_dof == 0) pulse_test_dof = 512
    !
    call create_pulse_test_grid
    !
  else if (abs(itestcase) == Shu_Osher_Test) then
    !
    if (shu_osher_test_dof == 0) shu_osher_test_dof = 512
    !
    call create_shu_osher_test_grid
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
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
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
subroutine reset_bc_input_for_tcf_grid(lx,rotation_angle,inner_wall_bc)
  !
  !.. Use Statements ..
  use ovar, only : bc_input,bc_input_t,periodic_bc_t
  use ovar, only : bc_names
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: lx,rotation_angle
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: inner_wall_bc
  !
  !.. Local Scalars ..
  integer :: i,n,nbc,bc_idx
  !
  type(periodic_bc_t) :: tmp_periodic_bc1
  type(periodic_bc_t) :: tmp_periodic_bc2
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: bcname
  !
  !.. Local Allocatable Arrays ..
  type(bc_input_t), allocatable :: tmp_bc_input(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reset_bc_input_for_tcf_grid"
  !
  logical(lk), parameter :: use_rotation_axis = fals
 !logical(lk), parameter :: use_rotation_axis = true
  !
continue
  !
  ! Create a default bc_input for the TCF flow case
  !
  ! F2003 AUTO-REALLOC
  tmp_bc_input = get_default_tcf_grid_bc_input(inner_wall_bc)
  !
  nbc = size(tmp_bc_input)
  !
  ! If bc_input is not allocated, use move_alloc to allocate
  ! bc_input to tmp_bc_input and return from the subroutine.
  !
  if (.not. allocated(bc_input)) then
    call move_alloc( from=tmp_bc_input , to=bc_input )
    return
  end if
  !
  ! Find the inflow BC in bc_input if it exists
  !
  bc_idx = 1
  n = 0
  tmp_periodic_bc1 = tmp_bc_input(bc_idx)%periodic
  find_inflow_loop: &
  do i = 1,size(bc_input)
    !
    bcname = trim(adjustl(uppercase(bc_input(i)%name))) ! F2003 AUTO-REALLOC
    if (bcname == "XMIN") n = i
    if (bcname == "INFLOW") n = i
    if (bc_input(i)%set_bc_value == bc_sem) n = i
    !
    if (n > 0) then
      tmp_bc_input(bc_idx) = bc_input(n)
      tmp_bc_input(bc_idx)%name = "INFLOW"
      exit find_inflow_loop
    end if
    !
  end do find_inflow_loop
  !
  ! Find the outflow bc in bc_input if it exists
  !
  bc_idx = 2
  n = 0
  tmp_periodic_bc2 = tmp_bc_input(bc_idx)%periodic
  find_outflow_loop: &
  do i = 1,size(bc_input)
    !
    bcname = trim(adjustl(uppercase(bc_input(i)%name))) ! F2003 AUTO-REALLOC
    if (bcname == "XMAX") n = i
    if (bcname == "OUTFLOW") n = i
    if (any(bc_input(i)%set_bc_value == bc_outflows(:))) n = i
    !
    if (n > 0) then
      tmp_bc_input(bc_idx) = bc_input(n)
      tmp_bc_input(bc_idx)%name = "OUTFLOW"
      exit find_outflow_loop
    end if
    !
  end do find_outflow_loop
  !
  ! If the inflow and outflow BCs are periodic, copy the default values for
  ! the periodic BCs back into tmp_bc_input
  !
  if (all(tmp_bc_input(1:2)%set_bc_value == bc_periodic)) then
    tmp_bc_input(1)%periodic = tmp_periodic_bc1
    tmp_bc_input(2)%periodic = tmp_periodic_bc2
  end if
  !
  ! Find the inner cylinder wall bc in bc_input if it exists
  !
  bc_idx = 3
  n = 0
  find_iwall_loop: &
  do i = 1,size(bc_input)
    !
    bcname = trim(adjustl(uppercase(bc_input(i)%name))) ! F2003 AUTO-REALLOC
    if (bcname == "INNER WALL") n = i
    if (bcname == "INNER_WALL") n = i
    if (bcname == "INNERWALL") n = i
    if (bcname == "BOTTOM WALL") n = i
    if (bcname == "BOTTOM_WALL") n = i
    if (bcname == "BOTTOMWALL") n = i
    if (bcname == "YMIN") n = i
    if (bcname == "RMIN") n = i
    !
    if (n > 0) then
      tmp_bc_input(bc_idx) = bc_input(n)
      tmp_bc_input(bc_idx)%name = "INNERWALL"
      exit find_iwall_loop
    end if
    !
  end do find_iwall_loop
  !
  ! Find the outer cylinder wall bc in bc_input if it exists
  !
  bc_idx = 4
  n = 0
  find_owall_loop: &
  do i = 1,size(bc_input)
    !
    bcname = trim(adjustl(uppercase(bc_input(i)%name))) ! F2003 AUTO-REALLOC
    if (bcname == "OUTER WALL") n = i
    if (bcname == "OUTER_WALL") n = i
    if (bcname == "OUTERWALL") n = i
    if (bcname == "TOP WALL") n = i
    if (bcname == "TOP_WALL") n = i
    if (bcname == "TOPWALL") n = i
    if (bcname == "YMAX") n = i
    if (bcname == "RMAX") n = i
    !
    if (n > 0) then
      tmp_bc_input(bc_idx) = bc_input(n)
      tmp_bc_input(bc_idx)%name = "OUTERWALL"
      exit find_owall_loop
    end if
    !
  end do find_owall_loop
  !
  ! Find the rotationally periodic BCs in bc_input
  ! if running only a sector of the TCF problem
  !
  if (nbc > 4) then
    !
    ! Find the theta-min rotationally periodic bc in bc_input if it exists
    !
    bc_idx = 5
    n = 0
    tmp_periodic_bc1 = tmp_bc_input(bc_idx)%periodic
    find_thetamin_loop: &
    do i = 1,size(bc_input)
      !
      bcname = trim(adjustl(uppercase(bc_input(i)%name))) ! F2003 AUTO-REALLOC
      if (bcname == "ZMIN") n = i
      if (bcname == "THETAMIN") n = i
      !
      if (n > 0) then
        tmp_bc_input(bc_idx) = bc_input(n)
        tmp_bc_input(bc_idx)%name = "THETAMIN"
        exit find_thetamin_loop
      end if
      !
    end do find_thetamin_loop
    !
    ! Find the theta-max rotationally periodic bc in bc_input if it exists
    !
    bc_idx = 6
    n = 0
    tmp_periodic_bc2 = tmp_bc_input(bc_idx)%periodic
    find_thetamax_loop: &
    do i = 1,size(bc_input)
      !
      bcname = trim(adjustl(uppercase(bc_input(i)%name))) ! F2003 AUTO-REALLOC
      if (bcname == "ZMAX") n = i
      if (bcname == "THETAMAX") n = i
      !
      if (n > 0) then
        tmp_bc_input(bc_idx) = bc_input(n)
        tmp_bc_input(bc_idx)%name = "THETAMAX"
        exit find_thetamax_loop
      end if
      !
    end do find_thetamax_loop
    !
    ! Copy the default values for the periodic BCs back into tmp_bc_input
    !
    tmp_bc_input(5)%periodic = tmp_periodic_bc1
    tmp_bc_input(6)%periodic = tmp_periodic_bc2
    !
  end if
  !
  ! Finally, use move_alloc to replace whatever might be in bc_input
  ! with the reconstructed complete bc_input contained in tmp_bc_input.
  !
  call move_alloc( from=tmp_bc_input , to=bc_input )
  !
  ! Reset bc_names to the names just assigned to bc_input for each BC
  !
  bc_names = bc_input(1:nbc)%name ! F2003 AUTO-REALLOC
  !
  contains
  !
  pure function get_default_tcf_grid_bc_input(iw_bc) result(return_value)
    !
    !.. Optional Arguments ..
    integer, optional, intent(in) :: iw_bc
    !
    !.. Function Result ..
    type(bc_input_t), allocatable :: return_value(:)
    !
    !.. Local Scalars ..
    integer :: n,ierr
    !
  continue
    !
    n = merge( 6 , 4 , rotation_angle > zero )
    !
    allocate ( return_value(1:n) , stat=ierr )
    !
    return_value(1)%name = "INFLOW"
    return_value(1)%set_bc_value = bc_periodic
    return_value(1)%periodic%matching_bc_name = "OUTFLOW"
    return_value(1)%periodic%matching_bc_idx  = 2
    return_value(1)%periodic%translation(1)   = lx
    !
    return_value(2)%name = "OUTFLOW"
    return_value(2)%set_bc_value = bc_periodic
    return_value(2)%periodic%matching_bc_name = "INFLOW"
    return_value(2)%periodic%matching_bc_idx  = 1
    return_value(2)%periodic%translation(1)   = -lx
    !
    return_value(3)%name = "INNERWALL"
    return_value(3)%set_bc_value = bc_isothermal_wall
    if (present(iw_bc)) then
      if (any(iw_bc == bc_heirarchy(:))) then
        ! Use iw_bc as the BC if iw_bc is a valid physical BC
        return_value(3)%set_bc_value = iw_bc
      end if
    end if
    if (return_value(3)%set_bc_value == bc_isothermal_wall) then
      return_value(3)%wall_temp = three*ten2 ! 300.0 K
    end if
    !
    return_value(4)%name = "OUTERWALL"
    return_value(4)%set_bc_value = bc_adiabatic_wall
    !
    ! Create the rotationally periodic boundary conditions
    ! if running only a sector of the TCF problem
    !
    if (n > 4) then
      !
      ! Create the rotationally periodic boundary conditions
      !
      return_value(5)%name = "THETAMIN"
      return_value(5)%set_bc_value = bc_periodic
      return_value(5)%periodic%matching_bc_name = "THETAMAX"
      return_value(5)%periodic%matching_bc_idx  = 6
      if (use_rotation_axis) then
        return_value(5)%periodic%RotationAxis(1) = one
        return_value(5)%periodic%RotationAxis(4) = rotation_angle
      else
        return_value(5)%periodic%RotationAngle(1) = rotation_angle
      end if
      !
      return_value(6)%name = "THETAMAX"
      return_value(6)%set_bc_value = bc_periodic
      return_value(6)%periodic%matching_bc_name = "THETAMIN"
      return_value(6)%periodic%matching_bc_idx  = 5
      if (use_rotation_axis) then
        return_value(6)%periodic%RotationAxis(1) = one
        return_value(6)%periodic%RotationAxis(4) = -rotation_angle
      else
        return_value(6)%periodic%RotationAngle(1) = -rotation_angle
      end if
      !
    end if
    !
  end function get_default_tcf_grid_bc_input
  !
end subroutine reset_bc_input_for_tcf_grid
!
!###############################################################################
!
subroutine create_taylor_couette_grid()
  !
  !.. Use Statements ..
  use order_mod,   only : n_order
  use ovar,        only : Lref,bc_input
  use geomvar_mod, only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geomvar_mod, only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geomvar_mod, only : grid,cell_geom,cell_order
  use ovar,        only : grid_format,bc_names,itestcase,internal_grid_order
  !
  use ovar, only : TCF_Inner_Radius,TCF_Inner_Wall_Spacing
  use ovar, only : TCF_Outer_Radius,TCF_Outer_Wall_Spacing
  use ovar, only : TCF_Aspect_Ratio,TCF_Axial_Cells
  use ovar, only : TCF_Radial_Cells,TCF_Azimuthal_Cells
  use ovar, only : TCF_Periodic_Sectors
  !
  use cgns_types_mod, only : HEXA_8,HEXA_27,HEXA_64,HEXA_125
  !
  use cgns_mod, only : get_face_elem
  !
  !.. Local Scalars ..
  integer  :: ierr,i,j,k,imax,jmax,azp1
  integer  :: nc,nf,n1,n2,bc1_idx,bc2_idx
  integer  :: ncx,ncr,nca,npx,npr,npa,npts
  integer  :: host_cell,host_side,elem_type
  integer  :: inner_wall_idx,outer_wall_idx
  integer  :: inflow_idx,outflow_idx
  integer  :: theta_min_idx,theta_max_idx
  integer  :: grid_order
  real(wp) :: lx,special_b_coeff,angle
  !
  logical(lk) :: rotationally_periodic
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
  ! Get the length of the cylinders
  !
  lx = TCF_Aspect_Ratio * abs(TCF_Outer_Radius-TCF_Inner_Radius)
  !
  ! Get the periodic rotation angle if running only a sector of the TCF problem
  ! NOTE: First make sure that TCF_Periodic_Sectors is one or greater
  !
  TCF_Periodic_Sectors = max(1,TCF_Periodic_Sectors)
  !
  angle = -huge(zero)
  rotationally_periodic = fals
  if (TCF_Periodic_Sectors > 1) then
    angle = three60/real(TCF_Periodic_Sectors,kind=wp)
    rotationally_periodic = true
  end if
  !
  ! Reset bc_input for this axially (and possibly rotationally) periodic grid
  ! so that the subroutine match_periodic_faces can correctly match all the
  ! periodic faces.
  !
  call reset_bc_input_for_tcf_grid(lx,angle)
  !
  ! Find the index in bc_input for the 4 (maybe 6) different boundary conditions
  !
  inner_wall_idx = -1
  outer_wall_idx = -1
  inflow_idx = -1
  outflow_idx = -1
  theta_min_idx = merge(-1,0,rotationally_periodic)
  theta_max_idx = merge(-1,0,rotationally_periodic)
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
    else if (trim(adjustl(uppercase(bc_names(i)))) == "THETAMIN") then
      theta_min_idx = i
    else if (trim(adjustl(uppercase(bc_names(i)))) == "THETAMAX") then
      theta_max_idx = i
    end if
  end do
  !
  ! If any of the above BC indices are still -1, there is a missing BC in
  ! bc_input needed for the TCF problem, so write out an error for that
  ! missing BC and stop the simulation.
  !
  if (any([inner_wall_idx,outer_wall_idx,inflow_idx, &
           outflow_idx,theta_min_idx,theta_max_idx] == -1)) then
    if (mypnum == glb_root) then
      write (iout,10)
      if (inner_wall_idx == -1) write (iout,11) "INNERWALL"
      if (outer_wall_idx == -1) write (iout,11) "OUTERWALL"
      if (inflow_idx == -1)     write (iout,11) "INFLOW"
      if (outflow_idx == -1)    write (iout,11) "OUTFLOW"
      if (theta_min_idx == -1)  write (iout,11) "THETAMIN"
      if (theta_max_idx == -1)  write (iout,11) "THETAMAX"
      write (iout,12)
    end if
    flush (iout)
    write (error_message,13)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
 !!
 !! I believe the following logic is no longer needed after adding
 !! the subroutine reset_bc_input_for_tcf_grid.
 !!
 !if (abs(itestcase) == Freestream_TCF) then
 !  inner_wall_bc = bc_input(inner_wall_idx)%set_bc_value
 !  outer_wall_bc = bc_input(outer_wall_idx)%set_bc_value
 !else
 !  inner_wall_bc = bc_isothermal_wall
 !  outer_wall_bc = bc_adiabatic_wall
 !end if
  !
  ! Reset the grid format input parameter to identify
  ! the grid as internally generated
  !
  grid_format = Internally_Generated_Grid
  !
  ! Get the grid order if specified for an internally
  ! generated grid, otherwise set it to the solution order.
  !
  grid_order = merge( internal_grid_order , n_order , internal_grid_order>0 )
  !
  ! For now we are limiting the grid order to
  ! the max order allowed by CGNS, which is 4
  !
  grid_order = min(grid_order,4)
  !
  ! If we are generating the grid on the fly, make the reference length 1
  !
  Lref = one
  !
  ! The Taylor-Couette flow problem is always 3D
  !
  nr = 3
  !
  ! Get the number of grid cells in each direction from the TCF inputs
  !
  ncr = TCF_Radial_Cells
  nca = max(1,TCF_Azimuthal_Cells/TCF_Periodic_Sectors)
  ncx = TCF_Axial_Cells
  !
  ! If the TCF grid is rotationally periodic, the grid will not be
  ! self-connecting and computing the number of nodes from the
  ! number of cells needs a +1, i.e.,
  !   azp1 = 1    if rotationally periodic in theta
  !   azp1 = 0    if self-connecting in theta (full 360 degree)
  !
  azp1 = merge(1,0,rotationally_periodic)
  !
  ! Get the number of grid nodes in each direction
  !
  npr = ncr + 1    ! radial nodes
  npa = nca + azp1 ! azimuthal nodes
  npx = ncx + 1    ! axial nodes
  !
  ncell = ncr*nca*ncx ! total number of grid cells
  ! The term azp1 used to compute the number of azimuthal nodes, npa, accounts
  ! for whether this grid is self-connecting (full 360 degrees) or only a
  ! sector using rotationally periodic BCs, so we can compute the total number
  ! of grid nodes in the normal way for structured grids.
  nnode = npr*npa*npx ! total number of grid nodes
  ! total number of boundary faces
  nfbnd = ncx*nca &      ! number of boundary faces for inner cylinder
        + ncx*nca &      ! number of boundary faces for outer cylinder
        + ncr*nca &      ! number of boundary faces for  inflow/periodic-1.A
        + ncr*nca &      ! number of boundary faces for outflow/periodic-1.B
        + ncx*ncr*azp1 & ! number of boundary faces for kmin/periodic-2.A
        + ncx*ncr*azp1   ! number of boundary faces for kmax/periodic-2.A
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
  call high_order_TCF_grid_points(ncx,ncr,nca,lx,angle,radial_pts,grid_order, &
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
  ! Loop through the 4 (possibly 6) grid boundaries
  ! NOTE: This corresponds to the grid boundaries in
  !       the order [InnerWall,OuterWall,Inflow,Outflow] which
  !       correspond to the [jmin,jmax,imin,imax] boundaries.
  ! NOTE: If running only a sector of the TCF problem, there will
  !       be kmin and kmax boundaries that are rotationally periodic.
  !       If running a full 360 degree version of the TCF problem,
  !       there are no kmin and kmax boundaries because
  !       the k-direction self-connects with itself.
  !
  ! InnerWall/j-min and OuterWall/j-max boundaries
  !
  n2 = 0
  !
  bc1_idx = inner_wall_idx
  bc2_idx = outer_wall_idx
  !
  n1 = n2           ! j-min (r-min) counter
  n2 = n2 + ncx*nca ! j-max (r-max) counter
  do k = 1,nca
    do i = 1,ncx
      !
      n1 = n1 + 1 ! j-min (r-min) face
      n2 = n2 + 1 ! j-max (r-max) face
      !
      bface( 1,n1) = bc_input(bc1_idx)%set_bc_value
      bface( 2,n1) = map3dto1d(i,  1,k,ncx,ncr,nca)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(3)
      bface( 5,n1) = bc1_idx ! group for boundary face
      bface( 6,n1) = bc1_idx ! index for bc_input for this BC
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_input(bc2_idx)%set_bc_value
      bface( 2,n2) = map3dto1d(i,ncr,k,ncx,ncr,nca)
      bface( 3,n2) = Geom_Hexa
      bface( 4,n2) = HostFace3D(4)
      bface( 5,n2) = bc2_idx ! group for boundary face
      bface( 6,n2) = bc2_idx ! index for bc_input for this BC
      bface(10,n2) = ncell + n2
      !
    end do
  end do
  !
  ! Min-X and Max-X boundary faces
  !
  bc1_idx = inflow_idx
  bc2_idx = outflow_idx
  !
  n1 = n2           ! min-x counter
  n2 = n2 + ncr*nca ! max-x counter
  do k = 1,nca
    do j = 1,ncr
      !
      n1 = n1 + 1 ! x-min face
      n2 = n2 + 1 ! x-max face
      !
      ! NOTE: The subroutine match_periodic_faces will set the values of
      !       bface([5,6,8,9],:), so we dont need to do that here.
      !
      bface( 1,n1) = bc_input(bc1_idx)%set_bc_value
      bface( 2,n1) = map3dto1d(  1,j,k,ncx,ncr,nca)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(1)
      bface( 5,n1) = bc1_idx ! group for boundary face
      bface( 6,n1) = bc1_idx ! index for bc_input needed by match_periodic_faces
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_input(bc2_idx)%set_bc_value
      bface( 2,n2) = map3dto1d(ncx,j,k,ncx,ncr,nca)
      bface( 3,n2) = Geom_Hexa
      bface( 4,n2) = HostFace3D(2)
      bface( 5,n2) = bc2_idx ! group for boundary face
      bface( 6,n2) = bc2_idx ! index for bc_input needed by match_periodic_faces
      bface(10,n2) = ncell + n2
      !
    end do
  end do
  !
  ! Min-Theta/K and Max-Theta/K rotationally periodic boundary faces
  !
  ! Create imax and jmax so that the following loop to create the rotationally
  ! periodic boundary faces in bface is skipped if running the full 360 degree
  ! version of the TCF problem (i.e., rotationally_periodic=F).
  !
  imax = merge(ncx,0,rotationally_periodic)
  jmax = merge(ncr,0,rotationally_periodic)
  !
  bc1_idx = theta_min_idx
  bc2_idx = theta_max_idx
  !
  n1 = n2           ! min-theta/k counter
  n2 = n2 + ncx*ncr ! max-theta/k counter
  do j = 1,jmax
    do i = 1,imax
      !
      n1 = n1 + 1 ! x-min face
      n2 = n2 + 1 ! x-max face
      !
      ! NOTE: The subroutine match_periodic_faces will set the values of
      !       bface([5,6,8,9],:), so we dont need to do that here.
      !
      bface( 1,n1) = bc_input(bc1_idx)%set_bc_value
      bface( 2,n1) = map3dto1d(i,j,  1,ncx,ncr,nca)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(5)
      bface( 5,n1) = bc1_idx ! group for boundary face
      bface( 6,n1) = bc1_idx ! index for bc_input needed by match_periodic_faces
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_input(bc2_idx)%set_bc_value
      bface( 2,n2) = map3dto1d(i,j,nca,ncx,ncr,nca)
      bface( 3,n2) = Geom_Hexa
      bface( 4,n2) = HostFace3D(6)
      bface( 5,n2) = bc2_idx ! group for boundary face
      bface( 6,n2) = bc2_idx ! index for bc_input needed by match_periodic_faces
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
  ! Only have the host roots create the grid derived type. The high-order
  ! data it will contain can take up a lot of memory. Only having the root
  ! process on each host/node store this data will save a lot of memory
  ! instead of having [# of cores per node] copies of the same data on each
  ! node. The partitioning code only expects the host root to have this
  ! information, so storing this information on non-root processes is just
  ! a waste of memory anyways.
  !
  if (i_am_host_root) then
    !
    allocate ( grid , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( grid%elem(1:ncell+nfbnd) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid%elem",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Create grid for the grid cells
    !
    if (grid_order == 1) then
      elem_type = int(HEXA_8)
      npts = 8
    else if (grid_order == 2) then
      elem_type = int(HEXA_27)
      npts = 27
    else if (grid_order == 3) then
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
      grid%elem(nc)%prop%order = grid_order
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
      !            cells %pts to put it back in CGNS order, then return the
      !            face elements %pts array in GFR ordering
      !
      if (grid%elem(host_cell)%prop%order <= 4) then
        grid%elem(nc) = get_face_elem(grid%elem(host_cell),host_side,host_cell)
      else
       !grid%elem(nc) = get_face_elem(grid%elem(host_cell),host_side,host_cell)
        write (error_message,1)
        call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
      end if
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
  end if
  !
 !if (mypnum == glb_root) then
 !  write (iout,2) minval( grid%elem(:)%prop%order ), &
 !                 maxval( grid%elem(:)%prop%order ), &
 !                 n_order
 !end if
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
  ! Format Statements
  !
  1 format ("ERROR: The internal grid generator for creating the ", &
                   "Taylor-Couette flow grid is currently limited to ", &
                   "creating P1-P4 order grids.")
  !
  2 format (/,4x,"After creating the Taylor-Couette flow grid:",/, &
             10x,"Min Grid Order = ",i0,/, &
             10x,"Max Grid Order = ",i0,/, &
             10x,"Solution Order = ",i0,/)
  !
  10 format (/,80("#"),//, &
               2x,"ERROR: There is a missing boundary condition that ", &
                  "is needed",/,9x,"to run the Taylor-Couette problem:")
  11 format (12x,"The ",a9," boundary condition is missing!")
  12 format (/,80("#"),/)
  13 format ("WHEN TRYING TO INTERNALLY GENERATE THE GRID FOR THE ", &
             "TAYLOR-COUETTE FLOW PROBLEM, THERE IS A MISSING BOUNDARY ", &
             "CONDITION THAT IS NEEDED TO BE ABLE TO RUN THIS PROBLEM!")
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
subroutine high_order_TCF_grid_points(ncx,ncr,nca,length,angle,radial_pts, &
                                      order,elem_nodes,elem_pts,xyz_nodes, &
                                      xyz_pts)
  !
  !.. Use Statements ..
  use ovar, only : TCF_Theta_Shift
  !
  !.. Formal Arguments ..
  integer,                intent(in) :: ncx
  integer,                intent(in) :: ncr
  integer,                intent(in) :: nca
  real(wp),               intent(in) :: length
  real(wp),               intent(in) :: angle
  real(wp),               intent(in) :: radial_pts(:)
  integer,                intent(in) :: order
  integer,  allocatable, intent(out) :: elem_nodes(:,:)
  integer,  allocatable, intent(out) :: elem_pts(:,:)
  real(wp), allocatable, intent(out) :: xyz_nodes(:,:)
  real(wp), allocatable, intent(out) :: xyz_pts(:,:)
  !
  !.. Local Scalars ..
  integer  :: nnx,nnr,nna,npx,npr,npa,p,azp1
  integer  :: i,j,k,ic,jc,kc,ip,jp,kp
  integer  :: n,nc,lcell,lnode,cpts,npts
  real(wp) :: y,z,theta,divp,dr,dlpa,dlpx,dlca
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: ijkmap(:,:,:)
  real(wp), allocatable :: rpts(:)
  real(wp), allocatable :: xpts(:)
  real(wp), allocatable :: apts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "high_order_TCF_grid_points"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! If angle > zero, the TCF grid is rotationally periodic which means the
  ! grid will not be self connecting and computing the number of nodes
  ! from the number of cells needs a +1.
  !
  azp1 = merge(1,0,angle>zero)
  !
  !
  ! Number of grid nodes in each direction
  !
  nnx = ncx + 1    ! axial     grid nodes
  nnr = ncr + 1    ! radial    grid nodes
  nna = nca + azp1 ! azimuthal grid nodes
 !!
 !! Get the order of the high order grid
 !! NOTE: For now we are limiting this to P4 to match the maximum allowed
 !!       by CGNS. One reason to do this is because we can then use the
 !!       functions in cgns_mod to extract the high-order connectivity for
 !!       each of the boundary face elements, which are needed for the
 !!       grid derived type.
 !!
 !p = min(order,4)
  p = order
  !
  ! Number of high order points in each direction
  !
  npx = ncx*p + 1    ! axial     high order points
  npr = ncr*p + 1    ! radial    high order points
  npa = nca*p + azp1 ! azimuthal high order points
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
  if (angle > zero) then
    dlpa = (angle*degrees2radians) / real(npa-1,kind=wp)
  else
    dlpa = two*pi / real(npa,kind=wp)
  end if
  ! Compute the uniform azimuthal cell spacing
  dlca = dlpa*real(p,kind=wp)
  !
  xpts = dlpx * real( intseq(0,npx-1) , kind=wp ) ! F2003 AUTO-REALLOC
  apts = dlpa * real( intseq(0,npa-1) , kind=wp ) ! F2003 AUTO-REALLOC
  !
  ! Shift the azimuthal points (only do this for full 360 grid version)
  !
  if (abs(TCF_Theta_Shift) > eps6) then
    ! NOTE: TCF_Theta_Shift gives the shift relative to the uniform cell size
    apts(:) = apts(:) + TCF_Theta_Shift*dlca
    ! Move any theta coordinates beyond 2*pi back to the interval [0,2*pi]
    where (apts > two*pi) apts = apts - two*pi
  end if
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
  ijkmap = init_array(npx,npr,npa+1-azp1,-1)
  !
  ! First create the indices for the linear grid nodes, and also
  ! fill in xyz_nodes for these nodes
  !
  n = 0
  do k = 1,nna
    !
    kp = (k-1)*p + 1
    !
    theta = apts(kp)
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
        xyz_nodes(1,n) = xpts(ip)
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
  ! Only the root process on each host/node will store the high-order
  ! grid information prior to partitioning, so we only want the host
  ! roots to create xyz_pts.
  !
  if (i_am_host_root) then
    !
    xyz_pts = init_array(3,npts,zero) ! F2003 AUTO-REALLOC
    !
    ! Copy xyz_nodes over to xyz_pts
    ! NOTE: n should equal lnode coming out of the above loops
    !
    xyz_pts(1:3,1:lnode) = xyz_nodes(1:3,1:lnode)
    !
  end if
  !
  ! Now create the indices for the high-order grid points, and also
  ! fill in xyz_pts for these points
  !
  n = lnode ! n should already be equal to lnode, but just in case
  do kp = 1,npa
    theta = apts(kp)
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
        ! NOTE: Only the host roots will store the high-order data
        !
        if (i_am_host_root) then
          xyz_pts(1,n) = xpts(ip)
          xyz_pts(2,n) = rpts(jp)*cos(theta)
          xyz_pts(3,n) = rpts(jp)*sin(theta)
        end if
        !
      end do
    end do
  end do
  !
  ! Finally, if not using a rotationally periodic grid, fill the
  ! self-connected kmax plane in ijkmap by copying over the indices
  ! from the kmin plane
  !
  if (angle <= zero) ijkmap(:,:,npa+1) = ijkmap(:,:,1)
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
  ! Only the root process on each host/node will store the high-order
  ! grid information prior to partitioning, so we only want the host
  ! roots to create elem_pts.
  !
  if (i_am_host_root) then
    elem_pts = init_array(cpts,lcell,0) ! F2003 AUTO-REALLOC
  end if
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
        ! NOTE: Only the host roots will store the high-order data
        !
        if (i_am_host_root) then
          n = 0
          do kp = k,k+p
            do jp = j,j+p
              do ip = i,i+p
                n = n + 1
                elem_pts(n,nc) = ijkmap(ip,jp,kp)
              end do
            end do
          end do
        end if
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
  use geomvar_mod, only : nr,grid,bface
  !
  !.. Formal Arguments ..
  integer, intent(in) :: lcell
  integer, intent(in) :: lfbnd
  !
  !.. Local Scalars ..
  integer :: nnod,npts,fnod,fpts
  integer :: nc,nf,i,l
  integer :: i1,i2,i3,i4,i5,i6,i7,i8
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
  1 format (i6,1x,"|",1x,100i8)
  2 format (i6,1x,"|",1x,i3,i2,1x,"|",1x,100i8)
  3 format (i7,2(5x,3f12.6))
  4 format (89("="),/,"Cell ",i0, &
            t21,"x",t33,"y",t45,"z",t62,"x",t74,"r",t86,"t",/, &
            89("-"))
  5 format (89("="),/,"Boundary Face ",i0, &
            t21,"x",t33,"y",t45,"z",t62,"x",t74,"r",t86,"t",/, &
            89("-"))
  6 format (1000(a1))
  !
  if (mypnum == glb_root) then
    !
    open (newunit=i1,file="p0-cell_conn.dat",status="replace",action="write")
    open (newunit=i2,file="ho-cell_conn.dat",status="replace",action="write")
    open (newunit=i3,file="p0-face_conn.dat",status="replace",action="write")
    open (newunit=i4,file="ho-face_conn.dat",status="replace",action="write")
    open (newunit=i5,file="p0-cell_xyz.dat",status="replace",action="write")
    open (newunit=i6,file="ho-cell_xyz.dat",status="replace",action="write")
    open (newunit=i7,file="p0-face_xyz.dat",status="replace",action="write")
    open (newunit=i8,file="ho-face_xyz.dat",status="replace",action="write")
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
    write (i1,1) 0,    (i, i=1,nnod)
    write (i3,2) 0,0,0,(i, i=1,fnod)
    write (i2,1) 0,    (i, i=1,npts)
    write (i4,2) 0,0,0,(i, i=1,fpts)
    !
    write (i2,6) ("-",i=1,npts*8+9)
    write (i4,6) ("-",i=1,fpts*8+17)
    !
    do nc = 1,lcell
      !
      nodes = grid%elem(nc)%nodes(:) ! F2003 AUTO-REALLOC
      pts   = grid%elem(nc)%pts(:)   ! F2003 AUTO-REALLOC
      !
      write (i1,1) nc,(nodes(i), i=1,size(nodes))
      write (i2,1) nc,(pts(i), i=1,size(pts))
      !
      write (i5,4) nc
      do i = 1,size(nodes)
        l = nodes(i)
        x = grid%xyz(1,l)
        y = grid%xyz(2,l)
        z = grid%xyz(3,l)
        r = sqrt(y*y+z*z)
        t = atan2(z,y)
        t = t + merge(two*pi,zero,t<zero)
        write (i5,3) l,x,y,z,x,r,t*radians2degrees
      end do
      !
      write (i6,4) nc
      do i = 1,size(pts)
        l = pts(i)
        x = grid%xyz(1,l)
        y = grid%xyz(2,l)
        z = grid%xyz(3,l)
        r = sqrt(y*y+z*z)
        t = atan2(z,y)
        t = t + merge(two*pi,zero,t<zero)
        write (i6,3) l,x,y,z,x,r,t*radians2degrees
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
      write (i3,2) nf,bface(1,nf),bface(4,nf),(nodes(i), i=1,size(nodes))
      write (i4,2) nf,bface(1,nf),bface(4,nf),(pts(i), i=1,size(pts))
      !
      write (i7,5) nf
      do i = 1,size(nodes)
        l = nodes(i)
        x = grid%xyz(1,l)
        y = grid%xyz(2,l)
        z = grid%xyz(3,l)
        r = sqrt(y*y+z*z)
        t = atan2(z,y)
        t = t + merge(two*pi,zero,t<zero)
        write (i7,3) l,x,y,z,x,r,t*radians2degrees
      end do
      !
      write (i8,5) nf
      do i = 1,size(pts)
        l = pts(i)
        x = grid%xyz(1,l)
        y = grid%xyz(2,l)
        z = grid%xyz(3,l)
        r = sqrt(y*y+z*z)
        t = atan2(z,y)
        t = t + merge(two*pi,zero,t<zero)
        write (i8,3) l,x,y,z,x,r,t*radians2degrees
      end do
      !
    end do
    !
    close (i1)
    close (i2)
    close (i3)
    close (i4)
    close (i5)
    close (i6)
    close (i7)
    close (i8)
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
  use order_mod,   only : n_order
  use ovar,        only : Lref,itestcase
  use ovar,        only : VortexGrid_DOF,VortexGrid_Lref
  use ovar,        only : VortexGrid_random_perturb
  use ovar,        only : VortexGrid_perturb_percent
  use geomvar_mod, only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geomvar_mod, only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geomvar_mod, only : cell_geom,cell_order
  use ovar,        only : grid_format
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
  use order_mod,   only : n_order
  use ovar,        only : VortexGrid_DOF,Lref
  use ovar,        only : VortexGrid_random_perturb
  use ovar,        only : VortexGrid_perturb_percent
  use geomvar_mod, only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geomvar_mod, only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geomvar_mod, only : cell_geom,cell_order
  use ovar,        only : grid_format
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
  use order_mod,   only : n_order,geom_solpts
  use ovar,        only : VortexGrid_DOF,Lref
  use ovar,        only : VortexGrid_random_perturb
  use ovar,        only : VortexGrid_perturb_percent
  use geomvar_mod, only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geomvar_mod, only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geomvar_mod, only : cell_geom,cell_order
  use ovar,        only : grid_format
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
subroutine reset_bc_input_for_internal_grid(Lref_x,Lref_y,Lref_z)
  !
  ! Reset bc_input for a fully periodic grid so that the subroutine
  ! match_periodic_faces can correctly match all the periodic faces.
  !
  !.. Use Statements ..
  use ovar,   only : bc_input
  use geomvar_mod, only : nr
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: Lref_x  !< x-translation distance (def=1.)
  real(wp), optional, intent(in) :: Lref_y  !< y-translation distance (def=1.)
  real(wp), optional, intent(in) :: Lref_z  !< z-translation distance (def=1.)
  !
  !.. Local Scalars ..
  integer  :: ierr
  real(wp) :: lx,ly,lz
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reset_bc_input_for_internal_grid"
  !
continue
  !
  lx = one
  ly = one
  lz = one
  !
  if (present(Lref_x)) lx = Lref_x
  if (present(Lref_y)) ly = Lref_y
  if (present(Lref_z)) lz = Lref_z
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
  bc_input(1)%periodic%translation(1)   = lx
  !
  bc_input(2)%name = "xmax"
  bc_input(2)%periodic%matching_bc_name = "xmin"
  bc_input(2)%periodic%matching_bc_idx  = 1
  bc_input(2)%periodic%translation(1)   = -lx
  !
  if (nr > 1) then
    bc_input(3)%name = "ymin"
    bc_input(3)%periodic%matching_bc_name = "ymax"
    bc_input(3)%periodic%matching_bc_idx  = 4
    bc_input(3)%periodic%translation(2)   = ly
    !
    bc_input(4)%name = "ymax"
    bc_input(4)%periodic%matching_bc_name = "ymin"
    bc_input(4)%periodic%matching_bc_idx  = 3
    bc_input(4)%periodic%translation(2)   = -ly
  end if
  !
  if (nr > 2) then
    bc_input(5)%name = "zmin"
    bc_input(5)%periodic%matching_bc_name = "zmax"
    bc_input(5)%periodic%matching_bc_idx  = 6
    bc_input(5)%periodic%translation(3)   = lz
    !
    bc_input(6)%name = "zmax"
    bc_input(6)%periodic%matching_bc_name = "zmin"
    bc_input(6)%periodic%matching_bc_idx  = 5
    bc_input(6)%periodic%translation(3)   = -lz
  end if
  !
end subroutine reset_bc_input_for_internal_grid
!
!###############################################################################
!
subroutine reset_bc_input_for_channel_grid(lx,lz)
  !
  !.. Use Statements ..
  use ovar, only : bc_input,bc_input_t
  use ovar, only : bc_names
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: lx,lz
  !
  !.. Local Scalars ..
  integer :: i,n,bc_idx
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: bcname
  !
  !.. Local Allocatable Arrays ..
  type(bc_input_t), allocatable :: tmp_bc_input(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reset_bc_input_for_channel_grid"
  !
continue
  !
  ! Create a default bc_input for the channel flow case
  !
  ! F2003 AUTO-REALLOC
  tmp_bc_input = get_default_channel_grid_bc_input()
  !
  ! If bc_input is not allocated, use move_alloc to allocate
  ! bc_input to tmp_bc_input and return from the subroutine.
  !
  if (.not. allocated(bc_input)) then
    call move_alloc( from=tmp_bc_input , to=bc_input )
    return
  end if
  !
  ! Find the inflow BC in bc_input if it exists
  !
  bc_idx = 1
  n = 0
  find_inflow_loop: &
  do i = 1,size(bc_input)
    !
    bcname = trim(adjustl(uppercase(bc_input(i)%name))) ! F2003 AUTO-REALLOC
    if (bcname == "XMIN") n = i
    if (bcname == "INFLOW") n = i
    if (bc_input(i)%set_bc_value == bc_sem) n = i
    !
    if (n > 0) then
      tmp_bc_input(bc_idx) = bc_input(n)
      tmp_bc_input(bc_idx)%name = "INFLOW"
      exit find_inflow_loop
    end if
    !
  end do find_inflow_loop
  !
  ! Find the outflow bc in bc_input if it exists
  !
  bc_idx = 2
  n = 0
  find_outflow_loop: &
  do i = 1,size(bc_input)
    !
    bcname = trim(adjustl(uppercase(bc_input(i)%name))) ! F2003 AUTO-REALLOC
    if (bcname == "XMAX") n = i
    if (bcname == "OUTFLOW") n = i
    if (any(bc_input(i)%set_bc_value == bc_outflows(:))) n = i
    !
    if (n > 0) then
      tmp_bc_input(bc_idx) = bc_input(n)
      tmp_bc_input(bc_idx)%name = "OUTFLOW"
      exit find_outflow_loop
    end if
    !
  end do find_outflow_loop
  !
  ! Find the bottom wall bc in bc_input if it exists
  !
  bc_idx = 3
  n = 0
  find_bwall_loop: &
  do i = 1,size(bc_input)
    !
    bcname = trim(adjustl(uppercase(bc_input(i)%name))) ! F2003 AUTO-REALLOC
    if (bcname == "BOTTOM WALL") n = i
    if (bcname == "BOTTOM_WALL") n = i
    if (bcname == "BOTTOMWALL") n = i
    if (bcname == "YMIN") n = i
    !
    if (n > 0) then
      tmp_bc_input(bc_idx) = bc_input(n)
      tmp_bc_input(bc_idx)%name = "BOTTOM WALL"
      exit find_bwall_loop
    end if
    !
  end do find_bwall_loop
  !
  ! Find the top wall bc in bc_input if it exists
  !
  bc_idx = 4
  n = 0
  find_twall_loop: &
  do i = 1,size(bc_input)
    !
    bcname = trim(adjustl(uppercase(bc_input(i)%name))) ! F2003 AUTO-REALLOC
    if (bcname == "TOP WALL") n = i
    if (bcname == "TOP_WALL") n = i
    if (bcname == "TOPWALL") n = i
    if (bcname == "YMAX") n = i
    !
    if (n > 0) then
      tmp_bc_input(bc_idx) = bc_input(n)
      tmp_bc_input(bc_idx)%name = "TOP WALL"
      exit find_twall_loop
    end if
    !
  end do find_twall_loop
  !
  ! Find the z-min periodic bc in bc_input if it exists
  !
  bc_idx = 5
  n = 0
  find_zmin_loop: &
  do i = 1,size(bc_input)
    !
    bcname = trim(adjustl(uppercase(bc_input(i)%name))) ! F2003 AUTO-REALLOC
    if (bcname == "ZMIN") n = i
    !
    if (n > 0) then
      tmp_bc_input(bc_idx) = bc_input(n)
      tmp_bc_input(bc_idx)%name = "ZMIN"
      exit find_zmin_loop
    end if
    !
  end do find_zmin_loop
  !
  ! Find the z-max periodic bc in bc_input if it exists
  !
  bc_idx = 6
  n = 0
  find_zmax_loop: &
  do i = 1,size(bc_input)
    !
    bcname = trim(adjustl(uppercase(bc_input(i)%name))) ! F2003 AUTO-REALLOC
    if (bcname == "ZMAX") n = i
    !
    if (n > 0) then
      tmp_bc_input(bc_idx) = bc_input(n)
      tmp_bc_input(bc_idx)%name = "ZMAX"
      exit find_zmax_loop
    end if
    !
  end do find_zmax_loop
  !
  ! Finally, use move_alloc to replace whatever might be in bc_input
  ! with the reconstructed complete bc_input contained in tmp_bc_input.
  !
  call move_alloc( from=tmp_bc_input , to=bc_input )
  !
  ! Reset bc_names to the names just assigned to bc_input for each BC
  !
  bc_names = bc_input(1:6)%name ! F2003 AUTO-REALLOC
  !
  contains
  !
  pure function get_default_channel_grid_bc_input() result(return_value)
    !
    !.. Function Result ..
    type(bc_input_t) :: return_value(1:6)
    !
  continue
    !
    return_value(1)%name = "INFLOW"
    return_value(1)%set_bc_value = bc_periodic
    return_value(1)%periodic%matching_bc_name = "OUTFLOW"
    return_value(1)%periodic%matching_bc_idx  = 2
    return_value(1)%periodic%translation(1)   = lx
    !
    return_value(2)%name = "OUTFLOW"
    return_value(2)%set_bc_value = bc_periodic
    return_value(2)%periodic%matching_bc_name = "INFLOW"
    return_value(2)%periodic%matching_bc_idx  = 1
    return_value(2)%periodic%translation(1)   = -lx
    !
    return_value(3)%name = "BOTTOM WALL"
    return_value(3)%set_bc_value = bc_adiabatic_wall
    !
    return_value(4)%name = "TOP WALL"
    return_value(4)%set_bc_value = bc_adiabatic_wall
    !
    return_value(5)%name = "ZMIN"
    return_value(5)%set_bc_value = bc_periodic
    return_value(5)%periodic%matching_bc_name = "ZMAX"
    return_value(5)%periodic%matching_bc_idx  = 6
    return_value(5)%periodic%translation(3)   = lz
    !
    return_value(6)%name = "ZMAX"
    return_value(6)%set_bc_value = bc_periodic
    return_value(6)%periodic%matching_bc_name = "ZMIN"
    return_value(6)%periodic%matching_bc_idx  = 5
    return_value(6)%periodic%translation(3)   = -lz
    !
  end function get_default_channel_grid_bc_input
  !
end subroutine reset_bc_input_for_channel_grid
!
!###############################################################################
!
subroutine create_channel_flow_grid()
  !
  !.. Use Statements ..
  use order_mod,      only : n_order
  use ovar,           only : Lref
  use ovar,           only : loc_solution_pts
  use ovar,           only : channel_grid_ave_xplus
  use ovar,           only : channel_grid_ave_yplus
  use ovar,           only : channel_grid_ave_zplus
  use ovar,           only : channel_grid_yplus_wall
  use ovar,           only : channel_yplus_of_one
  use geomvar_mod,    only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geomvar_mod,    only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geomvar_mod,    only : cell_geom,cell_order
  use ovar,           only : grid_format,bc_input,bc_names
  !
  use poly_mod,       only : legendre_gauss_quadrature
  use poly_mod,       only : legendre_gauss_lobatto_quadrature
  !
  !.. Local Scalars ..
  integer  :: i,j,k,n,n1,n2
  integer  :: ncx,ncy,ncz,bc1_idx,bc2_idx
  integer  :: nhy
  real(wp) :: xplus_ave,yplus_ave,zplus_ave,dx,dz,yplus1
  real(wp) :: ds1,ds2,dl1,dl2,a,b,sf,rnum,rden,sferr,val,ue,se
  real(wp) :: lx_plus,ly_plus,lz_plus,yplus_wall,pp1
  !
  integer :: inflow_idx,outflow_idx
  integer :: bottom_wall_idx,top_wall_idx
  integer :: zmin_idx,zmax_idx
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: xi(:)
  real(wp), allocatable :: x(:)
  real(wp), allocatable :: y(:)
  real(wp), allocatable :: z(:)
  real(wp), allocatable :: yhalf(:)
  real(wp), allocatable :: y1(:)
  real(wp), allocatable :: y2(:)
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
  ! Reset bc_input for this axially and spanwise periodic grid so that the
  ! subroutine match_periodic_faces can correctly match all the periodic faces.
  !
  call reset_bc_input_for_channel_grid(channel_grid_Lx,channel_grid_Lz)
  !
  ! Find the index in bc_input for the 6 different boundary conditions
  !
  inflow_idx = -1
  outflow_idx = -1
  zmin_idx = -1
  zmax_idx = -1
  bottom_wall_idx = -1
  top_wall_idx = -1
  !
  do i = 1,size(bc_names)
    if (trim(adjustl(uppercase(bc_names(i)))) == "INFLOW") then
      inflow_idx = i
    else if (trim(adjustl(uppercase(bc_names(i)))) == "OUTFLOW") then
      outflow_idx = i
    else if (trim(adjustl(uppercase(bc_names(i)))) == "ZMIN") then
      zmin_idx = i
    else if (trim(adjustl(uppercase(bc_names(i)))) == "ZMAX") then
      zmax_idx = i
    else if (trim(adjustl(uppercase(bc_names(i)))) == "BOTTOM WALL") then
      bottom_wall_idx = i
    else if (trim(adjustl(uppercase(bc_names(i)))) == "TOP WALL") then
      top_wall_idx = i
    end if
  end do
  !
  ! If any of the above BC indices are still -1, there is a missing BC in
  ! bc_input needed for the TCF problem, so write out an error for that
  ! missing BC and stop the simulation.
  !
  if (any([bottom_wall_idx,top_wall_idx,inflow_idx, &
           outflow_idx,zmin_idx,zmax_idx] == -1)) then
    if (mypnum == glb_root) then
      write (iout,10)
      if (bottom_wall_idx == -1) write (iout,11) "BOTTOM WALL"
      if (top_wall_idx == -1)    write (iout,11) "TOP WALL"
      if (inflow_idx == -1)      write (iout,11) "INFLOW"
      if (outflow_idx == -1)     write (iout,11) "OUTFLOW"
      if (zmin_idx == -1)        write (iout,11) "ZMIN"
      if (zmax_idx == -1)        write (iout,11) "ZMAX"
      write (iout,12)
    end if
    flush (iout)
    write (error_message,13)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Precompute P+1 as a real
  !
  pp1 = real(n_order+1,kind=wp)
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
  if (channel_grid_ave_xplus  == zero) channel_grid_ave_xplus = -three*ten
  if (channel_grid_ave_yplus  == zero) channel_grid_ave_yplus = -three*five
  if (channel_grid_ave_zplus  == zero) channel_grid_ave_zplus = -three*five
  if (channel_grid_yplus_wall == zero) channel_grid_yplus_wall = two
  !
  ! Save the channel grid parameters to local variables
  !
  yplus1 = channel_yplus_of_one
  xplus_ave = channel_grid_ave_xplus
  yplus_ave = channel_grid_ave_yplus
  zplus_ave = channel_grid_ave_zplus
  yplus_wall = channel_grid_yplus_wall
  !
  ! Get the length of the domain in each direction in y+ dimensions
  !
  lx_plus = channel_grid_Lx / yplus1
  ly_plus = channel_grid_Ly / yplus1
  lz_plus = channel_grid_Lz / yplus1
  !
  ! Compute the cell average x+,y+,z+ for the grid
  ! NOTE: We are assuming the y+ cell average is the y+ spacing in the middle
  !       of the channel and not the actual average y+ spacing including the
  !       near wall cells. This will tend towards finer y+ spacing in the
  !       middle of the channel, but it makes the input more directly
  !       comparable to the x+ and z+ inputs. Otherwise, if the average y+
  !       spacing includes the near wall cells in the average, the spacing
  !       in the middle of the channel will be larger than the equivalent
  !       x+ or z+ spacing given by the respective inputs.
  !
  ! If the inputs channel_grid_ave_[xyz]plus are:
  !   > 0 : the input specifies the average cell spacing
  !   < 0 : the abs(input) specifies the average DoF spacing
  !
  ! For a given coordinate direction (say [i]), the cell average [i]+ is
  !   if (channel_grid_ave_[i]plus > zero) then
  !     [i]plus_ave = channel_grid_ave_[i]plus
  !   else
  !     [i]plus_ave = abs(channel_grid_ave_[i]plus)*real(n_order+1,kind=wp)
  !   end if
  ! NOTE: We use (P+1) here because the delta plus value is the
  !       average distance between the DoF within a cell (or the
  !       subcells), and there are P+1 subcells in each cell (the
  !       subcells between the end points and the cell boundaries
  !       are kind of 1/2 a subcell each).
  !
  xplus_ave = abs(xplus_ave)*merge(one , pp1 , xplus_ave > zero)
  yplus_ave = abs(yplus_ave)*merge(one , pp1 , yplus_ave > zero)
  zplus_ave = abs(zplus_ave)*merge(one , pp1 , zplus_ave > zero)
  !
  ! Get the number of cells in the x and z directions
  !
  ncx = nint( lx_plus / xplus_ave )
  ncz = nint( lz_plus / zplus_ave )
  !
  ! Compute the actual (non-plus) x and z cell spacings
  ! NOTE: These will most likely not be equal to xplus_ave and zplus_ave
  !       since lx_plus and lz_plus will probably not be evenly divisible
  !       by these respective numbers.
  !
  dx = channel_grid_Lx / real(ncx,kind=wp)
  dz = channel_grid_Lz / real(ncz,kind=wp)
  !
  ! Save the x+ and z+ spacings to the module variables
  !
  xplus_cell_ave = dx / yplus1
  zplus_cell_ave = dz / yplus1
  !
  ! Initialize the x and z grid spacing arrays
  !
  x = init_array(ncx+1,zero) ! F2003 AUTO-REALLOC
  z = init_array(ncz+1,zero) ! F2003 AUTO-REALLOC
  !
  ! #######
  ! #######
  ! ## 2 ##
  ! #######
  ! #######
  !
  ! Get the grid cell spacings for x and z
  !
  ! ######
  ! # 2A #
  ! ######
  !
  ! Get the grid point locations for x
  !
  do n = 2,ncx+1
    x(n) = x(n-1) + dx
  end do
  !
  ! ######
  ! # 2B #
  ! ######
  !
  ! Get the grid point locations for z
  !
  do n = 2,ncz+1
    z(n) = z(n-1) + dz
  end do
  !
  ! ######
  ! # 2C #
  ! ######
  !
  ! Get the grid point locations for y
  !
  y = channel_y_spacing(yplus_wall,yplus_ave,ly_plus) ! F2003 AUTO-REALLOC
  ncy = size(y)-1
  !
#ifdef DEBUG_ON
!#define DEBUG_CHANNEL_Y_SPACING
#ifdef DEBUG_CHANNEL_Y_SPACING
  y1 = old_channel_y_spacing(yplus_wall,yplus_ave,ly_plus) ! F2003 AUTO-REALLOC
  !
  ! 2C.0 Get the location of the 1D solution points in reference space
  !
  xi = init_array(n_order+1,zero) ! F2003 AUTO-REALLOC
  !
  if (loc_solution_pts == Legendre_Gauss_Lobatto) then
    call legendre_gauss_lobatto_quadrature(xi)
  else
    call legendre_gauss_quadrature(xi)
  end if
  !
  ! Scale xi from [-1,1] to [0,1]
  !
  xi(:) = half*(xi(:) + one)
  !
  ! 2C.1 : Get the beginning and ending spacings needed
  !        to generate the tanh spacing
  !
  ! Need to get dl1/ds1 and dl2/ds2 to create tanh spacing for y spacing
  !
  if (yplus_wall < zero) then
    !
   !! OLD VERSION: The 1st DoF off the wall is located at abs(yplus_wall)
   !!
   !! OLD VERSION:
   !! Divide the abs(yplus_wall) by the location of the 1st solution point
   !! on the interval [0,1] to get the height of the first cell off the wall.
   !!
   !dl1 = abs(yplus_wall) / xi(1)
    !
    ! NEW VERSION: The nint(abs(yplus_wall))-th DoF
    !              off the wall is located at y+=1.
    ! NEW VERSION:
    ! The location of y+=1 is at the nint(abs(yplus_wall))-th DoF off the wall,
    ! so 1/xi(nint(abs(yplus_wall))) gives the y+ spacing of the first cell off
    ! the wall.
    ! NOTE: This assumes xi is on the interval [0,1].
    ! NOTE: We need to use min/max to make sure nint(abs(yplus_wall))
    !       is a valid index for the array xi.
    !
    n = max( 1 , min( size(xi) , nint(abs(yplus_wall)) ) )
    dl1 = one / xi(n)
    !
  else
    !
    ! The height of the first cell off the wall is yplus_wall
    !
    dl1 = yplus_wall
    !
  end if
  !
  ! Get the centerline y+ spacing
  !
  dl2 = yplus_ave
  !
  ! Map the beginning and ending grid cell spacings to the reference
  ! interval [0,1]. We are only worried about half the channel so we
  ! need to half of the channel height in plus coordinates.
  !
  ds1 = dl1 / (half*ly_plus)
  ds2 = dl2 / (half*ly_plus)
  !
  ! Compute the approximate number of cells (for half of the channel height)
  ! needed to get a constant stretching with end spacings ds1 and ds2.
  !
  nhy = nint( log(ds2/ds1) / log((ds1-one)/(ds2-one)) )
  !
  ! The total number of y cells is just twice nhy,
  ! since nhy is for half of the channel height
  !
  ncy = nhy*2
  !
  ! Save all the various y+ spacings to the module variables
  !
  yplus_cell_wall = dl1
  yplus_dof_wall = xi(:)*dl1 ! F2003 AUTO-REALLOC
  yplus_cell_centerline = dl2
  yplus_cell_ave = ly_plus / real(ncy,kind=wp)
  !
  ! Initialize the y grid spacing for half the channel
  !
  yhalf = init_array(nhy+1,zero) ! F2003 AUTO-REALLOC
  !
  ! 2C.2 : Get the a and b parameters to begin solving for the tanh y-spacing
  !
  a = sqrt(ds2/ds1)
  b = one / (real(nhy,kind=wp)*sqrt(ds1*ds2))
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
    ! Get the distribution yhalf(n)
    !
    do n = 0,nhy
      val = real(n,kind=wp)/real(nhy,kind=wp)
      ue = half*(one + tanh(sf*(val-half))/tanh(sf*half))
      se = ue / (a + (one-a)*ue)
      yhalf(n+1) = se * half*channel_grid_Ly
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
    ! Get the distribution yhalf(n)
    !
    do n = 0,nhy
      val = real(n,kind=wp)/real(nhy,kind=wp)
      ue = half*(one + tan(sf*(val-half))/tan(sf*half))
      se = ue / (a + (one-a)*ue)
      yhalf(n+1) = se * half*channel_grid_Ly
    end do
    !
  else
    !
    ! Both of the above Newton iterations can have numerical difficulties
    ! if the stretching factor is close to 1, therefore use this simpler
    ! stretching function that doesnt require a Newton iteration.
    !
    do n = 0,nhy
      val = real(n,kind=wp)/real(nhy,kind=wp)
      ue = val*(one + two*(b-one)*(val-half)*(one-val))
      se = ue / (a + (one-a)*ue)
      yhalf(n+1) = se * half*channel_grid_Ly
    end do
    !
  end if
  !
  ! Create the y grid spacing for the whole channel
  !
  y2 = init_array(ncy+1,zero) ! F2003 AUTO-REALLOC
  !
  ! Get the grid points for the bottom half of the channel
  do n = 1,nhy+1
    y2(n) = yhalf(n)
  end do
  ! Get the grid points for the top half of the channel
  do n = 1,nhy
    y2(ncy+2-n) = channel_grid_Ly - yhalf(n)
  end do
  if (mypnum == glb_root) then
    do n = 1,max(size(y),size(y1),size(y2))
      val = -huge(zero)
      a = val
      b = val
      if (n <= size(y)) val = y(n)
      if (n <= size(y1)) a = y1(n)
      if (n <= size(y2)) b = y2(n)
      write (9517,'(3es20.10)') val,a,b
    end do
  end if
#endif
#endif
  !
  ! #######
  ! #######
  ! ## 3 ##
  ! #######
  ! #######
  !
  ! We can finish creating the grid now that we have the x,y,z grid spacings
  !
  ncell = ncx*ncy*ncz                       ! total number of grid cells
  nnode = (ncx+1)*(ncy+1)*(ncz+1)           ! total number of grid nodes
  nfbnd = 2*ncx*ncy + 2*ncx*ncz + 2*ncy*ncz ! total number of boundary faces
  !
  ! Allocate all the arrays needed to define the grid
  !
  xyz_nodes = init_array(nr,nnode,zero) ! F2003 AUTO-REALLOC
  nodes_of_cell_ptr = init_array(ncell+1,0) ! F2003 AUTO-REALLOC
  nodes_of_cell = init_array(geom_nodes(Geom_Hexa)*ncell,0) ! F2003 AUTO-REALLOC
  bface = init_array(nbfai,nfbnd,0) ! F2003 AUTO-REALLOC
  cell_geom = init_array(ncell,Geom_Hexa) ! F2003 AUTO-REALLOC
  cell_order = init_array(ncell,n_order) ! F2003 AUTO-REALLOC
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
  do k = 1,ncz
    do j = 1,ncy
      do i = 1,ncx
        !
        n = nodes_of_cell_ptr( map3dto1d(i,j,k,ncx,ncy) )
        !
        nodes_of_cell(n+1) = map3dto1d(i  ,j  ,k  ,ncx+1,ncy+1)
        nodes_of_cell(n+2) = map3dto1d(i+1,j  ,k  ,ncx+1,ncy+1)
        nodes_of_cell(n+3) = map3dto1d(i+1,j+1,k  ,ncx+1,ncy+1)
        nodes_of_cell(n+4) = map3dto1d(i  ,j+1,k  ,ncx+1,ncy+1)
        nodes_of_cell(n+5) = map3dto1d(i  ,j  ,k+1,ncx+1,ncy+1)
        nodes_of_cell(n+6) = map3dto1d(i+1,j  ,k+1,ncx+1,ncy+1)
        nodes_of_cell(n+7) = map3dto1d(i+1,j+1,k+1,ncx+1,ncy+1)
        nodes_of_cell(n+8) = map3dto1d(i  ,j+1,k+1,ncx+1,ncy+1)
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
  do k = 1,ncz+1
    do j = 1,ncy+1
      do i = 1,ncx+1
        !
        n = map3dto1d(i,j,k,ncx+1,ncy+1)
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
  ! Loop through the 6 grid boundaries
  ! NOTE: This corresponds to the grid boundaries in
  !       the order [imin,imax,jmin,jmax,kmin,kmax]
  !
  n2 = 0
  !
  bc1_idx = inflow_idx
  bc2_idx = outflow_idx
  !
  n1 = n2
  n2 = n2 + ncy*ncz
  do k = 1,ncz
    do j = 1,ncy
      !
      n1 = n1 + 1
      n2 = n2 + 1
      !
      bface( 1,n1) = bc_input(bc1_idx)%set_bc_value
      bface( 2,n1) = map3dto1d(1 ,j,k,ncx,ncy)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(1)
      bface( 5,n1) = bc1_idx ! group for boundary face
      bface( 6,n1) = bc1_idx ! index for bc_input needed by match_periodic_faces
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_input(bc2_idx)%set_bc_value
      bface( 2,n2) = map3dto1d(ncx,j,k,ncx,ncy)
      bface( 3,n2) = Geom_Hexa
      bface( 4,n2) = HostFace3D(2)
      bface( 5,n2) = bc2_idx ! group for boundary face
      bface( 6,n2) = bc2_idx ! index for bc_input needed by match_periodic_faces
      bface(10,n2) = ncell + n2
      !
    end do
  end do
  !
  bc1_idx = bottom_wall_idx
  bc2_idx = top_wall_idx
  !
  n1 = n2
  n2 = n2 + ncx*ncz
  do k = 1,ncz
    do i = 1,ncx
      !
      n1 = n1 + 1
      n2 = n2 + 1
      !
      bface( 1,n1) = bc_input(bc1_idx)%set_bc_value
      bface( 2,n1) = map3dto1d(i,1 ,k,ncx,ncy)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(3)
      bface( 5,n1) = bc1_idx ! group for boundary face
      bface( 6,n1) = bc1_idx ! index for bc_input
      bface(10,n1) = ncell + n1
      !
      !
      bface( 1,n2) = bc_input(bc2_idx)%set_bc_value
      bface( 2,n2) = map3dto1d(i,ncy,k,ncx,ncy)
      bface( 3,n2) = Geom_Hexa
      bface( 4,n2) = HostFace3D(4)
      bface( 5,n2) = bc2_idx ! group for boundary face
      bface( 6,n2) = bc2_idx ! index for bc_input
      bface(10,n2) = ncell + n2
      !
    end do
  end do
  !
  bc1_idx = zmin_idx
  bc2_idx = zmax_idx
  !
  n1 = n2
  n2 = n2 + ncx*ncy
  do j = 1,ncy
    do i = 1,ncx
      !
      n1 = n1 + 1
      n2 = n2 + 1
      !
      bface( 1,n1) = bc_input(bc1_idx)%set_bc_value
      bface( 2,n1) = map3dto1d(i,j,1 ,ncx,ncy)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(5)
      bface( 5,n1) = bc1_idx ! group for boundary face
      bface( 6,n1) = bc1_idx ! index for bc_input needed by match_periodic_faces
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_input(bc2_idx)%set_bc_value
      bface( 2,n2) = map3dto1d(i,j,ncz,ncx,ncy)
      bface( 3,n2) = Geom_Hexa
      bface( 4,n2) = HostFace3D(6)
      bface( 5,n2) = bc2_idx ! group for boundary face
      bface( 6,n2) = bc2_idx ! index for bc_input needed by match_periodic_faces
      bface(10,n2) = ncell + n2
      !
    end do
  end do
  !
  if (mypnum == glb_root) then
    call write_channel_grid_info(iout,ncx,ncy,ncz)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  10 format (/,80("#"),//, &
               2x,"ERROR: There is a missing boundary condition that ", &
                  "is needed",/,9x,"to run the channel flow problem:")
  11 format (12x,"The ",a9," boundary condition is missing!")
  12 format (/,80("#"),/)
  13 format ("WHEN TRYING TO INTERNALLY GENERATE THE GRID FOR THE ", &
             "CHANNEL FLOW PROBLEM, THERE IS A MISSING BOUNDARY ", &
             "CONDITION THAT IS NEEDED TO BE ABLE TO RUN THIS PROBLEM!")
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
function channel_y_spacing(yplus_wall,yplus_ave,ly_plus,inner_wall_radius) &
                    result(return_value)
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use ovar,      only : loc_solution_pts
  !
  use poly_mod,  only : legendre_gauss_quadrature
  use poly_mod,  only : legendre_gauss_lobatto_quadrature
  
  !.. Formal Arguments ..
  real(wp), intent(in) :: yplus_wall
  real(wp), intent(in) :: yplus_ave
  real(wp), intent(in) :: ly_plus
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: inner_wall_radius
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer  :: n,nhy,ncy
  real(wp) :: dl1,dl2,ds1,ds2,ymin,ymax
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: xi(:)
  real(wp), allocatable :: yhalf(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "channel_y_spacing"
  !
continue
  !
  ! Get y location of the bottom wall (inner wall radius
  ! for pipe channel) and the y location of the top wall
  !
  ymin = zero
  if (present(inner_wall_radius)) then
    ymin = inner_wall_radius
  end if
  ymax = ymin + half*channel_grid_Ly
  !
  ! 2C.0 Get the location of the 1D solution points in reference space
  !
  xi = init_array(n_order+1,zero) ! F2003 AUTO-REALLOC
  !
  if (loc_solution_pts == Legendre_Gauss_Lobatto) then
    call legendre_gauss_lobatto_quadrature(xi)
  else
    call legendre_gauss_quadrature(xi)
  end if
  !
  ! Scale xi from [-1,1] to [0,1]
  !
  xi(:) = half*(xi(:) + one)
  !
  ! 2C.1 : Get the beginning and ending spacings needed
  !        to generate the tanh spacing
  !
  ! Need to get dl1/ds1 and dl2/ds2 to create tanh spacing for y spacing
  !
  if (yplus_wall < zero) then
    !
   !! OLD VERSION: The 1st DoF off the wall is located at abs(yplus_wall)
   !!
   !! OLD VERSION:
   !! Divide the abs(yplus_wall) by the location of the 1st solution point
   !! on the interval [0,1] to get the height of the first cell off the wall.
   !!
   !dl1 = abs(yplus_wall) / xi(1)
    !
    ! NEW VERSION: The nint(abs(yplus_wall))-th DoF
    !              off the wall is located at y+=1.
    ! NEW VERSION:
    ! The location of y+=1 is at the nint(abs(yplus_wall))-th DoF off the wall,
    ! so 1/xi(nint(abs(yplus_wall))) gives the y+ spacing of the first cell off
    ! the wall.
    ! NOTE: This assumes xi is on the interval [0,1].
    ! NOTE: We need to use min/max to make sure nint(abs(yplus_wall))
    !       is a valid index for the array xi.
    !
    n = max( 1 , min( size(xi) , nint(abs(yplus_wall)) ) )
    dl1 = one / xi(n)
    !
  else
    !
    ! The height of the first cell off the wall is yplus_wall
    !
    dl1 = yplus_wall
    !
  end if
  !
  ! Get the centerline y+ spacing
  !
  dl2 = yplus_ave
  !
  ! Map the beginning and ending grid cell spacings to the reference
  ! interval [0,1]. We are only worried about half the channel so we
  ! need to half of the channel height in plus coordinates.
  !
  ds1 = dl1 / (half*ly_plus)
  ds2 = dl2 / (half*ly_plus)
  !
  ! Compute the approximate number of cells (for half of the channel height)
  ! needed to get a constant stretching with end spacings ds1 and ds2.
  !
  nhy = nint( log(ds2/ds1) / log((ds1-one)/(ds2-one)) )
  !
  ! The total number of y cells is just twice nhy,
  ! since nhy is for half of the channel height
  !
  ncy = nhy*2
  !
  ! Save all the various y+ spacings to the module variables
  !
  yplus_cell_wall = dl1
  yplus_dof_wall = xi(:)*dl1 ! F2003 AUTO-REALLOC
  yplus_cell_centerline = dl2
  yplus_cell_ave = ly_plus / real(ncy,kind=wp)
  !
  ! The tanh_spacing functions require the end spacings in the same scale
  ! as the end coordinates (the spacings and coordinates will be scaled
  ! to the interval [0,1] inside the function), so unscale ds1 and ds2 from
  ! the interval [0,1] to the actual scale of the end coordinates.
  !
  ds1 = ds1 * (ymax-ymin)
  ds2 = ds2 * (ymax-ymin)
  !
  ! Get the y grid spacing for half the channel using a tanh spacing
  ! NOTE: The optional argument fixed_end_spacings determines if the tanh
  !       solver is either
  !         1) (when true) enforcing the end spacings ds1 and ds2 and solving
  !            for the np-3 interior spacings between points 2 and np-1, or
  !         2) (when false) not strictly enforcing the end spacings ds1 and ds2
  !            and solving for the np-1 interior spacings between the end points
  !
  yhalf = tanh_spacing(ymin,ymax,ds1,ds2,nhy+1, & ! F2003 AUTO-REALLOC
                       fixed_end_spacings=fals) 
  !
  ! Create the y grid spacing for the whole channel
  !
  return_value = init_array(ncy+1,zero) ! F2003 AUTO-REALLOC
  !
  ! Get the grid points for the bottom half of the channel
  do n = 1,nhy+1
    return_value(n) = yhalf(n)
  end do
  ! Get the grid points for the top half of the channel
  do n = 1,nhy
    return_value(ncy+2-n) = two*ymin + channel_grid_Ly - yhalf(n)
  end do
  !
end function channel_y_spacing
!
!###############################################################################
!
function old_channel_y_spacing(yplus_wall,yplus_ave,ly_plus,inner_wall_radius) &
                        result(return_value)
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use ovar,      only : loc_solution_pts
  !
  use poly_mod,  only : legendre_gauss_quadrature
  use poly_mod,  only : legendre_gauss_lobatto_quadrature
  
  !.. Formal Arguments ..
  real(wp), intent(in) :: yplus_wall
  real(wp), intent(in) :: yplus_ave
  real(wp), intent(in) :: ly_plus
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: inner_wall_radius
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer  :: n,nhy,ncy,i
  real(wp) :: dl1,dl2,ds1,ds2,a,b,sf,rnum,rden,sferr
  real(wp) :: val,ue,se,ymin,ymax
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: xi(:)
  real(wp), allocatable :: yhalf(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "old_channel_y_spacing"
  !
continue
  !
  ! Get y location of the bottom wall (inner wall radius
  ! for pipe channel) and the y location of the top wall
  !
  ymin = zero
  if (present(inner_wall_radius)) then
    ymin = inner_wall_radius
  end if
  ymax = ymin + half*channel_grid_Ly
  !
  ! 2C.0 Get the location of the 1D solution points in reference space
  !
  xi = init_array(n_order+1,zero) ! F2003 AUTO-REALLOC
  !
  if (loc_solution_pts == Legendre_Gauss_Lobatto) then
    call legendre_gauss_lobatto_quadrature(xi)
  else
    call legendre_gauss_quadrature(xi)
  end if
  !
  ! Scale xi from [-1,1] to [0,1]
  !
  xi(:) = half*(xi(:) + one)
  !
  ! 2C.1 : Get the beginning and ending spacings needed
  !        to generate the tanh spacing
  !
  ! Need to get dl1/ds1 and dl2/ds2 to create tanh spacing for y spacing
  !
  if (yplus_wall < zero) then
    !
   !! OLD VERSION: The 1st DoF off the wall is located at abs(yplus_wall)
   !!
   !! OLD VERSION:
   !! Divide the abs(yplus_wall) by the location of the 1st solution point
   !! on the interval [0,1] to get the height of the first cell off the wall.
   !!
   !dl1 = abs(yplus_wall) / xi(1)
    !
    ! NEW VERSION: The nint(abs(yplus_wall))-th DoF
    !              off the wall is located at y+=1.
    ! NEW VERSION:
    ! The location of y+=1 is at the nint(abs(yplus_wall))-th DoF off the wall,
    ! so 1/xi(nint(abs(yplus_wall))) gives the y+ spacing of the first cell off
    ! the wall.
    ! NOTE: This assumes xi is on the interval [0,1].
    ! NOTE: We need to use min/max to make sure nint(abs(yplus_wall))
    !       is a valid index for the array xi.
    !
    n = max( 1 , min( size(xi) , nint(abs(yplus_wall)) ) )
    dl1 = one / xi(n)
    !
  else
    !
    ! The height of the first cell off the wall is yplus_wall
    !
    dl1 = yplus_wall
    !
  end if
  !
  ! Get the centerline y+ spacing
  !
  dl2 = yplus_ave
  !
  ! Map the beginning and ending grid cell spacings to the reference
  ! interval [0,1]. We are only worried about half the channel so we
  ! need to half of the channel height in plus coordinates.
  !
  ds1 = dl1 / (half*ly_plus)
  ds2 = dl2 / (half*ly_plus)
  !
  ! Compute the approximate number of cells (for half of the channel height)
  ! needed to get a constant stretching with end spacings ds1 and ds2.
  !
  nhy = nint( log(ds2/ds1) / log((ds1-one)/(ds2-one)) )
  !
  ! The total number of y cells is just twice nhy,
  ! since nhy is for half of the channel height
  !
  ncy = nhy*2
  !
  ! Save all the various y+ spacings to the module variables
  !
  yplus_cell_wall = dl1
  yplus_dof_wall = xi(:)*dl1 ! F2003 AUTO-REALLOC
  yplus_cell_centerline = dl2
  yplus_cell_ave = ly_plus / real(ncy,kind=wp)
  !
  ! Initialize the y grid spacing for half the channel
  !
  yhalf = init_array(nhy+1,zero) ! F2003 AUTO-REALLOC
  !
  ! 2C.2 : Get the a and b parameters to begin solving for the tanh y-spacing
  !
  a = sqrt(ds2/ds1)
  b = one / (real(nhy,kind=wp)*sqrt(ds1*ds2))
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
    ! Get the distribution yhalf(n)
    !
    do n = 0,nhy
      val = real(n,kind=wp)/real(nhy,kind=wp)
      ue = half*(one + tanh(sf*(val-half))/tanh(sf*half))
      se = ue / (a + (one-a)*ue)
      yhalf(n+1) = se  * half*channel_grid_Ly
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
    ! Get the distribution yhalf(n)
    !
    do n = 0,nhy
      val = real(n,kind=wp)/real(nhy,kind=wp)
      ue = half*(one + tan(sf*(val-half))/tan(sf*half))
      se = ue / (a + (one-a)*ue)
      yhalf(n+1) = se  * half*channel_grid_Ly
    end do
    !
  else
    !
    ! Both of the above Newton iterations can have numerical difficulties
    ! if the stretching factor is close to 1, therefore use this simpler
    ! stretching function that doesnt require a Newton iteration.
    !
    do n = 0,nhy
      val = real(n,kind=wp)/real(nhy,kind=wp)
      ue = val*(one + two*(b-one)*(val-half)*(one-val))
      se = ue / (a + (one-a)*ue)
      yhalf(n+1) = se  * half*channel_grid_Ly
    end do
    !
  end if
  !
  ! Create the y grid spacing for the whole channel
  !
  return_value = init_array(ncy+1,zero) ! F2003 AUTO-REALLOC
  !
  ! Get the grid points for the bottom half of the channel
  do n = 1,nhy+1
    return_value(n) = yhalf(n)
  end do
  ! Get the grid points for the top half of the channel
  do n = 1,nhy
    return_value(ncy+2-n) = two*ymin + channel_grid_Ly - yhalf(n)
  end do
  !
end function old_channel_y_spacing
!
!###############################################################################
!
subroutine write_channel_grid_info(iunit,ncx,ncy,ncz)
  !
  !.. Use Statements ..
  use order_mod,   only : n_order
  use ovar,        only : yplus1 => channel_yplus_of_one
  use geomvar_mod, only : n_global_cell
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: iunit,ncx,ncy,ncz
  !
  !.. Local Scalars ..
  integer  :: i,io,icx,icy,icz
  real(wp) :: pp1
  real(wp) :: xplus_dof_ave
  real(wp) :: yplus_dof_ave
  real(wp) :: zplus_dof_ave
  real(wp) :: yplus_dof_centerline
  !
continue
  !
  io = iout
  if (present(iunit)) io = iunit
  !
  pp1 = real(n_order+1,kind=wp)
  !
  xplus_dof_ave = xplus_cell_ave / pp1
  yplus_dof_ave = yplus_cell_ave / pp1
  zplus_dof_ave = zplus_cell_ave / pp1
  yplus_dof_centerline = yplus_cell_centerline / pp1
  !
  if (present(ncx)) then
    icx = ncx
  else
    icx = nint( channel_grid_Lx / (xplus_cell_ave*yplus1) )
  end if
  !
  if (present(ncz)) then
    icz = ncz
  else
    icz = nint( channel_grid_Lz / (zplus_cell_ave*yplus1) )
  end if
  !
  if (present(ncy)) then
    icy = ncy
  else
    icy = n_global_cell / (icx*icz)
  end if
  !
  write (io,1)
  write (io,2) icx,icy,icz
  write (io,3) &
    "dx+",xplus_cell_ave,xplus_dof_ave, &
    "(centerline) dy+",yplus_cell_centerline,yplus_dof_centerline, &
    "(grid) dy+",yplus_cell_ave,yplus_dof_ave, &
    "dz+",zplus_cell_ave,zplus_dof_ave
  !
  if (allocated(yplus_dof_wall)) then
    write (io,4) yplus_cell_wall
    do i = 1,size(yplus_dof_wall)
      write (io,5) i,yplus_dof_wall(i)
    end do
  else
    write (io,6) yplus_cell_wall
  end if
  write (io,1)
  !
  ! Format Statements
  !
  1 format (/,80("="),/)
  2 format (2x,"Details of internally generated grid used for channel flow",/, &
            80("-"),//, &
            4x,"Cells [nx,ny,nz] = [",i5,",",i5,",",i5,"]")
  3 format (/,26x,"Cell Ave",7x,"DoF Ave",/, &
            4(a22,f12.2,2x,f12.2,/))
  4 format (4x,"y+ spacing  for first cell off wall",/,24x,f10.5,//, &
            4x,"y+ spacings for DoF in first cell off wall")
  5 format (14x,"xi+(",i2,") = ",f10.5)
  6 format (4x,"y+ spacing  for first cell off wall",/,24x,f10.5)
  !
end subroutine write_channel_grid_info
!
!###############################################################################
!
subroutine create_pipe_channel_grid()
  !
  ! This subroutine internally generates a grid for the pipe channel flow (PCF)
  ! problem. It is sort a combination of the grid generators for the TCF and
  ! channel flow problems, so it utilizes input variables from both problems.
  !
  !.. Use Statements ..
  use order_mod,   only : n_order
  use ovar,        only : Lref,bc_input
  use geomvar_mod, only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geomvar_mod, only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geomvar_mod, only : grid,cell_geom,cell_order
  use ovar,        only : grid_format,bc_names,internal_grid_order
  !
  use ovar, only : channel_grid_ave_xplus
  use ovar, only : channel_grid_ave_yplus
  use ovar, only : channel_grid_ave_zplus
  use ovar, only : channel_grid_yplus_wall
  use ovar, only : channel_yplus_of_one
  use ovar, only : pipe_channel_angle
  !
  use cgns_types_mod, only : HEXA_8,HEXA_27,HEXA_64,HEXA_125
  !
  use cgns_mod, only : get_face_elem
  !
  !.. Local Scalars ..
  integer  :: ierr,i,j,k
  integer  :: nc,nf,n1,n2,bc1_idx,bc2_idx
  integer  :: ncx,ncy,ncz,npts
  integer  :: host_cell,host_side,elem_type
  integer  :: inner_wall_idx,outer_wall_idx
  integer  :: inflow_idx,outflow_idx
  integer  :: theta_min_idx,theta_max_idx
  integer  :: grid_order
  real(wp) :: lx,ly,lz,angle,ri,ro
  real(wp) :: xplus_ave,yplus_ave,zplus_ave,dx,dz,yplus1
  real(wp) :: lx_plus,ly_plus,lz_plus,yplus_wall,pp1
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: elem_nodes(:,:)
  integer,  allocatable :: elem_pts(:,:)
  real(wp), allocatable :: y(:)
  real(wp), allocatable :: xyz_pts(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_pipe_channel_grid"
  !
 !logical(lk), parameter :: lz_is_inner_wall = true
  logical(lk), parameter :: lz_is_inner_wall = fals
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
  ! The Taylor-Couette flow problem is always 3D
  !
  nr = 3
  !
  ! Get the length of the cylinders
  !
  lx = channel_grid_Lx
  ly = channel_grid_Ly
  lz = channel_grid_Lz
  !
  angle = min(nine*ten,max(eps1,pipe_channel_angle))
  !
  ! Get the inner and outer radii for the pipe walls
  !
  if (lz_is_inner_wall) then
    ri = lz/(angle*degrees2radians)
    ro = ri + ly
  else
    ro = lz/(angle*degrees2radians)
    ri = ro - ly
  end if
  !
  ! Reset bc_input for this axially (and possibly rotationally) periodic grid
  ! so that the subroutine match_periodic_faces can correctly match all the
  ! periodic faces.
  !
  call reset_bc_input_for_tcf_grid(lx,angle,inner_wall_bc=bc_adiabatic_wall)
  !
  ! Find the index in bc_input for the 4 (maybe 6) different boundary conditions
  !
  inner_wall_idx = -1
  outer_wall_idx = -1
  inflow_idx = -1
  outflow_idx = -1
  theta_min_idx = -1
  theta_max_idx = -1
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
    else if (trim(adjustl(uppercase(bc_names(i)))) == "THETAMIN") then
      theta_min_idx = i
    else if (trim(adjustl(uppercase(bc_names(i)))) == "THETAMAX") then
      theta_max_idx = i
    end if
  end do
  !
  ! If any of the above BC indices are still -1, there is a missing BC in
  ! bc_input needed for the TCF problem, so write out an error for that
  ! missing BC and stop the simulation.
  !
  if (any([inner_wall_idx,outer_wall_idx,inflow_idx, &
           outflow_idx,theta_min_idx,theta_max_idx] == -1)) then
    if (mypnum == glb_root) then
      write (iout,10)
      if (inner_wall_idx == -1) write (iout,11) "INNERWALL"
      if (outer_wall_idx == -1) write (iout,11) "OUTERWALL"
      if (inflow_idx == -1)     write (iout,11) "INFLOW"
      if (outflow_idx == -1)    write (iout,11) "OUTFLOW"
      if (theta_min_idx == -1)  write (iout,11) "THETAMIN"
      if (theta_max_idx == -1)  write (iout,11) "THETAMAX"
      write (iout,12)
    end if
    flush (iout)
    write (error_message,13)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Get the grid order if specified for an internally
  ! generated grid, otherwise set it to the solution order.
  !
  grid_order = merge( internal_grid_order , n_order , internal_grid_order>0 )
  !
  ! For now we are limiting the grid order to
  ! the max order allowed by CGNS, which is 4
  !
  grid_order = min(grid_order,4)
  !
  ! Precompute P+1 as a real
  !
  pp1 = real(n_order+1,kind=wp)
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
  if (channel_grid_ave_xplus  == zero) channel_grid_ave_xplus = -three*ten
  if (channel_grid_ave_yplus  == zero) channel_grid_ave_yplus = -three*five
  if (channel_grid_ave_zplus  == zero) channel_grid_ave_zplus = -three*five
  if (channel_grid_yplus_wall == zero) channel_grid_yplus_wall = two
  !
  ! Save the channel grid parameters to local variables
  !
  yplus1 = channel_yplus_of_one
  xplus_ave = channel_grid_ave_xplus
  yplus_ave = channel_grid_ave_yplus
  zplus_ave = channel_grid_ave_zplus
  yplus_wall = channel_grid_yplus_wall
  !
  ! Get the length of the domain in each direction in y+ dimensions
  !
  lx_plus = lx / yplus1
  ly_plus = ly / yplus1
  lz_plus = lz / yplus1
  !
  ! Compute the cell average x+,y+,z+ for the grid
  ! NOTE: We are assuming the y+ cell average is the y+ spacing in the middle
  !       of the channel and not the actual average y+ spacing including the
  !       near wall cells. This will tend towards finer y+ spacing in the
  !       middle of the channel, but it makes the input more directly
  !       comparable to the x+ and z+ inputs. Otherwise, if the average y+
  !       spacing includes the near wall cells in the average, the spacing
  !       in the middle of the channel will be larger than the equivalent
  !       x+ or z+ spacing given by the respective inputs.
  !
  ! If the inputs channel_grid_ave_[xyz]plus are:
  !   > 0 : the input specifies the average cell spacing
  !   < 0 : the abs(input) specifies the average DoF spacing
  !
  ! For a given coordinate direction (say [i]), the cell average [i]+ is
  !   if (channel_grid_ave_[i]plus > zero) then
  !     [i]plus_ave = channel_grid_ave_[i]plus
  !   else
  !     [i]plus_ave = abs(channel_grid_ave_[i]plus)*real(n_order+1,kind=wp)
  !   end if
  ! NOTE: We use (P+1) here because the delta plus value is the
  !       average distance between the DoF within a cell (or the
  !       subcells), and there are P+1 subcells in each cell (the
  !       subcells between the end points and the cell boundaries
  !       are kind of 1/2 a subcell each).
  !
  xplus_ave = abs(xplus_ave)*merge(one , pp1 , xplus_ave > zero)
  yplus_ave = abs(yplus_ave)*merge(one , pp1 , yplus_ave > zero)
  zplus_ave = abs(zplus_ave)*merge(one , pp1 , zplus_ave > zero)
  !
  ! Get the number of cells in the x and z directions
  !
  ncx = nint( lx_plus / xplus_ave )
  ncz = nint( lz_plus / zplus_ave )
  !
  ! Compute the actual (non-plus) x and z cell spacings
  ! NOTE: These will most likely not be equal to xplus_ave and zplus_ave
  !       since lx_plus and lz_plus will probably not be evenly divisible
  !       by these respective numbers.
  !
  dx = channel_grid_Lx / real(ncx,kind=wp)
  dz = channel_grid_Lz / real(ncz,kind=wp)
  
  !
  ! Save the x+ and z+ spacings to the module variables
  !
  xplus_cell_ave = dx / yplus1
  zplus_cell_ave = dz / yplus1
  !
  ! #######
  ! #######
  ! ## 2 ##
  ! #######
  ! #######
  !
  ! Get the grid point locations for y
  !
  y = channel_y_spacing(yplus_wall,yplus_ave,ly_plus, & ! F2003 AUTO-REALLOC
                        inner_wall_radius=ri)
  ncy = size(y)-1
  !
  ! Create the linear grid cell connectivities (elem_nodes), the high-order
  ! grid cell connectivities (elem_pts), the linear grid node coordinates
  ! (xyz_nodes), and the high-order grid point coordinates (xyz_pts).
  !
  call high_order_TCF_grid_points(ncx,ncy,ncz,lx,angle,y,grid_order, &
                                  elem_nodes,elem_pts,xyz_nodes,xyz_pts)
  !
  ! Deallocate y since it is no longer needed
  !
  if (allocated(y)) then
    deallocate ( y , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"y",2,__LINE__,__FILE__,ierr,error_message)
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
  ncell = ncx*ncy*ncz                       ! total number of grid cells
  nnode = (ncx+1)*(ncy+1)*(ncz+1)           ! total number of grid nodes
  nfbnd = 2*ncx*ncy + 2*ncx*ncz + 2*ncy*ncz ! total number of boundary faces
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
  ! Create the nodes_of_cell and nodes_of_cell_ptr arrays. nodes_of_cell is
  ! just a 1D version of elem_nodes, and nodes_of_cell_ptr is an indexing
  ! array that gives the beginning/ending index in nodes_of_cell for each
  ! cell.
  !
  ! Create the nodes_of_cell_ptr array
  !
  nodes_of_cell_ptr = intseq(0,ncell)*geom_nodes(Geom_Hexa) ! F2003 AUTO-REALLOC
  !
  ! Create the nodes_of_cell array
  !
  ! F2003 AUTO-REALLOC
  nodes_of_cell = reshape( elem_nodes(:,:) , shape=[size(elem_nodes)] )
  !
  ! #######
  ! #######
  ! ## 4 ##
  ! #######
  ! #######
  !
  ! Create the bface array
  !
  ! Loop through the 6 grid boundaries
  ! NOTE: This corresponds to the grid boundaries in
  !       the order [InnerWall,OuterWall,Inflow,Outflow,ThetaMin,ThetaMax] which
  !       correspond to the [jmin,jmax,imin,imax,kmin,kmax] boundaries.
  !
  ! InnerWall/j-min and OuterWall/j-max boundaries
  !
  n2 = 0
  !
  bc1_idx = inner_wall_idx
  bc2_idx = outer_wall_idx
  !
  n1 = n2           ! j-min (r-min) counter
  n2 = n2 + ncx*ncz ! j-max (r-max) counter
  do k = 1,ncz
    do i = 1,ncx
      !
      n1 = n1 + 1 ! j-min (r-min) face
      n2 = n2 + 1 ! j-max (r-max) face
      !
      bface( 1,n1) = bc_input(bc1_idx)%set_bc_value
      bface( 2,n1) = map3dto1d(i,  1,k,ncx,ncy,ncz)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(3)
      bface( 5,n1) = bc1_idx ! group for boundary face
      bface( 6,n1) = bc1_idx ! index for bc_input for this BC
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_input(bc2_idx)%set_bc_value
      bface( 2,n2) = map3dto1d(i,ncy,k,ncx,ncy,ncz)
      bface( 3,n2) = Geom_Hexa
      bface( 4,n2) = HostFace3D(4)
      bface( 5,n2) = bc2_idx ! group for boundary face
      bface( 6,n2) = bc2_idx ! index for bc_input for this BC
      bface(10,n2) = ncell + n2
      !
    end do
  end do
  !
  ! Min-X and Max-X boundary faces
  !
  bc1_idx = inflow_idx
  bc2_idx = outflow_idx
  !
  n1 = n2           ! min-x counter
  n2 = n2 + ncy*ncz ! max-x counter
  do k = 1,ncz
    do j = 1,ncy
      !
      n1 = n1 + 1 ! x-min face
      n2 = n2 + 1 ! x-max face
      !
      ! NOTE: The subroutine match_periodic_faces will set the values of
      !       bface([5,6,8,9],:), so we dont need to do that here.
      !
      bface( 1,n1) = bc_input(bc1_idx)%set_bc_value
      bface( 2,n1) = map3dto1d(  1,j,k,ncx,ncy,ncz)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(1)
      bface( 5,n1) = bc1_idx ! group for boundary face
      bface( 6,n1) = bc1_idx ! index for bc_input needed by match_periodic_faces
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_input(bc2_idx)%set_bc_value
      bface( 2,n2) = map3dto1d(ncx,j,k,ncx,ncy,ncz)
      bface( 3,n2) = Geom_Hexa
      bface( 4,n2) = HostFace3D(2)
      bface( 5,n2) = bc2_idx ! group for boundary face
      bface( 6,n2) = bc2_idx ! index for bc_input needed by match_periodic_faces
      bface(10,n2) = ncell + n2
      !
    end do
  end do
  !
  ! Min-Theta/K and Max-Theta/K rotationally periodic boundary faces
  !
  bc1_idx = theta_min_idx
  bc2_idx = theta_max_idx
  !
  n1 = n2           ! min-theta/k counter
  n2 = n2 + ncx*ncy ! max-theta/k counter
  do j = 1,ncy
    do i = 1,ncx
      !
      n1 = n1 + 1 ! x-min face
      n2 = n2 + 1 ! x-max face
      !
      ! NOTE: The subroutine match_periodic_faces will set the values of
      !       bface([5,6,8,9],:), so we dont need to do that here.
      !
      bface( 1,n1) = bc_input(bc1_idx)%set_bc_value
      bface( 2,n1) = map3dto1d(i,j,  1,ncx,ncy,ncz)
      bface( 3,n1) = Geom_Hexa
      bface( 4,n1) = HostFace3D(5)
      bface( 5,n1) = bc1_idx ! group for boundary face
      bface( 6,n1) = bc1_idx ! index for bc_input needed by match_periodic_faces
      bface(10,n1) = ncell + n1
      !
      bface( 1,n2) = bc_input(bc2_idx)%set_bc_value
      bface( 2,n2) = map3dto1d(i,j,ncz,ncx,ncy,ncz)
      bface( 3,n2) = Geom_Hexa
      bface( 4,n2) = HostFace3D(6)
      bface( 5,n2) = bc2_idx ! group for boundary face
      bface( 6,n2) = bc2_idx ! index for bc_input needed by match_periodic_faces
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
  ! Only have the host roots create the grid derived type. The high-order
  ! data it will contain can take up a lot of memory. Only having the root
  ! process on each host/node store this data will save a lot of memory
  ! instead of having [# of cores per node] copies of the same data on each
  ! node. The partitioning code only expects the host root to have this
  ! information, so storing this information on non-root processes is just
  ! a waste of memory anyways.
  !
  if (i_am_host_root) then
    !
    allocate ( grid , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( grid%elem(1:ncell+nfbnd) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid%elem",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Create grid for the grid cells
    !
    if (grid_order == 1) then
      elem_type = int(HEXA_8)
      npts = 8
    else if (grid_order == 2) then
      elem_type = int(HEXA_27)
      npts = 27
    else if (grid_order == 3) then
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
      grid%elem(nc)%prop%order = grid_order
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
      !            cells %pts to put it back in CGNS order, then return the
      !            face elements %pts array in GFR ordering
      !
      if (grid%elem(host_cell)%prop%order <= 4) then
        grid%elem(nc) = get_face_elem(grid%elem(host_cell),host_side,host_cell)
      else
       !grid%elem(nc) = get_face_elem(grid%elem(host_cell),host_side,host_cell)
        write (error_message,1)
        call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
      end if
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
  end if
  !
 !if (mypnum == glb_root) then
 !  write (iout,2) minval( grid%elem(:)%prop%order ), &
 !                 maxval( grid%elem(:)%prop%order ), &
 !                 n_order
 !end if
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
  ! Format Statements
  !
  1 format ("ERROR: The internal grid generator for creating the ", &
                   "Taylor-Couette flow grid is currently limited to ", &
                   "creating P1-P4 order grids.")
  !
  2 format (/,4x,"After creating the Taylor-Couette flow grid:",/, &
             10x,"Min Grid Order = ",i0,/, &
             10x,"Max Grid Order = ",i0,/, &
             10x,"Solution Order = ",i0,/)
  !
  10 format (/,80("#"),//, &
               2x,"ERROR: There is a missing boundary condition that ", &
                  "is needed",/,9x,"to run the Taylor-Couette problem:")
  11 format (12x,"The ",a9," boundary condition is missing!")
  12 format (/,80("#"),/)
  13 format ("WHEN TRYING TO INTERNALLY GENERATE THE GRID FOR THE ", &
             "TAYLOR-COUETTE FLOW PROBLEM, THERE IS A MISSING BOUNDARY ", &
             "CONDITION THAT IS NEEDED TO BE ABLE TO RUN THIS PROBLEM!")
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
end subroutine create_pipe_channel_grid
!
!###############################################################################
!
subroutine create_flame_speed_test_grid()
  !
  !.. Use Statements ..
  use order_mod,   only : n_order
  use ovar,        only : Lref,itestcase
  use ovar,        only : flame_speed_test_dof
  use geomvar_mod, only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geomvar_mod, only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geomvar_mod, only : cell_geom,cell_order
  use ovar,        only : grid_format,bc_names
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
subroutine create_sod_test_grid()
  !
  !.. Use Statements ..
  use order_mod,   only : n_order
  use ovar,        only : Lref,itestcase
  use ovar,        only : sod_test_dof
  use geomvar_mod, only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geomvar_mod, only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geomvar_mod, only : cell_geom,cell_order
  use ovar,        only : grid_format,bc_names
  !
  !.. Local Scalars ..
  integer  :: n,nc,nf,np,ierr
  real(wp) :: dl
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_sod_test_grid"
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
  ! Make the reference length one meter for this test case
  !
  Lref = one
  !
  ! The Sod test case is always 1D
  !
  nr = 1
  !
  ! Estimate the number of cells required to get the requested DOF
  !
  if ( sod_test_dof > 0 ) then
    nc = nint( real(sod_test_dof,kind=wp) / real(n_order+1,kind=wp) )
  else
    nc = abs(sod_test_dof)  ! directly specify number of cells
  endif
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
! bface( 1,nf) = bc_slip_wall
  bface( 1,nf) = bc_sup_outflow
  bface( 2,nf) = 1      ! host cell index
  bface( 3,nf) = Geom_Edge
  bface( 4,nf) = 1      ! x-min face
  bface( 5,nf) = 1      ! this face belongs to boundary group 1
  bface( 6,nf) = 0      ! index of array bc_in(:) to get BC flow conditions
  bface( 7,nf) = 0      ! not used for physical boundary
  bface( 8,nf) = 0      ! not used for physical boundary
  bface( 9,nf) = 0      ! not used for physical boundary
  bface(10,nf) = ncell + nf ! index of ghost cell/face
  !
  ! Create the i-max boundary face
  !
  nf = 2
  !
! bface( 1,nf) = bc_slip_wall
  bface( 1,nf) = bc_sup_outflow
  bface( 2,nf) = ncell  ! host cell index
  bface( 3,nf) = Geom_Edge
  bface( 4,nf) = 2      ! x-max face
  bface( 5,nf) = 2      ! this face belongs to boundary group 2
  bface( 6,nf) = 0      ! index of array bc_in(:) to get BC flow conditions
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
! bc_names(1) = "SLIPWALL"
! bc_names(2) = "SLIPWALL"
  bc_names(1) = "SUP_OUTFLOW"
  bc_names(2) = "SUP_OUTFLOW"
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine create_sod_test_grid
!
!###############################################################################
!
subroutine create_pulse_test_grid()
  !
  !.. Use Statements ..
  use order_mod,   only : n_order
  use ovar,        only : Lref,itestcase
  use ovar,        only : pulse_test_dof, pulse_test_length
  use geomvar_mod, only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geomvar_mod, only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geomvar_mod, only : cell_geom,cell_order
  use ovar,        only : grid_format,bc_names
  !
  !.. Local Scalars ..
  integer  :: n,nc,nf,np,ierr
  real(wp) :: dl
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_pulse_test_grid"
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
  ! Make the reference length one meter for this test case
  !
  Lref = one
  !
  ! The pulse test case is always 1D
  !
  nr = 1
  !
  ! Reset bc_input for this fully periodic grid so that the subroutine
  ! match_periodic_faces can correctly match all the periodic faces.
  !
  call reset_bc_input_for_internal_grid( Lref_x = pulse_test_length )
  !
  ! Estimate the number of cells required to get the requested DOF
  !
  if ( pulse_test_dof > 0 ) then
    nc = nint( real(pulse_test_dof,kind=wp) / real(n_order+1,kind=wp) )
  else
    nc = abs(pulse_test_dof)  ! directly specify number of cells
  endif
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
  ! Compute the cell size for a Cartesian/uniform grid of total length
  !
  dl = pulse_test_length / real(ncell,kind=wp)
  !
  ! Create the grid node coordinates as a uniform 1D grid from 0:length
  !
  do n = 1,nnode
    xyz_nodes(1,n) = dl*real(n-1,kind=wp)
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
  ! NOTE: The subroutine match_periodic_faces will set the values of
  !       bface([5,6,8,9],:), so we dont need to do that here.
  !
  !
  ! Create the i-min boundary face
  !
  nf = 1
  !
  bface( 1,nf) = bc_periodic
  bface( 2,nf) = 1      ! host cell index
  bface( 3,nf) = Geom_Edge
  bface( 4,nf) = 1      ! x-min face
! bface( 5,nf) = 1      ! this face belongs to boundary group 1
  bface( 6,nf) = 1      ! index of array bc_input(:) to get BC flow conditions
  bface(10,nf) = ncell + nf ! index of ghost cell/face
  !
  ! Create the i-max boundary face
  !
  nf = 2
  !
  bface( 1,nf) = bc_periodic
  bface( 2,nf) = ncell  ! host cell index
  bface( 3,nf) = Geom_Edge
  bface( 4,nf) = 2      ! x-max face
! bface( 5,nf) = 2      ! this face belongs to boundary group 2
  bface( 6,nf) = 2      ! index of array bc_input(:) to get BC flow conditions
  bface(10,nf) = ncell + nf ! index of ghost cell/face
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine create_pulse_test_grid
!
!###############################################################################
!
subroutine create_shu_osher_test_grid()
  !
  !.. Use Statements ..
  use order_mod,   only : n_order
  use ovar,        only : Lref,itestcase
  use ovar,        only : shu_osher_test_dof
  use geomvar_mod, only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geomvar_mod, only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geomvar_mod, only : cell_geom,cell_order
  use ovar,        only : grid_format,bc_names
  !
  !.. Local Scalars ..
  integer  :: n,nc,nf,np,ierr
  integer  :: inflow_idx,outflow_idx
  real(wp) :: dl
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_shu_osher_test_grid"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Reset bc_input for the 1D Shu-Osher problem
  !
  call reset_bc_input_for_1d_grid(inflow=bc_sup_inflow)
  !
  ! Find the inflow and outflow indices in bc_input
  ! NOTE: The subroutine reset_bc_input_for_1d_grid should force the
  !       names of the boundary conditions to INFLOW and OUTFLOW.
  !
  inflow_idx = -1
  outflow_idx = -1
  do n = 1,size(bc_names)
    if (trim(adjustl(uppercase(bc_names(n)))) == "INFLOW") then
      inflow_idx = n
    else if (trim(adjustl(uppercase(bc_names(n)))) == "OUTFLOW") then
      outflow_idx = n
    end if
  end do
  !
  ! If any of the above BC indices are still -1, there is a missing BC in
  ! bc_input needed for the 1D Shu-Osher problem, so write out an error for
  ! that missing BC and stop the simulation.
  !
  if (any([inflow_idx,outflow_idx] == -1)) then
    if (mypnum == glb_root) then
      write (iout,1)
      if (inflow_idx == -1)  write (iout,2) "INFLOW"
      if (outflow_idx == -1) write (iout,2) "OUTFLOW"
      write (iout,3)
    end if
    flush (iout)
    write (error_message,4)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Reset the grid format input parameter to identify
  ! the grid as internally generated
  !
  grid_format = Internally_Generated_Grid
  !
  ! Make the reference length the domain width in meters for this test case
  !
  Lref = ten
  !
  ! The Sod test case is always 1D
  !
  nr = 1
  !
  ! Estimate the number of cells required to get the requested DOF
  !
  if ( shu_osher_test_dof > 0 ) then
    nc = nint( real(shu_osher_test_dof,kind=wp) / real(n_order+1,kind=wp) )
  else
    nc = abs(shu_osher_test_dof)  ! directly specify number of cells
  endif
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
  bface( 1,nf) = inflow_idx
  bface( 2,nf) = 1      ! host cell index
  bface( 3,nf) = Geom_Edge
  bface( 4,nf) = 1      ! x-min face
  bface( 5,nf) = 1      ! this face belongs to boundary group 1
  bface( 6,nf) = 0      ! index of array bc_in(:) to get BC flow conditions
  bface( 7,nf) = 0      ! not used for physical boundary
  bface( 8,nf) = 0      ! not used for physical boundary
  bface( 9,nf) = 0      ! not used for physical boundary
  bface(10,nf) = ncell + nf ! index of ghost cell/face
  !
  ! Create the i-max boundary face
  !
  nf = 2
  !
  bface( 1,nf) = outflow_idx
  bface( 2,nf) = ncell  ! host cell index
  bface( 3,nf) = Geom_Edge
  bface( 4,nf) = 2      ! x-max face
  bface( 5,nf) = 2      ! this face belongs to boundary group 2
  bface( 6,nf) = 0      ! index of array bc_in(:) to get BC flow conditions
  bface( 7,nf) = 0      ! not used for physical boundary
  bface( 8,nf) = 0      ! not used for physical boundary
  bface( 9,nf) = 0      ! not used for physical boundary
  bface(10,nf) = ncell + nf ! index of ghost cell/face
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/,80("#"),//, &
              2x,"ERROR: There is a missing boundary condition that ", &
                 "is needed",/,9x,"to run the 1D Shu-Osher problem:")
  2 format (12x,"The ",a9," boundary condition is missing!")
  3 format (/,80("#"),/)
  4 format ("WHEN TRYING TO INTERNALLY GENERATE THE GRID FOR THE ", &
            "1D SHU-OSHER FLOW PROBLEM, THERE IS A MISSING BOUNDARY ", &
            "CONDITION THAT IS NEEDED TO BE ABLE TO RUN THIS PROBLEM!")
  !
end subroutine create_shu_osher_test_grid
!
!###############################################################################
!
subroutine reset_bc_input_for_1d_grid(periodic_distance,inflow,outflow)
  !
  !.. Use Statements ..
  use ovar, only : bc_input,bc_input_t,periodic_bc_t
  use ovar, only : bc_names
  !
  !.. Optional Arguments ..
  real(wp), optional, intent(in) :: periodic_distance
  integer,  optional, intent(in) :: inflow
  integer,  optional, intent(in) :: outflow
  !
  !.. Local Scalars ..
  integer :: i,n,nbc,bc_idx
  !
  type(periodic_bc_t) :: tmp_periodic_bc1
  type(periodic_bc_t) :: tmp_periodic_bc2
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: bcname
  !
  !.. Local Allocatable Arrays ..
  type(bc_input_t), allocatable :: tmp_bc_input(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reset_bc_input_for_1d_grid"
  !
continue
  !
  ! Create a default bc_input for one of the 1D test cases
  !
  ! F2003 AUTO-REALLOC
  tmp_bc_input = get_default_1d_bc_input(periodic_distance,inflow,outflow)
  !
  nbc = size(tmp_bc_input)
  !
  ! If bc_input is not allocated, use move_alloc to allocate
  ! bc_input to tmp_bc_input and return from the subroutine.
  !
  if (.not. allocated(bc_input)) then
    call move_alloc( from=tmp_bc_input , to=bc_input )
    return
  end if
  !
  ! Find the inflow BC in bc_input if it exists
  !
  bc_idx = 1
  n = 0
  tmp_periodic_bc1 = tmp_bc_input(bc_idx)%periodic
  find_inflow_loop: &
  do i = 1,size(bc_input)
    !
    bcname = trim(adjustl(uppercase(bc_input(i)%name))) ! F2003 AUTO-REALLOC
    if (bcname == "XMIN") n = i
    if (bcname == "INFLOW") n = i
    if (any(bc_input(i)%set_bc_value == bc_inflows(:))) n = i
    !
    if (n > 0) then
      tmp_bc_input(bc_idx) = bc_input(n)
      tmp_bc_input(bc_idx)%name = "INFLOW"
      exit find_inflow_loop
    end if
    !
  end do find_inflow_loop
  !
  ! Find the outflow bc in bc_input if it exists
  !
  bc_idx = 2
  n = 0
  tmp_periodic_bc2 = tmp_bc_input(bc_idx)%periodic
  find_outflow_loop: &
  do i = 1,size(bc_input)
    !
    bcname = trim(adjustl(uppercase(bc_input(i)%name))) ! F2003 AUTO-REALLOC
    if (bcname == "XMAX") n = i
    if (bcname == "OUTFLOW") n = i
    if (any(bc_input(i)%set_bc_value == bc_outflows(:))) n = i
    !
    if (n > 0) then
      tmp_bc_input(bc_idx) = bc_input(n)
      tmp_bc_input(bc_idx)%name = "OUTFLOW"
      exit find_outflow_loop
    end if
    !
  end do find_outflow_loop
  !
  ! If the inflow and outflow BCs are periodic, copy the default values for
  ! the periodic BCs back into tmp_bc_input
  !
  if (all(tmp_bc_input(1:2)%set_bc_value == bc_periodic)) then
    tmp_bc_input(1)%periodic = tmp_periodic_bc1
    tmp_bc_input(2)%periodic = tmp_periodic_bc2
  end if
  !
  ! Finally, use move_alloc to replace whatever might be in bc_input
  ! with the reconstructed complete bc_input contained in tmp_bc_input.
  !
  call move_alloc( from=tmp_bc_input , to=bc_input )
  !
  ! Reset bc_names to the names just assigned to bc_input for each BC
  !
  bc_names = bc_input(1:nbc)%name ! F2003 AUTO-REALLOC
  !
  contains
  !
  pure function get_default_1d_bc_input(periodic_distance,inflow,outflow) &
                                 result(return_value)
    !
    !.. Optional Arguments ..
    real(wp), optional, intent(in) :: periodic_distance
    integer,  optional, intent(in) :: inflow
    integer,  optional, intent(in) :: outflow
    !
    !.. Function Result ..
    type(bc_input_t), allocatable :: return_value(:)
    !
    !.. Local Scalars ..
    integer :: ierr,inflow_bc,outflow_bc
    logical(lk) :: set_periodic
    !
  continue
    !
    set_periodic = fals
    if (present(periodic_distance)) then
      if (abs(periodic_distance) > eps6) then
        set_periodic = true
      end if
    end if
    !
    inflow_bc = bc_sup_outflow
    if (present(inflow)) then
      inflow_bc = inflow
    end if
    !
    outflow_bc = bc_sup_outflow
    if (present(outflow)) then
      outflow_bc = outflow
    end if
    !
    allocate ( return_value(1:2) , stat=ierr )
    !
    ! Default both boundary conditions to supersonic outflow
    ! NOTE: This has seemed to work best for non-periodic 1D problems
    !
    return_value(1)%name = "INFLOW"
    return_value(1)%set_bc_value = inflow_bc
    !
    return_value(2)%name = "OUTFLOW"
    return_value(2)%set_bc_value = outflow_bc
    !
    if (set_periodic) then
      !
      ! Set the boundary conditions to periodic if periodic_distance
      ! is present and has a magnitude greater than 1.0e-6.
      !
      return_value(1)%set_bc_value = bc_periodic
      return_value(1)%periodic%matching_bc_name = "OUTFLOW"
      return_value(1)%periodic%matching_bc_idx  = 2
      return_value(1)%periodic%translation(1)   = periodic_distance
      !
      return_value(2)%set_bc_value = bc_periodic
      return_value(2)%periodic%matching_bc_name = "INFLOW"
      return_value(2)%periodic%matching_bc_idx  = 1
      return_value(2)%periodic%translation(1)   = -periodic_distance
      !
    end if
    !
  end function get_default_1d_bc_input
  !
end subroutine reset_bc_input_for_1d_grid
!
!###############################################################################
!
subroutine create_postproc_debug_grid
  !
  !.. Use Statements ..
  use order_mod,   only : n_order
  use ovar,        only : Lref
  use geomvar_mod, only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geomvar_mod, only : bface,nodes_of_cell_ptr,nodes_of_cell
  use geomvar_mod, only : cell_geom,cell_order
  use ovar,        only : VortexGrid_random_perturb
  use ovar,        only : VortexGrid_perturb_percent
  use ovar,        only : grid_format
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
function tanh_spacing_1d(x1in,x2in,ds1in,ds2in,npts,special_b_coeff, &
                         fixed_end_spacings) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: x1in
  real(wp), intent(in) :: x2in
  real(wp), intent(in) :: ds1in
  real(wp), intent(in) :: ds2in
  integer,  intent(in) :: npts
  !
  !.. Optional Arguments ..
  real(wp),    optional, intent(in) :: special_b_coeff
  logical(lk), optional, intent(in) :: fixed_end_spacings
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:)
  !
  !.. Local Scalars ..
  integer  :: nl,np,n1,n2,l
  real(wp) :: dl,x,c0
  !
  logical(lk) :: enforce_end_spacings
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
  enforce_end_spacings = true
  if (present(fixed_end_spacings)) then
    enforce_end_spacings = fixed_end_spacings
  end if
  !
  if (enforce_end_spacings) then
    !
    ! We are going to enforce the end spacings ds1in and ds2in,
    ! so we only want to solve for the npts-3 interior spacings
    ! (meaning npts-2 points on the line).
    !
    np = npts-2
    !
    n1 = 2
    n2 = np+1
    !
    c0 = one
    !
  else
    !
    ! We are not enforcing the end spacings, so we need to solve for
    ! all npts-1 interior spacings (meaning npts points on the line).
    !
    np = npts
    !
    n1 = 1
    n2 = np
    !
    c0 = zero ! cancel the adjustment to x1i and x2i
    !
  end if
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
  ! NOTE: If not enforcing end spacings, c0=zero and
  !       this just sets [x1i,x2i] = [x1,x2]
  !
  do l = 1,nl
    x1i(l) = x1(l) + c0*(x2(l)-x1(l))*ds1in/dl
    x2i(l) = x2(l) - c0*(x2(l)-x1(l))*ds2in/dl
  end do
  !
  ! F2003 AUTO-REALLOC
  xp = get_tanh_spacing(x1i,x2i,ds1in,ds2in,np,special_b_coeff)
  !
  ! Copy the interior point spacings to the function result
  !
  return_value(n1:n2) = xp(1,1:np)
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
function tanh_spacing_md(x1in,x2in,ds1in,ds2in,npts,special_b_coeff, &
                         fixed_end_spacings) result(return_value)
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: x1in(:)
  real(wp), intent(in) :: x2in(:)
  real(wp), intent(in) :: ds1in
  real(wp), intent(in) :: ds2in
  integer,  intent(in) :: npts
  !
  !.. Optional Arguments ..
  real(wp),    optional, intent(in) :: special_b_coeff
  logical(lk), optional, intent(in) :: fixed_end_spacings
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: nl,np,n1,n2,l
  real(wp) :: dl,c0
  !
  logical(lk) :: enforce_end_spacings
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
  enforce_end_spacings = true
  if (present(fixed_end_spacings)) then
    enforce_end_spacings = fixed_end_spacings
  end if
  !
  if (enforce_end_spacings) then
    !
    ! We are going to enforce the end spacings ds1in and ds2in,
    ! so we only want to solve for the npts-3 interior spacings
    ! (meaning npts-2 points on the line).
    !
    np = npts-2
    !
    n1 = 2
    n2 = np+1
    !
    c0 = one
    !
  else
    !
    ! We are not enforcing the end spacings, so we need to solve for
    ! all npts-1 interior spacings (meaning npts points on the line).
    !
    np = npts
    !
    n1 = 1
    n2 = np
    !
    c0 = zero ! cancel the adjustment to x1i and x2i
    !
  end if
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
  ! NOTE: If not enforcing end spacings, c0=zero and
  !       this just sets [x1i,x2i] = [x1,x2]
  !
  do l = 1,nl
    x1i(l) = x1(l) + c0*(x2(l)-x1(l))*ds1in/dl
    x2i(l) = x2(l) - c0*(x2(l)-x1(l))*ds2in/dl
  end do
  !
  ! F2003 AUTO-REALLOC
  xp = get_tanh_spacing(x1i,x2i,ds1in,ds2in,np,special_b_coeff)
  !
  ! Copy the interior point spacings to the function result
  !
  return_value(1:nl,n1:n2) = xp(1:nl,1:np)
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
