module mappings_mod
  !
  !.. Global Use Statements ..
  use module_kind_types
  use matrix_types_mod, only : matrix
  use eqn_idx, only : nec,nmb,nme,nee,nq
  !
  implicit none
  !
  private
  !
#ifdef DEBUG_ON
!!!#define TIME_PARAMETRIC_MAPPING
!!!#define OUTPUT_GRID_MAPPINGS
#endif
  !
  !
  ! ##################################
  ! ###  Module Public Procedures  ###
  ! ##################################
  !
  public :: get_solution_points
  public :: mapping
  !
  type :: mapping_matrices_t
    integer :: min_order = 1
    integer :: max_order = 1
    integer :: min_edge_def = huge(0)
    integer :: max_edge_def = -huge(0)
    integer :: min_face_def = huge(0)
    integer :: max_face_def = -huge(0)
    type(matrix), allocatable :: MeshOrder(:)
    type(matrix), allocatable :: EdgeDefOrder(:)
    type(matrix), allocatable :: FaceDefOrder(:)
  end type mapping_matrices_t
  !
  type(mapping_matrices_t), save, target, allocatable :: mapping(:,:)
  !
  !
  !
 !logical(lk), parameter :: edge_defined_uses_endpts = fals
  logical(lk), parameter :: edge_defined_uses_endpts = true
  !
  ! Generic interfaces
  !
  interface map_edge_to_solpts
    module procedure equidistant_map_edge_to_solpts, &
                     parametric_map_edge_to_solpts
  end interface
  !
  interface map_tria_to_solpts
    module procedure equidistant_map_tria_to_solpts
  end interface
  !
  interface map_quad_to_solpts
    module procedure equidistant_map_quad_to_solpts, &
                     parametric_map_quad_to_solpts
  end interface
  !
  interface map_tetr_to_solpts
    module procedure equidistant_map_tetr_to_solpts
  end interface
  !
  interface map_pyra_to_solpts
    module procedure equidistant_map_pyra_to_solpts
  end interface
  !
  interface map_pris_to_solpts
    module procedure equidistant_map_pris_to_solpts
  end interface
  !
  interface map_hexa_to_solpts
    module procedure equidistant_map_hexa_to_solpts, &
                     parametric_map_hexa_to_solpts
  end interface
  !
contains
!
!###############################################################################
!
subroutine get_solution_points
  !
  !.. Use Statements ..
  use ovar, only : grid_format
  use ovar, only : parametric_agglomeration
  use ovar, only : plot3d_agglomerate_order
  use ovar, only : cgns_grid_needs_correction
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_solution_points"
  !
continue
  !
  if (cgns_grid_needs_correction) then
    !
    call custom_equidistant_mapping
    !
  else
    !
    if (grid_format == Plot3D_Format .and. &
        plot3d_agglomerate_order > 0 .and. &
        parametric_agglomeration /= 0) then
      !
      if (parametric_agglomeration == Axisymmetric_Parameterization) then
        write (error_message,1)
        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      end if
      !
      call parametrically_corrected_mapping
      !
    else
      !
      call assume_equidistant_mapping
      !
    end if
    !
  end if
  !
  ! Format Statements
  !
  1 format ("Axisymmetric parameterization to correct mappings for ", &
            "Plot3D agglomerated elements is not yet supported.")
  !
end subroutine get_solution_points
!
!###############################################################################
!
subroutine custom_equidistant_mapping()
  !
  use geovar,    only : grid
  use geovar,    only : nr,nfbnd,ncell,nnode,n_totcel,n_totnod
  use geovar,    only : n_solpts,n_totpts
  use geovar,    only : nodes_of_cell_ptr,nodes_of_cell
  use geovar,    only : bface,xyz_nodes,xyz
  use geovar,    only : cell_geom,cell_order
  !
  use order_mod, only : geom_solpts
  !
  use ovar,      only : rhoref,aref
  use ovar,      only : interpto_input_mode
  !
  use flowvar,   only : usp
  !
  !.. Local Scalars ..
  integer :: n,n0,n1,n2,np,p1,p2,ierr
  integer :: this_cell,mesh_order
  integer :: this_geom,this_order
  logical(lk) :: use_nodes_of_cell
  logical(lk) :: check_grid_for_pts
  logical(lk) :: grid_sol_is_available
  real(wp) :: rhoref_inv,rhoa_inv,rhoa2_inv
  !
  integer :: l,io,io1,io2
  !
  type(geom_family_t) :: elem_family
  !
  !.. Local Pointers ..
  type(matrix), pointer :: this_map
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: pts(:)
  real(wp), allocatable :: elem_xyz(:,:)
  real(wp), allocatable :: elem_sol(:,:)
  real(wp), allocatable :: elem_usp(:,:)
  real(wp), allocatable :: custom_map(:,:)
  real(wp), allocatable :: eqd_xyz(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "custom_equidistant_mapping"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  check_grid_for_pts = fals
  !
  rhoref_inv = one / rhoref
  rhoa_inv = rhoref_inv / aref
  rhoa2_inv = rhoa_inv / aref
  !
  ! We need to reset various grid size parameters in case partitioning
  ! was done for use with multiple processors. If only using with a
  ! single processor, this should just give the values we already have.
  !
  nfbnd = size(bface,dim=2)
  n_totcel = size(nodes_of_cell_ptr) - 1
  ncell = minval(bface(10,:)) - 1
  n_totnod = size(xyz_nodes,dim=2)
  nnode = maxval(nodes_of_cell(1:nodes_of_cell_ptr(ncell+1)))
  !
  ! Allocate the mapping array and the component
  ! arrays MeshOrder, EdgeDefOrder, and FaceDefOrder
  !
  call create_mapping_array(cell_geom,cell_order,grid)
  !
  ! Get the total number of interior solution points
  !
  n_solpts = 0
  do n = 1,ncell
    n_solpts = n_solpts + geom_solpts( cell_geom(n) , cell_order(n) )
  end do
  !
  ! Get the total number of interior plus ghost cell solution points
  !
  n_totpts = n_solpts
  do n = ncell+1,n_totcel
    n_totpts = n_totpts + geom_solpts( cell_geom(n) , cell_order(n) )
  end do
  !
  ! Now calculate the coordinates of the solution points within each cell
  !
  xyz = init_array(nr,n_solpts,zero) ! F2003 AUTO-REALLOC
  !$ACC UPDATE DEVICE(xyz)
  !
  ! Determine whether grid%sol is available
  !
  grid_sol_is_available = fals
  if (allocated(grid)) then
    if (allocated(grid%sol)) then
      grid_sol_is_available = true
      usp = init_array(n_solpts,nq,zero) ! F2003 AUTO-REALLOC
      !$ACC UPDATE DEVICE(usp)
    end if
  end if
  !
  ! Initialize the counting index to current
  ! storage location within the xyz array
  !
  n2 = 0
  !
  ! Loop over all the cells local to this processor and map the coordinates
  ! of the points defining each element that were given by the grid file to
  ! the coordinates of the solution points within each cell
  !
!#define OUTPUT_GRID_MAPPINGS
#ifdef OUTPUT_GRID_MAPPINGS
  open (newunit=io,file="xyz.dat",status="replace",action="write",iostat=ierr)
  open (newunit=io1,file="xyz.x",status="replace",action="write",iostat=ierr)
  open (newunit=io2,file="grid.x",status="replace",action="write",iostat=ierr)
  !
  write (io1,*) ncell
  write (io2,*) ncell
  do this_cell = 1,ncell
    this_order = cell_order(this_cell)
    write (io1,'(3i6)') this_order+1,this_order+1,this_order+1
    write (io2,'(3i6)') this_order+1,this_order+1,this_order+1
  end do
#endif
  !
  do this_cell = 1,ncell
    !
    this_geom = cell_geom(this_cell)
    this_order = cell_order(this_cell)
    !
    np = geom_solpts(this_geom,this_order) ! number of solution points
    !
    n0 = n2       ! offset to location within xyz for this cell
    n1 = n2 + 1   ! index of first point of this cell within xyz
    n2 = n2 + np  ! index of last  point of this cell within xyz
    !
    ! Get the mesh order and family type for this cell 
    !
    mesh_order = grid%elem(this_cell)%prop%order
    elem_family = grid%elem(this_cell)%prop%family
    !
    ! Abort if the family type for this cell is not complete
    if (elem_family /= Complete) then
      if (elem_family == Edge_Defined) then
        write (iout,1) this_cell,"Edge"
      else if (elem_family == Face_Defined) then
        write (iout,1) this_cell,"Face"
      else
        write (iout,1) this_cell,"Unknown"
      end if
      call stop_gfr(abort,pname,__LINE__,__FILE__)
    end if
              
    !
    ! Collect the coordinates for the points of
    ! this cell from grid into a local array
    !
    pts = grid%elem(this_cell)%pts(:) ! F2003 AUTO-REALLOC
    !
    elem_xyz = grid%xyz(1:nr,pts(:)) ! F2003 AUTO-REALLOC
    !
!#ifdef DEBUG_ON
!block
!    real(wp) :: x,y,z,r,t
!    if (mypnum == glb_root) then
!      if (this_cell == 1) then
!        write (iout,'(///)')
!        write (iout,184) "mesh_order = ",mesh_order
!        write (iout,184) "cell_order = ",this_order
!        do n = 1,size(elem_xyz,dim=2)
!          x = elem_xyz(1,n)
!          y = elem_xyz(2,n)
!          z = elem_xyz(3,n)
!          r = sqrt(y*y + z*z)
!          t = atan2(z,y)
!          if (t < zero) t = t + two*pi
!          write (iout,185) n,x,y,z,x,r,t*radians2degrees
!        end do
!        write (iout,'(///)')
!      end if
!    end if
!    184 format (a,i0)
!    185 format (i6,2(3x,3f12.6))
!end block
!#endif
    !
    ! Assign a local pointer to the mapping
    ! matrix as an alias to simplify the code
    !
    if (elem_family == Edge_Defined) then
      this_map => mapping(this_geom,this_order)%EdgeDefOrder(mesh_order)
    else if (elem_family == Face_Defined) then
      this_map => mapping(this_geom,this_order)%FaceDefOrder(mesh_order)
    else
      this_map => mapping(this_geom,this_order)%MeshOrder(mesh_order)
    end if
    !
    ! Create the mapping matrix for this combination of cell geometry and
    ! cell order if the mat component of this_map has not been allocated yet
    ! NOTE: We dont need this now since we are doing a custom mapping, so this
    !       is only done just in case the mapping array is used somewhere
    !       later outside this module. 
    !
    if (.not.allocated(this_map%mat)) then
      !
      ! F2003 AUTO-REALLOC
      this_map%mat = map_general_mesh_elem_to_solpts(this_geom, &
                                                     this_order, &
                                                     mesh_order)
      !
    end if
    !
    ! Create a custom mapping matrix for this cell that (possibly/hopefully)
    ! corrects for any errors with the high-order points from the grid file
    ! not being defined at equally distant points as defined by the CGNS
    ! standard.
    !
    ! F2003 AUTO-REALLOC
    custom_map = map_custom_mesh_elem_to_solpts(this_geom,this_order, &
                                                mesh_order,elem_xyz)
    !
    xyz(1:nr,n1:n2) = matmul( elem_xyz(:,:) , custom_map(:,:) )
    !
!#define USE_CUSTOM_MAPS_TO_FIX_PW
#ifdef USE_CUSTOM_MAPS_TO_FIX_PW
block
  real(wp) :: x,y,z,r,t,xmin,xmax,rmin,rmax
    if (ncpu == 1) then
      xmin = minval(elem_xyz(1,:))
      xmax = maxval(elem_xyz(1,:))
      rmin = minval( sqrt( (elem_xyz(2,:)**2) + (elem_xyz(3,:)**2) ) )
      rmax = maxval( sqrt( (elem_xyz(2,:)**2) + (elem_xyz(3,:)**2) ) )
     !if (this_cell == 26) then
     !if (xmin > -3.45_wp .and. xmax < -3.25_wp .and. &
     !    rmin > 0.88_wp .and. rmax < one) then
        write (36,725) this_cell
        eqd_xyz = matmul( elem_xyz(:,:) , this_map%mat(:,:) ) ! F03 AUTO-REALLOC
        do n = 1,size(eqd_xyz,dim=2)
          x = eqd_xyz(1,n)
          y = eqd_xyz(2,n)
          z = eqd_xyz(3,n)
          r = sqrt(y*y+z*z)
          t = atan2(z,y)
          t = (t + merge(two*pi,zero,t<zero)) * radians2degrees
          write (36,726) n,x,y,z,r,t
        end do
        write (37,725) this_cell
        do n = n1,n2
          x = xyz(1,n)
          y = xyz(2,n)
          z = xyz(3,n)
          r = sqrt(y*y+z*z)
          t = atan2(z,y)
          t = (t + merge(two*pi,zero,t<zero)) * radians2degrees
          write (37,726) n-n1+1,x,y,z,r,t
        end do
     !end if
    end if
    725 format (/,2x,"Cell # ",i0)
    726 format (2x,i4,5(2x,f13.8))
   !if (this_cell == 26) call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
end block
#endif
    !
#ifdef OUTPUT_GRID_MAPPINGS
    write (io,*)
    write (io,203) "ell",this_cell
    do n = 1,size(elem_xyz,dim=2)
      if (nr == 2) then
        write (io,204) (elem_xyz(l,n),l=1,nr)
      else if (nr == 3) then
        write (io,304) (elem_xyz(l,n),l=1,nr)
      end if
    end do
    if (nr == 2) then
      write (io,205) (elem_xyz(l,1),l=1,nr)
    else if (nr == 3) then
      write (io,305) (elem_xyz(l,1),l=1,nr)
    end if
    write (io,206)
    !
    write (io,203) "pts",this_cell
    do n = n1,n2-1
      if (nr == 2) then
        write (io,204) (xyz(l,n),l=1,nr)
      else if (nr == 3) then
        write (io,304) (xyz(l,n),l=1,nr)
      end if
    end do
    if (nr == 2) then
      write (io,205) (xyz(l,n2),l=1,nr)
    else if (nr == 3) then
      write (io,305) (xyz(l,n2),l=1,nr)
    end if
    write (io,206)
    203 format ("c",a3,i0,"={")
    204 format ("{",f7.4,",",f7.4,"},")
    205 format ("{",f7.4,",",f7.4,"}")
    206 format ("};")
    304 format ("{",f7.4,",",f7.4,",",f7.4,"},")
    305 format ("{",f7.4,",",f7.4,",",f7.4,"}")
    !
   !write (io,100) this_cell
   !write (io,*)
   !do n = 1,size(elem_xyz,dim=2)
   !  write (io,101) n,(elem_xyz(l,n),l=1,nr)
   !end do
   !write (io,*)
   !do n = n1,n2
   !  write (io,101) n,(xyz(l,n),l=1,nr)
   !end do
    100 format (/," Cell # ",i0)
    101 format (i8,3x,3f10.4)
    !
    do l = 1,nr
      write (io1,102) (xyz(l,n),n=n1,n2)
      write (io2,102) (elem_xyz(l,n),n=1,size(elem_xyz,dim=2))
    end do
    102 format (5es18.10)
#endif
    !
    if (grid_sol_is_available) then
      !
      elem_sol = grid%sol(1:nq,pts(:)) ! F2003 AUTO-REALLOC
      elem_usp = matmul( elem_sol(:,:), custom_map(:,:) ) ! F2003 AUTO-REALLOC
      !
      ! Copy the transposed solution back into usp
      ! NOTE: elem_usp should be sized elem_usp(1:nq,1:np)
      ! NOTE: grid%sol is expected to be dimensional and in metric units,
      !       and thus elem_sol and elem_usp should also be dimensional and
      !       in metrics units. While we are copying the transpose of elem_usp
      !       into usp, nondimensionalize the flow variables.
      !
      do n = 1,np
        usp(n0+n,  1:nec) = elem_usp(  1:nec,n) * rhoref_inv
        usp(n0+n,nmb:nme) = elem_usp(nmb:nme,n) * rhoa_inv
        usp(n0+n,    nee) = elem_usp(    nee,n) * rhoa2_inv
      end do
     !usp(n1:n2,1:nq) = transpose( elem_usp(:,:) )
      !
    end if
    !
    if (associated(this_map)) this_map => null()
    !
  end do
  !
#ifdef OUTPUT_GRID_MAPPINGS
  close (io)
  close (io1)
  close (io2)
#endif
#ifdef USE_CUSTOM_MAPS_TO_FIX_PW
   call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
#endif
  !
  ! Deallocate the local arrays before leaving
  !
  if (allocated(pts)) then
    deallocate ( pts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"pts",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(eqd_xyz)) then
    deallocate ( eqd_xyz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"eqd_xyz",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(elem_xyz)) then
    deallocate ( elem_xyz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"elem_xyz",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(elem_sol)) then
    deallocate ( elem_sol , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"elem_sol",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(elem_usp)) then
    deallocate ( elem_usp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"elem_usp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(mapping)) then
    deallocate ( mapping , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"mapping",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(custom_map)) then
    deallocate ( custom_map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"custom_map",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  !$ACC UPDATE DEVICE(xyz)
  if (grid_sol_is_available) then
    !$ACC UPDATE DEVICE(usp)
  end if
  !
  if (interpto_input_mode == 2) then
    if (ncpu > 1) then
      call mpi_bcast(xyz,size(xyz,kind=int_mpi),mpi_flttyp,glb_root, &
                     MPI_COMM_WORLD,mpierr)
      if (grid_sol_is_available) then
        call mpi_bcast(usp,size(usp,kind=int_mpi),mpi_flttyp,glb_root, &
                       MPI_COMM_WORLD,mpierr)
      end if
    end if
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/,"ERROR! Cell # ",i0," is defined using an unsupported",/, &
                  8x,"family type: ",a,"-Defined",/)
  !
end subroutine custom_equidistant_mapping
!
!###############################################################################
!
subroutine parametrically_corrected_mapping()
  !
  use geovar,           only : grid
  use geovar,           only : nr,nfbnd,ncell,nnode,n_totcel,n_totnod
  use geovar,           only : n_solpts,n_totpts
  use geovar,           only : nodes_of_cell_ptr,nodes_of_cell
  use geovar,           only : bface,xyz_nodes,xyz
  use geovar,           only : cell_geom,cell_order
  !
  use order_mod,        only : geom_solpts
  use order_mod,        only : n_min_geom,n_max_geom
  use order_mod,        only : n_min_order,n_max_order
  !
  use ovar,             only : rhoref,aref
  !
  use flowvar,          only : usp
  !
  use element_type_mod, only : std_elem
  !
  !.. Local Scalars ..
  integer :: l,n,n0,n1,n2,np,op,p1,p2,ierr,io,io1,io2
  integer :: this_cell,mesh_order
  integer :: this_geom,this_order
  logical(lk) :: grid_sol_is_available
  real(wp) :: rhoref_inv,rhoa_inv,rhoa2_inv
  !
  character(len=200) :: array_name
  !
  !.. Local Allocatable Scalars ..
  type(matrix), allocatable :: this_map
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: elem_xyz(:,:)
  real(wp), allocatable :: elem_sol(:,:)
  real(wp), allocatable :: elem_usp(:,:)
  real(wp), allocatable :: edge_pts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "parametrically_corrected_mapping"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
#ifdef TIME_PARAMETRIC_MAPPING
  call debug_timer(reset_accumulated_time,"parametric initialization",acc=1)
  call debug_timer(reset_accumulated_time,"create comp_coord",acc=2)
  call debug_timer(reset_accumulated_time,"rescale comp_coord",acc=3)
  call debug_timer(reset_accumulated_time,"create mapping matrix",acc=4)
  call debug_timer(reset_accumulated_time,"chopping mapping matrix",acc=5)
#endif
  !
  rhoref_inv = one / rhoref
  rhoa_inv = rhoref_inv / aref
  rhoa2_inv = rhoa_inv / aref
  !
  ! We need to reset various grid size parameters in case partitioning
  ! was done for use with multiple processors. If only using with a
  ! single processor, this should just give the values we already have.
  !
  nfbnd = size(bface,dim=2)
  n_totcel = size(nodes_of_cell_ptr) - 1
  ncell = minval(bface(10,:)) - 1
  n_totnod = size(xyz_nodes,dim=2)
  nnode = maxval(nodes_of_cell(1:nodes_of_cell_ptr(ncell+1)))
  !
  ! Get the total number of interior solution points
  !
  n_solpts = 0
  do n = 1,ncell
    n_solpts = n_solpts + geom_solpts( cell_geom(n) , cell_order(n) )
  end do
  !
  ! Get the total number of interior plus ghost cell solution points
  !
  n_totpts = n_solpts
  do n = ncell+1,n_totcel
    n_totpts = n_totpts + geom_solpts( cell_geom(n) , cell_order(n) )
  end do
  !
  ! Now calculate the coordinates of the solution points within each cell
  !
  allocate ( xyz(1:nr,1:n_solpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(xyz)
  !
  ! Determine whether grid%sol is available
  !
  grid_sol_is_available = fals
  if (allocated(grid)) then
    if (allocated(grid%sol)) then
      grid_sol_is_available = true
      allocate ( usp(1:n_solpts,1:nq) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"usp",1,__LINE__,__FILE__,ierr,error_message)
      !$ACC UPDATE DEVICE(usp)
    end if
  end if
  !
  ! Initialize the counting index to current
  ! storage location within the xyz array
  !
  n2 = 0
  !
  ! Loop over all the cells local to this processor and map the coordinates
  ! of the points defining each element that were given by the grid file to
  ! the coordinates of the solution points within each cell
  !
#ifdef OUTPUT_GRID_MAPPINGS
  open (newunit=io,file="xyz.dat",status="replace",action="write",iostat=ierr)
  open (newunit=io1,file="xyz.x",status="replace",action="write",iostat=ierr)
  open (newunit=io2,file="grid.x",status="replace",action="write",iostat=ierr)
  !
  write (io1,*) ncell
  write (io2,*) ncell
  do this_cell = 1,ncell
    this_order = cell_order(this_cell)
    write (io1,'(3i6)') this_order+1,this_order+1,this_order+1
    write (io2,'(3i6)') this_order+1,this_order+1,this_order+1
  end do
#endif
  !
  do this_cell = 1,ncell
    !
    this_geom = cell_geom(this_cell)
    this_order = cell_order(this_cell)
    !
    np = geom_solpts(this_geom,this_order) ! number of solution points
    !
    n0 = n2
    n1 = n2 + 1
    n2 = n2 + np
    !
    ! Using parametric equations to compute custom mappings for each
    ! cell assumes that high-order mesh information is available in
    ! grid%elem(:)%pts, so always use that information instead of
    ! nodes_of_cell_ptr and nodes_of_cell.
    !
    mesh_order = grid%elem(this_cell)%prop%order
    !
    ! Using F2003 auto-reallocation
    elem_xyz = grid%xyz(1:nr, grid%elem(this_cell)%pts(:) )
    !
    ! Deallocate this_map and then reallocate it to remove the custom
    ! mapping information it contains for the previous grid cell
    !
    if (allocated(this_map)) then
      deallocate ( this_map , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"this_map",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    allocate ( this_map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"this_map",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Create the custom mapping matrix for this cell
    !
    op = size(elem_xyz,dim=2)
    !
    ! Allocate the mat component of this_map
    !
    allocate ( this_map%mat(1:op,1:np) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    edge_pts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F2003 AUTO-REALLOCATION
    !
    ! This parametric mapping should only be used with Plot3D grids,
    ! so we only need to worry about structured grid elements.
    !
    select case (this_geom)
      !
      case (Geom_Edge)
        !
        this_map%mat = map_edge_to_solpts(np,edge_pts,elem_xyz)
        !
      case (Geom_Quad)
        !
        this_map%mat = map_quad_to_solpts(np,edge_pts,elem_xyz,mesh_order)
        !
      case (Geom_Hexa)
        !
        this_map%mat = map_hexa_to_solpts(np,edge_pts,elem_xyz,mesh_order)
        !
      case default
        !
        write (iout,2) this_cell,cell_geom(this_cell)
        call stop_gfr(abort,pname,__LINE__,__FILE__)
        !
    end select
    !
    ! Now use the mapping matrix to map the coordinates of the points for the
    ! element provided by the grid file to the coordinates of the solution
    ! points for the cell
    !
    xyz(1:nr,n1:n2) = matmul( elem_xyz(:,:) , this_map%mat(:,:) )
    !
#ifdef OUTPUT_GRID_MAPPINGS
    write (io,100) this_cell
    do n = n1,n2
      write (io,101) n,(xyz(l,n),l=1,nr)
    end do
    100 format (/," Cell # ",i0)
    101 format (i8,3x,3es23.15)
    !
    do l = 1,nr
      write (io1,102) (xyz(l,n),n=n1,n2)
      write (io2,102) (elem_xyz(l,n),n=1,size(elem_xyz,dim=2))
    end do
    102 format (5es18.10)
#endif
    !
    if (grid_sol_is_available) then
      ! USING F2003 AUTO-REALLOCATION
      elem_sol = grid%sol(1:nq, grid%elem(this_cell)%pts(:) )
      ! USING F2003 AUTO-REALLOCATION
      elem_usp = matmul( elem_sol(:,:), this_map%mat(:,:) )
      !
      ! Copy the transposed solution back into usp
      ! NOTE: elem_usp should be sized elem_usp(1:nq,1:np)
      ! NOTE: grid%sol is expected to be dimensional and in metric units,
      !       and thus elem_sol and elem_usp should also be dimensional and
      !       in metrics units. While we are copying the transpose of elem_usp
      !       into usp, nondimensionalize the flow variables.
      !
      do n = 1,np
        usp(n0+n,  1:nec) = elem_usp(  1:nec,n) * rhoref_inv
        usp(n0+n,nmb:nme) = elem_usp(nmb:nme,n) * rhoa_inv
        usp(n0+n,    nee) = elem_usp(    nee,n) * rhoa2_inv
      end do
     !usp(n1:n2,1:nq) = transpose( elem_usp(:,:) )
      !
    end if
    !
  end do
  !
#ifdef OUTPUT_GRID_MAPPINGS
  close (io)
  close (io1)
  close (io2)
#endif
  !
  ! Deallocate the local arrays before leaving
  !
  if (allocated(elem_xyz)) then
    deallocate ( elem_xyz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"elem_xyz",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(elem_sol)) then
    deallocate ( elem_sol , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"elem_sol",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(elem_usp)) then
    deallocate ( elem_usp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"elem_usp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(mapping)) then
    deallocate ( mapping , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"mapping",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  !$ACC UPDATE DEVICE(xyz)
  if (grid_sol_is_available) then
    !$ACC UPDATE DEVICE(usp)
  end if
  !
#ifdef TIME_PARAMETRIC_MAPPING
  call debug_timer(report_accumulated_time,"parametric initialization",acc=1)
  call debug_timer(report_accumulated_time,"create comp_coord",acc=2)
  call debug_timer(report_accumulated_time,"rescale comp_coord",acc=3)
  call debug_timer(report_accumulated_time,"create mapping matrix",acc=4)
  call debug_timer(report_accumulated_time,"chopping mapping matrix",acc=5)
#endif
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("this_map => mapping(",a,",",i0,")%",a,"(",i0,")%mat")
  2 format (" ERROR: Invalid cell geometry!",/, &
            "        Cell = ",i0,/, &
            "        Geom = ",i0,/)
  3 format ("ERROR: ",a," cells are currently not supported!")
  4 format (" std_elem(",a,",",i0,")%",a,1x,a," contiguous")
  !
end subroutine parametrically_corrected_mapping
!
!###############################################################################
!
subroutine assume_equidistant_mapping()
  !
  use geovar,           only : grid
  use geovar,           only : nr,nfbnd,ncell,nnode,n_totcel,n_totnod
  use geovar,           only : n_solpts,n_totpts
  use geovar,           only : nodes_of_cell_ptr,nodes_of_cell
  use geovar,           only : bface,xyz_nodes,xyz
  use geovar,           only : cell_geom,cell_order
  !
  use order_mod,        only : geom_solpts
  !
  use ovar,             only : rhoref,aref
  use ovar,             only : interpto_input_mode
  !
  use flowvar,          only : usp
  !
  !.. Local Scalars ..
  integer :: n,n0,n1,n2,np,p1,p2,ierr
  integer :: this_cell,mesh_order
  integer :: this_geom,this_order
  logical(lk) :: use_nodes_of_cell
  logical(lk) :: check_grid_for_pts
  logical(lk) :: grid_sol_is_available
  real(wp) :: rhoref_inv,rhoa_inv,rhoa2_inv
  !
  integer :: l,io,io1,io2
  !
  type(geom_family_t) :: elem_family
  !
  !.. Local Pointers ..
  type(matrix), pointer :: this_map
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: elem_xyz(:,:)
  real(wp), allocatable :: elem_sol(:,:)
  real(wp), allocatable :: elem_usp(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "assume_equidistant_mapping"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  check_grid_for_pts = fals
  !
  rhoref_inv = one / rhoref
  rhoa_inv = rhoref_inv / aref
  rhoa2_inv = rhoa_inv / aref
  !
  ! We need to reset various grid size parameters in case partitioning
  ! was done for use with multiple processors. If only using with a
  ! single processor, this should just give the values we already have.
  !
  nfbnd = size(bface,dim=2)
  n_totcel = size(nodes_of_cell_ptr) - 1
  ncell = minval(bface(10,:)) - 1
  n_totnod = size(xyz_nodes,dim=2)
  nnode = maxval(nodes_of_cell(1:nodes_of_cell_ptr(ncell+1)))
  !
  ! Allocate the mapping array and the component
  ! arrays MeshOrder, EdgeDefOrder, and FaceDefOrder
  !
  call create_mapping_array(cell_geom,cell_order,grid)
  !
  if (allocated(grid)) then
    check_grid_for_pts = (allocated(grid%elem) .and. allocated(grid%xyz))
  end if
  !
  ! Get the total number of interior solution points
  !
  n_solpts = 0
  do n = 1,ncell
    n_solpts = n_solpts + geom_solpts( cell_geom(n) , cell_order(n) )
  end do
  !
  ! Get the total number of interior plus ghost cell solution points
  !
  n_totpts = n_solpts
  do n = ncell+1,n_totcel
    n_totpts = n_totpts + geom_solpts( cell_geom(n) , cell_order(n) )
  end do
  !
  ! Now calculate the coordinates of the solution points within each cell
  !
  xyz = init_array(nr,n_solpts,zero) ! F2003 AUTO-REALLOC
  !$ACC UPDATE DEVICE(xyz)
  !
  ! Determine whether grid%sol is available
  !
  grid_sol_is_available = fals
  if (allocated(grid)) then
    if (allocated(grid%sol)) then
      grid_sol_is_available = true
      usp = init_array(n_solpts,nq,zero) ! F2003 AUTO-REALLOC
      !$ACC UPDATE DEVICE(usp)
    end if
  end if
  !
  ! Initialize the counting index to current
  ! storage location within the xyz array
  !
  n2 = 0
  !
  ! Loop over all the cells local to this processor and map the coordinates
  ! of the points defining each element that were given by the grid file to
  ! the coordinates of the solution points within each cell
  !
!#define OUTPUT_GRID_MAPPINGS
#ifdef OUTPUT_GRID_MAPPINGS
  open (newunit=io,file="xyz.dat",status="replace",action="write",iostat=ierr)
  open (newunit=io1,file="xyz.x",status="replace",action="write",iostat=ierr)
  open (newunit=io2,file="grid.x",status="replace",action="write",iostat=ierr)
  !
  write (io1,*) ncell
  write (io2,*) ncell
  do this_cell = 1,ncell
    this_order = cell_order(this_cell)
    write (io1,'(3i6)') this_order+1,this_order+1,this_order+1
    write (io2,'(3i6)') this_order+1,this_order+1,this_order+1
  end do
#endif
  !
  do this_cell = 1,ncell
    !
    this_geom = cell_geom(this_cell)
    this_order = cell_order(this_cell)
    !
    np = geom_solpts(this_geom,this_order) ! number of solution points
 !  write (iout,*) "this_cell = ",this_cell
 !  write (iout,*) "this_geom = ",this_geom
 !  write (iout,*) "this_order = ",this_order
 !  write (iout,*) "np = ",np
    !
    n0 = n2
    n1 = n2 + 1
    n2 = n2 + np
    !
    ! Reset use_nodes_of_cell to true for this cell
    !
    use_nodes_of_cell = true
    !
    ! Check to see if we will be able to extract mesh points for this
    ! cell using grid%elem(:)%pts.
    ! If we can, set use_nodes_of_cell to false so that we use the
    ! possibly higher-order mesh information in grid%elem and grid%xyz.
    !
    if (check_grid_for_pts) then
      if (allocated(grid%elem(this_cell)%pts)) then
        use_nodes_of_cell = fals
      end if
    end if
    !
    if (use_nodes_of_cell) then
      !
      mesh_order = 1
      elem_family = Complete
      !
      p1 = nodes_of_cell_ptr(this_cell)+1
      p2 = nodes_of_cell_ptr(this_cell+1)
      !
      ! Using F2003 auto-reallocation
      elem_xyz = xyz_nodes(1:nr, nodes_of_cell(p1:p2) )
      !
    else
      !
      mesh_order = grid%elem(this_cell)%prop%order
      elem_family = grid%elem(this_cell)%prop%family
      !
      ! Using F2003 auto-reallocation
      elem_xyz = grid%xyz(1:nr, grid%elem(this_cell)%pts(:) )
      !
    end if
!#ifdef DEBUG_ON
!block
!    real(wp) :: x,y,z,r,t
!    if (mypnum == glb_root) then
!      if (this_cell == 1) then
!        write (iout,'(///)')
!        write (iout,184) "mesh_order = ",mesh_order
!        write (iout,184) "cell_order = ",this_order
!        do n = 1,size(elem_xyz,dim=2)
!          x = elem_xyz(1,n)
!          y = elem_xyz(2,n)
!          z = elem_xyz(3,n)
!          r = sqrt(y*y + z*z)
!          t = atan2(z,y)
!          if (t < zero) t = t + two*pi
!          write (iout,185) n,x,y,z,x,r,t*radians2degrees
!        end do
!        write (iout,'(///)')
!      end if
!    end if
!    184 format (a,i0)
!    185 format (i6,2(3x,3f12.6))
!end block
!#endif
    !
    !
    ! Assign a local pointer to the mapping
    ! matrix as an alias to simplify the code
    !
    if (elem_family == Edge_Defined) then
      write (iout,1) this_cell,"Edge"
      call stop_gfr(abort,pname,__LINE__,__FILE__)
      this_map => mapping(this_geom,this_order)%EdgeDefOrder(mesh_order)
    else if (elem_family == Face_Defined) then
      write (iout,1) this_cell,"Face"
      call stop_gfr(abort,pname,__LINE__,__FILE__)
      this_map => mapping(this_geom,this_order)%FaceDefOrder(mesh_order)
    else
      this_map => mapping(this_geom,this_order)%MeshOrder(mesh_order)
    end if
    !
    ! Create the mapping matrix for this combination of cell geometry and
    ! cell order if the mat component of this_map has not been allocated yet
    !
    if (.not.allocated(this_map%mat)) then
      !
      ! F2003 AUTO-REALLOC
      this_map%mat = map_general_mesh_elem_to_solpts(this_geom, &
                                                     this_order, &
                                                     mesh_order)
      !
    end if
    !
    ! Now use the mapping matrix to map the coordinates of the points for the
    ! element provided by the grid file to the coordinates of the solution
    ! points for the cell
    !
    xyz(1:nr,n1:n2) = matmul( elem_xyz(:,:) , this_map%mat(:,:) )
!#ifdef DEBUG_ON
!block
!    real(wp) :: x,y,z,r,t
!    if (mypnum == glb_root) then
!      if (this_cell == 1) then
!        write (iout,'(///)')
!        write (iout,184) "mesh_order = ",mesh_order
!        write (iout,184) "cell_order = ",this_order
!        do n = n1,n2
!          x = xyz(1,n)
!          y = xyz(2,n)
!          z = xyz(3,n)
!          r = sqrt(y*y + z*z)
!          t = atan2(z,y)
!          if (t < zero) t = t + two*pi
!          write (iout,185) n,x,y,z,x,r,t*radians2degrees
!        end do
!        write (iout,'(///)')
!      end if
!    end if
!   !call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
!    184 format (a,i0)
!    185 format (i6,2(3x,3f12.6))
!end block
!#endif
    !
#ifdef OUTPUT_GRID_MAPPINGS
    write (io,*)
    write (io,203) "ell",this_cell
    do n = 1,size(elem_xyz,dim=2)
      if (nr == 2) then
        write (io,204) (elem_xyz(l,n),l=1,nr)
      else if (nr == 3) then
        write (io,304) (elem_xyz(l,n),l=1,nr)
      end if
    end do
    if (nr == 2) then
      write (io,205) (elem_xyz(l,1),l=1,nr)
    else if (nr == 3) then
      write (io,305) (elem_xyz(l,1),l=1,nr)
    end if
    write (io,206)
    !
    write (io,203) "pts",this_cell
    do n = n1,n2-1
      if (nr == 2) then
        write (io,204) (xyz(l,n),l=1,nr)
      else if (nr == 3) then
        write (io,304) (xyz(l,n),l=1,nr)
      end if
    end do
    if (nr == 2) then
      write (io,205) (xyz(l,n2),l=1,nr)
    else if (nr == 3) then
      write (io,305) (xyz(l,n2),l=1,nr)
    end if
    write (io,206)
    203 format ("c",a3,i0,"={")
    204 format ("{",f7.4,",",f7.4,"},")
    205 format ("{",f7.4,",",f7.4,"}")
    206 format ("};")
    304 format ("{",f7.4,",",f7.4,",",f7.4,"},")
    305 format ("{",f7.4,",",f7.4,",",f7.4,"}")
    !
   !write (io,100) this_cell
   !write (io,*)
   !do n = 1,size(elem_xyz,dim=2)
   !  write (io,101) n,(elem_xyz(l,n),l=1,nr)
   !end do
   !write (io,*)
   !do n = n1,n2
   !  write (io,101) n,(xyz(l,n),l=1,nr)
   !end do
    100 format (/," Cell # ",i0)
    101 format (i8,3x,3f10.4)
    !
    do l = 1,nr
      write (io1,102) (xyz(l,n),n=n1,n2)
      write (io2,102) (elem_xyz(l,n),n=1,size(elem_xyz,dim=2))
    end do
    102 format (5es18.10)
#endif
    !
    if (grid_sol_is_available) then
      ! USING F2003 AUTO-REALLOCATION
      elem_sol = grid%sol(1:nq, grid%elem(this_cell)%pts(:) )
      ! USING F2003 AUTO-REALLOCATION
      elem_usp = matmul( elem_sol(:,:), this_map%mat(:,:) )
      !
      ! Copy the transposed solution back into usp
      ! NOTE: elem_usp should be sized elem_usp(1:nq,1:np)
      ! NOTE: grid%sol is expected to be dimensional and in metric units,
      !       and thus elem_sol and elem_usp should also be dimensional and
      !       in metrics units. While we are copying the transpose of elem_usp
      !       into usp, nondimensionalize the flow variables.
      !
      do n = 1,np
        usp(n0+n,  1:nec) = elem_usp(  1:nec,n) * rhoref_inv
        usp(n0+n,nmb:nme) = elem_usp(nmb:nme,n) * rhoa_inv
        usp(n0+n,    nee) = elem_usp(    nee,n) * rhoa2_inv
      end do
     !usp(n1:n2,1:nq) = transpose( elem_usp(:,:) )
      !
    end if
    !
    if (associated(this_map)) this_map => null()
    !
  end do
  !
#ifdef OUTPUT_GRID_MAPPINGS
  close (io)
  close (io1)
  close (io2)
#endif
  !
  !
  !
#ifdef DEBUG_ON
!#define DEBUG_MRF_FREESTREAM_TCF
#ifdef DEBUG_MRF_FREESTREAM_TCF
block
  use ovar, only : itestcase
  if (abs(itestcase) == Freestream_TCF) then
    call exact_TCF_solution_points()
  end if
end block
  !
!block
!    real(wp) :: x,y,z,r,t
!    if (mypnum == glb_root) then
!      n2 = 0
!      do this_cell = 1,ncell
!        this_geom = cell_geom(this_cell)
!        this_order = cell_order(this_cell)
!        np = geom_solpts(this_geom,this_order) ! number of solution points
!        n0 = n2
!        n1 = n2 + 1
!        n2 = n2 + np
!        write (701,284) "this_cell = ",this_cell
!        do n = n1,n2
!          x = xyz(1,n)
!          y = xyz(2,n)
!          z = xyz(3,n)
!          r = sqrt(y*y + z*z)
!          t = atan2(z,y)
!          if (t < zero) t = t + two*pi
!          write (701,285) n,x,y,z,x,r,t*radians2degrees
!        end do
!      end do
!    end if
!    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
!    284 format (a,i0)
!    285 format (i6,3x,3f16.12,3x,2f16.12,f18.12)
!end block
#endif
#endif
  !
  ! Deallocate the local arrays before leaving
  !
  if (allocated(elem_xyz)) then
    deallocate ( elem_xyz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"elem_xyz",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(elem_sol)) then
    deallocate ( elem_sol , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"elem_sol",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(elem_usp)) then
    deallocate ( elem_usp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"elem_usp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(mapping)) then
    deallocate ( mapping , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"mapping",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  !$ACC UPDATE DEVICE(xyz)
  if (grid_sol_is_available) then
    !$ACC UPDATE DEVICE(usp)
  end if
  !
  if (interpto_input_mode == 2) then
    if (ncpu > 1) then
      call mpi_bcast(xyz,size(xyz,kind=int_mpi),mpi_flttyp,glb_root, &
                     MPI_COMM_WORLD,mpierr)
      if (grid_sol_is_available) then
        call mpi_bcast(usp,size(usp,kind=int_mpi),mpi_flttyp,glb_root, &
                       MPI_COMM_WORLD,mpierr)
      end if
    end if
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/,"ERROR! Cell # ",i0," is defined using an unsupported",/, &
                  8x,"family type: ",a,"-Defined",/)
  !
end subroutine assume_equidistant_mapping
!
!###############################################################################
!
subroutine exact_TCF_solution_points()
  !
  !.. Use Statements ..
  use geovar,           only : nr,ncell
  use geovar,           only : nodes_of_cell_ptr,nodes_of_cell
  use geovar,           only : xyz_nodes,xyz
  use geovar,           only : cell_order
  !
  use element_type_mod, only : std_elem
  !
  !.. Local Scalars ..
  integer  :: this_cell,this_order
  integer  :: i,j,k,n,np,p1,p2
  real(wp) :: y,z,r,t
  real(wp) :: x_min,x_max,r_min,r_max,t_min,t_max
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: elem_xrt(:,:)
  real(wp), allocatable :: edge_pts(:)
  real(wp), allocatable :: xpts(:)
  real(wp), allocatable :: rpts(:)
  real(wp), allocatable :: tpts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "exact_TCF_solution_points"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize the counting index to current
  ! storage location within the xyz array
  !
  n = 0
  !
  ! Loop over all the cells local to this processor and map the coordinates
  ! of the points defining each element that were given by the grid file to
  ! the coordinates of the solution points within each cell
  !
  do this_cell = 1,ncell
    !
    this_order = cell_order(this_cell)
    !
    p1 = nodes_of_cell_ptr(this_cell)+1
    p2 = nodes_of_cell_ptr(this_cell+1)
    !
    elem_xrt = xyz_nodes(1:nr, nodes_of_cell(p1:p2) ) ! F2003 AUTO-REALLOC
    !
    ! Convert elem_xrt from Cartesian to cylindrical coordinates
    !
    do k = 1,size(elem_xrt,dim=2)
      y = elem_xrt(2,k)
      z = elem_xrt(3,k)
      r = sqrt(y*y + z*z)
      t = atan2(z,y)
      elem_xrt(2,k) = r
      elem_xrt(3,k) = t + merge(two*pi,zero,t<zero)
    end do
    !
    x_min = elem_xrt(1,1)
    x_max = elem_xrt(1,2)
    !
    r_min = elem_xrt(2,1)
    r_max = elem_xrt(2,4)
    !
    t_min = elem_xrt(3,1)
    t_max = elem_xrt(3,5)
    if (t_max < t_min) t_max = two*pi
    !
    edge_pts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F2003 AUTO-REALLOC
    !
    edge_pts(:) = half*(edge_pts(:) + one)
    !
    xpts = x_min + edge_pts(:)*(x_max-x_min) ! F2003 AUTO-REALLOC
    rpts = r_min + edge_pts(:)*(r_max-r_min) ! F2003 AUTO-REALLOC
    tpts = t_min + edge_pts(:)*(t_max-t_min) ! F2003 AUTO-REALLOC
    !
    do k = 1,size(tpts)
      do j = 1,size(rpts)
        do i = 1,size(xpts)
          !
          n = n + 1
          !
          xyz(1,n) = xpts(i)
          xyz(2,n) = rpts(j)*cos(tpts(k))
          xyz(3,n) = rpts(j)*sin(tpts(k))
          !
        end do
      end do
    end do
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine exact_TCF_solution_points
!
!###############################################################################
!
subroutine create_mapping_array(cell_geom,cell_order,grid)
  !
  !.. Use Statements ..
  use geovar, only : grid_t
  !
  !.. Formal Arguments ..
  integer,                   intent(in) :: cell_geom(:)
  integer,                   intent(in) :: cell_order(:)
  type(grid_t), allocatable, intent(in) :: grid
  !
  !.. Local Scalars ..
  integer :: gmin,gmax,omin,omax,ierr
  integer :: this_cell,elem_order
  integer :: this_geom,this_order
  !
  character(len=200) :: array_name
  !
  !.. Local Pointers ..
  type(mapping_matrices_t), pointer :: this_map
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_mapping_array"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize the local pointers
  !
  this_map => null()
  !
  ! Initialize the minimum and maximum cell geometries and orders
  !
  gmin = minval(cell_geom)
  gmax = maxval(cell_geom)
  omin = minval(cell_order)
  omax = maxval(cell_order)
  !
  ! Allocate the mapping array
  !
  allocate ( mapping(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"mapping",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Find the minimum and maximum order of the complete and
  ! incomplete (i.e. edge-defined and face-defined families) elements
  ! from the grid file
  !
  if (allocated(grid)) then
    !
    if (allocated(grid%elem)) then
      !
      do this_cell = 1,size(grid%elem)
        !
        this_geom = cell_geom(this_cell)
        this_order = cell_order(this_cell)
        !
        elem_order = grid%elem(this_cell)%prop%order
        !
        ! Assign a local pointer to the mapping
        ! matrix as an alias to simplify the code
        !
        this_map => mapping(this_geom,this_order)
        !
        if (grid%elem(this_cell)%prop%family == Edge_Defined) then
          this_map%min_edge_def = min( this_map%min_edge_def, elem_order )
          this_map%max_edge_def = max( this_map%max_edge_def, elem_order )
        else if (grid%elem(this_cell)%prop%family == Face_Defined) then
          this_map%min_face_def = min( this_map%min_face_def, elem_order )
          this_map%max_face_def = max( this_map%max_face_def, elem_order )
        else
          this_map%min_order = min( this_map%min_order, elem_order )
          this_map%max_order = max( this_map%max_order, elem_order )
        end if
        !
        ! Remove the association of the local pointer before continuing
        !
        if (associated(this_map)) this_map => null()
        !
      end do
      !
    end if
    !
  end if
  !
  ! Allocate the MeshOrder, EdgeDefOrder, and FaceDefOrder components
  ! of the mapping array as needed
  !
  do this_geom = lbound(mapping,dim=1),ubound(mapping,dim=1)
    do this_order = lbound(mapping,dim=2),ubound(mapping,dim=2)
      !
      ! Assign a local pointer to the mapping
      ! matrix as an alias to simplify the code
      !
      this_map => mapping(this_geom,this_order)
      !
      omin = this_map%min_order
      omax = this_map%max_order
      !
      if (all([omin,omax] > 0)) then
        allocate ( this_map%MeshOrder(omin:omax) , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"MeshOrder"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message)
      end if
      !
      omin = this_map%min_edge_def
      omax = this_map%max_edge_def
      !
      if (all([omin,omax] > 0)) then
        allocate ( this_map%EdgeDefOrder(omin:omax) , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"EdgeDefOrder"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message)
      end if
      !
      omin = this_map%min_face_def
      omax = this_map%max_face_def
      !
      if (all([omin,omax] > 0)) then
        allocate ( this_map%FaceDefOrder(omin:omax) , &
                   stat=ierr , errmsg=error_message )
        write (array_name,1) Geom_Name(this_geom),this_order,"FaceDefOrder"
        call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                         error_message)
      end if
      !
      ! Remove the association of the local pointer before continuing
      !
      if (associated(this_map)) this_map => null()
      !
    end do
  end do
  !
  ! Just to be sure, make sure to remove any association of
  ! the local pointer before leaving
  !
  if (associated(this_map)) this_map => null()
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("mapping(",a,",",i0,")%",a)
  !
end subroutine create_mapping_array
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function map_general_mesh_elem_to_solpts(this_geom,this_order,mesh_order) &
                                  result(return_value)
  !
  !.. Use Statements ..
  use element_type_mod, only : std_elem
  use quadrature_mod,   only : get_general_quadrature_pts
  use vandermonde_mod,  only : General_Basis_Matrix
  use vandermonde_mod,  only : Vandermonde_Inverse
  !
  !.. Formal Arguments ..
  integer, intent(in) :: this_geom
  integer, intent(in) :: this_order
  integer, intent(in) :: mesh_order
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: grdpts(:,:) ! grid/mesh points
  real(wp), allocatable :: solpts(:,:) ! solution points
  real(wp), allocatable :: vmat_m2n(:,:)
  real(wp), allocatable :: vmat_n2m(:,:)
  real(wp), allocatable :: bmat(:,:)
  real(wp), allocatable :: mmat(:,:)
  !
continue
  !
  ! First get the reference coordinates for the equal distant points
  ! for the grid/mesh element of order mesh_order
  !
  ! F2003 AUTO-REALLOC
  grdpts = get_general_quadrature_pts(this_geom,mesh_order,Equi_Distant_Points)
  !
  ! First create a modal-to-nodal Vandermonde matrix
  ! for the equal distant mesh points
  !
  ! F2003 AUTO-REALLOC
  vmat_m2n = General_Basis_Matrix(this_geom,mesh_order,grdpts)
  !
  ! Create the nodal-to-modal Vandermonde matrix
  ! for the equal distant mesh points
  !
  vmat_n2m = Vandermonde_Inverse( vmat_m2n ) ! F2003 AUTO-REALLOC
  !
  ! Extract the solution points for this geom/order combo
  !
  solpts = std_elem(this_geom,this_order)%pts ! F2003 AUTO-REALLOC
  !
  ! Create the basis matrix of the mesh-order basis functions
  ! evaluated at the solution points
  !
  bmat = General_Basis_Matrix(this_geom,mesh_order,solpts) ! F2003 AUTO-REALLOC
  !
  ! Compute the mapping (interpolation) matrix
  !
  mmat = matmul(transpose(vmat_n2m),transpose(bmat)) ! F2003 AUTO-REALLOC
  !
  ! Use the chop function to copy mmat into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( real(mmat,kind=qp) ) ! F2003 AUTO-REALLOC
  !
end function map_general_mesh_elem_to_solpts
!
!###############################################################################
!
#ifndef DISABLE_PURE
!pure &
#endif
function map_custom_mesh_elem_to_solpts(this_geom,this_order,mesh_order, &
                                        elem_xyz) result(return_value)
  !
  !.. Use Statements ..
  use element_type_mod, only : std_elem
  use vandermonde_mod,  only : General_Basis_Matrix
  use vandermonde_mod,  only : Vandermonde_Inverse
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: this_geom
  integer,  intent(in) :: this_order
  integer,  intent(in) :: mesh_order
  real(wp), intent(in) :: elem_xyz(:,:)
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: grdpts(:,:) ! grid/mesh points
  real(wp), allocatable :: solpts(:,:) ! solution points
  real(wp), allocatable :: vmat_m2n(:,:)
  real(wp), allocatable :: vmat_n2m(:,:)
  real(wp), allocatable :: bmat(:,:)
  real(wp), allocatable :: mmat(:,:)
  !
continue
  !
  ! First get the reference coordinates for the equal distant points
  ! for the grid/mesh element of order mesh_order
  !
  ! F2003 AUTO-REALLOC
  grdpts = get_custom_quadrature_pts(this_geom,mesh_order,elem_xyz)
#ifdef USE_CUSTOM_MAPS_TO_FIX_PW
block
  integer :: i,j,k,n,np
  np = mesh_order+1
  n = 0
  do k = 1,np
    do j = 1,np
      do i = 1,np
        n = n + 1
        write (38,'(2x,i4,3(3x,es24.15))') n,grdpts(1,n),grdpts(2,n),grdpts(3,n)
      end do
    end do
  end do
  write (38,*)
end block
#endif
  !
  ! First create a modal-to-nodal Vandermonde matrix
  ! for the equal distant mesh points
  !
  ! F2003 AUTO-REALLOC
  vmat_m2n = General_Basis_Matrix(this_geom,mesh_order,grdpts)
  !
  ! Create the nodal-to-modal Vandermonde matrix
  ! for the equal distant mesh points
  !
  vmat_n2m = Vandermonde_Inverse( vmat_m2n ) ! F2003 AUTO-REALLOC
  !
  ! Extract the solution points for this geom/order combo
  !
  solpts = std_elem(this_geom,this_order)%pts ! F2003 AUTO-REALLOC
  !
  ! Create the basis matrix of the mesh-order basis functions
  ! evaluated at the solution points
  !
  bmat = General_Basis_Matrix(this_geom,mesh_order,solpts) ! F2003 AUTO-REALLOC
  !
  ! Compute the mapping (interpolation) matrix
  !
  mmat = matmul(transpose(vmat_n2m),transpose(bmat)) ! F2003 AUTO-REALLOC
  !
  ! Use the chop function to copy mmat into return_value while also
  ! setting any values that are less than working machine epsilon to zero.
  !
  return_value = chop( real(mmat,kind=qp) ) ! F2003 AUTO-REALLOC
  !
end function map_custom_mesh_elem_to_solpts
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_custom_quadrature_pts(this_geom,p,xyz) result(return_value)
  !
  !.. Use Statements ..
  use quadrature_mod, only : get_general_quadrature_pts
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: this_geom
  integer,  intent(in) :: p
  real(wp), intent(in) :: xyz(:,:)
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,n,np,i1,i2,j1,j2,k1,k2
  real(wp) :: dl
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: ds(:)
  real(wp), allocatable :: x(:,:,:)
  real(wp), allocatable :: y(:,:,:)
  real(wp), allocatable :: z(:,:,:)
  !
continue
  !
  if (this_geom /= Geom_Hexa) then
    ! F2003 AUTO-REALLOC
    return_value = get_general_quadrature_pts(this_geom,p,Equi_Distant_Points)
    return
  end if
  !
  np = p+1
  !
  x = init_array(np,np,np,zero) ! F2003 AUTO-REALLOC
  y = init_array(np,np,np,zero) ! F2003 AUTO-REALLOC
  z = init_array(np,np,np,zero) ! F2003 AUTO-REALLOC
  ds = init_array(p,zero) ! F2003 AUTO-REALLOC
  !
  return_value = init_array(3,np**3,zero) ! F2003 AUTO-REALLOC
  !
  ! Create the x/i/xi coordinates
  !
  do k = 1,np
    do j = 1,np
      !
      dl = zero
      do i = 1,p
        i1 = map3dto1d(i  ,j,k,np)
        i2 = map3dto1d(i+1,j,k,np)
        ds(i) = norm2( xyz(:,i2) - xyz(:,i1) )
        dl = dl + ds(i)
      end do
      !
      ! Normalize ds from the interval [0,dl] to [0,2]
      !
      ds(:) = two*ds(:)/dl
      !
      x(1,j,k) = -one
      do i = 1,p-1
        x(i+1,j,k) = x(i,j,k) + ds(i)
      end do
      x(np,j,k) = one
      !
    end do
  end do
  !
  ! Create the y/j/eta coordinates
  !
  do i = 1,np
    do k = 1,np
      !
      dl = zero
      do j = 1,p
        j1 = map3dto1d(i,j  ,k,np)
        j2 = map3dto1d(i,j+1,k,np)
        ds(j) = norm2( xyz(:,j2) - xyz(:,j1) )
        dl = dl + ds(j)
      end do
      !
      ! Normalize ds from the interval [0,dl] to [0,2]
      !
      ds(:) = two*ds(:)/dl
      !
      y(i,1,k) = -one
      do j = 1,p-1
        y(i,j+1,k) = y(i,j,k) + ds(j)
      end do
      y(i,np,k) = one
      !
    end do
  end do
  !
  ! Create the z/k/zeta coordinates
  !
  do j = 1,np
    do i = 1,np
      !
      dl = zero
      do k = 1,p
        k1 = map3dto1d(i,j,k  ,np)
        k2 = map3dto1d(i,j,k+1,np)
        ds(k) = norm2( xyz(:,k2) - xyz(:,k1) )
        dl = dl + ds(k)
      end do
      !
      ! Normalize ds from the interval [0,dl] to [0,2]
      !
      ds(:) = two*ds(:)/dl
      !
      z(i,j,1) = -one
      do k = 1,p-1
        z(i,j,k+1) = z(i,j,k) + ds(k)
      end do
      z(i,j,np) = one
      !
    end do
  end do
  !
  ! Copy x/y/z into return_value
  !
  do k = 1,np
    do j = 1,np
      do i = 1,np
        n = map3dto1d(i,j,k,np)
        return_value(1,n) = x(i,j,k)
        return_value(2,n) = y(i,j,k)
        return_value(3,n) = z(i,j,k)
      end do
    end do
  end do
  !
  contains
  !
  pure function map3dto1d(i,j,k,n) result(return_value)
    integer, intent(in) :: i,j,k,n
    integer :: return_value
    return_value = i + n*(j-1) + n*n*(k-1)
  end function map3dto1d
  !
end function get_custom_quadrature_pts
!
!###############################################################################
!
function old_mappings(this_geom,this_order,mesh_order,op,np,family,nc) &
               result(return_value)
  !
  !.. Use Statements ..
  use element_type_mod, only : std_elem
  !
  !.. Formal Arguments ..
  integer,             intent(in) :: this_geom
  integer,             intent(in) :: this_order
  integer,             intent(in) :: mesh_order
  integer,             intent(in) :: op
  integer,             intent(in) :: np
  type(geom_family_t), intent(in) :: family
  integer,             intent(in) :: nc
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: edge_pts(:)
  real(wp), allocatable :: tria_pts(:,:)
  real(wp), allocatable :: tetr_pts(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "old_mappings"
  !
continue
  !
  select case (this_geom)
    !
    case (Geom_Edge)
      !
      edge_pts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F2003 AUTO-REALLOC
      !
      ! F2003 AUTO-REALLOC
      return_value = map_edge_to_solpts( op , np , edge_pts , mesh_order )
      !
    case (Geom_Tria)
      !
      tria_pts = std_elem(Geom_Tria,this_order)%pts  ! F2003 AUTO-REALLOC
      !
      ! F2003 AUTO-REALLOC
      return_value = map_tria_to_solpts( tria_pts , mesh_order )
      !
    case (Geom_Quad)
      !
      edge_pts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F2003 AUTO-REALLOC
      !
      ! F2003 AUTO-REALLOC
      return_value = map_quad_to_solpts( op , np , edge_pts , &
                                         mesh_order , family )
      !
    case (Geom_Tetr)
      !
      tetr_pts = std_elem(Geom_Tetr,this_order)%pts ! F2003 AUTO-REALLOC
      !
      ! F2003 AUTO-REALLOC
      return_value = map_tetr_to_solpts( tetr_pts , mesh_order )
      !
    case (Geom_Pyra)
      !
      edge_pts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F2003 AUTO-REALLOC
      tria_pts = std_elem(Geom_Tria,this_order)%pts      ! F2003 AUTO-REALLOC
      !
      ! F2003 AUTO-REALLOC
      return_value = map_pyra_to_solpts( edge_pts , tria_pts , mesh_order )
      !
    case (Geom_Pris)
      !
      edge_pts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F2003 AUTO-REALLOC
      tria_pts = std_elem(Geom_Tria,this_order)%pts      ! F2003 AUTO-REALLOC
      !
      ! F2003 AUTO-REALLOC
      return_value = map_pris_to_solpts( edge_pts , tria_pts , mesh_order )
      !
    case (Geom_Hexa)
      !
      edge_pts = std_elem(Geom_Edge,this_order)%pts(1,:) ! F2003 AUTO-REALLOC
      !
      ! F2003 AUTO-REALLOC
      return_value = map_hexa_to_solpts( op , np , edge_pts , &
                                         mesh_order , family )
      !
    case default
      !
      write (iout,1) nc,this_geom
      call stop_gfr(abort,pname,__LINE__,__FILE__)
      !
  end select
  !
  ! Format Statements
  !
  1 format (" ERROR: Invalid cell geometry!",/, &
            "        Cell = ",i0,/, &
            "        Geom = ",i0,/)
  !
end function old_mappings
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function equidistant_map_edge_to_solpts(op,np,xi,mesh_order) &
                                 result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  use polynomial_mod, only : equal_distant_quadrature
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: op
  integer,  intent(in) :: np
  real(wp), intent(in) :: xi(:)
  integer,  intent(in) :: mesh_order
  !
  !.. Function Result ..
  real(wp) :: return_value(1:op,1:np)
  !
  !.. Local Scalars ..
  integer :: i,mi,ni,xpts,epts
  !
  !.. Local Arrays ..
  real(mp) :: xi_qp(1:size(xi))
  real(mp) :: eq_pts(1:mesh_order+1)
  real(mp) :: map_qp(1:op,1:np)
  !
continue
  !
  xpts = size(xi_qp)
  epts = size(eq_pts)
  !
  return_value(:,:) = zero
  !
  call equal_distant_quadrature(eq_pts)
  !
  xi_qp(:) = real(xi(:),kind=kind(xi_qp))
  !
  do mi = 1,xpts
    do ni = 1,epts
      map_qp(ni,mi) = eval_LagrangePoly([ni],[xi_qp(mi)],eq_pts)
    end do
  end do
  !
  ! Use the chop function to copy map_qp into return_value while also
  ! setting any values that are less that working machine epsilon to zero.
  !
  return_value = chop( map_qp )
  !
end function equidistant_map_edge_to_solpts
!
!###############################################################################
!
#ifndef DISABLE_PURE
!pure &
#endif
function equidistant_map_tria_to_solpts(rs,mesh_order) result(return_value)
  !
  ! RS : location of solution points on the reference triangle element
  ! MESH_ORDER : polynomial order of the grid file
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: rs(:,:)
  integer,  intent(in) :: mesh_order
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,ni,nj,nk
  integer :: pm,ps,mp,np
  !
continue
  !
  pm = mesh_order                 ! order of grid cells
  mp = geom_solpts(Geom_Tria,pm)  ! no. of pts defining a tria cell in grid file
  !
  np = size(rs,dim=2)             ! no. of solution points in a triangle cell
  ps = solpts2order(Geom_Tria,np) ! solution order
  !
  ! Initialize the allocatable arrays
  !
  return_value = init_array(mp,np,zero) ! F2003 AUTO-REALLOC
  !
  ! Loop over each of the solution points
  !
  nk = 0
  do nj = 0,ps
    do ni = 0,ps-nj
      !
      nk = nk + 1
      !
      ! Evaluate the Lagrange polynomial for each grid point at
      ! the reference coordinates of the current solution point
      !
      k = 0
      do j = 0,pm
        do i = 0,pm-j
          !
          ! Evaluate the Lagrange polynomial for grid point (i,j) at
          ! the reference space coordinates for solution point (ni,nj)
          !
          k = k + 1
          return_value(k,nk) = eval_triangle_lagrange(i,j,pm,rs(1:2,nk))
          !
        end do
      end do
      !
    end do
  end do
  !
end function equidistant_map_tria_to_solpts
!
!###############################################################################
!
pure &
function eval_triangle_lagrange(i,j,p,rs) result(return_value)
  !
  integer,  intent(in) :: i
  integer,  intent(in) :: j
  integer,  intent(in) :: p
  real(wp), intent(in) :: rs(:)
  !
  real(wp) :: return_value
  !
  integer  :: k,m
  real(wp) :: a,b,c
  real(wp) :: val_i,val_j,val_k
  real(wp) :: nrm_i,nrm_j,nrm_k
  real(wp) :: r,s,rn,sn,rp,ri,rj
  !
continue
  !
  a = zero
  b = one
  c = one
  !
  if (tria_is_biunit) then
    a = one
    b = two
    c = half
  end if
  !
  k = p - i - j
  !
  ! Initialize the three values to one
  !
  val_i = one
  val_j = one
  val_k = one
  !
  nrm_i = one
  nrm_j = one
  nrm_k = one
  !
  r = rs(1)
  s = rs(2)
  !
  rp = real(p,kind=wp)
  ri = real(i,kind=wp)
  rj = real(j,kind=wp)
  !
  rn = b*ri/rp - a
  sn = b*rj/rp - a
  !
  ! Compute the product series for the x/r coordinate
  !
  do m = 0,i-1
    val_i = val_i * (c*(r +a) - real(m,kind=wp)/rp)
    nrm_i = nrm_i * (c*(rn+a) - real(m,kind=wp)/rp)
  end do
  !
  ! Compute the product series for the y/s coordinate
  !
  do m = 0,j-1
    val_j = val_j * (c*(s +a) - real(m,kind=wp)/rp)
    nrm_j = nrm_j * (c*(sn+a) - real(m,kind=wp)/rp)
  end do
  !
  ! Compute the product series for the auxilary index k = p-i-j
  !
  do m = 0,k-1
    val_k = val_k * (one-c*(r +a)-c*(s +a) - real(m,kind=wp)/rp)
    nrm_k = nrm_k * (one-c*(rn+a)-c*(sn+a) - real(m,kind=wp)/rp)
  end do
  !
  ! Compute the value of the triangle lagrange at (r,s)
  !
  return_value = val_i*val_j*val_k / (nrm_i*nrm_j*nrm_k)
  !  !
end function eval_triangle_lagrange
!
!###############################################################################
!
pure &
function eval_tetrahedra_lagrange(i,j,k,p,rs) result(return_value)
  !
  integer,  intent(in) :: i
  integer,  intent(in) :: j
  integer,  intent(in) :: k
  integer,  intent(in) :: p
  real(wp), intent(in) :: rs(:)
  !
  real(wp) :: return_value
  !
  integer  :: l,m
  real(wp) :: a,b,c
  real(wp) :: val_i,val_j,val_k,val_l
  real(wp) :: nrm_i,nrm_j,nrm_k,nrm_l
  real(wp) :: r,s,t,rn,sn,tn,rp,ri,rj,rk
  !
continue
  !
  a = zero
  b = one
  c = one
  !
  if (tria_is_biunit) then
    a = one
    b = two
    c = half
  end if
  !
  l = p - i - j - k
  !
  ! Initialize the three values to one
  !
  val_i = one
  val_j = one
  val_k = one
  val_l = one
  !
  nrm_i = one
  nrm_j = one
  nrm_k = one
  nrm_l = one
  !
  r = rs(1)
  s = rs(2)
  t = rs(3)
  !
  rp = real(p,kind=wp)
  ri = real(i,kind=wp)
  rj = real(j,kind=wp)
  rk = real(k,kind=wp)
  !
  rn = b*ri/rp - a
  sn = b*rj/rp - a
  tn = b*rk/rp - a
  !
  ! Compute the product series for the x/r coordinate
  !
  do m = 0,i-1
    val_i = val_i * (c*(r +a) - real(m,kind=wp)/rp)
    nrm_i = nrm_i * (c*(rn+a) - real(m,kind=wp)/rp)
  end do
  !
  ! Compute the product series for the y/s coordinate
  !
  do m = 0,j-1
    val_j = val_j * (c*(s +a) - real(m,kind=wp)/rp)
    nrm_j = nrm_j * (c*(sn+a) - real(m,kind=wp)/rp)
  end do
  !
  ! Compute the product series for the z/t coordinate
  !
  do m = 0,k-1
    val_k = val_k * (c*(t +a) - real(m,kind=wp)/rp)
    nrm_k = nrm_k * (c*(tn+a) - real(m,kind=wp)/rp)
  end do
  !
  ! Compute the product series for the auxilary index l = p-i-j-k
  !
  do m = 0,l-1
    val_l = val_l * (one-c*(r +a)-c*(s +a)-c*(t +a) - real(m,kind=wp)/rp)
    nrm_l = nrm_l * (one-c*(rn+a)-c*(sn+a)-c*(tn+a) - real(m,kind=wp)/rp)
  end do
  !
  ! Compute the value of the tetrahedra lagrange at (r,s,t)
  !
  return_value = val_i*val_j*val_k*val_l / (nrm_i*nrm_j*nrm_k*nrm_l)
  !  !
end function eval_tetrahedra_lagrange
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function equidistant_map_quad_to_solpts(op,np,xi,mesh_order,family) &
                                 result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  use polynomial_mod, only : eval_LagrangePoly2D
  use polynomial_mod, only : equal_distant_quadrature
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: op
  integer,  intent(in) :: np
  real(wp), intent(in) :: xi(:)
  integer,  intent(in) :: mesh_order
  !
  type(geom_family_t), intent(in) :: family
  !
  !.. Function Result ..
  real(wp) :: return_value(1:op,1:np)
  !
  !.. Local Scalars ..
  integer  :: i,j,m,mi,mj,n,ni,nj,xpts,epts
  !
  !.. Local Arrays ..
  real(mp) :: xy_qp(1:2)
  real(mp) :: xi_qp(1:size(xi))
  real(mp) :: eq_pts(1:mesh_order+1)
  real(mp) :: map_qp(1:op,1:np)
  !
  !.. Local Parameters ..
  real(mp), parameter :: end_pts(1:2) = [-qone,qone]
  !
continue
  !
  xpts = size(xi_qp)
  epts = size(eq_pts)
  !
  return_value(:,:) = zero
  !
  call equal_distant_quadrature(eq_pts)
  !
  xi_qp(:) = real(xi(:),kind=kind(xi_qp))
  !
  if (family == Edge_Defined .or. family == Face_Defined) then
    !
    do mj = 1,xpts
      !
      xy_qp(2) = xi_qp(mj)
      !
      do mi = 1,xpts
        !
        xy_qp(1) = xi_qp(mi)
        !
        m = (mj-1)*xpts + mi
        !
        n = 0
        if (edge_defined_uses_endpts) then
          !
          do nj = 1,epts
            do ni = 1,epts,merge(1,epts-1,any(nj==[1,epts]))
              n = n + 1
              if (all(ni /= [1,epts])) then
                i = ni
                j = merge(2,1,nj==epts)
                map_qp(n,m) = eval_LagrangePoly2D([i,j],xy_qp,eq_pts,end_pts)
              else if (all(nj /= [1,epts])) then
                i = merge(2,1,ni==epts)
                j = nj
                map_qp(n,m) = eval_LagrangePoly2D([i,j],xy_qp,end_pts,eq_pts)
              else
                map_qp(n,m) = eval_LagrangePoly([ni,nj],xy_qp,eq_pts)
              end if
            end do
          end do
          !
        else
          !
          do nj = 1,epts
            do ni = 1,epts,merge(1,epts-1,any(nj==[1,epts]))
              n = n + 1
              map_qp(n,m) = eval_LagrangePoly([ni,nj],xy_qp,eq_pts)
            end do
          end do
          !
        end if
        !
      end do
    end do
    !
  else
    !
    do mj = 1,xpts
      !
      xy_qp(2) = xi_qp(mj)
      !
      do mi = 1,xpts
        !
        xy_qp(1) = xi_qp(mi)
        !
        m = (mj-1)*xpts + mi
        !
        do nj = 1,epts
          do ni = 1,epts
            !
            n = (nj-1)*epts + ni
            !
            map_qp(n,m) = eval_LagrangePoly([ni,nj],xy_qp,eq_pts)
            !
          end do
        end do
        !
      end do
    end do
    !
  end if
  !
  ! Use the chop function to copy map_qp into return_value while also
  ! setting any values that are less that working machine epsilon to zero.
  !
  return_value = chop( map_qp )
  !
end function equidistant_map_quad_to_solpts
!
!###############################################################################
!
#ifndef DISABLE_PURE
!pure &
#endif
function equidistant_map_tetr_to_solpts(rst,mesh_order) result(return_value)
  !
  ! RST : location of solution points on the reference tetrahedra element
  ! MESH_ORDER : polynomial order of the grid file
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: rst(:,:)
  integer,  intent(in) :: mesh_order
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,ni,nj,nk,nl
  integer :: pm,ps,mp,np
  !
continue
  !
  pm = mesh_order                 ! order of grid cells
  mp = geom_solpts(Geom_Tetr,pm)  ! no. of pts defining a tetr cell in grid file
  !
  np = size(rst,dim=2)             ! no. of solution points in a tetrahedra cell
  ps = solpts2order(Geom_Tetr,np) ! solution order
  !
  ! Initialize the allocatable arrays
  !
  return_value = init_array(mp,np,zero) ! F2003 AUTO-REALLOC
  !
  ! Loop over each of the solution points
  !
  nl = 0
  do nk = 0,ps
    do nj = 0,ps-nk
      do ni = 0,ps-nk-nj
        !
        nl = nl + 1
        !
        ! Evaluate the Lagrange polynomial for each grid point at
        ! the reference coordinates of the current solution point
        !
        l = 0
        do k = 0,pm
          do j = 0,pm-k
            do i = 0,pm-k-j
              !
              ! Evaluate the Lagrange polynomial for grid point (i,j,k) at
              ! the reference space coordinates for solution point (ni,nj,nk)
              !
              l = l + 1
              return_value(l,nl) = eval_tetrahedra_lagrange(i,j,k,pm,rst(1:3,nl))
              !
            end do
          end do
        end do
        !
      end do
    end do
  end do
  !
end function equidistant_map_tetr_to_solpts
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function equidistant_map_pyra_to_solpts(xi,rs,mesh_order) result(return_value)
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xi(:)
  real(wp), intent(in) :: rs(:,:)
  integer,  intent(in) :: mesh_order
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer :: n
  integer :: pm,ps,mp,np
  !
continue
  !
  pm = mesh_order                ! order of grid cells
  mp = geom_solpts(Geom_Pyra,pm) ! no. of pts defining a pyra cell in grid file
  !
  ps = size(xi) - 1              ! solution order
  np = geom_solpts(Geom_Pyra,ps) ! no. of solution points in a pyramid cell
  !
  return_value = init_array(mp,np,zero) ! F2003 AUTO-REALLOC
  !
end function equidistant_map_pyra_to_solpts
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function equidistant_map_pris_to_solpts(xi,rs,mesh_order) result(return_value)
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xi(:)
  real(wp), intent(in) :: rs(:,:)
  integer,  intent(in) :: mesh_order
  !
  !.. Function Result ..
  real(wp), allocatable :: return_value(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,j,k
  integer  :: pm,ps,mp,np
  real(wp) :: r,s,a,b,c,zetam1,zetap1
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(rs,1)) :: p1
  real(wp), dimension(1:size(rs,1)) :: p2
  real(wp), dimension(1:size(rs,1)) :: p3
  !
continue
  !
  pm = mesh_order                ! order of grid cells
  mp = geom_solpts(Geom_Pris,pm) ! no. of pts defining a pris cell in grid file
  !
  ps = size(xi) - 1              ! solution order
  np = geom_solpts(Geom_Pris,ps) ! no. of solution points in a prism cell
  !
  ! Initialize local allocatable array
  !
  return_value = init_array(mp,np,zero) ! F2003 AUTO-REALLOC
  !
  ! If tria_is_biunit is true:  a=0 ; b=1 ; c=1/2
  ! If tria_is_biunit is false: a=1 ; b=0 ; c=1
  !
  b = deltafun(tria_is_biunit)
  a = one - b
  c = one - half*b
  !
  k = 0
  !
  do j = 1,size(xi)
    !
    zetam1 = one - xi(j)
    zetap1 = one + xi(j)
    !
    ! Find the nodes of the triangle for the current xi-eta plane.
    !
    p1 = half * (zetam1*rs (:,1) + zetap1*rs (:,4))
    p2 = half * (zetam1*rs (:,2) + zetap1*rs (:,5))
    p3 = half * (zetam1*rs (:,3) + zetap1*rs (:,6))
    !
    do i = 1,size(rs,dim=2)
      !
      k = k + 1
      !
      r = rs(1,i)
      s = rs(2,i)
      !
      return_value(1,k) = c * ((a-r-s)*p1(1) + (r+b)*p2(1) + (s+b)*p3(1))
      return_value(2,k) = c * ((a-r-s)*p1(2) + (r+b)*p2(2) + (s+b)*p3(2))
      return_value(3,k) = xi(j)
      !
    end do
    !
  end do
  !
end function equidistant_map_pris_to_solpts
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function equidistant_map_hexa_to_solpts(op,np,xi,mesh_order,family) &
                                 result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  use polynomial_mod, only : eval_LagrangePoly3D
  use polynomial_mod, only : equal_distant_quadrature
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: op
  integer,  intent(in) :: np
  real(wp), intent(in) :: xi(:)
  integer,  intent(in) :: mesh_order
  !
  type(geom_family_t), intent(in) :: family
  !
  !.. Function Result ..
  real(wp) :: return_value(1:op,1:np)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,m,mi,mj,mk,n,ni,nj,nk,xpts,epts
  !
  !.. Local Arrays ..
  real(mp) :: xy_qp(1:3)
  real(mp) :: xi_qp(1:size(xi))
  real(mp) :: eq_pts(1:mesh_order+1)
  real(mp) :: map_qp(1:op,1:np)
  !
  !.. Local Parameters ..
  real(mp), parameter :: end_pts(1:2) = [-qone,qone]
  !
continue
  !
  xpts = size(xi_qp)
  epts = size(eq_pts)
  !
  return_value(:,:) = zero
  !
  call equal_distant_quadrature(eq_pts)
  !
  xi_qp(:) = real(xi(:),kind=kind(xi_qp))
  !
  if (family == Edge_Defined) then
    !
    do mk = 1,xpts
      !
      xy_qp(3) = xi_qp(mk)
      !
      do mj = 1,xpts
        !
        xy_qp(2) = xi_qp(mj)
        !
        do mi = 1,xpts
          !
          xy_qp(1) = xi_qp(mi)
          !
          m = (mj-1)*xpts + mi
          !
          n = 0
          if (edge_defined_uses_endpts) then
            !
            do nk = 1,epts
              do nj = 1,epts,merge(1,epts-1,any(nk==[1,epts]))
                do ni = 1,epts,merge(1,epts-1,any(nk==[1,epts]) .and. &
                                              any(nj==[1,epts]))
                  n = n + 1
                  if (all(ni /= [1,epts])) then
                    i = ni
                    j = merge(2,1,nj==epts)
                    k = merge(2,1,nk==epts)
                    map_qp(n,m) = eval_LagrangePoly3D([i,j,k],xy_qp, &
                                                      eq_pts,end_pts,end_pts)
                  else if (all(nj /= [1,epts])) then
                    i = merge(2,1,ni==epts)
                    j = nj
                    k = merge(2,1,nk==epts)
                    map_qp(n,m) = eval_LagrangePoly3D([i,j,k],xy_qp, &
                                                      end_pts,eq_pts,end_pts)
                  else if (all(nk /= [1,epts])) then
                    i = merge(2,1,ni==epts)
                    j = merge(2,1,nj==epts)
                    k = nk
                    map_qp(n,m) = eval_LagrangePoly3D([i,j,k],xy_qp, &
                                                      end_pts,end_pts,eq_pts)
                  else
                    map_qp(n,m) = eval_LagrangePoly([ni,nj,nk],xy_qp,eq_pts)
                  end if
                end do
              end do
            end do
            !
          else
            !
            do nk = 1,epts
              do nj = 1,epts,merge(1,epts-1,any(nk==[1,epts]))
                do ni = 1,epts,merge(1,epts-1,any(nk==[1,epts]) .and. &
                                              any(nj==[1,epts]))
                  n = n + 1
                  map_qp(n,m) = eval_LagrangePoly([ni,nj,nk],xy_qp,eq_pts)
                end do
              end do
            end do
            !
          end if
          !
        end do
      end do
    end do
    !
  else if (family == Face_Defined) then
    !
  else
    !
    do mk = 1,xpts
      !
      xy_qp(3) = xi_qp(mk)
      !
      do mj = 1,xpts
        !
        xy_qp(2) = xi_qp(mj)
        !
        do mi = 1,xpts
          !
          xy_qp(1) = xi_qp(mi)
          !
          m = (mk-1)*xpts*xpts + (mj-1)*xpts + mi
          !
          do nk = 1,epts
            do nj = 1,epts
              do ni = 1,epts
                !
                n = (nk-1)*epts*epts + (nj-1)*epts + ni
                !
                map_qp(n,m) = eval_LagrangePoly([ni,nj,nk],xy_qp,eq_pts)
                !
              end do
            end do
          end do
          !
        end do
      end do
    end do
    !
  end if
  !
  ! Use the chop function to copy map_qp into return_value while also
  ! setting any values that are less that working machine epsilon to zero.
  !
  return_value = chop( map_qp )
  !
end function equidistant_map_hexa_to_solpts
!
!###############################################################################
!###############################################################################
!###############################################################################
!
! CUSTOM MAPPINGS USING PARAMETERIZED EQUATIONS TO CORRECT PLOT3D AGGLOMERATIONS
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function parametric_map_edge_to_solpts(np,xi,elem_xyz) result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  use ovar,           only : parametric_agglomeration
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: np
  real(wp), intent(in) :: xi(:)
  real(wp), intent(in) :: elem_xyz(:,:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(elem_xyz,dim=2),1:np)
  !
  !.. Local Scalars ..
  integer  :: i,mi,ni,xpts,gpts
  real(mp) :: arc_len
  !
  !.. Local Arrays ..
  real(mp) :: xi_qp(1:size(xi))
  real(mp) :: gp_qp(1:size(elem_xyz,dim=2))
  real(mp) :: map_qp(1:size(elem_xyz,dim=2),1:np)
  !
continue
  !
  xpts = size(xi_qp)
  gpts = size(gp_qp)
  !
  return_value(:,:) = zero
  !
  xi_qp(:) = real(xi(:),kind=kind(xi_qp))
  !
  if (parametric_agglomeration == Cartesian_Parameterization) then
    !
    gp_qp(:) = real(elem_xyz(1,:),kind=kind(gp_qp))
    !
    arc_len = abs(gp_qp(gpts) - gp_qp(1))
    !
    if (gp_qp(gpts) > gp_qp(1)) then
      gp_qp(:) = (gp_qp(:) - gp_qp(1)) / arc_len
    else
      gp_qp(:) = (gp_qp(:) - gp_qp(gpts)) / arc_len
    end if
    !
    gp_qp(:) = qtwo*gp_qp(:) - qone
    !
  else if (parametric_agglomeration == Axisymmetric_Parameterization) then
    !
  end if
  !
  do mi = 1,xpts
    do ni = 1,gpts
      map_qp(ni,mi) = eval_LagrangePoly([ni],[xi_qp(mi)],gp_qp)
    end do
  end do
  !
  ! Use the chop function to copy map_qp into return_value while also
  ! setting any values that are less that working machine epsilon to zero.
  !
  return_value = chop( map_qp )
  !
end function parametric_map_edge_to_solpts
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function parametric_map_quad_to_solpts(np,xi,elem_xyz,mesh_order) &
                                result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly2D
  use ovar,           only : parametric_agglomeration
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: np
  real(wp), intent(in) :: xi(:)
  real(wp), intent(in) :: elem_xyz(:,:)
  integer,  intent(in) :: mesh_order
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(elem_xyz,dim=2),1:np)
  !
  !.. Local Scalars ..
  integer  :: i,j,l,m,mi,mj,n,ni,nj,xpts,gpts
  real(mp) :: arc_len
  !
  !.. Local Arrays ..
  real(mp) :: xy_qp(1:2)
  real(mp) :: xi_qp(1:size(xi))
  real(mp) :: gp_qp(1:mesh_order+1)
  real(mp) :: comp_coord(1:mesh_order+1,1:2,1:size(elem_xyz,dim=2))
  real(mp) :: map_qp(1:size(elem_xyz,dim=2),1:np)
  !
continue
  !
  xpts = size(xi_qp)
  gpts = mesh_order+1
  !
  return_value(:,:) = zero
  !
  xi_qp(:) = real(xi(:),kind=kind(xi_qp))
  !
  if (parametric_agglomeration == Cartesian_Parameterization) then
    !
    do mj = 1,gpts
      do mi = 1,gpts
        !
        m = (mj-1)*gpts + mi
        !
        do n = 1,gpts
          !
          ni = (mj-1)*gpts + n
          nj = ( n-1)*gpts + mi
          !
          ! Save the 1D x coordinates
          comp_coord(n,1,m) = real(elem_xyz(1,ni),kind=kind(comp_coord))
          ! Save the 1D y coordinates
          comp_coord(n,2,m) = real(elem_xyz(2,nj),kind=kind(comp_coord))
          !
        end do
        !
      end do
    end do
    !
    ! Rescale the 1D coordinates so they are on the interval [-1,+1]
    !
    do n = 1,size(comp_coord,dim=3)
      do l = 1,size(comp_coord,dim=2)
        !
        gp_qp(1:gpts) = comp_coord(1:gpts,l,n)
        !
        arc_len = abs(gp_qp(gpts) - gp_qp(1))
        !
        if (gp_qp(gpts) > gp_qp(1)) then
          gp_qp(:) = (gp_qp(:) - gp_qp(1)) / arc_len
        else
          gp_qp(:) = (gp_qp(:) - gp_qp(gpts)) / arc_len
        end if
        !
        comp_coord(1:gpts,l,n) = qtwo*gp_qp(1:gpts) - qone
        !
      end do
    end do
    !
  else if (parametric_agglomeration == Axisymmetric_Parameterization) then
    !
  end if
  !
  ! Create the mapping for this element
  !
  do mj = 1,xpts
    !
    xy_qp(2) = xi_qp(mj)
    !
    do mi = 1,xpts
      !
      xy_qp(1) = xi_qp(mi)
      !
      m = (mj-1)*xpts + mi
      !
      do nj = 1,gpts
        do ni = 1,gpts
          !
          n = (nj-1)*gpts + ni
          !
          map_qp(n,m) = eval_LagrangePoly2D([ni,nj],xy_qp, &
                                            comp_coord(1:gpts,1,n), &
                                            comp_coord(1:gpts,2,n))
          !
        end do
      end do
      !
    end do
  end do
  !
  ! Use the chop function to copy map_qp into return_value while also
  ! setting any values that are less that working machine epsilon to zero.
  !
  return_value = chop( map_qp )
  !
end function parametric_map_quad_to_solpts
!
!###############################################################################
!
#ifndef DISABLE_PURE
#ifndef TIME_PARAMETRIC_MAPPING
pure &
#endif
#endif
function parametric_map_hexa_to_solpts(np,xi,elem_xyz,mesh_order) &
                                result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly1D
  use polynomial_mod, only : eval_LagrangePoly3D
  use ovar,           only : parametric_agglomeration
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: np
  real(wp), intent(in) :: xi(:)
  real(wp), intent(in) :: elem_xyz(:,:)
  integer,  intent(in) :: mesh_order
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(elem_xyz,dim=2),1:np)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,m,mi,mj,mk,n,ni,nj,nk,xpts,gpts
  real(mp) :: arc_len,x1,x2,x3,rndoff_qp,d1,d2,d3
  !
  !.. Local Arrays ..
  integer  :: imask(1:mesh_order+1)
  real(mp) :: xy_qp(1:3)
  real(mp) :: xi_qp(1:size(xi))
  real(mp) :: gp_qp(1:mesh_order+1)
  real(mp) :: x1_pts(1:mesh_order+1)
  real(mp) :: x2_pts(1:mesh_order+1)
  real(mp) :: x3_pts(1:mesh_order+1)
  real(mp) :: comp_coord(1:mesh_order+1,1:3,1:size(elem_xyz,dim=2))
  real(mp) :: lp_qp(1:size(xi),1:3)
  real(mp) :: map_qp(1:size(elem_xyz,dim=2),1:np)
  !
continue
  !
#ifdef TIME_PARAMETRIC_MAPPING
  call debug_timer(start_accumulation,"parametric initialization",acc=1)
#endif
  rndoff_qp = real(rndoff,kind=kind(rndoff_qp))
  xpts = size(xi_qp)
  gpts = mesh_order+1
  imask(:) = intseq(1,mesh_order+1)
  !
  return_value(:,:) = zero
  !
  xi_qp(:) = real(xi(:),kind=kind(xi_qp))
#ifdef TIME_PARAMETRIC_MAPPING
  call debug_timer(stop_accumulation,"parametric initialization",acc=1)
#endif
  !
  if (parametric_agglomeration == Cartesian_Parameterization) then
    !
#ifdef TIME_PARAMETRIC_MAPPING
    call debug_timer(start_accumulation,"create comp_coord",acc=2)
#endif
    do mk = 1,gpts
      do mj = 1,gpts
        do mi = 1,gpts
          !
          m = (mk-1)*gpts*gpts + (mj-1)*gpts + mi
          !
          do n = 1,gpts
            !
            ni = (mk-1)*gpts*gpts + (mj-1)*gpts + n
            nj = (mk-1)*gpts*gpts + ( n-1)*gpts + mi
            nk = ( n-1)*gpts*gpts + (mj-1)*gpts + mi
            !
            ! Save the 1D x coordinates
            comp_coord(n,1,m) = real(elem_xyz(1,ni),kind=kind(comp_coord))
            ! Save the 1D y coordinates
            comp_coord(n,2,m) = real(elem_xyz(2,nj),kind=kind(comp_coord))
            ! Save the 1D z coordinates
            comp_coord(n,3,m) = real(elem_xyz(3,nk),kind=kind(comp_coord))
            !
          end do
          !
        end do
      end do
    end do
#ifdef TIME_PARAMETRIC_MAPPING
    call debug_timer(stop_accumulation,"create comp_coord",acc=2)
#endif
    !
    ! Rescale the 1D coordinates so they are on the interval [-1,+1]
    !
#ifdef TIME_PARAMETRIC_MAPPING
    call debug_timer(start_accumulation,"rescale comp_coord",acc=3)
#endif
    do n = 1,size(comp_coord,dim=3)
      do l = 1,size(comp_coord,dim=2)
        !
        gp_qp(1:gpts) = comp_coord(1:gpts,l,n)
        !
        arc_len = abs(gp_qp(gpts) - gp_qp(1))
        !
        if (gp_qp(gpts) > gp_qp(1)) then
          gp_qp(:) = (gp_qp(:) - gp_qp(1)) / arc_len
        else
          gp_qp(:) = (gp_qp(:) - gp_qp(gpts)) / arc_len
        end if
        !
        comp_coord(1:gpts,l,n) = qtwo*gp_qp(1:gpts) - qone
        !
      end do
    end do
#ifdef TIME_PARAMETRIC_MAPPING
    call debug_timer(stop_accumulation,"rescale comp_coord",acc=3)
#endif
    !
  else if (parametric_agglomeration == Axisymmetric_Parameterization) then
    !
  end if
  !
  ! Create the mapping for this element
  !
#ifdef TIME_PARAMETRIC_MAPPING
  call debug_timer(start_accumulation,"create mapping matrix",acc=4)
#endif
  do nk = 1,gpts
    do nj = 1,gpts
      do ni = 1,gpts
        !
        n = (nk-1)*gpts*gpts + (nj-1)*gpts + ni
        !
        x1_pts(1:gpts) = comp_coord(1:gpts,1,n)
        x2_pts(1:gpts) = comp_coord(1:gpts,2,n)
        x3_pts(1:gpts) = comp_coord(1:gpts,3,n)
        !
       !d1 = qone
       !d2 = qone
       !d3 = qone
       !do l = 1,gpts
       !  if (l /= ni) d1 = d1 * (x1_pts(ni) - x1_pts(l))
       !  if (l /= nj) d2 = d2 * (x2_pts(nj) - x2_pts(l))
       !  if (l /= nk) d3 = d3 * (x3_pts(nk) - x3_pts(l))
       !end do
       !d1 = qone / d1
       !d2 = qone / d2
       !d3 = qone / d3
        !
        d1 = qone / product( x1_pts(ni) - x1_pts , imask /= ni )
        d2 = qone / product( x2_pts(nj) - x2_pts , imask /= nj )
        d3 = qone / product( x3_pts(nk) - x3_pts , imask /= nk )
        !
       !lp_qp(:,1) = d1
       !lp_qp(:,2) = d2
       !lp_qp(:,3) = d3
       !do m = 1,xpts
       !  do l = 1,gpts
       !    if (l /= ni) lp_qp(m,1) = lp_qp(m,1) * (xi_qp(m) - x1_pts(l))
       !    if (l /= nj) lp_qp(m,2) = lp_qp(m,2) * (xi_qp(m) - x2_pts(l))
       !    if (l /= nk) lp_qp(m,3) = lp_qp(m,3) * (xi_qp(m) - x3_pts(l))
       !  end do
       !end do
        do m = 1,xpts
          lp_qp(m,1) = product(  xi_qp(m) - x1_pts , imask /= ni ) * d1
          lp_qp(m,2) = product(  xi_qp(m) - x2_pts , imask /= nj ) * d2
          lp_qp(m,3) = product(  xi_qp(m) - x3_pts , imask /= nk ) * d3
        end do
        !
        do mk = 1,xpts
          do mj = 1,xpts
            do mi = 1,xpts
              !
              m = (mk-1)*xpts*xpts + (mj-1)*xpts + mi
              !
              x1 = lp_qp(mi,1)
              x2 = lp_qp(mj,2)
              x3 = lp_qp(mk,3)
              !
              map_qp(n,m) = x1*x2*x3
              !
             !map_qp(n,m) = eval_LagrangePoly3D([ni,nj,nk],xy_qp, &
             !                                  x1_pts,x2_pts,x3_pts)
              !
            end do
          end do
        end do
        !
      end do
    end do
  end do
#ifdef TIME_PARAMETRIC_MAPPING
  call debug_timer(stop_accumulation,"create mapping matrix",acc=4)
#endif
  !
  ! Use the chop function to copy map_qp into return_value while also
  ! setting any values that are less that working machine epsilon to zero.
  !
#ifdef TIME_PARAMETRIC_MAPPING
  call debug_timer(start_accumulation,"chopping mapping matrix",acc=5)
#endif
  do m = 1,size(return_value,dim=2)
    do n = 1,size(return_value,dim=1)
      if (abs(map_qp(n,m)) > rndoff_qp) then
        return_value(n,m) = real(map_qp(n,m),kind=wp)
      end if
    end do
  end do
 !return_value = chop( map_qp )
#ifdef TIME_PARAMETRIC_MAPPING
  call debug_timer(stop_accumulation,"chopping mapping matrix",acc=5)
#endif
  !
end function parametric_map_hexa_to_solpts
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
! OLD MAPPING FUNCTIONS USING ONLY THE NODES OF THE CELLS
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_solpts_edge(xi,xy) result(return_value)
  !
  real(wp), dimension(:),   intent(in) :: xi
  real(wp), dimension(:,:), intent(in) :: xy
  !
  real(wp), dimension(1:size(xy,1),1:size(xi)) :: return_value
  !
  integer  :: i
  real(wp) :: xim1,xip1
  !
continue
  !
  do i = 1,size(xi)
    !
    xim1 = one - xi(i)
    xip1 = one + xi(i)
    !
    return_value(:,i) = half * ( xim1*xy(:,1) + xip1*xy(:,2) )
    !
  end do
  !
end function get_solpts_edge
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_solpts_tria(rs,xy) result(return_value)
  !
  real(wp), dimension(:,:), intent(in) :: rs
  real(wp), dimension(:,:), intent(in) :: xy
  !
  real(wp), dimension(1:size(rs,1),1:size(rs,2)) :: return_value
  !
  integer  :: i
  real(wp) :: r,s,a,b,c
  !
continue
  !
  ! If tria_is_biunit is true:  a=0 ; b=1 ; c=1/2
  ! If tria_is_biunit is false: a=1 ; b=0 ; c=1
  !
  b = deltafun(tria_is_biunit)
  a = one - b
  c = one - half*b
  !
  ! On entry, r and s both need to be defined on the interval [0,1]
  !
  do i = 1,size(rs,2)
    !
    r = rs(1,i)
    s = rs(2,i)
    !
    return_value(:,i) = c * ((a-r-s)*xy(:,1) + (r+b)*xy(:,2) + (s+b)*xy(:,3))
    !
  end do
  !
end function get_solpts_tria
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_solpts_quad(xi,xy) result(return_value)
  !
  real(wp), dimension(:),   intent(in) :: xi
  real(wp), dimension(:,:), intent(in) :: xy
  !
  real(wp), dimension(1:size(xy,1),1:size(xi)**2) :: return_value
  !
  integer  :: i,j,l
  real(wp) :: etam1,etap1,xim1,xip1
  !
continue
  !
  l = 1
  do j = 1,size(xi)
    !
    etam1 = one - xi(j)
    etap1 = one + xi(j)
    !
    do i = 1,size(xi)
      !
      xim1 = one - xi(i)
      xip1 = one + xi(i)
      !
      return_value(:,l) = one4 * ( etam1*xim1*xy(:,1) + &
                                   etam1*xip1*xy(:,2) + &
                                   etap1*xip1*xy(:,3) + &
                                   etap1*xim1*xy(:,4) )
      !
      l = l + 1
      !
    end do
    !
  end do
  !
end function get_solpts_quad
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_solpts_tetr(rs,xyz) result(return_value)
  !
  real(wp), dimension(:,:), intent(in) :: rs
  real(wp), dimension(:,:), intent(in) :: xyz
  !
  real(wp), allocatable :: return_value(:,:)
  !
  integer :: ndim,npts,ierr
  integer :: n
  !
continue
  !
  n = size(rs,dim=2)
  !
  ndim = size(xyz,dim=1)
  npts = n*(3+nint(sqrt(real(1+8*n,kind=wp))))/6
  !
  allocate ( return_value(1:ndim,1:npts) , source=zero , stat=ierr )
  !
end function get_solpts_tetr
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_solpts_pyra(xi,rs,xyz) result(return_value)
  !
  real(wp), dimension(:),   intent(in) :: xi
  real(wp), dimension(:,:), intent(in) :: rs
  real(wp), dimension(:,:), intent(in) :: xyz
  !
  real(wp), allocatable :: return_value(:,:)
  !
  integer :: ndim,npts,ierr
  integer :: n
  !
continue
  !
  n = size(xi)-1
  !
  ndim = size(xyz,dim=1)
  npts = (n+1)*(n+2)*(2*n+3)/6
  !
  allocate ( return_value(1:ndim,1:npts) , source=zero , stat=ierr )
  !
end function get_solpts_pyra
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_solpts_pris(xi,rs,xyz) result(return_value)
  !
  real(wp), dimension(:),   intent(in) :: xi
  real(wp), dimension(:,:), intent(in) :: rs
  real(wp), dimension(:,:), intent(in) :: xyz
  !
  real(wp), allocatable :: return_value(:,:)
  !
  integer  :: ndim,npts,ierr
  integer  :: i,j,k
  real(wp) :: r,s,a,b,c,zetam1,zetap1
  !
  real(wp), dimension(1:size(xyz,1)) :: p1
  real(wp), dimension(1:size(xyz,1)) :: p2
  real(wp), dimension(1:size(xyz,1)) :: p3
  !
continue
  !
  ndim = size(xyz,dim=1)
  npts = size(xi)*size(rs,dim=2)
  !
  allocate ( return_value(1:ndim,1:npts) , source=zero , stat=ierr )
  !
  ! If tria_is_biunit is true:  a=0 ; b=1 ; c=1/2
  ! If tria_is_biunit is false: a=1 ; b=0 ; c=1
  !
  b = deltafun(tria_is_biunit)
  a = one - b
  c = one - half*b
  !
  k = 0
  !
  do j = 1,size(xi)
    !
    zetam1 = one - xi(j)
    zetap1 = one + xi(j)
    !
    ! Find the nodes of the triangle for the current xi-eta plane.
    !
    p1 = half * (zetam1*xyz(:,1) + zetap1*xyz(:,4))
    p2 = half * (zetam1*xyz(:,2) + zetap1*xyz(:,5))
    p3 = half * (zetam1*xyz(:,3) + zetap1*xyz(:,6))
    !
    do i = 1,size(rs,dim=2)
      !
      k = k + 1
      !
      r = rs(1,i)
      s = rs(2,i)
      !
      return_value(:,k) = c * ((a-r-s)*p1 + (r+b)*p2 + (s+b)*p3)
      !
    end do
    !
  end do
  !
end function get_solpts_pris
!
!###############################################################################
!
#ifndef DISABLE_PURE
pure &
#endif
function get_solpts_hexa(xi,xyz) result(return_value)
  !
  real(wp), dimension(:),   intent(in) :: xi
  real(wp), dimension(:,:), intent(in) :: xyz
  !
  real(wp), dimension(1:size(xyz,1),1:size(xi)**3) :: return_value
  !
  integer  :: i,j,k,l
  real(wp) :: xim1,xip1,etam1,etap1,zetam1,zetap1
  !
continue
  !
  l = 1
  do k = 1,size(xi)
    !
    zetam1 = one - xi(k)
    zetap1 = one + xi(k)
    !
    do j = 1,size(xi)
      !
      etam1 = one - xi(j)
      etap1 = one + xi(j)
      !
      do i = 1,size(xi)
        !
        xim1 = one - xi(i)
        xip1 = one + xi(i)
        !
        return_value(:,l) = one8 * ( xim1*etam1*zetam1*xyz(:,1) + &
                                     xip1*etam1*zetam1*xyz(:,2) + &
                                     xip1*etap1*zetam1*xyz(:,3) + &
                                     xim1*etap1*zetam1*xyz(:,4) + &
                                     xim1*etam1*zetap1*xyz(:,5) + &
                                     xip1*etam1*zetap1*xyz(:,6) + &
                                     xip1*etap1*zetap1*xyz(:,7) + &
                                     xim1*etap1*zetap1*xyz(:,8) )
        !
        l = l + 1
        !
      end do
      !
    end do
    !
  end do
  !
end function get_solpts_hexa
!
!###############################################################################
!
subroutine move_mapping_to_device()
  !
  !.. Local Scalars ..
  integer :: i,i1,i2,j,j1,j2,k,k1,k2,m1,m2,n1,n2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "move_mapping_to_device"
  !
continue
  !
  if (allocated(mapping)) then
    !$ACC ENTER DATA COPYIN(mapping)
    i1 = lbound(mapping,dim=1)
    i2 = ubound(mapping,dim=1)
    j1 = lbound(mapping,dim=2)
    j2 = ubound(mapping,dim=2)
    !$ACC ENTER DATA COPYIN(mapping(i1:i2,j1:j2))
    do j = j1,j2
      do i = i1,i2
        !
        !---------------------------------------------------------------
        !
        if (allocated(mapping(i,j)%MeshOrder)) then
          !$ACC ENTER DATA COPYIN(mapping(i,j)%MeshOrder)
          k1 = lbound(mapping(i,j)%MeshOrder,dim=1)
          k2 = ubound(mapping(i,j)%MeshOrder,dim=1)
          !$ACC ENTER DATA COPYIN(mapping(i,j)%MeshOrder(k1:k2))
          do k = k1,k2
            !
            if (allocated(mapping(i,j)%MeshOrder(k)%mat)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%MeshOrder(k)%mat)
              m1 = lbound(mapping(i,j)%MeshOrder(k)%mat,dim=1)
              m2 = ubound(mapping(i,j)%MeshOrder(k)%mat,dim=1)
              n1 = lbound(mapping(i,j)%MeshOrder(k)%mat,dim=2)
              n2 = ubound(mapping(i,j)%MeshOrder(k)%mat,dim=2)
             !$ACC ENTER DATA COPYIN(mapping(i,j)%MeshOrder(k)%mat(m1:m2,n1:n2))
            end if
            !
            if (allocated(mapping(i,j)%MeshOrder(k)%row_ptr)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%MeshOrder(k)%row_ptr)
              m1 = lbound(mapping(i,j)%MeshOrder(k)%row_ptr,dim=1)
              m2 = ubound(mapping(i,j)%MeshOrder(k)%row_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(mapping(i,j)%MeshOrder(k)%row_ptr(m1:m2))
            end if
            !
            if (allocated(mapping(i,j)%MeshOrder(k)%row_idx)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%MeshOrder(k)%row_idx)
              m1 = lbound(mapping(i,j)%MeshOrder(k)%row_idx,dim=1)
              m2 = ubound(mapping(i,j)%MeshOrder(k)%row_idx,dim=1)
              !$ACC ENTER DATA COPYIN(mapping(i,j)%MeshOrder(k)%row_idx(m1:m2))
            end if
            !
            if (allocated(mapping(i,j)%MeshOrder(k)%row_val)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%MeshOrder(k)%row_val)
              m1 = lbound(mapping(i,j)%MeshOrder(k)%row_val,dim=1)
              m2 = ubound(mapping(i,j)%MeshOrder(k)%row_val,dim=1)
              !$ACC ENTER DATA COPYIN(mapping(i,j)%MeshOrder(k)%row_val(m1:m2))
            end if
            !
            if (allocated(mapping(i,j)%MeshOrder(k)%col_ptr)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%MeshOrder(k)%col_ptr)
              m1 = lbound(mapping(i,j)%MeshOrder(k)%col_ptr,dim=1)
              m2 = ubound(mapping(i,j)%MeshOrder(k)%col_ptr,dim=1)
              !$ACC ENTER DATA COPYIN(mapping(i,j)%MeshOrder(k)%col_ptr(m1:m2))
            end if
            !
            if (allocated(mapping(i,j)%MeshOrder(k)%col_idx)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%MeshOrder(k)%col_idx)
              m1 = lbound(mapping(i,j)%MeshOrder(k)%col_idx,dim=1)
              m2 = ubound(mapping(i,j)%MeshOrder(k)%col_idx,dim=1)
              !$ACC ENTER DATA COPYIN(mapping(i,j)%MeshOrder(k)%col_idx(m1:m2))
            end if
            !
            if (allocated(mapping(i,j)%MeshOrder(k)%col_val)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%MeshOrder(k)%col_val)
              m1 = lbound(mapping(i,j)%MeshOrder(k)%col_val,dim=1)
              m2 = ubound(mapping(i,j)%MeshOrder(k)%col_val,dim=1)
              !$ACC ENTER DATA COPYIN(mapping(i,j)%MeshOrder(k)%col_val(m1:m2))
            end if
            !
          end do
          !
        end if
        !
        !---------------------------------------------------------------
        !
        if (allocated(mapping(i,j)%EdgeDefOrder)) then
          !$ACC ENTER DATA COPYIN(mapping(i,j)%EdgeDefOrder)
          k1 = lbound(mapping(i,j)%EdgeDefOrder,dim=1)
          k2 = ubound(mapping(i,j)%EdgeDefOrder,dim=1)
          !$ACC ENTER DATA COPYIN(mapping(i,j)%EdgeDefOrder(k1:k2))
          do k = k1,k2
            !
            if (allocated(mapping(i,j)%EdgeDefOrder(k)%mat)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%EdgeDefOrder(k)%mat)
              m1 = lbound(mapping(i,j)%EdgeDefOrder(k)%mat,dim=1)
              m2 = ubound(mapping(i,j)%EdgeDefOrder(k)%mat,dim=1)
              n1 = lbound(mapping(i,j)%EdgeDefOrder(k)%mat,dim=2)
              n2 = ubound(mapping(i,j)%EdgeDefOrder(k)%mat,dim=2)
          !$ACC ENTER DATA COPYIN(mapping(i,j)%EdgeDefOrder(k)%mat(m1:m2,n1:n2))
            end if
            !
            if (allocated(mapping(i,j)%EdgeDefOrder(k)%row_ptr)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%EdgeDefOrder(k)%row_ptr)
              m1 = lbound(mapping(i,j)%EdgeDefOrder(k)%row_ptr,dim=1)
              m2 = ubound(mapping(i,j)%EdgeDefOrder(k)%row_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(mapping(i,j)%EdgeDefOrder(k)%row_ptr(m1:m2))
            end if
            !
            if (allocated(mapping(i,j)%EdgeDefOrder(k)%row_idx)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%EdgeDefOrder(k)%row_idx)
              m1 = lbound(mapping(i,j)%EdgeDefOrder(k)%row_idx,dim=1)
              m2 = ubound(mapping(i,j)%EdgeDefOrder(k)%row_idx,dim=1)
            !$ACC ENTER DATA COPYIN(mapping(i,j)%EdgeDefOrder(k)%row_idx(m1:m2))
            end if
            !
            if (allocated(mapping(i,j)%EdgeDefOrder(k)%row_val)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%EdgeDefOrder(k)%row_val)
              m1 = lbound(mapping(i,j)%EdgeDefOrder(k)%row_val,dim=1)
              m2 = ubound(mapping(i,j)%EdgeDefOrder(k)%row_val,dim=1)
            !$ACC ENTER DATA COPYIN(mapping(i,j)%EdgeDefOrder(k)%row_val(m1:m2))
            end if
            !
            if (allocated(mapping(i,j)%EdgeDefOrder(k)%col_ptr)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%EdgeDefOrder(k)%col_ptr)
              m1 = lbound(mapping(i,j)%EdgeDefOrder(k)%col_ptr,dim=1)
              m2 = ubound(mapping(i,j)%EdgeDefOrder(k)%col_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(mapping(i,j)%EdgeDefOrder(k)%col_ptr(m1:m2))
            end if
            !
            if (allocated(mapping(i,j)%EdgeDefOrder(k)%col_idx)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%EdgeDefOrder(k)%col_idx)
              m1 = lbound(mapping(i,j)%EdgeDefOrder(k)%col_idx,dim=1)
              m2 = ubound(mapping(i,j)%EdgeDefOrder(k)%col_idx,dim=1)
            !$ACC ENTER DATA COPYIN(mapping(i,j)%EdgeDefOrder(k)%col_idx(m1:m2))
            end if
            !
            if (allocated(mapping(i,j)%EdgeDefOrder(k)%col_val)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%EdgeDefOrder(k)%col_val)
              m1 = lbound(mapping(i,j)%EdgeDefOrder(k)%col_val,dim=1)
              m2 = ubound(mapping(i,j)%EdgeDefOrder(k)%col_val,dim=1)
            !$ACC ENTER DATA COPYIN(mapping(i,j)%EdgeDefOrder(k)%col_val(m1:m2))
            end if
            !
          end do
          !
        end if
        !
        !---------------------------------------------------------------
        !
        if (allocated(mapping(i,j)%FaceDefOrder)) then
          !$ACC ENTER DATA COPYIN(mapping(i,j)%FaceDefOrder)
          k1 = lbound(mapping(i,j)%FaceDefOrder,dim=1)
          k2 = ubound(mapping(i,j)%FaceDefOrder,dim=1)
          !$ACC ENTER DATA COPYIN(mapping(i,j)%FaceDefOrder(k1:k2))
          do k = k1,k2
            !
            if (allocated(mapping(i,j)%FaceDefOrder(k)%mat)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%FaceDefOrder(k)%mat)
              m1 = lbound(mapping(i,j)%FaceDefOrder(k)%mat,dim=1)
              m2 = ubound(mapping(i,j)%FaceDefOrder(k)%mat,dim=1)
              n1 = lbound(mapping(i,j)%FaceDefOrder(k)%mat,dim=2)
              n2 = ubound(mapping(i,j)%FaceDefOrder(k)%mat,dim=2)
          !$ACC ENTER DATA COPYIN(mapping(i,j)%FaceDefOrder(k)%mat(m1:m2,n1:n2))
            end if
            !
            if (allocated(mapping(i,j)%FaceDefOrder(k)%row_ptr)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%FaceDefOrder(k)%row_ptr)
              m1 = lbound(mapping(i,j)%FaceDefOrder(k)%row_ptr,dim=1)
              m2 = ubound(mapping(i,j)%FaceDefOrder(k)%row_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(mapping(i,j)%FaceDefOrder(k)%row_ptr(m1:m2))
            end if
            !
            if (allocated(mapping(i,j)%FaceDefOrder(k)%row_idx)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%FaceDefOrder(k)%row_idx)
              m1 = lbound(mapping(i,j)%FaceDefOrder(k)%row_idx,dim=1)
              m2 = ubound(mapping(i,j)%FaceDefOrder(k)%row_idx,dim=1)
            !$ACC ENTER DATA COPYIN(mapping(i,j)%FaceDefOrder(k)%row_idx(m1:m2))
            end if
            !
            if (allocated(mapping(i,j)%FaceDefOrder(k)%row_val)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%FaceDefOrder(k)%row_val)
              m1 = lbound(mapping(i,j)%FaceDefOrder(k)%row_val,dim=1)
              m2 = ubound(mapping(i,j)%FaceDefOrder(k)%row_val,dim=1)
            !$ACC ENTER DATA COPYIN(mapping(i,j)%FaceDefOrder(k)%row_val(m1:m2))
            end if
            !
            if (allocated(mapping(i,j)%FaceDefOrder(k)%col_ptr)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%FaceDefOrder(k)%col_ptr)
              m1 = lbound(mapping(i,j)%FaceDefOrder(k)%col_ptr,dim=1)
              m2 = ubound(mapping(i,j)%FaceDefOrder(k)%col_ptr,dim=1)
            !$ACC ENTER DATA COPYIN(mapping(i,j)%FaceDefOrder(k)%col_ptr(m1:m2))
            end if
            !
            if (allocated(mapping(i,j)%FaceDefOrder(k)%col_idx)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%FaceDefOrder(k)%col_idx)
              m1 = lbound(mapping(i,j)%FaceDefOrder(k)%col_idx,dim=1)
              m2 = ubound(mapping(i,j)%FaceDefOrder(k)%col_idx,dim=1)
            !$ACC ENTER DATA COPYIN(mapping(i,j)%FaceDefOrder(k)%col_idx(m1:m2))
            end if
            !
            if (allocated(mapping(i,j)%FaceDefOrder(k)%col_val)) then
              !$ACC ENTER DATA COPYIN(mapping(i,j)%FaceDefOrder(k)%col_val)
              m1 = lbound(mapping(i,j)%FaceDefOrder(k)%col_val,dim=1)
              m2 = ubound(mapping(i,j)%FaceDefOrder(k)%col_val,dim=1)
            !$ACC ENTER DATA COPYIN(mapping(i,j)%FaceDefOrder(k)%col_val(m1:m2))
            end if
            !
          end do
          !
        end if
        !
        !---------------------------------------------------------------
        !
      end do
    end do
  end if
  !
end subroutine move_mapping_to_device
!
!###############################################################################
!
include "Functions/solpts2order.f90"
!
!###############################################################################
!
end module mappings_mod
!
!###############################################################################
!
