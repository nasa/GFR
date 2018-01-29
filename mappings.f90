module mappings_mod
  !
  !.. Global Use Statements ..
  use module_kind_types
  use generic_types_mod, only : matrix
  !
  implicit none
  !
  private
  !
  !
  ! ##################################
  ! ###  Module Public Procedures  ###
  ! ##################################
  !
  public :: get_solution_points
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
contains
!
!###############################################################################
!
subroutine get_solution_points()
  !
  use geovar,         only : grid
  use geovar,         only : nr,nfbnd,ncell,nnode,n_totcel,n_totnod
  use geovar,         only : n_solpts,n_totpts
  use geovar,         only : nodes_of_cell_ptr,nodes_of_cell
  use geovar,         only : bface,xyz_nodes,xyz
  use geovar,         only : cell_geom,cell_order
  !
  use order_mod,      only : geom_solpts
  use order_mod,      only : n_min_geom,n_max_geom
  use order_mod,      only : n_min_order,n_max_order
  !
  use quadrature_mod, only : std_elem
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,np,mp,p1,p2,ierr
  integer :: this_cell,mesh_order
  integer :: this_geom,this_order
  logical(lk) :: use_nodes_of_cell
  logical(lk) :: check_grid_for_pts
  !
  type(geom_family_t) :: elem_family
  !
  character(len=200) :: array_name
  !
  !.. Local Pointers ..
  type(matrix), pointer             :: this_map
 !real(wp),     pointer, contiguous :: edge_pts(:)
  real(wp),     pointer :: edge_pts(:)
 !real(wp),     pointer, contiguous :: tria_pts(:,:)
  real(wp),     pointer :: tria_pts(:,:)
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: elem_xyz(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_solution_points"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  check_grid_for_pts = fals
  !
  ! Initialize the local pointers
  !
  edge_pts => null()
  tria_pts => null()
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
  allocate ( xyz(1:nr,1:n_solpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz",1,__LINE__,__FILE__,ierr,error_message)
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
  do this_cell = 1,ncell
    !
    this_geom = cell_geom(this_cell)
    this_order = cell_order(this_cell)
    !
    np = geom_solpts(this_geom,this_order) ! number of solution points
    !
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
    !
    ! Assign a local pointer to the mapping
    ! matrix as an alias to simplify the code
    !
    if (elem_family == Edge_Defined) then
      this_map => mapping(this_geom,this_order)%EdgeDefOrder(mesh_order)
      write (array_name,1) Geom_Name(this_geom),this_order, &
                           "EdgeDefOrder",mesh_order
    else if (elem_family == Face_Defined) then
      this_map => mapping(this_geom,this_order)%FaceDefOrder(mesh_order)
      write (array_name,1) Geom_Name(this_geom),this_order, &
                           "FaceDefOrder",mesh_order
    else
      this_map => mapping(this_geom,this_order)%MeshOrder(mesh_order)
      write (array_name,1) Geom_Name(this_geom),this_order, &
                           "MeshOrder",mesh_order
    end if
    !
    ! Create the mapping matrix for this combination of cell geometry and
    ! cell order if the mat component of this_map has not been allocated yet
    !
    if (.not.allocated(this_map%mat)) then
      !
      mp = size(elem_xyz,dim=2)
      !
      ! Allocate the mat component of this_map
      !
      allocate ( this_map%mat(1:mp,1:np) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      edge_pts => std_elem(Geom_Edge,this_order)%pts(1,:)
      tria_pts => std_elem(Geom_Tria,this_order)%pts
      !
      select case (this_geom)
        !
        case (Geom_Edge)
          !
          this_map%mat = map_edge_to_solpts( mp , np , edge_pts , mesh_order )
          !
        case (Geom_Tria)
          !
          this_map%mat = map_tria_to_solpts( mp , np , tria_pts , &
                                             mesh_order , elem_family )
          !
        case (Geom_Quad)
          !
          this_map%mat = map_quad_to_solpts( mp , np , edge_pts , &
                                             mesh_order , elem_family )
          !
        case (Geom_Tetr)
          !
          write (error_message,3) "Tetrahedral"
          call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
          !
          this_map%mat = map_tetr_to_solpts( mp , np , tria_pts , &
                                             mesh_order , elem_family )
          !
        case (Geom_Pyra)
          !
          write (error_message,3) "Pyramid"
          call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
          !
          this_map%mat = map_pyra_to_solpts( mp , np , edge_pts , tria_pts , &
                                             mesh_order , elem_family )
          !
        case (Geom_Pris)
          !
          write (error_message,3) "Prism"
          call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
          !
          this_map%mat = map_pris_to_solpts( mp , np , edge_pts , tria_pts , &
                                             mesh_order , elem_family )
          !
        case (Geom_Hexa)
          !
          this_map%mat = map_hexa_to_solpts( mp , np , edge_pts , &
                                             mesh_order , elem_family )
          !
        case default
          !
          write (iout,2) this_cell,cell_geom(this_cell)
          call stop_gfr(abort,pname,__LINE__,__FILE__)
          !
      end select
      !
      if (associated(edge_pts)) edge_pts => null()
      if (associated(tria_pts)) tria_pts => null()
      !
      ! Check the newly made matrix for sparsity
      !
      call this_map%check_for_sparsity
      !
    end if
    !
    ! Now use the mapping matrix to map the coordinates of the points for the
    ! element provided by the grid file to the coordinates of the solution
    ! points for the cell
    !
    xyz(1:nr,n1:n2) = this_map%mm_mult( elem_xyz(:,:) )
    !
    if (associated(this_map)) this_map => null()
    !
  end do
  !
  ! Deallocate the local arrays before leaving
  !
  if (allocated(elem_xyz)) then
    deallocate ( elem_xyz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"elem_xyz",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(mapping)) then
    deallocate ( mapping , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"mapping",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
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
end subroutine get_solution_points
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
pure function map_edge_to_solpts(mp,np,xi,mesh_order) &
                          result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: mp
  integer,  intent(in) :: np
  real(wp), intent(in) :: xi(:)
  integer,  intent(in) :: mesh_order
  !
  !.. Function Result ..
  real(wp) :: return_value(1:mp,1:np)
  !
  !.. Local Scalars ..
  integer  :: i,mi,ni,xpts,epts
  !
  !.. Local Arrays ..
  real(qp) :: xi_qp(1:size(xi))
  real(qp) :: eq_pts(1:mesh_order+1)
  real(qp) :: map_qp(1:mp,1:np)
  !
continue
  !
  xpts = size(xi_qp)
  epts = size(eq_pts)
  !
  return_value(:,:) = zero
  !
  eq_pts(:) = real(intseq(-mesh_order,mesh_order,2),kind=qp) / &
              real(max(mesh_order,1),kind=qp)
  !
  xi_qp(:) = real(xi(:),kind=qp)
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
end function map_edge_to_solpts
!
!###############################################################################
!
pure function map_tria_to_solpts(mp,np,rs,mesh_order,family) &
                          result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  use polynomial_mod, only : eval_LagrangePoly2D
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: mp
  integer,  intent(in) :: np
  real(wp), intent(in) :: rs(:,:)
  integer,  intent(in) :: mesh_order
  !
  type(geom_family_t), intent(in) :: family
  !
  !.. Function Result ..
  real(wp) :: return_value(1:mp,1:np)
  !
  !.. Local Scalars ..
  integer  :: i
  real(wp) :: r,s,a,b,c
  !
  !.. Local Arrays ..
  real(qp) :: eq_pts(1:mesh_order+1)
  real(qp) :: map_qp(1:mp,1:np)
  !
  !.. Local Parameters ..
  real(qp), parameter :: end_pts(1:2) = [-qone,qone]
  !
continue
  !
  return_value(:,:) = zero
  !
  eq_pts(:) = real(intseq(-mesh_order,mesh_order,2),kind=qp) / &
              real(max(mesh_order,1),kind=qp)
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
    return_value(:,i) = c * ((a-r-s)*rs(:,1) + (r+b)*rs(:,2) + (s+b)*rs(:,3))
    !
  end do
  !
  ! Use the chop function to copy map_qp into return_value while also
  ! setting any values that are less that working machine epsilon to zero.
  !
  return_value = chop( map_qp )
  !
end function map_tria_to_solpts
!
!###############################################################################
!
pure function map_quad_to_solpts(mp,np,xi,mesh_order,family) &
                          result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  use polynomial_mod, only : eval_LagrangePoly2D
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: mp
  integer,  intent(in) :: np
  real(wp), intent(in) :: xi(:)
  integer,  intent(in) :: mesh_order
  !
  type(geom_family_t), intent(in) :: family
  !
  !.. Function Result ..
  real(wp) :: return_value(1:mp,1:np)
  !
  !.. Local Scalars ..
  integer  :: i,j,m,mi,mj,n,ni,nj,xpts,epts
  !
  !.. Local Arrays ..
  real(qp) :: xy_qp(1:2)
  real(qp) :: xi_qp(1:size(xi))
  real(qp) :: eq_pts(1:mesh_order+1)
  real(qp) :: map_qp(1:mp,1:np)
  !
  !.. Local Parameters ..
  real(qp), parameter :: end_pts(1:2) = [-qone,qone]
  !
continue
  !
  xpts = size(xi_qp)
  epts = size(eq_pts)
  !
  return_value(:,:) = zero
  !
  eq_pts(:) = real(intseq(-mesh_order,mesh_order,2),kind=qp) / &
              real(max(mesh_order,1),kind=qp)
  !
  xi_qp(:) = real(xi(:),kind=qp)
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
end function map_quad_to_solpts
!
!###############################################################################
!
pure function map_tetr_to_solpts(mp,np,rs,mesh_order,family) &
                          result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  use polynomial_mod, only : eval_LagrangePoly3D
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: mp
  integer,  intent(in) :: np
  real(wp), intent(in) :: rs(:,:)
  integer,  intent(in) :: mesh_order
  !
  type(geom_family_t), intent(in) :: family
  !
  !.. Function Result ..
  real(wp) :: return_value(1:mp,1:np)
  !
  !.. Local Scalars ..
  integer :: n,ndim,npts
  !
  !.. Local Arrays ..
  real(qp) :: eq_pts(1:mesh_order+1)
  real(qp) :: map_qp(1:mp,1:np)
  real(qp) :: xi_qp(1:size(rs,dim=2))
  !
  !.. Local Parameters ..
  real(qp), parameter :: end_pts(1:2) = [-qone,qone]
  !
continue
  !
  return_value(:,:) = zero
  !
  eq_pts(:) = real(intseq(-mesh_order,mesh_order,2),kind=qp) / &
              real(max(mesh_order,1),kind=qp)
  !
  xi_qp(:) = real(rs(1,:),kind=qp)
  !
  n = size(rs,dim=2)
  !
 !ndim = size(xyz,dim=1)
  npts = n*(3+nint(sqrt(real(1+8*n,kind=wp))))/6
  !
  ! Use the chop function to copy map_qp into return_value while also
  ! setting any values that are less that working machine epsilon to zero.
  !
  return_value = chop( map_qp )
  !
end function map_tetr_to_solpts
!
!###############################################################################
!
pure function map_pyra_to_solpts(mp,np,xi,rs,mesh_order,family) &
                          result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  use polynomial_mod, only : eval_LagrangePoly3D
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: mp
  integer,  intent(in) :: np
  real(wp), intent(in) :: xi(:)
  real(wp), intent(in) :: rs(:,:)
  integer,  intent(in) :: mesh_order
  !
  type(geom_family_t), intent(in) :: family
  !
  !.. Function Result ..
  real(wp) :: return_value(1:mp,1:np)
  !
  !.. Local Scalars ..
  integer :: n,ndim,npts
  !
  !.. Local Arrays ..
  real(qp) :: xi_qp(1:size(xi))
  real(qp) :: eq_pts(1:mesh_order+1)
  real(qp) :: map_qp(1:mp,1:np)
  !
  !.. Local Parameters ..
  real(qp), parameter :: end_pts(1:2) = [-qone,qone]
  !
continue
  !
  return_value(:,:) = zero
  !
  eq_pts(:) = real(intseq(-mesh_order,mesh_order,2),kind=qp) / &
              real(max(mesh_order,1),kind=qp)
  !
  xi_qp(:) = real(xi(:),kind=qp)
  !
  n = size(xi)-1
  !
 !ndim = size(xyz,dim=1)
  npts = (n+1)*(n+2)*(2*n+3)/6
  !
  ! Use the chop function to copy map_qp into return_value while also
  ! setting any values that are less that working machine epsilon to zero.
  !
  return_value = chop( map_qp )
  !
end function map_pyra_to_solpts
!
!###############################################################################
!
pure function map_pris_to_solpts(mp,np,xi,rs,mesh_order,family) &
                          result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  use polynomial_mod, only : eval_LagrangePoly3D
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: mp
  integer,  intent(in) :: np
  real(wp), intent(in) :: xi(:)
  real(wp), intent(in) :: rs(:,:)
  integer,  intent(in) :: mesh_order
  !
  type(geom_family_t), intent(in) :: family
  !
  !.. Function Result ..
  real(wp) :: return_value(1:mp,1:np)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,ndim,npts
  real(wp) :: r,s,a,b,c,zetam1,zetap1
  !
  !.. Local Arrays ..
  real(qp) :: xi_qp(1:size(xi))
  real(qp) :: eq_pts(1:mesh_order+1)
  real(qp) :: map_qp(1:mp,1:np)
  !
  !.. Local Arrays ..
  real(wp), dimension(1:size(rs ,1)) :: p1
  real(wp), dimension(1:size(rs ,1)) :: p2
  real(wp), dimension(1:size(rs ,1)) :: p3
  !
  !.. Local Parameters ..
  real(qp), parameter :: end_pts(1:2) = [-qone,qone]
  !
continue
  !
  return_value(:,:) = zero
  !
  eq_pts(:) = real(intseq(-mesh_order,mesh_order,2),kind=qp) / &
              real(max(mesh_order,1),kind=qp)
  !
  xi_qp(:) = real(xi(:),kind=qp)
  !
  ndim = size(rs ,dim=1)
  npts = size(xi)*size(rs,dim=2)
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
      return_value(:,k) = c * ((a-r-s)*p1 + (r+b)*p2 + (s+b)*p3)
      !
    end do
    !
  end do
  !
  ! Use the chop function to copy map_qp into return_value while also
  ! setting any values that are less that working machine epsilon to zero.
  !
  return_value = chop( map_qp )
  !
end function map_pris_to_solpts
!
!###############################################################################
!
pure function map_hexa_to_solpts(mp,np,xi,mesh_order,family) &
                          result(return_value)
  !
  !.. Use Statements ..
  use polynomial_mod, only : eval_LagrangePoly
  use polynomial_mod, only : eval_LagrangePoly3D
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: mp
  integer,  intent(in) :: np
  real(wp), intent(in) :: xi(:)
  integer,  intent(in) :: mesh_order
  !
  type(geom_family_t), intent(in) :: family
  !
  !.. Function Result ..
  real(wp) :: return_value(1:mp,1:np)
  !
  !.. Local Scalars ..
  integer  :: i,j,k,m,mi,mj,mk,n,ni,nj,nk,xpts,epts
  !
  !.. Local Arrays ..
  real(qp) :: xy_qp(1:3)
  real(qp) :: xi_qp(1:size(xi))
  real(qp) :: eq_pts(1:mesh_order+1)
  real(qp) :: map_qp(1:mp,1:np)
  !
  !.. Local Parameters ..
  real(qp), parameter :: end_pts(1:2) = [-qone,qone]
  !
continue
  !
  xpts = size(xi_qp)
  epts = size(eq_pts)
  !
  return_value(:,:) = zero
  !
  eq_pts(:) = real(intseq(-mesh_order,mesh_order,2),kind=qp) / &
              real(max(mesh_order,1),kind=qp)
  !
  xi_qp(:) = real(xi(:),kind=qp)
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
end function map_hexa_to_solpts
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
pure function get_solpts_edge(xi,xy) result(return_value)
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
pure function get_solpts_tria(rs,xy) result(return_value)
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
pure function get_solpts_quad(xi,xy) result(return_value)
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
pure function get_solpts_tetr(rs,xyz) result(return_value)
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
pure function get_solpts_pyra(xi,rs,xyz) result(return_value)
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
pure function get_solpts_pris(xi,rs,xyz) result(return_value)
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
pure function get_solpts_hexa(xi,xyz) result(return_value)
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
end module mappings_mod
