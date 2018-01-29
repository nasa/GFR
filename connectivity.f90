module connectivity_mod
  !
  !.. Global Use Statements ..
  use module_kind_types
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
  public :: get_grid_connectivity
  public :: find_cells_surrounding_each_node
  public :: find_cells_surrounding_each_cell
  public :: find_nodes_on_each_edge
  public :: find_cells_surrounding_each_edge
  public :: create_face_derived_type
  !
  public :: get_unit_cell_subconnectivity
  !
  public :: GetFaceNodes_for_QuadFace
  public :: GetFaceNodes_for_TriaFace
  public :: RotateFluxPoints_Quad
  public :: RotateFluxPoints_Tria
  !
  !
  ! #########################################
  ! #####   MODULE GENERIC INTERFACES   #####
  ! #########################################
  !
  interface RotateFluxPoints_Edge
    module procedure RotateFluxPoints_Edge, &
                     FluxPointRotation_Edge
  end interface RotateFluxPoints_Edge
  !
  interface RotateFluxPoints_Quad
    module procedure RotateFluxPoints_Quad, &
                     FluxPointRotation_Quad
  end interface RotateFluxPoints_Quad
  !
  interface RotateFluxPoints_Tria
    module procedure RotateFluxPoints_Tria, &
                     FluxPointRotation_Tria
  end interface RotateFluxPoints_Tria
  !
contains
!
!###############################################################################
!
subroutine get_grid_connectivity(bface,cell_geom,cell_order,xyz_nodes, &
                                 nodes_of_cell,nodes_of_cell_ptr)
  !
  !.. Use Statements ..
  use ovar,   only : continuous_output
  use geovar, only : nface,face,fp
  !
  !.. Formal Arguments ..
  integer,  dimension(:,:), intent(in) :: bface
  integer,  dimension(:),   intent(in) :: cell_geom
  integer,  dimension(:),   intent(in) :: cell_order
  integer,  dimension(:),   intent(in) :: nodes_of_cell
  integer,  dimension(:),   intent(in) :: nodes_of_cell_ptr
  real(wp), dimension(:,:), intent(in) :: xyz_nodes
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:)   :: cells_with_node
  integer, allocatable, dimension(:)   :: cells_with_node_ptr
  integer, allocatable, dimension(:)   :: cells_surr_cell
  integer, allocatable, dimension(:)   :: cells_surr_cell_ptr
  integer, allocatable, dimension(:)   :: cells_with_edge
  integer, allocatable, dimension(:)   :: cells_with_edge_ptr
  integer, allocatable, dimension(:)   :: edges_with_node
  integer, allocatable, dimension(:)   :: edges_with_node_ptr
  integer, allocatable, dimension(:,:) :: nodes_on_edge
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_grid_connectivity"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Get the cells surrounding each node
  !
  call memory_pause("entering grid connectivity master subroutine", &
                    "finding cells surrounding each node")
  call find_cells_surrounding_each_node(bface, &
                                        nodes_of_cell,nodes_of_cell_ptr, &
                                        cells_with_node,cells_with_node_ptr)
  !
  ! Get the cells surrounding each cell
  !
  call memory_pause("finding cells surrounding each node", &
                    "finding the cells surrounding each cell")
  call find_cells_surrounding_each_cell(bface,cell_geom, &
                                        nodes_of_cell,nodes_of_cell_ptr, &
                                        cells_with_node,cells_with_node_ptr, &
                                        cells_surr_cell,cells_surr_cell_ptr)
  !
  ! If using the continuous output option, find the
  ! node and cell connectivity for all the grid edges
  !
  if (continuous_output) then
    !
    ! Get the nodes defining each edge as well as the edges
    ! surrounding each node
    !
    call find_nodes_on_each_edge(cell_geom, &
                                 nodes_of_cell,nodes_of_cell_ptr, &
                                 cells_with_node,cells_with_node_ptr, &
                                 nodes_on_edge)
                                !nodes_on_edge, &
                                !edges_with_node,edges_with_node_ptr)
    !
    ! Get the cells that surround each edge
    !
    call find_cells_surrounding_each_edge(nodes_on_edge, &
                                          nodes_of_cell,nodes_of_cell_ptr, &
                                          cells_with_node,cells_with_node_ptr, &
                                          cells_with_edge,cells_with_edge_ptr)
    !
  end if
  !
  ! Get the cells on each face
  !
  call memory_pause("finding the cells surrounding each cell", &
                    "finding the cells on each face")
  call create_face_derived_type(bface,cell_geom,cell_order,xyz_nodes, &
                                nodes_of_cell,nodes_of_cell_ptr, &
                                cells_surr_cell,cells_surr_cell_ptr, &
                                face,fp)
  !
  ! Save the total number of faces to nface
  !
  nface = size(face)
  !
  ! Get the indices for the flux points for each cell
  !
  call memory_pause("finding the cells on each face", &
                    "finding the flux points of each cell")
  call create_cell_derived_type(nodes_of_cell,nodes_of_cell_ptr, &
                                cells_with_node_ptr,nodes_on_edge, &
                                cells_with_edge,cells_with_edge_ptr)
  !
  ! If using the continuous output option, create the edge derived type array
  !
  if (continuous_output) then
    call create_edge_derived_type(cells_with_edge,cells_with_edge_ptr, &
                                  nodes_on_edge)
  end if
  !
  ! Get the unit cell subconnectivities for use when writing solution file
  !
  call memory_pause("finding the flux points of each cell", &
                    "getting the unit cell subconnectivities")
  call get_unit_cell_subconnectivity
  call memory_pause("getting the unit cell subconnectivities", &
                    "leaving the grid connectivity master subroutine")
  !
  ! Deallocate the local arrays before leaving
  !
  if (allocated(cells_with_node)) then
    deallocate ( cells_with_node , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cells_with_node",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (allocated(cells_with_node_ptr)) then
    deallocate ( cells_with_node_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cells_with_node_ptr",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (allocated(cells_surr_cell)) then
    deallocate ( cells_surr_cell , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cells_surr_cell",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (allocated(cells_surr_cell_ptr)) then
    deallocate ( cells_surr_cell_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cells_surr_cell_ptr",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (allocated(edges_with_node)) then
    deallocate ( edges_with_node , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edges_with_node",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (allocated(edges_with_node_ptr)) then
    deallocate ( edges_with_node_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edges_with_node_ptr",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (allocated(cells_with_edge)) then
    deallocate ( cells_with_edge , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cells_with_edge",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (allocated(cells_with_edge_ptr)) then
    deallocate ( cells_with_edge_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cells_with_edge_ptr",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  if (allocated(nodes_on_edge)) then
    deallocate ( nodes_on_edge , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodes_on_edge",2,__LINE__,__FILE__,ierr, &
                     error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_grid_connectivity
!
!###############################################################################
!
subroutine find_cells_surrounding_each_node(bface, &
                                          nodes_of_cell,nodes_of_cell_ptr, &
                                          cells_with_node,cells_with_node_ptr, &
                                          include_bnd_faces, &
                                          make_bnd_faces_negative)
  !
  ! Find the cells that surround each node
  !
  ! For node i_node,
  ! cells_with_node(ibeg:iend) gives all of the grid cells that contain i_node
  ! where : ibeg = cells_with_node_ptr(i_node)+1
  !         iend = cells_with_node_ptr(i_node+1)
  !
  !.. Formal Arguments ..
  integer,              dimension(:,:),  intent(in) :: bface
  integer,              dimension(:),    intent(in) :: nodes_of_cell
  integer,              dimension(:),    intent(in) :: nodes_of_cell_ptr
  integer, allocatable, dimension(:), intent(inout) :: cells_with_node
  integer, allocatable, dimension(:), intent(inout) :: cells_with_node_ptr
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: include_bnd_faces
  logical(lk), optional, intent(in) :: make_bnd_faces_negative
  !
  !.. Local Scalars ..
  integer :: i,ip,is,ibeg,iend,n,nc,ierr
  integer :: n1,n2,np,lnode,lfbnd,lcell
  integer :: bnd_face_coef
  logical(lk) :: include_bnd_face_ghost_cells
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:) :: cwn
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_cells_surrounding_each_node"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Get the grid size from the input arrays
  !
  lnode = maxval(nodes_of_cell)
  lfbnd = size(bface,dim=2)
  lcell = size(nodes_of_cell_ptr)-1 - lfbnd
  !
  ! Determine whether to include the boundary face ghost cells
  ! when finding the cells surrounding each node
  !
  include_bnd_face_ghost_cells = true
  if (present(include_bnd_faces)) then
    include_bnd_face_ghost_cells = include_bnd_faces
  end if
  !
  ! Determine the coefficient for the boundary face ghost cell indices
  ! based on the optional argument make_bnd_faces_negative
  !
  bnd_face_coef = 1
  if (present(make_bnd_faces_negative)) then
    bnd_face_coef = merge(-1,1,make_bnd_faces_negative)
  end if
  !
  ! Allocate cells_with_node_ptr
  !
  allocate ( cells_with_node_ptr(1:lnode+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_with_node_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! Loop over the cells, storing 'ahead'
  !
  do n = 1,lcell
    !
    ibeg = nodes_of_cell_ptr(n)+1
    iend = nodes_of_cell_ptr(n+1)
    !
    do i = ibeg,iend
      ip = nodes_of_cell(i) + 1
      cells_with_node_ptr(ip) = cells_with_node_ptr(ip) + 1
    end do
    !
  end do
  !
  ! Add in the nodes of the ghost faces
  !
  if (include_bnd_face_ghost_cells) then
    !
    do n = 1,lfbnd
      !
      nc = bface(10,n)
      !
      ibeg = nodes_of_cell_ptr(nc)+1
      iend = nodes_of_cell_ptr(nc+1)
      !
      do i = ibeg,iend
        ip = nodes_of_cell(i)
        cells_with_node_ptr(ip+1) = cells_with_node_ptr(ip+1) + 1
      end do
      !
    end do
    !
  end if
  !
  ! Reshuffle cells_with_node_ptr
  !
  do n = 2,size(cells_with_node_ptr)
    cells_with_node_ptr(n) = cells_with_node_ptr(n) + cells_with_node_ptr(n-1)
  end do
  !
  allocate ( cells_with_node(1:last(cells_with_node_ptr)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_with_node",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! Now store the cells surrounding each node in cells_with_node
  !
  do n = 1,lcell
    !
    ibeg = nodes_of_cell_ptr(n)+1
    iend = nodes_of_cell_ptr(n+1)
    !
    do i = ibeg,iend
      ip = nodes_of_cell(i)
      is = cells_with_node_ptr(ip) + 1
      cells_with_node_ptr(ip) = is
      cells_with_node(is) = n
    end do
    !
  end do
  !
  ! Add in the nodes of the ghost faces
  !
  if (include_bnd_face_ghost_cells) then
    !
    do n = 1,lfbnd
      !
      nc = bface(10,n)
      !
      ibeg = nodes_of_cell_ptr(nc)+1
      iend = nodes_of_cell_ptr(nc+1)
      !
      do i = ibeg,iend
        ip = nodes_of_cell(i)
        is = cells_with_node_ptr(ip) + 1
        cells_with_node_ptr(ip) = is
        !
        ! Make the index for this boundary face ghost cell negative
        ! if the optional argument make_bnd_faces_negative is present
        ! with the value .TRUE.
        !
        cells_with_node(is) = bnd_face_coef * nc
        !
      end do
      !
    end do
    !
  end if
  !
  ! Finally, reorder cells_with_node_ptr
  !
  do n = size(cells_with_node_ptr),2,-1
    cells_with_node_ptr(n) = cells_with_node_ptr(n-1)
  end do
  cells_with_node_ptr(1) = 0
  !
  ! Perform one last error check to make sure that no cells have been added
  ! more than once to a given node. This was known to happen when the ghost
  ! cells were treated as actual cells but shouldnt happen anymore now that
  ! ghost cells are actually ghost faces. I have removed the old code that
  ! prevented this from happening in the above loops since it required a mask
  ! array that had the potential to consume a lot of memory. This error check
  ! has been added JUST-IN-CASE I have overlooked a scenario in which this
  ! could possible happen.
  !
  nc = 0
  do n = 1,lnode
    nc = max( nc , cells_with_node_ptr(n+1) - cells_with_node_ptr(n) )
  end do
  !
  allocate ( cwn(1:nc) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cwn",1,__LINE__,__FILE__,ierr,error_message, &
                   skip_alloc_pause)
  !
  do n = 1,lnode
    !
    n1 = cells_with_node_ptr(n)+1
    n2 = cells_with_node_ptr(n+1)
    np = n2-n1+1
    !
    cwn(1:np) = cells_with_node(n1:n2)
    call make_unique(cwn(1:np),ip)
    if (ip /= np) then
      !
      write (iout,9023) n,(cells_with_node(i),i=n1,n2)
      9023 format ("Found a duplicate cell around node ",i0,/, &
                   "  Cell list: ",100(2x,i0))
                  !"  Cell list: ",*(2x,i0))
      !
      ! Compute the difference from the original number of cells found for
      ! this node and the new unique number of cells that surround this node
      nc = np - ip
      ! Adjust the pointer index for the remaining nodes
      cells_with_node_ptr(n+1:) = cells_with_node_ptr(n+1:) - nc
      ! Get the new value for n2
      n2 = n1 + ip - 1
      ! Copy the unique cell list back into cells_with_node
      cells_with_node(n1:n2) = cwn(1:ip)
      ! Perform end-off-shifts to move the remaining parts of cells_with_node
      ! to the end of the list of cells for the current node
      cells_with_node(n2+1:) = eoshift(cells_with_node(n2+1:),shift=nc)
      !
    end if
    !
  end do
  !
  if (allocated(cwn)) then
    deallocate ( cwn , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cwn",2,__LINE__,__FILE__,ierr,error_message, &
                     skip_alloc_pause)
  end if
  !
  ! Reallocate cells_with_node if it was changed
  !
  if (last(cells_with_node_ptr) /= size(cells_with_node)) then
    !
    n1 = last(cells_with_node_ptr) + 1
    n2 = ubound(cells_with_node,dim=1)
    !
    if (any(cells_with_node(n1:n2) /= 0)) then
      write (error_message,10)
      call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
    else
      call reallocate(cells_with_node,last(cells_with_node_ptr))
    end if
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  10 format ("Error: The end of cells_with_node should be zero!")
  !
contains
  !
  pure subroutine make_unique(cell_list,nc)
    !.. Formal Arguments ..
    integer,                 intent(out) :: nc
    integer, dimension(:), intent(inout) :: cell_list
    !.. Local Scalars ..
    integer :: i
    !.. Local Arrays ..
    integer, dimension(1:size(cell_list)) :: array_of_zeros
  continue
    !
    array_of_zeros = 0
    !
    nc = size(cell_list)
    do i = 1,nc-1
      where (cell_list(i+1:nc) == cell_list(i)) cell_list(i+1:nc) = 0
    end do
    !
    nc = count( cell_list /= 0 )
    cell_list = pack( cell_list , cell_list /= 0 , array_of_zeros )
    !
  end subroutine make_unique
  !
end subroutine find_cells_surrounding_each_node
!
!###############################################################################
!
subroutine find_cells_surrounding_each_cell(bface,cell_geom, &
                                          nodes_of_cell,nodes_of_cell_ptr, &
                                          cells_with_node,cells_with_node_ptr, &
                                          cells_surr_cell,cells_surr_cell_ptr)
  !
  ! Find the cells that surround each cell
  !
  ! For cell i_cell, ...
  ! cells_surr_cell(ibeg:iend) gives all of the grid cells
  !                            that are adjacent to i_cell
  ! where : ibeg = cells_surr_cell_ptr(i_cell)+1
  !         iend = cells_surr_cell_ptr(i_cell+1)
  !
  !.. Formal Arguments ..
  integer,              dimension(:,:),  intent(in) :: bface
  integer,              dimension(:),    intent(in) :: cell_geom
  integer,              dimension(:),    intent(in) :: nodes_of_cell
  integer,              dimension(:),    intent(in) :: nodes_of_cell_ptr
  integer,              dimension(:),    intent(in) :: cells_with_node
  integer,              dimension(:),    intent(in) :: cells_with_node_ptr
  integer, allocatable, dimension(:), intent(inout) :: cells_surr_cell
  integer, allocatable, dimension(:), intent(inout) :: cells_surr_cell_ptr
  !
  !.. Local Scalars ..
  integer :: i,i1,i2,ieadj,csc_ptr
  integer :: j,j1,j2,n,n1,n2,ierr
  integer :: this_geom,this_side,nfnods
  integer :: lnode,lcell
  !
  !.. Local Arrays ..
  integer, dimension(1:4) :: ip
  integer, dimension(1:4) :: face_nodes
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:) :: lpoin
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_cells_surrounding_each_cell"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  lnode = maxval(nodes_of_cell)
  lcell = size(nodes_of_cell_ptr)-1 - size(bface,dim=2)
  !
  allocate ( lpoin(1:lnode) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"lpoin",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cells_surr_cell_ptr(1:lcell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_surr_cell_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! Loop to set up cells_surr_cell_ptr
  !
  do n = 1,lcell
    cells_surr_cell_ptr(n+1) = cells_surr_cell_ptr(n) + geom_faces(cell_geom(n))
  end do
  !
  allocate ( cells_surr_cell(1:last(cells_surr_cell_ptr)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_surr_cell",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! Loop to find cells_surr_cell
  !
  do n = 1,lcell
    !
    ! Get the surrounding cells based on the geometry of the current cell
    !
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    !
    ! Get the pointer index of cells_surr_cell for this cell
    !
    csc_ptr = cells_surr_cell_ptr(n)
    !
    ! Geometry of the current cell
    !
    this_geom = cell_geom(n)
    !
    ! Loop over the number of cell faces for the current geometry
    !
    do this_side = 1,geom_faces(this_geom)
      !
      ! Find the number of nodes that need to be matched for this face
      !
      nfnods = num_face_nodes(this_side,this_geom)
      !
      ! Node offsets for this face of the current cell
      !
      ip(1:nfnods) = side_nodes(1:nfnods,this_side,this_geom)
      !
      ! Store nodes of the current face for the current cell
      !
      face_nodes(1:nfnods) = nodes_of_cell( n1 + ip(1:nfnods) )
      !
      ! Mark each of these nodes on the current face
      !
      lpoin( face_nodes(1:nfnods) ) = 1
      !
      ! Get the beginning and ending indices for the
      ! cells surrounding the first node on the face
      !
      i1 = cells_with_node_ptr( face_nodes(1) )+1
      i2 = cells_with_node_ptr( face_nodes(1) +1)
      !
      ! Initialize the adjacent face marker to zero
      !
      ieadj = 0
      !
      ! Loop over the cells that contain the first node on the face
      !
      find_cell: do i = i1,i2
        !
        ! Current cell containing this node
        ! NOTE: We need the absolute value of this because the cell index could
        !       be negative if it is for a boundary face ghost cell and the
        !       optional argument make_bnd_faces_negative was present with the
        !       value .TRUE. when subroutine find_cells_surrounding_each_node
        !       was called
        !
        j = abs(cells_with_node(i))
        !
        ! Skip to the next cell if this is the current cell we are working on
        !
        if (j == n) cycle find_cell
        !
        ! Get the beginning and ending indices of the nodes defining the
        ! possible adjacent cell
        !
        j1 = nodes_of_cell_ptr(j)+1
        j2 = nodes_of_cell_ptr(j+1)
        !
        ! If summing the lpoin array using the current nodes of the possible
        ! adjacent cell is equal to the number of nodes we need to match, then
        ! this is the adjacent cell to the current face.
        ! NOTE: We need to use the actual value of cells_with_node(i) instead
        !       of j because we want to keep the cell index negative if this
        !       is a boundary face ghost cell
        !
        if (sum(lpoin(nodes_of_cell(j1:j2))) == nfnods) then
          ieadj = cells_with_node(i)
          exit find_cell
        end if
        !
      end do find_cell
      !
      ! Store the adjacent cell index
      !
      cells_surr_cell(csc_ptr+this_side) = ieadj
      !
      ! Re-mark the nodes in lpoin to zero
      !
      lpoin( face_nodes(1:nfnods) ) = 0
      !
    end do
    !
  end do
  !
  deallocate ( lpoin , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"lpoin",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine find_cells_surrounding_each_cell
!
!###############################################################################
!
subroutine initialize_flxpts(fp)
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  use order_mod, only : n_min_order,n_max_order
  use geovar,    only : fp_t
  !
  !.. Formal Arguments ..
  type(fp_t), allocatable, intent(inout) :: fp
  !
  !.. Local Scalars ..
  integer :: gmin,gmax,omin,omax,ierr
  integer :: n,nep,nfp,this_order
  !
  character(len=100) :: array_name
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: fp_idx(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_flxpts"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  gmin = min(Geom_Edge,Geom_Tria,Geom_Quad)
  gmax = max(Geom_Edge,Geom_Tria,Geom_Quad)
  !
  omin = n_min_order
  omax = n_max_order
  !
  allocate ( fp , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fp",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( fp%geord(gmin:gmax,omin:omax) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"fp%geord",1,__LINE__,__FILE__,ierr,error_message)
  !
  do this_order = omin,omax
    !
    ! Edge face
    !
    nep = geom_solpts(Geom_Edge,this_order)
    !
    allocate ( fp%geord(Geom_Edge,this_order)%rot(0:2) , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "Geom_Edge",this_order
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( fp%geord(Geom_Edge,this_order)%rot(0)%pts(1:nep) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "Geom_Edge",this_order,1
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( fp%geord(Geom_Edge,this_order)%rot(1)%pts(1:nep) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "Geom_Edge",this_order,1
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( fp%geord(Geom_Edge,this_order)%rot(2)%pts(1:nep) , source=0 , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "Geom_Edge",this_order,2
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    fp%geord(Geom_Edge,this_order)%rot(0)%pts(1:nep) = 0
    fp%geord(Geom_Edge,this_order)%rot(1)%pts(1:nep) = intseq(0,nep-1)
    fp%geord(Geom_Edge,this_order)%rot(2)%pts(1:nep) = intseq(nep-1,0,-1)
    !
    ! Triangular Face
    !
    nfp = geom_solpts(Geom_Tria,this_order)
    !
    allocate ( fp%geord(Geom_Tria,this_order)%rot(0:6) , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "Geom_Tria",this_order
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    do n = lbound(fp%geord(Geom_Tria,this_order)%rot,dim=1), &
           ubound(fp%geord(Geom_Tria,this_order)%rot,dim=1)
      !
      allocate ( fp%geord(Geom_Tria,this_order)%rot(n)%pts(1:nfp) , source=0 , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) "Geom_Tria",this_order,n
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
    end do
    !
    ! Quadrilateral Face
    !
    nfp = geom_solpts(Geom_Quad,this_order)
    !
    allocate ( fp%geord(Geom_Quad,this_order)%rot(0:8) , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "Geom_Quad",this_order
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    do n = lbound(fp%geord(Geom_Quad,this_order)%rot,dim=1), &
           ubound(fp%geord(Geom_Quad,this_order)%rot,dim=1)
      !
      allocate ( fp%geord(Geom_Quad,this_order)%rot(n)%pts(1:nfp) , source=0 , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) "Geom_Quad",this_order,n
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
    end do
    !
    ! See the function RotateFluxPoints_Quad for a complete description of
    ! how these rotations are defined
    !
    if (allocated(fp_idx)) then
      deallocate ( fp_idx , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"fp_idx",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    ! Create fp_idx as a square 2D array containing
    ! the consecutive indices from 0 to nfp-1
    !
    allocate ( fp_idx(1:nep,1:nep) , &
               source=reshape( intseq(0,nfp-1) , [nep,nep] ) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"fp_idx",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Set all the flux point indices for the error rotation to 0
    !
    fp%geord(Geom_Quad,this_order)%rot(0)%pts(1:nfp) = 0
    !
    ! Create the flux point indices for the first four rotations where the
    ! "a" and "b" vectors are aligned in the same coordinate directions
    !
    fp%geord(Geom_Quad,this_order)%rot(1)%pts(1:nfp) = reshape( fp_idx , [nfp] )
    !
    fp%geord(Geom_Quad,this_order)%rot(2)%pts(1:nfp) = &
                                         reshape( fp_idx(:,nep:1:-1) , [nfp] )
    !
    fp%geord(Geom_Quad,this_order)%rot(3)%pts(1:nfp) = &
                                         reshape( fp_idx(nep:1:-1,:) , [nfp] )
    !
    fp%geord(Geom_Quad,this_order)%rot(4)%pts(1:nfp) = &
                                  reshape( fp_idx(nep:1:-1,nep:1:-1) , [nfp] )
    !
    ! Create the flux point indices for the last four rotations where the
    ! "a" and "b" vectors are aligned in opposite coordinate directions,
    ! which requires that we first transpose fp_idx
    !
    fp_idx = transpose(fp_idx)
    !
    fp%geord(Geom_Quad,this_order)%rot(5)%pts(1:nfp) = reshape( fp_idx , [nfp] )
    !
    fp%geord(Geom_Quad,this_order)%rot(6)%pts(1:nfp) = &
                                         reshape( fp_idx(:,nep:1:-1) , [nfp] )
    !
    fp%geord(Geom_Quad,this_order)%rot(7)%pts(1:nfp) = &
                                         reshape( fp_idx(nep:1:-1,:) , [nfp] )
    !
    fp%geord(Geom_Quad,this_order)%rot(8)%pts(1:nfp) = &
                                  reshape( fp_idx(nep:1:-1,nep:1:-1) , [nfp] )
    !
  end do
  !
  if (allocated(fp_idx)) then
    deallocate ( fp_idx , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"fp_idx",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("fp%geord(",a,",",i0,")%rot",:,"%(",i0,")%pts")
  !
end subroutine initialize_flxpts
!
!###############################################################################
!
subroutine create_face_derived_type(bface,cell_geom,cell_order,xyz_nodes, &
                                    nodes_of_cell,nodes_of_cell_ptr, &
                                    cells_surr_cell,cells_surr_cell_ptr, &
                                    face,fp,force_flxpt_offset_to_zero)
  !
  ! Procedure that creates the face array.
  ! Look in the module geovar for the definition of the derived type.
  !
  !.. Use Statements ..
  use order_mod, only : maxFP
  use order_mod, only : geom_solpts
  use order_mod, only : find_geom_max_interp_order
  use geovar,    only : face_t,fp_t
  !
  !.. Formal Arguments ..
  integer,                      intent(in) :: bface(:,:)
  integer,                      intent(in) :: cell_geom(:)
  integer,                      intent(in) :: cell_order(:)
  integer,                      intent(in) :: nodes_of_cell(:)
  integer,                      intent(in) :: nodes_of_cell_ptr(:)
  real(wp),                     intent(in) :: xyz_nodes(:,:)
  integer,                      intent(in) :: cells_surr_cell(:)
  integer,                      intent(in) :: cells_surr_cell_ptr(:)
  type(face_t), allocatable, intent(inout) :: face(:)
  type(fp_t),   allocatable, intent(inout) :: fp
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: force_flxpt_offset_to_zero
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,nf,nc,ierr
  integer :: i,ic,ifa,ip1,ip2,ip3,jc
  integer :: ibeg,iend,interval
  integer :: gmin,gmax,omin,omax
  integer :: face_counter,nfp,l1,r1,l2,r2
  integer :: left_cell,left_geom
  integer :: host_cell,ghost_cell,host_geom
  integer :: host_side,nfnods,ncnods
  integer :: right_cell,right_geom,ncf
  integer :: left_side,right_side
  integer :: face_geom,face_order,face_pts
  integer :: c1,c2,this_cell,this_geom,adj_cell
  integer :: left_order,right_order
  integer :: rotation
  character(len=200) :: array_name
  logical(lk) :: error_finding_rotation
  logical(lk) :: flxpt_offset_is_zero
  !
  integer :: nr,lnode,lfbnd,lface,lcell
  !
  real(wp) :: dir
  !
  !.. Local Arrays ..
  integer,  dimension(1:4) :: ip
  integer,  dimension(1:4) :: left_nodes
  integer,  dimension(1:4) :: right_nodes
 !integer,  dimension(1:4) :: face_nodes
  real(wp), dimension(1:3) :: rnrm
  real(wp), dimension(1:3) :: face_center
  real(wp), dimension(1:3) :: cell_center
  real(wp), dimension(1:3,1:4) :: xyz_face
  real(wp), dimension(1:3,1:8) :: xyz_cell
 !integer,  dimension(1:maxFP) :: rot_fpts
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:) :: lpoin
  integer, allocatable, dimension(:,:) :: faces_of_cell
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_face_derived_type"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Get the grid dimensions from the input arrays
  !
  nr = size(xyz_nodes,dim=1)
  lnode = size(xyz_nodes,dim=2)
  lfbnd = size(bface,dim=2)
  lcell = size(nodes_of_cell_ptr)-1 - lfbnd
  !
  ! Initialize the flxpts array
  !
  call initialize_flxpts(fp)
  !
  ! ###########
  ! ###########
  ! ###  1  ###
  ! ###########
  ! ###########
  call debug_timer(start_timer,"section 1")
  !
  ! Find the total number of faces
  ! and then allocate the face array
  !
  lface = lfbnd
  do n = 1,lcell
    !
    n1 = cells_surr_cell_ptr(n)+1
    n2 = cells_surr_cell_ptr(n+1)
    !
    lface = lface + count( cells_surr_cell(n1:n2)>n .and. &
                           cells_surr_cell(n1:n2)<=lcell )
    !
  end do
  !
  allocate ( face(1:lface) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( faces_of_cell(1:6,1:lcell) , source=-99 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"faces_of_cell",1,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(stop_timer,"section 1")
  ! ###########
  ! ###########
  ! ###  2  ###
  ! ###########
  ! ###########
  call debug_timer(start_timer,"section 2")
  !
  ! Get the cells on each boundary face
  !
  do nf = 1,lfbnd
    !
    host_cell = bface(2,nf)   ! Host cell on boundary face
    ghost_cell = bface(10,nf) ! Ghost cell on boundary face
    !
    ! Store the cells on this boundary face
    !
    face(nf)%left%cell  = host_cell
    face(nf)%right%cell = ghost_cell
    !
    host_geom = cell_geom(host_cell)             ! geometry of host cell
    host_side = bface(4,nf)                      ! face of host cell on boundary
    nfnods = num_face_nodes(host_side,host_geom) ! # of nodes defining bnd face
    !
    ! beginning index of nodes defining host cell
    n1 = nodes_of_cell_ptr(host_cell)+1
    ! ending    index of nodes defining host cell
    n2 = nodes_of_cell_ptr(host_cell+1)
    !
    ! Node offsets for this face
    ip(1:nfnods) = side_nodes(1:nfnods,host_side,host_geom)
    !
    ! Store the nodes on this boundary face
    !
    face(nf)%nodes(1:nfnods) = nodes_of_cell( n1 + ip(1:nfnods) )
    !
    ! Store the geometry of this face
    !
    face(nf)%geom = geom_of_face(nfnods)
    !
    ! Store the side of the left cell that is on this face since we have this
    ! information available right now. Store the side of the right/ghost cell
    ! as well since this will always be 1.
    !
    face(nf)%left%cell_face  = host_side
    face(nf)%right%cell_face = 1
    !
    ! Store this face for the host cell
    !
    faces_of_cell(host_side,host_cell) = nf
    !
  end do
  !
  call debug_timer(stop_timer,"section 2")
  ! ###########
  ! ###########
  ! ###  3  ###
  ! ###########
  ! ###########
  call debug_timer(start_timer,"section 3")
  !
  ! Get the cells for the interior faces
  !
  face_counter = lfbnd
  !
  do this_cell = 1,lcell
    !
    c1 = cells_surr_cell_ptr(this_cell)+1
    c2 = cells_surr_cell_ptr(this_cell+1)
    !
    n1 = nodes_of_cell_ptr(this_cell)+1
    n2 = nodes_of_cell_ptr(this_cell+1)
    !
    this_geom = cell_geom(this_cell)
    !
    do nc = c1,c2
      !
      adj_cell = cells_surr_cell(nc)
      !
      if (adj_cell>this_cell .and. adj_cell<=lcell) then
        !
        face_counter = face_counter + 1
        !
        nf = nc-c1+1
        !
        ! Store the cells on this interior face
        !
        face(face_counter)%left%cell  = this_cell
        face(face_counter)%right%cell = adj_cell
        !
        ! # of nodes defining this face
        nfnods = num_face_nodes(nf,this_geom)
        ! node offsets for this face
        ip(1:nfnods) = side_nodes(1:nfnods,nf,this_geom)
        !
        ! Store the nodes on this interior face
        !
        face(face_counter)%nodes(1:nfnods) = nodes_of_cell( n1 + ip(1:nfnods) )
        !
        ! Store the geometry of this face
        !
        face(face_counter)%geom = geom_of_face(nfnods)
        !
        ! Store the side of the left cell that is on this face since we have
        ! this information available right now. We will figure out the side of
        ! the right cell in section 5.
        !
        face(face_counter)%left%cell_face = nf
        !
        ! Store this face for the left cell
        !
        faces_of_cell(nf,this_cell) = face_counter
        !
      end if
      !
    end do
    !
  end do
  !
  call debug_timer(stop_timer,"section 3")
  ! ###########
  ! ###########
  ! ###  4  ###
  ! ###########
  ! ###########
  call debug_timer(start_timer,"section 4")
  !
  ! This is a test to make sure that face(nf)%left%cell is on the left side of
  ! this face and face(nf)%right%cell is on the right side
  !
  ic = 0
  do nf = lfbnd+1,lface
    !
    left_cell = face(nf)%left%cell ! left cell on this face
    face_geom = face(nf)%geom      ! geometry of this face
    left_geom = cell_geom(left_cell) ! geometry of the left cell
    !
    ! # of nodes defining the face and defining the left cell
    !
    nfnods = geom_nodes(face_geom)
    ncnods = geom_nodes(left_geom)
    !
    ! Coordinates of the nodes on this face
    !
    xyz_face(1:nr,1:nfnods) = xyz_nodes(1:nr,face(nf)%nodes(1:nfnods))
   !face_nodes(1:nfnods) = face(nf)%nodes(1:nfnods)
   !do n = 1,nfnods
   !  do i = 1,nr
   !    xyz_face(i,n) = xyz_nodes(i,face_nodes(n))
   !  end do
   !end do
    !
    ! Coordinates of the center of this face
    !
   !do i = 1,nr
   !  face_center(i) = zero
   !  do n = 1,nfnods
   !    face_center(i) = face_center(i) + xyz_face(i,n)
   !  end do
   !  face_center(i) = face_center(i) / real(nfnods,kind=wp)
   !end do
    face_center(1:nr) = sum(xyz_face(1:nr,1:nfnods),dim=2) &
                      / real(nfnods,kind=wp)
    !
    ! Get the coordinates of the center of the left cell on this face
    !
    n1 = nodes_of_cell_ptr(left_cell)+1 ! index of first node of left cell
    n2 = nodes_of_cell_ptr(left_cell+1) ! index of last  node of left cell
    !
    ! Coordinates of nodes defining left cell
    !
    xyz_cell(1:nr,1:ncnods) = xyz_nodes(1:nr,nodes_of_cell(n1:n2))
    !
    ! Coordinates of the center of the left cell
    !
    cell_center(1:nr) = sum(xyz_cell(1:nr,1:ncnods),dim=2) &
                      / real(ncnods,kind=wp)
    !
    ! Get the normal to this face using the nodes that define it
    !
    rnrm(1:nr) = face_normal( xyz_face(1:nr,1:nfnods) )
    !
    ! Compute the dot product of the vector from the center of
    ! the left cell to the center of the face with the face normal
    !
    dir = dot( rnrm(1:nr) , face_center(1:nr) - cell_center(1:nr) )
    !
    ! If dot is negative, then the cell stored in face(nf)%right%cell
    ! is actually the cell to the left of the face vector and needs
    ! to be stored in face(nf)%left%cell
    !
    if (dir < zero) then
      write (*,1) nf,(face(nf)%nodes(i),i=1,nfnods)
      ic = 1
    end if
    !
  end do
  if (ic == 1) call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  !
  call debug_timer(stop_timer,"section 4")
  ! ###########
  ! ###########
  ! ###  5  ###
  ! ###########
  ! ###########
  call debug_timer(start_timer,"section 5")
  !
  ! Allocate the array for marking the face nodes of a cell
  !
  allocate ( lpoin(1:lnode) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"lpoin",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Find the side of the right cell that is on each interior face
  !
  loop_right_cell_face: do nf = lfbnd+1,lface
    !
    nfnods = geom_nodes(face(nf)%geom) ! # of nodes defining this face
    !
    ! Mark the nodes on the current face
    !
    lpoin( face(nf)%nodes(1:nfnods) ) = 1
    !
    ! Information for the right cell on the face
    !
    right_cell = face(nf)%right%cell
    right_geom = cell_geom(right_cell)
    !
    n1 = nodes_of_cell_ptr(right_cell)+1
    n2 = nodes_of_cell_ptr(right_cell+1)
    !
    loop_cell_faces: do ncf = 1,geom_faces(right_geom)
      !
      ! Skip this cell face if the number of face nodes doesnt match
      !
      if (num_face_nodes(ncf,right_geom) /= nfnods) cycle loop_cell_faces
      !
      ! Store the nodes of the cell face in question
      !
      ip(1:nfnods) = nodes_of_cell(n1+side_nodes(1:nfnods,ncf,right_geom))
      !
      ! If this face subset of lpoin sums up to nfnods, this is the side
      ! of the right cell on the current face. Store the value, reset all
      ! of lpoin to zero, and skip to the next face.
      !
      if (sum(lpoin(ip(1:nfnods))) == nfnods) then
        face(nf)%right%cell_face = ncf
        faces_of_cell(ncf,right_cell) = nf
        lpoin( ip(1:nfnods) ) = 0
        cycle loop_right_cell_face
      end if
      !
    end do loop_cell_faces
    !
  end do loop_right_cell_face
  !
  deallocate ( lpoin , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"lpoin",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(stop_timer,"section 5")
  ! ###########
  ! ###########
  ! ###  6  ###
  ! ###########
  ! ###########
  call debug_timer(start_timer,"section 6")
  !
  ! Get the order for each face based on the max order of each cell on the face
  !
  do nf = 1,lface
    !
    left_cell = face(nf)%left%cell
    right_cell = face(nf)%right%cell
    !
    left_order = cell_order(left_cell)
    right_order = cell_order(right_cell)
    !
    face(nf)%order = max(left_order,right_order)
    !
  end do
  !
  call debug_timer(stop_timer,"section 6")
  ! ###########
  ! ###########
  ! ###  7  ###
  ! ###########
  ! ###########
  call debug_timer(start_timer,"section 7")
  !
  ! Before we start finding the face/flux points of the cells on each face,
  ! check the whether the optional argument force_flxpt_offset_to_zero is
  ! present and assign its value to flxpt_offset_is_zero.
  !
  ! Initialize flxpt_offset_is_zero to .FALSE.
  !
  flxpt_offset_is_zero = fals
  !
  if (present(force_flxpt_offset_to_zero)) then
    !
    ! If this is called from the parallel module, the flux point indices are
    ! used for communicating the packed up face data stored in the faceusp and
    ! facedusp arrays. These indices are needed for the displacement length
    ! offsets that are required to create MPI datatypes to communication this
    ! data. Since MPI requires these displacement length offsets to zero, this
    ! optional flag is used to make sure the indices of the flux point start at
    ! 0 instead of 1 as is the default. Additionally, the MPI offsets are for
    ! each face so the flux point indices will always be in the range of
    ! [0:nfp-1]. For the left cell, the indices will be exactly [0:nfp-1], but
    ! the indices for the right cell will depend on the side of the right cell
    ! on the interface and the rotation of the right cell relative to the left
    ! cell on the interface.
    !
    flxpt_offset_is_zero = force_flxpt_offset_to_zero
    !
  end if
  !
  ! For each cell on the face, get the indices for the face/flux points of the
  ! cell that are on this face. These indices are used in the correct and
  ! interp(:)%toFace matrices.
  !
  error_finding_rotation = fals
  !
  do nf = 1,lface
    !
    left_cell = face(nf)%left%cell
    right_cell = face(nf)%right%cell
    !
    left_side = face(nf)%left%cell_face
    right_side = face(nf)%right%cell_face
    !
    left_geom = cell_geom(left_cell)
    right_geom = cell_geom(right_cell)
    !
    left_order = cell_order(left_cell)
    right_order = cell_order(right_cell)
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    nfp = geom_solpts(face_geom,face_order)
    !
    nfnods = geom_nodes(face_geom)
    !
    l1 = nodes_of_cell_ptr(left_cell)+1
    l2 = nodes_of_cell_ptr(left_cell+1)
    !
    r1 = nodes_of_cell_ptr(right_cell)+1
    r2 = nodes_of_cell_ptr(right_cell+1)
    !
    ip(1:nfnods) = side_nodes(1:nfnods,left_side,left_geom)
    left_nodes(1:nfnods) = nodes_of_cell(l1 + ip(1:nfnods))
    !
    ! #############################
    ! #############################
    ! #############################
    ! ###                       ###
    ! ### LEFT CELL ON THE FACE ###
    ! ###                       ###
    ! #############################
    ! #############################
    ! #############################
    !
    ! Find the offset for the flux points of the left cell
    !
    if (flxpt_offset_is_zero) then
      !
      ! Set the offset for the flux points to zero for the left cell
      ! if the optional argument force_flxpt_offset_to_zero was present
      ! with the value .TRUE.
      !
      face_pts = 0
      !
    else
      !
      ! Initialize the flux point offset to 1
      !
      face_pts = 1
      !
      ! For the left cell, the starting index is the sum of all face points
      ! for all faces up to the current face + 1. The '+1' is why face_pts
      ! is initialized to 1. We need to use face_order when computing the
      ! number of flux points for each cell face before the current one.
      ! This is the index needed because this is used to get the correct
      ! location in the interp%toFace(face_order)%mat and correct%mat arrays.
      ! These arrays are set up with the indexing within mat assuming that
      ! the order of all faces is face_order.
      !
      do ncf = 1,left_side-1
        !
        ifa = faces_of_cell(ncf,left_cell)
        !
        face_pts = face_pts + geom_solpts(face(ifa)%geom,face_order)
        !
      end do
      !
    end if
    !
    ! Get the rotation index for this face
    ! NOTE: The rotation index is always 1 for the left cell on a face
    !
   !nc = 1
   !rotation = nc
    rotation = 1
    !
    ! Save the offset for the flux points
    !
    face(nf)%left%fp_offset = face_pts
    !
    ! Save the rotation index for this face
    !
    face(nf)%left%rotation = rotation
    !
    ! ##############################
    ! ##############################
    ! ##############################
    ! ###                        ###
    ! ### RIGHT CELL ON THE FACE ###
    ! ###                        ###
    ! ##############################
    ! ##############################
    ! ##############################
    !
    ! For the right cell, it depends on whether
    ! this is a boundary or interior face
    !
    ! Find the offset for the flux points of the left cell
    !
    if (flxpt_offset_is_zero) then
      !
      ! Set the offset for the flux points to zero for the left cell
      ! if the optional argument force_flxpt_offset_to_zero was present
      ! with the value .TRUE.
      !
      face_pts = 0
      !
    else
      !
      ! Initialize the flux point offset to 1
      !
      face_pts = 1
      !
      ! The offset is just 1 for a boundary face, so we only need to continue
      ! with computing the flux point offset if this is an interior face
      !
      if (nf > lfbnd) then
        !
        ! For the right cell, the starting index is the sum of all face points
        ! for all faces up to the current face + 1. The '+1' is why face_pts
        ! is initialized to 1. We need to use face_order when computing the
        ! number of flux points for each cell face before the current one.
        ! This is the index needed because this is used to get the correct
        ! location in the interp%toFace(face_order)%mat and correct%mat arrays.
        ! These arrays are set up with the indexing within mat assuming that
        ! the order of all faces is face_order.
        !
        do ncf = 1,right_side-1
          !
          ifa = faces_of_cell(ncf,right_cell)
          !
          face_pts = face_pts + geom_solpts(face(ifa)%geom,face_order)
          !
        end do
        !
      end if
      !
    end if
    !
    ! Get the local cell indices for the nodes of the right cell that are
    ! located on the current face
    !
    if (nf <= lfbnd) then
      !
      ! This is a boundary face so the nodes on this face are simply the nodes
      ! defining this face, i.e. nodes_of_cell(r1:r1-1+nfnods)
      !
      ip(1:nfnods) = intseq(0,nfnods-1)
      !
    else
      !
      ! This is an interior face so the indices for the nodes of the right cell
      ! located on this face depend on the geometry of the right cell and the
      ! side of the right cell thats on this face
      !
      ip(1:nfnods) = side_nodes(1:nfnods,right_side,right_geom)
      !
    end if
    !
    ! Use the node indices to collect the nodes of the right cell on this face
    !
    right_nodes(1:nfnods) = nodes_of_cell(r1 + ip(1:nfnods))
    !
    ! Get the flux point rotation for the right cell
    !
    if (face_geom == Geom_Edge) then
      !
      ! If this is a 1D edge between two 2D cells, the ordering of the flux
      ! points for the right cell is just the reverse ordering of the flux
      ! points for the left cell.
      !
      rotation = RotateFluxPoints_Edge(left_nodes(1:nfnods), &
                                       right_nodes(1:nfnods))
      !
      ! If the previous function could not find the rotation, mark that there
      ! was an error finding the correct face rotation between cells.
      !
      if (rotation == 0) error_finding_rotation = true
      !
      face(nf)%right%fp_offset = face_pts
      face(nf)%right%rotation = rotation
      !
    else if (face_geom == Geom_Quad) then
      !
      ! Get the nodes from the left and right cells that are on the current face
      ! and in the correct ABCD ordering needed for the RotateFluxPoints_Quad
      ! function.
      !
      call GetFaceNodes_for_QuadFace(nodes_of_cell(l1:l2),left_side, &
                                     nodes_of_cell(r1:r2),right_side, &
                                     left_nodes(1:nfnods),right_nodes(1:nfnods))
      !
      ! Find the face rotation between the two cells on the interface and use
      ! this rotation to get the correct ordering of the flux points for the
      ! right cell relative to the flux points of the left cell.
      !
      rotation = RotateFluxPoints_Quad(left_nodes(1:nfnods), &
                                       right_nodes(1:nfnods))
      !
      ! If the previous function could not find the rotation, mark that there
      ! was an error finding the correct face rotation between cells.
      !
      if (rotation == 0) error_finding_rotation = true
      !
      face(nf)%right%fp_offset = face_pts
      face(nf)%right%rotation = rotation
      !
    else if (face_geom == Geom_Tria) then
      !
      call GetFaceNodes_for_TriaFace(nodes_of_cell(l1:l2),left_side, &
                                     nodes_of_cell(r1:r2),right_side, &
                                     left_nodes(1:nfnods),right_nodes(1:nfnods))
      !
      rotation = RotateFluxPoints_Tria(left_nodes(1:nfnods), &
                                       right_nodes(1:nfnods))
      !
      ! If the previous function could not find the rotation, mark that there
      ! was an error finding the correct face rotation between cells.
      !
      if (rotation == 0) error_finding_rotation = true
      !
      face(nf)%right%fp_offset = face_pts
      face(nf)%right%rotation = rotation
      !
    end if
    !
  end do
  !
  if (error_finding_rotation) then
    write (iout,*) "There was an error rotating flux points!"
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  !
  call debug_timer(stop_timer,"section 7")
  ! ###########
  ! ###########
  ! ###  8  ###
  ! ###########
  ! ###########
  call debug_timer(start_timer,"section 8")
  !
  ! Now that we have finished creating the face array, use it to find the
  ! the maximum interpolation order required for each combination of
  ! cell geometry and cell order.
  !
  call find_geom_max_interp_order(cell_geom,cell_order,face(:)%left%cell, &
                                  face(:)%right%cell,face(:)%order)
  call debug_timer(stop_timer,"section 8")
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (" Error on face - ",i0,",  nodes ",i0,3(:,2x,i0))
  2 format ("face(",i0,")%",a,"%flxpts")
  !
end subroutine create_face_derived_type
!
!###############################################################################
!
subroutine GetFaceNodes_for_QuadFace(left_cell_nodes,left_side, &
                                     right_cell_nodes,right_side, &
                                     left_face_nodes,right_face_nodes)
  !
  !.. Formal Arguments ..
  integer, dimension(:),    intent(in) :: left_cell_nodes
  integer,                  intent(in) :: left_side
  integer, dimension(:),    intent(in) :: right_cell_nodes
  integer,                  intent(in) :: right_side
  integer, dimension(:), intent(inout) :: left_face_nodes
  integer, dimension(:), intent(inout) :: right_face_nodes
  !
  !.. Local Scalars ..
  integer :: i,left_geom,right_geom
  integer :: left_ncnods,right_ncnods
  !
  !.. Local Arrays ..
  integer, dimension(1:4) :: ip
  !
  !.. Local Parameters ..
  integer, parameter :: flip_left_hand_to_right_hand(1:4) = [1,4,3,2]
  integer, parameter :: min_face_ABCD(1:4) = [4,2,1,3]
  integer, parameter :: max_face_ABCD(1:4) = [2,4,1,3]
  integer, parameter :: ncnod_to_geom(4:8) = [Geom_Quad,Geom_Pyra,Geom_Pris, &
                                              Geom_Unknown,Geom_Hexa]
  !
  character(len=*), parameter :: pname = "GetFaceNodes_for_QuadFace"
  !
continue
  !
  ! Count the number of cell nodes for the left and right cells
  !
  left_ncnods  = size( left_cell_nodes)
  right_ncnods = size(right_cell_nodes)
  !
  ! Make sure the number of nodes for the left cell is between 5 and 8 because
  ! only pyramids, prisms, and hexahedra are valid geometries for the left cell
  !
  if ( left_ncnods < 5 .or. left_ncnods > 8 ) then
    write (error_message,1) ("left",i=1,3),left_ncnods
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Make sure the number of nodes for the right cell is between 4 and 8 because
  ! only pyramids, prisms, hexahedra and quadrilateral boundary faces are valid
  ! geometries for the right cell
  !
  if ( right_ncnods < 4 .or. right_ncnods > 8 ) then
    write (error_message,2) ("right",i=1,3),right_ncnods
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Get the geometry of the left and right cells from the number of nodes
  ! defining the cells
  !
  left_geom  = ncnod_to_geom( left_ncnods)
  right_geom = ncnod_to_geom(right_ncnods)
  !
  ! Get the nodes of the left cell that are on the interface.
  ! NOTE: We have to add one because the values stored in side_nodes are
  !       used as offset values,
  !       i.e. they range from [0:size(left_cell_nodes)-1] and we want them
  !              to range from [1:size(left_cell_nodes)].
  !
  ip = side_nodes(:,left_side,left_geom) + 1
  !
  ! Get the face nodes of the left cell in the correct order needed
  ! to find the rotation between the cells on the interface in the
  ! function RotateFluxPoints_Quad
  !
  select case (left_geom)
      !
    case (Geom_Pyra)
      !
    case (Geom_Pris)
      !
    case (Geom_Hexa)
      !
      ! If the side of the left cell on the interface is either a Y-min or Y-max
      ! face, we need to flip the second and fourth indices. This is because the
      ! vectors 'a' and 'b' for 'Y' faces dont follow the standard rotation of a
      ! right-handed coordinate system like the 'X' and 'Z' face do. To further
      ! illustrate:
      !               For a 'Z' face:   a -> 'X'   and   b -> 'Y'
      !               For a 'X' face:   a -> 'Y'   and   b -> 'Z'
      !
      !               However,
      !   ##CORRECT## For a 'Y' face:   a -> 'X'   and   b -> 'Z' ##CORRECT##
      !               instead of
      ! ##INCORRECT## For a 'Y' face:   a -> 'Z'   and   b -> 'X' ##INCORRECT##
      !
      ! since we want 'a' to correspond to the lowest coordinate direction on
      ! the face. By flipping the second and fourth indices, we can then alter
      ! the face indices as needed based on whether they are simply minimum or
      ! maximum faces.
      !
      if (any(left_side == [Ymin3D_Face,Ymax3D_Face])) then
        ip = ip( flip_left_hand_to_right_hand )
      end if
      !
      ! Rearrange the indices of the face nodes so that they are in the order
      ! [a,b,c,d] depending on whether the side of the left cell on the
      ! interface is a minimum or maximum cell face.
      !
      if (any(left_side == [Xmin3D_Face,Ymin3D_Face,Zmin3D_Face])) then
        ip = ip( min_face_ABCD )
      else
        ip = ip( max_face_ABCD )
      end if
      !
    case default
      !
      ! ERROR! This geometry is invalid
      !
      write (error_message,2) ("left",i=1,2),left_geom
      call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
      !
  end select
  !
  ! Extract the nodes of the left cell on the interface in the order [a,b,c,d]
  !
  left_face_nodes = left_cell_nodes( ip )
  !
  ! Get the nodes of the right cell that are on the interface.
  ! NOTE: We have to add one because the values stored in side_nodes are
  !       used as offset values,
  !       i.e. they range from [0:size(right_cell_nodes)-1] and we want them
  !              to range from [1:size(right_cell_nodes)].
  !
  ip = side_nodes(:,right_side,right_geom) + 1
  !
  ! Get the face nodes of the right cell in the correct order needed
  ! to find the rotation between the cells on the interface in the
  ! function RotateFluxPoints_Quad
  !
  select case (right_geom)
      !
    case (Geom_Quad)
      !
      ! The right cell can be a 2D Quad if this is for a boundary face. For a
      ! quad face, the face indices are always in a minimum face orientation.
      !
      ip = min_face_ABCD
      !
    case (Geom_Pyra)
      !
    case (Geom_Pris)
      !
    case (Geom_Hexa)
      !
      ! If the side of the right cell on the interface is either a Y-min or
      ! Y-max face, we need to flip the second and fourth indices. This is
      ! because the vectors 'a' and 'b' for 'Y' faces dont follow the standard
      ! rotation of a right-handed coordinate system like the 'X' and 'Z' face
      ! do. To further illustrate:
      !               For a 'Z' face:   a -> 'X'   and   b -> 'Y'
      !               For a 'X' face:   a -> 'Y'   and   b -> 'Z'
      !
      !               However,
      !   ##CORRECT## For a 'Y' face:   a -> 'X'   and   b -> 'Z' ##CORRECT##
      !               instead of
      ! ##INCORRECT## For a 'Y' face:   a -> 'Z'   and   b -> 'X' ##INCORRECT##
      !
      ! since we want 'a' to correspond to the lowest coordinate direction on
      ! the face. By flipping the second and fourth indices, we can then alter
      ! the face indices as needed based on whether they are simply minimum or
      ! maximum faces.
      !
      if (any(right_side == [Ymin3D_Face,Ymax3D_Face])) then
        ip = ip( flip_left_hand_to_right_hand )
      end if
      !
      ! Rearrange the indices of the face nodes so that they are in the order
      ! [a,b,c,d] depending on whether the side of the right cell on the
      ! interface is a minimum or maximum cell face.
      !
      if (any(right_side == [Xmin3D_Face,Ymin3D_Face,Zmin3D_Face])) then
        ip = ip( min_face_ABCD )
      else
        ip = ip( max_face_ABCD )
      end if
      !
    case default
      !
      ! ERROR! This geometry is invalid
      !
      write (error_message,2) ("right",i=1,2),right_geom
      call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
      !
  end select
  !
  ! Extract the nodes of the right cell on the interface in the order [a,b,c,d]
  !
  right_face_nodes = right_cell_nodes( ip )
  !
  ! Format Statements
  !
  1 format ("ERROR! The number of cell nodes contained in ",a,"_cell_nodes ", &
            "does not match a valid geometry for the ",a," cell on an ", &
            "interface. The number of cell nodes in ",a,"_cell_nodes is ", &
            i0,".")
  2 format ("ERROR! Encountered an invalid geometry when trying to find ", &
            "the nodes of the ",a," cell on the interface in the order ", &
            "[a,b,c,d]. The value of ",a,"_geom is ",i0".")
  !
end subroutine GetFaceNodes_for_QuadFace
!
!###############################################################################
!
subroutine GetFaceNodes_for_TriaFace(left_cell_nodes,left_side, &
                                     right_cell_nodes,right_side, &
                                     left_face_nodes,right_face_nodes)
  !
  !.. Formal Arguments ..
  integer, dimension(:),    intent(in) :: left_cell_nodes
  integer,                  intent(in) :: left_side
  integer, dimension(:),    intent(in) :: right_cell_nodes
  integer,                  intent(in) :: right_side
  integer, dimension(:), intent(inout) :: left_face_nodes
  integer, dimension(:), intent(inout) :: right_face_nodes
  !
  !.. Local Scalars ..
  integer :: i,left_geom,right_geom
  integer :: left_ncnods,right_ncnods
  !
  !.. Local Arrays ..
  integer, dimension(1:3) :: ip
  !
  !.. Local Parameters ..
  integer, parameter :: ncnod_to_geom(3:6) = [Geom_Tria,Geom_Tetr, &
                                              Geom_Pyra,Geom_Pris]
  !
  character(len=*), parameter :: pname = "GetFaceNodes_for_TriaFace"
  !
continue
  !
  ! Count the number of cell nodes for the left and right cells
  !
  left_ncnods  = size( left_cell_nodes)
  right_ncnods = size(right_cell_nodes)
  !
  ! Make sure the number of nodes for the left cell is between 4 and 6 because
  ! only tetrahedra, pyramids, and prisms are valid geometries for the left cell
  !
  if ( left_ncnods < 4 .or. left_ncnods > 6 ) then
    write (error_message,1) ("left",i=1,3),left_ncnods
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Make sure the number of nodes for the right cell is between 3 and 6 because
  ! only tetrahedra, pyramids, prisms, and triangular boundary faces are valid
  ! geometries for the right cell
  !
  if ( right_ncnods < 3 .or. right_ncnods > 6 ) then
    write (error_message,2) ("right",i=1,3),right_ncnods
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Get the geometry of the left and right cells from the number of nodes
  ! defining the cells
  !
  left_geom  = ncnod_to_geom( left_ncnods)
  right_geom = ncnod_to_geom(right_ncnods)
  !
  ! Get the nodes of the left cell that are on the interface.
  ! NOTE: We have to add one because the values stored in side_nodes are
  !       used as offset values,
  !       i.e. they range from [0:size(left_cell_nodes)-1] and we want them
  !              to range from [1:size(left_cell_nodes)].
  !
  ip = side_nodes(1:3,left_side,left_geom) + 1
  !
  ! Get the face nodes of the left cell in the correct order needed
  ! to find the rotation between the cells on the interface in the
  ! function RotateFluxPoints_Quad
  !
  select case (left_geom)
      !
    case (Geom_Tetr)
      !
    case (Geom_Pyra)
      !
    case (Geom_Pris)
      !
    case default
      !
      ! ERROR! This geometry is invalid
      !
      write (error_message,2) ("left",i=1,2),left_geom
      call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
      !
  end select
  !
  ! Extract the nodes of the left cell on the interface
  !
  left_face_nodes = left_cell_nodes( ip )
  !
  ! Get the nodes of the right cell that are on the interface.
  ! NOTE: We have to add one because the values stored in side_nodes are
  !       used as offset values,
  !       i.e. they range from [0:size(right_cell_nodes)-1] and we want them
  !              to range from [1:size(right_cell_nodes)].
  !
  ip = side_nodes(1:3,right_side,right_geom) + 1
  !
  ! Get the face nodes of the right cell in the correct order needed
  ! to find the rotation between the cells on the interface in the
  ! function RotateFluxPoints_Tria
  !
  select case (right_geom)
      !
    case (Geom_Tria)
      !
      ! The right cell can be a 2D triangle if this is for a boundary face. For
      ! a triangular face, the face indices are always
      !
    case (Geom_Tetr)
      !
    case (Geom_Pyra)
      !
    case (Geom_Pris)
      !
    case default
      !
      ! ERROR! This geometry is invalid
      !
      write (error_message,2) ("right",i=1,2),right_geom
      call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
      !
  end select
  !
  ! Extract the nodes of the right cell on the interface
  !
  right_face_nodes = right_cell_nodes( ip )
  !
  ! Format Statements
  !
  1 format ("ERROR! The number of cell nodes contained in ",a,"_cell_nodes ", &
            "does not match a valid geometry for the ",a," cell on an ", &
            "interface. The number of cell nodes in ",a,"_cell_nodes is ", &
            i0,".")
  2 format ("ERROR! Encountered an invalid geometry when trying to find ", &
            "the nodes of the ",a," cell on the interface in the order ", &
            "[a,b,c]. The value of ",a,"_geom is ",i0".")
  !
end subroutine GetFaceNodes_for_TriaFace
!
!###############################################################################
!
pure function FluxPointRotation_Edge(ab_left,ab_right) result(return_value)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: ab_left
  integer, dimension(:), intent(in) :: ab_right
 !integer,               intent(in) :: face_order
  !
  !.. Function Result ..
  integer :: return_value
  !
  !.. Local Scalars ..
 !integer :: nep
  !
  !.. Local Parameters ..
  integer, parameter :: rotation(1:2,1:2) = reshape( [1,2, 2,1] , [2,2] )
  !
continue
  !
 !nep = face_order+1
  !
  if (all( ab_left == ab_right(rotation(:,1)) )) then
    !
    ! The nodes from the left and right cells are in the same order
    !    a_left = a_right
    !    b_left = b_right
    !
    ! ################
    ! ## ROTATION 1 ##
    ! ################
    !
    ! Therefore, the flux point indices for the right cell are the exact same
    ! as those for the left cell.
    !
    return_value = 1
   !return_value = intseq(0,nep-1)
    !
  else if (all( ab_left == ab_right(rotation(:,2)) )) then
    !
    ! The nodes from the left and right cells are in the opposite order
    !    a_left = b_right
    !    b_left = a_right
    !
    ! ################
    ! ## ROTATION 2 ##
    ! ################
    !
    ! Therefore, the flux point indices for the right cell are the opposite
    ! order as those for the left cell.
    !
    return_value = 2
   !return_value = intseq(nep-1,0,-1)
    !
  else
    !
    ! Something was wrong here so just return the value 0.
    !
    return_value = 0
    !
  end if
  !
end function FluxPointRotation_Edge
!
!###############################################################################
!
pure function RotateFluxPoints_Edge(ab_left,ab_right,face_order) &
                             result(return_value)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: ab_left
  integer, dimension(:), intent(in) :: ab_right
  integer,               intent(in) :: face_order
  !
  !.. Function Result ..
  integer, dimension(1:face_order+1) :: return_value
  !
  !.. Local Scalars ..
  integer :: nep
  !
  !.. Local Parameters ..
  integer, parameter :: rotation(1:2,1:2) = reshape( [1,2, 2,1] , [2,2] )
  !
continue
  !
  nep = face_order+1
  !
  if (all( ab_left == ab_right(rotation(:,1)) )) then
    !
    ! The nodes from the left and right cells are in the same order
    !    a_left = a_right
    !    b_left = b_right
    !
    ! ################
    ! ## ROTATION 1 ##
    ! ################
    !
    ! Therefore, the flux point indices for the right cell are the exact same
    ! as those for the left cell.
    !
    return_value = intseq(0,nep-1)
    !
  else if (all( ab_left == ab_right(rotation(:,2)) )) then
    !
    ! The nodes from the left and right cells are in the opposite order
    !    a_left = b_right
    !    b_left = a_right
    !
    ! ################
    ! ## ROTATION 2 ##
    ! ################
    !
    ! Therefore, the flux point indices for the right cell are the opposite
    ! order as those for the left cell.
    !
    return_value = intseq(nep-1,0,-1)
    !
  else
    !
    ! Something was wrong here so just return the value 0.
    !
    return_value = 0
    !
  end if
  !
end function RotateFluxPoints_Edge
!
!###############################################################################
!
pure function FluxPointRotation_Quad(abcd_left,abcd_right) result(return_value)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: abcd_left
  integer, dimension(:), intent(in) :: abcd_right
 !integer,               intent(in) :: face_order
  !
  !.. Function Result ..
  integer :: return_value
  !
  !.. Local Scalars ..
 !integer :: nfp,nep
  !
  !.. Local Arrays ..
 !integer, dimension(1:face_order+1,1:face_order+1) :: fp_idx
  !
  !.. Local Parameters ..
  integer, parameter :: rotation(1:4,1:8) = &
                                         reshape( [1,2,3,4, 4,3,2,1, &
                                                   3,4,1,2, 2,1,4,3, &
                                                   2,1,3,4, 4,3,1,2, &
                                                   3,4,2,1, 1,2,4,3] , [4,8] )
  !
continue
  !
 !nep = face_order+1
 !nfp = nep*nep
  !
  ! Create the offset for the flux point indices for the left face
  !
 !fp_idx = reshape( intseq(0,nfp-1) , [nep,nep] )
  !
  ! If the left cell on the interface has vectors 'a_left' and 'b_left' with the
  ! following orientation:
  !
  !                   (a_left) A
  !                            |
  !                            |
  !                            -----> (b_left)
  !
  ! Then the following if statements give the possible orientations of the
  ! 'a_right' and 'b_right' vectors from the right cell on the interface.
  !
  if (all(abcd_left == abcd_right(rotation(:,1)))) then
    !
    ! The 'a_right' and 'b_right' vectors of
    ! the right cell are oriented such that:
    !
    !                a_left corresponds to +a_right
    !                            and
    !                b_left corresponds to +b_right
    !
    ! ################    (a_right) A
    ! ## ROTATION 1 ##              |
    ! ################              |
    !                               -----> (b_right)
    !
    ! Therefore, the flux point indices for the right cell are the exact same
    ! as those for the left cell.
    !
    return_value = 1
   !return_value = reshape( fp_idx , [nfp] )
    !
  else if (all(abcd_left == abcd_right(rotation(:,2)))) then
    !
    ! The 'a_right' and 'b_right' vectors of
    ! the right cell are oriented such that:
    !
    !                a_left corresponds to +a_right
    !                            and
    !                b_left corresponds to -b_right
    !
    ! ################                   A (a_right)
    ! ## ROTATION 2 ##                   |
    ! ################                   |
    !                     (b_right) <-----
    !
    ! Therefore, the flux point indices for the right cell are the same ordering
    ! for the first dimension and the opposite ordering for the second.
    !
    return_value = 2
   !return_value = reshape( fp_idx(:,nep:1:-1) , [nfp] )
    !
  else if (all(abcd_left == abcd_right(rotation(:,3)))) then
    !
    ! The 'a_right' and 'b_right' vectors of
    ! the right cell are oriented such that:
    !
    !                a_left corresponds to -a_right
    !                            and
    !                b_left corresponds to +b_right
    !
    ! ################              -----> (b_right)
    ! ## ROTATION 3 ##              |
    ! ################              |
    !                     (a_right) V
    !
    ! Therefore, the flux point indices for the right cell are the opposite
    ! ordering for the first dimension and the same ordering for the second.
    !
    return_value = 3
   !return_value = reshape( fp_idx(nep:1:-1,:) , [nfp] )
    !
  else if (all(abcd_left == abcd_right(rotation(:,4)))) then
    !
    ! The 'a_right' and 'b_right' vectors of
    ! the right cell are oriented such that:
    !
    !                a_left corresponds to -a_right
    !                            and
    !                b_left corresponds to -b_right
    !
    ! ################    (b_right) <-----
    ! ## ROTATION 4 ##                   |
    ! ################                   |
    !                                    V (a_right)
    !
    ! Therefore, the flux point indices for the right cell are the opposite
    ! ordering for both the first and second dimensions
    !
    return_value = 4
   !return_value = reshape( fp_idx(nep:1:-1,nep:1:-1) , [nfp] )
    !
  else if (all(abcd_left == abcd_right(rotation(:,5)))) then
    !
    ! The 'a_right' and 'b_right' vectors of
    ! the right cell are oriented such that:
    !
    !                a_left corresponds to +b_right
    !                            and
    !                b_left corresponds to +a_right
    !
    ! ################    (b_right) A
    ! ## ROTATION 5 ##              |
    ! ################              |
    !                               -----> (a_right)
    !
    ! Because the a and b vectors from the left cell are aligned with the
    ! opposite vectors from the right cell, we first need to transpose the
    ! fp_idx array.
    !
   !fp_idx = transpose(fp_idx)
    !
    ! Both aligned vectors have the same directional signs so keep the same
    ! orderings for both the first and second dimensions.
    !
    return_value = 5
   !return_value = reshape( fp_idx , [nfp] )
    !
  else if (all(abcd_left == abcd_right(rotation(:,6)))) then
    !
    ! The 'a_right' and 'b_right' vectors of
    ! the right cell are oriented such that:
    !
    !                a_left corresponds to +b_right
    !                            and
    !                b_left corresponds to -a_right
    !
    ! ################                   A (b_right)
    ! ## ROTATION 6 ##                   |
    ! ################                   |
    !                     (a_right) <-----
    !
    ! Because the a and b vectors from the left cell are aligned with the
    ! opposite vectors from the right cell, we first need to transpose the
    ! fp_idx array.
    !
   !fp_idx = transpose(fp_idx)
    !
    ! The second aligned vector has opposite sign so keep the same ordering for
    ! the first dimension and reverse the ordering of the second.
    !
    return_value = 6
   !return_value = reshape( fp_idx(:,nep:1:-1) , [nfp] )
    !
  else if (all(abcd_left == abcd_right(rotation(:,7)))) then
    !
    ! The 'a_right' and 'b_right' vectors of
    ! the right cell are oriented such that:
    !
    !                a_left corresponds to -b_right
    !                            and
    !                b_left corresponds to +a_right
    !
    ! ################              -----> (a_right)
    ! ## ROTATION 7 ##              |
    ! ################              |
    !                     (b_right) V
    !
    ! Because the a and b vectors from the left cell are aligned with the
    ! opposite vectors from the right cell, we first need to transpose the
    ! fp_idx array.
    !
   !fp_idx = transpose(fp_idx)
    !
    ! The first aligned vector has opposite sign so reverse the ordering of the
    ! first dimension and keep the same ordering for the second.
    !
    return_value = 7
   !return_value = reshape( fp_idx(nep:1:-1,:) , [nfp] )
    !
  else if (all(abcd_left == abcd_right(rotation(:,8)))) then
    !
    ! The 'a_right' and 'b_right' vectors of
    ! the right cell are oriented such that:
    !
    !                a_left corresponds to -b_right
    !                            and
    !                b_left corresponds to -a_right
    !
    ! ################    (a_right) <-----
    ! ## ROTATION 8 ##                   |
    ! ################                   |
    !                                    V (b_right)
    !
    ! Because the a and b vectors from the left cell are aligned with the
    ! opposite vectors from the right cell, we first need to transpose the
    ! fp_idx array.
    !
   !fp_idx = transpose(fp_idx)
    !
    ! Both of the aligned vectors have opposite signs so reverse the orderings
    ! of both the first and second dimensions.
    !
    return_value = 8
   !return_value = reshape( fp_idx(nep:1:-1,nep:1:-1) , [nfp] )
    !
  else
    !
    ! Something was wrong here so just return the value 0.
    !
    return_value = 0
    !
  end if
  !
end function FluxPointRotation_Quad
!
!###############################################################################
!
pure function RotateFluxPoints_Quad(abcd_left,abcd_right,face_order) &
                             result(return_value)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: abcd_left
  integer, dimension(:), intent(in) :: abcd_right
  integer,               intent(in) :: face_order
  !
  !.. Function Result ..
  integer, dimension(1:(face_order+1)**2) :: return_value
  !
  !.. Local Scalars ..
  integer :: nfp,nep
  !
  !.. Local Arrays ..
  integer, dimension(1:face_order+1,1:face_order+1) :: fp_idx
  !
  !.. Local Parameters ..
  integer, parameter :: rotation(1:4,1:8) = &
                                         reshape( [1,2,3,4, 4,3,2,1, &
                                                   3,4,1,2, 2,1,4,3, &
                                                   2,1,3,4, 4,3,1,2, &
                                                   3,4,2,1, 1,2,4,3] , [4,8] )
  !
continue
  !
  nep = face_order+1
  nfp = nep*nep
  !
  ! Create the offset for the flux point indices for the left face
  !
  fp_idx = reshape( intseq(0,nfp-1) , [nep,nep] )
  !
  ! If the left cell on the interface has vectors 'a_left' and 'b_left' with the
  ! following orientation:
  !
  !                   (a_left) A
  !                            |
  !                            |
  !                            -----> (b_left)
  !
  ! Then the following if statements give the possible orientations of the
  ! 'a_right' and 'b_right' vectors from the right cell on the interface.
  !
  if (all(abcd_left == abcd_right(rotation(:,1)))) then
    !
    ! The 'a_right' and 'b_right' vectors of
    ! the right cell are oriented such that:
    !
    !                a_left corresponds to +a_right
    !                            and
    !                b_left corresponds to +b_right
    !
    ! ################    (a_right) A
    ! ## ROTATION 1 ##              |
    ! ################              |
    !                               -----> (b_right)
    !
    ! Therefore, the flux point indices for the right cell are the exact same
    ! as those for the left cell.
    !
    return_value = reshape( fp_idx , [nfp] )
    !
  else if (all(abcd_left == abcd_right(rotation(:,2)))) then
    !
    ! The 'a_right' and 'b_right' vectors of
    ! the right cell are oriented such that:
    !
    !                a_left corresponds to +a_right
    !                            and
    !                b_left corresponds to -b_right
    !
    ! ################                   A (a_right)
    ! ## ROTATION 2 ##                   |
    ! ################                   |
    !                     (b_right) <-----
    !
    ! Therefore, the flux point indices for the right cell are the same ordering
    ! for the first dimension and the opposite ordering for the second.
    !
    return_value = reshape( fp_idx(:,nep:1:-1) , [nfp] )
    !
  else if (all(abcd_left == abcd_right(rotation(:,3)))) then
    !
    ! The 'a_right' and 'b_right' vectors of
    ! the right cell are oriented such that:
    !
    !                a_left corresponds to -a_right
    !                            and
    !                b_left corresponds to +b_right
    !
    ! ################              -----> (b_right)
    ! ## ROTATION 3 ##              |
    ! ################              |
    !                     (a_right) V
    !
    ! Therefore, the flux point indices for the right cell are the opposite
    ! ordering for the first dimension and the same ordering for the second.
    !
    return_value = reshape( fp_idx(nep:1:-1,:) , [nfp] )
    !
  else if (all(abcd_left == abcd_right(rotation(:,4)))) then
    !
    ! The 'a_right' and 'b_right' vectors of
    ! the right cell are oriented such that:
    !
    !                a_left corresponds to -a_right
    !                            and
    !                b_left corresponds to -b_right
    !
    ! ################    (b_right) <-----
    ! ## ROTATION 4 ##                   |
    ! ################                   |
    !                                    V (a_right)
    !
    ! Therefore, the flux point indices for the right cell are the opposite
    ! ordering for both the first and second dimensions
    !
    return_value = reshape( fp_idx(nep:1:-1,nep:1:-1) , [nfp] )
    !
  else if (all(abcd_left == abcd_right(rotation(:,5)))) then
    !
    ! The 'a_right' and 'b_right' vectors of
    ! the right cell are oriented such that:
    !
    !                a_left corresponds to +b_right
    !                            and
    !                b_left corresponds to +a_right
    !
    ! ################    (b_right) A
    ! ## ROTATION 5 ##              |
    ! ################              |
    !                               -----> (a_right)
    !
    ! Because the a and b vectors from the left cell are aligned with the
    ! opposite vectors from the right cell, we first need to transpose the
    ! fp_idx array.
    !
    fp_idx = transpose(fp_idx)
    !
    ! Both aligned vectors have the same directional signs so keep the same
    ! orderings for both the first and second dimensions.
    !
    return_value = reshape( fp_idx , [nfp] )
    !
  else if (all(abcd_left == abcd_right(rotation(:,6)))) then
    !
    ! The 'a_right' and 'b_right' vectors of
    ! the right cell are oriented such that:
    !
    !                a_left corresponds to +b_right
    !                            and
    !                b_left corresponds to -a_right
    !
    ! ################                   A (b_right)
    ! ## ROTATION 6 ##                   |
    ! ################                   |
    !                     (a_right) <-----
    !
    ! Because the a and b vectors from the left cell are aligned with the
    ! opposite vectors from the right cell, we first need to transpose the
    ! fp_idx array.
    !
    fp_idx = transpose(fp_idx)
    !
    ! The second aligned vector has opposite sign so keep the same ordering for
    ! the first dimension and reverse the ordering of the second.
    !
    return_value = reshape( fp_idx(:,nep:1:-1) , [nfp] )
    !
  else if (all(abcd_left == abcd_right(rotation(:,7)))) then
    !
    ! The 'a_right' and 'b_right' vectors of
    ! the right cell are oriented such that:
    !
    !                a_left corresponds to -b_right
    !                            and
    !                b_left corresponds to +a_right
    !
    ! ################              -----> (a_right)
    ! ## ROTATION 7 ##              |
    ! ################              |
    !                     (b_right) V
    !
    ! Because the a and b vectors from the left cell are aligned with the
    ! opposite vectors from the right cell, we first need to transpose the
    ! fp_idx array.
    !
    fp_idx = transpose(fp_idx)
    !
    ! The first aligned vector has opposite sign so reverse the ordering of the
    ! first dimension and keep the same ordering for the second.
    !
    return_value = reshape( fp_idx(nep:1:-1,:) , [nfp] )
    !
  else if (all(abcd_left == abcd_right(rotation(:,8)))) then
    !
    ! The 'a_right' and 'b_right' vectors of
    ! the right cell are oriented such that:
    !
    !                a_left corresponds to -b_right
    !                            and
    !                b_left corresponds to -a_right
    !
    ! ################    (a_right) <-----
    ! ## ROTATION 8 ##                   |
    ! ################                   |
    !                                    V (b_right)
    !
    ! Because the a and b vectors from the left cell are aligned with the
    ! opposite vectors from the right cell, we first need to transpose the
    ! fp_idx array.
    !
    fp_idx = transpose(fp_idx)
    !
    ! Both of the aligned vectors have opposite signs so reverse the orderings
    ! of both the first and second dimensions.
    !
    return_value = reshape( fp_idx(nep:1:-1,nep:1:-1) , [nfp] )
    !
  else
    !
    ! Something was wrong here so just return the value 0.
    !
    return_value = 0
    !
  end if
  !
end function RotateFluxPoints_Quad
!
!###############################################################################
!
pure function RotateFluxPoints_Tria(left_nodes,right_nodes,nfp) &
                             result(return_value)
  !
  !.. Formal Arguments ..
  integer,               intent(in) :: nfp
  integer, dimension(:), intent(in) :: left_nodes
  integer, dimension(:), intent(in) :: right_nodes
  !
  !.. Function Result ..
  integer, dimension(1:nfp) :: return_value
  !
continue
  !
  return_value = 0
  !
end function RotateFluxPoints_Tria
!
!###############################################################################
!
pure function FluxPointRotation_Tria(left_nodes,right_nodes) &
                              result(return_value)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: left_nodes
  integer, dimension(:), intent(in) :: right_nodes
  !
  !.. Function Result ..
  integer :: return_value
  !
continue
  !
  return_value = 0
  !
end function FluxPointRotation_Tria
!
!###############################################################################
!
subroutine create_cell_derived_type(nodes_of_cell,nodes_of_cell_ptr, &
                                    cells_with_node_ptr,nodes_on_edge, &
                                    cells_with_edge,cells_with_edge_ptr)
  !
  !.. Use Statements ..
  use geovar,    only : ncell,cell,nedge
  use geovar,    only : nfbnd,nface,face
  use geovar,    only : cell_geom
  use geovar,    only : cell_order
  use order_mod, only : geom_solpts
  use ovar,      only : continuous_output
  !
  !.. Formal Arguments ..
  integer,              dimension(:),   intent(in) :: nodes_of_cell
  integer,              dimension(:),   intent(in) :: nodes_of_cell_ptr
  integer, allocatable, dimension(:),   intent(in) :: cells_with_node_ptr
  integer, allocatable, dimension(:),   intent(in) :: cells_with_edge
  integer, allocatable, dimension(:),   intent(in) :: cells_with_edge_ptr
  integer, allocatable, dimension(:,:), intent(in) :: nodes_on_edge
  !
  !.. Local Scalars ..
  integer :: k,l,i,ibeg,iend,interval,ierr
  integer :: n,nc,nc1,nc2,np,np1,np2
  integer :: ncf,nf,nfp,nep,lbeg,rbeg
  integer :: num_faces,num_edges,num_nodes
  integer :: face_geom,face_order,bad_edges
  integer :: host_cell,left_cell,right_cell
  integer :: host_face,left_face,right_face
  integer :: this_cell,this_geom,this_order
  integer :: this_edge,this_node,cell_edge
  real(wp) :: ave
  character(len=200) :: array_name
  !
  !.. Local Arrays ..
  integer, dimension(1:2) :: cp,ep
  integer, dimension(1:8) :: ip
  !
  !.. Local Parameters ..
  integer, parameter :: left_side = 1
  integer, parameter :: right_side = 2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_cell_derived_type"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  !
  !  NOTE: CURRENTLY, EVERYTHING IN THIS ROUTINE HAS THE ARRAY CELL GIVING THE
  !        NUMBER OF FLUX POINTS ON EACH SIDE OF THE FACE BASED ON THE ORDER OF
  !        THE RESPECTIVE CELL ON THAT SIDE. THIS IS PERFECTLY FINE IF BOTH
  !        CELLS ON THE FACE ARE OF THE SAME ORDER. IF THE ORDER IS DIFFERENT
  !        BETWEEN THE CELLS ON A FACE, WE MUST BE CAREFUL WHEN USING THIS ARRAY
  !        TO FILL ARRAYS LIKE FACEUSP BY INTERPOLATING INFORMATION TO THE
  !        FACES. THOSE FACE ARRAYS CURRENTLY ASSUME THE NUMBER OF FLUX POINTS
  !        ON BOTH SIDES OF THE FACE IS THE SAME. THIS ASSUMPTION IS DESIRED AND
  !        THE NUMBER OF FLUX POINTS ON BOTH SIDES SHOULD BE ASSOCIATED WITH THE
  !        ORDER OF THE CELL THAT IS OF HIGHER ORDER. THIS MEANS THE CELL OF
  !        LOWER ORDER MUST INTERPOLATE ITS POLYNOMIAL SOLUTION TO THE HIGHER
  !        ORDER FLUX POINTS OF THE ADJACENT CELL. LONG STORY SHORT, THE
  !        ASSUMPTION MADE IN THIS ROUTINE WILL BE INCORRECT FOR ONE CELL ON A
  !        FACE IF THE ORDERS OF THE ATTACHED CELLS ARE DIFFERENT.
  !
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  ! ############################################################################
  !
  !
  ! Allocate the outer part of the cell array
  !
  allocate ( cell(1:ncell) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Go through the cells and allocate the face, edge, and node parts of cell.
  ! While we are doing this, store the geometry, solution order, first solution
  ! point, and last solution point of each cell within the cell array.
  !
  iend = 0
  !
  do nc = 1,ncell
    !
    this_geom = cell_geom(nc)
    this_order = cell_order(nc)
    !
    ibeg = iend + 1
    iend = iend + geom_solpts( this_geom , this_order )
    !
    cell(nc)%geom = this_geom
    cell(nc)%order = this_order
    cell(nc)%beg_sp = ibeg
    cell(nc)%end_sp = iend
    !
    num_faces = geom_faces(this_geom)
    !
    allocate ( cell(nc)%face(1:num_faces) , stat=ierr , errmsg=error_message )
    write (array_name,1) nc
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !
    if (continuous_output) then
      !
      num_edges = geom_edges(this_geom)
      !
      allocate ( cell(nc)%edge(1:num_edges) , stat=ierr , errmsg=error_message )
      write (array_name,2) nc
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                       error_message,skip_alloc_pause)
      !
      num_nodes = geom_nodes(this_geom)
      !
      allocate ( cell(nc)%node(1:num_nodes) , stat=ierr , errmsg=error_message )
      write (array_name,3) nc
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                       error_message,skip_alloc_pause)
      !
    end if
    !
  end do
  !
  ! Loop through the faces and create the flx-pt information for each
  ! cell attached to the face.
  !
  ! Start with only the boundary faces
  !
  do nf = 1,nfbnd
    !
    host_cell  = face(nf)%left%cell
    host_face  = face(nf)%left%cell_face
    face_geom  = face(nf)%geom
    face_order = face(nf)%order
    !
    ! Copy in the face information for the current cell face
    !
    cell(host_cell)%face(host_face)%idx = nf
    cell(host_cell)%face(host_face)%side = left_side
    cell(host_cell)%face(host_face)%geom = face_geom
    cell(host_cell)%face(host_face)%order = face_order
    !
    cell(host_cell)%face(host_face)%fp_offset = face(nf)%left%fp_offset
    cell(host_cell)%face(host_face)%rotation = face(nf)%left%rotation
    !
  end do
  !
  ! Now do the interior faces
  !
  do nf = nfbnd+1,nface
    !
    ! Information for this face and its attached cells
    !
    left_cell  = face(nf)%left%cell
    right_cell = face(nf)%right%cell
    left_face  = face(nf)%left%cell_face
    right_face = face(nf)%right%cell_face
    face_geom  = face(nf)%geom
    face_order = face(nf)%order
    !
    ! ####################################
    ! ####################################
    ! #####  Left cell on this face  #####
    ! ####################################
    ! ####################################
    !
    ! Copy in the face information for the left cell face.
    ! For now, initialize the face_idx component to [1:nfp]
    !
    cell(left_cell)%face(left_face)%idx = nf
    cell(left_cell)%face(left_face)%side = left_side
    cell(left_cell)%face(left_face)%geom = face_geom
    cell(left_cell)%face(left_face)%order = face_order
    !
    cell(left_cell)%face(left_face)%fp_offset = face(nf)%left%fp_offset
    cell(left_cell)%face(left_face)%rotation = face(nf)%left%rotation
    !
    ! #####################################
    ! #####################################
    ! #####  Right cell on this face  #####
    ! #####################################
    ! #####################################
    !
    ! Copy in the face information for the right cell face.
    ! For now, initialize the face_idx component to [1:nfp]
    !
    cell(right_cell)%face(right_face)%idx = nf
    cell(right_cell)%face(right_face)%side = right_side
    cell(right_cell)%face(right_face)%geom = face_geom
    cell(right_cell)%face(right_face)%order = face_order
    !
    cell(right_cell)%face(right_face)%fp_offset = face(nf)%right%fp_offset
    cell(right_cell)%face(right_face)%rotation = face(nf)%right%rotation
    !
  end do
  !
  if (.not. continuous_output) then
    call debug_timer(leaving_procedure,pname)
    return
  end if
  !
  ! Loop through the edges and create the edg-pt information for each
  ! cell attached to the edge
  !
  nedge = size(nodes_on_edge,dim=2)
  !
  bad_edges = 0
  !
  do this_edge = 1,size(nodes_on_edge,dim=2)
    !
    ep(1:2) = nodes_on_edge(1:2,this_edge)
    !
    nc1 = cells_with_edge_ptr(this_edge)+1
    nc2 = cells_with_edge_ptr(this_edge+1)
    ave = one / real(nc2-nc1+1,kind=wp)
    !
    cell_loop: do nc = nc1,nc2
      !
      this_cell = cells_with_edge(nc)
      !
      if (this_cell > ncell) cycle cell_loop
      !
      this_geom = cell(this_cell)%geom
      this_order = cell(this_cell)%order
      !
      np1 = nodes_of_cell_ptr(this_cell)+1
      np2 = nodes_of_cell_ptr(this_cell+1)
      np  = np2-np1+1
      !
      nep = geom_solpts(Geom_Edge,this_order)
      !
      ip(1:np) = nodes_of_cell(np1:np2)
      !
      ! Loop through the edges of this cell to find the edge with the same nodes
      !
      one_pass: do l = 1,1
        !
        find_cell_edge: do i = 1,geom_edges(this_geom)
          !
          cp(1:2) = ip( edge_nodes(1:2,i,this_geom) + 1 )
          !
          if (cp(1)==ep(1) .and. cp(2)==ep(2)) then
            !
            cell_edge = i
            ibeg = nep*(i-1) + 1
            iend = nep*i
            interval = 1
            !
            exit one_pass
            !
          else if (cp(1)==ep(2) .and. cp(2)==ep(1)) then
            !
            cell_edge = i
            ibeg = nep*i
            iend = nep*(i-1) + 1
            interval = -1
            !
            exit one_pass
            !
          else
            !
            ! This nodes for this cell edge do not match so cycle to the
            ! next cell edge
            !
            cycle find_cell_edge
            !
          end if
          !
        end do find_cell_edge
        !
        if (bad_edges == 0) write (iout,10)
        !
        bad_edges = bad_edges + 1
        write (iout,11) this_edge,this_cell
        !
      end do one_pass
      !
      !
      !
      cell(this_cell)%edge(cell_edge)%idx = this_edge
      cell(this_cell)%edge(cell_edge)%ave_coef = ave
      !
      cell(this_cell)%edge(cell_edge)%ep_offset = nep*(cell_edge-1)
      cell(this_cell)%edge(cell_edge)%rev_dir = (interval == -1)
      !
    end do cell_loop
    !
  end do
  !
  ! Kill the program with an error message if any bad edges were found
  !
  if (bad_edges /= 0) then
    write (error_message,12)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Loop through the cells and create the information for each node of the cell
  !
  do this_cell = 1,ncell
    !
    np1 = nodes_of_cell_ptr(this_cell)+1
    np2 = nodes_of_cell_ptr(this_cell+1)
    !
    do n = np1,np2
      !
      this_node = nodes_of_cell(n)
      !
      nc1 = cells_with_node_ptr(this_node)+1
      nc2 = cells_with_node_ptr(this_node+1)
      !
      l = n-np1+1
      !
      cell(this_cell)%node(l)%idx = this_node
      cell(this_cell)%node(l)%ave_coef = one / real(nc2-nc1+1,kind=wp)
      !
    end do
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("cell(",i0,")%face",:,"(",i0,")%",a)
  2 format ("cell(",i0,")%edge",:,"(",i0,")%",a)
  3 format ("cell(",i0,")%node",:,"(",i0,")%",a)
  !
  10 format (/,"ERROR! Bad edges detected!",/)
  11 format ("   BAD EDGE: this_edge = ",i8,";  this_cell = ",i8)
  12 format ("ERROR! Bad edges detected!")
  !
end subroutine create_cell_derived_type
!
!###############################################################################
!
subroutine get_unit_cell_subconnectivity()
  !
  !.. Use Statements ..
  use geovar,    only : pst_geom
  use geovar,    only : global_cell_order
  use ovar,      only : continuous_output
  use order_mod, only : o_order
  !
  !.. Local Scalars ..
  integer :: i,j,k,sk,ierr
  integer :: no1d,no2d,no3d,np1d
  integer :: p1,p2,p3,p4,p5,p6,p7,p8
  integer :: omin,omax,order
  integer :: add_points
  !
  !.. Local Arrays ..
  integer, allocatable, dimension(:,:) :: counter
  integer, allocatable, dimension(:,:,:) :: counter3
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_unit_cell_subconnectivity"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate the pst_geom array
  !
  omin = min( minval( global_cell_order ) , o_order )
  omax = max( maxval( global_cell_order ) , o_order )
  !
  allocate ( pst_geom(Geom_Min:Geom_Max,omin:omax) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"pst_geom",1,__LINE__,__FILE__,ierr, &
                   error_message,skip_alloc_pause)
  !
  ! Loop through all the possible cell orders to get the connectivity
  ! for each of them
  !
  add_points = merge(2,0,continuous_output)
  !
  loop_over_orders: do order = omin,omax
    !
    no1d = max(1,order) + add_points
    no2d = no1d**2
    no3d = no1d**3
    np1d = no1d + 1
    !
    ! Allocate counter and counter3 for this order
    !
    allocate ( counter(1:no1d+1,1:no1d+1) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"counter",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    allocate ( counter3(1:no1d+1,1:no1d+1,1:no1d+1) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"counter3",1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !
    ! ###########################
    ! ## IGNORE NODAL GEOMETRY ##
    ! ###########################
    !
    !
    ! ##################################################################
    ! ##################################################################
    ! ###                       1D GEOMETRIES                        ###
    ! ##################################################################
    ! ##################################################################
    !
    !
    ! ###################
    ! ## Edge Geometry ##
    ! ###################
    !
    allocate ( pst_geom(Geom_Edge,order)%connectivity(1:2,1:no1d) , &
               source=0 , stat=ierr , errmsg=error_message )
    !
    do i = 1,no1d
      pst_geom(Geom_Edge,order)%connectivity(1,i) = i
      pst_geom(Geom_Edge,order)%connectivity(2,i) = i + 1
    end do
    !
    ! Subtract one from the sub-connectivity so that it starts at zero
    !
    pst_geom(Geom_Edge,order)%connectivity = &
                                   pst_geom(Geom_Edge,order)%connectivity - 1
    !
    !
    ! ##################################################################
    ! ##################################################################
    ! ###                       2D GEOMETRIES                        ###
    ! ##################################################################
    ! ##################################################################
    !
    !
    ! #######################
    ! ## Triangle Geometry ##
    ! #######################
    !
    allocate ( pst_geom(Geom_Tria,order)%connectivity(1:4,1:no2d) , &
               source=0 , stat=ierr , errmsg=error_message )
    !
    ! Get the subconnectivity for the unit triangle
    !
    sk = 0
    do j = 1,no1d+1
      do i = 1,no1d+2-j
      sk = sk + 1
      counter(j,i) = sk
      end do
    end do
    !
    sk = 0
    do j = 1,no1d+1
      do i = 1,no1d+1-j
        !
        p1 = counter(j,i)
        p2 = counter(j,i+1)
        p3 = counter(j+1,i)
        p4 = counter(j+1,i+1)
        !
        sk = sk + 1
        pst_geom(Geom_Tria,order)%connectivity(1,sk) = p1
        pst_geom(Geom_Tria,order)%connectivity(2,sk) = p2
        pst_geom(Geom_Tria,order)%connectivity(3,sk) = p3
        pst_geom(Geom_Tria,order)%connectivity(4,sk) = p3
        !
        if (p4 /= 0) then
          sk = sk + 1
          pst_geom(Geom_Tria,order)%connectivity(1,sk) = p2
          pst_geom(Geom_Tria,order)%connectivity(2,sk) = p4
          pst_geom(Geom_Tria,order)%connectivity(3,sk) = p3
          pst_geom(Geom_Tria,order)%connectivity(4,sk) = p3
        end if
        !
      end do
    end do
    !
    ! Subtract one from the sub-connectivity so that it starts at zero
    !
    pst_geom(Geom_Tria,order)%connectivity = &
                                   pst_geom(Geom_Tria,order)%connectivity - 1
    !
    ! ############################
    ! ## Quadrilateral Geometry ##
    ! ############################
    !
    allocate ( pst_geom(Geom_Quad,order)%connectivity(1:4,1:no2d) , &
               source=0 , stat=ierr , errmsg=error_message )
    !
    ! Get the subconnectivity for the unit square
    !
    sk = 0
    do j = 1,no1d
      do i = 1,no1d
        !
        sk = sk + 1
        pst_geom(Geom_Quad,order)%connectivity(1,sk) = (no1d+1)*(j-1) + i
        pst_geom(Geom_Quad,order)%connectivity(2,sk) = (no1d+1)*(j-1) + i + 1
        pst_geom(Geom_Quad,order)%connectivity(3,sk) = (no1d+1)*(j  ) + i + 1
        pst_geom(Geom_Quad,order)%connectivity(4,sk) = (no1d+1)*(j  ) + i
        !
      end do
    end do
    !
    ! Subtract one from the sub-connectivity so that it starts at zero
    !
    pst_geom(Geom_Quad,order)%connectivity = &
                                   pst_geom(Geom_Quad,order)%connectivity - 1
    !
    !
    ! ##################################################################
    ! ##################################################################
    ! ###                       3D GEOMETRIES                        ###
    ! ##################################################################
    ! ##################################################################
    !
    !
    ! ##########################
    ! ## Tetrahedral Geometry ##
    ! ##########################
    !
    allocate ( pst_geom(Geom_Tetr,order)%connectivity(1:8,1:no3d) , &
               source=0 , stat=ierr , errmsg=error_message )
    !
    sk = 0
    do k = 1,np1d
      do j = 1,np1d+1-k
        do i = 1,np1d+2-k-j
          sk = sk + 1
          counter3(k,j,i) = sk
        end do
      end do
    end do
    !
    sk = 0
    sk = 0
    do k = 1,no1d
      do j = 1,no1d+1-k
        do i = 1,no1d+2-k-j
          !
          p1 = counter3(k  ,j  ,i  )
          p2 = counter3(k+1,j  ,i  )
          p3 = counter3(k+1,j+1,i  )
          p4 = counter3(k  ,j+1,i  )
          p5 = counter3(k  ,j  ,i+1)
          p6 = counter3(k+1,j  ,i+1)
          p7 = counter3(k+1,j+1,i+1)
          p8 = counter3(k  ,j+1,i+1)
          !
          sk = sk + 1
          pst_geom(Geom_Tetr,order)%connectivity(1:1,sk) = p1
          pst_geom(Geom_Tetr,order)%connectivity(2:2,sk) = p2
          pst_geom(Geom_Tetr,order)%connectivity(3:4,sk) = p4
          pst_geom(Geom_Tetr,order)%connectivity(5:8,sk) = p5
          !
          if (p6 /= 0) then
            !
            sk = sk + 1
            pst_geom(Geom_Tetr,order)%connectivity(1:1,sk) = p2
            pst_geom(Geom_Tetr,order)%connectivity(2:2,sk) = p4
            pst_geom(Geom_Tetr,order)%connectivity(3:4,sk) = p5
            pst_geom(Geom_Tetr,order)%connectivity(5:8,sk) = p6
            !
            sk = sk + 1
            pst_geom(Geom_Tetr,order)%connectivity(1:1,sk) = p4
            pst_geom(Geom_Tetr,order)%connectivity(2:2,sk) = p5
            pst_geom(Geom_Tetr,order)%connectivity(3:4,sk) = p6
            pst_geom(Geom_Tetr,order)%connectivity(5:8,sk) = p8
            !
            sk = sk + 1
            pst_geom(Geom_Tetr,order)%connectivity(1:1,sk) = p2
            pst_geom(Geom_Tetr,order)%connectivity(2:2,sk) = p3
            pst_geom(Geom_Tetr,order)%connectivity(3:4,sk) = p4
            pst_geom(Geom_Tetr,order)%connectivity(5:8,sk) = p6
            !
            sk = sk + 1
            pst_geom(Geom_Tetr,order)%connectivity(1:1,sk) = p4
            pst_geom(Geom_Tetr,order)%connectivity(2:2,sk) = p6
            pst_geom(Geom_Tetr,order)%connectivity(3:4,sk) = p3
            pst_geom(Geom_Tetr,order)%connectivity(5:8,sk) = p8
            !
            if (p7 /= 0) then
              sk = sk + 1
              pst_geom(Geom_Tetr,order)%connectivity(1:1,sk) = p3
              pst_geom(Geom_Tetr,order)%connectivity(2:2,sk) = p8
              pst_geom(Geom_Tetr,order)%connectivity(3:4,sk) = p6
              pst_geom(Geom_Tetr,order)%connectivity(5:8,sk) = p7
            end if
            !
          end if
          !
        end do
      end do
    end do
    !
    ! Subtract one from the sub-connectivity so that it starts at zero
    !
    pst_geom(Geom_Tetr,order)%connectivity = &
                                   pst_geom(Geom_Tetr,order)%connectivity - 1
    !
    ! ######################
    ! ## Pyramid Geometry ##
    ! ######################
    !
    allocate ( pst_geom(Geom_Pyra,order)%connectivity(1:8,1:no3d) , &
               source=0 , stat=ierr , errmsg=error_message )
    !
    ! Subtract one from the sub-connectivity so that it starts at zero
    !
    pst_geom(Geom_Pyra,order)%connectivity = &
                                   pst_geom(Geom_Pyra,order)%connectivity - 1
    !
    ! ####################
    ! ## Prism Geometry ##
    ! ####################
    !
    allocate ( pst_geom(Geom_Pris,order)%connectivity(1:8,1:no3d) , &
               source=0 , stat=ierr , errmsg=error_message )
    !
    counter3 = 0
    sk = 0
    do k = 1,np1d
      do j = 1,np1d
        do i = 1,np1d+1-k
          sk = sk + 1
          counter3(k,j,i) = sk
        end do
      end do
    end do
    !
    sk = 0
    do k = 1,no1d
      do j = 1,no1d
        do i = 1,no1d+1-k
          !
          p1 = counter3(k  ,j  ,i  )
          p2 = counter3(k  ,j  ,i+1)
          p3 = counter3(k  ,j+1,i+1)
          p4 = counter3(k  ,j+1,i  )
          p5 = counter3(k+1,j  ,i  )
          p6 = counter3(k+1,j  ,i+1)
          p7 = counter3(k+1,j+1,i+1)
          p8 = counter3(k+1,j+1,i  )
          !
          sk = sk + 1
          pst_geom(Geom_Pris,order)%connectivity(1:1,sk) = p1
          pst_geom(Geom_Pris,order)%connectivity(2:2,sk) = p2
          pst_geom(Geom_Pris,order)%connectivity(3:4,sk) = p5
          pst_geom(Geom_Pris,order)%connectivity(5:5,sk) = p4
          pst_geom(Geom_Pris,order)%connectivity(6:6,sk) = p3
          pst_geom(Geom_Pris,order)%connectivity(7:8,sk) = p8
          !
          if (p6 /= 0 .and. p7 /= 0) then
            sk = sk + 1
            pst_geom(Geom_Pris,order)%connectivity(1:1,sk) = p2
            pst_geom(Geom_Pris,order)%connectivity(2:2,sk) = p6
            pst_geom(Geom_Pris,order)%connectivity(3:4,sk) = p5
            pst_geom(Geom_Pris,order)%connectivity(5:5,sk) = p3
            pst_geom(Geom_Pris,order)%connectivity(6:6,sk) = p7
            pst_geom(Geom_Pris,order)%connectivity(7:8,sk) = p8
          end if
          !
        end do
      end do
    end do
    !
    ! Subtract one from the sub-connectivity so that it starts at zero
    !
    pst_geom(Geom_Pris,order)%connectivity = &
                                   pst_geom(Geom_Pris,order)%connectivity - 1
    !
    ! #########################
    ! ## Hexahedral Geometry ##
    ! #########################
    !
    allocate ( pst_geom(Geom_Hexa,order)%connectivity(1:8,1:no3d) , &
               source=0 , stat=ierr , errmsg=error_message )
    !
    sk = 0
    do k = 1,no1d
      do j = 1,no1d
        do i = 1,no1d
          !
          p1 = np1d*np1d*(k-1) + np1d*(j-1) + i
          p2 = np1d*np1d*(k-1) + np1d*(j-1) + i + 1
          p3 = np1d*np1d*(k-1) + np1d*(j  ) + i + 1
          p4 = np1d*np1d*(k-1) + np1d*(j  ) + i
          p5 = np1d*np1d*(k  ) + np1d*(j-1) + i
          p6 = np1d*np1d*(k  ) + np1d*(j-1) + i + 1
          p7 = np1d*np1d*(k  ) + np1d*(j  ) + i + 1
          p8 = np1d*np1d*(k  ) + np1d*(j  ) + i
          !
          sk = sk + 1
          pst_geom(Geom_Hexa,order)%connectivity(1,sk) = p1
          pst_geom(Geom_Hexa,order)%connectivity(2,sk) = p2
          pst_geom(Geom_Hexa,order)%connectivity(3,sk) = p3
          pst_geom(Geom_Hexa,order)%connectivity(4,sk) = p4
          pst_geom(Geom_Hexa,order)%connectivity(5,sk) = p5
          pst_geom(Geom_Hexa,order)%connectivity(6,sk) = p6
          pst_geom(Geom_Hexa,order)%connectivity(7,sk) = p7
          pst_geom(Geom_Hexa,order)%connectivity(8,sk) = p8
          !
        end do
      end do
    end do
    !
    ! Subtract one from the sub-connectivity so that it starts at zero
    !
    pst_geom(Geom_Hexa,order)%connectivity = &
                                   pst_geom(Geom_Hexa,order)%connectivity - 1
    !
    ! Deallocate counter and counter3 before going on to the next order
    !
    deallocate ( counter , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"counter",2,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    deallocate ( counter3 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"counter3",2,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !
  end do loop_over_orders
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_unit_cell_subconnectivity
!
!###############################################################################
!
subroutine find_nodes_on_each_edge(cell_geom, &
                                   nodes_of_cell,nodes_of_cell_ptr, &
                                   cells_with_node,cells_with_node_ptr, &
                                   nodes_on_edge, &
                                   edges_with_node,edges_with_node_ptr)
  !
  !.. Formal Arguments ..
  integer,              dimension(:),      intent(in) :: cell_geom
  integer,              dimension(:),      intent(in) :: nodes_of_cell
  integer,              dimension(:),      intent(in) :: nodes_of_cell_ptr
  integer,              dimension(:),      intent(in) :: cells_with_node
  integer,              dimension(:),      intent(in) :: cells_with_node_ptr
  integer, allocatable, dimension(:,:), intent(inout) :: nodes_on_edge
  !
  !.. Optional Arguments ..
  integer, optional, allocatable, intent(inout) :: edges_with_node(:)
  integer, optional, allocatable, intent(inout) :: edges_with_node_ptr(:)
  !
  !.. Local Scalars ..
  integer :: i,ibeg,iend,ierr,idx,np,n1,n2,lnode,ledge,lon
  integer :: this_node,this_cell,this_edge,this_geom
  integer :: other_node,unique_edges,next_edge_idx
  integer :: max_lon,max_ledge
  !
  !.. Local Arrays ..
  integer, dimension(1:8) :: ip
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: list_of_other_nodes(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_nodes_on_each_edge"
  !
continue
  !
  lnode = maxval(nodes_of_cell)
  !
  ! Before we start trying to find the edges, check to make sure that
  ! the grid contains 3D elements. If the grid is 2D, the edges are
  ! actually the faces so set the number of edges to zero. We will
  ! still allocate the edge arrays (to zero size) because they will
  ! be still be accessed later.
  !
  if (all(geom_dimen(cell_geom) <= 2)) then
    !
    ledge = 0
    !
    ! Allocate the nodes_on_edge array
    !
    allocate ( nodes_on_edge(1:2,1:ledge) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodes_on_edge",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    ! Allocate the edges_with_node_ptr and edges_with_node
    ! arrays if they are present
    !
    if (present(edges_with_node_ptr) .and. present(edges_with_node)) then
      !
      allocate ( edges_with_node_ptr(1:lnode+1) , source=0 , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"edges_with_node_ptr",1,__LINE__,__FILE__,ierr, &
                       error_message)
      !
      allocate ( edges_with_node(1:last(edges_with_node_ptr)) , source=0 , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"edges_with_node",1,__LINE__,__FILE__,ierr, &
                       error_message)
      !
    end if
    !
    return
    !
  end if
  !
  call debug_timer(entering_procedure,pname)
  !
  !##################################################################
  !##################################################################
  !##################################################################
  !#####  1. Count the total number of unique edges that exist  #####
  !#####     and create the pointer array needed to create the  #####
  !#####     list of edges attached to each grid node.          #####
  !##################################################################
  !##################################################################
  !##################################################################
  !
  ! Loop through the grid nodes and count the number of edges that
  ! contain each node
  !
 !unique_edges = 0
 !list_of_other_nodes(:) = 0
 !!
 !call debug_timer(start_timer,"counting number of edges")
 !node_count_loop: do this_node = 1,lnode
 !  !
 !  ! Beginning and ending indices for the cells surrounding the current node
 !  !
 !  ibeg = cells_with_node_ptr(this_node)+1
 !  iend = cells_with_node_ptr(this_node+1)
 !  !
 !  ! Loop over the cells that surround this node
 !  !
 !  lon = 0
 !  !
 !  cell_count_loop: do i = ibeg,iend
 !    !
 !    ! NOTE: We need the absolute value of this because the cell index could
 !    !       be negative if it is for a boundary face ghost cell and the
 !    !       optional argument make_bnd_faces_negative was present with the
 !    !       value .TRUE. when subroutine find_cells_surrounding_each_node
 !    !       was called
 !    !
 !    this_cell = abs(cells_with_node(i))
 !    this_geom = cell_geom(this_cell)
 !    !
 !    n1 = nodes_of_cell_ptr(this_cell)+1
 !    n2 = nodes_of_cell_ptr(this_cell+1)
 !    np = n2-n1+1
 !    !
 !    ! Get the indices for all nodes defining this cell
 !    !
 !    ip(1:np) = nodes_of_cell(n1:n2)
 !    !
 !    ! Get the local node index (i.e., 1<=idx<=geom_nodes(this_geom))
 !    ! of this_node for this_cell
 !    !
 !    idx = minloc( ip(1:np) , dim=1 , mask=(this_node==ip(1:np)) )
 !    !
 !    ! Loop through the cell edges connected to this node
 !    !
 !    edge_count_loop: do this_edge = 1,num_node_edges(idx,this_geom)
 !      !
 !      ! Get the index for the other node that is connected to
 !      ! this_node through this_edge
 !      !
 !      other_node = ip( node_nodes(this_edge,idx,this_geom)+1 )
 !      !
 !      ! Skip to the next edge of this cell if this edge has already been found
 !      !
 !      if ( any(list_of_other_nodes(1:lon) == other_node) ) then
 !        cycle edge_count_loop
 !      end if
 !      !
 !      ! Add this edge to the edge counter ONLY IF the index for the other
 !      ! node is GREATER than the index for the current node. This ensures
 !      ! that the counted edges remain unique, i.e. a-b and b-a are the same
 !      ! edge and arent counted twice.
 !      !
 !      if (other_node > this_node) then
 !        unique_edges = unique_edges + 1
 !      end if
 !      !
 !      ! Mark the edge between this_node and other_node as already found
 !      !
 !      lon = lon + 1
 !      if (lon > max_lon) then
 !        write (error_message,3) max_lon
 !        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
 !      end if
 !      list_of_other_nodes(lon) = other_node
 !      !
 !    end do edge_count_loop
 !    !
 !  end do cell_count_loop
 !  !
 !  ! Reset the list_of_other_nodes array for the next grid node
 !  !
 !  list_of_other_nodes(1:lon) = 0
 !  !
 !end do node_count_loop
 !call debug_timer(stop_timer,"counting number of edges")
 !!
 !! Save the number of unique edges in ledge
 !!
 !ledge = unique_edges
 !write (iout,2) mypnum,"ledge = ",ledge
  !
  !###################################################################
  !###################################################################
  !###################################################################
  !#####  2. Repeat step 1, but this time store the nodes that   #####
  !#####     define each unique edge in the array nodes_on_edge  #####
  !###################################################################
  !###################################################################
  !###################################################################
  !
  ! Allocate the nodes_on_edge array
  ! NOTE: We start by allocating the second dimension of nodes_on_edge to
  !       3*lnodes because the ratio of edges/nodes for a structured block
  !       approaches the limit of 3 as the number of nodes in the i, j, and k
  !       directions go to infinity. We will dynamically increase the size of
  !       this array as needed.
  !
  !                           number_of_edges(i,j,k)
  !             LIMIT         ---------------------- = 3
  !       i,j,k -> infinity   number_of_nodes(i,j,k)
  !
  max_ledge = 3*lnode
  allocate ( nodes_on_edge(1:2,1:max_ledge) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_on_edge",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate the list_of_other_nodes array
  ! NOTE : We start by allocating the list_of_other_nodes array to 100 since
  !        this array should be dimensioned to the maximum number of edges for
  !        a given node. Theoretically, this should be infinity, but a node
  !        with 100 edges would mean there are a lot of highly anisotropic
  !        elements containing this node and the grid quality should probably
  !        be questioned. We will dynamically increase the size of this array
  !        as needed.
  !
  max_lon = 100
  allocate ( list_of_other_nodes(1:max_lon) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"list_of_other_nodes",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! Loop back through the grid nodes, but this time store the edge information
  ! now that we have properly allocated the necessary arrays
  !
  unique_edges = 0
  !
 !call debug_timer(start_timer,"storing the edges")
  node_store_loop: do this_node = 1,lnode
    !
    ! Beginning and ending indices for the cells surrounding the current node
    !
    ibeg = cells_with_node_ptr(this_node)+1
    iend = cells_with_node_ptr(this_node+1)
    !
    ! Loop over the cells that surround this node
    !
    lon = 0
    !
    cell_store_loop: do i = ibeg,iend
      !
      ! NOTE: We need the absolute value of this because the cell index could
      !       be negative if it is for a boundary face ghost cell and the
      !       optional argument make_bnd_faces_negative was present with the
      !       value .TRUE. when subroutine find_cells_surrounding_each_node
      !       was called
      !
      this_cell = abs(cells_with_node(i))
      this_geom = cell_geom(this_cell)
      !
      n1 = nodes_of_cell_ptr(this_cell)+1
      n2 = nodes_of_cell_ptr(this_cell+1)
      np = n2-n1+1
      !
      ! Get the indices for all nodes defining this cell
      !
      ip(1:np) = nodes_of_cell(n1:n2)
      !
      ! Get the local node index (i.e., 1<=idx<=geom_nodes(this_geom))
      ! of this_node for this_cell
      !
      idx = minloc( ip(1:np) , dim=1 , mask=(this_node==ip(1:np)) )
      !
      ! Loop through the cell edges connected to this node
      !
      edge_store_loop: do this_edge = 1,num_node_edges(idx,this_geom)
        !
        ! Get the index for the other node that is connected to
        ! this_node through this_edge
        !
        other_node = ip( node_nodes(this_edge,idx,this_geom)+1 )
        !
        ! Skip to the next edge of this cell if this edge has already been found
        !
        if ( any(list_of_other_nodes(1:lon) == other_node) ) then
          cycle edge_store_loop
        end if
        !
        ! Add this edge to the edge counter ONLY IF the index for the other
        ! node is GREATER than the index for the current node. This ensures
        ! that the counted edges remain unique, i.e. a-b and b-a are the same
        ! edge and arent counted twice.
        !
        if (other_node > this_node) then
          !
          unique_edges = unique_edges + 1
          !
          ! Increase the size of nodes_on_edge if there is not enough space
          ! to store this edge in nodes_on_edge
          !
          if (unique_edges > max_ledge) then
            max_ledge = max_ledge + lnode
            call reallocate(nodes_on_edge,2,max_ledge)
          end if
          !
          nodes_on_edge(1,unique_edges) = this_node
          nodes_on_edge(2,unique_edges) = other_node
          !
        end if
        !
        ! Mark the edge between this_node and other_node as already found
        !
        lon = lon + 1
        !
        ! Increase the size of list_of_other_nodes if there is not enough
        ! space to add this other node to list_of_other_nodes
        !
        if (lon > max_lon) then
          max_lon = max_lon + 100
          call reallocate(list_of_other_nodes,max_lon)
        end if
        !
        list_of_other_nodes(lon) = other_node
        !
      end do edge_store_loop
      !
    end do cell_store_loop
    !
    ! Reset the list_of_other_nodes array for the next grid node
    !
    list_of_other_nodes(1:lon) = 0
    !
  end do node_store_loop
 !call debug_timer(stop_timer,"storing the edges")
  !
  ! Save the number of unique edges in ledge
  !
  ledge = unique_edges
 !write (iout,2) mypnum,"ledge = ",ledge
  !
  ! Reallocate nodes_on_edge to reduce its size so that it is exactly
  ! allocated to the number of edges in the grid
  !
  call reallocate(nodes_on_edge,2,ledge)
  !
  ! Deallocate list_of_other_nodes since it is no longer needed
  !
  deallocate ( list_of_other_nodes , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"list_of_other_nodes",2,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  !#########################################################################
  !#########################################################################
  !#########################################################################
  !#####  3. Use the newly created nodes_on_edge array to go back and  #####
  !#####     create the list of edges attached to each grid node.      #####
  !#########################################################################
  !#########################################################################
  !#########################################################################
  !
  if (present(edges_with_node_ptr) .and. present(edges_with_node)) then
    !
    ! Allocate the edges_with_node_ptr array
    !
    allocate ( edges_with_node_ptr(1:lnode+1) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edges_with_node_ptr",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    ! Create the edges_with_node_ptr array
    !
    do this_edge = 1,ledge
      do i = 1,2
        this_node = nodes_on_edge(i,this_edge) + 1
        edges_with_node_ptr(this_node) = edges_with_node_ptr(this_node) + 1
      end do
    end do
    !
    ! Reshuffle edges_with_node_ptr so that its cumulative
    !
    do this_node = 2,size(edges_with_node_ptr)
      edges_with_node_ptr(this_node) = edges_with_node_ptr(this_node) + &
                                       edges_with_node_ptr(this_node-1)
    end do
    !
    ! Allocate the edges_with_node array
    !
    allocate ( edges_with_node(1:last(edges_with_node_ptr)) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edges_with_node",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    ! Add each edge to the list of attached edges for the two nodes
    ! defining the edge
    !
    do this_edge = 1,ledge
      do i = 1,2
        this_node = nodes_on_edge(i,this_edge)
        next_edge_idx = edges_with_node_ptr(this_node) + 1
        edges_with_node(next_edge_idx) = this_edge
        edges_with_node_ptr(this_node) = next_edge_idx
      end do
    end do
    !
    ! The previous loop has shifted edges_with_node_ptr left 1 index so we
    ! need to shift it back to the right 1 index
    !
    edges_with_node_ptr = eoshift( edges_with_node_ptr , shift=-1 , boundary=0 )
    !
    if (last(edges_with_node_ptr) /= size(edges_with_node)) then
      write (error_message,1)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
    !
    write (iout,2) mypnum,"ledge",ledge
    write (iout,2) mypnum,"size(edges_with_node)",size(edges_with_node)
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("ERROR! The number of edges found with nodes is inconsistent!")
  2 format (" CPU-",i0," : ",a," = ",i0)
  3 format ("The parameter max_lon needs to be greater than ",i0)
  !
end subroutine find_nodes_on_each_edge
!
!###############################################################################
!
subroutine find_cells_surrounding_each_edge(nodes_on_edge, &
                                   nodes_of_cell,nodes_of_cell_ptr, &
                                   cells_with_node,cells_with_node_ptr, &
                                   cells_with_edge,cells_with_edge_ptr)
  !
 !use ifport, only : qsort,sizeof_size_t
  !
  !.. Formal Arguments ..
  integer,              dimension(:,:),  intent(in) :: nodes_on_edge
  integer,              dimension(:),    intent(in) :: nodes_of_cell
  integer,              dimension(:),    intent(in) :: nodes_of_cell_ptr
  integer,              dimension(:),    intent(in) :: cells_with_node
  integer,              dimension(:),    intent(in) :: cells_with_node_ptr
  integer, allocatable, dimension(:), intent(inout) :: cells_with_edge
  integer, allocatable, dimension(:), intent(inout) :: cells_with_edge_ptr
  !
  !.. Local Scalars ..
  integer :: i,ibeg,iend,ierr,n1,n2,ledge
  integer :: this_node,this_cell,this_edge
  integer :: other_node,next_cell_idx
 !integer, allocatable :: sorted_cwe(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_cells_surrounding_each_edge"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ledge = size(nodes_on_edge,dim=2)
  !
  !###################################################################
  !###################################################################
  !###################################################################
  !#####  1. Create the pointer array needed to create the list  #####
  !#####     of cells that contain each edge.                    #####
  !###################################################################
  !###################################################################
  !###################################################################
  !
  allocate ( cells_with_edge_ptr(1:ledge+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_with_edge_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! Loop through all the edges and count the number of cells that
  ! contain each edge
  !
  do this_edge = 1,ledge
    !
    ! Initialize cells_with_edge_ptr(this_edge+1) to the value in
    ! cells_with_edge_ptr(this_edge) to make the array cumulative
    !
    cells_with_edge_ptr(this_edge+1) = cells_with_edge_ptr(this_edge)
    !
    ! Get the two nodes defining this edge
    !
    this_node  = nodes_on_edge(1,this_edge)
    other_node = nodes_on_edge(2,this_edge)
    !
    ! Beginning and ending indices for the cells surrounding the current node
    !
    ibeg = cells_with_node_ptr(this_node)+1
    iend = cells_with_node_ptr(this_node+1)
    !
    ! Loop over the cells that surround this_node and find those cells that
    ! also contain other_node. As long as there is a 1-to-1 matching of all
    ! grid faces (e.g., C(0) continuous, and triangle and quad faces only
    ! match faces of the same geometry), any cell containing these two nodes
    ! also contains this edge.
    !
    do i = ibeg,iend
      !
      ! NOTE: We need the absolute value of this because the cell index could
      !       be negative if it is for a boundary face ghost cell and the
      !       optional argument make_bnd_faces_negative was present with the
      !       value .TRUE. when subroutine find_cells_surrounding_each_node
      !       was called
      !
      this_cell = abs(cells_with_node(i))
      !
      n1 = nodes_of_cell_ptr(this_cell)+1
      n2 = nodes_of_cell_ptr(this_cell+1)
      !
      ! If other_node matches any of the nodes in nodes_of_cell(n1:n2),
      ! increase the cell counter for this edge by one so that later we
      ! can add this cell to the list of cells containing the current edge
      !
      if (any(other_node == nodes_of_cell(n1:n2))) then
        cells_with_edge_ptr(this_edge+1) = cells_with_edge_ptr(this_edge+1) + 1
      end if
      !
    end do
    !
  end do
  !
  !#####################################################################
  !#####################################################################
  !#####################################################################
  !#####  2. Use the previously created integer pointer array      #####
  !#####     cells_with_edge_ptr to store the list of cells        #####
  !#####     containing each edge in the array cells_with_edge.    #####
  !#####################################################################
  !#####################################################################
  !#####################################################################
  !
  allocate ( cells_with_edge(1:last(cells_with_edge_ptr)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_with_edge",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! Loop back through the edges, but this time add the found cells
  ! to the list of cells containing each edge
  !
  do this_edge = 1,ledge
    !
    ! Get the two nodes defining this edge
    !
    this_node  = nodes_on_edge(1,this_edge)
    other_node = nodes_on_edge(2,this_edge)
    !
    ! Beginning and ending indices for the cells surrounding the current node
    !
    ibeg = cells_with_node_ptr(this_node)+1
    iend = cells_with_node_ptr(this_node+1)
    !
    ! Loop over the cells that surround this_node and find those cells that
    ! also contain other_node. As long as there is a 1-to-1 matching of all
    ! grid faces (e.g., C(0) continuous, and triangle and quad faces only
    ! match faces of the same geometry), any cell containing these two nodes
    ! also contains this edge.
    !
    do i = ibeg,iend
      !
      ! NOTE: We need the absolute value of this because the cell index could
      !       be negative if it is for a boundary face ghost cell and the
      !       optional argument make_bnd_faces_negative was present with the
      !       value .TRUE. when subroutine find_cells_surrounding_each_node
      !       was called
      !
      this_cell = abs(cells_with_node(i))
      !
      n1 = nodes_of_cell_ptr(this_cell)+1
      n2 = nodes_of_cell_ptr(this_cell+1)
      !
      ! If other_node matches any of the nodes in nodes_of_cell(n1:n2),
      ! add this cell to the list of cells containing the current edge
      ! NOTE: We need to use the actual value of cells_with_node(i) instead
      !       of this_cell because we want to keep the cell index negative if
      !       this is a boundary face ghost cell
      !
      if (any(other_node == nodes_of_cell(n1:n2))) then
        next_cell_idx = cells_with_edge_ptr(this_edge) + 1
        cells_with_edge(next_cell_idx) = cells_with_node(i)
        cells_with_edge_ptr(this_edge) = next_cell_idx
      end if
      !
    end do
    !
  end do
  !
  ! The previous loop has shifted cells_with_edge_ptr left 1 index so we
  ! need to shift it back to the right 1 index
  !
  cells_with_edge_ptr = eoshift( cells_with_edge_ptr , shift=-1 , boundary=0 )
  !
  if (last(cells_with_edge_ptr) /= size(cells_with_edge)) then
    write (error_message,1)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
 !allocate ( sorted_cwe(1:size(cells_with_edge)) , source=cells_with_edge , &
 !           stat=ierr , errmsg=error_message )
 !call alloc_error(pname,"sorted_cwe",1,__LINE__,__FILE__,ierr,error_message)
 !!
 !if (mypnum == 0) then
 !  call qsort(sorted_cwe,size(sorted_cwe,kind=sizeof_size_t), &
 !             int(wip,kind=sizeof_size_t),cmp_function)
 !  do i = 1,size(sorted_cwe)
 !    write (mypnum+130,'(i9)') sorted_cwe(i)
 !  end do
 !end if
 !call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("ERROR! The number of cells found with edges is inconsistent!")
  !
 !contains
 !integer(2) function cmp_function(a1, a2)
 !integer :: a1, a2
 !cmp_function=a1-a2
 !end function
end subroutine find_cells_surrounding_each_edge
!
!###############################################################################
!
subroutine create_edge_derived_type(cells_with_edge,cells_with_edge_ptr, &
                                    nodes_on_edge)
  !
  !.. Use Statements ..
  use geovar, only : edge,cell
  !
  !.. Formal Arguments ..
  integer, dimension(:),   intent(in) :: cells_with_edge
  integer, dimension(:),   intent(in) :: cells_with_edge_ptr
  integer, dimension(:,:), intent(in) :: nodes_on_edge
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,nc,ne,ierr
  integer :: ledge,lcell,this_order
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_edge_derived_type"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ledge = size(nodes_on_edge,dim=2)
  lcell = size(cell)
  !
  ! Allocate the edge array
  !
  allocate ( edge(1:ledge) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"edge",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Loop through the edges and save the required information in the edge array
  !
  do ne = 1,ledge
    !
    n1 = cells_with_edge_ptr(ne)+1
    n2 = cells_with_edge_ptr(ne+1)
    !
    ! Initialize the order for this edge to a large negative number
    !
    this_order = -huge( int(0,kind=kind(this_order)) )
    !
    ! Find the maximum order of all cells containing this edge
    !
    do n = n1,n2
      !
      ! The index for a cell containing this edge
      ! NOTE: We want to ignore all boundary faces within cells_with_edge
      !       since the order of the boundary face is determined by its host
      !       cell, and the host cell will already be included in the list
      !       of cells containing this edge thats stored in cells_with_edge.
      ! NOTE: We need the absolute value of this because the cell index could
      !       be negative if it is for a boundary face ghost cell and the
      !       optional argument make_bnd_faces_negative was present with the
      !       value .TRUE. when subroutine find_cells_surrounding_each_node
      !       was called
      !
      nc = abs(cells_with_edge(n))
      !
      if (nc <= lcell) then
        this_order = max( this_order , cell(nc)%order )
      end if
      !
    end do
    !
    ! Save this information in the edge array
    !
    edge(ne)%nodes(1) = nodes_on_edge(1,ne)
    edge(ne)%nodes(2) = nodes_on_edge(2,ne)
    edge(ne)%order    = this_order
    !
  end do
  !
  ! Go back through the cell derived type and save the order of the
  ! the edges into the derived type component cell(:)%edge(:)%order
  !
  do nc = 1,lcell
    do ne = 1,size(cell(nc)%edge)
      cell(nc)%edge(ne)%order = edge( cell(nc)%edge(ne)%idx )%order
    end do
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine create_edge_derived_type
!
!###############################################################################
!
include "Functions/last.f90"
!
!###############################################################################
!
end module connectivity_mod
