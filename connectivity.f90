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
  !
contains
!
!###############################################################################
!
subroutine get_grid_connectivity()
  !
  !.. Use Statements ..
  use geovar, only : bface
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:)   :: cells_with_node
  integer, allocatable, dimension(:)   :: cells_with_node_ptr
  integer, allocatable, dimension(:)   :: cells_surr_cell
  integer, allocatable, dimension(:)   :: cells_surr_cell_ptr
  !
continue
  !
  ! Get the cells surrounding each node
  !
  call find_cells_surrounding_each_node(cells_with_node,cells_with_node_ptr)
  !
  ! Get the cells surrounding each cell
  !
  call find_cells_surrounding_each_cell(cells_with_node,cells_with_node_ptr, &
                                        cells_surr_cell,cells_surr_cell_ptr)
  !
  ! Get the cells on each face
  !
  call find_cells_on_each_face(cells_surr_cell,cells_surr_cell_ptr)
  !
  ! Get the indices for the flux points for each cell
  !
  call find_flxpts_of_each_cell
  !
  ! Get the unit cell subconnectivities for use when writing solution file
  !
  call get_unit_cell_subconnectivity
  !
end subroutine get_grid_connectivity
!
!###############################################################################
!
subroutine find_cells_surrounding_each_node(cells_with_node,cells_with_node_ptr)
  !
  ! Find the cells that surround each node
  !
  ! For node i_node,
  ! cells_with_node(ibeg:iend) gives all of the grid cells that contain i_node
  ! where : ibeg = cells_with_node_ptr(i_node)+1
  !         iend = cells_with_node_ptr(i_node+1)
  !
  !.. Use Statements ..
  use geovar, only : nnode,ncell,n_totcel,nfbnd,bface
  use geovar, only : nodes_of_cell_ptr,nodes_of_cell
  !
  !.. Formal Arguments ..
  integer, allocatable, dimension(:), intent(inout) :: cells_with_node
  integer, allocatable, dimension(:), intent(inout) :: cells_with_node_ptr
  !
  !.. Local Scalars ..
  integer :: i,ip,is,ibeg,iend,n,nc,ierr
  !
  !.. Local Allocatable Arrays ..
  logical, allocatable, dimension(:,:) :: ghost_msk
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_cells_surrounding_each_node"
  !
continue
  !
  allocate ( cells_with_node_ptr(1:nnode+1) , source=0 , stat=ierr )
  call alloc_error(pname,"cells_with_node_ptr",1,ierr)
  !
  ! Loop over the cells, storing 'ahead'
  !
  do n = 1,ncell
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
  ! Add in the ghost cells. Make sure to only use the nodes
  ! of the ghost cell which lie on the boundary face.
  ! NOTE: We need a mask array here for the ghost cells to
  !       prevent this loop from finding a particular ghost
  !       cell more than once on an individual node.
  !
  allocate ( ghost_msk(ncell+1:n_totcel,1:nnode) , source=fals , stat=ierr )
  call alloc_error(pname,"ghost_msk",1,ierr)
  !
  do n = 1,nfbnd
    !
    nc = bface(10,n)
    !
    ibeg = nodes_of_cell_ptr(nc)+1
    iend = nodes_of_cell_ptr(nc+1)
    !
    do i = ibeg,iend
      ip = nodes_of_cell(i)
      if (ip <= nnode) then
        if (.not. ghost_msk(nc,ip)) then
          cells_with_node_ptr(ip+1) = cells_with_node_ptr(ip+1) + 1
          ghost_msk(nc,ip) = true
        end if
      end if
    end do
    !
  end do
  !
  ! Reshuffle cells_with_node_ptr
  !
  do n = 2,size(cells_with_node_ptr)
    cells_with_node_ptr(n) = cells_with_node_ptr(n) + cells_with_node_ptr(n-1)
  end do
  !
  allocate ( cells_with_node(1:last(cells_with_node_ptr)) , source=0 , stat=ierr )
  call alloc_error(pname,"cells_with_node",1,ierr)
  !
  ! Now store the cells surrounding each node in cells_with_node
  !
  do n = 1,ncell
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
  ! Add in the ghost cells. Make sure to only use the nodes
  ! of the ghost cell which lie on the boundary face. 
  ! NOTE: We need a mask array here for the ghost cells to
  !       prevent this loop from finding a particular ghost
  !       cell more than once on an individual node.
  !
  ghost_msk(:,:) = fals
  !
  do n = 1,nfbnd
    !
    nc = bface(10,n)
    !
    ibeg = nodes_of_cell_ptr(nc)+1
    iend = nodes_of_cell_ptr(nc+1)
    !
    do i = ibeg,iend
      ip = nodes_of_cell(i)
      if (ip <= nnode) then
        if (.not. ghost_msk(nc,ip)) then
          is = cells_with_node_ptr(ip) + 1
          cells_with_node_ptr(ip) = is
          cells_with_node(is) = nc
          ghost_msk(nc,ip) = true
        end if
      end if
    end do
    !
  end do
  !
  ! We are finished with ghost_msk so deallocate
  !
  deallocate ( ghost_msk , stat=ierr )
  call alloc_error(pname,"ghost_msk",2,ierr)
  !
  ! Finally, reorder cells_with_node_ptr
  !
  do n = size(cells_with_node_ptr),2,-1
    cells_with_node_ptr(n) = cells_with_node_ptr(n-1)
  end do
  cells_with_node_ptr(1) = 0
  !
end subroutine find_cells_surrounding_each_node
!
!###############################################################################
!
subroutine find_cells_surrounding_each_cell(cells_with_node,cells_with_node_ptr, &
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
  !.. Use Statements ..
  use geovar, only : cell_geom,n_totnod,ncell
  use geovar, only : nodes_of_cell_ptr,nodes_of_cell
  !
  !.. Formal Arguments ..
  integer, allocatable, dimension(:), intent(inout) :: cells_with_node
  integer, allocatable, dimension(:), intent(inout) :: cells_with_node_ptr
  integer, allocatable, dimension(:), intent(inout) :: cells_surr_cell
  integer, allocatable, dimension(:), intent(inout) :: cells_surr_cell_ptr
  !
  !.. Local Scalars ..
  integer :: i,i1,i2,ieadj,csc_ptr
  integer :: j,j1,j2,n,n1,n2,ierr
  integer :: this_geom,this_side,nfnods
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
  allocate ( lpoin(1:n_totnod) , source=0 , stat=ierr )
  call alloc_error(pname,"lpoin",1,ierr)
  allocate ( cells_surr_cell_ptr(1:ncell+1) , source=0 , stat=ierr )
  call alloc_error(pname,"cells_surr_cell_ptr",1,ierr)
  !
  ! Loop to set up cells_surr_cell_ptr
  !
  do n = 1,ncell
    cells_surr_cell_ptr(n+1) = cells_surr_cell_ptr(n) + geom_faces(cell_geom(n))
  end do
  !
  allocate ( cells_surr_cell(1:last(cells_surr_cell_ptr)) , source=0 , stat=ierr )
  call alloc_error(pname,"cells_surr_cell",1,ierr)
  !
  ! Loop to find cells_surr_cell
  !
  do n = 1,ncell
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
      do i = i1,i2
        !
        ! Current cell containing this node
        !
        j = cells_with_node(i)
        !
        ! Skip to the next cell if this is the current cell we are working on
        !
        if (j == n) cycle
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
        !
        if (sum(lpoin(nodes_of_cell(j1:j2))) == nfnods) then
          ieadj = j
          exit
        end if
        !
      end do
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
  deallocate ( lpoin , stat=ierr )
  call alloc_error(pname,"lpoin",2,ierr)
  !
  ! We are also finished with cells_with_node
  ! and cells_with_node_ptr so deallocate
  !
  deallocate ( cells_with_node , stat=ierr )
  call alloc_error(pname,"cells_with_node",2,ierr)
  deallocate ( cells_with_node_ptr , stat=ierr )
  call alloc_error(pname,"cells_with_node_ptr",2,ierr)
  !
end subroutine find_cells_surrounding_each_cell
!
!###############################################################################
!
subroutine find_cells_on_each_face(cells_surr_cell,cells_surr_cell_ptr)
  !
  ! For face i_face, ...
  ! nodes_on_face(1:2,i_face) : grid nodes defining i_face
  !
  ! For face i_face, with the face vector going from
  ! nodes_on_face(1,i_face) to nodes_on_face(2,i_face), ...
  ! cells_on_face(1,i_face) : attached cell to the left  of i_face
  ! cells_on_face(2,i_face) : attached cell to the right of i_face
  ! cells_on_face(3,i_face) : node of left  triangular cell not on i_face (quads = 0)
  ! cells_on_face(4,i_face) : node of right triangular cell not on i_face (quads = 0)
  ! cells_on_face(5,i_face) : index of correction function array needed for left  cell
  ! cells_on_face(6,i_face) : index of correction function array needed for right cell
  ! cells_on_face(7,i_face) : geometry of i_face
  ! cells_on_face(8,i_face) : order of i_face
  !
  !
  !.. Use Statements ..
  use order_mod, only : geom_solpts
  use geovar,    only : nr,nfbnd,nface,ncell,n_totnod
  use geovar,    only : xyz_nodes,bface,cell_geom,cell_order
  use geovar,    only : nodes_of_cell_ptr,nodes_of_cell
  use geovar,    only : cells_on_face,nodes_on_face,faces_of_cell
  !
  !.. Formal Arguments ..
  integer, allocatable, dimension(:),   intent(inout) :: cells_surr_cell
  integer, allocatable, dimension(:),   intent(inout) :: cells_surr_cell_ptr
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,nf,nc,ierr
  integer :: i,ic,ifa,ip1,ip2,ip3,jc
  integer :: face_counter
  integer :: left_cell,left_geom
  integer :: host_cell,ghost_cell,host_geom
  integer :: host_side,nfnods,ncnods
  integer :: right_cell,right_geom,ncf
  integer :: left_side,right_side
  integer :: face_geom,face_order,face_pts
  integer :: c1,c2,this_cell,this_geom,adj_cell
  integer :: left_order,right_order
  !
  real(wp) :: dir
  !
  !.. Local Arrays ..
  integer, dimension(1:4) :: lnode
  integer, dimension(1:4) :: ip
  real(wp), dimension(1:nr) :: rnrm
  real(wp), dimension(1:nr) :: face_center
  real(wp), dimension(1:nr) :: cell_center
  real(wp), dimension(1:3,1:4) :: xyz_face
  real(wp), dimension(1:3,1:8) :: xyz_cell
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:) :: lpoin
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_cells_on_each_face"
  !
continue
  !
  ! ###########
  ! ###########
  ! ###  1  ###
  ! ###########
  ! ###########
  !
  ! Find the total number of faces
  ! and then allocate cells_on_face
  !
  nface = nfbnd
  do n = 1,ncell
    !
    n1 = cells_surr_cell_ptr(n)+1
    n2 = cells_surr_cell_ptr(n+1)
    !
    nface = nface + count( cells_surr_cell(n1:n2)>n .and. &
                           cells_surr_cell(n1:n2)<=ncell )
    !
  end do
  !
  allocate ( cells_on_face(1:8,1:nface) , source=0 , stat=ierr )
  call alloc_error(pname,"cells_on_face",1,ierr)
  !
  allocate ( nodes_on_face(1:4,1:nface) , source=0 , stat=ierr )
  call alloc_error(pname,"nodes_on_face",1,ierr)
  !
  allocate ( faces_of_cell(1:6,1:ncell) , source=-99 , stat=ierr )
  call alloc_error(pname,"faces_of_cell",1,ierr)
  !
  ! ###########
  ! ###########
  ! ###  2  ###
  ! ###########
  ! ###########
  !
  ! Get the cells on each boundary face
  !
  do nf = 1,nfbnd
    !
    host_cell = bface(2,nf)   ! Host cell on boundary face
    ghost_cell = bface(10,nf) ! Ghost cell on boundary face
    !
    ! Store the cells on this boundary face
    !
    cells_on_face(1,nf) = host_cell
    cells_on_face(2,nf) = ghost_cell
    !
    host_geom = cell_geom(host_cell)             ! geometry of host cell
    host_side = bface(4,nf)                      ! face of host cell on boundary
    nfnods = num_face_nodes(host_side,host_geom) ! # of nodes defining bnd face
    !
    n1 = nodes_of_cell_ptr(host_cell)+1 ! beginning index of nodes defining host cell
    n2 = nodes_of_cell_ptr(host_cell+1) ! ending    index of nodes defining host cell
    !
    ! Node offsets for this face
    ip(1:nfnods) = side_nodes(1:nfnods,host_side,host_geom)
    !
    ! Store the nodes on this boundary face
    !
    nodes_on_face(1:nfnods,nf) = nodes_of_cell( n1 + ip(1:nfnods) )
    !
    ! Store the geometry of this face
    !
    cells_on_face(7,nf) = geom_of_face(nfnods)
    !
    ! Store the side of the left cell that is on this face since we have this
    ! information available right now. Store the side of the right/ghost cell
    ! as well since this will always be 1.
    !
    cells_on_face(3,nf) = host_side
    cells_on_face(4,nf) = 1
    !
    ! Store this face for the host cell
    !
    faces_of_cell(host_side,host_cell) = nf
    !
  end do
  !
  ! ###########
  ! ###########
  ! ###  3  ###
  ! ###########
  ! ###########
  !
  ! Get the cells for the interior faces
  !
  face_counter = nfbnd
  !
  do this_cell = 1,ncell
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
      if (adj_cell>this_cell .and. adj_cell<=ncell) then
        !
        face_counter = face_counter + 1
        !
        nf = nc-c1+1
        !
        ! Store the cells on this interior face
        !
        cells_on_face(1,face_counter) = this_cell
        cells_on_face(2,face_counter) = adj_cell
        !
        nfnods = num_face_nodes(nf,this_geom)             ! # of nodes defining this face
        ip(1:nfnods) = side_nodes(1:nfnods,nf,this_geom) ! node offsets for this face
        !
        ! Store the nodes on this interior face
        !
        nodes_on_face(1:nfnods,face_counter) = nodes_of_cell( n1 + ip(1:nfnods) )
        !
        ! Store the geometry of this face
        !
        cells_on_face(7,face_counter) = geom_of_face(nfnods)
        !
        ! Store the side of the left cell that is on this face since we have
        ! this information available right now. We will figure out the side of
        ! the right cell in section 5.
        !
        cells_on_face(3,face_counter) = nf
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
  ! We are finished with cells_surr_cell
  ! and cells_surr_cell_ptr so deallocate
  !
  deallocate ( cells_surr_cell , stat=ierr )
  call alloc_error(pname,"cells_surr_cell",2,ierr)
  !
  deallocate ( cells_surr_cell_ptr , stat=ierr )
  call alloc_error(pname,"cells_surr_cell_ptr",2,ierr)
  !
  ! ###########
  ! ###########
  ! ###  4  ###
  ! ###########
  ! ###########
  !
  ! For the face vector going from nodes_on_face(1,nf) to nodes_on_face(2,nf),
  ! this is a test to make sure that cells_on_face(1,nf) is to the left of
  ! this face vector and cells_on_face(2,nf) is on the right
  !
  ic = 0
  do nf = nfbnd+1,nface
    !
    left_cell = cells_on_face(1,nf)  ! left cell on this face
    face_geom = cells_on_face(7,nf)  ! geometry of this face
    left_geom = cell_geom(left_cell) ! geometry of the left cell
    !
    ! # of nodes defining the face and defining the left cell
    !
    nfnods = geom_nodes(face_geom)
    ncnods = geom_nodes(left_geom)
    !
    ! Coordinates of the nodes on this face
    !
    xyz_face(1:nr,1:nfnods) = xyz_nodes(1:nr,nodes_on_face(1:nfnods,nf))
    !
    ! Coordinates of the center of this face
    !
    face_center(1:nr) = sum(xyz_face(1:nr,1:nfnods),dim=2) / real(nfnods,kind=wp)
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
    cell_center(1:nr) = sum(xyz_cell(1:nr,1:ncnods),dim=2) / real(ncnods,kind=wp)
    !
    ! Get the normal to this face using the nodes that define it
    !
    rnrm = face_normal( xyz_face(1:nr,1:nfnods) )
    !
    ! Compute the dot product of the vector from the center of
    ! the left cell to the center of the face with the face normal
    !
    dir = dot( rnrm(1:nr) , face_center(1:nr) - cell_center(1:nr) )
    !
    ! If dot is negative, then the cell stored in cells_on_face(2,nf)
    ! is actually the cell to the left of the face vector and needs
    ! to be stored in cells_on_face(1,nf)
    !
    if (dir < zero) then
      write (*,1) nf,(nodes_on_face(i,nf),i=1,nfnods)
      ic = 1
    end if
    !
  end do
  if (ic == 1) call stop_recon(stop_mpi,pname)
  !
  ! ###########
  ! ###########
  ! ###  5  ###
  ! ###########
  ! ###########
  !
  allocate ( lpoin(1:n_totnod) , source=0 , stat=ierr )
  call alloc_error(pname,"lpoin",1,ierr)
  !
  ! Find the side of the right cell that is on each interior face
  !
 !do n = 1,size(faces_of_cell,dim=2)
 !  write (102,2) n,(faces_of_cell(i,n),i=1,size(faces_of_cell,dim=1))
 !end do
  !
  loop_right_cell_face: do nf = nfbnd+1,nface
    !
    nfnods = geom_nodes(cells_on_face(7,nf)) ! # of nodes defining this face
    !
    ! Mark the nodes on the current face
    !
    lpoin( nodes_on_face(1:nfnods,nf) ) = 1
    !
    ! Information for the right cell on the face
    !
    right_cell = cells_on_face(2,nf)
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
        cells_on_face(4,nf) = ncf
        faces_of_cell(ncf,right_cell) = nf
        lpoin( ip(1:nfnods) ) = 0
        cycle loop_right_cell_face
      end if
      !
    end do loop_cell_faces
    !
   !write (103,11) nf,(nodes_on_face(n,nf),n=1,nfnods),right_cell, &
   !                (nodes_of_cell(n),n=n1,n2)
   !11 format (" Error: nf=",i0,"; fn=",4(1x,i0),"; rc=",i0,"; cn=",8(1x,i0))
  end do loop_right_cell_face
  !
  deallocate ( lpoin , stat=ierr )
  call alloc_error(pname,"lpoin",2,ierr)
  !
  ! ###########
  ! ###########
  ! ###  6  ###
  ! ###########
  ! ###########
  !
  ! Get the order for each face based on the max order of each cell on the face
  !
  do nf = 1,nface
    !
    left_cell = cells_on_face(1,nf)
    right_cell = cells_on_face(2,nf)
    !
    left_order = cell_order(left_cell)
    right_order = cell_order(right_cell)
    !
    cells_on_face(8,nf) = max(left_order,right_order)
    !
  end do
  !
  ! ###########
  ! ###########
  ! ###  7  ###
  ! ###########
  ! ###########
  !
  ! Get the starting index used in the correction procedure for both sides of
  ! the faces
  !
 !do n = 1,size(faces_of_cell,dim=2)
 !  if (any(faces_of_cell(:,n)==-99)) then
 !    write (101,2) n,(faces_of_cell(i,n),i=1,size(faces_of_cell,dim=1))
 !  end if
 !end do
 !2 format (" Cell ",i4,": ",100i7)
  !
  do nf = 1,nface
    !
    left_cell = cells_on_face(1,nf)
    right_cell = cells_on_face(2,nf)
    !
    left_side = cells_on_face(3,nf)
    right_side = cells_on_face(4,nf)
    !
    left_order = cell_order(left_cell)
    right_order = cell_order(right_cell)
    !
    ! For the left cell, the starting index is the sum of all face points
    ! for all faces up to the current face + 1. The '+1' is why face_pts
    ! is initialized to 1.
    !
    face_pts = 1
    do ncf = 1,left_side-1
      !
      ifa = faces_of_cell(ncf,left_cell)
      !
      face_geom = cells_on_face(7,ifa)
      face_pts = face_pts + geom_solpts(face_geom,left_order)
      !
    end do
    !
    cells_on_face(5,nf) = face_pts
    !
    ! For the right cell, it depends on whether
    ! this is a boundary or interior face
    !
    if (nf <= nfbnd) then
      !
      ! This is a boundary face so the starting index is just
      ! the number of face points on the boundary face
      !
      face_geom = cells_on_face(7,ifa)
      face_pts = geom_solpts(face_geom,right_order)
      !
    else
      !
      ! This is an interior face so the starting index is the sum of all
      ! face points for all faces up to and including the current face.
      !
      face_pts = 0
      do ncf = 1,right_side
        !
        ifa = faces_of_cell(ncf,right_cell)
        !
        face_geom = cells_on_face(7,ifa)
        face_pts = face_pts + geom_solpts(face_geom,right_order)
        !
      end do
      !
    end if
    !
    cells_on_face(6,nf) = face_pts
    !
   !face_geom = cells_on_face(7,nf)
   !face_order = cells_on_face(8,nf)
   !!
   !face_pts = geom_solpts(face_geom,face_order)
   !!
   !! NOTE: THIS ASSUMES THAT THE NUMBER OF FACE POINTS IS THE SAME FOR ALL
   !!       FACES OF THE LEFT AND RIGHT CELLS (I.E., 2D, TETS, AND HEXES AND THE
   !!       CELL ORDER IS CONSTANT THROUGHOUT). THIS IS NOT THE CASE FOR PYRAMIDS
   !!       AND PRISMS OR IF THE CELL ORDER IS NOT THE SAME FOR ALL CELLS.
   !!
   !cells_on_face(5,nf) = face_pts*(left_side-1) + 1
   !cells_on_face(6,nf) = face_pts*right_side
    !
  end do
  !
 !write (110,*)
 !write (110,*) "Boundary Faces"
 !write (110,*)
 !do nf = 1,nface
 !  write (110,1000) nf, cells_on_face(1:5:2,nf),cells_on_face(2:6:2,nf)
 !  if (nf == nfbnd) then
 !    write (110,*)
 !    write (110,*) "Interior Faces"
 !    write (110,*)
 !  end if
 !end do
 !1000 format ("Face ",i4,": ",2(6x,i4,1x,i1,1x,i2))
  !
  ! Format Statements
  !
  1 format (" Error on face - ",i0,",  nodes ",i0,3(:,2x,i0))
  !
end subroutine find_cells_on_each_face
!
!###############################################################################
!
subroutine find_flxpts_of_each_cell()
  !
  !.. Use Statements ..
  use geovar, only : nfbnd,ncell,n_totcel,nface,nflxprm
  use geovar, only : flxpts_of_cell_ptr,flxpts_of_cell,bface
  use geovar, only : cells_on_face,nodes_of_cell_ptr
  use geovar, only : cell_geom,cell_order
  use order_mod, only : geom_solpts
  !
  !.. Local Scalars ..
  integer :: n,nf,ierr,k
  integer :: jl,lf1,lf2,lk1,lk,lfp,nlfp
  integer :: jr,rf1,rf2,rk1,rk,rfp,nrfp
  integer :: this_order,this_geom
  integer :: face_nodes,face_geom
  integer :: n_flxpts
  integer :: left_cell,left_order,nfp,right_cell,right_order
  !
  !.. Local Parameters ..
  integer, parameter :: left_side = 1
  integer, parameter :: right_side = 2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_flxpts_of_each_cell"
  !
continue
  !
  ! Create the pointer for the flxpts_of_cell array
  !
  allocate ( flxpts_of_cell_ptr(1:ncell+1) , source=0 , stat=ierr )
  call alloc_error(pname,"flxpts_of_cell_ptr",1,ierr)
  !
  ! Set up the flux points pointer array for each interior cell
  ! NOTE: This array is only used for interior cells so it does
  !       not need to account for ghost cells
  !
  do n = 1,ncell
    !
    this_order = cell_order(n)
    this_geom = cell_geom(n)
    !
    n_flxpts = 0
    do nf = 1,geom_faces(this_geom)
      !
      face_nodes = num_face_nodes(nf,this_geom)
      face_geom = geom_of_face(face_nodes)
      !
      n_flxpts = n_flxpts + geom_solpts(face_geom,this_order)
      !
    end do
    !
    flxpts_of_cell_ptr(n+1) = flxpts_of_cell_ptr(n) + n_flxpts
    !
  end do
  !
  ! Allocate the flxpts_of_cell array
  !
  allocate ( flxpts_of_cell(1:nflxprm,1:last(flxpts_of_cell_ptr)) , &
             source=0 , stat=ierr )
  call alloc_error(pname,"flxpts_of_cell",1,ierr)
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
  !  NOTE: CURRENTLY, EVERYTHING IN THIS ROUTINE HAS FLXPTS_OF_CELL GIVING THE
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
  ! Loop through the boundary faces and fill
  ! the left part of the flxpts_of_cell array
  !
  do nf = 1,nfbnd
    !
    left_cell = cells_on_face(1,nf) ! left  cell on face
    !
    lf1 = flxpts_of_cell_ptr(left_cell)+1 ! beginning flxpt index for left  cell
    lf2 = flxpts_of_cell_ptr(left_cell+1) ! ending    flxpt index for left  cell
    !
    face_geom = cells_on_face(7,nf) ! geometry of this face
    left_order = cell_order(left_cell) ! order of this face
    !
    nfp = geom_solpts(face_geom,left_order)
    !
    lk1 = cells_on_face(5,nf) ! beginning flxpt on current face for left  cell
    !
    ! Loop through each flux point on the current face and get
    ! the flux point information for both cells on this face
    !
    do k = 1,nfp
      !
      lk = lk1-1 + k ! idx of current flxpt for left  cell (from 1 to lf2-lf1+1)
      !
      lfp = lf1-1 + lk ! idx of current flxpt for left  cell (from lf1 to lf2)
      !
      flxpts_of_cell(1,lfp) = nf
      flxpts_of_cell(2,lfp) = k
      flxpts_of_cell(3,lfp) = lk
      flxpts_of_cell(4,lfp) = left_side
      !
    end do
    !
  end do
  !
  ! Loop through each interior face to fill the
  ! remaining parts of the flxpts_of_cell array
  !
  do nf = nfbnd+1,nface
    !
    ! Loop through each flux point on the current face and get
    ! the flux point information for both cells on this face
    !
    face_geom = cells_on_face(7,nf) ! geometry of this face
    !
    ! Left cell on face
    !
    left_cell = cells_on_face(1,nf) ! left  cell on face
    left_order = cell_order(left_cell) ! order of the left side of this face
    !
    lf1 = flxpts_of_cell_ptr(left_cell)+1 ! beginning flxpt index for left  cell
    lf2 = flxpts_of_cell_ptr(left_cell+1) ! ending    flxpt index for left  cell
    !
    nfp = geom_solpts(face_geom,left_order) ! # of flxpts on left side of this face
    !
    lk1 = cells_on_face(5,nf) ! beginning flxpt on current face for left  cell
    !
    do k = 1,nfp
      !
      lk = lk1-1 + k ! idx of current flxpt for left  cell (from 1 to lf2-lf1+1)
      lfp = lf1-1 + lk ! idx of current flxpt for left  cell (from lf1 to lf2)
      !
      flxpts_of_cell(1,lfp) = nf
      flxpts_of_cell(2,lfp) = k
      flxpts_of_cell(3,lfp) = lk
      flxpts_of_cell(4,lfp) = left_side
      !
    end do
    !
    ! Right cell on face
    !
    right_cell = cells_on_face(2,nf) ! right cell on face
    right_order = cell_order(right_cell) ! order of the right side of this face
    !
    rf1 = flxpts_of_cell_ptr(right_cell)+1 ! beginning flxpt index for right cell
    rf2 = flxpts_of_cell_ptr(right_cell+1) ! ending    flxpt index for right cell
    !
    nfp = geom_solpts(face_geom,right_order) ! # of flxpts on right side of this face
    !
    rk1 = cells_on_face(6,nf) ! beginning flxpt on current face for right cell
    !
    do k = 1,nfp
      !
      rk = rk1+1 - k ! idx of current flxpt for right cell (from 1 to rf2-rf1+1)
      rfp = rf1-1 + rk ! idx of current flxpt for right cell (from rf1 to rf2)
      !
      flxpts_of_cell(1,rfp) = nf
      flxpts_of_cell(2,rfp) = k
      flxpts_of_cell(3,rfp) = rk
      flxpts_of_cell(4,rfp) = right_side
      !
    end do
    !
  end do
  !
end subroutine find_flxpts_of_each_cell
!
!###############################################################################
!
subroutine get_unit_cell_subconnectivity()
  !
  !.. Use Statements ..
  use geovar, only : pst_geom
  use geovar, only : global_cell_order
  !
  !.. Local Scalars ..
  integer :: i,j,k,sk,ierr
  integer :: no1d,no2d,no3d,np1d
  integer :: p1,p2,p3,p4,p5,p6,p7,p8
  integer :: omin,omax,order
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
  ! Allocate the pst_geom array
  !
  omin = minval( global_cell_order )
  omax = maxval( global_cell_order )
  !
  allocate ( pst_geom(Geom_Min:Geom_Max,omin:omax) , stat=ierr )
  call alloc_error(pname,"pst_geom",1,ierr)
  !
  ! Loop through all the possible cell orders to get the connectivity
  ! for each of them
  !
  loop_over_orders: do order = omin,omax
    !
    no1d = max(1,order)
    no2d = no1d**2
    no3d = no1d**3
    np1d = no1d + 1
    !
    ! Allocate counter and counter3 for this order
    !
    allocate ( counter(1:no1d+1,1:no1d+1) , source=0 , stat=ierr )
    call alloc_error(pname,"counter",1,ierr)
    allocate ( counter3(1:no1d+1,1:no1d+1,1:no1d+1) , source=0 , stat=ierr )
    call alloc_error(pname,"counter3",1,ierr)
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
               source=0 , stat=ierr)
    !
    do i = 1,no1d
      pst_geom(Geom_Edge,order)%connectivity(1,i) = i
      pst_geom(Geom_Edge,order)%connectivity(2,i) = i + 1
    end do
    !
    ! Subtract one from the sub-connectivity so that it starts at zero
    !
    pst_geom(Geom_Edge,order)%connectivity = pst_geom(Geom_Edge,order)%connectivity - 1
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
               source=0 , stat=ierr)
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
    pst_geom(Geom_Tria,order)%connectivity = pst_geom(Geom_Tria,order)%connectivity - 1
    !
    ! ############################
    ! ## Quadrilateral Geometry ##
    ! ############################
    !
    allocate ( pst_geom(Geom_Quad,order)%connectivity(1:4,1:no2d) , &
               source=0 , stat=ierr)
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
    pst_geom(Geom_Quad,order)%connectivity = pst_geom(Geom_Quad,order)%connectivity - 1
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
               source=0 , stat=ierr)
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
    pst_geom(Geom_Tetr,order)%connectivity = pst_geom(Geom_Tetr,order)%connectivity - 1
    !
    ! ######################
    ! ## Pyramid Geometry ##
    ! ######################
    !
    allocate ( pst_geom(Geom_Pyra,order)%connectivity(1:8,1:no3d) , &
               source=0 , stat=ierr)
    !
    ! Subtract one from the sub-connectivity so that it starts at zero
    !
    pst_geom(Geom_Pyra,order)%connectivity = pst_geom(Geom_Pyra,order)%connectivity - 1
    !
    ! ####################
    ! ## Prism Geometry ##
    ! ####################
    !
    allocate ( pst_geom(Geom_Pris,order)%connectivity(1:8,1:no3d) , &
               source=0 , stat=ierr)
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
    pst_geom(Geom_Pris,order)%connectivity = pst_geom(Geom_Pris,order)%connectivity - 1
    !
    ! #########################
    ! ## Hexahedral Geometry ##
    ! #########################
    !
    allocate ( pst_geom(Geom_Hexa,order)%connectivity(1:8,1:no3d) , &
               source=0 , stat=ierr)
    !
    sk = 0
    do k = 1,no1d
      do j = 1,no1d
        do i = 1,no1d
          sk = sk + 1
          pst_geom(Geom_Hexa,order)%connectivity(1,sk) = np1d*np1d*(k-1) + np1d*(j-1) + i
          pst_geom(Geom_Hexa,order)%connectivity(2,sk) = np1d*np1d*(k-1) + np1d*(j-1) + i + 1
          pst_geom(Geom_Hexa,order)%connectivity(3,sk) = np1d*np1d*(k-1) + np1d*(j  ) + i + 1
          pst_geom(Geom_Hexa,order)%connectivity(4,sk) = np1d*np1d*(k-1) + np1d*(j  ) + i
          pst_geom(Geom_Hexa,order)%connectivity(5,sk) = np1d*np1d*(k  ) + np1d*(j-1) + i
          pst_geom(Geom_Hexa,order)%connectivity(6,sk) = np1d*np1d*(k  ) + np1d*(j-1) + i + 1
          pst_geom(Geom_Hexa,order)%connectivity(7,sk) = np1d*np1d*(k  ) + np1d*(j  ) + i + 1
          pst_geom(Geom_Hexa,order)%connectivity(8,sk) = np1d*np1d*(k  ) + np1d*(j  ) + i
        end do
      end do
    end do
    !
    ! Subtract one from the sub-connectivity so that it starts at zero
    !
    pst_geom(Geom_Hexa,order)%connectivity = pst_geom(Geom_Hexa,order)%connectivity - 1
    !
    ! Deallocate counter and counter3 before going on to the next order
    !
    deallocate ( counter , stat=ierr )
    call alloc_error(pname,"counter",2,ierr)
    deallocate ( counter3 , stat=ierr )
    call alloc_error(pname,"counter3",2,ierr)
    !
  end do loop_over_orders
  !
end subroutine get_unit_cell_subconnectivity
!
!###############################################################################
!
include "Functions/last.f90"
!
!###############################################################################
!
end module connectivity_mod
