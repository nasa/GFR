subroutine read_grid_file()
  !
  ! Routine that uses the variable grid_format to choose
  ! which routine to use in order to read the grid file.
  !        
continue
  !
  if (grid_format == patran_neutral) then
    call read_patran_neutral_gridfile
  else if (grid_format == gambit_neutral) then
    call read_gambit_neutral_gridfile
  else if (grid_format == gmsh_file) then
    call read_gmsh_gridfile
  end if
  !
  ! Check that all the nodes are in counter-clockwise
  ! ordering for each cell and adjust if needed
  !
  call check_node_ordering
  !
  ! Adjust the grid coordinates and nodes_of_cell
  ! connectivity to account for ghost cells
  !
  call adjust_for_ghost_cells
  !
end subroutine read_grid_file
!
!###############################################################################
!
subroutine read_gmsh_gridfile()
  !
  ! Routine to read in a grid file using GMSH formatting
  !
  ! Explanation of values given in array cell_nodes
  !
  ! 1  : 2-node line
  ! 2  : 3-node triangle
  ! 3  : 4-node quadrangle
  ! 4  : 4-node tetrahedron
  ! 5  : 8-node hexahedron
  ! 6  : 6-node prism
  ! 7  : 5-node pyramid
  ! 8  : 3-node second order line
  !        2 nodes associated with the vertices
  !        1 with the edge
  ! 9  : 6-node second order triangle
  !        3 nodes associated with the vertices
  !        3 with the edges
  ! 10 : 9-node second order quadrangle
  !        4 nodes associated with the vertices
  !        4 with the edges
  !        1 with the face
  ! 11 : 10-node second order tetrahedron
  !        4 nodes associated with the vertices
  !        6 with the edges
  ! 12 : 27-node second order hexahedron
  !        8 nodes associated with the vertices
  !        12 with the edges
  !        6 with the faces
  !        1 with the volume
  ! 13 : 18-node second order prism
  !        6 nodes associated with the vertices
  !        9 with the edges
  !        3 with the quadrangular faces
  ! 14 : 14-node second order pyramid
  !        5 nodes associated with the vertices
  !        8 with the edges
  !        1 with the quadrangular face
  ! 15 : 1-node point
  ! 16 : 8-node second order quadrangle
  !        4 nodes associated with the vertices
  !        4 with the edges
  ! 17 : 20-node second order hexahedron
  !        8 nodes associated with the vertices
  !        12 with the edges
  ! 18 : 15-node second order prism
  !        6 nodes associated with the vertices
  !        9 with the edges
  ! 19 : 13-node second order pyramid
  !        5 nodes associated with the vertices
  !        8 with the edges
  ! 20 : 9-node third order incomplete triangle
  !        3 nodes associated with the vertices
  !        6 with the edges
  ! 21 : 10-node third order triangle
  !        3 nodes associated with the vertices
  !        6 with the edges
  !        1 with the face
  ! 22 : 12-node fourth order incomplete triangle
  !        3 nodes associated with the vertices
  !        9 with the edges
  ! 23 : 15-node fourth order triangle
  !        3 nodes associated with the vertices
  !        9 with the edges
  !        3 with the face
  ! 24 : 15-node fifth order incomplete triangle
  !        3 nodes associated with the vertices
  !        12 with the edges
  ! 25 : 21-node fifth order complete triangle
  !        3 nodes associated with the vertices
  !        12 with the edges
  !        6 with the face
  ! 26 : 4-node third order edge
  !        2 nodes associated with the vertices
  !        2 internal to the edge
  ! 27 : 5-node fourth order edge
  !        2 nodes associated with the vertices
  !        3 internal to the edge
  ! 28 : 6-node fifth order edge 
  !        2 nodes associated with the vertices
  !        4 internal to the edge
  ! 29 : 20-node third order tetrahedron
  !        4 nodes associated with the vertices
  !        12 with the edges
  !        4 with the faces
  ! 30 : 35-node fourth order tetrahedron
  !        4 nodes associated with the vertices
  !        18 with the edges
  !        12 with the faces
  !        1 in the volume
  ! 31 : 56-node fifth order tetrahedron
  !        4 nodes associated with the vertices
  !        24 with the edges
  !        24 with the faces
  !        4 in the volume
  ! 92 : 64-node third order hexahedron
  !        8 nodes associated with the vertices
  !        24 with the edges
  !        24 with the faces
  !        8 in the volume
  ! 93 : 125-node fourth order hexahedron
  !        8 nodes associated with the vertices
  !        36 with the edges
  !        54 with the faces
  !        27 in the volume
  !
  integer, parameter, dimension(1:33) :: &
      cell_nodes = [ 2,  3,   4,  4,  8, 6,  5,  3,  6,  9, &
                    10, 27,  18, 14,  1, 8, 20, 15, 13,  9, &
                    10, 12,  15, 15, 21, 4,  5,  6, 20, 35, &
                    56, 64, 125]
continue
  !
end subroutine read_gmsh_gridfile
!
!###############################################################################
!
subroutine read_gambit_neutral_gridfile()
  !
  ! Routine to read in a grid file using Gambit Neutral formatting.
  !
  use geovar, only : nr,nbfai,nnode,ncell,nfbnd,xyz_nodes
  use geovar, only : bface,nodes_of_cell_ptr,nodes_of_cell
  use ovar,   only : iop_solpts,iop_flxpts
  !
  !.. Local Scalars ..
  integer   :: i,j,k,n,nf,n1,n2,i1,nc,pc,pf
  integer   :: icell,inode,iside,ierr,itype
  integer   :: ngrps,nbsets,ndfvl,nentry,nvalues
  !
  character(len=80) :: text
  !
  !.. Local Allocatable Arrays ..
  logical, allocatable, dimension(:) :: bc_done
  !
continue
  !
  ! Open the grid file
  !
  open (iogrd,file=gridfile,status="old",action="read", &
        position="rewind",iostat=ierr)
  call io_error("read_gambit_neutral_gridfile",trim(adjustl(gridfile)),1,ierr)
  !
  ! Skip the file header
  !
  do i = 1,6
    read (iogrd,1) text
  end do
  !
  ! Read the grid size and dimensions
  !
  read (iogrd,*) nnode,ncell,ngrps,nbsets,nr,ndfvl
  !
  allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , stat=ierr )
  call alloc_error("read_gambit_neutral_gridfile","nodes_of_cell_ptr",1,ierr)
  !
  ! Skip some lines
  !
  read (iogrd,1) text
  read (iogrd,1) text
  !
  ! Skip over the node coordinates
  !
  do i = 1,nnode
    read (iogrd,1) text
  end do
  !
  ! Skip some more lines
  !
  read (iogrd,1) text
  read (iogrd,1) text
  !
  ! Read through the cell connectivities to
  ! find the number of nodes defining each cell
  !
  do i = 1,ncell
    read (iogrd,*) n,itype,nodes_of_cell_ptr(n+1)
  end do
  !
  ! Reshuffle nodes_of_cell_ptr
  !
  do i = 2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(i) = nodes_of_cell_ptr(i) + nodes_of_cell_ptr(i-1)
  end do
  !
  ! Skip the group sections
  !
  read (iogrd,1) text
  do i = 1,ngrps
    read (iogrd,1) text
    do
      read (iogrd,1) text
      if (trim(adjustl(text)) == "ENDOFSECTION") exit
    end do
  end do
  !
  ! Cycle through the different boundary conditions
  ! to find the total number of boundary faces
  !
  nfbnd = 0
  do i = 1,nbsets
    read (iogrd,1) text
    read (iogrd,*) text,itype,nentry,nvalues
    do j = 1,nentry
      read (iogrd,1) text
    end do
    read (iogrd,1) text
    nfbnd = nfbnd + nentry
  end do
  !
  ! Now go back to the beginning and read in the grid data
  !
  rewind (iogrd)
  !
  ! Allocate the arrays needed to store this data
  !
  allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , stat=ierr )
  call alloc_error("read_gambit_neutral_gridfile","xyz_nodes",1,ierr)
  allocate ( nodes_of_cell(1:last(nodes_of_cell_ptr)) , source=0 , stat=ierr )
  call alloc_error("read_gambit_neutral_gridfile","nodes_of_cell",1,ierr)
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , stat=ierr )
  call alloc_error("read_gambit_neutral_gridfile","bface",1,ierr)
  allocate ( bc_done(1:nfbnd) , source=true , stat=ierr )
  call alloc_error("read_gambit_neutral_gridfile","bc_done",1,ierr)
  !
  ! Skip the header and other data that is already known
  !
  do i = 1,9
    read (iogrd,1) text
  end do
  !
  ! Read in the nodal coordinates
  !
  do i = 1,nnode
    read (iogrd,*) n,(xyz_nodes(k,n),k=1,nr)
  end do
  !
  ! Skip a couple lines
  !
  read (iogrd,1) text
  read (iogrd,1) text
  !
  ! Read the nodes defining each cell into nodes_of_cell
  !
  do i = 1,ncell
    read (iogrd,*) n,itype,j,(nodes_of_cell(k),k=nodes_of_cell_ptr(n)+1, &
                                                 nodes_of_cell_ptr(n+1))
  end do
  !
  ! Skip the group sections
  !
  read (iogrd,1) text
  do i = 1,ngrps
    read (iogrd,1) text
    do
      read (iogrd,1) text
      if (trim(adjustl(text)) == "ENDOFSECTION") exit
    end do
  end do
  !
  ! Read in all the boundary condition sets
  !
  n = 0
  do i = 1,nbsets
    !
    read (iogrd,1) text
    read (iogrd,*) text,itype,nentry,nvalues
    !
    if (trim(adjustl(text)) == "SubInflow") then
      itype = bc_sub_inflow
    else if (trim(adjustl(text)) == "SubOutflow") then
      itype = bc_sub_outflow
    else if (trim(adjustl(text)) == "SupInflow") then
      itype = bc_sup_inflow
    else if (trim(adjustl(text)) == "SupOutflow") then
      itype = bc_sup_outflow
    else if (trim(adjustl(text)) == "Characteristic") then
      itype = bc_characteristic
    else if (trim(adjustl(text)) == "Wall") then
      itype = bc_slip_wall
    else if (trim(adjustl(text)) == "FreeStream") then
      itype = bc_freestream
    else if (trim(adjustl(text)) == "Periodic") then
      itype = bc_periodic
      bc_done(n+1:n+nentry) = fals
    else
      itype = bc_unknown
    end if
    !
    do j = 1,nentry
      n = n + 1
      bface(1,n) = itype
      read (iogrd,*) (bface(k,n),k=2,4),(bface(4+k,n),k=1,nvalues)
    end do
    !
    read (iogrd,1) text
    !
  end do
  !
  close (iogrd,iostat=ierr)
  call io_error("read_gambit_neutral_gridfile",trim(adjustl(gridfile)),2,ierr)
  !
  ! Save the index for the ghost cell on the face in bface(10,:)
  !
  bface(10,:) = ncell + [1:nfbnd]
  !
  ! Lets do a quick modification to the periodic boundary conditions
  ! For each periodic boundary condition, we will find the opposite
  ! connecting periodic boundary condition and link the two together
  !
  do nf = 1,nfbnd
    !
    ! Skip the boundary faces that are not periodic or have been completed
    !
    if (bc_done(nf)) cycle
    !
    ! Cells that are included in this periodic bc
    !
    nc = bface(2,nf)
    pc = bface(5,nf)
    !
    ! Find the other periodic boundary face with these two cells
    !
    search_loop: do n = 1,nfbnd
      if (bc_done(n)) cycle
      if (bface(2,n)==pc .and. bface(5,n)==nc) then
        pf = n
        exit search_loop
      end if
    end do search_loop
    !
    ! Link the two boundary faces to each other
    !
    bface(5,nf) = pf          ! boundary face of connecting cell
    bface(6,nf) = pc          ! connecting cell
    bface(8,nf) = bface(3,pf) ! connecting cell geometry type
    bface(9,nf) = bface(4,pf) ! face of connecting cell on boundary
    !
    bface(5,pf) = nf          ! boundary face of connecting cell
    bface(6,pf) = nc          ! connecting cell
    bface(8,pf) = bface(3,nf) ! connecting cell geometry type
    bface(9,pf) = bface(4,nf) ! face of connecting cell on boundary
    !
    ! Mark these two periodic boundary faces as finished
    !
    bc_done(nf) = true
    bc_done(pf) = true
    !
  end do
  !
  deallocate ( bc_done , stat=ierr )
  call alloc_error("read_gambit_neutral_gridfile","bc_done",2,ierr)
  !
  ! Make sure that Lobatto nodes are requested if the grid contains triangle
  ! elements since Gauss nodes is not currently set up for triangles.
  !
  do n = 1,ncell
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    if (n2-n1+1 == 3) then
      write (*,3)
      call stop_recon(stop_mpi,"read_gambit_neutral_gridfile")
    end if
  end do
  if (iop_solpts /= Legendre_Gauss_Lobatto .or. &
      iop_flxpts /= Legendre_Gauss_Lobatto) then
    do n = 1,ncell
      n1 = nodes_of_cell_ptr(n)+1
      n2 = nodes_of_cell_ptr(n+1)
      if (n2-n1+1 == 3) then
        write (*,2)
        call stop_recon(stop_mpi,"read_gambit_neutral_gridfile")
      end if
    end do
  end if
  !
  ! Format statements
  !
  1 format (a)
  2 format (/," Quadrature points besides Gauss-Lobatto nodes are not",/, &
              " currently supported for a grid containing triangle cells.",//)
  3 format (///," Triangles are temporarily disabled due to significant",/, &
                " changes in the quadrilateral part of the code. These ",/, &
                " changes have yet to be applied to triangles.",///)
  !
end subroutine read_gambit_neutral_gridfile
!
!###############################################################################
!
subroutine read_patran_neutral_gridfile
  !
  use geovar
  !
  !  converts a patran neutral file generated by Gridgen using the
  !  NCC Analysis SoftWare selection to geometry output for our purpose
  !
  !.. Local Allocatable Arrays ..
  integer,   allocatable, dimension(:,:) :: edge,c_vtemp
  real(vbp), allocatable, dimension(:)   :: x,y,area_
  !
  !.. Local Scalars ..
  integer   :: ndum,i,k,ierr
  real(vbp) :: rdum
  !
continue
  !
  open(iogrd,file=trim(adjustl(gridfile)),status="old",action="read", &
             position="rewind",iostat=ierr)
  call io_error("read_patran_neutral_gridfile",trim(adjustl(gridfile)),1,ierr)
  !
 !open(ioxyd,file="x_y.data",      status="replace",action="write")
 !open(iocnd,file="cell_node.data",status="replace",action="write")
 !open(iobcd,file="b_c.data",      status="replace",action="write")
 !!
 !!  read patran neutral file to get nnode, ncell, and nfbnd
 !!
 !read (iogrd,*)
 !read (iogrd,*)
 !read (iogrd,*) ndum,ndum,ndum,ndum, nnode, ncell
 !read (iogrd,*)
 !!
 !!   allocate memory for x and y and c_vtemp
 !allocate (x(nnode))
 !allocate (y(nnode))
 !allocate (c_vtemp(3,ncell))
 !!
 !do i = 1,nnode
 !  read (iogrd,*)
 !  read (iogrd,11) x(i),y(i)
 !  read (iogrd,*)
 !end do
 !!
 !do i = 1,ncell
 !  read (iogrd,*)
 !  read (iogrd,*)
 !  read (iogrd,*) c_vtemp(1,i),c_vtemp(2,i),c_vtemp(3,i)
 !end do
 !!
 !nfbnd = 0
 !do
 !  read (iogrd,*,end = 20)
 !  read (iogrd,*)
 !  read (iogrd,*)
 !  nfbnd = nfbnd+1
 !end do
 !20 continue
 !!
 !deallocate (x, y, c_vtemp)
 !!
 !!.. with nnode, ncell, and nfbnd known, now allocate memory
 !!... for geovar
 !allocate (xyz_nodes(2,-nfbnd:nnode))
 !allocate (c_v(3,-nfbnd:ncell), c_e(3,-nfbnd:ncell),c_nbr(3,ncell))
 !allocate (be_v(2,nfbnd), be_ic(nfbnd), be_e(nfbnd), &
 !          be_pe(nfbnd), be_pc(nfbnd), id_be(nfbnd))
 !!
 !!now local memory
 !allocate (area_(-nfbnd:ncell))
 !!...
 !rewind(iogrd)
 !read (iogrd,*)
 !read (iogrd,*)
 !read (iogrd,*)
 !read (iogrd,*)
 !!
 !do i = 1,nnode
 !  read (iogrd,*)
 !  read (iogrd,11)xyz_nodes(1,i),xyz_nodes(2,i)
 !  read (iogrd,*)
 !end do
 !!
 !do i = 1,ncell
 !  read (iogrd,*)
 !  read (iogrd,*)
 !  read (iogrd,*) c_v(1,i),c_v(2,i),c_v(3,i)
 !end do
 !!
 !!... we do not use the value i = 0, but initialze anyway
 !xyz_nodes(1,0) = zero
 !xyz_nodes(2,0) = zero
 !!
 !!   allocate memory for edge
 !allocate (edge(3,nfbnd))
 !!
 !do i = 1,nfbnd
 !  read (iogrd,*) ndum,be_ic(i)
 !  read (iogrd,10) ndum,edge(1,i),edge(2,i),edge(3,i)
 !  read (iogrd,*) rdum
 !  id_be(i) = int(rdum)
 !end do
 !!
 !!   determine which edge is on boundary
 !do i = 1,nfbnd
 !  !
 !  if (edge(3,i) == 0) then
 !    be_v(1,i) = c_v(1,be_ic(i))
 !    be_v(2,i) = c_v(2,be_ic(i))
 !  else
 !    if (edge(2,i) == 0) then
 !      be_v(1,i) = c_v(3,be_ic(i))
 !      be_v(2,i) = c_v(1,be_ic(i))
 !    else
 !      if (edge(1,i) == 0) then
 !        be_v(1,i) = c_v(2,be_ic(i))
 !        be_v(2,i) = c_v(3,be_ic(i))
 !      end if
 !    end if
 !  end if
 !  !
 !end do
 !!
 !!  write output for x and y coordinate data; cell to vertex data;
 !!  and boundary condition data.
 !!
 !write (ioxyd,*) nnode
 !do i = 1,nnode
 !  write (ioxyd,15) xyz_nodes(1,i),xyz_nodes(2,i)
 !end do
 !!
 !write (iocnd,*) ncell
 !do i = 1,ncell
 !  write (iocnd,211) c_v(1,i),c_v(2,i),c_v(3,i)
 !end do
 !!
 !! reset id_be to match cese convention
 !do i = 1,nfbnd
 !  !
 !  if (id_be(i) == 1) then ! inflow
 !    id_be(i) = 1
 !  else if (id_be(i) == 10) then ! inside wall
 !    id_be(i) = 0
 !  else if (id_be(i) == 11) then ! outside wall
 !    id_be(i) = 3
 !  else if (id_be(i) == 12) then ! outflow
 !    id_be(i) = 2
 !  end if
 !  !
 !end do
 !!
 !write (iobcd,*) nfbnd,0
 !write (iobcd,222) (id_be(i),i = 1,nfbnd)
 !write (iobcd,222) (be_ic(i),i = 1,nfbnd)
 !write (iobcd,222) ((be_v(k,i),k = 1,2),i = 1,nfbnd)
 !!
 !deallocate (edge)
 !!
 !close (ioxyd)
 !close (iocnd)
 !close (iobcd)
  close (iogrd,iostat=ierr)
  call io_error("read_patran_neutral_gridfile",trim(adjustl(gridfile)),2,ierr)
  !
  11 format (2e16.8)
  10 format (i9,i1,i1,i1)
  15 format (2e16.8)
  211 format (1x,3i10)
  222 format (1x,10i8)
  !
end subroutine read_patran_neutral_gridfile
!
!###############################################################################
!
subroutine check_node_ordering()
  !
  !.. Use Statements ..
  use geovar, only : ncell,nfbnd,nodes_of_cell_ptr, &
                     xyz_nodes,nodes_of_cell,bface
  !
  !.. Local Scalars ..
  integer   :: i,n,nf,n1,n2,ip1,ip2,ip3,ip4
  real(vbp) :: x31,y31,x42,y42,ss
  !
continue
  !
  ! Lets do a quick test on the nodes defining each cell to
  ! make sure they are all in counter-clockwise order. If they
  ! arent, switch the ordering of the nodes so that they are.
  !
  ! First lets look at the cells on the grid boundaries since we
  ! also need to change bface(4) if the ordering is wrong
  !
  do nf = 1,nfbnd
    !
    n = bface(2,nf)
    !
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    !
    ip1 = nodes_of_cell(n1+0)
    ip2 = nodes_of_cell(n1+1)
    ip3 = nodes_of_cell(n1+2)
    ip4 = nodes_of_cell(n2  ) ! for triangles, ip4 = ip3 because n2 = n1+2
    !
    x31 = xyz_nodes(1,ip3) - xyz_nodes(1,ip1)
    y31 = xyz_nodes(2,ip3) - xyz_nodes(2,ip1)
    !
    x42 = xyz_nodes(1,ip4) - xyz_nodes(1,ip2) ! for triangles x42 = x3-x2
    y42 = xyz_nodes(2,ip4) - xyz_nodes(2,ip2) ! for triangles y42 = y3-y2
    !
    ! Take the cross product - x31 X x42
    !
    ss = x31*y42 - x42*y31
    !
    ! If the cross product is negative, this means the cell is in clockwise
    ! ordering.  Switch two nodes (2 and 3 for triangles; 2 and 4 for quads)
    ! to change the cell to counter-clockwise ordering.
    !
    if (ss < zero) then
      !
      nodes_of_cell(n1+1) = ip4 ! for triangles, ip4 = ip3
      nodes_of_cell(n2  ) = ip2
      !
      ! We also need to adjust the value in bface(4,nf) depending on cell type
      ! For triangles, faces need to go from 1,2,3 to 3,2,1 ordering.
      !                This is done using:  bface(4) = 4 - bface(4)
      ! For quad, faces need to go from 1,2,3,4 to 4,3,2,1 ordering.
      !                This is done using:  bface(4) = 5 - bface(4)
      !
      bface(4,nf) = merge(4,5,bface(3,nf)==3) - bface(4,nf)
      !
      ! Now we have to loop through the remaining boundary faces
      ! and apply this correction to bface(4) for any other boundary
      ! faces that have this same host cell. Otherwise, when it gets
      ! to the other boundary face > nf with this same host cell, the
      ! node ordering will already be fixed and it will not make it
      ! into this ss<zero conditional to apply a correction to bface(4).
      !
      do i = nf+1,nfbnd
        if (bface(2,i) == n) then
          bface(4,i) = merge(4,5,bface(3,i)==3) - bface(4,i)
        end if
      end do
      !
    end if
    !
  end do
  !
  ! Now lets do it over all the cells in the grid. This includes checking
  ! the cells on the grid boundaries a second time but they shouldnt be
  ! adjusted this time and the additional cost shouldnt hurt too much.
  !
  do n = 1,ncell
    !
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    !
    ip1 = nodes_of_cell(n1+0)
    ip2 = nodes_of_cell(n1+1)
    ip3 = nodes_of_cell(n1+2)
    ip4 = nodes_of_cell(n2  ) ! for triangles, ip4 = ip3 because n2 = n1+2
    !
    x31 = xyz_nodes(1,ip3) - xyz_nodes(1,ip1)
    y31 = xyz_nodes(2,ip3) - xyz_nodes(2,ip1)
    !
    x42 = xyz_nodes(1,ip4) - xyz_nodes(1,ip2) ! for triangles x42 = x3-x2
    y42 = xyz_nodes(2,ip4) - xyz_nodes(2,ip2) ! for triangles y42 = y3-y2
    !
    ! Take the cross product - x31 X x42
    !
    ss = x31*y42 - x42*y31
    !
    ! If the cross product is negative, this means the cell is in clockwise
    ! ordering.  Switch two nodes (2 and 3 for triangles; 2 and 4 for quads)
    ! to change the cell have counter-clockwise ordering.
    !
    if (ss < zero) then
      !
      nodes_of_cell(n1+1) = ip4 ! for triangles, ip4 = ip3
      nodes_of_cell(n2  ) = ip2
      !
    end if
    !
  end do
  !
end subroutine check_node_ordering
!
!###############################################################################
!
subroutine adjust_for_ghost_cells()
  !
  ! Routine to make some adjustments to the data from the grid file
  ! to account for ghost cells and new nodes associated with them 
  !
  !.. Use Statements ..
  use geovar, only : nr,ncell,nfbnd,nnode,n_totcel,n_totnod, &
                     bface,nodes_of_cell_ptr,nodes_of_cell,xyz_nodes
 !use ovar,   only : Period_Time,Period_Timesteps
 !use ovar,   only : Final_Time,Total_Periods
 !use ovar,   only : vxref,vyref,aref,dt_in,machref
  use ovar,   only : Grid_Min_X,Grid_Max_X,Grid_Length_X
  use ovar,   only : Grid_Min_Y,Grid_Max_Y,Grid_Length_Y
  !
  !.. Local Scalars ..
  integer   :: g1,g2,i,i1,i2,n1,n2,nf,ng,ni,np
  integer   :: gcell,icell,igeom,inode,iside
  integer   :: pf,pcell,p1,p2,pgeom,pside
  integer   :: np1,np2,ni1,ni2,ierr
  integer   :: Period_Timesteps_X,Period_Timesteps_Y
  real(vbp) :: Period_Time_X,Period_Time_Y
  real(vbp) :: theta,real_Total_Periods
  !
  !.. Local Arrays ..
  real(vbp), dimension(1:nr)     :: ptf,rnrm,ds
  real(vbp), dimension(1:nr)     :: ds1,ds2,pnrm,inrm
  real(vbp), dimension(1:nr,1:4) :: xyp
  !
  !.. Local Allocatable Arrays ..
  integer,   allocatable, dimension(:)   :: tmp_int
  real(vbp), allocatable, dimension(:,:) :: tmp_real
  !
  integer, parameter :: rotation = 0
  !
continue
  !
  ! Before we start adding ghost cells, lets get the extents of the grid domain
  ! NOTE: ASSUME GRID BOUNDARIES OF CONSTANT X AND Y
  !
  Grid_Min_X = minval(xyz_nodes(1,:))
  Grid_Max_X = maxval(xyz_nodes(1,:))
  Grid_Min_Y = minval(xyz_nodes(2,:))
  Grid_Max_Y = maxval(xyz_nodes(2,:))
  Grid_Length_X = Grid_Max_X - Grid_Min_X
  Grid_Length_Y = Grid_Max_Y - Grid_Min_Y
  !
 !if (vxref < eps3) then
 !  !
 !  Period_Time_Y = Grid_Length_Y / machref
 !  !
 !  Period_Timesteps_Y = nint(Period_Time_Y/dt_in)
 !  !
 !  Period_Timesteps = Period_Timesteps_Y
 !  !
 !else if (vyref < eps3) then
 !  !
 !  Period_Time_X = Grid_Length_X / machref
 !  !
 !  Period_Timesteps_X = nint(Period_Time_X/dt_in)
 !  !
 !  Period_Timesteps = Period_Timesteps_X
 !  !
 !else
 !  !
 !  Period_Time_X = Grid_Length_X / (vxref/aref)
 !  Period_Time_Y = Grid_Length_Y / (vyref/aref)
 !  !
 !  Period_Timesteps_X = nint(Period_Time_X/dt_in)
 !  Period_Timesteps_Y = nint(Period_Time_Y/dt_in)
 !  !
 !  if (Period_Timesteps_X == Period_Timesteps_Y) then
 !    Period_Timesteps = Period_Timesteps_X
 !  else
 !    Period_Timesteps = Period_Timesteps_X * Period_Timesteps_Y &
 !                     / mod( max(Period_Timesteps_X,Period_Timesteps_Y), &
 !                            min(Period_Timesteps_X,Period_Timesteps_Y) )
 !  end if
 !  !
 !end if
 !!
 !Period_Time = dt_in*real(Period_Timesteps,vbp)
 !!
 !if (Total_Periods /= 0) then
 !  Final_Time = Period_Time*real(Total_Periods,vbp)
 !else
 !  real_Total_Periods = Final_Time/Period_Time
 !  Total_Periods = nint(real_Total_Periods)
 !  if (abs(real_Total_Periods-real(Total_Periods,vbp)) > ten*eps6) then
 !    write (iout,2) real_Total_Periods
 !  end if
 !end if
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Total number of interior plus ghost cells
  !
  n_totcel = ncell + nfbnd
  !
  ! Reallocate nodes_of_cell_ptr to account for ghost cells
  !
  allocate ( tmp_int(1:n_totcel+1) , source=0 , stat=ierr )
  call alloc_error("adjust_for_ghost_cells","tmp_int",1,ierr)
  ! 
  tmp_int(1:size(nodes_of_cell_ptr)) = nodes_of_cell_ptr
  call move_alloc( from=tmp_int , to=nodes_of_cell_ptr )
  !
  ! Set the number of nodes for each ghost cell
  !
  do nf = 1,nfbnd
    !
    if (bface(1,nf) /= bc_periodic) then
      !
      ! Not a periodic boundary, number of nodes is the same as interior
      !
      inode = merge(3,4,bface(3,nf)==3)
      !
    else
      !
      pf = bface(5,nf)
      pcell = bface(2,pf)
      !
      n1 = nodes_of_cell_ptr(pcell)+1
      n2 = nodes_of_cell_ptr(pcell+1)
      !
      inode = merge(3,4,n2-n1+1==3)
      !
    end if
    !
    nodes_of_cell_ptr(ncell+nf+1) = inode
    !
  end do
  !
  ! Reshuffle the ghost cell portion of nodes_of_cell_ptr
  !
  do i = ncell+2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(i) = nodes_of_cell_ptr(i) + nodes_of_cell_ptr(i-1)
  end do
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Reallocate nodes_of_cell to account for ghost cells
  ! 
  allocate ( tmp_int(1:last(nodes_of_cell_ptr)) , source=0 , stat=ierr )
  call alloc_error("adjust_for_ghost_cells","tmp_int",1,ierr)
  !
  tmp_int(1:size(nodes_of_cell)) = nodes_of_cell
  call move_alloc( from=tmp_int , to=nodes_of_cell )
  !
  ! Loop over the boundary faces and find the nodes of each ghost cell.
  ! The first two nodes are on the boundary face while the other nodes
  ! are newly created and to be defined in step 3. The nodes are in
  ! counter-clockwise ordering.
  !
  inode = nnode
  !
  do nf = 1,nfbnd
    !
    gcell = nf + ncell
    !
    ! Node indices for the ghost cell
    !
    g1 = nodes_of_cell_ptr(gcell)+1
    g2 = nodes_of_cell_ptr(gcell+1)
    !
    ! Information about boundary face
    !
    icell = bface(2,nf)
    igeom = merge(3,4,bface(3,nf)==3)
    iside = bface(4,nf)
    !
    ! Starting node index for the interior cell
    !
    i1 = nodes_of_cell_ptr(icell)+1
    !
    ! Assign nodes on boundary from interior cell
    !
    nodes_of_cell(g1+0) = nodes_of_cell(i1+mod(iside,igeom))
    nodes_of_cell(g1+1) = nodes_of_cell(i1+iside-1)
    !
    ! Create new nodes for the remaining ghost cell nodes
    !
    do i = g1+2,g2
      inode = inode + 1
      nodes_of_cell(i) = inode
    end do
    !
  end do
  !
  ! Total number of nodes after the addition of
  ! the newly created nodes for the ghost cells
  !
  n_totnod = inode
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Reallocate xyz_nodes to account for ghost nodes
  !
  allocate ( tmp_real(1:nr,1:n_totnod) , source=zero , stat=ierr )
  call alloc_error("adjust_for_ghost_cells","tmp_real",1,ierr)
  !
  tmp_real(1:nr,1:nnode) = xyz_nodes
  call move_alloc( from=tmp_real , to=xyz_nodes )
  !
  ! Loop over boundary faces to compute the
  ! coordinates of the newly created ghost nodes
  !
  do nf = 1,nfbnd
    !
    ! Ghost cell and the corresponding node indices
    !
    gcell = nf + ncell
    !
    g1 = nodes_of_cell_ptr(gcell)+1
    g2 = nodes_of_cell_ptr(gcell+1)
    !
    ! Interior cell and the corresponding node indices
    !
    icell = bface(2,nf)
    !
    i1 = nodes_of_cell_ptr(icell)+1
    i2 = nodes_of_cell_ptr(icell+1)
    !
    ! Information about boundary face
    !
    igeom = merge(3,4,bface(3,nf)==3)
    iside = bface(4,nf)
    !
    ! Conditional to treat periodic boundary conditions separately
    !
    if (bface(1,nf) /= bc_periodic) then
      !
      ! This is not a periodic boundary face
      !
      ! Nodes on face
      !
      n1 = nodes_of_cell(g1)
      n2 = nodes_of_cell(g1+1)
      !
      ! Center of boundary face by averaging the coordinates of the two boundary nodes
      !
      ptf(:) = half*(xyz_nodes(:,n1)+xyz_nodes(:,n2))
      !
      ! Outward unit normal to boundary face. NOTE, we are using the
      ! ordering of the boundary nodes for the ghost cell which is
      ! opposite that of the interior cell so the normal is the negative
      ! of the regular normal equation. (We are essentially finding the
      ! inward unit normal of the boundary face within the ghost cell)
      !
      rnrm(1) = xyz_nodes(2,n1) - xyz_nodes(2,n2)
      rnrm(2) = xyz_nodes(1,n2) - xyz_nodes(1,n1)
      rnrm = unit_vector(rnrm)
      !
      ! Coordinates for nodes g1 and g1+1 of the ghost cell have already been defined
      ! so we need to define nodes g1+2 and g2 (g1+2==g2 for triangular cells).
      !
      if (igeom == 3) then
        !
        ! Ghost node 3 for triangle
        !
        ni = nodes_of_cell( i1 + mod(iside+1,igeom) )
        ng = nodes_of_cell( g1 + 2 )
        !
        ! Get the coordinates of ghost node 3
        !
        xyz_nodes(:,ng) = merge(rotate(xyz_nodes(:,n1),xyz_nodes(:,n2),xyz_nodes(:,ni)), &
                                reflect(xyz_nodes(:,ni),ptf,rnrm), &
                                rotation == 1)
        !
      else if (igeom == 4) then
        !
        ! Ghost node 3 for quadrilateral
        !
        ng = nodes_of_cell( g1 + 2 )
        ni = nodes_of_cell( i1 + mod(iside+2-max(0,min(1,rotation)),igeom) )
        !
        ! Get the coordinates of ghost node 3
        !
        xyz_nodes(:,ng) = merge(rotate(xyz_nodes(:,n1),xyz_nodes(:,n2),xyz_nodes(:,ni)), &
                                reflect(xyz_nodes(:,ni),ptf,rnrm), &
                                rotation == 1)
        !
        ! Ghost node 4 for quadrilateral
        !
        ng = nodes_of_cell( g1 + 3 )
        ni = nodes_of_cell( i1 + mod(iside+1+max(0,min(1,rotation)),igeom) )
        !
        ! Get the coordinates of ghost node 4
        !
        xyz_nodes(:,ng) = merge(rotate(xyz_nodes(:,n1),xyz_nodes(:,n2),xyz_nodes(:,ni)), &
                                reflect(xyz_nodes(:,ni),ptf,rnrm), &
                                rotation == 1)
        !
      end if
      !
    else
      !
      ! This is a periodic boundary face
      !
      ! Boundary face linked to this periodic boundary face
      !
      pf = bface(5,nf)
      !
      ! Linked cell and the corresponding node indices
      !
      pcell = bface(2,pf)
      !
      p1 = nodes_of_cell_ptr(pcell)+1
      p2 = nodes_of_cell_ptr(pcell+1)
      !
      ! Information about the linked face
      !
      pgeom = merge(3,4,bface(3,pf)==3)
      pside = bface(4,pf)
      !
      ! Node coordinates for periodic cell
      !
      xyp(:,1:pgeom) = xyz_nodes(:,nodes_of_cell(p1:p2))
      !
      ! Nodes on the boundary face for pcell
      !
      np1 = nodes_of_cell( p1 + pside - 1 )
      np2 = nodes_of_cell( p1 + mod(pside,pgeom) )
      !
      ! Nodes on the boundary face for icell
      !
      ni1 = nodes_of_cell( i1 + iside - 1)
      ni2 = nodes_of_cell( i1 + mod(iside,igeom) )
      !
      ! Center of boundary face for pcell
      !
      ptf(:) = half * (xyz_nodes(:,np1) + xyz_nodes(:,np2))
      !
      ! Outward unit normal for boundary face of pcell
      !
      pnrm(1) = xyz_nodes(2,np2) - xyz_nodes(2,np1)
      pnrm(2) = xyz_nodes(1,np1) - xyz_nodes(1,np2)
      pnrm = unit_vector(pnrm)
      !
      ! Outward unit normal for boundary face of icell
      !
      inrm(1) = xyz_nodes(2,ni2) - xyz_nodes(2,ni1)
      inrm(2) = xyz_nodes(1,ni1) - xyz_nodes(1,ni2)
      inrm = unit_vector(inrm)
      !
      ! Get the angle between the two outward unit normals
      !
      theta = acos(dot_product(pnrm,inrm))
      !
      ! Translation of the pcell nodes so that ptf is at the origin
      !
      do i = 1,pgeom
        xyp(:,i) = xyp(:,i) - ptf(:)
      end do
      !
      ! Rotate xyp about ptf theta+180 degrees
      !
      do i = 1,pgeom
        xyp(:,i) = rotate_2d( xyp(:,i) , theta+pi )
      end do
      !
      ! Translation of the pcell nodes back to original ptf
      !
      do i = 1,pgeom
        xyp(:,i) = xyp(:,i) + ptf(:)
      end do
      !
      ! Find the distances between the two corresponding nodes of
      ! the two boundary faces. The two boundary faces should now be
      ! parallel to each other but with opposite outward unit normals
      ! and therefore the distances should be equal. 
      !
      ds1(:) = xyz_nodes(:,ni1) - xyp(:,1+mod(pside,pgeom))
      ds2(:) = xyz_nodes(:,ni2) - xyp(:,pside)
      !
      if (abs(ds1(1)-ds2(1))>eps6 .or. abs(ds1(2)-ds2(2))>eps6) then
        write (*,1) icell,pcell,abs(ds1-ds2)
      end if
      !
      ! Move the interior nodes of pcell to the ghost nodes of icell
      !
      do i = 3,pgeom
        !
        ng = nodes_of_cell( g1 + i-1 )
        np = nodes_of_cell( p1 + mod(pside+i-2,pgeom) )
        !
        xyz_nodes(:,ng) = xyz_nodes(:,np) + ds1(:)
        !
      end do
      !
    end if
    !
  end do
  !
  ! If using multiple processors, save the global grid arrays
  !
  call save_global_grid
  !
  1 format ("ICELL = ",i0,"  PCELL = ",i0,"   (delta x,delta y) = ",2es24.11)
  2 format (/," ********************* WARNING!!! *********************",/, &
              "  IT APPEARS THE FINAL TIME DOES NOT CORRESPOND TO THE",/, &
              "  EXACT END OF A PERIOD FOR THE VORTEX PROBLEM. THE",/, &
              "  FINAL TIME CORRESPONDS TO",f12.5," TOTAL PERIODS.",/, &
              " ********************* WARNING!!! *********************",/)
  !
end subroutine adjust_for_ghost_cells
!
!###############################################################################
!
subroutine save_global_grid
  !
  !.. Use Statements ..
  use geovar
  use flowvar, only : global_usp,global_vsp
  !
  !.. Local Scalars ..
  integer :: n,np,n1,n2,ierr
  !
continue
  !
  n_global_node = nnode
  n_global_cell = ncell
  n_global_fbnd = nfbnd
  n_global_totnod = n_totnod
  n_global_totcel = n_totcel
  !
  n1 = size(xyz_nodes,dim=1)
  n2 = size(xyz_nodes,dim=2)
  allocate ( global_xyz_nodes(1:n1,1:n2) , source=xyz_nodes , stat=ierr )
  call alloc_error("save_global_grid","global_xyz_nodes",1,ierr)
  !
  n1 = size(nodes_of_cell_ptr)
  allocate ( global_nofc_ptr(1:n1) , source=nodes_of_cell_ptr , stat=ierr )
  call alloc_error("save_global_grid","global_nofc_ptr",1,ierr)
  !
  n1 = size(nodes_of_cell)
  allocate ( global_nofc(1:n1) , source=nodes_of_cell , stat=ierr )
  call alloc_error("save_global_grid","global_nofc",1,ierr)
  !
  ! Create the global pts_in_cell_ptr. This will be done later
  ! but a global version is never made before the global grid
  ! is partitioned and each processor creates its own separate
  ! local grid.
  !
  allocate ( global_pinc_ptr(1:n_totcel+1) , source=0 , stat=ierr )
  call alloc_error("save_global_grid","global_pinc_ptr",1,ierr)
  !
  do n = 1,n_totcel
    np = nodes_of_cell_ptr(n+1) - nodes_of_cell_ptr(n)
    global_pinc_ptr(n+1) = global_pinc_ptr(n) + merge(n_sp2dt,n_sp2dq,np==3)
  end do
  !
  n_global_solpts = global_pinc_ptr(ncell+1)
  n_global_totpts = last(global_pinc_ptr)
  !
  ! Now the only remaining global grid array to create is global_xyz
  ! We will have to do a collective communication of the local xyz
  ! arrays to accomplish this after they are created. We will allocate
  ! the array here though since we have the parameters needed to do so.
  !
  allocate ( global_xyz(1:nr,1:n_global_totpts) , source=zero , stat=ierr )
  call alloc_error("save_global_grid","global_xyz",1,ierr)
  !
  ! Also allocate the global solution arrays for use later
  !
  allocate ( global_usp(1:neqns,1:n_global_totpts) , source=zero , stat=ierr )
  call alloc_error("save_global_grid","global_usp",1,ierr)
  allocate ( global_vsp(1:neqns,1:n_global_totpts) , source=zero , stat=ierr )
  call alloc_error("save_global_grid","global_vsp",1,ierr)
  !
  ! Create the global cell offset array
  !
  allocate ( global_cell_offset(1:n_global_totcel) , source=0 , stat=ierr )
  call alloc_error("save_global_grid","global_cell_offset",1,ierr)
  !
  do n = 1,n_global_totcel-1
    np = global_nofc_ptr(n+1) - global_nofc_ptr(n)
    global_cell_offset(n+1) = global_cell_offset(n) + &
                              merge(n_sp2dt,n_sp2dq,np==3)
  end do
  !
end subroutine save_global_grid
!
!###############################################################################
!
