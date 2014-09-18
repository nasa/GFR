module module_gmsh
  !
  !.. Use Statements ..
  use module_kind_types
  !
  implicit none
  !
  private
  !
  ! ##################################
  ! ###  Module Public Procedures  ###
  ! ##################################
  !
  public :: read_gmsh_gridfile
  !
  ! ###########################
  ! ###  Module Parameters  ###
  ! ###########################
  !
  ! CELL_NODES: gives the number of nodes defining each type of grid cell
  ! CELL_DIMEN: gives the dimensions of the cell-type (i.e., [0,1,2,3]-D)
  ! CELL_ORDER: gives the order of the cell-type
  ! CELL_GEOM : gives the shape/geometry of the cell type
  !              0 -> node/point
  !              1 -> line/edge
  !              3 -> triangle
  !              2 -> quadrilateral
  !              4 -> tetrahedron
  !              5 -> pyramid
  !              6 -> prism
  !              8 -> hexahedron
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
  !
  integer, parameter, dimension(1:33) :: &
      cell_dimen = [ 1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3, 3, 3, 0, 2, 3, &
                     3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 3, 3, 3, 3, 3 ]
  !
  integer, parameter, dimension(1:33) :: &
      cell_order = [ 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, &
                     2, 2, 3, 3, 4, 4, 5, 5, 3, 4, 5, 3, 4, 5, 3, 4 ]
  !
  integer, parameter, dimension(1:33) :: &
      cell_geom  = [ 1, 3, 2, 4, 8, 6, 5, 1, 3, 2, 4, 8, 6, 5, 0, 2, 8, &
                     6, 5, 3, 3, 3, 3, 3, 3, 1, 1, 1, 4, 4, 4, 8, 8 ]
  !
  !                  1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,
  !                 18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33
  !
  ! ###################################
  ! ###  Module Derived Data Types  ###
  ! ###################################
  !
  type cell_t
    integer :: cell_type
    integer :: cell_geom
    integer :: num_tags
    integer :: num_nodes
    integer :: host_cell = 0
    integer, allocatable, dimension(:) :: nodes
    integer, allocatable, dimension(:) :: tags
  end type cell_t
  !
  type phys_t
    integer :: phys_dim
    integer :: phys_id
    integer :: ibc
    character(len=80) :: phys_nam
    integer, allocatable, dimension(:) :: cells
  end type phys_t
  !
  ! ################################
  ! ###  Module Local Variables  ###
  ! ################################
  !
  logical, save :: meshformat_successful
  logical, save :: nodes_successful
  logical, save :: elements_successful
  logical, save :: physicalnames_successful
  !
  integer, save :: num_phys_grps
  integer, save :: num_cell_data
  !
  ! ###################################
  ! ###  Module Allocatable Arrays  ###
  ! ###################################
  !
  type(cell_t), allocatable, dimension(:) :: cell_data
  type(phys_t), allocatable, dimension(:) :: phys_grps
  !
contains
!
!###############################################################################
!
subroutine read_gmsh_gridfile( gridfile )
  !
  ! Routine to read in a grid file using GMSH formatting
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: gridfile
  !
  !.. Local Scalars ..
  integer :: ierr
  logical :: end_of_file
  character(len=80) :: section_header
  !
  !.. Local Allocatable Arrays ..
  !
continue
  !
  meshformat_successful = fals
  nodes_successful = fals
  elements_successful = fals
  physicalnames_successful = fals
  !
  ! Open the grid file
  !
  open (iogrd,file=gridfile,status="old",action="read", &
       position="rewind",iostat=ierr)
  call io_error("read_gmsh_gridfile",trim(adjustl(gridfile)),1,ierr)
  !
  ! Reading the grid file sections, looping until reaching the end of the file
  !
  do
    !
    call find_gmsh_section(section_header,end_of_file)
    !
    ! If the end_of_file flag was set, exit the read loop
    !
    if (end_of_file) exit
    !
    ! Read the next section based on the section header
    !
    if (section_header == "MESHFORMAT") then
      call read_gmsh_section_meshformat
    else if (section_header == "NODES") then
      call read_gmsh_section_nodes
    else if (section_header == "ELEMENTS") then
      call read_gmsh_section_elements
    else if (section_header == "PHYSICALNAMES") then
      call read_gmsh_section_physicalnames
    else if (section_header == "PERIODIC") then
      call read_gmsh_section_periodic
    else if (section_header == "NODEDATA") then
      call read_gmsh_section_nodedata
    else if (section_header == "ELEMENTDATA") then
      call read_gmsh_section_elementdata
    else if (section_header == "ELEMENTNODEDATA") then
      call read_gmsh_section_elementnodedata
    else if (section_header == "INTERPOLATIONSCHEME") then
      call read_gmsh_section_interpolationscheme
    else
      !
      ! The section header wasnt recognized so continue reading
      ! the grid file until a valid section header is found
      !
      cycle
      !
    end if
    !
  end do
  !
  ! Check that all the necessary sections were successfully read
  !
  if ( .not. meshformat_successful ) then
    !
    if (mypnum == 0) write (iout,2)
    call stop_recon(stop_mpi,"read_gmsh_gridfile")
    !
  else if ( .not. nodes_successful ) then
    !
    if (mypnum == 0) write (iout,3)
    call stop_recon(stop_mpi,"read_gmsh_gridfile")
    !
  else if ( .not. elements_successful ) then
    !
    if (mypnum == 0) write (iout,4)
    call stop_recon(stop_mpi,"read_gmsh_gridfile")
    !
  else if ( .not. physicalnames_successful ) then
    !
    if (mypnum == 0) write (iout,5)
    call stop_recon(stop_mpi,"read_gmsh_gridfile")
    !
  else
    !
    ! We are done reading the grid file so close it
    !
    close (iogrd,iostat=ierr)
    call io_error("read_gmsh_gridfile",trim(adjustl(gridfile)),2,ierr)
    !
  end if
  !
  ! Finish creating the arrays defining the grid geometry.
  ! bface, nodes_of_cell_ptr, nodes_of_cell, etc.
  !
  call create_gmsh_grid_arrays
  !
  ! Format statements
  !
  1 format (a)
  2 format (/,"******************** ERROR! *********************",/, &
              "  Aborting because there was an error reading the",/, &
              "  'MeshFormat' section of the GMSH grid file!",/, &
              "******************** ERROR! *********************",/)
  3 format (/,"******************** ERROR! *********************",/, &
              "  Aborting because there was an error reading the",/, &
              "  'Nodes' section of the GMSH grid file!",/, &
              "******************** ERROR! *********************",/)
  4 format (/,"******************** ERROR! *********************",/, &
              "  Aborting because there was an error reading the",/, &
              "  'Elements' section of the GMSH grid file!",/, &
              "******************** ERROR! *********************",/)
  5 format (/,"******************** ERROR! *********************",/, &
              "  Aborting because there was an error reading the",/, &
              "  'PhysicalNames' section of the GMSH grid file!",/, &
              "******************** ERROR! *********************",/)
  !
end subroutine read_gmsh_gridfile
!
!###############################################################################
!
subroutine find_gmsh_section(text_out,end_of_file)
  !
  !.. Formal Arguments ..
  logical,          intent(out) :: end_of_file
  character(len=*), intent(out) :: text_out
  !
  !.. Local Scalars ..
  integer :: ierr,n
  character(len=len(text_out)) :: text
  !
continue
  !
  end_of_file = fals
  text_out = ""
  !
  ! Find the next line beginning with "$"
  !
  do
    !
    read (iogrd,1,iostat=ierr,end=100,err=101) text
    text = trim(adjustl(text))
    !
    if (text(1:1) == "$") then
      n = len_trim(text)
      text_out = uppercase(text(2:n))
      exit
    else
      cycle
    end if
    !
  end do
  !
  ! Return if successfully found a new section header
  !
  return
  !
  ! Read error continuation statements
  !
  100 end_of_file = true
      return
  !
  101 if (mypnum == 0) write (iout,2) ierr
      call stop_recon(stop_mpi,"find_gmsh_section")
  !
  ! Format Statements
  !
  1 format (a)
  2 format (/,"******************** ERROR! *********************",/, &
              "    Error trying to read the input grid file!", &
              "       IO Error = ",i0,/, &
              "******************** ERROR! *********************",/)
  !
end subroutine find_gmsh_section
!
!###############################################################################
!
subroutine read_gmsh_section_meshformat
  !
  !.. Local Scalars ..
  integer :: ierr,gmsh_file_type,gmsh_data_size
  logical :: end_of_file
  real(wp) :: gmsh_version_num
  character(len=80) :: text
  !
continue
  !
  ! Read the MeshFormat section of the GMSH grid file
  !
  read (iogrd,*,iostat=ierr,end=100,err=101) gmsh_version_num, &
                                             gmsh_file_type, &
                                             gmsh_data_size
  !
  ! Check that the GMSH version number is supported
  !
  gmsh_version_num = gmsh_version_num * ten
  !
  if (nint(gmsh_version_num) /= 22) then
    if (mypnum == 0) write (iout,1)
    call stop_recon(stop_mpi,"read_gmsh_section_meshformat")
  end if
  !
  ! Return if the line marking the end of the MeshFormat section is found
  !
  call find_gmsh_section(text,end_of_file)
  !
  if (text == "ENDMESHFORMAT") then
    meshformat_successful = true
    return
  else if (end_of_file) then
    if (mypnum == 0) write (iout,2)
  else
    if (mypnum == 0) write (iout,3)
  end if
  !
  call stop_recon(stop_mpi,"read_gmsh_section_meshformat")
  !
  ! Read error continuation statements
  !
  100 if (mypnum == 0) write (iout,2) ierr
      call stop_recon(stop_mpi,"read_gmsh_section_meshformat")
  !
  101 if (mypnum == 0) write (iout,4) ierr
      call stop_recon(stop_mpi,"read_gmsh_section_meshformat")
  !
  ! Format Statements
  !
  1 format (/,"******************** ERROR! *********************",/, &
              "  The version of the GMSH grid file is",/, &
              "  not supported! Currently, only version",/, &
              "  2.2 of the GMSH format is supported!",/, &
              "******************** ERROR! *********************",/)
  2 format (/,"******************** ERROR! *********************",/, &
              "  Pre-mature end-of-file while trying to read",/, &
              "  the 'MeshFormat' section of the GMSH grid file!",/, &
              "******************** ERROR! *********************",/)
  3 format (/,"*********************** ERROR! ************************",/, &
              "  While trying to read the 'MeshFormat' section of the",/, &
              "  GMSH grid file, an unknown and invalid section",/, &
              "  header was found before the required '$EndMeshFormat'",/, &
              "  line signifiying the end of the 'MeshFormat' section!",/, &
              "*********************** ERROR! ************************",/)
  4 format (/,"******************** ERROR! *********************",/, &
              "  Error trying to read the 'MeshFormat'",/, &
              "  section of the GMSH grid file!", &
              "       IO Error = ",i0,/, &
              "******************** ERROR! *********************",/)
  !
end subroutine read_gmsh_section_meshformat
!
!###############################################################################
!
subroutine read_gmsh_section_nodes
  !
  !.. Use Statements ..
  use geovar, only : nr,nnode,xyz_nodes
  !
  !.. Local Scalars ..
  integer :: i,j,n,ierr
  logical :: end_of_file
  real(wp) :: coord_ave
  character(len=80) :: text
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable, dimension(:,:) :: temp_nodes
  !
continue
  !
  ! Default the number of grid dimensions to 3D since GMSH grid files
  ! seem to always contain three coordinate values for each node.
  !
  nr = 3
  !
  ! Read the Nodes section of the GMSH grid file
  !
  read (iogrd,*,iostat=ierr,end=100,err=101) nnode
  !
  ! Allocate the temporary array to read in the coordinates of the grid nodes
  !
  allocate ( temp_nodes(1:nr,1:nnode) , source=zero , stat=ierr )
  call alloc_error("read_gmsh_section_nodes","temp_nodes",1,ierr)
  !
  ! Read the coordinates of the grid nodes
  !
  do n = 1,nnode
    read (iogrd,*,iostat=ierr,end=100,err=101) i, (temp_nodes(j,i), j=1,nr)
  end do
  !
  ! If all the nodes coordinates for the third dimensions
  ! are zero, reduce nr to 2 since this is a 2D grid.
  !
  coord_ave = sum( abs(temp_nodes(nr,1:nnode)) ) / real(nnode,kind=wp)
  !
  if (coord_ave <= eps14) nr = 2
  !
  ! Allocate the xyz_nodes array and use the temp_nodes array as the data source
  !
  allocate ( xyz_nodes(1:nr,1:nnode) , source=temp_nodes(1:nr,1:nnode) , stat=ierr )
  call alloc_error("read_gmsh_section_nodes","xyz_nodes",1,ierr)
  !
  ! Deallocate temp_nodes now that we no longer need it
  !
  deallocate ( temp_nodes , stat=ierr )
  call alloc_error("read_gmsh_section_nodes","temp_nodes",2,ierr)
  !
  ! Return if the line marking the end of the Nodes section is found
  !
  call find_gmsh_section(text,end_of_file)
  !
  if (text == "ENDNODES") then
    nodes_successful = true
    return
  else if (end_of_file) then
    if (mypnum == 0) write (iout,1)
  else
    if (mypnum == 0) write (iout,2)
  end if
  !
  call stop_recon(stop_mpi,"read_gmsh_section_nodes")
  !
  ! Read error continuation statements
  !
  100 if (mypnum == 0) write (iout,1) ierr
      call stop_recon(stop_mpi,"read_gmsh_section_nodes")
  !
  101 if (mypnum == 0) write (iout,3) ierr
      call stop_recon(stop_mpi,"read_gmsh_section_nodes")
  !
  ! Format Statements
  !
  1 format (/,"******************** ERROR! *********************",/, &
              "  Pre-mature end-of-file while trying to read",/, &
              "  the 'Nodes' section of the GMSH grid file!",/, &
              "******************** ERROR! *********************",/)
  2 format (/,"*********************** ERROR! ************************",/, &
              "  While trying to read the 'Nodes' section of the",/, &
              "  GMSH grid file, an unknown and invalid section",/, &
              "  header was found before the required '$EndNodes'",/, &
              "  line signifiying the end of the 'Nodes' section!",/, &
              "*********************** ERROR! ************************",/)
  3 format (/,"******************** ERROR! *********************",/, &
              "  Error trying to read the 'Nodes'",/, &
              "  section of the GMSH grid file!", &
              "       IO Error = ",i0,/, &
              "******************** ERROR! *********************",/)
  !
end subroutine read_gmsh_section_nodes
!
!###############################################################################
!
subroutine read_gmsh_section_elements
  !
  !.. Local Scalars ..
  integer :: i,j,m,n,ierr
  integer :: ntype,ntags,nnods,ngeom
  logical :: end_of_file
  character(len=80) :: text,array_name
  !
continue
  !
  ! Read the Elements section of the GMSH grid file
  !
  read (iogrd,*,iostat=ierr,end=100,err=101) num_cell_data
  !
  ! Allocate the array for the element definitions
  !
  allocate ( cell_data(1:num_cell_data) , stat=ierr )
  call alloc_error("read_gmsh_section_elements","cell_data",1,ierr)
  !
  ! Read the element definitions
  !
  do n = 1,num_cell_data
    !
    ! Read the first three values on the line for the current cell
    ! and do not advance to the next record. We need to allocate
    ! some memory based on the second two values and then we will
    ! read the remaining information on this line.
    !
    read (iogrd,*,iostat=ierr,end=100,err=101) i, ntype, ntags
    backspace (iogrd)
    !
    ! Corrections for Hex-64 and Hex-125 cells to change the cell_type
    ! variable to a value applicable to the cell_nodes array.
    !
    if (ntype == 92) then
      ntype = 32 ! Hex-64
    else if (ntype == 93) then
      ntype = 33 ! Hex-125
    end if
    !
    nnods = cell_nodes(ntype)
    ngeom = cell_geom(ntype)
    !
    ! Allocate the tags and nodes arrays within cell_data for the
    ! current cell and store the cell type and number of tags as well
    !
    cell_data(i)%cell_type = ntype
    cell_data(i)%cell_geom = ngeom
    cell_data(i)%num_tags = ntags
    cell_data(i)%num_nodes = nnods
    !
    write (array_name,4) "cell_data(", i, ")%nodes"
    allocate ( cell_data(i)%nodes(1:nnods) , source=0 , stat=ierr )
    call alloc_error("read_gmsh_section_elements",array_name,1,ierr)
    !
    write (array_name,4) "cell_data(", i, ")%tags "
    allocate ( cell_data(i)%tags(1:ntags) , source=0 , stat=ierr )
    call alloc_error("read_gmsh_section_elements",array_name,1,ierr)
    !
    ! Continue reading the remaining information for this cell
    !
    read (iogrd,*,iostat=ierr,end=100,err=101) m, m, m, &
                              (cell_data(i)%tags(j),  j=1,ntags), &
                              (cell_data(i)%nodes(j), j=1,nnods)
    !
  end do
  !
  ! Return if the line marking the end of the Elements section is found
  !
  call find_gmsh_section(text,end_of_file)
  !
  if (text == "ENDELEMENTS") then
    elements_successful = true
    return
  else if (end_of_file) then
    if (mypnum == 0) write (iout,1)
  else
    if (mypnum == 0) write (iout,2)
  end if
  !
  call stop_recon(stop_mpi,"read_gmsh_section_elements")
  !
  ! Read error continuation statements
  !
  100 if (mypnum == 0) write (iout,1) ierr
      call stop_recon(stop_mpi,"read_gmsh_section_elements")
  !
  101 if (mypnum == 0) write (iout,3) ierr
      call stop_recon(stop_mpi,"read_gmsh_section_elements")
  !
  ! Format Statements
  !
  1 format (/,"******************** ERROR! *********************",/, &
              "  Pre-mature end-of-file while trying to read",/, &
              "  the 'Elements' section of the GMSH grid file!",/, &
              "******************** ERROR! *********************",/)
  2 format (/,"*********************** ERROR! ************************",/, &
              "  While trying to read the 'Elements' section of the",/, &
              "  GMSH grid file, an unknown and invalid section",/, &
              "  header was found before the required '$EndElements'",/, &
              "  line signifiying the end of the 'Elements' section!",/, &
              "*********************** ERROR! ************************",/)
  3 format (/,"******************** ERROR! *********************",/, &
              "  Error trying to read the 'Elements'",/, &
              "  section of the GMSH grid file!", &
              "       IO Error = ",i0,/, &
              "******************** ERROR! *********************",/)
  4 format (a,i0,a)
  !
end subroutine read_gmsh_section_elements
!
!###############################################################################
!
subroutine read_gmsh_section_physicalnames
  !
  !.. Local Scalars ..
  integer :: i,n,ierr
  integer :: phys_dim,phys_id
  logical :: end_of_file
  character(len=80) :: text,phys_nam
  !
  !.. Local Parameters ..
  character(len=*), parameter :: delimiters(1:4) = (/"'",'"',",",achar(9)/)
  !
continue
  !
  ! Read the PhysicalNames section of the GMSH grid file
  !
  read (iogrd,*,iostat=ierr,end=100,err=101) num_phys_grps
  !
  ! Allocate the array for the element definitions
  !
  allocate ( phys_grps(1:num_phys_grps) , stat=ierr )
  call alloc_error("read_gmsh_section_physicalnames","phys_grps",1,ierr)
  !
  ! Read the element definitions
  !
  do n = 1,num_phys_grps
    !
    read (iogrd,*,iostat=ierr,end=100,err=101) phys_dim, phys_id, phys_nam
    !
    ! Change all delimiters in phys_name to spaces
    !
    forall ( i=1:len_trim(phys_nam) , any(phys_nam(i:i)==delimiters) ) &
      phys_nam(i:i) = " "
    phys_nam = uppercase( trim(adjustl(phys_nam)) )
    !
    ! Assign these values to the phys_grps array
    !
    phys_grps(n)%phys_nam = phys_nam
    phys_grps(n)%phys_dim = phys_dim
    phys_grps(n)%phys_id  = phys_id
    phys_grps(n)%ibc = bc_string_to_integer( phys_nam )
    !
  end do
  !
  ! Return if the line marking the end of the PhysicalNames section is found
  !
  call find_gmsh_section(text,end_of_file)
  !
  if (text == "ENDPHYSICALNAMES") then
    physicalnames_successful = true
    return
  else if (end_of_file) then
    if (mypnum == 0) write (iout,1)
  else
    if (mypnum == 0) write (iout,2)
  end if
  !
  call stop_recon(stop_mpi,"read_gmsh_section_physicalnames")
  !
  ! Read error continuation statements
  !
  100 if (mypnum == 0) write (iout,1) ierr
      call stop_recon(stop_mpi,"read_gmsh_section_physicalnames")
  !
  101 if (mypnum == 0) write (iout,3) ierr
      call stop_recon(stop_mpi,"read_gmsh_section_physicalnames")
  !
  ! Format Statements
  !
  1 format (/,"******************** ERROR! *********************",/, &
              "  Pre-mature end-of-file while trying to read",/, &
              "  the 'PhysicalNames' section of the GMSH grid file!",/, &
              "******************** ERROR! *********************",/)
  2 format (/,"*********************** ERROR! ************************",/, &
              "  While trying to read the 'PhysicalNames' section of the",/, &
              "  GMSH grid file, an unknown and invalid section",/, &
              "  header was found before the required '$EndPhysicalNames'",/, &
              "  line signifiying the end of the 'PhysicalNames' section!",/, &
              "*********************** ERROR! ************************",/)
  3 format (/,"******************** ERROR! *********************",/, &
              "  Error trying to read the 'PhysicalNames'",/, &
              "  section of the GMSH grid file!", &
              "       IO Error = ",i0,/, &
              "******************** ERROR! *********************",/)
  !
end subroutine read_gmsh_section_physicalnames
!
!###############################################################################
!
subroutine create_gmsh_grid_arrays
  !
  ! Procedure to create the remaining grid geometry arrays using
  ! the information stored in the arrays cell_data and phys_grps.
  !
  !.. Use Statements ..
  use geovar, only : nr,ncell,nnode,nfbnd,nbfai,bface
  use geovar, only : nodes_of_cell_ptr,nodes_of_cell
  use ovar,   only : iop_solpts,iop_flxpts
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,ifa,nf,nh,nb,h1,h2,b1,b2,i,j
  integer :: np,nc,nfp,igrps,bgrps,pid,pdm,inum,ierr
  logical :: error_found
  character(len=80) :: array_name
  !
  !.. Local Arrays ..
  integer, dimension(1:4) :: lnode
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:) :: int_grps
  integer, allocatable, dimension(:) :: bnd_grps
  logical, allocatable, dimension(:) :: mskp
  logical, allocatable, dimension(:) :: mske
  integer, allocatable, dimension(:) :: lpoin
  logical, allocatable, dimension(:) :: bc_done
  !
continue
  !
  error_found = fals
  !
  ! Initialize the number of interior cells and boundary faces to zero
  !
  ncell = 0
  nfbnd = 0
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Count the number of and collect the interior and boundary physical groups
  !
  igrps = count(phys_grps(:)%phys_dim == nr)
  bgrps = count(phys_grps(:)%phys_dim == nr-1)
  !
  allocate ( int_grps(1:igrps) , source=0 , stat=ierr)
  call alloc_error("create_gmsh_grid_arrays","int_grps",1,ierr)
  allocate ( bnd_grps(1:bgrps) , source=0 , stat=ierr)
  call alloc_error("create_gmsh_grid_arrays","bnd_grps",1,ierr)
  !
  int_grps(1:igrps) = pack( [1:num_phys_grps] , phys_grps(:)%phys_dim == nr )
  bnd_grps(1:bgrps) = pack( [1:num_phys_grps] , phys_grps(:)%phys_dim == nr-1 )
  !
  ! Check that the total number of interior and boundary groups is correct
  !
  if ( igrps+bgrps /= num_phys_grps ) then
    if (mypnum == 0) write (iout,1) igrps,bgrps,num_phys_grps
    call stop_recon(stop_mpi,"create_gmsh_grid_arrays")
  end if
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Loop through the physical groups and collect
  ! the elements for each physical group
  !
  do i = 1,num_phys_grps
    !
    pid = phys_grps(i)%phys_id
    pdm = phys_grps(i)%phys_dim
    !
    ! Count the number of elements that match this physical groups id number
    !
    inum = 0
    do j = 1,size(cell_data)
      if (allocated(cell_data(j)%tags)) then
        inum = inum + merge( 1 , 0 , cell_data(j)%tags(1) == pid )
      end if
    end do
   !inum = count( cell_data(:)%tags(1) == pid )
    !
    ! Allocate the cells array for this physical group
    !
    write (array_name,3) "phys_grps(", i, ")%cells "
    allocate ( phys_grps(i)%cells(1:inum) , source=0 , stat=ierr )
    call alloc_error("create_gmsh_grid_arrays",array_name,1,ierr)
    !
    ! Create the list of elements that are in this physical group
    !
    inum = 0
    do j = 1,size(cell_data)
      if (allocated(cell_data(j)%tags)) then
        if (cell_data(j)%tags(1) == pid) then
          inum = inum + 1
          phys_grps(i)%cells(inum) = j
        end if
      end if
    end do
   !phys_grps(i)%cells(1:inum) = pack( [1:num_cell_data] , &
   !                                   cell_data(:)%tags(1) == pid )
    !
  end do
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Count the total number of interior cells and boundary faces
  !
  ncell = 0
  do i = 1,igrps
    ncell = ncell + size( phys_grps(int_grps(i))%cells )
  end do
  write (iout,*) "ncell = ",ncell
  !
  nfbnd = 0
  do i = 1,bgrps
    nfbnd = nfbnd + size( phys_grps(bnd_grps(i))%cells )
  end do
  write (iout,*) "nfbnd = ",nfbnd
  !
  ! Check that the number of interior cells and boundary faces are correct
  !
  if ( ncell+nfbnd /= num_cell_data ) then
    if (mypnum == 0) write (iout,2) ncell,nfbnd,num_cell_data
    call stop_recon(stop_mpi,"create_gmsh_grid_arrays")
  end if
  !
  ! #######
  ! ## 4 ##
  ! #######
  !
  ! Allocate the nodes_of_cell_ptr array
  !
  allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , stat=ierr )
  call alloc_error("create_gmsh_grid_arrays","nodes_of_cell_ptr",1,ierr)
  !
  ! Loop through the cells for each interior physical
  ! group and use cell_data to create nodes_of_cell_ptr
  !
  do i = 1,igrps
    np = int_grps(i)
    do j = 1,size(phys_grps(np)%cells)
      n = phys_grps(np)%cells(j)
      nodes_of_cell_ptr(n+1) = cell_data(n)%num_nodes
    end do
  end do
  !
  ! Reshuffle nodes_of_cell_ptr
  !
  do n = 2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(n) = nodes_of_cell_ptr(n) + nodes_of_cell_ptr(n-1)
  end do
  !
  ! #######
  ! ## 5 ##
  ! #######
  !
  ! Allocate the nodes_of_cell array
  !
  allocate ( nodes_of_cell(1:last(nodes_of_cell_ptr)) , source=0 , stat=ierr )
  call alloc_error("create_gmsh_grid_arrays","nodes_of_cell",1,ierr)
  !
  ! Loop over the interior cells and fill in the node connectivity
  ! for that cell from the information in cell_data
  !
  do n = 1,ncell
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    nodes_of_cell(n1:n2) = cell_data(n)%nodes(:)
  end do
  !
  ! #######
  ! ## 6 ##
  ! #######
  !
  ! Before we create the bface array, we need to
  ! find the host cells for each boundary face
  !
  ! Allocate mask arrays for all the grid cells and nodes
  !
  allocate ( mskp(1:nnode) , source=fals , stat=ierr )
  call alloc_error("create_gmsh_grid_arrays","mskp",1,ierr)
  allocate ( mske(1:ncell) , source=fals , stat=ierr )
  call alloc_error("create_gmsh_grid_arrays","mske",1,ierr)
  !
  ! Identify the nodes that are on boundary faces
  !
  do i = 1,bgrps
    np = bnd_grps(i)
    do j = 1,size(phys_grps(np)%cells)
      n = phys_grps(np)%cells(j)
      mskp( cell_data(n)%nodes(:) ) = true
    end do
  end do
  !
  ! Identify all cells that contain a boundary node
  !
  do i = 1,igrps
    np = int_grps(i)
    do j = 1,size(phys_grps(np)%cells)
      n = phys_grps(np)%cells(j)
      mske(n) = any( mskp(cell_data(n)%nodes(:)) )
    end do
  end do
  !
  ! Deallocate mskp and int_grps since we dont need them anymore
  !
  deallocate ( mskp , stat=ierr )
  call alloc_error("create_gmsh_grid_arrays","mskp",2,ierr)
  deallocate ( int_grps , stat=ierr )
  call alloc_error("create_gmsh_grid_arrays","int_grps",2,ierr)
  !
  ! Allocate the point counting array, lpoin, and the host_cell array
  !
  allocate ( lpoin(1:nnode) , source=0 , stat=ierr )
  call alloc_error("create_gmsh_grid_arrays","lpoin",1,ierr)
  !
  ! Find the host cells for each boundary face
  !
  do i = 1,bgrps
    !
    np = bnd_grps(i)
    !
    face_loop: do j = 1,size(phys_grps(np)%cells)
      !
      n = phys_grps(np)%cells(j)
      nfp = cell_data(n)%num_nodes
      !
      ! Mark the nodes for this boundary face
      !
      lpoin(cell_data(n)%nodes(:)) = 1
      !
      ! Loop over the cells containing boundary nodes and find
      ! the cell with the same number of marked nodes as the number
      ! of nodes defining this boundary face
      !
      do nc = 1,ncell
        if (.not. mske(nc)) cycle
        n1 = nodes_of_cell_ptr(nc)+1
        n2 = nodes_of_cell_ptr(nc+1)
        if ( sum(lpoin(nodes_of_cell(n1:n2))) == nfp) then
          cell_data(n)%host_cell = nc
          lpoin(cell_data(n)%nodes(:)) = 0
          cycle face_loop
        end if
      end do
      !
    end do face_loop
    !
  end do
  !
  ! Deallocate mske and lpoin since we no longer need them
  !
  deallocate ( mske , stat=ierr )
  call alloc_error("create_gmsh_grid_arrays","mske",2,ierr)
  deallocate ( lpoin , stat=ierr )
  call alloc_error("create_gmsh_grid_arrays","lpoin",2,ierr)
  !
  ! #######
  ! ## 6 ##
  ! #######
  !
  ! Create and fill out bface with data from the
  ! PhysicalName section of the GMSH grid file
  !
  ! Make sure nbfai is large enough for high-order boundary faces
  !
  nbfai = 10
  do i = 1,bgrps
    np = bnd_grps(i)
    do j = 1,size(phys_grps(np)%cells)
      n = phys_grps(np)%cells(j)
      nbfai = max( nbfai , cell_data(n)%num_nodes+10 )
    end do
  end do
  !
  ! Allocate the bface array
  !
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , stat=ierr )
  call alloc_error("create_gmsh_grid_arrays","bface",1,ierr)
  !
  ! Loop through the faces for each boundary physical
  ! group and use cell_data to create bface
  !
  nf = 0
  boundary_phys_grps: do i = 1,bgrps
    !
    np = bnd_grps(i)
    !
    current_phys_grp: do j = 1,size(phys_grps(np)%cells)
      !
      nb = phys_grps(np)%cells(j)
      nh = cell_data(nb)%host_cell
      nf = nf + 1
      !
      bface( 1,nf) = phys_grps(np)%ibc       ! face boundary condition
      bface( 2,nf) = nh                      ! host cell of boundary face
      bface( 3,nf) = cell_data(nh)%cell_geom ! host cell geometry
      bface( 5,nf) = i                       ! group for boundary face
      bface(10,nf) = ncell + nf              ! ghost cell of boundary face
      !
      ! Find which face of the host cell is on the boundary
      !
     !b1 = nodes_of_cell_ptr(nb)+1 ! beginning index for boundary element nodes
     !b2 = nodes_of_cell_ptr(nb+1) ! ending    index for boundary element nodes
      !
      h1 = nodes_of_cell_ptr(nh)+1 ! beginning index for host cell nodes
      h2 = nodes_of_cell_ptr(nh+1) ! ending    index for host cell nodes
      !
      ! Find which face of the host cell is on the boundary
      ! depending on the element type of the host cell
      !
      ifa = 0
      select case (bface(3,nf))
        case (0)
          ifa = -1
        case (1)
          ifa = -2
        case (2)
         !ifa = host_face_quad( nodes_of_cell(h1:h2) , nodes_of_cell(b1:b2) )
          ifa = host_face_quad( nodes_of_cell(h1:h2) , cell_data(nb)%nodes(:) )
        case (3)
         !ifa = host_face_tria( nodes_of_cell(h1:h2) , nodes_of_cell(b1:b2) )
          ifa = host_face_tria( nodes_of_cell(h1:h2) , cell_data(nb)%nodes(:) )
        case (4)
         !ifa = host_face_tetr( nodes_of_cell(h1:h2) , nodes_of_cell(b1:b2) )
          ifa = host_face_tetr( nodes_of_cell(h1:h2) , cell_data(nb)%nodes(:) )
        case (5)
         !ifa = host_face_pyra( nodes_of_cell(h1:h2) , nodes_of_cell(b1:b2) )
          ifa = host_face_pyra( nodes_of_cell(h1:h2) , cell_data(nb)%nodes(:) )
        case (6)
         !ifa = host_face_pris( nodes_of_cell(h1:h2) , nodes_of_cell(b1:b2) )
          ifa = host_face_pris( nodes_of_cell(h1:h2) , cell_data(nb)%nodes(:) )
        case (8)
         !ifa = host_face_hexa( nodes_of_cell(h1:h2) , nodes_of_cell(b1:b2) )
          ifa = host_face_hexa( nodes_of_cell(h1:h2) , cell_data(nb)%nodes(:) )
        case default
          ifa = -10
      end select
      !
      ! Check for an error finding which face of the host cell is on the
      ! boundary, and if there is no error, store the value in bface
      !
      if (ifa < 0) then
        !
        if (.not. error_found) then
          if (mypnum == 0) write (iout,4)
          error_found = true
        end if
        !
        ! Write out the specific error type corresponding 
        ! to the negative value of ifa
        !
        if (mypnum == 0) then
          if (ifa == -1) then
            write (iout,5) "A 0-D 'Point'-type",nh,nb
          else if (ifa == -2) then
            write (iout,5) "A 1-D 'Edge'-type",nh,nb
          else if (ifa == -4) then
            write (iout,5) "A 3-D tetrahedron",nh,nb
          else if (ifa == -5) then
            write (iout,5) "A 3-D pyramid",nh,nb
          else if (ifa == -6) then
            write (iout,5) "A 3-D prism",nh,nb
          else if (ifa == -8) then
            write (iout,5) "A 3-D hexahedron",nh,nb
          else if (ifa == -10) then
            write (iout,5) "An unknown type of",nh,nb
          else if (ifa == -99) then
            write (iout,6) nh,nb
          end if
        end if
        !
      else
        !
        ! Store the host cell face that is on the boundary
        !
        bface(4,nf) = ifa
        !
      end if
      !
    end do current_phys_grp
    !
  end do boundary_phys_grps
  !
  ! Deallocate bnd_grps since we dont need it anymore
  !
  deallocate ( bnd_grps , stat=ierr )
  call alloc_error("create_gmsh_grid_arrays","bnd_grps",2,ierr)
  !
  ! #######
  ! ## 7 ##
  ! #######
  !
  ! Finish any remaining tasks for the grid geometry
  !
  ! ###############################################################
  ! ###############################################################
  ! ####################   !!! WARNING !!!   ######################
  ! ###############################################################
  ! ###############################################################
  ! ###                                                         ### 
  ! ###  PERIODIC BOUNDARY CONDITIONS ARE NOT CURRENTLY SET UP  ###
  ! ###  WHEN USING GMSH GRID FILES! GAMBIT NEUTRAL GRID FILES  ###
  ! ###    ARE REQUIRED TO USE PERIODIC BOUNDARY CONDITIONS!    ###
  ! ###                                                         ### 
  ! ###############################################################
  ! ###############################################################
  ! ####################   !!! WARNING !!!   ######################
  ! ###############################################################
  ! ###############################################################
  !
  ! If an error was found, kill the program
  !
  if (error_found) then
    if (mypnum == 0) write (iout,4)
    call stop_recon(stop_mpi,"create_gmsh_grid_arrays")
  end if
  !
  ! Make sure that Lobatto nodes are requested if the grid contains triangle
  ! elements since Gauss nodes is not currently set up for triangles.
  !
  do n = 1,ncell
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    if (n2-n1+1 == 3) then
      if (mypnum == 0) write (iout,8)
      call stop_recon(stop_mpi,"create_gmsh_grid_arrays")
    end if
  end do
  if (iop_solpts /= Legendre_Gauss_Lobatto .or. &
      iop_flxpts /= Legendre_Gauss_Lobatto) then
    do n = 1,ncell
      n1 = nodes_of_cell_ptr(n)+1
      n2 = nodes_of_cell_ptr(n+1)
      if (n2-n1+1 == 3) then
        if (mypnum == 0) write (iout,7)
        call stop_recon(stop_mpi,"create_gmsh_grid_arrays")
      end if
    end do
  end if
  !
  ! Format Statements
  !
  1 format (/,"************************ ERROR! **************************",/, &
              "  The number of interior and boundary physical groups is",/, &
              "  inconsistent with the total number of physical groups!",/, &
              "     number of interior groups = ",i0,/, &
              "     number of boundary groups = ",i0,/, &
              "        total number of groups = ",i0,/, &
              "************************ ERROR! **************************",/)
  2 format (/,"************************ ERROR! ************************",/, &
              "  The number of interior cells and boundary faces is",/, &
              "  inconsistent with the total number of grid elements!",/, &
              "          number of interior cells = ",i0,/, &
              "          number of boundary faces = ",i0,/, &
              "     total number of grid elements = ",i0,/, &
              "************************ ERROR! ************************",/)
  3 format (a,i0,a)
  4 format (/,"**********************************************************",/, &
              "************************ ERROR! **************************",/, &
              "**********************************************************",/)
  5 format ("  ",a," host element was found on a boundary!",/, &
            "     Host Cell = ",i0,";    Boundary Face = ",i0,/, &
            "  This element type is not currently supported!",/)
  6 format ("  Unable to find which face of host cell ",i0,/, &
            "  that makes up the boundary face ",i0," !",/)
  7 format (/," Quadrature points besides Gauss-Lobatto nodes are not",/, &
              " currently supported for a grid containing triangle cells.",//)
  8 format (///," Triangles are temporarily disabled due to significant",/, &
                " changes in the quadrilateral part of the code. These ",/, &
                " changes have yet to be applied to triangles.",///)
  !
end subroutine create_gmsh_grid_arrays
!
!###############################################################################
!###############################################################################
!###############################################################################
!
! The remaining GMSH sections have not been set up for use
!
!###############################################################################
!###############################################################################
!###############################################################################
!
subroutine read_gmsh_section_periodic
  !
  !.. Local Scalars ..
  logical :: end_of_file
  character(len=80) :: text
  !
continue
  !
  ! The code is not currently set up to read the gmsh
  ! Periodic section. Continue reading the grid file
  ! until the end of the Periodic section is found.
  !
  call find_gmsh_section(text,end_of_file)
  !
  if (text == "ENDPERIODIC") then
    return
  else if (end_of_file) then
    if (mypnum == 0) write (iout,1)
  else
    if (mypnum == 0) write (iout,2)
  end if
  !
  call stop_recon(stop_mpi,"read_gmsh_section_periodic")
  !
  ! Format Statements
  !
  1 format (/,"******************** ERROR! *********************",/, &
              "  Pre-mature end-of-file while trying to read",/, &
              "  the 'Periodic' section of the GMSH grid file!",/, &
              "******************** ERROR! *********************",/)
  2 format (/,"*********************** ERROR! ************************",/, &
              "  While trying to read the 'Periodic' section of the",/, &
              "  GMSH grid file, an unknown and invalid section",/, &
              "  header was found before the required '$EndPeriodic'",/, &
              "  line signifiying the end of the 'Periodic' section!",/, &
              "*********************** ERROR! ************************",/)
  !
end subroutine read_gmsh_section_periodic
!
!###############################################################################
!
subroutine read_gmsh_section_nodedata
  !
  !.. Local Scalars ..
  logical :: end_of_file
  character(len=80) :: text
  !
continue
  !
  ! The code is not currently set up to read the gmsh
  ! NodeData section. Continue reading the grid file
  ! until the end of the NodeData section is found.
  !
  call find_gmsh_section(text,end_of_file)
  !
  if (text == "ENDNODEDATA") then
    return
  else if (end_of_file) then
    if (mypnum == 0) write (iout,1)
  else
    if (mypnum == 0) write (iout,2)
  end if
  !
  call stop_recon(stop_mpi,"read_gmsh_section_nodedata")
  !
  ! Format Statements
  !
  1 format (/,"******************** ERROR! *********************",/, &
              "  Pre-mature end-of-file while trying to read",/, &
              "  the 'NodeData' section of the GMSH grid file!",/, &
              "******************** ERROR! *********************",/)
  2 format (/,"*********************** ERROR! ************************",/, &
              "  While trying to read the 'NodeData' section of the",/, &
              "  GMSH grid file, an unknown and invalid section",/, &
              "  header was found before the required '$EndNodeData'",/, &
              "  line signifiying the end of the 'NodeData' section!",/, &
              "*********************** ERROR! ************************",/)
  !
end subroutine read_gmsh_section_nodedata
!
!###############################################################################
!
subroutine read_gmsh_section_elementdata
  !
  !.. Local Scalars ..
  logical :: end_of_file
  character(len=80) :: text
  !
continue
  !
  ! The code is not currently set up to read the gmsh
  ! ElementData section. Continue reading the grid file
  ! until the end of the ElementData section is found.
  !
  call find_gmsh_section(text,end_of_file)
  !
  if (text == "ENDELEMENTDATA") then
    return
  else if (end_of_file) then
    if (mypnum == 0) write (iout,1)
  else
    if (mypnum == 0) write (iout,2)
  end if
  !
  call stop_recon(stop_mpi,"read_gmsh_section_elementdata")
  !
  ! Format Statements
  !
  1 format (/,"******************** ERROR! *********************",/, &
              "  Pre-mature end-of-file while trying to read",/, &
              "  the 'ElementData' section of the GMSH grid file!",/, &
              "******************** ERROR! *********************",/)
  2 format (/,"*********************** ERROR! ************************",/, &
              "  While trying to read the 'ElementData' section of the",/, &
              "  GMSH grid file, an unknown and invalid section",/, &
              "  header was found before the required '$EndElementData'",/, &
              "  line signifiying the end of the 'ElementData' section!",/, &
              "*********************** ERROR! ************************",/)
  !
end subroutine read_gmsh_section_elementdata
!
!###############################################################################
!
subroutine read_gmsh_section_elementnodedata
  !
  !.. Local Scalars ..
  logical :: end_of_file
  character(len=80) :: text
  !
continue
  !
  ! The code is not currently set up to read the gmsh
  ! ElementNodeData section. Continue reading the grid file
  ! until the end of the ElementNodeData section is found.
  !
  call find_gmsh_section(text,end_of_file)
  !
  if (text == "ENDELEMENTNODEDATA") then
    return
  else if (end_of_file) then
    if (mypnum == 0) write (iout,1)
  else
    if (mypnum == 0) write (iout,2)
  end if
  !
  call stop_recon(stop_mpi,"read_gmsh_section_elementnodedata")
  !
  ! Format Statements
  !
  1 format (/,"******************** ERROR! *********************",/, &
              "  Pre-mature end-of-file while trying to read",/, &
              "  the 'ElementNodeData' section of the GMSH grid file!",/, &
              "******************** ERROR! *********************",/)
  2 format (/,"*********************** ERROR! ************************",/, &
              "  While trying to read the 'ElementNodeData' section of the",/, &
              "  GMSH grid file, an unknown and invalid section",/, &
              "  header was found before the required '$EndElementNodeData'",/, &
              "  line signifiying the end of the 'ElementNodeData' section!",/, &
              "*********************** ERROR! ************************",/)
  !
end subroutine read_gmsh_section_elementnodedata
!
!###############################################################################
!
subroutine read_gmsh_section_interpolationscheme
  !
  !.. Local Scalars ..
  logical :: end_of_file
  character(len=80) :: text
  !
continue
  !
  ! The code is not currently set up to read the gmsh
  ! InterpolationScheme section. Continue reading the grid file
  ! until the end of the InterpolationScheme section is found.
  !
  call find_gmsh_section(text,end_of_file)
  !
  if (text == "ENDINTERPOLATIONSCHEME") then
    return
  else if (end_of_file) then
    if (mypnum == 0) write (iout,1)
  else
    if (mypnum == 0) write (iout,2)
  end if
  !
  call stop_recon(stop_mpi,"read_gmsh_section_interpolationscheme")
  !
  ! Format Statements
  !
  1 format (/,"******************** ERROR! *********************",/, &
              "  Pre-mature end-of-file while trying to read",/, &
              "  the 'InterpolationScheme' section of the GMSH grid file!",/, &
              "******************** ERROR! *********************",/)
  2 format (/,"*********************** ERROR! ************************",/, &
              "  While trying to read the 'InterpolationScheme' section of the",/, &
              "  GMSH grid file, an unknown and invalid section",/, &
              "  header was found before the required '$EndInterpolationScheme'",/, &
              "  line signifiying the end of the 'InterpolationScheme' section!",/, &
              "*********************** ERROR! ************************",/)
  !
end subroutine read_gmsh_section_interpolationscheme
!
!###############################################################################
!
pure function host_face_quad(int_nodes,bnd_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: int_nodes
  integer, dimension(:), intent(in) :: bnd_nodes
  !
  !.. Function Return Value
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Arrays ..
  integer, dimension(1:4) :: marker
  !
  !.. Local Parameters ..
  integer, parameter :: num_face_nodes = 2
  integer, parameter :: fnodes(1:2,1:4) = reshape( [ 1, 2, &
                                                     2, 3, &
                                                     3, 4, &
                                                     4, 1 ] , [2,4] )
  !
continue
  !
  ! Initialize return_value with the face not found value of -99
  !
  return_value = -99
  !
  ! Initialize the marker array to zero
  !
  marker(:) = 0
  !
  ! Mark each host cell node that is on the boundary
  !
  do i = 1,num_face_nodes
    where (int_nodes(1:4) == bnd_nodes(i)) marker = 1
  end do
  !
  ! Now find which face corresponds to the marked nodes
  !
  do i = 1,size(fnodes,dim=2)
    if ( sum( marker(fnodes(:,i)) ) == num_face_nodes ) then
      return_value = i
      exit
    end if
  end do
  !
end function host_face_quad
!
!###############################################################################
!
pure function host_face_tria(int_nodes,bnd_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: int_nodes
  integer, dimension(:), intent(in) :: bnd_nodes
  !
  !.. Function Return Value
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Arrays ..
  integer, dimension(1:3) :: marker
  !
  !.. Local Parameters ..
  integer, parameter :: num_face_nodes = 2
  integer, parameter :: fnodes(1:2,1:3) = reshape( [ 1, 2, &
                                                     2, 3, &
                                                     3, 1 ] , [2,3] )
  !
continue
  !
  ! Initialize return_value with the face not found value of -99
  !
  return_value = -99
  !
  ! Initialize the marker array to zero
  !
  marker(:) = 0
  !
  ! Mark each host cell node that is on the boundary
  !
  do i = 1,num_face_nodes
    where (int_nodes(1:3) == bnd_nodes(i)) marker = 1
  end do
  !
  ! Now find which face corresponds to the marked nodes
  !
  do i = 1,size(fnodes,dim=2)
    if ( sum( marker(fnodes(:,i)) ) == num_face_nodes ) then
      return_value = i
      exit
    end if
  end do
  !
end function host_face_tria
!
!###############################################################################
!
pure function host_face_tetr(int_nodes,bnd_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: int_nodes
  integer, dimension(:), intent(in) :: bnd_nodes
  !
  !.. Function Return Value
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Arrays ..
  integer, dimension(1:4) :: marker
  !
  !.. Local Parameters ..
  integer, parameter :: num_face_nodes = 3
  integer, parameter :: fnodes(1:3,1:4) = reshape( [ 2, 3, 4, &
                                                     3, 1, 4, &
                                                     4, 1, 2, &
                                                     1, 3, 2 ] , [3,4] )
  !
continue
  !
  ! Initialize return_value with the face not found value of -99
  !
  return_value = -99
  !
  ! Initialize the marker array to zero
  !
  marker(:) = 0
  !
  ! Mark each host cell node that is on the boundary
  !
  do i = 1,num_face_nodes
    where (int_nodes(1:4) == bnd_nodes(i)) marker = 1
  end do
  !
  ! Now find which face corresponds to the marked nodes
  !
  do i = 1,size(fnodes,dim=2)
    if ( sum( marker(fnodes(:,i)) ) == num_face_nodes ) then
      return_value = i
      exit
    end if
  end do
  !
  ! Tetrahedron elements are not currently set up so set
  ! return_value to the currently not set up error value of -4
  !
  return_value = -4
  !
end function host_face_tetr
!
!###############################################################################
!
pure function host_face_pyra(int_nodes,bnd_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: int_nodes
  integer, dimension(:), intent(in) :: bnd_nodes
  !
  !.. Function Return Value
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Arrays ..
  integer, dimension(1:5) :: marker
  !
  !.. Local Parameters ..
  integer, parameter :: num_face_nodes = 3
  integer, parameter :: fnodes(1:4,1:5) = reshape( [ 2, 3, 5, 0, &
                                                     3, 4, 5, 0, &
                                                     4, 1, 5, 0, &
                                                     1, 2, 5, 0, &
                                                     1, 4, 3, 2 ] , [4,5] )
  !
continue
  !
  ! Initialize return_value with the face not found value of -99
  !
  return_value = -99
  !
  ! Initialize the marker array to zero
  !
  marker(:) = 0
  !
  ! Mark each host cell node that is on the boundary
  ! NOTE: We only check the first three boundary nodes because those three nodes
  !       are enough to identify which of the five faces is on the boundary.
  !       This way we dont have to differentiate between the quad and triangle
  !       faces of the pyramid element.
  !
  do i = 1,num_face_nodes
    where (int_nodes(1:5) == bnd_nodes(i)) marker = 1
  end do
  !
  ! Now find which face corresponds to the marked nodes
  !
  do i = 1,size(fnodes,dim=2)
    if ( sum( marker(fnodes(:,i)) ) == num_face_nodes ) then
      return_value = i
      exit
    end if
  end do
  !
  ! Pyramid elements are not currently set up so set
  ! return_value to the currently not set up error value of -5
  !
  return_value = -5
  !
end function host_face_pyra
!
!###############################################################################
!
pure function host_face_pris(int_nodes,bnd_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: int_nodes
  integer, dimension(:), intent(in) :: bnd_nodes
  !
  !.. Function Return Value
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Arrays ..
  integer, dimension(1:6) :: marker
  !
  !.. Local Parameters ..
  integer, parameter :: num_face_nodes = 3
  integer, parameter :: fnodes(1:4,1:5) = reshape( [ 1, 3, 2, 0, &
                                                     4, 5, 6, 0, &
                                                     1, 2, 5, 4, &
                                                     2, 3, 6, 5, &
                                                     4, 6, 3, 1 ] , [4,5] )
  !
continue
  !
  ! Initialize return_value with the face not found value of -99
  !
  return_value = -99
  !
  ! Initialize the marker array to zero
  !
  marker(:) = 0
  !
  ! Mark each host cell node that is on the boundary
  ! NOTE: We only check the first three boundary nodes because those three nodes
  !       are enough to identify which of the five faces is on the boundary.
  !       This way we dont have to differentiate between the quad and triangle
  !       faces of the prism element.
  !
  do i = 1,num_face_nodes
    where (int_nodes(1:6) == bnd_nodes(i)) marker = 1
  end do
  !
  ! Now find which face corresponds to the marked nodes
  !
  do i = 1,size(fnodes,dim=2)
    if ( sum( marker(fnodes(:,i)) ) == num_face_nodes ) then
      return_value = i
      exit
    end if
  end do
  !
  ! Prism elements are not currently set up so set
  ! return_value to the currently not set up error value of -6
  !
  return_value = -6
  !
end function host_face_pris
!
!###############################################################################
!
pure function host_face_hexa(int_nodes,bnd_nodes) result(return_value)
  !
  !.. Formal Arguments ..
  integer, dimension(:), intent(in) :: int_nodes
  integer, dimension(:), intent(in) :: bnd_nodes
  !
  !.. Function Return Value
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Arrays ..
  integer, dimension(1:8) :: marker
  !
  !.. Local Parameters ..
  integer, parameter :: num_face_nodes = 4
  integer, parameter :: fnodes(1:4,1:6) = reshape( [ 1, 4, 3, 2, &
                                                     5, 6, 7, 8, &
                                                     2, 3, 7, 6, &
                                                     5, 8, 4, 1, &
                                                     1, 2, 6, 5, &
                                                     3, 4, 8, 7 ] , [4,6] )
  !
continue
  !
  ! Initialize return_value with the face not found value of -99
  !
  return_value = -99
  !
  ! Initialize the marker array to zero
  !
  marker(:) = 0
  !
  ! Mark each host cell node that is on the boundary
  !
  do i = 1,num_face_nodes
    where (int_nodes(1:8) == bnd_nodes(i)) marker = 1
  end do
  !
  ! Now find which face corresponds to the marked nodes
  !
  do i = 1,size(fnodes,dim=2)
    if ( sum( marker(fnodes(:,i)) ) == num_face_nodes ) then
      return_value = i
      exit
    end if
  end do
  !
  ! Hexahedron elements are not currently set up so set
  ! return_value to the currently not set up error value of -8
  !
  return_value = -8
  !
end function host_face_hexa
!
!###############################################################################
!
include "Functions/bc_string_to_integer.f90"
include "Functions/last.f90"
!
!###############################################################################
!
end module module_gmsh
!
!###############################################################################
!
