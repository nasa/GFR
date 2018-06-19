module module_gmsh
  !
  !.. Use Statements ..
  use module_kind_types
  use geovar, only : grid
  use geovar, only : elem_t
  use geovar, only : prop_t
  !
  implicit none
  !
  private
  !
  ! ##################################
  ! ###  Module Public Procedures  ###
  ! ##################################
  !
  ! Public Interface for Gmsh element_properties function
  !
  interface gmsh_element_properties
    module procedure element_properties
  end interface gmsh_element_properties
  !
  public :: read_gmsh_gridfile
  public :: gmsh_element_properties
  public :: gmsh_memory_usage
  !
  ! ###########################
  ! ###  Module Parameters  ###
  ! ###########################
  !
  ! Explanation of different element types created by Gmsh
  !
  !   1 : 2-node 1st order line
  !   2 : 3-node 1st order triangle
  !   3 : 4-node 1st order quadrangle
  !   4 : 4-node 1st order tetrahedron
  !   5 : 8-node 1st order hexahedron
  !   6 : 6-node 1st order prism
  !   7 : 5-node 1st order pyramid
  !   8 : 3-node 2nd order line
  !            2 nodes associated with the vertices
  !            1 with the edge
  !   9 : 6-node 2nd order triangle
  !            3 nodes associated with the vertices
  !            3 with the edges
  !  10 : 9-node 2nd order quadrangle
  !            4 nodes associated with the vertices
  !            4 with the edges
  !            1 with the interior
  !  11 : 10-node 2nd order tetrahedron
  !            4 nodes associated with the vertices
  !            6 with the edges
  !  12 : 27-node 2nd order hexahedron
  !            8 nodes associated with the vertices
  !           12 with the edges
  !            6 with the faces
  !            1 with the volume
  !  13 : 18-node 2nd order prism
  !            6 nodes associated with the vertices
  !            9 with the edges
  !            3 with the quadrangular faces
  !  14 : 14-node 2nd order pyramid
  !            5 nodes associated with the vertices
  !            8 with the edges
  !            1 with the quadrangular face
  !  15 : 1-node 1st order point
  !  16 : 8-node 2nd order EDGE-DEFINED quadrangle
  !            4 nodes associated with the vertices
  !            4 with the edges
  !  17 : 20-node 2nd order EDGE-DEFINED hexahedron
  !            8 nodes associated with the vertices
  !           12 with the edges
  !  18 : 15-node 2nd order EDGE-DEFINED prism
  !            6 nodes associated with the vertices
  !            9 with the edges
  !  19 : 13-node 2nd order EDGE-DEFINED pyramid
  !            5 nodes associated with the vertices
  !            8 with the edges
  !  20 : 9-node 3rd order EDGE-DEFINED triangle
  !            3 nodes associated with the vertices
  !            6 with the edges
  !  21 : 10-node 3rd order triangle
  !            3 nodes associated with the vertices
  !            6 with the edges
  !            1 with the interior
  !  22 : 12-node 4th order EDGE-DEFINED triangle
  !            3 nodes associated with the vertices
  !            9 with the edges
  !  23 : 15-node 4th order triangle
  !            3 nodes associated with the vertices
  !            9 with the edges
  !            3 with the interior
  !  24 : 15-node 5th order EDGE-DEFINED triangle
  !            3 nodes associated with the vertices
  !           12 with the edges
  !  25 : 21-node 5th order triangle
  !            3 nodes associated with the vertices
  !           12 with the edges
  !            6 with the interior
  !  26 : 4-node 3rd order edge
  !            2 nodes associated with the vertices
  !            2 internal to the edge
  !  27 : 5-node 4th order edge
  !            2 nodes associated with the vertices
  !            3 internal to the edge
  !  28 : 6-node 5th order edge
  !            2 nodes associated with the vertices
  !            4 internal to the edge
  !  29 : 20-node 3rd order tetrahedron
  !            4 nodes associated with the vertices
  !           12 with the edges
  !            4 with the faces
  !  30 : 35-node 4th order tetrahedron
  !            4 nodes associated with the vertices
  !           18 with the edges
  !           12 with the faces
  !            1 in the volume
  !  31 : 56-node 5th order tetrahedron
  !            4 nodes associated with the vertices
  !           24 with the edges
  !           24 with the faces
  !            4 in the volume
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###                                                                         ###
!###               NEWLY ADDED ELEMENTS FROM GMSH SOURCE CODE                ###
!###       THAT DO NOT SEEM TO HAVE ANY DOCUMENTATION RELATING TO THEM       ###
!###                                                                         ###
!###  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  ###
!###  VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV  ###
!###                                                                         ###
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
  !  32 : 22-node  4th order EDGE-DEFINED tetrahedron
  !            4 nodes associated with the vertices
  !           18 with the edges
  !  33 : 28-node  5th order EDGE-DEFINED tetrahedron
  !            4 nodes associated with the vertices
  !           24 with the edges
  !  34 : ?-node  1st order polygon
  !           * nodes associated with the vertices
  !           * with the edges
  !           * with the faces
  !           * in the volume
  !  35 : ?-node  1st order polyhedron
  !           * nodes associated with the vertices
  !           * with the edges
  !           * with the faces
  !           * in the volume
  !  36 : 16-node  3rd order quadrangle
  !            4 nodes associated with the vertices
  !            8 with the edges
  !            4 in the interior
  !  37 : 25-node  4th order quadrangle
  !            4 nodes associated with the vertices
  !           12 with the edges
  !            9 in the interior
  !  38 : 36-node  5th order quadrangle
  !            4 nodes associated with the vertices
  !           16 with the edges
  !           16 in the interior
  !  39 : 12-node  3rd order EDGE-DEFINED quadrangle
  !            4 nodes associated with the vertices
  !            8 with the edges
  !  40 : 16-node  4th order EDGE-DEFINED quadrangle
  !            4 nodes associated with the vertices
  !           12 with the edges
  !  41 : 20-node  5th order EDGE-DEFINED quadrangle
  !            4 nodes associated with the vertices
  !           16 with the edges
  !  42 : 28-node  6th order triangle
  !            3 nodes associated with the vertices
  !           15 with the edges
  !           10 in the interior
  !  43 : 36-node  7th order triangle
  !            3 nodes associated with the vertices
  !           18 with the edges
  !           15 in the interior
  !  44 : 45-node  8th order triangle
  !            3 nodes associated with the vertices
  !           21 with the edges
  !           21 in the interior
  !  45 : 55-node  9th order triangle
  !            3 nodes associated with the vertices
  !           24 with the edges
  !           28 in the interior
  !  46 : 66-node 10th order triangle
  !            3 nodes associated with the vertices
  !           27 with the edges
  !           36 in the interior
  !  47 : 49-node  6th order quadrangle
  !            4 nodes associated with the vertices
  !           20 with the edges
  !           25 in the interior
  !  48 : 64-node  7th order quadrangle
  !            4 nodes associated with the vertices
  !           24 with the edges
  !           36 in the interior
  !  49 : 81-node  8th order quadrangle
  !            4 nodes associated with the vertices
  !           28 with the edges
  !           49 in the interior
  !  50 : 100-node  9th order quadrangle
  !            4 nodes associated with the vertices
  !           32 with the edges
  !           64 in the interior
  !  51 : 121-node 10th order quadrangle
  !            4 nodes associated with the vertices
  !           36 with the edges
  !           81 in the interior
  !  52 : 18-node  6th order EDGE-DEFINED triangle
  !            3 nodes associated with the vertices
  !           15 with the edges
  !  53 : 21-node  7th order EDGE-DEFINED triangle
  !            3 nodes associated with the vertices
  !           18 with the edges
  !  54 : 24-node  8th order EDGE-DEFINED triangle
  !            3 nodes associated with the vertices
  !           21 with the edges
  !  55 : 27-node  9th order EDGE-DEFINED triangle
  !            3 nodes associated with the vertices
  !           24 with the edges
  !  56 : 30-node 10th order EDGE-DEFINED triangle
  !            3 nodes associated with the vertices
  !           27 with the edges
  !  57 : 24-node  6th order EDGE-DEFINED quadrangle
  !            4 nodes associated with the vertices
  !           20 with the edges
  !  58 : 28-node  7th order EDGE-DEFINED quadrangle
  !            4 nodes associated with the vertices
  !           24 with the edges
  !  59 : 32-node  8th order EDGE-DEFINED quadrangle
  !            4 nodes associated with the vertices
  !           28 with the edges
  !  60 : 36-node  9th order EDGE-DEFINED quadrangle
  !            4 nodes associated with the vertices
  !           32 with the edges
  !  61 : 40-node 10th order EDGE-DEFINED quadrangle
  !            4 nodes associated with the vertices
  !           36 with the edges
  !  62 : 7-node  6th order line
  !            2 nodes associated with the vertices
  !            5 on the edge
  !  63 : 8-node  7th order line
  !            2 nodes associated with the vertices
  !            6 on the edge
  !  64 : 9-node  8th order line
  !            2 nodes associated with the vertices
  !            7 on the edge
  !  65 : 10-node  9th order line
  !            2 nodes associated with the vertices
  !            8 on the edge
  !  66 : 11-node 10th order line
  !            2 nodes associated with the vertices
  !            9 on the edge
  !  67 : 2-node  1st order line border
  !            * nodes associated with the vertices
  !            * on the edge
  !  68 : 3-node  1st order triangle border
  !            * nodes associated with the vertices
  !            * with the edges
  !            * with the interior
  !  69 : ?-node 1st order polygon border
  !            * nodes associated with the vertices
  !            * with the edges
  !            * with the faces
  !            * in the volume
  !  70 : 2-node  1st order line child
  !            * nodes associated with the vertices
  !            * on the edge
  !  71 : 84-node  6th order tetrahedron
  !            4 nodes associated with the vertices
  !           30 with the edges
  !           40 with the faces
  !           10 in the volume
  !  72 : 120-node  7th order tetrahedron
  !            4 nodes associated with the vertices
  !           36 with the edges
  !           60 with the faces
  !           20 in the volume
  !  73 : 165-node  8th order tetrahedron
  !            4 nodes associated with the vertices
  !           42 with the edges
  !           84 with the faces
  !           35 in the volume
  !  74 : 220-node  9th order tetrahedron
  !            4 nodes associated with the vertices
  !           48 with the edges
  !          112 with the faces
  !           56 in the volume
  !  75 : 286-node 10th order tetrahedron
  !            4 nodes associated with the vertices
  !           54 with the edges
  !          144 with the faces
  !           84 in the volume
  !  79 : 34-node  6th order EDGE-DEFINED tetrahedron
  !            4 nodes associated with the vertices
  !           30 with the edges
  !  80 : 40-node  7th order EDGE-DEFINED tetrahedron
  !            4 nodes associated with the vertices
  !           36 with the edges
  !  81 : 46-node  8th order EDGE-DEFINED tetrahedron
  !            4 nodes associated with the vertices
  !           42 with the edges
  !  82 : 52-node  9th order EDGE-DEFINED tetrahedron
  !            4 nodes associated with the vertices
  !           48 with the edges
  !  83 : 58-node 10th order EDGE-DEFINED tetrahedron
  !            4 nodes associated with the vertices
  !           54 with the edges
  !  84 : 1-node  0th order line
  !            1 on the edge
  !  85 : 1-node  0th order triangle
  !            1 in the interior
  !  86 : 1-node  0th order quadrangle
  !            1 in the interior
  !  87 : 1-node  0th order tetrahedron
  !            1 in the volume
  !  88 : 1-node  0th order hexahedron
  !            1 in the volume
  !  89 : 1-node  0th order prism
  !            1 in the volume
  !  90 : 40-node  3rd order prism
  !            6 nodes associated with the vertices
  !           18 with the edges
  !            2 with the triangular faces
  !           12 with the quadrangular faces
  !            2 in the volume
  !  91 : 75-node  4th order prism
  !            6 nodes associated with the vertices
  !           27 with the edges
  !            6 with the triangular faces
  !           27 with the quadrangular faces
  !            9 in the volume
  !  92 : 64-node  3rd order hexahedron
  !            8 nodes associated with the vertices
  !           24 with the edges
  !           24 with the faces
  !            8 in the volume
  !  93 : 125-node  4th order hexahedron
  !            8 nodes associated with the vertices
  !           36 with the edges
  !           54 with the faces
  !           27 in the volume
  !  94 : 216-node  5th order hexahedron
  !            8 nodes associated with the vertices
  !           48 with the edges
  !           96 with the faces
  !           64 in the volume
  !  95 : 343-node  6th order hexahedron
  !            8 nodes associated with the vertices
  !           60 with the edges
  !          150 with the faces
  !          125 in the volume
  !  96 : 512-node  7th order hexahedron
  !            8 nodes associated with the vertices
  !           72 with the edges
  !          216 with the faces
  !          216 in the volume
  !  97 : 729-node  8th order hexahedron
  !            8 nodes associated with the vertices
  !           84 with the edges
  !          294 with the faces
  !          343 in the volume
  !  98 : 1000-node  9th order hexahedron
  !            8 nodes associated with the vertices
  !           96 with the edges
  !          384 with the faces
  !          512 in the volume
  !  99 : 32-node  3rd order EDGE-DEFINED hexahedron
  !            8 nodes associated with the vertices
  !           24 with the edges
  ! 100 : 44-node  4th order EDGE-DEFINED hexahedron
  !            8 nodes associated with the vertices
  !           36 with the edges
  ! 101 : 56-node  5th order EDGE-DEFINED hexahedron
  !            8 nodes associated with the vertices
  !           48 with the edges
  ! 102 : 68-node  6th order EDGE-DEFINED hexahedron
  !            8 nodes associated with the vertices
  !           60 with the edges
  ! 103 : 80-node  7th order EDGE-DEFINED hexahedron
  !            8 nodes associated with the vertices
  !           72 with the edges
  ! 104 : 92-node  8th order EDGE-DEFINED hexahedron
  !            8 nodes associated with the vertices
  !           84 with the edges
  ! 105 : 104-node  9th order EDGE-DEFINED hexahedron
  !            8 nodes associated with the vertices
  !           96 with the edges
  ! 106 : 126-node  5th order prism
  !            6 nodes associated with the vertices
  !           36 with the edges
  !           12 with the triangular faces
  !           48 with the quadrangular faces
  !           24 in the volume
  ! 107 : 196-node  6th order prism
  !            6 nodes associated with the vertices
  !           45 with the edges
  !           20 with the triangular faces
  !           75 with the quadrangular faces
  !           50 in the volume
  ! 108 : 288-node  7th order prism
  !            6 nodes associated with the vertices
  !           54 with the edges
  !           30 with the triangular faces
  !          108 with the quadrangular faces
  !           90 in the volume
  ! 109 : 405-node  8th order prism
  !            6 nodes associated with the vertices
  !           63 with the edges
  !           42 with the triangular faces
  !          147 with the quadrangular faces
  !          147 in the volume
  ! 110 : 550-node  9th order prism
  !            6 nodes associated with the vertices
  !           72 with the edges
  !           56 with the triangular faces
  !          192 with the quadrangular faces
  !          224 in the volume
  ! 111 : 24-node  3rd order EDGE-DEFINED prism
  !            6 nodes associated with the vertices
  !           18 with the edges
  ! 112 : 33-node  4th order EDGE-DEFINED prism
  !            6 nodes associated with the vertices
  !           27 with the edges
  ! 113 : 42-node  5th order EDGE-DEFINED prism
  !            6 nodes associated with the vertices
  !           36 with the edges
  ! 114 : 51-node  6th order EDGE-DEFINED prism
  !            6 nodes associated with the vertices
  !           45 with the edges
  ! 115 : 60-node  7th order EDGE-DEFINED prism
  !            6 nodes associated with the vertices
  !           54 with the edges
  ! 116 : 69-node  8th order EDGE-DEFINED prism
  !            6 nodes associated with the vertices
  !           63 with the edges
  ! 117 : 78-node  9th order EDGE-DEFINED prism
  !            6 nodes associated with the vertices
  !           72 with the edges
  ! 118 : 30-node  3rd order pyramid
  !            5 nodes associated with the vertices
  !           16 with the edges
  !            4 with the triangular faces
  !            4 with the quadrangular faces
  !            1 in the volume
  ! 119 : 55-node  4th order pyramid
  !            5 nodes associated with the vertices
  !           24 with the edges
  !           12 with the triangular faces
  !            9 with the quadrangular faces
  !            5 in the volume
  ! 120 : 91-node  5th order pyramid
  !            5 nodes associated with the vertices
  !           32 with the edges
  !           24 with the triangular faces
  !           16 with the quadrangular faces
  !           14 in the volume
  ! 121 : 140-node  6th order pyramid
  !            5 nodes associated with the vertices
  !           40 with the edges
  !           40 with the triangular faces
  !           25 with the quadrangular faces
  !           30 in the volume
  ! 122 : 204-node  7th order pyramid
  !            5 nodes associated with the vertices
  !           48 with the edges
  !           60 with the triangular faces
  !           36 with the quadrangular faces
  !           55 in the volume
  ! 123 : 285-node  8th order pyramid
  !            5 nodes associated with the vertices
  !           56 with the edges
  !           84 with the triangular faces
  !           49 with the quadrangular faces
  !           91 in the volume
  ! 124 : 385-node  9th order pyramid
  !            5 nodes associated with the vertices
  !           64 with the edges
  !          112 with the triangular faces
  !           64 with the quadrangular faces
  !          140 in the volume
  ! 125 : 21-node  3rd order EDGE-DEFINED pyramid
  !            5 nodes associated with the vertices
  !           16 with the edges
  ! 126 : 29-node  4th order EDGE-DEFINED pyramid
  !            5 nodes associated with the vertices
  !           24 with the edges
  ! 127 : 37-node  5th order EDGE-DEFINED pyramid
  !            5 nodes associated with the vertices
  !           32 with the edges
  ! 128 : 45-node  6th order EDGE-DEFINED pyramid
  !            5 nodes associated with the vertices
  !           40 with the edges
  ! 129 : 53-node  7th order EDGE-DEFINED pyramid
  !            5 nodes associated with the vertices
  !           48 with the edges
  ! 130 : 61-node  8th order EDGE-DEFINED pyramid
  !            5 nodes associated with the vertices
  !           56 with the edges
  ! 131 : 69-node  9th order EDGE-DEFINED pyramid
  !            5 nodes associated with the vertices
  !           64 with the edges
  ! 132 : 1-node  0th order pyramid
  !            1 in the volume
  ! 133 : 1-node  1st order point Xfem
  !            * nodes associated with the vertices
  !            * with the edges
  !            * with the faces
  !            * in the volume
  ! 134 : 2-node  1st order line Xfem
  !            * nodes associated with the vertices
  !            * on the edge
  ! 135 : 3-node  1st order triangle Xfem
  !            * nodes associated with the vertices
  !            * with the edges
  !            * with the interior
  ! 136 : 4-node  1st order tetrahedron Xfem
  !            * nodes associated with the vertices
  !            * with the edges
  !            * with the faces
  !            * in the volume
  ! 137 : 16-node  3rd order EDGE-DEFINED tetrahedron
  !            4 nodes associated with the vertices
  !           12 with the edges
  ! 138 : MINI-node  1st order triangle
  !            * nodes associated with the vertices
  !            * with the edges
  !            * with the interior
  ! 139 : MINI-node  1st order tetrahedron
  !            * nodes associated with the vertices
  !            * with the edges
  !            * with the faces
  !            * in the volume
  ! 140 : 4-node  1st order trihedron
  !            * nodes associated with the vertices
  !            * with the edges
  !            * with the faces
  !            * in the volume
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###                                                                         ###
!###  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  ###
!###  |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  ###
!###                                                                         ###
!###               NEWLY ADDED ELEMENTS FROM GMSH SOURCE CODE                ###
!###       THAT DO NOT SEEM TO HAVE ANY DOCUMENTATION RELATING TO THEM       ###
!###                                                                         ###
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
  !
  ! ###################################
  ! ###  Module Derived Data Types  ###
  ! ###################################
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
  logical(lk), save :: meshformat_successful
  logical(lk), save :: nodes_successful
  logical(lk), save :: elements_successful
  logical(lk), save :: physicalnames_successful
  !
  integer, save :: iogmsh
  !
  ! ###################################
  ! ###  Module Allocatable Arrays  ###
  ! ###################################
  !
 !type(elem_t), allocatable, dimension(:) :: cell_data
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
  logical(lk) :: end_of_file
  character(len=80) :: section_header
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_gmsh_gridfile"
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
  open (newunit=iogmsh,file=gridfile,status="old",action="read", &
                       position="rewind",iostat=ierr,iomsg=error_message)
  call io_error(pname,gridfile,1,__LINE__,__FILE__,ierr,error_message)
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
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
    !
  else if ( .not. nodes_successful ) then
    !
    if (mypnum == 0) write (iout,3)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
    !
  else if ( .not. elements_successful ) then
    !
    if (mypnum == 0) write (iout,4)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
    !
  else if ( .not. physicalnames_successful ) then
    !
    if (mypnum == 0) write (iout,5)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
    !
  else
    !
    ! We are done reading the grid file so close it
    !
    close (iogmsh,iostat=ierr,iomsg=error_message)
    call io_error(pname,gridfile,2,__LINE__,__FILE__,ierr,error_message)
    !
  end if
  !
  ! Finish creating the arrays defining the grid geometry.
  ! bface, nodes_of_cell_ptr, nodes_of_cell, etc.
  !
  call create_solver_geom_arrays_gmsh
  !
  ! Deallocate the allocatable module arrays cell_data and phys_grps
  ! before we leave this subroutine.
  !
 !if (allocated(cell_data)) then
 !  deallocate ( cell_data , stat=ierr , errmsg=error_message )
 !  call alloc_error(pname,"cell_data",2,__LINE__,__FILE__,ierr,error_message)
 !end if
  if (allocated(phys_grps)) then
    deallocate ( phys_grps , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"phys_grps",2,__LINE__,__FILE__,ierr,error_message)
  end if
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
  logical(lk),          intent(out) :: end_of_file
  character(len=*), intent(out) :: text_out
  !
  !.. Local Scalars ..
  integer :: ierr,n
  character(len=len(text_out)) :: text
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_gmsh_section"
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
    read (iogmsh,1,iostat=ierr,end=100,err=101) text
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
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
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
  logical(lk) :: end_of_file
  real(wp) :: gmsh_version_num
  character(len=80) :: text
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_gmsh_section_meshformat"
  !
continue
  !
  ! Read the MeshFormat section of the GMSH grid file
  !
  read (iogmsh,*,iostat=ierr,end=100,err=101) gmsh_version_num, &
                                              gmsh_file_type, &
                                              gmsh_data_size
  !
  ! Check that the GMSH version number is supported
  !
  gmsh_version_num = gmsh_version_num * ten
  !
  if (nint(gmsh_version_num) /= 22) then
    if (mypnum == 0) write (iout,1)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
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
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  !
  ! Read error continuation statements
  !
  100 if (mypnum == 0) write (iout,2) ierr
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  !
  101 if (mypnum == 0) write (iout,4) ierr
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
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
  use geovar, only : nr
  use ovar,   only : grid_scaling_factor
  !
  !.. Local Scalars ..
  integer :: i,j,n,npts,ierr
  logical(lk) :: end_of_file
  real(wp) :: minz,maxz
  character(len=80) :: text
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable, dimension(:,:) :: temp_nodes
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_gmsh_section_nodes"
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
  read (iogmsh,*,iostat=ierr,end=100,err=101) npts
  !
  ! Allocate the temporary array to read in the coordinates of the grid nodes
  !
  allocate ( temp_nodes(1:nr,1:npts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_nodes",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read the coordinates of the grid nodes
  !
  do n = 1,npts
    read (iogmsh,*,iostat=ierr,end=100,err=101) i, (temp_nodes(j,i), j=1,nr)
  end do
  !
  ! If the minimum and maximum values of the coordinate in the third
  ! dimension are equal, reduce nr to 2 since this is a 2D grid.
  !
  minz = minval(temp_nodes(nr,:))
  maxz = maxval(temp_nodes(nr,:))
  !
  if (abs(maxz-minz) <= eps6) nr = 2
  !
  ! Allocate grid if it has not yet been allocated
  !
  if (.not. allocated(grid)) then
    allocate ( grid , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Allocate the xyz component of the grid derived type
  ! and use the temp_nodes array as the data source
  !
  allocate ( grid%xyz(1:nr,1:npts) , source=temp_nodes(1:nr,1:npts) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid%xyz",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Apply the grid scaling factor to the coordinates of the grid nodes
  !
  do n = 1,size(grid%xyz,dim=2)
    do i = 1,size(grid%xyz,dim=1)
      grid%xyz(i,n) = grid_scaling_factor*grid%xyz(i,n)
    end do
  end do
  !
  ! Deallocate temp_nodes now that we no longer need it
  !
  deallocate ( temp_nodes , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_nodes",2,__LINE__,__FILE__,ierr,error_message)
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
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  !
  ! Read error continuation statements
  !
  100 if (mypnum == 0) write (iout,1) ierr
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  !
  101 if (mypnum == 0) write (iout,3) ierr
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
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
  integer :: i,j,m,n,ierr,nods
  integer :: ntype,ntags,npts,geom
  integer :: num_cell_data
  logical(lk) :: end_of_file
  character(len=80) :: text,array_name
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_gmsh_section_elements"
  !
continue
  !
  ! Read the Elements section of the GMSH grid file
  !
  read (iogmsh,*,iostat=ierr,end=100,err=101) num_cell_data
  !
  ! Allocate grid if it has not yet been allocated
  !
  if (.not. allocated(grid)) then
    allocate ( grid , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Allocate the elem component of the grid derived type
  ! to store the element definitions
  !
  allocate ( grid%elem(1:num_cell_data) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid%elem",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read the element definitions
  !
  do n = 1,size(grid%elem)
    !
    ! Read the first three values on the line for the current cell
    ! and do not advance to the next record. We need to allocate
    ! some memory based on the second two values and then we will
    ! read the remaining information on this line.
    !
    read (iogmsh,*,iostat=ierr,end=100,err=101) i, ntype, ntags
    backspace (iogmsh)
    !
    ! Get the properties for this type of element
    !
    grid%elem(i)%prop = element_properties(ntype)
    !
    ! Check for bad element types
    !
    if ( (grid%elem(i)%prop%npts == -1) .or. &
         (grid%elem(i)%prop%geom == Geom_Unknown) ) then
      !
      write (error_message,5) ntype
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
      !
   !else if (grid%elem(i)%prop%family /= Complete) then
   !  !
   !  write (error_message,6) ntype
   !  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
   !  !
    else if (grid%elem(i)%prop%warning) then
      !
      if (mypnum == 0) then
        write (iout,7) grid%elem(i)%prop%npts, &
                       grid%elem(i)%prop%order, &
                       trim(adjustl(Geom_String(grid%elem(i)%prop%geom)))
      end if
      !
    end if
    !
    npts = grid%elem(i)%prop%npts
    nods = geom_nodes( grid%elem(i)%prop%geom )
    !
    ! Allocate the pts array within grid%elem for the current cell
    !
    write (array_name,4) "grid%elem(", i, ")%pts"
    allocate ( grid%elem(i)%pts(1:npts) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Allocate the tags array within grid%elem for the current cell
    !
    write (array_name,4) "grid%elem(", i, ")%tags "
    allocate ( grid%elem(i)%tags(1:ntags) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Allocate the nodes array within grid%elem for the current cell
    !
    write (array_name,4) "grid%elem(", i, ")%nodes "
    allocate ( grid%elem(i)%nodes(1:nods) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Continue reading the remaining information for this cell
    !
    read (iogmsh,*,iostat=ierr,end=100,err=101) m, m, m, &
                              (grid%elem(i)%tags(j),  j=1,ntags), &
                              (grid%elem(i)%pts(j), j=1,npts)
    !
    ! The recursive point ordering format used by Gmsh means the nodes
    ! (i.e., corners) of each element are provided by grid%elem(i)%pts(1:nods).
    ! Copy the indices for the element nodes into the grid%elem(i)%nodes
    ! array since the pts array will get reordered next.
    !
    grid%elem(i)%nodes(1:nods) = grid%elem(i)%pts(1:nods)
    !
    ! Reorder the points for this grid element from the recursive format
    ! used by Gmsh to the i,j,k ordering used by the solver.
    !
    grid%elem(i) = reorder_gmsh_pts( grid%elem(i) )
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
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  !
  ! Read error continuation statements
  !
  100 if (mypnum == 0) write (iout,1) ierr
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  !
  101 if (mypnum == 0) write (iout,3) ierr
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
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
  5 format (" ERROR!!! ENCOUNTERED INVALID ELEMENT TYPE (",i0, &
            ") IN THE GMSH GRID FILE!")
  6 format (" ERROR!!! ENCOUNTERED INVALID ELEMENT TYPE (",i0, &
            ") IN THE GMSH GRID FILE! THIS IS AN INCOMPLETE", &
            " (EDGE-DEFINED OR SERENDIPITY-TYPE) ELEMENT", &
            " WHICH ARE NOT SUPPORTED AT THIS TIME!")
  7 format (" WARNING!!! ",i0,"-point, ",i0,"-order ",a, &
            " geometry could cause problems!")
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
  integer :: num_phys_grps
  logical(lk) :: end_of_file
  character(len=80) :: text,phys_nam
  !
  !.. Local Parameters ..
  character(len=*), parameter :: delimiters(1:4) = (/"'",'"',",",achar(9)/)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_gmsh_section_physicalnames"
  !
continue
  !
  ! Read the PhysicalNames section of the GMSH grid file
  !
  read (iogmsh,*,iostat=ierr,end=100,err=101) num_phys_grps
  !
  ! Allocate the array for the element definitions
  !
  allocate ( phys_grps(1:num_phys_grps) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"phys_grps",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read the element definitions
  !
  do n = 1,size(phys_grps)
    !
    read (iogmsh,*,iostat=ierr,end=100,err=101) phys_dim, phys_id, phys_nam
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
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  !
  ! Read error continuation statements
  !
  100 if (mypnum == 0) write (iout,1) ierr
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  !
  101 if (mypnum == 0) write (iout,3) ierr
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
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
subroutine create_solver_geom_arrays_gmsh
  !
  ! Procedure to create the remaining grid geometry arrays using
  ! the information stored in the arrays grid%elem and phys_grps.
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use geovar, only : nr,ncell,nnode,nfbnd,nbfai,bface
  use geovar, only : nodes_of_cell_ptr,nodes_of_cell
  use geovar, only : cell_geom,cell_order,xyz_nodes
  use ovar,   only : loc_solution_pts,loc_flux_pts
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,ifa,nf,h1,h2,i,j,ie
  integer :: np,nc,nfp,pid,pdm,inum,ierr
  integer :: igrps,bgrps,b0grps,b0cells
  integer :: host_cell,host_geom,host_side
  integer :: npts,ncnods,nfnods,bnd_cell
  logical(lk) :: error_found
  character(len=80) :: array_name
  !
  !.. Local Arrays ..
  integer, dimension(1:8) :: int_nodes
  integer, dimension(1:4) :: face_nodes
  integer, dimension(1:8) :: marker
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:) :: int_grps
  integer, allocatable, dimension(:) :: bnd_grps
  integer, allocatable, dimension(:) :: b0_grps
  integer, allocatable, dimension(:) :: mapp
  integer, allocatable, dimension(:) :: mapr
  integer, allocatable, dimension(:) :: lpoin
  integer, allocatable, dimension(:) :: cells_with_node
  integer, allocatable, dimension(:) :: cells_with_node_ptr
  !
  logical(lk), allocatable, dimension(:) :: mskp
  !
  type(elem_t), allocatable, dimension(:) :: tmp_cell_data
  type(phys_t), allocatable, dimension(:) :: tmp_phys_grps
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_solver_geom_arrays_gmsh"
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
  ! Count the number of interior physical groups
  !
  igrps = count(phys_grps(:)%phys_dim == nr)
  !
  ! Count the number of boundary physical groups
  ! NOTE: Only count the 0-dimensional physical groups if the interior is 1D
  !
  if (nr > 1) then
    bgrps = count( phys_grps(:)%phys_dim < nr .and. &
                   phys_grps(:)%phys_dim > 0 )
    b0grps = count(phys_grps(:)%phys_dim == 0)
  else
    bgrps = count(phys_grps(:)%phys_dim < nr)
    b0grps = 0
  end if
  !
  ! Allocate the arrays to identify the interior and boundary physical groups
  !
  allocate ( int_grps(1:igrps) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"int_grps",1,__LINE__,__FILE__,ierr,error_message)
  allocate ( bnd_grps(1:bgrps) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bnd_grps",1,__LINE__,__FILE__,ierr,error_message)
  allocate ( b0_grps(1:b0grps) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"b0_grps",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Collect the interior physical groups
  !
  int_grps(1:igrps) = pack( intseq(1,size(phys_grps)) , &
                            phys_grps(:)%phys_dim == nr )
  !
  ! Collect the boundary physical groups
  ! NOTE: Only collect the 0-dimensional physical groups if the interior is 1D
  !
  if (nr > 1) then
    bnd_grps(1:bgrps) = pack( intseq(1,size(phys_grps)) , &
                              phys_grps(:)%phys_dim < nr .and. &
                              phys_grps(:)%phys_dim > 0 )
    b0_grps(1:b0grps) = pack( intseq(1,size(phys_grps)) , &
                              phys_grps(:)%phys_dim == 0 )
  else
    bnd_grps(1:bgrps) = pack( intseq(1,size(phys_grps)) , &
                              phys_grps(:)%phys_dim < nr )
  end if
  !
  ! Check that the total number of interior and boundary groups is correct
  !
  if ( igrps+bgrps+b0grps /= size(phys_grps) ) then
    if (mypnum == 0) write (iout,1) igrps,bgrps,b0grps,size(phys_grps)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Reorder grid%elem so that the interior cells are first followed by
  ! the boundary faces
  !
  ! First, count the number of 0-D elements in grid%elem
  !
  b0cells = 0
  if (b0grps > 0) then
    do j = 1,size(grid%elem)
      if (allocated(grid%elem(j)%tags)) then
        if (any(grid%elem(j)%tags(1) == phys_grps(b0_grps(:))%phys_id)) then
          b0cells = b0cells + 1
        end if
      end if
    end do
  end if
  !
  ! Allocate the temporary array to store the reordered grid%elem
  !
  n = size(grid%elem) - b0cells
  allocate ( tmp_cell_data(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"tmp_cell_data",1,__LINE__,__FILE__,ierr,error_message)
  !
  n = 0
  !
  ! Find all the interior cells and put them at the front of tmp_cell_data
  !
  do j = 1,size(grid%elem)
    if (allocated(grid%elem(j)%tags)) then
      if (any(grid%elem(j)%tags(1) == phys_grps(int_grps(:))%phys_id)) then
        n = n + 1
        tmp_cell_data(n) = grid%elem(j)
      end if
    end if
  end do
  !
  ! Now append all the boundary faces to the end of tmp_cell_data
  !
  do j = 1,size(grid%elem)
    if (allocated(grid%elem(j)%tags)) then
      if (any(grid%elem(j)%tags(1) == phys_grps(bnd_grps(:))%phys_id)) then
        n = n + 1
        tmp_cell_data(n) = grid%elem(j)
      end if
    end if
  end do
  !
  if (n /= size(grid%elem) - b0cells) then
    if (mypnum == glb_root) then
      write (iout,4)
      write (iout,9) size(grid%elem),b0cells,size(grid%elem)-b0cells,n
      write (iout,4)
    end if
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Finally, use move_alloc to rebuild grid%elem with the new ordering
  ! from tmp_cell_data, and then deallocate tmp_cell_data
  !
  call move_alloc( from=tmp_cell_data , to=grid%elem )
  !
  ! Deallocate the b0_grps array if it is allocated
  !
  if (allocated(b0_grps)) then
    deallocate ( b0_grps , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"b0_grps",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Remove any 0-D physical groups if the interior is not 1D
  !
  if (b0grps > 0) then
    !
    ! Allocate the temporary physical derived type
    !
    n = size(phys_grps) - b0grps
    allocate ( tmp_phys_grps(1:n) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"tmp_phys_grps",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    ! Save all the non-0-D physical groups to
    ! the temporary physical group derived type
    !
    n = 0
    do i = 1,size(phys_grps)
      if (phys_grps(i)%phys_dim > 0) then
        n = n + 1
        tmp_phys_grps(n) = phys_grps(i)
      end if
    end do
    !
    if (n /= size(phys_grps) - b0grps) then
      write (error_message,10) size(phys_grps),b0grps,size(phys_grps)-b0grps,n
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
    !
    ! Use move_alloc to rebuild phys_grps without any 0-D physical groups
    !
    call move_alloc( from=tmp_phys_grps , to=phys_grps )
    !
    ! Finally, rebuild the int_grps and bnd_grps arrays since the
    ! index of these groups within phys_grps may have now changed.
    ! NOTE: The number of interior and boundary groups should not have changed
    !       so we do not need to reallocate these arrays to a different size
    !
    int_grps(1:igrps) = pack( intseq(1,size(phys_grps)) , &
                              phys_grps(:)%phys_dim == nr )
    bnd_grps(1:bgrps) = pack( intseq(1,size(phys_grps)) , &
                              phys_grps(:)%phys_dim < nr )
    !
  end if
  !
  ! #######
  ! ## 4 ##
  ! #######
  !
  ! Loop through the physical groups and collect
  ! the elements for each physical group
  !
  do i = 1,size(phys_grps)
    !
    pid = phys_grps(i)%phys_id
    pdm = phys_grps(i)%phys_dim
    !
    ! Count the number of elements that match this physical groups id number
    !
    inum = 0
    do j = 1,size(grid%elem)
      if (allocated(grid%elem(j)%tags)) then
        inum = inum + merge( 1 , 0 , grid%elem(j)%tags(1) == pid )
      end if
    end do
   !inum = count( grid%elem(:)%tags(1) == pid )
    !
    ! Allocate the cells array for this physical group
    !
    write (array_name,3) "phys_grps(", i, ")%cells"
    allocate ( phys_grps(i)%cells(1:inum) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Create the list of elements that are in this physical group
    !
    inum = 0
    do j = 1,size(grid%elem)
      if (allocated(grid%elem(j)%tags)) then
        if (grid%elem(j)%tags(1) == pid) then
          inum = inum + 1
          phys_grps(i)%cells(inum) = j
        end if
      end if
    end do
   !phys_grps(i)%cells(1:inum) = pack( intseq(1,size(grid%elem)) , &
   !                                   grid%elem(:)%tags(1) == pid )
    !
  end do
  !
  ! #######
  ! ## 5 ##
  ! #######
  !
  ! Count the total number of interior cells and boundary faces
  !
  ncell = 0
  do i = 1,size(int_grps)
    ncell = ncell + size( phys_grps(int_grps(i))%cells )
  end do
  !
  nfbnd = 0
  do i = 1,size(bnd_grps)
    nfbnd = nfbnd + size( phys_grps(bnd_grps(i))%cells )
  end do
  !
  ! Check that the number of interior cells and boundary faces are correct
  !
  if ( ncell+nfbnd /= size(grid%elem) ) then
    if (mypnum == 0) write (iout,2) ncell,nfbnd,size(grid%elem)
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  !
  ! #######
  ! ## 6 ##
  ! #######
  !
  ! Go through all the grid%elem(:)%nodes arrays and redo the node indices
  ! as if this was a linear grid with only the corner nodes defined. This
  ! will allow us to create the nodes_of_cell, nodes_of_cell_ptr, and xyz_nodes
  ! arrays later on in the way they will be needed by the rest of the code.
  !
  ! Find the corner nodes of all the grid elements
  !
  npts = size(grid%xyz,dim=2)
  !
  allocate ( mskp(1:npts) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"mskp",1,__LINE__,__FILE__,ierr,error_message)
  !
  do n = 1,size(grid%elem)
    mskp(grid%elem(n)%nodes(:)) = true
  end do
  !
  ! Count the number of corner nodes
  !
  nnode = count( mskp )
  !
  if (mypnum == 0) then
    write (iout,*)
    write (iout,3) "ncell = ",ncell
    write (iout,3) "nfbnd = ",nfbnd
    write (iout,3) "nnode = ",nnode
    write (iout,3) "npts  = ",npts
    write (iout,*)
  end if
  !
  ! Create the forward map (mapp) and the reverse map (mapr)
  !
  allocate ( mapp(1:nnode) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"mapp",1,__LINE__,__FILE__,ierr,error_message)
  allocate ( mapr(1:npts) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"mapr",1,__LINE__,__FILE__,ierr,error_message)
  !
  mapp(1:nnode) = pack( intseq(1,npts) , mask=mskp )
  mapr(1:npts) = unpack( intseq(1,nnode) , mask=mskp , field=0 )
  !
  ! Create the xyz_nodes array
  !
  allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_nodes",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Save the coordinates of the corner nodes into the xyz_array
  !
  do n = 1,nnode
    xyz_nodes(1:nr,n) = grid%xyz(1:nr,mapp(n))
  end do
 !do j = 1,17
 !  do i = 1,49
 !    n = (i-1)*17 + j
 !    write (182,'(3(1x,es14.6))') (xyz_nodes(n1,n),n1=1,nr)
 !  end do
 !end do
  !
  ! Convert the node indices stored in grid%elem(:)%nodes
  ! from the range 1:npts to the range 1:nnode
  !
  do n = 1,size(grid%elem)
    grid%elem(n)%nodes(:) = mapr( grid%elem(n)%nodes(:) )
  end do
  !
  deallocate ( mskp , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"mskp",2,__LINE__,__FILE__,ierr,error_message)
  deallocate ( mapp , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"mapp",2,__LINE__,__FILE__,ierr,error_message)
  deallocate ( mapr , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"mapr",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! #######
  ! ## 7 ##
  ! #######
  !
  ! Allocate the nodes_of_cell_ptr array
  !
  allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  allocate ( cell_geom(1:ncell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_geom",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( cell_order(1:ncell) , source=n_order , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_order",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Loop through the cells for each interior physical
  ! group and use grid%elem to create nodes_of_cell_ptr
  !
 !n = 0
 !do i = 1,size(int_grps)
 !  np = int_grps(i)
 !  do j = 1,size(phys_grps(np)%cells)
 !    n = n + 1
 !    nc = phys_grps(np)%cells(j)
 !    nodes_of_cell_ptr(n+1) = size(grid%elem(nc)%nodes)
 !    cell_geom(n) = grid%elem(nc)%prop%geom
 !  end do
 !end do
  !
  ! NOTE: THERE MIGHT BE A BUG HERE BECAUSE IT SEEMS TO ASSUME THAT THE
  !       ORDERING WITHIN NODES_OF_CELL_POINTER IS RELATED TO THE ORDERING
  !       IN GRID%ELEM. THE COMMENTED CODE ABOVE REMOVES THIS ASSUMPTION.
  !
  ! NOTE: THIS ONLY COUNTS THE NODES DEFINING THE LINEAR CELL GEOMETRY
  !       AND DOESNT INCLUDE ANY OF THE ADDITIONAL POINTS PROVIDED IF
  !       THIS CELL IS CURVILINEAR (>= 2ND ORDER).
  !
  do i = 1,size(int_grps)
    np = int_grps(i)
    do j = 1,size(phys_grps(np)%cells)
      n = phys_grps(np)%cells(j)
      cell_geom(n) = grid%elem(n)%prop%geom
      nodes_of_cell_ptr(n+1) = size(grid%elem(n)%nodes)
    end do
  end do
  !
  ! Reshuffle nodes_of_cell_ptr
  !
  do n = 2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(n) = nodes_of_cell_ptr(n) + nodes_of_cell_ptr(n-1)
  end do
  !
  ! Deallocate int_grps since we dont need it anymore
  !
  deallocate ( int_grps , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"int_grps",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! #######
  ! ## 8 ##
  ! #######
  !
  ! Allocate the nodes_of_cell array
  !
  allocate ( nodes_of_cell(1:last(nodes_of_cell_ptr)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Loop over the interior cells and fill in the node connectivity
  ! for that cell from the information in grid%elem.
  ! NOTE: THIS ONLY COPIES OVER THE NODES DEFINING THE LINEAR CELL GEOMETRY
  !       AND DOESNT INCLUDE ANY OF THE ADDITIONAL POINTS PROVIDED IF THIS
  !       CELL IS CURVILINEAR (>= 2ND ORDER).
  !
  do n = 1,ncell
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    nodes_of_cell(n1:n2) = grid%elem(n)%nodes(:)
  end do
  !
  ! #######
  ! ## 9 ##
  ! #######
  !
  ! Before we create the bface array, we need to
  ! find the host cells for each boundary face
  !
  ! Allocate cells_with_node_ptr
  !
  allocate ( cells_with_node_ptr(1:nnode+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_with_node_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! Loop over the cells, storing 'ahead'
  !
  do n = 1,ncell
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    do i = n1,n2
      np = nodes_of_cell(i) + 1
      cells_with_node_ptr(np) = cells_with_node_ptr(np) + 1
    end do
  end do
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
  do n = 1,ncell
    !
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    !
    do i = n1,n2
      j = nodes_of_cell(i)
      np = cells_with_node_ptr(j) + 1
      cells_with_node_ptr(j) = np
      cells_with_node(np) = n
    end do
    !
  end do
  !
  ! Finally, reorder cells_with_node_ptr
  !
  do n = size(cells_with_node_ptr),2,-1
    cells_with_node_ptr(n) = cells_with_node_ptr(n-1)
  end do
  cells_with_node_ptr(1) = 0
  !
  ! Allocate the point counting array, lpoin, and the host_cell array
  !
  allocate ( lpoin(1:npts) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"lpoin",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Find the host cells for each boundary face
  !
  do i = 1,size(bnd_grps)
    !
    np = bnd_grps(i)
    !
    face_loop: do j = 1,size(phys_grps(np)%cells)
      !
      n = phys_grps(np)%cells(j)
      nfp = size(grid%elem(n)%nodes)
      !
      face_nodes(1:nfp) = grid%elem(n)%nodes(1:nfp)
      !
      ! Mark the nodes for this boundary face
      !
      lpoin(face_nodes(1:nfp)) = 1
      !
      n1 = cells_with_node_ptr( face_nodes(1) )+1
      n2 = cells_with_node_ptr( face_nodes(1) +1)
      !
      do ie = n1,n2
        nc = cells_with_node(ie)
        if (sum( lpoin(grid%elem(nc)%nodes(:)) ) == nfp) then
          grid%elem(n)%host_cell = nc
          lpoin(face_nodes(1:nfp)) = 0
          cycle face_loop
        end if
      end do
      write (iout,1092) j,n
      !
    end do face_loop
    !
  end do
  !
  ! Deallocate lpoin, cells_with_node, and cells_with_node_ptr
  ! since we no longer need them
  !
  deallocate ( cells_with_node_ptr , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_with_node_ptr",2,__LINE__,__FILE__,ierr, &
                   error_message)
  deallocate ( cells_with_node , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_with_node",2,__LINE__,__FILE__,ierr, &
                   error_message)
  deallocate ( lpoin , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"lpoin",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! ########
  ! ## 10 ##
  ! ########
  !
  ! Create and fill out bface with data from the
  ! PhysicalName section of the GMSH grid file
  !
  ! Make sure nbfai is large enough for high-order boundary faces
  !
 !nbfai = 10
 !do i = 1,size(bnd_grps)
 !  np = bnd_grps(i)
 !  do j = 1,size(phys_grps(np)%cells)
 !    n = phys_grps(np)%cells(j)
 !    nbfai = max( nbfai , size(grid%elem(n)%pts)+10 )
 !  end do
 !end do
  !
  ! Allocate the bface array
  !
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Loop through the faces for each boundary physical
  ! group and use grid%elem to create bface
  !
  nf = 0
  boundary_phys_grps: do i = 1,size(bnd_grps)
    !
    np = bnd_grps(i)
    !
    current_phys_grp: do j = 1,size(phys_grps(np)%cells)
      !
      bnd_cell = phys_grps(np)%cells(j)
      host_cell = grid%elem(bnd_cell)%host_cell
      nf = nf + 1
      !
      bface( 1,nf) = phys_grps(np)%ibc       ! face boundary condition
      bface( 2,nf) = host_cell                      ! host cell of boundary face
      bface( 3,nf) = grid%elem(host_cell)%prop%geom ! host cell geometry
      bface( 5,nf) = i                       ! group for boundary face
      bface(10,nf) = ncell + nf              ! ghost cell of boundary face
      !
      host_cell = bface(2,nf) ! host cell of boundary face
      host_geom = bface(3,nf) ! geometry of host cell
      !
      ! If the cell geometry is a node or edge, change ifa to a special value.
      ! Otherwise, continue on to finding the side of the host cell that is
      ! on the boundary.
      !
      if (all(host_geom /= [Geom_2D,Geom_3D])) then
        !
        if (host_geom == Geom_Node) then
          ifa = -1
        else if (host_geom == Geom_Edge) then
          ifa = -2
        else
          ifa = -10
        end if
        !
      else
        !
        ! Find which face of the host cell is on the boundary
        !
        ! total number of host cell nodes
        ncnods = size(grid%elem(host_cell)%nodes)
        !
        ! Create a local copy of the interior nodes
        !
        int_nodes(1:ncnods) = grid%elem(host_cell)%nodes
        !
        marker = 0 ! Initialize the marker array to zero
        !
        ! Mark each host cell node that is on the boundary
        !
        do n = 1,size(grid%elem(bnd_cell)%nodes)
          where (int_nodes(1:ncnods) == grid%elem(bnd_cell)%nodes(n)) marker = 1
        end do
        !
        ! Loop over the number of faces for this geometry type and find which
        ! face corresponds to the marked nodes
        !
        ifa = -99 ! Initialize ifa to make sure that a side is actually found
        !
        do host_side = 1,geom_faces(host_geom)
          !
          ! Get the number of nodes defining this face
          !
          nfnods = num_face_nodes(host_side,host_geom)
          !
          ! Extract the cell nodes that define the current face
          ! NOTE: The +1 at the end is needed since the values stored in
          !       side_nodes are normally used as offsets, meaning that their
          !       values are in the range [0:ncnods-1]. We want the values
          !       in the range of [1:ncnods] hence the +1.
          !
          face_nodes(1:nfnods) = side_nodes(1:nfnods,host_side,host_geom) + 1
          !
          ! If all of these face nodes have been marked, then this is the
          ! side of the host cell on the boundary. Save the value into ifa
          ! and exit the face loop.
          !
          if ( sum( marker(face_nodes(1:nfnods)) ) == nfnods ) then
            ifa = host_side
            exit
          end if
          !
        end do
        !
      end if
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
            write (iout,5) "A 0-D 'Point'-type",host_cell,bnd_cell
          else if (ifa == -2) then
            write (iout,5) "A 1-D 'Edge'-type",host_cell,bnd_cell
          else if (ifa == -10) then
            write (iout,5) "An unknown type of",host_cell,bnd_cell
          else if (ifa == -99) then
            write (iout,6) host_cell,bnd_cell
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
  deallocate ( bnd_grps , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bnd_grps",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! ########
  ! ## 11 ##
  ! ########
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
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  end if
  !
  ! Make sure that Lobatto nodes are requested if the grid contains triangle
  ! elements since Gauss nodes is not currently set up for triangles.
  !
 !if (any(cell_geom == Geom_Tria)) then
 !  if (mypnum == 0) write (iout,8)
 !  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
 !end if
  !
  if (loc_solution_pts /= Legendre_Gauss_Lobatto .or. &
      loc_flux_pts /= Legendre_Gauss_Lobatto) then
    if (any(cell_geom == Geom_Tria)) then
      if (mypnum == 0) write (iout,7)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
    end if
  end if
  !
  ! Format Statements
  !
  1 format (/,"************************ ERROR! **************************",/, &
              "  The number of interior and boundary physical groups is",/, &
              "  inconsistent with the total number of physical groups!",/, &
              "     number of     interior groups = ",i0,/, &
              "     number of     boundary groups = ",i0,/, &
              "     number of 0-D boundary groups = ",i0,/, &
              "            total number of groups = ",i0,/, &
              "************************ ERROR! **************************",/)
  2 format (/,"************************ ERROR! ************************",/, &
              "  The number of interior cells and boundary faces is",/, &
              "  inconsistent with the total number of grid elements!",/, &
              "          number of interior cells = ",i0,/, &
              "          number of boundary faces = ",i0,/, &
              "     total number of grid elements = ",i0,/, &
              "************************ ERROR! ************************",/)
  3 format (a,i0,:,a)
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
  9 format ("  The final counter for tmp_cell_data does not",/, &
            "  match the size of the grid%elem array!",/, &
            "                             size(grid%elem) = ",i0,/, &
            "                 number of 0-D grid elements = ",i0,/, &
            "         size(grid%elem) - 0-D grid elements = ",i0,/, &
            "                 final tmp_cell_data counter = ",i0,/)
  10 format ("  The final counter for tmp_phys_grps does not",/, &
             "  match the size of the phys_grps array!",/, &
             "                             size(phys_grps) = ",i0,/, &
             "               number of 0-D physical groups = ",i0,/, &
             "       size(phys_grps) - 0-D physical groups = ",i0,/, &
             "                 final tmp_phys_grps counter = ",i0,/)
  !
  1092 format ("    ERROR finding host cell for boundary face ",i10,", ",i10)
  !
end subroutine create_solver_geom_arrays_gmsh
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
  logical(lk) :: end_of_file
  character(len=80) :: text
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_gmsh_section_periodic"
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
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
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
  logical(lk) :: end_of_file
  character(len=80) :: text
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_gmsh_section_nodedata"
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
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
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
  logical(lk) :: end_of_file
  character(len=80) :: text
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_gmsh_section_elementdata"
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
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
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
  logical(lk) :: end_of_file
  character(len=80) :: text
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_gmsh_section_elementnodedata"
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
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  !
  ! Format Statements
  !
  1 format (/, &
            "******************** ERROR! *********************",/, &
            "  Pre-mature end-of-file while trying to read",/, &
            "  the 'ElementNodeData' section of the GMSH grid file!",/, &
            "******************** ERROR! *********************",/)
  2 format (/, &
            "*********************** ERROR! ************************",/, &
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
  logical(lk) :: end_of_file
  character(len=80) :: text
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_gmsh_interpolationscheme"
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
    write (error_message,1)
  else
    write (error_message,2)
  end if
  !
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
  !
  ! Format Statements
  !
  1 format ("Pre-mature end-of-file while trying to read the ", &
            "'InterpolationScheme' section of the GMSH grid file!")
  2 format ("While trying to read the 'InterpolationScheme' section of the ", &
            "GMSH grid file, an unknown and invalid section header was ", &
            "found before the required '$EndInterpolationScheme' line ", &
            "signifiying the end of the 'InterpolationScheme' section!")
  !
end subroutine read_gmsh_section_interpolationscheme
!
!###############################################################################
!
pure function npts_for_element_type(elem_type) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: elem_type
  !
  !.. Function Result ..
  integer :: return_value
  !
  !.. Local Scalars ..
  type(prop_t) :: prop
  !
continue
  !
  prop = element_properties(elem_type)
  return_value = prop%npts
  !
end function npts_for_element_type
!
!###############################################################################
!
pure function element_properties(elem_type) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: elem_type
  !
  !.. Function Result ..
  type(prop_t) :: return_value
  !
continue
  !
  select case (elem_type)
    case (1)
     !return_value = prop_t(npts=2, &
     !                      order=1, &
     !                      geom=Geom_Edge, &
     !                      elem_type=elem_type)
      return_value%npts = 2
      return_value%order = 1
      return_value%geom = Geom_Edge
      return_value%elem_type = elem_type
    case (2)
     !return_value = prop_t(npts=3, &
     !                      order=1, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type)
      return_value%npts = 3
      return_value%order = 1
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
    case (3)
     !return_value = prop_t(npts=4, &
     !                      order=1, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type)
      return_value%npts = 4
      return_value%order = 1
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
    case (4)
     !return_value = prop_t(npts=4, &
     !                      order=1, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=elem_type)
      return_value%npts = 4
      return_value%order = 1
      return_value%geom = Geom_Tetr
      return_value%elem_type = elem_type
    case (5)
     !return_value = prop_t(npts=8, &
     !                      order=1, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=elem_type)
      return_value%npts = 8
      return_value%order = 1
      return_value%geom = Geom_Hexa
      return_value%elem_type = elem_type
    case (6)
     !return_value = prop_t(npts=6, &
     !                      order=1, &
     !                      geom=Geom_Pris, &
     !                      elem_type=elem_type)
      return_value%npts = 6
      return_value%order = 1
      return_value%geom = Geom_Pris
      return_value%elem_type = elem_type
    case (7)
     !return_value = prop_t(npts=5, &
     !                      order=1, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=elem_type)
      return_value%npts = 5
      return_value%order = 1
      return_value%geom = Geom_Pyra
      return_value%elem_type = elem_type
    case (8)
     !return_value = prop_t(npts=3, &
     !                      order=2, &
     !                      geom=Geom_Edge, &
     !                      elem_type=elem_type)
      return_value%npts = 3
      return_value%order = 2
      return_value%geom = Geom_Edge
      return_value%elem_type = elem_type
    case (9)
     !return_value = prop_t(npts=6, &
     !                      order=2, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type)
      return_value%npts = 6
      return_value%order = 2
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
    case (10)
     !return_value = prop_t(npts=9, &
     !                      order=2, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type)
      return_value%npts = 9
      return_value%order = 2
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
    case (11)
     !return_value = prop_t(npts=10, &
     !                      order=2, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=elem_type)
      return_value%npts = 10
      return_value%order = 2
      return_value%geom = Geom_Tetr
      return_value%elem_type = elem_type
    case (12)
     !return_value = prop_t(npts=27, &
     !                      order=2, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=elem_type)
      return_value%npts = 27
      return_value%order = 2
      return_value%geom = Geom_Hexa
      return_value%elem_type = elem_type
    case (13)
     !return_value = prop_t(npts=18, &
     !                      order=2, &
     !                      geom=Geom_Pris, &
     !                      elem_type=elem_type)
      return_value%npts = 18
      return_value%order = 2
      return_value%geom = Geom_Pris
      return_value%elem_type = elem_type
    case (14)
     !return_value = prop_t(npts=14, &
     !                      order=2, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=elem_type)
      return_value%npts = 14
      return_value%order = 2
      return_value%geom = Geom_Pyra
      return_value%elem_type = elem_type
    case (15)
     !return_value = prop_t(npts=1, &
     !                      order=1, &
     !                      geom=Geom_Node, &
     !                      elem_type=elem_type)
      return_value%npts = 1
      return_value%order = 1
      return_value%geom = Geom_Node
      return_value%elem_type = elem_type
    case (16)
     !return_value = prop_t(npts=8, &
     !                      order=2, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 8
      return_value%order = 2
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (17)
     !return_value = prop_t(npts=20, &
     !                      order=2, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 20
      return_value%order = 2
      return_value%geom = Geom_Hexa
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (18)
     !return_value = prop_t(npts=15, &
     !                      order=2, &
     !                      geom=Geom_Pris, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 15
      return_value%order = 2
      return_value%geom = Geom_Pris
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (19)
     !return_value = prop_t(npts=13, &
     !                      order=2, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 13
      return_value%order = 2
      return_value%geom = Geom_Pyra
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (20)
     !return_value = prop_t(npts=9, &
     !                      order=3, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 9
      return_value%order = 3
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (21)
     !return_value = prop_t(npts=10, &
     !                      order=3, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type)
      return_value%npts = 10
      return_value%order = 3
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
    case (22)
     !return_value = prop_t(npts=12, &
     !                      order=4, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 12
      return_value%order = 4
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (23)
     !return_value = prop_t(npts=15, &
     !                      order=4, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type)
      return_value%npts = 15
      return_value%order = 4
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
    case (24)
     !return_value = prop_t(npts=15, &
     !                      order=5, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 15
      return_value%order = 5
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (25)
     !return_value = prop_t(npts=21, &
     !                      order=5, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type)
      return_value%npts = 21
      return_value%order = 5
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
    case (26)
     !return_value = prop_t(npts=4, &
     !                      order=3, &
     !                      geom=Geom_Edge, &
     !                      elem_type=elem_type)
      return_value%npts = 4
      return_value%order = 3
      return_value%geom = Geom_Edge
      return_value%elem_type = elem_type
    case (27)
     !return_value = prop_t(npts=5, &
     !                      order=4, &
     !                      geom=Geom_Edge, &
     !                      elem_type=elem_type)
      return_value%npts = 5
      return_value%order = 4
      return_value%geom = Geom_Edge
      return_value%elem_type = elem_type
    case (28)
     !return_value = prop_t(npts=6, &
     !                      order=5, &
     !                      geom=Geom_Edge, &
     !                      elem_type=elem_type)
      return_value%npts = 6
      return_value%order = 5
      return_value%geom = Geom_Edge
      return_value%elem_type = elem_type
    case (29)
     !return_value = prop_t(npts=20, &
     !                      order=3, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=elem_type)
      return_value%npts = 20
      return_value%order = 3
      return_value%geom = Geom_Tetr
      return_value%elem_type = elem_type
    case (30)
     !return_value = prop_t(npts=35, &
     !                      order=4, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=elem_type)
      return_value%npts = 35
      return_value%order = 4
      return_value%geom = Geom_Tetr
      return_value%elem_type = elem_type
    case (31)
     !return_value = prop_t(npts=56, &
     !                      order=5, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=elem_type)
      return_value%npts = 56
      return_value%order = 5
      return_value%geom = Geom_Tetr
      return_value%elem_type = elem_type
    case (32)
     !return_value = prop_t(npts=22, &
     !                      order=4, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 22
      return_value%order = 4
      return_value%geom = Geom_Tetr
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (33)
     !return_value = prop_t(npts=28, &
     !                      order=5, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 28
      return_value%order = 5
      return_value%geom = Geom_Tetr
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (36)
     !return_value = prop_t(npts=16, &
     !                      order=3, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type)
      return_value%npts = 16
      return_value%order = 3
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
    case (37)
     !return_value = prop_t(npts=25, &
     !                      order=4, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type)
      return_value%npts = 25
      return_value%order = 4
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
    case (38)
     !return_value = prop_t(npts=36, &
     !                      order=5, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type)
      return_value%npts = 36
      return_value%order = 5
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
    case (39)
     !return_value = prop_t(npts=12, &
     !                      order=3, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 12
      return_value%order = 3
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (40)
     !return_value = prop_t(npts=16, &
     !                      order=4, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 16
      return_value%order = 4
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (41)
     !return_value = prop_t(npts=20, &
     !                      order=5, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 20
      return_value%order = 5
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (42)
     !return_value = prop_t(npts=28, &
     !                      order=6, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type)
      return_value%npts = 28
      return_value%order = 6
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
    case (43)
     !return_value = prop_t(npts=36, &
     !                      order=7, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type)
      return_value%npts = 36
      return_value%order = 7
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
    case (44)
     !return_value = prop_t(npts=45, &
     !                      order=8, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type)
      return_value%npts = 45
      return_value%order = 8
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
    case (45)
     !return_value = prop_t(npts=55, &
     !                      order=9, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type)
      return_value%npts = 55
      return_value%order = 9
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
    case (46)
     !return_value = prop_t(npts=66, &
     !                      order=10, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type)
      return_value%npts = 66
      return_value%order = 10
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
    case (47)
     !return_value = prop_t(npts=49, &
     !                      order=6, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type)
      return_value%npts = 49
      return_value%order = 6
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
    case (48)
     !return_value = prop_t(npts=64, &
     !                      order=7, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type)
      return_value%npts = 64
      return_value%order = 7
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
    case (49)
     !return_value = prop_t(npts=81, &
     !                      order=8, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type)
      return_value%npts = 81
      return_value%order = 8
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
    case (50)
     !return_value = prop_t(npts=100, &
     !                      order=9, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type)
      return_value%npts = 100
      return_value%order = 9
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
    case (51)
     !return_value = prop_t(npts=121, &
     !                      order=10, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type)
      return_value%npts = 121
      return_value%order = 10
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
    case (52)
     !return_value = prop_t(npts=18, &
     !                      order=6, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 18
      return_value%order = 6
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (53)
     !return_value = prop_t(npts=21, &
     !                      order=7, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 21
      return_value%order = 7
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (54)
     !return_value = prop_t(npts=24, &
     !                      order=8, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 24
      return_value%order = 8
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (55)
     !return_value = prop_t(npts=27, &
     !                      order=9, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 27
      return_value%order = 9
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (56)
     !return_value = prop_t(npts=30, &
     !                      order=10, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 30
      return_value%order = 10
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (57)
     !return_value = prop_t(npts=24, &
     !                      order=6, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 24
      return_value%order = 6
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (58)
     !return_value = prop_t(npts=28, &
     !                      order=7, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 28
      return_value%order = 7
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (59)
     !return_value = prop_t(npts=32, &
     !                      order=8, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 32
      return_value%order = 8
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (60)
     !return_value = prop_t(npts=36, &
     !                      order=9, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 36
      return_value%order = 9
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (61)
     !return_value = prop_t(npts=40, &
     !                      order=10, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 40
      return_value%order = 10
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (62)
     !return_value = prop_t(npts=7, &
     !                      order=6, &
     !                      geom=Geom_Edge, &
     !                      elem_type=elem_type)
      return_value%npts = 7
      return_value%order = 6
      return_value%geom = Geom_Edge
      return_value%elem_type = elem_type
    case (63)
     !return_value = prop_t(npts=8, &
     !                      order=7, &
     !                      geom=Geom_Edge, &
     !                      elem_type=elem_type)
      return_value%npts = 8
      return_value%order = 7
      return_value%geom = Geom_Edge
      return_value%elem_type = elem_type
    case (64)
     !return_value = prop_t(npts=9, &
     !                      order=8, &
     !                      geom=Geom_Edge, &
     !                      elem_type=elem_type)
      return_value%npts = 9
      return_value%order = 8
      return_value%geom = Geom_Edge
      return_value%elem_type = elem_type
    case (65)
     !return_value = prop_t(npts=10, &
     !                      order=9, &
     !                      geom=Geom_Edge, &
     !                      elem_type=elem_type)
      return_value%npts = 10
      return_value%order = 9
      return_value%geom = Geom_Edge
      return_value%elem_type = elem_type
    case (66)
     !return_value = prop_t(npts=11, &
     !                      order=10, &
     !                      geom=Geom_Edge, &
     !                      elem_type=elem_type)
      return_value%npts = 11
      return_value%order = 10
      return_value%geom = Geom_Edge
      return_value%elem_type = elem_type
    case (71)
     !return_value = prop_t(npts=84, &
     !                      order=6, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=elem_type)
      return_value%npts = 84
      return_value%order = 6
      return_value%geom = Geom_Tetr
      return_value%elem_type = elem_type
    case (72)
     !return_value = prop_t(npts=120, &
     !                      order=7, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=elem_type)
      return_value%npts = 120
      return_value%order = 7
      return_value%geom = Geom_Tetr
      return_value%elem_type = elem_type
    case (73)
     !return_value = prop_t(npts=165, &
     !                      order=8, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=elem_type)
      return_value%npts = 165
      return_value%order = 8
      return_value%geom = Geom_Tetr
      return_value%elem_type = elem_type
    case (74)
     !return_value = prop_t(npts=220, &
     !                      order=9, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=elem_type)
      return_value%npts = 220
      return_value%order = 9
      return_value%geom = Geom_Tetr
      return_value%elem_type = elem_type
    case (75)
     !return_value = prop_t(npts=286, &
     !                      order=10, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=elem_type)
      return_value%npts = 286
      return_value%order = 10
      return_value%geom = Geom_Tetr
      return_value%elem_type = elem_type
    case (84)
     !return_value = prop_t(npts=1, &
     !                      order=0, &
     !                      geom=Geom_Edge, &
     !                      elem_type=elem_type, &
     !                      warning=true)
      return_value%npts = 1
      return_value%order = 0
      return_value%geom = Geom_Edge
      return_value%elem_type = elem_type
      return_value%warning = true
    case (85)
     !return_value = prop_t(npts=1, &
     !                      order=0, &
     !                      geom=Geom_Tria, &
     !                      elem_type=elem_type, &
     !                      warning=true)
      return_value%npts = 1
      return_value%order = 0
      return_value%geom = Geom_Tria
      return_value%elem_type = elem_type
      return_value%warning = true
    case (86)
     !return_value = prop_t(npts=1, &
     !                      order=0, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type, &
     !                      warning=true)
      return_value%npts = 1
      return_value%order = 0
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
      return_value%warning = true
    case (87)
     !return_value = prop_t(npts=1, &
     !                      order=0, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=elem_type, &
     !                      warning=true)
      return_value%npts = 1
      return_value%order = 0
      return_value%geom = Geom_Tetr
      return_value%elem_type = elem_type
      return_value%warning = true
    case (88)
     !return_value = prop_t(npts=1, &
     !                      order=0, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=elem_type, &
     !                      warning=true)
      return_value%npts = 1
      return_value%order = 0
      return_value%geom = Geom_Hexa
      return_value%elem_type = elem_type
      return_value%warning = true
    case (89)
     !return_value = prop_t(npts=1, &
     !                      order=0, &
     !                      geom=Geom_Pris, &
     !                      elem_type=elem_type, &
     !                      warning=true)
      return_value%npts = 1
      return_value%order = 0
      return_value%geom = Geom_Pris
      return_value%elem_type = elem_type
      return_value%warning = true
    case (90)
     !return_value = prop_t(npts=40, &
     !                      order=3, &
     !                      geom=Geom_Pris, &
     !                      elem_type=elem_type)
      return_value%npts = 40
      return_value%order = 3
      return_value%geom = Geom_Pris
      return_value%elem_type = elem_type
    case (91)
     !return_value = prop_t(npts=75, &
     !                      order=4, &
     !                      geom=Geom_Pris, &
     !                      elem_type=elem_type)
      return_value%npts = 75
      return_value%order = 4
      return_value%geom = Geom_Pris
      return_value%elem_type = elem_type
    case (92)
     !return_value = prop_t(npts=64, &
     !                      order=3, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=elem_type)
      return_value%npts = 64
      return_value%order = 3
      return_value%geom = Geom_Hexa
      return_value%elem_type = elem_type
    case (93)
     !return_value = prop_t(npts=125, &
     !                      order=4, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=elem_type)
      return_value%npts = 125
      return_value%order = 4
      return_value%geom = Geom_Hexa
      return_value%elem_type = elem_type
    case (94)
     !return_value = prop_t(npts=216, &
     !                      order=5, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=elem_type)
      return_value%npts = 216
      return_value%order = 5
      return_value%geom = Geom_Hexa
      return_value%elem_type = elem_type
    case (95)
     !return_value = prop_t(npts=343, &
     !                      order=6, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=elem_type, &
     !                      warning=true)
      return_value%npts = 343
      return_value%order = 6
      return_value%geom = Geom_Hexa
      return_value%elem_type = elem_type
      return_value%warning = true
    case (96)
     !return_value = prop_t(npts=512, &
     !                      order=7, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=elem_type, &
     !                      warning=true)
      return_value%npts = 512
      return_value%order = 7
      return_value%geom = Geom_Hexa
      return_value%elem_type = elem_type
      return_value%warning = true
    case (97)
     !return_value = prop_t(npts=729, &
     !                      order=8, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=elem_type, &
     !                      warning=true)
      return_value%npts = 729
      return_value%order = 8
      return_value%geom = Geom_Hexa
      return_value%elem_type = elem_type
      return_value%warning = true
    case (98)
     !return_value = prop_t(npts=1000, &
     !                      order=9, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=elem_type, &
     !                      warning=true)
      return_value%npts = 1000
      return_value%order = 9
      return_value%geom = Geom_Hexa
      return_value%elem_type = elem_type
      return_value%warning = true
    case (106)
     !return_value = prop_t(npts=126, &
     !                      order=5, &
     !                      geom=Geom_Pris, &
     !                      elem_type=elem_type)
      return_value%npts = 126
      return_value%order = 5
      return_value%geom = Geom_Pris
      return_value%elem_type = elem_type
    case (107)
     !return_value = prop_t(npts=196, &
     !                      order=6, &
     !                      geom=Geom_Pris, &
     !                      elem_type=elem_type, &
     !                      warning=true)
      return_value%npts = 196
      return_value%order = 6
      return_value%geom = Geom_Pris
      return_value%elem_type = elem_type
      return_value%warning = true
    case (108)
     !return_value = prop_t(npts=288, &
     !                      order=7, &
     !                      geom=Geom_Pris, &
     !                      elem_type=elem_type, &
     !                      warning=true)
      return_value%npts = 288
      return_value%order = 7
      return_value%geom = Geom_Pris
      return_value%elem_type = elem_type
      return_value%warning = true
    case (109)
     !return_value = prop_t(npts=405, &
     !                      order=8, &
     !                      geom=Geom_Pris, &
     !                      elem_type=elem_type, &
     !                      warning=true)
      return_value%npts = 405
      return_value%order = 8
      return_value%geom = Geom_Pris
      return_value%elem_type = elem_type
      return_value%warning = true
    case (110)
     !return_value = prop_t(npts=550, &
     !                      order=9, &
     !                      geom=Geom_Pris, &
     !                      elem_type=elem_type, &
     !                      warning=true)
      return_value%npts = 550
      return_value%order = 9
      return_value%geom = Geom_Pris
      return_value%elem_type = elem_type
      return_value%warning = true
    case (118)
     !return_value = prop_t(npts=30, &
     !                      order=3, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=elem_type)
      return_value%npts = 30
      return_value%order = 3
      return_value%geom = Geom_Pyra
      return_value%elem_type = elem_type
    case (119)
     !return_value = prop_t(npts=55, &
     !                      order=4, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=elem_type)
      return_value%npts = 55
      return_value%order = 4
      return_value%geom = Geom_Pyra
      return_value%elem_type = elem_type
    case (120)
     !return_value = prop_t(npts=91, &
     !                      order=5, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=elem_type)
      return_value%npts = 91
      return_value%order = 5
      return_value%geom = Geom_Pyra
      return_value%elem_type = elem_type
    case (121)
     !return_value = prop_t(npts=140, &
     !                      order=6, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=elem_type)
      return_value%npts = 140
      return_value%order = 6
      return_value%geom = Geom_Pyra
      return_value%elem_type = elem_type
    case (122)
     !return_value = prop_t(npts=204, &
     !                      order=7, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=elem_type)
      return_value%npts = 204
      return_value%order = 7
      return_value%geom = Geom_Pyra
      return_value%elem_type = elem_type
    case (123)
     !return_value = prop_t(npts=285, &
     !                      order=8, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=elem_type)
      return_value%npts = 285
      return_value%order = 8
      return_value%geom = Geom_Pyra
      return_value%elem_type = elem_type
    case (124)
     !return_value = prop_t(npts=385, &
     !                      order=9, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=elem_type)
      return_value%npts = 385
      return_value%order = 9
      return_value%geom = Geom_Pyra
      return_value%elem_type = elem_type
    case (137)
     !return_value = prop_t(npts=16, &
     !                      order=3, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=elem_type, &
     !                      family=Edge_Defined)
      return_value%npts = 16
      return_value%order = 3
      return_value%geom = Geom_Tetr
      return_value%elem_type = elem_type
      return_value%family = Edge_Defined
      return_value%warning = true
    case (132)
     !return_value = prop_t(npts=1, &
     !                      order=0, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=elem_type, &
     !                      warning=true)
      return_value%npts = 1
      return_value%order = 0
      return_value%geom = Geom_Pyra
      return_value%elem_type = elem_type
      return_value%warning = true
    case default
     !return_value = prop_t(npts=-1, &
     !                      order=-1, &
     !                      geom=Geom_Unknown, &
     !                      elem_type=elem_type)
      return_value%npts = -1
      return_value%order = -1
      return_value%geom = Geom_Unknown
      return_value%elem_type = elem_type
  end select
  !
end function element_properties
!
!###############################################################################
!
elemental function reorder_gmsh_pts(elem) result(return_value)
  !
  !.. Formal Arguments ..
  type(elem_t), intent(in) :: elem
  !
  !.. Function Result ..
  type(elem_t) :: return_value
  !
continue
  !
  return_value = elem
  !
  if (elem%prop%family == Edge_Defined) then
    !
    ! This is an incomplete element in the edge-defined/serendipity family
    !
    select case (elem%prop%geom)
      case (Geom_Edge)
        return_value%pts = reorder_gmsh_edge_pts(elem%pts)
      case (Geom_Tria)
        return_value%pts = reorder_gmsh_edge_defined_tria_pts(elem%pts)
      case (Geom_Quad)
        return_value%pts = reorder_gmsh_edge_defined_quad_pts(elem%pts)
      case (Geom_Tetr)
        return_value%pts = reorder_gmsh_edge_defined_tetr_pts(elem%pts)
      case (Geom_Pyra)
        return_value%pts = reorder_gmsh_edge_defined_pyra_pts(elem%pts)
      case (Geom_Pris)
        return_value%pts = reorder_gmsh_edge_defined_pris_pts(elem%pts)
      case (Geom_Hexa)
        return_value%pts = reorder_gmsh_edge_defined_hexa_pts(elem%pts)
    end select
    !
  else
    !
    ! This is a complete element with internal node points
    !
    select case (elem%prop%geom)
      case (Geom_Edge)
        return_value%pts = reorder_gmsh_edge_pts(elem%pts)
      case (Geom_Tria)
        return_value%pts = reorder_gmsh_tria_pts(elem%pts)
      case (Geom_Quad)
        return_value%pts = reorder_gmsh_quad_pts(elem%pts)
      case (Geom_Tetr)
        return_value%pts = reorder_gmsh_tetr_pts(elem%pts)
      case (Geom_Pyra)
        return_value%pts = reorder_gmsh_pyra_pts(elem%pts)
      case (Geom_Pris)
        return_value%pts = reorder_gmsh_pris_pts(elem%pts)
      case (Geom_Hexa)
        return_value%pts = reorder_gmsh_hexa_pts(elem%pts)
    end select
    !
  end if
  !
end function reorder_gmsh_pts
!
!###############################################################################
!
pure function reorder_gmsh_edge_pts(pts) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Scalars ..
  integer :: n,i
continue
  !
  n = size(pts)
  !
  return_value(1) = pts(1)
  return_value(n) = pts(2)
  !
  do i = 2,n-1
    return_value(i) = pts(i+1)
  end do
  !
end function reorder_gmsh_edge_pts
!
!###############################################################################
!
pure function reorder_gmsh_tria_pts(pts) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: p1 = 3
  integer, parameter :: tria1(1:p1) = [1,2,3]
  !
  integer, parameter :: p2 = 6
  integer, parameter :: tria2(1:p2) = [1,4,2,6,5,3]
  !
  integer, parameter :: p3 = 10
  integer, parameter :: tria3(1:p3) = [1,4,5,2,9,10,6,8,7,3]
  !
  integer, parameter :: p4 = 15
  integer, parameter :: tria4(1:p4) = [1,4,5,6,2,12,13,14,7,11,15,8,10,9,3]
  !
  integer, parameter :: p5 = 21
  integer, parameter :: tria5(1:p5) = [ 1, 4, 5, 6, 7, 2,15, &
                                       16,19,17, 8,14,21,20, &
                                        9,13,18,10,12,11, 3]
  !
  integer, parameter :: p6 = 28
  integer, parameter :: tria6(1:p6)  = [ 1, 4, 5, 6, 7, 8, 2, &
                                        18,19,22,23,20, 9,17, &
                                        27,28,24,10,16,26,25, &
                                        11,15,21,12,14,13, 3]
  !
  integer, parameter :: p7 = 36
  integer, parameter :: tria7(1:p7)  = [ 1, 4, 5, 6, 7, 8, 9, 2,21,22,25,26, &
                                        27,23,10,20,33,34,35,28,11,19,32,36, &
                                        29,12,18,31,30,13,17,24,14,16,15, 3]
  !
  integer, parameter :: p8 = 45
  integer, parameter :: tria8(1:p8)  = [ 1, 4, 5, 6, 7, 8, 9,10, 2, &
                                        24,25,28,29,30,31,26,11,23, &
                                        39,40,43,41,32,12,22,38,45, &
                                        44,33,13,21,37,42,34,14,20, &
                                        36,35,15,19,27,16,18,17, 3]
  !
  integer, parameter :: p9 = 55
  integer, parameter :: tria9(1:p9)  = [ 1, 4, 5, 6, 7, 8, 9,10,11, 2,27, &
                                        28,31,32,33,34,35,29,12,26,45,46, &
                                        49,50,47,36,13,25,44,54,55,51,37, &
                                        14,24,43,53,52,38,15,23,42,48,39, &
                                        16,22,41,40,17,21,30,18,20,19, 3]
  !
  integer, parameter :: p10 = 66
  integer, parameter :: tria10(1:p10) = [ 1, 4, 5, 6, 7, 8, 9,10,11,12, 2, &
                                         30,31,34,35,36,37,38,39,32,13,29, &
                                         51,52,55,56,57,53,40,14,28,50,63, &
                                         64,65,58,41,15,27,49,62,66,59,42, &
                                         16,26,48,61,60,43,17,25,47,54,44, &
                                         18,24,46,45,19,23,33,20,22,21, 3]
continue
  !
  select case (size(pts))
    case (p1)
      return_value = pts( tria1(:) )
    case (p2)
      return_value = pts( tria2(:) )
    case (p3)
      return_value = pts( tria3(:) )
    case (p4)
      return_value = pts( tria4(:) )
    case (p5)
      return_value = pts( tria5(:) )
    case (p6)
      return_value = pts( tria6(:) )
    case (p7)
      return_value = pts( tria7(:) )
    case (p8)
      return_value = pts( tria8(:) )
    case (p9)
      return_value = pts( tria9(:) )
    case (p10)
      return_value = pts( tria10(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_gmsh_tria_pts
!
!###############################################################################
!
pure function reorder_gmsh_quad_pts(pts) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: p1 = 4
  integer, parameter :: quad1(1:p1) = [1,2, &
                                       4,3]
  !
  integer, parameter :: p2 = 9
  integer, parameter :: quad2(1:p2) = [1,5,2, &
                                       8,9,6, &
                                       4,7,3]
  !
  integer, parameter :: p3 = 16
  integer, parameter :: quad3(1:p3) = [ 1, 5, 6, 2, &
                                       12,13,14, 7, &
                                       11,16,15, 8, &
                                        4,10, 9, 3]
  !
  integer, parameter :: p4 = 25
  integer, parameter :: quad4(1:p4) = [ 1, 5, 6, 7, 2, &
                                       16,17,21,18, 8, &
                                       15,24,25,22, 9, &
                                       14,20,23,19,10, &
                                        4,13,12,11, 3]
  !
  integer, parameter :: p5 = 36
  integer, parameter :: quad5(1:p5) = [ 1, 5, 6, 7, 8, 2, &
                                       20,21,25,26,22, 9, &
                                       19,32,33,34,27,10, &
                                       18,31,36,35,28,11, &
                                       17,24,30,29,23,12, &
                                        4,16,15,14,13, 3]
  !
  integer, parameter :: p6 = 49
  integer, parameter :: quad6(1:p6) = [ 1, 5, 6, 7, 8, 9, 2, &
                                       24,25,29,30,31,26,10, &
                                       23,40,41,45,42,32,11, &
                                       22,39,48,49,46,33,12, &
                                       21,38,44,47,43,34,13, &
                                       20,28,37,36,35,27,14, &
                                        4,19,18,17,16,15, 3]
  !
  integer, parameter :: p7 = 64
  integer, parameter :: quad7(1:p7) = [ 1, 5, 6, 7, 8, 9,10, 2, &
                                       28,29,33,34,35,36,30,11, &
                                       27,48,49,53,54,50,37,12, &
                                       26,47,60,61,62,55,38,13, &
                                       25,46,59,64,63,56,39,14, &
                                       24,45,52,58,57,51,40,15, &
                                       23,32,44,43,42,41,31,16, &
                                        4,22,21,20,19,18,17, 3]
  !
  integer, parameter :: p8 = 81
  integer, parameter :: quad8(1:p8) = [ 1, 5, 6, 7, 8, 9,10,11, 2, &
                                       32,33,37,38,39,40,41,34,12, &
                                       31,56,57,61,62,63,58,42,13, &
                                       30,55,72,73,77,74,64,43,14, &
                                       29,54,71,80,81,78,65,44,15, &
                                       28,53,70,76,79,75,66,45,16, &
                                       27,52,60,69,68,67,59,46,17, &
                                       26,36,51,50,49,48,47,35,18, &
                                        4,25,24,23,22,21,20,19, 3]
  !
  integer, parameter :: p9 = 100
  integer, parameter :: quad9(1:p9) = [ 1, 5, 6, 7,  8, 9,10,11,12, 2, &
                                        36,37,41,42, 43,44,45,46,38,13, &
                                        35,64,65,69, 70,71,72,66,47,14, &
                                        34,63,84,85, 89,90,86,73,48,15, &
                                        33,62,83,96, 97,98,91,74,49,16, &
                                        32,61,82,95,100,99,92,75,50,17, &
                                        31,60,81,88, 94,93,87,76,51,18, &
                                        30,59,68,80, 79,78,77,67,52,19, &
                                        29,40,58,57, 56,55,54,53,39,20, &
                                         4,28,27,26, 25,24,23,22,21, 3]
  !
  integer, parameter :: p10 = 121
  integer, parameter :: quad10(1:p10) = &
    [  1,  5,  6,  7,  8,  9, 10, 11, 12, 13,  2, &
      40, 41, 45, 46, 47, 48, 49, 50, 51, 42, 14, &
      39, 72, 73, 77, 78, 79, 80, 81, 74, 52, 15, &
      38, 71, 96, 97,101,102,103, 98, 82, 53, 16, &
      37, 70, 95,112,113,117,114,104, 83, 54, 17, &
      36, 69, 94,111,120,121,118,105, 84, 55, 18, &
      35, 68, 93,110,116,119,115,106, 85, 56, 19, &
      34, 67, 92,100,109,108,107, 99, 86, 57, 20, &
      33, 66, 76, 91, 90, 89, 88, 87, 75, 58, 21, &
      32, 44, 65, 64, 63, 62, 61, 60, 59, 43, 22, &
       4, 31, 30, 29, 28, 27, 26, 25, 24, 23,  3]
continue
  !
  select case (size(pts))
    case (p1)
      return_value = pts( quad1(:) )
    case (p2)
      return_value = pts( quad2(:) )
    case (p3)
      return_value = pts( quad3(:) )
    case (p4)
      return_value = pts( quad4(:) )
    case (p5)
      return_value = pts( quad5(:) )
    case (p6)
      return_value = pts( quad6(:) )
    case (p7)
      return_value = pts( quad7(:) )
    case (p8)
      return_value = pts( quad8(:) )
    case (p9)
      return_value = pts( quad9(:) )
    case (p10)
      return_value = pts( quad10(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_gmsh_quad_pts
!
!###############################################################################
!
pure function reorder_gmsh_tetr_pts(pts) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: p1 = 4
  integer, parameter :: tetr1(1:p1) = [1,4,2,3]
  !
  integer, parameter :: p2 = 10
  integer, parameter :: tetr2(1:p2) = [1,8,4,5,10,2,7,9,6,3]
  !
  integer, parameter :: p3 = 20
  integer, parameter :: tetr3(1:p3) = [ 1,12,11, 4, 5,18,15, 6,16, 2, &
                                       10,19,13,17,20, 7, 9,14, 8, 3]
  !
  integer, parameter :: p4 = 35
  integer, parameter :: tetr4(1:p4) = [ 1,16,15,14, 4, 5,26, &
                                       28,20, 6,27,21, 7,22, &
                                        2,13,29,30,17,23,35, &
                                       32,25,33, 8,12,31,18, &
                                       24,34, 9,11,19,10, 3]
  !
  integer, parameter :: p5 = 56
  integer, parameter :: tetr5(1:p5) = &
    [ 1,20,19,18,17, 4, 5,35,40,37,25, 6,38,39, &
     26, 7,36,27, 8,28, 2,16,41,44,42,21,29,53, &
     56,47,34,54,50,31,48, 9,15,46,45,22,32,55, &
     52,33,51,10,14,43,23,30,49,11,13,24,12, 3]
  !
  integer, parameter :: p6 = 84
  integer, parameter :: tetr6(1:p6) = &
    [ 1,24,23,22,21,20, 4, 5,45,53,52,47,30, 6,48,54,51,31, 7,49,50, &
     32, 8,46,33, 9,34, 2,19,55,58,59,56,25,35,75,82,78,65,43,79,84, &
     68,42,76,69,37,66,10,18,63,64,60,26,38,81,83,73,44,80,74,41,70, &
     11,17,62,61,27,39,77,72,40,71,12,16,57,28,36,67,13,15,29,14, 3]
  !
  integer, parameter :: p7 = 120
  integer, parameter :: tetr7(1:p7) = &
    [  1, 28, 27, 26, 25, 24, 23,  4,  5, 56, 67, 66, 65, 58, 35, &
       6, 59, 68, 70, 64, 36,  7, 60, 69, 63, 37,  8, 61, 62, 38, &
       9, 57, 39, 10, 40,  2, 22, 71, 74, 75, 76, 72, 29, 41,101, &
     112,111,104, 86, 52,105,118,115, 89, 51,106,116, 90, 50,102, &
      91, 43, 87, 11, 21, 82, 83, 84, 77, 30, 44,110,119,113, 97, &
      53,117,120, 98, 55,107, 99, 49, 92, 12, 20, 81, 85, 78, 31, &
      45,109,114, 96, 54,108,100, 48, 93, 13, 19, 80, 79, 32, 46, &
     103, 95, 47, 94, 14, 18, 73, 33, 42, 88, 15, 17, 34, 16,  3]
  !
  integer, parameter :: p8 = 165
  integer, parameter :: tetr8(1:p8) = &
    [  1, 32, 31, 30, 29, 28, 27, 26,  4,  5, 68, 82, 81, 80, 79, &
      70, 40,  6, 71, 83, 88, 85, 78, 41,  7, 72, 86, 87, 77, 42, &
       8, 73, 84, 76, 43,  9, 74, 75, 44, 10, 69, 45, 11, 46,  2, &
      25, 89, 92, 93, 94, 95, 90, 33, 47,131,146,145,144,134,110, &
      61,135,156,158,150,113, 60,136,157,151,114, 59,137,152,115, &
      58,132,116, 49,111, 12, 24,103,104,107,105, 96, 34, 50,143, &
     159,160,147,124, 62,153,165,162,125, 67,155,163,128, 64,138, &
     126, 57,117, 13, 23,102,109,108, 97, 35, 51,142,161,148,123, &
      65,154,164,130, 66,139,129, 56,118, 14, 22,101,106, 98, 36, &
      52,141,149,122, 63,140,127, 55,119, 15, 21,100, 99, 37, 53, &
     133,121, 54,120, 16, 20, 91, 38, 48,112, 17, 19, 39, 18,  3]
  !
  integer, parameter :: p9 = 220
  integer, parameter :: tetr9(1:p9) = &
    [  1, 36, 35, 34, 33, 32, 31, 30, 29,  4, &
       5, 81, 98, 97, 96, 95, 94, 83, 45,  6, &
      84, 99,107,106,101, 93, 46,  7, 85,102, &
     108,105, 92, 47,  8, 86,103,104, 91, 48, &
       9, 87,100, 90, 49, 10, 88, 89, 50, 11, &
      82, 51, 12, 52,  2, 28,109,112,113,114, &
     115,116,110, 37, 53,165,184,183,182,181, &
     168,137, 70,169,199,204,201,189,140, 69, &
     170,202,203,190,141, 68,171,200,191,142, &
      67,172,192,143, 66,166,144, 55,138, 13, &
      27,126,127,130,131,128,117, 38, 56,180, &
     205,208,206,185,154, 71,193,217,220,211, &
     155, 79,198,218,214,158, 78,195,212,159, &
      73,173,156, 65,145, 14, 26,125,135,136, &
     132,118, 39, 57,179,210,209,186,153, 74, &
     196,219,216,163, 80,197,215,164, 77,174, &
     160, 64,146, 15, 25,124,134,133,119, 40, &
      58,178,207,187,152, 75,194,213,162, 76, &
     175,161, 63,147, 16, 24,123,129,120, 41, &
      59,177,188,151, 72,176,157, 62,148, 17, &
      23,122,121, 42, 60,167,150, 61,149, 18, &
      22,111, 43, 54,139, 19, 21, 44, 20,  3]
  !
  integer, parameter :: p10 = 286
  integer, parameter :: tetr10(1:p10) = &
    [  1, 40, 39, 38, 37, 36, 35, 34, 33, 32,  4,  5, 95, &
     115,114,113,112,111,110, 97, 50,  6, 98,116,127,126, &
     125,118,109, 51,  7, 99,119,128,130,124,108, 52,  8, &
     100,120,129,123,107, 53,  9,101,121,122,106, 54, 10, &
     102,117,105, 55, 11,103,104, 56, 12, 96, 57, 13, 58, &
       2, 31,131,134,135,136,137,138,139,132, 41, 59,203, &
     226,225,224,223,222,206,167, 79,207,247,255,254,249, &
     232,170, 78,208,250,256,253,233,171, 77,209,251,252, &
     234,172, 76,210,248,235,173, 75,211,236,174, 74,204, &
     175, 61,168, 14, 30,151,152,155,156,157,153,140, 42, &
      62,221,257,260,261,258,227,187, 80,237,277,284,280, &
     267,188, 91,245,281,286,270,191, 90,244,278,271,192, &
      89,239,268,193, 82,212,189, 73,176, 15, 29,150,163, &
     164,165,158,141, 43, 63,220,265,266,262,228,186, 83, &
     240,283,285,275,199, 92,246,282,276,200, 94,243,272, &
     201, 88,213,194, 72,177, 16, 28,149,162,166,159,142, &
      44, 64,219,264,263,229,185, 84,241,279,274,198, 93, &
     242,273,202, 87,214,195, 71,178, 17, 27,148,161,160, &
     143, 45, 65,218,259,230,184, 85,238,269,197, 86,215, &
     196, 70,179, 18, 26,147,154,144, 46, 66,217,231,183, &
      81,216,190, 69,180, 19, 25,146,145, 47, 67,205,182, &
      68,181, 20, 24,133, 48, 60,169, 21, 23, 49, 22,  3]
  !
continue
  !
  select case (size(pts))
    case (p1)
      return_value = pts( tetr1(:) )
    case (p2)
      return_value = pts( tetr2(:) )
    case (p3)
      return_value = pts( tetr3(:) )
    case (p4)
      return_value = pts( tetr4(:) )
    case (p5)
      return_value = pts( tetr5(:) )
    case (p6)
      return_value = pts( tetr6(:) )
    case (p7)
      return_value = pts( tetr7(:) )
    case (p8)
      return_value = pts( tetr8(:) )
    case (p9)
      return_value = pts( tetr9(:) )
    case (p10)
      return_value = pts( tetr10(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_gmsh_tetr_pts
!
!###############################################################################
!
pure function reorder_gmsh_pyra_pts(pts) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: p1 = 5
  integer, parameter :: pyra1(1:p1) = [1,2,4,3,5]
  !
  integer, parameter :: p2 = 14
  integer, parameter :: pyra2(1:p2) = [1,6,2,7,14,9,4,11,3,8,10,13,12,5]
  !
  integer, parameter :: p3 = 30
  integer, parameter :: pyra3(1:p3) = [ 1, 6, 7, 2, 8,26,29,12, 9,27, &
                                       28,13, 4,17,16, 3,10,22,14,23, &
                                       30,24,20,25,18,11,15,21,19, 5]
  !
  integer, parameter :: p4 = 55
  integer, parameter :: pyra4(1:p4) = [ 1, 6, 7, 8, 2, 9,42,49,45,15,10, &
                                       46,50,48,16,11,43,47,44,17, 4,23, &
                                       22,21, 3,12,30,31,18,34,51,52,36, &
                                       33,54,53,37,27,40,39,24,13,32,19, &
                                       35,55,38,28,41,25,14,20,29,26, 5]
  !
  integer, parameter :: p5 = 91
  integer, parameter :: p6 = 140
  integer, parameter :: p7 = 204
  integer, parameter :: p8 = 285
  integer, parameter :: p9 = 385
  !
continue
  !
  select case (size(pts))
    case (p1)
      return_value = pts( pyra1(:) )
    case (p2)
      return_value = pts( pyra2(:) )
    case (p3)
      return_value = pts( pyra3(:) )
    case (p4)
      return_value = pts( pyra4(:) )
   !case (p5)
   !  return_value = pts( pyra5(:) )
   !case (p6)
   !  return_value = pts( pyra6(:) )
   !case (p7)
   !  return_value = pts( pyra7(:) )
   !case (p8)
   !  return_value = pts( pyra8(:) )
   !case (p9)
   !  return_value = pts( pyra9(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_gmsh_pyra_pts
!
!###############################################################################
!
pure function reorder_gmsh_pris_pts(pts) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: p1 = 6
  integer, parameter :: pris1(1:p1) = [1,2,3,4,5,6]
  !
  integer, parameter :: p2 = 18
  integer, parameter :: pris2(1:p2) = [ 1, 7, 2, 8,10, 3, 9,16,11, &
                                       17,18,12, 4,13, 5,14,15, 6]
  !
  integer, parameter :: p3 = 40
  integer, parameter :: pris3(1:p3) = [ 1, 7, 8, 2, 9,25,13,10,14, 3, &
                                       11,27,28,15,31,39,35,34,36,17, &
                                       12,30,29,16,32,40,38,33,37,18, &
                                        4,19,20, 5,21,26,23,22,24, 6]
  !
  integer, parameter :: p4 = 75
  integer, parameter :: pris4(1:p4) = &
    [ 1, 7, 8, 9, 2,10,34,36,16,11,35,17,12,18, 3, &
     13,40,44,41,19,49,67,70,58,56,73,62,52,59,22, &
     14,47,48,45,20,53,69,72,65,57,75,66,55,63,23, &
     15,43,46,42,21,50,68,71,61,54,74,64,51,60,24, &
      4,25,26,27, 5,28,37,38,31,29,39,32,30,33, 6]
  !
  integer, parameter :: p5 = 126
  integer, parameter :: pris5(1:p5) = &
    [  1,  7,  8,  9, 10,  2, 11, 43, 48, 45, 19, 12, 46, 47, 20, 13, 44, 21, &
      14, 22,  3, 15, 55, 59, 60, 56, 23, 71,103,115,107, 87, 82,123,119, 91, &
      81,111, 92, 74, 88, 27, 16, 66, 67, 68, 61, 24, 75,105,117,109, 98, 83, &
     125,121, 99, 86,113,100, 80, 93, 28, 17, 65, 70, 69, 62, 25, 76,106,118, &
     110, 97, 84,126,122,102, 85,114,101, 79, 94, 29, 18, 58, 64, 63, 57, 26, &
      72,104,116,108, 90, 77,124,120, 96, 78,112, 95, 73, 89, 30,  4, 31, 32, &
      33, 34,  5, 35, 49, 52, 50, 39, 36, 54, 53, 40, 37, 51, 41, 38, 42,  6]
  !
  integer, parameter :: p6 = 196
  integer, parameter :: p7 = 288
  integer, parameter :: p8 = 405
  integer, parameter :: p9 = 550
  !
continue
  !
  select case (size(pts))
    case (p1)
      return_value = pts( pris1(:) )
    case (p2)
      return_value = pts( pris2(:) )
    case (p3)
      return_value = pts( pris3(:) )
    case (p4)
      return_value = pts( pris4(:) )
    case (p5)
      return_value = pts( pris5(:) )
   !case (p6)
   !  return_value = pts( pris6(:) )
   !case (p7)
   !  return_value = pts( pris7(:) )
   !case (p8)
   !  return_value = pts( pris8(:) )
   !case (p9)
   !  return_value = pts( pris9(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_gmsh_pris_pts
!
!###############################################################################
!
pure function reorder_gmsh_hexa_pts(pts) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: p1 = 8
  integer, parameter :: hexa1(1:p1) = [1,2,4,3,5,6,8,7]
  !
  integer, parameter :: p2 = 27
  integer, parameter :: hexa2(1:p2) = [ 1, 9, 2,10,21,12, 4,14, 3, &
                                       11,22,13,23,27,24,16,25,15, &
                                        5,17, 6,18,26,19, 8,20, 7]
  !
  integer, parameter :: p3 = 64
  integer, parameter :: hexa3(1:p3) = &
    [ 1, 9,10, 2,11,33,36,15,12,34,35,16, 4,20,19, 3, &
     13,37,38,17,41,57,58,45,44,60,59,46,23,50,49,21, &
     14,40,39,18,42,61,62,48,43,64,63,47,24,51,52,22, &
      5,25,26, 6,27,53,54,29,28,56,55,30, 8,32,31, 7]
  !
  integer, parameter :: p4 = 125
  integer, parameter :: hexa4(1:p4) = &
    [  1,  9, 10, 11,  2, 12, 45, 52, 48, 18, 13, 49, 53, 51, 19, 14, 46, 50, &
      47, 20,  4, 26, 25, 24,  3, 15, 54, 58, 55, 21, 63, 99,107,100, 72, 70, &
     108,119,110, 76, 66,102,112,101, 73, 30, 82, 85, 81, 27, 16, 61, 62, 59, &
      22, 67,109,120,111, 79, 71,121,125,122, 80, 69,114,123,113, 77, 31, 86, &
      89, 88, 28, 17, 57, 60, 56, 23, 64,103,115,104, 75, 68,116,124,117, 78, &
      65,106,118,105, 74, 32, 83, 87, 84, 29,  5, 33, 34, 35,  6, 36, 90, 94, &
      91, 39, 37, 97, 98, 95, 40, 38, 93, 96, 92, 41,  8, 44, 43, 42,  7]
  !
  integer, parameter :: p5 = 216
  integer, parameter :: hexa5(1:p5) = &
    [  1,  9, 10, 11, 12,  2, 13, 57, 68, 67, 60, 21, 14, 61, 69, 72, 66, 22, &
      15, 62, 70, 71, 65, 23, 16, 58, 63, 64, 59, 24,  4, 32, 31, 30, 29,  3, &
      17, 73, 77, 78, 74, 25, 89,153,161,162,154,105,100,163,185,188,167,109, &
      99,164,186,187,168,110, 92,156,172,171,155,106, 37,122,126,125,121, 33, &
      18, 84, 85, 86, 79, 26, 93,165,189,190,169,116,101,193,209,210,197,117, &
     104,196,212,211,198,118, 98,175,202,201,173,111, 38,127,134,133,132, 34, &
      19, 83, 88, 87, 80, 27, 94,166,192,191,170,115,102,194,213,214,200,120, &
     103,195,216,215,199,119, 97,176,203,204,174,112, 39,128,135,136,131, 35, &
      20, 76, 82, 81, 75, 28, 90,157,177,178,158,108, 95,179,205,206,181,114, &
      96,180,208,207,182,113, 91,160,184,183,159,107, 40,123,129,130,124, 36, &
       5, 41, 42, 43, 44,  6, 45,137,141,142,138, 49, 46,148,149,150,143, 50, &
      47,147,152,151,144, 51, 48,140,146,145,139, 52,  8, 56, 55, 54, 53,  7]
  !
  integer, parameter :: p6 = 343
  integer, parameter :: p7 = 512
  integer, parameter :: p8 = 729
  integer, parameter :: p9 = 1000
  !
continue
  !
  select case (size(pts))
    case (p1)
      return_value = pts( hexa1(:) )
    case (p2)
      return_value = pts( hexa2(:) )
    case (p3)
      return_value = pts( hexa3(:) )
    case (p4)
      return_value = pts( hexa4(:) )
    case (p5)
      return_value = pts( hexa5(:) )
   !case (p6)
   !  return_value = pts( hexa6(:) )
   !case (p7)
   !  return_value = pts( hexa7(:) )
   !case (p8)
   !  return_value = pts( hexa8(:) )
   !case (p9)
   !  return_value = pts( hexa9(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_gmsh_hexa_pts
!
!###############################################################################
!
pure function reorder_gmsh_edge_defined_tria_pts(pts) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: p1 = 3
  integer, parameter :: tria1(1:p1) = [1,2,3]
  !
  integer, parameter :: p2 = 6
  integer, parameter :: tria2(1:p2) = [1,4,2,6,5,3]
  !
  integer, parameter :: p3 = 9
  integer, parameter :: tria3(1:p3) = [1,4,5,2,9,6,8,7,3]
  !
  integer, parameter :: p4 = 12
  integer, parameter :: tria4(1:p4) = [ 1, 4, 5, 6, &
                                        2,12, 7,11, &
                                        8,10, 9, 3]
  !
  integer, parameter :: p5 = 15
  integer, parameter :: tria5(1:p5) = [ 1, 4, 5, 6, 7, &
                                        2,15, 8,14, 9, &
                                       13,10,12,11, 3]
  !
  integer, parameter :: p6 = 18
  integer, parameter :: tria6(1:p6) = [ 1, 4, 5, 6, 7, 8, &
                                        2,18, 9,17,10,16, &
                                       11,15,12,14,13, 3]
  !
  integer, parameter :: p7 = 21
  integer, parameter :: tria7(1:p7) = [ 1, 4, 5, 6, 7, 8, 9, &
                                        2,21,10,20,11,19,12, &
                                       18,13,17,14,16,15, 3]
  !
  integer, parameter :: p8 = 24
  integer, parameter :: tria8(1:p8) = [ 1, 4, 5, 6, 7, 8, 9,10, &
                                        2,24,11,23,12,22,13,21, &
                                       14,20,15,19,16,18,17, 3]
  !
  integer, parameter :: p9 = 27
  integer, parameter :: tria9(1:p9) = [ 1, 4, 5, 6, 7, 8, 9,10,11, &
                                        2,27,12,26,13,25,14,24,15, &
                                       23,16,22,17,21,18,20,19, 3]
  !
  integer, parameter :: p10 = 30
  integer, parameter :: tria10(1:p10) = [ 1, 4, 5, 6, 7, 8, 9,10,11,12, &
                                          2,30,13,29,14,28,15,27,16,26, &
                                         17,25,18,24,19,23,20,22,21, 3]
  !
continue
  !
  select case (size(pts))
    case (p1)
      return_value = pts( tria1(:) )
    case (p2)
      return_value = pts( tria2(:) )
    case (p3)
      return_value = pts( tria3(:) )
    case (p4)
      return_value = pts( tria4(:) )
    case (p5)
      return_value = pts( tria5(:) )
    case (p6)
      return_value = pts( tria6(:) )
    case (p7)
      return_value = pts( tria7(:) )
    case (p8)
      return_value = pts( tria8(:) )
    case (p9)
      return_value = pts( tria9(:) )
    case (p10)
      return_value = pts( tria10(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_gmsh_edge_defined_tria_pts
!
!###############################################################################
!
pure function reorder_gmsh_edge_defined_quad_pts(pts) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: p1 = 4
  integer, parameter :: quad1(1:p1) = [1,2,4,3]
  !
  integer, parameter :: p2 = 8
  integer, parameter :: quad2(1:p2) = [1,5,2,8,6,4,7,3]
  !
  integer, parameter :: p3 = 12
  integer, parameter :: quad3(1:p3) = [ 1, 5, 6, &
                                        2,12, 7, &
                                       11, 8, 4, &
                                       10, 9, 3]
  !
  integer, parameter :: p4 = 16
  integer, parameter :: quad4(1:p4) = [ 1, 5, 6, 7, &
                                        2,16, 8,15, &
                                        9,14,10, 4, &
                                       13,12,11, 3]
  !
  integer, parameter :: p5 = 20
  integer, parameter :: quad5(1:p5) = [ 1, 5, 6, 7, 8, &
                                        2,20, 9,19,10, &
                                       18,11,17,12, 4, &
                                       16,15,14,13, 3]
  !
  integer, parameter :: p6 = 24
  integer, parameter :: quad6(1:p6) = [ 1, 5, 6, 7, 8, 9, &
                                        2,24,10,23,11,22, &
                                       12,21,13,20,14, 4, &
                                       19,18,17,16,15, 3]
  !
  integer, parameter :: p7 = 28
  integer, parameter :: quad7(1:p7) = [ 1, 5, 6, 7, 8, 9,10, &
                                        2,28,11,27,12,26,13, &
                                       25,14,24,15,23,16, 4, &
                                       22,21,20,19,18,17, 3]
  !
  integer, parameter :: p8 = 32
  integer, parameter :: quad8(1:p8) = [ 1, 5, 6, 7, 8, 9,10,11, &
                                        2,32,12,31,13,30,14,29, &
                                       15,28,16,27,17,26,18, 4, &
                                       25,24,23,22,21,20,19, 3]
  !
  integer, parameter :: p9 = 36
  integer, parameter :: quad9(1:p9) = [ 1, 5, 6, 7, 8, 9,10,11,12, &
                                        2,36,13,35,14,34,15,33,16, &
                                       32,17,31,18,30,19,29,20, 4, &
                                       28,27,26,25,24,23,22,21, 3]
  !
  integer, parameter :: p10 = 40
  integer, parameter :: quad10(1:p10) = [ 1, 5, 6, 7, 8, 9,10,11,12,13, &
                                          2,40,14,39,15,38,16,37,17,36, &
                                         18,35,19,34,20,33,21,32,22, 4, &
                                         31,30,29,28,27,26,25,24,23, 3]
  !
continue
  !
  select case (size(pts))
    case (p1)
      return_value = pts( quad1(:) )
    case (p2)
      return_value = pts( quad2(:) )
    case (p3)
      return_value = pts( quad3(:) )
    case (p4)
      return_value = pts( quad4(:) )
    case (p5)
      return_value = pts( quad5(:) )
    case (p6)
      return_value = pts( quad6(:) )
    case (p7)
      return_value = pts( quad7(:) )
    case (p8)
      return_value = pts( quad8(:) )
    case (p9)
      return_value = pts( quad9(:) )
    case (p10)
      return_value = pts( quad10(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_gmsh_edge_defined_quad_pts
!
!###############################################################################
!
pure function reorder_gmsh_edge_defined_tetr_pts(pts) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: p1 = 4
  integer, parameter :: tetr1(1:p1) = [1,4,2,3]
  !
  integer, parameter :: p2 = 10
  integer, parameter :: tetr2(1:p2) = [1,8,4,5,10,2,7,9,6,3]
  !
  integer, parameter :: p3 = 16
  integer, parameter :: tetr3(1:p3) = [ 1,12,11, 4, 5,15, 6,16, &
                                        2,10,13, 7, 9,14, 8, 3]
  !
  integer, parameter :: p4 = 22
  integer, parameter :: tetr4(1:p4) = [ 1,16,15,14, 4, 5,20, 6,21, 7,22, &
                                        2,13,17, 8,12,18, 9,11,19,10, 3]
  !
  integer, parameter :: p5 = 28
  integer, parameter :: tetr5(1:p5) = [ 1,20,19,18,17, 4, 5, &
                                       25, 6,26, 7,27, 8,28, &
                                        2,16,21, 9,15,22,10, &
                                       14,23,11,13,24,12, 3]
  !
  integer, parameter :: p6 = 34
  integer, parameter :: tetr6(1:p6) = [ 1,24,23,22,21,20, 4, 5,30, 6,31, 7, &
                                       32, 8,33, 9,34, 2,19,25,10,18,26,11, &
                                       17,27,12,16,28,13,15,29,14, 3]
  !
  integer, parameter :: p7 = 40
  integer, parameter :: tetr7(1:p7) = [ 1,28,27,26,25,24,23, 4, 5,35, &
                                        6,36, 7,37, 8,38, 9,39,10,40, &
                                        2,22,29,11,21,30,12,20,31,13, &
                                       19,32,14,18,33,15,17,34,16, 3]
  !
  integer, parameter :: p8 = 46
  integer, parameter :: tetr8(1:p8) = [ 1,32,31,30,29,28,27,26, 4, 5,40, 6,41, &
                                        7,42, 8,43, 9,44,10,45,11,46, 2,25,33, &
                                       12,24,34,13,23,35,14,22,36,15,21,37,16, &
                                       20,38,17,19,39,18, 3]
  !
  integer, parameter :: p9 = 52
  integer, parameter :: tetr9(1:p9) = [ 1,36,35,34,33,32,31,30,29, 4, 5,45, 6, &
                                       46, 7,47, 8,48, 9,49,10,50,11,51,12,52, &
                                        2,28,37,13,27,38,14,26,39,15,25,40,16, &
                                       24,41,17,23,42,18,22,43,19,21,44,20, 3]
  !
  integer, parameter :: p10 = 58
  integer, parameter :: tetr10(1:p10) = [ 1,40,39,38,37,36,35,34,33,32, 4, 5, &
                                         50, 6,51, 7,52, 8,53, 9,54,10,55,11, &
                                         56,12,57,13,58, 2,31,41,14,30,42,15, &
                                         29,43,16,28,44,17,27,45,18,26,46,19, &
                                         25,47,20,24,48,21,23,49,22, 3]
  !
continue
  !
  select case (size(pts))
    case (p1)
      return_value = pts( tetr1(:) )
    case (p2)
      return_value = pts( tetr2(:) )
    case (p3)
      return_value = pts( tetr3(:) )
    case (p4)
      return_value = pts( tetr4(:) )
    case (p5)
      return_value = pts( tetr5(:) )
    case (p6)
      return_value = pts( tetr6(:) )
    case (p7)
      return_value = pts( tetr7(:) )
    case (p8)
      return_value = pts( tetr8(:) )
    case (p9)
      return_value = pts( tetr9(:) )
    case (p10)
      return_value = pts( tetr10(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_gmsh_edge_defined_tetr_pts
!
!###############################################################################
!
pure function reorder_gmsh_edge_defined_pyra_pts(pts) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: p1 = 5
  integer, parameter :: pyra1(1:p1) = [1,2,4,3,5]
  !
  integer, parameter :: p2 = 14
  integer, parameter :: pyra2(1:p2) = [1,14,6,2,7,9,4,11,3,8,10,13,12,5]
  !
  integer, parameter :: p3 = 21
  integer, parameter :: p4 = 29
  integer, parameter :: p5 = 37
  integer, parameter :: p6 = 45
  integer, parameter :: p7 = 53
  integer, parameter :: p8 = 61
  integer, parameter :: p9 = 69
  !
continue
  !
  select case (size(pts))
    case (p1)
      return_value = pts( pyra1(:) )
    case (p2)
      return_value = pts( pyra2(:) )
   !case (p3)
   !  return_value = pts( pyra3(:) )
   !case (p4)
   !  return_value = pts( pyra4(:) )
   !case (p5)
   !  return_value = pts( pyra5(:) )
   !case (p6)
   !  return_value = pts( pyra6(:) )
   !case (p7)
   !  return_value = pts( pyra7(:) )
   !case (p8)
   !  return_value = pts( pyra8(:) )
   !case (p9)
   !  return_value = pts( pyra9(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_gmsh_edge_defined_pyra_pts
!
!###############################################################################
!
pure function reorder_gmsh_edge_defined_pris_pts(pts) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: p1 = 6
  integer, parameter :: pris1(1:p1) = [1,2,3,4,5,6]
  !
  integer, parameter :: p2 = 15
  integer, parameter :: pris2(1:p2) = [1,7,2,8,10,3,9,11,12,4,13,5,14,15,6]
  !
  integer, parameter :: p3 = 24
  integer, parameter :: pris3(1:p3) = [ 1, 7, 8, 2, 9,13,10,14, 3,11,15,17, &
                                       12,16,18, 4,19,20, 5,21,23,22,24, 6]
  !
  integer, parameter :: p4 = 33
  integer, parameter :: pris4(1:p4) = [ 1, 7, 8, 9, 2,10,16,11,17,12,18, &
                                        3,13,19,22,14,20,23,15,21,24, 4, &
                                       25,26,27, 5,28,31,29,32,30,33, 6]
  !
  integer, parameter :: p5 = 42
  integer, parameter :: pris5(1:p5) = [ 1, 7, 8, 9,10, 2,11,19,12,20,13,21,14, &
                                       22, 3,15,23,27,16,24,28,17,25,29,18,26, &
                                       30, 4,31,32,33,34, 5,35,39,36,40,37,41, &
                                       38,42, 6]
  !
  integer, parameter :: p6 = 51
  integer, parameter :: pris6(1:p6) = [ 1, 7, 8, 9,10,11, 2,12,22,13,23,14,24, &
                                       15,25,16,26, 3,17,27,32,18,28,33,19,29, &
                                       34,20,30,35,21,31,36, 4,37,38,39,40,41, &
                                        5,42,47,43,48,44,49,45,50,46,51, 6]
  !
  integer, parameter :: p7 = 60
  integer, parameter :: pris7(1:p7) = [ 1, 7, 8, 9,10,11,12, 2,13,25,14,26, &
                                       15,27,16,28,17,29,18,30, 3,19,31,37, &
                                       20,32,38,21,33,39,22,34,40,23,35,41, &
                                       24,36,42, 4,43,44,45,46,47,48, 5,49, &
                                       55,50,56,51,57,52,58,53,59,54,60, 6]
  !
  integer, parameter :: p8 = 69
  integer, parameter :: pris8(1:p8) = [ 1, 7, 8, 9,10,11,12,13, 2,14,28,15,29, &
                                       16,30,17,31,18,32,19,33,20,34, 3,21,35, &
                                       42,22,36,43,23,37,44,24,38,45,25,39,46, &
                                       26,40,47,27,41,48, 4,49,50,51,52,53,54, &
                                       55, 5,56,63,57,64,58,65,59,66,60,67,61, &
                                       68,62,69, 6]
  !
  integer, parameter :: p9 = 78
  integer, parameter :: pris9(1:p9) = [ 1, 7, 8, 9,10,11,12,13,14, 2,15,31,16, &
                                       32,17,33,18,34,19,35,20,36,21,37,22,38, &
                                        3,23,39,47,24,40,48,25,41,49,26,42,50, &
                                       27,43,51,28,44,52,29,45,53,30,46,54, 4, &
                                       55,56,57,58,59,60,61,62, 5,63,71,64,72, &
                                       65,73,66,74,67,75,68,76,69,77,70,78, 6]
  !
continue
  !
  select case (size(pts))
    case (p1)
      return_value = pts( pris1(:) )
    case (p2)
      return_value = pts( pris2(:) )
    case (p3)
      return_value = pts( pris3(:) )
    case (p4)
      return_value = pts( pris4(:) )
    case (p5)
      return_value = pts( pris5(:) )
    case (p6)
      return_value = pts( pris6(:) )
    case (p7)
      return_value = pts( pris7(:) )
    case (p8)
      return_value = pts( pris8(:) )
    case (p9)
      return_value = pts( pris9(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_gmsh_edge_defined_pris_pts
!
!###############################################################################
!
pure function reorder_gmsh_edge_defined_hexa_pts(pts) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: p1 = 8
  integer, parameter :: hexa1(1:p1) = [1,2,4,3,5,6,8,7]
  !
  integer, parameter :: p2 = 20
  integer, parameter :: hexa2(1:p2) = [ 1, 9, 2,10,12, 4,14, 3,11,13, &
                                       16,15, 5,17, 6,18,19, 8,20, 7]
  !
  integer, parameter :: p3 = 32
  integer, parameter :: hexa3(1:p3) = [ 1, 9,10, 2,11,15,12,16, &
                                        4,20,19, 3,13,17,23,21, &
                                       14,18,24,22, 5,25,26, 6, &
                                       27,29,28,30, 8,32,31, 7]
  !
  integer, parameter :: p4 = 44
  integer, parameter :: hexa4(1:p4) = [ 1, 9,10,11, 2,12,18,13,19,14,20, &
                                        4,26,25,24, 3,15,21,30,27,16,22, &
                                       31,28,17,23,32,29, 5,33,34,35, 6, &
                                       36,39,37,40,38,41, 8,44,43,42, 7]
  !
  integer, parameter :: p5 = 56
  integer, parameter :: hexa5(1:p5) = [ 1, 9,10,11,12, 2,13,21, &
                                       14,22,15,23,16,24, 4,32, &
                                       31,30,29, 3,17,25,37,33, &
                                       18,26,38,34,19,27,39,35, &
                                       20,28,40,36, 5,41,42,43, &
                                       44, 6,45,49,46,50,47,51, &
                                       48,52, 8,56,55,54,53, 7]
  !
  integer, parameter :: p6 = 68
  integer, parameter :: hexa6(1:p6) = [ 1, 9,10,11,12,13, 2,14,24,15, &
                                       25,16,26,17,27,18,28, 4,38,37, &
                                       36,35,34, 3,19,29,44,39,20,30, &
                                       45,40,21,31,46,41,22,32,47,42, &
                                       23,33,48,43, 5,49,50,51,52,53, &
                                        6,54,59,55,60,56,61,57,62,58, &
                                       63, 8,68,67,66,65,64, 7]
  !
  integer, parameter :: p7 = 80
  integer, parameter :: hexa7(1:p7) = [ 1, 9,10,11,12,13,14, 2,15,27, &
                                       16,28,17,29,18,30,19,31,20,32, &
                                        4,44,43,42,41,40,39, 3,21,33, &
                                       51,45,22,34,52,46,23,35,53,47, &
                                       24,36,54,48,25,37,55,49,26,38, &
                                       56,50, 5,57,58,59,60,61,62, 6, &
                                       63,69,64,70,65,71,66,72,67,73, &
                                       68,74, 8,80,79,78,77,76,75, 7]
  !
  integer, parameter :: p8 = 92
  integer, parameter :: hexa8(1:p8) = [ 1, 9,10,11,12,13,14,15, 2,16,30,17, &
                                       31,18,32,19,33,20,34,21,35,22,36, 4, &
                                       50,49,48,47,46,45,44, 3,23,37,58,51, &
                                       24,38,59,52,25,39,60,53,26,40,61,54, &
                                       27,41,62,55,28,42,63,56,29,43,64,57, &
                                        5,65,66,67,68,69,70,71, 6,72,79,73, &
                                       80,74,81,75,82,76,83,77,84,78,85, 8, &
                                       92,91,90,89,88,87,86, 7]
  !
  integer, parameter :: p9 = 104
  integer, parameter :: hexa9(1:p9) = [  1, 9,10,11,12,13,14, 15, 16,  2, 17, &
                                        33,18,34,19,35,20,36, 21, 37, 22, 38, &
                                        23,39,24,40, 4,56,55, 54, 53, 52, 51, &
                                        50,49, 3,25,41,65,57, 26, 42, 66, 58, &
                                        27,43,67,59,28,44,68, 60, 29, 45, 69, &
                                        61,30,46,70,62,31,47, 71, 63, 32, 48, &
                                        72,64, 5,73,74,75,76, 77, 78, 79, 80, &
                                         6,81,89,82,90,83,91, 84, 92, 85, 93, &
                                        86,94,87,95,88,96, 8,104,103,102,101, &
                                       100,99,98,97, 7]
  !
continue
  !
  select case (size(pts))
    case (p1)
      return_value = pts( hexa1(:) )
    case (p2)
      return_value = pts( hexa2(:) )
    case (p3)
      return_value = pts( hexa3(:) )
    case (p4)
      return_value = pts( hexa4(:) )
    case (p5)
      return_value = pts( hexa5(:) )
    case (p6)
      return_value = pts( hexa6(:) )
    case (p7)
      return_value = pts( hexa7(:) )
    case (p8)
      return_value = pts( hexa8(:) )
    case (p9)
      return_value = pts( hexa9(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_gmsh_edge_defined_hexa_pts
!
!###############################################################################
!
subroutine gmsh_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
  !
  !.. Formal Arguments ..
  integer, optional, intent(in) :: iunit
  !
  !.. Local Scalars ..
  integer :: iou,i,j
  character(len=36) :: array_name
  type(memory) :: array,dt
  type(memory) :: total
  !
#ifndef DISABLE_DTIO
continue
  !
  iou = iout
  if (present(iunit)) iou = iunit
  !
  call total%reset
  !
  write (iou,1)
  !
  ! Check the memory of the arrays within this module
  !
  if (allocated(phys_grps)) then
    !
    call dt%set( storage_size(phys_grps,kind=inttype) , &
                         size(phys_grps,kind=inttype) )
    !
    do i = lbound(phys_grps,dim=1),ubound(phys_grps,dim=1)
      !
      call dt%add( storage_size(phys_grps(i)%phys_dim,kind=inttype) )
      call dt%add( storage_size(phys_grps(i)%phys_id,kind=inttype) )
      call dt%add( storage_size(phys_grps(i)%ibc,kind=inttype) )
      call dt%add( storage_size(phys_grps(i)%phys_nam,kind=inttype) )
      !
      call array%set( storage_size(phys_grps(i)%cells,kind=inttype) , &
                              size(phys_grps(i)%cells,kind=inttype) )
      !
      write (array_name,4) "phys_grps",i,"cells"
      write (iou,2) array_name,array%get_units()
      !
      call dt%add(array)
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "phys_grps"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "phys_grps"
    !
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module module_gmsh")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  4 format (a,"(",i0,")%",a)
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine gmsh_memory_usage
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
