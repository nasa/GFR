module module_cgns
  !
  !.. Use Statements
  use module_kind_types
  use module_cgns_types
  use geovar, only : grid
  use geovar, only : elem_t
  use geovar, only : prop_t
 !use geovar, only : cgns_idx_t
  !
  implicit none
  !
  private
  !
  ! ##################################
  ! ###  Module Public Procedures  ###
  ! ##################################
  !
  ! Public Interface for CGNS element_properties function
  !
  interface cgns_element_properties
    module procedure element_properties
  end interface cgns_element_properties
  !
  public :: read_cgns_gridfile
  public :: cgns_element_properties
  public :: cgns_memory_usage
  !
  ! ###########################################################################
  ! ###########################################################################
  ! ###                                                                     ###
  ! ### Explanation of different element types defined by the CGNS standard ###
  ! ###                                                                     ###
  ! ###########################################################################
  ! ###########################################################################
  !
  ! ----------------------------------------------------------------------------
  !
  ! NOTES: * The terms EDGE-DEFINED and SERENDIPITY are synonymous.
  !        * Gmsh exclusively uses the terms incomplete and serendipity
  !          interchangeably but the CGNS standard uses none of these terms.
  !        * I coined the terms EDGE-DEFINED and FACE-DEFINED for CGNS to
  !          differentiate between these two element types because they both are
  !          incomplete element types.
  !        * EDGE-DEFINED elements are equivalent to Gmsh SERENDIPITY elements
  !          and only define the high-order nodes on the edges of the element.
  !        * FACE-DEFINED elements go one step further and also define the
  !          high-order interior nodes on the faces of a volume element.
  !        * A COMPLETE element takes the final step and defines the high-order
  !          nodes within the interior of the area or volume element.
  !        * EDGE-DEFINED/SERENDIPITY elements are possible for both 2D and 3D
  !          element types.
  !        * FACE-DEFINED elements are only possible for 3D element types.
  !
  ! ----------------------------------------------------------------------------
  !
  ! NODE : 1-node 1st order node
  !
  ! BAR_2 : 2-node 1st order edge
  ! BAR_3 : 3-node 2nd order edge
  ! BAR_4 : 4-node 3rd order edge
  ! BAR_5 : 5-node 4th order edge
  !
  ! TRI_3  :  3-node 1st order     COMPLETE triangle
  ! TRI_6  :  6-node 2nd order     COMPLETE triangle
  ! TRI_9  :  9-node 3rd order EDGE-DEFINED triangle
  ! TRI_10 : 10-node 3rd order     COMPLETE triangle
  ! TRI_12 : 12-node 4th order EDGE-DEFINED triangle
  ! TRI_15 : 15-node 4th order     COMPLETE triangle
  !
  ! QUAD_4     :  4-node 1st order     COMPLETE quadrilateral
  ! QUAD_8     :  8-node 2nd order EDGE-DEFINED quadrilateral
  ! QUAD_9     :  9-node 2nd order     COMPLETE quadrilateral
  ! QUAD_12    : 12-node 3rd order EDGE-DEFINED quadrilateral
  ! QUAD_16    : 16-node 3rd order     COMPLETE quadrilateral
  ! QUAD_P4_16 : 16-node 4th order EDGE-DEFINED quadrilateral
  ! QUAD_25    : 25-node 4th order     COMPLETE quadrilateral
  !
  ! TETRA_4  :  4-node 1st order     COMPLETE tetrahedron
  ! TETRA_10 : 10-node 2nd order     COMPLETE tetrahedron
  ! TETRA_16 : 16-node 3rd order EDGE-DEFINED tetrahedron
  ! TETRA_20 : 20-node 3rd order     COMPLETE tetrahedron
  ! TETRA_22 : 22-node 4th order EDGE-DEFINED tetrahedron
  ! TETRA_34 : 34-node 4th order FACE-DEFINED tetrahedron
  ! TETRA_35 : 35-node 4th order     COMPLETE tetrahedron
  !
  ! PYRA_5     :  5-node 1st order     COMPLETE pyramid
  ! PYRA_13    : 13-node 2nd order EDGE-DEFINED pyramid
  ! PYRA_14    : 14-node 2nd order     COMPLETE pyramid
  ! PYRA_21    : 21-node 3rd order EDGE-DEFINED pyramid
  ! PYRA_29    : 29-node 3rd order FACE-DEFINED pyramid
  ! PYRA_30    : 30-node 3rd order     COMPLETE pyramid
  ! PYRA_P4_29 : 29-node 4th order EDGE-DEFINED pyramid
  ! PYRA_50    : 50-node 4th order FACE-DEFINED pyramid
  ! PYRA_55    : 55-node 4th order     COMPLETE pyramid
  !
  ! PENTA_6  :  6-node 1st order     COMPLETE prism
  ! PENTA_15 : 15-node 2nd order EDGE-DEFINED prism
  ! PENTA_18 : 18-node 2nd order     COMPLETE prism
  ! PENTA_24 : 24-node 3rd order EDGE-DEFINED prism
  ! PENTA_38 : 38-node 3rd order FACE-DEFINED prism
  ! PENTA_40 : 40-node 3rd order     COMPLETE prism
  ! PENTA_33 : 33-node 4th order EDGE-DEFINED prism
  ! PENTA_66 : 66-node 4th order FACE-DEFINED prism
  ! PENTA_75 : 75-node 4th order     COMPLETE prism
  !
  ! HEXA_8   :   8-node 1st order     COMPLETE hexahedron
  ! HEXA_20  :  20-node 2nd order EDGE-DEFINED hexahedron
  ! HEXA_27  :  27-node 2nd order     COMPLETE hexahedron
  ! HEXA_32  :  32-node 3rd order EDGE-DEFINED hexahedron
  ! HEXA_56  :  56-node 3rd order FACE-DEFINED hexahedron
  ! HEXA_64  :  64-node 3rd order     COMPLETE hexahedron
  ! HEXA_44  :  44-node 4th order EDGE-DEFINED hexahedron
  ! HEXA_98  :  98-node 4th order FACE-DEFINED hexahedron
  ! HEXA_125 : 125-node 4th order     COMPLETE hexahedron
  !
  ! MIXED : Arbitrary mixture of all possible element types 
  !         except NGON_n or NFACE_n
  !
  ! ################################################################
  ! #####   NON-COMPATIBLE ELEMENT TYPES, THEREFORE NOT USED   #####
  ! ################################################################
  !
  ! NGON_n     : Face elements used for defining arbitrary polyhedral elements
  ! NFACE_n    : Arbitrary polyhedral element defined by collection of
  !              face elements of NGON_n element type
  !
  !
  ! ###########################
  ! ###  Module Parameters  ###
  ! ###########################
  !
  logical(lk), parameter :: debug_cgns = fals
 !logical(lk), parameter :: debug_cgns = true
  !
  ! ###################################
  ! ###  Module Derived Data Types  ###
  ! ###################################
  !
  type :: cell_grps_t
    integer :: cell_dim
    integer :: grp_id
    integer :: ibc
    character(len=CGLEN) :: name
    integer, allocatable :: cells(:)
  end type cell_grps_t
  !
  type :: cg_t
    character(len=CGLEN) :: name
  end type cg_t
  !
  type, extends(cg_t) :: cg_sz_t
    integer(CST) :: nsize
  end type cg_sz_t
  !
  type, extends(cg_t) :: cg_dataarray_t
    real(wp), allocatable :: data(:)
  end type cg_dataarray_t
  !
  type, extends(cg_t) :: cg_bc_t
    integer(CET) :: boco_type, ptset_type, normal_datatype, location
    integer(CST) :: npnts, normal_list_size
    integer(CBT) :: normal_index(3), ndataset
    integer(CST), allocatable :: pnts(:)!, faces(:)
    real(wp),     allocatable :: normal_list(:)
  end type cg_bc_t
  !
  type, extends(cg_t) :: cg_coord_t
    real(wp), allocatable :: x(:)
  end type cg_coord_t
  !
  type, extends(cg_t) :: cg_grid_t
    type(cg_coord_t), allocatable :: coord(:)
  end type cg_grid_t
  !
  type, extends(cg_t) :: cg_section_t
    integer(CET) :: type
    integer(CST) :: start, end
    integer(CBT) :: nbndry, parent_flag
    integer(CST), allocatable :: elements(:), parentdata(:)
   !integer(CET), allocatable :: elemtypes(:)
  end type cg_section_t
  !
  type, extends(cg_sz_t) :: cg_zone_t
    integer(CST) :: nvertex, ncell, nboundvertex
    integer(CET) :: type
    integer(CBT) :: index_dim
    type(cg_grid_t),    allocatable :: grid(:)
    type(cg_section_t), allocatable :: section(:)
    type(cg_bc_t),      allocatable :: bc(:)
  end type cg_zone_t
  !
  type, extends(cg_sz_t) :: cg_base_t
    integer(CBT) :: cell_dim, phys_dim
    type(cg_zone_t), allocatable :: zone(:)
  end type cg_base_t
  !
  ! ################################
  ! ###  Module Local Variables  ###
  ! ################################
  !
  integer(CBT), save :: if_g  ! CGNS file index number
  integer(CBT), save :: ib_g  ! base index number
  integer(CBT), save :: iz_g  ! zone index number
  integer(CBT), save :: ig_g  ! grid index number
  integer(CBT), save :: is_g  ! element section index number
  integer(CBT), save :: ic_g  ! coordinate array index number
  integer(CBT), save :: ia_g  ! data array index number
  integer(CBT), save :: ibc_g ! boundary condition index number
  !
  character(len=150),  save :: gridfile
  character(len=1000), save :: array_name
  !
  ! Grid dimensions specified by the CGNS grid file
  !
  integer, save :: npnts_correct ! total # of grid points
  integer, save :: ncell_correct ! total # of interior cells
  integer, save :: nfbnd_correct ! total # of boundary faces
  !
  ! ##################################
  ! ###  Module Fixed-Size Arrays  ###
  ! ##################################
  !
  !    data_sizes( 1) = base(1)%cell_dim
  !    data_sizes( 2) = base(1)%phys_dim
  !    data_sizes( 3) = base(1)%zone(1)%nvertex
  !    data_sizes( 4) = base(1)%zone(1)%ncell
  !    data_sizes( 5) = base(1)%zone(1)%nboundvertex
  !    data_sizes( 6) = total # of grids across all zones
  !    data_sizes( 7) = size( grid_data )
  !    data_sizes( 8) = total # of element sections across all zones
  !    data_sizes( 9) = size( section_data )
  !    data_sizes(10) = total # of BCs across all zones
  !    data_sizes(11) = size( boco_data )
  !
  integer(wip), save :: data_sizes(1:11) = 0
  !
  ! ###################################
  ! ###  Module Allocatable Arrays  ###
  ! ###################################
  !
  type(cg_base_t), save, allocatable :: base(:)
  !
  real(wp),     save, allocatable :: grid_data(:,:)
  integer(wip), save, allocatable :: section_data(:)
  integer(wip), save, allocatable :: boco_data(:)
  !
  type(cell_grps_t), save, allocatable :: cell_grps(:)
  !
  ! ######################################
  ! ###  Generic Procedure Interfaces  ###
  ! ######################################
  !
  interface create_element
    module procedure create_element_i4, create_element_i8
  end interface create_element
  !
contains
!
!###############################################################################
!
subroutine read_cgns_gridfile( input_gridfile )
  !
  !.. Use Statements ..
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: input_gridfile
  !
  !.. Local Scalars ..
  integer      :: ierr
  integer(CBT) :: nbases,grid_precision
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_cgns_gridfile"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  gridfile = trim(adjustl(input_gridfile))
  !
  call cgp_mpi_comm_f(mpi_comm_self_val,cgierr)
 !call cgns_error(pname,gridfile,"cgp_mpi_comm_f",cgierr,__LINE__,__FILE__)
  !
  ! Open the grid file
  !
  call cgp_open_f(trim(adjustl(gridfile)),CG_MODE_READ,if_g,cgierr)
  if (debug_cgns) write (*,1) mypnum,"opened CGNS grid file"
  call cgns_error(pname,gridfile,"cgp_open_f",cgierr, &
                  __LINE__,__FILE__,"Mode = CG_MODE_READ")
  !
  ! Have only the root processor read in the data from the CGNS grid file
  !
  if (mypnum == 0) then
    !
    ! Get the precision of the grid file
    !
    call cg_precision_f(if_g,grid_precision,cgierr)
    call cgns_error(pname,gridfile,"cg_precision_f",cgierr, &
                    __LINE__,__FILE__,not_parallel=true)
    if (debug_cgns) then
      write (iout,1) mypnum, &
                     "The CGNS grid file was written using a precision of ", &
                     grid_precision
    end if
    !
    ! Read the number of bases
    !
    call cg_nbases_f(if_g,nbases,cgierr)
    call cgns_error(pname,gridfile,"cg_nbases_f",cgierr, &
                    __LINE__,__FILE__,not_parallel=true)
    !
    ! Allocate the base array
    !
    allocate ( base(1:nbases) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"base",1,__LINE__,__FILE__,ierr,error_message)
    !
    do ib_g = 1,size(base)
      if (debug_cgns) write (iout,1) mypnum,"reading CGNS base node #",ib_g
      call read_cgns_base( base(ib_g) )
    end do
    !
  end if
  !
  ! We are done reading the grid file so close it
  !
  call cgp_close_f(if_g,cgierr)
  if (debug_cgns) write (*,1) mypnum,"closed CGNS grid file"
  call cgns_error(pname,gridfile,"cgp_close_f",cgierr,__LINE__,__FILE__)
  !
  !
  !
  if (mypnum == 0) then
    !
    ! Check to see if there are any unsupported features in the CGNS grid file
    !
    call check_for_unsupported_cgns_file
    !
    ! Take the information read from the CGNS grid file and transfer it into
    ! the grid derived type and cell_grps arrays
    !
    call transfer_cgns_data_to_grid_dt
    !
    ! Create the arrays needed by the solver that define the grid
    ! geometry: bface, nodes_of_cell_ptr, nodes_of_cell, etc.
    !
    call create_solver_geom_arrays_cgns
    !
  end if
  !
  ! It should now be safe to deallocate the base derived-type array.
  ! We dont need it wasting a bunch of memory since we are about to
  ! give each MPI process copies of the grid_data, section_data, and
  ! boco_data arrays.
  !
  if (allocated(base)) then
    deallocate ( base , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"base",2,__LINE__,__FILE__,ierr,error_message)
  end if
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  ! Reset the MPI communicator for CGNS back to MPI_COMM_WORLD
  !
  call cgp_mpi_comm_f(mpi_comm_world_val,cgierr)
 !call cgns_error(pname,gridfile,"cgp_mpi_comm_f",cgierr,__LINE__,__FILE__)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  ! Broadcast the grid arrays to all processors
  !
  if (ncpu > 1) then
   !call broadcast_grid_arrays
    call broadcast_solver_geom_arrays
  else
    call create_serial_geom_arrays
  end if
  if (debug_cgns) then
    if (mypnum == 0) then
      write (*,*)
      write (*,*) "Leaving read_cgns_gridfile! Press any key to continue"
      read (*,*)
    end if
  end if
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
 !call stop_recon(stop_mpi,pname,__LINE__,__FILE__, &
 !                "TEMPORARY STOP AFTER READING CGNS GRID")
  !
  ! Finish create the arrays defining the grid geometry:
  ! bface, nodes_of_cell_ptr, nodes_of_cell, etc.
  !
 !call create_arrays_from_cgns_data
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (" CPU #",i0,": ",a,:,i0)
  !
end subroutine read_cgns_gridfile
!
!###############################################################################
!
subroutine check_for_unsupported_cgns_file
  !
  !.. Local Scalars ..
  integer :: ib,iz,ig,ic,is,n
  logical(lk) :: grid_requires_int64
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_for_unsupported_cgns_file"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Check whether 64-bit integers are required for this grid
  !
  call debug_timer(start_timer,"check_if_grid_requires_int64")
  grid_requires_int64 = check_if_grid_requires_int64()
  call debug_timer(stop_timer,"check_if_grid_requires_int64")
  !
  if (grid_requires_int64) then
    write (error_message,1)
    call stop_recon(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! FOR NOW: Abort if there is more than one base in the CGNS file
  !
  if (allocated(base)) then
    if (size(base) > 1) then
      write (error_message,2)
      call stop_recon(abort,pname,__LINE__,__FILE__,error_message)
    end if
  else
    write (error_message,12)
    call stop_recon(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! FOR NOW: Abort if there is more than one zone in the CGNS file for any base
  !
  do ib = 1,size(base)
    if (allocated(base(ib)%zone)) then
      if (size(base(ib)%zone) > 1) then
        write (error_message,3)
        call stop_recon(abort,pname,__LINE__,__FILE__,error_message)
      end if
    else
      write (error_message,13)
      call stop_recon(abort,pname,__LINE__,__FILE__,error_message)
    end if
  end do
  !
  ! Abort if there is a section that contains an
  ! unsupported element type (NGON_n or NFACE_n)
  !
  do ib = 1,size(base)
    do iz = 1,size(base(ib)%zone)
      if (allocated(base(ib)%zone(iz)%section)) then
        do is = 1,size(base(ib)%zone(iz)%section)
          if (any(base(ib)%zone(iz)%section(is)%type == [NGON_n,NFACE_n])) then
            write (error_message,4)
            call stop_recon(abort,pname,__LINE__,__FILE__,error_message)
          end if
        end do
      else
        write (error_message,14)
        call stop_recon(abort,pname,__LINE__,__FILE__,error_message)
      end if
    end do
  end do
  !
  ! FOR NOW: Abort if there is more than one set of grid coordinates
  !          for any zone
  !
  do ib = 1,size(base)
    do iz = 1,size(base(ib)%zone)
      if (allocated(base(ib)%zone(iz)%grid)) then
        if (size(base(ib)%zone(iz)%grid) > 1) then
          write (error_message,5)
          call stop_recon(abort,pname,__LINE__,__FILE__,error_message)
        end if
      else
        write (error_message,15)
        call stop_recon(abort,pname,__LINE__,__FILE__,error_message)
      end if
    end do
  end do
  !
  ! FOR NOW: Abort if the number of coordinate data is inconsistent within
  !          any of the grids
  !
  do ib = 1,size(base)
    do iz = 1,size(base(ib)%zone)
      do ig = 1,size(base(ib)%zone(iz)%grid)
        if (allocated(base(ib)%zone(iz)%grid(ig)%coord)) then
          n = size(base(ib)%zone(iz)%grid(ig)%coord(1)%x)
          do ic = 2,size(base(ib)%zone(iz)%grid(ig)%coord)
            if (n /= size(base(ib)%zone(iz)%grid(ig)%coord(ic)%x)) then
              write (error_message,6)
              call stop_recon(abort,pname,__LINE__,__FILE__,error_message)
            end if
          end do
        else
          write (error_message,16)
          call stop_recon(abort,pname,__LINE__,__FILE__,error_message)
        end if
      end do
    end do
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (" THE CGNS GRID FILE APPEARS TO REQUIRE 64-BIT INTEGERS", &
            " WHICH IS NOT CURRENTLY SUPPORTED IN THE SOLVER!!!")
  2 format (" THE CGNS GRID FILE APPEARS TO CONTAIN MORE THAN ONE CGNSBase_t", &
            " NODE WHICH IS NOT CURRENTLY SUPPORTED IN THE CGNS GRID READER!!!")
  3 format (" THE CGNS GRID FILE APPEARS TO CONTAIN MORE THAN ONE Zone_t", &
            " NODE WHICH IS NOT CURRENTLY SUPPORTED IN THE CGNS GRID READER!!!")
  4 format (" THE CGNS GRID FILE APPEARS TO CONTAIN AN Elements_t NODE", &
            " WITH AN UNSUPPORTED ELEMENT TYPE OF NGON_n or NFACE_n!!")
  5 format (" THE CGNS GRID FILE APPEARS TO CONTAIN MORE THAN ONE", &
            " GridCoordinates_t NODE WHICH IS NOT CURRENTLY SUPPORTED IN", &
            " THE CGNS GRID READER!!!")
  6 format (" THE CGNS GRID FILE APPEARS TO CONTAIN AN INCONSISTENCY WITHIN", &
            " ONE OF THE GridCoordinates_t NODES!!! THE SIZE OF ALL THE", &
            " DataArray_t NODES UNDER A GridCoordinate_t NODE SHOULD ALL", &
            " BE THE SAME SIZE!!!")
  !
  12 format(" THE base DERIVED-TYPE ARRAY USED TO STORE ALL THE DATA READ", &
            " IN FROM THE CGNS GRID FILE IS NOT ALLOCATED WHEN IT SHOULD BE!!!")
  13 format(" THE zone DERIVED-TYPE ARRAY FOR STORING THE Zone_t NODES UNDER", &
            " THE base(",i0,") CGNSBase_t NODE IS NOT ALLOCATED WHEN IT", &
            " SHOULD BE!!!")
  14 format(" THE section DERIVED-TYPE ARRAY FOR STORING THE Elements_t", &
            " NODES UNDER THE base(",i0,")%zone(",i0,") Zone_t NODE IS NOT", &
            " ALLOCATED WHEN IT SHOULD BE!!!")
  15 format(" THE grid DERIVED-TYPE ARRAY FOR STORING THE GridCoordinates_t", &
            " NODES UNDER THE base(",i0,")%zone(",i0,") Zone_t NODE IS NOT", &
            " ALLOCATED WHEN IT SHOULD BE!!!")
  16 format(" THE coord DERIVED-TYPE ARRAY FOR STORING THE DataArray_t", &
            " NODES UNDER THE base(",i0,")%zone(",i0,")%grid(",i0,")", &
            " GridCoordinates_t NODE IS NOT ALLOCATED WHEN IT SHOULD BE!!!")
  !
end subroutine check_for_unsupported_cgns_file
!
!###############################################################################
!
subroutine create_broadcast_arrays
  !
  !.. Local Scalars ..
  integer :: i,j,n,n1,n2,ne
  integer :: ib,iz,ig,is,ibc,ierr
  integer :: ngrids,nsections,nbocos
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_broadcast_arrays"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Create the initial array that will contain the sizes of 
  ! all the subsequent arrays that will be broadcast
  !
  data_sizes(1) = int( base(1)%cell_dim , kind=kind(data_sizes) )
  data_sizes(2) = int( base(1)%phys_dim , kind=kind(data_sizes) )
  data_sizes(3) = int( base(1)%zone(1)%nvertex , kind=kind(data_sizes) )
  data_sizes(4) = int( base(1)%zone(1)%ncell , kind=kind(data_sizes) )
  data_sizes(5) = int( base(1)%zone(1)%nboundvertex , kind=kind(data_sizes) )
  !
  ! Get the size needed for the grid_data array
  !
  ngrids = 0
  !
  n1 = data_sizes(2)
  n2 = 0
  !
  do ib = 1,size(base)
    do iz = 1,size(base(ib)%zone)
      do ig = 1,size(base(ib)%zone(iz)%grid)
        if (allocated(base(ib)%zone(iz)%grid(ig)%coord)) then
          ngrids = ngrids + 1
          n2 = n2 + 1 ! Add one to store the dimensions for this grid
          n2 = n2 + size(base(ib)%zone(iz)%grid(ig)%coord(1)%x)
        end if
      end do
    end do
  end do
  !
  ! Allocate the grid_data array
  !
  allocate ( grid_data(1:n1,1:n2) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid_data",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the grid_data array
  !
  n = 0
  do ib = 1,size(base)
    do iz = 1,size(base(ib)%zone)
      do ig = 1,size(base(ib)%zone(iz)%grid)
        if (allocated(base(ib)%zone(iz)%grid(ig)%coord)) then
          ! 
          n1 = size(base(ib)%zone(iz)%grid(ig)%coord)
          n2 = size(base(ib)%zone(iz)%grid(ig)%coord(1)%x)
          !
          n = n + 1
          !
          ! Transfer the bits of the integers n1 and n2 into (haha) a
          ! working-precision real and store in the grid_data array.
          ! The other processors will have to transfer these back into
          ! integers in order to get the coordinate array dimensions.
          ! NOTE: The constant zero is used just as a "mold" variable which
          !       means the integers will be transfered into a variable
          !       with the same kind and kind-type as the mold variable.
          !
          grid_data(1,n) = transfer( n1 , zero )
          grid_data(2,n) = transfer( n2 , zero )
          !
          do j = 1,n2
            n = n + 1
            do i = 1,n1
              grid_data(i,n) = base(ib)%zone(iz)%grid(ig)%coord(i)%x(j)
            end do
          end do
          !
        end if
      end do
    end do
  end do
  !
  ! Get the size needed for the section_data array
  !
  nsections = 0
  !
  n1 = 0
  !
  do ib = 1,size(base)
    do iz = 1,size(base(ib)%zone)
      do is = 1,size(base(ib)%zone(iz)%section)
        nsections = nsections + 1
        n1 = n1 + 5 ! account for type, start, end, and nbndry components,
        !         ! and size of elements component array
        n1 = n1 + size(base(ib)%zone(iz)%section(is)%elements)
        ! IGNORE PARENT DATA FOR NOW
      end do
    end do
  end do
  !
  ! Allocate the section_data array
  !
  allocate ( section_data(1:n1) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"section_data",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the section_data array
  !
  n = 0
  do ib = 1,size(base)
    do iz = 1,size(base(ib)%zone)
      do is = 1,size(base(ib)%zone(iz)%section)
        !
        ne = size(base(ib)%zone(iz)%section(is)%elements)
        !
        section_data(n+1) = base(ib)%zone(iz)%section(is)%type
        section_data(n+2) = int( base(ib)%zone(iz)%section(is)% &
                                 start , kind=kind(section_data) )
        section_data(n+3) = int( base(ib)%zone(iz)%section(is)% &
                                 end , kind=kind(section_data) )
        section_data(n+4) = base(ib)%zone(iz)%section(is)%nbndry
        section_data(n+5) = ne
        !
        n = n + 5 ! account for type, start, end, and nbndry components,
        !           and size of elements component array
        n1 = n + 1
        n2 = n + ne
        section_data(n1:n2) = int( base(ib)%zone(iz)%section(is)% &
                                   elements(1:ne) , kind=kind(section_data) )
        !
        n = n2
        ! IGNORE PARENT DATA FOR NOW
      end do
    end do
  end do
  !
  ! Get the size needed for the boco_data array
  !
  nbocos = 0
  !
  n1 = 0
  !
  do ib = 1,size(base)
    do iz = 1,size(base(ib)%zone)
      do ibc = 1,size(base(ib)%zone(iz)%bc)
        nbocos = nbocos + 1
        n1 = n1 + 3 ! account for boco_type, npnts, and location components
        n1 = n1 + size(base(ib)%zone(iz)%bc(ibc)%pnts)
        ! IGNORE NORMAL DATA FOR NOW
      end do
    end do
  end do
  !
  ! Allocate the boco_data array
  !
  allocate ( boco_data(1:n1) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"boco_data",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the boco_data array
  !
  n = 0
  !
  do ib = 1,size(base)
    do iz = 1,size(base(ib)%zone)
      do ibc = 1,size(base(ib)%zone(iz)%bc)
        !
        n1 = size( base(ib)%zone(iz)%bc(ibc)%pnts )
        !
        boco_data(n+1) = base(ib)%zone(iz)%bc(ibc)%boco_type
        boco_data(n+2) = base(ib)%zone(iz)%bc(ibc)%location
        boco_data(n+3) = n1
        !
        n = n + 3 ! account for boco_type, npnts, and location components
        boco_data(n+1:n+n1) = int( base(ib)%zone(iz)%bc(ibc)%pnts(1:n1) , &
                                   kind=kind(boco_data) )
        !
        n = n + n1
        ! IGNORE NORMAL DATA FOR NOW
      end do
    end do
  end do
  !
  ! Store the sizes of the three data arrays in data_sizes
  !
  data_sizes( 6) = ngrids
  data_sizes( 7) = size( grid_data )
  data_sizes( 8) = nsections
  data_sizes( 9) = size( section_data )
  data_sizes(10) = nbocos
  data_sizes(11) = size( boco_data )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine create_broadcast_arrays
!
!###############################################################################
!
subroutine broadcast_grid_arrays
  !
  !.. Local Scalars ..
  integer :: n1,n2,ierr
  integer(INT_MPI) :: isz
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "broadcast_grid_arrays"
  integer(INT_MPI), parameter :: root = 0_int_mpi
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Broadcast the data_sizes array
  !
 !write (iout,1) mypnum,wip,(data_sizes(n1),n1=1,size(data_sizes))
  call mpi_bcast(data_sizes,size(data_sizes,kind=int_mpi),mpi_inttyp, &
                 root,MPI_COMM_WORLD,mpierr)
 !write (iout,1) mypnum,wip,(data_sizes(n1),n1=1,size(data_sizes))
  !
  ! Have all the non-root processors allocate the
  ! grid_data, section_data, and boco_data arrays
  !
  if (mypnum /= root) then
    !
    ! Allocate the grid_data array
    !
    n1 = data_sizes(2)
    n2 = data_sizes(7) / n1
    allocate ( grid_data(1:n1,1:n2) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid_data",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Allocate the section_data array
    !
    n1 = data_sizes(9)
    allocate ( section_data(1:n1) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"section_data",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    ! Allocate the boco_data array
    !
    n1 = data_sizes(11)
    allocate ( boco_data(1:n1) , source=0 , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"boco_data",1,__LINE__,__FILE__,ierr,error_message)
    !
  end if
  !
  ! Broadcast the grid_data array
  !
  isz = int( data_sizes(7) , kind=int_mpi )
  call mpi_bcast(grid_data,isz,mpi_inttyp,root,MPI_COMM_WORLD,mpierr)
  !
  ! Broadcast the section_data array
  !
  isz = int( data_sizes(9) , kind=int_mpi )
  call mpi_bcast(section_data,isz,mpi_inttyp,root,MPI_COMM_WORLD,mpierr)
  !
  ! Broadcast the boco_data array
  !
  isz = int( data_sizes(11) , kind=int_mpi )
  call mpi_bcast(boco_data,isz,mpi_inttyp,root,MPI_COMM_WORLD,mpierr)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (" CPU #",i0,": data_sizes(:) = ",6(1x,i0),/,28x,*(6(1x,i0,:),/))
  !
end subroutine broadcast_grid_arrays
!
!###############################################################################
!
subroutine read_cgns_base(b)
  !
  !.. Formal Arguments ..
  type(cg_base_t), intent(inout) :: b
  !
  !.. Local Scalars ..
  integer      :: ierr
  integer(CBT) :: nzones
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_cgns_base"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Read the information for the current base
  !
  call cg_base_read_f(if_g,ib_g,b%name,b%cell_dim,b%phys_dim,cgierr)
  call cgns_error(pname,gridfile,"cg_base_read_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  ! Read the number of zones in the current base
  !
  call cg_nzones_f(if_g,ib_g,nzones,cgierr)
  call cgns_error(pname,gridfile,"cg_nzones_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  ! Allocate the zones component for the current base
  !
  allocate ( b%zone(1:nzones) , stat=ierr , errmsg=error_message )
  array_name = get_array_name()
  call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read each of the zones
  !
  do iz_g = 1,size(b%zone)
    if (debug_cgns) write (iout,1) mypnum,"reading CGNS zone node #",iz_g
    call read_cgns_zone( b%zone(iz_g) )
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (" CPU #",i0,": ",a,:,i0)
  !
end subroutine read_cgns_base
!
!###############################################################################
!
subroutine read_cgns_zone(z)
  !
  !.. Formal Arguments ..
  type(cg_zone_t), intent(inout) :: z
  !
  !.. Local Scalars ..
  integer      :: ierr
  integer(CBT) :: ngrids,nsections,nbocos
  !
  !.. Local Arrays ..
  integer(CST) :: zone_size(3*3)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_cgns_zone"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Read the information for the current zone
  !
  call cg_zone_read_f(if_g,ib_g,iz_g,z%name,zone_size,cgierr)
  call cgns_error(pname,gridfile,"cg_zone_read_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  z%nvertex = zone_size(1)
  z%ncell = zone_size(2)
  z%nboundvertex = zone_size(3)
  !
  call cg_zone_type_f(if_g,ib_g,iz_g,z%type,cgierr)
  call cgns_error(pname,gridfile,"cg_zone_type_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  call cg_index_dim_f(if_g,ib_g,iz_g,z%index_dim,cgierr)
  call cgns_error(pname,gridfile,"cg_index_dim_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  if (z%type /= UNSTRUCTURED .or. base(ib_g)%zone(iz_g)%index_dim /= 1) then
    write (error_message,1) ib_g,iz_g
    call stop_recon(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Read the number of grids that are in the current zone
  !
  call cg_ngrids_f(if_g,ib_g,iz_g,ngrids,cgierr)
  call cgns_error(pname,gridfile,"cg_ngrids_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  ! Allocate the grid component for the current zone
  !
  allocate ( z%grid(1:ngrids) , stat=ierr , errmsg=error_message )
  array_name = get_array_name("grid")
  call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read the number of element sections that are in the current zone
  !
  call cg_nsections_f(if_g,ib_g,iz_g,nsections,cgierr)
  call cgns_error(pname,gridfile,"cg_nsections_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  ! Allocate the section component for the current zone
  !
  allocate ( z%section(1:nsections) , stat=ierr , errmsg=error_message )
  array_name = get_array_name("section")
  call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read the number of boundary conditions that are in the current zone
  !
  call cg_nbocos_f(if_g,ib_g,iz_g,nbocos,cgierr)
  call cgns_error(pname,gridfile,"cg_nbocos_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  ! Allocate the bc component for the current zone
  !
  allocate ( z%bc(1:nbocos) , stat=ierr , errmsg=error_message )
  array_name = get_array_name("bc")
  call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read the grids for the current zone
  !
  do ig_g = 1,size(z%grid)
    if (debug_cgns) write (iout,2) mypnum,"reading CGNS grid node #",ig_g
    call read_cgns_grid( z%grid(ig_g) )
  end do
  !
  ! Read the sections for the current zone
  !
  do is_g = 1,size(z%section)
    if (debug_cgns) write (iout,2) mypnum,"reading CGNS section node #",is_g
    call read_cgns_section( z%section(is_g) )
  end do
  !
  ! Read the sections for the current zone
  !
  do ibc_g = 1,size(z%bc)
    if (debug_cgns) write (iout,2) mypnum,"reading CGNS bc node #",ibc_g
    call read_cgns_bc( z%bc(ibc_g) )
  end do
  !
  ! Fix the boundary condition lists
  !
 !call fix_boco_lists(z%nvertex,z%ncell,z%section,z%bc)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (" The CGNS grid file seems to contain a structured zone", &
            " which are not currently supported!  Error occured in", &
            " CGNS base # ",i0," , zone # ",i0," !")
  2 format (" CPU #",i0,": ",a,:,i0)
  !
end subroutine read_cgns_zone
!
!###############################################################################
!
subroutine read_cgns_grid(g)
  !
  !.. Formal Arguments ..
  type(cg_grid_t), intent(inout) :: g
  !
  !.. Local Scalars ..
  integer      :: ierr
  integer(CBT) :: ncoords,narrays
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_cgns_grid"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Get the name of the current grid
  !
  call cg_grid_read_f(if_g,ib_g,iz_g,ig_g,g%name,cgierr)
  call cgns_error(pname,gridfile,"cg_grid_read_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  if (ig_g == 1) then
    !
    ! This is for the original grid, so we need to use the cg_coord_* functions
    !
    ! Get the number of coordinates
    !
    call cg_ncoords_f(if_g,ib_g,iz_g,ncoords,cgierr)
    call cgns_error(pname,gridfile,"cg_ncoords_f",cgierr, &
                    __LINE__,__FILE__,not_parallel=true)
    !
    ! Allocate the coord component of the current grid
    !
    allocate ( g%coord(1:ncoords) , stat=ierr , errmsg=error_message )
    array_name = get_array_name("grid",ig_g,"coord")
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Read the coordinate data
    !
    do ic_g = 1,size(g%coord)
      if (debug_cgns) write (iout,1) mypnum,"reading CGNS coord node #",ic_g
      call read_cgns_coord( g%coord(ic_g) )
    end do
    !
  else
    !
    ! This is for an additional grid beyond the original
    ! so we need to use the cg_array_* functions
    !
    ! Because the grid data is stored in a generic DataArray_t, we need to
    ! manually move to the current GridCoordinates_t node to read the data
    !
    call cg_goto_f(if_g,ib_g,cgierr, &
                   "Zone_t",iz_g, &
                   "GridCoordinates_t",ig_g, &
                   CGNS_GOTO_END)
    call cgns_error(pname,gridfile,"cg_goto_f",cgierr, &
                    __LINE__,__FILE__,not_parallel=true)
    !
    ! Get the number of coordinate arrays for the current grid
    !
    call cg_narrays_f(narrays,cgierr)
    call cgns_error(pname,gridfile,"cg_narrays_f",cgierr, &
                    __LINE__,__FILE__,not_parallel=true)
    !
    ! Allocate the coord component of the current grid
    !
    allocate ( g%coord(1:narrays) , stat=ierr , errmsg=error_message )
    array_name = get_array_name("grid",ig_g,"coord")
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Read the coordinate data
    !
    do ia_g = 1,size(g%coord)
      if (debug_cgns) then
        write (iout,1) mypnum,"reading CGNS additional coord node #",ia_g
      end if
      call read_cgns_additional_coord( g%coord(ia_g) )
    end do
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (" CPU #",i0,": ",a,:,i0)
  !
end subroutine read_cgns_grid
!
!###############################################################################
!
subroutine read_cgns_coord(c)
  !
  !.. Formal Arguments ..
  type(cg_coord_t), intent(inout) :: c
  !
  !.. Local Scalars ..
  integer      :: ierr
  integer(CST) :: np,ibeg,iend
  integer(CET) :: stored_datatype
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_cgns_coord"
  integer(CET),     parameter :: read_datatype = merge(REALDOUBLE, &
                                                       REALSINGLE,wp==r8)
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Read the coordinate info
  !
  call cg_coord_info_f(if_g,ib_g,iz_g,ic_g,stored_datatype,c%name,cgierr)
  call cgns_error(pname,gridfile,"cg_coord_info_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  ! Allocate the data array to the total number of nodes in the grid
  !
  np = base(ib_g)%zone(iz_g)%nvertex
  !
  allocate ( c%x(1:np) , source=zero , stat=ierr , errmsg=error_message )
  array_name = get_array_name("grid",ig_g,"coord",ic_g,"x")
  call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read in the coordinate data
  !
  ibeg = 1_CST
  iend = np
  !
  call cg_coord_read_f(if_g,ib_g,iz_g,trim(adjustl(c%name)),read_datatype, &
                       ibeg,iend,c%x,cgierr)
  call cgns_error(pname,gridfile,"cg_coord_read_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_cgns_coord
!
!###############################################################################
!
subroutine read_cgns_additional_coord(c)
  !
  !.. Formal Arguments ..
  type(cg_coord_t), intent(inout) :: c
  !
  !.. Local Scalars ..
  integer      :: ierr
  integer(CST) :: np,ibeg,iend
  integer(CBT) :: data_dimension
  integer(CET) :: stored_datatype
  !
  !.. Local Arrays ..
  integer(CST) :: dimension_vector(1:12) ! NOTE: 12 is the maximum number of
  !                                              data dimensions for a
  !                                              DataArray_t node defined by the
  !                                              CGNS standard
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_cgns_additional_coord"
  integer(CET),     parameter :: read_datatype = merge(REALDOUBLE, &
                                                       REALSINGLE,wp==r8)
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Read the coordinate info
  !
  call cg_array_info_f(ia_g,c%name,stored_datatype,data_dimension, &
                       dimension_vector,cgierr)
  call cgns_error(pname,gridfile,"cg_array_info_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  ! Allocate the data array to the total number of nodes in the grid
  !
  np = product( dimension_vector(1:data_dimension) )
  !
  allocate ( c%x(1:np) , source=zero , stat=ierr , errmsg=error_message )
  array_name = get_array_name("grid",ig_g,"coord",ia_g,"x")
  call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read in the coordinate data
  !
  call cg_array_read_as_f(ia_g,read_datatype,c%x,cgierr)
  call cgns_error(pname,gridfile,"cg_array_read_as_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_cgns_additional_coord
!
!###############################################################################
!
subroutine read_cgns_section(s)
  !
  !.. Use Statements ..
  use iso_c_binding, only : c_null_ptr
  !
  !.. Formal Arguments ..
  type(cg_section_t), intent(inout) :: s
  !
  !.. Local Scalars ..
  integer      :: ierr
  integer(CST) :: np,nelem,element_data_size
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_cgns_section"
  !
 !logical(lk), parameter :: read_parent_data = true
  logical(lk), parameter :: read_parent_data = fals
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Read the info for this element section
  !
  call cg_section_read_f(if_g,ib_g,iz_g,is_g,s%name,s%type,s%start,s%end, &
                         s%nbndry,s%parent_flag,cgierr)
  call cgns_error(pname,gridfile,"cg_section_read_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  ! Get the size of the element data
  !
  call cg_elementdatasize_f(if_g,ib_g,iz_g,is_g,element_data_size,cgierr)
  call cgns_error(pname,gridfile,"cg_elementdatasize_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  np = element_data_size
  !
  allocate ( s%elements(1:np) , stat=ierr , errmsg=error_message )
  array_name = get_array_name("section",is_g,"elements")
  call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  !
  if (s%parent_flag == 1_CBT .and. read_parent_data) then
    !
    nelem = s%end - s%start + 1_CST
    np = nelem*4_CST
    !
    allocate ( s%parentdata(1:np) , stat=ierr , errmsg=error_message )
    array_name = get_array_name("section",is_g,"parentdata")
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    call cg_elements_read_f(if_g,ib_g,iz_g,is_g,s%elements,s%parentdata,cgierr)
    call cgns_error(pname,gridfile,"cg_elements_read_f",cgierr, &
                    __LINE__,__FILE__,not_parallel=true)
    !
  else
    !
    call cg_elements_read_f(if_g,ib_g,iz_g,is_g,s%elements,c_null_ptr,cgierr)
    call cgns_error(pname,gridfile,"cg_elements_read_f",cgierr, &
                    __LINE__,__FILE__,not_parallel=true)
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_cgns_section
!
!###############################################################################
!
subroutine read_cgns_bc(bc)
  !
  !.. Formal Arguments ..
  type(cg_bc_t), intent(inout) :: bc
  !
  !.. Local Scalars ..
  integer      :: ierr
  integer(CST) :: np
  !
  !.. Local Allocatable Arrays ..
  real(r4),     allocatable :: normal_list_sp(:)
  real(r8),     allocatable :: normal_list_dp(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_cgns_bc"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Read the info for this boundary condition
  !
  call cg_boco_info_f(if_g,ib_g,iz_g,ibc_g, &
                      bc%name,bc%boco_type,bc%ptset_type,bc%npnts, &
                      bc%normal_index,bc%normal_list_size,bc%normal_datatype, &
                      bc%ndataset,cgierr)
  call cgns_error(pname,gridfile,"cg_boco_info_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  ! Read the location of this boundary condition
  !
  call cg_boco_gridlocation_read_f(if_g,ib_g,iz_g,ibc_g,bc%location,cgierr)
  call cgns_error(pname,gridfile,"cg_boco_gridlocation_read_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
 !write (*,*)
 !write (*,'(a,i0)') "base(ib_g)%zone(iz_g)%index_dim = ", &
 !                    base(ib_g)%zone(iz_g)%index_dim
 !write (*,*)
 !write (*,'(a,a)')  "bc%name             = ",bc%name
 !write (*,'(a,i0)') "bc%boco_type        = ",bc%boco_type
 !write (*,'(a,i0)') "bc%ptset_type       = ",bc%ptset_type
 !write (*,'(a,i0)') "bc%npnts            = ",bc%npnts
 !write (*,'(a,*(1x,i0))') "bc%normal_index     = ",bc%normal_index
 !write (*,'(a,i0)') "bc%normal_list_size = ",bc%normal_list_size
 !write (*,'(a,i0)') "bc%ndataset         = ",bc%ndataset
 !write (*,'(a,i0)') "bc%location         = ",bc%location
 !write (*,'(a,i0)') "bc%normal_datatype  = ",bc%normal_datatype
  !
  ! Make some adjustments if ptset_type is using either of the depreciated
  ! PointSetType values of ElementRange or ElementList. If it is, switch
  ! ptset_type to the respective PointRange or PointList and switch location
  ! to either FaceCenter or EdgeCenter depending on base(ib_g)%cell_dim.
  !
  if (any(bc%ptset_type == [ELEMENTRANGE,ELEMENTLIST])) then
    if (bc%ptset_type == ELEMENTRANGE) then
      bc%ptset_type = POINTRANGE
    else if (bc%ptset_type == ELEMENTLIST) then
      bc%ptset_type = POINTLIST
    end if
    if (base(ib_g)%cell_dim == 3) then
      bc%location = FACECENTER
    else if (base(ib_g)%cell_dim == 2) then
      bc%location = EDGECENTER
    end if
  end if
  !
  ! FOR NOW: Abort if the BC location is not FaceCenter or EdgeCenter
  !
  if (base(ib_g)%cell_dim == 3 .and. bc%location /= FACECENTER) then
    write (error_message,1) 1,ibc_g,trim(adjustl(bc%name)),"Face"
    call stop_recon(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  if (base(ib_g)%cell_dim == 2 .and. bc%location /= EDGECENTER) then
    write (error_message,1) 2,ibc_g,trim(adjustl(bc%name)),"Edge"
    call stop_recon(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Allocate the pnts component for this boundary condition
  !
  np = base(ib_g)%zone(iz_g)%index_dim * bc%npnts
  !
  allocate ( bc%pnts(1:np) , source=0_CST , stat=ierr , errmsg=error_message )
  array_name = get_array_name("bc",ibc_g,"pnts")
  call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  !
  if (bc%normal_list_size > 0) then
    !
    ! Allocate the normal_list component for this boundary condition
    !
    np = base(ib_g)%phys_dim * bc%normal_list_size
    !
    allocate ( bc%normal_list(1:np), source=zero , &
               stat=ierr , errmsg=error_message )
    array_name = get_array_name("bc",ibc_g,"normal_list")
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Read the boundary condition data depending on the datatype of normal_list
    !
    if (bc%normal_datatype == REALSINGLE) then
      !
      allocate ( normal_list_sp(1:np) , source=0.0_r4 , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"normal_list_sp",1,__LINE__,__FILE__,ierr, &
                       error_message)
      !
      call cg_boco_read_f(if_g,ib_g,iz_g,ibc_g,bc%pnts,normal_list_sp,cgierr)
      call cgns_error(pname,gridfile,"cg_boco_read_f(SP)",cgierr, &
                      __LINE__,__FILE__,not_parallel=true)
      !
      bc%normal_list(:) = real( normal_list_sp(:) , kind=wp )
      !
      deallocate ( normal_list_sp , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"normal_list_sp",2,__LINE__,__FILE__,ierr, &
                       error_message)
      !
    else if (bc%normal_datatype == REALDOUBLE) then
      !
      allocate ( normal_list_dp(1:np) , source=0.0_r8 , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"normal_list_dp",1,__LINE__,__FILE__,ierr, &
                       error_message)
      !
      call cg_boco_read_f(if_g,ib_g,iz_g,ibc_g,bc%pnts,normal_list_dp,cgierr)
      call cgns_error(pname,gridfile,"cg_boco_read_f(DP)",cgierr, &
                      __LINE__,__FILE__,not_parallel=true)
      !
      bc%normal_list(:) = real( normal_list_dp(:) , kind=wp )
      !
      deallocate ( normal_list_dp , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"normal_list_dp",2,__LINE__,__FILE__,ierr, &
                       error_message)
      !
    end if
    !
  else
    !
    call cg_boco_read_f(if_g,ib_g,iz_g,ibc_g,bc%pnts,c_null_ptr,cgierr)
    call cgns_error(pname,gridfile,"cg_boco_read_f(null)",cgierr, &
                    __LINE__,__FILE__,not_parallel=true)
    !
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (" The CGNS grid file appears to be ",i0,"D but boundary", &
            " condition # ",i0,' ( "',a,'" ) is not specified to be', &
            " ",a,"Centered. Currently, the only supported locations", &
            " for boundary conditions are FaceCentered for 3D and", &
            " EdgeCentered for 2D.")
  !
end subroutine read_cgns_bc
!
!###############################################################################
!
pure function get_array_name(c1,n1,c2,n2,c3) result(return_value)
  !
  character(len=*), optional, intent(in) :: c1
  integer,          optional, intent(in) :: n1
  character(len=*), optional, intent(in) :: c2
  integer,          optional, intent(in) :: n2
  character(len=*), optional, intent(in) :: c3
  !
  character(len=len(array_name)) :: return_value
  !
  integer :: i
  !
continue
  !
  write (return_value,1) ib_g
  !
  one_pass_loop: do i = 1,1
    !
    if (present(c1)) then
      write (return_value,2) trim(adjustl(return_value)),iz_g, &
                             trim(adjustl(c1))
    else
      exit one_pass_loop
    end if
    !
    if (present(n1)) then
      write (return_value,3) trim(adjustl(return_value)),n1
    else
      exit one_pass_loop
    end if
    !
    if (present(c2)) then
      write (return_value,4) trim(adjustl(return_value)),trim(adjustl(c2))
    else
      exit one_pass_loop
    end if
    !
    if (present(n2)) then
      write (return_value,3) trim(adjustl(return_value)),n2
    else
      exit one_pass_loop
    end if
    !
    if (present(c3)) then
      write (return_value,4) trim(adjustl(return_value)),trim(adjustl(c3))
    else
      exit one_pass_loop
    end if
    !
  end do one_pass_loop
  !
  return_value = trim(adjustl(return_value))
  !
  ! Format Statements
  !
  1 format ("base(",i0,")%zone")
  2 format (a,"(",i0,")%",a)
  3 format (a,"(",i0,")")
  4 format (a,"%",a)
  !
end function get_array_name
!
!###############################################################################
!
subroutine transfer_cgns_data_to_grid_dt
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,n3,np,nc,nf
  integer :: i,j,j0,j1,j2
  integer :: ib,iz,ig,ic,ix,is,ibc
  integer :: elem_beg,elem_end
  integer :: mine,maxe,sume,ierr
  integer :: ngrids,nsections,nbocos
  integer(CET) :: sect_type,elem_type
  !
  character(len=CGLEN) :: cname
  character(len=CGLEN) :: bc_string
  !
  !.. Local Arrays ..
  integer :: idx(1:3)
  integer :: idxr(1:3)
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: mske(:)
  character(len=32), allocatable :: unsupported_bc(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "transfer_cgns_data_to_grid_dt"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! First, start by allocating the variable grid 
  !
  allocate ( grid , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Get the size needed for the grid%xyz array
  !
  call debug_timer(start_timer,"finding size for grid%xyz")
  !
  ngrids = 0
  n1 = 0
  n2 = 0
  !
  do ib = 1,size(base)
    n1 = max( n1 , base(ib)%phys_dim )
    do iz = 1,size(base(ib)%zone)
      do ig = 1,size(base(ib)%zone(iz)%grid)
        if (allocated(base(ib)%zone(iz)%grid(ig)%coord)) then
          ngrids = ngrids + 1
          n2 = n2 + size(base(ib)%zone(iz)%grid(ig)%coord(1)%x)
        end if
      end do
    end do
  end do
  !
  ! Create the grid%xyz array
  !
  call debug_timer(stop_timer,"finding size for grid%xyz")
  call debug_timer(start_timer,"creating grid%xyz")
  !
  allocate ( grid%xyz(1:n1,1:n2) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid%xyz",1,__LINE__,__FILE__,ierr,error_message)
  !
  n = 0
  do ib = 1,size(base)
    do iz = 1,size(base(ib)%zone)
      do ig = 1,size(base(ib)%zone(iz)%grid)
        if (allocated(base(ib)%zone(iz)%grid(ig)%coord)) then
          !
          idx(1:3) = [1,2,3]
          do ic = 1,size(base(ib)%zone(iz)%grid(ig)%coord)
            cname = trim(adjustl(base(ib)%zone(iz)%grid(ig)%coord(ic)%name))
            if (cname == "CoordinateX" .or. cname == "X") then
              idx(ic) = 1
              idxr(1) = ic
            else if (cname == "CoordinateY" .or. cname == "Y") then
              idx(ic) = 2
              idxr(2) = ic
            else if (cname == "CoordinateZ" .or. cname == "Z") then
              idx(ic) = 3
              idxr(3) = ic
            end if
          end do
          !
          n1 = size(base(ib)%zone(iz)%grid(ig)%coord)
          n2 = size(base(ib)%zone(iz)%grid(ig)%coord(1)%x)
          !
         !call debug_timer(start_timer,"creating grid%xyz, crazy order")
         !grid%xyz(1:n1,n+1:n+n2) = transpose( &
         !   reshape( [ (base(ib)%zone(iz)%grid(ig)%coord(idxr(ic))%x, &
         !               ic=1,size(base(ib)%zone(iz)%grid(ig)%coord)) ], &
         !            [n2,n1] ) )
         !call debug_timer(stop_timer,"creating grid%xyz, crazy order")
          !
         !call debug_timer(start_timer,"creating grid%xyz, old order")
          do ix = 1,n2
            n = n + 1
            do ic = 1,n1
              grid%xyz(idx(ic),n) = base(ib)%zone(iz)%grid(ig)%coord(ic)%x(ix)
            end do
          end do
         !call debug_timer(stop_timer,"creating grid%xyz, old order")
          !
        end if
      end do
    end do
  end do
  !
  ! Get the size needed for the grid%elem array
  !
  call debug_timer(stop_timer,"creating grid%xyz")
  call debug_timer(start_timer,"find size for grid%elem")
  !
  nsections = 0
  !
  n2 = 0
  !
  do ib = 1,size(base)
    do iz = 1,size(base(ib)%zone)
      n1 = 0
      do is = 1,size(base(ib)%zone(iz)%section)
        nsections = nsections + 1
        elem_end = int( base(ib)%zone(iz)%section(is)% &
                        end , kind=kind(elem_end) )
        n1 = max( n1 , elem_end )
        ! IGNORE PARENT DATA FOR NOW
      end do
      base(ib)%zone(iz)%nsize = int( n1 , kind=kind(base(ib)%zone(iz)%nsize) )
    end do
    base(ib)%nsize = sum( base(ib)%zone(:)%nsize )
    n2 = n2 + int( base(ib)%nsize , kind=kind(n2) )
  end do
  !
  ! Allocate the mske array
  !
  call debug_timer(stop_timer,"find size for grid%elem")
  call debug_timer(start_timer,"checking for only unique elements")
  !
  allocate ( mske(1:n2) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"mske",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Check to make sure that each element defined in a section is
  ! defined once and only once across all the CGNS bases and zones
  ! NOTE: Under a given zone, all elements (whether its a cell, a face,
  !       or an edge) are required by the CGNS standard to have a unique
  !       index within that zone, with the indexing starting at 1.
  !       Additionally, within a given element section of a zone, the
  !       element indices are required to be continuous.
  !
  n2 = 0
  !
  do ib = 1,size(base)
    do iz = 1,size(base(ib)%zone)
      do is = 1,size(base(ib)%zone(iz)%section)
        elem_beg = n2 + int( base(ib)%zone(iz)%section(is)% &
                             start , kind=kind(elem_end) )
        elem_end = n2 + int( base(ib)%zone(iz)%section(is)% &
                             end , kind=kind(elem_end) )
        do n = elem_beg,elem_end
          mske(n) = mske(n) + 1
        end do
      end do
      n2 = n2 + int( base(ib)%zone(iz)%nsize , kind=kind(n2) )
    end do
  end do
  !
  mine = minval(mske)
  maxe = maxval(mske)
  sume = sum(mske)
  !
  ! Check for any errors, and abort if any were found
  !
  if (mine /= 1 .or. maxe /= 1 .or. sume /= size(mske)) then
    if (mine /= 1 .or. maxe /= 1) then
      if (mypnum == 0) write (iout,11) mine,maxe
    end if
    if (sume /= size(mske)) then
      if (mypnum == 0) write (iout,12) sume,size(mske)
    end if
    write (error_message,13)
    call stop_recon(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Deallocate mske
  !
  deallocate ( mske , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"mske",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the grid%elem array
  !
  call debug_timer(stop_timer,"checking for only unique elements")
  call debug_timer(start_timer,"creating grid%elem")
  !
  ! Allocate the elem component of the grid derived type
  !
  allocate ( grid%elem(1:sume) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid%elem",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Store the CGNS indices (base, zone, section) for all the elements
  ! in the CGNS grid file for easy access later.
  !
  n2 = 0
  !
  do ib = 1,size(base)
    do iz = 1,size(base(ib)%zone)
      do is = 1,size(base(ib)%zone(iz)%section)
        !
        sect_type = base(ib)%zone(iz)%section(is)%type
        elem_beg = n2 + int( base(ib)%zone(iz)%section(is)% &
                             start , kind=kind(elem_end) )
        elem_end = n2 + int( base(ib)%zone(iz)%section(is)% &
                             end , kind=kind(elem_end) )
        !
        if (sect_type == MIXED) then
          !
          j2 = 0
          do n = elem_beg,elem_end
            !
            j0 = j2 + 1
            !
            elem_type = int( base(ib)%zone(iz)%section(is)% &
                             elements(j0) , kind=kind(elem_type) )
            np = npts_for_element_type(elem_type)
            !
            j1 = j0 + 1
            j2 = j0 + np
            !
            grid%elem(n) = &
                create_element(elem_type, &
                               base(ib)%zone(iz)%section(is)%elements(j1:j2))
            !
           !grid%elem(n)%idx = cgns_idx_t(ib,iz,is)
           !allocate ( grid%elem(n)%idx , source=cgns_idx_t(ib,iz,is) , &
           !           stat=ierr , errmsg=error_message )
           !write (array_name,1) n,"idx"
           !call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
           !                 error_message)
            !
            grid%elem(n)%tags = [0] ! using F2003 automatic reallocation
           !allocate ( grid%elem(n)%tags(1:1) , source=0 , &
           !           stat=ierr , errmsg=error_message )
           !write (array_name,1) n,"tags"
           !call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
           !                 error_message)
            !
          end do
          !
        else
          !
          elem_type = sect_type
          np = npts_for_element_type(elem_type)
          !
          j2 = 0
          do n = elem_beg,elem_end
            !
            j1 = j2 + 1
            j2 = j2 + np
            !
            grid%elem(n) = &
                create_element(elem_type, &
                               base(ib)%zone(iz)%section(is)%elements(j1:j2))
            !
           !grid%elem(n)%idx = cgns_idx_t(ib,iz,is)
           !allocate ( grid%elem(n)%idx , source=cgns_idx_t(ib,iz,is) , &
           !           stat=ierr , errmsg=error_message )
           !write (array_name,1) n,"idx"
           !call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
           !                 error_message)
            !
            grid%elem(n)%tags = [0] ! using F2003 automatic reallocation
           !allocate ( grid%elem(n)%tags(1:1) , source=0 , &
           !           stat=ierr , errmsg=error_message )
           !write (array_name,1) n,"tags"
           !call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
           !                 error_message)
            !
          end do
          !
        end if
        !
      end do
      n2 = n2 + int( base(ib)%zone(iz)%nsize , kind=kind(n2) )
    end do
  end do
  !
  ! Add the CGNS boco index to all the elements given a boundary condition
  !
  call debug_timer(stop_timer,"creating grid%elem")
  call debug_timer(start_timer,"adding BC to grid%elem")
  !
  n2 = 0
  nbocos = 0
  !
  do ib = 1,size(base)
    do iz = 1,size(base(ib)%zone)
      do ibc = 1,size(base(ib)%zone(iz)%bc)
        !
        nbocos = nbocos + 1
        !
        if (base(ib)%zone(iz)%bc(ibc)%ptset_type == POINTLIST) then
          !
          do i = 1,size(base(ib)%zone(iz)%bc(ibc)%pnts)
            !
            n = n2 + int( base(ib)%zone(iz)%bc(ibc)% &
                          pnts(i) , kind=kind(elem_end) )
            !
           !grid%elem(n)%idx%boco = ibc
            grid%elem(n)%tags(1) = nbocos
            !
          end do
          !
        else if (base(ib)%zone(iz)%bc(ibc)%ptset_type == POINTRANGE) then
          !
          elem_beg = n2 + int( base(ib)%zone(iz)%bc(ibc)% &
                               pnts(1) , kind=kind(elem_end) )
          elem_end = n2 + int( base(ib)%zone(iz)%bc(ibc)% &
                               pnts(2) , kind=kind(elem_end) )
          !
          do n = elem_beg,elem_end
           !grid%elem(n)%idx%boco = ibc
            grid%elem(n)%tags(1) = nbocos
          end do
          !
        end if
        !
      end do
      !
      n2 = n2 + int( base(ib)%zone(iz)%nsize , kind=kind(n2) )
      !
    end do
  end do
  !
  ! Create the cell_grps derived type
  !
  call debug_timer(stop_timer,"adding BC to grid%elem")
  call debug_timer(start_timer,"creating cell groups")
  !
  allocate ( cell_grps(1:1+nbocos) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_grps",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! All interior across all zones
  !
  n = 1
  cell_grps(n)%name     = "Interior"
  cell_grps(n)%cell_dim = maxval( base(:)%cell_dim )
  cell_grps(n)%grp_id   = 0
  cell_grps(n)%ibc = NOT_A_BC
  !
  ! Add all the BCs
  !
  do ib = 1,size(base)
    do iz = 1,size(base(ib)%zone)
      do ibc = 1,size(base(ib)%zone(iz)%bc)
        !
        n = n + 1
        !
        bc_string = BCTYPENAME( base(ib)%zone(iz)%bc(ibc)%boco_type )
        if (uppercase(bc_string) == 'NULL') then
          bc_string = trim(adjustl(base(ib)%zone(iz)%bc(ibc)%name))
        end if
        !
        !#################################################################
        !#################################################################
        !#################################################################
        !    TEMPORARY FOR POINTWISE CURVED MESH OF SMC000 NOZZLE
        !
        if (uppercase(bc_string) == 'FARFIELD') then
          bc_string = 'CHARACTERISTIC'
        else if (uppercase(bc_string) == 'INFLOW') then
          bc_string = 'REFERENCE'
        else if (uppercase(bc_string) == 'INFLOWSUBSONIC') then
          bc_string = 'SUBINFLOW'
        else if (uppercase(bc_string) == 'NOZZLEBASE') then
          bc_string = 'ADIABATIC'
        else if (uppercase(bc_string) == 'NOZZLEEXTERIOR') then
          bc_string = 'ADIABATIC'
        else if (uppercase(bc_string) == 'NOZZLEINTERIOR') then
          bc_string = 'ADIABATIC'
        else if (uppercase(bc_string) == 'OUTFLOW') then
          bc_string = 'CHARACTERISTIC'
        end if
        !#################################################################
        !#################################################################
        !#################################################################
        !
        cell_grps(n)%name = base(ib)%zone(iz)%bc(ibc)%name
        if (base(ib)%zone(iz)%bc(ibc)%location == FACECENTER) then
          cell_grps(n)%cell_dim = 2
        else if (base(ib)%zone(iz)%bc(ibc)%location == EDGECENTER) then
          cell_grps(n)%cell_dim = 1
        end if
        cell_grps(n)%grp_id = n-1
        cell_grps(n)%ibc = bc_string_to_integer( bc_string )
        if (cell_grps(n)%ibc == bc_unsupported) then
          if (allocated(unsupported_bc)) then
            unsupported_bc = [unsupported_bc,trim(adjustl(bc_string))]
          else
            unsupported_bc = [bc_string]
          end if
        end if
      end do
    end do
  end do
  !
  if (allocated(unsupported_bc)) then
    do ibc = 1,size(unsupported_bc)
      write (iout,201) trim(adjustl(unsupported_bc(ibc)))
    end do
  end if
  201 format (" UNSUPPORTED BC = '",a,"'")
  !
  ! Save the total number of interior grid cells, grid vertices, and
  ! boundary faces that the CGNS grid file says that there should be,
  ! so that we can compare these values to those found later when
  ! creating all the geometry arrays for the solver.
  !
  call debug_timer(stop_timer,"creating cell groups")
  call debug_timer(start_timer,"saving grid parameters")
  !
  npnts_correct = 0
  ncell_correct = 0
  nfbnd_correct = 0
  !
  do ib = 1,size(base)
    do iz = 1,size(base(ib)%zone)
      np = int( base(ib)%zone(iz)%nvertex , kind=kind(np) )
      nc = int( base(ib)%zone(iz)%ncell , kind=kind(nc) )
      nf = int( base(ib)%zone(iz)%nsize , kind=kind(nf) ) - nc
      npnts_correct = npnts_correct + np
      ncell_correct = ncell_correct + nc
      nfbnd_correct = nfbnd_correct + nf
    end do
  end do
  !
  write (iout,2) npnts_correct,ncell_correct,nfbnd_correct
  call debug_timer(stop_timer,"saving grid parameters")
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1  format ("grid%elem(",i0,")%",a)
  2  format (/,3x,"The total dimensions of the CGNS grid file are:",/, &
               10x,"total # of  grid vertices = ",i0,/, &
               10x,"total # of interior cells = ",i0,/, &
               10x,"total # of boundary faces = ",i0,/)
  !
  11 format (" ERROR!!! mine = ",i0,"; maxe = ",i0)
  12 format (" ERROR!!! sume = ",i0,"; size(mske) = ",i0)
  13 format (" THERE APPEARS TO BE AN ERROR WITH THE ELEMENT DEFINITIONS", &
             " WITHIN THE CGNS GRID FILE. ALL ELEMENTS SHOULD BE DEFINED", &
             " ONCE AND ONLY ONCE, BUT THIS DOES NOT SEEM TO BE THE CASE!!")
  !
 !101 format (a," : Finding the size needed for the grid%xyz array")
 !102 format (a," : Creating the grid%xyz array")
 !103 format (a," : Finding the size needed for the grid%elem array")
 !104 format (a," : Checking to make sure that all elements",/, &
 !            31x," are defined once and only once")
 !105 format (a," : Creating the grid%elem array")
 !106 format (a," : Adding BC information into grid%elem for boundary faces")
 !107 format (a," : Creating the cell groups")
 !108 format (a," : Saving the total number of grid cells, grid vertices",/, &
 !            31x," and boundary faces defined by the CGNS grid file")
  !
end subroutine transfer_cgns_data_to_grid_dt
!
!###############################################################################
!
subroutine create_solver_geom_arrays_cgns
  !
  ! Procedure to create the remaining grid geometry arrays using
  ! the information stored in the arrays grid%elem and cell_grps.
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use geovar,    only : nr,ncell,nnode,nfbnd,nbfai,bface
  use geovar,    only : nodes_of_cell_ptr,nodes_of_cell
  use geovar,    only : xyz_nodes
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
  integer :: int_nodes(1:8)
  integer :: face_nodes(1:4)
  integer :: marker(1:8)
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: int_grps(:)
  integer, allocatable :: bnd_grps(:)
  integer, allocatable :: b0_grps(:)
  integer, allocatable :: mapp(:)
  integer, allocatable :: mapr(:)
  integer, allocatable :: lpoin(:)
  integer, allocatable :: cells_with_node_ptr(:)
  integer, allocatable :: cells_with_node(:)
  !
  logical(lk), allocatable :: mskp(:)
  !
  type(elem_t), allocatable :: tmp_cell_data(:)
  type(cell_grps_t), allocatable :: tmp_cell_grps(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_solver_geom_arrays_cgns"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! NOTE: UNLIKE THE SIMILAR SUBROUTINE IN THE GMSH MODULE, THIS SUBROUTINE
  !       WILL ONLY BE CALLED BY THE ROOT PROCESSOR IN ORDER TO SAVE MEMORY.
  !       THE ROOT PROCESSOR WILL HAVE TO BROADCAST THE NECESSARY ARRAYS
  !       (e.g., nodes_of_cell_ptr, nodes_of_cell, xyz_nodes, bface, etc.)
  !       TO ALL THE OTHER PROCESSORS AFTER RETURNING FROM THIS SUBROUTINE.
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
#ifdef DEBUG_ON
  write (iout,101) pname
#endif
  !
  ! Get the physical directional dimensions from the CGNS bases
  !
  nr = maxval( base(:)%phys_dim )
  !
  ! Count the number of interior cell groups
  !
  igrps = count(cell_grps(:)%cell_dim == nr)
  !
  ! Count the number of boundary cell groups
  ! NOTE: Only count the 0-dimensional cell groups if the interior is 1D
  !
  if (nr > 1) then
    bgrps = count( cell_grps(:)%cell_dim < nr .and. &
                   cell_grps(:)%cell_dim > 0 )
    b0grps = count(cell_grps(:)%cell_dim == 0)
  else
    bgrps = count(cell_grps(:)%cell_dim < nr)
    b0grps = 0
  end if
  !
  ! Allocate the arrays to identify the interior and boundary cell groups
  !
  allocate ( int_grps(1:igrps) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"int_grps",1,__LINE__,__FILE__,ierr,error_message)
  allocate ( bnd_grps(1:bgrps) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bnd_grps",1,__LINE__,__FILE__,ierr,error_message)
  allocate ( b0_grps(1:b0grps) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"b0_grps",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Collect the interior cell groups
  !
  int_grps(1:igrps) = pack( intseq(1,size(cell_grps)) , &
                            cell_grps(:)%cell_dim == nr )
  !
  ! Collect the boundary cell groups
  ! NOTE: Only collect the 0-dimensional cell groups if the interior is 1D
  !
  if (nr > 1) then
    bnd_grps(1:bgrps) = pack( intseq(1,size(cell_grps)) , &
                              cell_grps(:)%cell_dim < nr .and. &
                              cell_grps(:)%cell_dim > 0 )
    b0_grps(1:b0grps) = pack( intseq(1,size(cell_grps)) , &
                              cell_grps(:)%cell_dim == 0 )
  else
    bnd_grps(1:bgrps) = pack( intseq(1,size(cell_grps)) , &
                              cell_grps(:)%cell_dim < nr )
  end if
  !
  ! Check that the total number of interior and boundary groups is correct
  !
  if ( igrps+bgrps+b0grps /= size(cell_grps) ) then
    write (iout,11) igrps,bgrps,b0grps,size(cell_grps)
    call stop_recon(abort,pname,__LINE__,__FILE__)
  end if
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
#ifdef DEBUG_ON
  write (iout,102) pname
#endif
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
        if (any(grid%elem(j)%tags(1) == cell_grps(b0_grps(:))%grp_id)) then
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
      if (any(grid%elem(j)%tags(1) == cell_grps(int_grps(:))%grp_id)) then
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
      if (any(grid%elem(j)%tags(1) == cell_grps(bnd_grps(:))%grp_id)) then
        n = n + 1
        tmp_cell_data(n) = grid%elem(j)
      end if
    end if
  end do
  !
  if (n /= size(grid%elem) - b0cells) then
    write (iout,21) size(grid%elem),b0cells,size(grid%elem)-b0cells,n
    call stop_recon(abort,pname,__LINE__,__FILE__)
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
  ! Remove any 0-D cell groups if the interior is not 1D
  !
  if (b0grps > 0) then
    !
#ifdef DEBUG_ON
    write (iout,103) pname
#endif
    !
    ! Allocate the temporary cell groups derived type
    !
    n = size(cell_grps) - b0grps
    allocate ( tmp_cell_grps(1:n) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"tmp_cell_grps",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    ! Save all the non-0-D cell groups to
    ! the temporary cell groups derived type
    !
    n = 0
    do i = 1,size(cell_grps)
      if (cell_grps(i)%cell_dim > 0) then
        n = n + 1
        tmp_cell_grps(n) = cell_grps(i)
      end if
    end do
    !
    if (n /= size(cell_grps) - b0grps) then
      write (iout,22) size(cell_grps),b0grps,size(cell_grps)-b0grps,n
      call stop_recon(abort,pname,__LINE__,__FILE__)
    end if
    !
    ! Use move_alloc to rebuild cell_grps without any 0-D cell groups
    !
    call move_alloc( from=tmp_cell_grps , to=cell_grps )
    !
    ! Finally, rebuild the int_grps and bnd_grps arrays since the
    ! index of these groups within cell_grps may have now changed.
    ! NOTE: The number of interior and boundary groups should not have changed
    !       so we do not need to reallocate these arrays to a different size
    !
    int_grps(1:igrps) = pack( intseq(1,size(cell_grps)) , &
                              cell_grps(:)%cell_dim == nr )
    bnd_grps(1:bgrps) = pack( intseq(1,size(cell_grps)) , &
                              cell_grps(:)%cell_dim < nr )
    !
  end if
  !
  ! #######
  ! ## 4 ##
  ! #######
  !
#ifdef DEBUG_ON
  write (iout,104) pname
#endif
  !
  ! Loop through the cell groups and collect
  ! the elements for each cell group
  !
  do i = 1,size(cell_grps)
    !
    pid = cell_grps(i)%grp_id
    pdm = cell_grps(i)%cell_dim
    !
    ! Count the number of elements that match the id number for this cell group
    !
    inum = 0
    do j = 1,size(grid%elem)
      if (allocated(grid%elem(j)%tags)) then
        inum = inum + merge( 1 , 0 , grid%elem(j)%tags(1) == pid )
      end if
    end do
   !inum = count( grid%elem(:)%tags(1) == pid )
    !
    ! Allocate the cells array for this cell group
    !
    write (array_name,1) "cell_grps(", i, ")%cells"
    allocate ( cell_grps(i)%cells(1:inum) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Create the list of elements that are in this cell group
    !
    inum = 0
    do j = 1,size(grid%elem)
      if (allocated(grid%elem(j)%tags)) then
        if (grid%elem(j)%tags(1) == pid) then
          inum = inum + 1
          cell_grps(i)%cells(inum) = j
        end if
      end if
    end do
   !cell_grps(i)%cells(1:inum) = pack( intseq(1,size(grid%elem)) , &
   !                                   grid%elem(:)%tags(1) == pid )
    !
  end do
  !
  ! #######
  ! ## 5 ##
  ! #######
  !
#ifdef DEBUG_ON
  write (iout,105) pname
#endif
  !
  ! Count the total number of interior cells and boundary faces
  !
  ncell = 0
  do i = 1,size(int_grps)
    ncell = ncell + size( cell_grps(int_grps(i))%cells )
  end do
  !
  if (ncell /= ncell_correct) then
    write (iout,12) ncell,ncell_correct
    call stop_recon(abort,pname,__LINE__,__FILE__)
  end if
  !
  nfbnd = 0
  do i = 1,size(bnd_grps)
    nfbnd = nfbnd + size( cell_grps(bnd_grps(i))%cells )
  end do
  !
  if (nfbnd /= nfbnd_correct) then
    write (iout,13) nfbnd,nfbnd_correct
    call stop_recon(abort,pname,__LINE__,__FILE__)
  end if
  !
  ! Check that the number of interior cells and boundary faces are correct
  !
  if ( ncell+nfbnd /= size(grid%elem) ) then
    write (iout,14) ncell,nfbnd,size(grid%elem)
    call stop_recon(abort,pname,__LINE__,__FILE__)
  end if
  !
  ! #######
  ! ## 6 ##
  ! #######
  !
#ifdef DEBUG_ON
  write (iout,106) pname
#endif
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
  write (iout,*)
  write (iout,1) "ncell = ",ncell
  write (iout,1) "nfbnd = ",nfbnd
  write (iout,1) "nnode = ",nnode
  write (iout,1) "npts  = ",npts
  write (iout,*)
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
#ifdef DEBUG_ON
  write (iout,107) pname
#endif
  !
  ! Allocate the nodes_of_cell_ptr array
  !
  allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__,ierr, &
                   error_message)
  !
  ! NOTE: Unlike the GMSH version of this subroutine, we will wait to allocate
  !       cell_geom and cell_order until the subroutine that broadcasts all the
  !       geometry arrays to all the non-root processors
  !
  ! Loop through the cells for each interior cell
  ! group and use grid%elem to create nodes_of_cell_ptr
  !
 !n = 0
 !do i = 1,size(int_grps)
 !  np = int_grps(i)
 !  do j = 1,size(cell_grps(np)%cells)
 !    n = n + 1
 !    nc = cell_grps(np)%cells(j)
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
    do j = 1,size(cell_grps(np)%cells)
      n = cell_grps(np)%cells(j)
     !cell_geom(n) = grid%elem(n)%prop%geom
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
#ifdef DEBUG_ON
  write (iout,108) pname
#endif
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
#ifdef DEBUG_ON
  write (iout,109) pname
#endif
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
  allocate ( lpoin(1:nnode) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"lpoin",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Find the host cells for each boundary face
  !
  do i = 1,size(bnd_grps)
    !
#ifdef DEBUG_ON
   !write (iout,1091) i 
#endif
    !
    np = bnd_grps(i)
    !
    face_loop: do j = 1,size(cell_grps(np)%cells)
      !
      n = cell_grps(np)%cells(j)
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
#ifdef DEBUG_ON
  write (iout,110) pname
#endif
  !
  ! Create and fill out bface with data from the
  ! ZoneBC_t node of the CGNS grid file
  !
  ! Allocate the bface array
  !
  allocate ( bface(1:nbfai,1:nfbnd) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Loop through the faces for each boundary cell
  ! group and use grid%elem to create bface
  !
  nf = 0
  boundary_cell_grps: do i = 1,size(bnd_grps)
    !
    np = bnd_grps(i)
    !
    current_cell_grp: do j = 1,size(cell_grps(np)%cells)
      !
      bnd_cell = cell_grps(np)%cells(j)
      host_cell = grid%elem(bnd_cell)%host_cell
      nf = nf + 1
      !
      bface( 1,nf) = cell_grps(np)%ibc       ! face boundary condition
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
          write (iout,99)
          error_found = true
        end if
        !
        ! Write out the specific error type corresponding
        ! to the negative value of ifa
        !
        if (ifa == -1) then
          write (iout,31) "A 0-D 'Point'-type",host_cell,bnd_cell
        else if (ifa == -2) then
          write (iout,31) "A 1-D 'Edge'-type",host_cell,bnd_cell
        else if (ifa == -10) then
          write (iout,31) "An unknown type of",host_cell,bnd_cell
        else if (ifa == -99) then
          write (iout,32) host_cell,bnd_cell
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
    end do current_cell_grp
    !
  end do boundary_cell_grps
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
#ifdef DEBUG_ON
  write (iout,111) pname
#endif
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
  ! ###  WHEN USING CGNS GRID FILES! GAMBIT NEUTRAL GRID FILES  ###
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
    write (iout,99)
    call stop_recon(abort,pname,__LINE__,__FILE__)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1  format (a,i0,:,a)
  !
  11 format (/,"************************ ERROR! **************************", &
             /,"  The number of interior and boundary cell groups is", &
             /,"  inconsistent with the total number of cell groups!", &
             /,"     number of     interior groups = ",i0, &
             /,"     number of     boundary groups = ",i0, &
             /,"     number of 0-D boundary groups = ",i0, &
             /,"            total number of groups = ",i0, &
             /,"************************ ERROR! **************************",/)
  12 format (/,"************************** ERROR! **************************", &
             /,"  The number of interior cells identified by the cell", &
             /,"  groups is inconsistent with the total number of interior", &
             /,"  cells that the CGNS grid file says should exist!", &
             /,"     # of interior cells from    cell groups = ",i0, &
             /,"     # of interior cells from CGNS grid file = ",i0, &
             /,"************************** ERROR! **************************",/)
  13 format (/,"************************** ERROR! **************************", &
             /,"  The number of boundary faces identified by the cell", &
             /,"  groups is inconsistent with the total number of boundary", &
             /,"  cells that the CGNS grid file says should exist!", &
             /,"     # of boundary faces from    cell groups = ",i0, &
             /,"     # of boundary faces from CGNS grid file = ",i0, &
             /,"************************** ERROR! **************************",/)
  14 format (/,"************************ ERROR! ************************", &
             /,"  The number of interior cells and boundary faces is", &
             /,"  inconsistent with the total number of grid elements!", &
             /,"          number of interior cells = ",i0, &
             /,"          number of boundary faces = ",i0, &
             /,"     total number of grid elements = ",i0, &
             /,"************************ ERROR! ************************",/)
  !
  31 format ("  ",a," host element was found on a boundary!",/, &
             "     Host Cell = ",i0,";    Boundary Face = ",i0,/, &
             "  This element type is not currently supported!",/)
  32 format ("  Unable to find which face of host cell ",i0,/, &
             "  that makes up the boundary face ",i0," !",/)
  !
  21 format ("  The final counter for tmp_cell_data does not",/, &
             "  match the size of the grid%elem array!",/, &
             "                             size(grid%elem) = ",i0,/, &
             "                 number of 0-D grid elements = ",i0,/, &
             "         size(grid%elem) - 0-D grid elements = ",i0,/, &
             "                 final tmp_cell_data counter = ",i0,/)
  22 format ("  The final counter for tmp_cell_grps does not",/, &
             "  match the size of the cell_grps array!",/, &
             "                             size(cell_grps) = ",i0,/, &
             "                   number of 0-D cell groups = ",i0,/, &
             "           size(cell_grps) - 0-D cell groups = ",i0,/, &
             "                 final tmp_cell_grps counter = ",i0,/)
  !
  99 format (/,"**********************************************************", &
             /,"************************ ERROR! **************************", &
             /,"**********************************************************",/)
  !
  101 format (a," : Separating out interior and boundary cell groups")
  102 format (a," : Reordering grid%elem so interior cells are first",/, &
              32x," followed by boundary faces")
  103 format (a," : Removing any 0-D cell groups since the interior is not 1D")
  104 format (a," : Collecting the elements for each cell group")
  105 format (a," : Counting the total number of interior cells",/, &
              32x," and boundary faces")
  106 format (a," : Building the node list for the linear grid")
  107 format (a," : Creating the nodes_of_cell_ptr array")
  108 format (a," : Creating the nodes_of_cell array")
  109 format (a," : Finding the host cells for each boundary face")
  110 format (a," : Creating the bface array")
  111 format (a," : Finishing up any remaining tasks for the grid geometry")
  !
  1091 format ("  Finding host cells for boundary group # ",i0)
  1092 format ("    ERROR finding host cell for boundary face ",i10,", ",i10)
  !
end subroutine create_solver_geom_arrays_cgns
!
!###############################################################################
!
subroutine broadcast_solver_geom_arrays
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use geovar,    only : nr,ncell,nnode,nfbnd,nbfai,bface
  use geovar,    only : nodes_of_cell_ptr,nodes_of_cell
  use geovar,    only : cell_geom,cell_order,xyz_nodes
  use ovar,      only : iop_solpts,iop_flxpts
  !
  !.. Local Scalars ..
  integer :: n,nocsz,ierr
  integer(INT_MPI) :: isz
  !
  !.. Local Arrays ..
  integer :: grid_size(1:5)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "broadcast_solver_geom_arrays"
  integer(INT_MPI), parameter :: root = 0_int_mpi
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Broadcast the grid dimensions
  !
  if (mypnum == 0) then
    grid_size(:) = [ nr , nnode , ncell , nfbnd , size(nodes_of_cell) ]
  end if
  !
  call mpi_bcast(grid_size,size(grid_size,kind=int_mpi),mpi_inttyp, &
                 root,MPI_COMM_WORLD,mpierr)
  if (debug_cgns) write (iout,1) mypnum,(grid_size(n),n=1,size(grid_size))
  !
  nr    = grid_size(1)
  nnode = grid_size(2)
  ncell = grid_size(3)
  nfbnd = grid_size(4)
  nocsz = grid_size(5)
  !
  ! Have all the non-root processors allocate the arrays
  !
  if (mypnum /= root) then
    !
    allocate ( bface(1:nbfai,1:nfbnd) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    allocate ( nodes_of_cell(1:nocsz) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_nodes",1,__LINE__,__FILE__,ierr,error_message)
    !
  end if
  !
  ! Broadcast the bface array
  !
  isz = int( nbfai*nfbnd , kind=int_mpi )
  call mpi_bcast(bface,isz,mpi_inttyp,root,MPI_COMM_WORLD,mpierr)
  !
  ! Broadcast the nodes_of_cell_ptr array
  !
  isz = int( ncell+1 , kind=int_mpi )
  call mpi_bcast(nodes_of_cell_ptr,isz,mpi_inttyp,root,MPI_COMM_WORLD,mpierr)
  !
  ! Broadcast the nodes_of_cell array
  !
  isz = int( nocsz , kind=int_mpi )
  call mpi_bcast(nodes_of_cell,isz,mpi_inttyp,root,MPI_COMM_WORLD,mpierr)
  !
  ! Broadcast the xyz_nodes array
  !
  isz = int( nr*nnode , kind=int_mpi )
  call mpi_bcast(xyz_nodes,isz,mpi_flttyp,root,MPI_COMM_WORLD,mpierr)
  !
  ! Get the geometry for all the interior cells
  !
  allocate ( cell_geom(1:ncell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_geom",1,__LINE__,__FILE__,ierr,error_message)
  !
  do n = 1,ncell
    cell_geom(n) = nodes_to_geom( nr , nodes_of_cell_ptr(n+1) - &
                                       nodes_of_cell_ptr(n) )
  end do
  !
  allocate ( cell_order(1:ncell) , source=n_order , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_order",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Make sure that Lobatto nodes are requested if the grid contains triangle
  ! elements since Gauss nodes is not currently set up for triangles.
  !
 !if (any(cell_geom == Geom_Tria)) then
 !  write (error_message,11)
 !  call stop_recon(stop_mpi,pname,__LINE__,__FILE__,error_message)
 !end if
  !
  if (iop_solpts /= Legendre_Gauss_Lobatto .or. &
      iop_flxpts /= Legendre_Gauss_Lobatto) then
    if (any(cell_geom == Geom_Tria)) then
      write (error_message,12)
      call stop_recon(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (" CPU #",i0,": grid_size(:) = ",5(1x,i0),:,/,27x,*(5(1x,i0),/))
  !
  11 format (" Triangles are temporarily disabled due to significant", &
             " changes in the quadrilateral part of the code. These ", &
             " changes have yet to be applied to triangles.")
  12 format (" Quadrature points besides Gauss-Lobatto nodes are not", &
             " currently supported for a grid containing triangle cells.")
  !
end subroutine broadcast_solver_geom_arrays
!
!###############################################################################
!
subroutine create_serial_geom_arrays
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use geovar,    only : nr,ncell
  use geovar,    only : nodes_of_cell_ptr
  use geovar,    only : cell_geom,cell_order
  use ovar,      only : iop_solpts,iop_flxpts
  !
  !.. Local Scalars ..
  integer :: n,ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_serial_geom_arrays"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Get the geometry for all the interior cells
  !
  allocate ( cell_geom(1:ncell) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_geom",1,__LINE__,__FILE__,ierr,error_message)
  !
  do n = 1,ncell
    cell_geom(n) = nodes_to_geom( nr , nodes_of_cell_ptr(n+1) - &
                                       nodes_of_cell_ptr(n) )
  end do
  !
  allocate ( cell_order(1:ncell) , source=n_order , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_order",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Make sure that Lobatto nodes are requested if the grid contains triangle
  ! elements since Gauss nodes is not currently set up for triangles.
  !
 !if (any(cell_geom == Geom_Tria)) then
 !  write (error_message,11)
 !  call stop_recon(stop_mpi,pname,__LINE__,__FILE__,error_message)
 !end if
  !
  if (iop_solpts /= Legendre_Gauss_Lobatto .or. &
      iop_flxpts /= Legendre_Gauss_Lobatto) then
    if (any(cell_geom == Geom_Tria)) then
      write (error_message,12)
      call stop_recon(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  11 format (" Triangles are temporarily disabled due to significant", &
             " changes in the quadrilateral part of the code. These ", &
             " changes have yet to be applied to triangles.")
  12 format (" Quadrature points besides Gauss-Lobatto nodes are not", &
             " currently supported for a grid containing triangle cells.")
  !
end subroutine create_serial_geom_arrays
!
!###############################################################################
!
subroutine create_arrays_from_cgns_data
!!!
!!!.. Local Parameters ..
!!character(len=*), parameter :: pname = "create_arrays_from_cgns_data"
!!!
continue
!!!
!!call debug_timer(entering_procedure,pname)
!!!
!!npts  = data_sizes(3)
!!ncell = data_sizes(4)
!!!
!!! Extract the integers n1 and n2 from grid data
!!!
!!n1 = transfer( grid_data(1,1) , n1 )
!!n2 = transfer( grid_data(2,1) , n2 )
!!!
!!! If the minimum and maximum values of the coordinate in the third
!!! dimension are equal, reduce nr to 2 since this is a 2D grid
!!!
!!nr = n1
!!!
!!do n = 1,size(grid_data,dim=1)
!!  !
!!  minz = minval(grid_data(nr,2:1+n2))
!!  maxz = maxval(grid_data(nr,2:1+n2))
!!  !
!!  if (abs(maxz-minz) > eps6) exit
!!  !
!!  nr = nr-1
!!  !
!!end do
!!!
!!if (n2 /= npts) then
!!  write (error_message,101) n2,npts
!!  call stop_recon(stop_mpi,pname,__LINE__,__FILE__,error_message)
!!end if
!!!
!!! Allocate the xyz component of the grid derived type and use the grid_data
!!! array as the data source
!!!
!!allocate ( grid%xyz(1:nr,1:n2) , source=grid_data(1:nr,2:1+n2) , &
!!           stat=ierr , errmsg=error_message )
!!call alloc_error(pname,"grid%xyz",1,__LINE__,__FILE__,ierr,error_message)
!!!
!!! Deallocate grid_data now that we no longer need it
!!!
!!deallocate ( grid_data , stat=ierr , errmsg=error_message )
!!call alloc_error(pname,"grid_data",2,__LINE__,__FILE__,ierr,error_message)
!!!
!!! Use the section_data array to create the
!!! elem component of the grid derived type
!!!
!!nsections = data_sizes(8)
!!!
!!n = 0
!!n2 = 0
!!!
!!count_elem_loop: do is = 1,nsections
!!  !
!!  elem_end = section_data(n+3)
!!  sect_end = section_data(n+5)
!!  !
!!  n2 = max(n2,elem_end)
!!  !
!!  n = n + 5 + sect_end
!!  !
!!  if (n == size(section_data)) then
!!    exit count_elem_loop
!!  else if (n > size(section_data)) then
!!    write (error_message,102) is,nsections,n,size(section_data)
!!    call stop_recon(stop_mpi,pname,__LINE__,__FILE__,error_message)
!!  end if
!!  !
!!end do count_elem_loop
!!!
!!! Allocate the elem component of the grid derived type
!!!
!!allocate ( grid%elem(1:n2) , stat=ierr , errmsg=error_message )
!!call alloc_error(pname,"grid%elem",1,__LINE__,__FILE__,ierr,error_message)
!!!
!!n = 0
!!!
!!extract_elem_loop: do is = 1,nsections
!!  !
!!  sect_type   = section_data(n+1)
!!  elem_beg    = section_data(n+2)
!!  elem_end    = section_data(n+3)
!!  elem_nbndry = section_data(n+4)
!!  sect_end    = section_data(n+5)
!!  n = n + 5
!!  !
!!  if (sect_type == MIXED) then
!!    !
!!    n2 = n
!!    do i = elem_beg,elem_end
!!      n0 = n2 + 1
!!      elem_type = section_data(n0)
!!      np = npts_for_element_type(elem_type)
!!      n1 = n0 + 1
!!      n2 = n0 + np
!!      grid%elem(i) = create_element(elem_type,section_data(n1:n2))
!!    end do
!!    n = n2
!!    !
!!  else
!!    !
!!    elem_type = sect_type
!!    np = npts_for_element_type(elem_type)
!!    n2 = n
!!    do i = elem_beg,elem_end
!!      n1 = n2 + 1
!!      n2 = n2 + np
!!      grid%elem(i) = create_element(elem_type,section_data(n1:n2))
!!    end do
!!    n = n2
!!    !
!!  end if
!!  !
!!end do extract_elem_loop
!!!
!!! Deallocate section_data
!!!
!!deallocate ( section_data , stat=ierr , errmsg=error_message )
!!call alloc_error(pname,"section_data",2,__LINE__,__FILE__,ierr,error_message)
!!!
!!! Check that the cell dimension of all the interior cells are the same
!!!
!!!
!!! Use the boco_data array to create the bface array
!!!
!!nbocos = data_sizes(10)
!!!
!!n = 0
!!n2 = 0
!!!
!!count_boco_loop: do ibc = 1,nbocos
!!  !
!!  boco_end = boco_data(n+3)
!!  !
!!  n2 = n2 + boco_end
!!  n = n + 3 + boco_end
!!  !
!!  if (n == size(boco_data)) then
!!    exit count_boco_loop
!!  else if (n > size(boco_data)) then
!!    write (error_message,103) ibc,nbocos,n,size(boco_data)
!!    call stop_recon(stop_mpi,pname,__LINE__,__FILE__,error_message)
!!  end if
!!  !
!!end do count_boco_loop
!!!
!!! Allocate the bface array
!!!
!!allocate ( bface(1:nbfai,1:n2) , source=0 , stat=ierr , errmsg=error_message )
!!call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
!!!
!!n = 0
!!nf = 0
!!!
!!extract_boco_loop: do ibc = 1,nbocos
!!  !
!!  boco_type = boco_data(n+1)
!!  boco_loc  = boco_data(n+2)
!!  boco_end  = boco_data(n+3)
!!  !
!!  n = n + 3
!!  !
!!  do i = 1,boco_end
!!    !
!!    nf = nf + 1
!!    nc = boco_data(n+i)
!!    !
!!    host_cell = grid%elem(nc)%host_cell
!!    !
!!    bface( 1,nf) = boco_type               ! face boundary condition
!!    bface( 2,nf) = host_cell                      ! host cell of boundary face
!!    bface( 3,nf) = grid%elem(host_cell)%prop%geom ! host cell geometry
!!    bface( 5,nf) = ibc                     ! group for boundary face
!!    bface(10,nf) = ncell + nf              ! ghost cell of boundary face
!!    !
!!  end do
!!  !
!!  n = n + boco_end
!!  !
!!end do extract_boco_loop
!!!
!!! Deallocate boco_data since we dont need it anymore
!!!
!!deallocate ( boco_data , stat=ierr , errmsg=error_message )
!!call alloc_error(pname,"boco_data",2,__LINE__,__FILE__,ierr,error_message)
!!!
!!! 
!!!
!!do nf = 1,size(bface,dim=2)
!!  !
!!    !
!!    host_cell = bface(2,nf) ! host cell of boundary face
!!    host_geom = bface(3,nf) ! geometry of host cell
!!    !
!!    ! If the cell geometry is a node or edge, change ifa to a special value.
!!    ! Otherwise, continue on to finding the side of the host cell that is
!!    ! on the boundary.
!!    !
!!    if (all(host_geom /= [Geom_2D,Geom_3D])) then
!!      !
!!      if (host_geom == Geom_Node) then
!!        ifa = -1
!!      else if (host_geom == Geom_Edge) then
!!        ifa = -2
!!      else
!!        ifa = -10
!!      end if
!!      !
!!    else
!!      !
!!      ! Find which face of the host cell is on the boundary
!!      !
!!      ! total number of host cell nodes
!!      ncnods = size(grid%elem(host_cell)%nodes)
!!      !
!!      ! Create a local copy of the interior nodes
!!      !
!!      int_nodes(1:ncnods) = grid%elem(host_cell)%nodes
!!      !
!!      marker = 0 ! Initialize the marker array to zero
!!      !
!!      ! Mark each host cell node that is on the boundary
!!      !
!!      do n = 1,size(grid%elem(bnd_cell)%nodes)
!!        where (int_nodes(1:ncnods) == grid%elem(bnd_cell)%nodes(n)) marker = 1
!!      end do
!!      !
!!      ! Loop over the number of faces for this geometry type and find which
!!      ! face corresponds to the marked nodes
!!      !
!!      ifa = -99 ! Initialize ifa to make sure that a side is actually found
!!      !
!!      do host_side = 1,geom_faces(host_geom)
!!        !
!!        ! Get the number of nodes defining this face
!!        !
!!        nfnods = num_face_nodes(host_side,host_geom)
!!        !
!!        ! Extract the cell nodes that define the current face
!!        ! NOTE: The +1 at the end is needed since the values stored in
!!        !       side_nodes are normally used as offsets, meaning that their
!!        !       values are in the range [0:ncnods-1]. We want the values
!!        !       in the range of [1:ncnods] hence the +1.
!!        !
!!        face_nodes(1:nfnods) = side_nodes(1:nfnods,host_side,host_geom) + 1
!!        !
!!        ! If all of these face nodes have been marked, then this is the
!!        ! side of the host cell on the boundary. Save the value into ifa
!!        ! and exit the face loop.
!!        !
!!        if ( sum( marker(face_nodes(1:nfnods)) ) == nfnods ) then
!!          ifa = host_side
!!          exit
!!        end if
!!        !
!!      end do
!!      !
!!    end if
!!    !
!!    ! Check for an error finding which face of the host cell is on the
!!    ! boundary, and if there is no error, store the value in bface
!!    !
!!    if (ifa < 0) then
!!      !
!!      if (.not. error_found) then
!!        if (mypnum == 0) write (iout,4)
!!        error_found = true
!!      end if
!!      !
!!      ! Write out the specific error type corresponding
!!      ! to the negative value of ifa
!!      !
!!      if (mypnum == 0) then
!!        if (ifa == -1) then
!!          write (iout,5) "A 0-D 'Point'-type",host_cell,bnd_cell
!!        else if (ifa == -2) then
!!          write (iout,5) "A 1-D 'Edge'-type",host_cell,bnd_cell
!!        else if (ifa == -10) then
!!          write (iout,5) "An unknown type of",host_cell,bnd_cell
!!        else if (ifa == -99) then
!!          write (iout,6) host_cell,bnd_cell
!!        end if
!!      end if
!!      !
!!    else
!!      !
!!      ! Store the host cell face that is on the boundary
!!      !
!!      bface(4,nf) = ifa
!!      !
!!    end if
!!    !
!!  end do current_cell_grp
!!  !
!!end do boundary_cell_grps

!!  !

!!!
!!! - Take VERTEX located ZoneBC data and create FACECENTERED BC Data
!!! 
!!call debug_timer(leaving_procedure,pname)
!!!
end subroutine create_arrays_from_cgns_data
!
!###############################################################################
!
!!!subroutine fix_boco_lists(nvertex,ncell,s,bc)
!!!  !
!!!  !.. Formal Arguments ..
!!!  integer(CST),       intent(in) :: nvertex
!!!  integer(CST),       intent(in) :: ncell
!!!  type(cg_section_t), intent(in) :: s(:)
!!!  type(cg_bc_t),   intent(inout) :: bc(:)
!!!  !
!!!  !.. Local Scalars ..
!!!  integer :: is,ibc,nd
!!!  !
!!!  !.. Local Parameters ..
!!!  character(len=*), parameter :: pname = "fix_boco_lists"
!!!  !
!!!  type :: bc_elem_t
!!!    logical(lk)  :: done = fals
!!!    integer(CET) :: face_bc = 0
!!!    integer, allocatable :: nodes(:)
!!!  end type bc_elem_t
!!!  !
!!!  type(bc_elem_t), allocatable :: bc_elem(:)
!!!  !
!!!  type :: bc_vertex_t
!!!    integer :: n = 0
!!!    integer(CET), allocatable :: bc(:)
!!!  end type bc_vertex_t
!!!  !
!!!  type(bc_vertex_t), allocatable :: bc_vertex(:)
!!!  !
!!!  integer(CET), parameter :: elem_loc(0:3) = [VERTEX, &
!!!                                              EDGECENTER, &
!!!                                              FACECENTER, &
!!!                                              CELLCENTER]
!!!  !
!!!continue
!!!  !
!!!  call debug_timer(entering_procedure,pname)
!!!  !
!!!  ! Check for boundary conditions that have a location of VERTEX
!!!  !
!!!  do ibc = 1,size(bc)
!!!    !
!!!    if (bc(ibc)%location == VERTEX) then
!!!      !
!!!      ! Check to see if the boundary condition ptset_type is using either of
!!!      ! the depreciated PointSetType values of ELEMENTRANGE or ELEMENTLIST.
!!!      ! If it is, switch ptset_type to the respective POINTRANGE OR POINTLIST
!!!      ! and switch location to either FACECENTER OR EDGECENTER depending on
!!!      ! the element types.
!!!      !
!!!      nd = 0
!!!      if (bc(ibc)%ptset_type == ELEMENTRANGE) then
!!!        bc%ptset_type = POINTRANGE
!!!        nd = get_elem_dim_for_boco_pointrange(bc(ibc)%pnts(1), &
!!!                                              bc(ibc)%pnts(2),s)
!!!      else if (bc(ibc)%ptset_type == ELEMENTLIST) then
!!!        bc%ptset_type = POINTLIST
!!!        nd = get_elem_dim_for_boco_pointlist(bc(ibc)%pnts,s)
!!!      end if
!!!      !
!!!      bc%location = elem_loc(nd)
!!!      !
!!!    end if
!!!    !
!!!  end do
!!!  !
!!!  n1 = maxval( s(:)%end )
!!!  !
!!!  allocate ( bc_elem(1:n1) , stat=ierr , errmsg=error_message )
!!!  !
!!!  ! Mark nodes on boundaries
!!!  !
!!!  allocate ( bc_vertex(1:nvertex) , stat=ierr , errmsg=error_message )
!!!  !
!!!  do ibc = 1,size(bc)
!!!    !
!!!    if (bc(ibc)%ptset_type == PointRange) then
!!!      do n = bc(ibc)%pnts(1),bc(ibc)%pnts(2)
!!!        bc_vertex(n)%n = bc_vertex(n)%n + 1
!!!      end do
!!!    else if (bc(ibc)%ptset_type == PointList) then
!!!      do i = 1,size(bc(ibc)%pnts)
!!!        n = bc(ibc)%pnts(i)
!!!        bc_vertex(n)%n = bc_vertex(n)%n + 1
!!!      end do
!!!    else if (bc(ibc)%ptset_type == ElementRange) then
!!!      do n = bc(ibc)%pnts(1),bc(ibc)%pnts(2)
!!!        bc_elem(n)%face_bc = bc(ibc)%boco_type
!!!        bc_elem(n)%done = true
!!!      end do
!!!    else if (bc(ibc)%ptset_type == ElementList) then
!!!      do i = 1,size(bc(ibc)%pnts)
!!!        n = bc(ibc)%pnts(i)
!!!        bc_elem(n)%face_bc = bc(ibc)%boco_type
!!!        bc_elem(n)%done = true
!!!      end do
!!!    end if
!!!    !
!!!  end do
!!!  !
!!!  do n = 1,size(bc_vertex)
!!!    if (bc_vertex(n)%n > 0) then
!!!      allocate ( bc_vertex(n)%bc(1:bc_vertex(n)%n) , source=0_CET , &
!!!                 stat=ierr , errmsg=error_message )
!!!      bc_vertex(n)%n = 0
!!!    end if
!!!  end do
!!!  !
!!!  do ibc = 1,size(bc)
!!!    !
!!!    if (bc(ibc)%ptset_type == PointRange) then
!!!      do i = bc(ibc)%pnts(1),bc(ibc)%pnts(2)
!!!        n = bc_vertex(n)%n + 1
!!!        bc_vertex(i)%n = n
!!!        bc_vertex(i)%bc(n) = bc(ibc)%boco_type
!!!      end do
!!!    else if (bc(ibc)%ptset_type == PointList) then
!!!      do j = 1,size(bc(ibc)%pnts)
!!!        i = bc(ibc)%pnts(j)
!!!        n = bc_vertex(i)%n + 1
!!!        bc_vertex(i)%n = n
!!!        bc_vertex(i)%bc(n) = bc(ibc)%boco_type
!!!      end do
!!!    end if
!!!    !
!!!  end do
!!!  !
!!!  ! Initial loop through all elements which assigns a boundary condition
!!!  ! to any elements that have one or more nodes that has only a single
!!!  ! boundary condition
!!!  !
!!!  do n = 1,size(bc_elem)
!!!    !
!!!    if (bc_elem(n)%done) cycle
!!!    !
!!!    if (sum(bc_vertex(bc_elem(n)%nodes(:))%n) == 0) then
!!!      bc_elem(n)%done = true
!!!      bc_elem(n)%face_bc = 0
!!!    end if
!!!    !
!!!  end do
!!!       
!!!
!!!
!!!
!!!
!!!  
!!!  if (any(bc%ptset_type == [ElementRange,ElementList])) then
!!!    if (
!!!  if (any(bc%ptset_type == [PointRange,ElementRange])) then
!!!    !
!!!    np = pnts(2) - pnts(1) + 1
!!!    !
!!!    bc%pnts(:) = intseq(pnts(1),pnts(2))
!!!    !
!!!  else
!!!    !
!!!    np = size( pnts , kind=CST )
!!!    !
!!!    allocate ( bc%pnts(1:np) , stat=ierr , errmsg=error_message )
!!!    array_name = get_array_name("bc",ibc_g,"pnts")
!!!    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
!!!    !
!!!    bc%pnts(:) = pnts(:)
!!!    !
!!!  end if
!!!  !
!!!  call debug_timer(leaving_procedure,pname)
!!!  !
!!!end subroutine fix_boco_lists
!
!###############################################################################
!
pure function get_elem_dim_for_boco_pointlist(pnts,s) result(return_value)
  !
  !.. Formal Arguments ..
  integer(CST),       intent(in) :: pnts(:)
  type(cg_section_t), intent(in) :: s(:)
  !
  !.. Function Result ..
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: np,is,ndim,i
  integer(CET) :: elem_type
  integer(CST) :: ibeg,iend
  !
  !.. Local Arrays ..
  integer     :: edim(1:size(pnts))
  logical(lk) :: msk(1:size(pnts))
  !
continue
  !
  ! Initialize return_value to 0 (VERTEX dimension)
  !
  return_value = 0
  !
  np = size(pnts)
  !
  ! Loop through the sections and find all
  ! the elements that belong to this section
  !
  msk(:) = fals
  edim(:) = 0
  !
  do is = 1,size(s)
    !
    ibeg = s(is)%start
    iend = s(is)%end
    !
    msk(:) = ( (pnts(:) >= ibeg) .and. (pnts(:) <= iend) )
    !
    if (s(is)%type /= MIXED) then
      !
      ndim = geom_dimen( geom_for_element_type(s(is)%type) )
      where (msk(:)) edim(:) = ndim
      !
    else
      !
     !do i = 1,np
     !  if (msk(i)) then
     !    elem_type = s(is)%elements(pnts(i))
     !    edim(i) = geom_dimen( geom_for_element_type(elem_type) )
     !  end if
     !end do
      !
    end if
    !
    msk(:) = fals
    !
  end do
  !
  if (minval(edim) == maxval(edim)) then
    return_value = minval(edim)
  end if
  !
end function get_elem_dim_for_boco_pointlist
!
!###############################################################################
!
pure function get_elem_dim_for_boco_pointrange(pbeg,pend,s) result(return_value)
  !
  !.. Formal Arguments ..
  integer(CST),       intent(in) :: pbeg
  integer(CST),       intent(in) :: pend
  type(cg_section_t), intent(in) :: s(:)
  !
  !.. Function Result ..
  integer :: return_value
  !
  !.. Local Arrays ..
  integer(CST) :: pnts(1:pend-pbeg+1)
  !
continue
  !
  pnts(:) = intseq(pbeg,pend)
  !
  return_value = get_elem_dim_for_boco_pointlist(pnts,s)
  !
end function get_elem_dim_for_boco_pointrange
!
!###############################################################################
!
pure function create_element_i4(elem_type,pts) result(return_value)
  !
  integer, parameter :: lp = i4
  !
  !.. Formal Arguments ..
  integer(CET),              intent(in) :: elem_type
  integer(lp), dimension(:), intent(in) :: pts
  !
  !.. Function Result ..
  type(elem_t) :: return_value
  !
  !.. Local Scalars ..
  integer :: ierr,nn,np
  !
continue
  !
  return_value%prop = element_properties(elem_type)
  !
  np = size(pts)
  nn = geom_nodes( return_value%prop%geom )
  !
  ! Using F2003 auto-reallocation
  !
  return_value%pts = int(pts(1:np),kind=kind(return_value%pts))
  return_value%nodes = int(pts(1:nn),kind=kind(return_value%nodes))
 !!
 !! Using traditional allocate statements
 !!
 !allocate ( return_value%pts(1:np) , &
 !           source=int(pts(1:np),kind=kind(return_value%pts)) , &
 !           stat=ierr )
 !!
 !allocate ( return_value%nodes(1:nn) , &
 !           source=int(pts(1:nn),kind=kind(return_value%nodes)) , &
 !           stat=ierr )
  !
  return_value = reorder_cgns_pts(return_value)
  !
end function create_element_i4
!
!###############################################################################
!
pure function create_element_i8(elem_type,pts) result(return_value)
  !
  integer, parameter :: lp = i8
  !
  !.. Formal Arguments ..
  integer(CET),              intent(in) :: elem_type
  integer(lp), dimension(:), intent(in) :: pts
  !
  !.. Function Result ..
  type(elem_t) :: return_value
  !
  !.. Local Scalars ..
  integer :: ierr,nn,np
  !
continue
  !
  return_value%prop = element_properties(elem_type)
  !
  np = size(pts)
  nn = geom_nodes( return_value%prop%geom )
  !
  ! Using F2003 auto-reallocation
  !
  return_value%pts = int(pts(1:np),kind=kind(return_value%pts))
  return_value%nodes = int(pts(1:nn),kind=kind(return_value%nodes))
 !!
 !! Using traditional allocate statements
 !!
 !allocate ( return_value%pts(1:np) , &
 !           source=int(pts(1:np),kind=kind(return_value%pts)) , &
 !           stat=ierr )
 !!
 !allocate ( return_value%nodes(1:nn) , &
 !           source=int(pts(1:nn),kind=kind(return_value%nodes)) , &
 !           stat=ierr )
  !
  return_value = reorder_cgns_pts(return_value)
  !
end function create_element_i8
!
!###############################################################################
!
pure function npts_for_element_type(elem_type) result(return_value)
  !
  !.. Formal Arguments ..
  integer(CET), intent(in) :: elem_type
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
pure function geom_for_element_type(elem_type) result(return_value)
  !
  !.. Formal Arguments ..
  integer(CET), intent(in) :: elem_type
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
  return_value = prop%geom
  !
end function geom_for_element_type
!
!###############################################################################
!
pure function element_properties(elem_type) result(return_value)
  !
  !.. Formal Arguments ..
  integer(CET), intent(in) :: elem_type
  !
  !.. Function Result ..
  type(prop_t) :: return_value
  !
continue
  !
  select case (elem_type)
    case (NODE)
      return_value = prop_t(npts=1, &
                            order=1, &
                            geom=Geom_Node, &
                            elem_type=int(elem_type))
    case (BAR_2)
      return_value = prop_t(npts=2, &
                            order=1, &
                            geom=Geom_Edge, &
                            elem_type=int(elem_type))
    case (BAR_3)
      return_value = prop_t(npts=3, &
                            order=2, &
                            geom=Geom_Edge, &
                            elem_type=int(elem_type))
    case (TRI_3)
      return_value = prop_t(npts=3, &
                            order=1, &
                            geom=Geom_Tria, &
                            elem_type=int(elem_type))
    case (TRI_6)
      return_value = prop_t(npts=6, &
                            order=2, &
                            geom=Geom_Tria, &
                            elem_type=int(elem_type))
    case (QUAD_4)
      return_value = prop_t(npts=4, &
                            order=1, &
                            geom=Geom_Quad, &
                            elem_type=int(elem_type))
    case (QUAD_8)
      return_value = prop_t(npts=8, &
                            order=2, &
                            geom=Geom_Quad, &
                            elem_type=int(elem_type), &
                            family=Edge_Defined)
    case (QUAD_9)
      return_value = prop_t(npts=9, &
                            order=2, &
                            geom=Geom_Quad, &
                            elem_type=int(elem_type))
    case (TETRA_4)
      return_value = prop_t(npts=4, &
                            order=1, &
                            geom=Geom_Tetr, &
                            elem_type=int(elem_type))
    case (TETRA_10)
      return_value = prop_t(npts=10, &
                            order=2, &
                            geom=Geom_Tetr, &
                            elem_type=int(elem_type))
    case (PYRA_5)
      return_value = prop_t(npts=5, &
                            order=1, &
                            geom=Geom_Pyra, &
                            elem_type=int(elem_type))
    case (PYRA_13)
      return_value = prop_t(npts=13, &
                            order=2, &
                            geom=Geom_Pyra, &
                            elem_type=int(elem_type), &
                            family=Edge_Defined)
    case (PYRA_14)
      return_value = prop_t(npts=14, &
                            order=2, &
                            geom=Geom_Pyra, &
                            elem_type=int(elem_type))
    case (PENTA_6)
      return_value = prop_t(npts=6, &
                            order=1, &
                            geom=Geom_Pris, &
                            elem_type=int(elem_type))
    case (PENTA_15)
      return_value = prop_t(npts=15, &
                            order=2, &
                            geom=Geom_Pris, &
                            elem_type=int(elem_type), &
                            family=Edge_Defined)
    case (PENTA_18)
      return_value = prop_t(npts=18, &
                            order=2, &
                            geom=Geom_Pris, &
                            elem_type=int(elem_type))
    case (HEXA_8)
      return_value = prop_t(npts=8, &
                            order=1, &
                            geom=Geom_Hexa, &
                            elem_type=int(elem_type))
    case (HEXA_20)
      return_value = prop_t(npts=20, &
                            order=2, &
                            geom=Geom_Hexa, &
                            elem_type=int(elem_type), &
                            family=Edge_Defined)
    case (HEXA_27)
      return_value = prop_t(npts=27, &
                            order=2, &
                            geom=Geom_Hexa, &
                            elem_type=int(elem_type))
   !case (MIXED)   ! MIXED TYPE SHOULD BE SEPARATED BY NOW
   !case (NGON_n)  ! UNSUPPORTED TYPE
   !case (NFACE_n) ! UNSUPPORTED TYPE
    case (BAR_4)
      return_value = prop_t(npts=4, &
                            order=3, &
                            geom=Geom_Edge, &
                            elem_type=int(elem_type))
    case (TRI_9)
      return_value = prop_t(npts=9, &
                            order=3, &
                            geom=Geom_Tria, &
                            elem_type=int(elem_type), &
                            family=Edge_Defined)
    case (TRI_10)
      return_value = prop_t(npts=10, &
                            order=3, &
                            geom=Geom_Tria, &
                            elem_type=int(elem_type))
    case (QUAD_12)
      return_value = prop_t(npts=12, &
                            order=3, &
                            geom=Geom_Quad, &
                            elem_type=int(elem_type), &
                            family=Edge_Defined)
    case (QUAD_16)
      return_value = prop_t(npts=16, &
                            order=3, &
                            geom=Geom_Quad, &
                            elem_type=int(elem_type))
    case (TETRA_16)
      return_value = prop_t(npts=16, &
                            order=3, &
                            geom=Geom_Tetr, &
                            elem_type=int(elem_type), &
                            family=Edge_Defined)
    case (TETRA_20)
      return_value = prop_t(npts=20, &
                            order=3, &
                            geom=Geom_Tetr, &
                            elem_type=int(elem_type))
    case (PYRA_21)
      return_value = prop_t(npts=21, &
                            order=3, &
                            geom=Geom_Pyra, &
                            elem_type=int(elem_type), &
                            family=Edge_Defined)
    case (PYRA_29)
      return_value = prop_t(npts=29, &
                            order=3, &
                            geom=Geom_Pyra, &
                            elem_type=int(elem_type), &
                            family=Face_Defined)
    case (PYRA_30)
      return_value = prop_t(npts=30, &
                            order=3, &
                            geom=Geom_Pyra, &
                            elem_type=int(elem_type))
    case (PENTA_24)
      return_value = prop_t(npts=24, &
                            order=3, &
                            geom=Geom_Pris, &
                            elem_type=int(elem_type), &
                            family=Edge_Defined)
    case (PENTA_38)
      return_value = prop_t(npts=38, &
                            order=3, &
                            geom=Geom_Pris, &
                            elem_type=int(elem_type), &
                            family=Face_Defined)
    case (PENTA_40)
      return_value = prop_t(npts=40, &
                            order=3, &
                            geom=Geom_Pris, &
                            elem_type=int(elem_type))
    case (HEXA_32)
      return_value = prop_t(npts=32, &
                            order=3, &
                            geom=Geom_Hexa, &
                            elem_type=int(elem_type), &
                            family=Edge_Defined)
    case (HEXA_56)
      return_value = prop_t(npts=56, &
                            order=3, &
                            geom=Geom_Hexa, &
                            elem_type=int(elem_type), &
                            family=Face_Defined)
    case (HEXA_64)
      return_value = prop_t(npts=64, &
                            order=3, &
                            geom=Geom_Hexa, &
                            elem_type=int(elem_type))
#ifdef CGNS_3_3
    case (BAR_5)
      return_value = prop_t(npts=5, &
                            order=4, &
                            geom=Geom_Edge, &
                            elem_type=int(elem_type))
    case (TRI_12)
      return_value = prop_t(npts=12, &
                            order=4, &
                            geom=Geom_Tria, &
                            elem_type=int(elem_type), &
                            family=Edge_Defined)
    case (TRI_15)
      return_value = prop_t(npts=15, &
                            order=4, &
                            geom=Geom_Tria, &
                            elem_type=int(elem_type))
    case (QUAD_P4_16)
      return_value = prop_t(npts=16, &
                            order=4, &
                            geom=Geom_Quad, &
                            elem_type=int(elem_type), &
                            family=Edge_Defined)
    case (QUAD_25)
      return_value = prop_t(npts=25, &
                            order=4, &
                            geom=Geom_Quad, &
                            elem_type=int(elem_type))
    case (TETRA_22)
      return_value = prop_t(npts=22, &
                            order=4, &
                            geom=Geom_Tetr, &
                            elem_type=int(elem_type), &
                            family=Edge_Defined)
    case (TETRA_34)
      return_value = prop_t(npts=34, &
                            order=4, &
                            geom=Geom_Tetr, &
                            elem_type=int(elem_type), &
                            family=Face_Defined)
    case (TETRA_35)
      return_value = prop_t(npts=35, &
                            order=4, &
                            geom=Geom_Tetr, &
                            elem_type=int(elem_type))
    case (PYRA_P4_29)
      return_value = prop_t(npts=29, &
                            order=4, &
                            geom=Geom_Pyra, &
                            elem_type=int(elem_type), &
                            family=Edge_Defined)
    case (PYRA_50)
      return_value = prop_t(npts=50, &
                            order=4, &
                            geom=Geom_Pyra, &
                            elem_type=int(elem_type), &
                            family=Face_Defined)
    case (PYRA_55)
      return_value = prop_t(npts=55, &
                            order=4, &
                            geom=Geom_Pyra, &
                            elem_type=int(elem_type))
    case (PENTA_33)
      return_value = prop_t(npts=33, &
                            order=4, &
                            geom=Geom_Pris, &
                            elem_type=int(elem_type), &
                            family=Edge_Defined)
    case (PENTA_66)
      return_value = prop_t(npts=66, &
                            order=4, &
                            geom=Geom_Pris, &
                            elem_type=int(elem_type), &
                            family=Face_Defined)
    case (PENTA_75)
      return_value = prop_t(npts=75, &
                            order=4, &
                            geom=Geom_Pris, &
                            elem_type=int(elem_type))
    case (HEXA_44)
      return_value = prop_t(npts=44, &
                            order=4, &
                            geom=Geom_Hexa, &
                            elem_type=int(elem_type), &
                            family=Edge_Defined)
    case (HEXA_98)
      return_value = prop_t(npts=98, &
                            order=4, &
                            geom=Geom_Hexa, &
                            elem_type=int(elem_type), &
                            family=Face_Defined)
    case (HEXA_125)
      return_value = prop_t(npts=125, &
                            order=4, &
                            geom=Geom_Hexa, &
                            elem_type=int(elem_type))
#endif
    case default
      return_value = prop_t(npts=-1, &
                            order=-1, &
                            geom=Geom_Unknown, &
                            elem_type=int(elem_type))
  end select
  !
end function element_properties
!
!###############################################################################
!
elemental function reorder_cgns_pts(elem) result(return_value)
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
    ! This is a complete element with internal node points
    !
  select case (elem%prop%geom)
    case (Geom_Edge)
      return_value%pts = reorder_cgns_edge_pts(elem%pts)
    case (Geom_Tria)
      return_value%pts = reorder_cgns_tria_pts(elem%pts,elem%prop%elem_type)
    case (Geom_Quad)
      return_value%pts = reorder_cgns_quad_pts(elem%pts,elem%prop%elem_type)
    case (Geom_Tetr)
      return_value%pts = reorder_cgns_tetr_pts(elem%pts,elem%prop%elem_type)
    case (Geom_Pyra)
      return_value%pts = reorder_cgns_pyra_pts(elem%pts,elem%prop%elem_type)
    case (Geom_Pris)
      return_value%pts = reorder_cgns_pris_pts(elem%pts,elem%prop%elem_type)
    case (Geom_Hexa)
      return_value%pts = reorder_cgns_hexa_pts(elem%pts,elem%prop%elem_type)
  end select
  !
end function reorder_cgns_pts
!
!###############################################################################
!
pure function reorder_cgns_edge_pts(pts) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Scalars ..
  integer :: n,i
  !
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
end function reorder_cgns_edge_pts
!
!###############################################################################
!
pure function reorder_cgns_tria_pts(pts,elem_type) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  integer, intent(in) :: elem_type
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: tria3(3) = [1,2,3]
  !
  integer, parameter :: tria6(6) = [1,4,2,6,5,3]
  !
  integer, parameter :: tria9(9) = [1,4,5,2,9,6,8,7,3]
  !
  integer, parameter :: tria10(10) = [1,4,5,2,9,10,6,8,7,3]
  !
  integer, parameter :: tria12(12) = [1,4,5,6,2,12,7,11,8,10,9,3]
  !
  integer, parameter :: tria15(15) = [1,4,5,6,2,12,13,14,7,11,15,8,10,9,3]
  !
continue
  !
  select case (elem_type)
    case (TRI_3)
      return_value = pts( tria3(:) )
    case (TRI_6)
      return_value = pts( tria6(:) )
    case (TRI_9)
      return_value = pts( tria9(:) )
    case (TRI_10)
      return_value = pts( tria10(:) )
    case (TRI_12)
      return_value = pts( tria12(:) )
    case (TRI_15)
      return_value = pts( tria15(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_cgns_tria_pts
!
!###############################################################################
!
pure function reorder_cgns_quad_pts(pts,elem_type) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  integer, intent(in) :: elem_type
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: quad4(4) = [1,2,4,3]
  !
  integer, parameter :: quad8(8) = [1,5,2,8,6,4,7,3]
  !
  integer, parameter :: quad9(9) = [1,5,2,8,9,6,4,7,3]
  !
  integer, parameter :: quad12(12) = [1,5,6,2,12,7,11,8,4,10,9,3]
  !
  integer, parameter :: quad3_16(16) = [ 1, 5, 6, 2, &
                                       12,13,14, 7, &
                                       11,16,15, 8, &
                                        4,10, 9, 3]
  !
  integer, parameter :: quad4_16(16) = [ 1, 5, 6, 7, &
                                        2,16, 8,15, &
                                        9,14,10, 4, &
                                       13,12,11, 3]
  !
  integer, parameter :: quad25(25) = [ 1, 5, 6, 7, 2, &
                                     16,17,18,19, 8, &
                                     15,24,25,20, 9, &
                                     14,23,22,21,10, &
                                      4,13,12,11, 3]
  !
continue
  !
  select case (elem_type)
    case (QUAD_4)
      return_value = pts( quad4(:) )
    case (QUAD_8)
      return_value = pts( quad8(:) )
    case (QUAD_9)
      return_value = pts( quad9(:) )
    case (QUAD_12)
      return_value = pts( quad12(:) )
    case (QUAD_16)
      return_value = pts( quad3_16(:) )
    case (QUAD_P4_16)
      return_value = pts( quad4_16(:) )
    case (QUAD_25)
      return_value = pts( quad25(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_cgns_quad_pts
!
!###############################################################################
!
pure function reorder_cgns_tetr_pts(pts,elem_type) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  integer, intent(in) :: elem_type
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: tetr4(4) = [1,4,2,3]
  !
  integer, parameter :: tetr10(10) = [1,8,4,5,9,2,7,10,6,3]
  !
  integer, parameter :: tetr16(16) = [ 1,11,12, 4, 5,14, 6,13, &
                                      2,10,16, 7, 9,15, 8, 3]
  !
  integer, parameter :: tetr20(20) = [ 1,11,12, 4, 5,18,14, 6,13, 2, &
                                     10,20,16,17,19, 7, 9,15, 8, 3]
  !
  integer, parameter :: tetr22(22) = [ 1,14,15,16, 4, 5,19, 6,18, 7,17, &
                                      2,13,22, 8,12,21, 9,11,20,10, 3]
  !
  integer, parameter :: tetr34(34) = [ 1,14,15,16, 4, 5,26, &
                                     28,19, 6,27,18, 7,17, &
                                      2,13,33,34,22,23,    &
                                     31,24,29, 8,12,32,21, &
                                     25,30, 9,11,20,10, 3]
  !
  integer, parameter :: tetr35(35) = [ 1,14,15,16, 4, 5,26, &
                                     28,19, 6,27,18, 7,17, &
                                      2,13,33,34,22,23,35, &
                                     31,24,29, 8,12,32,21, &
                                     25,30, 9,11,20,10, 3]
  !
continue
  !
  select case (elem_type)
    case (TETRA_4)
      return_value = pts( tetr4(:) )
    case (TETRA_10)
      return_value = pts( tetr10(:) )
    case (TETRA_16)
      return_value = pts( tetr16(:) )
    case (TETRA_20)
      return_value = pts( tetr20(:) )
    case (TETRA_22)
      return_value = pts( tetr22(:) )
    case (TETRA_34)
      return_value = pts( tetr34(:) )
    case (TETRA_35)
      return_value = pts( tetr35(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_cgns_tetr_pts
!
!###############################################################################
!
pure function reorder_cgns_pyra_pts(pts,elem_type) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  integer, intent(in) :: elem_type
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: pyra5(5) = [1,2,4,3,5]
  !
  integer, parameter :: pyra13(13) = [1,6,2,9,7,4,8,3,8,10,13,12,5]
  !
  integer, parameter :: pyra14(14) = [1,6,2,9,14,7,4,8,3,10,11,13,12,5]
  !
  integer, parameter :: pyra21(21) = [ 1, 6, 7, 2,13, 8,12, &
                                      9, 4,11,10, 3,14,16, &
                                     20,18,15,17,21,19, 5]
  !
  integer, parameter :: pyra3_29(29) = [ 1, 6, 7, 2,13,22,23, 8,12,25, &
                                       24, 9, 4,11,10, 3,14,26,16,29, &
                                       27,20,28,18,15,17,21,19, 5]
  !
  integer, parameter :: pyra4_29(29) = [ 1, 6, 7, 8, 2,17, 9,16,10,15, &
                                       11, 4,14,13,12, 3,18,21,27,24, &
                                       19,22,28,25,20,23,29,26, 5]
  !
  integer, parameter :: pyra30(30) = [ 1, 6, 7, 2,13,22,23, 8,12,25, &
                                     24, 9, 4,11,10, 3,14,26,16,29, &
                                     30,27,20,28,18,15,17,21,19, 5]
  !
  integer, parameter :: pyra50(50) = [ 1, 6, 7, 8, 2,17,30,31,32, 9, &
                                     16,37,38,33,10,15,36,35,34,11, &
                                      4,14,13,12, 3,18,39,40,21,49, &
                                     42,48,43,27,46,45,24,19,41,22, &
                                     50,44,28,47,25,20,23,29,26, 5]
  !
  integer, parameter :: pyra55(55) = [ 1, 6, 7, 8, 2,17,30,31,32, 9,16, &
                                     37,38,33,10,15,36,35,34,11, 4,14, &
                                     13,12, 3,18,39,40,21,49,51,52,42, &
                                     48,54,53,43,27,46,45,24,19,41,22, &
                                     50,55,44,28,47,25,20,23,29,26, 5]
  !
continue
  !
  select case (elem_type)
    case (PYRA_5)
      return_value = pts( pyra5(:) )
    case (PYRA_13)
      return_value = pts( pyra13(:) )
    case (PYRA_14)
      return_value = pts( pyra14(:) )
    case (PYRA_21)
      return_value = pts( pyra21(:) )
    case (PYRA_29)
      return_value = pts( pyra3_29(:) )
    case (PYRA_P4_29)
      return_value = pts( pyra4_29(:) )
    case (PYRA_30)
      return_value = pts( pyra30(:) )
    case (PYRA_50)
      return_value = pts( pyra50(:) )
    case (PYRA_55)
      return_value = pts( pyra55(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_cgns_pyra_pts
!
!###############################################################################
!
pure function reorder_cgns_pris_pts(pts,elem_type) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  integer, intent(in) :: elem_type
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: pris6(6) = [1,2,3,4,5,6]
  !
  integer, parameter :: pris15(15) = [ 1, 7, 2, 9, 8, 3, &
                                     10,   11,      12, &
                                      4,13, 5,15,14, 6]
  !
  integer, parameter :: pris18(18) = [ 1, 7, 2, 9, 8, 3, &
                                     10,16,11,18,17,12, &
                                      4,13, 5,15,14, 6]
  !
  integer, parameter :: pris24(24) = [ 1, 7, 8, 2,12, 9,11,10, 3, &
                                     13,      15,            17, &
                                     14,      16,            18, &
                                      4,19,20, 5,24,21,23,22, 6]
  !
  integer, parameter :: pris38(38) = [ 1, 7, 8, 2,12,25, 9,11,10, 3, &
                                     13,26,27,15,35,   30,34,31,17, &
                                     14,29,28,16,36,   33,37,32,18, &
                                      4,19,20, 5,24,38,21,23,22, 6]
  !
  integer, parameter :: pris40(40) = [ 1, 7, 8, 2,12,25, 9,11,10, 3, &
                                     13,26,27,15,35,39,30,34,31,17, &
                                     14,29,28,16,36,40,33,37,32,18, &
                                      4,19,20, 5,24,38,21,23,22, 6]
  !
  integer, parameter :: pris33(33) = &
    [ 1, 7, 8, 9, 2,15,10,14,11,13,12, 3, &
     16,         19,                  22, &
     17,         20,                  23, &
     18,         21,                  24, &
      4,25,26,27, 5,33,28,32,29,31,30, 6]
  !
  integer, parameter :: pris66(66) = &
    [ 1, 7, 8, 9, 2,15,34,35,10,14,36,11,13,12, 3, &
     16,37,38,39,19,57,      46,56,   47,55,48,22, &
     17,44,45,40,20,58,      53,63,   54,62,49,23, &
     18,43,42,41,21,59,      52,60,   51,61,50,24, &
      4,25,26,27, 5,33,64,65,28,32,66,29,31,30, 6]
  !
  integer, parameter :: pris75(75) = &
    [ 1, 7, 8, 9, 2,15,34,35,10,14,36,11,13,12, 3, &
     16,37,38,39,19,57,67,68,46,56,69,47,55,48,22, &
     17,44,45,40,20,58,70,71,53,63,72,54,62,49,23, &
     18,43,42,41,21,59,73,74,52,60,75,51,61,50,24, &
      4,25,26,27, 5,33,64,65,28,32,66,29,31,30, 6]
  !
continue
  !
  select case (elem_type)
    case (PENTA_6)
      return_value = pts( pris6(:) )
    case (PENTA_15)
      return_value = pts( pris15(:) )
    case (PENTA_18)
      return_value = pts( pris18(:) )
    case (PENTA_24)
      return_value = pts( pris24(:) )
    case (PENTA_33)
      return_value = pts( pris33(:) )
    case (PENTA_38)
      return_value = pts( pris38(:) )
    case (PENTA_40)
      return_value = pts( pris40(:) )
    case (PENTA_66)
      return_value = pts( pris66(:) )
    case (PENTA_75)
      return_value = pts( pris75(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_cgns_pris_pts
!
!###############################################################################
!
pure function reorder_cgns_hexa_pts(pts,elem_type) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  integer, intent(in) :: elem_type
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Parameters ..
  integer, parameter :: hexa8(8) = [1,2,4,3,5,6,8,7]
  !
  integer, parameter :: hexa20(20) = [ 1, 9, 2,12,   10, 4,11, 3, &
                                     13,   14,         16,   15, &
                                      5,17, 6,20,   18, 8,19, 7]
  !
  integer, parameter :: hexa27(27) = [ 1, 9, 2,12,21,10, 4,11, 3, &
                                     13,22,14,25,27,23,16,24,15, &
                                      5,17, 6,20,26,18, 8,19, 7]
  !
  integer, parameter :: hexa32(32) = &
    [ 1, 9,10, 2,16,      11,15,      12, 4,14,13, 3, &
     17,      19,                        23,      21, &
     18,      20,                        24,      22, &
      5,25,26, 6,32,      27,31,      28, 8,30,29, 7]
  !
  integer, parameter :: hexa56(56) = &
    [ 1, 9,10, 2,16,33,34,11,15,36,35,12, 4,14,13, 3, &
     17,37,38,19,50,      41,49,      42,23,46,45,21, &
     18,40,39,20,51,      44,52,      43,24,47,48,22, &
      5,25,26, 6,32,53,54,27,31,56,55,28, 8,30,29, 7]
  !
  integer, parameter :: hexa64(64) = &
    [ 1, 9,10, 2,16,33,34,11,15,36,35,12, 4,14,13, 3, &
     17,37,38,19,50,57,58,41,49,60,59,42,23,46,45,21, &
     18,40,39,20,51,61,62,44,52,64,63,43,24,47,48,22, &
      5,25,26, 6,32,53,54,27,31,56,55,28, 8,30,29, 7]
  !
  integer, parameter :: hexa44(44) = &
  [ 1, 9,10,11, 2,20,         12,19,         13,18,         14, 4,17,16,15, 3, &
   21,            24,                                          30,         27, &
   22,            25,                                          31,         28, &
   23,            26,                                          32,         29, &
    5,33,34,35, 6,44,         36,43,         37,42,         38, 8,41,40,39, 7]
  !
  integer, parameter :: hexa98(98) = &
  [ 1, 9,10,11, 2,20,45,46,47,12,19,52,53,48,13,18,51,50,49,14, 4,17,16,15, 3, &
   21,54,55,56,24,83,         63,82,         64,81,         65,30,74,73,72,27, &
   22,61,62,57,25,84,         70,89,         71,88,         66,31,75,80,79,28, &
   23,60,59,58,26,85,         69,86,         68,87,         67,32,76,77,78,29, &
    5,33,34,35, 6,44,90,91,92,36,43,97,98,93,37,42,96,95,94,38, 8,41,40,39, 7]
  !
  ! SAME ORDER AS ABOVE hexa44 AND hexa98 BUT COMPLETELY FILLED IN
  integer, parameter :: hexa125(125) = &
    [  1,  9, 10, 11,  2, 20, 45, 46, 47, 12, 19, 52, 53, 48, 13, 18, 51, 50, &
      49, 14,  4, 17, 16, 15,  3, 21, 54, 55, 56, 24, 83, 99,100,101, 63, 82, &
     106,107,102, 64, 81,105,104,103, 65, 30, 74, 73, 72, 27, 22, 61, 62, 57, &
      25, 84,108,109,110, 70, 89,115,116,111, 71, 88,114,113,112, 66, 31, 75, &
      80, 79, 28, 23, 60, 59, 58, 26, 85,117,118,119, 69, 86,124,125,120, 68, &
      87,123,122,121, 67, 32, 76, 77, 78, 29,  5, 33, 34, 35,  6, 44, 90, 91, &
      92, 36, 43, 97, 98, 93, 37, 42, 96, 95, 94, 38,  8, 41, 40, 39,  7]
  !
continue
  !
  select case (elem_type)
    case (HEXA_8)
      return_value = pts( hexa8(:) )
    case (HEXA_20)
      return_value = pts( hexa20(:) )
    case (HEXA_27)
      return_value = pts( hexa27(:) )
    case (HEXA_32)
      return_value = pts( hexa32(:) )
    case (HEXA_44)
      return_value = pts( hexa44(:) )
    case (HEXA_56)
      return_value = pts( hexa56(:) )
    case (HEXA_64)
      return_value = pts( hexa64(:) )
    case (HEXA_98)
      return_value = pts( hexa98(:) )
    case (HEXA_125)
      return_value = pts( hexa125(:) )
    case default
      return_value = pts(:)
  end select
  !
end function reorder_cgns_hexa_pts
!
!###############################################################################
!
pure &
function check_if_grid_requires_int64() result(return_value)
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Scalars ..
  integer      :: ii,b,z,g,s,ibc,c
  integer(i8)  :: isize,itemp
  integer(CST) :: ival
  !
  !.. Local Parameters ..
  integer(i8), parameter :: max_i4 = huge(0_i4)
  !
continue
  !
  return_value = fals
  !
  one_pass_loop: do ii = 1,1
    !
    base_loop: do b = 1,size(base)
      !
      zone_loop: do z = 1,size(base(b)%zone)
        !
        ival = max( base(b)%zone(z)%nvertex,     &
                    base(b)%zone(z)%ncell,       &
                    base(b)%zone(z)%nboundvertex )
        if (int(ival,kind=i8) >= max_i4) then
          return_value = true
          exit one_pass_loop
        end if
        !
        grid_loop: do g = 1,size(base(b)%zone(z)%grid)
          coord_loop: do c = 1,size(base(b)%zone(z)%grid(g)%coord)
            !
            isize = size( base(b)%zone(z)%grid(g)%coord(c)%x , kind=i8 )
            if (isize >= max_i4) then
              return_value = true
              exit one_pass_loop
            end if
            !
          end do coord_loop
        end do grid_loop
        !
        section_loop: do s = 1,size(base(b)%zone(z)%section)
          !
          ival = max( base(b)%zone(z)%section(s)%start, &
                      base(b)%zone(z)%section(s)%end    )
          isize = int( ival , kind=i8 )
          !
          if (allocated(base(b)%zone(z)%section(s)%elements)) then
            itemp = size( base(b)%zone(z)%section(s)%elements , kind=i8 )
            ! MAXVAL IS COSTLY
            ival = maxval( base(b)%zone(z)%section(s)%elements )
            ! MAXVAL IS COSTLY
            itemp = max( itemp , int(ival,kind=i8) )
          end if
          !
          isize = max( isize , itemp )
          !
          if (allocated(base(b)%zone(z)%section(s)%parentdata)) then
            itemp = size( base(b)%zone(z)%section(s)%parentdata , kind=i8 )
            ival = maxval( base(b)%zone(z)%section(s)%parentdata )
            itemp = max( itemp , int(ival,kind=i8) )
          end if
          !
          isize = max( isize , itemp )
          !
          if (isize >= max_i4) then
            return_value = true
            exit one_pass_loop
          end if
          !
        end do section_loop
        !
        bc_loop: do ibc = 1,size(base(b)%zone(z)%bc)
          !
          ival = max( base(b)%zone(z)%bc(ibc)%npnts,           &
                      base(b)%zone(z)%bc(ibc)%normal_list_size )
          isize = int( ival , kind=i8 )
          !
          if (allocated(base(b)%zone(z)%bc(ibc)%pnts)) then
            itemp = size( base(b)%zone(z)%bc(ibc)%pnts , kind=i8 )
            ival = maxval( base(b)%zone(z)%bc(ibc)%pnts )
            itemp = max( itemp , int(ival,kind=i8) )
          end if
          !
          isize = max( isize , itemp )
          !
          if (allocated(base(b)%zone(z)%bc(ibc)%normal_list)) then
            itemp = size( base(b)%zone(z)%bc(ibc)%normal_list , kind=i8 )
          end if
          !
          isize = max( isize , itemp )
          !
          if (isize >= max_i4) then
            return_value = true
            exit one_pass_loop
          end if
          !
        end do bc_loop
        !
      end do zone_loop
      !
    end do base_loop
    !
  end do one_pass_loop
  !
end function check_if_grid_requires_int64
!
!###############################################################################
!
subroutine cgns_memory_usage(iunit)
  !
  !.. Use Statements ..
  use module_memory
  !
  !.. Formal Arguments ..
  integer, optional, intent(in) :: iunit
  !
  !.. Local Scalars ..
  integer :: iou,i,j
  integer :: ib,iz,ig,is,ic,ibc
  character(len=36) :: array_name
  type(memory) :: array,dt
  type(memory) :: total
  type(memory) :: bdt,zdt,sdt,gdt,cdt,bcdt
  !
continue
  !
  iou = iout
  if (present(iunit)) iou = iunit
  !
  call total%reset
  !
  write (iou,1)
  !
  !-----------------------------------------------------------------------------
  ! 
  if (allocated(base)) then
    !
    call bdt%set( storage_size(base,kind=inttype) , &
                         size(base,kind=inttype) )
    !
    base_loop: do ib = lbound(base,dim=1), &
                       ubound(base,dim=1)
      !
      call bdt%add( storage_size(base(ib)%name,kind=inttype) )
      call bdt%add( storage_size(base(ib)%nsize,kind=inttype) )
      call bdt%add( storage_size(base(ib)%cell_dim,kind=inttype) )
      call bdt%add( storage_size(base(ib)%phys_dim,kind=inttype) )
      !
      if (allocated(base(ib)%zone)) then
        !
        rename_zone: associate ( zone => base(ib)%zone )
          !
          call zdt%set( storage_size(zone,kind=inttype) , &
                                size(zone,kind=inttype) )
          !
          zone_loop: do iz = lbound(zone,dim=1),ubound(zone,dim=1)
            !
            call zdt%add( storage_size(zone(iz)%name,kind=inttype) )
            call zdt%add( storage_size(zone(iz)%nsize,kind=inttype) )
            call zdt%add( storage_size(zone(iz)%nvertex,kind=inttype) )
            call zdt%add( storage_size(zone(iz)%ncell,kind=inttype) )
            call zdt%add( storage_size(zone(iz)%nboundvertex,kind=inttype) )
            call zdt%add( storage_size(zone(iz)%type,kind=inttype) )
            call zdt%add( storage_size(zone(iz)%index_dim,kind=inttype) )
            !
!###############################################################################
            if (allocated(zone(iz)%grid)) then
              !
              rename_grid: associate ( grid => zone(iz)%grid )
                !
                call gdt%set( storage_size(grid,kind=inttype) , &
                                      size(grid,kind=inttype) )
                !
                grid_loop: do ig = lbound(grid,dim=1),ubound(grid,dim=1)
                  !
                  call gdt%add( storage_size(grid(ig)%name,kind=inttype) )
                  !
                  if (allocated(grid(ig)%coord)) then
                    !
                    rename_coord: associate ( coord => grid(ig)%coord )
                      !
                      call cdt%set( storage_size(coord,kind=inttype) , &
                                            size(coord,kind=inttype) )
                      !
                      coord_loop: do ic = lbound(coord,dim=1), &
                                          ubound(coord,dim=1)
                        !
                        call cdt%add(storage_size(coord(ic)%name,kind=inttype))
                        !
                        if (allocated(coord(ic)%x)) then
                          !
                          call cdt%add(storage_size(coord(ic)%x,kind=inttype), &
                                               size(coord(ic)%x,kind=inttype))
                          !
                        end if
                        !
                      end do coord_loop
                      !
                    end associate rename_coord
                    !
                    write (array_name,11) ib,iz,"grid",ig,"coord"
                    write (iou,2) array_name,cdt%get_units()
                    !
                    call gdt%add(cdt)
                    !
                  else
                    write (array_name,11) ib,iz,"grid",ig,"coord"
                    write (iou,7) array_name
                  end if
                  !
                end do grid_loop
                !
              end associate rename_grid
              !
              write (array_name,11) ib,iz,"grid"
              write (iou,2) array_name,gdt%get_units()
              call zdt%add(gdt)
              !
            else
              write (array_name,11) ib,iz,"grid"
              write (iou,7) array_name
            end if
!###############################################################################
            !
!###############################################################################
            if (allocated(zone(iz)%section)) then
              !
              rename_section: associate ( section => zone(iz)%section )
                !
                call sdt%set( storage_size(section,kind=inttype) , &
                                      size(section,kind=inttype) )
                !
                section_loop: do is = lbound(section,dim=1), &
                                      ubound(section,dim=1)
                  !
                  call sdt%add( storage_size(section(is)%name,kind=inttype) )
                  call sdt%add( storage_size(section(is)%type,kind=inttype) )
                  call sdt%add( storage_size(section(is)%start,kind=inttype) )
                  call sdt%add( storage_size(section(is)%end,kind=inttype) )
                  call sdt%add( storage_size(section(is)%nbndry,kind=inttype) )
                  call sdt%add( storage_size(section(is)%parent_flag, &
                                             kind=inttype) )
                  if (allocated(section(is)%elements)) then
                    call sdt%add( &
                           storage_size(section(is)%elements,kind=inttype) , &
                                   size(section(is)%elements,kind=inttype) )
                  end if
                  if (allocated(section(is)%parentdata)) then
                    call sdt%add( &
                           storage_size(section(is)%parentdata,kind=inttype) , &
                                   size(section(is)%parentdata,kind=inttype) )
                  end if
                  !
                end do section_loop
                !
              end associate rename_section
              !
              write (array_name,11) ib,iz,"section"
              write (iou,2) array_name,sdt%get_units()
              call zdt%add(sdt)
              !
            else
              write (array_name,11) ib,iz,"section"
              write (iou,7) array_name
            end if
!###############################################################################
            !
!###############################################################################
            if (allocated(zone(iz)%bc)) then
              !
              rename_bc: associate ( bc => zone(iz)%bc )
                !
                call bcdt%set( storage_size(bc,kind=inttype) , &
                                       size(bc,kind=inttype) )
                !
                bc_loop: do ibc = lbound(bc,dim=1),ubound(bc,dim=1)
                  !
                  call bcdt%add( storage_size(bc(ibc)%name,kind=inttype) )
                  call bcdt%add( storage_size(bc(ibc)%boco_type,kind=inttype) )
                  call bcdt%add( storage_size(bc(ibc)%ptset_type,kind=inttype) )
                  call bcdt%add( storage_size(bc(ibc)%normal_datatype, &
                                              kind=inttype) )
                  call bcdt%add( storage_size(bc(ibc)%location,kind=inttype) )
                  call bcdt%add( storage_size(bc(ibc)%npnts,kind=inttype) )
                  call bcdt%add( storage_size(bc(ibc)%normal_list_size, &
                                              kind=inttype) )
                  call bcdt%add( &
                              storage_size(bc(ibc)%normal_index,kind=inttype) ,&
                                      size(bc(ibc)%normal_index,kind=inttype) )
                  call bcdt%add( storage_size(bc(ibc)%ndataset,kind=inttype) )
                  if (allocated(bc(ibc)%pnts)) then
                    call bcdt%add( storage_size(bc(ibc)%pnts,kind=inttype) , &
                                           size(bc(ibc)%pnts,kind=inttype) )
                  end if
                  if (allocated(bc(ibc)%normal_list)) then
                    call bcdt%add( &
                              storage_size(bc(ibc)%normal_list,kind=inttype) , &
                                      size(bc(ibc)%normal_list,kind=inttype) )
                  end if
                  !
                end do bc_loop
                !
              end associate rename_bc
              !
              write (array_name,11) ib,iz,"bc"
              write (iou,2) array_name,bcdt%get_units()
              call zdt%add(bcdt)
              !
            else
              write (array_name,11) ib,iz,"bc"
              write (iou,7) array_name
            end if
!###############################################################################
            !
          end do zone_loop
          !
          call bdt%add(zdt)
          write (array_name,10) ib
          write (iou,2) array_name,zdt%get_units()
          !
        end associate rename_zone
        !
      else
        write (iou,7) "base",ib,"zone"
      end if
      !
    end do base_loop
    !
    call total%add(bdt)
    write (array_name,5) "base"
    write (iou,2) array_name,bdt%get_units()
    !
  else
    !
    write (iou,7) "base"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  if (allocated(cell_grps)) then
    !
    call dt%set( storage_size(cell_grps,kind=inttype) , &
                         size(cell_grps,kind=inttype) )
    !
    do i = lbound(cell_grps,dim=1),ubound(cell_grps,dim=1)
      !
      call dt%add( storage_size(cell_grps(i)%cell_dim,kind=inttype) )
      call dt%add( storage_size(cell_grps(i)%grp_id,kind=inttype) )
      call dt%add( storage_size(cell_grps(i)%ibc,kind=inttype) )
      call dt%add( storage_size(cell_grps(i)%name,kind=inttype) )
      !
      if (allocated(cell_grps(i)%cells)) then
        call dt%add( storage_size(cell_grps(i)%cells,kind=inttype) , &
                             size(cell_grps(i)%cells,kind=inttype) )
      end if
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "cell_grps"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "cell_grps"
    !
  end if
  !
  !-----------------------------------------------------------------------------
  !
  !
  ! Check the memory of the arrays within this module
  !
  !---------------------------------------------------------------------
  if (allocated(grid_data)) then
    call array%set( storage_size(grid_data,kind=inttype) , &
                            size(grid_data,kind=inttype) )
    call total%add(array)
    write (array_name,5) "grid_data"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "grid_data"
  end if
  !---------------------------------------------------------------------
  if (allocated(section_data)) then
    call array%set( storage_size(section_data,kind=inttype) , &
                            size(section_data,kind=inttype) )
    call total%add(array)
    write (array_name,5) "section_data"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "section_data"
  end if
  !---------------------------------------------------------------------
  if (allocated(boco_data)) then
    call array%set( storage_size(boco_data,kind=inttype) , &
                            size(boco_data,kind=inttype) )
    call total%add(array)
    write (array_name,5) "boco_data"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "boco_data"
  end if
  !---------------------------------------------------------------------
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module module_cgns")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  4 format (a,"(",i0,")%",a)
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  10 format ("base(",i0,")%zone")
  11 format ("base(",i0,")%zone(",i0,")%",a,:,"(",i0,")%",a)
  !
end subroutine cgns_memory_usage
!
!###############################################################################
!
include "Functions/bc_string_to_integer.f90"
include "Functions/last.f90"
!
!###############################################################################
!
end module module_cgns
!
!###############################################################################
!
