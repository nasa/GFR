module module_cgns
  !
#define PROFILE_ON
!!!#define DEBUG_HOST_HANG
  !.. Use Statements
  use module_kind_types, global_debug_timer => debug_timer
  use module_cgns_types
  use geovar, only : grid_t
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
  logical(lk), parameter :: create_nozzle_interior_cgns_gridfile = fals
 !logical(lk), parameter :: create_nozzle_interior_cgns_gridfile = true
  !
  ! ###################################
  ! ###  Module Derived Data Types  ###
  ! ###################################
  !
  type :: cg_t
    character(len=CGLEN) :: name
  end type cg_t
  !
  type, extends(cg_t) :: cg_sz_t
   !character(len=CGLEN) :: name  ! INHERITED FROM CG_T
    integer(CST) :: nsize
  end type cg_sz_t
  !
  type, extends(cg_t) :: cg_dataarray_t
   !character(len=CGLEN) :: name  ! INHERITED FROM CG_T
    real(wp), allocatable :: data(:)
  end type cg_dataarray_t
  !
  type, extends(cg_t) :: cg_bc_t
   !character(len=CGLEN) :: name  ! INHERITED FROM CG_T
    integer(CET) :: boco_type, ptset_type, normal_datatype, location
    integer(CST) :: npnts, normal_list_size
    integer(CBT) :: normal_index(3), ndataset
    integer(CST), allocatable :: pnts(:)!, faces(:)
    real(wp),     allocatable :: normal_list(:)
  end type cg_bc_t
  !
  type, extends(cg_t) :: cg_coord_t
   !character(len=CGLEN) :: name  ! INHERITED FROM CG_T
    real(wp), allocatable :: x(:)
  end type cg_coord_t
  !
  type, extends(cg_t) :: cg_grid_t
   !character(len=CGLEN) :: name  ! INHERITED FROM CG_T
    type(cg_coord_t), allocatable :: coord(:)
  end type cg_grid_t
  !
  type, extends(cg_t) :: cg_section_t
   !character(len=CGLEN) :: name  ! INHERITED FROM CG_T
    integer(CET) :: type
    integer(CST) :: start, end
    integer(CBT) :: nbndry, parent_flag
    integer(CST), allocatable :: elements(:), parentdata(:)
   !integer(CET), allocatable :: elemtypes(:)
  end type cg_section_t
  !
  type, extends(cg_sz_t) :: cg_zone_t
   !character(len=CGLEN) :: name  ! INHERITED FROM CG_T
   !integer(CST) :: nsize         ! INHERITED FROM CG_SZ_T
    integer(CST) :: nvertex, ncell, nboundvertex
    integer(CET) :: type
    integer(CBT) :: index_dim
    type(cg_grid_t),    allocatable :: grid(:)
    type(cg_section_t), allocatable :: section(:)
    type(cg_bc_t),      allocatable :: bc(:)
  end type cg_zone_t
  !
  type, extends(cg_sz_t) :: cg_base_t
   !character(len=CGLEN) :: name  ! INHERITED FROM CG_T
   !integer(CST) :: nsize         ! INHERITED FROM CG_SZ_T
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
  type :: cell_grps_t
    integer :: cell_dim
    integer :: grp_id
    integer :: ibc
    character(len=CGLEN) :: name
    integer, allocatable :: cells(:)
  end type cell_grps_t
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
  use ovar, only : itestcase
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
  if (mypnum == glb_root) then
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
  ! Reset the MPI communicator for CGNS back to MPI_COMM_WORLD
  !
  call cgp_mpi_comm_f(mpi_comm_world_val,cgierr)
 !call cgns_error(pname,gridfile,"cgp_mpi_comm_f",cgierr,__LINE__,__FILE__)
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  !
  !
  if (mypnum == glb_root) then
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
    if (itestcase == Nozzle_Jet) then
      if (create_nozzle_interior_cgns_gridfile) then
        call extract_nozzle_grid
      else
        call find_elem_for_initial_jet_profile
      end if
    end if
    !
    ! Create the arrays needed by the solver that define the grid
    ! geometry: bface, nodes_of_cell_ptr, nodes_of_cell, etc.
    !
    call create_solver_geom_arrays_cgns
    !
  end if
  !
  ! Broadcast the grid derived type to all processors
  ! but save it on only the root processes of each host
  !
  call root_broadcasts_grid_dt
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
  ! Broadcast the grid arrays to all processors
  !
  if (ncpu > 1) then
   !call broadcast_grid_arrays
    call broadcast_solver_geom_arrays
  else
    call create_serial_geom_arrays
  end if
  if (debug_cgns) then
    if (mypnum == glb_root) then
      write (*,*)
      write (*,*) "Leaving read_cgns_gridfile! Press any key to continue"
      read (*,*)
    end if
  end if
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
 !call stop_gfr(stop_mpi,pname,__LINE__,__FILE__, &
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
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! FOR NOW: Abort if there is more than one base in the CGNS file
  !
  if (allocated(base)) then
    if (size(base) > 1) then
      write (error_message,2)
      call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
    end if
  else
    write (error_message,12)
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! FOR NOW: Abort if there is more than one zone in the CGNS file for any base
  !
  do ib = 1,size(base)
    if (allocated(base(ib)%zone)) then
      if (size(base(ib)%zone) > 1) then
        write (error_message,3)
        call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
      end if
    else
      write (error_message,13)
      call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
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
            call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
          end if
        end do
      else
        write (error_message,14)
        call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
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
          call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
        end if
      else
        write (error_message,15)
        call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
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
              call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
            end if
          end do
        else
          write (error_message,16)
          call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
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
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
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
  integer(CET) :: read_datatype
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_cgns_coord"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  read_datatype = merge(REALDOUBLE,REALSINGLE,wp==r8)
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
  integer(CET) :: read_datatype
  !
  !.. Local Arrays ..
  integer(CST) :: dimension_vector(1:12) ! NOTE: 12 is the maximum number of
  !                                              data dimensions for a
  !                                              DataArray_t node defined by the
  !                                              CGNS standard
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_cgns_additional_coord"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  read_datatype = merge(REALDOUBLE,REALSINGLE,wp==r8)
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
  !.. Use Statements ..
  use iso_c_binding, only : c_null_ptr
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
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  if (base(ib_g)%cell_dim == 2 .and. bc%location /= EDGECENTER) then
    write (error_message,1) 2,ibc_g,trim(adjustl(bc%name)),"Edge"
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
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
  !.. Use Statements ..
  use geovar, only : grid
  use ovar,   only : grid_scaling_factor
  use ovar,   only : governing_equations
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
  character(len=9)     :: wall_type
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
  wall_type = 'SLIPWALL'
  if (governing_equations == NavierStokes_Eqns) then
    wall_type = 'ADIABATIC'
  end if
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
  ! Apply the grid scaling factor to the coordinates of the grid nodes
  !
  do n = 1,size(grid%xyz,dim=2)
    do i = 1,size(grid%xyz,dim=1)
      grid%xyz(i,n) = grid_scaling_factor*grid%xyz(i,n)
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
      if (mypnum == glb_root) write (iout,11) mine,maxe
    end if
    if (sume /= size(mske)) then
      if (mypnum == glb_root) write (iout,12) sume,size(mske)
    end if
    write (error_message,13)
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
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
            grid%elem(n)%tags = [0,0] ! using F2003 automatic reallocation
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
            grid%elem(n)%tags = [0,0] ! using F2003 automatic reallocation
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
        j1 = 3
        bc_string = BCTYPENAME( base(ib)%zone(iz)%bc(ibc)%boco_type )
        if (uppercase(bc_string) == 'NULL') then
          j1 = 1
          bc_string = trim(adjustl(base(ib)%zone(iz)%bc(ibc)%name))
        end if
        j2 = len_trim(bc_string)
        !
        !#################################################################
        !#################################################################
        !#################################################################
        !    TEMPORARY FOR POINTWISE CURVED MESH OF SMC000 NOZZLE
        !
        if (uppercase(bc_string(j1:j2)) == 'FARFIELD') then
          bc_string = 'CHARACTERISTIC'
        else if (uppercase(bc_string(j1:j2)) == 'INFLOW') then
         !bc_string = 'REFERENCE'
          bc_string = 'CHARACTERISTIC'
        else if (uppercase(bc_string(j1:j2)) == 'INFLOWSUBSONIC') then
         !bc_string = 'SUBINFLOW'
          bc_string = 'CHARACTERISTIC'
        else if (uppercase(bc_string(j1:j2)) == 'NOZZLEBASE') then
          bc_string = trim(adjustl(wall_type))
        else if (uppercase(bc_string(j1:j2)) == 'NOZZLEEXTERIOR') then
          bc_string = trim(adjustl(wall_type))
        else if (uppercase(bc_string(j1:j2)) == 'NOZZLEINTERIOR') then
          bc_string = trim(adjustl(wall_type))
        else if (uppercase(bc_string(j1:j2)) == 'WALL') then
          bc_string = trim(adjustl(wall_type))
        else if (uppercase(bc_string(j1:j2)) == 'OUTFLOW') then
          bc_string = 'CHARACTERISTIC'
        else if (uppercase(bc_string(j1:j2)) == 'OUTFLOWSUBSONIC') then
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
  use geovar,    only : grid
  use ovar,      only : bc_names
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
    call stop_gfr(abort,pname,__LINE__,__FILE__)
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
    call stop_gfr(abort,pname,__LINE__,__FILE__)
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
      call stop_gfr(abort,pname,__LINE__,__FILE__)
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
    call stop_gfr(abort,pname,__LINE__,__FILE__)
  end if
  !
  nfbnd = 0
  do i = 1,size(bnd_grps)
    nfbnd = nfbnd + size( cell_grps(bnd_grps(i))%cells )
  end do
  !
  if (nfbnd /= nfbnd_correct) then
    write (iout,13) nfbnd,nfbnd_correct
    call stop_gfr(abort,pname,__LINE__,__FILE__)
  end if
  !
  ! Check that the number of interior cells and boundary faces are correct
  !
  if ( ncell+nfbnd /= size(grid%elem) ) then
    write (iout,14) ncell,nfbnd,size(grid%elem)
    call stop_gfr(abort,pname,__LINE__,__FILE__)
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
  allocate ( bc_names(1:size(bnd_grps)) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bc_names",1,__LINE__,__FILE__,ierr,error_message)
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
    bc_names(i) = trim(adjustl(cell_grps(np)%name))
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
      bface( 6,nf) = 0                       ! index for bc_in array
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
  ! Deallocate bnd_grps and cell_grps since we dont need them anymore
  !
  deallocate ( cell_grps , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cell_grps",2,__LINE__,__FILE__,ierr,error_message)
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
    call stop_gfr(abort,pname,__LINE__,__FILE__)
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
subroutine root_broadcasts_grid_dt()
  !
  !.. Use Statements ..
  use geovar, only : grid
  !
  !.. Local Scalars ..
  integer          :: n1,n2,ierr
  integer(INT_MPI) :: n
  !
  !.. Local Arrays ..
  integer :: array_sizes(1:4)
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: grid_elem(:)
  real(wp), allocatable :: grid_xyz(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "root_broadcasts_grid_dt"
 !logical(lk),      parameter :: use_all_processors = fals
  logical(lk),      parameter :: use_all_processors = true
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Have the global root process copy the information in the grid derived type
  ! to the grid_xyz and grid_elem arrays, and save the necessary size
  ! information into the array_sizes array.
  !
#ifdef DEBUG_HOST_HANG
  if (i_am_host_root) then
    write (iout,100) mypnum,"b) before root_packs_grid_dt, before mpi_barrier"
    call mpi_barrier(host_roots_comm,mpierr)
    write (iout,100) mypnum, &
            "c) before root_packs_grid_dt, after  mpi_barrier",mpierr
  end if
#endif
  100 format (" CPU # ",i4.4," : ",a,:,"  mpierr=",i0)
  101 format (" CPU # ",i4.4," : ",a,:,i0)
  102 format (a," : size=",i0,",   storage_size=",i0)
  if (mypnum == glb_root) then
    !
    call root_packs_grid_dt(grid_elem,grid_xyz,grid)
    !
    array_sizes(1) = size(grid_elem)
    array_sizes(2) = size(grid_xyz,dim=1)
    array_sizes(3) = size(grid_xyz,dim=2)
    array_sizes(4) = size(grid%elem)
    !
  end if
#ifdef DEBUG_HOST_HANG
  if (i_am_host_root) then
    if (mypnum == glb_root) then
      write (iout,102) "grid_elem",size(grid_elem),storage_size(grid_elem)
      write (iout,102) "grid_xyz ",size(grid_xyz),storage_size(grid_xyz)
      flush (iout)
    end if
    write (iout,100) mypnum,"d) after  root_packs_grid_dt, before mpi_barrier"
    call mpi_barrier(host_roots_comm,mpierr)
    write (iout,100) mypnum, &
            "e) after  mpi_barrier, before mpi_bast(array_sizes)",mpierr
  end if
#endif
  !
  ! Broadcast the array_sizes array to all host roots
  !
  n = size(array_sizes,kind=int_mpi)
 !if (use_all_processors) then
    call mpi_bcast(array_sizes,n,mpi_inttyp,glb_root,MPI_COMM_WORLD,mpierr)
 !else
 !  if (i_am_host_root) then
 !    call mpi_bcast(array_sizes,n,mpi_inttyp,glb_root,host_roots_comm,mpierr)
 !  end if
 !end if
#ifdef DEBUG_HOST_HANG
  if (i_am_host_root) then
    write (iout,100) mypnum, &
            "f) after  mpi_bast(array_sizes), before allocation",mpierr
  end if
#endif
  !
  ! Allocate the grid_elem array on the non-root processes
  ! so it can be received from the root process.
  !
  if (.not. allocated(grid_elem)) then
    n1 = array_sizes(1)
    allocate ( grid_elem(1:n1) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid_elem",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
#ifdef DEBUG_HOST_HANG
  if (i_am_host_root) then
    write (iout,100) mypnum,"g) after  allocation, before mpi_bast(grid_elem)"
  end if
#endif
  !
  ! Broadcast the grid_elem array to all host roots
  !
  n = size(grid_elem,kind=int_mpi)
 !if (mypnum == glb_root .or. mypnum == 1064 .or. mypnum == 1092) then
 !  write (iout,101) mypnum,"array_sizes(1) = ",array_sizes(1)
 !  write (iout,101) mypnum,"array_sizes(2) = ",array_sizes(2)
 !  write (iout,101) mypnum,"array_sizes(3) = ",array_sizes(3)
 !  write (iout,101) mypnum,"array_sizes(4) = ",array_sizes(4)
 !  write (iout,101) mypnum,"size(grid_elem) = ",size(grid_elem)
 !  write (iout,101) mypnum,"size(grid_xyz,dim=1) = ",size(grid_xyz,dim=1)
 !  write (iout,101) mypnum,"size(grid_xyz,dim=2) = ",size(grid_xyz,dim=2)
 !end if
  if (use_all_processors) then
    call mpi_bcast(grid_elem,n,mpi_inttyp,glb_root,MPI_COMM_WORLD,mpierr)
    if (.not. i_am_host_root) then
      deallocate ( grid_elem , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"grid_elem",2,__LINE__,__FILE__,ierr,error_message)
    end if
  else
    if (i_am_host_root) then
      call mpi_bcast(grid_elem,n,mpi_inttyp,glb_root,host_roots_comm,mpierr)
    end if
  end if
#ifdef DEBUG_HOST_HANG
  if (i_am_host_root) then
    write (iout,100) mypnum, &
            "h) after  mpi_bcast(grid_elem), before mpi_bast(grid_xyz)",mpierr
  end if
#endif
  !
  ! Allocate the grid_xyz array on the non-root processes
  ! so it can be received from the root process.
  !
  if (.not. allocated(grid_xyz)) then
    n1 = array_sizes(2)
    n2 = array_sizes(3)
    allocate ( grid_xyz(1:n1,1:n2) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid_xyz",1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Broadcast the grid_xyz array to all host roots
  !
  n = size(grid_xyz,kind=int_mpi)
  if (use_all_processors) then
    call mpi_bcast(grid_xyz,n,mpi_flttyp,glb_root,MPI_COMM_WORLD,mpierr)
    if (.not. i_am_host_root) then
      deallocate ( grid_xyz , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"grid_xyz",2,__LINE__,__FILE__,ierr,error_message)
    end if
  else
    call mpi_bcast(grid_xyz,n,mpi_flttyp,glb_root,host_roots_comm,mpierr)
#ifdef DEBUG_HOST_HANG
    if (i_am_host_root) then
      write (iout,100) mypnum, &
              "i) after  mpi_bcast(grid_xyz), before unpack_grid_dt",mpierr
    end if
#endif
  end if
  !
  ! Have all host roots unpack the grid_elem and grid_xyz arrays
  ! into the global grid derived type.
  ! NOTE: This also needs to be done on the global root process since move_alloc
  !       was used to copy grid%xyz over to grid_xyz, which leaves grid%xyz
  !       deallocated after the copy.
  !
  if (i_am_host_root) then
    call unpack_grid_dt(array_sizes(4),grid_elem,grid_xyz,grid)
#ifdef DEBUG_HOST_HANG
    write (iout,100) mypnum,"j) after unpack_grid_dt, before mpi_barrier"
    call mpi_barrier(host_roots_comm,mpierr)
    write (iout,100) mypnum,"k) after unpack_grid_dt, after mpi_barrier",mpierr
#endif
  end if
  !
  ! The local allocatable arrays grid_xyz and grid_elem should
  ! have been deallocated at the end of the previous routine.
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine root_broadcasts_grid_dt
!
!###############################################################################
!
subroutine root_packs_grid_dt(grid_elem,grid_xyz,grid)
  !
  !.. Use Statements ..
  use geovar, only : grid_t
  !
  !.. Formal Arguments ..
  integer,      allocatable, intent(inout) :: grid_elem(:)
  real(wp),     allocatable, intent(inout) :: grid_xyz(:,:)
  type(grid_t), allocatable, intent(inout) :: grid
  !
  !.. Local Scalars ..
  integer :: i,j,n,nc,nge,ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "root_packs_grid_dt"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Deallocate the grid_elem and grid_xyz arrays
  !
  if (allocated(grid_elem)) then
    deallocate ( grid_elem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid_elem",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(grid_xyz)) then
    deallocate ( grid_xyz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid_xyz",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Find the size needed for grid_elem array
  !
  call debug_timer(start_timer,"find the size needed for grid_elem")
  nge = 0
  do nc = 1,size(grid%elem)
    !
    i = 1       ! 1 for elem_type, size(tags) combo
    i = i + 1 ! 1 for host cell
    !
    ! Add the size of grid%elem(:)%pts to the counter for this element
    !
    if (allocated(grid%elem(nc)%pts)) then
      ! NOTE: Make sure that the counter is increased by at least 1 if the
      !       pts array is allocated but its size is 0.
      i = i + max( 1 , size(grid%elem(nc)%pts) )
    else
      i = i + 1
    end if
    !
    ! Add the size of grid%elem(:)%tags to the counter for this element
    ! if we are sending this data to the other processors
    !
    if (allocated(grid%elem(nc)%tags)) then
      i = i + size(grid%elem(nc)%tags)
    else
      ! We will know if tags data exists for this element from the
      ! elem_type, size(tags) combo, so we dont need to add 1 here
    end if
    !
    ! Add the size of grid%elem(:)%nodes to the counter for this element
    !
    if (allocated(grid%elem(nc)%nodes)) then
      ! NOTE: Make sure that the counter is increased by at least 1 if the
      !       nodes array is allocated but its size is 0.
      i = i + max( 1 , size(grid%elem(nc)%nodes) )
    else
      i = i + 1
    end if
    !
    nge = nge + i
    !
  end do
  call debug_timer(stop_timer,"find the size needed for grid_elem")
  call debug_timer(start_timer,"creating grid_xyz")
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Create grid_xyz
  !
  ! Copy all of grid%xyz into grid_xyz
  !
  ! USE MOVE_ALLOC TO PREVENT HAVING TO MAKE A COPY OF THE grid%xyz
  ! ARRAY JUST TO BROADCAST IT TO THE OTHER HOST ROOTS. ONCE THE BROADCAST
  ! HAS COMPLETED, USE ANOTHER MOVE_ALLOC TO MOVE grid_xyz BACK INTO
  ! grid%xyz.
  !
  call move_alloc( from=grid%xyz , to=grid_xyz )
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Create grid_elem
  !
  ! Allocate the grid_elem array
  !
  call debug_timer(stop_timer,"creating grid_xyz")
  call debug_timer(start_timer,"creating grid_elem")
  allocate ( grid_elem(1:nge) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid_elem",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Loop over the elements for this partition and store the data for
  ! each element in the grid_elem array
  !
  n = 0
  do nc = 1,size(grid%elem)
    !
    ! Since the values of grid%elem(glb_cell)%prop%elem_type and
    ! size(grid%elem(glb_cell)%tags) should both be at most 16-bit
    ! integers, combine them into a single default sized integer (probably
    ! 32-bit) to reduce the size needed for the grid_elem array.
    !
    ! FOR ELEM_TYPE, SIZE(TAGS) COMBO:
    ! combo = 1000*size(tags) + elem_type
    ! elem_type = mod(combo,1000)
    ! ntags = combo/1000
    !
    i = 1 ! 1 for elem_type, size(tags) combo
    j = grid%elem(nc)%prop%elem_type
    if (allocated(grid%elem(nc)%tags)) then
      j = j + 1000 * size(grid%elem(nc)%tags)
    end if
    grid_elem(n+i) = j
    !
    ! Add grid%elem(:)%host_cell if needed
    !
    i = i + 1 ! 1 for host cell
    grid_elem(n+i) = grid%elem(nc)%host_cell
    !
    ! Add grid%elem(:)%pts
    !
    ! Start by setting the next value of grid_elem to 0 in case the pts
    ! component isnt allocated or its allocated with a size of 0. This will
    ! indicate to the receiving processor that the grid%elem(:)%pts
    ! array didnt contain any data for this element.
    !
    i = i + 1
    grid_elem(n+i) = 0
    !
    ! Check to make sure grid%elem(nc)%pts is
    ! allocated before trying to add it to grid_elem
    !
    if (allocated(grid%elem(nc)%pts)) then
      !
      ! The grid%elem(:)%pts array is allocated, now make sure
      ! its size is greater than zero and contains data
      ! NOTE: We need to perform this check so that we can decrease the
      !       counter for this element before we start adding stuff to
      !       grid_elem.
      !
      if (size(grid%elem(nc)%pts) > 0) then
        !
        ! Its size is greater than 0, so first decrease the counter for this
        ! element by 1, so the previous assignment of 0 to what was the next
        ! value of grid_elem can be overwritten
        !
        i = i - 1
        !
        ! Now add grid%elem(:)%pts to grid_elem
        !
        do j = 1,size(grid%elem(nc)%pts)
          i = i + 1
          grid_elem(n+i) = grid%elem(nc)%pts(j)
        end do
        !
      end if
      !
    end if
    !
    ! Add grid%elem(:)%tags if its allocated and we
    ! are sending this data to the other processors
    !
    if (allocated(grid%elem(nc)%tags)) then
      do j = 1,size(grid%elem(nc)%tags)
        i = i + 1
        grid_elem(n+i) = grid%elem(nc)%tags(j)
      end do
    end if
    !
    ! Finally, add grid%elem(:)%nodes
    !
    ! Start by setting the next value of grid_elem to 0 in case the nodes
    ! component isnt allocated or its allocated with a size of 0. This will
    ! indicate to the receiving processor that the grid%elem(:)%nodes
    ! array didnt contain any data for this element.
    !
    i = i + 1
    grid_elem(n+i) = 0
    !
    ! Check to make sure grid%elem(nc)%nodes is
    ! allocated before trying to add it to grid_elem
    !
    if (allocated(grid%elem(nc)%nodes)) then
      !
      ! The grid%elem(:)%nodes array is allocated, now make sure
      ! its size is greater than zero and contains data
      ! NOTE: We need to perform this check so that we can decrease the
      !       counter for this element before we start adding stuff to
      !       grid_elem.
      !
      if (size(grid%elem(nc)%nodes) > 0) then
        !
        ! Its size is greater than 0, so first decrease the counter for this
        ! element by 1, so the previous assignment of 0 to what was the next
        ! value of grid_elem can be overwritten
        !
        i = i - 1
        !
        ! Now add grid%elem(:)%nodes to grid_elem
        !
        do j = 1,size(grid%elem(nc)%nodes)
          i = i + 1
          grid_elem(n+i) = grid%elem(nc)%nodes(j)
        end do
        !
      end if
      !
    end if
    !
    ! Add the final counter for this element to the grid_elem
    ! counter before going to the next element
    !
    n = n + i
    !
  end do
  call debug_timer(stop_timer,"creating grid_elem")
  !
  ! Check to make sure the final counter matches the size of the grid_elem array
  !
  if (n /= size(grid_elem)) then
    write (error_message,1) n,size(grid_elem)
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("ERROR CREATING grid_elem TO BROADCAST TO HOST ROOTS!!! ", &
            "THE FINAL COUNTER n=",i0,", WHEN IT SHOULD BE ", &
            "size(grid_elem)=",i0,"!!!")
  !
end subroutine root_packs_grid_dt
!
!###############################################################################
!
subroutine unpack_grid_dt(num_cells,grid_elem,grid_xyz,grid)
  !
  !.. Use Statements ..
  use geovar, only : grid_t
  !
  !.. Formal Arguments ..
  integer,                      intent(in) :: num_cells
  integer,      allocatable, intent(inout) :: grid_elem(:)
  real(wp),     allocatable, intent(inout) :: grid_xyz(:,:)
  type(grid_t), allocatable, intent(inout) :: grid
  !
  !.. Local Scalars ..
  integer :: n,nc,nn,np,n1,n2,ierr
  integer :: i,ntags,elem_type
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "unpack_grid_dt"
  !
continue
  !
  ! First allocate the grid derived type
  !
  if (allocated(grid)) then
    deallocate ( grid , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  allocate ( grid , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the grid%xyz array
  !
  ! USE MOVE_ALLOC (AVOIDS COPYING DATA TO NEW MEMORY)
  call move_alloc( from=grid_xyz , to=grid%xyz )
  ! USE STANDARD ALLOCATION (COPIES DATA TO NEW MEMORY)
 !allocate ( grid%xyz , source=grid_xyz , stat=ierr , errmsg=error_message )
 !call alloc_error(pname,"grid%xyz",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the grid%elem array
  !
  allocate ( grid%elem(1:num_cells) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid%elem",1,__LINE__,__FILE__,ierr,error_message)
  !
  n = 0
  do nc = 1,num_cells
    !
    ! Unpack the elem_type/ntags combo
    !
    i = 1
    elem_type = mod( grid_elem(n+i) , 1000 )
    ntags = grid_elem(n+i) / 1000
    !
    ! Get the element properties for this type of element
    !
    grid%elem(nc)%prop = element_properties(elem_type)
    !
    ! Unpack the host cell if this is being used
    !
    i = i + 1
    grid%elem(nc)%host_cell = grid_elem(n+i)
    !
    ! Unpack the point indices from grid_elem if the next value in grid_elem
    ! is not 0; otherwise skip this step and increase the counter for this
    ! element by 1
    !
    if (grid_elem(n+i+1) /= 0) then
      np = grid%elem(nc)%prop%npts
      n1 = n + i + 1
      n2 = n + i + np
      grid%elem(nc)%pts = grid_elem(n1:n2) ! using F2003 auto-reallocation
      i = i + np
    else
      i = i + 1
    end if
    !
    ! Unpack the tags from grid_elem
    !
    n1 = n + i + 1
    n2 = n + i + ntags
    grid%elem(nc)%tags = grid_elem(n1:n2) ! using F2003 auto-reallocation
    i = i + ntags
    !
    ! Unpack the node indices from grid_elem if the next value in grid_elem
    ! is not 0; otherwise skip this step and increase the counter for this
    !
    if (grid_elem(n+i+1) /= 0) then
      nn = geom_nodes( grid%elem(nc)%prop%geom )
      n1 = n + i + 1
      n2 = n + i + nn
      grid%elem(nc)%nodes = grid_elem(n1:n2) ! using F2003 auto-reallocation
      i = i + nn
    else
      i = i + 1
    end if
    !
    ! Add the final counter for this element to the grid_elem
    ! counter before going to the next element
    !
    n = n + i
    !
  end do
  !
  ! Check to make sure the final counter matches the size of the grid_elem array
  !
  if (n /= size(grid_elem)) then
    write (error_message,1) mypnum+1,n,size(grid_elem)
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ! Deallocate grid_elem and grid_xyz because they are no longer needed
  !
  if (allocated(grid_elem)) then
    deallocate ( grid_elem , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid_elem",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(grid_xyz)) then
    deallocate ( grid_xyz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid_xyz",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Format Statements
  !
  1 format ("ERROR CREATING GLOBAL grid%elem ON THE HOST ROOT #",i0," ", &
            "AFTER RECEIVING THE BROADCAST FROM THE GLOBAL ROOT PROCESS!!! ", &
            "THE FINAL COUNTER n=",i0,", WHEN IT SHOULD BE ", &
            "size(grid_elem)=",i0,"!!!")
  !
end subroutine unpack_grid_dt
!
!###############################################################################
!
subroutine extract_nozzle_grid
  !
  !.. Use Statements ..
  use geovar, only : grid
  use ovar,   only : grid_scaling_factor
  !
  !.. Local Scalars ..
  integer :: l,l1,l2,m,n,ierr
  integer :: ngrps,mintag,maxtag,nsubin
  integer :: nnodes,ngpts,npts,nelem,ncells
  integer :: interior_grp,subin_grp,nsize
  !
  character(len=CGLEN) :: grp_name
  character(len=150)   :: nozzle_grid_file
  !
  !.. Local Arrays ..
  real(wp) :: xyz_ave(1:3)
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: pts(:)
  integer, allocatable :: nodes(:)
  integer, allocatable :: pts_map(:)
  integer, allocatable :: pts_rmap(:)
  integer, allocatable :: elem_map(:)
  !
  logical(lk), allocatable :: grp_still_remains(:)
  logical(lk), allocatable :: nodes_mask(:)
  logical(lk), allocatable :: pts_mask(:)
  logical(lk), allocatable :: elem_mask(:)
  !
  type(grid_t), allocatable :: nozzle
  type(grid_t), allocatable :: temp_nozzle
  !
  type(cell_grps_t), allocatable :: nozzle_grps(:)
  !
  type(cg_base_t), allocatable :: b
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "extract_nozzle_grid"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  write (nozzle_grid_file,10) grid%elem(1)%prop%order
  10 format ("P",i0,".nozzle_only.cgns")
  write (iout,1) trim(adjustl(nozzle_grid_file))
  1 format (//," Extracting the interior of the nozzle to the", &
               " new grid file '",a,"' !!",//)
  !
  ngpts = size(grid%xyz,dim=2)
  !
  subin_grp = -huge(0)
  interior_grp = -huge(0)
  find_subin_grp: do n = 1,size(cell_grps)
    !
    grp_name = trim(adjustl(uppercase(cell_grps(n)%name)))
    if (grp_name == "INTERIOR") then
      interior_grp = cell_grps(n)%grp_id
    else if (grp_name == "INFLOWSUBSONIC") then
      subin_grp = cell_grps(n)%grp_id
    end if
    !
    if (all([subin_grp,interior_grp] /= -huge(0))) exit find_subin_grp
    !
  end do find_subin_grp
  !
  !
  allocate ( elem_mask(1:size(grid%elem)) , source=fals , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"elem_mask",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( nodes_mask(1:ngpts) , source=fals , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_mask",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( pts_mask(1:ngpts) , source=fals , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"pts_mask",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Mark all the interior and boundary elements that are within
  ! the interior of the nozzle
  !
  nsubin = 0
  ncells = 0
  !
  do n = 1,size(grid%elem)
    !
    if (allocated(grid%elem(n)%tags)) then
      if (grid%elem(n)%tags(1) == subin_grp) then
        nsubin = nsubin + 1
      end if
    end if
    !
    pts = grid%elem(n)%pts(:)     ! USING F2003 AUTO-REALLOCATION
    nodes = grid%elem(n)%nodes(:) ! USING F2003 AUTO-REALLOCATION
    !
    do l = 1,size(grid%xyz,dim=1)
      xyz_ave(l) = sum(grid%xyz(l,pts(:))) / real(size(pts),kind=kind(xyz_ave))
    end do
    !
    elem_mask(n) = point_is_in_smc000_nozzle(xyz_ave)
    !
    if (elem_mask(n)) then
      if (allocated(grid%elem(n)%tags)) then
        if (grid%elem(n)%tags(1) == interior_grp) then
          ncells = ncells + 1
        end if
      end if
      !
      pts_mask(pts(:)) = true ! Mark all nodes and points
      pts_mask(nodes(:)) = fals ! Unmark all nodes
      !
      nodes_mask(nodes(:)) = true ! Mark nodes only
      !
    end if
    !
  end do
  !
  if (allocated(nodes)) then
    deallocate ( nodes , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodes",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(pts)) then
    deallocate ( pts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"pts",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Get the size parameters for the nozzle part of the grid
  !
  nnodes = count(nodes_mask)
  npts = nnodes + count(pts_mask)
  nelem = count(elem_mask)
  !
  ! THE RESULTS SHOULD BE :
  ! ncells = 83,200
  ! nelem = 88,064
  ! nelem-ncells = 4,864
  ! nsubin = 1,664
  ! nelem+nsubin = 89,728
  write (iout,2) "               Number of interior cells = ",ncells
  write (iout,2) "            Number of linear grid nodes = ",nnodes
  write (iout,2) "            Total number of grid points = ",npts
  write (iout,2) "          Number of found grid elements = ",nelem
  write (iout,2) "Number of subsonic inflow face elements = ",nsubin
  write (iout,*)
  2 format (4x,a,i0)
  !
  ! Create the nozzle object containing the subset of the grid object
  ! for all elements and points inside the nozzle interior
  !
  allocate ( nozzle , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nozzle",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Create the mapping array for the grid elements
  !
  ! USING F2003 AUTO_REALLOCATION
  elem_map = pack(intseq(1,size(elem_mask)),mask=elem_mask)
  !
  if (allocated(elem_mask)) then
    deallocate ( elem_mask , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"elem_mask",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Create the elem component of the nozzle object
  ! NOTE: We need to add nsubin to the number of marked elements in elem_mask
  !       to account for the yet to be defined face elements that define the
  !       outflow boundary condition for the nozzle interior.
  !
 !nsize = nelem+nsubin
  nsize = nelem
  allocate ( nozzle%elem(1:nsize) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nozzle%elem",1,__LINE__,__FILE__,ierr,error_message)
  !
  nozzle%elem(1:nelem) = grid%elem( elem_map(:) )
  !
  if (allocated(elem_map)) then
    deallocate ( elem_map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"elem_map",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Create the forward mapping array for all the grid points
  !
  allocate ( pts_map(1:npts) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"pts_map",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the beginning of the mapping array for linear grid nodes
  !
  pts_map(1:nnodes) = pack(intseq(1,size(nodes_mask)),mask=nodes_mask)
  !
  if (allocated(nodes_mask)) then
    deallocate ( nodes_mask , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodes_mask",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Create the remaining of the mapping array for the high-order interior points
  !
  pts_map(nnodes+1:npts) = pack(intseq(1,size(pts_mask)),mask=pts_mask)
  !
  if (allocated(pts_mask)) then
    deallocate ( pts_mask , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"pts_mask",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Create the xyz component of the nozzle object
  !
  allocate ( nozzle%xyz(1:size(grid%xyz,dim=1),1:npts) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nozzle%xyz",1,__LINE__,__FILE__,ierr,error_message)
  !
  do n = 1,size(pts_map)
    nozzle%xyz(:,n) = grid%xyz(:,pts_map(n))
  end do
  !
  ! #######
  ! ## 3 ##
  ! #######
  !
  ! Create the reverse mapping array for all points
  ! before deallocating the forward mapping array
  !
  allocate ( pts_rmap(1:ngpts) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"pts_rmap",1,__LINE__,__FILE__,ierr,error_message)
  !
  do n = 1,size(pts_map)
    pts_rmap(pts_map(n)) = n
  end do
  !
  if (allocated(pts_map)) then
    deallocate ( pts_map , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"pts_map",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Create an array to mark the remaining group ids
  !
#ifdef SPECIAL_FOR_GCC
  mintag =  huge(mintag)
  maxtag = -huge(maxtag)
  do n = lbound(grid%elem,dim=1),ubound(grid%elem,dim=1)
    mintag = min( mintag , grid%elem(n)%tags(1) )
    maxtag = max( maxtag , grid%elem(n)%tags(1) )
  end do
#else
  mintag = minval(grid%elem(:)%tags(1))
  maxtag = maxval(grid%elem(:)%tags(1))
#endif
  !
  allocate ( grp_still_remains(mintag:maxtag) , source=fals , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grp_still_remains",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  ! Localize the point indices within the elem array and count
  ! how many of each of the remaining group ids there are
  !
  do n = 1,nelem
    !
    nozzle%elem(n)%pts(:) = pts_rmap( nozzle%elem(n)%pts(:) )
    nozzle%elem(n)%nodes(:) = pts_rmap( nozzle%elem(n)%nodes(:) )
    !
    grp_still_remains(nozzle%elem(n)%tags(1)) = true
    !
  end do
  !
  if (allocated(pts_rmap)) then
    deallocate ( pts_rmap , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"pts_rmap",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! #######
  ! ## 4 ##
  ! #######
  !
  ! Count the remaining number of element groups
  ! NOTE: Add 1 to the number of remaining element groups to include
  !       the outflow boundary condition that will need to be added.
  !
  ngrps = count(grp_still_remains(:)) + 1
 !ngrps = count(grp_still_remains(:))
  !
  ! Create a new nozzle_grps object that should be a subset of the
  ! cell_grps object for just the nozzle portion of the grid.
  !
  allocate ( nozzle_grps(1:ngrps) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nozzle_grps",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Copy the respective group information from cell_grps
  ! into the remaining element groups in nozzle_grps.
  !
  l = 0
  do n = 1,size(cell_grps)
    if (grp_still_remains(cell_grps(n)%grp_id)) then
      l = l + 1
      nozzle_grps(l) = cell_grps(n)
    end if
  end do
  !
  if (allocated(grp_still_remains)) then
    deallocate ( grp_still_remains , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grp_still_remains",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  ! Add the final group, the outflow boundary condition for the nozzle
  ! interior which was previously just a group of interior faces
  !
  nozzle_grps(ngrps)%cell_dim = 2
  nozzle_grps(ngrps)%grp_id = ngrps-1
  nozzle_grps(ngrps)%ibc = bc_characteristic
  nozzle_grps(ngrps)%name = "NozzleOutflow"
  !
  ! #######
  ! ## 5 ##
  ! #######
  !
  ! Create tags for the still undefined face elements at the outflow
  ! boundary so these future elements will be found and sorted correctly
  ! in the following loop that creates the temporary nozzle object.
  !
 !do n = nelem+1,size(nozzle%elem)
 !  ! USING F2003 AUTO-REALLOCATION
 !  nozzle%elem(n)%tags = [nozzle_grps(ngrps)%grp_id,0]
 !end do
  !
  ! Create a new temporary nozzle object to copy over the elements in
  ! the current nozzle object so they are ordered by their group id.
  !
  allocate ( temp_nozzle , source=nozzle , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_nozzle",1,__LINE__,__FILE__,ierr,error_message)
  !
  l2 = 0
  do m = 1,size(nozzle_grps)-1
    !
    l1 = l2 + 1
    !
    do n = 1,size(nozzle%elem)
      if (nozzle%elem(n)%tags(1) == nozzle_grps(m)%grp_id) then
        l2 = l2 + 1
        temp_nozzle%elem(l2) = nozzle%elem(n)
      end if
    end do
    !
    ! Reset the id for the current nozzle group
    nozzle_grps(m)%grp_id = 1-m
    ! Make sure the id for these elements matches the new group id
#ifdef SPECIAL_FOR_GCC
    do l = l1,l2
      temp_nozzle%elem(l)%tags(1) = 1-m
    end do
#else
    temp_nozzle%elem(l1:l2)%tags(1) = 1-m
#endif
    !
    nozzle_grps(m)%cells = intseq(l1,l2) ! USING F2003 AUTO-REALLOCATION
    !
  end do
  !
  ! Move the new sorted temporary nozzle object over to replace the old
  ! unsorted nozzle object, which will leave temp_nozzle deallocated.
  !
  call move_alloc( from=temp_nozzle , to=nozzle )
  !
  call create_nozzle_outflow_faces(ncells,nsubin,nozzle)
  !
  nozzle_grps(ngrps)%cells = intseq(nelem+1,size(nozzle%elem))
  !
  do n = nelem+1,size(nozzle%elem)
    nozzle%elem(n)%tags = [nozzle_grps(ngrps)%grp_id,0]
  end do
  !
  ! #######
  ! ## 6 ##
  ! #######
  !
 !call create_nozzle_outflow_faces(ncells,nsubin,nozzle)
  !
  ! Reorder all the element points back into CGNS ordering
  !
  do n = 1,size(nozzle%elem)
    nozzle%elem(n) = reorder_cgns_pts(nozzle%elem(n),reverse=true)
  end do
  !
  ! Remove the grid scaling factor from the grid coordinates
  !
  do n = 1,size(nozzle%xyz,dim=2)
    do l = 1,size(nozzle%xyz,dim=1)
      nozzle%xyz(l,n) = nozzle%xyz(l,n) / grid_scaling_factor
    end do
  end do
  !
  ! Create the new CGNS base object for the nozzle subset of the grid
  !
  call create_cgns_base_for_nozzle(nozzle,nozzle_grps,b)
  !
  if (allocated(nozzle)) then
    deallocate ( nozzle , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nozzle",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(nozzle_grps)) then
    deallocate ( nozzle_grps , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nozzle_grps",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Write out the nozzle subset of the grid
  !
  call write_cgns_gridfile(nozzle_grid_file,b)
  !
  if (allocated(b)) then
    deallocate ( b , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"b",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  write (error_message,100)
  call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  100 format ("After writing out a new CGNS grid file containing only ", &
              "the region of the grid within the interior of the nozzle.")
  !
end subroutine extract_nozzle_grid
!
!###############################################################################
!
subroutine create_nozzle_outflow_faces(lcells,lfaces,nozzle)
  !
  !.. Use Statements ..
  use geovar,           only : nbfai
  use connectivity_mod, only : find_cells_surrounding_each_node
  use connectivity_mod, only : find_cells_surrounding_each_cell
  !
  !.. Formal Arguments ..
  integer,                      intent(in) :: lcells
  integer,                      intent(in) :: lfaces
  type(grid_t), allocatable, intent(inout) :: nozzle
  !
  !.. Local Scalars ..
  integer :: c0,c1,c2,i,j,n,n1,n2,ierr
  integer :: host_cell,host_geom,host_side
  integer :: new_elem,lelem,lbfac,nfnods
  integer :: missing_elem
  !
  !.. Local Arrays ..
  integer :: ip(1:4)
  integer :: face_nodes(1:4)
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: bface(:,:)
  integer, allocatable :: noc(:)
  integer, allocatable :: noc_ptr(:)
  integer, allocatable :: cwn(:)
  integer, allocatable :: cwn_ptr(:)
  integer, allocatable :: csc(:)
  integer, allocatable :: csc_ptr(:)
  integer, allocatable :: cell_geom(:)
  integer, allocatable :: nodes(:)
  !
  type(grid_t), allocatable :: temp_nozzle
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_nozzle_outflow_faces"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  lelem = size(nozzle%elem)!- lfaces ! Total number of defined elements
  lbfac = lelem - lcells             ! Total number of defined boundary faces
  !
  ! Create a temporary bface array
  !
  allocate ( bface(1:nbfai,1:lbfac) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! For this specific case, we only need bface(10,:) to be defined
  ! and it just needs to be an integer sequence from lcells+1 to lelem
  !
  bface(10,:) = intseq(lcells+1,lelem)
  !
  ! Create temporary nodes_of_cell and nodes_of_cell_ptr arrays
  !
  allocate ( noc_ptr(1:lelem+1) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"noc_ptr",1,__LINE__,__FILE__,ierr,error_message)
  !
  do n = 1,lelem
    noc_ptr(n+1) = noc_ptr(n) + size(nozzle%elem(n)%nodes)
  end do
  !
  allocate ( noc(1:last(noc_ptr)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"noc",1,__LINE__,__FILE__,ierr,error_message)
  !
  do n = 1,lelem
    n1 = noc_ptr(n)+1
    n2 = noc_ptr(n+1)
    noc(n1:n2) = nozzle%elem(n)%nodes(:)
  end do
  !
  ! Create a temporary cell_geom array
  ! USING F2003 AUTO-REALLOCATION
  cell_geom = nozzle%elem(1:lelem)%prop%geom
  !
  ! Find the cells surrounding each node
  !
  call find_cells_surrounding_each_node(bface,noc,noc_ptr,cwn,cwn_ptr)
  !
  ! Find the cells surrounding each cell
  !
  call find_cells_surrounding_each_cell(bface,cell_geom,noc,noc_ptr, &
                                        cwn,cwn_ptr,csc,csc_ptr)
  !
  ! Deallocate the local arrays which are no longer needed
  !
  if (allocated(cwn)) then
    deallocate ( cwn , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cwn",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(cwn_ptr)) then
    deallocate ( cwn_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cwn_ptr",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(bface)) then
    deallocate ( bface , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bface",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(cell_geom)) then
    deallocate ( cell_geom , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cell_geom",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  ! Make sure the number of missing surrounding
  ! cells matches the value of lfaces as expected
  !
 !if (count(csc(:) == 0) /= lfaces) then
 !  write (error_message,100) count(csc(:) == 0),lfaces
 !  call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
 !end if
  missing_elem = count(csc(:) == 0)
  !
  allocate ( temp_nozzle , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_nozzle",1,__LINE__,__FILE__,ierr,error_message)
  !
  temp_nozzle%xyz = nozzle%xyz ! USING F2003 AUTO-REALLOCATION
  !
  allocate ( temp_nozzle%elem(1:lelem+missing_elem) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"temp_nozzle",1,__LINE__,__FILE__,ierr,error_message)
  !
  temp_nozzle%elem(1:lelem) = nozzle%elem(1:lelem)
  !
  ! Find all the cells that are missing surrounding cells
  ! and add the missing cell as a new face element
  !
  new_elem = lelem
  !
  csc_cell_loop: do host_cell = 1,lcells
    !
    c0 = csc_ptr(host_cell)
    c1 = csc_ptr(host_cell)+1
    c2 = csc_ptr(host_cell+1)
    !
    if (all(csc(c1:c2) /= 0)) cycle csc_cell_loop
    !
    n1 = noc_ptr(host_cell)+1
    n2 = noc_ptr(host_cell+1)
    !
    host_geom = temp_nozzle%elem(host_cell)%prop%geom
    !
    csc_side_loop: do host_side = 1,geom_faces(host_geom)
      !
      if (csc(c0+host_side) /= 0) cycle csc_side_loop
      !
      ! Create a new elem
      !
      new_elem = new_elem + 1
      !
      ! Create the new face element from the host element
      !
      temp_nozzle%elem(new_elem) = get_face_elem(temp_nozzle%elem(host_cell), &
                                            host_side,host_cell)
      !
      ! THE REMAINDER OF THIS LOOP IS JUST TO CHECK THAT THE get_face_elem
      ! FUNCTION SEEMS TO BE CREATING THE CORRECT FACE ELEMENT
      !
      ! Find the number of nodes that need to be matched for this face
      !
      nfnods = num_face_nodes(host_side,host_geom)
      !
      ! Node offsets for this face of the current cell
      !
      ip(1:nfnods) = side_nodes(1:nfnods,host_side,host_geom)
      !
      ! Get the nodes of the current face for the current cell
      !
      face_nodes(1:nfnods) = noc( n1 + ip(1:nfnods) )
      !
      if (allocated(temp_nozzle%elem(new_elem)%nodes)) then
        if (size(temp_nozzle%elem(new_elem)%nodes) == nfnods) then
          one_pass_loop: do j = 1,1
            do i = 0,nfnods-1
              nodes = cshift(temp_nozzle%elem(new_elem)%nodes,shift=i)
              if (all(nodes == face_nodes(1:nfnods))) then
                exit one_pass_loop
              end if
            end do
            write (iout,203) new_elem
            write (iout,204) (face_nodes(i),i=1,nfnods)
            write (iout,205) (temp_nozzle%elem(new_elem)%nodes(i),i=1,nfnods)
          end do one_pass_loop
        else
          write (iout,202) size(temp_nozzle%elem(new_elem)%nodes), &
                           new_elem,nfnods
        end if
      else
        write (iout,201) new_elem
      end if
      !
    end do csc_side_loop
    !
  end do csc_cell_loop
  !
  call move_alloc( from=temp_nozzle , to=nozzle )
  !
  ! Deallocate the local arrays before leaving
  !
  if (allocated(noc)) then
    deallocate ( noc , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"noc",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(noc_ptr)) then
    deallocate ( noc_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"noc_ptr",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(csc)) then
    deallocate ( csc , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"csc",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(csc_ptr)) then
    deallocate ( csc_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"csc_ptr",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  100 format ("The number of missing surrounding cells (",i0, &
              ") does not match the expected number of boundary faces (", &
              i0,") that need to be created!!")
  !
  201 format (" The nodes array for new face element ",i0, &
              " is not allocated!!!")
  202 format (" The size of the nodes array (",i0, &
              ") for new face element ",i0, &
              " does not match the expected size (",i0,") !!!")
  203 format (" The face nodes for the new face element ",i0, &
              " do not match the expected face nodes!!!")
  204 format ("   The extracted face nodes = [ ",100(i0,:,", "))
 !204 format ("   The extracted face nodes = [ ",*(i0,:,", "))
  205 format ("   The expected  face nodes = [ ",100(i0,:,", "))
 !205 format ("   The expected  face nodes = [ ",*(i0,:,", "))
  !
end subroutine create_nozzle_outflow_faces
!
!###############################################################################
!
subroutine find_elem_for_initial_jet_profile
  !
  !.. Use Statements ..
  use geovar, only : grid
  !
  !.. Local Scalars ..
  integer :: l,n,init
  !
  !.. Local Arrays ..
  real(wp) :: xyz_ave(1:3)
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: pts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_elem_for_initial_jet_profile"
  !
continue
  !
  xyz_ave(:) = zero
  !
  do n = 1,size(grid%elem)
    !
    if (grid%elem(n)%tags(1) == 0) then
      !
      pts = grid%elem(n)%pts(:) ! USING F2003 AUTO-REALLOCATION
      !
      do l = 1,size(grid%xyz,dim=1)
        xyz_ave(l) = sum(grid%xyz(l,pts(:))) / real(size(pts),kind=wp)
      end do
      !
      if (xyz_ave(1) > zero) then
        init = 1
      else
        init = merge(1,0,point_is_in_smc000_nozzle(xyz_ave))
      end if
      !
      if (size(grid%elem(n)%tags) < 2) then
        grid%elem(n)%tags = [grid%elem(n)%tags(1),init] ! F2003 AUTO-REALLOC
      else
        grid%elem(n)%tags(2) = init
      end if
      !
    end if
    !
  end do
  !
end subroutine find_elem_for_initial_jet_profile
!
!###############################################################################
!
subroutine create_cgns_base_for_nozzle(nozzle,nozzle_grps,b)
  !
  !.. Formal Arguments ..
  type(grid_t),                    intent(in) :: nozzle
  type(cell_grps_t),               intent(in) :: nozzle_grps(:)
  type(cg_base_t), allocatable, intent(inout) :: b
  !
  !.. Local Scalars ..
  integer :: n,n1,n2,l,l1,l2,ierr
  integer :: ncells,npts,ncon
  integer :: this_cell,this_type
  integer :: mintype,maxtype
  !
  character(len=CGLEN) :: grp_name
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_cgns_base_for_nozzle"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Create the CGNS base
  !
  allocate ( b , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"b",1,__LINE__,__FILE__,ierr,error_message)
  !
  b%name = "Base"
  b%cell_dim = 3
  b%phys_dim = 3
  !
  ! Create the CGNS zone
  !
  allocate ( b%zone(1:1) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"b%zone",1,__LINE__,__FILE__,ierr,error_message)
  !
  b%zone(1)%name = "SMC Nozzle Interior"
  b%zone(1)%nvertex = size(nozzle%xyz,dim=2)
  b%zone(1)%ncell = size(nozzle_grps(1)%cells)
  b%zone(1)%nboundvertex = 0
  b%zone(1)%type = UNSTRUCTURED
  !
  ! Create the GridCoordinates_t node for the zone
  !
  allocate ( b%zone(1)%grid(1:1) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"b%zone(1)%grid",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  allocate ( b%zone(1)%grid(1)%coord(1:size(nozzle%xyz,dim=1)) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"b%zone(1)%grid(1)%coord",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  b%zone(1)%grid(1)%coord(1)%name = "CoordinateX"
  b%zone(1)%grid(1)%coord(2)%name = "CoordinateY"
  b%zone(1)%grid(1)%coord(3)%name = "CoordinateZ"
  !
  b%zone(1)%grid(1)%coord(1)%x = nozzle%xyz(1,:) ! USING F2003 AUTO-REALLOCATION
  b%zone(1)%grid(1)%coord(2)%x = nozzle%xyz(2,:) ! USING F2003 AUTO-REALLOCATION
  b%zone(1)%grid(1)%coord(3)%x = nozzle%xyz(3,:) ! USING F2003 AUTO-REALLOCATION
  !
  ! Create an Elements_t node for each of the nozzle groups
  !
  allocate ( b%zone(1)%section(1:size(nozzle_grps)) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"b%zone(1)%section",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  ! Also create BC_t nodes for each of the boundary condition groups
  !
  allocate ( b%zone(1)%bc(1:size(nozzle_grps)-1) , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"b%zone(1)%bc",1,__LINE__,__FILE__,ierr,error_message)
  !
  n2 = 0
  do n = 1,size(nozzle_grps)
    !
    grp_name = trim(adjustl(nozzle_grps(n)%name))
    !
    ncells = size(nozzle_grps(n)%cells)
    !
    write (iout,101) trim(adjustl(grp_name)),ncells
    101 format (5x,"Number of elements in group '",a,"' = ",i0)
    !
    n1 = n2 + 1
    n2 = n2 + ncells
    !
    if (uppercase(grp_name) == "INTERIOR") then
      b%zone(1)%section(n)%name = "Hex_Elements"
    else
      b%zone(1)%section(n)%name = "quad_" // grp_name
    end if
    b%zone(1)%section(n)%start = n1
    b%zone(1)%section(n)%end = n2
    b%zone(1)%section(n)%nbndry = 0
    !
    mintype = minval(nozzle%elem( nozzle_grps(n)%cells(:) )%prop%elem_type)
    maxtype = maxval(nozzle%elem( nozzle_grps(n)%cells(:) )%prop%elem_type)
    !
    if (mintype == maxtype) then
      !
      b%zone(1)%section(n)%type = mintype
      !
      npts = size(nozzle%elem(nozzle_grps(n)%cells(1))%pts)
      ncon = npts*ncells
      !
      allocate ( b%zone(1)%section(n)%elements(1:ncon) , source=0_CST , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) n
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      l2 = 0
      do l = 1,ncells
        !
        this_cell = nozzle_grps(n)%cells(l)
        !
        l1 = l2 + 1
        l2 = l2 + npts
        !
        b%zone(1)%section(n)%elements(l1:l2) = nozzle%elem(this_cell)%pts(:)
        !
      end do
      !
    else
      !
      b%zone(1)%section(n)%type = MIXED
      !
      ncon = 0
      do l = 1,ncells
        ncon = ncon + 1 + size(nozzle%elem( nozzle_grps(n)%cells(l) )%pts)
      end do
      !
      allocate ( b%zone(1)%section(n)%elements(1:ncon) , source=0_CST , &
                 stat=ierr , errmsg=error_message )
      write (array_name,1) n
      call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
      !
      l2 = 0
      do l = 1,ncells
        !
        this_cell = nozzle_grps(n)%cells(l)
        this_type = nozzle%elem(this_cell)%prop%elem_type
        !
        l1 = l2 + 1
        l2 = l2 + 1 + size(nozzle%elem(this_cell)%pts)
        !
        b%zone(1)%section(n)%elements(l1:l2) = [this_type, &
                                                nozzle%elem(this_cell)%pts(:)]
        !
      end do
      !
    end if
    !
    if (n > 1) then
      !
      grp_name = trim(adjustl(nozzle_grps(n)%name))
      !
      b%zone(1)%bc(n-1)%name = grp_name
      b%zone(1)%bc(n-1)%ptset_type = PointRange
      b%zone(1)%bc(n-1)%npnts = 2
      b%zone(1)%bc(n-1)%pnts = [n1,n2] ! USING F2003 AUTO-REALLOCATION
      b%zone(1)%bc(n-1)%location = FaceCenter
      !
      if (uppercase(grp_name) == "NOZZLEINTERIOR") then
        b%zone(1)%bc(n-1)%boco_type = BCWall
      else if (uppercase(grp_name) == "NOZZLEOUTFLOW") then
        b%zone(1)%bc(n-1)%boco_type = BCOutflowSubsonic
      else if (uppercase(grp_name) == "INFLOWSUBSONIC") then
        b%zone(1)%bc(n-1)%boco_type = BCInflowSubsonic
      else
        b%zone(1)%bc(n-1)%boco_type = BCTypeNull
      end if
      !
    end if
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ("b%zone(1)%section(",i0,")%elements")
  !
end subroutine create_cgns_base_for_nozzle
!
!###############################################################################
!
subroutine write_cgns_gridfile(cgns_grid_file,b)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: cgns_grid_file
  type(cg_base_t),  intent(in) :: b ! CGNS base to write
  !
  !.. Local Scalars ..
  integer :: n,ierr
  !
  integer(CBT) :: if_n,ib_n,iz_n,is_n,ibc_n
  !
  !.. Local Arrays ..
  integer(CST) :: zone_size(1:3)
  integer(CBT) :: ixyz_n(1:3)
  !
  !.. Local Allocatable Arrays ..
 !real(CRT), allocatable :: var(:)
  real(r8), allocatable :: var(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_cgns_gridfile"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ibc_n = 0
  !
  ! Open the CGNS grid file
  !
 !write (iout,128) cgns_grid_file
 !128 format ("cgns_grid_file = '",a,"'")
  call cgp_open_f(trim(adjustl(cgns_grid_file)),CG_MODE_WRITE,if_n,cgierr)
  call cgns_error(pname,cgns_grid_file,"cg_open_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  ! Set the CGNS base (who knows what this is for)
  !
  call cg_base_write_f(if_n,trim(adjustl(b%name)), &
                       b%cell_dim,b%phys_dim,ib_n,cgierr)
  call cgns_error(pname,cgns_grid_file,"cg_base_write_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  ! Write the title for this solution
  !
  zone_size(1) = b%zone(1)%nvertex
  zone_size(2) = b%zone(1)%ncell
  zone_size(3) = b%zone(1)%nboundvertex
  !
  call cg_zone_write_f(if_n,ib_n,trim(adjustl(b%zone(1)%name)), &
                       zone_size,b%zone(1)%type,iz_n,cgierr)
  call cgns_error(pname,cgns_grid_file,"cg_zone_write_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  ! Create data nodes for the coordinates
  !
  do n = 1,size(b%zone(1)%grid(1)%coord)
    ! USING F2003 AUTO-REALLOCATION
    var = real( b%zone(1)%grid(1)%coord(n)%x(:) , kind(var) )
   !call cg_coord_write_f(if_n,ib_n,iz_n,CGNS_KIND_REAL, &
    call cg_coord_write_f(if_n,ib_n,iz_n,RealDouble, &
                          trim(adjustl(b%zone(1)%grid(1)%coord(n)%name)), &
                          var,ixyz_n(n),cgierr)
    call cgns_error(pname,cgns_grid_file,"cg_coord_write_f",cgierr, &
                    __LINE__,__FILE__,"Write node: Coordinate",n, &
                    not_parallel=true)
  end do
  !
  ! Deallocate the arrays var
  !
  if (allocated(var)) then
    deallocate ( var , stat=ierr )
    call alloc_error(pname,"var",2,__LINE__,__FILE__,ierr)
  end if
  !
  ! Write the grid connectivity for the interior cells
  !
  do n = 1,size(b%zone(1)%section)
    call cg_section_write_f(if_n,ib_n,iz_n, &
               trim(adjustl(b%zone(1)%section(n)%name)), &
                            b%zone(1)%section(n)%type, &
                            b%zone(1)%section(n)%start, &
                            b%zone(1)%section(n)%end, &
                            b%zone(1)%section(n)%nbndry, &
                            b%zone(1)%section(n)%elements, &
                            is_n,cgierr)
    call cgns_error(pname,cgns_grid_file,"cg_section_write_f",cgierr, &
                    __LINE__,__FILE__,"Writing section # ",n, &
                    not_parallel=true)
  end do
  !
  do n = 1,size(b%zone(1)%bc)
    !
    call cg_boco_write_f(if_n,ib_n,iz_n, &
            trim(adjustl(b%zone(1)%bc(n)%name)), &
                         b%zone(1)%bc(n)%boco_type, &
                         b%zone(1)%bc(n)%ptset_type, &
                         b%zone(1)%bc(n)%npnts, &
                         b%zone(1)%bc(n)%pnts, &
                         ibc_n,cgierr)
    call cgns_error(pname,cgns_grid_file,"cg_boco_write_f",cgierr, &
                    __LINE__,__FILE__,"Writing boundary condition #",n, &
                    not_parallel=true)
    !
    call cg_boco_gridlocation_write_f(if_n,ib_n,iz_n,ibc_n, &
                                      b%zone(1)%bc(n)%location,cgierr)
    call cgns_error(pname,cgns_grid_file,"cg_boco_gridlocation_write_f", &
                    cgierr,__LINE__,__FILE__,"Writing grid location for bc #", &
                    n,not_parallel=true)
    !
  end do
  !
  ! Close the CGNS FILE
  !
  call cgp_close_f(if_n,cgierr)
  call cgns_error(pname,cgns_grid_file,"cg_close_f",cgierr, &
                  __LINE__,__FILE__,not_parallel=true)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine write_cgns_gridfile
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
  use ovar,      only : loc_solution_pts,loc_flux_pts,bc_names
  !
  !.. Local Scalars ..
  integer :: n,nocsz,ierr,nbcnm
  integer(INT_MPI) :: isz
  !
  !.. Local Arrays ..
  integer :: grid_size(1:6)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "broadcast_solver_geom_arrays"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Broadcast the grid dimensions
  !
  if (mypnum == glb_root) then
    grid_size(:) = [ nr , nnode , ncell , nfbnd , &
                     size(nodes_of_cell) , size(bc_names) ]
  end if
  !
  call mpi_bcast(grid_size,size(grid_size,kind=int_mpi),mpi_inttyp, &
                 glb_root,MPI_COMM_WORLD,mpierr)
  if (debug_cgns) write (iout,1) mypnum,(grid_size(n),n=1,size(grid_size))
  !
  nr    = grid_size(1)
  nnode = grid_size(2)
  ncell = grid_size(3)
  nfbnd = grid_size(4)
  nocsz = grid_size(5)
  nbcnm = grid_size(6)
  !
  ! Have all the non-root processors allocate the arrays
  !
  if (mypnum /= glb_root) then
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
    allocate ( bc_names(1:nbcnm) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bc_names",1,__LINE__,__FILE__,ierr,error_message)
    !
  end if
  !
  ! Broadcast the bface array
  !
  isz = int( nbfai*nfbnd , kind=int_mpi )
  call mpi_bcast(bface,isz,mpi_inttyp,glb_root,MPI_COMM_WORLD,mpierr)
  !
  ! Broadcast the nodes_of_cell_ptr array
  !
  isz = int( ncell+1 , kind=int_mpi )
  call mpi_bcast(nodes_of_cell_ptr,isz,mpi_inttyp, &
                 glb_root,MPI_COMM_WORLD,mpierr)
  !
  ! Broadcast the nodes_of_cell array
  !
  isz = int( nocsz , kind=int_mpi )
  call mpi_bcast(nodes_of_cell,isz,mpi_inttyp,glb_root,MPI_COMM_WORLD,mpierr)
  !
  ! Broadcast the xyz_nodes array
  !
  isz = int( nr*nnode , kind=int_mpi )
  call mpi_bcast(xyz_nodes,isz,mpi_flttyp,glb_root,MPI_COMM_WORLD,mpierr)
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
  do n = 1,nbcnm
    isz = len(bc_names(n),kind=int_mpi)
    call mpi_bcast(bc_names(n),isz,MPI_CHARACTER,glb_root,MPI_COMM_WORLD,mpierr)
  end do
  !
  ! Make sure that Lobatto nodes are requested if the grid contains triangle
  ! elements since Gauss nodes is not currently set up for triangles.
  !
 !if (any(cell_geom == Geom_Tria)) then
 !  write (error_message,11)
 !  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
 !end if
  !
  if (loc_solution_pts /= Legendre_Gauss_Lobatto .or. &
      loc_flux_pts /= Legendre_Gauss_Lobatto) then
    if (any(cell_geom == Geom_Tria)) then
      write (error_message,12)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
    end if
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (" CPU #",i0,": grid_size(:) = ",5(1x,i0),:,/,27x,100(5(1x,i0),/))
 !1 format (" CPU #",i0,": grid_size(:) = ",5(1x,i0),:,/,27x,*(5(1x,i0),/))
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
  use ovar,      only : loc_solution_pts,loc_flux_pts
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
 !  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
 !end if
  !
  if (loc_solution_pts /= Legendre_Gauss_Lobatto .or. &
      loc_flux_pts /= Legendre_Gauss_Lobatto) then
    if (any(cell_geom == Geom_Tria)) then
      write (error_message,12)
      call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
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
pure function get_element_type(elem_geom,elem_order,elem_family) &
                        result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: elem_geom
  integer, intent(in) :: elem_order
  !
  !.. Optional Arguments ..
  type(geom_family_t), optional, intent(in) :: elem_family
  !
  !.. Function Result ..
  integer(CET) :: return_value
  !
  !.. Local Scalars ..
  type(geom_family_t) :: family
  !
continue
  !
  return_value = ElementTypeNull
  !
  family = Complete
  if (present(elem_family)) then
    family = elem_family
  end if
  !
  select case (elem_geom)
    !
    case (Geom_Node)
      !
      return_value = NODE
      !
    case (Geom_Edge)
      !
      select case (elem_order)
        case (1)
            return_value = BAR_2
        case (2)
            return_value = BAR_3
        case (3)
            return_value = BAR_4
        case (4)
            return_value = BAR_5
      end select
      !
    case (Geom_Tria)
      !
      select case (elem_order)
        case (1)
            return_value = TRI_3
        case (2)
            return_value = TRI_6
        case (3)
            if (family == Edge_Defined) then
              return_value = TRI_9
            else
              return_value = TRI_10
            end if
        case (4)
            if (family == Edge_Defined) then
              return_value = TRI_12
            else
              return_value = TRI_15
            end if
      end select
      !
    case (Geom_Quad)
      !
      select case (elem_order)
        case (1)
            return_value = QUAD_4
        case (2)
            if (family == Edge_Defined) then
              return_value = QUAD_8
            else
              return_value = QUAD_9
            end if
        case (3)
            if (family == Edge_Defined) then
              return_value = QUAD_12
            else
              return_value = QUAD_16
            end if
        case (4)
            if (family == Edge_Defined) then
              return_value = QUAD_P4_16
            else
              return_value = QUAD_25
            end if
      end select
      !
    case (Geom_Tetr)
      !
      select case (elem_order)
        case (1)
            return_value = TETRA_4
        case (2)
            return_value = TETRA_10
        case (3)
            if (family == Edge_Defined) then
              return_value = TETRA_16
            else
              return_value = TETRA_20
            end if
        case (4)
            if (family == Edge_Defined) then
              return_value = TETRA_22
            else if (family == Face_Defined) then
              return_value = TETRA_34
            else
              return_value = TETRA_35
            end if
      end select
      !
    case (Geom_Pyra)
      !
      select case (elem_order)
        case (1)
            return_value = PYRA_5
        case (2)
            if (family == Edge_Defined) then
              return_value = PYRA_13
            else
              return_value = PYRA_14
            end if
        case (3)
            if (family == Edge_Defined) then
              return_value = PYRA_21
            else if (family == Face_Defined) then
              return_value = PYRA_29
            else
              return_value = PYRA_30
            end if
        case (4)
            if (family == Edge_Defined) then
              return_value = PYRA_P4_29
            else if (family == Face_Defined) then
              return_value = PYRA_50
            else
              return_value = PYRA_55
            end if
      end select
      !
    case (Geom_Pris)
      !
      select case (elem_order)
        case (1)
            return_value = PENTA_6
        case (2)
            if (family == Edge_Defined) then
              return_value = PENTA_15
            else
              return_value = PENTA_18
            end if
        case (3)
            if (family == Edge_Defined) then
              return_value = PENTA_24
            else if (family == Face_Defined) then
              return_value = PENTA_38
            else
              return_value = PENTA_40
            end if
        case (4)
            if (family == Edge_Defined) then
              return_value = PENTA_33
            else if (family == Face_Defined) then
              return_value = PENTA_66
            else
              return_value = PENTA_75
            end if
      end select
      !
    case (Geom_Hexa)
      !
      select case (elem_order)
        case (1)
            return_value = HEXA_8
        case (2)
            if (family == Edge_Defined) then
              return_value = HEXA_20
            else
              return_value = HEXA_27
            end if
        case (3)
            if (family == Edge_Defined) then
              return_value = HEXA_32
            else if (family == Face_Defined) then
              return_value = HEXA_56
            else
              return_value = HEXA_64
            end if
        case (4)
            if (family == Edge_Defined) then
              return_value = HEXA_44
            else if (family == Face_Defined) then
              return_value = HEXA_98
            else
              return_value = HEXA_125
            end if
      end select
      !
  end select
  !
end function get_element_type
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
     !return_value = prop_t(npts=1, &
     !                      order=1, &
     !                      geom=Geom_Node, &
     !                      elem_type=int(elem_type)
      return_value%npts = 1
      return_value%order = 1
      return_value%geom = Geom_Node
      return_value%elem_type = int(elem_type)
    case (BAR_2)
     !return_value = prop_t(npts=2, &
     !                      order=1, &
     !                      geom=Geom_Edge, &
     !                      elem_type=int(elem_type)
      return_value%npts = 2
      return_value%order = 1
      return_value%geom = Geom_Edge
      return_value%elem_type = int(elem_type)
    case (BAR_3)
     !return_value = prop_t(npts=3, &
     !                      order=2, &
     !                      geom=Geom_Edge, &
     !                      elem_type=int(elem_type)
      return_value%npts = 3
      return_value%order = 2
      return_value%geom = Geom_Edge
      return_value%elem_type = int(elem_type)
    case (TRI_3)
     !return_value = prop_t(npts=3, &
     !                      order=1, &
     !                      geom=Geom_Tria, &
     !                      elem_type=int(elem_type)
      return_value%npts = 3
      return_value%order = 1
      return_value%geom = Geom_Tria
      return_value%elem_type = int(elem_type)
    case (TRI_6)
     !return_value = prop_t(npts=6, &
     !                      order=2, &
     !                      geom=Geom_Tria, &
     !                      elem_type=int(elem_type)
      return_value%npts = 6
      return_value%order = 2
      return_value%geom = Geom_Tria
      return_value%elem_type = int(elem_type)
    case (QUAD_4)
     !return_value = prop_t(npts=4, &
     !                      order=1, &
     !                      geom=Geom_Quad, &
     !                      elem_type=int(elem_type)
      return_value%npts = 4
      return_value%order = 1
      return_value%geom = Geom_Quad
      return_value%elem_type = int(elem_type)
    case (QUAD_8)
     !return_value = prop_t(npts=8, &
     !                      order=2, &
     !                      geom=Geom_Quad, &
     !                      elem_type=int(elem_type), &
     !                      family=Edge_Defined)
      return_value%npts = 8
      return_value%order = 2
      return_value%geom = Geom_Quad
      return_value%elem_type = int(elem_type)
      return_value%family = Edge_Defined
    case (QUAD_9)
     !return_value = prop_t(npts=9, &
     !                      order=2, &
     !                      geom=Geom_Quad, &
     !                      elem_type=int(elem_type)
      return_value%npts = 9
      return_value%order = 2
      return_value%geom = Geom_Quad
      return_value%elem_type = int(elem_type)
    case (TETRA_4)
     !return_value = prop_t(npts=4, &
     !                      order=1, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=int(elem_type)
      return_value%npts = 4
      return_value%order = 1
      return_value%geom = Geom_Tetr
      return_value%elem_type = int(elem_type)
    case (TETRA_10)
     !return_value = prop_t(npts=10, &
     !                      order=2, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=int(elem_type)
      return_value%npts = 10
      return_value%order = 2
      return_value%geom = Geom_Tetr
      return_value%elem_type = int(elem_type)
    case (PYRA_5)
     !return_value = prop_t(npts=5, &
     !                      order=1, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=int(elem_type)
      return_value%npts = 5
      return_value%order = 1
      return_value%geom = Geom_Pyra
      return_value%elem_type = int(elem_type)
    case (PYRA_13)
     !return_value = prop_t(npts=13, &
     !                      order=2, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=int(elem_type), &
     !                      family=Edge_Defined)
      return_value%npts = 13
      return_value%order = 2
      return_value%geom = Geom_Pyra
      return_value%elem_type = int(elem_type)
      return_value%family = Edge_Defined
    case (PYRA_14)
     !return_value = prop_t(npts=14, &
     !                      order=2, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=int(elem_type)
      return_value%npts = 14
      return_value%order = 2
      return_value%geom = Geom_Pyra
      return_value%elem_type = int(elem_type)
    case (PENTA_6)
     !return_value = prop_t(npts=6, &
     !                      order=1, &
     !                      geom=Geom_Pris, &
     !                      elem_type=int(elem_type)
      return_value%npts = 6
      return_value%order = 1
      return_value%geom = Geom_Pris
      return_value%elem_type = int(elem_type)
    case (PENTA_15)
     !return_value = prop_t(npts=15, &
     !                      order=2, &
     !                      geom=Geom_Pris, &
     !                      elem_type=int(elem_type), &
     !                      family=Edge_Defined)
      return_value%npts = 15
      return_value%order = 2
      return_value%geom = Geom_Pris
      return_value%elem_type = int(elem_type)
      return_value%family = Edge_Defined
    case (PENTA_18)
     !return_value = prop_t(npts=18, &
     !                      order=2, &
     !                      geom=Geom_Pris, &
     !                      elem_type=int(elem_type)
      return_value%npts = 18
      return_value%order = 2
      return_value%geom = Geom_Pris
      return_value%elem_type = int(elem_type)
    case (HEXA_8)
     !return_value = prop_t(npts=8, &
     !                      order=1, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=int(elem_type)
      return_value%npts = 8
      return_value%order = 1
      return_value%geom = Geom_Hexa
      return_value%elem_type = int(elem_type)
    case (HEXA_20)
     !return_value = prop_t(npts=20, &
     !                      order=2, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=int(elem_type), &
     !                      family=Edge_Defined)
      return_value%npts = 20
      return_value%order = 2
      return_value%geom = Geom_Hexa
      return_value%elem_type = int(elem_type)
      return_value%family = Edge_Defined
    case (HEXA_27)
     !return_value = prop_t(npts=27, &
     !                      order=2, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=int(elem_type)
      return_value%npts = 27
      return_value%order = 2
      return_value%geom = Geom_Hexa
      return_value%elem_type = int(elem_type)
   !case (MIXED)   ! MIXED TYPE SHOULD BE SEPARATED BY NOW
   !case (NGON_n)  ! UNSUPPORTED TYPE
   !case (NFACE_n) ! UNSUPPORTED TYPE
    case (BAR_4)
     !return_value = prop_t(npts=4, &
     !                      order=3, &
     !                      geom=Geom_Edge, &
     !                      elem_type=int(elem_type)
      return_value%npts = 4
      return_value%order = 3
      return_value%geom = Geom_Edge
      return_value%elem_type = int(elem_type)
    case (TRI_9)
     !return_value = prop_t(npts=9, &
     !                      order=3, &
     !                      geom=Geom_Tria, &
     !                      elem_type=int(elem_type), &
     !                      family=Edge_Defined)
      return_value%npts = 9
      return_value%order = 3
      return_value%geom = Geom_Tria
      return_value%elem_type = int(elem_type)
      return_value%family = Edge_Defined
    case (TRI_10)
     !return_value = prop_t(npts=10, &
     !                      order=3, &
     !                      geom=Geom_Tria, &
     !                      elem_type=int(elem_type)
      return_value%npts = 10
      return_value%order = 3
      return_value%geom = Geom_Tria
      return_value%elem_type = int(elem_type)
    case (QUAD_12)
     !return_value = prop_t(npts=12, &
     !                      order=3, &
     !                      geom=Geom_Quad, &
     !                      elem_type=int(elem_type), &
     !                      family=Edge_Defined)
      return_value%npts = 12
      return_value%order = 3
      return_value%geom = Geom_Quad
      return_value%elem_type = int(elem_type)
      return_value%family = Edge_Defined
    case (QUAD_16)
     !return_value = prop_t(npts=16, &
     !                      order=3, &
     !                      geom=Geom_Quad, &
     !                      elem_type=int(elem_type)
      return_value%npts = 16
      return_value%order = 3
      return_value%geom = Geom_Quad
      return_value%elem_type = int(elem_type)
    case (TETRA_16)
     !return_value = prop_t(npts=16, &
     !                      order=3, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=int(elem_type), &
     !                      family=Edge_Defined)
      return_value%npts = 16
      return_value%order = 3
      return_value%geom = Geom_Tetr
      return_value%elem_type = int(elem_type)
      return_value%family = Edge_Defined
    case (TETRA_20)
     !return_value = prop_t(npts=20, &
     !                      order=3, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=int(elem_type)
      return_value%npts = 20
      return_value%order = 3
      return_value%geom = Geom_Tetr
      return_value%elem_type = int(elem_type)
    case (PYRA_21)
     !return_value = prop_t(npts=21, &
     !                      order=3, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=int(elem_type), &
     !                      family=Edge_Defined)
      return_value%npts = 21
      return_value%order = 3
      return_value%geom = Geom_Pyra
      return_value%elem_type = int(elem_type)
      return_value%family = Edge_Defined
    case (PYRA_29)
     !return_value = prop_t(npts=29, &
     !                      order=3, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=int(elem_type), &
     !                      family=Face_Defined)
      return_value%npts = 29
      return_value%order = 3
      return_value%geom = Geom_Pyra
      return_value%elem_type = int(elem_type)
      return_value%family = Face_Defined
    case (PYRA_30)
     !return_value = prop_t(npts=30, &
     !                      order=3, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=int(elem_type)
      return_value%npts = 30
      return_value%order = 3
      return_value%geom = Geom_Pyra
      return_value%elem_type = int(elem_type)
    case (PENTA_24)
     !return_value = prop_t(npts=24, &
     !                      order=3, &
     !                      geom=Geom_Pris, &
     !                      elem_type=int(elem_type), &
     !                      family=Edge_Defined)
      return_value%npts = 24
      return_value%order = 3
      return_value%geom = Geom_Pris
      return_value%elem_type = int(elem_type)
      return_value%family = Edge_Defined
    case (PENTA_38)
     !return_value = prop_t(npts=38, &
     !                      order=3, &
     !                      geom=Geom_Pris, &
     !                      elem_type=int(elem_type), &
     !                      family=Face_Defined)
      return_value%npts = 38
      return_value%order = 3
      return_value%geom = Geom_Pris
      return_value%elem_type = int(elem_type)
      return_value%family = Face_Defined
    case (PENTA_40)
     !return_value = prop_t(npts=40, &
     !                      order=3, &
     !                      geom=Geom_Pris, &
     !                      elem_type=int(elem_type)
      return_value%npts = 40
      return_value%order = 3
      return_value%geom = Geom_Pris
      return_value%elem_type = int(elem_type)
    case (HEXA_32)
     !return_value = prop_t(npts=32, &
     !                      order=3, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=int(elem_type), &
     !                      family=Edge_Defined)
      return_value%npts = 32
      return_value%order = 3
      return_value%geom = Geom_Hexa
      return_value%elem_type = int(elem_type)
      return_value%family = Edge_Defined
    case (HEXA_56)
     !return_value = prop_t(npts=56, &
     !                      order=3, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=int(elem_type), &
     !                      family=Face_Defined)
      return_value%npts = 56
      return_value%order = 3
      return_value%geom = Geom_Hexa
      return_value%elem_type = int(elem_type)
      return_value%family = Face_Defined
    case (HEXA_64)
     !return_value = prop_t(npts=64, &
     !                      order=3, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=int(elem_type)
      return_value%npts = 64
      return_value%order = 3
      return_value%geom = Geom_Hexa
      return_value%elem_type = int(elem_type)
    case (BAR_5)
     !return_value = prop_t(npts=5, &
     !                      order=4, &
     !                      geom=Geom_Edge, &
     !                      elem_type=int(elem_type)
      return_value%npts = 5
      return_value%order = 4
      return_value%geom = Geom_Edge
      return_value%elem_type = int(elem_type)
    case (TRI_12)
     !return_value = prop_t(npts=12, &
     !                      order=4, &
     !                      geom=Geom_Tria, &
     !                      elem_type=int(elem_type), &
     !                      family=Edge_Defined)
      return_value%npts = 12
      return_value%order = 4
      return_value%geom = Geom_Tria
      return_value%elem_type = int(elem_type)
      return_value%family = Edge_Defined
    case (TRI_15)
     !return_value = prop_t(npts=15, &
     !                      order=4, &
     !                      geom=Geom_Tria, &
     !                      elem_type=int(elem_type)
      return_value%npts = 15
      return_value%order = 4
      return_value%geom = Geom_Tria
      return_value%elem_type = int(elem_type)
    case (QUAD_P4_16)
     !return_value = prop_t(npts=16, &
     !                      order=4, &
     !                      geom=Geom_Quad, &
     !                      elem_type=int(elem_type), &
     !                      family=Edge_Defined)
      return_value%npts = 16
      return_value%order = 4
      return_value%geom = Geom_Quad
      return_value%elem_type = int(elem_type)
      return_value%family = Edge_Defined
    case (QUAD_25)
     !return_value = prop_t(npts=25, &
     !                      order=4, &
     !                      geom=Geom_Quad, &
     !                      elem_type=int(elem_type)
      return_value%npts = 25
      return_value%order = 4
      return_value%geom = Geom_Quad
      return_value%elem_type = int(elem_type)
    case (TETRA_22)
     !return_value = prop_t(npts=22, &
     !                      order=4, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=int(elem_type), &
     !                      family=Edge_Defined)
      return_value%npts = 22
      return_value%order = 4
      return_value%geom = Geom_Tetr
      return_value%elem_type = int(elem_type)
      return_value%family = Edge_Defined
    case (TETRA_34)
     !return_value = prop_t(npts=34, &
     !                      order=4, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=int(elem_type), &
     !                      family=Face_Defined)
      return_value%npts = 34
      return_value%order = 4
      return_value%geom = Geom_Tetr
      return_value%elem_type = int(elem_type)
      return_value%family = Face_Defined
    case (TETRA_35)
     !return_value = prop_t(npts=35, &
     !                      order=4, &
     !                      geom=Geom_Tetr, &
     !                      elem_type=int(elem_type)
      return_value%npts = 35
      return_value%order = 4
      return_value%geom = Geom_Tetr
      return_value%elem_type = int(elem_type)
    case (PYRA_P4_29)
     !return_value = prop_t(npts=29, &
     !                      order=4, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=int(elem_type), &
     !                      family=Edge_Defined)
      return_value%npts = 29
      return_value%order = 4
      return_value%geom = Geom_Pyra
      return_value%elem_type = int(elem_type)
      return_value%family = Edge_Defined
    case (PYRA_50)
     !return_value = prop_t(npts=50, &
     !                      order=4, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=int(elem_type), &
     !                      family=Face_Defined)
      return_value%npts = 50
      return_value%order = 4
      return_value%geom = Geom_Pyra
      return_value%elem_type = int(elem_type)
      return_value%family = Face_Defined
    case (PYRA_55)
     !return_value = prop_t(npts=55, &
     !                      order=4, &
     !                      geom=Geom_Pyra, &
     !                      elem_type=int(elem_type)
      return_value%npts = 55
      return_value%order = 4
      return_value%geom = Geom_Pyra
      return_value%elem_type = int(elem_type)
    case (PENTA_33)
     !return_value = prop_t(npts=33, &
     !                      order=4, &
     !                      geom=Geom_Pris, &
     !                      elem_type=int(elem_type), &
     !                      family=Edge_Defined)
      return_value%npts = 33
      return_value%order = 4
      return_value%geom = Geom_Pris
      return_value%elem_type = int(elem_type)
      return_value%family = Edge_Defined
    case (PENTA_66)
     !return_value = prop_t(npts=66, &
     !                      order=4, &
     !                      geom=Geom_Pris, &
     !                      elem_type=int(elem_type), &
     !                      family=Face_Defined)
      return_value%npts = 66
      return_value%order = 4
      return_value%geom = Geom_Pris
      return_value%elem_type = int(elem_type)
      return_value%family = Face_Defined
    case (PENTA_75)
     !return_value = prop_t(npts=75, &
     !                      order=4, &
     !                      geom=Geom_Pris, &
     !                      elem_type=int(elem_type)
      return_value%npts = 75
      return_value%order = 4
      return_value%geom = Geom_Pris
      return_value%elem_type = int(elem_type)
    case (HEXA_44)
     !return_value = prop_t(npts=44, &
     !                      order=4, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=int(elem_type), &
     !                      family=Edge_Defined)
      return_value%npts = 44
      return_value%order = 4
      return_value%geom = Geom_Hexa
      return_value%elem_type = int(elem_type)
      return_value%family = Edge_Defined
    case (HEXA_98)
     !return_value = prop_t(npts=98, &
     !                      order=4, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=int(elem_type), &
     !                      family=Face_Defined)
      return_value%npts = 98
      return_value%order = 4
      return_value%geom = Geom_Hexa
      return_value%elem_type = int(elem_type)
      return_value%family = Face_Defined
    case (HEXA_125)
     !return_value = prop_t(npts=125, &
     !                      order=4, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=int(elem_type)
      return_value%npts = 125
      return_value%order = 4
      return_value%geom = Geom_Hexa
      return_value%elem_type = int(elem_type)
    case default
     !return_value = prop_t(npts=-1, &
     !                      order=-1, &
     !                      geom=Geom_Unknown, &
     !                      elem_type=int(elem_type)
      return_value%npts = -1
      return_value%order = -1
      return_value%geom = Geom_Unknown
      return_value%elem_type = int(elem_type)
  end select
  !
end function element_properties
!
!###############################################################################
!
elemental function reorder_cgns_pts(elem,reverse) result(return_value)
  !
  !.. Formal Arguments ..
  type(elem_t), intent(in) :: elem
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: reverse
  !
  !.. Function Result ..
  type(elem_t) :: return_value
  !
continue
  !
  return_value = elem
  !
  select case (elem%prop%geom)
    case (Geom_Edge)
      return_value%pts = reorder_cgns_edge_pts(elem%pts,reverse)
    case (Geom_Tria)
      return_value%pts = reorder_cgns_tria_pts(elem%pts,elem%prop%elem_type, &
                                               reverse)
    case (Geom_Quad)
      return_value%pts = reorder_cgns_quad_pts(elem%pts,elem%prop%elem_type, &
                                               reverse)
    case (Geom_Tetr)
      return_value%pts = reorder_cgns_tetr_pts(elem%pts,elem%prop%elem_type, &
                                               reverse)
    case (Geom_Pyra)
      return_value%pts = reorder_cgns_pyra_pts(elem%pts,elem%prop%elem_type, &
                                               reverse)
    case (Geom_Pris)
      return_value%pts = reorder_cgns_pris_pts(elem%pts,elem%prop%elem_type, &
                                               reverse)
    case (Geom_Hexa)
      return_value%pts = reorder_cgns_hexa_pts(elem%pts,elem%prop%elem_type, &
                                               reverse)
  end select
  !
end function reorder_cgns_pts
!
!###############################################################################
!
pure function reorder_cgns_edge_pts(pts,reverse) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: reverse
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Scalars ..
  integer     :: n,i
  logical(lk) :: reorder_to_cgns
  !
continue
  !
  reorder_to_cgns = fals
  if (present(reverse)) then
    reorder_to_cgns = reverse
  end if
  !
  n = size(pts)
  !
  if (reorder_to_cgns) then
   !return_value([1,n,(i+1,i=2,n-1)]) = pts(:)
    !
    return_value(1) = pts(1)
    return_value(2) = pts(n)
    !
    do i = 2,n-1
      return_value(i+1) = pts(i)
    end do
    !
  else
   !return_value(:) = pts([1,(i+1,i=2,n-1),2])
    !
    return_value(1) = pts(1)
    return_value(n) = pts(2)
    !
    do i = 2,n-1
      return_value(i) = pts(i+1)
    end do
    !
  end if
  !
end function reorder_cgns_edge_pts
!
!###############################################################################
!
pure function reorder_cgns_tria_pts(pts,elem_type,reverse) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  integer, intent(in) :: elem_type
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: reverse
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Scalars ..
  logical(lk) :: reorder_to_cgns
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
  reorder_to_cgns = fals
  if (present(reverse)) then
    reorder_to_cgns = reverse
  end if
  !
  if (reorder_to_cgns) then
    !
    select case (elem_type)
      case (TRI_3)
        return_value( tria3(:) ) = pts(:)
      case (TRI_6)
        return_value( tria6(:) ) = pts(:)
      case (TRI_9)
        return_value( tria9(:) ) = pts(:)
      case (TRI_10)
        return_value( tria10(:) ) = pts(:)
      case (TRI_12)
        return_value( tria12(:) ) = pts(:)
      case (TRI_15)
        return_value( tria15(:) ) = pts(:)
      case default
        return_value = pts(:)
    end select
    !
  else
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
  end if
  !
end function reorder_cgns_tria_pts
!
!###############################################################################
!
pure function reorder_cgns_quad_pts(pts,elem_type,reverse) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  integer, intent(in) :: elem_type
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: reverse
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Scalars ..
  logical(lk) :: reorder_to_cgns
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
  reorder_to_cgns = fals
  if (present(reverse)) then
    reorder_to_cgns = reverse
  end if
  !
  if (reorder_to_cgns) then
    !
    select case (elem_type)
      case (QUAD_4)
        return_value( quad4(:) ) = pts(:)
      case (QUAD_8)
        return_value( quad8(:) ) = pts(:)
      case (QUAD_9)
        return_value( quad9(:) ) = pts(:)
      case (QUAD_12)
        return_value( quad12(:) ) = pts(:)
      case (QUAD_16)
        return_value( quad3_16(:) ) = pts(:)
      case (QUAD_P4_16)
        return_value( quad4_16(:) ) = pts(:)
      case (QUAD_25)
        return_value( quad25(:) ) = pts(:)
      case default
        return_value = pts(:)
    end select
    !
  else
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
  end if
  !
end function reorder_cgns_quad_pts
!
!###############################################################################
!
pure function reorder_cgns_tetr_pts(pts,elem_type,reverse) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  integer, intent(in) :: elem_type
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: reverse
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Scalars ..
  logical(lk) :: reorder_to_cgns
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
  reorder_to_cgns = fals
  if (present(reverse)) then
    reorder_to_cgns = reverse
  end if
  !
  if (reorder_to_cgns) then
    !
    select case (elem_type)
      case (TETRA_4)
        return_value( tetr4(:) ) = pts(:)
      case (TETRA_10)
        return_value( tetr10(:) ) = pts(:)
      case (TETRA_16)
        return_value( tetr16(:) ) = pts(:)
      case (TETRA_20)
        return_value( tetr20(:) ) = pts(:)
      case (TETRA_22)
        return_value( tetr22(:) ) = pts(:)
      case (TETRA_34)
        return_value( tetr34(:) ) = pts(:)
      case (TETRA_35)
        return_value( tetr35(:) ) = pts(:)
      case default
        return_value = pts(:)
    end select
    !
  else
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
  end if
  !
end function reorder_cgns_tetr_pts
!
!###############################################################################
!
pure function reorder_cgns_pyra_pts(pts,elem_type,reverse) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  integer, intent(in) :: elem_type
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: reverse
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Scalars ..
  logical(lk) :: reorder_to_cgns
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
  reorder_to_cgns = fals
  if (present(reverse)) then
    reorder_to_cgns = reverse
  end if
  !
  if (reorder_to_cgns) then
    !
    select case (elem_type)
      case (PYRA_5)
        return_value( pyra5(:) ) = pts(:)
      case (PYRA_13)
        return_value( pyra13(:) ) = pts(:)
      case (PYRA_14)
        return_value( pyra14(:) ) = pts(:)
      case (PYRA_21)
        return_value( pyra21(:) ) = pts(:)
      case (PYRA_29)
        return_value( pyra3_29(:) ) = pts(:)
      case (PYRA_P4_29)
        return_value( pyra4_29(:) ) = pts(:)
      case (PYRA_30)
        return_value( pyra30(:) ) = pts(:)
      case (PYRA_50)
        return_value( pyra50(:) ) = pts(:)
      case (PYRA_55)
        return_value( pyra55(:) ) = pts(:)
      case default
        return_value = pts(:)
    end select
    !
  else
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
  end if
  !
end function reorder_cgns_pyra_pts
!
!###############################################################################
!
pure function reorder_cgns_pris_pts(pts,elem_type,reverse) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  integer, intent(in) :: elem_type
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: reverse
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Scalars ..
  logical(lk) :: reorder_to_cgns
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
  reorder_to_cgns = fals
  if (present(reverse)) then
    reorder_to_cgns = reverse
  end if
  !
  if (reorder_to_cgns) then
    !
    select case (elem_type)
      case (PENTA_6)
        return_value( pris6(:) ) = pts(:)
      case (PENTA_15)
        return_value( pris15(:) ) = pts(:)
      case (PENTA_18)
        return_value( pris18(:) ) = pts(:)
      case (PENTA_24)
        return_value( pris24(:) ) = pts(:)
      case (PENTA_33)
        return_value( pris33(:) ) = pts(:)
      case (PENTA_38)
        return_value( pris38(:) ) = pts(:)
      case (PENTA_40)
        return_value( pris40(:) ) = pts(:)
      case (PENTA_66)
        return_value( pris66(:) ) = pts(:)
      case (PENTA_75)
        return_value( pris75(:) ) = pts(:)
      case default
        return_value = pts(:)
    end select
    !
  else
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
  end if
  !
end function reorder_cgns_pris_pts
!
!###############################################################################
!
pure function reorder_cgns_hexa_pts(pts,elem_type,reverse) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  integer, intent(in) :: elem_type
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: reverse
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
  !.. Local Scalars ..
  logical(lk) :: reorder_to_cgns
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
  reorder_to_cgns = fals
  if (present(reverse)) then
    reorder_to_cgns = reverse
  end if
  !
  if (reorder_to_cgns) then
    !
    select case (elem_type)
      case (HEXA_8)
        return_value( hexa8(:) ) = pts(:)
      case (HEXA_20)
        return_value( hexa20(:) ) = pts(:)
      case (HEXA_27)
        return_value( hexa27(:) ) = pts(:)
      case (HEXA_32)
        return_value( hexa32(:) ) = pts(:)
      case (HEXA_44)
        return_value( hexa44(:) ) = pts(:)
      case (HEXA_56)
        return_value( hexa56(:) ) = pts(:)
      case (HEXA_64)
        return_value( hexa64(:) ) = pts(:)
      case (HEXA_98)
        return_value( hexa98(:) ) = pts(:)
      case (HEXA_125)
        return_value( hexa125(:) ) = pts(:)
      case default
        return_value = pts(:)
    end select
    !
  else
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
  end if
  !
end function reorder_cgns_hexa_pts
!
!###############################################################################
!
elemental function get_face_elem(elem,side,host_cell) result(return_value)
  !
  !.. Formal Arguments ..
  type(elem_t), intent(in) :: elem
  integer,      intent(in) :: side
  integer,      intent(in) :: host_cell
  !
  !.. Function Result ..
  type(elem_t) :: return_value
  !
  !.. Local Scalars ..
  integer :: nn,ierr
  !
  !.. Local Allocatable Scalars ..
  type(elem_t), allocatable :: tmp_elem
  !
continue
  !
  nn = size(elem%nodes)
  !
  ! Make sure the element points are in CGNS ordering
  ! before trying to extract the face points
  !
  if (all(elem%nodes == elem%pts(1:nn))) then
    tmp_elem = elem
  else
    tmp_elem = reorder_cgns_pts(elem,reverse=true)
  end if
  !
  ! Create a new face element from the requested side of the host element
  !
  select case (tmp_elem%prop%geom)
    case (Geom_Edge)
      return_value = get_edge_face_elem(tmp_elem,side)
    case (Geom_Tria)
      return_value = get_tria_face_elem(tmp_elem,side)
    case (Geom_Quad)
      return_value = get_quad_face_elem(tmp_elem,side)
    case (Geom_Tetr)
      return_value = get_tetr_face_elem(tmp_elem,side)
    case (Geom_Pyra)
      return_value = get_pyra_face_elem(tmp_elem,side)
    case (Geom_Pris)
      return_value = get_pris_face_elem(tmp_elem,side)
    case (Geom_Hexa)
      return_value = get_hexa_face_elem(tmp_elem,side)
  end select
  !
  ! Convert the new face points from CGNS ordering
  ! to the ordering needed by the solver
  !
  return_value = reorder_cgns_pts(return_value)
  !
  ! Save the host cell index to the new face element
  !
  return_value%host_cell = host_cell
  !
  ! Explicitly deallocate tmp_elem to make sure
  ! all of its components are deallocated
  !
  if (allocated(tmp_elem)) deallocate ( tmp_elem , stat=ierr )
  !
end function get_face_elem
!
!###############################################################################
!
pure function get_edge_face_elem(elem,side) result(return_value)
  !
  !.. Formal Arguments ..
  type(elem_t), intent(in) :: elem
  integer,      intent(in) :: side
  !
  !.. Function Result ..
  type(elem_t) :: return_value
  !
continue
  !
  return_value%prop = element_properties(NODE)
  return_value%tags = elem%tags ! USING F2003 AUTO-REALLOCATION
  !
  if (side == 1) then
    return_value%pts   = elem%pts(1:1)    ! USING F2003 AUTO-REALLOCATION
  else
    return_value%pts   = elem%pts(2:2)    ! USING F2003 AUTO-REALLOCATION
  end if
  !
  return_value%nodes = return_value%pts ! USING F2003 AUTO-REALLOCATION
  !
end function get_edge_face_elem
!
!###############################################################################
!
pure function get_tria_face_elem(elem,side) result(return_value)
  !
  !.. Formal Arguments ..
  type(elem_t), intent(in) :: elem
  integer,      intent(in) :: side
  !
  !.. Function Result ..
  type(elem_t) :: return_value
  !
  !.. Local Parameters ..
  integer, parameter :: tria_p1(2,3) = reshape( [1,2, &
                                                 2,3, &
                                                 3,1], &
                                                shape(tria_p1) )
  !
  integer, parameter :: tria_p2(3,3) = reshape( [1,2,4, &
                                                 2,3,5, &
                                                 3,1,6], &
                                                shape(tria_p2) )
  !
  integer, parameter :: tria_p3(4,3) = reshape( [1,2,4,5, &
                                                 2,3,6,7, &
                                                 3,1,8,9], &
                                                shape(tria_p3) )
  !
  integer, parameter :: tria_p4(5,3) = reshape( [1,2, 4, 5, 6, &
                                                 2,3, 7, 8, 9, &
                                                 3,1,10,11,12], &
                                                shape(tria_p4) )
  !
continue
  !
  return_value%tags = elem%tags ! USING F2003 AUTO-REALLOCATION
  !
  select case (elem%prop%elem_type)
    case (TRI_6)
      !
      return_value%prop = element_properties(BAR_3)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( tria_p2(:,side) )
      !
    case (TRI_9,TRI_10)
      !
      return_value%prop = element_properties(BAR_4)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( tria_p3(:,side) )
      !
    case (TRI_12,TRI_15)
      !
      return_value%prop = element_properties(BAR_5)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( tria_p4(:,side) )
      !
    case default
      !
      return_value%prop = element_properties(BAR_2)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( tria_p1(:,side) )
      !
  end select
  !
  ! USING F2003 AUTO-REALLOCATION
  return_value%nodes = return_value%pts(1:num_face_nodes(side,Geom_Tria))
  !
end function get_tria_face_elem
!
!###############################################################################
!
pure function get_quad_face_elem(elem,side) result(return_value)
  !
  !.. Formal Arguments ..
  type(elem_t), intent(in) :: elem
  integer,      intent(in) :: side
  !
  !.. Function Result ..
  type(elem_t) :: return_value
  !
  !.. Local Parameters ..
  integer, parameter :: quad_p1(2,4) = reshape( [1,2, &
                                                 2,3, &
                                                 3,4, &
                                                 4,1], &
                                                shape(quad_p1) )
  !
  integer, parameter :: quad_p2(3,4) = reshape( [1,2,5, &
                                                 2,3,6, &
                                                 3,4,7, &
                                                 4,1,8], &
                                                shape(quad_p2) )
  !
  integer, parameter :: quad_p3(4,4) = reshape( [1,2, 5, 6, &
                                                 2,3, 7, 8, &
                                                 3,4, 9,10, &
                                                 4,1,11,12], &
                                                shape(quad_p3) )
  !
  integer, parameter :: quad_p4(5,4) = reshape( [1,2, 5, 6, 7, &
                                                 2,3, 8, 9,10, &
                                                 3,4,11,12,13, &
                                                 4,1,14,15,16], &
                                                shape(quad_p4) )
  !
continue
  !
  return_value%tags = elem%tags ! USING F2003 AUTO-REALLOCATION
  !
  select case (elem%prop%elem_type)
    case (QUAD_8,QUAD_9)
      !
      return_value%prop = element_properties(BAR_3)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( quad_p2(:,side) )
      !
    case (QUAD_12,QUAD_16)
      !
      return_value%prop = element_properties(BAR_4)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( quad_p3(:,side) )
      !
    case (QUAD_P4_16,QUAD_25)
      !
      return_value%prop = element_properties(BAR_5)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( quad_p4(:,side) )
      !
    case default
      !
      return_value%prop = element_properties(BAR_2)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( quad_p1(:,side) )
      !
  end select
  !
  ! USING F2003 AUTO-REALLOCATION
  return_value%nodes = return_value%pts(1:num_face_nodes(side,Geom_Quad))
  !
end function get_quad_face_elem
!
!###############################################################################
!
pure function get_tetr_face_elem(elem,side) result(return_value)
  !
  !.. Formal Arguments ..
  type(elem_t), intent(in) :: elem
  integer,      intent(in) :: side
  !
  !.. Function Result ..
  type(elem_t) :: return_value
  !
  !.. Local Parameters ..
  integer, parameter :: tetr_p1(3,4) = reshape( [2,3,4, &
                                                 3,1,4, &
                                                 1,2,4, &
                                                 1,3,2], &
                                                shape(tetr_p1) )
  !
  integer, parameter :: tetr_p2(6,4) = reshape( [2,3,4,6,10, 9, &
                                                 3,1,4,7, 8,10, &
                                                 1,2,4,5, 9, 8, &
                                                 1,3,2,7, 6, 5], &
                                                shape(tetr_p2) )
  !
  integer, parameter :: tetr_e3(9,4) = reshape( [2,3,4, 7, 8,15,16,14,13, &
                                                 3,1,4, 9,10,11,12,16,15, &
                                                 1,2,4, 5, 6,13,14,12,11, &
                                                 1,3,2,10, 9, 8, 7, 6, 5], &
                                                shape(tetr_e3) )
  !
  integer, parameter :: tetr_p3(10,4) = reshape( [2,3,4, 7, 8,15,16,14,13,19, &
                                                  3,1,4, 9,10,11,12,16,15,20, &
                                                  1,2,4, 5, 6,13,14,12,11,18, &
                                                  1,3,2,10, 9, 8, 7, 6, 5,17], &
                                                shape(tetr_p3) )
  !
  integer, parameter :: tetr_e4(12,4) = reshape( &
    [2,3,4, 8, 9,10,20,21,22,19,18,17, &
     3,1,4,11,12,13,14,15,16,22,21,20, &
     1,2,4, 5, 6, 7,17,18,19,16,15,14, &
     1,3,2,13,12,11,10, 9, 8, 7, 6, 5], &
                                                shape(tetr_e4) )
  !
  integer, parameter :: tetr_p4(15,4) = reshape( &
    [2,3,4, 8, 9,10,20,21,22,19,18,17,29,30,31, &
     3,1,4,11,12,13,14,15,16,22,21,20,32,33,34, &
     1,2,4, 5, 6, 7,17,18,19,16,15,14,26,27,28, &
     1,3,2,13,12,11,10, 9, 8, 7, 6, 5,23,25,24], &
                                                shape(tetr_p4) )

  !
continue
  !
  return_value%tags = elem%tags ! USING F2003 AUTO-REALLOCATION
  !
  select case (elem%prop%elem_type)
    !
    case (TETRA_10)
      !
      return_value%prop = element_properties(TRI_6)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( tetr_p2(:,side) )
      !
    case (TETRA_16)
      !
      return_value%prop = element_properties(TRI_9)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( tetr_e3(:,side) )
      !
    case (TETRA_20)
      !
      return_value%prop = element_properties(TRI_10)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( tetr_p3(:,side) )
      !
    case (TETRA_22)
      !
      return_value%prop = element_properties(TRI_12)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( tetr_e4(:,side) )
      !
    case (TETRA_34,TETRA_35)
      !
      return_value%prop = element_properties(TRI_15)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( tetr_p4(:,side) )
      !
    case default
      !
      return_value%prop = element_properties(TRI_3)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( tetr_p1(:,side) )
      !
  end select
  !
  ! USING F2003 AUTO-REALLOCATION
  return_value%nodes = return_value%pts(1:num_face_nodes(side,Geom_Tetr))
  !
end function get_tetr_face_elem
!
!###############################################################################
!
pure function get_pyra_face_elem(elem,side) result(return_value)
  !
  !.. Formal Arguments ..
  type(elem_t), intent(in) :: elem
  integer,      intent(in) :: side
  !
  !.. Function Result ..
  type(elem_t) :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Parameters ..
  integer, parameter :: p1_pts(5) = [ 3, 3, 3, 3, 4]
  integer, parameter :: e2_pts(5) = [ 6, 6, 6, 6, 8]
  integer, parameter :: p2_pts(5) = [ 6, 6, 6, 6, 9]
  integer, parameter :: e3_pts(5) = [ 9, 9, 9, 9,12]
  integer, parameter :: p3_pts(5) = [10,10,10,10,16]
  integer, parameter :: e4_pts(5) = [12,12,12,12,16]
  integer, parameter :: p4_pts(5) = [15,15,15,15,25]
  !
  integer, parameter :: p1_face_type(5) = [(TRI_3 ,i=1,4),QUAD_4]
  integer, parameter :: e2_face_type(5) = [(TRI_6 ,i=1,4),QUAD_8]
  integer, parameter :: p2_face_type(5) = [(TRI_6 ,i=1,4),QUAD_9]
  integer, parameter :: e3_face_type(5) = [(TRI_9 ,i=1,4),QUAD_12]
  integer, parameter :: p3_face_type(5) = [(TRI_10,i=1,4),QUAD_16]
  integer, parameter :: e4_face_type(5) = [(TRI_12,i=1,4),QUAD_P4_16]
  integer, parameter :: p4_face_type(5) = [(TRI_15,i=1,4),QUAD_25]
  !
  integer, parameter :: pyra_p1(4,5) = reshape( [2,3,5,0, &
                                                 3,4,5,0, &
                                                 4,1,5,0, &
                                                 1,2,5,0, &
                                                 1,4,3,2], &
                                                shape(pyra_p1) )
  !
  integer, parameter :: pyra_e2(8,5) = reshape( [2,3,5,7,12,11,0,0, &
                                                 3,4,5,8,13,12,0,0, &
                                                 4,1,5,9,10,13,0,0, &
                                                 1,2,5,6,11,10,0,0, &
                                                 1,4,3,2, 9, 8,7,6], &
                                                shape(pyra_e2) )
  !
  integer, parameter :: pyra_p2(9,5) = reshape( [2,3,5,7,12,11,0,0, 0, &
                                                 3,4,5,8,13,12,0,0, 0, &
                                                 4,1,5,9,10,13,0,0, 0, &
                                                 1,2,5,6,11,10,0,0, 0, &
                                                 1,4,3,2, 9, 8,7,6,14], &
                                                shape(pyra_p2) )
  !
  integer, parameter :: pyra_e3(12,5) = reshape( &
    [2,3,5, 8, 9,18,19,17,16,0,0,0, &
     3,4,5,10,11,20,21,19,18,0,0,0, &
     4,1,5,12,13,14,15,21,20,0,0,0, &
     1,2,5, 6, 7,16,17,15,14,0,0,0, &
     1,4,3, 2,13,12,11,10, 9,8,7,6], &
                                                shape(pyra_e3) )
  !
  integer, parameter :: pyra_p3(16,5) = reshape( &
    [2,3,5, 8, 9,18,19,17,16,27,0,0, 0, 0, 0, 0, &
     3,4,5,10,11,20,21,19,18,28,0,0, 0, 0, 0, 0, &
     4,1,5,12,13,14,15,21,20,29,0,0, 0, 0, 0, 0, &
     1,2,5, 6, 7,16,17,15,14,26,0,0, 0, 0, 0, 0, &
     1,4,3, 2,13,12,11,10, 9, 8,7,6,22,25,24,23], &
                                                shape(pyra_p3) )
  !
  integer, parameter :: pyra_e4(16,5) = reshape( &
    [2,3,5, 9,10,11,24,25,26,23,22,21,0,0,0,0, &
     3,4,5,12,13,14,27,28,29,26,25,24,0,0,0,0, &
     4,1,5,15,16,17,18,19,20,29,28,27,0,0,0,0, &
     1,2,5, 6, 7, 8,21,22,23,20,19,18,0,0,0,0, &
     1,4,3, 2,17,16,15,14,13,12,11,10,9,8,7,6], &
                                                shape(pyra_e4) )
  !
  integer, parameter :: pyra_p4(25,5) = reshape( &
    [2,3,5, 9,10,11,24,25,26,23,22,21,42,43,44,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
     3,4,5,12,13,14,27,28,29,26,25,24,45,46,47,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
     4,1,5,15,16,17,18,19,20,29,28,27,48,49,50,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
     1,2,5, 6, 7, 8,21,22,23,20,19,18,39,40,41,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
     1,4,3, 2,17,16,15,14,13,12,11,10, 9, 8, 7,6,30,37,36,35,34,33,32,31,38], &
                                                shape(pyra_p4) )

  !
continue
  !
  return_value%tags = elem%tags ! USING F2003 AUTO-REALLOCATION
  !
  select case (elem%prop%elem_type)
    !
    case (PYRA_13)
      !
      return_value%prop = element_properties(e2_face_type(side))
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( pyra_e2(1:e2_pts(side),side) )
      !
    case (PYRA_14)
      !
      return_value%prop = element_properties(p2_face_type(side))
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( pyra_p2(1:p2_pts(side),side) )
      !
    case (PYRA_21)
      !
      return_value%prop = element_properties(e3_face_type(side))
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( pyra_e3(1:e3_pts(side),side) )
      !
    case (PYRA_29,PYRA_30)
      !
      return_value%prop = element_properties(p3_face_type(side))
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( pyra_p3(1:p3_pts(side),side) )
      !
    case (PYRA_P4_29)
      !
      return_value%prop = element_properties(e4_face_type(side))
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( pyra_e4(1:e4_pts(side),side) )
      !
    case (PYRA_50,PYRA_55)
      !
      return_value%prop = element_properties(p4_face_type(side))
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( pyra_p4(1:p4_pts(side),side) )
      !
    case default
      !
      return_value%prop = element_properties(p1_face_type(side))
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( pyra_p1(1:p1_pts(side),side) )
      !
  end select
  !
  ! USING F2003 AUTO-REALLOCATION
  return_value%nodes = return_value%pts(1:num_face_nodes(side,Geom_Pyra))
  !
end function get_pyra_face_elem
!
!###############################################################################
!
pure function get_pris_face_elem(elem,side) result(return_value)
  !
  !.. Formal Arguments ..
  type(elem_t), intent(in) :: elem
  integer,      intent(in) :: side
  !
  !.. Function Result ..
  type(elem_t) :: return_value
  !
  !.. Local Scalars ..
  integer :: i
  !
  !.. Local Parameters ..
  integer, parameter :: p1_pts(5) = [ 3, 3, 4, 4, 4]
  integer, parameter :: e2_pts(5) = [ 6, 6, 8, 8, 8]
  integer, parameter :: p2_pts(5) = [ 6, 6, 9, 9, 9]
  integer, parameter :: e3_pts(5) = [ 9, 9,12,12,12]
  integer, parameter :: p3_pts(5) = [10,10,16,16,16]
  integer, parameter :: e4_pts(5) = [12,12,16,16,16]
  integer, parameter :: p4_pts(5) = [15,15,25,25,25]
  !
  integer, parameter :: p1_face_type(5) = [(TRI_3 ,i=1,2),(QUAD_4    ,i=1,3)]
  integer, parameter :: e2_face_type(5) = [(TRI_6 ,i=1,2),(QUAD_8    ,i=1,3)]
  integer, parameter :: p2_face_type(5) = [(TRI_6 ,i=1,2),(QUAD_9    ,i=1,3)]
  integer, parameter :: e3_face_type(5) = [(TRI_9 ,i=1,2),(QUAD_12   ,i=1,3)]
  integer, parameter :: p3_face_type(5) = [(TRI_10,i=1,2),(QUAD_16   ,i=1,3)]
  integer, parameter :: e4_face_type(5) = [(TRI_12,i=1,2),(QUAD_P4_16,i=1,3)]
  integer, parameter :: p4_face_type(5) = [(TRI_15,i=1,2),(QUAD_25   ,i=1,3)]
  !
  integer, parameter :: pris_p1(4,5) = reshape( [1,3,2,0, &
                                                 4,5,6,0, &
                                                 1,2,5,4, &
                                                 2,3,6,5, &
                                                 3,1,4,6], &
                                                shape(pris_p1) )
  !
  integer, parameter :: pris_e2(8,5) = reshape( [1,3,2, 9, 8, 7, 0, 0, &
                                                 4,5,6,13,14,15, 0, 0, &
                                                 1,2,5, 4, 7,11,13,10, &
                                                 2,3,6, 5, 8,12,14,11, &
                                                 3,1,4, 6, 9,10,15,12], &
                                                shape(pris_e2) )
  !
  integer, parameter :: pris_p2(9,5) = reshape( [1,3,2, 9, 8, 7, 0, 0, 0, &
                                                 4,5,6,13,14,15, 0, 0, 0, &
                                                 1,2,5, 4, 7,11,13,10,16, &
                                                 2,3,6, 5, 8,12,14,11,17, &
                                                 3,1,4, 6, 9,10,15,12,18], &
                                                shape(pris_p2) )
  !
  integer, parameter :: pris_e3(12,5) = reshape( &
    [1,3,2,12,11,10, 9, 8, 7, 0, 0, 0, &
     4,5,6,19,20,21,22,23,24, 0, 0, 0, &
     1,2,5, 4, 7, 8,15,16,20,19,14,13, &
     2,3,6, 5, 9,10,17,18,22,21,16,15, &
     3,1,4, 6,11,12,13,14,24,23,18,17], &
                                                shape(pris_e3) )
  !
  integer, parameter :: pris_p3(16,5) = reshape( &
    [1,3,2,12,11,10, 9, 8, 7,25, 0, 0, 0, 0, 0, 0, &
     4,5,6,19,20,21,22,23,24,38, 0, 0, 0, 0, 0, 0, &
     1,2,5, 4, 7, 8,15,16,20,19,14,13,26,27,28,29, &
     2,3,6, 5, 9,10,17,18,22,21,16,15,30,31,32,33, &
     3,1,4, 6,11,12,13,14,24,23,18,17,34,35,36,37], &
                                                shape(pris_p3) )
  !
  integer, parameter :: pris_e4(16,5) = reshape( &
    [1,3,2,15,14,13,12,11,10, 9, 8, 7, 0, 0, 0, 0, &
     4,5,6,25,26,27,28,29,30,31,32,33, 0, 0, 0, 0, &
     1,2,5, 4, 7, 8, 9,19,20,21,27,26,25,18,17,16, &
     2,3,6, 5,10,11,12,22,23,24,30,29,28,21,20,19, &
     3,1,4, 6,13,14,15,16,17,18,33,32,31,24,23,22], &
                                                shape(pris_e4) )
  !
  integer, parameter :: pris_p4(25,5) = reshape( &
    [1,3,2,15,14,13,12,11,10, 9, 8, 7,34,36,35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
     4,5,6,25,26,27,28,29,30,31,32,33,64,65,66, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
     1,2,5, 4, 7, 8, 9,19,20,21,27,26,25,18,17,16,37,38,39,40,41,42,43,44,45, &
     2,3,6, 5,10,11,12,22,23,24,30,29,28,21,20,19,46,47,48,49,50,51,52,53,54, &
     3,1,4, 6,13,14,15,16,17,18,33,32,31,24,23,22,55,56,57,58,59,60,61,62,63], &
                                                shape(pris_p4) )
  !
continue
  !
  return_value%tags = elem%tags ! USING F2003 AUTO-REALLOCATION
  !
  select case (elem%prop%elem_type)
    !
    case (PENTA_15)
      !
      return_value%prop = element_properties(e2_face_type(side))
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( pris_e2(1:e2_pts(side),side) )
      !
    case (PENTA_18)
      !
      return_value%prop = element_properties(p2_face_type(side))
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( pris_p2(1:p2_pts(side),side) )
      !
    case (PENTA_24)
      !
      return_value%prop = element_properties(e3_face_type(side))
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( pris_e3(1:e3_pts(side),side) )
      !
    case (PENTA_38,PENTA_40)
      !
      return_value%prop = element_properties(p3_face_type(side))
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( pris_p3(1:p3_pts(side),side) )
      !
    case (PENTA_33)
      !
      return_value%prop = element_properties(e4_face_type(side))
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( pris_e4(1:e4_pts(side),side) )
      !
    case (PENTA_66,PENTA_75)
      !
      return_value%prop = element_properties(p4_face_type(side))
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( pris_p4(1:p4_pts(side),side) )
      !
    case default
      !
      return_value%prop = element_properties(p1_face_type(side))
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( pris_p1(1:p1_pts(side),side) )
      !
  end select
  !
  ! USING F2003 AUTO-REALLOCATION
  return_value%nodes = return_value%pts(1:num_face_nodes(side,Geom_Pris))
  !
end function get_pris_face_elem
!
!###############################################################################
!
pure function get_hexa_face_elem(elem,side) result(return_value)
  !
  !.. Formal Arguments ..
  type(elem_t), intent(in) :: elem
  integer,      intent(in) :: side
  !
  !.. Function Result ..
  type(elem_t) :: return_value
  !
  !.. Local Parameters ..
  integer, parameter :: hexa_p1(4,6) = reshape( [1,4,3,2, &
                                                 5,6,7,8, &
                                                 2,3,7,6, &
                                                 1,5,8,4, &
                                                 1,2,6,5, &
                                                 3,4,8,7], &
                                                shape(hexa_p1) )
  !
  integer, parameter :: hexa_e2(8,6) = reshape( [1,4,3,2,12,11,10, 9, &
                                                 5,6,7,8,17,18,19,20, &
                                                 2,3,7,6,10,15,18,14, &
                                                 1,5,8,4,13,20,16,12, &
                                                 1,2,6,5, 9,14,17,13, &
                                                 3,4,8,7,11,16,19,15], &
                                                shape(hexa_e2) )
  !
  integer, parameter :: hexa_p2(9,6) = reshape( [1,4,3,2,12,11,10, 9,21, &
                                                 5,6,7,8,17,18,19,20,26, &
                                                 2,3,7,6,10,15,18,14,23, &
                                                 1,5,8,4,13,20,16,12,25, &
                                                 1,2,6,5, 9,14,17,13,22, &
                                                 3,4,8,7,11,16,19,15,24], &
                                                shape(hexa_p2) )
  !
  integer, parameter :: hexa_e3(12,6) = reshape( &
    [1,4,3,2,16,15,14,13,12,11,10, 9, &
     5,6,7,8,25,26,27,28,29,30,31,32, &
     2,3,7,6,11,12,21,22,28,27,20,19, &
     1,5,8,4,17,18,32,31,24,23,15,16, &
     1,2,6,5, 9,10,19,20,26,25,18,17, &
     3,4,8,7,13,14,23,24,30,29,22,21], &
                                                shape(hexa_e3) )
  !
  integer, parameter :: hexa_p3(16,6) = reshape( &
    [1,4,3,2,16,15,14,13,12,11,10, 9,33,36,35,34, &
     5,6,7,8,25,26,27,28,29,30,31,32,53,54,55,56, &
     2,3,7,6,11,12,21,22,28,27,20,19,41,42,43,44, &
     1,5,8,4,17,18,32,31,24,23,15,16,50,51,52,49, &
     1,2,6,5, 9,10,19,20,26,25,18,17,37,38,39,40, &
     3,4,8,7,13,14,23,24,30,29,22,21,45,46,47,48], &
                                                shape(hexa_p3) )
  !
  integer, parameter :: hexa_e4(16,6) = reshape( &
    [1,4,3,2,20,19,18,17,16,15,14,13,12,11,10, 9, &
     5,6,7,8,33,34,35,36,37,38,39,40,41,42,43,44, &
     2,3,7,6,12,13,14,27,28,29,38,37,36,26,25,24, &
     1,5,8,4,21,22,23,44,43,42,32,31,30,18,19,20, &
     1,2,6,5, 9,10,11,24,25,26,35,34,33,23,22,21, &
     3,4,8,7,15,16,17,30,31,32,41,40,39,29,28,27], &
                                                shape(hexa_e4) )
  !
  integer, parameter :: hexa_p4(25,6) = reshape( &
    [1,4,3,2,20,19,18,17,16,15,14,13,12,11,10, 9,45,52,51,50,49,48,47,46,53, &
     5,6,7,8,33,34,35,36,37,38,39,40,41,42,43,44,90,91,92,93,94,95,96,97,98, &
     2,3,7,6,12,13,14,27,28,29,38,37,36,26,25,24,63,64,65,66,67,68,69,70,71, &
     1,5,8,4,21,22,23,44,43,42,32,31,30,18,19,20,83,84,85,86,87,88,81,82,89, &
     1,2,6,5, 9,10,11,24,25,26,35,34,33,23,22,21,54,55,56,57,58,59,60,61,62, &
     3,4,8,7,15,16,17,30,31,32,41,40,39,29,28,27,72,73,74,75,76,77,78,79,80], &
                                                shape(hexa_p4) )

  !
continue
  !
  return_value%tags = elem%tags ! USING F2003 AUTO-REALLOCATION
  !
  select case (elem%prop%elem_type)
    !
    case (HEXA_20)
      !
      return_value%prop = element_properties(QUAD_8)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( hexa_e2(:,side) )
      !
    case (HEXA_27)
      !
      return_value%prop = element_properties(QUAD_9)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( hexa_p2(:,side) )
      !
    case (HEXA_32)
      !
      return_value%prop = element_properties(QUAD_12)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( hexa_e3(:,side) )
      !
    case (HEXA_56,HEXA_64)
      !
      return_value%prop = element_properties(QUAD_16)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( hexa_p3(:,side) )
      !
    case (HEXA_44)
      !
      return_value%prop = element_properties(QUAD_P4_16)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( hexa_e4(:,side) )
      !
    case (HEXA_98,HEXA_125)
      !
      return_value%prop = element_properties(QUAD_25)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( hexa_p4(:,side) )
      !
    case default
      !
      return_value%prop = element_properties(QUAD_4)
      ! USING F2003 AUTO-REALLOCATION
      return_value%pts = elem%pts( hexa_p1(:,side) )
      !
  end select
  !
  ! USING F2003 AUTO-REALLOCATION
  return_value%nodes = return_value%pts(1:num_face_nodes(side,Geom_Hexa))
  !
end function get_hexa_face_elem
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
#endif
end subroutine cgns_memory_usage
!
!###############################################################################
!
subroutine debug_timer(time_flag,message,processor,part_num,acc)
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: time_flag
  !
  character(len=*), intent(in) :: message
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: processor
  integer, optional, intent(in) :: part_num
  integer, optional, intent(in) :: acc
  !
#ifdef PROFILE_ON
  !.. Local Scalars ..
  real(wp), save :: beg_wall_time
  real(wp)       :: end_wall_time
  real(wp)       :: wall_time
  !
  integer, save :: procedure_counter = 0
  integer       :: i,hours,minutes,seconds,total_sec
  integer       :: output_cpu
  !
  character(len=300) :: out_message
  !
  !.. Local Arrays ..
  real(wp), save :: beg_procedure_time(10)
  real(wp), save :: accumulated_time(10)
  real(wp), save :: accumulation_beg_time(10)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: leading_space = "                    "
  !
continue
  !
  output_cpu = 0
  if (present(processor)) output_cpu = processor
  !
  if (time_flag == start_timer) then
    !
    if (present(part_num)) then
      write (out_message,6) message,part_num
    else
      write (out_message,5) message
    end if
    !
    i = 3 + procedure_counter*2
    !
    beg_wall_time = mpi_wtime()
   !beg_wall_time = dclock()
   !call cpu_time( beg_wall_time )
    !
    if (mypnum == output_cpu) then
      write (iout,1) leading_space(1:i),trim(adjustl(out_message))
    end if
    !
  else if (time_flag == stop_timer) then
    !
    end_wall_time = mpi_wtime()
   !end_wall_time = dclock()
   !call cpu_time( end_wall_time )
    wall_time = end_wall_time - beg_wall_time
    !
    total_sec = nint(wall_time)
    !
    hours   = total_sec / 3600
    minutes = (total_sec - hours*3600) / 60
    seconds = total_sec - hours*3600 - minutes*60
    !
    if (present(part_num)) then
      write (out_message,6) message,part_num
    else
      write (out_message,5) message
    end if
    !
    i = 3 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,2) leading_space(1:i),wall_time,trim(adjustl(out_message))
    end if
    !
    ! Reset beg_wall_time in case this procedure is called twice
    ! in a row with time_flag set to stop_timer
    !
    beg_wall_time = mpi_wtime()
   !beg_wall_time = dclock()
   !call cpu_time( beg_wall_time )
    !
  else if (time_flag == reset_accumulated_time) then
    !
    if (present(acc)) then
      accumulated_time(acc) = zero
    end if
    !
  else if (time_flag == report_accumulated_time) then
    !
    if (present(part_num)) then
      write (out_message,6) message,part_num
    else
      write (out_message,5) message
    end if
    !
    i = 3 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,2) leading_space(1:i),accumulated_time(acc), &
                     trim(adjustl(out_message))
    end if
    !
  else if (time_flag == start_accumulation) then
    !
    if (present(acc)) then
      accumulation_beg_time(acc) = mpi_wtime()
    end if
    !
  else if (time_flag == stop_accumulation) then
    !
    if (present(acc)) then
      wall_time = mpi_wtime() - accumulation_beg_time(acc)
      accumulated_time(acc) = accumulated_time(acc) + wall_time
    end if
    !
  else if (time_flag == entering_procedure) then
    !
    i = 1 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,3) leading_space(1:i),trim(adjustl(message))
    end if
    !
    procedure_counter = procedure_counter + 1
    !
    beg_procedure_time(procedure_counter) = mpi_wtime()
   !beg_procedure_time(procedure_counter) = dclock()
   !call cpu_time( beg_procedure_time(procedure_counter) )
    !
  else if (time_flag == leaving_procedure) then
    !
    end_wall_time = mpi_wtime()
   !end_wall_time = dclock()
   !call cpu_time( end_wall_time )
    wall_time = end_wall_time - beg_procedure_time(procedure_counter)
    !
    procedure_counter = max( 0 , procedure_counter-1 )
    !
    i = 1 + procedure_counter*2
    !
    if (mypnum == output_cpu) then
      write (iout,4) leading_space(1:i),trim(adjustl(message)),wall_time
    end if
    !
  end if
  !
  flush (iout)
  !
  ! Format Statements
  !
  1 format (a,"## Starting timer for: ",a)
  2 format (a,"     -- TIME WAS",f12.6," seconds for ",a)
  3 format (a,"=> Entering procedure : ",a)
  4 format (a," <= Leaving procedure : ",a,"   after",f12.6," seconds")
  5 format (a)
  6 format (a," for partition # ",i0)
  !
#endif
  !
end subroutine debug_timer
!
!##############################################################################
!
include "Functions/bc_string_to_integer.f90"
include "Functions/last.f90"
include "Functions/smc000_interior_wall_radius.f90"
!
!###############################################################################
!
end module module_cgns
!
!###############################################################################
!
