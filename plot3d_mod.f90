module module_plot3d
  !
#ifdef DEBUG_ON
!!!#define DEBUG_CUTS
#endif
  !
  !.. Use Statements ..
  use module_kind_types
  use ovar, only : output_plot3d_to_tecplot
  !
  implicit none
  !
  private
  !
  ! ##################################
  ! ###  Module Public Procedures  ###
  ! ##################################
  !
  ! Public Interface for plot3d_element_properties function
  !
  interface plot3d_element_properties
    module procedure element_properties
  end interface plot3d_element_properties
  !
  public :: read_plot3d_gridfile
  public :: plot3d_memory_usage
  public :: plot3d_element_properties
  !
  ! ###########################
  ! ###  Module Parameters  ###
  ! ###########################
  !
  character(len=*), parameter :: bcs_default = "bcs.p3d"
  character(len=*), parameter :: cut_default = "cut.p3d"
  !
  character(len=*), parameter :: &
     errstr = " ######################  ERROR!  ######################"
  !
  ! ################################
  ! ###  Module Local Variables  ###
  ! ################################
  !
  logical(lk), save :: is_multiblock
  logical(lk), save :: is_formatted
  logical(lk), save :: is_dbl_precision
  logical(lk), save :: blk_header_is_dbl_precision
  logical(lk), save :: is_big_endian
  !
  integer, save :: nblk
  integer, save :: nbc
  integer, save :: ncut
  !
  logical(lk), save :: plot3d_has_periodic_faces = fals
  !
  character(len=150), save :: grid_folder
  !
  ! ###################################
  ! ###  Module Allocatable Arrays  ###
  ! ###################################
  !
  ! Arrays whose size is dependent on the number of blocks, nblk
  integer, save, allocatable, dimension(:,:) :: offset
  integer, save, allocatable, dimension(:,:) :: ijkpts
  !
  ! Arrays whose size is dependent on the number of boundary conditions
  !
  type condition
    integer :: block
    integer :: dir ! 1=i,2=j,3=k
    integer :: dir1
    integer :: dir2
    integer :: loc ! 1=imin,2=imax,3=jmin,4=jmax,5=kmin,6=kmax
    integer, dimension(1:3) :: begpt
    integer, dimension(1:3) :: endpt
  end type condition
  !
  type boundary_condition
    integer :: bc_type
    integer :: nfaces
    type(condition) :: bc
  end type boundary_condition
  !
  type cut_pair
    character(len=20), dimension(1:2) :: cut_name
    type(condition), dimension(1:2) :: side
    integer :: points
    integer, allocatable, dimension(:,:) :: nodes
  end type cut_pair
  !
  type(boundary_condition), save, allocatable, dimension(:) :: bcs_info
  type(cut_pair), save, allocatable, dimension(:) :: cut_info
  !
  interface map1dto3d
    module procedure map1dto3d_known_block,map1dto3d_unknown_block
  end interface map1dto3d
  !
contains
!
!###############################################################################
!
subroutine read_plot3d_gridfile ( gridfile , cutfile , bcsfile )
  !
  !.. Use Statements ..
  use geovar, only : grid
  !
  ! This routine reads in a multi-block grid in plot3d
  ! format and converts it to an unstructured grid.
  !
  ! GRIDFILE    : Gives the name of the Plot3D grid file to use
  ! BCSFILE     : Gives the name of the file that specifies the boundary
  !               conditions for block boundaries
  ! CUTFILE     : Gives the name of the file that specifies the cut conditions
  !               for a multiblock Plot3D file
  ! NOTE: Single block format means the first line of the grid file does not
  !       give the number of blocks and immediately gives the dimensions of
  !       the one and only block.
  !       Multi block format means the first line of the grid file does give
  !       the number of blocks. The second line (and however many more lines
  !       it takes) gives the dimensions of each block.
  !       HOWEVER, a single block grid can be written in a multi block format,
  !       this just means the first line would contain only the value "1".
  !
  ! NOTE: In order for this Plot3D grid to be used by the FR unstructured
  !       solver, this subroutine and the helper subroutines that it calls need
  !       to define the following arrays/variables in the module geovar_mod:
  !             NR : Number of grid dimensions
  !             NNODE : Total number of grid nodes
  !             NCELL : Total number of grid cells
  !             NFBND : Total number of boundary edges/faces
  !             XYZ_NODES : Cartesian coordinates of all the grid nodes
  !             BFACE : Boundary information for each boundary edge/face
  !             NODES_OF_CELL_PTR : Pointer array for array nodes of cells
  !             NODES_OF_CELL : Gives the nodes that define each grid cell
  !
  !.. Formal Arguments ..
  character(len=*),    intent(in) :: gridfile
  character(len=*), intent(inout) :: cutfile
  character(len=*), intent(inout) :: bcsfile
  !
  !.. Local Scalars ..
  integer :: ierr,itest
  integer :: iotst
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_plot3d_gridfile"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Have only the host roots read in the grid file
  !
  if (i_am_host_root) then
    !
    ! First, allocate the grid derived type so that it can be filled
    ! up as needed while processing the Plot3D grid file
    !
    if (allocated(grid)) then
      deallocate ( grid , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"grid",2,__LINE__,__FILE__,ierr,error_message)
    end if
    !
    allocate ( grid , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Default the Plot3D grid file to unformatted
    !
    is_formatted = fals
    !
    ! Reading a multi-block grid in Plot3d format
    !
    ! Check to see if the Plot3D grid file can be
    ! opened and read as a formatted file.
    !
    open (newunit=iotst,file=gridfile,form="formatted",status="old", &
                        action="read",position="rewind",iostat=ierr, &
                        iomsg=error_message)
    call io_error(pname,trim(adjustl(gridfile)),1,__LINE__,__FILE__, &
                  ierr,error_message)
    !
    ! Try to read the grid file as a formatted file.
    ! Mark the grid file as formated if this succeeds.
    !
    read (iotst,*,iostat=ierr) itest
    if (ierr == 0) then
      is_formatted = true
    end if
    !
    close (iotst,iostat=ierr,iomsg=error_message)
    call io_error(pname,trim(adjustl(gridfile)),2,__LINE__,__FILE__, &
                  ierr,error_message)
    !
    ! If the Plot3D grid file could not be read as a formatted file,
    ! check to see if it can be opened and read as an unformatted file.
    !
    if (.not.is_formatted) then
      !
      open (newunit=iotst,file=gridfile,form="unformatted",status="old", &
                          action="read",position="rewind",iostat=ierr, &
                          iomsg=error_message)
      call io_error(pname,trim(adjustl(gridfile)),1,__LINE__,__FILE__, &
                    ierr,error_message)
      !
      ! Try to read the grid file as an unformatted file. If this fails,
      ! something is wrong so throw an error and stop execution.
      !
      read (iotst,iostat=ierr,iomsg=error_message) itest
      if (ierr == 0) then
        is_formatted = fals
      else
        call io_error(pname,trim(adjustl(gridfile)),4,__LINE__,__FILE__, &
                      ierr,error_message)
      end if
      !
      close (iotst,iostat=ierr,iomsg=error_message)
      call io_error(pname,trim(adjustl(gridfile)),2,__LINE__,__FILE__, &
                    ierr,error_message)
      !
    end if
    !
    ! Now read the Plot3D grid file in what ever format was determined
    !
    call readplot3d(gridfile)
    !
    ! Reading the boundary condition information
    !
    call readbcsinfo(bcsfile)
    !
    ! Create grid cells from the plot3d nodes
    !
    call create_plot3d_grid_cells
    !
    ! If this is a multiblock grid, we need to do some additional
    ! work to deal with cut boundaries between blocks
    !
    if (nblk > 1) then
      !
      ! Read cut information if this is a multiblock grid
      !
      call readcutinfo(cutfile)
      !
      ! Find the points in the cuts
      !
      call findcutpoints
      !
    else
      !
      ncut = 0
      !
    end if
    !
    ! Go through the cut boundaries and remove all duplicate nodes.
    ! This shouldnt do anything if there are no cut boundaries.
    !
    call remove_duplicate_nodes
    !
    ! Reorder the grid points so that the corners of all the grid cells
    ! are first followed by any high-order cell points existing due
    ! to cell agglomeration.
    !
    call reorder_plot3d_nodes
    !
    ! Get the boundary face information
    !
    call getbface
    !
    ! Match the periodic boundaries if they exist
    !
    if (plot3d_has_periodic_faces) then
      call match_plot3d_periodic_boundaries
    end if
    !
  end if
  !
  ! Write the output file for tecplot
  !
  if (output_plot3d_to_tecplot) then
    if (mypnum == glb_root) call writetecplot
  end if
  !
 !if (mypnum == glb_root) then
 !  call create_grid_dt
 !end if
 !!
 !call root_broadcasts_grid_dt
  !
  if (ncpu > 1) then
    call broadcast_solver_geom_arrays
  else
    call create_serial_geom_arrays
  end if
  !
  ! Deallocate all module arrays before we leave
  ! since we wont be needing them anymore.
  !
  if (allocated(bcs_info)) then
    deallocate ( bcs_info , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bcs_info",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(cut_info)) then
    deallocate ( cut_info , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cut_info",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(ijkpts)) then
    deallocate ( ijkpts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ijkpts",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(offset)) then
    deallocate ( offset , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"offset",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine read_plot3d_gridfile
!
!###############################################################################
!
subroutine readplot3d(gridfile)
  !
  ! This sub reads a multi-block grid in Plot3d format
  !
  !     ijkpts(1,iblk) : number of points in the x-direction for iblk
  !     ijkpts(2,iblk) : number of points in the y-direction for iblk
  !     ijkpts(3,iblk) : number of points in the z-direction for iblk
  !
  !.. Use Statements ..
  use geovar, only : nr,xyz_nodes
  use geovar, only : grid
  use geovar, only : ncell,nnode
  use ovar,   only : plot3d_agglomerate_order
  use ovar,   only : plot3d_all_strides
  use ovar,   only : plot3d_i_stride
  use ovar,   only : plot3d_j_stride
  use ovar,   only : plot3d_k_stride
  use ovar,   only : grid_scaling_factor
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: gridfile
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,n,iblk,ierr,cpos,recl_max
  integer :: noff,npts,nnode_new,ncell_new
  integer :: ibeg,iend,jbeg,jend,kbeg,kend
  integer :: iogrd,tempval
  character(len=100) :: read_message
  !
  !.. Local Arrays ..
  integer, dimension(1:3) :: strides
  !
  !.. Local Allocatable Arrays ..
  real(r4), allocatable, dimension(:,:) :: r4_nodes
  real(wp), allocatable, dimension(:,:) :: temp_xyz
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "readplot3d"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Set all the coarsening strides to -1 to prevent coarsening
  ! if agglomerating cells to create high-order curved cells.
  !
  if (plot3d_agglomerate_order > 0) then
    plot3d_all_strides = -1
    plot3d_i_stride = -1
    plot3d_j_stride = -1
    plot3d_k_stride = -1
  end if
  !
  strides(:) = max(1,plot3d_all_strides)
  !
  if (plot3d_i_stride > 0) strides(1) = plot3d_i_stride
  if (plot3d_j_stride > 0) strides(2) = plot3d_j_stride
  if (plot3d_k_stride > 0) strides(3) = plot3d_k_stride
  !
  ! Before we start, extract the name of the
  ! folder that the grid file resides within
  !
  cpos = scan( gridfile(1:len_trim(gridfile)), '/' , back=true )
  grid_folder = gridfile(1:cpos)
  grid_folder = trim(adjustl(grid_folder))
  !
  ! Open the file for multi-block grid
  !
  if (is_formatted) then
    open (newunit=iogrd,file=gridfile,form="formatted",status="old", &
                        action="read",position="rewind",iostat=ierr, &
                        iomsg=error_message)
  else
    open (newunit=iogrd,file=gridfile,form="unformatted",status="old", &
                        action="read",position="rewind",iostat=ierr, &
                        iomsg=error_message)
  end if
  call io_error(pname,trim(adjustl(gridfile)),1,__LINE__,__FILE__, &
                ierr,error_message)
  !
  ! If this is a multiblock file, read in the number of blocks
  !
  if (is_formatted) then
    call read_p3d_formatted_header(iogrd,gridfile)
  else
    call read_p3d_unformatted_header(iogrd,gridfile)
  end if
  !
  ! Check that the requested plot3d node stride is acceptable
  !
  l = 0
  do iblk = 1,nblk
    do i = 1,nr
      if (mod(ijkpts(i,iblk)-1,strides(i)) /= 0) then
        if (mypnum == glb_root) then
          if (l == 0) write (iout,4)
          write (iout,5) iblk,i
        end if
        l = 1
      end if
    end do
  end do
  !
  ! If an error was found with the strides, only have the global root
  ! call stop_gfr with an abort flag. We cant use the stop_mpi flag
  ! here because only the host roots will call this routine, and the
  ! stop_mpi flag requires all processes to call stop_gfr together,
  ! which would result in the program hanging within stop_gfr. We
  ! only want the global root to call it with abort, to make sure that
  ! the error statements in the above loop are written to iout. If some
  ! other host root calls stop_gfr with the abort flag, there isnt a
  ! guarantee that those error statements will be written before the
  ! call to mpi_abort.
  !
  if (l == 1) then
    if (mypnum == glb_root) then
      call stop_gfr(abort,pname,__LINE__,__FILE__)
    end if
    ! Have all the host roots call mpi_barrier here so that the MPI processes
    ! that arent the global root do not continue past this point before the
    ! global root has a chance to call mpi_abort within the above conditional.
    call mpi_barrier(host_roots_comm,mpierr)
  end if
  !
  ! Allocate the node and cell offset array
  !
  allocate ( offset(1:2,1:nblk) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"offset",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Create the offset arrays for each block
  !
  offset(1,1) = 1 ! node offset
  offset(2,1) = 1 ! cell offset
  !
  do iblk = 2,nblk
    !
    offset(1,iblk) = offset(1,iblk-1) + product( ijkpts(1:nr,iblk-1) )
    offset(2,iblk) = offset(2,iblk-1) + product( ijkpts(1:nr,iblk-1)-1 )
    !
  end do
  !
  ! Get the total number of nodes and cells
  !
  nnode = sum( product( ijkpts(1:nr,1:nblk) , dim=1 ) )
  ncell = sum( product( ijkpts(1:nr,1:nblk)-1 , dim=1 ) )
  !
  ! Allocate the node coordinate array
  !
  allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"xyz_nodes",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read in the node coordinates
  !
  if (is_formatted) then
    !
    ! If the grid is formatted, this is straight
    ! forward so read in the node coordinates
    !
    do iblk = 1,nblk
      write (read_message,10) iblk
      call debug_timer(start_timer,read_message)
      noff = offset(1,iblk)-1
      npts = product(ijkpts(1:nr,iblk))
      read (iogrd,*) ((xyz_nodes(l,i+noff), i=1,npts), l=1,nr)
      call debug_timer(stop_timer,read_message)
    end do
    !
  else
    !
    ! If the grid is unformatted, we dont know if the grid is in
    ! single or double precision so we have to be careful about
    ! reading in the node coordinates.
    !
    ! In an unformatted Plot3D grid file, a separate record exists for
    ! the node coordinates of each block. Assume the grid file is in
    ! double precision and try reading in the node coordinates for the
    ! first block.
    !
    write (read_message,11) 1
    call debug_timer(start_timer,read_message)
    noff = offset(1,1)-1
    npts = product(ijkpts(1:nr,1))
    read (iogrd,iostat=ierr,iomsg=error_message) &
                        ((xyz_nodes(l,i+noff), i=1,npts), l=1,nr)
    call debug_timer(stop_timer,read_message)
    !
    if (ierr /= 0) then
      !
      write (iout,20) ierr,adjustl(trim(error_message))
      !
      ! If there was an error reading the coordinates in double precision,
      ! backspace and try to read the grid file in single precision,
      ! aborting if there is an error for single precision.
      !
      ! Allocate the single precision coordinate array
      !
      allocate ( r4_nodes(1:nr,1:nnode) , source=0.0_r4 , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"r4_nodes",1,__LINE__,__FILE__,ierr,error_message)
      !
     !backspace (iogrd,iostat=ierr,iomsg=error_message)
     !if (ierr /= 0) then
     !  write (iout,22) ierr,adjustl(trim(error_message))
     !end if
      rewind (iogrd,iostat=ierr,iomsg=error_message)
      if (ierr /= 0) then
        write (iout,21) ierr,adjustl(trim(error_message))
      end if
      read (iogrd) tempval
      read (iogrd) ((tempval,l=1,nr),i=1,nblk)
      !
     !do iblk = 1,nblk
     !  write (iout,'(*(2x,i0))') (ijkpts(l,iblk),l=1,nr)
     !end do
      do iblk = 1,nblk
        write (read_message,12) iblk
        call debug_timer(start_timer,read_message)
        noff = offset(1,iblk)-1
        npts = product(ijkpts(1:nr,iblk))
        read (iogrd,iostat=ierr,iomsg=error_message) &
                                      ((r4_nodes(l,i+noff), i=1,npts), l=1,nr)
       !do i = 1,npts
       !  write (110,'(*(es20.10))') (r4_nodes(l,i+noff),l=1,nr)
       !end do
        call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
        call debug_timer(stop_timer,read_message)
      end do
      !
      is_dbl_precision = fals
      !
      ! If the single precision coordinates were successfully read,
      ! covert them to working precision
      !
      xyz_nodes(1:nr,1:nnode) = real( r4_nodes(1:nr,1:nnode) , kind=wp )
      !
      ! Finish by deallocate the single precision array
      !
      deallocate ( r4_nodes , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"r4_nodes",2,__LINE__,__FILE__,ierr,error_message)
      !
    else
      !
      ! If there was no error reading the coordinates for the first
      ! block in double precision, continue reading the remaining blocks.
      !
      do iblk = 2,nblk
        write (read_message,11) iblk
        call debug_timer(start_timer,read_message)
        noff = offset(1,iblk)-1
        npts = product(ijkpts(1:nr,iblk))
        read (iogrd,iostat=ierr,iomsg=error_message) &
                                      ((xyz_nodes(l,i+noff), i=1,npts), l=1,nr)
        call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
        call debug_timer(stop_timer,read_message)
      end do
      !
      is_dbl_precision = true
      !
    end if
    !
  end if
  !
  ! Close the Plot3D grid file since we are finished reading it
  !
  close (iogrd,iostat=ierr,iomsg=error_message)
  call io_error(pname,trim(adjustl(gridfile)),2,__LINE__,__FILE__, &
                ierr,error_message)
  !
  if (mypnum == glb_root) then
    !
    write (iout,*)
    write (iout,1) " Finish reading the Plot3D file"
    write (iout,2) " nblk = ", nblk
    write (iout,2) " ncell = ", ncell
    write (iout,2) " nnode = ", nnode
    !
    write (iout,2) " Number of blocks = ",nblk
    do iblk = 1,nblk
      write (iout,3) iblk,(ijkpts(l,iblk),l=1,3)
    end do
    if (is_multiblock) then
      write (iout,*) "This is a multi block grid"
    else
      write (iout,*) "This is a single block grid"
    end if
    if (is_formatted) then
      write (iout,*) "This is a formatted grid"
    else
      if (is_dbl_precision) then
        write (iout,*) "This is an unformatted grid in double precision"
      else
        write (iout,*) "This is an unformatted grid in single precision"
      end if
     !if (is_big_endian) then
     !  write (iout,*) "The unformatted grid file is in big endian format"
     !else
     !  write (iout,*) "The unformatted grid file is in little endian format"
     !end if
    end if
    write (iout,*)
    !
  end if
  !
  ! Apply the grid scaling factor to the coordinates of the grid nodes
  !
  do n = 1,size(xyz_nodes,dim=2)
    do i = 1,size(xyz_nodes,dim=1)
      xyz_nodes(i,n) = grid_scaling_factor*xyz_nodes(i,n)
    end do
  end do
  !
  ! Coarsen the grid nodes if any of the coarsening strides are larger than 1
  !
  if (any(strides > 1)) then
    !
    nnode_new = 0
    do iblk = 1,size(ijkpts,dim=2)
      npts = 1
      do l = 1,size(ijkpts,dim=1)
        npts = npts * ( (ijkpts(l,iblk)-1)/strides(l) + 1 )
      end do
      nnode_new = nnode_new + npts
    end do
    !
    allocate ( temp_xyz(1:nr,1:nnode_new) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"temp_xyz",1,__LINE__,__FILE__,ierr,error_message)
    !
    n = 0
    !
    do iblk = 1,size(ijkpts,dim=2)
      !
      noff = offset(1,iblk)-1
      !
      ibeg = 1
      jbeg = 1
      kbeg = 1
      !
      iend = max( ibeg , ijkpts(1,iblk) )
      jend = max( jbeg , ijkpts(2,iblk) )
      kend = max( kbeg , ijkpts(3,iblk) )
      !
      do k = kbeg,kend,strides(3)
        do j = jbeg,jend,strides(2)
          do i = ibeg,iend,strides(1)
            !
            l = map3dto1d( [i,j,k] , iblk , 1 )
            n = n + 1
            !
            temp_xyz(1:nr,n) = xyz_nodes(1:nr,l+noff)
            !
          end do
        end do
      end do
      !
    end do
    !
    npts = n
    !
    do iblk = 1,size(ijkpts,dim=2)
      do l = 1,size(ijkpts,dim=1)
        ijkpts(l,iblk) = (ijkpts(l,iblk)-1) / strides(l) + 1
      end do
    end do
    !
    nnode_new = sum( product( ijkpts(1:nr,1:nblk)   , dim=1 ) )
    ncell_new = sum( product( ijkpts(1:nr,1:nblk)-1 , dim=1 ) )
    !
    if (npts /= nnode_new) then
      if (mypnum == glb_root) then
        write (iout,*) " ERROR WITH NEW NUMBER OF NODES"
        write (iout,*) " npts      = ",npts
        write (iout,*) " nnode_new = ",nnode_new
        call stop_gfr(abort,pname,__LINE__,__FILE__)
      end if
      ! Have all the host roots call mpi_barrier here so that the MPI processes
      ! that arent the global root do not continue past this point before the
      ! global root has a chance to call mpi_abort within the above conditional.
      call mpi_barrier(host_roots_comm,mpierr)
    end if
    !
    offset(1,1) = 1 ! node offset
    offset(2,1) = 1 ! cell offset
    !
    do iblk = 2,nblk
      !
      offset(1,iblk) = offset(1,iblk-1) + product( ijkpts(1:nr,iblk-1) )
      offset(2,iblk) = offset(2,iblk-1) + product( ijkpts(1:nr,iblk-1)-1 )
      !
    end do
    !
    if (mypnum == glb_root) then
      !
      write (iout,*)
      write (iout,1) " AFTER COARSENING USING STRIDES OF"
      write (iout,2) "    Stride in I-direction = ",strides(1)
      write (iout,2) "    Stride in J-direction = ",strides(2)
      write (iout,2) "    Stride in K-direction = ",strides(3)
      write (iout,*)
      write (iout,2) " nblk = ", nblk
      write (iout,2) " ncell = ", ncell
      write (iout,2) " nnode = ", nnode
      !
      write (iout,2) " Number of blocks = ",nblk
      do iblk = 1,nblk
        write (iout,3) iblk,(ijkpts(l,iblk),l=1,3)
      end do
      !
    end if
    !
    call move_alloc( from=temp_xyz , to=xyz_nodes )
    !
    nnode = nnode_new
    ncell = ncell_new
    !
  end if
  !
 !call stop_gfr(stop_mpi,pname,__LINE__,__FILE__, &
 !                "multiblock testing in readplot3d")
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (a)
  2 format (a,i0)
  3 format ("   Block #",i0," dimensions: ",i0,2(" x ",i0))
  4 format (/," ERROR TRYING TO USE COARSENING STRIDES WITH PLOT3D!!!",/, &
              " THE NUMBER OF GRID CELLS MUST BE EVENLY",/, &
              " DIVISIBLE BY THE COARSENING STRIDES!!!",/)
  5 format (3x,"Cells not evenly divisible in block ",i0," for direction ",i0)
  !
  10 format ("Reading formatted coordinates for block # ",i0)
  11 format ("Reading 64-bit unformatted coordinates for block # ",i0)
  12 format ("Reading 32-bit unformatted coordinates for block # ",i0)
  !
  20 format (/,1x,"read(xyz_nodes[64-bit]) iostat = ",i0,/, &
               1x,"read(xyz_nodes[64-bit])  iomsg = ",a,/)
  21 format (/,1x,"rewind(xyz_nodes[32-bit]) iostat = ",i0,/, &
               1x,"rewind(xyz_nodes[32-bit])  iomsg = ",a,/)
  22 format (/,1x,"backspace(xyz_nodes[32-bit]) iostat = ",i0,/, &
               1x,"backspace(xyz_nodes[32-bit])  iomsg = ",a,/)
  !
end subroutine readplot3d
!
!###############################################################################
!
subroutine read_p3d_formatted_header(iogrd,gridfile)
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: iogrd
  character(len=*), intent(in) :: gridfile
  !
  !.. Local Scalars ..
  integer :: ierr,num_val,n,l,lines
  character(len=100) :: text
  !
  !.. Local Allocatable Arrays ..
  character(len=50), allocatable, dimension(:) :: clist
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_p3d_formatted_header"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Read the first line as text to see how many values it contains
  !
  read (iogrd,1,iostat=ierr,iomsg=error_message) text
  call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Parse text to see how many values are contained on the first line
  !
  call parse_list(text,clist)
  !
  num_val = size(clist)
  !
  ! If the number of values on the first line is one, this
  ! is probably a multiblock grid. If it is greater than one,
  ! this is probably a single block grid file.
  !
  if (num_val > 1) then
    !
    ! This is probably a single block grid file. If this is a single
    ! block grid, num_val should be equal to the number of grid dimensions.
    ! Make sure the size of the third dimension is greater than one
    ! before declaring this a 3D grid.
    !
    nblk = 1
    !
    ! Allocate ijkpts, and initialize to 1 node in each direction
    !
    allocate ( ijkpts(1:3,1:nblk) , source=1 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ijkpts (SB)",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Make sure that num_val is equal to either 2 or 3
    !
    if (all(num_val /= [2,3])) goto 100 ! go to 100 error statement
    !
    ! Read the values in clist if num_val is equal to either 2 or 3
    !
    do l = 1,num_val
      read (clist(l),*,iostat=ierr,iomsg=error_message) ijkpts(l,1)
      call io_error(pname,"clist(l)",-1,__LINE__,__FILE__,ierr,error_message)
    end do
    !
    ! If we have reached this point, it probably safe
    ! to say this is a single block grid file
    !
    is_multiblock = fals
    !
  else if (num_val == 1) then
    !
    ! This is probably a multi block grid
    !
    read (clist(1),*,iostat=ierr,iomsg=error_message) nblk
    call io_error(pname,"clist(1)",-1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Make sure that nblk is a valid number of blocks
    !
    if (nblk <= 0) goto 101 ! go to 101 error statement
    !
    ! Allocate ijkpts, and initialize to 1 node in each direction
    !
    allocate ( ijkpts(1:3,1:nblk) , source=1 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ijkpts (MB)",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! We need to get the dimensions of each block but we dont know
    ! if this is a 3D or 2D grid. Therefore, we must continue reading
    ! lines until we find a line that does not contain integer values.
    ! To do this, we have to:
    !    1. Read each line as a character string.
    !    2. Parse this string to separate the individual values contained
    !       within the string into substrings.
    !    3. Test each substring to determine if it is a valid integer value.
    !    4. Once a non-integer value is found, this means we have reached the
    !       node coordinate section and everything before this was the block
    !       dimensions.
    ! NOTE: We have to do this because Fortran will not produce an error
    !       if a floating point value is read as an integer, and
    !       automatically converts the float to an integer.
    !
    l = 0
    blk_dims: do
      read (iogrd,1,iostat=ierr,iomsg=error_message) text
      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
      call parse_list(text,clist)
      num_val = size(clist)
      do n = 1,num_val
        if (string_is_integer(clist(n))) then
          l = l + 1
        else
          exit blk_dims
        end if
      end do
    end do blk_dims
    !
    if (mod(l,nblk) /= 0) then
      goto 102 ! go to 102 error statement
    else
      nr = l/nblk
      if (all(nr /= [2,3])) goto 100 ! go to 100 error statement
    end if
    !
    ! We have to rewind the grid file so that we
    ! can actually read in the block dimensions
    !
    rewind (iogrd)
    !
    ! Skip the first line containing the number of blocks
    !
    read (iogrd,1,iostat=ierr,iomsg=error_message) text
    call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Now read in the block dimensions
    !
    read (iogrd,*,iostat=ierr,iomsg=error_message) &
                                  ((ijkpts(l,n), l=1,nr), n=1,nblk)
    call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
    !
    ! If we have reached this point, it probably safe
    ! to say this is a single block grid file
    !
    is_multiblock = true
    !
  else
    !
    ! There was an unknown error reading the header of the Plot3D grid file
    !
    goto 103 ! go to 103 error statement
    !
  end if
  !
  ! Reading the header seemed to be successful. Make sure
  ! that all the block dimensions are valid values.
  !
 !if (mypnum == glb_root) then
 !  do n = 1,nblk
 !    write (iout,7) n,(ijkpts(l,n),l=1,3)
 !  end do
 !end if
  if (any(ijkpts(:,:) < 1)) goto 104 ! go to 104 error statement
  !
  ! If all the third dimensions are 1, this is a 2D grid
  !
  nr = merge( 2 , 3 , all(ijkpts(3,1:nblk) == 1) )
  !
  ! Deallocate clist before we leave if it is allocated
  !
  if (allocated(clist)) then
    deallocate ( clist , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"clist",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Return statement to ignore the remaining error crap
  !
  return
  !
  ! Error continuations
  !
  100 if (mypnum == glb_root) write (iout,2) errstr,errstr ; goto 999
  101 if (mypnum == glb_root) write (iout,3) errstr,errstr ; goto 999
  102 if (mypnum == glb_root) write (iout,4) errstr,errstr ; goto 999
  103 if (mypnum == glb_root) write (iout,5) errstr,errstr ; goto 999
  104 if (mypnum == glb_root) write (iout,6) errstr,errstr ; goto 999
  !
  999 continue
  if (mypnum == glb_root) then
    call stop_gfr(abort,pname,__LINE__,__FILE__)
  end if
  ! Have all the host roots call mpi_barrier here so that the MPI processes
  ! that arent the global root do not continue past this point before the
  ! global root has a chance to call mpi_abort within the above conditional.
  call mpi_barrier(host_roots_comm,mpierr)
  !
  ! Format Statements
  !
  1 format (a)
  2 format (/,a,/," ERROR READING THE FORMATTED PLOT3D GRID FILE HEADER!",/, &
                  " THE GRID FILE SEEMED TO BE IN SINGLE BLOCK FORMAT",/, &
                  " BUT THE NUMBER OF VALUES ON THE FIRST LINE OF THE",/, &
                  " GRID FILE WAS INVALID!",/,a,/)
  3 format (/,a,/," ERROR READING THE FORMATTED PLOT3D GRID FILE HEADER!",/, &
                  " THE GRID FILE SEEMED TO BE IN MULTIBLOCK FORMAT",/, &
                  " BUT THE NUMBER OF BLOCKS ON THE FIRST LINE OF THE",/, &
                  " GRID FILE WAS INVALID!",/,a,/)
  4 format (/,a,/, &
            " ERROR READING THE FORMATTED PLOT3D GRID FILE HEADER!",/, &
            " THE GRID FILE SEEMED TO BE IN MULTIBLOCK FORMAT BUT THE",/, &
            " NUMBER OF GRID DIMENSIONS WAS NOT EVENLY DIVISIBLE BY",/, &
            " THE NUMBER OF BLOCKS, TELLING US IF THE GRID IS 2D OR 3D!",/,a,/)
  5 format (/,a,/," ERROR READING THE FIRST LINE OF", &
                  " THE FORMATTED PLOT3D GRID FILE!",/,a,/)
  6 format (/,a,/," ERROR READING THE FORMATTED PLOT3D GRID FILE HEADER!",/, &
                  " THE MINIMUM NUMBER OF NODES IN ANY DIRECTION IS 1!",/,a,/)
  7 format ("Block #",i0," dimensions: ",i0,2(" x ",i0))
  !
end subroutine read_p3d_formatted_header
!
!###############################################################################
!
subroutine read_p3d_unformatted_header(iogrd,gridfile)
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  integer,          intent(in) :: iogrd
  character(len=*), intent(in) :: gridfile
  !
  !.. Local Scalars ..
  integer :: int_rec1,int_rec2
  integer :: ierr,n,l
  real(r4) :: rl4_rec2
  real(r8) :: rl8_rec2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_p3d_unformatted_header"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Whatever the Plot3D format, the first record
  ! should contain at least one integer
  !
  read (iogrd,iostat=ierr,iomsg=error_message) int_rec1
  call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
  !
  ! #######
  ! ## 1 ##
  ! #######
  !
  ! Determine whether this grid is in single or multi block format
  !
  ! Try reading the second record as an integer
  !
  read (iogrd,iostat=ierr) int_rec2
  !
  int_check: if (ierr == 0) then
    !
    ! This seems to be a multiblock grid file
    !
    is_multiblock = true
    !
  else
    !
    ! Failed to read the second record as an integer.
    ! This seems to indicate this is a single block grid file.
    !
    is_multiblock = fals
    !
    ! Backspace and try to read it again, this time as a single precision real.
    !
    backspace (iogrd)
    !
  end if int_check
  !
  ! #######
  ! ## 2 ##
  ! #######
  !
  ! Determine whether this is a 2D or 3D grid by reading the block dimensions
  !
  ! Rewind the grid file to start over
  !
  rewind (iogrd)
  !
  ! If multiblock read the number of blocks
  !
  if (is_multiblock) then
    !
    read (iogrd,iostat=ierr,iomsg=error_message) nblk
    call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Make sure that nblk is a valid number of blocks
    !
    if (nblk <= 0) goto 101 ! go to 101 error statement
    !
  else
    !
    nblk = 1
    !
  end if
  !
  ! Allocate ijkpts, and initialize to 1 node in each direction
  !
  allocate ( ijkpts(1:3,1:nblk) , source=1 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ijkpts",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Read the grid dimensions, assuming it is a 3D grid
  !
  read (iogrd,iostat=ierr) ((ijkpts(l,n), l=1,3), n=1,nblk)
  !
  ! If the read was not successful, try reading the
  ! header again but assuming it is a 2D grid.
  !
  if (ierr /= 0) then
    !
    ! Reinitialize ijkpts to 1 node in each direction
    !
    ijkpts(1:3,1:nblk) = 1
    !
    ! Try reading the block dimensions again but in 2D
    !
    backspace (iogrd)
    read (iogrd,iostat=ierr) ((ijkpts(l,n), l=1,2), n=1,nblk)
    !
    if (ierr /= 0) goto 102 ! go to 102 error statement
    !
  end if
  !
  ! Reading the header seemed to be successful. Make sure
  ! that all the block dimensions are valid values.
  !
 !if (mypnum == glb_root) then
 !  do n = 1,nblk
 !    write (iout,7) n,(ijkpts(l,n),l=1,3)
 !  end do
 !end if
  !
  if (any(ijkpts(:,:) < 1)) goto 103 ! go to 104 error statement
  !
  ! If all the third dimensions are 1, this is a 2D grid
  !
  nr = merge( 2 , 3 , all(ijkpts(3,1:nblk) == 1) )
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Return statement to ignore the remaining error crap
  !
  return
  !
  ! Error continuations
  !
  100 if (mypnum == glb_root) write (iout,2) errstr,errstr ; goto 999
  101 if (mypnum == glb_root) write (iout,3) errstr,errstr ; goto 999
  102 if (mypnum == glb_root) write (iout,4) errstr,errstr ; goto 999
  103 if (mypnum == glb_root) write (iout,5) errstr,errstr ; goto 999
  !
  999 continue
  if (mypnum == glb_root) then
    call stop_gfr(abort,pname,__LINE__,__FILE__)
  end if
  ! Have all the host roots call mpi_barrier here so that the MPI processes
  ! that arent the global root do not continue past this point before the
  ! global root has a chance to call mpi_abort within the above conditional.
  call mpi_barrier(host_roots_comm,mpierr)
  !
  ! Format Statements
  !
  2 format (/,a,/, &
            " AN ERROR OCCURED WITH THE UNFORMATTED PLOT3D GRID FILE",/, &
            " WHILE TRYING TO DETERMINE IF THE NODE COORDINATES ARE",/, &
            " STORED IN SINGLE OR DOUBLE PRECISION!",/,a,/)
  3 format (/,a,/, &
            " ERROR READING THE UNFORMATTED PLOT3D GRID FILE HEADER!",/, &
            " THE GRID FILE SEEMED TO BE IN MULTIBLOCK FORMAT",/, &
            " BUT THE NUMBER OF BLOCKS ON THE FIRST LINE OF THE",/, &
            " GRID FILE WAS INVALID!",/,a,/)
  4 format (/,a,/, &
            " ERROR READING THE UNFORMATTED PLOT3D GRID FILE HEADER!",/, &
            " THE GRID DIMENSIONS COULD NOT BE READ IN EITHER 2D OR 3D!",/,a,/)
  5 format (/,a,/, &
            " ERROR READING THE UNFORMATTED PLOT3D GRID FILE HEADER!",/, &
            " THE MINIMUM NUMBER OF NODES IN ANY DIRECTION IS 1!",/,a,/)
  7 format ("Block #",i0," dimensions: ",i0,2(" x ",i0))
  !
end subroutine read_p3d_unformatted_header
!
!###############################################################################
!
!subroutine read_p3d_unformatted_header_generic_endian(iogrd,gridfile)
!  !
!  !.. Use Statements ..
!  use geovar, only : nr
!  !
!  !.. Formal Arguments ..
!  integer,          intent(in) :: iogrd
!  character(len=*), intent(in) :: gridfile
!  !
!  !.. Local Scalars ..
!  integer :: int_rec1,int_rec2
!  integer :: ierr,n,l,npass
!  real(r4) :: rl4_rec2
!  real(r8) :: rl8_rec2
!  real(qp) :: max_int
!  real(qp) :: tot_pts_i4,tot_pts_i8
!  !
!  integer, dimension(1:20) :: ijktest
!  real(wp), dimension(1:100) :: xyztest
!  !
!  !.. Local Allocatable Arrays ..
!  integer(i8), allocatable :: ijkpts_i8(:,:)
!  !
!  !.. Local Parameters ..
!  character(len=*), parameter :: pname = "read_p3d_unformatted_header"
!  !
!continue
!  !
!  ! Determine the structure of the unformatted Plot3D file
!  !
!  call determine_unformatted_structure(gridfile)
!  !
!  ! Open the grid file and read the header
!  !
!  open (newunit=iogrd,file=gridfile,form="unformatted",status="old", &
!                      action="read",access="stream", &
!                      iostat=ierr,iomsg=error_message)
!  call io_error(pname,gridfile,1,__LINE__,__FILE__,ierr,error_message)
!  !
!  ! Initialize the file position
!  !
!  file_pos = 1
!  !
!  ! Read the first record if this is a multiblock grid file
!  !
!  if (is_multiblock) then
!    !
!    if (blk_header_is_dbl_precision) then
!      !
!      ! Read the opening record marker
!      !
!      read (iogrd,pos=1,iostat=ierr,iomsg=error_message) int_i8
!      error_message = "int_i8 "//trim(adjustl(error_message))
!      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!      !
!      ! Read the number of blocks
!      !
!      read (iogrd,iostat=ierr,iomsg=error_message) int_i8
!      error_message = "int_i8 "//trim(adjustl(error_message))
!      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!      !
!      if (is_big_endian) int_i8 = big_to_little_endian(int_i8)
!      !
!      nblk = int( int_i8 , kind=kind(nblk) )
!      !
!      ! Read the closing record marker
!      !
!      read (iogrd,iostat=ierr,iomsg=error_message) int_i8
!      error_message = "int_i8 "//trim(adjustl(error_message))
!      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!      !
!    else
!      !
!      ! Read the opening record marker
!      !
!      read (iogrd,pos=1,iostat=ierr,iomsg=error_message) int_i4
!      error_message = "int_i4 "//trim(adjustl(error_message))
!      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!      !
!      ! Read the number of blocks
!      !
!      read (iogrd,iostat=ierr,iomsg=error_message) int_i4
!      error_message = "int_i4 "//trim(adjustl(error_message))
!      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!      !
!      if (is_big_endian) int_i4 = big_to_little_endian(int_i4)
!      !
!      nblk = int( int_i4 , kind=kind(nblk) )
!      !
!      ! Read the closing record marker
!      !
!      read (iogrd,iostat=ierr,iomsg=error_message) int_i4
!      error_message = "int_i4 "//trim(adjustl(error_message))
!      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!      !
!    end if
!    !
!    ! Allocate ijkpts and initialize to 1 node in each direction
!    !
!    allocate ( ijkpts(1:3,1:nblk) , source=1 , &
!               stat=ierr , errmsg=error_message )
!    call alloc_error(pname,"ijkpts",1,__LINE__,__FILE__,ierr,error_message)
!    !
!    ! Get the current file position before continuing
!    !
!    inquire (iogrd,pos=file_pos)
!    !
!    ! Try reading the opening record marker for the block
!    ! dimensions using a 4-byte integer
!    !
!    read (iogrd,pos=file_pos,iostat=ierr,iomsg=error_message) int_i4
!    error_message = "int_i4 "//trim(adjustl(error_message))
!    call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!    !
!    if (is_big_endian) int_i4 = big_to_little_endian(int_i4)
!    !
!    nr = int( int_i4 / 4_i4 / int(nblk,kind=i4) , kind=kind(nr) )
!    !
!    if (any(nr == [2,3])) then
!      !
!      ! Read in the block dimensions
!      !
!
!      !
!    else
!      !
!      ! Try reading as an 8-byte integer
!      !
!      read (iogrd,pos=file_pos,iostat=ierr,iomsg=error_message) int_i8
!      error_message = "int_i8 "//trim(adjustl(error_message))
!      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!      !
!      if (is_big_endian) int_i8 = big_to_little_endian(int_i8)
!      !
!      nr = int( int_i8 / 8_i8 / int(nblk,kind=i8) , kind=kind(nr) )
!      !
!      if (all(nr /= [2,3])) then
!        write (error_message,20)
!        call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
!      end if
!      !
!    end if
!    !
!
!  else
!    !
!    ! This is a single block grid file so set nblk to 1
!    !
!    nblk = 1
!    !
!    ! Allocate ijkpts and initialize to 1 node in each direction
!    !
!    allocate ( ijkpts(1:3,1:nblk) , source=1 , &
!               stat=ierr , errmsg=error_message )
!    call alloc_error(pname,"ijkpts",1,__LINE__,__FILE__,ierr,error_message)
!    !
!    if (blk_header_is_dbl_precision) then
!      !
!      ! Allocate ijkpts_i8 and initialize to 1 node in each direction
!      !
!      allocate ( ijkpts_i8(1:3,1:nblk) , source=1_i8 , &
!                 stat=ierr , errmsg=error_message )
!      call alloc_error(pname,"ijkpts_i8",1,__LINE__,__FILE__,ierr, &
!                       error_message)
!      !
!      ! Read the opening record marker for the block dimensions
!      !
!      read (iogrd,pos=1,iostat=ierr,iomsg=error_message) int_i8
!      error_message = "int_i8 "//trim(adjustl(error_message))
!      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!      !
!      ! Read the block dimensions
!      !
!      read (iogrd,iostat=ierr,iomsg=error_message) (ijkpts_i8(l,1),l=1,nr)
!      error_message = "ijkpts_i8 "//trim(adjustl(error_message))
!      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!      !
!      if (is_big_endian) then
!        ijkpts_i8(1:nr,:) = big_to_little_endian(ijkpts_i8(1:nr,:))
!      end if
!      !
!      ! Read the closing record marker for the block dimensions
!      !
!      read (iogrd,iostat=ierr,iomsg=error_message) int_i8
!      error_message = "int_i8 "//trim(adjustl(error_message))
!      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!      !
!      ! Copy ijkpts_i8 into ijkpts
!      !
!      ijkpts = int( ijkpts_i8 , kind=kind(ijkpts) )
!      !
!      ! Deallocate ijkpts_i8
!      !
!      deallocate ( ijkpts_i8 , stat=ierr , errmsg=error_message )
!      call alloc_error(pname,"ijkpts_i8",2,__LINE__,__FILE__,ierr, &
!                       error_message)
!      !
!    else
!      !
!      ! Allocate ijkpts_i4 and initialize to 1 node in each direction
!      !
!      allocate ( ijkpts_i4(1:3,1:nblk) , source=1_i4 , &
!                 stat=ierr , errmsg=error_message )
!      call alloc_error(pname,"ijkpts_i4",1,__LINE__,__FILE__,ierr, &
!                       error_message)
!      !
!      ! Read the opening record marker for the block dimensions
!      !
!      read (iogrd,pos=1,iostat=ierr,iomsg=error_message) int_i4
!      error_message = "int_i4 "//trim(adjustl(error_message))
!      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!      !
!      ! Read the block dimensions
!      !
!      read (iogrd,iostat=ierr,iomsg=error_message) (ijkpts_i4(l,1),l=1,nr)
!      error_message = "ijkpts_i4 "//trim(adjustl(error_message))
!      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!      !
!      if (is_big_endian) then
!        ijkpts_i4(1:nr,:) = big_to_little_endian(ijkpts_i4(1:nr,:))
!      end if
!      !
!      ! Read the closing record marker for the block dimensions
!      !
!      read (iogrd,iostat=ierr,iomsg=error_message) int_i4
!      error_message = "int_i4 "//trim(adjustl(error_message))
!      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!      !
!      ! Copy ijkpts_i4 into ijkpts
!      !
!      ijkpts = int( ijkpts_i4 , kind=kind(ijkpts) )
!      !
!      ! Deallocate ijkpts_i4
!      !
!      deallocate ( ijkpts_i4 , stat=ierr , errmsg=error_message )
!      call alloc_error(pname,"ijkpts_i4",2,__LINE__,__FILE__,ierr, &
!                       error_message)
!      !
!    end if
!    !
!  end if
!  ! Try reading the marker using 4-byte integer
!  !
!  read (iogrd,pos=file_pos,iostat=ierr
!  ! Test to see if the
!
!  !
!  !
!  ! Loop through this procedure twice incase
!  ! there the grid file is in big endian
!  !
!  endian_check_loop: do npass = 1,2
!    !
!    rewind (iogrd)
!    !
!    ! Whatever the Plot3D format, the first record
!    ! should contain at least one integer
!    !
!    read (iogrd,iostat=ierr) int_rec1
!    call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,"int_rec1")
!    !
!    if (npass == 2) int_rec1 = big_to_little_endian( int_rec1 )
!    !
!    ! #######
!    ! ## 1 ##
!    ! #######
!    !
!    ! Determine whether this grid is in single or multi block format
!    !
!    ! Try reading the second record as an integer
!    !
!    read (iogrd,iostat=ierr) int_rec2
!    !
!    int_check: if (ierr == 0) then
!      !
!      ! This seems to be a multiblock grid file
!      !
!      is_multiblock = true
!      !
!    else
!      !
!      ! Failed to read the second record as an integer.
!      ! This seems to indicate this is a single block grid file.
!      !
!      is_multiblock = fals
!      !
!      ! Backspace and try to read it again,
!      ! this time as a single precision real.
!      !
!      backspace (iogrd)
!      !
!    end if int_check
!    !
!    ! #######
!    ! ## 2 ##
!    ! #######
!    !
!    ! Determine whether this is a 2D or 3D grid by reading the block dimensions
!    !
!    ! Rewind the grid file to start over
!    !
!    rewind (iogrd)
!    !
!    ! If multiblock read the number of blocks
!    !
!    if (is_multiblock) then
!      !
!      read (iogrd,iostat=ierr) nblk
!      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,"nblk")
!      !
!      if (npass == 2) nblk = big_to_little_endian( nblk )
!      !
!      ! Make sure that nblk is a valid number of blocks
!      !
!      if (nblk <= 0) goto 101 ! go to 101 error statement
!      !
!    else
!      !
!      nblk = 1
!      !
!    end if
!    !
!    ! Allocate ijkpts, and initialize to 1 node in each direction
!    !
!    if (allocated(ijkpts)) then
!      deallocate ( ijkpts , stat=ierr , errmsg=error_message )
!      call alloc_error(pname,"ijkpts",2,__LINE__,__FILE__,ierr,error_message)
!    end if
!    !
!    allocate ( ijkpts(1:3,1:nblk) , source=1 , &
!               stat=ierr , errmsg=error_message )
!    call alloc_error(pname,"ijkpts",1,__LINE__,__FILE__,ierr,error_message)
!    !
!    ! Read the grid dimensions, assuming it is a 3D grid
!    !
!   !read (iogrd,iostat=ierr) ((ijkpts(l,n), l=1,3), n=1,nblk)
!    ijktest = 0
!    read (iogrd,iostat=ierr) (ijktest(l), l=1,2)
!    if (npass == 2) ijktest = big_to_little_endian(ijktest)
!    write (iout,*) ierr,(ijktest(l),l=1,2)
!    xyztest = real(0.0,kind=kind(xyztest))
!    read (iogrd,iostat=ierr) (xyztest(l),l=1,size(xyztest))
!    if (npass == 2) xyztest = big_to_little_endian(xyztest)
!    write (iout,*) ierr,(xyztest(l),l=1,size(xyztest))
!    cycle endian_check_loop
!    !
!    if (npass == 2) ijkpts = big_to_little_endian( ijkpts )
!    !
!    ! If the read was not successful, try reading the
!    ! header again but assuming it is a 2D grid.
!    !
!    if (ierr /= 0) then
!      !
!      ! Reinitialize ijkpts to 1 node in each direction
!      !
!      ijkpts(1:3,1:nblk) = 1
!      !
!      ! Try reading the block dimensions again but in 2D
!      !
!      backspace (iogrd)
!      read (iogrd,iostat=ierr) ((ijkpts(l,n), l=1,2), n=1,nblk)
!      !
!      if (npass == 2) then
!        ijkpts(1:2,1:nblk) = big_to_little_endian( ijkpts(1:2,1:nblk) )
!      end if
!      !
!      if (ierr /= 0) goto 102 ! go to 102 error statement
!      !
!    end if
!    !
!    ! Reading the header seemed to be successful. Make sure
!    ! that all the block dimensions are valid values.
!    !
!    if (any(ijkpts(:,:) < 1)) goto 103 ! go to 104 error statement
!    !
!    ! If all the third dimensions are 1, this is a 2D grid
!    !
!    nr = merge( 2 , 3 , all(ijkpts(3,1:nblk) == 1) )
!    !
!    ! Make an i8 copy of ijkpts
!    !
!    if (allocated(ijkpts_i8)) then
!      deallocate ( ijkpts_i8 , stat=ierr , errmsg=error_message )
!      call alloc_error(pname,"ijkpts_i8",2,__LINE__,__FILE__,ierr, &
!      error_message)
!    end if
!    !
!    allocate ( ijkpts_i8(1:3,1:nblk) , source=int(ijkpts,kind=i8) , &
!               stat=ierr , errmsg=error_message )
!    call alloc_error(pname,"ijkpts_i8",1,__LINE__,__FILE__,ierr,error_message)
!    !
!    ! Get the total number of points in i4 and i8 integer types
!    !
!    tot_pts_i4 = sum( product( real( abs(ijkpts) , kind=qp ) , dim=1 ) )
!    tot_pts_i8 = sum( product( real( abs(ijkpts_i8) , kind=qp ) , dim=1 ) )
!    write (iout,*) "tot_pts_i4 = ",tot_pts_i4
!    write (iout,*) "tot_pts_i8 = ",tot_pts_i8
!    write (iout,*)
!    write (iout,*) "ijkpts"
!    do n = 1,nblk
!      write (iout,'(*(2x,i0))') (ijkpts(l,n),l=1,3)
!    end do
!    write (iout,*)
!    write (iout,*) "ijkpts_i8"
!    do n = 1,nblk
!      write (iout,'(*(2x,i0))') (ijkpts_i8(l,n),l=1,3)
!    end do
!    !
!    ! If these two values are equal, more than likely the Plot3D header
!    ! was successfully read
!    !
!    if (tot_pts_i4 < real(huge(0),kind=qp)) then
!      is_big_endian = (npass == 2)
!      exit endian_check_loop
!    else if (tot_pts_i8 < real(huge(0_i8),kind=qp)) then
!      is_big_endian = (npass == 2)
!      exit endian_check_loop
!    else if (npass == 2) then
!      goto 104
!    end if
!    !
!  end do endian_check_loop
!  !
!  if (allocated(ijkpts_i8)) then
!    deallocate ( ijkpts_i8 , stat=ierr , errmsg=error_message )
!    call alloc_error(pname,"ijkpts_i8",2,__LINE__,__FILE__,ierr,error_message)
!  end if
!  !
!  ! Return statement to ignore the remaining error crap
!  !
!  return
!  !
!  ! Error continuations
!  !
!  100 if (mypnum == 0) write (iout,2) errstr,errstr ; goto 999
!  101 if (mypnum == 0) write (iout,3) errstr,errstr ; goto 999
!  102 if (mypnum == 0) write (iout,4) errstr,errstr ; goto 999
!  103 if (mypnum == 0) write (iout,5) errstr,errstr ; goto 999
!  104 if (mypnum == 0) write (iout,6) errstr,errstr ; goto 999
!  !
!  999 call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
!  !
!  ! Format Statements
!  !
!  2 format (/,a,/, &
!            " AN ERROR OCCURED WITH THE UNFORMATTED PLOT3D GRID FILE",/, &
!            " WHILE TRYING TO DETERMINE IF THE NODE COORDINATES ARE",/, &
!            " STORED IN SINGLE OR DOUBLE PRECISION!",/,a,/)
!  3 format (/,a,/, &
!            " ERROR READING THE UNFORMATTED PLOT3D GRID FILE HEADER!",/, &
!                  " THE GRID FILE SEEMED TO BE IN MULTIBLOCK FORMAT",/, &
!                  " BUT THE NUMBER OF BLOCKS ON THE FIRST LINE OF THE",/, &
!                  " GRID FILE WAS INVALID!",/,a,/)
!  4 format (/,a,/, &
!            " ERROR READING THE UNFORMATTED PLOT3D GRID FILE HEADER!",/, &
!            " THE GRID DIMENSIONS COULD NOT BE READ IN EITHER 2D OR 3D!",/,a,/)
!  5 format (/,a,/, &
!            " ERROR READING THE UNFORMATTED PLOT3D GRID FILE HEADER!",/, &
!                  " THE MINIMUM NUMBER OF NODES IN ANY DIRECTION IS 1!",/,a,/)
!  6 format (/,a,/, &
!            " ERROR READING THE UNFORMATTED PLOT3D GRID FILE HEADER!",/, &
!            " FAILED TO READ THE HEADER IN LITTLE OR BIG ENDIAN FORMAT!",/,a,/)
!  !
!end subroutine read_p3d_unformatted_header_generic_endian
!
!###############################################################################
!
!!subroutine determine_unformatted_structure(gridfile)
!!  !
!!  character(len=*), intent(in) :: gridfile
!!  !
!!  !.. Local Scalars ..
!!  integer :: iogrd,n,ierr
!!  logical(lk) :: first_record_successfully_read
!!  integer(i4) :: rec_header_i4,int_i4,ncount_i4
!!  integer(i8) :: rec_header_i8,int_i8,ncount_i8
!!  character(len=200) :: error_message
!!  !
!!  !.. Local Parameters ..
!!  character(len=*), parameter :: pname = "determine_unformatted_structure"
!!  !
!!continue
!!  !
!!  first_record_successfully_read = fals
!!  !
!!  open (newunit=iogrd,file=gridfile,form="unformatted",status="old", &
!!                      action="read",access="stream", &
!!                      iostat=ierr,iomsg=error_message)
!!  call io_error(pname,gridfile,1,__LINE__,__FILE__,ierr,error_message)
!!  !
!!  ! Attempt to read record header using 4-byte integer
!!  !
!!  read (iogrd,pos=1,iostat=ierr,iomsg=error_message) rec_header_i4
!!  error_message = "rec_header_i4 "//trim(adjustl(error_message))
!!  call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!!  !
!! !write (iout,4) "rec_header_i4 = ",rec_header_i4
!!  !
!!  ! Skip the next value in the record
!!  !
!!  read (iogrd,iostat=ierr,iomsg=error_message) int_i4
!!  error_message = "int_i4, skip, "//trim(adjustl(error_message))
!!  call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!!  !
!! !write (iout,4) "Skip,   int_i4 = ",int_i4
!!  !
!!  ! Try to find the matching record header.
!!  ! NOTE: Use a single count do-loop so we can exit
!!  !       once we find the matching record header.
!!  ! NOTE: The record headers should match regardless of the files endianess.
!!  !
!!  i4_header_test: do n = 1,1
!!    !
!!    ! For a multi-block Plot3D file the first record should look like:
!!    !    beg_record_header , number_of_blocks , end_record_header
!!    ! so the matching record header should be in the third file position.
!!    !
!!    read (iogrd,iostat=ierr,iomsg=error_message) int_i4
!!    error_message = "int_i4, MB 3D, "//trim(adjustl(error_message))
!!    call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!!    !
!!   !write (iout,4) "MB, 3D, int_i4 = ",int_i4
!!    !
!!    if (int_i4 == rec_header_i4) then
!!      is_multiblock = true
!!      blk_header_is_dbl_precision = fals
!!      first_record_successfully_read = true
!!      ncount_i4 = 1_i4
!!      exit i4_header_test
!!    end if
!!    !
!!    ! For a signle-block 2D Plot3D file the first record should look like:
!!    !    beg_record_header , nx , ny , end_record_header
!!    ! so the matching record header should be in the fourth file position.
!!    !
!!    read (iogrd,iostat=ierr,iomsg=error_message) int_i4
!!    error_message = "int_i4, SB 2D, "//trim(adjustl(error_message))
!!    call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!!    !
!!   !write (iout,4) "SB, 2D, int_i4 = ",int_i4
!!    !
!!    if (int_i4 == rec_header_i4) then
!!      is_multiblock = fals
!!      blk_header_is_dbl_precision = fals
!!      first_record_successfully_read = true
!!      ncount_i4 = 2_i4
!!      nr = 2
!!      exit i4_header_test
!!    end if
!!    !
!!    ! For a single-block 3D Plot3D file the first record should look like:
!!    !    beg_record_header , nx , ny , nz , end_record_header
!!    ! so the matching record header should be in the fifth file position.
!!    !
!!    read (iogrd,iostat=ierr,iomsg=error_message) int_i4
!!    error_message = "int_i4, SB 3D, "//trim(adjustl(error_message))
!!    call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!!    !
!!   !write (iout,4) "SB, 3D, int_i4 = ",int_i4
!!    !
!!    if (int_i4 == rec_header_i4) then
!!      is_multiblock = fals
!!      blk_header_is_dbl_precision = fals
!!      first_record_successfully_read = true
!!      ncount_i4 = 3_i4
!!      nr = 3
!!      exit i4_header_test
!!    end if
!!    !
!!  end do i4_header_test
!!  !
!!  ! If the previous loop wasnt successful, try again using 8-byte integers
!!  !
!!  if (.not. first_record_successfully_read) then
!!    !
!!    ! Attempt to read record header using 8-byte integer
!!    !
!!    read (iogrd,pos=1,iostat=ierr,iomsg=error_message) rec_header_i8
!!    error_message = "rec_header_i8 "//trim(adjustl(error_message))
!!    call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!!    !
!!   !write (iout,4) "rec_header_i8 = ",rec_header_i8
!!    !
!!    ! Skip the next value in the record
!!    !
!!    read (iogrd,iostat=ierr,iomsg=error_message) int_i8
!!    error_message = "int_i8, skip, "//trim(adjustl(error_message))
!!    call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!!    !
!!   !write (iout,4) "Skip,   int_i8 = ",int_i8
!!    !
!!    i8_header_test: do n = 1,1
!!      !
!!      ! For a multi-block Plot3D file the first record should look like:
!!      !    beg_record_header , number_of_blocks , end_record_header
!!      ! so the matching record header should be in the third file position.
!!      !
!!      read (iogrd,iostat=ierr,iomsg=error_message) int_i8
!!      error_message = "int_i8, MB 3D, "//trim(adjustl(error_message))
!!      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!!      !
!!     !write (iout,4) "MB, 3D, int_i8 = ",int_i8
!!      !
!!      if (int_i8 == rec_header_i8) then
!!        is_multiblock = true
!!        blk_header_is_dbl_precision = true
!!        first_record_successfully_read = true
!!        ncount_i8 = 1_i8
!!        exit i8_header_test
!!      end if
!!      !
!!      ! For a signle-block 2D Plot3D file the first record should look like:
!!      !    beg_record_header , nx , ny , end_record_header
!!      ! so the matching record header should be in the fourth file position.
!!      !
!!      read (iogrd,iostat=ierr,iomsg=error_message) int_i8
!!      error_message = "int_i8, SB 2D, "//trim(adjustl(error_message))
!!      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!!      !
!!     !write (iout,4) "SB, 2D, int_i8 = ",int_i8
!!      !
!!      if (int_i8 == rec_header_i8) then
!!        is_multiblock = fals
!!        blk_header_is_dbl_precision = true
!!        first_record_successfully_read = true
!!        ncount_i8 = 2_i8
!!        nr = 2
!!        exit i8_header_test
!!      end if
!!      !
!!      ! For a single-block 3D Plot3D file the first record should look like:
!!      !    beg_record_header , nx , ny , nz , end_record_header
!!      ! so the matching record header should be in the fifth file position.
!!      !
!!      read (iogrd,iostat=ierr,iomsg=error_message) int_i8
!!      error_message = "int_i8, SB 3D, "//trim(adjustl(error_message))
!!      call io_error(pname,gridfile,-1,__LINE__,__FILE__,ierr,error_message)
!!      !
!!     !write (iout,4) "SB, 3D, int_i8 = ",int_i8
!!      !
!!      if (int_i8 == rec_header_i8) then
!!        is_multiblock = fals
!!        blk_header_is_dbl_precision = true
!!        first_record_successfully_read = true
!!        ncount_i8 = 3_i8
!!        nr = 3
!!        exit i8_header_test
!!      end if
!!      !
!!    end do i8_header_test
!!    !
!!  end if
!!  !
!!  ! If neither read attempts were successful, exit with an error
!!  !
!!  if (.not. first_record_successfully_read) then
!!    write (error_message,1)
!!    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
!!  end if
!!  !
!!  ! Determine the endianess of the file
!!  !
!!  if (blk_header_is_dbl_precision) then
!!    !
!!    if (rec_header_i8/ncount_i8 == 8) then
!!      ! File is in native endianess
!!      is_big_endian = fals
!!    else if (big_to_little_endian(rec_header_i8)/ncount_i8 == 8) then
!!      ! File is in non-native endianess
!!      is_big_endian = true
!!    else
!!      write (error_message,2)
!!    end if
!!    !
!!  else
!!    !
!!    if (rec_header_i4/ncount_i4 == 4) then
!!      ! File is in native endianess
!!      is_big_endian = fals
!!    else if (big_to_little_endian(rec_header_i4)/ncount_i4 == 4) then
!!      ! File is in non-native endianess
!!      is_big_endian = true
!!    else
!!      write (error_message,3)
!!    end if
!!    !
!!  end if
!! !write (iout,*) "is_big_endian = ",is_big_endian
!! !write (iout,*) "is_multiblock = ",is_multiblock
!! !write (iout,*) "blk_header_is_dbl_precision = ",blk_header_is_dbl_precision
!!  !
!!  ! This should be everything we need to successfully read the Plot3D file
!!  ! so close the file and return
!!  !
!!  close (iogrd,iostat=ierr,iomsg=error_message)
!!  call io_error(pname,gridfile,2,__LINE__,__FILE__,ierr,error_message)
!!  !
!!  ! Format Statements
!!  !
!!  1 format ("UNABLE TO READ THE FIRST RECORD OF", &
!!            "THE UNFORMATTED PLOT3D GRID FILE!")
!!  2 format ("THE FIRST RECORD OF THE UNFORMATTED PLOT3D GRID FILE SEEMS", &
!!            " TO BE AN 8-BYTE RECORD BUT THE ENDIANESS OF THE FILE", &
!!            " COULD NOT BE DETERMINED!")
!!  3 format ("THE FIRST RECORD OF THE UNFORMATTED PLOT3D GRID FILE SEEMS", &
!!            " TO BE A 4-BYTE RECORD BUT THE ENDIANESS OF THE FILE", &
!!            " COULD NOT BE DETERMINED!")
!!  4 format (a,i0)
!!  !
!!  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__)
!!  !
!!end subroutine determine_unformatted_structure
!
!###############################################################################
!
subroutine create_plot3d_grid_cells()
  !
  !.. Use Statements ..
  use geovar, only : nr,ncell,nnode,xyz_nodes,grid
  use geovar, only : nodes_of_cell_ptr,nodes_of_cell
  use ovar,   only : plot3d_agglomerate_order
  !
  !.. Local Scalars ..
  integer :: l,n,nb,nc,nf,nn,np,location,dir,ie
  integer :: i,ic,ip,ifa,ibeg,iend,npi
  integer :: j,jc,jp,jfa,jbeg,jend,npj
  integer :: k,kc,kp,kfa,kbeg,kend,npk
  integer :: d1,dir1,beg1,end1,n1,nc1
  integer :: d2,dir2,beg2,end2,n2,nc2
  integer :: elem_type,face_type,bc_type
  integer :: cell_pts,geom_pts,iblk,ierr
  integer :: lcell,lfbnd,lnode,stride
  integer :: tcell,tnode,tc,tn
  !
  !.. Local Arrays ..
  integer, dimension(1:3) :: ijk,ijk0,ijk1,ijk2,ijkf
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: ipg(:),ipc(:)
  integer, allocatable :: cells_with_node_ptr(:)
  integer, allocatable :: cells_with_node(:)
  integer, allocatable :: lpoin(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_plot3d_grid_cells"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  stride = max(1,plot3d_agglomerate_order)
  !
  ! Count the total number of cells across all blocks in this grid
  !
  call debug_timer(start_timer,"Count cells in all blocks")
  l = 0
  lcell = 0
  lnode = 0
  tcell = 0
  tnode = 0
  do iblk = 1,size(ijkpts,dim=2)
    !
    ! Check that the agglomeration order is acceptable for this block
    !
    nc = 1
    nn = 1
    tc = 1
    tn = 1
    do n = 1,size(ijkpts,1)
      !
      if (mod(ijkpts(n,iblk)-1,stride) /= 0) then
        if (mypnum == glb_root) then
          if (l == 0) write (iout,1)
          write (iout,2) iblk,n
        end if
        l = 1
      end if
      !
      ! Multiply the number of cells for this block
      ! by the number of cells for this direction.
      !
      nc = nc * max( 1 , (ijkpts(n,iblk)-1)/stride )
      nn = nn * ( (ijkpts(n,iblk)-1)/stride + 1 )
      !
      tc = tc * max( 1 , ijkpts(n,iblk)-1 )
      tn = tn * ijkpts(n,iblk)
      !
    end do
    !
    ! Add the number of cells for this block to
    ! the total number of cells for this grid
    !
    lcell = lcell + nc
    lnode = lnode + nn
    !
    tcell = tcell + tc
    tnode = tnode + tn
    !
  end do
  !
  ncell = lcell
  call debug_timer(stop_timer,"Count cells in all blocks")
  !
 !write (iout,'(a," = ",i0)') "lcell",lcell
 !write (iout,'(a," = ",i0)') "tcell",tcell
 !write (iout,'(a," = ",i0)') "lnode",lnode
 !write (iout,'(a," = ",i0)') "tnode",tnode
  !
  ! Count the total number of boundary faces across all blocks in this grid
  !
  call debug_timer(start_timer,"Counting boundary faces across all blocks")
  lfbnd = 0
  do nb = 1,size(bcs_info)
    !
    iblk = bcs_info(nb)%bc%block
    !
    dir1 = bcs_info(nb)%bc%dir1
    beg1 = bcs_info(nb)%bc%begpt(dir1)
    end1 = bcs_info(nb)%bc%endpt(dir1)
    !
    dir2 = bcs_info(nb)%bc%dir2
    beg2 = bcs_info(nb)%bc%begpt(dir2)
    end2 = bcs_info(nb)%bc%endpt(dir2)
    !
    nc1 = end1 - beg1
    nc2 = end2 - beg2
    !
    ! Check that the agglomeration order is acceptable for this block
    !
    if (mod(nc1,stride) /= 0) then
      if (mypnum == glb_root) then
        if (l == 0) write (iout,1)
        write (iout,3) nb,iblk,dir1
      end if
      l = 1
    end if
    ! NOTE: This should be fine if nr==2 where nc2==0
    if (mod(nc2,stride) /= 0) then
      if (mypnum == glb_root) then
        if (l == 0) write (iout,1)
        write (iout,3) nb,iblk,dir2
      end if
      l = 1
    end if
    !
    if (nc2 == 0) then
      nc2 = stride ! Make sure nc2 is equal to the stride if nr==2 where nc2==0
    end if
    !
    lfbnd = lfbnd + (nc1/stride)*(nc2/stride)
    !
  end do
 !write (iout,'(a," = ",i0)') "lfbnd",lfbnd
  call debug_timer(stop_timer,"Counting boundary faces across all blocks")
  !
  ! If an error was found with the strides, abort the simulation
  !
  if (l == 1) then
    if (mypnum == glb_root) then
      call stop_gfr(abort,pname,__LINE__,__FILE__)
    end if
    ! Have all the host roots call mpi_barrier here so that the MPI processes
    ! that arent the global root do not continue past this point before the
    ! global root has a chance to call mpi_abort within the above conditional.
    call mpi_barrier(host_roots_comm,mpierr)
  end if
  !
  ! Allocate the elem component of the grid derived type
  !
  allocate ( grid%elem(1:lcell+lfbnd) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid%elem",1,__LINE__,__FILE__,ierr,error_message)
  !
  geom_pts = 2**nr
  cell_pts = (stride+1)**nr
  elem_type = merge(Geom_Hexa,Geom_Quad,nr==3)
  face_type = merge(Geom_Quad,Geom_Edge,nr==3)
  !
  allocate ( ipg(1:geom_pts) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ipg",1,__LINE__,__FILE__,ierr,error_message)
  allocate ( ipc(1:cell_pts) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ipc",1,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(start_timer,"Creating grid%elem(1:ncell)")
  nc = 0
  !
  do iblk = 1,size(ijkpts,dim=2)
    !
    npi = ijkpts(1,iblk)
    npj = ijkpts(2,iblk)
    npk = ijkpts(3,iblk)
    !
    kc_loop: do kc = 1,max(1,npk-1),stride
      jc_loop: do jc = 1,max(1,npj-1),stride
        ic_loop: do ic = 1,max(1,npi-1),stride
          !
          nc = nc + 1
          !
          ipg(1) = map3dto1d( [ic,jc,kc] , iblk , 1 )
          ipg(2) = ipg(1) + stride
          ipg(3) = ipg(2) + stride*npi
          ipg(4) = ipg(3) - stride
          !
          np = 0
          if (nr == 3) then
            ipg(5:8) = ipg(1:4) + stride*npi*npj
            kp3_loop: do kp = 0,stride
              k = kc + kp
              jp3_loop: do jp = 0,stride
                j = jc + jp
                i = ic
                np = np + 1
                ipc(np) = map3dto1d( [i,j,k] , iblk , 1 )
                ip3_loop: do ip = 1,stride
                  np = np + 1
                  ipc(np) = ipc(np-1) + 1
                end do ip3_loop
               !ip3_loop: do ip = 0,stride
               !  i = ic + ip
               !  np = np + 1
               !  ipc(np) = map3dto1d( [i,j,k] , iblk , 1 )
               !end do ip3_loop
              end do jp3_loop
            end do kp3_loop
          else if (nr == 2) then
            k = 1
            jp2_loop: do jp = 0,stride
              j = jc + jp
              ip2_loop: do ip = 0,stride
                i = ic + ip
                np = np + 1
                ipc(np) = map3dto1d( [i,j,k] , iblk , 1 )
              end do ip2_loop
            end do jp2_loop
          end if
          !
          grid%elem(nc)%prop = element_properties(elem_type,stride)
          !
          ! USING F2003 AUTO-REALLOCATION
          !
          grid%elem(nc)%nodes = ipg
          grid%elem(nc)%pts = ipc
          grid%elem(nc)%tags = [0,0]
          !
        end do ic_loop
      end do jc_loop
    end do kc_loop
    !
  end do
  call debug_timer(stop_timer,"Creating grid%elem(1:ncell)")
  !
 !do nc = 1,lcell
 !  n1 = size(grid%elem(nc)%nodes)
 !  n2 = size(grid%elem(nc)%pts)
 !  write (153,153) nc,(grid%elem(nc)%nodes(n),n=1,n1)
 !  write (154,153) nc,(grid%elem(nc)%pts(n),n=1,n2)
 !end do
 !153 format ("nc = ",i4," : ",*(1x,i4,:))
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  call debug_timer(start_timer,"Creating nodes_of_cell")
  allocate ( nodes_of_cell_ptr(1:lcell+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  do nc = 2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(nc) = nodes_of_cell_ptr(nc-1) + &
                            size(grid%elem(nc-1)%nodes)
  end do
  !
  allocate ( nodes_of_cell(1:last(nodes_of_cell_ptr)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr,error_message)
  !
  do nc = 1,lcell
    n1 = nodes_of_cell_ptr(nc)+1
    n2 = nodes_of_cell_ptr(nc+1)
    nodes_of_cell(n1:n2) = grid%elem(nc)%nodes(:)
  end do
  call debug_timer(stop_timer,"Creating nodes_of_cell")
  !
 !do nc = 1,size(nodes_of_cell_ptr)-1
 !  n1 = nodes_of_cell_ptr(nc)+1
 !  n2 = nodes_of_cell_ptr(nc+1)
 !  write (151,151) nc,n1,n2,(nodes_of_cell(n),n=n1,n2)
 !end do
 !151 format ("nc = ",i4," (",i4,"->",i4,") : ",*(1x,i4,:))
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  call debug_timer(start_timer,"Creating cells_with_node")
 !allocate ( cells_with_node_ptr(1:lnode+1) , source=0 , &
  allocate ( cells_with_node_ptr(1:tnode+1) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_with_node_ptr",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  do nc = 1,lcell
    n1 = nodes_of_cell_ptr(nc)+1
    n2 = nodes_of_cell_ptr(nc+1)
    do n = n1,n2
      np = nodes_of_cell(n) + 1
      cells_with_node_ptr(np) = cells_with_node_ptr(np) + 1
    end do
  end do
  !
  do n = 2,size(cells_with_node_ptr)
    cells_with_node_ptr(n) = cells_with_node_ptr(n) + cells_with_node_ptr(n-1)
  end do
  !
  allocate ( cells_with_node(1:last(cells_with_node_ptr)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cells_with_node",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  do nc = 1,lcell
    n1 = nodes_of_cell_ptr(nc)+1
    n2 = nodes_of_cell_ptr(nc+1)
    do n = n1,n2
      j = nodes_of_cell(n)
      np = cells_with_node_ptr(j) + 1
      cells_with_node_ptr(j) = np
      cells_with_node(np) = nc
    end do
  end do
  !
  do n = size(cells_with_node_ptr),2,-1
    cells_with_node_ptr(n) = cells_with_node_ptr(n-1)
  end do
  cells_with_node_ptr(1) = 0
  call debug_timer(stop_timer,"Creating cells_with_node")
  !
 !do nc = 1,size(cells_with_node_ptr)-1
 !  n1 = cells_with_node_ptr(nc)+1
 !  n2 = cells_with_node_ptr(nc+1)
 !  write (152,152) nc,n1,n2,(cells_with_node(n),n=n1,n2)
 !end do
 !152 format ("nn = ",i4," (",i4,"->",i4,") : ",*(1x,i4,:))
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  if (allocated(ipc)) then
    deallocate ( ipc , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ipc",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(ipg)) then
    deallocate ( ipg , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ipg",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  geom_pts = 2**(nr-1)
  cell_pts = (stride+1)**(nr-1)
  !
  allocate ( ipg(1:geom_pts) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ipg",1,__LINE__,__FILE__,ierr,error_message)
  allocate ( ipc(1:cell_pts) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ipc",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate the point counting array lpoin
  !
 !allocate ( lpoin(1:lnode) , source=0 , stat=ierr , errmsg=error_message )
  allocate ( lpoin(1:tnode) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"lpoin",1,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(start_timer,"Creating grid%elem(ncell+1:ncell+nfbnd)")
  nc = lcell
  !
  do nb = 1,size(bcs_info)
    !
    ! Boundary condition type for this face
    !
    bc_type = bcs_info(nb)%bc_type
    !
   !if (bc_type == bc_periodic) then
   !  plot3d_has_periodic_faces = true
   !end if
    !
    ! Block number for this face
    !
    iblk = bcs_info(nb)%bc%block
    !
    ! Side/Direction (i,j,k) of this BC and its location (min,max)
    !   dir: [1,2,3]       = [i,j,k]
    !   loc: [1,2,3,4,5,6] = [imin,imax,jmin,jmax,kmin,kmax]
    !
    dir = bcs_info(nb)%bc%dir
    location = bcs_info(nb)%bc%loc
    !
    ijk0(:) = 0
    ijk0(dir) = mod(location+1,2)
    !
    ! Beginning index in each direction
    !
    ibeg = bcs_info(nb)%bc%begpt(1)
    jbeg = bcs_info(nb)%bc%begpt(2)
    kbeg = bcs_info(nb)%bc%begpt(3)
    !
    ! Ending index in each direction
    ! NOTE: The index in bcs_info(nb)%bc%endpt(:) gives the node index that
    !       the boundary condition ends at. We have to subtract one from
    !       this index to convert it to a cell index. We also have to make
    !       sure that an end index does not become 0 after subtracting one.
    !       This would happen for whatever index corresponds to the face the
    !       BC is on if it is a MIN face.
    !
    iend = max( 1 , bcs_info(nb)%bc%endpt(1)-1 )
    jend = max( 1 , bcs_info(nb)%bc%endpt(2)-1 )
    kend = max( 1 , bcs_info(nb)%bc%endpt(3)-1 )
    !
    ! If the BC is on a MAX face, we need to make sure that the beginning
    ! index for this face direction is also converted to a cell index. For
    ! all other directions, the beginning index will be less than the ending
    ! index so this will not affect those indices.
    ! NOTE: For whichever direction, dir, matches the face that the boundary
    !       condition is on:
    !          bcs_info(nb)%bc%begpt(dir) = bcs_info(nb)%bc%endpt(dir)
    !       Thus subtracting one from the endpt value to convert to cell
    !       indices will make one of the following i,j,k loops a zero-trip
    !       loop and nothing will be done. Therefore, we must make sure
    !       that the endpt values are not less than the begpt values.
    !
    ibeg = min(ibeg,iend)
    jbeg = min(jbeg,jend)
    kbeg = min(kbeg,kend)
   !!
   !! Output information about this boundary condition
   !!
   !if (mypnum == glb_root) then
   !  write (iout,1) nb,trim(bc_integer_to_string(bc_type)),iblk, &
   !                 cface(location),ibeg,iend,jbeg,jend,kbeg,kend
   !end if
    !
    ijk1(:) = 0
    ijk1( bcs_info(nb)%bc%dir1 ) = 1
    !
    ijk2(:) = 0
    ijk2( bcs_info(nb)%bc%dir2 ) = 1
    !
    ! Loop through the boundary faces of this boundary condition and
    ! store the necessary information in bface.
    !
    do kfa = kbeg,kend,stride
      do jfa = jbeg,jend,stride
        do ifa = ibeg,iend,stride
          !
          ijkf = [ifa,jfa,kfa] + ijk0
          !
          nc = nc + 1
          !
          ijk = ijkf
          ipg(1) = map3dto1d( ijk , iblk , 1 )
          !
          ijk = ijkf + stride*ijk1
          ipg(2) = map3dto1d( ijk , iblk , 1 )
          !
          np = 0
          if (nr == 3) then
            !
            ijk = ijkf + stride*(ijk1 + ijk2)
            ipg(3) = map3dto1d( ijk , iblk , 1 )
            !
            ijk = ijkf + stride*ijk2
            ipg(4) = map3dto1d( ijk , iblk , 1 )
            !
            do d2 = 0,stride
             !ijk = ijkf + d2*ijk2 - ijk1
              do d1 = 0,stride
               !ijk = ijk + ijk1
                ijk = ijkf + d1*ijk1 + d2*ijk2
                np = np + 1
                ipc(np) = map3dto1d( ijk , iblk , 1 )
              end do
            end do
            !
          else
            !
            do d1 = 0,stride
              ijk = ijkf + d1*ijk1
              np = np + 1
              ipc(np) = map3dto1d( ijk , iblk , 1 )
            end do
            !
          end if
          !
          grid%elem(nc)%prop = element_properties(face_type,stride)
          !
          ! USING F2003 AUTO-REALLOCATION
          !
          grid%elem(nc)%pts = ipc
          grid%elem(nc)%nodes = ipg
          grid%elem(nc)%tags = [nb,0]
          !
          ! Find the host cell for this face
          !
          ! Mark the nodes for this boundary face
          !
          lpoin(ipg(:)) = 1
          !
          n1 = cells_with_node_ptr( ipg(1) )+1
          n2 = cells_with_node_ptr( ipg(1) +1)
          !
          grid%elem(nc)%host_cell = 0
          do ie = n1,n2
            n = cells_with_node(ie)
            if (sum( lpoin(grid%elem(n)%nodes(:)) ) == geom_pts) then
              grid%elem(nc)%host_cell = n
              lpoin(ipg(:)) = 0
            end if
          end do
          if (grid%elem(nc)%host_cell == 0) write (iout,1092) nc-lcell,n,nc
          !
         !bface( 1,nf) = bc_type ! actual boundary condition for this face
         !bface( 2,nf) = host_cell
         !bface( 3,nf) = merge(Geom_Hexa,Geom_Quad,nr==3) ! geom of host cell
         !bface( 4,nf) = merge(HostFace3D(location),HostFace2D(location),nr==3)
         !bface( 5,nf) = nb ! boundary condition group for this boundary face
         !bface(10,nf) = ncell + nf ! ghost cell number
          !
        end do
      end do
    end do
    !
    ! End of loop over the boundary faces
    !
  end do
  call debug_timer(stop_timer,"Creating grid%elem(ncell+1:ncell+nfbnd)")
  !
 !do nc = lcell+1,size(grid%elem)
 !  nf = nc-lcell
 !  nb = grid%elem(nc)%tags(1)
 !  n1 = size(grid%elem(nc)%nodes)
 !  n2 = size(grid%elem(nc)%pts)
 !  write (155,155) nf,nb,bcs_info(nb)%bc%block, &
 !                  (grid%elem(nc)%nodes(n),n=1,n1)
 !  write (156,155) nf,nb,bcs_info(nb)%bc%block, &
 !                  (grid%elem(nc)%pts(n),n=1,n2)
 !end do
 !155 format ("nf = ",i4,", BC #",i2,", BLK #",i1," : ",*(1x,i4,:))
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
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," ERROR TRYING TO USE COARSENING STRIDES WITH PLOT3D!!!",/, &
              " THE NUMBER OF GRID CELLS MUST BE EVENLY",/, &
              " DIVISIBLE BY THE COARSENING STRIDES!!!",/)
  2 format (3x,"Cells not evenly divisible in block ",i0," for direction ",i0)
  3 format (3x,"Faces not evenly divisible for BC #",i0,/, &
            6x,"in block ",i0," for direction ",i0)
  1092 format (5x,"ERROR finding host cell for boundary face ",i10,", ",i10,i10)
  !
end subroutine create_plot3d_grid_cells
!
!###############################################################################
!
subroutine readbcsinfo( bcsfile )
  !
  ! This subs reads boundary conditions information by vulcan input file
  !
  !.. Use Statements ..
  use geovar, only : nr,nfbnd
  use ovar, only : bc_names
  !
  !.. Formal Arguments ..
  character(len=*), intent(inout) :: bcsfile
  !
  !.. Local Scalars ..
  integer :: nb,ierr,dir,dir1,dir2,beg1,end1,beg2,end2,idx1,location,iblk
  integer :: iobcs
  character(len=150) :: text
  !
  !.. Local Allocatable Arrays ..
  character(len=50), allocatable, dimension(:) :: clist
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "readbcsinfo"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Initialize the total number of boundary faces
  !
  nfbnd = 0
  !
  ! Check the name of the bcs file specified in the input file
  ! and alter it if necessary
  !
  call check_supplementary_file(bcsfile,bcs_default)
  !
  ! Open the input file for bcs information
  !
  open (newunit=iobcs,file=bcsfile,form="formatted",status="old", &
                      action="read",position="rewind",iostat=ierr, &
                      iomsg=error_message)
  call io_error(pname,trim(adjustl(bcsfile)),1,__LINE__,__FILE__, &
                ierr,error_message)
  !
  ! Read through the file and see how many lines of
  ! boundary conditions there are.
  ! NOTE: Initialize nbc to -1 to account for the
  !       first line which is just the file header
  !
  nbc = -1
  find_nbc: do
    read (iobcs,1,iostat=ierr,iomsg=error_message) text
    text = trim(adjustl(text))
    if (ierr < 0) then
      exit find_nbc
    else if (ierr > 0) then
      call io_error(pname,trim(adjustl(bcsfile)),-1,__LINE__,__FILE__, &
                    ierr,error_message)
    else if (text(1:1) == "!" .or. len_trim(text) == 0) then
      cycle find_nbc
    else
      nbc = nbc + 1
    end if
  end do find_nbc
  !
  ! Allocate the boundary condition array
  !
  allocate ( bcs_info(1:nbc) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bcs_info",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( bc_names(1:nbc) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bcs_info",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Rewind back to the beginning of the BC file and
  ! skip/read the header line.
  !
  rewind (iobcs)
  !
  skip_header: do
    read (iobcs,1,iostat=ierr,iomsg=error_message) text
    text = trim(adjustl(text))
    if (text(1:1) == "!" .or. len_trim(text) == 0) then
      cycle skip_header
    else
      exit skip_header
    end if
  end do skip_header
 !!
 !read (iobcs,1,iostat=ierr,iomsg=error_message) text
  !
  ! Now actually read the boundary condition information
  ! contained in the BC file
  !
  do nb = 1,nbc
    !
    read_bcs: do
      read (iobcs,1,iostat=ierr) text
      text = trim(adjustl(text))
      if (text(1:1) == "!" .or. len_trim(text) == 0) then
        cycle read_bcs
      else
        exit read_bcs
      end if
    end do read_bcs
   !!
   !read (iobcs,1,iostat=ierr) text
    !
    call parse_list(text,clist,list_min_size=11)
    !
    ! Convert the BC string to its respective integer value
    !
    bcs_info(nb)%bc_type = bc_string_to_integer(clist(1))
    !
    ! Read the block this BC belongs to
    !
    read (clist(2),*,iostat=ierr,iomsg=error_message) iblk
    call io_error(pname,"bcs_info(nb)%block",-1,__LINE__,__FILE__, &
                  ierr,error_message)
    !
    bcs_info(nb)%bc%block = iblk
    !
    ! Get the face location for this BC
    !
    dir = get_direction( clist(3) )
    bcs_info(nb)%bc%dir = dir
    !
    ! Find whether this BC is on a min or max face
    !
    if (clist(4) == "MIN") then
      location = 1
      bcs_info(nb)%bc%begpt(dir) = 1
      bcs_info(nb)%bc%endpt(dir) = 1
    else if (clist(4) == "MAX") then
      location = 2
      bcs_info(nb)%bc%begpt(dir) = ijkpts(dir,iblk)
      bcs_info(nb)%bc%endpt(dir) = ijkpts(dir,iblk)
    else
      read (clist(4),*,iostat=ierr,iomsg=error_message) location
      call io_error(pname,"bcs_info(nb)%bc%begpt(dir)",-1,__LINE__,__FILE__, &
                    ierr,error_message)
      if (location == 1) then
        bcs_info(nb)%bc%begpt(dir) = 1
        bcs_info(nb)%bc%endpt(dir) = 1
      else
        bcs_info(nb)%bc%begpt(dir) = ijkpts(dir,iblk)
        bcs_info(nb)%bc%endpt(dir) = ijkpts(dir,iblk)
      end if
    end if
    !
    bcs_info(nb)%bc%loc = 2*dir-1 + merge(0,1,location==1)
   !write (iout,'("nb = ",i2," ; location = ",i1," ; hostface = ",i1)') &
   !  nb,bcs_info(nb)%bc%loc, &
   !  merge(HostFace3D(bcs_info(nb)%bc%loc), &
   !        HostFace2D(bcs_info(nb)%bc%loc),nr==3)
    !
    ! Fill in the remaining directional information based on 2D/3D grid
    !
    if (nr == 2) then
      !
      ! This is a 2D grid so get the first direction
      !
      dir1 = get_direction( clist(5) )
      idx1 = 6
      !
      ! If the first direction is k, try the second direction
      !
      if (dir1 == 3) then
        dir1 = get_direction( clist(8) )
        idx1 = 9
      end if
      !
      if ( (dir==1 .and. dir1/=2) .or. (dir==2 .and. dir1/=1) ) then
        if (mypnum == glb_root) then
          write (iout,3) nb,dir,dir1
          call stop_gfr(abort,pname,__LINE__,__FILE__)
        end if
        ! Have all the host roots call mpi_barrier here so that the MPI
        ! processes that arent the global root do not continue past this
        ! point before the global root has a chance to call mpi_abort
        ! within the above conditional.
        call mpi_barrier(host_roots_comm,mpierr)
      end if
      !
      beg1 = get_dir_limit( clist(idx1) , ijkpts(dir1,iblk) )
      end1 = get_dir_limit( clist(idx1+1) , ijkpts(dir1,iblk) )
      !
      bcs_info(nb)%bc%dir1 = dir1
      bcs_info(nb)%bc%begpt(dir1) = beg1
      bcs_info(nb)%bc%endpt(dir1) = end1
      !
      ! Default for the second direction which shouldnt be used for 2D
      !
      bcs_info(nb)%bc%dir2 = 3
      bcs_info(nb)%bc%begpt(3) = 1
      bcs_info(nb)%bc%endpt(3) = 1
      !
      ! Store the total number of faces for this boundary condition
      ! and update the total number of boundary faces for the grid.
      !
      bcs_info(nb)%nfaces = end1 - beg1
      !
      nfbnd = nfbnd + bcs_info(nb)%nfaces
      !
    else
      !
      ! This is a 3D grid so get the first direction
      !
      dir1 = get_direction( clist(5) )
      beg1 = get_dir_limit( clist(6) , ijkpts(dir1,iblk) )
      end1 = get_dir_limit( clist(7) , ijkpts(dir1,iblk) )
      !
      bcs_info(nb)%bc%dir1 = dir1
      bcs_info(nb)%bc%begpt(dir1) = beg1
      bcs_info(nb)%bc%endpt(dir1) = end1
      !
      ! Get the second direction
      !
      dir2 = get_direction( clist(8) )
      beg2 = get_dir_limit( clist(9) , ijkpts(dir2,iblk) )
      end2 = get_dir_limit( clist(10) , ijkpts(dir2,iblk) )
      !
      bcs_info(nb)%bc%dir2 = dir2
      bcs_info(nb)%bc%begpt(dir2) = beg2
      bcs_info(nb)%bc%endpt(dir2) = end2
      !
      ! Store the total number of faces for this boundary condition
      ! and update the total number of boundary faces for the grid.
      !
      bcs_info(nb)%nfaces = (end1 - beg1)*(end2 - beg2)
      !
      nfbnd = nfbnd + bcs_info(nb)%nfaces
      !
    end if
    !
    bc_names(nb) = trim(clist(11))
   !if ( bc_names(nb) /= '1' ) bc_names(nb) = ''
    !
  end do
  !
  close (iobcs,iostat=ierr,iomsg=error_message)
  call io_error(pname,trim(adjustl(bcsfile)),2,__LINE__,__FILE__, &
                ierr,error_message)
  !
  ! Deallocate clist if it is still allocated
  !
  if (allocated(clist)) then
    deallocate ( clist , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"clist",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (a)
  3 format (" ERROR READING BCS FILE! THE DIRECTIONS GIVEN ARE NOT VALID!",/, &
            "        BOUNDARY CONDITION = ",i0,/, &
            "        MAIN DIRECTION = ",i0,/, &
            "        FIRST DIRECTION = ",i0,/)
  !
end subroutine readbcsinfo
!
!###############################################################################
!
subroutine readcutinfo( cutfile )
  !
  ! This subs reads cuts information by vulcan input file
  !
  !.. Use Statements ..
  use geovar, only : nr
  !
  !.. Formal Arguments ..
  character(len=*), intent(inout) :: cutfile
  !
  !.. Local Scalars ..
  integer :: n,m,ierr,dir,dir1,dir2,beg1,beg2,end1,end2,idx1,location,iblk
  integer :: iocut
  character(len=150) :: text
  !
  !.. Local Arrays ..
  integer, dimension(1:2) :: cut_points
  !
  !.. Local Allocatable Arrays ..
  character(len=50), allocatable, dimension(:) :: clist
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "readcutinfo"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Check the name of the cut file specified in the input file
  ! and alter it if necessary
  !
  call check_supplementary_file(cutfile,cut_default)
  !
  ! Open the input file for cut information
  !
  open (newunit=iocut,file=cutfile,form="formatted",status="old", &
                      action="read",position="rewind",iostat=ierr, &
                      iomsg=error_message)
  call io_error(pname,trim(adjustl(cutfile)),1,__LINE__,__FILE__, &
                ierr,error_message)
  !
  ! Read through the file and see how many lines of
  ! cut conditions there are.
  ! NOTE: Initialize ncut to -1 to account for the
  !       first line which is just the file header
  !
  ncut = -1
  find_ncut: do
    read (iocut,1,iostat=ierr,iomsg=error_message) text
    text = trim(adjustl(text))
    if (ierr < 0) then
      exit find_ncut
    else if (ierr > 0) then
      call io_error(pname,trim(adjustl(cutfile)),-1,__LINE__,__FILE__, &
                    ierr,error_message)
    else if (text(1:1) == "!" .or. len_trim(text) == 0) then
      cycle find_ncut
    else
      ncut = ncut + 1
    end if
  end do find_ncut
  !
  if (mod(ncut,2) /= 0) then
    write (error_message,4)
    call stop_gfr(abort,pname,__LINE__,__FILE__,error_message)
  end if
  !
  ncut = ncut/2
  !
  ! Allocate the cut condition array
  !
  allocate ( cut_info(1:ncut) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"cut_info",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Rewind back to the beginning of the cut file and
  ! skip/read the header line.
  !
  rewind (iocut)
  !
  skip_header: do
    read (iocut,1,iostat=ierr) text
    text = trim(adjustl(text))
    if (text(1:1) == "!" .or. len_trim(text) == 0) then
      cycle skip_header
    else
      exit skip_header
    end if
  end do skip_header
 !!
 !read (iocut,1,iostat=ierr,iomsg=error_message) text
  !
  ! Now actually read the cut condition information
  ! contained in the cut file
  !
  cut_loop: do n = 1,ncut
    !
    side_of_cut_loop: do m = 1,2
      !
      read_cut: do
        read (iocut,1,iostat=ierr) text
        text = trim(adjustl(text))
        if (text(1:1) == "!" .or. len_trim(text) == 0) then
          cycle read_cut
        else
          exit read_cut
        end if
      end do read_cut
     !!
     !read (iocut,1,iostat=ierr) text
      !
      call parse_list(text,clist,list_min_size=10)
      !
      ! Store the name of the cut
      !
      cut_info(n)%cut_name(m) = clist(1)
      !
      ! Read the block this cut belongs to
      !
      read (clist(2),*,iostat=ierr,iomsg=error_message) iblk
      call io_error(pname,"cut_info(n)%side(m)%block",-1,__LINE__,__FILE__, &
                    ierr,error_message)
      !
      cut_info(n)%side(m)%block = iblk
      !
      ! Get the face location for this cut
      !
      dir = get_direction( clist(3) )
      cut_info(n)%side(m)%dir = dir
      !
      ! Find whether this cut is on a min or max face
      !
      if (clist(4) == "MIN") then
        location = 1
        cut_info(n)%side(m)%begpt(dir) = 1
        cut_info(n)%side(m)%endpt(dir) = 1
      else if (clist(4) == "MAX") then
        location = 2
        cut_info(n)%side(m)%begpt(dir) = ijkpts(dir,iblk)
        cut_info(n)%side(m)%endpt(dir) = ijkpts(dir,iblk)
      else
        read (clist(4),*,iostat=ierr,iomsg=error_message) location
        call io_error(pname,"cut_info(n)%side(m)%begpt(dir)",-1, &
                      __LINE__,__FILE__,ierr,error_message)
        if (location == 1) then
          cut_info(n)%side(m)%begpt(dir) = 1
          cut_info(n)%side(m)%endpt(dir) = 1
        else
          cut_info(n)%side(m)%begpt(dir) = ijkpts(dir,iblk)
          cut_info(n)%side(m)%endpt(dir) = ijkpts(dir,iblk)
        end if
      end if
      !
      cut_info(n)%side(m)%loc = 2*dir-1 + merge(0,1,location==1)
      !
      ! Fill in the remaining directional information based on 2D/3D grid
      !
      if (nr == 2) then
        !
        ! This is a 2D grid so get the first direction
        !
        dir1 = get_direction( clist(5) )
        idx1 = 6
        !
        ! If the first direction is k, try the second direction
        !
        if (dir1 == 3) then
          dir1 = get_direction( clist(8) )
          idx1 = 9
        end if
        !
        if ( (dir==1 .and. dir1/=2) .or. (dir==2 .and. dir1/=1) ) then
          if (mypnum == glb_root) then
            write (iout,3) n,m,dir,dir1
            call stop_gfr(abort,pname,__LINE__,__FILE__)
          end if
          ! Have all the host roots call mpi_barrier here so that the MPI
          ! processes that arent the global root do not continue past this
          ! point before the global root has a chance to call mpi_abort
          ! within the above conditional.
          call mpi_barrier(host_roots_comm,mpierr)
        end if
        !
        beg1 = get_dir_limit( clist(idx1) , ijkpts(dir1,iblk) )
        end1 = get_dir_limit( clist(idx1+1) , ijkpts(dir1,iblk) )
        !
        ! Default for the second direction which shouldnt be used for 2D
        !
        dir2 = 3
        beg2 = 1
        end2 = 1
        !
      else
        !
        ! This is a 3D grid so get the first direction
        !
        dir1 = get_direction( clist(5) )
        beg1 = get_dir_limit( clist(6) , ijkpts(dir1,iblk) )
        end1 = get_dir_limit( clist(7) , ijkpts(dir1,iblk) )
        !
        ! Get the second direction
        !
        dir2 = get_direction( clist(8) )
        beg2 = get_dir_limit( clist(9) , ijkpts(dir2,iblk) )
        end2 = get_dir_limit( clist(10) , ijkpts(dir2,iblk) )
        !
      end if
      !
      ! Store the information for the first direction in cut_info
      !
      cut_info(n)%side(m)%dir1 = dir1
      cut_info(n)%side(m)%begpt(dir1) = beg1
      cut_info(n)%side(m)%endpt(dir1) = end1
      !
      ! Store the information for the second direction in cut_info
      !
      cut_info(n)%side(m)%dir2 = dir2
      cut_info(n)%side(m)%begpt(dir2) = beg2
      cut_info(n)%side(m)%endpt(dir2) = end2
      !
      ! Get the number of points on this side of the cut
      !
      cut_points(m) = (abs(end1-beg1)+1) * (abs(end2-beg2)+1)
      !
    end do side_of_cut_loop
    !
    ! Check that both sides of the cut have the same number of points
    !
    if (cut_points(1) /= cut_points(2)) then
      if (mypnum == glb_root) then
        write (iout,2) trim(cut_info(n)%cut_name(1)), &
                            cut_points(1),cut_points(2)
        call stop_gfr(abort,pname,__LINE__,__FILE__)
      end if
      ! Have all the host roots call mpi_barrier here so that the MPI
      ! processes that arent the global root do not continue past this
      ! point before the global root has a chance to call mpi_abort
      ! within the above conditional.
      call mpi_barrier(host_roots_comm,mpierr)
    else
      cut_info(n)%points = cut_points(1)
    end if
    !
  end do cut_loop
  !
  close (iocut,iostat=ierr,iomsg=error_message)
  call io_error(pname,trim(adjustl(cutfile)),2,__LINE__,__FILE__, &
                ierr,error_message)
  !
  ! Deallocate clist if it is still allocated
  !
  if (allocated(clist)) then
    deallocate ( clist , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"clist",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (a)
  2 format (/," *********************************", &
              "  ERROR!  *********************************",/, &
              " ***  The number of points on each side of cut '", &
                     a,"' do not match! ***",/, &
              " ***          Number of points on side 1 = ",i0,/, &
              " ***          Number of points on side 2 = ",i0,/, &
              " *********************************", &
              "  ERROR!  *********************************",/)
  3 format (" ERROR READING CUT FILE! THE DIRECTIONS GIVEN ARE NOT VALID!",/, &
            "        CUT CONDITION = ",i0,/, &
            "        SIDE OF CUT = ",i0,/, &
            "        MAIN DIRECTION = ",i0,/, &
            "        FIRST DIRECTION = ",i0,/)
  4 format (" ERROR READING CUT FILE! THERE SHOULD BE AN EVEN NUMBER OF CUT", &
            " CONDITIONS IN THE CUT FILE, ONE FOR EACH SIDE OF A CUT!")
  !
end subroutine readcutinfo
!
!###############################################################################
!
subroutine findcutpoints()
  !
  ! This subs reads cuts information by vulcan input file
  !
  !.. Use Statements ..
  use geovar, only : nr,xyz_nodes
  !
  !.. Local Scalars ..
  integer :: i,j,k,l,m,n,np,iblk,ierr
  integer :: ipcut,p1,p2,b1,b2
  integer :: location,dir,dir1,dir2
  integer :: beg1,beg2,end1,end2,stp1,stp2,n1,n2
  real(wp) :: diff
  character(len=100) :: cut_message
  !
  !.. Local Arrays ..
  integer, dimension(1:3) :: ijk,ijk1,ijk2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "findcutpoints"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Loop over the number of cuts
  !
  cut_loop: do n = 1,ncut
    !
    ! Get the number of points on this cut
    !
    np = cut_info(n)%points
    !
    ! Allocate the cut_info(n)%nodes array
    !
    allocate ( cut_info(n)%nodes(1:2,1:np) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"cut_info(n)%nodes",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    ! Loop over the two blocks
    !
    write (cut_message,10) n
    call debug_timer(start_timer,cut_message)
    two_blocks: do m = 1,2
      !
      ! Set the number of cut points
      !
      np = 0
      !
      ! Block numbers
      !
      iblk = cut_info(n)%side(m)%block
      !
      ! Get the block face the cut is on and the
      ! first and second directions of the cut
      !
      dir = cut_info(n)%side(m)%dir
      dir1 = cut_info(n)%side(m)%dir1
      dir2 = cut_info(n)%side(m)%dir2
      !
      ! Get the node index in the primary direction
      !
      ijk(dir) = cut_info(n)%side(m)%begpt(dir)
     !location = mod(cut_info(n)%side(m)%loc+1,2) + 1 ! converts [1,2,3,4,5,6]
     !!                                                       to [1,2,1,2,1,2]
     !ijk(dir) = merge(1,ijkpts(dir,iblk),location==1)
      !
      ! Get the extent of the cut in the first direction
      ! NOTE: stp1 adjusts the following loop for
      !       the situation when beg1 > end1
      !
      beg1 = cut_info(n)%side(m)%begpt(dir1)
      end1 = cut_info(n)%side(m)%endpt(dir1)
      stp1 = (end1-beg1)/abs(end1-beg1)
      !
      ! Get the extent of the cut in the second direction
      ! NOTE: stp2 adjusts the following loop for
      !       the situation when beg2 > end2
      !
      beg2 = cut_info(n)%side(m)%begpt(dir2)
      end2 = cut_info(n)%side(m)%endpt(dir2)
      if (nr == 2) then
        stp2 = 1
      else
        stp2 = (end2-beg2)/abs(end2-beg2)
      end if
      !
      do n1 = beg1,end1,stp1
        do n2 = beg2,end2,stp2
          np = np + 1
          ijk(dir1) = n1
          ijk(dir2) = n2
          cut_info(n)%nodes(m,np) = map3dto1d(ijk,iblk,1)
        end do
      end do
      !
    end do two_blocks
    call debug_timer(stop_timer,cut_message)
    !
    ! Check the consistency of data
    !
    ! Loop over the number of points in this cut
    !
    write (cut_message,11) n
    call debug_timer(start_timer,cut_message)
    do ipcut = 1,cut_info(n)%points
      !
      p1 = cut_info(n)%nodes(1,ipcut)
      p2 = cut_info(n)%nodes(2,ipcut)
      !
      b1 = cut_info(n)%side(1)%block
      b2 = cut_info(n)%side(2)%block
      !
      ijk1(1:3) = map1dto3d(p1,b1,1)
      ijk2(1:3) = map1dto3d(p2,b2,1)
      !
      diff = zero
      do l = 1,nr
        diff = diff + abs(xyz_nodes(l,p1)-xyz_nodes(l,p2)) / &
                     (abs(xyz_nodes(l,p1))+eps7)
      end do
      if (diff > eps2) then
        if (mypnum == glb_root) then
          if (nr == 2) then
            write (iout,1) ipcut, p1, (ijk1(l),l=1,nr), (xyz_nodes(l,p1),l=1,nr)
            write (iout,1) ipcut, p2, (ijk2(l),l=1,nr), (xyz_nodes(l,p2),l=1,nr)
          else
            write (iout,2) ipcut, p1, (ijk1(l),l=1,nr), (xyz_nodes(l,p1),l=1,nr)
            write (iout,2) ipcut, p2, (ijk2(l),l=1,nr), (xyz_nodes(l,p2),l=1,nr)
          end if
        end if
      end if
    end do
    call debug_timer(stop_timer,cut_message)
    !
    ! End of loops over the cuts
    !
  end do cut_loop
  !
  call debug_timer(leaving_procedure,pname)
  !
  1 format (2i8,2i5,4e14.6)
  2 format (2i8,3i5,6e14.6)
  !
  10 format ("Creating cut matches for cut # ",i0)
  11 format ("Checking the consistency of cut matches for cut # ",i0)
  !
end subroutine findcutpoints
!
!###############################################################################
!
subroutine remove_duplicate_nodes()
  !
  !.. Use Statements ..
  use geovar, only : ncell,xyz_nodes,grid,nodes_of_cell
  !
  !.. Local Scalars ..
  integer :: i,ip,ipcut,ip1,ip2,ip3,n,nc,ib1,ib2,ipn1,ipn2,ierr
  integer :: new_pt,old_pt,duplicate_node,l
  integer :: lpts_new,lpts_old
  logical(lk) :: first_line_not_output
  logical(lk) :: output_info
 !real(wp) :: diff
  !
  !.. Local Arrays ..
  integer  :: ijk1(1:4),ijk2(1:4),ijk3(1:4)
  real(wp) :: old_xyz(1:3),new_xyz(1:3)
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: lpoin(:)
  integer,  allocatable :: new2old(:)
  integer,  allocatable :: old2new(:)
  real(wp), allocatable :: temp_xyz(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "remove_duplicate_nodes"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  first_line_not_output = true
  output_info = fals
  !
  ! Save the old number of grid nodes lpts_old
  !
  lpts_old = size(xyz_nodes,dim=2)
  !
  ! Initialize lpoin
  ! lpoin is defined as:
  !    lpoin(ip) > 0  =>  ip is a keeper point
  !    lpoin(ip) < 0  =>  ip is a duplicate point to be removed
  !                       and abs(lpoin(ip)) is the keeper point for
  !                       this duplicate.
  !
  lpoin = intseq(1,lpts_old) ! USING F2003 AUTO-REALLOCATION
  !
  ! Remove the multiple defined points
  !
 !do n = ncut,1,-1
  do n = 1,ncut
    !
    ! Two blocks associated with this cut
    !
    ib1 = cut_info(n)%side(1)%block
    ib2 = cut_info(n)%side(2)%block
    !
    ! Number of grid points on this cut
    !
    do ipcut = 1,cut_info(n)%points
      !
      ! Points on the right and left blocks
      !
      ip1 = cut_info(n)%nodes(1,ipcut)
      ip2 = cut_info(n)%nodes(2,ipcut)
      !
#ifdef DEBUG_CUTS
      if (any(n == [7,8,9])) then
        ijk1 = map1dto3d(ip1,1)
        ijk2 = map1dto3d(ip2,1)
       !if ( (ijk1(3) ==   1 .and. ijk2(3) == 193) .or. &
       !     (ijk1(3) == 193 .and. ijk2(3) ==   1) ) then
       !  output_info = true
        if (n == 7) then
          if (ijk1(1) == 157 .and. ijk2(1) == 157) then
            if (ijk1(3) == 1 .and. ijk2(3) == 193) then
              output_info = true
            end if
          end if
        else if (n == 8) then
          if (ijk1(1) == 157 .and. ijk2(1) == 1) then
            if (any(ijk1(3) == [1,193]) .and. ijk1(3) == ijk2(3)) then
              output_info = true
            end if
          end if
        else if (n == 9) then
          if (ijk1(1) == 1 .and. ijk2(1) == 1) then
            if (ijk1(3) == 1 .and. ijk2(3) == 193) then
              if (ijk1(2) >= 133 .and. ijk1(2) == ijk2(2)) then
                output_info =true
              end if
            end if
          end if
        end if
      end if
#endif
      !
      if (output_info) then
        if (lpoin(ip2) > 0) then
          write (iout,200) ipcut,(ijk1(i),i=1,4),ip1,lpoin(ip1), &
                                 (ijk2(i),i=1,4),ip2,lpoin(ip2)
          write (120+n,202) ipcut,(ijk1(i),i=1,4),ip1,lpoin(ip1), &
                                  (ijk2(i),i=1,4),ip2,lpoin(ip2)
        else
          ip3 = -lpoin(ip2)
          ijk3 = map1dto3d(ip3,1)
          write (iout,200) ipcut,(ijk1(i),i=1,4),ip1,lpoin(ip1), &
                                 (ijk2(i),i=1,4),ip2,lpoin(ip2), &
                                 (ijk3(i),i=1,4),ip3,lpoin(ip3)
          write (120+n,202) ipcut,(ijk1(i),i=1,4),ip1,lpoin(ip1), &
                                  (ijk2(i),i=1,4),ip2,lpoin(ip2), &
                                  (ijk3(i),i=1,4),ip3,lpoin(ip3)
        end if
      end if
      !
      ! Point to the adjacent point
      !
      if (lpoin(ip2) > 0) then
        !
        if (lpoin(ip1) > 0) then
          ! ip2 is a duplicate of ip1 so keep which ever one
          ! is higher and have the other point to it
          if (ip2 > ip1) then
            lpoin(ip2) = -ip1
          else
            lpoin(ip1) = -ip2
          end if
        else
          ! ip1 has already been marked as a duplicate
          ! and already points to a keeper
          lpoin(ip2) = lpoin(ip1)
        end if
        !
      else
        !
        ! Remove the point connected by more than 4 cells
        !
        ip3 = -lpoin(ip2)
        if (ip3 == ip1) then
          if (lpoin(ip1) > 0) then
            lpoin(ip2) = -ip1
          else
            lpoin(ip2) = lpoin(ip1)
          end if
        else if (ip3 > ip1) then
          if (lpoin(ip1) > 0) then
            lpoin(ip3) = -ip1
            lpoin(ip2) = -ip1
          else ! lpoin(ip1) < 0
           !if (lpoin(ip3) > 0) then
           !  lpoin(ip1) = -ip3
           !  lpoin(ip2) = -ip3
           !else
              lpoin(ip3) = lpoin(ip1)
              lpoin(ip2) = lpoin(ip1)
           !end if
          end if
        else
          if (lpoin(ip3) > 0) then
            lpoin(ip1) = -ip3
            lpoin(ip2) = -ip3
          else
            lpoin(ip1) = lpoin(ip3)
            lpoin(ip2) = lpoin(ip3)
          end if
        end if
        !
      end if
      !
      if (output_info) then
        ip3 = abs(lpoin(ip2))
        ijk3 = map1dto3d(ip3,1)
        write (iout,201) ipcut,(ijk1(i),i=1,4),ip1,lpoin(ip1), &
                               (ijk2(i),i=1,4),ip2,lpoin(ip2), &
                               (ijk3(i),i=1,4),ip3,lpoin(ip3)
        write (120+n,203) ipcut,(ijk1(i),i=1,4),ip1,lpoin(ip1), &
                                (ijk2(i),i=1,4),ip2,lpoin(ip2), &
                                (ijk3(i),i=1,4),ip3,lpoin(ip3)
        flush (iout)
        flush (120+n)
      end if
      !
      output_info = fals
      !
    end do
    !
  end do
  200 format (/,"      Before (ipcut=",i0,") :",/, &
                "  ip1 [",3i4," (",i1,") ] : lpoin(",i8,") = ",i0,/, &
                "  ip2 [",3i4," (",i1,") ] : lpoin(",i8,") = ",i0,:,/, &
                "  ip3 [",3i4," (",i1,") ] : lpoin(",i8,") = ",i0)
  201 format (  "       After (ipcut=",i0,") :",/, &
                "  ip1 [",3i4," (",i1,") ] : lpoin(",i8,") = ",i0,/, &
                "  ip2 [",3i4," (",i1,") ] : lpoin(",i8,") = ",i0,/, &
                "  ip3 [",3i4," (",i1,") ] : lpoin(",i8,") = ",i0)
  202 format (/,"Before (ipcut=",i0,") =>",/, &
                " ip1 (",3(i3,","),i1,") : lpoin(",i8,") = ",i9,/, &
                " ip2 (",3(i3,","),i1,") : lpoin(",i8,") = ",i9,:,/, &
                " ip3 (",3(i3,","),i1,") : lpoin(",i8,") = ",i9)
  203 format (  " After (ipcut=",i0,") =>",/, &
                " ip1 (",3(i3,","),i1,") : lpoin(",i8,") = ",i9,/, &
                " ip2 (",3(i3,","),i1,") : lpoin(",i8,") = ",i9,/, &
                " ip3 (",3(i3,","),i1,") : lpoin(",i8,") = ",i9)
  !
  ! new2old(ip) = the old point number for a new point ip
  ! old2new(ip) = the new point number for a old point ip
  !
  lpts_new = count( lpoin(:) > 0 )
  !
  allocate ( new2old(1:lpts_new) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"new2old",1,__LINE__,__FILE__,ierr,error_message)
  allocate ( old2new(1:lpts_old) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"old2new",1,__LINE__,__FILE__,ierr,error_message)
  !
  new_pt = 0
  do old_pt = 1,lpts_old
    if (lpoin(old_pt) > 0) then
      new_pt = new_pt + 1
      new2old(new_pt) = old_pt
      old2new(old_pt) = new_pt
    end if
  end do
  !
  ! Get the new point number for all the old grid points
  !
  do old_pt = 1,lpts_old
    !
    if (lpoin(old_pt) < 0) then
      !
      duplicate_node = -lpoin(old_pt)
      new_pt = old2new(duplicate_node)
      old2new(old_pt) = new_pt
      !
    end if
    !
  end do
  !
  if (mypnum == glb_root) write (iout,1) lpts_old,lpts_new
  !
  ! Check the consistency
  !
  do n = 1,ncut
    !
    ! Number of grid points on this cut
    !
    do ipcut = 1,cut_info(n)%points
      !
      ! Points on the right and left blocks
      !
      ip1 = cut_info(n)%nodes(1,ipcut)
      ip2 = cut_info(n)%nodes(2,ipcut)
      ipn1 = old2new(ip1)
      ipn2 = old2new(ip2)
      !
      if (any([ipn1,ipn2] == 0)) then
        if (first_line_not_output) then
          if (mypnum == glb_root) write (iout,2)
          first_line_not_output = fals
        end if
        if (mypnum == glb_root) then
          write (iout,8) 1,n,ipcut,map1dto3d(ip1,1),ipn1
          write (iout,8) 2,n,ipcut,map1dto3d(ip2,1),ipn2
        end if
      else if (ipn1 /= ipn2) then
        if (first_line_not_output) then
          if (mypnum == glb_root) write (iout,2)
          first_line_not_output = fals
        end if
        if (mypnum == glb_root)  write (iout,9) n, ipcut, ip1, ip2, ipn1, ipn2
      end if
      !
    end do
  end do
  !
 !new_xyz = zero
 !old_xyz = zero
 !!
 !do ip = 1,lpts_old
 !  !
 !  new_pt = old2new(ip)
 !  old_pt = new2old(new_pt)
 !  !
 !  new_xyz(:) = xyz_nodes(:,new_pt)
 !  old_xyz(:) = xyz_nodes(:,old_pt)
 !  !
 !  write (44,3) ip, new_pt, old_pt, lpoin(ip)
 !  !
 !  diff = sum( abs(new_xyz-old_xyz) / (abs(new_xyz)+eps7) )
 !  !
 !  if (diff > eps2) then
 !    write (44,4) ip, old2new(ip), (new_xyz(n), old_xyz(n), n=1,size(new_xyz))
 !  end if
 !  !
 !end do
  !
  ! Deallocate lpoin since we are done with it
  !
  deallocate ( lpoin , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"lpoin",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Adjust xyz_nodes to reflect the nodal set without any duplicates
  !
  ! Allocate the temporary nodal array
  !
  temp_xyz = xyz_nodes(:,new2old(:)) ! F2003 AUTO-REALLOCATION
 !allocate( temp_xyz(1:size(xyz_nodes,dim=1),1:lpts_new) , source=zero , &
 !          stat=ierr , errmsg=error_message )
 !call alloc_error(pname,"temp_xyz",1,__LINE__,__FILE__,ierr,error_message)
 !!
 !! Copy over the non-duplicate nodes from xyz_nodes into temp_xyz
 !!
 !do n = 1,lpts_new
 !  do l = 1,size(xyz_nodes,dim=1)
 !    temp_xyz(l,n) = xyz_nodes(l,new2old(n))
 !  end do
 !end do
  !
  ! Use move_alloc to reallocate xyz_nodes. This does the following:
  !     1. Deallocates xyz_nodes
  !     2. Allocates xyz_nodes with the same type, type
  !        parameters, array bounds, and value as temp_xyz
  !     3. Deallocates temp_xyz
  !
  call move_alloc( from=temp_xyz , to=xyz_nodes )
  !
#ifdef DEBUG_CUTS
  if (mypnum == glb_root) then
    do new_pt = 1,lpts_new
      write (2000,200) new_pt,new2old(new_pt)
    end do
    200 format ("new_point ",i8," -> old_point ",i8)
    flush (2000)
    do old_pt = 1,lpts_old
      write (2001,201) old_pt,old2new(old_pt)
    end do
    201 format ("old_point ",i8," -> new_point ",i8)
    flush (2001)
  end if
  call mpi_barrier(host_roots_comm,mpierr)
#endif
  !
  ! Deallocate new2old
  !
  deallocate ( new2old , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"new2old",2,__LINE__,__FILE__,ierr,error_message)
 !!
 !! Create the nodes_of_cell_ptr and nodes_of_cell arrays
 !!
 !call create_nodes_of_cells(old2new)
  !
  ! Fix all the node information in the grid derived type
#ifdef DEBUG_CUTS
  if (mypnum == glb_root) then
    do nc = 1,size(grid%elem)
      if (allocated(grid%elem(nc)%nodes)) then
        write (1000+mypnum,100) nc,(grid%elem(nc)%nodes(n), &
                                    n=1,size(grid%elem(nc)%nodes))
        write (1002+mypnum,102) nc,(map1dto3d(grid%elem(nc)%nodes(n),1), &
                                    n=1,size(grid%elem(nc)%nodes))
      else
        write (1000+mypnum,101) nc
        write (1002+mypnum,101) nc
      end if
    end do
    flush (1000+mypnum)
    flush (1002+mypnum)
  end if
  call mpi_barrier(host_roots_comm,mpierr)
  100 format (1x,i6," : ",8i9)
  101 format (1x,i6," : NOT ALLOCATED")
  102 format (1x,i6," :",8(2x,3i4," (",i1,")"))
#endif
  !
  do nc = 1,size(grid%elem)
    if (allocated(grid%elem(nc)%nodes)) then
     !do n = 1,size(grid%elem(nc)%nodes)
     !  grid%elem(nc)%nodes(n) = old2new( grid%elem(nc)%nodes(n) )
     !end do
      grid%elem(nc)%nodes(:) = old2new( grid%elem(nc)%nodes(:) )
    end if
    if (allocated(grid%elem(nc)%pts)) then
     !do n = 1,size(grid%elem(nc)%pts)
     !  grid%elem(nc)%pts(n) = old2new( grid%elem(nc)%pts(n) )
     !end do
      grid%elem(nc)%pts(:) = old2new( grid%elem(nc)%pts(:) )
    end if
  end do
  !
  ! Fix the node information in nodes_of_cell
  !
  if (allocated(nodes_of_cell)) then
   !do n = 1,size(nodes_of_cell)
   !  nodes_of_cell(n) = old2new( nodes_of_cell(n) )
   !end do
    nodes_of_cell(:) = old2new( nodes_of_cell(:) )
  end if
  !
  ! Deallocate old2new
  !
  deallocate ( old2new , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"old2new",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," After removing duplicate nodes due to cut", &
              " interfaces of the Plot3D grid:",/, &
              "      old nnode = ",i0,/, &
              "      new nnode = ",i0)
  2 format (/," Error checking the consistency !!!!",/, &
              "   n, ipcut, ip1, ip2, ipn1, ipn2")
  3 format (4i10)
  4 format (i8,i10,6e14.6)
  8 format ("Zero found on side ",i1," of cut:", &
            2x,i4,i8,2x,3i5," (",i0,") => ",i0)
  9 format (i4,i8, 4i12)
  !
end subroutine remove_duplicate_nodes
!
!###############################################################################
!
subroutine reorder_plot3d_nodes
  !
  !.. Use Statements ..
  use geovar, only : ncell,nnode,xyz_nodes,grid
  use geovar, only : nodes_of_cell
  !
  !.. Local Scalars ..
  integer :: n,nc,lpts,ierr
  !
  !.. Local Allocatable Arrays ..
  logical(lk), allocatable :: msk(:)
  integer,     allocatable :: new2old(:)
  integer,     allocatable :: old2new(:)
  real(wp),    allocatable :: temp_xyz(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "reorder_plot3d_nodes"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  lpts = size(xyz_nodes,dim=2)
  !
  allocate ( msk(1:lpts) , source=fals , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",1,__LINE__,__FILE__,ierr,error_message)
  !
#ifdef DEBUG_CUTS
  if (mypnum == glb_root) then
    do nc = 1,size(grid%elem)
      if (allocated(grid%elem(nc)%nodes)) then
        write (1001+mypnum,100) nc,(grid%elem(nc)%nodes(n), &
                                    n=1,size(grid%elem(nc)%nodes))
      else
        write (1001+mypnum,101) nc
      end if
    end do
    flush (1001+mypnum)
  end if
  call mpi_barrier(host_roots_comm,mpierr)
  100 format (1x,i6," : ",8i9)
  101 format (1x,i6," : NOT ALLOCATED")
#endif
  !
  do nc = 1,ncell
    msk( grid%elem(nc)%nodes(:) ) = true
  end do
  !
  nnode = count(msk)
  !
  ! Initialize new2old and old2new
  !
  new2old = [ pack( intseq(1,lpts) , msk ) , & ! USING F2003 AUTO-REALLOCATION
              pack( intseq(1,lpts) , .not. msk ) ]
  !
  ! Deallocate the msk array
  !
  deallocate ( msk , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"msk",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Initialize old2new to new2old to allocate it to the right size
  old2new = new2old ! USING F2003 AUTO-REALLOCATION
  ! Now, create the inverse map of new2old
  old2new( new2old(:) ) = intseq(1,size(new2old))
  !
  call move_alloc( from=xyz_nodes , to=temp_xyz )
  !
  xyz_nodes = temp_xyz(:,new2old(1:nnode))
  !
  temp_xyz(:,:) = temp_xyz(:,new2old(:))
  !
  call move_alloc( from=temp_xyz , to=grid%xyz )
  !
  ! Deallocate new2old
  !
  deallocate ( new2old , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"new2old",2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Fix all the node information in the grid derived type
  !
  do nc = 1,size(grid%elem)
    if (allocated(grid%elem(nc)%nodes)) then
      do n = 1,size(grid%elem(nc)%nodes)
        grid%elem(nc)%nodes(n) = old2new( grid%elem(nc)%nodes(n) )
      end do
    end if
    if (allocated(grid%elem(nc)%pts)) then
      do n = 1,size(grid%elem(nc)%pts)
        grid%elem(nc)%pts(n) = old2new( grid%elem(nc)%pts(n) )
      end do
    end if
  end do
  !
  ! Fix the node information in nodes_of_cell
  !
  if (allocated(nodes_of_cell)) then
    do n = 1,size(nodes_of_cell)
      nodes_of_cell(n) = old2new( nodes_of_cell(n) )
    end do
  end if
  !
  ! Deallocate old2new
  !
  deallocate ( old2new , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"old2new",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine reorder_plot3d_nodes
!
!###############################################################################
!
subroutine create_nodes_of_cells(old2new)
  !
  !.. Use Statements ..
  use order_mod, only : n_order
  use geovar, only : ncell,nr
  use geovar, only : nodes_of_cell_ptr,nodes_of_cell
  use geovar, only : cell_geom,cell_order
  !
  !.. Formal Arguments ..
  integer, intent(in) :: old2new(:)
  !
  !.. Local Scalars ..
  integer :: n,cell_pts,ierr,iblk,nc,npx,npy,npz,i,j,k,n1,n2
  !
  !.. Local Arrays ..
  integer, dimension(1:8) :: ip
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_nodes_of_cells"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate the nodes_of_cell pointer array
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
  ! Set the cell geometry type for all the cells
  !
  cell_geom(:) = merge(Geom_Hexa,Geom_Quad,nr==3)
  !
  ! The number of nodes for all cells will be the same with it being
  ! either 4 (quads) for a 2D grid or 8 (hexes) for a 3D grid.
  !
  cell_pts = merge(8,4,nr==3)
  !
  ! Set up the pointer array for nodes_of_cell
  !
  do n = 2,size(nodes_of_cell_ptr)
    nodes_of_cell_ptr(n) = nodes_of_cell_ptr(n-1) + cell_pts
  end do
  !
  ! Allocate the nodes_of_cell array based on the last pointer value
  !
  allocate ( nodes_of_cell(1:last(nodes_of_cell_ptr)) , source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Initialize the cell counter
  !
  nc = 0
  !
  ! Loop each block and store the nodes for each cell within the block
  ! in the nodes_of_cell array.
  !
  block_loop: do iblk = 1,nblk
    !
    npx = ijkpts(1,iblk)
    npy = ijkpts(2,iblk)
    npz = ijkpts(3,iblk)
    !
    k_loop: do k = 1,merge(1,npz-1,nr==2) ! Make sure the k loop is evaluated at
      j_loop: do j = 1,(npy-1)            ! least once since npz=1 for a 2D grid
        i_loop: do i = 1,(npx-1)
          !
          nc = nc + 1
          !
          ! Get the first four nodes for this cell
          !
          ip(1) = map3dto1d( [i,j,k] , iblk , 1 )
          ip(2) = ip(1) + 1
          ip(3) = ip(2) + npx
          ip(4) = ip(1) + npx
          !
          ! Add the k direction if this is a 3D grid
          !
          if (nr == 3) then
            ip(5:8) = ip(1:4) + npx*npy
          end if
          !
          n1 = nodes_of_cell_ptr(nc)+1
          n2 = nodes_of_cell_ptr(nc+1)
          !
          nodes_of_cell(n1:n2) = old2new( ip(1:cell_pts) )
          !
        end do i_loop
      end do j_loop
    end do k_loop
    !
  end do block_loop
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine create_nodes_of_cells
!
!###############################################################################
!
subroutine getbface()
  !
  !.. Use Statements ..
  use geovar, only : nr,nbfai,ncell,nfbnd,bface,grid
  !
  !.. Local Scalars ..
  integer :: nf,nb,bc_type,iblk,ierr
  integer :: i,j,k,host_cell,dir,location
  integer :: ibeg,iend,jbeg,jend,kbeg,kend
  !
  !.. Local Arrays ..
  integer, dimension(1:3) :: ijk
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "getbface"
  character(len=*), parameter :: cface(1:6) = ["IMIN","IMAX", &
                                               "JMIN","JMAX", &
                                               "KMIN","KMAX"]
  !
  ! The variable loc is in the range 1:4 for 2D and 1:6 in 3D and gives
  ! the face of the block on the boundary corresponding to
  !       loc = [1,2,3,4,5,6] = [imin,imax,jmin,jmax,kmin,kmax]
  !
  ! hostface2d(loc) gives the face of the 2D quad host cell on the boundary
  ! hostface3d(loc) gives the face of the 3D hex  host cell on the boundary
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (allocated(grid)) then
    if (allocated(grid%elem)) then
      call get_bface_from_grid_dt
      call debug_timer(leaving_procedure,pname)
      return
    end if
  end if
  !
  ! Allocate the bface array
  !
  allocate ( bface(1:nbfai,1:nfbnd), source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Initialize the face counter
  !
  nf = 0
  !
  face_loop: do nb = 1,nbc
    !
    ! Boundary condition type for this face
    !
    bc_type = bcs_info(nb)%bc_type
    !
    if (bc_type == bc_periodic) then
      plot3d_has_periodic_faces = true
    end if
    !
    ! Block number for this face
    !
    iblk = bcs_info(nb)%bc%block
    !
    ! Side/Direction (i,j,k) of this BC and its location (min,max)
    !   dir: [1,2,3]       = [i,j,k]
    !   location: [1,2,3,4,5,6] = [imin,imax,jmin,jmax,kmin,kmax]
    !
    dir = bcs_info(nb)%bc%dir
    location = bcs_info(nb)%bc%loc
    !
    ! Beginning index in each direction
    !
    ibeg = bcs_info(nb)%bc%begpt(1)
    jbeg = bcs_info(nb)%bc%begpt(2)
    kbeg = bcs_info(nb)%bc%begpt(3)
    !
    ! Ending index in each direction
    ! NOTE: The index in bcs_info(nb)%bc%endpt(:) gives the node index that
    !       the boundary condition ends at. We have to subtract one from
    !       this index to convert it to a cell index. We also have to make
    !       sure that an end index does not become 0 after subtracting one.
    !       This would happen for whatever index corresponds to the face the
    !       BC is on if it is a MIN face.
    !
    iend = max( 1 , bcs_info(nb)%bc%endpt(1)-1 )
    jend = max( 1 , bcs_info(nb)%bc%endpt(2)-1 )
    kend = max( 1 , bcs_info(nb)%bc%endpt(3)-1 )
    !
    ! If the BC is on a MAX face, we need to make sure that the beginning
    ! index for this face direction is also converted to a cell index. For
    ! all other directions, the beginning index will be less than the ending
    ! index so this will not affect those indices.
    ! NOTE: For whichever direction, dir, matches the face that the boundary
    !       condition is on:
    !          bcs_info(nb)%bc%begpt(dir) = bcs_info(nb)%bc%endpt(dir)
    !       Thus subtracting one from the endpt value to convert to cell
    !       indices will make one of the following i,j,k loops a zero-trip
    !       loop and nothing will be done. Therefore, we must make sure
    !       that the endpt values are not less than the begpt values.
    !
    ibeg = min(ibeg,iend)
    jbeg = min(jbeg,jend)
    kbeg = min(kbeg,kend)
    !
    ! Output information about this boundary condition
    !
    if (mypnum == glb_root) then
      write (iout,1) nb,trim(bc_integer_to_string(bc_type)),iblk, &
                     cface(location),ibeg,iend,jbeg,jend,kbeg,kend
    end if
    !
    ! Loop through the boundary faces of this boundary condition and
    ! store the necessary information in bface.
    !
    do k = kbeg,kend
      do j = jbeg,jend
        do i = ibeg,iend
          !
          nf = nf + 1
          !
          host_cell = map3dto1d( [i,j,k] , iblk , 2 )
          !
          bface( 1,nf) = bc_type ! actual boundary condition for this face
          bface( 2,nf) = host_cell
          bface( 3,nf) = merge(Geom_Hexa,Geom_Quad,nr==3) ! geom of host cell
          bface( 4,nf) = merge(HostFace3D(location),HostFace2D(location),nr==3)
          bface( 5,nf) = nb ! boundary condition group for this boundary face
          bface(10,nf) = ncell + nf ! ghost cell number
          !
        end do
      end do
    end do
    !
    ! End of loop over the boundary faces
    !
  end do face_loop
  !
  if (nfbnd /= nf) then
    if (mypnum == glb_root) then
      write (iout,2) nfbnd,nf
      call stop_gfr(abort,pname,__LINE__,__FILE__)
    end if
    ! Have all the host roots call mpi_barrier here so that the MPI processes
    ! that arent the global root do not continue past this point before the
    ! global root has a chance to call mpi_abort within the above conditional.
    call mpi_barrier(host_roots_comm,mpierr)
  end if
  !
 !do nf = 1,nfbnd
 !  write (110,3) trim(bc_integer_to_string(bface(1,nf))), &
 !                (bface(i,nf),i=2,4),bface(10,nf)
 !end do
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format (/," BC #",i0," - ",a,": on Block ",i0," on the ",a," face.",/, &
              " IMIN/MAX = ",i0,"/",i0, &
             ", JMIN/MAX = ",i0,"/",i0, &
             ", KMIN/MAX = ",i0,"/",i0)
  2 format (/," Error creating the boundary face array!",/, &
              " The number of expected boundary faces = ",i0,/, &
              " The number of created  boundary faces = ",i0,/)
  3 format (a," - hc=",i0,", geo=",i1,", location=",i1,", gc=",i0)
  !
end subroutine getbface
!
!###############################################################################
!
subroutine get_bface_from_grid_dt()
  !
  !.. Use Statements ..
  use geovar, only : nr,nbfai,ncell,nfbnd,bface,grid
  !
  !.. Local Scalars ..
  integer :: nb,nc,nf,ierr,location,host_cell
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_bface_from_grid_dt"
  !
  ! hostface2d(loc) gives the face of the 2D quad host cell on the boundary
  ! hostface3d(loc) gives the face of the 3D hex  host cell on the boundary
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nfbnd = size(grid%elem) - ncell
  !
  ! Allocate the bface array
  !
  allocate ( bface(1:nbfai,1:nfbnd), source=0 , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
  !
  do nf = 1,nfbnd
    !
    nc = nf + ncell
    !
    host_cell = grid%elem(nc)%host_cell
    !
    nb = grid%elem(nc)%tags(1)
    !
    location = bcs_info(nb)%bc%loc
    !
    if (bcs_info(nb)%bc_type == bc_periodic) then
      plot3d_has_periodic_faces = true
    end if
    !
    bface( 1,nf) = bcs_info(nb)%bc_type
    bface( 2,nf) = host_cell
    bface( 3,nf) = grid%elem(host_cell)%prop%geom
    bface( 4,nf) = merge(HostFace3D(location),HostFace2D(location),nr==3)
    bface( 5,nf) = nb
    bface(10,nf) = nc
    !
    grid%elem(nc)%tags(1) = bface(1,nf)
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine get_bface_from_grid_dt
!
!###############################################################################
!
subroutine match_plot3d_periodic_boundaries()
  !
  !.. Use Statements ..
  use geovar, only : nr,nfbnd,bface,xyz_nodes
  use geovar, only : nodes_of_cell_ptr,nodes_of_cell
  !
  !.. Local Scalars ..
  integer :: n,npf,nf,pf,p1,p2,ierr
  integer :: host_cell,host_geom,host_side
  integer :: nmin,nmax,pmin,pmax,nfnods
  real(wp) :: dscoef
  !
  !.. Local Arrays ..
  integer, dimension(1:4) :: ip
  integer, dimension(1:6) :: hostface
  integer, dimension(1:6) :: pfaces
  integer, dimension(1:2,1:3) :: idx
  !
  !.. Local Allocatable Arrays ..
  logical(lk), allocatable :: face_msk(:,:)
  logical(lk), allocatable :: matched(:)
  integer, allocatable :: minfac(:)
  integer, allocatable :: maxfac(:)
  real(wp), allocatable :: ds(:)
  real(wp), allocatable :: xyz_min(:,:)
  real(wp), allocatable :: xyz_max(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "match_plot3d_periodic_boundaries"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nr == 3) then
    hostface = HostFace3D
    dscoef  = one
    idx(1:2,1:3) = reshape( [2,3,1,3,1,2] , [2,3] )
  else
    hostface = HostFace2D
    dscoef  = one / sqrt(two)
    idx(1:2,1:3) = reshape( [2,2,1,1,0,0] , [2,3] )
  end if
  !
  allocate ( face_msk(1:nfbnd,1:6) , source=fals , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_msk",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Count the number of periodic faces for each direction
  !
  do n = 1,6
    face_msk(:,n) = ( bface(1,:)==bc_periodic .and. bface(4,:)==hostface(n) )
  end do
  !
  pfaces(1:6) = count( face_msk(:,1:6) , dim=1 )
  !
  ! Check to make sure the number of faces in each direction is consistent
  !
  if (any( pfaces([1,3,5]) /= pfaces([2,4,6]) )) then
    if (mypnum == glb_root) then
      if (pfaces(1) /= pfaces(2)) write (iout,1) "X",pfaces(1),pfaces(2)
      if (pfaces(3) /= pfaces(4)) write (iout,1) "Y",pfaces(3),pfaces(4)
      if (pfaces(5) /= pfaces(6)) write (iout,1) "Z",pfaces(5),pfaces(6)
      call stop_gfr(abort,pname,__LINE__,__FILE__)
    end if
    ! Have all the host roots call mpi_barrier here so that the MPI processes
    ! that arent the global root do not continue past this point before the
    ! global root has a chance to call mpi_abort within the above conditional.
    call mpi_barrier(host_roots_comm,mpierr)
  end if
  !
  !
  !
  do n = 1,nr
    !
    nmin = (n-1)*2 + 1
    nmax = nmin + 1
    !
    npf = pfaces(nmin)
    !
    allocate ( minfac(1:npf) , source=pack(intseq(1,nfbnd),face_msk(:,nmin)) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"minfac",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( maxfac(1:npf) , source=pack(intseq(1,nfbnd),face_msk(:,nmax)) , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"maxfac",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( xyz_min(1:nr,1:npf) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_min",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( xyz_max(1:nr,1:npf) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_max",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( matched(1:npf) , source=fals , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"matched",1,__LINE__,__FILE__,ierr,error_message)
    !
    allocate ( ds(1:npf) , source=zero , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ds",1,__LINE__,__FILE__,ierr,error_message)
    !
    !
    !
    do nf = 1,npf
      !
      ! Get the coordinates for the centers of all
      ! the minimum faces in the current direction
      !
      pf = minfac(nf)
      !
      host_cell = bface(2,pf)
      host_geom = bface(3,pf)
      host_side = bface(4,pf)
      !
      nfnods = num_face_nodes(host_side,host_geom)
      !
      p1 = nodes_of_cell_ptr(host_cell)+1
      p2 = nodes_of_cell_ptr(host_cell+1)
      !
      ip(1:nfnods) = side_nodes(1:nfnods,host_side,host_geom)
      !
      ip(1:nfnods) = nodes_of_cell( p1 + ip(1:nfnods) )
      !
      xyz_min(1:nr,nf) = sum( xyz_nodes(1:nr,ip(1:nfnods)) , dim=2 )
      xyz_min(1:nr,nf) = xyz_min(1:nr,nf) / real(nfnods,kind=wp)
      !
      ! Now for all the max faces
      !
      pf = maxfac(nf)
      !
      host_cell = bface(2,pf)
      host_geom = bface(3,pf)
      host_side = bface(4,pf)
      !
      nfnods = num_face_nodes(host_side,host_geom)
      !
      p1 = nodes_of_cell_ptr(host_cell)+1
      p2 = nodes_of_cell_ptr(host_cell+1)
      !
      ip(1:nfnods) = side_nodes(1:nfnods,host_side,host_geom)
      !
      ip(1:nfnods) = nodes_of_cell( p1 + ip(1:nfnods) )
      !
      xyz_max(1:nr,nf) = sum( xyz_nodes(1:nr,ip(1:nfnods)) , dim=2 )
      xyz_max(1:nr,nf) = xyz_max(1:nr,nf) / real(nfnods,kind=wp)
      !
    end do
    !
    !
    !
    do nf = 1,npf
      !
      ds(:) = huge(zero)
      !
      do pf = 1,npf
        !
        if (matched(pf)) cycle
        !
        ds(pf) = dscoef * norm2( xyz_min(idx(:,n),nf) - &
                                 xyz_max(idx(:,n),pf) )
        !
      end do
      !
      pf = minloc( ds , dim=1 )
      matched(pf) = true
      !
      pmin = minfac( nf )
      pmax = maxfac( pf )
      !
      bface(5,pmin) = pmax
      bface(6,pmin) = bface(2,pmax)
      bface(8,pmin) = bface(3,pmax)
      bface(9,pmin) = bface(4,pmax)
      !
      bface(5,pmax) = pmin
      bface(6,pmax) = bface(2,pmin)
      bface(8,pmax) = bface(3,pmin)
      bface(9,pmax) = bface(4,pmin)
      !
    end do
    !
    !
    !
    deallocate ( minfac , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"minfac",2,__LINE__,__FILE__,ierr,error_message)
    !
    deallocate ( maxfac , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"maxfac",2,__LINE__,__FILE__,ierr,error_message)
    !
    deallocate ( xyz_min , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_min",2,__LINE__,__FILE__,ierr,error_message)
    !
    deallocate ( xyz_max , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_max",2,__LINE__,__FILE__,ierr,error_message)
    !
    deallocate ( matched , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"matched",2,__LINE__,__FILE__,ierr,error_message)
    !
    deallocate ( ds , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ds",2,__LINE__,__FILE__,ierr,error_message)
    !
  end do
  !
  deallocate ( face_msk , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_msk",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
  1 format (/," The number of periodic faces in the ",a, &
              " direction is not consistent!",/, &
              "     Periodic faces on MIN boundary: ",i0,/, &
              "     Periodic faces on MAX boundary: ",i0,/)
  !
end subroutine match_plot3d_periodic_boundaries
!
!###############################################################################
!
subroutine create_grid_dt()
  !
  !.. Use Statements ..
  use geovar, only : grid
  use geovar, only : ncell,cell_geom
  use geovar, only : nfbnd,bface
  use geovar, only : nodes_of_cell,nodes_of_cell_ptr
  use geovar, only : nr,xyz_nodes
  !
  !.. Local Scalars ..
  integer :: nc,nf,n1,n2,ierr
  integer :: host_cell,host_geom
  integer :: host_side,nfnods
  !
  !.. Local Arrays ..
  integer :: ip(1:4)
  !
  !.. Local Allocatable Arrays ..
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_grid_dt"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (allocated(grid)) then
    if (allocated(grid%elem)) then
      if (allocated(grid%xyz)) then
        call debug_timer(leaving_procedure,pname)
        return
      end if
    end if
  end if
  !
  grid%xyz = xyz_nodes(:,:) ! F2003 AUTO-REALLOCATION
  !
  allocate ( grid%elem(1:ncell+nfbnd) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"grid%elem",1,__LINE__,__FILE__,ierr,error_message)
  !
  do nc = 1,ncell
    !
    n1 = nodes_of_cell_ptr(nc)+1
    n2 = nodes_of_cell_ptr(nc+1)
    !
    grid%elem(nc) = create_element(cell_geom(nc),nodes_of_cell(n1:n2))
    !
    grid%elem(nc)%tags = [0,0] ! F2003 AUTO-REALLOCATION
    !
  end do
  !
  do nf = 1,nfbnd
    !
    nc = ncell + nf
    !
    host_cell = bface(2,nf)
    host_geom = bface(3,nf)
    host_side = bface(4,nf)
    nfnods = num_face_nodes(host_side,host_geom)
    !
    n1 = nodes_of_cell_ptr(host_cell)+1
    n2 = nodes_of_cell_ptr(host_cell+1)
    !
    ip(1:nfnods) = side_nodes(1:nfnods,host_side,host_geom)
    !
    grid%elem(nc) = create_element(merge(Geom_Quad,Geom_Edge,nr==3), &
                                   nodes_of_cell( n1 + ip(1:nfnods) ))
    !
    grid%elem(nc)%host_cell = host_cell
    !
    grid%elem(nc)%tags = [bface(1,nf),0] ! F2003 AUTO-REALLOCATION
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine create_grid_dt
!
!###############################################################################
!
pure function create_element(elem_type,pts) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : elem_t
  !
  !.. Formal Arguments ..
  integer, intent(in) :: elem_type
  integer, intent(in) :: pts(:)
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
  ! Reorder the nodes and points
  !
  return_value = reorder_plot3d_pts(return_value)
  !
end function create_element
!
!###############################################################################
!
elemental function reorder_plot3d_pts(elem,reverse) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : elem_t
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
      return_value%pts = reorder_plot3d_edge_pts(elem%pts)
    case (Geom_Quad)
      return_value%pts = reorder_plot3d_quad_pts(elem%pts,reverse)
    case (Geom_Hexa)
      return_value%pts = reorder_plot3d_hexa_pts(elem%pts,reverse)
  end select
  !
end function reorder_plot3d_pts
!
!###############################################################################
!
pure function reorder_plot3d_edge_pts(pts) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: pts(:)
  !
  !.. Function Result ..
  integer :: return_value(1:size(pts))
  !
continue
  !
  return_value(:) = pts(:)
  !
end function reorder_plot3d_edge_pts
!
!###############################################################################
!
pure function reorder_plot3d_quad_pts(pts,reverse) result(return_value)
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
  logical(lk) :: reverse_reorder
  !
  !.. Local Parameters ..
  integer, parameter :: quad4(4) = [1,2,4,3]
  !
continue
  !
  reverse_reorder = fals
  if (present(reverse)) then
    reverse_reorder = reverse
  end if
  !
  if (reverse_reorder) then
    return_value( quad4(:) ) = pts(:)
  else
    return_value = pts( quad4(:) )
  end if
  !
end function reorder_plot3d_quad_pts
!
!###############################################################################
!
pure function reorder_plot3d_hexa_pts(pts,reverse) result(return_value)
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
  logical(lk) :: reverse_reorder
  !
  !.. Local Parameters ..
  integer, parameter :: hexa8(8) = [1,2,4,3,5,6,8,7]
  !
continue
  !
  reverse_reorder = fals
  if (present(reverse)) then
    reverse_reorder = reverse
  end if
  !
  if (reverse_reorder) then
    return_value( hexa8(:) ) = pts(:)
  else
    return_value = pts( hexa8(:) )
  end if
  !
end function reorder_plot3d_hexa_pts
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
  call mpi_bcast(array_sizes,n,mpi_inttyp,glb_root,MPI_COMM_WORLD,mpierr)
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
  call debug_timer(entering_procedure,pname)
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
  call debug_timer(leaving_procedure,pname)
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
    !===========================================================================
    !
    if (allocated(bface)) then
      deallocate ( bface , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"bface",2,__LINE__,__FILE__,ierr,error_message)
    end if
    allocate ( bface(1:nbfai,1:nfbnd) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bface",1,__LINE__,__FILE__,ierr,error_message)
    !
    !===========================================================================
    !
    if (allocated(nodes_of_cell_ptr)) then
      deallocate ( nodes_of_cell_ptr , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"nodes_of_cell_ptr",2,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
    allocate ( nodes_of_cell_ptr(1:ncell+1) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodes_of_cell_ptr",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    !===========================================================================
    !
    if (allocated(nodes_of_cell)) then
      deallocate ( nodes_of_cell , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"nodes_of_cell",2,__LINE__,__FILE__, &
                       ierr,error_message)
    end if
    allocate ( nodes_of_cell(1:nocsz) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodes_of_cell",1,__LINE__,__FILE__,ierr, &
                     error_message)
    !
    !===========================================================================
    !
    if (allocated(xyz_nodes)) then
      deallocate ( xyz_nodes , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"xyz_nodes",2,__LINE__,__FILE__,ierr,error_message)
    end if
    allocate ( xyz_nodes(1:nr,1:nnode) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xyz_nodes",1,__LINE__,__FILE__,ierr,error_message)
    !
    !===========================================================================
    !
    if (allocated(bc_names)) then
      deallocate ( bc_names , stat=ierr , errmsg=error_message )
      call alloc_error(pname,"bc_names",2,__LINE__,__FILE__,ierr,error_message)
    end if
    allocate ( bc_names(1:nbcnm) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"bc_names",1,__LINE__,__FILE__,ierr,error_message)
    !
    !===========================================================================
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
pure function element_properties(elem_type,elem_order) result(return_value)
  !
  !.. Use Statements ..
  use geovar, only : prop_t
  !
  !.. Formal Arguments ..
  integer, intent(in) :: elem_type
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: elem_order
  !
  !.. Function Result ..
  type(prop_t) :: return_value
  !
  !.. Local Scalars ..
  integer :: np,order
  !
continue
  !
  order = 1
  if (present(elem_order)) order = elem_order
  !
  select case (elem_type)
    case (Geom_Edge)
      np = order + 1
     !return_value = prop_t(npts=np, &
     !                      order=order, &
     !                      geom=Geom_Edge, &
     !                      elem_type=elem_type)
      return_value%npts = np
      return_value%order = order
      return_value%geom = Geom_Edge
      return_value%elem_type = elem_type
    case (Geom_Quad)
      np = (order+1)**2
     !return_value = prop_t(npts=np, &
     !                      order=order, &
     !                      geom=Geom_Quad, &
     !                      elem_type=elem_type)
      return_value%npts = np
      return_value%order = order
      return_value%geom = Geom_Quad
      return_value%elem_type = elem_type
    case (Geom_Hexa)
      np = (order+1)**3
     !return_value = prop_t(npts=np, &
     !                      order=order, &
     !                      geom=Geom_Hexa, &
     !                      elem_type=elem_type)
      return_value%npts = np
      return_value%order = order
      return_value%geom = Geom_Hexa
      return_value%elem_type = elem_type
    case default
     !return_value = prop_t(npts=-1, &
     !                      order=-1, &
     !                      geom=Geom_Unknown, &
     !                      elem_type=int(elem_type))
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
subroutine writetecplot()
  !
  !.. Use Statements ..
  use geovar, only : nr,ncell,nnode,xyz_nodes
  use geovar, only : nodes_of_cell_ptr,nodes_of_cell
  !
  !.. Local Scalars ..
  integer :: l,n,n1,n2,ierr
  integer :: iotec
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "writetecplot"
  character(len=*), parameter :: fname = "unst_plot3d.tec.dat"
  !
continue
  !
  ! Open a file for an unstructured grid in tecplot format
  !
  open (newunit=iotec,file=fname,form="formatted",status="replace", &
                      action="write",position="rewind",iostat=ierr, &
                      iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Write the Tecplot file header
  !
  if (nr == 2) then
    write (iotec,10) nnode,ncell
  else
    write (iotec,11) nnode,ncell
  end if
  !
  ! Write the nodal coordinates
  !
  do l = 1,nr
    write (iotec,21) (xyz_nodes(l,n), n=1,nnode)
  end do
  !
  ! Write the connectivity matrix
  !
  do n = 1,ncell
    n1 = nodes_of_cell_ptr(n)+1
    n2 = nodes_of_cell_ptr(n+1)
    write (iotec,22) (nodes_of_cell(l), l=n1,n2)
  end do
  !
  ! Close the Tecplot file
  !
  close (iotec,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Format Statements
  !
  10 format ('TITLE = "Unstructured Grid from Plot2D Grid"',/, &
             'VARIABLES = "X", "Y"',/, &
             'ZONE T="Z1", N=',i0,', E=',i0, &
             ', DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL')
  11 format ('TITLE = "Unstructured Grid from Plot3D Grid"',/, &
             'VARIABLES = "X", "Y", "Z"',/, &
             'ZONE T="Z1", N=',i0,', E=',i0, &
             ', DATAPACKING=BLOCK, ZONETYPE=FEBRICK')
  21 format (5e18.10)
  22 format (i0,7(1x,i0))
  !
end subroutine writetecplot
!
!###############################################################################
!
subroutine check_supplementary_file(filename,file_default)
  !
  !.. Formal Arguments ..
  character(len=*),    intent(in) :: file_default
  character(len=*), intent(inout) :: filename
  !
  !.. Local Scalars ..
  character(len=len(grid_folder)+len(filename)) :: test_file
  logical(ldk) :: found
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "check_supplementary_file"
  !
continue
  !
  filename = trim(adjustl(filename))
  !
  ! The input argument filename will only contain a string if one was
  ! actually provided within the input file. If a file name was given,
  ! check to see if the input file exists first before trying other
  ! places.
  !
  if (len_trim(filename) > 0) then
    !
    inquire (file=filename,exist=found)
    !
    ! If the input file exists, return from this subroutine
    ! without making any alterations to the variable filename
    !
    if (found) return
    !
    ! The input file name doesnt seem to exist so try to find
    ! this file within the same folder as the grid file
    !
    test_file = trim(grid_folder) // trim(filename)
    !
    inquire (file=test_file,exist=found)
    !
    if (found) then
      !
      filename = trim(adjustl(test_file))
      !
      return
      !
    end if
    !
  end if
  !
  ! If we have reached this point, either no file was specified in the
  ! input file or the file specified could not be found as specified
  ! or within the same folder as the grid file. Either way, reset
  ! filename to the default value and look for the default file name
  ! in the grid directory or the current working directory.
  !
  ! First check the same directory containing the grid file
  !
  test_file = trim(grid_folder) // trim(file_default)
  !
  inquire (file=test_file,exist=found)
  !
  ! If this file exists, set filename to this
  ! file and return from this subroutine
  !
  if (found) then
    !
    filename = trim(adjustl(test_file))
    !
    return
    !
  end if
  !
  ! If that didnt work, just set filename to the default value and let
  ! the subroutine for reading the file spit out an error if it cant
  ! find the file.
  !
  filename = file_default
  !
end subroutine check_supplementary_file
!
!###############################################################################
!
subroutine parse_list(strng,clist,list_min_size)
  !
  ! Routine to parse a character string (STRNG) to extract a list of characters
  !
  ! STRNG : character string to be parsed
  ! CLIST : list of character strings
  ! LIST_MIN_SIZE : optional argument that makes sure the clist array contains
  !                 at least list_min_size number of array elements
  !
  !.. Formal Arguments ..
  character(len=*), intent(inout) :: strng
  character(len=*), allocatable, dimension(:), intent(inout) :: clist
  integer, optional, intent(in) :: list_min_size
  !
  !.. Local Scalars ..
  integer :: ibeg,is,m,n,num_str,ierr
  integer :: num_rpt,new_size,old_size
  logical(lk) :: repeats_exist
  character(len=len(clist)) :: substr_before,substr_after
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable, dimension(:) :: ast_loc
  !
  !.. Local Parameters ..
  character(len=*), parameter :: delimiters(1:4) = ["'",'"',",",achar(9)]
  character(len=*), parameter :: pname = "parse_list"
  !
continue
  !
  ! Replace all single quotes, double quotes, commas, and tabs with spaces
  !
  forall (n=1:len_trim(strng),any(strng(n:n)==delimiters)) strng(n:n) = " "
  !
  m = 0
  !
  if (len_trim(strng) > 0) then
    is = 1
    find_next_word: do
      !
      if (strng(is:is) /= " ") then
        !
        ibeg = is
        !
        find_end_of_word: do
          if (is == len_trim(strng)) then
            call reallocate( clist , m+1 )
            clist(m+1) = strng(ibeg:is)
          else if (strng(is:is) == " ") then
            m = m + 1
            call reallocate( clist , m )
            clist(m) = strng(ibeg:is-1)
            ibeg = is
            exit find_end_of_word
          end if
          is = is + 1
        end do find_end_of_word
        !
      else
        !
        is = is + 1
        !
      end if
      !
      if (is > len_trim(strng)) exit find_next_word
      !
    end do find_next_word
  end if
  !
  ! Get the final number of sub-strings found in strng
  !
  num_str = size(clist)
  !
  ! Make sure that everything in clist is adjusted left and trimmed.
  ! Probably not necessary but the extra overhead is basically zero
  ! and worth the extra piece of mind.
  !
  do n = 1,num_str
    clist(n) = trim(adjustl(clist(n)))
  end do
  !
  ! Do a quick scan through and see if we need to worry about repeat counts
  !
  repeats_exist = .false.
  do n = 1,num_str
    if (scan(trim(clist(n)),"*") > 0) then
      repeats_exist = .true.
      exit
    end if
  end do
  !
  ! If repeats exist, we need to expand them, adjusting clist accordingly
  !
  if (repeats_exist) then
    !
    ! Allocate the array giving the location of the asterisk
    ! in each clist element and the array giving the number
    ! of repeat counts for those elements
    !
    allocate ( ast_loc(1:num_str) , source=0 , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ast_loc",1,__LINE__,__FILE__,ierr,error_message)
    !
    ! Loop through the elements of clist and find the location of any asterisks
    !
    do n = 1,num_str
      ast_loc(n) = scan( trim(clist(n)) , "*" )
    end do
    !
    ! Now loop back through the list in reverse order to expand the
    ! repeats. By going in reverse order, we dont have to worry about
    ! the value of num_str increasing as we keep expanding the repeats.
    !
    new_size = num_str
    !
    expand_repeats: do n = num_str,1,-1
      !
      if (ast_loc(n) /= 0) then
        !
        ! This element of clist has an asterisk so extract the
        ! substrings before and after its location within the string.
        !
        substr_before = clist(n)(1:ast_loc(n)-1)
        substr_after = clist(n)(ast_loc(n)+1:len_trim(clist(n)))
        !
        ! Try to read the substring before the asterisk. If the substring
        ! cant be read as an integer, assume that this element of clist
        ! is intended to be of character type and this asterisk is not
        ! meant as a repeat counter.
        !
        read (substr_before,*,iostat=ierr) num_rpt
        !
        ! If the substring was successfully read, expand the repeat count.
        !
        if (ierr == 0) then
          !
          ! If num_rpt <= 0 ignore this repeat counter
          !
          if (num_rpt <= 0) cycle expand_repeats
          !
          ! We need to subtract one from the repeat count because one of
          ! the repeated values already exists within clist. Get the
          ! new size of the clist array and reallocate it accordingly.
          !
          num_rpt = num_rpt-1
          old_size = new_size
          new_size = old_size + num_rpt
          !
          call reallocate ( clist , new_size )
          !
          ! Move everything after clist(n) to
          ! the end of the newly expanded array
          !
          do m = old_size,n+1,-1
            clist(num_rpt+m) = clist(m)
          end do
          !
          ! Copy the substring that is to be repeated to the new array elements
          !
          clist(n:n+num_rpt) = substr_after
          !
        end if
        !
      end if
      !
    end do expand_repeats
    !
    ! Deallocate ast_loc since we are done expanding the repeat counts
    !
    deallocate ( ast_loc , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ast_loc",2,__LINE__,__FILE__,ierr,error_message)
    !
  end if
  !
  ! Update the number of elements in clist
  !
  num_str = size(clist)
  !
  ! Make sure that clist is big enough for the bcs and cut subroutines
  !
  if (present(list_min_size)) then
    if (num_str < list_min_size) then
      call reallocate( clist , list_min_size )
      do n = num_str+1,list_min_size
        clist(n) = "1"
      end do
    end if
  end if
  !
  ! Finished parsing a character string (STRNG) to extract a list of characters
  !
  1 format ("size(clist) = ",i0,"; new_size = ",i0)
  !
end subroutine parse_list
!
!###############################################################################
!
pure function string_is_integer(strng) result(return_value)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: strng
  !
  !.. Function Result ..
  logical(lk) :: return_value
  !
  !.. Local Scalars ..
  integer :: n,lenstr,ierr,ival
  character(len=len(strng)) :: loc_strng
  character(len=20) :: read_format
  !
  !.. Local Allocatable Arrays ..
  logical(lk), allocatable, dimension(:) :: is_a_number
  !
  !.. Local Parameters ..
  character(len=*), parameter :: nums_and_signs = "0123456789+-"
  character(len=*), parameter :: nums_only = "0123456789"
  !
continue
  !
  ! Ensure the input string has no leading blanks
  !
  loc_strng = trim(adjustl(strng))
  !
  ! ################
  ! ### Option 1 ###
  ! ################
  !
  write (read_format,'(a,i0,a)') '(i', len(loc_strng), ')'
  read (loc_strng,read_format,iostat=ierr) ival
  !
  return_value = (ierr == 0)
 !!
 !! ################
 !! ### Option 2 ###
 !! ################
 !!
 !! Get the length of the (adjusted) input string
 !!
 !lenstr = len_trim(loc_strng)
 !!
 !! First check that the second to last characters within loc_strng are numbers
 !!
 !return_value = (verify(loc_strng(2:lenstr),nums_only) == 0)
 !!
 !! If the previous test was true, make sure that the first character is
 !! a number or sign to finalize the verification that this is an integer
 !!
 !if (return_value) then
 !  return_value = (verify(loc_strng(1:lenstr),nums_and_signs) == 0)
 !end if
 !!
 !! ################
 !! ### Option 3 ###
 !! ################
 !!
 !! Allocate the is_a_number logical array
 !!
 !allocate ( is_a_number(1:lenstr) , source=fals , &
 !           stat=ierr , errmsg=error_message )
 !if (ierr /= 0) then
 !  return_value = fals
 !  return
 !end if
 !!
 !! Allow for the first character in string to contain a +,- sign
 !!
 !is_a_number(1) = (verify(loc_strng(1:1),nums_and_signs) == 0)
 !!
 !! Check that the remaining characters are numbers
 !!
 !do n = 2,lenstr
 !  is_a_number(n) = (verify(loc_strng(n:n),nums_only) == 0)
 !end do
 !!
 !! Finally make sure that all the elements of is_a_number are true
 !!
 !return_value = all(is_a_number)
 !!
 !! Deallocate is_a_number before leaving the function
 !!
 !deallocate ( is_a_number , stat=ierr , errmsg=error_message )
  !
  ! Ignore if there is an error deallocating is_a_number since
  ! the function result should have been correctly computed
  !
end function string_is_integer
!
!###############################################################################
!
function get_direction(strng) result(return_value)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: strng
  !
  !.. Function Result ..
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: dir,ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_direction"
  !
continue
  !
  dir = index( "IJK" , trim(uppercase(strng)) )
  !
  if (dir == 0) then
    read (strng,*,iostat=ierr,iomsg=error_message) dir
    call io_error(pname,"strng",-1,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  return_value = dir
  !
end function get_direction
!
!###############################################################################
!
function get_dir_limit(strng,mxpts) result(return_value)
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: strng
  integer,          intent(in) :: mxpts
  !
  !.. Function Result ..
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: ipt,ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_dir_limit"
  !
continue
  !
  if (uppercase(strng) == "MIN") then
    return_value = 1
  else if (uppercase(strng) == "MAX") then
    return_value = mxpts
  else
    read (strng,*,iostat=ierr,iomsg=error_message) ipt
    call io_error(pname,"strng",-1,__LINE__,__FILE__,ierr,error_message)
    if (ipt == 0) then
      return_value = mxpts
    else
      return_value = min(ipt,mxpts)
    end if
  end if
  !
end function get_dir_limit
!
!###############################################################################
!
pure function map3dto1d(ijk,n,offset_type) result(return_value)
  !
  ! Map 2D-i,j or 3D-i,j,k structured indices into 1D array space
  !
  !.. Formal Arguments ..
  integer,               intent(in) :: n   ! current block
  integer, dimension(:), intent(in) :: ijk ! i,j,k indices in current block
  integer,               intent(in) :: offset_type
  !
  !.. Function Result ..
  integer :: return_value
  !
  !.. Local Scalars ..
  integer :: l,nr,adj
  !
continue
  !
  nr = size(ijk)
  !
  ! offset_type determines whether we are offsetting or not,
  ! and if we are, whether the offset is for nodes or cells
  !     = 1 : offset for nodes
  !     = 2 : offset for cells
  !     = anything else : no offset
  !
  if (any(offset_type == [1,2])) then
    return_value = offset(offset_type,n)
  else
    return_value = 1
  end if
  !
  ! If we are offsetting for cells, we need to adjust the argument to
  ! the product intrinsic in the next loop.
  !
  if (offset_type == 2) then
    adj = 1
  else
    adj = 0
  end if
  !
  ! Compute the 1D index in the 1D array space given the
  ! 3D index-direction sizes, indices, and block offset
  !
  do l = 1,nr
    return_value = return_value + product(ijkpts(1:l-1,n)-adj) * (ijk(l)-1)
  end do
  !
end function map3dto1d
!
!###############################################################################
!
pure function map1dto3d_known_block(ijk1d,n,offset_type) result(return_value)
  !
  ! Map 2D-i,j or 3D-i,j,k structured indices into 1D array space
  !
  !.. Formal Arguments ..
  integer, intent(in) :: n     ! current block
  integer, intent(in) :: ijk1d ! 1D value of i,j,k indices
  integer, intent(in) :: offset_type
  !
  !.. Function Result ..
  integer, dimension(1:3) :: return_value
  !
  !.. Local Scalars ..
  integer :: blk_offset,ipts,jpts,numer,denom,i,j,k
  !
continue
  !
  ! offset_type determines whether we are offsetting or not,
  ! and if we are, whether the offset is for nodes or cells
  !     = 1 : offset for nodes
  !     = 2 : offset for cells
  !     = anything else : no offset
  !
  if (any(offset_type == [1,2])) then
    blk_offset = offset(offset_type,n)
  else
    blk_offset = 1
  end if
  !
  ! If we are offsetting for cells, we need to adjust the number of points
  ! in each direction to reflect the number of cells in each direction
  !
  ipts = ijkpts(1,n)
  jpts = ijkpts(2,n)
  if (offset_type == 2) then
    ipts = ipts - 1
    jpts = jpts - 1
  end if
  !
  ! k-index
  !
  numer = ijk1d - blk_offset + ipts*jpts
  denom = ipts*jpts
  k = floor( real(numer,kind=wp) / real(denom,kind=wp) )
  !
  ! j-index
  !
  numer = ijk1d - blk_offset + ipts*jpts*(1-k) + ipts
  denom = ipts
  j = floor( real(numer,kind=wp) / real(denom,kind=wp) )
  !
  ! i-index
  !
  i = ijk1d - blk_offset + ipts*jpts*(1-k) + ipts*(1-j) + 1
  !
  ! Put the i,j,k values in the function result
  !
  return_value = [i,j,k]
  !
end function map1dto3d_known_block
!
!###############################################################################
!
pure function map1dto3d_unknown_block(ijk1d,offset_type) result(return_value)
  !
  ! Map 2D-i,j or 3D-i,j,k structured indices into 1D array space
  !
  !.. Formal Arguments ..
  integer, intent(in) :: ijk1d ! 1D value of i,j,k indices
  integer, intent(in) :: offset_type
  !
  !.. Function Result ..
  integer, dimension(1:4) :: return_value
  !
  !.. Local Scalars ..
  integer :: blk_offset,ipts,jpts,numer,denom,i,j,k,l,m,n
  !
continue
  !
  m = offset_type
  if (all(offset_type /= [1,2])) m = 1
  !
  ! Find which block this point is in
  n = size(offset,dim=2)
  do l = 1,size(offset,dim=2)-1
    if (ijk1d < offset(m,l+1)) then
      n = l
      exit
    end if
  end do
  !
  ! offset_type determines whether we are offsetting or not,
  ! and if we are, whether the offset is for nodes or cells
  !     = 1 : offset for nodes
  !     = 2 : offset for cells
  !     = anything else : no offset
  !
  if (any(offset_type == [1,2])) then
    blk_offset = offset(offset_type,n)
  else
    blk_offset = 1
  end if
  !
  ! If we are offsetting for cells, we need to adjust the number of points
  ! in each direction to reflect the number of cells in each direction
  !
  ipts = ijkpts(1,n)
  jpts = ijkpts(2,n)
  if (offset_type == 2) then
    ipts = ipts - 1
    jpts = jpts - 1
  end if
  !
  ! k-index
  !
  numer = ijk1d - blk_offset + ipts*jpts
  denom = ipts*jpts
  k = floor( real(numer,kind=wp) / real(denom,kind=wp) )
  !
  ! j-index
  !
  numer = ijk1d - blk_offset + ipts*jpts*(1-k) + ipts
  denom = ipts
  j = floor( real(numer,kind=wp) / real(denom,kind=wp) )
  !
  ! i-index
  !
  i = ijk1d - blk_offset + ipts*jpts*(1-k) + ipts*(1-j) + 1
  !
  ! Put the i,j,k values in the function result
  !
  return_value = [i,j,k,n]
  !
end function map1dto3d_unknown_block
!
!###############################################################################
!
subroutine plot3d_memory_usage(iunit)
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
  if (allocated(bcs_info)) then
    !
    call dt%set( storage_size(bcs_info,kind=inttype) , &
                         size(bcs_info,kind=inttype) )
    !
    do i = lbound(bcs_info,dim=1),ubound(bcs_info,dim=1)
      !
      call dt%add( storage_size(bcs_info(i)%bc_type,kind=inttype) )
      call dt%add( storage_size(bcs_info(i)%nfaces,kind=inttype) )
      call dt%add( storage_size(bcs_info(i)%bc,kind=inttype) )
      call dt%add( storage_size(bcs_info(i)%bc%block,kind=inttype) )
      call dt%add( storage_size(bcs_info(i)%bc%dir,kind=inttype) )
      call dt%add( storage_size(bcs_info(i)%bc%dir1,kind=inttype) )
      call dt%add( storage_size(bcs_info(i)%bc%dir2,kind=inttype) )
      call dt%add( storage_size(bcs_info(i)%bc%loc,kind=inttype) )
      call dt%add( storage_size(bcs_info(i)%bc%begpt,kind=inttype) , &
                           size(bcs_info(i)%bc%begpt,kind=inttype) )
      call dt%add( storage_size(bcs_info(i)%bc%endpt,kind=inttype) , &
                           size(bcs_info(i)%bc%endpt,kind=inttype) )
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "bcs_info"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "bcs_info"
    !
  end if
  !
  if (allocated(cut_info)) then
    !
    call dt%set( storage_size(cut_info,kind=inttype) , &
                         size(cut_info,kind=inttype) )
    !
    do i = lbound(cut_info,dim=1),ubound(cut_info,dim=1)
      !
      call dt%add( storage_size(cut_info(i)%points,kind=inttype) )
      call dt%add( storage_size(cut_info(i)%cut_name,kind=inttype) , &
                           size(cut_info(i)%cut_name,kind=inttype) )
      call dt%add( storage_size(cut_info(i)%side,kind=inttype) , &
                           size(cut_info(i)%side,kind=inttype) )
      !
      do j = lbound(cut_info(i)%side,dim=1),ubound(cut_info(i)%side,dim=1)
        call dt%add( storage_size(cut_info(i)%side(j)%block,kind=inttype) )
        call dt%add( storage_size(cut_info(i)%side(j)%dir,kind=inttype) )
        call dt%add( storage_size(cut_info(i)%side(j)%dir1,kind=inttype) )
        call dt%add( storage_size(cut_info(i)%side(j)%dir2,kind=inttype) )
        call dt%add( storage_size(cut_info(i)%side(j)%loc,kind=inttype) )
        call dt%add( storage_size(cut_info(i)%side(j)%begpt,kind=inttype) , &
                             size(cut_info(i)%side(j)%begpt,kind=inttype) )
        call dt%add( storage_size(cut_info(i)%side(j)%endpt,kind=inttype) , &
                             size(cut_info(i)%side(j)%endpt,kind=inttype) )
      end do
      !
      call dt%add( storage_size(cut_info(i)%nodes,kind=inttype) , &
                           size(cut_info(i)%nodes,kind=inttype) )
      !
    end do
    !
    call total%add(dt)
    write (array_name,5) "cut_info"
    write (iou,2) array_name,dt%get_units()
    !
  else
    !
    write (iou,7) "cut_info"
    !
  end if
  !
  if (allocated(offset)) then
    call array%set( storage_size(offset,kind=inttype) , &
                            size(offset,kind=inttype) )
    call total%add(array)
    write (array_name,5) "offset"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "offset"
  end if
  if (allocated(ijkpts)) then
    call array%set( storage_size(ijkpts,kind=inttype) , &
                            size(ijkpts,kind=inttype) )
    call total%add(array)
    write (array_name,5) "ijkpts"
    write (iou,2) array_name,array%get_units()
  else
    write (iou,7) "ijkpts"
  end if
  !
  total = total%get_units()
  write (iou,3) total%rmem,trim(adjustl(total%units))
  !
  ! Format Statements
  !
  1 format (/,"Status of arrays in module module_plot3d")
  2 format (" Array: ",a,dt)
  3 format (" Total memory = ",f8.2,1x,a)
  4 format ("cut_info(",a,",",i0,")%nodes")
  5 format ("'",a,"'")
  7 format (" Array: ",a36," -       NOT ALLOCATED")
  !
#endif
end subroutine plot3d_memory_usage
!
!###############################################################################
!
include "Functions/bc_string_to_integer.f90"
include "Functions/last.f90"
!
!###############################################################################
!
end module module_plot3d
!
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
