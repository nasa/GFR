!###############################################################################
!> Module for point probes output.
!>
!> Usage:
!> - Call collect_probe_data to collect and occasionallcollect and occasionally.
!>   Will allocate/initialize local probes and open necessary output files.
!> - Call finalize_probes to close probe files and deallocate memory.
!###############################################################################
module probes_mod
  !
  !.. Use Statements ..
  use module_kind_types
  !
  implicit none
  private
  !
  ! .. Public Procedures ..
  public  :: initialize_probes
  public  :: collect_probe_data
  public  :: finalize_probes
  !
  ! .. Private Procedures ..
  private :: probe_file_open
  private :: probe_file_close
  !
  ! .. Module Parameters ..
  integer, parameter :: nvarmax = 7
  !
  ! .. Module Derived Types ..
  !
  ! PROBE_INFO_T : derived type containing index and coordinate info for
  !                an individual probe
  !
  type :: old_probe_info_t
    ! rake_idx: index of the rake from the input file
    integer :: rake_idx = 0
    ! cell_idx: index of the cell from the input file
    integer :: cell_idx = 0
    ! probe_idx:
    !     If rake_idx >  0, this is the index of the probe within the rake
    !     If cell_idx >  0, this is the index of the probe within the cell
    !     Otherwise,        this is the point_probe index from the input file
    integer :: probe_idx = 0
    ! xyz: coordinates for this probe in meters, must divide this by
    !      grid_scaling_units to get the user grid units which were given
    !      in the input file for this probe.
    real(wp) :: xyz(1:3) = [zero,zero,zero]
  end type old_probe_info_t
  !
  type :: probe_info_t
    ! probe_idx: main probe index based on the input probe arrays
    integer :: probe_idx = 0
    ! sub_idx: sub-index of the probe within a rake or cell probe
    integer :: sub_idx = 0
    ! probe_type: 1=point, 2=rake, 3=cell
    integer :: probe_type
    ! xyz: coordinates for this probe in meters, must divide this by
    !      grid_scaling_units to get the user grid units which were given
    !      in the input file for this probe.
    real(wp) :: xyz(1:3) = [zero,zero,zero]
  end type probe_info_t
  !
  ! MY_PROBES_T : Derived type containing info for an
  !               individual probe on the local partition
  ! 
  type :: my_probes_t
    type(probe_info_t) :: probe_info
    integer :: solpt_idx = 0 ! sol pt on local partition
    integer :: io_unit = 0
    logical(lk) :: ignore_this_probe = fals
  end type my_probes_t
  !
  ! MY_PROBES : Derived type array containing info for
  !             all probes on the local partition
  !
  type(my_probes_t), allocatable, save :: my_probes(:)
  !
  ! MY_SAVED_DATA_T : Derived type to store each data collection for a probe
  !                   so it can be written out at the probes_write_interval
  !
  type :: my_saved_data_t
    integer,  allocatable :: iter(:)
    real(wp), allocatable :: time(:)
    real(wp), allocatable :: sol(:,:)
  end type my_saved_data_t
  !
  ! MY_SAVED_DATA : Derived type array containing the saved data
  !                 collections for all probes on the local partition
  !
  type(my_saved_data_t), allocatable, save :: my_saved_data(:)
  !
  character(len=:), allocatable, save :: probes_dir
  !
  integer, save :: num_total_probes = 0
  integer, save :: num_point_probes = 0
  integer, save :: num_rake_probes  = 0
  integer, save :: num_cell_probes  = 0
  !
  integer, save :: current_save_idx = 0
  !
contains
!
!###############################################################################
!> Initialize memory and file io for local point probes.
!###############################################################################
subroutine initialize_probes()
  !
  use geomvar_mod, only : nr,xyz,cell
  use ovar,        only : point_probes, rake_probes, cell_probes
  use ovar,        only : probes_collection_interval, probes_write_interval
  use io_mod,      only : create_results_dir
  !
  !.. Local Variables ..
  integer  :: io,i,j,k,n,np,npts,n0
  integer  :: this_cell,this_sp,this_idx
  integer  :: my_num_probes
  real(wp) :: ds
  !
  real(wp), allocatable :: xyz_beg(:)
  real(wp), allocatable :: xyz_end(:)
  !
  type(probe_info_t), allocatable, save :: pp_info(:)
  type(probe_info_t), allocatable, save :: rp_info(:)
  type(probe_info_t), allocatable, save :: cp_info(:)
  !
  type(my_probes_t), allocatable :: my_point_probes(:)
  type(my_probes_t), allocatable :: my_rake_probes(:)
  type(my_probes_t), allocatable :: my_cell_probes(:)
  !
  type(my_probes_t), allocatable :: tmp_probes(:)
  !
  character(len=1000) :: temp_dir
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_probes"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  num_point_probes = 0
  num_rake_probes = 0
  num_cell_probes = 0
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  if (allocated(point_probes)) then
    if (size(point_probes) > 0) then
      num_point_probes = size(point_probes)
    end if
  end if
  !
  if (allocated(rake_probes)) then
    do i = 1,size(rake_probes)
      if (rake_probes(i)%num_probes > 0) then
        num_rake_probes = num_rake_probes + rake_probes(i)%num_probes
      end if
    end do
  end if
  !
  if (allocated(cell_probes)) then
    if (size(cell_probes) > 0) then
      num_cell_probes = size(cell_probes)
    end if
  end if
  !
  num_total_probes = num_point_probes + num_rake_probes + num_cell_probes
  !
  ! If there are no probes, set the probe collection interval to zero so
  ! it will stop trying to collect/write probe data and exit this subroutine
  !
  if (num_total_probes < 1) then
    num_total_probes = 0
    probes_collection_interval = 0
    call debug_timer(leaving_procedure,pname)
    return
  end if
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Continuing on since there is at least one probe
  !
  ! Create directory to write all the probe data
  !
  probes_dir = "" ! F2003 AUTO-REALLOC
  if (size(point_probes) > 0) then
    call create_results_dir("Probe_files",temp_dir,skip_lustre=true)
    probes_dir = trim(adjustl(temp_dir)) ! F2003 AUTO-REALLOC
  end if
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  allocate ( pp_info(1:num_point_probes) )
  !
  ! Add any point probes to pp_info
  !
  do n = 1,num_point_probes
    pp_info(n)%probe_idx = n
    pp_info(n)%probe_type = 1
    pp_info(n)%xyz(1:nr) = point_probes(n)%xyz(1:nr)
  end do
  !
  call find_local_probes(pp_info,my_point_probes)
  !
  if (allocated(pp_info)) deallocate ( pp_info )
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  allocate ( rp_info(1:num_rake_probes) )
  !
  ! Add the info for each individual probe from each rake to rp_info
  !
  n = 0
  loop_over_rakes: &
  do i = 1,size(rake_probes)
    !
    np = rake_probes(i)%num_probes
    !
    xyz_beg = rake_probes(i)%xyz_beg(1:nr) ! F2003 AUTO-REALLOC
    xyz_end = rake_probes(i)%xyz_end(1:nr) ! F2003 AUTO-REALLOC
    !
    if (np > 1) then
      !
      ! Loop over the probes on this rake to get the coordinates for each one
      !
      do j = 1,np
        n = n + 1
        ds = real(j-1,kind=wp)/real(np-1,kind=wp)
        rp_info(n)%probe_idx = i
        rp_info(n)%sub_idx = j
        rp_info(n)%probe_type = 2
        rp_info(n)%xyz(1:nr) = xyz_beg(1:nr) + ds*(xyz_end(1:nr)-xyz_beg(1:nr))
      end do
      !
    else if (np == 1) then
      !
      n = n + 1
      !
      ! If only one probe is requested, average the beginning and end
      ! coordinates for the rake to get the probe coordinates.
      !
      rp_info(n)%probe_idx = i
      rp_info(n)%sub_idx = 1
      rp_info(n)%probe_type = 2
      rp_info(n)%xyz(1:nr) = half*(xyz_beg(1:nr) + xyz_end(1:nr))
      !
   !NOTE: skip if np <= 0
    end if
    !
  end do loop_over_rakes
  !
  if (allocated(xyz_beg)) deallocate ( xyz_beg )
  if (allocated(xyz_end)) deallocate ( xyz_end )
  !
  call find_local_probes(rp_info,my_rake_probes)
  !
  if (allocated(rp_info)) deallocate ( rp_info )
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  allocate ( cp_info(1:num_cell_probes) )
  !
  ! Add any cell probes to cp_info
  !
  do n = 1,num_cell_probes
    cp_info(n)%probe_idx = n
    cp_info(n)%probe_type = 3
    cp_info(n)%xyz(1:nr) = cell_probes(n)%xyz(1:nr)
  end do
  !
  call find_local_probes(cp_info,my_cell_probes)
  !
  if (allocated(cp_info)) deallocate ( cp_info )
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Find any duplicate point probes and mark them to be ignored
  !
  do n = 2,size(my_point_probes)
    this_sp = my_point_probes(n)%solpt_idx
    if (any(this_sp == my_point_probes(1:n-1)%solpt_idx)) then
      my_point_probes(n)%ignore_this_probe = true
      write (iout,1) my_point_probes(n)%probe_info%probe_idx
    end if
  end do
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Find any duplicate rake probes and mark them to be ignored
  !
  do n = 2,size(my_rake_probes)
    this_sp = my_rake_probes(n)%solpt_idx
    if (any(this_sp == my_rake_probes(1:n-1)%solpt_idx)) then
      my_rake_probes(n)%ignore_this_probe = true
      write (iout,2) my_rake_probes(n)%probe_info%sub_idx, &
                     my_rake_probes(n)%probe_info%probe_idx
    end if
  end do
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Find any duplicate cell probes and mark them to be ignored
  !
  npts = 0
  do n = 1,size(my_cell_probes)
    this_idx = my_cell_probes(n)%probe_info%probe_idx
    call get_cell_info(my_cell_probes(n)%solpt_idx,this_cell,np)
    my_cell_probes(n)%probe_info%probe_idx = this_cell
    if (any(this_cell == my_cell_probes(1:n-1)%probe_info%probe_idx)) then
      my_cell_probes(n)%ignore_this_probe = true
      write (iout,3) this_idx
    else
      npts = npts + np
    end if
  end do
  !
  ! Now expand out individual point probes for all the solution points
  ! of the cell probes
  !
  allocate ( tmp_probes(1:npts) )
  !
  j = 0
  do n = 1,size(my_cell_probes)
    !
    if (my_cell_probes(n)%ignore_this_probe) cycle
    !
    this_cell = my_cell_probes(n)%probe_info%probe_idx
    n0 = cell(this_cell)%beg_sp - 1
    np = cell(this_cell)%end_sp - n0
    !
    do k = 1,np
      j = j + 1
      tmp_probes(j)%probe_info%probe_idx = this_cell
      tmp_probes(j)%probe_info%sub_idx = k
      tmp_probes(j)%probe_info%probe_type = 3
      tmp_probes(j)%probe_info%xyz(1:nr) = xyz(1:nr,n0+k)
      tmp_probes(j)%solpt_idx = n0+k
    end do
    !
  end do
  !
  call move_alloc( from=tmp_probes , to=my_cell_probes )
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  my_point_probes(:)%probe_info%probe_type = 1
  my_rake_probes(:)%probe_info%probe_type = 2
  my_cell_probes(:)%probe_info%probe_type = 3
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Count the total number of probes
  !
  my_num_probes = size(my_point_probes) + &
                  size(my_rake_probes) + &
                  size(my_cell_probes)
  !
  ! Collect the details of all local point, rake, and cell probes to my_probes
  !
  allocate ( my_probes(1:my_num_probes) )
  !
  j = 0
  do n = 1,size(my_point_probes)
    j = j + 1
    my_probes(j) = my_point_probes(n)
  end do
  !
  do n = 1,size(my_rake_probes)
    j = j + 1
    my_probes(j) = my_rake_probes(n)
  end do
  !
  do n = 1,size(my_cell_probes)
    j = j + 1
    my_probes(j) = my_cell_probes(n)
  end do
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Open all the probe files
  !
  do j = 1,size(my_probes)
    if (my_probes(j)%ignore_this_probe) cycle
    call probe_file_open(j,io)
    my_probes(j)%io_unit = io
  end do
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Allocate the my_saved_data array so it can be used to save the probe
  ! data at the collection interval for writing out at the write interval
  !
  current_save_idx = 0
  !
  allocate ( my_saved_data(1:size(my_probes)) )
  !
  ! NOTE: EVERYTHING IN THIS LOOP IS DOING AN F2003 AUTO-REALLOC
  do i = 1,size(my_saved_data)
    my_saved_data(i)%iter = init_array(probes_write_interval,0)
    my_saved_data(i)%time = init_array(probes_write_interval,zero)
    my_saved_data(i)%sol  = init_array(nr+2,probes_write_interval,zero)
  end do
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1 format ( 2x,"Warning: POINT_PROBES(",i0,") will be ignored because the",/, &
            11x,"solution point identified for this probe was",/, &
            11x,"already identified by another point probe.")
  2 format ( 2x,"Warning: Probe #",i0," of RAKE_PROBES(",i0,") will be ", &
                "ignored because",/, &
            11x,"the solution point identified for this probe was already",/, &
            11x,"identified by another rake probe.")
  3 format ( 2x,"Warning: CELL_PROBES(",i0,") will be ignored because the ", &
                "cell identified for",/, &
            11x,"this probe was already identified by another cell probe.")
  !
end subroutine initialize_probes
!
!###############################################################################
!
subroutine find_local_probes(probe_info,local_probes)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr,n_solpts,xyz
  use ovar,        only : grid_scaling_factor
  !
  !.. Formal Arguments ..
  type(probe_info_t),              intent(in) :: probe_info(:)
  type(my_probes_t), allocatable, intent(out) :: local_probes(:)
  !
  !.. Local Scalars ..
  integer  :: i,j,n,np,num_probes,min_idx
  real(wp) :: min_dist2,dist2
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: prb_xyz(:)
  integer,  allocatable :: local_sp_idx(:)
  real(wp), allocatable :: pairs(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "find_local_probes"
  !
continue
  !
  num_probes = size(probe_info)
  !
  if (num_probes < 1) return
  !
  ! Allocate the arrays needed to find the closest solution point for each
  ! probe and which partition that solution point belongs to.
  !
  local_sp_idx = init_array(num_probes,0)    ! F2003 AUTO-REALLOC
  pairs = init_array(2,num_probes,-99.0_wp)  ! F2003 AUTO-REALLOC
  !
  ! Find local closest sol pt to each probe.
  !
  do i = 1,size(probe_info)
    !
    min_idx = 0            ! sol pt index with min_dist2
    min_dist2 = huge(one)  ! min distance^2 from probe to sol pt
    !
    prb_xyz = probe_info(i)%xyz(1:nr)*grid_scaling_factor ! F2003 AUTO-REALLOC
    !
    do n = 1, n_solpts
      ! Convert probe location from input grid units to meters.
      dist2 = norm2( prb_xyz(1:nr) - xyz(1:nr,n) )
      if ( dist2 < min_dist2 ) then
        min_dist2 = dist2
        min_idx = n
      end if
    end do
    local_sp_idx(i) = min_idx  ! min solpt index
    pairs(1,i) = min_dist2               ! min dist^2
    pairs(2,i) = real(mypnum,kind=wp)    ! mpi process
    !
  end do
  !
  ! Find global closest sol pt to each probe and its mpi process
  !
  if (ncpu > 1) then
    call mpi_allreduce(MPI_IN_PLACE,pairs,size(pairs,dim=2),mpi_2flttyp, &
                       MPI_MINLOC,MPI_COMM_WORLD,mpierr)
  end if
  !
  np = count( nint(pairs(2,:)) == mypnum )
  !
  allocate ( local_probes(1:np) )
  !
  j = 0
  do i = 1,num_probes
    if ( nint(pairs(2,i)) == mypnum ) then
      j = j + 1
      local_probes(j)%probe_info = probe_info(i)
      local_probes(j)%solpt_idx = local_sp_idx(i)
    end if
  end do
  !
  if (allocated(prb_xyz))      deallocate ( prb_xyz )
  if (allocated(local_sp_idx)) deallocate ( local_sp_idx )
  if (allocated(pairs))        deallocate ( pairs )
  !
end subroutine find_local_probes
!
!###############################################################################
!
subroutine get_cell_info(this_solpt,this_cell,npts)
  !
  !.. Use Statements ..
  use geomvar_mod, only : ncell,cell
  !
  !.. Formal Arguments ..
  integer, intent(in)  :: this_solpt
  integer, intent(out) :: this_cell
  integer, intent(out) :: npts
  !
  !.. Local Scalars ..
  integer :: n,nc,k
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_cell_info"
  !
continue
  !
  find_cell_loop: &
  do nc = 1,ncell
    if (this_solpt <= cell(nc)%end_sp) then
      this_cell = nc
      npts = cell(nc)%end_sp - cell(nc)%beg_sp + 1
      exit find_cell_loop
    end if
  end do find_cell_loop
  !
end subroutine get_cell_info
!
!###############################################################################
!> Finalize memory and file io for local point probes.
!###############################################################################
subroutine finalize_probes()
  !
  !.. Local Variables ..
  integer :: i
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "finalize_probes"
  !
continue
  !
  if (num_total_probes < 1) return
  !
  call debug_timer(entering_procedure,pname)
  !
  if (current_save_idx > 0) then
    call write_saved_probe_data
  end if
  !
  do i = 1, size(my_probes)
    if (my_probes(i)%ignore_this_probe) cycle
    call probe_file_close(i)
  end do
  !
  deallocate( my_probes )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine finalize_probes
!
!###############################################################################
!> Collect data for all local probes.
!>
!> If this interval of collecting data is equal to probes_write_interval,
!> write out the saved probe data.
!> Must call finalize_probes to close all probe output files.
!###############################################################################
subroutine collect_probe_data(iter,time,usp)
  !
  !.. Use Statements ..
  use eqn_idx_mod, only : nec,nmb,nme,nee
  use geomvar_mod, only : nr
  use ovar,        only : time_ref,rhoref,aref
  use ovar,        only : probes_write_interval
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: iter      !< solution iteration number
  real(wp), intent(in) :: time      !< solution time [nondim]
  real(wp), intent(in) :: usp(:,:)  !< solution conserved variables [nondim]
  !
  !.. Local Variables ..
  integer  :: i, j, io, ipt, n
  real(wp) :: rhoa, rhoa2
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "collect_probe_data"
  !
continue
  !
  if (num_total_probes < 1) return
  !
  call debug_timer(entering_procedure,pname)
  !
  ! non-dim reference values
 !rhoref                 ! density
  rhoa  = rhoref * aref  ! momentum
  rhoa2 = rhoa * aref    ! energy
  !
  current_save_idx = current_save_idx + 1
  !
  do i = 1,size(my_probes)
    !
    if (my_probes(i)%ignore_this_probe) cycle
    !
    ipt = my_probes(i)%solpt_idx
    !
    my_saved_data(i)%iter(current_save_idx) = iter
    my_saved_data(i)%time(current_save_idx) = time*time_ref
    !
    my_saved_data(i)%sol(1     ,current_save_idx) = usp(ipt,nec)*rhoref
    my_saved_data(i)%sol(2:nr+1,current_save_idx) = usp(ipt,nmb:nme)*rhoa
    my_saved_data(i)%sol(nr+2  ,current_save_idx) = usp(ipt,nee)*rhoa2
    !
  end do
  !
  ! Write out all the saved probe data if the current save index is equal
  ! to the write interval. After writing the data for each probe, reset the
  ! my_saved_data components to zero. After writing the data for all probes,
  ! reset current_save_idx to 0.
  !
  if (current_save_idx == probes_write_interval) then
    call write_saved_probe_data
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine collect_probe_data
!
!###############################################################################
!
! Subroutine to write out all the saved probe data in my_saved_data.
! NOTE: This subroutine will also reset the saved probe data to zero after it
!       has been written out. It will also reset the current save index counter
!       to 0.
subroutine write_saved_probe_data()
  !
  !.. Local Scalars ..
  integer :: i,j,n,ni,np,io
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_saved_probe_data"
  !
continue
  !
  loop_over_probes: &
  do i = 1,size(my_probes)
    !
    if (my_probes(i)%ignore_this_probe) cycle loop_over_probes
    !
    io = my_probes(i)%io_unit
    !
    np = size(my_saved_data(i)%sol,dim=1)
    ni = size(my_saved_data(i)%sol,dim=2)
    !
    write_probe_data: &
    do j = 1,ni
      write (io,901) my_saved_data(i)%iter(j), &
                     my_saved_data(i)%time(j), &
                    (my_saved_data(i)%sol(n,j), n=1,np)
    end do write_probe_data
    !
    flush (io)
    !
    ! Reset the saved data for this probe
    !
    my_saved_data(i)%iter(:) = 0
    my_saved_data(i)%time(:) = zero
    my_saved_data(i)%sol(:,:) = zero
    !
  end do loop_over_probes
  !
  ! Reset the current save index to 0
  !
  current_save_idx = 0
  !
  ! Format Statements
  !
  ! NOTE: concern with low-precision output of time when small time steps
 !901 format(i10,es26.16e3,99es26.16e3)  ! real(dp) needs es24.16e3
  901 format(i10,es26.16e3,99es22.12e2)  ! custom
  !
end subroutine write_saved_probe_data
!
!###############################################################################
!> Open output file for a local probe number and write header.
!###############################################################################
subroutine probe_file_open(i,io)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr,n_solpts,xyz,cell_map
  use ovar, only : grid_scaling_factor
  !
  !.. Formal Arguments ..
  integer, intent(in)  :: i   !< local probe index
  integer, intent(out) :: io  !< io unit number
  !
  !.. Local Variables ..
  ! max number of output variables
  !
  character (len=:), allocatable :: file_prefix
  character (len=:), allocatable :: zone_prefix
  character (len=:), allocatable :: filename
  character (len=:), allocatable :: zonename
  character (len=30) :: varnames(nvarmax)
  !
  character (len=100) :: tmptxt
  !
  character(len=:), allocatable :: text_rakenum
  character(len=:), allocatable :: text_probenum
  !
  logical :: exists
  integer :: ipt, n, nvar
  !
  real(wp), allocatable :: prb_xyz(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "probe_file_open"
  character(len=*), parameter :: char_xyz = "xyz"
  !
continue
  !
  ! Set actual number of output variables.
  nvar = nr + 4
  !
  ! Set probe output variable names.
  varnames(1:nvarmax) = [ 'Iteration              ', &
                          'TimeValues             ', &
                          'Density                ', &
                          'MomentumX              ', &
                          'MomentumY              ', &
                          'MomentumZ              ', &
                          'EnergyStagnationDensity'  ]
  ! Adjust energy output position for 1d, 2d cases.
  varnames(nvar) = 'EnergyStagnationDensity'
  !
  ! Nearest solution point index
  ipt = my_probes(i)%solpt_idx
  !
  if (my_probes(i)%probe_info%probe_type == 3) then
    !
    ! F2003 AUTO-REALLOC
    text_probenum = get_cell_idx_string(my_probes(i)%probe_info%probe_idx, &
                                        my_probes(i)%probe_info%sub_idx)
    !
    file_prefix = "cell" ! F2003 AUTO-REALLOC
    zone_prefix = "Cell" ! F2003 AUTO-REALLOC
    !
  else if (my_probes(i)%probe_info%probe_type == 2) then
    !
    ! F2003 AUTO-REALLOC
    text_probenum = get_rake_idx_string(my_probes(i)%probe_info%probe_idx, &
                                        my_probes(i)%probe_info%sub_idx)
    !
    file_prefix = "rake" ! F2003 AUTO-REALLOC
    zone_prefix = "Rake" ! F2003 AUTO-REALLOC
    !
   !write (tmptxt,101) my_probes(i)%probe_info%rake_idx
   !text_rakenum = trim(adjustl(tmptxt)) ! F2003 AUTO-REALLOC
   !!
   !write (tmptxt,101) my_probes(i)%probe_info%probe_idx
   !text_probenum = trim(adjustl(tmptxt)) ! F2003 AUTO-REALLOC
   !!
   !filename = probes_dir // 'rake-probe.' // & ! F2003 AUTO-REALLOC
   !           text_rakenum // '-' // text_probenum // '.tec'
   !zonename = 'Rake-Probe ' // &               ! F2003 AUTO-REALLOC
   !           text_rakenum // '-' // text_probenum
    !
  else
    !
    ! F2003 AUTO-REALLOC
    text_probenum = get_probe_idx_string(my_probes(i)%probe_info%probe_idx)
    !
    file_prefix = "point" ! F2003 AUTO-REALLOC
    zone_prefix = "Point" ! F2003 AUTO-REALLOC
    !
   !write (tmptxt,101) my_probes(i)%probe_info%probe_idx
   !text_probenum = trim(adjustl(tmptxt)) ! F2003 AUTO-REALLOC
   !!
   !filename = probes_dir // 'point-probe.' // & ! F2003 AUTO-REALLOC
   !           text_probenum // '.tec'
   !zonename = 'Point-Probe ' // text_probenum ! F2003 AUTO-REALLOC
    !
  end if
  !
  filename = probes_dir // file_prefix // '-probe.' // & ! F2003 AUTO-REALLOC
             text_probenum // '.tec'
  zonename = zone_prefix // '-Probe ' // text_probenum ! F2003 AUTO-REALLOC
  !
  ! Test if file exists. Used to determine which header to write.
  inquire (file=trim(filename),exist=exists)
  !
  open (newunit=io,file=trim(filename),form='formatted',position='append')
  !
  ! Write file header comment
  write (io,900)
  !
  ! Write variable header
  if ( exists ) then
    write (io,901) '# ',(trim(varnames(n)),n=1,nvar)
  else
    write (io,901)   '',(trim(varnames(n)),n=1,nvar)
  end if
  !
  ! Write zone header
  write (io,902) trim(zonename)
  write (io,903)
  !
  ! Write AuxData
  !
  prb_xyz = my_probes(i)%probe_info%xyz(1:nr) ! F2003 AUTO-REALLOC
  !
  ! Probe xyz values in input grid units.
  !
  do n = 1, nr
    write (io,904) 'probe_' // char_xyz(n:n) // '_input',prb_xyz(n)
  end do
  !
  ! Probe xyz values in meters.
  !
  prb_xyz(:) = prb_xyz(:) * grid_scaling_factor
  do n = 1, nr
    write (io,904) 'probe_' // char_xyz(n:n) // '_meters',prb_xyz(n)
  end do
  !
  ! Nearest solution point xyz values in meters.
  !
  do n = 1, nr
    write (io,904) 'found_' // char_xyz(n:n) // '_meters',xyz(n,ipt)
  end do
  !
  ! Flush this io buffer before returning
  !
  flush (io)
  !
  ! Format Statements
  !
  101 format(i0)
  900 format('# All quantities in SI units with meters and seconds.')
  901 format(a,'Variables=',99('"',a,'"',:,','))
  902 format('Zone T="',a,'"')
  903 format('F=POINT')
! 904 format('AuxData ',a,' = "',es13.6e2,'"')  ! real(sp) needs es13.6e2
  904 format('AuxData ',a,' = "',es24.16e3,'"')  ! real(dp) needs es24.16e3
  !
end subroutine probe_file_open
!
!###############################################################################
!
pure &
function get_cell_idx_string(cell_idx,probe_idx) result(return_value)
  !
  !.. Use Statements ..
  use geomvar_mod, only : cell,cell_map
  !
  !.. Formal Arguments ..
  integer, intent(in) :: cell_idx
  integer, intent(in) :: probe_idx
  !
  !.. Function Result ..
  character(len=:), allocatable :: return_value
  !
  !.. Local Scalars ..
  integer :: npts,glb_cell
  character(len=:), allocatable :: cell_strng
  character(len=:), allocatable :: probe_strng
  !
continue
  !
  glb_cell = cell_map(mypnum+1)%loc_to_glb(cell_idx)
  !
  ! Get the string for the cell index
  !
  ! F2003 AUTO-REALLOC
  cell_strng = get_zero_padded_int_string(glb_cell,glb_cell)
  !
  ! Get the string for the probe index
  !
  npts = cell(cell_idx)%end_sp - cell(cell_idx)%beg_sp + 1
  probe_strng = get_zero_padded_int_string(probe_idx,npts)
  !
  ! Combine the cell and probe strings to get the final string
  !
  return_value = cell_strng // '-' // probe_strng ! F2003 AUTO-REALLOC
  !
end function get_cell_idx_string
!
!###############################################################################
!
pure &
function get_rake_idx_string(rake_idx,probe_idx) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : rake_probes
  !
  !.. Formal Arguments ..
  integer, intent(in) :: rake_idx
  integer, intent(in) :: probe_idx
  !
  !.. Function Result ..
  character(len=:), allocatable :: return_value
  !
  !.. Local Scalars ..
  character(len=:), allocatable :: rake_strng
  character(len=:), allocatable :: probe_strng
  !
continue
  !
  ! Get the string for the rake index
  !
  ! F2003 AUTO-REALLOC
  rake_strng = get_zero_padded_int_string(rake_idx,size(rake_probes))
  !
  ! Get the string for the probe index
  !
  probe_strng = get_zero_padded_int_string(probe_idx, &
                                           rake_probes(rake_idx)%num_probes)
  !
  ! Combine the rake and probe strings to get the final string
  !
  return_value = rake_strng // '-' // probe_strng ! F2003 AUTO-REALLOC
  !
end function get_rake_idx_string
!
!###############################################################################
!
pure &
function get_probe_idx_string(probe_idx) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: probe_idx
  !
  !.. Function Result ..
  character(len=:), allocatable :: return_value
  !
continue
  !
  ! F2003 AUTO-REALLOC
  return_value = get_zero_padded_int_string(probe_idx,num_point_probes)
  !
end function get_probe_idx_string
!
!###############################################################################
!
pure &
function get_zero_padded_int_string(val,max_val) result(return_value)
  !
  !.. Formal Arguments ..
  integer, intent(in) :: val
  integer, intent(in) :: max_val
  !
  !.. Function Result ..
  character(len=:), allocatable :: return_value
  !
  !.. Local Scalars ..
  integer :: n
  !
  character(len=10) :: cnum
  character(len=20) :: frmt
  character(len=50) :: strng
  !
continue
  !
  write (cnum,1) max_val
  !
  n = len_trim(adjustl(cnum))
  !
  write (frmt,2) n,n
  !
  write (strng,frmt) val
  !
  return_value = trim(adjustl(strng)) ! F2003 AUTO-REALLOC
  !
  ! Format Statements
  !
  1 format (i0)
  2 format ('(i',i0,'.',i0,')')
  !
end function get_zero_padded_int_string
!
!###############################################################################
!> Close output file for probe number i.
!###############################################################################
subroutine probe_file_close(i)
  !
  integer, intent(in) :: i  !< global probe index
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "probe_file_close"
  !
continue
  !
  close (my_probes(i)%io_unit)
  my_probes(i)%io_unit = 0
  !
end subroutine probe_file_close
!
!###############################################################################
!
end module probes_mod
