!###############################################################################
!> Module for point probes output.
!>
!> Usage:
!> - Call local_probes_write to write data.
!>   Will allocate/initialize local probes and open necessary output files.
!> - Call local_probes_finalize to close probe files and deallocate memory.
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
  public  :: local_probes_write
  public  :: local_probes_finalize
  !
  ! .. Private Procedures ..
  private :: local_probes_initialize
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
  type :: probe_info_t
    ! rake_idx: index of the rake from the input file
    integer :: rake_idx = 0
    ! probe_idx:
    !     If rake_idx >  0, this is the index of the probe within the rake
    !     If rake_idx <= 0, this is the point_probe index from the input file
    integer :: probe_idx = 0
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
  end type my_probes_t
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
  ! LOCAL_PROBE : Derived type containing probe info for local partition
  !
  type :: local_probe_t
    integer :: global_point_probe_idx = 0
    integer :: global_rake_probe_idx = 0
    integer :: solpt_idx = 0          ! sol pt on local partition
    integer :: io_unit = 0
  end type local_probe_t
  !
  ! MY_PROBES : Derived type array containing info for
  !             all probes on the local partition
  !
  type(my_probes_t), allocatable, save :: my_probes(:)
  !
  ! MY_SAVED_DATA : Derived type array containing the saved data
  !                 collections for all probes on the local partition
  !
  type(my_saved_data_t), allocatable, save :: my_saved_data(:)
  !
  type(local_probe_t), allocatable, save :: local_probes(:)
  !
  character(len=:), allocatable, save :: probes_dir
  !
  logical(lk), save :: local_probes_initialized = fals
  !
  integer, save :: num_total_probes = 0
  integer, save :: num_point_probes = 0
  integer, save :: num_rake_probes  = 0
  !
  integer, save :: current_save_idx = 0
  !
contains
!
!###############################################################################
!> Initialize memory and file io for local point probes.
!###############################################################################
subroutine initialize_my_probes()
  !
  use geomvar_mod, only : nr, n_solpts, xyz
  use ovar,        only : grid_scaling_factor, point_probes, rake_probes
  use ovar,        only : probes_collection_interval, probes_write_interval
  use io_mod,      only : create_results_dir
  !
  !.. Local Variables ..
  integer  :: io, i, j, n, np
  integer  :: min_idx
  real(wp) :: ds, dist2, min_dist2
  !
  integer,  allocatable :: local_probes_solpt_idx(:)  !   (num_total_probes)
  real(wp), allocatable :: pairs(:,:)                 ! (2,num_total_probes)
  !
  real(wp), allocatable :: prb_xyz(:)
  real(wp), allocatable :: xyz_beg(:)
  real(wp), allocatable :: xyz_end(:)
  !
  type(probe_info_t), allocatable :: probe_info(:)
  !
  character(len=1000) :: temp_dir
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "initialize_my_probes"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  num_point_probes = 0
  num_rake_probes = 0
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
  num_total_probes = num_point_probes + num_rake_probes
  !
  ! If there are no probes, set the probe collection interval to zero so
  ! it will stop trying to collect/write probe data and exit this subroutine
  !
  if (num_total_probes < 1) then
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
  ! Create an array to hold the coordinates for the all point and rake probes
  ! and their respective indices pointing back to where they came from
  !
  allocate ( probe_info(1:num_total_probes) )
  !
  ! Add any point probes to probe_info
  !
  do n = 1,num_point_probes
    probe_info(n)%rake_idx = 0 ! not a rake
    probe_info(n)%probe_idx = n
    probe_info(n)%xyz(1:nr) = point_probes(n)%xyz(1:nr)
  end do
  !
  ! Add the info for each individual probe from each rake to probe_info
  !
  n = num_point_probes
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
        probe_info(n)%rake_idx = i
        probe_info(n)%probe_idx = j
        probe_info(n)%xyz(1:nr) = xyz_beg(1:nr) + &
                                  ds*(xyz_end(1:nr)-xyz_beg(1:nr))
      end do
      !
    else if (np == 1) then
      !
      n = n + 1
      !
      ! If only one probe is requested, average the beginning and end
      ! coordinates for the rake to get the probe coordinates.
      !
      probe_info(n)%rake_idx = i
      probe_info(n)%probe_idx = 1
      probe_info(n)%xyz(1:nr) = half*(xyz_beg(1:nr) + xyz_end(1:nr))
      !
   !NOTE: skip if np <= 0
    end if
    !
  end do loop_over_rakes
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Allocate the arrays needed to find the closest solution point for each
  ! probe and which partition that solution point belongs to.
  !
  local_probes_solpt_idx = init_array(num_total_probes,0) ! F2003 AUTO-REALLOC
  pairs = init_array(2,num_total_probes,-99.0_wp)         ! F2003 AUTO-REALLOC
  !
  ! Find local closest sol pt to each probe.
  !
  do i = 1,num_total_probes
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
    local_probes_solpt_idx(i) = min_idx  ! min solpt index
    pairs(1,i) = min_dist2               ! min dist^2
    pairs(2,i) = real(mypnum,kind=wp)    ! mpi process
  end do
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Find global closest sol pt to each probe and its mpi process
  !
  call mpi_allreduce(MPI_IN_PLACE,pairs,size(pairs,dim=2),mpi_2flttyp, &
                     MPI_MINLOC,MPI_COMM_WORLD,mpierr)
  !
  !#############################################################################
  !#############################################################################
  !#############################################################################
  !
  ! Create the array my_probes to hold all the info for point and rake probes
  ! that are located on this partition
  !
  n = count( nint(pairs(2,:)) == mypnum )
  !
  allocate( my_probes(1:n) )
  !
  j = 0  ! local point_probes index
  do i = 1,num_total_probes
    if ( nint(pairs(2,i)) == mypnum ) then
      j = j + 1
      my_probes(j)%probe_info = probe_info(i)
      my_probes(j)%solpt_idx = local_probes_solpt_idx(i)
      call probe_file_open(j,io)
      my_probes(j)%io_unit = io
    end if
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
  ! Mark this module as being initialized
  !
  local_probes_initialized = true
  !
  deallocate ( local_probes_solpt_idx )
  deallocate ( pairs )
  deallocate ( prb_xyz )
  deallocate ( xyz_beg )
  deallocate ( xyz_end )
  deallocate ( probe_info )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine initialize_my_probes
!
!###############################################################################
!> Initialize memory and file io for local point probes.
!###############################################################################
subroutine local_probes_initialize()
  !
  use geomvar_mod, only : nr, n_solpts, xyz
  use ovar,        only : grid_scaling_factor, point_probes
  use io_mod,      only : create_results_dir
  !
  !.. Local Variables ..
  integer  :: io, i, j, n
  integer  :: min_idx
  real(wp) :: dist2, min_dist2
  !
  integer,  allocatable :: local_probes_solpt_idx(:)  !   (size(point_probes))
  real(wp), allocatable :: pairs(:,:)                 ! (2,size(point_probes))
  !
  character(len=1000) :: temp_dir
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "local_probes_initialize"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  probes_dir = "" ! F2003 AUTO-REALLOC
  if (size(point_probes) > 0) then
    call create_results_dir("Probe_files",temp_dir,skip_lustre=true)
    probes_dir = trim(adjustl(temp_dir)) ! F2003 AUTO-REALLOC
  end if
  !
  allocate( local_probes_solpt_idx(size(point_probes)) )
  allocate( pairs(2,size(point_probes)) )
  !
  ! Find local closest sol pt to each probe.
  !
  local_probes_solpt_idx(:) = 0
  do i = 1, size(point_probes)
    min_idx = 0               ! sol pt index with min_dist2
    min_dist2 = huge(1.0_wp)  ! min distance^2 from probe to sol pt
    do n = 1, n_solpts
      ! Convert probe location from input grid units to meters.
      dist2 = sum( ( point_probes(i)%xyz(1:nr) * grid_scaling_factor &
                   - xyz(1:nr,n) )**2 )
      if ( dist2 < min_dist2 ) then
        min_dist2 = dist2
        min_idx = n
      end if
    end do
    local_probes_solpt_idx(i) = min_idx  ! min solpt index
    pairs(1,i) = min_dist2               ! min dist^2
    pairs(2,i) = real(mypnum,wp)         ! mpi process
  end do
  !
  ! Find global closest sol pt to each probe and its mpi process
  !
  call mpi_allreduce(MPI_IN_PLACE,pairs,size(pairs,dim=2),mpi_2flttyp, &
                     MPI_MINLOC,MPI_COMM_WORLD,mpierr)
  !
  ! Fill the local point_probes data for this process and open output files
  !
  n = count( int(pairs(2,:)) == mypnum )
  allocate( local_probes(n) )
  j = 0  ! local point_probes index
  do i = 1, size(point_probes)  ! global point_probes index
    if ( int(pairs(2,i)) == mypnum ) then
      j = j + 1
      local_probes(j)%global_point_probe_idx = i
      local_probes(j)%solpt_idx = local_probes_solpt_idx(i)
      call probe_file_open_og(j,io)
      local_probes(j)%io_unit = io
    end if
  end do
  !
  local_probes_initialized = .true.
  !
  deallocate( local_probes_solpt_idx )
  deallocate( pairs )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine local_probes_initialize
!
!###############################################################################
!> Finalize memory and file io for local point probes.
!###############################################################################
subroutine local_probes_finalize()
  !
  !.. Local Variables ..
  integer :: i
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "local_probes_finalize"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (current_save_idx > 0) then
    call write_saved_probe_data
  end if
  !
  do i = 1, size(my_probes)
    call probe_file_close(i)
  end do
  !
  deallocate( my_probes )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine local_probes_finalize
!
!###############################################################################
!> Collect data for all local probes.
!>
!> Will call local_probes_initialize automatically when necessary,
!> which will open the probe output files.
!> If this interval of collecting data is equal to probes_write_interval,
!> write out the saved probe data.
!> Must call local_probes_finalize to close all probe output files.
!###############################################################################
subroutine local_probes_write(iter,time,usp)
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
  character(len=*), parameter :: pname = "local_probes_write"
  !
continue
  !
  if ( .not. local_probes_initialized ) then
    call initialize_my_probes()
  endif
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
end subroutine local_probes_write
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
  use geomvar_mod, only : nr, n_solpts, xyz
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
  if (my_probes(i)%probe_info%rake_idx > 0) then
    !
    ! F2003 AUTO-REALLOC
    text_probenum = get_rake_idx_string(my_probes(i)%probe_info%rake_idx, &
                                        my_probes(i)%probe_info%probe_idx)
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
!> Open output file for a local probe number and write header.
!###############################################################################
subroutine probe_file_open_og(i,io)
  !
  !.. Use Statements ..
  use geomvar_mod, only : nr, n_solpts, xyz
  use ovar, only : grid_scaling_factor, point_probes
  !
  !.. Formal Arguments ..
  integer, intent(in)  :: i   !< local probe index
  integer, intent(out) :: io  !< io unit number
  !
  !.. Local Variables ..
  ! max number of output variables
  integer, parameter :: nvarmax = 7
  !
  character (len=:), allocatable :: filename
  character (len=:), allocatable :: zonename
  character (len=5)  :: textnum
  character (len=1)  :: char_xyz(3) = ['x','y','z']
  character (len=30) :: varnames(nvarmax)
  !
  logical :: exists
  integer :: ig, ipt, n, nvar
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "probe_file_open_og"
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
  ipt = local_probes(i)%solpt_idx
  !
  n = local_probes(i)%global_point_probe_idx
  write(textnum,101) n
  filename = trim(probes_dir)//'probe.'//trim(textnum)//'.tec'
  zonename = 'Probe '//trim(textnum)
  !
  ! Test if file exists. Used to determine which header to write.
  inquire(file=trim(filename),exist=exists)
  !
  open(newunit=io,file=trim(filename),form='formatted',position='append')
  !
  ! Write file header comment
  write(io,900)
  !
  ! Write variable header
  if ( exists ) then
    write(io,901) '# ',(trim(varnames(n)),n=1,nvar)
  else
    write(io,901)   '',(trim(varnames(n)),n=1,nvar)
  end if
  !
  ! Write zone header
  write(io,902) trim(zonename)
  write(io,903)
  !
  ! Write AuxData
  ig = local_probes(i)%global_point_probe_idx
  ! Probe xyz values in input grid units.
  do n = 1, nr
    write(io,904) 'probe_' // char_xyz(n:n) // &
                  '_input',point_probes(ig)%xyz(n)
  end do
  ! Probe xyz values in meters.
  do n = 1, nr
    write(io,904) 'probe_' // char_xyz(n:n) // &
                  '_meters',point_probes(ig)%xyz(n) * grid_scaling_factor
  end do
  ! Nearest solution point xyz values in meters.
  do n = 1, nr
    write(io,904) 'found_' // char_xyz(n:n) // '_meters',xyz(n,ipt)
  end do
  !
  101 format(i0)
  900 format('# All quantities in SI units with meters and seconds.')
  901 format(a,'Variables=',99('"',a,'"',:,','))
  902 format('Zone T="',a,'"')
  903 format('F=POINT')
! 904 format('AuxData ',a,' = "',es13.6e2,'"')  ! real(sp) needs es13.6e2
  904 format('AuxData ',a,' = "',es24.16e3,'"')  ! real(dp) needs es24.16e3
  !
  flush(io)
  !
end subroutine probe_file_open_og
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
