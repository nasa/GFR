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
  use ovar, only : probes_interval, probes
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
  ! LOCAL_PROBE : Derived type containing probe info for local partition
  !
  type, public :: local_probe_t
    integer :: global_probe_idx = 0
    integer :: solpt_idx = 0          ! sol pt on local partition
    integer :: io_unit = 0
  end type local_probe_t
  !
  type(local_probe_t), allocatable :: local_probes(:)
  !
  character(len=:), allocatable, save :: probes_dir
  !
  logical :: local_probes_initialized = .false.
  !
contains
!
!###############################################################################
!> Initialize memory and file io for local point probes.
!###############################################################################
subroutine local_probes_initialize()
  !
  use geomvar_mod, only : nr, n_solpts, xyz
  use ovar,        only : grid_scaling_factor, probes
  use io_mod,      only : create_results_dir
  !
  !.. Local Variables ..
  integer  :: io, i, j, n
  integer  :: min_idx
  real(wp) :: dist2, min_dist2
  !
  integer,  allocatable :: local_probes_solpt_idx(:)  !   (size(probes))
  real(wp), allocatable :: pairs(:,:)                 ! (2,size(probes))
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
  if (size(probes) > 0) then
    call create_results_dir("Probe_files",temp_dir,skip_lustre=true)
    probes_dir = trim(adjustl(temp_dir)) ! F2003 AUTO-REALLOC
  end if
  !
  allocate( local_probes_solpt_idx(size(probes)) )
  allocate( pairs(2,size(probes)) )
  !
  ! Find local closest sol pt to each probe.
  !
  local_probes_solpt_idx(:) = 0
  do i = 1, size(probes)
    min_idx = 0               ! sol pt index with min_dist2
    min_dist2 = huge(1.0_wp)  ! min distance^2 from probe to sol pt
    do n = 1, n_solpts
      ! Convert probe location from input grid units to meters.
      dist2 = sum( ( probes(i)%xyz(1:nr) * grid_scaling_factor &
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
  ! Fill the local probes data for this process and open output files
  !
  n = count( int(pairs(2,:)) == mypnum )
  allocate( local_probes(n) )
  j = 0  ! local probes index
  do i = 1, size(probes)  ! global probes index
    if ( int(pairs(2,i)) == mypnum ) then
      j = j + 1
      local_probes(j)%global_probe_idx = i
      local_probes(j)%solpt_idx = local_probes_solpt_idx(i)
      call probe_file_open(j,io)
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
  do i = 1, size(local_probes)
    call probe_file_close(i)
  end do
  !
  deallocate( local_probes )
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine local_probes_finalize
!
!###############################################################################
!> Write data for all local probes.
!>
!> Will call local_probes_initialize automatically when necessary,
!> which will open the probe output files.
!> Must call local_probes_finalize to close all probe output files.
!###############################################################################
subroutine local_probes_write(iter,time,usp)
  !
  !.. Use Statements ..
  use eqn_idx_mod, only : nec,nmx,nmy,nmz,nmb,nme,nee
  use ovar,        only : time_ref,rhoref,aref
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: iter      !< solution iteration number
  real(wp), intent(in) :: time      !< solution time [nondim]
  real(wp), intent(in) :: usp(:,:)  !< solution conserved variables [nondim]
  !
  !.. Local Variables ..
  integer :: i, io, ipt, n
  real(wp) :: ref(5)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "local_probes_write"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if ( .not. local_probes_initialized ) then
    call local_probes_initialize()
  endif
  !
  ! non-dim reference values
  ref(nec)     = rhoref            ! density
  ref(nmb:nme) = rhoref * aref     ! momentum
  ref(nee)     = rhoref * aref**2  ! energy
  !
  do i = 1, size(local_probes)
    io  = local_probes(i)%io_unit
    ipt = local_probes(i)%solpt_idx
    !
    write(io,901) iter,time*time_ref,(usp(ipt,n)*ref(n),n=nec,nee)
    flush(io)
    ! concern with low-precision output of time when small time steps
!   901 format(i10,es26.16e3,99es15.6e2)   ! real(sp) needs es13.6e2
!   901 format(i10,es26.16e3,99es26.16e3)  ! real(dp) needs es24.16e3
    901 format(i10,es26.16e3,99es22.12e2)  ! custom
  enddo
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine local_probes_write
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
  integer, parameter :: nvarmax = 7
  !
  character (len=:), allocatable :: filename
  character (len=:), allocatable :: zonename
  character (len=5)  :: textnum
  character (len=1)  :: char_xyz(3) = ['x','y','z']
  character (len=2)  :: prefix
  character (len=30) :: varnames(nvarmax)
  !
  logical :: exists
  integer :: ig, ipt, n, nvar
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "probe_file_open"
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
  n = local_probes(i)%global_probe_idx
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
  ig = local_probes(i)%global_probe_idx
  ! Probe xyz values in input grid units.
  do n = 1, nr
    write(io,904) 'probe_' // char_xyz(n:n) // &
                  '_input',probes(ig)%xyz(n)
  end do
  ! Probe xyz values in meters.
  do n = 1, nr
    write(io,904) 'probe_' // char_xyz(n:n) // &
                  '_meters',probes(ig)%xyz(n) * grid_scaling_factor
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
end subroutine probe_file_open
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
  close(local_probes(i)%io_unit)
  local_probes(i)%io_unit = 0
  !
end subroutine probe_file_close
!
!###############################################################################
!
end module probes_mod
