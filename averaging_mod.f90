module averaging_mod
  !
  use module_kind_types
  !
  implicit none
  !
  private
  !
  public :: average_restart_files
  !
  character(len=100), save, allocatable :: restart_file_list(:)
  real(wp), save, allocatable :: restart_times(:)
  !
contains
!
!###############################################################################
!
subroutine average_restart_files()
  !
  !.. Use Statements ..
  use ovar,        only : restart_list_file,restart_file
  use ovar,        only : output_time_averaging,read_time_ave_restart
  use ovar,        only : ave_start_time,prev_ave_time
  use ovar,        only : time,itcur,tmrst,itrst
  use flowvar,     only : uavesp
  use time_mod,    only : initialize_time_ave,update_time_ave
  use restart_mod, only : read_restart_file,write_restart_file
  use io_mod,      only : write_parallel_cgns
  !
  !.. Local Scalars ..
  integer :: ierr,io,n
  !
  character(len=150) :: fname
  character(len=150) :: text
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "average_restart_files"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! First deallocate unnecessary arrays that were allocated for the solver but
  ! wont be needed since execution will abort at the end of this subroutine
  !
  call deallocate_solver_arrays
  !
  ! Read in the list of restart files that will be averaged
  !
  fname = trim(adjustl(restart_list_file))
  !
  open (newunit=io,file=fname,status="old",action="read", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  n = 0
  do
    read (io,1,iostat=ierr,iomsg=error_message) text
    if (ierr < 0) then
      exit
    else if (ierr > 0) then
      call io_error(pname,fname,-1,__LINE__,__FILE__,ierr,error_message)
    else
      n = n + 1
    end if
  end do
  !
  rewind (io)
  !
  allocate ( restart_times(1:n) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"restart_times",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( restart_file_list(1:n) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"restart_file_list",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  do n = 1,size(restart_file_list)
    read (io,1,iostat=ierr,iomsg=error_message) restart_file_list(n)
    call io_error(pname,fname,-1,__LINE__,__FILE__,ierr,error_message)
  end do
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
  !
  read_time_ave_restart = fals
  !
  restart_file = trim(adjustl(restart_file_list(1)))
  !
  if (mypnum == glb_root) then
    write (iout,2) trim(adjustl(restart_file))
  end if
  !
  call read_restart_file
  !
  time = tmrst
  itcur = itrst
  !
  restart_times(1) = time
  ave_start_time = time
  prev_ave_time = time
  !
  if (mypnum == glb_root) write (iout,3)
  !
  call initialize_time_ave
  !
  ! Initialize uavesp to zero
  !
  uavesp(:,:) = zero
  !
  do n = 2,size(restart_file_list)
    !
    restart_file = trim(adjustl(restart_file_list(n)))
    !
    if (mypnum == glb_root) then
      write (iout,4) trim(adjustl(restart_file))
    end if
    !
    call read_restart_file
    !
    time = tmrst
    itcur = itrst
    !
    restart_times(n) = time
    !
    if (mypnum == glb_root) write (iout,5)
    call update_time_ave
    !
  end do
  !
  output_time_averaging = true
  !
  if (mypnum == glb_root) write (iout,7) pname
  call write_restart_file(itcur)
  !
  if (mypnum == glb_root) write (iout,8) pname
  call write_parallel_cgns(open_CGNS)
  !
  call debug_timer(leaving_procedure,pname)
  !
  write (error_message,6)
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  !
  ! Format Statements
  !
  1 format (a)
  2 format (/,4x,"Reading initial restart file '",a,"'")
  3 format (/,4x,"Initializing the time-averaged variables")
  4 format (/,4x,"Working on restart file '",a,"'",/, &
              8x,"Reading the restart file")
  5 format (8x,"Updating the time-averaged variables")
  6 format ("STOPPING EXECUTION AFTER TIME-AVERAGING ", &
            "THE LIST OF RESTART FILES!")
  7 format (/,4x,"Calling write_restart_file from ",a)
  8 format (/,4x,"Calling write_parallel_cgns from ",a)
  !
end subroutine average_restart_files
!
!###############################################################################
!
subroutine deallocate_solver_arrays()
  !
  !.. Use Statements ..
  use geovar,      only : face,grid
  use flowvar,     only : faceusp,facedusp,faceflx,facexyz
  use flowvar,     only : residual,dtsp,uoldsp,usp_sum,rssprk
  use flowvar,     only : nodeusp,nodedusp,edgeusp,edgedusp,edgexyz
  use flowvar,     only : global_usp,global_var
  use metrics_mod, only : metrics_dt,geofa
  !
  !.. Local Scalars ..
  integer :: ierr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "deallocate_solver_arrays"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (allocated(face)) then
    deallocate ( face , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"face",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(grid)) then
    deallocate ( grid , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"grid",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(faceusp)) then
    deallocate ( faceusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"faceusp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(facedusp)) then
    deallocate ( facedusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facedusp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(faceflx)) then
    deallocate ( faceflx , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"faceflx",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(facexyz)) then
    deallocate ( facexyz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"facexyz",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(residual)) then
    deallocate ( residual , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"residual",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(dtsp)) then
    deallocate ( dtsp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"dtsp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(uoldsp)) then
    deallocate ( uoldsp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"uoldsp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(usp_sum)) then
    deallocate ( usp_sum , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"usp_sum",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(rssprk)) then
    deallocate ( rssprk , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"rssprk",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(nodeusp)) then
    deallocate ( nodeusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodeusp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(nodedusp)) then
    deallocate ( nodedusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"nodedusp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edgeusp)) then
    deallocate ( edgeusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edgeusp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edgedusp)) then
    deallocate ( edgedusp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edgedusp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(edgexyz)) then
    deallocate ( edgexyz , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"edgexyz",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(global_usp)) then
    deallocate ( global_usp , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_usp",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(global_var)) then
    deallocate ( global_var , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"global_var",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(metrics_dt)) then
    deallocate ( metrics_dt , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"metrics_dt",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(geofa)) then
    deallocate ( geofa , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"geofa",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine deallocate_solver_arrays
!
!###############################################################################
!
end module averaging_mod
