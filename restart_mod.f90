module restart_mod
  !
  use module_kind_types
  !
  implicit none
  !
  private
  !
  public :: read_restart_file
  public :: write_restart_file
  !
contains
!
!###############################################################################
!
subroutine read_restart_file
  !
  !.. Use Statements ..
  use geovar,       only : ncell,pts_in_cell_ptr
  use geovar,       only : n_global_cell,global_nofc_ptr,global_pinc_ptr
  use ovar,         only : itrst,tmrst,dtrst
  use flowvar,      only : usp,global_usp
  use parallel_mod, only : cell_map
  !
  !.. Local Scalars ..
  integer :: n,nc,n1,n2,np,gp,ierr
  integer :: iorst
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "read_restart_file"
  character(len=*), parameter :: fname = "out.rst"
  !
continue
  !
  open(newunit=iorst,file=fname,status="old",action="read", &
                     form="unformatted",iostat=ierr)
  call io_error(pname,fname,1,ierr)
  !
  read (iorst) itrst, tmrst, dtrst
  read (iorst) global_usp
  !
  close (iorst,iostat=ierr)
  call io_error(pname,fname,2,ierr)
  !
  ! If using multiple processors, copy the solution data from 
  ! only those cells that are local to this processor into usp
  !
  if (ncpu > 1) then
    !
    do nc = 1,ncell
      !
      ! Indices of the datat for this cell in the usp array
      !
      n1 = pts_in_cell_ptr(nc)+1
      n2 = pts_in_cell_ptr(nc+1)
      np = n2-n1+1
      !
      ! Starting point of the data for this cell in the global_usp array
      !
      gp = global_pinc_ptr( cell_map(mypnum+1)%loc_to_glb(nc) )
      !
      ! Copy the solution data from the global_usp array for this cell
      !
      usp(:,n1:n2) = global_usp(:,gp+1:gp+np)
      !
    end do
    !
  else
    !
    ! For a serial job, just copy the global_usp array into usp
    !
    usp(:,:) = global_usp(:,:)
    !
  end if
  !
end subroutine read_restart_file
!
!###############################################################################
!
subroutine write_restart_file(itnum)
  !
  !.. Use Statements ..
  use ovar,    only : itcur,time,d_t
  use flowvar, only : global_usp
  use parallel_mod, only : collect_global_solution
  !
  !.. Formal Arguments ..
  integer, intent(in) :: itnum
  !
  !.. Local Scalars ..
  integer :: ierr,iorst
  character(len=40) :: fname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_restart_file"
  !
continue
  !
  ! Collect the global solution
  !
  call collect_global_solution
  !
  ! Only use the root partition to output the restart file
  !
  if (mypnum == 0) then
    !
    if (itnum > 0) then
      write (fname,'("out.",i0.8,".rst")') itcur
    else
      write (fname,'("out.",i0.8,".done.rst")') itcur
    end if
    !
    open (newunit=iorst,file=trim(adjustl(fname)),status="replace", &
                        action="write",form='unformatted',iostat=ierr)
    call io_error(pname,trim(adjustl(fname)),1,ierr)
    !
    write (iorst) itcur, time, d_t
    write (iorst) global_usp
    !
    close (iorst,iostat=ierr)
    call io_error(pname,trim(adjustl(fname)),2,ierr)
    !
  end if
  !
  call mpi_barrier(mpi_comm_world,mpierr)
  !
end subroutine write_restart_file
!
!###############################################################################
!
end module restart_mod
