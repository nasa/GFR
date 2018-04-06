module channel_mod
  !
  ! Module containing procedures for the channel flow problem
  !
  !.. Use Statements ..
  use module_kind_types
  use eqn_idx, only : nec,nmx,nee,nq
  !
  implicit none
  !
  private
  !
  !.. Public procedures
  public :: read_channel_init_file
  public :: channel_initialization
  public :: channel_source
  public :: spatial_ave_channel_flow
  !
  !.. Local Module Variables ..
  real(wp), allocatable, save :: init_sol(:,:)
  !
contains
!
!###############################################################################
!
subroutine read_channel_init_file()
  !
  !.. Use Statements ..
  use ovar, only : channel_init_file
  !
  !.. Local Parameters..
  character(len=*), parameter :: pname = "read_channel_init_file"
  !
continue
  !
end subroutine read_channel_init_file
!
!###############################################################################
!
pure function channel_initialization( xyz ) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : bc_in
  !
  !.. Formal Arguments ..
  real(wp), dimension(:,:), intent(in) :: xyz
  !
  !.. Function Result ..
  real(wp), dimension(1:size(xyz,dim=2),1:nq) :: return_value
  !
  !.. Local Scalars ..
  integer :: k,np
  !
continue
  !
  np = size(xyz,dim=2)
  !
  if (allocated(init_sol)) then
    !
    !
  else
    !
    do k = 1,np
      return_value(k,:) = bc_in(0)%cv(:)
    end do
    !
  end if
  !
end function channel_initialization
!
!###############################################################################
!
subroutine channel_source()
  !
  !.. Use Statements ..
  use geovar,      only : n_solpts
  use ovar,        only : channel_body_force
  use ovar,        only : rhoref,aref
  use flowvar,     only : residual,usp
  use metrics_mod, only : jac
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: bfx,hc,vx
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "channel_source"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  bfx = channel_body_force / (rhoref*aref*aref)
  !
  do n = 1,n_solpts
    !
    vx = usp(n,nmx) / usp(n,nec)
    !
    hc = -bfx * jac(n)
    !
    residual(n,nmx) = residual(n,nmx) + hc
    residual(n,nee) = residual(n,nee) + hc * vx
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine channel_source
!
!###############################################################################
!
subroutine spatial_ave_channel_flow()
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts,ncell,cell,xyz
  use flowvar, only : uavesp,time_ave_variables
  use ovar,    only : rhoref,aref,tref,tkeref,tlsref
  !
  !.. Local Scalars ..
  integer  :: l,m,n,nc,n1,n2,ierr,io
  integer  :: nvar,n_yloc
  !
  character(len=1)  :: ta_char
  character(len=1)  :: vmult
  character(len=50) :: varname
  !
  !.. Local Arrays ..
  real(wp) :: y_vals(1:10000)
  real(wp) :: var_dim(1:size(uavesp,dim=2))
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: y_ptr(:)
  integer,  allocatable :: y_counts(:)
  real(wp), allocatable :: y_data(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "spatial_ave_channel_flow"
  character(len=*), parameter :: avefile = "spatial_ave.dat"
  !
continue
  !
  nvar = size(uavesp,dim=2)
  !
  var_dim(1:nvar) = one
  !
  allocate ( y_ptr(1:n_solpts) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"y_ptr",1,__LINE__,__FILE__,ierr,error_message)
  !
  !
  !
  y_vals(:) = zero
  !
  n_yloc = 1
  y_vals(1) = xyz(2,1)
  y_ptr(1) = 1
  !
  solpts_loop: do n = 2,n_solpts
    !
    do l = 1,n_yloc
      if (abs(y_vals(l)-xyz(2,n)) <= eps6) then
        y_ptr(n) = l
        cycle solpts_loop
      end if
    end do
    !
    n_yloc = n_yloc + 1
    y_vals(n_yloc) = xyz(2,n)
    y_ptr(n) = n_yloc
    !
  end do solpts_loop
  !
  !
  !
  allocate ( y_data(1:nvar,1:n_yloc) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"y_data",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( y_counts(1:n_yloc) , source=0 , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"y_counts",1,__LINE__,__FILE__,ierr,error_message)
  !
  !
  !
  do n = 1,n_solpts
    !
    l = y_ptr(n)
    !
    y_counts(l) = y_counts(l) + 1
    !
    do m = 1,nvar
      y_data(m,l) = y_data(m,l) + uavesp(n,m)
    end do
    !
  end do
  !
  write (iout,2) minval(y_counts),maxval(y_counts)
  !
  !
  !
  do l = 1,n_yloc
    do m = 1,nvar
      y_data(m,l) = y_data(m,l) / real(y_counts(l),kind=wp)
    end do
  end do
  !
  !
  open (newunit=io,file=avefile,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,avefile,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (io,10)
  !
  do m = 1,nvar
    !
    vmult = ""
    varname = ""
    do l = 1,len_trim(time_ave_variables(m))
      ta_char = uppercase( time_ave_variables(m)(l:l) )
      select case (ta_char)
        case ("D")
          varname = trim(adjustl(varname)) // trim(vmult) // "Rho"
          var_dim(m) = var_dim(m) * rhoref
        case ("X")
          varname = trim(adjustl(varname)) // trim(vmult) // "Vx"
          var_dim(m) = var_dim(m) * aref
        case ("Y")
          varname = trim(adjustl(varname)) // trim(vmult) // "Vy"
          var_dim(m) = var_dim(m) * aref
        case ("Z")
          varname = trim(adjustl(varname)) // trim(vmult) // "Vz"
          var_dim(m) = var_dim(m) * aref
        case ("E")
          varname = trim(adjustl(varname)) // trim(vmult) // "E"
          var_dim(m) = var_dim(m) * rhoref * aref * aref
        case ("P")
          varname = trim(adjustl(varname)) // trim(vmult) // "P"
          var_dim(m) = var_dim(m) * rhoref * aref * aref
        case ("T")
          varname = trim(adjustl(varname)) // trim(vmult) // "T"
          var_dim(m) = var_dim(m) * tref
        case ("R")
          varname = trim(adjustl(varname)) // trim(vmult) // "Vr"
          var_dim(m) = var_dim(m) * aref
        case ("F")
          varname = trim(adjustl(varname)) // trim(vmult) // "Vf"
          var_dim(m) = var_dim(m) * aref
        case ("K")
          varname = trim(adjustl(varname)) // trim(vmult) // "TKE"
          var_dim(m) = var_dim(m) * tkeref
        case ("W")
          varname = trim(adjustl(varname)) // trim(vmult) // "TLS"
          var_dim(m) = var_dim(m) * tlsref
      end select
      vmult = "*"
    end do
    !
    if (m == nvar) then
      write (io,11) trim(adjustl(varname))
    else
      write (io,11) trim(adjustl(varname)),","
    end if
    !
  end do
  !
  write (io,12) n_yloc
  !
  do n = 1,n_yloc
    write (io,13) y_vals(n),(y_data(m,n)*var_dim(m), m=1,nvar)
  end do
  !
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,avefile,2,__LINE__,__FILE__,ierr,error_message)
  !
  !
  write (error_message,1)
  call stop_recon(stop_mpi,pname,__LINE__,__FILE__,error_message)
  !
  1 format ("After spatial averaging a channel flow restart file.")
  2 format (/,"  Points used for averaging: Min = ",i0,/, &
              "                             Max = ",i0,/)
  !
  10 format ('VARIABLES = "Y",')
  11 format ('"',a,'"',:,a)
  12 format ('ZONE T="Spatial Average", I=',i0,', DATAPACKING=POINT')
  13 format (100es14.6)
  !
end subroutine spatial_ave_channel_flow
!
!###############################################################################
!
end module channel_mod
!
!###############################################################################
!
