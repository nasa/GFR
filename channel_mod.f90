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
  public :: spatial_ave_z_planes
  public :: spatial_ave_x_and_z_planes
  !
  !.. Local Module Variables ..
  real(wp), allocatable, save :: init_sol(:,:)
  !
  !.. Local Module Parameters ..
  real(wp), parameter :: match_tolerance = eps9
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
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "spatial_ave_channel_flow"
  !
continue
  !
  call spatial_ave_z_planes(called_from_channel_mod=true)
  !
  call spatial_ave_x_and_z_planes(called_from_channel_mod=true)
  !
end subroutine spatial_ave_channel_flow
!
!###############################################################################
!
subroutine spatial_ave_x_and_z_planes(called_from_channel_mod)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts,xyz
  use flowvar, only : uavesp
  use ovar,    only : itcur,time,time_ref,time_scaling_factor
  use ovar,    only : ave_start_time,output_time_averaging
  use ovar,    only : grid_scaling_factor
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: called_from_channel_mod
  !
  !.. Local Scalars ..
  integer  :: j,m,n,ierr,io,ip
  integer  :: nvar,ny
  real(wp) :: out_time,rnum,ymax
  !
  integer  :: nvx,nvxvx,nvy,nvyvy,nvz,nvzvz,nvxvy,nrho,nt
  !
  logical(lk) :: averaging_channel_flow
  !
  character(len=50)  :: varname
  character(len=100) :: avefile
  character(len=100) :: proffile
  !
  !.. Local Arrays ..
  real(wp) :: y_vals(1:10000)
  real(wp) :: var_dim(1:size(uavesp,dim=2))
  real(wp) :: aux_data(1:17)
  !
  character(len=50) :: varnames(1:size(uavesp,dim=2))
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: y_ptr(:)
  integer,  allocatable :: y_counts(:)
  real(wp), allocatable :: y_data(:,:)
  real(wp), allocatable :: y_prof(:,:)
  real(wp), allocatable :: y_dudy(:)
  real(wp), allocatable :: dudy(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "spatial_ave_x_and_z_planes"
  !
continue
  !
  averaging_channel_flow = fals
  if (present(called_from_channel_mod)) then
    averaging_channel_flow = called_from_channel_mod
  end if
  !
  nvar = size(uavesp,dim=2)
  !
  !#############################################################################
  !
  ! Get the variable names
  !
  call get_variable_names(varnames,var_dim,nvx,nvxvx,nvy,nvyvy, &
                          nvz,nvzvz,nvxvy,nrho,nt)
  !
  !#############################################################################
  !
  ! Get du/dy
  !
  dudy = init_array(n_solpts,zero) ! F2003 AUTO-REALLOC
  !
  call get_dudy(uavesp(:,nvx),dudy)
  !
  !#############################################################################
  !
  ! Find the unique y coordinates
  !
  y_ptr = init_array(n_solpts,0) ! F2003 AUTO-REALLOC
  !
  y_vals(:) = zero
  !
  ny = 1
  y_vals(1) = xyz(2,1)
  y_ptr(1) = 1
  !
  solpts_loop: do n = 2,n_solpts
    !
    ! If the y-coordinate for the current solution point is equal to any of
    ! the previously found unique y coordinates in y_vals, set y_ptr for the
    ! current solution point to the index within y_vals of the matching
    ! y coordinate value, then continue on to the next solution point.
    !
    do j = 1,ny
      if (abs(y_vals(j)-xyz(2,n)) <= match_tolerance) then
        y_ptr(n) = j
        cycle solpts_loop
      end if
    end do
    !
    ! If we are here, that means the y-coordinate for the current solution
    ! point doesnt exist yet in y_vals, so we need to update the counter
    ! ny for number of unique y coordinate values, save the unique y coordinate
    ! value to y_vals, and set y_ptr for this solution point to the new
    ! counter value.
    !
    ny = ny + 1
    y_vals(ny) = xyz(2,n)
    y_ptr(n) = ny
    !
  end do solpts_loop
  !
  ymax = maxval(y_vals)
  !
  !#############################################################################
  !
  ! Collect the averaged data into an array
  !
  y_data   = init_array(nvar,ny,zero) ! F2003 AUTO-REALLOC
  y_counts = init_array(ny,0)         ! F2003 AUTO-REALLOC
  y_dudy   = init_array(ny,zero)      ! F2003 AUTO-REALLOC
  y_prof   = init_array(13,ny,zero)    ! F2003 AUTO-REALLOC
  !
  do n = 1,n_solpts
    !
    j = y_ptr(n)
    !
    y_counts(j) = y_counts(j) + 1
    !
    do m = 1,nvar
      y_data(m,j) = y_data(m,j) + uavesp(n,m)
    end do
    y_dudy(j) = y_dudy(j) + dudy(n)
    !
  end do
  !
  write (iout,2) minval(y_counts),maxval(y_counts)
  !
  !#############################################################################
  !
  ! Average the data by dividing through by the number of averaged points
  !
  do j = 1,ny
    rnum = one / real(y_counts(j),kind=wp)
    do m = 1,nvar
      y_data(m,j) = y_data(m,j) * rnum
    end do
    y_dudy(j) = y_dudy(j) * rnum
  end do
  !
  if (averaging_channel_flow) then
    rnum = one/ymax
  else
    rnum = one/grid_scaling_factor
  end if
  !
  !#############################################################################
  !
  ! Get the normalized profile
  !
  call get_normalized_profile(y_vals(1:ny), &
                              y_data(nvx  ,1:ny), &
                              y_data(nvxvx,1:ny), &
                              y_data(nvy  ,1:ny), &
                              y_data(nvyvy,1:ny), &
                              y_data(nvz  ,1:ny), &
                              y_data(nvzvz,1:ny), &
                              y_data(nvxvy,1:ny), &
                              y_data(nrho ,1:ny), &
                              y_data(nt   ,1:ny), &
                              y_dudy(1:ny), &
                              y_prof(:,1:ny), &
                              aux_data(:), &
                              averaging_channel_flow,1)
  !
  !#############################################################################
  !
  ! Output the spatially-averaged time-averaged solution
  !
  write (avefile,20) itcur
  avefile = trim(adjustl(avefile))
  !
  open (newunit=io,file=avefile,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,avefile,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (proffile,40) itcur
  proffile = trim(adjustl(proffile))
  !
  open (newunit=ip,file=proffile,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,proffile,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (io,10)
  write (ip,41)
  !
  do m = 1,nvar
    !
    varname = trim(adjustl(varnames(m)))
    !
    if (m == nvar) then
      write (io,11) trim(adjustl(varname))
    else
      write (io,11) trim(adjustl(varname)),","
    end if
    !
  end do
  !
  ! Write the averaged y-profile
  !
  write (io,12) ny
  !
  do j = 1,ny
    write (io,13) y_vals(j)*rnum,(y_data(m,j)*var_dim(m), m=1,nvar)
  end do
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,avefile,2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Write out the normalized y-profile
  !
  write (ip,12) ny
  write (ip,51) "r_avg",         aux_data(1)
  write (ip,52) "T_avg",         aux_data(2)
  write (ip,53) "visc_avg",      aux_data(3)
  write (ip,53) "dudy_wall",     aux_data(4)
  write (ip,51) "tau_wall",      aux_data(5)
  write (ip,51) "u_CL",          aux_data(6)
  write (ip,51) "u_tau",         aux_data(7)
  write (ip,51) "u_Bulk",        aux_data(8)
  write (ip,53) "dstar",         aux_data(9)
  write (ip,53) "theta",         aux_data(10)
  write (ip,54) "Re_Bulk",  nint(aux_data(11))
  write (ip,54) "Re_dstar", nint(aux_data(12))
  write (ip,54) "Re_theta", nint(aux_data(13))
  write (ip,54) "Re_CL",    nint(aux_data(14))
  write (ip,54) "Re_tau",   nint(aux_data(15))
  write (ip,53) "Cf_Bulk",       aux_data(16)
  write (ip,53) "Cf_CL",         aux_data(17)
  !
  do j = 1,ny
    write (ip,13) (y_prof(m,j), m=1,size(y_prof,dim=1))
  end do
  !
  close (ip,iostat=ierr,iomsg=error_message)
  call io_error(pname,proffile,2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Output the current iteration and how many
  ! flow-through-times have been averaged
  !
  if (mypnum == glb_root) then
    !
    out_time = time*time_ref*time_scaling_factor
    if (output_time_averaging) then
      out_time = out_time - ave_start_time*time_ref*time_scaling_factor
    end if
    !
    write (iout,21) itcur,out_time
    !
  end if
  !
  !
  write (error_message,1)
  call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
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
  20 format ("spatial_ave.",i0,".dat")
  21 format (/,2x,"Number of iterations completed = ",i0,/, &
               2x,"Number of flow-through-times averaged = ",f9.4,/)
  !
  40 format ("spatial_ave.profile.",i0,".dat")
  41 format ('VARIABLES = "y+",',/, &
             '"y/<greek>d</greek>",',/, &
             '"u/u<sub>cl</sub>",',/, &
             '"u+",',/, &
             '"u<sub>rms</sub>/u<sub><greek>t</greek></sub>",',/, &
             '"v<sub>rms</sub>/u<sub><greek>t</greek></sub>",',/, &
             '"w<sub>rms</sub>/u<sub><greek>t</greek></sub>",',/, &
             '"u<math>\\034</math>v<math>\\034\</math>/', &
             'u<sub><greek>t</greek></sub><sup>2</sup>",',/, &
             '"<greek>t</greek><sub>xx</sub>",',/, &
             '"<greek>t</greek><sub>yy</sub>",',/, &
             '"<greek>t</greek><sub>zz</sub>",',/, &
             '"<greek>t</greek><sub>xy</sub>",',/, &
             '"u"')
  !
  51 format ('AUXDATA ',a,'="',f9.6,'"')
  52 format ('AUXDATA ',a,'="',f8.3,'"')
  53 format ('AUXDATA ',a,'="',es14.6,'"')
  54 format ('AUXDATA ',a,'="',i0,'"')
  !
end subroutine spatial_ave_x_and_z_planes
!
!###############################################################################
!
subroutine spatial_ave_z_planes(called_from_channel_mod)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts,xyz
  use flowvar, only : uavesp
  use ovar,    only : itcur,time,time_ref,time_scaling_factor
  use ovar,    only : ave_start_time,output_time_averaging
  use ovar,    only : grid_scaling_factor
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: called_from_channel_mod
  !
  !.. Local Scalars ..
  integer  :: i,j,m,n,ierr,io,ia1,ia2,ip
  integer  :: nvar,ny,nx,nix
  real(wp) :: out_time,rnum,xmin,xmax,ymax
  !
  integer  :: nvx,nvxvx,nvy,nvyvy,nvz,nvzvz,nvxvy,nrho,nt
  !
  logical(lk) :: averaging_channel_flow
  !
  character(len=50)  :: varname
  character(len=100) :: avefile
  character(len=100) :: axialfile
  character(len=100) :: axfullfile
  character(len=100) :: proffile
  character(len=100) :: auxfile
  !
  !.. Local Arrays ..
  real(wp) :: x_vals(1:10000)
  real(wp) :: y_vals(1:10000)
  real(wp) :: var_dim(1:size(uavesp,dim=2))
  !
  character(len=50) :: varnames(1:size(uavesp,dim=2))
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: x_ptr(:)
  integer,  allocatable :: y_ptr(:)
  integer,  allocatable :: ix_ptr(:)
  integer,  allocatable :: xy_counts(:,:)
  real(wp), allocatable :: xy_data(:,:,:)
  real(wp), allocatable :: xy_prof(:,:,:)
  real(wp), allocatable :: xy_dudy(:,:)
  real(wp), allocatable :: dudy(:)
  real(wp), allocatable :: aux_data(:,:)
  real(wp), allocatable :: sem_profile(:,:)
  real(wp), allocatable :: aux_sem_data(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "spatial_ave_z_planes"
  !
continue
  !
  averaging_channel_flow = fals
  if (present(called_from_channel_mod)) then
    averaging_channel_flow = called_from_channel_mod
  end if
  !
  nvar = size(uavesp,dim=2)
  !
  !#############################################################################
  !
  ! Get the variable names
  !
  call get_variable_names(varnames,var_dim,nvx,nvxvx,nvy,nvyvy, &
                          nvz,nvzvz,nvxvy,nrho,nt)
  !
  !#############################################################################
  !
  ! Get du/dy
  !
  dudy = init_array(n_solpts,zero) ! F2003 AUTO-REALLOC
  !
  call get_dudy(uavesp(:,nvx),dudy)
  !
  !#############################################################################
  !
  ! Find the unique x coordinates
  !
  x_ptr = init_array(n_solpts,0) ! F2003 AUTO-REALLOC
  !
  x_vals(:) = zero
  !
  nx = 1
  x_vals(1) = xyz(1,1)
  x_ptr(1) = 1
  !
  find_unique_x_loc_loop: &
  do n = 2,n_solpts
    !
    do i = 1,nx
      if (abs(x_vals(i)-xyz(1,n)) <= match_tolerance) then
        x_ptr(n) = i
        cycle find_unique_x_loc_loop
      end if
    end do
    !
    nx = nx + 1
    x_vals(nx) = xyz(1,n)
    x_ptr(n) = nx
    !
  end do find_unique_x_loc_loop
  !
  xmin = minval(x_vals)
  xmax = maxval(x_vals)
  !
  !#############################################################################
  !
  ! Find the unique y coordinates
  !
  y_ptr = init_array(n_solpts,0) ! F2003 AUTO-REALLOC
  !
  y_vals(:) = zero
  !
  ny = 1
  y_vals(1) = xyz(2,1)
  y_ptr(1) = 1
  !
  find_unique_y_loc_loop: &
  do n = 2,n_solpts
    !
    do j = 1,ny
      if (abs(y_vals(j)-xyz(2,n)) <= match_tolerance) then
        y_ptr(n) = j
        cycle find_unique_y_loc_loop
      end if
    end do
    !
    ny = ny + 1
    y_vals(ny) = xyz(2,n)
    y_ptr(n) = ny
    !
  end do find_unique_y_loc_loop
  !
  ymax = maxval(y_vals)
  !
  !#############################################################################
  !
  ! Collect the averaged data into an x-y structured array
  !
  xy_data   = init_array(nvar,ny,nx,zero) ! F2003 AUTO-REALLOC
  xy_counts = init_array(ny,nx,0)         ! F2003 AUTO-REALLOC
  xy_dudy   = init_array(ny,nx,zero)      ! F2003 AUTO-REALLOC
  xy_prof   = init_array(13,ny,nx,zero)    ! F2003 AUTO-REALLOC
  aux_data  = init_array(17,nx,zero)       ! F2003 AUTO-REALLOC
  !
  do n = 1,n_solpts
    !
    i = x_ptr(n)
    j = y_ptr(n)
    !
    xy_counts(j,i) = xy_counts(j,i) + 1
    !
    do m = 1,nvar
      xy_data(m,j,i) = xy_data(m,j,i) + uavesp(n,m)
    end do
    xy_dudy(j,i) = xy_dudy(j,i) + dudy(n)
    !
  end do
  !
  write (iout,2) minval(xy_counts),maxval(xy_counts)
  !
  !#############################################################################
  !
  ! Average the data by dividing through by the number of averaged points
  !
  do i = 1,nx
    do j = 1,ny
      rnum = one / real(xy_counts(j,i),kind=wp)
      do m = 1,nvar
        xy_data(m,j,i) = xy_data(m,j,i) * rnum
      end do
      xy_dudy(j,i) = xy_dudy(j,i) * rnum
    end do
  end do
  !
  !#############################################################################
  !
  ! Get the spatially averaged profiles at axial locations
  ! equal to integer multiples of the channel height.
  !
  if (averaging_channel_flow) then
    rnum = one/ymax
  else
    rnum = one/grid_scaling_factor
  end if
  !
  nix = nint( abs(xmax-xmin)*rnum ) + 1
  !
  ix_ptr = init_array(nix,0) ! F2003 AUTO-REALLOC
  !
  do n = 1,nix
    ix_ptr(n) = minloc( abs(abs(x_vals(1:nx)-xmin)*rnum - real(n-1,kind=wp)) , dim=1 )
  end do
  write (iout,*) "xmin = ",xmin
  write (iout,*) "xmax = ",xmax
  write (iout,*) "abs(xmax-xmin)*rnum = ",abs(xmax-xmin)*rnum
  write (iout,*) "nint(abs(xmax-xmin)*rnum) = ",nint(abs(xmax-xmin)*rnum)
  do n = 1,nix
    write (iout,201) n,ix_ptr(n),x_vals(ix_ptr(n)),x_vals(ix_ptr(n))*rnum
  end do
  201 format (2x,i2,": x_vals(",i0,") = ",f10.5,2x,f10.5)
  !
  !#############################################################################
  !
  ! Get the normalized profiles
  !
  do i = 1,nx
    call get_normalized_profile(y_vals(1:ny), &
                                xy_data(nvx  ,1:ny,i), &
                                xy_data(nvxvx,1:ny,i), &
                                xy_data(nvy  ,1:ny,i), &
                                xy_data(nvyvy,1:ny,i), &
                                xy_data(nvz  ,1:ny,i), &
                                xy_data(nvzvz,1:ny,i), &
                                xy_data(nvxvy,1:ny,i), &
                                xy_data(nrho ,1:ny,i), &
                                xy_data(nt   ,1:ny,i), &
                                xy_dudy(1:ny,i), &
                                xy_prof(:,1:ny,i), &
                                aux_data(:,i), &
                                averaging_channel_flow,ix=i)
  end do
  !
  ! Get the SEM profile
  !
  call get_sem_profile(sem_profile,aux_sem_data,averaging_channel_flow)
  !
  !#############################################################################
  !
  ! Output the spatially-averaged time-averaged solution
  !
  write (avefile,20) itcur
  avefile = trim(adjustl(avefile))
  !
  open (newunit=io,file=avefile,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,avefile,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (axialfile,30) "",itcur
  axialfile = trim(adjustl(axialfile))
  !
  open (newunit=ia1,file=axialfile,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,axialfile,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (axfullfile,30) "-full",itcur
  axfullfile = trim(adjustl(axfullfile))
  !
  open (newunit=ia2,file=axfullfile,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,axfullfile,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (proffile,40) itcur
  proffile = trim(adjustl(proffile))
  !
  open (newunit=ip,file=proffile,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,proffile,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (io,10)
  write (ia1,31)
  write (ia2,31)
  write (ip,41)
  !
  do m = 1,nvar
    !
    varname = trim(adjustl(varnames(m)))
    !
    if (m == nvar) then
      write (io,11) trim(adjustl(varname))
      write (ia1,11) trim(adjustl(varname))
      write (ia2,11) trim(adjustl(varname))
    else
      write (io,11) trim(adjustl(varname)),","
      write (ia1,11) trim(adjustl(varname)),","
      write (ia2,11) trim(adjustl(varname)),","
    end if
    !
  end do
  !
  do i = 1,nx
    do j = 1,ny
      do m = 1,nvar
        xy_data(m,j,i) = xy_data(m,j,i)*var_dim(m)
      end do
    end do
  end do
  !
  ! Write the averaged x-y planes file
  !
  write (io,12) nx,ny
  !
  do j = 1,ny
    do i = 1,nx
      write (io,13) x_vals(i)*rnum, &
                    y_vals(j)*rnum, &
                    (xy_data(m,j,i), m=1,nvar)
    end do
  end do
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,avefile,2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Write the averaged y-profiles at each axial location equal to
  ! an integer interval of the channel height
  !
  do n = 1,nix
    !
    i = ix_ptr(n)
    !
    if (averaging_channel_flow) then
      write (ia1,133) nint(x_vals(i)*rnum),ny
    else
      write (ia1,33) nint(x_vals(i)*rnum),ny
    end if
    !
    do j = 1,ny
      write (ia1,13) y_vals(j)*rnum,(xy_data(m,j,i), m=1,nvar)
    end do
    !
  end do
  !
  close (ia1,iostat=ierr,iomsg=error_message)
  call io_error(pname,axialfile,2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Write the averaged y-profiles at every axial location
  !
  do i = 1,nx
    !
    if (averaging_channel_flow) then
      write (ia2,132) x_vals(i)*rnum,ny
    else
      write (ia2,32) x_vals(i)*rnum,ny
    end if
    !
    do j = 1,ny
      write (ia2,13) y_vals(j)*rnum,(xy_data(m,j,i), m=1,nvar)
    end do
    !
  end do
  !
  close (ia2,iostat=ierr,iomsg=error_message)
  call io_error(pname,axfullfile,2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Write the normalized y-profile at each axial location equal to
  ! an integer interval of the channel height
  !
  do n = 1,nix
    !
    i = ix_ptr(n)
    !
    if (averaging_channel_flow) then
      write (ip,133) nint(x_vals(i)*rnum),ny
    else
      write (ip,33) nint(x_vals(i)*rnum),ny
    end if
    write (ip,51) "r_avg",         aux_data(1,i)
    write (ip,52) "T_avg",         aux_data(2,i)
    write (ip,53) "visc_avg",      aux_data(3,i)
    write (ip,53) "dudy_wall",     aux_data(4,i)
    write (ip,53) "tau_wall",      aux_data(5,i)
    write (ip,52) "u_CL",          aux_data(6,i)
    write (ip,52) "u_tau",         aux_data(7,i)
    write (ip,52) "u_Bulk",        aux_data(8,i)
    write (ip,53) "dstar",         aux_data(9,i)
    write (ip,53) "theta",         aux_data(10,i)
    write (ip,54) "Re_Bulk",  nint(aux_data(11,i))
    write (ip,54) "Re_dstar", nint(aux_data(12,i))
    write (ip,54) "Re_theta", nint(aux_data(13,i))
    write (ip,54) "Re_CL",    nint(aux_data(14,i))
    write (ip,54) "Re_tau",   nint(aux_data(15,i))
    write (ip,53) "Cf_Bulk",       aux_data(16,i)
    write (ip,53) "Cf_CL",         aux_data(17,i)
    !
    do j = 1,ny
      write (ip,13) (xy_prof(m,j,i), m=1,size(xy_prof,dim=1))
    end do
    !
  end do
  !
  if (allocated(sem_profile) .and. allocated(aux_sem_data)) then
    !
    ! Write the SEM profile
    !
    write (ip,34) size(sem_profile,dim=2)
    write (ip,51) "r_avg",         aux_sem_data(1)
    write (ip,52) "T_avg",         aux_sem_data(2)
    write (ip,53) "visc_avg",      aux_sem_data(3)
    write (ip,53) "dudy_wall",     aux_sem_data(4)
    write (ip,53) "tau_wall",      aux_sem_data(5)
    write (ip,52) "u_CL",          aux_sem_data(6)
    write (ip,52) "u_tau",         aux_sem_data(7)
    write (ip,52) "u_Bulk",        aux_sem_data(8)
    write (ip,53) "dstar",         aux_sem_data(9)
    write (ip,53) "theta",         aux_sem_data(10)
    write (ip,54) "Re_Bulk",  nint(aux_sem_data(11))
    write (ip,54) "Re_dstar", nint(aux_sem_data(12))
    write (ip,54) "Re_theta", nint(aux_sem_data(13))
    write (ip,54) "Re_CL",    nint(aux_sem_data(14))
    write (ip,54) "Re_tau",   nint(aux_sem_data(15))
    write (ip,53) "Cf_Bulk",       aux_sem_data(16)
    write (ip,53) "Cf_CL",         aux_sem_data(17)
    !
    do j = 1,size(sem_profile,dim=2)
      write (ip,13) (sem_profile(m,j), m=1,size(sem_profile,dim=1))
    end do
    !
  end if
  !
  close (ip,iostat=ierr,iomsg=error_message)
  call io_error(pname,proffile,2,__LINE__,__FILE__,ierr,error_message)
  !
  !
  write (auxfile,50) itcur
  auxfile = trim(adjustl(auxfile))
  !
  open (newunit=io,file=auxfile,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,auxfile,1,__LINE__,__FILE__,ierr,error_message)
  !
  if (averaging_channel_flow) then
    write (io,160) nx
  else
    write (io,60) nx
  end if
  !
  do i = 1,nx
    write (io,13) x_vals(i)*rnum,(aux_data(m,i), m=1,size(aux_data,dim=1))
  end do
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,auxfile,2,__LINE__,__FILE__,ierr,error_message)
  !
  if (size(xy_prof,dim=2) == size(sem_profile,dim=2)) then
    x_vals(:) = x_vals(:)*rnum
    call write_profile_errors_file(x_vals,ix_ptr,xy_prof,sem_profile, &
                                   aux_data,aux_sem_data)
  end if
  !
  !#############################################################################
  !
  ! Output the current iteration and how many
  ! flow-through-times have been averaged
  !
  if (mypnum == glb_root) then
    !
    out_time = time*time_ref*time_scaling_factor
    if (output_time_averaging) then
      out_time = out_time - ave_start_time*time_ref*time_scaling_factor
    end if
    !
    write (iout,21) itcur,out_time
    !
  end if
  !
  !#############################################################################
  !
  if (allocated(x_ptr)) then
    deallocate ( x_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"x_ptr",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(y_ptr)) then
    deallocate ( y_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"y_ptr",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(ix_ptr)) then
    deallocate ( ix_ptr , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"ix_ptr",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(xy_counts)) then
    deallocate ( xy_counts , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xy_counts",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(xy_data)) then
    deallocate ( xy_data , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xy_data",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(xy_prof)) then
    deallocate ( xy_prof , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xy_prof",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(xy_dudy)) then
    deallocate ( xy_dudy , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"xy_dudy",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(dudy)) then
    deallocate ( dudy , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"dudy",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(aux_data)) then
    deallocate ( aux_data , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"aux_data",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(sem_profile)) then
    deallocate ( sem_profile , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"sem_profile",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  if (allocated(aux_sem_data)) then
    deallocate ( aux_sem_data , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"aux_sem_data",2,__LINE__,__FILE__, &
                     ierr,error_message)
  end if
  !
  1 format ("After spatial averaging a channel flow restart file.")
  2 format (/,"  Points used for averaging: Min = ",i0,/, &
              "                             Max = ",i0,/)
  !
  10 format ('VARIABLES = "X","Y",')
  11 format ('"',a,'"',:,a)
  12 format ('ZONE T="Spatial Average", I=',i0,', J=',i0,', DATAPACKING=POINT')
  13 format (100es14.6)
  !
  20 format ("spatial_ave.z-planes.",i0,".dat")
  21 format (/,2x,"Number of iterations completed = ",i0,/, &
               2x,"Number of flow-through-times averaged = ",f9.4,/)
  !
  30 format ("spatial_ave.axial-locations",a,".",i0,".dat")
  31 format ('VARIABLES = "Y",')
  32 format ('ZONE T="Spatial average profile at x =',f6.3, &
             '", I=',i0,', DATAPACKING=POINT')
  33 format ('ZONE T="Spatial average profile at x = ',i0, &
             '", I=',i0,', DATAPACKING=POINT')
  34 format ('ZONE T="SEM Inflow Profile", I=',i0,', DATAPACKING=POINT')
  !
  40 format ("spatial_ave.axial-profiles.",i0,".dat")
  41 format ('VARIABLES = "y+",',/, &
             '"y/<greek>d</greek>",',/, &
             '"u/u<sub>cl</sub>",',/, &
             '"u+",',/, &
             '"u<sub>rms</sub>/u<sub><greek>t</greek></sub>",',/, &
             '"v<sub>rms</sub>/u<sub><greek>t</greek></sub>",',/, &
             '"w<sub>rms</sub>/u<sub><greek>t</greek></sub>",',/, &
             '"u<math>\\034</math>v<math>\\034\</math>/', &
             'u<sub><greek>t</greek></sub><sup>2</sup>",',/, &
             '"<greek>t</greek><sub>xx</sub>",',/, &
             '"<greek>t</greek><sub>yy</sub>",',/, &
             '"<greek>t</greek><sub>zz</sub>",',/, &
             '"<greek>t</greek><sub>xy</sub>",',/, &
             '"u"')
  !
  50 format ("spatial_ave.aux-data.",i0,".dat")
  51 format ('AUXDATA ',a,'="',f9.6,'"')
  52 format ('AUXDATA ',a,'="',f8.3,'"')
  53 format ('AUXDATA ',a,'="',es14.6,'"')
  54 format ('AUXDATA ',a,'="',i0,'"')
  !
  60 format ('VARIABLES = "x",',/, &
                         '"r_avg",',/, &
                         '"T_avg",',/, &
                         '"visc_avg",',/, &
                         '"dudy_wall",',/, &
                         '"tau_wall",',/, &
                         '"u_CL",',/, &
                         '"u_tau",',/, &
                         '"u_Bulk",',/, &
                         '"dstar",',/, &
                         '"theta",',/, &
                         '"Re_Bulk",',/, &
                         '"Re_dstar",',/, &
                         '"Re_theta",',/, &
                         '"Re_CL",',/, &
                         '"Re_tau",',/, &
                         '"Cf_Bulk",',/, &
                         '"Cf_CL"',/, &
             'ZONE T="Axial Profile", I=',i0,', DATAPACKING=POINT')
  !
  132 format ('ZONE T="Spatial average profile at x =',f6.3, &
              'H", I=',i0,', DATAPACKING=POINT')
  133 format ('ZONE T="Spatial average profile at x = ',i0, &
              ' H", I=',i0,', DATAPACKING=POINT')
  160 format ('VARIABLES = "x/H",',/, &
                          '"r_avg",',/, &
                          '"T_avg",',/, &
                          '"visc_avg",',/, &
                          '"dudy_wall",',/, &
                          '"tau_wall",',/, &
                          '"u_CL",',/, &
                          '"u_tau",',/, &
                          '"u_Bulk",',/, &
                          '"dstar",',/, &
                          '"theta",',/, &
                          '"Re_Bulk",',/, &
                          '"Re_dstar",',/, &
                          '"Re_theta",',/, &
                          '"Re_CL",',/, &
                          '"Re_tau",',/, &
                          '"Cf_Bulk",',/, &
                          '"Cf_CL"',/, &
              'ZONE T="Axial Profile", I=',i0,', DATAPACKING=POINT')
  !
end subroutine spatial_ave_z_planes
!
!###############################################################################
!
pure &
subroutine get_variable_names(varnames,var_dim,nvx,nvxvx,nvy,nvyvy, &
                              nvz,nvzvz,nvxvy,nrho,nt)
  !
  !.. Use Statements ..
  use ovar,    only : rhoref,aref,tref,tkeref,tlsref
  use flowvar, only : time_ave_variables
  !
  !.. Formal Arguments ..
  integer,          intent(inout) :: nvx
  integer,          intent(inout) :: nvxvx
  integer,          intent(inout) :: nvy
  integer,          intent(inout) :: nvyvy
  integer,          intent(inout) :: nvz
  integer,          intent(inout) :: nvzvz
  integer,          intent(inout) :: nvxvy
  integer,          intent(inout) :: nrho
  integer,          intent(inout) :: nt
  character(len=*), intent(inout) :: varnames(:)
  real(wp),         intent(inout) :: var_dim(:)
  !
  !.. Local Scalars ..
  integer :: l,m,nvar
  !
  character(len=1) :: ta_char
  character(len=1) :: vmult
  !
  character(len=len(varnames)) :: varname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_variable_names"
  !
continue
  !
  nvx = 0
  nvxvx = 0
  nvy = 0
  nvyvy = 0
  nvz = 0
  nvzvz = 0
  nvxvy = 0
  nrho = 0
  nt = 0
  !
  nvar = size(var_dim)
  !
  var_dim(1:nvar) = one
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
    select case (uppercase(trim(adjustl(varname))))
      case ("VX")
        nvx = m
      case ("VX*VX")
        nvxvx = m
      case ("VY")
        nvy = m
      case ("VY*VY")
        nvyvy = m
      case ("VZ")
        nvz = m
      case ("VZ*VZ")
        nvzvz = m
      case ("VX*VY")
        nvxvy = m
      case ("RHO")
        nrho = m
      case ("T")
        nt = m
    end select
    !
    varnames(m) = trim(adjustl(varname))
    !
  end do
  !
end subroutine get_variable_names
!
!###############################################################################
!
pure &
subroutine get_dudy(vx,dudy)
  !
  !.. Use Statements ..
  use geovar,          only : ncell,cell
  use metrics_mod,     only : metrics_array
  use derivatives_mod, only : deriv
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: vx(:)
  real(wp), intent(inout) :: dudy(:)
  !
  !.. Local Scalars ..
  integer  :: i,j,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,dudi,dudj,dudk
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_dudy"
  !
continue
  !
  do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp-1  ! beginning index for solpts in cell nc
    np = cell(nc)%end_sp-n0 ! number of solution points in this cell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Compute the derivative of vx with respect to y
    !
    do n = 1,np
      !
      dudi = zero
      dudj = zero
      dudk = zero
      !
      ! Compute du/dxi
      !
      do i = deriv(this_geom,this_order)%d(1)%row_ptr(n)+1, &
             deriv(this_geom,this_order)%d(1)%row_ptr(n+1)
        !
        j = deriv(this_geom,this_order)%d(1)%row_idx(i)
        a = deriv(this_geom,this_order)%d(1)%row_val(i)
        !
        dudi = dudi + a * vx(n0+j)
        !
      end do
      !
      ! Compute du/deta
      !
      do i = deriv(this_geom,this_order)%d(2)%row_ptr(n)+1, &
             deriv(this_geom,this_order)%d(2)%row_ptr(n+1)
        !
        j = deriv(this_geom,this_order)%d(2)%row_idx(i)
        a = deriv(this_geom,this_order)%d(2)%row_val(i)
        !
        dudj = dudj + a * vx(n0+j)
        !
      end do
      !
      ! Compute du/dzeta
      !
      do i = deriv(this_geom,this_order)%d(3)%row_ptr(n)+1, &
             deriv(this_geom,this_order)%d(3)%row_ptr(n+1)
        !
        j = deriv(this_geom,this_order)%d(3)%row_idx(i)
        a = deriv(this_geom,this_order)%d(3)%row_val(i)
        !
        dudk = dudk + a * vx(n0+j)
        !
      end do
      !
      ! Finally, compute d(vx)/d(y) by transforming the reference space
      ! derivatives computed above and the metrics terms
      !
      !   du/dy = (du/dxi)*(dxi/dy) + 
      !           (du/deta)*(deta/dy) + 
      !           (du/dzeta)*(dzeta/dy)
      !
      dudy(n0+n) = dudi*metrics_array(n0+n,1,2) + &
                   dudj*metrics_array(n0+n,2,2) + &
                   dudk*metrics_array(n0+n,3,2)
      !
    end do
    !
  end do
  !
end subroutine get_dudy
!
!###############################################################################
!
!pure &
subroutine get_normalized_profile(y_vals,u,uu,v,vv,w,ww,uv,rho,t,dudy, &
                                  prof,aux_data,averaging_channel_flow,ix, &
                                  opt_debug_output)
  !
  !.. Use Statements ..
  use order_mod,        only : n_order
  use ovar,             only : rhoref,aref,tref
  use ovar,             only : grid_scaling_factor
  use element_type_mod, only : std_elem
  use transport_mod,    only : compute_viscosity_cpg
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: y_vals(:)
  real(wp),    intent(in) :: u(:)
  real(wp),    intent(in) :: uu(:)
  real(wp),    intent(in) :: v(:)
  real(wp),    intent(in) :: vv(:)
  real(wp),    intent(in) :: w(:)
  real(wp),    intent(in) :: ww(:)
  real(wp),    intent(in) :: uv(:)
  real(wp),    intent(in) :: rho(:)
  real(wp),    intent(in) :: t(:)
  real(wp),    intent(in) :: dudy(:)
  real(wp), intent(inout) :: prof(:,:)
  real(wp), intent(inout) :: aux_data(:)
  logical(lk), intent(in) :: averaging_channel_flow
  integer,     intent(in) :: ix
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: opt_debug_output
  !
  !.. Local Scalars ..
  integer  :: n,nc,n0,npts,n1d
  real(wp) :: rnum,height,half_height,dudy_wall,gsf
  real(wp) :: rho_avg,t_avg,visc_avg,u_cl,u_tau,u_bulk
  real(wp) :: xi_min,xi_max,y1,y2,tau_wall
  real(wp) :: dy,wt,udt,umt,udtm,udtp,umtm,umtp
  real(wp) :: displacement_thickness,momentum_thickness
  !
  logical(lk) :: debug_output
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_normalized_profile"
  !
continue
  !
  debug_output = fals
  if (present(opt_debug_output)) debug_output = opt_debug_output
  !
  n1d = n_order+1
  !
  npts = size(u)
  rnum = one / real(npts,kind=wp)
  !
  ! Get the channel height
  ! NOTE: The min y-value is slightly larger than zero, and the max y-value
  !       is slightly less than the actually channel height. Since these two
  !       y-locations are the same distance to their respective walls, we can
  !       just add these two values together to get the true channel height.
  !
  if (averaging_channel_flow) then
    !
    height = y_vals(npts) + y_vals(1)
    half_height = half * height
    gsf = one
    !
    rho_avg = sum(rho(:)) * rnum
    t_avg = sum(t(:)) * rnum
    !
    if (mod(npts,2) == 0) then
      u_cl = max( u(npts/2) , u((npts/2)+1) )
    else
      u_cl = u((npts/2)+1)
    end if
    !
    dudy_wall = half*(abs(dudy(1)) + abs(dudy(npts)))
    !
  else
    !
    height = abs(y_vals(npts) - y_vals(1))
    half_height = one
    gsf = grid_scaling_factor
    !
    rho_avg = rho(1)
    t_avg = t(1)
    !
    u_cl = u(npts)
    !
    dudy_wall = dudy(1)
    !
  end if
  !
  call compute_viscosity_cpg(t_avg,visc_avg)
  !
  tau_wall = visc_avg * dudy_wall
  !
  u_tau = sqrt( tau_wall / rho_avg )
  !
  do n = 1,npts
    prof( 1,n) = y_vals(n) * u_tau * rho_avg / visc_avg ! y+
    prof( 2,n) = y_vals(n) / half_height / gsf          ! y/b
    prof( 3,n) = u(n) / u_cl                            ! u/u_cl
    prof( 4,n) = u(n) / u_tau                           ! u+
    prof( 5,n) = sqrt(max(zero,uu(n) - u(n)*u(n))) / u_tau        ! u'/u_tau
    prof( 6,n) = sqrt(max(zero,vv(n) - v(n)*v(n))) / u_tau        ! v'/u_tau
    prof( 7,n) = sqrt(max(zero,ww(n) - w(n)*w(n))) / u_tau        ! w'/u_tau
    prof( 8,n) = (uv(n) - u(n)*v(n)) / (u_tau*u_tau)    ! u'v'/(u_tau**2)
    prof( 9,n) = (uu(n) - u(n)*u(n))*aref*aref          ! tauxx
    prof(10,n) = (vv(n) - v(n)*v(n))*aref*aref          ! tauyy
    prof(11,n) = (ww(n) - w(n)*w(n))*aref*aref          ! tauzz
    prof(12,n) = (uv(n) - u(n)*v(n))*aref*aref          ! tauxy
    prof(13,n) = u(n)*aref                              ! u
  end do
  !
  ! Compute the bulk velocity (Eq. 7.3 from Turbulent Flows by Pope)
  !
  u_bulk = zero
  momentum_thickness = zero
  displacement_thickness = zero
  !
  if (mod(npts,n1d) == 0) then
    !
    ! Use quadrature to integrate the bulk velocity
    !
    xi_min = std_elem(Geom_Edge,n_order)%pts(1,1)
    xi_max = std_elem(Geom_Edge,n_order)%pts(1,n1d)
    !
    do nc = 1,npts/n1d
      !
      n0 = (nc-1)*n1d
      !
      y1 = y_vals(n0+1)
      y2 = y_vals(n0+n1d)
      !
      rnum = (y2 - y1) / (xi_max - xi_min) ! Simplified using Mathematica
      !
      do n = 1,n1d
        wt = rnum*std_elem(Geom_Edge,n_order)%wts(n)
        u_bulk = u_bulk + u(n0+n)*wt
        udt = (one - u(n0+n)/u_cl)
        umt = udt*u(n0+n)/u_cl
        displacement_thickness = displacement_thickness + udt*wt
        momentum_thickness = momentum_thickness + umt*wt
      end do
      !
    end do
    !
  else
    !
    ! Use trapezoid rule to integrate bulk velocity
    !
    do n = 1,npts-1
      dy = abs(y_vals(n+1)-y_vals(n))
      u_bulk = u_bulk + half*(u(n) + u(n+1))*dy
      udtm = (one - u(n  )/u_cl)
      udtp = (one - u(n+1)/u_cl)
      umtm = udtm*u(n  )/u_cl
      umtp = udtp*u(n+1)/u_cl
      displacement_thickness = displacement_thickness + half*(udtm+udtp)*dy
      momentum_thickness = momentum_thickness + half*(umtm+umtp)*dy
    end do
    !
  end if
  !
  u_bulk = u_bulk / height
  !
  ! Save the auxiliary data for this profile
  !
  aux_data(1) = rho_avg * rhoref
  aux_data(2) = t_avg * tref
  aux_data(3) = visc_avg * rhoref * aref
  aux_data(4) = dudy_wall * aref
  aux_data(5) = tau_wall * rhoref * aref * aref
  aux_data(6) = u_cl * aref
  aux_data(7) = u_tau * aref
  aux_data(8) = u_bulk * aref
  aux_data(9) = displacement_thickness
  aux_data(10) = momentum_thickness
  ! Re = Reynolds number based on bulk velocity (Pope eq. 7.1)
  aux_data(11) = rho_avg * u_bulk * height / visc_avg
  ! Retheta = Reynolds number based on displacement thickness (Pope eq. 7.126)
  aux_data(12) = rho_avg * u_cl * displacement_thickness / visc_avg
  ! Retheta = Reynolds number based on momentum thickness (Pope eq. 7.127)
  aux_data(13) = rho_avg * u_cl * momentum_thickness / visc_avg
  ! Re0 = Reynolds number based on center-line velocity (Pope eq. 7.2)
  aux_data(14) = rho_avg * u_cl * half_height / visc_avg
  ! Re_tau = Reynolds number based on frication velocity
  aux_data(15) = rho_avg * u_tau * half_height / visc_avg
  ! Skin friction based on bulk velocity (Pope eq. 7.15)
  aux_data(16) = tau_wall / (half * rho_avg * u_bulk * u_bulk)
  ! Skin friction based on center-line velocity (Pope eq. 7.14)
  aux_data(17) = tau_wall / (half * rho_avg * u_cl * u_cl)
  !
  if (debug_output) then
    if (mod(npts,n1d) == 0) then
      write (iout,*)
      write (iout,801) "u_Bulk",u_bulk*aref
      write (iout,901) " dstar",displacement_thickness
      write (iout,901) " theta",momentum_thickness
      write (iout,*)
    end if
    write (iout,*)
  end if
  801 format (2x,"Using     Quadrature:  ",a6," = ",f11.6," m/s")
  901 format (2x,"Using     Quadrature:  ",a6," = ",es14.6)
  !
end subroutine get_normalized_profile
!
!###############################################################################
!
!pure &
subroutine get_sem_profile(prof,aux_data,averaging_channel_flow)
  !
  !.. Use Statements ..
  use ovar,             only : rhoref,aref,tref
  use ovar,             only : grid_scaling_factor
  use order_mod,        only : n_order
  use transport_mod,    only : compute_viscosity_cpg
  use element_type_mod, only : std_elem
  !
  !.. Formal Arguments ..
  real(wp), allocatable, intent(out) :: prof(:,:)
  real(wp), allocatable, intent(out) :: aux_data(:)
  logical(lk),            intent(in) :: averaging_channel_flow
  !
  !.. Local Scalars ..
  integer  :: n,npts,ierr
  real(wp) :: rnum,height,half_height,dudy_wall,gsf
  real(wp) :: rho_avg,t_avg,visc_avg,u_cl,u_tau,u_bulk
  real(wp) :: tau_wall,yc,yh,ywall
  real(wp) :: dy,udtm,udtp,umtm,umtp
  real(wp) :: displacement_thickness,momentum_thickness
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: y(:)
  real(wp), allocatable :: u(:)
  real(wp), allocatable :: r(:)
  real(wp), allocatable :: t(:)
  real(wp), allocatable :: tau11(:)
  real(wp), allocatable :: tau22(:)
  real(wp), allocatable :: tau33(:)
  real(wp), allocatable :: tau21(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_sem_profile"
  !
continue
  !
  ierr = 0
  !
  call extract_sem_profile(y,u,r,t,tau11,tau22,tau33,tau21,ierr)
  !
  if (ierr /= 0) return
  !
  npts = size(y)
  !
  rnum = one / real(npts,kind=wp)
  !
  prof     = init_array(13,npts,zero)         ! F2003 AUTO-REALLOC
  aux_data = init_array(17,zero)              ! F2003 AUTO-REALLOC
  !
  ! Get the channel height
  ! NOTE: The min y-value is slightly larger than zero, and the max y-value
  !       is slightly less than the actually channel height. Since these two
  !       y-locations are the same distance to their respective walls, we can
  !       just add these two values together to get the true channel height.
  !
  if (averaging_channel_flow) then
    !
    height = y(npts) + y(1)
    half_height = half * height
    gsf = one
    !
    rho_avg = sum( r(:) ) * rnum
    t_avg   = sum( t(:) ) * rnum
    !
    if (mod(npts,2) == 0) then
      u_cl = max( u(npts/2) , u((npts/2)+1) )
    else
      u_cl = u((npts/2)+1)
    end if
    !
  else
    !
    height = abs(y(npts) - y(1))
    half_height = one
    gsf = grid_scaling_factor
    !
    rho_avg = r(1)
    t_avg = t(1)
    !
    u_cl = u(npts)
    !
  end if
  !
  ! The follow is needed in case the wall is not at y=0
  !
  ! Find center of first cell off the wall
  if (mod(n_order+1,2) == 0) then
    yc = half*( y( (n_order+1)/2 ) + y( (n_order+1)/2 + 1 ) )
  else
    yc = y( (n_order/2) + 1 )
  end if
  ! Find height of first cell off the wall
  yh = two*( (y(1)-yc) / std_elem(Geom_Edge,n_order)%pts(1,1) )
  ! Get the coordinates of the wall
  ywall = yc - half*yh
  ! Get the distance of the first point off the wall
  dy = y(1)-ywall
  !
  if (abs(dy) <= eps10) then
    dudy_wall = abs(u(1)/y(1))
  else
    dudy_wall = u(1)/dy
  end if
  !
  call compute_viscosity_cpg(t_avg,visc_avg)
  !
  tau_wall = visc_avg * dudy_wall
  !
  u_tau = sqrt( tau_wall / rho_avg )
  !
  do n = 1,npts
    prof( 1,n) = y(n) * u_tau * rho_avg / visc_avg  ! y+
    prof( 2,n) = y(n) / half_height / gsf           ! y/b
    prof( 3,n) = u(n) / u_cl                        ! u/u_cl
    prof( 4,n) = u(n) / u_tau                       ! u+
    prof( 5,n) = sqrt(tau11(n)) / u_tau             ! u'/u_tau
    prof( 6,n) = sqrt(tau22(n)) / u_tau             ! v'/u_tau
    prof( 7,n) = sqrt(tau33(n)) / u_tau             ! w'/u_tau
    prof( 8,n) = tau21(n) / (u_tau*u_tau)           ! u'v'/(u_tau**2)
    prof( 9,n) = tau11(n)*aref*aref                 ! tauxx
    prof(10,n) = tau22(n)*aref*aref                 ! tauyy
    prof(11,n) = tau33(n)*aref*aref                 ! tauzz
    prof(12,n) = tau21(n)*aref*aref                 ! tauxy
    prof(13,n) = u(n)*aref                          ! u
  end do
  !
  ! Compute the bulk velocity (Eq. 7.3 from Turbulent Flows by Pope)
  !
  ! Use trapezoid rule to integrate bulk velocity
  !
  u_bulk = zero
  momentum_thickness = zero
  displacement_thickness = zero
  !
  do n = 1,npts-1
    dy = abs(y(n+1)-y(n))
    u_bulk = u_bulk + half*(u(n) + u(n+1))*dy
    udtm = (one - u(n  )/u_cl)
    udtp = (one - u(n+1)/u_cl)
    umtm = udtm*u(n  )/u_cl
    umtp = udtp*u(n+1)/u_cl
    displacement_thickness = displacement_thickness + half*(udtm+udtp)*dy
    momentum_thickness = momentum_thickness + half*(umtm+umtp)*dy
  end do
  !
  u_bulk = u_bulk / height
  !
  ! Save the auxiliary data for this profile
  !
  aux_data(1) = rho_avg * rhoref
  aux_data(2) = t_avg * tref
  aux_data(3) = visc_avg * rhoref * aref
  aux_data(4) = dudy_wall * aref
  aux_data(5) = tau_wall * rhoref * aref * aref
  aux_data(6) = u_cl * aref
  aux_data(7) = u_tau * aref
  aux_data(8) = u_bulk * aref
  aux_data(9) = displacement_thickness
  aux_data(10) = momentum_thickness
  ! Re = Reynolds number based on bulk velocity (Pope eq. 7.1)
  aux_data(11) = rho_avg * u_bulk * height / visc_avg
  ! Retheta = Reynolds number based on displacement thickness (Pope eq. 7.126)
  aux_data(12) = rho_avg * u_cl * displacement_thickness / visc_avg
  ! Retheta = Reynolds number based on momentum thickness (Pope eq. 7.127)
  aux_data(13) = rho_avg * u_cl * momentum_thickness / visc_avg
  ! Re0 = Reynolds number based on center-line velocity (Pope eq. 7.2)
  aux_data(14) = rho_avg * u_cl * half_height / visc_avg
  ! Re_tau = Reynolds number based on frication velocity
  aux_data(15) = rho_avg * u_tau * half_height / visc_avg
  ! Skin friction based on bulk velocity (Pope eq. 7.15)
  aux_data(16) = tau_wall / (half * rho_avg * u_bulk * u_bulk)
  ! Skin friction based on center-line velocity (Pope eq. 7.14)
  aux_data(17) = tau_wall / (half * rho_avg * u_cl * u_cl)
  !
end subroutine get_sem_profile
!
!###############################################################################
!
pure &
subroutine extract_sem_profile(y,uavg,ravg,tavg,tau11,tau22,tau33,tau21,ierr)
  !
  !.. Use Statements ..
  use synthetic_eddy_mod, only : my_sem_pts
  !
  !.. Formal Arguments ..
  real(wp), allocatable, intent(out) :: y(:)
  real(wp), allocatable, intent(out) :: uavg(:)
  real(wp), allocatable, intent(out) :: ravg(:)
  real(wp), allocatable, intent(out) :: tavg(:)
  real(wp), allocatable, intent(out) :: tau11(:)
  real(wp), allocatable, intent(out) :: tau22(:)
  real(wp), allocatable, intent(out) :: tau33(:)
  real(wp), allocatable, intent(out) :: tau21(:)
  integer,               intent(out) :: ierr
  !
  !.. Local Scalars ..
  integer  :: j,n,ny,npts
  real(wp) :: rnum
  !
  !.. Local Arrays ..
  real(wp) :: y_vals(1:10000)
  !
  !.. Local Allocatable Arrays ..
  integer, allocatable :: y_ptr(:)
  integer, allocatable :: y_counts(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "extract_sem_profile"
  !
continue
  !
  ierr = 0
  !
  if (allocated(my_sem_pts)) then
    if (.not. allocated(my_sem_pts(1)%sem_fp)) ierr = 1
  else
    ierr = 1
  end if
  !
  if (ierr /= 0) return
  !
  npts = my_sem_pts(1)%num_fp
  !
  !#############################################################################
  !
  ! Find the unique y coordinates
  !
  y_ptr = init_array(npts,0) ! F2003 AUTO-REALLOC
  !
  y_vals(:) = zero
  !
  ny = 1
  y_vals(1) = my_sem_pts(1)%sem_fp(1)%y
  y_ptr(1) = 1
  !
  find_unique_y_loc_loop: &
  do n = 2,npts
    !
    do j = 1,ny
      if (abs(y_vals(j)-my_sem_pts(1)%sem_fp(n)%y) <= match_tolerance) then
        y_ptr(n) = j
        cycle find_unique_y_loc_loop
      end if
    end do
    !
    ny = ny + 1
    y_vals(ny) = my_sem_pts(1)%sem_fp(n)%y
    y_ptr(n) = ny
    !
  end do find_unique_y_loc_loop
  !
  !#############################################################################
  !
  ! Collect the averaged data into an array
  !
  y_counts = init_array(ny,0)    ! F2003 AUTO-REALLOC
  y        = init_array(ny,zero) ! F2003 AUTO-REALLOC
  uavg     = init_array(ny,zero) ! F2003 AUTO-REALLOC
  ravg     = init_array(ny,zero) ! F2003 AUTO-REALLOC
  tavg     = init_array(ny,zero) ! F2003 AUTO-REALLOC
  tau11    = init_array(ny,zero) ! F2003 AUTO-REALLOC
  tau22    = init_array(ny,zero) ! F2003 AUTO-REALLOC
  tau33    = init_array(ny,zero) ! F2003 AUTO-REALLOC
  tau21    = init_array(ny,zero) ! F2003 AUTO-REALLOC
  !
  do n = 1,npts
    !
    j = y_ptr(n)
    !
    y_counts(j) = y_counts(j) + 1
    !
    y(j)     = y(j)     + my_sem_pts(1)%sem_fp(n)%y
    uavg(j)  = uavg(j)  + my_sem_pts(1)%sem_fp(n)%uavg
    ravg(j)  = ravg(j)  + my_sem_pts(1)%sem_fp(n)%ravg
    tavg(j)  = tavg(j)  + my_sem_pts(1)%sem_fp(n)%tavg
    tau11(j) = tau11(j) + my_sem_pts(1)%sem_fp(n)%r11
    tau22(j) = tau22(j) + my_sem_pts(1)%sem_fp(n)%r22
    tau33(j) = tau33(j) + my_sem_pts(1)%sem_fp(n)%r33
    tau21(j) = tau21(j) + my_sem_pts(1)%sem_fp(n)%r21
    !
  end do
  !
  !#############################################################################
  !
  ! Average the data by dividing through by the number of averaged points
  !
  do j = 1,ny
    rnum = one / real(y_counts(j),kind=wp)
    y(j)     = y(j)     * rnum
    uavg(j)  = uavg(j)  * rnum
    ravg(j)  = ravg(j)  * rnum
    tavg(j)  = tavg(j)  * rnum
    tau11(j) = tau11(j) * rnum
    tau22(j) = tau22(j) * rnum
    tau33(j) = tau33(j) * rnum
    tau21(j) = tau21(j) * rnum
  end do
  !
end subroutine extract_sem_profile
!
!###############################################################################
!
subroutine write_profile_errors_file(x_vals,ix_ptr,xy_prof,sem_profile, &
                                     aux_data,aux_sem_data)
  !
  !.. Use Statements ..
  use ovar, only : itcur
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: x_vals(:)
  integer,  intent(in) :: ix_ptr(:)
  real(wp), intent(in) :: xy_prof(:,:,:)
  real(wp), intent(in) :: sem_profile(:,:)
  real(wp), intent(in) :: aux_data(:,:)
  real(wp), intent(in) :: aux_sem_data(:)
  !
  !.. Local Scalars ..
  integer  :: i,j,m,n,npts,nvar,io,ierr
  real(wp) :: dq
  !
  character(len=100) :: error_file
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: prof(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_profile_errors_file"
  !
continue
  !
  nvar = size(xy_prof,dim=1)
  npts = size(xy_prof,dim=2)
  !
  prof = init_array(nvar,npts,zero) ! F2003 AUTO-REALLOC
  !
  write (error_file,1) itcur
  error_file = trim(adjustl(error_file))
  !
  open (newunit=io,file=error_file,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,error_file,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (io,2)
  !
  do n = 1,size(ix_ptr)
    !
    if (n /= 1 .and. n /= 7) cycle
    !
    i = ix_ptr(n)
    !
    do j = 1,npts
      prof(1,j) = xy_prof(1,j,i)
      prof(2,j) = xy_prof(2,j,i)
      do m = 3,nvar
       !dq = sem_profile(m,j) - xy_prof(m,j,i)
        prof(m,j) = log10(max(eps12,abs(sem_profile(m,j)-xy_prof(m,j,i))))
      end do
    end do
    !
    write (io,3) nint(x_vals(i)),npts
    !
    write (io,11) "r_avg",         aux_data(1,i)
    write (io,12) "T_avg",         aux_data(2,i)
    write (io,13) "visc_avg",      aux_data(3,i)
    write (io,13) "dudy_wall",     aux_data(4,i)
    write (io,11) "tau_wall",      aux_data(5,i)
    write (io,11) "u_CL",          aux_data(6,i)
    write (io,11) "u_tau",         aux_data(7,i)
    write (io,11) "u_Bulk",        aux_data(8,i)
    write (io,13) "dstar",         aux_data(9,i)
    write (io,13) "theta",         aux_data(10,i)
    write (io,14) "Re_Bulk",  nint(aux_data(11,i))
    write (io,14) "Re_dstar", nint(aux_data(12,i))
    write (io,14) "Re_theta", nint(aux_data(13,i))
    write (io,14) "Re_CL",    nint(aux_data(14,i))
    write (io,14) "Re_tau",   nint(aux_data(15,i))
    write (io,13) "Cf_Bulk",       aux_data(16,i)
    write (io,13) "Cf_CL",         aux_data(17,i)
    !
    do j = 1,npts
      write (io,4) (prof(m,j), m=1,nvar)
    end do
    !
  end do
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,error_file,2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Format Statements
  !
  1 format ("spatial_ave.axial-profiles-error.",i0,".dat")
  2 format ('VARIABLES = "y+",',/, &
            '"y/<greek>d</greek>",',/, &
            '"u/u<sub>cl</sub>",',/, &
            '"u+",',/, &
            '"u<sub>rms</sub>/u<sub><greek>t</greek></sub>",',/, &
            '"v<sub>rms</sub>/u<sub><greek>t</greek></sub>",',/, &
            '"w<sub>rms</sub>/u<sub><greek>t</greek></sub>",',/, &
            '"u<math>\\034</math>v<math>\\034\</math>/', &
            'u<sub><greek>t</greek></sub><sup>2</sup>",',/, &
            '"<greek>t</greek><sub>xx</sub>",',/, &
            '"<greek>t</greek><sub>yy</sub>",',/, &
            '"<greek>t</greek><sub>zz</sub>",',/, &
            '"<greek>t</greek><sub>xy</sub>",',/, &
            '"u"')
  3 format ('ZONE T="Error in spatial average profile at x = ',i0, &
            ' H", I=',i0,', DATAPACKING=POINT')
  4 format (100es14.6)
  !
  11 format ('AUXDATA ',a,'="',f9.6,'"')
  12 format ('AUXDATA ',a,'="',f8.3,'"')
  13 format ('AUXDATA ',a,'="',es14.6,'"')
  14 format ('AUXDATA ',a,'="',i0,'"')
  !
end subroutine write_profile_errors_file
!
!###############################################################################
!
end module channel_mod
!
!###############################################################################
!
