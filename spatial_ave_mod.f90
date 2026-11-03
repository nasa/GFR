module spatial_ave_mod
  !
  ! Module containing procedures for spatially averaging the time-averaged
  ! solution read in from a restart file
  !
  !.. Use Statements ..
  use module_kind_types
  !
  implicit none
  !
  private
  !
  !.. Public procedures
  public :: spatial_ave_channel_flow
  public :: spatial_ave_z_planes
  public :: spatial_ave_x_and_z_planes
  !
  !.. Module Parameters ..
  real(wp), parameter :: match_tolerance = eps9
  !
  integer, parameter :: n_prof_data = 13
  integer, parameter :: n_aux_data = 19
  !
  !.. Module Allocatable Scalars ..
  character(len=:), allocatable, save :: spatial_ave_dir
  !
contains
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
  call spatial_ave_z_planes_with_constant_y(called_from_spatial_ave_mod=true)
  !
  call spatial_ave_x_and_z_planes(called_from_spatial_ave_mod=true)
  !
end subroutine spatial_ave_channel_flow
!
!###############################################################################
!
subroutine spatial_ave_x_and_z_planes(called_from_spatial_ave_mod)
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_solpts,xyz
  use flowvar_mod, only : uavesp
  use ovar,        only : itcur,time,time_ref,time_scaling_factor
  use ovar,        only : ave_start_time,output_time_averaging
  use ovar,        only : grid_scaling_factor
  use io_mod,      only : create_results_dir
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: called_from_spatial_ave_mod
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
  character(len=300) :: dir_result
  character(len=300) :: fname
  !
  !.. Local Arrays ..
  real(wp) :: y_vals(1:10000)
  real(wp) :: var_dim(1:size(uavesp,dim=2))
  real(wp) :: aux_data(1:n_aux_data)
  !
  character(len=50) :: varnames(1:size(uavesp,dim=2))
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: avefile
  character(len=:), allocatable :: proffile
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
  if (present(called_from_spatial_ave_mod)) then
    averaging_channel_flow = called_from_spatial_ave_mod
  end if
  !
  !#############################################################################
  !
  ! Create the Spatial_Ave_files directory if it that hasnt been done yet
  !
  if (.not.allocated(spatial_ave_dir)) then
    call create_results_dir("Spatial_Ave_files",dir_result,skip_lustre=true)
    spatial_ave_dir = trim(adjustl(dir_result)) ! F2003 AUTO-REALLOC
  end if
  !
  !#############################################################################
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
    rnum = one/(half*ymax)
  else
    rnum = one/grid_scaling_factor
  end if
  !
  !#############################################################################
  !
  ! Get the normalized profile
  !
  call get_normalized_profile(one,one, &
                              y_vals(1:ny), &
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
                              averaging_channel_flow,1, &
                              opt_debug_output=true)
  !
  !#############################################################################
  !
  ! Output the spatially-averaged time-averaged solution
  !
  write (fname,20) itcur
  avefile = spatial_ave_dir // trim(adjustl(fname)) ! F2003 AUTO-REALLOC
  !
  open (newunit=io,file=avefile,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,avefile,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (fname,40) itcur
  proffile = spatial_ave_dir // trim(adjustl(fname)) ! F2003 AUTO-REALLOC
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
  write (ip,53) "delta99",       aux_data(18)
  write (ip,53) "ywall",         aux_data(19)
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
subroutine spatial_ave_z_planes_with_constant_y(called_from_spatial_ave_mod)
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_solpts,xyz
  use flowvar_mod, only : uavesp
  use ovar,        only : itcur,time,time_ref,time_scaling_factor
  use ovar,        only : ave_start_time,output_time_averaging
  use ovar,        only : grid_scaling_factor
  use io_mod,      only : create_results_dir
  !
  !.. Optional Arguments ..
  logical(lk), optional, intent(in) :: called_from_spatial_ave_mod
  !
  !.. Local Scalars ..
  integer  :: i,j,m,n,ierr,io,ia1,ia2,ip
  integer  :: nvar,ny,nx,nix
  real(wp) :: out_time,rnum,xmin,xmax,ymax
  real(wp) :: sem_profile_height,sem_u_cl
  !
  integer  :: nvx,nvxvx,nvy,nvyvy,nvz,nvzvz,nvxvy,nrho,nt
  !
  logical(lk) :: averaging_channel_flow
  !
  character(len=50)  :: varname
  character(len=300) :: dir_result
  character(len=300) :: fname
  !
  !.. Local Arrays ..
  real(wp) :: x_vals(1:10000)
  real(wp) :: y_vals(1:10000)
  real(wp) :: var_dim(1:size(uavesp,dim=2))
  !
  character(len=50) :: varnames(1:size(uavesp,dim=2))
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: avefile
  character(len=:), allocatable :: axialfile
  character(len=:), allocatable :: axfullfile
  character(len=:), allocatable :: proffile
  character(len=:), allocatable :: auxfile
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
  character(len=*), parameter :: pname = "spatial_ave_z_planes_with_constant_y"
  !
continue
  !
  averaging_channel_flow = fals
  if (present(called_from_spatial_ave_mod)) then
    averaging_channel_flow = called_from_spatial_ave_mod
  end if
  !
  !#############################################################################
  !
  ! Create the Spatial_Ave_files directory if it that hasnt been done yet
  !
  if (.not.allocated(spatial_ave_dir)) then
    call create_results_dir("Spatial_Ave_files",dir_result,skip_lustre=true)
    spatial_ave_dir = trim(adjustl(dir_result)) ! F2003 AUTO-REALLOC
  end if
  !
  !#############################################################################
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
  xy_data   = init_array(nvar,ny,nx,zero)        ! F2003 AUTO-REALLOC
  xy_counts = init_array(ny,nx,0)                ! F2003 AUTO-REALLOC
  xy_dudy   = init_array(ny,nx,zero)             ! F2003 AUTO-REALLOC
  xy_prof   = init_array(n_prof_data,ny,nx,zero) ! F2003 AUTO-REALLOC
  aux_data  = init_array(n_aux_data,nx,zero)     ! F2003 AUTO-REALLOC
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
  ! Get the SEM profile
  !
  call get_sem_profile(sem_profile,aux_sem_data, &
                       sem_profile_height,sem_u_cl, &
                       averaging_channel_flow)
  !
  !#############################################################################
  !
  ! Get the spatially averaged profiles at axial locations
  ! equal to integer multiples of the channel half height.
  !
  if (averaging_channel_flow) then
    rnum = one/(half*ymax)
  else
    rnum = one/sem_profile_height
  end if
  !
  nix = nint( abs(xmax-xmin)*rnum ) + 1
  !
  ix_ptr = init_array(nix,0) ! F2003 AUTO-REALLOC
  !
  do n = 1,nix
    ix_ptr(n) = minloc( abs(abs(x_vals(1:nx)-xmin)*rnum - real(n-1,kind=wp)) , &
                        dim=1 )
  end do
  write (iout,*) "xmin = ",xmin
  write (iout,*) "xmax = ",xmax
  write (iout,*) "abs(xmax-xmin)*rnum = ",abs(xmax-xmin)*rnum
  write (iout,*) "nint(abs(xmax-xmin)*rnum) = ",nint(abs(xmax-xmin)*rnum)
  do n = 1,nix
    write (iout,201) n,ix_ptr(n),x_vals(ix_ptr(n)), &
                     (x_vals(ix_ptr(n))-xmin)*rnum
  end do
  201 format (2x,i2,": x_vals(",i0,") = ",f10.5,2x,f10.5)
  !
  !#############################################################################
  !
  ! Get the normalized profiles
  !
  do i = 1,nx
    call get_normalized_profile(sem_profile_height,sem_u_cl, &
                                y_vals(1:ny), &
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
  !#############################################################################
  !
  ! Output the spatially-averaged time-averaged solution
  !
  write (fname,20) itcur
  avefile = spatial_ave_dir // trim(adjustl(fname)) ! F2003 AUTO-REALLOC
  !
  open (newunit=io,file=avefile,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,avefile,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (fname,30) "",itcur
  axialfile = spatial_ave_dir // trim(adjustl(fname)) ! F2003 AUTO-REALLOC
  !
  open (newunit=ia1,file=axialfile,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,axialfile,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (fname,30) "-full",itcur
  axfullfile = spatial_ave_dir // trim(adjustl(fname)) ! F2003 AUTO-REALLOC
  !
  open (newunit=ia2,file=axfullfile,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,axfullfile,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (fname,40) itcur
  proffile = spatial_ave_dir // trim(adjustl(fname)) ! F2003 AUTO-REALLOC
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
    write (ip,53) "delta99",       aux_data(18,i)
    write (ip,53) "ywall",         aux_data(19,i)
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
    write (ip,53) "delta99",       aux_sem_data(18)
    write (ip,53) "ywall",         aux_sem_data(19)
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
  write (fname,50) itcur
  auxfile = spatial_ave_dir // trim(adjustl(fname)) ! F2003 AUTO-REALLOC
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
                         '"Cf_CL",',/, &
                         '"delta99",',/, &
                         '"ywall"',/, &
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
                          '"Cf_CL",',/, &
                          '"delta99",',/, &
                          '"ywall"',/, &
              'ZONE T="Axial Profile", I=',i0,', DATAPACKING=POINT')
  !
end subroutine spatial_ave_z_planes_with_constant_y
!
!###############################################################################
!
subroutine spatial_ave_z_planes()
  !
  !.. Use Statements ..
  use geomvar_mod, only : n_solpts,xyz
  use flowvar_mod, only : uavesp
  use ovar,        only : itcur,time,time_ref,time_scaling_factor
  use ovar,        only : ave_start_time,output_time_averaging
  use ovar,        only : grid_scaling_factor
  use io_mod,      only : create_results_dir
  !
  !.. Local Scalars ..
  integer  :: i,j,k,m,n,nj,ierr,io,ia1,ia2,ip
  integer  :: nvar,ny,nx,nz,nix
  real(wp) :: out_time,rnum,xmin,xmax,zmin,zmax
  real(wp) :: sem_profile_height,sem_u_cl
  !
  integer  :: nvx,nvxvx,nvy,nvyvy,nvz,nvzvz,nvxvy,nrho,nt
  !
  character(len=50)  :: varname
  character(len=300) :: dir_result
  character(len=300) :: fname
  !
  !.. Local Arrays ..
  real(wp) :: x_vals(1:10000)
  real(wp) :: z_vals(1:10000)
  real(wp) :: var_dim(1:size(uavesp,dim=2))
  !
  character(len=50) :: varnames(1:size(uavesp,dim=2))
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: proffile
  character(len=:), allocatable :: auxfile
  !
  !.. Local Allocatable Arrays ..
  integer,  allocatable :: x_ptr(:)
  integer,  allocatable :: y_ptr(:)
  integer,  allocatable :: z_ptr(:)
  integer,  allocatable :: ix_ptr(:)
  integer,  allocatable :: xy_counts(:,:)
  real(wp), allocatable :: j_vals(:)
  real(wp), allocatable :: y_vals(:,:)
  real(wp), allocatable :: xy_data(:,:,:)
  real(wp), allocatable :: xy_ta_data(:,:,:)
  real(wp), allocatable :: xy_time_ave(:,:,:)
  real(wp), allocatable :: xy_prof(:,:,:)
  real(wp), allocatable :: xy_dudy(:,:)
  real(wp), allocatable :: dudy(:)
  real(wp), allocatable :: aux_data(:,:)
  real(wp), allocatable :: sem_profile(:,:)
  real(wp), allocatable :: aux_sem_data(:)
  real(wp), allocatable :: x_nd(:)
  real(wp), allocatable :: y_nd(:,:)
  !
  character(len=32), allocatable :: ta_var_names(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "spatial_ave_z_planes"
  !
continue
  !
  !#############################################################################
  !
  ! Create the Spatial_Ave_files directory if it that hasnt been done yet
  !
  if (.not.allocated(spatial_ave_dir)) then
    call create_results_dir("Spatial_Ave_files",dir_result,skip_lustre=true)
    spatial_ave_dir = trim(adjustl(dir_result)) ! F2003 AUTO-REALLOC
  end if
  !
  !#############################################################################
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
  ! Find the unique z coordinates
  !
  z_ptr = init_array(n_solpts,0) ! F2003 AUTO-REALLOC
  !
  z_vals(:) = zero
  !
  nz = 1
  z_vals(1) = xyz(3,1)
  z_ptr(1) = 1
  !
  find_unique_z_loc_loop: &
  do n = 2,n_solpts
    !
    do k = 1,nz
      if (abs(z_vals(k)-xyz(3,n)) <= match_tolerance) then
        z_ptr(n) = k
        cycle find_unique_z_loc_loop
      end if
    end do
    !
    nz = nz + 1
    z_vals(nz) = xyz(3,n)
    z_ptr(n) = nz
    !
  end do find_unique_z_loc_loop
  !
  zmin = minval(z_vals(1:nz))
  zmax = maxval(z_vals(1:nz))
  !
  write (iout,300) "z",nz,"z",zmin,"z",zmax
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
  xmin = minval(x_vals(1:nx))
  xmax = maxval(x_vals(1:nx))
  !
  write (iout,300) "x",nx,"x",xmin,"x",xmax
  !
  !#############################################################################
  !
  ! Check to make sure this is a pure structured grid by seeing if the number
  ! of solution points is evenly divisible by the number of unique x coordinates
  ! and z coordinates
  !
  if (mod(n_solpts,nx*nz) > 0) then
    write (iout,999) n_solpts, nx, nz
    return
  else
    ny = n_solpts/(nx*nz)
    write (iout,998) n_solpts,nx,ny,nz
  end if
  !
  !#############################################################################
  !
  ! Find the y coordinates for each unique x coordinate
  !
  y_ptr = init_array(n_solpts,0)  ! F2003 AUTO-REALLOC
  y_vals = init_array(ny,nx,zero) ! F2003 AUTO-REALLOC
  j_vals = init_array(ny,zero)    ! F2003 AUTO-REALLOC
  xy_counts = init_array(ny,nx,0)  ! F2003 AUTO-REALLOC
  !
  do i = 1,nx
    !
    nj = 0
    j_vals(:) = zero
    !
    find_unique_y_loc_loop: &
    do n = 1,n_solpts
      !
      if (x_ptr(n) /= i) cycle
      !
      if (nj == 0) then
        !
        nj = 1
        y_ptr(n) = nj
        j_vals(nj) = xyz(2,n)
        y_vals(nj,i) = xyz(2,n)
        xy_counts(nj,i) = xy_counts(nj,i) + 1
        !
      else
        !
        do j = 1,nj
          if (abs(j_vals(j)-xyz(2,n)) <= match_tolerance) then
            y_ptr(n) = j
            xy_counts(j,i) = xy_counts(j,i) + 1
            cycle find_unique_y_loc_loop
          end if
        end do
        !
        nj = nj + 1
        y_ptr(n) = nj
        j_vals(nj) = xyz(2,n)
        y_vals(nj,i) = xyz(2,n)
        xy_counts(nj,i) = xy_counts(nj,i) + 1
        !
      end if
      !
    end do find_unique_y_loc_loop
    !
  end do
  !
  write (iout,301) nz,minval(xy_counts),maxval(xy_counts)
  !
  !#############################################################################
  !
  ! Collect the averaged data into an x-y structured array
  !
  xy_data   = init_array(nvar,ny,nx,zero)        ! F2003 AUTO-REALLOC
  xy_counts = init_array(ny,nx,0)                ! F2003 AUTO-REALLOC
  xy_dudy   = init_array(ny,nx,zero)             ! F2003 AUTO-REALLOC
  xy_prof   = init_array(n_prof_data,ny,nx,zero) ! F2003 AUTO-REALLOC
  aux_data  = init_array(n_aux_data,nx,zero)     ! F2003 AUTO-REALLOC
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
  ! Get the time averaged variables from the time averaged restart data
  !
  call get_time_ave_variables(xy_data,xy_time_ave,ta_var_names)
  !
  !#############################################################################
  !
  ! Get the SEM profile
  !
  call get_sem_profile(sem_profile,aux_sem_data, &
                       sem_profile_height,sem_u_cl, &
                       averaging_channel_flow=fals)
  !
  !#############################################################################
  !
  ! Get the spatially averaged profiles at axial locations
  ! equal to integer multiples of the grid scaling factor
  !
  rnum = one/sem_profile_height
  !
  nix = nint( abs(xmax-xmin)*rnum ) + 1
  !
  ix_ptr = init_array(nix,0) ! F2003 AUTO-REALLOC
  !
  do n = 1,nix
    ix_ptr(n) = minloc( abs(abs(x_vals(1:nx)-xmin)*rnum - real(n-1,kind=wp)) , &
                        dim=1 )
  end do
  write (iout,*) "xmin = ",xmin
  write (iout,*) "xmax = ",xmax
  write (iout,*) "abs(xmax-xmin)*rnum = ",abs(xmax-xmin)*rnum
  write (iout,*) "nint(abs(xmax-xmin)*rnum) = ",nint(abs(xmax-xmin)*rnum)
  do n = 1,nix
    write (iout,201) n,ix_ptr(n),x_vals(ix_ptr(n)), &
                     (x_vals(ix_ptr(n))-xmin)*rnum
  end do
  !
  !#############################################################################
  !
  ! Get the normalized profiles
  !
  do i = 1,nx
    call get_normalized_profile(sem_profile_height,sem_u_cl, &
                                y_vals(1:ny,i), &
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
                                averaging_channel_flow=fals,ix=i)
  end do
  !
  !#############################################################################
  !
  ! Dimensionalize the time-averaged restart data in xy_data
  !
  xy_ta_data = init_array(ny,nx,nvar,zero) ! F2003 AUTO-REALLOC
  !
  do i = 1,nx
    do j = 1,ny
      do m = 1,nvar
        xy_ta_data(j,i,m) = xy_data(m,j,i)*var_dim(m)
      end do
    end do
  end do
  !
  ! Save the x and y coordinates nondimensionalized by rnum
  !
  x_nd = x_vals(1:nx)*rnum   ! F2003 AUTO-REALLOC
  y_nd = xy_prof(2,1:ny,1:nx) ! F2003 AUTO-REALLOC
  !
  !#############################################################################
  !
  ! Output the spatial-average of the quantities in the
  ! time-averaged restart file
  !
  call write_spatial_ave_ta_file(varnames,"ta_data", &
                                 x_vals,x_nd,y_vals,y_nd, &
                                 xy_ta_data,aux_data,ix_ptr)
  !
  !#############################################################################
  !
  ! Output the spatial-average of the fluctuating quantities that are
  ! derived from the time-averaged data
  !
  call write_spatial_ave_ta_file(ta_var_names,"time_ave", &
                                 x_vals,x_nd,y_vals,y_nd, &
                                 xy_time_ave,aux_data,ix_ptr)
  !
  !#############################################################################
  !
  ! Write the spatial_ave.axial-profiles.#.dat file
  !
  write (fname,40) itcur
  proffile = spatial_ave_dir // trim(adjustl(fname)) ! F2003 AUTO-REALLOC
  !
  open (newunit=ip,file=proffile,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,proffile,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (ip,41)
  !
  ! Write the normalized y-profile at each axial location equal to
  ! an integer interval of the channel height
  !
  do n = 1,nix
    !
    i = ix_ptr(n)
    !
    write (ip,35) nint((x_vals(i)-xmin)*rnum),ny
    !
    call write_auxdata_section(ip,aux_data(:,i),x_vals(i))
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
    !
    call write_auxdata_section(ip,aux_sem_data,xmin)
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
  !#############################################################################
  !
  ! Write the spatial_ave.aux-data.#.dat file
  !
  write (fname,50) itcur
  auxfile = spatial_ave_dir // trim(adjustl(fname)) ! F2003 AUTO-REALLOC
  !
  open (newunit=io,file=auxfile,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,auxfile,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (io,60) nx
  !
  do i = 1,nx
    write (io,13) x_vals(i),x_nd(i), &
                  x_vals(i)/aux_data(18,i), &
                  (aux_data(m,i), m=1,size(aux_data,dim=1))
  end do
  !
  close (io,iostat=ierr,iomsg=error_message)
  call io_error(pname,auxfile,2,__LINE__,__FILE__,ierr,error_message)
  !
  !#############################################################################
  !
  ! Write the profile errors file
  !
  if (size(xy_prof,dim=2) == size(sem_profile,dim=2)) then
    call write_profile_errors_file(x_nd,ix_ptr,xy_prof,sem_profile, &
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
  13 format (100es14.6)
  !
  21 format (/,2x,"Number of iterations completed = ",i0,/, &
               2x,"Number of flow-through-times averaged = ",f9.4,/)
  !
  34 format ('ZONE T="SEM Inflow Profile", I=',i0,', DATAPACKING=POINT')
  35 format ('ZONE T="Spatial average profile at x/d = ',i0, &
             '", I=',i0,', DATAPACKING=POINT')
  !
  40 format ("spatial_ave.nondim-profiles.",i0,".dat")
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
  !
  60 format ('VARIABLES = "x",',/, &
                         '"x/deltaSEM",',/, &
                         '"x/delta99",',/, &
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
                         '"Cf_CL",',/, &
                         '"delta99",',/, &
                         '"ywall"',/, &
             'ZONE T="Axial Profile", I=',i0,', DATAPACKING=POINT')
  !
  201 format (2x,i2,": x_vals(",i0,") = ",f10.5,2x,f10.5)
  !
  300 format (/,2x,"Number of unique ",a1,"-coordinate locations = ",i0,/, &
                36x,a1,"-min = ",es22.15,/, &
                36x,a1,"-max = ",es22.15,/)
  301 format (/,21x,"Actual number of z-planes = ",i0,/, &
                2x,"Minimum # of z-planes found for any x-y pair = ",i0,/, &
                2x,"Maximum # of z-planes found for any x-y pair = ",i0,/)
  !
  998 format (/,4x,"n_solpts = ",i0,/, &
                4x,"      nx = ",i0,/, &
                4x,"      ny = ",i0,/, &
                4x,"      nz = ",i0,/)
  !
  999 format (/,80("#"),/, &
              4x,"ERROR ATTEMPTING TO SPATIAL AVERAGE THE Z-PLANE DATA!",/, &
              4x,"THE GRID DOES NOT APPEAR TO BE TRUE I-J STRUCTURED GRID.",/, &
              4x,"BECAUSE THE NUMBER OF SOLUTION POINTS IS NOT EVENLY",/, &
              4x,"DIVISIBLE BY THE NUMBER OF UNIQUE X-COORDINATES AND",/, &
              4x,"Z-COORDINTES.",/, &
              8x,"     NUMBER OF SOLUTION POINTS = ",i0,/, &
              8x,"NUMBER OF UNIQUE X COORDINATES = ",i0,/, &
              8x,"NUMBER OF UNIQUE Z COORDINATES = ",i0,/, &
              80("#"),//)
  !
end subroutine spatial_ave_z_planes
!
!###############################################################################
!
subroutine write_spatial_ave_ta_file(vnames,cdata,xdim,x,ydim,y, &
                                     xy_data,aux_data,ix_ptr)
  !
  !.. Use Statements ..
  use ovar, only : itcur
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: vnames(:)
  character(len=*), intent(in) :: cdata
  real(wp),         intent(in) :: xdim(:)
  real(wp),         intent(in) :: x(:)
  real(wp),         intent(in) :: ydim(:,:)
  real(wp),         intent(in) :: y(:,:)
  real(wp),         intent(in) :: xy_data(:,:,:)
  real(wp),         intent(in) :: aux_data(:,:)
  integer,          intent(in) :: ix_ptr(:)
  !
  !.. Local Scalars ..
  integer :: i,j,n,m,nx,ny,nvar,ierr
  integer :: io_zp,io_ai,io_af
  character(len=300) :: fname
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: zplane_file
  character(len=:), allocatable :: axial_file_intervals
  character(len=:), allocatable :: axial_file_full
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: varname
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "write_spatial_ave_ta_file"
  !
continue
  !
  nx = size(x)
  ny = size(y,dim=1)
  nvar = size(vnames)
  !
  !#############################################################################
  !
  ! Write the averaged x-y planes file
  !
  write (fname,1) cdata,itcur
  zplane_file = spatial_ave_dir // trim(adjustl(fname)) ! F2003 AUTO-REALLOC
  !
  open (newunit=io_zp,file=zplane_file,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,zplane_file,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (io_zp,10) '"X/d",','"X",'
  !
  do m = 1,nvar-1
    write (io_zp,11) trim(adjustl(vnames(m))),","
  end do
  write (io_zp,11) trim(adjustl(vnames(nvar)))
  !
  write (io_zp,21) nx,ny
  !
  do j = 1,ny
    do i = 1,nx
      write (io_zp,30) x(i), y(j,i), xdim(i), ydim(j,i), &
                       (xy_data(j,i,m), m=1,nvar)
    end do
  end do
  !
  close (io_zp,iostat=ierr,iomsg=error_message)
  call io_error(pname,zplane_file,2,__LINE__,__FILE__,ierr,error_message)
  !
  !#############################################################################
  !
  ! Write the averaged y-profiles at each axial location equal to
  ! an integer interval of the channel height
  !
  write (fname,2) "",cdata,itcur
  ! F2003 AUTO-REALLOC
  axial_file_intervals = spatial_ave_dir // trim(adjustl(fname))
  !
  open (newunit=io_ai,file=axial_file_intervals, &
                      status="replace",action="write", &
                      iostat=ierr,iomsg=error_message)
  call io_error(pname,axial_file_intervals,1,__LINE__,__FILE__, &
                ierr,error_message)
  !
  write (io_ai,10) '',''
  !
  do m = 1,nvar-1
    write (io_ai,11) trim(adjustl(vnames(m))),","
  end do
  write (io_ai,11) trim(adjustl(vnames(nvar)))
  !
  do n = 1,size(ix_ptr)
    !
    i = ix_ptr(n)
    !
    write (io_ai,22) nint(x(i)),ny
    !
    call write_auxdata_section(io_ai,aux_data(:,i),xdim(i))
    !
    do j = 1,ny
      write (io_ai,30) y(j,i), ydim(j,i), (xy_data(j,i,m), m=1,nvar)
    end do
    !
  end do
  !
  close (io_ai,iostat=ierr,iomsg=error_message)
  call io_error(pname,axial_file_intervals,2,__LINE__,__FILE__, &
                ierr,error_message)
  !
  !#############################################################################
  !
  ! Write the averaged y-profiles at every axial location
  !
  write (fname,2) "-full",cdata,itcur
  axial_file_full = spatial_ave_dir // trim(adjustl(fname)) ! F2003 AUTO-REALLOC
  !
  open (newunit=io_af,file=axial_file_full,status="replace",action="write", &
                   iostat=ierr,iomsg=error_message)
  call io_error(pname,axial_file_full,1,__LINE__,__FILE__,ierr,error_message)
  !
  write (io_af,10) '',''
  !
  do m = 1,nvar-1
    write (io_af,11) trim(adjustl(vnames(m))),","
  end do
  write (io_af,11) trim(adjustl(vnames(nvar)))
  !
  do i = 1,nx
    !
    write (io_af,23) x(i),ny
    !
    call write_auxdata_section(io_af,aux_data(:,i),xdim(i))
    !
    do j = 1,ny
      write (io_af,30) y(j,i), ydim(j,i), (xy_data(j,i,m), m=1,nvar)
    end do
    !
  end do
  !
  close (io_af,iostat=ierr,iomsg=error_message)
  call io_error(pname,axial_file_full,2,__LINE__,__FILE__,ierr,error_message)
  !
  ! Format Statements
  !
  1  format ("spatial_ave.z-plane.",a,".",i0,".dat")
  2  format ("spatial_ave.dimen-profiles",a,".",a,".",i0,".dat")
  !
  10 format ('VARIABLES = ',a,'"Y/d",',a,'"Y",')
  11 format ('"',a,'"',:,a)
  !
  21 format ('ZONE T="Spatial Average", I=',i0,', J=',i0,', DATAPACKING=POINT')
  22 format ('ZONE T="Spatial average profile at X/d = ',i0, &
             '", I=',i0,', DATAPACKING=POINT')
  23 format ('ZONE T="Spatial average profile at X/d =',f6.3, &
             '", I=',i0,', DATAPACKING=POINT')
  !
  30 format (100es14.6)
  !
end subroutine write_spatial_ave_ta_file
!
!###############################################################################
!
subroutine write_auxdata_section(io,aux_data,x)
  !
  !.. Formal Arguments ..
  integer,  intent(in) :: io
  real(wp), intent(in) :: aux_data(:)
  real(wp), intent(in) :: x
  !
continue
  !
  write (io,1) "r_avg",         aux_data(1)
  write (io,2) "T_avg",         aux_data(2)
  write (io,3) "visc_avg",      aux_data(3)
  write (io,3) "dudy_wall",     aux_data(4)
  write (io,3) "tau_wall",      aux_data(5)
  write (io,2) "u_CL",          aux_data(6)
  write (io,2) "u_tau",         aux_data(7)
  write (io,2) "u_Bulk",        aux_data(8)
  write (io,3) "dstar",         aux_data(9)
  write (io,3) "theta",         aux_data(10)
  write (io,4) "Re_Bulk",  nint(aux_data(11))
  write (io,4) "Re_dstar", nint(aux_data(12))
  write (io,4) "Re_theta", nint(aux_data(13))
  write (io,4) "Re_CL",    nint(aux_data(14))
  write (io,4) "Re_tau",   nint(aux_data(15))
  write (io,3) "Cf_Bulk",       aux_data(16)
  write (io,3) "Cf_CL",         aux_data(17)
  write (io,3) "delta99",       aux_data(18)
  write (io,3) "ywall",         aux_data(19)
  write (io,3) "x_loc",         x
  !
  1 format ('AUXDATA ',a,'="',f9.6,'"')
  2 format ('AUXDATA ',a,'="',f8.3,'"')
  3 format ('AUXDATA ',a,'="',es14.6,'"')
  4 format ('AUXDATA ',a,'="',i0,'"')
  !
end subroutine write_auxdata_section
!
!###############################################################################
!
subroutine get_time_ave_variables(ta_data,ta_var,ta_var_names)
  !
  !.. Use Statements ..
  use geomvar_mod,    only : nr
  use cgns_types_mod, only : CGLEN
  use write_cgns_mod, only : get_cgns_time_ave_variables
  use write_cgns_mod, only : count_output_variables
  !
  !.. Formal Arguments ..
  real(wp),                       intent(in) :: ta_data(:,:,:)
  real(wp), allocatable,         intent(out) :: ta_var(:,:,:)
  character(len=*), allocatable, intent(out) :: ta_var_names(:)
  !
  !.. Local Scalars ..
  integer :: l,m,n,nvar,nx,ny,ierr
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: varname
  !
  !.. Local Allocatable Arrays ..
#ifdef SPECIAL_FOR_INTEL
  character(len=:), allocatable :: variable_list(:)
#else
  character(len=CGLEN), allocatable :: variable_list(:)
#endif
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_time_ave_variables"
  !
  character(len=*), parameter :: char_xyz(1:3) = ["X","Y","Z"]
  character(len=*), parameter :: char_cyl(1:3) = ["R    ","Theta","Phi  "]
  !
continue
  !
  call get_cgns_time_ave_variables(variable_list)
  !
  nvar = count_output_variables(variable_list)
  !
  ny = size(ta_data,dim=2)
  nx = size(ta_data,dim=3)
  !
  ta_var = init_array(ny,nx,nvar,zero) ! F2003 AUTO-REALLOC
  !
  allocate ( ta_var_names(1:nvar) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"ta_var_names",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !
  ! Get the complete list of variables that need to be output, which includes
  ! expanding certain variables like Momentum and Velocity into their
  ! individual vector strings.
  !
  n = 0
  do l = 1,size(variable_list)
    !
    select case(uppercase(trim(adjustl(variable_list(l)))))
      !
      case ("MOMENTUM")
        !
        do m = 1,nr
          n = n + 1
          varname = "Momentum"//char_xyz(m) ! F2003 AUTO-REALLOC
          call compute_time_ave_variable(varname,ta_data,ta_var(:,:,n))
          ta_var_names(n) = varname
        end do
        !
      case ("VELOCITY")
        !
        do m = 1,nr
          n = n + 1
          varname = "Velocity"//char_xyz(m) ! F2003 AUTO-REALLOC
          call compute_time_ave_variable(varname,ta_data,ta_var(:,:,n))
          ta_var_names(n) = varname
        end do
        !
      case ("CYLINDRICALVELOCITY")
        !
        do m = 1,2
          n = n + 1
          varname = "Velocity"//trim(adjustl(char_cyl(m))) ! F2003 AUTO-REALLOC
          call compute_time_ave_variable(varname,ta_data,ta_var(:,:,n))
          ta_var_names(n) = varname
        end do
        !
      case default
        !
        n = n + 1
        varname = trim(adjustl(variable_list(l))) ! F2003 AUTO-REALLOC
        call compute_time_ave_variable(varname,ta_data,ta_var(:,:,n))
        ta_var_names(n) = varname
        !
    end select
    !
  end do
  !
end subroutine get_time_ave_variables
!
!###############################################################################
!
subroutine compute_time_ave_variable(varname,ta_data,var)
  !
  !.. Use Statements ..
  use ovar, only : rhoref,aref,tref
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: varname
  real(wp),         intent(in) :: ta_data(:,:,:)
  real(wp),      intent(inout) :: var(:,:)
  !
  !.. Local Scalars ..
  integer  :: idx1,idx2,idx3
  real(wp) :: ref_coef
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: vname
  !
  !.. Local Saved Scalars ..
  integer,     save :: n_d = 0
  integer,     save :: n_x = 0
  integer,     save :: n_y = 0
  integer,     save :: n_z = 0
  integer,     save :: n_r = 0
  integer,     save :: n_f = 0
  integer,     save :: n_t = 0
  logical(lk), save :: sub_needs_initialization = true
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_time_ave_variable"
  !
continue
  !
  if (sub_needs_initialization) then
    !
    call get_ta_index("d",n_d)
    call get_ta_index("x",n_x)
    call get_ta_index("y",n_y)
    call get_ta_index("z",n_z)
    call get_ta_index("r",n_r)
    call get_ta_index("f",n_f)
    call get_ta_index("t",n_t)
    !
    sub_needs_initialization = fals
    !
  end if
  !
  idx1 = 0
  idx2 = -99
  idx3 = -99
  !
  vname = trim(adjustl(varname)) ! F2003 AUTO-REALLOC
  !
  select case (vname)
    !
    case ("Density")
      idx1 = n_d
      !
    case ("VelocityX")
      idx1 = n_x
    case ("VelocityY")
      idx1 = n_y
    case ("VelocityZ")
      idx1 = n_z
    case ("VelocityR")
      idx1 = n_r
    case ("VelocityTheta")
      idx1 = n_f
      !
    case ("MomentumX")
      call get_ta_index("dx",idx1)
    case ("MomentumY")
      call get_ta_index("dy",idx1)
    case ("MomentumZ")
      call get_ta_index("dz",idx1)
    case ("MomentumR")
      call get_ta_index("dr",idx1)
    case ("MomentumTheta")
      call get_ta_index("df",idx1)
      !
    case ("Temperature")
      idx1 = n_t
    case ("Pressure")
      call get_ta_index("p",idx1)
      !
    case ("T0")
      call get_ta_index("T",idx1)
    case ("P0")
      call get_ta_index("P",idx1)
      !
    case ("uu")
      call get_ta_index("xx",idx1) ; idx2 = n_x
    case ("vv")
      call get_ta_index("yy",idx1) ; idx2 = n_y
    case ("ww")
      call get_ta_index("zz",idx1) ; idx2 = n_z
    case ("rr")
      call get_ta_index("rr",idx1) ; idx2 = n_r
    case ("ff")
      call get_ta_index("ff",idx1) ; idx2 = n_f
      !
    case ("uv")
      call get_ta_index("xy",idx1) ; idx2 = n_x ; idx3 = n_y
    case ("uw")
      call get_ta_index("xz",idx1) ; idx2 = n_x ; idx3 = n_z
    case ("vw")
      call get_ta_index("yz",idx1) ; idx2 = n_y ; idx3 = n_z
    case ("ur")
      call get_ta_index("xr",idx1) ; idx2 = n_x ; idx3 = n_r
    case ("uf")
      call get_ta_index("xf",idx1) ; idx2 = n_x ; idx3 = n_f
    case ("rf")
      call get_ta_index("rf",idx1) ; idx2 = n_r ; idx3 = n_f
      !
    case ("TT")
      call get_ta_index("tt",idx1) ; idx2 = n_t
    case ("uT")
      call get_ta_index("xt",idx1) ; idx2 = n_x ; idx3 = n_t
    case ("vT")
      call get_ta_index("yt",idx1) ; idx2 = n_y ; idx3 = n_t
    case ("wT")
      call get_ta_index("zt",idx1) ; idx2 = n_z ; idx3 = n_t
    case ("rT")
      call get_ta_index("rt",idx1) ; idx2 = n_r ; idx3 = n_t
    case ("fT")
      call get_ta_index("ft",idx1) ; idx2 = n_f ; idx3 = n_t
      !
    case ("dd")
      call get_ta_index("dd",idx1) ; idx2 = n_d
    case ("du")
      call get_ta_index("dx",idx1) ; idx2 = n_d ; idx3 = n_x
    case ("dv")
      call get_ta_index("dy",idx1) ; idx2 = n_d ; idx3 = n_y
    case ("dw")
      call get_ta_index("dz",idx1) ; idx2 = n_d ; idx3 = n_z
    case ("dr")
      call get_ta_index("dr",idx1) ; idx2 = n_d ; idx3 = n_r
    case ("df")
      call get_ta_index("df",idx1) ; idx2 = n_d ; idx3 = n_f
    case ("dT")
      call get_ta_index("dt",idx1) ; idx2 = n_d ; idx3 = n_t
      !
    case ("PP")
      call get_ta_index("pp",idx1)
      call get_ta_index("p",idx2)
      !
    case default
      idx1 = n_d
      !
  end select
  !
  ! Get the coefficient to dimensionalize the variable
  !
  select case (vname(1:2))
    case ("De")
      ref_coef = rhoref
    case ("Ve")
      ref_coef = aref
    case ("Mo","du","dv","dw","dr","df")
      ref_coef = rhoref*aref
    case ("Te","T0")
      ref_coef = tref
    case ("Pr","P0")
      ref_coef = rhoref*aref*aref
    case ("uu","vv","ww","rr","ff","uv","uw","vw","ur","uf","rf")
      ref_coef = aref*aref
    case ("TT")
      ref_coef = tref**2
    case ("uT","vT","wT","rT","fT")
      ref_coef = aref*tref
    case ("dd")
      ref_coef = rhoref*rhoref
    case ("dT")
      ref_coef = rhoref*tref
    case ("PP")
      ref_coef = (rhoref*aref*aref)**2
    case default
      ref_coef = rhoref
  end select
  !
  ! Now use the found indices to actually compute this time averaged variable
  !
  if (all([idx2,idx3] /= -99)) then
    call compute_ta_variable(ta_data,var,ref_coef,idx1,idx2,idx3)
  else if (idx2 /= -99 .and. idx3 == -99) then
   !call compute_ta_variable(ta_data,var,idx1,idx2)
    call compute_ta_variable(ta_data,var,ref_coef,idx1,idx2,idx2)
  else
    call compute_ta_variable(ta_data,var,ref_coef,idx1)
  end if
  !
end subroutine compute_time_ave_variable
!
!###############################################################################
!
pure &
subroutine get_ta_index(string,idx)
  !
  !.. Use Statements ..
  use flowvar_mod, only : time_ave_variables
  !
  !.. Formal Arguments ..
  character(len=*), intent(in) :: string
  integer,       intent(inout) :: idx
  !
  !.. Local Scalars ..
  integer :: m
  !
continue
  !
  if (idx == 0) then
    !
    do m = 1,size(time_ave_variables)
      if (trim(adjustl(time_ave_variables(m))) == string) then
        idx = m
        exit
      end if
    end do
    !
  end if
  !
end subroutine get_ta_index
!
!###############################################################################
!
pure &
subroutine compute_ta_variable(ta_data,var,ref_coef,index1,index2,index3)
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: ta_data(:,:,:)
  real(wp), intent(inout) :: var(:,:)
  real(wp),    intent(in) :: ref_coef
  integer,     intent(in) :: index1
  !
  !.. Optional Arguments ..
  integer, optional, intent(in) :: index2
  integer, optional, intent(in) :: index3
  !
  !.. Local Scalars ..
  integer  :: i,j,n,lb,ub
  integer  :: idx1,idx2,idx3
  real(wp) :: coef,val
  !
continue
  !
  lb = lbound(ta_data,dim=1)
  ub = ubound(ta_data,dim=1)
  !
  var(:,:) = zero
  !
  one_pass_loop: &
  do n = 1,1
    !
    if (index1 < lb .or. index1 > ub) then
      ! exit the subroutine leaving var=0 if index1 isnt valid
      exit one_pass_loop
    else
      idx1 = index1
    end if
    !
    if (present(index2)) then
      if (index2 < lb .or. index2 > ub) then
        ! exit the subroutine leaving var=0 if index2 isnt valid
        exit one_pass_loop
      else
        idx2 = index2
        coef = one
      end if
    else
      ! if index2 isnt present, set idx2 to a valid index, but make the
      ! leading coefficient zero
      idx2 = lb
      coef = zero
    end if
    !
    if (present(index3)) then
      ! NOTE: If index2 is not present, index3 will get zeroed out by coef
      !       even if it is a valid index.
      if (index3 < lb .or. index3 > ub) then
        ! exit the subroutine leaving var=0 if index3 isnt valid
        exit one_pass_loop
      else
        idx3 = index3
      end if
    else
      ! set the idx3 the idx2 in case index2 is present and valid,
      ! and it is intended to be the square of that index
      idx3 = idx2
    end if
    !
    do j = 1,size(var,dim=2)
      do i = 1,size(var,dim=1)
        val = ta_data(idx1,i,j) - coef*ta_data(idx2,i,j)*ta_data(idx3,i,j)
        var(i,j) = ref_coef*val
      end do
    end do
    !
  end do one_pass_loop
  !
end subroutine compute_ta_variable
!
!###############################################################################
!
pure &
subroutine get_variable_names(varnames,var_dim,nvx,nvxvx,nvy,nvyvy, &
                              nvz,nvzvz,nvxvy,nrho,nt)
  !
  !.. Use Statements ..
  use ovar,        only : rhoref,aref,tref,tkeref,tlsref
  use flowvar_mod, only : time_ave_variables
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
      ta_char = time_ave_variables(m)(l:l)
      select case (ta_char)
        case ("d")
          varname = trim(adjustl(varname)) // trim(vmult) // "Rho"
          var_dim(m) = var_dim(m) * rhoref
        case ("x")
          varname = trim(adjustl(varname)) // trim(vmult) // "Vx"
          var_dim(m) = var_dim(m) * aref
        case ("y")
          varname = trim(adjustl(varname)) // trim(vmult) // "Vy"
          var_dim(m) = var_dim(m) * aref
        case ("z")
          varname = trim(adjustl(varname)) // trim(vmult) // "Vz"
          var_dim(m) = var_dim(m) * aref
        case ("e")
          varname = trim(adjustl(varname)) // trim(vmult) // "E"
          var_dim(m) = var_dim(m) * rhoref * aref * aref
        case ("p")
          varname = trim(adjustl(varname)) // trim(vmult) // "P"
          var_dim(m) = var_dim(m) * rhoref * aref * aref
        case ("t")
          varname = trim(adjustl(varname)) // trim(vmult) // "T"
          var_dim(m) = var_dim(m) * tref
        case ("P")
          varname = trim(adjustl(varname)) // trim(vmult) // "P0"
          var_dim(m) = var_dim(m) * rhoref * aref * aref
        case ("T")
          varname = trim(adjustl(varname)) // trim(vmult) // "T0"
          var_dim(m) = var_dim(m) * tref
        case ("r")
          varname = trim(adjustl(varname)) // trim(vmult) // "Vr"
          var_dim(m) = var_dim(m) * aref
        case ("f")
          varname = trim(adjustl(varname)) // trim(vmult) // "Vf"
          var_dim(m) = var_dim(m) * aref
        case ("k")
          varname = trim(adjustl(varname)) // trim(vmult) // "TKE"
          var_dim(m) = var_dim(m) * tkeref
        case ("w")
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
  use geomvar_mod,     only : ncell,cell
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
subroutine get_normalized_profile(profile_height,profile_uref, &
                                  y_vals,u,uu,v,vv,w,ww,uv,rho,t,dudy, &
                                  prof,aux_data,averaging_channel_flow,ix, &
                                  opt_debug_output)
  !
  !.. Use Statements ..
  use order_mod,        only : n_order
  use ovar,             only : rhoref,aref,tref
  use ovar,             only : grid_scaling_factor
  use element_type_mod, only : std_elem
  use transport_mod,    only : compute_viscosity_cpg
  use interpolation_mod, only : interp,General_Interp_Matrix
  !
  !.. Formal Arguments ..
  real(wp),    intent(in) :: profile_height
  real(wp),    intent(in) :: profile_uref
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
  real(wp) :: xi_min,xi_max,y1,y2,yc,yh,ywall,tau_wall
  real(wp) :: dy,wt,udt,umt,udtm,udtp,umtm,umtp
  real(wp) :: displacement_thickness,momentum_thickness
  real(wp) :: u_bulk_tr,displacement_thickness_tr,momentum_thickness_tr
  real(wp) :: delta99,u99
  !
  logical(lk) :: debug_output
  !
  real(wp) :: ymin,ymax,length
  real(wp) :: dsplt,mmntt,ublk
  real(wp), allocatable :: ycell(:)
  real(wp), allocatable :: xifp(:,:)
  real(wp), allocatable :: mat(:,:)
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
  ymin = y_vals(1)
  !
  ! The following is needed in case the wall is not at y=0
  !
  ! Find center of first cell off the wall
  if (mod(n_order+1,2) == 0) then
    yc = half*( y_vals( (n_order+1)/2 ) + y_vals( (n_order+1)/2 + 1 ) )
  else
    yc = y_vals( (n_order/2) + 1 )
  end if
  ! Find height of first cell off the wall
  yh = two*( (y_vals(1)-yc) / std_elem(Geom_Edge,n_order)%pts(1,1) )
  ! Get the coordinates of the wall
  ywall = yc - half*yh
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
   !u_cl = profile_uref
    !
    delta99 = height
    u99 = u_cl
    do n = 1,npts
      if (u(n) >= 0.99_wp*u_cl) then
        delta99 = y_vals(n)-ywall
        u99 = u(n)
        exit
      end if
    end do
    !
    half_height = delta99
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
    prof( 1,n) = (y_vals(n)-ywall) * u_tau * rho_avg / visc_avg ! y+
    prof( 2,n) = (y_vals(n)-ywall) / half_height / gsf          ! y/d
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
  ! Use trapezoid rule to integrate bulk velocity
  !
  u_bulk_tr = zero
  displacement_thickness_tr = zero
  momentum_thickness_tr = zero
  do n = 1,npts-1
    dy = abs(y_vals(n+1)-y_vals(n))
    u_bulk_tr = u_bulk_tr + half*(u(n) + u(n+1))*dy
    udtm = (one - u(n  )/u_cl)
    udtp = (one - u(n+1)/u_cl)
    umtm = udtm*u(n  )/u_cl
    umtp = udtp*u(n+1)/u_cl
    displacement_thickness_tr = displacement_thickness_tr + half*(udtm+udtp)*dy
    momentum_thickness_tr = momentum_thickness_tr + half*(umtm+umtp)*dy
  end do
  !
  ! Use quadrature to integrate the bulk velocity
  !
  u_bulk = zero
  momentum_thickness = zero
  displacement_thickness = zero
  !
  ublk = zero
  dsplt = zero
  mmntt = zero
  !
  if (mod(npts,n1d) == 0) then
    !
    xifp = std_elem(Geom_Edge,n1d-1)%face_pts(:,:) ! F2003 AUTO-REALLOC
    xifp(1,1) = -one
    xifp(1,2) =  one
    mat = General_Interp_Matrix(Geom_Edge,n1d-1,xifp) ! F2003 AUTO-REALLOC
    !
    xi_min = std_elem(Geom_Edge,n_order)%pts(1,1)
    xi_max = std_elem(Geom_Edge,n_order)%pts(1,n1d)
    !
    do nc = 1,npts/n1d
      !
      n0 = (nc-1)*n1d
      !
      ycell = y_vals(n0+1:n0+n1d) ! F2003 AUTO-REALLOC
      !
      ymin = dot_product(ycell,mat(:,1))
      ymax = dot_product(ycell,mat(:,2))
      length = abs(ymax-ymin)
      !
      y1 = y_vals(n0+1)
      y2 = y_vals(n0+n1d)
      !
      rnum = (y2 - y1) / (xi_max - xi_min) ! Simplified using Mathematica
      !
      do n = 1,n1d
        !
        wt = half*length*std_elem(Geom_Edge,n_order)%wts(n)
        u_bulk = u_bulk + u(n0+n)*wt
        udt = (one - u(n0+n)/u_cl)
        umt = udt*u(n0+n)/u_cl
        displacement_thickness = displacement_thickness + udt*wt
        momentum_thickness = momentum_thickness + umt*wt
        !
        wt = rnum*std_elem(Geom_Edge,n_order)%wts(n)
        ublk = ublk + u(n0+n)*wt
        udt = (one - u(n0+n)/u_cl)
        umt = udt*u(n0+n)/u_cl
        dsplt = dsplt + udt*wt
        mmntt = mmntt + umt*wt
        !
      end do
      !
    end do
    !
  else
    !
    u_bulk = u_bulk_tr
    displacement_thickness = displacement_thickness_tr
    momentum_thickness = momentum_thickness_tr
    !
  end if
  !
  u_bulk = u_bulk / height
  u_bulk_tr = u_bulk_tr / height
  ublk = ublk / height
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
  ! BL delta99 or channel half-height
  aux_data(18) = half_height
  ! y location of the wall (useful when it is nonzero)
  aux_data(19) = ywall
  !
  if (debug_output) then
    if (mod(npts,n1d) == 0) then
      write (iout,*)
      write (iout,801) "u_Bulk",u_bulk*aref
      write (iout,901) " dstar",displacement_thickness
      write (iout,901) " theta",momentum_thickness
      write (iout,*)
      write (iout,803) "u_Bulk",ublk*aref
      write (iout,903) " dstar",dsplt
      write (iout,903) " theta",mmntt
    end if
    write (iout,*)
    write (iout,802) "u_Bulk",u_bulk_tr*aref
    write (iout,902) " dstar",displacement_thickness_tr
    write (iout,902) " theta",momentum_thickness_tr
  end if
  801 format (2x,"Using     Quadrature:  ",a6," = ",f11.6," m/s")
  901 format (2x,"Using     Quadrature:  ",a6," = ",es14.6)
  802 format (2x,"Using Trapezoid Rule:  ",a6," = ",f11.6," m/s")
  902 format (2x,"Using Trapezoid Rule:  ",a6," = ",es14.6)
  803 format (2x,"Using Old Quadrature:  ",a6," = ",f11.6," m/s")
  903 format (2x,"Using Old Quadrature:  ",a6," = ",es14.6)
  !
end subroutine get_normalized_profile
!
!###############################################################################
!
!pure &
subroutine get_sem_profile(prof,aux_data,profile_height,profile_uref, &
                           averaging_channel_flow)
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
  real(wp),            intent(inout) :: profile_height
  real(wp),            intent(inout) :: profile_uref
  logical(lk),            intent(in) :: averaging_channel_flow
  !
  !.. Local Scalars ..
  integer  :: n,npts,ierr
  real(wp) :: rnum,height,half_height,dudy_wall,gsf
  real(wp) :: rho_avg,t_avg,visc_avg,u_cl,u_tau,u_bulk
  real(wp) :: tau_wall,yc,yh,ywall,ymin
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
  profile_height = one/grid_scaling_factor
  !
  call extract_sem_profile(y,u,r,t,tau11,tau22,tau33,tau21,ierr)
  !
  if (ierr /= 0) return
  !
  npts = size(y)
  !
  rnum = one / real(npts,kind=wp)
  !
  prof     = init_array(n_prof_data,npts,zero) ! F2003 AUTO-REALLOC
  aux_data = init_array(n_aux_data,zero)       ! F2003 AUTO-REALLOC
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
    profile_height = half_height
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
    profile_height = height
    !
    rho_avg = r(1)
    t_avg = t(1)
    !
    u_cl = u(npts)
    !
  end if
  !
  profile_uref = u_cl
  !
  !
  ! The following is needed in case the wall is not at y=0
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
  write (iout,*)
  write (iout,2) "     y(1)",y(1)
  write (iout,2) "    ywall",ywall
  write (iout,2) "       dy",dy
  write (iout,2) "     u(1)",u(1)
  write (iout,2) "dudy_wall",dudy_wall
  write (iout,*)
  write (iout,1) "y loc of wall",ywall
  write (iout,1) "y loc for center of first cell",yc
  write (iout,1) "max y loc for sem face points",y(npts)
  write (iout,*)
  write (iout,1) "height of first cell",yh
  write (iout,1) "wall-adjusted y loc for center of first cell",yc
  write (iout,1) "wall-adjusted max y loc for sem face points",y(npts)
  write (iout,*)
  1 format (a," = ",f15.10)
  2 format (a," = ",f20.10)
  !
  write (iout,*)
  do n = 1,npts
    write (iout,10) n,y(n)-ywall,y(n),u(n),u(n)*aref
  end do
  write (iout,*)
  10 format ("y(",i4,") = ",f15.10,5x,f15.10,5x,f17.10,5x,f17.10)
  !
  call compute_viscosity_cpg(t_avg,visc_avg)
  !
  tau_wall = visc_avg * dudy_wall
  !
  u_tau = sqrt( tau_wall / rho_avg )
  !
  do n = 1,npts
    prof( 1,n) = (y(n)-ywall) * u_tau * rho_avg / visc_avg  ! y+
    prof( 2,n) = (y(n)-ywall) / profile_height              ! y/d
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
  aux_data(14) = rho_avg * u_cl * profile_height / visc_avg
  ! Re_tau = Reynolds number based on frication velocity
  aux_data(15) = rho_avg * u_tau * profile_height / visc_avg
  ! Skin friction based on bulk velocity (Pope eq. 7.15)
  aux_data(16) = tau_wall / (half * rho_avg * u_bulk * u_bulk)
  ! Skin friction based on center-line velocity (Pope eq. 7.14)
  aux_data(17) = tau_wall / (half * rho_avg * u_cl * u_cl)
  ! BL delta99 or channel half-height
  aux_data(18) = profile_height
  ! y location of the wall (useful when it is nonzero)
  aux_data(19) = ywall
  !
end subroutine get_sem_profile
!
!###############################################################################
!
!pure &
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
  integer,     allocatable :: y_ptr(:)
  integer,     allocatable :: y_counts(:)
  integer,     allocatable :: y_sort(:)
  logical(lk), allocatable :: not_sorted(:)
  real(wp),    allocatable :: tmp_data(:,:)
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
  write (iout,*)
  write (iout,*) "npts = ",npts
  write (iout,*) "minval(y_counts) = ",minval(y_counts,dim=1)
  write (iout,*) "maxval(y_counts) = ",maxval(y_counts,dim=1)
  write (iout,*)
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
  !#############################################################################
  !
  ! Create an indexing array that gives the indices of the y values
  ! sorted from smallest to largest
  !
  y_sort = init_array(ny,0)        ! F2003 AUTO-REALLOC
  not_sorted = init_array(ny,true) ! F2003 AUTO-REALLOC
  !
  do j = 1,ny
    ! Find the minimum value in array y masked by the true values
    ! of the not_sorted array
    n = minloc(y,mask=not_sorted,dim=1)
    ! Save this index as the minimum remaining not-sorted y value
    y_sort(j) = n
    ! Mark this location of the y array as not-not-sorted (i.e., sorted)
    not_sorted(n) = fals
  end do
  !
  !#############################################################################
  !
  ! Finally, sort the arrays
  !
  tmp_data = init_array(ny,8,zero) ! F2003 AUTO-REALLOC
  !
  do j = 1,ny
    !
    n = y_sort(j)
    !
    tmp_data(j,1) = y(n)
    tmp_data(j,2) = uavg(n)
    tmp_data(j,3) = ravg(n)
    tmp_data(j,4) = tavg(n)
    tmp_data(j,5) = tau11(n)
    tmp_data(j,6) = tau22(n)
    tmp_data(j,7) = tau33(n)
    tmp_data(j,8) = tau21(n)
    !
  end do
  !
  y(1:ny) = tmp_data(1:ny,1)
  uavg(1:ny) = tmp_data(1:ny,2)
  ravg(1:ny) = tmp_data(1:ny,3)
  tavg(1:ny) = tmp_data(1:ny,4)
  tau11(1:ny) = tmp_data(1:ny,5)
  tau22(1:ny) = tmp_data(1:ny,6)
  tau33(1:ny) = tmp_data(1:ny,7)
  tau21(1:ny) = tmp_data(1:ny,8)
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
  character(len=300) :: fname
  !
  !.. Local Allocatable Scalars ..
  character(len=:), allocatable :: error_file
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
  write (fname,1) itcur
  error_file = spatial_ave_dir // trim(adjustl(fname)) ! F2003 AUTO-REALLOC
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
    write (io,13) "delta99",       aux_data(18,i)
    write (io,13) "ywall",         aux_data(19,i)
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
end module spatial_ave_mod
!
!###############################################################################
!
